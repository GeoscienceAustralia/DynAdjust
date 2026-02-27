//============================================================================
// Name         : dnadiff.cpp
// Author       : Dale Roberts <dale.o.roberts@gmail.com>
// Copyright    : Copyright 2025 Geoscience Australia
//
//                Licensed under the Apache License, Version 2.0 (the
//                "License"); you may not use this file except in compliance
//                with the License. You may obtain a copy of the License at
//
//                http ://www.apache.org/licenses/LICENSE-2.0
//
//                Unless required by applicable law or agreed to in writing,
//                software distributed under the License is distributed on an
//                "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND,
//                either express or implied. See the License for the specific
//                language governing permissions and limitations under the
//                License.
//
// Description  : DynAdjust file comparison with tolerances
//============================================================================

#include "dnadiff.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <regex>
#include <sstream>
#include <vector>

namespace dynadjust {
namespace dnadiff {

DnaDiff::DnaDiff(): difference_count_(0), total_comparisons_(0) {}

DnaDiff::~DnaDiff() {}

bool DnaDiff::IsNumeric(const std::string& str) const {
    if (str.empty()) return false;

    std::regex numeric_pattern(R"(^[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$)");
    return std::regex_match(str, numeric_pattern);
}

double DnaDiff::ParseDouble(const std::string& str) const {
    try {
        return std::stod(str);
    } catch (...) { return 0.0; }
}

bool DnaDiff::CompareValues(double val1, double val2) const {
    return std::abs(val1 - val2) <= options_.tolerance;
}

std::string DnaDiff::NormalisePath(const std::string& path) const {
    try {
        namespace fs = std::filesystem;
        fs::path p(path);
        return fs::weakly_canonical(p).lexically_normal().string();
    } catch (...) {
        return path;
    }
}

bool DnaDiff::CompareTokens(const std::string& token1,
                            const std::string& token2) const {
    // First check if tokens are identical
    if (token1 == token2) return true;

    // If both are numeric, compare with tolerance
    if (IsNumeric(token1) && IsNumeric(token2)) {
        double val1 = ParseDouble(token1);
        double val2 = ParseDouble(token2);
        return CompareValues(val1, val2);
    }

    // For non-numeric tokens, they must match exactly
    return false;
}

bool DnaDiff::ShouldSkipLine(const std::string& line) const {
    // Skip empty lines
    if (line.empty())
        return false;  // Don't skip - both files should have same empty lines

    // Always skip timestamp and version lines that will differ
    if (line.find("File created:") != std::string::npos) return true;
    if (line.find("Build:") != std::string::npos) return true;
    if (line.find("Version:") != std::string::npos) return true;
    if (line.find("time") != std::string::npos) return true;

    // Always ignore paths-related lines
    if (line.find("File name:") != std::string::npos) return true;
    if (line.find("Input files:") != std::string::npos) return true;
    if (line.find("Output folder:") != std::string::npos) return true;
    if (line.find("Input folder:") != std::string::npos) return true;
    if (line.find("Command line arguments:") != std::string::npos)
        return true;

    // Skip thread-related lines
    if (line.find("threads") != std::string::npos) return true;

    // Skip station correction lines
    if (line.find("Maximum station correction") != std::string::npos) return true;

    // Skip coordinate adjustment lines
    if (line.find("(e, n, up)") != std::string::npos) return true;

    return false;
}



bool DnaDiff::SkipToMarker(std::ifstream& f, const std::string& marker) const {
    std::string line;
    while (std::getline(f, line)) {
        if (line.find(marker) != std::string::npos)
            return true;
    }
    return false;
}

bool DnaDiff::CompareResultsOnly(const std::string& file1,
                                 const std::string& file2) {
    std::ifstream f1(file1);
    std::ifstream f2(file2);
    std::string line1, line2;
    int line_count = 0;

    if (!options_.skip_to_marker.empty()) {
        if (!SkipToMarker(f1, options_.skip_to_marker)) {
            std::cerr << "Error: Marker \"" << options_.skip_to_marker
                      << "\" not found in " << file1_name_ << std::endl;
            return false;
        }
        if (!SkipToMarker(f2, options_.skip_to_marker)) {
            std::cerr << "Error: Marker \"" << options_.skip_to_marker
                      << "\" not found in " << file2_name_ << std::endl;
            return false;
        }
    }

    // Skip headers in both files
    for (int i = 0; i < options_.skip_header_lines; ++i) {
        if (f1.good()) std::getline(f1, line1);
        if (f2.good()) std::getline(f2, line2);
        line_count++;
    }

    // Compare line by line
    while (std::getline(f1, line1) && std::getline(f2, line2)) {
        line_count++;

        // Skip lines that should be ignored
        if (ShouldSkipLine(line1) && ShouldSkipLine(line2)) {
            continue;
        }

        // Skip empty lines
        if (line1.empty() && line2.empty()) continue;

        total_comparisons_++;

        // Compare tokens on each line
        std::istringstream iss1(line1);
        std::istringstream iss2(line2);
        std::string token1, token2;
        std::vector<std::string> tokens1, tokens2;

        // Collect all tokens
        while (iss1 >> token1) tokens1.push_back(token1);
        while (iss2 >> token2) tokens2.push_back(token2);

        bool line_matches = true;

        // Check if this line might contain a file path
        bool likely_path_line = (line1.find("model:") != std::string::npos ||
                                 line1.find("file:") != std::string::npos ||
                                 line1.find(".gsb") != std::string::npos ||
                                 line1.find(".dat") != std::string::npos ||
                                 line1.find(".msr") != std::string::npos ||
                                 line1.find(".asl") != std::string::npos ||
                                 line1.find(".xyz") != std::string::npos ||
                                 line1.find("./") != std::string::npos ||
                                 line1.find("/") != std::string::npos);

        // Compare tokens
        size_t min_tokens = std::min(tokens1.size(), tokens2.size());
        for (size_t i = 0; i < min_tokens; ++i) {
            // If both are numeric, compare with tolerance
            if (IsNumeric(tokens1[i]) && IsNumeric(tokens2[i])) {
                double val1 = ParseDouble(tokens1[i]);
                double val2 = ParseDouble(tokens2[i]);
                if (!CompareValues(val1, val2)) {
                    line_matches = false;
                    if (options_.verbose) {
                        std::cout << "Line " << line_count << " - numeric difference: "
                                  << val1 << " vs " << val2
                                  << " (diff=" << (val1 - val2) << ")\n";
                    }
                    break;  // Stop at first difference in this line
                }
            } else if (tokens1[i] != tokens2[i]) {
                // Try normalising as paths if on a path line
                bool paths_match = false;
                if (likely_path_line) {
                    std::string norm1 = NormalisePath(tokens1[i]);
                    std::string norm2 = NormalisePath(tokens2[i]);
                    paths_match = (norm1 == norm2);
                }

                if (!paths_match) {
                    // Non-numeric tokens must match exactly
                    line_matches = false;
                    if (options_.verbose) {
                        std::cout << "Line " << line_count << " - text difference: "
                                  << tokens1[i] << " vs " << tokens2[i] << "\n";
                    }
                    break;  // Stop at first difference in this line
                }
            }
        }

        // Check if one line has more tokens
        if (tokens1.size() != tokens2.size()) {
            line_matches = false;
            if (options_.verbose) {
                std::cout << "Line " << line_count << " - different number of tokens\n";
            }
        }

        if (!line_matches) {
            difference_count_++;
            if (options_.verbose) {
                std::cout << "  " << file1_name_ << ": " << line1 << "\n";
                std::cout << "  " << file2_name_ << ": " << line2 << "\n";
            }

            // Early exit for 3-token lines with colon (like "Reference frame: GDA2020")
            if (tokens1.size() == 3 && tokens2.size() == 3 &&
                line1.find(':') != std::string::npos) {
                if (options_.verbose) {
                    std::cout << "\nFundamental difference found. Stopping comparison.\n";
                }
                f1.close();
                f2.close();
                return false;  // Files differ
            }
        }
    }

    // Check if one file has more lines
    if (std::getline(f1, line1) || std::getline(f2, line2)) {
        difference_count_++;
        if (options_.verbose) {
            std::cout << "Files have different number of lines\n";
        }
    }

    f1.close();
    f2.close();

    return difference_count_ == 0;
}



bool DnaDiff::CompareFiles(const std::string& file1, const std::string& file2,
                           const DiffOptions& options) {
    options_ = options;
    difference_count_ = 0;
    total_comparisons_ = 0;
    differences_.clear();

    namespace fs = std::filesystem;

    fs::path path1(file1);
    fs::path path2(file2);

    if (!fs::exists(path1)) {
        std::cerr << "Error: File does not exist: " << path1 << std::endl;
        return false;
    }

    if (!fs::exists(path2)) {
        std::cerr << "Error: File does not exist: " << path2 << std::endl;
        return false;
    }

    file1_name_ = path1.filename().string();
    file2_name_ = path2.filename().string();

    // Only compare measurement results
    return CompareResultsOnly(path1.string(), path2.string());
}

void DnaDiff::PrintSummary() const {
    std::cout << "\nComparison Summary:\n";
    std::cout << "  Total lines compared: " << total_comparisons_ << "\n";
    std::cout << "  Differences found: " << difference_count_ << "\n";

    if (difference_count_ > 0 && !differences_.empty()) {
        std::cout << "\nFirst "
                  << std::min(10, static_cast<int>(differences_.size()))
                  << " differences:\n";
        for (size_t i = 0; i < std::min(size_t(10), differences_.size()); ++i) {
            std::cout << "  " << differences_[i] << "\n";
        }
    }

    std::cout << "\nTolerance used: " << options_.tolerance << "\n";
}

bool ParseCommandLine(int argc, char* argv[], std::string& file1,
                      std::string& file2, DiffOptions& options) {
    if (argc < 3) { return false; }

    file1 = argv[1];
    file2 = argv[2];

    for (int i = 3; i < argc; ++i) {
        std::string arg = argv[i];

        if ((arg == "--tolerance" || arg == "--tol" || arg == "-t") &&
            i + 1 < argc) {
            options.tolerance = std::stod(argv[++i]);
        } else if (arg == "--skip-headers" && i + 1 < argc) {
            options.skip_header_lines = std::stoi(argv[++i]);
        } else if (arg == "--skip-to-marker" && i + 1 < argc) {
            options.skip_to_marker = argv[++i];
        } else if (arg == "--verbose" || arg == "-v") {
            options.verbose = true;
        } else if (arg == "--help" || arg == "-h") {
            return false;
        }
    }

    return true;
}

void PrintUsage(const std::string& program_name) {
    std::cout << "Usage: " << program_name << " file1 file2 [options]\n\n";
    std::cout << "Compare two DynAdjust output files with tolerance for "
                 "numerical values.\n\n";
    std::cout << "Options:\n";
    std::cout << "  --tolerance, --tol, -t <value>  Tolerance for numeric "
                 "comparisons (default: 1e-6)\n";
    std::cout << "  --skip-headers <n>              Skip first n lines "
                 "(default: 0)\n";
    std::cout << "  --skip-to-marker <text>         Skip lines in each file "
                 "until marker text is found\n";
    std::cout << "  --verbose, -v                   Show detailed differences\n";
    std::cout << "  --help, -h                      Show this help message\n";
}

}  // namespace dnadiff
}  // namespace dynadjust

int main(int argc, char* argv[]) {
    using namespace dynadjust::dnadiff;

    std::string file1, file2;
    DiffOptions options;

    if (!ParseCommandLine(argc, argv, file1, file2, options)) {
        PrintUsage(argv[0]);
        return 1;
    }

    DnaDiff differ;
    bool files_match = differ.CompareFiles(file1, file2, options);

    if (options.verbose || !files_match) { differ.PrintSummary(); }

    if (files_match) {
        std::cout << "Files match within tolerance.\n";
        return 0;
    } else {
        std::cout << "Files differ beyond tolerance.\n";
        return 1;
    }
}
