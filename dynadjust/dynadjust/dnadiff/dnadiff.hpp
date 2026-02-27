//============================================================================
// Name         : dnadiff.hpp
// Author       : Dale Roberts <dale.o.roberts@gmail.com>
// Copyright    : Copyright 2025 Geoscience Australia
//
//                Licensed under the Apache License, Version 2.0 (the "License");
//                you may not use this file except in compliance with the License.
//                You may obtain a copy of the License at
//
//                http ://www.apache.org/licenses/LICENSE-2.0
//
//                Unless required by applicable law or agreed to in writing, software
//                distributed under the License is distributed on an "AS IS" BASIS,
//                WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//                See the License for the specific language governing permissions and
//                limitations under the License.
//
// Description  : DynAdjust file comparison with tolerances
//============================================================================

#ifndef DNADIFF_HPP_
#define DNADIFF_HPP_

#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

namespace dynadjust {
namespace dnadiff {

struct DiffOptions {
    double tolerance = 1e-6;  // Single tolerance for all numeric comparisons
    int skip_header_lines = 0;
    std::string skip_to_marker;
    bool verbose = false;
};

class DnaDiff {
   public:
    DnaDiff();
    ~DnaDiff();

    // Main comparison function
    bool CompareFiles(const std::string& file1, const std::string& file2, const DiffOptions& options);

    // Get summary of differences
    void PrintSummary() const;
    int GetDifferenceCount() const { return difference_count_; }

   private:
    // Helper functions
    bool IsNumeric(const std::string& str) const;
    double ParseDouble(const std::string& str) const;
    bool CompareValues(double val1, double val2) const;
    bool CompareTokens(const std::string& token1, const std::string& token2) const;
    std::string NormalisePath(const std::string& path) const;

    // Skip non-comparable content
    bool ShouldSkipLine(const std::string& line) const;
    bool SkipToMarker(std::ifstream& f, const std::string& marker) const;

    // Main comparison function
    bool CompareResultsOnly(const std::string& file1, const std::string& file2);

    // Member variables
    DiffOptions options_;
    int difference_count_;
    int total_comparisons_;
    std::vector<std::string> differences_;
    std::string file1_name_;
    std::string file2_name_;
};

// Command line parsing
bool ParseCommandLine(int argc, char* argv[], std::string& file1, std::string& file2, DiffOptions& options);
void PrintUsage(const std::string& program_name);

}  // namespace dnadiff
}  // namespace dynadjust

#endif  // DNADIFF_HPP_
