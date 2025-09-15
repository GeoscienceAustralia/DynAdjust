//============================================================================
// Name         : asl_file.cpp
// Author       : Roger Fraser
// Contributors : Dale Roberts <dale.o.roberts@gmail.com>
// Copyright    : Copyright 2017-2025 Geoscience Australia
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
// Description  : DynAdjust associated station file io operations
//============================================================================

#include <include/io/asl_file.hpp>

namespace dynadjust {
namespace iostreams {

AslFile::AslFile(const std::filesystem::path& filename)
    : DynadjustFile(), path_(filename) {}

AslFile::AslFile(const AslFile& asl)
    : DynadjustFile(asl), path_(asl.path_) {}

AslFile& AslFile::operator=(const AslFile& rhs) {
  if (this == &rhs) {
    return *this;
  }
  DynadjustFile::operator=(rhs);
  path_ = rhs.path_;
  return *this;
}

AslLoadResult AslFile::Load() {
  AslLoadResult result;
  result.count = LoadLegacy(&result.stations, &result.free_stations);
  return result;
}

std::optional<AslLoadResult> AslFile::TryLoad() {
  try {
    return Load();
  } catch (const std::runtime_error&) {
    return std::nullopt;
  }
}

std::uint64_t AslFile::LoadLegacy(measurements::vASL* binary_asl,
                                  vUINT32* free_stn) {
  std::ifstream asl_file;
  std::ostringstream os;
  os << "LoadLegacy(): An error was encountered when opening " << path_
     << "." << std::endl;

  try {
    file_opener(asl_file, path_.string(), std::ios::in | std::ios::binary,
                binary, true);
  } catch (const std::runtime_error& e) {
    os << e.what();
    throw std::runtime_error(os.str());
  } catch (...) {
    throw std::runtime_error(os.str());
  }

  os.str("");
  os << "LoadLegacy(): An error was encountered when reading from " <<
      path_ << "." << std::endl;

  std::uint64_t stn_count;

  try {
    ReadFileInfo(asl_file);

    asl_file.read(reinterpret_cast<char*>(&stn_count), sizeof(std::uint64_t));

    initialiseIncrementingIntegerVector(free_stn,
                                        static_cast<UINT32>(stn_count));

    binary_asl->clear();
    binary_asl->resize(stn_count);

    for (auto& asl_entry : *binary_asl) {
      asl_file >> &asl_entry;
    }
  } catch (const std::ios_base::failure& f) {
    os << f.what();
    throw std::runtime_error(os.str());
  } catch (const std::runtime_error& e) {
    os << e.what();
    throw std::runtime_error(os.str());
  } catch (...) {
    throw std::runtime_error(os.str());
  }

  asl_file.close();
  return stn_count;
}

void AslFile::Write(const measurements::vASLPtr& binary_asl) {
  std::ofstream asl_file;
  std::ostringstream os;
  os << "Write(): An error was encountered when opening " << path_
     << "." << std::endl;

  try {
    file_opener(asl_file, path_.string(), std::ios::out | std::ios::binary,
                binary);
  } catch (const std::runtime_error& e) {
    os << e.what();
    throw std::runtime_error(os.str());
  } catch (...) {
    throw std::runtime_error(os.str());
  }

  os.str("");
  os << "Write(): An error was encountered when writing to " << path_
     << "." << std::endl;

  std::uint64_t asl_count = static_cast<std::uint64_t>(binary_asl.size());

  try {
    WriteFileInfo(asl_file);

    asl_file.write(reinterpret_cast<char*>(&asl_count), sizeof(std::uint64_t));
    for (const auto& asl_entry : binary_asl) {
      if (asl_entry.get()) {
        asl_file << asl_entry.get();
      }
    }
  } catch (const std::ios_base::failure& f) {
    os << f.what();
    throw std::runtime_error(os.str());
  } catch (const std::runtime_error& e) {
    os << e.what();
    throw std::runtime_error(os.str());
  } catch (...) {
    throw std::runtime_error(os.str());
  }

  asl_file.close();
}

void AslFile::WriteText(const measurements::vASLPtr& binary_asl,
                       const measurements::vdnaStnPtr& stations) {
  std::ofstream asl_file;
  std::ostringstream os;
  os << "WriteText(): An error was encountered when opening " <<
      path_ << "." << std::endl;

  try {
    file_opener(asl_file, path_.string());
  } catch (const std::runtime_error& e) {
    os << e.what();
    throw std::runtime_error(os.str());
  } catch (...) {
    throw std::runtime_error(os.str());
  }

  os.str("");
  os << "WriteText(): An error was encountered when writing to " <<
      path_ << "." << std::endl;

  size_t asl_count = binary_asl.size();

  std::stringstream ss_asl;
  ss_asl << asl_count << " stations";
  asl_file << std::left << std::setw(STATION) << ss_asl.str();
  asl_file << std::setw(HEADER_20) << std::right << "No. connected msrs";
  asl_file << std::setw(STATION) << std::right << "AML index";
  asl_file << std::setw(STATION) << std::right << "Unused?" << std::endl;

  vUINT32 asl_ptrs(binary_asl.size());
  initialiseIncrementingIntegerVector<UINT32>(
      asl_ptrs, static_cast<UINT32>(binary_asl.size()));

  CompareMeasCount2<measurements::ASLPtr, UINT32>
      msr_count_compare_func(&binary_asl);
  std::sort(asl_ptrs.begin(), asl_ptrs.end(), msr_count_compare_func);

  binary_asl.at(0).get()->GetAMLStnIndex();

  try {
    for (const auto& asl_ptr : asl_ptrs) {
      if (!binary_asl.at(asl_ptr)) {
        continue;
      }

      auto it_pasl = binary_asl.begin() + asl_ptr;

      asl_file << std::setw(STATION) << std::left
               << stations.at(asl_ptr)->GetName()
               << std::setw(HEADER_20) << std::right
               << it_pasl->get()->GetAssocMsrCount();

      if (it_pasl->get()->GetAssocMsrCount() == 0) {
        asl_file << std::setw(STATION) << "-";
      } else {
        asl_file << std::setw(STATION) << std::right
                 << it_pasl->get()->GetAMLStnIndex();
      }

      asl_file << std::setw(STATION) << std::right
               << (it_pasl->get()->IsValid() ? " " : "*") << std::endl;
    }
  } catch (const std::ios_base::failure& f) {
    os << f.what();
    throw std::runtime_error(os.str());
  } catch (const std::runtime_error& e) {
    os << e.what();
    throw std::runtime_error(os.str());
  } catch (...) {
    throw std::runtime_error(os.str());
  }

  asl_file.close();
}

}  // namespace iostreams
}  // namespace dynadjust
