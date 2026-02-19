//============================================================================
// Name         : bms_file_loader.cpp
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
// Description  : DynAdjust binary measurement file io operations
//============================================================================

#include <include/io/bms_file_loader.hpp>

#include <algorithm>
#include <fstream>
#include <ios>
#include <iostream>
#include <map>
#include <sstream>

#include <include/functions/dnaiostreamfuncs.hpp>

namespace dynadjust {
namespace iostreams {

BmsFileLoader& BmsFileLoader::operator=(const BmsFileLoader& rhs) {
  if (this == &rhs) {
    return *this;
  }
  dna_io_base::operator=(rhs);
  return *this;
}

std::uint64_t BmsFileLoader::CreateMsrInputFileMeta(
    vifm_t& vinput_file_meta,
    input_file_meta_t** input_file_meta) {
  std::uint64_t msr_file_count(0);
  std::stringstream ss;
  ss << "CreateMsrInputFileMeta(): An error was encountered when creating "
     << std::endl
     << "  the binary measurement file metadata." << std::endl;

  try {
    std::for_each(vinput_file_meta.begin(), vinput_file_meta.end(),
                  [&msr_file_count](input_file_meta_t& ifm) {
                    if (ifm.datatype == msr_data ||
                        ifm.datatype == stn_msr_data) {
                      msr_file_count++;
                    }
                  });

    if (*input_file_meta != nullptr) {
      delete[] * input_file_meta;
    }

    (*input_file_meta) = new input_file_meta_t[msr_file_count];

    msr_file_count = 0;
    std::for_each(vinput_file_meta.begin(), vinput_file_meta.end(),
                  [&msr_file_count, &input_file_meta](input_file_meta_t& ifm) {
                    if (ifm.datatype == msr_data ||
                        ifm.datatype == stn_msr_data) {
                      (*input_file_meta)[msr_file_count] = ifm;
                      msr_file_count++;
                    }
                  });
  } catch (...) {
    throw std::runtime_error(ss.str());
  }

  return msr_file_count;
}

void BmsFileLoader::LoadFileMeta(const std::string& bms_filename,
                                 binary_file_meta_t& bms_meta) {
  std::ifstream bms_file;
  std::stringstream ss;
  ss << "LoadFileMeta(): An error was encountered when opening " << bms_filename
     << "." << std::endl;

  try {
    file_opener(bms_file, bms_filename, std::ios::in | std::ios::binary,
                binary, true);
  } catch (const std::runtime_error& e) {
    ss << e.what();
    throw std::runtime_error(ss.str());
  } catch (...) {
    throw std::runtime_error(ss.str());
  }

  ss.str("");
  ss << "LoadFileMeta(): An error was encountered when reading from "
     << bms_filename << "." << std::endl;

  try {
    ReadFileInfo(bms_file);
    ReadFileMetadata(bms_file, bms_meta);
  } catch (const std::ios_base::failure& f) {
    ss << f.what();
    throw std::runtime_error(ss.str());
  } catch (const std::runtime_error& e) {
    ss << e.what();
    throw std::runtime_error(ss.str());
  } catch (...) {
    throw std::runtime_error(ss.str());
  }

  bms_file.close();
}

std::uint64_t BmsFileLoader::LoadFile(const std::string& bms_filename,
                                      measurements::pvmsr_t vbinary_msr,
                                      binary_file_meta_t& bms_meta) {
  std::ifstream bms_file;
  std::stringstream ss;
  ss << "LoadFile(): An error was encountered when opening " << bms_filename
     << "." << std::endl;

  try {
    file_opener(bms_file, bms_filename, std::ios::in | std::ios::binary,
                binary, true);
  } catch (const std::runtime_error& e) {
    ss << e.what();
    throw std::runtime_error(ss.str());
  } catch (...) {
    throw std::runtime_error(ss.str());
  }

  measurements::measurement_t measRecord;
  std::uint64_t msr;

  ss.str("");
  ss << "LoadFile(): An error was encountered when reading from "
     << bms_filename << "." << std::endl;

  try {
    ReadFileInfo(bms_file);
    ReadFileMetadata(bms_file, bms_meta);

    vbinary_msr->reserve(bms_meta.binCount);
    for (msr = 0; msr < bms_meta.binCount; msr++) {
      bms_file.read(reinterpret_cast<char*>(&measRecord),
                    sizeof(measurements::measurement_t));
      vbinary_msr->push_back(measRecord);
    }
  } catch (const std::ios_base::failure& f) {
    ss << f.what();
    throw std::runtime_error(ss.str());
  } catch (const std::runtime_error& e) {
    ss << e.what();
    throw std::runtime_error(ss.str());
  } catch (...) {
    throw std::runtime_error(ss.str());
  }

  bms_file.close();
  return bms_meta.binCount;
}

std::optional<std::uint64_t> BmsFileLoader::LoadWithOptional(
    const std::string& bms_filename, measurements::pvmsr_t vbinary_msr,
    binary_file_meta_t& bms_meta) {
  try {
    return LoadFile(bms_filename, vbinary_msr, bms_meta);
  } catch (const std::runtime_error&) {
    return std::nullopt;
  }
}

void BmsFileLoader::WriteFile(const std::string& bms_filename,
                              measurements::pvmsr_t vbinary_msr,
                              binary_file_meta_t& bms_meta) {
  std::ofstream bms_file;
  std::stringstream ss;
  ss << "WriteFile(): An error was encountered when opening " << bms_filename
     << "." << std::endl;

  try {
    file_opener(bms_file, bms_filename, std::ios::out | std::ios::binary,
                binary);
  } catch (const std::runtime_error& e) {
    ss << e.what();
    throw std::runtime_error(ss.str());
  } catch (...) {
    throw std::runtime_error(ss.str());
  }

  ss.str("");
  ss << "WriteFile(): An error was encountered when writing to " << bms_filename
     << "." << std::endl;

  try {
    WriteFileInfo(bms_file);
    WriteFileMetadata(bms_file, bms_meta);

    measurements::it_vmsr_t it_msr(vbinary_msr->begin());
    for (it_msr = vbinary_msr->begin(); it_msr != vbinary_msr->end(); ++it_msr) {
      bms_file.write(reinterpret_cast<char*>(&(*it_msr)),
                     sizeof(measurements::measurement_t));
    }
  } catch (const std::ios_base::failure& f) {
    ss << f.what();
    throw std::runtime_error(ss.str());
  } catch (const std::runtime_error& e) {
    ss << e.what();
    throw std::runtime_error(ss.str());
  } catch (...) {
    throw std::runtime_error(ss.str());
  }

  bms_file.close();
}

void BmsFileLoader::WriteFile(const std::string& bms_filename,
                              measurements::vdnaMsrPtr* vMeasurements,
                              binary_file_meta_t& bms_meta) {
  std::ofstream bms_file;
  std::stringstream ss;
  ss << "WriteFile(): An error was encountered when opening " << bms_filename
     << "." << std::endl;

  try {
    file_opener(bms_file, bms_filename, std::ios::out | std::ios::binary,
                binary);
  } catch (const std::runtime_error& e) {
    ss << e.what();
    throw std::runtime_error(ss.str());
  } catch (...) {
    throw std::runtime_error(ss.str());
  }

  ss.str("");
  ss << "WriteFile(): An error was encountered when writing to " << bms_filename
     << "." << std::endl;

  measurements::_it_vdnamsrptr it_msr;
  UINT32 msrIndex(0);

  try {
    // Build unique source file index map
    std::map<std::string, UINT32> sourceFileMap;
    std::vector<std::string> sourceFileList;
    for (it_msr = vMeasurements->begin(); it_msr != vMeasurements->end(); ++it_msr) {
      std::string src = it_msr->get()->GetSource();
      if (sourceFileMap.find(src) == sourceFileMap.end()) {
        sourceFileMap[src] = static_cast<UINT32>(sourceFileList.size());
        sourceFileList.push_back(src);
      }
      it_msr->get()->SetSourceFileIndex(sourceFileMap[src]);
    }

    // Populate source file metadata
    bms_meta.sourceFileCount = sourceFileList.size();
    if (bms_meta.sourceFileMeta != nullptr)
      delete[] bms_meta.sourceFileMeta;
    bms_meta.sourceFileMeta = new source_file_meta_t[bms_meta.sourceFileCount];
    for (std::uint64_t i(0); i < bms_meta.sourceFileCount; ++i) {
      memset(bms_meta.sourceFileMeta[i].filename, '\0', sizeof(bms_meta.sourceFileMeta[i].filename));
      snprintf(bms_meta.sourceFileMeta[i].filename, sizeof(bms_meta.sourceFileMeta[i].filename),
        "%s", sourceFileList[i].c_str());
    }

    WriteFileInfo(bms_file);
    WriteFileMetadata(bms_file, bms_meta);

    for (it_msr = vMeasurements->begin(); it_msr != vMeasurements->end();
         ++it_msr) {
      it_msr->get()->WriteBinaryMsr(&bms_file, &msrIndex);
    }
  } catch (const std::ios_base::failure& f) {
    ss << f.what();
    throw std::runtime_error(ss.str());
  } catch (const std::runtime_error& e) {
    ss << e.what();
    throw std::runtime_error(ss.str());
  } catch (...) {
    throw std::runtime_error(ss.str());
  }

  bms_file.close();
}

}  // namespace iostreams
}  // namespace dynadjust