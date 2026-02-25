//============================================================================
// Name         : bst_file_loader.cpp
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
// Description  : DynAdjust binary station file io operations
//============================================================================

#include <include/io/bst_file_loader.hpp>

#include <algorithm>
#include <fstream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>

#include <include/functions/dnaiostreamfuncs.hpp>

namespace dynadjust {
namespace iostreams {

namespace {

[[nodiscard]] std::ifstream OpenBstForRead(const std::string& filename) {
  std::ifstream stream{filename, std::ios::binary};
  stream.exceptions(std::ios::badbit | std::ios::failbit);
  if (!stream.is_open()) {
    std::ostringstream os;
    os << "LoadFile(): could not open file '" << filename << "' for reading.";
    throw std::runtime_error(os.str());
  }
  return stream;  // move-return
}

[[nodiscard]] std::ofstream OpenBstForWrite(const std::string& filename) {
  std::ofstream stream{filename, std::ios::binary | std::ios::trunc};
  stream.exceptions(std::ios::badbit | std::ios::failbit);
  if (!stream.is_open()) {
    std::ostringstream os;
    os << "WriteFile(): could not open file '" << filename << "' for writing.";
    throw std::runtime_error(os.str());
  }
  return stream;  // move-return
}

}  // anonymous namespace

BstFileLoader& BstFileLoader::operator=(const BstFileLoader& rhs) {
  if (this != &rhs) {
    DynadjustFile::operator=(rhs);
  }
  return *this;
}

std::uint64_t BstFileLoader::CreateStnInputFileMeta(
    vifm_t& vinput_file_meta, input_file_meta_t** input_file_meta) {
  try {
    // Allocate temporary buffer large enough, copy while counting
    auto temp = std::make_unique<input_file_meta_t[]>(vinput_file_meta.size());

    auto* out = temp.get();
    std::uint64_t cnt = 0;
    for (const auto& ifm : vinput_file_meta) {
      if (ifm.datatype == stn_data || ifm.datatype == stn_msr_data) {
        out[cnt++] = ifm;
      }
    }

    // Free previous buffer (legacy interface)
    delete[] * input_file_meta;
    *input_file_meta = nullptr;

    if (cnt > 0) {
      // Shrink-to-fit: allocate exact-sized array expected by callers
      auto exact = std::make_unique<input_file_meta_t[]>(cnt);
      std::copy_n(temp.get(), cnt, exact.get());
      *input_file_meta = exact.release();
    } else {
      // Legacy interface expects non-null pointer even for empty arrays
      *input_file_meta = new input_file_meta_t[1];  // minimal allocation
    }

    return cnt;

  } catch (const std::exception& e) {
    std::ostringstream os;
    os << "CreateStnInputFileMeta(): " << e.what();
    throw std::runtime_error(os.str());
  }
}

void BstFileLoader::LoadFileMeta(const std::string& bst_filename,
                                 binary_file_meta_t& bst_meta) {
  try {
    auto bst_file = OpenBstForRead(bst_filename);

    // May throw ios::failure
    ReadFileInfo(bst_file);
    ReadFileMetadata(bst_file, bst_meta);

  } catch (const std::exception& e) {
    std::ostringstream os;
    os << "LoadFileMeta(): " << e.what();
    throw std::runtime_error(os.str());
  }
}

std::uint64_t BstFileLoader::LoadFile(const std::string& bst_filename,
                                      pvstn_t vbinary_stn,
                                      binary_file_meta_t& bst_meta) {
  try {
    auto bst_file = OpenBstForRead(bst_filename);

    ReadFileInfo(bst_file);
    ReadFileMetadata(bst_file, bst_meta);

    vbinary_stn->resize(bst_meta.binCount);
    static_assert(std::is_trivially_copyable_v<station_t>,
                  "station_t must be trivially copyable for bulk I/O");

    bst_file.read(reinterpret_cast<char*>(vbinary_stn->data()),
                  sizeof(station_t) * bst_meta.binCount);

    return bst_meta.binCount;

  } catch (const std::exception& e) {
    std::ostringstream os;
    os << "LoadFile(): " << e.what();
    throw std::runtime_error(os.str());
  }
}

std::optional<std::uint64_t> BstFileLoader::LoadWithOptional(
    const std::string& bst_filename, pvstn_t vbinary_stn,
    binary_file_meta_t& bst_meta) {
  try {
    return LoadFile(bst_filename, vbinary_stn, bst_meta);
  } catch (const std::runtime_error&) {
    return std::nullopt;
  }
}

void BstFileLoader::WriteFile(const std::string& bst_filename,
                              pvstn_t vbinary_stn,
                              binary_file_meta_t& bst_meta) {
  try {
    auto bst_file = OpenBstForWrite(bst_filename);

    WriteFileInfo(bst_file);
    WriteFileMetadata(bst_file, bst_meta);

    static_assert(std::is_trivially_copyable_v<station_t>,
                  "station_t must be trivially copyable for bulk I/O");

    bst_file.write(reinterpret_cast<const char*>(vbinary_stn->data()),
                   sizeof(station_t) * vbinary_stn->size());
  } catch (const std::exception& e) {
    std::ostringstream os;
    os << "WriteFile(): " << e.what();
    throw std::runtime_error(os.str());
  }
}

bool BstFileLoader::WriteFile(const std::string& bst_filename,
                              measurements::vdnaStnPtr* vStations,
                              pvstring vUnusedStns,
                              binary_file_meta_t& bst_meta,
                              bool flagUnused) {
  try {
    auto bst_file = OpenBstForWrite(bst_filename);

    WriteFileInfo(bst_file);
    WriteFileMetadata(bst_file, bst_meta);

    if (flagUnused) {
      std::sort(vUnusedStns->begin(), vUnusedStns->end());
    }

    for (auto& stnPtr : *vStations) {
      if (flagUnused &&
          std::binary_search(vUnusedStns->begin(), vUnusedStns->end(),
                             stnPtr->GetName())) {
        stnPtr->SetStationUnused();
        stnPtr->WriteBinaryStn(&bst_file, true);
      } else {
        stnPtr->WriteBinaryStn(&bst_file);
      }
    }
    return true;
  } catch (const std::exception& e) {
    std::ostringstream os;
    os << "WriteFile(): " << e.what();
    throw std::runtime_error(os.str());
  }
}

}  // namespace iostreams
}  // namespace dynadjust