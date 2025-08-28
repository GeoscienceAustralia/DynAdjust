//============================================================================
// Name         : asl_file.hpp
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

#ifndef DYNADJUST_ASL_FILE_H_
#define DYNADJUST_ASL_FILE_H_

#if defined(_MSC_VER)
#if defined(LIST_INCLUDES_ON_BUILD)
#pragma message("  " __FILE__)
#endif
#endif

#include <cstdint>
#include <fstream>
#include <iomanip>
#include <ios>
#include <iostream>
#include <sstream>
#include <filesystem>
#include <optional>

#include <include/io/dynadjust_file.hpp>
#include <include/functions/dnaintegermanipfuncs.hpp>
#include <include/functions/dnaiostreamfuncs.hpp>
#include <include/functions/dnatemplatefuncs.hpp>
#include <include/functions/dnatemplatestnmsrfuncs.hpp>
#include <include/measurement_types/dnastation.hpp>


namespace dynadjust {
namespace iostreams {

// Structure to hold the results of loading an ASL file
struct AslLoadResult {
    measurements::vASL stations;
    vUINT32 free_stations;
    std::uint64_t count;
};

class AslFile : public DynadjustFile {
   public:
    AslFile() = delete;  // No default constructor
    explicit AslFile(const std::filesystem::path& filename);
    AslFile(const AslFile& asl);
    AslFile(AslFile&& asl) noexcept = default;
    virtual ~AslFile() = default;

    AslFile& operator=(const AslFile& rhs);
    AslFile& operator=(AslFile&& rhs) noexcept = default;

    // Modern C++ interface - returns structured data
    AslLoadResult Load();
    std::optional<AslLoadResult> TryLoad();

    // Legacy interface for backward compatibility
    std::uint64_t LoadLegacy(measurements::vASL* binary_asl, vUINT32* free_stn);

    void Write(const measurements::vASLPtr& binary_asl);

    void WriteText(const measurements::vASLPtr& binary_asl, const measurements::vdnaStnPtr& stations);

    const std::filesystem::path& GetPath() const { return path_; }
    bool Exists() const { return std::filesystem::exists(path_); }

   protected:
    std::filesystem::path path_;
};

}  // namespace iostreams
}  // namespace dynadjust

#endif  // DYNADJUST_ASL_FILE_H_
