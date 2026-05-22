//============================================================================
// Name         : test_dynadjust_file.cpp
// Author       : Dale Roberts <dale.o.roberts@gmail.com>
// Contributors :
// Copyright    : Copyright 2017-2025 Geoscience Australia
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
// Description  : Unit tests for DynadjustFile (issue #325 back-compat)
//============================================================================

#define TESTING_MAIN
#ifndef __BINARY_NAME__
#define __BINARY_NAME__ "test_dynadjust_file"
#endif
#ifndef __BINARY_DESC__
#define __BINARY_DESC__ "Unit tests for DynadjustFile v1.1 back-compat"
#endif

#include "testing.hpp"

#include "../dynadjust/include/io/dynadjust_file.hpp"
#include <cstdint>
#include <filesystem>
#include <fstream>

using namespace dynadjust::iostreams;

namespace {
const std::string TEMP_V11_FILE = "temp_test_v11.bin";
const std::string TEMP_V12_FILE = "temp_test_v12.bin";

void cleanup(const std::string& path) {
    if (std::filesystem::exists(path))
        std::filesystem::remove(path);
}

void write_identifier_field(std::ofstream& os, const std::string& value) {
    char buf[identifier_field_width + 1];
    snprintf(buf, sizeof(buf), "%*s", identifier_field_width,
             value.substr(0, identifier_field_width).c_str());
    os.write(buf, identifier_field_width);
}

// Writes the DynadjustFile file-info header block with the supplied version
// string (no class setter forces the on-disk version independent of __FILE_VERSION__).
void write_file_info(std::ofstream& os, const std::string& version_value) {
    os.write(version_header, identifier_field_width);
    write_identifier_field(os, version_value);
    os.write(create_date_header, identifier_field_width);
    write_identifier_field(os, "2024-01-15");
    os.write(create_by_header, identifier_field_width);
    write_identifier_field(os, "DNA1.4");
}

// Writes metadata in the pre-v1.2 layout: omits the observation_epoch fields
// from both the top-level meta and each input_file_meta entry.
void write_v11_metadata(std::ofstream& os,
                        const std::string& epsg,
                        const std::string& epoch,
                        std::uint64_t input_file_count,
                        const std::string& input_file_epoch) {
    std::uint64_t bin_count = 0;
    os.write(reinterpret_cast<const char*>(&bin_count), sizeof(bin_count));

    bool reduced = false;
    os.write(reinterpret_cast<const char*>(&reduced), sizeof(reduced));

    char modified_by[MOD_NAME_WIDTH] = {0};
    snprintf(modified_by, sizeof(modified_by), "test_v11");
    os.write(modified_by, MOD_NAME_WIDTH);

    char epsg_buf[STN_EPSG_WIDTH] = {0};
    snprintf(epsg_buf, sizeof(epsg_buf), "%s", epsg.c_str());
    os.write(epsg_buf, STN_EPSG_WIDTH);

    char epoch_buf[STN_EPOCH_WIDTH] = {0};
    snprintf(epoch_buf, sizeof(epoch_buf), "%s", epoch.c_str());
    os.write(epoch_buf, STN_EPOCH_WIDTH);

    // v1.1: no observation_epoch field at this position.

    bool reftran = false;
    os.write(reinterpret_cast<const char*>(&reftran), sizeof(reftran));

    bool geoid = false;
    os.write(reinterpret_cast<const char*>(&geoid), sizeof(geoid));

    os.write(reinterpret_cast<const char*>(&input_file_count), sizeof(input_file_count));
    for (std::uint64_t i = 0; i < input_file_count; ++i) {
        char file_name[FILE_NAME_WIDTH] = {0};
        snprintf(file_name, sizeof(file_name), "input_%llu.xml",
                 static_cast<unsigned long long>(i));
        os.write(file_name, FILE_NAME_WIDTH);

        char input_epsg_buf[STN_EPSG_WIDTH] = {0};
        snprintf(input_epsg_buf, sizeof(input_epsg_buf), "%s", epsg.c_str());
        os.write(input_epsg_buf, STN_EPSG_WIDTH);

        char input_epoch_buf[STN_EPOCH_WIDTH] = {0};
        snprintf(input_epoch_buf, sizeof(input_epoch_buf), "%s", input_file_epoch.c_str());
        os.write(input_epoch_buf, STN_EPOCH_WIDTH);

        // v1.1: no observation_epoch per input file.

        std::uint16_t filetype = 1;
        os.write(reinterpret_cast<const char*>(&filetype), sizeof(filetype));

        std::uint16_t datatype = 1;
        os.write(reinterpret_cast<const char*>(&datatype), sizeof(datatype));
    }

    // sourceFileCount is v1.1+, so it is present (set to zero for simplicity).
    std::uint64_t source_file_count = 0;
    os.write(reinterpret_cast<const char*>(&source_file_count), sizeof(source_file_count));
}
} // namespace

// Regression test for issue #325: v1.1 segment/project files (no
// observation_epoch on disk) must load with observation_epoch defaulted to
// epoch so legacy files remain readable with identical semantics.
TEST_CASE("DynadjustFile v1.1 loads with observation_epoch defaulted to epoch",
          "[DynadjustFile][backcompat]") {
    cleanup(TEMP_V11_FILE);

    const std::string epsg = "7843";
    const std::string epoch = "01.01.2020";
    const std::string input_epoch = "01.01.2020";

    {
        std::ofstream os(TEMP_V11_FILE, std::ios::out | std::ios::binary);
        REQUIRE(os.is_open());
        write_file_info(os, "1.1");
        write_v11_metadata(os, epsg, epoch, 2, input_epoch);
        os.close();
    }

    DynadjustFile file;
    std::ifstream is(TEMP_V11_FILE, std::ios::in | std::ios::binary);
    REQUIRE(is.is_open());

    file.ReadFileInfo(is);
    REQUIRE(file.GetVersion() == "1.1");

    binary_file_meta_t meta;
    file.ReadFileMetadata(is, meta);
    is.close();

    REQUIRE(std::string(meta.epsgCode) == epsg);
    REQUIRE(std::string(meta.epoch) == epoch);
    // Back-compat invariant: observation_epoch falls back to epoch.
    REQUIRE(std::string(meta.observation_epoch) == std::string(meta.epoch));

    REQUIRE(meta.inputFileCount == 2);
    for (std::uint64_t i = 0; i < meta.inputFileCount; ++i) {
        REQUIRE(std::string(meta.inputFileMeta[i].epoch) == input_epoch);
        REQUIRE(std::string(meta.inputFileMeta[i].observation_epoch)
                == std::string(meta.inputFileMeta[i].epoch));
    }

    cleanup(TEMP_V11_FILE);
}

// v1.2 round-trip: writing via DynadjustFile and reading back must preserve a
// distinct observation_epoch independent of reference-frame epoch.
TEST_CASE("DynadjustFile v1.2 round-trip preserves distinct observation_epoch",
          "[DynadjustFile][roundtrip]") {
    cleanup(TEMP_V12_FILE);

    binary_file_meta_t write_meta;
    write_meta.binCount = 0;
    write_meta.reduced = false;
    write_meta.reftran = false;
    write_meta.geoid = false;
    snprintf(write_meta.modifiedBy, sizeof(write_meta.modifiedBy), "test_v12");
    snprintf(write_meta.epsgCode, sizeof(write_meta.epsgCode), "7843");
    snprintf(write_meta.epoch, sizeof(write_meta.epoch), "01.01.2020");
    snprintf(write_meta.observation_epoch, sizeof(write_meta.observation_epoch),
             "15.06.2015");
    write_meta.inputFileCount = 0;
    write_meta.sourceFileCount = 0;

    {
        DynadjustFile file;
        std::ofstream os(TEMP_V12_FILE, std::ios::out | std::ios::binary);
        REQUIRE(os.is_open());
        file.WriteFileInfo(os);
        file.WriteFileMetadata(os, write_meta);
        os.close();
    }

    DynadjustFile file;
    std::ifstream is(TEMP_V12_FILE, std::ios::in | std::ios::binary);
    REQUIRE(is.is_open());

    file.ReadFileInfo(is);
    REQUIRE(file.GetVersion() == "1.2");

    binary_file_meta_t read_meta;
    file.ReadFileMetadata(is, read_meta);
    is.close();

    REQUIRE(std::string(read_meta.epoch) == "01.01.2020");
    REQUIRE(std::string(read_meta.observation_epoch) == "15.06.2015");

    cleanup(TEMP_V12_FILE);
}
