//============================================================================
// Name         : test_bst_file.cpp
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
// Description  : Unit tests
//============================================================================

#define TESTING_MAIN
#ifndef __BINARY_NAME__
#define __BINARY_NAME__ "test_BstFile"
#endif
#ifndef __BINARY_DESC__
#define __BINARY_DESC__ "Unit tests for BstFile class"
#endif

#include "testing.hpp"

#include "../dynadjust/include/io/bst_file.hpp"
#include <filesystem>

using namespace dynadjust::iostreams;
using namespace dynadjust::measurements;

namespace {
const std::string TEST_BST_FILE = "data/test.bst";
const std::string TEMP_BST_FILE = "temp_test.bst";
const std::string NONEXISTENT_FILE = "does_not_exist.bst";

void create_test_input_file_meta(vifm_t& vinput_file_meta) {
    input_file_meta_t meta1;
    snprintf(meta1.filename, sizeof(meta1.filename), "station_file1.xml");
    snprintf(meta1.epsgCode, sizeof(meta1.epsgCode), "7843");
    snprintf(meta1.epoch, sizeof(meta1.epoch), "01.01.2020");
    meta1.filetype = 1;
    meta1.datatype = stn_data;

    input_file_meta_t meta2;
    snprintf(meta2.filename, sizeof(meta2.filename), "measurement_file1.xml");
    snprintf(meta2.epsgCode, sizeof(meta2.epsgCode), "7843");
    snprintf(meta2.epoch, sizeof(meta2.epoch), "01.01.2020");
    meta2.filetype = 1;
    meta2.datatype = msr_data;

    input_file_meta_t meta3;
    snprintf(meta3.filename, sizeof(meta3.filename), "combined_file1.xml");
    snprintf(meta3.epsgCode, sizeof(meta3.epsgCode), "7843");
    snprintf(meta3.epoch, sizeof(meta3.epoch), "01.01.2020");
    meta3.filetype = 1;
    meta3.datatype = stn_msr_data;

    vinput_file_meta.push_back(meta1);
    vinput_file_meta.push_back(meta2);
    vinput_file_meta.push_back(meta3);
}

void create_test_station_data(vstn_t& stations) {
    station_t stn1;
    snprintf(stn1.stationName, sizeof(stn1.stationName), "TEST001");
    snprintf(stn1.stationNameOrig, sizeof(stn1.stationNameOrig), "TEST001_20200101");
    // -35.30812534548754, 149.1242571561052
    stn1.initialLatitude = -35.308;
    stn1.currentLatitude = -35.308125345;
    stn1.initialLongitude = 149.124;
    stn1.currentLongitude = 149.124257156;
    stn1.initialHeight = 100.123;
    stn1.currentHeight = 100.123;
    stn1.zone = 55;
    stn1.fileOrder = 1;
    stn1.nameOrder = 1;
    snprintf(stn1.description, sizeof(stn1.description), "Test station 1");

    station_t stn2;
    snprintf(stn2.stationName, sizeof(stn2.stationName), "TEST002");
    snprintf(stn2.stationNameOrig, sizeof(stn2.stationNameOrig), "TEST002_20200101");
    // -35.34254365271435, 149.16111145439325
    stn2.initialLatitude = -35.342;
    stn2.currentLatitude = -35.342543652;
    stn2.initialLongitude = 149.161;
    stn2.currentLongitude = 149.161111454;
    stn2.initialHeight = 150.456;
    stn2.currentHeight = 150.456;
    stn2.zone = 55;
    stn2.fileOrder = 2;
    stn2.nameOrder = 2;
    snprintf(stn2.description, sizeof(stn2.description), "Test station 2");

    stations.push_back(stn1);
    stations.push_back(stn2);
}

void create_test_binary_meta(binary_file_meta_t& meta, UINT32 binCount) {
    meta.binCount = binCount;
    meta.reduced = false;
    meta.reftran = false;
    meta.geoid = false;
    snprintf(meta.modifiedBy, sizeof(meta.modifiedBy), "test_program");
    snprintf(meta.epsgCode, sizeof(meta.epsgCode), "7843");
    snprintf(meta.epoch, sizeof(meta.epoch), "01.01.2020");
    meta.inputFileCount = 0;
    meta.inputFileMeta = nullptr;
}

void cleanup_temp_files() {
    if (std::filesystem::exists(TEMP_BST_FILE)) {
        std::filesystem::remove(TEMP_BST_FILE);
    }
}
} // namespace

// Basic constructor and destructor tests (known to work)
TEST_CASE("BstFile constructor", "[BstFile][basic]") {
    BstFile bst_loader;
    REQUIRE(true);
}

TEST_CASE("BstFile copy constructor", "[BstFile][basic]") {
    BstFile bst_loader1;
    BstFile bst_loader2(bst_loader1);
    REQUIRE(true);
}

TEST_CASE("BstFile assignment operator", "[BstFile][basic]") {
    BstFile bst_loader1;
    BstFile bst_loader2;
    bst_loader2 = bst_loader1;
    REQUIRE(true);
}

TEST_CASE("create_stn_input_file_meta empty input", "[BstFile]") {
    BstFile bst_loader;
    vifm_t vinput_file_meta;
    input_file_meta_t* input_file_meta = nullptr;

    UINT16 count = bst_loader.CreateStnInputFileMeta(vinput_file_meta, &input_file_meta);

    REQUIRE(count == 0);
    REQUIRE(input_file_meta != nullptr);

    delete[] input_file_meta;
}

TEST_CASE("create_stn_input_file_meta with mixed file types", "[BstFile]") {
    BstFile bst_loader;
    vifm_t vinput_file_meta;
    input_file_meta_t* input_file_meta = nullptr;

    create_test_input_file_meta(vinput_file_meta);

    UINT16 count = bst_loader.CreateStnInputFileMeta(vinput_file_meta, &input_file_meta);

    REQUIRE(count == 2); // stn_data and stn_msr_data should be included
    REQUIRE(input_file_meta != nullptr);

    // Verify the first entry (stn_data)
    REQUIRE(input_file_meta[0].datatype == stn_data);
    REQUIRE(std::string(input_file_meta[0].filename) == "station_file1.xml");

    // Verify the second entry (stn_msr_data)
    REQUIRE(input_file_meta[1].datatype == stn_msr_data);
    REQUIRE(std::string(input_file_meta[1].filename) == "combined_file1.xml");

    delete[] input_file_meta;
}

TEST_CASE("Error handling for non-existent files", "[BstFile]") {
    BstFile bst_loader;
    binary_file_meta_t bst_meta;

    bool threw_exception = false;
    try {
        bst_loader.LoadFileMeta(NONEXISTENT_FILE, bst_meta);
    } catch (const std::runtime_error&) { threw_exception = true; }
    REQUIRE(threw_exception);
}

TEST_CASE("Load metadata from test file if available", "[BstFile]") {
    BstFile bst_loader;
    binary_file_meta_t bst_meta;

    if (std::filesystem::exists(TEST_BST_FILE)) {
        bst_loader.LoadFileMeta(TEST_BST_FILE, bst_meta);
        REQUIRE(bst_meta.binCount > 0);
        REQUIRE(std::strlen(bst_meta.epsgCode) > 0);
    } else {
        REQUIRE(true); // Skip if file doesn't exist
    }
}

TEST_CASE("Write and read back synthetic data", "[BstFile]") {
    BstFile bst_loader;
    vstn_t stations;
    binary_file_meta_t bst_meta;

    cleanup_temp_files();

    create_test_station_data(stations);
    create_test_binary_meta(bst_meta, static_cast<UINT32>(stations.size()));

    bst_loader.WriteFile(TEMP_BST_FILE, &stations, bst_meta);
    REQUIRE(std::filesystem::exists(TEMP_BST_FILE));

    vstn_t loaded_stations;
    binary_file_meta_t loaded_meta;
    UINT32 count = bst_loader.LoadFile(TEMP_BST_FILE, &loaded_stations, loaded_meta);

    REQUIRE(count == 2);
    REQUIRE(loaded_stations.size() == 2);
    REQUIRE(std::string(loaded_stations[0].stationName) == "TEST001");
    REQUIRE(std::string(loaded_stations[1].stationName) == "TEST002");

    cleanup_temp_files();
}
