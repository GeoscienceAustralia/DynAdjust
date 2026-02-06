//============================================================================
// Name         : test_bms_file.cpp
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
#define __BINARY_NAME__ "test_BmsFile"
#define __BINARY_DESC__ "Unit tests for BmsFile class"

#include "testing.hpp"

#ifndef __BINARY_NAME__
#define __BINARY_NAME__ "test_BmsFile"
#endif
#ifndef __BINARY_DESC__
#define __BINARY_DESC__ "Unit tests for BmsFile class"
#endif

#include "../dynadjust/include/io/bms_file.hpp"
#include "../dynadjust/include/measurement_types/dnameasurement.hpp"
#include <filesystem>

using namespace dynadjust::iostreams;
using namespace dynadjust::measurements;

namespace {
const std::string TEST_BMS_FILE = "data/test.bms";
const std::string TEMP_BMS_FILE = "temp_test.bms";
const std::string NONEXISTENT_FILE = "does_not_exist.bms";

void create_test_measurement_records(vmsr_t& measurements) {
    // Create measurement record 1: GPS Baseline
    measurement_t msr1 = {};
    msr1.measType = 'G';
    msr1.measStart = 0;
    msr1.measurementStations = 2;
    snprintf(msr1.epsgCode, sizeof(msr1.epsgCode), "7843");
    snprintf(msr1.epoch, sizeof(msr1.epoch), "01.01.2020");
    snprintf(msr1.coordType, sizeof(msr1.coordType), "LLH");
    msr1.ignore = false;
    msr1.station1 = 0;
    msr1.station2 = 1;
    msr1.station3 = 0;
    msr1.vectorCount1 = 1;
    msr1.vectorCount2 = 0;
    msr1.clusterID = 1;
    msr1.fileOrder = 1;
    msr1.term1 = 1000.123;
    msr1.term2 = 2000.456;
    msr1.term3 = 100.789;
    msr1.scale1 = 1.0;
    msr1.scale2 = 1.0;
    msr1.scale3 = 1.0;

    // Create measurement record 2: Distance
    measurement_t msr2 = {};
    msr2.measType = 'S';
    msr2.measStart = 0;
    msr2.measurementStations = 2;
    snprintf(msr2.epsgCode, sizeof(msr2.epsgCode), "7843");
    snprintf(msr2.epoch, sizeof(msr2.epoch), "01.01.2020");
    snprintf(msr2.coordType, sizeof(msr2.coordType), "LLH");
    msr2.ignore = false;
    msr2.station1 = 1;
    msr2.station2 = 2;
    msr2.station3 = 0;
    msr2.vectorCount1 = 0;
    msr2.vectorCount2 = 0;
    msr2.clusterID = 0;
    msr2.fileOrder = 2;
    msr2.term1 = 1500.789;
    msr2.scale1 = 1.0;

    // Create measurement record 3: Ignored level difference
    measurement_t msr3 = {};
    msr3.measType = 'L';
    msr3.measStart = 0;
    msr3.measurementStations = 2;
    snprintf(msr3.epsgCode, sizeof(msr3.epsgCode), "7843");
    snprintf(msr3.epoch, sizeof(msr3.epoch), "01.01.2020");
    snprintf(msr3.coordType, sizeof(msr3.coordType), "LLH");
    msr3.ignore = true;
    msr3.station1 = 0;
    msr3.station2 = 2;
    msr3.station3 = 0;
    msr3.vectorCount1 = 0;
    msr3.vectorCount2 = 0;
    msr3.clusterID = 0;
    msr3.fileOrder = 3;
    msr3.term1 = 5.123;
    msr3.scale1 = 1.0;

    measurements.push_back(msr1);
    measurements.push_back(msr2);
    measurements.push_back(msr3);
}

void create_test_binary_meta(binary_file_meta_t& meta, std::uint64_t binCount) {
    meta.binCount = binCount;
    meta.reduced = false;
    meta.reftran = false;
    meta.geoid = false;
    snprintf(meta.modifiedBy, sizeof(meta.modifiedBy), "test_program");
    snprintf(meta.epsgCode, sizeof(meta.epsgCode), "7843");
    snprintf(meta.epoch, sizeof(meta.epoch), "01.01.2020");
    meta.inputFileCount = 0;
    meta.inputFileMeta = nullptr;
    meta.sourceFileCount = 0;
    meta.sourceFileMeta = nullptr;
}

void create_test_input_file_meta(vifm_t& input_files) {
    input_file_meta_t meta1;
    snprintf(meta1.filename, sizeof(meta1.filename), "measurement_file1.xml");
    snprintf(meta1.epsgCode, sizeof(meta1.epsgCode), "7843");
    snprintf(meta1.epoch, sizeof(meta1.epoch), "01.01.2020");
    meta1.filetype = 1;
    meta1.datatype = msr_data;

    input_file_meta_t meta2;
    snprintf(meta2.filename, sizeof(meta2.filename), "combined_file1.xml");
    snprintf(meta2.epsgCode, sizeof(meta2.epsgCode), "7843");
    snprintf(meta2.epoch, sizeof(meta2.epoch), "01.01.2020");
    meta2.filetype = 1;
    meta2.datatype = stn_msr_data;

    input_file_meta_t meta3;
    snprintf(meta3.filename, sizeof(meta3.filename), "station_file1.xml");
    snprintf(meta3.epsgCode, sizeof(meta3.epsgCode), "7843");
    snprintf(meta3.epoch, sizeof(meta3.epoch), "01.01.2020");
    meta3.filetype = 1;
    meta3.datatype = stn_data;

    input_files.push_back(meta1);
    input_files.push_back(meta2);
    input_files.push_back(meta3);
}

void cleanup_temp_files() {
    if (std::filesystem::exists(TEMP_BMS_FILE)) {
        std::filesystem::remove(TEMP_BMS_FILE);
    }
}
} // namespace

// Basic constructor and destructor tests
TEST_CASE("BmsFile constructor", "[BmsFile][basic]") {
    BmsFile bms_loader;
    REQUIRE(true);
}

TEST_CASE("BmsFile copy constructor", "[BmsFile][basic]") {
    BmsFile bms_loader1;
    BmsFile bms_loader2(bms_loader1);
    REQUIRE(true);
}

TEST_CASE("BmsFile assignment operator", "[BmsFile][basic]") {
    BmsFile bms_loader1;
    BmsFile bms_loader2;
    bms_loader2 = bms_loader1;
    REQUIRE(true);
}

// Test create_msr_input_file_meta function
TEST_CASE("create_msr_input_file_meta empty input", "[BmsFile]") {
    BmsFile bms_loader;
    vifm_t input_files;
    input_file_meta_t* file_meta = nullptr;

    std::uint64_t count = bms_loader.CreateMsrInputFileMeta(input_files, &file_meta);

    REQUIRE(count == 0);
    REQUIRE(file_meta != nullptr);

    delete[] file_meta;
}

TEST_CASE("create_msr_input_file_meta with mixed file types", "[BmsFile]") {
    BmsFile bms_loader;
    vifm_t input_files;
    input_file_meta_t* file_meta = nullptr;

    create_test_input_file_meta(input_files);

    std::uint64_t count = bms_loader.CreateMsrInputFileMeta(input_files, &file_meta);

    REQUIRE(count == 2); // msr_data and stn_msr_data should be included
    REQUIRE(file_meta != nullptr);

    // Verify the first entry (msr_data)
    REQUIRE(file_meta[0].datatype == msr_data);
    REQUIRE(std::string(file_meta[0].filename) == "measurement_file1.xml");

    // Verify the second entry (stn_msr_data)
    REQUIRE(file_meta[1].datatype == stn_msr_data);
    REQUIRE(std::string(file_meta[1].filename) == "combined_file1.xml");

    delete[] file_meta;
}

// Error handling tests
TEST_CASE("Error handling for non-existent files", "[BmsFile][error]") {
    BmsFile bms_loader;
    binary_file_meta_t bms_meta;

    bool threw_exception = false;
    try {
        bms_loader.LoadFileMeta(NONEXISTENT_FILE, bms_meta);
    } catch (const std::runtime_error&) { threw_exception = true; }
    REQUIRE(threw_exception);
}

TEST_CASE("Error handling for write to invalid path", "[BmsFile][error]") {
    BmsFile bms_loader;
    vmsr_t measurements;
    binary_file_meta_t bms_meta;

    create_test_measurement_records(measurements);
    create_test_binary_meta(bms_meta, static_cast<std::uint64_t>(measurements.size()));

    bool threw_exception = false;
    try {
        bms_loader.WriteFile("/invalid/path/test.bms", &measurements, bms_meta);
    } catch (const std::runtime_error&) { threw_exception = true; }
    REQUIRE(threw_exception);
}

// Load metadata from test file if available
TEST_CASE("Load metadata from test file if available", "[BmsFile][file]") {
    BmsFile bms_loader;
    binary_file_meta_t bms_meta;

    if (std::filesystem::exists(TEST_BMS_FILE)) {
        bms_loader.LoadFileMeta(TEST_BMS_FILE, bms_meta);
        REQUIRE(bms_meta.binCount > 0);
        REQUIRE(std::strlen(bms_meta.epsgCode) > 0);
    } else {
        REQUIRE(true); // Skip if file doesn't exist
    }
}

// Load full BMS file if available
TEST_CASE("Load BMS file if available", "[BmsFile][file]") {
    BmsFile bms_loader;
    vmsr_t measurements;
    binary_file_meta_t bms_meta;

    if (std::filesystem::exists(TEST_BMS_FILE)) {
        std::uint64_t count = bms_loader.LoadFile(TEST_BMS_FILE, &measurements, bms_meta);
        REQUIRE(count > 0);
        REQUIRE(measurements.size() == count);
        REQUIRE(measurements.size() == bms_meta.binCount);

        // Verify some basic measurement data
        for (const auto& msr : measurements) {
            REQUIRE(msr.measType != 0);        // Should have a measurement type
            REQUIRE(msr.station1 != UINT_MAX); // Should have valid station indices
        }
    } else {
        REQUIRE(true); // Skip if file doesn't exist
    }
}

// Write and read back synthetic BMS data
TEST_CASE("Write and read back synthetic BMS data", "[BmsFile][roundtrip]") {
    BmsFile bms_loader;
    vmsr_t original_measurements;
    binary_file_meta_t original_meta;

    cleanup_temp_files();

    create_test_measurement_records(original_measurements);
    create_test_binary_meta(original_meta, static_cast<std::uint64_t>(original_measurements.size()));

    // Write the BMS file
    bms_loader.WriteFile(TEMP_BMS_FILE, &original_measurements, original_meta);
    REQUIRE(std::filesystem::exists(TEMP_BMS_FILE));

    // Read it back
    vmsr_t loaded_measurements;
    binary_file_meta_t loaded_meta;
    std::uint64_t count = bms_loader.LoadFile(TEMP_BMS_FILE, &loaded_measurements, loaded_meta);

    // Verify data integrity
    REQUIRE(count == 3);
    REQUIRE(loaded_measurements.size() == 3);
    REQUIRE(loaded_meta.binCount == 3);

    // Verify metadata
    REQUIRE(std::string(loaded_meta.epsgCode) == "7843");
    REQUIRE(std::string(loaded_meta.epoch) == "01.01.2020");
    REQUIRE(std::string(loaded_meta.modifiedBy) == "test_program");
    REQUIRE(loaded_meta.reduced == false);
    REQUIRE(loaded_meta.reftran == false);
    REQUIRE(loaded_meta.geoid == false);

    // Verify measurement data
    REQUIRE(loaded_measurements[0].measType == 'G');
    REQUIRE(loaded_measurements[0].station1 == 0);
    REQUIRE(loaded_measurements[0].station2 == 1);
    REQUIRE(loaded_measurements[0].ignore == false);
    REQUIRE(loaded_measurements[0].term1 == 1000.123);

    REQUIRE(loaded_measurements[1].measType == 'S');
    REQUIRE(loaded_measurements[1].station1 == 1);
    REQUIRE(loaded_measurements[1].station2 == 2);
    REQUIRE(loaded_measurements[1].ignore == false);
    REQUIRE(loaded_measurements[1].term1 == 1500.789);

    REQUIRE(loaded_measurements[2].measType == 'L');
    REQUIRE(loaded_measurements[2].station1 == 0);
    REQUIRE(loaded_measurements[2].station2 == 2);
    REQUIRE(loaded_measurements[2].ignore == true);
    REQUIRE(loaded_measurements[2].term1 == 5.123);

    cleanup_temp_files();
}

// Test empty BMS data
TEST_CASE("Write and read empty BMS data", "[BmsFile][empty]") {
    BmsFile bms_loader;
    vmsr_t empty_measurements;
    binary_file_meta_t empty_meta;

    cleanup_temp_files();

    create_test_binary_meta(empty_meta, 0);

    // Write empty BMS file
    bms_loader.WriteFile(TEMP_BMS_FILE, &empty_measurements, empty_meta);
    REQUIRE(std::filesystem::exists(TEMP_BMS_FILE));

    // Read it back
    vmsr_t loaded_measurements;
    binary_file_meta_t loaded_meta;
    std::uint64_t count = bms_loader.LoadFile(TEMP_BMS_FILE, &loaded_measurements, loaded_meta);

    REQUIRE(count == 0);
    REQUIRE(loaded_measurements.size() == 0);
    REQUIRE(loaded_meta.binCount == 0);

    cleanup_temp_files();
}

// Test large BMS data (uint64_t usage)
TEST_CASE("Handle large BMS data", "[BmsFile][large]") {
    BmsFile bms_loader;
    vmsr_t large_measurements;
    binary_file_meta_t large_meta;

    cleanup_temp_files();

    // Create a reasonably large number of measurements (1000)
    const std::size_t num_measurements = 1000;
    large_measurements.reserve(num_measurements);

    for (std::size_t i = 0; i < num_measurements; ++i) {
        measurement_t msr = {};
        msr.measType = 'S';
        msr.station1 = static_cast<UINT32>(i % 100);
        msr.station2 = static_cast<UINT32>((i + 1) % 100);
        msr.ignore = (i % 50 == 0); // Every 50th measurement is ignored
        msr.fileOrder = static_cast<UINT32>(i);
        msr.term1 = 1000.0 + static_cast<double>(i);
        snprintf(msr.epsgCode, sizeof(msr.epsgCode), "7843");
        snprintf(msr.epoch, sizeof(msr.epoch), "01.01.2020");
        large_measurements.push_back(msr);
    }

    create_test_binary_meta(large_meta, static_cast<std::uint64_t>(large_measurements.size()));

    bms_loader.WriteFile(TEMP_BMS_FILE, &large_measurements, large_meta);
    REQUIRE(std::filesystem::exists(TEMP_BMS_FILE));

    // Read it back
    vmsr_t loaded_measurements;
    binary_file_meta_t loaded_meta;
    std::uint64_t count = bms_loader.LoadFile(TEMP_BMS_FILE, &loaded_measurements, loaded_meta);

    REQUIRE(count == num_measurements);
    REQUIRE(loaded_measurements.size() == num_measurements);
    REQUIRE(loaded_meta.binCount == num_measurements);

    // Verify some entries
    REQUIRE(loaded_measurements[0].term1 == 1000.0);
    REQUIRE(loaded_measurements[999].term1 == 1999.0);
    REQUIRE(loaded_measurements[0].ignore == true);  // i=0, 0%50==0
    REQUIRE(loaded_measurements[1].ignore == false); // i=1, 1%50!=0
    REQUIRE(loaded_measurements[50].ignore == true); // i=50, 50%50==0

    cleanup_temp_files();
}

// Test metadata-only operations
TEST_CASE("Load metadata only", "[BmsFile][metadata]") {
    BmsFile bms_loader;
    vmsr_t measurements;
    binary_file_meta_t write_meta;

    cleanup_temp_files();

    create_test_measurement_records(measurements);
    create_test_binary_meta(write_meta, static_cast<std::uint64_t>(measurements.size()));

    // Write a BMS file
    bms_loader.WriteFile(TEMP_BMS_FILE, &measurements, write_meta);
    REQUIRE(std::filesystem::exists(TEMP_BMS_FILE));

    // Load only metadata
    binary_file_meta_t read_meta;
    bms_loader.LoadFileMeta(TEMP_BMS_FILE, read_meta);

    // Verify metadata matches
    REQUIRE(read_meta.binCount == write_meta.binCount);
    REQUIRE(std::string(read_meta.epsgCode) == std::string(write_meta.epsgCode));
    REQUIRE(std::string(read_meta.epoch) == std::string(write_meta.epoch));
    REQUIRE(std::string(read_meta.modifiedBy) == std::string(write_meta.modifiedBy));
    REQUIRE(read_meta.reduced == write_meta.reduced);
    REQUIRE(read_meta.reftran == write_meta.reftran);
    REQUIRE(read_meta.geoid == write_meta.geoid);

    cleanup_temp_files();
}

// Test measurement type variety
TEST_CASE("Handle various measurement types", "[BmsFile][types]") {
    BmsFile bms_loader;
    vmsr_t measurements;
    binary_file_meta_t meta;

    cleanup_temp_files();

    // Create measurements of different types
    const char types[] = {'A', 'B', 'C', 'D', 'E', 'G', 'H', 'I', 'J', 'K',
                          'L', 'M', 'P', 'Q', 'R', 'S', 'V', 'X', 'Y', 'Z'};
    const std::size_t num_types = sizeof(types) / sizeof(types[0]);

    for (std::size_t i = 0; i < num_types; ++i) {
        measurement_t msr = {};
        msr.measType = types[i];
        msr.station1 = static_cast<UINT32>(i);
        msr.station2 = static_cast<UINT32>((i + 1) % num_types);
        msr.ignore = false;
        msr.fileOrder = static_cast<UINT32>(i);
        msr.term1 = 1000.0 + static_cast<double>(i);
        snprintf(msr.epsgCode, sizeof(msr.epsgCode), "7843");
        snprintf(msr.epoch, sizeof(msr.epoch), "01.01.2020");
        measurements.push_back(msr);
    }

    create_test_binary_meta(meta, static_cast<std::uint64_t>(measurements.size()));

    // Write and read back
    bms_loader.WriteFile(TEMP_BMS_FILE, &measurements, meta);

    vmsr_t loaded_measurements;
    binary_file_meta_t loaded_meta;
    std::uint64_t count = bms_loader.LoadFile(TEMP_BMS_FILE, &loaded_measurements, loaded_meta);

    REQUIRE(count == num_types);
    REQUIRE(loaded_measurements.size() == num_types);

    // Verify all measurement types are preserved
    for (std::size_t i = 0; i < num_types; ++i) {
        REQUIRE(loaded_measurements[i].measType == types[i]);
        REQUIRE(loaded_measurements[i].station1 == i);
        REQUIRE(loaded_measurements[i].term1 == 1000.0 + static_cast<double>(i));
    }

    cleanup_temp_files();
}

TEST_CASE("Source file index round-trip", "[BmsFile][roundtrip][source]") {
    BmsFile bms_loader;
    vmsr_t measurements;
    binary_file_meta_t meta;

    cleanup_temp_files();

    // Create measurements with different source file indices
    measurement_t msr1 = {};
    msr1.measType = 'S';
    msr1.measurementStations = 2;
    snprintf(msr1.epsgCode, sizeof(msr1.epsgCode), "7843");
    snprintf(msr1.epoch, sizeof(msr1.epoch), "01.01.2020");
    msr1.station1 = 0;
    msr1.station2 = 1;
    msr1.fileOrder = 0;
    msr1.term1 = 100.0;
    msr1.sourceFileIndex = 0;

    measurement_t msr2 = {};
    msr2.measType = 'L';
    msr2.measurementStations = 2;
    snprintf(msr2.epsgCode, sizeof(msr2.epsgCode), "7843");
    snprintf(msr2.epoch, sizeof(msr2.epoch), "01.01.2020");
    msr2.station1 = 1;
    msr2.station2 = 2;
    msr2.fileOrder = 1;
    msr2.term1 = 5.0;
    msr2.sourceFileIndex = 1;

    measurement_t msr3 = {};
    msr3.measType = 'H';
    msr3.measurementStations = 2;
    snprintf(msr3.epsgCode, sizeof(msr3.epsgCode), "7843");
    snprintf(msr3.epoch, sizeof(msr3.epoch), "01.01.2020");
    msr3.station1 = 2;
    msr3.station2 = 0;
    msr3.fileOrder = 2;
    msr3.term1 = 50.0;
    msr3.sourceFileIndex = 0;

    measurements.push_back(msr1);
    measurements.push_back(msr2);
    measurements.push_back(msr3);

    create_test_binary_meta(meta, measurements.size());

    meta.sourceFileCount = 2;
    meta.sourceFileMeta = new source_file_meta_t[2]();
    snprintf(meta.sourceFileMeta[0].filename, sizeof(meta.sourceFileMeta[0].filename), "network_a.xml");
    snprintf(meta.sourceFileMeta[1].filename, sizeof(meta.sourceFileMeta[1].filename), "network_b.xml");

    bms_loader.WriteFile(TEMP_BMS_FILE, &measurements, meta);
    REQUIRE(std::filesystem::exists(TEMP_BMS_FILE));

    vmsr_t loaded_measurements;
    binary_file_meta_t loaded_meta;
    std::uint64_t count = bms_loader.LoadFile(TEMP_BMS_FILE, &loaded_measurements, loaded_meta);

    REQUIRE(count == 3);
    REQUIRE(loaded_measurements[0].sourceFileIndex == 0);
    REQUIRE(loaded_measurements[1].sourceFileIndex == 1);
    REQUIRE(loaded_measurements[2].sourceFileIndex == 0);

    REQUIRE(loaded_meta.sourceFileCount == 2);
    REQUIRE(loaded_meta.sourceFileMeta != nullptr);
    REQUIRE(std::string(loaded_meta.sourceFileMeta[0].filename) == "network_a.xml");
    REQUIRE(std::string(loaded_meta.sourceFileMeta[1].filename) == "network_b.xml");

    cleanup_temp_files();
}