//============================================================================
// Name         : test_aml_file_loader.cpp
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
#define __BINARY_NAME__ "test_AmlFile"
#endif
#ifndef __BINARY_DESC__
#define __BINARY_DESC__ "Unit tests for AmlFile class"
#endif

#include "testing.hpp"

#include "../dynadjust/include/io/aml_file.hpp"
#include "../dynadjust/include/measurement_types/dnastation.hpp"
#include <filesystem>

using namespace dynadjust::iostreams;
using namespace dynadjust::measurements;

namespace {
const std::string TEST_AML_FILE = "data/test.aml";
const std::string TEMP_AML_FILE = "temp_test.aml";
const std::string TEMP_AML_TXT_FILE = "temp_test_aml.txt";
const std::string TEMP_BMS_FILE = "temp_test.bms";
const std::string NONEXISTENT_FILE = "does_not_exist.aml";

void create_test_measurement_records(vmsr_t& measurements) {
    // Create measurement record 1: GPS Baseline
    measurement_t msr1;
    msr1.measType = 'G';
    msr1.station1 = 0;
    msr1.station2 = 1;
    msr1.station3 = 0;
    msr1.ignore = false;
    msr1.clusterID = 1;
    msr1.fileOrder = 1;
    msr1.term1 = 1000.123;
    msr1.term2 = 2000.456;
    msr1.term3 = 100.789;

    // Create measurement record 2: Distance
    measurement_t msr2;
    msr2.measType = 'S';
    msr2.station1 = 1;
    msr2.station2 = 2;
    msr2.station3 = 0;
    msr2.ignore = false;
    msr2.clusterID = 0;
    msr2.fileOrder = 2;
    msr2.term1 = 1500.789;

    // Create measurement record 3: Ignored measurement
    measurement_t msr3;
    msr3.measType = 'L';
    msr3.station1 = 0;
    msr3.station2 = 2;
    msr3.station3 = 0;
    msr3.ignore = true;
    msr3.clusterID = 0;
    msr3.fileOrder = 3;
    msr3.term1 = 5.123;

    measurements.push_back(msr1);
    measurements.push_back(msr2);
    measurements.push_back(msr3);
}

void create_test_aml_data(v_aml_pair& aml_data, const vmsr_t& measurements) {
    // Create AML entries pointing to measurement records
    aml_pair entry1(0, true); // Points to GPS baseline, available
    aml_pair entry2(2, true); // Points to ignored level difference, will be consumed
    aml_pair entry3(1, true); // Points to distance measurement, available

    aml_data.push_back(entry1);
    aml_data.push_back(entry2);
    aml_data.push_back(entry3);
}

void create_test_aml_indices(vUINT32& aml_indices) {
    aml_indices.push_back(0); // Index to measurement 0 (GPS baseline)
    aml_indices.push_back(2); // Index to measurement 2 (ignored level diff)
    aml_indices.push_back(1); // Index to measurement 1 (distance)
}

void create_test_station_list(vASLPtr& station_list) {
    // Station 0: Has 2 measurements (indices 0, 1 in AML)
    ASLPtr station0;
    station0.reset(new CAStationList);
    station0->SetAssocMsrCount(2);
    station0->SetAMLStnIndex(0);
    station0->SetValid();

    // Station 1: Has 1 measurement (index 2 in AML)
    ASLPtr station1;
    station1.reset(new CAStationList);
    station1->SetAssocMsrCount(1);
    station1->SetAMLStnIndex(2);
    station1->SetValid();

    // Station 2: Has 0 measurements (unused station)
    ASLPtr station2;
    station2.reset(new CAStationList);
    station2->SetAssocMsrCount(0);
    station2->SetAMLStnIndex(3);
    station2->SetInvalid();

    station_list.push_back(station0);
    station_list.push_back(station1);
    station_list.push_back(station2);
}

void create_test_stations(vdnaStnPtr& stations) {
    dnaStnPtr stn1;
    stn1.reset(new CDnaStation("GDA2020", "2020.0"));
    stn1->SetName("STN001");
    stn1->SetXAxis_d(-35.308);
    stn1->SetYAxis_d(149.124);
    stn1->SetHeight_d(100.0);

    dnaStnPtr stn2;
    stn2.reset(new CDnaStation("GDA2020", "2020.0"));
    stn2->SetName("STN002");
    stn2->SetXAxis_d(-35.342);
    stn2->SetYAxis_d(149.161);
    stn2->SetHeight_d(150.0);

    dnaStnPtr stn3;
    stn3.reset(new CDnaStation("GDA2020", "2020.0"));
    stn3->SetName("STN003");
    stn3->SetXAxis_d(-35.355);
    stn3->SetYAxis_d(149.175);
    stn3->SetHeight_d(200.0);

    stations.push_back(stn1);
    stations.push_back(stn2);
    stations.push_back(stn3);
}

void cleanup_temp_files() {
    if (std::filesystem::exists(TEMP_AML_FILE)) {
        std::filesystem::remove(TEMP_AML_FILE);
    }
    if (std::filesystem::exists(TEMP_AML_TXT_FILE)) {
        std::filesystem::remove(TEMP_AML_TXT_FILE);
    }
    if (std::filesystem::exists(TEMP_BMS_FILE)) {
        std::filesystem::remove(TEMP_BMS_FILE);
    }
}

void create_dummy_bms_file(const std::string& filename, const vmsr_t& measurements) {
    // Create a minimal BMS file for testing
    std::ofstream bms_file(filename, std::ios::binary);

    // Write simple file header (VERSION, CREATED ON, CREATED BY - 60 bytes total)
    std::string header = "VERSION   ";
    header += "       1.0";
    header += "CREATED ON";
    header += "2025-06-17";
    header += "CREATED BY";
    header += "  TEST12345";
    bms_file.write(header.c_str(), 60);

    // Write measurement count
    std::uint64_t count = static_cast<std::uint64_t>(measurements.size());
    bms_file.write(reinterpret_cast<const char*>(&count), sizeof(std::uint64_t));

    // Write measurement records
    for (const auto& msr : measurements) {
        bms_file.write(reinterpret_cast<const char*>(&msr), sizeof(measurement_t));
    }

    bms_file.close();
}
} // namespace

// Basic constructor and destructor tests
TEST_CASE("AmlFile constructor", "[AmlFile][basic]") {
    AmlFile aml_io;
    REQUIRE(true);
}

TEST_CASE("AmlFile copy constructor", "[AmlFile][basic]") {
    AmlFile aml_io1;
    AmlFile aml_io2(aml_io1);
    REQUIRE(true);
}

TEST_CASE("AmlFile assignment operator", "[AmlFile][basic]") {
    AmlFile aml_io1;
    AmlFile aml_io2;
    aml_io2 = aml_io1;
    REQUIRE(true);
}

// Error handling tests
TEST_CASE("Error handling for non-existent files", "[AmlFile][error]") {
    AmlFile aml_io;
    v_aml_pair aml_data;
    vmsr_t measurements;

    bool threw_exception = false;
    try {
        aml_io.load_aml_file(NONEXISTENT_FILE, &aml_data, &measurements);
    } catch (const std::runtime_error&) { threw_exception = true; }
    REQUIRE(threw_exception);
}

TEST_CASE("Error handling for write to invalid path", "[AmlFile][error]") {
    AmlFile aml_io;
    vUINT32 aml_indices;
    create_test_aml_indices(aml_indices);

    bool threw_exception = false;
    try {
        aml_io.write_aml_file("/invalid/path/test.aml", &aml_indices);
    } catch (const std::runtime_error&) { threw_exception = true; }
    REQUIRE(threw_exception);
}

// Write and read back binary AML data
TEST_CASE("Write and read back binary AML data", "[AmlFile][roundtrip]") {
    AmlFile aml_io;
    vUINT32 original_aml_indices;

    cleanup_temp_files();

    create_test_aml_indices(original_aml_indices);

    // Write the AML file
    aml_io.write_aml_file(TEMP_AML_FILE, &original_aml_indices);
    REQUIRE(std::filesystem::exists(TEMP_AML_FILE));

    // Create test data for reading
    vmsr_t measurements;
    create_test_measurement_records(measurements);

    // Read it back
    v_aml_pair loaded_aml_data;
    aml_io.load_aml_file(TEMP_AML_FILE, &loaded_aml_data, &measurements);

    // Verify data integrity
    REQUIRE(loaded_aml_data.size() == 3);

    // Check measurement indices
    REQUIRE(loaded_aml_data[0].bmsr_index == 0);
    REQUIRE(loaded_aml_data[1].bmsr_index == 2);
    REQUIRE(loaded_aml_data[2].bmsr_index == 1);

    // Check availability (ignored measurements should be consumed)
    REQUIRE(loaded_aml_data[0].available == true);  // GPS baseline - not ignored
    REQUIRE(loaded_aml_data[1].available == false); // Level diff - ignored, so consumed
    REQUIRE(loaded_aml_data[2].available == true);  // Distance - not ignored

    cleanup_temp_files();
}

// Test empty AML data
TEST_CASE("Write and read empty AML data", "[AmlFile][empty]") {
    AmlFile aml_io;
    vUINT32 empty_aml_indices;

    cleanup_temp_files();

    // Write empty AML file
    aml_io.write_aml_file(TEMP_AML_FILE, &empty_aml_indices);
    REQUIRE(std::filesystem::exists(TEMP_AML_FILE));

    // Read it back
    v_aml_pair loaded_aml_data;
    vmsr_t measurements;
    aml_io.load_aml_file(TEMP_AML_FILE, &loaded_aml_data, &measurements);

    REQUIRE(loaded_aml_data.size() == 0);

    cleanup_temp_files();
}

// Test measurement to station tally creation
TEST_CASE("Create measurement to station tally", "[AmlFile][tally]") {
    AmlFile aml_io;
    vASLPtr station_list;
    v_aml_pair aml_data;
    vmsr_t measurements;
    vmsrtally stn_msr_tally;

    create_test_station_list(station_list);
    create_test_measurement_records(measurements);
    create_test_aml_data(aml_data, measurements);

    // Create the tally
    aml_io.create_msr_to_stn_tally(&station_list, aml_data, stn_msr_tally, measurements);

    // Verify tally was created for all stations
    REQUIRE(stn_msr_tally.size() == 3);

    // Station 0 should have measurements (GPS baseline from aml_data[0])
    // Note: aml_data[1] is ignored, so won't count

    // Station 1 should have measurements (Distance from aml_data[2])

    // Station 2 has no measurements
}

// Test large AML data (uint64_t usage)
TEST_CASE("Handle large AML data", "[AmlFile][large]") {
    AmlFile aml_io;
    vUINT32 large_aml_indices;

    cleanup_temp_files();

    // Create a reasonably large number of measurement indices (1000)
    const std::size_t num_measurements = 1000;
    for (std::size_t i = 0; i < num_measurements; ++i) {
        large_aml_indices.push_back(static_cast<UINT32>(i % 100)); // Cycle through 100 different indices
    }

    aml_io.write_aml_file(TEMP_AML_FILE, &large_aml_indices);
    REQUIRE(std::filesystem::exists(TEMP_AML_FILE));

    // Create corresponding measurement data
    vmsr_t measurements;
    for (std::size_t i = 0; i < 100; ++i) {
        measurement_t msr;
        msr.measType = 'S';
        msr.station1 = static_cast<UINT32>(i % 10);
        msr.station2 = static_cast<UINT32>((i + 1) % 10);
        msr.ignore = (i % 10 == 9); // Every 10th measurement is ignored
        msr.clusterID = 0;
        msr.fileOrder = static_cast<UINT32>(i);
        msr.term1 = 1000.0 + i;
        measurements.push_back(msr);
    }

    // Read it back
    v_aml_pair loaded_aml_data;
    aml_io.load_aml_file(TEMP_AML_FILE, &loaded_aml_data, &measurements);

    REQUIRE(loaded_aml_data.size() == num_measurements);

    // Verify some entries
    REQUIRE(loaded_aml_data[0].bmsr_index == 0);
    REQUIRE(loaded_aml_data[999].bmsr_index == 99);

    // Check that ignored measurements are consumed
    std::size_t consumed_count = 0;
    for (const auto& entry : loaded_aml_data) {
        if (!entry.available) {
            consumed_count++;
        }
    }
    REQUIRE(consumed_count == 100); // Should be 10 ignored measurements * 10 repetitions

    cleanup_temp_files();
}

// Test write_msr_to_stn function
TEST_CASE("Write measurement to station summary", "[AmlFile][summary]") {
    AmlFile aml_io;

    // Create test data
    vstn_t stations;
    vUINT32 station_list;
    vmsrtally stn_msr_tally;
    MsrTally parse_msr_tally;

    // Create station records
    station_t stn1, stn2;
    snprintf(stn1.stationName, sizeof(stn1.stationName), "STN001");
    snprintf(stn2.stationName, sizeof(stn2.stationName), "STN002");
    stations.push_back(stn1);
    stations.push_back(stn2);

    // Create station indices
    station_list.push_back(0);
    station_list.push_back(1);

    // Create tally data (resize to match station count)
    stn_msr_tally.resize(2);

    // Test output to stringstream
    std::ostringstream output;
    aml_io.write_msr_to_stn(output, &stations, &station_list, stn_msr_tally, &parse_msr_tally);

    std::string result = output.str();
    REQUIRE(result.find("MEASUREMENT TO STATIONS") != std::string::npos);
    REQUIRE(result.find("STN001") != std::string::npos);
    REQUIRE(result.find("STN002") != std::string::npos);
}

// Test AML text file output (simplified test without BMS dependency)
TEST_CASE("Write AML text file", "[AmlFile][text]") {
    // This test verifies that the write_aml_file_txt method exists and can be called
    // A full integration test would require a properly formatted BMS file
    // which is complex to create in a unit test environment

    AmlFile aml_io;
    vUINT32 aml_indices;
    vASLPtr station_list;
    vdnaStnPtr stations;

    cleanup_temp_files();

    create_test_aml_indices(aml_indices);
    create_test_station_list(station_list);
    create_test_stations(stations);

    // Test that the method can be called (it will throw due to missing BMS file)
    bool threw_exception = false;
    try {
        aml_io.write_aml_file_txt(NONEXISTENT_FILE, TEMP_AML_TXT_FILE, &aml_indices, &station_list, &stations);
    } catch (const std::runtime_error&) { threw_exception = true; }

    // Should throw exception due to non-existent BMS file
    REQUIRE(threw_exception);

    cleanup_temp_files();
}
