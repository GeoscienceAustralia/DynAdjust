//============================================================================
// Name         : test_asl_file.cpp
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
#define __BINARY_NAME__ "test_AslFile"
#endif
#ifndef __BINARY_DESC__
#define __BINARY_DESC__ "Unit tests for AslFile class"
#endif

#include "testing.hpp"

#include "../dynadjust/include/io/asl_file.hpp"
#include "../dynadjust/include/measurement_types/dnastation.hpp"
#include <filesystem>

using namespace dynadjust::iostreams;
using namespace dynadjust::measurements;

namespace {
const std::string TEST_ASL_FILE = "data/test.asl";
const std::string TEMP_ASL_FILE = "temp_test.asl";
const std::string TEMP_ASL_TXT_FILE = "temp_test.txt";
const std::string NONEXISTENT_FILE = "does_not_exist.asl";

void create_test_asl_data(vASL& asl_data) {
    // Since CAStationList copy constructor is private, resize and use references
    asl_data.resize(3);

    asl_data[0].SetAssocMsrCount(5);
    asl_data[0].SetAMLStnIndex(0);
    asl_data[0].SetValid();

    asl_data[1].SetAssocMsrCount(3);
    asl_data[1].SetAMLStnIndex(1);
    asl_data[1].SetValid();

    asl_data[2].SetAssocMsrCount(0);
    asl_data[2].SetAMLStnIndex(2);
    asl_data[2].SetInvalid();
}

void create_test_asl_ptr_data(pvASLPtr vbinary_asl) {
    ASLPtr asl1;
    asl1.reset(new CAStationList);
    asl1->SetAssocMsrCount(5);
    asl1->SetAMLStnIndex(0);
    asl1->SetValid();

    ASLPtr asl2;
    asl2.reset(new CAStationList);
    asl2->SetAssocMsrCount(3);
    asl2->SetAMLStnIndex(1);
    asl2->SetValid();

    ASLPtr asl3;
    asl3.reset(new CAStationList);
    asl3->SetAssocMsrCount(0);
    asl3->SetAMLStnIndex(2);
    asl3->SetInvalid();

    vbinary_asl->push_back(asl1);
    vbinary_asl->push_back(asl2);
    vbinary_asl->push_back(asl3);
}

void create_test_station_data(vdnaStnPtr& stations) {
    dnaStnPtr stn1;
    stn1.reset(new CDnaStation("GDA2020", "2020.0"));
    stn1->SetName("TEST001");
    stn1->SetXAxis_d(-35.308125345);
    stn1->SetYAxis_d(149.124257156);
    stn1->SetHeight_d(100.123);

    dnaStnPtr stn2;
    stn2.reset(new CDnaStation("GDA2020", "2020.0"));
    stn2->SetName("TEST002");
    stn2->SetXAxis_d(-35.342543652);
    stn2->SetYAxis_d(149.161111454);
    stn2->SetHeight_d(150.456);

    dnaStnPtr stn3;
    stn3.reset(new CDnaStation("GDA2020", "2020.0"));
    stn3->SetName("TEST003");
    stn3->SetXAxis_d(-35.355);
    stn3->SetYAxis_d(149.175);
    stn3->SetHeight_d(200.789);

    stations.push_back(stn1);
    stations.push_back(stn2);
    stations.push_back(stn3);
}

void cleanup_temp_files() {
    if (std::filesystem::exists(TEMP_ASL_FILE)) {
        std::filesystem::remove(TEMP_ASL_FILE);
    }
    if (std::filesystem::exists(TEMP_ASL_TXT_FILE)) {
        std::filesystem::remove(TEMP_ASL_TXT_FILE);
    }
}
} // namespace

// Basic constructor and destructor tests
TEST_CASE("AslFile constructor", "[asl_file][basic]") {
    AslFile asl_loader("test.asl");
    REQUIRE(asl_loader.GetPath() == "test.asl");
}

TEST_CASE("AslFile copy constructor", "[asl_file][basic]") {
    AslFile asl_loader1("test.asl");
    AslFile asl_loader2(asl_loader1);
    REQUIRE(asl_loader2.GetPath() == "test.asl");
}

TEST_CASE("AslFile assignment operator", "[asl_file][basic]") {
    AslFile asl_loader1("test1.asl");
    AslFile asl_loader2("test2.asl");
    asl_loader2 = asl_loader1;
    REQUIRE(asl_loader2.GetPath() == "test1.asl");
}

// Error handling tests
TEST_CASE("Error handling for non-existent files", "[asl_file][error]") {
    AslFile asl_loader(NONEXISTENT_FILE);

    auto result = asl_loader.TryLoad();
    REQUIRE(result == std::nullopt);
}

TEST_CASE("Error handling for write to invalid path", "[asl_file][error]") {
    AslFile asl_loader("/invalid/path/test.asl");
    vASLPtr asl_ptr_data;
    create_test_asl_ptr_data(&asl_ptr_data);

    bool threw_exception = false;
    try {
        asl_loader.Write(asl_ptr_data);
    } catch (const std::runtime_error&) { threw_exception = true; }
    REQUIRE(threw_exception);
}

// Load test file if available
TEST_CASE("Load ASL from test file if available", "[asl_file][file]") {
    if (std::filesystem::exists(TEST_ASL_FILE)) {
        AslFile asl_loader(TEST_ASL_FILE);

        auto result = asl_loader.Load();
        REQUIRE(result.count > 0);
        REQUIRE(result.stations.size() == result.count);
        REQUIRE(result.free_stations.size() == result.count);

        // Verify free stations are initialized as incrementing vector
        for (std::size_t i = 0; i < result.free_stations.size(); ++i) {
            REQUIRE(result.free_stations[i] == i);
        }
    } else {
        REQUIRE(true); // Skip if file doesn't exist
    }
}

// Write and read back binary ASL data
TEST_CASE("Write and read back binary ASL data", "[asl_file][roundtrip]") {
    cleanup_temp_files();

    vASLPtr original_asl_data;
    create_test_asl_ptr_data(&original_asl_data);

    // Write the ASL file
    AslFile asl_writer(TEMP_ASL_FILE);
    asl_writer.Write(original_asl_data);
    REQUIRE(std::filesystem::exists(TEMP_ASL_FILE));

    // Read it back
    AslFile asl_reader(TEMP_ASL_FILE);
    auto loaded_result = asl_reader.Load();

    // Verify data integrity
    REQUIRE(loaded_result.count == 3);
    REQUIRE(loaded_result.stations.size() == 3);
    REQUIRE(loaded_result.free_stations.size() == 3);

    // Verify station data
    REQUIRE(loaded_result.stations[0].GetAssocMsrCount() == 5);
    REQUIRE(loaded_result.stations[0].GetAMLStnIndex() == 0);
    REQUIRE(loaded_result.stations[0].IsValid() == true);

    REQUIRE(loaded_result.stations[1].GetAssocMsrCount() == 3);
    REQUIRE(loaded_result.stations[1].GetAMLStnIndex() == 1);
    REQUIRE(loaded_result.stations[1].IsValid() == true);

    REQUIRE(loaded_result.stations[2].GetAssocMsrCount() == 0);
    REQUIRE(loaded_result.stations[2].GetAMLStnIndex() == 2);
    REQUIRE(loaded_result.stations[2].IsValid() == false);

    // Verify free stations list is correctly initialized
    for (std::size_t i = 0; i < loaded_result.free_stations.size(); ++i) {
        REQUIRE(loaded_result.free_stations[i] == i);
    }

    cleanup_temp_files();
}

// Test empty ASL data
TEST_CASE("Write and read empty ASL data", "[asl_file][empty]") {
    cleanup_temp_files();

    vASLPtr empty_asl_data;

    // Write empty ASL file
    AslFile asl_writer(TEMP_ASL_FILE);
    asl_writer.Write(empty_asl_data);
    REQUIRE(std::filesystem::exists(TEMP_ASL_FILE));

    // Read it back
    AslFile asl_reader(TEMP_ASL_FILE);
    auto loaded_result = asl_reader.Load();

    REQUIRE(loaded_result.count == 0);
    REQUIRE(loaded_result.stations.size() == 0);
    REQUIRE(loaded_result.free_stations.size() == 0);

    cleanup_temp_files();
}

// Test ASL text output
TEST_CASE("Write ASL text file", "[asl_file][text]") {
    cleanup_temp_files();

    vASLPtr asl_ptr_data;
    vdnaStnPtr stations;
    create_test_asl_ptr_data(&asl_ptr_data);
    create_test_station_data(stations);

    // Write text ASL file
    AslFile asl_writer(TEMP_ASL_TXT_FILE);
    asl_writer.WriteText(asl_ptr_data, stations);
    REQUIRE(std::filesystem::exists(TEMP_ASL_TXT_FILE));

    // Verify file has content
    std::filesystem::file_size(TEMP_ASL_TXT_FILE);
    REQUIRE(std::filesystem::file_size(TEMP_ASL_TXT_FILE) > 0);

    // Read and verify basic content structure
    std::ifstream txt_file(TEMP_ASL_TXT_FILE);
    std::string first_line;
    std::getline(txt_file, first_line);
    REQUIRE(first_line.find("3 stations") != std::string::npos);

    cleanup_temp_files();
}

// Test with null pointers in ASL data
TEST_CASE("Handle null pointers in ASL data", "[asl_file][nulls]") {
    cleanup_temp_files();

    vASLPtr asl_ptr_data;

    // Create ASL data with some null entries
    ASLPtr asl1;
    asl1.reset(new CAStationList);
    asl1->SetAssocMsrCount(5);
    asl1->SetAMLStnIndex(0);
    asl1->SetValid();

    asl_ptr_data.push_back(asl1);
    asl_ptr_data.push_back(nullptr); // null entry
    asl_ptr_data.push_back(asl1);    // valid entry

    // Should write only non-null entries
    AslFile asl_writer(TEMP_ASL_FILE);
    asl_writer.Write(asl_ptr_data);
    REQUIRE(std::filesystem::exists(TEMP_ASL_FILE));

    cleanup_temp_files();
}

// Test large station count (uint64_t usage)
TEST_CASE("Handle large station counts", "[asl_file][large]") {
    cleanup_temp_files();

    vASL large_asl_data;

    // Create a reasonably large number of stations (1000)
    const std::size_t num_stations = 1000;

    // Work directly with pointer format since copy constructor is private
    vASLPtr asl_ptr_data;
    for (std::size_t i = 0; i < num_stations; ++i) {
        ASLPtr asl_ptr;
        asl_ptr.reset(new CAStationList);
        asl_ptr->SetAssocMsrCount(static_cast<UINT32>(i % 10));
        asl_ptr->SetAMLStnIndex(static_cast<UINT32>(i));
        if (i % 2 == 0) {
            asl_ptr->SetValid();
        } else {
            asl_ptr->SetInvalid();
        }
        asl_ptr_data.push_back(asl_ptr);
    }

    // Resize large_asl_data for reading back
    large_asl_data.resize(num_stations);

    AslFile asl_writer(TEMP_ASL_FILE);
    asl_writer.Write(asl_ptr_data);
    REQUIRE(std::filesystem::exists(TEMP_ASL_FILE));

    // Read it back
    AslFile asl_reader(TEMP_ASL_FILE);
    auto loaded_result = asl_reader.Load();

    REQUIRE(loaded_result.count == num_stations);
    REQUIRE(loaded_result.stations.size() == num_stations);
    REQUIRE(loaded_result.free_stations.size() == num_stations);

    // Verify a few entries
    REQUIRE(loaded_result.stations[0].GetAssocMsrCount() == 0);
    REQUIRE(loaded_result.stations[0].GetAMLStnIndex() == 0);
    REQUIRE(loaded_result.stations[0].IsValid() == true);

    REQUIRE(loaded_result.stations[999].GetAssocMsrCount() == 9);
    REQUIRE(loaded_result.stations[999].GetAMLStnIndex() == 999);
    REQUIRE(loaded_result.stations[999].IsValid() == false);

    cleanup_temp_files();
}