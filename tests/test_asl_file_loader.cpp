#define TESTING_MAIN
#define __BINARY_NAME__ "test_AslFileLoader"
#define __BINARY_DESC__ "Unit tests for AslFileLoader class"

#include "testing.hpp"

#ifndef __BINARY_NAME__
#define __BINARY_NAME__ "test_AslFileLoader"
#endif
#ifndef __BINARY_DESC__
#define __BINARY_DESC__ "Unit tests for AslFileLoader class"
#endif

#include "../dynadjust/include/io/asl_file_loader.hpp"
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
TEST_CASE("AslFileLoader constructor", "[asl_file_loader][basic]") {
    AslFileLoader asl_loader;
    REQUIRE(true);
}

TEST_CASE("AslFileLoader copy constructor", "[asl_file_loader][basic]") {
    AslFileLoader asl_loader1;
    AslFileLoader asl_loader2(asl_loader1);
    REQUIRE(true);
}

TEST_CASE("AslFileLoader assignment operator", "[asl_file_loader][basic]") {
    AslFileLoader asl_loader1;
    AslFileLoader asl_loader2;
    asl_loader2 = asl_loader1;
    REQUIRE(true);
}

// Error handling tests
TEST_CASE("Error handling for non-existent files", "[asl_file_loader][error]") {
    AslFileLoader asl_loader;
    vASL asl_data;
    vUINT32 free_stations;

    auto result = asl_loader.LoadWithOptional(NONEXISTENT_FILE, &asl_data, &free_stations);
    REQUIRE(result == std::nullopt);
}

TEST_CASE("Error handling for write to invalid path", "[asl_file_loader][error]") {
    AslFileLoader asl_loader;
    vASLPtr asl_ptr_data;
    create_test_asl_ptr_data(&asl_ptr_data);

    bool threw_exception = false;
    try {
        asl_loader.WriteFile("/invalid/path/test.asl", &asl_ptr_data);
    } catch (const std::runtime_error&) { threw_exception = true; }
    REQUIRE(threw_exception);
}

// Load test file if available
TEST_CASE("Load ASL from test file if available", "[AslFileLoader][file]") {
    AslFileLoader asl_loader;
    vASL asl_data;
    vUINT32 free_stations;

    if (std::filesystem::exists(TEST_ASL_FILE)) {
        std::uint64_t count = asl_loader.LoadFile(TEST_ASL_FILE, &asl_data, &free_stations);
        REQUIRE(count > 0);
        REQUIRE(asl_data.size() == count);
        REQUIRE(free_stations.size() == count);

        // Verify free stations are initialized as incrementing vector
        for (std::size_t i = 0; i < free_stations.size(); ++i) {
            REQUIRE(free_stations[i] == i);
        }
    } else {
        REQUIRE(true); // Skip if file doesn't exist
    }
}

// Write and read back binary ASL data
TEST_CASE("Write and read back binary ASL data", "[AslFileLoader][roundtrip]") {
    AslFileLoader asl_loader;
    vASLPtr original_asl_data;

    cleanup_temp_files();

    create_test_asl_ptr_data(&original_asl_data);

    // Write the ASL file
    asl_loader.WriteFile(TEMP_ASL_FILE, &original_asl_data);
    REQUIRE(std::filesystem::exists(TEMP_ASL_FILE));

    // Read it back
    vASL loaded_asl_data;
    vUINT32 free_stations;
    std::uint64_t count = asl_loader.LoadFile(TEMP_ASL_FILE, &loaded_asl_data, &free_stations);

    // Verify data integrity
    REQUIRE(count == 3);
    REQUIRE(loaded_asl_data.size() == 3);
    REQUIRE(free_stations.size() == 3);

    // Verify station data
    REQUIRE(loaded_asl_data[0].GetAssocMsrCount() == 5);
    REQUIRE(loaded_asl_data[0].GetAMLStnIndex() == 0);
    REQUIRE(loaded_asl_data[0].IsValid() == true);

    REQUIRE(loaded_asl_data[1].GetAssocMsrCount() == 3);
    REQUIRE(loaded_asl_data[1].GetAMLStnIndex() == 1);
    REQUIRE(loaded_asl_data[1].IsValid() == true);

    REQUIRE(loaded_asl_data[2].GetAssocMsrCount() == 0);
    REQUIRE(loaded_asl_data[2].GetAMLStnIndex() == 2);
    REQUIRE(loaded_asl_data[2].IsValid() == false);

    // Verify free stations list is correctly initialized
    for (std::size_t i = 0; i < free_stations.size(); ++i) {
        REQUIRE(free_stations[i] == i);
    }

    cleanup_temp_files();
}

// Test empty ASL data
TEST_CASE("Write and read empty ASL data", "[AslFileLoader][empty]") {
    AslFileLoader asl_loader;
    vASLPtr empty_asl_data;

    cleanup_temp_files();

    // Write empty ASL file
    asl_loader.WriteFile(TEMP_ASL_FILE, &empty_asl_data);
    REQUIRE(std::filesystem::exists(TEMP_ASL_FILE));

    // Read it back
    vASL loaded_asl_data;
    vUINT32 free_stations;
    std::uint64_t count = asl_loader.LoadFile(TEMP_ASL_FILE, &loaded_asl_data, &free_stations);

    REQUIRE(count == 0);
    REQUIRE(loaded_asl_data.size() == 0);
    REQUIRE(free_stations.size() == 0);

    cleanup_temp_files();
}

// Test ASL text output
TEST_CASE("Write ASL text file", "[AslFileLoader][text]") {
    AslFileLoader asl_loader;
    vASLPtr asl_ptr_data;
    vdnaStnPtr stations;

    cleanup_temp_files();

    create_test_asl_ptr_data(&asl_ptr_data);
    create_test_station_data(stations);

    // Write text ASL file
    asl_loader.WriteTextFile(TEMP_ASL_TXT_FILE, &asl_ptr_data, &stations);
    REQUIRE(std::filesystem::exists(TEMP_ASL_TXT_FILE));

    // Verify file has content
    std::filesystem::file_size(TEMP_ASL_TXT_FILE);
    REQUIRE(std::filesystem::file_size(TEMP_ASL_TXT_FILE) > 0);

    // Read and verify basic content structure
    {
        std::ifstream txt_file(TEMP_ASL_TXT_FILE);
        std::string first_line;
        std::getline(txt_file, first_line);
        REQUIRE(first_line.find("3 stations") != std::string::npos);
    }

    cleanup_temp_files();
}

// Test with null pointers in ASL data
TEST_CASE("Handle null pointers in ASL data", "[AslFileLoader][nulls]") {
    AslFileLoader asl_loader;
    vASLPtr asl_ptr_data;

    cleanup_temp_files();

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
    asl_loader.WriteFile(TEMP_ASL_FILE, &asl_ptr_data);
    REQUIRE(std::filesystem::exists(TEMP_ASL_FILE));

    cleanup_temp_files();
}

// Test large station count (uint64_t usage)
TEST_CASE("Handle large station counts", "[AslFileLoader][large]") {
    AslFileLoader asl_loader;
    vASL large_asl_data;

    cleanup_temp_files();

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

    asl_loader.WriteFile(TEMP_ASL_FILE, &asl_ptr_data);
    REQUIRE(std::filesystem::exists(TEMP_ASL_FILE));

    // Read it back
    vASL loaded_asl_data;
    vUINT32 free_stations;
    std::uint64_t count = asl_loader.LoadFile(TEMP_ASL_FILE, &loaded_asl_data, &free_stations);

    REQUIRE(count == num_stations);
    REQUIRE(loaded_asl_data.size() == num_stations);
    REQUIRE(free_stations.size() == num_stations);

    // Verify a few entries
    REQUIRE(loaded_asl_data[0].GetAssocMsrCount() == 0);
    REQUIRE(loaded_asl_data[0].GetAMLStnIndex() == 0);
    REQUIRE(loaded_asl_data[0].IsValid() == true);

    REQUIRE(loaded_asl_data[999].GetAssocMsrCount() == 9);
    REQUIRE(loaded_asl_data[999].GetAMLStnIndex() == 999);
    REQUIRE(loaded_asl_data[999].IsValid() == false);

    cleanup_temp_files();
}
