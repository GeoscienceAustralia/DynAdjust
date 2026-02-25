#define TESTING_MAIN
#define __BINARY_NAME__ "test_snx_file_writer"
#define __BINARY_DESC__ "Unit tests for SNX file writer and SINEX date formatting"

#include "testing.hpp"

#ifndef __BINARY_NAME__
#define __BINARY_NAME__ "test_snx_file_writer"
#endif
#ifndef __BINARY_DESC__
#define __BINARY_DESC__ "Unit tests for SNX file writer and SINEX date formatting"
#endif

#include "../dynadjust/include/io/dnaiosnx.hpp"
#include "../dynadjust/include/functions/dnatemplatedatetimefuncs.hpp"
#include "../dynadjust/include/functions/dnachronutils.hpp"
#include "../dynadjust/include/math/dnamatrix_contiguous.hpp"
#include "../dynadjust/include/parameters/dnadatum.hpp"
#include "../dynadjust/include/config/dnaoptions.hpp"
#include "../dynadjust/include/functions/dnastrmanipfuncs.hpp"
#include "../dynadjust/include/functions/dnaiostreamfuncs.hpp"
#include <filesystem>
#include <fstream>
#include <sstream>
#include <cstring>

using namespace dynadjust::iostreams;
using namespace dynadjust::math;
using namespace dynadjust::datum_parameters;

namespace {
const std::string TEMP_SNX_FILE = "temp_test.snx";
const std::string TEMP_WARNING_FILE = "temp_test_warnings.txt";

// Helper to create test stations
void create_test_stations(vstn_t& stations) {
    station_t stn1;
    strcpy(stn1.stationName, "STA1");
    strcpy(stn1.description, "Test Station 1");
    stn1.currentLatitude = -0.6108652381980153;  // ~35 degrees S
    stn1.currentLongitude = 2.5307274153917776;  // ~145 degrees E
    stn1.currentHeight = 100.5;
    
    station_t stn2;
    strcpy(stn2.stationName, "LONGNAME");  // Deliberately long name
    strcpy(stn2.description, "Test Station with Long Name");
    stn2.currentLatitude = -0.6283185307179586;  // ~36 degrees S
    stn2.currentLongitude = 2.5481807079117210;  // ~146 degrees E
    stn2.currentHeight = 200.7;
    
    station_t stn3;
    strcpy(stn3.stationName, "STA3");
    strcpy(stn3.description, "Test Station 3");
    stn3.currentLatitude = -0.5934119456780721;  // ~34 degrees S
    stn3.currentLongitude = 2.5132741228718345;  // ~144 degrees E
    stn3.currentHeight = 50.2;
    
    stations.push_back(stn1);
    stations.push_back(stn2);
    stations.push_back(stn3);
}

// Helper to create test estimates and variances
void create_test_matrices(matrix_2d& estimates, matrix_2d& variances, size_t num_stations) {
    size_t dim = num_stations * 3;
    
    // Initialize estimates (X, Y, Z for each station)
    estimates.redim(dim, 1);
    for (size_t i = 0; i < num_stations; ++i) {
        estimates.put(i * 3, 0, -4000000.0 + i * 1000);      // X
        estimates.put(i * 3 + 1, 0, 2500000.0 + i * 1000);   // Y
        estimates.put(i * 3 + 2, 0, -3700000.0 + i * 1000);  // Z
    }
    
    // Initialize variances (diagonal matrix with some covariances)
    variances.redim(dim, dim);
    variances.zero();
    for (size_t i = 0; i < dim; ++i) {
        variances.put(i, i, 0.0001 + i * 0.00001);  // Diagonal elements
        // Add some off-diagonal elements
        if (i > 0) {
            variances.put(i, i - 1, 0.00001);
            variances.put(i - 1, i, 0.00001);
        }
    }
}

// Helper to create test metadata
void create_test_metadata(binary_file_meta_t& bst_meta, binary_file_meta_t& bms_meta) {
    // BST metadata
    strcpy(bst_meta.modifiedBy, "test_program");
    strcpy(bst_meta.epsgCode, "7844");
    strcpy(bst_meta.epoch, "01.01.2020");
    bst_meta.inputFileCount = 1;
    bst_meta.inputFileMeta = new input_file_meta_t[1];
    strcpy(bst_meta.inputFileMeta[0].filename, "test_stations.xml");
    
    // BMS metadata
    strcpy(bms_meta.modifiedBy, "test_program");
    strcpy(bms_meta.epsgCode, "7844");
    strcpy(bms_meta.epoch, "01.01.2020");
    bms_meta.inputFileCount = 1;
    bms_meta.inputFileMeta = new input_file_meta_t[1];
    strcpy(bms_meta.inputFileMeta[0].filename, "test_measurements.xml");
}

// Helper to cleanup metadata - prevents double free
void cleanup_metadata(binary_file_meta_t& bst_meta, binary_file_meta_t& bms_meta) {
    if (bst_meta.inputFileMeta != nullptr) {
        delete[] bst_meta.inputFileMeta;
        bst_meta.inputFileMeta = nullptr;
        bst_meta.inputFileCount = 0;
    }
    if (bms_meta.inputFileMeta != nullptr) {
        delete[] bms_meta.inputFileMeta;
        bms_meta.inputFileMeta = nullptr;
        bms_meta.inputFileCount = 0;
    }
}

// Helper to create block station mappings
void create_block_mappings(uint32_uint32_map& blockStationsMap, vUINT32& blockStations, size_t num_stations) {
    for (size_t i = 0; i < num_stations; ++i) {
        blockStationsMap[i] = i;
        blockStations.push_back(i);
    }
}

// Cleanup function
void cleanup_temp_files() {
    if (std::filesystem::exists(TEMP_SNX_FILE)) {
        std::filesystem::remove(TEMP_SNX_FILE);
    }
    if (std::filesystem::exists(TEMP_WARNING_FILE)) {
        std::filesystem::remove(TEMP_WARNING_FILE);
    }
}

// Helper to read file into string
std::string read_file_contents(const std::string& filename) {
    std::ifstream file(filename);
    std::stringstream buffer;
    buffer << file.rdbuf();
    return buffer.str();
}

// Helper to verify SINEX format line
bool verify_sinex_line_format(const std::string& line, const std::string& expected_start) {
    return line.substr(0, expected_start.length()) == expected_start;
}

} // namespace

// Test dateSINEXFormat function with various dates
TEST_CASE("dateSINEXFormat standard dates", "[sinex][date]") {
    std::stringstream ss;
    using namespace dynadjust;
    
    // Test 1: April 30, 1995 should produce "95:120:00000"
    {
        ss.str("");
        auto d1 = MakeDate(1995, 4, 30);
        DateSINEXFormat(&ss, d1);
        REQUIRE(ss.str() == "95:120:00000");
    }
    
    // Test 2: January 1, 2000 should produce "00:001:00000"
    {
        ss.str("");
        auto d2 = MakeDate(2000, 1, 1);
        DateSINEXFormat(&ss, d2);
        REQUIRE(ss.str() == "00:001:00000");
    }
    
    // Test 3: December 31, 1999 should produce "99:365:00000"
    {
        ss.str("");
        auto d3 = MakeDate(1999, 12, 31);
        DateSINEXFormat(&ss, d3);
        REQUIRE(ss.str() == "99:365:00000");
    }
}

TEST_CASE("dateSINEXFormat leap year handling", "[sinex][date]") {
    std::stringstream ss;
    using namespace dynadjust;
    
    // Test 1: Feb 29, 2020 (leap year) should produce "20:060:00000"
    {
        ss.str("");
        auto d1 = MakeDate(2020, 2, 29);
        DateSINEXFormat(&ss, d1);
        REQUIRE(ss.str() == "20:060:00000");
    }
    
    // Test 2: Dec 31, 2020 (leap year) should produce "20:366:00000"
    {
        ss.str("");
        auto d2 = MakeDate(2020, 12, 31);
        DateSINEXFormat(&ss, d2);
        REQUIRE(ss.str() == "20:366:00000");
    }
    
    // Test 3: Feb 28, 2021 (non-leap year) should produce "21:059:00000"
    {
        ss.str("");
        auto d3 = MakeDate(2021, 2, 28);
        DateSINEXFormat(&ss, d3);
        REQUIRE(ss.str() == "21:059:00000");
    }
    
    // Test 4: Dec 31, 2021 (non-leap year) should produce "21:365:00000"
    {
        ss.str("");
        auto d4 = MakeDate(2021, 12, 31);
        DateSINEXFormat(&ss, d4);
        REQUIRE(ss.str() == "21:365:00000");
    }
}

TEST_CASE("dateSINEXFormat year transitions", "[sinex][date]") {
    std::stringstream ss;
    using namespace dynadjust;
    
    // Test year 2-digit formatting across century
    {
        ss.str("");
        auto d1 = MakeDate(1998, 1, 1);
        DateSINEXFormat(&ss, d1);
        REQUIRE(ss.str() == "98:001:00000");
    }
    
    {
        ss.str("");
        auto d2 = MakeDate(2001, 1, 1);
        DateSINEXFormat(&ss, d2);
        REQUIRE(ss.str() == "01:001:00000");
    }
    
    {
        ss.str("");
        auto d3 = MakeDate(2099, 12, 31);
        DateSINEXFormat(&ss, d3);
        REQUIRE(ss.str() == "99:365:00000");
    }
}

TEST_CASE("dateSINEXFormat with seconds calculation", "[sinex][date]") {
    std::stringstream ss;
    using namespace dynadjust;
    auto d = MakeDate(2020, 6, 15);
    
    // Test with seconds calculation enabled
    // Note: This will calculate seconds from the date to current time
    // We can't test exact value, but we can verify format
    DateSINEXFormat(&ss, d, true);
    std::string result = ss.str();
    
    // Should be in format YY:DOY:SSSSS where SSSSS is at least 5 characters
    REQUIRE(result.substr(0, 3) == "20:");
    REQUIRE(result.substr(3, 4) == "167:");  // June 15 is day 167 in leap year
    REQUIRE(result.length() >= 12);  // YY:DOY:SSSSS = 2+1+3+1+5 = 12
    
    // Verify the seconds field contains only digits or spaces
    std::string seconds_field = result.substr(7);
    for (char c : seconds_field) {
        REQUIRE((std::isdigit(c) || c == ' ' || c == '-'));
    }
}

// Test date helper functions
TEST_CASE("Date helper functions", "[sinex][date]") {
    using namespace dynadjust;
    auto d = MakeDate(2020, 6, 15);
    
    REQUIRE(DateYear<int>(d) == 2020);
    REQUIRE(DateMonth<int>(d) == 6);
    REQUIRE(DateDay<int>(d) == 15);
    REQUIRE(DateDOY<int>(d) == 167);  // June 15 is day 167 in leap year 2020
}

TEST_CASE("yearFraction calculations", "[sinex][date]") {
    using namespace dynadjust;
    // Test 1: January 1 should give near 0 fraction
    {
        auto d1 = MakeDate(2020, 1, 1);
        double fraction = YearFraction<double>(d1);
        // (1 - 0.5) / 366 = 0.5 / 366 ≈ 0.00137
        REQUIRE(fraction > 0.001);
        REQUIRE(fraction < 0.002);
    }
    
    // Test 2: Mid-year (approximately)
    {
        auto d2 = MakeDate(2020, 7, 1);  // Day 183 in leap year
        double fraction = YearFraction<double>(d2);
        // (183 - 0.5) / 366 ≈ 0.499
        REQUIRE(fraction > 0.498);
        REQUIRE(fraction < 0.501);
    }
    
    // Test 3: Last day of year
    {
        auto d3 = MakeDate(2020, 12, 31);  // Day 366
        double fraction = YearFraction<double>(d3);
        // (366 - 0.5) / 366 ≈ 0.9986
        REQUIRE(fraction > 0.998);
        REQUIRE(fraction < 1.0);
    }
}

TEST_CASE("referenceEpoch calculations", "[sinex][date]") {
    using namespace dynadjust;
    // Test reference epoch calculation (year + fraction)
    {
        auto d1 = MakeDate(2020, 1, 1);
        double epoch = ReferenceEpoch<double>(d1);
        REQUIRE(epoch > 2020.0);
        REQUIRE(epoch < 2020.01);
    }
    
    {
        auto d2 = MakeDate(2020, 12, 31);
        double epoch = ReferenceEpoch<double>(d2);
        REQUIRE(epoch > 2020.99);
        REQUIRE(epoch < 2021.0);
    }
}

TEST_CASE("elapsedTime calculations", "[sinex][date]") {
    using namespace dynadjust;
    auto d1 = MakeDate(2020, 1, 1);
    auto d2 = MakeDate(2020, 12, 31);
    
    double elapsed = ElapsedTime<double>(d2, d1);
    // Should be close to 1 year (365/366 ≈ 0.9973)
    REQUIRE(elapsed > 0.99);
    REQUIRE(elapsed < 1.0);
}

// Test DnaIoSnx class functionality
TEST_CASE("DnaIoSnx constructor", "[sinex][basic]") {
    DnaIoSnx snx_writer;
    REQUIRE(true);
}

TEST_CASE("DnaIoSnx warning system", "[sinex][warnings]") {
    cleanup_temp_files();
    
    DnaIoSnx snx_writer;
    vstn_t stations;
    create_test_stations(stations);
    
    // Create necessary data structures
    project_settings p;
    p.g.network_name = "TEST_NET";
    
    binary_file_meta_t bst_meta, bms_meta;
    create_test_metadata(bst_meta, bms_meta);
    
    matrix_2d estimates, variances;
    create_test_matrices(estimates, variances, stations.size());
    
    uint32_uint32_map blockStationsMap;
    vUINT32 blockStations;
    create_block_mappings(blockStationsMap, blockStations, stations.size());
    
    CDnaDatum datum;
    datum.SetEpoch(boost::gregorian::date(2020, 1, 1));
    
    UINT32 measurementParams = 30;
    UINT32 unknownParams = 9;
    double sigmaZero = 1.234567;
    
    // Write SINEX file - this should generate a warning for "LONGNAME"
    std::ofstream snx_file(TEMP_SNX_FILE);
    snx_writer.SerialiseSinex(&snx_file, &stations, bst_meta, bms_meta,
                               &estimates, &variances, p,
                               measurementParams, unknownParams, sigmaZero,
                               &blockStationsMap, &blockStations,
                               1, 0, &datum);
    snx_file.close();
    
    // Check that warnings were generated
    REQUIRE(snx_writer.WarningsExist() == true);
    
    // Write warnings to file
    std::ofstream warning_file(TEMP_WARNING_FILE);
    snx_writer.PrintWarnings(&warning_file, TEMP_SNX_FILE);
    warning_file.close();
    
    // Verify warning file contains expected warning
    std::string warnings = read_file_contents(TEMP_WARNING_FILE);
    REQUIRE(warnings.find("LONGNAME") != std::string::npos);
    REQUIRE(warnings.find("exceeds four characters") != std::string::npos);
    
    // Cleanup
    cleanup_metadata(bst_meta, bms_meta);
    cleanup_temp_files();
}

TEST_CASE("Write complete SINEX file", "[sinex][integration]") {
    cleanup_temp_files();
    
    DnaIoSnx snx_writer;
    vstn_t stations;
    create_test_stations(stations);
    
    project_settings p;
    p.g.network_name = "TEST_NET";
    
    binary_file_meta_t bst_meta, bms_meta;
    create_test_metadata(bst_meta, bms_meta);
    
    matrix_2d estimates, variances;
    create_test_matrices(estimates, variances, stations.size());
    
    uint32_uint32_map blockStationsMap;
    vUINT32 blockStations;
    create_block_mappings(blockStationsMap, blockStations, stations.size());
    
    CDnaDatum datum;
    datum.SetEpoch(boost::gregorian::date(2020, 1, 1));
    
    UINT32 measurementParams = 30;
    UINT32 unknownParams = 9;
    double sigmaZero = 1.234567;
    
    // Write SINEX file
    std::ofstream snx_file(TEMP_SNX_FILE);
    snx_writer.SerialiseSinex(&snx_file, &stations, bst_meta, bms_meta,
                               &estimates, &variances, p,
                               measurementParams, unknownParams, sigmaZero,
                               &blockStationsMap, &blockStations,
                               1, 0, &datum);
    snx_file.close();
    
    // Read and verify file structure
    std::string contents = read_file_contents(TEMP_SNX_FILE);
    
    // Check header
    REQUIRE(contents.find("%=SNX 2.00 DNA") != std::string::npos);
    
    // Check blocks exist
    REQUIRE(contents.find("+FILE/REFERENCE") != std::string::npos);
    REQUIRE(contents.find("-FILE/REFERENCE") != std::string::npos);
    REQUIRE(contents.find("+FILE/COMMENT") != std::string::npos);
    REQUIRE(contents.find("-FILE/COMMENT") != std::string::npos);
    REQUIRE(contents.find("+SITE/ID") != std::string::npos);
    REQUIRE(contents.find("-SITE/ID") != std::string::npos);
    REQUIRE(contents.find("+SOLUTION/STATISTICS") != std::string::npos);
    REQUIRE(contents.find("-SOLUTION/STATISTICS") != std::string::npos);
    REQUIRE(contents.find("+SOLUTION/ESTIMATE") != std::string::npos);
    REQUIRE(contents.find("-SOLUTION/ESTIMATE") != std::string::npos);
    REQUIRE(contents.find("+SOLUTION/MATRIX_ESTIMATE L COVA") != std::string::npos);
    REQUIRE(contents.find("-SOLUTION/MATRIX_ESTIMATE L COVA") != std::string::npos);
    REQUIRE(contents.find("%ENDSNX") != std::string::npos);
    
    // Check epoch formatting
    REQUIRE(contents.find("20:001:00000") != std::string::npos);  // 2020-01-01
    
    // Check statistics
    REQUIRE(contents.find("NUMBER OF OBSERVATIONS") != std::string::npos);
    REQUIRE(contents.find("30") != std::string::npos);  // measurementParams
    REQUIRE(contents.find("NUMBER OF UNKNOWNS") != std::string::npos);
    REQUIRE(contents.find("9") != std::string::npos);   // unknownParams
    REQUIRE(contents.find("VARIANCE FACTOR") != std::string::npos);
    REQUIRE(contents.find("1.234567") != std::string::npos);  // sigmaZero
    
    // Cleanup
    cleanup_metadata(bst_meta, bms_meta);
    cleanup_temp_files();
}

TEST_CASE("SINEX site ID formatting", "[sinex][formatting]") {
    cleanup_temp_files();
    
    DnaIoSnx snx_writer;
    vstn_t stations;
    
    // Create station with specific coordinates for testing
    station_t stn;
    strcpy(stn.stationName, "TEST");
    strcpy(stn.description, "Test Station Description");
    stn.currentLatitude = -0.6108652381980153;   // -35 degrees
    stn.currentLongitude = 2.5307274153917776;   // 145 degrees
    stn.currentHeight = 123.4;
    stations.push_back(stn);
    
    project_settings p;
    p.g.network_name = "TEST_NET";
    
    binary_file_meta_t bst_meta, bms_meta;
    create_test_metadata(bst_meta, bms_meta);
    
    matrix_2d estimates, variances;
    create_test_matrices(estimates, variances, 1);
    
    uint32_uint32_map blockStationsMap;
    vUINT32 blockStations;
    blockStationsMap[0] = 0;
    blockStations.push_back(0);
    
    CDnaDatum datum;
    datum.SetEpoch(boost::gregorian::date(2020, 1, 1));
    
    UINT32 measurementParams = 10;
    UINT32 unknownParams = 3;
    double sigmaZero = 1.0;
    
    std::ofstream snx_file(TEMP_SNX_FILE);
    snx_writer.SerialiseSinex(&snx_file, &stations, bst_meta, bms_meta,
                               &estimates, &variances, p,
                               measurementParams, unknownParams, sigmaZero,
                               &blockStationsMap, &blockStations,
                               1, 0, &datum);
    snx_file.close();
    
    std::string contents = read_file_contents(TEMP_SNX_FILE);
    
    // Find SITE/ID section
    size_t site_id_pos = contents.find("+SITE/ID");
    REQUIRE(site_id_pos != std::string::npos);
    
    // Check for proper formatting of station entry
    // The actual format has extra space after TEST
    REQUIRE(contents.find(" TEST  A TEST      P Test Station Descripti") != std::string::npos);
    
    // Cleanup
    cleanup_metadata(bst_meta, bms_meta);
    cleanup_temp_files();
}

TEST_CASE("SINEX solution estimates formatting", "[sinex][formatting]") {
    cleanup_temp_files();
    
    DnaIoSnx snx_writer;
    vstn_t stations;
    
    station_t stn;
    strcpy(stn.stationName, "STA1");
    strcpy(stn.description, "Station 1");
    stations.push_back(stn);
    
    project_settings p;
    p.g.network_name = "TEST_NET";
    
    binary_file_meta_t bst_meta, bms_meta;
    create_test_metadata(bst_meta, bms_meta);
    
    // Create specific estimates for testing
    matrix_2d estimates(3, 1);
    estimates.put(0, 0, -3999999.12345678901234);  // X
    estimates.put(1, 0, 2500000.98765432109876);   // Y
    estimates.put(2, 0, -3700000.11111111111111);  // Z
    
    matrix_2d variances(3, 3);
    variances.zero();
    variances.put(0, 0, 0.0001);  // X variance
    variances.put(1, 1, 0.0004);  // Y variance
    variances.put(2, 2, 0.0009);  // Z variance
    
    uint32_uint32_map blockStationsMap;
    vUINT32 blockStations;
    blockStationsMap[0] = 0;
    blockStations.push_back(0);
    
    CDnaDatum datum;
    datum.SetEpoch(boost::gregorian::date(2020, 1, 1));
    
    UINT32 measurementParams = 10;
    UINT32 unknownParams = 3;
    double sigmaZero = 1.0;
    
    std::ofstream snx_file(TEMP_SNX_FILE);
    snx_writer.SerialiseSinex(&snx_file, &stations, bst_meta, bms_meta,
                               &estimates, &variances, p,
                               measurementParams, unknownParams, sigmaZero,
                               &blockStationsMap, &blockStations,
                               1, 0, &datum);
    snx_file.close();
    
    std::string contents = read_file_contents(TEMP_SNX_FILE);
    
    // Check scientific notation formatting
    REQUIRE(contents.find("E+") != std::string::npos);  // Scientific notation
    
    // Check standard deviations (sqrt of variances)
    REQUIRE(contents.find("1.00000E-02") != std::string::npos);  // sqrt(0.0001) = 0.01
    REQUIRE(contents.find("2.00000E-02") != std::string::npos);  // sqrt(0.0004) = 0.02
    REQUIRE(contents.find("3.00000E-02") != std::string::npos);  // sqrt(0.0009) = 0.03
    
    // Cleanup
    cleanup_metadata(bst_meta, bms_meta);
    cleanup_temp_files();
}

TEST_CASE("SINEX variance matrix formatting", "[sinex][formatting]") {
    cleanup_temp_files();
    
    DnaIoSnx snx_writer;
    vstn_t stations;
    
    // Create 2 stations for 6x6 variance matrix
    station_t stn1, stn2;
    strcpy(stn1.stationName, "STA1");
    strcpy(stn2.stationName, "STA2");
    stations.push_back(stn1);
    stations.push_back(stn2);
    
    project_settings p;
    p.g.network_name = "TEST_NET";
    
    binary_file_meta_t bst_meta, bms_meta;
    create_test_metadata(bst_meta, bms_meta);
    
    matrix_2d estimates(6, 1);
    estimates.zero();
    
    // Create a 6x6 variance matrix with specific pattern
    matrix_2d variances(6, 6);
    variances.zero();
    
    // Fill lower triangular with test values
    for (size_t i = 0; i < 6; ++i) {
        for (size_t j = 0; j <= i; ++j) {
            variances.put(i, j, 0.0001 * (i + 1) * (j + 1));
        }
    }
    
    uint32_uint32_map blockStationsMap;
    vUINT32 blockStations;
    blockStationsMap[0] = 0;
    blockStationsMap[1] = 1;
    blockStations.push_back(0);
    blockStations.push_back(1);
    
    CDnaDatum datum;
    UINT32 measurementParams = 20;
    UINT32 unknownParams = 6;
    double sigmaZero = 1.0;
    
    std::ofstream snx_file(TEMP_SNX_FILE);
    snx_writer.SerialiseSinex(&snx_file, &stations, bst_meta, bms_meta,
                               &estimates, &variances, p,
                               measurementParams, unknownParams, sigmaZero,
                               &blockStationsMap, &blockStations,
                               1, 0, &datum);
    snx_file.close();
    
    std::string contents = read_file_contents(TEMP_SNX_FILE);
    
    // Check for matrix block
    REQUIRE(contents.find("+SOLUTION/MATRIX_ESTIMATE L COVA") != std::string::npos);
    
    // Check matrix indices format (1-based)
    REQUIRE(contents.find("     1     1 ") != std::string::npos);  // First element
    REQUIRE(contents.find("     6     1 ") != std::string::npos);  // Row 6, column 1
    
    // Cleanup
    cleanup_metadata(bst_meta, bms_meta);
    cleanup_temp_files();
}

TEST_CASE("Handle empty data", "[sinex][edge]") {
    cleanup_temp_files();
    
    DnaIoSnx snx_writer;
    vstn_t stations;  // Empty stations
    
    project_settings p;
    p.g.network_name = "EMPTY_NET";
    
    binary_file_meta_t bst_meta, bms_meta;
    create_test_metadata(bst_meta, bms_meta);
    
    matrix_2d estimates(0, 0);
    matrix_2d variances(0, 0);
    
    uint32_uint32_map blockStationsMap;
    vUINT32 blockStations;
    
    CDnaDatum datum;
    UINT32 measurementParams = 0;
    UINT32 unknownParams = 0;
    double sigmaZero = 0.0;
    
    // Should handle empty data gracefully
    std::ofstream snx_file(TEMP_SNX_FILE);
    snx_writer.SerialiseSinex(&snx_file, &stations, bst_meta, bms_meta,
                               &estimates, &variances, p,
                               measurementParams, unknownParams, sigmaZero,
                               &blockStationsMap, &blockStations,
                               1, 0, &datum);
    snx_file.close();
    
    // Verify file was created and has basic structure
    std::string contents = read_file_contents(TEMP_SNX_FILE);
    REQUIRE(contents.find("%=SNX 2.00") != std::string::npos);
    REQUIRE(contents.find("%ENDSNX") != std::string::npos);
    
    // Cleanup
    cleanup_metadata(bst_meta, bms_meta);
    cleanup_temp_files();
}

TEST_CASE("Multi-block SINEX output", "[sinex][multiblock]") {
    cleanup_temp_files();
    
    DnaIoSnx snx_writer;
    vstn_t stations;
    create_test_stations(stations);
    
    project_settings p;
    p.g.network_name = "MULTIBLOCK_NET";
    
    binary_file_meta_t bst_meta, bms_meta;
    create_test_metadata(bst_meta, bms_meta);
    
    matrix_2d estimates, variances;
    create_test_matrices(estimates, variances, stations.size());
    
    uint32_uint32_map blockStationsMap;
    vUINT32 blockStations;
    create_block_mappings(blockStationsMap, blockStations, stations.size());
    
    CDnaDatum datum;
    UINT32 measurementParams = 30;
    UINT32 unknownParams = 9;
    double sigmaZero = 1.0;
    
    // Test with block 2 of 5
    std::ofstream snx_file(TEMP_SNX_FILE);
    snx_writer.SerialiseSinex(&snx_file, &stations, bst_meta, bms_meta,
                               &estimates, &variances, p,
                               measurementParams, unknownParams, sigmaZero,
                               &blockStationsMap, &blockStations,
                               5, 1, &datum);  // 5 blocks total, this is block 2 (0-indexed)
    snx_file.close();
    
    std::string contents = read_file_contents(TEMP_SNX_FILE);
    
    // Check for multi-block comments
    REQUIRE(contents.find("block 2 of") != std::string::npos);
    REQUIRE(contents.find("segmented") != std::string::npos);
    REQUIRE(contents.find("5 blocks") != std::string::npos);
    
    // Cleanup
    cleanup_metadata(bst_meta, bms_meta);
    cleanup_temp_files();
}

// Test edge cases in date formatting
TEST_CASE("Date edge cases", "[sinex][date][edge]") {
    std::stringstream ss;
    using namespace dynadjust;
    
    // Test earliest supported date
    {
        ss.str("");
        auto d1 = MakeDate(1900, 1, 1);
        DateSINEXFormat(&ss, d1);
        REQUIRE(ss.str() == "00:001:00000");  // Year 1900 shows as 00
    }
    
    // Test century year (divisible by 100 but not 400 - not a leap year)
    {
        ss.str("");
        auto d2 = MakeDate(1900, 12, 31);
        DateSINEXFormat(&ss, d2);
        REQUIRE(ss.str() == "00:365:00000");  // Not 366
    }
    
    // Test century leap year (divisible by 400)
    {
        ss.str("");
        auto d3 = MakeDate(2000, 2, 29);
        DateSINEXFormat(&ss, d3);
        REQUIRE(ss.str() == "00:060:00000");  // Feb 29 exists
    }
    
    // Test day padding
    {
        ss.str("");
        auto d4 = MakeDate(2020, 1, 9);  // Day 9 - should be zero-padded
        DateSINEXFormat(&ss, d4);
        REQUIRE(ss.str() == "20:009:00000");
    }
}

// Test date parsing functions
TEST_CASE("Date parsing from string", "[sinex][date][parsing]") {
    using namespace dynadjust;
    // Test dateFromString_safe with dd.mm.yyyy format
    {
        std::string dateStr = "15.06.2020";
        auto d = DateFromStringSafe(dateStr);
        REQUIRE(DateYear<int>(d) == 2020);
        REQUIRE(DateMonth<int>(d) == 6);
        REQUIRE(DateDay<int>(d) == 15);
    }
    
    // Test with single digits
    {
        std::string dateStr = "1.1.2020";
        auto d = DateFromStringSafe(dateStr);
        REQUIRE(DateYear<int>(d) == 2020);
        REQUIRE(DateMonth<int>(d) == 1);
        REQUIRE(DateDay<int>(d) == 1);
    }
    
    // Test invalid date handling
    {
        bool threw_exception = false;
        try {
            std::string dateStr = "invalid.date.string";
            DateFromStringSafe(dateStr);
        } catch (const std::runtime_error&) {
            threw_exception = true;
        }
        REQUIRE(threw_exception);
    }
}

// Test special date calculations
TEST_CASE("Special date calculations", "[sinex][date][calculations]") {
    using namespace dynadjust;
    // Test year_doy_Average function
    {
        UINT32 start_year = 2020, end_year = 2020;
        UINT32 start_doy = 100, end_doy = 200;
        UINT32 avg_year, avg_doy;
        
        YearDoyAverage(start_year, end_year, start_doy, end_doy, avg_year, avg_doy);
        
        REQUIRE(avg_year == 2020);
        REQUIRE(avg_doy == 150);  // Average of 100 and 200
    }
    
    // Test year crossover
    {
        UINT32 start_year = 2020, end_year = 2021;
        UINT32 start_doy = 350, end_doy = 20;
        UINT32 avg_year, avg_doy;
        
        YearDoyAverage(start_year, end_year, start_doy, end_doy, avg_year, avg_doy);
        
        // average(350, 350 + 20) = average(350, 370) = 360
        // Since 360 < 365, year stays as start_year (2020)
        REQUIRE(avg_year == 2020);
        REQUIRE(avg_doy == 360);
    }
    
    // Test year crossover with result in next year
    {
        UINT32 start_year = 2020, end_year = 2021;
        UINT32 start_doy = 360, end_doy = 20;
        UINT32 avg_year, avg_doy;
        
        YearDoyAverage(start_year, end_year, start_doy, end_doy, avg_year, avg_doy);
        
        // average(360, 360 + 20) = average(360, 380) = 370
        // Since 370 > 365, year = 2021, doy = 370 - 365 = 5
        REQUIRE(avg_year == 2021);
        REQUIRE(avg_doy == 5);
    }
}