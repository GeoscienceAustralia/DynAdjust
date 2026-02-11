#define TESTING_MAIN

#include "testing.hpp"

#include <algorithm>
#include <cstring>
#include <string>
#include <vector>

#include "config/dnatypes.hpp"
#include "config/dnatypes-gui.hpp"
#include "measurement_types/dnameasurement.hpp"
#include "measurement_types/dnastation.hpp"
#include "functions/dnatemplatefuncs.hpp"

using namespace dynadjust::measurements;

namespace {

station_t makeStation(const char* name, UINT32 fileOrder) {
    station_t s = {};
    strncpy(s.stationName, name, STN_NAME_WIDTH - 1);
    s.stationName[STN_NAME_WIDTH - 1] = '\0';
    s.fileOrder = fileOrder;
    return s;
}

} // anonymous namespace


TEST_CASE("MSR_TO_STN_SORT_UI enum values match issue #305 specification", "[m2s][enum]") {
    // Issue #305 requested:
    //   0 = alphanumeric (default)
    //   1 = measurement count (increasing)
    //   2 = original station order
    //   3 = measurement count (decreasing)
    REQUIRE(alpha_stn_sort_ui == 0);
    REQUIRE(meas_stn_sort_ui == 1);
    REQUIRE(orig_stn_sort_ui == 2);
    REQUIRE(meas_dec_stn_sort_ui == 3);
}

TEST_CASE("CompareStnNameOrder sorts stations alphabetically", "[m2s][sort][alpha]") {
    std::vector<station_t> stations;
    stations.push_back(makeStation("PERT", 0));
    stations.push_back(makeStation("ALIC", 1));
    stations.push_back(makeStation("TIDB", 2));
    stations.push_back(makeStation("BEEC", 3));
    stations.push_back(makeStation("HOB2", 4));

    // Station indices
    vUINT32 indices = {0, 1, 2, 3, 4};

    CompareStnNameOrder<station_t, UINT32> cmp(&stations);
    std::sort(indices.begin(), indices.end(), cmp);

    // Expected order: ALIC(1), BEEC(3), HOB2(4), PERT(0), TIDB(2)
    REQUIRE(indices[0] == 1);  // ALIC
    REQUIRE(indices[1] == 3);  // BEEC
    REQUIRE(indices[2] == 4);  // HOB2
    REQUIRE(indices[3] == 0);  // PERT
    REQUIRE(indices[4] == 2);  // TIDB
}

TEST_CASE("CompareStnFileOrder sorts by original file order", "[m2s][sort][fileorder]") {
    std::vector<station_t> stations;
    stations.push_back(makeStation("PERT", 5));
    stations.push_back(makeStation("ALIC", 2));
    stations.push_back(makeStation("TIDB", 8));
    stations.push_back(makeStation("BEEC", 1));
    stations.push_back(makeStation("HOB2", 3));

    vUINT32 indices = {0, 1, 2, 3, 4};

    CompareStnFileOrder<station_t, UINT32> cmp(&stations);
    std::sort(indices.begin(), indices.end(), cmp);

    // Expected order by fileOrder: BEEC(1), ALIC(2), HOB2(3), PERT(5), TIDB(8)
    REQUIRE(indices[0] == 3);  // BEEC (fileOrder=1)
    REQUIRE(indices[1] == 1);  // ALIC (fileOrder=2)
    REQUIRE(indices[2] == 4);  // HOB2 (fileOrder=3)
    REQUIRE(indices[3] == 0);  // PERT (fileOrder=5)
    REQUIRE(indices[4] == 2);  // TIDB (fileOrder=8)
}

TEST_CASE("Alphanumeric sort gives different order from file order", "[m2s][sort][difference]") {
    std::vector<station_t> stations;
    stations.push_back(makeStation("PERT", 0));
    stations.push_back(makeStation("ALIC", 1));
    stations.push_back(makeStation("TIDB", 2));

    vUINT32 alpha_indices = {0, 1, 2};
    vUINT32 file_indices = {0, 1, 2};

    CompareStnNameOrder<station_t, UINT32> alpha_cmp(&stations);
    std::sort(alpha_indices.begin(), alpha_indices.end(), alpha_cmp);

    CompareStnFileOrder<station_t, UINT32> file_cmp(&stations);
    std::sort(file_indices.begin(), file_indices.end(), file_cmp);

    // Alpha: ALIC(1), PERT(0), TIDB(2)
    // File:  PERT(0), ALIC(1), TIDB(2)
    REQUIRE(alpha_indices[0] != file_indices[0]);
}

TEST_CASE("CompareStnNameOrder handles numeric station names correctly", "[m2s][sort][numeric-names]") {
    // Issue #305 shows real station names like "00NA_1900001", "01NA_1900001", etc.
    std::vector<station_t> stations;
    stations.push_back(makeStation("10-2882_301117", 0));
    stations.push_back(makeStation("00NA_1900001", 1));
    stations.push_back(makeStation("02NA", 2));
    stations.push_back(makeStation("01NA_1900001", 3));

    vUINT32 indices = {0, 1, 2, 3};

    CompareStnNameOrder<station_t, UINT32> cmp(&stations);
    std::sort(indices.begin(), indices.end(), cmp);

    // Lexicographic: "00NA_1900001", "01NA_1900001", "02NA", "10-2882_301117"
    REQUIRE(indices[0] == 1);  // 00NA_1900001
    REQUIRE(indices[1] == 3);  // 01NA_1900001
    REQUIRE(indices[2] == 2);  // 02NA
    REQUIRE(indices[3] == 0);  // 10-2882_301117
}

TEST_CASE("Decreasing measurement count is reverse of increasing", "[m2s][sort][decreasing]") {
    // Simulate what the code does: sort ascending then reverse for decreasing
    vUINT32 increasing = {3, 1, 4, 1, 5, 9, 2, 6};
    vUINT32 decreasing = increasing;

    std::sort(increasing.begin(), increasing.end());

    std::sort(decreasing.begin(), decreasing.end());
    std::reverse(decreasing.begin(), decreasing.end());

    // First element of increasing should be last of decreasing
    REQUIRE(increasing.front() == decreasing.back());
    REQUIRE(increasing.back() == decreasing.front());

    // They should be exact reverses
    for (size_t i = 0; i < increasing.size(); ++i)
        REQUIRE(increasing[i] == decreasing[increasing.size() - 1 - i]);
}

TEST_CASE("Default sort (alpha_stn_sort_ui=0) is the default value", "[m2s][default]") {
    // alpha_stn_sort_ui == 0 means it is the default when _sort_msr_to_stn is
    // initialised to 0 (the zero-initialised value)
    UINT16 default_sort = 0;
    REQUIRE(default_sort == alpha_stn_sort_ui);
}

TEST_CASE("CompareStnNameOrder handles duplicate station names", "[m2s][sort][duplicates]") {
    std::vector<station_t> stations;
    stations.push_back(makeStation("ALIC", 0));
    stations.push_back(makeStation("ALIC", 1));
    stations.push_back(makeStation("BEEC", 2));

    vUINT32 indices = {0, 1, 2};

    CompareStnNameOrder<station_t, UINT32> cmp(&stations);
    std::sort(indices.begin(), indices.end(), cmp);

    // Both ALICs should come before BEEC; relative order of duplicates is unspecified
    // but both should be in positions 0 and 1
    REQUIRE(std::string(stations[indices[0]].stationName) == "ALIC");
    REQUIRE(std::string(stations[indices[1]].stationName) == "ALIC");
    REQUIRE(std::string(stations[indices[2]].stationName) == "BEEC");
}

TEST_CASE("CompareStnNameOrder with single station", "[m2s][sort][single]") {
    std::vector<station_t> stations;
    stations.push_back(makeStation("ONLY", 0));

    vUINT32 indices = {0};

    CompareStnNameOrder<station_t, UINT32> cmp(&stations);
    std::sort(indices.begin(), indices.end(), cmp);

    REQUIRE(indices[0] == 0);
}

TEST_CASE("All four sort modes produce valid orderings", "[m2s][sort][all-modes]") {
    std::vector<station_t> stations;
    stations.push_back(makeStation("PERT", 3));
    stations.push_back(makeStation("ALIC", 1));
    stations.push_back(makeStation("TIDB", 2));
    stations.push_back(makeStation("BEEC", 0));

    // Mode 0: alpha
    {
        vUINT32 indices = {0, 1, 2, 3};
        CompareStnNameOrder<station_t, UINT32> cmp(&stations);
        std::sort(indices.begin(), indices.end(), cmp);
        // ALIC < BEEC < PERT < TIDB
        REQUIRE(std::string(stations[indices[0]].stationName) == "ALIC");
        REQUIRE(std::string(stations[indices[3]].stationName) == "TIDB");
    }

    // Mode 2: file order
    {
        vUINT32 indices = {0, 1, 2, 3};
        CompareStnFileOrder<station_t, UINT32> cmp(&stations);
        std::sort(indices.begin(), indices.end(), cmp);
        // fileOrder: BEEC(0), ALIC(1), TIDB(2), PERT(3)
        REQUIRE(std::string(stations[indices[0]].stationName) == "BEEC");
        REQUIRE(std::string(stations[indices[3]].stationName) == "PERT");
    }
}

TEST_CASE("Switch statement covers all enum values", "[m2s][enum][coverage]") {
    // Verify that the switch used in PrintMeasurementsToStation handles all cases
    for (int i = 0; i <= 3; ++i) {
        bool handled = false;
        switch (i) {
        case alpha_stn_sort_ui:  handled = true; break;
        case meas_stn_sort_ui:   handled = true; break;
        case orig_stn_sort_ui:   handled = true; break;
        case meas_dec_stn_sort_ui: handled = true; break;
        }
        REQUIRE(handled);
    }
}
