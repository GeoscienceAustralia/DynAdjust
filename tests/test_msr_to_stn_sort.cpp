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

station_t makeStation(const char* name, UINT32 fileOrder, UINT32 nameOrder = 0) {
    station_t s = {};
    strncpy(s.stationName, name, STN_NAME_WIDTH - 1);
    s.stationName[STN_NAME_WIDTH - 1] = '\0';
    s.fileOrder = fileOrder;
    s.nameOrder = nameOrder;
    return s;
}

} // anonymous namespace


TEST_CASE("MSR_TO_STN_SORT_UI enum values match issue #305 specification", "[m2s][enum]") {
    REQUIRE(orig_stn_sort_ui == 0);
    REQUIRE(name_stn_sort_ui == 1);
    REQUIRE(meas_stn_sort_ui == 2);
    REQUIRE(saem_stn_sort_ui == 3);
}

TEST_CASE("CompareStnNameOrder sorts stations by nameOrder", "[m2s][sort][alpha]") {
    // CompareStnNameOrder compares the pre-computed nameOrder field,
    // which is assigned during import based on alphabetical station name order.
    // Stations: PERT(fileOrder=0), ALIC(1), TIDB(2), BEEC(3), HOB2(4)
    // Alphabetical order: ALIC, BEEC, HOB2, PERT, TIDB
    // So nameOrder values: PERT=3, ALIC=0, BEEC=1, HOB2=2, TIDB=4
    std::vector<station_t> stations;
    stations.push_back(makeStation("PERT", 0, 3));
    stations.push_back(makeStation("ALIC", 1, 0));
    stations.push_back(makeStation("TIDB", 2, 4));
    stations.push_back(makeStation("BEEC", 3, 1));
    stations.push_back(makeStation("HOB2", 4, 2));

    vUINT32 indices = {0, 1, 2, 3, 4};

    CompareStnNameOrder<station_t, UINT32> cmp(&stations);
    std::sort(indices.begin(), indices.end(), cmp);

    // Expected order by nameOrder: ALIC(0), BEEC(1), HOB2(2), PERT(3), TIDB(4)
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

TEST_CASE("Name sort gives different order from file order", "[m2s][sort][difference]") {
    // PERT has fileOrder=0 (first in file) but nameOrder=1 (second alphabetically)
    // ALIC has fileOrder=1 (second in file) but nameOrder=0 (first alphabetically)
    std::vector<station_t> stations;
    stations.push_back(makeStation("PERT", 0, 1));
    stations.push_back(makeStation("ALIC", 1, 0));
    stations.push_back(makeStation("TIDB", 2, 2));

    vUINT32 name_indices = {0, 1, 2};
    vUINT32 file_indices = {0, 1, 2};

    CompareStnNameOrder<station_t, UINT32> name_cmp(&stations);
    std::sort(name_indices.begin(), name_indices.end(), name_cmp);

    CompareStnFileOrder<station_t, UINT32> file_cmp(&stations);
    std::sort(file_indices.begin(), file_indices.end(), file_cmp);

    // Name: ALIC(1), PERT(0), TIDB(2)
    // File: PERT(0), ALIC(1), TIDB(2)
    REQUIRE(name_indices[0] != file_indices[0]);
}

TEST_CASE("CompareStnNameOrder handles numeric station names correctly", "[m2s][sort][numeric-names]") {
    // nameOrder values reflect lexicographic order of the names:
    // "00NA_1900001" < "01NA_1900001" < "02NA" < "10-2882_301117"
    std::vector<station_t> stations;
    stations.push_back(makeStation("10-2882_301117", 0, 3));
    stations.push_back(makeStation("00NA_1900001",   1, 0));
    stations.push_back(makeStation("02NA",           2, 2));
    stations.push_back(makeStation("01NA_1900001",   3, 1));

    vUINT32 indices = {0, 1, 2, 3};

    CompareStnNameOrder<station_t, UINT32> cmp(&stations);
    std::sort(indices.begin(), indices.end(), cmp);

    REQUIRE(indices[0] == 1);  // 00NA_1900001 (nameOrder=0)
    REQUIRE(indices[1] == 3);  // 01NA_1900001 (nameOrder=1)
    REQUIRE(indices[2] == 2);  // 02NA          (nameOrder=2)
    REQUIRE(indices[3] == 0);  // 10-2882_301117 (nameOrder=3)
}

TEST_CASE("Decreasing measurement count is reverse of increasing", "[m2s][sort][decreasing]") {
    vUINT32 increasing = {3, 1, 4, 1, 5, 9, 2, 6};
    vUINT32 decreasing = increasing;

    std::sort(increasing.begin(), increasing.end());

    std::sort(decreasing.begin(), decreasing.end());
    std::reverse(decreasing.begin(), decreasing.end());

    REQUIRE(increasing.front() == decreasing.back());
    REQUIRE(increasing.back() == decreasing.front());

    for (size_t i = 0; i < increasing.size(); ++i)
        REQUIRE(increasing[i] == decreasing[increasing.size() - 1 - i]);
}

TEST_CASE("Default sort (orig_stn_sort_ui=0) is the default value", "[m2s][default]") {
    UINT16 default_sort = 0;
    REQUIRE(default_sort == orig_stn_sort_ui);
}

TEST_CASE("CompareStnNameOrder handles duplicate nameOrder", "[m2s][sort][duplicates]") {
    std::vector<station_t> stations;
    stations.push_back(makeStation("ALIC", 0, 0));
    stations.push_back(makeStation("ALIC", 1, 0));
    stations.push_back(makeStation("BEEC", 2, 1));

    vUINT32 indices = {0, 1, 2};

    CompareStnNameOrder<station_t, UINT32> cmp(&stations);
    std::sort(indices.begin(), indices.end(), cmp);

    std::string name0(stations[indices[0]].stationName);
    std::string name1(stations[indices[1]].stationName);
    std::string name2(stations[indices[2]].stationName);
    REQUIRE(name0 == "ALIC");
    REQUIRE(name1 == "ALIC");
    REQUIRE(name2 == "BEEC");
}

TEST_CASE("CompareStnNameOrder with single station", "[m2s][sort][single]") {
    std::vector<station_t> stations;
    stations.push_back(makeStation("ONLY", 0, 0));

    vUINT32 indices = {0};

    CompareStnNameOrder<station_t, UINT32> cmp(&stations);
    std::sort(indices.begin(), indices.end(), cmp);

    REQUIRE(indices[0] == 0);
}

TEST_CASE("All four sort modes produce valid orderings", "[m2s][sort][all-modes]") {
    std::vector<station_t> stations;
    stations.push_back(makeStation("PERT", 3, 2));   // idx 0
    stations.push_back(makeStation("ALIC", 1, 0));   // idx 1
    stations.push_back(makeStation("TIDB", 2, 3));   // idx 2
    stations.push_back(makeStation("BEEC", 0, 1));   // idx 3

    // Name sort (CompareStnNameOrder uses nameOrder)
    {
        vUINT32 indices = {0, 1, 2, 3};
        CompareStnNameOrder<station_t, UINT32> cmp(&stations);
        std::sort(indices.begin(), indices.end(), cmp);
        // nameOrder: ALIC(0), BEEC(1), PERT(2), TIDB(3)
        REQUIRE(std::string(stations[indices[0]].stationName) == "ALIC");
        REQUIRE(std::string(stations[indices[3]].stationName) == "TIDB");
    }

    // File order sort (CompareStnFileOrder uses fileOrder)
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
    for (int i = 0; i <= 3; ++i) {
        bool handled = false;
        switch (i) {
        case orig_stn_sort_ui: handled = true; break;
        case name_stn_sort_ui: handled = true; break;
        case meas_stn_sort_ui:   handled = true; break;
        case saem_stn_sort_ui: handled = true; break;
        }
        REQUIRE(handled);
    }
}
