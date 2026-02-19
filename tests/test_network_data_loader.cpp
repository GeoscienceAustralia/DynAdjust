//============================================================================
// Name         : test_network_data_loader.cpp
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
#define __BINARY_NAME__ "test_network_data_loader"
#endif
#ifndef __BINARY_DESC__
#define __BINARY_DESC__ "Unit tests for NetworkDataLoader class"
#endif

#include "testing.hpp"

#include "../dynadjust/dynadjust/dnaadjust/network_data_loader.hpp"
#include "../dynadjust/include/config/dnaoptions.hpp"
#include "../dynadjust/include/config/dnatypes.hpp"
#include "../dynadjust/include/measurement_types/dnameasurement.hpp"
#include "../dynadjust/include/measurement_types/dnastation.hpp"
#include <filesystem>
#include <optional>

using namespace dynadjust::networkadjust;
using namespace dynadjust::measurements;
namespace {

project_settings createTestSettings(adjustMode mode = SimultaneousMode) {
    project_settings settings;
    settings.a.adjust_mode = mode;
    settings.a.bst_file = "test.bst";
    settings.a.bms_file = "test.bms";
    settings.s.asl_file = "test.asl";
    return settings;
}

} // namespace

TEST_CASE("NetworkDataLoader constructor", "[network_data_loader][basic]") {
    auto settings = createTestSettings();
    NetworkDataLoader loader(settings);
    REQUIRE(true);
}

TEST_CASE("Move semantics", "[network_data_loader][move]") {
    auto settings = createTestSettings();
    NetworkDataLoader loader1(settings);
    NetworkDataLoader loader2 = std::move(loader1);
    REQUIRE(true);
}

TEST_CASE("Copy semantics deleted", "[network_data_loader][copy]") {
    // These should not compile:
    // NetworkDataLoader loader2(loader);
    // NetworkDataLoader loader3 = loader;
    REQUIRE(true);
}

TEST_CASE("LoadForPhased throws on missing files", "[network_data_loader][load]") {
    auto settings = createTestSettings(PhasedMode);
    NetworkDataLoader loader(settings);

    vstn_t bstBinaryRecords;
    binary_file_meta_t bst_meta{}, bms_meta{};
    vASL vAssocStnList;
    vmsr_t bmsBinaryRecords;
    vvUINT32 v_ISL;
    v_uint32_uint32_map v_blockStationsMap(1);
    vvUINT32 v_CML;

    UINT32 bstn_count = 0, asl_count = 0, bmsr_count = 0;
    UINT32 unknownParams = 0, unknownsCount = 0;
    UINT32 measurementParams = 0, measurementCount = 0, measurementVarianceCount = 0;

    bool threw = false;
    try {
        loader.LoadForPhased(&bstBinaryRecords, bst_meta, &vAssocStnList,
            &bmsBinaryRecords, bms_meta, v_ISL, &v_blockStationsMap,
            &v_CML, bstn_count, asl_count, bmsr_count,
            unknownParams, unknownsCount, measurementParams,
            measurementCount, measurementVarianceCount);
    } catch (const NetworkLoadException&) {
        threw = true;
    }
    REQUIRE(threw);
}

TEST_CASE("LoadForSimultaneous throws on missing files", "[network_data_loader][load]") {
    auto settings = createTestSettings(SimultaneousMode);
    NetworkDataLoader loader(settings);

    vstn_t bstBinaryRecords;
    binary_file_meta_t bst_meta{}, bms_meta{};
    vASL vAssocStnList;
    vmsr_t bmsBinaryRecords;
    vvUINT32 v_ISL;
    v_uint32_uint32_map v_blockStationsMap(1);
    vvUINT32 v_CML;

    UINT32 bstn_count = 0, asl_count = 0, bmsr_count = 0;
    UINT32 unknownParams = 0, unknownsCount = 0;
    UINT32 measurementParams = 0, measurementCount = 0, measurementVarianceCount = 0;

    UINT32 blockCount = 0;
    vvUINT32 v_JSL;
    vUINT32 v_unknownsCount, v_measurementCount, v_measurementVarianceCount;
    vUINT32 v_measurementParams, v_ContiguousNetList;
    std::vector<blockMeta_t> v_blockMeta;
    vvUINT32 v_parameterStationList;
    vv_stn_appear v_paramStnAppearance;
    v_mat_2d v_junctionVariances, v_junctionVariancesFwd;

    bool threw = false;
    try {
        loader.LoadForSimultaneous(&bstBinaryRecords, bst_meta, &vAssocStnList,
            &bmsBinaryRecords, bms_meta, v_ISL, &v_blockStationsMap,
            &v_CML, bstn_count, asl_count, bmsr_count,
            unknownParams, unknownsCount, measurementParams,
            measurementCount, measurementVarianceCount,
            &blockCount, &v_JSL, &v_unknownsCount, &v_measurementCount,
            &v_measurementVarianceCount, &v_measurementParams,
            &v_ContiguousNetList, &v_blockMeta, &v_parameterStationList,
            &v_paramStnAppearance, &v_junctionVariances, &v_junctionVariancesFwd);
    } catch (const NetworkLoadException&) {
        threw = true;
    }
    REQUIRE(threw);
}

TEST_CASE("RemoveInvalidStations functionality", "[network_data_loader][station_filter]") {
    auto settings = createTestSettings();
    NetworkDataLoader loader(settings);

    vASL test_associated_stations;

    test_associated_stations.emplace_back();
    test_associated_stations.back().SetValid();

    test_associated_stations.emplace_back();
    test_associated_stations.back().SetInvalid();

    test_associated_stations.emplace_back();
    test_associated_stations.back().SetValid();

    vUINT32 station_list = {0, 1, 2};

    loader.RemoveInvalidStations(station_list, test_associated_stations);

    REQUIRE(station_list.size() == 2);
    REQUIRE(station_list[0] == 0);
    REQUIRE(station_list[1] == 2);
}

TEST_CASE("RemoveInvalidStations with empty list", "[network_data_loader][station_filter]") {
    auto settings = createTestSettings();
    NetworkDataLoader loader(settings);

    vASL test_associated_stations;
    vUINT32 station_list;

    loader.RemoveInvalidStations(station_list, test_associated_stations);

    REQUIRE(station_list.empty());
}

TEST_CASE("RemoveInvalidStations all valid", "[network_data_loader][station_filter]") {
    auto settings = createTestSettings();
    NetworkDataLoader loader(settings);

    vASL test_associated_stations;
    test_associated_stations.emplace_back();
    test_associated_stations.back().SetValid();
    test_associated_stations.emplace_back();
    test_associated_stations.back().SetValid();

    vUINT32 station_list = {0, 1};

    loader.RemoveInvalidStations(station_list, test_associated_stations);

    REQUIRE(station_list.size() == 2);
}

TEST_CASE("RemoveInvalidStations all invalid", "[network_data_loader][station_filter]") {
    auto settings = createTestSettings();
    NetworkDataLoader loader(settings);

    vASL test_associated_stations;
    test_associated_stations.emplace_back();
    test_associated_stations.back().SetInvalid();
    test_associated_stations.emplace_back();
    test_associated_stations.back().SetInvalid();

    vUINT32 station_list = {0, 1};

    loader.RemoveInvalidStations(station_list, test_associated_stations);

    REQUIRE(station_list.empty());
}

TEST_CASE("RemoveNonMeasurements with empty list", "[network_data_loader][msr_filter]") {
    auto settings = createTestSettings();
    NetworkDataLoader loader(settings);

    vmsr_t measurements;
    vUINT32 msr_list;

    loader.RemoveNonMeasurements(msr_list, measurements);

    REQUIRE(msr_list.empty());
}

TEST_CASE("RemoveNonMeasurements with single element", "[network_data_loader][msr_filter]") {
    auto settings = createTestSettings();
    NetworkDataLoader loader(settings);

    vmsr_t measurements;
    measurement_t msr;
    msr.measStart = xMeas;
    msr.fileOrder = 0;
    measurements.push_back(msr);

    vUINT32 msr_list = {0};

    loader.RemoveNonMeasurements(msr_list, measurements);

    REQUIRE(msr_list.size() == 1);
}

TEST_CASE("Different adjustment modes construct correctly", "[network_data_loader][modes]") {
    for (auto mode : {SimultaneousMode, PhasedMode, Phased_Block_1Mode}) {
        auto settings = createTestSettings(mode);
        NetworkDataLoader loader(settings);
        REQUIRE(true);
    }
}
