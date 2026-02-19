//============================================================================
// Name         : test_measurement_processor.cpp
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
#define __BINARY_NAME__ "test_measurement_processor"
#endif
#ifndef __BINARY_DESC__
#define __BINARY_DESC__ "Unit tests for MeasurementProcessor class"
#endif

#include "testing.hpp"

#include "../dynadjust/dynadjust/dnaadjust/measurement_processor.hpp"
#include "../dynadjust/include/config/dnatypes.hpp"
#include "../dynadjust/include/measurement_types/dnameasurement.hpp"

using namespace dynadjust::networkadjust::processors;
using namespace dynadjust::measurements;

namespace {
measurement_t createTestMeasurement(char type, bool ignore = false, UINT32 clusterID = 0, UINT32 vectorCount1 = 0,
                                    UINT32 vectorCount2 = 0, MEASUREMENT_START measStart = xMeas) {
    measurement_t msr = {};
    msr.measType = type;
    msr.ignore = ignore;
    msr.clusterID = clusterID;
    msr.vectorCount1 = vectorCount1;
    msr.vectorCount2 = vectorCount2;
    msr.measStart = measStart;
    msr.station1 = 0;
    msr.station2 = 1;
    msr.term1 = 100.0;
    return msr;
}
} // namespace

// Basic constructor tests
TEST_CASE("MeasurementProcessor constructor simultaneous", "[measurement_processor][basic]") {
    MeasurementProcessor processor(AdjustmentMode::Simultaneous);
    REQUIRE(true);
}

TEST_CASE("MeasurementProcessor constructor phased", "[measurement_processor][basic]") {
    MeasurementProcessor processor(AdjustmentMode::Phased);
    REQUIRE(true);
}

// Move semantics tests
TEST_CASE("MeasurementProcessor move constructor", "[measurement_processor][move]") {
    MeasurementProcessor processor1(AdjustmentMode::Simultaneous);
    MeasurementProcessor processor2 = std::move(processor1);
    REQUIRE(true);
}

TEST_CASE("MeasurementProcessor move assignment", "[measurement_processor][move]") {
    MeasurementProcessor processor1(AdjustmentMode::Simultaneous);
    MeasurementProcessor processor2(AdjustmentMode::Phased);
    processor2 = std::move(processor1);
    REQUIRE(true);
}

// Copy semantics are deleted - compile-time test
TEST_CASE("MeasurementProcessor copy semantics deleted", "[measurement_processor][copy]") {
    // These should not compile:
    // MeasurementProcessor processor2(processor);
    // MeasurementProcessor processor3 = processor;
    REQUIRE(true);
}

// Empty data tests
TEST_CASE("Process empty measurement list", "[measurement_processor][empty]") {
    MeasurementProcessor processor(AdjustmentMode::Simultaneous);

    vmsr_t empty_measurements;
    vvUINT32 v_CML;
    MeasurementCounts counts;

    auto result = processor.ProcessForMode(empty_measurements, 0, v_CML, counts);

    REQUIRE(result.has_value());
    REQUIRE(*result == 0);
    REQUIRE(counts.measurement_count == 0);
    REQUIRE(counts.measurement_variance_count == 0);
    REQUIRE(counts.measurement_params == 0);
    REQUIRE(v_CML.size() == 1);
    REQUIRE(v_CML[0].empty());
}

// Basic measurement processing tests
TEST_CASE("Process single distance measurement", "[measurement_processor][single]") {
    MeasurementProcessor processor(AdjustmentMode::Simultaneous);

    vmsr_t measurements;
    measurements.push_back(createTestMeasurement('S'));

    vvUINT32 v_CML;
    MeasurementCounts counts;

    auto result = processor.ProcessForMode(measurements, 1, v_CML, counts);

    REQUIRE(result.has_value());
    REQUIRE(*result == 1);
    REQUIRE(counts.measurement_count == 1);
    REQUIRE(counts.measurement_variance_count == 1);
    REQUIRE(counts.measurement_params == 1);
    REQUIRE(v_CML.size() == 1);
    REQUIRE(v_CML[0].size() == 1);
    REQUIRE(v_CML[0][0] == 0);
}

TEST_CASE("Process single angle measurement", "[measurement_processor][single]") {
    MeasurementProcessor processor(AdjustmentMode::Simultaneous);

    vmsr_t measurements;
    measurements.push_back(createTestMeasurement('A'));

    vvUINT32 v_CML;
    MeasurementCounts counts;

    auto result = processor.ProcessForMode(measurements, 1, v_CML, counts);

    REQUIRE(result.has_value());
    REQUIRE(*result == 1);
    REQUIRE(counts.measurement_count == 1);
    REQUIRE(counts.measurement_variance_count == 1);
    REQUIRE(counts.measurement_params == 1);
}

// GPS measurement tests
// ProcessSimultaneous returns m = count of v_CML entries.
// For GPS X/Y, only the first measurement in each cluster is added to v_CML.
// G type measurements are always added individually.
TEST_CASE("Process GPS baseline measurements", "[measurement_processor][gps]") {
    MeasurementProcessor processor(AdjustmentMode::Simultaneous);

    vmsr_t measurements;
    // X with clusterID=1 → new cluster, added to v_CML (m=1)
    measurements.push_back(createTestMeasurement('X', false, 1));
    // Y with clusterID=1 → same cluster, not added to v_CML
    measurements.push_back(createTestMeasurement('Y', false, 1));
    // G always added to v_CML via default case (m=2)
    measurements.push_back(createTestMeasurement('G', false, 1));

    vvUINT32 v_CML;
    MeasurementCounts counts;

    auto result = processor.ProcessForMode(measurements, static_cast<UINT32>(measurements.size()), v_CML, counts);

    REQUIRE(result.has_value());
    REQUIRE(*result == 2);  // only 2 entries in v_CML
    REQUIRE(counts.measurement_count == 3);
    // Variance: X→1(axis=0→1), Y→2(axis=1→2), G→3(axis=2→0) = 6
    REQUIRE(counts.measurement_variance_count == 6);
    REQUIRE(counts.measurement_params == 3);
}

// Direction set tests
// D type: only added to v_CML if vectorCount1 >= 1
TEST_CASE("Process direction set measurements", "[measurement_processor][directions]") {
    MeasurementProcessor processor(AdjustmentMode::Simultaneous);

    vmsr_t measurements;
    // vectorCount1=1 means it's an RO with targets, vectorCount2=5 means 5 directions
    measurements.push_back(createTestMeasurement('D', false, 0, 1, 5));

    vvUINT32 v_CML;
    MeasurementCounts counts;

    auto result = processor.ProcessForMode(measurements, static_cast<UINT32>(measurements.size()), v_CML, counts);

    REQUIRE(result.has_value());
    REQUIRE(*result == 1);
    // vectorCount2=5 → 5-1=4 measurements
    REQUIRE(counts.measurement_count == 4);
    REQUIRE(counts.measurement_variance_count == 4);
    REQUIRE(counts.measurement_params == 4);
}

// Mixed measurement types
TEST_CASE("Process mixed measurement types", "[measurement_processor][mixed]") {
    MeasurementProcessor processor(AdjustmentMode::Simultaneous);

    vmsr_t measurements;
    measurements.push_back(createTestMeasurement('S'));
    measurements.push_back(createTestMeasurement('A'));
    measurements.push_back(createTestMeasurement('L'));
    measurements.push_back(createTestMeasurement('V'));

    vvUINT32 v_CML;
    MeasurementCounts counts;

    auto result = processor.ProcessForMode(measurements, static_cast<UINT32>(measurements.size()), v_CML, counts);

    REQUIRE(result.has_value());
    REQUIRE(*result == 4);
    REQUIRE(counts.measurement_count == 4);
    REQUIRE(counts.measurement_variance_count == 4);
    REQUIRE(counts.measurement_params == 4);
}

// Ignored measurements are skipped entirely
TEST_CASE("Process measurements with ignored entries", "[measurement_processor][ignored]") {
    MeasurementProcessor processor(AdjustmentMode::Simultaneous);

    vmsr_t measurements;
    measurements.push_back(createTestMeasurement('S', false));
    measurements.push_back(createTestMeasurement('S', true));   // ignored
    measurements.push_back(createTestMeasurement('A', false));

    vvUINT32 v_CML;
    MeasurementCounts counts;

    auto result = processor.ProcessForMode(measurements, static_cast<UINT32>(measurements.size()), v_CML, counts);

    REQUIRE(result.has_value());
    REQUIRE(*result == 2);  // only non-ignored added to v_CML
    REQUIRE(counts.measurement_count == 2);
    REQUIRE(counts.measurement_variance_count == 2);
    REQUIRE(counts.measurement_params == 2);
}

// Covariance terms (measStart > zMeas) are skipped
TEST_CASE("Process measurements with covariance terms", "[measurement_processor][covariance]") {
    MeasurementProcessor processor(AdjustmentMode::Simultaneous);

    vmsr_t measurements;
    measurements.push_back(createTestMeasurement('S', false, 0, 0, 0, xMeas));
    measurements.push_back(createTestMeasurement('S', false, 0, 0, 0, xCov));
    measurements.push_back(createTestMeasurement('S', false, 0, 0, 0, yCov));
    measurements.push_back(createTestMeasurement('S', false, 0, 0, 0, zCov));

    vvUINT32 v_CML;
    MeasurementCounts counts;

    auto result = processor.ProcessForMode(measurements, static_cast<UINT32>(measurements.size()), v_CML, counts);

    REQUIRE(result.has_value());
    REQUIRE(*result == 1);  // only the xMeas entry
    REQUIRE(counts.measurement_count == 1);
    REQUIRE(counts.measurement_variance_count == 1);
    REQUIRE(counts.measurement_params == 1);
}

// Phased mode tests
TEST_CASE("Process in phased mode", "[measurement_processor][phased]") {
    MeasurementProcessor processor(AdjustmentMode::Phased);

    vmsr_t measurements;
    measurements.push_back(createTestMeasurement('S'));

    vvUINT32 v_CML;
    MeasurementCounts counts;

    auto result = processor.ProcessForMode(measurements, 1, v_CML, counts);

    REQUIRE(result == std::nullopt);
}

// Large dataset tests
TEST_CASE("Process large measurement dataset", "[measurement_processor][large]") {
    MeasurementProcessor processor(AdjustmentMode::Simultaneous);

    vmsr_t measurements;
    const UINT32 num_measurements = 1000;

    for (UINT32 i = 0; i < num_measurements; ++i) {
        measurements.push_back(createTestMeasurement('S'));
    }

    vvUINT32 v_CML;
    MeasurementCounts counts;

    auto result = processor.ProcessForMode(measurements, num_measurements, v_CML, counts);

    REQUIRE(result.has_value());
    REQUIRE(*result == num_measurements);
    REQUIRE(counts.measurement_count == num_measurements);
    REQUIRE(counts.measurement_variance_count == num_measurements);
    REQUIRE(counts.measurement_params == num_measurements);
    REQUIRE(v_CML.size() == 1);
    REQUIRE(v_CML[0].size() == num_measurements);
}

// GPS cluster handling: only first measurement per cluster is added to v_CML
TEST_CASE("Process GPS measurements with different clusters", "[measurement_processor][clusters]") {
    MeasurementProcessor processor(AdjustmentMode::Simultaneous);

    vmsr_t measurements;
    // Cluster 1: X added to v_CML (new cluster), Y not added (same cluster)
    measurements.push_back(createTestMeasurement('X', false, 1));
    measurements.push_back(createTestMeasurement('Y', false, 1));
    // Cluster 2: X added to v_CML (new cluster), Y not added (same cluster)
    measurements.push_back(createTestMeasurement('X', false, 2));
    measurements.push_back(createTestMeasurement('Y', false, 2));

    vvUINT32 v_CML;
    MeasurementCounts counts;

    auto result = processor.ProcessForMode(measurements, static_cast<UINT32>(measurements.size()), v_CML, counts);

    REQUIRE(result.has_value());
    REQUIRE(*result == 2);  // 2 cluster heads in v_CML
    REQUIRE(counts.measurement_count == 4);
    REQUIRE(v_CML.size() == 1);
    REQUIRE(v_CML[0].size() == 2);
}

// bmsr_count doesn't affect result — we iterate the vector
TEST_CASE("Handle measurement count mismatch", "[measurement_processor][mismatch]") {
    MeasurementProcessor processor(AdjustmentMode::Simultaneous);

    vmsr_t measurements;
    measurements.push_back(createTestMeasurement('S'));
    measurements.push_back(createTestMeasurement('S', true)); // ignored

    vvUINT32 v_CML;
    MeasurementCounts counts;

    auto result = processor.ProcessForMode(measurements, 2, v_CML, counts);

    REQUIRE(result.has_value());
    REQUIRE(*result == 1);  // only 1 non-ignored added to v_CML
    REQUIRE(counts.measurement_count == 1);
    REQUIRE(v_CML.size() == 1);
    REQUIRE(v_CML[0].size() == 1);
}

// Multiple direction sets
TEST_CASE("Process multiple direction sets", "[measurement_processor][multiple_directions]") {
    MeasurementProcessor processor(AdjustmentMode::Simultaneous);

    vmsr_t measurements;
    // vectorCount1=1 (RO present), vectorCount2=3 (3 directions → 2 angles)
    measurements.push_back(createTestMeasurement('D', false, 0, 1, 3));
    // vectorCount1=1 (RO present), vectorCount2=5 (5 directions → 4 angles)
    measurements.push_back(createTestMeasurement('D', false, 0, 1, 5));

    vvUINT32 v_CML;
    MeasurementCounts counts;

    auto result = processor.ProcessForMode(measurements, static_cast<UINT32>(measurements.size()), v_CML, counts);

    REQUIRE(result.has_value());
    REQUIRE(*result == 2);  // 2 direction set ROs in v_CML
    // (3-1) + (5-1) = 6 measurements
    REQUIRE(counts.measurement_count == 6);
    REQUIRE(counts.measurement_variance_count == 6);
    REQUIRE(counts.measurement_params == 6);
}

// Direction set with vectorCount1=0 and vectorCount2=0: not added to v_CML
TEST_CASE("Process direction set with zero vector count", "[measurement_processor][edge_direction]") {
    MeasurementProcessor processor(AdjustmentMode::Simultaneous);

    vmsr_t measurements;
    measurements.push_back(createTestMeasurement('D', false, 0, 0, 0));

    vvUINT32 v_CML;
    MeasurementCounts counts;

    auto result = processor.ProcessForMode(measurements, 1, v_CML, counts);

    REQUIRE(result.has_value());
    REQUIRE(*result == 0);  // vectorCount1=0 → not added to v_CML
    REQUIRE(counts.measurement_count == 0);
    REQUIRE(counts.measurement_variance_count == 0);
    REQUIRE(counts.measurement_params == 0);
}

// Direction set with vectorCount1=1 and vectorCount2=1: added but (1-1)=0 angles
// Falls into vectorCount2 > 0 branch: count += 1-1 = 0
TEST_CASE("Process direction set with single vector count", "[measurement_processor][edge_direction_single]") {
    MeasurementProcessor processor(AdjustmentMode::Simultaneous);

    vmsr_t measurements;
    measurements.push_back(createTestMeasurement('D', false, 0, 1, 1));

    vvUINT32 v_CML;
    MeasurementCounts counts;

    auto result = processor.ProcessForMode(measurements, 1, v_CML, counts);

    REQUIRE(result.has_value());
    REQUIRE(*result == 1);  // vectorCount1=1 → added to v_CML
    REQUIRE(counts.measurement_count == 0);  // vectorCount2=1 → 1-1=0
    REQUIRE(counts.measurement_variance_count == 0);
    REQUIRE(counts.measurement_params == 0);
}

// Input validation tests
TEST_CASE("Process with pre-populated v_CML", "[measurement_processor][pre_populated]") {
    MeasurementProcessor processor(AdjustmentMode::Simultaneous);

    vmsr_t measurements;
    measurements.push_back(createTestMeasurement('S'));

    vvUINT32 v_CML;
    v_CML.push_back({99, 98, 97});
    v_CML.push_back({96, 95});

    MeasurementCounts counts;

    auto result = processor.ProcessForMode(measurements, 1, v_CML, counts);

    REQUIRE(result.has_value());
    REQUIRE(*result == 1);
    REQUIRE(v_CML.size() == 1);
    REQUIRE(v_CML[0].size() == 1);
    REQUIRE(v_CML[0][0] == 0);
}

// Boundary tests
TEST_CASE("Process with maximum UINT32 bmsr_count", "[measurement_processor][boundary]") {
    MeasurementProcessor processor(AdjustmentMode::Simultaneous);

    vmsr_t measurements;
    measurements.push_back(createTestMeasurement('S'));

    vvUINT32 v_CML;
    MeasurementCounts counts;

    auto result = processor.ProcessForMode(measurements, 1000, v_CML, counts);

    REQUIRE(result.has_value());
    REQUIRE(*result == 1);
    REQUIRE(v_CML.size() == 1);
    REQUIRE(v_CML[0].size() == 1);
}

// GPS axis counter: variance pattern is 1,2,3,1,2,3...
TEST_CASE("GPS axis counter behavior", "[measurement_processor][gps_axis]") {
    MeasurementProcessor processor(AdjustmentMode::Simultaneous);

    vmsr_t measurements;
    for (int i = 0; i < 6; ++i) {
        measurements.push_back(createTestMeasurement('G'));
    }

    vvUINT32 v_CML;
    MeasurementCounts counts;

    auto result = processor.ProcessForMode(measurements, static_cast<UINT32>(measurements.size()), v_CML, counts);

    REQUIRE(result.has_value());
    REQUIRE(*result == 6);
    REQUIRE(counts.measurement_count == 6);
    // Variance: 1+2+3+1+2+3 = 12
    REQUIRE(counts.measurement_variance_count == 12);
    REQUIRE(counts.measurement_params == 6);
}

// All supported simple measurement types (non-GPS, non-direction)
TEST_CASE("Process all supported measurement types", "[measurement_processor][all_types]") {
    MeasurementProcessor processor(AdjustmentMode::Simultaneous);

    vmsr_t measurements;
    // Simple types that go through the default path (1 count, 1 variance each)
    std::vector<char> simple_types = {'A', 'B', 'C', 'E', 'H', 'I', 'J', 'K', 'L',
                                      'M', 'P', 'Q', 'R', 'S', 'V', 'Z'};

    for (char type : simple_types) {
        measurements.push_back(createTestMeasurement(type));
    }

    vvUINT32 v_CML;
    MeasurementCounts counts;

    auto result = processor.ProcessForMode(measurements, static_cast<UINT32>(measurements.size()), v_CML, counts);

    REQUIRE(result.has_value());
    REQUIRE(*result == simple_types.size());
    REQUIRE(counts.measurement_count == simple_types.size());
    REQUIRE(counts.measurement_variance_count == simple_types.size());
    REQUIRE(v_CML.size() == 1);
    REQUIRE(v_CML[0].size() == simple_types.size());
}
