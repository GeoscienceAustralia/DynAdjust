//============================================================================
// Name         : test_dnaparser_jsonl.cpp
// Author       : Dale Roberts <dale.o.roberts@gmail.com>
// Copyright    : Copyright 2026 Geoscience Australia
//
//                Licensed under the Apache License, Version 2.0 (the "License");
//                you may not use this file except in compliance with the License.
//                You may obtain a copy of the License at
//
//                http://www.apache.org/licenses/LICENSE-2.0
//
// Description  : JSONL parser regression tests.
//============================================================================

#define TESTING_MAIN
#ifndef __BINARY_NAME__
#define __BINARY_NAME__ "test_dnaparser_jsonl"
#endif
#ifndef __BINARY_DESC__
#define __BINARY_DESC__ "Regression tests for JSONL parser"
#endif

#include "testing.hpp"

#include <filesystem>
#include <fstream>
#include <string>
#include <system_error>

#include <dynadjust/dnaimport/dnaparser_jsonl.hpp>
#include <include/measurement_types/dnadirection.hpp>
#include <include/measurement_types/dnagpsbaseline.hpp>
#include <include/measurement_types/dnagpspoint.hpp>
#include <include/measurement_types/dnameasurement.hpp>
#include <include/measurement_types/dnastntally.hpp>

using dynadjust::import::JsonlParseContext;
using dynadjust::import::ParseJsonlFile;
using dynadjust::measurements::MsrTally;
using dynadjust::measurements::StnTally;
using dynadjust::measurements::vdnaMsrPtr;
using dynadjust::measurements::vdnaStnPtr;

MsrTally g_parsemsr_tally;
StnTally g_parsestn_tally;
UINT32 g_fileOrder = 0;

namespace {

struct TempJsonlFile {
  explicit TempJsonlFile(const char* name)
      : path(std::filesystem::temp_directory_path() / name) {}

  ~TempJsonlFile() {
    std::error_code ec;
    std::filesystem::remove(path, ec);
  }

  std::filesystem::path path;
};

JsonlParseContext DefaultContext() {
  JsonlParseContext ctx{};
  ctx.default_frame = "GDA2020";
  ctx.default_epoch = "01.01.2020";
  ctx.first_file = true;
  ctx.user_supplied_frame = false;
  ctx.user_supplied_epoch = false;
  ctx.override_input_frame = false;
  ctx.prefer_single_x_as_g = false;
  ctx.cluster_id = 0;
  return ctx;
}

}  // namespace

TEST_CASE("JSONL parser ignores DnaStatistics adjustment metadata records") {
  TempJsonlFile tmp("dynadjust-jsonl-statistics-regression.jsonl");
  {
    std::ofstream out(tmp.path);
    out << "{\"DnaAdjustmentReport\":{\"type\":\"Adjustment\","
           "\"report\":\"adj\",\"software\":\"DynAdjust 1.4.0\","
           "\"referenceframe\":\"GDA2020\",\"epoch\":\"01.01.2020\"}}\n";
    out << "{\"DnaStatistics\":{\"iteration\":1,\"dof\":288,"
           "\"sigma_zero\":1.16,\"chisq_test\":\"passed\"}}\n";
    out << "{\"DnaStation\":{\"Name\":\"STN1\",\"Constraints\":\"FFF\","
           "\"Type\":\"LLH\",\"StationCoord\":{\"Name\":\"STN1\","
           "\"XAxis\":-36.3348253511,\"YAxis\":145.5741006918,"
           "\"Height\":172.1735},\"Description\":\"Parser regression\"}}\n";
  }

  vdnaStnPtr stations;
  vdnaMsrPtr measurements;
  JsonlParseContext ctx = DefaultContext();

  ParseJsonlFile(tmp.path.string(), &stations, &measurements, ctx);

  REQUIRE(ctx.stn_count == 1);
  REQUIRE(ctx.msr_count == 0);
  REQUIRE(stations.size() == 1);
  REQUIRE(measurements.empty());
}

TEST_CASE("JSONL parser rejects invalid header reference frame") {
  TempJsonlFile tmp("dynadjust-jsonl-invalid-header-frame.jsonl");
  {
    std::ofstream out(tmp.path);
    out << "{\"DnaXmlFormat\":{\"type\":\"Measurement File\","
           "\"referenceframe\":\"GDA202O\",\"epoch\":\"01.01.2020\"}}\n";
  }

  vdnaStnPtr stations;
  vdnaMsrPtr measurements;
  JsonlParseContext ctx = DefaultContext();

  bool threw = false;
  std::string message;
  try {
    ParseJsonlFile(tmp.path.string(), &stations, &measurements, ctx);
  } catch (const std::exception& e) {
    threw = true;
    message = e.what();
  }

  REQUIRE(threw);
  REQUIRE(message.find("JSONL line 1") != std::string::npos);
}

TEST_CASE("JSONL parser re-ingests adjustment measurement records") {
  TempJsonlFile tmp("dynadjust-jsonl-measurement-regression.jsonl");
  {
    std::ofstream out(tmp.path);
    out << "{\"DnaAdjustmentReport\":{\"type\":\"Adjustment\","
           "\"report\":\"adj\",\"software\":\"DynAdjust 1.4.0\","
           "\"referenceframe\":\"GDA2020\",\"epoch\":\"01.01.2020\"}}\n";
    out << "{\"DnaMeasurement\":{\"Type\":\"L\",\"First\":\"STN1\","
           "\"Second\":\"STN2\",\"Value\":1234.565,\"StdDev\":0.004,"
           "\"Adjusted\":1234.567,\"Correction\":0.002,"
           "\"AdjustedPrecision\":0.001,\"ResidualPrecision\":0.001,"
           "\"NStat\":0.8,\"TStat\":1.1,\"PelzerRel\":0.92}}\n";
    out << "{\"DnaMeasurement\":{\"Type\":\"D\",\"First\":\"STN1\","
           "\"Second\":\"STN2\",\"Value\":45.0,\"StdDev\":0.01,"
           "\"Total\":1,\"Directions\":[{\"Target\":\"STN2\","
           "\"Value\":45.0,\"StdDev\":0.01}],\"Adjusted\":45.1,"
           "\"Correction\":0.1}}\n";
    out << "{\"DnaMeasurement\":{\"Type\":\"G\",\"First\":\"STN1\","
           "\"Second\":\"STN2\",\"Total\":1,\"GPSBaseline\":[{"
           "\"First\":\"STN1\",\"Second\":\"STN2\",\"X\":1.099,"
           "\"Y\":2.198,\"Z\":3.297,\"SigmaXX\":0.0001,"
           "\"SigmaXY\":0.00001,\"SigmaXZ\":0.00002,"
           "\"SigmaYY\":0.0002,\"SigmaYZ\":0.00003,"
           "\"SigmaZZ\":0.0003}],\"Adjusted\":{\"X\":1.1,"
           "\"Y\":2.2,\"Z\":3.3}}}\n";
    out << "{\"DnaMeasurement\":{\"Type\":\"Y\",\"First\":\"STN1\","
           "\"Total\":1,\"Coords\":\"XYZ\",\"Clusterpoint\":[{"
           "\"First\":\"STN1\",\"X\":4.099,\"Y\":5.198,\"Z\":6.297,"
           "\"SigmaXX\":0.0004,\"SigmaXY\":0.00004,"
           "\"SigmaXZ\":0.00005,\"SigmaYY\":0.0005,"
           "\"SigmaYZ\":0.00006,\"SigmaZZ\":0.0006}],"
           "\"Adjusted\":{\"X\":4.1,\"Y\":5.2,\"Z\":6.3}}}\n";
  }

  vdnaStnPtr stations;
  vdnaMsrPtr measurements;
  JsonlParseContext ctx = DefaultContext();

  ParseJsonlFile(tmp.path.string(), &stations, &measurements, ctx);

  REQUIRE(ctx.stn_count == 0);
  REQUIRE(ctx.msr_count == 8);
  REQUIRE(stations.empty());
  REQUIRE(measurements.size() == 4);
}

TEST_CASE("JSONL parser preserves multi-item GNSS cluster covariances") {
  TempJsonlFile tmp("dynadjust-jsonl-gnss-cluster-covariance.jsonl");
  {
    std::ofstream out(tmp.path);
    out << "{\"DnaAdjustmentReport\":{\"type\":\"Adjustment\","
           "\"report\":\"adj\",\"software\":\"DynAdjust 1.4.0\","
           "\"referenceframe\":\"GDA2020\",\"epoch\":\"01.01.2020\"}}\n";
    out << "{\"DnaMeasurement\":{\"Type\":\"X\",\"First\":\"STN1\","
           "\"Second\":\"STN2\",\"Total\":2,\"GPSBaseline\":[{"
           "\"First\":\"STN1\",\"Second\":\"STN2\",\"X\":1.0,"
           "\"Y\":2.0,\"Z\":3.0,\"SigmaXX\":0.0001,"
           "\"SigmaXY\":0.00001,\"SigmaXZ\":0.00002,"
           "\"SigmaYY\":0.0002,\"SigmaYZ\":0.00003,"
           "\"SigmaZZ\":0.0003,\"GPSCovariance\":[{"
           "\"m11\":0.000001,\"m12\":0.000002,\"m13\":0.000003,"
           "\"m21\":0.000004,\"m22\":0.000005,\"m23\":0.000006,"
           "\"m31\":0.000007,\"m32\":0.000008,\"m33\":0.000009}]},"
           "{\"First\":\"STN1\",\"Second\":\"STN3\",\"X\":4.0,"
           "\"Y\":5.0,\"Z\":6.0,\"SigmaXX\":0.0004,"
           "\"SigmaXY\":0.00004,\"SigmaXZ\":0.00005,"
           "\"SigmaYY\":0.0005,\"SigmaYZ\":0.00006,"
           "\"SigmaZZ\":0.0006}]}}\n";
    out << "{\"DnaMeasurement\":{\"Type\":\"Y\",\"First\":\"STN1\","
           "\"Total\":2,\"Coords\":\"XYZ\",\"Clusterpoint\":[{"
           "\"First\":\"STN1\",\"X\":7.0,\"Y\":8.0,\"Z\":9.0,"
           "\"SigmaXX\":0.0007,\"SigmaXY\":0.00007,"
           "\"SigmaXZ\":0.00008,\"SigmaYY\":0.0008,"
           "\"SigmaYZ\":0.00009,\"SigmaZZ\":0.0009,"
           "\"PointCovariance\":[{\"m11\":0.000011,"
           "\"m12\":0.000012,\"m13\":0.000013,"
           "\"m21\":0.000014,\"m22\":0.000015,"
           "\"m23\":0.000016,\"m31\":0.000017,"
           "\"m32\":0.000018,\"m33\":0.000019}]},"
           "{\"First\":\"STN2\",\"X\":10.0,\"Y\":11.0,\"Z\":12.0,"
           "\"SigmaXX\":0.0010,\"SigmaXY\":0.00010,"
           "\"SigmaXZ\":0.00011,\"SigmaYY\":0.0011,"
           "\"SigmaYZ\":0.00012,\"SigmaZZ\":0.0012}]}}\n";
  }

  vdnaStnPtr stations;
  vdnaMsrPtr measurements;
  JsonlParseContext ctx = DefaultContext();

  ParseJsonlFile(tmp.path.string(), &stations, &measurements, ctx);

  REQUIRE(ctx.stn_count == 0);
  REQUIRE(ctx.msr_count == 12);
  REQUIRE(stations.empty());
  REQUIRE(measurements.size() == 2);

  auto* baselines = measurements[0]->GetBaselines_ptr();
  REQUIRE(baselines != nullptr);
  REQUIRE(baselines->size() == 2);
  REQUIRE(baselines->at(0).GetCovariances_ptr()->size() == 1);
  REQUIRE(baselines->at(1).GetCovariances_ptr()->empty());
  REQUIRE(baselines->at(0).GetCovariances_ptr()->at(0).GetM11() == 0.000001);

  auto* points = measurements[1]->GetPoints_ptr();
  REQUIRE(points != nullptr);
  REQUIRE(points->size() == 2);
  REQUIRE(points->at(0).GetCovariances_ptr()->size() == 1);
  REQUIRE(points->at(1).GetCovariances_ptr()->empty());
  REQUIRE(points->at(0).GetCovariances_ptr()->at(0).GetM11() == 0.000011);
}

TEST_CASE("JSONL parser preserves measurement observation epochs") {
  TempJsonlFile tmp("dynadjust-jsonl-observation-epoch.jsonl");
  {
    std::ofstream out(tmp.path);
    out << "{\"DnaAdjustmentReport\":{\"type\":\"Adjustment\","
           "\"report\":\"adj\",\"software\":\"DynAdjust 1.4.0\","
           "\"referenceframe\":\"GDA2020\",\"epoch\":\"01.01.2020\"}}\n";
    out << "{\"DnaMeasurement\":{\"Type\":\"L\",\"First\":\"STN1\","
           "\"Second\":\"STN2\",\"Epoch\":\"01.01.2020\","
           "\"EpochOfObservation\":\"15.07.2021\","
           "\"Value\":10.0,\"StdDev\":0.01}}\n";
    out << "{\"DnaMeasurement\":{\"Type\":\"D\",\"First\":\"STN1\","
           "\"Second\":\"STN2\",\"EpochOfObservation\":\"16.07.2021\","
           "\"Value\":44.0,\"StdDev\":0.01,\"Total\":1,"
           "\"Directions\":[{\"Target\":\"STN3\",\"Value\":45.0,"
           "\"StdDev\":0.01}]}}\n";
    out << "{\"DnaMeasurement\":{\"Type\":\"X\",\"First\":\"STN1\","
           "\"Second\":\"STN2\",\"Epoch\":\"01.01.2020\","
           "\"EpochOfObservation\":\"17.07.2021\",\"Total\":1,"
           "\"GPSBaseline\":[{\"First\":\"STN1\",\"Second\":\"STN2\","
           "\"X\":1.0,\"Y\":2.0,\"Z\":3.0,\"SigmaXX\":0.0001,"
           "\"SigmaXY\":0.00001,\"SigmaXZ\":0.00002,"
           "\"SigmaYY\":0.0002,\"SigmaYZ\":0.00003,"
           "\"SigmaZZ\":0.0003}]}}\n";
    out << "{\"DnaMeasurement\":{\"Type\":\"Y\",\"First\":\"STN1\","
           "\"Epoch\":\"01.01.2020\","
           "\"EpochOfObservation\":\"18.07.2021\",\"Total\":1,"
           "\"Coords\":\"XYZ\",\"Clusterpoint\":[{\"First\":\"STN1\","
           "\"X\":4.0,\"Y\":5.0,\"Z\":6.0,\"SigmaXX\":0.0004,"
           "\"SigmaXY\":0.00004,\"SigmaXZ\":0.00005,"
           "\"SigmaYY\":0.0005,\"SigmaYZ\":0.00006,"
           "\"SigmaZZ\":0.0006}]}}\n";
  }

  vdnaStnPtr stations;
  vdnaMsrPtr measurements;
  JsonlParseContext ctx = DefaultContext();

  ParseJsonlFile(tmp.path.string(), &stations, &measurements, ctx);

  REQUIRE(measurements.size() == 4);
  REQUIRE(measurements[0]->GetObservationEpoch() == "15.07.2021");
  REQUIRE(measurements[1]->GetObservationEpoch() == "16.07.2021");
  REQUIRE(measurements[2]->GetObservationEpoch() == "17.07.2021");
  REQUIRE(measurements[3]->GetObservationEpoch() == "18.07.2021");

  auto* directions = measurements[1]->GetDirections_ptr();
  REQUIRE(directions != nullptr);
  REQUIRE(directions->at(0).GetObservationEpoch() == "16.07.2021");

  auto* baselines = measurements[2]->GetBaselines_ptr();
  REQUIRE(baselines != nullptr);
  REQUIRE(baselines->at(0).GetObservationEpoch() == "17.07.2021");

  auto* points = measurements[3]->GetPoints_ptr();
  REQUIRE(points != nullptr);
  REQUIRE(points->at(0).GetObservationEpoch() == "18.07.2021");
}

TEST_CASE("JSONL parser preserves grouped adjustment direction sets") {
  TempJsonlFile tmp("dynadjust-jsonl-direction-set-regression.jsonl");
  {
    std::ofstream out(tmp.path);
    out << "{\"DnaAdjustmentReport\":{\"type\":\"Adjustment\","
           "\"report\":\"adj\",\"software\":\"DynAdjust 1.4.0\","
           "\"referenceframe\":\"GDA2020\",\"epoch\":\"01.01.2020\"}}\n";
    out << "{\"DnaMeasurement\":{\"Type\":\"D\",\"First\":\"STN1\","
           "\"Second\":\"STN2\",\"Value\":44.0,\"StdDev\":0.01,"
           "\"Total\":2,\"Directions\":[{\"Target\":\"STN3\","
           "\"Value\":45.0,\"StdDev\":0.01,\"Adjusted\":45.1,"
           "\"Correction\":0.1},{\"Target\":\"STN4\","
           "\"Value\":46.0,\"StdDev\":0.01,\"Ignore\":true,"
           "\"Adjusted\":46.2,\"Correction\":0.2}]}}\n";
  }

  vdnaStnPtr stations;
  vdnaMsrPtr measurements;
  JsonlParseContext ctx = DefaultContext();

  ParseJsonlFile(tmp.path.string(), &stations, &measurements, ctx);

  REQUIRE(ctx.stn_count == 0);
  REQUIRE(ctx.msr_count == 2);
  REQUIRE(stations.empty());
  REQUIRE(measurements.size() == 1);

  auto* directions = measurements[0]->GetDirections_ptr();
  REQUIRE(directions != nullptr);
  REQUIRE(directions->size() == 2);
  REQUIRE(directions->at(0).GetTarget() == "STN3");
  REQUIRE(directions->at(1).GetTarget() == "STN4");
  REQUIRE(directions->at(1).GetIgnore() == true);
}

TEST_CASE("JSONL parser treats boolean Ignore false as active") {
  TempJsonlFile tmp("dynadjust-jsonl-ignore-false-regression.jsonl");
  {
    std::ofstream out(tmp.path);
    out << "{\"DnaAdjustmentReport\":{\"type\":\"Adjustment\","
           "\"report\":\"adj\",\"software\":\"DynAdjust 1.4.0\","
           "\"referenceframe\":\"GDA2020\",\"epoch\":\"01.01.2020\"}}\n";
    out << "{\"DnaMeasurement\":{\"Type\":\"L\",\"First\":\"STN1\","
           "\"Second\":\"STN2\",\"Value\":10.0,\"StdDev\":0.01,"
           "\"Ignore\":false}}\n";
    out << "{\"DnaMeasurement\":{\"Type\":\"D\",\"First\":\"STN1\","
           "\"Second\":\"STN2\",\"Value\":44.0,\"StdDev\":0.01,"
           "\"Ignore\":false,\"Total\":2,\"Directions\":[{"
           "\"Target\":\"STN3\",\"Value\":45.0,\"StdDev\":0.01,"
           "\"Ignore\":false},{\"Target\":\"STN4\",\"Value\":46.0,"
           "\"StdDev\":0.01,\"Ignore\":true}]}}\n";
  }

  vdnaStnPtr stations;
  vdnaMsrPtr measurements;
  JsonlParseContext ctx = DefaultContext();

  ParseJsonlFile(tmp.path.string(), &stations, &measurements, ctx);

  REQUIRE(ctx.stn_count == 0);
  REQUIRE(ctx.msr_count == 3);
  REQUIRE(stations.empty());
  REQUIRE(measurements.size() == 2);
  REQUIRE(measurements[0]->GetIgnore() == false);
  REQUIRE(measurements[1]->GetIgnore() == false);

  auto* directions = measurements[1]->GetDirections_ptr();
  REQUIRE(directions != nullptr);
  REQUIRE(directions->size() == 2);
  REQUIRE(directions->at(0).GetIgnore() == false);
  REQUIRE(directions->at(1).GetIgnore() == true);
}
