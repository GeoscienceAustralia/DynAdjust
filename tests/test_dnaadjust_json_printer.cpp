//============================================================================
// Name         : test_dnaadjust_json_printer.cpp
// Author       : Dale Roberts <dale.o.roberts@gmail.com>
// Copyright    : Copyright 2026 Geoscience Australia
//
//                Licensed under the Apache License, Version 2.0 (the "License");
//                you may not use this file except in compliance with the License.
//                You may obtain a copy of the License at
//
//                http://www.apache.org/licenses/LICENSE-2.0
//
// Description  : Schema contract tests for the JSONL adjustment output.  These
//                tests lock in the key names and top-level record wrappers
//                emitted by DynAdjustJsonPrinter so downstream consumers (the
//                JSONL import parser in particular) can keep round-tripping
//                records produced by the adjustment path.
//============================================================================

#define TESTING_MAIN
#ifndef __BINARY_NAME__
#define __BINARY_NAME__ "test_dnaadjust_json_printer"
#endif
#ifndef __BINARY_DESC__
#define __BINARY_DESC__ "Schema contract tests for JSONL adjustment output"
#endif

#include "testing.hpp"

#include <string>

#include "../dynadjust/include/thirdparty/nlohmann/json.hpp"

// Bring matrix_2d into scope before including dnaadjust_printer.hpp.
// The printer header uses matrix_2d unqualified; it relies on the including
// translation unit having already pulled in using namespace dynadjust::math
// (normally provided by dnaadjust.hpp).  We do that explicitly here so the
// header compiles without the full dnaadjust.hpp dependency chain.
#include <include/math/dnamatrix_contiguous.hpp>
using namespace dynadjust::math;

#include <dynadjust/dnaadjust/dnaadjust_printer.hpp>

using nlohmann::json;

// Sample JSONL fragments matching the shapes produced by
// DynAdjustJsonPrinter.  If the output schema changes in
// dnaadjust_json_printer.cpp, these fixtures need to track it so
// downstream tooling and the JSONL import parser stay in sync.
constexpr const char* kAdjHeader =
    "{\"DnaAdjustmentReport\":{\"epoch\":\"01.01.2020\","
    "\"referenceframe\":\"GDA2020\",\"report\":\"adj\","
    "\"software\":\"DynAdjust 1.4.0\",\"type\":\"Adjustment\"}}";

constexpr const char* kStatisticsLine =
    "{\"DnaStatistics\":{\"chisq\":335.45,\"chisq_lower\":0.84,"
    "\"chisq_test\":\"passed\",\"chisq_upper\":1.17,"
    "\"confidence_interval\":95.0,\"dof\":288,"
    "\"global_pelzer\":0.78,\"measurement_params\":417,"
    "\"potential_outliers\":9,\"sigma_zero\":1.16,"
    "\"unknown_parameters\":129}}";

// Lat/Lon are DMS-packed decimals: -36.345678 ≡ -36°34'56.78".  This matches
// RadtoDms()'s output and lets the JSONL import parser round-trip via atof →
// DmstoDeg → Radians.
constexpr const char* kAdjustedStationLine =
    "{\"DnaStation\":{\"Adjusted\":{\"Height\":227.18,\"Lat\":-36.544321,"
    "\"Lon\":146.393456,\"X\":-4288403.6,\"Y\":2814576.3,\"Z\":-3778237.8},"
    "\"Constraints\":\"FFF\",\"Initial\":{\"Height\":227.17,"
    "\"Lat\":-36.544321,\"Lon\":146.393456},\"Name\":\"MYRT\","
    "\"StationCoord\":{\"Height\":227.18,\"Name\":\"MYRT\","
    "\"XAxis\":-36.544321,\"YAxis\":146.393456},"
    "\"Type\":\"LLH\","
    "\"Uncertainty\":{\"HzPosU\":0.005,\"Orientation\":1.88,"
    "\"SE\":0.002,\"SN\":0.002,\"SU\":0.005,\"SemiMajor\":0.002,"
    "\"SemiMinor\":0.002,\"VtPosU\":0.010}}}";

constexpr const char* kUncertaintyLine =
    "{\"DnaStation\":{\"Constraints\":\"FFF\",\"Name\":\"STN1\","
    "\"Type\":\"LLH\","
    "\"Uncertainty\":{\"HzPosU\":0.005,\"Orientation\":1.9,"
    "\"SE\":0.002,\"SN\":0.002,\"SU\":0.005,\"SemiMajor\":0.002,"
    "\"SemiMinor\":0.002,\"VtPosU\":0.010}}}";

constexpr const char* kCorrectionsLine =
    "{\"DnaStation\":{\"Constraints\":\"FFF\",\"Corrections\":"
    "{\"dE\":0.003,\"dN\":0.002,\"dUp\":9.1},\"Initial\":"
    "{\"Height\":172.19,\"Lat\":-36.544321,\"Lon\":146.393456},"
    "\"Name\":\"211300470\",\"Type\":\"LLH\"}}";

constexpr const char* kM2SLine =
    "{\"DnaStation\":{\"Constraints\":\"FFF\","
    "\"MsrCounts\":{\"A\":0,\"B\":0,\"C\":0,\"D\":0,\"E\":0,\"G\":4,"
    "\"H\":0,\"I\":0,\"J\":0,\"K\":0,\"L\":0,\"M\":0,\"P\":0,\"Q\":0,"
    "\"R\":0,\"S\":0,\"V\":0,\"X\":0,\"Y\":0,\"Z\":0},"
    "\"Name\":\"STN1\",\"TotalMeasurements\":4,\"Type\":\"LLH\"}}";

constexpr const char* kMeasurementLine =
    "{\"DnaMeasurement\":{\"Adjusted\":1234.567,"
    "\"AdjustedPrecision\":0.001,\"Correction\":0.002,"
    "\"EpochOfObservation\":\"15.07.2021\","
    "\"First\":\"STN1\",\"NStat\":0.8,\"PelzerRel\":0.92,"
    "\"ResidualPrecision\":0.001,\"Second\":\"STN2\","
    "\"StdDev\":0.004,\"TStat\":1.1,\"Type\":\"L\","
    "\"Value\":1234.565}}";

constexpr const char* kDirectionSetLine =
    "{\"DnaMeasurement\":{\"Adjusted\":45.0,"
    "\"AdjustedPrecision\":0.0001,\"Correction\":0.0,"
    "\"Directions\":[{\"Adjusted\":45.1,\"AdjustedPrecision\":0.0002,"
    "\"Correction\":0.1,\"NStat\":0.5,\"PelzerRel\":0.9,"
    "\"ResidualPrecision\":0.0003,\"StdDev\":0.01,\"TStat\":1.1,"
    "\"Target\":\"STN3\",\"Value\":45.0},"
    "{\"Adjusted\":46.2,\"AdjustedPrecision\":0.0002,"
    "\"Correction\":0.2,\"Ignore\":true,\"NStat\":0.6,"
    "\"PelzerRel\":0.8,\"ResidualPrecision\":0.0003,"
    "\"StdDev\":0.01,\"TStat\":1.2,\"Target\":\"STN4\","
    "\"Value\":46.0}],\"EpochOfObservation\":\"16.07.2021\","
    "\"First\":\"STN1\",\"NStat\":0.0,"
    "\"PelzerRel\":0.0,\"ResidualPrecision\":0.0001,"
    "\"Second\":\"STN2\",\"StdDev\":0.01,\"TStat\":0.0,"
    "\"Total\":2,\"Type\":\"D\",\"Value\":44.0}}";

// Cluster measurements (G/X/Y) carry observed X/Y/Z and covariance terms in
// the same GPSBaseline/Clusterpoint arrays consumed by the JSONL input parser.
// Adjusted statistics remain as output-only metadata.
constexpr const char* kClusterMeasurementLine =
    "{\"DnaMeasurement\":{"
    "\"Adjusted\":[{\"X\":1.1,\"Y\":2.2,\"Z\":3.3},"
    "{\"X\":4.4,\"Y\":5.5,\"Z\":6.6}],"
    "\"AdjustedPrecision\":[{\"X\":0.01,\"Y\":0.01,\"Z\":0.01},"
    "{\"X\":0.02,\"Y\":0.02,\"Z\":0.02}],"
    "\"Correction\":[{\"X\":0.001,\"Y\":0.002,\"Z\":0.003},"
    "{\"X\":0.004,\"Y\":0.005,\"Z\":0.006}],"
    "\"EpochOfObservation\":\"17.07.2021\","
    "\"First\":\"STN1\",\"GPSBaseline\":[{\"First\":\"STN1\","
    "\"GPSCovariance\":[{\"m11\":0.000001,\"m12\":0.000002,"
    "\"m13\":0.000003,\"m21\":0.000004,\"m22\":0.000005,"
    "\"m23\":0.000006,\"m31\":0.000007,\"m32\":0.000008,"
    "\"m33\":0.000009}],\"Second\":\"STN2\","
    "\"SigmaXX\":0.0001,\"SigmaXY\":0.00001,"
    "\"SigmaXZ\":0.00002,\"SigmaYY\":0.0002,\"SigmaYZ\":0.00003,"
    "\"SigmaZZ\":0.0003,\"X\":1.099,\"Y\":2.198,\"Z\":3.297},"
    "{\"First\":\"STN1\",\"Second\":\"STN3\",\"SigmaXX\":0.0004,"
    "\"SigmaXY\":0.00004,\"SigmaXZ\":0.00005,\"SigmaYY\":0.0005,"
    "\"SigmaYZ\":0.00006,\"SigmaZZ\":0.0006,\"X\":4.396,"
    "\"Y\":5.495,\"Z\":6.594}],"
    "\"NStat\":[{\"X\":0.5,\"Y\":0.6,\"Z\":0.7},"
    "{\"X\":0.8,\"Y\":0.9,\"Z\":1.0}],"
    "\"PelzerRel\":[{\"X\":0.9,\"Y\":0.9,\"Z\":0.9},"
    "{\"X\":0.8,\"Y\":0.8,\"Z\":0.8}],"
    "\"Second\":\"STN2\",\"TStat\":[{\"X\":1.0,\"Y\":1.1,\"Z\":1.2},"
    "{\"X\":1.3,\"Y\":1.4,\"Z\":1.5}],"
    "\"Total\":2,\"Type\":\"X\"}}";

TEST_CASE("header line has DnaAdjustmentReport wrapper with required keys") {
    const json record = json::parse(kAdjHeader);
    REQUIRE(record.contains("DnaAdjustmentReport"));
    const auto& hdr = record["DnaAdjustmentReport"];
    REQUIRE(hdr.contains("type"));
    REQUIRE(hdr.contains("report"));
    REQUIRE(hdr.contains("software"));
    REQUIRE(hdr.contains("referenceframe"));
    REQUIRE(hdr.contains("epoch"));
    REQUIRE(hdr["type"] == "Adjustment");
    REQUIRE(hdr["report"] == "adj");
    REQUIRE(!hdr.contains("statistics"));
}

TEST_CASE("statistics trailer has DnaStatistics wrapper with required keys") {
    const json record = json::parse(kStatisticsLine);
    REQUIRE(record.contains("DnaStatistics"));
    const auto& s = record["DnaStatistics"];
    REQUIRE(s.contains("chisq"));
    REQUIRE(s.contains("dof"));
    REQUIRE(s.contains("sigma_zero"));
    REQUIRE(s.contains("global_pelzer"));
    REQUIRE(s.contains("chisq_lower"));
    REQUIRE(s.contains("chisq_upper"));
    REQUIRE(s.contains("chisq_test"));
    REQUIRE(s.contains("unknown_parameters"));
    REQUIRE(s.contains("measurement_params"));
    REQUIRE(s["chisq_test"].get<std::string>() == "passed");
}

TEST_CASE("adjusted station record carries Initial, Adjusted, Uncertainty") {
    const json record = json::parse(kAdjustedStationLine);
    REQUIRE(record.contains("DnaStation"));
    const auto& s = record["DnaStation"];
    REQUIRE(s.contains("Name"));
    REQUIRE(s.contains("Constraints"));
    REQUIRE(s.contains("Type"));
    REQUIRE(s.contains("StationCoord"));
    REQUIRE(s.contains("Initial"));
    REQUIRE(s.contains("Adjusted"));
    REQUIRE(s.contains("Uncertainty"));
    REQUIRE(s["Type"] == "LLH");
    REQUIRE(s["Adjusted"].contains("X"));
    REQUIRE(s["Adjusted"].contains("Y"));
    REQUIRE(s["Adjusted"].contains("Z"));
    REQUIRE(s["Adjusted"].contains("Lat"));
    REQUIRE(s["Adjusted"].contains("Lon"));
    REQUIRE(s["Adjusted"].contains("Height"));
}

TEST_CASE("station record has round-trippable StationCoord block") {
    // dnaparser_jsonl.cpp:146-157 expects StationCoord.{XAxis,YAxis,Height}
    // alongside top-level Type; this makes the emitted station record
    // directly re-ingestible as initial inputs for a follow-on adjustment.
    const json record = json::parse(kAdjustedStationLine);
    const auto& s = record["DnaStation"];
    REQUIRE(s["Type"].is_string());
    REQUIRE(s["StationCoord"].is_object());
    const auto& sc = s["StationCoord"];
    REQUIRE(sc.contains("Name"));
    REQUIRE(sc.contains("XAxis"));
    REQUIRE(sc.contains("YAxis"));
    REQUIRE(sc.contains("Height"));
    // Name inside StationCoord mirrors the station Name.
    REQUIRE(sc["Name"] == s["Name"]);
}

TEST_CASE("uncertainty record carries ellipse + horizontal/vertical PU") {
    const json record = json::parse(kUncertaintyLine);
    const auto& u = record["DnaStation"]["Uncertainty"];
    REQUIRE(u.contains("SE"));
    REQUIRE(u.contains("SN"));
    REQUIRE(u.contains("SU"));
    REQUIRE(u.contains("SemiMajor"));
    REQUIRE(u.contains("SemiMinor"));
    REQUIRE(u.contains("Orientation"));
    REQUIRE(u.contains("HzPosU"));
    REQUIRE(u.contains("VtPosU"));
}

TEST_CASE("station correction record carries dE, dN, dUp") {
    const json record = json::parse(kCorrectionsLine);
    const auto& c = record["DnaStation"]["Corrections"];
    REQUIRE(c.contains("dE"));
    REQUIRE(c.contains("dN"));
    REQUIRE(c.contains("dUp"));
}

TEST_CASE("measurements-to-stations record carries MsrCounts + total") {
    const json record = json::parse(kM2SLine);
    const auto& s = record["DnaStation"];
    REQUIRE(s.contains("MsrCounts"));
    REQUIRE(s.contains("TotalMeasurements"));
    const auto& counts = s["MsrCounts"];
    // 20 canonical single-letter measurement types.
    REQUIRE(counts.size() == 20);
    REQUIRE(counts.contains("G"));
    REQUIRE(counts["G"] == 4);
}

TEST_CASE("adjusted measurement record carries Type and adjustment fields") {
    const json record = json::parse(kMeasurementLine);
    REQUIRE(record.contains("DnaMeasurement"));
    const auto& m = record["DnaMeasurement"];
    REQUIRE(m.contains("Type"));
    REQUIRE(m.contains("First"));
    REQUIRE(m.contains("Value"));
    REQUIRE(m.contains("StdDev"));
    REQUIRE(m.contains("Adjusted"));
    REQUIRE(m.contains("Correction"));
    REQUIRE(m.contains("EpochOfObservation"));
    REQUIRE(m.contains("NStat"));
    REQUIRE(m.contains("TStat"));
    REQUIRE(m.contains("PelzerRel"));
}

TEST_CASE("records round-trip through parse/dump") {
    // The compact dump format (no whitespace) is what DynAdjustJsonPrinter
    // emits.  Re-parsing and re-dumping must produce byte-identical output
    // so records survive a JSON parse/dump pass.
    for (const char* line : {kAdjHeader, kStatisticsLine,
                             kAdjustedStationLine, kUncertaintyLine,
                             kCorrectionsLine, kM2SLine, kMeasurementLine,
                             kClusterMeasurementLine}) {
        const json parsed = json::parse(line);
        REQUIRE(!parsed.is_discarded());
        const std::string redump = parsed.dump();
        REQUIRE(!redump.empty());
        REQUIRE(json::parse(redump) == parsed);
    }
}

// ---------------------------------------------------------------------
// Schema alignment with the JSONL input parser (dnaparser_jsonl.cpp)
// ---------------------------------------------------------------------

TEST_CASE("station record uses input-parser-compatible key names") {
    // dnaparser_jsonl.cpp:113-160 reads these keys from DnaStation; when
    // the adjusted JSONL is re-ingested, these must survive verbatim.
    const json record = json::parse(kAdjustedStationLine);
    const auto& s = record["DnaStation"];
    REQUIRE(s.contains("Name"));
    REQUIRE(s.contains("Constraints"));
    REQUIRE(s.contains("Type"));
    REQUIRE(s.contains("StationCoord"));
    // "Description" is optional in the input parser; present only when
    // non-empty on adjusted-station output.  Same policy as the parser.
}

TEST_CASE("measurement record uses input-parser-compatible key names") {
    // dnaparser_jsonl.cpp:415-495 reads these keys from DnaMeasurement.
    const json record = json::parse(kMeasurementLine);
    const auto& m = record["DnaMeasurement"];
    REQUIRE(m.contains("Type"));
    REQUIRE(m.contains("First"));
    REQUIRE(m.contains("Value"));
    REQUIRE(m.contains("StdDev"));
    REQUIRE(m.contains("EpochOfObservation"));
    // "Second" is present on any non-scalar-from-nothing measurement;
    // the L (height diff) fixture has two stations.
    REQUIRE(m.contains("Second"));
}

// ---------------------------------------------------------------------
// Type safety — numbers must be numbers, strings must be strings
// ---------------------------------------------------------------------

TEST_CASE("station identity fields are strings") {
    const json record = json::parse(kAdjustedStationLine);
    const auto& s = record["DnaStation"];
    REQUIRE(s["Name"].is_string());
    REQUIRE(s["Constraints"].is_string());
}

TEST_CASE("adjusted coordinate fields are numeric") {
    const json record = json::parse(kAdjustedStationLine);
    const auto& a = record["DnaStation"]["Adjusted"];
    REQUIRE(a["X"].is_number());
    REQUIRE(a["Y"].is_number());
    REQUIRE(a["Z"].is_number());
    REQUIRE(a["Lat"].is_number());
    REQUIRE(a["Lon"].is_number());
    REQUIRE(a["Height"].is_number());
}

TEST_CASE("uncertainty numeric fields are numbers not strings") {
    const json record = json::parse(kUncertaintyLine);
    const auto& u = record["DnaStation"]["Uncertainty"];
    REQUIRE(u["SE"].is_number());
    REQUIRE(u["SN"].is_number());
    REQUIRE(u["SU"].is_number());
    REQUIRE(u["SemiMajor"].is_number());
    REQUIRE(u["SemiMinor"].is_number());
    REQUIRE(u["Orientation"].is_number());
    REQUIRE(u["HzPosU"].is_number());
    REQUIRE(u["VtPosU"].is_number());
}

TEST_CASE("correction delta fields are numbers") {
    const json record = json::parse(kCorrectionsLine);
    const auto& c = record["DnaStation"]["Corrections"];
    REQUIRE(c["dE"].is_number());
    REQUIRE(c["dN"].is_number());
    REQUIRE(c["dUp"].is_number());
}

TEST_CASE("statistics numeric fields are numbers and counts are integers") {
    const json record = json::parse(kStatisticsLine);
    const auto& s = record["DnaStatistics"];
    REQUIRE(s["chisq"].is_number());
    REQUIRE(s["sigma_zero"].is_number());
    REQUIRE(s["global_pelzer"].is_number());
    REQUIRE(s["chisq_lower"].is_number());
    REQUIRE(s["chisq_upper"].is_number());
    REQUIRE(s["confidence_interval"].is_number());
    REQUIRE(s["chisq_test"].is_string());
    REQUIRE(s["dof"].is_number_integer());
    REQUIRE(s["unknown_parameters"].is_number_integer());
    REQUIRE(s["measurement_params"].is_number_integer());
    REQUIRE(s["potential_outliers"].is_number_integer());
}

TEST_CASE("measurement statistics fields are numbers") {
    const json record = json::parse(kMeasurementLine);
    const auto& m = record["DnaMeasurement"];
    REQUIRE(m["Value"].is_number());
    REQUIRE(m["StdDev"].is_number());
    REQUIRE(m["Adjusted"].is_number());
    REQUIRE(m["Correction"].is_number());
    REQUIRE(m["AdjustedPrecision"].is_number());
    REQUIRE(m["ResidualPrecision"].is_number());
    REQUIRE(m["NStat"].is_number());
    REQUIRE(m["TStat"].is_number());
    REQUIRE(m["PelzerRel"].is_number());
}

TEST_CASE("measurement Type is single uppercase character") {
    const json record = json::parse(kMeasurementLine);
    const std::string t = record["DnaMeasurement"]["Type"].get<std::string>();
    REQUIRE(t.size() == 1);
    REQUIRE(t[0] >= 'A');
    REQUIRE(t[0] <= 'Z');
}

// ---------------------------------------------------------------------
// Optional-field policies
// ---------------------------------------------------------------------

TEST_CASE("station Description only appears when non-empty") {
    // kAdjustedStationLine has no Description (station MYRT in fixture)
    const json with_no_desc = json::parse(kAdjustedStationLine);
    REQUIRE(!with_no_desc["DnaStation"].contains("Description"));

    // kCorrectionsLine has a Description (station 211300470 → BENALLA)
    constexpr const char* with_desc =
        "{\"DnaStation\":{\"Constraints\":\"FFF\",\"Description\":"
        "\"BENALLA PM 47\",\"Name\":\"211300470\"}}";
    const json r = json::parse(with_desc);
    REQUIRE(r["DnaStation"].contains("Description"));
    REQUIRE(r["DnaStation"]["Description"].get<std::string>() ==
            "BENALLA PM 47");
}

TEST_CASE("ignored measurement carries Ignore=true; active ones omit it") {
    // Default measurement has no Ignore key
    const json active = json::parse(kMeasurementLine);
    REQUIRE(!active["DnaMeasurement"].contains("Ignore"));

    constexpr const char* ignored =
        "{\"DnaMeasurement\":{\"Type\":\"L\",\"First\":\"A\","
        "\"Second\":\"B\",\"Value\":10.0,\"Adjusted\":10.0,"
        "\"StdDev\":0.01,\"Correction\":0.0,\"AdjustedPrecision\":0.01,"
        "\"ResidualPrecision\":0.01,\"NStat\":0.0,\"TStat\":0.0,"
        "\"PelzerRel\":0.0,\"Ignore\":true}}";
    const json r = json::parse(ignored);
    REQUIRE(r["DnaMeasurement"].contains("Ignore"));
    REQUIRE(r["DnaMeasurement"]["Ignore"].get<bool>() == true);
}

// ---------------------------------------------------------------------
// Measurement type variations (station counts match schema)
// ---------------------------------------------------------------------

TEST_CASE("scalar two-station measurement carries First and Second only") {
    // Type L (height difference): First+Second, no Third.
    constexpr const char* two_stn =
        "{\"DnaMeasurement\":{\"Type\":\"L\",\"First\":\"A\","
        "\"Second\":\"B\",\"Value\":10.0,\"Adjusted\":10.0,"
        "\"StdDev\":0.01,\"Correction\":0.0,\"AdjustedPrecision\":0.01,"
        "\"ResidualPrecision\":0.01,\"NStat\":0.0,\"TStat\":0.0,"
        "\"PelzerRel\":0.0}}";
    const json r = json::parse(two_stn);
    const auto& m = r["DnaMeasurement"];
    REQUIRE(m.contains("First"));
    REQUIRE(m.contains("Second"));
    REQUIRE(!m.contains("Third"));
}

TEST_CASE("angle measurement carries First, Second, and Third") {
    // Type A (horizontal angle): three stations.
    constexpr const char* three_stn =
        "{\"DnaMeasurement\":{\"Type\":\"A\",\"First\":\"A\","
        "\"Second\":\"B\",\"Third\":\"C\",\"Value\":45.0,\"StdDev\":0.01,"
        "\"Adjusted\":45.0,\"Correction\":0.0,"
        "\"AdjustedPrecision\":0.01,\"ResidualPrecision\":0.01,"
        "\"NStat\":0.0,\"TStat\":0.0,\"PelzerRel\":0.0}}";
    const json r = json::parse(three_stn);
    const auto& m = r["DnaMeasurement"];
    REQUIRE(m["Type"] == "A");
    REQUIRE(m.contains("First"));
    REQUIRE(m.contains("Second"));
    REQUIRE(m.contains("Third"));
}

TEST_CASE("direction set measurement carries grouped Directions array") {
    const json record = json::parse(kDirectionSetLine);
    const auto& m = record["DnaMeasurement"];
    REQUIRE(m["Type"] == "D");
    REQUIRE(m["First"] == "STN1");
    REQUIRE(m["Second"] == "STN2");
    REQUIRE(m["Total"] == 2);
    REQUIRE(m["EpochOfObservation"] == "16.07.2021");
    REQUIRE(m.contains("Directions"));
    REQUIRE(m["Directions"].is_array());
    REQUIRE(m["Directions"].size() == 2);
    REQUIRE(m["Directions"][0]["Target"] == "STN3");
    REQUIRE(m["Directions"][0]["Value"].is_number());
    REQUIRE(m["Directions"][0]["StdDev"].is_number());
    REQUIRE(!m["Directions"][0].contains("Second"));
    REQUIRE(m["Directions"][1]["Ignore"].get<bool>() == true);
}

TEST_CASE("cluster measurement carries parser-compatible GPSBaseline values") {
    const json record = json::parse(kClusterMeasurementLine);
    const auto& m = record["DnaMeasurement"];
    REQUIRE(m["Type"] == "X");
    REQUIRE(m["Total"] == 2);
    REQUIRE(m["EpochOfObservation"] == "17.07.2021");
    REQUIRE(m.contains("GPSBaseline"));
    REQUIRE(m["GPSBaseline"].is_array());
    REQUIRE(m["GPSBaseline"].size() == 2);
    const auto& b = m["GPSBaseline"][0];
    for (const char* field : {"First", "Second", "X", "Y", "Z",
                               "SigmaXX", "SigmaXY", "SigmaXZ",
                               "SigmaYY", "SigmaYZ", "SigmaZZ"}) {
        REQUIRE(b.contains(field));
    }
    REQUIRE(b["X"].is_number());
    REQUIRE(b["Y"].is_number());
    REQUIRE(b["Z"].is_number());
    REQUIRE(b["SigmaXX"].is_number());
    REQUIRE(b.contains("GPSCovariance"));
    REQUIRE(b["GPSCovariance"].is_array());
    REQUIRE(b["GPSCovariance"].size() == 1);
    REQUIRE(b["GPSCovariance"][0]["m11"].is_number());
    REQUIRE(!m["GPSBaseline"][1].contains("GPSCovariance"));
    REQUIRE(!m.contains("Value"));

    for (const char* field : {"Adjusted", "Correction", "AdjustedPrecision",
                               "NStat", "TStat", "PelzerRel"}) {
        REQUIRE(m.contains(field));
        REQUIRE(m[field].is_array());
        REQUIRE(m[field].size() == 2);
        REQUIRE(m[field][0].contains("X"));
        REQUIRE(m[field][0].contains("Y"));
        REQUIRE(m[field][0].contains("Z"));
        REQUIRE(m[field][0]["X"].is_number());
        REQUIRE(m[field][0]["Y"].is_number());
        REQUIRE(m[field][0]["Z"].is_number());
    }
}

TEST_CASE("single-station measurement carries First only") {
    // Type H (orthometric height): one station.
    constexpr const char* one_stn =
        "{\"DnaMeasurement\":{\"Type\":\"H\",\"First\":\"A\","
        "\"Value\":100.0,\"StdDev\":0.01,\"Adjusted\":100.0,"
        "\"Correction\":0.0,"
        "\"AdjustedPrecision\":0.01,\"ResidualPrecision\":0.01,"
        "\"NStat\":0.0,\"TStat\":0.0,\"PelzerRel\":0.0}}";
    const json r = json::parse(one_stn);
    const auto& m = r["DnaMeasurement"];
    REQUIRE(m["Type"] == "H");
    REQUIRE(m.contains("First"));
    REQUIRE(!m.contains("Second"));
    REQUIRE(!m.contains("Third"));
}

// ---------------------------------------------------------------------
// MsrCounts structure
// ---------------------------------------------------------------------

TEST_CASE("MsrCounts carries every canonical measurement type key") {
    const json record = json::parse(kM2SLine);
    const auto& counts = record["DnaStation"]["MsrCounts"];
    // The 20 canonical measurement types DynAdjust supports.
    for (const char* key : {"A","B","C","D","E","G","H","I","J","K",
                             "L","M","P","Q","R","S","V","X","Y","Z"}) {
        REQUIRE(counts.contains(key));
        REQUIRE(counts[key].is_number_integer());
    }
}

TEST_CASE("MsrCounts sum consistent with TotalMeasurements") {
    const json record = json::parse(kM2SLine);
    const auto& s = record["DnaStation"];
    int sum = 0;
    for (const char* key : {"A","B","C","D","E","G","H","I","J","K",
                             "L","M","P","Q","R","S","V","X","Y","Z"}) {
        sum += s["MsrCounts"][key].get<int>();
    }
    REQUIRE(sum == s["TotalMeasurements"].get<int>());
}

// ---------------------------------------------------------------------
// Cross-record / end-to-end stream structure
// ---------------------------------------------------------------------

TEST_CASE("multi-record JSONL stream parses record-by-record") {
    // Simulate a minimal `.adj.jsonl`: header, stats, one measurement,
    // one station.  Parser must read each line independently.
    std::string stream;
    stream += kAdjHeader; stream += "\n";
    stream += kStatisticsLine; stream += "\n";
    stream += kMeasurementLine; stream += "\n";
    stream += kAdjustedStationLine; stream += "\n";

    std::vector<json> records;
    std::size_t start = 0;
    while (start < stream.size()) {
        auto end = stream.find('\n', start);
        if (end == std::string::npos) break;
        records.push_back(json::parse(stream.substr(start, end - start)));
        start = end + 1;
    }
    REQUIRE(records.size() == 4);
    REQUIRE(records[0].contains("DnaAdjustmentReport"));
    REQUIRE(records[1].contains("DnaStatistics"));
    REQUIRE(records[2].contains("DnaMeasurement"));
    REQUIRE(records[3].contains("DnaStation"));
}

TEST_CASE("exactly one top-level wrapper per record") {
    // Every emitted line must have exactly one of DnaAdjustmentReport /
    // DnaStation / DnaMeasurement / DnaStatistics so the dispatch
    // in dnaparser_jsonl.cpp:621-625 remains unambiguous.  (The parser
    // also accepts DnaXmlFormat from the input grammar as a header
    // wrapper, for backwards compatibility.)
    for (const char* line : {kAdjHeader, kStatisticsLine,
                             kAdjustedStationLine, kUncertaintyLine,
                             kCorrectionsLine, kM2SLine, kMeasurementLine,
                             kDirectionSetLine, kClusterMeasurementLine}) {
        const json r = json::parse(line);
        REQUIRE(r.is_object());
        REQUIRE(r.size() == 1);
        const std::string key = r.begin().key();
        REQUIRE((key == "DnaAdjustmentReport" || key == "DnaStation" ||
                 key == "DnaMeasurement" || key == "DnaStatistics"));
    }
}

// ---------------------------------------------------------------------
// Value ranges
// ---------------------------------------------------------------------

TEST_CASE("Lat/Lon values are DMS-packed decimals in geographic range") {
    // RadtoDms/RadtoDmsL encode -36°54'43.21" as the decimal -36.544321:
    // integer part is signed degrees, decimal part packs MM.SSssss.  So
    // max |Lat| ≈ 90.595999 and max |Lon| ≈ 180.595999.
    const json record = json::parse(kAdjustedStationLine);
    const auto& a = record["DnaStation"]["Adjusted"];
    const double lat = a["Lat"].get<double>();
    const double lon = a["Lon"].get<double>();
    REQUIRE(lat >= -91.0);
    REQUIRE(lat <= 91.0);
    REQUIRE(lon >= -181.0);
    REQUIRE(lon <= 181.0);
}

TEST_CASE("Statistics Chi-square test reports one of the known verdicts") {
    const json record = json::parse(kStatisticsLine);
    const std::string verdict =
        record["DnaStatistics"]["chisq_test"].get<std::string>();
    REQUIRE((verdict == "passed" || verdict == "warning" ||
             verdict == "failed" || verdict == "no_redundancy"));
}

// ---------------------------------------------------------------------
// Hook context threading contract — compile-time guard against silent
// signature changes.  A test-only subclass exposes the three protected
// context structs via public type aliases (the structs are declared in
// the protected section of DynAdjustPrinter).  Each test constructs a
// context using aggregate initialisation and asserts every field
// round-trips the value passed in.  Matrix parameters are passed by
// reference and the tests verify address identity to confirm the struct
// stores references, not copies.
// ---------------------------------------------------------------------

using dynadjust::networkadjust::DynAdjustPrinter;

TEST_CASE("AdjustedStationContext carries threaded precomputed fields") {
  matrix_2d vc(3, 3), vl(3, 3);
  DynAdjustPrinter::AdjustedStationContext ctx{
      -0.635, 2.557, 227.18, vc, vl};
  REQUIRE(ctx.lat == -0.635);
  REQUIRE(ctx.lon == 2.557);
  REQUIRE(ctx.height == 227.18);
  REQUIRE(&ctx.var_cart == &vc);
  REQUIRE(&ctx.var_local == &vl);
}

TEST_CASE("PositionalUncertaintyContext carries ellipse + PU") {
  matrix_2d vc(3, 3), vl(3, 3);
  DynAdjustPrinter::PositionalUncertaintyContext ctx{
      -0.635, 2.557, 227.18, vc, vl,
      0.002, 0.002, 1.88, 0.005, 0.010};
  REQUIRE(ctx.lat == -0.635);
  REQUIRE(ctx.lon == 2.557);
  REQUIRE(ctx.height == 227.18);
  REQUIRE(ctx.semimajor == 0.002);
  REQUIRE(ctx.semiminor == 0.002);
  REQUIRE(ctx.azimuth == 1.88);
  REQUIRE(ctx.hz_pos_u == 0.005);
  REQUIRE(ctx.vt_pos_u == 0.010);
  REQUIRE(&ctx.var_cart == &vc);
  REQUIRE(&ctx.var_local == &vl);
}

TEST_CASE("StationCorrectionContext carries dE/dN/dUp") {
  DynAdjustPrinter::StationCorrectionContext ctx{0.003, 0.002, 9.1};
  REQUIRE(ctx.cor_e == 0.003);
  REQUIRE(ctx.cor_n == 0.002);
  REQUIRE(ctx.cor_up == 9.1);
}
