//============================================================================
// Name         : dnaparser_jsonl.cpp
// Author       : Dale Roberts <dale.o.roberts@gmail.com>
// Copyright    : Copyright 2025 Geoscience Australia
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
// Description  : JSONL parser for DynaML-equivalent station and measurement
//                data.  One JSON object per line, using the same element
//                names and hierarchy as the DynaML XML schema.
//============================================================================

#include <fstream>
#include <exception>
#include <sstream>
#include <string>

#include <include/thirdparty/nlohmann/json.hpp>

#include <dynadjust/dnaimport/dnaparser_jsonl.hpp>
#include <dynadjust/dnaimport/dnaparser_validators.hpp>
#include <include/exception/dnaexception.hpp>
#include <include/parameters/dnaepsg.hpp>
#include <include/measurement_types/dnastntally.hpp>

using json = nlohmann::json;
using dynadjust::exception::XMLInteropException;
using dynadjust::measurements::CDnaCovariance;
using dynadjust::measurements::CDnaDirection;
using dynadjust::measurements::CDnaGpsBaseline;
using dynadjust::measurements::CDnaGpsPoint;
using dynadjust::measurements::CDnaStation;
using dynadjust::measurements::MsrTally;
using dynadjust::measurements::StnTally;
using dynadjust::measurements::dnaMsrPtr;
using dynadjust::measurements::dnaStnPtr;
using dynadjust::measurements::vdnaMsrPtr;
using dynadjust::measurements::vdnaStnPtr;

extern MsrTally g_parsemsr_tally;
extern StnTally g_parsestn_tally;
extern UINT32 g_fileOrder;

namespace dynadjust {
namespace import {
namespace {

// Helper to get a string value from a JSON object, returning empty string
// if the key is absent.  Handles both JSON strings ("1.23") and JSON
// numbers (1.23) so that callers don't need to care about quoting.
std::string jstr(const json& obj, const char* key) {
  auto it = obj.find(key);
  if (it == obj.end() || it->is_null())
    return "";
  if (it->is_string())
    return it->get<std::string>();
  if (it->is_number()) {
    // dump() produces a full-precision decimal without quotes.
    return it->dump();
  }
  return it->dump();  // fallback for booleans etc.
}

bool jignore(const json& obj, const char* key) {
  auto it = obj.find(key);
  if (it == obj.end() || it->is_null())
    return false;
  if (it->is_boolean())
    return it->get<bool>();
  if (it->is_string())
    return !it->get<std::string>().empty();
  return !it->dump().empty();
}

// Throw with line number context.
[[noreturn]] void throwLine(const std::string& msg, UINT32 line_no) {
  std::stringstream ss;
  ss << "JSONL line " << line_no << ": " << msg;
  throw XMLInteropException(ss.str(), 0);
}

// -------------------------------------------------------------------
// Header
// -------------------------------------------------------------------

void parseHeader(const json& obj, JsonlParseContext& ctx, UINT32 line_no) {
  const auto& hdr = obj.contains("DnaAdjustmentReport")
                        ? obj.at("DnaAdjustmentReport")
                        : obj.at("DnaXmlFormat");
  ctx.file_epsg = "";
  ctx.file_epoch = "";

  std::string rf = jstr(hdr, "referenceframe");
  std::string ep = jstr(hdr, "epoch");

  if (!rf.empty()) {
    ctx.file_specified_frame = true;
    try {
      ctx.file_epsg = epsg::epsgStringFromName<std::string>(rf);
    } catch (const std::exception& e) {
      throwLine(e.what(), line_no);
    }
  }

  if (!ep.empty()) {
    ctx.file_specified_epoch = true;
    ctx.file_epoch = ep;
  }
}

// -------------------------------------------------------------------
// Station
// -------------------------------------------------------------------

void parseStation(const json& obj, vdnaStnPtr* vStations,
                  JsonlParseContext& ctx, UINT32 line_no) {
  const auto& stn = obj.at("DnaStation");

  std::string name = jstr(stn, "Name");
  ValidateStationName(name);

  std::string frame, epoch;
  if (ctx.override_input_frame)
    frame = ctx.default_frame;
  else if (!ctx.file_epsg.empty())
    frame = epsg::datumFromEpsgString<std::string>(ctx.file_epsg);
  else
    frame = ctx.default_frame;

  if (ctx.override_input_frame)
    epoch = ctx.default_epoch;
  else if (!ctx.file_epoch.empty())
    epoch = ctx.file_epoch;
  else
    epoch = ctx.default_epoch;

  dnaStnPtr station(new CDnaStation(frame, epoch));
  station->SetfileOrder(g_fileOrder++);
  station->SetName(name);

  std::string constraints = DefaultConstraints(jstr(stn, "Constraints"));
  station->SetConstraints(constraints);
  g_parsestn_tally.addstation(constraints);

  station->SetCoordType(jstr(stn, "Type"));

  // StationCoord sub-object
  auto it_coord = stn.find("StationCoord");
  if (it_coord != stn.end()) {
    const auto& coord = *it_coord;
    // Name inside StationCoord is the same; we already set it.
    station->SetXAxis(jstr(coord, "XAxis"));
    station->SetYAxis(jstr(coord, "YAxis"));
    station->SetHeight(jstr(coord, "Height"));

    std::string hz = jstr(coord, "HemisphereZone");
    if (!hz.empty())
      station->SetHemisphereZone(hz);
  }

  station->SetDescription(SanitiseDescription(jstr(stn, "Description")));

  vStations->push_back(station);
  ctx.stn_count++;
}

// -------------------------------------------------------------------
// GPS Covariance  (shared between X-baselines and Y-points)
// -------------------------------------------------------------------

void parseGPSCovariance(const json& cov_obj, CDnaGpsBaseline& bsl) {
  CDnaCovariance cov;
  cov.SetType("X");
  cov.SetClusterID(bsl.GetClusterID());
  cov.SetM11(jstr(cov_obj, "m11"));
  cov.SetM12(jstr(cov_obj, "m12"));
  cov.SetM13(jstr(cov_obj, "m13"));
  cov.SetM21(jstr(cov_obj, "m21"));
  cov.SetM22(jstr(cov_obj, "m22"));
  cov.SetM23(jstr(cov_obj, "m23"));
  cov.SetM31(jstr(cov_obj, "m31"));
  cov.SetM32(jstr(cov_obj, "m32"));
  cov.SetM33(jstr(cov_obj, "m33"));
  bsl.AddGpsCovariance(&cov);
}

void parsePointCovariance(const json& cov_obj, CDnaGpsPoint& pt) {
  CDnaCovariance cov;
  cov.SetType("Y");
  cov.SetClusterID(pt.GetClusterID());
  cov.SetM11(jstr(cov_obj, "m11"));
  cov.SetM12(jstr(cov_obj, "m12"));
  cov.SetM13(jstr(cov_obj, "m13"));
  cov.SetM21(jstr(cov_obj, "m21"));
  cov.SetM22(jstr(cov_obj, "m22"));
  cov.SetM23(jstr(cov_obj, "m23"));
  cov.SetM31(jstr(cov_obj, "m31"));
  cov.SetM32(jstr(cov_obj, "m32"));
  cov.SetM33(jstr(cov_obj, "m33"));
  pt.AddPointCovariance(&cov);
}

// -------------------------------------------------------------------
// GPS Baselines  (G and X types)
// -------------------------------------------------------------------

void parseGPSBaselines(const json& msr_obj, dnaMsrPtr& msr,
                       JsonlParseContext& ctx, UINT32 line_no) {
  auto it = msr_obj.find("GPSBaseline");
  if (it == msr_obj.end())
    throwLine("GPS measurement requires \"GPSBaseline\" array.", line_no);

  const auto& baselines = *it;
  if (!baselines.is_array())
    throwLine("\"GPSBaseline\" must be an array.", line_no);

  char type_c = msr->GetTypeC();
  UINT32 total = static_cast<UINT32>(baselines.size());

  // Reserve space in the cluster
  msr->ReserveGpsBaselinesCount(total);

  for (UINT32 i = 0; i < total; ++i) {
    const auto& bsl_obj = baselines[i];

    CDnaGpsBaseline bsl;
    bsl.SetType(msr->GetType());

    // First and Second: each baseline can override the parent's values
    // (for X-type clusters, each baseline may connect different stations)
    std::string first = jstr(bsl_obj, "First");
    bsl.SetFirst(first.empty() ? msr->GetFirst() : first);

    std::string second = jstr(bsl_obj, "Second");
    bsl.SetTarget(second.empty() ? msr->GetTarget() : second);

    bsl.SetReferenceFrame(msr->GetReferenceFrame());
    bsl.SetEpoch(msr->GetEpoch());
    bsl.SetObservationEpoch(msr->GetObservationEpoch());
    bsl.SetClusterID(msr->GetClusterID());
    bsl.SetClusterDBID(msr->GetClusterDBID(), msr->GetClusterDBIDset());

    // For X-type clusters with >1 baselines, reserve covariance space
    if (type_c == 'X' && total > 1)
      bsl.ReserveGpsCovariancesCount(total - 1);

    // Set X, Y, Z components — each increments tally 1x
    bsl.SetX(jstr(bsl_obj, "X"));
    ctx.msr_count++;
    if (type_c == 'G')
      g_parsemsr_tally.G++;
    else
      g_parsemsr_tally.X++;

    bsl.SetY(jstr(bsl_obj, "Y"));
    ctx.msr_count++;
    if (type_c == 'G')
      g_parsemsr_tally.G++;
    else
      g_parsemsr_tally.X++;

    bsl.SetZ(jstr(bsl_obj, "Z"));
    ctx.msr_count++;
    if (type_c == 'G')
      g_parsemsr_tally.G++;
    else
      g_parsemsr_tally.X++;

    // Per-baseline MeasurementID
    std::string bsl_dbid = jstr(bsl_obj, "MeasurementID");
    if (!bsl_dbid.empty())
      bsl.SetMeasurementDBID(bsl_dbid);

    // Sigma values
    bsl.SetSigmaXX(jstr(bsl_obj, "SigmaXX"));
    bsl.SetSigmaXY(jstr(bsl_obj, "SigmaXY"));
    bsl.SetSigmaXZ(jstr(bsl_obj, "SigmaXZ"));
    bsl.SetSigmaYY(jstr(bsl_obj, "SigmaYY"));
    bsl.SetSigmaYZ(jstr(bsl_obj, "SigmaYZ"));
    bsl.SetSigmaZZ(jstr(bsl_obj, "SigmaZZ"));

    // Covariances (only for X-type clusters)
    auto it_cov = bsl_obj.find("GPSCovariance");
    if (it_cov != bsl_obj.end() && it_cov->is_array()) {
      for (const auto& cov_obj : *it_cov) {
        parseGPSCovariance(cov_obj, bsl);
      }
    }

    msr->AddGpsBaseline(&bsl);
  }
}

// -------------------------------------------------------------------
// GPS Points  (Y type)
// -------------------------------------------------------------------

void parseClusterpoints(const json& msr_obj, dnaMsrPtr& msr,
                        JsonlParseContext& ctx, UINT32 line_no) {
  auto it = msr_obj.find("Clusterpoint");
  if (it == msr_obj.end())
    throwLine("Y measurement requires \"Clusterpoint\" array.", line_no);

  const auto& points = *it;
  if (!points.is_array())
    throwLine("\"Clusterpoint\" must be an array.", line_no);

  UINT32 total = static_cast<UINT32>(points.size());

  // Reserve space in the cluster
  msr->ReserveGpsPointsCount(total);

  for (UINT32 i = 0; i < total; ++i) {
    const auto& pt_obj = points[i];

    CDnaGpsPoint pt;
    pt.SetType(msr->GetType());
    pt.SetCoordType(msr->GetCoordType());

    // Each clusterpoint can have its own First station
    std::string first = jstr(pt_obj, "First");
    pt.SetFirst(first.empty() ? msr->GetFirst() : first);

    pt.SetClusterID(msr->GetClusterID());
    pt.SetClusterDBID(msr->GetClusterDBID(), msr->GetClusterDBIDset());
    pt.SetReferenceFrame(msr->GetReferenceFrame());
    pt.SetEpoch(msr->GetEpoch());
    pt.SetObservationEpoch(msr->GetObservationEpoch());

    // Set X, Y, Z components — each increments tally 1x
    pt.SetX(jstr(pt_obj, "X"));
    ctx.msr_count++;
    g_parsemsr_tally.Y++;

    pt.SetY(jstr(pt_obj, "Y"));
    ctx.msr_count++;
    g_parsemsr_tally.Y++;

    pt.SetZ(jstr(pt_obj, "Z"));
    ctx.msr_count++;
    g_parsemsr_tally.Y++;

    // Per-point MeasurementID
    std::string pt_dbid = jstr(pt_obj, "MeasurementID");
    if (!pt_dbid.empty())
      pt.SetMeasurementDBID(pt_dbid);

    // Sigma values
    pt.SetSigmaXX(jstr(pt_obj, "SigmaXX"));
    pt.SetSigmaXY(jstr(pt_obj, "SigmaXY"));
    pt.SetSigmaXZ(jstr(pt_obj, "SigmaXZ"));
    pt.SetSigmaYY(jstr(pt_obj, "SigmaYY"));
    pt.SetSigmaYZ(jstr(pt_obj, "SigmaYZ"));
    pt.SetSigmaZZ(jstr(pt_obj, "SigmaZZ"));

    // Point covariances
    auto it_cov = pt_obj.find("PointCovariance");
    if (it_cov != pt_obj.end() && it_cov->is_array()) {
      for (const auto& cov_obj : *it_cov) {
        parsePointCovariance(cov_obj, pt);
      }
    }

    msr->AddGpsPoint(&pt);
  }
}

// -------------------------------------------------------------------
// Directions  (D type)
// -------------------------------------------------------------------

void parseDirections(const json& msr_obj, dnaMsrPtr& msr,
                     JsonlParseContext& ctx, UINT32 line_no) {
  auto it = msr_obj.find("Directions");
  if (it == msr_obj.end())
    throwLine("D measurement requires \"Directions\" array.", line_no);

  const auto& directions = *it;
  if (!directions.is_array())
    throwLine("\"Directions\" must be an array.", line_no);

  UINT32 total = msr->GetTotal();

  for (const auto& dir_obj : directions) {
    CDnaDirection dir;
    dir.SetClusterID(msr->GetClusterID());
    dir.SetClusterDBID(msr->GetClusterDBID(), msr->GetClusterDBIDset());
    dir.SetEpoch(msr->GetEpoch());
    dir.SetObservationEpoch(msr->GetObservationEpoch());

    // Ignore for individual directions
    bool dir_ignored = jignore(dir_obj, "Ignore");
    dir.SetIgnore(dir_ignored);

    dir.SetTarget(jstr(dir_obj, "Target"));
    dir.SetValue(jstr(dir_obj, "Value"));
    dir.SetStdDev(jstr(dir_obj, "StdDev"));

    // Direction Value always increments measurement count
    ctx.msr_count++;

    msr->AddDirection(&dir);

    // Tally counting: only count non-ignored directions
    // (matching pimpl behaviour in Directions_pimpl::post_Directions)
    if (msr->GetIgnore())
      continue;

    if (!dir_ignored) {
      g_parsemsr_tally.D++;
    } else {
      ValidateDirectionInSet(msr, true, total);
    }
  }
}

// -------------------------------------------------------------------
// Measurement dispatch
// -------------------------------------------------------------------

void parseMeasurement(const json& obj, vdnaMsrPtr* vMeasurements,
                      JsonlParseContext& ctx, UINT32 line_no) {
  const auto& msr_obj = obj.at("DnaMeasurement");

  std::string type_str = jstr(msr_obj, "Type");
  char type_c = ValidateMeasurementType(type_str);

  // Determine frame and epoch for this measurement
  std::string frame, epoch;
  std::string msr_rf = jstr(msr_obj, "ReferenceFrame");
  std::string msr_ep = jstr(msr_obj, "Epoch");
  std::string msr_obs_ep = jstr(msr_obj, "EpochOfObservation");

  bool rf_supplied = !msr_rf.empty();
  bool ep_supplied = !msr_ep.empty();
  bool obs_ep_supplied = !msr_obs_ep.empty();

  if (rf_supplied || ctx.override_input_frame)
    frame = ctx.override_input_frame ? ctx.default_frame : msr_rf;
  else if (!ctx.file_epsg.empty())
    frame = epsg::datumFromEpsgString<std::string>(ctx.file_epsg);
  else
    frame = ctx.default_frame;

  if (ep_supplied || ctx.override_input_frame)
    epoch = ctx.override_input_frame ? ctx.default_epoch : msr_ep;
  else if (!ctx.file_epoch.empty())
    epoch = ctx.file_epoch;
  else
    epoch = ctx.default_epoch;

  // Create measurement (increments tally for scalar types)
  dnaMsrPtr msr = CreateMeasurement(type_c, &ctx.cluster_id,
      ctx.msr_count, g_parsemsr_tally, frame, epoch);
  msr->SetType(type_str);

  // Ignore
  if (jignore(msr_obj, "Ignore")) {
    msr->SetIgnore(true);
    g_parsemsr_tally.ignored++;
  } else {
    msr->SetIgnore(false);
  }

  // Source, MeasurementID, ClusterID
  msr->SetSource(jstr(msr_obj, "Source"));

  std::string msr_dbid = jstr(msr_obj, "MeasurementID");
  if (!msr_dbid.empty())
    msr->SetMeasurementDBID(msr_dbid);

  std::string cluster_dbid = jstr(msr_obj, "ClusterID");
  if (!cluster_dbid.empty())
    msr->SetClusterDBID(cluster_dbid);

  // First station
  std::string first = jstr(msr_obj, "First");
  ValidateFirst(first);
  msr->SetFirst(first);

  // Second station
  std::string second = jstr(msr_obj, "Second");
  ValidateSecond(second, type_c, first);
  if (TypeRequiresSecond(type_c))
    msr->SetTarget(second);

  // Third station (only for angles)
  std::string third = jstr(msr_obj, "Third");
  ValidateThird(third, type_c, first);
  if (!third.empty())
    msr->SetTarget2(third);

  // Value and StdDev (scalar types)
  std::string value = jstr(msr_obj, "Value");
  ValidateValue(value, type_c, first);
  if (!value.empty())
    msr->SetValue(value);

  std::string stddev = jstr(msr_obj, "StdDev");
  ValidateStdDev(stddev, type_c, first);
  if (!stddev.empty())
    msr->SetStdDev(stddev);

  // Instrument and target heights (S, V, Z types)
  switch (type_c) {
  case 'S': case 'V': case 'Z': {
    std::string ih = jstr(msr_obj, "InstHeight");
    if (!ih.empty())
      msr->SetInstrumentHeight(ih);
    std::string th = jstr(msr_obj, "TargHeight");
    if (!th.empty())
      msr->SetTargetHeight(th);
    break;
  }
  }

  // Reference frame and epoch (for GPS types, already set via
  // CreateMeasurement but may need per-measurement override)
  if (rf_supplied && !ctx.override_input_frame)
    msr->SetReferenceFrame(msr_rf);

  if (ep_supplied && !ctx.override_input_frame)
    msr->SetEpoch(msr_ep);

  if (obs_ep_supplied)
    msr->SetObservationEpoch(msr_obs_ep);

  // Scales
  std::string vscale = jstr(msr_obj, "Vscale");
  if (!vscale.empty())
    msr->SetVscale(DefaultScale(vscale));

  std::string pscale = jstr(msr_obj, "Pscale");
  if (!pscale.empty())
    msr->SetPscale(DefaultScale(pscale));

  std::string lscale = jstr(msr_obj, "Lscale");
  if (!lscale.empty())
    msr->SetLscale(DefaultScale(lscale));

  std::string hscale = jstr(msr_obj, "Hscale");
  if (!hscale.empty())
    msr->SetHscale(DefaultScale(hscale));

  // Total (for cluster/set types D, X, Y; G auto-sets from baseline count)
  std::string total = jstr(msr_obj, "Total");
  if (type_c != 'G') {
    ValidateTotal(total, type_c);
    if (!total.empty())
      msr->SetTotal(total);
  }

  // Coords (for Y type)
  std::string coords = jstr(msr_obj, "Coords");
  if (!coords.empty())
    msr->SetCoordType(coords);

  // Dispatch to sub-parsers for complex types
  switch (type_c) {
  case 'G':
  case 'X':
    parseGPSBaselines(msr_obj, msr, ctx, line_no);
    // For G type, set Total from actual baseline count
    if (type_c == 'G')
      msr->SetTotal(static_cast<UINT32>(msr->GetBaselines_ptr()->size()));
    break;
  case 'Y':
    parseClusterpoints(msr_obj, msr, ctx, line_no);
    break;
  case 'D':
    parseDirections(msr_obj, msr, ctx, line_no);
    break;
  }

  // Ignored empty direction set — skip silently (matching pimpl)
  if (type_c == 'D' && msr->GetTotal() == 0 && msr->GetIgnore())
    return;

  // Finalise and add
  vMeasurements->push_back(msr);
  FinaliseMeasurement(msr, ctx.prefer_single_x_as_g);
}

}  // anonymous namespace

// -------------------------------------------------------------------
// Public entry point
// -------------------------------------------------------------------

void ParseJsonlFile(const std::string& filename,
                    vdnaStnPtr* vStations,
                    vdnaMsrPtr* vMeasurements,
                    JsonlParseContext& ctx) {
  std::ifstream ifs(filename);
  if (!ifs.is_open()) {
    std::stringstream ss;
    ss << "ParseJsonlFile(): Could not open " << filename;
    throw XMLInteropException(ss.str(), 0);
  }

  // Initialise tallies
  g_parsemsr_tally.initialise();
  g_parsestn_tally.initialise();

  ctx.stn_count = 0;
  ctx.msr_count = 0;
  ctx.file_specified_frame = false;
  ctx.file_specified_epoch = false;
  ctx.file_epsg = "";
  ctx.file_epoch = "";

  std::string line;
  UINT32 line_no = 0;
  while (std::getline(ifs, line)) {
    ++line_no;

    // Skip empty lines
    if (line.empty() || line.find_first_not_of(" \t\r\n") == std::string::npos)
      continue;

    json obj;
    try {
      obj = json::parse(line);
    } catch (const json::parse_error& e) {
      throwLine(std::string("JSON parse error: ") + e.what(), line_no);
    }

    if (!obj.is_object() || obj.empty())
      throwLine("Expected a JSON object.", line_no);

    try {
      // Accept both header wrappers: DnaXmlFormat is the original name shared
      // with the XML schema; DnaAdjustmentReport is the adjustment-output
      // header emitted by the JSONL side-channel printer.
      if (obj.contains("DnaXmlFormat") ||
          obj.contains("DnaAdjustmentReport")) {
        parseHeader(obj, ctx, line_no);
      } else if (obj.contains("DnaStatistics")) {
        // Adjustment-report JSONL includes statistics records that are useful
        // to consumers but have no station/measurement input equivalent.
        // Treat them as metadata so .adj.jsonl streams remain re-ingestible.
      } else if (obj.contains("DnaStation")) {
        parseStation(obj, vStations, ctx, line_no);
      } else if (obj.contains("DnaMeasurement")) {
        parseMeasurement(obj, vMeasurements, ctx, line_no);
      } else {
        throwLine("Unrecognised top-level key.", line_no);
      }
    } catch (const XMLInteropException&) {
      throw;  // re-throw as-is
    } catch (const std::exception& e) {
      throwLine(e.what(), line_no);
    }
  }

  // Build success message
  std::ostringstream ss;
  if (ctx.stn_count > 0)
    ss << "Loaded " << ctx.stn_count << " stations";
  if (ctx.msr_count > 0) {
    if (ctx.stn_count > 0)
      ss << ", ";
    ss << "Loaded " << ctx.msr_count << " measurements";
  }
  ctx.message = ss.str();
}

}  // namespace import
}  // namespace dynadjust
