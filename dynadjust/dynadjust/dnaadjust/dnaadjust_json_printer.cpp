//============================================================================
// Name         : dnaadjust_json_printer.cpp
// Author       : Dale Roberts <dale.o.roberts@gmail.com>
// Copyright    : Copyright 2026 Geoscience Australia
//
//                Licensed under the Apache License, Version 2.0 (the "License");
//                you may not use this file except in compliance with the License.
//                You may obtain a copy of the License at
//
//                http://www.apache.org/licenses/LICENSE-2.0
//
//                Unless required by applicable law or agreed to in writing, software
//                distributed under the License is distributed on an "AS IS" BASIS,
//                WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//                See the License for the specific language governing permissions and
//                limitations under the License.
//
// Description  : JSONL output implementation.  Builds nlohmann::json objects
//                from the same in-memory data consumed by the text printer
//                and emits one JSON object per line.  Key names mirror the
//                DynaML / JSONL input schema so records round-trip through
//                dnaimport.
//============================================================================

#include <cmath>
#include <mutex>
#include <string>

#include <include/config/dnaversion.hpp>
#include <include/functions/dnaiostreamfuncs.hpp>
#include <include/functions/dnastrmanipfuncs.hpp>
#include <include/functions/dnatemplatecalcfuncs.hpp>
#include <include/functions/dnatemplatematrixfuncs.hpp>
#include <include/measurement_types/dnameasurement.hpp>
#include <include/parameters/dnadatum.hpp>
#include <include/parameters/dnaellipsoid.hpp>
#include <include/thirdparty/nlohmann/json.hpp>

#include <dynadjust/dnaadjust/dnaadjust_json_printer.hpp>
#include <dynadjust/dnaadjust/dnaadjust.hpp>

namespace dynadjust {
namespace networkadjust {

namespace {

using nlohmann::json;

// Trim trailing space/\r/\n from a fixed-width station field. Fields are
// right-padded with spaces on disk and null-terminated in memory, so the
// std::string(const char*) conversion plus trimstr yields the clean name.
std::string Trimmed(const char* buf) {
  return trimstr<std::string>(std::string(buf));
}

json HeaderRecord(const CDnaDatum& datum, const char* report_kind) {
  json hdr;
  hdr["type"] = "Adjustment";
  hdr["report"] = report_kind;
  hdr["software"] = std::string("DynAdjust ") + __BINARY_VERSION__;
  hdr["referenceframe"] = datum.GetName();
  hdr["epoch"] = datum.GetEpoch_s();
  // Output wrapper distinct from DnaXmlFormat (used by DynaML XML / JSONL
  // input) so readers can tell an adjustment report apart from input data.
  // The JSONL input parser recognises both wrappers as a header line.
  return json{{"DnaAdjustmentReport", hdr}};
}

json StationIdentity(const it_vstn_t& stn_it) {
  json s;
  s["Name"] = Trimmed(stn_it->stationName);
  s["Constraints"] = Trimmed(stn_it->stationConst);
  // Round-trip StationCoord blocks on .xyz records carry Lat/Lon in
  // DMS-packed form, so the type tag is always LLH for adjustment output.
  s["Type"] = "LLH";
  const std::string desc = Trimmed(stn_it->description);
  if (!desc.empty())
    s["Description"] = desc;
  return s;
}

// Geographic representation of an adjusted Cartesian position, computed once
// per station so OnAdjustedStation does not run CartToGeo twice.
struct AdjustedPosition {
  double x, y, z;
  double lat, lon, height;  // lat/lon in radians, height in metres
};

AdjustedPosition PositionFromContext(
    const matrix_2d* estimates, const UINT32& mat_idx,
    const DynAdjustPrinter::AdjustedStationContext& ctx) {
  // x/y/z read straight from the estimate matrix; lat/lon/h were already
  // computed by PrintAdjStation and threaded through ctx, so CartToGeo
  // does not need to run again.
  AdjustedPosition p;
  p.x = estimates->get(mat_idx, 0);
  p.y = estimates->get(mat_idx + 1, 0);
  p.z = estimates->get(mat_idx + 2, 0);
  p.lat = ctx.lat;
  p.lon = ctx.lon;
  p.height = ctx.height;
  return p;
}

// Lat and Lon come out of the adjustment as radians but are written in DynaML
// XML output as DMS-packed decimal (31.483180 ≡ 31°48'31.80").  Emitting the
// same encoding keeps the JSONL round-trip through SetXAxis/SetYAxis, whose
// DmstoDeg() converts the decimal back to radians.
json InitialCoords(const it_vstn_t& stn_it) {
  json init;
  init["Lat"] = RadtoDms(stn_it->initialLatitude);
  init["Lon"] = RadtoDmsL(stn_it->initialLongitude);
  init["Height"] = stn_it->initialHeight;
  return init;
}

json AdjustedCoords(const AdjustedPosition& p) {
  json adj;
  adj["X"] = p.x;
  adj["Y"] = p.y;
  adj["Z"] = p.z;
  adj["Lat"] = RadtoDms(p.lat);
  adj["Lon"] = RadtoDmsL(p.lon);
  adj["Height"] = p.height;
  return adj;
}

// Round-trippable StationCoord block using the DynaML LLH convention so the
// JSONL input parser can re-ingest an adjustment report as a station file.
json StationCoordLLH(const it_vstn_t& stn_it, const AdjustedPosition& p) {
  json c;
  c["Name"] = Trimmed(stn_it->stationName);
  c["XAxis"] = RadtoDms(p.lat);   // Lat, DMS-packed decimal
  c["YAxis"] = RadtoDmsL(p.lon);  // Lon, DMS-packed decimal
  c["Height"] = p.height;
  return c;
}

json UncertaintyBlockFromCart(const it_vstn_t& stn_it,
                              const matrix_2d& var_cart,
                              const matrix_2d& var_local_with_geoid) {
  // var_local already has geoid variance in (2,2) per AdjustedStationContext
  // contract; ellipse + PU are derived from it directly so the numeric
  // output matches the previous path byte-for-byte.
  double semimajor = 0., semiminor = 0., azimuth = 0.;
  ErrorEllipseParameters<double>(var_local_with_geoid,
                                 semimajor, semiminor, azimuth);
  double hz_pos_u = 0., vt_pos_u = 0.;
  PositionalUncertainty<double>(semimajor, semiminor,
      std::sqrt(var_local_with_geoid.get(2, 2)), hz_pos_u, vt_pos_u);

  (void)stn_it;  // reserved for future fields; currently all values come
                 // from ctx-derived matrices.

  json u;
  u["SE"] = std::sqrt(var_local_with_geoid.get(0, 0));
  u["SN"] = std::sqrt(var_local_with_geoid.get(1, 1));
  u["SU"] = std::sqrt(var_local_with_geoid.get(2, 2));
  u["SemiMajor"] = semimajor;
  u["SemiMinor"] = semiminor;
  u["Orientation"] = azimuth;
  u["HzPosU"] = hz_pos_u;
  u["VtPosU"] = vt_pos_u;
  u["VarianceLocal"] = {
      {var_local_with_geoid.get(0, 0), var_local_with_geoid.get(0, 1),
       var_local_with_geoid.get(0, 2)},
      {var_local_with_geoid.get(1, 0), var_local_with_geoid.get(1, 1),
       var_local_with_geoid.get(1, 2)},
      {var_local_with_geoid.get(2, 0), var_local_with_geoid.get(2, 1),
       var_local_with_geoid.get(2, 2)}};
  u["VarianceCart"] = {
      {var_cart.get(0, 0), var_cart.get(0, 1), var_cart.get(0, 2)},
      {var_cart.get(1, 0), var_cart.get(1, 1), var_cart.get(1, 2)},
      {var_cart.get(2, 0), var_cart.get(2, 1), var_cart.get(2, 2)}};
  return u;
}

// Applies the geoid model variance to var_local[2,2] and re-derives the
// ellipse + positional uncertainty terms.  Used by OnPositionalUncertainty
// to keep .apu.jsonl numerics identical to the pre-refactor behaviour:
// text .apu excludes geoid from ellipse/PU, but the JSON sibling has
// always included it and downstream consumers depend on that.
json UncertaintyBlockWithGeoid(const it_vstn_t& stn_it,
                               const matrix_2d& var_cart,
                               const matrix_2d& var_local_no_geoid) {
  matrix_2d var_local = var_local_no_geoid;  // copy; keeps base var_local clean
  var_local.elementadd(2, 2,
      static_cast<double>(stn_it->geoidSepUnc) * stn_it->geoidSepUnc);
  return UncertaintyBlockFromCart(stn_it, var_cart, var_local);
}

json CorrectionsBlockFromContext(
    const DynAdjustPrinter::StationCorrectionContext& ctx) {
  json c;
  c["dE"] = ctx.cor_e;
  c["dN"] = ctx.cor_n;
  c["dUp"] = ctx.cor_up;
  return c;
}

bool ScalarUsesAngularInputUnits(const char type) {
  switch (type) {
    case 'A':
    case 'B':
    case 'D':
    case 'I':
    case 'J':
    case 'K':
    case 'P':
    case 'Q':
    case 'V':
    case 'Z':
      return true;
    default:
      return false;
  }
}

double ScalarInputValue(const it_vmsr_t& it_msr) {
  return ScalarUsesAngularInputUnits(it_msr->measType)
      ? RadtoDms(it_msr->term1)
      : it_msr->term1;
}

double ScalarInputStdDev(const it_vmsr_t& it_msr) {
  const double stddev = std::sqrt(it_msr->term2);
  return ScalarUsesAngularInputUnits(it_msr->measType)
      ? Seconds(stddev)
      : stddev;
}

void AddObservationEpoch(json& m, const measurement_t& msr) {
  const std::string observation_epoch = Trimmed(msr.observation_epoch);
  if (!observation_epoch.empty())
    m["EpochOfObservation"] = observation_epoch;
}

// Populate the scalar DnaMeasurement fields common to every type.  Used for
// non-cluster types (A/B/C/E/H/I/J/K/L/M/P/Q/R/S/V/Z) where the bms iterator
// points to a single record whose adjusted value is the scalar term1/measAdj.
json ScalarMeasurement(const vstn_t& stations, const it_vmsr_t& it_msr) {
  json m;
  m["Type"] = std::string(1, it_msr->measType);
  AddObservationEpoch(m, *it_msr);
  if (it_msr->ignore)
    m["Ignore"] = true;
  m["First"] = Trimmed(stations.at(it_msr->station1).stationName);
  if (it_msr->measurementStations >= 2)
    m["Second"] = Trimmed(stations.at(it_msr->station2).stationName);
  if (it_msr->measurementStations >= 3)
    m["Third"] = Trimmed(stations.at(it_msr->station3).stationName);
  m["Value"] = ScalarInputValue(it_msr);
  m["StdDev"] = ScalarInputStdDev(it_msr);
  m["Adjusted"] = it_msr->measAdj;
  m["Correction"] = it_msr->measCorr;
  m["AdjustedPrecision"] = it_msr->measAdjPrec;
  m["ResidualPrecision"] = it_msr->residualPrec;
  m["NStat"] = it_msr->NStat;
  m["TStat"] = it_msr->TStat;
  m["PelzerRel"] = it_msr->PelzerRel;
  return m;
}

json DirectionEntry(const vstn_t& stations, const it_vmsr_t& it_dir) {
  json direction;
  direction["Target"] = Trimmed(stations.at(it_dir->station2).stationName);
  direction["Value"] = ScalarInputValue(it_dir);
  direction["StdDev"] = ScalarInputStdDev(it_dir);
  if (it_dir->ignore)
    direction["Ignore"] = true;

  direction["Adjusted"] = it_dir->measAdj;
  direction["Correction"] = it_dir->measCorr;
  direction["AdjustedPrecision"] = it_dir->measAdjPrec;
  direction["ResidualPrecision"] = it_dir->residualPrec;
  direction["NStat"] = it_dir->NStat;
  direction["TStat"] = it_dir->TStat;
  direction["PelzerRel"] = it_dir->PelzerRel;
  return direction;
}

json DirectionSetJson(const vstn_t& stations, const it_vmsr_t& it_set) {
  json m = ScalarMeasurement(stations, it_set);
  json directions = json::array();

  const UINT32 direction_count =
      it_set->vectorCount1 > 0 ? it_set->vectorCount1 - 1 : 0;
  it_vmsr_t it_dir = it_set + 1;
  for (UINT32 dir = 0; dir < direction_count; ++dir, ++it_dir)
    directions.push_back(DirectionEntry(stations, it_dir));

  m["Total"] = direction_count;
  m["Directions"] = directions;
  return m;
}

json ComponentTriple(const it_vmsr_t& it_x, double measurement_t::*field) {
  const it_vmsr_t it_y = it_x + 1;
  const it_vmsr_t it_z = it_x + 2;
  return json{{"X", (*it_x).*field},
              {"Y", (*it_y).*field},
              {"Z", (*it_z).*field}};
}

json CovarianceComponent(const it_vmsr_t& it_cov_x) {
  const it_vmsr_t it_cov_y = it_cov_x + 1;
  const it_vmsr_t it_cov_z = it_cov_x + 2;

  return json{{"m11", it_cov_x->term1},
              {"m12", it_cov_x->term2},
              {"m13", it_cov_x->term3},
              {"m21", it_cov_y->term1},
              {"m22", it_cov_y->term2},
              {"m23", it_cov_y->term3},
              {"m31", it_cov_z->term1},
              {"m32", it_cov_z->term2},
              {"m33", it_cov_z->term3}};
}

// Build one parser-compatible GPSBaseline or Clusterpoint entry from three
// consecutive bms records (X/Y/Z).  Any covariance rows immediately following
// the triplet are attached to the component using the input grammar.
json ClusterComponent(const vstn_t& stations, const it_vmsr_t& it_x) {
  const it_vmsr_t it_y = it_x + 1;
  const it_vmsr_t it_z = it_x + 2;
  json component;

  component["First"] = Trimmed(stations.at(it_x->station1).stationName);
  if (it_x->measType != 'Y')
    component["Second"] = Trimmed(stations.at(it_x->station2).stationName);
  component["X"] = it_x->term1;
  component["Y"] = it_y->term1;
  component["Z"] = it_z->term1;
  component["SigmaXX"] = it_x->term2;
  component["SigmaXY"] = it_y->term2;
  component["SigmaXZ"] = it_z->term2;
  component["SigmaYY"] = it_y->term3;
  component["SigmaYZ"] = it_z->term3;
  component["SigmaZZ"] = it_z->term4;

  const UINT32 covariance_count = it_x->vectorCount2;
  if (covariance_count > 0) {
    json covariances = json::array();
    it_vmsr_t it_cov = it_x + 3;
    for (UINT32 cov = 0; cov < covariance_count; ++cov) {
      covariances.push_back(CovarianceComponent(it_cov));
      it_cov += 3;
    }
    component[it_x->measType == 'Y' ? "PointCovariance" : "GPSCovariance"] =
        covariances;
  }

  return component;
}

json CollapseSingletonArray(const json& values) {
  return values.size() == 1 ? values.at(0) : values;
}

// Build a full G/X/Y cluster record. The iterator must point to the first X
// component of the cluster. Each element is encoded in the same array shape the
// JSONL parser consumes, preserving Total and inter-component covariance rows.
json ClusterMeasurement(const vstn_t& stations, const it_vmsr_t& it_cluster) {
  json m;
  m["Type"] = std::string(1, it_cluster->measType);
  AddObservationEpoch(m, *it_cluster);
  if (it_cluster->ignore)
    m["Ignore"] = true;
  m["First"] = Trimmed(stations.at(it_cluster->station1).stationName);
  // Y-clusters store a single station per point; G/X baselines reference a
  // Second station on the far end of the first baseline. Child entries carry
  // their own station names, which is what the parser uses for X clusters.
  if (it_cluster->measType != 'Y')
    m["Second"] = Trimmed(stations.at(it_cluster->station2).stationName);

  const UINT32 cluster_count = it_cluster->vectorCount1;
  m["Total"] = cluster_count;

  json components = json::array();
  json adjusted = json::array();
  json correction = json::array();
  json adjusted_precision = json::array();
  json nstat = json::array();
  json tstat = json::array();
  json pelzer = json::array();

  it_vmsr_t it_x = it_cluster;
  for (UINT32 cluster_msr = 0; cluster_msr < cluster_count; ++cluster_msr) {
    components.push_back(ClusterComponent(stations, it_x));
    adjusted.push_back(ComponentTriple(it_x, &measurement_t::measAdj));
    correction.push_back(ComponentTriple(it_x, &measurement_t::measCorr));
    adjusted_precision.push_back(
        ComponentTriple(it_x, &measurement_t::measAdjPrec));
    nstat.push_back(ComponentTriple(it_x, &measurement_t::NStat));
    tstat.push_back(ComponentTriple(it_x, &measurement_t::TStat));
    pelzer.push_back(ComponentTriple(it_x, &measurement_t::PelzerRel));

    it_x += 3 + (it_x->vectorCount2 * 3);
  }

  if (it_cluster->measType == 'Y') {
    const std::string coords = Trimmed(it_cluster->coordType);
    if (!coords.empty())
      m["Coords"] = coords;
    m["Clusterpoint"] = components;
  } else {
    m["GPSBaseline"] = components;
  }

  m["Adjusted"] = CollapseSingletonArray(adjusted);
  m["Correction"] = CollapseSingletonArray(correction);
  m["AdjustedPrecision"] = CollapseSingletonArray(adjusted_precision);
  m["NStat"] = CollapseSingletonArray(nstat);
  m["TStat"] = CollapseSingletonArray(tstat);
  m["PelzerRel"] = CollapseSingletonArray(pelzer);
  return m;
}

}  // namespace

DynAdjustJsonPrinter::DynAdjustJsonPrinter(dna_adjust& adjust_instance)
    : DynAdjustPrinter(adjust_instance) {}

DynAdjustJsonPrinter::~DynAdjustJsonPrinter() { CloseStreams(); }

// Serialise the record outside the lock so concurrent hooks can dump in
// parallel; only the stream write is serialised.
static void WriteRecord(DynAdjustJsonPrinter::JsonStream& js,
                        const json& record) {
  std::string line = record.dump();
  line.push_back('\n');
  std::lock_guard<std::mutex> lock(js.mu);
  if (js.stream.is_open())
    js.stream.write(line.data(), static_cast<std::streamsize>(line.size()));
}

void DynAdjustJsonPrinter::OpenStreams() {
  const auto& o = adjust_.projectSettings_.o;

  // Route file_opener failure through SignalExceptionAdjustment so a JSONL
  // open error surfaces via the same channel as text-report open errors.
  // SignalExceptionAdjustment throws, so control never returns past the
  // catch; any already-opened streams are torn down by ~DynAdjustJsonPrinter.
  auto open_with_header = [this](JsonStream& js,
                                 const std::string& path,
                                 const char* report_kind) {
    if (path.empty())
      return;
    try {
      file_opener(js.stream, path);
    } catch (const std::runtime_error& e) {
      adjust_.SignalExceptionAdjustment(e.what(), 0);
    }
    WriteRecord(js, HeaderRecord(adjust_.datum_, report_kind));
  };

  open_with_header(adj_, o._adj_json_file, "adj");
  open_with_header(xyz_, o._xyz_json_file, "xyz");
  open_with_header(apu_, o._apu_json_file, "apu");
  open_with_header(cor_, o._cor_json_file, "cor");
  open_with_header(m2s_, o._m2s_json_file, "m2s");
}

void DynAdjustJsonPrinter::CloseStreams() {
  auto close_if_open = [](JsonStream& js) {
    if (js.stream.is_open())
      js.stream.close();
  };
  close_if_open(adj_);
  close_if_open(xyz_);
  close_if_open(apu_);
  close_if_open(cor_);
  close_if_open(m2s_);
}

void DynAdjustJsonPrinter::OnAdjustedStation(ReportKind kind,
    const it_vstn_t& stn_it, const matrix_2d* estimates,
    const matrix_2d* variances, const UINT32& mat_idx,
    const AdjustedStationContext& ctx) {
  (void)variances;
  const AdjustedPosition p = PositionFromContext(estimates, mat_idx, ctx);

  json station = StationIdentity(stn_it);
  // StationCoord carries the adjusted Lat/Lon/H in the same encoding used
  // by the JSONL input parser, so the emitted record can be re-ingested as
  // initial inputs for a follow-on adjustment.
  station["StationCoord"] = StationCoordLLH(stn_it, p);
  station["Initial"] = InitialCoords(stn_it);
  station["Adjusted"] = AdjustedCoords(p);
  station["Uncertainty"] =
      UncertaintyBlockFromCart(stn_it, ctx.var_cart, ctx.var_local);

  const json record{{"DnaStation", station}};

  // Route based on the caller-declared kind. Works under staged-mode where
  // the text sink is a stringstream rather than the real file stream.
  switch (kind) {
    case ReportKind::kAdj: WriteRecord(adj_, record); break;
    case ReportKind::kXyz: WriteRecord(xyz_, record); break;
  }
}

void DynAdjustJsonPrinter::OnAdjustedMeasurement(const it_vmsr_t& it_msr) {
  // Cluster types (G single baseline, X baseline cluster, Y point cluster)
  // contain one or more X/Y/Z triplets plus covariance rows. Emit one JSONL
  // DnaMeasurement for the full cluster. Everything else is a scalar.
  const bool is_cluster =
      it_msr->measType == 'G' || it_msr->measType == 'X' ||
      it_msr->measType == 'Y';
  const json body =
      is_cluster ? ClusterMeasurement(adjust_.bstBinaryRecords_, it_msr)
                 : (it_msr->measType == 'D'
                        ? DirectionSetJson(adjust_.bstBinaryRecords_, it_msr)
                        : ScalarMeasurement(adjust_.bstBinaryRecords_,
                                            it_msr));
  WriteRecord(adj_, json{{"DnaMeasurement", body}});
}

void DynAdjustJsonPrinter::OnPositionalUncertainty(const it_vstn_t& stn_it,
    const matrix_2d* variances, const UINT32& mat_idx,
    const PositionalUncertaintyContext& ctx) {
  (void)variances;
  (void)mat_idx;
  json station = StationIdentity(stn_it);
  station["Uncertainty"] =
      UncertaintyBlockWithGeoid(stn_it, ctx.var_cart, ctx.var_local);
  WriteRecord(apu_, json{{"DnaStation", station}});
}

void DynAdjustJsonPrinter::OnStationCorrection(const UINT32& block,
    const it_vstn_t& stn_it, const matrix_2d* estimates,
    const UINT32& mat_idx, const StationCorrectionContext& ctx) {
  (void)estimates;
  (void)mat_idx;
  // Corrections are computed per-block; drop the record if the block index
  // is out of range rather than emitting a malformed DnaStation with no
  // Corrections payload.
  if (block >= adjust_.v_originalStations_.size())
    return;

  json station = StationIdentity(stn_it);
  station["Initial"] = InitialCoords(stn_it);
  station["Corrections"] = CorrectionsBlockFromContext(ctx);
  WriteRecord(cor_, json{{"DnaStation", station}});
}

void DynAdjustJsonPrinter::OnM2SRecord(const it_vstn_t& stn_it,
    MsrTally& tally) {
  json station = StationIdentity(stn_it);

  // Enumerate the 20 canonical MsrTally fields directly rather than going
  // through FillMsrList + MeasurementCount's switch-per-character.
  station["MsrCounts"] = json{
      {"A", tally.A}, {"B", tally.B}, {"C", tally.C}, {"D", tally.D},
      {"E", tally.E}, {"G", tally.G}, {"H", tally.H}, {"I", tally.I},
      {"J", tally.J}, {"K", tally.K}, {"L", tally.L}, {"M", tally.M},
      {"P", tally.P}, {"Q", tally.Q}, {"R", tally.R}, {"S", tally.S},
      {"V", tally.V}, {"X", tally.X}, {"Y", tally.Y}, {"Z", tally.Z}};
  // totalCount is populated by CreateMsrToStnTally before this hook fires,
  // so the cached value is current.
  station["TotalMeasurements"] = tally.totalCount;
  WriteRecord(m2s_, json{{"DnaStation", station}});
}

void DynAdjustJsonPrinter::OnStatistics() {
  json stats;
  // When _adj_stat_iteration is enabled, PrintStatistics also fires after
  // each iteration before the terminal call, so include the iteration
  // number with each record; consumers can pick the largest for the final
  // outcome and reconstruct the convergence history from the rest.
  stats["iteration"] = adjust_.CurrentIteration();
  stats["unknown_parameters"] = adjust_.unknownParams_;
  stats["measurement_params"] = adjust_.measurementParams_;
  stats["potential_outliers"] = adjust_.potentialOutlierCount_;
  stats["dof"] = adjust_.degreesofFreedom_;
  stats["chisq"] = adjust_.chiSquared_;
  stats["sigma_zero"] = adjust_.sigmaZero_;
  stats["global_pelzer"] = adjust_.globalPelzerReliability_;
  stats["chisq_lower"] = adjust_.chiSquaredLowerLimit_;
  stats["chisq_upper"] = adjust_.chiSquaredUpperLimit_;
  stats["confidence_interval"] =
      adjust_.projectSettings_.a.confidence_interval;

  const char* result = "unknown";
  if (adjust_.degreesofFreedom_ < 1)
    result = "no_redundancy";
  else {
    switch (adjust_.passFail_) {
      case test_stat_pass:
        result = "passed";
        break;
      case test_stat_warning:
        result = "warning";
        break;
      case test_stat_fail:
        result = "failed";
        break;
    }
  }
  stats["chisq_test"] = result;
  WriteRecord(adj_, json{{"DnaStatistics", stats}});

  // Flush every open stream so a subsequent crash or signal does not lose
  // JSONL lines buffered by the C++ stream layer.  We flush on every
  // statistics fire (terminal + per-iteration when enabled) — flush is
  // cheap relative to the adjustment iteration that just completed.
  for (JsonStream* js : {&adj_, &xyz_, &apu_, &cor_, &m2s_}) {
    std::lock_guard<std::mutex> lock(js->mu);
    if (js->stream.is_open())
      js->stream.flush();
  }
}

}  // namespace networkadjust
}  // namespace dynadjust
