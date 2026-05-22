//============================================================================
// Name         : dnaparser_validators.cpp
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
// Description  : Shared validation, factory and finalisation functions for
//                DynaML parsing.
//============================================================================

#include <sstream>
#include <string>

#include <dynadjust/dnaimport/dnaparser_validators.hpp>
#include <include/exception/dnaexception.hpp>
#include <include/functions/dnastringfuncs.hpp>

using dynadjust::exception::XMLInteropException;
using dynadjust::measurements::CDnaAngle;
using dynadjust::measurements::CDnaAzimuth;
using dynadjust::measurements::CDnaCoordinate;
using dynadjust::measurements::CDnaDirection;
using dynadjust::measurements::CDnaDirectionSet;
using dynadjust::measurements::CDnaDistance;
using dynadjust::measurements::CDnaGpsBaselineCluster;
using dynadjust::measurements::CDnaGpsPointCluster;
using dynadjust::measurements::CDnaGpsBaseline;
using dynadjust::measurements::CDnaGpsPoint;
using dynadjust::measurements::CDnaHeight;
using dynadjust::measurements::CDnaHeightDifference;

namespace dynadjust {
namespace import {

bool TypeRequiresSecond(char type_c) {
  switch (type_c) {
  case 'A': case 'B': case 'K':
  case 'C': case 'E': case 'M': case 'S':
  case 'D': case 'G': case 'L':
  case 'V': case 'Z': case 'X':
    return true;
  default:
    return false;
  }
}

bool TypeRequiresValueStdDev(char type_c) {
  switch (type_c) {
  case 'A': case 'B': case 'K':
  case 'C': case 'E': case 'M': case 'S':
  case 'D':
  case 'H': case 'R':
  case 'I': case 'J': case 'P': case 'Q':
  case 'L': case 'V': case 'Z':
    return true;
  default:
    return false;
  }
}

char ValidateMeasurementType(const std::string& type) {
  if (type.empty()) {
    throw XMLInteropException("\"Type\" element cannot be empty.", 0);
  }
  char c = type[0];
  switch (c) {
  case 'A': case 'B': case 'C': case 'D': case 'E':
  case 'G': case 'H': case 'I': case 'J': case 'K':
  case 'L': case 'M': case 'P': case 'Q': case 'R':
  case 'S': case 'V': case 'X': case 'Y': case 'Z':
    return c;
  default:
    std::stringstream ss;
    ss << "Unknown measurement type: " << type;
    throw XMLInteropException(ss.str(), 0);
  }
}

void ValidateFirst(const std::string& first) {
  if (first.empty()) {
    throw XMLInteropException("\"First\" element cannot be empty.", 0);
  }
}

void ValidateSecond(const std::string& second, char type_c,
                    const std::string& first) {
  if (TypeRequiresSecond(type_c) && second.empty()) {
    std::stringstream ss;
    ss << "\"Second\" element cannot be empty for measurement type "
       << type_c << " (first station: " << first << ").";
    throw XMLInteropException(ss.str(), 0);
  }
}

void ValidateThird(const std::string& third, char type_c,
                   const std::string& first) {
  if (type_c == 'A' && third.empty()) {
    std::stringstream ss;
    ss << "\"Third\" element cannot be empty for measurement type "
       << type_c << " (first station: " << first << ").";
    throw XMLInteropException(ss.str(), 0);
  }
}

void ValidateValue(const std::string& value, char type_c,
                   const std::string& first) {
  if (TypeRequiresValueStdDev(type_c) && value.empty()) {
    std::stringstream ss;
    ss << "\"Value\" element cannot be empty for measurement type "
       << type_c << " (first station: " << first << ").";
    throw XMLInteropException(ss.str(), 0);
  }
}

void ValidateStdDev(const std::string& stddev, char type_c,
                    const std::string& first) {
  if (!TypeRequiresValueStdDev(type_c)) {
    return;
  }
  if (stddev.empty()) {
    std::stringstream ss;
    ss << "\"StdDev\" element cannot be empty for " << type_c
       << " measurements (first station: " << first << ").";
    throw XMLInteropException(ss.str(), 0);
  }
  if (DoubleFromString<double>(stddev) < PRECISION_1E25) {
    std::stringstream ss;
    ss << "\"StdDev\" element cannot contain a zero or negative value for "
       << type_c << " measurements (first station: " << first << ").";
    throw XMLInteropException(ss.str(), 0);
  }
}

void ValidateTotal(const std::string& total, char type_c) {
  switch (type_c) {
  case 'D': case 'G': case 'X': case 'Y':
    if (total.empty()) {
      throw XMLInteropException("\"Total\" element cannot be empty.", 0);
    }
  }
}

void ValidateStationName(const std::string& name) {
  if (name.empty()) {
    throw XMLInteropException(
        "DnaStation \"Name\" element cannot be empty.", 0);
  }
}

measurements::dnaMsrPtr CreateMeasurement(
    char type_c, PUINT32 cluster_id, UINT32& msr_count,
    measurements::MsrTally& tally, const std::string& frame,
    const std::string& epoch) {
  measurements::dnaMsrPtr msr;

  switch (type_c) {
  case 'A':
    tally.A++;
    msr_count++;
    msr.reset(new CDnaAngle);
    break;
  case 'B':
    tally.B++;
    msr_count++;
    msr.reset(new CDnaAzimuth);
    break;
  case 'C':
    tally.C++;
    msr_count++;
    msr.reset(new CDnaDistance);
    break;
  case 'D':
    // D tallies are incremented per-direction in post_Directions.
    msr.reset(new CDnaDirectionSet(++(*cluster_id)));
    break;
  case 'E':
    tally.E++;
    msr_count++;
    msr.reset(new CDnaDistance);
    break;
  case 'G':
  case 'X':
    // G/X tallies are incremented per-component in GPSBaseline X/Y/Z.
    msr.reset(new CDnaGpsBaselineCluster(++(*cluster_id), frame, epoch));
    break;
  case 'H':
    tally.H++;
    msr_count++;
    msr.reset(new CDnaHeight);
    break;
  case 'I':
    tally.I++;
    msr_count++;
    msr.reset(new CDnaCoordinate);
    break;
  case 'J':
    tally.J++;
    msr_count++;
    msr.reset(new CDnaCoordinate);
    break;
  case 'K':
    tally.K++;
    msr_count++;
    msr.reset(new CDnaAzimuth);
    break;
  case 'L':
    tally.L++;
    msr_count++;
    msr.reset(new CDnaHeightDifference);
    break;
  case 'M':
    tally.M++;
    msr_count++;
    msr.reset(new CDnaDistance);
    break;
  case 'P':
    tally.P++;
    msr_count++;
    msr.reset(new CDnaCoordinate);
    break;
  case 'Q':
    tally.Q++;
    msr_count++;
    msr.reset(new CDnaCoordinate);
    break;
  case 'R':
    tally.R++;
    msr_count++;
    msr.reset(new CDnaHeight);
    break;
  case 'S':
    tally.S++;
    msr_count++;
    msr.reset(new CDnaDistance);
    break;
  case 'V':
    tally.V++;
    msr_count++;
    msr.reset(new CDnaDirection);
    break;
  case 'Y':
    // Y tallies are incremented per-component in Clusterpoint X/Y/Z.
    msr.reset(new CDnaGpsPointCluster(++(*cluster_id), frame, epoch));
    break;
  case 'Z':
    tally.Z++;
    msr_count++;
    msr.reset(new CDnaDirection);
    break;
  }

  return msr;
}

void FinaliseMeasurement(measurements::dnaMsrPtr& msr,
                         bool prefer_single_x_as_g) {
  UINT32 total = 0;
  UINT32 found = 0;

  switch (msr->GetTypeC()) {
  case 'D': {
    total = msr->GetTotal();

    // Ignored empty direction set — skip silently.
    if (total == 0 && msr->GetIgnore()) {
      return;
    }

    found = static_cast<UINT32>(msr->GetDirections_ptr()->size());
    if (total != found) {
      std::stringstream ss;
      ss << "Direction set declares total of " << total
         << " but found " << found << " directions.";
      throw XMLInteropException(ss.str(), 0);
    }

    // Count non-ignored directions.
    found = 0;
    for (CDnaDirection& d : *msr->GetDirections_ptr()) {
      if (d.NotIgnored()) {
        found++;
      }
    }

    if (found == 0 && msr->NotIgnored()) {
      std::stringstream ss;
      ss << "Direction set declares total of " << total
         << " but there aren't any non-ignored directions in the set.";
      throw XMLInteropException(ss.str(), 0);
    }

    msr->SetNonIgnoredDirns(found);

    // Propagate parent properties to child directions.
    for (CDnaDirection& d : *msr->GetDirections_ptr()) {
      d.SetType(msr->GetType());
      d.SetFirst(msr->GetFirst());
      d.SetRecordedTotal(0);
    }
    break;
  }

  case 'G':
  case 'X': {
    total = msr->GetTotal();
    found = static_cast<UINT32>(msr->GetBaselines_ptr()->size());
    if (total != found) {
      std::stringstream ss;
      ss << "GPS baseline cluster declares total of " << total
         << " but found " << found << " baselines.";
      throw XMLInteropException(ss.str(), 0);
    }

    // Promote single-baseline X to G if preferred.
    if (msr->GetTypeC() == 'X' && total < 2 && prefer_single_x_as_g) {
      msr->SetType("G");
    }

    // Propagate parent properties to child baselines.
    for (CDnaGpsBaseline& b : *msr->GetBaselines_ptr()) {
      b.SetIgnore(msr->GetIgnore());
      b.SetRecordedTotal(msr->GetTotal());
      b.SetPscale(msr->GetPscale());
      b.SetLscale(msr->GetLscale());
      b.SetHscale(msr->GetHscale());
      b.SetVscale(msr->GetVscale());
    }
    break;
  }

  case 'Y': {
    total = msr->GetTotal();
    found = static_cast<UINT32>(msr->GetPoints_ptr()->size());
    if (total != found) {
      std::stringstream ss;
      ss << "GPS point cluster declares total of " << total
         << " but found " << found << " clusterpoints.";
      throw XMLInteropException(ss.str(), 0);
    }

    // Propagate parent properties to child points.
    for (CDnaGpsPoint& p : *msr->GetPoints_ptr()) {
      p.SetIgnore(msr->GetIgnore());
      p.SetCoordType(msr->GetCoordType());
      p.SetRecordedTotal(msr->GetTotal());
      p.SetPscale(msr->GetPscale());
      p.SetLscale(msr->GetLscale());
      p.SetHscale(msr->GetHscale());
      p.SetVscale(msr->GetVscale());
    }
    break;
  }
  }
}

void ValidateDirectionInSet(const measurements::dnaMsrPtr& direction_set,
                            bool direction_ignored,
                            UINT32 declared_total) {
  // If the parent set is ignored, no further checks.
  if (direction_set->GetIgnore()) {
    return;
  }

  // If the direction itself is not ignored, no problem.
  if (!direction_ignored) {
    return;
  }

  // Ignoring this direction — check it doesn't leave the set empty.
  UINT32 remaining = direction_set->GetTotal() - 1;
  if (remaining == 0) {
    std::stringstream ss;
    ss << "Direction set declares total of " << declared_total
       << " but there aren't any non-ignored directions in the set.";
    throw XMLInteropException(ss.str(), 0);
  }
}

std::string DefaultConstraints(const std::string& constraints) {
  return constraints.empty() ? "FFF" : constraints;
}

std::string DefaultScale(const std::string& scale) {
  return scale.empty() ? "1" : scale;
}

std::string SanitiseDescription(const std::string& desc) {
  return findandreplace(desc, std::string("&"), std::string("and"));
}

}  // namespace import
}  // namespace dynadjust
