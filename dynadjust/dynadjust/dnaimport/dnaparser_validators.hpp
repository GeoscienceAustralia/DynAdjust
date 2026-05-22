//============================================================================
// Name         : dnaparser_validators.hpp
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
//                DynaML parsing.  These are used by the XML (pimpl) parser
//                and can be reused by any future parser (e.g. JSON).
//============================================================================

#ifndef DYNADJUST_DNAIMPORT_DNAPARSER_VALIDATORS_HPP_
#define DYNADJUST_DNAIMPORT_DNAPARSER_VALIDATORS_HPP_

#include <string>

#include <include/measurement_types/dnameasurement_types.hpp>

namespace dynadjust {
namespace import {

// Validates type is non-empty and a recognised measurement letter.
// Returns the type character.
char ValidateMeasurementType(const std::string& type);

// First station must be non-empty.
void ValidateFirst(const std::string& first);

// Second station must be non-empty for types that require it.
void ValidateSecond(const std::string& second, char type_c,
                    const std::string& first);

// Third station must be non-empty for horizontal angles.
void ValidateThird(const std::string& third, char type_c,
                   const std::string& first);

// Value must be non-empty for scalar measurement types.
void ValidateValue(const std::string& value, char type_c,
                   const std::string& first);

// StdDev must be non-empty and positive for scalar measurement types.
void ValidateStdDev(const std::string& stddev, char type_c,
                    const std::string& first);

// Total must be non-empty for cluster/set measurement types (D, G, X, Y).
void ValidateTotal(const std::string& total, char type_c);

// Station name must be non-empty.
void ValidateStationName(const std::string& name);

// Returns true if the measurement type requires a Second station.
bool TypeRequiresSecond(char type_c);

// Returns true if the measurement type is a scalar type that
// requires Value and StdDev fields.
bool TypeRequiresValueStdDev(char type_c);

// Measurement factory — creates the correct subclass from a type
// character.  Increments tallies and measurement count for scalar
// types.  Cluster types (G, X, D, Y) have their tallies updated
// per-component by the caller.
measurements::dnaMsrPtr CreateMeasurement(
    char type_c, PUINT32 cluster_id, UINT32& msr_count,
    measurements::MsrTally& tally, const std::string& frame,
    const std::string& epoch);

// Post-parse finalisation — validates child counts match declared
// totals, propagates parent properties (ignore, scales, type) to
// child elements.  Call after all children (directions / baselines /
// clusterpoints) have been added.
void FinaliseMeasurement(measurements::dnaMsrPtr& msr,
                         bool prefer_single_x_as_g);

// Direction-in-set post-validation.  Called after parsing each
// Directions sub-element to check that ignoring a direction doesn't
// leave the set empty.
void ValidateDirectionInSet(const measurements::dnaMsrPtr& direction_set,
                            bool direction_ignored,
                            UINT32 declared_total);

// Returns "FFF" if empty, otherwise returns input unchanged.
std::string DefaultConstraints(const std::string& constraints);

// Returns "1" if empty, otherwise returns input unchanged.
std::string DefaultScale(const std::string& scale);

// Replaces "&" with "and" in station descriptions.
std::string SanitiseDescription(const std::string& description);

}  // namespace import
}  // namespace dynadjust

#endif  // DYNADJUST_DNAIMPORT_DNAPARSER_VALIDATORS_HPP_
