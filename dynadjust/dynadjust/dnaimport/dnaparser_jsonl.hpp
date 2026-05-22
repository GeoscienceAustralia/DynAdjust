//============================================================================
// Name         : dnaparser_jsonl.hpp
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
//                data.  One JSON object per line, using the same element names
//                and hierarchy as the DynaML XML schema.
//============================================================================

#ifndef DYNADJUST_DNAIMPORT_DNAPARSER_JSONL_HPP_
#define DYNADJUST_DNAIMPORT_DNAPARSER_JSONL_HPP_

#include <string>

#include <include/measurement_types/dnameasurement_types.hpp>

namespace dynadjust {
namespace import {

struct JsonlParseContext {
  // Inputs
  std::string default_frame;
  std::string default_epoch;
  bool first_file;
  bool user_supplied_frame;
  bool user_supplied_epoch;
  bool override_input_frame;
  bool prefer_single_x_as_g;

  // Outputs
  UINT32 stn_count;
  UINT32 msr_count;
  UINT32 cluster_id;
  std::string file_epsg;
  std::string file_epoch;
  std::string message;
  bool file_specified_frame;
  bool file_specified_epoch;
};

// Parse a JSONL file containing DynaML-equivalent station and measurement
// data.  Populates vStations and vMeasurements; updates context tallies,
// counts and metadata.
void ParseJsonlFile(const std::string& filename,
                    measurements::vdnaStnPtr* vStations,
                    measurements::vdnaMsrPtr* vMeasurements,
                    JsonlParseContext& ctx);

}  // namespace import
}  // namespace dynadjust

#endif  // DYNADJUST_DNAIMPORT_DNAPARSER_JSONL_HPP_
