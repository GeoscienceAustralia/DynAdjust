//============================================================================
// Name         : network_data_loader.cpp
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
// Description  : Network Data Loader implementation
//============================================================================

#include "network_data_loader.hpp"

#include <filesystem>
#include <sstream>
#include <system_error>
#include <include/functions/dnatemplatefuncs.hpp>
#include <include/functions/dnastrmanipfuncs.hpp>
#include <include/functions/dnatemplatestnmsrfuncs.hpp>
#include <include/measurement_types/dnastation.hpp>

namespace dynadjust {
namespace networkadjust {

namespace {

bool FileExists(const std::string& filename) {
    if (filename.empty()) {
        return false;
    }

    std::error_code ec;
    const bool exists = std::filesystem::exists(filename, ec);
    return !ec && exists;
}

}  // namespace

NetworkDataLoader::NetworkDataLoader(const project_settings &settings)
    : settings_(settings),
      bst_loader_(std::make_unique<dynadjust::iostreams::BstFile>()),
      bms_loader_(std::make_unique<dynadjust::iostreams::BmsFile>()),
      map_loader_(std::make_unique<dynadjust::iostreams::MapFile>()),
      measurement_processor_(std::make_unique<processors::MeasurementProcessor>(
          settings.a.adjust_mode == SimultaneousMode ? processors::AdjustmentMode::Simultaneous
                                                     : processors::AdjustmentMode::Phased)) {}

bool NetworkDataLoader::LoadCommon(vstn_t *bstBinaryRecords, binary_file_meta_t &bst_meta, vASL *vAssocStnList,
                                   vmsr_t *bmsBinaryRecords, binary_file_meta_t &bms_meta, vUINT32 &v_ISLTemp,
                                   UINT32 &bstn_count, UINT32 &asl_count, UINT32 &bmsr_count, UINT32 &unknownParams,
                                   UINT32 &unknownsCount) {
    if (!LoadStations(bstBinaryRecords, bst_meta, bstn_count)) return false;

    // Apply constraints internally
    if (bstBinaryRecords && !settings_.a.station_constraints.empty()) { ApplyConstraints(*bstBinaryRecords); }

    if (!LoadAssociatedStations(vAssocStnList, v_ISLTemp, asl_count)) return false;

    // Remove invalid stations internally
    if (vAssocStnList) { RemoveInvalidStations(v_ISLTemp, *vAssocStnList); }

    unknownParams = unknownsCount = static_cast<UINT32>(v_ISLTemp.size() * 3);

    if (!LoadMeasurements(bmsBinaryRecords, bms_meta, bmsr_count)) { return false; }

    return true;
}

bool NetworkDataLoader::LoadForPhased(vstn_t *bstBinaryRecords, binary_file_meta_t &bst_meta, vASL *vAssocStnList,
                                      vmsr_t *bmsBinaryRecords, binary_file_meta_t &bms_meta, vvUINT32 &v_ISL,
                                      v_uint32_uint32_map *v_blockStationsMap, vvUINT32 *v_CML, UINT32 &bstn_count,
                                      UINT32 &asl_count, UINT32 &bmsr_count, UINT32 &unknownParams,
                                      UINT32 &unknownsCount, UINT32 &measurementParams, UINT32 &measurementCount,
                                      UINT32 &measurementVarianceCount) {
    vUINT32 v_ISLTemp;
    if (!LoadCommon(bstBinaryRecords, bst_meta, vAssocStnList, bmsBinaryRecords, bms_meta, v_ISLTemp, bstn_count,
                    asl_count, bmsr_count, unknownParams, unknownsCount)) {
        return false;
    }

    // For phased mode, we just need the measurements loaded for segmentation
    // file processing but we don't process them the same way as simultaneous mode
    measurementParams = measurementCount = bmsr_count;
    measurementVarianceCount = bmsr_count;

    // v_ISL is not used in phased mode, but we need to ensure it's initialized
    // to maintain compatibility
    if (v_ISL.empty()) { v_ISL.push_back(v_ISLTemp); }

    return true;
}

bool NetworkDataLoader::LoadForSimultaneous(
    vstn_t *bstBinaryRecords, binary_file_meta_t &bst_meta, vASL *vAssocStnList, vmsr_t *bmsBinaryRecords,
    binary_file_meta_t &bms_meta, vvUINT32 &v_ISL, v_uint32_uint32_map *v_blockStationsMap, vvUINT32 *v_CML,
    UINT32 &bstn_count, UINT32 &asl_count, UINT32 &bmsr_count, UINT32 &unknownParams, UINT32 &unknownsCount,
    UINT32 &measurementParams, UINT32 &measurementCount, UINT32 &measurementVarianceCount, UINT32 *blockCount,
    vvUINT32 *v_JSL, vUINT32 *v_unknownsCount, vUINT32 *v_measurementCount, vUINT32 *v_measurementVarianceCount,
    vUINT32 *v_measurementParams, vUINT32 *v_ContiguousNetList, std::vector<blockMeta_t> *v_blockMeta,
    vvUINT32 *v_parameterStationList, vv_stn_appear *v_paramStnAppearance, v_mat_2d *v_junctionVariances,
    v_mat_2d *v_junctionVariancesFwd) {
    vUINT32 v_ISLTemp;
    if (!LoadCommon(bstBinaryRecords, bst_meta, vAssocStnList, bmsBinaryRecords, bms_meta, v_ISLTemp, bstn_count,
                    asl_count, bmsr_count, unknownParams, unknownsCount)) {
        return false;
    }

    ProcessSimultaneousMode(v_ISLTemp, v_ISL, v_blockStationsMap, unknownParams, unknownsCount);

    ProcessMeasurements(bmsBinaryRecords, bmsr_count, v_CML, measurementParams, measurementCount,
                        measurementVarianceCount);

    // Apply non-measurement removal for simultaneous mode
    if (v_CML && !v_CML->empty()) { RemoveNonMeasurements((*v_CML)[0], *bmsBinaryRecords); }

    // Initialize additional vectors for simultaneous mode
    InitializeSimultaneousModeVectors(
        v_ISL, unknownsCount, measurementParams, measurementCount, measurementVarianceCount, blockCount, v_JSL,
        v_unknownsCount, v_measurementCount, v_measurementVarianceCount, v_measurementParams, v_ContiguousNetList,
        v_blockMeta, v_parameterStationList, v_paramStnAppearance, v_junctionVariances, v_junctionVariancesFwd);

    return true;
}

bool NetworkDataLoader::LoadStations(vstn_t *bstBinaryRecords, binary_file_meta_t &bst_meta, UINT32 &bstn_count) {
    if (!FileExists(settings_.a.bst_file)) {
        throw StationLoadException("Failed to load binary station file: " + settings_.a.bst_file);
    }

    auto result = bst_loader_->LoadWithOptional(settings_.a.bst_file, bstBinaryRecords, bst_meta);
    if (!result) { throw StationLoadException("Failed to load binary station file: " + settings_.a.bst_file); }
    bstn_count = static_cast<UINT32>(*result);
    return true;
}

bool NetworkDataLoader::LoadAssociatedStations(vASL *vAssocStnList, vUINT32 &v_ISLTemp, UINT32 &asl_count) {
    if (!FileExists(settings_.s.asl_file)) {
        throw StationLoadException("Failed to load associated station list file: " + settings_.s.asl_file);
    }

    dynadjust::iostreams::AslFile asl_loader(settings_.s.asl_file);
    auto result = asl_loader.TryLoad();
    if (!result) { throw StationLoadException("Failed to load associated station list file: " + settings_.s.asl_file); }

    // Copy results to output parameters
    *vAssocStnList = std::move(result->stations);
    v_ISLTemp = std::move(result->free_stations);
    asl_count = static_cast<UINT32>(result->count);
    return true;
}

bool NetworkDataLoader::LoadMeasurements(vmsr_t *bmsBinaryRecords, binary_file_meta_t &bms_meta, UINT32 &bmsr_count) {
    if (!FileExists(settings_.a.bms_file)) {
        throw MeasurementLoadException("Failed to load binary measurement file: " + settings_.a.bms_file);
    }

    auto result = bms_loader_->LoadWithOptional(settings_.a.bms_file, bmsBinaryRecords, bms_meta);
    if (!result) { throw MeasurementLoadException("Failed to load binary measurement file: " + settings_.a.bms_file); }
    bmsr_count = static_cast<UINT32>(*result);
    return true;
}

void NetworkDataLoader::ProcessSimultaneousMode(const vUINT32 &v_ISLTemp, vvUINT32 &v_ISL,
                                                v_uint32_uint32_map *v_blockStationsMap, UINT32 &unknownParams,
                                                UINT32 &unknownsCount) {
    v_ISL.clear();
    v_ISL.push_back(v_ISLTemp);

    // Ensure v_blockStationsMap has at least one map
    if (v_blockStationsMap->empty()) {
        throw NetworkLoadException("ProcessSimultaneousMode: blockStationsMap is empty");
    }

    try {
        // Insert station mappings
        for (UINT32 i = 0; i < v_ISL.at(0).size(); ++i) (*v_blockStationsMap)[0][v_ISL.at(0).at(i)] = i;
    } catch (const std::bad_alloc &e) {
        std::stringstream ss;
        ss << "ProcessSimultaneousMode: Could not allocate memory for blockStationsMap: " << e.what();
        throw NetworkLoadException(ss.str());
    }
}

void NetworkDataLoader::ProcessMeasurements(vmsr_t *bmsBinaryRecords, UINT32 bmsr_count, vvUINT32 *v_CML,
                                            UINT32 &measurementParams, UINT32 &measurementCount,
                                            UINT32 &measurementVarianceCount) {
    processors::MeasurementCounts counts;
    auto result = measurement_processor_->ProcessForMode(*bmsBinaryRecords, bmsr_count, *v_CML, counts);

    if (!result) {
        throw MeasurementLoadException(
            "No measurements were found. "
            "If measurements were successfully loaded on import, "
            "ensure that all measurements have not been ignored.");
    }

    measurementParams = measurementCount = counts.measurement_count;
    measurementVarianceCount = counts.measurement_variance_count;
}

void NetworkDataLoader::ApplyConstraints(vstn_t &stations, std::string_view station_map_file) {
    if (settings_.a.station_constraints.empty()) { return; }

    // Load station map
    v_string_uint32_pair station_map;
    std::string map_file{station_map_file};
    if (map_file.empty()) { map_file = settings_.g.input_folder + FOLDER_SLASH + settings_.g.network_name + ".map"; }
    LoadStationMap(&station_map, map_file);

    // Parse constraints from comma-delimited string
    vstring constraint_tokens;
    try {
        SplitDelimitedString<std::string>(settings_.a.station_constraints, std::string(","), &constraint_tokens);
    } catch (...) {
        return;  // Invalid constraint format, silently ignore
    }

    if (apply_discontinuities_) { AddDiscontinuitySites(constraint_tokens, stations); }

    // Process station constraints in pairs (station_name, constraint)
    for (auto it = constraint_tokens.begin(); it != constraint_tokens.end(); ++it) {
        const std::string &station_name = *it;

        // Get the constraint (next token)
        if (++it == constraint_tokens.end()) {
            break;  // No constraint found for this station
        }

        std::string constraint = *it;
        std::transform(constraint.begin(), constraint.end(), constraint.begin(), ::toupper);

        // Find station in map
        auto station_range =
            std::equal_range(station_map.begin(), station_map.end(), station_name, StationNameIDCompareName());

        if (station_range.first == station_range.second) {
            if (apply_discontinuities_) {
                continue;  // Station may have been renamed during discontinuity processing
            }

            std::ostringstream ss;
            ss << "The supplied constraint station '" << station_name << "' is not in the stations map";
            throw ConstraintException(ss.str());
        }

        // Validate constraint
        if (!CDnaStation::IsValidConstraint(constraint)) {
            throw ConstraintException("Invalid station constraint: '" + constraint + "'");
        }

        // Apply constraint to station
        const UINT32 station_index = station_range.first->second;
        if (station_index < stations.size()) {
            snprintf(stations[station_index].stationConst, STN_CONST_WIDTH, "%s", constraint.c_str());
        }
    }
}

void NetworkDataLoader::RemoveNonMeasurements(vUINT32 &measurement_list, const vmsr_t &measurements) {
    if (measurement_list.size() < 2) { return; }

    // Modern C++ approach using lambda instead of custom functor
    auto is_non_measurement = [&measurements](UINT32 index) -> bool {
        return index >= measurements.size() || measurements[index].measStart != xMeas;
    };

    // Remove non-measurements using modern algorithm
    auto new_end = std::remove_if(measurement_list.begin(), measurement_list.end(), is_non_measurement);
    measurement_list.erase(new_end, measurement_list.end());

    // Sort by file order if we have measurements remaining
    if (!measurement_list.empty()) {
        std::sort(measurement_list.begin(), measurement_list.end(), [&measurements](UINT32 a, UINT32 b) -> bool {
            if (a >= measurements.size() || b >= measurements.size()) { return a < b; }
            return measurements[a].fileOrder < measurements[b].fileOrder;
        });
    }
}

void NetworkDataLoader::RemoveInvalidStations(vUINT32 &station_list, const vASL &associated_stations) {
    if (station_list.empty()) { return; }

    // Modern C++ approach using remove_if + erase and lambda
    // Remove stations that are invalid (equivalent to original FALSE check)
    auto new_end =
        std::remove_if(station_list.begin(), station_list.end(), [&associated_stations](UINT32 index) -> bool {
            return index >= associated_stations.size() ||
                   associated_stations[index].Validity() == 0;  // 0 = INVALID_STATION
        });
    station_list.erase(new_end, station_list.end());

    // Sort remaining valid stations by index
    std::sort(station_list.begin(), station_list.end());
}

void NetworkDataLoader::LoadStationMap(const pv_string_uint32_pair station_map, std::string_view map_file) {
    try {
        map_loader_->LoadFile(std::string{map_file}, station_map);
    } catch (const std::runtime_error &e) {
        throw StationLoadException("Failed to load station map file '" + std::string{map_file} + "': " + e.what());
    }
}

void NetworkDataLoader::AddDiscontinuitySites(vstring &constraint_stations, vstn_t &stations) {
    // Sort stations by original name to enable efficient searching
    std::sort(stations.begin(), stations.end(), [](const station_t &lhs, const station_t &rhs) {
        return std::string(lhs.stationNameOrig) < std::string(rhs.stationNameOrig);
    });

    vstring constraint_discontinuity_stations;

    // Process constraint stations in pairs (station_name, constraint_code)
    for (auto const_it = constraint_stations.begin(); const_it != constraint_stations.end(); ++const_it) {
        // Find this constraint station in the station records by original name
        // Create a comparator that works with heterogeneous lookup
        auto comparator = [](const station_t &lhs, const station_t &rhs) {
            return std::string(lhs.stationNameOrig) < std::string(rhs.stationNameOrig);
        };

        // Create a temporary station object for comparison
        station_t temp_station;
        std::strncpy(temp_station.stationNameOrig, const_it->c_str(), STN_NAME_WIDTH);
        temp_station.stationNameOrig[STN_NAME_WIDTH - 1] = '\0';

        auto station_range = std::equal_range(stations.begin(), stations.end(), temp_station, comparator);

        if (station_range.first == station_range.second) {
            // Station not found in records, skip to next pair
            continue;
        }

        // This constraint station is one of the discontinuity sites.
        // Add the current "discontinuity name" (stationName, not stationNameOrig)
        constraint_discontinuity_stations.push_back(station_range.first->stationName);

        // Add the constraint code (next token)
        if (++const_it == constraint_stations.end()) { break; }
        constraint_discontinuity_stations.push_back(*const_it);
    }

    // Add the constraint discontinuity sites to the original constraint list
    constraint_stations.insert(constraint_stations.end(), constraint_discontinuity_stations.begin(),
                               constraint_discontinuity_stations.end());

    // Restore original sort order by regular station name
    std::sort(stations.begin(), stations.end(), [](const station_t &lhs, const station_t &rhs) {
        return std::string(lhs.stationName) < std::string(rhs.stationName);
    });
}

void NetworkDataLoader::InitializeSimultaneousModeVectors(
    const vvUINT32 &v_ISL, UINT32 unknownsCount, UINT32 measurementParams, UINT32 measurementCount,
    UINT32 measurementVarianceCount, UINT32 *blockCount, vvUINT32 *v_JSL, vUINT32 *v_unknownsCount,
    vUINT32 *v_measurementCount, vUINT32 *v_measurementVarianceCount, vUINT32 *v_measurementParams,
    vUINT32 *v_ContiguousNetList, std::vector<blockMeta_t> *v_blockMeta, vvUINT32 *v_parameterStationList,
    vv_stn_appear *v_paramStnAppearance, v_mat_2d *v_junctionVariances, v_mat_2d *v_junctionVariancesFwd) {
    // Only initialize if we're in simultaneous mode
    if (settings_.a.adjust_mode != SimultaneousMode) { return; }

    if (blockCount) { *blockCount = 1; }

    // Initialize simple vectors
    InitializeIfEmpty(v_JSL, 1);
    InitializeIfEmpty(v_paramStnAppearance, 1);
    InitializeIfEmpty(v_junctionVariances, 1);
    InitializeIfEmpty(v_junctionVariancesFwd, 1);

    // Initialize vectors with specific values
    if (v_unknownsCount && v_unknownsCount->empty()) {
        v_unknownsCount->resize(1);
        v_unknownsCount->at(0) = unknownsCount;
    }

    if (v_measurementCount && v_measurementCount->empty()) {
        v_measurementCount->resize(1);
        v_measurementCount->at(0) = measurementCount;
    }

    if (v_measurementVarianceCount && v_measurementVarianceCount->empty()) {
        v_measurementVarianceCount->resize(1);
        v_measurementVarianceCount->at(0) = measurementVarianceCount;
    }

    if (v_measurementParams && v_measurementParams->empty()) {
        v_measurementParams->resize(1);
        v_measurementParams->at(0) = measurementParams;
    }

    if (v_ContiguousNetList && v_ContiguousNetList->empty()) {
        v_ContiguousNetList->resize(1);
        v_ContiguousNetList->at(0) = 0;
    }

    // Initialize block metadata vectors
    if (v_blockMeta && v_blockMeta->empty()) {
        v_blockMeta->resize(1);
        v_blockMeta->at(0)._blockFirst = true;
        v_blockMeta->at(0)._blockLast = true;
        v_blockMeta->at(0)._blockIntermediate = false;
        v_blockMeta->at(0)._blockIsolated = false;
    }

    if (v_parameterStationList && v_parameterStationList->empty()) {
        v_parameterStationList->resize(1);
        // For simultaneous mode, all stations are parameters
        v_parameterStationList->at(0) = v_ISL.at(0);
    }
}

}  // namespace networkadjust
}  // namespace dynadjust
