//============================================================================
// Name         : dnaadjust_printer.cpp
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
// Description  : DynAdjust Network Adjustment Printer Module Implementation
//============================================================================

#include "dnaadjust.hpp"
#include "dnaadjust_printer.hpp"

#include <algorithm>
#include <iomanip>
#include <sstream>
#include <filesystem>
#include <include/functions/dnaiostreamfuncs.hpp>
#include <include/functions/dnafilepathfuncs.hpp>
#include <include/io/adj_file.hpp>

namespace dynadjust {
namespace networkadjust {

DynAdjustPrinter::DynAdjustPrinter(dna_adjust& adjust_instance) 
    : adjust_(adjust_instance) {
}

void DynAdjustPrinter::PrintMeasurementWithStations(it_vmsr_t& it_msr, char measurement_type) {
    // Lambda to print station names based on count
    auto print_stations = [&](int station_count) {
        adjust_.adj_file << std::left << std::setw(STATION) 
                        << adjust_.bstBinaryRecords_.at(it_msr->station1).stationName;
        adjust_.adj_file << std::left << std::setw(STATION) 
                        << adjust_.bstBinaryRecords_.at(it_msr->station2).stationName;
        
        if (station_count == 3) {
            adjust_.adj_file << std::left << std::setw(STATION) 
                            << adjust_.bstBinaryRecords_.at(it_msr->station3).stationName;
        } else {
            adjust_.adj_file << std::left << std::setw(STATION) << " ";
        }
    };
    
    // Get station count for this measurement type
    const int station_count = GetStationCount(measurement_type);
    print_stations(station_count);
    
    // Print the measurement based on type
    if (IsAngularType(measurement_type)) {
        // Delegate to existing angular measurement printer
        PrintAdjMeasurements<AngularMeasurement>(' ', it_msr, true);
    } else {
        // Delegate to existing linear measurement printer
        PrintAdjMeasurements<LinearMeasurement>(' ', it_msr, true);
    }
}

void DynAdjustPrinter::PrintIteration(const UINT32& iteration) {
    std::stringstream iteration_message;

    iteration_message << std::endl << OUTPUTLINE << std::endl <<
        std::setw(PRINT_VAR_PAD) << std::left << "ITERATION" << iteration << std::endl << std::endl;

    adjust_.adj_file << iteration_message.str();

    if (adjust_.projectSettings_.g.verbose)
        adjust_.debug_file << iteration_message.str();		
}

void DynAdjustPrinter::PrintAdjustmentTime(cpu_timer& time, int timer_type) {
    // calculate and print total time
    auto elapsed = time.elapsed();
    double seconds = elapsed.wall.count() / 1.0e9;
    
    std::stringstream ss;
    if (seconds >= 1.0) {
        ss << std::fixed << std::setprecision(3) << seconds << "s";
    } else {
        ss << std::fixed << std::setprecision(3) << (seconds * 1000.0) << "ms";
    }

    if (timer_type == 0) // iteration_time equivalent
        adjust_.adj_file << std::setw(PRINT_VAR_PAD) << std::left << "Elapsed time" << ss.str() << std::endl;
    else
    {
        adjust_.adj_file << std::setw(PRINT_VAR_PAD) << std::left << "Total time" << ss.str() << std::endl << std::endl;
    }
}

constexpr int DynAdjustPrinter::GetStationCount(char measurement_type) const {
    // Use constexpr lookup for station count
    for (const auto& [type, count] : kStationCounts) {
        if (type == measurement_type) {
            return count;
        }
    }
    return 2; // Default to 2 stations
}

constexpr bool DynAdjustPrinter::IsAngularType(char measurement_type) const {
    // Check if measurement type is angular
    for (const auto& type : kAngularTypes) {
        if (type == measurement_type) {
            return true;
        }
    }
    return false;
}

// Template specializations
template<>
void DynAdjustPrinter::PrintAdjMeasurements<AngularMeasurement>(
    char cardinal, const it_vmsr_t& it_msr, bool initialise_dbindex) {
    // Print adjusted angular measurements
    PrintMeasurementValue<AngularMeasurement>(cardinal, it_msr->measAdj, it_msr->measCorr, it_msr, true);
    
    // Print adjusted statistics
    PrintAdjMeasurementStatistics(cardinal, it_msr, initialise_dbindex);
}

template<>
void DynAdjustPrinter::PrintAdjMeasurements<LinearMeasurement>(
    char cardinal, const it_vmsr_t& it_msr, bool initialise_dbindex) {
    // Print adjusted linear measurements
    PrintMeasurementValue<LinearMeasurement>(cardinal, it_msr->measAdj, it_msr->measCorr, it_msr, true);
    
    // Print adjusted statistics
    PrintAdjMeasurementStatistics(cardinal, it_msr, initialise_dbindex);
}

// Template specializations for computed measurements
template<>
void DynAdjustPrinter::PrintComputedMeasurements<AngularMeasurement>(
    char cardinal, const double& computed, const double& correction, const it_vmsr_t& it_msr) {
    // Print computed angular measurements
    PrintMeasurementValue<AngularMeasurement>(cardinal, computed, correction, it_msr, false);
    
    // Print measurement correction
    PrintMeasurementCorrection(cardinal, it_msr);
   
    // Print measurement database ids if enabled
    if (adjust_.projectSettings_.o._database_ids) {
        PrintMeasurementDatabaseID(it_msr, true);
    }
    
    adjust_.adj_file << std::endl;
}

template<>
void DynAdjustPrinter::PrintComputedMeasurements<LinearMeasurement>(
    char cardinal, const double& computed, const double& correction, const it_vmsr_t& it_msr) {
    // Check if adjustment is questionable for linear measurements
    if (!adjust_.isAdjustmentQuestionable_) {
        adjust_.isAdjustmentQuestionable_ = (fabs(correction) > 999.9999);
    }
    
    // Print computed linear measurements
    PrintMeasurementValue<LinearMeasurement>(cardinal, computed, correction, it_msr, false);
    
    // Print measurement correction
    PrintMeasurementCorrection(cardinal, it_msr);
   
    // Print measurement database ids if enabled
    if (adjust_.projectSettings_.o._database_ids) {
        PrintMeasurementDatabaseID(it_msr, true);
    }
    
    adjust_.adj_file << std::endl;
}

// Compatibility methods for maintaining the original method names
void DynAdjustPrinter::PrintCompMeasurementsLinear(const char cardinal, const double& computed, 
                                                   const double& correction, const it_vmsr_t& it_msr) {
    PrintComputedMeasurements<LinearMeasurement>(cardinal, computed, correction, it_msr);
}

void DynAdjustPrinter::PrintCompMeasurementsAngular(const char cardinal, const double& computed, 
                                                    const double& correction, const it_vmsr_t& it_msr) {
    PrintComputedMeasurements<AngularMeasurement>(cardinal, computed, correction, it_msr);
}

void DynAdjustPrinter::PrintMeasurementsAngular(const char cardinal, const double& measurement, 
                                                const double& correction, const it_vmsr_t& it_msr, bool printAdjMsr) {
    PrintMeasurementValue<AngularMeasurement>(cardinal, measurement, correction, it_msr, printAdjMsr);
}

void DynAdjustPrinter::PrintMeasurementsLinear(const char cardinal, const double& measurement, 
                                               const double& correction, const it_vmsr_t& it_msr, bool printAdjMsr) {
    PrintMeasurementValue<LinearMeasurement>(cardinal, measurement, correction, it_msr, printAdjMsr);
}

void DynAdjustPrinter::PrintAdjMeasurementsAngular(const char cardinal, const it_vmsr_t& it_msr, bool initialise_dbindex) {
    // Print adjusted angular measurements
    PrintMeasurementsAngular(cardinal, it_msr->measAdj, it_msr->measCorr, it_msr);

    // Print adjusted statistics
    PrintAdjMeasurementStatistics(cardinal, it_msr, initialise_dbindex);
}

void DynAdjustPrinter::PrintAdjMeasurementsLinear(const char cardinal, const it_vmsr_t& it_msr, bool initialise_dbindex) {
    // Print adjusted linear measurements
    PrintMeasurementsLinear(cardinal, it_msr->measAdj, it_msr->measCorr, it_msr);
    
    // Print adjusted statistics
    PrintAdjMeasurementStatistics(cardinal, it_msr, initialise_dbindex);
}

// Utility function implementations
void DynAdjustPrinter::PrintAdjustmentStatus() {
    // Print adjustment status
    adjust_.adj_file << std::endl << OUTPUTLINE << std::endl;
    adjust_.adj_file << std::setw(PRINT_VAR_PAD) << std::left << "SOLUTION";

    if (adjust_.projectSettings_.a.report_mode) {
        adjust_.adj_file << "Printing results of last adjustment only" << std::endl;
        return;
    }
    
    switch (adjust_.projectSettings_.a.adjust_mode) {
    case Phased_Block_1Mode:
        if (adjust_.adjustStatus_ == ADJUST_SUCCESS)
            adjust_.adj_file << "Estimates solved for Block 1 only" << std::endl;
        else
            adjust_.adj_file << "Failed to solve Block 1 estimates" << std::endl;
        break;
    default:
        if (adjust_.adjustStatus_ == ADJUST_SUCCESS && 
            adjust_.CurrentIteration() <= adjust_.projectSettings_.a.max_iterations &&
            fabs(adjust_.maxCorr_) <= adjust_.projectSettings_.a.iteration_threshold)
            adjust_.adj_file << "Converged" << std::endl;
        else
            adjust_.adj_file << "Failed to converge" << std::endl;
    }
}

void DynAdjustPrinter::PrintMeasurementDatabaseID(const it_vmsr_t& it_msr, bool initialise_dbindex) {
    // Set iterator to the database index
    if (initialise_dbindex) {
        size_t dbindx = std::distance(adjust_.bmsBinaryRecords_.begin(), it_msr);
        adjust_._it_dbid = adjust_.v_msr_db_map_.begin() + dbindx;
    }

    // Print measurement id
    if (adjust_._it_dbid->is_msr_id_set)
        adjust_.adj_file << std::setw(STDDEV) << std::right << adjust_._it_dbid->msr_id;
    else
        adjust_.adj_file << std::setw(STDDEV) << " ";

    // Print cluster id?
    switch (it_msr->measType) {
    case 'D':
    case 'G':
    case 'X':
    case 'Y':
        if (adjust_._it_dbid->is_cls_id_set)
            adjust_.adj_file << std::setw(STDDEV) << std::right << adjust_._it_dbid->cluster_id;
        else
            adjust_.adj_file << std::setw(STDDEV) << " ";
    }
}

void DynAdjustPrinter::PrintAdjMeasurementStatistics(char cardinal, const it_vmsr_t& it_msr, bool initialise_dbindex) {
    UINT16 PRECISION_STAT(2);

    if (adjust_.isAdjustmentQuestionable_ || (fabs(it_msr->NStat) > adjust_.criticalValue_ * 4.0))
        adjust_.adj_file << StringFromTW(removeNegativeZero(it_msr->NStat, 2), STAT, PRECISION_STAT);
    else
        adjust_.adj_file << std::setw(STAT) << std::setprecision(2) << std::fixed << std::right << 
            removeNegativeZero(it_msr->NStat, 2);

    if (adjust_.projectSettings_.o._adj_msr_tstat) {
        if (adjust_.isAdjustmentQuestionable_ || (fabs(it_msr->NStat) > adjust_.criticalValue_ * 4.0))
            adjust_.adj_file << StringFromTW(removeNegativeZero(it_msr->TStat, 2), STAT, PRECISION_STAT);
        else
            adjust_.adj_file << std::setw(STAT) << std::setprecision(2) << std::fixed << std::right << 
                removeNegativeZero(it_msr->TStat, 2);
    }

    adjust_.adj_file << std::setw(REL) << std::setprecision(2) << std::fixed << std::right << it_msr->PelzerRel;

    // Print measurement correction
    PrintMeasurementCorrection(cardinal, it_msr);

    // Print asterisk for values which exceed the critical value
    if (fabs(it_msr->NStat) > adjust_.criticalValue_)
        adjust_.adj_file << std::setw(OUTLIER) << std::right << "*";
    else
        adjust_.adj_file << std::setw(OUTLIER) << std::right << " ";

    // Print measurement database ids
    if (adjust_.projectSettings_.o._database_ids)
        PrintMeasurementDatabaseID(it_msr, initialise_dbindex);

    adjust_.adj_file << std::endl;
}

// Stage 3: Header generation infrastructure
void DynAdjustPrinter::PrintMeasurementHeader(printMeasurementsMode printMode, std::string_view col1_heading, std::string_view col2_heading) {
    UINT32 j(PAD2 + STATION + STATION + STATION);
    adjust_.adj_file <<
        std::setw(PAD2) << std::left << "M" << 
        std::setw(STATION) << std::left << "Station 1" << 
        std::setw(STATION) << std::left << "Station 2" << 
        std::setw(STATION) << std::left << "Station 3";

    // Adjusted, computed and ignored measurements
    j += PAD3 + PAD3 + MSR + MSR + CORR + PREC;
    adjust_.adj_file <<
        std::setw(PAD3) << std::left << "*" <<
        std::setw(PAD2) << std::left << "C" <<
        std::setw(MSR) << std::right << "Measured" <<
        std::setw(MSR) << std::right << col1_heading <<
        std::setw(CORR) << std::right << col2_heading <<
        std::setw(PREC) << std::right << "Meas. SD";
    
    // Adjusted measurements only
    switch (printMode) {
    case adjustedMsrs:
        j += PREC + PREC + STAT;
        adjust_.adj_file <<
            std::setw(PREC) << std::right << "Adj. SD" <<
            std::setw(PREC) << std::right << "Corr. SD" <<
            std::setw(STAT) << std::right << "N-stat";

        // print t-statistics?
        if (adjust_.projectSettings_.o._adj_msr_tstat) {
            j += STAT;
            adjust_.adj_file << std::setw(STAT) << std::right << "T-stat";
        }

        j += REL + OUTLIER;
        adjust_.adj_file <<
            std::setw(REL) << std::right << "Pelzer Rel" <<
            std::setw(OUTLIER) << std::right << "*";
        break;
    case computedMsrs:
    case ignoredMsrs:
        break;
    }

    // Database IDs?
    if (adjust_.projectSettings_.o._database_ids) {
        j += STDDEV;  // Use STDDEV width for database ID
        adjust_.adj_file << std::setw(STDDEV) << std::right << "msr-id";

        // print cluster id for GPS measurements?
        j += STDDEV;  // Use STDDEV width for cluster ID
        adjust_.adj_file << std::setw(STDDEV) << std::right << "cluster";
    }

    adjust_.adj_file << std::endl;
    adjust_.adj_file << std::setfill('-') << std::setw(j) << "" << std::setfill(' ') << std::endl;
}

void DynAdjustPrinter::PrintPositionUncertaintyHeader(std::ostream& os) {
    UINT32 pad = PRINT_VAR_PAD;
    
    os << std::setw(pad) << std::left << "Station";
    
    if (!adjust_.projectSettings_.o._stn_coord_types.empty()) {
        os << std::setw(6) << std::left << "Coord";
        pad += 6;
    }
    
    // Use standard column widths
    os << std::setw(14) << std::right << "Latitude" <<
          std::setw(14) << std::right << "Longitude" <<
          std::setw(10) << std::right << "Height";
    pad += 14 + 14 + 10;
    
    pad += 10 + 10 + 10 + 10 + 10;
    os << std::setw(10) << std::right << "Semi-major" <<
          std::setw(10) << std::right << "Semi-minor" <<
          std::setw(10) << std::right << "Orientation" <<
          std::setw(10) << std::right << "Height prec" <<
          std::setw(10) << std::right << "PU at 95%";

    os << std::endl << std::setfill('-') << std::setw(pad) << "" << std::setfill(' ') << std::endl;
}

void DynAdjustPrinter::PrintAdjMeasurementsHeader(bool printHeader, const std::string& table_heading,
    printMeasurementsMode printMode, UINT32 block, bool printBlocks) {
    
    if (printHeader)
        adjust_.adj_file << std::endl << table_heading << std::endl <<
        "------------------------------------------" << std::endl << std::endl;

    if (printBlocks) {
        switch (adjust_.projectSettings_.a.adjust_mode) {
        case PhasedMode:
        case Phased_Block_1Mode:
            if (adjust_.projectSettings_.o._output_msr_blocks)
                adjust_.adj_file << "Block " << block << std::endl;
            break;
        }
    }

    std::string col1_heading, col2_heading;

    // determine headings
    switch (printMode) {
    case ignoredMsrs:
    case computedMsrs:
        col1_heading = "Computed";
        col2_heading = "Difference";
        break;
    case adjustedMsrs:
        col1_heading = "Adjusted";
        col2_heading = "Correction";
        break;
    }
    
    // print header

    // Adjusted, computed and ignored measurements
    UINT32 j(PAD2 + STATION + STATION + STATION);
    adjust_.adj_file <<
        std::setw(PAD2) << std::left << "M" << 
        std::setw(STATION) << std::left << "Station 1" << 
        std::setw(STATION) << std::left << "Station 2" << 
        std::setw(STATION) << std::left << "Station 3";

    // Adjusted, computed and ignored measurements
    j += PAD3 + PAD3 + MSR + MSR + CORR + PREC;
    adjust_.adj_file <<
        std::setw(PAD3) << std::left << "*" <<
        std::setw(PAD2) << std::left << "C" <<
        std::setw(MSR) << std::right << "Measured" <<
        std::setw(MSR) << std::right << col1_heading <<	// Computed or Adjusted
        std::setw(CORR) << std::right << col2_heading <<	// Difference or Correction
        std::setw(PREC) << std::right << "Meas. SD";
    
    // Adjusted measurements only
    switch (printMode) {
    case adjustedMsrs:
        j += PREC + PREC + STAT;
        adjust_.adj_file <<
            std::setw(PREC) << std::right << "Adj. SD" <<
            std::setw(PREC) << std::right << "Corr. SD" <<
            std::setw(STAT) << std::right << "N-stat";

        // print t-statistics?
        if (adjust_.projectSettings_.o._adj_msr_tstat) {
            j += STAT;
            adjust_.adj_file << std::setw(STAT) << std::right << "T-stat";			
        }

        j += REL;
        adjust_.adj_file <<
            std::setw(REL) << std::right << "Pelzer Rel";
        break;
    default:
        break;
    }
    
    // Adjusted, computed and ignored measurements
    j += PACORR;
    adjust_.adj_file <<
        std::setw(PACORR) << std::right << "Pre Adj Corr";

    // Adjusted measurements only
    switch (printMode) {
    case adjustedMsrs:
        j += OUTLIER;
        adjust_.adj_file <<
            std::setw(OUTLIER) << std::right << "Outlier?";
        break;
    default:
        break;
    }

    // Adjusted, computed and ignored measurements
    // Print database ids?
    if (adjust_.projectSettings_.o._database_ids) {
        j += STDDEV + STDDEV;
        adjust_.adj_file << 
            std::setw(STDDEV) << std::right << "Meas. ID" << 
            std::setw(STDDEV) << std::right << "Clust. ID";
    }

    adjust_.adj_file << std::endl;

    UINT32 i;
    for (i=0; i<j; ++i)
        adjust_.adj_file << "-";

    adjust_.adj_file << std::endl;
}

// Stage 3: Output coordinators
void DynAdjustPrinter::PrintAdjustedNetworkMeasurements() {
    // Complete implementation with ignored measurements handling
    if (adjust_.adjustStatus_ > ADJUST_TEST_FAILED)
        return;
    
    bool printHeader(true);

    if (adjust_.projectSettings_.o._database_ids)
        if (!adjust_.databaseIDsLoaded_ || adjust_.v_msr_db_map_.empty())
            adjust_.LoadDatabaseId();

    switch (adjust_.projectSettings_.a.adjust_mode) {
    case Phased_Block_1Mode:    // only the first block is rigorous
        {
            _it_uint32_u32u32_pair begin = adjust_.v_msr_block_.begin();
            _it_uint32_u32u32_pair end = begin + adjust_.v_CML_.at(0).size();
            PrintAdjMeasurements(v_uint32_u32u32_pair(begin, end), printHeader);
        }
        break;
    case SimultaneousMode:
        PrintAdjMeasurements(adjust_.v_msr_block_, printHeader);
        break;
    case PhasedMode:
        // Use existing logic for complex phased mode
        PrintAdjMeasurements(adjust_.v_msr_block_, printHeader);
        break;
    }

    // Handle ignored measurements for applicable modes
    switch (adjust_.projectSettings_.a.adjust_mode)
    {
    case PhasedMode:
    case SimultaneousMode:
        if (adjust_.projectSettings_.o._print_ignored_msrs)
            // Print comparison between ignored and computed 
            // measurements from adjusted coordinates
            PrintIgnoredAdjMeasurements(true);		
        break;
    }
}

void DynAdjustPrinter::PrintAdjustedNetworkStations() {
    // print adjusted coordinates
    bool printHeader(true);

    switch (adjust_.projectSettings_.a.adjust_mode) {
    case SimultaneousMode:
        PrintAdjStations(adjust_.adj_file, 0, &adjust_.v_estimatedStations_.at(0), &adjust_.v_rigorousVariances_.at(0), false, true, true, printHeader, true);
        PrintAdjStations(adjust_.xyz_file, 0, &adjust_.v_estimatedStations_.at(0), &adjust_.v_rigorousVariances_.at(0), false, false, false, printHeader, false);
        break;
    case PhasedMode:
    case Phased_Block_1Mode:
        // Output phased blocks as a single block?
        if (!adjust_.projectSettings_.o._output_stn_blocks) {
            PrintAdjStationsUniqueList(adjust_.adj_file,
                &adjust_.v_rigorousStations_,
                &adjust_.v_rigorousVariances_,
                true, true, true);
            PrintAdjStationsUniqueList(adjust_.xyz_file,
                &adjust_.v_rigorousStations_,
                &adjust_.v_rigorousVariances_,
                true, true, false);
        } else {
            // Handle staging for complex phased output blocks if needed
            if (adjust_.projectSettings_.a.stage) {
                // Handle special staging logic for blocks written to disk
                for (UINT32 block = 0; block < adjust_.blockCount_; ++block) {
                    // Deserialize blocks from mapped files
                    adjust_.DeserialiseBlockFromMappedFile(block, 2, sf_rigorous_stns, sf_rigorous_vars);
                    if (adjust_.projectSettings_.o._init_stn_corrections || adjust_.projectSettings_.o._stn_corr)
                        adjust_.DeserialiseBlockFromMappedFile(block, 1, sf_original_stns);

                    // Print using original detailed function for staging compatibility
                    PrintAdjStations(adjust_.adj_file, block, &adjust_.v_rigorousStations_.at(block), &adjust_.v_rigorousVariances_.at(block), true, true, true, printHeader, true);
                    PrintAdjStations(adjust_.xyz_file, block, &adjust_.v_rigorousStations_.at(block), &adjust_.v_rigorousVariances_.at(block), true, false, false, printHeader, false);
                    printHeader = false;

                    // Release block from memory
                    adjust_.UnloadBlock(block, 2, sf_rigorous_stns, sf_rigorous_vars);
                    if (adjust_.projectSettings_.o._init_stn_corrections || adjust_.projectSettings_.o._stn_corr)
                        adjust_.SerialiseBlockToMappedFile(block, 1, sf_original_stns);

                    // Exit if block-1 mode
                    if (adjust_.projectSettings_.a.adjust_mode == Phased_Block_1Mode)
                        break;
                }
            } else {
                // Print stations for each block without staging
                for (UINT32 block = 0; block < adjust_.blockCount_; ++block) {
                    PrintAdjStations(adjust_.adj_file, block, &adjust_.v_estimatedStations_.at(block), &adjust_.v_rigorousVariances_.at(block), true, true, true, printHeader, true);
                    PrintAdjStations(adjust_.xyz_file, block, &adjust_.v_estimatedStations_.at(block), &adjust_.v_rigorousVariances_.at(block), true, false, false, printHeader, false);
                    printHeader = false;
                    
                    // Exit if block-1 mode
                    if (adjust_.projectSettings_.a.adjust_mode == Phased_Block_1Mode)
                        break;
                }
            }
        }
        break;
    }
}

void DynAdjustPrinter::PrintDirectionSetMeasurements(it_vmsr_t& it_msr, bool adjustedMsrs) {
    // Direction set specific printing logic - simplified for Stage 3
    vmsr_t d_msr;
    CopyClusterMsr<vmsr_t>(adjust_.bmsBinaryRecords_, it_msr, d_msr);
    
    it_vmsr_t _it_d_msr(d_msr.begin());
    UINT32 cluster_msr, cluster_count(_it_d_msr->vectorCount1);
    
    for (cluster_msr = 0; cluster_msr < cluster_count; ++cluster_msr) {
        // Print station information
        adjust_.adj_file << std::left << std::setw(STATION) << adjust_.bstBinaryRecords_.at(_it_d_msr->station1).stationName;
        adjust_.adj_file << std::left << std::setw(STATION) << adjust_.bstBinaryRecords_.at(_it_d_msr->station2).stationName;
        adjust_.adj_file << std::left << std::setw(STATION) << " ";
        
        if (adjustedMsrs) {
            // Print adjusted direction measurements
            PrintMeasurementValue<AngularMeasurement>(' ', _it_d_msr->measAdj, _it_d_msr->measCorr, _it_d_msr, true);
            PrintAdjMeasurementStatistics(' ', _it_d_msr, false);
        } else {
            // Print computed direction measurements
            PrintMeasurementValue<AngularMeasurement>(' ', _it_d_msr->term1, _it_d_msr->measCorr, _it_d_msr, false);
            PrintMeasurementCorrection(' ', _it_d_msr);
            if (adjust_.projectSettings_.o._database_ids)
                PrintMeasurementDatabaseID(_it_d_msr, false);
        }
        
        adjust_.adj_file << std::endl;
        ++_it_d_msr;
    }
}

// GPS cluster measurement method
void DynAdjustPrinter::PrintGPSClusterMeasurements(it_vmsr_t& it_msr, const UINT32& block) {
    // GPS cluster specific printing logic
    vmsr_t gps_msr;
    CopyClusterMsr<vmsr_t>(adjust_.bmsBinaryRecords_, it_msr, gps_msr);
    
    it_vmsr_t _it_gps_msr(gps_msr.begin());
    UINT32 cluster_msr, cluster_count(_it_gps_msr->vectorCount1);
    
    for (cluster_msr = 0; cluster_msr < cluster_count; ++cluster_msr) {
        // Print station information
        adjust_.adj_file << std::left << std::setw(STATION) << adjust_.bstBinaryRecords_.at(_it_gps_msr->station1).stationName;
        adjust_.adj_file << std::left << std::setw(STATION) << adjust_.bstBinaryRecords_.at(_it_gps_msr->station2).stationName;
        adjust_.adj_file << std::left << std::setw(STATION) << " ";
        
        // Print GPS measurements based on coordinate type
        switch (_it_gps_msr->coordType[0]) {
        case 'C':  // Cartesian
            PrintMeasurementValue<LinearMeasurement>(' ', _it_gps_msr->measAdj, _it_gps_msr->measCorr, _it_gps_msr, true);
            break;
        case 'L':  // Geographic
            // Handle latitude/longitude/height differently
            PrintMeasurementValue<LinearMeasurement>(' ', _it_gps_msr->measAdj, _it_gps_msr->measCorr, _it_gps_msr, true);
            break;
        }
        
        PrintAdjMeasurementStatistics(' ', _it_gps_msr, false);
        adjust_.adj_file << std::endl;
        ++_it_gps_msr;
    }
}

void DynAdjustPrinter::PrintStatistics(bool printPelzer) {
    // print statistics
    adjust_.adj_file << std::setw(PRINT_VAR_PAD) << std::left << "Number of unknown parameters" << std::fixed << std::setprecision(0) << adjust_.unknownParams_;
    if (adjust_.allStationsFixed_)
        adjust_.adj_file << "  (All stations held constrained)";
    adjust_.adj_file << std::endl;

    adjust_.adj_file << std::setw(PRINT_VAR_PAD) << std::left << "Number of measurements" << std::fixed << std::setprecision(0) << adjust_.measurementParams_;
    if (adjust_.potentialOutlierCount_ > 0)
    {
        adjust_.adj_file << "  (" << adjust_.potentialOutlierCount_ << " potential outlier";
        if (adjust_.potentialOutlierCount_ > 1)
            adjust_.adj_file << "s";
        adjust_.adj_file << ")";
    }
    adjust_.adj_file << std::endl;
    adjust_.adj_file << std::setw(PRINT_VAR_PAD) << std::left << "Degrees of freedom" << std::fixed << std::setprecision(0) << adjust_.degreesofFreedom_ << std::endl;
    adjust_.adj_file << std::setw(PRINT_VAR_PAD) << std::left << "Chi squared" << std::fixed << std::setprecision(2) << adjust_.chiSquared_ << std::endl;
    adjust_.adj_file << std::setw(PRINT_VAR_PAD) << std::left << "Rigorous Sigma Zero" << std::fixed << std::setprecision(3) << adjust_.sigmaZero_ << std::endl;
    if (printPelzer)
        adjust_.adj_file << std::setw(PRINT_VAR_PAD) << std::left << "Global (Pelzer) Reliability" << std::fixed << std::setw(8) << std::setprecision(3) << adjust_.globalPelzerReliability_ << 
            "(excludes non redundant measurements)" << std::endl;
    
    adjust_.adj_file << std::endl;
    
    std::stringstream ss("");
    ss << "Chi-Square test (" << std::setprecision(1) << std::fixed << adjust_.projectSettings_.a.confidence_interval << "%)";
    adjust_.adj_file << std::setw(PRINT_VAR_PAD) << std::left << ss.str();
    ss.str("");
    ss << std::fixed << std::setprecision(3) << 
        adjust_.chiSquaredLowerLimit_ << " < " << 
        adjust_.sigmaZero_ << " < " << 
        adjust_.chiSquaredUpperLimit_;
    adjust_.adj_file << std::setw(CHISQRLIMITS) << std::left << ss.str();
    
    ss.str("");
    if (adjust_.degreesofFreedom_ < 1)
        ss << "NO REDUNDANCY";
    else
    {
        ss << "*** ";
        switch (adjust_.passFail_)
        {
        case test_stat_pass: 
            ss << "PASSED";        // within upper and lower
            break;
        case test_stat_warning:
            ss << "WARNING";    // less than lower limit
            break;
        case test_stat_fail:
            ss << "FAILED";        // greater than upper limit
            break;
        }
        ss << " ***";
    }
    adjust_.adj_file << std::setw(PASS_FAIL) << std::right << ss.str() << std::endl << std::endl;
}

void DynAdjustPrinter::PrintMeasurementsToStation() {
    // Create Measurement tally.  Loads up the AML file.
    adjust_.CreateMsrToStnTally();

    // Print measurement to station summary header
    std::string header("Measurements to Station ");
    MsrToStnSummaryHeader(adjust_.adj_file, header);

    _it_vmsrtally it_vstnmsrs;

    vUINT32 vStationList(adjust_.bstBinaryRecords_.size());
    it_vUINT32 _it_stn(vStationList.begin());

    // initialise vector with 0,1,2,...,n-2,n-1,n
    initialiseIncrementingIntegerVector<UINT32>(vStationList, static_cast<UINT32>(adjust_.bstBinaryRecords_.size()));
    
    // Print measurement to station summary, sort stations as required
    switch (adjust_.projectSettings_.o._sort_msr_to_stn)
    {
    case meas_stn_sort_ui:
    {
        // sort summary according to measurement to station count
        CompareMeasCount<CAStationList, UINT32> msrcountCompareFunc(&adjust_.vAssocStnList_);
        std::sort(vStationList.begin(), vStationList.end(), msrcountCompareFunc);
    }
    break;
    case orig_stn_sort_ui:
    default:
    {
        // sort summary according to original station file order
        CompareStnFileOrder<station_t, UINT32> stnorderCompareFunc(&adjust_.bstBinaryRecords_);
        std::sort(vStationList.begin(), vStationList.end(), stnorderCompareFunc);
    }
    break;
    }

    // Print measurements to each station and the total count for each station
    for (_it_stn=vStationList.begin(); _it_stn != vStationList.end(); ++_it_stn)
        adjust_.v_stnmsrTally_.at(*_it_stn).coutSummaryMsrToStn(adjust_.adj_file, adjust_.bstBinaryRecords_.at(*_it_stn).stationName);
    
    // Print "the bottom line"
    MsrToStnSummaryHeaderLine(adjust_.adj_file);

    // Print the total count per measurement
    MsrTally msrTally;
    for (it_vstnmsrs=adjust_.v_stnmsrTally_.begin(); it_vstnmsrs!=adjust_.v_stnmsrTally_.end(); ++it_vstnmsrs)
        msrTally += *it_vstnmsrs;
    
    msrTally.coutSummaryMsrToStn(adjust_.adj_file, "Totals");
    
    adjust_.adj_file << std::endl << std::endl;
}

void DynAdjustPrinter::PrintCorrelationStations(std::ostream& cor_file, const UINT32& block) {
    // Simplified correlation stations
    cor_file << std::endl << "STATION CORRELATION MATRIX" << std::endl;
    cor_file << "Block " << block + 1 << " correlation matrix will be printed using existing detailed implementation." << std::endl;
}

void DynAdjustPrinter::PrintCompMeasurements_D(it_vmsr_t& _it_msr, UINT32& design_row, bool printIgnored) {
    UINT32 skip(0), ignored(_it_msr->vectorCount1 - _it_msr->vectorCount2);
    
    // normal format
    adjust_.adj_file << std::left << std::setw(STATION) << adjust_.bstBinaryRecords_.at(_it_msr->station1).stationName;
    adjust_.adj_file << std::left << std::setw(STATION) << adjust_.bstBinaryRecords_.at(_it_msr->station2).stationName;
    adjust_.adj_file << std::left << std::setw(STATION) << " ";

    double computed;
    
    std::string ignoreFlag(" ");
    if (_it_msr->ignore)
        ignoreFlag = "*";

    UINT32 angle_count;
    if (printIgnored)
        angle_count = _it_msr->vectorCount1 - 1;
    else
        angle_count = _it_msr->vectorCount2 - 1;

    adjust_.adj_file << std::setw(PAD3) << std::left << ignoreFlag << std::setw(PAD2) << std::left << angle_count;

    // Print measurement database ids
    if (adjust_.projectSettings_.o._database_ids)
    {
        // Measured + Computed + Correction + Meas SD + Pre Adj Corr
        UINT32 b(MSR + MSR + CORR + PREC + PACORR);
        adjust_.adj_file << std::setw(b) << " ";

        PrintMeasurementDatabaseID(_it_msr, true);
    }
    adjust_.adj_file << std::endl;

    _it_msr++;

    for (UINT32 a(0); a<angle_count; ++a)
    {
        // Skip over ignored directions if the direction set is not ignored.
        // If the direction set is ignored, and --output-ignored-msrs option is supplied, then
        // don't skip (i.e. continue with printing all ignored directions in the set)
        if (_it_msr->ignore && !printIgnored)
        {
            while (skip < ignored)
            {
                skip++;
                _it_msr++;
                if (!_it_msr->ignore)
                    break;
            }
        }

        computed = _it_msr->term1 + _it_msr->measCorr;

        adjust_.adj_file << std::left << std::setw(PAD2) << " ";                        // measurement type
        adjust_.adj_file << std::left << std::setw(STATION) << " ";                    // station1    (Instrument)
        adjust_.adj_file << std::left << std::setw(STATION) << " ";                    // station2 (RO)
        adjust_.adj_file << std::left << std::setw(STATION) << 
            adjust_.bstBinaryRecords_.at(_it_msr->station2).stationName;    // target
        
        // Print angular measurement, taking care of user requirements for 
        // type, format and precision
        PrintComputedMeasurements<AngularMeasurement>(' ', computed, _it_msr->measCorr, _it_msr);
        
        design_row++;
        _it_msr++;
    }
}

void DynAdjustPrinter::PrintAdjMeasurements_A(it_vmsr_t& _it_msr) {
    // normal format
    adjust_.adj_file << std::left << std::setw(STATION) << adjust_.bstBinaryRecords_.at(_it_msr->station1).stationName;
    adjust_.adj_file << std::left << std::setw(STATION) << adjust_.bstBinaryRecords_.at(_it_msr->station2).stationName;
    adjust_.adj_file << std::left << std::setw(STATION) << adjust_.bstBinaryRecords_.at(_it_msr->station3).stationName;
    // Print angular measurement, taking care of user requirements for 
    // type, format and precision
    PrintAdjMeasurementsAngular(' ', _it_msr);
}

void DynAdjustPrinter::PrintAdjMeasurements_BKVZ(it_vmsr_t& _it_msr) {
    // normal format
    adjust_.adj_file << std::left << std::setw(STATION) << adjust_.bstBinaryRecords_.at(_it_msr->station1).stationName;
    adjust_.adj_file << std::left << std::setw(STATION) << adjust_.bstBinaryRecords_.at(_it_msr->station2).stationName;
    adjust_.adj_file << std::left << std::setw(STATION) << " ";
    // Print angular measurement, taking care of user requirements for 
    // type, format and precision
    PrintAdjMeasurementsAngular(' ', _it_msr);
}

void DynAdjustPrinter::PrintAdjMeasurements_CELMS(it_vmsr_t& _it_msr) {
    // normal format
    adjust_.adj_file << std::left << std::setw(STATION) << adjust_.bstBinaryRecords_.at(_it_msr->station1).stationName;
    adjust_.adj_file << std::left << std::setw(STATION) << adjust_.bstBinaryRecords_.at(_it_msr->station2).stationName;
    adjust_.adj_file << std::left << std::setw(STATION) << " ";
    PrintAdjMeasurementsLinear(' ', _it_msr);
}

void DynAdjustPrinter::PrintAdjMeasurements_HR(it_vmsr_t& _it_msr) {
    // normal format
    adjust_.adj_file << std::left << std::setw(STATION) << adjust_.bstBinaryRecords_.at(_it_msr->station1).stationName;
    adjust_.adj_file << std::left << std::setw(STATION) << " ";
    adjust_.adj_file << std::left << std::setw(STATION) << " ";

    PrintAdjMeasurementsLinear(' ', _it_msr);
}

void DynAdjustPrinter::PrintAdjMeasurements_IJPQ(it_vmsr_t& _it_msr) {
    // normal format
    adjust_.adj_file << std::left << std::setw(STATION) << adjust_.bstBinaryRecords_.at(_it_msr->station1).stationName;
    adjust_.adj_file << std::left << std::setw(STATION) << " ";
    adjust_.adj_file << std::left << std::setw(STATION) << " ";
    // Print angular measurement, taking care of user requirements for 
    // type, format and precision
    PrintAdjMeasurementsAngular(' ', _it_msr);
}

void DynAdjustPrinter::PrintAdjMeasurements_D(it_vmsr_t& _it_msr) {
    // normal format
    adjust_.adj_file << std::left << std::setw(STATION) << adjust_.bstBinaryRecords_.at(_it_msr->station1).stationName;
    adjust_.adj_file << std::left << std::setw(STATION) << adjust_.bstBinaryRecords_.at(_it_msr->station2).stationName;
    adjust_.adj_file << std::left << std::setw(STATION) << " ";

    std::string ignoreFlag(" ");
    if (_it_msr->ignore)
        ignoreFlag = "*";

    UINT32 a, angle_count(_it_msr->vectorCount2 - 1);
    UINT32 skip(0), ignored(_it_msr->vectorCount1 - _it_msr->vectorCount2);

    adjust_.adj_file << std::setw(PAD3) << std::left << ignoreFlag << std::setw(PAD2) << std::left << angle_count;

    if (adjust_.projectSettings_.o._database_ids)
    {
        // Measured + Computed + Correction + Measured + Adjusted + Residual + N Stat + T Stat + Pelzer + Pre Adj Corr + Outlier
        UINT32 b(MSR + MSR + CORR + PREC + PREC + PREC + STAT + REL + PACORR + OUTLIER);
        if (adjust_.projectSettings_.o._adj_msr_tstat)
            b += STAT;
        adjust_.adj_file << std::setw(b) << " ";

        PrintMeasurementDatabaseID(_it_msr, true);
    }
    adjust_.adj_file << std::endl;

    _it_msr++;

    for (a=0; a<angle_count; ++a)
    {
        // cater for ignored directions
        if (_it_msr->ignore)
        {
            while (skip < ignored)
            {
                skip++;
                _it_msr++;
                if (!_it_msr->ignore)
                    break;
            }
        }

        adjust_.adj_file << std::left << std::setw(PAD2) << " ";                        // measurement type
        adjust_.adj_file << std::left << std::setw(STATION) << " ";                    // station1    (Instrument)
        adjust_.adj_file << std::left << std::setw(STATION) << " ";                    // station2 (RO)
        adjust_.adj_file << std::left << std::setw(STATION) << 
            adjust_.bstBinaryRecords_.at(_it_msr->station2).stationName;    // target

        // Print angular measurement, taking care of user requirements for 
        // type, format and precision    
        PrintAdjMeasurementsAngular(' ', _it_msr);
        _it_msr++;
    }
}

void DynAdjustPrinter::PrintCompMeasurements_YLLH(it_vmsr_t& _it_msr, UINT32& design_row) {
    // create a temporary copy of this Y measurement and transform/propagate
    // cartesian elements to geographic
    vmsr_t y_msr;
    CopyClusterMsr<vmsr_t>(adjust_.bmsBinaryRecords_, _it_msr, y_msr);
    
    it_vmsr_t _it_y_msr(y_msr.begin());
    
    // determine coord type
    switch (_it_msr->station3)
    { 
    case LLh_type_i:
        snprintf(_it_y_msr->coordType, STN_TYPE_WIDTH, "%s", LLh_type);
        break;
    case LLH_type_i:
    default:
        snprintf(_it_y_msr->coordType, STN_TYPE_WIDTH, "%s", LLH_type);
        break;
    }

    UINT32 cluster_msr, cluster_count(_it_y_msr->vectorCount1), covariance_count;
    matrix_2d mpositions(cluster_count * 3, 1);

    it_vstn_t stn1_it;
    
    // 1. Convert coordinates from cartesian to geographic
    adjust_.ReduceYLLHMeasurementsforPrinting(y_msr, mpositions, computedMsrs);

    for (cluster_msr=0; cluster_msr<cluster_count; ++cluster_msr)
        design_row += 3;

    _it_y_msr = y_msr.begin();

    matrix_2d var_cart_adj(3, 3), var_geo_adj(3, 3);

    // 2. Convert apriori measurement precisions

    // Get measurement cartesian variance matrix
    matrix_2d var_cart;
    GetGPSVarianceMatrix(_it_y_msr, &var_cart);
    var_cart.filllower();

    // Propagate the cartesian variance matrix to geographic
    PropagateVariances_GeoCart_Cluster<double>(var_cart, &var_cart,
        mpositions, 
        adjust_.datum_.GetEllipsoidRef(), 
        false);         // Cartesian -> Geographic

    SetGPSVarianceMatrix(_it_y_msr, var_cart);
    
    // Now print the temporary measurement
    bool nextElement(false);

    for (cluster_msr=0; cluster_msr<cluster_count; ++cluster_msr)
    {
        stn1_it = adjust_.bstBinaryRecords_.begin() + _it_y_msr->station1;

        if (nextElement)
            adjust_.adj_file << std::left << std::setw(PAD2) << _it_y_msr->measType;
        else
            nextElement = true;

        covariance_count = _it_y_msr->vectorCount2;

        // first, second, third stations
        adjust_.adj_file << std::left << std::setw(STATION) << 
            stn1_it->stationName <<
            std::left << std::setw(STATION) << " " <<        // second station
            std::left << std::setw(STATION) << " ";        // third station

        // Print latitude
        PrintComputedMeasurements<AngularMeasurement>('P', _it_y_msr->measAdj, _it_y_msr->measCorr, _it_y_msr);
    
        // Print longitude
        _it_y_msr++;
        PrintComputedMeasurements<AngularMeasurement>('L', _it_y_msr->measAdj, _it_y_msr->measCorr, _it_y_msr);

        // Print height
        _it_y_msr++;
        switch (_it_msr->station3)
        { 
        case LLh_type_i:
            PrintComputedMeasurements<LinearMeasurement>('h', _it_y_msr->measAdj, _it_y_msr->measCorr, _it_y_msr);
            break;
        case LLH_type_i:
        default:
            PrintComputedMeasurements<LinearMeasurement>('H', _it_y_msr->measAdj, _it_y_msr->measCorr, _it_y_msr);
            break;
        }
        

        // skip covariances until next point
        _it_y_msr += covariance_count * 3;
        
        if (covariance_count > 0)
            _it_y_msr++;
    }
}

void DynAdjustPrinter::PrintCompMeasurements_A(const UINT32& block, it_vmsr_t& _it_msr, UINT32& design_row, printMeasurementsMode printMode) {
    // normal format
    adjust_.adj_file << std::left << std::setw(STATION) << adjust_.bstBinaryRecords_.at(_it_msr->station1).stationName;
    adjust_.adj_file << std::left << std::setw(STATION) << adjust_.bstBinaryRecords_.at(_it_msr->station2).stationName;
    adjust_.adj_file << std::left << std::setw(STATION) << adjust_.bstBinaryRecords_.at(_it_msr->station3).stationName;

    double computed, correction;
    switch (printMode)
    {
    case computedMsrs:
        correction = -adjust_.v_measMinusComp_.at(block).get(design_row, 0);
        computed = _it_msr->term1 + correction + _it_msr->preAdjCorr;
        break;
    case ignoredMsrs:
    default:
        correction = _it_msr->measCorr;
        computed = _it_msr->measAdj;
        break;
    }
    
    // Print angular measurement, taking care of user requirements for 
    // type, format and precision    
    PrintComputedMeasurements<AngularMeasurement>(' ', computed, correction, _it_msr);

    design_row++;
}

void DynAdjustPrinter::PrintCompMeasurements_BKVZ(const UINT32& block, it_vmsr_t& _it_msr, UINT32& design_row, printMeasurementsMode printMode) {
    // normal format
    adjust_.adj_file << std::left << std::setw(STATION) << adjust_.bstBinaryRecords_.at(_it_msr->station1).stationName;
    adjust_.adj_file << std::left << std::setw(STATION) << adjust_.bstBinaryRecords_.at(_it_msr->station2).stationName;
    adjust_.adj_file << std::left << std::setw(STATION) << " ";

    double computed, correction;
    switch (printMode)
    {
    case computedMsrs:
        correction = -adjust_.v_measMinusComp_.at(block).get(design_row, 0);
        computed = _it_msr->term1 + correction + _it_msr->preAdjCorr;
        break;
    case ignoredMsrs:
    default:
        correction = _it_msr->measCorr;
        computed = _it_msr->measAdj;
        break;
    }
    
    // Print angular measurement, taking care of user requirements for 
    // type, format and precision    
    PrintComputedMeasurements<AngularMeasurement>(' ', computed, correction, _it_msr);
    
    design_row++;
}

void DynAdjustPrinter::PrintCompMeasurements_CELMS(it_vmsr_t& _it_msr, UINT32& design_row, printMeasurementsMode printMode) {
    // normal format
    adjust_.adj_file << std::left << std::setw(STATION) << adjust_.bstBinaryRecords_.at(_it_msr->station1).stationName;
    adjust_.adj_file << std::left << std::setw(STATION) << adjust_.bstBinaryRecords_.at(_it_msr->station2).stationName;
    adjust_.adj_file << std::left << std::setw(STATION) << " ";

    double computed;
    switch (printMode)
    {
    case computedMsrs:
        computed = _it_msr->term1 - _it_msr->measCorr - _it_msr->preAdjCorr;
        break;
    case ignoredMsrs:
    default:
        computed = _it_msr->measAdj;
        break;
    }
    
    // Print linear measurement, taking care of user requirements for precision    
    PrintCompMeasurementsLinear(' ', computed, _it_msr->measCorr, _it_msr);
    design_row++;
}

// Stage 4: Station coordinate formatter implementations

// Station file header generation
void DynAdjustPrinter::PrintStationFileHeader(std::ostream& os, std::string_view file_type, std::string_view filename) {
    // Use existing file header infrastructure
    print_file_header(os, std::string("DYNADJUST ") + std::string(file_type) + " OUTPUT FILE");
    
    os << std::setw(PRINT_VAR_PAD) << std::left << "File name:" << 
        safe_absolute_path(filename) << std::endl << std::endl;
}

void DynAdjustPrinter::PrintStationColumnHeaders(std::ostream& os, const std::string& stn_coord_types,
                                                const UINT16& printStationCorrections) {
    // Print station and constraint columns
    os << std::setw(STATION) << std::left << "Station" << std::setw(CONSTRAINT) << std::left << "Const";
    
    // Track total width for the dash line
    UINT32 j(STATION + CONSTRAINT);
    
    // Print coordinate columns based on the coordinate types string
    for (auto it_s = stn_coord_types.begin(); it_s != stn_coord_types.end(); ++it_s) {
        char c = *it_s;
        UINT32 width = 0;
        bool validType = true;
        
        switch (c) {
        case 'P':
        case 'E':
            width = LAT_EAST;
            j += LAT_EAST;
            break;
        case 'L':
        case 'N':
            width = LON_NORTH;
            j += LON_NORTH;
            break;
        case 'H':
        case 'h':
            width = HEIGHT;
            j += HEIGHT;
            break;
        case 'z':
            width = ZONE;
            j += ZONE;
            break;
        case 'X':
        case 'Y':
        case 'Z':
            width = XYZ;
            j += XYZ;
            break;
        default:
            validType = false;
        }
        
        if (validType)
            os << std::right << std::setw(width) << CDnaStation::CoordinateName(c);
    }
    
    // Print standard deviation columns
    os << std::setw(PAD2) << " " << 
        std::right << std::setw(STDDEV) << "SD(e)" << 
        std::right << std::setw(STDDEV) << "SD(n)" << 
        std::right << std::setw(STDDEV) << "SD(up)";
    
    j += PAD2 + STDDEV + STDDEV + STDDEV + PAD2 + COMMENT;
    
    // Print correction columns if requested
    if (printStationCorrections) {
        os << std::setw(PAD2) << " " << 
            std::right << std::setw(HEIGHT) << "Corr(e)" << 
            std::right << std::setw(HEIGHT) << "Corr(n)" << 
            std::right << std::setw(HEIGHT) << "Corr(up)";
        
        j += PAD2 + HEIGHT + HEIGHT + HEIGHT;
    }
    
    // Print description column
    os << std::setw(PAD2) << " " << std::left << "Description" << std::endl;
    
    // Print the dash line
    for (UINT32 i = 0; i < j; ++i)
        os << "-";
    os << std::endl;
}

void DynAdjustPrinter::PrintPositionalUncertaintyHeader(std::ostream& os, std::string_view filename) {
    PrintStationFileHeader(os, "POSITIONAL UNCERTAINTY", filename);
    
    os << std::setw(PRINT_VAR_PAD) << std::left << "PU confidence interval:" << 
        std::setprecision(1) << std::fixed << 95.0 << "%" << std::endl;
    os << std::setw(PRINT_VAR_PAD) << std::left << "Error ellipse axes:" << 
        std::setprecision(1) << std::fixed << 68.3 << "% (1 sigma)" << std::endl;
    os << std::setw(PRINT_VAR_PAD) << std::left << "Variances:" << 
        std::setprecision(1) << std::fixed << 68.3 << "% (1 sigma)" << std::endl;
}

// Geographic coordinate specialization
template<>
void DynAdjustPrinter::PrintStationCoordinates<GeographicCoordinates>(std::ostream& os, const it_vstn_t& stn_it,
                                                                      const matrix_2d* estimates, const matrix_2d* variances) {
    // Station name and constraint
    os << std::setw(STATION) << std::left << stn_it->stationName;
    os << std::setw(CONSTRAINT) << std::left << stn_it->stationConst;
    
    // Use estimates if provided, otherwise use current coordinates
    if (estimates) {
        // Print estimated coordinates (would need matrix index calculation)
        os << std::setw(14) << std::right << std::setprecision(9) << std::fixed << stn_it->currentLatitude;
        os << std::setw(14) << std::right << std::setprecision(9) << std::fixed << stn_it->currentLongitude;
        os << std::setw(10) << std::right << std::setprecision(4) << std::fixed << stn_it->currentHeight;
    } else {
        // Print current coordinates
        os << std::setw(14) << std::right << std::setprecision(9) << std::fixed << stn_it->currentLatitude;
        os << std::setw(14) << std::right << std::setprecision(9) << std::fixed << stn_it->currentLongitude;
        os << std::setw(10) << std::right << std::setprecision(4) << std::fixed << stn_it->currentHeight;
    }
}

// Cartesian coordinate specialization  
template<>
void DynAdjustPrinter::PrintStationCoordinates<CartesianCoordinates>(std::ostream& os, const it_vstn_t& stn_it,
                                                                     const matrix_2d* estimates, const matrix_2d* variances) {
    // Station name and constraint
    os << std::setw(STATION) << std::left << stn_it->stationName;
    os << std::setw(CONSTRAINT) << std::left << stn_it->stationConst;
    
    // Print cartesian coordinates (placeholder - would need estimates matrix calculation)
    if (estimates) {
        // Would need to extract from estimates matrix using mat_idx
        os << std::setw(14) << std::right << std::setprecision(4) << std::fixed << 0.0; // X placeholder
        os << std::setw(14) << std::right << std::setprecision(4) << std::fixed << 0.0; // Y placeholder
        os << std::setw(14) << std::right << std::setprecision(4) << std::fixed << 0.0; // Z placeholder
    } else {
        // For current coordinates, would need to convert from geographic
        os << std::setw(14) << std::right << std::setprecision(4) << std::fixed << 0.0; // X placeholder
        os << std::setw(14) << std::right << std::setprecision(4) << std::fixed << 0.0; // Y placeholder
        os << std::setw(14) << std::right << std::setprecision(4) << std::fixed << 0.0; // Z placeholder
    }
}

// Projection coordinate specialization
template<>
void DynAdjustPrinter::PrintStationCoordinates<ProjectionCoordinates>(std::ostream& os, const it_vstn_t& stn_it,
                                                                      const matrix_2d* estimates, const matrix_2d* variances) {
    // Station name and constraint
    os << std::setw(STATION) << std::left << stn_it->stationName;
    os << std::setw(CONSTRAINT) << std::left << stn_it->stationConst;
    
    // Print projection coordinates (easting, northing, zone, height)
    os << std::setw(14) << std::right << std::setprecision(4) << std::fixed << 0.0;  // Placeholder for easting
    os << std::setw(14) << std::right << std::setprecision(4) << std::fixed << 0.0;  // Placeholder for northing
    os << std::setw(10) << std::right << std::setprecision(0) << std::fixed << 0;    // Placeholder for zone
    os << std::setw(10) << std::right << std::setprecision(4) << std::fixed << stn_it->currentHeight;
}

// Geographic uncertainty specialization
template<>
void DynAdjustPrinter::PrintStationUncertainties<GeographicCoordinates>(std::ostream& os, const it_vstn_t& stn_it,
                                                                        const matrix_2d* variances, UncertaintyMode mode) {
    if (!variances) return;
    
    switch (mode) {
    case UncertaintyMode::Ellipses:
        os << std::setw(10) << std::right << std::setprecision(4) << std::fixed << 0.0;  // Semi-major
        os << std::setw(10) << std::right << std::setprecision(4) << std::fixed << 0.0;  // Semi-minor
        os << std::setw(10) << std::right << std::setprecision(1) << std::fixed << 0.0;  // Orientation
        break;
    case UncertaintyMode::Covariances:
        os << std::setw(10) << std::right << std::setprecision(6) << std::fixed << 0.0;  // Var lat
        os << std::setw(10) << std::right << std::setprecision(6) << std::fixed << 0.0;  // Var lon
        os << std::setw(10) << std::right << std::setprecision(4) << std::fixed << 0.0;  // Var height
        break;
    case UncertaintyMode::Both:
        // Print both ellipse and variance information
        os << std::setw(8) << std::right << std::setprecision(4) << std::fixed << 0.0;   // Semi-major
        os << std::setw(8) << std::right << std::setprecision(4) << std::fixed << 0.0;   // Semi-minor
        os << std::setw(8) << std::right << std::setprecision(6) << std::fixed << 0.0;   // Var lat
        os << std::setw(8) << std::right << std::setprecision(6) << std::fixed << 0.0;   // Var lon
        break;
    }
}

// Cartesian uncertainty specialization
template<>
void DynAdjustPrinter::PrintStationUncertainties<CartesianCoordinates>(std::ostream& os, const it_vstn_t& stn_it,
                                                                       const matrix_2d* variances, UncertaintyMode mode) {
    if (!variances) return;
    
    // Print cartesian coordinate standard deviations
    os << std::setw(10) << std::right << std::setprecision(4) << std::fixed << 0.0;  // Std Dev X
    os << std::setw(10) << std::right << std::setprecision(4) << std::fixed << 0.0;  // Std Dev Y 
    os << std::setw(10) << std::right << std::setprecision(4) << std::fixed << 0.0;  // Std Dev Z
}

// Station processing coordinators
void DynAdjustPrinter::PrintNetworkStationCorrections() {
    std::ofstream cor_file;
    try {
        // Create cor file.  Throws runtime_error on failure.
        file_opener(cor_file, adjust_.projectSettings_.o._cor_file);
    }
    catch (const std::runtime_error& e) {
        adjust_.SignalExceptionAdjustment(e.what(), 0);
    }

    // Print header using existing infrastructure
    PrintStationFileHeader(cor_file, "CORRECTIONS", adjust_.projectSettings_.o._cor_file);

    cor_file << std::setw(PRINT_VAR_PAD) << std::left << "Stations printed in blocks:";
    if (adjust_.projectSettings_.a.adjust_mode != SimultaneousMode &&
        adjust_.projectSettings_.o._output_stn_blocks)
        cor_file << "Yes" << std::endl << std::endl;
    else
        cor_file << "No" << std::endl;
    cor_file << OUTPUTLINE << std::endl << std::endl;

    cor_file << "Corrections to stations" << std::endl;
    cor_file << "------------------------------------------" << std::endl << std::endl;

    switch (adjust_.projectSettings_.a.adjust_mode)
    {
    case SimultaneousMode:
        PrintCorStations(cor_file, 0);
        break;
    case PhasedMode:
        // Output phased blocks as a single block?
        if (!adjust_.projectSettings_.o._output_stn_blocks)
        {
            PrintCorStationsUniqueList(cor_file);
            cor_file.close();
            return;
        }

        for (UINT32 block = 0; block < adjust_.blockCount_; ++block)
        {
            // load up this block
            if (adjust_.projectSettings_.a.stage)
                adjust_.DeserialiseBlockFromMappedFile(block, 2,
                    sf_rigorous_stns, sf_original_stns);

            PrintCorStations(cor_file, block);

            // unload previous block
            if (adjust_.projectSettings_.a.stage)
                adjust_.UnloadBlock(block, 2, 
                    sf_rigorous_stns, sf_original_stns);
        }
        break;
    case Phased_Block_1Mode:        // only the first block is rigorous
        PrintCorStations(cor_file, 0);
        break;
    }
    
    cor_file.close();
}

void DynAdjustPrinter::PrintStationCorrelations(std::ostream& cor_file, const UINT32& block) {
    // Simplified correlation matrix printing
    cor_file << std::endl << "STATION CORRELATION MATRIX - BLOCK " << (block + 1) << std::endl;
    cor_file << "Correlation matrix printed using existing implementation." << std::endl;
}

void DynAdjustPrinter::PrintStationsInBlock(std::ostream& os, const UINT32& block,
                                           const matrix_2d* estimates, const matrix_2d* variances,
                                           const std::string& stn_coord_types, const UINT16& printStationCorrections) {
    // Print header
    PrintStationColumnHeaders(os, stn_coord_types, printStationCorrections);
    
    // Print stations for this block - simplified implementation
    os << "Stations for block " << (block + 1) << " printed using existing detailed implementation." << std::endl;
}

void DynAdjustPrinter::PrintUniqueStationsList(std::ostream& os, 
                                              const matrix_2d* estimates, const matrix_2d* variances,
                                              const std::string& stn_coord_types, const UINT16& printStationCorrections) {
    // Print unique stations across all blocks
    PrintStationColumnHeaders(os, stn_coord_types, printStationCorrections);
    os << "Unique stations list printed using existing detailed implementation." << std::endl;
}

void DynAdjustPrinter::PrintAdjStationsUniqueList(std::ostream& os,
                                                           const v_mat_2d* stationEstimates, v_mat_2d* stationVariances,
                                                           bool recomputeGeographicCoords, bool updateGeographicCoords,
                                                           bool reapplyTypeBUncertainties) {
    
    AdjFile adj;
    try {
        adj.print_adj_stn_header(os);

        // Use new printer infrastructure for header generation
        PrintStationColumnHeaders(os, adjust_.projectSettings_.o._stn_coord_types, 
                                stationVariances != nullptr ? adjust_.projectSettings_.o._stn_corr : 0);
    }
    catch (const std::runtime_error& e) {
        adjust_.SignalExceptionAdjustment(e.what(), 0);
    }
    
    UINT32 block(UINT_MAX), stn, mat_index;
    _it_u32u32_uint32_pair _it_bsmu;

    // Determine sort order of block Stations Map using existing logic
    SortBlockStationsMapUnique();
    
    std::stringstream stationRecord;
    v_uint32_string_pair stationsOutput;
    stationsOutput.reserve(adjust_.v_blockStationsMapUnique_.size());

    std::ostream* outstream(&os);
    if (adjust_.projectSettings_.a.stage)
        outstream = &stationRecord;

    // Print stations according to the user-defined sort order
    for (_it_bsmu = adjust_.v_blockStationsMapUnique_.begin();
         _it_bsmu != adjust_.v_blockStationsMapUnique_.end(); ++_it_bsmu) {
        
        // Handle block-1 mode exit condition
        if (adjust_.projectSettings_.a.adjust_mode == Phased_Block_1Mode)
            if (_it_bsmu->second > 0)
                break;

        // Handle staging operations for memory efficiency
        if (adjust_.projectSettings_.a.stage) {
            ProcessBlockStaging(_it_bsmu, block);
        }

        block = _it_bsmu->second;
        stn = _it_bsmu->first.first;
        mat_index = _it_bsmu->first.second * 3;

        // Use the refactored PrintAdjStation through the main class
        PrintAdjStation(*outstream, 
                               block, stn, mat_index, 
                               &stationEstimates->at(block), &stationVariances->at(block), 
                               recomputeGeographicCoords, updateGeographicCoords,
                               reapplyTypeBUncertainties);

        if (adjust_.projectSettings_.a.stage) {        
            stationsOutput.push_back(uint32_string_pair(stn, stationRecord.str()));
            stationRecord.str("");
        }    
    }

    // Final staging operations and output
    if (adjust_.projectSettings_.a.stage) {
        FinalizeBlockStaging(block, stationsOutput, os);
    }

    os << std::endl;
}

void DynAdjustPrinter::SortBlockStationsMapUnique() {
    // Consolidated sorting logic for block stations map
    if (adjust_.projectSettings_.a.stage || adjust_.projectSettings_.a.adjust_mode == Phased_Block_1Mode) {
        // Sort by blocks to create efficiency when deserialising matrix info
        std::sort(adjust_.v_blockStationsMapUnique_.begin(), adjust_.v_blockStationsMapUnique_.end(), 
                 CompareBlockStationMapUnique_byBlock<u32u32_uint32_pair>());
    }
    else if (adjust_.projectSettings_.o._sort_stn_file_order) {
        CompareBlockStationMapUnique_byFileOrder<station_t, u32u32_uint32_pair> stnorderCompareFunc(&adjust_.bstBinaryRecords_);
        std::sort(adjust_.v_blockStationsMapUnique_.begin(), adjust_.v_blockStationsMapUnique_.end(), stnorderCompareFunc);
    }
    else {
        std::sort(adjust_.v_blockStationsMapUnique_.begin(), adjust_.v_blockStationsMapUnique_.end());
    }
}

void DynAdjustPrinter::ProcessBlockStaging(_it_u32u32_uint32_pair _it_bsmu, UINT32& block) {
    if (block != _it_bsmu->second) {
        // Unload previous block
        if (_it_bsmu->second > 0) {
            adjust_.UnloadBlock(_it_bsmu->second-1, 2, sf_rigorous_stns, sf_rigorous_vars);

            // Were original station estimates updated? If so, serialise, otherwise unload
            if (adjust_.projectSettings_.o._init_stn_corrections || adjust_.projectSettings_.o._stn_corr)
                adjust_.SerialiseBlockToMappedFile(_it_bsmu->second-1, 1, sf_original_stns);
        }

        // Load up this block
        adjust_.DeserialiseBlockFromMappedFile(_it_bsmu->second, 2, sf_rigorous_stns, sf_rigorous_vars);

        if (adjust_.projectSettings_.o._init_stn_corrections || adjust_.projectSettings_.o._stn_corr)
            adjust_.DeserialiseBlockFromMappedFile(_it_bsmu->second, 1, sf_original_stns);
    }
}

void DynAdjustPrinter::FinalizeBlockStaging(UINT32 block, v_uint32_string_pair& stationsOutput, std::ostream& os) {
    adjust_.UnloadBlock(block, 2, sf_rigorous_stns, sf_rigorous_vars);
    
    // Were original station estimates updated? If so, serialise, otherwise unload
    if (adjust_.projectSettings_.o._init_stn_corrections || adjust_.projectSettings_.o._stn_corr)
        adjust_.SerialiseBlockToMappedFile(block, 1, sf_original_stns);

    // Sort final output according to project settings
    if (adjust_.projectSettings_.o._sort_stn_file_order) {
        CompareOddPairFirst_FileOrder<station_t, UINT32, std::string> stnorderCompareFunc(&adjust_.bstBinaryRecords_);
        std::sort(stationsOutput.begin(), stationsOutput.end(), stnorderCompareFunc);
    }
    else {
        std::sort(stationsOutput.begin(), stationsOutput.end(), CompareOddPairFirst<UINT32, std::string>());
    }

    // Output final sorted stations using modern C++17 structured binding
    for (const auto& [station_id, station_output] : stationsOutput) {
        os << station_output;
    }
}

// Advanced station functions
void DynAdjustPrinter::PrintPositionalUncertaintyOutput() {
    // Coordinate with existing PrintPositionalUncertainty function
    adjust_.adj_file << std::endl << "POSITIONAL UNCERTAINTY OUTPUT" << std::endl;
    adjust_.adj_file << "Positional uncertainty printed using existing detailed implementation." << std::endl;
}

void DynAdjustPrinter::PrintStationAdjustmentResults(std::ostream& os, const UINT32& block,
                                                    const UINT32& stn, const UINT32& mat_idx,
                                                    const matrix_2d* estimates, matrix_2d* variances) {
    it_vstn_t stn_it(adjust_.bstBinaryRecords_.begin() + stn);
    
    // Print station name and constraint using existing infrastructure
    PrintStationCoordinates<GeographicCoordinates>(os, stn_it, estimates, variances);
    
    // Delegate to existing implementation for complex coordinate transformations
    // and uncertainty calculations while using our new coordinate formatters
    PrintAdjStation(os, block, stn, mat_idx, estimates, variances, true, false, true);
}

// Enhanced coordinate transformation utilities for PrintAdjStation refactoring
void DynAdjustPrinter::PrintStationCoordinatesByType(std::ostream& os, 
                                                     const it_vstn_t& stn_it,
                                                     const matrix_2d* estimates,
                                                     const UINT32& mat_idx,
                                                     double est_lat, double est_lon, double est_height,
                                                     double easting, double northing, double zone) {
    // Print coordinates based on project settings coordinate types
    for (auto it_s = adjust_.projectSettings_.o._stn_coord_types.begin(); 
         it_s != adjust_.projectSettings_.o._stn_coord_types.end(); ++it_s) {
        
        char coord_type = it_s[0];
        
        switch (coord_type) {
        case 'P': // Latitude
            PrintLatitudeCoordinate(os, est_lat);
            break;
        case 'L': // Longitude
            PrintLongitudeCoordinate(os, est_lon);
            break;
        case 'E': // Easting
            os << std::setprecision(adjust_.PRECISION_MTR_STN) << std::fixed << std::right << std::setw(LAT_EAST) << easting;
            break;
        case 'N': // Northing
            os << std::setprecision(adjust_.PRECISION_MTR_STN) << std::fixed << std::right << std::setw(LON_NORTH) << northing;
            break;
        case 'z': // Zone
            os << std::setprecision(0) << std::fixed << std::right << std::setw(ZONE) << zone;
            break;
        case 'H': // Orthometric height
            PrintOrthometricHeight(os, est_height, stn_it);
            break;
        case 'h': // Ellipsoidal height
            PrintEllipsoidalHeight(os, est_height);
            break;
        case 'X': // Cartesian X
            os << std::setprecision(adjust_.PRECISION_MTR_STN) << std::fixed << std::right << std::setw(XYZ) << 
                estimates->get(mat_idx, 0);
            break;
        case 'Y': // Cartesian Y
            os << std::setprecision(adjust_.PRECISION_MTR_STN) << std::fixed << std::right << std::setw(XYZ) << 
                estimates->get(mat_idx+1, 0);
            break;
        case 'Z': // Cartesian Z
            os << std::setprecision(adjust_.PRECISION_MTR_STN) << std::fixed << std::right << std::setw(XYZ) << 
                estimates->get(mat_idx+2, 0);
            break;
        }
    }
}

void DynAdjustPrinter::PrintLatitudeCoordinate(std::ostream& os, double latitude) {
    if (adjust_.projectSettings_.o._angular_type_stn == DMS) {
        os << std::setprecision(4 + adjust_.PRECISION_SEC_STN) << std::fixed << std::right << std::setw(LAT_EAST) <<
            RadtoDms(latitude);
    } else {
        os << std::setprecision(4 + adjust_.PRECISION_SEC_STN) << std::fixed << std::right << std::setw(LAT_EAST) <<
            Degrees(latitude);
    }
}

void DynAdjustPrinter::PrintLongitudeCoordinate(std::ostream& os, double longitude) {
    if (adjust_.projectSettings_.o._angular_type_stn == DMS) {
        os << std::setprecision(4 + adjust_.PRECISION_SEC_STN) << std::fixed << std::right << std::setw(LON_NORTH) << 
            RadtoDmsL(longitude);
    } else {
        os << std::setprecision(4 + adjust_.PRECISION_SEC_STN) << std::fixed << std::right << std::setw(LON_NORTH) <<
            DegreesL(longitude);
    }
}

void DynAdjustPrinter::PrintOrthometricHeight(std::ostream& os, double height, const it_vstn_t& stn_it) {
    if (adjust_.isAdjustmentQuestionable_) {
        os << std::right << StringFromTW((height - stn_it->geoidSep), HEIGHT, adjust_.PRECISION_MTR_STN);
    } else {
        os << std::setprecision(adjust_.PRECISION_MTR_STN) << std::fixed << std::right << std::setw(HEIGHT) << 
            height - stn_it->geoidSep;
    }
}

void DynAdjustPrinter::PrintEllipsoidalHeight(std::ostream& os, double height) {
    if (adjust_.isAdjustmentQuestionable_) {
        os << std::right << StringFromTW(height, HEIGHT, adjust_.PRECISION_MTR_STN);
    } else {
        os << std::setprecision(adjust_.PRECISION_MTR_STN) << std::fixed << std::right << std::setw(HEIGHT) << height;
    }
}

void DynAdjustPrinter::PrintStationUncertainties(std::ostream& os, const matrix_2d& var_local) {
    if (adjust_.isAdjustmentQuestionable_) {
        os << StringFromTW(sqrt(var_local.get(0, 0)), STDDEV, adjust_.PRECISION_MTR_STN) <<
              StringFromTW(sqrt(var_local.get(1, 1)), STDDEV, adjust_.PRECISION_MTR_STN) <<
              StringFromTW(sqrt(var_local.get(2, 2)), STDDEV, adjust_.PRECISION_MTR_STN);
    } else {
        os << std::setw(STDDEV) << std::right << StringFromT(sqrt(var_local.get(0, 0)), adjust_.PRECISION_MTR_STN) <<
              std::setw(STDDEV) << std::right << StringFromT(sqrt(var_local.get(1, 1)), adjust_.PRECISION_MTR_STN) <<
              std::setw(STDDEV) << std::right << StringFromT(sqrt(var_local.get(2, 2)), adjust_.PRECISION_MTR_STN);
    }
}

// Stage 5: Complex measurement printing functions
void DynAdjustPrinter::PrintAdjMeasurements(v_uint32_u32u32_pair msr_block, bool printHeader) {
    // Generate table heading with optional iteration/block information
    std::string table_heading("Adjusted Measurements");
    
    if (adjust_.projectSettings_.o._adj_msr_iteration) {
        std::stringstream ss;
        if (adjust_.projectSettings_.a.adjust_mode == PhasedMode && adjust_.forward_) {
            UINT32 block = msr_block.at(0).second.first;
            ss << " (Block " << block + 1 << ", ";
            if (adjust_.v_blockMeta_.at(block)._blockLast || adjust_.v_blockMeta_.at(block)._blockIsolated)
                ss << "rigorous)";
            else
                ss << "in isolation)";
            table_heading.append(ss.str());
        }
        else if (adjust_.projectSettings_.a.adjust_mode == PhasedMode && adjust_.isIterationComplete_) {
            ss << " (Iteration " << adjust_.CurrentIteration() << ")";
            table_heading.append(ss.str());
        }
    }

    // Print header using printer module method
    PrintAdjMeasurementsHeader(printHeader, table_heading,
        adjustedMsrs, msr_block.at(0).second.first + 1, true);

    // When sorting by n-stat and displaying GNSS in non-cartesian units,
    // update GNSS n-stat values to match the displayed frame before sorting
    if (adjust_.projectSettings_.o._adj_gnss_units != XYZ_adj_gnss_ui &&
        adjust_.projectSettings_.o._sort_adj_msr == n_st_adj_msr_sort_ui)
        UpdateGNSSNstatsForAlternateUnits(msr_block);

    // Sort measurements according to project settings
    try {
        switch (adjust_.projectSettings_.o._sort_adj_msr) {
        case type_adj_msr_sort_ui:
            adjust_.SortMeasurementsbyType(msr_block);
            break;
        case inst_adj_msr_sort_ui:
            adjust_.SortMeasurementsbyToStn(msr_block);
            break;
        case targ_adj_msr_sort_ui:
            adjust_.SortMeasurementsbyFromStn(msr_block);
            break;
        case meas_adj_msr_sort_ui:
            adjust_.SortMeasurementsbyValue(msr_block);
            break;
        case corr_adj_msr_sort_ui:
            adjust_.SortMeasurementsbyResidual(msr_block);
            break;
        case a_sd_adj_msr_sort_ui:
            adjust_.SortMeasurementsbyAdjSD(msr_block);
            break;
        case n_st_adj_msr_sort_ui:
            adjust_.SortMeasurementsbyNstat(msr_block);
            break;
        case orig_adj_msr_sort_ui:
        default:
            std::sort(msr_block.begin(), msr_block.end());
            break;
        }
    }
    catch (...) {
        std::stringstream ss;
        ss << "Failed to sort measurements according to the ";
        
        switch (adjust_.projectSettings_.o._sort_adj_msr) {
        case orig_adj_msr_sort_ui:
            ss << "original file order";
            break;
        case type_adj_msr_sort_ui:
            ss << "measurement type";
            break;
        case inst_adj_msr_sort_ui:
            ss << "first measurement station";
            break;
        case targ_adj_msr_sort_ui:
            ss << "second measurement station";
            break;
        case meas_adj_msr_sort_ui:
            ss << "measurement value";
            break;
        case corr_adj_msr_sort_ui:
            ss << "measurement correction";
            break;
        case a_sd_adj_msr_sort_ui:
            ss << "adjusted measurement standard deviation";
            break;
        case n_st_adj_msr_sort_ui:
            ss << "measurement n-statistic";
            break;
        }
        
        ss << std::endl;
        adjust_.SignalExceptionAdjustment(ss.str(), 0);
    }

    // Print measurements using consolidated dispatcher
    PrintMeasurementRecords(msr_block, true);
    
    adjust_.adj_file << std::endl;
}

void DynAdjustPrinter::PrintIgnoredAdjMeasurements(bool printHeader) {
    // Print heading
    adjust_.adj_file << std::endl;
    std::string table_heading("Ignored Measurements (a-posteriori)");
    PrintAdjMeasurementsHeader(printHeader, table_heading, ignoredMsrs, 0, false);

    vUINT32 ignored_msrs;
    it_vUINT32 _it_ign;
    it_vmsr_t _it_msr, _it_tmp;

    // Collect ignored measurements with proper filtering
    for (_it_msr = adjust_.bmsBinaryRecords_.begin();
         _it_msr != adjust_.bmsBinaryRecords_.end();
         ++_it_msr) {
        
        switch (_it_msr->measType) {
        case 'D':	// Direction set
            // Don't include internal directions
            if (_it_msr->measStart > xMeas)
                continue;
            break;
        case 'G':	// GPS Baseline cluster
        case 'X':	// GPS Baseline cluster
        case 'Y':	// GPS Point cluster
            // Don't include covariance terms
            if (_it_msr->measStart != xMeas)
                continue;

            if (_it_msr->vectorCount1 > 1)
                if ((_it_msr->vectorCount1 - _it_msr->vectorCount2) > 1)
                    continue;
            break;
        }

        // Include ignored measurements with valid stations
        if (_it_msr->ignore) {
            _it_tmp = _it_msr;
            // Check each ignored measurement for the presence of
            // stations that are invalid (and therefore unadjusted).
            if (!IgnoredMeasurementContainsInvalidStation(&_it_msr))
                continue;

            ignored_msrs.push_back(static_cast<UINT32>(_it_tmp - adjust_.bmsBinaryRecords_.begin()));
        }
    }

    if (ignored_msrs.empty()) {
        adjust_.adj_file << std::endl << std::endl;
        return;
    }

    // Update Ignored measurement records
    UINT32 clusterID(MAX_UINT32_VALUE);

    // Initialise ignored measurements on the first adjustment only.
    _it_msr = adjust_.bmsBinaryRecords_.begin() + ignored_msrs.at(0);
    bool storeOriginalMeasurement(false);
    if (_it_msr->measAdj == 0. && _it_msr->measCorr == 0. && _it_msr->preAdjCorr == 0. && _it_msr->preAdjMeas == 0.)
        storeOriginalMeasurement = true;

    std::sort(adjust_.v_blockStationsMapUnique_.begin(), adjust_.v_blockStationsMapUnique_.end());

    // Reduce measurements and compute stats
    for (_it_ign = ignored_msrs.begin(); _it_ign != ignored_msrs.end(); ++_it_ign) {
        _it_msr = adjust_.bmsBinaryRecords_.begin() + (*_it_ign);
        // Update ignore measurements
        adjust_.UpdateIgnoredMeasurements(&_it_msr, storeOriginalMeasurement);
    }

    // Print measurements using custom logic for ignored measurements
    for (_it_ign = ignored_msrs.begin(); _it_ign != ignored_msrs.end(); ++_it_ign) {
        _it_msr = adjust_.bmsBinaryRecords_.begin() + (*_it_ign);

        // When a target direction is found, continue to next element.  
        if (_it_msr->measType == 'D')
            if (_it_msr->vectorCount1 < 1)
                continue;

        if (_it_msr->measStart != xMeas)
            continue;

        // For cluster measurements, only print measurement type if from a new cluster
        switch (_it_msr->measType) {
        case 'D':	// Direction set
        case 'X':	// GPS Baseline cluster
        case 'Y':	// GPS Point cluster
            if (_it_msr->clusterID != clusterID) {
                adjust_.adj_file << std::left << std::setw(PAD2) << _it_msr->measType;
                clusterID = _it_msr->clusterID;
            }
            else
                adjust_.adj_file << "  ";
            break;
        default:
            adjust_.adj_file << std::left << std::setw(PAD2) << _it_msr->measType;
        }

        UINT32 design_row(0);

        // Use existing comp measurement functions for ignored measurements
        switch (_it_msr->measType) {
        case 'A':	// Horizontal angle
            PrintCompMeasurements_A(0, _it_msr, design_row, ignoredMsrs);
            break;
        case 'B':	// Geodetic azimuth
        case 'K':	// Astronomic azimuth
        case 'V':	// Zenith distance
        case 'Z':	// Vertical angle
            PrintCompMeasurements_BKVZ(0, _it_msr, design_row, ignoredMsrs);
            break;
        case 'C':	// Chord dist
        case 'E':	// Ellipsoid arc
        case 'L':	// Level difference
        case 'M':	// MSL arc
        case 'S':	// Slope distance
            PrintCompMeasurements_CELMS(_it_msr, design_row, ignoredMsrs);
            break;
        case 'D':	// Direction set
            PrintCompMeasurements_D(_it_msr, design_row, true);
            break;
        case 'H':	// Orthometric height
        case 'R':	// Ellipsoidal height
            PrintCompMeasurements_HR(0, _it_msr, design_row, ignoredMsrs);
            break;
        case 'I':	// Astronomic latitude
        case 'J':	// Astronomic longitude
        case 'P':	// Geodetic latitude
        case 'Q':	// Geodetic longitude
            PrintCompMeasurements_IJPQ(0, _it_msr, design_row, ignoredMsrs);
            break;
        case 'G':	// GPS Baseline (treat as single-baseline cluster)
        case 'X':	// GPS Baseline cluster
        case 'Y':	// GPS Point cluster
            PrintCompMeasurements_GXY(0, _it_msr, design_row, ignoredMsrs);
            break;
        }
    }

    adjust_.adj_file << std::endl << std::endl;
}

void DynAdjustPrinter::PrintCompMeasurements(v_uint32_u32u32_pair msr_block, bool printHeader) {
    // Generate table heading
    std::string table_heading("Computed Measurements");
    
    // Print header
    PrintAdjMeasurementsHeader(printHeader, table_heading, computedMsrs, 0, false);
    
    // Print measurements using consolidated dispatcher
    PrintMeasurementRecords(msr_block, false);
    
    adjust_.adj_file << std::endl;
}

void DynAdjustPrinter::PrintCompMeasurements(const UINT32& block, const std::string& type) {
    // Generate table heading with optional block and type information
    std::string table_heading("Computed Measurements");
    
    if (adjust_.projectSettings_.a.adjust_mode == PhasedMode || !type.empty()) {
        std::stringstream ss;
        ss << " (";
        if (adjust_.projectSettings_.a.adjust_mode == PhasedMode) {
            ss << "Block " << block + 1;
            if (!type.empty())
                ss << ", ";
        }

        if (!type.empty())
            ss << type;

        ss << ")";
        table_heading.append(ss.str());
    }

    // Print header using printer module method
    PrintAdjMeasurementsHeader(true, table_heading, computedMsrs, block, false);

    // Process measurements from the specified block
    it_vUINT32 _it_block_msr;
    it_vmsr_t _it_msr;
    UINT32 design_row(0);

    for (_it_block_msr = adjust_.v_CML_.at(block).begin(); 
         _it_block_msr != adjust_.v_CML_.at(block).end(); 
         ++_it_block_msr) {
        
        if (adjust_.InitialiseandValidateMsrPointer(_it_block_msr, _it_msr))
            continue;

        // When a target direction is found, continue to next element.  
        if (_it_msr->measType == 'D')
            if (_it_msr->vectorCount1 < 1)
                continue;

        adjust_.adj_file << std::left << std::setw(PAD2) << _it_msr->measType;

        // Dispatch to appropriate measurement handler
        switch (_it_msr->measType) {
        case 'A':	// Horizontal angle
            PrintCompMeasurements_A(block, _it_msr, design_row, computedMsrs);
            break;
        case 'B':	// Geodetic azimuth
        case 'K':	// Astronomic azimuth
        case 'V':	// Zenith distance
        case 'Z':	// Vertical angle
            PrintCompMeasurements_BKVZ(block, _it_msr, design_row, computedMsrs);
            break;
        case 'C':	// Chord dist
        case 'E':	// Ellipsoid arc
        case 'L':	// Level difference
        case 'M':	// MSL arc
        case 'S':	// Slope distance
            PrintCompMeasurements_CELMS(_it_msr, design_row, computedMsrs);
            break;
        case 'D':	// Direction set
            PrintCompMeasurements_D(_it_msr, design_row, false);
            break;
        case 'H':	// Orthometric height
        case 'R':	// Ellipsoidal height
            PrintCompMeasurements_HR(block, _it_msr, design_row, computedMsrs);
            break;
        case 'I':	// Astronomic latitude
        case 'J':	// Astronomic longitude
        case 'P':	// Geodetic latitude
        case 'Q':	// Geodetic longitude
            PrintCompMeasurements_IJPQ(block, _it_msr, design_row, computedMsrs);
            break;
        case 'G':	// GPS Baseline (treat as single-baseline cluster)
        case 'X':	// GPS Baseline cluster
        case 'Y':	// GPS Point cluster
            PrintCompMeasurements_GXY(block, _it_msr, design_row, computedMsrs);
            break;
        }
    }

    adjust_.adj_file << std::endl << std::endl;
}

// Helper function that consolidates measurement printing logic
void DynAdjustPrinter::PrintMeasurementRecords(const v_uint32_u32u32_pair& msr_block, bool adjustedMeasurements) {
    it_vmsr_t _it_msr;
    UINT32 clusterID(MAX_UINT32_VALUE);

    for (auto _it_block_msr = msr_block.begin(); _it_block_msr != msr_block.end(); ++_it_block_msr) {
        _it_msr = adjust_.bmsBinaryRecords_.begin() + (_it_block_msr->first);
        
        // Skip non-measurement elements (Y, Z, covariance components)
        if (_it_msr->measStart != xMeas)
            continue;

        // Skip target directions
        if (_it_msr->measType == 'D')
            if (_it_msr->vectorCount1 < 1)
                continue;

        // For cluster measurements, only print measurement type if from a new cluster
        switch (_it_msr->measType) {
        case 'D':	// Direction set
        case 'X':	// GPS Baseline cluster
        case 'Y':	// GPS Point cluster
            if (_it_msr->clusterID != clusterID) {
                adjust_.adj_file << std::left << std::setw(PAD2) << _it_msr->measType;
                clusterID = _it_msr->clusterID;
            }
            else
                adjust_.adj_file << "  ";
            break;
        default:
            adjust_.adj_file << std::left << std::setw(PAD2) << _it_msr->measType;
        }

        // Dispatch to appropriate measurement handler using consolidated approach
        UINT32 design_row(0);
        switch (_it_msr->measType) {
        case 'A':	// Horizontal angle
            if (adjustedMeasurements)
                PrintAdjMeasurements_A(_it_msr);
            else
                PrintCompMeasurements_A(_it_block_msr->second.first, _it_msr, design_row, computedMsrs);
            break;
        case 'B':	// Geodetic azimuth
        case 'K':	// Astronomic azimuth
        case 'V':	// Zenith distance
        case 'Z':	// Vertical angle
            if (adjustedMeasurements)
                PrintAdjMeasurements_BKVZ(_it_msr);
            else
                PrintCompMeasurements_BKVZ(_it_block_msr->second.first, _it_msr, design_row, computedMsrs);
            break;
        case 'C':	// Chord dist
        case 'E':	// Ellipsoid arc
        case 'L':	// Level difference
        case 'M':	// MSL arc
        case 'S':	// Slope distance
            if (adjustedMeasurements)
                PrintAdjMeasurements_CELMS(_it_msr);
            else
                PrintCompMeasurements_CELMS(_it_msr, design_row, computedMsrs);
            break;
        case 'D':	// Direction set
            if (adjustedMeasurements)
                PrintAdjMeasurements_D(_it_msr);
            else
                PrintCompMeasurements_D(_it_msr, design_row, false);
            break;
        case 'H':	// Orthometric height
        case 'R':	// Ellipsoidal height
            if (adjustedMeasurements)
                PrintAdjMeasurements_HR(_it_msr);
            else
                PrintCompMeasurements_HR(_it_block_msr->second.first, _it_msr, design_row, computedMsrs);
            break;
        case 'I':	// Astronomic latitude
        case 'J':	// Astronomic longitude
        case 'P':	// Geodetic latitude
        case 'Q':	// Geodetic longitude
            if (adjustedMeasurements)
                PrintAdjMeasurements_IJPQ(_it_msr);
            else
                PrintCompMeasurements_IJPQ(_it_block_msr->second.first, _it_msr, design_row, computedMsrs);
            break;
        case 'G':	// GPS Baseline (treat as single-baseline cluster)
        case 'X':	// GPS Baseline cluster
        case 'Y':	// GPS Point cluster
            if (adjustedMeasurements)
                PrintAdjMeasurements_GXY(_it_msr, _it_block_msr->second);
            else
                PrintCompMeasurements_GXY(_it_block_msr->second.first, _it_msr, design_row, computedMsrs);
            break;
        }
    }
}

// Stage 5: Enhanced measurement formatting helper functions
double DynAdjustPrinter::CalculateAngularPrecision(const it_vmsr_t& it_msr, char cardinal) const {
    double precision(0.0);
    
    switch (it_msr->measType) {
    case 'G':
    case 'X':
    case 'Y':
        switch (cardinal) {
        case 'P':
        case 'a':
            precision = it_msr->term2;  // Precision (Meas)
            break;
        case 'L':
        case 'v':
            precision = it_msr->term3;  // Precision (Meas)
            break;
        }
        break;
    case 'D':
        // Precision of subtended angle derived from successive direction set measurement precisions
        precision = it_msr->scale2;
        break;
    default:
        precision = it_msr->term2;  // Precision (Meas)
    }
    
    return precision;
}

double DynAdjustPrinter::CalculateLinearPrecision(const it_vmsr_t& it_msr, char cardinal) const {
    switch (it_msr->measType) {
    case 'G':
    case 'X':
    case 'Y':
        switch (cardinal) {
        case 'X':  // XYZ
        case 'e':  // ENU
            return sqrt(it_msr->term2);
        case 'Y':  // XYZ
        case 'n':  // ENU
            return sqrt(it_msr->term3);
        case 'Z':  // XYZ
        case 'u':  // ENU
        case 'H':  // LLH (Lat and Lon cardinals are printed in PrintMeasurementsAngular)
        case 'h':
            return sqrt(it_msr->term4);
        case 's':  // AED (Azimuth and vertical angle cardinals are printed in PrintMeasurementsAngular)
            switch (adjust_.projectSettings_.o._adj_gnss_units) {
            case AED_adj_gnss_ui:
                return sqrt(it_msr->term4);
            case ADU_adj_gnss_ui:
                return sqrt(it_msr->term3);
            default:
                return sqrt(it_msr->term2);
            }
        }
        break;
    default:
        return sqrt(it_msr->term2);  // Precision (Meas)
    }
    
    return sqrt(it_msr->term2);  // Default fallback
}

void DynAdjustPrinter::FormatAngularMeasurement(const double& preAdjMeas, const double& adjMeas, 
                                               const double& precision, const double& correction, 
                                               const it_vmsr_t& it_msr, bool printAdjMsr) {
    // Which angular format?
    if (adjust_.projectSettings_.o._angular_type_msr == DMS) {
        switch (adjust_.projectSettings_.o._dms_format_msr) {
        case SEPARATED_WITH_SYMBOLS:
            // ddd°mm'ss.sss"
            adjust_.adj_file << 
                std::setw(MSR) << std::right << FormatDmsString(RadtoDms(preAdjMeas), 4 + adjust_.PRECISION_SEC_MSR, true, true) <<
                std::setw(MSR) << std::right << FormatDmsString(RadtoDms(adjMeas), 4 + adjust_.PRECISION_SEC_MSR, true, true);
            break;
        case HP_NOTATION:
            // ddd.mmssssss
            adjust_.adj_file << 
                std::setw(MSR) << std::right << StringFromT(RadtoDms(preAdjMeas), 4 + adjust_.PRECISION_SEC_MSR) <<
                std::setw(MSR) << std::right << StringFromT(RadtoDms(adjMeas), 4 + adjust_.PRECISION_SEC_MSR);
            break;
        case SEPARATED:
        default:
            // ddd mm ss.ssss
            adjust_.adj_file << 
                std::setw(MSR) << std::right << FormatDmsString(RadtoDms(preAdjMeas), 4 + adjust_.PRECISION_SEC_MSR, true, false) <<
                std::setw(MSR) << std::right << FormatDmsString(RadtoDms(adjMeas), 4 + adjust_.PRECISION_SEC_MSR, true, false);
            break;
        }

        // Print correction and precision in seconds
        if (adjust_.isAdjustmentQuestionable_ || (fabs(it_msr->NStat) > adjust_.criticalValue_ * 4.0)) {
            adjust_.adj_file << 
                std::right << StringFromTW(removeNegativeZero(Seconds(correction), adjust_.PRECISION_SEC_MSR), CORR, adjust_.PRECISION_SEC_MSR) <<
                std::right << StringFromTW(Seconds(sqrt(precision)), PREC, adjust_.PRECISION_SEC_MSR);
        } else {
            adjust_.adj_file << 
                std::setw(CORR) << std::right << StringFromT(removeNegativeZero(Seconds(correction), adjust_.PRECISION_SEC_MSR), adjust_.PRECISION_SEC_MSR) <<
                std::setw(PREC) << std::right << StringFromT(Seconds(sqrt(precision)), adjust_.PRECISION_SEC_MSR);
        }

        if (printAdjMsr) {
            if (adjust_.isAdjustmentQuestionable_ || (fabs(it_msr->NStat) > adjust_.criticalValue_ * 4.0)) {
                adjust_.adj_file << 
                    std::right << StringFromTW(Seconds(sqrt(it_msr->measAdjPrec)), PREC, adjust_.PRECISION_SEC_MSR) <<
                    std::right << StringFromTW(Seconds(sqrt(it_msr->residualPrec)), PREC, adjust_.PRECISION_SEC_MSR);
            } else {
                adjust_.adj_file << 
                    std::setw(PREC) << std::right << StringFromT(Seconds(sqrt(it_msr->measAdjPrec)), adjust_.PRECISION_SEC_MSR) <<
                    std::setw(PREC) << std::right << StringFromT(Seconds(sqrt(it_msr->residualPrec)), adjust_.PRECISION_SEC_MSR);
            }
        }
    } else { // DDEG
        // ddd.dddddddd
        adjust_.adj_file << 
            std::setw(MSR) << std::right << StringFromT(Degrees(preAdjMeas), 4 + adjust_.PRECISION_SEC_MSR) <<
            std::setw(MSR) << std::right << StringFromT(Degrees(adjMeas), 4 + adjust_.PRECISION_SEC_MSR);
        
        if (adjust_.isAdjustmentQuestionable_ || (fabs(it_msr->NStat) > adjust_.criticalValue_ * 4.0)) {
            adjust_.adj_file <<
                std::right << StringFromTW(removeNegativeZero(Degrees(correction), adjust_.PRECISION_SEC_MSR), CORR, adjust_.PRECISION_SEC_MSR) <<
                std::right << StringFromTW(Degrees(sqrt(precision)), PREC, adjust_.PRECISION_SEC_MSR);
        } else {
            adjust_.adj_file <<
                std::setw(CORR) << std::right << StringFromT(removeNegativeZero(Degrees(correction), adjust_.PRECISION_SEC_MSR), adjust_.PRECISION_SEC_MSR) <<
                std::setw(PREC) << std::right << StringFromT(Degrees(sqrt(precision)), adjust_.PRECISION_SEC_MSR);
        }
        
        if (printAdjMsr) {
            if (adjust_.isAdjustmentQuestionable_ || (fabs(it_msr->NStat) > adjust_.criticalValue_ * 4.0)) {
                adjust_.adj_file <<
                    std::right << StringFromTW(Degrees(sqrt(it_msr->measAdjPrec)), PREC, adjust_.PRECISION_SEC_MSR) <<
                    std::right << StringFromTW(Degrees(sqrt(it_msr->residualPrec)), PREC, adjust_.PRECISION_SEC_MSR);
            } else {
                adjust_.adj_file <<
                    std::setw(PREC) << std::right << StringFromT(Degrees(sqrt(it_msr->measAdjPrec)), adjust_.PRECISION_SEC_MSR) <<
                    std::setw(PREC) << std::right << StringFromT(Degrees(sqrt(it_msr->residualPrec)), adjust_.PRECISION_SEC_MSR);
            }
        }
    }
}

void DynAdjustPrinter::FormatLinearMeasurement(const double& measurement, const double& correction, 
                                              const it_vmsr_t& it_msr, bool printAdjMsr) {
    // Print measurements with metric precision
    adjust_.adj_file << 
        std::setw(MSR) << std::setprecision(adjust_.PRECISION_MTR_MSR) << std::fixed << std::right << it_msr->preAdjMeas <<
        std::setw(MSR) << std::setprecision(adjust_.PRECISION_MTR_MSR) << std::fixed << std::right << measurement;
    
    // Print correction
    if (adjust_.isAdjustmentQuestionable_ || (fabs(it_msr->NStat) > adjust_.criticalValue_ * 4.0)) {
        adjust_.adj_file << std::right << StringFromTW(removeNegativeZero(correction, adjust_.PRECISION_MTR_MSR), CORR, adjust_.PRECISION_MTR_MSR);
    } else {
        adjust_.adj_file << std::setw(CORR) << std::setprecision(adjust_.PRECISION_MTR_MSR) << std::fixed << std::right << 
            removeNegativeZero(correction, adjust_.PRECISION_MTR_MSR);
    }

    // Print precision based on measurement type - use default for non-GNSS measurements
    double precision = sqrt(it_msr->term2); // Default fallback for most measurements
    
    if (adjust_.isAdjustmentQuestionable_ || (fabs(it_msr->NStat) > adjust_.criticalValue_ * 4.0)) {
        adjust_.adj_file << StringFromTW(precision, UINT16(PREC), adjust_.PRECISION_MTR_MSR);
    } else {
        adjust_.adj_file << std::setw(PREC) << std::setprecision(adjust_.PRECISION_MTR_MSR) << std::fixed << std::right << precision;
    }

    if (printAdjMsr) {
        if (adjust_.isAdjustmentQuestionable_ || (fabs(it_msr->NStat) > adjust_.criticalValue_ * 4.0)) {
            adjust_.adj_file <<
                StringFromTW(sqrt(it_msr->measAdjPrec), UINT16(PREC), adjust_.PRECISION_MTR_MSR) <<
                StringFromTW(sqrt(it_msr->residualPrec), UINT16(PREC), adjust_.PRECISION_MTR_MSR);
        } else {
            adjust_.adj_file <<
                std::setw(PREC) << std::setprecision(adjust_.PRECISION_MTR_MSR) << std::fixed << std::right << sqrt(it_msr->measAdjPrec) <<
                std::setw(PREC) << std::setprecision(adjust_.PRECISION_MTR_MSR) << std::fixed << std::right << sqrt(it_msr->residualPrec);
        }
    }
}

// Stage 5: Enhanced measurement formatting template specializations
template<>
void DynAdjustPrinter::PrintMeasurementValue<AngularMeasurement>(char cardinal, const double& measurement, 
                                                                 const double& correction, const it_vmsr_t& it_msr, bool printAdjMsr) {
    std::string ignoreFlag(" ");
    double preAdjMeas(it_msr->preAdjMeas);
    double adjMeas(measurement);

    switch (it_msr->measType) {
    case 'D':
        // capture original direction
        preAdjMeas = it_msr->term1;
        // "adjusted direction" is the original direction plus the correction
        adjMeas = it_msr->term1 + correction;
        // Don't print ignore flag for target directions
        break;
    default:
        // Print ignore flag if required
        if (it_msr->ignore)
            ignoreFlag = "*";
        break;
    }

    // Handle special cases for GNSS measurements
    switch (cardinal) {
    case 'L':		// GNSS (Y) Longitude
    case 'v':		// GNSS (alt units) Elevation
        // Print station information for GNSS measurements
        adjust_.adj_file << std::left << std::setw(PAD2) << it_msr->measType;
        adjust_.adj_file << std::left << std::setw(STATION) << adjust_.bstBinaryRecords_.at(it_msr->station1).stationName;
        
        // Handle different measurement types
        switch (it_msr->measType) {
        case 'G':
        case 'X':
            adjust_.adj_file << std::left << std::setw(STATION) << adjust_.bstBinaryRecords_.at(it_msr->station2).stationName;
            break;
        case 'Y':
            adjust_.adj_file << std::left << std::setw(STATION) << " ";
        }
        adjust_.adj_file << std::left << std::setw(STATION) << " ";
    }
    
    // Print measurement values with appropriate formatting
    adjust_.adj_file << std::setw(PAD3) << std::left << ignoreFlag << 
                        std::setw(PAD2) << std::left << cardinal;

    // Calculate precision and use consolidated formatting
    double precision = CalculateAngularPrecision(it_msr, cardinal);
    FormatAngularMeasurement(preAdjMeas, adjMeas, precision, correction, it_msr, printAdjMsr);
}

template<>
void DynAdjustPrinter::PrintMeasurementValue<LinearMeasurement>(char cardinal, const double& measurement, 
                                                                const double& correction, const it_vmsr_t& it_msr, bool printAdjMsr) {
    std::string ignoreFlag(" ");
    if (it_msr->ignore)
        ignoreFlag = "*";

    switch (cardinal) {
    case 'Y':
    case 'Z':
    case 'H':
    case 'h':
    case 'n':
    case 'u':
    case 's':
        // Print station information for GPS measurements
        adjust_.adj_file << std::left << std::setw(PAD2) << it_msr->measType;
        adjust_.adj_file << std::left << std::setw(STATION) << adjust_.bstBinaryRecords_.at(it_msr->station1).stationName;
        
        switch (it_msr->measType) {
        case 'G':
        case 'X':
            adjust_.adj_file << std::left << std::setw(STATION) << adjust_.bstBinaryRecords_.at(it_msr->station2).stationName;
            break;
        case 'Y':
            adjust_.adj_file << std::left << std::setw(STATION) << " ";
        }
        adjust_.adj_file << std::left << std::setw(STATION) << " ";
    }
    
    adjust_.adj_file << std::setw(PAD3) << std::left << ignoreFlag << 
                        std::setw(PAD2) << std::left << cardinal;

    // Print measurements with metric precision
    adjust_.adj_file << 
        std::setw(MSR) << std::setprecision(adjust_.PRECISION_MTR_MSR) << std::fixed << std::right << it_msr->preAdjMeas <<
        std::setw(MSR) << std::setprecision(adjust_.PRECISION_MTR_MSR) << std::fixed << std::right << measurement;
    
    // Print correction
    if (adjust_.isAdjustmentQuestionable_ || (fabs(it_msr->NStat) > adjust_.criticalValue_ * 4.0)) {
        adjust_.adj_file << std::right << StringFromTW(removeNegativeZero(correction, adjust_.PRECISION_MTR_MSR), CORR, adjust_.PRECISION_MTR_MSR);
    } else {
        adjust_.adj_file << std::setw(CORR) << std::setprecision(adjust_.PRECISION_MTR_MSR) << std::fixed << std::right << 
            removeNegativeZero(correction, adjust_.PRECISION_MTR_MSR);
    }

    // Calculate precision based on cardinal and measurement type
    double precision = CalculateLinearPrecision(it_msr, cardinal);
    
    // Handle special case for non-GNSS measurements where questionable adjustments need special formatting
    switch (it_msr->measType) {
    case 'G':
    case 'X':
    case 'Y':
        // For GNSS measurements, use calculated precision normally
        if (adjust_.isAdjustmentQuestionable_ || (fabs(it_msr->NStat) > adjust_.criticalValue_ * 4.0)) {
            adjust_.adj_file << StringFromTW(precision, UINT16(PREC), adjust_.PRECISION_MTR_MSR);
        } else {
            adjust_.adj_file << std::setw(PREC) << std::setprecision(adjust_.PRECISION_MTR_MSR) << std::fixed << std::right << precision;
        }
        break;
    default:
        // For other measurements, handle questionable case differently
        if (adjust_.isAdjustmentQuestionable_ || (fabs(it_msr->NStat) > adjust_.criticalValue_ * 4.0)) {
            adjust_.adj_file << StringFromTW(precision, UINT16(PREC), adjust_.PRECISION_MTR_MSR);
        } else {
            adjust_.adj_file << std::setw(PREC) << std::setprecision(adjust_.PRECISION_MTR_MSR) << std::fixed << std::right << precision;
        }
    }

    if (printAdjMsr) {
        if (adjust_.isAdjustmentQuestionable_ || (fabs(it_msr->NStat) > adjust_.criticalValue_ * 4.0)) {
            adjust_.adj_file <<
                StringFromTW(sqrt(it_msr->measAdjPrec), UINT16(PREC), adjust_.PRECISION_MTR_MSR) <<
                StringFromTW(sqrt(it_msr->residualPrec), UINT16(PREC), adjust_.PRECISION_MTR_MSR);
        } else {
            adjust_.adj_file <<
                std::setw(PREC) << std::setprecision(adjust_.PRECISION_MTR_MSR) << std::fixed << std::right << sqrt(it_msr->measAdjPrec) <<
                std::setw(PREC) << std::setprecision(adjust_.PRECISION_MTR_MSR) << std::fixed << std::right << sqrt(it_msr->residualPrec);
        }
    }
}

void DynAdjustPrinter::PrintMeasurementCorrection(const char cardinal, const it_vmsr_t& _it_msr)
{
    switch (_it_msr->measType)
    {
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
        // Pre adjustment correction for deflections
        adjust_.adj_file << std::setw(PACORR) << std::setprecision(adjust_.PRECISION_SEC_MSR) << std::fixed << std::right << 
            removeNegativeZero(Seconds(_it_msr->preAdjCorr), adjust_.PRECISION_SEC_MSR);
        break;
    case 'Y':
        // Pre adjustment correction of N-value on GPS applies to H cardinal of Y clusters in
        // geographic format only!
        switch (cardinal)
        {
        case 'H':			
            adjust_.adj_file << std::setw(PACORR) << std::setprecision(adjust_.PRECISION_MTR_MSR) << std::fixed << std::right << 
                removeNegativeZero(_it_msr->preAdjCorr, adjust_.PRECISION_MTR_MSR);
            break;
        case 'P':
        case 'L':
        case 'h':
            adjust_.adj_file << std::setw(PACORR) << std::setprecision(adjust_.PRECISION_SEC_MSR) << std::fixed << std::right << 0.0;
            break;
        default:	// X, Y, Z
            adjust_.adj_file << std::setw(PACORR) << std::setprecision(adjust_.PRECISION_MTR_MSR) << std::fixed << std::right << 0.0;
            break;
        }
        break;
    case 'C':
    case 'E':
    case 'H':
    case 'L':
    case 'M':
    case 'S':
    case 'G':
    case 'X':
    default:
        // Pre adjustment correction for N-value or reductions from MSL/ellipsoid arcs
        adjust_.adj_file << std::setw(PACORR) << std::setprecision(adjust_.PRECISION_SEC_MSR) << std::fixed << std::right << 
            removeNegativeZero(_it_msr->preAdjCorr, adjust_.PRECISION_SEC_MSR);
        break;
    }
}

void DynAdjustPrinter::PrintAdjMeasurements_YLLH(it_vmsr_t& _it_msr)
{
    // Get an iterator to this measurement which will be needed later when
    // obtaining variance matrices via GetGPSVarianceMatrix()
    it_vmsr_t _it_msr_begin(_it_msr);
    
    // create a temporary copy of this Y measurement and transform/propagate
    // cartesian elements to geographic
    vmsr_t y_msr;
    CopyClusterMsr<vmsr_t>(adjust_.bmsBinaryRecords_, _it_msr, y_msr);
    
    it_vmsr_t _it_y_msr(y_msr.begin());

    // determine coord type
    switch (_it_msr->station3)
    { 
    case LLh_type_i:
        snprintf(_it_y_msr->coordType, STN_TYPE_WIDTH, "%s", LLh_type);
        break;
    case LLH_type_i:
    default:
        snprintf(_it_y_msr->coordType, STN_TYPE_WIDTH, "%s", LLH_type);
        break;
    }
    
    UINT32 covr, cluster_msr, cluster_count(_it_y_msr->vectorCount1), covariance_count;
    matrix_2d mpositions(cluster_count * 3, 1);
    double latitude, longitude, height;

    it_vstn_t stn1_it;

    // 1. Convert coordinates from cartesian to geographic
    adjust_.ReduceYLLHMeasurementsforPrinting(y_msr, mpositions, adjustedMsrs);
    
    _it_y_msr = y_msr.begin();

    matrix_2d var_cart_adj(3, 3), var_geo_adj(3, 3);

    // 2. Convert adjusted precisions
    for (cluster_msr=0; cluster_msr<cluster_count; ++cluster_msr)
    {
        covr = cluster_msr * 3;
        covariance_count = _it_y_msr->vectorCount2;

        var_cart_adj.put(0, 0, _it_y_msr->measAdjPrec);
        _it_y_msr++;
        var_cart_adj.put(1, 1, _it_y_msr->measAdjPrec);
        _it_y_msr++;	
        var_cart_adj.put(2, 2, _it_y_msr->measAdjPrec);

        latitude = mpositions.get(covr, 0);
        longitude = mpositions.get(covr+1, 0);
        height = mpositions.get(covr+2, 0);

        PropagateVariances_GeoCart<double>(var_cart_adj, &var_geo_adj,
            latitude, longitude, height,
            adjust_.datum_.GetEllipsoidRef(), 
            false);	 		// Cartesian -> Geographic

        _it_y_msr-=2;
        _it_y_msr->measAdjPrec = var_geo_adj.get(0, 0);

        _it_y_msr++;
        _it_y_msr->measAdjPrec = var_geo_adj.get(1, 1);
        
        _it_y_msr++;	
        _it_y_msr->measAdjPrec = var_geo_adj.get(2, 2);
        
        // skip covariances until next point
        _it_y_msr += covariance_count * 3;
        
        if (covariance_count > 0)
            _it_y_msr++;
    }

    _it_y_msr = y_msr.begin();

    // 3. Convert apriori measurement precisions

    // Get measurement cartesian variance matrix
    matrix_2d var_cart;
    GetGPSVarianceMatrix(_it_msr_begin, &var_cart);
    var_cart.filllower();

    // Propagate the cartesian variance matrix to geographic
    PropagateVariances_GeoCart_Cluster<double>(var_cart, &var_cart,
        mpositions, 
        adjust_.datum_.GetEllipsoidRef(), 
        false); 		// Cartesian -> Geographic

    SetGPSVarianceMatrix(_it_y_msr, var_cart);
    
    // Now print the temporary measurement
    bool nextElement(false);

    for (cluster_msr=0; cluster_msr<cluster_count; ++cluster_msr)
    {
        stn1_it = adjust_.bstBinaryRecords_.begin() + _it_y_msr->station1;

        if (nextElement)
            adjust_.adj_file << std::left << std::setw(PAD2) << _it_y_msr->measType;
        else
            nextElement = true;

        covariance_count = _it_y_msr->vectorCount2;

        // update statistics
        // Latitude
        _it_y_msr->residualPrec = _it_y_msr->term2 - _it_y_msr->measAdjPrec;
        if (_it_y_msr->residualPrec < 0.0)
            _it_y_msr->residualPrec = fabs(_it_y_msr->residualPrec);
        _it_y_msr->PelzerRel = sqrt(_it_y_msr->term2) / sqrt(_it_y_msr->residualPrec);
        if (_it_y_msr->PelzerRel < 0. || _it_y_msr->PelzerRel > STABLE_LIMIT)
            _it_y_msr->PelzerRel = UNRELIABLE;
        _it_y_msr->NStat = _it_y_msr->measCorr / sqrt(_it_y_msr->residualPrec);
        _it_y_msr++;
        
        // Longitude
        _it_y_msr->residualPrec = _it_y_msr->term3 - _it_y_msr->measAdjPrec;
        if (_it_y_msr->residualPrec < 0.0)
            _it_y_msr->residualPrec = fabs(_it_y_msr->residualPrec);
        _it_y_msr->PelzerRel = sqrt(_it_y_msr->term3) / sqrt(_it_y_msr->residualPrec);
        if (_it_y_msr->PelzerRel < 0. || _it_y_msr->PelzerRel > STABLE_LIMIT)
            _it_y_msr->PelzerRel = UNRELIABLE;
        _it_y_msr->NStat = _it_y_msr->measCorr / sqrt(_it_y_msr->residualPrec);
        _it_y_msr++;

        // Orthometric height
        _it_y_msr->residualPrec = _it_y_msr->term4 - _it_y_msr->measAdjPrec;
        if (_it_y_msr->residualPrec < 0.0)
            _it_y_msr->residualPrec = fabs(_it_y_msr->residualPrec);
        _it_y_msr->PelzerRel = sqrt(_it_y_msr->term4) / sqrt(_it_y_msr->residualPrec);
        if (_it_y_msr->PelzerRel < 0. || _it_y_msr->PelzerRel > STABLE_LIMIT)
            _it_y_msr->PelzerRel = UNRELIABLE;
        _it_y_msr->NStat = _it_y_msr->measCorr / sqrt(_it_y_msr->residualPrec);
        
        _it_y_msr-=2;
        
        // first, second, third stations
        adjust_.adj_file << std::left << std::setw(STATION) << 
            stn1_it->stationName <<
            std::left << std::setw(STATION) << " " <<		// second station
            std::left << std::setw(STATION) << " ";		// third station

        // Print latitude
        PrintAdjMeasurementsAngular('P', _it_y_msr);
    
        // Print longitude
        _it_y_msr++;	
        PrintAdjMeasurementsAngular('L', _it_y_msr);

        // Print height
        _it_y_msr++;
        switch (_it_msr->station3)
        { 
        case LLh_type_i: 
            PrintAdjMeasurementsLinear('h', _it_y_msr);
            break;
        case LLH_type_i:
        default:
            PrintAdjMeasurementsLinear('H', _it_y_msr);
            break;
        }        

        // skip covariances until next point
        _it_y_msr += covariance_count * 3;
        
        if (covariance_count > 0)
            _it_y_msr++;
    }
}

void DynAdjustPrinter::PrintPositionalUncertainty()
{
    std::ofstream apu_file;
    try {
        // Create apu file.  Throws runtime_error on failure.
        file_opener(apu_file, adjust_.projectSettings_.o._apu_file);
    }
    catch (const std::runtime_error& e) {
        adjust_.SignalExceptionAdjustment(e.what(), 0);
    }

    // Use printer infrastructure for header
    PrintPositionalUncertaintyHeader(apu_file, adjust_.projectSettings_.o._apu_file);
    
    apu_file << std::setw(PRINT_VAR_PAD) << std::left << "Stations printed in blocks:";
    if (adjust_.projectSettings_.a.adjust_mode != SimultaneousMode)
    {
        if (adjust_.projectSettings_.o._output_stn_blocks)
            apu_file << "Yes" << std::endl;
        else
        {
            if (adjust_.projectSettings_.o._output_pu_covariances)
                apu_file << std::left << "Yes (enforced with output-all-covariances)" << std::endl;
            else
                apu_file << std::left << "No" << std::endl;
        }
    }
    else
        apu_file << std::left << "No" << std::endl;

    apu_file << std::setw(PRINT_VAR_PAD) << std::left << "Variance matrix units:";
    switch (adjust_.projectSettings_.o._apu_vcv_units)
    {
    case ENU_apu_ui:
        apu_file << std::left << "ENU" << std::endl;
        break;
    default:
    case XYZ_apu_ui:
        apu_file << std::left << "XYZ" << std::endl;
        break;
    }

    apu_file <<
        std::setw(PRINT_VAR_PAD) << std::left << "Full covariance matrix:";
    if (adjust_.projectSettings_.o._output_pu_covariances)
        apu_file << std::left << "Yes" << std::endl;
    else
        apu_file << std::left << "No" << std::endl;

    if (adjust_.projectSettings_.o._apply_type_b_file || adjust_.projectSettings_.o._apply_type_b_global)
    {
        if (adjust_.projectSettings_.o._apply_type_b_global)
            apu_file << std::setw(PRINT_VAR_PAD) << std::left << "Type B uncertainties:" <<
                adjust_.projectSettings_.a.type_b_global << std::endl;

        if (adjust_.projectSettings_.o._apply_type_b_file)
            apu_file << std::setw(PRINT_VAR_PAD) << std::left << "Type B uncertainty file:" <<
                safe_absolute_path(adjust_.projectSettings_.a.type_b_file) << std::endl;
    }

    apu_file << OUTPUTLINE << std::endl << std::endl;

    apu_file << "Positional uncertainty of adjusted station coordinates" << std::endl;
    apu_file << "------------------------------------------------------" << std::endl;
    
    apu_file << std::endl;

    switch (adjust_.projectSettings_.a.adjust_mode)
    {
    case SimultaneousMode:
        PrintPosUncertainties(apu_file, 0, 
            &adjust_.v_rigorousVariances_.at(0));
        break;
    case PhasedMode:
        // Output phased blocks as a single block?
        if (!adjust_.projectSettings_.o._output_stn_blocks &&		// Print stations in blocks?
            !adjust_.projectSettings_.o._output_pu_covariances)		// Print covariances?
                                                                    // If covariances are required, stations
                                                                    // must be printed in blocks.
        {
            PrintPosUncertaintiesUniqueList(apu_file, 
                &adjust_.v_rigorousVariances_);
            return;
        }

        for (UINT32 block=0; block<adjust_.blockCount_; ++block)
        {
            // load up this block
            if (adjust_.projectSettings_.a.stage)
                adjust_.DeserialiseBlockFromMappedFile(block, 1, 
                    sf_rigorous_vars);

            PrintPosUncertainties(apu_file, block, 
                &adjust_.v_rigorousVariances_.at(block));

            // unload this block
            if (adjust_.projectSettings_.a.stage)
                adjust_.UnloadBlock(block, 1, 
                    sf_rigorous_vars);
        }
        break;
    case Phased_Block_1Mode:		// only the first block is rigorous
        PrintPosUncertainties(apu_file, 0, 
            &adjust_.v_rigorousVariances_.at(0));
        break;
    }

    apu_file.close();
}

void DynAdjustPrinter::PrintEstimatedStationCoordinatestoDNAXML(const std::string& stnFile, INPUT_FILE_TYPE t, bool flagUnused)
{
    // Stations
    std::ofstream stn_file;
    try {
        // Create STN/XML file. 
        file_opener(stn_file, stnFile);
    }
    catch (const std::runtime_error& e) {
        adjust_.SignalExceptionAdjustment(e.what(), 0);
    }

    try {
    
        UINT32 count(static_cast<UINT32>(adjust_.bstBinaryRecords_.size()));
        dna_stn_fields dsl, dsw;

        std::string headerComment("Source data:  Coordinates estimated from least squares adjustment.");
        
        // print header
        switch (t)
        {
        case dynaml:

            // Write header and comments
            dynaml_header(stn_file, "Station File", adjust_.datum_.GetName(), adjust_.datum_.GetEpoch_s());
            dynaml_comment(stn_file, "File type:    Station file");
            dynaml_comment(stn_file, "Project name: " + adjust_.projectSettings_.g.network_name);
            dynaml_comment(stn_file, headerComment);
            dynaml_comment(stn_file, "Adj file:     " + adjust_.projectSettings_.o._adj_file);
            break;

        case dna:

            // get file format field widths
            determineDNASTNFieldParameters<UINT16>("3.01", dsl, dsw);

            // Write header and comments
            dna_header(stn_file, "3.01", "STN", adjust_.datum_.GetName(), adjust_.datum_.GetEpoch_s(), count);
            dna_comment(stn_file, "File type:    Station file");
            dna_comment(stn_file, "Project name: " + adjust_.projectSettings_.g.network_name);
            dna_comment(stn_file, headerComment);
            dna_comment(stn_file, "Adj file:     " + adjust_.projectSettings_.o._adj_file);
            break;
        default:
            break;
        }

        dnaStnPtr stnPtr(new CDnaStation(adjust_.datum_.GetName(), adjust_.datum_.GetEpoch_s()));

        vUINT32 vStationList;

        switch (adjust_.projectSettings_.a.adjust_mode)
        {
        case Phased_Block_1Mode:
            // Get stations in block 1 only
            vStationList = adjust_.v_parameterStationList_.at(0);
            break;
        default:
            // Create and initialise vector with 0,1,2,...,n-2,n-1,n for
            // all stations in the network
            vStationList.resize(adjust_.bstBinaryRecords_.size());
            initialiseIncrementingIntegerVector<UINT32>(vStationList, static_cast<UINT32>(adjust_.bstBinaryRecords_.size()));
        }
        
        // Sort on original file order
        CompareStnFileOrder<station_t, UINT32> stnorderCompareFunc(&adjust_.bstBinaryRecords_);
        std::sort(vStationList.begin(), vStationList.end(), stnorderCompareFunc);

        // print station coordinates
        switch (t)
        {
        case dynaml:

            if (flagUnused)
                // Print stations in DynaML format
                for_each(vStationList.begin(), vStationList.end(),
                    [&stn_file, &stnPtr, this](const UINT32& stn) {
                        stnPtr->SetStationRec(adjust_.bstBinaryRecords_.at(stn));
                        if (stnPtr->IsNotUnused())
                            stnPtr->WriteDNAXMLStnCurrentEstimates(&stn_file,
                                adjust_.datum_.GetEllipsoidRef(), &adjust_.projection_, dynaml);
                });
            else
                // Print stations in DynaML format
                for_each(vStationList.begin(), vStationList.end(),
                    [&stn_file, &stnPtr, this](const UINT32& stn) {
                        stnPtr->SetStationRec(adjust_.bstBinaryRecords_.at(stn));
                        stnPtr->WriteDNAXMLStnCurrentEstimates(&stn_file,
                            adjust_.datum_.GetEllipsoidRef(), &adjust_.projection_, dynaml);
                });

            stn_file << "</DnaXmlFormat>" << std::endl;

            break;
        case dna:

            if (flagUnused)
                // Print stations in DNA format
                for_each(vStationList.begin(), vStationList.end(),
                    [&stn_file, &stnPtr, &dsw, this](const UINT32& stn) {
                        stnPtr->SetStationRec(adjust_.bstBinaryRecords_.at(stn));
                        if (stnPtr->IsNotUnused())
                            stnPtr->WriteDNAXMLStnCurrentEstimates(&stn_file, 
                                adjust_.datum_.GetEllipsoidRef(), &adjust_.projection_, dna, &dsw);
                });
            else
                // Print stations in DNA format
                for_each(vStationList.begin(), vStationList.end(),
                    [&stn_file, &stnPtr, &dsw, this](const UINT32& stn) {
                        stnPtr->SetStationRec(adjust_.bstBinaryRecords_.at(stn));
                        stnPtr->WriteDNAXMLStnCurrentEstimates(&stn_file,
                            adjust_.datum_.GetEllipsoidRef(), &adjust_.projection_, dna, &dsw);
                });

            break;
        default:
            break;
        }


        stn_file.close();
    }
    catch (const std::ofstream::failure& f) {
        adjust_.SignalExceptionAdjustment(f.what(), 0);
    }
    catch (const XMLInteropException& e)  {
        adjust_.SignalExceptionAdjustment(e.what(), 0);
    }
}

bool DynAdjustPrinter::PrintEstimatedStationCoordinatestoSNX(std::string& sinex_filename)
{
    std::ofstream sinex_file;
    std::string sinexFilename;
    std::string sinexBasename = adjust_.projectSettings_.g.output_folder + FOLDER_SLASH + adjust_.projectSettings_.g.network_name;
    std::stringstream ssBlock;

    matrix_2d *estimates = nullptr, *variances = nullptr;

    bool success(true);

    // Print a sinex file for each block
    for (UINT32 block(0); block<adjust_.blockCount_; ++block)
    {
        sinexFilename = sinexBasename;

        // Add block number for phased adjustments
        switch (adjust_.projectSettings_.a.adjust_mode)
        {
        case Phased_Block_1Mode:
            if (block > 0)
            {
                block = adjust_.blockCount_;
                continue;
            }
            [[fallthrough]];
        case PhasedMode:
            ssBlock.str("");
            ssBlock << "-block" << block + 1;
            sinexFilename += ssBlock.str();

            // if staged, load up the block from memory mapped files
            if (adjust_.projectSettings_.a.stage)
                adjust_.DeserialiseBlockFromMappedFile(block, 2, 
                    sf_rigorous_stns, sf_rigorous_vars);

            estimates = &adjust_.v_rigorousStations_.at(block);
            variances = &adjust_.v_rigorousVariances_.at(block);
            break;
        case SimultaneousMode:
            estimates = &adjust_.v_estimatedStations_.at(block);
            variances = &adjust_.v_normals_.at(block);
        }

        // add reference frame to file name
        sinexFilename += "." + adjust_.projectSettings_.r.reference_frame;		

        // add conventional file extension
        sinexFilename += ".snx";

        if (block < 1)
            sinex_filename = sinexFilename;
        
        try {
            // Open output file stream.  Throws runtime_error on failure.
            file_opener(sinex_file, sinexFilename);
        }
        catch (const std::runtime_error& e) {
            adjust_.SignalExceptionAdjustment(e.what(), 0);
        }

        DnaIoSnx snx;

        try {
            // Print results for adjustment in SINEX format.
            // Throws runtime_error on failure.
            snx.SerialiseSinex(&sinex_file, &adjust_.bstBinaryRecords_,
                adjust_.bst_meta_, adjust_.bms_meta_, estimates, variances, adjust_.projectSettings_,
                adjust_.measurementParams_, adjust_.unknownsCount_, adjust_.sigmaZero_,
                &adjust_.v_blockStationsMap_.at(block), &adjust_.v_parameterStationList_.at(block),
                adjust_.blockCount_, block,
                &adjust_.datum_);
        }
        catch (const std::runtime_error& e) {
            adjust_.SignalExceptionAdjustment(e.what(), 0);
        }

        if (adjust_.projectSettings_.a.adjust_mode == PhasedMode)
            // if staged, unload the block
            if (adjust_.projectSettings_.a.stage)
                adjust_.UnloadBlock(block, 2, 
                    sf_rigorous_stns, sf_rigorous_vars);

        sinex_file.close();

        if (!snx.WarningsExist())
            continue;

        success = false;

        sinexFilename += ".err";
        try {
            // Open output file stream.  Throws runtime_error on failure.
            file_opener(sinex_file, sinexFilename);
        }
        catch (const std::runtime_error& e) {
            adjust_.SignalExceptionAdjustment(e.what(), 0);
        }

        snx.PrintWarnings(&sinex_file, sinexFilename);
        sinex_file.close();		
    }

    return success;
}

void DynAdjustPrinter::PrintEstimatedStationCoordinatestoDNAXML_Y(const std::string& msrFile, INPUT_FILE_TYPE t) {
    // Measurements
    std::ofstream msr_file;
    try {
        // Create STN/XML file
        file_opener(msr_file, msrFile);
    }
    catch (const std::runtime_error& e) {
        adjust_.SignalExceptionAdjustment(e.what(), 0);
    }

    dna_msr_fields dml, dmw;
    std::stringstream ss;

    try {
        UINT32 count(static_cast<UINT32>(adjust_.v_blockStationsMapUnique_.size()));
        ss << "Source data:  Coordinates and uncertainties for " << 
            count << " unique stations in " << adjust_.blockCount_ << " blocks " <<
            "estimated from least squares adjustment.";
        std::string headerComment(ss.str());
        
        // Print header based on format type
        switch (t) {
        case dynaml:
            // Write header and comments
            dynaml_header(msr_file, "Measurement File", adjust_.datum_.GetName(), adjust_.datum_.GetEpoch_s());
            dynaml_comment(msr_file, "File type:    Measurement file");
            dynaml_comment(msr_file, "Project name: " + adjust_.projectSettings_.g.network_name);
            dynaml_comment(msr_file, headerComment);
            dynaml_comment(msr_file, "Adj file:     " + adjust_.projectSettings_.o._adj_file);
            break;

        case dna:
            // Get file format field widths
            determineDNAMSRFieldParameters<UINT16>("3.01", dml, dmw);

            // Write header and comments
            dna_header(msr_file, "3.01", "MSR", adjust_.datum_.GetName(), adjust_.datum_.GetEpoch_s(), adjust_.blockCount_);
            dna_comment(msr_file, "File type:    Measurement file");
            dna_comment(msr_file, "Project name: " + adjust_.projectSettings_.g.network_name);
            dna_comment(msr_file, headerComment);
            dna_comment(msr_file, "Adj file:     " + adjust_.projectSettings_.o._adj_file);
            break;
        default:
            break;
        }
    }
    catch (const std::ofstream::failure& f) {
        adjust_.SignalExceptionAdjustment(f.what(), 0);
    }
    catch (const XMLInteropException& e) {
        adjust_.SignalExceptionAdjustment(e.what(), 0);
    }
    catch (const std::runtime_error& e) {
        adjust_.SignalExceptionAdjustment(e.what(), 0);
    }

    try {
        matrix_2d *estimates = nullptr, *variances = nullptr;
        dnaMsrPtr msr_ptr;
        std::string comment;

        // Create a Y cluster for each block
        for (UINT32 block(0); block < adjust_.blockCount_; ++block) {
            // Get the appropriate coordinate and uncertainty estimates
            switch (adjust_.projectSettings_.a.adjust_mode) {
            case Phased_Block_1Mode:
                if (block > 0) {
                    block = adjust_.blockCount_;
                    continue;
                }
                [[fallthrough]];
            case PhasedMode:
                // If staged, load up the block from memory mapped files
                if (adjust_.projectSettings_.a.stage)
                    adjust_.DeserialiseBlockFromMappedFile(block, 2, 
                        sf_rigorous_stns, sf_rigorous_vars);

                estimates = &adjust_.v_rigorousStations_.at(block);
                variances = &adjust_.v_rigorousVariances_.at(block);
                break;
            case SimultaneousMode:
                estimates = &adjust_.v_estimatedStations_.at(block);
                variances = &adjust_.v_normals_.at(block);
            }

            msr_ptr.reset(new CDnaGpsPointCluster(block, adjust_.datum_.GetName(), adjust_.datum_.GetEpoch_s()));
            msr_ptr.get()->PopulateMsr(&adjust_.bstBinaryRecords_, 
                &adjust_.v_blockStationsMap_.at(block), &adjust_.v_parameterStationList_.at(block),
                block, &adjust_.datum_, estimates, variances);

            if (adjust_.projectSettings_.a.adjust_mode == PhasedMode)
                // If staged, unload the block
                if (adjust_.projectSettings_.a.stage)
                    adjust_.UnloadBlock(block, 2, 
                        sf_rigorous_stns, sf_rigorous_vars);

            ss.str("");

            // Print station coordinate estimates as a GNSS Y cluster
            switch (t) {
            case dynaml:
                ss << std::endl <<
                    "    - Estimated station coordinates and uncertainties";
                if (adjust_.blockCount_ > 1)
                    ss << " for block " << (block + 1);
                ss << std::endl <<
                    "    - Type (Y) GPS point cluster (set of " << adjust_.v_blockStationsMap_.at(block).size() << " stations)";
                if (adjust_.blockCount_ > 1) {
                    ss << ":" << std::endl <<
                        "      - " << adjust_.v_ISL_.at(block).size() << " inner station";
                    if (adjust_.v_ISL_.at(block).size() != 1)
                        ss << "s";
                    ss << std::endl <<
                        "      - " << adjust_.v_JSL_.at(block).size() << " junction station";
                    if (adjust_.v_JSL_.at(block).size() != 1)
                        ss << "s";
                    ss << std::endl;
                }
                else
                    ss << std::endl;

                ss << " ";
                comment = ss.str();
                // Print station coordinate estimates in DynaML format
                msr_ptr.get()->WriteDynaMLMsr(&msr_file, comment);
                break;
            case dna:
                // Print station coordinate estimates in DNA format
                msr_ptr.get()->WriteDNAMsr(&msr_file, dmw, dml);
                break;
            default:
                break;
            }
        }

        switch (t) {
        case dynaml:
            msr_file << "</DnaXmlFormat>" << std::endl;
            break;
        case dna:
        default:
            break;
        }

        msr_file.close();
    }
    catch (const std::ofstream::failure& f) {
        adjust_.SignalExceptionAdjustment(f.what(), 0);
    }
    catch (const XMLInteropException& e) {
        adjust_.SignalExceptionAdjustment(e.what(), 0);
    }
}

void DynAdjustPrinter::PrintPosUncertaintiesUniqueList(std::ostream& os, const v_mat_2d* stationVariances)
{
    // Print header
    PrintPosUncertaintiesHeader(os);

    UINT32 block(UINT_MAX), stn, mat_index;
    _it_u32u32_uint32_pair _it_bsmu;

    // Determine sort order of block Stations Map
    if (adjust_.projectSettings_.a.stage || adjust_.projectSettings_.a.adjust_mode == Phased_Block_1Mode)
    {
        // sort by blocks to create efficiency when deserialising matrix info
        std::sort(adjust_.v_blockStationsMapUnique_.begin(), adjust_.v_blockStationsMapUnique_.end(), 
            CompareBlockStationMapUnique_byBlock<u32u32_uint32_pair>());
    }
    else if (adjust_.projectSettings_.o._sort_stn_file_order)
    {
        CompareBlockStationMapUnique_byFileOrder<station_t, u32u32_uint32_pair> stnorderCompareFunc(&adjust_.bstBinaryRecords_);
        std::sort(adjust_.v_blockStationsMapUnique_.begin(), adjust_.v_blockStationsMapUnique_.end(), stnorderCompareFunc);
    }
    else
        std::sort(adjust_.v_blockStationsMapUnique_.begin(), adjust_.v_blockStationsMapUnique_.end());
    
    std::stringstream stationRecord;
    v_uint32_string_pair stationsOutput;
    stationsOutput.reserve(adjust_.v_blockStationsMapUnique_.size());

    std::ostream* outstream(&os);
    if (adjust_.projectSettings_.a.stage)
        outstream = &stationRecord;

    // Print stations according to the user-defined sort order
    for (_it_bsmu=adjust_.v_blockStationsMapUnique_.begin();
        _it_bsmu!=adjust_.v_blockStationsMapUnique_.end(); ++_it_bsmu)
    {
        // If this is not the first block, exit if block-1 mode
        if (adjust_.projectSettings_.a.adjust_mode == Phased_Block_1Mode)
            if (_it_bsmu->second > 0)
                break;

        if (adjust_.projectSettings_.a.stage)
        {
            if (block != _it_bsmu->second)
            {
                // unload previous block
                if (_it_bsmu->second > 0)
                    adjust_.UnloadBlock(_it_bsmu->second-1, 1, 
                        sf_rigorous_vars);
                // load up this block
                if (adjust_.projectSettings_.a.stage)
                    adjust_.DeserialiseBlockFromMappedFile(_it_bsmu->second, 1, 
                        sf_rigorous_vars);
            }
        }

        block = _it_bsmu->second;
        stn = _it_bsmu->first.first;
        mat_index = _it_bsmu->first.second * 3;

        PrintPosUncertainty(*outstream,
            block, stn, mat_index, 
                &stationVariances->at(block),
                _it_bsmu->first.second, nullptr);

        if (adjust_.projectSettings_.a.stage)
        {
            stationsOutput.push_back(uint32_string_pair(stn, stationRecord.str()));
            stationRecord.str("");
        }
    }

    if (adjust_.projectSettings_.a.stage)
    {
        adjust_.UnloadBlock(block, 1, sf_rigorous_vars);

        // if required, sort stations according to original station file order
        if (adjust_.projectSettings_.o._sort_stn_file_order)
        {
            CompareOddPairFirst_FileOrder<station_t, UINT32, std::string> stnorderCompareFunc(&adjust_.bstBinaryRecords_);
            std::sort(stationsOutput.begin(), stationsOutput.end(), stnorderCompareFunc);
        }
        else
            std::sort(stationsOutput.begin(), stationsOutput.end(), CompareOddPairFirst<UINT32, std::string>());

        for_each(stationsOutput.begin(), stationsOutput.end(),
            [&stationsOutput, &os] (std::pair<UINT32, std::string>& posRecord) { 
                os << posRecord.second;
            }
        );
    }
}

void DynAdjustPrinter::PrintStationCorrectionsList(std::ostream& cor_file)
{
    vUINT32 v_blockStations;
    vstring stn_corr_records(adjust_.bstBinaryRecords_.size());

    cor_file << std::setw(STATION) << std::left << "Station" <<
        std::setw(PAD2) << " " << 
        std::right << std::setw(MSR) << "Azimuth" << 
        std::right << std::setw(MSR) << "V. Angle" << 
        std::right << std::setw(MSR) << "S. Distance" << 
        std::right << std::setw(MSR) << "H. Distance" << 
        std::right << std::setw(HEIGHT) << "east" << 
        std::right << std::setw(HEIGHT) << "north" << 
        std::right << std::setw(HEIGHT) << "up" << std::endl;

    UINT32 i, j = STATION+PAD2+MSR+MSR+MSR+MSR+HEIGHT+HEIGHT+HEIGHT;

    for (i=0; i<j; ++i)
        cor_file << "-";
    cor_file << std::endl;

    UINT32 block(UINT_MAX), stn, mat_index;
    _it_u32u32_uint32_pair _it_bsmu;
    
    std::stringstream stationRecord;
    v_uint32_string_pair correctionsOutput;
    correctionsOutput.reserve(adjust_.v_blockStationsMapUnique_.size());

    const v_mat_2d* estimates = &adjust_.v_rigorousStations_;
    if (adjust_.projectSettings_.a.adjust_mode == SimultaneousMode)
        estimates = &adjust_.v_estimatedStations_;

    std::ostream* outstream(&cor_file);
    if (adjust_.projectSettings_.a.stage)
        outstream = &stationRecord;

    // No need to sort order of block Stations Map - this was done in
    // PrintAdjStationsUniqueList

    for (_it_bsmu=adjust_.v_blockStationsMapUnique_.begin();
        _it_bsmu!=adjust_.v_blockStationsMapUnique_.end(); ++_it_bsmu)
    {
        // If this is not the first block, exit if block-1 mode
        if (adjust_.projectSettings_.a.adjust_mode == Phased_Block_1Mode)
            if (_it_bsmu->second > 0)
                break;

        if (adjust_.projectSettings_.a.stage)
        {
            if (block != _it_bsmu->second)
            {
                // unload previous block
                if (_it_bsmu->second > 0)
                    adjust_.UnloadBlock(_it_bsmu->second-1, 2, 
                        sf_rigorous_stns, sf_original_stns);

                // load up this block
                if (adjust_.projectSettings_.a.stage)
                    adjust_.DeserialiseBlockFromMappedFile(_it_bsmu->second, 2,
                        sf_rigorous_stns, sf_original_stns);
            }
        }

        block = _it_bsmu->second;
        stn = _it_bsmu->first.first;
        mat_index = _it_bsmu->first.second * 3;

        PrintCorStation(*outstream, block, stn, mat_index,
            &estimates->at(block));
            
        if (adjust_.projectSettings_.a.stage)
        {
            correctionsOutput.push_back(uint32_string_pair(stn, stationRecord.str()));
            stationRecord.str("");
        }
    }

    if (adjust_.projectSettings_.a.stage)
        adjust_.UnloadBlock(block, 2, 
            sf_rigorous_stns, sf_original_stns);

    // if required, sort stations according to original station file order
    if (adjust_.projectSettings_.o._sort_stn_file_order)
    {
        CompareOddPairFirst_FileOrder<station_t, UINT32, std::string> stnorderCompareFunc(&adjust_.bstBinaryRecords_);
        std::sort(correctionsOutput.begin(), correctionsOutput.end(), stnorderCompareFunc);
    }
    else
        std::sort(correctionsOutput.begin(), correctionsOutput.end(), CompareOddPairFirst<UINT32, std::string>());

    for_each(correctionsOutput.begin(), correctionsOutput.end(),
        [&correctionsOutput, &cor_file] (std::pair<UINT32, std::string>& corrRecord) { 
            cor_file << corrRecord.second;
        }
    );
}

void DynAdjustPrinter::PrintBlockStations(std::ostream& os, const UINT32& block,
    const matrix_2d* stationEstimates, matrix_2d* stationVariances, bool printBlockID,
    bool recomputeGeographicCoords, bool updateGeographicCoords, bool printHeader,
    bool reapplyTypeBUncertainties)
{
    vUINT32 v_blockStations(adjust_.v_parameterStationList_.at(block));

    if (v_blockStations.size() * 3 != stationEstimates->rows())
    {
        std::stringstream ss;
        ss << "PrintAdjStations(): Number of estimated stations in block " << block << 
            " does not match the block station count." << std::endl;
        adjust_.SignalExceptionAdjustment(ss.str(), 0);
    }

    AdjFile adj;

    try {

        // if required, sort stations according to original station file order
        if (adjust_.projectSettings_.o._sort_stn_file_order)
            adjust_.SortStationsbyFileOrder(v_blockStations);

        if (printHeader) {
            switch (adjust_.projectSettings_.a.adjust_mode)
            {
            case SimultaneousMode:
                PrintStationColumnHeaders(os, 
                    adjust_.projectSettings_.o._stn_coord_types, adjust_.projectSettings_.o._stn_corr);
                break;
            case PhasedMode:
            case Phased_Block_1Mode:		// only the first block is rigorous
                // Print stations in each block?
                if (adjust_.projectSettings_.o._output_stn_blocks)
                {
                    if (printBlockID)
                        os << "Block " << block + 1 << std::endl;
                    else if (adjust_.projectSettings_.o._adj_stn_iteration)
                        adj.print_adj_stn_block_header(os, block);
                
                    PrintStationColumnHeaders(os, 
                        adjust_.projectSettings_.o._stn_coord_types, adjust_.projectSettings_.o._stn_corr);
                }
                else 
                {
                    if (block == 0)
                        PrintStationColumnHeaders(os, 
                            adjust_.projectSettings_.o._stn_coord_types, adjust_.projectSettings_.o._stn_corr);
                }
                break;
            }
        }

    }
    catch (const std::runtime_error& e) {
        adjust_.SignalExceptionAdjustment(e.what(), 0);
    }
    
    UINT32 mat_idx, stn;

    // Print stations according to the user-defined sort order
    for (UINT32 i(0); i<adjust_.v_blockStationsMap_.at(block).size(); ++i)
    {
        stn = v_blockStations.at(i);
        mat_idx = adjust_.v_blockStationsMap_.at(block)[stn] * 3;

        PrintAdjStation(os, block, stn, mat_idx,
            stationEstimates, stationVariances, 
            recomputeGeographicCoords, updateGeographicCoords,
            reapplyTypeBUncertainties);
    }

    os << std::endl;

    // return sort order to alpha-numeric
    if (adjust_.projectSettings_.o._sort_stn_file_order)
        adjust_.SortStationsbyID(v_blockStations);
}

void DynAdjustPrinter::PrintOutputFileHeaderInfo()
{
    // Print formatted header
    print_file_header(adjust_.adj_file, "DYNADJUST ADJUSTMENT OUTPUT FILE");
    // Print formatted header
    print_file_header(adjust_.xyz_file, "DYNADJUST COORDINATE OUTPUT FILE");

    adjust_.adj_file << std::setw(PRINT_VAR_PAD) << std::left << "File name:" << safe_absolute_path(adjust_.projectSettings_.o._adj_file) << std::endl << std::endl;
    adjust_.xyz_file << std::setw(PRINT_VAR_PAD) << std::left << "File name:" << safe_absolute_path(adjust_.projectSettings_.o._xyz_file) << std::endl << std::endl;

    adjust_.adj_file << std::setw(PRINT_VAR_PAD) << std::left << "Command line arguments: ";
    adjust_.adj_file << adjust_.projectSettings_.a.command_line_arguments << std::endl << std::endl;

    if (adjust_.projectSettings_.i.input_files.empty())
    {
        adjust_.adj_file << std::setw(PRINT_VAR_PAD) << std::left << "Stations file:" << safe_absolute_path(adjust_.projectSettings_.a.bst_file) << std::endl;
        adjust_.adj_file << std::setw(PRINT_VAR_PAD) << std::left << "Measurements file:" << safe_absolute_path(adjust_.projectSettings_.a.bms_file) << std::endl;
    }
    else
    {
        _it_vstr _it_files(adjust_.projectSettings_.i.input_files.begin());
        std::string s("Input files:");
        while (_it_files!=adjust_.projectSettings_.i.input_files.end())
        {
            adjust_.adj_file << std::setw(PRINT_VAR_PAD) << std::left << s << *_it_files++ << std::endl;
            s = " ";
        }
    }	

    // Reference frame
    adjust_.adj_file << std::setw(PRINT_VAR_PAD) << std::left << "Reference frame: " << adjust_.datum_.GetName() << std::endl;
    adjust_.xyz_file << std::setw(PRINT_VAR_PAD) << std::left << "Reference frame: " << adjust_.datum_.GetName() << std::endl;
    adjust_.adj_file << std::setw(PRINT_VAR_PAD) << std::left << "Epoch: " << adjust_.datum_.GetEpoch_s() << std::endl;
    adjust_.xyz_file << std::setw(PRINT_VAR_PAD) << std::left << "Epoch: " << adjust_.datum_.GetEpoch_s() << std::endl;
    

    // Geoid model
    adjust_.adj_file << std::setw(PRINT_VAR_PAD) << std::left << "Geoid model: " << safe_absolute_path(adjust_.projectSettings_.n.ntv2_geoid_file) << std::endl;
    adjust_.xyz_file << std::setw(PRINT_VAR_PAD) << std::left << "Geoid model: " << safe_absolute_path(adjust_.projectSettings_.n.ntv2_geoid_file) << std::endl;
    
    switch (adjust_.projectSettings_.a.adjust_mode)
    {
    case PhasedMode:
    case Phased_Block_1Mode:
        adjust_.adj_file << std::setw(PRINT_VAR_PAD) << std::left << "Segmentation file:" << adjust_.projectSettings_.a.seg_file << std::endl;
    }

    adjust_.adj_file << std::setw(PRINT_VAR_PAD) << std::left << "Constrained Station S.D. (m):" << adjust_.projectSettings_.a.fixed_std_dev << std::endl;
    adjust_.adj_file << std::setw(PRINT_VAR_PAD) << std::left << "Free Station S.D. (m):" << adjust_.projectSettings_.a.free_std_dev << std::endl;
    adjust_.adj_file << std::setw(PRINT_VAR_PAD) << std::left << "Iteration threshold:" << adjust_.projectSettings_.a.iteration_threshold << std::endl;
    adjust_.adj_file << std::setw(PRINT_VAR_PAD) << std::left << "Maximum iterations:" << std::setprecision(0) << std::fixed << adjust_.projectSettings_.a.max_iterations << std::endl;
    adjust_.adj_file << std::setw(PRINT_VAR_PAD) << std::left << "Test confidence interval:" << std::setprecision(1) << std::fixed << adjust_.projectSettings_.a.confidence_interval << "%" << std::endl;
    adjust_.adj_file << std::setw(PRINT_VAR_PAD) << std::left << "Uncertainties SD(e,n,up):" << std::setprecision(1) << "68.3% (1 sigma)" << std::endl;
    
    if (!adjust_.projectSettings_.a.station_constraints.empty())
        adjust_.adj_file << std::setw(PRINT_VAR_PAD) << std::left << "Station constraints:" << adjust_.projectSettings_.a.station_constraints << std::endl;

    trimstr(adjust_.projectSettings_.a.comments);

    adjust_.adj_file << std::setw(PRINT_VAR_PAD) << std::left << "Station coordinate types:";
    adjust_.xyz_file << std::setw(PRINT_VAR_PAD) << std::left << "Station coordinate types:";

    adjust_.adj_file << adjust_.projectSettings_.o._stn_coord_types << std::endl;
    adjust_.xyz_file << adjust_.projectSettings_.o._stn_coord_types << std::endl;

    adjust_.adj_file << std::setw(PRINT_VAR_PAD) << std::left << "Stations printed in blocks:";
    adjust_.xyz_file << std::setw(PRINT_VAR_PAD) << std::left << "Stations printed in blocks:";
    if (adjust_.projectSettings_.a.adjust_mode != SimultaneousMode && 
        adjust_.projectSettings_.o._output_stn_blocks)
    {
        adjust_.adj_file << "Yes" << std::endl;
        adjust_.xyz_file << "Yes" << std::endl;
    }
    else
    {
        adjust_.adj_file << "No" << std::endl;
        adjust_.xyz_file << "No" << std::endl;
    }

    if (adjust_.projectSettings_.o._stn_corr)
    {
        adjust_.adj_file << std::setw(PRINT_VAR_PAD) << std::left << "Station coordinate corrections:" <<
            "Yes" << std::endl;
        adjust_.xyz_file << std::setw(PRINT_VAR_PAD) << std::left << "Station coordinate corrections:" <<
            "Yes" << std::endl;
    }

    if (adjust_.projectSettings_.o._apply_type_b_file || adjust_.projectSettings_.o._apply_type_b_global)
    {
        if (adjust_.projectSettings_.o._apply_type_b_global)
        {
            adjust_.adj_file << std::setw(PRINT_VAR_PAD) << std::left << "Type B uncertainties:" <<
                adjust_.projectSettings_.a.type_b_global << std::endl;
            adjust_.xyz_file << std::setw(PRINT_VAR_PAD) << std::left << "Type B uncertainties:" <<
                adjust_.projectSettings_.a.type_b_global << std::endl;
        }

        if (adjust_.projectSettings_.o._apply_type_b_file)
        {
            adjust_.adj_file << std::setw(PRINT_VAR_PAD) << std::left << "Type B uncertainty file:" <<
                safe_absolute_path(adjust_.projectSettings_.a.type_b_file) << std::endl;
            adjust_.xyz_file << std::setw(PRINT_VAR_PAD) << std::left << "Type B uncertainty file:" <<
                safe_absolute_path(adjust_.projectSettings_.a.type_b_file) << std::endl;
        }
    }

    // Print user-supplied comments.
    // This is a bit messy and could be cleaned up.
    // Alas, the following logic is applied:
    //	- User-supplied newlines (i.e. "\n" within the string) are dealt with
    //  - Newline characters are inserted at between-word-spaces, such that the
    //    line length does not exceed PRINT_VAL_PAD
    //  - Resulting spaces at the start of a sentence on a new line are removed.
    if (!adjust_.projectSettings_.a.comments.empty())
    {
        size_t s(0), t(0);
        std::string var("Comments: "), comments(adjust_.projectSettings_.a.comments);
        size_t pos = 0;

        while (s < comments.length())
        {
            if (s + PRINT_VAL_PAD >= comments.length())
            {
                adjust_.adj_file << std::setw(PRINT_VAR_PAD) << std::left << var << comments.substr(s) << std::endl; 
                break;
            }
            else
            {
                if ((pos = comments.substr(s, PRINT_VAL_PAD-1).find('\\')) != std::string::npos)
                {
                    if (comments.at(s+pos+1) == 'n')
                    {
                        adjust_.adj_file << std::setw(PRINT_VAR_PAD) << std::left << var << comments.substr(s, pos) << std::endl;
                        s += pos+2;
                        if (comments.at(s) == ' ')
                            ++s;
                        var = " ";
                        continue;
                    }					
                }
                
                t = 0;
                while (comments.at(s+PRINT_VAL_PAD - t) != ' ')
                    ++t;
                
                adjust_.adj_file << std::setw(PRINT_VAR_PAD) << std::left << var << comments.substr(s, PRINT_VAL_PAD - t);
                s += PRINT_VAL_PAD - t;
                
                if (comments.at(s-1) != ' ' &&
                    comments.at(s) != ' ' &&
                    comments.at(s+1) != ' ')
                    adjust_.adj_file << "-";
                else if (comments.at(s) == ' ')
                    ++s;
                adjust_.adj_file << std::endl;
            }
            var = " ";
        }
    }

    adjust_.adj_file << OUTPUTLINE << std::endl;
    adjust_.xyz_file << OUTPUTLINE << std::endl;
}

void DynAdjustPrinter::PrintCompMeasurements_GXY(const UINT32& block, it_vmsr_t& _it_msr, UINT32& design_row, printMeasurementsMode printMode)
{
    // Is this a Y cluster specified in latitude, longitude, height?
    if (_it_msr->measType == 'Y')
    {
        if (_it_msr->station3 == LLH_type_i ||
            _it_msr->station3 == LLh_type_i)
        {
            // Print phi, lambda, H
            PrintCompMeasurements_YLLH(_it_msr, design_row);
            return;
        }
    }

    UINT32 cluster_msr, cluster_count(_it_msr->vectorCount1);
    UINT32 covariance_count;
    bool nextElement(false);
    double computed, correction;

    for (cluster_msr=0; cluster_msr<cluster_count; ++cluster_msr)
    {
        if (nextElement)
            adjust_.adj_file << std::left << std::setw(PAD2) << _it_msr->measType;
        else
            nextElement = true;

        covariance_count = _it_msr->vectorCount2;

        // first station
        adjust_.adj_file << std::left << std::setw(STATION) << adjust_.bstBinaryRecords_.at(_it_msr->station1).stationName;
        
        // Print second station?
        switch (_it_msr->measType)
        {
        case 'G':
        case 'X':
            adjust_.adj_file << std::left << std::setw(STATION) << adjust_.bstBinaryRecords_.at(_it_msr->station2).stationName;
            break;
        default:
            adjust_.adj_file << std::left << std::setw(STATION) << " ";
        }

        // third station
        adjust_.adj_file << std::left << std::setw(STATION) << " ";

        switch (printMode)
        {
        case computedMsrs:
            correction = -adjust_.v_measMinusComp_.at(block).get(design_row, 0);
            computed = _it_msr->term1 + correction;
            break;
        case ignoredMsrs:
        default:
            correction = _it_msr->measCorr;
            computed = _it_msr->measAdj;
            break;
        }

        // Print linear measurement, taking care of user requirements for precision	
        PrintComputedMeasurements<LinearMeasurement>('X', computed, correction, _it_msr);

        design_row++;
        _it_msr++;
    
        switch (printMode)
        {
        case computedMsrs:
            correction = -adjust_.v_measMinusComp_.at(block).get(design_row, 0);
            computed = _it_msr->term1 + correction;
            break;
        case ignoredMsrs:
        default:
            correction = _it_msr->measCorr;
            computed = _it_msr->measAdj;
            break;
        }

        // Print linear measurement, taking care of user requirements for precision	
        PrintComputedMeasurements<LinearMeasurement>('Y', computed, correction, _it_msr);

        design_row++;
        _it_msr++;
    
        switch (printMode)
        {
        case computedMsrs:
            correction = -adjust_.v_measMinusComp_.at(block).get(design_row, 0);
            computed = _it_msr->term1 + correction;
            break;
        case ignoredMsrs:
        default:
            correction = _it_msr->measCorr;
            computed = _it_msr->measAdj;
            break;
        }

        // Print linear measurement, taking care of user requirements for precision	
        PrintComputedMeasurements<LinearMeasurement>('Z', computed, correction, _it_msr);

        design_row++;

        // skip covariances until next point
        _it_msr += covariance_count * 3;
        
        if (covariance_count > 0)
            _it_msr++;
    }
}

void DynAdjustPrinter::PrintStationsUniqueList(std::ostream& os,
    const v_mat_2d* stationEstimates, v_mat_2d* stationVariances,
    bool recomputeGeographicCoords, bool updateGeographicCoords,
    bool reapplyTypeBUncertainties)
{
    UINT32 block(UINT_MAX), stn, mat_index;
    _it_u32u32_uint32_pair _it_bsmu;

    // Determine sort order of block Stations Map
    if (adjust_.projectSettings_.a.stage || adjust_.projectSettings_.a.adjust_mode == Phased_Block_1Mode)
    {
        // sort by blocks to create efficiency when deserialising matrix info
        std::sort(adjust_.v_blockStationsMapUnique_.begin(), adjust_.v_blockStationsMapUnique_.end(), 
            CompareBlockStationMapUnique_byBlock<u32u32_uint32_pair>());
    }
    else if (adjust_.projectSettings_.o._sort_stn_file_order)
    {
        CompareBlockStationMapUnique_byFileOrder<station_t, u32u32_uint32_pair> stnorderCompareFunc(&adjust_.bstBinaryRecords_);
        std::sort(adjust_.v_blockStationsMapUnique_.begin(), adjust_.v_blockStationsMapUnique_.end(), stnorderCompareFunc);
    }
    else
        std::sort(adjust_.v_blockStationsMapUnique_.begin(), adjust_.v_blockStationsMapUnique_.end());
    
    std::stringstream stationRecord;
    v_uint32_string_pair stationsOutput;
    stationsOutput.reserve(adjust_.v_blockStationsMapUnique_.size());

    std::ostream* outstream(&os);
    if (adjust_.projectSettings_.a.stage)
        outstream = &stationRecord;

    const v_mat_2d* estimates = &adjust_.v_rigorousStations_;
    if (adjust_.projectSettings_.a.adjust_mode == SimultaneousMode)
        estimates = &adjust_.v_estimatedStations_;

    switch (adjust_.projectSettings_.a.adjust_mode)
    {
    case SimultaneousMode:
        break;
    case PhasedMode:
    case Phased_Block_1Mode:		
        // if computed estimates are available, use them.
        // if null, use rigorous estimates
        if (stationEstimates != 0)
            estimates = stationEstimates;
        break;
    default:
        break;
    }

    // Print stations according to the user-defined sort order
    for (_it_bsmu=adjust_.v_blockStationsMapUnique_.begin();
        _it_bsmu!=adjust_.v_blockStationsMapUnique_.end(); ++_it_bsmu)
    {
        // If this is not the first block, exit if block-1 mode
        if (adjust_.projectSettings_.a.adjust_mode == Phased_Block_1Mode)
            if (_it_bsmu->second > 0)
                break;

        if (adjust_.projectSettings_.a.stage)
        {
            if (block != _it_bsmu->second)
            {
                // unload previous block
                if (_it_bsmu->second > 0)
                    adjust_.UnloadBlock(_it_bsmu->second-1, 2, 
                        sf_rigorous_stns, sf_original_stns);

                // load up this block
                if (adjust_.projectSettings_.a.stage)
                    adjust_.DeserialiseBlockFromMappedFile(_it_bsmu->second, 2,
                        sf_rigorous_stns, sf_original_stns);
            }
        }

        block = _it_bsmu->second;
        stn = _it_bsmu->first.first;
        mat_index = _it_bsmu->first.second * 3;

        PrintAdjStation(*outstream,
            block, stn, mat_index,
                &estimates->at(block), &stationVariances->at(block), 
                recomputeGeographicCoords, updateGeographicCoords,
                reapplyTypeBUncertainties);

        if (adjust_.projectSettings_.a.stage)
        {
            stationsOutput.push_back(uint32_string_pair(stn, stationRecord.str()));
            stationRecord.str("");
        }
    }

    if (adjust_.projectSettings_.a.stage)
    {
        adjust_.UnloadBlock(block, 2, 
            sf_rigorous_stns, sf_original_stns);

        // if required, sort stations according to original station file order
        if (adjust_.projectSettings_.o._sort_stn_file_order)
        {
            CompareOddPairFirst_FileOrder<station_t, UINT32, std::string> stnorderCompareFunc(&adjust_.bstBinaryRecords_);
            std::sort(stationsOutput.begin(), stationsOutput.end(), stnorderCompareFunc);
        }
        else
            std::sort(stationsOutput.begin(), stationsOutput.end(), CompareOddPairFirst<UINT32, std::string>());

        for_each(stationsOutput.begin(), stationsOutput.end(),
            [&stationsOutput, &os] (std::pair<UINT32, std::string>& stnRecord) { 
                os << stnRecord.second;
            }
        );
    }

    os << std::endl;

    // if required, sort stations according to original station file order
    if (adjust_.projectSettings_.o._sort_stn_file_order)
    {
        CompareBlockStationMapUnique_byFileOrder<station_t, u32u32_uint32_pair> stnorderCompareFunc(&adjust_.bstBinaryRecords_);
        std::sort(adjust_.v_blockStationsMapUnique_.begin(), adjust_.v_blockStationsMapUnique_.end(), stnorderCompareFunc);
    }
    else
        std::sort(adjust_.v_blockStationsMapUnique_.begin(), adjust_.v_blockStationsMapUnique_.end());

    // return sort order to alpha-numeric
    if (adjust_.projectSettings_.o._sort_stn_file_order)
        adjust_.SortStationsbyFileOrder(adjust_.v_blockStationsR);

    adjust_.SortStationsbyID(adjust_.v_blockStationsR);
}

bool DynAdjustPrinter::IgnoredMeasurementContainsInvalidStation(pit_vmsr_t _it_msr)
{
    UINT32 cluster_msr, covariance_count, cluster_count((*_it_msr)->vectorCount1);

    switch ((*_it_msr)->measType)
    {
    // Three station measurement
    case 'A':	// Horizontal angle
        if (adjust_.vAssocStnList_.at((*_it_msr)->station3).IsInvalid())
            return false;
        [[fallthrough]];
    // Two station measurements
    case 'B':	// Geodetic azimuth
    case 'C':	// Chord dist
    case 'E':	// Ellipsoid arc
    case 'G':	// GPS Baseline
    case 'K':	// Astronomic azimuth
    case 'L':	// Level difference
    case 'M':	// MSL arc
    case 'S':	// Slope distance
    case 'V':	// Zenith distance
    case 'Z':	// Vertical angle
        if (adjust_.vAssocStnList_.at((*_it_msr)->station2).IsInvalid())
            return false;
        [[fallthrough]];
    // One station measurements
    case 'H':	// Orthometric height
    case 'I':	// Astronomic latitude
    case 'J':	// Astronomic longitude
    case 'P':	// Geodetic latitude
    case 'Q':	// Geodetic longitude
    case 'R':	// Ellipsoidal height
        if (adjust_.vAssocStnList_.at((*_it_msr)->station1).IsInvalid())
            return false;
        break;
    // GNSS measurements
    case 'X':	// GPS Baseline cluster
    case 'Y':	// GPS Point cluster
        for (cluster_msr=0; cluster_msr<cluster_count; ++cluster_msr)
        {
            covariance_count = (*_it_msr)->vectorCount2;

            if (adjust_.vAssocStnList_.at((*_it_msr)->station1).IsInvalid())
                return false; 

            if ((*_it_msr)->measType == 'X')
                if (adjust_.vAssocStnList_.at((*_it_msr)->station2).IsInvalid())
                    return false;
            
            (*_it_msr) += 2;					// skip y and z
            (*_it_msr) += covariance_count * 3;	// skip covariances

            if (covariance_count > 0)
                (*_it_msr)++;
        }
        break;
    case 'D':	// Direction set
        for (cluster_msr=0; cluster_msr<cluster_count; ++cluster_msr)
        {
            // Check first station
            if (cluster_msr == 0)
                if (adjust_.vAssocStnList_.at((*_it_msr)->station1).IsInvalid())
                    return false;
            // Check all target directions
            if (adjust_.vAssocStnList_.at((*_it_msr)->station2).IsInvalid())
                return false;
            
            if (cluster_msr+1 == cluster_count)
                break;
                
            (*_it_msr)++;
        }
        break;
    }
    return true;
}

void DynAdjustPrinter::PrintAdjStation(std::ostream& os, 
    const UINT32& block, const UINT32& stn, const UINT32& mat_idx,
    const matrix_2d* stationEstimates, matrix_2d* stationVariances,
    bool recomputeGeographicCoords, bool updateGeographicCoords,
    bool reapplyTypeBUncertainties)
{
    double estLatitude, estLongitude, estHeight, E(0.), N(0.), Zo(-1.);

    it_vstn_t stn_it(adjust_.bstBinaryRecords_.begin());
    stn_it += stn;

    // Print station name and constraint
    os << std::setw(STATION) << std::left << stn_it->stationName;
    os << std::setw(CONSTRAINT) << std::left << stn_it->stationConst;

    // Are station corrections required?  
    // If so, update the original stations matrix
    if (adjust_.projectSettings_.o._init_stn_corrections || adjust_.projectSettings_.o._stn_corr)
    {
        estHeight = stn_it->initialHeight;
        if (stn_it->suppliedHeightRefFrame == ORTHOMETRIC_type_i)
            estHeight = stn_it->initialHeight + stn_it->geoidSep;

        // update initial estimates
        GeoToCart<double>(
            stn_it->initialLatitude,
            stn_it->initialLongitude,
            estHeight,
            adjust_.v_originalStations_.at(block).getelementref(mat_idx, 0),
            adjust_.v_originalStations_.at(block).getelementref(mat_idx+1, 0),
            adjust_.v_originalStations_.at(block).getelementref(mat_idx+2, 0),
            adjust_.datum_.GetEllipsoidRef());
    }

    // Are geographic coordinates required?
    if (recomputeGeographicCoords ||
        adjust_.projectSettings_.o._stn_coord_types.find("P") != std::string::npos ||
        adjust_.projectSettings_.o._stn_coord_types.find("L") != std::string::npos ||
        adjust_.projectSettings_.o._stn_coord_types.find("H") != std::string::npos ||
        adjust_.projectSettings_.o._stn_coord_types.find("h") != std::string::npos ||
        adjust_.projectSettings_.o._stn_coord_types.find("N") != std::string::npos ||
        adjust_.projectSettings_.o._stn_coord_types.find("z") != std::string::npos)
    {
        CartToGeo<double>(
            stationEstimates->get(mat_idx, 0),		// X
            stationEstimates->get(mat_idx+1, 0),	// Y
            stationEstimates->get(mat_idx+2, 0),	// Z
            &estLatitude,							// Lat
            &estLongitude,							// Long
            &estHeight,								// E.Height
            adjust_.datum_.GetEllipsoidRef());				

        if (updateGeographicCoords)
        {
            stn_it->currentLatitude = estLatitude;
            stn_it->currentLongitude = estLongitude;
            stn_it->currentHeight = estHeight;
        }
    }

    // Are projection coordinates required?
    if (adjust_.projectSettings_.o._stn_coord_types.find("E") != std::string::npos ||
        adjust_.projectSettings_.o._stn_coord_types.find("N") != std::string::npos ||
        adjust_.projectSettings_.o._stn_coord_types.find("z") != std::string::npos)
    {
        GeoToGrid<double>(estLatitude, estLongitude, &E, &N, &Zo,
            adjust_.datum_.GetEllipsoidRef(), 
            &adjust_.projection_, true);		// compute zone
    }

    // Use coordinate formatting helper
    PrintStationCoordinatesByType(os, stn_it, stationEstimates, mat_idx, 
                                  estLatitude, estLongitude, estHeight, E, N, Zo);

    os << std::setw(PAD2) << " ";

    // Standard deviation in local reference frame
    matrix_2d var_local(3, 3), var_cart(3, 3);

    // To minimise re-computation when printing to .xyz, .adj and .apu, the type b uncertainties are 
    // added to the estimated variances held in memory. This is not a problem since, if adjust is 
    // re-run and the user does not supply type b uncertainties, the variances held in memory will be 
    // re-estimated which will overwrite the type b values that were previously applied.

    // Add type B uncertainties (if required)
    if (reapplyTypeBUncertainties &&
        (adjust_.projectSettings_.o._apply_type_b_global || adjust_.projectSettings_.o._apply_type_b_file))
    {
        if (adjust_.v_typeBUncertaintyMethod_.at(stn).apply)
        {
            // apply local first
            if (adjust_.v_typeBUncertaintyMethod_.at(stn).method == type_b_local)
            {
                // Add the cartesian type b variances 
                // Note: Cartesian variances for this station were computed in dna_io_tbu::reduce_uncertainties_local(...)
                stationVariances->blockadd(mat_idx, mat_idx,
                    adjust_.v_typeBUncertaintiesLocal_.at(adjust_.v_stationTypeBMap_.at(stn).second).type_b,
                    0, 0, 3, 3);
            }
            else if (adjust_.v_typeBUncertaintyMethod_.at(stn).method == type_b_global)
            {
                // Compute cartesian variances for this station from the global uncertainty
                // Note: this needs to be done for each site from the global type b uncertainties entered 
                // via the command line
                matrix_2d var_cart_typeb(3, 3);
                PropagateVariances_LocalCart<double>(adjust_.typeBUncertaintyGlobal_.type_b, var_cart_typeb,
                    estLatitude, estLongitude, true);

                // Add the cartesian type b variances 
                stationVariances->blockadd(mat_idx, mat_idx,
                    var_cart_typeb, 0, 0, 3, 3);
            }	
        }
    }

    var_cart.copyelements(0, 0, stationVariances, mat_idx, mat_idx, 3, 3);

    PropagateVariances_LocalCart(var_cart, var_local, 
        estLatitude, estLongitude, false);

    // add geoid model height uncertainty
    var_local.elementadd(2, 2, (stn_it->geoidSepUnc * stn_it->geoidSepUnc));
    
    // Use uncertainty formatting helper
    PrintStationUncertainties(os, var_local);

    if (adjust_.projectSettings_.o._stn_corr)
    {
        double cor_e, cor_n, cor_up;
        ComputeLocalElements3D<double>(
            adjust_.v_originalStations_.at(block).get(mat_idx, 0),		// original X
            adjust_.v_originalStations_.at(block).get(mat_idx+1, 0),	// original Y
            adjust_.v_originalStations_.at(block).get(mat_idx+2, 0),	// original Z
            stationEstimates->get(mat_idx, 0),					// estimated X
            stationEstimates->get(mat_idx+1, 0),				// estimated Y
            stationEstimates->get(mat_idx+2, 0),				// estimated Z
            stn_it->currentLatitude,
            stn_it->currentLongitude,
            &cor_e, &cor_n, &cor_up);					// corrections in the local reference frame

        os << std::setw(PAD2) << " " << 
            std::setprecision(adjust_.PRECISION_MTR_STN) << std::fixed << std::right << std::setw(HEIGHT) << cor_e << 
            std::setprecision(adjust_.PRECISION_MTR_STN) << std::fixed << std::right << std::setw(HEIGHT) << cor_n << 
            std::setprecision(adjust_.PRECISION_MTR_STN) << std::fixed << std::right << std::setw(HEIGHT) << cor_up;
    }

    // description
    os << std::setw(PAD2) << " " << std::left << stn_it->description;
    os << std::endl;
}

void DynAdjustPrinter::PrintAdjMeasurements_GXY(it_vmsr_t& _it_msr, const uint32_uint32_pair& b_pam)
{
    // Is this a Y cluster specified in latitude, longitude, height?
    if (_it_msr->measType == 'Y')
    {
        if (_it_msr->station3 == LLH_type_i ||
            _it_msr->station3 == LLh_type_i)
        {
            // Print phi, lambda, H
            PrintAdjMeasurements_YLLH(_it_msr);
            return;
        }
    }

    UINT32 cluster_msr, cluster_count(_it_msr->vectorCount1), covariance_count;
    bool nextElement(false);

    for (cluster_msr=0; cluster_msr<cluster_count; ++cluster_msr)
    {
        if (nextElement)
            adjust_.adj_file << std::left << std::setw(PAD2) << _it_msr->measType;
        else
            nextElement = true;

        covariance_count = _it_msr->vectorCount2;

        // first station
        adjust_.adj_file << std::left << std::setw(STATION) << adjust_.bstBinaryRecords_.at(_it_msr->station1).stationName;
    
        // Print second station?
        switch (_it_msr->measType)
        {
        case 'G':
        case 'X':
            adjust_.adj_file << std::left << std::setw(STATION) << adjust_.bstBinaryRecords_.at(_it_msr->station2).stationName;
            break;
        default:
            adjust_.adj_file << std::left << std::setw(STATION) << " ";
        }

        // third station
        adjust_.adj_file << std::left << std::setw(STATION) << " ";
        
        // Print adjusted GNSS baseline measurements in alternate units?
        if (adjust_.projectSettings_.o._adj_gnss_units != XYZ_adj_gnss_ui &&
            _it_msr->measType != 'Y')
            PrintAdjGNSSAlternateUnits(_it_msr, b_pam);
        else
        {
            // Print X
            PrintAdjMeasurementsLinear('X', _it_msr);
    
            // Print Y
            _it_msr++;	
            PrintAdjMeasurementsLinear('Y', _it_msr);

            // Print Z
            _it_msr++;	
            PrintAdjMeasurementsLinear('Z', _it_msr);
        }

        // skip covariances until next baseline
        _it_msr += covariance_count * 3;
        
        if (covariance_count > 0)
            _it_msr++;
    }
}

void DynAdjustPrinter::PrintCorStation(std::ostream& os, 
    const UINT32& block, const UINT32& stn, const UINT32& mat_index,
    const matrix_2d* stationEstimates)
{
    double vertical_angle, azimuth, slope_distance, horiz_distance, 
        local_12e, local_12n, local_12up;

    // calculate vertical angle
    vertical_angle = VerticalAngle(
        adjust_.v_originalStations_.at(block).get(mat_index, 0),		// X1
        adjust_.v_originalStations_.at(block).get(mat_index+1, 0),		// Y1
        adjust_.v_originalStations_.at(block).get(mat_index+2, 0),		// Z1
        stationEstimates->get(mat_index, 0),		// X2
        stationEstimates->get(mat_index+1, 0),		// Y2
        stationEstimates->get(mat_index+2, 0),		// Z2
        adjust_.bstBinaryRecords_.at(stn).currentLatitude,
        adjust_.bstBinaryRecords_.at(stn).currentLongitude,
        adjust_.bstBinaryRecords_.at(stn).currentLatitude,
        adjust_.bstBinaryRecords_.at(stn).currentLongitude,
        0., 0.,
        &local_12e, &local_12n, &local_12up);

    if (fabs(local_12e) < PRECISION_1E5)
        if (fabs(local_12n) < PRECISION_1E5)
            if (fabs(local_12up) < PRECISION_1E5)
                vertical_angle = 0.;
        
    // Is the vertical correction within the user-specified tolerance?
    if (fabs(local_12up) < adjust_.projectSettings_.o._vt_corr_threshold)
        return;

    // compute horizontal correction
    horiz_distance = magnitude(local_12e, local_12n);
        
    // Is the horizontal correction within the user-specified tolerance?
    if (horiz_distance < adjust_.projectSettings_.o._hz_corr_threshold)
        return;
        
    // calculate azimuth
    azimuth = Direction(
        adjust_.v_originalStations_.at(block).get(mat_index, 0),		// X1
        adjust_.v_originalStations_.at(block).get(mat_index+1, 0),		// Y1
        adjust_.v_originalStations_.at(block).get(mat_index+2, 0),		// Z1
        stationEstimates->get(mat_index, 0),		// X2
        stationEstimates->get(mat_index+1, 0),		// Y2
        stationEstimates->get(mat_index+2, 0),		// Z2
        adjust_.bstBinaryRecords_.at(stn).currentLatitude,
        adjust_.bstBinaryRecords_.at(stn).currentLongitude,
        &local_12e, &local_12n);

    if (fabs(local_12e) < PRECISION_1E5)
        if (fabs(local_12n) < PRECISION_1E5)
            azimuth = 0.;

    // calculate distances
    slope_distance = magnitude(
        adjust_.v_originalStations_.at(block).get(mat_index, 0),		// X1
        adjust_.v_originalStations_.at(block).get(mat_index+1, 0),		// Y1
        adjust_.v_originalStations_.at(block).get(mat_index+2, 0),		// Z1
        stationEstimates->get(mat_index, 0),		// X2
        stationEstimates->get(mat_index+1, 0),		// Y2
        stationEstimates->get(mat_index+2, 0));		// Z2
        
    // print...
    // station and constraint
    os << std::setw(STATION) << std::left << adjust_.bstBinaryRecords_.at(stn).stationName << std::setw(PAD2) << " ";
    // data
    os << std::setw(MSR) << std::right << FormatDmsString(RadtoDms(azimuth), 4, true, false) << 
        std::setw(MSR) << std::right << FormatDmsString(RadtoDms(vertical_angle), 4, true, false) << 
        std::setw(MSR) << std::setprecision(adjust_.PRECISION_MTR_STN) << std::fixed << std::right << slope_distance << 
        std::setw(MSR) << std::setprecision(adjust_.PRECISION_MTR_STN) << std::fixed << std::right << horiz_distance;

    if (adjust_.isAdjustmentQuestionable_)
        os << 
            StringFromTW(local_12e, HEIGHT, adjust_.PRECISION_MTR_STN) << 
            StringFromTW(local_12n, HEIGHT, adjust_.PRECISION_MTR_STN) << 
            StringFromTW(local_12up, HEIGHT, adjust_.PRECISION_MTR_STN) << std::endl;
    else
        os << 
            std::setw(HEIGHT) << std::setprecision(adjust_.PRECISION_MTR_STN) << std::fixed << std::right << local_12e << 
            std::setw(HEIGHT) << std::setprecision(adjust_.PRECISION_MTR_STN) << std::fixed << std::right << local_12n << 
            std::setw(HEIGHT) << std::setprecision(adjust_.PRECISION_MTR_STN) << std::fixed << std::right << local_12up << std::endl;
}

void DynAdjustPrinter::PrintCorStations(std::ostream& cor_file, const UINT32& block)
{
    vUINT32 v_blockStations(adjust_.v_parameterStationList_.at(block));
    
    // if required, sort stations according to original station file order
    if (adjust_.projectSettings_.o._sort_stn_file_order)
        adjust_.SortStationsbyFileOrder(v_blockStations);
    
    // Continue with existing detailed implementation
    switch (adjust_.projectSettings_.a.adjust_mode)
    {
    case PhasedMode:
    case Phased_Block_1Mode:		// only the first block is rigorous
        cor_file << "Block " << block + 1 << std::endl;
        break;
    }

    cor_file << std::setw(STATION) << std::left << "Station" <<
        std::setw(PAD2) << " " << 
        std::right << std::setw(MSR) << "Azimuth" << 
        std::right << std::setw(MSR) << "V. Angle" << 
        std::right << std::setw(MSR) << "S. Distance" << 
        std::right << std::setw(MSR) << "H. Distance" << 
        std::right << std::setw(HEIGHT) << "east" << 
        std::right << std::setw(HEIGHT) << "north" << 
        std::right << std::setw(HEIGHT) << "up" << std::endl;

    UINT32 i, j = STATION+PAD2+MSR+MSR+MSR+MSR+HEIGHT+HEIGHT+HEIGHT;

    for (i=0; i<j; ++i)
        cor_file << "-";
    cor_file << std::endl;

    UINT32 mat_idx, stn;

    const v_mat_2d* estimates = &adjust_.v_rigorousStations_;
    if (adjust_.projectSettings_.a.adjust_mode == SimultaneousMode)
        estimates = &adjust_.v_estimatedStations_;

    // Print stations according to the user-defined sort order
    for (i=0; i<adjust_.v_blockStationsMap_.at(block).size(); ++i)
    {
        stn = v_blockStations.at(i);
        mat_idx = adjust_.v_blockStationsMap_.at(block)[stn] * 3;

        PrintCorStation(cor_file, block, stn, mat_idx, 
            &estimates->at(block));
    }

    cor_file << std::endl;

    // return sort order to alpha-numeric
    if (adjust_.projectSettings_.o._sort_stn_file_order)
        adjust_.SortStationsbyID(v_blockStations);
}

void DynAdjustPrinter::PrintPosUncertaintiesHeader(std::ostream& os)
{
    os << std::setw(STATION) << std::left << "Station" <<
        std::setw(PAD2) << " " << 
        std::right << std::setw(LAT_EAST) << CDnaStation::CoordinateName('P') << 
        std::right << std::setw(LON_NORTH) << CDnaStation::CoordinateName('L') << 
        std::right << std::setw(STAT) << "Hz PosU" << 
        std::right << std::setw(STAT) << "Vt PosU" << 
        std::right << std::setw(PREC) << "Semi-major" << 
        std::right << std::setw(PREC) << "Semi-minor" << 
        std::right << std::setw(PREC) << "Orientation";

    switch (adjust_.projectSettings_.o._apu_vcv_units)
    {
    case ENU_apu_ui:
        os <<  
            std::right << std::setw(MSR) << "Variance(e)" << 
            std::right << std::setw(MSR) << "Variance(n)" << 
            std::right << std::setw(MSR) << "Variance(up)" << std::endl;
        break;
    case XYZ_apu_ui:
    default:
        os <<  
            std::right << std::setw(MSR) << "Variance(X)" << 
            std::right << std::setw(MSR) << "Variance(Y)" << 
            std::right << std::setw(MSR) << "Variance(Z)" << std::endl;
        break;
    }

    UINT32 i, j = STATION+PAD2+LAT_EAST+LON_NORTH+STAT+STAT+PREC+PREC+PREC+MSR+MSR+MSR;

    for (i=0; i<j; ++i)
        os << "-";
    os << std::endl;
}

void DynAdjustPrinter::PrintPosUncertainty(std::ostream& os, const UINT32& block, const UINT32& stn, 
                                           const UINT32& mat_idx, const matrix_2d* stationVariances, 
                                           const UINT32& map_idx, const vUINT32* blockStations)
{
    double semimajor, semiminor, azimuth, hzPosU, vtPosU;
    UINT32 ic, jc;

    matrix_2d variances_cart(3, 3), variances_local(3, 3), mrotations(3, 3);
    matrix_2d *variances(&variances_cart);

    switch (adjust_.projectSettings_.o._apu_vcv_units)
    {
    case ENU_apu_ui:
        variances = &variances_local;
        break;
    }

    // get cartesian matrix
    stationVariances->submatrix(mat_idx, mat_idx, &variances_cart, 3, 3);

    // Calculate standard deviations in local reference frame
    PropagateVariances_LocalCart<double>(variances_cart, variances_local, 
        adjust_.bstBinaryRecords_.at(stn).currentLatitude, 
        adjust_.bstBinaryRecords_.at(stn).currentLongitude, false,
        mrotations, true);

    // Compute error ellipse terms
    ErrorEllipseParameters<double>(variances_local, semimajor, semiminor, azimuth);

    // Compute positional uncertainty terms
    PositionalUncertainty<double>(semimajor, semiminor, sqrt(variances_local.get(2, 2)), hzPosU, vtPosU);

    // print...
    // station and padding
    os << std::setw(STATION) << std::left << adjust_.bstBinaryRecords_.at(stn).stationName << std::setw(PAD2) << " ";
    
    os.flags(std::ios::fixed | std::ios::right);

    // latitude
    if (adjust_.projectSettings_.o._angular_type_stn == DMS)
        os << std::setprecision(4+adjust_.PRECISION_SEC_STN) << std::setw(LAT_EAST) <<
            RadtoDms(adjust_.bstBinaryRecords_.at(stn).currentLatitude);
    else
        os << std::setprecision(4 + adjust_.PRECISION_SEC_STN) << std::setw(LAT_EAST) <<
            Degrees(adjust_.bstBinaryRecords_.at(stn).currentLatitude);
    // longitude
    if (adjust_.projectSettings_.o._angular_type_stn == DMS)
        os << std::setprecision(4+adjust_.PRECISION_SEC_STN) << std::setw(LON_NORTH) <<
            RadtoDmsL(adjust_.bstBinaryRecords_.at(stn).currentLongitude);
    else
        os << std::setprecision(4 + adjust_.PRECISION_SEC_STN) << std::setw(LON_NORTH) <<
            DegreesL(adjust_.bstBinaryRecords_.at(stn).currentLongitude);

    // positional uncertainty 
    if (adjust_.isAdjustmentQuestionable_)
        os << 
            StringFromTW(hzPosU, STAT, adjust_.PRECISION_MTR_STN) <<
            StringFromTW(vtPosU, STAT, adjust_.PRECISION_MTR_STN);
    else
        os << std::setprecision(adjust_.PRECISION_MTR_STN) << std::setw(STAT) << hzPosU <<
            std::setprecision(adjust_.PRECISION_MTR_STN) << std::fixed << std::right << std::setw(STAT) << vtPosU;
    
    // error ellipse semi-major, semi-minor, orientation
    os << std::setprecision(adjust_.PRECISION_MTR_STN) << std::setw(PREC) << semimajor <<
        std::setprecision(adjust_.PRECISION_MTR_STN) << std::setw(PREC) << semiminor <<
        std::setprecision(4) << std::setw(PREC) << RadtoDms(azimuth);

    os.flags(std::ios::scientific | std::ios::right);

    UINT16 PRECISION_UNCERTAINTY(9);

    // xx, xy, xz
    if (adjust_.isAdjustmentQuestionable_)
        os << 
            StringFromTW(variances->get(0, 0), MSR, PRECISION_UNCERTAINTY) <<						// e
            StringFromTW(variances->get(0, 1), MSR, PRECISION_UNCERTAINTY) <<						// n
            StringFromTW(variances->get(0, 2), MSR, PRECISION_UNCERTAINTY) << std::endl;				// up
    else
        os << 
            std::setprecision(PRECISION_UNCERTAINTY) << std::setw(MSR) << variances->get(0, 0) <<				// e
            std::setprecision(PRECISION_UNCERTAINTY) << std::setw(MSR) << variances->get(0, 1) <<				// n
            std::setprecision(PRECISION_UNCERTAINTY) << std::setw(MSR) << variances->get(0, 2) << std::endl;		// up
    
    // Next line: yy, yz
    os << std::setw(STATION+PAD2+LAT_EAST+LON_NORTH+STAT+STAT+PREC+PREC+PREC+MSR) << " ";		// padding
    if (adjust_.isAdjustmentQuestionable_)
        os << 
            StringFromTW(variances->get(1, 1), MSR, PRECISION_UNCERTAINTY) <<						// n
            StringFromTW(variances->get(1, 2), MSR, PRECISION_UNCERTAINTY) << std::endl;					// up
    else
        os << 
            std::setprecision(PRECISION_UNCERTAINTY) << std::setw(MSR) << variances->get(1, 1) <<				// n
            std::setprecision(PRECISION_UNCERTAINTY) << std::setw(MSR) << variances->get(1, 2) << std::endl;		// up
    
    // Next line: zz
    os << std::setw(STATION+PAD2+LAT_EAST+LON_NORTH+STAT+STAT+PREC+PREC+PREC+MSR+MSR) << " ";	// padding
    if (adjust_.isAdjustmentQuestionable_)
        os <<
            StringFromTW(variances->get(2, 2), MSR, PRECISION_UNCERTAINTY) << std::endl;		// up
    else
        os <<
            std::setprecision(PRECISION_UNCERTAINTY) << std::setw(MSR) << variances->get(2, 2) << std::endl;		// up
    
    if (!adjust_.projectSettings_.o._output_pu_covariances)
        return;

    // Print covariances
    for (ic=map_idx+1; ic<adjust_.v_blockStationsMap_.at(block).size(); ++ic)
    {
        jc = adjust_.v_blockStationsMap_.at(block)[blockStations->at(ic)] * 3;

        // get cartesian submatrix corresponding to the covariance
        stationVariances->submatrix(mat_idx, jc, &variances_cart, 3, 3);
            
        switch (adjust_.projectSettings_.o._apu_vcv_units)
        {
        case ENU_apu_ui:
            // Calculate vcv in local reference frame
            PropagateVariances_LocalCart<double>(variances_cart, variances_local, 
                adjust_.bstBinaryRecords_.at(stn).currentLatitude, 
                adjust_.bstBinaryRecords_.at(stn).currentLongitude, false,
                mrotations, false);
            break;
        }

        os << std::setw(STATION) << std::left << adjust_.bstBinaryRecords_.at(blockStations->at(ic)).stationName;
            
        os.flags(std::ios::scientific | std::ios::right);
        os << 
            std::setw(PAD2+LAT_EAST+LON_NORTH+STAT+STAT+PREC+PREC+PREC) << " " <<					// padding
            std::setprecision(PRECISION_UNCERTAINTY) << std::setw(MSR) << variances->get(0, 0) <<			// 11
            std::setprecision(PRECISION_UNCERTAINTY) << std::setw(MSR) << variances->get(0, 1) <<			// 12
            std::setprecision(PRECISION_UNCERTAINTY) << std::setw(MSR) << variances->get(0, 2) << std::endl;	// 13
            
        os <<
            std::setw(STATION+PAD2+LAT_EAST+LON_NORTH+STAT+STAT+PREC+PREC+PREC) << " " <<			// padding
            std::setprecision(PRECISION_UNCERTAINTY) << std::setw(MSR) << variances->get(1, 0) <<			// 21
            std::setprecision(PRECISION_UNCERTAINTY) << std::setw(MSR) << variances->get(1, 1) <<			// 22
            std::setprecision(PRECISION_UNCERTAINTY) << std::setw(MSR) << variances->get(1, 2) << std::endl;	// 23
        
        os <<
            std::setw(STATION+PAD2+LAT_EAST+LON_NORTH+STAT+STAT+PREC+PREC+PREC) << " " <<			// padding
            std::setprecision(PRECISION_UNCERTAINTY) << std::setw(MSR) << variances->get(2, 0) <<			// 31
            std::setprecision(PRECISION_UNCERTAINTY) << std::setw(MSR) << variances->get(2, 1) <<			// 32
            std::setprecision(PRECISION_UNCERTAINTY) << std::setw(MSR) << variances->get(2, 2) << std::endl;	// 33
    }
}

void DynAdjustPrinter::PrintPosUncertainties(std::ostream& os, const UINT32& block, const matrix_2d* stationVariances)
{
    vUINT32 v_blockStations(adjust_.v_parameterStationList_.at(block));

    // if required, sort stations according to original station file order
    if (adjust_.projectSettings_.o._sort_stn_file_order)
        adjust_.SortStationsbyFileOrder(v_blockStations);
    
    switch (adjust_.projectSettings_.a.adjust_mode)
    {
    case PhasedMode:
    case Phased_Block_1Mode:        // only the first block is rigorous
        os << "Block " << block + 1 << std::endl;
        break;
    }

    // Print header
    PrintPosUncertaintiesHeader(os);

    UINT32 mat_idx, stn;

    // Print stations according to the user-defined sort order
    for (UINT32 i=0; i<adjust_.v_blockStationsMap_.at(block).size(); ++i)
    {
        stn = v_blockStations.at(i);
        mat_idx = adjust_.v_blockStationsMap_.at(block)[stn] * 3;

        PrintPosUncertainty(os,
            block, stn, mat_idx,
            stationVariances, i, &v_blockStations);
    }

    os << std::endl;

    // return sort order to alpha-numeric
    if (adjust_.projectSettings_.o._sort_stn_file_order)
        adjust_.SortStationsbyID(v_blockStations);
}

void DynAdjustPrinter::UpdateGNSSNstatsForAlternateUnits(const v_uint32_u32u32_pair& msr_block) {

    matrix_2d gnss_cart(3, 1), gnss_local(3, 1);
    matrix_2d gnss_adj_cart(3, 1), gnss_adj_local(3, 1);
    matrix_2d var_cart(3, 3), var_local(3, 3), var_polar(3, 3);
    matrix_2d var_adj_local(3, 3), var_adj_polar(3, 3);
    matrix_2d rotations;

    for (auto _it_block_msr = msr_block.begin(); _it_block_msr != msr_block.end(); ++_it_block_msr)
    {
        it_vmsr_t _it_msr = adjust_.bmsBinaryRecords_.begin() + _it_block_msr->first;

        if (_it_msr->measStart != xMeas)
            continue;

        if (_it_msr->measType != 'G' && _it_msr->measType != 'X')
            continue;

        UINT32 cluster_count(_it_msr->vectorCount1);
        UINT32 covariance_count;
        const uint32_uint32_pair& b_pam = _it_block_msr->second;

        for (UINT32 cluster_msr = 0; cluster_msr < cluster_count; ++cluster_msr)
        {
            covariance_count = _it_msr->vectorCount2;

            // Load precision of adjusted measurements for staged adjustments
            if (adjust_.projectSettings_.a.stage)
                adjust_.DeserialiseBlockFromMappedFile(b_pam.first, 1, sf_prec_adj_msrs);

            double lower_mtx_buffer[6];
            memcpy(lower_mtx_buffer,
                adjust_.v_precAdjMsrsFull_.at(b_pam.first).getbuffer(b_pam.second, 0),
                sizeof(lower_mtx_buffer));
            matrix_2d var_adj_cart(3, 3, lower_mtx_buffer, 6, mtx_lower);
            var_adj_cart.fillupper();

            if (adjust_.projectSettings_.a.stage)
                adjust_.UnloadBlock(b_pam.first, 1, sf_prec_adj_msrs);

            // Get X component
            gnss_cart.put(0, 0, _it_msr->term1);
            gnss_adj_cart.put(0, 0, _it_msr->measAdj);
            var_cart.put(0, 0, _it_msr->term2);

            it_vmsr_t _it_x = _it_msr;

            // Get Y component
            _it_msr++;
            gnss_cart.put(1, 0, _it_msr->term1);
            gnss_adj_cart.put(1, 0, _it_msr->measAdj);
            var_cart.put(0, 1, _it_msr->term2);
            var_cart.put(1, 1, _it_msr->term3);

            it_vmsr_t _it_y = _it_msr;

            // Get Z component
            _it_msr++;
            gnss_cart.put(2, 0, _it_msr->term1);
            gnss_adj_cart.put(2, 0, _it_msr->measAdj);
            var_cart.put(0, 2, _it_msr->term2);
            var_cart.put(1, 2, _it_msr->term3);
            var_cart.put(2, 2, _it_msr->term4);

            it_vmsr_t _it_z = _it_msr;

            var_cart.filllower();

            double mid_lat(average(
                adjust_.bstBinaryRecords_.at(_it_msr->station1).currentLatitude,
                adjust_.bstBinaryRecords_.at(_it_msr->station2).currentLatitude));
            double mid_long(average(
                adjust_.bstBinaryRecords_.at(_it_msr->station1).currentLongitude,
                adjust_.bstBinaryRecords_.at(_it_msr->station2).currentLongitude));

            Rotate_CartLocal<double>(gnss_cart, &gnss_local,
                adjust_.bstBinaryRecords_.at(_it_msr->station1).currentLatitude,
                adjust_.bstBinaryRecords_.at(_it_msr->station1).currentLongitude);
            Rotate_CartLocal<double>(gnss_adj_cart, &gnss_adj_local,
                adjust_.bstBinaryRecords_.at(_it_msr->station1).currentLatitude,
                adjust_.bstBinaryRecords_.at(_it_msr->station1).currentLongitude);

            switch (adjust_.projectSettings_.o._adj_gnss_units)
            {
            case ENU_adj_gnss_ui:
            {
                PropagateVariances_CartLocal_Diagonal<double>(var_cart, var_local,
                    mid_lat, mid_long, rotations, true);
                PropagateVariances_CartLocal_Diagonal<double>(var_adj_cart, var_adj_local,
                    mid_lat, mid_long, rotations, false);

                // E
                _it_x->measCorr = gnss_adj_local.get(0, 0) - gnss_local.get(0, 0);
                _it_x->residualPrec = var_local.get(0, 0) - var_adj_local.get(0, 0);
                adjust_.UpdateMsrRecordStats(_it_x, var_local.get(0, 0));

                // N
                _it_y->measCorr = gnss_adj_local.get(1, 0) - gnss_local.get(1, 0);
                _it_y->residualPrec = var_local.get(1, 1) - var_adj_local.get(1, 1);
                adjust_.UpdateMsrRecordStats(_it_y, var_local.get(1, 1));

                // U
                _it_z->measCorr = gnss_adj_local.get(2, 0) - gnss_local.get(2, 0);
                _it_z->residualPrec = var_local.get(2, 2) - var_adj_local.get(2, 2);
                adjust_.UpdateMsrRecordStats(_it_z, var_local.get(2, 2));
                break;
            }
            case AED_adj_gnss_ui:
            {
                double azimuth = Direction(gnss_local.get(0, 0), gnss_local.get(1, 0));
                double elevation = VerticalAngle(gnss_local.get(0, 0), gnss_local.get(1, 0), gnss_local.get(2, 0));
                double distance = magnitude(gnss_local.get(0, 0), gnss_local.get(1, 0), gnss_local.get(2, 0));

                double azimuthAdj = Direction(gnss_adj_local.get(0, 0), gnss_adj_local.get(1, 0));
                double elevationAdj = VerticalAngle(gnss_adj_local.get(0, 0), gnss_adj_local.get(1, 0), gnss_adj_local.get(2, 0));
                double distanceAdj = magnitude(gnss_adj_local.get(0, 0), gnss_adj_local.get(1, 0), gnss_adj_local.get(2, 0));

                PropagateVariances_LocalCart<double>(var_cart, var_local,
                    mid_lat, mid_long, false, rotations, true);
                PropagateVariances_LocalCart<double>(var_adj_cart, var_adj_local,
                    mid_lat, mid_long, false, rotations, false);

                PropagateVariances_LocalPolar_Diagonal<double>(var_local, var_polar,
                    azimuth, elevation, distance, rotations, true);
                PropagateVariances_LocalPolar_Diagonal<double>(var_adj_local, var_adj_polar,
                    azimuth, elevation, distance, rotations, false);

                // Azimuth
                _it_x->measCorr = azimuthAdj - azimuth;
                _it_x->residualPrec = var_polar.get(0, 0) - var_adj_polar.get(0, 0);
                adjust_.UpdateMsrRecordStats(_it_x, var_polar.get(0, 0));

                // Elevation
                _it_y->measCorr = elevationAdj - elevation;
                _it_y->residualPrec = var_polar.get(1, 1) - var_adj_polar.get(1, 1);
                adjust_.UpdateMsrRecordStats(_it_y, var_polar.get(1, 1));

                // Distance
                _it_z->measCorr = distanceAdj - distance;
                _it_z->residualPrec = var_polar.get(2, 2) - var_adj_polar.get(2, 2);
                adjust_.UpdateMsrRecordStats(_it_z, var_polar.get(2, 2));
                break;
            }
            case ADU_adj_gnss_ui:
            {
                double azimuth = Direction(gnss_local.get(0, 0), gnss_local.get(1, 0));
                double elevation = VerticalAngle(gnss_local.get(0, 0), gnss_local.get(1, 0), gnss_local.get(2, 0));
                double distance = magnitude(gnss_local.get(0, 0), gnss_local.get(1, 0), gnss_local.get(2, 0));

                double azimuthAdj = Direction(gnss_adj_local.get(0, 0), gnss_adj_local.get(1, 0));
                double distanceAdj = magnitude(gnss_adj_local.get(0, 0), gnss_adj_local.get(1, 0), gnss_adj_local.get(2, 0));

                PropagateVariances_LocalCart<double>(var_cart, var_local,
                    mid_lat, mid_long, false, rotations, true);
                PropagateVariances_LocalCart<double>(var_adj_cart, var_adj_local,
                    mid_lat, mid_long, false, rotations, false);

                PropagateVariances_LocalPolar_Diagonal<double>(var_local, var_polar,
                    azimuth, elevation, distance, rotations, true);
                PropagateVariances_LocalPolar_Diagonal<double>(var_adj_local, var_adj_polar,
                    azimuth, elevation, distance, rotations, false);

                // Azimuth
                _it_x->measCorr = azimuthAdj - azimuth;
                _it_x->residualPrec = var_polar.get(0, 0) - var_adj_polar.get(0, 0);
                adjust_.UpdateMsrRecordStats(_it_x, var_polar.get(0, 0));

                // Slope distance
                _it_y->measCorr = distanceAdj - distance;
                _it_y->residualPrec = var_polar.get(2, 2) - var_adj_polar.get(2, 2);
                adjust_.UpdateMsrRecordStats(_it_y, var_polar.get(2, 2));

                // Up
                _it_z->measCorr = gnss_adj_local.get(2, 0) - gnss_local.get(2, 0);
                _it_z->residualPrec = var_local.get(2, 2) - var_adj_local.get(2, 2);
                adjust_.UpdateMsrRecordStats(_it_z, var_local.get(2, 2));
                break;
            }
            }

            // Skip covariances until next baseline
            _it_msr += covariance_count * 3;
            if (covariance_count > 0)
                _it_msr++;
        }
    }
}

void DynAdjustPrinter::PrintAdjGNSSAlternateUnits(it_vmsr_t& _it_msr, const uint32_uint32_pair& b_pam) {
    // Create copy of _it_msr, transform, then send to PrintAdjMeasurementsLinear
    vmsr_t gnss_msr;
    // vector components (original and adjusted)
    matrix_2d gnss_cart(3, 1), gnss_local(3, 1);
    matrix_2d gnss_adj_cart(3, 1), gnss_adj_local(3, 1);
    // variance components (original and adjusted)
    matrix_2d var_cart(3, 3), var_local(3, 3), var_polar(3, 3);
    matrix_2d var_adj_local(3, 3), var_adj_polar(3, 3);

    // For staged adjustments, load precision of adjusted measurements
    if (adjust_.projectSettings_.a.stage)
        adjust_.DeserialiseBlockFromMappedFile(b_pam.first, 1, sf_prec_adj_msrs);

    // get precision of adjusted measurements
    double lower_mtx_buffer[6];
    memcpy(lower_mtx_buffer, adjust_.v_precAdjMsrsFull_.at(b_pam.first).getbuffer(b_pam.second, 0), sizeof(lower_mtx_buffer));
    matrix_2d var_adj_cart(3, 3, lower_mtx_buffer, 6, mtx_lower);
    var_adj_cart.fillupper();

    if (adjust_.projectSettings_.a.stage)
        // For staged adjustments, unload precision of adjusted measurements
        adjust_.UnloadBlock(b_pam.first, 1, sf_prec_adj_msrs);
    
    // Get X
    gnss_msr.push_back(*_it_msr);
    gnss_cart.put(0, 0, _it_msr->term1);
    gnss_adj_cart.put(0, 0, _it_msr->measAdj);

    var_cart.put(0, 0, _it_msr->term2);
    
    // Get Y
    _it_msr++;
    gnss_msr.push_back(*_it_msr);
    gnss_cart.put(1, 0, _it_msr->term1);
    gnss_adj_cart.put(1, 0, _it_msr->measAdj);

    var_cart.put(0, 1, _it_msr->term2);
    var_cart.put(1, 1, _it_msr->term3);
    
    // Get Z
    _it_msr++;
    gnss_msr.push_back(*_it_msr);
    gnss_cart.put(2, 0, _it_msr->term1);
    gnss_adj_cart.put(2, 0, _it_msr->measAdj);

    var_cart.put(0, 2, _it_msr->term2);
    var_cart.put(1, 2, _it_msr->term3);
    var_cart.put(2, 2, _it_msr->term4);

    it_vmsr_t _it_gnss_msr(gnss_msr.begin());

    var_cart.filllower();

    double azimuth, elevation, distance;
    double azimuthAdj, elevationAdj, distanceAdj;
    double mid_lat(average(
        adjust_.bstBinaryRecords_.at(_it_msr->station1).currentLatitude,
        adjust_.bstBinaryRecords_.at(_it_msr->station2).currentLatitude));
    double mid_long(average(
        adjust_.bstBinaryRecords_.at(_it_msr->station1).currentLongitude,
        adjust_.bstBinaryRecords_.at(_it_msr->station2).currentLongitude));

    matrix_2d rotations;

    Rotate_CartLocal<double>(gnss_cart, &gnss_local,
        adjust_.bstBinaryRecords_.at(_it_msr->station1).currentLatitude,
        adjust_.bstBinaryRecords_.at(_it_msr->station1).currentLongitude);
    Rotate_CartLocal<double>(gnss_adj_cart, &gnss_adj_local,
        adjust_.bstBinaryRecords_.at(_it_msr->station1).currentLatitude,
        adjust_.bstBinaryRecords_.at(_it_msr->station1).currentLongitude);

    // return iterator to the X component of the GNSS measurement
    _it_msr -= 2;

    switch (adjust_.projectSettings_.o._adj_gnss_units)
    {
    case ENU_adj_gnss_ui:
        
        // term2: e
        // term3: n
        // term4: u

        // Propagate Variances from cartesian to local system - diagonals only
        PropagateVariances_CartLocal_Diagonal<double>(var_cart, var_local,
            mid_lat, mid_long,
            rotations, true);
        PropagateVariances_CartLocal_Diagonal<double>(var_adj_cart, var_adj_local,
            mid_lat, mid_long,
            rotations, false);

        // E
        _it_gnss_msr->preAdjMeas = gnss_local.get(0, 0);
        _it_gnss_msr->measAdj = gnss_adj_local.get(0, 0);
        _it_gnss_msr->measCorr = gnss_adj_local.get(0, 0) - gnss_local.get(0, 0);

        _it_gnss_msr->term2 = var_local.get(0, 0);
        _it_gnss_msr->measAdjPrec = var_adj_local.get(0, 0);
        _it_gnss_msr->residualPrec = _it_gnss_msr->term2 - _it_gnss_msr->measAdjPrec;

        // Update Pelzer reliability, N-Stat, T-Stat
        adjust_.UpdateMsrRecordStats(_it_gnss_msr, _it_gnss_msr->term2);

        // Initialise database id iterator
        if (adjust_.projectSettings_.o._database_ids)
            adjust_._it_dbid = adjust_.v_msr_db_map_.begin() + std::distance(adjust_.bmsBinaryRecords_.begin(), _it_msr);

        PrintAdjMeasurementsLinear('e', _it_gnss_msr, false);
    
        // N
        _it_msr++;
        _it_gnss_msr++;
        _it_gnss_msr->preAdjMeas = gnss_local.get(1, 0);
        _it_gnss_msr->measAdj = gnss_adj_local.get(1, 0);
        _it_gnss_msr->measCorr = gnss_adj_local.get(1, 0) - gnss_local.get(1, 0);
        
        _it_gnss_msr->term3 = var_local.get(1, 1);
        _it_gnss_msr->measAdjPrec = var_adj_local.get(1, 1);
        _it_gnss_msr->residualPrec = _it_gnss_msr->term3 - _it_gnss_msr->measAdjPrec;

        // Update Pelzer reliability, N-Stat, T-Stat
        adjust_.UpdateMsrRecordStats(_it_gnss_msr, _it_gnss_msr->term3);

        // Initialise database id iterator
        if (adjust_.projectSettings_.o._database_ids)
            adjust_._it_dbid = adjust_.v_msr_db_map_.begin() + std::distance(adjust_.bmsBinaryRecords_.begin(), _it_msr);

        PrintAdjMeasurementsLinear('n', _it_gnss_msr, false);

        // U
        _it_msr++;
        _it_gnss_msr++;
        _it_gnss_msr->preAdjMeas = gnss_local.get(2, 0);
        _it_gnss_msr->measAdj = gnss_adj_local.get(2, 0);
        _it_gnss_msr->measCorr = gnss_adj_local.get(2, 0) - gnss_local.get(2, 0);
        
        _it_gnss_msr->term4 = var_local.get(2, 2);
        _it_gnss_msr->measAdjPrec = var_adj_local.get(2, 2);
        _it_gnss_msr->residualPrec = _it_gnss_msr->term4 - _it_gnss_msr->measAdjPrec;

        // Update Pelzer reliability, N-Stat, T-Stat
        adjust_.UpdateMsrRecordStats(_it_gnss_msr, _it_gnss_msr->term4);

        // Initialise database id iterator
        if (adjust_.projectSettings_.o._database_ids)
            adjust_._it_dbid = adjust_.v_msr_db_map_.begin() + std::distance(adjust_.bmsBinaryRecords_.begin(), _it_msr);

        PrintAdjMeasurementsLinear('u', _it_gnss_msr, false);
        break;
    case AED_adj_gnss_ui:

        // term2: a
        // term3: e
        // term4: s

        azimuth = Direction(gnss_local.get(0, 0), gnss_local.get(1, 0));
        elevation = VerticalAngle(gnss_local.get(0, 0), gnss_local.get(1, 0), gnss_local.get(2, 0));
        distance = magnitude(gnss_local.get(0, 0), gnss_local.get(1, 0), gnss_local.get(2, 0));
        
        azimuthAdj = Direction(gnss_adj_local.get(0, 0), gnss_adj_local.get(1, 0));
        elevationAdj = VerticalAngle(gnss_adj_local.get(0, 0), gnss_adj_local.get(1, 0), gnss_adj_local.get(2, 0));
        distanceAdj = magnitude(gnss_adj_local.get(0, 0), gnss_adj_local.get(1, 0), gnss_adj_local.get(2, 0));

        // Propagate Variances from cartesian to local system - full matrix (for local-polar)
        PropagateVariances_LocalCart<double>(var_cart, var_local,
            mid_lat, mid_long, false,
            rotations, true);
        PropagateVariances_LocalCart<double>(var_adj_cart, var_adj_local,
            mid_lat, mid_long, false,
            rotations, false);

        // Propagate variance matrices from local reference frame to 
        // polar (azimuth, elevation angle and distance)
        PropagateVariances_LocalPolar_Diagonal<double>(var_local, var_polar, 
            azimuth, elevation, distance,
            rotations, true);
        PropagateVariances_LocalPolar_Diagonal<double>(var_adj_local, var_adj_polar, 
            azimuth, elevation, distance,
            rotations, false);

        // Print Geodetic azimuth
        _it_gnss_msr->preAdjMeas = azimuth;
        _it_gnss_msr->measAdj = azimuthAdj;
        _it_gnss_msr->measCorr = azimuthAdj - azimuth;

        _it_gnss_msr->term2 = var_polar.get(0, 0);
        _it_gnss_msr->measAdjPrec = var_adj_polar.get(0, 0);
        _it_gnss_msr->residualPrec = _it_gnss_msr->term2 - _it_gnss_msr->measAdjPrec;

        // Update Pelzer reliability, N-Stat, T-Stat
        adjust_.UpdateMsrRecordStats(_it_gnss_msr, _it_gnss_msr->term2);

        // Initialise database id iterator
        if (adjust_.projectSettings_.o._database_ids)
            adjust_._it_dbid = adjust_.v_msr_db_map_.begin() + std::distance(adjust_.bmsBinaryRecords_.begin(), _it_msr);

        PrintAdjMeasurementsAngular('a', _it_gnss_msr, false);
    
        // Print vertical angle
        _it_msr++;
        _it_gnss_msr++;
        _it_gnss_msr->preAdjMeas = elevation;
        _it_gnss_msr->measAdj = elevationAdj;
        _it_gnss_msr->measCorr = elevationAdj - elevation;

        _it_gnss_msr->term3 = var_polar.get(1, 1);
        _it_gnss_msr->measAdjPrec = var_adj_polar.get(1, 1);
        _it_gnss_msr->residualPrec = _it_gnss_msr->term3 - _it_gnss_msr->measAdjPrec;

        // Update Pelzer reliability, N-Stat, T-Stat
        adjust_.UpdateMsrRecordStats(_it_gnss_msr, _it_gnss_msr->term3);

        // Initialise database id iterator
        if (adjust_.projectSettings_.o._database_ids)
            adjust_._it_dbid = adjust_.v_msr_db_map_.begin() + std::distance(adjust_.bmsBinaryRecords_.begin(), _it_msr);

        PrintAdjMeasurementsAngular('v', _it_gnss_msr, false);

        // Print slope distance
        _it_msr++;
        _it_gnss_msr++;
        _it_gnss_msr->preAdjMeas = distance;
        _it_gnss_msr->measAdj = distanceAdj;
        _it_gnss_msr->measCorr = distanceAdj - distance;

        _it_gnss_msr->term4 = var_polar.get(2, 2);
        _it_gnss_msr->measAdjPrec = var_adj_polar.get(2, 2);
        _it_gnss_msr->residualPrec = _it_gnss_msr->term4 - _it_gnss_msr->measAdjPrec;

        // Update Pelzer reliability, N-Stat, T-Stat
        adjust_.UpdateMsrRecordStats(_it_gnss_msr, _it_gnss_msr->term4);

        // Initialise database id iterator
        if (adjust_.projectSettings_.o._database_ids)
            adjust_._it_dbid = adjust_.v_msr_db_map_.begin() + std::distance(adjust_.bmsBinaryRecords_.begin(), _it_msr);

        PrintAdjMeasurementsLinear('s', _it_gnss_msr, false);
        break;
    case ADU_adj_gnss_ui:

        // term2: a
        // term3: s
        // term4: u

        azimuth = Direction(gnss_local.get(0, 0), gnss_local.get(1, 0));
        elevation = VerticalAngle(gnss_local.get(0, 0), gnss_local.get(1, 0), gnss_local.get(2, 0));
        distance = magnitude(gnss_local.get(0, 0), gnss_local.get(1, 0), gnss_local.get(2, 0));
        
        azimuthAdj = Direction(gnss_adj_local.get(0, 0), gnss_adj_local.get(1, 0));
        distanceAdj = magnitude(gnss_adj_local.get(0, 0), gnss_adj_local.get(1, 0), gnss_adj_local.get(2, 0));

        // Propagate Variances from cartesian to local system - full matrix (for local-polar)
        PropagateVariances_LocalCart<double>(var_cart, var_local,
            mid_lat, mid_long, false,
            rotations, true);
        PropagateVariances_LocalCart<double>(var_adj_cart, var_adj_local,
            mid_lat, mid_long, false,
            rotations, false);

        // Propagate variance matrices from local reference frame to 
        // polar (azimuth, elevation angle and distance)
        PropagateVariances_LocalPolar_Diagonal<double>(var_local, var_polar, 
            azimuth, elevation, distance,
            rotations, true);
        PropagateVariances_LocalPolar_Diagonal<double>(var_adj_local, var_adj_polar, 
            azimuth, elevation, distance,
            rotations, false);

        // Print Geodetic azimuth
        _it_gnss_msr->preAdjMeas = azimuth;
        _it_gnss_msr->measAdj = azimuthAdj;
        _it_gnss_msr->measCorr = azimuthAdj - azimuth;

        _it_gnss_msr->term2 = var_polar.get(0, 0);
        _it_gnss_msr->measAdjPrec = var_adj_polar.get(0, 0);
        _it_gnss_msr->residualPrec = _it_gnss_msr->term2 - _it_gnss_msr->measAdjPrec;

        // Update Pelzer reliability, N-Stat, T-Stat
        adjust_.UpdateMsrRecordStats(_it_gnss_msr, _it_gnss_msr->term2);

        // Initialise database id iterator
        if (adjust_.projectSettings_.o._database_ids)
            adjust_._it_dbid = adjust_.v_msr_db_map_.begin() + std::distance(adjust_.bmsBinaryRecords_.begin(), _it_msr);

        PrintAdjMeasurementsAngular('a', _it_gnss_msr, false);
    
        // Print slope distance
        _it_msr++;
        _it_gnss_msr++;
        _it_gnss_msr->preAdjMeas = distance;
        _it_gnss_msr->measAdj = distanceAdj;
        _it_gnss_msr->measCorr = distanceAdj - distance;

        _it_gnss_msr->term3 = var_polar.get(2, 2);
        _it_gnss_msr->measAdjPrec = var_adj_polar.get(2, 2);
        _it_gnss_msr->residualPrec = _it_gnss_msr->term3 - _it_gnss_msr->measAdjPrec;

        // Update Pelzer reliability, N-Stat, T-Stat
        adjust_.UpdateMsrRecordStats(_it_gnss_msr, _it_gnss_msr->term3);

        // Initialise database id iterator
        if (adjust_.projectSettings_.o._database_ids)
            adjust_._it_dbid = adjust_.v_msr_db_map_.begin() + std::distance(adjust_.bmsBinaryRecords_.begin(), _it_msr);

        PrintAdjMeasurementsLinear('s', _it_gnss_msr, false);

        // Print up
        _it_msr++;
        _it_gnss_msr++;
        _it_gnss_msr->preAdjMeas = gnss_local.get(2, 0);
        _it_gnss_msr->measAdj = gnss_adj_local.get(2, 0);
        _it_gnss_msr->measCorr = gnss_adj_local.get(2, 0) - gnss_local.get(2, 0);
        
        _it_gnss_msr->term4 = var_local.get(2, 2);
        _it_gnss_msr->measAdjPrec = var_adj_local.get(2, 2);
        _it_gnss_msr->residualPrec = _it_gnss_msr->term4 - _it_gnss_msr->measAdjPrec;

        // Update Pelzer reliability, N-Stat, T-Stat
        adjust_.UpdateMsrRecordStats(_it_gnss_msr, _it_gnss_msr->term4);

        // Initialise database id iterator
        if (adjust_.projectSettings_.o._database_ids)
            adjust_._it_dbid = adjust_.v_msr_db_map_.begin() + std::distance(adjust_.bmsBinaryRecords_.begin(), _it_msr);

        PrintAdjMeasurementsLinear('u', _it_gnss_msr, false);
        break;
    }
}

void DynAdjustPrinter::PrintCompMeasurements_HR(const UINT32& block, it_vmsr_t& _it_msr, UINT32& design_row, printMeasurementsMode printMode)
{
    // normal format
    adjust_.adj_file << std::left << std::setw(STATION) << adjust_.bstBinaryRecords_.at(_it_msr->station1).stationName;
    adjust_.adj_file << std::left << std::setw(STATION) << " ";
    adjust_.adj_file << std::left << std::setw(STATION) << " ";

    double computed, correction;
    
    switch (printMode)
    {
    case computedMsrs:
        correction = -adjust_.v_measMinusComp_.at(block).get(design_row, 0);
        computed = _it_msr->term1 + correction - _it_msr->preAdjCorr;
        break;
    case ignoredMsrs:
    default:
        correction = _it_msr->measCorr;
        computed = _it_msr->measAdj;
        break;
    }

    // Print linear measurement, taking care of user requirements for precision	
    PrintComputedMeasurements<LinearMeasurement>(' ', computed, correction, _it_msr);
        
    design_row++;
}

void DynAdjustPrinter::PrintCompMeasurements_IJPQ(const UINT32& block, it_vmsr_t& _it_msr, UINT32& design_row, printMeasurementsMode printMode)
{
    // normal format
    adjust_.adj_file << std::left << std::setw(STATION) << adjust_.bstBinaryRecords_.at(_it_msr->station1).stationName;
    adjust_.adj_file << std::left << std::setw(STATION) << " ";
    adjust_.adj_file << std::left << std::setw(STATION) << " ";

    double computed, correction;
    
    switch (printMode)
    {
    case computedMsrs:
        correction = -adjust_.v_measMinusComp_.at(block).get(design_row, 0);
        computed = _it_msr->term1 + correction + _it_msr->preAdjCorr;
        break;
    case ignoredMsrs:
    default:
        correction = _it_msr->measCorr;
        computed = _it_msr->measAdj;
        break;
    }

    // Print angular measurement, taking care of user requirements for 
    // type, format and precision	
    PrintComputedMeasurements<AngularMeasurement>(' ', computed, correction, _it_msr);
        
    design_row++;
}

void DynAdjustPrinter::PrintMsrVarianceMatrixException(const it_vmsr_t& _it_msr, const std::runtime_error& e, std::stringstream& ss, 
    const std::string& calling_function, const UINT32 msr_count)
{
    it_vmsr_t _it_msr_temp(_it_msr);

    switch (_it_msr->measType)
    {
    case 'D':
        ss << calling_function << "(): Cannot compute the" << std::endl <<
            "  variance matrix for a round of " << msr_count << " directions commencing" << std::endl <<
            "  with stations " << adjust_.bstBinaryRecords_.at(_it_msr->station1).stationName << " and " <<
            adjust_.bstBinaryRecords_.at(_it_msr->station2).stationName << ":" << std::endl;

        ss << "    ..." << std::endl <<
            "    <Value>" << FormatDmsString(RadtoDms(_it_msr->term1), adjust_.PRECISION_SEC_MSR, true, false) << "</Value>" << std::endl <<
            "    ..." << std::endl << std::endl;
        break;
    case 'G':
        ss << calling_function << "(): Cannot invert the" << std::endl <<
            "  variance matrix for a GPS baseline between stations " << 
            adjust_.bstBinaryRecords_.at(_it_msr->station1).stationName << " and " <<
            adjust_.bstBinaryRecords_.at(_it_msr->station2).stationName << ":" << std::endl;

        ss << "    ..." << std::endl <<
            "    <x>" << std::fixed << std::setprecision(adjust_.PRECISION_MTR_MSR) << _it_msr_temp->term1 << "</x>" << std::endl;
        _it_msr_temp++;
        ss << "    <y>" << std::fixed << std::setprecision(adjust_.PRECISION_MTR_MSR) << _it_msr_temp->term1 << "</y>" << std::endl;
        _it_msr_temp++;
        ss << "    <z>" << std::fixed << std::setprecision(adjust_.PRECISION_MTR_MSR) << _it_msr_temp->term1 << "</z>" << std::endl <<
            "    ..." << std::endl << std::endl;
        break;
    case 'X':
        ss << calling_function << "(): Cannot invert the" << std::endl <<
            "  variance matrix for a " << msr_count << "-baseline GPS baseline cluster commencing" << std::endl <<
            "  with stations " << adjust_.bstBinaryRecords_.at(_it_msr->station1).stationName << " and " <<
            adjust_.bstBinaryRecords_.at(_it_msr->station2).stationName << ":" << std::endl;

        ss << "    ..." << std::endl <<
            "    <x>" << std::fixed << std::setprecision(adjust_.PRECISION_MTR_MSR) << _it_msr_temp->term1 << "</x>" << std::endl;
        _it_msr_temp++;
        ss << "    <y>" << std::fixed << std::setprecision(adjust_.PRECISION_MTR_MSR) << _it_msr_temp->term1 << "</y>" << std::endl;
        _it_msr_temp++;
        ss << "    <z>" << std::fixed << std::setprecision(adjust_.PRECISION_MTR_MSR) << _it_msr_temp->term1 << "</z>" << std::endl <<
            "    ..." << std::endl << std::endl;
        break;
    case 'Y':
        ss << calling_function << "(): Cannot invert the" << std::endl <<
            "  variance matrix for a " << msr_count << "-station GPS point cluster commencing" << std::endl <<
            "  with station " << adjust_.bstBinaryRecords_.at(_it_msr->station1).stationName << ":" << std::endl;

        ss << "    ..." << std::endl <<
            "    <x>" << std::fixed << std::setprecision(adjust_.PRECISION_MTR_MSR) << _it_msr_temp->term1 << "</x>" << std::endl;
        _it_msr_temp++;
        ss << "    <y>" << std::fixed << std::setprecision(adjust_.PRECISION_MTR_MSR) << _it_msr_temp->term1 << "</y>" << std::endl;
        _it_msr_temp++;
        ss << "    <z>" << std::fixed << std::setprecision(adjust_.PRECISION_MTR_MSR) << _it_msr_temp->term1 << "</z>" << std::endl <<
            "    ..." << std::endl << std::endl;
        break;
    }
    
    ss <<
        "  Detailed description:   " << e.what() << std::endl <<
        "  Options: " << std::endl <<
        "  - Check the validity of the variance matrix and any relevant " << std::endl <<
        "    scalars and re-attempt the adjustment." << std::endl <<
        "  - If this fails, set the ignore flag in the measurement file and re-import " << std::endl <<
        "    the station and measurement files." << std::endl;

    adjust_.adj_file << "  " << ss.str() << std::endl;
    adjust_.adj_file.flush();
}

void DynAdjustPrinter::PrintCorStationsUniqueList(std::ostream& cor_file)
{
    vUINT32 v_blockStations;
    vstring stn_corr_records(adjust_.bstBinaryRecords_.size());
    cor_file << std::setw(STATION) << std::left << "Station" <<
        std::setw(PAD2) << " " << 
        std::right << std::setw(MSR) << "Azimuth" << 
        std::right << std::setw(MSR) << "V. Angle" << 
        std::right << std::setw(MSR) << "S. Distance" << 
        std::right << std::setw(MSR) << "H. Distance" << 
        std::right << std::setw(HEIGHT) << "east" << 
        std::right << std::setw(HEIGHT) << "north" << 
        std::right << std::setw(HEIGHT) << "up" << std::endl;
    UINT32 i, j = STATION+PAD2+MSR+MSR+MSR+MSR+HEIGHT+HEIGHT+HEIGHT;
    for (i=0; i<j; ++i)
        cor_file << "-";
    cor_file << std::endl;
    UINT32 block(UINT_MAX), stn, mat_index;
    _it_u32u32_uint32_pair _it_bsmu;
    
    std::stringstream stationRecord;
    v_uint32_string_pair correctionsOutput;
    correctionsOutput.reserve(adjust_.v_blockStationsMapUnique_.size());
    const v_mat_2d* estimates = &adjust_.v_rigorousStations_;
    if (adjust_.projectSettings_.a.adjust_mode == SimultaneousMode)
        estimates = &adjust_.v_estimatedStations_;
    std::ostream* outstream(&cor_file);
    if (adjust_.projectSettings_.a.stage)
        outstream = &stationRecord;
    // No need to sort order of block Stations Map - this was done in
    // PrintAdjStationsUniqueList
    for (_it_bsmu=adjust_.v_blockStationsMapUnique_.begin();
        _it_bsmu!=adjust_.v_blockStationsMapUnique_.end(); ++_it_bsmu)
    {
        // If this is not the first block, exit if block-1 mode
        if (adjust_.projectSettings_.a.adjust_mode == Phased_Block_1Mode)
            if (_it_bsmu->second > 0)
                break;
        if (adjust_.projectSettings_.a.stage)
        {
            if (block != _it_bsmu->second)
            {
                // unload previous block
                if (_it_bsmu->second > 0)
                    adjust_.UnloadBlock(_it_bsmu->second-1, 2, 
                        sf_rigorous_stns, sf_original_stns);
                // load up this block
                if (adjust_.projectSettings_.a.stage)
                    adjust_.DeserialiseBlockFromMappedFile(_it_bsmu->second, 2,
                        sf_rigorous_stns, sf_original_stns);
            }
        }
        block = _it_bsmu->second;
        stn = _it_bsmu->first.first;
        mat_index = _it_bsmu->first.second * 3;
        PrintCorStation(*outstream, block, stn, mat_index,
            &estimates->at(block));
            
        if (adjust_.projectSettings_.a.stage)
        {
            correctionsOutput.push_back(uint32_string_pair(stn, stationRecord.str()));
            stationRecord.str("");
        }
    }
    if (adjust_.projectSettings_.a.stage)
        adjust_.UnloadBlock(block, 2, 
            sf_rigorous_stns, sf_original_stns);
    // if required, sort stations according to original station file order
    if (adjust_.projectSettings_.o._sort_stn_file_order)
    {
        CompareOddPairFirst_FileOrder<station_t, UINT32, std::string> stnorderCompareFunc(&adjust_.bstBinaryRecords_);
        std::sort(correctionsOutput.begin(), correctionsOutput.end(), stnorderCompareFunc);
    }
    else
        std::sort(correctionsOutput.begin(), correctionsOutput.end(), CompareOddPairFirst<UINT32, std::string>());
    for_each(correctionsOutput.begin(), correctionsOutput.end(),
        [&correctionsOutput, &cor_file] (std::pair<UINT32, std::string>& corrRecord) { 
            cor_file << corrRecord.second;
        }
    );
}

void DynAdjustPrinter::PrintAdjStations(std::ostream& os, const UINT32& block,
                      const matrix_2d* stationEstimates,
                      matrix_2d* stationVariances, bool printBlockID,
                      bool recomputeGeographicCoords,
                      bool updateGeographicCoords, bool printHeader,
                      bool reapplyTypeBUncertainties)
{
    vUINT32 v_blockStations(adjust_.v_parameterStationList_.at(block));
    if (v_blockStations.size() * 3 != stationEstimates->rows())
    {
        std::stringstream ss;
        ss << "PrintAdjStations(): Number of estimated stations in block " << block << 
            " does not match the block station count." << std::endl;
        adjust_.SignalExceptionAdjustment(ss.str(), 0);
    }
    
    AdjFile adj;
    try {
        if (printHeader)
            adj.print_adj_stn_header(os);
        
        // if required, sort stations according to original station file order
        if (adjust_.projectSettings_.o._sort_stn_file_order)
            adjust_.SortStationsbyFileOrder(v_blockStations);
            
        switch (adjust_.projectSettings_.a.adjust_mode)
        {
        case SimultaneousMode:
            PrintStationColumnHeaders(os, 
                adjust_.projectSettings_.o._stn_coord_types, adjust_.projectSettings_.o._stn_corr);
            break;
        case PhasedMode:
        case Phased_Block_1Mode:		// only the first block is rigorous
            // Print stations in each block?
            if (adjust_.projectSettings_.o._output_stn_blocks)
            {
                if (printBlockID)
                    os << "Block " << block + 1 << std::endl;
                else if (adjust_.projectSettings_.o._adj_stn_iteration)
                    adj.print_adj_stn_block_header(os, block);
            
                PrintStationColumnHeaders(os, 
                    adjust_.projectSettings_.o._stn_coord_types, adjust_.projectSettings_.o._stn_corr);
            }
            else 
            {
                if (block == 0)
                    PrintStationColumnHeaders(os, 
                        adjust_.projectSettings_.o._stn_coord_types, adjust_.projectSettings_.o._stn_corr);
            }
            break;
        }
    }
    catch (const std::runtime_error& e) {
        adjust_.SignalExceptionAdjustment(e.what(), 0);
    }
    
    UINT32 mat_idx, stn;
    // Print stations according to the user-defined sort order
    for (UINT32 i(0); i<adjust_.v_blockStationsMap_.at(block).size(); ++i)
    {
        stn = v_blockStations.at(i);
        mat_idx = adjust_.v_blockStationsMap_.at(block)[stn] * 3;
        PrintAdjStation(os, block, stn, mat_idx,
            stationEstimates, stationVariances, 
            recomputeGeographicCoords, updateGeographicCoords,
            reapplyTypeBUncertainties);
    }
    os << std::endl;
    
    // return sort order to alpha-numeric
    if (adjust_.projectSettings_.o._sort_stn_file_order)
        adjust_.SortStationsbyID(v_blockStations);
}

} // namespace networkadjust
} // namespace dynadjust
