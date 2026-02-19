//============================================================================
// Name         : dnaadjust_printer.hpp
// Author       : Dale Roberts <dale.o.roberts@gmail.com>
// Contributors : 
// Copyright    : Copyright 2017-2025 Geoscience Australia
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
// Description  : DynAdjust Network Adjustment Printer Module
//============================================================================

#ifndef DNAADJUST_PRINTER_H_
#define DNAADJUST_PRINTER_H_

#if defined(_MSC_VER)
#if defined(LIST_INCLUDES_ON_BUILD)
#pragma message("  " __FILE__)
#endif
#endif

#include <array>
#include <iostream>
#include <string_view>
#include <type_traits>

#include <include/config/dnaexports.hpp>

#include <include/config/dnaconsts-iostream.hpp>
#include <include/config/dnaconsts.hpp>
#include <include/config/dnaoptions-interface.hpp>
#include <include/config/dnatypes.hpp>
#include <include/measurement_types/dnameasurement.hpp>
#include <include/measurement_types/dnastation.hpp>
#include <include/measurement_types/dnagpspoint.hpp>
#include <include/io/dnaiodna.hpp>
#include <include/functions/dnaiostreamfuncs.hpp>
#include <include/functions/dnatimer.hpp>

using namespace dynadjust::datum_parameters;
using namespace dynadjust::measurements;

namespace dynadjust {
namespace networkadjust {

// Forward declarations
class dna_adjust;

// Print context structure containing all necessary references for printing
struct PrintContext {
    std::ostream& adj_file;
    std::ostream& debug_file;
    const project_settings& settings;
    const vstn_t& binary_stations;
    const vmsr_t& binary_measurements;

    PrintContext(std::ostream& adj, std::ostream& debug, const project_settings& proj_settings,
                 const vstn_t& bst_records, const vmsr_t& bms_records)
        : adj_file(adj), debug_file(debug), settings(proj_settings), binary_stations(bst_records),
          binary_measurements(bms_records) {}
};

// Tag types for measurement classification
struct AngularMeasurement {};
struct LinearMeasurement {};
struct PrintAdjMeasurementsHeader {};
struct DirectionSetMeasurement {};

// Stage 4: Coordinate type tags for station formatting
struct GeographicCoordinates {};
struct CartesianCoordinates {};
struct ProjectionCoordinates {};

// Stage 4: Output mode enums
enum class UncertaintyMode {
    Ellipses,
    Covariances,
    Both
};

// Type traits for measurement classification
template <typename T> struct IsAngularMeasurement : std::false_type {};

template <> struct IsAngularMeasurement<AngularMeasurement> : std::true_type {};

template <typename T> constexpr bool kIsAngularMeasurementV = IsAngularMeasurement<T>::value;

// Printer class for DynAdjust output formatting
class DNAADJUST_API DynAdjustPrinter {
  public:
    // Constructor taking dna_adjust reference for direct access to members
    explicit DynAdjustPrinter(dna_adjust& adjust_instance);

    // Template-based measurement printing
    template <typename MeasurementType>
    void PrintAdjMeasurements(char cardinal, const it_vmsr_t& it_msr, bool initialise_dbindex = true);

    // Template-based computed measurement printing
    template <typename MeasurementType>
    void PrintComputedMeasurements(char cardinal, const double& computed, const double& correction, const it_vmsr_t& it_msr);
    
    // Compatibility methods - delegate to template implementations
    void PrintCompMeasurementsLinear(const char cardinal, const double& computed, const double& correction, const it_vmsr_t& it_msr);
    void PrintCompMeasurementsAngular(const char cardinal, const double& computed, const double& correction, const it_vmsr_t& it_msr);
    void PrintMeasurementsAngular(const char cardinal, const double& measurement, const double& correction, const it_vmsr_t& it_msr, bool printAdjMsr = true);
    void PrintMeasurementsLinear(const char cardinal, const double& measurement, const double& correction, const it_vmsr_t& it_msr, bool printAdjMsr = true);
    void PrintAdjMeasurementsAngular(const char cardinal, const it_vmsr_t& it_msr, bool initialise_dbindex = true);
    void PrintAdjMeasurementsLinear(const char cardinal, const it_vmsr_t& it_msr, bool initialise_dbindex = true);

    // Stage 5: Enhanced measurement formatting templates
    template <typename MeasurementType>
    void PrintMeasurementValue(char cardinal, const double& measurement, const double& correction, const it_vmsr_t& it_msr, bool printAdjMsr = true);

    // Consolidated station formatter for measurement types A, BKVZ, CELMS, HR, IJPQ
    void PrintMeasurementWithStations(it_vmsr_t& it_msr, char measurement_type);

    // Utility functions
    void PrintIteration(const UINT32& iteration);
    void PrintAdjustmentTime(cpu_timer& time, int timer_type);
    void PrintAdjustmentStatus();
    void PrintMeasurementDatabaseID(const it_vmsr_t& it_msr, bool initialise_dbindex = false);
    void PrintAdjMeasurementStatistics(char cardinal, const it_vmsr_t& it_msr, bool initialise_dbindex);
    
    // Stage 3: Header generation infrastructure
    void PrintMeasurementHeader(printMeasurementsMode printMode, std::string_view col1_heading, std::string_view col2_heading);
    void PrintAdjMeasurementsHeader(bool printHeader, const std::string& table_heading, printMeasurementsMode printMode, UINT32 block, bool printBlocks = false);
    void PrintPositionUncertaintyHeader(std::ostream& os);
    
    // Stage 3: Output coordinators
    void PrintAdjustedNetworkMeasurements();
    void PrintAdjustedNetworkStations();
    void PrintNetworkStationCorrections();
    
    // Stage 5: Complex measurement handlers
    void PrintAdjMeasurements(v_uint32_u32u32_pair msr_block, bool printHeader);
    void PrintIgnoredAdjMeasurements(bool printHeader);
    void PrintCompMeasurements(v_uint32_u32u32_pair msr_block, bool printHeader);
    void PrintMeasurementCorrection(const char cardinal, const it_vmsr_t& _it_msr);
    void PrintAdjMeasurements_YLLH(it_vmsr_t& _it_msr);
    void PrintPositionalUncertainty();
    void PrintEstimatedStationCoordinatestoDNAXML(const std::string& stnFile, INPUT_FILE_TYPE t, bool flagUnused = false);
    bool PrintEstimatedStationCoordinatestoSNX(std::string& sinex_filename);
    void PrintCompMeasurements(const UINT32& block, const std::string& type);
    void PrintPosUncertaintiesUniqueList(std::ostream& os, const v_mat_2d* stationVariances);
    void PrintStationCorrectionsList(std::ostream& cor_file);
    void PrintBlockStations(std::ostream& os, const UINT32& block, const matrix_2d* stationEstimates, 
                           matrix_2d* stationVariances, bool printBlockID, bool recomputeGeographicCoords, 
                           bool updateGeographicCoords, bool printHeader, bool reapplyTypeBUncertainties);
    void PrintOutputFileHeaderInfo();
    void PrintCompMeasurements_GXY(const UINT32& block, it_vmsr_t& _it_msr, UINT32& design_row, printMeasurementsMode printMode);
    void PrintAdjGNSSAlternateUnits(it_vmsr_t& _it_msr, const uint32_uint32_pair& b_pam);
    void UpdateGNSSNstatsForAlternateUnits(const v_uint32_u32u32_pair& msr_block);
    void PrintStationsUniqueList(std::ostream& os, const v_mat_2d* stationEstimates, v_mat_2d* stationVariances, 
                                 bool recomputeGeographicCoords, bool updateGeographicCoords, bool reapplyTypeBUncertainties);
    void PrintCompMeasurements_D(it_vmsr_t& _it_msr, UINT32& design_row, bool printIgnored);
    void PrintAdjMeasurements_D(it_vmsr_t& _it_msr);
    void PrintAdjMeasurements_A(it_vmsr_t& _it_msr);
    void PrintAdjMeasurements_BKVZ(it_vmsr_t& _it_msr);
    void PrintAdjMeasurements_CELMS(it_vmsr_t& _it_msr);
    void PrintAdjMeasurements_HR(it_vmsr_t& _it_msr);
    void PrintAdjMeasurements_IJPQ(it_vmsr_t& _it_msr);
    void PrintCompMeasurements_YLLH(it_vmsr_t& _it_msr, UINT32& design_row);
    void PrintCompMeasurements_A(const UINT32& block, it_vmsr_t& _it_msr, UINT32& design_row, printMeasurementsMode printMode);
    void PrintCompMeasurements_BKVZ(const UINT32& block, it_vmsr_t& _it_msr, UINT32& design_row, printMeasurementsMode printMode);
    void PrintCompMeasurements_CELMS(it_vmsr_t& _it_msr, UINT32& design_row, printMeasurementsMode printMode);
    void PrintCompMeasurements_HR(const UINT32& block, it_vmsr_t& _it_msr, UINT32& design_row, printMeasurementsMode printMode);
    void PrintCompMeasurements_IJPQ(const UINT32& block, it_vmsr_t& _it_msr, UINT32& design_row, printMeasurementsMode printMode);
    void PrintMsrVarianceMatrixException(const it_vmsr_t& _it_msr, const std::runtime_error& e, std::stringstream& ss, 
                                        const std::string& calling_function, const UINT32 msr_count);
    
    // Stage 3: Specialized measurement handlers
    void PrintGPSClusterMeasurements(it_vmsr_t& it_msr, const UINT32& block = 0);
    void PrintDirectionSetMeasurements(it_vmsr_t& it_msr, bool adjustedMsrs = true);
    
    // Stage 3: Statistical and summary generators
    void PrintStatistics(bool printPelzer = true);
    void PrintMeasurementsToStation();
    void PrintCorrelationStations(std::ostream& cor_file, const UINT32& block);

    // Stage 4: Station coordinate formatters
    template<typename CoordinateType>
    void PrintStationCoordinates(std::ostream& os, const it_vstn_t& stn_it, 
                                const matrix_2d* estimates = nullptr, 
                                const matrix_2d* variances = nullptr);
                                
    template<typename CoordinateType>
    void PrintStationUncertainties(std::ostream& os, const it_vstn_t& stn_it,
                                  const matrix_2d* variances, UncertaintyMode mode);

    // Stage 4: Station file headers
    void PrintStationFileHeader(std::ostream& os, std::string_view file_type, 
                               std::string_view filename);
    void PrintStationColumnHeaders(std::ostream& os, const std::string& stn_coord_types,
                                  const UINT16& printStationCorrections);
    void PrintPositionalUncertaintyHeader(std::ostream& os, std::string_view filename);

    // Stage 4: Station processing coordinators  
    // void PrintNetworkStationCorrections(); - Already declared above
    void PrintStationCorrelations(std::ostream& cor_file, const UINT32& block);
    void PrintStationsInBlock(std::ostream& os, const UINT32& block,
                             const matrix_2d* estimates, const matrix_2d* variances,
                             const std::string& stn_coord_types, const UINT16& printStationCorrections);
    void PrintUniqueStationsList(std::ostream& os, 
                                const matrix_2d* estimates, const matrix_2d* variances,
                                const std::string& stn_coord_types, const UINT16& printStationCorrections);
                                
    // Enhanced unique stations list processing
    void PrintAdjStationsUniqueList(std::ostream& os,
                                              const v_mat_2d* stationEstimates, v_mat_2d* stationVariances,
                                              bool recomputeGeographicCoords, bool updateGeographicCoords,
                                              bool reapplyTypeBUncertainties);

    // Stage 4: Advanced station functions
    void PrintPositionalUncertaintyOutput();
    void PrintStationAdjustmentResults(std::ostream& os, const UINT32& block,
                                      const UINT32& stn, const UINT32& mat_idx,
                                      const matrix_2d* estimates, matrix_2d* variances);

    // Stage 6: Export functions
    void PrintEstimatedStationCoordinatestoDNAXML_Y(const std::string& msrFile, INPUT_FILE_TYPE t);
    
    // Stage 7: Validation functions
    bool IgnoredMeasurementContainsInvalidStation(pit_vmsr_t _it_msr);
    
    // Enhanced station formatting
    void PrintAdjStation(std::ostream& os, const UINT32& block, const UINT32& stn, const UINT32& mat_idx,
                        const matrix_2d* stationEstimates, matrix_2d* stationVariances,
                        bool recomputeGeographicCoords, bool updateGeographicCoords, bool reapplyTypeBUncertainties);
    
    // GPS cluster measurement printing
    void PrintAdjMeasurements_GXY(it_vmsr_t& _it_msr, const uint32_uint32_pair& b_pam);
    
    // Station correction printing
    void PrintCorStation(std::ostream& os, const UINT32& block, const UINT32& stn, const UINT32& mat_index,
                        const matrix_2d* stationEstimates);
    void PrintCorStations(std::ostream& cor_file, const UINT32& block);
    void PrintCorStationsUniqueList(std::ostream& cor_file);
    void PrintAdjStations(std::ostream& os, const UINT32& block,
                          const matrix_2d* stationEstimates,
                          matrix_2d* stationVariances, bool printBlockID,
                          bool recomputeGeographicCoords,
                          bool updateGeographicCoords, bool printHeader,
                          bool reapplyTypeBUncertainties);
    void PrintPosUncertaintiesHeader(std::ostream& os);
    void PrintPosUncertainty(std::ostream& os, const UINT32& block, const UINT32& stn, 
                            const UINT32& mat_idx, const matrix_2d* stationVariances, 
                            const UINT32& map_idx, const vUINT32* blockStations);
    void PrintPosUncertainties(std::ostream& os, const UINT32& block, const matrix_2d* stationVariances);

    // Stage 4: Enhanced coordinate formatting utilities for PrintAdjStation refactoring
    void PrintStationCoordinatesByType(std::ostream& os, const it_vstn_t& stn_it,
                                      const matrix_2d* estimates, const UINT32& mat_idx,
                                      double est_lat, double est_lon, double est_height,
                                      double easting, double northing, double zone);
    void PrintLatitudeCoordinate(std::ostream& os, double latitude);
    void PrintLongitudeCoordinate(std::ostream& os, double longitude);
    void PrintOrthometricHeight(std::ostream& os, double height, const it_vstn_t& stn_it);
    void PrintEllipsoidalHeight(std::ostream& os, double height);
    void PrintStationUncertainties(std::ostream& os, const matrix_2d& var_local);

    // Stage 4: Enhanced unique stations list processing utilities
    void SortBlockStationsMapUnique();
    void ProcessBlockStaging(_it_u32u32_uint32_pair _it_bsmu, UINT32& block);
    void FinalizeBlockStaging(UINT32 block, v_uint32_string_pair& stationsOutput, std::ostream& os);

  private:
    dna_adjust& adjust_;

    // Helper functions
    constexpr int GetStationCount(char measurement_type) const;
    constexpr bool IsAngularType(char measurement_type) const;
    void PrintMeasurementRecords(const v_uint32_u32u32_pair& msr_block, bool adjustedMeasurements);
    
    // Stage 5: Enhanced measurement formatting helpers
    double CalculateAngularPrecision(const it_vmsr_t& it_msr, char cardinal) const;
    double CalculateLinearPrecision(const it_vmsr_t& it_msr, char cardinal) const;
    void FormatAngularMeasurement(const double& preAdjMeas, const double& adjMeas, const double& precision, const double& correction, const it_vmsr_t& it_msr, bool printAdjMsr);
    void FormatLinearMeasurement(const double& measurement, const double& correction, const it_vmsr_t& it_msr, bool printAdjMsr);

    // Station count lookup table
    static constexpr std::array kStationCounts{
        std::pair{'A', 3}, std::pair{'B', 2}, std::pair{'K', 2}, std::pair{'V', 2},
        std::pair{'Z', 2}, std::pair{'C', 2}, std::pair{'E', 2}, std::pair{'L', 2},
        std::pair{'M', 2}, std::pair{'S', 2}, std::pair{'H', 2}, std::pair{'R', 2},
        std::pair{'I', 2}, std::pair{'J', 2}, std::pair{'P', 2}, std::pair{'Q', 2}};

    // Angular measurement types
    static constexpr std::array kAngularTypes{'A', 'B', 'K', 'V', 'Z'};
};

// Template implementation - will be specialized for Angular/Linear types
template <typename MeasurementType>
void DynAdjustPrinter::PrintAdjMeasurements(char cardinal, const it_vmsr_t& it_msr, bool initialise_dbindex) {
    // Default implementation delegates to existing print functions
    // This will be replaced by explicit specializations
    static_assert(sizeof(MeasurementType) == 0, "Must use specialization");
}

template <typename MeasurementType>
void DynAdjustPrinter::PrintComputedMeasurements(char cardinal, const double& computed, const double& correction, const it_vmsr_t& it_msr) {
    // Default implementation - will be replaced by explicit specializations
    static_assert(sizeof(MeasurementType) == 0, "Must use specialization");
}

// Stage 4: Template implementations for station coordinate formatting
template <typename CoordinateType>
void DynAdjustPrinter::PrintStationCoordinates(std::ostream& os, const it_vstn_t& stn_it,
                                               const matrix_2d* estimates, const matrix_2d* variances) {
    // Default implementation - will be replaced by explicit specializations
    static_assert(sizeof(CoordinateType) == 0, "Must use specialization");
}

template <typename CoordinateType>
void DynAdjustPrinter::PrintStationUncertainties(std::ostream& os, const it_vstn_t& stn_it,
                                                 const matrix_2d* variances, UncertaintyMode mode) {
    // Default implementation - will be replaced by explicit specializations
    static_assert(sizeof(CoordinateType) == 0, "Must use specialization");
}

// Template specializations - defined in .cpp file
template <>
void DynAdjustPrinter::PrintAdjMeasurements<AngularMeasurement>(char cardinal, const it_vmsr_t& it_msr,
                                                                     bool initialise_dbindex);

template <>
void DynAdjustPrinter::PrintAdjMeasurements<LinearMeasurement>(char cardinal, const it_vmsr_t& it_msr,
                                                                    bool initialise_dbindex);

template <>
void DynAdjustPrinter::PrintComputedMeasurements<AngularMeasurement>(char cardinal, const double& computed, 
                                                                        const double& correction, const it_vmsr_t& it_msr);

template <>
void DynAdjustPrinter::PrintComputedMeasurements<LinearMeasurement>(char cardinal, const double& computed, 
                                                                       const double& correction, const it_vmsr_t& it_msr);

// Stage 4: Station coordinate formatting specializations
template <>
void DynAdjustPrinter::PrintStationCoordinates<GeographicCoordinates>(std::ostream& os, const it_vstn_t& stn_it,
                                                                      const matrix_2d* estimates, const matrix_2d* variances);

template <>
void DynAdjustPrinter::PrintStationCoordinates<CartesianCoordinates>(std::ostream& os, const it_vstn_t& stn_it,
                                                                     const matrix_2d* estimates, const matrix_2d* variances);

template <>
void DynAdjustPrinter::PrintStationCoordinates<ProjectionCoordinates>(std::ostream& os, const it_vstn_t& stn_it,
                                                                      const matrix_2d* estimates, const matrix_2d* variances);

template <>
void DynAdjustPrinter::PrintStationUncertainties<GeographicCoordinates>(std::ostream& os, const it_vstn_t& stn_it,
                                                                        const matrix_2d* variances, UncertaintyMode mode);

template <>
void DynAdjustPrinter::PrintStationUncertainties<CartesianCoordinates>(std::ostream& os, const it_vstn_t& stn_it,
                                                                       const matrix_2d* variances, UncertaintyMode mode);

// Stage 5: Enhanced measurement formatting template specializations
template <>
void DynAdjustPrinter::PrintMeasurementValue<AngularMeasurement>(char cardinal, const double& measurement, 
                                                                 const double& correction, const it_vmsr_t& it_msr, bool printAdjMsr);

template <>
void DynAdjustPrinter::PrintMeasurementValue<LinearMeasurement>(char cardinal, const double& measurement, 
                                                                const double& correction, const it_vmsr_t& it_msr, bool printAdjMsr);

} // namespace networkadjust
} // namespace dynadjust

#endif // DNAADJUST_PRINTER_H_
