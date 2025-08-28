//============================================================================
// Name         : dnaplot.hpp
// Author       : Roger Fraser
// Contributors : Dale Roberts <dale.o.roberts@gmail.com>
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
// Description  : DynAdjust Plot library
//============================================================================

#ifndef DNAPLOT_H_
#define DNAPLOT_H_

#if defined(_MSC_VER)
	#if defined(LIST_INCLUDES_ON_BUILD) 
		#pragma message("  " __FILE__) 
	#endif
#endif

/// \cond
#include <exception>
#include <stdexcept>
#include <iostream>
#include <string>
#include <utility>
#include <vector>
#include <map>
#include <cstdarg>
#include <math.h>

#include <memory>
#include <filesystem>
#include <thread>

#include <include/config/dnaexports.hpp>
#include <include/config/dnaversion.hpp>
#include <include/config/dnatypes-gui.hpp>

#include <include/exception/dnaexception.hpp>

#include <include/functions/dnatemplatefuncs.hpp>
#include <include/functions/dnatemplatestnmsrfuncs.hpp>
#include <include/functions/dnaiostreamfuncs.hpp>
#include <include/functions/dnastringfuncs.hpp>
#include <include/functions/dnaintegermanipfuncs.hpp>
#include <include/functions/dnafilepathfuncs.hpp>
#include <include/functions/dnatemplategeodesyfuncs.hpp>
#include <include/functions/dnatemplatematrixfuncs.hpp>
#include <include/functions/dnastrmanipfuncs.hpp>
#include <include/functions/dnaprocessfuncs.hpp>

#include <include/config/dnaconsts.hpp>
#include <include/config/dnatypes.hpp>

#include <include/parameters/dnaepsg.hpp>
#include <include/parameters/dnadatum.hpp>
#include <include/parameters/dnaprojection.hpp>

#include <include/measurement_types/dnadirection.hpp>
#include <include/measurement_types/dnagpsbaseline.hpp>

#include <include/math/dnamatrix_contiguous.hpp>

#include <include/io/bst_file.hpp>
#include <include/io/bms_file.hpp>
#include <include/io/map_file.hpp>
#include <include/io/seg_file.hpp>
#include <include/io/dnaiotpb.hpp>

using namespace dynadjust::measurements;
using namespace dynadjust::epsg;
using namespace dynadjust::exception;
using namespace dynadjust::datum_parameters;
using namespace dynadjust::iostreams;
using namespace dynadjust::math;

namespace dynadjust {
namespace networkplot {


#ifdef _MSC_VER
class DNAPLOT_API dna_plot {
#else
class dna_plot {
#endif

public:
	dna_plot();
	virtual ~dna_plot();

private:
	// Disallow use of compiler generated functions. See dnaadjust.hpp
	dna_plot(const dna_plot&);
	dna_plot& operator=(const dna_plot&);	

public:
	void LoadNetworkFiles(const project_settings& projectSettings);
	
	// GMT maps
	void CreateGMTPlotEnvironment(project_settings* pprj);
	void CreateGMTPlot();

	// gnuplot graphs
	void CreategnuplotGraphEnvironment(project_settings* pprj, const plotGraphMode& graphMode);
	void CreateSegmentationGraph(const plotGraphMode& graphMode);

	inline UINT32 blockCount() const { return blockCount_; }
	inline UINT32 stationCount() const { return stationCount_; }
	inline UINT32 measurementCount() const { return measurementCount_; }
	inline UINT32 measurementCategories() const { return measurementCategories_; }
	inline UINT32 blockThreshold() const { return blockThreshold_; }
	inline UINT32 minInnerStns() const { return minInnerStns_; }

private:

	CDnaDatum			datum_;
	CDnaProjection		projection_;
	
	// Gnuplot plotting methods
	void PlotGnuplotDatFileStns();
	void PlotGnuplotDatFileMsrs();
	void PrintGnuplotCommandFile(const std::string& gnuplot_cmd_file, const plotGraphMode& graphMode);
	void PrintGnuplotCommandFileStns(const UINT32& fontSize);
	void PrintGnuplotCommandFileMsrs(const UINT32& fontSize);

	void LoadBinaryFiles();
	void LoadStationMap();
	void LoadSegmentationFile();
	void LoadPosUncertaintyFile();

	void InitialiseAppsandSystemCommands();	

	void InitialiseGMTParameters();
	void FinaliseGMTParameters();
	void InitialiseGMTFilenames();
	void CreateGMTInputFiles();
	void CreateExtraInputFiles();
	void CreateGMTCommandFiles();
	void CreateGMTCommandFile(const UINT32& block);
	bool InitialiseandValidateStartingBlock(UINT32& block);

	void InvokeGnuplot();
	void CleanupGnuplotFiles(const plotGraphMode& graphMode);

	void InvokeGMT();
	void AggregateGMTPDFs();
	void CleanupGMTFiles();

	void ComputeStationCorrections();
	void ComputeStationCorrection(it_vstn_t_const _it_stn, stationCorrections_t& stnCor,
		matrix_2d& currentEstimates, matrix_2d& initialEstimates);
	void LoadCorrectionsFile();
	
	//void SortandMapStations();
	void ComputeMeasurementCount();

	// GMT plotting methods
	void FormGMTDataFileNames(const UINT32& block=0);
	void PrintGMTParameters();
	
	void PrintStationDataFile(std::ostream& os, it_vstn_t_const _it_stn);
	void PrintStationsDataFile();
	void PrintStationsDataFileBlock(const UINT32& block);
	
	void PrintStationLabel(std::ostream& os, it_vstn_t_const _it_stn);
	void PrintStationLabels();
	void PrintStationLabelsBlock(const UINT32& block);
	
	void PrintPlateBoundaries(const UINT32& block);
	void PrintCorrectionArrows(const UINT32& block);
	void PrintPositionalUncertainty(const UINT32& block);
	void PrintErrorEllipses(const UINT32& block);
	
	void PrintMeasurementsDatFiles();
	void PrintMeasurementsDatFilesBlock(const UINT32& block);

	void PrintMeasurementsDatFile(char msrType, std::ofstream* msr_file_stream);
	void PrintMeasurementsDatFileBlock(const UINT32& block, char msrType, std::ofstream* msr_file_stream);
	
	//void PrintGMTPlotCoords(vdnaStnPtr* vStations, vdnaMsrPtr* vMeasurements, ostream& osStn, ostream& osStn2, ostream& osMsr); 
	//void PrintGMTPlotCoords_D(vdnaStnPtr* vStations, std::vector<CDnaDirection>* vDirections, ostream& osMsr); 
	//void PrintGMTPlotCoords_GX(vdnaStnPtr* vStations, std::vector<CDnaGpsBaseline>* vBaselines, ostream& osMsr); 

	void CalculateLimitsFromStation();
	void CalculateLimitsFromPoint();
	void DetermineBoundingBox();

	void SignalExceptionPlot(const std::string& msg, const int& line_no, const char *streamType, ...);

	void ReduceStationCoordinates(station_t* stationRecord);

	void NormaliseScale(double& scale);
	void NormaliseScaleBar(double& scale_bar_width);
	void NormaliseGraticule(double& graticule_width, UINT32& graticule_width_precision);
	void SelectCoastlineResolution(const double& dDimension, std::string& coastResolution, plot_settings* plotCriteria);

	bool WithinLimits(const double& latitude, const double& longitude);

	void BuildParameterStationList();

	vvUINT32				v_parameterStationList_;

	std::ofstream			gmtbat_file_;
	std::ofstream			gnuplotbat_file_;

	MsrTally				parsemsrTally_;		// total network tally
	std::vector<MsrTally>		v_msr_tally_;		// per block tally

	vstn_t					bstBinaryRecords_;
	vmsr_t					bmsBinaryRecords_;
	binary_file_meta_t		bst_meta_;
	binary_file_meta_t		bms_meta_;
	std::string					output_folder_;
	std::string					network_name_;
	v_string_uint32_pair	stnsMap_;		// Station Name Map sorted on name (string)
	_PLOT_STATUS_			plotStatus_;
	project_settings		projectSettings_;
	project_settings*		pprj_;
	std::string					reference_frame_;
	
	double					lowerDeg_;
	double					leftDeg_;
	double					upperDeg_;
	double					rightDeg_;
	double					centre_width_;
	double					centre_height_;

	double					mapScale_;
	double					arrow_legend_long_;
	double					arrow_legend_lat_;
	double					average_correction_;
	double					arrow_legend_length_;

	double					uncertainty_legend_long_;
	double					uncertainty_legend_lat_;
	double					largest_uncertainty_;
	double					uncertainty_legend_length_;

	double 					graticule_width_;
	UINT32 					graticule_width_precision_;
	
	_it_pair_vchar			it_vchar_range;
	it_pair_string_vUINT32	it_stnmap_range;
	vchar					_combined_msr_list;

	double					default_paper_width_;
	double					default_paper_height_;

	std::string					seg_stn_graph_file_;
	std::string					seg_msr_graph_file_;

	vstring					v_isl_pts_file_;
	vstring					v_isl_const_file_;
	vstring					v_isl_lbl_file_;
	vstring					v_tectonic_plate_file_;
	vstring					v_stn_cor_file_;
	vstring					v_stn_apu_file_;
	vstring					v_stn_err_file_;
	vstring					v_jsl_pts_file_;
	vstring					v_jsl_const_file_;
	vstring					v_jsl_lbl_file_;
	vstring					v_msr_def_file_;
	vstring					v_msr_all_file_;

	vstring					v_gmt_cmd_filenames_;
	vstring					v_gmt_pdf_filenames_;
	
	vv_string_string_pair	v_msr_file_;

	vstnCor_t				v_stn_corrs_;
	vstnPU_t				v_stn_pu_;

	bool					isLandscape;
	bool					default_limits_;
	bool					plotConstraints_;

	// segmentation block variables
	bool					plotBlocks_;
	UINT32					stationCount_;
	UINT32					measurementCount_;
	UINT32					measurementCategories_;
	UINT32					blockCount_;
	UINT32					blockThreshold_;
	UINT32					minInnerStns_;
	vUINT32					YClusterStations_;
	vvUINT32				v_ISL_;				// Inner stations
	vvUINT32				v_JSL_;				// Junction stations
	vvUINT32				v_CML_;				// Measurements

	vUINT32					v_measurementCount_;		// number of raw measurements and constrained stations
	vUINT32					v_unknownsCount_;			// number of all stations (constrained and free)

	std::string					_APP_GMTSET_;
	std::string					_APP_PSCOAST_;
	std::string					_APP_PSCONVERT_;
	std::string					_APP_PSTEXT_;
	std::string					_APP_PSVELO_;
	std::string					_APP_PSXY_;
	std::string					_APP_PSLEGEND_;
	std::string					_PDF_AGGREGATE_;
	std::string					_CMD_EXT_;
	std::string 					_ECHO_CMD_;
	std::string 					_LEGEND_ECHO_;
	std::string 					_LEGEND_CMD_1_;
	std::string 					_LEGEND_CMD_2_;

	std::string 					_COMMENT_PREFIX_;

	std::string					_CMD_HEADER_;
	std::string 					_DELETE_CMD_;
	std::string 					_CHMOD_CMD_;
	std::string 					_COPY_CMD_;
	std::string 					_MOVE_CMD_;
	std::string					_NULL_OUTPUT_;
	std::string 					_MAKEDIR_CMD_;
	std::string 					_RMDIR_CMD_;
	std::string 					_MAKEENV_CMD_;
	std::string 					_MAKETEMP_CMD_;
	std::string 					_ENV_GMT_TMP_DIR_;
	std::string 					_GMT_TMP_DIR_;

	double					dWidth_;
	double					dHeight_;
	double 					dDimension_;
	double 					dBuffer_;
	double 					dScaleLat_;
	double 					title_block_height_;
	double 					page_width_;
	double 					avg_data_scale_;
	double 					scale_;
	double 					scale_bar_width_;
	int 					scale_precision_;

	double 					line_width_;
	double 					circle_radius_;
	double 					circle_radius_2_;
	double 					circle_line_width_;

	std::string 					coastResolution_;
	
};


}	// namespace mathcomp
}	// namespace dynadjust



#endif /* DNAPLOT_H_ */
