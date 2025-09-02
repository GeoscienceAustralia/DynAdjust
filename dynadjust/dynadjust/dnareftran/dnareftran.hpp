//============================================================================
// Name         : dnareftran.hpp
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
// Description  : Reference Frame Transformation library
//============================================================================

#ifndef DNAREFTRAN_H_
#define DNAREFTRAN_H_

#if defined(_MSC_VER)
	#if defined(LIST_INCLUDES_ON_BUILD) 
		#pragma message("  " __FILE__) 
	#endif
#endif

/// \cond
#include <exception>
#include <stdexcept>
#include <memory>
#include <filesystem>
/// \endcond

#include <include/ide/trace.hpp>

#include <include/io/bst_file.hpp>
#include <include/io/bms_file.hpp>
#include <include/io/dnaiotpb.hpp>
#include <include/io/dnaiodna.hpp>

#include <include/config/dnaconsts.hpp>
#include <include/config/dnaversion.hpp>
#include <include/config/dnaexports.hpp>
#include <include/config/dnaversion-stream.hpp>
#include <include/config/dnaoptions.hpp>
#include <include/config/dnaoptions-interface.hpp>
#include <include/exception/dnaexception.hpp>

#include <include/functions/dnastrutils.hpp>
#include <include/functions/dnatemplatedatetimefuncs.hpp>
#include <include/functions/dnaiostreamfuncs.hpp>
#include <include/functions/dnastringfuncs.hpp>
#include <include/functions/dnastrmanipfuncs.hpp>
#include <include/functions/dnatemplatematrixfuncs.hpp>
#include <include/functions/dnatransparamfuncs.hpp>
#include <include/functions/dnafilepathfuncs.hpp>

#include <include/parameters/dnaepsg.hpp>
#include <include/parameters/dnaellipsoid.hpp>
#include <include/parameters/dnaprojection.hpp>
#include <include/parameters/dnadatum.hpp>
#include <include/parameters/dnatransformationparameters.hpp>
#include <include/parameters/dnaplatemotionmodels.hpp>
#include <include/parameters/dnaframesubstitutions.hpp>

#include <include/math/dnamatrix_contiguous.hpp>
#include <include/measurement_types/dnameasurement.hpp>
#include <include/measurement_types/dnageometries.hpp>

using namespace dynadjust::measurements;
using namespace dynadjust::exception;
using namespace dynadjust::iostreams;
using namespace dynadjust::epsg;
using namespace dynadjust::datum_parameters;
using namespace dynadjust::pmm_parameters;
using namespace dynadjust::math;
using namespace dynadjust::geometries;
using namespace dynadjust::frame_substitutions;

namespace dynadjust {
namespace referenceframe {

// This class is exported from the dnaReftran.dll
#ifdef _MSC_VER
class DNAREFTRAN_API dna_reftran {
#else
class dna_reftran {
#endif

public:
	dna_reftran();
	dna_reftran(const project_settings& p, std::ofstream* f_out);
	virtual ~dna_reftran();

private:
	// Disallow use of compiler generated functions. See dnaadjust.hpp
	dna_reftran(const dna_reftran&);
	dna_reftran& operator=(const dna_reftran&);	

public:
	void TransformBinaryFiles(const std::string& bstFile, const std::string& bmsFile, const std::string& newFrame, const std::string& newEpoch="");
	
	// Returns the file progress
	//inline int ReturnFileProgress() const { return (int)m_dPercentComplete; }

	// Returns the byte offset
	inline int GetByteOffset() const { return m_iBytesRead; }

	// Sets the byte offset
	inline void SetByteOffset() { m_iBytesRead = 0; }
	inline UINT32 StationsTransformed() const { return m_stnsTransformed; }
	inline UINT32 StationsNotTransformed() const { return m_stnsNotTransformed; }
	inline UINT32 MeasurementsTransformed() const { return m_msrsTransformed; }
	inline UINT32 MeasurementsNotTransformed() const { return m_msrsNotTransformed; }

	// file handling
	void SerialiseDNA(const std::string& stnfilename, const std::string& msrfilename, bool flagUnused=false);
	void SerialiseDynaML(const std::string& xmlfilename, bool flagUnused=false);
	void SerialiseDynaML(const std::string& stnfilename, const std::string& msrfilename, bool flagUnused=false);

	void SerialiseDynaMLStn(std::ofstream* xml_file, CDnaProjection& projection, bool flagUnused=false);
	void SerialiseDynaMLMsr(std::ofstream* xml_file);

	//bool PrintTransformedStationCoordinatestoSNX();

	void LoadTectonicPlateParameters(const std::string& pltfileName, const std::string& pmmfileName);

	//void LoadFrameSubstitutions(const std::string& frxfileName);
	void LoadWGS84FrameSubstitutions();

	inline void InitialiseSettings(const project_settings& p) {projectSettings_ = p;}

private:

	void LoadBinaryStationFile(const std::string& bstfileName);
	void LoadBinaryMeasurementFile(const std::string& bmsfileName);

	void WriteBinaryStationFile(const std::string& bstfileName);
	void WriteBinaryMeasurementFile(const std::string& bmsfileName);

	void TransformStationRecords(const std::string& newFrame, const std::string& newEpoch);
	void TransformMeasurementRecords(const std::string& newFrame, const std::string& newEpoch);

	double DetermineElapsedTime(const CDnaDatum& datumFrom, const CDnaDatum& datumTo,
		transformation_parameter_set& transParams, transformationType transType);
	
	void ObtainHelmertParameters(const CDnaDatum& datumFrom, const CDnaDatum& datumTo, 
		transformation_parameter_set& transParams, double& timeElapsed, transformationType transType);
	
	UINT32 DetermineTectonicPlate(const std::string& plate);

	void ObtainPlateMotionParameters(it_vstn_t& stn_it, double* reduced_parameters, const CDnaDatum& datumFrom, const CDnaDatum& datumTo,
		transformation_parameter_set& transformParameters, double& timeElapsed);

	void JoinTransformationParameters(it_vstn_t& stn_it, double* reduced_parameters, const CDnaDatum& datumFrom, const CDnaDatum& datumTo,
		transformation_parameter_set& transformParameters, transformationType transType, const matrix_2d& coordinates);

	void Transform(it_vstn_t& stn_it, const matrix_2d& coordinates, matrix_2d& coordinates_mod, 
		const CDnaDatum& datumFrom, transformation_parameter_set& transformParameters);

	void TransformDynamic(it_vstn_t& stn_it, const matrix_2d& coordinates, matrix_2d& coordinates_mod,
		const CDnaDatum& datumFrom, const CDnaDatum& datumTo, transformation_parameter_set& transformParameters,
		transformationType transType);

	void TransformFrames_Join(it_vstn_t& stn_it, const matrix_2d& coordinates, matrix_2d& coordinates_mod,
		const CDnaDatum& datumFrom, const CDnaDatum& datumTo, transformation_parameter_set& transformParameters,
		transformationType transType);

	void TransformFrames_PlateMotionModel(it_vstn_t& stn_it, const matrix_2d& coordinates, matrix_2d& coordinates_mod,
		const CDnaDatum& datumFrom, const CDnaDatum& datumTo, transformation_parameter_set& transformParameters);

	void TransformFrames_WithoutPlateMotionModel(it_vstn_t& stn_it, const matrix_2d& coordinates, matrix_2d& coordinates_mod,
		const CDnaDatum& datumFrom, const CDnaDatum& datumTo, transformation_parameter_set& transformParameters,
		transformationType transType);

	void TransformEpochs_PlateMotionModel(it_vstn_t& stn_it, const matrix_2d& coordinates, matrix_2d& coordinates_mod,
		const CDnaDatum& datumFrom, const CDnaDatum& datumTo);

	void TransformStation(it_vstn_t& stn_it, const CDnaDatum& datumFrom,
		transformation_parameter_set& transformParameters);
	
	void TransformMeasurement(it_vmsr_t& msr_it, const CDnaDatum& datumFrom,
		transformation_parameter_set& transformParameters);
	void TransformMeasurement_GX(it_vmsr_t& msr_it, const CDnaDatum& datumFrom, 
		transformation_parameter_set& transformParameters);
	void TransformMeasurement_Y(it_vmsr_t& msr_it, const CDnaDatum& datumFrom, 
		transformation_parameter_set& transformParameters);

	void IdentifyStationPlate();
	
	void CalculateRotations();

	void LogFrameSubstitutions(std::vector<string_string_pair>& substitutions, const std::string& type);
	void ApplyStationFrameSubstitutions();
	void ApplyMeasurementFrameSubstitutions();

	bool IsolateandApplySubstitute(const std::string& epsgCode, const std::string& epoch, std::string& epsgSubstitute);

	//double							m_dPercentComplete;			// percentage of bytes read from file
	int								m_iBytesRead;				// bytes read from file
	UINT32							m_stnsTransformed;
	UINT32							m_stnsNotTransformed;
	UINT32							m_msrsTransformed;
	UINT32							m_msrsNotTransformed;

	binary_file_meta_t				bst_meta_;
	binary_file_meta_t				bms_meta_;

	vstn_t							bstBinaryRecords_;
	vmsr_t							bmsBinaryRecords_;
	
	project_settings				projectSettings_;
	CDnaDatum						datumTo_;
	CDnaDatum						datumFrom_;

	bool							transformationPerformed_;

	v_string_v_doubledouble_pair	global_plates_;				// Tectonic plate boundaries
	v_plate_motion_eulers			plate_motion_eulers_;		// Euler parameters corresponding to each plate
	v_plate_motion_cartesians		plate_motion_cartesians_;	// Helmert parameters computed from Euler parameters

	vframeSubsPtr					_frameSubstitutions;		// Reference frame substitutions
	std::vector<string_string_pair>		_v_stn_substitutions;		// station substitutions made
	std::vector<string_string_pair>		_v_msr_substitutions;		// station substitutions made

	v_string_uint32_pair 			vplateMap_;					// Plate Map index sorted on plate ID

	std::ofstream* 					rft_file;
	_INPUT_DATA_TYPE_ 				data_type_;					// stations or measurements

	// Database management
	v_msr_database_id_map	v_msr_db_map_;
	bool					databaseIDsLoaded_;
	bool					databaseIDsSet_;

	void LoadDatabaseId();
};

}	// namespace referenceframe
}	// namespace dynadjust


#endif /* DNAREFTRAN_H_ */
