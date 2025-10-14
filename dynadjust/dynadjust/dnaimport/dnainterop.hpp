//============================================================================
// Name         : dnainterop.hpp
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
// Description  : DynAdjust Interoperability library
//============================================================================

#ifndef DNAINTEROP_H_
#define DNAINTEROP_H_

#if defined(_MSC_VER)
#if defined(LIST_INCLUDES_ON_BUILD) 
#pragma message("  " __FILE__) 
#endif
#endif

/// \cond
#include <system_error>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdarg>
#include <math.h>
#include <mutex>

#include <memory>
#include <filesystem>
/// \endcond

#include <include/config/dnaexports.hpp>
#include <include/config/dnaversion.hpp>
#include <include/config/dnaconsts.hpp>
#include <include/config/dnatypes.hpp>
#include <include/config/dnatypes-gui.hpp>
#include <include/exception/dnaexception.hpp>

#include <include/parameters/dnaepsg.hpp>
#include <include/parameters/dnadatum.hpp>

#include <dynadjust/dnaimport/dnaparser_pimpl.hxx>

#include <include/io/dnaiodna.hpp>
#include <include/io/bst_file.hpp>
#include <include/io/bms_file.hpp>
#include <include/io/aml_file.hpp>
#include <include/io/asl_file.hpp>
#include <include/io/map_file.hpp>
#include <include/io/seg_file.hpp>
#include <include/io/dnaiosnx.hpp>
#include <include/io/dnaioscalar.hpp>

#include <include/functions/dnafilepathfuncs.hpp>
#include <include/functions/dnastrutils.hpp>
#include <include/functions/dnatemplatefuncs.hpp>
#include <include/functions/dnatemplatestnmsrfuncs.hpp>
#include <include/functions/dnaiostreamfuncs.hpp>
#include <include/functions/dnastringfuncs.hpp>

using namespace dynadjust::measurements;
using namespace dynadjust::math;
using namespace dynadjust::exception;
using namespace dynadjust::datum_parameters;
using namespace dynadjust::iostreams;

extern UINT32 g_fileOrder;

namespace dynadjust {
namespace dynamlinterop {

#ifdef _MSC_VER
class DNAIMPORT_API dna_import {
#else
class dna_import {
#endif
public:
	dna_import();
	virtual ~dna_import();

private:
	// Disallow use of compiler generated functions. See dnaadjust.hpp
	dna_import(const dna_import&);	
	dna_import& operator=(const dna_import&);

public:
	// Parse an xml file
	_PARSE_STATUS_ ParseInputFile(const std::string& filename, vdnaStnPtr* vStations, PUINT32 stnCount, 
		vdnaMsrPtr* vMeasurements, PUINT32 msrCount, 
		PUINT32 clusterID, input_file_meta_t* input_file_meta, bool firstFILE,
		std::string* success_msg, project_settings* p);

	_PARSE_STATUS_ LoadDNAGeoidFile(const std::string& fileName, vdnaStnPtr* vStations);

	void RemoveIgnoredMeasurements(vdnaMsrPtr* vMeasurements, MsrTally* parsemsrTally);
	void IncludeMeasurementTypes(const std::string& includeMsrs, vdnaMsrPtr* vMeasurements, MsrTally* parsemsrTally);
	void ExcludeMeasurementTypes(const std::string& excludeMsrs, vdnaMsrPtr* vMeasurements, MsrTally* parsemsrTally);
	void IgnoreInsufficientMeasurements(vdnaStnPtr* vStations, vdnaMsrPtr* vMeasurements, pvstring vPoorlyConstrainedStns);

	void ExcludeAllOutsideBoundingBox(vdnaStnPtr* vStations, vdnaMsrPtr* vMeasurements, 
		StnTally* parsestnTally, MsrTally* parsemsrTally, pvstring pvUnusedStns, const project_settings& p,
		bool& splitXmsrs, bool& splitYmsrs);
	void ExtractStnsAndAssociatedMsrs(const std::string& stnListInclude, const std::string& stnListExclude, vdnaStnPtr* vStations, vdnaMsrPtr* vMeasurements, 
		StnTally* parsestnTally, MsrTally* parsemsrTally, pvstring vExcludedStns, const project_settings& p,
		bool& splitXmsrs, bool& splitYmsrs);
	void ImportStnsMsrsFromBlock(vdnaStnPtr* vStations, vdnaMsrPtr* vMeasurements, const project_settings& p);
	void ImportStnsMsrsFromNetwork(vdnaStnPtr* vStations, vdnaMsrPtr* vMeasurements, const project_settings& p);

	UINT32 RemoveDuplicateStations(vdnaStnPtr* vStations, pvstring vduplicateStations, pv_stringstring_doubledouble_pair vnearbyStations);
	UINT32 FindSimilarMeasurements(vdnaMsrPtr* vMeasurements, vdnaMsrPtr* vSimilarMeasurements);
	UINT32 FindSimilarGXMeasurements(vdnaMsrPtr* vMeasurements, vdnaMsrPtr* vSimilarMeasurements);
	void FullSortandMapStations(vdnaStnPtr* vStations, pv_string_uint32_pair vStnsMap_sortName);
	//void SortandMapStations(vdnaStnPtr* vStations, pv_string_uint32_pair vStnsMap_sortName);
	void SortStationsForExport(vdnaStnPtr* vStations);
	void ReduceStations(vdnaStnPtr* vStations, const CDnaProjection& projection);
	void ReduceStations(vdnaStnPtr* vStations, const CDnaProjection& projection, const UINT32& cores);
	void RenameStations(vdnaStnPtr* vStations, vdnaMsrPtr* vMeasurements, project_settings* p);
	void ApplyDiscontinuities(vdnaMsrPtr* vMeasurements);
	void TrackDiscontinuitySite(const std::string& site, const std::string& site_renamed);
	void ApplyDiscontinuitiesMeasurements(vdnaMsrPtr* vMeasurements);

	void ApplyDiscontinuitiesMeasurements_GX(std::vector<CDnaGpsBaseline>* vGpsBaselines);
	void ApplyDiscontinuitiesMeasurements_Y(std::vector<CDnaGpsPoint>* vGpsPoints);
	void ApplyDiscontinuitiesMeasurements_D(std::vector<CDnaDirection>* vDirections, const boost::gregorian::date& site_date);
	
	void AddDiscontinuityStations(vdnaStnPtr* vstationsTotal);

	void EditGNSSMsrScalars(vdnaMsrPtr* vMeasurements, project_settings* p);
	void ApplyGNSSMsrScalar(std::vector<CDnaGpsBaseline>::iterator& _it_bsl, vscl_t& bslScalars);

	void MapMeasurementStations(vdnaMsrPtr* vMeasurements, pvASLPtr vAssocStnList, PUINT32 lMapCount, pvstring vUnusedStns, pvUINT32 vIgnoredMsrs);
	UINT32 ComputeMeasurementCount(vdnaMsrPtr* vMeasurements, const vUINT32& vIgnoredMsrs);

	// file handling
	void SerialiseDNA(vdnaStnPtr* vStations, vdnaMsrPtr* vMeasurements, const std::string& stnfilename, const std::string& msrfilename, const project_settings& p, vifm_t* vinput_file_meta, bool flagUnused=false);
	void SerialiseMSR(vdnaMsrPtr* vMeasurements, const std::string& msrfilename, const project_settings& p, vifm_t* vinput_file_meta);
	void SerialiseSTN(vdnaStnPtr* vStations, const std::string& stnfilename, const project_settings& p, vifm_t* vinput_file_meta, bool flagUnused = false);
	//void SerialiseDynaMLfromBinary(const std::string& outfilename, const project_settings& p, vifm_t* vinput_file_meta, bool flagUnused=false);
	void SerialiseDynaMLfromMemory(vdnaStnPtr* vStations, vdnaMsrPtr* vMeasurements, const std::string& outfilename, const project_settings& p, vifm_t* vinput_file_meta, bool flagUnused=false);
	//void SerialiseDynaMLSepfromBinary(const std::string& stnfilename, const std::string& msrfilename, const project_settings& p, vifm_t* vinput_file_meta, bool flagUnused=false);
	void SerialiseDynaMLSepfromMemory(vdnaStnPtr* vStations, vdnaMsrPtr* vMeasurements, const std::string& stnfilename, const std::string& msrfilename, const project_settings& p, vifm_t* vinput_file_meta, bool flagUnused=false);
	//void SerialiseGeoidData(vdnaStnPtr* vStations, const std::string& geofilename);
	void SimulateMSR(vdnaStnPtr* vStations, vdnaMsrPtr* vMeasurements, const std::string& msrfilename, const project_settings& p);

	void SerialiseBst(const std::string& bst_filename, vdnaStnPtr* vStations, pvstring vUnusedStns, vifm_t& vinput_file_meta, bool flagUnused=false);
	void SerialiseBms(const std::string& bms_filename, vdnaMsrPtr* vMeasurements, vifm_t& vinput_file_meta);
	void SerialiseMap(const std::string& stnmap_file);
	void SerialiseAml(const std::string& filename, pvUINT32 vAssocMsrList);
	void SerialiseAsl(const std::string& filename, pvASLPtr vAssocStnList);
	void SerialiseMapTextFile(const std::string& stnmap_file);
	void SerialiseDiscontTextFile(const std::string& discont_file);
	void SerialiseAmlTextFile(const std::string& bms_filename, const std::string& aml_filename, pvUINT32 vAML, pvASLPtr vAssocStnList, vdnaStnPtr* vStations);
	void SerialiseAslTextFile(const std::string& filename, pvASLPtr vAssocStnList, vdnaStnPtr* vStations);

	void SerialiseDatabaseId(const std::string& dbid_filename, pvdnaMsrPtr vMeasurements);

	void CompleteAssociationLists(vdnaMsrPtr* vMeasurements, pvASLPtr vAssocStnList, pvUINT32 vAssocMsrList, const _AML_TYPE_ aml_type = str_msr);

	double GetProgress();
	inline bool IsProcessing() const { return isProcessing_; }
	inline MsrTally& GetMsrTally() const { return *p_parsemsr_tally; }
	inline StnTally& GetStnTally() const { return *p_parsestn_tally; }

	inline _PARSE_STATUS_ GetStatus() const { return parseStatus_; }

	inline void ResetFileOrder() const { g_fileOrder = 0; }
	inline bool filespecifiedReferenceFrame() const { return _filespecifiedreferenceframe; }
	inline bool filespecifiedEpoch() const { return _filespecifiedepoch; }
	void InitialiseDatum(const std::string& reference_frame, const std::string epoch="");
	
	void PrintMeasurementsToStations(std::string& m2s_file, MsrTally* parsemsrTally,
		std::string& bst_file, std::string& bms_file, std::string& aml_file, pvASLPtr vAssocStnList);

	// Discontonuity file
	void ParseDiscontinuities(const std::string& fileName);

private:

	// DynaML files
	bool CheckAndWarnInvalidSchemaLocation(const std::string& fileName);

	void ParseXML(const std::string& fileName, vdnaStnPtr* vStations, PUINT32 stnCount,
							   vdnaMsrPtr* vMeasurements, PUINT32 msrCount, PUINT32 clusterID,
							   std::string& fileEpsg, std::string& fileEpoch, bool firstFile, std::string* success_msg);
	
	// SINEX files
	void ParseSNX(const std::string& fileName, vdnaStnPtr* vStations, PUINT32 stnCount, 
							   vdnaMsrPtr* vMeasurements, PUINT32 msrCount, PUINT32 clusterID);
	
	// DNA Ascii files
	//void ParseDNAVersion(const INPUT_DATA_TYPE& idt);
	void ParseDNA(const std::string& fileName, vdnaStnPtr* vStations, PUINT32 stnCount, 
							   vdnaMsrPtr* vMeasurements, PUINT32 msrCount, PUINT32 clusterID, 
							   std::string& fileEpsg, std::string& fileEpoch, bool firstFile);
	void ParseDNASTN(vdnaStnPtr* vStations, PUINT32 stnCount,
								const std::string& fileEpsg, const std::string& fileEpoch);
	void ParseDNAMSR(pvdnaMsrPtr vMeasurements, PUINT32 msrCount, PUINT32 clusterID,
								const std::string& fileEpsg, const std::string& fileEpoch);

	//void SetDefaultReferenceFrame(vdnaStnPtr* vStations, vdnaMsrPtr* vMeasurements);

	void ParseDNAMSRAngular(const std::string& sBuf, dnaMsrPtr& msr_ptr);
	void ParseDNAMSRCoordinate(const std::string& sBuf, dnaMsrPtr& msr_ptr);
	UINT32 ParseDNAMSRDirections(std::string& sBuf, dnaMsrPtr& msr_ptr, bool ignoreMsr);
	void ParseDNAMSRGPSBaselines(std::string& sBuf, dnaMsrPtr& msr_ptr, bool ignoreMsr);
	void ParseDNAMSRGPSPoints(std::string& sBuf, dnaMsrPtr& msr_ptr, bool ignoreMsr);
	void ParseDNAMSRLinear(const std::string& sBuf, dnaMsrPtr& msr_ptr);
	void ParseDNAMSRCovariance(CDnaCovariance& cov);

	std::string ParseAngularValue(const std::string& sBuf, const std::string& calling_function);
	std::string ParseLinearValue(const std::string& sBuf, const std::string& msrName, const std::string& calling_function);
	std::string ParseGPSMsrValue(const std::string& sBuf, const std::string& element, const std::string& calling_function);
	std::string ParseGPSVarValue(const std::string& sBuf, const std::string& element, const UINT32 location, const UINT32 width, const std::string& calling_function);
	std::string ParseInstrumentValue(const std::string& sBuf, const std::string& calling_function);
	std::string ParseTargetValue(const std::string& sBuf, const std::string& calling_function);
	std::string ParseTarget2Value(const std::string& sBuf, const std::string& calling_function);
	std::string ParseStdDevValue(const std::string& sBuf, const std::string& calling_function);
	std::string ParseInstHeightValue(const std::string& sBuf, const std::string& calling_function);
	std::string ParseTargHeightValue(const std::string& sBuf, const std::string& calling_function);
	std::string ParseMsrCountValue(const std::string& sBuf, UINT32& msrCount, const std::string& calling_function);
	std::string ParseScaleVValue(const std::string& sBuf, const std::string& calling_function);
	std::string ParseScalePValue(const std::string& sBuf, const std::string& calling_function);
	std::string ParseScaleLValue(const std::string& sBuf, const std::string& calling_function);
	std::string ParseScaleHValue(const std::string& sBuf, const std::string& calling_function);
	std::string ParseRefFrameValue(const std::string& sBuf, const std::string& calling_function);
	std::string ParseEpochValue(const std::string& sBuf, const std::string& calling_function);

	void ParseDatabaseIds(const std::string& sBuf, const std::string& calling_function, const char msrType);
	void ParseDatabaseClusterId(const std::string& sBuf, const std::string& calling_function);
	void ParseDatabaseMsrId(const std::string& sBuf, const std::string& calling_function);

	void LoadNetworkFiles(pvstn_t binaryStn, pvmsr_t binaryMsr, const project_settings& projectSettings, bool loadSegmentFile);
	void LoadBinaryFiles(pvstn_t binaryStn, pvmsr_t binaryMsr);
	void LoadSegmentationFile(pvmsr_t binaryMsr);

	void SplitClusterMsrsConnectedToStns(vdnaMsrPtr* vClusterMsrs, 
		pvstring pvIncludedStns, pvstring pvExcludedStns, 
		bool& splitXmsrs, bool& splitYmsrs);
	void SplitClusterMsrs(vdnaMsrPtr& msrsConnectedToStns, 
		pvstring pvIncludedStns, pvstring pvExcludedStns, 
		vdnaMsrPtr* vMeasurements, bool& splitXmsrs, bool& splitYmsrs);
	void ExtractAssociatedMsrsConnectedToStns(vdnaMsrPtr* vMeasurements, MsrTally* parsemsrTally, 
		pvstring pvUsedStns, pvstring pvUnusedStns, const project_settings& p,
		bool& splitXmsrs, bool& splitYmsrs);
	void ExtractAssociatedMsrsBoundingBox(vdnaMsrPtr* vMeasurements, MsrTally* parsemsrTally, 
		pvstring pvUsedStns, pvstring pvUnusedStns, const project_settings& p,
		bool& splitXmsrs, bool& splitYmsrs);
	void ExtractAssociatedStns(vdnaMsrPtr* vMeasurements, pvstring pvUsedStns);
	void ExtractAssociatedStns_GX(std::vector<CDnaGpsBaseline>* vGpsBaselines, pvstring pvUsedStns);
	void ExtractAssociatedStns_Y(std::vector<CDnaGpsPoint>* vGpsPoints, pvstring pvUsedStns);
	void ExtractAssociatedStns_D(std::vector<CDnaDirection>* vDirections, pvstring pvUsedStns);

	void SortandMapStations(vdnaStnPtr* vStations);

	void MapMeasurementStationsBsl(std::vector<CDnaGpsBaseline>* vGpsBaselines,
		pvASLPtr vAssocStnList, PUINT32 lMapCount);
	void MapMeasurementStationsPnt(std::vector<CDnaGpsPoint>* vGpsPoints,
		pvASLPtr vAssocStnList, PUINT32 lMapCount);
	void MapMeasurementStationsDir(std::vector<CDnaDirection>* vDirections,
		pvASLPtr vAssocStnList, PUINT32 lMapCount);
	
	//void RenameStationsMsr(CDnaMeasurement* msr, v_string_string_pair& variables);
	void RenameStationsBsl(std::vector<CDnaGpsBaseline>* vGpsBaselines, v_string_vstring_pair& stnRenaming);
	void RenameStationsPnt(std::vector<CDnaGpsPoint>* vGpsPoints, v_string_vstring_pair& stnRenaming);
	void RenameStationsDir(std::vector<CDnaDirection>* vDirections, v_string_vstring_pair& stnRenaming);
	
	void CompleteASLDirections(_it_vdnamsrptr _it_msr, std::vector<CDnaDirection>* vDirections, pvASLPtr vAssocStnList, 
		pvUINT32 vAssocMsrList, PUINT32 currentBmsFileIndex, const _AML_TYPE_ aml_type);
	void CompleteASLGpsBaselines(std::vector<CDnaGpsBaseline>* vGpsBaselines, pvASLPtr vAssocStnList, 
		pvUINT32 vAssocMsrList, PUINT32 currentBmsFileIndex, const _AML_TYPE_ aml_type);
	void CompleteASLGpsPoints(std::vector<CDnaGpsPoint>* vGpsPoints, pvASLPtr vAssocStnList, 
		pvUINT32 vAssocMsrList, PUINT32 currentBmsFileIndex, const _AML_TYPE_ aml_type);

	//void FindUnusedStationsInIgnoredMeasurements(vdnaMsrPtr* vMeasurements, pvASLPtr vAssocStnList, pvUINT32 vAssocMsrList, pvstring vUnusedStns, pvUINT32 vIgnoredMsrs);

	void SerialiseXmlStn(std::ifstream* ifs_stns, std::ofstream* ofs_dynaml);
	void SerialiseXmlMsr(std::ifstream* ifs_stns, std::ifstream* ifs_msrs, std::ofstream* ofs_dynaml);

	void InitialiseDynaMLSepStationFile(const project_settings& p, vifm_t* vinput_file_meta,
		const std::string& stnfilename, std::ofstream* dynaml_stn_file);
	void InitialiseDynaMLSepMeasurementFile(const project_settings& p, vifm_t* vinput_file_meta,
		const std::string& msrfilename, std::ofstream* dynaml_msr_file);

	void InitialiseDynaMLFile(const project_settings& p, vifm_t* vinput_file_meta,
		const std::string& outfilename, std::ofstream* dynaml_file);

	void DetermineBoundingBox();
	void BuildExtractStationsList(const std::string& stnList, pvstring vstnList);
	
	void RemoveNonMeasurements(const UINT32& block, pvmsr_t binaryMsr);

	void SignalComplete();
	void SignalExceptionParseDNA(const std::string& message, const std::string& sBuf, const int& column_no);
	void SignalExceptionParse(std::string msg, int i);
	void SignalExceptionInterop(std::string msg, int i, const char *streamType, ...);

	
	project_settings		projectSettings_;
	binary_file_meta_t		bst_meta_;
	binary_file_meta_t		bms_meta_;
	CDnaDatum				datum_;

	dna_stn_fields			dsl_, dsw_;
	dna_msr_fields			dml_, dmw_;
	
	volatile double percentComplete_;

	std::ifstream*	ifsInputFILE_;
	size_t		sifsFileSize_;
	bool		isProcessing_;
	
	double		bbox_upperLat_;
	double		bbox_upperLon_;
	double		bbox_lowerLat_;
	double		bbox_lowerLon_;

	PUINT32		_pclusterID;

	MsrTally*	p_parsemsr_tally;
	StnTally*	p_parsestn_tally;

	UINT32		m_binaryRecordCount;
	UINT32		m_dbidRecordCount;
	UINT32		m_lineNo;
	UINT32		m_columnNo;

	std::string		m_strProjectDefaultEpsg;
	std::string		m_strProjectDefaultEpoch;
	std::string		m_msrComments;

	vvUINT32	v_ISL_;				// Inner stations
	vvUINT32	v_JSL_;				// Junction stations
	vvUINT32	v_CML_;				// Measurements

	it_pair_string_vUINT32	it_stnmap_range;
	
	v_string_uint32_pair vStnsMap_sortName_;		// Station Name Map sorted on name (string)
	
	_PARSE_STATUS_ parseStatus_;

	_INPUT_FILE_TYPE_ m_ift;
	_INPUT_DATA_TYPE_ m_idt;

	msr_database_id_map		m_msr_db_map;

	v_discontinuity_tuple	stn_discontinuities_;
	bool					m_discontsSortedbyName;

	v_string_string_pair	stn_renamed_;

	// Database management
	v_msr_database_id_map	v_msr_db_map_;
	bool					databaseIDsLoaded_;
	bool					databaseIDsSet_;

	void LoadDatabaseId();

	vUINT32 v_measurementCount_, v_unknownsCount_, v_ContiguousNetList_, v_parameterStationCount_;

	bool _filespecifiedreferenceframe;
	bool _filespecifiedepoch;

};

}	// namespace dynamlinterop
}	// namespace dynadjust

#endif /* DNAINTEROP_H_ */
