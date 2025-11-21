//============================================================================
// Name         : dnainterop.hpp
// Author       : Roger Fraser
// Contributors :
// Version      : 1.00
// Copyright    : Copyright 2017 Geoscience Australia
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

// MKL/BLAS/LAPACK is included via dnamatrix_contiguous.hpp
// #include <mkl.h>

#include <exception>
#include <system_error>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <string>
#include <utility>
#include <vector>
#include <map>
#include <cstdarg>
#include <math.h>

#include <boost/shared_ptr.hpp>
#include <boost/bind/bind.hpp>
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/thread/mutex.hpp>

#include <include/exception/dnaexception.hpp>
#include <include/config/dnaexports.hpp>
#include <include/config/dnaversion.hpp>
#include <include/config/dnaconsts.hpp>
#include <include/config/dnatypes.hpp>
#include <include/config/dnatypes-gui.hpp>

#include <include/math/dnamatrix_contiguous.hpp>

#include <include/parameters/dnaepsg.hpp>
#include <include/parameters/dnadatum.hpp>

#include <dynadjust/dnaimport/dnaparser_pimpl.hxx>

#include <include/io/dnaiodna.hpp>
#include <include/io/dnaiobst.hpp>
#include <include/io/dnaiobms.hpp>
#include <include/io/dnaioaml.hpp>
#include <include/io/dnaioasl.hpp>
#include <include/io/dnaiomap.hpp>
#include <include/io/dnaioseg.hpp>
#include <include/io/dnaiosnx.hpp>
#include <include/io/dnaioscalar.hpp>

#include <include/functions/dnatemplatefuncs.hpp>
#include <include/functions/dnatemplatestnmsrfuncs.hpp>
#include <include/functions/dnaiostreamfuncs.hpp>
#include <include/functions/dnastringfuncs.hpp>

using namespace std;
using namespace boost;
using namespace boost::filesystem;

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
	_PARSE_STATUS_ ParseInputFile(const string& filename, vdnaStnPtr* vStations, PUINT32 stnCount, 
		vdnaMsrPtr* vMeasurements, PUINT32 msrCount, 
		PUINT32 clusterID, input_file_meta_t* input_file_meta,
		string* success_msg, project_settings* p);

	_PARSE_STATUS_ LoadDNAGeoidFile(const string& fileName, vdnaStnPtr* vStations);

	void RemoveIgnoredMeasurements(vdnaMsrPtr* vMeasurements, MsrTally* parsemsrTally);
	void IncludeMeasurementTypes(const string& includeMsrs, vdnaMsrPtr* vMeasurements, MsrTally* parsemsrTally);
	void ExcludeMeasurementTypes(const string& excludeMsrs, vdnaMsrPtr* vMeasurements, MsrTally* parsemsrTally);
	void ExcludeAllOutsideBoundingBox(vdnaStnPtr* vStations, vdnaMsrPtr* vMeasurements, 
		StnTally* parsestnTally, MsrTally* parsemsrTally, pvstring pvUnusedStns, const project_settings& p,
		bool& splitXmsrs, bool& splitYmsrs);
	void ExtractStnsAndAssociatedMsrs(const string& stnListInclude, const string& stnListExclude, vdnaStnPtr* vStations, vdnaMsrPtr* vMeasurements, 
		StnTally* parsestnTally, MsrTally* parsemsrTally, pvstring vExcludedStns, const project_settings& p,
		bool& splitXmsrs, bool& splitYmsrs);
	void ImportStnsMsrsFromBlock(vdnaStnPtr* vStations, vdnaMsrPtr* vMeasurements, const project_settings& p);

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
	void TrackDiscontinuitySite(const string& site, const string& site_renamed);
	void ApplyDiscontinuitiesMeasurements(vdnaMsrPtr* vMeasurements);

	void ApplyDiscontinuitiesMeasurements_GX(vector<CDnaGpsBaseline>* vGpsBaselines);
	void ApplyDiscontinuitiesMeasurements_Y(vector<CDnaGpsPoint>* vGpsPoints);
	void ApplyDiscontinuitiesMeasurements_D(vector<CDnaDirection>* vDirections, const date& site_date);
	
	void AddDiscontinuityStations(vdnaStnPtr* vstationsTotal);

	void EditGNSSMsrScalars(vdnaMsrPtr* vMeasurements, project_settings* p);
	void ApplyGNSSMsrScalar(vector<CDnaGpsBaseline>::iterator& _it_bsl, vscl_t& bslScalars);

	void MapMeasurementStations(vdnaMsrPtr* vMeasurements, pvASLPtr vAssocStnList, PUINT32 lMapCount, pvstring vUnusedStns, pvUINT32 vIgnoredMsrs);
	UINT32 ComputeMeasurementCount(vdnaMsrPtr* vMeasurements, const vUINT32& vIgnoredMsrs);

	// file handling
	void SerialiseDNA(vdnaStnPtr* vStations, vdnaMsrPtr* vMeasurements, const string& stnfilename, const string& msrfilename, const project_settings& p, vifm_t* vinput_file_meta, bool flagUnused=false);
	//void SerialiseDynaMLfromBinary(const string& outfilename, const project_settings& p, vifm_t* vinput_file_meta, bool flagUnused=false);
	void SerialiseDynaMLfromMemory(vdnaStnPtr* vStations, vdnaMsrPtr* vMeasurements, const string& outfilename, const project_settings& p, vifm_t* vinput_file_meta, bool flagUnused=false);
	//void SerialiseDynaMLSepfromBinary(const string& stnfilename, const string& msrfilename, const project_settings& p, vifm_t* vinput_file_meta, bool flagUnused=false);
	void SerialiseDynaMLSepfromMemory(vdnaStnPtr* vStations, vdnaMsrPtr* vMeasurements, const string& stnfilename, const string& msrfilename, const project_settings& p, vifm_t* vinput_file_meta, bool flagUnused=false);
	//void SerialiseGeoidData(vdnaStnPtr* vStations, const string& geofilename);
	void SimulateMSR(vdnaStnPtr* vStations, vdnaMsrPtr* vMeasurements, const string& msrfilename, const project_settings& p);

	void SerialiseBst(const string& bst_filename, vdnaStnPtr* vStations, pvstring vUnusedStns, vifm_t& vinput_file_meta, bool flagUnused=false);
	void SerialiseBms(const string& bms_filename, vdnaMsrPtr* vMeasurements, vifm_t& vinput_file_meta);
	void SerialiseMap(const string& stnmap_file);
	void SerialiseAml(const string& filename, pvUINT32 vAssocMsrList);
	void SerialiseAsl(const string& filename, pvASLPtr vAssocStnList);
	void SerialiseMapTextFile(const string& stnmap_file);
	void SerialiseDiscontTextFile(const string& discont_file);
	void SerialiseAmlTextFile(const string& bms_filename, const string& aml_filename, pvUINT32 vAML, pvASLPtr vAssocStnList, vdnaStnPtr* vStations);
	void SerialiseAslTextFile(const string& filename, pvASLPtr vAssocStnList, vdnaStnPtr* vStations);

	void SerialiseDatabaseId(const string& dbid_filename, pvdnaMsrPtr vMeasurements);

	void CompleteAssociationLists(vdnaMsrPtr* vMeasurements, pvASLPtr vAssocStnList, pvUINT32 vAssocMsrList);

	double GetProgress();
	inline bool IsProcessing() const { return isProcessing_; }
	inline MsrTally& GetMsrTally() const { return *p_parsemsr_tally; }
	inline StnTally& GetStnTally() const { return *p_parsestn_tally; }

	inline _PARSE_STATUS_ GetStatus() const { return parseStatus_; }

	void inline ResetFileOrder() const { g_fileOrder = 0; }
	void InitialiseDatum(const string& reference_frame);
	void UpdateEpoch(const vifm_t* vinput_file_meta);
	
	void PrintMeasurementsToStations(string& m2s_file, MsrTally* parsemsrTally,
		string& bst_file, string& bms_file, string& aml_file, pvASLPtr vAssocStnList);

	// Discontonuity file
	void ParseDiscontinuities(const string& fileName);

private:
	
	// DynaML files
	void ParseXML(const string& fileName, vdnaStnPtr* vStations, PUINT32 stnCount, 
							   vdnaMsrPtr* vMeasurements, PUINT32 msrCount, PUINT32 clusterID, 
							   string& fileEpsg, string& fileEpoch, string* success_msg);
	
	// SINEX files
	void ParseSNX(const string& fileName, vdnaStnPtr* vStations, PUINT32 stnCount, 
							   vdnaMsrPtr* vMeasurements, PUINT32 msrCount, PUINT32 clusterID);
	
	// DNA Ascii files
	//void ParseDNAVersion(const INPUT_DATA_TYPE& idt);
	void ParseDNA(const string& fileName, vdnaStnPtr* vStations, PUINT32 stnCount, 
							   vdnaMsrPtr* vMeasurements, PUINT32 msrCount, PUINT32 clusterID, 
							   string& fileEpsg, string& fileEpoch);
	void ParseDNASTN(vdnaStnPtr* vStations, PUINT32 stnCount);
	void ParseDNAMSR(pvdnaMsrPtr vMeasurements, PUINT32 msrCount, PUINT32 clusterID);

	//void SetDefaultReferenceFrame(vdnaStnPtr* vStations, vdnaMsrPtr* vMeasurements);

	void ParseDNAMSRAngular(const string& sBuf, dnaMsrPtr& msr_ptr);
	void ParseDNAMSRCoordinate(const string& sBuf, dnaMsrPtr& msr_ptr);
	UINT32 ParseDNAMSRDirections(string& sBuf, dnaMsrPtr& msr_ptr, bool ignoreMsr);
	void ParseDNAMSRGPSBaselines(string& sBuf, dnaMsrPtr& msr_ptr, bool ignoreMsr);
	void ParseDNAMSRGPSPoints(string& sBuf, dnaMsrPtr& msr_ptr, bool ignoreMsr);
	void ParseDNAMSRLinear(const string& sBuf, dnaMsrPtr& msr_ptr);
	void ParseDNAMSRCovariance(CDnaCovariance& cov);

	string ParseAngularValue(const string& sBuf, const string& calling_function);
	string ParseLinearValue(const string& sBuf, const string& msrName, const string& calling_function);
	string ParseGPSMsrValue(const string& sBuf, const string& element, const string& calling_function);
	string ParseGPSVarValue(const string& sBuf, const string& element, const UINT32 location, const UINT32 width, const string& calling_function);
	string ParseInstrumentValue(const string& sBuf, const string& calling_function);
	string ParseTargetValue(const string& sBuf, const string& calling_function);
	string ParseTarget2Value(const string& sBuf, const string& calling_function);
	string ParseStdDevValue(const string& sBuf, const string& calling_function);
	string ParseInstHeightValue(const string& sBuf, const string& calling_function);
	string ParseTargHeightValue(const string& sBuf, const string& calling_function);
	string ParseMsrCountValue(const string& sBuf, UINT32& msrCount, const string& calling_function);
	string ParseScaleVValue(const string& sBuf, const string& calling_function);
	string ParseScalePValue(const string& sBuf, const string& calling_function);
	string ParseScaleLValue(const string& sBuf, const string& calling_function);
	string ParseScaleHValue(const string& sBuf, const string& calling_function);
	string ParseRefFrameValue(const string& sBuf, const string& calling_function);
	string ParseEpochValue(const string& sBuf, const string& calling_function);

	void ParseDatabaseIds(const string& sBuf, const string& calling_function, const char msrType);
	void ParseDatabaseClusterId(const string& sBuf, const string& calling_function);
	void ParseDatabaseMsrId(const string& sBuf, const string& calling_function);

	void LoadNetworkFiles(pvstn_t binaryStn, pvmsr_t binaryMsr, const project_settings& projectSettings, bool loadSegmentFile);
	void LoadBinaryFiles(pvstn_t binaryStn, pvmsr_t binaryMsr);
	void LoadSegmentationFile();

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
	void ExtractAssociatedStns_GX(vector<CDnaGpsBaseline>* vGpsBaselines, pvstring pvUsedStns);
	void ExtractAssociatedStns_Y(vector<CDnaGpsPoint>* vGpsPoints, pvstring pvUsedStns);
	void ExtractAssociatedStns_D(vector<CDnaDirection>* vDirections, pvstring pvUsedStns);

	void MapMeasurementStationsBsl(vector<CDnaGpsBaseline>* vGpsBaselines,
		pvASLPtr vAssocStnList, PUINT32 lMapCount);
	void MapMeasurementStationsPnt(vector<CDnaGpsPoint>* vGpsPoints,
		pvASLPtr vAssocStnList, PUINT32 lMapCount);
	void MapMeasurementStationsDir(vector<CDnaDirection>* vDirections,
		pvASLPtr vAssocStnList, PUINT32 lMapCount);
	
	//void RenameStationsMsr(CDnaMeasurement* msr, v_string_string_pair& variables);
	void RenameStationsBsl(vector<CDnaGpsBaseline>* vGpsBaselines, v_string_vstring_pair& stnRenaming);
	void RenameStationsPnt(vector<CDnaGpsPoint>* vGpsPoints, v_string_vstring_pair& stnRenaming);
	void RenameStationsDir(vector<CDnaDirection>* vDirections, v_string_vstring_pair& stnRenaming);
	
	void CompleteASLDirections(_it_vdnamsrptr _it_msr, vector<CDnaDirection>* vDirections, pvASLPtr vAssocStnList, 
		pvUINT32 vAssocMsrList, PUINT32 currentBmsFileIndex);
	void CompleteASLGpsBaselines(vector<CDnaGpsBaseline>* vGpsBaselines, pvASLPtr vAssocStnList, 
		pvUINT32 vAssocMsrList, PUINT32 currentBmsFileIndex);
	void CompleteASLGpsPoints(vector<CDnaGpsPoint>* vGpsPoints, pvASLPtr vAssocStnList, 
		pvUINT32 vAssocMsrList, PUINT32 currentBmsFileIndex);

	//void FindUnusedStationsInIgnoredMeasurements(vdnaMsrPtr* vMeasurements, pvASLPtr vAssocStnList, pvUINT32 vAssocMsrList, pvstring vUnusedStns, pvUINT32 vIgnoredMsrs);

	void SerialiseXmlStn(std::ifstream* ifs_stns, std::ofstream* ofs_dynaml);
	void SerialiseXmlMsr(std::ifstream* ifs_stns, std::ifstream* ifs_msrs, std::ofstream* ofs_dynaml);

	void InitialiseDynaMLSepStationFile(const project_settings& p, vifm_t* vinput_file_meta,
		const string& stnfilename, std::ofstream* dynaml_stn_file);
	void InitialiseDynaMLSepMeasurementFile(const project_settings& p, vifm_t* vinput_file_meta,
		const string& msrfilename, std::ofstream* dynaml_msr_file);

	void InitialiseDynaMLFile(const project_settings& p, vifm_t* vinput_file_meta,
		const string& outfilename, std::ofstream* dynaml_file);

	void DetermineBoundingBox();
	void BuildExtractStationsList(const string& stnList, pvstring vstnList);
	
	void RemoveNonMeasurements(const UINT32& block, pvmsr_t binaryMsr);

	void SignalComplete();
	void SignalExceptionParseDNA(const string& message, const string& sBuf, const int& column_no);
	void SignalExceptionParse(string msg, int i);
	void SignalExceptionInterop(string msg, int i, const char *streamType, ...);

	
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

	string		m_strProjectDefaultEpsg;
	string		m_strProjectDefaultEpoch;
	string		m_msrComments;

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

	//v_msr_database_id_map	v_msr_db_map;
	//UINT32					m_bmsIndex;
	//char					m_msrStart;
};

}	// namespace dynamlinterop
}	// namespace dynadjust

#endif /* DNAINTEROP_H_ */
