//============================================================================
// Name         : dnasegment.hpp
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
// Description  : DynAdjust Network Segmentation library
//============================================================================

#ifndef DNASEGMENT_H_
#define DNASEGMENT_H_

#if defined(_MSC_VER)
	#if defined(LIST_INCLUDES_ON_BUILD) 
		#pragma message("  " __FILE__) 
	#endif
#endif

#include <exception>
#include <stdexcept>
#include <iostream>
#include <string>
#include <utility>
#include <vector>
#include <map>
#include <cstdarg>
#include <math.h>

#include <boost/shared_ptr.hpp>
#include <boost/chrono/time_point.hpp>
#include <boost/timer/timer.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/filesystem.hpp>

using namespace std;
using namespace boost;
using namespace boost::posix_time;
using namespace boost::filesystem;

#include <include/config/dnaexports.hpp>
#include <include/config/dnaversion.hpp>
#include <include/exception/dnaexception.hpp>
#include <include/functions/dnatemplatefuncs.hpp>
#include <include/functions/dnatemplatestnmsrfuncs.hpp>
#include <include/functions/dnatemplatedatetimefuncs.hpp>
#include <include/functions/dnaiostreamfuncs.hpp>
#include <include/functions/dnastringfuncs.hpp>
#include <include/config/dnatypes.hpp>
#include <include/config/dnaconsts.hpp>

#include <include/io/dnaioasl.hpp>
#include <include/io/dnaioaml.hpp>
#include <include/io/dnaiobst.hpp>
#include <include/io/dnaiobms.hpp>
#include <include/io/dnaiomap.hpp>
#include <include/io/dnaioseg.hpp>

using namespace dynadjust::measurements;
using namespace dynadjust::exception;
using namespace dynadjust::iostreams;

namespace dynadjust {
namespace networksegment {

// This class is exported from the dnaSegment.dll
#ifdef _MSC_VER
class DNASEGMENT_API dna_segment {
#else
class dna_segment {
#endif
public:
	dna_segment();
	dna_segment(const dna_segment& orig);
	virtual ~dna_segment();

	void PrepareSegmentation(project_settings* projectSettings);

	// Segment a memory-based network and store blocked lists in memory
	_SEGMENT_STATUS_ SegmentNetwork(project_settings* projectSettings, string* success_msg);

	void InitialiseSegmentation();
	
	double GetProgress() const;
	inline UINT32 currentBlock() const { return currentBlock_; }
	inline size_t blockCount() const { return vISL_.size(); }
	inline double averageblockSize() const { return averageBlockSize_; }
	inline UINT32 stationSolutionCount() const { return stationSolutionCount_; }
	inline UINT32 maxBlockSize() const { return maxBlockSize_; }
	inline UINT32 minBlockSize() const { return minBlockSize_; }
	
	static void coutVersion();
	void coutSummary() const;
	void coutCurrentBlockSummary(ostream &os);


	void LoadNetFile();
	void LoadBinaryFiles(const string& bstrfileName, const string& bmsrfileName);
	void LoadAssociationFiles(const string& aslfileName, const string& amlfileName);
	void LoadStationMap(const string& stnmap_file);

	void WriteSegmentedNetwork(const string& segfileName);
	void WriteFreeStnListSortedbyASLMsrCount();

	void VerifyStationConnections();

	string DefaultStartingStation();
	vstring StartingStations();
	
	inline _SEGMENT_STATUS_ GetStatus() const { return segmentStatus_; }

private:
	void SignalExceptionSerialise(const string& msg, const int& i, const char *streamType, ...);
	void ParseStartingStations();
	void BuildFirstBlock();
	void BuildNextBlock();
	void FinaliseBlock();

	void SortbyMeasurementCount(pvUINT32 vStnList);

	UINT32 SelectInner();
	void SelectJunction();
	
	it_vUINT32 MoveStation(vUINT32& fromList, it_vUINT32 it_from, vUINT32& toList, const UINT32& stn_index);
	void MoveFreeStn(it_vUINT32 it_freeisl, vUINT32& toList, const UINT32& stn_index, const string& type);
	void MoveFreeStnToInnerList(it_vUINT32 it_freeisl, const UINT32& stn_index);
	void MoveFreeStnToJunctionList(it_vUINT32 it_freeisl, const UINT32& stn_index);

	bool IncrementNextAvailableAMLIndex(UINT32& amlIndex, const UINT32& lastamlIndex);
	bool IncrementNextAvailableAMLIndex(it_aml_pair& _it_aml, const it_aml_pair& _it_lastaml);
	void GetInnerMeasurements(const UINT32& innerStation);
	void IdentifyLowestStationAssociation(pvUINT32 vStnList, vUINT32& totalStations,const int currentLevel, const int maxLevel, pvUINT32 vStnCount);
	void IdentifyInnerMsrsandAssociatedStns(const UINT32& innerStation, vUINT32& totalStations);
	void AddtoJunctionStnList(const vUINT32& msrStations);
	void AddtoJunctionStnList_sortedbyMsrCount(measurement_t* msr, const UINT32& innerStation);
	void AddtoJunctionStnList_sortedbyName(measurement_t* msr, const UINT32& innerStation);
	bool AddtoCurrentMsrList(const UINT32& amlindex, const vUINT32& msrStations);
	UINT32 GetAvailableMsrCount(const UINT32& stn_index);
	void FindCommonMeasurements();
	void MoveJuncttoInnerStnList();
	void MoveJunctiontoISL();
	void SetAvailableMsrCount();
	void RemoveDuplicateStations(pvstring vStations);
	void RemoveInvalidFreeStations();
	void RemoveNonMeasurements();
	void RemoveIgnoredMeasurements();
	void BuildFreeStationAvailabilityList();
	void BuildFreeStationAvailabilityList(vASL& assocStnList, v_freestn_pair& freeStnList);
	void BuildStationAppearanceList();
	void CalculateAverageBlockSize();


	void VerifyStationConnections_Block(const UINT32& block);
	void VerifyStationsandBuildBlock(bool validationOnly=false);

	project_settings	projectSettings_;

	_SEGMENT_STATUS_	segmentStatus_;
	bool				isProcessing_;
	UINT32				currentBlock_;
	UINT32				currentNetwork_;			// current contiguous network ID

	string				network_name_;				// network name
	UINT16				debug_level_;

	string				output_folder_;
	
	// functors
	CompareClusterID<measurement_t, UINT32> clusteridCompareFunc_;

	// block vectors
	vstring					vinitialStns_;
	vUINT32					v_ContiguousNetList_;	// vector of contiguous network IDs (corresponding to each block)
	vUINT32					vCurrJunctStnList_;
	vUINT32					vCurrInnerStnList_;
	vUINT32					vCurrMeasurementList_;
	
	vstn_t					bstBinaryRecords_;
	vmsr_t					bmsBinaryRecords_;
	vASL					vAssocStnList_;
	vUINT32					vASLCount_;
	vUINT32					vAssocMsrList_;
	v_aml_pair				vAssocFreeMsrList_;
	vUINT32					vfreeStnList_;
	vUINT32					vfreeMsrList_;			// vAssocMsrList_, less non-measurements and duplicate measurement references, sorted by ClusterID.
	v_string_uint32_pair	stnsMap_;

	v_freestn_pair			vfreeStnAvailability_;
	
	vvUINT32				vJSL_;
	vvUINT32				vISL_;
	vvUINT32				vCML_;

	std::ofstream				debug_file;
	std::ofstream				trace_file;

	double					averageBlockSize_;
	UINT32					stationSolutionCount_;
	UINT32					minBlockSize_;
	UINT32					maxBlockSize_;
};

}	// namespace networksegment
}	// namespace dynadjust


#endif /* DNASEGMENT_H_ */
