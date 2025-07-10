//============================================================================
// Name         : dnaioseg.hpp
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
// Description  : DynAdjust segmentation file io operations
//============================================================================

#ifndef DNAIOSEG_H_
#define DNAIOSEG_H_

#if defined(_MSC_VER)
	#if defined(LIST_INCLUDES_ON_BUILD) 
		#pragma message("  " __FILE__) 
	#endif
#endif

#include <include/io/dynadjust_file.hpp>
#include <include/config/dnatypes.hpp>
#include <include/functions/dnaiostreamfuncs.hpp>
#include <include/functions/dnatemplatestnmsrfuncs.hpp>
#include <include/measurement_types/dnameasurement.hpp>


using namespace dynadjust::measurements;

namespace dynadjust {
namespace iostreams {

class dna_io_seg : public DynadjustFile
{
public:
	dna_io_seg(void) {};
	dna_io_seg(const dna_io_seg& seg) : DynadjustFile(seg) {};
	virtual ~dna_io_seg(void) {};

	dna_io_seg& operator=(const dna_io_seg& rhs);

	void load_seg_file_header(const std::string& seg_filename, std::istream& seg_file, UINT32& blockCount, 
		UINT32& blockThreshold, UINT32& minInnerStns);

	void load_seg_file_header_f(const std::string& seg_filename, UINT32& blockCount, 
		UINT32& blockThreshold, UINT32& minInnerStns);

	void load_seg_file(const std::string& seg_filename, UINT32& blockCount, 
		UINT32& blockThreshold, UINT32& minInnerStns,
		vvUINT32& v_ISL, vvUINT32& v_JSL, vvUINT32& v_CML,
		bool loadMetrics,
		pvmsr_t bmsBinaryRecords, pvUINT32 v_measurementCount, 
		pvUINT32 v_unknownsCount, pvUINT32 v_ContiguousNetList,
		pvUINT32 v_parameterStationCount);

	void create_stn_appearance_list(vv_stn_appear& v_paramStnAppearance,
		const vvUINT32& paramStationList,
		vASL& assocStnList);

	void write_seg_block(std::ostream &os, 
		const vUINT32& vISL, const vUINT32& vJSL, const vUINT32& vCML, 
		const UINT32& currentBlock, 
		const vstn_t* bstBinaryRecords, const vmsr_t* bmsBinaryRecords, 
		bool PRINT_NAMES=false);

	void write_seg_file(const std::string& seg_filename, const std::string& bst_filename, const std::string& bms_filename,
		const UINT32& min_inner_stns, const UINT32& max_block_stns,
		const std::string& seg_starting_stns, const vstring& vinitialStns,
		const std::string& command_line_arguments,
		vvUINT32& v_ISL, vvUINT32& v_JSL, vvUINT32& v_CML,
		vUINT32& v_ContiguousNetList, const pvstn_t bstBinaryRecords, const pvmsr_t bmsBinaryRecords);

	void build_free_stn_availability(vASL& assocStnList, v_freestn_pair& freeStnList);

	void write_stn_appearance(const std::string& sap_filename, const v_stn_block_map& stnAppearance);
	void load_stn_appearance(const std::string& sap_filename, v_stn_block_map& stnAppearance);

protected:

};

}	// namespace measurements
}	// namespace dynadjust

#endif
