//============================================================================
// Name         : dnaiodna.hpp
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
// Description  : DNA file io helper
//============================================================================

#ifndef DNAIODNA_H_
#define DNAIODNA_H_

#if defined(_MSC_VER)
	#if defined(LIST_INCLUDES_ON_BUILD) 
		#pragma message("  " __FILE__) 
	#endif
#endif

#include <include/config/dnatypes-fwd.hpp>
#include <include/config/dnaconsts-iostream.hpp>
#include <include/io/dynadjust_file.hpp>
#include <include/io/dnaiodnatypes.hpp>
#include <include/functions/dnaintegermanipfuncs.hpp>
#include <include/functions/dnatemplatestnmsrfuncs.hpp>
#include <include/functions/dnaiostreamfuncs.hpp>
#include <include/functions/dnastringfuncs.hpp>
#include <include/functions/dnatemplatefuncs.hpp>
#include <include/measurement_types/dnastation.hpp>
#include <include/measurement_types/dnameasurement.hpp>
#include <include/parameters/dnadatum.hpp>

namespace dynadjust {
namespace iostreams {

class dna_io_dna : public DynadjustFile
{
public:
	dna_io_dna(void)
		: pv_msr_db_map_(0)
		, m_databaseIDsSet_(false)
		, m_filespecifiedReferenceFrame_(false)
		, m_filespecifiedEpoch_(false) {
	};
	dna_io_dna(const dna_io_dna& dna) 
		: DynadjustFile(dna) 
		, pv_msr_db_map_(0)
		, m_databaseIDsSet_(false)
		, m_filespecifiedReferenceFrame_(false)
		, m_filespecifiedEpoch_(false) {
	};
	virtual ~dna_io_dna(void) {};

	dna_io_dna& operator=(const dna_io_dna& rhs);

	void read_ren_file(const std::string& filename, pv_string_vstring_pair stnRenaming);

	void write_dna_files(vdnaStnPtr* vStations, vdnaMsrPtr* vMeasurements, 
		const std::string& stnfilename, const std::string& msrfilename, const std::string& networkname,
		const CDnaDatum& datum, const CDnaProjection& projection, bool flagUnused,
		const std::string& stn_comment, const std::string& msr_comment);
	
	void write_dna_files(pvstn_t vbinary_stn, pvmsr_t vbinary_msr, 
		const std::string& stnfilename, const std::string& msrfilename, const std::string& networkname,
		const CDnaDatum& datum, const CDnaProjection& projection, bool flagUnused,
		const std::string& stn_comment, const std::string& msr_comment);

	// CDnaStation
	void write_stn_file(vdnaStnPtr* vStations, const std::string& stnfilename, const std::string& networkname,
		const CDnaDatum& datum, const CDnaProjection& projection, bool flagUnused,
		const std::string& comment);
	// station_t
	void write_stn_file(pvstn_t vbinary_stn, const std::string& stnfilename, const std::string& networkname,
		const CDnaDatum& datum, const CDnaProjection& projection, bool flagUnused,
		const std::string& comment);
	
	// CDnaMeasurement
	void write_msr_file(vdnaMsrPtr* vMeasurements, const std::string& msrfilename, const std::string& networkname,
		const CDnaDatum& datum, const std::string& comment);

	// measurement_t
	void write_msr_file(const vstn_t& binaryStn, pvmsr_t vbinary_msr, const std::string& msrfilename, const std::string& networkname,
		const CDnaDatum& datum, const std::string& comment);

	void read_dna_header(std::ifstream* ptr, std::string& version, INPUT_DATA_TYPE& idt,
		CDnaDatum& referenceframe, 
		std::string& fileEpsg, std::string& fileEpoch, std::string& geoidversion, UINT32& count);

	inline const dna_stn_fields	dna_stn_positions() { return dsl_; }
	inline const dna_stn_fields	dna_stn_widths() { return dsw_; }
	inline const dna_msr_fields	dna_msr_positions() { return dml_; }
	inline const dna_msr_fields	dna_msr_widths() { return dmw_; }

	void set_dbid_ptr(pv_msr_database_id_map pv_msr_db_map);
	void set_source_file_meta_ptr(const source_file_meta_t* sfm, std::uint64_t count) {
		sourceFileMeta_ = sfm;
		sourceFileCount_ = count;
	}

	inline bool filespecifiedReferenceFrame() { return m_filespecifiedReferenceFrame_; }
	inline bool filespecifiedEpoch() { return m_filespecifiedEpoch_; }

protected:

private:

	void create_file_stn(std::ofstream* ptr, const std::string& filename);
	void create_file_msr(std::ofstream* ptr, const std::string& filename);
	void create_file_pointer(std::ofstream* ptr, const std::string& filename);
	void open_file_pointer(std::ifstream* ptr, const std::string& filename);

	void write_stn_header_data(std::ofstream* ptr, const std::string& networkname, const std::string& datum,
		const std::string& epoch, const size_t& count, const std::string& comment);
	void write_stn_header(std::ofstream* ptr, vdnaStnPtr* vStations, const std::string& networkname,
		const CDnaDatum& datum, bool flagUnused, const std::string& comment);
	void write_stn_header(std::ofstream* ptr, pvstn_t vbinary_stn, const std::string& networkname,
		const CDnaDatum& datum, bool flagUnused, const std::string& comment);

	void write_msr_header_data(std::ofstream* ptr, const std::string& networkname, const std::string& datum,
		const std::string& epoch, const size_t& count, const std::string& comment);
	void write_msr_header(std::ofstream* ptr, vdnaMsrPtr* vMeasurements, const std::string& networkname,
		const CDnaDatum& datum, const std::string& comment);
	void write_msr_header(std::ofstream* ptr, pvmsr_t vbinary_msrn, const std::string& networkname,
		const CDnaDatum& datum, const std::string& comment);

	void read_ren_data(std::ifstream* ptr, pv_string_vstring_pair stnRenaming);
	
	void prepare_sort_list(const UINT32 count);

	dna_stn_fields			dsl_, dsw_;
	dna_msr_fields			dml_, dmw_;

	vUINT32					vStationList_;

	pv_msr_database_id_map	pv_msr_db_map_;
	bool					m_databaseIDsSet_;
	bool                    m_filespecifiedReferenceFrame_;
	bool					m_filespecifiedEpoch_;

	const source_file_meta_t*	sourceFileMeta_ = nullptr;
	std::uint64_t				sourceFileCount_ = 0;
};

}	// namespace iostreams
}	// namespace dynadjust
	
#endif // DNAIODNA_H_
