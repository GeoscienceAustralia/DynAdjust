//============================================================================
// Name         : dnatypes-structs.hpp
// Author       : Roger Fraser
// Contributors : Dale Roberts <dale.o.roberts@gmail.com>
//              : Dale Roberts
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
// Description  : Custom structures and complex types for DynAdjust
//============================================================================

#pragma once
#ifndef DNATYPES_STRUCTS_H_
#define DNATYPES_STRUCTS_H_

#include <cstring>		// memset
#include <stdio.h>		// snprintf
#include <string>
#include <vector>
#include <include/config/dnatypes-basic.hpp>
#include <include/config/dnatypes-containers.hpp>

/////////////////////////////////////////////////////////////
// Custom pair type to manage appearances of a station
template <class T1=UINT32, class T2=bool>
struct stn_appearance_t
{
	typedef T1 station_type;
	typedef T2 appearance_type;

	T1 station_id;
	T2 first_appearance_fwd;
	T2 first_appearance_rev;

	void set_id(const T1& id) { station_id = id; }
	void first_fwd() { first_appearance_fwd = true; }
	void first_rev() { first_appearance_rev = true; }

	stn_appearance_t()
		: station_id(0), first_appearance_fwd(false), first_appearance_rev(false) {}
	stn_appearance_t(const T1& id, const T2& f, const T2& r)
		: station_id(id), first_appearance_fwd(f), first_appearance_rev(r) {}
	template <class U, class V>
	stn_appearance_t (const stn_appearance_t<U, V> &p)
		: station_id(p.station_id)
		, first_appearance_fwd(p.first_appearance_fwd)
		, first_appearance_rev(p.first_appearance_rev) {}
};

typedef stn_appearance_t<UINT32, bool> stn_appear;
typedef std::vector<stn_appear> v_stn_appear;
typedef v_stn_appear::iterator it_vstn_appear;
typedef std::vector<v_stn_appear> vv_stn_appear;
typedef vv_stn_appear::iterator it_vvstn_appear;
/////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////
// Custom pair type to map the location of stations in blocks
template <class T=UINT32, class U=bool>
struct stn_block_map_t
{
	T block_no;
	U first_appearance_fwd;
	U first_appearance_rev;
	U valid_stn;

	stn_block_map_t()
		: block_no(0)
		, first_appearance_fwd(false), first_appearance_rev(false)
		, valid_stn(false) {}
	stn_block_map_t(const T& block, const U fwd, const U rev, const U validity)
		: block_no(block)
		, first_appearance_fwd(fwd), first_appearance_rev(rev)
		, valid_stn(validity) {}
	stn_block_map_t (const stn_block_map_t<T> &p)
		: block_no(p.block_no)
		, first_appearance_fwd(p.first_appearance_fwd)
		, first_appearance_rev(p.first_appearance_rev)
		, valid_stn(p.valid_stn) {}
	void firstAppearanceFwd(const T& block) {
		block_no = block;
		first_appearance_fwd = true;
	}
	void firstAppearanceRev() {
		first_appearance_rev = true;
	}
};

typedef stn_block_map_t<UINT32> stn_block_map;
typedef std::vector<stn_block_map> v_stn_block_map;
typedef v_stn_block_map::iterator it_vstn_block_map;
/////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////
// Custom pair type for "free stations" vector -default value is free
template <class T1=UINT32, class T2=bool>
struct freestnpair_t
{
	typedef T1 stnindex_type;
	typedef T2 free_type;

	T1 stn_index;
	T2 available;
	void consume() { available = false; }
	const T2 isfree() const { return available; }
	freestnpair_t() : stn_index(T1()), available(true) {}
	freestnpair_t(const T1& x, const T2& y) : stn_index(x), available(y) {}
	template <class U, class V>
	freestnpair_t (const freestnpair_t <U, V> &p) : stn_index(p.stn_index), available(p.available) { }
};

typedef freestnpair_t<UINT32, bool> freestn_pair;
typedef std::vector<freestn_pair> v_freestn_pair;
typedef v_freestn_pair::iterator it_freestn_pair;
/////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////
// Custom pair type for AML vector -default value is free
template <class T1=UINT32, class T2=bool>
struct amlpair_t
{
	typedef T1 bmsindex_type;
	typedef T2 free_type;

	T1 bmsr_index;
	T2 available;
	void consume() { available = false; }
	amlpair_t() : bmsr_index(T1()), available(true) {}
	amlpair_t(const T1& x, const T2& y) : bmsr_index(x), available(y) {}
	template <class U, class V>
	amlpair_t (const amlpair_t<U, V> &p) : bmsr_index(p.bmsr_index), available(p.available) { }
};

typedef amlpair_t<UINT32, bool> aml_pair;
typedef std::vector<aml_pair> v_aml_pair;
typedef v_aml_pair::iterator it_aml_pair;
/////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////
// Custom pair type for MT adjustments -default value is not adjusted (false)
template <class T1=UINT32, class T2=bool>
struct sequential_adj_t
{
	typedef T1 blockindex_type;
	typedef T2 adjusted_type;

	T1 block_index;
	T2 adjusted;
	void solution() { adjusted = true; }
	void nosolution() { adjusted = false; }
	sequential_adj_t()
		: block_index(T1()), adjusted(false) {}
	sequential_adj_t(const T1& x, const T2& y)
		: block_index(x), adjusted(y) {}
	template <class U, class V>
	sequential_adj_t (const sequential_adj_t<U, V> &p)
		: block_index(p.block_index), adjusted(p.adjusted) { }
};

typedef sequential_adj_t<UINT32, bool> sequential_adj;
typedef std::vector<sequential_adj> v_sequential_adj;
typedef v_sequential_adj::iterator it_sequential_adj;
/////////////////////////////////////////////////////////////

// Scalar struct
typedef struct scl_t {
    scl_t(): station1(""), station2(""), v_scale(1.), p_scale(1.), l_scale(1.), h_scale(1.) {}

    std::string station1;
    std::string station2;
    double v_scale; // phi, n or X scalar
    double p_scale; // lambda, e or Y scalar
    double l_scale; // height, up or Z scalar
    double h_scale; // matrix scalar
} scalar_t;

typedef std::vector<scalar_t> vscl_t, *pvscl_t;
typedef vscl_t::iterator it_vscl_t, *pit_vscl_t;
typedef vscl_t::const_iterator it_vscl_t_const;

// Statistics summary
typedef struct {
	double	_fwdChiSquared;
	double	_revChiSquared;
	double	_rigChiSquared;
	int	_degreesofFreedom;
} statSummary_t;

// Station corrections
typedef struct stationCorrections {
	stationCorrections(const std::string& stn="")
		: _station(stn), _azimuth(0.), _vAngle(0.), _sDistance(0.)
		, _hDistance(0.), _east(0.), _north(0.), _up(0.) {}

	std::string	_station;
	double	_azimuth;
	double	_vAngle;
	double	_sDistance;
	double	_hDistance;
	double	_east;
	double	_north;
	double	_up;
} stationCorrections_t;

typedef std::vector<stationCorrections_t> vstnCor_t, *pvstnCor_t;
typedef std::vector<vstnCor_t> vvstnCor_t;
typedef vstnCor_t::iterator it_vstnCor_t;

// Station position uncertainty
typedef struct {
	std::string	_station;
	double	_latitude;
	double	_longitude;
	double	_hzPosU;
	double	_vtPosU;
	double	_semimMajor;
	double	_semimMinor;
	double	_orientation;
	double	_xx;
	double	_xy;
	double	_xz;
	double	_yy;
	double	_yz;
	double	_zz;
} stationPosUncertainty_t;

typedef std::vector<stationPosUncertainty_t> vstnPU_t, *pvstnPU_t;
typedef std::vector<vstnPU_t> vvstnPU_t;
typedef vstnPU_t::iterator it_vstnPU_t;

// Segment parameters
typedef struct {
	std::string	_networkName;
	vstring	_initialStns;
	UINT32	_minInnerStns;
	UINT32	_maxTotalStns;
	UINT16	_sortStnsByMsrs;
	bool	_quiet;
	UINT16	_verbose;			// Give detailed information about what dnainterop is doing.
								// 0: No information (default)
								// 1: Helpful information
								// 2: Extended information
								// 3: Debug level information
} segmentParam_t;

// Block metadata
typedef struct block_meta {
	block_meta() : _blockIsolated(false), _blockFirst(false), _blockLast(false), _blockIntermediate(false) {}

	bool	_blockIsolated;		// Does this block comprise a single contiguous network (an isolated block)
								//		True:	A single block in isolation
								//		False:	One of two or more segmented blocks
	bool	_blockFirst;		// Is this the first block in a list, or an only block?
	bool	_blockLast;			// Is this the last block in a list, or an only block?
	bool	_blockIntermediate;	// Not the first or last, but an intermediate block of a list of three or more
} blockMeta_t;

// Station structure for binary station file
typedef struct stn_t {
	stn_t(const short& u=0)
		: suppliedStationType(LLH_type_i), initialLatitude(0.), currentLatitude(0.), initialLongitude(0.), currentLongitude(0.)
		, initialHeight(0.), currentHeight(0.), suppliedHeightRefFrame(ELLIPSOIDAL_type_i)
		, geoidSep(0.), geoidSepUnc(0.), meridianDef(0.), verticalDef(0.), zone(u)
		, fileOrder(0), nameOrder(0), clusterID(0), unusedStation(FALSE)
	{
		memset(stationName, '\0', sizeof(stationName));
		memset(stationNameOrig, '\0', sizeof(stationNameOrig));
		memset(stationConst, '\0', sizeof(stationConst));
		memset(stationType, '\0', sizeof(stationType));
		memset(description, '\0', sizeof(description));
		// GDA2020, lat, long, height
		memset(epsgCode, '\0', sizeof(epsgCode));
		snprintf(epsgCode, sizeof(epsgCode), "7843");
		memset(epoch, '\0', sizeof(epoch));
		memset(plate, '\0', sizeof(plate));
	}

	char	stationName[STN_NAME_WIDTH];			// 30 characters
	char    stationNameOrig[STN_NAME_ORIG_WIDTH];	// 39 characters: 30 (name) + 1 (_) + 8 (date)
	char	stationConst[STN_CONST_WIDTH];			// constraint: lat, long, height
	char	stationType[STN_TYPE_WIDTH];			// type: LLH, UTM, XYZ
	UINT16	suppliedStationType;					// type supplied by the user (required for adjust)
	double	initialLatitude;						// initial estimate
	double	currentLatitude;						// current estimate
	double	initialLongitude;
	double	currentLongitude;
	double	initialHeight;					// initialHeight and currentHeight are always assumed to be
	double	currentHeight;					// ellipsoidal (ELLIPSOIDAL_type_i).  If the user runs geoid using
											// the convert option, then suppliedHeightRefFrame is set to
											// ORTHOMETRIC_type_i and currentHeight is set to (initialHeight + N).
	UINT16	suppliedHeightRefFrame;			// Used to signify which reference frame supplied height refers to
	float	geoidSep; 					 	// ellipsoid / geoid separation
	float   geoidSepUnc;                    // ellipsoid / geoid separation uncertainty (std deviation in m)
	double	meridianDef;					// deflection in meridian (N/S)
	double	verticalDef;					// deflection in vertical (E/W)
	short	zone;
	char	description[STN_DESC_WIDTH];	// 128 characters
	UINT32	fileOrder;						// original file order
	UINT32	nameOrder;						// station name sorted position
	UINT32	clusterID;						// cluster ID (which cluster this station belongs to)
	UINT16	unusedStation;					// is this station unused?
	char	epsgCode[STN_EPSG_WIDTH];		// epsg ID, i.e. NNNNN (where NNNNN is in the range 0-32767)
	char	epoch[STN_EPOCH_WIDTH];			// date, i.e. "DD.MM.YYYY" (10 chars)
											// if datum is dynamic, Epoch is YYYY MM DD
											// if datum is static, Epoch is ignored
	char	plate[STN_PLATE_WIDTH];			// Tectonic plate identifier. Typically two characters.
} station_t;

typedef std::vector<station_t> vstn_t, *pvstn_t;
typedef vstn_t::iterator it_vstn_t;
typedef vstn_t::const_iterator it_vstn_t_const;
typedef std::vector<statSummary_t> vsummary_t, *pvsummary_t;

typedef std::pair <station_t, std::string> stn_t_string_pair;
typedef std::vector<stn_t_string_pair> v_stn_string;
typedef v_stn_string::iterator it_stn_string;

// Input file metadata
typedef struct input_file_meta {
	char	filename[FILE_NAME_WIDTH+1];	// Input file path
	char	epsgCode[STN_EPSG_WIDTH+1];		// Input file epsg ID, i.e. NNNNN (where NNNNN is in the range 0-32767). "Mixed" if stations are on different reference frames
	char	epoch[STN_EPOCH_WIDTH+1];		// Input file epoch
	UINT16	filetype;						// Input file type (geodesyml, dynaml, dna, csv, sinex)
	UINT16	datatype;						// Input data type (station, measurement, both)
} input_file_meta_t;

// Source file metadata for measurement provenance
typedef struct source_file_meta {
	char	filename[FILE_NAME_WIDTH+1];
} source_file_meta_t;

// Binary file metadata
typedef struct binary_file_meta {
	binary_file_meta ()
		: binCount(0), reduced(false), reftran(false), geoid(false)
		, inputFileCount(0), inputFileMeta(NULL)
		, sourceFileCount(0), sourceFileMeta(nullptr) {}
	binary_file_meta (const std::string& app_name)
		: binCount(0), reduced(false), reftran(false), geoid(false)
		, inputFileCount(0), inputFileMeta(NULL)
		, sourceFileCount(0), sourceFileMeta(nullptr) {
            snprintf(modifiedBy, sizeof(modifiedBy), "%s", app_name.c_str());
	}
	~binary_file_meta() {
		if (inputFileMeta != NULL)
			delete []inputFileMeta;
		if (sourceFileMeta != nullptr)
			delete []sourceFileMeta;
	}

	binary_file_meta(const binary_file_meta&) = delete;
	binary_file_meta& operator=(const binary_file_meta&) = delete;

	binary_file_meta(binary_file_meta&& rhs) noexcept
		: binCount(rhs.binCount), reduced(rhs.reduced)
		, reftran(rhs.reftran), geoid(rhs.geoid)
		, inputFileCount(rhs.inputFileCount), inputFileMeta(rhs.inputFileMeta)
		, sourceFileCount(rhs.sourceFileCount), sourceFileMeta(rhs.sourceFileMeta)
	{
		memcpy(modifiedBy, rhs.modifiedBy, sizeof(modifiedBy));
		memcpy(epsgCode, rhs.epsgCode, sizeof(epsgCode));
		memcpy(epoch, rhs.epoch, sizeof(epoch));
		rhs.inputFileMeta = nullptr;
		rhs.sourceFileMeta = nullptr;
	}

	binary_file_meta& operator=(binary_file_meta&& rhs) noexcept
	{
		if (this != &rhs)
		{
			delete []inputFileMeta;
			delete []sourceFileMeta;
			binCount = rhs.binCount;
			reduced = rhs.reduced;
			reftran = rhs.reftran;
			geoid = rhs.geoid;
			memcpy(modifiedBy, rhs.modifiedBy, sizeof(modifiedBy));
			memcpy(epsgCode, rhs.epsgCode, sizeof(epsgCode));
			memcpy(epoch, rhs.epoch, sizeof(epoch));
			inputFileCount = rhs.inputFileCount;
			inputFileMeta = rhs.inputFileMeta;
			sourceFileCount = rhs.sourceFileCount;
			sourceFileMeta = rhs.sourceFileMeta;
			rhs.inputFileMeta = nullptr;
			rhs.sourceFileMeta = nullptr;
		}
		return *this;
	}
    std::uint64_t		binCount;						// number of records in the binary file
	bool				reduced;						// indicates whether the data is reduced(true) or raw(false)
	char				modifiedBy[MOD_NAME_WIDTH+1];	// the program that modified this file
	char				epsgCode[STN_EPSG_WIDTH+1];		// epsg ID, i.e. NNNNN (where NNNNN is in the range 0-32767). "Mixed" if stations are on different reference frames
	char				epoch[STN_EPOCH_WIDTH+1];		// date, i.e. "DD.MM.YYYY" (10 chars)
	bool				reftran;						// the data has been transformed to another frame and/or epoch
	bool				geoid;							// geoid separation values have been obtained
    std::uint64_t		inputFileCount;					// Number of source file metadata elements
	input_file_meta_t*	inputFileMeta;					// Source file metadata
	std::uint64_t		sourceFileCount;				// Number of unique source files for measurements
	source_file_meta_t*	sourceFileMeta;					// Source file metadata for measurement provenance
} binary_file_meta_t;

typedef std::vector<input_file_meta_t> vifm_t;
typedef vifm_t::iterator it_vifm_t;
typedef std::vector<source_file_meta_t> vsfm_t;
typedef std::vector<binary_file_meta_t> vbfm_t;
typedef vbfm_t::iterator it_vbfm_t;

// Template functions
template <typename S>
S formatStnMsrFileSourceString(const vifm_t* vfile_meta, const size_t& file_type)
{
	std::string source_files("");
	bool this_file;

	for (std::size_t i(0); i<vfile_meta->size(); ++i)
	{
		this_file= false;
		switch (file_type)
		{
		case stn_data:
			if (vfile_meta->at(i).datatype == stn_data ||
				vfile_meta->at(i).datatype == stn_msr_data)
				this_file= true;
			break;
		case msr_data:
			if (vfile_meta->at(i).datatype == msr_data ||
				vfile_meta->at(i).datatype == stn_msr_data)
				this_file= true;
			break;
		case stn_msr_data:
			if (vfile_meta->at(i).datatype == stn_data ||
				vfile_meta->at(i).datatype == msr_data ||
				vfile_meta->at(i).datatype == stn_msr_data)
				this_file= true;
			break;
		}

		if (this_file)
		{
			source_files += vfile_meta->at(i).filename;
			source_files += " ";
		}
	}

	return source_files;
}

template <typename S = std::string>
S FormatFileType(const size_t& file_type)
{
	switch (file_type)
	{
	case dna:
		return "DNA";
		break;
	case sinex:
		return "SNX";
		break;
	case geodesyml:
		return "GML";
		break;
	case dynaml:
		return "XML";
		break;
	case csv:
		return "CSV";
		break;
	}
	return "unknown";
}

#endif // DNATYPES_STRUCTS_H_