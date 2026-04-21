//============================================================================
// Name         : dnaiodnatypes.hpp
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
// Description  : DNA file consts
//============================================================================

#ifndef DNAIODNATYPES_H_
#define DNAIODNATYPES_H_

#if defined(_MSC_VER)
	#if defined(LIST_INCLUDES_ON_BUILD) 
		#pragma message("  " __FILE__) 
	#endif
#endif

/// \cond
#include <sstream>
/// \endcond

#include <include/exception/dnaexception.hpp>
#include <include/functions/dnastrutils.hpp>

using namespace dynadjust::exception;

// msr_database_t is a struct intended to be contained within a vector
// whose size is identical to the binary measurement vector. Having a
// 1:1 relationship allows for direct referencing of database ID values
// without having to store the index to the measurement in either binary 
// measurement vector or the database ID vector.
typedef struct msr_database_t {
	UINT32	msr_id;
	UINT32	cluster_id;
	bool	is_msr_id_set;
	bool	is_cls_id_set;

	msr_database_t()
		: msr_id(0), cluster_id(0)
		, is_msr_id_set(false), is_cls_id_set(false) {}

} msr_database_id_map;

typedef std::vector<msr_database_id_map> v_msr_database_id_map, *pv_msr_database_id_map;
typedef v_msr_database_id_map::iterator it_vdbid_t;

typedef struct {
	UINT32 stn_name;
	UINT32 stn_const;
	UINT32 stn_type;
	UINT32 stn_e_phi_x;
	UINT32 stn_n_lam_y;
	UINT32 stn_ht_z;
	UINT32 stn_hemi_zo;
	UINT32 stn_desc;
} dna_stn_fields;

typedef struct {
	UINT32 msr_type;
	UINT32 msr_ignore;
	UINT32 msr_inst;
	UINT32 msr_targ1;
	UINT32 msr_targ2;
	UINT32 msr_linear;
	UINT32 msr_gps;
	UINT32 msr_gps_vcv_1;
	UINT32 msr_gps_vcv_2;
	UINT32 msr_gps_vcv_3;
	UINT32 msr_gps_vscale;
	UINT32 msr_gps_pscale;
	UINT32 msr_gps_lscale;
	UINT32 msr_gps_hscale;
	UINT32 msr_gps_reframe;
	UINT32 msr_gps_epoch;
	UINT32 msr_ang_d;
	UINT32 msr_ang_m;
	UINT32 msr_ang_s;
	UINT32 msr_stddev;
	UINT32 msr_inst_ht;
	UINT32 msr_targ_ht;
	UINT32 msr_id_msr;
	UINT32 msr_id_cluster;
	UINT32 msr_gps_obs_epoch;	// DNA v3.02+ column; 0 if file predates v3.02
} dna_msr_fields;

// DNA version 1.00
// Stations
template <class U>
struct _dna_stn_fields_100_
{
	static const U _locations_[8];
	static const U _widths_[8];
};

template <class U>
const U _dna_stn_fields_100_<U>::_locations_[8] = 
{
	0,		// station name
	10,		// constraints
	14,		// type
	17,		// e / phi / x
	31,		// n / lambda / y
	45,		// height / z
	59,		// hemisphere-zone
	63		// description
};

template <class U>
const U _dna_stn_fields_100_<U>::_widths_[8] = 
{
	10,		// station name
	3,		// constraints
	3,		// type
	14,		// e / phi / x
	14,		// n / lambda / y
	14,		// height / z
	3,		// hemisphere-zone
	128		// description
};

class dna_stn_fields_100
	: public _dna_stn_fields_100_<UINT16>
{
public:
};

// Measurements
template <class U>
struct _dna_msr_fields_100_
{
	static const U _locations_[24];
	static const U _widths_[24];
};

template <class U>
const U _dna_msr_fields_100_<U>::_locations_[24] = 
{
	0,		// msr type
	1,		// ignore 
	2,		// instrument name
	12,		// target 1
	22,		// target 2
	32,		// linear measurement
	32,		// gps measurement
	46,		// gps vcv 1
	57,		// gps vcv 2
	68,		// gps vcv 3
	32,		// gps vscale
	46,		// gps pscale
	57,		// gps lscale
	68,		// gps hscale
	79,		// gps ref frame
	89,		// gps epoch
	46,		// angular d
	50,		// angular m
	54,		// angular s
	60,		// standard deviation
	69,		// inst height
	76,		// target height
	83,		// msr id
	93		// cluster id
};

template <class U>
const U _dna_msr_fields_100_<U>::_widths_[24] = 
{
	1,		// msr type
	1,		// ignore 
	10,		// instrument name
	10,		// target 1
	10,		// target 2
	14,		// linear measurement
	14,		// gps measurement
	11,		// gps vcv 1
	11,		// gps vcv 2
	11,		// gps vcv 3
	14,		// gps vscale
	11,		// gps pscale
	11,		// gps lscale
	11,		// gps hscale
	10,		// gps ref frame
	10,		// gps epoch
	4,		// angular d
	4,		// angular m
	6,		// angular s
	9,		// standard deviation
	7,		// inst height
	7,		// target height
	10,		// msr id
	10		// cluster id
};

class dna_msr_fields_100
	: public _dna_msr_fields_100_<UINT16>
{
public:
};


// DNA version 3.00
// Stations
template <class U>
struct _dna_stn_fields_300_
{
	static const U _locations_[8];
	static const U _widths_[8];
};

template <class U>
const U _dna_stn_fields_300_<U>::_locations_[8] = 
{
	0,		// station name
	20,		// constraints
	24,		// type
	27,		// e / phi / x
	41,		// n / lambda / y
	55,		// height / z
	69,		// hemisphere-zone
	73		// description
};

template <class U>
const U _dna_stn_fields_300_<U>::_widths_[8] = 
{
	20,		// station name
	3,		// constraints
	3,		// type
	14,		// e / phi / x
	14,		// n / lambda / y
	14,		// height / z
	3,		// hemisphere-zone
	128		// description
};

class dna_stn_fields_300
	: public _dna_stn_fields_300_<UINT16>
{
public:
};


// Measurements
template <class U>
struct _dna_msr_fields_300_
{
	static const U _locations_[24];
	static const U _widths_[24];
};

template <class U>
const U _dna_msr_fields_300_<U>::_locations_[24] = 
{
	0,		// msr type
	1,		// ignore 
	2,		// instrument name
	22,		// target 1
	42,		// target 2
	62,		// linear measurement
	62,		// gps measurement
	76,		// gps vcv 1
	90,		// gps vcv 2
	104,	// gps vcv 3
	62,		// gps vscale
	69,		// gps pscale
	76,		// gps lscale
	83,		// gps hscale
	90,		// gps ref frame
	104,	// gps epoch
	76,		// angular d
	80,		// angular m
	82,		// angular s
	90,		// standard deviation
	99,		// inst height
	106,	// target height
	118,	// msr id
	128		// cluster id
};

template <class U>
const U _dna_msr_fields_300_<U>::_widths_[24] = 
{
	1,		// msr type
	1,		// ignore 
	20,		// instrument name
	20,		// target 1
	20,		// target 2
	14,		// linear measurement
	14,		// gps measurement
	14,		// gps vcv 1
	14,		// gps vcv 2
	14,		// gps vcv 3
	7,		// gps vscale
	7,		// gps pscale
	7,		// gps lscale
	7,		// gps hscale
	14,		// gps ref frame
	14,		// gps epoch
	4,		// angular d
	2,		// angular m
	8,		// angular s
	9,		// standard deviation
	7,		// inst height
	7,		// target height
	10,		// msr id
	10		// cluster id
};

class dna_msr_fields_300
	: public _dna_msr_fields_300_<UINT16>
{
public:
};


// DNA version 3.01
// Stations
template <class U>
struct _dna_stn_fields_301_
{
	static const U _locations_[8];
	static const U _widths_[8];
};

template <class U>
const U _dna_stn_fields_301_<U>::_locations_[8] = 
{
	0,		// station name
	20,		// constraints
	24,		// type
	27,		// e / phi / x
	47,		// n / lambda / y
	67,		// height / z
	87,		// hemisphere-zone
	91		// description
};

template <class U>
const U _dna_stn_fields_301_<U>::_widths_[8] = 
{
	20,		// station name
	3,		// constraints
	3,		// type
	20,		// e / phi / x
	20,		// n / lambda / y
	20,		// height / z
	3,		// hemisphere-zone
	128		// description
};

class dna_stn_fields_301
	: public _dna_stn_fields_301_<UINT16>
{
public:
};


// Measurements
template <class U>
struct _dna_msr_fields_301_
{
	static const U _locations_[24];
	static const U _widths_[24];
};

template <class U>
const U _dna_msr_fields_301_<U>::_locations_[24] = 
{
	0,		// msr type
	1,		// ignore 
	2,		// instrument name
	22,		// target 1
	42,		// target 2
	62,		// linear measurement
	62,		// gps measurement
	82,		// gps vcv 1
	102,	// gps vcv 2
	122,	// gps vcv 3
	62,		// gps vscale
	72,		// gps pscale
	82,		// gps lscale
	92,		// gps hscale
	102,	// gps ref frame
	122,	// gps epoch
	76,		// angular d
	80,		// angular m
	82,		// angular s
	90,		// standard deviation
	99,		// inst height
	106,	// target height
	142,	// msr id
	152		// cluster id
};

template <class U>
const U _dna_msr_fields_301_<U>::_widths_[24] = 
{
	1,		// msr type
	1,		// ignore 
	20,		// instrument name
	20,		// target 1
	20,		// target 2
	14,		// linear measurement
	20,		// gps measurement
	20,		// gps vcv 1
	20,		// gps vcv 2
	20,		// gps vcv 3
	10,		// gps vscale
	10,		// gps pscale
	10,		// gps lscale
	10,		// gps hscale
	20,		// gps ref frame
	20,		// gps epoch
	4,		// angular d
	2,		// angular m
	8,		// angular s
	9,		// standard deviation
	7,		// inst height
	7,		// target height
	10,		// msr id
	10		// cluster id
};

class dna_msr_fields_301
	: public _dna_msr_fields_301_<UINT16>
{
public:
};


// DNA version 3.02 - adds a 25th column for the observation epoch.
// Columns 0-23 are identical to version 3.01.
template <class U>
struct _dna_msr_fields_302_
{
	static const U _locations_[25];
	static const U _widths_[25];
};

template <class U>
const U _dna_msr_fields_302_<U>::_locations_[25] =
{
	0,		// msr type
	1,		// ignore
	2,		// instrument name
	22,		// target 1
	42,		// target 2
	62,		// linear measurement
	62,		// gps measurement
	82,		// gps vcv 1
	102,	// gps vcv 2
	122,	// gps vcv 3
	62,		// gps vscale
	72,		// gps pscale
	82,		// gps lscale
	92,		// gps hscale
	102,	// gps ref frame
	122,	// gps epoch (reference-frame epoch)
	76,		// angular d
	80,		// angular m
	82,		// angular s
	90,		// standard deviation
	99,		// inst height
	106,	// target height
	142,	// msr id
	152,	// cluster id
	162		// gps observation epoch
};

template <class U>
const U _dna_msr_fields_302_<U>::_widths_[25] =
{
	1,		// msr type
	1,		// ignore
	20,		// instrument name
	20,		// target 1
	20,		// target 2
	14,		// linear measurement
	20,		// gps measurement
	20,		// gps vcv 1
	20,		// gps vcv 2
	20,		// gps vcv 3
	10,		// gps vscale
	10,		// gps pscale
	10,		// gps lscale
	10,		// gps hscale
	20,		// gps ref frame
	20,		// gps epoch (reference-frame epoch)
	4,		// angular d
	2,		// angular m
	8,		// angular s
	9,		// standard deviation
	7,		// inst height
	7,		// target height
	10,		// msr id
	10,		// cluster id
	20		// gps observation epoch
};

class dna_msr_fields_302
	: public _dna_msr_fields_302_<UINT16>
{
public:
};



template <typename U>
void assignDNASTNFieldParameters(const UINT16* locs, const UINT16* widths, 
	dna_stn_fields& dflocs, dna_stn_fields& dfwidths)
{
	dflocs.stn_name = locs[0];
	dflocs.stn_const = locs[1];
	dflocs.stn_type = locs[2];
	dflocs.stn_e_phi_x = locs[3];
	dflocs.stn_n_lam_y = locs[4];
	dflocs.stn_ht_z = locs[5];
	dflocs.stn_hemi_zo = locs[6];
	dflocs.stn_desc = locs[7];
		
	dfwidths.stn_name = widths[0];
	dfwidths.stn_const = widths[1];
	dfwidths.stn_type = widths[2];
	dfwidths.stn_e_phi_x = widths[3];
	dfwidths.stn_n_lam_y = widths[4];
	dfwidths.stn_ht_z = widths[5];
	dfwidths.stn_hemi_zo = widths[6];
	dfwidths.stn_desc = widths[7];
}
	

template <typename U>
void determineDNASTNFieldParameters(const std::string& version, 
	dna_stn_fields& dflocs, dna_stn_fields& dfwidths)
{
	if (iequals(version, "3.01"))
	{
		assignDNASTNFieldParameters<U>(dna_stn_fields_301::_locations_,
			dna_stn_fields_301::_widths_,
			dflocs, dfwidths);
		return;
	}

	if (iequals(version, "3.00"))
	{
		assignDNASTNFieldParameters<U>(dna_stn_fields_300::_locations_,
			dna_stn_fields_300::_widths_,
			dflocs, dfwidths);
		return;
	}

	if (iequals(version, "1.00"))
	{
		assignDNASTNFieldParameters<U>(dna_stn_fields_100::_locations_,
			dna_stn_fields_100::_widths_,
			dflocs, dfwidths);
		return;
	}
}


template <typename U>
void assignDNAMSRFieldParameters(const UINT16* locs, const UINT16* widths, 
	dna_msr_fields& dflocs, dna_msr_fields& dfwidths)
{
	dflocs.msr_type = locs[0];
	dflocs.msr_ignore = locs[1];
	dflocs.msr_inst = locs[2];
	dflocs.msr_targ1 = locs[3];
	dflocs.msr_targ2 = locs[4];
	dflocs.msr_linear = locs[5];
	dflocs.msr_gps = locs[6];
	dflocs.msr_gps_vcv_1 = locs[7];
	dflocs.msr_gps_vcv_2 = locs[8];
	dflocs.msr_gps_vcv_3 = locs[9];
	dflocs.msr_gps_vscale = locs[10];
	dflocs.msr_gps_pscale = locs[11];
	dflocs.msr_gps_lscale = locs[12];
	dflocs.msr_gps_hscale = locs[13];
	dflocs.msr_gps_reframe = locs[14];
	dflocs.msr_gps_epoch = locs[15];
	dflocs.msr_ang_d = locs[16];
	dflocs.msr_ang_m = locs[17];
	dflocs.msr_ang_s = locs[18];
	dflocs.msr_stddev = locs[19];
	dflocs.msr_inst_ht = locs[20];
	dflocs.msr_targ_ht = locs[21];
	dflocs.msr_id_msr = locs[22];
	dflocs.msr_id_cluster = locs[23];
	// Default the 3.02-introduced field to 0; a v3.02 branch overrides it below.
	dflocs.msr_gps_obs_epoch = 0;

	dfwidths.msr_type = widths[0];
	dfwidths.msr_ignore = widths[1];
	dfwidths.msr_inst = widths[2];
	dfwidths.msr_targ1 = widths[3];
	dfwidths.msr_targ2 = widths[4];
	dfwidths.msr_linear = widths[5];
	dfwidths.msr_gps = widths[6];
	dfwidths.msr_gps_vcv_1 = widths[7];
	dfwidths.msr_gps_vcv_2 = widths[8];
	dfwidths.msr_gps_vcv_3 = widths[9];
	dfwidths.msr_gps_vscale = widths[10];
	dfwidths.msr_gps_pscale = widths[11];
	dfwidths.msr_gps_lscale = widths[12];
	dfwidths.msr_gps_hscale = widths[13];
	dfwidths.msr_gps_reframe = widths[14];
	dfwidths.msr_gps_epoch = widths[15];
	dfwidths.msr_ang_d = widths[16];
	dfwidths.msr_ang_m = widths[17];
	dfwidths.msr_ang_s = widths[18];
	dfwidths.msr_stddev = widths[19];
	dfwidths.msr_inst_ht = widths[20];
	dfwidths.msr_targ_ht = widths[21];
	dfwidths.msr_id_msr = widths[22];
	dfwidths.msr_id_cluster = widths[23];
	dfwidths.msr_gps_obs_epoch = 0;
}
	

template <typename U>
void determineDNAMSRFieldParameters(const std::string& version,
	dna_msr_fields& dflocs, dna_msr_fields& dfwidths, const U u=0)
{
	if (iequals(version, "3.02"))
	{
		assignDNAMSRFieldParameters<U>(dna_msr_fields_302::_locations_,
			dna_msr_fields_302::_widths_,
			dflocs, dfwidths);
		// The 3.02-introduced observation-epoch column is beyond the shared
		// 24-entry assigner, so set it explicitly here.
		dflocs.msr_gps_obs_epoch = dna_msr_fields_302::_locations_[24];
		dfwidths.msr_gps_obs_epoch = dna_msr_fields_302::_widths_[24];
		return;
	}

	if (iequals(version, "3.01"))
	{
		assignDNAMSRFieldParameters<U>(dna_msr_fields_301::_locations_,
			dna_msr_fields_301::_widths_,
			dflocs, dfwidths);
		return;
	}

	if (iequals(version, "3.00"))
	{
		assignDNAMSRFieldParameters<U>(dna_msr_fields_300::_locations_,
			dna_msr_fields_300::_widths_,
			dflocs, dfwidths);
		return;
	}

	if (iequals(version, "1.00"))
	{
		assignDNAMSRFieldParameters<U>(dna_msr_fields_100::_locations_,
			dna_msr_fields_100::_widths_,
			dflocs, dfwidths);
		return;
	}

	std::stringstream ss;
	ss << " Unknown file version: " << version << std::endl;
	throw XMLInteropException(ss.str(), u);
}


#endif // DNAIODNATYPES_H_
