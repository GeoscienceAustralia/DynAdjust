//============================================================================
// Name         : dnatypes-basic.hpp
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
// Description  : Basic data types and fundamental definitions for DynAdjust
//============================================================================

#pragma once
#ifndef DNATYPES_BASIC_H_
#define DNATYPES_BASIC_H_

#include <cstdint>

// Basic type definitions
#ifdef UINT32
#undef UINT32
#endif

#ifdef UINT16
#undef UINT16
#endif

typedef unsigned int	UINT32, *PUINT32;
typedef const PUINT32	PUINT32_const;
typedef unsigned short	UINT16, *PUINT16;

typedef uint32_t        index_t;

// Basic constants
#ifndef INT_MIN
#define INT_MIN -2147483648
#endif

#ifndef TRUE
#define TRUE 1
#endif

#ifndef FALSE
#define FALSE 0
#endif

#ifndef DNATYPES_FWD_H_
const char LOWER_TRIANGLE = 'L';
const char UPPER_TRIANGLE = 'U';

const UINT16 VALID_STATION = 1;
const UINT16 INVALID_STATION = 0;
#endif

// Station dimensions constants
const UINT16 STN_NAME_WIDTH(31);
const UINT16 STN_NAME_ORIG_WIDTH(40);
const UINT16 STN_DESC_WIDTH(129);
const UINT16 STN_CONST_WIDTH(4);
const UINT16 STN_TYPE_WIDTH(4);
const UINT16 STN_EPSG_WIDTH(7);
const UINT16 STN_EPOCH_WIDTH(12);
const UINT16 STN_PLATE_WIDTH(3);

const UINT32 MOD_NAME_WIDTH(20);
const UINT32 FILE_NAME_WIDTH(256);

// Enumerations
typedef enum _SIGMA_ZERO_STAT_PASS_
{
	test_stat_pass = 0,
	test_stat_warning = 1,
	test_stat_fail = 2
} SIGMA_ZERO_STAT_PASS;

typedef enum _INPUT_FILE_TYPE_
{
	geodesyml = 0,
	dynaml = 1,
	dna = 2,
	csv = 3,
	sinex = 4,
	jsonl = 5
} INPUT_FILE_TYPE;

typedef enum _INPUT_DATA_TYPE_
{
	stn_data = 0,
	msr_data = 1,
	stn_msr_data = 2,
	geo_data = 3,
	ren_data = 4,
	tbu_data = 5,
	unknown = 6
} INPUT_DATA_TYPE;

typedef enum _TIMER_TYPE_
{
	iteration_time = 0,
	total_time = 1
} TIMER_TYPE;

typedef enum _ANGULAR_TYPE_
{
	DMS = 0,
	//#define DMIN 1
	DDEG = 1
} ANGULAR_TYPE;

typedef enum _DMS_FORMAT_
{
	SEPARATED = 0,
	SEPARATED_WITH_SYMBOLS = 1,
	HP_NOTATION = 2
} DMS_FORMAT;

typedef enum _COORD_TYPE_
{
	XYZ_type_i = 0,		// Cartesian
	LLh_type_i = 1,		// Geographic (ellipsoid height)
	LLH_type_i = 2,		// Geographic (orthometric height)
	UTM_type_i = 3,		// Projection
	ENU_type_i = 4,		// Local
	AED_type_i = 5		// Azimuth, elevation and distance
} COORD_TYPE;

typedef enum _STATION_ELEM_
{
	station_1 = 0,
	station_2 = 1,
	station_3 = 2
} STATION_ELEM;

typedef enum _CART_ELEM_
{
	x_element = 0,
	y_element = 1,
	z_element = 2
} CART_ELEM;

typedef enum _COORDINATE_TYPES_
{
	latitude_t = 0,
	longitude_t = 1,
	easting_t = 2,
	northing_t = 3
} COORDINATE_TYPES;

typedef enum _HEIGHT_SYSTEM_
{
	ORTHOMETRIC_type_i = 0,
	ELLIPSOIDAL_type_i = 1
} HEIGHT_SYSTEM;

typedef enum _MEASUREMENT_STATIONS_ {
	ONE_STATION = 1,
	TWO_STATION = 2,
	THREE_STATION = 3,
	UNKNOWN_TYPE = -1
} MEASUREMENT_STATIONS;

typedef enum _MEASUREMENT_START_ {
	xMeas = 0,
	yMeas = 1,
	zMeas = 2,
	xCov = 3,
	yCov = 4,
	zCov = 5
} MEASUREMENT_START;

typedef enum _CONSTRAINT_TYPE_ {
	free_3D = 0,			// FFF - positional or 3d adjustment
	constrained_3D = 1,		// CCC
	free_2D = 2,			// FFC - horizontal or 2d adjustment. Implies vertical is tightly constrained (or fixed).
	free_1D = 3,			// CCF - vertical or 1d adjustment. Implies horizontal is tightly constrained (or fixed).
	custom_constraint = 4   // A custom constraint type that isn't necessarily
} CONSTRAINT_TYPE;

typedef enum _AML_TYPE_ {
	str_msr = 0,	// msr_t struct
	cls_msr = 1		// CDnaMeasurement class
} AML_TYPE;

enum iosMode
{
	binary = 0,
	ascii = 1
};

enum mtxType
{
	// mtx_full is intended for non-square matrices which are
	// commonly fully populated, such as a vector of estimates.
	// In this case, the whole buffer is copied/stored.
	mtx_full = 0,
	// mtx_lower is intended for square matrices, such as
	// full variance matrices and normal equations.
	// In this case, data is stored in columns, from left to right
	mtx_lower = 1,
	// mtx_sparse is intended for matrices which have many
	// zeros, such as design and AtV-1 matrices.
	// In this case, each element is stored with its index
	mtx_sparse = 2
};

typedef enum _STAGE_FILE_
{
	sf_normals = 0,
	sf_normals_r = 1,
	sf_atvinv = 2,
	sf_design = 3,
	sf_meas_minus_comp = 4,
	sf_estimated_stns = 5,
	sf_original_stns = 6,
	sf_rigorous_stns = 7,
	sf_junction_vars = 8,
	sf_junction_vars_f = 9,
	sf_junction_ests_f = 10,
	sf_junction_ests_r = 11,
	sf_rigorous_vars = 12,
	sf_prec_adj_msrs = 13,
	sf_corrections = 14
} STAGE_FILE;

#endif // DNATYPES_BASIC_H_