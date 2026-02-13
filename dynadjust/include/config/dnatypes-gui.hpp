//============================================================================
// Name         : dnatypes-gui.hpp
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
// Description  : DynAdjust data types for enums used in the gui include file
//============================================================================

#pragma once
#ifndef DNATYPES_GUI_H_
#define DNATYPES_GUI_H_

#if defined(_MSC_VER)
	#if defined(LIST_INCLUDES_ON_BUILD) 
		#pragma message("  " __FILE__) 
	#endif
#endif

typedef enum APU_UNITS_UI
{
	XYZ_apu_ui = 0,		// Cartesian
	ENU_apu_ui = 1,		// Local
	LLH_apu_ui = 2,		// Polar (LLH)
} _APU_UNITS_UI_;

typedef enum ADJ_GNSS_UNITS_UI
{
	XYZ_adj_gnss_ui = 0,		// Cartesian
	ENU_adj_gnss_ui = 1,		// Local
	AED_adj_gnss_ui = 2,		// Polar (Az, elev, dist)
	ADU_adj_gnss_ui = 3			// Polar (Az, dist, ht)
} _ADJ_GNSS_UNITS_UI_;

typedef enum ADJ_MSR_SORT_UI
{
	orig_adj_msr_sort_ui = 0,		// original file order
	type_adj_msr_sort_ui = 1,		// measurement type
	inst_adj_msr_sort_ui = 2,		// instrument station
	targ_adj_msr_sort_ui = 3,
	meas_adj_msr_sort_ui = 4,
	corr_adj_msr_sort_ui = 5,
	a_sd_adj_msr_sort_ui = 6,
	n_st_adj_msr_sort_ui = 7
} _ADJ_MSR_SORT_UI_;

typedef enum MSR_TO_STN_SORT_UI
{
	orig_stn_sort_ui = 0,		// original station order
    name_stn_sort_ui = 1,		// station name
	meas_stn_sort_ui = 2,		// measurement count (ascending)
    saem_stn_sort_ui = 3		// measurement count (descending)
} _MSR_TO_STN_SORT_UI_;


#endif









