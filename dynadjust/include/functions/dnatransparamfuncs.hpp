//============================================================================
// Name         : dnatransparamfuncs.hpp
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
// Description  : Functions to determine the correct 7/14 parameters between
//                input and output frames.
//============================================================================

#ifndef DNATRANS_PARAM_FUNCS_H_
#define DNATRANS_PARAM_FUNCS_H_

#if defined(_MSC_VER)
	#if defined(LIST_INCLUDES_ON_BUILD) 
		#pragma message("  " __FILE__)
	#endif
#endif

/// \cond
#include <string>
#include <iomanip>
#include <sstream>
/// \endcond

#include <include/parameters/dnaepsg.hpp>
#include <include/parameters/dnatransformationparameters.hpp>

using namespace dynadjust::epsg;

namespace dynadjust {
namespace datum_parameters {

///////////////////////////////////////////////////////////////////////////////////////////
template <class S, class U>
S message_parameters_undefined(U from_epsg_code, U to_epsg_code)
{
	// No direct parameters exist!
	std::stringstream ss;
	ss << datumFromEpsgCode<S, U>(from_epsg_code) << " to " <<
		datumFromEpsgCode<S, U>(to_epsg_code) << " parameters have not been defined.";
	return ss.str();
}

///////////////////////////////////////////////////////////////////////////////////////////

template <typename U>
void determineGDA94Parameters(transformation_parameter_set& tParam)
{
	// Which 'to' datum?
	switch (tParam.from_to_.second)
	{
	///////////////////////////////////////////////////////////////////////////////////////////////
	// GDA94-to-GDA2020
	case GDA2020_i_xyz:
	case GDA2020_i_2d:
	case GDA2020_i:
		memcpy(&tParam.parameters_, GDA94_GDA2020<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = GDA94_GDA2020<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = GDA94_GDA2020<double, UINT32>::reference_frame;
		break;
	// GDA94-to-GDA94
	case GDA94_i_xyz:
	case GDA94_i_2d:
	case GDA94_i:
		// no transformation - set zero
		memset(&tParam.parameters_, 0, sizeof(double) * 14);
		break;
	///////////////////////////////////////////////////////////////////////////////////////////////
	// GDA94-to-ITRF1988
	case ITRF1988_i_xyz:
	case ITRF1988_i:
	// GDA94-to-ITRF1989
	case ITRF1989_i_xyz:
	case ITRF1989_i:
	// GDA94-to-ITRF1990
	case ITRF1990_i_xyz:
	case ITRF1990_i:
	// GDA94-to-ITRF1991
	case ITRF1991_i_xyz:
	case ITRF1991_i:
	// GDA94-to-ITRF1992
	case ITRF1992_i_xyz:
	case ITRF1992_i:
	// GDA94-to-ITRF1993
	case ITRF1993_i_xyz:
	case ITRF1993_i:
	// GDA94-to-ITRF1994
	case ITRF1994_i_xyz:
	case ITRF1994_i:
		// No direct parameters exist!
		throw RefTranException(message_parameters_undefined<std::string, UINT32>(tParam.from_to_.first, tParam.from_to_.second),
			REFTRAN_DIRECT_PARAMS_UNAVAILABLE);
	// GDA94-to-ITRF1996
	case ITRF1996_i_xyz:
	case ITRF1996_i:
		memcpy(&tParam.parameters_, ITRF1996_GDA94<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF1996_GDA94<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF1996_GDA94<double, UINT32>::reference_frame;
		tParam.paramDirection_ = __paramReverse__;
		tParam.reverse();
		break;
	// GDA94-to-ITRF1997
	case ITRF1997_i_xyz:
	case ITRF1997_i:
		memcpy(&tParam.parameters_, ITRF1997_GDA94<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF1997_GDA94<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF1997_GDA94<double, UINT32>::reference_frame;
		tParam.paramDirection_ = __paramReverse__;
		tParam.reverse();
		break;
	// GDA94-to-ITRF2000
	case ITRF2000_i_xyz:
	case ITRF2000_i:
		memcpy(&tParam.parameters_, ITRF2000_GDA94<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2000_GDA94<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2000_GDA94<double, UINT32>::reference_frame;
		tParam.paramDirection_ = __paramReverse__;
		tParam.reverse();
		break;
	// GDA94-to-ITRF2005
	case ITRF2005_i_xyz:
	case ITRF2005_i:
		memcpy(&tParam.parameters_, ITRF2005_GDA94<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2005_GDA94<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2005_GDA94<double, UINT32>::reference_frame;
		tParam.paramDirection_ = __paramReverse__;
		tParam.reverse();
		break;
	// GDA94-to-ITRF2008
	case ITRF2008_i_xyz:
	case ITRF2008_i:
		memcpy(&tParam.parameters_, ITRF2008_GDA94<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2008_GDA94<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2008_GDA94<double, UINT32>::reference_frame;
		tParam.paramDirection_ = __paramReverse__;
		tParam.reverse();
		break;
	// GDA94-to-ITRF2014
	case ITRF2014_i_xyz:
	case ITRF2014_i:
		memcpy(&tParam.parameters_, ITRF2014_GDA94<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2014_GDA94<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2014_GDA94<double, UINT32>::reference_frame;
		tParam.paramDirection_ = __paramReverse__;
		tParam.reverse();
		break;
	// GDA94-to-ITRF2020
	case ITRF2020_i_xyz:
	case ITRF2020_i:
		// No direct parameters exist!
		throw RefTranException(message_parameters_undefined<std::string, UINT32>(tParam.from_to_.first, tParam.from_to_.second),
			REFTRAN_DIRECT_PARAMS_UNAVAILABLE);
	// GDA94-to-WGS84
	case WGS84_i_xyz:
	case WGS84_i:
	case WGS84_transit_i_xyz:
	case WGS84_transit_i:
	case WGS84_G730_i_xyz:
	case WGS84_G730_i:
	case WGS84_G873_i_xyz:
	case WGS84_G873_i:
	case WGS84_G1150_i_xyz:
	case WGS84_G1150_i:
	case WGS84_G1674_i_xyz:
	case WGS84_G1674_i:
	case WGS84_G1762_i_xyz:
	case WGS84_G1762_i:
	case WGS84_G2139_i_xyz:
	case WGS84_G2139_i:
	default:
		throw std::runtime_error("determineGDA94Parameters: ");
	}
}
	

template <typename U>
void determineGDA2020Parameters(transformation_parameter_set& tParam)
{
	// Which 'to' datum?
	switch (tParam.from_to_.second)
	{
	///////////////////////////////////////////////////////////////////////////////////////////////
	// GDA2020-to-GDA94
	case GDA94_i_xyz:
	case GDA94_i_2d:
	case GDA94_i:
		memcpy(&tParam.parameters_, GDA94_GDA2020<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = GDA94_GDA2020<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = GDA94_GDA2020<double, UINT32>::reference_frame;
		tParam.paramDirection_ = __paramReverse__;
		tParam.reverse();
		break;
	// GDA2020-to-GDA2020
	case GDA2020_i_xyz:
	case GDA2020_i_2d:
	case GDA2020_i:
		// no transformation - set zero
		memset(&tParam.parameters_, 0, sizeof(double) * 14);
		break;
	///////////////////////////////////////////////////////////////////////////////////////////////
	// GDA2020-to-ITRF1988
	case ITRF1988_i_xyz:
	case ITRF1988_i:
	// GDA2020-to-ITRF1989
	case ITRF1989_i_xyz:
	case ITRF1989_i:
	// GDA2020-to-ITRF1990
	case ITRF1990_i_xyz:
	case ITRF1990_i:
	// GDA2020-to-ITRF1991
	case ITRF1991_i_xyz:
	case ITRF1991_i:
	// GDA2020-to-ITRF1992
	case ITRF1992_i_xyz:
	case ITRF1992_i:
	// GDA2020-to-ITRF1993
	case ITRF1993_i_xyz:
	case ITRF1993_i:
	// GDA2020-to-ITRF1994
	case ITRF1994_i_xyz:
	case ITRF1994_i:
		// No direct parameters exist!
		throw RefTranException(message_parameters_undefined<std::string, UINT32>(tParam.from_to_.first, tParam.from_to_.second),
			REFTRAN_DIRECT_PARAMS_UNAVAILABLE);
	// GDA2020-to-ITRF1996
	case ITRF1996_i_xyz:
	case ITRF1996_i:
		memcpy(&tParam.parameters_, ITRF1996_GDA2020<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF1996_GDA2020<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF1996_GDA2020<double, UINT32>::reference_frame;
		tParam.paramDirection_ = __paramReverse__;
		tParam.reverse();
		break;
	// GDA2020-to-ITRF1997
	case ITRF1997_i_xyz:
	case ITRF1997_i:
		memcpy(&tParam.parameters_, ITRF1997_GDA2020<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF1997_GDA2020<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF1997_GDA2020<double, UINT32>::reference_frame;
		tParam.paramDirection_ = __paramReverse__;
		tParam.reverse();
		break;
	// GDA2020-to-ITRF2000
	case ITRF2000_i_xyz:
	case ITRF2000_i:
		memcpy(&tParam.parameters_, ITRF2000_GDA2020<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2000_GDA2020<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2000_GDA2020<double, UINT32>::reference_frame;
		tParam.paramDirection_ = __paramReverse__;
		tParam.reverse();
		break;
	// GDA2020-to-ITRF2005
	case ITRF2005_i_xyz:
	case ITRF2005_i:
		memcpy(&tParam.parameters_, ITRF2005_GDA2020<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2005_GDA2020<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2005_GDA2020<double, UINT32>::reference_frame;
		tParam.paramDirection_ = __paramReverse__;
		tParam.reverse();
		break;
	// GDA2020-to-ITRF2008
	case ITRF2008_i_xyz:
	case ITRF2008_i:
		memcpy(&tParam.parameters_, ITRF2008_GDA2020<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2008_GDA2020<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2008_GDA2020<double, UINT32>::reference_frame;
		tParam.paramDirection_ = __paramReverse__;
		tParam.reverse();
		break;
	// GDA2020-to-ITRF2014
	case ITRF2014_i_xyz:
	case ITRF2014_i:
		memcpy(&tParam.parameters_, ITRF2014_GDA2020<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2014_GDA2020<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2014_GDA2020<double, UINT32>::reference_frame;
		tParam.paramDirection_ = __paramReverse__;
		tParam.reverse();
		break;
	// GDA2020-to-ITRF2020
	case ITRF2020_i_xyz:
	case ITRF2020_i:
		memcpy(&tParam.parameters_, ITRF2020_GDA2020<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2020_GDA2020<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2020_GDA2020<double, UINT32>::reference_frame;
		tParam.paramDirection_ = __paramReverse__;
		tParam.reverse();
		break;
	// GDA2020-to-WGS84
	case WGS84_i_xyz:
	case WGS84_i:
	case WGS84_transit_i_xyz:
	case WGS84_transit_i:
	case WGS84_G730_i_xyz:
	case WGS84_G730_i:
	case WGS84_G873_i_xyz:
	case WGS84_G873_i:
	case WGS84_G1150_i_xyz:
	case WGS84_G1150_i:
	case WGS84_G1674_i_xyz:
	case WGS84_G1674_i:
	case WGS84_G1762_i_xyz:
	case WGS84_G1762_i:
	case WGS84_G2139_i_xyz:
	case WGS84_G2139_i:
	default:
		throw std::runtime_error("determineGDA2020Parameters: ");
	}
}
	

template <typename U>
void determineITRF1988Parameters(transformation_parameter_set& tParam)
{
	// Which 'to' datum?
	switch (tParam.from_to_.second)
	{
	//////////////////////////////////
	// ITRF1988-to-GDA2020
	case GDA2020_i_xyz:
	case GDA2020_i_2d:
	case GDA2020_i:
	// ITRF1988-to-GDA94
	case GDA94_i_xyz:
	case GDA94_i_2d:
	case GDA94_i:
	//////////////////////////////////
	// ITRF1988-to-ITRF1989
	case ITRF1989_i_xyz:
	case ITRF1989_i:
	// ITRF1988-to-ITRF1990
	case ITRF1990_i_xyz:
	case ITRF1990_i:
	// ITRF1988-to-ITRF1991
	case ITRF1991_i_xyz:
	case ITRF1991_i:
	// ITRF1988-to-ITRF1992
	case ITRF1992_i_xyz:
	case ITRF1992_i:
	// ITRF1988-to-ITRF1993
	case ITRF1993_i_xyz:
	case ITRF1993_i:
	// ITRF1988-to-ITRF1994
	case ITRF1994_i_xyz:
	case ITRF1994_i:
	// ITRF1988-to-ITRF1996
	case ITRF1996_i_xyz:
	case ITRF1996_i:
	// ITRF1988-to-ITRF1997
	case ITRF1997_i_xyz:
	case ITRF1997_i:
	// ITRF1988-to-ITRF2005
	case ITRF2005_i_xyz:
	case ITRF2005_i:
		// No direct parameters exist!
		throw RefTranException(message_parameters_undefined<std::string, UINT32>(tParam.from_to_.first, tParam.from_to_.second),
			REFTRAN_DIRECT_PARAMS_UNAVAILABLE);
	// ITRF1988-to-ITRF1988
	case ITRF1988_i_xyz:
	case ITRF1988_i:
		// not a transformation of frame, but a shift in epoch only
		// Not defined yet - requires interpolation from the ITRF1992 velocity model
		throw RefTranException("A transformation between different epochs on the same frame is required.",
			REFTRAN_TRANS_ON_PLATE_REQUIRED);
	// ITRF1988-to-ITRF2000
	case ITRF2000_i_xyz:
	case ITRF2000_i:
		memcpy(&tParam.parameters_, ITRF2000_ITRF1988<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2000_ITRF1988<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2000_ITRF1988<double, UINT32>::reference_frame;
		tParam.paramDirection_ = __paramReverse__;
		tParam.reverse();
		break;
	// ITRF1988-to-ITRF2008
	case ITRF2008_i_xyz:
	case ITRF2008_i:
		memcpy(&tParam.parameters_, ITRF2008_ITRF1988<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2008_ITRF1988<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2008_ITRF1988<double, UINT32>::reference_frame;
		tParam.paramDirection_ = __paramReverse__;
		tParam.reverse();
		break;
	// ITRF1988-to-ITRF2014
	case ITRF2014_i_xyz:
	case ITRF2014_i:
		memcpy(&tParam.parameters_, ITRF2014_ITRF1988<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2014_ITRF1988<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2014_ITRF1988<double, UINT32>::reference_frame;
		tParam.paramDirection_ = __paramReverse__;
		tParam.reverse();
		break;
	// ITRF1988-to-ITRF2020
	case ITRF2020_i_xyz:
	case ITRF2020_i:
		memcpy(&tParam.parameters_, ITRF2020_ITRF1988<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2020_ITRF1988<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2020_ITRF1988<double, UINT32>::reference_frame;
		tParam.paramDirection_ = __paramReverse__;
		tParam.reverse();
		break;
	// ITRF1988-to-WGS84
	case WGS84_i_xyz:
	case WGS84_i:
	case WGS84_transit_i_xyz:
	case WGS84_transit_i:
	case WGS84_G730_i_xyz:
	case WGS84_G730_i:
	case WGS84_G873_i_xyz:
	case WGS84_G873_i:
	case WGS84_G1150_i_xyz:
	case WGS84_G1150_i:
	case WGS84_G1674_i_xyz:
	case WGS84_G1674_i:
	case WGS84_G1762_i_xyz:
	case WGS84_G1762_i:
	case WGS84_G2139_i_xyz:
	case WGS84_G2139_i:
	default:
		throw std::runtime_error("determineITRF1988Parameters: ");
	}
}

template <typename U>
void determineITRF1989Parameters(transformation_parameter_set& tParam)
{
	// Which 'to' datum?
	switch (tParam.from_to_.second)
	{
	//////////////////////////////////
	// ITRF1989-to-GDA2020
	case GDA2020_i_xyz:
	case GDA2020_i_2d:
	case GDA2020_i:
	// ITRF1989-to-GDA94
	case GDA94_i_xyz:
	case GDA94_i_2d:
	case GDA94_i:
	//////////////////////////////////
	// ITRF1989-to-ITRF1988
	case ITRF1988_i_xyz:
	case ITRF1988_i:
	// ITRF1989-to-ITRF1990
	case ITRF1990_i_xyz:
	case ITRF1990_i:
	// ITRF1989-to-ITRF1991
	case ITRF1991_i_xyz:
	case ITRF1991_i:
	// ITRF1989-to-ITRF1992
	case ITRF1992_i_xyz:
	case ITRF1992_i:
	// ITRF1989-to-ITRF1993
	case ITRF1993_i_xyz:
	case ITRF1993_i:
	// ITRF1989-to-ITRF1994
	case ITRF1994_i_xyz:
	case ITRF1994_i:
	// ITRF1989-to-ITRF1996
	case ITRF1996_i_xyz:
	case ITRF1996_i:
	// ITRF1989-to-ITRF1997
	case ITRF1997_i_xyz:
	case ITRF1997_i:
	// ITRF1989-to-ITRF2005
	case ITRF2005_i_xyz:
	case ITRF2005_i:
		// No direct parameters exist!
		throw RefTranException(message_parameters_undefined<std::string, UINT32>(tParam.from_to_.first, tParam.from_to_.second),
			REFTRAN_DIRECT_PARAMS_UNAVAILABLE);
	// ITRF1989-to-ITRF1989
	case ITRF1989_i_xyz:
	case ITRF1989_i:
		// not a transformation of frame, but a shift in epoch only
		// Not defined yet - requires interpolation from the ITRF1992 velocity model
		throw RefTranException("A transformation between different epochs on the same frame is required.", 
			REFTRAN_TRANS_ON_PLATE_REQUIRED);
	// ITRF1989-to-ITRF2000
	case ITRF2000_i_xyz:
	case ITRF2000_i:
		memcpy(&tParam.parameters_, ITRF2000_ITRF1989<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2000_ITRF1989<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2000_ITRF1989<double, UINT32>::reference_frame;
		tParam.paramDirection_ = __paramReverse__;
		tParam.reverse();
		break;
	// ITRF1989-to-ITRF2008
	case ITRF2008_i_xyz:
	case ITRF2008_i:
		memcpy(&tParam.parameters_, ITRF2008_ITRF1989<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2008_ITRF1989<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2008_ITRF1989<double, UINT32>::reference_frame;
		tParam.paramDirection_ = __paramReverse__;
		tParam.reverse();
		break;
	// ITRF1989-to-ITRF2014
	case ITRF2014_i_xyz:
	case ITRF2014_i:
		memcpy(&tParam.parameters_, ITRF2014_ITRF1989<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2014_ITRF1989<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2014_ITRF1989<double, UINT32>::reference_frame;
		tParam.paramDirection_ = __paramReverse__;
		tParam.reverse();
		break;
	// ITRF1989-to-ITRF2020
	case ITRF2020_i_xyz:
	case ITRF2020_i:
		memcpy(&tParam.parameters_, ITRF2020_ITRF1989<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2020_ITRF1989<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2020_ITRF1989<double, UINT32>::reference_frame;
		tParam.paramDirection_ = __paramReverse__;
		tParam.reverse();
		break;
	// ITRF1989-to-WGS84
	case WGS84_i_xyz:
	case WGS84_i:
	case WGS84_transit_i_xyz:
	case WGS84_transit_i:
	case WGS84_G730_i_xyz:
	case WGS84_G730_i:
	case WGS84_G873_i_xyz:
	case WGS84_G873_i:
	case WGS84_G1150_i_xyz:
	case WGS84_G1150_i:
	case WGS84_G1674_i_xyz:
	case WGS84_G1674_i:
	case WGS84_G1762_i_xyz:
	case WGS84_G1762_i:
	case WGS84_G2139_i_xyz:
	case WGS84_G2139_i:
	default:		
		throw std::runtime_error("determineITRF1989Parameters: ");
	}
}

template <typename U>
void determineITRF1990Parameters(transformation_parameter_set& tParam)
{
	// Which 'to' datum?
	switch (tParam.from_to_.second)
	{
	//////////////////////////////////
	// ITRF1990-to-GDA2020
	case GDA2020_i_xyz:
	case GDA2020_i_2d:
	case GDA2020_i:
	// ITRF1990-to-GDA94
	case GDA94_i_xyz:
	case GDA94_i_2d:
	case GDA94_i:
	//////////////////////////////////
	// ITRF1990-to-ITRF1988
	case ITRF1988_i_xyz:
	case ITRF1988_i:
	// ITRF1990-to-ITRF1989
	case ITRF1989_i_xyz:
	case ITRF1989_i:
	// ITRF1990-to-ITRF1991
	case ITRF1991_i_xyz:
	case ITRF1991_i:
	// ITRF1990-to-ITRF1992
	case ITRF1992_i_xyz:
	case ITRF1992_i:
	// ITRF1990-to-ITRF1993
	case ITRF1993_i_xyz:
	case ITRF1993_i:
	// ITRF1990-to-ITRF1994
	case ITRF1994_i_xyz:
	case ITRF1994_i:
	// ITRF1990-to-ITRF1996
	case ITRF1996_i_xyz:
	case ITRF1996_i:
	// ITRF1990-to-ITRF1997
	case ITRF1997_i_xyz:
	case ITRF1997_i:
	// ITRF1990-to-ITRF2005
	case ITRF2005_i_xyz:
	case ITRF2005_i:
		// No direct parameters exist!
		throw RefTranException(message_parameters_undefined<std::string, UINT32>(tParam.from_to_.first, tParam.from_to_.second),
			REFTRAN_DIRECT_PARAMS_UNAVAILABLE);
	// ITRF1990-to-ITRF1990
	case ITRF1990_i_xyz:
	case ITRF1990_i:
		// not a transformation of frame, but a shift in epoch only
		// Not defined yet - requires interpolation from the ITRF1990 velocity model
		throw RefTranException("A transformation between different epochs on the same frame is required.", 
			REFTRAN_TRANS_ON_PLATE_REQUIRED);
	// ITRF1990-to-ITRF2000
	case ITRF2000_i_xyz:
	case ITRF2000_i:
		memcpy(&tParam.parameters_, ITRF2000_ITRF1990<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2000_ITRF1990<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2000_ITRF1990<double, UINT32>::reference_frame;
		tParam.paramDirection_ = __paramReverse__;
		tParam.reverse();
		break;
	// ITRF1990-to-ITRF2008
	case ITRF2008_i_xyz:
	case ITRF2008_i:
		memcpy(&tParam.parameters_, ITRF2008_ITRF1990<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2008_ITRF1990<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2008_ITRF1990<double, UINT32>::reference_frame;
		tParam.paramDirection_ = __paramReverse__;
		tParam.reverse();
		break;
	// ITRF1990-to-ITRF2014
	case ITRF2014_i_xyz:
	case ITRF2014_i:
		memcpy(&tParam.parameters_, ITRF2014_ITRF1990<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2014_ITRF1990<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2014_ITRF1990<double, UINT32>::reference_frame;
		tParam.paramDirection_ = __paramReverse__;
		tParam.reverse();
		break;
	// ITRF1990-to-ITRF2020
	case ITRF2020_i_xyz:
	case ITRF2020_i:
		memcpy(&tParam.parameters_, ITRF2020_ITRF1990<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2020_ITRF1990<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2020_ITRF1990<double, UINT32>::reference_frame;
		tParam.paramDirection_ = __paramReverse__;
		tParam.reverse();
		break;
	// ITRF1990-to-WGS84
	case WGS84_i_xyz:
	case WGS84_i:
	case WGS84_transit_i_xyz:
	case WGS84_transit_i:
	case WGS84_G730_i_xyz:
	case WGS84_G730_i:
	case WGS84_G873_i_xyz:
	case WGS84_G873_i:
	case WGS84_G1150_i_xyz:
	case WGS84_G1150_i:
	case WGS84_G1674_i_xyz:
	case WGS84_G1674_i:
	case WGS84_G1762_i_xyz:
	case WGS84_G1762_i:
	case WGS84_G2139_i_xyz:
	case WGS84_G2139_i:
	default:
		throw std::runtime_error("determineITRF1990Parameters: ");
	}
}

template <typename U>
void determineITRF1991Parameters(transformation_parameter_set& tParam)
{
	// Which 'to' datum?
	switch (tParam.from_to_.second)
	{
	//////////////////////////////////
	// ITRF1991-to-GDA2020
	case GDA2020_i_xyz:
	case GDA2020_i_2d:
	case GDA2020_i:
	// ITRF1991-to-GDA94
	case GDA94_i_xyz:
	case GDA94_i_2d:
	case GDA94_i:
	//////////////////////////////////
	// ITRF1991-to-ITRF1988
	case ITRF1988_i_xyz:
	case ITRF1988_i:
	// ITRF1991-to-ITRF1989
	case ITRF1989_i_xyz:
	case ITRF1989_i:
	// ITRF1991-to-ITRF1990
	case ITRF1990_i_xyz:
	case ITRF1990_i:
	// ITRF1991-to-ITRF1992
	case ITRF1992_i_xyz:
	case ITRF1992_i:
	// ITRF1991-to-ITRF1993
	case ITRF1993_i_xyz:
	case ITRF1993_i:
	// ITRF1991-to-ITRF1994
	case ITRF1994_i_xyz:
	case ITRF1994_i:
	// ITRF1991-to-ITRF1996
	case ITRF1996_i_xyz:
	case ITRF1996_i:
	// ITRF1991-to-ITRF1997
	case ITRF1997_i_xyz:
	case ITRF1997_i:
	// ITRF1991-to-ITRF2005
	case ITRF2005_i_xyz:
	case ITRF2005_i:
		// No direct parameters exist!
		throw RefTranException(message_parameters_undefined<std::string, UINT32>(tParam.from_to_.first, tParam.from_to_.second),
			REFTRAN_DIRECT_PARAMS_UNAVAILABLE);
	// ITRF1991-to-ITRF1991
	case ITRF1991_i_xyz:
	case ITRF1991_i:
		// not a transformation of frame, but a shift in epoch only
		// Not defined yet - requires interpolation from the ITRF1991 velocity model
		throw RefTranException("A transformation between different epochs on the same frame is required.", 
			REFTRAN_TRANS_ON_PLATE_REQUIRED);
	// ITRF1991-to-ITRF2000
	case ITRF2000_i_xyz:
	case ITRF2000_i:
		memcpy(&tParam.parameters_, ITRF2000_ITRF1991<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2000_ITRF1991<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2000_ITRF1991<double, UINT32>::reference_frame;
		tParam.paramDirection_ = __paramReverse__;
		tParam.reverse();
		break;
	// ITRF1991-to-ITRF2008
	case ITRF2008_i_xyz:
	case ITRF2008_i:
		memcpy(&tParam.parameters_, ITRF2008_ITRF1991<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2008_ITRF1991<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2008_ITRF1991<double, UINT32>::reference_frame;
		tParam.paramDirection_ = __paramReverse__;
		tParam.reverse();
		break;
	// ITRF1991-to-ITRF2014
	case ITRF2014_i_xyz:
	case ITRF2014_i:
		memcpy(&tParam.parameters_, ITRF2014_ITRF1991<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2014_ITRF1991<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2014_ITRF1991<double, UINT32>::reference_frame;
		tParam.paramDirection_ = __paramReverse__;
		tParam.reverse();
		break;
	// ITRF1991-to-ITRF2020
	case ITRF2020_i_xyz:
	case ITRF2020_i:
		memcpy(&tParam.parameters_, ITRF2020_ITRF1991<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2020_ITRF1991<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2020_ITRF1991<double, UINT32>::reference_frame;
		tParam.paramDirection_ = __paramReverse__;
		tParam.reverse();
		break;
	// ITRF1991-to-WGS84
	case WGS84_i_xyz:
	case WGS84_i:
	case WGS84_transit_i_xyz:
	case WGS84_transit_i:
	case WGS84_G730_i_xyz:
	case WGS84_G730_i:
	case WGS84_G873_i_xyz:
	case WGS84_G873_i:
	case WGS84_G1150_i_xyz:
	case WGS84_G1150_i:
	case WGS84_G1674_i_xyz:
	case WGS84_G1674_i:
	case WGS84_G1762_i_xyz:
	case WGS84_G1762_i:
	case WGS84_G2139_i_xyz:
	case WGS84_G2139_i:
	default:		
		throw std::runtime_error("determineITRF1991Parameters: ");
	}
}

template <typename U>
void determineITRF1992Parameters(transformation_parameter_set& tParam)
{
	// Which 'to' datum?
	switch (tParam.from_to_.second)
	{
	//////////////////////////////////
	// ITRF1992-to-GDA2020
	case GDA2020_i_xyz:
	case GDA2020_i_2d:
	case GDA2020_i:
	// ITRF1992-to-GDA94
	case GDA94_i_xyz:
	case GDA94_i_2d:
	case GDA94_i:
	//////////////////////////////////
	// ITRF1992-to-ITRF1988
	case ITRF1988_i_xyz:
	case ITRF1988_i:
	// ITRF1992-to-ITRF1989
	case ITRF1989_i_xyz:
	case ITRF1989_i:
	// ITRF1992-to-ITRF1990
	case ITRF1990_i_xyz:
	case ITRF1990_i:
	// ITRF1992-to-ITRF1991
	case ITRF1991_i_xyz:
	case ITRF1991_i:
	// ITRF1992-to-ITRF1993
	case ITRF1993_i_xyz:
	case ITRF1993_i:
	// ITRF1992-to-ITRF1994
	case ITRF1994_i_xyz:
	case ITRF1994_i:
	// ITRF1992-to-ITRF1996
	case ITRF1996_i_xyz:
	case ITRF1996_i:
	// ITRF1992-to-ITRF1997
	case ITRF1997_i_xyz:
	case ITRF1997_i:
	// ITRF1992-to-ITRF2005
	case ITRF2005_i_xyz:
	case ITRF2005_i:
		// No direct parameters exist!
		throw RefTranException(message_parameters_undefined<std::string, UINT32>(tParam.from_to_.first, tParam.from_to_.second),
			REFTRAN_DIRECT_PARAMS_UNAVAILABLE);
	// ITRF1992-to-ITRF1992
	case ITRF1992_i_xyz:
	case ITRF1992_i:
		// not a transformation of frame, but a shift in epoch only
		// Not defined yet - requires interpolation from the ITRF1992 velocity model
		throw RefTranException("A transformation between different epochs on the same frame is required.", 
			REFTRAN_TRANS_ON_PLATE_REQUIRED);
	// ITRF1992-to-ITRF2000
	case ITRF2000_i_xyz:
	case ITRF2000_i:
		memcpy(&tParam.parameters_, ITRF2000_ITRF1992<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2000_ITRF1992<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2000_ITRF1992<double, UINT32>::reference_frame;
		tParam.paramDirection_ = __paramReverse__;
		tParam.reverse();
		break;
	// ITRF1992-to-ITRF2008
	case ITRF2008_i_xyz:
	case ITRF2008_i:
		memcpy(&tParam.parameters_, ITRF2008_ITRF1992<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2008_ITRF1992<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2008_ITRF1992<double, UINT32>::reference_frame;
		tParam.paramDirection_ = __paramReverse__;
		tParam.reverse();
		break;
	// ITRF1992-to-ITRF2014
	case ITRF2014_i_xyz:
	case ITRF2014_i:
		memcpy(&tParam.parameters_, ITRF2014_ITRF1992<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2014_ITRF1992<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2014_ITRF1992<double, UINT32>::reference_frame;
		tParam.paramDirection_ = __paramReverse__;
		tParam.reverse();
		break;
	// ITRF1992-to-ITRF2020
	case ITRF2020_i_xyz:
	case ITRF2020_i:
		memcpy(&tParam.parameters_, ITRF2020_ITRF1992<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2020_ITRF1992<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2020_ITRF1992<double, UINT32>::reference_frame;
		tParam.paramDirection_ = __paramReverse__;
		tParam.reverse();
		break;
	// ITRF1992-to-WGS84
	case WGS84_i_xyz:
	case WGS84_i:
	case WGS84_transit_i_xyz:
	case WGS84_transit_i:
	case WGS84_G730_i_xyz:
	case WGS84_G730_i:
	case WGS84_G873_i_xyz:
	case WGS84_G873_i:
	case WGS84_G1150_i_xyz:
	case WGS84_G1150_i:
	case WGS84_G1674_i_xyz:
	case WGS84_G1674_i:
	case WGS84_G1762_i_xyz:
	case WGS84_G1762_i:
	case WGS84_G2139_i_xyz:
	case WGS84_G2139_i:
	default:
		throw std::runtime_error("determineITRF1992Parameters: ");
	}
}

template <typename U>
void determineITRF1993Parameters(transformation_parameter_set& tParam)
{
	// Which 'to' datum?
	switch (tParam.from_to_.second)
	{
	//////////////////////////////////
	// ITRF1993-to-GDA2020
	case GDA2020_i_xyz:
	case GDA2020_i_2d:
	case GDA2020_i:
	// ITRF1993-to-GDA94
	case GDA94_i_xyz:
	case GDA94_i_2d:
	case GDA94_i:
	//////////////////////////////////
	// ITRF1993-to-ITRF1988
	case ITRF1988_i_xyz:
	case ITRF1988_i:
	// ITRF1993-to-ITRF1989
	case ITRF1989_i_xyz:
	case ITRF1989_i:
	// ITRF1993-to-ITRF1990
	case ITRF1990_i_xyz:
	case ITRF1990_i:
	// ITRF1993-to-ITRF1991
	case ITRF1991_i_xyz:
	case ITRF1991_i:
	// ITRF1993-to-ITRF1992
	case ITRF1992_i_xyz:
	case ITRF1992_i:
	// ITRF1993-to-ITRF1994
	case ITRF1994_i_xyz:
	case ITRF1994_i:
	// ITRF1993-to-ITRF1996
	case ITRF1996_i_xyz:
	case ITRF1996_i:
	// ITRF1993-to-ITRF1997
	case ITRF1997_i_xyz:
	case ITRF1997_i:
	// ITRF1993-to-ITRF2005
	case ITRF2005_i_xyz:
	case ITRF2005_i:
		// No direct parameters exist!
		throw RefTranException(message_parameters_undefined<std::string, UINT32>(tParam.from_to_.first, tParam.from_to_.second),
			REFTRAN_DIRECT_PARAMS_UNAVAILABLE);
	// ITRF1993-to-ITRF1993
	case ITRF1993_i_xyz:
	case ITRF1993_i:
		// not a transformation of frame, but a shift in epoch only
		// Not defined yet - requires interpolation from the ITRF1993 velocity model
		throw RefTranException("A transformation between different epochs on the same frame is required.", 
			REFTRAN_TRANS_ON_PLATE_REQUIRED);
	// ITRF1993-to-ITRF2000
	case ITRF2000_i_xyz:
	case ITRF2000_i:
		memcpy(&tParam.parameters_, ITRF2000_ITRF1993<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2000_ITRF1993<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2000_ITRF1993<double, UINT32>::reference_frame;
		tParam.paramDirection_ = __paramReverse__;
		tParam.reverse();
		break;
	// ITRF1993-to-ITRF2008
	case ITRF2008_i_xyz:
	case ITRF2008_i:
		memcpy(&tParam.parameters_, ITRF2008_ITRF1993<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2008_ITRF1993<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2008_ITRF1993<double, UINT32>::reference_frame;
		tParam.paramDirection_ = __paramReverse__;
		tParam.reverse();
		break;
	// ITRF1993-to-ITRF2014
	case ITRF2014_i_xyz:
	case ITRF2014_i:
		memcpy(&tParam.parameters_, ITRF2014_ITRF1993<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2014_ITRF1993<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2014_ITRF1993<double, UINT32>::reference_frame;
		tParam.paramDirection_ = __paramReverse__;
		tParam.reverse();
		break;
	// ITRF1993-to-ITRF2020
	case ITRF2020_i_xyz:
	case ITRF2020_i:
		memcpy(&tParam.parameters_, ITRF2020_ITRF1993<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2020_ITRF1993<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2020_ITRF1993<double, UINT32>::reference_frame;
		tParam.paramDirection_ = __paramReverse__;
		tParam.reverse();
		break;
	// ITRF1993-to-WGS84
	case WGS84_i_xyz:
	case WGS84_i:
	case WGS84_transit_i_xyz:
	case WGS84_transit_i:
	case WGS84_G730_i_xyz:
	case WGS84_G730_i:
	case WGS84_G873_i_xyz:
	case WGS84_G873_i:
	case WGS84_G1150_i_xyz:
	case WGS84_G1150_i:
	case WGS84_G1674_i_xyz:
	case WGS84_G1674_i:
	case WGS84_G1762_i_xyz:
	case WGS84_G1762_i:
	case WGS84_G2139_i_xyz:
	case WGS84_G2139_i:
	default:
		throw std::runtime_error("determineITRF1993Parameters: ");
	}
}

template <typename U>
void determineITRF1994Parameters(transformation_parameter_set& tParam)
{
	// Which 'to' datum?
	switch (tParam.from_to_.second)
	{
	//////////////////////////////////
	// ITRF1994-to-GDA2020
	case GDA2020_i_xyz:
	case GDA2020_i_2d:
	case GDA2020_i:
	// ITRF1994-to-GDA94
	case GDA94_i_xyz:
	case GDA94_i_2d:
	case GDA94_i:
	//////////////////////////////////
	// ITRF1994-to-ITRF1988
	case ITRF1988_i_xyz:
	case ITRF1988_i:
	// ITRF1994-to-ITRF1989
	case ITRF1989_i_xyz:
	case ITRF1989_i:
	// ITRF1994-to-ITRF1990
	case ITRF1990_i_xyz:
	case ITRF1990_i:
	// ITRF1994-to-ITRF1991
	case ITRF1991_i_xyz:
	case ITRF1991_i:
	// ITRF1994-to-ITRF1992
	case ITRF1992_i_xyz:
	case ITRF1992_i:
	// ITRF1994-to-ITRF1993
	case ITRF1993_i_xyz:
	case ITRF1993_i:
	// ITRF1994-to-ITRF1996
	case ITRF1996_i_xyz:
	case ITRF1996_i:
	// ITRF1994-to-ITRF1997
	case ITRF1997_i_xyz:
	case ITRF1997_i:
	// ITRF1994-to-ITRF2005
	case ITRF2005_i_xyz:
	case ITRF2005_i:
		// No direct parameters exist!
		throw RefTranException(message_parameters_undefined<std::string, UINT32>(tParam.from_to_.first, tParam.from_to_.second),
			REFTRAN_DIRECT_PARAMS_UNAVAILABLE);
	// ITRF1994-to-ITRF1994
	case ITRF1994_i_xyz:
	case ITRF1994_i:
		// not a transformation of frame, but a shift in epoch only
		// Not defined yet - requires interpolation from the ITRF1994 velocity model
		throw RefTranException("A transformation between different epochs on the same frame is required.", 
			REFTRAN_TRANS_ON_PLATE_REQUIRED);
	// ITRF1994-to-ITRF2000
	case ITRF2000_i_xyz:
	case ITRF2000_i:
		memcpy(&tParam.parameters_, ITRF2000_ITRF1994<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2000_ITRF1994<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2000_ITRF1994<double, UINT32>::reference_frame;
		tParam.paramDirection_ = __paramReverse__;
		tParam.reverse();
		break;
	// ITRF1994-to-ITRF2008
	case ITRF2008_i_xyz:
	case ITRF2008_i:
		memcpy(&tParam.parameters_, ITRF2008_ITRF1994<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2008_ITRF1994<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2008_ITRF1994<double, UINT32>::reference_frame;
		tParam.paramDirection_ = __paramReverse__;
		tParam.reverse();
		break;
	// ITRF1994-to-ITRF2014
	case ITRF2014_i_xyz:
	case ITRF2014_i:
		memcpy(&tParam.parameters_, ITRF2014_ITRF1994<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2014_ITRF1994<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2014_ITRF1994<double, UINT32>::reference_frame;
		tParam.paramDirection_ = __paramReverse__;
		tParam.reverse();
		break;
	// ITRF1994-to-ITRF2020
	case ITRF2020_i_xyz:
	case ITRF2020_i:
		memcpy(&tParam.parameters_, ITRF2020_ITRF1994<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2020_ITRF1994<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2020_ITRF1994<double, UINT32>::reference_frame;
		tParam.paramDirection_ = __paramReverse__;
		tParam.reverse();
		break;
	// ITRF1994-to-WGS84
	case WGS84_i_xyz:
	case WGS84_i:
	case WGS84_transit_i_xyz:
	case WGS84_transit_i:
	case WGS84_G730_i_xyz:
	case WGS84_G730_i:
	case WGS84_G873_i_xyz:
	case WGS84_G873_i:
	case WGS84_G1150_i_xyz:
	case WGS84_G1150_i:
	case WGS84_G1674_i_xyz:
	case WGS84_G1674_i:
	case WGS84_G1762_i_xyz:
	case WGS84_G1762_i:
	case WGS84_G2139_i_xyz:
	case WGS84_G2139_i:
	default:
		throw std::runtime_error("determineITRF1994Parameters: ");
	}
}

template <typename U>
void determineITRF1996Parameters(transformation_parameter_set& tParam)
{
	// Which 'to' datum?
	switch (tParam.from_to_.second)
	{
	//////////////////////////////////
	// ITRF1996-to-GDA2020
	case GDA2020_i_xyz:
	case GDA2020_i_2d:
	case GDA2020_i:
		memcpy(&tParam.parameters_, ITRF1996_GDA2020<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF1996_GDA2020<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF1996_GDA2020<double, UINT32>::reference_frame;
		break;
	// ITRF1996-to-GDA94
	case GDA94_i_xyz:
	case GDA94_i_2d:
	case GDA94_i:
		memcpy(&tParam.parameters_, ITRF1996_GDA94<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF1996_GDA94<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF1996_GDA94<double, UINT32>::reference_frame;
		break;
	//////////////////////////////////
	// ITRF1996-to-ITRF1988
	case ITRF1988_i_xyz:
	case ITRF1988_i:
	// ITRF1996-to-ITRF1989
	case ITRF1989_i_xyz:
	case ITRF1989_i:
	// ITRF1996-to-ITRF1990
	case ITRF1990_i_xyz:
	case ITRF1990_i:
	// ITRF1996-to-ITRF1991
	case ITRF1991_i_xyz:
	case ITRF1991_i:
	// ITRF1996-to-ITRF1992
	case ITRF1992_i_xyz:
	case ITRF1992_i:
	// ITRF1996-to-ITRF1993
	case ITRF1993_i_xyz:
	case ITRF1993_i:
	// ITRF1996-to-ITRF1994
	case ITRF1994_i_xyz:
	case ITRF1994_i:
	// ITRF1996-to-ITRF1997
	case ITRF1997_i_xyz:
	case ITRF1997_i:
	// ITRF1996-to-ITRF2005
	case ITRF2005_i_xyz:
	case ITRF2005_i:
		// No direct parameters exist!
		throw RefTranException(message_parameters_undefined<std::string, UINT32>(tParam.from_to_.first, tParam.from_to_.second),
			REFTRAN_DIRECT_PARAMS_UNAVAILABLE);
	// ITRF1996-to-ITRF1996
	case ITRF1996_i_xyz:
	case ITRF1996_i:
		// not a transformation of frame, but a shift in epoch only
		// Not defined yet - requires interpolation from the ITRF1996 velocity model
		throw RefTranException("A transformation between different epochs on the same frame is required.", 
			REFTRAN_TRANS_ON_PLATE_REQUIRED);
	// ITRF1996-to-ITRF2000
	case ITRF2000_i_xyz:
	case ITRF2000_i:
		memcpy(&tParam.parameters_, ITRF2000_ITRF1996<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2000_ITRF1996<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2000_ITRF1996<double, UINT32>::reference_frame;
		tParam.paramDirection_ = __paramReverse__;
		tParam.reverse();
		break;
	// ITRF1996-to-ITRF2008
	case ITRF2008_i_xyz:
	case ITRF2008_i:
		memcpy(&tParam.parameters_, ITRF2008_ITRF1996<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2008_ITRF1996<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2008_ITRF1996<double, UINT32>::reference_frame;
		tParam.paramDirection_ = __paramReverse__;
		tParam.reverse();
		break;
	// ITRF1996-to-ITRF2014
	case ITRF2014_i_xyz:
	case ITRF2014_i:
		memcpy(&tParam.parameters_, ITRF2014_ITRF1996<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2014_ITRF1996<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2014_ITRF1996<double, UINT32>::reference_frame;
		tParam.paramDirection_ = __paramReverse__;
		tParam.reverse();
		break;
	// ITRF1996-to-ITRF2020
	case ITRF2020_i_xyz:
	case ITRF2020_i:
		memcpy(&tParam.parameters_, ITRF2020_ITRF1996<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2020_ITRF1996<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2020_ITRF1996<double, UINT32>::reference_frame;
		tParam.paramDirection_ = __paramReverse__;
		tParam.reverse();
		break;
	// ITRF1996-to-WGS84
	case WGS84_i_xyz:
	case WGS84_i:
	case WGS84_transit_i_xyz:
	case WGS84_transit_i:
	case WGS84_G730_i_xyz:
	case WGS84_G730_i:
	case WGS84_G873_i_xyz:
	case WGS84_G873_i:
	case WGS84_G1150_i_xyz:
	case WGS84_G1150_i:
	case WGS84_G1674_i_xyz:
	case WGS84_G1674_i:
	case WGS84_G1762_i_xyz:
	case WGS84_G1762_i:
	case WGS84_G2139_i_xyz:
	case WGS84_G2139_i:
	default:
		throw std::runtime_error("determineITRF1996Parameters: ");
	}
}

template <typename U>
void determineITRF1997Parameters(transformation_parameter_set& tParam)
{
	// Which 'to' datum?
	switch (tParam.from_to_.second)
	{
	//////////////////////////////////
	// ITRF1997-to-GDA2020
	case GDA2020_i_xyz:
	case GDA2020_i_2d:
	case GDA2020_i:
		memcpy(&tParam.parameters_, ITRF1997_GDA2020<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF1997_GDA2020<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF1997_GDA2020<double, UINT32>::reference_frame;
		break;
	// ITRF1997-to-GDA94
	case GDA94_i_xyz:
	case GDA94_i_2d:
	case GDA94_i:
		memcpy(&tParam.parameters_, ITRF1997_GDA94<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF1997_GDA94<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF1997_GDA94<double, UINT32>::reference_frame;
		break;
	//////////////////////////////////
	// ITRF1997-to-ITRF1988
	case ITRF1988_i_xyz:
	case ITRF1988_i:
	// ITRF1997-to-ITRF1989
	case ITRF1989_i_xyz:
	case ITRF1989_i:
	// ITRF1997-to-ITRF1990
	case ITRF1990_i_xyz:
	case ITRF1990_i:
	// ITRF1997-to-ITRF1991
	case ITRF1991_i_xyz:
	case ITRF1991_i:
	// ITRF1997-to-ITRF1992
	case ITRF1992_i_xyz:
	case ITRF1992_i:
	// ITRF1997-to-ITRF1993
	case ITRF1993_i_xyz:
	case ITRF1993_i:
	// ITRF1997-to-ITRF1994
	case ITRF1994_i_xyz:
	case ITRF1994_i:
	// ITRF1997-to-ITRF1996
	case ITRF1996_i_xyz:
	case ITRF1996_i:
	// ITRF1997-to-ITRF2005
	case ITRF2005_i_xyz:
	case ITRF2005_i:
		// No direct parameters exist!
		throw RefTranException(message_parameters_undefined<std::string, UINT32>(tParam.from_to_.first, tParam.from_to_.second),
			REFTRAN_DIRECT_PARAMS_UNAVAILABLE);
	// ITRF1997-to-ITRF1997
	case ITRF1997_i_xyz:
	case ITRF1997_i:
		// not a transformation of frame, but a shift in epoch only
		// Not defined yet - requires interpolation from the ITRF1997 velocity model
		throw RefTranException("A transformation between different epochs on the same frame is required.", 
			REFTRAN_TRANS_ON_PLATE_REQUIRED);
	// ITRF1997-to-ITRF2000
	case ITRF2000_i_xyz:
	case ITRF2000_i:
		memcpy(&tParam.parameters_, ITRF2000_ITRF1997<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2000_ITRF1997<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2000_ITRF1997<double, UINT32>::reference_frame;
		tParam.paramDirection_ = __paramReverse__;
		tParam.reverse();
		break;
	// ITRF1997-to-ITRF2008
	case ITRF2008_i_xyz:
	case ITRF2008_i:
		memcpy(&tParam.parameters_, ITRF2008_ITRF1997<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2008_ITRF1997<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2008_ITRF1997<double, UINT32>::reference_frame;
		tParam.paramDirection_ = __paramReverse__;
		tParam.reverse();
		break;
	// ITRF1997-to-ITRF2014
	case ITRF2014_i_xyz:
	case ITRF2014_i:
		memcpy(&tParam.parameters_, ITRF2014_ITRF1997<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2014_ITRF1997<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2014_ITRF1997<double, UINT32>::reference_frame;
		tParam.paramDirection_ = __paramReverse__;
		tParam.reverse();
		break;
	// ITRF1997-to-ITRF2020
	case ITRF2020_i_xyz:
	case ITRF2020_i:
		memcpy(&tParam.parameters_, ITRF2020_ITRF1997<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2020_ITRF1997<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2020_ITRF1997<double, UINT32>::reference_frame;
		tParam.paramDirection_ = __paramReverse__;
		tParam.reverse();
		break;
	// ITRF1997-to-WGS84
	case WGS84_i_xyz:
	case WGS84_i:
	case WGS84_transit_i_xyz:
	case WGS84_transit_i:
	case WGS84_G730_i_xyz:
	case WGS84_G730_i:
	case WGS84_G873_i_xyz:
	case WGS84_G873_i:
	case WGS84_G1150_i_xyz:
	case WGS84_G1150_i:
	case WGS84_G1674_i_xyz:
	case WGS84_G1674_i:
	case WGS84_G1762_i_xyz:
	case WGS84_G1762_i:
	case WGS84_G2139_i_xyz:
	case WGS84_G2139_i:
	default:		
		throw std::runtime_error("determineITRF1997Parameters: ");
	}
}

template <typename U>
void determineITRF2000Parameters(transformation_parameter_set& tParam)
{
	// Which 'to' datum?
	switch (tParam.from_to_.second)
	{
	//////////////////////////////////
	// ITRF2000-to-GDA2020
	case GDA2020_i_xyz:
	case GDA2020_i_2d:
	case GDA2020_i:
		memcpy(&tParam.parameters_, ITRF2000_GDA2020<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2000_GDA2020<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2000_GDA2020<double, UINT32>::reference_frame;
		break;
	// ITRF2000-to-GDA94
	case GDA94_i_xyz:
	case GDA94_i_2d:
	case GDA94_i:
		memcpy(&tParam.parameters_, ITRF2000_GDA94<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2000_GDA94<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2000_GDA94<double, UINT32>::reference_frame;
		break;
	//////////////////////////////////
	// ITRF2000-to-ITRF1988
	case ITRF1988_i_xyz:
	case ITRF1988_i:
		memcpy(&tParam.parameters_, ITRF2000_ITRF1988<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2000_ITRF1988<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2000_ITRF1988<double, UINT32>::reference_frame;
		break;
	// ITRF2000-to-ITRF1989
	case ITRF1989_i_xyz:
	case ITRF1989_i:
		memcpy(&tParam.parameters_, ITRF2000_ITRF1989<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2000_ITRF1989<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2000_ITRF1989<double, UINT32>::reference_frame;
		break;
	// ITRF2000-to-ITRF1990
	case ITRF1990_i_xyz:
	case ITRF1990_i:
		memcpy(&tParam.parameters_, ITRF2000_ITRF1990<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2000_ITRF1990<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2000_ITRF1990<double, UINT32>::reference_frame;
		break;
	// ITRF2000-to-ITRF1991
	case ITRF1991_i_xyz:
	case ITRF1991_i:
		memcpy(&tParam.parameters_, ITRF2000_ITRF1991<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2000_ITRF1991<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2000_ITRF1991<double, UINT32>::reference_frame;
		break;
	// ITRF2000-to-ITRF1992
	case ITRF1992_i_xyz:
	case ITRF1992_i:
		memcpy(&tParam.parameters_, ITRF2000_ITRF1992<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2000_ITRF1992<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2000_ITRF1992<double, UINT32>::reference_frame;
		break;
	// ITRF2000-to-ITRF1993
	case ITRF1993_i_xyz:
	case ITRF1993_i:
		memcpy(&tParam.parameters_, ITRF2000_ITRF1993<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2000_ITRF1993<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2000_ITRF1993<double, UINT32>::reference_frame;
		break;
	// ITRF2000-to-ITRF1994
	case ITRF1994_i_xyz:
	case ITRF1994_i:
		memcpy(&tParam.parameters_, ITRF2000_ITRF1994<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2000_ITRF1994<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2000_ITRF1994<double, UINT32>::reference_frame;
		break;
	// ITRF2000-to-ITRF1996
	case ITRF1996_i_xyz:
	case ITRF1996_i:
		memcpy(&tParam.parameters_, ITRF2000_ITRF1996<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2000_ITRF1996<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2000_ITRF1996<double, UINT32>::reference_frame;
		break;
	// ITRF2000-to-ITRF1997
	case ITRF1997_i_xyz:
	case ITRF1997_i:
		memcpy(&tParam.parameters_, ITRF2000_ITRF1997<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2000_ITRF1997<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2000_ITRF1997<double, UINT32>::reference_frame;
		break;
	// ITRF2000-to-ITRF2000
	case ITRF2000_i_xyz:
	case ITRF2000_i:
		// not a transformation of frame, but a shift in epoch only
		// Not defined yet - requires interpolation from the ITRF2005 velocity model
		throw RefTranException("A transformation between different epochs on the same frame is required.", 
			REFTRAN_TRANS_ON_PLATE_REQUIRED);
	// ITRF2000-to-ITRF2005
	case ITRF2005_i_xyz:
	case ITRF2005_i:
		memcpy(&tParam.parameters_, ITRF2005_ITRF2000<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2005_ITRF2000<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2005_ITRF2000<double, UINT32>::reference_frame;
		tParam.paramDirection_ = __paramReverse__;
		tParam.reverse();
		break;
	// ITRF2000-to-ITRF2008
	case ITRF2008_i_xyz:
	case ITRF2008_i:
		memcpy(&tParam.parameters_, ITRF2008_ITRF2000<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2008_ITRF2000<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2008_ITRF2000<double, UINT32>::reference_frame;
		tParam.paramDirection_ = __paramReverse__;
		tParam.reverse();
		break;
	// ITRF2000-to-ITRF2014
	case ITRF2014_i_xyz:
	case ITRF2014_i:
		memcpy(&tParam.parameters_, ITRF2014_ITRF2000<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2014_ITRF2000<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2014_ITRF2000<double, UINT32>::reference_frame;
		tParam.paramDirection_ = __paramReverse__;
		tParam.reverse();
		break;
	// ITRF2000-to-ITRF2020
	case ITRF2020_i_xyz:
	case ITRF2020_i:
		memcpy(&tParam.parameters_, ITRF2020_ITRF2000<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2020_ITRF2000<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2020_ITRF2000<double, UINT32>::reference_frame;
		tParam.paramDirection_ = __paramReverse__;
		tParam.reverse();
		break;
	// ITRF2000-to-WGS84
	case WGS84_i_xyz:
	case WGS84_i:
	case WGS84_transit_i_xyz:
	case WGS84_transit_i:
	case WGS84_G730_i_xyz:
	case WGS84_G730_i:
	case WGS84_G873_i_xyz:
	case WGS84_G873_i:
	case WGS84_G1150_i_xyz:
	case WGS84_G1150_i:
	case WGS84_G1674_i_xyz:
	case WGS84_G1674_i:
	case WGS84_G1762_i_xyz:
	case WGS84_G1762_i:
	case WGS84_G2139_i_xyz:
	case WGS84_G2139_i:
	default:
		throw std::runtime_error("determineITRF2000Parameters: ");
	}
}

template <typename U>
void determineITRF2005Parameters(transformation_parameter_set& tParam)
{
	// Which 'to' datum?
	switch (tParam.from_to_.second)
	{
	//////////////////////////////////
	// ITRF2005-to-GDA2020
	case GDA2020_i_xyz:
	case GDA2020_i_2d:
	case GDA2020_i:
		memcpy(&tParam.parameters_, ITRF2005_GDA2020<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2005_GDA2020<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2005_GDA2020<double, UINT32>::reference_frame;
		break;
	// ITRF2005-to-GDA94
	case GDA94_i_xyz:
	case GDA94_i_2d:
	case GDA94_i:
		memcpy(&tParam.parameters_, ITRF2005_GDA94<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2005_GDA94<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2005_GDA94<double, UINT32>::reference_frame;
		break;
	//////////////////////////////////
	// ITRF2005-to-ITRF1988
	case ITRF1988_i_xyz:
	case ITRF1988_i:
	// ITRF2005-to-ITRF1989
	case ITRF1989_i_xyz:
	case ITRF1989_i:
	// ITRF2005-to-ITRF1990
	case ITRF1990_i_xyz:
	case ITRF1990_i:
	// ITRF2005-to-ITRF1991
	case ITRF1991_i_xyz:
	case ITRF1991_i:
	// ITRF2005-to-ITRF1992
	case ITRF1992_i_xyz:
	case ITRF1992_i:
	// ITRF2005-to-ITRF1993
	case ITRF1993_i_xyz:
	case ITRF1993_i:
	// ITRF2005-to-ITRF1994
	case ITRF1994_i_xyz:
	case ITRF1994_i:
	// ITRF2005-to-ITRF1996
	case ITRF1996_i_xyz:
	case ITRF1996_i:
	// ITRF2005-to-ITRF1997
	case ITRF1997_i_xyz:
	case ITRF1997_i:
		// No direct parameters exist!
		throw RefTranException(message_parameters_undefined<std::string, UINT32>(tParam.from_to_.first, tParam.from_to_.second),
			REFTRAN_DIRECT_PARAMS_UNAVAILABLE);
	// ITRF2005-to-ITRF2000
	case ITRF2000_i_xyz:
	case ITRF2000_i:
		memcpy(&tParam.parameters_, ITRF2005_ITRF2000<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2005_ITRF2000<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2005_ITRF2000<double, UINT32>::reference_frame;
		break;
	// ITRF2005-to-ITRF2005
	case ITRF2005_i_xyz:
	case ITRF2005_i:
		// not a transformation of frame, but a shift in epoch only
		// Not defined yet - requires interpolation from the ITRF2005 velocity model
		throw RefTranException("A transformation between different epochs on the same frame is required.",
			REFTRAN_TRANS_ON_PLATE_REQUIRED);
	// ITRF2005-to-ITRF2008
	case ITRF2008_i_xyz:
	case ITRF2008_i:
		memcpy(&tParam.parameters_, ITRF2008_ITRF2005<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2008_ITRF2005<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2008_ITRF2005<double, UINT32>::reference_frame;
		tParam.paramDirection_ = __paramReverse__;
		tParam.reverse();
		break;
	// ITRF2005-to-ITRF2014
	case ITRF2014_i_xyz:
	case ITRF2014_i:
		memcpy(&tParam.parameters_, ITRF2014_ITRF2005<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2014_ITRF2005<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2014_ITRF2005<double, UINT32>::reference_frame;
		tParam.paramDirection_ = __paramReverse__;
		tParam.reverse();
		break;
	// ITRF2005-to-ITRF2020
	case ITRF2020_i_xyz:
	case ITRF2020_i:
		memcpy(&tParam.parameters_, ITRF2020_ITRF2005<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2020_ITRF2005<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2020_ITRF2005<double, UINT32>::reference_frame;
		tParam.paramDirection_ = __paramReverse__;
		tParam.reverse();
		break;
	// ITRF2005-to-WGS84
	case WGS84_i_xyz:
	case WGS84_i:
	case WGS84_transit_i_xyz:
	case WGS84_transit_i:
	case WGS84_G730_i_xyz:
	case WGS84_G730_i:
	case WGS84_G873_i_xyz:
	case WGS84_G873_i:
	case WGS84_G1150_i_xyz:
	case WGS84_G1150_i:
	case WGS84_G1674_i_xyz:
	case WGS84_G1674_i:
	case WGS84_G1762_i_xyz:
	case WGS84_G1762_i:
	case WGS84_G2139_i_xyz:
	case WGS84_G2139_i:
	default:
		throw std::runtime_error("determineITRF2005Parameters: ");
	}
}

template <typename U>
void determineITRF2008Parameters(transformation_parameter_set& tParam)
{
	// Which 'to' datum?
	switch (tParam.from_to_.second)
	{
	//////////////////////////////////
	// ITRF2008-to-GDA2020
	case GDA2020_i_xyz:
	case GDA2020_i_2d:
	case GDA2020_i:
		memcpy(&tParam.parameters_, ITRF2008_GDA2020<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2008_GDA2020<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2008_GDA2020<double, UINT32>::reference_frame;
		break;
	// ITRF2008-to-GDA94
	case GDA94_i_xyz:
	case GDA94_i_2d:
	case GDA94_i:
		memcpy(&tParam.parameters_, ITRF2008_GDA94<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2008_GDA94<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2008_GDA94<double, UINT32>::reference_frame;
		break;
	//////////////////////////////////
	// ITRF2008-to-ITRF1988
	case ITRF1988_i_xyz:
	case ITRF1988_i:
		memcpy(&tParam.parameters_, ITRF2008_ITRF1988<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2008_ITRF1988<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2008_ITRF1988<double, UINT32>::reference_frame;
		break;
	// ITRF2008-to-ITRF1989
	case ITRF1989_i_xyz:
	case ITRF1989_i:
		memcpy(&tParam.parameters_, ITRF2008_ITRF1989<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2008_ITRF1989<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2008_ITRF1989<double, UINT32>::reference_frame;
		break;
	// ITRF2008-to-ITRF1990
	case ITRF1990_i_xyz:
	case ITRF1990_i:
		memcpy(&tParam.parameters_, ITRF2008_ITRF1990<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2008_ITRF1990<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2008_ITRF1990<double, UINT32>::reference_frame;
		break;
	// ITRF2008-to-ITRF1991
	case ITRF1991_i_xyz:
	case ITRF1991_i:
		memcpy(&tParam.parameters_, ITRF2008_ITRF1991<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2008_ITRF1991<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2008_ITRF1991<double, UINT32>::reference_frame;
		break;
	// ITRF2008-to-ITRF1992
	case ITRF1992_i_xyz:
	case ITRF1992_i:
		memcpy(&tParam.parameters_, ITRF2008_ITRF1992<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2008_ITRF1992<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2008_ITRF1992<double, UINT32>::reference_frame;
		break;
	// ITRF2008-to-ITRF1993
	case ITRF1993_i_xyz:
	case ITRF1993_i:
		memcpy(&tParam.parameters_, ITRF2008_ITRF1993<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2008_ITRF1993<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2008_ITRF1993<double, UINT32>::reference_frame;
		break;
	// ITRF2008-to-ITRF1994
	case ITRF1994_i_xyz:
	case ITRF1994_i:
		memcpy(&tParam.parameters_, ITRF2008_ITRF1994<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2008_ITRF1994<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2008_ITRF1994<double, UINT32>::reference_frame;
		break;
	// ITRF2008-to-ITRF1996
	case ITRF1996_i_xyz:
	case ITRF1996_i:
		memcpy(&tParam.parameters_, ITRF2008_ITRF1996<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2008_ITRF1996<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2008_ITRF1996<double, UINT32>::reference_frame;
		break;
	// ITRF2008-to-ITRF1997
	case ITRF1997_i_xyz:
	case ITRF1997_i:
		memcpy(&tParam.parameters_, ITRF2008_ITRF1997<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2008_ITRF1997<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2008_ITRF1997<double, UINT32>::reference_frame;
		break;
	// ITRF2008-to-ITRF2000
	case ITRF2000_i_xyz:
	case ITRF2000_i:
		memcpy(&tParam.parameters_, ITRF2008_ITRF2000<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2008_ITRF2000<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2008_ITRF2000<double, UINT32>::reference_frame;
		break;
	// ITRF2008-to-ITRF2005
	case ITRF2005_i_xyz:
	case ITRF2005_i:
		memcpy(&tParam.parameters_, ITRF2008_ITRF2005<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2008_ITRF2005<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2008_ITRF2005<double, UINT32>::reference_frame;
		break;
	// ITRF2008-to-ITRF2008
	case ITRF2008_i_xyz:
	case ITRF2008_i:
		// not a transformation of frame, but a shift in epoch only
		// Not defined yet - requires interpolation from the ITRF2008 velocity model
		throw RefTranException("A transformation between different epochs on the same frame is required.", 
			REFTRAN_TRANS_ON_PLATE_REQUIRED);
	// ITRF2008-to-ITRF2014
	case ITRF2014_i_xyz:
	case ITRF2014_i:
		memcpy(&tParam.parameters_, ITRF2014_ITRF2008<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2014_ITRF2008<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2014_ITRF2008<double, UINT32>::reference_frame;
		tParam.paramDirection_ = __paramReverse__;
		tParam.reverse();
		break;
	// ITRF2008-to-ITRF2020
	case ITRF2020_i_xyz:
	case ITRF2020_i:
		memcpy(&tParam.parameters_, ITRF2020_ITRF2008<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2020_ITRF2008<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2020_ITRF2008<double, UINT32>::reference_frame;
		tParam.paramDirection_ = __paramReverse__;
		tParam.reverse();
		break;
	// ITRF2008-to-WGS84
	case WGS84_i_xyz:
	case WGS84_i:
	case WGS84_transit_i_xyz:
	case WGS84_transit_i:
	case WGS84_G730_i_xyz:
	case WGS84_G730_i:
	case WGS84_G873_i_xyz:
	case WGS84_G873_i:
	case WGS84_G1150_i_xyz:
	case WGS84_G1150_i:
	case WGS84_G1674_i_xyz:
	case WGS84_G1674_i:
	case WGS84_G1762_i_xyz:
	case WGS84_G1762_i:
	case WGS84_G2139_i_xyz:
	case WGS84_G2139_i:
	default:
		throw std::runtime_error("determineITRF2008Parameters: ");
	}
}
	

template <typename U>
void determineITRF2014Parameters(transformation_parameter_set& tParam)
{
	// Which 'to' datum?
	switch (tParam.from_to_.second)
	{
	//////////////////////////////////
	// ITRF2014-to-GDA2020
	case GDA2020_i_xyz:
	case GDA2020_i_2d:
	case GDA2020_i:
		memcpy(&tParam.parameters_, ITRF2014_GDA2020<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2014_GDA2020<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2014_GDA2020<double, UINT32>::reference_frame;
		break;
	// ITRF2014-to-GDA94
	case GDA94_i_xyz:
	case GDA94_i_2d:
	case GDA94_i:
		memcpy(&tParam.parameters_, ITRF2014_GDA94<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2014_GDA94<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2014_GDA94<double, UINT32>::reference_frame;
		break;
	//////////////////////////////////
	// ITRF2014-to-ITRF1988
	case ITRF1988_i_xyz:
	case ITRF1988_i:
		memcpy(&tParam.parameters_, ITRF2014_ITRF1989<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2014_ITRF1989<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2014_ITRF1989<double, UINT32>::reference_frame;
		break;
	// ITRF2014-to-ITRF1989
	case ITRF1989_i_xyz:
	case ITRF1989_i:
		memcpy(&tParam.parameters_, ITRF2014_ITRF1989<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2014_ITRF1989<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2014_ITRF1989<double, UINT32>::reference_frame;
		break;
	// ITRF2014-to-ITRF1990
	case ITRF1990_i_xyz:
	case ITRF1990_i:
		memcpy(&tParam.parameters_, ITRF2014_ITRF1990<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2014_ITRF1990<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2014_ITRF1990<double, UINT32>::reference_frame;
		break;
	// ITRF2014-to-ITRF1991
	case ITRF1991_i_xyz:
	case ITRF1991_i:
		memcpy(&tParam.parameters_, ITRF2014_ITRF1991<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2014_ITRF1991<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2014_ITRF1991<double, UINT32>::reference_frame;
		break;
	// ITRF2014-to-ITRF1992
	case ITRF1992_i_xyz:
	case ITRF1992_i:
		memcpy(&tParam.parameters_, ITRF2014_ITRF1992<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2014_ITRF1992<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2014_ITRF1992<double, UINT32>::reference_frame;
		break;
	// ITRF2014-to-ITRF1993
	case ITRF1993_i_xyz:
	case ITRF1993_i:
		memcpy(&tParam.parameters_, ITRF2014_ITRF1993<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2014_ITRF1993<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2014_ITRF1993<double, UINT32>::reference_frame;
		break;
	// ITRF2014-to-ITRF1994
	case ITRF1994_i_xyz:
	case ITRF1994_i:
		memcpy(&tParam.parameters_, ITRF2014_ITRF1994<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2014_ITRF1994<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2014_ITRF1994<double, UINT32>::reference_frame;
		break;
	// ITRF2014-to-ITRF1996
	case ITRF1996_i_xyz:
	case ITRF1996_i:
		memcpy(&tParam.parameters_, ITRF2014_ITRF1996<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2014_ITRF1996<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2014_ITRF1996<double, UINT32>::reference_frame;
		break;
	// ITRF2014-to-ITRF1997
	case ITRF1997_i_xyz:
	case ITRF1997_i:
		memcpy(&tParam.parameters_, ITRF2014_ITRF1997<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2014_ITRF1997<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2014_ITRF1997<double, UINT32>::reference_frame;
		break;
	// ITRF2014-to-ITRF2000
	case ITRF2000_i_xyz:
	case ITRF2000_i:
		memcpy(&tParam.parameters_, ITRF2014_ITRF2000<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2014_ITRF2000<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2014_ITRF2000<double, UINT32>::reference_frame;
		break;
	// ITRF2014-to-ITRF2005
	case ITRF2005_i_xyz:
	case ITRF2005_i:
		memcpy(&tParam.parameters_, ITRF2014_ITRF2005<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2014_ITRF2005<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2014_ITRF2005<double, UINT32>::reference_frame;
		break;
	// ITRF2014-to-ITRF2008
	case ITRF2008_i_xyz:
	case ITRF2008_i:
		memcpy(&tParam.parameters_, ITRF2014_ITRF2008<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2014_ITRF2008<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2014_ITRF2008<double, UINT32>::reference_frame;
		break;
	// ITRF2014-to-ITRF2014
	case ITRF2014_i_xyz:
	case ITRF2014_i:
		// not a transformation of frame, but a shift in epoch only
		// Not defined yet - requires interpolation from the ITRF2014 velocity model
		throw RefTranException("A transformation between different epochs on the same frame is required.",
			REFTRAN_TRANS_ON_PLATE_REQUIRED);
	// ITRF2014-to-ITRF2020
	case ITRF2020_i_xyz:
	case ITRF2020_i:
		memcpy(&tParam.parameters_, ITRF2020_ITRF2014<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2020_ITRF2014<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2020_ITRF2014<double, UINT32>::reference_frame;
		tParam.paramDirection_ = __paramReverse__;
		tParam.reverse();
		break;
	// ITRF2014-to-WGS84
	case WGS84_i_xyz:
	case WGS84_i:
	case WGS84_transit_i_xyz:
	case WGS84_transit_i:
	case WGS84_G730_i_xyz:
	case WGS84_G730_i:
	case WGS84_G873_i_xyz:
	case WGS84_G873_i:
	case WGS84_G1150_i_xyz:
	case WGS84_G1150_i:
	case WGS84_G1674_i_xyz:
	case WGS84_G1674_i:
	case WGS84_G1762_i_xyz:
	case WGS84_G1762_i:
	case WGS84_G2139_i_xyz:
	case WGS84_G2139_i:
	default:
		throw std::runtime_error("determineITRF2014Parameters: ");
	}
}

template <typename U>
void determineITRF2020Parameters(transformation_parameter_set& tParam)
{
	// Which 'to' datum?
	switch (tParam.from_to_.second)
	{
	//////////////////////////////////
	// ITRF2020-to-GDA2020
	case GDA2020_i_xyz:
	case GDA2020_i_2d:
	case GDA2020_i:
		memcpy(&tParam.parameters_, ITRF2020_GDA2020<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2020_GDA2020<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2020_GDA2020<double, UINT32>::reference_frame;
		break;
	// ITRF2020-to-GDA94
	case GDA94_i_xyz:
	case GDA94_i_2d:
	case GDA94_i:
		// No direct parameters exist!
		throw RefTranException(message_parameters_undefined<std::string, UINT32>(tParam.from_to_.first, tParam.from_to_.second),
			REFTRAN_DIRECT_PARAMS_UNAVAILABLE);
		break;
	//////////////////////////////////
	// ITRF2020-to-ITRF1988
	case ITRF1988_i_xyz:
	case ITRF1988_i:
		memcpy(&tParam.parameters_, ITRF2020_ITRF1988<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2020_ITRF1988<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2020_ITRF1988<double, UINT32>::reference_frame;
		break;
	// ITRF2020-to-ITRF1989
	case ITRF1989_i_xyz:
	case ITRF1989_i:
		memcpy(&tParam.parameters_, ITRF2020_ITRF1989<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2020_ITRF1989<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2020_ITRF1989<double, UINT32>::reference_frame;
		break;
	// ITRF2020-to-ITRF1990
	case ITRF1990_i_xyz:
	case ITRF1990_i:
		memcpy(&tParam.parameters_, ITRF2020_ITRF1990<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2020_ITRF1990<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2020_ITRF1990<double, UINT32>::reference_frame;
		break;
	// ITRF2020-to-ITRF1991
	case ITRF1991_i_xyz:
	case ITRF1991_i:
		memcpy(&tParam.parameters_, ITRF2020_ITRF1991<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2020_ITRF1991<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2020_ITRF1991<double, UINT32>::reference_frame;
		break;
	// ITRF2020-to-ITRF1992
	case ITRF1992_i_xyz:
	case ITRF1992_i:
		memcpy(&tParam.parameters_, ITRF2020_ITRF1992<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2020_ITRF1992<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2020_ITRF1992<double, UINT32>::reference_frame;
		break;
	// ITRF2020-to-ITRF1993
	case ITRF1993_i_xyz:
	case ITRF1993_i:
		memcpy(&tParam.parameters_, ITRF2020_ITRF1993<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2020_ITRF1993<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2020_ITRF1993<double, UINT32>::reference_frame;
		break;
	// ITRF2020-to-ITRF1994
	case ITRF1994_i_xyz:
	case ITRF1994_i:
		memcpy(&tParam.parameters_, ITRF2020_ITRF1994<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2020_ITRF1994<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2020_ITRF1994<double, UINT32>::reference_frame;
		break;
	// ITRF2020-to-ITRF1996
	case ITRF1996_i_xyz:
	case ITRF1996_i:
		memcpy(&tParam.parameters_, ITRF2020_ITRF1996<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2020_ITRF1996<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2020_ITRF1996<double, UINT32>::reference_frame;
		break;
	// ITRF2020-to-ITRF1997
	case ITRF1997_i_xyz:
	case ITRF1997_i:
		memcpy(&tParam.parameters_, ITRF2020_ITRF1997<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2020_ITRF1997<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2020_ITRF1997<double, UINT32>::reference_frame;
		break;
	// ITRF2020-to-ITRF2000
	case ITRF2000_i_xyz:
	case ITRF2000_i:
		memcpy(&tParam.parameters_, ITRF2020_ITRF2000<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2020_ITRF2000<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2020_ITRF2000<double, UINT32>::reference_frame;
		break;
	// ITRF2020-to-ITRF2005
	case ITRF2005_i_xyz:
	case ITRF2005_i:
		memcpy(&tParam.parameters_, ITRF2020_ITRF2005<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2020_ITRF2005<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2020_ITRF2005<double, UINT32>::reference_frame;
		break;
	// ITRF2020-to-ITRF2008
	case ITRF2008_i_xyz:
	case ITRF2008_i:
		memcpy(&tParam.parameters_, ITRF2020_ITRF2008<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2020_ITRF2008<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2020_ITRF2008<double, UINT32>::reference_frame;
		break;
	// ITRF2020-to-ITRF2014
	case ITRF2014_i_xyz:
	case ITRF2014_i:
		memcpy(&tParam.parameters_, ITRF2020_ITRF2014<double, UINT32>::transformationParameters, sizeof(double) * 14);
		tParam.reference_epoch_ = ITRF2020_ITRF2014<double, UINT32>::reference_epoch;
		tParam.reference_frame_ = ITRF2020_ITRF2014<double, UINT32>::reference_frame;
		break;
	// ITRF2020-to-ITRF2020
	case ITRF2020_i_xyz:
	case ITRF2020_i:
		// not a transformation of frame, but a shift in epoch only
		// Not defined yet - requires interpolation from the ITRF2014 velocity model
		throw RefTranException("A transformation between different epochs on the same frame is required.",
			REFTRAN_TRANS_ON_PLATE_REQUIRED);
	// ITRF2014-to-WGS84
	case WGS84_i_xyz:
	case WGS84_i:
	case WGS84_transit_i_xyz:
	case WGS84_transit_i:
	case WGS84_G730_i_xyz:
	case WGS84_G730_i:
	case WGS84_G873_i_xyz:
	case WGS84_G873_i:
	case WGS84_G1150_i_xyz:
	case WGS84_G1150_i:
	case WGS84_G1674_i_xyz:
	case WGS84_G1674_i:
	case WGS84_G1762_i_xyz:
	case WGS84_G1762_i:
	case WGS84_G2139_i_xyz:
	case WGS84_G2139_i:
	default:
		throw std::runtime_error("determineITRF2020Parameters: ");
	}
}

// This function returns the Australian plate motion model parameters, which are selected when the user
// chooses to use the Australian PMM as default for all transformations.  
template <typename T, typename U>
void getAustralianPlateMotionModelParameters(transformation_parameter_set& tParam)
{
	tParam.paramDirection_ = __paramForward__;
	memcpy(&tParam.parameters_, AUS_PLATE_MOTION_MODEL<T, U>::transformationParameters, sizeof(T) * 14);
	tParam.reference_epoch_ = AUS_PLATE_MOTION_MODEL<T, U>::reference_epoch;
	tParam.reference_frame_ = AUS_PLATE_MOTION_MODEL<T, U>::reference_frame;
}

// This function returns plate motion model parameters according to pre-defined rotations defined
// from Euler rotations.  
template <typename T, typename U>
void setDefinedPlateMotionModelParameters(transformation_parameter_set& tParam, const T& x, const T& y, const T& z)
{
	tParam.paramDirection_ = __paramForward__;

	tParam.parameters_[0] = 0.0;	// x translation (millimetres)
	tParam.parameters_[1] = 0.0;	// y translation (millimetres)
	tParam.parameters_[2] = 0.0;	// z translation (millimetres)
	tParam.parameters_[3] = 0.0;	// scale (ppb)
	tParam.parameters_[4] = 0.0;	// x rotation (milli-arc-seconds)
	tParam.parameters_[5] = 0.0;	// y rotation (milli-arc-seconds)
	tParam.parameters_[6] = 0.0;	// z rotation (milli-arc-seconds)
	tParam.parameters_[7] = 0.0;	// x translation rate (millimetres p/yr)
	tParam.parameters_[8] = 0.0;	// y translation rate (millimetres p/yr)
	tParam.parameters_[9] = 0.0;	// z translation rate (millimetres p/yr)
	tParam.parameters_[10] = 0.0;	// scale rate (ppb p/yr)
	tParam.parameters_[11] = x;		// x rotation rate (milli-arc-seconds p/yr)
	tParam.parameters_[12] = y;		// y rotation rate (milli-arc-seconds p/yr)
	tParam.parameters_[13] = z;		// z rotation rate (milli-arc-seconds p/yr)
}

template <typename U>
void determineHelmertParameters(transformation_parameter_set& tParam)
{
	tParam.paramDirection_ = __paramForward__;

	try 
	{
		// Which 'from' datum?
		switch (tParam.from_to_.first)
		{
		// GDA94
		case GDA94_i_xyz:
		case GDA94_i_2d:
		case GDA94_i:
			determineGDA94Parameters<U>(tParam);
			break;
		// GDA2020
		case GDA2020_i_xyz:
		case GDA2020_i_2d:
		case GDA2020_i:
			determineGDA2020Parameters<U>(tParam);
			break;
		// ITRF1988
		case ITRF1988_i_xyz:
		case ITRF1988_i:
			determineITRF1988Parameters<U>(tParam);
			break;
		// ITRF1989
		case ITRF1989_i_xyz:
		case ITRF1989_i:
			determineITRF1989Parameters<U>(tParam);
			break;
		// ITRF1990
		case ITRF1990_i_xyz:
		case ITRF1990_i:
			determineITRF1990Parameters<U>(tParam);
			break;
		// ITRF1991
		case ITRF1991_i_xyz:
		case ITRF1991_i:
			determineITRF1991Parameters<U>(tParam);
			break;
		// ITRF1992
		case ITRF1992_i_xyz:
		case ITRF1992_i:
			determineITRF1992Parameters<U>(tParam);
			break;
		// ITRF1993
		case ITRF1993_i_xyz:
		case ITRF1993_i:
			determineITRF1993Parameters<U>(tParam);
			break;
		// ITRF1994
		case ITRF1994_i_xyz:
		case ITRF1994_i:
			determineITRF1994Parameters<U>(tParam);
			break;
		// ITRF1996
		case ITRF1996_i_xyz:
		case ITRF1996_i:
			determineITRF1996Parameters<U>(tParam);
			break;
		// ITRF1997
		case ITRF1997_i_xyz:
		case ITRF1997_i:
			determineITRF1997Parameters<U>(tParam);
			break;
		// ITRF2000
		case ITRF2000_i_xyz:
		case ITRF2000_i:
			determineITRF2000Parameters<U>(tParam);
			break;
		// ITRF2005
		case ITRF2005_i_xyz:
		case ITRF2005_i:
			determineITRF2005Parameters<U>(tParam);
			break;
		// ITRF2008
		case ITRF2008_i_xyz:
		case ITRF2008_i:
			determineITRF2008Parameters<U>(tParam);
			break;
		// ITRF2014
		case ITRF2014_i_xyz:
		case ITRF2014_i:
			determineITRF2014Parameters<U>(tParam);
			break;
		// ITRF2020
		case ITRF2020_i_xyz:
		case ITRF2020_i:
			determineITRF2020Parameters<U>(tParam);
			break;
		// WGS84
		case WGS84_i_xyz:
		case WGS84_i:
		case WGS84_transit_i_xyz:
		case WGS84_transit_i:
		case WGS84_G730_i_xyz:
		case WGS84_G730_i:
		case WGS84_G873_i_xyz:
		case WGS84_G873_i:
		case WGS84_G1150_i_xyz:
		case WGS84_G1150_i:
		case WGS84_G1674_i_xyz:
		case WGS84_G1674_i:
		case WGS84_G1762_i_xyz:
		case WGS84_G1762_i:
		case WGS84_G2139_i_xyz:
		case WGS84_G2139_i:
			throw RefTranException("Transformation parameters from WGS 84 are unavailable.",
				REFTRAN_WGS84_TRANS_UNSUPPORTED);
			break;
		default:
			throw std::runtime_error("");
		}
	}
	catch (std::runtime_error& e)
	{
		std::stringstream ss;
		ss << "determineHelmertParameters():" << std::endl << 
			"  " << e.what() << " Parameters for transforming between" << std::endl <<
			"  " << datumFromEpsgCode<std::string, UINT32>(tParam.from_to_.first) << " and " << 
			datumFromEpsgCode<std::string, UINT32>(tParam.from_to_.second) <<
			" have not been defined yet." << std::endl;
		throw std::runtime_error(ss.str());
	}
}

}	// namespace datum_parameters
}	// namespace dynadjust

#endif /* DNATRANSFORM_PARAM_H_ */
