//============================================================================
// Name         : dnadatumprojectionparam.hpp
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
// Description  : DynAdjust Datum and Projection parameters
//============================================================================

#ifndef DNADATUMPROJECTIONPARAM_H_
#define DNADATUMPROJECTIONPARAM_H_

#if defined(_MSC_VER)
	#if defined(LIST_INCLUDES_ON_BUILD) 
		#pragma message("  " __FILE__) 
	#endif
#endif

/// \cond
#include <math.h>
/// \endcond

#include <include/config/dnatypes-fwd.hpp>

// GRS80 parameters
const double GRS80_a = 6378137.0;				// Semi major axis (a)
const double GRS80_inv_f = 298.257222101;		// Inverse flattening (1/f)

// WGS84 parameters
const double WGS84_a = 6378137.0;				// Semi major axis (a)
const double WGS84_inv_f = 298.25722360;		// Inverse flattening (1/f)

// ANS parameters
const double ANS_a = 6378160.0;					// Semi major axis (a)
const double ANS_inv_f = 298.25;				// Inverse flattening (1/f)

// International 1924 parameters
const double International24_a = 6378388.0;     // Semi major axis (a)
const double International24_inv_f = 297.0;		// Inverse flattening (1/f)

// UTM parameters
const double FALSE_E = 500000.0;				// False Easting
const double FALSE_N = 10000000.0;				// False Northing
const double K0 = 0.9996;						// Central scale factor
const double ZW = 6.0;							// Zone width
const double LCMZ1 = -177.0;					// Longitude of central meridian of zone 1
const double LWEZ0 = -186.0;					// Longitude of western edge of zone 0
const double LCMZ0 = -183.0;					// Longitude of central meridian of zone 0

// nu
template <class T>
T primeVertical_(const T& semiMajor, const T& e1Sqd, const T& latitude) 
{
	return (semiMajor / sqrt(1.0 - e1Sqd * (sin(latitude)*sin(latitude))));
}

// rho
template <class T>
double primeMeridian_(const T& semiMajor, const T& e1Sqd, const T& latitude)
{
	T dDel = sqrt(1.0 - e1Sqd * (sin(latitude)*sin(latitude)));
	return (semiMajor * ((1.0 - e1Sqd) / (dDel * dDel * dDel)) );
}

// nu and rho
template <class T>
void primeVerticalandMeridian_(const T& semiMajor, const T& e1Sqd, const T& latitude, T& nu, T& rho) 
{
	T dDel = sqrt(1.0 - e1Sqd * (sin(latitude)*sin(latitude)));

	nu = semiMajor / dDel;
	rho = semiMajor * ((1.0 - e1Sqd) / (dDel * dDel * dDel));
}

// average radius of curvature
template <class T>
T averageRadiusofCurvature_(const T& semiMajor, const T& e1Sqd, const T& latitude)
{
	T nu, rho;
	primeVerticalandMeridian_(semiMajor, e1Sqd, latitude, nu, rho);
	return sqrt(nu * rho);
}

typedef enum _DATUM_TYPE_
{
	STATIC_DATUM = 0,
	DYNAMIC_DATUM = 1
} DATUM_TYPE;

typedef enum _ELLIPSOID_TYPE_
{
	EPSG_DEFINED_ELLIPSOID = 0,
	USER_DEFINED_ELLIPSOID = 1
} ELLIPSOID_TYPE;

typedef enum _PROJECTION_NAME_
{
	UTM = 0,
	USER_DEFINED_PROJECTION = 2
} PROJECTION_NAME;

namespace dynadjust {
namespace datum_parameters {



}	// namespace datum_parameters
}	// namespace dynadjust

#endif /* DNADATUMPROJECTIONPARAM_H_ */
