//============================================================================
// Name         : dnatemplatecalcfuncs.hpp
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
// Description  : Advanced and common calculation functions using standard
//				  data types
//============================================================================

#ifndef DNATEMPLATECALCFUNCS_H_
#define DNATEMPLATECALCFUNCS_H_

#if defined(_MSC_VER)
	#if defined(LIST_INCLUDES_ON_BUILD) 
		#pragma message("  " __FILE__) 
	#endif
#endif

/// \cond
#include <numeric>
#include <math.h>
#include <memory>
/// \endcond

#include <include/config/dnatypes.hpp>
#include <include/config/dnaconsts.hpp>

template <class U>
std::size_t sumOfConsecutiveIntegers(const U& max, const U& min = 1)
{
    std::size_t sum(0);
	for (U u=min; u<=max; u++)
		sum += u;
	return sum;	
}

template <typename T, typename U>
T removeNegativeZero(const T& t, const U& precision)
{
	if (t < 0.0 || t == T(-0.0))
	{
		T v(floor((t * pow(10.0, precision)) + 0.5));
		if (fabs(v) > 0.)
			// t is non-zero at precision
			return t;
		// t is a negative zero
		return 0.;
	}
	return t;	
}

template <typename T>
bool are_floats_equal(const T& lhs, const T& rhs, const T epsilon=PRECISION_1E6)
{
	return (fabs(lhs - rhs) < epsilon);
}

template <typename T, typename U, typename iterator>
T average(const iterator begin, const iterator end, U& sum)
{
	U n(static_cast<UINT32>(std::distance(begin, end)));
	sum = std::accumulate(begin, end, 0);
	return static_cast<T>(sum) / n;
}

template <class T, class U>
T average(const U& a, const U& b)
{
	T t(a + b);
	return t / 2.;
}

template <class T>
T average(const T& a, const T& b)
{
	T t(a + b);
	return (t / (T)(2.));
}

template <class T, typename U>
T minVal(const T& lhs, const U& rhs)
{
	if (lhs < rhs)
		return lhs;
	else 
		return rhs;
}

template <class T, typename U>
T maxVal(const T& lhs, const U& rhs)
{
	if (lhs > rhs)
		return lhs;
	else 
		return rhs;
}

// use this for longitude values only, where 180 degrees
// marks the boundary between east (positive) and west (negative)
template <class T>
T Rad180Mod(const T &dValue)
{
	if (dValue < -PI)
		 return dValue + (PI + PI);
	else if (dValue > PI)
		 return dValue - (PI + PI);
	return dValue;
}

template <class T>
T degrees_to_radians_(T& degrees)
{
	return (degrees * DEG_TO_RAD);
}

template <class T>
T Radians(T& degrees)
{
	return degrees_to_radians_(degrees);
}

template <class T>
T Radians(const T& degrees)
{
	return degrees_to_radians_(degrees);
}

template <class T>
void Radians(T* degrees)
{
	*degrees = degrees_to_radians_(*degrees);
}

template <class T>
T radians_to_degrees_(T& radians)
{
	return (radians * RAD_TO_DEG);
}

template <class T>
T Degrees(T& radians)
{
	return radians_to_degrees_(radians);
}

template <class T>
T Degrees(const T& radians)
{
	return radians_to_degrees_(radians);
}

template <class T>
void Degrees(T* radians)
{
	*radians = radians_to_degrees_(*radians);
}

// use DegreesL(const T& radians) for longitude values only, where 180 degrees
// marks the boundary between east (positive) and west (negative)
template <class T>
T DegreesL(const T& radians)
{
	return (Rad180Mod(radians) * RAD_TO_DEG);
}

// use DegreesL(T* dValue) for longitude values only, where 180 degrees
// marks the boundary between east (positive) and west (negative)
template <class T>
void DegreesL(T* dValue)
{
	*dValue = (Rad180Mod(*dValue) * RAD_TO_DEG);
}

// Seconds from Radians
template <class T>
T Seconds(const T& radians)
{
	return radians * RAD_TO_SEC;
}

template <class T>
T Seconds(T* radians)
{
	return *radians * RAD_TO_SEC;
}

template <class T>
T SecondstoRadians(const T& seconds)
{
	return seconds / static_cast<T>(RAD_TO_SEC);
}

template <class T, typename U>
void DegtoDms(const T& dDegrees, U* dDegMinSec)
{
	T d, m, s;

	*dDegMinSec = fabs(dDegrees);		 // retain original value
	d = floor(*dDegMinSec);
	m = floor((((*dDegMinSec) - d) * 60.0));
	s = ((*dDegMinSec) - d - (m/60.0)) * 3600.0;
	if (fabs(s - 60.0) < 0.000000001)
	{
		s = 0.0;
		m += 1.0;
	}
	*dDegMinSec = d + (m/100.0) + (s/10000.0);
	if (dDegrees < 0.0)
		*dDegMinSec *= -1;
}


// DegtoDms helper
template <class T>
T DegtoDms(const T& dDegrees)
{
	T dDegMinSec;
	DegtoDms(dDegrees, &dDegMinSec);
	return dDegMinSec;
}


template <class T, typename U>
void DmstoDeg(const T& dDegMinSec, U* dDegrees)
{
	T dh, dm, ds;

	dh = fabs(dDegMinSec);		 // retain original value
	*dDegrees = floor(dh);
	dm = floor(((dh - (*dDegrees)) * 100.0) + 0.0001);
	ds = (((dh - (*dDegrees)) * 100.0) - dm) * 100.0;
	*dDegrees += ((dm / 60.0) + (ds / 3600.0));
	if (dDegMinSec < 0.0)
		(*dDegrees) *= -1;
}

// DmstoDeg helper
template <class T>
T DmstoDeg(const T& dDegMinSec)
{
	T dDegrees;
	DmstoDeg(dDegMinSec, &dDegrees);
	return dDegrees;
}

template <class T, typename U>
void DmstoRad(const T& dDegMinSec, U* dRadians)
{
	T dh, dm, ds;

	dh = fabs(dDegMinSec);		 // retain original value
	*dRadians = floor(dh);
	dm = floor(((dh - (*dRadians)) * 100.0) + 0.0001);
	ds = (((dh - (*dRadians)) * 100.0) - dm) * 100.0;
	*dRadians += ((dm / 60.0) + (ds / 3600.0));
	if (dDegMinSec < 0.0)
		(*dRadians) *= -1;
	Radians(dRadians);
}

// DmstoRad helper
template <class T>
T DmstoRad(const T& dDegMinSec)
{
	T dRadians;
	DmstoRad(dDegMinSec, &dRadians);
	return dRadians;
}

template <class T>
T RadtoDms(const T& dRadians)
{
	T dDms;
	DegtoDms(Degrees(dRadians), &dDms);
	return dDms;
}

template <class T>
T RadtoDmsL(const T& dRadians)
{
	T dDms;
	DegtoDms(DegreesL(dRadians), &dDms);
	return dDms;
}

template <class T, typename U>
void DmintoDeg(const T& dDegMin, U* dDegrees)
{
	T dh;

	dh = fabs(dDegMin);	// retain original value
	*dDegrees = floor(dh);
	*dDegrees += (dh - *dDegrees) * 100.0 / 60.0;
	if (dDegMin < 0.0)
		*dDegrees *= -1;
}

// DmintoDeg helper
template <class T>
T DmintoDeg(const T& dDegMin)
{
	T dDegrees;
	DmintoDeg(dDegMin, &dDegrees);
	return dDegrees;
}

template <class T, typename U>
void DegtoDmin(const T& dDegrees, U* dDegMin)
{
	T d;
	
	*dDegMin = fabs(dDegrees);	// retain original value
	d = floor(*dDegMin);
	*dDegMin = d + ((*dDegMin - d) * 0.6);
	
	if (dDegrees < 0.0)
		*dDegMin *= -1;
}

// DegDmin helper
template <class T>
T DegtoDmin(const T& dDegrees)
{
	T dDegMin;
	DegtoDin(dDegrees, &dDegMin);
	return dDegMin;
}


template <class T>
//            |
//     4th    |   1st
//            |
//	----------------------
//            |
//     3rd    |   2nd
//            |
T atan_2(const T& x, const T& y)
{
	T theta(atan(x / y));			// first quadrant (default)
	if (y < 0)						// second or third quadrant
		return theta + PI;
	else
	{
		if (x > 0)					// first quadrant 
			return theta;
		else						// fourth quadrant
			return theta + TWO_PI;
	}
  return 0.;
}

#endif /* DNATEMPLATECALCFUNCS_H_ */
