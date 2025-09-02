//============================================================================
// Name         : dnatemplategeodesyfuncs.hpp
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
// Description  : Basic Geodetic Functions
//============================================================================

#ifndef DNATEMPLATEGEODESYFUNCS_H_
#define DNATEMPLATEGEODESYFUNCS_H_

#if defined(_MSC_VER)
	#if defined(LIST_INCLUDES_ON_BUILD) 
		#pragma message("  " __FILE__) 
	#endif
#endif

/// \cond
#include <algorithm>
#include <functional>
#include <sstream>
#include <string>
#include <vector>
#include <math.h>
#include <iostream>
#include <memory>
/// \endcond

#include <include/parameters/dnaellipsoid.hpp>
#include <include/parameters/dnaprojection.hpp>
#include <include/config/dnatypes.hpp>

using namespace dynadjust::datum_parameters;

// nu helper
template <class T>
T primeVertical(const CDnaEllipsoid* ellipsoid, const T& latitude) 
{
	return primeVertical_(ellipsoid->GetSemiMajor(), ellipsoid->GetE1sqd(), latitude);
}

// rho helper
template <class T>
double primeMeridian(const CDnaEllipsoid* ellipsoid, const T& latitude)
{
	return primeMeridian_(ellipsoid->GetSemiMajor(), ellipsoid->GetE1sqd(), latitude);
}

// nu and rho helper
template <class T>
void primeVerticalandMeridian(const CDnaEllipsoid* ellipsoid, const T& latitude, T& nu, T& rho) 
{
	primeVerticalandMeridian_(ellipsoid->GetSemiMajor(), ellipsoid->GetE1sqd(), latitude, nu, rho);
}

// average radius of curvature helper
template <class T>
T averageRadiusofCurvature(const CDnaEllipsoid* ellipsoid, const T& latitude)
{
	return averageRadiusofCurvature_(ellipsoid->GetSemiMajor(), ellipsoid->GetE1sqd(), latitude);
}


template <class T>
void GeoToCart(const T& Latitude, const T& Longitude, const T& Height, T* X, T* Y, T* Z, 
			   const CDnaEllipsoid* ellipsoid)
{
	// copy variables in case the caller overwrites original values
	T latitude(Latitude), longitude(Longitude), height(Height);
	
	// calculate prime vertical (once) 
	T Nu(primeVertical(ellipsoid, latitude));
			
	*X = (Nu + height) * cos(latitude) * cos(longitude);
	*Y = (Nu + height) * cos(latitude) * sin(longitude);
	*Z = ((Nu * (1. - ellipsoid->GetE1sqd())) + height) * sin(latitude);
}


template <class T>
void CartToGeo_SimpleIteration(const T& X, const T& Y, const T& Z,
			   T* latitude, T* longitude, T* height,
			   const CDnaEllipsoid* ellipsoid) // ellipsoid parameters
{
	// copy variables in case the caller overwrites original values
	T x(X), y(Y), z(Z);

	T f(ellipsoid->GetFlattening());
	T e2(ellipsoid->GetE1sqd());
	T p(pow(((x * x) + (y * y)), 0.5));

	// "Cartesian to geographic" conversion problem resides in the 
	// equation to compute latitude (phi)...latitude is required on both
	// sides of the equation. A simple iterative method is used to determine 
	// the latitude, whereby the starting latitude is computed from 
	// tan(phi) = z / p

	// Compute starting lat value from second eccentricity squared
	T lat1(atan(z / p));
	*latitude = atan2((z + (e2 * ellipsoid->GetSemiMajor() * sin(lat1))), p);
	T Nu(primeVertical(ellipsoid, *latitude));

	for (UINT16 i(0); i<16; i++)
	{
		lat1 = atan2((z + (e2 * Nu * sin(*latitude))), p);
		//Nu = ellipsoid->PrimeVertical(*latitude);
		Nu = primeVertical(ellipsoid, *latitude);
		*latitude = atan2((z + (e2 * Nu * sin(lat1))), p);

		if (fabs(lat1 - *latitude) < PRECISION_1E16)
			break;
	}
	
	// Compute longitude
	*longitude = atan(y / x);

	// determine correct quadrant and apply negative long accordingly
	if (x < 0.0 && y > 0.0)
		*longitude += PI;
	else if (x < 0.0 && y < 0.0)
		*longitude = -(PI - *longitude);
	
	// Compute height
	*height = (p / cos(*latitude)) - Nu;
}



template <class T>
// K Lin and J Wang's method based on Newton's iteration
// 	double dXAxis(-3563081.362), dYAxis(-2057145.984), dZAxis(-4870449.482), dHeight(0.);
//	CDnaEllipsoid e;
//	CartToGeo<double>(dXAxis, dYAxis, dZAxis, &dXAxis, &dYAxis, &dHeight, &e);
//	std::stringstream ss;
//	ss << std::setw(MSR) << std::right << FormatDmsString(RadtoDms(dXAxis), 5, true, false) << ", ";
//	ss << std::setw(MSR) << std::right << FormatDmsString(RadtoDms(dYAxis), 5, true, false) << ", ";
//	ss << std::setw(MSR) << std::setprecision(4) << std::fixed << std::right << dHeight;
//	std::string comp(ss.str());
//	comp should equal  "-50 00 00.0000, -150 00 00.0000, 10000.000"
//
void CartToGeo(const T& X, const T& Y, const T& Z,
			   T* latitude, T* longitude, T* height,
			   const CDnaEllipsoid* ellipsoid) // ellipsoid parameters
{
	// copy variables in case the caller overwrites original values
	T x(X), y(Y), z(Z);

	T p2((x * x) + (y * y));
	T p(sqrt(p2));
	T a2(ellipsoid->GetSemiMajor() * ellipsoid->GetSemiMajor());
	T b2(ellipsoid->GetSemiMinor() * ellipsoid->GetSemiMinor());
	T Z2(z * z);
	T a2Z2(a2 * Z2);
	T b2p2(b2 * p2);
	T A(a2Z2 + b2p2);

	// Compute initial approximation of m (Lin and Wang 1995, eq. 9, p. 301)
	T m0((ellipsoid->GetSemiMajor() * ellipsoid->GetSemiMinor() * sqrt(A) * A - a2 * b2 * A) / (2. *
		((a2 * a2Z2) + (b2 * b2p2))));

	T twom, a2twom, b2twom, f, df, m;
	
	// Generally converges after one iteration and 
	// so shouldn't need more than five iterations.  
	for (UINT16 i(0); i<5; ++i)
	{
		m = m0;
		twom = m * 2.;
		a2twom = a2 + twom;
		b2twom = b2 + twom;
		
		f = (a2 * p2 / (a2twom * a2twom)) + (b2 * Z2 / (b2twom * b2twom)) - 1.;
		
		// if f is sufficiently close to zero, break.
		if (fabs(f) < PRECISION_1E12)
			break;
		
		df = -4. * ((a2 * p2 / (a2twom * a2twom * a2twom)) + (b2 * Z2 / (b2twom * b2twom * b2twom)));
		
		// recompute new value for m
		m0 = m - (f / df);
		m = m0;	
	}

	twom = m * 2.;

	T p_E(a2 * p / (a2 + twom));
	T Z_E(b2 * z / (b2 + twom));

	// Compute latitude
	*latitude = atan(a2 * Z_E / (b2 * p_E));

	// Compute longitude
	*longitude = atan(y / x);

	// determine correct quadrant and apply negative long accordingly
	if (x < 0.0 && y > 0.0)
		*longitude += PI;
	else if (x < 0.0 && y < 0.0)
		*longitude = -(PI - *longitude);

	// The following line causes an issue for west longitudes, which by nature are negative!
	// Not sure why this was introduced.  Removing the conditional absolute has no adverse impact on
	// positions which are located in the eastern hemisphere,
	//if (*longitude < 0.)
	//	*longitude += TWO_PI;
	
	// Compute height
	*height = sqrt(((p - p_E) * (p - p_E)) + ((z - Z_E) * (z - Z_E)));
	if ((p + fabs(z)) < (p_E + fabs(Z_E)))
		*height *= -1.;
}


template <class T>
// K Lin and J Wang's method based on Newton's iteration
T CartToLat(const T& X, const T& Y, const T& Z, const CDnaEllipsoid* ellipsoid)
{
	// copy variables in case the caller overwrites original values
	T x(X), y(Y), z(Z);

	T p2((x * x) + (y * y));
	T p(sqrt(p2));
	T a2(ellipsoid->GetSemiMajor() * ellipsoid->GetSemiMajor());
	T b2(ellipsoid->GetSemiMinor() * ellipsoid->GetSemiMinor());
	T Z2(z * z);
	T a2Z2(a2 * Z2);
	T b2p2(b2 * p2);
	T A(a2Z2 + b2p2);

	// Compute initial approximation of m (Lin and Wang 1995, eq. 9, p. 301)
	T m0((ellipsoid->GetSemiMajor() * ellipsoid->GetSemiMinor() * sqrt(A) * A - a2 * b2 * A) / (2. *
		((a2 * a2Z2) + (b2 * b2p2))));

	T twom, a2twom, b2twom, f, df, m;
	
	// Generally converges after one iteration and 
	// so shouldn't need more than five iterations.  
	for (UINT16 i(0); i<5; ++i)
	{
		m = m0;
		twom = m * 2.;
		a2twom = a2 + twom;
		b2twom = b2 + twom;
		
		f = (a2 * p2 / (a2twom * a2twom)) + (b2 * Z2 / (b2twom * b2twom)) - 1.;
		
		// if f is sufficiently close to zero, break.
		if (fabs(f) < PRECISION_1E12)
			break;
		
		df = -4. * ((a2 * p2 / (a2twom * a2twom * a2twom)) + (b2 * Z2 / (b2twom * b2twom * b2twom)));
		
		// recompute new value for m
		m0 = m - (f / df);
		m = m0;	
	}

	twom = m * 2.;

	T p_E(a2 * p / (a2 + twom));
	T Z_E(b2 * z / (b2 + twom));

	// Compute latitude
	return atan(a2 * Z_E / (b2 * p_E));
}

template <class T>
T PartialD_Latitude(const T& X, const T& Y, const T& Z,
			   const _CART_ELEM_& element,  const T& latitude, const CDnaEllipsoid* ellipsoid)
{
	if (element > z_element)
		return 0.;

	const T small_inc = PRECISION_1E4;

	// Compute the partial derivative.
	// 1. add small increment to the required element
	T cart[3] = { X, Y, Z };
	cart[element] += small_inc;

	// 2. f(x + small_inc)
	T fx_small_inc(CartToLat(
		cart[x_element],				// X1
		cart[y_element],				// Y1
		cart[z_element],				// Z1
		ellipsoid));

	//			  f(x + small_inc) - f(x) 
	// 3. f'(x) = -----------------------
	//					 small_inc
	return (fx_small_inc - latitude) / small_inc;
}


template <class T>
T PartialD_Latitude_F(const T& X, const T& Y, const T& Z,
			   const _CART_ELEM_& element,  T* latitude, const CDnaEllipsoid* ellipsoid)
{
	if (element > z_element)
		return 0.;

	// compute the new latitude
	*latitude = CartToLat(X, Y, Z, ellipsoid);

	return PartialD_Latitude(X, Y, Z, element, *latitude, ellipsoid);
}

template <class T>
T PartialD_HorizAngle(const T X1, const T Y1, const T Z1,
				 const T X2, const T Y2, const T Z2, 
				 const T X3, const T Y3, const T Z3, 
				 const T currentLatitude, const T currentLongitude,
				 const _STATION_ELEM_& station, const _CART_ELEM_& element, const T angle)
{
	if (element > z_element)
		return 0.;

	const T small_inc = PRECISION_1E4;

	// Compute the partial derivative.
	// 1. add small increment to the required element
	T cart[3][3] = { X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3 };
	cart[station][element] += small_inc;

	// Temporary variables
	T dir12, dir13, loc12e, loc12n, loc13e, loc13n;

	// 2. f(x + small_inc)
	T fx_small_inc(HorizontalAngle(
		cart[station_1][x_element],				// X1
		cart[station_1][y_element],				// Y1
		cart[station_1][z_element],				// Z1
		cart[station_2][x_element],				// X2
		cart[station_2][y_element],				// Y2
		cart[station_2][z_element],				// Z2
		cart[station_3][x_element],				// X3
		cart[station_3][y_element],				// Y3ZONE
		cart[station_3][z_element],				// Z3
		currentLatitude, currentLongitude,
		&dir12, &dir13, &loc12e, &loc12n, &loc13e, &loc13n));

	//			  f(x + small_inc) - f(x) 
	// 3. f'(x) = -----------------------
	//					 small_inc
	return (fx_small_inc - angle) / small_inc;
}


template <class T>
// latitude and longitude in radians
void GeoToGrid(const T& Latitude, const T& Longitude, T* easting, T* northing, T* zone, 
			   const CDnaEllipsoid* ellipsoid, const CDnaProjection* projection, bool COMPUTE_ZONE)
{
	T latitude(Latitude);
	T longitude(Longitude);

	// Compute Zone if not previously specified
	if (COMPUTE_ZONE)
		*zone = floor((Degrees(longitude) - projection->GetLongWesternEdgeZone0()) / projection->GetZoneWidth());

	// Compute Geodetic Longitude of the Central Meridian of pGeoValue->dNum3 (the UTM Zone)
	T CMeridian((*zone * projection->GetZoneWidth()) + projection->GetLongCentralMeridianZone0());

	// Compute diff in Longitude between CMeridian and Longitude
	T omega(longitude - Radians(CMeridian));

	T e2(ellipsoid->GetE1sqd());
	T e2_2(pow(e2, 2.0));
	T e2_3(pow(e2, 3.0));

	T Nu, Rho, Nu_div_Rho;
	primeVerticalandMeridian(ellipsoid, latitude, Nu, Rho);
	Nu_div_Rho = Nu / Rho;

	T A0(1.0 - (e2 / 4.0) - (3.0 * e2_2 / 64.0) - (5.0 * e2_3 / 256.0));
	T A2(3.0 / 8.0 * (e2 + (e2_2 / 4.0) + (15.0 * e2_3 / 128.0)));
	T A4(15.0 / 256.0 * (e2_2 + (3.0 * e2_3 / 4.0)));
	T A6(35.0 * e2_3 / 3072.0);

	T m(ellipsoid->GetSemiMajor() * ((A0 * latitude) - 
		(A2 * sin(2.0 * latitude)) +
		(A4 * sin(4.0 * latitude)) -
		(A6 * sin(6.0 * latitude))));

	T cos_lat(cos(latitude));
	T sin_lat(sin(latitude));
	T Num1(K0 * Nu * omega * cos_lat);

	T tan_lat_2(pow(tan(latitude), 2.0));
	T tan_lat_4(pow(tan(latitude), 4.0));

	// Compute Easting
	T Term1((pow(omega, 2.0) / 6.0) * (pow(cos_lat, 2.0)) * (Nu_div_Rho - tan_lat_2));
	T Term2((pow(omega, 4.0) / 120) * (pow(cos_lat, 4.0)) *
		(((4.0 * pow(Nu_div_Rho, 3.0)) * (1.0 - (6.0 * tan_lat_2))) +
		((pow(Nu_div_Rho, 2.0)) * (1.0 + (8.0 * tan_lat_2))) -
		(Nu_div_Rho * 2.0 * tan_lat_2) + tan_lat_4));
	T Term3((pow(omega, 6.0) / 5040) * (pow(cos_lat, 6.0)) *
		(61.0 - (479.0 * tan_lat_2) + (179.0 * tan_lat_4) - (tan_lat_4)));

	*easting = (Num1 * (1.0 + Term1 + Term2 + Term3)) + FALSE_E;

	// Compute Northing
	Term1 = (pow(omega, 2.0) / 2.0 * Nu * sin_lat * cos_lat);
	Term2 = (pow(omega, 4.0) / 24.0 * Nu * sin_lat * pow(cos_lat, 3.0));
	Term2 *= ((4.0 * pow(Nu_div_Rho, 2.0)) + Nu_div_Rho - tan_lat_2);
	Term3 = (pow(omega, 6.0) / 720.0 * Nu * sin_lat * pow(cos_lat, 5.0));
	Term3 *= ((8.0 * pow(Nu_div_Rho, 4.0) * (11.0 - (24.0 * tan_lat_2))) - 
		(28.0 * pow(Nu_div_Rho, 3.0) * (1.0 - (6.0 * tan_lat_2))) +
		(pow(Nu_div_Rho, 2.0) * (1.0 - (32.0 * tan_lat_2))) -
		(Nu_div_Rho * 2.0 * tan_lat_2) + tan_lat_4);
	
	T Term4(pow(omega, 8.0) / 40320 * Nu * sin_lat * pow(cos_lat, 7.0));
	Term4 *= (1385.0 - (3111.0 * tan_lat_2) + (543.0 * tan_lat_4) - (pow(tan(latitude), 6.0)));

	*northing = (K0 * (m + Term1 + Term2 + Term3 + Term4)) + FALSE_N;
}

// returns lat/long values in radians
template <class T, typename U>
void GridToGeo(const T& easting, const T& northing, const U& zone, T* latitude, T* longitude, 
			   const T& a, const T& inv_f, // ellipsoid parameters
			   const T& FALSE_E, const T& FALSE_N, const T& kO, const T& lcmZ1, const T& zW)	// projecton parameters
{
	T f = 1 / inv_f;
	T b = a * (1 - f);
	T e2 = (2 * f) - (f * f);
	//T e = sqrt(e2);
	//T Seconde2 = e2 / (1 - e2);
	//T Seconde = sqrt(Seconde2);
	T n = (a - b) / (a + b);
	T n2 = pow(n, 2.0);
	T n3 = pow(n, 3.0);
	T n4 = pow(n, 4.0);
	T G = a * (1 - n) * (1 - n2);
	G *= (1 + (9 * n2 / 4) + (225 * n4 / 64));
	G *= (PI / 180.);

	T ePrime = easting - FALSE_E;
	T nPrime = northing - FALSE_N;
	T m = nPrime / kO;
	T sigma = (m * PI) / (180 * G);
	
	T latPrime = sigma + (( (3 * n / 2) - (27 * n3 / 32) ) * sin(2 * sigma));
	latPrime += ( (21 * n2 / 16) - (55 * n4 / 32) ) * sin(4 * sigma);
	latPrime += (151 * n3 / 96) * sin(6 * sigma);
	latPrime += (1097 * n4 / 512) * sin(8 * sigma);

	T rho = a * (1 - e2) / pow( (1 - (e2 * pow(sin(latPrime), 2.0))), 1.5);
	T nu = a / pow( (1 - (e2 * pow(sin(latPrime), 2.0))), 0.5);
	T num1 = tan(latPrime) / (kO * rho);
	T x = ePrime / (kO * nu);
	
	T term1 = num1 * x * ePrime / 2;
	T term2 = num1 * ePrime * pow(x, 3.0) / 24;
	term2 *= ( (-4 * pow((nu / rho), 2.0)) + (9 * nu / rho * (1 - pow((tan(latPrime)), 2.0))) + (12 * pow((tan(latPrime)), 2.0))  );
	T term3 = num1 * ePrime * pow(x, 5.0) / 720;
	term3 *= ( 
		(8 * pow((nu / rho), 4.0) * (11 - (24 * pow((tan(latPrime)), 2.0)))) - 
		(12 * pow((nu / rho), 3.0) * (21 - (71 * pow((tan(latPrime)), 2.0)))) +
		(15 * pow((nu / rho), 2.0) * (15 - (98 * pow((tan(latPrime)), 2.0)) + (15 * pow((tan(latPrime)), 4.0)))) +
		(180 * (nu / rho) * ((5 * pow((tan(latPrime)), 2.0))-(3 * pow((tan(latPrime)), 4.0))))+
		(360 * pow((tan(latPrime)), 4.0))
		);
	T term4 = num1 * ePrime * pow(x, 7.0) / 40320;
	term4 *= (1385 +
		(3633 * pow((tan(latPrime)), 2.0)) +
		(4095 * pow((tan(latPrime)), 4.0)) +
		(1575 * pow((tan(latPrime)), 6.0))
		);

	// Store radians value of Geodetic Latitude
	*latitude = latPrime - term1 + term2 - term3 + term4;

	// Compute Geodetic Longitude of the Central Meridian of pGridValue->dNum3 (the UTM Zone) in radians
	T centralMeridian = ((zone * zW) + lcmZ1 - zW) * PI / 180;
	
	num1 = 1 / (cos(latPrime));

	term1 = x * num1;
	term2 = ( pow(x, 3.0) / 6 * num1 * ((nu / rho) + (2. * tan(latPrime) * tan(latPrime))) );
	term3 = ( pow(x, 5.0) / 120 * num1 * (
		(-4 * pow((nu / rho), 3.0) * (1 - (6 * tan(latPrime) * tan(latPrime)))) +
		(pow((nu / rho), 2.0) * (9 - (68 * tan(latPrime) * tan(latPrime)))) +
		(72 * (nu / rho) * tan(latPrime) * tan(latPrime)) +
		(24 * tan(latPrime) * tan(latPrime) * tan(latPrime) * tan(latPrime))
		));
	term4 = ( pow(x, 7.0) / 5040 * num1 * (
		61 + (662 * tan(latPrime) * tan(latPrime)) +
		(1320 * tan(latPrime) * tan(latPrime) * tan(latPrime) * tan(latPrime)) +
		(720 * tan(latPrime) * tan(latPrime) * tan(latPrime) * tan(latPrime) * tan(latPrime) * tan(latPrime))
		));

	// Store radians value of Geodetic Longitude
	*longitude = centralMeridian + term1 - term2 + term3 - term4;

}

// Great Circle Distance
template <class T>
T GreatCircleDistance(const T& dLatitudeAT, const T& dLongitudeAT, const T& dLatitudeTO, const T& dLongitudeTO)
{
	T deltaLatitude(dLatitudeTO - dLatitudeAT);
	T deltaLongitude(dLongitudeTO - dLongitudeAT);
	T a(sin(deltaLatitude / 2) * sin(deltaLatitude / 2) + cos(dLatitudeAT) * cos(dLatitudeTO) * sin(deltaLongitude / 2) * sin(deltaLongitude / 2));
	T c(2 * atan2(sqrt(a), sqrt(1 - a)));
	return c * T(6372797.);
}

// Rigorous Geodesic via Robbins' formula
// Robbins, A. R. (1962). �Long lines on the spheroid.� Surv. Rev., XVI(125), 301�309.
template <class T>
T RobbinsReverse(const T& dLatitudeAT, const T& dLongitudeAT, const T& dLatitudeTO, const T& dLongitudeTO, T* pAzimuth, const CDnaEllipsoid* ellipsoid)
{
	T s, z, c, g, h, h2, chi;
	T sinsigma, sigma, sigma2, sigma3, sigma4, sigma5;
		
	T dPVertA(primeVertical(ellipsoid, dLatitudeAT));
	T dPVertB(primeVertical(ellipsoid, dLatitudeTO));
	T tanzeta2 = (1 - ellipsoid->GetE1sqd()) * tan(dLatitudeTO) + ellipsoid->GetE1sqd() * dPVertA * sin(dLatitudeAT) / (dPVertB * cos(dLatitudeTO));
	T tau1 = cos(dLatitudeAT) * tanzeta2 - sin(dLatitudeAT) * cos(dLongitudeTO - dLongitudeAT);
	T tanazimuthAB;
	if (fabs (dLongitudeTO - dLongitudeAT) < PRECISION_1E15)
		tanazimuthAB = 0.0;
	else
		tanazimuthAB = sin (dLongitudeTO - dLongitudeAT) / tau1;

	// Compute the azimuth, check for the correct sign, quadrant etc. - 
	*pAzimuth = atan(tanazimuthAB);
	if ((*pAzimuth) < 0.0)
		*pAzimuth += PI;
	if (dLongitudeTO < dLongitudeAT)
		*pAzimuth += PI;
	if ((fabs (*pAzimuth) < PRECISION_1E15) && (dLatitudeTO < dLatitudeAT))
		*pAzimuth += PI;

	// Check here for the sign of the computed azimuth - eg to add pi or 2pi etc. 
	s = sin (*pAzimuth);
	z = atan (tanzeta2);
	c = cos (z);

	// If sin (alpha12) is close to zero then chi is calculated one way 
	// otherwise it is calculated differently.  An arbitrary 
	// "boundary value" of 0.2 has been used. 
	if (fabs (s) < 0.2)
		chi = tau1 / cos(*pAzimuth);
	else
		chi = sin(dLongitudeTO - dLongitudeAT) / s;
	sinsigma = chi * c;
	sigma = asin (sinsigma);
	sigma2 = sigma * sigma;
	sigma3 = sigma2 * sigma;
	sigma4 = sigma2 * sigma2;
	sigma5 = sigma3 * sigma2;
	g = ellipsoid->GetE2() * sin (dLatitudeAT);
	h = ellipsoid->GetE2() * cos (dLatitudeAT) * cos (*pAzimuth);
	h2 = h * h;
	return (dPVertA * sigma * (1 - sigma2 * h2 * (1 - h2) / 6 + sigma3 * g * h * (1 - 2 * h2) / 8 + sigma4 * (h2 * (4 - 7 * h2) - 3 * g * g * (1 - 7 * h2)) / 120 - sigma5 * g * h / 48));
}


template <class T>
void VincentyDirect(const T& dLatitudeAT, const T& dLongitudeAT, const T& dAzimuth, const T& dDistance, 
					T *dLatitudeTO, T *dLongitudeTO, const CDnaEllipsoid* ellipsoid)
{
	// calculate fundamentals
	T f = ellipsoid->GetFlattening();
	T b = ellipsoid->GetSemiMinor();

	// parametric latitude of P'
	T tanUI = (1.0 - f) * tan(dLatitudeAT);		
	// angular distance
	T tanSigma1 = tanUI / cos(dAzimuth);		
	// parametric latitude of the geodesic vertex, or
	// azimuth of the geodesic at the equator
	T sinAlpha = cos(atan(tanUI)) * sin(dAzimuth);
	T Alpha = asin(sinAlpha);
	T cosAlpha = cos(Alpha);
	// geodesic constant
	T u2 = pow(cosAlpha, 2.0) * (pow(ellipsoid->GetSemiMajor(), 2.0) - pow(b, 2.0)) / pow(b, 2.0);
	
	// Vincenty's constants A' and B'
	T A = 1.0 + (u2/16384.0) * (4096.0 + (u2 * (-768.0 + (u2 * (320.0 - (175.0 * u2))))));
	T B = (u2/1024.0) * (256.0 + (u2 * (-128.0 + (u2 * (74.0 - (47.0 * u2))))));
	
	T Sigma = dDistance / (b * A);
	T twoSigmam, deltaSigma, SigmaDiff(99.);
	
	// iterate until no signigicant change in sigma
	for (UINT16 i(0); i<10; i++)
	{
		twoSigmam = (2.0 * atan(tanSigma1)) + Sigma;
		deltaSigma = B * sin(Sigma) * (cos(twoSigmam) + (B / 4.0 * ((cos(Sigma) * (-1.0 + (2.0 * pow(cos(twoSigmam), 2.0)))) - (B / 6.0 * cos(twoSigmam) * (-3.0 + (4.0 * pow(sin(Sigma), 2.0))) * ((-3.0 + (4.0 * pow(cos(twoSigmam), 2.0))))))));
		SigmaDiff = Sigma;
		Sigma = (dDistance / (b * A)) + deltaSigma;
		SigmaDiff -= Sigma;

		if (fabs(SigmaDiff) < PRECISION_1E16)
			break;
	}
	
	// latitude of new position
	*dLatitudeTO = atan2(((sin(atan(tanUI)) * cos(Sigma)) + (cos(atan(tanUI)) * sin(Sigma) * cos(dAzimuth))), ((1.0 - f) * pow((pow(sinAlpha, 2.0) + pow(((sin(atan(tanUI)) * sin(Sigma)) - (cos(atan(tanUI)) * cos(Sigma) * cos(dAzimuth))), 2.0)), 0.5)));
	
	T Lambda = atan2((sin(Sigma) * sin(dAzimuth)), ((cos(atan(tanUI)) * cos(Sigma)) - (sin(atan(tanUI))*sin(Sigma)*cos(dAzimuth))));
	T C = (f / 16.0) * pow(cosAlpha, 2.0) * (4.0 + (f * (4.0 - (3.0 * pow(cosAlpha, 2.0)))));
	T Omega = Lambda - ((1.0 - C) * f * sinAlpha * (Sigma + (C * sin(Sigma) * (cos(twoSigmam) + (C * cos(Sigma) * (-1 + (2 * pow(cos(twoSigmam), 2.0))))))));
	
	// longitude of new position
	*dLongitudeTO = dLongitudeAT + Omega;
}

template <class T>
void ComputeLocalElements3D(const T X1, const T Y1, const T Z1,
			const T X2, const T Y2, const T Z2, 
			const T currentLatitude, const T currentLongitude,
			T* local_12e, T* local_12n, T* local_12up)
{
	// 1->2
	T dX12(X2 - X1);
	T dY12(Y2 - Y1);
	T dZ12(Z2 - Z1);

	// helpers
	T sin_lat(sin(currentLatitude));
	T cos_lat(cos(currentLatitude));
	T sin_long(sin(currentLongitude));
	T cos_long(cos(currentLongitude));


	*local_12e = -sin_long * dX12 + cos_long * dY12;
	*local_12n = -sin_lat * cos_long * dX12 - 
		sin_lat * sin_long * dY12 +
		cos_lat * dZ12;
	*local_12up = cos_lat * cos_long * dX12 +
		cos_lat * sin_long * dY12 +
		sin_lat * dZ12;
}

template <class T>
void ComputeLocalElements2D(const T X1, const T Y1, const T Z1,
	const T X2, const T Y2, const T Z2,
	const T currentLatitude, const T currentLongitude,
	T* local_12e, T* local_12n)
{
	// 1->2
	T dX12(X2 - X1);
	T dY12(Y2 - Y1);
	T dZ12(Z2 - Z1);

	// helpers
	T sin_lat(sin(currentLatitude));
	T cos_lat(cos(currentLatitude));
	T sin_long(sin(currentLongitude));
	T cos_long(cos(currentLongitude));


	*local_12e = -sin_long * dX12 + cos_long * dY12;
	*local_12n = -sin_lat * cos_long * dX12 -
		sin_lat * sin_long * dY12 +
		cos_lat * dZ12;
}

template <class T>
T Direction(const T local_12e, const T local_12n)
{
	// "computed" direction 1->2
	T direction12;

	if (fabs(local_12e) < fabs(local_12n))
		direction12 = atan_2(local_12e, local_12n);
	else
		direction12 = HALF_PI - atan_2(local_12n, local_12e);

	if (direction12 < 0)
		direction12 += TWO_PI;
	
	return direction12;
}

template <class T>
T Direction(const T X1, const T Y1, const T Z1,
			const T X2, const T Y2, const T Z2, 
			const T currentLatitude, const T currentLongitude,
			T* local_12e, T* local_12n)
{
	ComputeLocalElements2D(X1, Y1, Z1, X2, Y2, Z2, currentLatitude, currentLongitude,
		local_12e, local_12n);

	return Direction(*local_12e, *local_12n);
}

template <class T>
// helper function
T Direction(const T X1, const T Y1, const T Z1,
			const T X2, const T Y2, const T Z2, 
			const T currentLatitude, const T currentLongitude)
{
	T local_12e, local_12n;
	
	return Direction(X1, Y1, Z1,
		X2, Y2, Z2, 
		currentLatitude, currentLongitude,
		&local_12e, &local_12n);
}

template <class T>
T HorizontalAngle(const T X1, const T Y1, const T Z1, 
				  const T X2, const T Y2, const T Z2, 
				  const T X3, const T Y3, const T Z3, 
				  const T currentLatitude, const T currentLongitude,
				  T* direction12, T* direction13,
				  T* local_12e, T* local_12n, T* local_13e, T* local_13n)
{
	// compute vectors [1->2] & [1->3] in the local reference frame
	//
	// 1->2
	*direction12 = Direction(X1, Y1, Z1, X2, Y2, Z2, currentLatitude, currentLongitude, local_12e, local_12n);
	*direction13 = Direction(X1, Y1, Z1, X3, Y3, Z3, currentLatitude, currentLongitude, local_13e, local_13n);
	
	if (*direction12 > *direction13)
		*direction13 += TWO_PI;

	// angle 123
	T angle = *direction13 - *direction12;

	return angle;
}

template <class T>
// helper function
T HorizontalAngle(const T X1, const T Y1, const T Z1, 
				  const T X2, const T Y2, const T Z2, 
				  const T X3, const T Y3, const T Z3, 
				  const T currentLatitude, const T currentLongitude,
				  T* direction12, T* direction13)
{
	T local_12e, local_12n, local_13e, local_13n;

	return HorizontalAngle(X1, Y1, Z1, 
		X2, Y2, Z2, 
		X3, Y3, Z3, 
		currentLatitude, currentLongitude,
		direction12, direction13,
		&local_12e, &local_12n, &local_13e, &local_13n);
}

template <class T>
void CartesianElementsFromInstrumentHeight(const T height, T* dX, T* dY, T* dZ, 
				  const T Latitude, const T Longitude)
{
	// Use rotation matrix for local vector -> cartesian vector, whereby
	// local elements for e and n are zero
	*dX = cos(Latitude) * cos(Longitude) * height;
	*dY = cos(Latitude) * sin(Longitude) * height;
	*dZ = sin(Latitude) * height;
}
	
template <class T>
// The return value is the true vertical between the (local) horizontal plane
// and the instrument-target vector.  The local_12e/n/up elements represent 
// the geometric difference between the two stations
T VerticalAngle(const T& local_12e, const T& local_12n, const T& local_12up)
{
	return atan2(local_12up, sqrt((local_12e * local_12e) + (local_12n * local_12n)));
}

template <class T>
// The return value is the true vertical between the (local) horizontal plane
// and the instrument-target vector.  The local_12e/n/up elements represent 
// the geometric difference between the two stations
T VerticalAngle(const T X1, const T Y1, const T Z1, 
				  const T X2, const T Y2, const T Z2, 
				  const T latitude1, const T longitude1,
				  const T latitude2, const T longitude2,
				  const T instrumentHeight, const T targetHeight,
				  T* local_12e, T* local_12n, T* local_12up)
{
	// helpers
	T sin_lat1(sin(latitude1));
	T cos_lat1(cos(latitude1));
	T sin_long1(sin(longitude1));
	T cos_long1(cos(longitude1));

	// Compute cartesian elements dX, dY, dZ for instrument to target
	// First, compute cartesian vector for both instrument and target
	T dXih, dYih, dZih, dXth, dYth, dZth;
	CartesianElementsFromInstrumentHeight(instrumentHeight,
		&dXih, &dYih, &dZih, latitude1, longitude1);
	CartesianElementsFromInstrumentHeight(targetHeight,
		&dXth, &dYth, &dZth, latitude2, longitude2);

	T dX12(X2 - X1 + dXth - dXih);
	T dY12(Y2 - Y1 + dYth - dYih);
	T dZ12(Z2 - Z1 + dZth - dZih);

	// compute local reference frame elements (station1 to station2)
	*local_12e = -sin_long1 * dX12 + cos_long1 * dY12;
	*local_12n = -sin_lat1 * cos_long1 * dX12 -
					sin_lat1 * sin_long1 * dY12 +
					cos_lat1 * dZ12;
	*local_12up = cos_lat1 * cos_long1 * dX12 +
					cos_lat1 * sin_long1 * dY12 +
					sin_lat1 * dZ12;

	// compute angle (instrument to target)
	return VerticalAngle(*local_12e, *local_12n, *local_12up);
	//return atan2((*local_12up), sqrt(((*local_12e) * (*local_12e)) + ((*local_12n) * (*local_12n))));
	//////////////////////////////////////////////////////
}

template <class T>
// helper function
T VerticalAngle(const T X1, const T Y1, const T Z1, 
				  const T X2, const T Y2, const T Z2, 
				  const T latitude1, const T longitude1,
				  const T latitude2, const T longitude2,
				  const T instrumentHeight, const T targetHeight)
{
	T local_12e, local_12n, local_12up;

	return VerticalAngle(
		X1, Y1, Z1, 
		X2, Y2, Z2, 
		latitude1, longitude1,
		latitude2, longitude2,
		instrumentHeight, targetHeight,
		&local_12e, &local_12n, &local_12up);
}


template <class T>
// The return value is the true zenith distance between the ellipsoid normal
// and the instrument-target vector.  The local_12e/n/up elements represent 
// the geometric difference between the two stations
T ZenithDistance(const T X1, const T Y1, const T Z1, 
				  const T X2, const T Y2, const T Z2, 
				  const T latitude1, const T longitude1,
				  const T latitude2, const T longitude2,
				  const T instrumentHeight, const T targetHeight,
				  T* local_12e, T* local_12n, T* local_12up)
{
	// helpers
	T sin_lat1(sin(latitude1));
	T cos_lat1(cos(latitude1));
	T sin_long1(sin(longitude1));
	T cos_long1(cos(longitude1));

	// Compute cartesian elements dX, dY, dZ for instrument to target
	// First, compute cartesian vector for both instrument and target
	T dXih, dYih, dZih, dXth, dYth, dZth;
	CartesianElementsFromInstrumentHeight(instrumentHeight,
		&dXih, &dYih, &dZih, latitude1, longitude1);
	CartesianElementsFromInstrumentHeight(targetHeight,
		&dXth, &dYth, &dZth, latitude2, longitude2);

	T dX12(X2 - X1 + dXth - dXih);
	T dY12(Y2 - Y1 + dYth - dYih);
	T dZ12(Z2 - Z1 + dZth - dZih);

	// compute local reference frame elements (station1 to station2)
	*local_12e = -sin_long1 * dX12 + cos_long1 * dY12;
	*local_12n = -sin_lat1 * cos_long1 * dX12 -
					sin_lat1 * sin_long1 * dY12 +
					cos_lat1 * dZ12;
	*local_12up = cos_lat1 * cos_long1 * dX12 +
					cos_lat1 * sin_long1 * dY12 +
					sin_lat1 * dZ12;

	// compute angle (instrument to target)
	return atan2(sqrt((*local_12e) * (*local_12e) + (*local_12n) * (*local_12n)), (*local_12up));
	//////////////////////////////////////////////////////
}

template <class T>
// helper function
T ZenithDistance(const T X1, const T Y1, const T Z1, 
				  const T X2, const T Y2, const T Z2, 
				  const T latitude1, const T longitude1,
				  const T latitude2, const T longitude2,
				  const T instrumentHeight, const T targetHeight)
{
	T local_12e, local_12n, local_12up;

	return ZenithDistance(
		X1, Y1, Z1, X2, Y2, Z2, 
		latitude1, longitude1,
		latitude2, longitude2,
		instrumentHeight, targetHeight,
		&local_12e, &local_12n, &local_12up);
}
	

template <class T>
T EllipsoidHeight(const T X, const T Y, const T Z, 
				  const T latitude, T* nu, T* Zn,
				  const CDnaEllipsoid* ellipsoid)
{
	*nu = primeVertical(ellipsoid, latitude);
	// Zn is the z coordinate element of the point on the z-axis 
	// which intersects with the the normal at the given Latitude
	*Zn = ellipsoid->GetE1sqd() * (*nu) * sin(latitude);

	return sqrt(X*X + Y*Y + pow(Z+(*Zn), (int)2)) - (*nu);
}
	
template <class T>
T EllipsoidHeightDifference(const T X1, const T Y1, const T Z1, 
				  const T X2, const T Y2, const T Z2, 
				  const T Latitude1, const T Latitude2,
				  T* h1, T* h2, T* nu1, T* nu2, T* Zn1, T* Zn2, 
				  const CDnaEllipsoid* ellipsoid)
{
	return ((*h2 = EllipsoidHeight(X2, Y2, Z2, Latitude2, nu2, Zn2, ellipsoid)) - 
		(*h1 = EllipsoidHeight(X1, Y1, Z1, Latitude1, nu1, Zn1, ellipsoid)));
}

template <class T>
T magnitude(const T X1, const T Y1, const T Z1, 
				  const T X2, const T Y2, const T Z2)
{
	return sqrt(((X2 - X1) * (X2 - X1)) + ((Y2 - Y1) * (Y2 - Y1)) + ((Z2 - Z1) * (Z2 - Z1)));
}

template <class T>
T magnitude(const T dX, const T dY, const T dZ)
{
	return sqrt((dX * dX) + (dY * dY) + (dZ * dZ));
}

template <class T>
T magnitude(const T X1, const T N1, const T X2, const T N2)
{
	return sqrt(((X2 - X1) * (X2 - X1)) + ((N2 - N1) * (N2 - N1)));
}

template <class T>
T magnitude(const T dX, const T dN)
{
	return sqrt((dX * dX) + (dN * dN));
}

template <class T>
T EllipsoidChordDistance(const T X1, const T Y1, const T Z1, 
				  const T X2, const T Y2, const T Z2, 
				  const T latitude1, const T latitude2,
				  const T Height1, const T Height2,
				  T* dX, T* dY, T* dZ, 
				  const CDnaEllipsoid* ellipsoid)
{
	T nu1(primeVertical(ellipsoid, latitude1));
	T nu2(primeVertical(ellipsoid, latitude2));

	T scale1(nu1 / (nu1 + Height1));
	T scale2(nu2 / (nu2 + Height2));

	// Zn1,2 is the z coordinate element of the point on the z-axis 
	// which intersects with the the normal at the given Latitude
	T Zn1(ellipsoid->GetE1sqd() * nu1 * sin(latitude1));
	T Zn2(ellipsoid->GetE1sqd() * nu2 * sin(latitude2));

	// station 1
	T x1(X1 * scale1);
	T y1(Y1 * scale1);
	T z1((Z1 + Zn1) * scale1 - Zn1);
 
	// station 2
	T x2(X2 * scale2);
	T y2(Y2 * scale2);
	T z2((Z2 + Zn2) * scale2 - Zn2);

	*dX = x2 - x1;
	*dY = y2 - y1;
	*dZ = z2 - z1;

	return magnitude(*dX, *dY, *dZ);
}

template <class T>
T RadiusCurvatureInChordDirection(const T X1, const T Y1, const T Z1, 
			 const T X2, const T Y2, const T Z2, 				  
			 const T latitude1, const T longitude1, const T latitude2,
			 const CDnaEllipsoid* ellipsoid)
{
	T nu, rho;
	primeVerticalandMeridian(ellipsoid, average(latitude1, latitude2), nu, rho);

	T local_12e, local_12n;
	T direction12(Direction(X1, Y1, Z1, X2, Y2, Z2,
		latitude1, longitude1, &local_12e, &local_12n));				  
	T cos_dir(cos(direction12));
	T sin_dir(sin(direction12));
	return  rho * nu / ((nu * cos_dir * cos_dir) + (rho * sin_dir * sin_dir));
}

template <class T>
T EllipsoidArctoEllipsoidChord(const T arc, 
			 const T X1, const T Y1, const T Z1, 
			 const T X2, const T Y2, const T Z2, 				  
			 const T Latitude1, const T Longitude1, const T Latitude2,
			 const CDnaEllipsoid* ellipsoid)
{
	T r(RadiusCurvatureInChordDirection(X1, Y1, Z1, X2, Y2, Z2, Latitude1, Longitude1, Latitude2, ellipsoid));
	return 2.0 * r * sin(arc / 2.0 / r);
}


template <class T>
T EllipsoidChordtoEllipsoidArc(const T chord, 
			 const T X1, const T Y1, const T Z1, 
			 const T X2, const T Y2, const T Z2, 				  
			 const T Latitude1, const T Longitude1, const T Latitude2,
			 const CDnaEllipsoid* ellipsoid)
{
	T r(RadiusCurvatureInChordDirection(X1, Y1, Z1, X2, Y2, Z2, Latitude1, Longitude1, Latitude2, ellipsoid));
	return asin(chord / 2.0 / r) * 2.0 * r;
}

template <class T>
T EllipsoidArcDistance(
			 const T X1, const T Y1, const T Z1, 
			 const T X2, const T Y2, const T Z2, 				  
			 const T Latitude1, const T Longitude1, const T Latitude2,
			 const T Height1, const T Height2, 
			 const CDnaEllipsoid* ellipsoid)
{
	T dx, dy, dz, ellipsoid_chord;
	ellipsoid_chord = EllipsoidChordDistance<T>(
		X1, Y1, Z1,
		X2, Y2, Z2,
		Latitude1, Latitude2,
		Height1, Height2, &dx, &dy, &dz, ellipsoid);

	return EllipsoidChordtoEllipsoidArc<T>(
		ellipsoid_chord, 
		X1, Y1, Z1, 
		X2, Y2, Z2, 				  
		Latitude1, Longitude1, Latitude2,
		ellipsoid);
}

template <class T>
T MSLChordtoMSLArc(const T chord, 
			 const T latitude1, const T latitude2,
			 const T N1, const T N2,
			 const CDnaEllipsoid* ellipsoid)
{
	T nu, rho;
	primeVerticalandMeridian(ellipsoid, average(latitude1, latitude2), nu, rho);

	T r(sqrt(nu*rho) + average(N1, N2));
	return asin(chord / 2.0 / r) * 2.0 * r;
}


template <class T>
T MSLArctoMSLChord(const T arc, 
			 const T latitude1, const T latitude2,
			 const T N1, const T N2,
			 const CDnaEllipsoid* ellipsoid)
{
	T nu, rho;
	primeVerticalandMeridian(ellipsoid, average(latitude1, latitude2), nu, rho);

	T r(sqrt(nu*rho) + average(N1, N2));
	return 2.0 * r * sin(arc / 2.0 / r);
}


template <class T>
T MSLChordtoEllipsoidChord(const T msl_chord, 
						   const T Latitude1, const T Latitude2,
						   const T N1, const T N2,
						   const CDnaEllipsoid* ellipsoid)
{
	T ellipsoid_chord(msl_chord * msl_chord);
	ellipsoid_chord -= pow(N2 - N1, (int)2);

	T meanLat(average(Latitude1, Latitude2));
	ellipsoid_chord /= 1. + N1 / averageRadiusofCurvature(ellipsoid, meanLat);
	ellipsoid_chord /= 1. + N2 / averageRadiusofCurvature(ellipsoid, meanLat);
	return sqrt(ellipsoid_chord);
}

template <class T>
T MSLArctoEllipsoidChord(const T msl_arc, 
						   const T Latitude1, const T Latitude2,
						   const T N1, const T N2,
						   const CDnaEllipsoid* ellipsoid)
{
	// 1. Convert MSL Arc -> MSL Chord
	T msl_chord = MSLArctoMSLChord<T>(msl_arc, 
		Latitude1, Latitude2,
		N1, N2,
		ellipsoid);
	
	// 2. Convert MSL Chord -> Ellipsoid Chord
	return MSLChordtoEllipsoidChord<T>(msl_chord, 
		Latitude1, Latitude2,
		N1, N2,
		ellipsoid);
}

template <class T>
T EllipsoidChordtoMSLChord(const T ellipsoid_chord, 
						   const T Latitude1, const T Latitude2,
						   const T N1, const T N2,
						   const CDnaEllipsoid* ellipsoid)
{
	T msl_chord(ellipsoid_chord * ellipsoid_chord);
	T meanLat(average(Latitude1, Latitude2));
	
	msl_chord *= 1. + N1 / averageRadiusofCurvature(ellipsoid, meanLat);
	msl_chord *= 1. + N2 / averageRadiusofCurvature(ellipsoid, meanLat);
	msl_chord += pow(N2 - N1, (int)2);
	
	return sqrt(msl_chord);
}

template <class T>
T EllipsoidChordtoMSLArc(const T ellipsoid_chord, 
						   const T Latitude1, const T Latitude2,
						   const T N1, const T N2,
						   const CDnaEllipsoid* ellipsoid)
{
	// 1. Convert Ellipsoid Chord -> MSL Chord
	T msl_chord = EllipsoidChordtoMSLChord<T>(
		ellipsoid_chord, Latitude1, Latitude2,
		N1, N2, ellipsoid);

	// 2. Convert MSL Chord to MSL Arc
	return MSLChordtoMSLArc<T>(
		msl_chord, Latitude1, Latitude2,
		N1, N2, ellipsoid);
}

template <class T>
T MSLArcDistance(
			 const T X1, const T Y1, const T Z1, 
			 const T X2, const T Y2, const T Z2, 				  
			 const T Latitude1, const T Latitude2,
			 const T Height1, const T Height2, 
			 const T N1, const T N2,
			 const CDnaEllipsoid* ellipsoid)
{
	// 1. Calculate Ellipsoid Chord
	T dx, dy, dz, chord;
	chord = EllipsoidChordDistance<T>(
		X1, Y1, Z1,
		X2, Y2, Z2,
		Latitude1, Latitude2,
		Height1, Height2, &dx, &dy, &dz, ellipsoid);

	// 2. Convert Ellipsoid Chord -> MSL Chord
	chord = EllipsoidChordtoMSLChord<T>(
		chord, Latitude1, Latitude2,
		N1, N2, ellipsoid);

	// 3. Convert MSL Chord -> MSL Arc
	return MSLChordtoMSLArc<T>(
		chord, Latitude1, Latitude2,
		N1, N2, ellipsoid);
}
	

template <class T>
T LaplaceCorrection(const T azimuth, const T zenith,
					const T deflPrimeV, const T deflPrimeM,
					const T Latitude)
{
	return deflPrimeV * tan(Latitude) + ((deflPrimeM * sin(azimuth) - deflPrimeV * cos(azimuth)) / tan(zenith));	// cot(z) = 1/tan(z)
}

template <class T>
T ZenithDeflectionCorrection(const T azimuth, const T deflPrimeV, const T deflPrimeM)
{
	return deflPrimeM * cos(azimuth) + deflPrimeV * sin(azimuth);
}

template <class T>
T DirectionDeflectionCorrection(const T azimuth, const T zenith,
									  const T deflPrimeV, const T deflPrimeM)
{
	return (deflPrimeM * sin(azimuth) - deflPrimeV * cos(azimuth)) / tan(zenith);	// cot(z) = 1/tan(z)
}

template <class T>
T HzAngleDeflectionCorrection(const T azimuth12, const T zenith12,
									  const T azimuth13, const T zenith13,
									  const T deflPrimeV, const T deflPrimeM)
{
	return DirectionDeflectionCorrection(azimuth13, zenith13, deflPrimeV, deflPrimeM) -
		DirectionDeflectionCorrection(azimuth12, zenith12, deflPrimeV, deflPrimeM);
}

template <class T>
T HzAngleDeflectionCorrections(const T azimuth12, const T zenith12,
	const T azimuth13, const T zenith13,
	const T deflPrimeV, const T deflPrimeM, T& correction12, T& correction13)
{
	return (correction13 = DirectionDeflectionCorrection(azimuth13, zenith13, deflPrimeV, deflPrimeM)) -
		(correction12 = DirectionDeflectionCorrection(azimuth12, zenith12, deflPrimeV, deflPrimeM));
}
#endif /* DNATEMPLATEGEODESYFUNCS_H_ */
