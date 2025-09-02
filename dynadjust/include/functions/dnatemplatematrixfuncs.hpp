//============================================================================
// Name         : dnatemplatematrixfuncs.hpp
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
// Description  : Common functions involving matrix_2d operations
//============================================================================

#ifndef DNATEMPLATEMATRIXFUNCS_H_
#define DNATEMPLATEMATRIXFUNCS_H_

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
#include <include/config/dnatypes.hpp>
#include <include/math/dnamatrix_contiguous.hpp>
#include <include/measurement_types/dnameasurement.hpp>

using namespace dynadjust::datum_parameters;
using namespace dynadjust::math;

// msr_t_Iterator = std::vector<measurement_t>::iterator
template<typename msr_t_Iterator>
// Fills upper triangle
void GetDirectionsVarianceMatrix(msr_t_Iterator begin, matrix_2d* vmat)
{
	msr_t_Iterator bmsRecord(begin);
	UINT32 a, angle_count(bmsRecord->vectorCount2 - 1);		// number of directions excluding the RO
	UINT32 skip(0), ignored(bmsRecord->vectorCount1 - bmsRecord->vectorCount2);

	vmat->zero();
	vmat->redim(angle_count, angle_count);

	bmsRecord++;

	for (a=0; a<angle_count; ++a)
	{
		// Ignored direction?
		if (bmsRecord->ignore)
		{
			while (skip < ignored)
			{
				skip++;
				bmsRecord++;
				if (!bmsRecord->ignore)
					break;
			}
		}

		vmat->put(a, a, bmsRecord->scale2);				// derived angle variance
		if (a+1 < angle_count)
			vmat->put(a, a+1, bmsRecord->scale3);		// derived angle covariance
		bmsRecord++;
	}
}


// M = measurement_t, U = matrix_2d
// Fills upper triangle
template<typename msr_t_Iterator>
void GetGPSVarianceMatrix(const msr_t_Iterator begin, matrix_2d* vmat) 
{
	msr_t_Iterator bmsRecord(begin);
	UINT32 variance_dim(bmsRecord->vectorCount1 * 3), covariance_dim, cov;
	vmat->zero();
	vmat->redim(variance_dim, variance_dim);	

	for (UINT32 var(0), cov_elem; var<variance_dim; var+=3)
	{
		covariance_dim = bmsRecord->vectorCount2 * 3;

		vmat->put(var, var, (bmsRecord++)->term2);				// XX

		vmat->put(var, var+1, bmsRecord->term2);				// XY
		vmat->put(var+1, var+1, (bmsRecord++)->term3);			// YY
		
		vmat->put(var, var+2, bmsRecord->term2);				// XZ
		vmat->put(var+1, var+2, bmsRecord->term3);				// YZ
		vmat->put(var+2, var+2, (bmsRecord++)->term4);			// ZZ
		
		for (cov_elem=0; cov_elem<covariance_dim; cov_elem+=3)
		{
			cov = var + 3 + cov_elem;
			vmat->put(var, cov, bmsRecord->term1);				// m11
			vmat->put(var, cov+1, bmsRecord->term2);			// m12
			vmat->put(var, cov+2, (bmsRecord++)->term3);		// m13
			
			vmat->put(var+1, cov, bmsRecord->term1);			// m21
			vmat->put(var+1, cov+1, bmsRecord->term2);			// m22
			vmat->put(var+1, cov+2, (bmsRecord++)->term3);		// m23
			
			vmat->put(var+2, cov, bmsRecord->term1);			// m31
			vmat->put(var+2, cov+1, bmsRecord->term2);			// m32
			vmat->put(var+2, cov+2, (bmsRecord++)->term3);		// m33
		}
	}
}


// msr_t_Iterator = std::vector<measurement_t>::iterator
template<typename msr_t_Iterator>
// Sets values based on upper triangle
void SetDirectionsVarianceMatrix(msr_t_Iterator begin, const matrix_2d& vmat)
{
	msr_t_Iterator bmsRecord(begin);
	UINT32 a, angle_count(bmsRecord->vectorCount2 - 1);		// number of directions excluding the RO
	UINT32 skip(0), ignored(bmsRecord->vectorCount1 - bmsRecord->vectorCount2);

	bmsRecord->scale2 = 0.;	// variance (angle)
	bmsRecord->scale3 = 0.;	// covariance (angle)

	bmsRecord++;

	for (a=0; a<angle_count; ++a)
	{
		// Ignored direction?
		if (bmsRecord->ignore)
		{
			while (skip < ignored)
			{
				skip++;
				bmsRecord++;
				if (!bmsRecord->ignore)
					break;
			}
		}

		bmsRecord->scale2 = vmat.get(a, a);				// derived angle variance
		if (a+1 < angle_count)
			bmsRecord->scale3 = vmat.get(a, a+1);		// derived angle covariance
		else
			bmsRecord->scale3 = 0.;						// not necessary, but in the interest of 
														// preventing confusion, set to zero
		bmsRecord++;
	}
}

// msr_t_Iterator = std::vector<measurement_t>::iterator
template<typename msr_t_Iterator>
// Sets values based on upper triangle
void SetGPSVarianceMatrix(msr_t_Iterator begin, const matrix_2d& vmat)
{
	msr_t_Iterator bmsRecord(begin);
	UINT32 variance_dim(bmsRecord->vectorCount1 * 3), covariance_dim, cov;

	for (UINT32 var(0), cov_elem; var<variance_dim; var+=3)
	{
		covariance_dim = bmsRecord->vectorCount2 * 3;

		(bmsRecord++)->term2 = vmat.get(var, var);				// XX

		bmsRecord->term2 = vmat.get(var, var+1);				// XY
		(bmsRecord++)->term3 = vmat.get(var+1, var+1);			// YY
		
		bmsRecord->term2 = vmat.get(var, var+2);				// XZ
		bmsRecord->term3 = vmat.get(var+1, var+2);				// YZ
		(bmsRecord++)->term4 = vmat.get(var+2, var+2);			// ZZ
		
		for (cov_elem=0; cov_elem<covariance_dim; cov_elem+=3)
		{
			cov = var + 3 + cov_elem;
			bmsRecord->term1 = vmat.get(var, cov);				// m11
			bmsRecord->term2 = vmat.get(var, cov+1);			// m12
			(bmsRecord++)->term3 = vmat.get(var, cov+2);		// m13
			
			bmsRecord->term1 = vmat.get(var+1, cov);			// m21
			bmsRecord->term2 = vmat.get(var+1, cov+1);			// m22
			(bmsRecord++)->term3 = vmat.get(var+1, cov+2);		// m23
			
			bmsRecord->term1 = vmat.get(var+2, cov);			// m31
			bmsRecord->term2 = vmat.get(var+2, cov+1);			// m32
			(bmsRecord++)->term3 = vmat.get(var+2, cov+2);		// m33
		}
	}
}

template <class T>
void FormCarttoGeoRotationMatrix(const T& latitude, const T& longitude, const T& height, 
	matrix_2d& mrotations, const CDnaEllipsoid* ellipsoid, bool CLUSTER=false, const UINT32& n=0)
{
	if (!CLUSTER)
		mrotations.redim(3, 3);

	T coslat(cos(latitude));
	T sinlat(sin(latitude));
	T coslon(cos(longitude));
	T sinlon(sin(longitude));

	T term1_a(ellipsoid->GetSemiMajor() * ellipsoid->GetE1sqd());
	T one_minus_esq(1. - ellipsoid->GetE1sqd());
	
	T nu_plus_h(primeVertical(ellipsoid, latitude) + height);
	T nu_1minuse2_plus_h((primeVertical(ellipsoid, latitude) * (one_minus_esq) + height));
	T term1_b(term1_a * sinlat * coslat);
	T term1_c(pow((1. - ellipsoid->GetE1sqd() * sinlat * sinlat), 1.5));
	
	// set up Rotation matrix (geo to cart)
	mrotations.put(n,   n,   (term1_b * coslat * coslon / term1_c) - ((nu_plus_h) * sinlat * coslon));
	mrotations.put(n,   n+1, -(nu_plus_h) * coslat * sinlon);
	mrotations.put(n,   n+2, coslat * coslon);
	mrotations.put(n+1, n,   (term1_b * coslat * sinlon / term1_c) - ((nu_plus_h) * sinlat * sinlon));
	mrotations.put(n+1, n+1, (nu_plus_h) * coslat * coslon);
	mrotations.put(n+1, n+2, coslat * sinlon);
	mrotations.put(n+2, n,   (term1_b * one_minus_esq * sinlat / term1_c) + (nu_1minuse2_plus_h * coslat));
	mrotations.put(n+2, n+1, 0.);
	mrotations.put(n+2, n+2, sinlat);
}
	

template <class T>
void FormCarttoGeoRotationMatrix_Cluster(const matrix_2d& mpositions, matrix_2d& mrotations,
									  const CDnaEllipsoid* ellipsoid)
{
	mrotations.redim(mpositions.rows(), mpositions.rows());
	
	for (UINT32 i(0); i<mpositions.rows(); i+=3)
	{
		FormCarttoGeoRotationMatrix<T>(
			mpositions.get(i, 0),
			mpositions.get(i+1, 0),
			mpositions.get(i+2, 0),
			mrotations, ellipsoid, true, i);
	}
}
	
// Assumes design elements for station 1 and station 2 are always
// -1 and 1 respectively.
template <class T>
void Precision_Adjusted_GNSS_bsl(const matrix_2d& mvariances,
	const UINT32& stn1, const UINT32& stn2,
	matrix_2d* mvariances_mod, bool FILLLOWER=true)
{
	matrix_2d tmp(3, 6);
	mvariances_mod->zero();
	UINT32 i, j, k;

	// 1. Form A * V
	for (i=0, j=0; i<3; ++i, ++j)
	{
		// variance-covariance
		for (j=0; j<3; ++j)
		{
			// stn 11 variance
			tmp.elementadd(i, j, -mvariances.get(stn1+i, stn1+j));
			// stn 21 covariance
			tmp.elementadd(i, j, mvariances.get(stn2+i, stn1+j));			
		}
		k=j;
		// variance-covariance
		for (j=0; j<3; ++j, ++k)
		{
			// stn 12 covariance
			tmp.elementadd(i, k, -mvariances.get(stn1+i, stn2+j));
			// stn 22 variance
			tmp.elementadd(i, k, mvariances.get(stn2+i, stn2+j));			
		}
	}

	//tmp.trace("A x V", "%.16G ");

	// 2. Form AV * AT (upper triangular variance-covariance)
	for (i=0; i<3; ++i)		
		for (j=i; j<3; ++j)
			// Sum the variances & covariances
			mvariances_mod->put(i, j, tmp.get(i, j+3) - tmp.get(i, j));
		
	if (FILLLOWER)
		mvariances_mod->filllower();

	//mvariances_mod->trace("AV X At", "%.16G ");
}

template <class T>
void Prpagate_Variances_Geo_Cart(const matrix_2d& mvariances, matrix_2d mrotations, matrix_2d* mvariances_mod, bool FORWARD=true)
{
	// the rotation matrix is in the direction geo to cart (forward)
	// so to go cart to geo, perform inverse
	if (!FORWARD)
		mrotations = mrotations.sweepinverse();		// allows negative diagonal terms
	
	matrix_2d mV(mrotations);
	//mV.multiply(mvariances);		// original variance matrix
	mV.multiply("N", mvariances, "N");		// original variance matrix
	//mvariances_mod->multiply_square_t(mV, mrotations);
	mvariances_mod->multiply(mV, "N", mrotations, "T");
}


template <class T>
void PropagateVariances_GeoCart(const matrix_2d mvariances, matrix_2d* mvariances_mod, 
								const T& latitude, const T& longitude, const T& height, 
								matrix_2d& mrotations, 
								const CDnaEllipsoid* ellipsoid, bool GEO_TO_CART, bool CALCULATE_ROTATIONS)
{

	if (CALCULATE_ROTATIONS)
		FormCarttoGeoRotationMatrix<T>(latitude, longitude, height, mrotations, ellipsoid);

	Prpagate_Variances_Geo_Cart<T>(mvariances, mrotations, mvariances_mod, GEO_TO_CART);
}
	

template <class T>
void PropagateVariances_GeoCart(const matrix_2d mvariances, matrix_2d* mvariances_mod, 
								const T& latitude, const T& longitude, const T& height, 
								const CDnaEllipsoid* ellipsoid, bool GEO_TO_CART)
{

	matrix_2d mrotations;
	PropagateVariances_GeoCart<T>(mvariances, mvariances_mod, 
		latitude, longitude, height, 
		mrotations, 
		ellipsoid, GEO_TO_CART, 
		true);		// Calculate rotations
}
	

template <class T>
void PropagateVariances_GeoCart_Cluster(const matrix_2d& mvariances, matrix_2d* mvariances_mod, 
								const matrix_2d& mpositions, matrix_2d& mrotations, 
								const CDnaEllipsoid* ellipsoid, bool GEO_TO_CART, bool CALCULATE_ROTATIONS)
{
	if (CALCULATE_ROTATIONS)
		FormCarttoGeoRotationMatrix_Cluster<T>(mpositions, mrotations, ellipsoid);
	
	Prpagate_Variances_Geo_Cart<T>(mvariances, mrotations, mvariances_mod, GEO_TO_CART);
}
	

template <class T>
void PropagateVariances_GeoCart_Cluster(const matrix_2d& mvariances, matrix_2d* mvariances_mod, 
								const matrix_2d& mpositions_rad, 
								const CDnaEllipsoid* ellipsoid, bool GEO_TO_CART)
{
	matrix_2d mrotations;
	FormCarttoGeoRotationMatrix_Cluster<T>(mpositions_rad, mrotations, ellipsoid);
	
	Prpagate_Variances_Geo_Cart<T>(mvariances, mrotations, mvariances_mod, GEO_TO_CART);
}
	
template <class T>
void ScaleMatrix(const matrix_2d mvariances, matrix_2d* mvariances_mod, const matrix_2d& scalars)
{
	//matrix_2d mV(mvariances_mod->multiply(scalars, mvariances));
	matrix_2d mV(mvariances_mod->multiply(scalars, "N", mvariances, "N"));
	//mvariances_mod->multiply_square_t(mV, scalars);
	mvariances_mod->multiply(mV, "N", scalars, "T");
}


template <class T>
void ScaleGPSVCV(const matrix_2d& mvariances, matrix_2d* mvariances_mod, 
						const T& latitude, const T& longitude, const T& height, 
						const CDnaEllipsoid* ellipsoid, 
						const T& pScale, const T& lScale, const T& hScale)
{
	matrix_2d mrotations(3, 3);
	PropagateVariances_GeoCart<T>(mvariances, mvariances_mod,
		latitude, longitude, height, mrotations, ellipsoid,
		false, 		// Geographic -> Cartesian ?
		true);		// create the rotation matrix

	matrix_2d var_scalars(3, 3);
	var_scalars.put(0, 0, sqrt(pScale));
	var_scalars.put(1, 1, sqrt(lScale));
	var_scalars.put(2, 2, sqrt(hScale));
	ScaleMatrix<T>(*mvariances_mod, mvariances_mod, var_scalars);

	PropagateVariances_GeoCart<T>(*mvariances_mod, mvariances_mod,
		latitude, longitude, height, mrotations, ellipsoid,
		true,		// Geographic -> Cartesian ?
		false);		// don't create a rotation matrix
}

template <class T>
// coordType is passed so that ScaleGPSVCV_Cluster knows whether
// mvariances needs to be propagated to geographic first.  Hence,
// if coordType == LLH_type_i, no propagation is undertaken
void ScaleGPSVCV_Cluster(const matrix_2d& mvariances, matrix_2d* mvariances_mod, 
						const matrix_2d& mpositions, const CDnaEllipsoid* ellipsoid, 
						const T& pScale, const T& lScale, const T& hScale, _COORD_TYPE_ coordType=XYZ_type_i)
{
	matrix_2d mrotations;
	// form rotation matrix
	FormCarttoGeoRotationMatrix_Cluster<T>(mpositions, mrotations, ellipsoid);

	// Don't propagate if already in geographic
	if (coordType == XYZ_type_i)
		// propagate variances in cartesian system to geographic 
		PropagateVariances_GeoCart_Cluster<T>(mvariances, mvariances_mod,
			mpositions, mrotations, ellipsoid,
			false, 		// Cartesian -> Geographic
			false);		// don't create a rotation matrix

	// scale matrix
	matrix_2d var_scalars(mvariances.rows(), mvariances.columns());

	for (UINT32 r(0); r<var_scalars.rows(); r+=3)
	{	
		var_scalars.put(r, r, sqrt(pScale));
		var_scalars.put(r+1, r+1, sqrt(lScale));
		var_scalars.put(r+2, r+2, sqrt(hScale));
	}

	// perform the scaling
	ScaleMatrix<T>(*mvariances_mod, mvariances_mod, var_scalars);

	// propagate variances in geographic system to cartesian
	PropagateVariances_GeoCart_Cluster<T>(*mvariances_mod, mvariances_mod,
		mpositions, mrotations, ellipsoid, 
		true, 		// Geographic -> Cartesian
		false);		// don't create a rotation matrix
}

template <class T>
void FormLocaltoCartRotationMatrix(const T& latitude, const T& longitude,
	matrix_2d& mrotations, bool LOCAL_TO_CART=true)
{
	mrotations.redim(3, 3);

	T coslat(cos(latitude));
	T sinlat(sin(latitude));
	T coslon(cos(longitude));
	T sinlon(sin(longitude));

	if (LOCAL_TO_CART)
	{
		// set up Rotation matrix for local to cart
		mrotations.put(0, 0, -sinlon);
		mrotations.put(0, 1, -sinlat*coslon);
		mrotations.put(0, 2, coslat*coslon);
		mrotations.put(1, 0, coslon);
		mrotations.put(1, 1, -sinlat*sinlon);
		mrotations.put(1, 2, coslat*sinlon);
		mrotations.put(2, 0, 0.);
		mrotations.put(2, 1, coslat);
		mrotations.put(2, 2, sinlat);
	}
	else
	{
		// set up Rotation matrix for cart to local, which is 
		// transpose of local to cart!
		mrotations.put(0, 0, -sinlon); 
		mrotations.put(0, 1, coslon); 
		mrotations.put(0, 2, 0.); 
		mrotations.put(1, 0, -sinlat*coslon); 
		mrotations.put(1, 1, -sinlat*sinlon); 
		mrotations.put(1, 2, coslat);
		mrotations.put(2, 0, coslat*coslon); 
		mrotations.put(2, 1, coslat*sinlon); 
		mrotations.put(2, 2, sinlat);
	}
}
	
template <class T>
void FormLocaltoPolarRotationMatrix(const T& azimuth, const T& elevation, const T& distance,
	matrix_2d& mrotations, bool LOCAL_TO_POLAR=true)
{
	mrotations.redim(3, 3);

	T cos_azimuth(cos(azimuth));
	T sin_azimuth(sin(azimuth));
	T cos_elevation(cos(elevation));
	T sin_elevation(sin(elevation));
	
	if (LOCAL_TO_POLAR)
	{
		// set up Jacobian matrix for local to polar
		mrotations.put(0, 0, cos_azimuth/distance); 
		mrotations.put(0, 1, -sin_azimuth/distance); 
		mrotations.put(0, 2, 0.); 
		mrotations.put(1, 0, -sin_azimuth*sin_elevation/distance); 
		mrotations.put(1, 1, -cos_azimuth*sin_elevation/distance); 
		mrotations.put(1, 2, cos_elevation/distance);
		mrotations.put(2, 0, sin_azimuth*cos_elevation); 
		mrotations.put(2, 1, cos_azimuth*cos_elevation); 
		mrotations.put(2, 2, sin_elevation);
	}
	else
	{
		// Set up Jacobian matrix for polar to local, which is 
		// transpose of local to polar!
		mrotations.put(0, 0, cos_azimuth/distance);
		mrotations.put(0, 1, -sin_azimuth*sin_elevation/distance);
		mrotations.put(0, 2, sin_azimuth*cos_elevation);
		mrotations.put(1, 0, -sin_azimuth/distance);
		mrotations.put(1, 1, -cos_azimuth*sin_elevation/distance);
		mrotations.put(1, 2, cos_azimuth*cos_elevation);
		mrotations.put(2, 0, 0.); 
		mrotations.put(2, 1, cos_elevation/distance); 
		mrotations.put(2, 2, sin_elevation);
	}
}

template <class T>
void PropagateVariances_CartLocal_Diagonal(const matrix_2d& mvariances, matrix_2d& mvariances_mod, 
								  const T& latitude, const T& longitude,
								  matrix_2d& mrotations, bool CALCULATE_ROTATIONS=false)
{
	if (CALCULATE_ROTATIONS)
		FormLocaltoCartRotationMatrix<T>(latitude, longitude, mrotations);

	//matrix_2d mrotations_T(mrotations.rows(), mrotations.columns());
	matrix_2d rtv(mrotations.rows(), mrotations.columns());

	//mrotations_T.transpose(mrotations);
	
	mvariances_mod.redim(3, 3);
	
	// RtV
	//rtv.multiply(mrotations_T, mvariances);
	rtv.multiply(mrotations, "T", mvariances, "N");
	
	UINT32 row, col, i;
	
	// RtVR
	for (row=0; row<3; ++row) {
		for (col=0; col<3; ++col) {
			mvariances_mod.put(row, col, 0.0);
			// diagonals only
			if (row != col)
				continue;
			for (i=0; i<3; ++i)
				mvariances_mod.elementadd(row, col, rtv.get(row, i) * mrotations.get(i, col));
		}
	}
}


template <class T>
void PropagateVariances_LocalPolar_Diagonal(const matrix_2d& mvariances, matrix_2d& mvariances_mod, 
								  const T& azimuth, const T& elevation, const T& distance,
								  matrix_2d& mrotations, bool CALCULATE_ROTATIONS=false)
{
	if (CALCULATE_ROTATIONS)
		FormLocaltoPolarRotationMatrix<T>(azimuth, elevation, distance, mrotations);

	matrix_2d mrotations_T(mrotations.rows(), mrotations.columns());
	matrix_2d rtv(mrotations.rows(), mrotations.columns());

	mrotations_T.transpose(mrotations);
	
	mvariances_mod.redim(3, 3);
	
	// RtV
	//rtv.multiply(mrotations, mvariances);
	rtv.multiply(mrotations, "N", mvariances, "N");

	UINT32 row, col, i;
	
	// RtVR
	for (row=0; row<3; ++row) {
		for (col=0; col<3; ++col) {
			mvariances_mod.put(row, col, 0.0);
			// diagonals only
			if (row != col)
				continue;
			for (i=0; i<3; ++i)
				mvariances_mod.elementadd(row, col, rtv.get(row, i) * mrotations_T.get(i, col));
		}
	}
}


template <class T>
void PropagateVariances_LocalCart(const matrix_2d& mvariances, matrix_2d& mvariances_mod, 
								  const T& latitude, const T& longitude, bool LOCAL_TO_CART,
								  matrix_2d& mrotations, bool CALCULATE_ROTATIONS=false)
{
	if (CALCULATE_ROTATIONS)
		FormLocaltoCartRotationMatrix<T>(latitude, longitude, mrotations);

	// form transpose, from either passed in matrix or newly formed matrix
	matrix_2d mrotations_T(mrotations.rows(), mrotations.columns());
	mrotations_T.transpose(mrotations);

	// the rotation matrix is in the direction local to cart (forward)
	if (LOCAL_TO_CART)
	{
		// Vc = R * Vl * RT
		//matrix_2d mV(mvariances_mod.multiply(mrotations, mvariances));
		matrix_2d mV(mvariances_mod.multiply(mrotations, "N", mvariances, "N"));
		//mvariances_mod.multiply_square(mV, mrotations_T);
		mvariances_mod.multiply(mV, "N", mrotations, "T");
	}
	else
	{
		// Vc = R-1 * Vc * [R-1]T
		//    = RT * Vc * R (since R is orthogonal)
		//matrix_2d mV(mvariances_mod.multiply(mrotations_T, mvariances));
		matrix_2d mV(mvariances_mod.multiply(mrotations, "T", mvariances, "N"));
		//mvariances_mod.multiply_square(mV, mrotations);
		mvariances_mod.multiply(mV, "N", mrotations, "N");
	}
}

template <class T>
void PropagateVariances_LocalCart(const matrix_2d& mvariances, matrix_2d& mvariances_mod, 
								  const T& latitude, const T& longitude, bool LOCAL_TO_CART)
{
	matrix_2d mrotations;

	PropagateVariances_LocalCart<T>(mvariances, mvariances_mod, 
		latitude, longitude, LOCAL_TO_CART,
		mrotations, true);	// calculate rotations
}

template <class T>
void Rotate_LocalCart(const matrix_2d mvector, matrix_2d* mvector_mod, 
								  const T& latitude, const T& longitude)
{
	matrix_2d mrotations(3, 3);
	FormLocaltoCartRotationMatrix<T>(latitude, longitude, mrotations);

	mvector_mod->redim(3, 1);
	//mvector_mod->multiply(mrotations, mvector);
	mvector_mod->multiply(mrotations, "N", mvector, "N");
}

template <class T>
void Rotate_CartLocal(const matrix_2d mvector, matrix_2d* mvector_mod, 
								  const T& latitude, const T& longitude)
{
	// Helps
	T sin_lat(sin(latitude));
	T cos_lat(cos(latitude));
	T sin_long(sin(longitude));
	T cos_long(cos(longitude));

	mvector_mod->redim(3, 1);
	
	mvector_mod->put(0, 0, -sin_long * mvector.get(0,0) + cos_long * mvector.get(1,0));
	mvector_mod->put(1, 0, -sin_lat * cos_long * mvector.get(0,0) -
		sin_lat * sin_long * mvector.get(1,0) +
		cos_lat * mvector.get(2,0));
	mvector_mod->put(2, 0, cos_lat * cos_long * mvector.get(0,0) +
		cos_lat * sin_long * mvector.get(1,0) +
		sin_lat * mvector.get(2,0));
}


template <class T>
void Rotate_LocalPolar(const matrix_2d mvector, matrix_2d* mvector_mod, 
								  const T& azimuth, const T& elevation, const T& distance)
{
	// Helps
	T sin_azimuth(sin(azimuth));
	T cos_azimuth(cos(azimuth));
	T sin_elevation(sin(elevation));
	T cos_elevation(cos(elevation));
	
	mvector_mod->redim(3, 1);
	
	mvector_mod->put(0, 0, -cos_azimuth * mvector.get(0,0) - sin_azimuth * mvector.get(1,0));
	mvector_mod->put(1, 0, -sin_elevation * sin_azimuth * mvector.get(0,0) -
		sin_elevation * cos_azimuth * mvector.get(1,0) +
		cos_elevation * mvector.get(2,0));
	mvector_mod->put(2, 0, cos_elevation * sin_azimuth * mvector.get(0,0) +
		cos_elevation * cos_azimuth * mvector.get(1,0) +
		sin_elevation * mvector.get(2,0));
}


template <class T>
// rotx is x rotation in radians
// roty is y rotation in radians
// rotz is z rotation in radians
void FormHelmertRotationMatrix(const T& rotx, const T& roty, const T& rotz, matrix_2d& mrotations, bool RIGOROUS=false)
{
	mrotations.redim(3, 3);

	// Which rotation matrix is required?
	// Convert to seconds for the test
	if (RIGOROUS)
	{
		// rigorous formula for large (> 10 seconds) rotations
		mrotations.put(0, 0, cos(roty) * cos(rotz));
		mrotations.put(0, 1, cos(roty) * sin(rotz));
		mrotations.put(0, 2, -sin(roty));
		mrotations.put(1, 0, (sin(rotx) * sin(roty) * cos(rotz)) - (cos(rotx) * sin(rotz)));
		mrotations.put(1, 1, (sin(rotx) * sin(roty) * sin(rotz)) + (cos(rotx) * cos(rotz)));
		mrotations.put(1, 2, sin(rotx) * cos(roty));
		mrotations.put(2, 0, (cos(rotx) * sin(roty) * cos(rotz)) + (sin(rotx) * sin(rotz)));
		mrotations.put(2, 1, (cos(rotx) * sin(roty) * sin(rotz)) - (sin(rotx) * cos(rotz)));
		mrotations.put(2, 2, cos(rotx) * cos(roty));
	}
	else
	{
		mrotations.put(0, 0, 1.);
		mrotations.put(0, 1, rotz);
		mrotations.put(0, 2, -roty);
		mrotations.put(1, 0, -rotz);
		mrotations.put(1, 1, 1.);
		mrotations.put(1, 2, rotx);
		mrotations.put(2, 0, roty);
		mrotations.put(2, 1, -rotx);
		mrotations.put(2, 2, 1.);
	}
}
	

template <class T>
void ReduceParameters(const T* parameters, T* reduced_parameters, const T& elapsedTime, bool DYNAMIC=true)
{
	// translations (reduce to metres)
	reduced_parameters[0] = parameters[0] / 1000.;
	reduced_parameters[1] = parameters[1] / 1000.;
	reduced_parameters[2] = parameters[2] / 1000.;
	// scale
	reduced_parameters[3] = parameters[3] / 1E9;
	// rotations
	reduced_parameters[4] = parameters[4];
	reduced_parameters[5] = parameters[5];
	reduced_parameters[6] = parameters[6];

	if (DYNAMIC)
	{
		// apply rates to translations
		reduced_parameters[0] += parameters[7] / 1000. * elapsedTime;
		reduced_parameters[1] += parameters[8] / 1000. * elapsedTime;
		reduced_parameters[2] += parameters[9] / 1000. * elapsedTime;
		// apply rate to scale
		reduced_parameters[3] += parameters[10] / 1E9 * elapsedTime;
		// apply rates to rotations
		reduced_parameters[4] += parameters[11] * elapsedTime;
		reduced_parameters[5] += parameters[12] * elapsedTime;
		reduced_parameters[6] += parameters[13] * elapsedTime;
	}

	// reduce rotations from milli-arc-seconds to radians
	reduced_parameters[4] = SecondstoRadians(reduced_parameters[4]) / 1000.;
	reduced_parameters[5] = SecondstoRadians(reduced_parameters[5]) / 1000.;
	reduced_parameters[6] = SecondstoRadians(reduced_parameters[6]) / 1000.; 
}

template <class T>
// No check or safe guard in place to test if mcoordinates_mod is a reference to mcoordinates
void TransformCartesian(const matrix_2d& mcoordinates, matrix_2d& mcoordinates_mod, 
	const matrix_2d& parameters, const matrix_2d& mrotations)
{
	// Add rotation and scale contributions
	for (UINT16 i(0), j; i<3; ++i)
	{
		// Initialise 'to datum' matrix
		mcoordinates_mod.put(i, 0, 0.0);

		for (j=0; j<3; ++j)		// For each column in the row
			mcoordinates_mod.elementadd(i, 0, mrotations.get(i, j) * mcoordinates.get(j, 0));	// Sum partial products

		// Apply scale
		mcoordinates_mod.elementmultiply(i, 0, 1.0 + parameters.get(3, 0));

		// Add translations
		mcoordinates_mod.elementadd(i, 0, parameters.get(i, 0));
	}
}

template <class T>
void Transform_7parameter(const matrix_2d& mcoordinates, matrix_2d& mcoordinates_mod, const T parameters[])
{
	mcoordinates_mod.redim(3, 1);

	bool RIGOROUS(false);
	if (parameters[4] > 10. || parameters[5] > 10. || parameters[6] > 10.)
		RIGOROUS = true;

	// Form rotation matrix
	matrix_2d mrotations(3, 3);
	FormHelmertRotationMatrix<T>(parameters[4], parameters[5], parameters[6], mrotations, RIGOROUS);
	
	// Put translations and scale into reducedParameters
	matrix_2d mtrans_scale(4, 1, parameters, 4);

	// Transform
	TransformCartesian<T>(mcoordinates, mcoordinates_mod, 
		mtrans_scale, mrotations);

}
	

template <typename T>
void PositionalUncertainty(const T& semimajor, const T& semiminor, const T& sdHt,
						   T& hzPosU_Radius, T& vtPosU_Radius)
{
	hzPosU_Radius = vtPosU_Radius = -1.;
	if (semimajor < 0.0)
		return;
	if (semiminor < 0.0)
		return;

	// Horizontal
	T c(semiminor / semimajor);
	T K(HPOS_UNCERT_Q0 + (HPOS_UNCERT_Q1 * c) + (HPOS_UNCERT_Q2 * (c * c)) + (HPOS_UNCERT_Q3 * (c * c * c)));
	T R(semimajor * K);
		
	hzPosU_Radius = R;
		
	// Vertical
	vtPosU_Radius = sdHt * 1.96;
}

template <typename T>
T PedalVariance(const matrix_2d& mvariance, const T& direction)
{
	T cos_theta(cos(direction));
	T sin_theta(sin(direction));
	return mvariance.get(0, 0) * cos_theta * cos_theta +
		mvariance.get(1, 1) * sin_theta * sin_theta +
		2. * mvariance.get(0, 1) * cos_theta * sin_theta;
}
	

template <typename T>
void ErrorEllipseParameters(const matrix_2d& mvariance, T& semimajor, T& semiminor, T& azimuth)
{
	semimajor = semiminor = azimuth = -1.;

	if (mvariance.rows() < 2)
		return;
	if (mvariance.columns() < 2)
		return;

	T e2(mvariance.get(0, 0));
	T n2(mvariance.get(1, 1));
	T en(mvariance.get(0, 1));
	T e2_plus_n2(e2 + n2);
	T e2_minus_n2(e2 - n2);
	T n2_minus_e2(n2 - e2);
	
	T W((e2_minus_n2 * e2_minus_n2) + (4. * en * en));
	
	if (W < 0.0)
	{
		if (fabs(W) > PRECISION_1E15)
			return;								// temp term cannot be negative!!!
		else
			W = 0.0;
	}

	T a_sqd(0.5 * (e2_plus_n2 + sqrt(W)));		// semi-major
	T b_sqd(0.5 * (e2_plus_n2 - sqrt(W)));		// semi-minor
	
	if (a_sqd < 0.0)
		return;									// lamda2 term cannot be negative!!!
	if (b_sqd < 0.0)
		return;									// lamda2 term cannot be negative!!!

	semimajor = sqrt(a_sqd);
	semiminor = sqrt(b_sqd);

	// Compute the azimuth of the semi-major axis
	if (fabs(e2 - n2) < PRECISION_1E25)
	{
		if (en < PRECISION_1E25)
			azimuth = 0.;			// ellipse is a circle
		else
			azimuth = PI / 4.;		// azimuth = 45
	}
	else
	{
		T x(en + en);
		azimuth = 0.5 * atan_2(x, n2_minus_e2);
	}
}

#endif /* DNATEMPLATEMATRIXFUNCS_H_ */
