//============================================================================
// Name         : dnatemplatestnmsrfuncs.hpp
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
// Description  : Common calculation functions using predefined station and
//				  measurement types
//============================================================================

#ifndef DNATEMPLATESTNMSRFUNCS_H_
#define DNATEMPLATESTNMSRFUNCS_H_

#if defined(_MSC_VER)
	#if defined(LIST_INCLUDES_ON_BUILD) 
		#pragma message("  " __FILE__) 
	#endif
#endif

/// \cond
#include <stdio.h>

#include <iostream>
#include <algorithm>
#include <vector>
#include <string>
#include <cctype>

#include <memory>

/// \endcond

#include <include/config/dnaexports.hpp>
#include <include/functions/dnatemplatefuncs.hpp>
#include <include/functions/dnastrutils.hpp>
#include <include/config/dnatypes.hpp>
#include <include/measurement_types/dnameasurement_types.hpp>


using namespace dynadjust::measurements;

template <typename T>
// Get all the statons asociated with a measurement
void GetGXMsrStations(std::vector<CDnaGpsBaseline>* vgpsBsls, std::vector<T>& msrStations)
{
	msrStations.clear();
			
	// Step 1. Set ignore flag for all baselines to true
	for_each(vgpsBsls->begin(), vgpsBsls->end(), 
		[&msrStations](CDnaGpsBaseline& bsl) {
			msrStations.push_back(bsl.GetFirst());
			msrStations.push_back(bsl.GetTarget());
	});

	// Strip duplicate entries
	strip_duplicates(msrStations);
}


template <typename T>
// Get all the statons asociated with a measurement
void GetMsrStations(const vmsr_t& binaryMsrs, const T& bmsIndex, std::vector<T>& msrStations)
{
	msrStations.clear();
	T clusterID(binaryMsrs.at(bmsIndex).clusterID);

	it_vmsr_t_const _it_msr_start = binaryMsrs.begin() + bmsIndex;
	it_vmsr_t_const _it_msr_prev = _it_msr_start;

	bool bCluster(false);

	// I don't think this is necessary for X measruements, since
	// bmsIndex points to the first element in a cluster.  Perhaps keep
	// this here to safe guard against cases where bmsIndex is not the
	// first???
	switch (binaryMsrs.at(bmsIndex).measType)
	{
	case 'D':
	case 'X':
	case 'Y':
		// move to the beginning of the cluster
		while (_it_msr_start != binaryMsrs.begin())
		{
			if ((--_it_msr_prev)->clusterID != clusterID)
				break;
			_it_msr_start--;
		}
		bCluster = true;
		break;
	}

	while (_it_msr_start->clusterID == clusterID)		// Will be true for non-cluster measurements
	{	
		if (_it_msr_start->measType != 'D')
		{
			if (_it_msr_start->measStart != xMeas)
			{
				if (++_it_msr_start == binaryMsrs.end())
					break;
				continue;
			}
		}		
		
		// add the relevant stations
		// Station 1 is in all measurement types
		msrStations.push_back(_it_msr_start->station1);
		// Station 2 is in two- and three-station measurements
		if (MsrTally::Stations(_it_msr_start->measType) >= TWO_STATION)
			msrStations.push_back(_it_msr_start->station2);
		// Station 3 is only in three-station measurements
		if (MsrTally::Stations(_it_msr_start->measType) == THREE_STATION)
			msrStations.push_back(_it_msr_start->station3);
	
		if (!bCluster)
			break;

		if (++_it_msr_start == binaryMsrs.end())
			break;
	}

	// Strip duplicate entries
	strip_duplicates(msrStations);
}

template <typename T>
// Get all the binary measurement indices involved in a measurement
// On return, msrIndices will contain only one index.  For clusters
// and directions, msrIndices will be more than one.
void GetMsrIndices(const vmsr_t& binaryMsrs, const T& bmsIndex, std::vector<T>& msrIndices)
{
	msrIndices.clear();

	// First, sift out conventional measurement types and return immediately.
	switch (binaryMsrs.at(bmsIndex).measType)
	{
	case 'D':
	case 'X':
	case 'Y':
		// Clusters dealt with after this switch statement
		break;
	default:
		// basic measurement type for which there is only one index, so
		// add the index and return
		msrIndices.push_back(bmsIndex);
		return;
	}

	// Ok, now deal with cluster measurements
	UINT32 clusterID(binaryMsrs.at(bmsIndex).clusterID);

	it_vmsr_t_const _it_msr_start(binaryMsrs.begin() + bmsIndex);
	it_vmsr_t_const _it_msr_prev(_it_msr_start);
	
	// move to the beginning of the cluster
	while (_it_msr_start != binaryMsrs.begin())
	{
		if ((--_it_msr_prev)->clusterID != clusterID)
			break;
		--_it_msr_start;
	}


	while (_it_msr_start->clusterID == clusterID && _it_msr_start != binaryMsrs.end())
	{
		if (_it_msr_start->ignore)
		{
			if (++_it_msr_start == binaryMsrs.end())
				break;
			continue;
		}

		if (_it_msr_start->measType != 'D')
		{
			if (_it_msr_start->measStart > xMeas)
			{
				if (++_it_msr_start == binaryMsrs.end())
					break;
				continue;
			}
		}		

		msrIndices.push_back(static_cast<T>(std::distance(binaryMsrs.begin(), _it_msr_start)));
	
		if (++_it_msr_start == binaryMsrs.end())
			break;
	}

	// Sort and strip duplicate entries
	strip_duplicates(msrIndices);
}

template <typename T>
// Get the first binary measurement index for a measurement
T GetFirstMsrIndex(const vmsr_t& binaryMsrs, const T& bmsIndex)
{
	// First, sift out conventional measurement types and return immediately.
	switch (binaryMsrs.at(bmsIndex).measType)
	{
	case 'D':
	case 'X':
	case 'Y':
		// Clusters dealt with after this switch statement
		break;
	default:
		// basic measurement type for which there is only one index, so
		// return bmsIndex
		return bmsIndex;
	}

	// Ok, now deal with cluster measurements
	UINT32 clusterID(binaryMsrs.at(bmsIndex).clusterID);

	it_vmsr_t_const _it_msr_start(binaryMsrs.begin() + bmsIndex);
	it_vmsr_t_const _it_msr_prev(_it_msr_start);
	
	// move to the beginning of the cluster
	while (_it_msr_start != binaryMsrs.begin())
	{
		if ((--_it_msr_prev)->clusterID != clusterID)
			break;
		--_it_msr_start;
	}

	// return the index of the measurement's first record in the binary measurement file 
	return static_cast<UINT32>(std::distance(binaryMsrs.begin(), _it_msr_start));
}

//template <typename T, typename Iter>
//void RenameStationsMsr(T* msr, Iter begin, Iter end)
//{
//	for (_it_string_vstring_pair it = begin;
//		it != end;
//		++it)
//	{
//
//	it_string_pair it;
//	it = binary_search_index_pair(begin, end, msr->GetFirst());
//	
//	if (it != end)
//		msr->SetFirst(it->second);
//		
//	// Does this measurement have a second station?
//	if (msr->m_MSmeasurementStations == ONE_STATION)
//		return;
//
//	// Okay, at least two stations
//	it = binary_search_index_pair(begin, end, msr->GetTarget());
//	if (it != end)
//		msr->SetTarget(it->second);
//			
//	// Does this measurement have a third station?
//	if (msr->m_MSmeasurementStations == TWO_STATION)
//		return;
//			
//	// Okay, three stations
//	it = binary_search_index_pair(begin, end, msr->GetTarget2());
//	if (it != end)
//		msr->SetTarget2(it->second);
//}
	

template <typename T, typename Iter>
void RenameStationsMsr(T* msr, Iter begin, Iter end)
{
	_it_string_vstring_pair it;
	Iter original = begin;
	
	// Search the aliases for the first station name
	// Emulate for_each behaviour
	while (begin != end) 
	{		
		if (binary_search(begin->second.begin(), begin->second.end(), msr->GetFirst()))
		{
			msr->SetFirst(begin->first);
			break;
		}
		++begin;
	}

	// Is this measurement a one-station measurement?
	if (msr->m_MSmeasurementStations == ONE_STATION)
		return;

	// Search the aliases for the second station name
	// Emulate for_each behaviour
	begin = original;
	while (begin != end) 
	{		
		if (binary_search(begin->second.begin(), begin->second.end(), msr->GetTarget()))
		{
			msr->SetTarget(begin->first);
			break;
		}
		++begin;
	}

	// Is this measurement a two-station measurement?
	if (msr->m_MSmeasurementStations == TWO_STATION)
		return;

	// Search the aliases for the third station name
	// Emulate for_each behaviour
	begin = original;
	while (begin != end) 
	{		
		if (binary_search(begin->second.begin(), begin->second.end(), msr->GetTarget2()))
		{
			msr->SetTarget2(begin->first);
			break;
		}
		++begin;
	}
}

template <typename T, typename Iter>
void IgnoreGXMeasurements(T* msr, Iter begin, Iter end)
{
	
	// Search the aliases for the first station name
	// Emulate for_each behaviour	
	if (binary_search(begin, end, msr->GetClusterID()))
	{
		msr->SetIgnore(true);	
		for_each(msr->GetBaselines_ptr()->begin(),
			msr->GetBaselines_ptr()->end(),
			[] (CDnaGpsBaseline& bsl) {
				bsl.SetIgnore(true);
		});
	}
}


template <typename T, typename msriterator>
// Copy the cluster measurement 
void CopyClusterMsr(T& cluster, const msriterator _it_msr, T& clusterCopy)
{
	clusterCopy.clear();

	// First, sift out conventional measurement types and return immediately.
	switch (_it_msr->measType)
	{
	case 'D':
	case 'X':
	case 'Y':
		// Clusters dealt with after this switch statement
		break;
	default:
		// basic measurement type for which there is only one element
		clusterCopy.push_back(*_it_msr);
		return;
	}

	// Ok, now deal with cluster measurements
	UINT32 clusterID(_it_msr->clusterID);

	msriterator _it_msr_start(_it_msr);
	msriterator _it_msr_temp(_it_msr_start);
	
	while (_it_msr_start != cluster.begin())
	{
		if ((--_it_msr_temp)->clusterID != clusterID)
			break;
		--_it_msr_start;
	}

	msriterator _it_msr_last(_it_msr);

	while (_it_msr_last != cluster.end())
	{
        if ((++_it_msr_temp) == cluster.end()) 
		{
            ++_it_msr_last;
            break;
        }
		if (_it_msr_temp->measType != 'Y')
		{
			++_it_msr_last;
			break;
		}
		++_it_msr_last;
	}

	clusterCopy.insert(clusterCopy.begin(), _it_msr_start, _it_msr_last);
}


// Comparison functions
template <typename T = dnaStnPtr, typename S = std::string>
class CompareStationName {
public:
	bool operator()(const T& left, const T& right) {
		if (left->GetName() == right->GetName())
			return (left->GetfileOrder() < right->GetfileOrder());
		return left->GetName() < right->GetName();
	}
	bool operator()(const T& left, const S& right) {
		return left->GetName() < right;
	}
	bool operator()(const S& left, const T& right) {
		return left < right->GetName();
	}
};

template <typename T = dnaStnPtr, typename S = std::string>
class EqualStationNameSaveDuplicates {
public:
	EqualStationNameSaveDuplicates(std::vector<S>* stns) : _stns(stns) {}

	bool operator()(const T& left, const T& right) {
		if (equals(left->GetName(), right->GetName()))
			_stns->push_back(right->GetName());

		return (left->GetName() == right->GetName());
	}

	std::vector<S>* _stns;
};


template <typename T = dnaStnPtr>
class EqualStationName {
public:
	bool operator()(const T& left, const T& right) {
		return (left->GetName() == right->GetName());
	}
};


template <typename T = dnaStnPtr, typename S = std::string>
class EqualStationName_CaseInsensitive {
public:
	EqualStationName_CaseInsensitive(std::vector<S>* stns) : _stns(stns) {}

	bool operator()(const T& left, const T& right) {
		if (iequals(left->GetName(), right->GetName()))
		{
			_stns->push_back(right->GetName());
			return true;
		}
		return false;
	}
	std::vector<S>* _stns;
};


template <typename T = dnaStnPtr, typename S = std::string>
class TestEqualStationName {
public:
	TestEqualStationName(std::vector<S>* usedStns, std::vector<S>* unusedStns) 
		: _usedStns(usedStns)
		, _unusedStns(unusedStns) {}

	bool operator()(const T& stn) {
		// Is stn in the list of stations to keep?  If so, return false
		if (binary_search(_usedStns->begin(), _usedStns->end(), stn->GetName()))
			return true;
		
		// Station is not in the list, so add it to the unused stations
		_unusedStns->push_back(stn->GetName());
		return false;		
	}
	std::vector<S>* _usedStns;
	std::vector<S>* _unusedStns;
};


template <typename T = dnaStnPtr, typename S = std::string>
class TestNotEqualStationName {
public:
	TestNotEqualStationName(std::vector<S>* usedStns, std::vector<S>* unusedStns) 
		: _usedStns(usedStns)
		, _unusedStns(unusedStns) {}

	bool operator()(const T& stn) {
		// Is stn in the list of stations to keep?  If so, return false
		if (binary_search(_usedStns->begin(), _usedStns->end(), stn->GetName()))
			return false;
		
		// Station is not in the list, so add it to the unused stations
		_unusedStns->push_back(stn->GetName());
		return true;		
	}
	std::vector<S>* _usedStns;
	std::vector<S>* _unusedStns;
};


// uses low accuracy spherical formula
template <typename T = dnaStnPtr, typename S = stringstring_doubledouble_pair, typename U = double>
class NearbyStation_LowAcc {
public:
	NearbyStation_LowAcc (const U& tolerance, std::vector<S>* stns) 
		: _tolerance(tolerance), _stns(stns), _dist(0.) {}

	bool operator()(const T& left, const T& right) {
		if ((_dist = GreatCircleDistance(
			left->GetLatitude(), 
			left->GetLongitude(), 
			right->GetLatitude(), 
			right->GetLongitude())) < _tolerance)
		{
			_stns->push_back(stringstring_doubledouble_pair(
				string_string_pair(left->GetName(), right->GetName()),
				doubledouble_pair(_dist, right->GetHeight()-left->GetHeight())));
			return true;
		}
		return false;
	}

	U			_tolerance;
	std::vector<S>*	_stns;
	double		_dist;
};


// uses high accuracy spherical formula
template <typename T = dnaStnPtr, typename U = double, typename S = stringstring_doubledouble_pair, typename E = CDnaEllipsoid>
class NearbyStation_HighAcc {
public:
	NearbyStation_HighAcc(const U& tolerance, std::vector<S>* stns, const E& ellipsoid) 
		: _tolerance(tolerance), _stns(stns), _dAzimuth(0.), _dist(0.), _ellipsoid(ellipsoid) {}

	bool operator()(const T& left, const T& right) {
		if ((_dist = RobbinsReverse(
			left->GetLatitude(), 
			left->GetLongitude(), 
			right->GetLatitude(), 
			right->GetLongitude(), 
			&_dAzimuth, &_ellipsoid)) < _tolerance)
		{
			_stns->push_back(stringstring_doubledouble_pair(
				string_string_pair(left->GetName(), right->GetName()),
				doubledouble_pair(_dist, right->GetHeight()-left->GetHeight())));
			return true;
		}
		return false;
	}

	U			_tolerance;
	std::vector<S>*	_stns;
	U			_dAzimuth;
	double		_dist;
	E			_ellipsoid;		// GDA by default
};


// T = double/float
template<typename T = double>
class FindStnsWithinBoundingBox{
public:
	FindStnsWithinBoundingBox(const T& upperLat, const T& upperLon, const T& lowerLat, const T& lowerLon, pvstring stns)
		: _upperLat(upperLat), _upperLon(upperLon), _lowerLat(lowerLat), _lowerLon(lowerLon)
		, _stns(stns) {}

	bool operator()(dnaStnPtr s) {
		if (s->GetXAxis() > _upperLat)
		{
			_stns->push_back(s->GetName());
			return true;
		}
		if (s->GetXAxis() < _lowerLat)
		{
			_stns->push_back(s->GetName());
			return true;
		}
		if (s->GetYAxis() > _upperLon)
		{
			_stns->push_back(s->GetName());
			return true;
		}
		if (s->GetYAxis() < _lowerLon)
		{
			_stns->push_back(s->GetName());
			return true;		
		}
		return false;
	}

private:
	T _upperLat;
	T _upperLon;
	T _lowerLat;
	T _lowerLon;
	pvstring _stns;
};

// requires stns to be sorted
template<typename U = std::string>
class FindMsrsConnectedToStns_GX{
public:
	FindMsrsConnectedToStns_GX(const pvstring stns)
		: _stns(stns) {}

	bool operator()(const CDnaGpsBaseline& m) {
		if (binary_search(_stns->begin(), _stns->end(), m.GetFirst()))
			return true;
		if (binary_search(_stns->begin(), _stns->end(), m.GetTarget()))
			return true;
		return false;
	}
private:
	pvstring _stns;
};

// requires stns to be sorted
template<typename U = std::string>
class FindMsrsConnectedToStns_Y{
public:
	FindMsrsConnectedToStns_Y(const pvstring stns)
		: _stns(stns) {}

	bool operator()(const CDnaGpsPoint& m) {
		if (binary_search(_stns->begin(), _stns->end(), m.GetFirst()))
			return true;
		return false;
	}	
private:
	pvstring _stns;
};

// requires stns to be sorted
template<typename U = std::string>
class FindMsrsConnectedToStns_D{
public:
	FindMsrsConnectedToStns_D(const pvstring stns)
		: _stns(stns) {}

	bool operator()(const CDnaDirection& m) {
		if (binary_search(_stns->begin(), _stns->end(), m.GetTarget()))
			return true;
		return false;
	}	
private:
	pvstring _stns;
};

// T = double/float
// requires stns to be sorted
template<typename T = pvstring>
class FindMsrsConnectedToStns{
public:
	FindMsrsConnectedToStns(const T stns)
		: _stns(stns) {}

	bool operator()(dnaMsrPtr m) {

		FindMsrsConnectedToStns_GX<std::string> gpsbslFunc(_stns);
		FindMsrsConnectedToStns_Y<std::string> gpspntFunc(_stns);
		FindMsrsConnectedToStns_D<std::string> dirnFunc(_stns);
		switch (m->GetTypeC())
		{
		
		// vector-type measurements
		case 'G':
		case 'X':			
			return find_if(m->GetBaselines_ptr()->begin(), m->GetBaselines_ptr()->end(), gpsbslFunc) != m->GetBaselines_ptr()->end();
		case 'Y':			
			return find_if(m->GetPoints_ptr()->begin(), m->GetPoints_ptr()->end(), gpspntFunc) != m->GetPoints_ptr()->end();
		}

		if (binary_search(_stns->begin(), _stns->end(), m->GetFirst()))
			return true;

		//
		// single station measurements
		//
		switch (m->GetTypeC())
		{
		case 'H':	// Orthometric height
		case 'I':	// Astronomic latitude
		case 'J':	// Astronomic longitude
		case 'P':	// Geodetic latitude
		case 'Q':	// Geodetic longitude
		case 'R':	// Ellipsoidal height
			return false;
		}

		//
		// dual (or more) station measurements
		//
		if (binary_search(_stns->begin(), _stns->end(), m->GetTarget()))
			return true;

		switch (m->GetTypeC())
		{
		case 'D':	// Direction set
			return find_if(
				m->GetDirections_ptr()->begin(), 
				m->GetDirections_ptr()->end(), 
				dirnFunc) != m->GetDirections_ptr()->end();
		case 'C':	// Chord dist
		case 'E':	// Ellipsoid arc
		case 'M':	// MSL arc
		case 'S':	// Slope distance
		case 'L':	// Level difference
		case 'B':	// Geodetic azimuth
		case 'K':	// Astronomic azimuth
		case 'V':	// Zenith distance
		case 'Z':	// Vertical angle
			return false;
		}
		
		//
		// triple station measurements
		//
		if (binary_search(_stns->begin(), _stns->end(), m->GetTarget2()))
			return true;

		//switch (m->GetTypeC())
		//{
		//case 'A':	// Horizontal angles
		//default:
		//	return false;
		//}
		
		return false;
	}

private:
	T _stns;
};
	

template<typename T>
void ResetMeasurementPtr(dnaMsrPtr* msrPtr, const T& cType)
{
	switch (cType)
	{
	case 'A': // Horizontal angle
		msrPtr->reset(new CDnaAngle);
		break;
	case 'B': // Geodetic azimuth
	case 'K': // Astronomic azimuth
		msrPtr->reset(new CDnaAzimuth);
		break;
	case 'C': // Chord dist
	case 'E': // Ellipsoid arc
	case 'M': // MSL arc
	case 'S': // Slope distance
		msrPtr->reset(new CDnaDistance);
		break;
	case 'D': // Direction set
		msrPtr->reset(new CDnaDirectionSet);
		break;
	case 'G': // GPS Baseline (treat as single-baseline cluster)
	case 'X': // GPS Baseline cluster
		msrPtr->reset(new CDnaGpsBaselineCluster);
		break;
	case 'H': // Orthometric height
		msrPtr->reset(new CDnaHeight);
		break;
	case 'I': // Astronomic latitude
	case 'J': // Astronomic longitude
	case 'P': // Geodetic latitude
	case 'Q': // Geodetic longitude
		msrPtr->reset(new CDnaCoordinate);
		break;
	case 'L': // Level difference
		msrPtr->reset(new CDnaHeightDifference);
		break;
	case 'R': // Ellipsoidal height
		msrPtr->reset(new CDnaHeight);
		break;
	case 'V': // Zenith distance
	case 'Z': // Vertical angle
		msrPtr->reset(new CDnaDirection);
		break;
	case 'Y': // GPS point cluster
		msrPtr->reset(new CDnaGpsPointCluster);
		break;
	}
}
	
// T = double/float
// requires stns to be sorted
template<typename T = CDnaDirection>
class SortDirectionsObsClockwise{
public:
	SortDirectionsObsClockwise() {};

	bool operator()(const T& left, const T& right) {
		if (left.GetValue() == right.GetValue())
			return left.m_strTarget < right.m_strTarget;
		return left.GetValue() < right.GetValue();
	}
};
	
template <typename T = dnaStnPtr>
class CompareLatitude {
public:
	bool operator()(const T& left, const T& right) {
		return left->GetXAxis() < right->GetXAxis();
	}
};


template <typename T = dnaStnPtr>
class CompareLongitude {
public:
	bool operator()(const T& left, const T& right)	{
		return left->GetYAxis() < right->GetYAxis();
	}
};

// used to sort measurements
template <typename M = CDnaMeasurement>
class CompareMsr {
public:
	bool operator()(const std::shared_ptr<M> left, const std::shared_ptr<M> right) {
		if (left->GetIgnore() == right->GetIgnore()) {
			if (iequals(left->GetType(), right->GetType()))
			{
				switch (left->GetTypeC())
				{
				case 'A':	// Horizontal angles
					return *(dynamic_cast<const CDnaAngle*>(&(*left))) < 
						*(dynamic_cast<const CDnaAngle*>(&(*right)));
				
				case 'C':	// Chord dist
				case 'E':	// Ellipsoid arc
				case 'M':	// MSL arc
				case 'S':	// Slope distance
					return *(dynamic_cast<const CDnaDistance*>(&(*left))) < 
						*(dynamic_cast<const CDnaDistance*>(&(*right)));
				
				case 'D':	// Direction set
					return *(dynamic_cast<const CDnaDirectionSet*>(&(*left))) < 
						*(dynamic_cast<const CDnaDirectionSet*>(&(*right)));
				
				case 'B':	// Geodetic azimuth
				case 'K':	// Astronomic azimuth
					return *(dynamic_cast<const CDnaAzimuth*>(&(*left))) < 
						*(dynamic_cast<const CDnaAzimuth*>(&(*right)));
				
				case 'I':	// Astronomic latitude
				case 'J':	// Astronomic longitude
				case 'P':	// Geodetic latitude
				case 'Q':	// Geodetic longitude
					return *(dynamic_cast<const CDnaCoordinate*>(&(*left))) < 
						*(dynamic_cast<const CDnaCoordinate*>(&(*right)));
				
				case 'H':	// Orthometric height
				case 'R':	// Ellipsoidal height
					return *(dynamic_cast<const CDnaHeight*>(&(*left))) < 
						*(dynamic_cast<const CDnaHeight*>(&(*right)));
				
				case 'L':	// Level difference
					return *(dynamic_cast<const CDnaHeightDifference*>(&(*left))) < 
						*(dynamic_cast<const CDnaHeightDifference*>(&(*right)));
				
				case 'V':	// Zenith distance
				case 'Z':	// Vertical angle
					return *(dynamic_cast<const CDnaDirection*>(&(*left))) < 
						*(dynamic_cast<const CDnaDirection*>(&(*right)));
				
				case 'G':	// GPS Baseline (treat as single-baseline cluster)
				case 'X':	// GPS Baseline cluster
					return *(dynamic_cast<const CDnaGpsBaselineCluster*>(&(*left))) < 
						*(dynamic_cast<const CDnaGpsBaselineCluster*>(&(*right)));
				
				case 'Y':	// GPS point cluster				
					return *(dynamic_cast<const CDnaGpsPointCluster*>(&(*left))) < 
						*(dynamic_cast<const CDnaGpsPointCluster*>(&(*right)));
				}
				
				return left < right;				
			}
			else
				return left->GetType() < right->GetType(); }		
		return left->GetIgnore() < right->GetIgnore();
	}
};

	

template <class CharT, class Traits>
std::basic_istream<CharT, Traits>& operator>>(std::basic_istream<CharT, Traits>& is, CAStationList t)
{
	is.read(reinterpret_cast<char *>(t.GetAssocMsrCountPtr()), sizeof(UINT32));
	is.read(reinterpret_cast<char *>(t.GetAMLStnIndexPtr()), sizeof(UINT32));
	is.read(reinterpret_cast<char *>(t.ValidityPtr()), sizeof(UINT16));
	return is;
}

template <class CharT, class Traits>
std::basic_istream<CharT, Traits>& operator>>(std::basic_istream<CharT, Traits>& is, CAStationList* t)
{
	is.read(reinterpret_cast<char *>(t->GetAssocMsrCountPtr()), sizeof(UINT32));
	is.read(reinterpret_cast<char *>(t->GetAMLStnIndexPtr()), sizeof(UINT32));
	is.read(reinterpret_cast<char *>(t->ValidityPtr()), sizeof(UINT16));
	return is;
}

template <class CharT, class TraitsT>
std::basic_ostream<CharT, TraitsT>& operator<<(std::basic_ostream<CharT, TraitsT>& os, CAStationList t)
{
	os.write(reinterpret_cast<char *>(t.GetAssocMsrCountPtr()), sizeof(UINT32));
	os.write(reinterpret_cast<char *>(t.GetAMLStnIndexPtr()), sizeof(UINT32));
	os.write(reinterpret_cast<char *>(t.ValidityPtr()), sizeof(UINT16));
	return os;
}

template <class CharT, class TraitsT>
std::basic_ostream<CharT, TraitsT>& operator<<(std::basic_ostream<CharT, TraitsT>& os, CAStationList* t)
{
	os.write(reinterpret_cast<char *>(t->GetAssocMsrCountPtr()), sizeof(UINT32));
	os.write(reinterpret_cast<char *>(t->GetAMLStnIndexPtr()), sizeof(UINT32));
	os.write(reinterpret_cast<char *>(t->ValidityPtr()), sizeof(UINT16));
	return os;
}


class StationNameIDCompareName {						// class for station name and ID comparison
public: //functions
	// comparison func for sorting
	bool operator()(const string_uint32_pair& lhs, const string_uint32_pair& rhs) const {
		return keyLess(lhs.first, rhs.first);
	}
	// comparison func for lookups
	bool operator()(const string_uint32_pair& lhs, const string_uint32_pair::first_type& k) const {
		return keyLess(lhs.first, k);
	}
	// comparison func for lookups
	bool operator()(const string_uint32_pair::first_type& k, const string_uint32_pair& rhs) const {
		return keyLess(k, rhs.first);
	}
private:
	// the "real" comparison function
	bool keyLess(const string_uint32_pair::first_type& k1, const string_uint32_pair::first_type& k2) const {
		return k1 < k2;
	}
};

class StationNameIDCompareID {						// class for station name and ID comparison
public: //functions
	// comparison func for sorting
	bool operator()(const string_uint32_pair& lhs, const string_uint32_pair& rhs) const {
		return keyLess(lhs.second, rhs.second);
	}
	// comparison func for lookups
	bool operator()(const string_uint32_pair& lhs, const string_uint32_pair::second_type& k) const {
		return keyLess(lhs.second, k);
	}
	// comparison func for lookups
	bool operator()(const string_uint32_pair::second_type& k, const string_uint32_pair& rhs) const {
		return keyLess(k, rhs.second);
	}
private:
	// the "real" comparison function
	bool keyLess(const string_uint32_pair::second_type& k1, const string_uint32_pair::second_type& k2) const {
		return k1 < k2;
	}
};


class ParamStnMsrCompareStn {						// class for station name and ID comparison
public: //functions
	// comparison func for sorting
	bool operator()(const uint32_string_pair& lhs, const uint32_string_pair& rhs) const {
		return keyLess(lhs.first, rhs.first);
	}
	// comparison func for lookups
	bool operator()(const uint32_string_pair& lhs, const uint32_string_pair::first_type& k) const {
		return keyLess(lhs.first, k);
	}
	// comparison func for lookups
	bool operator()(const uint32_string_pair::first_type& k, const uint32_string_pair& rhs) const {
		return keyLess(k, rhs.first);
	}
private:
	// the "real" comparison function
	bool keyLess(const uint32_string_pair::first_type& k1, const uint32_string_pair::first_type& k2) const {
		return k1 < k2;
	}
};



// A = CAStationList, U = UINT32
// Example use:
// CompareMeasCount<CAStationList, UINT32> msrcountCompareFunc(&vAssocStnList_, vAssocStnList_.at(_it_stnmap->second).GetAssocMsrCount());
// std::shared_ptr<UINT32> stnID(new UINT32(_it_stnmap->second));




// M = measurement_t, U = UINT32
// M = measurement_t, U = UINT32
// M = measurement_t, U = UINT32

// S = station_t, U = u32u32_double_pair
template <typename S, typename U = u32u32_double_pair>
class CompareBlockStationMapUnique_byFileOrder {
public:
	CompareBlockStationMapUnique_byFileOrder(std::vector<S>* s)
		:  _s(s) {}
	bool operator()(const U& left, const U& right) {
		return _s->at(left.first.first).fileOrder < _s->at(right.first.first).fileOrder;
	}
private:
	std::vector<S>*	_s;
};

// M = measurement_t, U = UINT32
template<typename M, typename U>
class CompareMeasType_PairFirst
{
public:
	CompareMeasType_PairFirst(std::vector<M>* m)
		:  _m(m) {}
	bool operator()(const std::pair<U, std::pair<U, U> >& lhs, const std::pair<U, std::pair<U, U> >& rhs) {
		if (_m->at(lhs.first).measType == _m->at(rhs.first).measType)
		{
			if (_m->at(lhs.first).station1 == _m->at(rhs.first).station1)
			{
				if (_m->at(lhs.first).station2 == _m->at(rhs.first).station2)
					return _m->at(lhs.first).term1 < _m->at(rhs.first).term1;
				else
					return _m->at(lhs.first).station2 < _m->at(rhs.first).station2;	
			}
			return _m->at(lhs.first).station1 < _m->at(rhs.first).station1;
		}
		return _m->at(lhs.first).measType < _m->at(rhs.first).measType;
	}
private:
	std::vector<M>*	_m;
};


// M = measurement_t, U = UINT32
template<typename M, typename U>
class CompareMeasFromStn_PairFirst
{
public:
	CompareMeasFromStn_PairFirst(std::vector<M>* m)
		:  _m(m) {}
	bool operator()(const std::pair<U, std::pair<U, U> >& lhs, const std::pair<U, std::pair<U, U> >& rhs) {
		if (_m->at(lhs.first).station1 == _m->at(rhs.first).station1)
		{
			if (_m->at(lhs.first).measType == _m->at(rhs.first).measType)
			{
				if (_m->at(lhs.first).station2 == _m->at(rhs.first).station2)
					return _m->at(lhs.first).term1 < _m->at(rhs.first).term1;
				else
					return _m->at(lhs.first).station2 < _m->at(rhs.first).station2;
			}
			return _m->at(lhs.first).measType < _m->at(rhs.first).measType;
		}
		return _m->at(lhs.first).station1 < _m->at(rhs.first).station1;
	}
private:
	std::vector<M>*	_m;
};


// M = measurement_t, U = UINT32
template<typename M, typename U>
class CompareMeasToStn_PairFirst
{
public:
	CompareMeasToStn_PairFirst(std::vector<M>* m)
		:  _m(m) {}
	bool operator()(const std::pair<U, std::pair<U, U> >& lhs, const std::pair<U, std::pair<U, U> >& rhs) {
		if (_m->at(lhs.first).station2 == _m->at(rhs.first).station2)
		{
			if (_m->at(lhs.first).measType == _m->at(rhs.first).measType)
			{
				if (_m->at(lhs.first).station1 == _m->at(rhs.first).station1)
					return _m->at(lhs.first).term1 < _m->at(rhs.first).term1;
				else
					return _m->at(lhs.first).station1 < _m->at(rhs.first).station1;
			}
			return _m->at(lhs.first).measType < _m->at(rhs.first).measType;
		}
		return _m->at(lhs.first).station2 < _m->at(rhs.first).station2;
	}
private:
	std::vector<M>*	_m;
};


// m = measurement_t
template<typename m>
bool isCompoundMeas(const m& msrType)
{
	switch (msrType)
	{
	case 'G':
	case 'X':
	case 'Y':
		return true;
	}
	return false;
}

template<typename m>
bool notCompoundMeas(const m& msrType)
{
	switch (msrType)
	{
	case 'G':
	case 'X':
	case 'Y':
		return false;
	}
	return true;
}

// m = measurement_t
template<typename m>
bool isCompoundMeasAll(const m& msrType)
{
	switch (msrType)
	{
	case 'D':
	case 'G':
	case 'X':
	case 'Y':
		return true;
	}
	return false;
}

template<typename m>
bool notCompoundMeasAll(const m& msrType)
{
	switch (msrType)
	{
	case 'D':
	case 'G':
	case 'X':
	case 'Y':
		return false;
	}
	return true;
}

// M = measurement_t, U = UINT32
template<typename M, typename U>
class CompareMeasValue_PairFirst
{
public:
	CompareMeasValue_PairFirst(std::vector<M>* m)
		:  _m(m) {}
	bool operator()(const std::pair<U, std::pair<U, U> >& lhs, const std::pair<U, std::pair<U, U> >& rhs) {
		if (isCompoundMeasAll(_m->at(lhs.first).measType) && notCompoundMeasAll(_m->at(rhs.first).measType))
		{
			double lhsValue = 0.0;
			U increment(0);
			UINT32 vector_count(_m->at(lhs.first).vectorCount1), covariance_count;

			// Get the largest LHS value from the cluster
			switch (_m->at(lhs.first).measType)
			{
			case 'G':
			case 'X':
			case 'Y':
				for (UINT32 g(0); g < vector_count; ++g)
				{
					covariance_count = _m->at(lhs.first + increment).vectorCount2;
					if (fabs(_m->at(lhs.first + increment).term1) > lhsValue)				// X
						lhsValue = fabs(_m->at(lhs.first + increment).term1);
					if (fabs(_m->at(lhs.first + increment + 1).term1) > lhsValue)			// Y
						lhsValue = fabs(_m->at(lhs.first + increment + 1).term1);
					if (fabs(_m->at(lhs.first + increment + 2).term1) > lhsValue)			// Z
						lhsValue = fabs(_m->at(lhs.first + increment + 2).term1);
					increment += 3;							// move to covariances
					increment += (covariance_count * 3);	// skip over covariances
				}
				break;
			case 'D':
				//TRACE("%.9f\n", radians_to_degrees_(_m->at(lhs.first).term1));
				lhsValue = fabs(_m->at(lhs.first).term1);
				for (UINT32 d(1); d < vector_count; ++d)
				{
					//TRACE("%.9f\n", _m->at(lhs.first+d).term1);
					if (fabs(_m->at(lhs.first + d).term1) > lhsValue)
						lhsValue = fabs(_m->at(lhs.first + d).term1);
				}
				break;
			}
			//TRACE("LHS: %.2f; RHS: %.2f\n", fabs(lhsValue), fabs(_m->at(rhs.first).term1));
			return fabs(lhsValue) > fabs(_m->at(rhs.first).term1);
		}
		else if (notCompoundMeasAll(_m->at(lhs.first).measType) && isCompoundMeasAll(_m->at(rhs.first).measType))
		{
			double rhsValue = 0.0;
			U increment(0);
			UINT32 vector_count(_m->at(rhs.first).vectorCount1), covariance_count;

			// Get the largest RHS value from the cluster
			switch (_m->at(rhs.first).measType)
			{
			case 'G':
			case 'X':
			case 'Y':
				for (UINT32 g(0); g < vector_count; ++g)
				{
					covariance_count = _m->at(rhs.first + increment).vectorCount2;
					if (fabs(_m->at(rhs.first + increment).term1) > rhsValue)				// X
						rhsValue = fabs(_m->at(rhs.first + increment).term1);
					if (fabs(_m->at(rhs.first + increment + 1).term1) > rhsValue)			// Y
						rhsValue = fabs(_m->at(rhs.first + increment + 1).term1);
					if (fabs(_m->at(rhs.first + increment + 2).term1) > rhsValue)			// Z
						rhsValue = fabs(_m->at(rhs.first + increment + 2).term1);
					increment += 3;							// move to covariances
					increment += (covariance_count * 3);	// skip over covariances
				}
				break;
			case 'D':
				//TRACE("%.9f\n", radians_to_degrees_(_m->at(rhs.first).term1));
				rhsValue = fabs(_m->at(rhs.first).term1);
				for (UINT32 d(1); d < vector_count; ++d)
				{
					//TRACE("%.9f\n", _m->at(lhs.first+d).term1);
					if (fabs(_m->at(rhs.first + d).term1) > rhsValue)
						rhsValue = fabs(_m->at(rhs.first + d).term1);
				}
				break;
			}

			return fabs(_m->at(lhs.first).term1) > fabs(rhsValue);
		}
		else if (isCompoundMeasAll(_m->at(lhs.first).measType) && isCompoundMeasAll(_m->at(rhs.first).measType))
		{
			double lhsValue = 0.0;
			U increment(0);
			UINT32 vector_count(_m->at(lhs.first).vectorCount1), covariance_count;

			// Get the largest LHS value from the cluster
			switch (_m->at(lhs.first).measType)
			{
			case 'G':
			case 'X':
			case 'Y':
				for (UINT32 g(0); g < vector_count; ++g)
				{
					covariance_count = _m->at(lhs.first + increment).vectorCount2;
					if (fabs(_m->at(lhs.first + increment).term1) > lhsValue)		// X
						lhsValue = fabs(_m->at(lhs.first + increment).term1);
					if (fabs(_m->at(lhs.first + increment + 1).term1) > lhsValue)	// Y
						lhsValue = fabs(_m->at(lhs.first + increment + 1).term1);
					if (fabs(_m->at(lhs.first + increment + 2).term1) > lhsValue)	// Z
						lhsValue = fabs(_m->at(lhs.first + increment + 2).term1);
					increment += 3;							// move to covariances
					increment += (covariance_count * 3);	// skip over covariances
				}
				break;
			case 'D':
				//TRACE("%.9f\n", radians_to_degrees_(_m->at(lhs.first).term1));
				lhsValue = fabs(_m->at(lhs.first).term1);
				for (UINT32 d(1); d < vector_count; ++d)
				{
					//TRACE("%.9f\n", _m->at(lhs.first+d).term1);
					if (fabs(_m->at(lhs.first + d).term1) > lhsValue)
						lhsValue = fabs(_m->at(lhs.first + d).term1);
				}
				break;
			}

			double rhsValue = 0.0;
			increment = 0;
			vector_count = _m->at(rhs.first).vectorCount1;

			// Get the largest RHS value from the cluster
			switch (_m->at(rhs.first).measType)
			{
			case 'G':
			case 'X':
			case 'Y':
				for (UINT32 g(0); g < vector_count; ++g)
				{
					covariance_count = _m->at(rhs.first + increment).vectorCount2;
					if (fabs(_m->at(rhs.first + increment).term1) > rhsValue)				// X
						rhsValue = fabs(_m->at(rhs.first + increment).term1);
					if (fabs(_m->at(rhs.first + increment + 1).term1) > rhsValue)			// Y
						rhsValue = fabs(_m->at(rhs.first + increment + 1).term1);
					if (fabs(_m->at(rhs.first + increment + 2).term1) > rhsValue)			// Z
						rhsValue = fabs(_m->at(rhs.first + increment + 2).term1);
					increment += 3;							// move to covariances
					increment += (covariance_count * 3);	// skip over covariances
				}
				break;
			case 'D':
				//TRACE("%.9f\n", radians_to_degrees_(_m->at(rhs.first).term1));
				rhsValue = fabs(_m->at(rhs.first).term1);
				for (UINT32 d(1); d < vector_count; ++d)
				{
					//TRACE("%.9f\n", _m->at(lhs.first+d).term1);
					if (fabs(_m->at(rhs.first + d).term1) > rhsValue)
						rhsValue = fabs(_m->at(rhs.first + d).term1);
				}
			}

			return fabs(lhsValue) > fabs(rhsValue);
		}
		else
			return fabs(_m->at(lhs.first).term1) > fabs(_m->at(rhs.first).term1);
	}
private:
	std::vector<M>*	_m;
};


// M = measurement_t, U = UINT32
template<typename M, typename U>
class CompareMeasResidual_PairFirst
{
public:
	CompareMeasResidual_PairFirst(std::vector<M>* m)
		: _m(m) {}
	bool operator()(const std::pair<U, std::pair<U, U> >& lhs, const std::pair<U, std::pair<U, U> >& rhs) {
		if (isCompoundMeasAll(_m->at(lhs.first).measType) && notCompoundMeasAll(_m->at(rhs.first).measType))
		{
			double lhsValue = 0.0;
			U increment(0);
			UINT32 vector_count(_m->at(lhs.first).vectorCount1), covariance_count;

			// Get the largest LHS value from the cluster
			switch (_m->at(lhs.first).measType)
			{
			case 'G':
			case 'X':
			case 'Y':
				for (UINT32 g(0); g < vector_count; ++g)
				{
					covariance_count = _m->at(lhs.first + increment).vectorCount2;
					if (fabs(_m->at(lhs.first + increment).measCorr) > lhsValue)				// X
						lhsValue = fabs(_m->at(lhs.first + increment).measCorr);
					if (fabs(_m->at(lhs.first + increment + 1).measCorr) > lhsValue)			// Y
						lhsValue = fabs(_m->at(lhs.first + increment + 1).measCorr);
					if (fabs(_m->at(lhs.first + increment + 2).measCorr) > lhsValue)			// Z
						lhsValue = fabs(_m->at(lhs.first + increment + 2).measCorr);
					increment += 3;							// move to covariances
					increment += (covariance_count * 3);	// skip over covariances
				}
				break;
			case 'D':
				//TRACE("%.9f\n", radians_to_degrees_(_m->at(lhs.first).term1));
				lhsValue = fabs(_m->at(lhs.first).measCorr);
				for (UINT32 d(1); d < vector_count; ++d)
				{
					//TRACE("%.9f\n", _m->at(lhs.first+d).term1);
					if (fabs(_m->at(lhs.first + d).measCorr) > lhsValue)
						lhsValue = fabs(_m->at(lhs.first + d).measCorr);
				}
				break;
			}
			//TRACE("LHS: %.2f; RHS: %.2f\n", fabs(lhsValue), fabs(_m->at(rhs.first).measCorr));
			return fabs(lhsValue) > fabs(_m->at(rhs.first).measCorr);
		}
		else if (notCompoundMeasAll(_m->at(lhs.first).measType) && isCompoundMeasAll(_m->at(rhs.first).measType))
		{
			double rhsValue = 0.0;
			U increment(0);
			UINT32 vector_count(_m->at(rhs.first).vectorCount1), covariance_count;

			// Get the largest RHS value from the cluster
			switch (_m->at(rhs.first).measType)
			{
			case 'G':
			case 'X':
			case 'Y':
				for (UINT32 g(0); g < vector_count; ++g)
				{
					covariance_count = _m->at(rhs.first + increment).vectorCount2;
					if (fabs(_m->at(rhs.first + increment).measCorr) > rhsValue)				// X
						rhsValue = fabs(_m->at(rhs.first + increment).measCorr);
					if (fabs(_m->at(rhs.first + increment + 1).measCorr) > rhsValue)			// Y
						rhsValue = fabs(_m->at(rhs.first + increment + 1).measCorr);
					if (fabs(_m->at(rhs.first + increment + 2).measCorr) > rhsValue)			// Z
						rhsValue = fabs(_m->at(rhs.first + increment + 2).measCorr);
					increment += 3;							// move to covariances
					increment += (covariance_count * 3);	// skip over covariances
				}
				break;
			case 'D':
				//TRACE("%.9f\n", radians_to_degrees_(_m->at(rhs.first).term1));
				rhsValue = fabs(_m->at(rhs.first).measCorr);
				for (UINT32 d(1); d < vector_count; ++d)
				{
					//TRACE("%.9f\n", _m->at(lhs.first+d).term1);
					if (fabs(_m->at(rhs.first + d).measCorr) > rhsValue)
						rhsValue = fabs(_m->at(rhs.first + d).measCorr);
				}
				break;
			}

			return fabs(_m->at(lhs.first).measCorr) > fabs(rhsValue);
		}
		else if (isCompoundMeasAll(_m->at(lhs.first).measType) && isCompoundMeasAll(_m->at(rhs.first).measType))
		{
			double lhsValue = 0.0;
			U increment(0);
			UINT32 vector_count(_m->at(lhs.first).vectorCount1), covariance_count;

			// Get the largest LHS value from the cluster
			switch (_m->at(lhs.first).measType)
			{
			case 'G':
			case 'X':
			case 'Y':
				for (UINT32 g(0); g < vector_count; ++g)
				{
					covariance_count = _m->at(lhs.first + increment).vectorCount2;
					if (fabs(_m->at(lhs.first + increment).measCorr) > lhsValue)		// X
						lhsValue = fabs(_m->at(lhs.first + increment).measCorr);
					if (fabs(_m->at(lhs.first + increment + 1).measCorr) > lhsValue)	// Y
						lhsValue = fabs(_m->at(lhs.first + increment + 1).measCorr);
					if (fabs(_m->at(lhs.first + increment + 2).measCorr) > lhsValue)	// Z
						lhsValue = fabs(_m->at(lhs.first + increment + 2).measCorr);
					increment += 3;							// move to covariances
					increment += (covariance_count * 3);	// skip over covariances
				}
				break;
			case 'D':
				//TRACE("%.9f\n", radians_to_degrees_(_m->at(lhs.first).term1));
				lhsValue = fabs(_m->at(lhs.first).measCorr);
				for (UINT32 d(1); d < vector_count; ++d)
				{
					//TRACE("%.9f\n", _m->at(lhs.first+d).term1);
					if (fabs(_m->at(lhs.first + d).measCorr) > lhsValue)
						lhsValue = fabs(_m->at(lhs.first + d).measCorr);
				}
				break;
			}

			double rhsValue = 0.0;
			increment = 0;
			vector_count = _m->at(rhs.first).vectorCount1;

			// Get the largest RHS value from the cluster
			switch (_m->at(rhs.first).measType)
			{
			case 'G':
			case 'X':
			case 'Y':
				for (UINT32 g(0); g < vector_count; ++g)
				{
					covariance_count = _m->at(rhs.first + increment).vectorCount2;
					if (fabs(_m->at(rhs.first + increment).measCorr) > rhsValue)				// X
						rhsValue = fabs(_m->at(rhs.first + increment).measCorr);
					if (fabs(_m->at(rhs.first + increment + 1).measCorr) > rhsValue)			// Y
						rhsValue = fabs(_m->at(rhs.first + increment + 1).measCorr);
					if (fabs(_m->at(rhs.first + increment + 2).measCorr) > rhsValue)			// Z
						rhsValue = fabs(_m->at(rhs.first + increment + 2).measCorr);
					increment += 3;							// move to covariances
					increment += (covariance_count * 3);	// skip over covariances
				}
				break;
			case 'D':
				//TRACE("%.9f\n", radians_to_degrees_(_m->at(rhs.first).term1));
				rhsValue = fabs(_m->at(rhs.first).measCorr);
				for (UINT32 d(1); d < vector_count; ++d)
				{
					//TRACE("%.9f\n", _m->at(lhs.first+d).term1);
					if (fabs(_m->at(rhs.first + d).measCorr) > rhsValue)
						rhsValue = fabs(_m->at(rhs.first + d).measCorr);
				}
			}

			return fabs(lhsValue) > fabs(rhsValue);
		}
		else
			return fabs(_m->at(lhs.first).measCorr) > fabs(_m->at(rhs.first).measCorr);
	}
private:
	std::vector<M>* _m;
};


// M = measurement_t, U = UINT32
template<typename M, typename U>
class CompareMeasAdjSD_PairFirst
{
public:
	CompareMeasAdjSD_PairFirst(std::vector<M>* m)
		: _m(m) {}
	bool operator()(const std::pair<U, std::pair<U, U> >& lhs, const std::pair<U, std::pair<U, U> >& rhs) {
		if (isCompoundMeasAll(_m->at(lhs.first).measType) && notCompoundMeasAll(_m->at(rhs.first).measType))
		{
			double lhsValue = 0.0;
			U increment(0);
			UINT32 vector_count(_m->at(lhs.first).vectorCount1), covariance_count;

			// Get the largest LHS value from the cluster
			switch (_m->at(lhs.first).measType)
			{
			case 'G':
			case 'X':
			case 'Y':
				for (UINT32 g(0); g < vector_count; ++g)
				{
					covariance_count = _m->at(lhs.first + increment).vectorCount2;
					if (fabs(_m->at(lhs.first + increment).measAdjPrec) > lhsValue)				// X
						lhsValue = fabs(_m->at(lhs.first + increment).measAdjPrec);
					if (fabs(_m->at(lhs.first + increment + 1).measAdjPrec) > lhsValue)			// Y
						lhsValue = fabs(_m->at(lhs.first + increment + 1).measAdjPrec);
					if (fabs(_m->at(lhs.first + increment + 2).measAdjPrec) > lhsValue)			// Z
						lhsValue = fabs(_m->at(lhs.first + increment + 2).measAdjPrec);
					increment += 3;							// move to covariances
					increment += (covariance_count * 3);	// skip over covariances
				}
				break;
			case 'D':
				//TRACE("%.9f\n", radians_to_degrees_(_m->at(lhs.first).term1));
				lhsValue = fabs(_m->at(lhs.first).measAdjPrec);
				for (UINT32 d(1); d < vector_count; ++d)
				{
					//TRACE("%.9f\n", _m->at(lhs.first+d).term1);
					if (fabs(_m->at(lhs.first + d).measAdjPrec) > lhsValue)
						lhsValue = fabs(_m->at(lhs.first + d).measAdjPrec);
				}
				break;
			}
			//TRACE("LHS: %.2f; RHS: %.2f\n", fabs(lhsValue), fabs(_m->at(rhs.first).measAdjPrec));
			return fabs(lhsValue) > fabs(_m->at(rhs.first).measAdjPrec);
		}
		else if (notCompoundMeasAll(_m->at(lhs.first).measType) && isCompoundMeasAll(_m->at(rhs.first).measType))
		{
			double rhsValue = 0.0;
			U increment(0);
			UINT32 vector_count(_m->at(rhs.first).vectorCount1), covariance_count;

			// Get the largest RHS value from the cluster
			switch (_m->at(rhs.first).measType)
			{
			case 'G':
			case 'X':
			case 'Y':
				for (UINT32 g(0); g < vector_count; ++g)
				{
					covariance_count = _m->at(rhs.first + increment).vectorCount2;
					if (fabs(_m->at(rhs.first + increment).measAdjPrec) > rhsValue)				// X
						rhsValue = fabs(_m->at(rhs.first + increment).measAdjPrec);
					if (fabs(_m->at(rhs.first + increment + 1).measAdjPrec) > rhsValue)			// Y
						rhsValue = fabs(_m->at(rhs.first + increment + 1).measAdjPrec);
					if (fabs(_m->at(rhs.first + increment + 2).measAdjPrec) > rhsValue)			// Z
						rhsValue = fabs(_m->at(rhs.first + increment + 2).measAdjPrec);
					increment += 3;							// move to covariances
					increment += (covariance_count * 3);	// skip over covariances
				}
				break;
			case 'D':
				//TRACE("%.9f\n", radians_to_degrees_(_m->at(rhs.first).term1));
				rhsValue = fabs(_m->at(rhs.first).measAdjPrec);
				for (UINT32 d(1); d < vector_count; ++d)
				{
					//TRACE("%.9f\n", _m->at(lhs.first+d).term1);
					if (fabs(_m->at(rhs.first + d).measAdjPrec) > rhsValue)
						rhsValue = fabs(_m->at(rhs.first + d).measAdjPrec);
				}
				break;
			}

			return fabs(_m->at(lhs.first).measAdjPrec) > fabs(rhsValue);
		}
		else if (isCompoundMeasAll(_m->at(lhs.first).measType) && isCompoundMeasAll(_m->at(rhs.first).measType))
		{
			double lhsValue = 0.0;
			U increment(0);
			UINT32 vector_count(_m->at(lhs.first).vectorCount1), covariance_count;

			// Get the largest LHS value from the cluster
			switch (_m->at(lhs.first).measType)
			{
			case 'G':
			case 'X':
			case 'Y':
				for (UINT32 g(0); g < vector_count; ++g)
				{
					covariance_count = _m->at(lhs.first + increment).vectorCount2;
					if (fabs(_m->at(lhs.first + increment).measAdjPrec) > lhsValue)		// X
						lhsValue = fabs(_m->at(lhs.first + increment).measAdjPrec);
					if (fabs(_m->at(lhs.first + increment + 1).measAdjPrec) > lhsValue)	// Y
						lhsValue = fabs(_m->at(lhs.first + increment + 1).measAdjPrec);
					if (fabs(_m->at(lhs.first + increment + 2).measAdjPrec) > lhsValue)	// Z
						lhsValue = fabs(_m->at(lhs.first + increment + 2).measAdjPrec);
					increment += 3;							// move to covariances
					increment += (covariance_count * 3);	// skip over covariances
				}
				break;
			case 'D':
				//TRACE("%.9f\n", radians_to_degrees_(_m->at(lhs.first).term1));
				lhsValue = fabs(_m->at(lhs.first).measAdjPrec);
				for (UINT32 d(1); d < vector_count; ++d)
				{
					//TRACE("%.9f\n", _m->at(lhs.first+d).term1);
					if (fabs(_m->at(lhs.first + d).measAdjPrec) > lhsValue)
						lhsValue = fabs(_m->at(lhs.first + d).measAdjPrec);
				}
				break;
			}

			double rhsValue = 0.0;
			increment = 0;
			vector_count = _m->at(rhs.first).vectorCount1;

			// Get the largest RHS value from the cluster
			switch (_m->at(rhs.first).measType)
			{
			case 'G':
			case 'X':
			case 'Y':
				for (UINT32 g(0); g < vector_count; ++g)
				{
					covariance_count = _m->at(rhs.first + increment).vectorCount2;
					if (fabs(_m->at(rhs.first + increment).measAdjPrec) > rhsValue)				// X
						rhsValue = fabs(_m->at(rhs.first + increment).measAdjPrec);
					if (fabs(_m->at(rhs.first + increment + 1).measAdjPrec) > rhsValue)			// Y
						rhsValue = fabs(_m->at(rhs.first + increment + 1).measAdjPrec);
					if (fabs(_m->at(rhs.first + increment + 2).measAdjPrec) > rhsValue)			// Z
						rhsValue = fabs(_m->at(rhs.first + increment + 2).measAdjPrec);
					increment += 3;							// move to covariances
					increment += (covariance_count * 3);	// skip over covariances
				}
				break;
			case 'D':
				//TRACE("%.9f\n", radians_to_degrees_(_m->at(rhs.first).term1));
				rhsValue = fabs(_m->at(rhs.first).measAdjPrec);
				for (UINT32 d(1); d < vector_count; ++d)
				{
					//TRACE("%.9f\n", _m->at(lhs.first+d).term1);
					if (fabs(_m->at(rhs.first + d).measAdjPrec) > rhsValue)
						rhsValue = fabs(_m->at(rhs.first + d).measAdjPrec);
				}
			}

			return fabs(lhsValue) > fabs(rhsValue);
		}
		else
			return fabs(_m->at(lhs.first).measAdjPrec) > fabs(_m->at(rhs.first).measAdjPrec);
	}
private:
	std::vector<M>* _m;
};


// M = measurement_t, U = UINT32
template<typename M, typename U>
class CompareMeasNstat_PairFirst
{
public:
	CompareMeasNstat_PairFirst(std::vector<M>* m)
		:  _m(m) {}
	bool operator()(const std::pair<U, std::pair<U, U> >& lhs, const std::pair<U, std::pair<U, U> >& rhs) {
		if (isCompoundMeasAll(_m->at(lhs.first).measType) && notCompoundMeasAll(_m->at(rhs.first).measType))
		{
			double lhsValue = 0.0;
			U increment(0);			
			UINT32 vector_count(_m->at(lhs.first).vectorCount1), covariance_count;

			// Get the largest LHS value from the cluster
			switch (_m->at(lhs.first).measType)
			{
			case 'G':
			case 'X':
			case 'Y':				
				for (UINT32 g(0); g < vector_count; ++g)
				{					
					covariance_count = _m->at(lhs.first + increment).vectorCount2;
					if (fabs(_m->at(lhs.first + increment).NStat) > lhsValue)				// X
						lhsValue = fabs(_m->at(lhs.first + increment).NStat);				
					if (fabs(_m->at(lhs.first + increment + 1).NStat) > lhsValue)			// Y
						lhsValue = fabs(_m->at(lhs.first + increment + 1).NStat);
					if (fabs(_m->at(lhs.first + increment + 2).NStat) > lhsValue)			// Z
						lhsValue = fabs(_m->at(lhs.first + increment + 2).NStat);
					increment += 3;							// move to covariances
					increment += (covariance_count * 3);	// skip over covariances
				}
				break;
			case 'D':
				//TRACE("%.9f\n", radians_to_degrees_(_m->at(lhs.first).term1));
				lhsValue = fabs(_m->at(lhs.first).NStat);
				for (UINT32 d(1); d < vector_count; ++d)
				{
					//TRACE("%.9f\n", _m->at(lhs.first+d).term1);
					if (fabs(_m->at(lhs.first + d).NStat) > lhsValue)
						lhsValue = fabs( _m->at(lhs.first + d).NStat);
				}
				break;
			}
			//TRACE("LHS: %.2f; RHS: %.2f\n", fabs(lhsValue), fabs(_m->at(rhs.first).NStat));
			return fabs(lhsValue) > fabs(_m->at(rhs.first).NStat);
		}
		else if (notCompoundMeasAll(_m->at(lhs.first).measType) && isCompoundMeasAll(_m->at(rhs.first).measType))
		{			
			double rhsValue = 0.0;			
			U increment(0);			
			UINT32 vector_count(_m->at(rhs.first).vectorCount1), covariance_count;

			// Get the largest RHS value from the cluster
			switch (_m->at(rhs.first).measType)
			{
			case 'G':
			case 'X':
			case 'Y':				
				for (UINT32 g(0); g < vector_count; ++g)
				{
					covariance_count = _m->at(rhs.first + increment).vectorCount2;
					if (fabs(_m->at(rhs.first + increment).NStat) > rhsValue)				// X
						rhsValue = fabs(_m->at(rhs.first + increment).NStat);				
					if (fabs(_m->at(rhs.first + increment + 1).NStat) > rhsValue)			// Y
						rhsValue = fabs(_m->at(rhs.first + increment + 1).NStat);
					if (fabs(_m->at(rhs.first + increment + 2).NStat) > rhsValue)			// Z
						rhsValue = fabs(_m->at(rhs.first + increment + 2).NStat);
					increment += 3;							// move to covariances
					increment += (covariance_count * 3);	// skip over covariances
				}
				break;
			case 'D':
				//TRACE("%.9f\n", radians_to_degrees_(_m->at(rhs.first).term1));
				rhsValue = fabs(_m->at(rhs.first).NStat);
				for (UINT32 d(1); d < vector_count; ++d)
				{
					//TRACE("%.9f\n", _m->at(lhs.first+d).term1);
					if (fabs(_m->at(rhs.first + d).NStat) > rhsValue)
						rhsValue =  fabs(_m->at(rhs.first + d).NStat);
				}
				break;
			}

			return fabs(_m->at(lhs.first).NStat) > fabs(rhsValue);
		}
		else if (isCompoundMeasAll(_m->at(lhs.first).measType) && isCompoundMeasAll(_m->at(rhs.first).measType))
		{
			double lhsValue = 0.0;
			U increment(0);
			UINT32 vector_count(_m->at(lhs.first).vectorCount1), covariance_count;

			// Get the largest LHS value from the cluster
			switch (_m->at(lhs.first).measType)
			{
			case 'G':
			case 'X':
			case 'Y':
				for (UINT32 g(0); g < vector_count; ++g)
				{
					covariance_count = _m->at(lhs.first + increment).vectorCount2;
					if (fabs(_m->at(lhs.first + increment).NStat) > lhsValue)		// X
						lhsValue = fabs(_m->at(lhs.first + increment).NStat);			
					if (fabs(_m->at(lhs.first + increment + 1).NStat) > lhsValue)	// Y
						lhsValue = fabs(_m->at(lhs.first + increment + 1).NStat);
					if (fabs(_m->at(lhs.first + increment + 2).NStat) > lhsValue)	// Z
						lhsValue = fabs(_m->at(lhs.first + increment + 2).NStat);
					increment += 3;							// move to covariances
					increment += (covariance_count * 3);	// skip over covariances
				}
				break;
			case 'D':
				//TRACE("%.9f\n", radians_to_degrees_(_m->at(lhs.first).term1));
				lhsValue = fabs(_m->at(lhs.first).NStat);
				for (UINT32 d(1); d < vector_count; ++d)
				{
					//TRACE("%.9f\n", _m->at(lhs.first+d).term1);
					if (fabs(_m->at(lhs.first+d).NStat) > lhsValue)
						lhsValue =  fabs(_m->at(lhs.first+d).NStat);
				}
				break;
			}

			double rhsValue = 0.0;
			increment = 0;
			vector_count = _m->at(rhs.first).vectorCount1;

			// Get the largest RHS value from the cluster
			switch (_m->at(rhs.first).measType)
			{
			case 'G':
			case 'X':
			case 'Y':				
				for (UINT32 g(0); g < vector_count; ++g)
				{
					covariance_count = _m->at(rhs.first + increment).vectorCount2;
					if (fabs(_m->at(rhs.first + increment).NStat) > rhsValue)				// X
						rhsValue = fabs(_m->at(rhs.first + increment).NStat);
					if (fabs(_m->at(rhs.first + increment + 1).NStat) > rhsValue)			// Y
						rhsValue = fabs(_m->at(rhs.first + increment + 1).NStat);
					if (fabs(_m->at(rhs.first + increment + 2).NStat) > rhsValue)			// Z
						rhsValue = fabs(_m->at(rhs.first + increment + 2).NStat);
					increment += 3;							// move to covariances
					increment += (covariance_count * 3);	// skip over covariances
				}
				break;
			case 'D':
				//TRACE("%.9f\n", radians_to_degrees_(_m->at(rhs.first).term1));
				rhsValue = fabs(_m->at(rhs.first).NStat);
				for (UINT32 d(1); d < vector_count; ++d)
				{
					//TRACE("%.9f\n", _m->at(lhs.first+d).term1);
					if (fabs(_m->at(rhs.first+d).NStat) > rhsValue)
						rhsValue =  fabs(_m->at(rhs.first+d).NStat);
				}
			}

			return fabs(lhsValue) > fabs(rhsValue);
		}
		else
			return fabs(_m->at(lhs.first).NStat) > fabs(_m->at(rhs.first).NStat);
	}
private:
	std::vector<M>*	_m;
};


// M = measurement_t
template<typename M, typename C = char>
class CompareMeasTypeT
{
public:
	CompareMeasTypeT(const C& t) :  _type(t) {}
	inline void SetComparand(const char& t) { _type = t; }
	bool operator()(const M& msr) {
		return msr.measType == _type;
	}

private:
	C _type;
};
	

// M = measurement_t
template<typename M>
class CompareValidMeasTypeT
{
public:
	CompareValidMeasTypeT(const char& t) :  _type(t) {}
	inline void SetComparand(const char& t) { _type = t; }
	bool operator()(const M& msr) {
		return msr.ignore == false &&
			msr.measType == _type;
	}

private:
	char _type;
};


// M = measurement_t, U = UINT32

// M = CDnaMeasurement, U = UINT32
// M = CDnaGpsPoint or CDnaGpsBaseline (derived from CDnaMeasurement), U = UINT32
// M = CDnaMeasurement, U = UINT32

// M = CDnaGpsPoint or CDnaGpsBaseline (derived from CDnaMeasurement), U = UINT32
// M = CDnaGpsPoint or CDnaGpsBaseline (derived from CDnaMeasurement), U = UINT32
// M = CDnaGpsPoint or CDnaGpsBaseline (derived from CDnaMeasurement), U = UINT32

// M = CDnaMeasurement
// M = CDnaMeasurement

// M = measurement_t, U = UINT32
// Compare functor - searches for all stations that appear in the list
template<typename U, typename M, typename C>
class CompareFreeClusterAllStns
{
public:
	CompareFreeClusterAllStns(std::vector<U>* u, std::vector<M>* m, const C& c) :  _u(u), _m(m), _c(c) {}
	bool operator()(const U& amlindex) {
		// one-station measurement types
		// is this station on the list?
		// Is station 1 on inner or junction lists?
		switch (_c)
		{
		case 'H':	// Orthometric height
		case 'R':	// Ellipsoidal height
		case 'I':	// Astronomic latitude
		case 'J':	// Astronomic longitude
		case 'P':	// Geodetic latitude
		case 'Q':	// Geodetic longitude
		case 'Y':	// GPS point cluster
			if (binary_search(_u->begin(), _u->end(), _m->at(amlindex).station1))
				return true;
			return false;
		}
		// two-station measurement types
		// is this station on the list?
		switch (_c)
		{
		case 'D':	// Direction set
		case 'B':	// Geodetic azimuth
		case 'K':	// Astronomic azimuth
		case 'C':	// Chord dist
		case 'E':	// Ellipsoid arc
		case 'M':	// MSL arc
		case 'S':	// Slope distance
		case 'L':	// Level difference
		case 'V':	// Zenith distance
		case 'Z':	// Vertical angle
		case 'G':	// GPS Baseline (treat as single-baseline cluster)
		case 'X':	// GPS Baseline cluster
			if (binary_search(_u->begin(), _u->end(), _m->at(amlindex).station1) &&
				binary_search(_u->begin(), _u->end(), _m->at(amlindex).station2))
				return true;
			return false;
		}
		// three-station measurement types
		// is this station on the list?
		if (binary_search(_u->begin(), _u->end(), _m->at(amlindex).station1) &&
			binary_search(_u->begin(), _u->end(), _m->at(amlindex).station2) &&
			binary_search(_u->begin(), _u->end(), _m->at(amlindex).station3))
			return true;
		return false;
	}
private:
	std::vector<U>*	_u;
	std::vector<M>*	_m;
	char		_c;
};


#endif /* DNATEMPLATESTNMSRFUNCS_H_ */
