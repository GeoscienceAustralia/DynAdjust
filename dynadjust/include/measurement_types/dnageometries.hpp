//============================================================================
// Name         : dnageometries.hpp
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
// Description  : Basic geometry point 
//============================================================================

#ifndef DNAGEOMETRIES_H_
#define DNAGEOMETRIES_H_

#if defined(_MSC_VER)
	#if defined(LIST_INCLUDES_ON_BUILD) 
		#pragma message("  " __FILE__) 
	#endif
#endif

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/geometries/register/point.hpp>
#include <boost/geometry/geometries/multi_polygon.hpp>

using namespace std;
using namespace boost::geometry;

namespace dynadjust {
namespace geometries {

template <class T>
class dnaGeometryPoint
{
public:
	// dnaGeometryPoint()
	// 	: m_east_long(0.), m_north_lat(0.) 
	// {}

	// dnaGeometryPoint(const T& east_long, const T& north_lat)
	// 	: m_east_long(east_long), m_north_lat(north_lat)
	// {}

	T get_east_long() const {
		return m_east_long;
	}
	
	T get_north_lat() const {
		return m_north_lat;
	}
	
	void set_east_long(T east_long) {
		m_east_long = east_long;
	}
	
	void set_north_lat(T north_lat) {
		m_north_lat = north_lat;
	}
	
private:
	
	T m_east_long;
	T m_north_lat;
};

typedef boost::geometry::model::polygon<dnaGeometryPoint <double>> dnaGeometryPolygon;
typedef boost::geometry::strategy::centroid::bashein_detmer<dnaGeometryPoint<double>, dnaGeometryPoint<double>> dnaBasheinStrategy;
typedef boost::geometry::strategy::centroid::average<dnaGeometryPoint<double>, dnaGeometryPoint<double>> dnaAverageStrategy;

}	// namespace geometries
}	// namespace dynadjust


BOOST_GEOMETRY_REGISTER_POINT_2D_GET_SET(
	dynadjust::geometries::dnaGeometryPoint<double>,
	double,
	boost::geometry::cs::spherical_equatorial<boost::geometry::degree>,
	dynadjust::geometries::dnaGeometryPoint<double>::get_east_long,
	dynadjust::geometries::dnaGeometryPoint<double>::get_north_lat,
	dynadjust::geometries::dnaGeometryPoint<double>::set_east_long,
	dynadjust::geometries::dnaGeometryPoint<double>::set_north_lat)

#endif /* DNAGEOMETRIES_H_ */
