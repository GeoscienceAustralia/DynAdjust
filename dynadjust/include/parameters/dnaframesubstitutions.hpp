//============================================================================
// Name         : dnaframesubstitutions.hpp
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
// Description  : Reference frame substitutions file
//============================================================================

#ifndef DNAIOFRX_H_
#define DNAIOFRX_H_

#if defined(_MSC_VER)
	#if defined(LIST_INCLUDES_ON_BUILD) 
		#pragma message("  " __FILE__) 
	#endif
#endif

#include <include/config/dnatypes-fwd.hpp>
#include <include/functions/dnatemplatedatetimefuncs.hpp>
#include <include/parameters/dnaconsts-datums.hpp>
#include <include/parameters/dnaepsg.hpp>

namespace dynadjust {
namespace frame_substitutions {

/////////////////////////////////////////////////////////////
// Custom type to manage frame substitutions
template <class T1 = std::string, class T2 = UINT32, class T3 = double>
class frame_substitutions_t
{
public:
	frame_substitutions_t() {
		frame_name = WGS84_s;
		frame_epsg = WGS84_i_xyz;
		frame_desc = "";
		substitute_name = ITRF2014_s;
		substitute_epsg = ITRF2014_i_xyz;

		parameters_[0] = 0.0;
		parameters_[1] = 0.0;
		parameters_[2] = 0.0;
		parameters_[3] = 0.0;
		parameters_[4] = 0.0;
		parameters_[5] = 0.0;
		parameters_[6] = 0.0;
	}

	const boost::gregorian::date getFromEpoch() const { return from_epoch; }
	const boost::gregorian::date getToEpoch() const { return to_epoch; }
	const T1 getFrameName() const { return frame_name; }
	const T1 getSubstituteName() const { return substitute_name; }

protected:
	T1 frame_name;              // frame name
    T2 frame_epsg;              // frame epsg code
    T1 frame_desc;              // frame description
    T1 substitute_name;         // substitution frame name
    T2 substitute_epsg;         // substitution frame epsg code
    boost::gregorian::date from_epoch;            // start date of substitution
    boost::gregorian::date to_epoch;              // end date of substitution
    T3 parameters_[7];		    // transformation parameters (if any)
};

typedef frame_substitutions_t<std::string, UINT32, double> frame_substitutions;

typedef std::shared_ptr< frame_substitutions_t<std::string, UINT32, double> > frameSubsPtr;
typedef std::vector<frameSubsPtr> vframeSubsPtr, *pvframeSubsPtr;
typedef vframeSubsPtr::iterator _it_vframesubptr;


// WGS84 (Transit) to ITRF90
template <class T1 = std::string, class T2 = UINT32, class T3 = double>
class WGS84_TRANSIT_ITRF90 : public frame_substitutions_t<T1, T2, T3>
{
public:
	WGS84_TRANSIT_ITRF90() 
	{
		frame_substitutions::frame_name = WGS84_transit_s;
		frame_substitutions::frame_epsg = WGS84_transit_i_xyz;
		frame_substitutions::frame_desc = "";
		frame_substitutions::substitute_name = ITRF1990_s_brief;
		frame_substitutions::substitute_epsg = ITRF1990_i_xyz;

		frame_substitutions::parameters_[0] = 0.060;
		frame_substitutions::parameters_[1] = -0.517;
		frame_substitutions::parameters_[2] = -0.223;
		frame_substitutions::parameters_[3] = -0.011;
		frame_substitutions::parameters_[4] = 0.0183;
		frame_substitutions::parameters_[5] = -0.0003;
		frame_substitutions::parameters_[6] = 0.0070;
		
		frame_substitutions::from_epoch = dateFromString<boost::gregorian::date>("01.01.1987");
		frame_substitutions::to_epoch = dateFromString<boost::gregorian::date>("01.01.1994");

	};
	~WGS84_TRANSIT_ITRF90() {}
};


// WGS84 to ITRF90
template <class T1 = std::string, class T2 = UINT32, class T3 = double>
class WGS84_ITRF90 : public frame_substitutions_t<T1, T2, T3>
{
public:
	WGS84_ITRF90()
	{
		frame_substitutions::frame_name = WGS84_s;
		frame_substitutions::frame_epsg = WGS84_i_xyz;
		frame_substitutions::frame_desc = "";
		frame_substitutions::substitute_name = ITRF1990_s_brief;
		frame_substitutions::substitute_epsg = ITRF1990_i_xyz;

		frame_substitutions::parameters_[0] = 0.060;
		frame_substitutions::parameters_[1] = -0.517;
		frame_substitutions::parameters_[2] = -0.223;
		frame_substitutions::parameters_[3] = -0.011;
		frame_substitutions::parameters_[4] = 0.0183;
		frame_substitutions::parameters_[5] = -0.0003;
		frame_substitutions::parameters_[6] = 0.0070;

		frame_substitutions::from_epoch = dateFromString<boost::gregorian::date>("01.01.1987");
		frame_substitutions::to_epoch = dateFromString<boost::gregorian::date>("01.01.1994");

	};
	~WGS84_ITRF90() {}
};


// WGS84 (G730) to ITRF91
template <class T1 = std::string, class T2 = UINT32, class T3 = double>
class WGS84_G730_ITRF91 : public frame_substitutions_t<T1, T2, T3>
{
public:
	WGS84_G730_ITRF91() 
	{
		frame_substitutions::frame_name = WGS84_G730_s;
		frame_substitutions::frame_epsg = WGS84_G730_i_xyz;
		frame_substitutions::frame_desc = "";
		frame_substitutions::substitute_name = ITRF1991_s_brief;
		frame_substitutions::substitute_epsg = ITRF1991_i_xyz;

		frame_substitutions::from_epoch = dateFromString<boost::gregorian::date>("02.01.1994");
		frame_substitutions::to_epoch = dateFromString<boost::gregorian::date>("28.09.1996");

	};
	virtual ~WGS84_G730_ITRF91() {}
};


// WGS84 to ITRF91
template <class T1 = std::string, class T2 = UINT32, class T3 = double>
class WGS84_ITRF91 : public frame_substitutions_t<T1, T2, T3>
{
public:
	WGS84_ITRF91()
	{
		frame_substitutions::frame_name = WGS84_s;
		frame_substitutions::frame_epsg = WGS84_i_xyz;
		frame_substitutions::frame_desc = "";
		frame_substitutions::substitute_name = ITRF1991_s_brief;
		frame_substitutions::substitute_epsg = ITRF1991_i_xyz;

		frame_substitutions::from_epoch = dateFromString<boost::gregorian::date>("02.01.1994");
		frame_substitutions::to_epoch = dateFromString<boost::gregorian::date>("28.09.1996");

	};
	virtual ~WGS84_ITRF91() {}
};


// WGS84 (G873) to ITRF94
template <class T1 = std::string, class T2 = UINT32, class T3 = double>
class WGS84_G873_ITRF94 : public frame_substitutions_t<T1, T2, T3>
{
public:
	WGS84_G873_ITRF94() 
	{
		frame_substitutions::frame_name = WGS84_G873_s;
		frame_substitutions::frame_epsg = WGS84_G873_i_xyz;
		frame_substitutions::frame_desc = "";
		frame_substitutions::substitute_name = ITRF1994_s_brief;
		frame_substitutions::substitute_epsg = ITRF1994_i_xyz;

		frame_substitutions::from_epoch = dateFromString<boost::gregorian::date>("29.09.1996");
		frame_substitutions::to_epoch = dateFromString<boost::gregorian::date>("19.01.2002");

	};
	virtual ~WGS84_G873_ITRF94() {}
};

// WGS84 to ITRF94
template <class T1 = std::string, class T2 = UINT32, class T3 = double>
class WGS84_ITRF94 : public frame_substitutions_t<T1, T2, T3>
{
public:
	WGS84_ITRF94()
	{
		frame_substitutions::frame_name = WGS84_s;
		frame_substitutions::frame_epsg = WGS84_i_xyz;
		frame_substitutions::frame_desc = "";
		frame_substitutions::substitute_name = ITRF1994_s_brief;
		frame_substitutions::substitute_epsg = ITRF1994_i_xyz;

		frame_substitutions::from_epoch = dateFromString<boost::gregorian::date>("29.09.1996");
		frame_substitutions::to_epoch = dateFromString<boost::gregorian::date>("19.01.2002");

	};
	virtual ~WGS84_ITRF94() {}
};

// WGS84 (G1150) to ITRF2000
template <class T1 = std::string, class T2 = UINT32, class T3 = double>
class WGS84_G1150_ITRF2000 : public frame_substitutions_t<T1, T2, T3>
{
public:
	WGS84_G1150_ITRF2000()
	{
		frame_substitutions::frame_name = WGS84_G1150_s;
		frame_substitutions::frame_epsg = WGS84_G1150_i_xyz;
		frame_substitutions::frame_desc = "";
		frame_substitutions::substitute_name = ITRF2000_s;
		frame_substitutions::substitute_epsg = ITRF2000_i_xyz;

		frame_substitutions::from_epoch = dateFromString<boost::gregorian::date>("20.01.2002");
		frame_substitutions::to_epoch = dateFromString<boost::gregorian::date>("06.05.2012");

	};
	virtual ~WGS84_G1150_ITRF2000() {}
};

// WGS84 to ITRF2000
template <class T1 = std::string, class T2 = UINT32, class T3 = double>
class WGS84_ITRF2000 : public frame_substitutions_t<T1, T2, T3>
{
public:
	WGS84_ITRF2000()
	{
		frame_substitutions::frame_name = WGS84_s;
		frame_substitutions::frame_epsg = WGS84_i_xyz;
		frame_substitutions::frame_desc = "";
		frame_substitutions::substitute_name = ITRF2000_s;
		frame_substitutions::substitute_epsg = ITRF2000_i_xyz;

		frame_substitutions::from_epoch = dateFromString<boost::gregorian::date>("20.01.2002");
		frame_substitutions::to_epoch = dateFromString<boost::gregorian::date>("06.05.2012");

	};
	virtual ~WGS84_ITRF2000() {}
};

// WGS84 (G1674) to ITRF2008
template <class T1 = std::string, class T2 = UINT32, class T3 = double>
class WGS84_G1674_ITRF2008 : public frame_substitutions_t<T1, T2, T3>
{
public:
	WGS84_G1674_ITRF2008() 
	{
		frame_substitutions::frame_name = WGS84_G1674_s;
		frame_substitutions::frame_epsg = WGS84_G1674_i_xyz;
		frame_substitutions::frame_desc = "";
		frame_substitutions::substitute_name = ITRF2008_s;
		frame_substitutions::substitute_epsg = ITRF2008_i_xyz;

		frame_substitutions::from_epoch = dateFromString<boost::gregorian::date>("07.05.2012");
		frame_substitutions::to_epoch = dateFromString<boost::gregorian::date>("15.10.2013");

	};
	virtual ~WGS84_G1674_ITRF2008() {}
};

// WGS84 to ITRF2008
template <class T1 = std::string, class T2 = UINT32, class T3 = double>
class WGS84_ITRF2008_1 : public frame_substitutions_t<T1, T2, T3>
{
public:
	WGS84_ITRF2008_1()
	{
		frame_substitutions::frame_name = WGS84_s;
		frame_substitutions::frame_epsg = WGS84_i_xyz;
		frame_substitutions::frame_desc = "";
		frame_substitutions::substitute_name = ITRF2008_s;
		frame_substitutions::substitute_epsg = ITRF2008_i_xyz;

		frame_substitutions::from_epoch = dateFromString<boost::gregorian::date>("07.05.2012");
		frame_substitutions::to_epoch = dateFromString<boost::gregorian::date>("15.10.2013");

	};
	virtual ~WGS84_ITRF2008_1() {}
};

// WGS84 (G1762) to ITRF2008
template <class T1 = std::string, class T2 = UINT32, class T3 = double>
class WGS84_G1762_ITRF2008 : public frame_substitutions_t<T1, T2, T3>
{
public:
	WGS84_G1762_ITRF2008() 
	{
		frame_substitutions::frame_name = WGS84_G1762_s;
		frame_substitutions::frame_epsg = WGS84_G1762_i_xyz;
		frame_substitutions::frame_desc = "";
		frame_substitutions::substitute_name = ITRF2008_s;
		frame_substitutions::substitute_epsg = ITRF2008_i_xyz;

		frame_substitutions::from_epoch = dateFromString<boost::gregorian::date>("16.10.2013");
		frame_substitutions::to_epoch = dateFromString<boost::gregorian::date>("02.01.2021");

	};
	virtual ~WGS84_G1762_ITRF2008() {}
};

// WGS84 to ITRF2008
template <class T1 = std::string, class T2 = UINT32, class T3 = double>
class WGS84_ITRF2008_2 : public frame_substitutions_t<T1, T2, T3>
{
public:
	WGS84_ITRF2008_2()
	{
		frame_substitutions::frame_name = WGS84_s;
		frame_substitutions::frame_epsg = WGS84_i_xyz;
		frame_substitutions::frame_desc = "";
		frame_substitutions::substitute_name = ITRF2008_s;
		frame_substitutions::substitute_epsg = ITRF2008_i_xyz;

		frame_substitutions::from_epoch = dateFromString<boost::gregorian::date>("16.10.2013");
		frame_substitutions::to_epoch = dateFromString<boost::gregorian::date>("02.01.2021");

	};
	virtual ~WGS84_ITRF2008_2() {}
};

// WGS84 (G2139) to ITRF2014
template <class T1 = std::string, class T2 = UINT32, class T3 = double>
class WGS84_G2139_ITRF2014 : public frame_substitutions_t<T1, T2, T3>
{
public:
	WGS84_G2139_ITRF2014() 
	{
		frame_substitutions::frame_name = WGS84_G2139_s;
		frame_substitutions::frame_epsg = WGS84_G2139_i_xyz;
		frame_substitutions::frame_desc = "";
		frame_substitutions::substitute_name = ITRF2014_s;
		frame_substitutions::substitute_epsg = ITRF2014_i_xyz;

		frame_substitutions::from_epoch = dateFromString<boost::gregorian::date>("02.01.2021");
		frame_substitutions::to_epoch = boost::gregorian::day_clock::local_day() + boost::gregorian::years(100);

	};
	virtual ~WGS84_G2139_ITRF2014() {}
};

// WGS84 to ITRF2014
template <class T1 = std::string, class T2 = UINT32, class T3 = double>
class WGS84_ITRF2014 : public frame_substitutions_t<T1, T2, T3>
{
public:
	WGS84_ITRF2014()
	{
		frame_substitutions::frame_name = WGS84_s;
		frame_substitutions::frame_epsg = WGS84_i_xyz;
		frame_substitutions::frame_desc = "";
		frame_substitutions::substitute_name = ITRF2014_s;
		frame_substitutions::substitute_epsg = ITRF2014_i_xyz;

		frame_substitutions::from_epoch = dateFromString<boost::gregorian::date>("02.01.2021");
		frame_substitutions::to_epoch = boost::gregorian::day_clock::local_day() + boost::gregorian::years(100);

	};
	virtual ~WGS84_ITRF2014() {}
};

// SIRGAS95 to ITRF94
template <class T1 = std::string, class T2 = UINT32, class T3 = double>
class SIRGAS95_ITRF94 : public frame_substitutions_t<T1, T2, T3> {
   public:
    SIRGAS95_ITRF94() {
        frame_substitutions::frame_name = SIRGAS_1995_s;
        frame_substitutions::frame_epsg = SIRGAS_1995_i_xyz;
        frame_substitutions::frame_desc = "";
        frame_substitutions::substitute_name = ITRF1994_s;
        frame_substitutions::substitute_epsg = ITRF1994_i_xyz;

        frame_substitutions::from_epoch = dateFromString<boost::gregorian::date>("01.01.1900");
        frame_substitutions::to_epoch = boost::gregorian::day_clock::local_day() + boost::gregorian::years(100);
    };
    virtual ~SIRGAS95_ITRF94() {}
};

// SIRGAS2000 to ITRF2000
template <class T1 = std::string, class T2 = UINT32, class T3 = double>
class SIRGAS2000_ITRF2000 : public frame_substitutions_t<T1, T2, T3> {
   public:
    SIRGAS2000_ITRF2000() {
        frame_substitutions::frame_name = SIRGAS_2000_s;
        frame_substitutions::frame_epsg = SIRGAS_2000_i_xyz;
        frame_substitutions::frame_desc = "";
        frame_substitutions::substitute_name = ITRF2000_s;
        frame_substitutions::substitute_epsg = ITRF2000_i_xyz;

        frame_substitutions::from_epoch = dateFromString<boost::gregorian::date>("01.01.1900");
        frame_substitutions::to_epoch = boost::gregorian::day_clock::local_day() + boost::gregorian::years(100);
    };
    virtual ~SIRGAS2000_ITRF2000() {}
};

// FUNCTIONS

// Used to sort site tuples by WGS84 starting epoch
template <typename T>
class CompareSubstitutionByEpoch {
public:
	bool operator()(const T& left, const T& right) {
		// if the starting epochs are equal
		if (left.from_epoch == right.from_epoch)
			// sort on ending epochs
			return left.to_epoch < right.to_epoch;
		// sort on starting epoch
		return left.from_epoch < right.from_epoch;
	}
};

template<typename T, typename S>
struct CompareSubstituteOnFrameName
{
	//public:
	bool operator()(const std::shared_ptr<T>& lhs, const std::shared_ptr<T>& rhs) const {
		if (lhs->getFrameName() == rhs->getFrameName())
			return lhs->getFromEpoch() < rhs->getFromEpoch();
		return pair_firstless(lhs->getFrameName(), rhs->getFrameName());
	}
	bool operator()(const std::shared_ptr<T>& lhs, const S& rhs) const {
		return pair_firstless(lhs->getFrameName(), rhs);
	}
	bool operator()(const S& lhs, const std::shared_ptr<T>& rhs) const {
		return pair_firstless(lhs, rhs->getFrameName());
	}
	//private:
	bool pair_firstless(const S& s1, const S& s2) const {
		return s1 < s2;
	}
};


// tweak the binary search so it returns the iterator of the object found
template<typename T, typename S>
T binary_search_substitution(T begin, T end, S value)
{
	T i = lower_bound(begin, end, value, CompareSubstituteOnFrameName< frame_substitutions_t<std::string, UINT32, double> , std::string>());
	if (i != end && i->get()->getFrameName() == value)
		return i;
	else
		return end;
}
/////////////////////////////////////////////////////////////




}	// namespace frame_substitutions
}	// namespace dynadjust

#endif
