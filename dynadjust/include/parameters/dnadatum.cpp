//============================================================================
// Name         : dnadatum.hpp
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
// Description  : DynAdjust Datum Library
//============================================================================

#include <include/parameters/dnadatum.hpp>
#include <include/parameters/dnaepsg.hpp>

using namespace dynadjust::epsg;

namespace dynadjust {
namespace datum_parameters {

// default to GDA2020
CDnaDatum::CDnaDatum(void)
	// "GDA2020 (3D)" - geographic
	: epsgCode_(DEFAULT_EPSG_U)
	, epoch_(boost::gregorian::date(2020, 1, 1))
	, datumType_(STATIC_DATUM)
	, datumName_(DEFAULT_DATUM)
{
	ellipsoid_.SetEllipsoid(epsgCode_);
}

CDnaDatum::CDnaDatum(const UINT32& epsgCode)
{
	// Set datum with default epoch of today
	SetDatum(epsgCode);
}

CDnaDatum::CDnaDatum(const UINT32& epsgCode, const boost::gregorian::date& epoch)
{
	SetDatum(epsgCode, epoch);
}

CDnaDatum::CDnaDatum(const std::string& epsgCode, const std::string& epoch)
{
	SetDatumFromEpsg(epsgCode, epoch);
}
	
//// datumName corresponds to the name given in epsg
//// epoch is a space delimited string: YYY MM DD
//CDnaDatum::CDnaDatum(const std::string& datumName, const std::string& epoch)
//{
//	datumName_ = datumName;
//	SetDatumFromName(datumName, epoch);
//}

//CDnaDatum::CDnaDatum(const CDnaDatum& newDatum)
//{
//	epsgCode_ = newDatum.epsgCode_;
//	ellipsoid_ = newDatum.ellipsoid_;
//	epoch_ = newDatum.epoch_;
//	datumType_ = newDatum.datumType_;
//	datumName_ = newDatum.datumName_;
//}

CDnaDatum& CDnaDatum::operator=(const CDnaDatum& rhs)
{
	if (this == &rhs)	// check for assignment to self!
		return *this;

	epsgCode_ = rhs.epsgCode_;
	ellipsoid_ = rhs.ellipsoid_;
	epoch_ = rhs.epoch_;
	datumType_ = rhs.datumType_;
	datumName_ = rhs.datumName_;

	return *this;
}

bool CDnaDatum::operator==(const CDnaDatum& rhs) const
{
	return (
		epsgCode_ == rhs.epsgCode_ &&
		ellipsoid_ == rhs.ellipsoid_ &&
		epoch_ == rhs.epoch_ &&
		datumType_ == rhs.datumType_
		);
}

//bool CDnaDatum::isSameFrame(const CDnaDatum& rhs) const
//{
//	return (
//		epsgCode_ == rhs.epsgCode_ &&
//		ellipsoid_ == rhs.ellipsoid_ &&
//		datumType_ == rhs.datumType_
//		);
//}

void CDnaDatum::initialiseDatumFromEpsgCode()
{
	// Initialise ellipsoid using epsg code and set epoch if static
	ellipsoid_.SetEllipsoid(epsgCode_);
	
	// datum type?
	if (isEpsgDatumStatic(epsgCode_))
		datumType_ = STATIC_DATUM;
	else
		datumType_ = DYNAMIC_DATUM;
	
	// Set name
	datumName_ = datumFromEpsgCode<std::string, UINT32>(epsgCode_);
}
	
// returns:
//  - dd.mm.yyyy for datums that have a reference epoch
//  - "" for ensembles that do not have an epoch
std::string CDnaDatum::GetEpoch_s() const
{
	if (isEpsgWGS84Ensemble(epsgCode_) && isTimeImmemorial(epoch_))
		return "";
	else if (isTimeImmemorial(epoch_))
		return "";
	return stringFromDate<boost::gregorian::date>(epoch_);
}
	
std::string CDnaDatum::GetEpsgCode_s() const 
{
	char epsgCode[8];
	snprintf(epsgCode, sizeof(epsgCode), "%d", epsgCode_);
	return std::string(epsgCode); 
}

void CDnaDatum::SetDatum(const UINT32& epsgCode)
{
	epsgCode_ = epsgCode;
	epoch_ = dateFromString<boost::gregorian::date>(referenceepochFromEpsgCode(epsgCode_));

	// set datum parameters using epsg code
	initialiseDatumFromEpsgCode();
}

void CDnaDatum::SetDatum(const UINT32& epsgCode, const boost::gregorian::date& epoch)
{
	epsgCode_ = epsgCode;
	epoch_ = epoch;

	// set datum parameters using epsg code
	initialiseDatumFromEpsgCode();
}

void CDnaDatum::SetDatum(const std::string& epsgCode)
{
	// Reuse
	SetDatum(LongFromString<UINT32>(trimstr(epsgCode)));
}

//void CDnaDatum::SetDatum(const std::string& epsgCode, const boost::gregorian::date& epoch)
//{
//	// Reuse
//	SetDatum(LongFromString<UINT32>(trimstr(epsgCode)), epoch);
//}

void CDnaDatum::SetEpoch(const std::string& epoch)
{
	// Parse epoch
	if (epoch.empty())
	{
		if (isEpsgWGS84Ensemble(epsgCode_))
			epoch_ = timeImmemorial<boost::gregorian::date>();
		else
			epoch_ = dateFromString<boost::gregorian::date>(referenceepochFromEpsgCode(epsgCode_));
	}
	else
		epoch_ = dateFromString<boost::gregorian::date>(epoch);
}

//void CDnaDatum::SetEpoch(const double& decimal_year)
//{
//	epoch_ = dateFromDouble_doy_year<boost::gregorian::date, double>(decimal_year);
//
//}

void CDnaDatum::SetDatumFromName(const std::string& datumName, const std::string& epoch)
{
	// Get epsg code
	epsgCode_ = epsgCodeFromName<UINT32, std::string>(datumName);

	// Parse epoch
	SetEpoch(epoch);

	// set datum parameters using epsg code
	initialiseDatumFromEpsgCode();
}

void CDnaDatum::SetDatumFromEpsg(const std::string& epsgCode, const std::string& epoch)
{
	// Get epsg code
	epsgCode_ = LongFromString<UINT32>(trimstr(epsgCode));

	// Check valid epsg code (throws an exception on unknown code)
	validateEpsgCode(epsgCode_);

	// Parse epoch
	SetEpoch(epoch);

	// set datum parameters using epsg code
	initialiseDatumFromEpsgCode();
}

}	// namespace datum_parameters
}	// namespace dynadjust
