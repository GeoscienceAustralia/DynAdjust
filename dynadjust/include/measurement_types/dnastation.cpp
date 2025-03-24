//============================================================================
// Name         : dnastation.cpp
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
// Description  : DynAdjust CDnaStation implementation file
//============================================================================

#include <include/measurement_types/dnastation.hpp>
#include <include/parameters/dnaepsg.hpp>
#include <include/exception/dnaexception.hpp>

using namespace dynadjust::epsg;
using namespace dynadjust::datum_parameters;
using namespace dynadjust::exception;

namespace dynadjust {
namespace measurements {

// Comparison functions
// m_strName
bool operator<(const CDnaStation& left, const CDnaStation& right)
{
	return left.m_strName < right.m_strName;
}

bool operator<(const boost::shared_ptr<CDnaStation>& left, const boost::shared_ptr<CDnaStation>& right)
{
	return left.get()->m_strName < right.get()->m_strName;
}

CAStationList::CAStationList()
	: availMsrCount_ (0)
	, assocMsrCount_ (0)
	, amlStnIndex_ (0)
	, validStation_(VALID_STATION)
{
}

CAStationList::CAStationList(bool validStation)
	: availMsrCount_ (0)
	, assocMsrCount_ (0)
	, amlStnIndex_ (0)
	, validStation_(validStation ? VALID_STATION : INVALID_STATION)
{
}

CAStationList::~CAStationList()
{

}

// move constructor
CAStationList::CAStationList(const CAStationList&& s)
{
	availMsrCount_ = s.availMsrCount_;
	assocMsrCount_ = s.assocMsrCount_;
	amlStnIndex_ = s.amlStnIndex_;
	validStation_ = s.validStation_;
}


// move assignment operator
CAStationList& CAStationList::operator =(CAStationList&& rhs)
{
	// check for assignment to self!
	if (this == &rhs)
		return *this;

	availMsrCount_ = rhs.availMsrCount_;
	assocMsrCount_ = rhs.assocMsrCount_;
	amlStnIndex_ = rhs.amlStnIndex_;
	validStation_ = rhs.validStation_;

	return *this;
}















//<CDnaStation>//
CDnaStation::CDnaStation(const std::string& referenceframe, const std::string& epoch)
	: m_strName(""), m_strOriginalName("")
	, m_dXAxis(0.), m_dYAxis(0.), m_dZAxis(0.), m_dHeight(0.)
	, m_dStdDevX(0.), m_dStdDevY(0.), m_dStdDevZ(0.), m_dStdDevHt(0.)
	, m_strConstraints("FFF"), m_strType(""), m_strHemisphereZone(""), m_strDescription(""), m_strComment("")
	, m_cLatConstraint('F'), m_cLonConstraint('F'), m_cHtConstraint('F')
	, m_ctType(LLH_type_i), m_ctTypeSupplied(LLH_type_i), m_htType(ORTHOMETRIC_type_i)
	, m_dcurrentLatitude(0.), m_dcurrentLongitude(0.), m_dcurrentHeight(0.)
	, m_fgeoidSep(0.), m_dmeridianDef(0.), m_dverticalDef(0.)
	, m_lfileOrder(0), m_lnameOrder(0)
	, m_zone(0), m_unusedStation(INVALID_STATION)
	, m_referenceFrame(referenceframe), m_epoch(epoch)
	, m_constraintType(free_3D)
{
	m_epsgCode = epsgStringFromName<std::string>(referenceframe);
}

CDnaStation::~CDnaStation()
{
}

// copy constructors
CDnaStation::CDnaStation(const CDnaStation& newStation)
{
	m_strName = newStation.m_strName;
	m_strOriginalName = newStation.m_strOriginalName;
	m_strConstraints = newStation.m_strConstraints;
	m_strType = newStation.m_strType;
	m_dXAxis = newStation.m_dXAxis;
	m_dYAxis = newStation.m_dYAxis;
	m_dZAxis = newStation.m_dZAxis;
	m_dHeight = newStation.m_dHeight;
	m_dStdDevX = newStation.m_dStdDevX;
	m_dStdDevY = newStation.m_dStdDevY;
	m_dStdDevZ = newStation.m_dStdDevZ;
	m_dStdDevHt = newStation.m_dStdDevHt;
	m_strHemisphereZone = newStation.m_strHemisphereZone;
	m_strDescription = newStation.m_strDescription;
	m_cLatConstraint = newStation.m_cLatConstraint;
	m_cLonConstraint = newStation.m_cLonConstraint;
	m_cHtConstraint = newStation.m_cHtConstraint;
	m_strComment = newStation.m_strComment;

	m_ctTypeSupplied = newStation.m_ctTypeSupplied;
	m_ctType = newStation.m_ctType;
	m_htType = newStation.m_htType;

	m_dcurrentLatitude = newStation.m_dcurrentLatitude;
	m_dcurrentLongitude = newStation.m_dcurrentLongitude;
	m_dcurrentHeight = newStation.m_dcurrentHeight;
	m_fgeoidSep = newStation.m_fgeoidSep;
	m_dmeridianDef = newStation.m_dmeridianDef;
	m_dverticalDef = newStation.m_dverticalDef;
	m_lfileOrder = newStation.m_lfileOrder;
	m_lnameOrder = newStation.m_lnameOrder;
	m_zone = newStation.m_zone;
	m_unusedStation = newStation.m_unusedStation;

	m_referenceFrame = newStation.m_referenceFrame;
	m_epsgCode = newStation.m_epsgCode;
	m_epoch = newStation.m_epoch;

	m_constraintType = newStation.m_constraintType;
}

CDnaStation::CDnaStation(const std::string& strName, const std::string& strConstraints, 
	const std::string& strType, const double& dXAxis, const double& dYAxis, const double& dZAxis,
	const double& dHeight, const std::string& strHemisphereZone, const std::string& strDescription, 
	const std::string& strComment)
{
	m_strName = strName;

	m_strOriginalName = "";
	
	SetConstraints(strConstraints);
	SetCoordType(strType);

	m_dXAxis = dXAxis;
	m_dYAxis = dYAxis;
	m_dZAxis = dZAxis;
	m_dHeight = dHeight;

	m_dStdDevX = PRECISION_1E5;
	m_dStdDevY = PRECISION_1E5;
	m_dStdDevZ = PRECISION_1E5;
	m_dStdDevHt = PRECISION_1E4;

	if (strHemisphereZone.empty())
		m_zone = 0;
	else
		SetHemisphereZone(strHemisphereZone);

	m_strDescription = strDescription;
	m_strComment = strComment;
	m_unusedStation = INVALID_STATION;

	m_dcurrentLatitude = 0.;
	m_dcurrentLongitude = 0.;
	m_dcurrentHeight = 0.;
	m_fgeoidSep = 0.;
	m_dmeridianDef = 0.;
	m_dverticalDef = 0.;

	m_referenceFrame = DEFAULT_DATUM;
	m_epsgCode = DEFAULT_EPSG_S;
	m_epoch = "";
}

CDnaStation& CDnaStation::operator =(const CDnaStation& rhs)
{
	// check for assignment to self!
	if (this == &rhs)
		return *this;

	m_strName = rhs.m_strName;
	m_strConstraints = rhs.m_strConstraints;
	m_strType = rhs.m_strType;
	m_dXAxis = rhs.m_dXAxis;
	m_dYAxis = rhs.m_dYAxis;
	m_dZAxis = rhs.m_dZAxis;
	m_dHeight = rhs.m_dHeight;
	m_strHemisphereZone = rhs.m_strHemisphereZone;
	m_strDescription = rhs.m_strDescription;
	m_cLatConstraint = rhs.m_cLatConstraint;
	m_cLonConstraint = rhs.m_cLonConstraint;
	m_cHtConstraint = rhs.m_cHtConstraint;
	m_strComment = rhs.m_strComment;

	m_ctTypeSupplied = rhs.m_ctTypeSupplied;
	m_ctType = rhs.m_ctType;
	m_htType = rhs.m_htType;

	// measurement file extras
	m_dStdDevX = rhs.m_dStdDevX;
	m_dStdDevY = rhs.m_dStdDevY;
	m_dStdDevHt = rhs.m_dStdDevHt;

	m_dcurrentLatitude = rhs.m_dcurrentLatitude;
	m_dcurrentLongitude = rhs.m_dcurrentLongitude;
	m_dcurrentHeight = rhs.m_dcurrentHeight;
	m_fgeoidSep = rhs.m_fgeoidSep;
	m_dmeridianDef = rhs.m_dmeridianDef;
	m_dverticalDef = rhs.m_dverticalDef;
	m_lfileOrder = rhs.m_lfileOrder;
	m_lnameOrder = rhs.m_lnameOrder;
	m_zone = rhs.m_zone;

	m_unusedStation = rhs.m_unusedStation;

	m_referenceFrame = rhs.m_referenceFrame;
	m_epsgCode = rhs.m_epsgCode;
	m_epoch = rhs.m_epoch;

	m_constraintType = rhs.m_constraintType;
	
	return *this;
}


void CDnaStation::SetConstraints(const std::string& sConstraints)
{
	// capture string, trim whitespace
	m_strConstraints = trimstr(sConstraints);
	
	// No string provided?  Fill with FFF
	if (sConstraints.empty())
		m_strConstraints = "FFF";

	// Greater than 3 characters?  Trim to 3 characters
	if (m_strConstraints.length() > 3)
		m_strConstraints = m_strConstraints.substr(0, 3);
	
	// Less than 3 characters?  Pad with 'F'
	if (m_strConstraints.length() < 3)
		m_strConstraints.append(size_t(3 - m_strConstraints.length()), 'F');
	
	m_cLatConstraint = (char)(*m_strConstraints.substr(0, 1).c_str());
	m_cLonConstraint = (char)(*m_strConstraints.substr(1, 1).c_str());
	m_cHtConstraint = (char)(*m_strConstraints.substr(2, 1).c_str());

	// Free in all 3 dimensions
	if (boost::iequals(m_strConstraints, "FFF"))
		m_constraintType = free_3D;
	// Constrained in all 3 dimensions
	else if (boost::iequals(m_strConstraints, "CCC"))
		m_constraintType = constrained_3D;
	// Horizontal or 2D adjustment
	else if (boost::iequals(m_strConstraints, "FFC"))
		m_constraintType = free_2D;
	// Vertical or 1D adjustment
	else if (boost::iequals(m_strConstraints, "CCF"))
		m_constraintType = free_1D;
	else
		m_constraintType = custom_constraint;
}
	

void CDnaStation::SetCoordType(const std::string& sType) {
	m_strType = trimstr(sType);
	m_ctTypeSupplied = m_ctType = GetMyCoordTypeC();	
}

// X, Easting, Latitude
void CDnaStation::SetXAxis(const std::string& str)
{
	switch (m_ctType)
	{
	case LLH_type_i:
	case LLh_type_i:
		FromDmsString(&m_dXAxis, trimstr(str));
		m_dXAxis = Radians(m_dXAxis);
		break;
	default:
		// All other types will be converted by dna_import::ReduceStations_LLH()
		DoubleFromString(m_dXAxis, trimstr(str));
	}
}


// Y, Northing, Longitude
void CDnaStation::SetYAxis(const std::string& str)
{
	switch (m_ctType)
	{
	case LLH_type_i:
	case LLh_type_i:
		FromDmsString(&m_dYAxis, trimstr(str));
		m_dYAxis = Radians(m_dYAxis);
		break;
	default:
		DoubleFromString(m_dYAxis, trimstr(str));
	}
}


// Z
void CDnaStation::SetZAxis(const std::string& str)
{
	DoubleFromString(m_dZAxis, trimstr(str));
}


// Height
void CDnaStation::SetHeight(const std::string& str)
{
	if (GetMyCoordTypeC() == XYZ_type_i)
		SetZAxis(str);
	else
		DoubleFromString(m_dHeight, trimstr(str));
}

// Hemisphere zone
void CDnaStation::SetHemisphereZone(const std::string& sHemisphereZone)
{
	m_strHemisphereZone = trimstr(sHemisphereZone);
	if (m_strHemisphereZone.empty())
		return;
	if (isdigit(m_strHemisphereZone.substr(0).c_str()[0]))					// is the first character a decimal digit?
		m_zone = LongFromString<UINT32>(m_strHemisphereZone.substr(0));		// extract zone from first position
	else
		m_zone = LongFromString<UINT32>(m_strHemisphereZone.substr(1));		// extract zone after hemisphere
}

// Reduce from Cartesian or Projection to Geographic
// m_dXAxis always contains X, Easting and Latitude
// m_dYAxis always contains Y, Northing and Longitude
void CDnaStation::ReduceStations_LLH(const CDnaEllipsoid* m_eEllipsoid, const CDnaProjection* m_pProjection)
{
	switch (m_ctType)
	{
	case XYZ_type_i:
		// Convert from cartesian to geographic (radians)
		CartToGeo<double>(m_dXAxis, m_dYAxis, m_dZAxis, &m_dXAxis, &m_dYAxis, &m_dHeight, m_eEllipsoid);
		m_dZAxis= 0.;
		m_htType = ELLIPSOIDAL_type_i;

		// Force the current type to be geographic (ellipsoid height).  
		// The user-supplied type is retained in m_ctTypeSupplied
		m_ctType = LLh_type_i;
		m_strType = LLh_type;
		break;
	case UTM_type_i:
		// Convert from projection to geographic (radians)
		GridToGeo<double, UINT32>(m_dXAxis, m_dYAxis, m_zone, &m_dXAxis, &m_dYAxis, 
			m_eEllipsoid->GetSemiMajor(), m_eEllipsoid->GetInverseFlattening(),		// ellipsoid parameters
			m_pProjection->GetFalseEasting(), m_pProjection->GetFalseNorthing(),	// projection parameters
			m_pProjection->GetCentralScaleFactor(), m_pProjection->GetLongCentralMeridianZone1(), 
			m_pProjection->GetZoneWidth());

		// Force the current type to be geographic (orthometric height).
		// The user-supplied type is retained in m_ctTypeSupplied
		m_ctType = LLH_type_i;
		m_strType = LLH_type;
		break;
	default:
		break;
	}
}

// Reduce from Geographic or Projection to Cartesian
// m_dXAxis always contains X, Easting and Latitude
// m_dYAxis always contains Y, Northing and Longitude
void CDnaStation::ReduceStations_XYZ(const CDnaEllipsoid* m_eEllipsoid, const CDnaProjection* m_pProjection)
{
	switch (m_ctType)
	{
	case UTM_type_i:
		// Convert from projection to geographic (radians)
		GridToGeo<double, UINT32>(m_dXAxis, m_dYAxis, m_zone, &m_dXAxis, &m_dYAxis, 
			m_eEllipsoid->GetSemiMajor(), m_eEllipsoid->GetInverseFlattening(),	// ellipsoid parameters
			m_pProjection->GetFalseEasting(), m_pProjection->GetFalseNorthing(),	// projection parameters
			m_pProjection->GetCentralScaleFactor(), m_pProjection->GetLongCentralMeridianZone1(), 
			m_pProjection->GetZoneWidth());
		[[fallthrough]];
	case LLH_type_i:
	case LLh_type_i:
		m_dcurrentLatitude = m_dXAxis;
		m_dcurrentLongitude = m_dYAxis;
		m_dcurrentHeight = m_dHeight;
		if (m_ctType == LLH_type_i)
			m_dcurrentHeight += m_fgeoidSep;
		// Convert from geographic (radians) to cartesian
		// Assumes height is ellipsoidal
		GeoToCart<double>(m_dcurrentLatitude, m_dcurrentLongitude, m_dcurrentHeight, &m_dXAxis, &m_dYAxis, &m_dZAxis, m_eEllipsoid);
		break;
	default:
		break;
	}

	// Force the current type to be geographic.  The user-supplied type is retained in m_ctTypeSupplied
	m_ctType = XYZ_type_i;
	m_strType = XYZ_type;
}

bool CDnaStation::IsValidConstraint(const std::string& sConst)
{
	if (boost::iequals(sConst, "CCC"))
		return true;
	if (boost::iequals(sConst, "CCF"))
		return true;
	if (boost::iequals(sConst, "CFF"))
		return true;
	if (boost::iequals(sConst, "FFF"))
		return true;
	if (boost::iequals(sConst, "FFC"))
		return true;
	if (boost::iequals(sConst, "FCC"))
		return true;
	if (boost::iequals(sConst, "FCF"))
		return true;
	if (boost::iequals(sConst, "CFC"))
		return true;
	return false;
}
	

_COORD_TYPE_ CDnaStation::GetCoordTypeC(const std::string& sType)
{
	// case insensitive
	if (boost::iequals(sType, XYZ_type))
		return XYZ_type_i;
	else if (boost::iequals(sType, UTM_type))
		return UTM_type_i;				// height is assumed to be orthometric
	else if (boost::iequals(sType, ENU_type))
		return ENU_type_i;

	// case sensitive
	else if (boost::equals(sType, LLh_type))
		return LLh_type_i;				// ellipsoid height
	
	// default
	else if (boost::equals(sType, LLH_type))
		return LLH_type_i;					// orthometric height (default)

	// If this point is reached, sType is an unknown coordinate type, so throw!
	std::stringstream ss;
	ss << "  '" << sType << "' is not a recognised coordinate type." << std::endl;
	throw boost::enable_current_exception(std::runtime_error(ss.str()));

	return LLH_type_i;
}


_COORD_TYPE_ CDnaStation::GetMyCoordTypeC() const
{
	return GetCoordTypeC(m_strType);
}


_HEIGHT_SYSTEM_ CDnaStation::GetHeightSystemC(const std::string& sType) const
{
	if (sType.compare(ELLIPSOIDAL_type) == 0)
		return ELLIPSOIDAL_type_i;
	return ORTHOMETRIC_type_i;		// default
}


_HEIGHT_SYSTEM_ CDnaStation::GetMyHeightSystemC() const
{
	if (m_strType.compare(LLh_type) == 0)
		return ELLIPSOIDAL_type_i;
	return ORTHOMETRIC_type_i;	// default
}

// SetHeightSystem called by void Height_pimpl::system(const ::std::string& system)
// where system is an attribute of the element Height (either "ellipsoidal" or "orthometric")
void CDnaStation::SetHeightSystem(const std::string& sType)
{
	SetHeightSystem(GetHeightSystemC(sType));
}
	
void CDnaStation::SetHeightSystem(const HEIGHT_SYSTEM& type_i)
{
	m_htType = type_i;

	switch (m_htType)
	{
	case ELLIPSOIDAL_type_i:
		// Convert to orthometric
		if (m_ctType == LLH_type_i)
			m_ctType = LLh_type_i;
		if (m_strType.compare(LLH_type))
			m_strType = LLh_type;
		break;
	case ORTHOMETRIC_type_i:
		// Convert to ellipsoidal
		if (m_ctType == LLh_type_i)
			m_ctType = LLH_type_i;
		if (m_strType.compare(LLh_type))
			m_strType = LLH_type;
		break;
	}
}

//void CDnaStation::coutStationData(std::ostream &os, ostream &os2, const UINT16& uType) const
//{
//	UINT32 precision = 3;
//	if (m_ctType == LLH_type_i)
//		precision = 10;
//	std::stringstream ss;
//	std::string str;
//	size_t dot;
//
//	switch (uType)
//	{
//	case DNA_COUT:
//	case GEOLAB_COUT:
//		
//		os << "+ " << std::setw(16) << m_strName;
//		os << std::setw(4) << m_strConstraints;
//		os << std::setw(4) << m_strType;
//		os << std::setw(20) << std::setprecision(precision) << std::fixed << (m_ctType == LLH_type_i ? RadtoDms(m_dXAxis): m_dXAxis);
//		os << std::setw(20) << std::setprecision(precision) << std::fixed << (m_ctType == LLH_type_i ? RadtoDmsL(m_dYAxis): m_dYAxis);
//		//os << std::setw(16) << vStations[s].GetZAxis();
//		os << std::setw(13) << std::setprecision(4) << std::fixed << m_dHeight;
//		//os << std::setw(10) << vStations[s].GetRedHeight();
//		os << std::setw(5) << m_strHemisphereZone;
//		os << std::setw(6) << m_lnameOrder;
//		os << std::setw(19) << m_strDescription;
//		//os << std::setw(10) << vStations[s].GetComment();
//		os << std::endl;
//		break;
//	case NEWGAN_COUT:
//		os << std::setw(3) << "4";
//		os << std::right << std::setw(12) << m_strName;
//		os << " ";
//		os << std::left << std::setw(23) << m_strDescription.substr(0, 23);
//		ss.str("");
//		ss << std::setprecision(precision) << std::fixed << (m_ctType == LLH_type_i ? RadtoDms(m_dXAxis): m_dXAxis);
//		str = trimstr(ss.str());
//		dot = str.find(".");
//		str.replace(dot, 1, " ");
//		str.insert(dot+5, ".");
//		if (m_dXAxis < 0)
//		{
//			str.replace(0, 1, " ");		// replace negative sign
//			str = trimstr(str);
//			os << std::left << "S" << std::setw(precision + 4) << std::right << str;
//		}
//		else
//			os << std::left << "N" << std::setw(precision + 4) << std::right << str;
//		
//		ss.str("");
//		ss << std::setprecision(precision) << std::fixed << (m_ctType == LLH_type_i ? RadtoDmsL(m_dYAxis): m_dYAxis);
//		str = trimstr(ss.str());
//		dot = str.find(".");
//		str.replace(dot, 1, " ");
//		str.insert(dot+5, ".");
//		if (m_dYAxis < 0)
//		{
//			str.replace(0, 1, " ");		// replace negative sign
//			str = trimstr(str);
//			os << std::left << "W" << std::setw(precision + 5) << std::right << str;
//		}
//		else
//			os << std::left << "E" << std::setw(precision + 5) << std::right << str;
//		
//		
//		os << std::right << std::setw(9) << std::setprecision(3) << std::fixed << m_dHeight;
//		os << std::endl;
//		break;
//	case GMT_OUT:
//		os << std::setprecision(precision) << std::fixed << (m_ctType == LLH_type_i ? DegreesL(m_dYAxis): m_dYAxis);
//		os << "  ";
//		os << std::setprecision(precision) << std::fixed << (m_ctType == LLH_type_i ? Degrees(m_dXAxis): m_dXAxis);
//		os << std::endl;
//
//		os2 << std::setprecision(precision) << std::fixed << (m_ctType == LLH_type_i ? DegreesL(m_dYAxis): m_dYAxis);
//		os2 << "  ";
//		os2 << std::setprecision(precision) << std::fixed << (m_ctType == LLH_type_i ? Degrees(m_dXAxis): m_dXAxis);
//		os2 << "  6 0 0 LM " << m_strName << std::endl;		
//		break;
//	}
//	
//
//}


void CDnaStation::WriteBinaryStn(std::ofstream* binary_stream, const UINT16 bUnused)
{
	station_t stationRecord;
	strcpy(stationRecord.stationName, m_strName.substr(0, STN_NAME_WIDTH).c_str());
	strcpy(stationRecord.stationNameOrig, m_strOriginalName.substr(0, STN_NAME_WIDTH).c_str());
	strcpy(stationRecord.stationConst, m_strConstraints.substr(0, STN_CONST_WIDTH).c_str());
	strcpy(stationRecord.stationType, m_strType.substr(0, STN_TYPE_WIDTH).c_str());

	stationRecord.suppliedStationType = m_ctTypeSupplied;

	switch (m_ctTypeSupplied)
	{
	case LLH_type_i:
	case UTM_type_i:
		stationRecord.suppliedHeightRefFrame = ORTHOMETRIC_type_i;
		break;
	default:
		break;
	}

	stationRecord.initialLatitude = m_dXAxis;
	if (fabs(m_dcurrentLatitude - 0.) < PRECISION_1E15)
		stationRecord.currentLatitude = m_dXAxis;
	else
		stationRecord.currentLatitude = m_dcurrentLatitude;
	
	stationRecord.initialLongitude = m_dYAxis;
	if (fabs(m_dcurrentLongitude - 0.) < PRECISION_1E15)
		stationRecord.currentLongitude = m_dYAxis;
	else
		stationRecord.currentLongitude = m_dcurrentLongitude;

	stationRecord.initialHeight = m_dHeight;	
	if (fabs(m_dcurrentHeight - 0.) < PRECISION_1E15)
		stationRecord.currentHeight = m_dHeight;
	else
		stationRecord.currentHeight = m_dcurrentHeight;
		
	stationRecord.geoidSep = m_fgeoidSep;
	stationRecord.meridianDef = m_dmeridianDef;
	stationRecord.verticalDef = m_dverticalDef;
	stationRecord.zone = m_zone;
	strcpy(stationRecord.description, m_strDescription.substr(0, STN_DESC_WIDTH).c_str());
	stationRecord.fileOrder = m_lfileOrder;
	stationRecord.nameOrder = m_lnameOrder;
	stationRecord.unusedStation = bUnused;

	strcpy(stationRecord.epoch, m_epoch.substr(0, STN_EPOCH_WIDTH).c_str());
	strcpy(stationRecord.epsgCode, m_epsgCode.substr(0, STN_EPSG_WIDTH).c_str());

	binary_stream->write(reinterpret_cast<char *>(&stationRecord), sizeof(station_t));
}


void CDnaStation::WriteDNAXMLStnCurrentEstimates(std::ofstream* dna_ofstream, 
	const CDnaEllipsoid* ellipsoid, const CDnaProjection* projection,
	INPUT_FILE_TYPE t, const dna_stn_fields* dsw)
{
	// Get the latest estimates
	double lat_east_x(m_dcurrentLatitude);
	double lon_north_y(m_dcurrentLongitude);
	double ht_zone_z(m_dcurrentHeight);

	// Write latest station estimates
	WriteDNAXMLStn(dna_ofstream,
		ellipsoid, projection,
		lat_east_x, lon_north_y, ht_zone_z,
		t, dsw);
}
	

void CDnaStation::WriteDNAXMLStnInitialEstimates(std::ofstream* dna_ofstream, 
	const CDnaEllipsoid* ellipsoid, const CDnaProjection* projection,
	INPUT_FILE_TYPE t, const dna_stn_fields* dsw)
{
	// Write initial station estimates
	WriteDNAXMLStn(dna_ofstream,
		ellipsoid, projection,
		m_dXAxis, m_dYAxis, m_dHeight,
		t, dsw);
}


void CDnaStation::WriteDNAXMLStn(std::ofstream* dna_ofstream,
	const CDnaEllipsoid* ellipsoid, const CDnaProjection* projection,
	const double& lat_east_x, const double& lon_north_y, const double& ht_zone_z,
	INPUT_FILE_TYPE t, const dna_stn_fields* dsw)
{
	//m_ctTypeSupplied = LLh_type_i;

	std::string hemisphereZone(m_strHemisphereZone);
	std::string coordinateType(LLH_type);

	switch (m_ctTypeSupplied)
	{
	case XYZ_type_i:
		coordinateType = XYZ_type;
		break;
	case LLh_type_i:
		coordinateType = LLh_type;
		break;
	case UTM_type_i:
		coordinateType = UTM_type;
		break;
	case ENU_type_i:	// Not supported yet, so revert to Lat, Long, Height
		coordinateType = LLH_type;
		m_ctTypeSupplied = LLH_type_i;
		break;
	case LLH_type_i:
		coordinateType = LLH_type;
		break;
	default:
		break;
	}

	// Convert height to orthometric?
	// Decisions for printing (estimated) station height for various input coordinate types
	// and geoid model.  Estimated stations are in LLh.
	// h = H + Na
	//
	// Input Coord type		Geoid loaded?          Geoid not loaded?          Result
	// ---------------------------------------------------------------------------------
	// LLH				    ht_z - m_fgeoidSep     ht_z (m_fgeoidSep = 0)     H (orthometric)
	// UTM					ht_z - m_fgeoidSep     ht_z (m_fgeoidSep = 0)     H (orthometric)
	// ENU					ht_z - m_fgeoidSep     ht_z (m_fgeoidSep = 0)     H (orthometric)
	// LLh          	    ht_z			       ht_z					      h (ellipsoid)
	// XYZ          	    LLh->XYZ		       LLh->XYZ					  X, Y, Z
	double lat_east_x_mod(lat_east_x);
	double lon_north_y_mod(lon_north_y);
	double ht_zone_z_mod(ht_zone_z);

	switch (m_ctTypeSupplied)
	{
	case LLH_type_i:
	case UTM_type_i:
		// reduce to orthometric.  
		// If geoid hasn't been loaded, m_fgeoidSep will be 0.0
		ht_zone_z_mod -= m_fgeoidSep;
		break;
	default:
		break;
	}

	// Convert coordinates to cartesian or utm?
	// Convert radians values to degrees, minutes and seconds?
	std::stringstream ss;
	double zone;
	switch (m_ctTypeSupplied)
	{
	case XYZ_type_i:
		GeoToCart<double>(lat_east_x, lon_north_y, ht_zone_z, 
			&lat_east_x_mod, &lon_north_y_mod, &ht_zone_z_mod, ellipsoid);
		break;
	
	case UTM_type_i:
		GeoToGrid<double>(lat_east_x, lon_north_y, 
			&lat_east_x_mod, &lon_north_y_mod, &zone, ellipsoid, projection, true);
		
		ss << std::fixed << std::setprecision(0) << zone;
		hemisphereZone = ss.str();	
		break;

	case LLh_type_i:
	case LLH_type_i:
		lat_east_x_mod = RadtoDms(lat_east_x);
		lon_north_y_mod = RadtoDmsL(lon_north_y);
		break;
	default:
		break;
	}

	switch (t)
	{
	case dna:
	
		// Write DNA station file
		WriteDNAStn(dna_ofstream, coordinateType,
			lat_east_x_mod, lon_north_y_mod, ht_zone_z_mod, 
			hemisphereZone, *dsw);
		break;

	case dynaml:

		// Write DynAdjust XML station file
		WriteDynaMLStn(dna_ofstream, coordinateType,
			lat_east_x_mod, lon_north_y_mod, ht_zone_z_mod,
			hemisphereZone);
		break;
	default:
		break;
	}
}

void CDnaStation::WriteDNAStn(std::ofstream* dna_ofstream, const std::string& coordinateType, 
	const double& lat_east_x, const double& lon_north_y, const double& ht_zone_z,
	std::string& hemisphereZone, const dna_stn_fields& dsw)
{	
	UINT32 LEX_precision(4), LNY_precision(4), HZ_precision(4);	
	
	switch (m_ctTypeSupplied)
	{
	case LLh_type_i:
	case LLH_type_i:
		LEX_precision = LNY_precision = 10;
		break;
	default:
		break;
	}

	(*dna_ofstream) << std::left << std::setw(dsw.stn_name) << m_strName;
	(*dna_ofstream) << std::left << std::setw(dsw.stn_const) << m_strConstraints;
	(*dna_ofstream) << " ";
	(*dna_ofstream) << std::left << std::setw(dsw.stn_type) << coordinateType;
	(*dna_ofstream) << std::right << std::setw(dsw.stn_e_phi_x) << std::setprecision(LEX_precision) << std::fixed << lat_east_x;
	(*dna_ofstream) << std::right << std::setw(dsw.stn_n_lam_y) << std::setprecision(LNY_precision) << std::fixed << lon_north_y;
	(*dna_ofstream) << std::right << std::setw(dsw.stn_ht_z) << std::setprecision(HZ_precision) << std::fixed << ht_zone_z;	
	if (m_ctTypeSupplied == UTM_type_i)
		(*dna_ofstream) << std::right << std::setw(dsw.stn_hemi_zo) << hemisphereZone;
	else
		(*dna_ofstream) << std::right << std::setw(dsw.stn_hemi_zo) << " ";
	(*dna_ofstream) << " ";
	(*dna_ofstream) << std::left << m_strDescription << std::endl;
}
	

void CDnaStation::WriteDynaMLStn(std::ofstream* xml_ofstream, const std::string& coordinateType, 
	const double& lat_east_x, const double& lon_north_y, const double& ht_zone_z,
	std::string& hemisphereZone)
{
	UINT32 LEX_precision(4), LNY_precision(4), HZ_precision(4);
	
	// Convert radians values to degrees, minutes and seconds
	switch (m_ctTypeSupplied)
	{
	case LLh_type_i:
	case LLH_type_i:
		LEX_precision = LNY_precision = 10;
		break;
	default:
		break;
	}
	
	(*xml_ofstream) << "  <DnaStation>" << std::endl;
	(*xml_ofstream) << "    <Name>" << m_strName << "</Name>" << std::endl;
	(*xml_ofstream) << "    <Constraints>" << m_strConstraints << "</Constraints>" << std::endl;
	(*xml_ofstream) << "    <Type>" << coordinateType << "</Type>" << std::endl;
	(*xml_ofstream) << "    <StationCoord>" << std::endl;
	(*xml_ofstream) << "      <Name>" << m_strName << "</Name>" << std::endl;
	(*xml_ofstream) << "      <XAxis>" << std::setprecision(LEX_precision) << std::fixed << lat_east_x << "</XAxis>" << std::endl;
	(*xml_ofstream) << "      <YAxis>" << std::setprecision(LNY_precision) << std::fixed << lon_north_y << "</YAxis>" << std::endl;
	(*xml_ofstream) << "      <Height>" << std::setprecision(HZ_precision) << std::fixed << ht_zone_z << "</Height>" << std::endl;
	
	// Convert radians values to degrees, minutes and seconds
	if (m_ctTypeSupplied == UTM_type_i)
		(*xml_ofstream) << "      <HemisphereZone>" << hemisphereZone << "</HemisphereZone>" << std::endl;
	
	(*xml_ofstream) << "    </StationCoord>" << std::endl;
	(*xml_ofstream) << "    <Description>"  << m_strDescription << "</Description>" << std::endl;
	(*xml_ofstream) << "  </DnaStation>" << std::endl;
}
	

void CDnaStation::WriteGeoidfile(std::ofstream* geo_ofstream)
{
	(*geo_ofstream) << std::setw(44) << std::left << m_strName <<
		std::setw(15) << std::setprecision(3) << std::fixed << std::left << m_fgeoidSep <<
		std::setw(10) << std::setprecision(3) << std::fixed << std::right << m_dmeridianDef <<
		std::setw(10) << std::setprecision(3) << std::fixed << std::right << m_dverticalDef << std::endl;
}
	

void CDnaStation::SetStationRec(const station_t& stationRecord)
{
	SetName(stationRecord.stationName);
	SetOriginalName(stationRecord.stationNameOrig);
	SetConstraints(stationRecord.stationConst);
	SetCoordType(stationRecord.stationType);
	
	switch (stationRecord.suppliedStationType)
	{
	case UTM_type_i:
		m_ctTypeSupplied = UTM_type_i;
		break;
	case XYZ_type_i:
		m_ctTypeSupplied = XYZ_type_i;
		break;
	case ENU_type_i:
		m_ctTypeSupplied = ENU_type_i;
		break;
	case LLh_type_i:
		m_ctTypeSupplied = LLh_type_i;
		break;
	default:
	case LLH_type_i:
		m_ctTypeSupplied = LLH_type_i;
		break;
	}
	
	m_dXAxis = stationRecord.initialLatitude;
	m_dYAxis = stationRecord.initialLongitude;
	m_dHeight = stationRecord.initialHeight;
	m_dcurrentLatitude = stationRecord.currentLatitude;
	m_dcurrentLongitude = stationRecord.currentLongitude;
	m_dcurrentHeight = stationRecord.currentHeight;
	m_fgeoidSep = stationRecord.geoidSep;
	m_dmeridianDef = stationRecord.meridianDef;
	m_dverticalDef = stationRecord.verticalDef;
	m_zone = stationRecord.zone;
	if (stationRecord.initialLatitude < 0.)
		m_strHemisphereZone = "S";
	else
		m_strHemisphereZone = "N";
	char zone[3];
	snprintf(zone, sizeof(zone), "%d", m_zone);
	m_strHemisphereZone += zone;
	SetDescription(stationRecord.description);
	m_lfileOrder = stationRecord.fileOrder;
	m_lnameOrder = stationRecord.nameOrder;
	m_unusedStation = (stationRecord.unusedStation == VALID_STATION ? true : false);

	m_epoch = stationRecord.epoch;
	m_epsgCode = stationRecord.epsgCode;
	m_referenceFrame = datumFromEpsgCode<std::string, UINT32>(LongFromString<UINT32>(m_epsgCode));
}

std::string CDnaStation::CoordinateName(const char& c)
{
	switch (c)
	{
	case 'P':
		return "Latitude";
	case 'L':
		return "Longitude";
	case 'H':
		return "H(Ortho)";
	case 'h':
		return "h(Ellipse)";
	case 'E':
		return "Easting";
	case 'N':
		return "Northing";
	case 'z':
		return "Zone";
	case 'X':
		return "X";
	case 'Y':
		return "Y";
	case 'Z':
		return "Z";
	}
	return "";
}

}	// namespace measurements
}	// namespace dynadjust
