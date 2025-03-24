//============================================================================
// Name         : dnameasurement.cpp
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
// Description  : CDnaMeasurement and CDnaCovariance implementation file
//============================================================================

#include <include/exception/dnaexception.hpp>
#include <include/measurement_types/dnameasurement.hpp>

using namespace dynadjust::exception;
using namespace dynadjust::math;

namespace dynadjust {
namespace measurements {

CDnaCovariance::CDnaCovariance(void)
	: m_bIgnore(false)
	, m_lstn1Index(0)
	, m_lstn2Index(0)
	, m_strType("")
	, m_dM11(0.)
	, m_dM12(0.)
	, m_dM13(0.)
	, m_dM21(0.)
	, m_dM22(0.)
	, m_dM23(0.)
	, m_dM31(0.)
	, m_dM32(0.)
	, m_dM33(0.)
	, m_lclusterID(0)
{
}


CDnaCovariance::~CDnaCovariance(void)
{

}

// move constructor
CDnaCovariance::CDnaCovariance(CDnaCovariance&& c)
{
	m_bIgnore = c.m_bIgnore;
	m_lstn1Index = c.m_lstn1Index;
	m_lstn2Index = c.m_lstn2Index;
	m_strType = c.m_strType;
	m_dM11 = c.m_dM11;
	m_dM12 = c.m_dM12;
	m_dM13 = c.m_dM13;
	m_dM21 = c.m_dM21;
	m_dM22 = c.m_dM22;
	m_dM23 = c.m_dM23;
	m_dM31 = c.m_dM31;
	m_dM32 = c.m_dM32;
	m_dM33 = c.m_dM33;
	m_lclusterID = c.m_lclusterID;
}


// move assignment operator
CDnaCovariance& CDnaCovariance::operator= (CDnaCovariance&& rhs)
{
	// check for assignment to self!
	if (this == &rhs)
		return *this;

	m_bIgnore = rhs.m_bIgnore;
	m_lstn1Index = rhs.m_lstn1Index;
	m_lstn2Index = rhs.m_lstn2Index;
	m_strType = rhs.m_strType;
	m_dM11 = rhs.m_dM11;
	m_dM12 = rhs.m_dM12;
	m_dM13 = rhs.m_dM13;
	m_dM21 = rhs.m_dM21;
	m_dM22 = rhs.m_dM22;
	m_dM23 = rhs.m_dM23;
	m_dM31 = rhs.m_dM31;
	m_dM32 = rhs.m_dM32;
	m_dM33 = rhs.m_dM33;
	m_lclusterID = rhs.m_lclusterID;

	return *this;
}
	

bool CDnaCovariance::operator== (const CDnaCovariance& rhs) const
{
	return (
		m_bIgnore == rhs.m_bIgnore &&
		m_lstn1Index == rhs.m_lstn1Index &&
		m_lstn2Index == rhs.m_lstn2Index &&
		m_strType == rhs.m_strType &&
		m_dM11 == rhs.m_dM11 &&
		m_dM12 == rhs.m_dM12 &&
		m_dM13 == rhs.m_dM13 &&
		m_dM21 == rhs.m_dM21 &&
		m_dM22 == rhs.m_dM22 &&
		m_dM23 == rhs.m_dM23 &&
		m_dM31 == rhs.m_dM31 &&
		m_dM32 == rhs.m_dM32 &&
		m_dM33 == rhs.m_dM33
	);
}	

void CDnaCovariance::WriteDynaMLMsr(std::ofstream* dynaml_stream) const
{
	if (GetTypeC() == 'X')
		*dynaml_stream << "      <GPSCovariance>" << std::endl;
	else
		*dynaml_stream << "      <PointCovariance>" << std::endl;
		
	*dynaml_stream << "        <m11>" << std::scientific << std::setprecision(13) << m_dM11 << "</m11>" << std::endl;
	*dynaml_stream << "        <m12>" << m_dM12 << "</m12>" << std::endl;
	*dynaml_stream << "        <m13>" << m_dM13 << "</m13>" << std::endl;
	*dynaml_stream << "        <m21>" << m_dM21 << "</m21>" << std::endl;
	*dynaml_stream << "        <m22>" << m_dM22 << "</m22>" << std::endl;
	*dynaml_stream << "        <m23>" << m_dM23 << "</m23>" << std::endl;
	*dynaml_stream << "        <m31>" << m_dM31 << "</m31>" << std::endl;
	*dynaml_stream << "        <m32>" << m_dM32 << "</m32>" << std::endl;
	*dynaml_stream << "        <m33>" << m_dM33 << "</m33>" << std::endl;
	
	if (GetTypeC() == 'X')
		*dynaml_stream << "      </GPSCovariance>" << std::endl;
	else
		*dynaml_stream << "      </PointCovariance>" << std::endl;
}
	

void CDnaCovariance::WriteDNAMsr(std::ofstream* dna_stream, 
	const dna_msr_fields& dmw, const dna_msr_fields&) const
{
	UINT32 pad(dmw.msr_type + dmw.msr_ignore + dmw.msr_inst + dmw.msr_targ1 + dmw.msr_targ2 + dmw.msr_gps);
	// X
	*dna_stream << std::setw(pad) << " ";
	*dna_stream <<
		std::right << std::setw(dmw.msr_gps_vcv_1) << std::scientific << std::setprecision(13) << m_dM11 <<
		std::right << std::setw(dmw.msr_gps_vcv_2) << m_dM12 <<
		std::right << std::setw(dmw.msr_gps_vcv_3) << m_dM13;
	*dna_stream << std::endl;
		
	// Y
	*dna_stream << std::setw(pad) << " ";
	*dna_stream << 
		std::right << std::setw(dmw.msr_gps_vcv_1) << std::scientific << std::setprecision(13) << m_dM21 <<
		std::right << std::setw(dmw.msr_gps_vcv_2) << m_dM22 <<
		std::right << std::setw(dmw.msr_gps_vcv_3) << m_dM23;
	*dna_stream << std::endl;

	// Z
	*dna_stream << std::setw(pad) << " ";
	*dna_stream << 
		std::right << std::setw(dmw.msr_gps_vcv_1) << std::scientific << std::setprecision(13) << m_dM31 <<
		std::right << std::setw(dmw.msr_gps_vcv_2) << m_dM32 <<
		std::right << std::setw(dmw.msr_gps_vcv_3) << m_dM33;
	*dna_stream << std::endl;
}
	

void CDnaCovariance::SimulateMsr(vdnaStnPtr*, const CDnaEllipsoid*)
{
	m_dM11 = 0.0;
	m_dM12 = 0.0;
	m_dM13 = 0.0;
	m_dM21 = 0.0;
	m_dM22 = 0.0;
	m_dM23 = 0.0;
	m_dM31 = 0.0;
	m_dM32 = 0.0;
	m_dM33 = 0.0;
}
	

UINT32 CDnaCovariance::SetMeasurementRec(const vstn_t&, it_vmsr_t& it_msr)
{
	// get data relating to Z, sigmaZZ, sigmaXZ, sigmaYZ
	it_msr++;
	m_lstn1Index = it_msr->station1;
	m_lstn2Index = it_msr->station2;
	m_strType = it_msr->measType;
	m_dM11 = it_msr->term1;
	m_dM12 = it_msr->term2;
	m_dM13 = it_msr->term3;

	m_lclusterID = it_msr->clusterID;

	// get data relating to Z, sigmaZZ, sigmaXZ, sigmaYZ
	it_msr++;
	
	m_lstn1Index = it_msr->station1;
	m_lstn2Index = it_msr->station2;
	m_dM21 = it_msr->term1;
	m_dM22 = it_msr->term2;
	m_dM23 = it_msr->term3;

	// get data relating to Z, sigmaZZ, sigmaXZ, sigmaYZ
	it_msr++;
	
	m_lstn1Index = it_msr->station1;
	m_lstn2Index = it_msr->station2;
	m_dM31 = it_msr->term1;
	m_dM32 = it_msr->term2;
	m_dM33 = it_msr->term3;
	
	return 0;
}
	

void CDnaCovariance::WriteBinaryMsr(std::ofstream *binary_stream, PUINT32 msrIndex, const std::string& epsgCode, const std::string& epoch) const
{
	*msrIndex += 3;
	measurement_t measRecord;

	// Common
	measRecord.measType = GetTypeC();
	measRecord.station1 = m_lstn1Index;	
	measRecord.station2 = m_lstn2Index;	
	measRecord.clusterID = m_lclusterID;
	
	snprintf(measRecord.epsgCode, sizeof(measRecord.epsgCode), "%s", epsgCode.substr(0, STN_EPSG_WIDTH).c_str());
	snprintf(measRecord.epoch, sizeof(measRecord.epoch), "%s", epoch.substr(0, STN_EPOCH_WIDTH).c_str());

	// X
	measRecord.measStart = xCov;
	measRecord.term1 = m_dM11;
	measRecord.term2 = m_dM12;
	measRecord.term3 = m_dM13;

	binary_stream->write(reinterpret_cast<char *>(&measRecord), sizeof(measurement_t));

	// Y
	measRecord.measStart = yCov;
	measRecord.term1 = m_dM21;
	measRecord.term2 = m_dM22;
	measRecord.term3 = m_dM23;

	binary_stream->write(reinterpret_cast<char *>(&measRecord), sizeof(measurement_t));

	// Z
	measRecord.measStart = zCov;
	measRecord.term1 = m_dM31;
	measRecord.term2 = m_dM32;
	measRecord.term3 = m_dM33;
	
	binary_stream->write(reinterpret_cast<char *>(&measRecord), sizeof(measurement_t));
}

void CDnaCovariance::SerialiseDatabaseMap(std::ofstream* os, const msr_database_id_map& dbid)
{
	UINT16 msr, cls;
	msr = static_cast<UINT16>(dbid.is_msr_id_set);
	cls = static_cast<UINT16>(dbid.is_cls_id_set);
	
	// X
	os->write(reinterpret_cast<const char *>(&dbid.msr_id), sizeof(UINT32));
	os->write(reinterpret_cast<const char *>(&dbid.cluster_id), sizeof(UINT32));
	os->write(reinterpret_cast<const char *>(&msr), sizeof(UINT16));
	os->write(reinterpret_cast<const char *>(&cls), sizeof(UINT16));
	// Y
	os->write(reinterpret_cast<const char *>(&dbid.msr_id), sizeof(UINT32));
	os->write(reinterpret_cast<const char *>(&dbid.cluster_id), sizeof(UINT32));
	os->write(reinterpret_cast<const char *>(&msr), sizeof(UINT16));
	os->write(reinterpret_cast<const char *>(&cls), sizeof(UINT16));
	// Z
	os->write(reinterpret_cast<const char *>(&dbid.msr_id), sizeof(UINT32));
	os->write(reinterpret_cast<const char *>(&dbid.cluster_id), sizeof(UINT32));
	os->write(reinterpret_cast<const char *>(&msr), sizeof(UINT16));
	os->write(reinterpret_cast<const char *>(&cls), sizeof(UINT16));
}

void CDnaCovariance::SetM11(const std::string& str)
{
	DoubleFromString(m_dM11, trimstr(str));
}

void CDnaCovariance::SetM12(const std::string& str)
{
	DoubleFromString(m_dM12, trimstr(str));
}

void CDnaCovariance::SetM13(const std::string& str)
{
	DoubleFromString(m_dM13, trimstr(str));
}

void CDnaCovariance::SetM21(const std::string& str)
{
	DoubleFromString(m_dM21, trimstr(str));
}

void CDnaCovariance::SetM22(const std::string& str)
{
	DoubleFromString(m_dM22, trimstr(str));
}

void CDnaCovariance::SetM23(const std::string& str)
{
	DoubleFromString(m_dM23, trimstr(str));
}

void CDnaCovariance::SetM31(const std::string& str)
{
	DoubleFromString(m_dM31, trimstr(str));
}

void CDnaCovariance::SetM32(const std::string& str)
{
	DoubleFromString(m_dM32, trimstr(str));
}

void CDnaCovariance::SetM33(const std::string& str)
{
	DoubleFromString(m_dM33, trimstr(str));
}



// Comparison functions
// m_strFirst
bool operator<(const CDnaMeasurement& left, const CDnaMeasurement& right)
{
	return left.m_strFirst < right.m_strFirst;
}

bool operator<(const dnaMsrPtr& left, const dnaMsrPtr& right)
{
	return left.get()->m_strFirst < right.get()->m_strFirst;
}



CDnaMeasurement::CDnaMeasurement()
	: m_strFirst("")
	, m_MSmeasurementStations(ONE_STATION)
	, m_bIgnore(false)
	, m_lmeasurementIndex(0)
	, m_lstn1Index(0)
	, m_lstn2Index(0)
	, m_lstn3Index(0)
	, m_measAdj(0)
	, m_measCorr(0)
	, m_measAdjPrec(0)
	, m_residualPrec(0)
	, m_preAdjCorr(0)
	, m_epsgCode(DEFAULT_EPSG_S)
	, m_epoch("")
	, m_bInsufficient(false)
{
}
	

CDnaMeasurement::~CDnaMeasurement()
{
}

// move constructor
CDnaMeasurement::CDnaMeasurement(CDnaMeasurement&& m)
{
	m_strType = m.m_strType;
	m_bIgnore = m.m_bIgnore;
	m_strFirst = m.m_strFirst;
	m_MSmeasurementStations = m.m_MSmeasurementStations;
	m_lmeasurementIndex = m.m_lmeasurementIndex;
	m_lstn1Index = m.m_lstn1Index;
	m_lstn2Index = m.m_lstn2Index;
	m_lstn3Index = m.m_lstn3Index;

	m_measAdj = m.m_measAdj;
	m_measCorr = m.m_measCorr;
	m_measAdjPrec = m.m_measAdjPrec;
	m_residualPrec = m.m_residualPrec;
	m_preAdjCorr = m.m_preAdjCorr;

	m_epoch = m.m_epoch;
	m_epsgCode = m.m_epsgCode;

	m_msr_db_map = m.m_msr_db_map;

	m_bInsufficient = m.m_bInsufficient;
}

// move assignment operator
CDnaMeasurement& CDnaMeasurement::operator= (CDnaMeasurement&& rhs)
{
	// check for assignment to self!
	if (this == &rhs)
		return *this;

	m_strFirst = rhs.m_strFirst;
	m_MSmeasurementStations = rhs.m_MSmeasurementStations;
	m_strType = rhs.m_strType;
	m_bIgnore = rhs.m_bIgnore;
	m_lmeasurementIndex = rhs.m_lmeasurementIndex;
	m_lstn1Index = rhs.m_lstn1Index;
	m_lstn2Index = rhs.m_lstn2Index;
	m_lstn3Index = rhs.m_lstn3Index;

	m_measAdj = rhs.m_measAdj;
	m_measCorr = rhs.m_measCorr;
	m_measAdjPrec = rhs.m_measAdjPrec;
	m_residualPrec = rhs.m_residualPrec;
	m_preAdjCorr = rhs.m_preAdjCorr;

	m_epoch = rhs.m_epoch;
	m_epsgCode = rhs.m_epsgCode;

	m_msr_db_map = rhs.m_msr_db_map;

	m_bInsufficient = rhs.m_bInsufficient;

	return *this;
}
	

void CDnaMeasurement::coutMeasurement(std::ostream& os) const
{
	os << std::setw(2) << std::left << "+ " << std::setw(2) << m_strType;
}

void CDnaMeasurement::SetMeasurementDBID(const std::string& str)
{
	if (str.empty())
	{
		m_msr_db_map.msr_id = 0;
		m_msr_db_map.is_msr_id_set = false;
	}
	else
	{
		m_msr_db_map.msr_id = LongFromString<UINT32>(str);
		m_msr_db_map.is_msr_id_set = true;
	}
}

void CDnaMeasurement::SetClusterDBID(const std::string& str)
{
	if (str.empty())
	{
		m_msr_db_map.cluster_id = 0;
		m_msr_db_map.is_cls_id_set = false;
	}
	else
	{
		m_msr_db_map.cluster_id = LongFromString<UINT32>(str);
		m_msr_db_map.is_cls_id_set = true;
	}
}

void CDnaMeasurement::SetClusterDBID(const UINT32& u, bool s) 
{
	m_msr_db_map.cluster_id = u;
	m_msr_db_map.is_cls_id_set = s;
}

void CDnaMeasurement::SetType(const std::string& str)
{ 
	m_strType = trimstr(str);
	str_toupper<int>(m_strType);
}

void CDnaMeasurement::SetDatabaseMap(const msr_database_id_map& dbidmap) 
{
	m_msr_db_map = dbidmap;
}

void CDnaMeasurement::SetDatabaseMaps(it_vdbid_t& it_dbidmap)
{
	m_msr_db_map = *it_dbidmap;
}
	

void CDnaMeasurement::SerialiseDatabaseMap(std::ofstream* os)
{
	UINT16 val;

	os->write(reinterpret_cast<const char*>(&m_msr_db_map.msr_id), sizeof(UINT32));
	os->write(reinterpret_cast<const char*>(&m_msr_db_map.cluster_id), sizeof(UINT32));
	val = static_cast<UINT16>(m_msr_db_map.is_msr_id_set);
	os->write(reinterpret_cast<const char*>(&val), sizeof(UINT16));
	val = static_cast<UINT16>(m_msr_db_map.is_cls_id_set);
	os->write(reinterpret_cast<const char*>(&val), sizeof(UINT16));
}

void CDnaMeasurement::SetEpoch(const std::string& epoch)
{
	m_epoch = epoch;
}

}	// namespace measurements
}	// namespace dynadjust


