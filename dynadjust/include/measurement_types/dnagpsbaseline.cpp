//============================================================================
// Name         : dnagpsbaseline.cpp
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
// Description  : CDnaGpsBaseline and CDnaGpsBaselineCluster implementation file
//============================================================================

#include <include/exception/dnaexception.hpp>
#include <include/parameters/dnaepsg.hpp>
#include <include/measurement_types/dnagpsbaseline.hpp>

using namespace dynadjust::exception;
using namespace dynadjust::epsg;

namespace dynadjust {
namespace measurements {

CDnaGpsBaseline::CDnaGpsBaseline(void)
	: m_strTarget("")
	, m_lRecordedTotal(0)
	, m_dX(0.)
	, m_dY(0.)
	, m_dZ(0.)
	, m_dSigmaXX(0.)
	, m_dSigmaXY(0.)
	, m_dSigmaXZ(0.)
	, m_dSigmaYY(0.)
	, m_dSigmaYZ(0.)
	, m_dSigmaZZ(0.)
	, m_dPscale(1.)
	, m_dLscale(1.)
	, m_dHscale(1.)
	, m_dVscale(1.)
	, m_referenceFrame(DEFAULT_DATUM)
	, m_lclusterID(0)
{
	SetEpoch(DEFAULT_EPOCH);
	SetEpsg(epsgStringFromName<std::string>(m_referenceFrame));

	m_vGpsCovariances.clear();
	m_vGpsCovariances.reserve(0);
	m_MSmeasurementStations = TWO_STATION;
}


CDnaGpsBaseline::~CDnaGpsBaseline(void)
{

}

// move constructor
CDnaGpsBaseline::CDnaGpsBaseline(CDnaGpsBaseline&& g)
{
	m_strType = g.m_strType;
	m_strFirst = g.m_strFirst;
	m_strTarget = g.m_strTarget;
	m_bIgnore = g.m_bIgnore;
	m_lRecordedTotal = g.m_lRecordedTotal;
	
	m_referenceFrame = g.m_referenceFrame;
	m_epsgCode = g.m_epsgCode;
	m_epoch = g.m_epoch;
	
	m_dX = g.m_dX;
	m_dY = g.m_dY;
	m_dZ = g.m_dZ;
	m_dSigmaXX = g.m_dSigmaXX;
	m_dSigmaXY = g.m_dSigmaXY;
	m_dSigmaXZ = g.m_dSigmaXZ;
	m_dSigmaYY = g.m_dSigmaYY;
	m_dSigmaYZ = g.m_dSigmaYZ;
	m_dSigmaZZ = g.m_dSigmaZZ;
	
	m_dPscale = g.m_dPscale;
	m_dLscale = g.m_dLscale;
	m_dHscale = g.m_dHscale;
	m_dVscale = g.m_dVscale;
	
	m_lclusterID = g.m_lclusterID;
	m_MSmeasurementStations = g.m_MSmeasurementStations;
	
	m_vGpsCovariances = std::move(g.m_vGpsCovariances);
	
	m_msr_db_map = g.m_msr_db_map;
}

// move assignment operator
CDnaGpsBaseline& CDnaGpsBaseline::operator= (CDnaGpsBaseline&& rhs)
{
	// check for assignment to self!
	if (this == &rhs)
		return *this;

	CDnaMeasurement::operator=(std::move(rhs));
	m_strTarget = rhs.m_strTarget;
	m_lRecordedTotal = rhs.m_lRecordedTotal;

	m_referenceFrame = rhs.m_referenceFrame;
	m_epsgCode = rhs.m_epsgCode;
	m_epoch = rhs.m_epoch;

	m_dX = rhs.m_dX;
	m_dY = rhs.m_dY;
	m_dZ = rhs.m_dZ;
	m_dSigmaXX = rhs.m_dSigmaXX;
	m_dSigmaXY = rhs.m_dSigmaXY;
	m_dSigmaXZ = rhs.m_dSigmaXZ;
	m_dSigmaYY = rhs.m_dSigmaYY;
	m_dSigmaYZ = rhs.m_dSigmaYZ;
	m_dSigmaZZ = rhs.m_dSigmaZZ;

	m_dPscale = rhs.m_dPscale;
	m_dLscale = rhs.m_dLscale;
	m_dHscale = rhs.m_dHscale;
	m_dVscale = rhs.m_dVscale;

	m_lclusterID = rhs.m_lclusterID;
	m_MSmeasurementStations = rhs.m_MSmeasurementStations;

	m_vGpsCovariances = std::move(rhs.m_vGpsCovariances);

	m_msr_db_map = rhs.m_msr_db_map;

	return *this;
}
	

bool CDnaGpsBaseline::operator== (const CDnaGpsBaseline& rhs) const
{
	return (
		m_strFirst == rhs.m_strFirst &&
		m_strTarget == rhs.m_strTarget &&
		m_strType == rhs.m_strType &&
		m_lRecordedTotal == rhs.m_lRecordedTotal &&
		m_bIgnore == rhs.m_bIgnore &&
		fabs(m_dX - rhs.m_dX) < PRECISION_1E4 &&
		fabs(m_dY - rhs.m_dY) < PRECISION_1E4 &&
		fabs(m_dZ - rhs.m_dZ) < PRECISION_1E4 &&
		boost::iequals(m_referenceFrame, rhs.m_referenceFrame) &&
		m_epoch == rhs.m_epoch
		);
}

bool CDnaGpsBaseline::operator< (const CDnaGpsBaseline& rhs) const
{
	if (m_strFirst == rhs.m_strFirst) {
		if (m_strTarget == rhs.m_strTarget) {
			if (m_strType == rhs.m_strType) {
				if (m_lRecordedTotal == rhs.m_lRecordedTotal) {
					if (m_bIgnore == rhs.m_bIgnore) {
						if (m_epoch == rhs.m_epoch) {
							if (fabs(m_dX - rhs.m_dX) < PRECISION_1E4) {
								if (fabs(m_dY - rhs.m_dY) < PRECISION_1E4) {
									return m_dZ < rhs.m_dZ; }
								else
									return m_dY < rhs.m_dY; }
							else
								return m_dX < rhs.m_dX; }
						else
							return m_epoch < rhs.m_epoch; }
					else
						return m_bIgnore < rhs.m_bIgnore; }
				else
					return m_lRecordedTotal < rhs.m_lRecordedTotal;	}
			else
				return m_strType < rhs.m_strType; }
		else
			return m_strTarget < rhs.m_strTarget; }
	else
		return m_strFirst < rhs.m_strFirst;
}

void CDnaGpsBaseline::ReserveGpsCovariancesCount(const UINT32& size)
{
	m_vGpsCovariances.reserve(size);
}


void CDnaGpsBaseline::ResizeGpsCovariancesCount(const UINT32& size)
{
	m_vGpsCovariances.resize(size);
}


void CDnaGpsBaseline::AddGpsCovariance(const CDnaCovariance* pGpsCovariance)
{
	m_vGpsCovariances.push_back(std::move((CDnaCovariance&)*pGpsCovariance));
}
	

UINT32 CDnaGpsBaseline::CalcBinaryRecordCount() const
{
	UINT32 RecordCount = 3;
	std::vector<CDnaCovariance>::const_iterator _it_cov = m_vGpsCovariances.begin();
	for (; _it_cov!=m_vGpsCovariances.end(); ++_it_cov)
		RecordCount += _it_cov->CalcBinaryRecordCount();
	return RecordCount;
}
	

void CDnaGpsBaseline::WriteDynaMLMsr(std::ofstream* dynaml_stream, const std::string&, bool) const
{
	*dynaml_stream << "    <First>" << m_strFirst << "</First>" << std::endl;
	*dynaml_stream << "    <Second>" << m_strTarget << "</Second>" << std::endl;

	UINT32 precision = 4;
	
	*dynaml_stream << "    <GPSBaseline>" << std::endl;
	*dynaml_stream << "      <X>" << std::fixed << std::setprecision(precision) << m_dX << "</X>" << std::endl;
	*dynaml_stream << "      <Y>" << m_dY << "</Y>" << std::endl;
	*dynaml_stream << "      <Z>" << m_dZ << "</Z>" << std::endl;

	if (m_msr_db_map.is_msr_id_set)
		*dynaml_stream << "      <MeasurementID>" << m_msr_db_map.msr_id << "</MeasurementID>" << std::endl;

	*dynaml_stream << "      <SigmaXX>" << std::scientific << std::setprecision(13) << m_dSigmaXX << "</SigmaXX>" << std::endl;
	*dynaml_stream << "      <SigmaXY>" << m_dSigmaXY << "</SigmaXY>" << std::endl;
	*dynaml_stream << "      <SigmaXZ>" << m_dSigmaXZ << "</SigmaXZ>" << std::endl;
	*dynaml_stream << "      <SigmaYY>" << m_dSigmaYY << "</SigmaYY>" << std::endl;
	*dynaml_stream << "      <SigmaYZ>" << m_dSigmaYZ << "</SigmaYZ>" << std::endl;
	*dynaml_stream << "      <SigmaZZ>" << m_dSigmaZZ << "</SigmaZZ>" << std::endl;

	// write GPSPoint covariances
	std::vector<CDnaCovariance>::const_iterator _it_cov = m_vGpsCovariances.begin();
	for (; _it_cov!=m_vGpsCovariances.end(); ++_it_cov)
		_it_cov->WriteDynaMLMsr(dynaml_stream);

	*dynaml_stream << "    </GPSBaseline>" << std::endl;

}


void CDnaGpsBaseline::WriteDNAMsr(std::ofstream* dna_stream, const dna_msr_fields& dmw, const dna_msr_fields& dml, bool) const
{
	*dna_stream << std::setw(dmw.msr_type) << m_strType;
	if (m_bIgnore)
		*dna_stream << std::setw(dmw.msr_ignore) << "*";
	else
		*dna_stream << std::setw(dmw.msr_ignore) << " ";

	*dna_stream << std::left << std::setw(dmw.msr_inst) << m_strFirst;
	*dna_stream << std::left << std::setw(dmw.msr_targ1) << m_strTarget;

	// Print header for G baseline and first X cluster baseline
	bool printHeader(true);
	
	if (GetTypeC() == 'X')
	{
		if (m_lRecordedTotal != m_vGpsCovariances.size() + 1)
		{
			printHeader = false;

			// print database ids
			if (m_msr_db_map.is_msr_id_set)
				*dna_stream << std::setw(dml.msr_id_msr - dml.msr_targ2) << " " <<
					std::right << std::setw(dmw.msr_id_msr) << m_msr_db_map.msr_id;
			if (m_msr_db_map.is_cls_id_set)
				*dna_stream << std::setw(dmw.msr_id_cluster) << m_msr_db_map.cluster_id;			
		}
	}

	if (printHeader)
	{
		if (GetTypeC() == 'X')
			*dna_stream << std::left << std::setw(dmw.msr_targ2) << m_lRecordedTotal;
		else
			*dna_stream << std::right << std::setw(dmw.msr_targ2) << " ";
		
		// print scaling
		*dna_stream << 
			std::fixed << std::setprecision(2) << 
			std::right << std::setw(dmw.msr_gps_vscale) << double_string_width<double, UINT32, std::string>(m_dVscale, dmw.msr_gps_vscale) <<
			std::right << std::setw(dmw.msr_gps_pscale) << double_string_width<double, UINT32, std::string>(m_dPscale, dmw.msr_gps_vscale) <<
			std::right << std::setw(dmw.msr_gps_lscale) << double_string_width<double, UINT32, std::string>(m_dLscale, dmw.msr_gps_vscale) <<
			std::right << std::setw(dmw.msr_gps_hscale) << double_string_width<double, UINT32, std::string>(m_dHscale, dmw.msr_gps_vscale);

		// print reference frame and epoch
		*dna_stream <<
			std::right << std::setw(dmw.msr_gps_reframe) << m_referenceFrame <<
			std::right << std::setw(dmw.msr_gps_epoch) << m_epoch;

		// print database ids
		if (m_msr_db_map.is_msr_id_set)
			*dna_stream << std::setw(dmw.msr_id_msr) << m_msr_db_map.msr_id;
		if (m_msr_db_map.is_cls_id_set)
			*dna_stream << std::setw(dmw.msr_id_cluster) << m_msr_db_map.cluster_id;
	}

	*dna_stream << std::endl;

	UINT32 pad(dmw.msr_type + dmw.msr_ignore + dmw.msr_inst + dmw.msr_targ1 + dmw.msr_targ2);

	// X
	*dna_stream << std::setw(pad) << " ";
	*dna_stream << std::right << std::setw(dmw.msr_gps) << std::fixed << std::setprecision(4) << m_dX;
	*dna_stream << std::right << std::setw(dmw.msr_gps_vcv_1) << std::scientific << std::setprecision(13) << m_dSigmaXX;
	*dna_stream << std::endl;
		
	// Y
	*dna_stream << std::setw(pad) << " ";
	*dna_stream << std::right << std::setw(dmw.msr_gps) << std::fixed << std::setprecision(4) << m_dY;
	*dna_stream << 
		std::right << std::setw(dmw.msr_gps_vcv_1) << std::scientific << std::setprecision(13) << m_dSigmaXY << 
		std::right << std::setw(dmw.msr_gps_vcv_2) << m_dSigmaYY;
	*dna_stream << std::endl;

	// Z
	*dna_stream << std::setw(pad) << " ";
	*dna_stream << std::right << std::setw(dmw.msr_gps) << std::fixed << std::setprecision(4) << m_dZ;
	*dna_stream << 
		std::right << std::setw(dmw.msr_gps_vcv_1) << std::scientific << std::setprecision(13) << m_dSigmaXZ <<
		std::right << std::setw(dmw.msr_gps_vcv_2) << m_dSigmaYZ << 
		std::right << std::setw(dmw.msr_gps_vcv_3) << m_dSigmaZZ;
	*dna_stream << std::endl;

	// write GPSBaseline covariances (not supported by DNA format)
	std::vector<CDnaCovariance>::const_iterator _it_cov = m_vGpsCovariances.begin();
	for (_it_cov=m_vGpsCovariances.begin(); _it_cov!=m_vGpsCovariances.end(); ++_it_cov)
		_it_cov->WriteDNAMsr(dna_stream, dmw, dml);
}
	

void CDnaGpsBaseline::SimulateMsr(vdnaStnPtr* vStations, const CDnaEllipsoid* ellipsoid)
{
	m_dX = (vStations->at(m_lstn2Index).get()->GetXAxis() - vStations->at(m_lstn1Index).get()->GetXAxis());
	m_dY = (vStations->at(m_lstn2Index).get()->GetYAxis() - vStations->at(m_lstn1Index).get()->GetYAxis());
	m_dZ = (vStations->at(m_lstn2Index).get()->GetZAxis() - vStations->at(m_lstn1Index).get()->GetZAxis());

	m_dSigmaXX = 4.022E-05;                    
	m_dSigmaXY = -1.369E-05;
	m_dSigmaXZ = 3.975E-05;
	m_dSigmaYY = 1.487E-05;
	m_dSigmaYZ = -2.035E-05;
	m_dSigmaZZ = 6.803E-05;

	std::vector<CDnaCovariance>::iterator _it_cov = m_vGpsCovariances.begin();
	for (_it_cov=m_vGpsCovariances.begin(); _it_cov!=m_vGpsCovariances.end(); ++_it_cov)
		_it_cov->SimulateMsr(vStations, ellipsoid);
}
	

UINT32 CDnaGpsBaseline::SetMeasurementRec(const vstn_t& binaryStn, it_vmsr_t& it_msr, it_vdbid_t& dbidmap)
{
	// first station
	m_lstn1Index = it_msr->station1;
	m_strFirst = binaryStn.at(it_msr->station1).stationName;

	// target station
	m_lstn2Index = it_msr->station2;
	m_strTarget = binaryStn.at(it_msr->station2).stationName;
	
	m_dPscale = it_msr->scale1;
	m_dLscale = it_msr->scale2;
	m_dHscale = it_msr->scale3;
	m_dVscale = it_msr->scale4;

	m_epoch = it_msr->epoch;
	m_epsgCode = it_msr->epsgCode;
	m_referenceFrame = datumFromEpsgString<std::string>(it_msr->epsgCode);

	m_lclusterID = it_msr->clusterID;
	m_MSmeasurementStations = (MEASUREMENT_STATIONS)it_msr->measurementStations;

	// X, sigmaXX
	m_bIgnore = it_msr->ignore;
	m_strType = it_msr->measType;
	m_dX = it_msr->term1;
	m_dSigmaXX = it_msr->term2;
	m_lRecordedTotal = it_msr->vectorCount1;

	it_msr++;
	
	m_dY = it_msr->term1;
	m_dSigmaXY = it_msr->term2;
	m_dSigmaYY = it_msr->term3;

	it_msr++;

	m_dZ = it_msr->term1;
	m_dSigmaXZ = it_msr->term2;
	m_dSigmaYZ = it_msr->term3;
	m_dSigmaZZ = it_msr->term4;
	m_lRecordedTotal = it_msr->vectorCount1;

	m_vGpsCovariances.clear();
	m_vGpsCovariances.resize(it_msr->vectorCount2);

	CDnaMeasurement::SetDatabaseMap(*dbidmap);
	dbidmap += 2;
	dbidmap += (it_msr->vectorCount2 * 3);

	// now covariances
	std::vector<CDnaCovariance>::iterator _it_cov = m_vGpsCovariances.begin();
	for (; _it_cov!=m_vGpsCovariances.end(); ++_it_cov)
		_it_cov->SetMeasurementRec(binaryStn, it_msr);

	return it_msr->vectorCount1;
}


void CDnaGpsBaseline::WriteBinaryMsr(std::ofstream* binary_stream, PUINT32 msrIndex) const
{
	measurement_t measRecord;
	measRecord.measType = GetTypeC();
	measRecord.measStart = xMeas;
	measRecord.ignore = m_bIgnore;
	measRecord.station1 = m_lstn1Index;
	measRecord.station2 = m_lstn2Index;
	measRecord.vectorCount1 = m_lRecordedTotal;									// number of GpsBaselines in the parent cluster
	measRecord.vectorCount2 = static_cast<UINT32>(m_vGpsCovariances.size());	// number of Covariances in the GpsBaseline

	measRecord.clusterID = m_lclusterID;
	measRecord.measurementStations = m_MSmeasurementStations;

	measRecord.scale1 = m_dPscale;
	measRecord.scale2 = m_dLscale;
	measRecord.scale3 = m_dHscale;
	measRecord.scale4 = m_dVscale;

	snprintf(measRecord.epsgCode, sizeof(measRecord.epsgCode), "%s", m_epsgCode.substr(0, STN_EPSG_WIDTH).c_str());
	snprintf(measRecord.epoch, sizeof(measRecord.epoch), "%s", m_epoch.substr(0, STN_EPOCH_WIDTH).c_str());

	// X
	measRecord.measAdj = m_measAdj;
	measRecord.measCorr = m_measCorr;
	measRecord.measAdjPrec = m_measAdjPrec;
	measRecord.residualPrec = m_residualPrec;
	measRecord.preAdjCorr = m_preAdjCorr;
	measRecord.term1 = m_dX;
	measRecord.term2 = m_dSigmaXX;	// already a variance
	measRecord.term3 = 0.;
	measRecord.term4 = 0.;
	measRecord.fileOrder = ((*msrIndex)++);

	binary_stream->write(reinterpret_cast<char *>(&measRecord), sizeof(measurement_t));

	// Y
	measRecord.measStart = yMeas;
	measRecord.measAdj = m_measAdj;
	measRecord.measCorr = m_measCorr;
	measRecord.measAdjPrec = m_measAdjPrec;
	measRecord.residualPrec = m_residualPrec;
	measRecord.preAdjCorr = m_preAdjCorr;
	measRecord.term1 = m_dY;
	measRecord.term2 = m_dSigmaXY;	// already a variance
	measRecord.term3 = m_dSigmaYY;
	measRecord.term4 = 0.;
	measRecord.fileOrder = ((*msrIndex)++);

	binary_stream->write(reinterpret_cast<char *>(&measRecord), sizeof(measurement_t));

	// Z
	measRecord.measStart = zMeas;
	measRecord.measAdj = m_measAdj;
	measRecord.measCorr = m_measCorr;
	measRecord.measAdjPrec = m_measAdjPrec;
	measRecord.residualPrec = m_residualPrec;
	measRecord.preAdjCorr = m_preAdjCorr;
	measRecord.term1 = m_dZ;
	measRecord.term2 = m_dSigmaXZ;		// already a variance
	measRecord.term3 = m_dSigmaYZ;
	measRecord.term4 = m_dSigmaZZ;
	measRecord.fileOrder = ((*msrIndex)++);

	binary_stream->write(reinterpret_cast<char *>(&measRecord), sizeof(measurement_t));

	// now write covariance elements
	std::vector<CDnaCovariance>::const_iterator _it_cov;
	for (_it_cov=m_vGpsCovariances.begin(); _it_cov!=m_vGpsCovariances.end(); ++_it_cov)
		_it_cov->WriteBinaryMsr(binary_stream, msrIndex, m_epsgCode, m_epoch);
}

void CDnaGpsBaseline::SerialiseDatabaseMap(std::ofstream* os)
{
	// X
	CDnaMeasurement::SerialiseDatabaseMap(os);
	
	// Y
	CDnaMeasurement::SerialiseDatabaseMap(os);
	
	// Z
	CDnaMeasurement::SerialiseDatabaseMap(os);

	for_each(m_vGpsCovariances.begin(), m_vGpsCovariances.end(),
		[this, os](const CDnaCovariance& cov) {
		((CDnaCovariance*)&cov)->SerialiseDatabaseMap(os, m_msr_db_map);
	});
}



void CDnaGpsBaseline::SetX(const std::string& str)
{
	DoubleFromString(m_dX, trimstr(str));
}

void CDnaGpsBaseline::SetY(const std::string& str)
{
	DoubleFromString(m_dY, trimstr(str));
}

void CDnaGpsBaseline::SetZ(const std::string& str)
{
	DoubleFromString(m_dZ, trimstr(str));
}

void CDnaGpsBaseline::SetSigmaXX(const std::string& str)
{
	if (DoubleFromString_ZeroCheck(m_dSigmaXX, trimstr(str)))
	{
		std::stringstream ss;
		ss << "SetSigmaXX(): Variances cannot be zero: " << str << ".";
		throw XMLInteropException(ss.str(), 0);
	}
}

void CDnaGpsBaseline::SetSigmaXY(const std::string& str)
{
	DoubleFromString(m_dSigmaXY, trimstr(str));
}

void CDnaGpsBaseline::SetSigmaXZ(const std::string& str)
{
	DoubleFromString(m_dSigmaXZ, trimstr(str));
}

void CDnaGpsBaseline::SetSigmaYY(const std::string& str)
{
	if (DoubleFromString_ZeroCheck(m_dSigmaYY, trimstr(str)))
	{
		std::stringstream ss;
		ss << "SetSigmaYY(): Variances cannot be zero: " << str << ".";
		throw XMLInteropException(ss.str(), 0);
	}
}

void CDnaGpsBaseline::SetSigmaYZ(const std::string& str)
{
	DoubleFromString(m_dSigmaYZ, trimstr(str));
}

void CDnaGpsBaseline::SetSigmaZZ(const std::string& str)
{
	if (DoubleFromString_ZeroCheck(m_dSigmaZZ, trimstr(str)))
	{
		std::stringstream ss;
		ss << "SetSigmaZZ(): Variances cannot be zero: " << str << ".";
		throw XMLInteropException(ss.str(), 0);
	}
}
	

void CDnaGpsBaseline::SetEpsg(const std::string& epsg) 
{
	m_epsgCode = epsg;
}


void CDnaGpsBaseline::SetReferenceFrame(const std::string& refFrame) 
{
	m_referenceFrame = refFrame;
	SetEpsg(epsgStringFromName<std::string>(refFrame));
}


void CDnaGpsBaseline::SetPscale(const std::string& str)
{
	DoubleFromString(m_dPscale, trimstr(str));
}

void CDnaGpsBaseline::SetLscale(const std::string& str)
{
	DoubleFromString(m_dLscale, trimstr(str));
}

void CDnaGpsBaseline::SetHscale(const std::string& str)
{
	DoubleFromString(m_dHscale, trimstr(str));
}

void CDnaGpsBaseline::SetVscale(const std::string& str)
{
	DoubleFromString(m_dVscale, trimstr(str));
}









CDnaGpsBaselineCluster::CDnaGpsBaselineCluster(void)
	: m_strTarget("")
	, m_lRecordedTotal(1)
	, m_dPscale(1.)
	, m_dLscale(1.)
	, m_dHscale(1.)
	, m_dVscale(1.)
	, m_referenceFrame(DEFAULT_DATUM)
	, m_lclusterID(0)
{
	SetEpoch(DEFAULT_EPOCH);
	SetEpsg(epsgStringFromName<std::string>(m_referenceFrame));

	m_vGpsBaselines.clear();
	m_vGpsBaselines.reserve(1);

	// m_strType will be set by the parsing function as "X" or "G"
	m_MSmeasurementStations = TWO_STATION;
}


CDnaGpsBaselineCluster::~CDnaGpsBaselineCluster(void)
{

}

// move constructor
CDnaGpsBaselineCluster::CDnaGpsBaselineCluster(CDnaGpsBaselineCluster&& g)
{
	m_strType = g.m_strType;
	m_bIgnore = g.m_bIgnore;
	m_strTarget = g.m_strTarget;
	m_lRecordedTotal = g.m_lRecordedTotal;
	m_vGpsBaselines = std::move(g.m_vGpsBaselines);

	m_dPscale = g.m_dPscale;
	m_dLscale = g.m_dLscale;
	m_dHscale = g.m_dHscale;
	m_dVscale = g.m_dVscale;
	m_lclusterID = g.m_lclusterID;
	m_MSmeasurementStations = g.m_MSmeasurementStations;
	
	m_referenceFrame = g.m_referenceFrame;
	m_epsgCode = g.m_epsgCode;
	m_epoch = g.m_epoch;

	m_msr_db_map = g.m_msr_db_map;

	m_dbidmap = g.m_dbidmap;
}

// move assignment operator 
CDnaGpsBaselineCluster& CDnaGpsBaselineCluster::operator= (CDnaGpsBaselineCluster&& rhs)
{
	// check for assignment to self!
	if (this == &rhs)
		return *this;

	CDnaMeasurement::operator=(std::move(rhs));
	m_strTarget = rhs.m_strTarget;
	m_lRecordedTotal = rhs.m_lRecordedTotal;

	m_referenceFrame = rhs.m_referenceFrame;
	m_epsgCode = rhs.m_epsgCode;
	m_epoch = rhs.m_epoch;

	m_dPscale = rhs.m_dPscale;
	m_dLscale = rhs.m_dLscale;
	m_dHscale = rhs.m_dHscale;
	m_dVscale = rhs.m_dVscale;

	m_lclusterID = rhs.m_lclusterID;
	m_MSmeasurementStations = rhs.m_MSmeasurementStations;

	m_vGpsBaselines = std::move(rhs.m_vGpsBaselines);

	m_msr_db_map = rhs.m_msr_db_map;

	m_dbidmap = rhs.m_dbidmap;

	return *this;
}
	

CDnaGpsBaselineCluster::CDnaGpsBaselineCluster(const UINT32 lclusterID, const std::string& referenceframe, const std::string& epoch)
	: m_strTarget("")
	, m_lRecordedTotal(1)
	, m_dPscale(1.)
	, m_dLscale(1.)
	, m_dHscale(1.)
	, m_dVscale(1.)
	, m_referenceFrame(referenceframe)
	, m_lclusterID(lclusterID)
{
	SetEpoch(epoch);
	SetEpsg(epsgStringFromName<std::string>(referenceframe));

	m_vGpsBaselines.clear();
	m_vGpsBaselines.reserve(1);
	// m_strType will be set by the parsing function as "X" or "G"
	m_MSmeasurementStations = TWO_STATION;
}
	

bool CDnaGpsBaselineCluster::operator== (const CDnaGpsBaselineCluster& rhs) const
{
	return (
		m_strType == rhs.m_strType &&
		m_strFirst == rhs.m_strFirst &&
		m_bIgnore == rhs.m_bIgnore &&
		m_strTarget == rhs.m_strTarget &&
		m_lRecordedTotal == rhs.m_lRecordedTotal &&
		m_dPscale == rhs.m_dPscale &&
		m_dLscale == rhs.m_dLscale &&
		m_dHscale == rhs.m_dHscale &&
		m_dVscale == rhs.m_dVscale &&
		m_vGpsBaselines == rhs.m_vGpsBaselines &&
		boost::iequals(m_referenceFrame, rhs.m_referenceFrame) &&
		m_epoch == rhs.m_epoch
		);
}

bool CDnaGpsBaselineCluster::operator< (const CDnaGpsBaselineCluster& rhs) const
{
	if (m_strFirst == rhs.m_strFirst) {
		if (m_strType == rhs.m_strType) {	// could be G, X
			if (m_bIgnore == rhs.m_bIgnore) {
				if (m_epoch == rhs.m_epoch) {
					if (m_strTarget == rhs.m_strTarget) {
						if (m_lRecordedTotal == rhs.m_lRecordedTotal) 
							return m_vGpsBaselines < rhs.m_vGpsBaselines;
						else
							return m_lRecordedTotal < rhs.m_lRecordedTotal; }
					else
						return m_strTarget < rhs.m_strTarget; }
				else
					return m_epoch < rhs.m_epoch; }
			else
				return m_bIgnore < rhs.m_bIgnore; }
		else
			return m_strType < rhs.m_strType; }
	else
		return m_strFirst < rhs.m_strFirst;
}


void CDnaGpsBaselineCluster::ReserveGpsBaselinesCount(const UINT32& size)
{
	m_vGpsBaselines.reserve(size);
}

void CDnaGpsBaselineCluster::AddGpsBaseline(const CDnaMeasurement* pGpsBaseline)
{
	m_vGpsBaselines.push_back(std::move((CDnaGpsBaseline&)*pGpsBaseline));
}


void CDnaGpsBaselineCluster::ClearBaselines()
{
	m_vGpsBaselines.clear();
}
	

void CDnaGpsBaselineCluster::SerialiseDatabaseMap(std::ofstream* os)
{
	for_each(m_vGpsBaselines.begin(), m_vGpsBaselines.end(),
		[this, os](const CDnaGpsBaseline& bsl) {
		((CDnaGpsBaseline*)&bsl)->SerialiseDatabaseMap(os);
	});
	
}
	

void CDnaGpsBaselineCluster::SetDatabaseMaps(it_vdbid_t& dbidmap) 
{
	m_dbidmap = dbidmap;

	CDnaMeasurement::SetDatabaseMap(*m_dbidmap);

	for_each(m_vGpsBaselines.begin(), m_vGpsBaselines.end(),
		[this](const CDnaGpsBaseline& bsl) {
			((CDnaGpsBaseline*)&bsl)->SetDatabaseMap(*m_dbidmap);
			m_dbidmap += ((CDnaGpsBaseline*)&bsl)->GetCovariances_ptr()->size();
		});
	
}
	

UINT32 CDnaGpsBaselineCluster::CalcBinaryRecordCount() const
{
	UINT32 recordCount(0);
	for_each(m_vGpsBaselines.begin(), m_vGpsBaselines.end(),
		[&recordCount](const CDnaGpsBaseline& bsl) {
			recordCount += bsl.CalcBinaryRecordCount();
	});
	return recordCount;
}


void CDnaGpsBaselineCluster::WriteDynaMLMsr(std::ofstream* dynaml_stream, const std::string& comment, bool) const
{
	const size_t bslCount = m_vGpsBaselines.size();
	char cType = GetTypeC();

	if (comment.empty())
	{
		*dynaml_stream << "  <!-- Type " << measurement_name<char, std::string>(GetTypeC());
		if (cType == 'X')
		{
			if (bslCount > 1)
				*dynaml_stream << " (set of " << bslCount << ")" << std::endl;
			else
				*dynaml_stream << "  (single)" << std::endl;
		}
		*dynaml_stream << " -->" << std::endl;
	}
	else
		*dynaml_stream << "  <!-- " << comment << " -->" << std::endl;
	
	*dynaml_stream << "  <DnaMeasurement>" << std::endl;
	*dynaml_stream << "    <Type>" << cType << "</Type>" << std::endl;
	// Source file from which the measurement came
	*dynaml_stream << "    <Source>" << m_sourceFile << "</Source>" << std::endl;
	if (m_bIgnore)
		*dynaml_stream << "    <Ignore>*</Ignore>" << std::endl;
	else
		*dynaml_stream << "    <Ignore/>" << std::endl;
	
	// Reference frame and epoch
	*dynaml_stream << "    <ReferenceFrame>" << m_referenceFrame << "</ReferenceFrame>" << std::endl;
	*dynaml_stream << "    <Epoch>" << m_epoch << "</Epoch>" << std::endl;

	// Scalars
	*dynaml_stream << "    <Vscale>" << std::fixed << std::setprecision(3) << m_dVscale << "</Vscale>" << std::endl;
	*dynaml_stream << "    <Pscale>" << m_dPscale << "</Pscale>" << std::endl;
	*dynaml_stream << "    <Lscale>" << m_dLscale << "</Lscale>" << std::endl;
	*dynaml_stream << "    <Hscale>" << m_dHscale << "</Hscale>" << std::endl;

	if (m_msr_db_map.is_cls_id_set)
		*dynaml_stream << "    <ClusterID>" << m_msr_db_map.cluster_id << "</ClusterID>" << std::endl;

	// baseline count
	if (cType == 'X')
		*dynaml_stream << "    <Total>" << bslCount << "</Total>" << std::endl;

	// write GpsBaselines
	std::vector<CDnaGpsBaseline>::const_iterator _it_bsl;
	for (_it_bsl=m_vGpsBaselines.begin(); _it_bsl!=m_vGpsBaselines.end(); ++_it_bsl)
		_it_bsl->WriteDynaMLMsr(dynaml_stream, comment, true);
	
	*dynaml_stream << "  </DnaMeasurement>" << std::endl;

}
	

void CDnaGpsBaselineCluster::WriteDNAMsr(std::ofstream* dna_stream, const dna_msr_fields& dmw, const dna_msr_fields& dml, bool) const
{
	// write GpsBaselines
	std::vector<CDnaGpsBaseline>::const_iterator _it_bsl;
	for (_it_bsl=m_vGpsBaselines.begin(); _it_bsl!=m_vGpsBaselines.end(); ++_it_bsl)
		_it_bsl->WriteDNAMsr(dna_stream, dmw, dml, true);
}
	

void CDnaGpsBaselineCluster::SimulateMsr(vdnaStnPtr* vStations, const CDnaEllipsoid* ellipsoid)
{
	std::vector<CDnaGpsBaseline>::iterator _it_bsl;
	for (_it_bsl=m_vGpsBaselines.begin(); _it_bsl!=m_vGpsBaselines.end(); ++_it_bsl)
		_it_bsl->SimulateMsr(vStations, ellipsoid);
}
	

UINT32 CDnaGpsBaselineCluster::SetMeasurementRec(const vstn_t& binaryStn, it_vmsr_t& it_msr, it_vdbid_t& dbidmap)
{
	m_lclusterID = it_msr->clusterID;
	m_MSmeasurementStations = (MEASUREMENT_STATIONS)it_msr->measurementStations;
	m_lRecordedTotal = it_msr->vectorCount1;
	m_vGpsBaselines.clear();
	m_vGpsBaselines.resize(m_lRecordedTotal);

	m_referenceFrame = datumFromEpsgString<std::string>(it_msr->epsgCode);
	m_epoch = it_msr->epoch;
	m_epsgCode = it_msr->epsgCode;

	m_dPscale = it_msr->scale1;
	m_dLscale = it_msr->scale2;
	m_dHscale = it_msr->scale3;
	m_dVscale = it_msr->scale4;

	UINT32 measrecordCount(it_msr->vectorCount1);

	CDnaMeasurement::SetDatabaseMap(*dbidmap);

	// read remaining GpsPoint data and all Covariances from file
	for (UINT32 i=0; i<m_lRecordedTotal; i++)
	{
		if (i > 0)
		{
			it_msr++;
			dbidmap++;
		}

		m_strType = it_msr->measType;
		m_bIgnore = it_msr->ignore;
		m_vGpsBaselines.at(i).SetMeasurementRec(binaryStn, it_msr, dbidmap);
	}
	return measrecordCount - 1;	
}


void CDnaGpsBaselineCluster::WriteBinaryMsr(std::ofstream* binary_stream, PUINT32 msrIndex) const
{
	std::vector<CDnaGpsBaseline>::const_iterator _it_bsl;
	for (_it_bsl=m_vGpsBaselines.begin(); _it_bsl!=m_vGpsBaselines.end(); ++_it_bsl)
		_it_bsl->WriteBinaryMsr(binary_stream, msrIndex);
}


void CDnaGpsBaselineCluster::SetEpsg(const std::string& epsg) 
{
	// m_epsgCode is a member of CDnaMeasurement
	m_epsgCode = epsg;
}


void CDnaGpsBaselineCluster::SetReferenceFrame(const std::string& refFrame) 
{
	// m_referenceFrame is a member of CDnaGpsBaselineCluster
	m_referenceFrame = refFrame;
	SetEpsg(epsgStringFromName<std::string>(refFrame));
}


void CDnaGpsBaselineCluster::SetTotal(const std::string& str)
{
	m_lRecordedTotal = LongFromString<UINT32>(trimstr(str));
	m_vGpsBaselines.reserve(m_lRecordedTotal);
}

void CDnaGpsBaselineCluster::SetPscale(const std::string& str)
{
	DoubleFromString(m_dPscale, trimstr(str));
}

void CDnaGpsBaselineCluster::SetLscale(const std::string& str)
{
	DoubleFromString(m_dLscale, trimstr(str));
}

void CDnaGpsBaselineCluster::SetHscale(const std::string& str)
{
	DoubleFromString(m_dHscale, trimstr(str));
}

void CDnaGpsBaselineCluster::SetVscale(const std::string& str)
{
	DoubleFromString(m_dVscale, trimstr(str));
}




}	// namespace measurements
}	// namespace dynadjust
