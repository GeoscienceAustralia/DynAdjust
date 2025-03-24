//============================================================================
// Name         : dnaangle.cpp
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
// Description  : CDnaAngle implementation file
//============================================================================

#include <include/measurement_types/dnaangle.hpp>

namespace dynadjust {
namespace measurements {

CDnaAngle::CDnaAngle(void)
	: m_strTarget("")
	, m_strTarget2("")
	, m_drValue(0.)
	, m_dStdDev(0.)
{
	// from abstract class
	m_strType = "A";
	m_MSmeasurementStations = THREE_STATION;
}

CDnaAngle::~CDnaAngle(void)
{
}
	

bool CDnaAngle::operator== (const CDnaAngle& rhs) const
{
	return (
		m_strFirst == rhs.m_strFirst &&
		m_strType == rhs.m_strType &&
		m_bIgnore == rhs.m_bIgnore &&
		m_strTarget == rhs.m_strTarget &&
		m_strTarget2 == rhs.m_strTarget2 &&
		m_drValue == rhs.m_drValue &&
		m_dStdDev == rhs.m_dStdDev &&
		m_epoch == rhs.m_epoch
		);
}

bool CDnaAngle::operator< (const CDnaAngle& rhs) const
{
	if (m_strFirst == rhs.m_strFirst) {		
		if (m_strTarget == rhs.m_strTarget) {
			if (m_strTarget2 == rhs.m_strTarget2) {
				if (m_bIgnore == rhs.m_bIgnore) {
					if (m_epoch == rhs.m_epoch) {
						if (m_drValue == rhs.m_drValue)
							return m_dStdDev < rhs.m_dStdDev;
						else
							return m_drValue < rhs.m_drValue; }
					else
						return m_epoch < rhs.m_epoch; }
				else
					return m_bIgnore < rhs.m_bIgnore; }		
			else
				return m_strTarget2 < rhs.m_strTarget2; }
		else
			return m_strTarget < rhs.m_strTarget; }		
	else
		return m_strFirst < rhs.m_strFirst;
}
	

void CDnaAngle::WriteDynaMLMsr(std::ofstream* dynaml_stream, const std::string& comment, bool) const
{
	if (comment.empty())
		*dynaml_stream << "  <!-- Type " << measurement_name<char, std::string>(GetTypeC()) << " -->" << std::endl;
	else
		*dynaml_stream << "  <!-- " << comment << " -->" << std::endl;

	*dynaml_stream << "  <DnaMeasurement>" << std::endl;
	*dynaml_stream << "    <Type>" << m_strType << "</Type>" << std::endl;
	// Source file from which the measurement came
	*dynaml_stream << "    <Source>" << m_sourceFile << "</Source>" << std::endl;
	if (m_bIgnore)
		*dynaml_stream << "    <Ignore>*</Ignore>" << std::endl;
	else
		*dynaml_stream << "    <Ignore/>" << std::endl;

	if (m_epoch.empty())
		*dynaml_stream << "    <Epoch/>" << std::endl;
	else
		*dynaml_stream << "    <Epoch>" << m_epoch << "</Epoch>" << std::endl;

	*dynaml_stream << "    <First>" << m_strFirst << "</First>" << std::endl;
	*dynaml_stream << "    <Second>" << m_strTarget << "</Second>" << std::endl;
	*dynaml_stream << "    <Third>" << m_strTarget2 << "</Third>" << std::endl;
	*dynaml_stream << "    <Value>" << std::setprecision(8) << std::fixed << RadtoDms(m_drValue) << "</Value>" << std::endl;
	*dynaml_stream << "    <StdDev>" << std::scientific << std::setprecision(6) << Seconds(m_dStdDev) << "</StdDev>" << std::endl;
	
	if (m_msr_db_map.is_msr_id_set)
		*dynaml_stream << "    <MeasurementID>" << m_msr_db_map.msr_id << "</MeasurementID>" << std::endl;
	
	*dynaml_stream << "  </DnaMeasurement>" << std::endl;
}
	

void CDnaAngle::WriteDNAMsr(std::ofstream* dna_stream, const dna_msr_fields& dmw, const dna_msr_fields& dml, bool) const
{
	*dna_stream << std::setw(dmw.msr_type) << m_strType;
	if (m_bIgnore)
		*dna_stream << std::setw(dmw.msr_ignore) << "*";
	else
		*dna_stream << std::setw(dmw.msr_ignore) << " ";
	*dna_stream << std::left << std::setw(dmw.msr_inst) << m_strFirst;
	*dna_stream << std::left << std::setw(dmw.msr_targ1) << m_strTarget;
	*dna_stream << std::left << std::setw(dmw.msr_targ2) << m_strTarget2;
	*dna_stream << std::setw(dmw.msr_linear) << " ";	// linear measurement value
	*dna_stream << std::setw(dmw.msr_ang_d + dmw.msr_ang_m + dmw.msr_ang_s) << 
		std::right << FormatDnaDmsString(RadtoDms(m_drValue), 8);
	
	UINT32 m_stdDevPrec(3);
	*dna_stream << std::setw(dmw.msr_stddev) << StringFromTW(Seconds(m_dStdDev), dmw.msr_stddev, m_stdDevPrec);
	//*dna_stream << std::setw(dmw.msr_stddev) << std::fixed << std::setprecision(3) << Seconds(m_dStdDev);
	
	*dna_stream << std::setw(dml.msr_gps_epoch - dml.msr_inst_ht) << " ";
	*dna_stream << std::setw(dmw.msr_gps_epoch) << m_epoch;

	if (m_msr_db_map.is_msr_id_set)
		*dna_stream << std::setw(dmw.msr_id_msr) << m_msr_db_map.msr_id;
	
	*dna_stream << std::endl;
}
	

void CDnaAngle::SimulateMsr(vdnaStnPtr* vStations, const CDnaEllipsoid*)
{
	_it_vdnastnptr stn1_it(vStations->begin() + m_lstn1Index);
	_it_vdnastnptr stn2_it(vStations->begin() + m_lstn2Index);
	_it_vdnastnptr stn3_it(vStations->begin() + m_lstn3Index);

	double direction12, direction13;

	// calculated angle
	m_drValue = HorizontalAngle(
		stn1_it->get()->GetXAxis(),
		stn1_it->get()->GetYAxis(),
		stn1_it->get()->GetZAxis(),
		stn2_it->get()->GetXAxis(),
		stn2_it->get()->GetYAxis(),
		stn2_it->get()->GetZAxis(),
		stn3_it->get()->GetXAxis(),
		stn3_it->get()->GetYAxis(),
		stn3_it->get()->GetZAxis(),
		stn1_it->get()->GetcurrentLatitude(),
		stn1_it->get()->GetcurrentLongitude(),
		&direction12, &direction13);

	// Deflections available?
	if (fabs(stn1_it->get()->GetverticalDef()) > E4_SEC_DEFLECTION || fabs(stn1_it->get()->GetmeridianDef()) > E4_SEC_DEFLECTION)
	{
		/////////////////////////////////////////////////////////////////////////////////
		// Angles (observed or derived from directions) must be corrected for deflection 
		// of the vertical via "Laplace correction".  This correction requires zenith 
		// distance (zenith12, zenith13) and geodetic azimuth (direction12, direction13), 
		// both of which must be computed from coordinates.
		
		// 1. Compute zenith distance 1 -> 2
		double zenith12(ZenithDistance<double>(
			stn1_it->get()->GetXAxis(),
			stn1_it->get()->GetYAxis(),
			stn1_it->get()->GetZAxis(),
			stn2_it->get()->GetXAxis(),
			stn2_it->get()->GetYAxis(),
			stn2_it->get()->GetZAxis(),
			stn1_it->get()->GetcurrentLatitude(),
			stn1_it->get()->GetcurrentLongitude(),
			stn2_it->get()->GetcurrentLatitude(),
			stn2_it->get()->GetcurrentLongitude(),
			0.0,									// instrument height
			0.0));									// target height

		// 2. Compute zenith distance 1 -> 3
		double zenith13(ZenithDistance<double>(
			stn1_it->get()->GetXAxis(),
			stn1_it->get()->GetYAxis(),
			stn1_it->get()->GetZAxis(),
			stn3_it->get()->GetXAxis(),
			stn3_it->get()->GetYAxis(),
			stn3_it->get()->GetZAxis(),
			stn1_it->get()->GetcurrentLatitude(),
			stn1_it->get()->GetcurrentLongitude(),
			stn3_it->get()->GetcurrentLatitude(),
			stn3_it->get()->GetcurrentLongitude(),
			0.0,									// instrument height
			0.0));									// target height

		// 3. Laplace correction 1 -> 2 -> 3
		m_preAdjCorr = HzAngleDeflectionCorrection<double>(		
			direction12,							// geodetic azimuth 1 -> 2
			zenith12,								// zenith distance 1 -> 2
			direction13,							// geodetic azimuth 1 -> 3
			zenith13,								// zenith distance 1 -> 3
			stn1_it->get()->GetverticalDef(),		// deflection in prime vertical
			stn1_it->get()->GetmeridianDef());		// deflection in prime meridian

		// 4. apply deflection correction
		m_drValue += m_preAdjCorr;
	}

	// compute standard deviation in radians
	m_dStdDev = SecondstoRadians(0.01);

	m_epoch = "01.10.1985";
}
	

UINT32 CDnaAngle::SetMeasurementRec(const vstn_t& binaryStn, it_vmsr_t& it_msr, it_vdbid_t& dbidmap)
{
	m_strType = it_msr->measType;
	m_bIgnore = it_msr->ignore;
	m_MSmeasurementStations = (MEASUREMENT_STATIONS)it_msr->measurementStations;
	
	// first station
	m_lstn1Index = it_msr->station1;
	m_strFirst = binaryStn.at(it_msr->station1).stationName;

	// target station
	m_lstn2Index = it_msr->station2;
	m_strTarget = binaryStn.at(it_msr->station2).stationName;
	
	// target2 station
	m_lstn3Index = it_msr->station3;
	m_strTarget2 = binaryStn.at(it_msr->station3).stationName;
	
	m_measAdj = it_msr->measAdj;
	m_measCorr = it_msr->measCorr;
	m_measAdjPrec = it_msr->measAdjPrec;
	m_residualPrec = it_msr->residualPrec;
	m_preAdjCorr = it_msr->preAdjCorr;
	m_drValue = it_msr->term1;
	m_dStdDev = sqrt(it_msr->term2);
	
	m_epoch = it_msr->epoch;

	CDnaMeasurement::SetDatabaseMap(*dbidmap);

	return 0;
}

void CDnaAngle::WriteBinaryMsr(std::ofstream* binary_stream, PUINT32 msrIndex) const
{
	measurement_t measRecord;
	measRecord.measType = GetTypeC();
	measRecord.measStart = xMeas;
	measRecord.ignore = m_bIgnore;
	measRecord.station1 = m_lstn1Index;
	measRecord.station2 = m_lstn2Index;
	measRecord.station3 = m_lstn3Index;
	measRecord.measAdj = m_measAdj;
	measRecord.measCorr = m_measCorr;
	measRecord.measAdjPrec = m_measAdjPrec;
	measRecord.residualPrec = m_residualPrec;
	measRecord.preAdjCorr = m_preAdjCorr;
	measRecord.term1 = m_drValue;
	measRecord.term2 = m_dStdDev * m_dStdDev;	// convert to variance
	measRecord.measurementStations = m_MSmeasurementStations;
	measRecord.fileOrder = ((*msrIndex)++);

	snprintf(measRecord.epoch, sizeof(measRecord.epoch), "%s", m_epoch.substr(0, STN_EPOCH_WIDTH).c_str());

	binary_stream->write(reinterpret_cast<char *>(&measRecord), sizeof(measurement_t));
}

void CDnaAngle::SetValue(const std::string& str)
{
	// convert to radians
	RadFromDmsString(&m_drValue, trimstr(str));
}

void CDnaAngle::SetStdDev(const std::string& str)
{
	// convert seconds to radians
	RadFromSecondsString(&m_dStdDev, trimstr(str));
}

}	// namespace measurements
}	// namespace dynadjust
