//============================================================================
// Name         : dnacoordinate.cpp
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
// Description  : CDnaCoordinate implementation file
//============================================================================

#include <include/measurement_types/dnacoordinate.hpp>

namespace dynadjust {
namespace measurements {

CDnaCoordinate::CDnaCoordinate(void)
	: m_drValue(0.)
	, m_dStdDev(0.)
{
	m_MSmeasurementStations = ONE_STATION;
}


CDnaCoordinate::~CDnaCoordinate(void)
{

}
	

bool CDnaCoordinate::operator== (const CDnaCoordinate& rhs) const
{
	return (
		m_strFirst == rhs.m_strFirst &&
		m_strType == rhs.m_strType &&
		m_bIgnore == rhs.m_bIgnore &&
		m_drValue == rhs.m_drValue &&
		m_dStdDev == rhs.m_dStdDev &&
		m_epoch == rhs.m_epoch
		);
}

bool CDnaCoordinate::operator< (const CDnaCoordinate& rhs) const
{
	if (m_strFirst == rhs.m_strFirst) {
		if (m_strType == rhs.m_strType) {	// could be P, Q, J, K
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
			return m_strType < rhs.m_strType; }
	else
		return m_strFirst < rhs.m_strFirst;
}
	

void CDnaCoordinate::WriteDynaMLMsr(std::ofstream* dynaml_stream, const std::string& comment, bool) const
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
	
	*dynaml_stream << "    <Value>" << std::fixed << std::setprecision(10);
	switch (GetTypeC())
	{
	case 'P':
	case 'I':		// Latitude requires Degrees(...)
		*dynaml_stream << RadtoDms(m_drValue);
		break;
	case 'Q':
	case 'J':		// Longitude requires DegreesL(...)
		*dynaml_stream << RadtoDmsL(m_drValue);
		break;
	}

	*dynaml_stream << "</Value>" << std::endl;
	
	*dynaml_stream << "    <StdDev>" << std::scientific << std::setprecision(6) << Seconds(m_dStdDev) << "</StdDev>" << std::endl;
		
	if (m_msr_db_map.is_msr_id_set)
		*dynaml_stream << "    <MeasurementID>" << m_msr_db_map.msr_id << "</MeasurementID>" << std::endl;
	
	*dynaml_stream << "  </DnaMeasurement>" << std::endl;
}
	

void CDnaCoordinate::WriteDNAMsr(std::ofstream* dna_stream, const dna_msr_fields& dmw, const dna_msr_fields& dml, bool) const
{
	*dna_stream << std::setw(dmw.msr_type) << m_strType;
	if (m_bIgnore)
		*dna_stream << std::setw(dmw.msr_ignore) << "*";
	else
		*dna_stream << std::setw(dmw.msr_ignore) << " ";
	*dna_stream << std::left << std::setw(dmw.msr_inst) << m_strFirst;
	*dna_stream << std::setw(dmw.msr_targ1) << " ";	// first target
	*dna_stream << std::setw(dmw.msr_targ2) << " ";	// second target
	*dna_stream << std::right << std::setw(dmw.msr_linear) << std::fixed << std::setprecision(9) << RadtoDms(m_drValue);	// linear measurement value
	*dna_stream << std::setw(dmw.msr_ang_d + dmw.msr_ang_m + dmw.msr_ang_s) << " ";

	UINT32 m_stdDevPrec(6);
	*dna_stream << std::setw(dmw.msr_stddev) << StringFromTW(Seconds(m_dStdDev), dmw.msr_stddev, m_stdDevPrec);
	//*dna_stream << std::setw(dmw.msr_stddev) << std::fixed << std::setprecision(6) << Seconds(m_dStdDev);

	*dna_stream << std::setw(dml.msr_gps_epoch - dml.msr_inst_ht) << " ";
	*dna_stream << std::setw(dmw.msr_gps_epoch) << m_epoch;

	if (m_msr_db_map.is_msr_id_set)
		*dna_stream << std::setw(dmw.msr_id_msr) << m_msr_db_map.msr_id;

	*dna_stream << std::endl;
}
	

void CDnaCoordinate::SimulateMsr(vdnaStnPtr* vStations, const CDnaEllipsoid*)
{
	switch (GetTypeC())
	{
	case 'P':
	case 'I':
		m_drValue = vStations->at(m_lstn1Index).get()->GetcurrentLatitude();
		break;
	case 'Q':
	case 'J':
		m_drValue = vStations->at(m_lstn1Index).get()->GetcurrentLongitude();
		break;
	}

	switch (GetTypeC())
	{
	case 'I':
		// convert to astro (requires geoid)
		if (fabs(vStations->at(m_lstn1Index).get()->GetmeridianDef()) > E4_SEC_DEFLECTION)
			m_drValue = vStations->at(m_lstn1Index).get()->GetcurrentLatitude() +
				vStations->at(m_lstn1Index).get()->GetmeridianDef();
		break;
	case 'J':
		// convert to astro (requires geoid)
		if (fabs(vStations->at(m_lstn1Index).get()->GetverticalDef()) > E4_SEC_DEFLECTION)
			m_drValue = vStations->at(m_lstn1Index).get()->GetcurrentLongitude() +
				vStations->at(m_lstn1Index).get()->GetverticalDef() / 
				cos(vStations->at(m_lstn1Index).get()->GetcurrentLatitude());
		break;
	}

	// compute standard deviation in radians
	m_dStdDev = SecondstoRadians(0.021);

	m_epoch = "01.10.1985";
}
	

UINT32 CDnaCoordinate::SetMeasurementRec(const vstn_t& binaryStn, it_vmsr_t& it_msr, it_vdbid_t& dbidmap)
{
	m_strType = it_msr->measType;
	m_bIgnore = it_msr->ignore;
	m_MSmeasurementStations = (MEASUREMENT_STATIONS)it_msr->measurementStations;
	
	m_lstn1Index = it_msr->station1;
	m_strFirst = binaryStn.at(it_msr->station1).stationName;
	
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

void CDnaCoordinate::WriteBinaryMsr(std::ofstream* binary_stream, PUINT32 msrIndex) const
{
	measurement_t measRecord;
	measRecord.measType = GetTypeC();
	measRecord.measStart = xMeas;
	measRecord.ignore = m_bIgnore;
	measRecord.station1 = m_lstn1Index;
	measRecord.measAdj = m_measAdj;
	measRecord.measCorr = m_measCorr;
	measRecord.measAdjPrec = m_measAdjPrec;
	measRecord.residualPrec = m_residualPrec;
	measRecord.preAdjCorr = m_preAdjCorr;
	measRecord.term1 = m_drValue;
	measRecord.term2 = m_dStdDev * m_dStdDev;	// convert to variance
	measRecord.term3 = measRecord.term4 = 0.;
	measRecord.measurementStations = m_MSmeasurementStations;
	measRecord.fileOrder = ((*msrIndex)++);

	snprintf(measRecord.epoch, sizeof(measRecord.epoch), "%s", m_epoch.substr(0, STN_EPOCH_WIDTH).c_str());

	binary_stream->write(reinterpret_cast<char *>(&measRecord), sizeof(measurement_t));
}


void CDnaCoordinate::SetValue(const std::string& str)
{
	RadFromDmsString(&m_drValue, trimstr(str));
}

void CDnaCoordinate::SetStdDev(const std::string& str)
{
	RadFromSecondsString(&m_dStdDev, trimstr(str));			// convert to radians
}


}	// namespace measurements
}	// namespace dynadjust
