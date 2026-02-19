//============================================================================
// Name         : dnadistance.cpp
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
// Description  : CDnaDistance implementation file
//============================================================================

#include <include/measurement_types/dnadistance.hpp>

namespace dynadjust {
namespace measurements {

CDnaDistance::CDnaDistance(void)
	: m_strTarget("")
	, m_dValue(0.)
	, m_dStdDev(0.)
	, m_fInstHeight(0.)
	, m_fTargHeight(0.)
{
	m_strType = "S";		// default is Slope distance, but could also be 'C' or 'E' or 'M'
	m_MSmeasurementStations = TWO_STATION;
}


CDnaDistance::~CDnaDistance(void)
{

}
	

bool CDnaDistance::operator== (const CDnaDistance& rhs) const
{
	return (
		m_strFirst == rhs.m_strFirst &&
		m_bIgnore == rhs.m_bIgnore &&
		m_strTarget == rhs.m_strTarget &&
		m_dValue == rhs.m_dValue &&
		m_dStdDev == rhs.m_dStdDev &&
		m_strType == rhs.m_strType &&
		m_fInstHeight == rhs.m_fInstHeight &&
		m_fTargHeight == rhs.m_fTargHeight &&
		m_epoch == rhs.m_epoch
		);
}


bool CDnaDistance::operator< (const CDnaDistance& rhs) const
{
	if (m_strFirst == rhs.m_strFirst) {
		if (m_strType == rhs.m_strType) {	// could be C, E, M, S
			if (m_bIgnore == rhs.m_bIgnore) {
				if (m_epoch == rhs.m_epoch) {
					if (m_strTarget == rhs.m_strTarget) {
						if (m_dValue == rhs.m_dValue) {
							if (m_dStdDev == rhs.m_dStdDev) {
								if (m_fInstHeight == rhs.m_fInstHeight)
									return m_fTargHeight < rhs.m_fTargHeight;
								else
									return m_fInstHeight < rhs.m_fInstHeight; }
							else
								return m_dStdDev < rhs.m_dStdDev; }
						else
							return m_dValue < rhs.m_dValue; }
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
	

void CDnaDistance::WriteDynaMLMsr(std::ofstream* dynaml_stream, const std::string& comment, bool) const
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
	*dynaml_stream << "    <Value>" << std::fixed << std::setprecision(4) << m_dValue << "</Value>" << std::endl;
	*dynaml_stream << "    <StdDev>" << std::fixed << std::setprecision(4) << m_dStdDev << "</StdDev>" << std::endl;

	switch (GetTypeC())
	{
	case 'S':
		*dynaml_stream << "    <InstHeight>" << std::fixed << std::setprecision(3) << m_fInstHeight << "</InstHeight>" << std::endl;
		*dynaml_stream << "    <TargHeight>" << std::fixed << std::setprecision(3) << m_fTargHeight << "</TargHeight>" << std::endl;
		break;
	}
		
	if (m_msr_db_map.is_msr_id_set)
		*dynaml_stream << "    <MeasurementID>" << m_msr_db_map.msr_id << "</MeasurementID>" << std::endl;
	
	*dynaml_stream << "  </DnaMeasurement>" << std::endl;
}
	

void CDnaDistance::WriteDNAMsr(std::ofstream* dna_stream, const dna_msr_fields& dmw, const dna_msr_fields& dml, bool) const
{
	*dna_stream << std::setw(dmw.msr_type) << m_strType;
	if (m_bIgnore)
		*dna_stream << std::setw(dmw.msr_ignore) << "*";
	else
		*dna_stream << std::setw(dmw.msr_ignore) << " ";

	*dna_stream << std::left << std::setw(dmw.msr_inst) << m_strFirst;
	*dna_stream << std::left << std::setw(dmw.msr_targ1) << m_strTarget;
	*dna_stream << std::setw(dmw.msr_targ2) << " ";
	*dna_stream << std::right << std::setw(dmw.msr_linear) << std::fixed << std::setprecision(4) << m_dValue;
	*dna_stream << std::setw(dmw.msr_ang_d + dmw.msr_ang_m + dmw.msr_ang_s) << " ";

	UINT32 m_stdDevPrec(3);
	*dna_stream << std::setw(dmw.msr_stddev) << StringFromTW(m_dStdDev, dmw.msr_stddev, m_stdDevPrec);
	//*dna_stream << std::setw(dmw.msr_stddev) << std::fixed << std::setprecision(3) << m_dStdDev;

	// database id width
	UINT32 width(dml.msr_gps_epoch - dml.msr_inst_ht);
	
	// Slope distance
	switch (GetTypeC())
	{
	case 'S':
		*dna_stream << std::setw(dmw.msr_inst_ht) << std::fixed << std::setprecision(3) << m_fInstHeight;
		*dna_stream << std::setw(dmw.msr_targ_ht) << std::fixed << std::setprecision(3) << m_fTargHeight;
		width = dml.msr_gps_epoch - dml.msr_targ_ht - dmw.msr_targ_ht;
		break;
	}
	
	*dna_stream << std::setw(width) << " ";
	*dna_stream << std::setw(dmw.msr_gps_epoch) << m_epoch;

	if (m_msr_db_map.is_msr_id_set)
		*dna_stream << std::setw(dmw.msr_id_msr) << m_msr_db_map.msr_id;

	*dna_stream << std::endl;
}
	

void CDnaDistance::SimulateMsr(vdnaStnPtr* vStations, const CDnaEllipsoid* ellipsoid)
{
	double dX, dY, dZ;
	double dXih, dYih, dZih, dXth, dYth, dZth;

	_it_vdnastnptr stn1_it(vStations->begin() + m_lstn1Index);
	_it_vdnastnptr stn2_it(vStations->begin() + m_lstn2Index);

	m_fInstHeight = float(0.000);
	m_fTargHeight = float(0.000);

	switch (GetTypeC())
	{
	case 'S':
		CartesianElementsFromInstrumentHeight<double>(m_fInstHeight,				// instrument height
			&dXih, &dYih, &dZih, 
			stn1_it->get()->GetcurrentLatitude(),
			stn1_it->get()->GetcurrentLongitude());
		CartesianElementsFromInstrumentHeight<double>(m_fTargHeight,				// target height
			&dXth, &dYth, &dZth,
			stn1_it->get()->GetcurrentLatitude(),
			stn1_it->get()->GetcurrentLongitude());
		dX = (stn2_it->get()->GetXAxis() - stn1_it->get()->GetXAxis() + dXth - dXih);
		dY = (stn2_it->get()->GetYAxis() - stn1_it->get()->GetYAxis() + dYth - dYih);
		dZ = (stn2_it->get()->GetZAxis() - stn1_it->get()->GetZAxis() + dZth - dZih);

		m_dValue = magnitude(dX, dY, dZ);
		break;
	case 'C':
		m_dValue = EllipsoidChordDistance<double>(
			stn1_it->get()->GetXAxis(),
			stn1_it->get()->GetYAxis(),
			stn1_it->get()->GetZAxis(),
			stn2_it->get()->GetXAxis(),
			stn2_it->get()->GetYAxis(),
			stn2_it->get()->GetZAxis(),
			stn1_it->get()->GetcurrentLatitude(),
			stn2_it->get()->GetcurrentLatitude(),
			stn1_it->get()->GetcurrentHeight(),
			stn2_it->get()->GetcurrentHeight(),
			&dX, &dY, &dZ, ellipsoid);
		break;
	case 'E':		
		m_dValue = EllipsoidArcDistance<double>(
			stn1_it->get()->GetXAxis(),
			stn1_it->get()->GetYAxis(),
			stn1_it->get()->GetZAxis(),
			stn2_it->get()->GetXAxis(),
			stn2_it->get()->GetYAxis(),
			stn2_it->get()->GetZAxis(),
			stn1_it->get()->GetcurrentLatitude(),
			stn1_it->get()->GetcurrentLongitude(),
			stn2_it->get()->GetcurrentLatitude(),
			stn1_it->get()->GetcurrentHeight(),
			stn2_it->get()->GetcurrentHeight(),
			ellipsoid);
		break;
	case 'M':		
		m_dValue = MSLArcDistance<double>(
			stn1_it->get()->GetXAxis(),
			stn1_it->get()->GetYAxis(),
			stn1_it->get()->GetZAxis(),
			stn2_it->get()->GetXAxis(),
			stn2_it->get()->GetYAxis(),
			stn2_it->get()->GetZAxis(),
			stn1_it->get()->GetcurrentLatitude(),
			stn2_it->get()->GetcurrentLatitude(),
			stn1_it->get()->GetcurrentHeight(),
			stn2_it->get()->GetcurrentHeight(),
			stn1_it->get()->GetgeoidSep(),
			stn2_it->get()->GetgeoidSep(),
			ellipsoid);
		break;
	}

	m_dStdDev = 3.0 * sqrt(m_dValue / 1000.0) / 100.0;

	m_epoch = "01.10.1985";
}
	

UINT32 CDnaDistance::SetMeasurementRec(const vstn_t& binaryStn, it_vmsr_t& it_msr, it_vdbid_t& dbidmap)
{
	m_bIgnore = it_msr->ignore;
	m_MSmeasurementStations = (MEASUREMENT_STATIONS)it_msr->measurementStations;
	
	m_strType = it_msr->measType;
	
	// first station
	m_lstn1Index = it_msr->station1;
	m_strFirst = binaryStn.at(it_msr->station1).stationName;
	
	if (GetTypeC() == 'S')		// slope distance
	{
		// inst and targ heights
		m_fInstHeight = static_cast<float> (it_msr->term3);
		m_fTargHeight = static_cast<float> (it_msr->term4);
	}

	// target station
	m_lstn2Index = it_msr->station2;
	m_strTarget = binaryStn.at(it_msr->station2).stationName;
	
	m_measAdj = it_msr->measAdj;
	m_measCorr = it_msr->measCorr;
	m_measAdjPrec = it_msr->measAdjPrec;
	m_residualPrec = it_msr->residualPrec;
	m_preAdjCorr = it_msr->preAdjCorr;
	m_dValue = it_msr->term1;
	m_dStdDev = sqrt(it_msr->term2);

	m_epoch = it_msr->epoch;
	m_sourceFileIndex = it_msr->sourceFileIndex;

	CDnaMeasurement::SetDatabaseMap(*dbidmap);

	return 0;
}


void CDnaDistance::WriteBinaryMsr(std::ofstream* binary_stream, PUINT32 msrIndex) const
{
	measurement_t measRecord;
	measRecord.measType = GetTypeC();
	measRecord.measStart = xMeas;
	measRecord.ignore = m_bIgnore;
	measRecord.station1 = m_lstn1Index;
	measRecord.station2 = m_lstn2Index;
	measRecord.measAdj = m_measAdj;
	measRecord.measCorr = m_measCorr;
	measRecord.measAdjPrec = m_measAdjPrec;
	measRecord.residualPrec = m_residualPrec;
	measRecord.preAdjCorr = m_preAdjCorr;
	measRecord.term1 = m_dValue;
	measRecord.term2 = m_dStdDev * m_dStdDev;	// convert to variance
	measRecord.term3 = m_fInstHeight;
	measRecord.term4 = m_fTargHeight;
	measRecord.measurementStations = m_MSmeasurementStations;
	measRecord.sourceFileIndex = m_sourceFileIndex;
	measRecord.fileOrder = ((*msrIndex)++);

	snprintf(measRecord.epoch, sizeof(measRecord.epoch), "%s", m_epoch.substr(0, STN_EPOCH_WIDTH).c_str());

	binary_stream->write(reinterpret_cast<char *>(&measRecord), sizeof(measurement_t));
}


void CDnaDistance::SetValue(const std::string& str)
{
	DoubleFromString(m_dValue, trimstr(str));
}

void CDnaDistance::SetStdDev(const std::string& str)
{
	DoubleFromString(m_dStdDev, trimstr(str));
}

// Instrument and target heights only make sense for 
// slope distances, vertical angles and zenith distances.
// Note: slope distances are derived from CDnaDistance, so
// these methods are provided here, but don't need 
// instrument height for 'C', 'E' and 'M'  measurements.
void CDnaDistance::SetInstrumentHeight(const std::string& str)
{
	FloatFromString(m_fInstHeight, trimstr(str));
}

void CDnaDistance::SetTargetHeight(const std::string& str)
{
	DoubleFromString(m_fTargHeight, trimstr(str));
}


}	// namespace measurements
}	// namespace dynadjust
