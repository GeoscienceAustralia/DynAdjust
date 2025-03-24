//============================================================================
// Name         : dnaheightdifference.cpp
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
// Description  : CDnaHeightDifference implementation file
//============================================================================

#include <include/measurement_types/dnaheightdifference.hpp>

//extern boost::random::mt19937 rng;
//extern boost::random::uniform_real_distribution<double> stdev;
//extern boost::random::uniform_real_distribution<double> pertu;

namespace dynadjust {
namespace measurements {

CDnaHeightDifference::CDnaHeightDifference(void)
	: m_strTarget("")
	, m_dValue(0.)
	, m_dStdDev(0.)
{
	m_strType = "L";		// default is Slope distance, but could also be 'C' or 'E' or 'M'
	m_MSmeasurementStations = TWO_STATION;
}


CDnaHeightDifference::~CDnaHeightDifference(void)
{

}
	

bool CDnaHeightDifference::operator== (const CDnaHeightDifference& rhs) const
{
	return (
		m_strFirst == rhs.m_strFirst &&
		m_strTarget == rhs.m_strTarget &&
		m_strType == rhs.m_strType &&
		m_bIgnore == rhs.m_bIgnore &&
		m_dValue == rhs.m_dValue &&
		m_dStdDev == rhs.m_dStdDev &&
		m_epoch == rhs.m_epoch
		);
}


bool CDnaHeightDifference::operator< (const CDnaHeightDifference& rhs) const
{
	if (m_strFirst == rhs.m_strFirst) {
		if (m_strType == rhs.m_strType) {	// don't think this is needed
			if (m_bIgnore == rhs.m_bIgnore) {
				if (m_epoch == rhs.m_epoch) {
					if (m_strTarget == rhs.m_strTarget) {
						if (m_dValue == rhs.m_dValue)
							return m_dStdDev < rhs.m_dStdDev; 
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
	

void CDnaHeightDifference::WriteDynaMLMsr(std::ofstream* dynaml_stream, const std::string& comment, bool) const
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
	*dynaml_stream << "    <StdDev>" << std::fixed << std::setprecision(6) << m_dStdDev << "</StdDev>" << std::endl;
	
	if (m_msr_db_map.is_msr_id_set)
		*dynaml_stream << "    <MeasurementID>" << m_msr_db_map.msr_id << "</MeasurementID>" << std::endl;

	*dynaml_stream << "  </DnaMeasurement>" << std::endl;
}


void CDnaHeightDifference::WriteDNAMsr(std::ofstream* dna_stream, const dna_msr_fields& dmw, const dna_msr_fields& dml, bool) const
{
	*dna_stream << std::setw(dmw.msr_type) << m_strType;
	if (m_bIgnore)
		*dna_stream << std::setw(dmw.msr_ignore) << "*";
	else
		*dna_stream << std::setw(dmw.msr_ignore) << " ";

	*dna_stream << std::left << std::setw(dmw.msr_inst) << m_strFirst;
	*dna_stream << std::left << std::setw(dmw.msr_targ1) << m_strTarget;
	*dna_stream << std::setw(dmw.msr_targ2) << " ";
	*dna_stream << std::right << std::setw(dmw.msr_linear) << std::fixed << std::setprecision(4) << m_dValue;	// linear measurement value
	*dna_stream << std::setw(dmw.msr_ang_d + dmw.msr_ang_m + dmw.msr_ang_s) << " ";
	
	UINT32 m_stdDevPrec(3); // note change from 6 decimal places to 3
	*dna_stream << std::setw(dmw.msr_stddev) << StringFromTW(m_dStdDev, dmw.msr_stddev, m_stdDevPrec);
	//*dna_stream << std::setw(dmw.msr_stddev) << std::fixed << std::setprecision(6) << m_dStdDev;

	*dna_stream << std::setw(dml.msr_gps_epoch - dml.msr_inst_ht) << " ";
	*dna_stream << std::setw(dmw.msr_gps_epoch) << m_epoch;

	if (m_msr_db_map.is_msr_id_set)
		*dna_stream << std::setw(dmw.msr_id_msr) << m_msr_db_map.msr_id;

	*dna_stream << std::endl;
}
	

void CDnaHeightDifference::SimulateMsr(vdnaStnPtr* vStations, const CDnaEllipsoid* ellipsoid)
{
	// Zn is the z coordinate element of the point on the z-axis 
	// which intersects with the the normal at the given Latitude
	double h1, h2, nu1, Zn1, nu2, Zn2;

	// calculated height
	m_dValue = EllipsoidHeightDifference<double>(
		vStations->at(m_lstn1Index).get()->GetXAxis(), 
		vStations->at(m_lstn1Index).get()->GetYAxis(),
		vStations->at(m_lstn1Index).get()->GetZAxis(),
		vStations->at(m_lstn2Index).get()->GetXAxis(), 
		vStations->at(m_lstn2Index).get()->GetYAxis(),
		vStations->at(m_lstn2Index).get()->GetZAxis(),
		vStations->at(m_lstn1Index).get()->GetcurrentLatitude(),
		vStations->at(m_lstn2Index).get()->GetcurrentLatitude(),
		&h1, &h2, &nu1, &Zn1, &nu2, &Zn2, ellipsoid);

	double distance = GreatCircleDistance<double>(
		vStations->at(m_lstn1Index).get()->GetcurrentLatitude(),
		vStations->at(m_lstn1Index).get()->GetcurrentLongitude(),
		vStations->at(m_lstn2Index).get()->GetcurrentLatitude(),
		vStations->at(m_lstn2Index).get()->GetcurrentLongitude());

	m_dStdDev = 3.0 * sqrt(distance / 1000.0) / 100.0;

	m_epoch = "01.10.1985";

	// N value available?
	if (fabs(vStations->at(m_lstn1Index).get()->GetgeoidSep()) > PRECISION_1E4 ||
		fabs(vStations->at(m_lstn2Index).get()->GetgeoidSep()) > PRECISION_1E4)
	{
		// reduce to orthometric height difference
		m_preAdjCorr = vStations->at(m_lstn2Index).get()->GetgeoidSep() - vStations->at(m_lstn1Index).get()->GetgeoidSep();
		m_dValue -= m_preAdjCorr;
	}

    // TODO - add option to perturb all measurements and standard deviations by a small (random) amount
	//      - Below is an example of how this might be achieved
    //
	//// perturb standard deviation
	//m_dStdDev = static_cast<double>(stdev(rng)) * sqrt(distance / 1000.0) / 100.0;
	//
	//// perturb measurement
	//double whichway = static_cast<double>(pertu(rng));
	//if (whichway > 0.5)
	//	m_dValue -= static_cast<double>(pertu(rng))/10.0;
	//else
	//	m_dValue += static_cast<double>(pertu(rng))/10.0;
}
	

UINT32 CDnaHeightDifference::SetMeasurementRec(const vstn_t& binaryStn, it_vmsr_t& it_msr, it_vdbid_t& dbidmap)
{
	m_bIgnore = it_msr->ignore;
	m_MSmeasurementStations = (MEASUREMENT_STATIONS)it_msr->measurementStations;
	
	m_strType = it_msr->measType;
	
	// first station
	m_lstn1Index = it_msr->station1;
	m_strFirst = binaryStn.at(it_msr->station1).stationName;

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

	CDnaMeasurement::SetDatabaseMap(*dbidmap);

	return 0;
}
	

void CDnaHeightDifference::WriteBinaryMsr(std::ofstream* binary_stream, PUINT32 msrIndex) const
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
	measRecord.measurementStations = m_MSmeasurementStations;
	measRecord.fileOrder = ((*msrIndex)++);

	snprintf(measRecord.epoch, sizeof(measRecord.epoch), "%s", m_epoch.substr(0, STN_EPOCH_WIDTH).c_str());

	binary_stream->write(reinterpret_cast<char *>(&measRecord), sizeof(measurement_t));
}


void CDnaHeightDifference::SetValue(const std::string& str)
{
	DoubleFromString(m_dValue, trimstr(str));
}

void CDnaHeightDifference::SetStdDev(const std::string& str)
{
	DoubleFromString(m_dStdDev, trimstr(str));
}

}	// namespace measurements
}	// namespace dynadjust
