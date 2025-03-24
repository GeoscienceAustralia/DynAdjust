//============================================================================
// Name         : dnaheight.cpp
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
// Description  : CDnaHeight implementation file
//============================================================================

#include <include/measurement_types/dnaheight.hpp>

namespace dynadjust {
namespace measurements {

CDnaHeight::CDnaHeight(void)
	: m_dValue(0.)
	, m_dStdDev(0.)
{
	m_strType = "H";		// default is orthometric distance, but could also be 'R'
	m_MSmeasurementStations = ONE_STATION;
}


CDnaHeight::~CDnaHeight(void)
{

}
	

bool CDnaHeight::operator== (const CDnaHeight& rhs) const
{
	return (
		m_strFirst == rhs.m_strFirst &&
		m_strType == rhs.m_strType &&
		m_bIgnore == rhs.m_bIgnore &&
		m_dValue == rhs.m_dValue &&
		m_dStdDev == rhs.m_dStdDev &&
		m_epoch == rhs.m_epoch
		);
}

bool CDnaHeight::operator< (const CDnaHeight& rhs) const
{
	if (m_strFirst == rhs.m_strFirst) {
		if (m_strType == rhs.m_strType) {	// don't think this is needed
			if (m_bIgnore == rhs.m_bIgnore) {
				if (m_epoch == rhs.m_epoch) {
					if (m_dValue == rhs.m_dValue)
						return m_dStdDev < rhs.m_dStdDev; 
					else
						return m_dValue < rhs.m_dValue; }
				else
					return m_epoch < rhs.m_epoch; }
			else
				return m_bIgnore < rhs.m_bIgnore; }
		else
			return m_strType < rhs.m_strType; }
	else
		return m_strFirst < rhs.m_strFirst;
}
	

void CDnaHeight::WriteDynaMLMsr(std::ofstream* dynaml_stream, const std::string& comment, bool) const
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
	*dynaml_stream << "    <Value>" << std::fixed << std::setprecision(4) << m_dValue << "</Value>" << std::endl;
	*dynaml_stream << "    <StdDev>" << std::scientific << std::setprecision(6) << m_dStdDev << "</StdDev>" << std::endl;
		
	if (m_msr_db_map.is_msr_id_set)
		*dynaml_stream << "    <MeasurementID>" << m_msr_db_map.msr_id << "</MeasurementID>" << std::endl;
	
	*dynaml_stream << "  </DnaMeasurement>" << std::endl;
}
	

void CDnaHeight::WriteDNAMsr(std::ofstream* dna_stream, const dna_msr_fields& dmw, const dna_msr_fields& dml, bool) const
{
	*dna_stream << std::setw(dmw.msr_type) << m_strType;
	if (m_bIgnore)
		*dna_stream << std::setw(dmw.msr_ignore) << "*";
	else
		*dna_stream << std::setw(dmw.msr_ignore) << " ";

	*dna_stream << std::left << std::setw(dmw.msr_inst) << m_strFirst;
	*dna_stream << std::setw(dmw.msr_targ1) << " ";
	*dna_stream << std::setw(dmw.msr_targ2) << " ";
	*dna_stream << std::right << std::setw(dmw.msr_linear) << std::fixed << std::setprecision(4) << m_dValue;	// linear measurement value
	*dna_stream << std::setw(dmw.msr_ang_d + dmw.msr_ang_m + dmw.msr_ang_s) << " ";
	
	UINT32 m_stdDevPrec(3);	// note change from 6 decimal places to 3
	*dna_stream << std::setw(dmw.msr_stddev) << StringFromTW(m_dStdDev, dmw.msr_stddev, m_stdDevPrec);
	//*dna_stream << std::setw(dmw.msr_stddev) << std::fixed << std::setprecision(6) << m_dStdDev;

	*dna_stream << std::setw(dml.msr_gps_epoch - dml.msr_inst_ht) << " ";
	*dna_stream << std::setw(dmw.msr_gps_epoch) << m_epoch;

	if (m_msr_db_map.is_msr_id_set)
		*dna_stream << std::setw(dmw.msr_id_msr) << m_msr_db_map.msr_id;

	*dna_stream << std::endl;
}
	

void CDnaHeight::SimulateMsr(vdnaStnPtr* vStations, const CDnaEllipsoid* ellipsoid)
{
	// Zn is the z coordinate element of the point on the z-axis 
	// which intersects with the the normal at the given Latitude
	double nu1, Zn1;

	// calculated height
	m_dValue = EllipsoidHeight<double>(
		vStations->at(m_lstn1Index).get()->GetXAxis(), 
		vStations->at(m_lstn1Index).get()->GetYAxis(),
		vStations->at(m_lstn1Index).get()->GetZAxis(),
		vStations->at(m_lstn1Index).get()->GetcurrentLatitude(),
		&nu1, &Zn1, ellipsoid);	

	m_dStdDev = 0.024;

	switch (GetTypeC())
	{
	case 'H':
		// N value available?
		if (fabs(vStations->at(m_lstn1Index).get()->GetgeoidSep()) > PRECISION_1E4)
		{
			// reduce to orthometric height
			m_preAdjCorr = vStations->at(m_lstn1Index).get()->GetgeoidSep();
			m_dValue -= m_preAdjCorr;
		}
		break;
	}

	m_epoch = "01.10.1985";
}
	

UINT32 CDnaHeight::SetMeasurementRec(const vstn_t& binaryStn, it_vmsr_t& it_msr, it_vdbid_t& dbidmap)
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
	m_dValue = it_msr->term1;
	m_dStdDev = sqrt(it_msr->term2);

	m_epoch = it_msr->epoch;

	CDnaMeasurement::SetDatabaseMap(*dbidmap);

	return 0;
}
	

void CDnaHeight::WriteBinaryMsr(std::ofstream* binary_stream, PUINT32 msrIndex) const
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
	measRecord.term1 = m_dValue;
	measRecord.term2 = m_dStdDev * m_dStdDev;	// convert to variance
	measRecord.term3 = 0.;
	measRecord.term4 = 0.;
	measRecord.measurementStations = m_MSmeasurementStations;
	measRecord.fileOrder = ((*msrIndex)++);

	snprintf(measRecord.epoch, sizeof(measRecord.epoch), "%s", m_epoch.substr(0, STN_EPOCH_WIDTH).c_str());

	binary_stream->write(reinterpret_cast<char *>(&measRecord), sizeof(measurement_t));
}


void CDnaHeight::SetValue(const std::string& str)
{
	DoubleFromString(m_dValue, trimstr(str));
}

void CDnaHeight::SetStdDev(const std::string& str)
{
	DoubleFromString(m_dStdDev, trimstr(str));
}


}	// namespace measurements
}	// namespace dynadjust
