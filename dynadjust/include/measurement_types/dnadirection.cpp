//============================================================================
// Name         : dnadirection.cpp
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
// Description  : CDnaDirection implementation file
//============================================================================

#include <include/exception/dnaexception.hpp>
#include <include/measurement_types/dnadirection.hpp>

using namespace dynadjust::exception;

namespace dynadjust {
namespace measurements {

CDnaDirection::CDnaDirection(void)
	: m_strTarget("")
	, m_drValue(0.)
	, m_dStdDev(0.)
	, m_fInstHeight(0.)
	, m_fTargHeight(0.)
	, m_lRecordedTotal(0)
	, m_lsetID(0)
{
	m_strType = "D";
	m_MSmeasurementStations = TWO_STATION;
}

CDnaDirection::~CDnaDirection(void)
{
}

// move constructor
CDnaDirection::CDnaDirection(CDnaDirection&& d)
{
	m_strFirst = d.m_strFirst;
	m_bIgnore = d.m_bIgnore;
	m_bInsufficient = d.m_bInsufficient;
	m_strTarget = d.m_strTarget;
	m_drValue = d.m_drValue;
	m_dStdDev = d.m_dStdDev;
	m_fInstHeight = d.m_fInstHeight;
	m_fTargHeight = d.m_fTargHeight;
	m_lRecordedTotal = d.m_lRecordedTotal;
	m_lsetID = d.m_lsetID;
	m_MSmeasurementStations = d.m_MSmeasurementStations;

	m_strType = "D";

	m_msr_db_map = d.m_msr_db_map;

	m_epoch = d.m_epoch;
}

// move assignment operator 
CDnaDirection& CDnaDirection::operator= (CDnaDirection&& rhs)
{
	// check for assignment to self!
	if (this == &rhs)
		return *this;

	CDnaMeasurement::operator=(std::move(rhs));
	m_bInsufficient = rhs.m_bInsufficient;
	m_strTarget = rhs.m_strTarget;
	m_drValue = rhs.m_drValue;
	m_dStdDev = rhs.m_dStdDev;
	m_fInstHeight = rhs.m_fInstHeight;
	m_fTargHeight = rhs.m_fTargHeight;
	m_lRecordedTotal = rhs.m_lRecordedTotal;
	m_lsetID = rhs.m_lsetID;
	m_MSmeasurementStations = rhs.m_MSmeasurementStations;

	m_msr_db_map = rhs.m_msr_db_map;

	return *this;
}






//// copy constructor (disabled)
//CDnaDirection::CDnaDirection(const CDnaDirection& newDirection)
//{
//	m_strFirst = newDirection.m_strFirst;
//	m_bIgnore = newDirection.m_bIgnore;
//	m_strTarget = newDirection.m_strTarget;
//	m_drValue = newDirection.m_drValue;
//	m_dStdDev = newDirection.m_dStdDev;
//	m_fInstHeight = newDirection.m_fInstHeight;
//	m_fTargHeight = newDirection.m_fTargHeight;
//	m_lRecordedTotal = newDirection.m_lRecordedTotal;
//	m_lsetID = newDirection.m_lsetID;
//	m_MSmeasurementStations = newDirection.m_MSmeasurementStations;
//
//	m_strType = "D";
//
//	m_msr_db_map = newDirection.m_msr_db_map;
//}


//CDnaDirection::CDnaDirection(const bool bIgnore, const std::string& strFirst, const std::string& strTarget, const double& drValue, const double& dStdDev)
//{
//	m_strType = "D";
//	m_strFirst = strFirst;
//	m_bIgnore = bIgnore;
//	m_strTarget = strTarget;
//	m_drValue = drValue;
//	m_dStdDev = dStdDev;
//	m_lRecordedTotal = 0;
//	m_lsetID = 0;
//}

// assignment operator (disabled)
//CDnaDirection& CDnaDirection::operator= (const CDnaDirection& rhs)
//{
//	// check for assignment to self!
//	if (this == &rhs)
//		return *this;
//
//	CDnaMeasurement::operator=(rhs);
//	m_strTarget = rhs.m_strTarget;
//	m_drValue = rhs.m_drValue;
//	m_dStdDev = rhs.m_dStdDev;
//	m_fInstHeight = rhs.m_fInstHeight;
//	m_fTargHeight = rhs.m_fTargHeight;
//	m_lRecordedTotal = rhs.m_lRecordedTotal;
//	m_lsetID = rhs.m_lsetID;
//	m_MSmeasurementStations = rhs.m_MSmeasurementStations;
//
//	m_msr_db_map = rhs.m_msr_db_map;
//
//	return *this;
//}


bool CDnaDirection::operator== (const CDnaDirection& rhs) const
{
	return (
		m_strFirst == rhs.m_strFirst &&
		m_bIgnore == rhs.m_bIgnore &&
		m_strTarget == rhs.m_strTarget &&
		m_drValue == rhs.m_drValue &&
		m_dStdDev == rhs.m_dStdDev &&
		m_strType == rhs.m_strType &&
		m_fInstHeight == rhs.m_fInstHeight &&
		m_fTargHeight == rhs.m_fTargHeight &&
		m_epoch == rhs.m_epoch
		);
}


bool CDnaDirection::operator< (const CDnaDirection& rhs) const
{
	if (m_strFirst == rhs.m_strFirst) {
		if (m_strType == rhs.m_strType) {	// could be B, D, V, Z
			if (m_bIgnore == rhs.m_bIgnore) {
				if (m_epoch == rhs.m_epoch) {
					if (m_strTarget == rhs.m_strTarget) {
						if (m_drValue == rhs.m_drValue) {
							if (m_dStdDev == rhs.m_dStdDev) {
								if (m_fInstHeight == rhs.m_fInstHeight)
									return m_fTargHeight < rhs.m_fTargHeight;
								else
									return m_fInstHeight < rhs.m_fInstHeight; }
							else
								return m_dStdDev < rhs.m_dStdDev; }
						else
							return m_drValue < rhs.m_drValue; }
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
	

void CDnaDirection::WriteDynaMLMsr(std::ofstream* dynaml_stream, const std::string& comment, bool bSubMeasurement /*= false*/) const
{
	if (bSubMeasurement)
	{
		// Direction set
		*dynaml_stream << "    <Directions>" << std::endl;
		if (m_bIgnore)
			*dynaml_stream << "      <Ignore>*</Ignore>" << std::endl;
		else
			*dynaml_stream << "      <Ignore/>" << std::endl;
		*dynaml_stream << "      <Target>" << m_strTarget << "</Target>" << std::endl;
		*dynaml_stream << "      <Value>" << std::setprecision(8) << std::fixed 
			<< RadtoDms(m_drValue) << "</Value>" << std::endl;
		*dynaml_stream << "      <StdDev>" << std::scientific << std::setprecision(6) << Seconds(m_dStdDev) << "</StdDev>" << std::endl;	
		
		if (m_msr_db_map.is_msr_id_set)
			*dynaml_stream << "      <MeasurementID>" << m_msr_db_map.msr_id << "</MeasurementID>" << std::endl;
		
		*dynaml_stream << "    </Directions>" << std::endl;
	}
	else
	{
		// Zenith distance, vertical angle
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
		*dynaml_stream << "    <Value>" << std::setprecision(8) << std::fixed 
			<< RadtoDms(m_drValue) << "</Value>" << std::endl;
		*dynaml_stream << "    <StdDev>" << std::scientific << std::setprecision(6) << Seconds(m_dStdDev) << "</StdDev>" << std::endl;	

		// Zenith distance, vertical angle
		switch (GetTypeC())
		{
		case 'Z':
		case 'V':
			*dynaml_stream << "    <InstHeight>" << std::fixed << std::setprecision(3) << m_fInstHeight << "</InstHeight>" << std::endl;
			*dynaml_stream << "    <TargHeight>" << std::fixed << std::setprecision(3) << m_fTargHeight << "</TargHeight>" << std::endl;
			break;
		}
		
		if (m_msr_db_map.is_msr_id_set)
			*dynaml_stream << "    <MeasurementID>" << m_msr_db_map.msr_id << "</MeasurementID>" << std::endl;
		
		*dynaml_stream << "  </DnaMeasurement>" << std::endl;
	}
}
	
void CDnaDirection::WriteDNAMsr(std::ofstream* dna_stream, const dna_msr_fields& dmw, const dna_msr_fields& dml, bool bSubMeasurement /*= false*/) const
{
	*dna_stream << std::setw(dmw.msr_type) << m_strType;
	if (m_bIgnore)
		*dna_stream << std::setw(dmw.msr_ignore) << "*";
	else
		*dna_stream << std::setw(dmw.msr_ignore) << " ";

	// database id width
	UINT32 width(dml.msr_id_msr - dml.msr_stddev - dmw.msr_stddev);

	if (bSubMeasurement)
	{
		// Direction set
		*dna_stream << std::left << std::setw(dmw.msr_inst) << " ";
		*dna_stream << std::setw(dmw.msr_targ1) << " ";
		*dna_stream << std::left << std::setw(dmw.msr_targ2) << m_strTarget;
		*dna_stream << std::setw(dmw.msr_linear) << " ";	// linear measurement value
		*dna_stream << std::setw(dmw.msr_ang_d + dmw.msr_ang_m + dmw.msr_ang_s) << 
			std::right << FormatDnaDmsString(RadtoDms(m_drValue), 8);
		
		UINT32 m_stdDevPrec(3);
		*dna_stream << std::setw(dmw.msr_stddev) << StringFromTW(Seconds(m_dStdDev), dmw.msr_stddev, m_stdDevPrec);
		//*dna_stream << std::setw(dmw.msr_stddev) << std::fixed << std::setprecision(3) << Seconds(m_dStdDev);
	}
	else
	{
		// Azimuth, zenith distance, vertical angle
		*dna_stream << std::left << std::setw(dmw.msr_inst) << m_strFirst;
		*dna_stream << std::left << std::setw(dmw.msr_targ1) << m_strTarget;
		*dna_stream << std::setw(dmw.msr_targ2) << " ";
		*dna_stream << std::setw(dmw.msr_linear) << " ";	// linear measurement value
		*dna_stream << std::setw(dmw.msr_ang_d + dmw.msr_ang_m + dmw.msr_ang_s) << 
			std::right << FormatDnaDmsString(RadtoDms(m_drValue), 8);
		
		UINT32 m_stdDevPrec(3);
		*dna_stream << std::setw(dmw.msr_stddev) << StringFromTW(Seconds(m_dStdDev), dmw.msr_stddev, m_stdDevPrec);
		//*dna_stream << std::setw(dmw.msr_stddev) << std::fixed << std::setprecision(3) << Seconds(m_dStdDev);
		width = dml.msr_gps_epoch - dml.msr_inst_ht;

		// Zenith distance, vertical angle
		switch (GetTypeC())
		{
		case 'Z':
		case 'V':
			*dna_stream << std::setw(dmw.msr_inst_ht) << std::fixed << std::setprecision(3) << m_fInstHeight;
			*dna_stream << std::setw(dmw.msr_targ_ht) << std::fixed << std::setprecision(3) << m_fTargHeight;
			width = dml.msr_gps_epoch - dml.msr_targ_ht - dmw.msr_targ_ht;
			break;
		}

		*dna_stream << std::setw(width) << " ";
		*dna_stream << std::setw(dmw.msr_gps_epoch) << m_epoch;
		width = 0;
	}

	if (m_msr_db_map.is_msr_id_set)
	{
		if (width)
			*dna_stream << std::setw(width) << " ";
		*dna_stream << std::setw(dmw.msr_id_msr) << m_msr_db_map.msr_id;

		if (bSubMeasurement)
			if (m_msr_db_map.is_cls_id_set)
				*dna_stream << std::setw(dmw.msr_id_cluster) << m_msr_db_map.cluster_id;
	}

	*dna_stream << std::endl;
}
	
void CDnaDirection::SimulateMsr(vdnaStnPtr* vStations, const CDnaEllipsoid*)
{
	_it_vdnastnptr stn1_it(vStations->begin() + m_lstn1Index);
	_it_vdnastnptr stn2_it(vStations->begin() + m_lstn2Index);
	
	double zenithDistance;
	double direction;

	direction = Direction(
		stn1_it->get()->GetXAxis(), 
		stn1_it->get()->GetYAxis(),
		stn1_it->get()->GetZAxis(),
		stn2_it->get()->GetXAxis(), 
		stn2_it->get()->GetYAxis(),
		stn2_it->get()->GetZAxis(),
		stn1_it->get()->GetcurrentLatitude(),
		stn1_it->get()->GetcurrentLongitude());

	switch (GetTypeC())
	{
	
	// direction
	case 'D':	
		// assign direction
		m_drValue = direction;

		// Deflections available?
		if (fabs(stn1_it->get()->GetverticalDef()) > E4_SEC_DEFLECTION || fabs(stn1_it->get()->GetmeridianDef()) > E4_SEC_DEFLECTION)
		{
			// 1. compute zenith distance 
			zenithDistance = ZenithDistance<double>(
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
				0.0,				
				0.0);	

			// 2. Compute deflection
			m_preAdjCorr = DirectionDeflectionCorrection<double>(
				m_drValue,
				zenithDistance,
				stn1_it->get()->GetverticalDef(),			// deflection in prime vertical
				stn1_it->get()->GetmeridianDef());			// deflection in prime meridian
			
			// 3. apply deflection correction
			m_drValue += m_preAdjCorr;
		}
		break;
	
	// geodetic azimuth
	case 'B':	
		// assign direction
		m_drValue = direction;
		m_epoch = "01.10.1985";
		break;
	
	// astronomic azimuth
	case 'K':	
		// assign direction
		m_drValue = direction;

		// Deflections available?
		if (fabs(stn1_it->get()->GetverticalDef()) > E4_SEC_DEFLECTION || fabs(stn1_it->get()->GetmeridianDef()) > E4_SEC_DEFLECTION)
		{
			// 1. compute zenith distance 
			zenithDistance = ZenithDistance<double>(
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
				0.0,				
				0.0);	

			// 2. Compute correction
			m_preAdjCorr = LaplaceCorrection<double>(		// Laplace correction
				m_drValue,									// geodetic azimuth
				zenithDistance,								// zenith distance
				stn1_it->get()->GetverticalDef(),			// deflection in prime vertical
				stn1_it->get()->GetmeridianDef(),			// deflection in prime meridian
				stn1_it->get()->GetcurrentLatitude());

			// 3. apply deflection correction
			m_drValue += m_preAdjCorr;
		}
		m_epoch = "01.10.1985";
		break;
	
	// vertical angle
	case 'Z':
		m_fInstHeight = float(1.650);
		m_fTargHeight = float(1.651);

		// compute vertical angle 
		m_drValue = VerticalAngle<double>(
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
			m_fInstHeight,				
			m_fTargHeight);

		// Deflections available?
		if (fabs(stn1_it->get()->GetverticalDef()) > E4_SEC_DEFLECTION || fabs(stn1_it->get()->GetmeridianDef()) > E4_SEC_DEFLECTION)
		{
			// 1. compute bearing from estimated coordinates
			// see above
			
			// 2. Compute correction
			m_preAdjCorr = ZenithDeflectionCorrection<double>(	// Correction to vertical angle for deflection of vertical
				direction,										// geodetic azimuth
				stn1_it->get()->GetverticalDef(),				// deflection in prime vertical
				stn1_it->get()->GetmeridianDef());				// deflection in prime meridian

			// 3. apply deflection correction (in reverse 
			// to get the simulated observation)
			m_drValue += m_preAdjCorr;
		}
		m_epoch = "01.10.1985";
		break;
	
	// zenith distance
	case 'V':
		m_fInstHeight = float(1.650);
		m_fTargHeight = float(1.651);

		// compute zenith distance 
		m_drValue = ZenithDistance<double>(
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
			m_fInstHeight,
			m_fTargHeight);

		// Deflections available?
		if (fabs(stn1_it->get()->GetverticalDef()) > E4_SEC_DEFLECTION || fabs(stn1_it->get()->GetmeridianDef()) > E4_SEC_DEFLECTION)
		{
			// 1. compute bearing from estimated coordinates
			// see above

			// 2. Compute correction
			m_preAdjCorr = ZenithDeflectionCorrection<double>(	// Correction to zenith distance for deflection of vertical
				direction,										// geodetic azimuth
				stn1_it->get()->GetverticalDef(),				// deflection in prime vertical
				stn1_it->get()->GetmeridianDef());				// deflection in prime meridian

			// 3. apply deflection correction (in reverse 
			// to get the simulated observation)
			m_drValue -= m_preAdjCorr;
		}
	}

	m_dStdDev = SecondstoRadians(0.010);

	m_epoch = "01.10.1985";
}
	

//UINT32 CDnaDirection::SetMeasurementRec(std::ifstream* ifs_stns, std::ifstream* ifs_msrs, measurement_t* measRecord)
//{
//	char stationName[STN_NAME_WIDTH];
//
//	if (GetTypeC() == 'D')		// part of a DnaDirectionSet
//	{
//		// get data relating to each direction
//		if (ifs_msrs->eof() || !ifs_msrs->good())
//			throw XMLInteropException("SetMeasurementRec(): Errors were encountered when reading from the binary measurement file.", 0);
//		ifs_msrs->read(reinterpret_cast<char *>(measRecord), sizeof(measurement_t));
//	}
//
//	m_bIgnore = measRecord->ignore;
//	m_lsetID = measRecord->clusterID;
//	m_MSmeasurementStations = (MEASUREMENT_STATIONS)measRecord->measurementStations;
//	// measRecord holds the full number of measurement blocks, which is 
//	// the numberofdirections in the vector plus one for the RO
//	m_lRecordedTotal = measRecord->vectorCount1;
//
//	if (GetTypeC() != 'D')		// not part of a DnaDirectionSet (e.g. vertical angle, zenith distance)
//	{
//		m_strType = measRecord->measType;
//		
//		// first station
//		m_lstn1Index = measRecord->station1;
//		ifs_stns->seekg(sizeof(UINT32) + measRecord->station1 * sizeof(station_t), std::ios::beg);
//		ifs_stns->read(reinterpret_cast<char *>(&stationName), sizeof(stationName));
//		m_strFirst = stationName;
//
//		m_fInstHeight = static_cast<float> (measRecord->term3);
//		m_fTargHeight = static_cast<float> (measRecord->term4);
//	}
//
//	// target station
//	m_lstn2Index = measRecord->station2;
//	ifs_stns->seekg(sizeof(UINT32) + measRecord->station2 * sizeof(station_t), std::ios::beg);
//	ifs_stns->read(reinterpret_cast<char *>(&stationName), sizeof(stationName));
//	m_strTarget = stationName;
//	
//	m_measAdj = measRecord->measAdj;
//	m_measCorr = measRecord->measCorr;
//	m_measAdjPrec = measRecord->measAdjPrec;
//	m_residualPrec = measRecord->residualPrec;
//	m_preAdjCorr = measRecord->preAdjCorr;
//	m_drValue = measRecord->term1;
//	m_dStdDev = sqrt(measRecord->term2);
//	
//	return 1;
//}
	

UINT32 CDnaDirection::SetMeasurementRec(const vstn_t& binaryStn, it_vmsr_t& it_msr, it_vdbid_t& dbidmap)
{
	m_bIgnore = it_msr->ignore;
	m_lsetID = it_msr->clusterID;
	m_MSmeasurementStations = (MEASUREMENT_STATIONS)it_msr->measurementStations;
	// it_msr holds the full number of measurement blocks, which is 
	// the numberofdirections in the vector plus one for the RO
	m_lRecordedTotal = it_msr->vectorCount1;

	if (it_msr->measType != 'D')		// not part of a DnaDirectionSet (e.g. vertical angle, zenith distance)
	{
		m_strType = it_msr->measType;
		
		// first station
		m_lstn1Index = it_msr->station1;
		m_strFirst = binaryStn.at(it_msr->station1).stationName;
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
	m_drValue = it_msr->term1;
	m_dStdDev = sqrt(it_msr->term2);

	m_epoch = it_msr->epoch;

	CDnaMeasurement::SetDatabaseMap(*dbidmap);
	
	return 0;
}
	

void CDnaDirection::WriteBinaryMsr(std::ofstream* binary_stream, PUINT32 msrIndex) const
{
	measurement_t measRecord;
	measRecord.measType = GetTypeC();
	if (m_strType == "D")
		measRecord.measStart = yMeas;
	else
		measRecord.measStart = xMeas;
	measRecord.ignore = m_bIgnore;
	measRecord.station1 = m_lstn1Index;
	measRecord.station2 = m_lstn2Index;
	measRecord.measAdj = m_measAdj;
	measRecord.measCorr = m_measCorr;
	measRecord.measAdjPrec = m_measAdjPrec;
	measRecord.residualPrec = m_residualPrec;
	measRecord.preAdjCorr = m_preAdjCorr;
	measRecord.term1 = m_drValue;
	measRecord.term2 = m_dStdDev * m_dStdDev;	// convert to variance
	measRecord.term3 = m_fInstHeight;
	measRecord.term4 = m_fTargHeight;

	measRecord.clusterID = m_lsetID;
	// number of Directions in the parent cluster including the first
	measRecord.vectorCount1 = m_lRecordedTotal;
	measRecord.measurementStations = m_MSmeasurementStations;
	measRecord.fileOrder = ((*msrIndex)++);

	snprintf(measRecord.epoch, sizeof(measRecord.epoch), "%s", m_epoch.substr(0, STN_EPOCH_WIDTH).c_str());

	binary_stream->write(reinterpret_cast<char *>(&measRecord), sizeof(measurement_t));
}


//void CDnaDirection::coutDirectionData(std::ostream &os) const
//{
//	os << std::setw(16) << " " << std::setw(MSR2_WIDTH) << m_strTarget;
//	os << std::setw(3) << (m_bIgnore ? "*" : " ") << std::setw(MEAS_WIDTH) << m_drValue;
//	os << std::setw(VAR_WIDTH) << m_dStdDev;
//	os << std::setw(7) << " ";	// m_fInstHeight
//	switch (GetTypeC())
//	{
//	case 'Z':
//	case 'V':
//		os << std::setw(7) << std::setprecision(3) << std::fixed << m_fTargHeight;
//		break;
//	default:
//		os << std::setw(7) << " ";
//		break;
//	}
//	os << std::endl;
//}

void CDnaDirection::SetValue(const std::string& str)
{
	RadFromDmsString(&m_drValue, trimstr(str));
}

void CDnaDirection::SetStdDev(const std::string& str)
{
	// convert seconds to radians
	RadFromSecondsString(&m_dStdDev, trimstr(str));
}

// Instrument and target heights only make sense for 
// slope distances, vertical angles and zenith distances.
// Note: vertical angles and zenith distances are derived
// from CDnaDirection, so these methods are provided here,
// but don't need instrument height for 'D' measurements.
void CDnaDirection::SetInstrumentHeight(const std::string& str)
{
	FloatFromString<float>(m_fInstHeight, trimstr(str));
}

void CDnaDirection::SetTargetHeight(const std::string& str)
{
	FloatFromString<float>(m_fTargHeight, trimstr(str));
}












CDnaAzimuth::CDnaAzimuth(void)
{
	m_strType = "";		// could be B or K
	m_strTarget = "";
	m_drValue = 0.;
	m_dStdDev = 0.;
	m_fInstHeight = 0.;
	m_fTargHeight = 0.;
	m_lRecordedTotal = 0;
	m_lsetID = 0;
	m_MSmeasurementStations = TWO_STATION;

}

CDnaAzimuth::~CDnaAzimuth(void)
{
}

// copy constructors
//CDnaAzimuth::CDnaAzimuth(const CDnaAzimuth& newAzimuth)
//{
//	m_strType = newAzimuth.m_strType;
//	m_strFirst = newAzimuth.m_strFirst;
//	m_bIgnore = newAzimuth.m_bIgnore;
//	m_strTarget = newAzimuth.m_strTarget;
//	m_drValue = newAzimuth.m_drValue;
//	m_dStdDev = newAzimuth.m_dStdDev;
//	m_MSmeasurementStations = newAzimuth.m_MSmeasurementStations;
//
//}


//CDnaAzimuth::CDnaAzimuth(const bool bIgnore, const std::string& strFirst, const std::string& strTarget, const double& drValue, const double& dStdDev)
//{
//	m_bIgnore = bIgnore;
//	m_strTarget = strTarget;
//	m_drValue = drValue;
//	m_dStdDev = dStdDev;
//}


//CDnaAzimuth& CDnaAzimuth::operator= (const CDnaAzimuth& rhs)
//{
//	// check for assignment to self!
//	if (this == &rhs)
//		return *this;
//
//	CDnaMeasurement::operator=(rhs);
//	m_strTarget = rhs.m_strTarget;
//	m_drValue = rhs.m_drValue;
//	m_dStdDev = rhs.m_dStdDev;
//	m_lRecordedTotal = rhs.m_lRecordedTotal;
//	m_lsetID = rhs.m_lsetID;
//	m_MSmeasurementStations = rhs.m_MSmeasurementStations;
//
//	return *this;
//}


bool CDnaAzimuth::operator== (const CDnaAzimuth& rhs) const
{
	return (
		m_strFirst == rhs.m_strFirst &&
		m_bIgnore == rhs.m_bIgnore &&
		m_strTarget == rhs.m_strTarget &&
		m_drValue == rhs.m_drValue &&
		m_dStdDev == rhs.m_dStdDev &&
		m_strType == rhs.m_strType &&
		m_epoch == rhs.m_epoch
		);
}

}	// namespace measurements
}	// namespace dynadjust
