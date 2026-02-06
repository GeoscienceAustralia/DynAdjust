//============================================================================
// Name         : dnadirectionset.cpp
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
// Description  : CDnaDirectionSet implementation file
//============================================================================

#include <include/functions/dnatemplatefuncs.hpp>
#include <include/measurement_types/dnadirection.hpp>
#include <include/measurement_types/dnadirectionset.hpp>

namespace dynadjust {
namespace measurements {

CDnaDirectionSet::CDnaDirectionSet(void)
	: m_strTarget("")
	, m_drValue(0.)
	, m_dStdDev(0.)
	, m_lRecordedTotal(0)
	, m_lNonIgnoredDirns(0)
	, m_lsetID(0)
{
	m_strType = "D";
	m_MSmeasurementStations = TWO_STATION;

	m_vTargetDirections.clear();
	m_vTargetDirections.reserve(15);
}

CDnaDirectionSet::~CDnaDirectionSet(void)
{
	m_vTargetDirections.clear();
}
	
CDnaDirectionSet::CDnaDirectionSet(const UINT32 lsetID)
	: m_strTarget("")
	, m_drValue(0.)
	, m_dStdDev(0.)
	, m_lRecordedTotal(0)
	, m_lNonIgnoredDirns(0)
	, m_lsetID(lsetID)
{
	m_strType = "D";
	m_MSmeasurementStations = TWO_STATION;

	m_vTargetDirections.clear();
	m_vTargetDirections.reserve(15);
}

// move constructors
CDnaDirectionSet::CDnaDirectionSet(CDnaDirectionSet&& d)
{
	m_bIgnore = d.m_bIgnore;
	m_bInsufficient = d.m_bInsufficient;
	m_strFirst = d.m_strFirst;
	m_strTarget = d.m_strTarget;
	m_drValue = d.m_drValue;
	m_dStdDev = d.m_dStdDev;
	m_lRecordedTotal = d.m_lRecordedTotal;
	m_lNonIgnoredDirns = d.m_lNonIgnoredDirns;
	m_MSmeasurementStations = d.m_MSmeasurementStations;

	m_lsetID = d.m_lsetID;

	m_strType = "D";
	m_vTargetDirections = std::move(d.m_vTargetDirections);

	m_msr_db_map = d.m_msr_db_map;

	m_dbidmap = d.m_dbidmap;
}

// move assignment operator 
CDnaDirectionSet& CDnaDirectionSet::operator= (CDnaDirectionSet&& rhs)
{
	// check for assignment to self!
	if (this == &rhs)
		return *this;

	CDnaMeasurement::operator=(std::move(rhs));
	m_strTarget = rhs.m_strTarget;
	m_drValue = rhs.m_drValue;
	m_dStdDev = rhs.m_dStdDev;
	m_lRecordedTotal = rhs.m_lRecordedTotal;
	m_lNonIgnoredDirns = rhs.m_lNonIgnoredDirns;
	
	m_lsetID = rhs.m_lsetID;
	m_MSmeasurementStations = rhs.m_MSmeasurementStations;
	
	m_vTargetDirections = std::move(rhs.m_vTargetDirections);
	
	m_msr_db_map = rhs.m_msr_db_map;

	m_dbidmap = rhs.m_dbidmap;

	return *this;
}
	

bool CDnaDirectionSet::operator== (const CDnaDirectionSet& rhs) const
{
	return (
		m_strFirst == rhs.m_strFirst &&
		m_strType == rhs.m_strType &&
		m_bIgnore == rhs.m_bIgnore &&
		m_strTarget == rhs.m_strTarget &&
		m_drValue == rhs.m_drValue &&
		m_dStdDev == rhs.m_dStdDev &&
		m_lRecordedTotal == rhs.m_lRecordedTotal &&
		m_lNonIgnoredDirns == rhs.m_lNonIgnoredDirns &&
		m_vTargetDirections == rhs.m_vTargetDirections &&
		m_epoch == rhs.m_epoch
		);
}


bool CDnaDirectionSet::operator< (const CDnaDirectionSet& rhs) const
{
	if (m_strFirst == rhs.m_strFirst) {
		if (m_strType == rhs.m_strType) {	// I don't think this is needed
			if (m_bIgnore == rhs.m_bIgnore) {
				if (m_epoch == rhs.m_epoch) {
					if (m_strTarget == rhs.m_strTarget) {
						if (m_drValue == rhs.m_drValue) {
							if (m_dStdDev == rhs.m_dStdDev)
								return m_vTargetDirections < rhs.m_vTargetDirections;
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


void CDnaDirectionSet::AddDirection(const CDnaMeasurement* pDirection)
{
	m_vTargetDirections.push_back(std::move((CDnaDirection&)*pDirection));
}


void CDnaDirectionSet::ClearDirections()
{
	m_vTargetDirections.clear();
}


//bool CDnaDirectionSet::IsRepeatedDirection(std::string strTarget)
//{
//	if (m_vTargetDirections.empty())
//		return false;
//
//	if (find_if(m_vTargetDirections.begin(), m_vTargetDirections.end(),
//		bind2nd(operator_on_mem(&CDnaDirection::m_strTarget, std::equal_to<std::string>()), strTarget)) != m_vTargetDirections.end())
//		return true;
//	return false;
//}
	

void CDnaDirectionSet::SerialiseDatabaseMap(std::ofstream* os)
{
	CDnaMeasurement::SerialiseDatabaseMap(os);

	for_each(m_vTargetDirections.begin(), m_vTargetDirections.end(),
		[this, os](const CDnaDirection& dir) {
			((CDnaMeasurement*)&dir)->SerialiseDatabaseMap(os);
	});
}
	

void CDnaDirectionSet::SetDatabaseMaps(it_vdbid_t& dbidmap)
{
	m_dbidmap = dbidmap;

	CDnaMeasurement::SetDatabaseMap(*m_dbidmap);

	for_each(m_vTargetDirections.begin(), m_vTargetDirections.end(),
		[this](const CDnaDirection& dir) {
			((CDnaMeasurement*)&dir)->SetDatabaseMap(*m_dbidmap);
	m_dbidmap++;
		});
	
}


UINT32 CDnaDirectionSet::CalcBinaryRecordCount() const
{
	// Direction set has 1 RO direction and n directions
	UINT32 recordCount(1);
	for_each(m_vTargetDirections.begin(), m_vTargetDirections.end(),
		[&recordCount](const CDnaDirection& dir) {
			recordCount += dir.CalcBinaryRecordCount();
	});
	return recordCount;
}
	

void CDnaDirectionSet::WriteDynaMLMsr(std::ofstream* dynaml_stream, const std::string& comment, bool) const
{
	const size_t dirCount(GetNumDirections());

	if (comment.empty())
	{
		*dynaml_stream << "  <!-- Type " << measurement_name<char, std::string>(GetTypeC());
		if (dirCount > 1)
			*dynaml_stream << " (set of " << dirCount << ")";
		else
			*dynaml_stream << " (single)";
		*dynaml_stream << " -->" << std::endl;
	}
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
	*dynaml_stream << "    <Value>" << std::setprecision(8) << std::fixed << RadtoDms(m_drValue) << "</Value>" << std::endl;
	*dynaml_stream << "    <StdDev>" << std::scientific << std::setprecision(6) << Seconds(m_dStdDev) << "</StdDev>" << std::endl;
	
	if (m_msr_db_map.is_msr_id_set)
		*dynaml_stream << "    <MeasurementID>" << m_msr_db_map.msr_id << "</MeasurementID>" << std::endl;
	if (m_msr_db_map.is_cls_id_set)
		*dynaml_stream << "    <ClusterID>" << m_msr_db_map.cluster_id << "</ClusterID>" << std::endl;
	
	*dynaml_stream << "    <Total>" << dirCount << "</Total>" << std::endl;
	
	// write directions
	std::vector<CDnaDirection>::const_iterator _it_dir = m_vTargetDirections.begin();
	for (_it_dir = m_vTargetDirections.begin(); _it_dir!=m_vTargetDirections.end(); ++_it_dir)
		_it_dir->WriteDynaMLMsr(dynaml_stream, comment, true);
	
	*dynaml_stream << "  </DnaMeasurement>" << std::endl;
}
	

void CDnaDirectionSet::WriteDNAMsr(std::ofstream* dna_stream, const dna_msr_fields& dmw, const dna_msr_fields& dml, bool) const
{
	const size_t dirCount(GetNumDirections());
	
	*dna_stream << std::setw(dmw.msr_type) << m_strType;
	if (m_bIgnore)
		*dna_stream << std::setw(dmw.msr_ignore) << "*";
	else
		*dna_stream << std::setw(dmw.msr_ignore) << " ";
	
	*dna_stream << std::left << std::setw(dmw.msr_inst) << m_strFirst;
	*dna_stream << std::left << std::setw(dmw.msr_targ1) << m_strTarget;
	*dna_stream << std::left << std::setw(dmw.msr_targ2) << dirCount;
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
	if (m_msr_db_map.is_cls_id_set)
		*dna_stream << std::setw(dmw.msr_id_cluster) << m_msr_db_map.cluster_id;
	
	*dna_stream << std::endl;
	
	// write directions
	std::vector<CDnaDirection>::const_iterator _it_dir = m_vTargetDirections.begin();
	for (_it_dir = m_vTargetDirections.begin(); _it_dir!=m_vTargetDirections.end(); ++_it_dir)
		_it_dir->WriteDNAMsr(dna_stream, dmw, dml, true);
}
	

void CDnaDirectionSet::SimulateMsr(vdnaStnPtr* vStations, const CDnaEllipsoid* ellipsoid)
{
	_it_vdnastnptr stn1_it(vStations->begin() + m_lstn1Index);
	_it_vdnastnptr stn2_it(vStations->begin() + m_lstn2Index);
	
	// compute direction from instument to first target
	m_drValue = Direction(
		stn1_it->get()->GetXAxis(), 
		stn1_it->get()->GetYAxis(),
		stn1_it->get()->GetZAxis(),
		vStations->at(m_lstn2Index).get()->GetXAxis(), 
		vStations->at(m_lstn2Index).get()->GetYAxis(),
		vStations->at(m_lstn2Index).get()->GetZAxis(),
		stn1_it->get()->GetcurrentLatitude(),
		stn1_it->get()->GetcurrentLongitude());

	// Deflections available?
	if (fabs(stn1_it->get()->GetverticalDef()) > E4_SEC_DEFLECTION || fabs(stn1_it->get()->GetmeridianDef()) > E4_SEC_DEFLECTION)
	{
		// 1. compute zenith distance 
		double zenithDistance = ZenithDistance<double>(
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

	// compute standard deviation in radians
	m_dStdDev = SecondstoRadians(0.01);

	m_epoch = "01.10.1985";
	
	std::vector<CDnaDirection>::iterator _it_dir = m_vTargetDirections.begin();
	for (_it_dir = m_vTargetDirections.begin(); _it_dir!=m_vTargetDirections.end(); ++_it_dir)
		_it_dir->SimulateMsr(vStations, ellipsoid);
}
	

UINT32 CDnaDirectionSet::SetMeasurementRec(const vstn_t& binaryStn, it_vmsr_t& it_msr, it_vdbid_t& dbidmap)
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
	
	m_measAdj = it_msr->measAdj;
	m_measCorr = it_msr->measCorr;
	m_measAdjPrec = it_msr->measAdjPrec;
	m_residualPrec = it_msr->residualPrec;
	m_preAdjCorr = it_msr->preAdjCorr;
	m_drValue = it_msr->term1;
	m_dStdDev = sqrt(it_msr->term2);
	
	// it_msr holds the full number of measurement blocks, which is 
	// the number of directions in the vector plus one for the RO
	m_lRecordedTotal = it_msr->vectorCount1 - 1;
	// Calculate afresh the number of non-ignored directions
	if (it_msr->vectorCount2 > 0)
		m_lNonIgnoredDirns = it_msr->vectorCount2 - 1;
	
	m_lsetID = it_msr->clusterID;

	m_epoch = it_msr->epoch;
	m_sourceFileIndex = it_msr->sourceFileIndex;

	m_vTargetDirections.clear();
	m_vTargetDirections.resize(m_lRecordedTotal);

	UINT32 it_msrCount = m_lRecordedTotal;	

	CDnaMeasurement::SetDatabaseMap(*dbidmap);

	// now covariances
	std::vector<CDnaDirection>::iterator _it_dir = m_vTargetDirections.begin();
	for (; _it_dir!=m_vTargetDirections.end(); ++_it_dir)
	{
		it_msr++;
		dbidmap++;
		_it_dir->SetType(m_strType);
		_it_dir->SetFirst(m_strFirst);
		_it_dir->SetMeasurementRec(binaryStn, it_msr, dbidmap);
	}

	return it_msrCount;
}
	

void CDnaDirectionSet::WriteBinaryMsr(std::ofstream* binary_stream, PUINT32 msrIndex) const
{
	measurement_t measRecord;
	measRecord.measType = GetTypeC();
	measRecord.measStart = xMeas;
	measRecord.ignore = m_bIgnore;
	measRecord.station1 = m_lstn1Index;
	measRecord.station2 = m_lstn2Index;
	
	// number of Directions in the parent cluster including the first
	measRecord.vectorCount1 = static_cast<UINT32>(m_vTargetDirections.size()) + 1;	
	// number of non-ignored directions in the cluster including the first
	if (!m_bIgnore && m_lNonIgnoredDirns > 0)
		measRecord.vectorCount2 = m_lNonIgnoredDirns + 1;

	measRecord.measAdj = m_measAdj;
	measRecord.measCorr = m_measCorr;
	measRecord.measAdjPrec = m_measAdjPrec;
	measRecord.residualPrec = m_residualPrec;
	measRecord.preAdjCorr = m_preAdjCorr;
	measRecord.term1 = m_drValue;
	measRecord.term2 = m_dStdDev * m_dStdDev;	// convert to variance

	measRecord.clusterID = m_lsetID;
	measRecord.measurementStations = m_MSmeasurementStations;
	measRecord.sourceFileIndex = m_sourceFileIndex;
	measRecord.fileOrder = ((*msrIndex)++);

	snprintf(measRecord.epoch, sizeof(measRecord.epoch), "%s", m_epoch.substr(0, STN_EPOCH_WIDTH).c_str());

	binary_stream->write(reinterpret_cast<char *>(&measRecord), sizeof(measurement_t));

	std::vector<CDnaDirection>::const_iterator _it_dir;
	for (_it_dir = m_vTargetDirections.begin(); _it_dir!=m_vTargetDirections.end(); ++_it_dir)
		_it_dir->WriteBinaryMsr(binary_stream, msrIndex);
}

void CDnaDirectionSet::SetValue(const std::string& str)
{
	RadFromDmsString(&m_drValue, trimstr(str));
}

void CDnaDirectionSet::SetStdDev(const std::string& str)
{
	RadFromSecondsString(&m_dStdDev, trimstr(str));
}

void CDnaDirectionSet::SetTotal(const std::string& str)
{
	m_lRecordedTotal = LongFromString<UINT32>(trimstr(str));
	m_vTargetDirections.reserve(m_lRecordedTotal);
}




}	// namespace measurements
}	// namespace dynadjust
