//============================================================================
// Name         : dnagpspoint.hpp
// Author       : Roger Fraser
// Contributors : Dale Roberts <dale.o.roberts@gmail.com>
// Copyright    : Copyright 2017-2025 Geoscience Australia
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
// Description  : Interface for the CDnaGpsPoint and CDnaGpsPointCluster classes
//============================================================================

#ifndef DNAGPSPOINT_H_
#define DNAGPSPOINT_H_

#if defined(_MSC_VER)
	#if defined(LIST_INCLUDES_ON_BUILD) 
		#pragma message("  " __FILE__) 
	#endif
#endif

#include <include/measurement_types/dnameasurement.hpp>
#include <include/config/dnatypes-fwd.hpp>

namespace dynadjust {
namespace measurements {

class CDnaGpsPoint : public CDnaMeasurement
{
public:
	CDnaGpsPoint(void);
	virtual ~CDnaGpsPoint(void);

	// move constructor and move assignment operator
	CDnaGpsPoint(CDnaGpsPoint&& g);
	CDnaGpsPoint& operator=(CDnaGpsPoint&& rhs);
	
private:
	// disallowed in CDnaMeasurement
	//CDnaGpsPoint(const CDnaGpsPoint&);
	//CDnaGpsPoint& operator=(const CDnaGpsPoint& rhs);

public:
	//CDnaGpsPoint(const bool bIgnore, const std::string& strType, const std::string& strFirstStation);

	//virtual inline CDnaGpsPoint* clone() const { return new CDnaGpsPoint(*this); }
	bool operator==(const CDnaGpsPoint& rhs) const;
	bool operator<(const CDnaGpsPoint& rhs) const;

	//inline CDnaGpsPoint& operator[](int iIndex) { return this[iIndex]; }

	void AddPointCovariance(const CDnaCovariance* pGpsCovariance);

	inline UINT32 GetClusterID() const { return m_lclusterID; }
	inline std::string GetCoordType() const { return m_strCoordType; }
	
	inline std::vector<CDnaCovariance>* GetCovariances_ptr() { return &m_vPointCovariances; }

	inline void SetClusterID(const UINT32& id) { m_lclusterID = id; }
	void SetX(const std::string& str);
	void SetY(const std::string& str);
	void SetZ(const std::string& str);
	void SetSigmaXX(const std::string& str);
	void SetSigmaXY(const std::string& str);
	void SetSigmaXZ(const std::string& str);
	void SetSigmaYY(const std::string& str);
	void SetSigmaYZ(const std::string& str);
	void SetSigmaZZ(const std::string& str);
	
	void SetReferenceFrame(const std::string& refFrame);
	//void SetEpoch(const std::string& epoch);
	void SetEpsg(const std::string& epsg);	
	inline std::string GetReferenceFrame() const { return m_referenceFrame; }
	//inline std::string GetEpoch() const { return m_epoch; }
		
	//void SetPscale(const std::string& str);
	//void SetLscale(const std::string& str);
	//void SetHscale(const std::string& str);
	//void SetVscale(const std::string& str);

	inline void SetTotal(const UINT32& l) { m_lRecordedTotal = l; }

	inline void SetSigmaXX(const double& dbl) { m_dSigmaXX = dbl; }
	inline void SetSigmaXY(const double& dbl) { m_dSigmaXY = dbl; }
	inline void SetSigmaXZ(const double& dbl) { m_dSigmaXZ = dbl; }
	inline void SetSigmaYY(const double& dbl) { m_dSigmaYY = dbl; }
	inline void SetSigmaYZ(const double& dbl) { m_dSigmaYZ = dbl; }
	inline void SetSigmaZZ(const double& dbl) { m_dSigmaZZ = dbl; }

	inline void SetXAxis(const double& dbl) { m_dX = dbl; }
	inline void SetYAxis(const double& dbl) { m_dY = dbl; }
	inline void SetZAxis(const double& dbl) { m_dZ = dbl; }
	
	inline void SetPscale(const double& dbl) { m_dPscale = dbl; }
	inline void SetLscale(const double& dbl) { m_dLscale = dbl; }
	inline void SetHscale(const double& dbl) { m_dHscale = dbl; }
	inline void SetVscale(const double& dbl) { m_dVscale = dbl; }
	void SetCoordType(const std::string& str);

	inline void SetRecordedTotal(const UINT32& total) { m_lRecordedTotal = total; }

	_COORD_TYPE_ GetMyCoordTypeC();

	void ReserveGpsCovariancesCount(const UINT32& size);
	void ResizeGpsCovariancesCount(const UINT32& size = 0);

	virtual UINT32 CalcBinaryRecordCount() const;
	//void coutPointData(std::ostream &os) const;
	virtual void WriteBinaryMsr(std::ofstream* binary_stream, PUINT32 msrIndex) const;
	virtual UINT32 SetMeasurementRec(const vstn_t& binaryStn, it_vmsr_t& it_msr, it_vdbid_t& dbidmap);
	virtual void WriteDynaMLMsr(std::ofstream* dynaml_stream, const std::string& comment, bool) const;
	virtual void WriteDNAMsr(std::ofstream* dna_stream, const dna_msr_fields& dmw, const dna_msr_fields& dml, bool) const;
	virtual void SimulateMsr(vdnaStnPtr* vStations, const CDnaEllipsoid* ellipsoid);
	virtual void PopulateMsr(pvstn_t bstRecords, uint32_uint32_map* blockStationsMap, vUINT32* blockStations,
		const UINT32& stn, const CDnaDatum* datum, math::matrix_2d* estimates, math::matrix_2d* variances);

	virtual void SerialiseDatabaseMap(std::ofstream* os);

	inline double GetVscale() const { return m_dVscale; }
	inline double GetPscale() const { return m_dPscale; }
	inline double GetLscale() const { return m_dLscale; }
	inline double GetHscale() const { return m_dHscale; }

	inline double GetValue() const { return sqrt(m_dX*m_dX + m_dY*m_dY + m_dZ*m_dZ); }			// Virtual magnitude
	inline double GetStdDev() const { return sqrt(m_dSigmaXX + m_dSigmaYY + m_dSigmaZZ); }		// RMS

	//virtual inline void SetDatabaseMap_bmsIndex(const UINT32& bmsIndex) { m_msr_db_map.bms_index = bmsIndex; }
protected:

	UINT32 m_lRecordedTotal;

	double m_dX;
	double m_dY;
	double m_dZ;
	double m_dSigmaXX;
	double m_dSigmaXY;
	double m_dSigmaXZ;
	double m_dSigmaYY;
	double m_dSigmaYZ;
	double m_dSigmaZZ;

	double m_dPscale;
	double m_dLscale;
	double m_dHscale;
	double m_dVscale;
	std::string m_strCoordType;

	COORD_TYPE m_ctType;

	std::string	m_referenceFrame;
	//string	m_epoch;

	UINT32 m_lclusterID;

	std::vector<CDnaCovariance> m_vPointCovariances;
};


class CDnaGpsPointCluster : public CDnaMeasurement
{
public:
	CDnaGpsPointCluster(void);
	virtual ~CDnaGpsPointCluster(void);

	// move constructor and move assignment operator
	CDnaGpsPointCluster(CDnaGpsPointCluster&& p);
	CDnaGpsPointCluster& operator=(CDnaGpsPointCluster&& rhs);

private:
	// disallow copying
	CDnaGpsPointCluster(const CDnaGpsPointCluster&);
	CDnaGpsPointCluster& operator=(const CDnaGpsPointCluster& rhs);

public:
	//CDnaGpsPointCluster(const bool bIgnore, const std::string& strType, const std::string& strFirstStation);
	CDnaGpsPointCluster(const UINT32 lclusterID, const std::string& referenceframe, const std::string& epoch);

	//virtual inline CDnaGpsPointCluster* clone() const { return new CDnaGpsPointCluster(*this); }
	bool operator==(const CDnaGpsPointCluster& rhs) const;
	bool operator<(const CDnaGpsPointCluster& rhs) const;

	//inline CDnaGpsPointCluster& operator[](int iIndex) { return this[iIndex]; }

	//inline UINT32 GetNumPoints() { return m_vGpsPoints.size(); }
	inline std::vector<CDnaGpsPoint>& GetPoints() { return m_vGpsPoints; }
	inline std::vector<CDnaGpsPoint>* GetPoints_ptr() override { return &m_vGpsPoints; }

	inline UINT32 GetClusterID() const override { return m_lclusterID; }
	inline std::string GetCoordType() const override { return m_strCoordType; }
	inline UINT32 GetTotal() const override { return m_lRecordedTotal; }
	inline double GetPscale() const override { return m_dPscale; }
	inline double GetLscale() const override { return m_dLscale; }
	inline double GetHscale() const override { return m_dHscale; }
	inline double GetVscale() const override { return m_dVscale; }

	void SetCoordType(const std::string& str) override;

	_COORD_TYPE_ GetMyCoordTypeC();

	void SetReferenceFrame(const std::string& refFrame) override;
	//void SetEpoch(const std::string& epoch);
	void SetEpsg(const std::string& epsg);
	inline std::string GetReferenceFrame() const override { return m_referenceFrame; }
	//inline std::string GetEpoch() const { return m_epoch; }

	//inline void SetPoints(const std::vector<CDnaGpsPoint>& d) { m_vGpsPoints = d; }
	void SetTotal(const std::string& str) override;
	void SetPscale(const std::string& str) override;
	void SetLscale(const std::string& str) override;
	void SetHscale(const std::string& str) override;
	void SetVscale(const std::string& str) override;

	inline void SetTotal(const UINT32& u) override { m_lRecordedTotal = u; }

	inline void SetPscale(const double& dbl) override { m_dPscale = dbl; }
	inline void SetLscale(const double& dbl) override { m_dLscale = dbl; }
	inline void SetHscale(const double& dbl) override { m_dHscale = dbl; }
	inline void SetVscale(const double& dbl) override { m_dVscale = dbl; }

	void AddGpsPoint(const CDnaMeasurement* pGpsPoint) override;
	//void ClearPoints();

	void SetSourceFileIndex(const UINT32& idx) override;

	UINT32 CalcBinaryRecordCount() const override;
	void WriteBinaryMsr(std::ofstream* binary_stream, PUINT32 msrIndex) const override;
	UINT32 SetMeasurementRec(const vstn_t& binaryStn, it_vmsr_t& it_msr, it_vdbid_t& dbidmap) override;
	void WriteDynaMLMsr(std::ofstream* dynaml_stream, const std::string& comment, bool) const override;
	void WriteDNAMsr(std::ofstream* dna_stream, const dna_msr_fields& dmw, const dna_msr_fields& dml, bool) const override;
	void SimulateMsr(vdnaStnPtr* vStations, const CDnaEllipsoid* ellipsoid) override;
	void PopulateMsr(pvstn_t bstRecords, uint32_uint32_map* blockStationsMap, vUINT32* blockStations,
		const UINT32& block, const CDnaDatum* datum, math::matrix_2d* estimates, math::matrix_2d* variances) override;

	void SerialiseDatabaseMap(std::ofstream* os) override;

	void SetDatabaseMaps(it_vdbid_t& dbidmap) override;

protected:

	UINT32 m_lRecordedTotal;
	double m_dPscale;
	double m_dLscale;
	double m_dHscale;
	double m_dVscale;
	std::string m_strCoordType;

	COORD_TYPE m_ctType;

	std::vector<CDnaGpsPoint> m_vGpsPoints;

	std::string	m_referenceFrame;
	//string	m_epoch;

	UINT32 m_lclusterID;

	it_vdbid_t m_dbidmap;
};

}	// namespace measurements
}	// namespace dynadjust

#endif /* DNAGPSPOINT_H_ */
