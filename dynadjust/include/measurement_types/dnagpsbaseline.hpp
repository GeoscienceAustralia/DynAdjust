//============================================================================
// Name         : dnagpsbaseline.hpp
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
// Description  : Interface for the CDnaGpsBaseline and CDnaGpsBaselineCluster 
//                classes
//============================================================================

#ifndef DNAGPSBASELINE_H_
#define DNAGPSBASELINE_H_

#if defined(_MSC_VER)
	#if defined(LIST_INCLUDES_ON_BUILD) 
		#pragma message("  " __FILE__) 
	#endif
#endif

#include <include/measurement_types/dnameasurement.hpp>

namespace dynadjust {
namespace measurements {

class CDnaGpsBaseline : public CDnaMeasurement
{
public:
	CDnaGpsBaseline(void);
	virtual ~CDnaGpsBaseline(void);

	// move constructor and move assignment operator
	CDnaGpsBaseline(CDnaGpsBaseline&& g);
	CDnaGpsBaseline& operator=(CDnaGpsBaseline&& rhs);

private:
	// disallowed in CDnaMeasurement
	//CDnaGpsBaseline(const CDnaGpsBaseline&);
	//CDnaGpsBaseline& operator=(const CDnaGpsBaseline& rhs);

public:
	//CDnaGpsBaseline(const bool bIgnore, const std::string& strType, const std::string& strFirstStation, const std::string& strSecondStation);

	//virtual inline CDnaGpsBaseline* clone() const { return new CDnaGpsBaseline(*this); }
	bool operator==(const CDnaGpsBaseline& rhs) const;
	bool operator<(const CDnaGpsBaseline& rhs) const;

	//inline CDnaGpsBaseline& operator[](int iIndex) { return this[iIndex]; }

	inline UINT32 GetClusterID() const { return m_lclusterID; }
	inline std::string GetTarget() const { return m_strTarget; }
	inline std::vector<CDnaCovariance>* GetCovariances_ptr() { return &m_vGpsCovariances; }

	void AddGpsCovariance(const CDnaCovariance* pGpsCovariance);

	void ReserveGpsCovariancesCount(const UINT32& size);
	void ResizeGpsCovariancesCount(const UINT32& size = 0);

	inline void SetClusterID(const UINT32& id) { m_lclusterID = id; }
	inline void SetTarget(const std::string& str) { m_strTarget = trimstr(str); }
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
	//void SetEpoch(const std::string& epoch); //moved to CDnaMeasurement
	void SetEpsg(const std::string& epsg);	
	inline std::string GetReferenceFrame() const { return m_referenceFrame; }
	//inline std::string GetEpoch() const { return m_epoch; } // moved to CDnaMeasurement
		
	void SetPscale(const std::string& str);
	void SetLscale(const std::string& str);
	void SetHscale(const std::string& str);
	void SetVscale(const std::string& str);

	inline void SetPscale(const double& dbl) { m_dPscale = dbl; }
	inline void SetLscale(const double& dbl) { m_dLscale = dbl; }
	inline void SetHscale(const double& dbl) { m_dHscale = dbl; }
	inline void SetVscale(const double& dbl) { m_dVscale = dbl; }

	inline void SetRecordedTotal(const UINT32& total) { m_lRecordedTotal = total; }
	inline void SetTotal(const UINT32& l) { m_lRecordedTotal = l; }

	virtual UINT32 CalcBinaryRecordCount() const;
	//void coutBaselineData(std::ostream &os, const int& pad, const UINT16& uType = 0) const;
	virtual void WriteBinaryMsr(std::ofstream* binary_stream, PUINT32 msrIndex) const;
	virtual UINT32 SetMeasurementRec(const vstn_t& binaryStn, it_vmsr_t& it_msr, it_vdbid_t& dbidmap);
	virtual void WriteDynaMLMsr(std::ofstream* dynaml_stream, const std::string& comment, bool) const;
	virtual void WriteDNAMsr(std::ofstream* dna_stream, const dna_msr_fields& dmw, const dna_msr_fields& dml, bool) const;
	virtual void SimulateMsr(vdnaStnPtr* vStations, const CDnaEllipsoid* ellipsoid);

	virtual void SerialiseDatabaseMap(std::ofstream* os);

	std::string m_strTarget;

	inline double GetVscale() const { return m_dVscale; }
	inline double GetPscale() const { return m_dPscale; }
	inline double GetLscale() const { return m_dLscale; }
	inline double GetHscale() const { return m_dHscale; }

	inline double GetValue() const { return sqrt(m_dX*m_dX + m_dY*m_dY + m_dZ*m_dZ); }			// Magnitude
	inline double GetStdDev() const { return sqrt(m_dSigmaXX + m_dSigmaYY + m_dSigmaZZ); }		// RMS

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

	std::string	m_referenceFrame;
	//string	m_epoch;

	UINT32 m_lclusterID;

	std::vector<CDnaCovariance> m_vGpsCovariances;
};

// used for Type G or X
class CDnaGpsBaselineCluster : public CDnaMeasurement
{
public:
	CDnaGpsBaselineCluster(void);
	virtual ~CDnaGpsBaselineCluster(void);

	// move constructor and move assignment operator
	CDnaGpsBaselineCluster(CDnaGpsBaselineCluster&& g);
	CDnaGpsBaselineCluster& operator=(CDnaGpsBaselineCluster&& rhs);

private:
	// disallow copying
	CDnaGpsBaselineCluster(const CDnaGpsBaselineCluster&);
	CDnaGpsBaselineCluster& operator=(const CDnaGpsBaselineCluster& rhs);

public:
	//CDnaGpsBaselineCluster(const bool bIgnore, const std::string& strType, const std::string& strFirstStation);
	CDnaGpsBaselineCluster(const UINT32 lclusterID, const std::string& referenceframe, const std::string& epoch);

	//virtual inline CDnaGpsBaselineCluster* clone() const { return new CDnaGpsBaselineCluster(*this); }
	bool operator==(const CDnaGpsBaselineCluster& rhs) const;
	bool operator<(const CDnaGpsBaselineCluster& rhs) const;

	//inline CDnaGpsBaselineCluster& operator[](int iIndex) { return this[iIndex]; }

	//inline UINT32 GetNumBaselinens() { return m_vGpsBaselines.size(); }
	inline std::vector<CDnaGpsBaseline>& GetBaselines() { return m_vGpsBaselines; }
	inline std::vector<CDnaGpsBaseline>* GetBaselines_ptr() override { return &m_vGpsBaselines; }

	inline UINT32 GetClusterID() const override { return m_lclusterID; }
	inline std::string GetTarget() const override { return m_strTarget; }
	inline UINT32 GetTotal() const override { return m_lRecordedTotal; }
	inline double GetPscale() const override { return m_dPscale; }
	inline double GetLscale() const override { return m_dLscale; }
	inline double GetHscale() const override { return m_dHscale; }
	inline double GetVscale() const override { return m_dVscale; }

	void SetReferenceFrame(const std::string& refFrame) override;
	//void SetEpoch(const std::string& epoch);
	void SetEpsg(const std::string& epsg);
	inline std::string GetReferenceFrame() const override { return m_referenceFrame; }
	//inline std::string GetEpoch() const { return m_epoch; }

	inline void SetTarget(const std::string& str) override { m_strTarget = trimstr(str); }
	inline void SetTotal(const UINT32& l) override { m_lRecordedTotal = l; }

	void SetTotal(const std::string& str) override;
	void SetPscale(const std::string& str) override;
	void SetLscale(const std::string& str) override;
	void SetHscale(const std::string& str) override;
	void SetVscale(const std::string& str) override;

	inline void SetPscale(const double& dbl) override { m_dPscale = dbl; }
	inline void SetLscale(const double& dbl) override { m_dLscale = dbl; }
	inline void SetHscale(const double& dbl) override { m_dHscale = dbl; }
	inline void SetVscale(const double& dbl) override { m_dVscale = dbl; }

	void AddGpsBaseline(const CDnaMeasurement* pGpsBaseline) override;

	void ClearBaselines();

	void ReserveGpsBaselinesCount(const UINT32& size) override;

	void SetSourceFileIndex(const UINT32& idx) override;

	UINT32 CalcBinaryRecordCount() const override;
	void WriteBinaryMsr(std::ofstream* binary_stream, PUINT32 msrIndex) const override;
	UINT32 SetMeasurementRec(const vstn_t& binaryStn, it_vmsr_t& it_msr, it_vdbid_t& dbidmap) override;
	void WriteDynaMLMsr(std::ofstream* dynaml_stream, const std::string& comment, bool) const override;
	void WriteDNAMsr(std::ofstream* dna_stream, const dna_msr_fields& dmw, const dna_msr_fields& dml, bool) const override;
	void SimulateMsr(vdnaStnPtr* vStations, const CDnaEllipsoid* ellipsoid) override;

	void SerialiseDatabaseMap(std::ofstream* os) override;

	void SetDatabaseMaps(it_vdbid_t& dbidmap) override;

protected:

	std::string m_strTarget;
	UINT32 m_lRecordedTotal;
	double m_dPscale;
	double m_dLscale;
	double m_dHscale;
	double m_dVscale;

	std::vector<CDnaGpsBaseline> m_vGpsBaselines;

	std::string	m_referenceFrame;
	//std::string	m_epoch;

	UINT32 m_lclusterID;

	it_vdbid_t m_dbidmap;

};

}	// namespace measurements
}	// namespace dynadjust

#endif /* DNAGPSBASELINE_H_ */
