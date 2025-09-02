//============================================================================
// Name         : dnadirectionset.hpp
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
// Description  : Interface for the CDnaDirectionSet class
//============================================================================

#ifndef DNADIRECTIONSET_H_
#define DNADIRECTIONSET_H_

#if defined(_MSC_VER)
	#if defined(LIST_INCLUDES_ON_BUILD) 
		#pragma message("  " __FILE__) 
	#endif
#endif

/// \cond
#include <memory>

#include <boost/shared_ptr.hpp>
/// \endcond

#include <include/measurement_types/dnastation.hpp>
#include <include/measurement_types/dnadirection.hpp>
#include <include/measurement_types/dnameasurement.hpp>


namespace dynadjust {
namespace measurements {

class CDnaDirectionSet : public CDnaMeasurement
{
public:
	CDnaDirectionSet(void);
	virtual ~CDnaDirectionSet(void);

	// move constructor and move assignment operator
	CDnaDirectionSet(CDnaDirectionSet&& d);
	CDnaDirectionSet& operator=(CDnaDirectionSet&& rhs);

private:
	// disallow copying
	//CDnaDirectionSet(const CDnaDirectionSet&);
	//CDnaDirectionSet& operator=(const CDnaDirectionSet& rhs);

public:
	CDnaDirectionSet(const UINT32 lsetID);

	//CDnaDirectionSet(bool bIgnore,
	//		const std::string& strFirst, const std::string& strTarget,
	//		const double& drValue, const double& dStdDev,
	//		const float& fInstrHeight, const float& fTargetHeight);

	bool operator==(const CDnaDirectionSet& rhs) const;
	bool operator<(const CDnaDirectionSet& rhs) const;

	//inline CDnaDirectionSet& operator[](int iIndex) { return this[iIndex]; }

	void LoadDirectionSet(const char* const, const int&, const std::string&, const std::string&, bool, const int&);

	inline UINT32 GetClusterID() const { return m_lsetID; }
	inline std::string GetTarget() const { return m_strTarget; }
	inline UINT32 GetTotal() const { return m_lRecordedTotal; }
	inline double GetValue() const { return m_drValue; }
	inline double GetStdDev() const { return m_dStdDev; }
	
	inline size_t GetNumDirections() const { return m_vTargetDirections.size(); }
	inline std::vector<CDnaDirection>* GetDirections_ptr() { return &m_vTargetDirections; }

	inline void SetClusterID(const UINT32& id) { m_lsetID = id; }
	inline void SetTarget(const std::string& str) { m_strTarget = trimstr(str); }
	inline void SetTotal(const UINT32& l) { m_lRecordedTotal = l; }
	inline void SetNonIgnoredDirns(const UINT32& n) { m_lNonIgnoredDirns = n; }

	void SetTotal(const std::string& str);
	void SetValue(const std::string& str);
	void SetStdDev(const std::string& str);
	
	void AddDirection(const CDnaMeasurement* pDirection);
	void ClearDirections();
	//bool IsRepeatedDirection(string);

	virtual UINT32 CalcBinaryRecordCount() const;
	virtual void WriteBinaryMsr(std::ofstream* binary_stream, PUINT32 msrIndex) const;
	virtual UINT32 SetMeasurementRec(const vstn_t& binaryStn, it_vmsr_t& it_msr, it_vdbid_t& dbidmap);
	virtual void WriteDynaMLMsr(std::ofstream* dynaml_stream, const std::string& comment, bool) const;
	virtual void WriteDNAMsr(std::ofstream* dna_stream, const dna_msr_fields& dmw, const dna_msr_fields& dml, bool) const;
	virtual void SimulateMsr(vdnaStnPtr* vStations, const CDnaEllipsoid* ellipsoid);

	virtual void SerialiseDatabaseMap(std::ofstream* os);

	void SetDatabaseMaps(it_vdbid_t& dbidmap);

	std::string m_strTarget;

protected:
	double m_drValue;
	double m_dStdDev;
	UINT32 m_lRecordedTotal;
	UINT32 m_lNonIgnoredDirns;
	std::vector<CDnaDirection> m_vTargetDirections;
	UINT32 m_lsetID;

	it_vdbid_t m_dbidmap;
};

}	// namespace measurements
}	// namespace dynadjust

#endif /* DNADIRECTIONSET_H_ */
