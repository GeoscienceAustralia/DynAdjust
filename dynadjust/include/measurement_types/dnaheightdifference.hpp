//============================================================================
// Name         : dnaheightdifference.hpp
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
// Description  : Interface for the CDnaHeightDifference class
//============================================================================

#ifndef DNAHEIGHTDIFFERENCE_H_
#define DNAHEIGHTDIFFERENCE_H_

#if defined(_MSC_VER)
	#if defined(LIST_INCLUDES_ON_BUILD) 
		#pragma message("  " __FILE__) 
	#endif
#endif

#include <include/measurement_types/dnaheight.hpp>

namespace dynadjust {
namespace measurements {

class CDnaHeightDifference : public CDnaHeight
{
public:
	CDnaHeightDifference(void);
	virtual ~CDnaHeightDifference(void);

private:
	// disallowed in CDnaMeasurement
	//CDnaHeightDifference(const CDnaHeightDifference&) {};
	//CDnaHeightDifference& operator=(const CDnaHeightDifference& rhs);

public:
	//CDnaHeightDifference(const bool bIgnore, const std::string& strType, const std::string& strFirst, const std::string& strTarget, const double& dValue, const double& dStdDev);

	//virtual inline CDnaHeightDifference* clone() const { return new CDnaHeightDifference(*this); }
	bool operator==(const CDnaHeightDifference& rhs) const;
	bool operator<(const CDnaHeightDifference& rhs) const;

	//inline CDnaHeightDifference& operator[](int iIndex) { return this[iIndex]; }

	inline std::string GetTarget() const { return m_strTarget; }
	inline double GetValue() const { return m_dValue; }
	inline double GetStdDev() const { return m_dStdDev; }
	
	inline void SetTarget(const std::string& str) { m_strTarget = trimstr(str); }

	void SetValue(const std::string& str);
	void SetStdDev(const std::string& str);
	
	inline virtual UINT32 CalcBinaryRecordCount() const { return 1; }
	virtual void WriteBinaryMsr(std::ofstream* binary_stream, PUINT32 msrIndex) const;
	virtual UINT32 SetMeasurementRec(const vstn_t& binaryStn, it_vmsr_t& it_msr, it_vdbid_t& dbidmap);
	virtual void WriteDynaMLMsr(std::ofstream* dynaml_stream, const std::string& comment, bool) const;
	virtual void WriteDNAMsr(std::ofstream* dna_stream, const dna_msr_fields& dmw, const dna_msr_fields& dml, bool) const;
	virtual void SimulateMsr(vdnaStnPtr* vStations, const CDnaEllipsoid* ellipsoid);

public:
	std::string 	m_strTarget;
	double 	m_dValue;
	double 	m_dStdDev;
};

}	// namespace measurements
}	// namespace dynadjust


#endif /* DNAHEIGHTDIFFERENCE_H_ */
