//============================================================================
// Name         : dnadirection.hpp
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
// Description  : Interface for the CDnaDirection class
//============================================================================


#ifndef DNADIRECTION_H_
#define DNADIRECTION_H_

#if defined(_MSC_VER)
	#if defined(LIST_INCLUDES_ON_BUILD) 
		#pragma message("  " __FILE__) 
	#endif
#endif

#include <include/measurement_types/dnameasurement.hpp>
#include <include/functions/dnatemplatefuncs.hpp>

namespace dynadjust {
namespace measurements {

class CDnaDirection : public CDnaMeasurement
{
public:
	CDnaDirection(void);
	virtual ~CDnaDirection(void);

	// move constructor and move assignment operator
	CDnaDirection(CDnaDirection&& d);
	CDnaDirection& operator=(CDnaDirection&& rhs);

private:
	// disallowed in CDnaMeasurement
	//CDnaDirection(const CDnaDirection&);
	//CDnaDirection& operator=(const CDnaDirection& rhs);

public:
	//CDnaDirection(const bool bIgnore, const std::string& strFirst, const std::string& strTarget, const double& dValue, const double& dStdDev);

	//virtual inline CDnaDirection* clone() const { return new CDnaDirection(*this); }
	bool operator==(const CDnaDirection& rhs) const;
	bool operator<(const CDnaDirection& rhs) const;

	//inline CDnaDirection& operator[](int iIndex) { return this[iIndex]; }

	inline std::string GetTarget() const { return m_strTarget; }
	inline double GetValue() const { return m_drValue; }
	inline double GetStdDev() const { return m_dStdDev; }
	inline float GetInstrHeight() const { return m_fInstHeight; }
	inline float GetTargetHeight() const { return m_fTargHeight; }

	inline void SetRecordedTotal(const UINT32& l) { m_lRecordedTotal = l; }
	inline void SetClusterID(const UINT32& id) { m_lsetID = id; }
	inline void SetTarget(const std::string& str) { m_strTarget = trimstr(str); }
	void SetValue(const std::string& str);
	void SetStdDev(const std::string& str);
	void SetInstrumentHeight(const std::string& str);
	void SetTargetHeight(const std::string& str);

	inline virtual UINT32 CalcBinaryRecordCount() const { return 1; }
	//void coutDirectionData(std::ostream &os) const;
	virtual void WriteBinaryMsr(std::ofstream* binary_stream, PUINT32 msrIndex) const;
	virtual UINT32 SetMeasurementRec(const vstn_t& binaryStn, it_vmsr_t& it_msr, it_vdbid_t& dbidmap);
	virtual void WriteDynaMLMsr(std::ofstream* dynaml_stream, const std::string& comment, bool bSubMeasurement = false) const;
	virtual void WriteDNAMsr(std::ofstream* dna_stream, const dna_msr_fields& dmw, const dna_msr_fields& dml, bool bSubMeasurement = false) const;
	virtual void SimulateMsr(vdnaStnPtr* vStations, const CDnaEllipsoid* ellipsoid);

	//virtual inline void SetDatabaseMap_bmsIndex(const UINT32& bmsIndex) { m_msr_db_map.bms_index = bmsIndex; }
public:
	std::string	m_strTarget;

protected:
	double	m_drValue;
	double	m_dStdDev;
	float	m_fInstHeight;
	float	m_fTargHeight;
	UINT32	m_lRecordedTotal;
	UINT32	m_lsetID;
};





class CDnaAzimuth : public CDnaDirection
{
public:
	CDnaAzimuth(void);
	virtual ~CDnaAzimuth(void);

private:
	// disallowed in CDnaMeasurement
	//CDnaAzimuth(const CDnaAzimuth&);
	//CDnaAzimuth& operator=(const CDnaAzimuth& rhs);

public:

	//CDnaAzimuth(const bool bIgnore, const std::string& strFirst, const std::string& strTarget, const double& drValue, const double& dStdDev);

	bool operator==(const CDnaAzimuth& rhs) const;

	//inline CDnaAzimuth& operator[](int iIndex) { return this[iIndex]; }

protected:

};

}	// namespace measurements
}	// namespace dynadjust


#endif /* DNADIRECTION_H_ */
