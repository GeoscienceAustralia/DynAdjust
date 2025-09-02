//============================================================================
// Name         : dnadatum.hpp
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
// Description  : DynAdjust Datum Library
//============================================================================

#ifndef DNADATUM_H_
#define DNADATUM_H_

#if defined(_MSC_VER)
	#if defined(LIST_INCLUDES_ON_BUILD) 
		#pragma message("  " __FILE__) 
	#endif
#endif

/// \cond
#include <math.h>
#include <string>
#include <sstream>
/// \endcond

#include <include/config/dnatypes.hpp>
#include <include/config/dnaexports.hpp>
#include <include/functions/dnastrmanipfuncs.hpp>
#include <include/functions/dnatemplatedatetimefuncs.hpp>
#include <include/parameters/dnadatumprojectionparam.hpp>
#include <include/parameters/dnaellipsoid.hpp>

namespace dynadjust {
namespace datum_parameters {

class DNATYPE_API CDnaDatum
{
public:

	CDnaDatum(void);
	CDnaDatum(const UINT32& epsgCode);
	CDnaDatum(const UINT32& epsgCode, const boost::gregorian::date& epoch);
	CDnaDatum(const std::string& epsgCode, const std::string& epoch);
	virtual inline ~CDnaDatum(void) {}

private:	
	// Disallow use of compiler generated functions. 
	CDnaDatum(const CDnaDatum& newDatum);

	
	//virtual inline CDnaDatum* clone() const { 
	//	return new CDnaDatum(*this); 
	//}

public:
	CDnaDatum& operator=(const CDnaDatum& rhs);
	bool operator==(const CDnaDatum& rhs) const;
	
	//inline CDnaDatum& operator[](int iIndex) { return this[iIndex]; }

	inline CDnaEllipsoid GetEllipsoid() const { return ellipsoid_; }
	inline const CDnaEllipsoid* GetEllipsoidRef() const { return &ellipsoid_; }

	std::string GetName() const { return datumName_; }
	std::string GetEpsgCode_s() const;
	std::string GetEpoch_s() const;
	inline UINT32 GetEpsgCode_i() const { return epsgCode_; }
	inline boost::gregorian::date GetEpoch() const { return epoch_; }
	inline DATUM_TYPE GetType() const { return datumType_; }
	inline bool isStatic() const { return datumType_ == STATIC_DATUM; }
	inline bool isDynamic() const { return datumType_ == DYNAMIC_DATUM; }

	inline void SetEpoch(const boost::gregorian::date& epoch) { epoch_ = epoch; }
	void SetEpoch(const std::string& epoch);
	//void SetEpoch(const double& decimal_year);

	void SetDatum(const UINT32& epsgCode);
	void SetDatum(const UINT32& epsgCode, const boost::gregorian::date& epoch);
	void SetDatum(const std::string& epsgCode);
	//void SetDatum(const std::string& epsgCode, const boost::gregorian::date& epoch);

	void SetDatumFromName(const std::string& datumName, const std::string& epoch);
	void SetDatumFromEpsg(const std::string& epsgCode, const std::string& epoch);

	//bool isSameFrame(const CDnaDatum& rhs) const;

private:
	void initialiseDatumFromEpsgCode();
	
	UINT32			epsgCode_;
	CDnaEllipsoid	ellipsoid_;
	boost::gregorian::date	epoch_;
	DATUM_TYPE		datumType_;
	std::string			datumName_;

};


}	// namespace datum_parameters
}	// namespace dynadjust

#endif /* DNADATUMPROJECTION_H_ */
