//============================================================================
// Name         : dnameasurement.hpp
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
// Description  : Interface for the msr_t struct and CDnaMeasurement, 
//                CDnaCovariance and MsrTally classes
//============================================================================

#ifndef DNAMEASUREMENT_H
#define DNAMEASUREMENT_H

#if defined(_MSC_VER)
	#if defined(LIST_INCLUDES_ON_BUILD) 
		#pragma message("  " __FILE__) 
	#endif
#endif

/// \cond
#include <stdio.h>
#include <string.h>

#include <vector>
#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>

#include <memory>
/// \endcond

#include <include/functions/dnastrmanipfuncs.hpp>
#include <include/config/dnatypes.hpp>
#include <include/measurement_types/dnastation.hpp>
#include <include/math/dnamatrix_contiguous.hpp>
#include <include/parameters/dnadatum.hpp>


//#include <boost/random.hpp>

enum {
	INST_WIDTH = 12,
	TARG_WIDTH = 12,
	STN_WIDTH = 16,
	MEAS_WIDTH = 18,
	MSR2_WIDTH = 12,
	VAR_WIDTH = 13,
	COV_WIDTH = 13
};

namespace dynadjust { namespace measurements {

template<typename C, typename S>
S measurement_name(const C& type)
{
	switch (type)
	{
	case 'A':
		return "(A) Horizontal angle";
	case 'B':
		return "(B) Geodetic azimuth";
	case 'K':
		return "(K) Astronomic azimuth";
	case 'C':
		return "(C) Chord distance";
	case 'E':
		return "(E) Ellipsoid arc";
	case 'M':
		return "(M) Mean sea level arc";
	case 'S':
		return "(S) Slope distance";
	case 'D':
		return "(D) Direction set";
	case 'G':
		return "(G) GPS baseline";
	case 'X':
		return "(X) GPS baseline cluster";
	case 'H':
		return "(H) Orthometric height";
	case 'I':
		return "(I) Astronomic latitude";
	case 'J':
		return "(J) Astronomic longitude";
	case 'P':
		return "(P) Geodetic latitude";
	case 'Q':
		return "(Q) Geodetic longitude";
	case 'L':
		return "(L) Level difference";
	case 'R':
		return "(R) Ellipsoidal height";
	case 'V':
		return "(V) Zenith distance";
	case 'Z':
		return "(Z) Vertical angle";
	case 'Y':
		return "(Y) GPS point cluster";
	}
	return "";
}


// forward declarations
class CDnaGpsBaseline;
class CDnaGpsPoint;
class CDnaDirection;
class CDnaMeasurement;
class CDnaCovariance;

// measurement types
typedef std::shared_ptr<CDnaMeasurement> dnaMsrPtr;
typedef std::vector<dnaMsrPtr> vdnaMsrPtr, *pvdnaMsrPtr;						// vector of dnaMsrPtr
typedef vdnaMsrPtr::iterator _it_vdnamsrptr;
typedef vdnaMsrPtr::const_iterator _it_vdnamsrptr_const;

typedef std::shared_ptr<CDnaCovariance> dnaCovariancePtr;
typedef std::shared_ptr< std::vector<CDnaCovariance> > vecCovariancePtr;

// data struct for storing measurement information to binary measurement file
typedef struct msr_t {
	msr_t()
		: measType('\0'), measStart(0), measurementStations(1), ignore(false), station1(0), station2(0)
		, station3(0), vectorCount1(0), vectorCount2(0), clusterID(0), fileOrder(0)
		, term1(0.), term2(0.), term3(0.), term4(0.), scale1(1.), scale2(1.), scale3(1.), scale4(1.)
		, measAdj(0.), measCorr(0.), measAdjPrec(0.), residualPrec(0.)
		, NStat(0.), TStat(0.), PelzerRel(0.)
		, preAdjCorr(0.), preAdjMeas(0.) 
	{
			memset(coordType, '\0', sizeof(coordType));
			// GDA94, lat, long, height
			memset(epsgCode, '\0', sizeof(epsgCode));
			snprintf(epsgCode, sizeof(epsgCode), DEFAULT_EPSG_S);
			memset(epoch, '\0', sizeof(epoch));
	}

	char	measType;				// 'A', 'S', 'X', ... , etc.
	char	measStart;				// Start of a measurement (0=start or X, 1=Y, 2=Z, 3=covX, 4=covY, 5=covZ)
	char	measurementStations;	// One-, two- or three-station measurement
	char	epsgCode[7];			// epsg ID, i.e. NNNNN (where NNNNN is in the range 0-32767)
	char	epoch[12];				// date, i.e. "DD.MM.YYYY" (10 chars)
									// if datum is dynamic, Epoch is YYYY MM DD
									// if datum is static, Epoch is ignored
	char	coordType[4];			// "LLH", "UTM", ... , etc.
	bool	ignore;
	UINT32	station1;        		// stations 1, 2 and 3 are indices to
	UINT32	station2;       		// record numbers in the station
	UINT32	station3;				// record file.
	UINT32	vectorCount1;			// number of directions, GpsPoints or GpsBaselines
	UINT32	vectorCount2;			// number of covariances for GpsPoint or GpsBaseline
									// number of non-ignored directions
	UINT32	clusterID;				// cluster ID (which cluster this measurement belongs to)
	UINT32	fileOrder;				// original file order
	double	term1;					// measurement, X, Y, Z, dX, dY, dZ value
									// direction
	double	term2;					// measurement, XX, XY or XZ variance
									// direction variance
	double	term3;					// instrument height, YY or YZ variance
	double	term4;					// target height or ZZ variance
	double	scale1;					// phi, n or X scalar
									// derived angle corrected for deflection of the vertical
	double	scale2;					// lambda, e or Y scalar
									// derived angle variance
	double	scale3;					// height, up or Z scalar
									// derived angle covariance
	double	scale4;					// matrix scalar
	double	measAdj;
	double	measCorr;
	double	measAdjPrec;
	double	residualPrec;
	double	NStat;
	double	TStat;
	double	PelzerRel;
	double	preAdjCorr;
	double	preAdjMeas;
} measurement_t;

typedef std::vector<measurement_t> vmsr_t, *pvmsr_t;
typedef vmsr_t::iterator it_vmsr_t, *pit_vmsr_t;
typedef vmsr_t::const_iterator it_vmsr_t_const;

class CDnaCovariance
{
public:
	CDnaCovariance(void);
	virtual ~CDnaCovariance();

	// move constructor and move assignment operator
	CDnaCovariance(CDnaCovariance&& c);
	CDnaCovariance& operator=(CDnaCovariance&& rhs);

private:
	// disallow copying
	CDnaCovariance(const CDnaCovariance& newCovariance);
	CDnaCovariance& operator=(const CDnaCovariance& rhs);
	
public:
	//virtual inline CDnaCovariance* clone() const { return new CDnaCovariance(*this); }
	bool operator==(const CDnaCovariance& rhs) const;

	//inline CDnaCovariance& operator[](int iIndex) { return this[iIndex]; }

	inline void SetType(const std::string& str) { m_strType = trimstr(str); }
	inline char GetTypeC() const { return (m_strType.c_str())[0]; }

	// m_bIgnore used only to 'split' cluster measurements
	inline void SetIgnore(const bool bval) { m_bIgnore = bval; }
	inline bool GetIgnore() const { return m_bIgnore; }

	inline virtual UINT32 CalcBinaryRecordCount() const { return 3; }
	void WriteBinaryMsr(std::ofstream* binary_stream, PUINT32 msrIndex, const std::string& epsgCode, const std::string& epoch) const;
	virtual UINT32 SetMeasurementRec(const vstn_t& binaryStn, it_vmsr_t& it_msr);
	virtual void WriteDynaMLMsr(std::ofstream* dynaml_stream) const;
	virtual void WriteDNAMsr(std::ofstream* dna_stream, 
		const dna_msr_fields& dmw, const dna_msr_fields& dml) const;
	virtual void SimulateMsr(vdnaStnPtr*, const CDnaEllipsoid*);
	
	void SerialiseDatabaseMap(std::ofstream* os, const msr_database_id_map& dbid);

	inline void SetClusterID(const UINT32& id) { m_lclusterID = id; }
	inline void SetStn1Index(const UINT32& stn) { m_lstn1Index = stn; }
	inline void SetStn2Index(const UINT32& stn) { m_lstn2Index = stn; }
	void SetM11(const std::string& str);
	void SetM12(const std::string& str);
	void SetM13(const std::string& str);
	void SetM21(const std::string& str);
	void SetM22(const std::string& str);
	void SetM23(const std::string& str);
	void SetM31(const std::string& str);
	void SetM32(const std::string& str);
	void SetM33(const std::string& str);

	inline void SetM11(const double& dbl) { m_dM11 = dbl; }
	inline void SetM12(const double& dbl) { m_dM12 = dbl; }
	inline void SetM13(const double& dbl) { m_dM13 = dbl; }
	inline void SetM21(const double& dbl) { m_dM21 = dbl; }
	inline void SetM22(const double& dbl) { m_dM22 = dbl; }
	inline void SetM23(const double& dbl) { m_dM23 = dbl; }
	inline void SetM31(const double& dbl) { m_dM31 = dbl; }
	inline void SetM32(const double& dbl) { m_dM32 = dbl; }
	inline void SetM33(const double& dbl) { m_dM33 = dbl; }

	inline UINT32 GetClusterID() const { return m_lclusterID; }
	inline UINT32 GetStn1Index() const { return m_lstn1Index; }
	inline UINT32 GetStn2Index() const { return m_lstn2Index; }
	inline double GetM11() const { return m_dM11; }
	inline double GetM12() const { return m_dM12; }
	inline double GetM13() const { return m_dM13; }
	inline double GetM21() const { return m_dM21; }
	inline double GetM22() const { return m_dM22; }
	inline double GetM23() const { return m_dM23; }
	inline double GetM31() const { return m_dM31; }
	inline double GetM32() const { return m_dM32; }
	inline double GetM33() const { return m_dM33; }

protected:
	bool	m_bIgnore;
	UINT32	 m_lstn1Index;        			// This is an index to the record number in the station file
	UINT32	 m_lstn2Index;
	std::string	 m_strType;
	double	 m_dM11;
	double	 m_dM12;
	double	 m_dM13;
	double	 m_dM21;
	double	 m_dM22;
	double	 m_dM23;
	double	 m_dM31;
	double	 m_dM32;
	double	 m_dM33;
	UINT32	m_lclusterID;
};



// Abstract class for all measurement types
class CDnaMeasurement
{
public:
	CDnaMeasurement();
	virtual ~CDnaMeasurement();

	// move constructor and move assignment operator
	CDnaMeasurement(CDnaMeasurement&& m);
	CDnaMeasurement& operator=(CDnaMeasurement&& rhs);

private:
	// disallow copying
	CDnaMeasurement(const CDnaMeasurement&);
	CDnaMeasurement& operator=(const CDnaMeasurement& rhs);

public:
	//virtual CDnaMeasurement* clone() const = 0;  // The Virtual (Copy) Constructor

	inline std::string GetType() const { return m_strType; }
	inline char GetTypeC() const { return (m_strType.c_str())[0]; }
	inline bool GetIgnore() const { return m_bIgnore; }
	inline bool NotIgnored() const { return m_bIgnore == false; }
	inline std::string GetFirst() const { return m_strFirst; }
	inline MEASUREMENT_STATIONS GetMsrStnCount() const { return m_MSmeasurementStations; }

	inline std::string GetEpsg() const { return m_epsgCode; }
	inline std::string GetSource() const { return m_sourceFile; }

	inline bool GetInsufficient() const { return m_bInsufficient; }

	void SetType(const std::string& str);
	inline void SetIgnore(const bool bval) { m_bIgnore = bval; }
	inline void SetFirst(const std::string& str) { m_strFirst = trimstr(str); }
	
	inline void SetEpsg(const std::string& e) { m_epsgCode = trimstr(e); }
	inline void SetSource(const std::string& source) { m_sourceFile = source; }
	
	inline void SetInsufficient(const bool bval) { m_bInsufficient = bval; }

	inline void SetMsrIndex(const UINT32& order) { m_lmeasurementIndex = order; }
	inline void SetStn1Index(const UINT32& order) { m_lstn1Index = order; }
	inline void SetStn2Index(const UINT32& order) { m_lstn2Index = order; }
	inline void SetStn3Index(const UINT32& order) { m_lstn3Index = order; }

	inline UINT32 GetMsrIndex() const { return m_lmeasurementIndex; }
	//inline UINT32& GetBMSFIndex() const { return m_lmeasurementIndex * sizeof(measurement_t); }
	inline UINT32 GetStn1Index() const { return m_lstn1Index; }
	inline UINT32 GetStn2Index() const { return m_lstn2Index; }
	inline UINT32 GetStn3Index() const { return m_lstn3Index; }

	inline double GetMeasAdj() const { return m_measAdj; }
	inline double GetMeasCorr() const { return m_measCorr; }
	inline double GetMeasAdjPrec() const { return m_measAdjPrec; }
	inline double GetResidualPrec() const { return m_residualPrec; }
	inline double GetPreAdjCorr() const { return m_preAdjCorr; }

	// pure virtual functions overridden by specialised classes
	virtual UINT32 CalcBinaryRecordCount() const = 0;
	virtual void WriteBinaryMsr(std::ofstream* binary_stream, PUINT32 msrIndex) const = 0;
	virtual UINT32 SetMeasurementRec(const vstn_t&, it_vmsr_t& it_msr, it_vdbid_t& dbidmap) = 0;
	virtual void WriteDynaMLMsr(std::ofstream* dynaml_stream, const std::string& comment, bool bSubMeasurement = false) const = 0;
	virtual void WriteDNAMsr(std::ofstream* dna_stream, const dna_msr_fields& dmw, const dna_msr_fields& dml, bool bSubMeasurement = false) const = 0;
	virtual void SimulateMsr(vdnaStnPtr* vStations, const CDnaEllipsoid* ellipsoid) = 0;

	// A function used by CDnaGpsPoint and CDnaGpsPointCluster only.
	// Used to export latest station coordinate and variance estimates to Y measurement
	virtual void PopulateMsr(pvstn_t, uint32_uint32_map*, vUINT32*,
		const UINT32&, const CDnaDatum*, math::matrix_2d*, math::matrix_2d*) {}

	// virtual functions overridden by specialised classes
	virtual inline UINT32 GetClusterID() const { return 0; }
	virtual inline std::string GetCoordType() const { return m_strType; }
	virtual inline float GetInstrHeight() const { return 0; }
	virtual inline double GetVscale() const { return 0; }
	virtual inline double GetPscale() const { return 0; }
	virtual inline double GetLscale() const { return 0; }
	virtual inline double GetHscale() const { return 0; }
	virtual inline double GetStdDev() const { return 0; }
	virtual inline std::string GetTarget() const { return ""; }
	virtual inline std::string GetTarget2() const { return ""; }
	virtual inline float GetTargetHeight() const { return 0; }
	virtual inline UINT32 GetTotal() const { return 0; }
	virtual inline double GetValue() const { return 0; }

	virtual inline std::string GetReferenceFrame() const { return ""; }
	inline std::string GetEpoch() const { return m_epoch; }
	
	virtual inline std::vector<CDnaGpsBaseline>* GetBaselines_ptr() { return 0; }
	virtual inline std::vector<CDnaDirection>* GetDirections_ptr() { return 0; }
	virtual inline std::vector<CDnaGpsPoint>* GetPoints_ptr() { return 0; }
	virtual inline std::vector<CDnaCovariance>* GetCovariances_ptr() { return 0; }

	virtual void AddDirection(const CDnaMeasurement*) {}
	virtual void AddGpsBaseline(const CDnaMeasurement*) {}
	virtual void AddGpsPoint(const CDnaMeasurement*) {}
	virtual void AddGpsCovariance(const CDnaCovariance*) {}
	virtual void AddPointCovariance(const CDnaCovariance*) {}

	virtual void ReserveDirectionsCount(const UINT32&) {}
	virtual void ReserveGpsBaselinesCount(const UINT32&) {}
	virtual void ReserveGpsPointsCount(const UINT32&) {}
	virtual void ReserveGpsCovariancesCount(const UINT32&) {}
	virtual void ResizeGpsCovariancesCount(const UINT32&) {}

	virtual void SetRecordedTotal(const UINT32&) {}
	virtual void SetNonIgnoredDirns(const UINT32&) {}
	virtual void SetClusterID(const UINT32&) {}
	virtual void SetCoordType(const std::string&) {}
	virtual void SetHscale(const std::string&) {}
	virtual void SetHscale(const double&) {}
	virtual void SetInstrumentHeight(const std::string&) {}
	
	virtual void SetReferenceFrame(const std::string&) {}
	void SetEpoch(const std::string& epoch);
	
	virtual void SetLscale(const std::string&) {}
	virtual void SetLscale(const double&) {}
	virtual void SetPscale(const std::string&) {}
	virtual void SetPscale(const double&) {}
	virtual void SetSigmaXX(const std::string&) {}
	virtual void SetSigmaXY(const std::string&) {}
	virtual void SetSigmaXZ(const std::string&) {}
	virtual void SetSigmaYY(const std::string&) {}
	virtual void SetSigmaYZ(const std::string&) {}
	virtual void SetSigmaZZ(const std::string&) {}
	virtual void SetStdDev(const std::string&) {}
	virtual void SetTarget(const std::string&) {}
	virtual void SetTarget2(const std::string&) {}
	virtual void SetTargetHeight(const std::string&) {}
	virtual void SetTotal(const std::string&) {}
	virtual void SetValue(const std::string&) {}
	virtual void SetVscale(const std::string&) {}
	virtual void SetVscale(const double&) {}
	virtual void SetX(const std::string&) {}
	virtual void SetY(const std::string&) {}
	virtual void SetZ(const std::string&) {}

	virtual void SetXAxis(const double&) {}
	virtual void SetYAxis(const double&) {}
	virtual void SetZAxis(const double&) {}
	virtual void SetSigmaXX(const double&) {}
	virtual void SetSigmaXY(const double&) {}
	virtual void SetSigmaXZ(const double&) {}
	virtual void SetSigmaYY(const double&) {}
	virtual void SetSigmaYZ(const double&) {}
	virtual void SetSigmaZZ(const double&) {}
	virtual void SetTotal(const UINT32&) {}
	
	virtual void PreferGMeasurements() {}

	//virtual void coutBaselineData(std::ostream &os, const int& pad, const UINT16& uType = 0) {}

	void SetMeasurementDBID(const std::string& str);
	void SetClusterDBID(const std::string& str);
	
	void SetClusterDBID(const UINT32& u, bool s);

	inline UINT32 GetClusterDBID() { return m_msr_db_map.cluster_id; }
	inline bool GetClusterDBIDset() { return m_msr_db_map.is_cls_id_set; }
	
	void SetDatabaseMap(const msr_database_id_map& dbidmap);
	virtual void SetDatabaseMaps(it_vdbid_t& it_dbidmap);

	virtual void SerialiseDatabaseMap(std::ofstream* os);

protected:
	void coutMeasurement(std::ostream &os) const;

public:
	std::string	m_strFirst;
	MEASUREMENT_STATIONS m_MSmeasurementStations;
	
protected:
	std::string	m_strType;
	bool	m_bIgnore;

	UINT32 	m_lmeasurementIndex;
	UINT32	m_lstn1Index;        			// stations 1, 2 and 3 are indices to
	UINT32	m_lstn2Index;       			// record numbers in the station
	UINT32	m_lstn3Index;					// record file.

	double	m_measAdj;
	double	m_measCorr;
	double	m_measAdjPrec;
	double	m_residualPrec;
	double	m_preAdjCorr;

	std::string	m_epsgCode;
	std::string	m_sourceFile;

	std::string	m_epoch;
	
	msr_database_id_map		m_msr_db_map;

	bool	m_bInsufficient;
};

// In the event a new measurement type is added, ensure SUPPORTED_MSR_COUNT is
// updated accordingly
class DNATYPE_API MsrTally
{
public:
	MsrTally();

	inline UINT32 SupportedMsrTypes() {
		return SUPPORTED_MSR_COUNT;
	}

	static void FillMsrList(vchar& msr_list);
	static std::string GetMsrName(const char& c);
	
	void initialise();

	MsrTally& operator+=(const MsrTally& rhs);
	//MsrTally& operator-=(const MsrTally& rhs);
	MsrTally operator+(const MsrTally& rhs) const;
	//MsrTally operator-(const MsrTally& rhs) const;
	UINT32 TotalCount();
	void coutSummary(std::ostream &os, const std::string& title);
	UINT32 MeasurementCount(const char& msrType);
	
	void CreateTally(const vdnaMsrPtr& vMeasurements);
	void CreateTally(const vmsr_t& vMeasurements, const vUINT32& CML);
	UINT32 CreateTally(const vmsr_t& vMeasurements, bool countValidOnly=false);
	
	void IncrementMsrType(const char& msrType, const UINT32& count=1);
	
	void coutSummaryMsrToStn(std::ostream &os, const std::string& station);
	void coutSummaryMsrToStn_Compressed(std::ostream &os, const std::string& station);

	//bool GPSOnly();	
	inline bool ContainsNonGPS() { return containsNonGPS; }

	static _MEASUREMENT_STATIONS_ Stations(const char& measType);
	
	bool containsNonGPS;

	UINT32 A, B, C, D, E, G, H, I, J, K, L, M, P, Q, R, S, V, X, Y, Z;
	UINT32 ignored;
	UINT32 SUPPORTED_MSR_COUNT;
	UINT32 totalCount;
};

typedef std::vector<MsrTally> vmsrtally;
typedef vmsrtally::iterator _it_vmsrtally;

template <typename T>
void MsrToStnSummaryHeaderLine(
	T& stream)
{
	UINT32 i, j(STATION + (NUMERIC_WIDTH * 20) + STAT);
	for (i=0; i<j; ++i)
		stream << "-";

	stream << std::endl;
}

template <typename T>
void MsrToStnSummaryHeader(
	T& stream, std::string& header)
{
	stream << std::endl << header << std::endl;
	stream << "------------------------------------------" << std::endl << std::endl;

	stream << std::setw(STATION) << std::left << "Station" <<
		std::setw(NUMERIC_WIDTH) << std::right <<	"A" <<
		std::setw(NUMERIC_WIDTH) << std::right <<	"B" <<
		std::setw(NUMERIC_WIDTH) << std::right <<	"C" <<
		std::setw(NUMERIC_WIDTH) << std::right <<	"D" <<
		std::setw(NUMERIC_WIDTH) << std::right <<	"E" <<
		std::setw(NUMERIC_WIDTH) << std::right << "G" <<
		std::setw(NUMERIC_WIDTH) << std::right <<	"H" <<
		std::setw(NUMERIC_WIDTH) << std::right <<	"I" <<
		std::setw(NUMERIC_WIDTH) << std::right <<	"J" <<
		std::setw(NUMERIC_WIDTH) << std::right <<	"K" <<
		std::setw(NUMERIC_WIDTH) << std::right <<	"L" <<
		std::setw(NUMERIC_WIDTH) << std::right <<	"M" <<
		std::setw(NUMERIC_WIDTH) << std::right <<	"P" <<
		std::setw(NUMERIC_WIDTH) << std::right <<	"Q" <<
		std::setw(NUMERIC_WIDTH) << std::right <<	"R" <<
		std::setw(NUMERIC_WIDTH) << std::right <<	"S" <<
		std::setw(NUMERIC_WIDTH) << std::right <<	"V" <<
		std::setw(NUMERIC_WIDTH) << std::right <<	"X" <<
		std::setw(NUMERIC_WIDTH) << std::right <<	"Y" <<
		std::setw(NUMERIC_WIDTH) << std::right <<	"Z" <<
		// Total
		std::setw(STAT) << std::right <<	"Total" <<
	 	std::endl;

	MsrToStnSummaryHeaderLine(stream);

}




}	// namespace measurements
}	// namespace dynadjust


#endif // ifndef DNAMEASUREMENT_H
