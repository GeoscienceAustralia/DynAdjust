// Not copyrighted - public domain.
//
// This sample parser implementation was generated by CodeSynthesis XSD,
// an XML Schema to C++ data binding compiler. You may use it in your
// programs without any restrictions.
//

#ifndef DYNA_ML_PIMPL_HXX
#define DYNA_ML_PIMPL_HXX

#if defined(_MSC_VER)
#if defined(LIST_INCLUDES_ON_BUILD) 
#pragma message("  " __FILE__) 
#endif
#endif

#include <dynadjust/dnaimport/dnaparser_pskel.hxx>
#include <include/measurement_types/dnameasurement_types.hpp>
#include <include/functions/dnatemplatestnmsrfuncs.hpp>

using namespace dynadjust::measurements;

class Clusterpoint_pimpl: public virtual Clusterpoint_pskel
{
public:
	virtual void pre ();
	virtual void X (const ::std::string&);
	virtual void Y (const ::std::string&);
	virtual void Z (const ::std::string&);
	virtual void SigmaXX (const ::std::string&);
	virtual void SigmaXY (const ::std::string&);
	virtual void SigmaXZ (const ::std::string&);
	virtual void SigmaYY (const ::std::string&);
	virtual void SigmaYZ (const ::std::string&);
	virtual void SigmaZZ (const ::std::string&);
	virtual void PointCovariance ();
	virtual void MeasurementDBID (const ::std::string&);
	virtual void post_Clusterpoint ();

protected:
	std::string _X;
	std::string _Y;
	std::string _Z;
	std::string _SigmaXX;
	std::string _SigmaXY;
	std::string _SigmaXZ;
	std::string _SigmaYY;
	std::string _SigmaYZ;
	std::string _SigmaZZ;
};

class Directions_pimpl: public virtual Directions_pskel
{
public:
	virtual void pre ();
	virtual void Ignore (const ::std::string&);
	virtual void Target (const ::std::string&);
	virtual void Value (const ::std::string&);
	virtual void StdDev (const ::std::string&);
	virtual void MeasurementDBID (const ::std::string&);
	virtual void post_Directions (const UINT32&);

protected:
	std::string _Ignore;
	std::string _Target;
	std::string _Value;
	std::string _StdDev;

};

class DnaMeasurement_pimpl: public virtual DnaMeasurement_pskel
{
public:
	virtual void pre ();
	virtual void Type (const ::std::string&);
	virtual void Ignore (const ::std::string&);
	virtual void First (const ::std::string&);
	virtual void Second (const ::std::string&);
	virtual void Third (const ::std::string&);
	virtual void Value (const ::std::string&);
	virtual void StdDev (const ::std::string&);
	virtual void InstHeight (const ::std::string&);
	virtual void TargHeight (const ::std::string&);
	virtual void Total (const ::std::string&);
	virtual void Directions ();
	virtual void Vscale (const ::std::string&);
	virtual void Epoch (const ::std::string&);
	virtual void ReferenceFrame (const ::std::string&);
	virtual void GPSBaseline ();
	virtual void Hscale (const ::std::string&);
	virtual void Lscale (const ::std::string&);
	virtual void Pscale (const ::std::string&);
	virtual void Clusterpoint ();
	virtual void Coords (const ::std::string&);
	virtual void Source (const ::std::string&);
	virtual void MeasurementDBID (const ::std::string&);
	virtual void ClusterDBID (const ::std::string&);
	virtual void post_DnaMeasurement ();
};

class DnaStation_pimpl: public virtual DnaStation_pskel
{
public:
	virtual void pre ();
	virtual void Name (const ::std::string&);
	virtual void Constraints (const ::std::string&);
	virtual void Type (const ::std::string&);
	virtual void StationCoord ();
	virtual void Description (const ::std::string&);
	virtual void post_DnaStation ();
};

class DnaXmlFormat_pimpl: public virtual DnaXmlFormat_pskel
{
public:
	DnaXmlFormat_pimpl(std::ifstream* is, PUINT32 clusterID, 
		const std::string& referenceframe, const std::string& epoch,
		bool firstFile,
		bool userspecifiedreferenceframe, bool userspecifiedepoch,
		bool overridereferenceframe);
	virtual void pre ();
	virtual void DnaStation ();
	virtual void DnaMeasurement ();
	virtual void type ();
	virtual void referenceframe ();
	virtual void epoch ();
	virtual void post_DnaXmlFormat (vdnaStnPtr* vStations, vdnaMsrPtr* vMeasurements);
	inline bool filespecifiedreferenceframe() { return _filespecifiedreferenceframe; }
	inline bool filespecifiedepoch() { return _filespecifiedepoch; }

	std::string DnaXmlParseMessage() { return _parse_msg; }

	inline UINT32 NumStationsRead() const { return _station_count; }
	inline UINT32 NumMeasurementsRead() const { return _measurement_count; }
	inline UINT32 CurrentClusterID() const { return _clusterID; }

	inline std::string FileEpsg() const { return _fileEpsg; }
	inline std::string FileEpoch() const { return _fileEpoch; }

protected:
	std::string _parse_msg;

private:
	std::ifstream* is_;
};

class GPSBaseline_pimpl: public virtual GPSBaseline_pskel
{
public:
	virtual void pre ();
	virtual void X (const ::std::string&);
	virtual void Y (const ::std::string&);
	virtual void Z (const ::std::string&);
	virtual void MeasurementDBID (const ::std::string&);
	virtual void SigmaXX (const ::std::string&);
	virtual void SigmaXY (const ::std::string&);
	virtual void SigmaXZ (const ::std::string&);
	virtual void SigmaYY (const ::std::string&);
	virtual void SigmaYZ (const ::std::string&);
	virtual void SigmaZZ (const ::std::string&);
	virtual void GPSCovariance ();
	virtual void post_GPSBaseline ();
};

class GPSCovariance_pimpl: public virtual GPSCovariance_pskel
{
public:
	virtual void pre ();
	virtual void m11 (const ::std::string&);
	virtual void m12 (const ::std::string&);
	virtual void m13 (const ::std::string&);
	virtual void m21 (const ::std::string&);
	virtual void m22 (const ::std::string&);
	virtual void m23 (const ::std::string&);
	virtual void m31 (const ::std::string&);
	virtual void m32 (const ::std::string&);
	virtual void m33 (const ::std::string&);
	virtual void post_GPSCovariance ();
};

class Height_pimpl: public virtual Height_pskel, public ::xml_schema::string_pimpl
{
public:
	virtual void pre ();
	virtual void system (const ::std::string&);
	virtual void post_Height ();
};

class PointCovariance_pimpl: public virtual PointCovariance_pskel
{
public:
	virtual void pre ();
	virtual void m11 (const ::std::string&);
	virtual void m12 (const ::std::string&);
	virtual void m13 (const ::std::string&);
	virtual void m21 (const ::std::string&);
	virtual void m22 (const ::std::string&);
	virtual void m23 (const ::std::string&);
	virtual void m31 (const ::std::string&);
	virtual void m32 (const ::std::string&);
	virtual void m33 (const ::std::string&);
	virtual void post_PointCovariance ();
};

class StationCoord_pimpl: public virtual StationCoord_pskel
{
public:
	virtual void pre ();
	virtual void Name (const ::std::string&);
	virtual void XAxis (const ::std::string&);
	virtual void YAxis (const ::std::string&);
	virtual void Height (const ::std::string&);
	virtual void GeoidModel ();
	virtual void HemisphereZone (const ::std::string&);
	virtual void post_StationCoord ();
};

class GeoidModel_pimpl: public virtual GeoidModel_pskel
{
public:
	virtual void pre ();
	virtual void NValue (const ::std::string&);
	virtual void DeflPM (const ::std::string&);
	virtual void DeflPV (const ::std::string&);
	virtual void post_GeoidModel ();
};

class type_pimpl: public virtual type_pskel,
	public ::xml_schema::string_pimpl
{
public:
	virtual void pre ();
	virtual void post_type ();
};

class referenceframe_pimpl: public virtual referenceframe_pskel,
	public ::xml_schema::string_pimpl
{
public:
	virtual void pre ();
	virtual void post_type (std::string& referenceframe, std::string& fileEpsg,
		bool userspecifiedreferenceframe, bool firstFile);
	inline const std::string str() { return _referenceframe; }

protected:
	std::string _referenceframe;
	bool _overridereferenceframe;
};

class epoch_pimpl: public virtual epoch_pskel,
	public ::xml_schema::string_pimpl
{
public:
	virtual void pre ();
	virtual void post_type (std::string& epoch, std::string& fileEpoch,
		bool userspecifiedreferenceframe, bool userspecifiedepoch, bool firstFile);
	inline const std::string str() { return _epoch; }

protected:
	std::string _epoch;
};

class system_pimpl: public virtual system_pskel,
	public ::xml_schema::string_pimpl
{
public:
	virtual void pre ();
	virtual void post_system ();
};

#endif // DYNA_ML_PIMPL_HXX
