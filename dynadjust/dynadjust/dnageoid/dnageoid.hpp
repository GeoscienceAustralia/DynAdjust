//============================================================================
// Name         : dnageoid.hpp
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
// Description  : AusGeoid Grid File (NTv2) Interpolation library
//============================================================================

#ifndef DNAGEOID_H_
#define DNAGEOID_H_

#if defined(_MSC_VER)
	#if defined(LIST_INCLUDES_ON_BUILD) 
		#pragma message("  " __FILE__) 
	#endif
#endif

/// \cond
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <exception>
#include <new>
#include <math.h>

/// \endcond

#include <include/config/dnaversion.hpp>
#include <include/io/bst_file_loader.hpp>
#include <include/io/bst_file.hpp>
#include <include/config/dnaconsts.hpp>
#include <include/config/dnaexports.hpp>
#include <include/config/dnaoptions.hpp>
#include <include/config/dnatypes.hpp>
#include <include/config/dnaversion-stream.hpp>

#include <include/parameters/dnaepsg.hpp>

#include <include/measurement_types/dnastation.hpp>

#include <include/exception/dnaexception.hpp>

#include <include/functions/dnastrutils.hpp>
#include <include/functions/dnastrmanipfuncs.hpp>
#include <include/functions/dnatemplatecalcfuncs.hpp>
#include <include/functions/dnaiostreamfuncs.hpp>
#include <include/functions/dnastringfuncs.hpp>
#include <include/functions/dnatemplatedatetimefuncs.hpp>
#include <include/functions/dnatemplatefuncs.hpp>

// When building dnageoid.dll with standard C API:
//  1. uncomment to define BUILD_DNAGEOID_C_DLL 

	   //#define BUILD_DNAGEOID_C_DLL

//
//  2. add def file to project
//     .\dynadjust_1_00_00\dynadjust\dnageoid\dnageoid_ext.def
//     dnageoid > Properties > Configuration Properties > Linker > Input > Module Definition File = dnageoid_ext.def
//
//     Note that this def file doesn't need to be included in order to build driver, but
//     should be included to permit backward compatibility.  That is, by using the def file, new versions 
//     of dnageoid can be rebuilt without having to rebuild the driver.


#include <dynadjust/dnageoid/bicubic_interpolation.hpp>
#include <dynadjust/dnageoid/dnageoid_ext.hpp>

#include <filesystem>

using namespace dynadjust::measurements;
using namespace dynadjust::exception;
using namespace dynadjust::iostreams;

// grid node format identifiers for ASCII NTv2 grid file
const char* const SHIFTS = "%f%f%f%f";
// grid node line length for ASCII NTv2 grid file
const UINT32 ASCII_LINE_LENGTH = 40;
	
#define GSB		"gsb"
#define ASC		"asc"
#define CSV		"csv"
#define DAT		"dat"
#define TXT		"txt"
#define PRN		"prn"

#define TYPE_ASC	0
#define TYPE_GSB	1
#define TYPE_CSV	2
#define TYPE_DAT	3

#define AGD84 0
#define GDA94 1
#define GEOGRAPHIC 0
#define PROJECTED  1

#define HP_NOTATION 1
#define NULL_ACCURACY -1000.0

#define ANS 0
#define GRS80 1
#define ORTHOMETRIC 0
#define ELLIPSOIDAL 1

namespace dynadjust {
namespace geoidinterpolation {

// This class is exported from the dnageoid.dll
#ifdef _MSC_VER

// Is a dll with standard C API being built?
#ifdef BUILD_DNAGEOID_C_DLL
class dna_geoid_interpolation {
#else
class DNAGEOID_API dna_geoid_interpolation {
#endif

#else
	// Non MS Windows
class dna_geoid_interpolation {
#endif

public:
	dna_geoid_interpolation();
	virtual ~dna_geoid_interpolation();

private:
	// Disallow use of compiler generated functions. See dnaadjust.hpp
	dna_geoid_interpolation(const dna_geoid_interpolation&);
	dna_geoid_interpolation& operator=(const dna_geoid_interpolation&);	

public:
	// Interpolates AusGeoid values using bi linear interpolation 
	void FileTransformation(const char* fileIn, const char* fileOut, const int& method, 
		const int& intEllipsoidtoOrtho, const int& intDmsFlag, 
		bool exportDnaGeoidFile, const char* geoFile=NULL);

	void PopulateBinaryStationFile(const std::string& bstnFile, const int& method, 
		bool convertHeights, 
		bool exportDnaGeoidFile, const char* geoFile=NULL);

	// Interpolates AusGeoid values using bi linear interpolation and applies the N value to the supplied height
	void BiLinearTransformation(geoid_point* apPoint);

	// Interpolates AusGeoid values using bi cubic interpolation and applies the N value to the supplied height
	void BiCubicTransformation(geoid_point* apPoint);

	// Reads the gridfile header information into memory
	void CreateGridIndex(const char* fileName, const char* fileType);

	// Creates a NTv2 Grid file from conventional DAT file format
	void CreateNTv2File(const char* datFile, const n_file_par* grid);

	// Exports an Binary grid file to a Ascii format
	void ExportToAscii(const char* inputGrid, const char* gridType, const char* shiftType, const char* outputGrid, int* IO_Status);

	// Exports an Ascii grid file to a Binary format
	void ExportToBinary(const char* inputGrid, const char* gridType, const char* shiftType, const char* outputGrid, int* IO_Status);

	// Prints the gridfile header information to a file_par struct
	void ReportGridProperties(const char* fileName, const char* fileType, n_file_par* gridProperties);

	// Returns the version of the currently opened grid file
	void ReportGridVersion(char* version);

	// Returns the file progress
	inline int ReturnFileProgress() const { return (int)m_dPercentComplete; }

	// Returns the byte offset
	inline int GetByteOffset() const { return m_iBytesRead; }

	// Sets the byte offset
	inline void SetByteOffset() { m_iBytesRead = 0; }

	inline int PointsInterpolated() const { return m_pointsInterpolated; }
	inline int PointsNotInterpolated() const { return m_pointsNotInterpolated; }

	std::string ErrorString(const int& error, const std::string& data="");
	std::string ErrorCaption(const int& error);

	inline void SetInputCoordinates(const std::string& inputCoordinates) { m_inputCoordinates = inputCoordinates; } 

	std::string ReturnBadStationRecords();

private:

	// NTv2 File open/close
	void ClearGridFileMemory();
	int OpenGridFile(const char* filename, const char* filetype, 
		n_file_par* pGridfile, bool isTEMP=true);
	
	// NTv2 File creation
	void PrintGridHeaderInfoAscii(std::ofstream* f_out);
	void PrintGridHeaderInfoBinary(std::ofstream* f_out, n_file_par* pGridfile);
	void PrintSubGridHeaderInfoAscii(std::ofstream* f_out, n_gridfileindex* m_gfIndex);
	void PrintSubGridHeaderInfoBinary(std::ofstream* f_out, n_gridfileindex* m_gfIndex);
	void ScanDatFileValues(char* szLine, float* n_value, char* c_northsouth, int* lat_deg, int* lat_min, float* lat_sec, char* c_eastwest, int* lon_deg, int* lon_min, float* lon_sec, float* defl_meridian, float* defl_primev);
	void ScanNodeLocations(char* szLine, double* latitude, double* longitude);
	void ComputeLatLong(double* dlat_initial, const char& c_northsouth, const int& lat_deg, const int& lat_min, const double& lat_sec, double* dlon_initial, const char& c_eastwest, const int& lon_deg, const int& lon_min, const double& lon_sec);
	void WriteBinaryRecords(std::ofstream* f_out, float n_value, float defl_meridian, float defl_primev);

	// NTv2 interpolation
	void ApplyAusGeoidGrid(geoid_point* agCoord, const int& method);
	int FindSubGrid();
	bool IsWithinUpperLatitudeGridInterval(n_gridfileindex* gfIndex);
	bool IsWithinLowerLatitudeGridInterval(n_gridfileindex* gfIndex);
	bool IsWithinUpperLongitudeGridInterval(n_gridfileindex* gfIndex);
	bool IsWithinLowerLongitudeGridInterval(n_gridfileindex* gfIndex);
	int InterpolateNvalue_BiLinear(geoid_point* dInterpPoint);
	int InterpolateNvalue_BiCubic(geoid_point* dInterpPoint);
	int DetermineFileType(const char* cType);
	bool ReadAsciiShifts(geoid_values* pNShifts[], int iNodeIndex, long lNode);
	bool ReadBinaryShifts(geoid_values* pNShifts[], int iNodeIndex, long lNode);

	void LoadBinaryStationFile(const std::string& bstnfileName);
	void WriteBinaryStationFile(const std::string& bstnfileName);
	void PopulateStationRecords(const int& method, bool convertHeights);
	
	void WriteDNA1GeoidFile(const std::string& geofileName);
	void PrintDNA1GeoidRecord(std::ofstream& f_out, const std::string& station, const double& nValue, const double& meridianDef, const double& verticalDef);

	void ProcessCsvFile(std::ifstream* f_in, std::ofstream* f_out, std::ofstream* f_dnageo, const int& method, const int& intEllipsoidtoOrtho, const int& intDmsFlag);
	void ProcessDatFile(std::ifstream* f_in, std::ofstream* f_out, std::ofstream* f_dnageo, const int& method, const int& intEllipsoidtoOrtho, const int& intDmsFlag);

	//void OpenOutFile(std::ofstream& f_out, const std::string& fileOut, std::ios::open_mode mode, const std::string callingFunc);

	bool				m_isReading;
	bool				m_isWriting;

	// grid file to be held in memory
	std::ifstream		m_pGfileptr;
	n_file_par*			m_pGridfile;				// Holds the parameters used for NTv2 grid file

	double				m_dPercentComplete;			// percentage of bytes read from file
	int					m_iBytesRead;				// bytes read from file
	int					m_Grid_Success;				// Integer success value ( 0 = success; non zero = failure)
	
	binary_file_meta_t	bst_meta_;
	
	vstn_t				bstBinaryRecords_;
	v_stn_string		bstBadPoints_;

	UINT32				m_pointsInterpolated;
	UINT32				m_pointsNotInterpolated;

	bool				m_exportDNAGeoidFile;

	bool				m_fileMode;

	std::string				m_inputCoordinates;

	bool				m_isRadians;
	
	UINT32 				m_lineLength;
};

}	// namespace geoidinterpolation
}	// namespace dynadjust


#endif /* DNAGEOID_H_ */
