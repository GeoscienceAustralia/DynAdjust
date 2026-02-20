//============================================================================
// Name         : dnaversion.hpp
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
// Description  : DynAdjust version include file
//============================================================================

#pragma once
#ifndef DNAVERSION_HPP
#define DNAVERSION_HPP

#if defined(_MSC_VER)
#if defined(LIST_INCLUDES_ON_BUILD) 
#pragma message("  " __FILE__) 
#endif
#endif

#ifdef _MSC_VER
#ifdef _DEBUG
#define _MSDEBUG		// enables extra code to retrieve station names, variable
//#define _MSTRACE		// trace _DEBUG values to MSVC IDE debug output window
#endif // _DEBUG

// Force cross compatibility between GCC and MSVC
#ifdef NULL
#undef NULL					// #def'd in <stdio.h>
#define NULL 0
#endif

#endif	// _MSC_VER

/* definition to expand macro then apply to pragma message */
#define VALUE_TO_STRING(x) #x
#define VALUE(x) VALUE_TO_STRING(x)
#define VAR_NAME_VALUE(var) #var "=" VALUE(var)

// Force WINVER and _WIN32_WINNT to be a minimum of 0x0500, which 
// enables DynAdjust to run on Windows 2000. 
// Other "minimum" OS values include:
//  - 0x0501 for Windows XP
//  - 0x0502 for Windows Server 2003
//  - 0x0600 for Windows Vista
//  - 0x0601 for Windows 7
#ifdef WINVER
#undef WINVER
#endif

#define WINVER			0x0500

#ifdef _WIN32_WINNT
#undef _WIN32_WINNT
#endif

#define _WIN32_WINNT	0x0500

//
//_MSC_VER 
//	defines the compiler version for the versions supported by the current version of MKS Toolkit. Possible values include:
//
//		Microsoft Visual C++ 7.1	_MSC_VER = 1310
//		Microsoft Visual C++ 7.0	_MSC_VER = 1300
//		Microsoft Visual C++ 6.0	_MSC_VER = 1200
//		Microsoft Visual C++ 5.0	_MSC_VER = 1100
//
//_WIN32 
//	is defined for Win32 applications and is always defined as 1.
//
//_M_IX86 
//	defines the processor. Possible values include:
//
//		Blend		_M_IX86 = 500 
//		Pentium		_M_IX86 = 500
//		Pentium Pro 	_M_IX86 = 600
//		80386		_M_IX86 = 300
//		80486		_M_IX86 = 400
//
//

#if defined(_M_IX86)
#define __HARDWARE__ "Win 32"
#elif defined(_M_X64)
#define __HARDWARE__ "Win 64"
#elif defined(__linux)
#define __HARDWARE__ "Linux"
#elif defined(sun)
#define __HARDWARE__ "Sun"
#elif defined(__unix__)
#define __HARDWARE__ "UNIX"
#elif defined(__APPLE__)
#define __HARDWARE__ "Apple"

#endif

const char* const __geoidint_app_name__ = "GeoidInt";

// Define Release or Debug mode
// 
// The documentation for the predefined macros says:
//  _WIN32: Defined for applications for Win32 and Win64. Always defined.
//  _WIN64: Defined for applications for Win64.
//
// So not only should _WIN32 always be defined, it does not cause any 
// problem in 64-bit applications.

// http://stackoverflow.com/questions/6679396/should-i-define-both-win32-and-win64-in-64bit-build

const char* const __the_app_name__ = "DynAdjust";

const char* const __metadata_app_name__ = "metadata";
const char* const __dynadjust_app_name__ = "dynadjust";

#if defined(USE_MKL)
    #define __LINEAR_ALGEBRA_LIB__ "Intel MKL"
#elif defined(__APPLE__)
    #define __LINEAR_ALGEBRA_LIB__ "Apple Accelerate Framework"
#else
    #define __LINEAR_ALGEBRA_LIB__ "OpenBLAS"
#endif

#if defined(_WIN32) || defined(_Wp64)
  #ifdef _MSDEBUG
    #if defined(_M_IX86)
      #define __BINARY_BUILDTYPE__ "Debug (32-bit)"
    #elif defined(_M_X64)
      #define __BINARY_BUILDTYPE__ "Debug (64-bit)"
    #endif
  #else
    #if defined(_M_IX86)
      #define __BINARY_BUILDTYPE__ "Release (32-bit)"
    #elif defined(_M_X64)
      #define __BINARY_BUILDTYPE__ "Release (64-bit)"
    #endif
  #endif

const char* const __import_app_name__ = "import";
const char* const __reftran_app_name__ = "reftran";
const char* const __geoid_app_name__ = "geoid";
const char* const __segment_app_name__ = "segment";
const char* const __adjust_app_name__ = "adjust";
const char* const __plot_app_name__ = "plot";

const char* const __import_dll_name__ = "dnaInterop.dll";
const char* const __reftran_dll_name__ = "dnaRefTran.dll";
const char* const __geoid_dll_name__ = "dnaGeoid.dll";
const char* const __segment_dll_name__ = "dnaSegment.dll";
const char* const __adjust_dll_name__ = "dnaAdjust.dll";
const char* const __plot_dll_name__ = "dnaPlot.dll";

#elif defined(__linux) || defined(sun) || defined(__unix__) || defined(__APPLE__)
  #ifdef NDEBUG
    #define __BINARY_BUILDTYPE__ "Release"
  #else
    #define __BINARY_BUILDTYPE__ "Debug"
#endif

const char* const __import_app_name__ = "dnaimport";
const char* const __reftran_app_name__ = "dnareftran";
const char* const __geoid_app_name__ = "dnageoid";
const char* const __segment_app_name__ = "dnasegment";
const char* const __adjust_app_name__ = "dnaadjust";
const char* const __plot_app_name__ = "dnaplot";

const char* const __import_dll_name__ = "libdnaimport.so";
const char* const __reftran_dll_name__ = "libdnareftran.so";
const char* const __geoid_dll_name__ = "libdnageoid.so";
const char* const __segment_dll_name__ = "libdnasegment.so";
const char* const __adjust_dll_name__ = "libdnaadjust.so";
const char* const __plot_dll_name__ = "libdnaplot.so";

#endif

#define __CONTACT_EMAIL__   "geodesy@ga.gov.au"
#define __COPYRIGHT_YEAR__  "2025"
#define __COPYRIGHT_OWNER__ "Geoscience Australia"
#define __COPYRIGHT_MSG__   "This software is released under the Apache License."

// Version format is AABBCCDD where
//
// AA - major version number
// BB - minor version number
// CC - bugfix version number
// DD - alpha / beta (DD + 50) version number
//
// When DD is not 00, 1 is subtracted from AABBCC. For example:
//
// Version     AABBCCDD
// 2.0.0       02000000
// 2.1.0       02010000
// 2.1.1       02010100
// 2.2.0.a1    02020001		Alpha release 
// 2.2.0.b2    02020002		Beta release
// 2.2.0.rc2   02020003		Release candidate 1 (unless bug notification is received, this is taken to be the final release)
#define __BINARY_VERSION__ "1.3.0"
#define __SHORT_VERSION__ "10300"	// used to record DynAdjust version in binary file header

// define executable name
#define __GLOBAL_BINARY_NAME__ __dynadjust_app_name__
#define __GLOBAL_BINARY_DESC__ "Geodetic network adjustment software"

#ifdef BUILD_DNAIMPORT_DLL
#define __BINARY_NAME__ __import_dll_name__
#define __BINARY_DESC__ "File import and format exchange library"
#endif

// define BUILD_DNASEGMENT_DLL when building libdnasegment.dll
#ifdef BUILD_DNASEGMENT_DLL
#define __BINARY_NAME__ __segment_dll_name__
#define __BINARY_DESC__ "Automated network segmentation library"
#endif

// define BUILD_DNAADJUST_DLL when building libdnasegment.dll
#ifdef BUILD_DNAADJUST_DLL
#define __BINARY_NAME__ __adjust_dll_name__
#define __BINARY_DESC__ "Geodetic network adjustment library"
#endif

// define BUILD_DNAPLOT_DLL when building libdnaplot.dll
#ifdef BUILD_DNAPLOT_DLL
#define __BINARY_NAME__ __plot_dll_name__
#define __BINARY_DESC__ "Geodetic network plotting library"
#endif

// define BUILD_DNAGEOID_DLL when building libdnageoid.dll
#ifdef BUILD_DNAGEOID_DLL
#define __BINARY_NAME__ __geoid_dll_name__
#define __BINARY_DESC__ "Geoid grid file interpolation library"
#endif

// define BUILD_DNAREFTRAN_DLL when building libdnareftran.dll
#ifdef BUILD_DNAREFTRAN_DLL
#define __BINARY_NAME__ __reftran_dll_name__
#define __BINARY_DESC__ "Reference frame transformation library"
#endif

// define executable name
#ifdef BUILD_DYNADJUST_EXE
#define __BINARY_NAME__ __dynadjust_app_name__
#define __BINARY_DESC__ "Geodetic network adjustment software"
#endif

// define executable name
#ifdef BUILD_PLOT_EXE
#define __BINARY_NAME__ __plot_app_name__
#define __BINARY_DESC__ "Geodetic network plotting software"
#endif

// define executable name
#ifdef BUILD_IMPORT_EXE
#define __BINARY_NAME__ __import_app_name__
#define __BINARY_DESC__ "File import and format exchange software"
#endif

// define executable name
#ifdef BUILD_SEGMENT_EXE
#define __BINARY_NAME__ __segment_app_name__
#define __BINARY_DESC__ "Automated network segmentation software"
#endif

// define executable name
#ifdef BUILD_ADJUST_EXE
#define __BINARY_NAME__ __adjust_app_name__
#define __BINARY_DESC__ "Geodetic network adjustment software"
#endif

// define executable name
#ifdef BUILD_GEOID_EXE
#define __BINARY_NAME__ __geoid_app_name__
#define __BINARY_DESC__ "Geoid grid file interpolation software"
#endif

// define executable name
#ifdef BUILD_GEOIDINT_EXE
#define __BINARY_NAME__ __geoidint_app_name__
#define __BINARY_DESC__ "Geoid Interpolation software"
#endif

// define executable name
#ifdef BUILD_REFTRAN_EXE
#define __BINARY_NAME__ __reftran_app_name__
#define __BINARY_DESC__ "Reference frame transformation software"
#endif

// define executable name
#ifdef BUILD_METADATA_EXE
#define __BINARY_NAME__ __metadata_app_name__
#define __BINARY_DESC__ "Metadata manipulation software"
#endif


#if defined(_MSC_VER) && (_MSC_VER > 1100)
#define _MS_COMPILER_ _MSC_VER
//#include <boost/preprocessor/stringize.hpp>
//#pragma message("_MSC_VER=" BOOST_PP_STRINGIZE(_MSC_VER))
#endif


#if defined(__clang__)                                   // Clang compiler
  #define __COMPILER__ "Clang"
  #define __COMPILER_VERSION__ __VERSION__
#elif defined(__GNUC__) || defined(__GNUG__)				// GNU GCC
  #define __COMPILER__ "GNU GCC"
  #define __COMPILER_VERSION__ __VERSION__
#elif defined(__SUNPRO_CC)								// Oracle Solaris
  #define __COMPILER__ "Solaris"
  #define __COMPILER_VERSION__ __SUNPRO_CC
#elif defined(__ICC) || defined(__INTEL_COMPILER)		// Intel compiler
  #define __COMPILER__ "Intel ICC"
  #define __COMPILER_VERSION__ std::to_string(__INTEL_COMPILER)
#elif defined(__HP_aCC)								// Oracle Solaris
  #define __COMPILER__ "HP C++"
  #define __COMPILER_VERSION__ __HP_aCC
#elif defined(_MSC_VER)
  // https://learn.microsoft.com/en-us/cpp/overview/compiler-versions?view=msvc-170#version-macros
  #define __COMPILER__ "MSVC++"
  #if (_MSC_VER == 1100)	
  #define __COMPILER_VERSION__ "5.0"
  #elif (_MSC_VER == 1200)
  #define __COMPILER_VERSION__ "6.0"
  #elif (_MSC_VER == 1300)
  #define __COMPILER_VERSION__ "7.0"
  #elif (_MSC_VER == 1310)
  #define __COMPILER_VERSION__ "7.1, VS2003"
  #elif (_MSC_VER == 1400)
  #define __COMPILER_VERSION__ "8.0, VS2005"
  #elif (_MSC_VER == 1500)
  #define __COMPILER_VERSION__ "9.0, VS2008"
  #elif (_MSC_VER == 1600)
  #define __COMPILER_VERSION__ "10.0, VS2010"
  #elif (_MSC_VER == 1700)
  #define __COMPILER_VERSION__ "11.0, VS2012"
  #elif (_MSC_VER == 1800)
  #define __COMPILER_VERSION__ "12.0, VS2013"
  #elif (_MSC_VER == 1900)
  #define __COMPILER_VERSION__ "14.0, VS2015"
  #elif (_MSC_VER == 1910)
  #define __COMPILER_VERSION__ "14.1, VS2017"
  #elif (_MSC_VER == 1911)
  #define __COMPILER_VERSION__ "15.3, VS2017"
  #elif (_MSC_VER == 1912)
  #define __COMPILER_VERSION__ "15.5, VS2017"
  #elif (_MSC_VER == 1913)
  #define __COMPILER_VERSION__ "15.6, VS2017"
  #elif (_MSC_VER == 1914)
  #define __COMPILER_VERSION__ "15.7, VS2017"
  #elif (_MSC_VER == 1915)
  #define __COMPILER_VERSION__ "15.8, VS2017"
  #elif (_MSC_VER == 1916)
  #define __COMPILER_VERSION__ "15.9, VS2017"
  #elif (_MSC_VER == 1920)
  #define __COMPILER_VERSION__ "16.0, VS2019"
  #elif (_MSC_VER == 1921)
  #define __COMPILER_VERSION__ "16.1, VS2019"
  #elif (_MSC_VER == 1922)
  #define __COMPILER_VERSION__ "16.2, VS2019"
  #elif (_MSC_VER == 1923)
  #define __COMPILER_VERSION__ "16.3, VS2019"
  #elif (_MSC_VER == 1924)
  #define __COMPILER_VERSION__ "16.4, VS2019"
  #elif (_MSC_VER == 1925)
  #define __COMPILER_VERSION__ "16.5, VS2019"
  #elif (_MSC_VER == 1926)
  #define __COMPILER_VERSION__ "16.6, VS2019"
  #elif (_MSC_VER == 1927)
  #define __COMPILER_VERSION__ "16.7, VS2019"
  #elif (_MSC_VER == 1928)
  #define __COMPILER_VERSION__ "16.8/9, VS2019"
  #elif (_MSC_VER == 1929)
  #define __COMPILER_VERSION__ "16.10/11, VS2019"
  #elif (_MSC_VER == 1930)
  #define __COMPILER_VERSION__ "17.0, VS2022"
  #elif (_MSC_VER == 1931)
  #define __COMPILER_VERSION__ "17.1, VS2022"
  #elif (_MSC_VER == 1932)
  #define __COMPILER_VERSION__ "17.2, VS2022"
  #elif (_MSC_VER == 1933)
  #define __COMPILER_VERSION__ "17.3, VS2022"
  #elif (_MSC_VER == 1934)
  #define __COMPILER_VERSION__ "17.4, VS2022"
  #elif (_MSC_VER == 1935)
  #define __COMPILER_VERSION__ "17.5, VS2022"
  #elif (_MSC_VER == 1936)
  #define __COMPILER_VERSION__ "17.6, VS2022"
  #elif (_MSC_VER == 1937)
  #define __COMPILER_VERSION__ "17.7, VS2022"
  #elif (_MSC_VER == 1938)
  #define __COMPILER_VERSION__ "17.8, VS2022"
  #elif (_MSC_VER == 1939)
  #define __COMPILER_VERSION__ "17.9, VS2022"
  #elif (_MSC_VER == 1940)
  #define __COMPILER_VERSION__ "17.10, VS2022"
  #elif (_MSC_VER == 1941)
  #define __COMPILER_VERSION__ "17.11, VS2022"
  #elif (_MSC_VER == 1942)
  #define __COMPILER_VERSION__ "17.12, VS2022"
  #elif (_MSC_VER == 1943)
  #define __COMPILER_VERSION__ "17.13, VS2022"
  #elif (_MSC_VER == 1944)
  #define __COMPILER_VERSION__ "17.14, VS2022"
  #elif (_MSC_VER == 1950)
  #define __COMPILER_VERSION__ "18.0, VS2026"
  
  #elif (_MSC_VER > 1950)
    #define __COMPILER_VERSION__ _MSC_VER
    #pragma message("Warning: A later than recognised VC++ compiler version is being used: " VAR_NAME_VALUE(_MSC_VER))
    #pragma message("         Time to update " VALUE(__FILE__))
  #endif
#else
  #define __COMPILER__ "Unknown"
#endif

#endif  // DNAVERSION_HPP
