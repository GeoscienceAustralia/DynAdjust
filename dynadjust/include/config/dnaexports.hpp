//============================================================================
// Name         : dnaexports.hpp
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
// Description  : DynAdjust exported symbols include file
//============================================================================

#ifndef DNAEXPORTS_H_
#define DNAEXPORTS_H_

#if defined(_MSC_VER)
	#if defined(LIST_INCLUDES_ON_BUILD) 
		#pragma message("  " __FILE__) 
	#endif
#endif

// define BUILD_DNAIMPORT_DLL when building libdnaimport.dll
#if defined(_MSC_VER) && (!defined(BUILD_STATIC))
    #ifdef BUILD_DNAIMPORT_DLL
    /* DLL export */
    #define DNAIMPORT_API __declspec(dllexport)
    #else
    /* EXE import */
    #define DNAIMPORT_API __declspec(dllimport)
    #endif
#else
    /* static library or non-Windows */
    #define DNAIMPORT_API
#endif
	
// define BUILD_DNASEGMENT_DLL when building libdnasegment.dll
#if defined(_MSC_VER) && (!defined(BUILD_STATIC))
    #ifdef BUILD_DNASEGMENT_DLL
    /* DLL export */
    #define DNASEGMENT_API __declspec(dllexport)
    #else
    /* EXE import */
    #define DNASEGMENT_API __declspec(dllimport)
    #endif
#else
    /* static library or non-Windows */
    #define DNASEGMENT_API
#endif
	
// define BUILD_DNAADJUST_DLL when building libdnaadjust.dll
#if defined(_MSC_VER) && (!defined(BUILD_STATIC))
    #ifdef BUILD_DNAADJUST_DLL
    /* DLL export */
    #define DNAADJUST_API __declspec(dllexport)
    #else
    /* EXE import */
    #define DNAADJUST_API __declspec(dllimport)
    #endif
#else
    /* static library or non-Windows */
    #define DNAADJUST_API
#endif
	
// define BUILD_DNAMATH_DLL when building libdnamath.dll
#if defined(_MSC_VER) && (!defined(BUILD_STATIC))
    #ifdef BUILD_DNAMATH_DLL
    /* DLL export */
    #define DNAMATHCOMP_API __declspec(dllexport)
    #else
    /* EXE import */
    #define DNAMATHCOMP_API __declspec(dllimport)
    #endif
#else
    /* static library or non-Windows */
    #define DNAMATHCOMP_API
#endif
	
// define BUILD_EXPORTDNATYPES when building dlls that use Dna Measurement types
#if defined(_MSC_VER) && (!defined(BUILD_STATIC))
    #ifdef BUILD_EXPORTDNATYPES
    /* DLL export */
    #define DNATYPE_API __declspec(dllexport)
    #else
    /* EXE import */
    #define DNATYPE_API __declspec(dllimport)
    #endif
#else
    /* static library or non-Windows */
    #define DNATYPE_API
#endif

// define BUILD_DNAPLOT_DLL when building libdnaplot.dll
#if defined(_MSC_VER) && (!defined(BUILD_STATIC))
    #ifdef BUILD_DNAPLOT_DLL
    /* DLL export */
    #define DNAPLOT_API __declspec(dllexport)
    #else
    /* EXE import */
    #define DNAPLOT_API __declspec(dllimport)
    #endif
#else
    /* static library or non-Windows */
    #define DNAPLOT_API
#endif

// define BUILD_DNAGEOID_DLL when building libdnageoid.dll
#if defined(_MSC_VER) && (!defined(BUILD_STATIC))
    #ifdef BUILD_DNAGEOID_DLL
    /* DLL export */
    #define DNAGEOID_API __declspec(dllexport)
    #else
    /* EXE import */
    #define DNAGEOID_API __declspec(dllimport)
    #endif
#else
    /* static library or non-Windows */
    #define DNAGEOID_API
#endif
	
// define BUILD_DNAREFTRAN_DLL when building libdnareftran.dll
#if defined(_MSC_VER) && (!defined(BUILD_STATIC))
    #ifdef BUILD_DNAREFTRAN_DLL
    /* DLL export */
    #define DNAREFTRAN_API __declspec(dllexport)
    #else
    /* EXE import */
    #define DNAREFTRAN_API __declspec(dllimport)
    #endif
#else
    /* static library or non-Windows */
    #define DNAREFTRAN_API
#endif
	
#endif // DNAEXPORTS_H_
