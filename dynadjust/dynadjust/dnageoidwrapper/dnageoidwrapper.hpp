//============================================================================
// Name         : dnageoidwrapper.hpp
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
// Description  : AusGeoid Grid File (NTv2) Interpolation library Executable
//============================================================================

#ifndef DNAGEOIDWRAPPER_H_
#define DNAGEOIDWRAPPER_H_

#if defined(_MSC_VER)
	#if defined(LIST_INCLUDES_ON_BUILD) 
		#pragma message("  " __FILE__) 
	#endif
#endif

/// \cond
#include <iomanip>

#include <iostream>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <sstream>
#include <string>
#include <time.h>
#include <string.h>
#include <set>


#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <filesystem>
/// \endcond
#include <include/config/dnaexports.hpp>
#include <include/config/dnaversion.hpp>
#include <include/config/dnaconsts.hpp>
#include <include/config/dnaprojectfile.hpp>
#include <include/config/dnaoptions-interface.hpp>
#include <include/config/dnatypes-gui.hpp>
#include <include/config/dnaversion-stream.hpp>

#include <include/functions/dnafilepathfuncs.hpp>
#include <include/functions/dnatemplatefuncs.hpp>
#include <include/functions/dnastringfuncs.hpp>
#include <include/functions/dnafilepathfuncs.hpp>
#include <include/functions/dnatemplatedatetimefuncs.hpp>
#include <include/functions/dnastrmanipfuncs.hpp>
#include <include/functions/dnatemplatefuncs.hpp>
#include <include/functions/dnatimer.hpp>

#include <dynadjust/dnageoid/dnageoid.hpp>

using namespace dynadjust::geoidinterpolation;

bool CreateNTv2Grid(dna_geoid_interpolation* g, const char* dat_gridfilePath, const n_file_par* grid);
bool createGridIndex(dna_geoid_interpolation* g, const char* gridfilePath, const char* gridfileType, const int& quiet);
bool reportGridProperties(dna_geoid_interpolation* g, const char* gridfilePath, const char* gridfileType);
bool InterpolateGridPoint(dna_geoid_interpolation* g, const char* gridfilePath, geoid_point* apInterpolant, 
	const int& method, const int& coordinate_format, const std::string& inputLatitude, const std::string& inputLongitude);
bool InterpolateGridPointFile(dna_geoid_interpolation* g, const char* inputfilePath, 
	const int& method, const int EllipsoidtoOrtho, const int& coordinate_format, 
	bool exportDnaGeoidFile, const char* dnageofilePath, std::string& outputfilePath);
bool InterpolateGridBinaryStationFile(dna_geoid_interpolation* g, const std::string& bstnfilePath,
	const int& method, bool convertHeights, bool exportDnaGeoidFile, const char* dnageofilePath);
int ParseCommandLineOptions(const int& argc, char* argv[], const boost::program_options::variables_map& vm, project_settings& p);

#endif
