//============================================================================
// Name         : dnastringfuncs.hpp
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
// Description  : Common String Functions
//============================================================================

#ifndef DNASTRINGFUNCS_HPP_
#define DNASTRINGFUNCS_HPP_

#if defined(_MSC_VER)
	#if defined(LIST_INCLUDES_ON_BUILD) 
		#pragma message("  " __FILE__) 
	#endif
#endif

/// \cond
#include <stdarg.h>
#include <sstream>
#include <iomanip>
/// \endcond

#include <include/config/dnaconsts.hpp>
#include <include/config/dnaversion.hpp>

// Returns the number of fields with valid data
int GetFields(char *line, char delim, bool multiple_delim_as_one, const char *fmt, ...);
void fileproc_help_header(std::string* msg);
void dynaml_header(std::ostream& os, const std::string& fileType, const std::string& referenceFrame, const std::string& epoch);
void dynaml_footer(std::ostream& os);
void dynaml_comment(std::ostream& os, const std::string& comment);
void dna_header(std::ostream& os, const std::string& fileVersion, const std::string& fileType, const std::string& reference_frame, const std::string& epoch_version, const size_t& count);
void dna_comment(std::ostream& os, const std::string& comment);
void dnaproj_header(std::ostream& os, const std::string& comment);
void dnaproj_comment(std::ostream& os, const std::string& comment);

std::string snx_softwarehardware_text();

#endif // DNASTRINGFUNCS_HPP_