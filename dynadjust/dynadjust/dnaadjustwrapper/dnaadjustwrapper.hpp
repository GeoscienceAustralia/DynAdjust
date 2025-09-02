//============================================================================
// Name         : dnaadjustwrapper.hpp
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
// Description  : DynAdjust Adjustment library Executable
//============================================================================

#ifndef DNAADJUSTWRAPPER_H_
#define DNAADJUSTWRAPPER_H_

#if defined(_MSC_VER)
	#if defined(LIST_INCLUDES_ON_BUILD) 
		#pragma message("  " __FILE__) 
	#endif
#endif

/// \cond
#include <iostream>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <sstream>
#include <string>
#include <time.h>
#include <memory>
#include <thread>
#include <filesystem>

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
/// \endcond

#include <include/functions/dnastrutils.hpp>
#include <include/functions/dnafilepathfuncs.hpp>
#include <include/functions/dnastrmanipfuncs.hpp>

#include <include/config/dnaprojectfile.hpp>
#include <include/config/dnaoptions-interface.hpp>

#include <dynadjust/dnaadjust/dnaadjust.hpp>
#include <dynadjust/dnaadjust/dnaadjust_printer.hpp>
#include <dynadjust/dnaadjustwrapper/dnaadjustprogress.hpp>
#include <dynadjust/dnaadjustwrapper/threading_init.hpp>

using namespace dynadjust::networkadjust;
using namespace dynadjust::exception;
using namespace dynadjust::iostreams;

#endif
