// precompile.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#pragma once

// prevent conflict with std::min(...) std::max(...)
#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#endif

/// \cond
#include <windows.h>

#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/algorithm/string/predicate.hpp>

/// \endcond

#include <include/config/dnaversion.hpp>
#include <include/config/dnaconsts.hpp>
#include <include/config/dnaconsts-interface.hpp>
#include <include/config/dnaoptions.hpp>
#include <include/config/dnaoptions-interface.hpp>

#include <include/exception/dnaexception.hpp>

#include <include/functions/dnastringfuncs.hpp>
#include <include/functions/dnaiostreamfuncs.hpp>
#include <include/functions/dnafilepathfuncs.hpp>
#include <include/functions/dnastrmanipfuncs.hpp>
#include <include/functions/dnatemplatecalcfuncs.hpp>
#include <include/functions/dnatemplatedatetimefuncs.hpp>

#include <include/parameters/dnadatum.hpp>
#include <include/parameters/dnadatumprojectionparam.hpp>
#include <include/parameters/dnaellipsoid.hpp>