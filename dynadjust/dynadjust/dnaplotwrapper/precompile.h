// precompile.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#pragma once

// See https://github.com/boostorg/process/issues/161
#define _WIN32_WINNT 0x0501

#include <include/config/dnaconsts.hpp>
#include <include/config/dnaversion.hpp>
#include <include/config/dnatypes.hpp>
#include <include/config/dnatypes-gui.hpp>
#include <include/config/dnaexports.hpp>
#include <include/config/dnaoptions-interface.hpp>
#include <include/config/dnaprojectfile.hpp>

#include <include/exception/dnaexception.hpp>

#include <include/functions/dnastringfuncs.hpp>
#include <include/functions/dnaprocessfuncs.hpp>
#include <include/functions/dnafilepathfuncs.hpp>

#include <dynadjust/dnaplot/dnaplot.hpp>

#include <dynadjust/dnaplotwrapper/dnaplotwrapper.hpp>
