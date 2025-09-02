// precompile.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#pragma once

// See https://github.com/boostorg/process/issues/161
#define _WIN32_WINNT 0x0501

// Support MKL inverse oly if compiler is Intel
#if defined(__ICC) || defined(__INTEL_COMPILER) // Intel compiler
#include <mkl.h>
#endif

#include <include/config/dnaconsts.hpp>
#include <include/config/dnaversion.hpp>
#include <include/config/dnatypes.hpp>
#include <include/config/dnatypes-gui.hpp>
#include <include/config/dnaexports.hpp>

#include <include/exception/dnaexception.hpp>

#include <include/io/bst_file.hpp>
#include <include/io/bms_file.hpp>
#include <include/io/map_file.hpp>
#include <include/io/seg_file.hpp>
#include <include/io/dnaiopdf.hpp>

#include <include/functions/dnatemplatefuncs.hpp>
#include <include/functions/dnatemplatestnmsrfuncs.hpp>
#include <include/functions/dnaiostreamfuncs.hpp>
#include <include/functions/dnastringfuncs.hpp>
#include <include/functions/dnaintegermanipfuncs.hpp>
#include <include/functions/dnafilepathfuncs.hpp>
#include <include/functions/dnatemplategeodesyfuncs.hpp>
#include <include/functions/dnatemplatematrixfuncs.hpp>

#include <include/parameters/dnaepsg.hpp>
#include <include/parameters/dnadatum.hpp>
#include <include/parameters/dnaprojection.hpp>

#include <include/measurement_types/dnadirection.hpp>
#include <include/measurement_types/dnagpsbaseline.hpp>
#include <include/measurement_types/dnagpspoint.hpp>

#include <include/math/dnamatrix_contiguous.hpp>

#include <include/memory/dnamemory_handler.hpp>

#include <dynadjust/dnaplot/dnaplot.hpp>
