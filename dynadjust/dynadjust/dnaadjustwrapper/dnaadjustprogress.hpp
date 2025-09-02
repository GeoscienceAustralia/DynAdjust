//============================================================================
// Name         : dnaadjustprogress.hpp
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

#ifndef DNAADJUSTPROGRESS_H_
#define DNAADJUSTPROGRESS_H_

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
#include <thread>
#include <chrono>
/// \endcond

#include <include/config/dnaconsts-interface.hpp>

#include <dynadjust/dnaadjust/dnaadjust.hpp>

using namespace dynadjust::networkadjust;
using namespace dynadjust::exception;

class dna_adjust_thread
{
public:
	dna_adjust_thread(dna_adjust* dnaAdj, project_settings* p,
		_ADJUST_STATUS_* adjustStatus);
	void operator()();
	
private:
	bool prepareAdjustment();
	bool processAdjustment();

	void handlePrepareAdjustError(const std::string& error_msg);
	void handleProcessAdjustError(const std::string& error_msg);
	void printErrorMsg(const std::string& error_msg);
	void coutMessage(std::stringstream& message);
	void coutMessage(const std::string& message);

	dna_adjust*		_dnaAdj;
	project_settings*	_p;
	_ADJUST_STATUS_*	_adjustStatus;
};

class dna_adjust_progress_thread
{
public:
	dna_adjust_progress_thread(dna_adjust* dnaAdj, project_settings* p);
	void operator()();

private:
	void prepareAdjustment();
	void processAdjustment();

	void coutMessage(const std::string& message);

	dna_adjust*		_dnaAdj;
	project_settings*	_p;
};

#endif //DNAADJUSTPROGRESS_H_
