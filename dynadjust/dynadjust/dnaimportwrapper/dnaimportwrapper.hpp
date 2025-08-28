//============================================================================
// Name         : dnaimportwrapper.hpp
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
// Description  : DynAdjust Interoperability library Executable
//============================================================================

#ifndef DNAIMPORTWRAPPER_H_
#define DNAIMPORTWRAPPER_H_

#if defined(_MSC_VER)
	#if defined(LIST_INCLUDES_ON_BUILD) 
		#pragma message("  " __FILE__) 
	#endif
#endif

#if defined(_WIN32) || defined(__WIN32__)
	#include <io.h>
#elif defined(__linux) || defined(sun) || defined(__unix__) || defined(__APPLE__)
	#include <unistd.h>
#endif

/// \cond
#include <iostream>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <sstream>
#include <string>
#include <time.h>
#include <set>
#include <mutex>
#include <thread>
#include <chrono>

#include <memory>
#include <filesystem>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <include/config/dnaversion.hpp>
#include <include/config/dnaconsts.hpp>
#include <include/config/dnaconsts-interface.hpp>
#include <include/config/dnaoptions-interface.hpp>
#include <include/config/dnaprojectfile.hpp>

#include <include/functions/dnaiostreamfuncs.hpp>
#include <include/functions/dnafilepathfuncs.hpp>
#include <include/functions/dnastrmanipfuncs.hpp>
#include <include/functions/dnatimer.hpp>

#include <dynadjust/dnaimport/dnainterop.hpp>

using namespace dynadjust::dynamlinterop;
using namespace dynadjust::exception;
using namespace dynadjust::datum_parameters;

extern bool running;
extern std::mutex cout_mutex;

using dynadjust::cpu_timer;

class dna_import_thread
{
public:
	dna_import_thread(dna_import* dnaParse, project_settings* p, const std::string& filename,
		vdnaStnPtr* vStations, PUINT32 stnCount, vdnaMsrPtr* vMeasurements, PUINT32 msrCount,
		PUINT32 clusterID, input_file_meta_t* input_file_meta, bool firstFile, std::string* status_msg,
		boost::posix_time::milliseconds* ms)
		: _dnaParse(dnaParse), _p(p), _filename(filename)
		, _vStations(vStations), _stnCount(stnCount), _vMeasurements(vMeasurements), _msrCount(msrCount)
		, _clusterID(clusterID), _input_file_meta(input_file_meta), _firstFile(firstFile)
		, _status_msg(status_msg), _ms(ms) {};
	void operator()()
	{
		cpu_timer time;	// constructor of cpu_timer calls start()
		try {
			_dnaParse->ParseInputFile(_filename, 
				_vStations, _stnCount, 
				_vMeasurements, _msrCount, 
				_clusterID, _input_file_meta, _firstFile,
				_status_msg, _p);
			*_ms = boost::posix_time::milliseconds(std::chrono::duration_cast<std::chrono::milliseconds>(time.elapsed().wall).count());
		} 
		catch (const XMLInteropException& e) {
			running = false;
			std::this_thread::sleep_for(std::chrono::milliseconds(50));
			std::stringstream err_msg;
			cout_mutex.lock();
			err_msg << std::endl << "- Error: " << e.what() << std::endl;
			std::cout << err_msg.str();
			*_status_msg = err_msg.str();
			cout_mutex.unlock();
			return;
		}
		running = false;		
	}
	inline void SetFile(const std::string& file) { _filename = file; }

private:
	dna_import*			_dnaParse;
	project_settings*	_p;
	std::string				_filename;
	vdnaStnPtr*			_vStations;
	PUINT32				_stnCount;
	vdnaMsrPtr*			_vMeasurements;
	PUINT32				_msrCount;
	PUINT32				_clusterID;
	input_file_meta_t*	_input_file_meta;
	bool				_firstFile;
	std::string*				_status_msg;
	boost::posix_time::milliseconds*		_ms;
};
	
class dna_import_progress_thread
{
public:
	dna_import_progress_thread(dna_import* dnaParse, project_settings* p)
		: _dnaParse(dnaParse), _p(p) {};
	void operator()()
	{
		double percentComplete(0.);
		std::ostringstream ss;

		int is_terminal(isatty(fileno(stdout)));

		percentComplete = 0.;
		while (running)
		{
			if (percentComplete > 100)
				percentComplete = 0.;
			else if (percentComplete < 0.)
				percentComplete = 0.;
			if (is_terminal && !_p->g.quiet)
			{
				ss.str("");
				ss << std::setw(3) << std::fixed << std::setprecision(0) << std::right << percentComplete << "%";
				cout_mutex.lock();
				std::cout << PROGRESS_BACKSPACE_04 << std::setw(PROGRESS_PERCENT_04) << ss.str();
				std::cout.flush();
				cout_mutex.unlock();
			}
			std::this_thread::sleep_for(std::chrono::milliseconds(10));
			percentComplete = _dnaParse->GetProgress();
		}
	}
private:
	dna_import*		_dnaParse;
	project_settings*	_p;
};

#endif
