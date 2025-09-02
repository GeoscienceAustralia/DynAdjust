//============================================================================
// Name         : dynadjust.cpp
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
// Description  : DynAdjust implementation
//============================================================================

/// \cond
#include <iostream>
#include <fstream>
#include <algorithm>
#include <iomanip>
#include <sstream>
#include <string>
#include <time.h>
#include <mutex>
#include <thread>
#include <chrono>
#include <memory>
#include <filesystem>

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
/// \endcond

std::mutex cout_mutex;

#include <include/config/dnaconsts.hpp>
#include <include/config/dnaversion.hpp>
#include <include/config/dnaoptions.hpp>
#include <include/config/dnaoptions-interface.hpp>
#include <include/config/dnaprojectfile.hpp>

#include <include/functions/dnastringfuncs.hpp>
#include <include/functions/dnaprocessfuncs.hpp>
#include <include/functions/dnaiostreamfuncs.hpp>
#include <include/functions/dnafilepathfuncs.hpp>
#include <include/functions/dnatemplatedatetimefuncs.hpp>
#include <include/functions/dnastrmanipfuncs.hpp>
#include <include/functions/dnastrutils.hpp>

using namespace dynadjust;

int CloseLogandReturn(std::ofstream& os, const int status, const std::string& output_file="")
{
	if (status == EXIT_FAILURE)
	{
		os << std::setw(25) << std::left << formattedDateTimeString<std::string>() << 
			"Failed. " << std::endl << std::endl << 
			"+ DynAdjust ended prematurely. ";

		if (!output_file.empty())
		{
			std::ifstream f(output_file);
			os << " Contents of last output file:" << std::endl << std::endl << f.rdbuf();
		}
		else
			 os << std::endl;
	}
	else
		os << std::endl << std::endl << "+ DynAdjust finished successfully." << std::endl;

	os.close();
	return status;
}

void PrintAppStartTimeMessage(std::ofstream& os, std::string&& app)
{
	// App start time
	os << std::setw(15) << std::left << std::string("+ " + app) << std::setw(25) << 
		formattedDateTimeString<std::string>();
}

void PrintSuccessStatusMessage(std::ofstream& os)
{
	// end time
	os << std::setw(25) << std::left << formattedDateTimeString<std::string>() << 
		"Ended successfully." << std::endl;
}

int main(int argc, char* argv[])
{
	std::stringstream ss;

	std::string cmd_line_banner;
	fileproc_help_header(&cmd_line_banner);

	project_settings p;
	boost::program_options::variables_map vm;
	boost::program_options::positional_options_description positional_options;

	boost::program_options::options_description standard_options("+ " + std::string(ALL_MODULE_STDOPT), PROGRAM_OPTIONS_LINE_LENGTH);
	boost::program_options::options_description generic_options("+ " + std::string(ALL_MODULE_GENERIC), PROGRAM_OPTIONS_LINE_LENGTH);

	std::string cmd_line_usage("+ ");
	cmd_line_usage.append(__BINARY_NAME__).append(" usage:  ").append(__BINARY_NAME__).append(" ").append(" [options]");
	boost::program_options::options_description allowable_options(cmd_line_usage, PROGRAM_OPTIONS_LINE_LENGTH);

	try {
		standard_options.add_options()
			(PROJECT_FILE_P, boost::program_options::value<std::string>(&p.g.project_file),
				"Project file containing all user options for import, segment, geoid, reftran, adjust and plot.")
			(NETWORK_NAME_N, boost::program_options::value<std::string>(&p.g.network_name),
				"Network name. User defined name for all input and output files. Default is \"network#\".")
			(RUN_IMPORT,
				"Run import - DynAdjust file exchange software.")
			(RUN_GEOID,
				"Run geoid - geoid model interpolation software to determine geoid-ellipsoid separation and deflection values.")
			(RUN_REFTRAN,
				"Run reftran - reference frame transformation software to transform stations and measurements.")
			(RUN_SEGMENT,
				"Run segment - automated segmentation software to partition a geodetic network into smaller sized blocks.")
			(RUN_ADJUST,
				"Run adjust - geodetic network adjustment software.")
			;

		generic_options.add_options()
			(VERBOSE, boost::program_options::value<UINT16>(&p.g.verbose),
				std::string("Give detailed information about what ").append(__BINARY_NAME__).append(" is doing.\n  0: No information (default)\n  1: Helpful information\n  2: Extended information\n  3: Debug level information").c_str())
			(QUIET,
				std::string("Suppresses all explanation of what ").append(__BINARY_NAME__).append(" is doing unless an error occurs.").c_str())
			(VERSION_V, "Display the current program version.")
			(HELP_H, "Show this help message.")
			(HELP_MODULE_H, boost::program_options::value<std::string>(),
				"Provide help for a specific help category.")
			;

		allowable_options.add(standard_options).add(generic_options);

		// add "positional options" to handle command line tokens which have no option name
		positional_options.add(NETWORK_NAME, -1);
		
		boost::program_options::command_line_parser parser(argc, argv);
		store(parser.options(allowable_options).positional(positional_options).run(), vm);
		notify(vm);
	} 
	catch (const std::exception& e) {
		cout_mutex.lock();
		std::cout << "- Error: " << e.what() << std::endl;
		std::cout << cmd_line_banner << allowable_options << std::endl;
		cout_mutex.unlock();
		return EXIT_FAILURE;
	}
	
	if (argc < 2)
	{
		std::cout << std::endl << "- Nothing to do - no options provided. " << std::endl << std::endl;  
		std::cout << cmd_line_banner << allowable_options << std::endl;
		return EXIT_FAILURE;
	}

	if (vm.count(VERSION))
	{
		std::cout << cmd_line_banner << std::endl;
		return EXIT_SUCCESS;
	}

	if (vm.count(HELP))
	{
		std::cout << cmd_line_banner << allowable_options << std::endl;
		return EXIT_SUCCESS;
	}

	if (vm.count(HELP_MODULE))
	{
		std::cout << cmd_line_banner;
		std::string original_text = vm[HELP_MODULE].as<std::string>();
		std::string help_text = str_upper<std::string>(original_text);

		if (str_upper<std::string, char>(ALL_MODULE_STDOPT).find(help_text) != std::string::npos) {
			std::cout << standard_options << std::endl;
		}
		else if (str_upper<std::string, char>(ALL_MODULE_GENERIC).find(help_text) != std::string::npos) {
			std::cout << generic_options << std::endl;
		}
		else {
			std::cout << std::endl << "- Error: Help module '" <<
				original_text << "' is not in the list of options." << std::endl;
			return EXIT_FAILURE;
		}

		return EXIT_SUCCESS;
	}

	if (vm.count(QUIET))
		p.g.quiet = 1;
	else
		std::cout << std::endl << cmd_line_banner;

	if (!vm.count(PROJECT_FILE) && !vm.count(NETWORK_NAME))
	{
		std::cout << std::endl << "- Nothing to do - no project file or network name specified. " << std::endl << std::endl;  
		return EXIT_FAILURE;
	}

	if (!vm.count(RUN_IMPORT) && 
		!vm.count(RUN_GEOID) && 
		!vm.count(RUN_REFTRAN) && 
		!vm.count(RUN_SEGMENT) && 
		!vm.count(RUN_ADJUST))
	{
		std::cout << std::endl << "- Nothing to do - no programs were specified. " << std::endl << std::endl;  
		return EXIT_FAILURE;
	}

	if (vm.count(PROJECT_FILE) && vm.count(NETWORK_NAME))
	{
		if (equals(std::filesystem::path(p.g.project_file).stem().string(), p.g.network_name))
		{
			std::cout << std::endl << "- Error: project file name doesn't match network name.  Provide" << std::endl;  
			std::cout << std::endl << "         either a project file path or the network name. " << std::endl << std::endl;  
			return EXIT_FAILURE;
		}
	}
	
	// If a name was supplied, form full file path for project file using the current folder
	if (vm.count(NETWORK_NAME))
		p.g.project_file = formPath<std::string>(".", p.g.network_name, "dnaproj");
	
	if (!std::filesystem::exists(p.g.project_file))
	{
		std::cout << std::endl << 
			"- Error: Project file  " << p.g.project_file <<
				"  does not exist. " << std::endl << std::endl;  
		return EXIT_FAILURE;
	}

	// Load project file, obtain output folder and form log file path
	CDnaProjectFile projectFile;
	projectFile.LoadProjectFile(p.g.project_file);
	p = projectFile.GetSettings();
	
	std::string dynadjustLogFilePath(formPath<std::string>(p.g.output_folder, "dynadjust", "log"));
	
	p.g.log_file = dynadjustLogFilePath;
	projectFile.UpdateSettingsGeneral(p.g);
	projectFile.PrintProjectFile();

	std::ofstream dynadjust_log;

	ss << "- Error: Could not open dynadjust.log for writing." << std::endl;
	try {
		// Create dynadjust log file.  Exit failure.
		file_opener(dynadjust_log, dynadjustLogFilePath);
	}
	catch (const std::runtime_error& e) {
		ss << e.what();
		std::cout << std::endl << ss.str() << std::endl << std::endl;
		return EXIT_FAILURE;
	}
	catch (...) {
		std::cout << std::endl << ss.str() << std::endl << std::endl;
		return EXIT_FAILURE;
	}

	dynadjust_log << cmd_line_banner << std::endl;

	dynadjust_log << std::left << std::setw(15) << "+ Executable" << 
		std::setw(25) << "Start date and time" << 
		std::setw(25) << "End date and time" << 
		std::setw(25) << "Exit status" << std::endl;

	for (UINT32 d=0, linelen(15+25+25+25); d<linelen; ++d)
		dynadjust_log << "-";
	dynadjust_log << std::endl;
	
	// Run import (optional)
	if (vm.count(RUN_IMPORT))
	{
		std::stringstream cmd;
		cmd << __import_app_name__ << " -p " << p.g.project_file;
		
		if (p.g.quiet)
			cmd << " --quiet";
		
		// start time
		PrintAppStartTimeMessage(dynadjust_log, __import_app_name__);
		
		if (!run_command(cmd.str().c_str()) && !p.g.quiet)
		{
			p.i.imp_file = formPath<std::string>(p.g.output_folder, p.g.network_name, "imp");
			return CloseLogandReturn(dynadjust_log, EXIT_FAILURE, p.i.imp_file);
		}
		
		// end time
		PrintSuccessStatusMessage(dynadjust_log);

		std::this_thread::sleep_for(std::chrono::milliseconds(40));
	}
	
	// Run reftran (optional)
	if (vm.count(RUN_REFTRAN))
	{
		std::stringstream cmd;
		cmd << __reftran_app_name__ << " -p " << p.g.project_file;
		
		if (p.g.quiet)
			cmd << " --quiet";
		
		// start time
		PrintAppStartTimeMessage(dynadjust_log, __reftran_app_name__);
		
		if (!run_command(cmd.str().c_str()) && !p.g.quiet)
			return CloseLogandReturn(dynadjust_log, EXIT_FAILURE);
		
		// end time
		PrintSuccessStatusMessage(dynadjust_log);

		std::this_thread::sleep_for(std::chrono::milliseconds(40));
	}
	
	// Run geoid (optional)
	if (vm.count(RUN_GEOID))
	{
		std::stringstream cmd;
		cmd << __geoid_app_name__ << " -p " << p.g.project_file;
		
		if (p.g.quiet)
			cmd << " --quiet";
		
		// start time
		PrintAppStartTimeMessage(dynadjust_log, __geoid_app_name__);
		
		if (!run_command(cmd.str().c_str()) && !p.g.quiet)
			return CloseLogandReturn(dynadjust_log, EXIT_FAILURE);
	
		// end time
		PrintSuccessStatusMessage(dynadjust_log);

		std::this_thread::sleep_for(std::chrono::milliseconds(40));
	}
		
	// Run segment (optional)
	if (vm.count(RUN_SEGMENT))
	{
		std::stringstream cmd;
		cmd << __segment_app_name__ << " -p " << p.g.project_file;
		
		if (p.g.quiet)
			cmd << " --quiet";
		
		// start time
		PrintAppStartTimeMessage(dynadjust_log, __segment_app_name__);
		
		if (!run_command(cmd.str().c_str()) && !p.g.quiet)
			return CloseLogandReturn(dynadjust_log, EXIT_FAILURE);
	
		// end time
		PrintSuccessStatusMessage(dynadjust_log);

		std::this_thread::sleep_for(std::chrono::milliseconds(40));
	}
	
	// Run adjust (optional)
	if (vm.count(RUN_ADJUST))
	{		
		std::stringstream cmd;
		cmd << __adjust_app_name__ << " -p " << p.g.project_file;
		
		if (p.g.quiet)
			cmd << " --quiet";
		
		// start time
		PrintAppStartTimeMessage(dynadjust_log, __adjust_app_name__);
		
		if (!run_command(cmd.str().c_str()) && !p.g.quiet)
			return CloseLogandReturn(dynadjust_log, EXIT_FAILURE);

		// end time
		PrintSuccessStatusMessage(dynadjust_log);

		std::this_thread::sleep_for(std::chrono::milliseconds(40));
	}

	return CloseLogandReturn(dynadjust_log, EXIT_SUCCESS);
}

