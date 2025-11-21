// dynadjust.cpp : Defines the entry point for the console application.
//
#include <iostream>
#include <fstream>
#include <algorithm>
#include <iomanip>
#include <sstream>
#include <string>
#include <time.h>

#include <boost/timer/timer.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/thread.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/filesystem.hpp>
#include <libgen.h>
#include <boost/algorithm/string/predicate.hpp>

using namespace std;
using namespace boost;
using namespace boost::filesystem;
using namespace boost::posix_time;
using namespace boost::program_options;

boost::mutex cout_mutex;

#include <include/config/dnaversion.hpp>
#include <include/config/dnaconsts.hpp>
#include <include/config/dnaoptions.hpp>
#include <include/config/dnaoptions-interface.hpp>

#include <include/functions/dnastringfuncs.hpp>
#include <include/functions/dnaprocessfuncs.hpp>
#include <include/functions/dnaiostreamfuncs.hpp>
#include <include/functions/dnafilepathfuncs.hpp>
#include <include/functions/dnafilepathfuncs.hpp>
#include <include/functions/dnatemplatedatetimefuncs.hpp>
#include <include/functions/dnastrmanipfuncs.hpp>

#include <include/config/dnaprojectfile.hpp>

using namespace dynadjust;

int CloseLogandReturn(std::ofstream& os, const int status, const string& output_file="")
{
	if (status == EXIT_FAILURE)
	{
		os << setw(25) << left << formattedDateTimeString<string>() << 
			"Failed. " << endl << endl << 
			"+ DynAdjust ended prematurely. ";

		if (!output_file.empty())
		{
			std::ifstream f(output_file);
			os << " Contents of last output file:" << endl << endl << f.rdbuf();
		}
		else
			 os << endl;
	}
	else
		os << endl << endl << "+ DynAdjust finished successfully." << endl;

	os.close();
	return status;
}

void PrintAppStartTimeMessage(std::ofstream& os, string&& app)
{
	// App start time
	os << setw(15) << left << string("+ " + app) << setw(25) << 
		formattedDateTimeString<string>();
}

void PrintSuccessStatusMessage(std::ofstream& os)
{
	// end time
	os << setw(25) << left << formattedDateTimeString<string>() << 
		"Ended successfully." << endl;
}

int main(int argc, char* argv[])
{
	stringstream ss;

	string cmd_line_banner;
	fileproc_help_header(&cmd_line_banner);

	project_settings p;
	variables_map vm;
	positional_options_description positional_options;

	options_description standard_options("+ " + string(ALL_MODULE_STDOPT), PROGRAM_OPTIONS_LINE_LENGTH);
	options_description generic_options("+ " + string(ALL_MODULE_GENERIC), PROGRAM_OPTIONS_LINE_LENGTH);

	string cmd_line_usage("+ ");
	cmd_line_usage.append(__BINARY_NAME__).append(" usage:  ").append(__BINARY_NAME__).append(" ").append(" [options]");
	options_description allowable_options(cmd_line_usage, PROGRAM_OPTIONS_LINE_LENGTH);

	try {
		standard_options.add_options()
			(PROJECT_FILE_P, value<string>(&p.g.project_file),
				"Project file containing all user options for import, segment, geoid, reftran, adjust and plot.")
			(NETWORK_NAME_N, value<string>(&p.g.network_name),
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
			(VERBOSE, value<UINT16>(&p.g.verbose),
				string("Give detailed information about what ").append(__BINARY_NAME__).append(" is doing.\n  0: No information (default)\n  1: Helpful information\n  2: Extended information\n  3: Debug level information").c_str())
			(QUIET,
				string("Suppresses all explanation of what ").append(__BINARY_NAME__).append(" is doing unless an error occurs.").c_str())
			(VERSION_V, "Display the current program version.")
			(HELP_H, "Show this help message.")
			(HELP_MODULE, value<string>(),
				"Provide help for a specific help category.")
			;

		allowable_options.add(standard_options).add(generic_options);

		// add "positional options" to handle command line tokens which have no option name
		positional_options.add(NETWORK_NAME, -1);
		
		command_line_parser parser(argc, argv);
		store(parser.options(allowable_options).positional(positional_options).run(), vm);
		notify(vm);
	} 
	catch (const std::exception& e) {
		cout_mutex.lock();
		cout << "- Error: " << e.what() << endl;
		cout << cmd_line_banner << allowable_options << endl;
		cout_mutex.unlock();
		return EXIT_FAILURE;
	}
	
	if (argc < 2)
	{
		cout << endl << "- Nothing to do - no options provided. " << endl << endl;  
		cout << cmd_line_banner << allowable_options << endl;
		return EXIT_FAILURE;
	}

	if (vm.count(VERSION))
	{
		cout << cmd_line_banner << endl;
		return EXIT_SUCCESS;
	}

	if (vm.count(HELP))
	{
		cout << cmd_line_banner << allowable_options << endl;
		return EXIT_SUCCESS;
	}

	if (vm.count(HELP_MODULE))
	{
		cout << cmd_line_banner;
		string original_text = vm[HELP_MODULE].as<string>();
		string help_text = str_upper<string>(original_text);

		if (str_upper<string, char>(ALL_MODULE_STDOPT).find(help_text) != string::npos) {
			cout << standard_options << endl;
		}
		else if (str_upper<string, char>(ALL_MODULE_GENERIC).find(help_text) != string::npos) {
			cout << generic_options << endl;
		}
		else {
			cout << endl << "- Error: Help module '" <<
				original_text << "' is not in the list of options." << endl;
			return EXIT_FAILURE;
		}

		return EXIT_SUCCESS;
	}

	if (vm.count(QUIET))
		p.g.quiet = 1;
	else
		cout << endl << cmd_line_banner;

	if (!vm.count(PROJECT_FILE) && !vm.count(NETWORK_NAME))
	{
		cout << endl << "- Nothing to do - no project file or network name specified. " << endl << endl;  
		return EXIT_FAILURE;
	}

	if (!vm.count(RUN_IMPORT) && 
		!vm.count(RUN_GEOID) && 
		!vm.count(RUN_REFTRAN) && 
		!vm.count(RUN_SEGMENT) && 
		!vm.count(RUN_ADJUST))
	{
		cout << endl << "- Nothing to do - no programs were specified. " << endl << endl;  
		return EXIT_FAILURE;
	}

	if (vm.count(PROJECT_FILE) && vm.count(NETWORK_NAME))
	{
		if (equals(path(p.g.project_file).stem().string(), p.g.network_name))
		{
			cout << endl << "- Error: project file name doesn't match network name.  Provide" << endl;  
			cout << endl << "         either a project file path or the network name. " << endl << endl;  
			return EXIT_FAILURE;
		}
	}
	
	// If a name was supplied, form full file path for project file using the current folder
	if (vm.count(NETWORK_NAME))
		p.g.project_file = formPath<string>(".", p.g.network_name, "dnaproj");
	
	if (!exists(p.g.project_file))
	{
		cout << endl << 
			"- Error: Project file  " << p.g.project_file <<
				"  does not exist. " << endl << endl;  
		return EXIT_FAILURE;
	}

	// Load project file, obtain output folder and form log file path
	CDnaProjectFile projectFile;
	projectFile.LoadProjectFile(p.g.project_file);
	p = projectFile.GetSettings();
	
	string dynadjustLogFilePath(formPath<string>(p.g.output_folder, "dynadjust", "log"));
	
	p.g.log_file = dynadjustLogFilePath;
	projectFile.UpdateSettingsGeneral(p.g);
	projectFile.PrintProjectFile();

	std::ofstream dynadjust_log;

	ss << "- Error: Could not open dynadjust.log for writing." << endl;
	try {
		// Create dynadjust log file.  Exit failure.
		file_opener(dynadjust_log, dynadjustLogFilePath);
	}
	catch (const runtime_error& e) {
		ss << e.what();
		cout << endl << ss.str() << endl << endl;
		return EXIT_FAILURE;
	}
	catch (...) {
		cout << endl << ss.str() << endl << endl;
		return EXIT_FAILURE;
	}

	dynadjust_log << cmd_line_banner << endl;

	dynadjust_log << left << setw(15) << "+ Executable" << 
		setw(25) << "Start date and time" << 
		setw(25) << "End date and time" << 
		setw(25) << "Exit status" << endl;

	for (UINT32 d=0, linelen(15+25+25+25); d<linelen; ++d)
		dynadjust_log << "-";
	dynadjust_log << endl;
	
	// Run import (optional)
	if (vm.count(RUN_IMPORT))
	{
		stringstream cmd;
		cmd << __import_app_name__ << " -p " << p.g.project_file;
		
		if (p.g.quiet)
			cmd << " --quiet";
		
		// start time
		PrintAppStartTimeMessage(dynadjust_log, __import_app_name__);
		
		if (!run_command(cmd.str().c_str()), p.g.quiet)
		{
			p.i.imp_file = formPath<string>(p.g.output_folder, p.g.network_name, "imp");
			return CloseLogandReturn(dynadjust_log, EXIT_FAILURE, p.i.imp_file);
		}
		
		// end time
		PrintSuccessStatusMessage(dynadjust_log);

		boost::this_thread::sleep(milliseconds(40));
	}
	
	// Run reftran (optional)
	if (vm.count(RUN_REFTRAN))
	{
		stringstream cmd;
		cmd << __reftran_app_name__ << " -p " << p.g.project_file;
		
		if (p.g.quiet)
			cmd << " --quiet";
		
		// start time
		PrintAppStartTimeMessage(dynadjust_log, __reftran_app_name__);
		
		if (!run_command(cmd.str().c_str()), p.g.quiet)
			return CloseLogandReturn(dynadjust_log, EXIT_FAILURE);
		
		// end time
		PrintSuccessStatusMessage(dynadjust_log);

		boost::this_thread::sleep(milliseconds(40));
	}
	
	// Run geoid (optional)
	if (vm.count(RUN_GEOID))
	{
		stringstream cmd;
		cmd << __geoid_app_name__ << " -p " << p.g.project_file;
		
		if (p.g.quiet)
			cmd << " --quiet";
		
		// start time
		PrintAppStartTimeMessage(dynadjust_log, __geoid_app_name__);
		
		if (!run_command(cmd.str().c_str()), p.g.quiet)
			return CloseLogandReturn(dynadjust_log, EXIT_FAILURE);
	
		// end time
		PrintSuccessStatusMessage(dynadjust_log);

		boost::this_thread::sleep(milliseconds(40));
	}
		
	// Run segment (optional)
	if (vm.count(RUN_SEGMENT))
	{
		stringstream cmd;
		cmd << __segment_app_name__ << " -p " << p.g.project_file;
		
		if (p.g.quiet)
			cmd << " --quiet";
		
		// start time
		PrintAppStartTimeMessage(dynadjust_log, __segment_app_name__);
		
		if (!run_command(cmd.str().c_str()), p.g.quiet)
			return CloseLogandReturn(dynadjust_log, EXIT_FAILURE);
	
		// end time
		PrintSuccessStatusMessage(dynadjust_log);

		boost::this_thread::sleep(milliseconds(40));
	}
	
	// Run adjust (optional)
	if (vm.count(RUN_ADJUST))
	{		
		stringstream cmd;
		cmd << __adjust_app_name__ << " -p " << p.g.project_file;
		
		if (p.g.quiet)
			cmd << " --quiet";
		
		// start time
		PrintAppStartTimeMessage(dynadjust_log, __adjust_app_name__);
		
		if (!run_command(cmd.str().c_str()), p.g.quiet)
			return CloseLogandReturn(dynadjust_log, EXIT_FAILURE);

		// end time
		PrintSuccessStatusMessage(dynadjust_log);

		boost::this_thread::sleep(milliseconds(40));
	}

	return CloseLogandReturn(dynadjust_log, EXIT_SUCCESS);
}

