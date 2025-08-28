//============================================================================
// Name         : dnasegmentwrapper.cpp
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
// Description  : DynAdjust Network Segmentation Executable
//============================================================================

#include <dynadjust/dnasegmentwrapper/dnasegmentwrapper.hpp>
#include <mutex>
#include <thread>

using namespace dynadjust;

bool running;
std::mutex cout_mutex;

int ParseCommandLineOptions(const int& argc, char* argv[], const boost::program_options::variables_map& vm, project_settings& p)
{
	// capture command line arguments
	for (int cmd_arg(0); cmd_arg<argc; ++cmd_arg)
	{
		 p.s.command_line_arguments += argv[cmd_arg];
		 p.s.command_line_arguments += " ";
	}

	if (vm.count(PROJECT_FILE))
	{
		if (std::filesystem::exists(p.g.project_file))
		{
			try {
				CDnaProjectFile projectFile(p.g.project_file, segmentSetting);
				p = projectFile.GetSettings();
			}
			catch (const std::runtime_error& e) {
				std::cout << std::endl << "- Error: " << e.what() << std::endl;
				return EXIT_FAILURE;
			}
			
			return EXIT_SUCCESS;
		}

		std::cout << std::endl << "- Error: project file " << p.g.project_file << " does not exist." << std::endl << std::endl;
		return EXIT_FAILURE;
	}

	if (!vm.count(NETWORK_NAME))
	{
		std::cout << std::endl << "- Nothing to do - no network name specified. " << std::endl << std::endl;  
		return EXIT_FAILURE;
	}
	
	p.g.project_file = formPath<std::string>(p.g.output_folder, p.g.network_name, "dnaproj");

	// input files
	p.s.asl_file = formPath<std::string>(p.g.input_folder, p.g.network_name, "asl");		// associated stations list
	p.s.aml_file = formPath<std::string>(p.g.input_folder, p.g.network_name, "aml");		// associated measurements list
	p.s.map_file = formPath<std::string>(p.g.input_folder, p.g.network_name, "map");		// station names map
	
	if (vm.count(NET_FILE))
		p.s.net_file = formPath<std::string>(p.g.input_folder, p.g.network_name, "net");		// Starting stations file
	
	// binary station file location (input)
	if (vm.count(BIN_STN_FILE))
		p.s.bst_file = formPath<std::string>(p.g.input_folder, p.s.bst_file);
	else
		p.s.bst_file = formPath<std::string>(p.g.output_folder, p.g.network_name, "bst");

	// binary station file location (input)
	if (vm.count(BIN_MSR_FILE))
		p.s.bms_file = formPath<std::string>(p.g.input_folder, p.s.bms_file);
	else
		p.s.bms_file = formPath<std::string>(p.g.output_folder, p.g.network_name, "bms");

	if (!std::filesystem::exists(p.s.bst_file) || !std::filesystem::exists(p.s.bms_file))
	{
		std::cout << std::endl << "- Nothing to do: ";  
			
		if (p.g.network_name.empty())
			std::cout << std::endl << "network name has not been specified specified, and " << std::endl << "               ";  
		std::cout << p.s.bst_file << " and " << p.s.bms_file << " do not exist." << std::endl << std::endl;  
		return EXIT_FAILURE;
	}

	// output files
	// User supplied segmentation file
	if (vm.count(SEG_FILE))
	{
		// Does it exist?
		if (!std::filesystem::exists(p.s.seg_file))
		{
			// Look for it in the input folder
			p.s.seg_file = formPath<std::string>(p.g.input_folder, leafStr<std::string>(p.s.seg_file));

			if (!std::filesystem::exists(p.s.seg_file))
			{
				std::cout << std::endl << "- Error: " <<
					"Segmentation file " << leafStr<std::string>(p.s.seg_file) << " does not exist." << std::endl << std::endl;  
				return EXIT_FAILURE;
			}
		}
	}
	else
		p.s.seg_file = formPath<std::string>(p.g.input_folder, p.g.network_name, "seg");

	// Station appearance file
	p.s.sap_file = formPath<std::string>(p.g.input_folder, p.g.network_name, "sap");
	
	if (vm.count(TEST_INTEGRITY))
		p.i.test_integrity = 1;

	//if (vm.count(SEG_FORCE_CONTIGUOUS))
	//	p.s.force_contiguous_blocks = 1;

	return EXIT_SUCCESS;
}

int main(int argc, char* argv[])
{	
	// create banner message
	std::string cmd_line_banner;
	fileproc_help_header(&cmd_line_banner);

	std::string stnfilename, msrfilename;
	
	project_settings p;

	boost::program_options::variables_map vm;
	boost::program_options::positional_options_description positional_options;

	boost::program_options::options_description standard_options("+ " + std::string(ALL_MODULE_STDOPT), PROGRAM_OPTIONS_LINE_LENGTH);
	boost::program_options::options_description config_options("+ " + std::string(SEGMENT_MODULE_CONFIG), PROGRAM_OPTIONS_LINE_LENGTH);
	boost::program_options::options_description generic_options("+ " + std::string(ALL_MODULE_GENERIC), PROGRAM_OPTIONS_LINE_LENGTH);
	
	std::string cmd_line_usage("+ ");
	cmd_line_usage.append(__BINARY_NAME__).append(" usage:  ").append(__BINARY_NAME__).append(" ").append(NETWORK_NAME).append(" [options]");
	boost::program_options::options_description allowable_options(cmd_line_usage, PROGRAM_OPTIONS_LINE_LENGTH);

	try {
		// Declare a group of options that will be 
		// allowed only on command line		
		standard_options.add_options()
			(PROJECT_FILE_P, boost::program_options::value<std::string>(&p.g.project_file),
				"Project file name. Full path to project file. If none specified, a new file is created using input-folder and network-name.")
			(NETWORK_NAME_N, boost::program_options::value<std::string>(&p.g.network_name),
				"Network name. User defined name for all input and output files. Default is \"network#\".")
			(INPUT_FOLDER_I, boost::program_options::value<std::string>(&p.g.input_folder),
				"Path containing all input files")
			(OUTPUT_FOLDER_O, boost::program_options::value<std::string>(&p.g.output_folder),		// default is ./,
				"Path for all output files")
			(BIN_STN_FILE, boost::program_options::value<std::string>(&p.s.bst_file),
				"Binary station file name. Overrides network name.")
			(BIN_MSR_FILE, boost::program_options::value<std::string>(&p.s.bms_file),
				"Binary measurement file name. Overrides network name.")
			(SEG_FILE, boost::program_options::value<std::string>(&p.s.seg_file),
				"Segmentation output file name. Overrides network name.")
			;

		// Declare a group of options that will be 
		// allowed both on command line and in
		// config options        
		config_options.add_options()
			(NET_FILE,
				"Look for a .net file containing stations to be incorporated within the first block.")
			(SEG_STARTING_STN, boost::program_options::value<std::string>(&p.s.seg_starting_stns),
				"Additional stations to be incorporated within the first block. arg is a comma delimited string \"stn1, stn 2,stn3 , stn 4\".")
			(SEG_MIN_INNER_STNS, boost::program_options::value<UINT32>(&p.s.min_inner_stations),
				(std::string("Minimum number of inner stations within each block. Default is ")+
					StringFromT(p.s.min_inner_stations)+std::string(".")).c_str())
			(SEG_THRESHOLD_STNS, boost::program_options::value<UINT32>(&p.s.max_total_stations),
				(std::string("Threshold limit for maximum number of stations per block. Default is ")+
					StringFromT(p.s.max_total_stations)+std::string(".")).c_str())
			(SEG_FORCE_CONTIGUOUS, boost::program_options::value<UINT16>(&p.s.force_contiguous_blocks),
				(std::string("Treatment of isolated networks:\n")+
				std::string("  0: Isolated networks as individual blocks ")+
				(p.s.force_contiguous_blocks==0 ? "(default)\n" : "\n")+
				std::string("  1: Force production of contiguous blocks ")+
				(p.s.force_contiguous_blocks==1 ? "(default)" : "")
				).c_str())
			(SEG_SEARCH_LEVEL, boost::program_options::value<UINT16>(&p.s.seg_search_level),
				"Level to which searches should be conducted to find stations with the lowest measurement count. Default is 0.")
			(TEST_INTEGRITY,
				"Test the integrity of all output files.")
			;

		generic_options.add_options()
			(VERBOSE, boost::program_options::value<UINT16>(&p.g.verbose),
				std::string("Give detailed information about what ").append(__BINARY_NAME__).append(" is doing.\n  0: No information (default)\n  1: Helpful information\n  2: Extended information\n  3: Debug level information").c_str())
				(QUIET,
					std::string("Suppresses all explanation of what ").append(__BINARY_NAME__).append(" is doing unless an error occurs").c_str())
					(VERSION_V, "Display the current program version")
			(HELP_H, "Show this help message")
			(HELP_MODULE_H, boost::program_options::value<std::string>(),
				"Provide help for a specific help category.")
			;

		allowable_options.add(standard_options).add(config_options).add(generic_options);

		// add "positional options" to handle command line tokens which have no option name
		positional_options.add(NETWORK_NAME, -1);
		
		boost::program_options::command_line_parser parser(argc, argv);
		store(parser.options(allowable_options).positional(positional_options).run(), vm);
		notify(vm);
	} 
	catch (const std::exception& e) {
		std::cout << "- Error: " << e.what() << std::endl <<
			cmd_line_banner << allowable_options << std::endl;
		return EXIT_FAILURE;
	}

	if (argc < 2)
	{
		std::cout << std::endl << "- Nothing to do - no options provided. " << std::endl << std::endl <<
			cmd_line_banner << allowable_options << std::endl;
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
		else if (str_upper<std::string, char>(SEGMENT_MODULE_CONFIG).find(help_text) != std::string::npos) {
			std::cout << config_options << std::endl;
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

	bool userSuppliedSegFile(false);
	if (!p.s.seg_file.empty())
		userSuppliedSegFile = true;
	bool userSuppliedBstFile(false);
	if (!p.s.bst_file.empty())
		userSuppliedBstFile = true;
	bool userSuppliedBmsFile(false);
	if (!p.s.bms_file.empty())
		userSuppliedBmsFile = true;

	if (ParseCommandLineOptions(argc, argv, vm, p) != EXIT_SUCCESS)
		return EXIT_FAILURE;

	// Has the user supplied a maximum block size that is less 
	// than the minimum block size?
	if (p.s.min_inner_stations > p.s.max_total_stations)
		p.s.min_inner_stations = p.s.max_total_stations;

	if (vm.count(QUIET))
		p.g.quiet = 1;
	
	if (!p.g.quiet)
	{
		std::cout << std::endl << cmd_line_banner;
		
		std::cout << "+ Options:" << std::endl; 
		std::cout << std::setw(PRINT_VAR_PAD) << std::left << "  Network name: " <<  p.g.network_name << std::endl;
		std::cout << std::setw(PRINT_VAR_PAD) << std::left << "  Input folder: " << p.g.input_folder << std::endl;
		std::cout << std::setw(PRINT_VAR_PAD) << std::left << "  Output folder: " << p.g.output_folder << std::endl;
		std::cout << std::setw(PRINT_VAR_PAD) << std::left << "  Associated station file: " << p.s.asl_file << std::endl;
		std::cout << std::setw(PRINT_VAR_PAD) << std::left << "  Associated measurement file: " << p.s.aml_file << std::endl;
		std::cout << std::setw(PRINT_VAR_PAD) << std::left << "  Binary station file: " << p.s.bst_file << std::endl;
		std::cout << std::setw(PRINT_VAR_PAD) << std::left << "  Binary measurement file: " << p.s.bms_file << std::endl;
		std::cout << std::setw(PRINT_VAR_PAD) << std::left << "  Segmentation output file: " << p.s.seg_file << std::endl;
		if (!p.s.net_file.empty())
			std::cout << std::setw(PRINT_VAR_PAD) << std::left << "  Block 1 stations file: " << p.s.net_file << std::endl;
		std::cout << std::setw(PRINT_VAR_PAD) << std::left << "  Minimum inner stations: " <<  p.s.min_inner_stations << std::endl;
		std::cout << std::setw(PRINT_VAR_PAD) << std::left << "  Block size threshold: " <<  p.s.max_total_stations << std::endl;
		if (!p.s.seg_starting_stns.empty())
			std::cout << std::setw(PRINT_VAR_PAD) << std::left << "  Additional Block 1 stations: " << p.s.seg_starting_stns << std::endl;
		else if (p.s.net_file.empty())
			std::cout << "  No initial station specified. The first station will be used." << std::endl;
		
		std::cout << std::endl;
	}

	// Should segment look for a net file?
	if (!p.s.net_file.empty())
	{
		if (!std::filesystem::exists(p.s.net_file))
		{
			cout_mutex.lock();
			std::cout << std::endl <<
				"- Error: " << p.s.net_file << " does not exist." << std::endl <<
				"  A file named " << p.g.network_name << ".net must exist in the input folder\n  in order to use this option." << std::endl << std::endl;
			cout_mutex.unlock();
			return EXIT_FAILURE;
		}
	}

	dna_segment netSegment;
	boost::posix_time::milliseconds elapsed_time(boost::posix_time::milliseconds(0));
	_SEGMENT_STATUS_ segmentStatus;
	std::string status_msg;

	try {
		netSegment.InitialiseSegmentation();
		running = true;

		// segment blocks using threads
		std::vector<std::thread> ui_segment_threads;
		if (!p.g.quiet)
			ui_segment_threads.emplace_back(dna_segment_progress_thread(&netSegment, &p));
		ui_segment_threads.emplace_back(dna_segment_thread(&netSegment, &p, &segmentStatus, &elapsed_time, &status_msg));
		for (auto& t : ui_segment_threads) {
			t.join();
		}
		
		switch (netSegment.GetStatus())
		{
		case SEGMENT_EXCEPTION_RAISED:
			running = false;
			return EXIT_FAILURE;
		default:
			break;
		}

		//// print station appearance file
		//if (!p.g.quiet)
		//	std::cout << "+ Printing station appearance list... ";
		//netSegment.WriteStationAppearanceList(p.s.sap_file);
		//if (!p.g.quiet)
		//	std::cout << "done." << std::endl;

		if (p.g.verbose > 1 && !p.g.quiet)
			netSegment.coutSummary(); 

		if (segmentStatus != SEGMENT_SUCCESS)
			std::cout << status_msg << std::endl;

		if (!p.g.quiet)
		{
			std::cout << "+ Segmentation statistics:" << std::endl;
			std::cout << std::endl << std::left << "  " <<
				std::setw(STATION) << "No. blocks" << 
				std::setw(STAT) << "Max size" << 
				std::setw(STAT) << "Min size" << 
				std::setw(STAT) << "Average" <<
				std::setw(STATION) << "Total size" << std::endl;
			std::cout << "  ";
			for (UINT32 i(0), j(STATION + STAT*3 + STATION); i<j; ++i)
				std::cout << "-";
			std::cout << std::endl << "  " << std::left <<
				std::setw(STATION) << netSegment.blockCount() << 
				std::setw(STAT) << std::setprecision(0) << netSegment.maxBlockSize() << 
				std::setw(STAT) << std::setprecision(0) << netSegment.minBlockSize() << 
				std::setw(STAT) << std::setprecision(2) << std::fixed << netSegment.averageblockSize() <<
				std::setw(STATION) << std::setprecision(0) << netSegment.stationSolutionCount() << std::endl;

			std::cout << std::endl;
		
		}

		if (!p.g.quiet)
			std::cout << "+ Verifying station connections... ";
		netSegment.VerifyStationConnections();
		if (!p.g.quiet)
			std::cout << "done." << std::endl;

		// print network segmentation block file
		if (!p.g.quiet)
			std::cout << "+ Printing blocks to " << leafStr<std::string>(p.s.seg_file) << "... ";
		netSegment.WriteSegmentedNetwork(p.s.seg_file);
		if (!p.g.quiet)
			std::cout << "done." << std::endl;

		
	} 
	catch (const NetSegmentException& e) {
		std::cout << std::endl << "- Error: " << e.what() << std::endl;
		return EXIT_FAILURE;
	} 
	catch (const std::runtime_error& e) {
		std::cout << "+ Exception of unknown type: " << e.what();
		return EXIT_FAILURE;
	}

	if (!userSuppliedSegFile)
		p.s.seg_file = "";
	if (!userSuppliedBstFile)
		p.s.bst_file = "";
	if (!userSuppliedBmsFile)
		p.s.bms_file = "";

	// Look for a project file.  If it exists, open and load it.
	// Update the import settings.
	// Print the project file. If it doesn't exist, it will be created.
	CDnaProjectFile projectFile;
	if (std::filesystem::exists(p.g.project_file))
		projectFile.LoadProjectFile(p.g.project_file);
	
	projectFile.UpdateSettingsSegment(p);
	projectFile.PrintProjectFile();

	if (p.g.quiet)
		return EXIT_SUCCESS;

	std::cout << std::endl << formatedElapsedTime<std::string>(&elapsed_time, "+ Network segmentation took ") << std::endl;
	std::cout << "+ " << p.g.network_name << " is now ready for sequential phased adjustment." << std::endl << std::endl;
	
	return EXIT_SUCCESS;
}
