//============================================================================
// Name         : dnaimportwrapper.cpp
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
// Description  : DynAdjust Interoperability library Executable
//============================================================================

#include <filesystem>
#include <dynadjust/dnaimportwrapper/dnaimportwrapper.hpp>
#include <include/parameters/dnaepsg.hpp>
#include <include/functions/dnastrutils.hpp>
#include <mutex>
#include <thread>

using namespace dynadjust;
using namespace dynadjust::epsg;

bool running;
std::mutex cout_mutex;

// Helper function to convert std::filesystem::file_time_type to time_t
time_t file_time_to_time_t(const std::filesystem::file_time_type& ftime) {
    using namespace std::chrono;
    auto sctp = time_point_cast<system_clock::duration>(ftime - std::filesystem::file_time_type::clock::now() + system_clock::now());
    return system_clock::to_time_t(sctp);
}

void PrintOutputFileHeaderInfo(std::ofstream* f_out, const std::string& out_file, project_settings* p, const std::string& header)
{
	// Print formatted header
	print_file_header(*f_out, header);

	*f_out << std::setw(PRINT_VAR_PAD) << std::left << "File name:" << std::filesystem::absolute(out_file).string() << std::endl << std::endl;
	
	*f_out << std::setw(PRINT_VAR_PAD) << std::left << "Command line arguments: ";
	*f_out << p->i.command_line_arguments << std::endl << std::endl;

	*f_out << std::setw(PRINT_VAR_PAD) << std::left << "Network name:" <<  p->g.network_name << std::endl;
	*f_out << std::setw(PRINT_VAR_PAD) << std::left << "Input folder: " << std::filesystem::absolute(p->g.input_folder).string() << std::endl;
	*f_out << std::setw(PRINT_VAR_PAD) << std::left << "Output folder: " << std::filesystem::absolute(p->g.output_folder).string() << std::endl;
	*f_out << std::setw(PRINT_VAR_PAD) << std::left << "Stations file:" << std::filesystem::absolute(p->i.bst_file).string() << std::endl;
	*f_out << std::setw(PRINT_VAR_PAD) << std::left << "Measurements file:" << std::filesystem::absolute(p->i.bms_file).string() << std::endl;
	*f_out << std::setw(PRINT_VAR_PAD) << std::left << "Associated station file:" << std::filesystem::absolute(p->i.asl_file).string() << std::endl;
	*f_out << std::setw(PRINT_VAR_PAD) << std::left << "Associated measurement file:" << std::filesystem::absolute(p->i.aml_file).string() << std::endl;
	*f_out << std::setw(PRINT_VAR_PAD) << std::left << "Duplicate stations output file:" <<  std::filesystem::absolute(p->i.dst_file).string() << std::endl;
	*f_out << std::setw(PRINT_VAR_PAD) << std::left << "Similar measurement output file:" <<  std::filesystem::absolute(p->i.dms_file).string() << std::endl;

	if (!p->i.input_files.empty())
	{
		_it_vstr _it_files(p->i.input_files.begin());
		std::string str("Input files:");
		while (_it_files!=p->i.input_files.end())
		{
			*f_out << std::setw(PRINT_VAR_PAD) << std::left << str << *_it_files++ << std::endl;
			str = " ";
		}
	}

	// If a reference frame has been supplied, report it.  
	// If not, the assumption is, the project frame will be adopted from the first file and
	// in this case, it will be reported later
	if (p->i.user_supplied_frame)
		*f_out << std::setw(PRINT_VAR_PAD) << std::left << "Project reference frame:" << p->i.reference_frame << " (user supplied)" << std::endl;
	else
		*f_out << std::setw(PRINT_VAR_PAD) << std::left << "Project reference frame:" << "To be adopted from the first input file" << std::endl;

	if (p->i.override_input_rfame)
		*f_out << std::setw(PRINT_VAR_PAD) << std::left << "Override input file ref frame:" << yesno_string(p->i.override_input_rfame) << std::endl;

	UINT32 epsgCode(epsgCodeFromName<UINT32, std::string>(p->i.reference_frame));

	if (isEpsgDatumStatic(epsgCode) && p->i.user_supplied_frame)
		*f_out << std::setw(PRINT_VAR_PAD) << std::left << "Project epoch:" << p->i.epoch << " (adopted reference epoch of " << p->i.reference_frame << ")" << std::endl;
	else
	{
		if (p->i.user_supplied_epoch)
			*f_out << std::setw(PRINT_VAR_PAD) << std::left << "Project epoch:" << p->i.epoch << " (user supplied)" << std::endl;
		else
			*f_out << std::setw(PRINT_VAR_PAD) << std::left << "Project epoch:" << "To be adopted from the first input file" << std::endl;
	}
	
	if (!p->i.include_msrs.empty())
		*f_out << std::setw(PRINT_VAR_PAD) << std::left << "Strip all measurements except:" << p->i.include_msrs << std::endl;
	else if (!p->i.exclude_msrs.empty())
		*f_out << std::setw(PRINT_VAR_PAD) << std::left << "Strip measurement types:" << p->i.exclude_msrs << std::endl;

	if (p->i.ignore_insufficient_msrs == 1)
		*f_out << std::setw(PRINT_VAR_PAD) << std::left << "Ignore insufficient measurements:" << "yes" << std::endl;
	
	if (p->i.search_nearby_stn)
		*f_out << std::setw(PRINT_VAR_PAD) << std::left << "Search for nearby stations:" << "tolerance = " << p->i.search_stn_radius << "m" << std::endl;
	if (p->i.search_similar_msr)
		*f_out << std::setw(PRINT_VAR_PAD) << std::left << "Search for similar measurements:" << "yes" << std::endl;
	if (p->i.search_similar_msr_gx)
		*f_out << std::setw(PRINT_VAR_PAD) << std::left << "Search for similar GNSS measurements:" << "yes" << std::endl;
	
	if (!p->i.bounding_box.empty())
	{
		*f_out << std::setw(PRINT_VAR_PAD) << std::left << "Bounding box: " << p->i.bounding_box << std::endl;
		if (p->i.split_clusters)
			*f_out << std::setw(PRINT_VAR_PAD) << std::left << "Split GNSS clusters: " << (p->i.split_clusters ? "Yes" : "No") << std::endl;
	}
	else 
	{
		if (!p->i.stn_associated_msr_include.empty())
			*f_out << std::setw(PRINT_VAR_PAD) << std::left << "Stations to include: " << p->i.stn_associated_msr_include << std::endl;
		if (!p->i.stn_associated_msr_exclude.empty())
			*f_out << std::setw(PRINT_VAR_PAD) << std::left << "Stations to exclude: " << p->i.stn_associated_msr_exclude << std::endl;

		if (p->i.split_clusters)
			*f_out << std::setw(PRINT_VAR_PAD) << std::left << "Split GNSS clusters: " << (p->i.split_clusters ? "Yes" : "No") << std::endl;
	}

	if (!p->i.seg_file.empty())
		*f_out << std::setw(PRINT_VAR_PAD) << std::left << "Segmentation file:" << std::filesystem::absolute(p->i.seg_file).string() << std::endl;

	if (p->i.import_block)
		*f_out << std::setw(PRINT_VAR_PAD) << std::left << "Import stns & msrs from block: " << p->i.import_block_number << std::endl;
	else if (p->i.import_network)
		*f_out << std::setw(PRINT_VAR_PAD) << std::left << "Import stns & msrs from network: " << p->i.import_network_number << std::endl;
	
	*f_out << OUTPUTLINE << std::endl << std::endl;
}
	
int ParseCommandLineOptions(const int& argc, char* argv[], const boost::program_options::variables_map& vm, project_settings& p)
{
	// capture command line arguments
	for (int cmd_arg(0); cmd_arg<argc; ++cmd_arg)
	{
		 p.i.command_line_arguments += argv[cmd_arg];
		 p.i.command_line_arguments += " ";
	}

	if (vm.count(PROJECT_FILE))
	{
		if (std::filesystem::exists(p.g.project_file))
		{
			try {
				CDnaProjectFile projectFile(p.g.project_file, importSetting);
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

	if (!vm.count(IMPORT_FILE) && !vm.count(IMPORT_SEG_BLOCK) && !vm.count(IMPORT_CONTIG_NET) && !vm.count(SEG_FILE))
	{
		std::cout << std::endl << "- Nothing to do - no files specified. " << std::endl << std::endl;  
		return EXIT_FAILURE;
	}

	if (vm.count(EPOCH) && !vm.count(REFERENCE_FRAME))
	{
		std::cout << std::endl << "- A reference frame must be provided when providing an epoch. " << std::endl << std::endl;
		return EXIT_FAILURE;
	}

	// Normalise files using input folder
	for_each(p.i.input_files.begin(), p.i.input_files.end(),
		[&p] (std::string& file) { 
			formPath<std::string>(p.g.input_folder, file);
		}
	);

	//////////////////////////////////////////////////////////////////////////////
	// General options and file paths
	// Network name
	if (p.g.network_name == "network1")
	{
		// Iterate through network1, network2, network3, etc 
		// until the first name not used is found
		std::stringstream netname_ss;
		std::string networkASL;
		UINT32 networkID(1);
		
		// This loop terminates when a file name cannot be found
		while (true)
		{
			// 1. Form ASL file and see if it exists
			networkASL = formPath<std::string>(p.g.output_folder, p.g.network_name, "asl");

			// 2. Does this network exist?
			if (!std::filesystem::exists(networkASL))
				break;

			// 3. Flush and look for the next network
			netname_ss.str("");
			netname_ss << "network" << ++networkID;
			p.g.network_name = netname_ss.str();			
		}
	}

	p.g.project_file = formPath<std::string>(p.g.output_folder, p.g.network_name, "dnaproj");

	// binary station file location (output)
	if (vm.count(BIN_STN_FILE))
		p.i.bst_file = formPath<std::string>(p.g.output_folder, p.i.bst_file);
	else
		p.i.bst_file = formPath<std::string>(p.g.output_folder, p.g.network_name, "bst");
	
	// binary station file location (output)
	if (vm.count(BIN_MSR_FILE))
		p.i.bms_file = formPath<std::string>(p.g.output_folder, p.i.bms_file);
	else
		p.i.bms_file = formPath<std::string>(p.g.output_folder, p.g.network_name, "bms");

	if (vm.count(IMPORT_GEO_FILE))
		p.i.import_geo_file = 1;

	// output files
	p.i.asl_file = formPath<std::string>(p.g.output_folder, p.g.network_name, "asl");	// associated stations list
	p.i.aml_file = formPath<std::string>(p.g.output_folder, p.g.network_name, "aml");	// associated measurements list
	p.i.map_file = formPath<std::string>(p.g.output_folder, p.g.network_name, "map");	// station names map
	p.i.dst_file = formPath<std::string>(p.g.output_folder, p.g.network_name, "dst");	// duplicate stations
	p.i.dms_file = formPath<std::string>(p.g.output_folder, p.g.network_name, "dms");	// duplicate measurements
	p.i.imp_file = formPath<std::string>(p.g.output_folder, p.g.network_name, "imp");	// log
	
	//////////////////////////////////////////////////////////////////////////////
	// Ref frame options
	if (vm.count(OVERRIDE_INPUT_FRAME))
		p.i.override_input_rfame = 1;
	if (vm.count(REFERENCE_FRAME))
	{
		p.i.reference_frame = str_upper<std::string>(p.i.reference_frame);
		p.i.user_supplied_frame = 1;
	}

	if (vm.count(EPOCH))
	{
		// Get today's date?
		if (iequals(p.i.epoch, "today"))
			p.i.epoch = stringFromToday<boost::gregorian::date>();
		// Has the user supplied the year only?
		else if (p.i.epoch.rfind(".") == std::string::npos)
			p.i.epoch.insert(0, "01.01.");
		
		if (p.i.epoch.length() < 10)
		{
			std::string dateStr = FormatDateString(p.i.epoch);			
			if (dateStr.empty())
			{
				std::cout << std::endl << "- Error: Cannot parse date '" << p.i.epoch << "'." << std::endl <<
					"  Please supply date in the format dd.mm.yyyy" << std::endl << std::endl;
				return EXIT_FAILURE;
			}
			p.i.epoch = dateStr;
		}

		p.i.user_supplied_epoch = 1;
	}	

	//////////////////////////////////////////////////////////////////////////////
	// Data screening options
	if (vm.count(GET_MSRS_TRANSCENDING_BOX))
		p.i.include_transcending_msrs = 1;

	if (vm.count(SPLIT_CLUSTERS))
		p.i.split_clusters = 1;

	if (vm.count(IMPORT_SEG_BLOCK) && vm.count(IMPORT_CONTIG_NET))
	{
		std::cout << std::endl << "- Error: Cannot import stations and measurements using both options" << std::endl << "  --" << IMPORT_SEG_BLOCK << " and --" << IMPORT_CONTIG_NET << "." << std::endl <<
			"  Please supply only one option." << std::endl << std::endl;
		return EXIT_FAILURE;
	}

	// Import stations and measurements from a block
	if (vm.count(IMPORT_SEG_BLOCK))
	{
		p.i.import_block = 1;
		if (p.i.import_block_number < 1)
			p.i.import_block_number = 1;
	}

	// Import stations and measurements from a contiguous network
	if (vm.count(IMPORT_CONTIG_NET))
		p.i.import_network = 1;
	
	// User supplied segmentation file
	if (vm.count(SEG_FILE))
	{
		// Does it exist?
		if (!std::filesystem::exists(p.i.seg_file))
			// Look for it in the input folder
			p.i.seg_file = formPath<std::string>(p.g.input_folder, leafStr<std::string>(p.i.seg_file));
	}

	// convert all single X measurements to G measurements?
	if (vm.count(PREFER_X_MSR_AS_G))
		p.i.prefer_single_x_as_g = 1;

	if (vm.count(TEST_NEARBY_STNS))
		p.i.search_nearby_stn = 1;

	if (vm.count(TEST_SIMILAR_MSRS))
		p.i.search_similar_msr = 1;

	if (vm.count(TEST_SIMILAR_GNSS_MSRS))
		p.i.search_similar_msr_gx = 1;

	if (vm.count(IGNORE_SIMILAR_MSRS))
		p.i.ignore_similar_msr = 1;

	if (vm.count(IGNORE_INSUFFICIENT_MSRS))
		p.i.ignore_insufficient_msrs = 1;

	if (vm.count(REMOVE_IGNORED_MSRS))
		p.i.remove_ignored_msr = 1;

	if (vm.count(FLAG_UNUSED_STNS))
		p.i.flag_unused_stn = 1;

	if (vm.count(TEST_INTEGRITY))
		p.i.test_integrity = 1;
	
	//////////////////////////////////////////////////////////////////////////////
	// GNSS scaling
	if (vm.count(VSCALE) ||
		vm.count(PSCALE) ||
		vm.count(LSCALE) ||
		vm.count(HSCALE) ||
		vm.count(SCALAR_FILE))
		p.i.apply_scaling = 1;

	if (vm.count(SCALAR_FILE))
		if (!std::filesystem::exists(p.i.scalar_file))			// does it exist?
			// No.  Assume it is a filename contained in the input folder.  import will throw
			// an exception if it cannot be found.
			p.i.scalar_file = formPath<std::string>(p.g.input_folder, p.i.scalar_file);	
	
	//////////////////////////////////////////////////////////////////////////////
	// Export options
	
	if (vm.count(OUTPUT_MSR_TO_STN))
		p.o._msr_to_stn = 1;

	p.o._m2s_file = formPath<std::string>(p.g.output_folder, p.g.network_name, "m2s");	// measurement to stations table

	// Create file name based on the provided block number or contiguous network number
	std::string fileName(p.g.network_name);
	if (p.i.import_block)
	{
		std::stringstream blk("");
		blk << ".block-" << p.i.import_block_number;
		fileName += blk.str();
	}
	else if (p.i.import_network)
	{
		std::stringstream net("");
		net << ".network-" << p.i.import_network_number;
		fileName += net.str();
	}

	// Export to dynaml?
	if (vm.count(EXPORT_XML_FILES))
	{
		p.i.export_dynaml = 1;
	
		//// Export from binary
		//if (vm.count(EXPORT_FROM_BINARY))
		//	p.i.export_from_bfiles = 1;

		// single file for both stations and measurements
		if (vm.count(EXPORT_SINGLE_XML_FILE))
		{
			p.i.export_single_xml_file = 1;
			p.i.xml_outfile = formPath<std::string>(p.g.output_folder, 
				fileName, "xml");
		}
		// unique files for stations and measurements
		else
		{
			p.i.export_single_xml_file = 0;
			p.i.xml_stnfile = formPath<std::string>(p.g.output_folder, 
				fileName + "stn", "xml");
			p.i.xml_msrfile = formPath<std::string>(p.g.output_folder, 
				fileName + "msr", "xml");
		}
	}

	// Export dna files
	if (vm.count(EXPORT_DNA_FILES))
	{
		p.i.export_dna_files = 1;		
		p.i.dna_stnfile = formPath<std::string>(p.g.output_folder, 
			fileName, "stn");
		p.i.dna_msrfile = formPath<std::string>(p.g.output_folder, 
			fileName, "msr");
	}

	if (vm.count(EXPORT_ASL_FILE))
		p.i.export_asl_file = 1;
	if (vm.count(EXPORT_AML_FILE))
		p.i.export_aml_file = 1;
	if (vm.count(EXPORT_MAP_FILE))
		p.i.export_map_file = 1;
	if (vm.count(EXPORT_DISCONT_FILE))
		p.i.export_discont_file = 1;
	
	// Simulate
	if (vm.count(SIMULATE_MSR_FILE))
	{
		p.i.simulate_measurements = 1;
		p.i.simulate_msrfile = formPath<std::string>(p.g.output_folder,
			p.g.network_name, "simulated.msr");
	}

	// Station renaming file
	if (vm.count(STATION_RENAMING_FILE))
		p.i.rename_stations = 1;

	return EXIT_SUCCESS;
}

void PrintMeasurementstoStations(MsrTally* parsemsrTally, dna_import* parserDynaML, 
	project_settings* p, pvASLPtr vAssocStnList)
{
	// Print measurements to stations table
	if (p->o._msr_to_stn)
	{
		if (!p->g.quiet)
		{
			std::cout << "+ Printing summary of measurements connected to each station...";
			std::cout.flush();
		}
		parserDynaML->PrintMeasurementsToStations(p->o._m2s_file, parsemsrTally, 
			p->i.bst_file, p->i.bms_file, p->i.aml_file, vAssocStnList);
		if (!p->g.quiet)
			std::cout << " done." << std::endl;
	}
}

int SearchForSimilarMeasurements(dna_import* parserDynaML, project_settings* p, std::ofstream* imp_file,
	vdnaMsrPtr* vmeasurementsTotal)
{
	std::ofstream dms_file;
	UINT32 msr;
	vdnaMsrPtr vSimilarMeasurements;
	std::string comment("");

	try {
		if (!p->g.quiet)
		{
			std::cout << "+ Searching for similar measurements... ";
			std::cout.flush();
		}
		*imp_file << "+ Searching for similar measurements... ";

		// at this stage, perform only one search option, giving preference to GNSS searching 
		// if both are provided
		if (p->i.search_similar_msr_gx)
			msr = parserDynaML->FindSimilarGXMeasurements(vmeasurementsTotal, &vSimilarMeasurements);
		else
			msr = parserDynaML->FindSimilarMeasurements(vmeasurementsTotal, &vSimilarMeasurements);

		if (!vSimilarMeasurements.empty())
		{
			std::ostringstream ss_msg;
			if (!p->i.ignore_similar_msr)
				ss_msg << std::endl << "- Warning: ";
			else
				ss_msg << std::endl << "- Note: ";

			ss_msg << msr << " measurement" << (msr > 1 ? "s were" : " was") << 
				" found to be very similar (if not identical)\n  to " << 
				(msr > 1 ? "other measurements" : "another measurement") << 
				".  See " << p->i.dms_file << " for details." << std::endl;

			*imp_file << ss_msg.str();
			if (!p->g.quiet)
				std::cout << ss_msg.str();
			ss_msg.str("");

			try {
				// Create duplicate measurements file
				file_opener(dms_file, p->i.dms_file);
			}
			catch (const std::ios_base::failure& f) {
				std::stringstream ss;
				ss << "- Error: Could not open " << p->i.dms_file << "." << std::endl;
				ss << "  Check that the file exists and that the file is not already opened." << std::endl << f.what();
				if (!p->g.quiet)
					std::cout << ss.str();
				*imp_file << ss.str();
				imp_file->close();
				return EXIT_FAILURE;
			}

			PrintOutputFileHeaderInfo(&dms_file, p->i.dms_file, p, "DUPLICATE MEASUREMENTS FILE");

			// output message
			dms_file << std::endl << "- " << msr << " measurement" << (msr > 1 ? "s were" : " was") << 
				" found to be very similar (if not identical)\n  to " << 
				(msr > 1 ? "other measurements." : "another measurement.  ");

			if (p->i.ignore_similar_msr)
				dms_file << std::endl << "+ These measurements have been ignored.";
			dms_file << std::endl << std::endl;

			// dump measurements to dms file
			for_each (vSimilarMeasurements.begin(), vSimilarMeasurements.end(),
				[&dms_file, &comment] (dnaMsrPtr& m) {
					m->WriteDynaMLMsr(&dms_file, comment);
			});

			dms_file.close();

			if (p->i.ignore_similar_msr)
				ss_msg << "  These measurements have been ignored." << std::endl;
			else
				ss_msg << std::endl <<
				"  If the listed measurements are true duplicates, either remove each duplicate " << std::endl <<
				"  from the measurement file and re-run " << __BINARY_NAME__ << ", or re-run " << __BINARY_NAME__ << " with the" << std::endl <<
				"  --" << IGNORE_SIMILAR_MSRS << " option.  Alternatively, if each measurement " << std::endl <<
				"  is unique, then call " << __BINARY_NAME__ << " without the --" << TEST_SIMILAR_MSRS << " option." << std::endl << std::endl;

			if (!p->g.quiet)
				std::cout << ss_msg.str();
			*imp_file << ss_msg.str();

			if (!p->i.ignore_similar_msr)
			{
				//imp_file->close();
				return EXIT_SUCCESS;
			}
		}
		else
		{
			*imp_file << "Done. ";
			if (!p->g.quiet)
			{
				std::cout << "Done. ";
				std::cout.flush();
			}
		}

		*imp_file << std::endl;
		if (!p->g.quiet)
		{
			std::cout << std::endl;
			std::cout.flush();
		}

	} 
	catch (const XMLInteropException& e) {
		std::cout << std::endl << "- Error: " << e.what() << std::endl;
		*imp_file << std::endl << "- Error: " << e.what() << std::endl;
		imp_file->close();
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;

}

void ReportExportedData(const project_settings& p, std::ostream& out, const size_t& file_type,
	const UINT32& stnCount, const UINT32& msrCount)
{
	std::string dataExported;
	if (stnCount > 0 && msrCount == 0)
		dataExported = "stations";
	else if (stnCount == 0 && msrCount > 0)
		dataExported = "measurements";
	else if (stnCount > 0 && msrCount > 0)
		dataExported = "stations and measurements";

	std::string stn_file, msr_file, single_file;
	bool print_single_xml = (p.i.export_single_xml_file && p.i.export_dynaml);

	switch (file_type)
	{
	case dna:
		stn_file =  leafStr<std::string>(p.i.dna_stnfile);
		msr_file =  leafStr<std::string>(p.i.dna_msrfile);
		break;
	case dynaml:
		if (print_single_xml)
			single_file = leafStr<std::string>(p.i.xml_outfile);
		else
		{
			stn_file = leafStr<std::string>(p.i.xml_stnfile);
			msr_file = leafStr<std::string>(p.i.xml_msrfile);
		}
		break;
	default:
		break;
	}

	out << "+ Exporting " << dataExported << " to:" << std::endl;
	if (file_type == dynaml && print_single_xml)
		out << "  - " << single_file << "... ";
	else
	{
		if (stnCount > 0)
		{
			out << "  - " << stn_file;
			if (msrCount > 0)
				out << std::endl;
		}
		if (msrCount > 0)
			out << "  - " << msr_file;
		out << "... ";
	}
	out.flush();
}
	

void ExportStationsandMeasurements(dna_import* parserDynaML, const project_settings& p, std::ofstream* imp_file, vifm_t* vinput_file_meta,
	vdnaStnPtr* vstationsTotal, vdnaMsrPtr* vmeasurementsTotal, const UINT32& stnCount, const UINT32& msrCount)
{
	std::stringstream ssEpsgWarning;
	bool displayEpsgWarning(false);
	std::string epsgCode(epsgStringFromName<std::string>(p.i.reference_frame));

	// Check inconsistent reference frames
	if ((p.i.export_dynaml || p.i.export_dna_files) && !p.i.override_input_rfame)
	{
		for (UINT32 i(0); i<vinput_file_meta->size(); ++i)
		{
			if (!iequals(epsgCode, vinput_file_meta->at(i).epsgCode))
			{
				std::string inputFrame(datumFromEpsgString<std::string>(vinput_file_meta->at(i).epsgCode));
				ssEpsgWarning << std::endl << "- Warning: The default reference frame (" << p.i.reference_frame  << ")" << 
					" used for all exported" << std::endl <<
					"  files does not match the reference frame of one or more input files." << std::endl <<
					"  To suppress this warning, override the default reference frame using" << std::endl <<
					"  --reference-frame, or provide --override-input-ref-frame.";
				displayEpsgWarning = true;
				break;
			}
		}
	}

	// Sort on original file order
	parserDynaML->SortStationsForExport(vstationsTotal);

	// DynaML file format
	if (p.i.export_dynaml && (stnCount > 0 || msrCount > 0)) 
	{
		if (!p.g.quiet)
			ReportExportedData(p, std::cout, dynaml, stnCount, msrCount);
		ReportExportedData(p, *imp_file, dynaml, stnCount, msrCount);
		
		if (p.i.export_single_xml_file)
		{
			// Single output file
			parserDynaML->SerialiseDynaMLfromMemory(
				vstationsTotal, vmeasurementsTotal, 
				p.i.xml_outfile, p, vinput_file_meta, (p.i.flag_unused_stn ? true : false));
		}
		else
		{
			// Separate output files (default)
			parserDynaML->SerialiseDynaMLSepfromMemory(
				vstationsTotal, vmeasurementsTotal,
				p.i.xml_stnfile, p.i.xml_msrfile, p, vinput_file_meta, (p.i.flag_unused_stn ? true : false));
		}
		if (!p.g.quiet)
		{
			std::cout << "Done." << std::endl;
			std::cout.flush();
		}
		*imp_file << "Done." << std::endl;		
	}

	// DNA file format
	if (p.i.export_dna_files && (stnCount > 0 || msrCount > 0)) 
	{
		// Separate output files (default)
		if (!p.g.quiet)
			ReportExportedData(p, std::cout, dna, stnCount, msrCount);
		ReportExportedData(p, *imp_file, dna, stnCount, msrCount);

		parserDynaML->SerialiseDNA(
			vstationsTotal, vmeasurementsTotal,
			p.i.dna_stnfile, p.i.dna_msrfile, p, vinput_file_meta, (p.i.flag_unused_stn ? true : false));

		if (!p.g.quiet)
			std::cout << "Done." << std::endl;
		*imp_file << "Done." << std::endl;
	}

	if (displayEpsgWarning)
	{
		std::cout << ssEpsgWarning.str() << std::endl;
		*imp_file << ssEpsgWarning.str() << std::endl;
		std::cout.flush();
	}
}

int PrepareImportSegmentedData(project_settings& p, bool& userSuppliedSegFile)
{
	// Form default seg file path
	userSuppliedSegFile = false;

	// Has the user provided a segmentation file?
	if (!p.i.seg_file.empty())
		userSuppliedSegFile = true;
	else
		p.i.seg_file = formPath<std::string>(p.g.input_folder, p.g.network_name, "seg");

	if (!std::filesystem::exists(p.i.seg_file))
	{
		std::cout << std::endl << "- Error: The required segmentation file does not exist:" << std::endl <<
			"         " << p.i.seg_file << std::endl << std::endl <<
			"  Run  'segment " << p.g.network_name << "' to create a segmentation file" << std::endl << std::endl;
		return EXIT_FAILURE;
	}

	if (!std::filesystem::exists(p.i.bst_file))
	{
		std::cout << std::endl << "- Error: The required binary station file does not exist:" << std::endl <<
			"         " << p.i.bst_file << std::endl << std::endl;
		return EXIT_FAILURE;
	}

	if (!std::filesystem::exists(p.i.bms_file))
	{
		std::cout << std::endl << "- Error: The required binary measurement file does not exist:" << std::endl <<
			"         " << p.i.bms_file << std::endl << std::endl;
		return EXIT_FAILURE;
	}

	// If the user has not provided a seg file, check the meta of the default file
	if (!userSuppliedSegFile)
	{
		if (std::filesystem::last_write_time(p.i.seg_file) < std::filesystem::last_write_time(p.i.bst_file) ||
			std::filesystem::last_write_time(p.i.seg_file) < std::filesystem::last_write_time(p.i.bms_file))
		{
			// Has import been run after the segmentation file was created?
			binary_file_meta_t bst_meta, bms_meta;
			BstFile bst;
			BmsFile bms;
			bst.LoadFileMeta(p.i.bst_file, bst_meta);
			bms.LoadFileMeta(p.i.bms_file, bms_meta);

			bool bst_meta_import(iequals(bst_meta.modifiedBy, __import_app_name__) ||
				iequals(bst_meta.modifiedBy, __import_dll_name__));
			bool bms_meta_import(iequals(bms_meta.modifiedBy, __import_app_name__) ||
				iequals(bms_meta.modifiedBy, __import_dll_name__));

			if ((bst_meta_import && (std::filesystem::last_write_time(p.i.seg_file) < std::filesystem::last_write_time(p.i.bst_file))) ||
				(bms_meta_import && (std::filesystem::last_write_time(p.i.seg_file) < std::filesystem::last_write_time(p.i.bms_file))))
			{

				std::cout << std::endl << std::endl <<
					"- Error: The raw stations and measurements have been imported after" << std::endl <<
					"  the segmentation file was created:" << std::endl;

				time_t t_bst = file_time_to_time_t(std::filesystem::last_write_time(p.i.bst_file));
				time_t t_bms = file_time_to_time_t(std::filesystem::last_write_time(p.i.bms_file));
				time_t t_seg = file_time_to_time_t(std::filesystem::last_write_time(p.i.seg_file));

				std::cout << "   " << leafStr<std::string>(p.i.bst_file) << "  last modified on  " << ctime(&t_bst);
				std::cout << "   " << leafStr<std::string>(p.i.bms_file) << "  last modified on  " << ctime(&t_bms) << std::endl;
				std::cout << "   " << leafStr<std::string>(p.i.seg_file) << "  created on  " << ctime(&t_seg) << std::endl;
				std::cout << "  Run 'segment " << p.g.network_name << " [options]' to re-create the segmentation file, or re-run" << std::endl <<
					"  the import using the " << SEG_FILE << " option if this segmentation file must\n  be used." << std::endl << std::endl;
				return EXIT_FAILURE;
			}
		}
	}

	return EXIT_SUCCESS;
}

int ImportSegmentedBlock(dna_import& parserDynaML, vdnaStnPtr* vStations, vdnaMsrPtr* vMeasurements, 
		StnTally* parsestnTally, MsrTally* parsemsrTally, project_settings& p)
{	
	// Form default seg file path
	bool userSuppliedSegFile(false);

	if (PrepareImportSegmentedData(p, userSuppliedSegFile) != EXIT_SUCCESS)
		return EXIT_FAILURE;

	// Import stations and measurements from a particular block
	if (!p.g.quiet)
		std::cout << std::endl << "+ Importing stations and measurements from block " << p.i.import_block_number << " of\n  " << p.i.seg_file << "... ";

	try {
		parserDynaML.ImportStnsMsrsFromBlock(vStations, vMeasurements, p);
	}
	catch (const XMLInteropException& e) {
		std::stringstream ss;
		ss << std::endl << std::endl << "- Error: " << e.what();
		std::cout << ss.str() << std::endl;
		return EXIT_FAILURE;
	}

	*parsestnTally += parserDynaML.GetStnTally();
	*parsemsrTally += parserDynaML.GetMsrTally();
	if (!p.g.quiet)
		std::cout << "Done. " << std::endl;

	// Restore seg_file to null
	if (!userSuppliedSegFile)
		p.i.seg_file = "";

	return EXIT_SUCCESS;
}

int ImportContiguousNetwork(dna_import& parserDynaML, vdnaStnPtr* vStations, vdnaMsrPtr* vMeasurements,
	StnTally* parsestnTally, MsrTally* parsemsrTally, project_settings& p)
{
	// Form default seg file path
	bool userSuppliedSegFile(false);

	if (PrepareImportSegmentedData(p, userSuppliedSegFile) != EXIT_SUCCESS)
		return EXIT_FAILURE;

	// Import stations and measurements from a particular block
	if (!p.g.quiet)
		std::cout << std::endl << "+ Importing stations and measurements from contiguous network " << p.i.import_network_number << " of\n  " << p.i.seg_file << "... ";

	try {
		parserDynaML.ImportStnsMsrsFromNetwork(vStations, vMeasurements, p);
	}
	catch (const XMLInteropException& e) {
		std::stringstream ss;
		ss << std::endl << std::endl << "- Error: " << e.what();
		std::cout << ss.str() << std::endl;
		return EXIT_FAILURE;
	}

	*parsestnTally += parserDynaML.GetStnTally();
	*parsemsrTally += parserDynaML.GetMsrTally();
	if (!p.g.quiet)
		std::cout << "Done. " << std::endl;

	// Restore seg_file to null
	if (!userSuppliedSegFile)
		p.i.seg_file = "";

	return EXIT_SUCCESS;
}
	

int ImportDataFiles(dna_import& parserDynaML, vdnaStnPtr* vStations, vdnaMsrPtr* vMeasurements,
		vdnaStnPtr* vstationsTotal, vdnaMsrPtr* vmeasurementsTotal,
		std::ofstream* imp_file, vifm_t* vinput_file_meta, StnTally* parsestnTally, MsrTally* parsemsrTally, 
		UINT32& errorCount, project_settings& p)
{
	// For consideration:  
	//	- All input files could be read concurrently using multi-thread for faster input
	//	- The number of files that can be processed concurrently should be restricted to
	//	  the number of available cores
	//	- If multi-threading is to be used, then the progress thread will need to be redesigned.
	//	  Currently, the progress of each file (as a percentage) is written to cout as each file
	//	  is processed. If files are processed concurrently, progress reporting will need to be
	//    redesigned if multiple progress threads are writing to cout.
	//	- There is significant gain in using multi-thread.
	//	- The catch will be the need to carefully design one progress thread that captures progress
	//	  of all other "read" threads (rather than one progress for each), and somehow write to 
	//	  cout the respective progress.  For example, loading four files could be displayed as 
	//	  follows:
	//
	//	  + Parsing file1.xml (32%), file2.xml (16%), file3.xml (86%), file4.xml (56%)
	//	    Done: Loaded 190 stations from file3.xml in 00.032s
	//	  + Parsing file1.xml (72%), file2.xml (48%), file4.xml (91%)
	//	    Done: Loaded 236 stations from file4.xml in 00.050s
	//	  + Parsing file1.xml (87%), file2.xml (76%)
	//	    Done: Loaded 1596 measurements from file1.xml in 01.380s
	//	    Done: Loaded 3142 measurements from file2.xml in 02.034s
	//
	UINT32 stnCount(0), msrCount(0), clusterID(0);
	
	size_t pos = std::string::npos;
	size_t strlen_arg = 0;
	for_each(p.i.input_files.begin(), p.i.input_files.end(),
		[&strlen_arg](std::string& file) {
			if (leafStr<std::string>(file).length() > strlen_arg)
				strlen_arg = leafStr<std::string>(file).length();
	});

	strlen_arg += (6 + PROGRESS_PERCENT_04);

	size_t i, nfiles(p.i.input_files.size());		// for each file...
	std::string input_file, ss, status_msg;
	std::ostringstream ss_time, ss_msg;
	input_file_meta_t input_file_meta;
	boost::posix_time::milliseconds elapsed_time(boost::posix_time::milliseconds(0));
	boost::posix_time::ptime pt;

	if (!p.g.quiet)
		std::cout << "+ Parsing: " << std::endl;
	*imp_file << "+ Parsing " << std::endl;

	bool firstFile;

	// obtain the (default) project reference frame epsg code
	std::string projectEpsgCode(epsgStringFromName<std::string>(p.i.reference_frame));

	for (i=0; i<nfiles; i++)
	{
		stnCount = msrCount = 0;
		input_file = p.i.input_files.at(i);
		if (!std::filesystem::exists(input_file))
		{
			input_file = formPath<std::string>(p.g.input_folder, input_file);
			if (!std::filesystem::exists(input_file))
			{	
				std::cout << "- Error:  " << input_file << " does not exist" << std::endl;
				return EXIT_FAILURE;
			}
		}
			
		ss = leafStr<std::string>(p.i.input_files.at(i)) + "... ";
		if (!p.g.quiet)
			std::cout << "  " << std::setw(strlen_arg) << std::left << ss;
		*imp_file << "  " << std::setw(strlen_arg) << std::left << ss;

		running = true;
		firstFile = bool(i == 0);

		std::vector<std::thread> ui_interop_threads;
		if (!p.g.quiet)
			ui_interop_threads.emplace_back(dna_import_progress_thread(&parserDynaML, &p));
		ui_interop_threads.emplace_back(dna_import_thread(&parserDynaML, &p, input_file,
			vStations, &stnCount, vMeasurements, &msrCount,
			&clusterID, &input_file_meta, firstFile, &status_msg,
			&elapsed_time));
		for (auto& t : ui_interop_threads) {
			t.join();
		}

		switch (parserDynaML.GetStatus())
		{
		case PARSE_EXCEPTION_RAISED:
			*imp_file << std::endl << status_msg;
			running = false;
			return EXIT_FAILURE;
			break;
		case PARSE_UNRECOGNISED_FILE:
			*imp_file << status_msg << std::endl;
			errorCount++;
			continue;
		case PARSE_SUCCESS:
			running = false;
			break;
		default:
			errorCount++;
			std::cout << std::endl;
			continue;
		}

		vinput_file_meta->push_back(input_file_meta);

		ss_time.str("");

		if (stnCount > 0 || msrCount > 0) // stations or measurements only
		{			
			ss_time << "  Done. Loaded ";
			if (stnCount)
				ss_time << stnCount << " stations";
			if (stnCount && msrCount)
				ss_time << " and ";
			if (msrCount)
				ss_time << msrCount << " measurements";
			ss_time << " in ";
				
			pt = boost::posix_time::ptime(boost::gregorian::day_clock::local_day(), elapsed_time);
			if (elapsed_time < boost::posix_time::seconds(3))
			{
				boost::posix_time::time_facet* facet(new boost::posix_time::time_facet("%s"));
				ss_time.imbue(std::locale(ss_time.getloc(), facet));
				ss_time << pt << "s";			
			}
			else if (elapsed_time < boost::posix_time::seconds(61))
			{		
				boost::posix_time::time_facet* facet(new boost::posix_time::time_facet("%S"));
				ss_time.imbue(std::locale(ss_time.getloc(), facet));
				ss_time << pt << "s";
			}
			else
				ss_time << elapsed_time;

			std::string time_message = ss_time.str();
			while ((pos = time_message.find("0s")) != std::string::npos)
				time_message = time_message.substr(0, pos) + "s";

			if ((pos = time_message.find(" 00.")) != std::string::npos)
				time_message = time_message.replace(pos, 4, " 0.");
			if ((pos = time_message.find(" 0.s")) != std::string::npos)
				time_message = time_message.replace(pos, 4, " 0s");

			if (!p.g.quiet)
			{
				if (isatty(fileno(stdout)))
					std::cout << PROGRESS_BACKSPACE_04;
				std::cout << time_message << std::endl;
			}
			*imp_file << time_message << std::endl;

			// Capture the input file's default reference frame and set the
			// project reference frame (if not specified on the command-line)
			std::string inputFileEpsg(""), inputFileDatum(""), inputFileEpoch("");
			UINT32 inputFileEpsgi;
			try {
				inputFileDatum = datumFromEpsgString<std::string>(input_file_meta.epsgCode);
				inputFileEpsg = input_file_meta.epsgCode;
				inputFileEpsgi = LongFromString<UINT32>(inputFileEpsg);
				inputFileEpoch = input_file_meta.epoch;
			}
			catch (...) 
			{
				// Do nothing, revert to defaults
				// p.i.reference_frame
				// p.i.epoch
			}

			// Is this the first file?  If so, attempt to set the default datum from
			// the input file (if applicable).
			
			bool referenceframeChanged(false);
			bool epochChanged(false);
			
			if (firstFile)
			{
				bool isStaticFrame(false);

				// If the user has not provided a reference frame, then inspect the file reference frame
				// If the file does not contain a reference frame (e.g. SNX) or the user has left it blank, 
				// fileEpsg will be empty. Thus, no change will occur (retain the default GDA2020).
				if (!inputFileEpsg.empty() && !p.i.user_supplied_frame)
				{
					// Set the project defaults
					referenceframeChanged = true;
					
					p.i.reference_frame = inputFileDatum;
					p.r.reference_frame = inputFileDatum;
					projectEpsgCode = epsgStringFromName<std::string>(p.i.reference_frame);
				}

				isStaticFrame = isEpsgDatumStatic(LongFromString<UINT32>(projectEpsgCode));
				
				// Has the user supplied a static reference frame?  
				// If so, set the epoch from the user-supplied static frame
				if (p.i.user_supplied_frame && isStaticFrame)
				{
					p.i.epoch = referenceepochFromEpsgString<std::string>(projectEpsgCode);
					p.r.epoch = p.i.epoch;
				}
				// Has the user supplied a SINEX file, but no reference frame?
				// If so, set the epoch from the DynAdjust default datum
				else if (!p.i.user_supplied_frame && input_file_meta.filetype == sinex)
				{
					// revert to epoch of the file frame
					p.i.epoch = referenceepochFromEpsgString<std::string>(projectEpsgCode);
					p.r.epoch = p.i.epoch;
				}
				// At this point, one of the following is true:
				//  - a dynamic reference frame has been supplied
				//  - a frame has not been supplied, in which case the file frame will be takenelse if (!p.i.user_supplied_epoch)
				else if (!p.i.user_supplied_epoch)
				{
					if (inputFileEpoch.empty())
					{						
						// revert to epoch of the file frame
						p.i.epoch = referenceepochFromEpsgCode<UINT32>(inputFileEpsgi);						
						p.r.epoch = p.i.epoch;
					}
					else
					{
						p.i.epoch = inputFileEpoch;
						p.r.epoch = p.i.epoch;
					}
				}
							

				try {
					// Initialise the 'default' datum (frame and epoch) for the project, from the first file, unless the
					// frame and epoch have been set by the user, in which case InitialiseDatum has already initialised the datum.
					parserDynaML.InitialiseDatum(p.i.reference_frame, p.i.epoch);
				}
				catch (const XMLInteropException& e) {
					std::stringstream ss;
					ss << "- Error: ";
					std::cout << ss.str() << e.what() << std::endl;
					return EXIT_FAILURE;
				}

				if (p.i.import_block || p.i.import_network)
					continue;

				// When the user has not supplied a frame on the command line
				if (!p.i.user_supplied_frame)
				{
					// reference frame
					std::stringstream datumSource;
					switch (input_file_meta.filetype)
					{
					case sinex:
						datumSource << ". DynAdjust default (frame not present within SNX file)";
						break;
					default:
						if (parserDynaML.filespecifiedReferenceFrame())
							datumSource << ". Taken from " << FormatFileType<std::string>(input_file_meta.filetype) << " header.";
						else
							datumSource << ". Reference frame not supplied in " << FormatFileType<std::string>(input_file_meta.filetype) << " header.";
					}

					if (!p.g.quiet)
						std::cout << std::left << "  Project reference frame set to " << p.i.reference_frame << datumSource.str() << std::endl;
					*imp_file << std::left << "  Project reference frame set to " << p.i.reference_frame << datumSource.str() << std::endl;
				}				
				// When the user has supplied a frame on the command line, and the input file datum field is blank
				else if (!parserDynaML.filespecifiedReferenceFrame())
				{
					std::stringstream ssEpsgWarning;
					switch (input_file_meta.filetype)
					{
					case sinex:
						ssEpsgWarning << "  - Warning: SNX reference frame set to " << p.i.reference_frame << ".";
						break;
					default:
						ssEpsgWarning << "  - Warning: Input file reference frame not supplied. Adopting " << inputFileDatum << ".";
					}				

					if (!p.g.quiet)
						std::cout << ssEpsgWarning.str() << std::endl;
					*imp_file << ssEpsgWarning.str() << std::endl;
				}
				// when reference frame has been supplied and the value in the file doesn't match
				else if (p.i.user_supplied_frame)
				{
					if (!iequals(projectEpsgCode, input_file_meta.epsgCode))
					{
						std::stringstream ssEpsgWarning;
						ssEpsgWarning << "  - Warning: File reference frame (" << inputFileDatum << ") differs from project frame (" << p.i.reference_frame << ").";
						if (!p.g.quiet)
							std::cout << ssEpsgWarning.str() << std::endl;
						*imp_file << ssEpsgWarning.str() << std::endl;
					}
				}

				// When the user has not supplied an epoch on the command line
				if (!p.i.user_supplied_epoch)
				{
					// epoch
					std::stringstream epochSource;
					if (isStaticFrame)
					{
						epochSource << "  - Warning: File epoch (" << inputFileEpoch << ") differs from project epoch (" << p.i.epoch << ")." << std::endl;
						if (!p.g.quiet)
							std::cout << epochSource.str();
						*imp_file << epochSource.str();
					}
					else
					{
						if (parserDynaML.filespecifiedEpoch())
						{
							if (isEpsgDatumStatic(inputFileEpsgi))
								epochSource << " (adopted reference epoch of " << inputFileDatum << ").";
							else
								epochSource << ". Taken from " << FormatFileType<std::string>(input_file_meta.filetype) << " header.";
						}
						else
							epochSource << ". Epoch not supplied in " << FormatFileType<std::string>(input_file_meta.filetype) << " header.";
						if (!p.g.quiet)
							std::cout << std::left << "  Project epoch set to " << p.i.epoch << epochSource.str() << std::endl;
						*imp_file << std::left << "  Project epoch set to " << p.i.epoch << epochSource.str() << std::endl;
						epochChanged = true;
					}
				}				
				// When the user has supplied an epoch on the command line, and the input file epoch field is blank
				else if (!parserDynaML.filespecifiedEpoch())
				{
					std::stringstream ssEpochWarning;
					ssEpochWarning << "  - Warning: Input file epoch not supplied. Adopting " << p.i.epoch << ".";

					if (!p.g.quiet)
						std::cout << ssEpochWarning.str() << std::endl;
					*imp_file << ssEpochWarning.str() << std::endl;
				}
				// when epoch has been supplied and the value in the file doesn't match
				else if (p.i.user_supplied_epoch)
				{
					if (!iequals(p.i.epoch, input_file_meta.epoch))
					{
						std::stringstream ssEpochWarning;
						ssEpochWarning << "  - Warning: File epoch (" << inputFileEpoch << ") differs from project epoch (" << p.i.epoch << ").";
						if (!p.g.quiet)
							std::cout << ssEpochWarning.str() << std::endl;
						*imp_file << ssEpochWarning.str() << std::endl;
					}
				}

				if (referenceframeChanged || epochChanged)
				{
					if (!p.g.quiet)
						std::cout << std::endl;
					*imp_file << std::endl;
				}
			}
			else
			{
				if (p.i.import_block || p.i.import_network)
					continue;

				// Was the datum field empty in the file?
				if (!parserDynaML.filespecifiedReferenceFrame())
				{					
					std::stringstream ssEpsgWarning;
					ssEpsgWarning << "  - Warning: Input file reference frame not supplied. Adopting " << inputFileDatum << ".";
					if (!p.g.quiet)
						std::cout << ssEpsgWarning.str() << std::endl;
					*imp_file << ssEpsgWarning.str() << std::endl;					
				}
				// Is the datum in the file different to the project datum?
				else if (!iequals(projectEpsgCode, input_file_meta.epsgCode))
				{
					std::stringstream ssEpsgWarning;
					if (referenceframeChanged)
					{
						ssEpsgWarning << "  - Warning: The project reference frame has been set to the default" << std::endl <<
							"    file datum of " << leafStr<std::string>(p.i.input_files.at(i)) << " (" << inputFileDatum << ").";

						// set the project reference frame epsg code
						projectEpsgCode = epsgStringFromName<std::string>(p.i.reference_frame);
					}
					else
					{
						ssEpsgWarning << "  - Warning: File reference frame (" << inputFileDatum << ") differs from project frame (" << p.i.reference_frame << ").";
					}

					if (!p.g.quiet)
						std::cout << ssEpsgWarning.str() << std::endl;
					*imp_file << ssEpsgWarning.str() << std::endl;
				}

				// Was the epoch field empty in the file?
				if (!parserDynaML.filespecifiedEpoch())
				{
					std::stringstream ssEpochWarning;
					ssEpochWarning << "  - Warning: Input file epoch not supplied. Adopting " << p.i.epoch << ".";
					if (!p.g.quiet)
						std::cout << ssEpochWarning.str() << std::endl;
					*imp_file << ssEpochWarning.str() << std::endl;
				}
				// Is the epoch in the file different to the project epoch?
				else if (!iequals(p.i.epoch, input_file_meta.epoch))
				{
					std::stringstream ssEpochWarning;
					
					ssEpochWarning << "  - Warning: File epoch (" << inputFileEpoch << ") differs from project epoch (" << p.i.epoch << ").";

					if (!p.g.quiet)
						std::cout << ssEpochWarning.str() << std::endl;
					*imp_file << ssEpochWarning.str() << std::endl;
				}
			}
		}

		// Handle discontinuities for non-sinex files, if and only if a discontinuity file has been loaded.
		// NOTE: If a station is marked as a discontinuity site, and that same station is in the renaming
		// file, station renaming for that site will not occur.  This is because the site name will be changed
		// to <site-name>_yyyymmdd.
		// If renaming must still occur, then it might be possible to rename using the first part of the name.
		

		if (stnCount > 0) // stations only
		{
			//vstationsTotal.reserve(vstationsTotal.size() + stnCount);
			// combine stations and station tally
			vstationsTotal->insert(vstationsTotal->end(), vStations->begin(), vStations->end());
			*parsestnTally += parserDynaML.GetStnTally();
			vStations->clear();
		}
		if (msrCount > 0) // measurements only
		{		
			//vmeasurementsTotal.reserve(vmeasurementsTotal.size() + msrCount);
			// combine measurements
			vmeasurementsTotal->insert(vmeasurementsTotal->end(), vMeasurements->begin(), vMeasurements->end());
			*parsemsrTally += parserDynaML.GetMsrTally();
			vMeasurements->clear();
		}
	}

	if (errorCount == nfiles || (vstationsTotal->empty() && vmeasurementsTotal->empty()))
		return EXIT_FAILURE;

	running = false;
	return EXIT_SUCCESS;
}

int main(int argc, char* argv[])
{
	// create banner message
	std::string cmd_line_banner;
	fileproc_help_header(&cmd_line_banner);

	project_settings p;

	boost::program_options::variables_map vm;
	boost::program_options::positional_options_description positional_options;

	boost::program_options::options_description standard_options("+ " + std::string(ALL_MODULE_STDOPT), PROGRAM_OPTIONS_LINE_LENGTH);
	boost::program_options::options_description ref_frame_options("+ " + std::string(IMPORT_MODULE_FRAME), PROGRAM_OPTIONS_LINE_LENGTH);
	boost::program_options::options_description data_screening_options("+ " + std::string(IMPORT_MODULE_SCREEN), PROGRAM_OPTIONS_LINE_LENGTH);
	boost::program_options::options_description gnss_scaling_options("+ " + std::string(IMPORT_MODULE_GNSS_VAR), PROGRAM_OPTIONS_LINE_LENGTH);
	boost::program_options::options_description simulation_options("+ " + std::string(IMPORT_MODULE_SIMULATE), PROGRAM_OPTIONS_LINE_LENGTH);
	boost::program_options::options_description output_options("+ " + std::string(ALL_MODULE_OUTPUT), PROGRAM_OPTIONS_LINE_LENGTH);
	boost::program_options::options_description export_options("+ " + std::string(ALL_MODULE_EXPORT), PROGRAM_OPTIONS_LINE_LENGTH);
	boost::program_options::options_description generic_options("+ " + std::string(ALL_MODULE_GENERIC), PROGRAM_OPTIONS_LINE_LENGTH);

	std::string cmd_line_usage("+ ");
	cmd_line_usage.append(__BINARY_NAME__).append(" usage:  ").append(__BINARY_NAME__).append(" [options] [files]...");
	boost::program_options::options_description allowable_options(cmd_line_usage, PROGRAM_OPTIONS_LINE_LENGTH);

	try {
		// Declare a group of options that will be 
		// allowed only on command line		
		standard_options.add_options()
			(PROJECT_FILE_P, boost::program_options::value<std::string>(&p.g.project_file),
				"Project file name. Full path to project file. If none specified, a new file is created using input-folder and network-name.  If a project file exists, the file is used to set all command line options.")
			(NETWORK_NAME_N, boost::program_options::value<std::string>(&p.g.network_name), 
				"Network name. User defined name for all input and output files. Default is \"network#\", where # is incremented until the first available network is reached.")
			(IMPORT_FILE_F, boost::program_options::value< vstring >(&p.i.input_files), 
				"Station and measurement input file(s). Switch is not required.")
			(IMPORT_GEO_FILE_G, boost::program_options::value<std::string>(&p.i.geo_file),
				"Import DNA geoid file.")
			(INPUT_FOLDER_I, boost::program_options::value<std::string>(&p.g.input_folder),
				"Path containing all input files.")
			(OUTPUT_FOLDER_O, boost::program_options::value<std::string>(&p.g.output_folder),		// default is ./,
				"Path for all output files.")
			(BIN_STN_FILE_S, boost::program_options::value<std::string>(&p.i.bst_file),
				"Binary station output file name. Overrides network name.")
			(BIN_MSR_FILE_M, boost::program_options::value<std::string>(&p.i.bms_file),
				"Binary measurement output file name. Overrides network name.")
			;

		ref_frame_options.add_options()
			(REFERENCE_FRAME_R, boost::program_options::value<std::string>(&p.i.reference_frame), 
			(std::string("Project reference frame for all stations, measurements, and preliminary reductions on the ellipsoid when input files do not specify a reference frame. Default is ") +
				p.i.reference_frame + ".").c_str())
			(EPOCH_E, boost::program_options::value<std::string>(&p.i.epoch),
			(std::string("Project epoch for all stations and measurements when input files do not specify an epoch. Default is ") +
				p.i.epoch + ".").c_str())
			(OVERRIDE_INPUT_FRAME,
				(std::string("Replace the reference frame specified in the input files with the reference frame specified by arg in --") +
				std::string(REFERENCE_FRAME)).c_str())
			;

		data_screening_options.add_options()
			(BOUNDING_BOX, boost::program_options::value<std::string>(&p.i.bounding_box),
				"Import stations and measurements within bounding box. arg is a comma delimited string \"lat1,lon1,lat2,lon2\" (in dd.mmss) defining the upper-left and lower-right limits.")
			(GET_MSRS_TRANSCENDING_BOX,
				"Include measurements which transcend bounding box, including associated stations.")
			(INCLUDE_STN_ASSOC_MSRS, boost::program_options::value<std::string>(&p.i.stn_associated_msr_include),
				"Include stations and all associated measurements. arg is a comma delimited string \"stn 1,stn 2,stn 3,...,stn N\" of the stations to include.")
			(EXCLUDE_STN_ASSOC_MSRS, boost::program_options::value<std::string>(&p.i.stn_associated_msr_exclude),
				"Exclude stations and all associated measurements. arg is a comma delimited string \"stn 1,stn 2,stn 3,...,stn N\" of the stations to exclude.")
			(SPLIT_CLUSTERS,
				"Allow bounding-box or get-stns-and-assoc-msrs to split GNSS point and baseline cluster measurements.")
			(IMPORT_SEG_BLOCK, boost::program_options::value<UINT32>(&p.i.import_block_number),
				"Extract stations and measurements from this block.")
			(IMPORT_CONTIG_NET, boost::program_options::value<UINT32>(&p.i.import_network_number),
				"Extract stations and measurements from this contiguous network.")
			(SEG_FILE, boost::program_options::value<std::string>(&p.i.seg_file),
				"Network segmentation input file. Filename overrides network name.")
			(PREFER_X_MSR_AS_G,
				"Import single baseline cluster measurements (X) as single baselines (G).")
			(INCLUDE_MSRS, boost::program_options::value<std::string>(&p.i.include_msrs),
				"Import the specified measurement types. arg is a non-delimited string of measurement types (eg \"GXY\").")
			(EXCLUDE_MSRS, boost::program_options::value<std::string>(&p.i.exclude_msrs),
				"Exclude the specified measurement types. arg is a non-delimited string of measurement typs (eg \"IJK\").")
			(STATION_RENAMING_FILE, boost::program_options::value<std::string>(&p.i.stn_renamingfile),
				"Station renaming file")
			(STATION_DISCONTINUITY_FILE, boost::program_options::value<std::string>(&p.i.stn_discontinuityfile),
				"Station discontinuity file.  Applies discontinuity dates to station names in station and measurement files.")
			(TEST_NEARBY_STNS,
				"Search for nearby stations.")
			(TEST_NEARBY_STN_DIST, boost::program_options::value<double>(&p.i.search_stn_radius),
				(std::string("Specify the radius of the circle within which to search for nearby stations.  Default is ")+
				StringFromT(STN_SEARCH_RADIUS)+std::string("m")).c_str())
			(TEST_SIMILAR_GNSS_MSRS,
				"Search and provide warnings for GNSS baselines (G) and baseline clusters (X) which appear to have been derived from the same source data.")
			(TEST_SIMILAR_MSRS,
				"Search and provide warnings for similar measurements.")
			(IGNORE_INSUFFICIENT_MSRS,
				"Ignore measurements which do not sufficiently constrain a station in two dimensions.")
			(IGNORE_SIMILAR_MSRS,
				"Ignore similar measurements.")
			(REMOVE_IGNORED_MSRS,
				"Remove ignored measurements.")
			(FLAG_UNUSED_STNS,
				"Mark unused stations in binary file.  Stations marked will be excluded from any further processing.")
			(TEST_INTEGRITY,
				"Test the integrity of the association lists and binary files.")
			;

		gnss_scaling_options.add_options()
			(VSCALE, boost::program_options::value<double>(&p.i.vscale),
				(std::string("Global variance (v) matrix scalar for all GNSS measurements.  Replaces existing scalar.  Default is ")+
				StringFromT(p.i.vscale)+std::string(".")).c_str())
			(PSCALE, boost::program_options::value<double>(&p.i.pscale),
				(std::string("Latitude (p=phi) variance matrix scalar for all GNSS measurements.  Replaces existing scalar.  Default is ")+
				StringFromT(p.i.pscale)+std::string(".")).c_str())
			(LSCALE, boost::program_options::value<double>(&p.i.lscale),
				(std::string("Longitude (l=lambda) variance matrix scalar for all GNSS measurements.  Replaces existing scalar.  Default is ")+
				StringFromT(p.i.lscale)+std::string(".")).c_str())
			(HSCALE, boost::program_options::value<double>(&p.i.hscale),
				(std::string("Height (h) variance matrix scalar for all GNSS measurements.  Replaces existing scalar.  Default is ")+
				StringFromT(p.i.hscale)+std::string(".")).c_str())
			(SCALAR_FILE, boost::program_options::value<std::string>(&p.i.scalar_file),
				"File containing v, p, l and h scalars for GNSS baseline measurements between specific station pairs.  Scalar file values do not apply to GNSS point or baseline clusters.")
			;

		output_options.add_options()
			(OUTPUT_MSR_TO_STN,
				"Output summary of measurements connected to each station.")
			(OUTPUT_MSR_TO_STN_SORTBY, boost::program_options::value<UINT16>(&p.o._sort_msr_to_stn),
				std::string("Sort order for measurement to stations summary.\n  " +
					StringFromT(orig_stn_sort_ui) + ": Original station order (default)\n  " +
					StringFromT(meas_stn_sort_ui) + ": Measurement count").c_str())
			;

		export_options.add_options()
			(EXPORT_XML_FILES,
				"Export stations and measurements to DynaML (DynAdjust XML) format.")
			(EXPORT_SINGLE_XML_FILE,
				"Create a single DynaML file for stations and measurements.")
			// Nice, but somewhat redundant functionality that offers no
			// benefit to the user
			//(EXPORT_FROM_BINARY,
			//	"Create DynaML output file using binary files. Default option uses internal memory.")
			(EXPORT_DNA_FILES,
				"Export stations and measurements to DNA STN and MSR format.")
			(EXPORT_ASL_FILE,
				"Export the ASL file as raw text.")
			(EXPORT_AML_FILE,
				"Export the AML file as raw text.")
			(EXPORT_MAP_FILE,
				"Export the MAP file as raw text.")
			(EXPORT_DISCONT_FILE,
				"Export discontinuity information as raw text.")
			;

		simulation_options.add_options()
			(SIMULATE_MSR_FILE,
				"Simulate exact measurements corresponding to the input measurements using the coordinates in the station file. To apply geoid--ellipsoid separations and deflections of the vertical to the simulated measurements, introduce a geoid file using the --geo-file option.")
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

		allowable_options.add(standard_options).add(ref_frame_options).add(data_screening_options).add(gnss_scaling_options).add(simulation_options).add(output_options).add(export_options).add(generic_options);

		// add "positional options" to handle command line tokens which have no option name
		positional_options.add(IMPORT_FILE, -1);
		
		boost::program_options::command_line_parser parser(argc, argv);
		store(parser.options(allowable_options).positional(positional_options).run(), vm);
		notify(vm);
	}
	catch(const std::exception& e) 
	{
		std::cout << "- Error: " << e.what() << std::endl;
		std::cout << cmd_line_banner << allowable_options << std::endl;
		return EXIT_FAILURE;
	}
	catch (...) 
	{
		std::cout << "+ Exception of unknown type!\n";
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
		else if (str_upper<std::string, char>(IMPORT_MODULE_FRAME).find(help_text) != std::string::npos) {
			std::cout << ref_frame_options << std::endl;
		} 
		else if (str_upper<std::string, char>(IMPORT_MODULE_SCREEN).find(help_text) != std::string::npos) {
			std::cout << data_screening_options << std::endl;
		} 
		else if (str_upper<std::string, char>(IMPORT_MODULE_GNSS_VAR).find(help_text) != std::string::npos) {
			std::cout << gnss_scaling_options << std::endl;
		} 
		else if (str_upper<std::string, char>(IMPORT_MODULE_SIMULATE).find(help_text) != std::string::npos) {
			std::cout << simulation_options << std::endl;
		} 
		else if (str_upper<std::string, char>(ALL_MODULE_OUTPUT).find(help_text) != std::string::npos) {
			std::cout << output_options << std::endl;
		} 
		else if (str_upper<std::string, char>(ALL_MODULE_EXPORT).find(help_text) != std::string::npos) {
			std::cout << export_options << std::endl;
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
	if (!p.i.seg_file.empty())
		userSuppliedSegFile = true;
	bool userSuppliedBstFile(false);
	if (!p.i.bst_file.empty())
		userSuppliedBstFile = true;
	bool userSuppliedBmsFile(false);
	if (!p.i.bms_file.empty())
		userSuppliedBmsFile = true;

	if (ParseCommandLineOptions(argc, argv, vm, p) != EXIT_SUCCESS)
		return EXIT_FAILURE;

	UINT32 errorCount(0);
	bool stn_map_created = false, measurements_mapped = false;

	std::string input_file;
	vstring input_files;
	std::string status_msg;

	std::ofstream imp_file;
	try {
		// Create import log file.  Throws runtime_error on failure.
		file_opener(imp_file, p.i.imp_file);
	}
	catch (const std::runtime_error& e) {
		std::stringstream ss;
		ss << "- Error: Could not open " << p.i.imp_file << ". \n  Check that the file exists and that the file is not already opened." << std::endl;
		std::cout << ss.str() << e.what() << std::endl;
		return EXIT_FAILURE;
	}

	vdnaStnPtr vStations, vstationsTotal;
	vdnaMsrPtr vMeasurements, vmeasurementsTotal;
	vdnaMsrPtr::iterator _itm;
	vASLPtr associatedSL;
	vUINT32 associatedML;
	v_string_uint32_pair vStnsMap_sortName;	// Station Name Map sorted on name (string)
	
	vstationsTotal.clear();
	vmeasurementsTotal.clear();

	dna_import parserDynaML;

	////////////////////////////////////////////////////////////////////////////////////
	// First things first!
	// Set the 'default' reference frame for the binary station and measurement files
	// At this point, the frame may be the hard-coded default, or a user-specified 
	// frame via:   -r [--reference-frame] arg
	// See comments in InitialiseDatum()
	try {
		// Initialise the 'default' datum for the project.
		parserDynaML.InitialiseDatum(p.i.reference_frame, p.i.epoch);
	}
	catch (const XMLInteropException& e) {

		std::cout << std::endl << cmd_line_banner;
		imp_file << std::endl << cmd_line_banner;
		
		std::stringstream ss;
		ss << "- Error: ";
		std::cout << ss.str() << e.what() << std::endl;
		imp_file << ss.str() << e.what() << std::endl;
		return EXIT_FAILURE;
	}

	// obtain the project reference frame
	UINT32 epsgCode(epsgCodeFromName<UINT32>(p.i.reference_frame));

	// set the default output reference frame and epoch, so that
	// if adjust is called without reftran, it reflects the datum 
	// supplied on import
	p.r.reference_frame = p.i.reference_frame;

	// if the reference frame supplied is static, then force the epoch
	// to be the reference epoch
	if (isEpsgDatumStatic(epsgCode))
	{
		p.i.epoch = referenceepochFromEpsgCode<UINT32>(epsgCode);
		p.r.epoch = p.i.epoch;
	}
	else
	{
		// Dynamic reference frame.  For this code to be reached, which
		// is prior to the input files been parsed, the user must have
		// supplied a dynamic reference frame on the command line.
		// 
		// Has the user not supplied an epoch?
		if (!p.i.user_supplied_epoch)
			// No, so set the "default" epoch for the project to be the
			// reference epoch of the supplied dynamic reference frame
			p.i.epoch = referenceepochFromEpsgCode<UINT32>(epsgCode);
		
		// Set the project epoch for reftran also
		p.r.epoch = p.i.epoch;
	}

	////////////////////////////////////////////////////////////////////////////////////

	if (vm.count(QUIET))
		p.g.quiet = 1;

	if (!p.g.quiet)
	{
		std::cout << std::endl << cmd_line_banner;

		std::cout << "+ Options:" << std::endl; 
		std::cout << std::setw(PRINT_VAR_PAD) << std::left << "  Network name: " <<  p.g.network_name << std::endl;
		std::cout << std::setw(PRINT_VAR_PAD) << std::left << "  Input folder: " << p.g.input_folder << std::endl;
		std::cout << std::setw(PRINT_VAR_PAD) << std::left << "  Output folder: " << p.g.output_folder << std::endl;
		std::cout << std::setw(PRINT_VAR_PAD) << std::left << "  Associated station file: " << p.i.asl_file << std::endl;
		std::cout << std::setw(PRINT_VAR_PAD) << std::left << "  Associated measurement file: " << p.i.aml_file << std::endl;
		std::cout << std::setw(PRINT_VAR_PAD) << std::left << "  Binary station output file: " << p.i.bst_file << std::endl;
		std::cout << std::setw(PRINT_VAR_PAD) << std::left << "  Binary measurement output file: " << p.i.bms_file << std::endl;
		
		// If a reference frame and epoch have been supplied, report them.
		// If not, the assumption is, the project frame and epoch will be adopted from the first file and
		// in this case, it will be reported later
		if (p.i.user_supplied_frame)
			std::cout << std::setw(PRINT_VAR_PAD) << std::left << "  Project reference frame:" << p.i.reference_frame << " (user supplied)" << std::endl;
		else
			std::cout << std::setw(PRINT_VAR_PAD) << std::left << "  Project reference frame:" << "To be adopted from the first input file" << std::endl;
		
		if (p.i.override_input_rfame)
			std::cout << std::setw(PRINT_VAR_PAD) << std::left << "  Override input file ref frame:" << yesno_string(p.i.override_input_rfame) << std::endl;

		if (isEpsgDatumStatic(epsgCode) && p.i.user_supplied_frame)
			std::cout << std::setw(PRINT_VAR_PAD) << std::left << "  Project epoch:" << p.i.epoch << " (adopted reference epoch of " << p.i.reference_frame << ")" << std::endl;
		else
		{
			if (p.i.user_supplied_epoch)
				std::cout << std::setw(PRINT_VAR_PAD) << std::left << "  Project epoch:" << p.i.epoch << " (user supplied)" << std::endl;
			else
				std::cout << std::setw(PRINT_VAR_PAD) << std::left << "  Project epoch:" << "To be adopted from the first input file" << std::endl;
		}

		if (p.i.export_dynaml)
		{
			if (p.i.export_single_xml_file)
				std::cout << std::setw(PRINT_VAR_PAD) << std::left << "  DynaML output file: " << p.i.xml_outfile << std::endl;
			else
			{
				std::cout << std::setw(PRINT_VAR_PAD) << std::left << "  DynaML station file: " << p.i.xml_stnfile << std::endl;
				std::cout << std::setw(PRINT_VAR_PAD) << std::left << "  DynaML measurement file: " << p.i.xml_msrfile << std::endl;
			}				
		}
		if (p.i.export_dna_files)
		{
			std::cout << std::setw(PRINT_VAR_PAD) << std::left << "  DNA station file: " << p.i.dna_stnfile << std::endl;
			std::cout << std::setw(PRINT_VAR_PAD) << std::left << "  DNA measurement file: " << p.i.dna_msrfile << std::endl;
		}

		if (p.i.simulate_measurements)
		{
			std::cout << std::setw(PRINT_VAR_PAD) << std::left << "  DNA simulated msr file: " << p.i.simulate_msrfile << std::endl;
		}
		
		if (!p.i.bounding_box.empty())
		{
			std::cout << std::setw(PRINT_VAR_PAD) << std::left << "  Bounding box: " << p.i.bounding_box << std::endl;
			if (p.i.split_clusters)
				std::cout << std::setw(PRINT_VAR_PAD) << std::left << "  Split GNSS clusters: " << (p.i.split_clusters ? "Yes" : "No") << std::endl;
		}
		else
		{
			if (!p.i.stn_associated_msr_include.empty())
				std::cout << std::setw(PRINT_VAR_PAD) << std::left << "  Stations to include: " << p.i.stn_associated_msr_include << std::endl;
			if (!p.i.stn_associated_msr_exclude.empty())
				std::cout << std::setw(PRINT_VAR_PAD) << std::left << "  Stations to exclude: " << p.i.stn_associated_msr_exclude << std::endl;
			
			if (p.i.split_clusters)
				std::cout << std::setw(PRINT_VAR_PAD) << std::left << "  Split GNSS clusters: " << (p.i.split_clusters ? "Yes" : "No") << std::endl;
		}
		
		if (p.i.import_block)
		{
			if (!p.i.seg_file.empty())
				std::cout << std::setw(PRINT_VAR_PAD) << std::left << "  Segmentation file: " << p.i.seg_file << std::endl;
			std::cout << std::setw(PRINT_VAR_PAD) << std::left << "  Import stns & msrs from block: " << p.i.import_block_number << std::endl;
		}
		else if (p.i.import_network)
		{
			if (!p.i.seg_file.empty())
				std::cout << std::setw(PRINT_VAR_PAD) << std::left << "  Segmentation file: " << p.i.seg_file << std::endl;
			std::cout << std::setw(PRINT_VAR_PAD) << std::left << "  Import stns & msrs from network: " << p.i.import_network_number << std::endl;
		}

		if (!p.i.scalar_file.empty())
			std::cout << std::setw(PRINT_VAR_PAD) << std::left << "  GNSS baseline scalar file: " << p.i.scalar_file << std::endl;
		
		std::cout << std::endl;
	}

	PrintOutputFileHeaderInfo(&imp_file, p.i.imp_file, &p, "DYNADJUST IMPORT LOG FILE");

	MsrTally parsemsrTally;
	StnTally parsestnTally;

	CDnaProjection projection(UTM);

	vifm_t vinput_file_meta;

	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	// start "total" time
	dynadjust::cpu_timer time;
	
	
	// Import discontinuity file and apply discontinuities
	// Due to the structure and format of SINEX files, it is essential that
	// discontinuities be parsed prior to reading any SINEX files.
	if (vm.count(STATION_DISCONTINUITY_FILE))
	{
		p.i.apply_discontinuities = true;
		
		// Does it exist?
		if (!std::filesystem::exists(p.i.stn_discontinuityfile))
		{
			std::filesystem::path discontPath(p.i.stn_discontinuityfile);
			std::stringstream ss;
			ss << "- Warning: The station discontinuity file " << discontPath.filename().string() << " does not exist... ignoring discontinuity input." << std::endl;
			imp_file << std::endl << ss.str();
		}
		else
		{
			if (!p.g.quiet)
			{
				std::cout << "+ Importing station discontinuities from " << p.i.stn_discontinuityfile << "... ";
				std::cout.flush();
			}
			imp_file << "+ Importing station discontinuities from " << p.i.stn_discontinuityfile << "... ";

			parserDynaML.ParseDiscontinuities(p.i.stn_discontinuityfile);

			if (!p.g.quiet)
				std::cout << "Done." << std::endl;
			imp_file << "Done." << std::endl;
		}

		if (p.i.export_discont_file)
		{
			if (!p.g.quiet)
			{
				std::cout << "+ Exporting discontinuity information to text file... ";
				std::cout.flush();
			}
			imp_file << "+ Exporting discontinuity information to text file... ";
			parserDynaML.SerialiseDiscontTextFile(p.i.stn_discontinuityfile);
			if (!p.g.quiet)
			{
				std::cout << "Done." << std::endl;
				std::cout.flush();
			}
			imp_file << "Done." << std::endl;
		}
	}

	// Now, set the 'default' epoch in the binary station and measurement files
	std::string default_datum = p.i.reference_frame;
	
	// Import network information based on a segmentation block?
	if (p.i.import_block)
	{
		imp_file << "+ Extracting stations and measurements from segmented block " << 
			p.i.import_block_number << "... ";
		
		if (ImportSegmentedBlock(parserDynaML, &vstationsTotal, &vmeasurementsTotal, 
			&parsestnTally, &parsemsrTally, p) != EXIT_SUCCESS)
			return EXIT_FAILURE;

		imp_file << "Done." << std::endl;
	}
	// Import network information based on a contiguous network?
	else if (p.i.import_network)
	{
		imp_file << "+ Extracting stations and measurements from contiguous network " <<
			p.i.import_network_number << "... ";
		
		if (ImportContiguousNetwork(parserDynaML, &vstationsTotal, &vmeasurementsTotal,
			&parsestnTally, &parsemsrTally, p) != EXIT_SUCCESS)
			return EXIT_FAILURE;

		imp_file << "Done." << std::endl;
	}
	// Import data as normal
	else
	{
		// Change current directory to the import folder
		// A hack to circumvent the problem caused by importing DynaML files in 
		// different directories to where import is run from, causing errors 
		// because DynaML.xsd cannot be found.
		std::filesystem::path currentPath(std::filesystem::current_path());
		std::filesystem::current_path(std::filesystem::path(p.g.input_folder));

		// Import all data as-is.
		// All filtering is performed later below
		if (ImportDataFiles(parserDynaML, &vStations, &vMeasurements, &vstationsTotal, &vmeasurementsTotal,
			&imp_file, &vinput_file_meta, &parsestnTally, &parsemsrTally, errorCount, p) != EXIT_SUCCESS)
			return EXIT_FAILURE;

		current_path(currentPath);
	}

	epsgCode = epsgCodeFromName<UINT32>(p.i.reference_frame);

	if (!p.g.quiet)
		std::cout << std::endl;
	imp_file << std::endl;

	vstring vPoorlyConstrainedStns;

	// Ignore measurements that do not sufficiently (in themselves) allow for a station to be estimated
	if (p.i.ignore_insufficient_msrs == 1)
	{
		size_t msrCount = vmeasurementsTotal.size();
		if (!p.g.quiet)
		{
			std::cout << "+ Identifying stations with insufficient measurements";
			std::cout.flush();
		}
		imp_file << "+ Identifying stations with insufficient measurements...";
		parserDynaML.IgnoreInsufficientMeasurements(&vstationsTotal, &vmeasurementsTotal, &vPoorlyConstrainedStns);
		if (!p.g.quiet)
		{
			std::cout << " Done. " << std::endl;
			if (msrCount > vmeasurementsTotal.size())
				std::cout << "+ Removed " << (msrCount - vmeasurementsTotal.size()) << " measurements which alone do not sufficiently allow for" << std::endl <<
					"  the estimation of 2D coordinates." << std::endl;

		}
		imp_file << " Done. " << std::endl;
		if (msrCount > vmeasurementsTotal.size())
			imp_file << "+ Removed " << (msrCount - vmeasurementsTotal.size()) << " measurements which alone do not sufficiently allow for" << std::endl <<
				"  the estimation of 2D coordinates." << std::endl;

		if (!vPoorlyConstrainedStns.empty())
		{
			if (vPoorlyConstrainedStns.size() == 1)
			{
				if (p.g.verbose > 2)
				{
					if (!p.g.quiet)
						std::cout << "- Warning: station " << vPoorlyConstrainedStns.at(0) << " is not associated with sufficient measurements." << std::endl;
					imp_file << "- Warning: station " << vPoorlyConstrainedStns.at(0) << " is not associated with sufficient measurements." << std::endl;
				}
				else
				{
					if (!p.g.quiet)
						std::cout << "- Warning: 1 station is not associated with sufficient measurements." << std::endl;
					imp_file << "- Warning: 1 station is not associated with sufficient measurements." << std::endl;
				}
			}
			else
			{
				if (p.g.verbose > 2)
				{
					if (!p.g.quiet)
						std::cout << "- Warning: The following " << vPoorlyConstrainedStns.size() << " stations are not associated with sufficient measurements." << std::endl;
					imp_file << "- Warning: The following " << vPoorlyConstrainedStns.size() << " stations are not associated with sufficient measurements." << std::endl;
					_it_vstr poorly;
					for (poorly = vPoorlyConstrainedStns.begin(); poorly != vPoorlyConstrainedStns.end(); poorly++)
					{
						if (!p.g.quiet)
							outputObject(std::string("  - " + *poorly + "\n"), std::cout);
						outputObject(std::string("  - " + *poorly + "\n"), imp_file);
					}
				}
				else
				{
					if (!p.g.quiet)
						std::cout << "- Warning: " << vPoorlyConstrainedStns.size() << " stations are not associated with sufficient measurements." << std::endl;
					imp_file << "- Warning: " << vPoorlyConstrainedStns.size() << " stations are not associated with sufficient measurements." << std::endl;
				}
			}
		}
	}

	// Remove ignored measurements (if supplied)
	if (p.i.remove_ignored_msr)
	{
		if (!p.g.quiet)
		{
			std::cout << "+ Removing ignored measurements... ";
			std::cout.flush();
		}
		imp_file << "+ Removing ignored measurements... ";
		parserDynaML.RemoveIgnoredMeasurements(&vmeasurementsTotal, &parsemsrTally);
		if (!p.g.quiet)
			std::cout << "Done. " << std::endl;
		imp_file << "Done. " << std::endl;
	}

	// Strip all measurements except required measurements (if supplied)
	if (!p.i.include_msrs.empty())
	{
		if (!p.g.quiet)
		{
			std::cout << "+ Stripping all measurements except types " << p.i.include_msrs << "... ";
			std::cout.flush();
		}
		imp_file << "+ Stripping all measurements except types " << p.i.include_msrs << "... ";
		parserDynaML.IncludeMeasurementTypes(p.i.include_msrs, &vmeasurementsTotal, &parsemsrTally);
		if (!p.g.quiet)
			std::cout << "Done. " << std::endl;
		imp_file << "Done. " << std::endl;
	}
	
	// Strip all unwanted measurements (if supplied)
	if (!p.i.exclude_msrs.empty())
	{
		if (!p.g.quiet)
		{
			std::cout << "+ Stripping measurement types " << p.i.exclude_msrs << "... ";
			std::cout.flush();
		}
		imp_file << "+ Stripping measurement types " << p.i.exclude_msrs << "... ";
		parserDynaML.ExcludeMeasurementTypes(p.i.exclude_msrs, &vmeasurementsTotal, &parsemsrTally);
		if (!p.g.quiet)
			std::cout << "Done. " << std::endl;
		imp_file << "Done. " << std::endl;
	}

	// Reduce stations.
	// But, only reduce stations if not importing for a segmentation block or network
	// The reason stations are not reduced if IMPORT_SEG_BLOCK or IMPORT_CONTIG_NET has been set is
	// because stations in the binary file will have already been reduced
	// on last import!
	if (p.i.import_block == 0 && p.i.import_network == 0 && vstationsTotal.size())
	{
		// reduce stations (e.g. convert from UTM to LLH)
		try {
			if (!p.g.quiet)
			{
				std::cout << "+ Reducing stations... ";
				std::cout.flush();
			}
			imp_file << "+ Reducing stations... ";
			parserDynaML.ReduceStations(&vstationsTotal, projection);
			if (!p.g.quiet)
				std::cout << "Done." << std::endl;
			imp_file << "Done." << std::endl;
		} 
		catch (const XMLInteropException& e) {
			std::cout << std::endl << "- Error: " << e.what() << std::endl;
			imp_file << std::endl << "- Error: " << e.what() << std::endl;
			imp_file.close();
			return EXIT_FAILURE;
		}
	}
	
	/////////////////////////////////////////////////////////////////////////
	// Add discontinuity sites to vStations
	//
	// WARNING: If discontinuity sites are present in the renaming file, and
	// the renaming file changes the name from the four-character id in the 
	// discontinuity file to another name, this function call will not work
	// properly. Hence, sites in the discontinuity file must be consistently
	// named in the station and measurement files and not renamed to another
	// name.
	if (p.i.apply_discontinuities && !vstationsTotal.empty())
		parserDynaML.AddDiscontinuityStations(&vstationsTotal);

	UINT32 stn;

	// Extract user-defined stations and all connected measurements
	if (vstationsTotal.size() &&	// cannot be empty
		(!p.i.stn_associated_msr_include.empty() || !p.i.stn_associated_msr_exclude.empty()))
	{
		vstring vUnusedStns;
		bool splitXmsrs(false), splitYmsrs(false);

		try {
			if (!p.g.quiet)
			{
				std::cout << "+ Extracting stations and associated measurements... ";
				std::cout.flush();
			}
			imp_file << "+ Extracting stations and associated measurements... ";

			// this method reforms asl and aml
			parserDynaML.ExtractStnsAndAssociatedMsrs(p.i.stn_associated_msr_include, p.i.stn_associated_msr_exclude, &vstationsTotal, &vmeasurementsTotal, 
				&parsestnTally, &parsemsrTally, &vUnusedStns, p, splitXmsrs, splitYmsrs);
		} 
		catch (const XMLInteropException& e) {
			std::cout << "- Error: " << e.what() << std::endl;
			imp_file << std::endl << "- Error: " << e.what() << std::endl;
			imp_file.close();
			return EXIT_FAILURE;
		}

		if (!p.g.quiet)
			std::cout << "Done. " << std::endl;
		imp_file << "Done." << std::endl;

		//msrCount = vmeasurementsTotal.size();
	}	

	// Now strip all stations outside the bounding box (if supplied)
	if (vstationsTotal.size() &&	// cannot be empty
		!p.i.bounding_box.empty())
	{
		vstring vUnusedStns;
		bool splitXmsrs(false), splitYmsrs(false);

		try {
			if (!p.g.quiet)
			{
				std::cout << "+ Stripping stations and measurements outside the bounding box... ";
				std::cout.flush();
			}
			imp_file << "+ Stripping stations and measurements outside the bounding box... ";
			parserDynaML.ExcludeAllOutsideBoundingBox(&vstationsTotal, &vmeasurementsTotal, 
				&parsestnTally, &parsemsrTally, &vUnusedStns, p, splitXmsrs, splitYmsrs);
		} 
		catch (const XMLInteropException& e) {
			std::cout << "- Error: " << e.what() << std::endl;
			imp_file << std::endl << "- Error: " << e.what() << std::endl;
			imp_file.close();
			return EXIT_FAILURE;
		}

		if (!p.g.quiet)
			std::cout << "Done. " << std::endl;
		imp_file << "Done." << std::endl;

		if (!vUnusedStns.empty())
		{
			std::ostringstream ss_msg;
			ss_msg << "  - " << vstationsTotal.size() << (vstationsTotal.size() == 1 ? " station is" : " stations are") << 
				" within the bounding box." << std::endl;
			imp_file << ss_msg.str();
			if (!p.g.quiet)
				std::cout << ss_msg.str();
			ss_msg.str("");
			ss_msg << "  - ";
			if (vUnusedStns.size() > 1)
			{
				ss_msg << vUnusedStns.size() << " stations were";
			}
			else
				ss_msg << "Station " << vUnusedStns.at(0) << " was";

			ss_msg << " found outside the bounding box and " << (vUnusedStns.size() > 1 ? "have been" : "has been") << std::endl <<
				"    removed together with the corresponding measurements." << std::endl;
			if (p.i.split_clusters && (splitXmsrs || splitYmsrs))
			{
				ss_msg << "  - Note: GPS ";
				if (splitXmsrs && splitYmsrs)
					ss_msg << "point and baseline";
				else if (splitXmsrs)
					ss_msg << "baseline";
				else //if (splitYmsrs)
					ss_msg << "point";
				ss_msg << " cluster measurements straddling the limits of" << std::endl << 
				"    the bounding box have been split." << std::endl;
			}
			imp_file << ss_msg.str() << "  - Excluded stations:" << std::endl;
			if (!p.g.quiet)
				std::cout << ss_msg.str() << std::endl;
			
			for (stn=0; stn<vUnusedStns.size(); ++stn)
				imp_file << "     " << vUnusedStns.at(stn) << std::endl;
			imp_file << std::endl;
		}
		else
		{
			imp_file << "  - No stations were found outside the bounding box." << std::endl;
			if (!p.g.quiet)
				std::cout << "  - No stations were found outside the bounding box." << std::endl;
		}
	}

	UINT32 msrRead(parsemsrTally.TotalCount());
	UINT32 stnCount(static_cast<UINT32>(vstationsTotal.size()));
	UINT32 msrCount(static_cast<UINT32>(vmeasurementsTotal.size()));

	if (!p.g.quiet)
		std::cout << std::endl;
	imp_file << std::endl;

	///////////////////////////////////////////////////////////////////////
	// Ok, now that unwanted stations and measurements have been stripped,
	// provide station and measurement Summary
	//
	if ((stnCount + msrCount) > 0)
	{
		if (p.i.import_block || p.i.import_network)
		{
			if (!p.g.quiet)
				std::cout << "+ Binary file ";
			imp_file << "+ Binary file ";
		}
		else
		{
			if (!p.g.quiet)
				std::cout << "+ File ";
			imp_file << "+ File ";
		}
		
		if (!p.g.quiet)
			std::cout << "parsing summary:" << std::endl << std::endl;
		imp_file << "parsing summary:" << std::endl << std::endl;
	}
	
	//
	// Station summary
	if (stnCount)
	{
		if (!p.g.quiet)
		{
			parsestnTally.coutSummary(std::cout, std::string("  Read"));
			std::cout << std::endl;
		}
		parsestnTally.coutSummary(imp_file, std::string("  Read"));
		imp_file << std::endl;
	}
	
	//
	// Measurement summary
	if (msrCount)
	{
		if (!p.g.quiet)
		{
			parsemsrTally.coutSummary(std::cout, std::string("  Read"));
			std::cout << std::endl;
		}
		parsemsrTally.coutSummary(imp_file, std::string("  Read"));
		imp_file << std::endl;
	}

	////////////////////////////////////////////////////////////////////////
	// Can we proceed?
	if (stnCount < 1)
	{
		imp_file << "- No further processing can be done as no stations were loaded." << std::endl;
		if (!p.g.quiet)
			std::cout << "- No further processing can be done as no stations were loaded." << std::endl;
		imp_file.close();
		return PARSE_SUCCESS;
	}

	/////////////////////////////////////////////////////////////////////////
	// Rename stations
	if (p.i.rename_stations && (stnCount > 0 || msrCount > 0))
	{
		// Does it exist?
		if (!std::filesystem::exists(p.i.stn_renamingfile))
			// Look for it in the input folder
			p.i.stn_renamingfile = formPath<std::string>(p.g.input_folder, leafStr<std::string>(p.i.stn_renamingfile));

		// Apply renaming
		try {
			if (!p.g.quiet)
			{
				std::cout << "+ Renaming stations... ";
				std::cout.flush();
			}
			imp_file << "+ Renaming stations... ";
			parserDynaML.RenameStations(&vstationsTotal, &vmeasurementsTotal, &p);
			if (!p.g.quiet)
				std::cout << "Done." << std::endl;
			imp_file << "Done." << std::endl;
		} 
		catch (const XMLInteropException& e) {
			std::cout << std::endl << "- Error: " << e.what() << std::endl;
			imp_file << std::endl << "- Error: " << e.what() << std::endl;
			imp_file.close();
			return EXIT_FAILURE;
		}
	}

	// Apply GNSS scaling (if required)
	if (p.i.import_block == 0 && p.i.import_network == 0 && p.i.apply_scaling)
	{
		// Apply scaling
		try {
			if (!p.g.quiet)
			{
				std::cout << "+ Applying scalars to GNSS measurements... ";
				std::cout.flush();
			}
			imp_file << "+ Applying scalars to GNSS measurements... ";
			parserDynaML.EditGNSSMsrScalars(&vmeasurementsTotal, &p);
			if (!p.g.quiet)
				std::cout << "Done." << std::endl;
			imp_file << "Done." << std::endl;
		} 
		catch (const XMLInteropException& e) {
			std::cout << std::endl << "- Error: " << e.what() << std::endl;
			imp_file << std::endl << "- Error: " << e.what() << std::endl;
			imp_file.close();
			return EXIT_FAILURE;
		}
	}

	std::ofstream dst_file;

	// flush duplicate stations/measurements files
	try {
		if (std::filesystem::exists(p.i.dst_file))
			std::filesystem::remove(p.i.dst_file);
		if (std::filesystem::exists(p.i.dms_file))
			std::filesystem::remove(p.i.dms_file);
	}
	catch (const std::ios_base::failure& f) { 
		// do nothing on failure.
		imp_file << std::endl << "- Warning: " << f.what() << std::endl;
	}

	// Prepare file names if importing from a segmentation block or contiguous network
	if (p.i.import_block == 1 || p.i.import_network == 1)
	{
		std::stringstream modifier("");
		if (p.i.import_block == 1)
			modifier << ".block-" << p.i.import_block_number;
		else if (p.i.import_network == 1)
			modifier << ".network-" << p.i.import_network_number;

		// create new output file names based on block number
		// reform file name for each so as to preserve full path for each file
		std::stringstream ss("");
		ss << formPath<std::string>(std::filesystem::path(p.i.bst_file).parent_path().generic_string(), std::filesystem::path(p.i.bst_file).stem().generic_string());
		ss << modifier.str() << ".bst";
		p.i.bst_file = ss.str();

		ss.str("");
		ss << formPath<std::string>(std::filesystem::path(p.i.bms_file).parent_path().generic_string(), std::filesystem::path(p.i.bms_file).stem().generic_string());
		ss << modifier.str() << ".bms";
		p.i.bms_file = ss.str();

		ss.str("");
		ss << formPath<std::string>(std::filesystem::path(p.i.asl_file).parent_path().generic_string(), std::filesystem::path(p.i.asl_file).stem().generic_string());
		ss << modifier.str() << ".asl";
		p.i.asl_file = ss.str();

		ss.str("");
		ss << formPath<std::string>(std::filesystem::path(p.i.aml_file).parent_path().generic_string(), std::filesystem::path(p.i.aml_file).stem().generic_string());
		ss << modifier.str() << ".aml";
		p.i.aml_file = ss.str();

		ss.str("");
		ss << formPath<std::string>(std::filesystem::path(p.i.map_file).parent_path().generic_string(), std::filesystem::path(p.i.map_file).stem().generic_string());
		ss << modifier.str() << ".map";
		p.i.map_file = ss.str();

		ss.str("");
		ss << formPath<std::string>(std::filesystem::path(p.i.map_file).parent_path().generic_string(), std::filesystem::path(p.i.map_file).stem().generic_string());
		ss << ".dbid";
		p.i.dbid_file = ss.str();
	}

	// Remove duplicate stations. If required, test nearby stations
	if (stnCount > 0) {
		
		vstring vduplicateStns;
		v_stringstring_doubledouble_pair vnearbyStns;

		// Remove duplicates and, if required, identify station pairs 
		// separated by distances less than search_stn_radius
		try {
			std::ostringstream ss_msg;
			if (p.i.search_nearby_stn)
				ss_msg << "+ Testing for duplicate and nearby stations... ";
			else
				ss_msg << "+ Testing for duplicate stations... ";
			
			if (!p.g.quiet)
			{
				std::cout << ss_msg.str();
				std::cout.flush();
			}

			imp_file << ss_msg.str();			
			
			stn = parserDynaML.RemoveDuplicateStations(&vstationsTotal, &vduplicateStns, &vnearbyStns);
			
			if (!p.g.quiet)
				std::cout << "Done. ";
			imp_file << "Done. ";

			if (stn > 0)
			{
				try {
					// Create duplicate station file
					file_opener(dst_file, p.i.dst_file);
				}
				catch (const std::runtime_error& e) {
					ss_msg << "- Error: Could not open " << p.i.dst_file << ". \n  Check that the file exists and that the file is not already opened." << 
					 	std::endl << e.what() << std::endl;
					if (!p.g.quiet)
						std::cout << ss_msg.str();
					imp_file << ss_msg.str();
					imp_file.close();
					return EXIT_FAILURE;
				}
						
				PrintOutputFileHeaderInfo(&dst_file, p.i.dst_file, &p, "DUPLICATE STATION FILE");

				if (!vduplicateStns.empty())
				{
					ss_msg.str("");
					ss_msg << "Removed " << vduplicateStns.size() << " duplicate station" << (vduplicateStns.size() > 1 ? "s" : "");
					
					// print message to .dst file
					dst_file << ss_msg.str() << ":" << std::endl;
					for (stn=0; stn<vduplicateStns.size(); ++stn)
						dst_file << "  - " << vduplicateStns.at(stn) << std::endl;
					dst_file << std::endl;

					ss_msg << "." << std::endl << "  See " << p.i.dst_file << " for details." << std::endl;

					imp_file << ss_msg.str();
					if (!p.g.quiet)
						std::cout << std::endl << "- Warning: " << ss_msg.str();					
				}				

				if (!vnearbyStns.empty())
				{
					ss_msg.str("");
					ss_msg << vnearbyStns.size() << (vnearbyStns.size() > 1 ? " pairs of stations were" : " pair of station was") << 
						" found to be separated by less than " << std::setprecision(3) << p.i.search_stn_radius << "m.";

					imp_file << std::endl << "- Warning: " << ss_msg.str() << std::endl << 
						"  See " << p.i.dst_file << " for details." << std::endl << std::endl;
					
					if (!p.g.quiet)
						std::cout << std::endl << "- Warning: " << ss_msg.str() << std::endl << 
						"  See " << p.i.dst_file << " for details." << std::endl << std::endl;
					ss_msg.str("");

					// output message
					dst_file << "Nearby station search results:" << std::endl << ss_msg.str() << std::endl << std::endl;

					dst_file <<  
						std::setw(HEADER_20) << std::left << "First station" << 
						std::setw(HEADER_20) << "Nearby station" << 
						std::setw(HEADER_20) << std::right << "Separation (m)" << 
						std::setw(HEADER_20) << "Diff height (m)" << 
					 	std::endl;

					for (UINT32 i(0); i<(HEADER_20*4); ++i)
						dst_file << "-";
					dst_file << std::endl;

					// dump nearby stations to dst file
					for (stn=0; stn<vnearbyStns.size(); ++stn)
					{
						dst_file <<
							std::setw(HEADER_20) << std::left << vnearbyStns.at(stn).first.first << 									// First
							std::setw(HEADER_20) << vnearbyStns.at(stn).first.second <<									// Nearby
							std::setw(HEADER_20) << std::setprecision(3) << std::fixed << std::right << vnearbyStns.at(stn).second.first <<		// Separation (m)
							std::setw(HEADER_20) << std::setprecision(3) << std::fixed << vnearbyStns.at(stn).second.second <<		// Diff height (m)
						 	std::endl;
					}

					ss_msg <<
						"  If the names in each pair refer to the same station, then update the " << std::endl <<
						"  station and measurement files with the correct station name and re-run " << __BINARY_NAME__ << "." << std::endl <<
						"  Alternatively, if the names in each pair are unique, either call " << __BINARY_NAME__ << std::endl <<
						"  without the --" << TEST_NEARBY_STNS << " option, or decrease the radial search " << std::endl <<
						"  distance using --" << TEST_NEARBY_STN_DIST << "." << std::endl << std::endl;

					if (!p.g.quiet)
						std::cout << ss_msg.str();
					imp_file << ss_msg.str();
					imp_file.close();
					dst_file.close();
					return EXIT_FAILURE;
				}
				else if (p.g.verbose == 3)
				{
					if (!p.g.quiet)
						std::cout << "+ Total number of unique stations is " << vstationsTotal.size() << std::endl;
					imp_file << "+ Total number of unique stations is " << vstationsTotal.size() << std::endl;
				}

				// If this line is reached, then there are no nearby stations
				dst_file.close();
			}
			if (!p.g.quiet)
				std::cout << std::endl;
			imp_file << std::endl;
		} 
		catch (const XMLInteropException& e) {
			std::cout << std::endl << "- Error: " << e.what() << std::endl;
			imp_file << std::endl << "- Error: " << e.what() << std::endl;
			imp_file.close();
			return EXIT_FAILURE;
		}		

		/////////////////////////////////////////////////////////////////
		// Now commence sorting and mapping
		// 1. Sort stations
		// 2. Create station map
		try {
			if (!p.g.quiet)
			{
				std::cout << "+ Sorting stations... ";
				std::cout.flush();
			}
			imp_file << "+ Sorting stations... ";
			parserDynaML.FullSortandMapStations((vdnaStnPtr*) &vstationsTotal, &vStnsMap_sortName);
			if (!p.g.quiet)
			{
				std::cout << "Done." << std::endl << "+ Serialising station map... ";
				std::cout.flush();
			}
			imp_file << "Done." << std::endl << "+ Serialising station map... ";			
			parserDynaML.SerialiseMap(p.i.map_file);		// parserDynaML keeps the map in memory
			stn_map_created = true;
			if (!p.g.quiet)
			{
				std::cout << "Done." << std::endl;
				std::cout.flush();
			}
			imp_file << "Done." << std::endl;

			if (p.i.export_map_file)
			{
				if (!p.g.quiet)
				{
					std::cout << "+ Exporting station map to text file... ";
					std::cout.flush();
				}
				imp_file << "+ Exporting station map to text file... ";
				parserDynaML.SerialiseMapTextFile(p.i.map_file);
				if (!p.g.quiet)
				{
					std::cout << "Done." << std::endl;
					std::cout.flush();
				}
				imp_file << "Done." << std::endl;				
			}
		} 
		catch (const XMLInteropException& e) {
			std::cout << std::endl << "- Error: " << e.what() << std::endl;
			imp_file << std::endl << "- Error: " << e.what() << std::endl;
			imp_file.close();
			return EXIT_FAILURE;
		}	
	}

	// Search for similar measurements
	if (msrCount > 0 && (p.i.search_similar_msr || p.i.search_similar_msr_gx || p.i.ignore_similar_msr)) 
	{
		if (SearchForSimilarMeasurements(&parserDynaML, &p, &imp_file,
			&vmeasurementsTotal) != EXIT_SUCCESS)
			return EXIT_FAILURE;	
	}

	// Import DNA geoid file
	if (p.i.import_geo_file && p.i.import_block == 0 && p.i.import_network == 0 && vstationsTotal.size())
	{
		if (!std::filesystem::exists(p.i.geo_file))
		{
			std::filesystem::path geoPath(p.i.geo_file);
			std::stringstream ss;
			ss << "- Error: The geoid file " << geoPath.filename().string() << " does not exist." << std::endl;
			std::cout << std::endl << ss.str();
			imp_file << std::endl << ss.str();
			return EXIT_FAILURE;
		}

		if (!p.g.quiet)
		{
			std::cout << "+ Importing geoid information from " << p.i.geo_file << "... ";
			std::cout.flush();
		}
		imp_file << "+ Importing geoid information from " << p.i.geo_file << "... ";

		parserDynaML.LoadDNAGeoidFile(p.i.geo_file, &vstationsTotal);

		if (!p.g.quiet)
			std::cout << "Done." << std::endl;
		imp_file << "Done." << std::endl;
	}

	UINT32 ignMsrCount, mapCount;
	vstring vunusedStations;
	vUINT32 vignoredMeasurements;

	// Map measurements to stations
	if (stn_map_created) 
	{	
		try {
			if (!p.g.quiet)
			{
				std::cout << "+ Mapping measurements to stations... ";
				std::cout.flush();
			}
			imp_file << "+ Mapping measurements to stations... ";
			
			parserDynaML.MapMeasurementStations((vdnaMsrPtr*) &vmeasurementsTotal,	&associatedSL, &mapCount, 
				&vunusedStations, &vignoredMeasurements);
			
			if (!p.g.quiet)
				std::cout << "Done." << std::endl;
			imp_file << "Done. ";

			if (msrCount > 0)
			{
				measurements_mapped = true;
				imp_file << "Mapped " << mapCount << " measurements to " << vStnsMap_sortName.size() - vunusedStations.size() << " stations." << std::endl;
			}
			else
				imp_file << std::endl;

			if (!vunusedStations.empty())
			{
				if (vunusedStations.size() == 1)
				{
					if (p.g.verbose > 2)
					{
						if (!p.g.quiet)
							std::cout << "- Warning: station " << vunusedStations.at(0) << " was not associated with any measurements." << std::endl;
						imp_file << "- Warning: station " << vunusedStations.at(0) << " was not associated with any measurements." << std::endl;
					}
					else
					{
						if (!p.g.quiet)
							std::cout << "- Warning: " << vunusedStations.size() << " station was not associated with any measurements." << std::endl;
						imp_file << "- Warning: " << vunusedStations.size() << " station was not associated with any measurements." << std::endl;
					}
				}
				else
				{
					if (p.g.verbose > 2)
					{
						if (!p.g.quiet)
							std::cout << "- Warning: The following " << vunusedStations.size() << " stations were not associated with any measurements." << std::endl;
						imp_file << "- Warning: The following " << vunusedStations.size() << " stations were not associated with any measurements." << std::endl;
						_it_vstr unused;
						for (unused = vunusedStations.begin(); unused!=vunusedStations.end(); unused++)
						{
							if (!p.g.quiet)
								outputObject(std::string("  - " + *unused + "\n"), std::cout);
							outputObject(std::string("  - " + *unused + "\n"), imp_file);
						}
					}
					else
					{
						if (!p.g.quiet)
							std::cout << "- Warning: " << vunusedStations.size() << " stations were not associated with any measurements." << std::endl;
						imp_file << "- Warning: " << vunusedStations.size() << " stations were not associated with any measurements." << std::endl;
					}
				}
			}

			if (msrRead < mapCount && vignoredMeasurements.empty())
			{
				if (!p.g.quiet)
					std::cout << "- Warning: Not all measurements were mapped: " << msrRead << " msrs read vs. " << mapCount
						<< " msrs mapped." << std::endl;
			}
			else if (msrRead != mapCount && !vignoredMeasurements.empty())
			{
				ignMsrCount = parserDynaML.ComputeMeasurementCount(&vmeasurementsTotal, vignoredMeasurements);	
				if (!p.g.quiet)
				{
					std::cout << "- Warning: " << ignMsrCount << " ignored measurements were not mapped." << std::endl;
					if ((msrRead - mapCount) != ignMsrCount)
						std::cout << "-          " << msrRead << " m.read vs. " << mapCount << " m.mapped." << std::endl;
				}
			}
			if (!p.g.quiet)
				std::cout.flush();
		} 
		catch (const XMLInteropException& e) {
			std::cout << std::endl << "- Error: " << e.what() << std::endl;
			imp_file << std::endl << "- Error: " << e.what() << std::endl;
			imp_file.close();
			return EXIT_FAILURE;
		}
	}

	stnCount = static_cast<UINT32>(vstationsTotal.size());
	if (stnCount > static_cast<UINT32>(associatedSL.size()))
		stnCount = static_cast<UINT32>(associatedSL.size());
	
	// Does the user want to simulate measurements?
	if (p.i.simulate_measurements && stnCount > 0 && msrCount > 0) 
	{
		// No need for the following code if the outer 'if' statement checks 
		// for non-zero stnCount value
		//
		// if (stnCount == 0)
		// {
		// 	if (!p.g.quiet)
		// 		std::cout << "- Error: there are no stations from which to simulate measurements." << std::endl;
		// 	imp_file << "- Error: there are no stations from which to simulate measurements." << std::endl;
		// 	imp_file.close();
		// 	return EXIT_FAILURE;
		// }	
	
		try {
			// Simulate measurements
			if (!p.g.quiet)
			{
				std::cout << "+ Simulating and exporting measurements to " << leafStr<std::string>(p.i.simulate_msrfile) << "... ";
				std::cout.flush();
			}
			imp_file << "+ Simulating and exporting measurements to " << leafStr<std::string>(p.i.simulate_msrfile) << "... ";
			parserDynaML.SimulateMSR(
				((vdnaStnPtr*) &vstationsTotal), 
				((vdnaMsrPtr*) &vmeasurementsTotal), 
				p.i.simulate_msrfile, p);
			if (!p.g.quiet)
				std::cout << "Done." << std::endl;
			imp_file << "Done." << std::endl;
		}
		catch (const XMLInteropException& e) {
			std::cout.flush();
			std::cout << std::endl << "- Error: " << e.what() << std::endl;
			imp_file << std::endl << "- Error: " << e.what() << std::endl;
			imp_file.close();
			return EXIT_FAILURE;
		}
	}

	// Create association lists
	if (measurements_mapped) 
	{
		try {
			if (!p.g.quiet)
			{
				std::cout << "+ Creating association lists... ";
				std::cout.flush();
			}
			imp_file << "+ Creating association lists... ";
			parserDynaML.CompleteAssociationLists(&vmeasurementsTotal, &associatedSL, &associatedML);
			if (!p.g.quiet)
			{
				std::cout << "Done." << std::endl;
				std::cout.flush();
			}
			imp_file << "Done." << std::endl;
		} 
		catch (const XMLInteropException& e) {
			std::cout << "- Error: " << e.what() << std::endl;
			imp_file << std::endl << "- Error: " << e.what() << std::endl;
			imp_file.close();
			return EXIT_FAILURE;
		}

		try {
			if (!p.g.quiet)
			{
				std::cout << "+ Serialising association lists... ";
				std::cout.flush();
			}
			imp_file << "+ Serialising association lists... ";
			parserDynaML.SerialiseAsl(p.i.asl_file, &associatedSL);
			parserDynaML.SerialiseAml(p.i.aml_file, &associatedML);
			if (!p.g.quiet)
			{
				std::cout << "Done." << std::endl;
				std::cout.flush();
			}
			imp_file << "Done." << std::endl;
		} 
		catch (const XMLInteropException& e) {
			std::cout << "- Error: " << e.what() << std::endl;
			imp_file << std::endl << "- Error: " << e.what() << std::endl;
			imp_file.close();
			return EXIT_FAILURE;
		}
	}

	// Create binary station file.
	// Binary station file can be serialised without measurements
	if (stnCount > 0)
	{
		try {
			if (!vstationsTotal.empty())
			{
				if (!p.g.quiet)
				{
					std::cout << "+ Serialising binary station file " << leafStr<std::string>(p.i.bst_file) << "... ";
					std::cout.flush();
				}
				imp_file << "+ Serialising binary station file " << leafStr<std::string>(p.i.bst_file) << "... ";
				
				parserDynaML.SerialiseBst(
					p.i.bst_file, ((vdnaStnPtr*) &vstationsTotal), &vunusedStations,
					vinput_file_meta, (p.i.flag_unused_stn ? true : false));
				if (!p.g.quiet)
				{
					std::cout << "Done." << std::endl;
					std::cout.flush();
				}
				imp_file << "Done." << std::endl;
			}
		} 
		catch (const XMLInteropException& e) {
			std::cout << "- Error: " << e.what() << std::endl;
			imp_file << std::endl << "- Error: " << e.what() << std::endl;
			imp_file.close();
			return EXIT_FAILURE;
		}
	}

	// Create binary measurement file.
	// Binary measurement file can be serialised without stations
	if (msrCount > 0) 
	{
		try {
			if (!vmeasurementsTotal.empty())
			{
				if (!p.g.quiet)
				{
					std::cout << "+ Serialising binary measurement file " << leafStr<std::string>(p.i.bms_file) << "... ";
					std::cout.flush();
				}
				imp_file << "+ Serialising binary measurement file " << leafStr<std::string>(p.i.bms_file) << "... ";

				parserDynaML.SerialiseBms(
					p.i.bms_file, ((vdnaMsrPtr*) &vmeasurementsTotal),
					vinput_file_meta);
				if (!p.g.quiet)
				{
					std::cout << "Done." << std::endl;
					std::cout.flush();
				}
				imp_file << "Done." << std::endl;
			}
		} 
		catch (const XMLInteropException& e) {
			std::cout << "- Error: " << e.what() << std::endl;
			imp_file << std::endl << "- Error: " << e.what() << std::endl;
			imp_file.close();
			return EXIT_FAILURE;
		}
	}

	// Export ASL and AML to text if required
	if (measurements_mapped) 
	{
		try {
			if (p.i.export_asl_file)
			{
				if (!p.g.quiet)
				{
					std::cout << "+ Exporting associated station list to text file... ";
					std::cout.flush();
				}
				imp_file << "+ Exporting associated station list to text file... ";
				parserDynaML.SerialiseAslTextFile(p.i.asl_file, &associatedSL, (vdnaStnPtr*) &vstationsTotal);
				if (!p.g.quiet)
				{
					std::cout << "Done." << std::endl;
					std::cout.flush();
				}
				imp_file << "Done." << std::endl;				
			}
		} 
		catch (const XMLInteropException& e) {
			std::cout << "- Error: " << e.what() << std::endl;
			imp_file << std::endl << "- Error: " << e.what() << std::endl;
			imp_file.close();
			return EXIT_FAILURE;
		}

		try {
			if (p.i.export_aml_file)
			{
				if (!p.g.quiet)
				{
					std::cout << "+ Exporting associated measurement list to text file... ";
					std::cout.flush();
				}
				imp_file << "+ Exporting associated measurement list to text file... ";
				parserDynaML.SerialiseAmlTextFile(p.i.bms_file, p.i.aml_file, &associatedML, &associatedSL, (vdnaStnPtr*) &vstationsTotal);
				if (!p.g.quiet)
				{
					std::cout << "Done." << std::endl;
					std::cout.flush();
				}
				imp_file << "Done." << std::endl;				
			}
		}
		catch (const XMLInteropException& e) {
			std::cout << "- Error: " << e.what() << std::endl;
			imp_file << std::endl << "- Error: " << e.what() << std::endl;
			imp_file.close();
			return EXIT_FAILURE;
		}

		if (p.i.test_integrity)
		{
			if (!p.g.quiet)
			{
				std::cout << "+ Testing internal integrity of ASL, AML and binary files... ";
				std::cout.flush();
			}
			std::ifstream binaryMS(p.i.bms_file.c_str(), std::ios::in | std::ios::binary | std::ifstream::ate);	/// Open and seek to end immediately after opening.
			if (!binaryMS.good())
			{
				std::cout << std::endl << "- Could not open binary file for reading." << std::endl;
				imp_file << std::endl << "- Could not open binary file for reading." << std::endl;
				imp_file.close();
				return EXIT_FAILURE;
			}
			// get size, then go back to beginning
			size_t sFileSize = (size_t)binaryMS.tellg();
			binaryMS.seekg(0, std::ios::beg);

			measurement_t measRecord;
			
			UINT32 amlcount, amlindex, msr;
			for (stn=0; stn<stnCount; stn++)
			{
				amlcount = associatedSL.at(stn).get()->GetAssocMsrCount();
				
				// The following is to be used when reading from binary file
				for (msr=0; msr<amlcount; msr++)
				{
					amlindex = associatedSL.at(stn).get()->GetAMLStnIndex() + msr;
					if (associatedML.at(amlindex) * sizeof(measurement_t) >= sFileSize)
					{
						std::cout << "Error: index " << associatedML.at(amlindex) << " is out of range for the binary file." << std::endl;
						continue;
					}
					binaryMS.seekg(sizeof(UINT32) + associatedML.at(amlindex) * sizeof(measurement_t), std::ios::beg);
					binaryMS.read(reinterpret_cast<char *>(&measRecord), sizeof(measurement_t));
				}
			}
			binaryMS.close();
			if (!p.g.quiet)
				std::cout << "OK." << std::endl;
		}
	}

	try {
		// Print measurements to stations table
		PrintMeasurementstoStations(&parsemsrTally, &parserDynaML, &p, &associatedSL);
	}
	catch (const XMLInteropException& e) {
		std::cout.flush();
		std::cout << std::endl << "- Error: " << e.what() << std::endl;
		imp_file << std::endl << "- Error: " << e.what() << std::endl;
		imp_file.close();
		return EXIT_FAILURE;
	}


	// Serialise database ids
	if (p.i.import_block == 1 || p.i.import_network == 1)
	{
		// dbid_file filename formed earlier
		parserDynaML.SerialiseDatabaseId(p.i.dbid_file, &vmeasurementsTotal);
	}
	else
	{
		p.i.dbid_file = formPath<std::string>(p.g.output_folder, p.g.network_name, "dbid");
		parserDynaML.SerialiseDatabaseId(p.i.dbid_file, &vmeasurementsTotal);
	}

	// Export stations and measurements
	try {
		if (p.i.export_dynaml || p.i.export_dna_files)
			ExportStationsandMeasurements(&parserDynaML, p, &imp_file, &vinput_file_meta, 
				&vstationsTotal, &vmeasurementsTotal, stnCount, msrCount);
	}
	catch (const XMLInteropException& e) {
		std::cout.flush();
		std::cout << std::endl << "- Error: " << e.what() << std::endl;
		imp_file << std::endl << "- Error: " << e.what() << std::endl;
		imp_file.close();
		return EXIT_FAILURE;
	}
	
	if (!userSuppliedSegFile)
		p.i.seg_file = "";
	if (!userSuppliedBstFile)
		p.i.bst_file = "";
	if (!userSuppliedBmsFile)
		p.i.bms_file = "";

	// Look for a project file.  If it exists, open and load it.
	// Update the import settings.
	// Print the project file. If it doesn't exist, it will be created.
	CDnaProjectFile projectFile;
	if (std::filesystem::exists(p.g.project_file))
		projectFile.LoadProjectFile(p.g.project_file);
	
	projectFile.UpdateSettingsImport(p);
	projectFile.UpdateSettingsReftran(p);
	projectFile.UpdateSettingsOutput(p);
	projectFile.PrintProjectFile();

	if (msrCount == 0)
	{
		std::cout << "- Warning: there are no measurements to process." << std::endl;
		imp_file << "- Warning: there are no measurements to process." << std::endl;
	}
	
	if (errorCount)
	{
		if (!p.g.quiet)
			std::cout << "- Warning: some files were not parsed - please read the log file for more details." << std::endl;
		imp_file << "- Warning: some files were not parsed - please read the log file for more details." << std::endl;
	}

	// Produce a warning if an ensemble is set as the default reference frame
	if (isEpsgWGS84Ensemble(epsgCode))
	{
		std::stringstream ssEnsembleWarning;
		ssEnsembleWarning << std::endl <<
			"- Warning:  The '" << p.i.reference_frame << "' reference frame set for this project refers to the" << std::endl <<
			"  \"World Geodetic System 1984 (WGS 84) ensemble\".  The WGS 84 ensemble is" << std::endl <<
			"  only suitable for low accuracy (metre level) positioning and does not" << std::endl <<
			"  provide for precise transformations to other well-known reference frames." << std::endl <<
			"  To achieve reliable adjustment results from data on WGS 84, please refer" << std::endl <<
			"  to \"Configuring import options\" in the DynAdjust User's Guide." << std::endl;
		if (!p.g.quiet)
			std::cout << ssEnsembleWarning.str();
		imp_file << ssEnsembleWarning.str();
	}
	
	boost::posix_time::milliseconds elapsed_time(std::chrono::duration_cast<std::chrono::milliseconds>(time.elapsed().wall).count());
	std::string time_message = formatedElapsedTime<std::string>(&elapsed_time, "+ Total file handling process took ");

	if (!p.g.quiet)
		std::cout << std::endl << time_message << std::endl;
	imp_file << std::endl << time_message << std::endl;
	
	if (stnCount > 0 && msrCount > 0)
	{
		if (!p.g.quiet)
			std::cout << "+ Binary station and measurement files are now ready for processing." << std::endl << std::endl;
		imp_file << "+ Binary station and measurement files are now ready for processing." << std::endl << std::endl;
	}

	imp_file.close();
	return PARSE_SUCCESS;
}
