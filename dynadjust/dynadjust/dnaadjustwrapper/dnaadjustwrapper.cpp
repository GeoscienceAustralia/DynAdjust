//============================================================================
// Name         : dnaadjustwrapper.cpp
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

#include <dynadjust/dnaadjustwrapper/dnaadjustwrapper.hpp>

/// \cond
#include <time.h>
#include <algorithm>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <thread>
#include <vector>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
/// \endcond

#include <include/config/dnaoptions-interface.hpp>
#include <include/config/dnaprojectfile.hpp>
#include <include/config/dnatypes-gui.hpp>
#include <include/functions/dnafilepathfuncs.hpp>
#include <include/functions/dnastrmanipfuncs.hpp>
#include <include/functions/dnastrutils.hpp>
#include <include/parameters/dnaepsg.hpp>
#include <include/io/bst_file.hpp>
#include <include/io/bms_file.hpp>
#include <include/io/map_file.hpp>
#include <include/io/asl_file.hpp>
#include <dynadjust/dnaadjust/dnaadjust.hpp>
#include <dynadjust/dnaadjust/dnaadjust_printer.hpp>
#include <dynadjust/dnaadjustwrapper/dnaadjustprogress.hpp>
#include <dynadjust/dnaadjustwrapper/threading_init.hpp>

using namespace dynadjust::networkadjust;
using namespace dynadjust::exception;
using namespace dynadjust::iostreams;

extern bool running;
extern std::mutex cout_mutex;

using namespace dynadjust;
using namespace dynadjust::epsg;

// Helper function to convert std::filesystem::file_time_type to time_t
time_t file_time_to_time_t(const std::filesystem::file_time_type& ftime) {
    using namespace std::chrono;
    auto sctp = time_point_cast<system_clock::duration>(ftime - std::filesystem::file_time_type::clock::now() + system_clock::now());
    return system_clock::to_time_t(sctp);
}

void PrintSummaryMessage(dna_adjust* netAdjust, const project_settings* p, boost::posix_time::milliseconds *elapsed_time)
{
	if (p->g.quiet)
		return;

	cout_mutex.lock();
	std::cout.flush();
	UINT32 currentIteration(0);

	// any messages left
	while (netAdjust->NewMessagesAvailable())
	{
		if (!netAdjust->GetMessageIteration(currentIteration))
			break;
		std::stringstream ss("");
		ss << "  Iteration " << std::right << std::setw(2) << std::fixed << std::setprecision(0) << currentIteration;
		ss << ", max station corr: " << std::right << std::setw(12) << netAdjust->GetMaxCorrection(currentIteration) << std::endl;
		std::cout << PROGRESS_BACKSPACE_28 << std::setw(28) << std::left << ss.str();
	}

	if (p->a.report_mode)
	{
		std::cout << "+ Printing results of last adjustment only" << std::endl;
		cout_mutex.unlock();
		return;
	}

	std::cout << std::left << "+ Done." << std::endl;

	UINT32 block_count(netAdjust->blockCount());
	std::string block_str(" block");
	if (block_count > 1)
		block_str.append("s");
	block_str.append(".");
	
	switch (p->a.adjust_mode)
	{
	case Phased_Block_1Mode:
	case PhasedMode:
		if (netAdjust->GetStatus() == ADJUST_SUCCESS)
			std::cout << "+ Successfully adjusted " << block_count << block_str;
		else
			std::cout << "+ Attempted to adjust " << netAdjust->blockCount() << block_str;
		std::cout << std::endl;
	}
	
	
	std::cout << "+ Solution: ";

	if (netAdjust->GetStatus() != ADJUST_SUCCESS)
	{	
		std::cout << "failed to converge after ";
		if (p->a.adjust_mode == Phased_Block_1Mode ||
			p->a.max_iterations == 1)
			std::cout << "one iteration." << std::endl;
		else
			std::cout << p->a.max_iterations << " iterations." << std::endl;

		if (netAdjust->GetStatus() > ADJUST_THRESHOLD_EXCEEDED)
		{
			std::cout << std::endl << "+ Open " << leafStr<std::string>(p->o._adj_file) << " to view the adjustment details." << std::endl << std::endl;
			cout_mutex.unlock();
			return;
		}		
	}
	else
	{
		switch (p->a.adjust_mode)
		{
		case Phased_Block_1Mode:
			std::cout << "estimates solved for Block 1 only." << std::endl <<
			 	std::endl << 
				"- Warning: Depending on the quality of the apriori station estimates, further" << std::endl <<
				"  iterations may be needed. --block1-phased mode should only be used once" << std::endl <<
				"  rigorous estimates have been produced for the entire network." << std::endl << std::endl;
			break;
		default:
			std::cout << "converged after " << netAdjust->CurrentIteration() << " iteration"; 
			if (netAdjust->CurrentIteration() > 1)
				std::cout << "s";
			std::cout << "." << std::endl;
		}
	}

	std::cout << formatedElapsedTime<std::string>(elapsed_time, "+ Network adjustment took ") << std::endl;
	cout_mutex.unlock();
	
}

void SerialiseVarianceMatrices(dna_adjust* netAdjust, const project_settings* p)
{
	// No need to facilitate serialising if network adjustment is in stage,
	// as this will already be taken care of
	if (p->a.stage)
		return;

	if (!p->g.quiet)
	{
		std::cout << "+ Serialising adjustment matrices... ";
		std::cout.flush();
	}

	netAdjust->SerialiseAdjustedVarianceMatrices();

	if (!p->g.quiet)
	{
		std::cout << "done." << std::endl;
		std::cout.flush();
	}

	
}
	

void DeserialiseVarianceMatrices(dna_adjust* netAdjust, const project_settings* p)
{
	// No need to facilitate serialising if network adjustment is in stage,
	// as this will already be taken care of
	if (p->a.stage)
		return;

	netAdjust->DeSerialiseAdjustedVarianceMatrices();
}
	

void GenerateStatistics(dna_adjust* netAdjust, const project_settings* p)
{
	// Generate statistics
	// Don't produce statistics only for block 1 only adjustments
	if (p->a.adjust_mode != Phased_Block_1Mode)
	{
		if (!p->g.quiet)
		{
			std::cout << "+ Generating statistics...";
			std::cout.flush();
		}
		netAdjust->GenerateStatistics();
		if (!p->g.quiet)
		{
			std::cout << " done." << std::endl;

			std::cout << "+ Adjustment results:" << std::endl << std::endl;
			std::cout << "+" << OUTPUTLINE << std::endl;
			std::cout << std::setw(PRINT_VAR_PAD) << std::left << "  Number of unknown parameters" << std::fixed << std::setprecision(0) << netAdjust->GetUnknownsCount();
			if (netAdjust->GetAllFixed())
				std::cout << "  (All stations held constrained)";
			std::cout << std::endl;

			std::cout << std::setw(PRINT_VAR_PAD) << std::left << "  Number of measurements" << std::fixed << std::setprecision(0) << netAdjust->GetMeasurementCount();
			
			if (netAdjust->GetPotentialOutlierCount() > 0)
			{
				std::cout << "  (" << netAdjust->GetPotentialOutlierCount() << " potential outlier";
				if (netAdjust->GetPotentialOutlierCount() > 1)
					std::cout << "s";
				std::cout << ")";
			}
			std::cout << std::endl;
			std::cout << std::setw(PRINT_VAR_PAD) << std::left << "  Degrees of freedom" << std::fixed << std::setprecision(0) << netAdjust->GetDegreesOfFreedom() << std::endl;
			std::cout << std::setw(PRINT_VAR_PAD) << std::left << "  Chi squared" << std::fixed << std::setprecision(2) << netAdjust->GetChiSquared() << std::endl;
			std::cout << std::setw(PRINT_VAR_PAD) << std::left << "  Rigorous sigma zero" << std::fixed << std::setprecision(3) << netAdjust->GetSigmaZero() << std::endl;
			std::cout << std::setw(PRINT_VAR_PAD) << std::left << "  Global (Pelzer) Reliability" << std::fixed << std::setw(8) << std::setprecision(3) << netAdjust->GetGlobalPelzerRel() << "(excludes non redundant measurements)" << std::endl << std::endl;
		
			std::stringstream ss("");
			ss << std::left << "  Chi-Square test (" << std::setprecision(1) << std::fixed << p->a.confidence_interval << "%)";
			std::cout << std::setw(PRINT_VAR_PAD) << std::left << ss.str();
			ss.str("");
			ss << std::fixed << std::setprecision(3) << 
				netAdjust->GetChiSquaredLowerLimit() << " < " << 
				netAdjust->GetSigmaZero() << " < " <<
				netAdjust->GetChiSquaredUpperLimit();
			std::cout << std::setw(CHISQRLIMITS) << std::left << ss.str();
			ss.str("");

			if (netAdjust->GetDegreesOfFreedom() < 1)
				ss << "NO REDUNDANCY";
			else
			{
				ss << "*** ";
				switch (netAdjust->GetTestResult())
				{
				case test_stat_pass: 
					ss << "PASSED";		// within upper and lower
					break;
				case test_stat_warning:
					ss << "WARNING";	// less than lower limit
					break;
				case test_stat_fail:
					ss << "FAILED";		// greater than upper limit
					break;
				}
				ss << " ***";
			}

			std::cout << std::setw(PASS_FAIL) << std::right << ss.str() << std::endl;	
			std::cout << "+" << OUTPUTLINE << std::endl << std::endl;
		}
	}
	else
		std::cout << std::endl;
}

void PrintAdjustedMeasurements(dna_adjust* netAdjust, const project_settings* p)
{
	if (p->o._adj_msr_final)
	{
		if (!p->g.quiet)
		{
			std::cout << "+ Printing adjusted measurements...";
			std::cout.flush();
		}

		netAdjust->GetPrinter()->PrintAdjustedNetworkMeasurements();
		if (!p->g.quiet)
			std::cout << " done." << std::endl;
	}
}

void PrintMeasurementstoStations(dna_adjust* netAdjust, const project_settings* p)
{
	// Print measurements to stations table
	if (p->o._msr_to_stn)
	{
		if (!p->g.quiet)
		{
			std::cout << "+ Printing summary of measurements connected to each station...";
			std::cout.flush();
		}
		netAdjust->GetPrinter()->PrintMeasurementsToStation();
		if (!p->g.quiet)
			std::cout << " done." << std::endl;
	}
}

void PrintAdjustedNetworkStations(dna_adjust* netAdjust, const project_settings* p)
{
	// Print adjusted stations to ADJ file
	if (!p->g.quiet)
	{
		std::cout << "+ Printing adjusted station coordinates...";
		std::cout.flush();
	}
	netAdjust->GetPrinter()->PrintAdjustedNetworkStations();
	if (!p->g.quiet)
		std::cout << " done." << std::endl;
}

void PrintPositionalUncertainty(dna_adjust* netAdjust, const project_settings* p)
{
	// Print positional uncertainty
	if (p->o._positional_uncertainty)
	{
		if (!p->g.quiet)
		{
			std::cout << "+ Printing positional uncertainty of adjusted coordinates...";
			std::cout.flush();
		}
		// Print correlations as required
		netAdjust->GetPrinter()->PrintPositionalUncertainty();
		if (!p->g.quiet)
			std::cout << " done." << std::endl;
	}
}

void PrintStationCorrections(dna_adjust* netAdjust, const project_settings* p)
{
	// Print corrections
	if (p->o._init_stn_corrections)
	{
		if (!p->g.quiet)
		{
			std::cout << "+ Printing corrections to initial station coordinates...";
			std::cout.flush();
		}
		netAdjust->GetPrinter()->PrintNetworkStationCorrections();
		if (!p->g.quiet)
			std::cout << " done." << std::endl;
	}
}

void UpdateBinaryFiles(dna_adjust* netAdjust, const project_settings* p)
{
	// Update bst and bms files with adjustment results
	if (!p->g.quiet)
	{
		std::cout << "+ Updating binary station and measurement files...";
		std::cout.flush();
	}
	netAdjust->UpdateBinaryFiles();
	if (!p->g.quiet)
		std::cout << " done." << std::endl;
}

void ExportDynaML(dna_adjust* netAdjust, project_settings* p)
{
	// Output adjustment as XML stn
	if (p->o._export_xml_stn_file)
	{
		// single file for both stations and measurements
		p->o._xml_file = p->o._adj_file + ".stn.xml";
				
		if (!p->g.quiet)
			std::cout << "+ Serializing estimated coordinates to " << leafStr<std::string>(p->o._xml_file) << "... ";
				
		// Export Stations file
		netAdjust->GetPrinter()->PrintEstimatedStationCoordinatestoDNAXML(p->o._xml_file, dynaml, 
			(p->i.flag_unused_stn ? true : false));

		if (!p->g.quiet)
			std::cout << "Done." << std::endl;
	}

	// Output adjustment as XML msr
	if (p->o._export_xml_msr_file)
	{
		// single file for both stations and measurements
		p->o._xml_file = p->o._adj_file + ".msr.xml";
				
		if (!p->g.quiet)
			std::cout << "+ Serializing estimated coordinates and uncertainties to " << leafStr<std::string>(p->o._xml_file) << "... ";
				
		// Export Measurements file (exclude unused stations given 
		// they will not have been estimated)
		netAdjust->GetPrinter()->PrintEstimatedStationCoordinatestoDNAXML_Y(p->o._xml_file, dynaml);

		if (!p->g.quiet)
			std::cout << "Done." << std::endl;
	}
}

void ExportDNA(dna_adjust* netAdjust, project_settings* p)
{
	// Print adjusted stations and measurements to DNA stn
	if (p->o._export_dna_stn_file)
	{
		std::string stnfilename(p->o._adj_file + ".stn");
		
		if (!p->g.quiet)
			std::cout << "+ Serializing estimated coordinates to " << leafStr<std::string>(stnfilename) << "... ";
					
		// Export Station file
		netAdjust->GetPrinter()->PrintEstimatedStationCoordinatestoDNAXML(stnfilename, dna, 
			(p->i.flag_unused_stn ? true : false));

		if (!p->g.quiet)
			std::cout << "Done." << std::endl;
	}

	// Print adjusted stations and measurements to DNA msr
	if (p->o._export_dna_msr_file)
	{
		std::string msrfilename(p->o._adj_file + ".msr");
		
		if (!p->g.quiet)
			std::cout << "+ Serializing estimated coordinates and uncertainties to " << leafStr<std::string>(msrfilename) << "... ";
					
		// Export Measurements file (exclude unused stations given 
		// they will not have been estimated)
		netAdjust->GetPrinter()->PrintEstimatedStationCoordinatestoDNAXML_Y(msrfilename, dna);

		if (!p->g.quiet)
			std::cout << "Done." << std::endl;
	}
}

void ExportSinex(dna_adjust* netAdjust, const project_settings* p)
{
	// Print adjusted stations and measurements to SINEX
	if (p->o._export_snx_file)
	{
		std::string sinex_file;
		// Export to SINEX
		if (!p->g.quiet)
			std::cout << "+ Printing station estimates and uncertainties to SINEX...";
		bool success(netAdjust->GetPrinter()->PrintEstimatedStationCoordinatestoSNX(sinex_file));

		// SomeFunc()
		if (!p->g.quiet)
			std::cout << " done." << std::endl;

		if (!success)
		{
			std::cout << "- Warning: The SINEX export process produced some warnings." << std::endl;
			switch (p->a.adjust_mode)
			{
			case PhasedMode:
				sinex_file = findandreplace(sinex_file, std::string("-block1"), std::string("-block*"));
			}

			std::cout << "  See " << leafStr<std::string>(sinex_file) << ".err for details." << std::endl; 
		}
	}
}

int ParseCommandLineOptions(const int& argc, char* argv[], const boost::program_options::variables_map& vm, project_settings& p)
{
	// capture command line arguments
	for (int cmd_arg(0); cmd_arg<argc; ++cmd_arg)
	{
		 p.a.command_line_arguments += argv[cmd_arg];
		 p.a.command_line_arguments += " ";
	}

	if (vm.count(PROJECT_FILE))
	{
		if (std::filesystem::exists(p.g.project_file))
		{
			try {
				CDnaProjectFile projectFile(p.g.project_file, adjustSetting);
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

	if (std::filesystem::exists(p.g.project_file))
	{
		// update import settings from dnaproj file
		try {
			CDnaProjectFile projectFile(p.g.project_file, importSetting);
			p.i = projectFile.GetSettings().i;
		}
		catch (...) {
			// do nothing
		}

		// update geoid file name from dnaproj file (blank if geoid was not executed)
		try {
			CDnaProjectFile projectFile(p.g.project_file, geoidSetting);
			p.n = projectFile.GetSettings().n;
		}
		catch (...) {
			// do nothing
		}

		// update reftran settings from dnaproj file
		// Note, if reftran was not executed, reference frame and epoch will be set to 
		// the frame and epoch captured on import.
		try {
			CDnaProjectFile projectFile(p.g.project_file, reftranSetting);
			p.r = projectFile.GetSettings().r;
		}
		catch (...) {
			// do nothing
		}
	}

	// binary station file location (output)
	if (vm.count(BIN_STN_FILE))
		p.a.bst_file = formPath<std::string>(p.g.input_folder, p.a.bst_file);
	else
		p.a.bst_file = formPath<std::string>(p.g.output_folder, p.g.network_name, "bst");
	
	// binary station file location (output)
	if (vm.count(BIN_MSR_FILE))
		p.a.bms_file = formPath<std::string>(p.g.input_folder, p.a.bms_file);
	else
		p.a.bms_file = formPath<std::string>(p.g.output_folder, p.g.network_name, "bms");

	if (!std::filesystem::exists(p.a.bst_file) || !std::filesystem::exists(p.a.bms_file))
	{
		cout_mutex.lock();
		std::cout << std::endl << "- Nothing to do: ";  
			
		if (p.g.network_name.empty())
			std::cout << std::endl << "network name has not been specified specified, and " << std::endl << "               ";  
		std::cout << p.a.bst_file << " and " << p.a.bms_file << " do not exist." << std::endl << std::endl;  
		cout_mutex.unlock();
		return EXIT_FAILURE;
	}

	// output settings
	if (vm.count(OUTPUT_ADJ_STN_ITER))
		p.o._adj_stn_iteration = 1;
	if (vm.count(OUTPUT_ADJ_MSR_ITER))
		p.o._adj_msr_iteration = 1;
	if (vm.count(OUTPUT_CMP_MSR_ITER))
		p.o._cmp_msr_iteration = 1;
	if (vm.count(OUTPUT_ADJ_STAT_ITER))
		p.o._adj_stat_iteration = 1;
	if (vm.count(OUTPUT_ADJ_MSR) ||				// print adjusted measurements?
		vm.count(OUTPUT_ADJ_GNSS_UNITS) ||		// print alternative units for adjusted GNSS measurements?
		vm.count(OUTPUT_ADJ_MSR_SORTBY) ||		// print sort adjusted measurements?
		vm.count(OUTPUT_ADJ_MSR_TSTAT))			// print t-statistic for adjusted measurements?
		p.o._adj_msr_final = 1;

	if (vm.count(OUTPUT_ADJ_STN_BLOCKS))
		p.o._output_stn_blocks = 1;
	if (vm.count(OUTPUT_ADJ_MSR_BLOCKS))
		p.o._output_msr_blocks = 1;
	if (vm.count(OUTPUT_ADJ_STN_SORT_ORDER))
		p.o._sort_stn_file_order = 1;
	
	if (vm.count(MODE_SIMULTANEOUS))
		p.a.adjust_mode = SimultaneousMode;		// default
	else if (vm.count(MODE_PHASED_BLOCK1))
		p.a.adjust_mode = Phased_Block_1Mode;
	else if (vm.count(MODE_PHASED_MT))
	{
		p.a.multi_thread = 1;
		p.a.adjust_mode = PhasedMode;
	}
	else if (vm.count(MODE_PHASED))
		p.a.adjust_mode = PhasedMode;
	else if (vm.count(MODE_SIMULATION))
		p.a.adjust_mode = SimulationMode;

	// Report mode?
	if (vm.count(MODE_ADJ_REPORT))
	{
		p.a.report_mode = true;
		p.a.max_iterations = 0;
	}
	else if (p.a.max_iterations < 1)
		p.a.report_mode = true;

	if (vm.count(STAGED_ADJUSTMENT))
	{
		p.a.stage = true;
		p.a.multi_thread = false;
		p.a.adjust_mode = PhasedMode;
		//p.o._output_stn_blocks = true;
	}

	// Force inverse method for measurement variances to be the same as that which
	// is used for the inversion of the normals
	//if (vm.count(LSQ_INVERSE_METHOD))
	//	p.a.inverse_method_msr = p.a.inverse_method_lsq;
	if (vm.count(SCALE_NORMAL_UNITY))
		p.a.scale_normals_to_unity = 1;
	if (vm.count(OUTPUT_ADJ_MSR_TSTAT))
		p.o._adj_msr_tstat = 1;
	if (vm.count(OUTPUT_ADJ_MSR_DBID))
		p.o._database_ids = 1;
	if (vm.count(OUTPUT_IGNORED_MSRS))
		p.o._print_ignored_msrs = 1;
	if (vm.count(PURGE_STAGE_FILES))
		p.a.purge_stage_files = 1;
	if (vm.count(RECREATE_STAGE_FILES))
		p.a.recreate_stage_files = 1;
	if (vm.count(TYPE_B_GLOBAL))
		p.o._apply_type_b_global = 1;
	if (vm.count(TYPE_B_FILE))
		p.o._apply_type_b_file = 1;
	p.s.asl_file = formPath<std::string>(p.g.output_folder, p.g.network_name, "asl");	// associated stations list
	p.s.aml_file = formPath<std::string>(p.g.output_folder, p.g.network_name, "aml");	// associated measurements list
	p.a.map_file = formPath<std::string>(p.g.output_folder, p.g.network_name, "map");	// station names map
	
	// has a seg file name been specified?
	if (vm.count(SEG_FILE))
		p.a.seg_file = formPath<std::string>(p.g.input_folder, p.a.seg_file);
	else
		p.a.seg_file = formPath<std::string>(p.g.output_folder, p.g.network_name, "seg");
	
	if (vm.count(OUTPUT_APU_CORRELATIONS))
		p.o._output_pu_covariances = 1;

	// Set up file names dependent on adjustment mode
	p.o._xyz_file = p.o._adj_file = 
		formPath<std::string>(p.g.output_folder, p.g.network_name);

	if (vm.count(OUTPUT_POS_UNCERTAINTY))
	{
		p.o._positional_uncertainty = 1;
		p.o._apu_file = p.o._adj_file;
	}

	if (vm.count(OUTPUT_STN_COR_FILE))
		p.o._cor_file = p.o._adj_file;

	switch (p.a.adjust_mode)
	{
	case Phased_Block_1Mode:
	case PhasedMode:

		p.o._adj_file += ".phased";
		p.o._xyz_file += ".phased";

		if (vm.count(OUTPUT_POS_UNCERTAINTY))
			p.o._apu_file += ".phased";

		if (vm.count(OUTPUT_STN_COR_FILE))
			p.o._cor_file += ".phased";
		
		if (p.a.adjust_mode == Phased_Block_1Mode)
		{
			p.o._adj_file += "-block1";
			p.o._xyz_file += "-block1";
			
			if (vm.count(OUTPUT_POS_UNCERTAINTY))
				p.o._apu_file += "-block1";

			if (vm.count(OUTPUT_STN_COR_FILE))
				p.o._cor_file += "-block1";

		}
		else if (p.a.stage)
		{
			p.o._adj_file += "-stage";
			p.o._xyz_file += "-stage";

			if (vm.count(OUTPUT_POS_UNCERTAINTY))
				p.o._apu_file += "-stage";

			if (vm.count(OUTPUT_STN_COR_FILE))
				p.o._cor_file += "-stage";

		}
		else if (p.a.multi_thread)
		{
			p.o._adj_file += "-mt";
			p.o._xyz_file += "-mt";
			
			if (vm.count(OUTPUT_POS_UNCERTAINTY))
				p.o._apu_file += "-mt";

			if (vm.count(OUTPUT_STN_COR_FILE))
				p.o._cor_file += "-mt";
		}
		break;
	case SimultaneousMode:
		p.o._adj_file += ".simult";
		p.o._xyz_file += ".simult";
		
		if (vm.count(OUTPUT_POS_UNCERTAINTY))
			p.o._apu_file += ".simult";
		if (vm.count(OUTPUT_STN_COR_FILE))
			p.o._cor_file += ".simult";
		break;
	}

	p.o._adj_file += ".adj";
	p.o._xyz_file += ".xyz";

	if (vm.count(OUTPUT_POS_UNCERTAINTY))
		p.o._apu_file += ".apu";

	if (vm.count(OUTPUT_STN_COR_FILE))
		p.o._cor_file += ".cor";

	if (vm.count(OUTPUT_STN_COR_FILE))
		p.o._init_stn_corrections = 1;

	if (vm.count(OUTPUT_STN_CORR))
		p.o._stn_corr = 1;

	if (vm.count(OUTPUT_MSR_TO_STN))
		p.o._msr_to_stn = 1;
	
	if (vm.count(TEST_INTEGRITY))
		p.i.test_integrity = 1;

	if (vm.count(EXPORT_XML_STN_FILE))
		p.o._export_xml_stn_file = 1;

	if (vm.count(EXPORT_XML_MSR_FILE))
		p.o._export_xml_msr_file = 1;

	if (vm.count(EXPORT_DNA_STN_FILE))
		p.o._export_dna_stn_file = 1;

	if (vm.count(EXPORT_DNA_MSR_FILE))
		p.o._export_dna_msr_file = 1;

	if (vm.count(EXPORT_SNX_FILE))
		p.o._export_snx_file = 1;

	return EXIT_SUCCESS;
}

void LoadBinaryMeta(binary_file_meta_t& bst_meta, binary_file_meta_t& bms_meta,
	const project_settings& p, bool& bst_meta_import, bool& bms_meta_import)
{
	BstFile bst;
	BmsFile bms;

	bst.LoadFileMeta(p.a.bst_file, bst_meta);
	bms.LoadFileMeta(p.a.bms_file, bms_meta);

	bst_meta_import = (iequals(bst_meta.modifiedBy, __import_app_name__) ||
		iequals(bst_meta.modifiedBy, __import_dll_name__));
	bms_meta_import = (iequals(bms_meta.modifiedBy, __import_app_name__) ||
		iequals(bms_meta.modifiedBy, __import_dll_name__));
}

int main(int argc, char* argv[])
{
	// create banner message
	std::string cmd_line_banner, stnfilename, msrfilename;	
	fileproc_help_header(&cmd_line_banner);

	project_settings p;

	boost::program_options::variables_map vm;
	boost::program_options::positional_options_description positional_options;

	boost::program_options::options_description standard_options("+ " + std::string(ALL_MODULE_STDOPT), PROGRAM_OPTIONS_LINE_LENGTH);
	boost::program_options::options_description adj_mode_options("+ " + std::string(ADJUST_MODULE_MODE), PROGRAM_OPTIONS_LINE_LENGTH);
	boost::program_options::options_description phased_adj_options("+ " + std::string(ADJUST_MODULE_PHASED), PROGRAM_OPTIONS_LINE_LENGTH);
	boost::program_options::options_description adj_config_options("+ " + std::string(ADJUST_MODULE_CONFIG), PROGRAM_OPTIONS_LINE_LENGTH);
	boost::program_options::options_description staged_adj_options("+ " + std::string(ADJUST_MODULE_STAGE), PROGRAM_OPTIONS_LINE_LENGTH);
	boost::program_options::options_description output_options("+ " + std::string(ALL_MODULE_OUTPUT), PROGRAM_OPTIONS_LINE_LENGTH);
	boost::program_options::options_description export_options("+ " + std::string(ALL_MODULE_EXPORT), PROGRAM_OPTIONS_LINE_LENGTH);
	boost::program_options::options_description generic_options("+ " + std::string(ALL_MODULE_GENERIC), PROGRAM_OPTIONS_LINE_LENGTH);

	std::string cmd_line_usage("+ ");
	cmd_line_usage.append(__BINARY_NAME__).append(" usage:  ").append(__BINARY_NAME__).append(" ").append(NETWORK_NAME).append(" [options]");
	boost::program_options::options_description allowable_options(cmd_line_usage, PROGRAM_OPTIONS_LINE_LENGTH);
	
	try {
		standard_options.add_options()
			(PROJECT_FILE_P, boost::program_options::value<std::string>(&p.g.project_file),
				"Project file containing all user options. If specified, all other options are ignored.")
			(NETWORK_NAME_N, boost::program_options::value<std::string>(&p.g.network_name),
				"Network name. User defined name for all input and output files. Default is \"network#\".")
			(INPUT_FOLDER_I, boost::program_options::value<std::string>(&p.g.input_folder),
				"Path containing all input files")
			(OUTPUT_FOLDER_O, boost::program_options::value<std::string>(&p.g.output_folder),		// default is ./,
				"Path for all output files")
			(BIN_STN_FILE, boost::program_options::value<std::string>(&p.a.bst_file),
				"Binary station file name. Overrides network name.")
			(BIN_MSR_FILE, boost::program_options::value<std::string>(&p.a.bms_file),
				"Binary measurement file name. Overrides network name.")
			(SEG_FILE, boost::program_options::value<std::string>(&p.a.seg_file),
				"Network segmentation file name. Overrides network name.")
			(COMMENTS, boost::program_options::value<std::string>(&p.a.comments),
				"Comments about the adjustment. All comments are printed to the adj file.")
			;

		adj_mode_options.add_options()
			(MODE_SIMULTANEOUS,
				"Simultaneous adjustment mode. The default mode.")
			(MODE_PHASED,
				"Sequential phased adjustment mode.")
			//(MODE_SIMULATION,
			//	"Adjustment simulation mode.")
			(MODE_ADJ_REPORT,
				"Reproduce the adjustment output files without performing an adjustment.")
			;

		phased_adj_options.add_options()
			(STAGED_ADJUSTMENT,
				"Store adjustment matrices in memory mapped files instead of retaining data in memory.  This option decreases efficiency but may be required if there is insufficient RAM to hold an adjustment in memory.")
			(MODE_PHASED_MT,
				"Process forward, reverse and combination adjustments concurrently using all available CPU cores.")
			(MODE_PHASED_BLOCK1,
				"Sequential phased adjustment mode resulting in rigorous estimates for block 1 only.")
			;
		
		adj_config_options.add_options()
			(CONF_INTERVAL, boost::program_options::value<float>(&p.a.confidence_interval),
				(std::string("Confidence interval for testing the least squares solution and measurement corrections. Default is ")+
				StringFromT(p.a.confidence_interval, 1)+std::string("%.")).c_str())
			(ITERATION_THRESHOLD, boost::program_options::value<float>(&p.a.iteration_threshold),
				(std::string("Least squares iteration threshold. Default is ")+
				StringFromT(p.a.iteration_threshold, 4)+std::string("m.")).c_str())
			(MAX_ITERATIONS, boost::program_options::value<UINT16>(&p.a.max_iterations),
				(std::string("Maximum number of iterations. Default is ")+
				StringFromT(p.a.max_iterations)+std::string(".")).c_str())
			(STN_CONSTRAINTS, boost::program_options::value<std::string>(&p.a.station_constraints),
				"Station constraints. arg is a comma delimited string \"stn1,CCC,stn2,CCF\" defining specific station constraints. These constraints override those contained in the station file.")
			(FREE_STN_SD, boost::program_options::value<double>(&p.a.free_std_dev),
				(std::string("A-priori standard deviation for free stations. Default is ")+
				StringFromT(p.a.free_std_dev)+std::string("m.")).c_str())
			(FIXED_STN_SD, boost::program_options::value<double>(&p.a.fixed_std_dev),
				(std::string("A-priori standard deviation for fixed stations. Default is ")+
				StringFromT(p.a.fixed_std_dev, 6)+std::string("m.")).c_str())
			(SCALE_NORMAL_UNITY,
				"Scale adjustment normal matrices to unity prior to computing inverse to minimise loss of precision caused by tight variances placed on constraint stations.")
			(TYPE_B_GLOBAL, boost::program_options::value<std::string>(&p.a.type_b_global),
				"Type b uncertainties to be added to each computed uncertainty. arg is a comma delimited string that provides 1D, 2D or 3D uncertainties in the local reference frame (e.g. \"up\" or \"e,n\" or \"e,n,up\").")
			(TYPE_B_FILE, boost::program_options::value<std::string>(&p.a.type_b_file),
				"Type b uncertainties file name. Full path to a file containing Type b uncertainties to be added to the computed uncertainty for specific sites.")
			(RELAXATION, boost::program_options::value<double>(&p.a.relaxation),
				"Initial relaxation factor (0 < alpha <= 1). Adaptive oscillation-aware damping: automatically reduced when sign-alternating oscillation is detected, slowly recovered during monotonic convergence. Default is 1.0 (start with full Newton steps).")
			;

		staged_adj_options.add_options()
			(RECREATE_STAGE_FILES,
				"Recreate memory mapped files.")
			(PURGE_STAGE_FILES,
				"Purge memory mapped files from disk upon adjustment completion.")
			;

		output_options.add_options()
			(OUTPUT_MSR_TO_STN,
				"Output summary of measurements connected to each station.")
			(OUTPUT_MSR_TO_STN_SORTBY, boost::program_options::value<UINT16>(&p.o._sort_msr_to_stn),
				std::string("Sort order for measurement to stations summary.\n  " +
					StringFromT(orig_stn_sort_ui) + ": Original station order (default)\n  " +
					StringFromT(name_stn_sort_ui) + ": Station name\n  " +
					StringFromT(meas_stn_sort_ui) + ": Measurement count (ascending)\n  " +
					StringFromT(saem_stn_sort_ui) + ": Measurement count (descending)").c_str())
			(OUTPUT_ADJ_STN_ITER,
				"Output adjusted station coordinates on each iteration.")
			(OUTPUT_ADJ_STAT_ITER,
				"Output statistical summary on each iteration.")
			(OUTPUT_ADJ_MSR_ITER,
				"Output adjusted measurements on each iteration.")
			(OUTPUT_CMP_MSR_ITER,
				"Output computed measurements on each iteration.")
			(OUTPUT_ADJ_MSR,
				"Output final adjusted measurements.")
			(OUTPUT_ADJ_GNSS_UNITS, boost::program_options::value<UINT16>(&p.o._adj_gnss_units),
				std::string("Units for adjusted GNSS baseline measurements in the .adj file.\n  " + 
					StringFromT(XYZ_adj_gnss_ui) + ": As measured (default)\n  " + 
					StringFromT(ENU_adj_gnss_ui) + ": Local [east, north, up]\n  " + 
					StringFromT(AED_adj_gnss_ui) + ": Polar [azimuth, vert. angle, slope dist]\n  " + 
					StringFromT(ADU_adj_gnss_ui) + ": Polar [azimuth, slope dist, up]").c_str())
			(OUTPUT_ADJ_MSR_TSTAT,
				"Output t-statistics for adjusted measurements.")
			(OUTPUT_ADJ_MSR_DBID,
				"Output measurement and cluster ids for database mapping.")
			(OUTPUT_IGNORED_MSRS,
				"Output adjusted measurement statistics for ignored measurements.")
			(OUTPUT_ADJ_MSR_SORTBY, boost::program_options::value<UINT16>(&p.o._sort_adj_msr),
				std::string("Sort order for adjusted measurements.\n  " + 
					StringFromT(orig_adj_msr_sort_ui) + ": Original input file order (default)\n  " + 
					StringFromT(type_adj_msr_sort_ui) + ": Measurement type\n  " + 
					StringFromT(inst_adj_msr_sort_ui) + ": Station 1\n  " + 
					StringFromT(targ_adj_msr_sort_ui) + ": Station 2\n  " + 
					StringFromT(meas_adj_msr_sort_ui) + ": Measurement value\n  " + 
					StringFromT(corr_adj_msr_sort_ui) + ": Correction\n  " + 
					StringFromT(a_sd_adj_msr_sort_ui) + ": Adjusted std. dev.\n  " + 
					StringFromT(n_st_adj_msr_sort_ui) + ": N-statistic").c_str())
			(OUTPUT_ADJ_STN_BLOCKS,
				"For phased adjustments, output adjusted coordinates according to each block.")
			(OUTPUT_ADJ_MSR_BLOCKS,
				"For phased adjustments, output adjusted measurements according to each block.")
			(OUTPUT_ADJ_STN_SORT_ORDER,
				"Output station information using the station order in the original station file. By default, stations are output in alpha-numeric order.")
			(OUTPUT_STN_COORD_TYPES, boost::program_options::value<std::string>(&p.o._stn_coord_types),
				(std::string("Output station coordinate types. arg is a case-sensitive string of chars \"ENzPLHhXYZ\" defining the specific types to be printed. Default is ").append(
				p.o._stn_coord_types).append(
				".")).c_str())
			(OUTPUT_ANGULAR_TYPE_STN, boost::program_options::value<UINT16>(&p.o._angular_type_stn),
				std::string("Output type for angular station coordinates.\n"
					"  0: Degrees, minutes and seconds (default)\n"
					"  1: Decimal degrees").c_str())
			(OUTPUT_STN_CORR,
				"Output station corrections with adjusted station coordinates.")
			(OUTPUT_PRECISION_METRES_STN, boost::program_options::value<UINT16>(&p.o._precision_metres_stn),
				(std::string("Output precision for linear station coordinates in metres. Default is ")+
				StringFromT(p.o._precision_metres_stn, 0)).c_str())
			(OUTPUT_PRECISION_SECONDS_STN, boost::program_options::value<UINT16>(&p.o._precision_seconds_stn),
				(std::string("Output precision for angular station coordinates. For values in degrees, minutes and seconds, precision relates to seconds. For values in decimal degrees, precision relates to degrees. Default is ")+
				StringFromT(p.o._precision_seconds_stn, 0)).c_str())
			(OUTPUT_PRECISION_METRES_MSR, boost::program_options::value<UINT16>(&p.o._precision_metres_msr),
				(std::string("Output precision for linear measurements in metres. Default is ")+
				StringFromT(p.o._precision_metres_msr, 0)).c_str())
			(OUTPUT_PRECISION_SECONDS_MSR, boost::program_options::value<UINT16>(&p.o._precision_seconds_msr),
				(std::string("Output precision for angular measurements. For values in degrees, minutes and seconds, precision relates to seconds. For values in decimal degrees, precision relates to degrees. Default is ")+
				StringFromT(p.o._precision_seconds_msr, 0)).c_str())
			(OUTPUT_ANGULAR_TYPE_MSR, boost::program_options::value<UINT16>(&p.o._angular_type_msr),
				std::string("Output type for angular measurements.\n"
					"  0: Degrees, minutes and seconds (default)\n"
					"  1: Decimal degrees").c_str())
			(OUTPUT_DMS_FORMAT_MSR, boost::program_options::value<UINT16>(&p.o._dms_format_msr),
				std::string("Output format for angular (dms) measurements.\n"
					"  0: Separated fields (default)\n"
					"  1: Separated fields with symbols\n"
					"  2: HP notation").c_str())
			;

		export_options.add_options()
			(OUTPUT_POS_UNCERTAINTY,
				"Output positional uncertainty and variances of adjusted station coordinates to .apu file.")
			(OUTPUT_APU_CORRELATIONS,
				"Output covariances between adjusted station coordinates to .apu file.")
			(OUTPUT_APU_UNITS, boost::program_options::value<UINT16>(&p.o._apu_vcv_units),
				std::string("Variance matrix units in the .apu file.\n  " +
					StringFromT(XYZ_apu_ui) + ": Cartesian [X,Y,Z] (default)\n  " + 
					//StringFromT(LLH_apu_ui) + ": Geographic [Lat,Lon,ht]\n  " + 
					StringFromT(ENU_apu_ui) + ": Local [e,n,up]").c_str())
			(OUTPUT_STN_COR_FILE,
				"Output corrections (azimuth, distance, e, n, up) to initial station coordinates to .cor file.")
			(HZ_CORR_THRESHOLD, boost::program_options::value<double>(&p.o._hz_corr_threshold),
				(std::string("Minimum horizontal threshold by which to restrict output of station corrections to .cor file. Default is ")+
				StringFromT(p.o._hz_corr_threshold, 1)+std::string("m")).c_str())
			(VT_CORR_THRESHOLD, boost::program_options::value<double>(&p.o._vt_corr_threshold),
				(std::string("Minimum vertical threshold by which to restrict output of station corrections to .cor file. Default is ")+
				StringFromT(p.o._vt_corr_threshold, 1)+std::string("m")).c_str())
			//(UPDATE_ORIGINAL_STN_FILE,
			//	"Update original station file with adjusted station coordinates.")
			(EXPORT_XML_STN_FILE,
				"Export estimated station coordinates to DynaML (DynAdjust XML) station file.")
			(EXPORT_XML_MSR_FILE,
				"Export estimated station coordinates and uncertainties to DynaML (DynAdjust XML) measurement file as a GNSS Y cluster.")
			(EXPORT_DNA_STN_FILE,
				"Export estimated station coordinates to DNA station file.")
			(EXPORT_DNA_MSR_FILE,
				"Export estimated station coordinates and uncertainties to DNA measurement file as a GNSS Y cluster.")
			(EXPORT_SNX_FILE,
				"Export estimated station coordinates and full variance matrix to SINEX file. Note: station names will be truncated to four characters as per the SINEX standard.")
			;

		// Declare a group of options that will be 
		// allowed only on command line		
		generic_options.add_options()
			(VERBOSE, boost::program_options::value<UINT16>(&p.g.verbose),
				std::string("Give detailed information about what ").append(__BINARY_NAME__).append(" is doing.\n"
					"  0: No information (default)\n"
					"  1: Helpful information\n"
					"  2: Extended information\n"
					"  3: Debug level information").c_str())
			(QUIET,
				std::string("Suppresses all explanation of what ").append(__BINARY_NAME__).append(" is doing unless an error occurs").c_str())
			(VERSION_V, "Display the current program version")
			(HELP_H, "Show this help message")
			(HELP_MODULE_H, boost::program_options::value<std::string>(),
				"Provide help for a specific help category.")
			;

		allowable_options.add(standard_options).add(adj_mode_options).add(phased_adj_options).add(adj_config_options).add(staged_adj_options).add(output_options).add(export_options).add(generic_options);

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
		bool module_found(false);

		if (str_upper<std::string, char>(ALL_MODULE_STDOPT).find(help_text) != std::string::npos) {
			std::cout << standard_options << std::endl;
			module_found = true;
		}

		if (str_upper<std::string, char>(ADJUST_MODULE_MODE).find(help_text) != std::string::npos) {
			std::cout << adj_mode_options << std::endl;
			module_found = true;
		} 

		if (str_upper<std::string, char>(ADJUST_MODULE_PHASED).find(help_text) != std::string::npos) {
			std::cout << phased_adj_options << std::endl;
			module_found = true;
		} 

		if (str_upper<std::string, char>(ADJUST_MODULE_CONFIG).find(help_text) != std::string::npos) {
			std::cout << adj_config_options << std::endl;
			module_found = true;
		} 

		if (str_upper<std::string, char>(ADJUST_MODULE_STAGE).find(help_text) != std::string::npos) {
			std::cout << staged_adj_options << std::endl;
			module_found = true;
		} 

		if (str_upper<std::string, char>(ALL_MODULE_OUTPUT).find(help_text) != std::string::npos) {
			std::cout << output_options << std::endl;
			module_found = true;
		} 

		if (str_upper<std::string, char>(ALL_MODULE_EXPORT).find(help_text) != std::string::npos) {
			std::cout << export_options << std::endl;
			module_found = true;
		}

		if (str_upper<std::string, char>(ALL_MODULE_GENERIC).find(help_text) != std::string::npos) {
			std::cout << generic_options << std::endl;
			module_found = true;
		}

		if (!module_found) {
			std::cout << std::endl << "- Error: Help module '" <<
				original_text << "' is not in the list of options." << std::endl;
			return EXIT_FAILURE;
		}

		return EXIT_SUCCESS;
	}

	bool userSuppliedSegFile(false);
	if (!p.a.seg_file.empty())
		userSuppliedSegFile = true;
	bool userSuppliedBstFile(false);
	if (!p.a.bst_file.empty())
		userSuppliedBstFile = true;
	bool userSuppliedBmsFile(false);
	if (!p.a.bms_file.empty())
		userSuppliedBmsFile = true;

	if (ParseCommandLineOptions(argc, argv, vm, p) != EXIT_SUCCESS)
		return EXIT_FAILURE;

	// Create an instance of the dna_adjust object exposed by the dnaadjust dll
	dna_adjust netAdjust;

	// Capture binary file metadata
	binary_file_meta_t bst_meta, bms_meta;
	bool bst_meta_import, bms_meta_import;
	LoadBinaryMeta(bst_meta, bms_meta, p, bst_meta_import, bms_meta_import);

	// Capture datum set within project file
	CDnaDatum datum;
	
	// Inspect if reftran has been executed.  If so, select the appropriate 
	// reference frame label
	if (bst_meta.reftran)
		datum.SetDatumFromName(p.r.reference_frame, p.r.epoch);
	else
		datum.SetDatumFromName(p.i.reference_frame, p.i.epoch);

	if (vm.count(QUIET))
		p.g.quiet = 1;
	
	if (!p.g.quiet)
	{
		cout_mutex.lock();
		std::cout << std::endl << cmd_line_banner;

		std::cout << "+ Options:" << std::endl;
		std::cout << std::setw(PRINT_VAR_PAD) << std::left << "  Network name: " <<  p.g.network_name << std::endl;
		std::cout << std::setw(PRINT_VAR_PAD) << std::left << "  Input folder: " << p.g.input_folder << std::endl;
		std::cout << std::setw(PRINT_VAR_PAD) << std::left << "  Output folder: " << p.g.output_folder << std::endl;
		std::cout << std::setw(PRINT_VAR_PAD) << std::left << "  Associated station file: " << p.s.asl_file << std::endl;
		std::cout << std::setw(PRINT_VAR_PAD) << std::left << "  Associated measurement file: " << p.s.aml_file << std::endl;
		std::cout << std::setw(PRINT_VAR_PAD) << std::left << "  Binary station file: " << p.a.bst_file << std::endl;
		std::cout << std::setw(PRINT_VAR_PAD) << std::left << "  Binary measurement file: " << p.a.bms_file << std::endl;
		if (p.a.adjust_mode == PhasedMode || p.a.adjust_mode == Phased_Block_1Mode)
			std::cout << std::setw(PRINT_VAR_PAD) << std::left << "  Segmentation file: " << p.a.seg_file << std::endl;
		std::cout << std::setw(PRINT_VAR_PAD) << std::left << "  Adjustment output file: " << p.o._adj_file << std::endl;
		std::cout << std::setw(PRINT_VAR_PAD) << std::left << "  Coordinate output file: " << p.o._xyz_file << std::endl;
		if (p.o._init_stn_corrections)
			std::cout << std::setw(PRINT_VAR_PAD) << std::left << "  Corrections output file: " << p.o._cor_file << std::endl;
		
		if (p.a.stage)
		{
			std::cout << std::setw(PRINT_VAR_PAD) << std::left << "  Stage using hard disk: " << "yes" << std::endl;
			if (p.a.recreate_stage_files)
				std::cout << std::setw(PRINT_VAR_PAD) << std::left << "  Recreate mapped stage files: " << "yes" << std::endl;
			if (p.a.purge_stage_files)
				std::cout << std::setw(PRINT_VAR_PAD) << std::left << "  Purge mapped stage files: " << "yes" << std::endl;
		}		
		
		std::cout << std::setw(PRINT_VAR_PAD) << std::left << "  Reference frame: " << datum.GetName() << std::endl;
		std::cout << std::setw(PRINT_VAR_PAD) << std::left << "  Epoch: " << datum.GetEpoch_s() << std::endl;
		
		std::cout << std::setw(PRINT_VAR_PAD) << std::left << "  Geoid model: " << safe_absolute_path(p.n.ntv2_geoid_file) << std::endl;

		if (p.a.scale_normals_to_unity)
			std::cout << std::setw(PRINT_VAR_PAD) << std::left << "  Scale normals to unity: " << "yes" << std::endl;
		if (p.a.relaxation < 1.0)
			std::cout << std::setw(PRINT_VAR_PAD) << std::left << "  Initial relaxation: "
				<< std::fixed << std::setprecision(3) << p.a.relaxation << std::endl;
		if (!p.a.station_constraints.empty())
		{
			std::cout << std::setw(PRINT_VAR_PAD) << std::left << "  Station constraints: " << p.a.station_constraints << std::endl;
			if (p.i.apply_discontinuities)
				std::cout << std::setw(PRINT_VAR_PAD) << std::left << "  Apply discontinuities: " << "yes" << std::endl;
		}

		switch (p.a.adjust_mode)
		{
		case Phased_Block_1Mode:
		case PhasedMode:

			if (!std::filesystem::exists(p.a.seg_file))
			{
				std::cout << std::endl << std::endl << 
					"- Error: The required segmentation file does not exist:" << std::endl;  
				std::cout << "         " << p.a.seg_file << std::endl << std::endl;
				std::cout << "  Run  'segment " << p.g.network_name << "' to create a segmentation file" << std::endl << std::endl;
				cout_mutex.unlock();
				return EXIT_FAILURE;
			}

			// Load the segmentation file parameters into the dll, mainly
			// the segmentation block count, which is used by the progress 
			// thread before the adjustment begins
			netAdjust.LoadSegmentationFileParameters(p.a.seg_file);
		}
	
		std::cout << std::endl;
		std::cout << std::setw(PRINT_VAR_PAD) << std::left;
		switch (p.a.adjust_mode)
		{
		case SimultaneousMode:
			std::cout << "+ Simultaneous adjustment mode" << std::endl;
			break;
		case PhasedMode:
			std::cout << "+ Rigorous sequential phased adjustment mode";
			if (p.a.stage)
				std::cout << " (staged)";

			// If the user has not provided a seg file, check the meta of the default file
			if (!userSuppliedSegFile)
			{
				if (std::filesystem::last_write_time(p.a.seg_file) < std::filesystem::last_write_time(p.a.bst_file) ||
					std::filesystem::last_write_time(p.a.seg_file) < std::filesystem::last_write_time(p.a.bms_file))
				{
					// Has import been run after the segmentation file was created?
					if ((bst_meta_import && (std::filesystem::last_write_time(p.a.seg_file) < std::filesystem::last_write_time(p.a.bst_file))) || 
						(bms_meta_import && (std::filesystem::last_write_time(p.a.seg_file) < std::filesystem::last_write_time(p.a.bms_file))))
					{
						std::cout << std::endl << std::endl << 
							"- Error: The raw stations and measurements have been imported after" << std::endl <<
							"  the segmentation file was created:" << std::endl;

						time_t t_bst = file_time_to_time_t(std::filesystem::last_write_time(p.a.bst_file));
						time_t t_bms = file_time_to_time_t(std::filesystem::last_write_time(p.a.bms_file));
						time_t t_seg = file_time_to_time_t(std::filesystem::last_write_time(p.a.seg_file));

						std::cout << "   " << leafStr<std::string>(p.a.bst_file) << "  last modified on  " << ctime(&t_bst);
						std::cout << "   " << leafStr<std::string>(p.a.bms_file) << "  last modified on  " << ctime(&t_bms) << std::endl;
						std::cout << "   " << leafStr<std::string>(p.a.seg_file) << "  created on  " << ctime(&t_seg) << std::endl;
						std::cout << "  Run 'segment " << p.g.network_name << " [options]' to re-create the segmentation file, or re-run" << std::endl << 
							"  adjust using the --" << SEG_FILE << " option if the file " << 
							std::filesystem::path(p.a.seg_file).stem() << " must\n  be used." << std::endl << std::endl;
						cout_mutex.unlock();
						return EXIT_FAILURE;
					}
				}
			}

			// Has import been run after a staged adjustment was run, and this adjustment intends
			// to reuse memory mapped files created in the previous staged adjustment?
			if (p.a.stage && !p.a.recreate_stage_files)
			{
				// Simply test one file - the estimated stations file
				std::string est_mmapfile_name =
					p.g.output_folder + FOLDER_SLASH + 
					p.g.network_name + "-est.mtx";
				std::string est_mmapfile_wildcard =
					p.g.output_folder + FOLDER_SLASH + 
					p.g.network_name + "-*.mtx";
				if (std::filesystem::exists(est_mmapfile_name))
				{
					// Has import been run after the segmentation file was created?
					
					// This warning is here for the following reasons:
					//  1. import recreates the binary station and measurement files.  At this time, the reduced flag in the
					//     metadata record is set to false.
					//	2. When adjust is run for the first time, adjust reduces raw measurements to the ellipsoid (i.e. 
					//     applies n-values and deflections), applies scaling to GPS variance matrices, and then updates the 
					//     binary measurement file. adjust then sets the reduced flag in the metadata record to true.
					//  3. If a user decides to re-run import and segment, then attempts to run 
					//     'adjust network-name --stage' (i.e. without recreating the stage files) after import or segment, 
					//     then adjust will attempt to load memory map files using the same parameters from the first import
					//     and segment.
					//  Hence, force the user to run adjust with the --create-stage-files option.
					if ((bst_meta_import && (std::filesystem::last_write_time(est_mmapfile_name) < std::filesystem::last_write_time(p.a.bst_file))) ||
						(bms_meta_import && (std::filesystem::last_write_time(est_mmapfile_name) < std::filesystem::last_write_time(p.a.bms_file))))
					{
						std::cout << std::endl << std::endl << 
							"- Error: The raw stations and measurements have been imported after" << std::endl <<
							"  a staged adjustment created the memory map files:" << std::endl;
						
						time_t t_bst = file_time_to_time_t(std::filesystem::last_write_time(p.a.bst_file));
						time_t t_bms = file_time_to_time_t(std::filesystem::last_write_time(p.a.bms_file));
						time_t t_mtx = file_time_to_time_t(std::filesystem::last_write_time(est_mmapfile_name));

						std::cout << "   " << leafStr<std::string>(p.a.bst_file) << "  last modified on  " << ctime(&t_bst);
						std::cout << "   " << leafStr<std::string>(p.a.bms_file) << "  last modified on  " << ctime(&t_bms) << std::endl;
						std::cout << "   " << leafStr<std::string>(est_mmapfile_wildcard) << "  created on  " << ctime(&t_mtx) << std::endl;
						std::cout << "  To readjust this network, re-run adjust using the " << RECREATE_STAGE_FILES << " option." << std::endl;
						cout_mutex.unlock();
						return EXIT_FAILURE;
					}				
				}
			}
			
			if (p.a.multi_thread)
			{
				std::cout << std::endl << "+ Optimised for concurrent processing via multi-threading." << std::endl << std::endl;
				std::cout << "+ The active CPU supports the execution of " << std::thread::hardware_concurrency() << " concurrent threads.";
			}
			std::cout << std::endl;
			break;
		case Phased_Block_1Mode:
			std::cout << "+ Sequential phased adjustment resulting in rigorous estimates for Block 1 only" << std::endl;
			
			break;
		case SimulationMode:
			std::cout << "+ Adjustment simulation only" << std::endl;
			break;
		}
		std::cout << std::endl;
		
		if (p.a.report_mode)
		{
			std::cout << "+ Report last adjustment results" << std::endl;

			// Has report mode been requested as well as an argument to recreate stage files?
			// If so, return an error message and exit as this will lead to reporting of 
			// incorrect (zero!) results
			if (p.a.recreate_stage_files)
			{
				std::cout << std::endl <<
					"- Error: The option --" << RECREATE_STAGE_FILES << " cannot be used in Report results mode" << std::endl <<
					"  as it will erase the results from the latest adjustment and create new stage" << std::endl <<
					"  files initialised to zero." << std::endl << std::endl;
				cout_mutex.unlock();
				return EXIT_FAILURE;
			}
		}
		
		cout_mutex.unlock();
	}
	
	boost::posix_time::milliseconds elapsed_time(boost::posix_time::milliseconds(0));
	
	_ADJUST_STATUS_ adjustStatus;
	
	try {
		running = true;

        int nthreads_la = init_linear_algebra_threads();
        std::thread progress(dna_adjust_progress_thread(&netAdjust, &p));

        // Do adjustment using linear algebra threads
        dna_adjust_thread(&netAdjust, &p, &adjustStatus)(); 

        progress.join();
	
		switch (adjustStatus)
		{
		case ADJUST_EXCEPTION_RAISED:
			running = false;
			return EXIT_FAILURE;
		default:
			break;
		}

		if (p.a.report_mode)
			// Load variance matrices into memory
			DeserialiseVarianceMatrices(&netAdjust, &p);

		elapsed_time = boost::posix_time::milliseconds(netAdjust.adjustTime().count());

		// Print summary message
		PrintSummaryMessage(&netAdjust, &p, &elapsed_time);

		if (netAdjust.GetStatus() > ADJUST_THRESHOLD_EXCEEDED)
			return ADJUST_SUCCESS;

		// Generate statistics
		GenerateStatistics(&netAdjust, &p);

		if (p.a.max_iterations > 0)
			// Write variance matrices to disk
			SerialiseVarianceMatrices(&netAdjust, &p);

		// Print adjusted measurements to ADJ file
		PrintAdjustedMeasurements(&netAdjust, &p);

		// Print measurements to stations table
		PrintMeasurementstoStations(&netAdjust, &p);

		// Print adjusted stations to adj and xyz files
		PrintAdjustedNetworkStations(&netAdjust, &p);

		// close adj and xyz files
		netAdjust.CloseOutputFiles();

		// Print positional uncertainty
		PrintPositionalUncertainty(&netAdjust, &p);

		// Print station coordinates
		PrintStationCorrections(&netAdjust, &p);

		// Update bst and bms files with adjustment results
		UpdateBinaryFiles(&netAdjust, &p);

		// Print adjusted stations and measurements to DynaML
		ExportDynaML(&netAdjust, &p);

		// Print adjusted stations and measurements to DNA stn and msr
		ExportDNA(&netAdjust, &p);

		// Print adjusted stations and measurements to SINEX
		ExportSinex(&netAdjust, &p);
	}
	catch (const NetAdjustException& e) {
		cout_mutex.lock();
		std::cout << std::endl << 
			"- Error: " << e.what() << std::endl;
		cout_mutex.unlock();
		return EXIT_FAILURE;
	}

	if (!p.g.quiet)
		std::cout << std::endl << "+ Open " << leafStr<std::string>(p.o._adj_file) << " to view the adjustment details." << std::endl << std::endl;
	
	if (!userSuppliedSegFile)
		p.a.seg_file = "";
	if (!userSuppliedBstFile)
		p.a.bst_file = "";
	if (!userSuppliedBmsFile)
		p.a.bms_file = "";

	// Look for a project file.  If it exists, open and load it.
	// Update the import settings.
	// Print the project file. If it doesn't exist, it will be created.
	CDnaProjectFile projectFile;
	if (std::filesystem::exists(p.g.project_file))
		projectFile.LoadProjectFile(p.g.project_file);
	
	projectFile.UpdateSettingsAdjust(p);
	projectFile.UpdateSettingsOutput(p);
	projectFile.PrintProjectFile();	

	return ADJUST_SUCCESS;
}

