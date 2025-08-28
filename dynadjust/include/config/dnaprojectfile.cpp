//============================================================================
// Name         : dnaprojectfile.cpp
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
// Description  : DynAdjust Project class
//============================================================================

#include <include/config/dnatypes-gui.hpp>
#include <include/config/dnaprojectfile.hpp>
#include <include/functions/dnaiostreamfuncs.hpp>
#include <fstream>
#include <include/functions/dnatemplatefuncs.hpp>
#include <include/functions/dnastringfuncs.hpp>
#include <include/functions/dnastrutils.hpp>
#include <include/functions/dnastrmanipfuncs.hpp>
#include <include/functions/dnafilepathfuncs.hpp>
#include <include/config/dnaoptions-helper.hpp>

namespace dynadjust {

CDnaProjectFile::CDnaProjectFile(void)
{
	
}

// Create instance and initialise settings based on context
CDnaProjectFile::CDnaProjectFile(const std::string& projectFile, const UINT16& verifyContext)
{
	LoadProjectFile(projectFile);

	switch (verifyContext)
	{
	case importSetting:
		InitialiseImportSettings();
		break;
	case reftranSetting:
		InitialiseReftranSettings();
		break;
	case geoidSetting:
		InitialiseGeoidSettings();
		break;
	case segmentSetting:
		InitialiseSegmentSettings();
		break;
	case adjustSetting:
		InitialiseAdjustSettings();
		break;
	default:
		InitialiseGeneralSettings();
	}
}


//CDnaProjectFile::CDnaProjectFile(const project_settings& project)
//	: settings_(project)
//{
//}

//CDnaProjectFile::CDnaProjectFile(const CDnaProjectFile& newProject)
//	: settings_(newProject.settings_)
//{	
//}

//CDnaProjectFile& CDnaProjectFile::operator=(const CDnaProjectFile& rhs)
//{
//	if (this == &rhs)	// check for assignment to self!
//		return *this;
//
//	settings_ = rhs.settings_;
//
//	return *this;
//}

//bool CDnaProjectFile::operator==(const CDnaProjectFile& rhs) const
//{
//	return (settings_ == rhs.settings_);
//}
	

void CDnaProjectFile::LoadProjectFile(const std::string& projectFile)
{
	// load project file
	if (std::filesystem::exists(projectFile))
	{
		settings_.g.project_file = projectFile;
		LoadProjectFile();
	}
	else
	{
		std::stringstream err_msg;
		err_msg << "LoadProjectFile(): Project file " << 
		projectFile << " does not exist." << std::endl;	
		throw std::runtime_error(err_msg.str());
	}
}


void CDnaProjectFile::LoadProjectFile()
{
	// load project file using p.g.project_file
	std::ifstream dnaproj_file;

	std::stringstream err_msg;
	err_msg << "LoadProjectFile(): An error was encountered when opening " << 
		settings_.g.project_file << "." << std::endl;

	try {
		// create binary aml file.  Throws runtime_error on failure.
		file_opener(dnaproj_file, settings_.g.project_file,
			std::ios::in, ascii, true);
	}
	catch (const std::runtime_error& e) {
		err_msg << e.what();
		throw std::runtime_error(err_msg.str());
	}
	catch (...) {
		throw std::runtime_error(err_msg.str());
	}

	err_msg.str("");
	err_msg << "LoadProjectFile(): An error was encountered when loading from " << 
		settings_.g.project_file << "." << std::endl;

	std::string str, line, var, val;
	std::stringstream ss;

	settingMode mSetting(unknownSetting);
	
	while (!dnaproj_file.eof())
	{
		try {
			getline(dnaproj_file, line, '\n');
		}
		catch (std::ifstream::failure& f)
		{
			if (dnaproj_file.eof())
				break;
			
			err_msg << f.what() << std::endl;
			throw std::runtime_error(err_msg.str());
		}
		catch (...)
		{
			if (dnaproj_file.eof())
				break;

			err_msg << "Could not read file." << std::endl;
			throw std::runtime_error(err_msg.str());
		}


		// Blank line?
		if (line.empty())
			continue;

		// Line too short?
		if (line.length() <= PRINT_VAR_PAD)
			continue;

		// #variables
		if (line.find(section_variables) != std::string::npos)
		{
			mSetting = switchSetting;
			continue;
		}

		// #general
		if (line.find(section_general) != std::string::npos)
		{
			mSetting = generalSetting;
			continue;
		}

		// #import
		if (line.find(section_import) != std::string::npos)
		{
			mSetting = importSetting;
			continue;
		}

		// #reftran
		if (line.find(section_reftran) != std::string::npos)
		{
			mSetting = reftranSetting;
			continue;
		}

		// #general
		if (line.find(section_geoid) != std::string::npos)
		{
			mSetting = geoidSetting;
			continue;
		}

		// #segment
		if (line.find(section_segment) != std::string::npos)
		{
			mSetting = segmentSetting;
			continue;
		}

		// #adjust
		if (line.find(section_adjust) != std::string::npos)
		{
			mSetting = adjustSetting;
			continue;
		}

		// #output
		if (line.find(section_output) != std::string::npos)
		{
			mSetting = outputSetting;
			continue;
		}

		// #plot
		if (line.find(section_plot) != std::string::npos)
		{
			mSetting = plotSetting;
			continue;
		}

		// #display
		if (line.find(section_display) != std::string::npos)
		{
			mSetting = displaySetting;
			continue;
		}
		
		// Line?
		if (line.find("----------") != std::string::npos)
			continue;

		// Now get the variables and their value relating to the setting
		var = trimstr(line.substr(0, PRINT_VAR_PAD));
		if (line.length() > PRINT_VAR_PAD)
			val = trimstr(line.substr(PRINT_VAR_PAD));
		else
			val.clear();

		try {
			switch (mSetting)
			{
			case generalSetting:
				LoadSettingGeneral(mSetting, var, val);
				break;
			case importSetting:
				LoadSettingImport(mSetting, var, val);
				break;
			case reftranSetting:
				LoadSettingReftran(var, val);
				break;
			case geoidSetting:
				LoadSettingGeoid(var, val);
				break;
			case segmentSetting:
				LoadSettingSegment(var, val);
				break;
			case adjustSetting:
				LoadSettingAdjust(var, val);
				LoadSettingOutput(var, val);
				break;
			case outputSetting:
				LoadSettingOutput(var, val);
				break;
			case plotSetting:
				LoadSettingOutput(var, val);
				break;
			case displaySetting:
				break;
			default:
  				break;
			}
		}
		catch (const std::runtime_error& e) {
			err_msg << e.what();
			throw std::runtime_error(err_msg.str());
		}
		catch (...) {
			throw std::runtime_error(err_msg.str());
		}
	}

	dnaproj_file.clear();
	dnaproj_file.close();
}

void CDnaProjectFile::InitialiseGeneralSettings()
{
	if (settings_.g.network_name.empty())
	{
		std::stringstream ss;
		ss << "The project file does not provide the network name. " << std::endl << std::endl;
		throw std::runtime_error(ss.str());
	}
}

void CDnaProjectFile::InitialiseImportSettings()
{
	InitialiseGeneralSettings();

	if (settings_.i.input_files.empty() && 
		settings_.i.seg_file.empty())
	{
		std::stringstream ss;
		ss << "Nothing to do - the project file does not list any files to import. " << std::endl << std::endl;
		throw std::runtime_error(ss.str());
	}

	std::string firstPart(settings_.g.output_folder + FOLDER_SLASH + settings_.g.network_name + ".");
	
	if (!settings_.i.bst_file.empty())
	{
		if (!std::filesystem::exists(settings_.i.bst_file))
			// Does the file exist?  No.  
			// Assume it is a filename contained in the input folder.
			// import will throw an exception if it cannot be found.
			settings_.i.bst_file = formPath<std::string>(settings_.g.input_folder, 
				leafStr<std::string>(settings_.i.bst_file));
	}
	else
		settings_.i.bst_file = firstPart + "bst";	// binary stations file

	if (!settings_.i.bms_file.empty())
	{
		if (!std::filesystem::exists(settings_.i.bms_file))
			// Does the file exist?  No.  
			// Assume it is a filename contained in the input folder.
			// import will throw an exception if it cannot be found.
			settings_.i.bms_file = formPath<std::string>(settings_.g.input_folder, 
				leafStr<std::string>(settings_.i.bms_file));
	}
	else
		settings_.i.bms_file = firstPart + "bms";	// binary measurements file

	settings_.i.asl_file = firstPart + "asl";	// associated stations list
	settings_.i.aml_file = firstPart + "aml";	// associated measurements list
	settings_.i.map_file = firstPart + "map";	// station names map
	settings_.i.dst_file = firstPart + "dst";	// fuplicate stations
	settings_.i.dms_file = firstPart + "dms";	// duplicate measurements
	settings_.i.imp_file = firstPart + "imp";	// log

	if (settings_.o._msr_to_stn == 1)
		settings_.o._m2s_file = firstPart + "m2s";	// measurement to stations table

	if (!settings_.i.seg_file.empty())
	{
		if (!std::filesystem::exists(settings_.i.seg_file))
			// Does the file exist?  No.  
			// Assume it is a filename contained in the input folder.
			// import will throw an exception if it cannot be found.
			settings_.i.seg_file = formPath<std::string>(settings_.g.input_folder, 
				leafStr<std::string>(settings_.i.seg_file));
	}

	//////////////////////////////////////////////////////////////////////////////
	// GNSS scaling
	if (fabs(settings_.i.pscale - 1.0) > PRECISION_1E5 ||
		fabs(settings_.i.lscale - 1.0) > PRECISION_1E5 ||
		fabs(settings_.i.hscale - 1.0) > PRECISION_1E5 ||
		fabs(settings_.i.vscale - 1.0) > PRECISION_1E5 ||
		!settings_.i.scalar_file.empty())
		settings_.i.apply_scaling = 1;
	
	if (!settings_.i.scalar_file.empty())
	{
		if (!std::filesystem::exists(settings_.i.scalar_file))			
			// Does the file exist?  No.  
			// Assume it is a filename contained in the input folder.
			// import will throw an exception if it cannot be found.
			settings_.i.scalar_file = formPath<std::string>(settings_.g.input_folder, 
				leafStr<std::string>(settings_.i.scalar_file));	
	}

	// Create file name based on the provided block
	std::string fileName(settings_.g.network_name);
	if (settings_.i.import_block)
	{
		std::stringstream blk("");
		blk << ".block-" << settings_.i.import_block_number;
		fileName += blk.str();
	}

	if (settings_.i.export_dynaml)
	{
		if (settings_.i.export_single_xml_file)
			settings_.i.xml_outfile = formPath<std::string>(settings_.g.output_folder, 
				fileName, "xml");
		else
		{
			settings_.i.xml_stnfile = formPath<std::string>(settings_.g.output_folder, 
				fileName + "stn", "xml");
			settings_.i.xml_msrfile = formPath<std::string>(settings_.g.output_folder, 
				fileName + "msr", "xml");
		}
	}

	if (settings_.i.export_dna_files)
	{
		settings_.i.dna_stnfile = formPath<std::string>(settings_.g.output_folder, 
			fileName, "stn");
		settings_.i.dna_msrfile = formPath<std::string>(settings_.g.output_folder, 
			fileName, "msr");		
	}

	if (settings_.i.simulate_measurements)
		settings_.i.simulate_msrfile = formPath<std::string>(settings_.g.output_folder, 
			settings_.g.network_name, "simulated.msr");

}

void CDnaProjectFile::InitialiseReftranSettings()
{
	InitialiseGeneralSettings();

	if (settings_.r.reference_frame.empty())
	{
		std::stringstream ss;
		ss << "Nothing to do - the project file does not specify a reference frame." << std::endl << std::endl;
		throw std::runtime_error(ss.str());
	}

	std::string firstPart(settings_.g.output_folder + FOLDER_SLASH + settings_.g.network_name + ".");
	settings_.r.bst_file = firstPart + "bst";	// binary stations file
	settings_.r.bms_file = firstPart + "bms";	// binary measurements file
	settings_.r.rft_file = firstPart + "rft";	// log

	if (settings_.i.export_dynaml == 1)
	{
		if (settings_.i.export_single_xml_file == 1)
			settings_.i.xml_outfile = firstPart + "xml";		// single XML file
		else
		{
			settings_.i.xml_stnfile = formPath<std::string>(settings_.g.output_folder,
				settings_.g.network_name + "stn", "xml");
			settings_.i.xml_msrfile = formPath<std::string>(settings_.g.output_folder,
				settings_.g.network_name + "msr", "xml");;
		}
	}

	if (settings_.i.export_dna_files == 1)
	{
		settings_.i.dna_stnfile = firstPart + "stn";	// station file
		settings_.i.dna_msrfile = firstPart + "msr";	// measurement file
	}
	
	if (settings_.o._export_snx_file == 1)
		settings_.o._snx_file = firstPart + "snx";		// SINEX file

}

void CDnaProjectFile::InitialiseGeoidSettings()
{
	InitialiseGeneralSettings();

	// Expected behaviour: when running geoid using a project file, network name must exist
	if (settings_.g.network_name.empty())
	{
		std::stringstream ss;
		ss << "Nothing to do - the project file does not provide a network name. " << std::endl << std::endl;
		throw std::runtime_error(ss.str());
	}
	
	settings_.n.bst_file = formPath<std::string>(settings_.g.input_folder, settings_.g.network_name, "bst");		// binary stations file
	settings_.n.geo_file = formPath<std::string>(settings_.g.output_folder, settings_.g.network_name, "geo");	// dna geo file
	settings_.n.file_mode = 1;

	if (settings_.n.ntv2_geoid_file.empty())
	{
		std::stringstream ss;
		ss << "Nothing to do - the project file does not provide an NTv2 grid file path. " << std::endl << std::endl;
		throw std::runtime_error(ss.str());
	}
}

void CDnaProjectFile::InitialiseSegmentSettings()
{
	InitialiseGeneralSettings();

	std::string firstPart(settings_.g.output_folder + FOLDER_SLASH + settings_.g.network_name + ".");
	settings_.s.bst_file = firstPart + "bst";	// binary stations file
	settings_.s.bms_file = firstPart + "bms";	// binary measurements file
	settings_.s.asl_file = firstPart + "asl";	// associated stations list
	settings_.s.aml_file = firstPart + "aml";	// associated measurements list
	settings_.s.map_file = firstPart + "map";	// station names map
	settings_.s.seg_file = firstPart + "seg";	// segmentation file
	settings_.s.sap_file = firstPart + "sap";	// station appearance file
}

void CDnaProjectFile::InitialiseAdjustSettings()
{
	InitialiseGeneralSettings();

	std::string firstPart(settings_.g.output_folder + FOLDER_SLASH + settings_.g.network_name + ".");
	settings_.a.bst_file = firstPart + "bst";	// binary stations file
	settings_.a.bms_file = firstPart + "bms";	// binary measurements file

	if (!std::filesystem::exists(settings_.a.bst_file) || !std::filesystem::exists(settings_.a.bms_file))
	{
		std::stringstream ss;
		ss << "- Nothing to do: ";  
			
		if (settings_.g.network_name.empty())
			ss << std::endl << "network name has not been specified specified, and " << std::endl << "               ";  
		ss << settings_.a.bst_file << " and " << settings_.a.bms_file << " do not exist." << std::endl << std::endl;  
		throw std::runtime_error(ss.str());
	}
	
	// Station appearance file
	settings_.a.seg_file = firstPart + "seg";	// segmentation file
	settings_.s.sap_file = firstPart + "sap";	// station appearance file
	settings_.s.asl_file = firstPart + "asl";	// associated stations list
	settings_.s.aml_file = firstPart + "aml";	// associated measurements list
	settings_.a.map_file = firstPart + "map";	// station names map

	// Set up file names dependent on adjustment mode
	settings_.o._xyz_file = settings_.o._adj_file = 
		formPath<std::string>(settings_.g.output_folder, settings_.g.network_name);

	if (settings_.o._positional_uncertainty)
		settings_.o._apu_file = settings_.o._adj_file;

	if (settings_.o._init_stn_corrections)
		settings_.o._cor_file = settings_.o._adj_file;

	switch (settings_.a.adjust_mode)
	{
	case Phased_Block_1Mode:
	case PhasedMode:

		settings_.o._adj_file += ".phased";
		settings_.o._xyz_file += ".phased";

		if (settings_.o._positional_uncertainty)
			settings_.o._apu_file += ".phased";

		if (settings_.o._init_stn_corrections)
			settings_.o._cor_file += ".phased";
		
		if (settings_.a.adjust_mode == Phased_Block_1Mode)
		{
			settings_.o._adj_file += "-block1";
			settings_.o._xyz_file += "-block1";
			
			if (settings_.o._positional_uncertainty)
				settings_.o._apu_file += "-block1";

			if (settings_.o._init_stn_corrections)
				settings_.o._cor_file += "-block1";

		}
		else if (settings_.a.stage)
		{
			settings_.o._adj_file += "-stage";
			settings_.o._xyz_file += "-stage";

			if (settings_.o._positional_uncertainty)
				settings_.o._apu_file += "-stage";

			if (settings_.o._init_stn_corrections)
				settings_.o._cor_file += "-stage";

		}
		else if (settings_.a.multi_thread)
		{
			settings_.o._adj_file += "-mt";
			settings_.o._xyz_file += "-mt";
			
			if (settings_.o._positional_uncertainty)
				settings_.o._apu_file += "-mt";

			if (settings_.o._init_stn_corrections)
				settings_.o._cor_file += "-mt";
		}
		break;
	case SimultaneousMode:
		settings_.o._adj_file += ".simult";
		settings_.o._xyz_file += ".simult";
		
		if (settings_.o._positional_uncertainty)
			settings_.o._apu_file += ".simult";
		if (settings_.o._init_stn_corrections)
			settings_.o._cor_file += ".simult";
		break;
	}

	settings_.o._adj_file += ".adj";
	settings_.o._xyz_file += ".xyz";

	if (settings_.o._positional_uncertainty)
		settings_.o._apu_file += ".apu";

	if (settings_.o._init_stn_corrections)
		settings_.o._cor_file += ".cor";

	if (settings_.a.fixed_std_dev < 0.)
		settings_.a.fixed_std_dev = PRECISION_1E6;
	else if (settings_.a.fixed_std_dev < PRECISION_1E10)
		settings_.a.fixed_std_dev = PRECISION_1E6;

	if (settings_.a.free_std_dev < 0.)
		settings_.a.free_std_dev = 10.0;
	else if (settings_.a.free_std_dev < PRECISION_1E10)
		settings_.a.free_std_dev = 10.0;

}

void CDnaProjectFile::LoadSettingGeneral(const settingMode mSetting, const std::string& var, std::string& val)
{
	// general settings
	if (iequals(var, QUIET))
	{
		if (val.empty())
			return;			
		settings_.g.quiet = yesno_uint<UINT16, std::string>(val);
	}
	else if (iequals(var, VERBOSE))
	{
		if (val.empty())
			return;
		settings_.g.verbose = lexical_cast<UINT16, std::string>(val);
	}
	else if (iequals(var, PROJECT_FILE))
	{
		if (val.empty())
			return;			
		settings_.g.project_file = val;
	}
	else if (iequals(var, DYNADJUST_LOG_FILE))
	{
		if (val.empty())
			return;			
		settings_.g.log_file = val;
	}
	else if (iequals(var, NETWORK_NAME))
	{
		if (val.empty())
			return;
		string_string_pair t;
		switch (mSetting)
		{
		case switchSetting:
			t.first = var.c_str();
			t.second = val.substr(0,4).c_str();
			t.second = trimstr(t.second);
			if (t.second.substr(0,2) == "${" && t.second.at(3) == '}')
			{
				settings_.g.variables.push_back(t);
				std::sort(settings_.g.variables.begin(), settings_.g.variables.end(), PairCompareFirst<std::string, std::string>());
			}
			break;
		case generalSetting:
			settings_.g.network_name = val;
			break;
		default:
  			break;
		}			
	}
	else if (iequals(var, INPUT_FOLDER))
	{
		if (val.empty())
			return;			
		settings_.g.input_folder = val;
	}
	else if (iequals(var, OUTPUT_FOLDER))
	{
		if (val.empty())
			return;			
		settings_.g.output_folder = val;
	}
		
}
	
void CDnaProjectFile::LoadSettingImport(const settingMode mSetting, const std::string& var, std::string& val)
{
	// Import settings
	if (iequals(var, IMPORT_FILE))
	{
		if (val.empty())
			return;			
		settings_.i.input_files.push_back(formPath<std::string>(settings_.g.input_folder, val));
	}
	else if (iequals(var, REFERENCE_FRAME))
	{
		if (val.empty())
			return;
		settings_.i.reference_frame = val;
	}
	else if (iequals(var, EPOCH))
	{
		if (val.empty())
			return;
		settings_.i.epoch = val;
	}
	else if (iequals(var, OVERRIDE_INPUT_FRAME))
	{
		if (val.empty())
			return;
		settings_.i.override_input_rfame = yesno_uint<UINT16, std::string>(val);;
	}
	else if (iequals(var, IMPORT_GEO_FILE))
	{
		if (val.empty())
			return;
		settings_.i.geo_file = formPath<std::string>(settings_.g.input_folder, val);
		settings_.i.import_geo_file = 1;
	}
	else if (iequals(var, BIN_STN_FILE))
	{
		string_string_pair t;
		switch (mSetting)
		{
		case switchSetting:
			if (val.empty())
				return;
			t.first = var.c_str();
			t.second = val.substr(0,4).c_str();
			t.second = trimstr(t.second);
			if (t.second.substr(0,2) == "${" && t.second.at(3) == '}')
			{
				settings_.g.variables.push_back(t);
				std::sort(settings_.g.variables.begin(), settings_.g.variables.end(), PairCompareFirst<std::string, std::string>());
			}
			break;
		case importSetting:
			if (val.empty())
				// default export bst file - create in output folder
				settings_.i.bst_file = formPath<std::string>(settings_.g.output_folder, settings_.g.network_name, "bst");	// binary stations file
			else
				// user specified bst file - look in input folder
				settings_.i.bst_file = formPath<std::string>(settings_.g.input_folder, val);
			break;
		default:
  			break;
		}
	}
	else if (iequals(var, BIN_MSR_FILE))
	{
		string_string_pair t;
		switch (mSetting)
		{
		case switchSetting:
			if (val.empty())
				return;
			t.first = var.c_str();
			t.second = val.substr(0,4).c_str();
			t.second = trimstr(t.second);
			if (t.second.substr(0,2) == "${" && t.second.at(3) == '}')
			{
				settings_.g.variables.push_back(t);
				std::sort(settings_.g.variables.begin(), settings_.g.variables.end(), PairCompareFirst<std::string, std::string>());
			}
			break;
		case importSetting:
			if (val.empty())
				settings_.i.bms_file = formPath<std::string>(settings_.g.output_folder, settings_.g.network_name, "bms");	// binary stations file
			else
				settings_.i.bms_file = formPath<std::string>(settings_.g.input_folder, val);
			break;
		default:
  			break;
		}
	}
	else if (iequals(var, DNA_MAP_FILE))
	{
		if (val.empty())
			return;
		string_string_pair t;
		switch (mSetting)
		{
		case switchSetting:
			t.first = var.c_str();
			t.second = val.substr(0,4).c_str();
			t.second = trimstr(t.second);
			if (t.second.substr(0,2) == "${" && t.second.at(3) == '}')
			{
				settings_.g.variables.push_back(t);
				std::sort(settings_.g.variables.begin(), settings_.g.variables.end(), PairCompareFirst<std::string, std::string>());
			}
			break;
		case importSetting:
			settings_.i.map_file = formPath<std::string>(settings_.g.output_folder, val);
			break;
		default:
  			break;
		}			
	}
	else if (iequals(var, ASL_FILE))
	{
		if (val.empty())
			return;
		string_string_pair t;
		switch (mSetting)
		{
		case switchSetting:
			t.first = var.c_str();
			t.second = val.substr(0,4).c_str();
			t.second = trimstr(t.second);
			if (t.second.substr(0,2) == "${" && t.second.at(3) == '}')
			{
				settings_.g.variables.push_back(t);
				std::sort(settings_.g.variables.begin(), settings_.g.variables.end(), PairCompareFirst<std::string, std::string>());
			}
			break;
		case importSetting:
			settings_.i.asl_file = formPath<std::string>(settings_.g.output_folder, val);
			break;
		default:
  			break;
		}
	}
	else if (iequals(var, AML_FILE))
	{
		if (val.empty())
			return;
		string_string_pair t;
		switch (mSetting)
		{
		case switchSetting:
			t.first = var.c_str();
			t.second = val.substr(0,4).c_str();
			t.second = trimstr(t.second);
			if (t.second.substr(0,2) == "${" && t.second.at(3) == '}')
			{
				settings_.g.variables.push_back(t);
				std::sort(settings_.g.variables.begin(), settings_.g.variables.end(), PairCompareFirst<std::string, std::string>());
			}
			break;
		case importSetting:
			settings_.i.aml_file = formPath<std::string>(settings_.g.output_folder, val);
			break;
		default:
  			break;
		}
	}
	else if (iequals(var, DST_FILE))
	{
		if (val.empty())
			return;
		settings_.i.dst_file = formPath<std::string>(settings_.g.output_folder, val);
	}
	else if (iequals(var, DMS_FILE))
	{
		if (val.empty())
			return;
		settings_.i.dms_file = formPath<std::string>(settings_.g.output_folder, val);
	}
	else if (iequals(var, BOUNDING_BOX))
	{
		if (val.empty())
			return;
		settings_.i.bounding_box = val;
	}
	else if (iequals(var, GET_MSRS_TRANSCENDING_BOX))
	{
		if (val.empty())
			return;
		settings_.i.include_transcending_msrs = yesno_uint<UINT16, std::string>(val);
	}
	else if (iequals(var, INCLUDE_STN_ASSOC_MSRS))
	{
		if (val.empty())
			return;
		settings_.i.stn_associated_msr_include = valorno_string<std::string>(val);
	}
	else if (iequals(var, EXCLUDE_STN_ASSOC_MSRS))
	{
		if (val.empty())
			return;
		settings_.i.stn_associated_msr_exclude = valorno_string<std::string>(val);
	}
	else if (iequals(var, SPLIT_CLUSTERS))
	{
		if (val.empty())
			return;
		settings_.i.split_clusters = yesno_uint<UINT16, std::string>(val);
	}
	else if (iequals(var, IMPORT_SEG_BLOCK))
	{
		if (val.empty())
			return;
		settings_.i.import_block_number = valorno_uint<UINT16, std::string>(val, settings_.i.import_block);
	}
	else if (iequals(var, SEG_FILE))
	{
		if (val.empty())
			return;
		settings_.i.seg_file = formPath<std::string>(settings_.g.input_folder, val);
	}
	else if (iequals(var, PREFER_X_MSR_AS_G))
	{
		if (val.empty())
			return;
		settings_.i.prefer_single_x_as_g = yesno_uint<UINT16, std::string>(val);
	}
	else if (iequals(var, INCLUDE_MSRS))
	{
		if (val.empty())
			return;
		settings_.i.include_msrs = val;
	}
	else if (iequals(var, EXCLUDE_MSRS))
	{
		if (val.empty())
			return;
		settings_.i.exclude_msrs = val;
	}
	else if (iequals(var, STATION_RENAMING_FILE))
	{
		if (val.empty())
			return;
		settings_.i.stn_renamingfile = formPath<std::string>(settings_.g.input_folder, val);
	}
	else if (iequals(var, STATION_DISCONTINUITY_FILE))
	{
		if (val.empty())
			return;
		settings_.i.stn_discontinuityfile = formPath<std::string>(settings_.g.input_folder, val);
		settings_.i.apply_discontinuities = true;
	}
	else if (iequals(var, TEST_NEARBY_STNS))
	{
		if (val.empty())
			return;
		settings_.i.search_nearby_stn = yesno_uint<UINT16, std::string>(val);
	}
	else if (iequals(var, TEST_NEARBY_STN_DIST))
	{
		if (val.empty())
			return;
		settings_.i.search_stn_radius = DoubleFromString<double>(val);
	}
	else if (iequals(var, TEST_SIMILAR_MSRS))
	{
		if (val.empty())
			return;
		settings_.i.search_similar_msr = yesno_uint<UINT16, std::string>(val);
	}
	else if (iequals(var, TEST_SIMILAR_GNSS_MSRS))
	{
		if (val.empty())
			return;
		settings_.i.search_similar_msr_gx = yesno_uint<UINT16, std::string>(val);
	}
	else if (iequals(var, IGNORE_INSUFFICIENT_MSRS))
	{
		if (val.empty())
			return;
		settings_.i.ignore_insufficient_msrs = yesno_uint<UINT16, std::string>(val);
		}
	else if (iequals(var, IGNORE_SIMILAR_MSRS))
	{
		if (val.empty())
			return;
		settings_.i.ignore_similar_msr = yesno_uint<UINT16, std::string>(val);
	}
	else if (iequals(var, REMOVE_IGNORED_MSRS))
	{
		if (val.empty())
			return;
		settings_.i.remove_ignored_msr = yesno_uint<UINT16, std::string>(val);
	}
	else if (iequals(var, FLAG_UNUSED_STNS))
	{
		if (val.empty())
			return;
		settings_.i.flag_unused_stn = yesno_uint<UINT16, std::string>(val);
	}
	else if (iequals(var, TEST_INTEGRITY))
	{
		if (val.empty())
			return;			
		settings_.i.test_integrity = yesno_uint<UINT16, std::string>(val);
	}
	else if (iequals(var, VSCALE))
	{
		if (val.empty())
			return;
		settings_.i.vscale = DoubleFromString<double>(val);
	}
	else if (iequals(var, PSCALE))
	{
		if (val.empty())
			return;
		settings_.i.pscale = DoubleFromString<double>(val);
	}
	else if (iequals(var, LSCALE))
	{
		if (val.empty())
			return;
		settings_.i.lscale = DoubleFromString<double>(val);
	}
	else if (iequals(var, HSCALE))
	{
		if (val.empty())
			return;
		settings_.i.hscale = DoubleFromString<double>(val);
	}
	else if (iequals(var, SCALAR_FILE))
	{
		if (val.empty())
			return;
		settings_.i.scalar_file = formPath<std::string>(settings_.g.input_folder, val);
	}
	else if (iequals(var, EXPORT_XML_FILES))
	{
		if (val.empty())
			return;			
		settings_.i.export_dynaml = yesno_uint<UINT16, std::string>(val);
	}
	else if (iequals(var, EXPORT_SINGLE_XML_FILE))
	{
		if (val.empty())
			return;			
		settings_.i.export_single_xml_file = yesno_uint<UINT16, std::string>(val);
	}
	else if (iequals(var, EXPORT_FROM_BINARY))
	{
		if (val.empty())
			return;			
		settings_.i.export_from_bfiles = yesno_uint<UINT16, std::string>(val);
	}
	else if (iequals(var, EXPORT_DNA_FILES))
	{
		if (val.empty())
			return;			
		settings_.i.export_dna_files = yesno_uint<UINT16, std::string>(val);
	}
	else if (iequals(var, EXPORT_ASL_FILE))
	{
		if (val.empty())
			return;			
		settings_.i.export_asl_file = yesno_uint<UINT16, std::string>(val);
	}
	else if (iequals(var, EXPORT_AML_FILE))
	{
		if (val.empty())
			return;			
		settings_.i.export_aml_file = yesno_uint<UINT16, std::string>(val);
	}
	else if (iequals(var, EXPORT_MAP_FILE))
	{
		if (val.empty())
			return;			
		settings_.i.export_map_file = yesno_uint<UINT16, std::string>(val);
	}
	else if (iequals(var, SIMULATE_MSR_FILE))
	{
		if (val.empty())
			return;			
		settings_.i.simulate_measurements = yesno_uint<UINT16, std::string>(val);
	}
	//else if (iequals(var, VERIFY_COORDS))
	//{
	//	if (val.empty())
	//		return;			
	//	settings_.i.verify_coordinates = lexical_cast<UINT16, std::string>(val);
	//}
	
}
	
void CDnaProjectFile::LoadSettingReftran(const std::string& var, std::string& val)
{
	if (iequals(var, REFERENCE_FRAME))
	{
		if (val.empty())
			return;
		settings_.r.reference_frame = val;
	}
	else if (iequals(var, EPOCH))
	{
		if (val.empty())
			return;
		settings_.r.epoch = val;
	}
	else if (iequals(var, TECTONIC_PLATE_MODEL_OPTION))
	{
		if (val.empty())
			return;
		settings_.r.plate_model_option = lexical_cast<UINT16, std::string>(val);
	}
	else if (iequals(var, TECTONIC_PLATE_BDY_FILE))
	{
		if (val.empty())
			return;
		settings_.r.tpb_file = val;

		if (!std::filesystem::exists(val))
			if (std::filesystem::exists(formPath<std::string>(settings_.g.input_folder, val)))
				settings_.r.tpb_file = formPath<std::string>(settings_.g.input_folder, val);
	}
	else if (iequals(var, TECTONIC_PLATE_POLE_FILE))
	{
		if (val.empty())
			return;
		settings_.r.tpp_file = val;

		if (!std::filesystem::exists(val))
			if (std::filesystem::exists(formPath<std::string>(settings_.g.input_folder, val)))
				settings_.r.tpp_file = formPath<std::string>(settings_.g.input_folder, val);
	}
}
	
void CDnaProjectFile::LoadSettingGeoid(const std::string& var, std::string& val)
{
	if (iequals(var, DAT_FILEPATH))
	{
		if (val.empty())
			return;
		settings_.n.rdat_geoid_file = formPath<std::string>(settings_.g.input_folder, val);
	}
	else if (iequals(var, NTV2_FILEPATH))
	{
		if (val.empty())
			return;
		settings_.n.ntv2_geoid_file = val;

		if (!std::filesystem::exists(val))
			if (std::filesystem::exists(formPath<std::string>(settings_.g.input_folder, val)))
				settings_.n.ntv2_geoid_file = formPath<std::string>(settings_.g.input_folder, val);
	}
	else if (iequals(var, INPUT_FILE))
	{
		if (val.empty())
			return;
		settings_.n.input_file = formPath<std::string>(settings_.g.input_folder, val);
	}
	//else if (iequals(var, BIN_STN_FILE))
	//{
	//	if (val.empty())
	//		return;
	//	settings_.n.bst_file = formPath<std::string>(settings_.g.input_folder, val);
	//}
	else if (iequals(var, METHOD))
	{
		if (val.empty())
			return;
		settings_.n.interpolation_method = lexical_cast<UINT16, std::string>(val);
	}
	else if (iequals(var, DDEG_FORMAT))
	{
		if (val.empty())
			return;
		settings_.n.coordinate_format = lexical_cast<UINT16, std::string>(val);
	}
	else if (iequals(var, DIRECTION))
	{
		if (val.empty())
			return;
		settings_.n.ellipsoid_to_ortho = lexical_cast<UINT16, std::string>(val);
	}
	else if (iequals(var, CONVERT_BST_HT))
	{
		if (val.empty())
			return;
		settings_.n.convert_heights = yesno_uint<UINT16, std::string>(val);
	}
	else if (iequals(var, EXPORT_GEO_FILE))
	{
		if (val.empty())
			return;
		settings_.n.export_dna_geo_file = yesno_uint<UINT16, std::string>(val);
	}
}
	
void CDnaProjectFile::LoadSettingSegment(const std::string& var, std::string& val)
{
	if (iequals(var, NET_FILE))
	{
		if (val.empty())
			return;
		settings_.s.net_file = formPath<std::string>(settings_.g.input_folder, val);
	}
	else if (iequals(var, SEG_FILE))
	{
		if (val.empty())
			return;
		settings_.s.seg_file = val;

		if (!std::filesystem::exists(val))
			if (std::filesystem::exists(formPath<std::string>(settings_.g.input_folder, val)))
				settings_.s.seg_file = formPath<std::string>(settings_.g.input_folder, val);
	}
	else if (iequals(var, BIN_STN_FILE))
	{
		if (val.empty())
			return;
		settings_.s.bst_file = val;

		if (!std::filesystem::exists(val))
			if (std::filesystem::exists(formPath<std::string>(settings_.g.input_folder, val)))
				settings_.s.bst_file = formPath<std::string>(settings_.g.input_folder, val);
	}
	else if (iequals(var, BIN_MSR_FILE))
	{
		if (val.empty())
			return;
		settings_.s.bms_file = val;

		if (!std::filesystem::exists(val))
			if (std::filesystem::exists(formPath<std::string>(settings_.g.input_folder, val)))
				settings_.s.bms_file = formPath<std::string>(settings_.g.input_folder, val);
	}
	else if (iequals(var, SEG_MIN_INNER_STNS))
	{
		if (val.empty())
			return;
		settings_.s.min_inner_stations = lexical_cast<UINT32, std::string>(val);
	}
	else if (iequals(var, SEG_THRESHOLD_STNS))
	{
		if (val.empty())
			return;
		settings_.s.max_total_stations = lexical_cast<UINT32, std::string>(val);
	}
	else if (iequals(var, SEG_FORCE_CONTIGUOUS))
	{
		if (val.empty())
			return;
		settings_.s.force_contiguous_blocks = yesno_uint<UINT16, std::string>(val);
	}
	else if (iequals(var, SEG_STARTING_STN))
	{
		if (val.empty())
			return;
		settings_.s.seg_starting_stns = val;
	}
	else if (iequals(var, SEG_SEARCH_LEVEL))
	{
		if (val.empty())
			return;
		settings_.s.seg_search_level = lexical_cast<UINT16, std::string>(val);
	}
}
	
void CDnaProjectFile::LoadSettingAdjust(const std::string& var, std::string& val)
{
	if (iequals(var, BIN_STN_FILE))
	{
		if (val.empty())
			return;
		settings_.a.bst_file = val;

		if (!std::filesystem::exists(val))
			if (std::filesystem::exists(formPath<std::string>(settings_.g.input_folder, val)))
				settings_.a.bst_file = formPath<std::string>(settings_.g.input_folder, val);
	}
	else if (iequals(var, BIN_MSR_FILE))
	{
		if (val.empty())
			return;
		settings_.a.bms_file = val;

		if (!std::filesystem::exists(val))
			if (std::filesystem::exists(formPath<std::string>(settings_.g.input_folder, val)))
				settings_.a.bms_file = formPath<std::string>(settings_.g.input_folder, val);
	}
	else if (iequals(var, SEG_FILE))
	{
		if (val.empty())
			return;
		settings_.a.seg_file = val;

		if (!std::filesystem::exists(val))
			if (std::filesystem::exists(formPath<std::string>(settings_.g.input_folder, val)))
				settings_.a.seg_file = formPath<std::string>(settings_.g.input_folder, val);
	}
	else if (iequals(var, COMMENTS))
	{
		if (val.empty())
			return;
		settings_.a.comments = val;		
	}
	if (iequals(var, ADJUSTMENT_MODE))
	{
		if (val.empty())
			return;
		if (iequals(val, MODE_SIMULTANEOUS))
		{
			settings_.a.adjust_mode = SimultaneousMode;
			settings_.a.multi_thread = false;
			settings_.a.stage = false;
		}
		else if (iequals(val, MODE_PHASED))
		{
			settings_.a.adjust_mode = PhasedMode;
		}
		else if (iequals(val, MODE_PHASED_BLOCK1))
		{
			settings_.a.adjust_mode = Phased_Block_1Mode;
			settings_.a.multi_thread = false;
		}
		else if (iequals(val, MODE_SIMULATION))
		{
			settings_.a.adjust_mode = SimulationMode;
			settings_.a.multi_thread = false;
			settings_.a.stage = false;
		}
	}

	else if (iequals(var, MODE_PHASED_MT))
	{
		if (val.empty())
			return;
		settings_.a.multi_thread = yesno_uint<UINT16, std::string>(val);
	}
	else if (iequals(var, STAGED_ADJUSTMENT))
	{
		if (val.empty())
			return;
		settings_.a.stage = yesno_uint<UINT16, std::string>(val);
	}
	
	else if (iequals(var, CONF_INTERVAL))
	{
		if (val.empty())
			return;
		settings_.a.confidence_interval = FloatFromString<float>(val);
	}
	else if (iequals(var, ITERATION_THRESHOLD))
	{
		if (val.empty())
			return;
		settings_.a.iteration_threshold = FloatFromString<float>(val);
	}
	else if (iequals(var, MAX_ITERATIONS))
	{
		if (val.empty())
			return;
		settings_.a.max_iterations = lexical_cast<UINT16, std::string>(val);
	}
	else if (iequals(var, STN_CONSTRAINTS))
	{
		if (val.empty())
			return;
		settings_.a.station_constraints = val;
	}
	else if (iequals(var, FREE_STN_SD))
	{
		if (val.empty())
			return;
		settings_.a.free_std_dev = DoubleFromString<double>(val);
	}
	else if (iequals(var, FIXED_STN_SD))
	{
		if (val.empty())
			return;
		settings_.a.fixed_std_dev = DoubleFromString<double>(val);
	}
	//else if (iequals(var, LSQ_INVERSE_METHOD))
	//{
	//	if (val.empty())
	//		return;
	//	settings_.a.inverse_method_lsq = lexical_cast<UINT16, std::string>(val);
	//}
	else if (iequals(var, SCALE_NORMAL_UNITY))
	{
		if (val.empty())
			return;
		settings_.a.scale_normals_to_unity = yesno_uint<UINT16, std::string>(val);
	}
	else if (iequals(var, RECREATE_STAGE_FILES))
	{
		if (val.empty())
			return;
		settings_.a.recreate_stage_files = yesno_uint<UINT16, std::string>(val);
	}
	else if (iequals(var, PURGE_STAGE_FILES))
	{
		if (val.empty())
			return;
		settings_.a.purge_stage_files = yesno_uint<UINT16, std::string>(val) == 1;
	}
	else if (iequals(var, TYPE_B_GLOBAL))
	{
		if (val.empty())
			return;
		settings_.a.type_b_global = val;
	}
	else if (iequals(var, TYPE_B_FILE))
	{
		if (val.empty())
			return;
		settings_.a.type_b_file = val;

		if (!std::filesystem::exists(val))
			if (std::filesystem::exists(formPath<std::string>(settings_.g.input_folder, val)))
				settings_.a.type_b_file = formPath<std::string>(settings_.g.input_folder, val);
	}
}
	
void CDnaProjectFile::LoadSettingOutput(const std::string& var, std::string& val)
{
	if (iequals(var, OUTPUT_MSR_TO_STN))
	{
		if (val.empty())
			return;
		settings_.o._msr_to_stn = yesno_uint<UINT16, std::string>(val);
	}
	else if (iequals(var, OUTPUT_MSR_TO_STN_SORTBY))
	{
		if (val.empty())
			return;
		settings_.o._sort_msr_to_stn = lexical_cast<UINT16, std::string>(val);
	}
	else if (iequals(var, OUTPUT_ADJ_STN_ITER))
	{
		if (val.empty())
			return;
		settings_.o._adj_stn_iteration = yesno_uint<UINT16, std::string>(val);
	}
	else if (iequals(var, OUTPUT_ADJ_STAT_ITER))
	{
		if (val.empty())
			return;
		settings_.o._adj_stat_iteration = yesno_uint<UINT16, std::string>(val);
	}
	else if (iequals(var, OUTPUT_ADJ_MSR_ITER))
	{
		if (val.empty())
			return;
		settings_.o._adj_msr_iteration = yesno_uint<UINT16, std::string>(val);
	}
	else if (iequals(var, OUTPUT_CMP_MSR_ITER))
	{
		if (val.empty())
			return;
		settings_.o._cmp_msr_iteration = yesno_uint<UINT16, std::string>(val);
	}
	else if (iequals(var, OUTPUT_ADJ_MSR))
	{
		if (val.empty())
			return;
		settings_.o._adj_msr_final = yesno_uint<UINT16, std::string>(val);
	}
	else if (iequals(var, OUTPUT_ADJ_GNSS_UNITS))
	{
		if (val.empty())
			return;
		settings_.o._adj_gnss_units = lexical_cast<UINT16, std::string>(val);
	}
	else if (iequals(var, OUTPUT_ADJ_MSR_TSTAT))
	{
		if (val.empty())
			return;
		settings_.o._adj_msr_tstat = yesno_uint<UINT16, std::string>(val);
	}
	else if (iequals(var, OUTPUT_ADJ_MSR_SORTBY))
	{
		if (val.empty())
			return;
		settings_.o._sort_adj_msr = lexical_cast<UINT16, std::string>(val);
	}
	else if (iequals(var, OUTPUT_ADJ_MSR_DBID))
	{
		if (val.empty())
			return;
		settings_.o._database_ids = yesno_uint<UINT16, std::string>(val);
	}
	else if (iequals(var, OUTPUT_ADJ_STN_BLOCKS))
	{
		if (val.empty())
			return;
		settings_.o._output_stn_blocks = yesno_uint<UINT16, std::string>(val);
	}
	else if (iequals(var, OUTPUT_ADJ_MSR_BLOCKS))
	{
		if (val.empty())
			return;
		settings_.o._output_msr_blocks = yesno_uint<UINT16, std::string>(val);
	}
	else if (iequals(var, OUTPUT_ADJ_STN_SORT_ORDER))
	{
		if (val.empty())
			return;
		settings_.o._sort_stn_file_order = yesno_uint<UINT16, std::string>(val);
	}
	else if (iequals(var, OUTPUT_STN_COORD_TYPES))
	{
		if (val.empty())
			return;
		settings_.o._stn_coord_types = val;
	}
	else if (iequals(var, OUTPUT_ANGULAR_TYPE_STN))
	{
		if (val.empty())
			return;
		settings_.o._angular_type_stn = lexical_cast<UINT16, std::string>(val);
	}
	else if (iequals(var, OUTPUT_STN_CORR))
	{
		if (val.empty())
			return;
		settings_.o._stn_corr = yesno_uint<UINT16, std::string>(val);
	}
	else if (iequals(var, OUTPUT_PRECISION_METRES_STN))
	{
		if (val.empty())
			return;
		settings_.o._precision_metres_stn = lexical_cast<UINT16, std::string>(val);
	}
	else if (iequals(var, OUTPUT_PRECISION_SECONDS_STN))
	{
		if (val.empty())
			return;
		settings_.o._precision_seconds_stn = lexical_cast<UINT16, std::string>(val);
	}
	else if (iequals(var, OUTPUT_PRECISION_METRES_MSR))
	{
		if (val.empty())
			return;
		settings_.o._precision_metres_msr = lexical_cast<UINT16, std::string>(val);
	}
	else if (iequals(var, OUTPUT_PRECISION_SECONDS_MSR))
	{
		if (val.empty())
			return;
		settings_.o._precision_seconds_msr = lexical_cast<UINT16, std::string>(val);
	}
	else if (iequals(var, OUTPUT_ANGULAR_TYPE_MSR))
	{
		if (val.empty())
			return;
		settings_.o._angular_type_msr = lexical_cast<UINT16, std::string>(val);
	}
	else if (iequals(var, OUTPUT_DMS_FORMAT_MSR))
	{
		if (val.empty())
			return;
		settings_.o._dms_format_msr = lexical_cast<UINT16, std::string>(val);
	}
	else if (iequals(var, OUTPUT_POS_UNCERTAINTY))
	{
		if (val.empty())
			return;
		settings_.o._positional_uncertainty = yesno_uint<UINT16, std::string>(val);
	}
	else if (iequals(var, OUTPUT_APU_CORRELATIONS))
	{
		if (val.empty())
			return;
		settings_.o._output_pu_covariances = yesno_uint<UINT16, std::string>(val);
	}
	else if (iequals(var, OUTPUT_APU_UNITS))
	{
		if (val.empty())
			return;
		settings_.o._apu_vcv_units = lexical_cast<UINT16, std::string>(val);
	}
	else if (iequals(var, OUTPUT_STN_COR_FILE))
	{
		if (val.empty())
			return;
		settings_.o._init_stn_corrections = yesno_uint<UINT16, std::string>(val);
	}
	else if (iequals(var, HZ_CORR_THRESHOLD))
	{
		if (val.empty())
			return;
		settings_.o._hz_corr_threshold = DoubleFromString<double>(val);
	}
	else if (iequals(var, VT_CORR_THRESHOLD))
	{
		if (val.empty())
			return;
		settings_.o._vt_corr_threshold = DoubleFromString<double>(val);
	}
	//else if (iequals(var, UPDATE_ORIGINAL_STN_FILE))
	//{
	//	if (val.empty())
	//		return;
	//	settings_.o. = yesno_uint<UINT16, std::string>(val);
	//}
	else if (iequals(var, EXPORT_XML_STN_FILE))
	{
		if (val.empty())
			return;
		settings_.o._export_xml_stn_file= yesno_uint<UINT16, std::string>(val);
	}
	else if (iequals(var, EXPORT_DNA_STN_FILE))
	{
		if (val.empty())
			return;
		settings_.o._export_dna_stn_file = yesno_uint<UINT16, std::string>(val);
	}
	else if (iequals(var, EXPORT_SNX_FILE))
	{
		if (val.empty())
			return;
		settings_.o._export_snx_file = yesno_uint<UINT16, std::string>(val);
	}
	// ???
	// EXPORT_DISCONT_FILE = "export-discont-file";
	
}
	
//void CDnaProjectFile::LoadSettingPlot(const std::string& var, std::string& val)
//{
//}
//	
//void CDnaProjectFile::LoadSettingDisplay(const std::string& var, std::string& val)
//{
//}

template <typename T>
void CDnaProjectFile::AddOptionValue(std::ostream& os, const char* const option, const T& value) {
	os << " --" << option << " " << value;
}

template <typename T>
void CDnaProjectFile::AddDefaultValue(std::ostream& os, const T& value)
{
	os << " " << value;
}

template <typename T, typename U>
void CDnaProjectFile::PrintRecord(std::ostream& os, const T& variable, const U& value) {

	PrintVariable(os, variable);
	PrintValue(os, value);
}
	

template <typename T>
void CDnaProjectFile::PrintVariable(std::ostream& os, const T& variable) {
	// width of 35
	os << std::setw(PRINT_VAR_PAD) << std::left << variable;
}
	

template <typename T>
void CDnaProjectFile::PrintValue(std::ostream& os, const T& value) {
	// width of 45
	os << std::setw(PRINT_VAL_PAD) << std::left << value << std::endl;
}
	

std::string CDnaProjectFile::FormCommandLineOptionsStringImport()
{
	std::stringstream options;
	//////////////////////////////////////////////////////
	// General settings
	AddDefaultValue(options, FormCommandLineOptionsStringGeneral());
	
	//////////////////////////////////////////////////////
	// Import settings
	// files
	for (_it_vstr_const _it_tmp(settings_.i.input_files.begin());
		_it_tmp!=settings_.i.input_files.end(); ++_it_tmp)
		AddDefaultValue<std::string>(options, 
			leafStr<std::string>(trimstr<std::string>(_it_tmp->c_str())));
	

	return options.str();
}
	
std::string CDnaProjectFile::FormCommandLineOptionsStringGeneral()
{
	std::stringstream options;
	//////////////////////////////////////////////////////
	// General settings
	// project file
	AddOptionValue<std::string>(options, PROJECT_FILE, settings_.g.project_file);
	// network name
	AddOptionValue<std::string>(options, NETWORK_NAME, settings_.g.network_name);
	// verbose
	AddOptionValue<UINT16>(options, VERBOSE, settings_.g.verbose);
	// verbose
	AddOptionValue<UINT16>(options, QUIET, settings_.g.quiet);

	return options.str();
}
	
void CDnaProjectFile::PrintProjectFile()
{
	std::ofstream dnaproj_file;

	std::stringstream err_msg;
	err_msg << "PrintProjectFile(): An error was encountered when opening " << 
		settings_.g.project_file << "." << std::endl;

	try {
		// create binary aml file.  Throws runtime_error on failure.
		file_opener(dnaproj_file, settings_.g.project_file,
			std::ios::out, ascii);
	}
	catch (const std::runtime_error& e) {
		err_msg << e.what();
		throw std::runtime_error(err_msg.str());
	}
	catch (...) {
		throw std::runtime_error(err_msg.str());
	}

	err_msg.str("");
	err_msg << "PrintProjectFile(): An error was encountered when writing to " << 
		settings_.g.project_file << "." << std::endl;

	std::stringstream ss;

	// Write header line
	dnaproj_header(dnaproj_file, settings_.g.network_name + " project file");
	dnaproj_file << std::endl << std::endl;

	////////////////////////////////////////////////////////////////////////////////////////////////
	// #general
	ss.str("");
	ss << section_general << " (" << PRINT_VAR_PAD << ")";
	PrintRecord(dnaproj_file, ss.str(), "VALUE");
	dnaproj_file << OUTPUTLINE << std::endl;
	
	PrintRecord(dnaproj_file, NETWORK_NAME, settings_.g.network_name);								// network name
	PrintRecord(dnaproj_file, INPUT_FOLDER, safe_absolute_path(settings_.g.input_folder));	// Path containing all input files
	PrintRecord(dnaproj_file, OUTPUT_FOLDER, safe_absolute_path(settings_.g.output_folder));	// Path for all output files
	PrintRecord(dnaproj_file, VERBOSE, settings_.g.verbose);										// Give detailed information about what dnainterop is doing.
																									// 0: No information (default)
																									// 1: Helpful information
																									// 2: Extended information\n3: Debug level information
	PrintRecord(dnaproj_file, QUIET, yesno_string(settings_.g.quiet));								// Run quietly?
	PrintRecord(dnaproj_file, PROJECT_FILE, safe_absolute_path(settings_.g.project_file));	// project file
	PrintRecord(dnaproj_file, DYNADJUST_LOG_FILE, safe_absolute_path(settings_.g.log_file));	// dynadjust log file
	
	dnaproj_file << std::endl;

	////////////////////////////////////////////////////////////////////////////////////////////////
	// #import
	ss.str("");
	ss << section_import << " (" << PRINT_VAR_PAD << ")";
	PrintRecord(dnaproj_file, ss.str(), "VALUE");
	dnaproj_file << OUTPUTLINE << std::endl;
	
	// Input files
	for_each (settings_.i.input_files.begin(), settings_.i.input_files.end(),
		[&dnaproj_file, this] (std::string file) {
			PrintRecord(dnaproj_file, IMPORT_FILE, leafStr<std::string>(file.c_str()));
	});

	// geoid file
	PrintRecord(dnaproj_file, IMPORT_GEO_FILE, 
		(settings_.i.geo_file.empty() ? " " : leafStr<std::string>(settings_.i.geo_file)));

	// reference frame settings
	PrintRecord(dnaproj_file, REFERENCE_FRAME, settings_.i.reference_frame);
	PrintRecord(dnaproj_file, EPOCH, settings_.i.epoch);
	PrintRecord(dnaproj_file, OVERRIDE_INPUT_FRAME, 
		yesno_string(settings_.i.override_input_rfame));

	// data screening
	PrintRecord(dnaproj_file, BOUNDING_BOX, settings_.i.bounding_box);
	PrintRecord(dnaproj_file, GET_MSRS_TRANSCENDING_BOX, 
		yesno_string(settings_.i.include_transcending_msrs));
	
	PrintRecord(dnaproj_file, INCLUDE_STN_ASSOC_MSRS,
		(settings_.i.stn_associated_msr_include.empty() ? "no" : settings_.i.stn_associated_msr_include));
	PrintRecord(dnaproj_file, EXCLUDE_STN_ASSOC_MSRS,
		(settings_.i.stn_associated_msr_exclude.empty() ? "no" : settings_.i.stn_associated_msr_exclude));
	PrintRecord(dnaproj_file, SPLIT_CLUSTERS,
		yesno_string(settings_.i.split_clusters));	
	PrintRecord(dnaproj_file, IMPORT_SEG_BLOCK,
		(settings_.i.import_block ? val_uint<std::string, UINT32>(settings_.i.import_block_number) : "no"));
	PrintRecord(dnaproj_file, SEG_FILE, 
		(settings_.i.seg_file.empty() ? " " : leafStr<std::string>(settings_.i.seg_file)));
	PrintRecord(dnaproj_file, PREFER_X_MSR_AS_G, 
		yesno_string(settings_.i.prefer_single_x_as_g));	
	
	PrintRecord(dnaproj_file, INCLUDE_MSRS, settings_.i.include_msrs);
	PrintRecord(dnaproj_file, EXCLUDE_MSRS, settings_.i.exclude_msrs);							
	PrintRecord(dnaproj_file, STATION_RENAMING_FILE,
		leafStr<std::string>(settings_.i.stn_renamingfile));						
	PrintRecord(dnaproj_file, STATION_DISCONTINUITY_FILE,
		leafStr<std::string>(settings_.i.stn_discontinuityfile));						
	PrintRecord(dnaproj_file, TEST_NEARBY_STNS,
		yesno_string(settings_.i.search_nearby_stn));						
	PrintRecord(dnaproj_file, TEST_NEARBY_STN_DIST, settings_.i.search_stn_radius);					
	PrintRecord(dnaproj_file, TEST_SIMILAR_GNSS_MSRS,
		yesno_string(settings_.i.search_similar_msr_gx));						
	PrintRecord(dnaproj_file, TEST_SIMILAR_MSRS,
		yesno_string(settings_.i.search_similar_msr));
	PrintRecord(dnaproj_file, IGNORE_INSUFFICIENT_MSRS,
		yesno_string(settings_.i.ignore_insufficient_msrs));
	PrintRecord(dnaproj_file, IGNORE_SIMILAR_MSRS,
		yesno_string(settings_.i.ignore_similar_msr));					
	PrintRecord(dnaproj_file, REMOVE_IGNORED_MSRS,
		yesno_string(settings_.i.remove_ignored_msr));					
	PrintRecord(dnaproj_file, FLAG_UNUSED_STNS,
		yesno_string(settings_.i.flag_unused_stn));						
	PrintRecord(dnaproj_file, TEST_INTEGRITY,
		yesno_string(settings_.i.test_integrity));						
	
	// GNSS variance matrix scaling
	PrintRecord(dnaproj_file, VSCALE, settings_.i.vscale);				
	PrintRecord(dnaproj_file, PSCALE, settings_.i.pscale);				
	PrintRecord(dnaproj_file, LSCALE, settings_.i.lscale);				
	PrintRecord(dnaproj_file, HSCALE, settings_.i.hscale);				
	PrintRecord(dnaproj_file, SCALAR_FILE, leafStr<std::string>(settings_.i.scalar_file));

	// export options
	PrintRecord(dnaproj_file, EXPORT_XML_FILES, 
		yesno_string(settings_.i.export_dynaml));						// Create DynaML output file
	PrintRecord(dnaproj_file, EXPORT_SINGLE_XML_FILE, 
		yesno_string(settings_.i.export_single_xml_file));				// Create single station and measurement DynaML output files
	PrintRecord(dnaproj_file, EXPORT_DNA_FILES, 
		yesno_string(settings_.i.export_dna_files));
	PrintRecord(dnaproj_file, EXPORT_ASL_FILE,
		yesno_string(settings_.i.export_asl_file));
	PrintRecord(dnaproj_file, EXPORT_AML_FILE, 
		yesno_string(settings_.i.export_aml_file));
	PrintRecord(dnaproj_file, EXPORT_MAP_FILE,
		yesno_string(settings_.i.export_map_file));
	PrintRecord(dnaproj_file, SIMULATE_MSR_FILE, 
		yesno_string(settings_.i.simulate_measurements));

	dnaproj_file << std::endl;
	
	////////////////////////////////////////////////////////////////////////////////////////////////
	// #reftran
	ss.str("");
	ss << section_reftran << " (" << PRINT_VAR_PAD << ")";
	PrintRecord(dnaproj_file, ss.str(), "VALUE");
	dnaproj_file << OUTPUTLINE << std::endl;
	
	PrintRecord(dnaproj_file, REFERENCE_FRAME, settings_.r.reference_frame);
	PrintRecord(dnaproj_file, EPOCH, settings_.r.epoch);

	PrintRecord(dnaproj_file, TECTONIC_PLATE_MODEL_OPTION, settings_.r.plate_model_option);					// Plate motion model option
	PrintRecord(dnaproj_file, TECTONIC_PLATE_BDY_FILE, leafStr<std::string>(settings_.r.tpb_file));				// Tectonic plate boundary file
	PrintRecord(dnaproj_file, TECTONIC_PLATE_POLE_FILE, leafStr<std::string>(settings_.r.tpp_file));				// Tectonic plate pole file
	
	dnaproj_file << std::endl;
	
	////////////////////////////////////////////////////////////////////////////////////////////////
	// #geoid
	ss.str("");
	ss << section_geoid << " (" << PRINT_VAR_PAD << ")";
	PrintRecord(dnaproj_file, ss.str(), "VALUE");
	dnaproj_file << OUTPUTLINE << std::endl;

	// Output configured to populate binary station files
	PrintRecord(dnaproj_file, NTV2_FILEPATH, safe_absolute_path(settings_.n.ntv2_geoid_file));					// Full file path to geoid file
	
	PrintRecord(dnaproj_file, METHOD, settings_.n.interpolation_method);
	PrintRecord(dnaproj_file, DDEG_FORMAT, settings_.n.coordinate_format);
	PrintRecord(dnaproj_file, DIRECTION, settings_.n.ellipsoid_to_ortho);
	//PrintRecord(dnaproj_file, CONVERT_BST_HT,
	//	yesno_string(settings_.n.convert_heights));
	PrintRecord(dnaproj_file, EXPORT_GEO_FILE,
		yesno_string(settings_.n.export_dna_geo_file));
	
	dnaproj_file << std::endl;

	////////////////////////////////////////////////////////////////////////////////////////////////
	// #segment
	ss.str("");
	ss << section_segment << " (" << PRINT_VAR_PAD << ")";
	PrintRecord(dnaproj_file, ss.str(), "VALUE");
	dnaproj_file << OUTPUTLINE << std::endl;
	
	PrintRecord(dnaproj_file, NET_FILE, leafStr<std::string>(settings_.s.net_file));				// Starting stations file
	PrintRecord(dnaproj_file, SEG_FILE, leafStr<std::string>(settings_.s.seg_file));				// Segmentation output file
	PrintRecord(dnaproj_file, SEG_MIN_INNER_STNS, settings_.s.min_inner_stations);			// Minimum number of inner stations per block
	PrintRecord(dnaproj_file, SEG_THRESHOLD_STNS, settings_.s.max_total_stations);			// Maximum number of total stations per block
	PrintRecord(dnaproj_file, SEG_FORCE_CONTIGUOUS,
		yesno_string(settings_.s.force_contiguous_blocks));

	// Stations to be incorporated within the first block.
	PrintRecord(dnaproj_file, SEG_STARTING_STN, settings_.s.seg_starting_stns);	

	dnaproj_file << std::endl;

	////////////////////////////////////////////////////////////////////////////////////////////////
	// #adjust
	ss.str("");
	ss << section_adjust << " (" << PRINT_VAR_PAD << ")";
	PrintRecord(dnaproj_file, ss.str(), "VALUE");
	dnaproj_file << OUTPUTLINE << std::endl;

	PrintRecord(dnaproj_file, SEG_FILE, leafStr<std::string>(settings_.a.seg_file));				// Starting stations file
	PrintRecord(dnaproj_file, COMMENTS, settings_.a.comments);								// Starting stations file
	
	PrintRecord(dnaproj_file, ADJUSTMENT_MODE, adjustmentMode<std::string, UINT32>(settings_.a.adjust_mode));

	PrintRecord(dnaproj_file, MODE_PHASED_MT, 
		yesno_string(settings_.a.multi_thread));
	PrintRecord(dnaproj_file, STAGED_ADJUSTMENT, 
		yesno_string(settings_.a.stage));
	
	PrintRecord(dnaproj_file, CONF_INTERVAL, settings_.a.confidence_interval);				// Confidence interval
	PrintRecord(dnaproj_file, ITERATION_THRESHOLD, settings_.a.iteration_threshold);		// Threshold to halt iterations
	PrintRecord(dnaproj_file, MAX_ITERATIONS, settings_.a.max_iterations);					// Number of iteration before a solution is adopted
	PrintRecord(dnaproj_file, STN_CONSTRAINTS, settings_.a.station_constraints);			
	
	ss.str("");
	ss << std::fixed << std::setprecision(3) << settings_.a.free_std_dev;
	PrintRecord(dnaproj_file, FREE_STN_SD, ss.str());										// SD for free stations
	ss.str("");
	ss << std::scientific << std::setprecision(4) << settings_.a.fixed_std_dev;
	PrintRecord(dnaproj_file, FIXED_STN_SD, ss.str());										// SD for fixed stations
	
	//PrintRecord(dnaproj_file, LSQ_INVERSE_METHOD, settings_.a.inverse_method_lsq);			// Least squares inverse method

	PrintRecord(dnaproj_file, SCALE_NORMAL_UNITY, 
		yesno_string(settings_.a.scale_normals_to_unity));									// Scale normals to unity before inversion
	PrintRecord(dnaproj_file, RECREATE_STAGE_FILES, 
		yesno_string(settings_.a.recreate_stage_files));									// Recreate stage files
	PrintRecord(dnaproj_file, PURGE_STAGE_FILES, 
		yesno_string(settings_.a.purge_stage_files));										// Purge stage files

	PrintRecord(dnaproj_file, TYPE_B_GLOBAL, settings_.a.type_b_global);					// Global Type B uncertainties
	PrintRecord(dnaproj_file, TYPE_B_FILE, leafStr<std::string>(settings_.a.type_b_file));		// Type B uncertainty file
	
	dnaproj_file << std::endl;

	////////////////////////////////////////////////////////////////////////////////////////////////
	// #output
	ss.str("");
	ss << section_output << " (" << PRINT_VAR_PAD << ")";
	PrintRecord(dnaproj_file, ss.str(), "VALUE");
	dnaproj_file << OUTPUTLINE << std::endl;

	PrintRecord(dnaproj_file, OUTPUT_MSR_TO_STN, 
		yesno_string(settings_.o._msr_to_stn));
	PrintRecord(dnaproj_file, OUTPUT_MSR_TO_STN_SORTBY, settings_.o._sort_msr_to_stn);
	PrintRecord(dnaproj_file, OUTPUT_ADJ_STN_ITER, 
		yesno_string(settings_.o._adj_stn_iteration));
	PrintRecord(dnaproj_file, OUTPUT_ADJ_STAT_ITER, 
		yesno_string(settings_.o._adj_stat_iteration));
	PrintRecord(dnaproj_file, OUTPUT_ADJ_MSR_ITER, 
		yesno_string(settings_.o._adj_msr_iteration));
	PrintRecord(dnaproj_file, OUTPUT_CMP_MSR_ITER, 
		yesno_string(settings_.o._cmp_msr_iteration));
	PrintRecord(dnaproj_file, OUTPUT_ADJ_MSR, 
		yesno_string(settings_.o._adj_msr_final));
	PrintRecord(dnaproj_file, OUTPUT_ADJ_GNSS_UNITS, 
		settings_.o._adj_gnss_units);	
	PrintRecord(dnaproj_file, OUTPUT_ADJ_MSR_TSTAT, 
		yesno_string(settings_.o._adj_msr_tstat));

	PrintRecord(dnaproj_file, OUTPUT_ADJ_MSR_SORTBY, settings_.o._sort_adj_msr);

	PrintRecord(dnaproj_file, OUTPUT_ADJ_MSR_DBID, 
		yesno_string(settings_.o._database_ids));			
	
	PrintRecord(dnaproj_file, OUTPUT_ADJ_MSR_BLOCKS, 
		yesno_string(settings_.o._output_msr_blocks));				
	PrintRecord(dnaproj_file, OUTPUT_ADJ_STN_SORT_ORDER, 
		yesno_string(settings_.o._sort_stn_file_order));			
	PrintRecord(dnaproj_file, OUTPUT_STN_COORD_TYPES, settings_.o._stn_coord_types);
	PrintRecord(dnaproj_file, OUTPUT_ANGULAR_TYPE_STN, settings_.o._angular_type_stn);
	PrintRecord(dnaproj_file, OUTPUT_STN_CORR, 
		yesno_string(settings_.o._stn_corr));
	PrintRecord(dnaproj_file, OUTPUT_PRECISION_METRES_STN, settings_.o._precision_metres_stn);	
	PrintRecord(dnaproj_file, OUTPUT_PRECISION_SECONDS_STN, settings_.o._precision_seconds_stn);
	PrintRecord(dnaproj_file, OUTPUT_PRECISION_METRES_MSR, settings_.o._precision_metres_msr);	
	PrintRecord(dnaproj_file, OUTPUT_PRECISION_SECONDS_MSR, settings_.o._precision_seconds_msr);
	PrintRecord(dnaproj_file, OUTPUT_ANGULAR_TYPE_MSR, settings_.o._angular_type_msr);
	PrintRecord(dnaproj_file, OUTPUT_DMS_FORMAT_MSR, settings_.o._dms_format_msr);				
	PrintRecord(dnaproj_file, OUTPUT_POS_UNCERTAINTY, 
		yesno_string(settings_.o._positional_uncertainty));
	PrintRecord(dnaproj_file, OUTPUT_APU_CORRELATIONS, 
		yesno_string(settings_.o._output_pu_covariances));		
	
	PrintRecord(dnaproj_file, OUTPUT_APU_UNITS, 
		settings_.o._apu_vcv_units);

	PrintRecord(dnaproj_file, OUTPUT_STN_COR_FILE, 
		yesno_string(settings_.o._init_stn_corrections));
	
	ss.str("");
	ss << std::fixed << std::setprecision(3) << settings_.o._hz_corr_threshold;
	PrintRecord(dnaproj_file, HZ_CORR_THRESHOLD, ss.str());				
	ss.str("");
	ss << std::fixed << std::setprecision(3) << settings_.o._vt_corr_threshold;
	PrintRecord(dnaproj_file, VT_CORR_THRESHOLD, ss.str());				
	
	PrintRecord(dnaproj_file, EXPORT_XML_STN_FILE, 
		yesno_string(settings_.o._export_xml_stn_file));					
	PrintRecord(dnaproj_file, EXPORT_DNA_STN_FILE, 
		yesno_string(settings_.o._export_dna_stn_file));					
	PrintRecord(dnaproj_file, EXPORT_SNX_FILE, 
		yesno_string(settings_.o._export_snx_file));					
	
	dnaproj_file << std::endl;

	////////////////////////////////////////////////////////////////////////////////////////////////
	// #plot
	ss.str("");
	ss << section_plot << " (" << PRINT_VAR_PAD << ")";
	PrintRecord(dnaproj_file, ss.str(), "VALUE");
	dnaproj_file << OUTPUTLINE << std::endl;
	
	dnaproj_file << std::endl;

	dnaproj_file.close();
}



}	// namespace dynadjust
