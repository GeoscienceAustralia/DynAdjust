//============================================================================
// Name         : dnageoidwrapper.cpp
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
// Description  : AusGeoid Grid File (NTv2) Interpolation library Executable
//============================================================================

#include <dynadjust/dnageoidwrapper/dnageoidwrapper.hpp>

using namespace dynadjust;

bool CreateNTv2Grid(dna_geoid_interpolation* g, const char* dat_gridfilePath, const n_file_par* grid)
{
	// example:
	// geoid -d ausgeoid09_gda94_v1.01_clip_1x1.dat -c -g ausgeoid_clip_1.0.1.0.gsb --grid-shift radians --grid-version 1.0.1.0 --system-fr ___GDA94 --system-to ___AHD71 --sub-grid-n 1D-grid --creation 21.04.2021 --update 22.04.2021
	//
	
	std::cout << "+ Creating NTv2 geoid grid file from WINTER DAT file format..." << std::endl;

	try {
		g->CreateNTv2File(dat_gridfilePath, grid);
	}
	catch (const NetGeoidException& e) {
		std::cout << std::endl << "- Error: " << e.what() << std::endl;
		return false;
	}
	std::cout << std::endl;

	// Open the new grid file and print its properties
	if (!reportGridProperties(g, grid->filename, grid->filetype))
		return false;

	return true;
}

	
bool ExportNTv2GridToAscii(dna_geoid_interpolation* g, const char* dat_gridfilePath, const char* gridfileType, const char* gridshiftType, const char* exportfileType)
{
	// example:
	// geoid -g ausgeoid_clip_1.0.1.0.gsb --grid-shift radians --export-ntv2-asc
	//

	std::filesystem::path asciiGridFile(dat_gridfilePath);
	std::string outfile = asciiGridFile.filename().string() + "." + exportfileType;

	int ioStatus;
	
	std::cout << std::endl << "+ Exporting NTv2 geoid grid file to " << leafStr<std::string>(outfile) << "... ";

	try {
		g->ExportToAscii(dat_gridfilePath, gridfileType, gridshiftType, outfile.c_str(), &ioStatus);
	}
	catch (const NetGeoidException& e) {
		std::cout << std::endl << "- Error: " << e.what() << std::endl;
		return false;
	}

	std::cout << "done." << std::endl << std::endl;

	// Open the new grid file and print its properties
	if (!reportGridProperties(g, outfile.c_str(), exportfileType))
		return false;

	return true;
}


	
bool ExportNTv2GridToBinary(dna_geoid_interpolation* g, const char* dat_gridfilePath, const char* gridfileType, const char* gridshiftType, const char* exportfileType)
{
	std::filesystem::path asciiGridFile(dat_gridfilePath);
	std::string outfile = asciiGridFile.filename().string() + "." + exportfileType;

	int ioStatus;
	
	std::cout << std::endl << "+ Exporting geoid file to " << leafStr<std::string>(outfile) << "... ";

	try {
		g->ExportToBinary(dat_gridfilePath, gridfileType, gridshiftType, outfile.c_str(), &ioStatus);
	}
	catch (const NetGeoidException& e) {
		std::cout << std::endl << "- Error: " << e.what() << std::endl;
		return false;
	}

	std::cout << "done." << std::endl << std::endl;

	// Open the grid file and print its properties
	if (!reportGridProperties(g, outfile.c_str(), exportfileType))
		return false;

	return true;
}

void ReturnBadStationRecords(dna_geoid_interpolation* g, project_settings& p)
{
	std::string records, filename(p.g.network_name);
	std::string badpointsPath(formPath<std::string>(p.g.output_folder, filename, "int"));
	std::stringstream ss;
	std::ofstream badpoints_log;

	ss << "- Error: Could not open " << filename << " for writing." << std::endl;
	try {
		// Create dynadjust log file.  Throws runtime_error on failure.
		file_opener(badpoints_log, badpointsPath);
	}
	catch (const std::runtime_error& e) {
		ss << e.what();
		throw std::runtime_error(ss.str());
	}
	catch (...) {
		throw std::runtime_error(ss.str());
	}

	// Print formatted header
	print_file_header(badpoints_log, "DYNADJUST GEOID INTERPOLATION LOG FILE");

	badpoints_log << std::setw(PRINT_VAR_PAD) << std::left << "File name:" << badpointsPath << std::endl << std::endl;
	
	badpoints_log << std::setw(PRINT_VAR_PAD) << std::left << "Command line arguments: ";
	badpoints_log << p.n.command_line_arguments << std::endl << std::endl;

	badpoints_log << std::setw(PRINT_VAR_PAD) << std::left << "Network name:" <<  p.g.network_name << std::endl;
	badpoints_log << std::setw(PRINT_VAR_PAD) << std::left << "Stations file:" << safe_absolute_path(p.n.bst_file) << std::endl;
	badpoints_log << std::setw(PRINT_VAR_PAD) << std::left << "Geoid model: " << safe_absolute_path(p.n.ntv2_geoid_file) << std::endl << std::endl;
	badpoints_log << std::setw(PRINT_VAR_PAD) << std::left << "Stations not interpolated:" << g->PointsNotInterpolated() << std::endl;
	badpoints_log << OUTPUTLINE << std::endl << std::endl;
	
	records = g->ReturnBadStationRecords();

	badpoints_log << records << std::endl;

	if (p.g.verbose > 1)
	{
		std::string data("<file record or filename>");
		badpoints_log << std::endl << std::endl <<
			"DYNADJUST GEOID INTERPOLARION ERROR CODES" << std::endl << std::endl <<
			std::setw(PAD) << "Code" << "Description (short and long)" << std::endl <<
			"------------------------------------------------------" << std::endl;
		for (int i=ERR_AUS_BINARY; i<=ERR_INTERPOLATION_TYPE; ++i)
		{
			badpoints_log << 
				std::setw(PAD) << i << 
				std::setw(ZONE) << "Short: " << g->ErrorCaption(i) << std::endl <<
				std::setw(PAD) << " " <<
				std::setw(ZONE) << "Long:  " << g->ErrorString(i, data) << std::endl;
		}
	}
	

	badpoints_log.close();

	std::cout << std::endl << "  See " << badpointsPath << " to view the list of stations for which an" << std::endl <<
		"  N-value could not be interpolated." << std::endl;

}


	
bool createGridIndex(dna_geoid_interpolation* g, const char* gridfilePath, const char* gridfileType, const int& quiet)
{
	if (!quiet)
		std::cout << "+ Opening grid file... ";
	try {
		g->CreateGridIndex(gridfilePath, gridfileType);
	}
	catch (const NetGeoidException& e) {
		std::cout << std::endl << "- Error: " << e.what() << std::endl;
		return false;
	}

	if (!quiet)
		std::cout << "done." << std::endl;

	return true;
}
	

bool reportGridProperties(dna_geoid_interpolation* g, const char* gridfilePath, const char* gridfileType)
{
	n_file_par grid_properties;
	
	try {
		g->ReportGridProperties(gridfilePath, gridfileType, &grid_properties);
	}
	catch (const NetGeoidException& e) {
		std::cout << std::endl << "- Error: " << e.what() << std::endl;
		return false;
	}

	bool isRadians(false);
	std::string shiftType(grid_properties.chGs_type);
	if (iequals(trimstr(shiftType), "radians"))
		isRadians = true;

	std::string formattedLimit;

	std::cout << "+ Grid properties for " << gridfilePath << ":" << std::endl;
	std::cout << "  - GS_TYPE  = " << grid_properties.chGs_type << std::endl;							// grid shift type (GS_TYPE)
	std::cout << "  - VERSION  = " << grid_properties.chVersion << std::endl;							// grid file version (VERSION)
	std::cout << "  - SYSTEM_F = " << grid_properties.chSystem_f << std::endl;						// reference system (SYSTEM_F)
	std::cout << "  - SYSTEM_T = " << grid_properties.chSystem_t << std::endl;						// reference system (SYSTEM_T)
	std::cout << "  - MAJOR_F  = " << std::setprecision(3) << std::fixed << grid_properties.daf << std::endl;	// semi major of from system (MAJOR_F)
	std::cout << "  - MAJOR_T  = " << std::setprecision(3) << std::fixed << grid_properties.dat << std::endl;	// semi major of to system (MAJOR_T)
	std::cout << "  - MINOR_F  = " << std::setprecision(3) << std::fixed << grid_properties.dbf << std::endl;	// semi minor of from system (MINOR_F)
	std::cout << "  - MINOR_T  = " << std::setprecision(3) << std::fixed << grid_properties.dbt << std::endl;	// semi minor of to system (MINOR_T)
	std::cout << "  - NUM_OREC = " << grid_properties.iH_info << std::endl;							// Number of header identifiers (NUM_OREC)
	std::cout << "  - NUM_SREC = " << grid_properties.iSubH_info << std::endl;						// Number of sub-header idents (NUM_SREC)
	std::cout << "  - NUM_FILE = " << grid_properties.iNumsubgrids << std::endl;						// number of subgrids in file (NUM_FILE)

	for (int i=0; i<grid_properties.iNumsubgrids; ++i)
	{
		std::string formattedLimit;
		std::cout << "  - SUBGRID " << i << ":" << std::endl;
		std::cout << "    - SUB_NAME = " << grid_properties.ptrIndex[i].chSubname << std::endl;  		// name of subgrid (SUB_NAME)
		std::cout << "    - PARENT   = " << grid_properties.ptrIndex[i].chParent << std::endl;		// name of parent grid (PARENT)
		std::cout << "    - CREATED  = " << grid_properties.ptrIndex[i].chCreated << std::endl;		// date of creation (CREATED)
		std::cout << "    - UPDATED  = " << grid_properties.ptrIndex[i].chUpdated << std::endl;		// date of last file update (UPDATED)
		
		// In ptrIndex, all values for the limits of a grid are held in seconds, despite
		// whether the grid node records are in radians.

		// lower latitude (S_LAT)
		std::cout << "    - S_LAT    = " << std::right << std::setw(isRadians ? ZONE : REL) << grid_properties.ptrIndex[i].dSlat;						
		formattedLimit = "(" + FormatDmsString(DegtoDms(grid_properties.ptrIndex[i].dSlat / DEG_TO_SEC), 6, true, false) + ")";
		std::cout << std::right << std::setw(MEASR) << formattedLimit << std::endl;								
		
		// upper latitude (N_LAT)
		std::cout << "    - N_LAT    = " << std::right << std::setw(isRadians ? ZONE : REL) << grid_properties.ptrIndex[i].dNlat;
		formattedLimit = "(" + FormatDmsString(DegtoDms(grid_properties.ptrIndex[i].dNlat / DEG_TO_SEC), 6, true, false) + ")";
		std::cout << std::right << std::setw(MEASR) << formattedLimit << std::endl;								
		
		// lower longitude (E_LONG)
		std::cout << "    - E_LONG   = " << std::right << std::setw(isRadians ? ZONE : REL) << grid_properties.ptrIndex[i].dElong;
		formattedLimit = "(" + FormatDmsString(DegtoDms(grid_properties.ptrIndex[i].dElong / DEG_TO_SEC), 6, true, false) + ")";
		std::cout << std::right << std::setw(MEASR) << formattedLimit << std::endl;

		// upper longitude (W_LONG)
		std::cout << "    - W_LONG   = " << std::right << std::setw(isRadians ? ZONE : REL) << grid_properties.ptrIndex[i].dWlong;
		formattedLimit = "(" + FormatDmsString(DegtoDms(grid_properties.ptrIndex[i].dWlong / DEG_TO_SEC), 6, true, false) + ")";
		std::cout << std::right << std::setw(MEASR) << formattedLimit << std::endl;
		
		std::cout << "    - LAT_INC  = " << std::fixed << grid_properties.ptrIndex[i].dLatinc << std::endl;		// latitude interval (LAT_INC)
		std::cout << "    - LONG_INC = " << std::fixed << grid_properties.ptrIndex[i].dLonginc << std::endl;		// longitude interval (LONG_INC)
		std::cout << "    - GS_COUNT = " << std::fixed << std::setprecision(0) << grid_properties.ptrIndex[i].lGscount;						// number of nodes (GS_COUNT)
	}
	std::cout << std::endl;
	return true;
}

bool InterpolateGridPoint(dna_geoid_interpolation* g, geoid_point* apInterpolant, 
	const int& method, const int& coordinate_format, const std::string& inputLatitude, const std::string& inputLongitude)
{
	try {
		if (method == BICUBIC)
			g->BiCubicTransformation(apInterpolant);
		else
			g->BiLinearTransformation(apInterpolant);
	}
	catch (const NetGeoidException& e) {
		std::cout << std::endl << "- Error: " << e.what() << std::endl;
		return false;
	}

	if (apInterpolant->cVar.IO_Status != ERR_TRANS_SUCCESS)
		return false;
	
	std::cout << "+ Interpolation results for ";
	std::cout << inputLatitude << ", " << inputLongitude;
	if (coordinate_format == DMS)
		//cout << std::fixed << std::setprecision(6) << DegtoDms<double>(apInterpolant->cVar.dLatitude) << ", " << DegtoDms<double>(apInterpolant->cVar.dLongitude) << " (ddd.mmssss):" << std::endl;
		std::cout << " (ddd.mmssss):" << std::endl;
	else
		//cout << std::fixed << std::setprecision(6) << apInterpolant->cVar.dLatitude << ", " << apInterpolant->cVar.dLongitude << " (ddd.dddddd):" << std::endl;
		std::cout << " (ddd.dddddd):" << std::endl;
	
	std::cout << std::endl;

	std::cout << "  N value          = " << std::setw(6) << 
		std::right << std::setprecision(3) << apInterpolant->gVar.dN_value << " metres" << std::endl;			// N value
	std::cout << "  Deflections:" << std::endl;
	std::cout << "  - Prime meridian = " << std::setw(6) << 
		std::right << std::fixed << std::setprecision(2) << apInterpolant->gVar.dDefl_meridian << " seconds" << std::endl;			// N value
	std::cout << "  - Prime vertical = " << std::setw(6) << 
		std::right << apInterpolant->gVar.dDefl_primev << " seconds" << std::endl;			// N value
	std::cout << std::endl;
	return true;

} // InterpolateGridPoint
	

bool InterpolateGridPointFile(dna_geoid_interpolation* g, const char* inputfilePath, 
	const int& method, const int EllipsoidtoOrtho, const int& coordinate_format, 
	bool exportDnaGeoidFile, const char* dnageofilePath, std::string& outputfilePath)
{
	std::filesystem::path inputFile(inputfilePath);
	if (inputFile.has_extension())
		outputfilePath = inputFile.stem().string() + "_out" + inputFile.extension().string();
	else
	{
		outputfilePath = inputfilePath; 
		outputfilePath.append("_out");
	}

	char outfilePath[601];
	strcpy(outfilePath, outputfilePath.c_str());

	try {
		g->FileTransformation(inputfilePath, outfilePath, 
			method, EllipsoidtoOrtho, coordinate_format, 
			exportDnaGeoidFile, dnageofilePath);
	}
	catch (const NetGeoidException& e) {
		std::cout << std::endl << "- Error: " << e.what() << std::endl;
		return false;
	}

	return true;

} // InterpolateGridPointFile


bool InterpolateGridBinaryStationFile(dna_geoid_interpolation* g, const std::string& bstnfilePath,
	const int& method, bool convertHeights, 
	bool exportDnaGeoidFile, const char* dnageofilePath)
{
	try {
		g->PopulateBinaryStationFile(bstnfilePath, method, convertHeights, 
			exportDnaGeoidFile, dnageofilePath);
	}
	catch (const NetGeoidException& e) {
		std::cout << std::endl << "- Error: " << e.what() << std::endl;
		return false;
	}

	return true;

} // InterpolateGridBinaryStationFile


std::string GetFileType(const std::string inputfilePath)
{
	std::filesystem::path inputFile(inputfilePath);
	if (inputFile.has_extension())
		return inputFile.extension().string();
	else
		return "";
}

int ParseCommandLineOptions(const int& argc, char* argv[], const boost::program_options::variables_map& vm, project_settings& p)
{
	// capture command line arguments
	for (int cmd_arg(0); cmd_arg<argc; ++cmd_arg)
	{
		 p.n.command_line_arguments += argv[cmd_arg];
		 p.n.command_line_arguments += " ";
	}

	// Has the user supplied a project file?
	if (vm.count(PROJECT_FILE))
	{
		if (std::filesystem::exists(p.g.project_file))
		{
			try {
				CDnaProjectFile projectFile(p.g.project_file, geoidSetting);
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

	if (!vm.count(NETWORK_NAME) &&				// User wants to populate a binary station file for a particular DynAdjust project
		!vm.count(INTERACTIVE) &&				// User wants to interpolate geoid information from the command line
		!vm.count(CREATE_NTV2) &&				// User wants to create a NTv2 grid file
		!vm.count(SUMMARY) &&					// User wants to print to the screen the details of the grid file
		!vm.count(INPUT_FILE) &&				// User wants to interpolate geoid information in text file mode
		!vm.count(DAT_FILEPATH) &&				// User supplied file path to WINTER DAT file
		!vm.count(EXPORT_NTV2_ASCII_FILE) &&	// User wants to export a geoid file to ASCII
		!vm.count(EXPORT_NTV2_BINARY_FILE))		// User wants to export a geoid file to binary
	{
		std::cout << std::endl << "- Nothing to do - no standard options specified. " << std::endl << std::endl;
		return EXIT_FAILURE;
	}

	// check mandatory arguments for NTv2 file creation
	if (vm.count(CREATE_NTV2) && !vm.count(DAT_FILEPATH))
	{
		std::cout << std::endl << "- Error: no dat file specified. " << std::endl << std::endl;
		return EXIT_FAILURE;
	}

	// check mandatory arguments for interactive mode
	if (vm.count(INTERACTIVE) && !vm.count(LATITUDE))
	{
		std::cout << std::endl << "- Error: Interpolation latitide not specified. " << std::endl << std::endl;
		return EXIT_FAILURE;
	}

	// check mandatory arguments for interactive mode
	if (vm.count(INTERACTIVE) && !vm.count(LONGITUDE))
	{
		std::cout << std::endl << "- Error: Interpolation longitude not specified. " << std::endl << std::endl;
		return EXIT_FAILURE;
	}

	// check mandatory arguments for text file mode
	if (!vm.count(INPUT_FILE) && !vm.count(NTV2_FILEPATH) && !vm.count(CREATE_NTV2))
	{
		std::cout << std::endl << "- Error: Interpolation input file not specified. " << std::endl << std::endl;
		return EXIT_FAILURE;
	}

	// check mandatory arguments for NTv2 file summary
	if (vm.count(SUMMARY) && !vm.count(NTV2_FILEPATH))
	{
		std::cout << std::endl << "- Error: No NTv2 grid file specified. " << std::endl << std::endl;
		return EXIT_FAILURE;
	}

	// check mandatory arguments for DynAdjust station file population
	if (vm.count(NETWORK_NAME) && !vm.count(NTV2_FILEPATH))
	{
		std::cout << std::endl << "- Error: No NTv2 grid file specified. " << std::endl << std::endl;
		return EXIT_FAILURE;
	}

	// binary station file location (input)
	if (vm.count(NETWORK_NAME))
	{
		p.g.project_file = formPath<std::string>(p.g.output_folder, p.g.network_name, "dnaproj");

		// define bst file name
		p.n.bst_file = formPath<std::string>(p.g.input_folder, p.g.network_name, "bst");
		p.n.file_mode = 1;
		if (!std::filesystem::exists(p.n.bst_file))
		{
			// Look for it in the input folder
			p.n.bst_file = formPath<std::string>(p.g.input_folder, leafStr<std::string>(p.n.bst_file));

			if (!std::filesystem::exists(p.n.bst_file))
			{
				std::cout << std::endl << "- Error: ";  
				std::cout << "Binary station file " << p.n.bst_file << " does not exist." << std::endl << std::endl;  
				return EXIT_FAILURE;
			}
		}
	}

	// Geoid DAT grid file file location (input)
	if (vm.count(DAT_FILEPATH))
	{
		if (!std::filesystem::exists(p.n.rdat_geoid_file))
		{
			// Look for it in the input folder
			p.n.rdat_geoid_file = formPath<std::string>(p.g.input_folder, leafStr<std::string>(p.n.rdat_geoid_file));

			if (!std::filesystem::exists(p.n.rdat_geoid_file))
			{
				std::cout << std::endl << "- Error: ";  
				std::cout << "WINTER DAT grid file " << p.n.rdat_geoid_file << " does not exist." << std::endl << std::endl;  
				return EXIT_FAILURE;
			}
		}
	}

	// Is geoid to run in file mode?
	if (vm.count(INPUT_FILE))
	{
		p.n.file_mode = 1;

		// Geoid DAT grid file file location (input)
		if (!std::filesystem::exists(p.n.input_file))
		{
			// Look for it in the input folder
			p.n.input_file = formPath<std::string>(p.g.input_folder, leafStr<std::string>(p.n.input_file));

			if (!std::filesystem::exists(p.n.input_file))
			{
				std::cout << std::endl << "- Error: ";  
				std::cout << "Input coordinates text file " << leafStr<std::string>(p.n.input_file) << " does not exist." << std::endl << std::endl;  
				return EXIT_FAILURE;
			}
		}
	}

	if (vm.count(EXPORT_GEO_FILE))
	{
		if (vm.count(NETWORK_NAME))
			p.n.geo_file = formPath<std::string>(p.g.output_folder, p.g.network_name, "geo");	// dna geoid file
		p.n.export_dna_geo_file = 1;
	}

	if (vm.count(DDEG_FORMAT))
		p.n.coordinate_format = DDEG;
	
	// Change behaviour of geoid to always convert orthometric heights to ellipsoid
	// deprecate this command
	//if (vm.count(CONVERT_BST_HT))
	p.n.convert_heights = 1;
	//else
	//	p.n.convert_heights = 0;

	return EXIT_SUCCESS;
}

int main(int argc, char* argv[])
{
	boost::program_options::variables_map vm;
	boost::program_options::positional_options_description positional_options;
	
	boost::program_options::options_description standard_options("+ " + std::string(ALL_MODULE_STDOPT), PROGRAM_OPTIONS_LINE_LENGTH);
	boost::program_options::options_description ntv2_options("+ " + std::string(GEOID_MODULE_NTV2), PROGRAM_OPTIONS_LINE_LENGTH);
	boost::program_options::options_description interpolate_options("+ " + std::string(GEOID_MODULE_INTERPOLATE), PROGRAM_OPTIONS_LINE_LENGTH);
	boost::program_options::options_description interactive_options("+ " + std::string(GEOID_MODULE_INTERACTIVE), PROGRAM_OPTIONS_LINE_LENGTH);
	boost::program_options::options_description file_interpolate_options("+ " + std::string(GEOID_MODULE_FILE), PROGRAM_OPTIONS_LINE_LENGTH);
	boost::program_options::options_description export_options("+ " + std::string(ALL_MODULE_EXPORT), PROGRAM_OPTIONS_LINE_LENGTH);
	boost::program_options::options_description generic_options("+ " + std::string(ALL_MODULE_GENERIC), PROGRAM_OPTIONS_LINE_LENGTH);

	std::string cmd_line_usage("+ ");
	cmd_line_usage.append(__BINARY_NAME__).append(" usage:  ").append(__BINARY_NAME__).append(" [options]");
	boost::program_options::options_description allowable_options(cmd_line_usage, PROGRAM_OPTIONS_LINE_LENGTH);

	n_file_par ntv2;
	ntv2.ptrIndex = new n_gridfileindex[1];
	geoid_point apInterpolant;

	project_settings p;

	std::string cmd_line_banner, gs_type("seconds"), version("1.0.0.0"), system_f("GDA94   "), system_t("AHD_1971");
	std::string subgridname("AUSGEOID"), parent(""), created(""), updated("");

	std::string inputLatitude, inputLongitude;
	
	fileproc_help_header(&cmd_line_banner);
	p.g.project_file = "";

	try {
		// Declare a group of options that will be 
		// allowed only on command line		
		standard_options.add_options()
			(PROJECT_FILE_P, boost::program_options::value<std::string>(&p.g.project_file),
				"Project file name. Full path to project file. If none specified, a new file is created using input-folder and network-name.")
			(NETWORK_NAME_N, boost::program_options::value<std::string>(&p.g.network_name), 
				std::string("Network name. If [" + std::string(NETWORK_NAME) + "].bst exists, all records within the binary station file will be populated with N value and deflections of the vertical.").c_str())
			(INPUT_FOLDER_I, boost::program_options::value<std::string>(&p.g.input_folder),
				"Path containing all input files.")
			(OUTPUT_FOLDER_O, boost::program_options::value<std::string>(&p.g.output_folder),		// default is ./,
				"Path for all output files.")
			;

		interpolate_options.add_options()
			(INTERACTIVE_E, "Interpolate geoid information using coordinates provided on the command line.")
			(INPUT_FILE_T, boost::program_options::value<std::string>(&p.n.input_file),
				"Interpolate geoid information using coordinates contained in a text file. "
				"arg is the path of the input text file. "
				"The supported text file formats include formatted text (*.txt) and comma separated values (*.csv) files. "
				"Refer to the User's Guide for file format information.")
			(METHOD_M, boost::program_options::value<UINT16>(&p.n.interpolation_method),
				"Interpolation method.\n  0  Bi-linear\n  1  Bi-cubic (default)")
			(CONVERT_BST_HT, 
				"DEPRECATED. If a user-supplied height in the binary file is orthometric, the height will be converted to ellipsoidal automatically.")
			(DDEG_FORMAT, "Specify input coordinates in decimal degrees (dd.dddddd).  Default is degrees, minutes and seconds (dd.mmssss).")
			(CREATE_NTV2_C, "Create NTv2 grid file from standard DAT file.")
			(NTV2_FILEPATH_G, boost::program_options::value<std::string>(&p.n.ntv2_geoid_file), "Full file path of the NTv2 grid file.")
			(SUMMARY_U, "Print a summary of the grid file.")
			;

		// Declare a group of options that will be 
		// allowed both on command line and in
		// config file        
		ntv2_options.add_options()
			(DAT_FILEPATH_D, boost::program_options::value<std::string>(&p.n.rdat_geoid_file), 
				"File path of the WINTER DAT grid file.")
			(NTV2_GS_TYPE, boost::program_options::value<std::string>(&gs_type),
				"Units in which the grid parameters and deflections of the vertical will be stored. arg is either 'seconds' or 'radians'. Default is seconds.")
			(NTV2_VERSION, boost::program_options::value<std::string>(&version),
				"Grid file version. Default is 1.0.0.0.")
			(NTV2_SYSTEM_F, boost::program_options::value<std::string>(&system_f),
				"The 'From' reference system. Default is GDA94.")
			(NTV2_SYSTEM_T, boost::program_options::value<std::string>(&system_t),
				"The 'To' reference system. Default is AHD_1971")
			(NTV2_MAJOR_F, boost::program_options::value<double>(&ntv2.daf),
				"Semi major of 'From' system. Default is 6378137.000")
			(NTV2_MAJOR_T, boost::program_options::value<double>(&ntv2.dat),
				"Semi major of 'To' system. Default is 6378137.000")
			(NTV2_MINOR_F, boost::program_options::value<double>(&ntv2.dbf),
				"Semi minor of 'From' system. Default is 6356752.314")
			(NTV2_MINOR_T, boost::program_options::value<double>(&ntv2.dbt),
				"Semi minor of 'To' system. Default is 6356752.314")
			(NTV2_SUB_NAME, boost::program_options::value<std::string>(&subgridname),
				"The name of the sub-grid. Default is AUSGEOID")
			(NTV2_CREATED, boost::program_options::value<std::string>(&created), 
				"Date of geoid model creation. arg is a dot delimited string \"dd.mm.yyyy\". Default is today's date if no value is supplied.")
			(NTV2_UPDATED, boost::program_options::value<std::string>(&updated), 
				"Date of last file update. arg is a dot delimited string \"dd.mm.yyyy\". Default is today's date if no value is supplied.")
			;

		interactive_options.add_options()
			(LATITUDE, boost::program_options::value<std::string>(&inputLatitude),
				"Latitude of the interpolant. Default is degrees, minutes and seconds (dd.mmssss).")
			(LONGITUDE, boost::program_options::value<std::string>(&inputLongitude),
				"Longitude of the interpolant. Default is degrees, minutes and seconds (dd.mmssss).")
			;

		file_interpolate_options.add_options()
			(DIRECTION_R, boost::program_options::value<UINT16>(&p.n.ellipsoid_to_ortho), 
				"Conversion of heights:\n  0  Orthometric to ellipsoid (default)\n  1  Ellipsoid to orthometric")
			;

		// Declare a group of options that will be 
		// allowed both on command line and in
		// config file        
		export_options.add_options()
			(EXPORT_GEO_FILE, 
				"Create a DNA geoid file from interpolated geoid information.")
			(EXPORT_NTV2_ASCII_FILE,
				"Export a binary NTv2 geoid file to ASCII (.asc) format.")
			(EXPORT_NTV2_BINARY_FILE,
				"Export an ASCII NTv2 geoid file to binary (.gsb) format.")
			;

		generic_options.add_options()
			(VERBOSE, boost::program_options::value<UINT16>(&p.g.verbose),
				std::string("When importing geoid information into a project, print the stations for which an N-value could not be interpolated to a log (*.int) file.\n  0: No information (default)\n  1: Helpful information\n  2: Extended information").c_str())
				(QUIET,
					std::string("Suppresses all explanation of what ").append(__BINARY_NAME__).append(" is doing unless an error occurs").c_str())
					(VERSION_V, "Display the current program version")
			(HELP_H, "Show this help message")
			(HELP_MODULE_H, boost::program_options::value<std::string>(),
				"Provide help for a specific help category.")
			;

		allowable_options.add(standard_options).add(interpolate_options).add(ntv2_options).add(interactive_options).add(file_interpolate_options).add(export_options).add(generic_options);

		// add "positional options" to handle command line tokens which have no option name
		positional_options.add(NETWORK_NAME, -1);

		boost::program_options::command_line_parser parser(argc, argv);
		store(parser.options(allowable_options).positional(positional_options).run(), vm);
		notify(vm);
	} 
	catch (const std::exception& e) {
		std::cout << std::endl << "- Error: " << e.what() << std::endl;
		std::cout << cmd_line_banner << allowable_options << std::endl;
		return EXIT_FAILURE;
	}
	catch (...) 
	{
		std::cout << std::endl << "- Exception of unknown type!\n";
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
		else if (str_upper<std::string, char>(GEOID_MODULE_NTV2).find(help_text) != std::string::npos) {
			std::cout << ntv2_options << std::endl;
		}
		else if (str_upper<std::string, char>(GEOID_MODULE_INTERPOLATE).find(help_text) != std::string::npos) {
			std::cout << interpolate_options << std::endl;
		}
		else if (str_upper<std::string, char>(GEOID_MODULE_INTERACTIVE).find(help_text) != std::string::npos) {
			std::cout << interactive_options << std::endl;
		}
		else if (str_upper<std::string, char>(GEOID_MODULE_FILE).find(help_text) != std::string::npos) {
			std::cout << file_interpolate_options << std::endl;
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

	if (ParseCommandLineOptions(argc, argv, vm, p) != EXIT_SUCCESS)
		return EXIT_FAILURE;

	// Has the user supplied a project file?
	bool userSuppliedProjectFile(false);
	if (!p.g.project_file.empty())
		userSuppliedProjectFile = true;

	// Capture NTv2 file path and extension
	// grid file path not supplied.  Generate name from dat file
	if (p.n.ntv2_geoid_file.empty())
	{
		std::filesystem::path gsbFile(p.n.rdat_geoid_file);
		p.n.ntv2_geoid_file = gsbFile.stem().string() + gsbFile.extension().string() + ".gsb";
		strcpy(ntv2.filename, p.n.ntv2_geoid_file.c_str());
		strcpy(ntv2.filetype, GSB);
	}
	else
	{
		std::string extension(GetFileType(p.n.ntv2_geoid_file));
		if (extension.empty())
		{
			std::cout << std::endl << "- Error: NTv2 grid file type cannot be determined from a file without a file extension. " << std::endl << std::endl;
			return EXIT_FAILURE;
		}

		if (!iequals(extension.substr(1), GSB) &&
			!iequals(extension.substr(1), ASC))
		{
			std::cout << std::endl << "- Error: NTv2 grid file type cannot be determined from file extension \"" << extension << "\"." << std::endl << 
				"         Supported types are ." << GSB << " and ." << ASC << " only." << std::endl << std::endl;
			return EXIT_FAILURE;
		}

		if (vm.count(EXPORT_NTV2_ASCII_FILE) &&
			iequals(extension.substr(1), ASC))
		{
			std::cout << std::endl << "- Error: Export to ASCII NTv2 grid file option only supported for " << GSB " grid files." << std::endl << std::endl;
			return EXIT_FAILURE;
		}

		if (vm.count(EXPORT_NTV2_BINARY_FILE) &&
			iequals(extension.substr(1), GSB))
		{
			std::cout << std::endl << "- Error: Export to Binary NTv2 grid file option only supported for " << ASC " grid files." << std::endl << std::endl;
			return EXIT_FAILURE;
		}

		strcpy(ntv2.filename, p.n.ntv2_geoid_file.c_str());
		strcpy(ntv2.filetype, extension.substr(1).c_str());
	}	
	
	if (vm.count(QUIET))
		p.g.quiet = 1;

	if (vm.count(NTV2_GS_TYPE))
	{
		gs_type = trimstr(gs_type);
		// Unknown type?
		if (!iequals(gs_type, "seconds") && !iequals(gs_type, "radians"))
			gs_type = "seconds";
		str_toupper<int>(gs_type);
	}
	
	if (!p.g.quiet)
	{
		std::cout << std::endl << cmd_line_banner;

		std::cout << "+ Options:" << std::endl; 
		
		if (vm.count(NETWORK_NAME))
			std::cout << std::setw(PRINT_VAR_PAD) << std::left << "  Network name: " <<  p.g.network_name << std::endl;

		if (p.n.file_mode || vm.count(CREATE_NTV2))
		{
			std::cout << std::setw(PRINT_VAR_PAD) << std::left << "  Input folder: " << p.g.input_folder << std::endl;
			std::cout << std::setw(PRINT_VAR_PAD) << std::left << "  Output folder: " << p.g.output_folder << std::endl;

			if (!p.n.bst_file.empty())
			{
				std::cout << std::setw(PRINT_VAR_PAD) << std::left << "  Binary station file: " << p.n.bst_file << std::endl;
				//cout << std::setw(PRINT_VAR_PAD) << std::left << "  Convert orthometric heights: ";
				//if (p.n.convert_heights)
				//	std::cout << "Yes" << std::endl;
				//else
				//	std::cout << "No" << std::endl;
			}

			if (!p.n.input_file.empty())
				std::cout << std::setw(PRINT_VAR_PAD) << std::left << "  ASCII file: " << p.n.input_file << std::endl;
		}
		
		std::cout << std::setw(PRINT_VAR_PAD) << std::left << "  Geoid grid file: " <<  ntv2.filename << std::endl;
		
		// Not applicable for project file use
		if (vm.count(CREATE_NTV2))
		{
			std::cout << std::setw(PRINT_VAR_PAD) << std::left << "  WINTER DAT file: " << leafStr<std::string>(p.n.rdat_geoid_file) << std::endl;

			if (vm.count(NTV2_GS_TYPE))
				std::cout << std::setw(PRINT_VAR_PAD) << std::left << "  Grid shift type: " << gs_type.c_str() << std::endl;
			if (vm.count(NTV2_VERSION))
				std::cout << std::setw(PRINT_VAR_PAD) << std::left << "  Grid file version: " << version.c_str() << std::endl;
			if (vm.count(NTV2_SYSTEM_F))
				std::cout << std::setw(PRINT_VAR_PAD) << std::left << "  From reference system: " << system_f.c_str() << std::endl;
			if (vm.count(NTV2_SYSTEM_T))
				std::cout << std::setw(PRINT_VAR_PAD) << std::left << "  To reference system: " << system_t.c_str() << std::endl;
			if (vm.count(NTV2_MAJOR_F))
				std::cout << std::setw(PRINT_VAR_PAD) << std::left << "  From semi-major: " << std::fixed << std::setprecision(3) << ntv2.daf << std::endl;
			if (vm.count(NTV2_MAJOR_T))
				std::cout << std::setw(PRINT_VAR_PAD) << std::left << "  To semi-major: " << ntv2.dat << std::endl;
			if (vm.count(NTV2_MINOR_F))
				std::cout << std::setw(PRINT_VAR_PAD) << std::left << "  From semi-minor: " << ntv2.dbf << std::endl;
			if (vm.count(NTV2_MINOR_T))
				std::cout << std::setw(PRINT_VAR_PAD) << std::left << "  To semi-minor: " << ntv2.dbt << std::endl;
			if (vm.count(NTV2_SUB_NAME))
				std::cout << std::setw(PRINT_VAR_PAD) << std::left << "  Sub-grid name: " << subgridname.c_str() << std::endl;
			
			boost::gregorian::date creationDate, updateDate;
			if (created.empty())
				created = "today";
			
			if (updated.empty())
				updated = "today";
			
			creationDate = dateFromString<boost::gregorian::date>(created);
			updateDate = dateFromString<boost::gregorian::date>(updated);

			// Print dates by default
			//if (vm.count(NTV2_CREATED))
			std::cout << std::setw(PRINT_VAR_PAD) << std::left << "  Date of file creation: " << 
				stringFromDate<boost::gregorian::date>(creationDate, "%d %B %Y") << std::endl;
			//if (vm.count(NTV2_UPDATED))
			std::cout << std::setw(PRINT_VAR_PAD) << std::left << "  Date of file update: " << 
				stringFromDate<boost::gregorian::date>(updateDate, "%d %B %Y") << std::endl;

			created = stringFromDate<boost::gregorian::date>(creationDate, "%d%m%Y");
			updated = stringFromDate<boost::gregorian::date>(updateDate, "%d%m%Y");

		}

		if (vm.count(EXPORT_NTV2_ASCII_FILE) ||
			vm.count(EXPORT_NTV2_BINARY_FILE))
		{
			std::cout << std::setw(PRINT_VAR_PAD) << std::left << "  Export NTv2 grid file to: ";
			if (vm.count(EXPORT_NTV2_ASCII_FILE))
				std::cout << "ASCII (." << ASC << ")" << std::endl;
			if (vm.count(EXPORT_NTV2_BINARY_FILE))
				std::cout << "Binary (." << GSB << ")" << std::endl;

			if (vm.count(NTV2_GS_TYPE))
				std::cout << std::setw(PRINT_VAR_PAD) << std::left << "  Grid shift type: " << gs_type.c_str() << std::endl;
		}

		if (p.n.file_mode || vm.count(INTERACTIVE))
		{
			std::cout << std::setw(PRINT_VAR_PAD) << std::left << "  Interpolation method: ";
			if (p.n.interpolation_method == BICUBIC)
				std::cout << "Bi-cubic" << std::endl;
			else
				std::cout << "Bi-linear" << std::endl;

			std::cout << std::setw(PRINT_VAR_PAD) << std::left << "  Input coordinate format: ";
			if (p.n.coordinate_format == DDEG)
				std::cout << "Decimal degrees" << std::endl;
			else
				std::cout << "Degrees minutes seconds" << std::endl;
			
			if (!vm.count(INTERACTIVE))
			{
				std::cout << std::setw(PRINT_VAR_PAD) << std::left << "  Transformation direction: ";
				if (p.n.ellipsoid_to_ortho == 0)
					std::cout << "Orthometric to ellipsoid" << std::endl;
				else
					std::cout << "Ellipsoid to orthometric" << std::endl;
			}
		}

		if (p.n.export_dna_geo_file)
			std::cout << std::setw(PRINT_VAR_PAD) << std::left << "  Export to DNA geoid file: " << "Yes" << std::endl;

		std::cout << std::endl;

		if (vm.count(CONVERT_BST_HT))
			std::cout << "- Warning: The '--" << CONVERT_BST_HT << "' option has been deprecated. Orthometric" << std::endl <<
				"  heights in the binary file will be converted to ellipsoidal by default, " << std::endl <<
			    "  unless the transformation direction has been modified by supplying the " << std::endl << 
			    "  '--" << DIRECTION << "' option with an argument of 1." << std::endl;

		std::cout << std::endl;
	
		// File interpolation mode...
		if (p.n.file_mode)
		{
			if (!p.n.bst_file.empty())
				std::cout << "+ Binary station file interpolation mode." << std::endl << std::endl;
			else
				std::cout << "+ ASCII file interpolation mode." << std::endl << std::endl;
		}			
	}

	dna_geoid_interpolation g;
	
	// Not applicable for project file use
	if (vm.count(CREATE_NTV2))
	{
		if (vm.count(NTV2_GS_TYPE))
			strcpy(ntv2.chGs_type, gs_type.c_str());
		if (vm.count(NTV2_VERSION))
			strcpy(ntv2.chVersion, version.c_str());
		if (vm.count(NTV2_SYSTEM_F))
			strcpy(ntv2.chSystem_f, system_f.c_str());
		if (vm.count(NTV2_SYSTEM_T))
			strcpy(ntv2.chSystem_t, system_t.c_str());
		
		// subgrid header string handlers
		if (vm.count(NTV2_SUB_NAME))
			strcpy(ntv2.ptrIndex[0].chSubname, subgridname.c_str());
		strcpy(ntv2.ptrIndex[0].chParent, "NONE    ");		// don't give the user the option - must be "NONE    " !!!
		strcpy(ntv2.ptrIndex[0].chCreated, created.c_str());
		strcpy(ntv2.ptrIndex[0].chUpdated, updated.c_str());

		// Create a new grid file from AusGeoid DAT file
		if (!CreateNTv2Grid(&g, p.n.rdat_geoid_file.c_str(), &ntv2))
		{
			if (ntv2.ptrIndex)
				delete [] ntv2.ptrIndex;
			return EXIT_FAILURE;
		}

		if (ntv2.ptrIndex)
			delete [] ntv2.ptrIndex;

		std::cout << std::endl << "+ Geoid file creation completed successfully." << std::endl << std::endl;

		return EXIT_SUCCESS;
	}
	else if (vm.count(EXPORT_NTV2_ASCII_FILE))
	{
		if (p.n.ntv2_geoid_file.empty())
		{
			std::cout << std::endl << "- Error: No NTv2 grid file specified. " << std::endl << std::endl;
			return EXIT_FAILURE;
		}

		if (vm.count(NTV2_GS_TYPE))
			strcpy(ntv2.chGs_type, gs_type.c_str());

		if (!ExportNTv2GridToAscii(&g, ntv2.filename, ntv2.filetype, ntv2.chGs_type, ASC))
			return EXIT_FAILURE;
		
		std::cout << "+ Geoid file creation completed successfully." << std::endl << std::endl;

		return EXIT_SUCCESS;
	}
	else if (vm.count(EXPORT_NTV2_BINARY_FILE))
	{
		if (p.n.ntv2_geoid_file.empty())
		{
			std::cout << std::endl << "- Error: No NTv2 grid file specified. " << std::endl << std::endl;
			return EXIT_FAILURE;
		}

		if (vm.count(NTV2_GS_TYPE))
			strcpy(ntv2.chGs_type, gs_type.c_str());

		if (!ExportNTv2GridToBinary(&g, ntv2.filename, ntv2.filetype, ntv2.chGs_type, GSB))
			return EXIT_FAILURE;

		std::cout << "+ Geoid file creation completed successfully." << std::endl << std::endl;

		return EXIT_SUCCESS;
	}

	// Not applicable for project file use
	if (vm.count(SUMMARY))
	{
		// Open the grid file and print its properties
		if (!reportGridProperties(&g, ntv2.filename, ntv2.filetype))
			return EXIT_FAILURE;

		return EXIT_SUCCESS;
	}

	// Not applicable for project file use
	if (vm.count(INTERACTIVE))
	{
		apInterpolant.cVar.dLatitude = DoubleFromString<double>(inputLatitude);
		apInterpolant.cVar.dLongitude = DoubleFromString<double>(inputLongitude);

		std::stringstream ssInput;
		ssInput << inputLatitude << ", " << inputLongitude;
		g.SetInputCoordinates(ssInput.str());

		if  (p.n.coordinate_format == DMS)
		{
			apInterpolant.cVar.dLatitude = DmstoDeg<double>(apInterpolant.cVar.dLatitude);
			apInterpolant.cVar.dLongitude = DmstoDeg<double>(apInterpolant.cVar.dLongitude);
		}

		// Get geoid values for a point
		if (!createGridIndex(&g, ntv2.filename, ntv2.filetype, p.g.quiet))
			return EXIT_FAILURE;

		if (!InterpolateGridPoint(&g, &apInterpolant, p.n.interpolation_method,
			p.n.coordinate_format, inputLatitude, inputLongitude))
		{
			std::cout << std::endl;
			if (apInterpolant.cVar.IO_Status == ERR_FINDSUBGRID_OUTSIDE)
				reportGridProperties(&g, ntv2.filename, ntv2.filetype);

			return EXIT_FAILURE;
		}
		return EXIT_SUCCESS;
	}
	 
	// Should a file be processed?
	if (!p.n.file_mode)
		return EXIT_SUCCESS;
	
	// Get geoid values for a point
	if (!createGridIndex(&g, ntv2.filename, ntv2.filetype, p.g.quiet))
		return EXIT_FAILURE;

	if (!p.g.quiet)
	{
		std::cout << "+ Interpolating geoid components";
		if (!p.n.bst_file.empty() && p.n.convert_heights)
			std::cout << " and reducing" << std::endl <<
				"  heights to the ellipsoid";
		std::cout << "... ";
	}
	
	cpu_timer time;

	char dnageoFile[601], *geoFileptr;
	geoFileptr = NULL;
	memset(dnageoFile, '\0', sizeof(dnageoFile));
	std::string outputfilePath;
		
	if (!p.n.geo_file.empty())
	{
		snprintf(dnageoFile, sizeof(dnageoFile), "%s", p.n.geo_file.c_str());
		geoFileptr = dnageoFile;
	}
		
	if (!p.n.bst_file.empty())
	{
		// populate binary station file with geoid separation and deflections
		if (!InterpolateGridBinaryStationFile(&g, p.n.bst_file.c_str(), p.n.interpolation_method, 
			p.n.convert_heights ? true : false, 
			p.n.export_dna_geo_file? true : false, geoFileptr))
				return EXIT_FAILURE;
	}
	else if (!p.n.input_file.empty())
	{
		// Ascii text file (similar to old winter format)
		if (!InterpolateGridPointFile(&g, p.n.input_file.c_str(), p.n.interpolation_method, 
			p.n.ellipsoid_to_ortho, p.n.coordinate_format,
			p.n.export_dna_geo_file? true : false, geoFileptr,
			outputfilePath))
				return EXIT_FAILURE;
	}
	else
		return EXIT_SUCCESS;
		
	if (!p.g.quiet)
	{
		std::cout << "done." << std::endl;
		std::cout << "+ Interpolated data for " << g.PointsInterpolated();
		if (g.PointsInterpolated() == 1)
			std::cout << " point." << std::endl;
		else
			std::cout << " points." << std::endl;

		if (g.PointsNotInterpolated() > 0)
		{
			std::cout << "- Warning: Data for " << g.PointsNotInterpolated();

			// Is this wrapper being called to update DynAdjust station file?
			if (!p.n.bst_file.empty())
			{
				if (g.PointsNotInterpolated() == 1)
					std::cout << " station";
				else
					std::cout << " stations";

				std::cout << " could not be interpolated.  ";
					
				if (p.g.verbose > 0)
				{
					try {
						ReturnBadStationRecords(&g, p);
					}
					catch (const std::runtime_error& e) {
						// print error message, and continue
						std::cout << "- Error: " << e.what() << std::endl;
					}
				}
				else
					std::cout << "To view the list of stations " << std::endl <<
						"  for which an N-value could not be interpolated, call " << __BINARY_NAME__ << " with --verbose-level 1." << std::endl;
			}
			// If this point is reached, then this wrapper must have been called
			// to interpolate points in interactive mode.
			else //if (!p.n.input_file.empty())
			{
				if (g.PointsNotInterpolated() == 1)
					std::cout << " point";
				else
					std::cout << " points";

				std::cout << " could not be interpolated." << std::endl;
				std::cout << "  See " << outputfilePath << " for more information." << std::endl;
			}
		}
	}

	// Look for a project file (if the user has specified a network name.  
	// If it exists, open it, load it and update the geoid settings.
	if (userSuppliedProjectFile)
	{
		CDnaProjectFile projectFile;
		if (std::filesystem::exists(p.g.project_file))
			projectFile.LoadProjectFile(p.g.project_file);

		// Print the project file. If it doesn't exist, it will be created.
		projectFile.UpdateSettingsGeoid(p);
		projectFile.PrintProjectFile();
	}

	if (p.g.quiet)
		return EXIT_SUCCESS;

	// wall time is in nanoseconds
	// cout << time.elapsed().wall << std::endl << std::endl;
	auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(time.elapsed().wall).count();
	boost::posix_time::milliseconds elapsed_time(elapsed_ms);
	std::cout << std::endl << formatedElapsedTime<std::string>(&elapsed_time, "+ Geoid file interpolation took ") << std::endl << std::endl;
	
	return EXIT_SUCCESS;
}
