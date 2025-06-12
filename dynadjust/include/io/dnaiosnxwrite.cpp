//============================================================================
// Name         : dnaiosnxwrite.cpp
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
// Description  : DynAdjust SINEX file write operations
//============================================================================

#include <include/io/dnaiosnx.hpp>
#include <include/functions/dnastringfuncs.hpp>
#include <include/functions/dnaiostreamfuncs.hpp>
#include <include/functions/dnatemplatefuncs.hpp>

namespace dynadjust { 
namespace iostreams {

void dna_io_snx::serialise_sinex(std::ofstream* snx_file, pvstn_t bstRecords,
					binary_file_meta_t& bst_meta, binary_file_meta_t& bms_meta,
					matrix_2d* estimates, matrix_2d* variances, const project_settings& p,
					UINT32& measurementParams, UINT32& unknownParams, double& sigmaZero,
					uint32_uint32_map* blockStationsMap, vUINT32* blockStations,
					const UINT32& blockCount, const UINT32& block,
					const CDnaDatum* datum)
{
	warningMessages_.clear();

	blockCount_ = blockCount;
	block_ = block;

	unknownParams_ = unknownParams;
	measurementParams_ = measurementParams;
	sigmaZero_ = sigmaZero;
	blockStationsMap_ = blockStationsMap;
	blockStations_ = blockStations;

	serialise_meta(snx_file, bst_meta, bms_meta, p, datum);
	serialise_site_id(snx_file, bstRecords);
	serialise_statistics(snx_file);
	serialise_solution_estimates(snx_file, bstRecords, estimates, variances, datum);
	serialise_solution_variances(snx_file, variances);

	*snx_file << "%ENDSNX" << std::endl;
}

void dna_io_snx::print_line(std::ofstream* snx_file)
{
	*snx_file << 
		"*-------------------------------------------------------------------------------" << 
	 	std::endl;
}
	

void dna_io_snx::serialise_meta(std::ofstream* snx_file, 
	binary_file_meta_t& bst_meta, binary_file_meta_t& bms_meta, 
	const project_settings& p, const CDnaDatum* datum)
{
	
	*snx_file << 
		"%=SNX 2.00 " <<				// SINEX version
		"DNA "; 						// Agency creating this file

	// Creation time of this SINEX file
	dateSINEXFormat(snx_file, boost::gregorian::day_clock::local_day(), true);

	// the agency providing the data in the SINEX file
	*snx_file << " DNA ";
	
	// adjustment epoch
	dateSINEXFormat(snx_file, datum->GetEpoch());		// start
	*snx_file << " ";
	dateSINEXFormat(snx_file, datum->GetEpoch());		// end = start
	*snx_file << " ";
	
	std::stringstream numberofparams;
	numberofparams << std::right << std::setw(5) << unknownParams_;
	std::string numberofparamsstr(numberofparams.str());
	numberofparamsstr = findandreplace<std::string>(numberofparamsstr, " ", "0");

	*snx_file <<
		"P " <<							// Technique(s) used to generate the SINEX solution
		std::left << std::setw(6) << 
			numberofparamsstr <<		// Number of parameters estimated in the SINEX file
		"0 " <<							// Single character indicating the constraint in the SINEX solution
		"S           " <<				// Solution types contained in this SINEX file.
										// S � all station parameters (i.e. station coordinates, station velocities, biases, geocenter)
										// O - Orbits
										// E - Earth Orientation Parameter
										// T � Troposphere
										// C � Celestial Reference Frame
										// A � Antenna parameters
										// BLANK
	 	std::endl;
	
	// FILE/REFERENCE
	print_line(snx_file);
	std::stringstream ss;

	// description: Organization(s) gathering/altering the file contents.
	*snx_file << "+FILE/REFERENCE" << std::endl <<
		"*INFO_TYPE_________ INFO________________________________________________________" << std::endl <<
		" DESCRIPTION        " << std::left << "Network " << p.g.network_name << std::endl;
	
	// output: Description of the file contents.
	ss.str("");
	if (blockCount_ > 1)
		ss << "Phased adjustment results. Block " << block_ + 1 << " of " << blockCount_;
	else
		ss << "Simultaneous adjustment results.";
	*snx_file << " OUTPUT             " << std::left << std::setw(60) << ss.str() << std::endl;

	// contact: address of the relevant contact email
	//*snx_file << " OUTPUT             " << std::left << std::setw(60) << __COPYRIGHT_OWNER__ << std::endl;

	// software/hardware
	*snx_file << snx_softwarehardware_text();

	// input files.  first, create a unique list of filenames
	std::vector<std::string> files;
	for (std::uint64_t i(0); i<bms_meta.inputFileCount; ++i)
		files.push_back(bms_meta.inputFileMeta[i].filename);
	for (std::uint64_t i(0); i<bst_meta.inputFileCount; ++i)
		files.push_back(bst_meta.inputFileMeta[i].filename);
	strip_duplicates(files);
	// second, print
	for_each(
		files.begin(), files.end(),
		[this, snx_file] (std::string& file) { 
			*snx_file << " INPUT              " << std::left << std::setw(60) << file << std::endl;	
	});
	
	*snx_file << "-FILE/REFERENCE" << std::endl;

	// FILE/COMMENTS
	print_line(snx_file);
	*snx_file << "+FILE/COMMENT" << std::endl;

	// print explanation of block results
	if (blockCount_ > 1)
	{
		*snx_file << 
			" This file contains the rigorous estimates for block " << block_ + 1 << 
			" of a segmented" << std::endl << 
			" network comprised of " << blockCount_ << " blocks. Due to the way in which junction stations" << std::endl << 
			" are carried through successive blocks, stations appearing in this file" << std::endl << 
			" may also be found in other SINEX files relating to this network, such as" << std::endl <<
			" " << p.g.network_name << "-block1.snx, " << p.g.network_name << "-block2.snx, etc." << std::endl;
	}
	*snx_file <<
		"-FILE/COMMENT" << std::endl;
}
	

void dna_io_snx::add_warning(const std::string& message, SINEX_WARN_TYPE warning)
{
	std::stringstream ss;

	switch (warning)
	{
	case excessive_name_chars:
		ss << "Station name " << message << " exceeds four characters.";
		break;
	default:
		ss << message;
	}
	
	warningMessages_.push_back(ss.str());
}
	

void dna_io_snx::print_warnings(std::ofstream* warning_file, const std::string& fileName)
{
	// Print formatted header
	print_file_header(*warning_file, "DYNADJUST SINEX OUTPUT WARNINGS FILE");
	*warning_file << std::setw(PRINT_VAR_PAD) << std::left << "File name:" << boost::filesystem::system_complete(fileName).string() << std::endl;

	*warning_file << OUTPUTLINE << std::endl << std::endl;

	for_each(
		warningMessages_.begin(), warningMessages_.end(),
		[warning_file] (std::string warning) {
			*warning_file << warning << std::endl;
	});
}
	

void dna_io_snx::serialise_site_id(std::ofstream* snx_file, pvstn_t bstRecords)
{

	print_line(snx_file);

	*snx_file << "+SITE/ID" << std::endl;
	*snx_file << "*CODE PT __DOMES__ T _STATION DESCRIPTION__ APPROX_LON_ APPROX_LAT_ _APP_H_" << std::endl;

	const station_t* stn;
	std::string stationName;

	for (UINT32 i(0); i<blockStationsMap_->size(); ++i)
	{
		stn = &(bstRecords->at(blockStations_->at(i)));
		stationName = stn->stationName;

		if (stationName.length() > 4)
			add_warning(stationName, excessive_name_chars);

		*snx_file << " " <<
			std::left  << std::setw(4) << stationName.substr(0, 4) << " " <<
			std::right << std::setw(2) << "A" << " " <<
			std::left << std::setw(9) << stationName.substr(0, 9) << " " <<
			std::right << std::setw(1) << "P" << " " <<
			std::left << std::setw(22) << std::string(stn->description).substr(0, 22) << " " <<
			std::right << std::setw(11) << FormatDmsString(RadtoDms(stn->currentLongitude), 5, true, false) << " " <<
			std::right << std::setw(11) << FormatDmsString(RadtoDms(stn->currentLatitude), 5, true, false) << " " <<
			std::right << std::setw(7) << std::setprecision(1) << std::fixed << stn->currentHeight << std::endl;	
	}

	*snx_file << "-SITE/ID" << std::endl;
}
	

void dna_io_snx::serialise_statistics(std::ofstream* snx_file)
{
	print_line(snx_file);

	*snx_file << "+SOLUTION/STATISTICS" << std::endl;
	*snx_file << "*_STATISTICAL PARAMETER________ __VALUE(S)____________" << std::endl;

	*snx_file << " " <<
		std::left << std::setw(30) << "NUMBER OF OBSERVATIONS" << " " << 
		std::right << std::setw(22) << measurementParams_ << std::endl;
	*snx_file << " " <<
		std::left << std::setw(30) << "NUMBER OF UNKNOWNS" << " " <<
		std::right << std::setw(22) << unknownParams_ << std::endl;
	*snx_file << " " <<
		std::left << std::setw(30) << "NUMBER OF DEGREES OF FREEDOM" << " " <<
		std::right << std::setw(22) << (measurementParams_ - unknownParams_) << std::endl;
	*snx_file << " " <<
		std::left << std::setw(30) << "VARIANCE FACTOR" << " " << 
		std::right << std::setw(22) << std::fixed << std::setprecision(6) << sigmaZero_ << std::endl;
	
	*snx_file << "-SOLUTION/STATISTICS" << std::endl;
}
	

void dna_io_snx::serialise_solution_estimates(std::ofstream* snx_file, pvstn_t bstRecords,
				matrix_2d* estimates, matrix_2d* variances, const CDnaDatum* datum)
{
	print_line(snx_file);

	*snx_file << "+SOLUTION/ESTIMATE" << std::endl;
	*snx_file << "*INDEX TYPE__ CODE PT SOLN _REF_EPOCH__ UNIT S __ESTIMATED VALUE____ _STD_DEV___" << std::endl;

	UINT32 i, j, index(1);
	std::string floating_value;
	std::stringstream ss;

	// Print stations
	for (i=0; i<blockStationsMap_->size(); ++i)
	{		
		j = (*blockStationsMap_)[blockStations_->at(i)] * 3;

		// estimated parameter (X)
		*snx_file << " " <<
			// parameter index
			std::right << std::setw(5) << index++ << " " << 
			// parameter type
			"STAX   " <<
			// 4 character site code
			std::left << std::setw(4) << std::string(bstRecords->at(blockStations_->at(i)).stationName).substr(0, 4) << " " <<
			// Point code
			std::right << std::setw(2) << "A" << " " <<
			// Solution id
			"0001 ";
		// epoch
		dateSINEXFormat(snx_file, datum->GetEpoch());

		*snx_file << " " <<
			// Parameter units
			std::left << std::setw(4) << "m" << " " <<
			// Constraint code
			"0" << " ";

		// Parameter estimate
		ss.str("");
		ss << std::setiosflags(std::ios_base::uppercase | std::ios_base::scientific) << std::setprecision(14) << estimates->get(j, 0);
		*snx_file << std::right << std::setw(21) << ss.str() << " ";
		// standard deviation
		ss.str("");
		ss << std::setprecision(5) << sqrt(variances->get(j, j));
		*snx_file << std::right << std::setw(11) << ss.str() << std::endl;

		// estimated parameter (Y)
		*snx_file << " " <<
			// parameter index
			std::right << std::setw(5) << index++ << " " << 
			// parameter type
			"STAY   " <<
			// 4 character site code
			std::left << std::setw(4) << std::string(bstRecords->at(blockStations_->at(i)).stationName).substr(0, 4) << " " <<
			// Point code
			std::right << std::setw(2) << "A" << " " <<
			// Solution id
			"0001 ";
		// epoch
		dateSINEXFormat(snx_file, datum->GetEpoch());

		*snx_file << " " <<
			// Parameter units
			std::left << std::setw(4) << "m" << " " <<
			// Constraint code
			"0" << " ";

		// Parameter estimate
		ss.str("");
		ss << std::setprecision(14) << estimates->get(j+1, 0);
		*snx_file << std::right << std::setw(21) << ss.str() << " ";
		// standard deviation
		ss.str("");
		ss << std::setprecision(5) << sqrt(variances->get(j+1, j+1));
		*snx_file << std::right << std::setw(11) << ss.str() << std::endl;

		// estimated parameter (Z)
		*snx_file << " " <<
			// parameter index
			std::right << std::setw(5) << index++ << " " << 
			// parameter type
			"STAZ   " <<
			// 4 character site code
			std::left << std::setw(4) << std::string(bstRecords->at(blockStations_->at(i)).stationName).substr(0, 4) << " " <<
			// Point code
			std::right << std::setw(2) << "A" << " " <<
			// Solution id
			"0001 ";
		// epoch
		dateSINEXFormat(snx_file, datum->GetEpoch());

		*snx_file << " " <<
			// Parameter units
			std::left << std::setw(4) << "m" << " " <<
			// Constraint code
			"0" << " ";

		// Parameter estimate
		ss.str("");
		ss << std::setprecision(14) << estimates->get(j+2, 0);
		*snx_file << std::right << std::setw(21) << ss.str() << " ";
		// standard deviation
		ss.str("");
		ss << std::setprecision(5) << sqrt(variances->get(j+2, j+2));
		*snx_file << std::right << std::setw(11) << ss.str() << std::endl;
	}

	*snx_file << "-SOLUTION/ESTIMATE" << std::endl;
}
	

void dna_io_snx::print_matrix_index(std::ofstream* snx_file, const UINT32& row, const UINT32& col)
{
	*snx_file << " " <<
		std::right << std::setw(5) << row + 1 << " " <<
		std::right << std::setw(5) << col + 1 << " ";
}

void dna_io_snx::serialise_solution_variances(std::ofstream* snx_file, matrix_2d* variances)
{
	print_line(snx_file);

	*snx_file << "+SOLUTION/MATRIX_ESTIMATE L COVA" << std::endl;
	*snx_file << "*PARA1 PARA2 ____PARA2+0__________ ____PARA2+1__________ ____PARA2+2__________" << std::endl;

	UINT32 row, col, max_dimension(variances->rows());
	std::string floating_value;
	std::stringstream ss;
	ss << std::setiosflags(std::ios_base::uppercase | std::ios_base::scientific);

	bool newRecord(true);
	UINT16 field(1);

	// Print stations
	for (row=0; row<max_dimension; ++row)
	{
		field = 1;
		for (col=0; col<row+1; ++col)
		{
			if (newRecord)
			{
				print_matrix_index(snx_file, row, col);
				newRecord = false;
			}

			// variance
			ss.str("");
			ss << std::setprecision(14) << variances->get(row, col);
			*snx_file << std::right << std::setw(21) << ss.str() << " ";
			
			if (row == col || ++field > 3)
			{
				*snx_file << std::endl;
				newRecord = true;
				field = 1;
			}
		}
	}

	*snx_file << "-SOLUTION/MATRIX_ESTIMATE L COVA" << std::endl;
}
	


} // dnaiostreams
} // dynadjust
