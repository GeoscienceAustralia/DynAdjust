//============================================================================
// Name         : dnaiomap.cpp
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
// Description  : DynAdjust station map file io operations
//============================================================================

#include <fstream>
/// \endcond

#include <include/io/dnaiomap.hpp>
#include <include/functions/dnaiostreamfuncs.hpp>
#include <include/functions/dnastrmanipfuncs.hpp>

namespace dynadjust { 
namespace iostreams {

void dna_io_map::load_map_file(const std::string& map_filename, pv_string_uint32_pair stnsMap) 
{	
	std::ifstream map_file;
	std::stringstream ss;
	ss << "load_map_file(): An error was encountered when opening " << map_filename << "." << std::endl;

	try {
		// open stations map file.  Throws runtime_error on failure.
		file_opener(map_file, map_filename, std::ios::in | std::ios::binary, binary, true);
	}
	catch (const std::runtime_error& e) {
		ss << e.what();
		throw boost::enable_current_exception(std::runtime_error(ss.str()));
	}
	catch (...) {
		throw boost::enable_current_exception(std::runtime_error(ss.str()));
	}
	
	ss.str("");
	ss << "load_map_file(): An error was encountered when reading from " << map_filename << "." << std::endl;

	stnsMap->clear();
	UINT32 mapsize;
	string_uint32_pair stnID;
	char stationName[STN_NAME_WIDTH];	// 56 characters
	
	try {
		// read the file information
		ReadFileInfo(map_file);
		
		// read the number of records
		map_file.read(reinterpret_cast<char *>(&mapsize), sizeof(UINT32));
		stnsMap->reserve(mapsize);
		
		// read the records
		for (UINT32 i=0; i<mapsize; i++)
		{
			map_file.read(const_cast<char *>(stationName), sizeof(stationName));
			stnID.first = stationName;
			map_file.read(reinterpret_cast<char *>(&stnID.second), sizeof(UINT32));
			stnsMap->push_back(stnID);
		}
	}
	catch (const std::ios_base::failure& f) {
		ss << f.what();
		throw boost::enable_current_exception(std::runtime_error(ss.str()));
	}
	catch (const std::runtime_error& e) {
		ss << e.what();
		throw boost::enable_current_exception(std::runtime_error(ss.str()));
	}
	catch (...) {
		throw boost::enable_current_exception(std::runtime_error(ss.str()));
	}

	map_file.close();

	if (stnsMap->size() < mapsize)
		throw boost::enable_current_exception(
			std::runtime_error("load_map_file(): Could not allocate sufficient memory for the Station map."));
}

void dna_io_map::write_map_file(const std::string& map_filename, pv_string_uint32_pair stnsMap)
{
	std::ofstream map_file;
	std::stringstream ss;
	ss << "write_map_file(): An error was encountered when opening " << map_filename << "." << std::endl;

	try {
		// Open station map file.  Throws runtime_error on failure.
		file_opener(map_file, map_filename, 
			std::ios::out | std::ios::binary, binary);
	}
	catch (const std::runtime_error& e) {
		ss << e.what();
		throw boost::enable_current_exception(std::runtime_error(ss.str()));
	}
	catch (...) {
		throw boost::enable_current_exception(std::runtime_error(ss.str()));
	}

	ss.str("");
	ss << "write_map_file(): An error was encountered when writing to " << map_filename << "." << std::endl;

	v_string_uint32_pair::const_iterator _it_stnmap(stnsMap->begin());
	UINT32 mapval(static_cast<UINT32>(stnsMap->size()));
	try {
		// write the file information
		WriteFileInfo(map_file);
		
		// write the data
		map_file.write(reinterpret_cast<char *>(&mapval), sizeof(UINT32));
		char stationName[STN_NAME_WIDTH];
		memset(stationName, '\0', sizeof(stationName));
		for (; _it_stnmap!=stnsMap->end(); ++_it_stnmap)
		{
			strcpy(stationName, _it_stnmap->first.c_str());
			mapval = _it_stnmap->second;
			map_file.write(stationName, sizeof(stationName));
			map_file.write(reinterpret_cast<char *>(&mapval), sizeof(UINT32));
		}
	}
	catch (const std::ios_base::failure& f) {
		ss << f.what();
		throw boost::enable_current_exception(std::runtime_error(ss.str()));
	}
	catch (const std::runtime_error& e) {
		ss << e.what();
		throw boost::enable_current_exception(std::runtime_error(ss.str()));
	}
	catch (...) {
		throw boost::enable_current_exception(std::runtime_error(ss.str()));
	}

	map_file.close();
}
	

void dna_io_map::write_map_file_txt(const std::string& map_filename, pv_string_uint32_pair stnsMap)
{
	std::ofstream map_file;
	std::stringstream ss;
	ss << "write_map_file_txt(): An error was encountered when opening " << map_filename << "." << std::endl;

	try {
		// Create station map text file.  Throws runtime_error on failure.
		file_opener(map_file, map_filename);
	}
	catch (const std::runtime_error& e) {
		ss << e.what();
		throw boost::enable_current_exception(std::runtime_error(ss.str()));
	}
	catch (...) {
		throw boost::enable_current_exception(std::runtime_error(ss.str()));
	}

	ss.str("");
	ss << "write_map_file_txt(): An error was encountered when writing to " << map_filename << "." << std::endl;

	try {
		std::stringstream ss_map;
		ss_map << stnsMap->size() << " stations";
		map_file << std::left << std::setw(STATION) << ss_map.str();
		map_file << std::setw(HEADER_20) << std::right << "Stn. index" << std::endl;
		v_string_uint32_pair::const_iterator _it_stnmap(stnsMap->begin());
		for (_it_stnmap=stnsMap->begin(); _it_stnmap!=stnsMap->end(); ++_it_stnmap)
		{
			map_file << std::setw(STATION) << std::left << 
				_it_stnmap->first.c_str() <<
				std::setw(HEADER_20) << std::right << _it_stnmap->second << std::endl;
		}
	}
	catch (const std::ios_base::failure& f) {
		ss << f.what();
		throw boost::enable_current_exception(std::runtime_error(ss.str()));
	}
	catch (const std::runtime_error& e) {
		ss << e.what();
		throw boost::enable_current_exception(std::runtime_error(ss.str()));
	}
	catch (...) {
		throw boost::enable_current_exception(std::runtime_error(ss.str()));
	}

	map_file.close();
}

//void dna_io_map::load_renaming_file(const std::string& renaming_filename, pv_string_string_pair stnRenaming)
//{
//	std::ifstream renaming_file;
//
//	std::stringstream ss;
//	ss << "load_renaming_file(): An error was encountered when opening " << renaming_filename << "." << std::endl;
//
//	// The contents of the renaming file is as follows:
//	//
//	//  --------------------------------------------------------------------------------
//	//  DYNADJUST STATION RENAMING FILE
//	//  
//	//  Station count                      3 
//	//  Station name width                 20
//	//  --------------------------------------------------------------------------------
//	//  
//	//  STATION NAMES
//	//  
//	//  OLD NAME             NEW NAME
//	//  -------------------- --------------------
//	//  GABO                 262600060
//	//  BEEC                 209901750
//	//  PTLD                 341404410
//
//	//
//	try {
//		// Load renaming file.  Throws runtime_error on failure.
//		file_opener(renaming_file, renaming_filename, 
//			std::ios::in, ascii, true);
//	}
//	catch (const std::runtime_error& e) {
//		ss << e.what();
//		throw boost::enable_current_exception(std::runtime_error(ss.str()));
//	}
//	catch (...) {
//		throw boost::enable_current_exception(std::runtime_error(ss.str()));
//	}
//	
//	ss.str("");
//	ss << "load_renaming_file(): An error was encountered when reading from " << renaming_filename << "." << std::endl;
//
//	stnRenaming->clear();
//
//	string_string_pair stnNames;
//	std::string sBuf("");
//	UINT32 stationWidth(STATION), stationCount(100);
//	UINT32 line(0);
//	
//	try {
//		
//		// continue until "STATION NAMES" is found
//		do 
//		{
//			line++;
//			getline(renaming_file, sBuf);
//
//			if (boost::iequals(trimstr(sBuf), "STATION NAMES"))
//				break;
//
//			// blank or whitespace?
//			if (trimstr(sBuf).empty())		
//				continue;
//
//			if (boost::iequals(trimstr(sBuf.substr(0, 13)), "Station count"))
//			{
//				stationCount = boost::lexical_cast<UINT16, std::string>(trimstr(sBuf.substr(PRINT_VAR_PAD)));
//				continue;
//			}
//
//			if (boost::iequals(trimstr(sBuf.substr(0, 18)), "Station name width"))
//			{
//				stationWidth = boost::lexical_cast<UINT16, std::string>(trimstr(sBuf.substr(PRINT_VAR_PAD)));
//				continue;
//			}
//
//		}
//		while (!boost::iequals(trimstr(sBuf), "STATION NAMES"));
//		
//		stnRenaming->reserve(stationCount);
//		
//		// Okay, now get the data
//		while (!renaming_file.eof())
//		{
//			getline(renaming_file, sBuf);
//
//			// blank or whitespace?
//			if (trimstr(sBuf).empty())			
//				continue;
//
//			if (trimstr(sBuf).length() < stationWidth)
//				continue;
//
//			if (boost::iequals(trimstr(sBuf.substr(0, 8)), "OLD NAME"))
//				continue;
//
//			if (boost::iequals(trimstr(sBuf.substr(0, 3)), "---"))
//				continue;
//
//			// Ignore lines with blank station name
//			if (trimstr(sBuf.substr(0, stationWidth)).empty())			
//				continue;
//
//			// Ignore lines with blank substitute name
//			if (trimstr(sBuf.substr(stationWidth+1)).empty())			
//				continue;
//
//			// initialise
//			stnNames.first = "";
//			stnNames.second = "";
//
//			// get the names
//			stnNames.first = trimstr(sBuf.substr(0, stationWidth));	
//			stnNames.second = trimstr(sBuf.substr(stationWidth+1));
//			stnRenaming->push_back(stnNames);
//		}
//	}
//	catch (const std::ios_base::failure& f) {
//		if (renaming_file.eof())
//		{
//			renaming_file.close();
//			return;
//		}
//		ss << f.what();
//		throw boost::enable_current_exception(std::runtime_error(ss.str()));
//	}
//	catch (const std::runtime_error& e) {
//		if (renaming_file.eof())
//		{
//			renaming_file.close();
//			return;
//		}
//		ss << e.what();
//		throw boost::enable_current_exception(std::runtime_error(ss.str()));
//	}
//	catch (...) {
//		if (renaming_file.eof())
//		{
//			renaming_file.close();
//			return;
//		}
//		throw boost::enable_current_exception(std::runtime_error(ss.str()));
//	}
//
//	renaming_file.close();
//}


} // dnaiostreams
} // dynadjust
