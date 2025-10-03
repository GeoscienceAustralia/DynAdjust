//============================================================================
// Name         : dnaiotbu.cpp
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
// Description  : Type B uncertainty file io and helps
//============================================================================

#include <include/io/dnaiotbu.hpp>
#include <include/functions/dnastrutils.hpp>

namespace dynadjust {
namespace iostreams {

void dna_io_tbu::assign_typeb_values_global(const vstring& typeBUncertainties, type_b_uncertainty& type_b)
{
	double e(0.), n(0.), u(0.);
	UINT32 i(0);
	std::stringstream ss;

	if (typeBUncertainties.size() > 2)
	{
		// 3 dimensions (east, north, up)
		if (!typeBUncertainties.at(i).empty())
			e = DoubleFromString<double>(typeBUncertainties.at(i));
		if (!typeBUncertainties.at(++i).empty())
			n = DoubleFromString<double>(typeBUncertainties.at(i));
		if (!typeBUncertainties.at(++i).empty())
			u = DoubleFromString<double>(typeBUncertainties.at(i));		
	}
	else if (typeBUncertainties.size() > 1) 
	{
		// 2 dimensions (east, north)
		if (!typeBUncertainties.at(i).empty())
			e = DoubleFromString<double>(typeBUncertainties.at(i));
		if (!typeBUncertainties.at(++i).empty())
			n = DoubleFromString<double>(typeBUncertainties.at(i));
	}	
	else if (typeBUncertainties.size() > 0) 
	{
		// 1 dimension (up)
		if (!typeBUncertainties.at(i).empty())
			u = DoubleFromString<double>(typeBUncertainties.at(i));
	}
	else {
		ss << "  No Type b uncertainties provided." << std::endl;
		throw std::runtime_error(ss.str());
	}

	// assign
	type_b.set_typeb_values_3d(e*e, n*n, u*u);
}
	

void dna_io_tbu::assign_typeb_values_local(const vstring& typeBUncertainties, type_b_uncertainty& type_b)
{
	double e(0.), n(0.), u(0.);
	UINT32 i(0);

	// 3 dimensions (east, north, up)
	if (!typeBUncertainties.at(i).empty())
		e = DoubleFromString<double>(typeBUncertainties.at(i));
	if (!typeBUncertainties.at(++i).empty())
		n = DoubleFromString<double>(typeBUncertainties.at(i));
	if (!typeBUncertainties.at(++i).empty())
		u = DoubleFromString<double>(typeBUncertainties.at(i));
	
	// assign
	type_b.set_typeb_values_3d(e*e, n*n, u*u);
}
	

void dna_io_tbu::validate_typeb_values(const std::string& argument, vstring& typeBUncertainties)
{
	std::stringstream ss;

	if (typeBUncertainties.size() == 0)
	{
		ss << "  No Type b uncertainties provided:" << std::endl <<
			"    " << argument << std::endl;
		throw std::runtime_error(ss.str());
	}

	UINT32 i(0), j;

	for (j = i; j < typeBUncertainties.size(); j++)
	{
		if (!is_floatingpoint<std::string>(typeBUncertainties.at(i)))
		{
			ss << "  Type b uncertainty '" << typeBUncertainties.at(i) <<
				"' is not a number:" << std::endl <<
				"    " << argument << std::endl;
			throw std::runtime_error(ss.str());
		}
		
		i++;
	}
}

void dna_io_tbu::read_tbu_header(std::ifstream* ptr, std::string& version, INPUT_DATA_TYPE& idt)
{
	std::string sBuf;
	getline((*ptr), sBuf);
	sBuf = trimstr(sBuf);

	// Set the default version
	version = "1.00";

	// Attempt to get the file's version
	try {
		if (iequals("!#=DNA", sBuf.substr(0, 6)))
			version = trimstr(sBuf.substr(6, 6));
	}
	catch (const std::runtime_error& e) {
		throw std::runtime_error(e.what());
	}

	std::string type;
	std::stringstream ssError;
	ssError << "  File type has not been provided in the header:" << std::endl <<
		"  " << sBuf << std::endl;

	if (sBuf.length() < 15)
		throw std::runtime_error(ssError.str());

	// Attempt to get the file's type
	try {
		type = trimstr(sBuf.substr(12, 3));
	}
	catch (const std::runtime_error& e) {
		ssError << "  " << e.what() << std::endl;
		throw std::runtime_error(ssError.str());
	}

	// Check this is a Type B file
	if (iequals(type, "tbu"))
		idt = tbu_data;
	else
	{
		idt = unknown;
		std::stringstream ssError;
		ssError << "  The supplied filetype '" << type << "' is not recognised:" << std::endl <<
			"  " << sBuf << std::endl;
		throw std::runtime_error(ssError.str());
	}
}

void dna_io_tbu::load_tbu_argument(const std::string& argument, type_b_uncertainty& type_b_uncertainties)
{
	vstring typeBUncertainties;
	typeBUncertainties.resize(3);

	// Extract constraints from comma delimited string
	SplitDelimitedString<std::string>(
		argument,		// the comma delimited string
		std::string(","),							// the delimiter
		&typeBUncertainties);					// the respective values

	// validate values. Throws on exception
	validate_typeb_values(argument, typeBUncertainties);

	// Assign the type b uncertainties to be applied to all stations, except
	// when site-specific type b uncertainties have been loaded from a file
	assign_typeb_values_global(typeBUncertainties, type_b_uncertainties);
}
	

void dna_io_tbu::identify_station_id(const std::string& stn_str, UINT32& stn_id, v_string_uint32_pair& vStnsMap)
{
	it_pair_string_vUINT32 it_stnmap_range;
	std::stringstream ss;

	// find this station in the station map
	it_stnmap_range = equal_range(vStnsMap.begin(), vStnsMap.end(), stn_str, StationNameIDCompareName());
	if (it_stnmap_range.first == it_stnmap_range.second)
	{
		ss << "  Station '" << stn_str <<
			"' is not included in the network." << std::endl;
		throw std::runtime_error(ss.str());
	}

	stn_id = it_stnmap_range.first->second;

}
	

void dna_io_tbu::load_tbu_file(const std::string& tbu_filename, v_type_b_uncertainty& type_b_uncertainties, v_string_uint32_pair& vStnsMap)
{
	// Type B uncertainty file structure is as follows.
	// Note - uncertainties are in metres, in the local reference frame (e,n,u)
	// and are given at 1 sigma (68%). 
	//
	//   station id (20 chars)  east uncertainty (13) north uncertainty (13) up uncertainty (13)
	//      ...
	//   EOF

	std::ifstream tbu_file;
	std::stringstream ss;
	ss << "load_tbu_file(): An error was encountered when opening " << tbu_filename << "." << std::endl;
	
	INPUT_DATA_TYPE idt;
	std::string version;

	try {
		// open ascii plate boundaries file.  Throws runtime_error on failure.
		file_opener(tbu_file, tbu_filename, std::ios::in, ascii, true);

		// read header information
		read_tbu_header(&tbu_file, version, idt);

	}
	catch (const std::runtime_error& e) {
		ss << e.what();
		throw std::runtime_error(ss.str());
	}
	catch (...) {
		throw std::runtime_error(ss.str());
	}

	ss.str("");
	ss << "load_tbu_file(): An error was encountered when reading from " << tbu_filename << "." << std::endl;

	UINT32 stn_id;
	std::string record, stn_str;
	vstring typeb_values;
	const UINT16 STDDEVn(STATION + STDDEV);
	const UINT16 STDDEVu(STATION + STDDEV + STDDEV);

	type_b_uncertainty typeBUncertainties;

	type_b_uncertainties.clear();

	try {

		typeb_values.resize(3);

		while (!tbu_file.eof())			// while EOF not found
		{
			typeb_values.at(0) = "";
			typeb_values.at(1) = "";
			typeb_values.at(2) = "";

			// get the plate identifier
			getline(tbu_file, record);

			// blank or whitespace?
			if (trimstr(record).empty())
				continue;

			// Ignore lines with comments
			if ((record.compare(0, 1, "*") == 0))
				continue;

			// obtain station name and unique ID
			stn_str = trimstr(record.substr(0, STATION));
			try {
				identify_station_id(stn_str, stn_id, vStnsMap);
			}
			catch (const std::runtime_error&) {
				// Can't find the station, ignore.
				continue;
			}

			typeBUncertainties.set_station_id(stn_id);
			
			// east uncertainty (may be blank)
			if (record.length() > STATION)
			{
				if (record.length() > STDDEVn)
					typeb_values.at(0) = trimstr(record.substr(STATION, STDDEV));
				else
					typeb_values.at(0) = trimstr(record.substr(STATION));
			}
			else
				// no uncertainties!
				continue;

			// north uncertainty (may be blank)
			if (record.length() > STDDEVn)
			{
				if (record.length() > STDDEVu)
					typeb_values.at(1) = trimstr(record.substr(STDDEVn, STDDEV));
				else
					typeb_values.at(1) = trimstr(record.substr(STDDEVn));
			}
			
			// up uncertainty (may be blank)
			if (record.length() > STDDEVu)
				typeb_values.at(2) = trimstr(record.substr(STDDEVu));

			try {
				// validate values. Throws on exception.  Print error and continue 
				// to the next set of values
				validate_typeb_values(record, typeb_values);
			}
			catch (const std::runtime_error&) {
				continue;
			}

			// Assign the site-specific type b uncertainties
			assign_typeb_values_local(typeb_values, typeBUncertainties);

			// add to the list
			type_b_uncertainties.push_back(typeBUncertainties);

		}

		tbu_file.close();

	}
	catch (const std::ios_base::failure& f) {
		if (tbu_file.eof())
		{
			tbu_file.close();
			return;
		}
		ss << f.what();
		throw std::runtime_error(ss.str());
	}
	catch (const std::runtime_error& e) {
		ss << e.what();
		throw std::runtime_error(ss.str());
	}
	catch (...) {
		if (tbu_file.eof())
		{
			tbu_file.close();
			return;
		}
		throw std::runtime_error(ss.str());
	}	

	return;
}

void dna_io_tbu::reduce_uncertainties_global(type_b_uncertainty& type_b, matrix_2d& var_cart, stn_t& bstBinaryRecord)
{

	matrix_2d var_local(type_b.type_b);
	PropagateVariances_LocalCart<double>(var_local, var_cart,
		bstBinaryRecord.currentLatitude, bstBinaryRecord.currentLongitude, true);
}
	

void dna_io_tbu::reduce_uncertainties_local(v_type_b_uncertainty& type_b, vstn_t& bstBinaryRecords)
{
	it_vstn_t stn_it(bstBinaryRecords.begin());
	it_type_b_uncertainty _it_tbu;
	matrix_2d var_local;

	for (_it_tbu = type_b.begin(); _it_tbu != type_b.end(); ++_it_tbu)
	{
		var_local = _it_tbu->type_b;
		stn_it = bstBinaryRecords.begin() + _it_tbu->station_id;

		PropagateVariances_LocalCart<double>(var_local, _it_tbu->type_b,
			stn_it->currentLatitude, stn_it->currentLongitude, true);
	}
}



} // dnaiostreams
} // dynadjust

