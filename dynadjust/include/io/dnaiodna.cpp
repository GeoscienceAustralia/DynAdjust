//============================================================================
// Name         : dnaiodna.cpp
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
// Description  : DNA file io helper
//============================================================================

#include <include/io/dnaiodna.hpp>
#include <include/parameters/dnaepsg.hpp>
#include <include/measurement_types/dnameasurement_types.hpp>
#include <include/functions/dnastrutils.hpp>

using namespace dynadjust::measurements;
using namespace dynadjust::epsg;

namespace dynadjust { 
namespace iostreams {

void dna_io_dna::write_dna_files(vdnaStnPtr* vStations, vdnaMsrPtr* vMeasurements, 
			const std::string& stnfilename, const std::string& msrfilename, const std::string& networkname,
			const CDnaDatum& datum, const CDnaProjection& projection, bool flagUnused,
			const std::string& stn_comment, const std::string& msr_comment)
{

	write_stn_file(vStations, stnfilename, networkname, datum, projection, flagUnused, stn_comment);
	write_msr_file(vMeasurements, msrfilename, networkname, datum, msr_comment);
}

void dna_io_dna::write_dna_files(pvstn_t vbinary_stn, pvmsr_t vbinary_msr, 
			const std::string& stnfilename, const std::string& msrfilename, const std::string& networkname,
			const CDnaDatum& datum, const CDnaProjection& projection, bool flagUnused,
			const std::string& stn_comment, const std::string& msr_comment)
{
	write_stn_file(vbinary_stn, stnfilename, networkname, datum, projection, flagUnused, stn_comment);
	write_msr_file(*vbinary_stn, vbinary_msr, msrfilename, networkname, datum, msr_comment);
}
	
void dna_io_dna::create_file_stn(std::ofstream* ptr, const std::string& filename)
{
	determineDNASTNFieldParameters<UINT16>("3.01", dsl_, dsw_);
	create_file_pointer(ptr, filename);
}
	
void dna_io_dna::create_file_msr(std::ofstream* ptr, const std::string& filename)
{
	determineDNAMSRFieldParameters<UINT16>("3.01", dml_, dmw_);
	create_file_pointer(ptr, filename);
}

void dna_io_dna::prepare_sort_list(const UINT32 count)
{
	vStationList_.resize(count);

	// initialise vector with 0,1,2,...,n-2,n-1,n
	initialiseIncrementingIntegerVector<UINT32>(vStationList_, count);
}
	
void dna_io_dna::create_file_pointer(std::ofstream* ptr, const std::string& filename)
{
	try {
		// Create file pointer to DNA file. 
		file_opener(*ptr, filename);
	}
	catch (const std::runtime_error& e) {
		throw std::runtime_error(e.what());
	}
}
	
void dna_io_dna::open_file_pointer(std::ifstream* ptr, const std::string& filename)
{
	try {
		// Open DNA file.
		file_opener(*ptr, filename, std::ios::in, ascii, true);
	}
	catch (const std::runtime_error& e) {
		throw std::runtime_error(e.what());
	}
}

void dna_io_dna::write_stn_header_data(std::ofstream* ptr, const std::string& networkname, const std::string& datum,
	const std::string& epoch, const size_t& count, const std::string& comment)
{
	// Write version line
	dna_header(*ptr, "3.01", "STN", datum, epoch, count);

	// Write comments
	dna_comment(*ptr, "File type:    Station file");
	dna_comment(*ptr, "Project name: " + networkname);

	// Write header comment line, e.g. "Coordinates exported from reftran."
	dna_comment(*ptr, comment);
}

void dna_io_dna::write_stn_header(std::ofstream* ptr, vdnaStnPtr* vStations, const std::string& networkname,
			const CDnaDatum& datum, bool flagUnused, const std::string& comment)
{
	// print stations
	// Has the user specified --flag-unused-stations, in which case, do not
	// print stations marked unused?
	size_t count(0);
	if (flagUnused) 
	{
		for_each(vStations->begin(), vStations->end(),
			[&count](const dnaStnPtr& stn) {
				if (stn.get()->IsNotUnused())
					count++;
		});
	}
	else
		count = vStations->size();

	// write header
	write_stn_header_data(ptr, networkname, datum.GetName(), datum.GetEpoch_s(), count, comment);	
}
	

void dna_io_dna::write_stn_header(std::ofstream* ptr, pvstn_t vbinary_stn, const std::string& networkname,
			const CDnaDatum& datum, bool flagUnused, const std::string& comment)
{
	// print stations
	// Has the user specified --flag-unused-stations, in which case, do not
	// print stations marked unused?
	size_t count(0);
	if (flagUnused) 
	{
		for_each(vbinary_stn->begin(), vbinary_stn->end(),
			[&count](const station_t& stn) {
				if (stn.unusedStation == 0)
					count++;
		});
	}
	else
		count = vbinary_stn->size();

	// write header
	write_stn_header_data(ptr, networkname, datum.GetName(), datum.GetEpoch_s(), count, comment);
}
	

void dna_io_dna::read_ren_file(const std::string& filename, pv_string_vstring_pair stnRenaming)
{
	std::ifstream renaming_file;

	// create file pointer
	open_file_pointer(&renaming_file, filename);

	std::string version;
	INPUT_DATA_TYPE idt;
	CDnaDatum datum;
	std::string geoidVersion, fileEpsg, fileEpoch;
	UINT32 count(0);

	// read header information
	read_dna_header(&renaming_file, version, idt, datum, fileEpsg, fileEpoch, geoidVersion, count);

	read_ren_data(&renaming_file, stnRenaming);

	renaming_file.close();
}

void dna_io_dna::read_ren_data(std::ifstream* ptr, pv_string_vstring_pair stnRenaming)
{
	std::string sBuf;

	string_vstring_pair stnNames;
	vstring altStnNames;

	while (!ptr->eof())			// while EOF not found
	{
		try {
			getline((*ptr), sBuf);
		}
		catch (...) {
			if (ptr->eof())
				return;
			std::stringstream ss;
			ss << "read_ren_data(): Could not read from the renaming file." << std::endl;
			std::runtime_error(ss.str());
		}

		// blank or whitespace?
		if (trimstr(sBuf).empty())			
			continue;

		// Ignore lines with blank station name
		if (trimstr(sBuf.substr(dsl_.stn_name, dsw_.stn_name)).empty())			
			continue;
		
		// Ignore lines with comments
		if (sBuf.compare(0, 1, "*") == 0)
			continue;
		
		std::string tmp;
		
		// Add the preferred name
		try {
			tmp = trimstr(sBuf.substr(dsl_.stn_name, dsw_.stn_name));
			stnNames.first = tmp;
		}
		catch (...) {
			std::stringstream ss;
			ss << "read_ren_data(): Could not extract station name from the record:  " << 
			 	std::endl << "    " << sBuf << std::endl;
			std::runtime_error(ss.str());
		}		

		// Alternative names
		altStnNames.clear();
		UINT32 positionEndName(dsw_.stn_name+dsw_.stn_name);
		UINT32 positionNextName(dsw_.stn_name);
		while (sBuf.length() > positionNextName)
		{
			if (sBuf.length() > positionEndName)
				tmp = trimstr(sBuf.substr(positionNextName, dsw_.stn_name));
			else
				tmp = trimstr(sBuf.substr(positionNextName));

			if (!tmp.empty())
				altStnNames.push_back(tmp);

			positionNextName = positionEndName;
			positionEndName += dsw_.stn_name;
		}
		
		// sort for faster searching
		std::sort(altStnNames.begin(), altStnNames.end());

		// Add the alternate names
		stnNames.second = altStnNames;

		// Okay, add the pair to the list
		stnRenaming->push_back(stnNames);
	}
}
	

void dna_io_dna::read_dna_header(std::ifstream* ptr, std::string& version, INPUT_DATA_TYPE& idt, 
	CDnaDatum& referenceframe,
	std::string& fileEpsg, std::string& fileEpoch, std::string& geoidversion, UINT32& count)
{
	std::string sBuf;
	getline((*ptr), sBuf);
	sBuf = trimstr(sBuf);
	m_filespecifiedReferenceFrame_ = false;
	m_filespecifiedEpoch_ = false;

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

	// Attempt to get the file's type
	try {
		determineDNASTNFieldParameters<UINT16>(version, dsl_, dsw_);
		determineDNAMSRFieldParameters<UINT16>(version, dml_, dmw_);
	}
	catch (const std::runtime_error& e) {
		std::stringstream ssError;
		ssError << "- Error: Unable to determine DNA file version." << std::endl <<
			sBuf << std::endl << " " << e.what() << std::endl;
		throw std::runtime_error(ssError.str());
	}

	// Version 1
	if (iequals(version, "1.00"))
	{
		idt = unknown;							// could be stn or msr
		count = 100;							// set default 100 stations
		return;
	}

	std::string type;
	// Attempt to get the file's type
	try {
		type = trimstr(sBuf.substr(12, 3));
	}
	catch (const std::runtime_error& e) {
		std::stringstream ssError;
		ssError << "- Error: File type has not been provided in the header" << std::endl <<
			sBuf << std::endl << e.what() << std::endl;
		throw std::runtime_error(ssError.str());
	}

	// Station file
	if (iequals(type, "stn"))
		idt = stn_data;
	// Measurement file
	else if (iequals(type, "msr"))
		idt = msr_data;
	// Geoid information file
	else if (iequals(type, "geo"))
		idt = geo_data;
	// Station renaming
	else if (iequals(type, "ren"))
		idt = ren_data;
	else
	{
		idt = unknown;
		std::stringstream ssError;
		ssError << "The supplied filetype '" << type << "' is not recognised" << std::endl;
		throw std::runtime_error(ssError.str());
	}

	if (sBuf.length() < 29)
	{
		count = 100;
		return;
	}

	std::string file_referenceframe;
	// Attempt to get the default reference frame
	try {
		if (sBuf.length() > 42)
			file_referenceframe = trimstr(sBuf.substr(29, 14));
		else if (sBuf.length() > 29)
			file_referenceframe = trimstr(sBuf.substr(29));
		else
			file_referenceframe = "";
	}
	catch (...) {
		// datum not provided?
		file_referenceframe = "";
	}

	if (!file_referenceframe.empty())
		m_filespecifiedReferenceFrame_ = true;

	std::string epoch_version;
	// Attempt to get the default epoch / geoid version
	try {
		if (sBuf.length() > 56)
			epoch_version = trimstr(sBuf.substr(43, 14));
		else if (sBuf.length() > 43)
			epoch_version = trimstr(sBuf.substr(43));
		else
			epoch_version = "";
	}
	catch (...) {
		// epoch / version not provided?
		epoch_version = "";
	}

	if (!epoch_version.empty())
		m_filespecifiedEpoch_ = true;

	// Try to get the reference frame
	// For the first file, and unless the user has specified a reference frame,
	// the reference frame will be used to set the frame for all other files.
	try {

		switch (idt)
		{
		case stn_data:
		case msr_data:
			
			// Capture file epsg
			if (file_referenceframe.empty())
				// Empty? Get the epsg code from the project datum
				fileEpsg = referenceframe.GetEpsgCode_s();
			else
				fileEpsg = epsgStringFromName<std::string>(file_referenceframe);

			// capture file epoch
			if (epoch_version.empty())
			{
				// Empty?  Get the epoch of the nominated epsg (set above)
				fileEpoch = referenceepochFromEpsgString<std::string>(fileEpsg);
			}
			else
				fileEpoch = epoch_version;

			break;
		case geo_data:
			geoidversion = epoch_version;
			break;
		default:
			break;
		}
	}
	catch (const std::runtime_error& e) {
		std::stringstream ssError;
		ssError << "The supplied frame (" << file_referenceframe << ") is not recognised" << std::endl <<
			e.what() << std::endl;
		throw std::runtime_error(ssError.str());
	}	

	// Attempt to get the data count
	try {
		if (sBuf.length() > 66)
			count = val_uint<UINT32, std::string>(trimstr(sBuf.substr(57, 10)));
		else if (sBuf.length() > 57)
			count = val_uint<UINT32, std::string>(trimstr(sBuf.substr(57)));
		else 
			count = 100;
	}
	catch (...) {
		count = 100;
	}
}


void dna_io_dna::write_stn_file(pvstn_t vbinary_stn, const std::string& stnfilename, const std::string& networkname,
				const CDnaDatum& datum, const CDnaProjection& projection, bool flagUnused,
				const std::string& comment)
{
	// Print DNA stations from a vector of stn_t
	std::ofstream dna_stn_file;

	try {
		// Create file pointer
		create_file_stn(&dna_stn_file, stnfilename);

		// Write header and comment
		write_stn_header(&dna_stn_file, vbinary_stn, networkname, datum, flagUnused, comment);

		// Sort on original file order
		prepare_sort_list(static_cast<UINT32>(vbinary_stn->size()));
		CompareStnFileOrder<station_t, UINT32> stnorderCompareFunc(vbinary_stn);
		std::sort(vStationList_.begin(), vStationList_.end(), stnorderCompareFunc);

		dnaStnPtr stnPtr(new CDnaStation(datum.GetName(), datum.GetEpoch_s()));

		// print stations
		// Has the user specified --flag-unused-stations, in wich case, do not
		// print stations marked unused?
		if (flagUnused) 
		{
			// Print stations in DNA format
			for_each(vStationList_.begin(), vStationList_.end(),
				[&dna_stn_file, &stnPtr, &projection, &datum, &vbinary_stn, this](const UINT32& stn) {
					if (stnPtr->IsNotUnused())
					{
						stnPtr->SetStationRec(vbinary_stn->at(stn));
						stnPtr->WriteDNAXMLStnCurrentEstimates(&dna_stn_file, datum.GetEllipsoidRef(), 
							const_cast<CDnaProjection*>(&projection), dna, &dsw_);
					}
			});
		}
		else
		{
			// Print stations in DNA format
			for_each(vStationList_.begin(), vStationList_.end(),
				[&dna_stn_file, &stnPtr, &projection, &datum, &vbinary_stn, this](const UINT32& stn) {
					stnPtr->SetStationRec(vbinary_stn->at(stn));
					stnPtr->WriteDNAXMLStnCurrentEstimates(&dna_stn_file, datum.GetEllipsoidRef(), 
						const_cast<CDnaProjection*>(&projection), dna, &dsw_);
			});
		}

		dna_stn_file.close();
		
	}
	catch (const std::ifstream::failure& f) {
		throw std::runtime_error(f.what());
	}
}
	

void dna_io_dna::write_stn_file(vdnaStnPtr* vStations, const std::string& stnfilename, const std::string& networkname,
			const CDnaDatum& datum, const CDnaProjection& projection, bool flagUnused,
			const std::string& comment)
{
	// Print DNA stations from a vector of dnaStnPtr
	std::ofstream dna_stn_file;

	try {
		// Create file pointer
		create_file_stn(&dna_stn_file, stnfilename);

		// Write header and comment
		write_stn_header(&dna_stn_file, vStations, networkname, datum, flagUnused, comment);

		// Sort on original file order
		std::sort(vStations->begin(), vStations->end(), CompareStnFileOrder_CDnaStn<CDnaStation>());

		// print stations
		// Has the user specified --flag-unused-stations, in wich case, do not
		// print stations marked unused?
		if (flagUnused) 
		{
			// Print stations in DNA format
			for_each(vStations->begin(), vStations->end(),
				[&dna_stn_file, &projection, &datum, this](const dnaStnPtr& stn) {
					if (stn.get()->IsNotUnused())
						stn.get()->WriteDNAXMLStnCurrentEstimates(&dna_stn_file, datum.GetEllipsoidRef(), 
							const_cast<CDnaProjection*>(&projection), dna, &dsw_);
			});
		}
		else
		{
			// Print stations in DNA format
			for_each(vStations->begin(), vStations->end(),
				[&dna_stn_file, &projection, &datum, this](const dnaStnPtr& stn) {
					stn.get()->WriteDNAXMLStnCurrentEstimates(&dna_stn_file, datum.GetEllipsoidRef(), 
						const_cast<CDnaProjection*>(&projection), dna, &dsw_);
			});
		}

		dna_stn_file.close();
		
	}
	catch (const std::ifstream::failure& f) {
		throw std::runtime_error(f.what());
	}	
}

void dna_io_dna::write_msr_header_data(std::ofstream* ptr, const std::string& networkname, const std::string& datum,
	const std::string& epoch, const size_t& count, const std::string& comment)
{
	// Write version line
	dna_header(*ptr, "3.01", "MSR", datum, epoch, count);

	// Write comments
	dna_comment(*ptr, "File type:    Measurement file");
	dna_comment(*ptr, "Project name: " + networkname);

	// Write header comment line
	dna_comment(*ptr, comment);

}


void dna_io_dna::write_msr_header(std::ofstream* ptr, vdnaMsrPtr* vMeasurements, const std::string& networkname,
	const CDnaDatum& datum, const std::string& comment)
{
	// write header
	write_msr_header_data(ptr, networkname, datum.GetName(), datum.GetEpoch_s(), vMeasurements->size(), comment);
}
	

void dna_io_dna::write_msr_header(std::ofstream* ptr, pvmsr_t vbinary_msrn, const std::string& networkname,
	const CDnaDatum& datum, const std::string& comment)
{
	// write header
	write_msr_header_data(ptr, networkname, datum.GetName(), datum.GetEpoch_s(), vbinary_msrn->size(), comment);
}


void dna_io_dna::write_msr_file(const vstn_t& vbinary_stn, pvmsr_t vbinary_msr, const std::string& msrfilename, const std::string& networkname,
	const CDnaDatum& datum, const std::string& comment)
{
	// Print DNA measurements	

	std::ofstream dna_msr_file;

	it_vmsr_t _it_msr(vbinary_msr->begin());
	dnaMsrPtr msrPtr;
	size_t dbindex;
	it_vdbid_t _it_dbid;

	try {
		// Create file pointer
		create_file_msr(&dna_msr_file, msrfilename);

		// Write header and comment
		write_msr_header(&dna_msr_file, vbinary_msr, networkname, datum, comment);

		// print measurements
		for (_it_msr=vbinary_msr->begin(); _it_msr!=vbinary_msr->end(); ++_it_msr)
		{
			ResetMeasurementPtr<char>(&msrPtr, _it_msr->measType);

			// Database IDs
			if (m_databaseIDsSet_)
			{
				dbindex = std::distance(vbinary_msr->begin(), _it_msr);
				_it_dbid = pv_msr_db_map_->begin() + dbindex;
			}

			msrPtr->SetMeasurementRec(vbinary_stn, _it_msr, _it_dbid);
			msrPtr->ResolveSourceFile(sourceFileMeta_, sourceFileCount_);
			msrPtr->WriteDNAMsr(&dna_msr_file, dmw_, dml_);
		}

		dna_msr_file.close();

	}
	catch (const std::ifstream::failure& f) {
		throw std::runtime_error(f.what());
	}
}
	
void dna_io_dna::write_msr_file(vdnaMsrPtr* vMeasurements, const std::string& msrfilename, const std::string& networkname,
	const CDnaDatum& datum,	const std::string& comment)
{
	// Print DNA measurements	
	
	std::ofstream dna_msr_file;
	
	_it_vdnamsrptr _it_msr(vMeasurements->begin());
	size_t dbindex;
	it_vdbid_t _it_dbid;
	
	try {
		// Create file pointer
		create_file_msr(&dna_msr_file, msrfilename);

		// Write header and comment
		write_msr_header(&dna_msr_file, vMeasurements, networkname, datum, comment);

		// print measurements
		for (_it_msr = vMeasurements->begin(); _it_msr != vMeasurements->end(); _it_msr++)
		{
			// Get database IDs
			dbindex = std::distance(vMeasurements->begin(), _it_msr);
			_it_dbid = pv_msr_db_map_->begin() + dbindex;
			_it_msr->get()->SetDatabaseMaps(_it_dbid);

			// Write
			_it_msr->get()->WriteDNAMsr(&dna_msr_file, dmw_, dml_);
		}
		
		dna_msr_file.close();
		
	}
	catch (const std::ifstream::failure& f) {
		throw std::runtime_error(f.what());
	}
	

}

void dna_io_dna::set_dbid_ptr(pv_msr_database_id_map pv_msr_db_map)
{
	pv_msr_db_map_ = pv_msr_db_map;
	if (pv_msr_db_map->size() > 0)
		m_databaseIDsSet_ = true;
	else
		m_databaseIDsSet_ = false;
}


} // dnaiostreams
} // dynadjust

