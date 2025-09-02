//============================================================================
// Name         : dynadjust_file.cpp
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
// Description  : DynAdjust file io operations
//============================================================================

/// \cond
#include <boost/date_time/gregorian/gregorian.hpp>
/// \endcond

#include <include/config/dnaversion.hpp>
#include <include/io/dynadjust_file.hpp>
#include <include/functions/dnastrmanipfuncs.hpp>

namespace dynadjust {
namespace iostreams {

DynadjustFile::DynadjustFile(void)
	: version_(__FILE_VERSION__)
	, date_("")
	, app_name_("DNA" + std::string(__SHORT_VERSION__))
{
}

// copy constructors
DynadjustFile::DynadjustFile(const DynadjustFile& newFile)
	: version_(newFile.version_)
	, date_(newFile.date_)
	, app_name_(newFile.app_name_)
{
}
	

DynadjustFile::~DynadjustFile(void)
{
}
	

DynadjustFile& DynadjustFile::operator= (const DynadjustFile& rhs)
{
	// check for assignment to self!
	if (this == &rhs)
		return *this;

	version_ = rhs.version_;
	date_ = rhs.date_;
	app_name_ = rhs.app_name_;

	return *this;
}

void DynadjustFile::WriteFileInfo(std::ofstream& file_stream)
{
	WriteVersion(file_stream);
	WriteDate(file_stream);
	WriteApp(file_stream);
}
	

void DynadjustFile::ReadFileInfo(std::ifstream& file_stream)
{
	ReadVersion(file_stream);
	ReadDate(file_stream);
	ReadApp(file_stream);
}
	

void DynadjustFile::WriteFileMetadata(std::ofstream& file_stream, binary_file_meta_t& file_meta)
{
	// Write the metadata
	file_stream.write(reinterpret_cast<char *>(&file_meta.binCount), sizeof(std::uint64_t)); 
	file_stream.write(reinterpret_cast<char *>(&file_meta.reduced), sizeof(bool)); 
	file_stream.write(reinterpret_cast<char *>(file_meta.modifiedBy), MOD_NAME_WIDTH); 
	
	// Write the epsg code and epoch
	file_stream.write(reinterpret_cast<char *>(file_meta.epsgCode), STN_EPSG_WIDTH);
	file_stream.write(reinterpret_cast<char *>(file_meta.epoch), STN_EPOCH_WIDTH);

	file_stream.write(reinterpret_cast<char*>(&file_meta.reftran), sizeof(bool));
	file_stream.write(reinterpret_cast<char*>(&file_meta.geoid), sizeof(bool));

	// Write file count and file meta
	file_stream.write(reinterpret_cast<char *>(&file_meta.inputFileCount), sizeof(std::uint64_t)); 
	for (std::uint64_t i(0); i<file_meta.inputFileCount; ++i)
	{
		file_stream.write(reinterpret_cast<char *>(file_meta.inputFileMeta[i].filename), FILE_NAME_WIDTH); 
		file_stream.write(reinterpret_cast<char *>(file_meta.inputFileMeta[i].epsgCode), STN_EPSG_WIDTH);
		file_stream.write(reinterpret_cast<char *>(file_meta.inputFileMeta[i].epoch), STN_EPOCH_WIDTH);
		file_stream.write(reinterpret_cast<char *>(&file_meta.inputFileMeta[i].filetype), sizeof(UINT16)); 
		file_stream.write(reinterpret_cast<char *>(&file_meta.inputFileMeta[i].datatype), sizeof(UINT16)); 
	}
}
	

void DynadjustFile::ReadFileMetadata(std::ifstream& file_stream, binary_file_meta_t& file_meta)
{
	// Read the metadata
	file_stream.read(reinterpret_cast<char *>(&file_meta.binCount), sizeof(std::uint64_t)); 
	file_stream.read(reinterpret_cast<char *>(&file_meta.reduced), sizeof(bool)); 
	file_stream.read(reinterpret_cast<char *>(file_meta.modifiedBy), MOD_NAME_WIDTH);

	// Read the epsg code and epoch
	file_stream.read(reinterpret_cast<char *>(file_meta.epsgCode), STN_EPSG_WIDTH);
	file_stream.read(reinterpret_cast<char *>(file_meta.epoch), STN_EPOCH_WIDTH);

	file_stream.read(reinterpret_cast<char*>(&file_meta.reftran), sizeof(bool));
	file_stream.read(reinterpret_cast<char*>(&file_meta.geoid), sizeof(bool));

	// Read file count and file meta
	file_stream.read(reinterpret_cast<char *>(&file_meta.inputFileCount), sizeof(std::uint64_t));
	if (file_meta.inputFileMeta != NULL)
		delete []file_meta.inputFileMeta;

	file_meta.inputFileMeta = new input_file_meta_t[file_meta.inputFileCount];
		
	for (std::uint64_t i(0); i<file_meta.inputFileCount; ++i)
	{
		file_stream.read(reinterpret_cast<char *>(file_meta.inputFileMeta[i].filename), FILE_NAME_WIDTH); 
		file_stream.read(reinterpret_cast<char *>(file_meta.inputFileMeta[i].epsgCode), STN_EPSG_WIDTH);
		file_stream.read(reinterpret_cast<char *>(file_meta.inputFileMeta[i].epoch), STN_EPOCH_WIDTH);
		file_stream.read(reinterpret_cast<char *>(&file_meta.inputFileMeta[i].filetype), sizeof(UINT16)); 
		file_stream.read(reinterpret_cast<char *>(&file_meta.inputFileMeta[i].datatype), sizeof(UINT16)); 
	}
}
	

void DynadjustFile::WriteVersion(std::ofstream& file_stream)
{
	char versionField[identifier_field_width+1];
	
	// write version field name
	versionField[identifier_field_width] = '\0';
	file_stream.write(const_cast<char *>(version_header), identifier_field_width); 
	
	// write version, with safeguard against overwriting by substr
	snprintf(versionField, sizeof(versionField), "%*s", identifier_field_width, version_.substr(0, identifier_field_width).c_str());
	file_stream.write(reinterpret_cast<char *>(versionField), identifier_field_width); 
		
}
	

void DynadjustFile::ReadVersion(std::ifstream& file_stream)
{
	char versionField[identifier_field_width+1];

	// read version field name
	versionField[identifier_field_width] = '\0';
	file_stream.read(reinterpret_cast<char *>(versionField), identifier_field_width);
	
	// read version
	file_stream.read(reinterpret_cast<char *>(versionField), identifier_field_width); 
	version_ = versionField;
	version_ = trimstr(version_);
}
	

void DynadjustFile::WriteDate(std::ofstream& file_stream)
{
	char dateField[identifier_field_width+1];

	// Form date string
	boost::gregorian::date today(boost::gregorian::day_clock::local_day());
	std::stringstream date_string;
	date_string << std::right << to_iso_extended_string(today);
	date_ = date_string.str();
	
	// write creation date field name
	dateField[identifier_field_width] = '\0';
	file_stream.write(const_cast<char *>(create_date_header), identifier_field_width); 
	
	// write creation date, with safeguard against overwriting by substr
	snprintf(dateField, sizeof(dateField), "%*s", identifier_field_width, date_.substr(0, identifier_field_width).c_str());
	file_stream.write(reinterpret_cast<char *>(dateField), identifier_field_width); 
		
}
	

void DynadjustFile::ReadDate(std::ifstream& file_stream)
{
	char dateField[identifier_field_width+1];

	// read creation date field name
	dateField[identifier_field_width] = '\0';
	file_stream.read(reinterpret_cast<char *>(dateField), identifier_field_width);
	
	// read creation date
	file_stream.read(reinterpret_cast<char *>(dateField), identifier_field_width); 
	date_ = dateField;
	date_ = trimstr(date_);
}
	

void DynadjustFile::WriteApp(std::ofstream& file_stream)
{
	char appField[identifier_field_width+1];
	
	// write application field name
	appField[identifier_field_width] = '\0';
	file_stream.write(const_cast<char *>(create_by_header), identifier_field_width); 
	
	// write application, with safeguard against overwriting by substr
	snprintf(appField, sizeof(appField), "%*s", identifier_field_width, app_name_.substr(0, identifier_field_width).c_str());
	file_stream.write(reinterpret_cast<char *>(appField), identifier_field_width); 
		
}
	

void DynadjustFile::ReadApp(std::ifstream& file_stream)
{
	char appField[identifier_field_width+1];

	// read application field name
	appField[identifier_field_width] = '\0';
	file_stream.read(reinterpret_cast<char *>(appField), identifier_field_width);
	
	// read application
	file_stream.read(reinterpret_cast<char *>(appField), identifier_field_width); 
	app_name_ = appField;
	app_name_ = trimstr(app_name_);
}

}	// namespace measurements
}	// namespace dynadjust
