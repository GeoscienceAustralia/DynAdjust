//============================================================================
// Name         : dynadjust_file.hpp
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

#ifndef DYNADJUST_FILE_H_
#define DYNADJUST_FILE_H_

#if defined(_MSC_VER)
	#if defined(LIST_INCLUDES_ON_BUILD) 
		#pragma message("  " __FILE__) 
	#endif
#endif

/// \cond
#include <fstream>
#include <string>
/// \endcond

#include <include/config/dnaexports.hpp>
#include <include/config/dnatypes-fwd.hpp>
#include <include/config/dnatypes-structs.hpp>

#define __FILE_VERSION__ "1.1"

namespace dynadjust {
namespace iostreams {


// File Info data: 
// Version identifier (10 character string), e.g. "VERSION   "
// Version value (10 character string), e.g. "       0.1"
// Created identifier (10 character string), e.g. "CREATED ON"
// Created value (10 character string), e.g. "2012-11-27"
// Application identifier (10 character string), e.g. "CREATED BY"
// Application value (10 character string), e.g. "    GEOMAP" or " DYNADJUST"
//
// For example, the beginning of the binary stream would have:
//  "VERSION          0.1CREATED ON2012-11-27CREATED BY DYNADJUST"
//  OR
//  "VERSION          0.1CREATED ON2013-03-01CREATED BY   THE_APP"

const UINT16 identifier_field_width(10);
const char* const version_header     = "VERSION   ";	// 10 characters
const char* const create_date_header = "CREATED ON";	// 10 characters
const char* const create_by_header   = "CREATED BY";	// 10 characters

class DNATYPE_API DynadjustFile {
public:
	DynadjustFile(void);
	DynadjustFile(const DynadjustFile&);
	virtual ~DynadjustFile(void);

	DynadjustFile& operator=(const DynadjustFile& rhs);

	inline std::string GetVersion() const { return version_; }
	inline void SetVersion(const std::string& version) { version_ = version; }

	void WriteFileInfo(std::ofstream& file_stream);
	void ReadFileInfo(std::ifstream& file_stream);

	void WriteFileMetadata(std::ofstream& file_stream, binary_file_meta_t& file_meta);
	void ReadFileMetadata(std::ifstream& file_stream, binary_file_meta_t& file_meta);

protected:

	bool versionAtLeast(int major, int minor) const;

	void WriteVersion(std::ofstream& file_stream);
	void ReadVersion(std::ifstream& file_stream);
	void WriteDate(std::ofstream& file_stream);
	void ReadDate(std::ifstream& file_stream);
	void WriteApp(std::ofstream& file_stream);
	void ReadApp(std::ifstream& file_stream);

	std::string	version_;
	std::string	date_;
	std::string	app_name_;
};

}	// namespace measurements
}	// namespace dynadjust

#endif  // DYNADJUST_FILE_H_
