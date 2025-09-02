//============================================================================
// Name         : dnaiostreamfuncs.hpp
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
// Description  : Common file stream Functions
//============================================================================

#ifndef DNAIOSTREAMFUNCS_HPP_
#define DNAIOSTREAMFUNCS_HPP_

#if defined(_MSC_VER)
	#if defined(LIST_INCLUDES_ON_BUILD) 
		#pragma message("  " __FILE__) 
	#endif
#endif

/// \cond
#include <iostream>
#include <string>
#include <chrono>
#include <iomanip>
#include <ctime>
#include <filesystem>
#include <boost/date_time/posix_time/posix_time.hpp>

using std::ostringstream;
using std::locale;
/// \endcond

#include <include/config/dnaversion.hpp>
#include <include/config/dnaversion-stream.hpp>
#include <include/config/dnaconsts-iostream.hpp>

template <typename T>
T real_line_length_ascii(const T& line_length_ascii)
{
	// return the real length of a line based on the system definition 
	// of a "new line".  On Windows, a new line is comprised of two characters
	// namely, carriage return (CR, '\r') and line feed (LF, '\n'). On Unix
	// and Unix-like systems (Linux, mac OS, BSD, etc.), a new line is 
	// comprised of a line feed (LF, '\n') only.
	//

#if defined(__linux) || defined(sun) || defined(__unix__) || defined(__APPLE__)
	// "real" line length = line_length_ascii + "\n"
	return line_length_ascii + 1;

#else // #if defined(_WIN32) || defined(__WIN32__)
	// "real" line length = line_length_ascii + "\r\n"
	return line_length_ascii + 2;
	
#endif
}


template <typename T>
void file_opener(
	T* stream,
	const std::string& str, 
	std::ios_base::openmode mode=std::ios_base::out,	// default: output
	const iosMode type=ascii,				// default: ascii
	bool fileMustExist=false)				// default: no need for file to exist
{
	file_opener(*stream, str, mode, type, fileMustExist);
}

template <typename T>
void file_opener(
	T& stream,
	const std::string& str, 
	std::ios_base::openmode mode=std::ios_base::out,	// default: output
	const iosMode type=ascii,				// default: ascii
	bool fileMustExist=false)				// default: no need for file to exist
{
	try {
		stream.open(str.c_str(), mode);
		stream.exceptions (/*std::ios_base::eofbit | */std::ios_base::badbit | std::ios_base::failbit);
		stream.iword(0) = type;
	}
	catch (const std::ios_base::failure& f) {
		std::stringstream ss;
		if (fileMustExist && !std::filesystem::exists(str.c_str()))
			ss << "file_opener(): Can't find " << str << ".";
		else
			ss << "file_opener(): An error was encountered when opening " << 
				str << ". \n  Check that the file is not already opened.";
		ss << std::endl << f.what();
		throw std::runtime_error(ss.str());
	}

	if (!stream.good()) {
		std::stringstream ss;
		ss << "file_opener(): An error was encountered when opening " << 
			str << ". \n  Check that the file is not already opened.";
		throw std::runtime_error(ss.str());
	}
}

template <typename T>
void print_file_header(
	T& stream,
	const std::string& header)
{
	// Print formatted header
	stream << OUTPUTLINE << std::endl;
	stream << header << std::endl << std::endl;

	// version
	output_version(stream, true);
	stream << std::endl;

	// build
	output_build(stream, true);
	stream << std::endl;

	// File creation time
	stream << std::setw(PRINT_VAR_PAD) << std::left << "File created:";
	std::ostringstream datetime_ss;
	boost::posix_time::time_facet* p_time_output = new boost::posix_time::time_facet;
	std::locale special_locale (std::locale(""), p_time_output);
	// special_locale takes ownership of the p_time_output facet
	datetime_ss.imbue (special_locale);
	(*p_time_output).format("%A, %d %B %Y, %X");
	datetime_ss << boost::posix_time::second_clock::local_time();
	stream << datetime_ss.str().c_str() << std::endl;
}


#endif // DNAIOSTREAMFUNCS_HPP_
