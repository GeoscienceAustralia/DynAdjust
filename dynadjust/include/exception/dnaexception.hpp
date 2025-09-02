//============================================================================
// Name         : dnaexception.hpp
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
// Description  : DynAdjust exception include file
//============================================================================

#ifndef DNAEXCEPTION_HPP
#define DNAEXCEPTION_HPP

#if defined(_MSC_VER)
	#if defined(LIST_INCLUDES_ON_BUILD) 
		#pragma message("  " __FILE__) 
	#endif
#endif

/// \cond
#include <exception>
#include <stdexcept>
/// \endcond

#include <include/config/dnatypes.hpp>

typedef enum _PARSE_STATUS_ {
	PARSE_SUCCESS = 0,
	PARSE_MSR_MAP = 1,
	PARSE_UNRECOGNISED_FILE = 2,
	PARSE_EXCEPTION_RAISED = 3
} PARSE_STATUS;

typedef enum _SEGMENT_STATUS_ {
	SEGMENT_SUCCESS = 0,
	SEGMENT_BLOCK_ERROR = 1,
	SEGMENT_EXCEPTION_RAISED = 2
} SEGMENT_STATUS;

// order of severity
typedef enum _ADJUST_STATUS_ {
	ADJUST_SUCCESS = 0,
	ADJUST_MAX_ITERATIONS_EXCEEDED = 1,
	ADJUST_THRESHOLD_EXCEEDED = 2,
	ADJUST_TEST_FAILED = 3,
	ADJUST_BLOCK_ERROR = 4,
	ADJUST_EXCEPTION_RAISED = 5,
	ADJUST_CANCELLED = 6
} ADJUST_STATUS;

typedef enum _PLOT_STATUS_ {
	PLOT_SUCCESS = 0,
	PLOT_BLOCK_ERROR = 1,
	PLOT_EXCEPTION_RAISED = 2
} PLOT_STATUS;

typedef enum _REFTRAN_STATUS_ {
	REFTRAN_SUCCESS = 0,
	REFTRAN_ERROR = 1,
	REFTRAN_EXCEPTION_RAISED = 2,
	REFTRAN_DIRECT_PARAMS_UNAVAILABLE = 3,
	REFTRAN_TRANS_ON_PLATE_REQUIRED = 4,
	REFTRAN_WGS84_TRANS_UNSUPPORTED = 5
} REFTRAN_STATUS;

namespace dynadjust { namespace exception {

// Exception handler
class XMLInteropException : public std::exception
{
public:
	explicit XMLInteropException(const std::string& what, const UINT32& line_no)
		: what_(what), line_no_(line_no) {}
	virtual const char* what() const throw() { return what_.c_str(); }
	virtual ~XMLInteropException() throw() {}

private:
	std::string what_;
	UINT32 line_no_;
};

// Exception handler
class NetPlotException : public std::exception
{
public:
	explicit NetPlotException(const std::string& what, const UINT32& line_no)
		: what_(what), line_no_(line_no) {}
	virtual const char* what() const throw() { return what_.c_str(); }
	virtual ~NetPlotException() throw() {}

private:
	std::string what_;
	UINT32 line_no_;
};

// Exception handler
class NetSegmentException : public std::exception
{
public:
	explicit NetSegmentException(const std::string& what, const UINT32& line_no)
		: what_(what), line_no_(line_no) {}
	virtual const char* what() const throw() { return what_.c_str(); }
	virtual ~NetSegmentException() throw() {}

private:
	std::string what_;
	UINT32 line_no_;
};

// Exception handler
class NetAdjustException : public std::exception
{
public:
	explicit NetAdjustException(const std::string& what, const UINT32 block_no)
		: what_(what), block_no_(block_no) {}
	virtual const char* what() const throw() { return what_.c_str(); }
	virtual ~NetAdjustException() throw() {}

private:
	std::string what_;
	UINT32 block_no_;
};

// Exception handler
class NetMemoryException : public std::exception
{
public:
	explicit NetMemoryException(const std::string& what)
		: what_(what) {}
	virtual const char* what() const throw() { return what_.c_str(); }
	virtual ~NetMemoryException() throw() {}

private:
	std::string what_;
};

// Exception handler
class NetGeoidException : public std::exception
{
public:
	explicit NetGeoidException(const std::string& what, const int& error_no, const UINT32& line_no = 0)
		: what_(what), error_no_(error_no), line_no_(line_no) {}
	virtual const char* what() const throw() { return what_.c_str(); }
	virtual ~NetGeoidException() throw() {}

   inline int error_no () const { return error_no_; }
   inline UINT32 line_no () const { return line_no_; }

private:
	std::string what_;
	int error_no_;
	UINT32 line_no_;
};

// Exception handler
class RefTranException : public std::exception
{
public:
	explicit RefTranException(const std::string& what, size_t exception_type = REFTRAN_EXCEPTION_RAISED)
		: what_(what)
		, exception_type_(exception_type) {}
	virtual const char* what() const throw() { return what_.c_str(); }
	virtual const size_t exception_type() const throw() { return exception_type_; }
	virtual ~RefTranException() throw() {}

private:
	std::string what_;
	size_t exception_type_;
};


}	// namespace exception
}	// namespace dynadjust


#endif  // DNAEXCEPTION_HPP
