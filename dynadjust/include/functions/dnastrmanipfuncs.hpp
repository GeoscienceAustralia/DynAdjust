//============================================================================
// Name         : dnastrmanipfuncs.hpp
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
// Description  : String Manipulation Functions
//============================================================================

#ifndef DNASTRMANIPFUNCS_H_
#define DNASTRMANIPFUNCS_H_

#if defined(_MSC_VER)
	#if defined(LIST_INCLUDES_ON_BUILD) 
		#pragma message("  " __FILE__) 
	#endif
#endif

/// \cond
#include <stdio.h>
#include <stdarg.h>
#include <string>
#include <algorithm>
#include <functional>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <iomanip>
#include <math.h>
#include <iostream>
#include <charconv>
#include <stdexcept>

#include <boost/tokenizer.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/phoenix/core.hpp>
#include <boost/phoenix/operator.hpp>
/// \endcond

#include <include/functions/dnastrutils.hpp>
#include <include/functions/dnatemplatecalcfuncs.hpp>

using boost::spirit::qi::double_;
using boost::spirit::qi::float_;
using boost::spirit::qi::parse;

// Helper template to convert between types (replacement for boost::lexical_cast)
template <typename T, typename U>
T lexical_cast(const U& value)
{
	std::stringstream ss;
	ss << value;
	T result;
	if (!(ss >> result) || !ss.eof()) {
		throw std::runtime_error("lexical_cast failed");
	}
	return result;
}

// Specialization for string to string (no-op)
template <>
inline std::string lexical_cast<std::string, std::string>(const std::string& value)
{
	return value;
}

template <class T, class U>
T val_uint(const U& value)
{
	return lexical_cast<T, U>(value);
}

template <class T, class U>
T valorno_uint(const U& value, T& var)
{
	if (iequals(value, "no"))
		return var = 0;

	var = 1;
	return lexical_cast<T, U>(value);
}

template <class T>
T valorno_string(const T& value)
{
	if (iequals(value, "no"))
		return "";
	return value;
}

template <class T, class U>
T yesno_uint(const U& value)
{
	if (iequals(value, "yes"))
		return 1;
	//if (iequals(value, "no"))
		return 0;

}

template <class T>
std::string yesno_string(const T value)
{
	if (value == 1)
		return std::string("yes");
	else
		return std::string("no");
}

template <class T>
std::string bool_string(const T value)
{
	if (value)
		return std::string("true");
	else
		return std::string("false");
}

template <class T>
T findandreplace(const T& source, const T& find, const T& replace)
{
	size_t index;
	T str = source;
	while ((index = str.find(find)) != T::npos)
		str.replace(index, find.length(), replace);
	return str;
}

// convert to upper case
template <class T>
void str_toupper(std::string& strtext)
{
	transform(strtext.begin(), strtext.end(), strtext.begin(), (T(*)(T))toupper);
}

// convert to upper case
template <class S, class T>
S str_upper(const T* const strtext)
{
	S strtext_(strtext);
	str_toupper<int>(strtext_);
	return strtext_;
}

// convert to upper case
template <class S>
S str_upper(S& strtext)
{
	str_toupper<int>(strtext);
	return strtext;
}

template <class T>
T trimstr_(const T& Src, const T& c)
{
	size_t p2 = Src.find_last_not_of(c);
	if (p2 == std::string::npos)
		return std::string();
	size_t p1 = Src.find_first_not_of(c);
	if (p1 == std::string::npos)
		p1 = 0;
	return Src.substr(p1, (p2-p1)+1);
}

template <class T>
T trimstrleft_(const T& Src, const T& c)
{
	size_t p2 = Src.length();
	size_t p1 = Src.find_first_not_of(c);
	if (p1 == std::string::npos)
		p1 = 0;
	return Src.substr(p1, (p2-p1)+1);
}

template <class T>
T trimstrright_(const T& Src, const T& c)
{
	size_t p2 = Src.find_last_not_of(c);
	if (p2 == std::string::npos)
		return std::string();
	return Src.substr(0, p2+1);
}

template <class T>
T trimstr(T& Src) {
	return trimstr_(Src, static_cast<T>(" \r\n"));
}
template <class T>
T trimstr(const T& Src) {
	return trimstr_(Src, static_cast<T>(" \r\n"));
}
template <class T>
T trimstrleft(const T& Src) {
	return trimstrleft_(Src, static_cast<T>(" \r\n"));
}
template <class T>
T trimstrright(const T& Src) {
	return trimstrright_(Src, static_cast<T>(" \r\n"));
}

template <class T, class U>
std::string StringFromTW(const T& t, const U& width, const U& precision=0)
{
	std::stringstream ss;

	// Assume number at precision prints within width
	ss << std::setw(width) << std::fixed << std::right << std::setprecision(precision) << t;
	int trim = (int)trimstr(ss.str()).length() - width;

	// Formatted string length is less than or equal to the fixed width
	if (trim <= 0)
		return ss.str();

	ss.str("");

	// Formatted string length exceeds fixed width by trim.
	// Use scientific notation if the overflow is 
	// greater than the precision
	if (trim > 0)
	{
		// scientific notation requires 6 or 5
		// characters (e.g. -1e+00 or 1e+00), not
		// including the fractional part
		if (width < UINT32(t < 0. ? 6 : 5))
		{
			// insufficient space for scientific notation!!!
			ss << std::setw(width) << std::string(width, '#');
			return ss.str();
		}
		
		int prec1 = width - (t < 0. ? 6 : 5);
		// if precision is required, an extra space is needed for
		// the decimal point, so reduce by one.
		if (prec1 > 0)
			prec1--;
		
		int prec = std::min<int>(precision, prec1);
		
		ss << std::setw(width) << std::scientific << std::right << std::setprecision(prec) << t;
	}
	else
	{
		int prec = precision - trim;
		if (prec < 0)
			prec = 0;
		ss << std::setw(width) << std::fixed << std::right << std::setprecision(prec) << t;
	}

	return ss.str();
}

template <class T>
std::string StringFromT(const T& t, const int& precision=-1)
{
	std::stringstream ss;
	if (precision < 0)
		ss << t;
	else
		ss << std::setprecision(precision) << std::fixed << t;
	return ss.str();
}

template <class T, class U, class S>
S double_string_width(const T value, const U& width)
{
	S s;
	char* cvalue = new char[width+1];
#if defined(_WIN32) || defined(__WIN32__)
	_snprintf(cvalue, width, "%.2f", value);
#else
	snprintf(cvalue, width, "%.2f", value);
#endif
	s = cvalue;
	delete[] cvalue;
	return trimstr(s);
}




// See here for boost::lexical_cast vs sscanf performance 
// http://www.boost.org/doc/libs/1_48_0/doc/html/boost_lexical_cast/performance.html
//
// An interesting article on conversion speeds
// http://tinodidriksen.com/2011/05/28/cpp-convert-string-to-double-speed/

template <class T>
T DoubleFromString(const std::string& str)
{
	T t(0.);
	parse(str.begin(), str.end(), double_, t);
	return t;
}

template <class T>
T FloatFromString(const std::string& str)
{
	T t(0.);
	parse(str.begin(), str.end(), float_, t);
	return t;
}

template <class T>
void DoubleFromString(T& t, const std::string& str)
{
	parse(str.begin(), str.end(), double_, t);
}

template <class T>
void FloatFromString(T& t, const std::string& str)
{
	parse(str.begin(), str.end(), float_, t);
}

template <class T>
bool DoubleFromString_ZeroCheck(T& t, const std::string& str)
{
	parse(str.begin(), str.end(), double_, t);
	if (fabs(t) < PRECISION_1E35)
		return true;
	return false;
}

template <class T>
// signed
typename std::enable_if<std::is_signed<T>::value, T>::type LongFromString(const std::string& str)
{
	char* end(NULL);
	T t(strtol(str.c_str(), &end,  10));
	
	if (*end)
	{
		std::stringstream ss;
		ss << "String to long conversion error: non-convertible part: " << end;
		throw std::runtime_error(ss.str());
	}

	return t;
}

template <class T>
// unsigned
typename std::enable_if<std::is_unsigned<T>::value, T>::type LongFromString(const std::string& str)
{
	char* end(NULL);
	T t(strtoul(str.c_str(), &end,  10));
	
	if (*end)
	{
		std::stringstream ss;
		ss << "String to long conversion error: non-convertible part: " << end;
		throw std::runtime_error(ss.str());
	}

	return t;
}

//template <class T>
//T LongFromString(const std::string& str)
//{
//	char* end(NULL);
//	T t;
//	std::string s(typeid(t).name());
//
//	if (iequals(s, "unsigned int") || iequals(s, "unsigned long") ||	// msvc
//	    iequals(s, "j") || iequals(s, "m"))								// gcc
//		t = strtoul(str.c_str(), &end,  10);
//	else if (iequals(s, "int") || iequals(s, "long") ||					// msvc
//	    iequals(s, "i") || iequals(s, "l"))								// gcc
//		t = strtol(str.c_str(), &end,  10);
//	else
//	{
//		std::stringstream ss;
//		ss << "String to long conversion error: " << s << " is not an integer." << end;
//		throw std::runtime_error(ss.str());
//	}
//
//	if (*end)
//	{
//		std::stringstream ss;
//		ss << "String to long conversion error: non-convertible part: " << end;
//		throw std::runtime_error(ss.str());
//	}
//
//	return t;
//}

//template <class T>
//void LongFromString(T *l, const std::string& str)
//{
//	//*l = atol(str.c_str());
//
//	*l = LongFromString<T>(str);
//}
//
template <class T>
void FromDmsString(T *d, const std::string& str)
{
	DmstoDeg(atof(str.c_str()), d);
}

template <class T>
T FromDmsString(const std::string& str)
{
	return DmstoDeg(atof(str.c_str()));
}

template <class T>
void RadFromDmsString(T *d, const std::string& str)
{
	DmstoDeg(atof(str.c_str()), d);
	Radians(d);
}

template <class T>
void RadFromSecondsString(T *d, const std::string& str)
{
	*d = SecondstoRadians(atof(str.c_str()));
}

template <class T>
void SplitDelimitedString(const T& str, const T& separator, std::vector<T>* tokenList)
{
	// use boost to extract stations from comma delimited string
	typedef boost::tokenizer<boost::char_separator<char> > btokenizer;
	boost::char_separator<char> sepa(separator.c_str());
	btokenizer tokens(str, sepa);
	tokenList->clear();
	try {
		for_each(tokens.begin(), tokens.end(), 
			[&tokenList](const std::string& s) {
				if (!trimstr(s).empty())
					tokenList->push_back(trimstr(s));
		});
	}
	catch (...) {
		return;
	}	
}

#endif //DNASTRMANIPFUNCS_H_
