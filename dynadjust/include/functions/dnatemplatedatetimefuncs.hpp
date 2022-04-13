//============================================================================
// Name         : dnatemplatedatetimefuncs.hpp
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
// Description  : Advanced and common calculation functions involving date and
//				  time
//============================================================================

#ifndef DNATEMPLATEDATETIMEFUNCS_H_
#define DNATEMPLATEDATETIMEFUNCS_H_

#if defined(_MSC_VER)
	#if defined(LIST_INCLUDES_ON_BUILD) 
		#pragma message("  " __FILE__) 
	#endif
#endif

#include <boost/timer/timer.hpp>
#include <boost/date_time/local_time/local_time.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <include/functions/dnastrmanipfuncs.hpp>

const UINT32 TIME_IMMEMORIAL = 1900;

using namespace std;
using namespace boost::timer;
using namespace boost::posix_time;
using namespace boost::gregorian;

// Boost implementation of gregorian calendar supports dates in the range 1400-Jan-01 to 9999-Dec-31.
// Integer representation of these dates is 2232400 and 5373484, which appears to be the Julian Date.
// GPS time begins on 6 January 1980 (2444245)

//date::date_int_type dit(2444245);
////date d(dit);
//date gpsday(1980, 1, 6);
//date today(day_clock::local_day());
//UINT32 id = gpsday.day_number();
//UINT32 itoday = today.day_number();

//stringstream str;
//const std::locale fmt(locale::classic(), new date_facet("%m/%d/%Y"));
//str.imbue(fmt);
//str << gpsday << endl;
//str << "Today:   " << today.day_number() << endl;
//str <<
//	static_cast<short>(today.year()) << endl <<
//	static_cast<short>(today.month()) << endl <<
//	static_cast<short>(today.day()) << endl << endl;
//str << "GPS day: " << gpsday.day_number() << endl;
//str <<
//	static_cast<short>(gpsday.year()) << endl <<
//	static_cast<short>(gpsday.month()) << endl <<
//	static_cast<short>(gpsday.day()) << endl;

//cout << str.str().c_str() << endl;

//TRACE("%s", str.str().c_str());

//days duration_since_gps(today- gpsday);
//
//str.str("");
//str << duration_since_gps;
//TRACE("Today is %d days since the beginning of GPS time\n", str.str().c_str());

template <class T>
void dateSINEXFormat(T* os, const date& theDate, bool calculateSeconds=false)
{
	// Time stamps must be in the following format:  YY:DOY:SECOD 
	// where:
	//	 YY     : year
	//	 DOY    : day of year
	//	 SECOD  : sec of day (e.g. the epoch)
	// Example:
	//		95:120:86399 denotes April 30, 1995 (23:59:59UT).
	//

	stringstream year;
	time_facet* facet(new time_facet("%y"));

	year.imbue(locale(year.getloc(), facet));
	ptime dateTime(theDate);
	year << dateTime;

	stringstream dayofyear;
	dayofyear << right << setw(3) << theDate.day_of_year();
	string dayofyearstr(dayofyear.str());
	dayofyearstr = findandreplace<string>(dayofyearstr, " ", "0");

	*os <<
		right << setw(2) <<
			year.str() << ":" <<		// year (YY)
		right << setw(3) << 
			dayofyearstr << ":";		// day of year

	if (calculateSeconds)
	{
		ptime currentTime(second_clock::local_time());
		time_duration timeInterval(currentTime - dateTime);

		stringstream secondsofday;
		secondsofday << left << setw(5) << timeInterval.total_seconds();
		*os << left << setw(5) << secondsofday.str();
		return;
	}

	*os << "00000";
}

template <class U>
U dateYear(const date& theDate)
{
	return theDate.year();
}

template <class U>
U dateDOY(const date& theDate)
{
	return theDate.day_of_year();
}

template <class U>
U dateMonth(const date& theDate)
{
	return theDate.month();
}

template <class U>
U dateDay(const date& theDate)
{
	return theDate.day();
}

template <class U>
U dateDaynumber(const date& theDate)
{
	return theDate.day_number();
}

template <class U>
U todayYear()
{
	return dateYear<U>(day_clock::local_day());
}

template <class U>
U todayMonth()
{
	return dateMonth<U>(day_clock::local_day());
}

template <class U>
U todayDay()
{
	return dateDay<U>(day_clock::local_day());
}

template <class U>
U todayDaynumber()
{
	return dateDaynumber<U>(day_clock::local_day());
}

template <class U, class T>
T integerToDate(const U& integer)
{
	return date(integer);
}

template <class U, class T>
T dateToInteger(const T& theDate)
{
	return dateDaynumber(theDate);
}

template <class U, class T>
T ymdToDate(const U& year, const U& month, const U& day)
{
	return date(year, month, day);
}

template <class T>
T timeImmemorial()
{
	return T(TIME_IMMEMORIAL, 1, 1);
}

template <class T>
bool isTimeImmemorial(const T& theDate)
{
	return theDate < T(TIME_IMMEMORIAL, 1, 2);
}

// cumulative days in a year required for determining decimal year
template <class P>
struct _cumulative_days_in_year_
{
	static const short cumulativeDays[2][12];
};

template <class P>
short const _cumulative_days_in_year_<P>::cumulativeDays[2][12] =
	// normal year
	{{
		0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334
	},
	// leap year
	{
		0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335
	}};

class CUMULATIVE_DAY_YEAR
	: public _cumulative_days_in_year_<void>
{
public:
};


// calendar month names
template <class P>
struct _month_names_
{
	static const char monthNames[12][12];
};

template <class P>
const char _month_names_<P>::monthNames[12][12] =
	{ 
		{"January"},
		{"February"},
		{"March"},
		{"April"},
		{"May"},
		{"June"},
		{"July"},
		{"August"},
		{"September"},
		{"October"},
		{"November"},
		{"December"}
	};

class MONTH_NAMES
	: public _month_names_<void>
{
public:
};

template <class T>
// calculate average and test for year cross over
void year_doy_Average(const T& start_year, const T& end_year,
	const T& start_doy, const T& end_doy,
	T& average_year, T& average_doy)
{
	if (start_year != end_year)
	{
		average_doy = average(start_doy, start_doy + end_doy);
		if (average_doy > DAYS_PER_YEAR)
		{
			average_year = end_year;
			average_doy -= static_cast<T>(DAYS_PER_YEAR);
		}
		else
			average_year = start_year;
		return;
	}
	
	average_year = start_year;
	
	if (start_doy != end_doy)
		average_doy = average(start_doy, end_doy);
	else
		average_doy = start_doy;
}

template <class T>
// epochFrom -> epochTo
T yearFraction(const date& epoch)
{
	// Day of year, less 0.5 day correction to 
	// bring to middle of day, divided by number 
	// of days in a year
	//return (dateDOY<T>(epoch) - 0.5) / DAYS_PER_YEAR;
	
	short leap(0);
	T days_in_year(365);
	
	if (fmod((T)epoch.year(), 400.) == 0)
	{
		leap = 1;
		days_in_year = 366;
	}
	else if (fmod((T)epoch.year(), 100.) == 0)
		leap = 0;
	else if (fmod((T)epoch.year(), 4.) == 0)
	{
		leap = 1;
		days_in_year = 366;
	}
	else
		leap = 0;

	T day_of_year(CUMULATIVE_DAY_YEAR::cumulativeDays[leap][epoch.month() - 1] + epoch.day());

	return (day_of_year - (T)0.5) / days_in_year;
}

template <class T>
// Return as a double or float the decimal year
T referenceEpoch(const date& theDate)
{
	// Epoch of interest
	return yearFraction<T>(theDate) + theDate.year();
}

template <class T>
// elapsedTime = current_epoch - reference_epoch
T elapsedTime(const T& current_epoch, const T& reference_epoch)
{
	// Elapsed time
	return (current_epoch - reference_epoch);
}

template <class T>
// elapsedTime = current_epoch - reference_epoch
T elapsedTime(const T& current_epoch, const date& reference_epoch)
{
	return elapsedTime(current_epoch, referenceEpoch<T>(reference_epoch));
}

template <class T>
// elapsedTime = current_epoch - reference_epoch
T elapsedTime(const date& current_epoch, const T& reference_epoch)
{
	return elapsedTime(referenceEpoch<T>(current_epoch), reference_epoch);
}

template <class T>
// elapsedTime = current_epoch - reference_epoch
T elapsedTime(const date& current_epoch, const date& reference_epoch)
{
	return elapsedTime(referenceEpoch<T>(current_epoch), referenceEpoch<T>(reference_epoch));
}

template <class T>
// dateString = dd.mm.yyyy
date dateFromString_safe(const string& dateString)
{
	vector<string> vdateList;
	SplitDelimitedString<string>(dateString, string("."), &vdateList);

	try {
		// boost date requires input in the following order: year, month, day
		return date(
			LongFromString<UINT32>(vdateList.at(2)),	// year
			LongFromString<UINT32>(vdateList.at(1)),	// month
			LongFromString<UINT32>(vdateList.at(0)));	// day
	}
	catch (runtime_error& e)
	{
		stringstream ss;
		ss << "dateFromString_safe(): Could not parse the date string \"" << dateString << "\"" << endl <<
			"  Details: " << e.what() << endl;
		throw boost::enable_current_exception(runtime_error(ss.str()));
	}
	catch (...)
	{
		stringstream ss;
		ss << "dateFromString_safe(): Could not parse the date string \"" << dateString << "\"." << endl <<
			"  Check that the epoch is formatted as dd.mm.yyyy" << endl;
		throw boost::enable_current_exception(runtime_error(ss.str()));
	}
}

template <class T, class S>
// dateString = ddd yyyy
T dateFromStringstream_doy_year(S& date)
{
	T parsedDate;
	// From boost_1_54_0/doc/html/date_time/date_time_io.html
	// %j :	Day of year as decimal from 001 to 366 for leap years, 001 - 365 for non-leap years.
	//		"060" => Feb-29
	// %y : Two digit year
	//		"05" => 2005
	// %Y : Four digit year
	// 		"2005" => 2005
	// But! 94 => 2094!
	// So use %Y"
	date.imbue(locale(date.getloc(), new date_input_facet("%j %Y")));
	date >> parsedDate;
	return parsedDate;	
}

template <class T, class S>
// dateString = ddd yyyy
T dateFromString_doy_year(S& date)
{
	T parsedDate;
	// From boost_1_54_0/doc/html/date_time/date_time_io.html
	// %j :	Day of year as decimal from 001 to 366 for leap years, 001 - 365 for non-leap years.
	//		"060" => Feb-29
	// %y : Two digit year
	//		"05" => 2005
	// %Y : Four digit year
	// 		"2005" => 2005
	// But! 94 => 2094!
	// So use %Y"
	stringstream ss_date;
	ss_date << date;
	ss_date.imbue(locale(ss_date.getloc(), new date_input_facet("%j %Y")));
	ss_date >> parsedDate;
	return parsedDate;	
}

template <class T, class D>
// date = yyyy.yyyyy
T dateFromDouble_doy_year(const D& date)
{
	T parsedDate;
	D year(floor(date));
	D year_fraction(date - year);
	D days_in_year(365.0);

	if (fmod(year, 400.) == 0)
		days_in_year += 1.0;
	else if (fmod(year, 4.) == 0)
		days_in_year += 1.0;
	
	// Calculate integer day from floating point.
	// If this rounds up or down, the net result is a day before
	// or a day after.  For some applications, this might cause 
	// one some concern.  But for reference frame geodesy and plate
	// tectonics, not a lot happens in a day so this is quite acceptable.
	D doy = year_fraction * days_in_year;
	size_t idoy = static_cast<size_t>(doy);
	if (idoy < 1)
		idoy = 1;

	// From boost_1_54_0/doc/html/date_time/date_time_io.html
	// %j :	Day of year as decimal from 001 to 366 for leap years, 001 - 365 for non-leap years.
	//		"060" => Feb-29
	// %y : Two digit year
	//		"05" => 2005
	// %Y : Four digit year
	// 		"2005" => 2005
	// But! 94 => 2094!
	// So use %Y"
	stringstream strdate;
	if (idoy < 10.)
		strdate << "00";
	else if (idoy < 100.)
		strdate << "0";
	strdate << idoy << " " << fixed << setprecision(0) << year;
	strdate.imbue(locale(strdate.getloc(), new date_input_facet("%j %Y")));
	strdate >> parsedDate;
	return parsedDate;	
}

template <class T>
// dateString = dd.mm.yyyy
string stringFromDate(const T& date)
{
	stringstream ss;
	ss.imbue(locale(ss.getloc(), new date_facet("%d.%m.%Y")));

	ss << date;
	return ss.str();	
}

template <class T>
// dateString = dd.mm.yyyy
string stringFromToday()
{
	stringstream ss;
	ss.imbue(locale(ss.getloc(), new date_facet("%d.%m.%Y")));

	ss << date(day_clock::local_day());
	return ss.str();	
}

template <class T>
// dateString = dd.mm.yyyy
string stringFromDate(const T& date, const string& format)
{
	stringstream ss;
	ss.imbue(locale(ss.getloc(), new date_facet(format.c_str())));

	ss << date;
	return ss.str();	
}

template <class T>
// dateString = dd.mm.yyyy
T dateFromString(const string& dateString)
{
	// The following works well when dateString is "03.02.2013", but fails
	// if dateString is "3.2.2013".  
	//
	//istringstream ss(dateString);
	//ss.imbue(locale(std::locale::classic(), new gregorian::date_input_facet("%e.%m.%Y")));
	//T date;
	//ss >> date;
	//return date;
	//
	// The Boost.Gregorian documentation specifies that:
	//  %d      Day of the month as decimal 01 to 31
	//  %e #    Like %d, the day of the month as a decimal number, but a leading zero is replaced by a space
	//	The top part of the documentation states the following warning:
	//	The flags marked with a hash sign (#) are implemented by system locale and are known to be missing on some platforms

	// Is today's date required?
	if (iequals(dateString, "today"))
		return dateFromString_safe<T>(stringFromDate<T>(T(day_clock::local_day())));

	// Parse manually.	
	return dateFromString_safe<T>(dateString);
	
}

template <class T>
// dateString = dd.mm.yyyy
string formattedDateStringFromNumericString(const string& dateString)
{
	date date(dateFromString<T>(dateString));
	
	stringstream ss;
	ss.imbue(locale(ss.getloc(), new date_facet("%A, %d %B %Y")));

	ss << date;
	return ss.str();
}

template <typename T>
T formattedDateTimeString()
{
	stringstream datetime_ss, stream;
	time_facet* p_time_output = new time_facet;
	locale special_locale (locale(""), p_time_output);
	
	// special_locale takes ownership of the p_time_output facet
	datetime_ss.imbue (special_locale);
	(*p_time_output).format("%d-%m-%Y %X");
	datetime_ss << second_clock::local_time();
	stream << datetime_ss.str().c_str();

	return stream.str();
}

template <typename S>
S formatedElapsedTime(milliseconds* elapsed_time, S app_message)
{
	ostringstream ss_time;
	ptime pt(ptime(gregorian::day_clock::local_day(), *elapsed_time));

	if (*elapsed_time < seconds(3))
	{
		time_facet* facet(new time_facet("%s"));
		ss_time.imbue(locale(ss_time.getloc(), facet));
		ss_time.str("");
		ss_time << pt << "s";			
	}
	else if (*elapsed_time < seconds(61))
	{		
		time_facet* facet(new time_facet("%S"));
		ss_time.imbue(locale(ss_time.getloc(), facet));
		ss_time.str("");
		ss_time << pt << "s";
	}
	else
		ss_time << seconds(static_cast<long>(elapsed_time->total_seconds()));

	size_t pos = string::npos;
	string time_message = ss_time.str();
	while ((pos = time_message.find("0s")) != string::npos)
		time_message = time_message.substr(0, pos) + "s";

	time_message = app_message + time_message + "."; 

	if ((pos = time_message.find(" 00.")) != string::npos)
		time_message = time_message.replace(pos, 4, " 0.");
	if ((pos = time_message.find(" 0.s")) != string::npos)
		time_message = time_message.replace(pos, 4, " 0s");

	return time_message;
}


#endif /* DNATEMPLATEDATETIMEFUNCS_H_ */
