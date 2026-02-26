//============================================================================
// Name         : dnachronutils.hpp
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
// Description  : Modern C++ date utilities using std::chrono
//============================================================================

#ifndef DNACHRONUTILS_HPP
#define DNACHRONUTILS_HPP

#include <chrono>
#include <ctime>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include <stdexcept>
#include <cmath>

namespace dynadjust {

// Type alias for date representation
using DateT = std::chrono::system_clock::time_point;

// Convert std::chrono time_point to struct tm
inline std::tm ToTm(const DateT& tp) {
    std::time_t tt = std::chrono::system_clock::to_time_t(tp);
    std::tm tm = *std::localtime(&tt);
    return tm;
}

// Create date from year, month, day
inline DateT MakeDate(int year, int month, int day) {
    std::tm tm = {0};
    tm.tm_year = year - 1900;
    tm.tm_mon = month - 1;
    tm.tm_mday = day;
    tm.tm_isdst = -1;
    std::time_t tt = std::mktime(&tm);
    return std::chrono::system_clock::from_time_t(tt);
}

// Get current local date
inline DateT Today() {
    return std::chrono::system_clock::now();
}

// Get year from date
template <typename T = int>
T DateYear(const DateT& date) {
    std::tm tm = ToTm(date);
    return static_cast<T>(tm.tm_year + 1900);
}

// Get month from date (1-12)
template <typename T = int>
T DateMonth(const DateT& date) {
    std::tm tm = ToTm(date);
    return static_cast<T>(tm.tm_mon + 1);
}

// Get day from date (1-31)
template <typename T = int>
T DateDay(const DateT& date) {
    std::tm tm = ToTm(date);
    return static_cast<T>(tm.tm_mday);
}

// Calculate day of year (1-366)
template <typename T = int>
T DateDOY(const DateT& date) {
    std::tm tm = ToTm(date);
    return static_cast<T>(tm.tm_yday + 1);
}

// Check if year is leap year
inline bool IsLeapYear(int year) {
    return (year % 400 == 0) || (year % 4 == 0 && year % 100 != 0);
}

// Format date in SINEX format: YY:DOY:SSSSS
template <typename T>
void DateSINEXFormat(T* os, const DateT& theDate, bool calculateSeconds = false) {
    // Get two-digit year
    int year = DateYear<int>(theDate) % 100;
    
    // Get day of year (1-based)
    int doy = DateDOY<int>(theDate);
    
    // Format YY:DOY:
    *os << std::right << std::setfill('0') << std::setw(2) << year << ":"
        << std::right << std::setfill('0') << std::setw(3) << doy << ":";
    
    if (calculateSeconds) {
        // Calculate seconds since midnight
        auto now = std::chrono::system_clock::now();
        auto time_since_epoch = now.time_since_epoch();
        auto today_duration = time_since_epoch % std::chrono::hours(24);
        auto seconds = std::chrono::duration_cast<std::chrono::seconds>(today_duration).count();

        *os << std::right << std::setw(5) << seconds;
    } else {
        *os << "00000";
    }

    // Restore default fill character so callers aren't affected
    *os << std::setfill(' ');
}

// Calculate year fraction (0.0 to 1.0)
template <typename T = double>
T YearFraction(const DateT& date) {
    int year = DateYear<int>(date);
    int doy = DateDOY<int>(date);
    T days_in_year = IsLeapYear(year) ? 366.0 : 365.0;
    return (static_cast<T>(doy) - 0.5) / days_in_year;
}

// Calculate reference epoch (year.fraction)
template <typename T = double>
T ReferenceEpoch(const DateT& date) {
    return DateYear<T>(date) + YearFraction<T>(date);
}

// Calculate elapsed time between dates
template <typename T = double>
T ElapsedTime(const DateT& current, const DateT& reference) {
    return ReferenceEpoch<T>(current) - ReferenceEpoch<T>(reference);
}

// Parse date from string "dd.mm.yyyy"
inline DateT DateFromStringSafe(const std::string& dateString) {
    std::vector<std::string> parts;
    std::stringstream ss(dateString);
    std::string part;
    
    while (std::getline(ss, part, '.')) {
        parts.push_back(part);
    }
    
    if (parts.size() != 3) {
        throw std::runtime_error("DateFromStringSafe(): Invalid date format. Expected dd.mm.yyyy");
    }
    
    try {
        int day = std::stoi(parts[0]);
        int month = std::stoi(parts[1]);
        int year = std::stoi(parts[2]);
        return MakeDate(year, month, day);
    } catch (const std::exception& e) {
        throw std::runtime_error(std::string("DateFromStringSafe(): Could not parse date string \"") + 
                                 dateString + "\": " + e.what());
    }
}

// Average year and day of year
template <typename T>
void YearDoyAverage(const T& start_year, const T& end_year,
                      const T& start_doy, const T& end_doy,
                      T& average_year, T& average_doy) {
    if (start_year != end_year) {
        average_doy = (start_doy + start_doy + end_doy) / 2;
        if (average_doy > 365) {
            average_year = end_year;
            average_doy -= 365;
        } else {
            average_year = start_year;
        }
    } else {
        average_year = start_year;
        average_doy = (start_doy + end_doy) / 2;
    }
}

// Convert boost::gregorian::date to std::chrono
inline DateT FromBoostDate(int year, int month, int day) {
    return MakeDate(year, month, day);
}

} // namespace dynadjust

#endif // DNACHRONUTILS_HPP