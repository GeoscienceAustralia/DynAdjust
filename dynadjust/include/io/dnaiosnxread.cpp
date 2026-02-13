//============================================================================
// Name         : dnaiosnxread.cpp
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
// Description  : DynAdjust SINEX file read operations
//============================================================================

#include <include/io/dnaiosnx.hpp>
#include <include/functions/dnastringfuncs.hpp>
#include <include/functions/dnatemplatedatetimefuncs.hpp>
#include <include/functions/dnastrutils.hpp>
#include <include/measurement_types/dnagpspoint.hpp>

namespace dynadjust { 
namespace iostreams {

void DnaIoSnx::ParseSinex(std::ifstream** snx_file, const std::string& fileName, vdnaStnPtr* vStations, PUINT32 stnCount, 
					vdnaMsrPtr* vMeasurements, PUINT32 msrCount, PUINT32 clusterID,
					StnTally& parsestn_tally, MsrTally& parsemsr_tally, UINT32& fileOrder,
					CDnaDatum& datum, bool applyDiscontinuities, 
					v_discontinuity_tuple* stn_discontinuities, bool& m_discontsSortedbyName,
					UINT32& lineNo, UINT32& columnNo, _PARSE_STATUS_& parseStatus)
{
	parseStatus = PARSE_SUCCESS;
	
	(*stnCount) = 0;
	(*msrCount) = 0;
	parsestn_tally.initialise();
	parsemsr_tally.initialise();

	lineNo = 0;
	columnNo = 0;

	containsVelocities_ = false;
	applyDiscontinuities_ = (applyDiscontinuities && !stn_discontinuities->empty());
	containsDiscontinuities_ = false;

	// read header line and extract epoch
	if (!ParseSinexHeader(snx_file, datum, lineNo))
	{
		std::stringstream ss;
		ss << "ParseSinex(): The SINEX file  " << fileName << " did not contain a header record." << std::endl;
		throw std::runtime_error(ss.str());
	}

	try {
		
		// read data
		ParseSinexData(snx_file, vStations, stnCount, vMeasurements, msrCount, clusterID,
					parsestn_tally, parsemsr_tally, datum, 
					stn_discontinuities, m_discontsSortedbyName,
					fileOrder, lineNo, columnNo);
	}
	catch (const std::ios_base::failure& f) {
		if ((*snx_file)->eof())
			return;
		std::stringstream ss;
		ss << "ParseSinex(): An error was encountered when attempting to read SINEX file  " << fileName << "." << std::endl << "  " << f.what();
		throw std::runtime_error(ss.str());
	}
	catch (const std::runtime_error& f) {
		std::stringstream ss;
		ss << "  - line " << lineNo;
		ss << ", column " <<  columnNo << std::endl;
		ss << "  - " << f.what();
		throw std::runtime_error(ss.str());
	}
}
	

void DnaIoSnx::ParseDiscontinuityFile(std::ifstream* snx_file, const std::string& fileName,
	v_discontinuity_tuple* stn_discontinuities, bool& m_discontsSortedbyName,
	UINT32& lineNo, UINT32& columnNo, _PARSE_STATUS_& parseStatus)
{
	parseStatus = PARSE_SUCCESS;
	
	lineNo = 0;
	columnNo = 0;

	CDnaDatum datum;

	containsVelocities_ = false;
	containsDiscontinuities_ = false;

	std::string sBuf, tmp;

	try {
		// As the purpose of the discontinuity file is to provide dates on when discontinuities occur,
		// there is no real imperative to capture the date from the header line
		// Check for it anyway
		ParseSinexHeader(&snx_file, datum, lineNo);

		while (!snx_file->eof())			// while EOF not found
		{
			lineNo++;
			getline((*snx_file), sBuf);

			// End of data?
			if (iequals(sBuf.substr(0, 7), ENDSNX))
				break;

			try {

				if (sBuf.substr(0, 1) == "+")
					if (sBuf.substr(0, 23) == "+SOLUTION/DISCONTINUITY")
						ParseSinexDiscontinuities(snx_file, stn_discontinuities, m_discontsSortedbyName,
							lineNo, columnNo);
			}
			catch (...) {
				std::stringstream ss;
				ss << "ParseSinexData(): Could not extract data from the record:  " << std::endl << "    " << sBuf << ".";
				throw std::runtime_error(ss.str());
			}		
		}
	
	}
	catch (const std::ios_base::failure& f) {
		if (snx_file->eof())
			return;
		std::stringstream ss;
		ss << "parse_discontinuity_file(): An error was encountered when attempting to read SINEX file  " << fileName << "." << std::endl << "  " << f.what();
		throw std::runtime_error(ss.str());
	}
	catch (const std::runtime_error& f) {
		std::stringstream ss;
		ss << "  - line " << lineNo;
		ss << ", column " <<  columnNo << std::endl;
		ss << "  - " << f.what();
		throw std::runtime_error(ss.str());
	}
}
	

// Expects yy as two digit year, and doy as three digit day of year
std::stringstream DnaIoSnx::ParseDateFromYyDoy(const UINT32& yy, const UINT32& doy, DATE_FORMAT_TYPE date_type, const std::string& separator)
{
	std::stringstream ss_doy, ss_year, ss_date;

	// format the day of year so that boost date_input_facet is happy
	if (doy < 100)
		ss_doy << "0";
	if (doy < 10)
		ss_doy << "0";

	ss_doy << doy;

	// create four-character year?
	switch (date_type)
	{
	case yyyy_doy:
	case doy_yyyy:
		//	SINEX time specification:
		//	 _______________________________________________________________ 
		//	| Time:	| YY:DDD:SSSSS. "UTC"					|	I2.2,		|
		//	|		| YY = last 2 digits of the year,		|	1H:,I3.3,	|
		//	|		| if YY <= 50 implies 21-st century,	|	1H:,I5.5	|
		//	|		| if YY > 50 implies 20-th century,		|				|
		//	|		| DDD = 3-digit day in year,			|				|
		//	|		| SSSSS = 5-digit seconds in day.		|				|
		//	|_______|_______________________________________|_______________|
		if (yy <= 50)
			ss_year << "20";
		else //(average_year > 50)
			ss_year << "19";
		break;
	default:
		break;
	}

	if (yy < 10)
		ss_year << "0";
	ss_year << yy;

	switch (date_type)
	{
	case yyyy_doy:
	case yy_doy:
		ss_date << ss_year.str() << separator << ss_doy.str();
		break;
	case doy_yy:
	case doy_yyyy:
		ss_date << ss_doy.str() << separator << ss_year.str();
		break;
	}

	return ss_date;
}

// date is a string "YY:DDD:SSSSS"
std::stringstream DnaIoSnx::ParseDateFromString(const std::string& date_str, DATE_FORMAT_TYPE date_type, DATE_TERMINAL_TYPE terminal_type, const std::string& separator)
{
	// Year
	UINT32 yy = LongFromString<UINT32>(date_str.substr(0, 2));
	
	// Day of year (assume a colon separates the two values)
	UINT32 doy = LongFromString<UINT32>(date_str.substr(3, 3));

	if (yy == 0 && doy == 0)
	{
		boost::gregorian::date today(boost::gregorian::day_clock::local_day());
		std::stringstream date_ss;

		switch (terminal_type)
		{
		case date_from:
			// Set a date before the advent of GPS
			date_ss << "001 " << TIME_IMMEMORIAL;
			return date_ss;
		case date_to:
			// Set a date 100 years into the future
			date_ss << "365 " << (today.year() + 100);
			return date_ss;
		}
	}

	return ParseDateFromYyDoy(yy, doy, date_type, separator);
}

bool DnaIoSnx::ParseSinexHeader(std::ifstream** snx_file, CDnaDatum& datum, UINT32& lineNo)
{
	/*
	%=SNX 2.00 AUS 11:182:57860 IGS 11:157:00000 11:157:86370 P 00054 0 S          
	*/
	char c = static_cast<char>((*snx_file)->peek());
	if (c != '%')
	{
		// The SINEX file doesn't have a header.
		return false;
	}	

	std::string sBuf, sinex_code;
	
	UINT32 average_year, average_doy;
	std::stringstream ss;

	lineNo++;
	try {
		getline((**snx_file), sBuf);
	}
	catch (...) {
		ss << "ParseSinexHeader(): Could not read from the SINEX file.  " << std::endl << ".";
		throw std::runtime_error(ss.str());
	}

	try {
		// capture epoch of data, calculate average and test for year cross over
		year_doy_Average(
			LongFromString<UINT32>(sBuf.substr(32, 2)), 
			LongFromString<UINT32>(sBuf.substr(45, 2)), 
			LongFromString<UINT32>(sBuf.substr(35, 3)), 
			LongFromString<UINT32>(sBuf.substr(48, 3)),
			average_year, 
			average_doy);

		ss = ParseDateFromYyDoy(average_year, average_doy, doy_yyyy, std::string(" "));
		datum.SetEpoch(dateFromStringstream_doy_year<boost::gregorian::date, std::stringstream>(ss));
	}
	catch (...)
	{
		ss << "ParseSinexHeader(): Could not extract date and time information from the header record:  " << std::endl << "    " << sBuf << ".";
		throw std::runtime_error(ss.str());
	}

	return true;
}
	

void DnaIoSnx::ParseSinexData(std::ifstream** snx_file, vdnaStnPtr* vStations, PUINT32 stnCount, 
					vdnaMsrPtr* vMeasurements, PUINT32 msrCount, PUINT32 clusterID,
					StnTally& parsestn_tally, MsrTally& parsemsr_tally, CDnaDatum& datum, 
					v_discontinuity_tuple* stn_discontinuities, bool& m_discontsSortedbyName,
					UINT32& fileOrder, UINT32& lineNo, UINT32& columnNo)
{
	std::string sBuf, tmp;
	
	vStations->clear();
	vMeasurements->clear();

	siteIDsRead_ = false;
	solutionEpochsRead_ = false;

	while (!(*snx_file)->eof())			// while EOF not found
	{
		lineNo++;
		getline((**snx_file), sBuf);

		// End of data?
		if (iequals(sBuf.substr(0, 7), ENDSNX))
			break;

		try {

			if (sBuf.substr(0, 1) == "+")
				ParseSinexBlock(snx_file, sBuf.c_str(), vStations, stnCount, vMeasurements, msrCount, clusterID,
					parsestn_tally, parsemsr_tally, datum, fileOrder, lineNo, columnNo);
		}
		catch (std::runtime_error& e) {
			std::stringstream ss;
			ss << "ParseSinexData(): Error parsing SINEX file:  " << std::endl << 
				"    " << e.what();
			throw std::runtime_error(ss.str());
		}
		catch (...) {
			std::stringstream ss;
			ss << "ParseSinexData(): Could not extract data from the record:  " << std::endl << "    " << sBuf << ".";
			throw std::runtime_error(ss.str());
		}		
	}

	// Update station names in station and measurement vectors using discontinuity file 
	if (applyDiscontinuities_)
		// true if  --discontinuity-file sinex_discontinuities.snx
		FormatStationNames(stn_discontinuities, m_discontsSortedbyName, 
			vStations, vMeasurements);
}
	

void DnaIoSnx::FormatStationNames(v_discontinuity_tuple* stn_discontinuities, bool& m_discontsSortedbyName,
	vdnaStnPtr* vStations, vdnaMsrPtr* vMeasurements)
{
	UINT32 site, doy, year;
	std::string site_name;

#ifdef _MSDEBUG
	std::string date_str;
#endif

	_it_vdiscontinuity_tuple _it_discont(stn_discontinuities->end()), _it_discont_site;

	boost::gregorian::date site_date;

	if (!m_discontsSortedbyName)
	{
		std::sort(stn_discontinuities->begin(), stn_discontinuities->end(),
			CompareSiteTuplesByName<discontinuity_tuple>());
		m_discontsSortedbyName = true;
	}
	std::sort(siteOccurrence_.begin(), siteOccurrence_.end(),
		CompareSiteTuplesByName<site_id>());

	_it_discont = stn_discontinuities->begin();

	for (site=0; site<siteOccurrence_.size(); ++site)
	{
		site_name = siteOccurrence_.at(site).site_name;

		// At this point, _it_discont will point to either:
		//  1. the first element in the vector
		//  2. the first occurrence of discontinuity for site_name
		// To prevent searching through the entire array again on 2. , check if 
		// _it_discont->site_name = site_name
		if (!equals(site_name, _it_discont->site_name))
		{
			// Save the iterator
			_it_discont_site = _it_discont;

			if ((_it_discont = binary_search_discontinuity_site(
				_it_discont, 
				stn_discontinuities->end(), 
				site_name)) == stn_discontinuities->end())
			{
				// Site not found in discontinuity file, reset iterator to last saved state
				_it_discont = _it_discont_site;
				continue;
			}
		}
		
		// Is this site a discontinuity site?
		if (!_it_discont->discontinuity_exists)
			continue;

		// Capture the (start date of the site)
		site_date = dateFromString_doy_year<boost::gregorian::date, std::string>(siteOccurrence_.at(site).formatted_date);
		
#ifdef _MSDEBUG
		//date_str = stringFromDate(_it_discont->date_start, "%y:%j:00000");
#endif
		_it_discont_site = _it_discont;

		// find the next occurrence of this site
		while (equals(site_name, _it_discont_site->site_name))
		{
			// Test if the start epoch of this site is within this discontinuity window
			if (site_date >= _it_discont_site->date_start &&
				site_date < _it_discont_site->date_end)
//			if (siteOccurrence_.at(site).solution_id == _it_discont_site->solution_id)
			{
				std::stringstream ss;
				doy = _it_discont_site->date_start.day_of_year();
				year = _it_discont_site->date_start.year();
				ss << siteOccurrence_.at(site).site_name << "_";

				// format date using the start of the window defining the 
				// period up to which a discontinuity has been identified
				ss << year;
				if  (doy < 100)
					ss << "0";
				if (doy < 10)
					ss << "0";
				ss << doy;
				
				siteOccurrence_.at(site).formatted_name = ss.str();
				//TRACE("Site %s\n", siteOccurrence_.at(site).formatted_name.c_str());
				break;
			}
			_it_discont_site++;
		}
	}

	if (siteOccurrence_.size() != vStations->size() ||
		siteOccurrence_.size() != vMeasurements->at(0)->GetTotal())
	{
		std::stringstream ss;
		ss << "format_station_names(): Inconsistent number of stations in SINEX file.";
		throw std::runtime_error(ss.str());
	}

	_it_vsite_id_tuple _it_site(siteOccurrence_.begin());

	// sort on file index (default sort)
	std::sort(siteOccurrence_.begin(), siteOccurrence_.end(),
		CompareSiteTuples<site_id>());

	// Now update station vector
	for_each(vStations->begin(), vStations->end(),
		[&_it_site](const dnaStnPtr& stn) {
			if (!equals(_it_site->site_name, _it_site->formatted_name))
			{
				stn->SetOriginalName();
				stn->SetName(_it_site->formatted_name);
			}
			_it_site++;
	});

	// re-initialise iterator
	_it_site = siteOccurrence_.begin();
	
	_it_vdnamsrptr _it_msr(vMeasurements->begin());
	std::vector<CDnaGpsPoint>* vgpsPnts(_it_msr->get()->GetPoints_ptr());

	// Now update measurement vector
	for_each(vgpsPnts->begin(), vgpsPnts->end(),
		[&_it_site](CDnaGpsPoint& pnt) {
			if (!equals(_it_site->site_name, _it_site->formatted_name))
				pnt.SetFirst(_it_site->formatted_name);
			_it_site++;
	});
}



	
void DnaIoSnx::ParseSinexBlock(std::ifstream** snx_file, const char* sinexRec, vdnaStnPtr* vStations, PUINT32 stnCount, 
					vdnaMsrPtr* vMeasurements, PUINT32 msrCount, PUINT32 clusterID,
					StnTally& parsestn_tally, MsrTally& parsemsr_tally, CDnaDatum& datum, 
					UINT32& fileOrder, UINT32& lineNo, UINT32& columnNo)
{
	// write comments
	if (strncmp(sinexRec, "+FILE/REFERENCE", 15) == 0 || 
		strncmp(sinexRec, "+SITE/RECEIVER", 14) == 0 || 
		strncmp(sinexRec, "+SITE/ANTENNA", 13) == 0 || 
		strncmp(sinexRec, "+SITE/GPS_PHASE_CENTER", 22) == 0 || 
		strncmp(sinexRec, "+SITE/ECCENTRICITY", 18) == 0)
	{
		// nothing to do.
		return;
	}

	// According to the SINEX standard:
	// - BLOCKS start with "+" in the first column followed by standardized block labels.
	// - Each BLOCK ends with "-" and the block label. 
	// - BLOCK data starts in column 2 or higher. 
	// - BLOCKS can be in any order, provided that they start with (+) and 
	//   end with (-) BLOCK labels.
	//
	// Therefore, don't assume the SOLUTION/* block order is followed.
	// However, assume SITE/ID comes first, so that discontinuities can be parsed in an
	// intelligent way.

	if (strncmp(sinexRec, "+SITE/ID", 8) == 0)
	{
		// OK, found general site list
		//parse_sinex_sites(snx_file, lineNo, columnNo);
		return;
	}

	if (strncmp(sinexRec, "+SOLUTION/EPOCHS", 16) == 0)
	{
		// OK, found epochs associated with stations
		ParseSinexEpochs(snx_file, lineNo, columnNo);
		return;
	}

	if (strncmp(sinexRec, "+SOLUTION/ESTIMATE", 18) == 0)
	{
		// OK, found some stations
		ParseSinexStn(snx_file, sinexRec, vStations, stnCount,
			parsestn_tally, datum, fileOrder, lineNo, columnNo);
		return;
	}

	if (strncmp(sinexRec, "+SOLUTION/MATRIX_ESTIMATE", 25) == 0)
	{
		// Ok, found some measurements
		ParseSinexMsr(snx_file, sinexRec, vStations, vMeasurements, clusterID, msrCount,
			parsemsr_tally, datum, lineNo);
		return;
	}
}
	

void DnaIoSnx::ParseSinexDiscontinuities(std::ifstream* snx_file, 
					v_discontinuity_tuple* stn_discontinuities, bool& m_discontsSortedbyName,
					UINT32& lineNo, UINT32& columnNo)
{
	std::string sBuf("+"), stn, stn_prev(""), model;
	UINT32 file_rec(0), solution_id;
	boost::gregorian::date discont_start_date, discont_end_date;

	std::stringstream ss, discont_from, discont_to;

	stn_discontinuities->clear();
	m_discontsSortedbyName = false;

	while (sBuf.at(0) != '-')
	{
		lineNo++;
		getline((*snx_file), sBuf);

		// A comment
		if (sBuf.at(0) == '*')
			continue;

		if (sBuf.at(0) == '-')
			break;

		try {
			// station
			stn = trimstr(sBuf.substr(1, 4));				// station name

			// solution ID
			solution_id = val_uint<UINT32, std::string>(
				trimstr(sBuf.substr(9, 4)));				// solution id

			// model
			// Commonly used models for time series analysis are:
			// "M" Models Codes for coordinates time series Analysis:
			//      P = Position
			//	    V = Velocity
			//	    A = Annual Periodicity
			//	    S = Semi-Annual Periodicity
			//	    E = Exponential Decay
			//	    X = Exclude/Delete
			model = trimstr(sBuf.substr(42, 1));				// model

			// Is this a discontinuity in position?  If not, continue to next
			if (!iequals(model, "P"))
				continue;

			// If there are multiple solutions, then there must be a discontinuity, but only if
			// the model is the same.  For this version of DynAdjust, work with position only (i.e. disregard velocity)
			if (solution_id > 1)
				containsDiscontinuities_ = true;

			////////////////////////////////////////////////////////////////////////////////////////////
			// Dates FROM and TO define the temporal extents of time series divided by discontinuities.
			// To understand how to apply these dates, the following rules apply:
			//   - If FROM and TO both equal 00:000:00000, there is no known (or recorded) discontinuity 
			//     at this site (in which case, this site does not need renaming)
			//   - If FROM equals 00:000:00000 and TO is a realistic date, then this is the first leg
			//     in a sequence of time series, in which case solution_id should equal 1
			//   - If FROM and TO are realistic dates, then this is an intermediate leg in a sequence
			//   - If FROM is a realistic date and TO equals 00:000:00000, then this is the last leg in
			//     the sequence of time series.

			// start epoch (yy:doy:sssss)
			discont_from = ParseDateFromString(trimstr(sBuf.substr(16, 6)), doy_yyyy, date_from, std::string(" "));
			discont_start_date = dateFromStringstream_doy_year<boost::gregorian::date, std::stringstream>(discont_from);

			// end epoch (yy:doy:sssss)
			discont_to = ParseDateFromString(trimstr(sBuf.substr(29, 6)), doy_yyyy, date_to, std::string(" "));
			discont_end_date = dateFromStringstream_doy_year<boost::gregorian::date, std::stringstream>(discont_to);

			////////////////////////////////////////////////////////////////////////////////////////////

			stn_discontinuities->push_back(
				discontinuity_tuple_t<UINT32, std::string, boost::gregorian::date, bool>(
				file_rec,			// file_index 
				solution_id,		// solution_id
				stn,				// site_name
				discont_start_date,	// discontinuity start
				discont_end_date,	// discontinuity end
				false) 				// discontinuity_exists
				);

			file_rec++;
		}
		catch (...) {
			ss.str("");
			columnNo = 1;
			ss << "ParseSinexDiscontinuities(): Could not extract station name from the record:  " << std::endl << "    " << sBuf << ".";
			throw std::runtime_error(ss.str());
		}			
	}

	// No discontinuities to handle?
	if (!containsDiscontinuities_)
		return;	

	// sort the station discontinuities on station name
	std::sort(stn_discontinuities->begin(), stn_discontinuities->end(), 
		CompareSiteTuplesByName<discontinuity_tuple>());
	
	// Identify if discontinuity exists (by equal station names) and 
	// set flag to indicate that this station is a known discontinuity station
	_it_vdiscontinuity_tuple it_discont, it_discont_next;
	for (it_discont = stn_discontinuities->begin();
		it_discont != stn_discontinuities->end(); it_discont++)
	{
		it_discont_next =  it_discont + 1;

		if (it_discont_next == stn_discontinuities->end())
			break;
		
		// Are station names equal?
		if (equals(it_discont->site_name, it_discont_next->site_name))
		{
			// Then this is a discontinuity site!
			it_discont->discontinuity_exists = true;
			it_discont_next->discontinuity_exists = true;
		}
	}

	// re-sort the station discontinuities back to the original file order
	std::sort(stn_discontinuities->begin(), stn_discontinuities->end(),
		CompareSiteTuples<discontinuity_tuple>());
	m_discontsSortedbyName = false;
}
	

// Read sites
void DnaIoSnx::ParseSinexSites(std::ifstream** snx_file, UINT32& lineNo, UINT32& columnNo)
{
	std::string sBuf("+"), stn, domes;
	
	// Has the SOLUTION/EPOCHS block already been read?  If so, no need
	// to continue with this block - just skip to the end.
	if (solutionEpochsRead_)
	{
		// Skip to the end of the SITE/ID block
		while (sBuf.at(0) != '-')
		{
			lineNo++;
			getline((**snx_file), sBuf);
		}
		siteIDsRead_ = true;
		return;
	}
	
	siteOccurrence_.clear();

	UINT32 file_rec(0);
	std::stringstream ss;

	while (sBuf.at(0) != '-')
	{
		lineNo++;
		getline((**snx_file), sBuf);

		// A comment
		if (sBuf.at(0) == '*')
			continue;

		if (sBuf.at(0) == '-')
			break;

		try {
			// station
			stn = trimstr(sBuf.substr(1, 4));

			// domes
			domes = trimstr(sBuf.substr(9, 9));

			// set default values (amended later)
			siteOccurrence_.push_back(
				site_id_tuple_t<UINT32, std::string, bool>(
					file_rec,		// file_index, 
					1,				// solution_id (set to 1 since SITE/ID doesn't hold solution_id)
					stn,			// station name
					stn,			// formatted_name
					"",				// formatted_date
					false) 			// false (amended later to be true if last occurrence)			
				);

			file_rec++;
		}
		catch (...) {
			ss.str("");
			columnNo = 1;
			ss << "ParseSinexDiscontinuities(): Could not extract station name from the record:  " << std::endl << "    " << sBuf << ".";
			throw std::runtime_error(ss.str());
		}			
	}

	siteIDsRead_ = true;

	ReduceSinexSites();
}
	

// Read station epochs
void DnaIoSnx::ParseSinexEpochs(std::ifstream** snx_file, UINT32& lineNo, UINT32& columnNo)
{
	std::string sBuf("+"), stn, site, date_start;
	UINT32 file_rec(0), solution_id;

	std::stringstream ss;

	if (!siteIDsRead_)
		siteOccurrence_.clear();


	while (sBuf.at(0) != '-')
	{
		lineNo++;
		getline((**snx_file), sBuf);

		// A comment
		if (sBuf.at(0) == '*')
			continue;

		if (sBuf.at(0) == '-')
			break;

		try {
			// station
			stn = trimstr(sBuf.substr(1, 4));
			// solution ID
			solution_id = val_uint<UINT32, std::string>(trimstr(sBuf.substr(9, 4)));
			// Get start epoch (yy:doy:sssss) of data window used to process this site
			date_start = ParseDateFromString(trimstr(sBuf.substr(16, 6)), doy_yyyy, date_from, std::string(" ")).str();
		}
		catch (...) {
			ss.str("");
			columnNo = 1;
			ss << "ParseSinexEpochs(): Could not extract station name from the record:  " << std::endl << "    " << sBuf << ".";
			throw std::runtime_error(ss.str());
		}

		// Perform some file consistency checks
		if (siteIDsRead_)
		{
			try {
				// Is the count the same?  This will cause an exception of file_rec is equal to or greater than 
				// siteOccurrence_.size()
				site = siteOccurrence_.at(file_rec).site_name;
			}
			catch (...) {
				ss.str("");
				columnNo = 1;
				ss << "ParseSinexEpochs(): The number of sites in SITE/ID and SOLUTION/EPOCHS" << std::endl <<
					"    is inconsistent:  " << std::endl << 
					"    " << "SITE/ID block has " << siteOccurrence_.size() << " sites" << std::endl <<
					"    " << "SOLUTION/EPOCHS block contains additional sites not listed in SITE/ID.  Next record is:" << std::endl <<
					"    " << sBuf << ".";
				throw std::runtime_error(ss.str());
			}
		
			if (!equals(site, stn))
			{
				ss.str("");
				columnNo = 1;
				ss << "ParseSinexEpochs(): The order of sites in SITE/ID and SOLUTION/EPOCHS" << std::endl <<
					"    is inconsistent:  " << std::endl << 
					"    " << "Index " << file_rec << " in SITE/ID block is " << siteOccurrence_.at(file_rec).site_name << std::endl <<
					"    " << "Index " << file_rec << " in SOLUTION/EPOCHS block is " << stn << ".";
				throw std::runtime_error(ss.str());			
			}
		}

		try {	
			// Has the SITE/ID block been read before the SOLUTION/EPOCHS block?  
			// If so, simply update the records
			if (siteIDsRead_)
			{
				siteOccurrence_.at(file_rec).solution_id = solution_id;
				// The mean date is only necessary for handling discontinuities
				siteOccurrence_.at(file_rec).formatted_date = date_start;	
			}
			else
			{
				// Create new elements and add to the vector
				siteOccurrence_.push_back(
					site_id_tuple_t<UINT32, std::string, bool>(
					file_rec,		// file_index, 
					solution_id,	// solution_id
					stn,			// station name
					stn,			// formatted_name
					date_start,		// formatted_date
					false) 			// false (amended later to be true if last occurrence)			
					);
			}			

			file_rec++;
		}
		catch (...) {
			ss.str("");
			columnNo = 1;
			ss << "ParseSinexEpochs(): Failed to add a new record to the site occurrence list." << std::endl <<
				"    " << sBuf << ".";
			throw std::runtime_error(ss.str());
		}
	}

	solutionEpochsRead_ = true;

	if (!siteIDsRead_)
		ReduceSinexSites();
}
	
void DnaIoSnx::ReduceSinexSites()
{
	// sort the sites on station name
	std::sort(siteOccurrence_.begin(), siteOccurrence_.end(), 
		CompareSiteTuplesByName<site_id>());

	// Identify if discontinuity exists from internal evidence (i.e. equal station names).  If so,
	// increment the station occurrence to indicate that this station is duplicated because
	// of discontinuity.
	_it_vsite_id_tuple it_site, it_site_next;
	for (it_site = siteOccurrence_.begin();
		it_site != siteOccurrence_.end(); it_site++)
	{
		it_site_next =  it_site + 1;

		if (it_site_next == siteOccurrence_.end())
		{
			// this is the last occurrence
			it_site->last_occurrence = true;
			break;
		}

		// if this station occurs again:
		//   - increment the next occurrence of it
		//   - set the index of the first occurrence
		if (equals(it_site->site_name, it_site_next->site_name))
		{
			// increment next occurrence.  Note, if SOLUTION/EPOCHS block has been read, 
			// this will have already been set
			if (!solutionEpochsRead_)
				it_site_next->solution_id = it_site->solution_id + 1;

			// If this code is reached, then there must be a discontinuities in the file!
			containsDiscontinuities_ = true;
		}
		else
			// this is the last occurrence
			it_site->last_occurrence = true;

	}

	// sort on file index (default sort)
	std::sort(siteOccurrence_.begin(), siteOccurrence_.end(),
		CompareSiteTuples<site_id>());
}
	
void DnaIoSnx::ParseSinexStn(std::ifstream** snx_file, const char* sinexRec, vdnaStnPtr* vStations, PUINT32 stnCount,
					StnTally& parsestn_tally, CDnaDatum& datum,
					UINT32& fileOrder, UINT32& lineNo, UINT32& columnNo)
{
	std::stringstream ss;

	if (applyDiscontinuities_ && !solutionEpochsRead_)
	{
		ss.str("");
		ss << "parse_sinex_stn(): Cannot apply discontinuities to the stations" << std::endl <<
			"    if the station epochs have not been loaded beforehand.  To rectify this problem," << std::endl <<
			"    reformat the SINEX file so that the +SOLUTION/EPOCHS block appears before the +SOLUTION/ESTIMATE block.";
		throw std::runtime_error(ss.str());
	}

	std::string sBuf(sinexRec), site, tmp, stn;

	dnaStnPtr stn_ptr;	
	UINT32 yy, doy, file_rec(0);

	//fileOrder = 0;
	uniqueStationCount_ = 0;

	while (sBuf.at(0) != '-')
	{
		lineNo++;

		getline((**snx_file), sBuf);

		// A comment
		if (sBuf.at(0) == '*')
			continue;
		
		if (sBuf.at(0) == '-')
			break;

		if (sBuf.substr(7, 3) == "VEL")
			containsVelocities_ = true;

		if (sBuf.substr(7, 3) != "STA")
			continue;

		// Capture the X coordinate
		if (sBuf.substr(7, 4) == "STAX")
		{
			try {
				// get the measurement station name (will be 4 characters)
				stn = trimstr(sBuf.substr(14, 4));

				// reset CDnaStation instance (stn_ptr) with new station name
				stn_ptr.reset(new CDnaStation(
					stn,									// station name
					"FFF",									// Constraints
					XYZ_type, 0.0, 0.0, 0.0, 0.0, "",		// Type, X, Y, Z, Height, HemisphereZOne
					trimstr(sBuf.substr(14, 4)),			// description
					ss.str()));								// comment

				stn_ptr->SetReferenceFrame(datum.GetName());
				stn_ptr->SetEpsg(datum.GetEpsgCode_s());

				yy = LongFromString<UINT32>(sBuf.substr(27, 2)), 
				doy = LongFromString<UINT32>(sBuf.substr(30, 3)), 
				ss = ParseDateFromYyDoy(yy, doy, doy_yyyy, std::string(" "));
				stn_ptr->SetEpoch(stringFromDate(dateFromStringstream_doy_year<boost::gregorian::date, std::stringstream>(ss)));

				stn_ptr->SetfileOrder(fileOrder++);
				stn_ptr->SetXAxis_d(DoubleFromString<double>(trimstr(sBuf.substr(47, 21))));
				stn_ptr->SetXAxisStdDev_d(DoubleFromString<double>(trimstr(sBuf.substr(68, 12))));
			}
			catch (const std::runtime_error& f) {
				ss.str("");
				ss << "  - line " << lineNo;
				ss << ", column " <<  columnNo << std::endl;
				ss << "  - " << f.what();
				throw std::runtime_error(ss.str());
			}
			catch (...) {
				ss.str("");
				columnNo = 47;
				ss << "parse_sinex_stn(): Could not extract X coordinate estimate from the record:  " << std::endl << "    " << sBuf << ".";
				throw std::runtime_error(ss.str());
			}
			continue;
		}

		// Perform some file consistency checks
		try {
			// Is the count the same?  This will cause an exception of file_rec is equal to or greater than 
			// siteOccurrence_.size()
			site = siteOccurrence_.at(file_rec).site_name;
		}
		catch (...) {
			ss.str("");
			columnNo = 1;
			ss << "parse_sinex_stn(): The number of sites in SOLUTION/EPOCHS and SOLUTION/ESTIMATE" << std::endl <<
				"    is inconsistent:  " << std::endl << 
				"    " << "SOLUTION/EPOCHS block has " << siteOccurrence_.size() << " sites" << std::endl <<
				"    " << "SOLUTION/ESTIMATE block contains additional sites not listed in SOLUTION/EPOCHS.  Next record is:" << std::endl <<
				"    " << sBuf << ".";
			throw std::runtime_error(ss.str());
		}

		if (!equals(site, stn))
		{
			ss.str("");
			columnNo = 1;
			ss << "parse_sinex_stn(): The order of sites in SOLUTION/EPOCHS and SOLUTION/ESTIMATE" << std::endl <<
				"    is inconsistent:  " << std::endl << 
				"    " << "Index " << file_rec << " in SOLUTION/EPOCHS block is " << siteOccurrence_.at(file_rec).site_name << std::endl <<
				"    " << "Index " << file_rec << " in SOLUTION/ESTIMATE block is " << stn << ".";
			throw std::runtime_error(ss.str());			
		}

		// Capture the Y coordinate
		if (sBuf.substr(7, 4) == "STAY")
		{
			try {
				stn_ptr->SetYAxis_d(DoubleFromString<double>(trimstr(sBuf.substr(47, 21))));
				stn_ptr->SetYAxisStdDev_d(DoubleFromString<double>(trimstr(sBuf.substr(68, 12))));
			}
			catch (...) {
				ss.str("");
				columnNo = 47;
				ss << "parse_sinex_stn(): Could not extract Y coordinate estimate from the record:  " << std::endl << "    " << sBuf << ".";
				throw std::runtime_error(ss.str());
			}
			continue;
		}

		// Capture the Z coordinate
		if (sBuf.substr(7, 4) == "STAZ")
		{
			try {
				stn_ptr->SetZAxis_d(DoubleFromString<double>(trimstr(sBuf.substr(47, 21))));
				stn_ptr->SetZAxisStdDev_d(DoubleFromString<double>(trimstr(sBuf.substr(68, 12))));
			}
			catch (...) {
				ss.str("");
				columnNo = 47;
				ss << "parse_sinex_stn(): Could not extract Z coordinate estimate from the record:  " << std::endl << "    " << sBuf << ".";
				throw std::runtime_error(ss.str());
			}				

			// Add this station to the vStations vector
			(*stnCount)++;
			file_rec++;
			vStations->push_back(stn_ptr);
			parsestn_tally.addstation("FFF");
		}			
	}
	
	// override header epoch (which is really only the start epoch of the data used) with the epoch from
	// the estimates. As per the SINEX standard, this is the epoch at which the estimated parameters are valid.
	if (!vStations->empty())
		datum.SetEpoch(vStations->at(0)->GetEpoch());

	// Do we need to remove sites in the case of unwanted discontinuities?
	if (containsDiscontinuities_ && !applyDiscontinuities_)
	{
		// Go in reverse direction so that stnCount can be used
		for (ptrdiff_t i = static_cast<ptrdiff_t>(*stnCount); i>0; --i)
		{
			if (!siteOccurrence_.at(i-1).last_occurrence)
			{
				vStations->erase(vStations->begin() + (i - 1));
				parsestn_tally.removestation("FFF");
			}
		}
		(*stnCount) = static_cast<UINT32>(vStations->size());
	}
}
	

void DnaIoSnx::ParseSinexMsr(std::ifstream** snx_file, const char* sinexRec, vdnaStnPtr* vStations, vdnaMsrPtr* vMeasurements, PUINT32 clusterID, PUINT32 msrCount,
					MsrTally& parsemsr_tally, CDnaDatum& datum, UINT32& lineNo)
{
	std::stringstream ss;

	if (applyDiscontinuities_ && !solutionEpochsRead_)
	{
		ss.str("");
		ss << "parse_sinex_msr(): Cannot apply discontinuities to the measurements" << std::endl <<
			"    if the station epochs have not been loaded beforehand.  To rectify this problem," << std::endl <<
			"    reformat the SINEX file so that the +SOLUTION/EPOCHS block appears before the +SOLUTION/ESTIMATE block.";
			throw std::runtime_error(ss.str());
	}

	std::string sBuf(sinexRec), tmp;
	parsemsr_tally.Y = static_cast<UINT32>(vStations->size() * 3);
	(*msrCount) = static_cast<UINT32>(vStations->size() * 3);

	UINT32 i, count, param1, param2, dimension(static_cast<UINT32>(vStations->size() * 3));
	double dparam[3];

	UINT32 dimension_all = dimension;

	if (containsDiscontinuities_)
		dimension_all = static_cast<UINT32>(siteOccurrence_.size() * 3);

	if (containsVelocities_)
		dimension_all *= 2;

	// If this file contains discontinuity sites and the user wishes to retain them,
	// then the size of the variance matrix will be the size of the stn discontinuity
	// vector, which was parsed in parse_sinex_discontinuities(...)
	if (containsDiscontinuities_ && applyDiscontinuities_)
		dimension = static_cast<UINT32>(siteOccurrence_.size() * 3);

	// default covariance matrix assumes 
	matrix_2d covariance_all(dimension_all, dimension_all);
	char cBuf[MAX_RECORD_LENGTH];
	char* p;

	// Get all covariances (whether discontinuities or velocities exist or not) and store
	// in covariance_all.  The unwanted elements will be stripped later
	while (sBuf.at(0) != '-')
	{
		lineNo++;
		(*snx_file)->getline(cBuf, MAX_RECORD_LENGTH);
		
		// Comments
		if (cBuf[0] == '*')
			continue;
		
		if (cBuf[0] == '-')
			break;

		p = cBuf;
		while (*p == ' ')
			p++;

		if ((count = GetFields(p, ' ', true, "ddfff", &param1, &param2, &dparam[0], &dparam[1], &dparam[2])) < 3)
		{
			ss.str("");
			ss << "parse_sinex_msr: Failed to read covariance elements from the record  " << cBuf << ".";
			throw std::runtime_error(ss.str());
		}

		// if the number of covariance elements has reached the number of stations, break out.
		if (param1 > dimension_all)
			break;

		// Fill covariance matrix.
		for (i=0; i<(count-2); i++)
			covariance_all.put(param1-1, param2-1+i, dparam[i]);			
	}

	// Remember, DynAdjust requires the upper diagonal part
	if (sinexRec[26] == LOWER_TRIANGLE)
		covariance_all.fillupper();

	matrix_2d covariance(dimension, dimension);

	UINT32 stn_discont;
	
	// At this point, velocities and/or discontinuities may be present in covariance_all.
	// Hence, the purpose of this block is to remove:
	//	1. velocity elements (DynAdjust doesn't handle velocities yet)
	//  2. discontinuities if they exist and the user doesn't want them (on account of not 
	//     providing a discontinuities file)
	UINT32 row_v, col_v, row, col;
	if (containsVelocities_ || (containsDiscontinuities_ && !applyDiscontinuities_))
	{
		for (row_v=0, row=0; row_v<dimension_all; ++row_v)
		{
			// Do we need to skip this site (in the case of unwanted discontinuities?
			if (containsDiscontinuities_ && !applyDiscontinuities_)
			{
				stn_discont = (UINT32)floor(row_v / 3.);
				if (containsVelocities_)
					stn_discont /= 2;

				// Is this the last occurrence?  If not, skip it
				if (!siteOccurrence_.at(stn_discont).last_occurrence)
				{
					// skip to z component
					row_v += 2;
					// If this site also includes velocities, then skip to 
					// the z component of the velocity
					if (containsVelocities_)
						row_v += 3;
#ifdef _MSDEBUG
					//TRACE("Skip to row:   %d\n", row_v+1);
#endif
					continue;
				}
#ifdef _MSDEBUG
				//else
					//TRACE("Keep this row: %d\n", row_v);
#endif
			}

			// At this point, unnecessary rows (discontinuity sites with or without velocities) have been skipped
			// Now copy the data
			// X row
			for (col_v=row_v, col=row; col_v<dimension_all; ++col_v)
			{
				// Do we need to skip this site (in the case of unwanted discontinuities?
				if (containsDiscontinuities_ && !applyDiscontinuities_)
				{
					stn_discont = (UINT32)floor((col_v+1) / 3.);
					if (containsVelocities_)
						stn_discont /= 2;
					
					// Is this the last occurrence?  If not, skip it
					if (!siteOccurrence_.at(stn_discont).last_occurrence)
					{
						// skip to z component
						col_v += 2;
						// If this site also includes velocities, then skip to 
						// the z component of the velocity
						if (containsVelocities_)
							col_v += 3;
						continue;
					}
				}

				covariance.put(row, col++, covariance_all.get(row_v, col_v++));
				covariance.put(row, col++, covariance_all.get(row_v, col_v++));
				covariance.put(row, col++, covariance_all.get(row_v, col_v));

				// skip velocity elements
				if (containsVelocities_)
					col_v += 3;
			}

			row_v++;
			row++;


			// Y row
			for (col_v=row_v, col=row; col_v<dimension_all; ++col_v)
			{
				// Do we need to skip this site (in the case of unwanted discontinuities?
				if (containsDiscontinuities_ && !applyDiscontinuities_)
				{
					stn_discont = (UINT32)floor((col_v+1) / 3.);
					if (containsVelocities_)
						stn_discont /= 2;

					// Is this the last occurrence?  If not, skip it
					if (!siteOccurrence_.at(stn_discont).last_occurrence)
					{
						// skip to z component
						col_v += 2;
						// If this site also includes velocities, then skip to 
						// the z component of the velocity
						if (containsVelocities_)
							col_v += 3;
						continue;
					}
				}

				covariance.put(row, col++, covariance_all.get(row_v, col_v++));
				covariance.put(row, col++, covariance_all.get(row_v, col_v));
				if (col-row > 2)
					covariance.put(row, col++, covariance_all.get(row_v, ++col_v));
				
				// skip velocity elements
				if (containsVelocities_)
					col_v += 3;
			}

			row_v++;
			row++;

			// Z row
			for (col_v=row_v, col=row; col_v<dimension_all; ++col_v)
			{
				// Do we need to skip this site (in the case of unwanted discontinuities?
				if (containsDiscontinuities_ && !applyDiscontinuities_)
				{
					stn_discont = (UINT32)floor((col_v+1) / 3.);
					if (containsVelocities_)
						stn_discont /= 2;

					if (stn_discont < siteOccurrence_.size())
					{
						// Is this the last occurrence?  If not, skip it
						if (!siteOccurrence_.at(stn_discont).last_occurrence)
						{
							// skip to z component
							col_v += 2;
							// If this site also includes velocities, then skip to 
							// the z component of the velocity
							if (containsVelocities_)
								col_v += 3;
							continue;
						}
					}
				}

				covariance.put(row, col++, covariance_all.get(row_v, col_v));
				if (col-row > 1)
					covariance.put(row, col++, covariance_all.get(row_v, ++col_v));
				if (col-row > 2)
					covariance.put(row, col++, covariance_all.get(row_v, ++col_v));
			
				// skip velocity elements
				if (containsVelocities_)
					col_v += 3;
			}

			// skip velocity elements
			if (containsVelocities_)
				row_v += 3;
			row++;
		}
	}
	else
		covariance = covariance_all;

	dnaMsrPtr dnaGpsPoint;
	dnaMsrPtr dnaGpsPointCluster;
	dnaCovariancePtr dnaPointCovariance;

	dnaGpsPointCluster.reset(new CDnaGpsPointCluster(++(*clusterID), datum.GetName(), datum.GetEpoch_s()));

	dnaGpsPointCluster->SetFirst(vStations->at(0)->GetName());
	dnaGpsPointCluster->SetTarget("");
	dnaGpsPointCluster->SetIgnore(false);
	dnaGpsPointCluster->SetTotal(static_cast<UINT32>(vStations->size()));
	dnaGpsPointCluster->SetCoordType(XYZ_type);
	dnaGpsPointCluster->SetPscale(1.);
	dnaGpsPointCluster->SetLscale(1.);
	dnaGpsPointCluster->SetHscale(1.);
	dnaGpsPointCluster->SetVscale(1.);

	dnaGpsPointCluster->SetReferenceFrame(datum.GetName());
	dnaGpsPointCluster->SetEpsg(datum.GetEpsgCode_s());
	dnaGpsPointCluster->SetEpoch(vStations->at(0)->GetEpoch());

	UINT32 cov_count, ci;

	for (UINT32 k, c, p(0); p<vStations->size(); ++p)
	{
		dnaGpsPoint.reset(new CDnaGpsPoint);
		dnaGpsPoint->SetType("Y");
	
		dnaGpsPoint->SetIgnore(dnaGpsPointCluster->GetIgnore());
		
		dnaGpsPoint->SetFirst(vStations->at(p)->GetName());
		dnaGpsPoint->SetTarget("");
		dnaGpsPoint->SetCoordType("XYZ");
		dnaGpsPoint->SetRecordedTotal(static_cast<UINT32>(vStations->size()));

		dnaGpsPoint->SetReferenceFrame(datum.GetName());
		dnaGpsPoint->SetEpsg(datum.GetEpsgCode_s());
		dnaGpsPoint->SetEpoch(vStations->at(p)->GetEpoch());

		dnaGpsPoint->SetPscale(dnaGpsPointCluster->GetPscale());
		dnaGpsPoint->SetLscale(dnaGpsPointCluster->GetLscale());
		dnaGpsPoint->SetHscale(dnaGpsPointCluster->GetHscale());
		dnaGpsPoint->SetVscale(dnaGpsPointCluster->GetVscale());
		dnaGpsPoint->SetClusterID(dnaGpsPointCluster->GetClusterID());

		dnaGpsPoint->SetXAxis(vStations->at(p)->GetXAxis());
		dnaGpsPoint->SetYAxis(vStations->at(p)->GetYAxis());
		dnaGpsPoint->SetZAxis(vStations->at(p)->GetZAxis());

		k = p * 3;

		dnaGpsPoint->SetSigmaXX(covariance.get(k,k));
		dnaGpsPoint->SetSigmaXY(covariance.get(k,k+1));
		dnaGpsPoint->SetSigmaXZ(covariance.get(k,k+2));
		dnaGpsPoint->SetSigmaYY(covariance.get(k+1,k+1));
		dnaGpsPoint->SetSigmaYZ(covariance.get(k+1,k+2));
		dnaGpsPoint->SetSigmaZZ(covariance.get(k+2,k+2));

		cov_count = static_cast<UINT32>(vStations->size()) - p - 1;
		
		// add covariances
		for (c=0; c<cov_count; ++c)
		{
			ci = 3 + k + (c * 3);

			dnaPointCovariance.reset(new CDnaCovariance);
			dnaPointCovariance->SetType("Y");

			dnaPointCovariance->SetM11(covariance.get(k, ci));
			dnaPointCovariance->SetM12(covariance.get(k, ci+1));
			dnaPointCovariance->SetM13(covariance.get(k, ci+2));
			dnaPointCovariance->SetM21(covariance.get(k+1, ci));
			dnaPointCovariance->SetM22(covariance.get(k+1, ci+1));
			dnaPointCovariance->SetM23(covariance.get(k+1, ci+2));
			dnaPointCovariance->SetM31(covariance.get(k+2, ci));
			dnaPointCovariance->SetM32(covariance.get(k+2, ci+1));
			dnaPointCovariance->SetM33(covariance.get(k+2, ci+2));

			dnaPointCovariance->SetClusterID(dnaGpsPointCluster->GetClusterID());

			dnaGpsPoint->AddPointCovariance(dnaPointCovariance.get());
		}
		dnaGpsPointCluster->AddGpsPoint(dnaGpsPoint.get());
	}
	vMeasurements->push_back(dnaGpsPointCluster);
}


} // dnaiostreams
} // dynadjust
