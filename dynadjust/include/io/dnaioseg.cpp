//============================================================================
// Name         : dnaioseg.cpp
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
// Description  : DynAdjust segmentation file io operations
//============================================================================

#include <include/io/dnaioseg.hpp>
#include <include/functions/dnaiostreamfuncs.hpp>
#include <include/functions/dnatemplatestnmsrfuncs.hpp>

namespace dynadjust { 
namespace iostreams {

void dna_io_seg::load_seg_file_header_f(const string& seg_filename, UINT32& blockCount, 
						UINT32& blockThreshold, UINT32& minInnerStns) 
{	
	std::ifstream seg_file;
	stringstream ss_err;

	try {
		
		// open stations seg file.  Throws runtime_error on failure.
		ss_err << "load_seg_file_f(): An error was encountered when opening " << seg_filename << "." << endl;
		file_opener(seg_file, seg_filename, ios::in, ascii, true);
			
		ss_err.str("");
		ss_err << "load_seg_file_f(): An error was encountered when reading " << seg_filename << "." << endl;
		load_seg_file_header(seg_filename, seg_file, blockCount, 
			blockThreshold, minInnerStns);

		seg_file.close();
	}
	catch (const runtime_error& e) {
		ss_err << e.what();
		throw boost::enable_current_exception(runtime_error(ss_err.str()));
	}
	catch (...) {
		throw boost::enable_current_exception(runtime_error(ss_err.str()));
	}
}
	

void dna_io_seg::load_seg_file_header(const string& seg_filename, istream& seg_file, UINT32& blockCount, 
	UINT32& blockThreshold, UINT32& minInnerStns) 
{	

	stringstream ss_err;	
	ss_err << "load_seg_file_header(): An error was encountered when reading from " << seg_filename << "." << endl;

	char line[PRINT_LINE_LENGTH];
	memset(line, '\0', PRINT_LINE_LENGTH);
		
	try {
		// read header line
		seg_file.ignore(PRINT_LINE_LENGTH, '\n');		// ------------------------
		seg_file.getline(line, PRINT_LINE_LENGTH);		//            DYNADJUST SEGME...
		seg_file.ignore(PRINT_LINE_LENGTH, '\n');		// 
		seg_file.ignore(PRINT_LINE_LENGTH, '\n');		// Version
		seg_file.ignore(PRINT_LINE_LENGTH, '\n');		// Build
		seg_file.ignore(PRINT_LINE_LENGTH, '\n');		// File created
		seg_file.ignore(PRINT_LINE_LENGTH, '\n');		// File name
		seg_file.ignore(PRINT_LINE_LENGTH, '\n');		// 
		seg_file.ignore(PRINT_LINE_LENGTH, '\n');		// Command line arguments
		seg_file.ignore(PRINT_LINE_LENGTH, '\n');		// 
		seg_file.ignore(PRINT_LINE_LENGTH, '\n');		// Stations file
		seg_file.ignore(PRINT_LINE_LENGTH, '\n');		// Measurements file
		seg_file.ignore(PRINT_LINE_LENGTH, '\n');		// 

		seg_file.getline(line, PRINT_LINE_LENGTH);		// Minimum inner stations
		minInnerStns = strtoul(line+PRINT_VAR_PAD, NULL, 0);
		
		seg_file.getline(line, PRINT_LINE_LENGTH);		// Block size threshold
		blockThreshold = strtoul(line+PRINT_VAR_PAD, NULL, 0);
		
		seg_file.ignore(PRINT_LINE_LENGTH, '\n');		// Starting stations
		seg_file.getline(line, PRINT_LINE_LENGTH);		// ???
		while (line[0]!= '-')
			seg_file.getline(line, PRINT_LINE_LENGTH);	// ------------------------
		
		seg_file.ignore(PRINT_LINE_LENGTH, '\n');		// 		
		seg_file.ignore(PRINT_LINE_LENGTH, '\n');		// SEGMENTATION SUMMARY
		seg_file.ignore(PRINT_LINE_LENGTH, '\n');		// 		
		
		// read block count
		seg_file.getline(line, PRINT_LINE_LENGTH);		// No. blocks produced 
		blockCount = strtoul(line+PRINT_VAR_PAD, NULL, 0);
	}
	catch (const ios_base::failure& f) {
		ss_err << f.what();
		throw boost::enable_current_exception(runtime_error(ss_err.str()));
	}
	catch (const runtime_error& e) {
		ss_err << e.what();
		throw boost::enable_current_exception(runtime_error(ss_err.str()));
	}
	catch (...) {
		throw boost::enable_current_exception(runtime_error(ss_err.str()));
	}
}

void dna_io_seg::load_seg_file(const string& seg_filename, UINT32& blockCount, 
	UINT32& blockThreshold, UINT32& minInnerStns,
	vvUINT32& v_ISL, vvUINT32& v_JSL, vvUINT32& v_CML,
	bool loadMetrics,
	pvmsr_t bmsBinaryRecords, pvUINT32 v_measurementCount, 
	pvUINT32 v_unknownsCount, pvUINT32 v_ContiguousNetList,
	pvUINT32 v_parameterStationCount) 
{	
	std::ifstream seg_file;
	stringstream ss_err;
	ss_err << "load_seg_file(): An error was encountered when opening " << seg_filename << "." << endl;

	try {
		// open stations seg file.  Throws runtime_error on failure.
		file_opener(seg_file, seg_filename, ios::in, ascii, true);
	}
	catch (const runtime_error& e) {
		ss_err << e.what();
		throw boost::enable_current_exception(runtime_error(ss_err.str()));
	}
	catch (...) {
		throw boost::enable_current_exception(runtime_error(ss_err.str()));
	}
	
	ss_err.str("");
	ss_err << "load_seg_file(): An error was encountered when reading from " << seg_filename << "." << endl;

	UINT32 b, blk, c, blkCount(0), i, j, m, 
		netID(0), jslCount(0), islCount(0), msrCount(0), stnCount(0);	
	
	char line[PRINT_LINE_LENGTH], format_spec[41], format_spec_netid[6], 
		format_spec_junct[6], format_spec_inner[6], format_spec_measr[6];
		
	try {
		
		load_seg_file_header(seg_filename, seg_file, blockCount, 
			blockThreshold, minInnerStns);
		
		// Resize vectors based upon block count
		v_ISL.resize(blockCount);
		v_JSL.resize(blockCount);
		v_CML.resize(blockCount);

		if (loadMetrics)
		{
			v_ContiguousNetList->resize(blockCount);
			v_measurementCount->resize(blockCount);
			v_unknownsCount->resize(blockCount);
			v_parameterStationCount->resize(blockCount);
		}

		// skip header info
		seg_file.ignore(PRINT_LINE_LENGTH, '\n');		// ------------------------
		seg_file.ignore(PRINT_LINE_LENGTH, '\n');		//   Block   Junction stns  Inner stns  Measurements  Total stns  

		sprintf(format_spec, "%%%d%s%%%d%s%%%d%s%%%d%s%%%d%s%%%d%s", BLOCK, "lu", NETID, "lu", JUNCT, "lu", INNER, "lu", MEASR, "lu", MEASR, "lu");
		sprintf(format_spec_netid, "%%%d%s", NETID, "lu");
		sprintf(format_spec_junct, "%%%d%s", JUNCT, "lu");
		sprintf(format_spec_inner, "%%%d%s", INNER, "lu");
		sprintf(format_spec_measr, "%%%d%s", MEASR, "lu");
		
		
		UINT32 t;
		// read block sizes
		for (t=0; t<blockCount; ++t)
		{
			seg_file.getline(line, PRINT_LINE_LENGTH);

			if (strncmp(line, "--------------------", 20) == 0)
				throw boost::enable_current_exception(runtime_error("  Segmentation file is corrupt."));
		
			if (sscanf(line, format_spec, &blkCount, &netID, &jslCount, &islCount, &msrCount, &stnCount) < 6)
				throw boost::enable_current_exception(runtime_error("  Segmentation file is corrupt."));

			if (stnCount != islCount + jslCount)
				throw boost::enable_current_exception(runtime_error("  Segmentation file is corrupt."));

			v_JSL.at(t) = vUINT32(jslCount);
			v_ISL.at(t) = vUINT32(islCount);
			v_CML.at(t) = vUINT32(msrCount);
		
			if (loadMetrics)
			{
				v_ContiguousNetList->at(t) = netID;
				v_measurementCount->at(t) = v_unknownsCount->at(t) = 0;
				v_parameterStationCount->at(t) = stnCount;
			}
		}
		

		// skip header info
		seg_file.ignore(PRINT_LINE_LENGTH, '\n');			// ------------------------
		seg_file.ignore(PRINT_LINE_LENGTH, '\n');			// 
		seg_file.ignore(PRINT_LINE_LENGTH, '\n');			// Individual block data (ASL and AML indices):
		seg_file.ignore(PRINT_LINE_LENGTH, '\n');			// ------------------------
		
		// initialise block ID
		netID = 999999;

		// read block data
		for (b=0; b<blockCount; ++b)
		{
			seg_file.ignore(PRINT_LINE_LENGTH, '\n');		// 
			seg_file.getline(line, PRINT_LINE_LENGTH);		// Block #
			sscanf(line+5, "%ud", &blk);
			if (b+1 != blk)
				throw boost::enable_current_exception(runtime_error("load_seg_file(): segmentation file is corrupt."));
			
			seg_file.ignore(PRINT_LINE_LENGTH, '\n');		// ------------------------
			seg_file.ignore(PRINT_LINE_LENGTH, '\n');		// Junction stns:
			seg_file.ignore(PRINT_LINE_LENGTH, '\n');		// Inner stns:   
			seg_file.ignore(PRINT_LINE_LENGTH, '\n');		// Measurements:
			seg_file.ignore(PRINT_LINE_LENGTH, '\n');		// Total stns:  
			seg_file.ignore(PRINT_LINE_LENGTH, '\n');		// 
			seg_file.ignore(PRINT_LINE_LENGTH, '\n');		// Inners      Junctions      Measurements      
			seg_file.ignore(PRINT_LINE_LENGTH, '\n');		// ------------------------

			seg_file.getline(line, PRINT_LINE_LENGTH);		// Ok, now the data!
			
			jslCount = static_cast<UINT32>(v_JSL.at(b).size());
			islCount = static_cast<UINT32>(v_ISL.at(b).size());
			msrCount = static_cast<UINT32>(v_CML.at(b).size());
			
			c = 0;
			while (strncmp(line, "--------------------", 20) != 0)
			{
				if (c < islCount)
				{
					sscanf(line, format_spec_inner, &i);
					v_ISL.at(b).at(c) = i;
					if (loadMetrics)
						v_unknownsCount->at(b) += 3;
				}
				if (c < jslCount)
				{
					sscanf(line+12, format_spec_junct, &j);
					v_JSL.at(b).at(c) = j;
					if (loadMetrics)
						v_unknownsCount->at(b) += 3;
				}
				if (c < msrCount)
				{
					sscanf(line+24, format_spec_measr, &m);
					v_CML.at(b).at(c) = m;

					if (loadMetrics)
					{
						// Calculate number of 'measurements' in this measurement
						// No need to check for ignored measurement since the segmentation 
						// algorithm only includes measurements with the ignore flag cleared
						switch (bmsBinaryRecords->at(m).measType)
						{
						case 'G':	// GPS Baseline
							v_measurementCount->at(b) += 3;		// three measurements (x, y, z)
							break;
						case 'X':	// GPS Baseline cluster
						case 'Y':	// GPS point cluster
							v_measurementCount->at(b) += bmsBinaryRecords->at(m).vectorCount1 * 3;
							break;
						case 'D':	// Direction set	
							v_measurementCount->at(b) += bmsBinaryRecords->at(m).vectorCount1 - 1;
							break;
						case 'A':	// Horizontal angle
						case 'B':	// Geodetic azimuth
						case 'C':	// Chord dist
						case 'E':	// Ellipsoid arc
						case 'H':	// Orthometric height
						case 'I':	// Astronomic latitude
						case 'J':	// Astronomic longitude
						case 'K':	// Astronomic azimuth
						case 'L':	// Level difference
						case 'P':	// Geodetic latitude
						case 'Q':	// Geodetic longitude
						case 'R':	// Ellipsoidal height
						case 'M':	// MSL arc
						case 'S':	// Slope distance
						case 'V':	// Zenith angle
						case 'Z':	// Vertical angle
						default:
							v_measurementCount->at(b)++;		// single measurement quantity
							break;
						}
					}
				}
				c++;
				seg_file.getline(line, PRINT_LINE_LENGTH);		// get more data
			}
		}
	}
	catch (const ios_base::failure& f) {
		ss_err << f.what();
		throw boost::enable_current_exception(runtime_error(ss_err.str()));
	}
	catch (const runtime_error& e) {
		ss_err << e.what();
		throw boost::enable_current_exception(runtime_error(ss_err.str()));
	}
	catch (...) {
		throw boost::enable_current_exception(runtime_error(ss_err.str()));
	}

	seg_file.close();

}

void dna_io_seg::build_free_stn_availability(vASL& assocStnList, v_freestn_pair& freeStnList)
{
	freeStnList.clear();

	// Resize free station availability list.  Note - the
	// vfreeStnAvailability_ constructor sets all stations to 'free'
	freeStnList.resize(assocStnList.size());

	UINT32 stn_id(0);
	for_each(
		freeStnList.begin(), 
		freeStnList.end(), 
		[&stn_id, &assocStnList](freestn_pair& stn){
			// consume invalid stations
			if (assocStnList.at(stn_id).IsInvalid())
				stn.consume();
			// Set ID
			stn.stn_index = stn_id++;
		}
	);
}

void dna_io_seg::create_stn_appearance_list(vv_stn_appear& stnAppearance,
	const vvUINT32& paramStationList,
	vASL& assocStnList)
{
	// build free station list
	v_freestn_pair freeStnList;
	build_free_stn_availability(assocStnList, freeStnList);

	// Make a copy for the reverse pass
	v_freestn_pair revStnList(freeStnList);

	UINT32 block, blockCount(static_cast<UINT32>(stnAppearance.size()));

	it_vstn_appear _it_appear;
	it_vUINT32_const _it_const;
	
	// 1. Create the station 'appearances' map for the forward pass
	for (block=0; block<blockCount; ++block)
	{
		for (_it_const=paramStationList.at(block).begin(),
			_it_appear=stnAppearance.at(block).begin();
			_it_appear!=stnAppearance.at(block).end();
			++_it_const, ++_it_appear)
		{
			// assign id
			_it_appear->set_id(*_it_const);

			if (freeStnList.at(_it_appear->station_id).isfree())
			{
				// assign first appearance flag and consume free station
				_it_appear->first_fwd();

				freeStnList.at(_it_appear->station_id).consume();
			}
		}
	}

	UINT32 blockR(blockCount- 1);
	
	// 2. Create the station 'appearances' map for the reverse pass
	for (block=0; block<blockCount; ++block, --blockR)
	{
		for (_it_appear=stnAppearance.at(blockR).begin();
			_it_appear!=stnAppearance.at(blockR).end(); 
			++_it_appear)
		{			
			if (revStnList.at(_it_appear->station_id).isfree())
			{
				// assign first appearance flag and consume free station
				_it_appear->first_rev();

				revStnList.at(_it_appear->station_id).consume();
			}
		}
	}
}

void dna_io_seg::write_seg_block(ostream &os, 
	const vUINT32& vISL, const vUINT32& vJSL, const vUINT32& vCML, 
	const UINT32& currentBlock, 
	const vstn_t* bstBinaryRecords, const vmsr_t* bmsBinaryRecords, 
	bool PRINT_NAMES)
{
	char dash;

	os << endl << "Block " << currentBlock << endl;
	for (dash=INNER+JUNCT+MEASR+PAD; dash>0; dash--)
		os << "-";
	os << endl;

	//os << setw(JUNCT) << left << "Block:" << setw(BLOCK) << right << currentBlock_ << endl;
	os << setw(JUNCT) << left << "Junction stns:" << setw(BLOCK) << vJSL.size() << endl;
	os << setw(JUNCT) << left << "Inner stns:" << setw(BLOCK) << vISL.size() << endl;
	os << setw(JUNCT) << left << "Measurements:" << setw(BLOCK) << vCML.size() << endl;
	os << setw(JUNCT) << left << "Total stns:" << setw(BLOCK) << (vJSL.size() + vISL.size()) << endl << endl;

	it_vUINT32_const _it_isl(vISL.begin());
	it_vUINT32_const _it_jsl(vJSL.begin());
	it_vUINT32_const _it_msr(vCML.begin());
	
	os << setw(INNER) << left << "Inner stns" << setw(JUNCT) << left << "Junction stns" << 
		setw(MEASR) << "Measurements" << setw(PAD) << "Type" << endl;
	
	for (dash=INNER+JUNCT+MEASR+PAD; dash>0; dash--)
		os << "-";
	os << endl;
	ostringstream tmp;
	
	vUINT32 msrStations;

	while (1)
	{
		tmp.str("");
		if (_it_isl != vISL.end())
		{
			if (PRINT_NAMES) {
				os << setw(INNER) << left << ((ostringstream&)(tmp << bstBinaryRecords->at(*_it_isl).stationName)).str();
				tmp.str("");
			}
			else
				os << setw(INNER) << left << *_it_isl;
			_it_isl++;
		}
		else
			os << setw(INNER) << left << " ";
		if (_it_jsl != vJSL.end())
		{
			if (PRINT_NAMES) {
				os << setw(JUNCT) << left << ((ostringstream&)(tmp << bstBinaryRecords->at(*_it_jsl).stationName)).str();
				tmp.str("");
			}			
			else
				os << setw(JUNCT) << left << *_it_jsl;
			_it_jsl++;
		}
		else
			os << setw(JUNCT) << left << " ";
		if (_it_msr != vCML.end())
		{
			if (PRINT_NAMES) {

				// Get all the stations associated with this measurement
				GetMsrStations(*bmsBinaryRecords, *_it_msr, msrStations);

				os << left << bmsBinaryRecords->at(*_it_msr).measType << " (";
				
				for_each(msrStations.begin(), msrStations.end(),
					[&os, &bstBinaryRecords](UINT32 stn){
					os << bstBinaryRecords->at(stn).stationName << " ";
				});

				os << ") ";
				
				if (bmsBinaryRecords->at(*_it_msr).clusterID > 0)
					os << "cluster " << bmsBinaryRecords->at(*_it_msr).clusterID;
				
				os << endl;
			}
			else
				os << setw(MEASR) << left << *_it_msr << setw(PAD) << left << bmsBinaryRecords->at(*_it_msr).measType << endl;
			_it_msr++;
		}
		else
			os << endl;

		// continue until all lists are exhausted
		if (_it_msr==vCML.end() &&
			_it_jsl==vJSL.end() &&
			_it_isl==vISL.end())
			break;
	}

	for (dash=INNER+JUNCT+MEASR+PAD; dash>0; dash--)
		os << "-";
	os << endl;

}

void dna_io_seg::write_seg_file(const string& seg_filename, const string& bst_filename, const string& bms_filename,
	const UINT32& min_inner_stns, const UINT32& max_block_stns,
	const string& seg_starting_stns, const vstring& vinitialStns,
	const string& command_line_arguments,
	vvUINT32& v_ISL, vvUINT32& v_JSL, vvUINT32& v_CML,
	vUINT32& v_ContiguousNetList, const pvstn_t bstBinaryRecords, const pvmsr_t bmsBinaryRecords)
{
	std::ofstream seg_file;
	stringstream ss;
	ss << "write_seg_file(): An error was encountered when opening " << seg_filename << "." << endl;

	try {
		// Create segmentation file.  Throws runtime_error on failure.
		file_opener(seg_file, seg_filename);
	}
	catch (const runtime_error& e) {
		ss << e.what();
		throw boost::enable_current_exception(runtime_error(ss.str()));
	}
	catch (...) {
		throw boost::enable_current_exception(runtime_error(ss.str()));
	}

	// Print formatted header
	print_file_header(seg_file, "DYNADJUST SEGMENTATION OUTPUT FILE");
	
	seg_file << setw(PRINT_VAR_PAD) << left << "File name:" << system_complete(seg_filename).string() << endl << endl;

	seg_file << setw(PRINT_VAR_PAD) << left << "Command line arguments: ";
	seg_file << command_line_arguments << endl << endl;

	seg_file << setw(PRINT_VAR_PAD) << left << "Stations file:" << system_complete(bst_filename).string() << endl;
	seg_file << setw(PRINT_VAR_PAD) << left << "Measurements file:" << system_complete(bms_filename).string() << endl;	

	UINT32 b = 1;
	
	seg_file << endl << setw(PRINT_VAR_PAD) << left << "Minimum inner stations" << min_inner_stns << endl;
	seg_file << setw(PRINT_VAR_PAD) << left << "Block size threshold" << max_block_stns << endl;
	if (seg_starting_stns.size() == 1)
		seg_file << setw(PRINT_VAR_PAD) << left << "Starting station" << seg_starting_stns.at(0) << endl;
	else
	{
		string s("Starting station(s)");
		for (b=0; b<vinitialStns.size(); ++b)
		{
			if (b > 0)
				 seg_file << endl;
			seg_file << setw(PRINT_VAR_PAD) << left << s <<
				vinitialStns.at(b);
			s = " ";
		}
		seg_file << endl;
	}

	seg_file << OUTPUTLINE << endl << endl;
	seg_file << setw(PRINT_VAR_PAD) << left << "SEGMENTATION SUMMARY" << endl << endl;
	seg_file << setw(PRINT_VAR_PAD) << left << "No. blocks produced" << v_ISL.size() << endl;
	
	char dash;
	for (dash=BLOCK+NETID+INNER+JUNCT+TOTAL+MEASR; dash>2; dash--)
		seg_file << "-";
	
	seg_file << endl;
	seg_file << setw(BLOCK) << left << "  Block" << setw(NETID) << left << "Network ID" 
		<< setw(JUNCT) << left << "Junction stns" << setw(INNER) << left << "Inner stns" 
		<< setw(MEASR) << left << "Measurements" << setw(TOTAL) << left << "Total stns"  << endl;

	it_vUINT32 _it_net(v_ContiguousNetList.begin());
	it_vvUINT32 _it_isl(v_ISL.begin()), _it_jsl(v_JSL.begin()), _it_cml(v_CML.begin());
	
	b=1;
	
	for (; _it_isl!=v_ISL.end(); _it_isl++)
	{
		// block
		seg_file << "  " << setw(BLOCK-2) << left << b;

		// network id
		seg_file << setw(NETID) << left << *_it_net;

		// junction stns
		if (_it_jsl!=v_JSL.end())
			seg_file << setw(JUNCT) << left << _it_jsl->size();
		else
			seg_file << setw(JUNCT) << " ";
		
		// inner stns
		seg_file << setw(INNER) << left << _it_isl->size();
		
		// total measurements
		if (_it_cml!=v_CML.end())
			seg_file << setw(MEASR) << left << _it_cml->size();
		else
			seg_file << setw(MEASR) << " ";

		// total stns
		if (_it_jsl!=v_JSL.end())
			seg_file << setw(TOTAL) << left << (_it_isl->size() + _it_jsl->size());
		else
			seg_file << setw(TOTAL) << left << _it_isl->size();
		
		seg_file << endl;

		_it_jsl++;
		_it_cml++;
		_it_net++;
		b++;
	}

	for (dash=BLOCK+NETID+INNER+JUNCT+TOTAL+MEASR; dash>2; dash--)
		seg_file << "-";
	seg_file << endl << endl << "INDIVIDUAL BLOCK DATA" << endl;
	for (dash=BLOCK+NETID+INNER+JUNCT+TOTAL+MEASR; dash>2; dash--)
		seg_file << "-";
	seg_file << endl;

	_it_isl = v_ISL.begin();
	_it_jsl = v_JSL.begin();
	_it_cml = v_CML.begin();
	UINT32 i(1);
	while (_it_isl!=v_ISL.end())
	{
		write_seg_block(seg_file, *_it_isl, *_it_jsl, *_it_cml, i++, bstBinaryRecords, bmsBinaryRecords);
		_it_isl++;
		_it_jsl++;
		_it_cml++;
	}

	seg_file << endl;

	seg_file.close();
}

void dna_io_seg::write_stn_appearance(const string& sap_filename, const v_stn_block_map& stnAppearance)
{
	std::ofstream sap_file;
	stringstream ss;
	ss << "write_sap_file(): An error was encountered when opening " << sap_filename << "." << endl;

	try {
		// Create segmentation file.  Throws runtime_error on failure.
		file_opener(sap_file, sap_filename, ios::out | ios::binary, binary);
	}
	catch (const runtime_error& e) {
		ss << e.what();
		throw boost::enable_current_exception(runtime_error(ss.str()));
	}
	catch (...) {
		throw boost::enable_current_exception(runtime_error(ss.str()));
	}

	// write the list
	for_each(
		stnAppearance.begin(), 
		stnAppearance.end(), 
		[this, &sap_file](const stn_block_map& stn){
			sap_file.write(reinterpret_cast<const char *>(&stn.block_no), sizeof(UINT32));
			sap_file.write(reinterpret_cast<const char *>(&stn.first_appearance_fwd), sizeof(bool));
			sap_file.write(reinterpret_cast<const char *>(&stn.first_appearance_rev), sizeof(bool));
			sap_file.write(reinterpret_cast<const char *>(&stn.valid_stn), sizeof(bool));
		}
	);

	sap_file.close();
}

void dna_io_seg::load_stn_appearance(const string& sap_filename, v_stn_block_map& stnAppearance)
{
	std::ifstream sap_file;
	stringstream ss;
	ss << "load_sap_file(): An error was encountered when opening " << sap_filename << "." << endl;

	try {
		// Create segmentation file.  Throws runtime_error on failure.
		file_opener(sap_file, sap_filename, ios::in | ios::binary, binary);
	}
	catch (const runtime_error& e) {
		ss << e.what();
		throw boost::enable_current_exception(runtime_error(ss.str()));
	}
	catch (...) {
		throw boost::enable_current_exception(runtime_error(ss.str()));
	}

	// write the list
	for_each(
		stnAppearance.begin(), 
		stnAppearance.end(), 
		[this, &sap_file](stn_block_map& stn){
			sap_file.read(reinterpret_cast<char *>(&stn.block_no), sizeof(UINT32));
			sap_file.read(reinterpret_cast<char *>(&stn.first_appearance_fwd), sizeof(bool));
			sap_file.read(reinterpret_cast<char *>(&stn.first_appearance_rev), sizeof(bool));
			sap_file.read(reinterpret_cast<char *>(&stn.valid_stn), sizeof(bool));
		}
	);

	sap_file.close();
}

} // dnaiostreams
} // dynadjust
