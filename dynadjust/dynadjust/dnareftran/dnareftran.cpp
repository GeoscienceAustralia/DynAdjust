//============================================================================
// Name         : dnareftran.cpp
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
// Description  : Reference Frame Transformation library
//============================================================================

#include <dynadjust/dnareftran/dnareftran.hpp>

namespace dynadjust { 
namespace referenceframe {

dna_reftran::dna_reftran()
	: databaseIDsLoaded_(false)
	, databaseIDsSet_(false)
{
#ifdef _MSC_VER
#if (_MSC_VER < 1900)
	{
		// this function is obsolete in MS VC++ 14.0, VS2015
		// Set scientific format to print two places for the exponent
		_set_output_format(_TWO_DIGIT_EXPONENT);
	}
#endif
#endif
	rft_file = 0;

	LoadWGS84FrameSubstitutions();
}

dna_reftran::dna_reftran(const project_settings& p, std::ofstream* f_out) 
	: databaseIDsLoaded_(false)
	, databaseIDsSet_(false)
{
	InitialiseSettings(p);
	// reftran log file pointer
	rft_file = f_out;

	LoadWGS84FrameSubstitutions();
}


dna_reftran::~dna_reftran()
{

}

void dna_reftran::TransformBinaryFiles(const std::string& bstFile, const std::string& bmsFile, const std::string& newFrame, const std::string& newEpoch)
{
	// TODO - Would it be faster to use memory mapping instead of reading
	// binary files into memory?  Not sure.
	// Either way, using memory mapped files would provide a greater capacity
	// to run extremely large adjustments.

	// 1. Load binary files
	// load the binary station file into memory
	LoadBinaryStationFile(bstFile);

	// Identify and apply any substitutions for WGS84 in the list of stations
	ApplyStationFrameSubstitutions();

	if (projectSettings_.r.plate_model_option == 1)
		IdentifyStationPlate();

	// load the binary measurement file into memory
	LoadBinaryMeasurementFile(bmsFile);

	// Identify and apply any substitutions for WGS84 in the list of measurements
	ApplyMeasurementFrameSubstitutions();

	datumTo_.SetDatumFromName(newFrame, newEpoch);

	// 2. Transform measurements first (because pre-transformed station 
	// coordinates are required)
	TransformMeasurementRecords(newFrame, newEpoch);

	// write the binary measurement file
	WriteBinaryMeasurementFile(bmsFile);

	// 3. Transform stations
	TransformStationRecords(newFrame, newEpoch);

	// write the binary station file
	WriteBinaryStationFile(bstFile);
}
	
// Obtain the two-character ID for the tectonic plate on which each station lies, and assign
// to the binary station record.  Uses boost 'point in polygon' search.
// Then, when plate motion model parameters are required for a site, the appropriate parameters
// can be selected using he two-character ID.
void dna_reftran::IdentifyStationPlate()
{
	it_vstn_t stn_it;

	dnaGeometryPolygon platePolygon;
	v_doubledouble_pair::reverse_iterator _it_boundary_points; 

	it_v_string_v_doubledouble_pair _it_plates;
	dnaGeometryPoint<double> point;

	//cout << std::endl;

	size_t plateCount(global_plates_.size());
	string_uint32_pair stnPlate;
	UINT32 p(0);

	vplateMap_.clear();
	vplateMap_.reserve(plateCount);	

	for (_it_plates=global_plates_.begin(); _it_plates!=global_plates_.end(); ++_it_plates)
	{
		boost::geometry::clear(platePolygon);

		// To avoid this step, could global_plates_ be changed from a vector of (double pairs) to
		// dnaGeometryPoints?  It's not a large burden to load these as-is.
		// Note that coordinate tuples for tectonic plate boundaries are loaded in 
		// anticlockwise order.
		for (_it_boundary_points=_it_plates->second.rbegin();
			_it_boundary_points!=_it_plates->second.rend();
			++_it_boundary_points)
		{
			point.set_east_long(_it_boundary_points->first);
			point.set_north_lat(_it_boundary_points->second);
			boost::geometry::append(platePolygon, point);
		}

		stnPlate.first = _it_plates->first;
		stnPlate.second = p++;
		vplateMap_.push_back(stnPlate);		

		//if (iequals(_it_plates->first, "AU"))
		//	std::cout << std::endl << _it_plates->first << std::endl << 
		//		boost::geometry::wkt(platePolygon) << std::endl;

		for (stn_it=bstBinaryRecords_.begin(); stn_it!=bstBinaryRecords_.end(); ++stn_it)
		{
			if (stn_it->plate[0] >= 'A' && stn_it->plate[0] <= 'Z')
				continue;
				
			point.set_east_long(radians_to_degrees_(stn_it->currentLongitude));
			point.set_north_lat(radians_to_degrees_(stn_it->currentLatitude));

			if (boost::geometry::within(point, platePolygon))
			{
				snprintf(stn_it->plate, sizeof(stn_it->plate), "%s", _it_plates->first.c_str());
				//cout << "Station " << stn_it->stationName << " is on plate " << _it_plates->first << std::endl;
			}	
		}
	}

	std::sort(vplateMap_.begin(), vplateMap_.end());
}

	
void dna_reftran::LoadTectonicPlateParameters(const std::string& pltfileName, const std::string& pmmfileName)
{
	dna_io_tpb tpb;
	std::stringstream ss;
	ss << "LoadTectonicPlateParameters(): An error was encountered when loading" << std::endl <<
		"  tectonic plate information." << std::endl;
	
	projectSettings_.r.plate_model_option = 1;

	try {
	
		// Load tectonic plate parameters.  Throws runtime_error on failure.
		tpb.load_tpp_file(pmmfileName, plate_motion_eulers_);
		std::sort(plate_motion_eulers_.begin(), plate_motion_eulers_.end());

		// Load tectonic plate polygons.  Throws runtime_error on failure.
		tpb.load_tpb_file(pltfileName, global_plates_);
		std::sort(global_plates_.begin(), global_plates_.end());		
	}
	catch (const std::runtime_error& e) {
		ss << e.what();
		throw std::runtime_error(ss.str());
	}

	std::string message;
	if (!tpb.validate_plate_files(global_plates_, plate_motion_eulers_, message))
	{
		ss << "         " << message << std::endl;
		throw std::runtime_error(ss.str());
	}

	try {
		CalculateRotations();
	}
	catch (const std::runtime_error& e) {
		ss << e.what();
		throw std::runtime_error(ss.str());
	}
}
	
void dna_reftran::CalculateRotations()
{
	plate_motion_cartesians_.clear();

	plate_motion_cartesian pmm;
	double r_lat, r_lon, r_rot;
	UINT32 j(0);

	if (projectSettings_.g.verbose > 1)
	{
		j = PAD + (HEADER_18 * 6) + PAD + COMMENT;
		*rft_file << std::endl << std::endl << "Euler pole rotation parameters" << std::endl <<
			"-------------------------------------------" << std::endl << std::endl;
		*rft_file << std::setw(PAD) << std::left << "Plate" << 
			std::setw(HEADER_18) << std::right << "Pole Latitude" <<
			std::setw(HEADER_18) << std::right << "Pole Longitude" <<
			std::setw(HEADER_18) << std::right << "Euler Rot. Rate" <<
			std::setw(HEADER_18) << std::right << "X Rot. Rate" <<
			std::setw(HEADER_18) << std::right << "Y Rot. Rate" <<
			std::setw(HEADER_18) << std::right << "Z Rot. Rate" <<
			std::setw(PAD) << " " <<
			std::setw(COMMENT) << std::left << "Reference" << std::endl;
		UINT32 i;
		for (i=0; i<j; ++i)
			*rft_file << "-";
		*rft_file << std::endl;
	}
	
	for (it_plate_motion_euler i = plate_motion_eulers_.begin();
		i != plate_motion_eulers_.end(); ++i)
	{
		pmm.plate_name = i->plate_name;
		pmm.pole_param_author = i->pole_param_author;
		
		r_rot = Radians<double>(i->pole_rotation_rate);		// rotation rate
		r_lat = Radians<double>(i->pole_latitude);			// pole latitude
		r_lon = Radians<double>(i->pole_longitude);			// pole longitude
		
		pmm.x_rotation = r_rot * cos(r_lat) * cos(r_lon) * RAD_TO_SEC / 1000.;		// cos(r_lat) * cos(r_lon)
		pmm.y_rotation = r_rot * cos(r_lat) * sin(r_lon) * RAD_TO_SEC / 1000.;		// cos(r_lat) * sin(r_lon)
		pmm.z_rotation = r_rot * sin(r_lat) * RAD_TO_SEC / 1000.;		            // sin(r_lat) * r_rot

		plate_motion_cartesians_.push_back(pmm);

		if (projectSettings_.g.verbose > 1)
			*rft_file << std::setw(PAD) << std::left << i->plate_name << 
				std::setw(HEADER_18) << std::right << std::fixed << std::setprecision(4) << i->pole_latitude <<
				std::setw(HEADER_18) << std::right << std::fixed << std::setprecision(4) << i->pole_longitude <<
				std::setw(HEADER_18) << std::right << std::fixed << std::setprecision(4) << i->pole_rotation_rate <<
				std::setw(HEADER_18) << std::right << std::fixed << std::setprecision(6) << pmm.x_rotation <<
				std::setw(HEADER_18) << std::right << std::fixed << std::setprecision(6) << pmm.y_rotation <<
				std::setw(HEADER_18) << std::right << std::fixed << std::setprecision(6) << pmm.z_rotation <<
				std::setw(PAD) << " " <<
				std::left << i->pole_param_author << std::endl;
	}
}

// Load substitutions for WGS84 and WGS84 (...)
void dna_reftran::LoadWGS84FrameSubstitutions()
{
	frameSubsPtr frameSubstitution;

	_frameSubstitutions.clear();
	
	// WGS84 (transit) and WGS84 to ITRF90
	frameSubstitution.reset(new WGS84_TRANSIT_ITRF90<std::string, UINT32, double>);
	_frameSubstitutions.push_back(frameSubstitution);
	frameSubstitution.reset(new WGS84_ITRF90<std::string, UINT32, double>);
	_frameSubstitutions.push_back(frameSubstitution);

	// WGS84 (G730) and WGS84 to ITRF91
	frameSubstitution.reset(new WGS84_G730_ITRF91<std::string, UINT32, double>);
	_frameSubstitutions.push_back(frameSubstitution);
	frameSubstitution.reset(new WGS84_ITRF91<std::string, UINT32, double>);
	_frameSubstitutions.push_back(frameSubstitution);

	// WGS84 (G873) and WGS84 to ITRF94
	frameSubstitution.reset(new WGS84_G873_ITRF94<std::string, UINT32, double>);
	_frameSubstitutions.push_back(frameSubstitution);
	frameSubstitution.reset(new WGS84_ITRF94<std::string, UINT32, double>);
	_frameSubstitutions.push_back(frameSubstitution);

	// WGS84 (G1150) and WGS84 to ITRF2000
	frameSubstitution.reset(new WGS84_G1150_ITRF2000<std::string, UINT32, double>);
	_frameSubstitutions.push_back(frameSubstitution);
	frameSubstitution.reset(new WGS84_ITRF2000<std::string, UINT32, double>);
	_frameSubstitutions.push_back(frameSubstitution);

	// WGS84 (G1674) and WGS84 to ITRF2008
	frameSubstitution.reset(new WGS84_G1674_ITRF2008<std::string, UINT32, double>);
	_frameSubstitutions.push_back(frameSubstitution);
	frameSubstitution.reset(new WGS84_ITRF2008_1<std::string, UINT32, double>);
	_frameSubstitutions.push_back(frameSubstitution);

	// WGS84 (G1762) and WGS84 to ITRF2008
	frameSubstitution.reset(new WGS84_G1762_ITRF2008<std::string, UINT32, double>);
	_frameSubstitutions.push_back(frameSubstitution);
	frameSubstitution.reset(new WGS84_ITRF2008_2<std::string, UINT32, double>);
	_frameSubstitutions.push_back(frameSubstitution);

	// WGS84 (G2139) and WGS84 to ITRF2014
	frameSubstitution.reset(new WGS84_G2139_ITRF2014<std::string, UINT32, double>);
	_frameSubstitutions.push_back(frameSubstitution);
	frameSubstitution.reset(new WGS84_ITRF2014<std::string, UINT32, double>);
	_frameSubstitutions.push_back(frameSubstitution);

	std::sort(_frameSubstitutions.begin(), _frameSubstitutions.end(), 
		CompareSubstituteOnFrameName< frame_substitutions_t<std::string, UINT32, double>, std::string>());

}

void dna_reftran::LogFrameSubstitutions(std::vector<string_string_pair>& substitutions, const std::string& type)
{
	// Sort, count and remove duplicates
	std::vector<string_string_pair>::iterator _it_sub_1, _it_sub_2, _it_sub_newend;
	UINT32 i(0);
	vUINT32 subs;

	// sort
	std::sort(substitutions.begin(), substitutions.end());

	// count unique pairs
	for (_it_sub_1 = substitutions.begin();
		_it_sub_1 != substitutions.end();)
	{
		subs.push_back(1);
		for (_it_sub_2 = _it_sub_1 + 1; _it_sub_2 != substitutions.end(); ++_it_sub_2)
		{
			if (*_it_sub_1 == *_it_sub_2)
			{
				subs.at(i)++;
				++_it_sub_1;
			}
			else
				break;
		}
		i++;
		if (_it_sub_2 == substitutions.end())
			break;
		_it_sub_1 = _it_sub_2;
	}

	// remove duplicates
	_it_sub_newend = unique(substitutions.begin(), substitutions.end());
	UINT32 subs_count = static_cast<UINT32>(_it_sub_newend - substitutions.begin());
	if (_it_sub_newend != substitutions.end())
		substitutions.resize(_it_sub_newend - substitutions.begin());

	if (subs.size() != substitutions.size())
		// something went wrong
		return;

	// print
	vUINT32::iterator _it_subs = subs.begin();
	std::stringstream ss1, ss2;
	ss1 << type << " reference frame substitutions";
	ss2 << "(" << subs_count << ")";
	
	*rft_file << std::endl << std::endl <<
		std::setw(PRINT_VAL_PAD) << std::left << ss1.str() <<
		std::setw(NUMERIC_WIDTH) << std::right << ss2.str() << std::endl;
	*rft_file << std::string(PRINT_VAL_PAD + NUMERIC_WIDTH, '-') << std::endl;

	for_each(
		substitutions.begin(),
		substitutions.end(),
		[this, &_it_subs](string_string_pair& substitution) {
			*rft_file << std::setw(BLOCK) << std::left << substitution.first <<
				" --> " <<
				std::setw(BLOCK) << substitution.second <<
				std::setw(HEADER_20) << std::right << *(_it_subs++) << std::endl;
		}
	);
	
}

	
void dna_reftran::ApplyStationFrameSubstitutions()
{
	// loop through binary station records and replace occurrences of 
	// the frame to be replaced with a substitute
	it_vstn_t stn_it;
	std::string epsgSubstitute;

	_v_stn_substitutions.clear();
	
	for (stn_it = bstBinaryRecords_.begin(); stn_it != bstBinaryRecords_.end(); ++stn_it)
	{
		try {
			if (IsolateandApplySubstitute(stn_it->epsgCode, stn_it->epoch, epsgSubstitute))
			{
				_v_stn_substitutions.push_back(string_string_pair(
					datumFromEpsgString(std::string(stn_it->epsgCode)),
					datumFromEpsgString(epsgSubstitute)));
				strcpy(stn_it->epsgCode, epsgSubstitute.c_str());
			}
		}
		catch (const RefTranException& e) 
		{
			std::stringstream error_msg;
			error_msg << std::endl <<
				"    - Station:          " << stn_it->stationName << std::endl <<
				"    - Frame and epoch:  " << datumFromEpsgString<std::string>(stn_it->epsgCode) << 
				" (no epoch)" << std::endl;

			switch (e.exception_type())
			{
			case REFTRAN_WGS84_TRANS_UNSUPPORTED:
			{
				std::stringstream throw_msg;
				throw_msg << e.what() << error_msg.str() << std::endl;
				throw RefTranException(throw_msg.str(), REFTRAN_WGS84_TRANS_UNSUPPORTED);
				break;
			}
			default:
				throw RefTranException(e.what());
				break;
			}
		}
	}

	if (_v_stn_substitutions.empty())
		return;
	if (projectSettings_.g.verbose < 2)
		return;

	LogFrameSubstitutions(_v_stn_substitutions, "Station");
}
	
void dna_reftran::ApplyMeasurementFrameSubstitutions()
{
	// loop through binary station records and replace occurrences of 
	// the frame to be replaced with a substitute
	it_vmsr_t msr_it;
	std::string epsgSubstitute;

	_v_msr_substitutions.clear();

	for (msr_it = bmsBinaryRecords_.begin(); msr_it != bmsBinaryRecords_.end(); ++msr_it)
	{
		try {
			if (IsolateandApplySubstitute(msr_it->epsgCode, msr_it->epoch, epsgSubstitute))
			{
				if (msr_it->measStart == xMeas)
					_v_msr_substitutions.push_back(string_string_pair(
						datumFromEpsgString(std::string(msr_it->epsgCode)),
						datumFromEpsgString(epsgSubstitute)));
				strcpy(msr_it->epsgCode, epsgSubstitute.c_str());
			}
		}
		catch (const RefTranException& e)
		{
			std::stringstream error_msg;
			error_msg << std::endl <<
				"    - Measurement type: " << measurement_name<char, std::string>(msr_it->measType) << std::endl <<
				"    - From:             " << bstBinaryRecords_.at(msr_it->station1).stationName << std::endl <<
				"    - To:               " << bstBinaryRecords_.at(msr_it->station2).stationName << std::endl <<
				"    - Frame and epoch:  " << datumFromEpsgString<std::string>(msr_it->epsgCode) << 
				" (no epoch)" << std::endl;

			switch (e.exception_type())
			{
			case REFTRAN_WGS84_TRANS_UNSUPPORTED:
			{
				std::stringstream throw_msg;
				throw_msg << e.what() << error_msg.str() << std::endl;
				throw RefTranException(throw_msg.str(), REFTRAN_WGS84_TRANS_UNSUPPORTED);
				break;
			}
			default:
				throw RefTranException(e.what());
				break;
			}
		}
	}

	if (_v_msr_substitutions.empty())
		return;
	if (projectSettings_.g.verbose < 2)
		return;

	LogFrameSubstitutions(_v_msr_substitutions, "Measurement");
}
	

bool dna_reftran::IsolateandApplySubstitute(const std::string& epsgCode, const std::string& stnEpoch, std::string& epsgSubstitute)
{
	_it_vframesubptr _it_subst = _frameSubstitutions.begin();

	std::string frame;
	frame = datumFromEpsgCode<std::string, UINT32>(LongFromString<UINT32>(epsgCode));

	// first, find the first occurrence of the substitute in _frameSubstitutions 
	if ((_it_subst = binary_search_substitution(
		_it_subst,
		_frameSubstitutions.end(),
		frame)) == _frameSubstitutions.end())
	{
		// Frame not found in substitutions
		return false;
	}

	epsgSubstitute = "";

	if (isEpsgStringWGS84Ensemble(epsgCode))
	{
		if (stnEpoch.empty())
		{
			std::stringstream throw_msg;
			throw_msg << " Cannot perform a reference frame substitution for data on '" << frame << "'" << std::endl <<
				"  without a valid epoch.  '" << frame << "' refers to the \"World Geodetic System 1984" << std::endl <<
				"  (WGS 84) ensemble\".  When transforming stations and measurements from the" << std::endl <<
				"  WGS 84 ensemble, each record must be accompanied with an epoch.  Refer to" << std::endl <<
				"  the DynAdjust User's Guide (\"Configuring import options\") for information" << std::endl <<
				"  on how to achieve reliable transformation results using WGS 84." << std::endl;
			throw RefTranException(throw_msg.str(), REFTRAN_WGS84_TRANS_UNSUPPORTED);
		}

		// In this case, use the epoch to identify the correct substitution
		boost::gregorian::date epoch = dateFromString<boost::gregorian::date>(stnEpoch);

		while (_it_subst != _frameSubstitutions.end())
		{
			if (epoch >= _it_subst->get()->getFromEpoch() &&
				epoch <= _it_subst->get()->getToEpoch())
			{
				epsgSubstitute = _it_subst->get()->getSubstituteName();
				break;
			}
			_it_subst++;
		}
	}
	else
		epsgSubstitute = _it_subst->get()->getSubstituteName();	

	if (epsgSubstitute.empty())
		return false;

	epsgSubstitute = epsgStringFromName<std::string>(epsgSubstitute);

	return true;
}


void dna_reftran::LoadBinaryStationFile(const std::string& bstfileName)
{
	try {
		// Load binary stations data.  Throws runtime_error on failure.
		BstFile bst;
		bst.LoadFile(bstfileName, &bstBinaryRecords_, bst_meta_);
	}
	catch (const std::runtime_error& e) {
		throw RefTranException(e.what());
	}
}


void dna_reftran::WriteBinaryStationFile(const std::string& bstfileName)
{
	std::string strEpsg(datumTo_.GetEpsgCode_s());
	std::string strEpoch(datumTo_.GetEpoch_s());

	// update binary file meta
	snprintf(bst_meta_.modifiedBy, sizeof(bst_meta_.modifiedBy), "%s", __BINARY_NAME__);
	snprintf(bst_meta_.epsgCode, sizeof(bst_meta_.epsgCode), "%s", strEpsg.substr(0, STN_EPSG_WIDTH).c_str());
	snprintf(bst_meta_.epoch, sizeof(bst_meta_.epoch), "%s", strEpoch.substr(0, STN_EPOCH_WIDTH).c_str());
	bst_meta_.reftran = true;

	try {
		// write binary stations data.  Throws runtime_error on failure.
		BstFile bst;
		bst.WriteFile(bstfileName, &bstBinaryRecords_, bst_meta_);
	}
	catch (const std::runtime_error& e) {
		throw RefTranException(e.what());
	}
}

void dna_reftran::LoadBinaryMeasurementFile(const std::string& bmsfileName)
{
	try {
		// Load binary measurements data.  Throws runtime_error on failure.
		BmsFile bms;
		bms.LoadFile(bmsfileName, &bmsBinaryRecords_, bms_meta_);
	}
	catch (const std::runtime_error& e) {
		throw RefTranException(e.what());
	}
}
	

void dna_reftran::WriteBinaryMeasurementFile(const std::string& bmsfileName)
{
	std::string strEpsg(datumTo_.GetEpsgCode_s());
	std::string strEpoch(datumTo_.GetEpoch_s());

	// update binary file meta
	snprintf(bms_meta_.modifiedBy, sizeof(bms_meta_.modifiedBy), "%s", __BINARY_NAME__);
	snprintf(bms_meta_.epsgCode, sizeof(bms_meta_.epsgCode), "%s", strEpsg.substr(0, STN_EPSG_WIDTH).c_str());
	snprintf(bms_meta_.epoch, sizeof(bms_meta_.epoch), "%s", strEpoch.substr(0, STN_EPOCH_WIDTH).c_str());
	bms_meta_.reftran = true;

	try {
		// write binary measurement data.  Throws runtime_error on failure.
		BmsFile bms;
		bms.WriteFile(bmsfileName, &bmsBinaryRecords_, bms_meta_);
	}
	catch (const std::runtime_error& e) {
		throw RefTranException(e.what());
	}
}
	

UINT32 dna_reftran::DetermineTectonicPlate(const std::string& plate)
{
	it_pair_string_vUINT32 it_plate = equal_range(vplateMap_.begin(), vplateMap_.end(), 
		plate, StationNameIDCompareName());

	if (it_plate.first == it_plate.second)
	{
		std::stringstream error_msg;
		error_msg << "An attempt to find plate motion model parameters failed for " << 
			plate << ".";
		throw RefTranException(error_msg.str());
	}

	return it_plate.first->second;
}
	

// Use plate motion model to project coordinates between epochs
// This method selects plate motion parameters based on the 
// tectonic plate within which the point lies.
void dna_reftran::ObtainPlateMotionParameters(it_vstn_t& stn_it, double* reduced_parameters,
	const CDnaDatum& datumFrom, const CDnaDatum& datumTo, transformation_parameter_set& transformParameters, double& timeElapsed)
{
	transformationType transformation_type = __plate_motion_model__;

	if (projectSettings_.r.plate_model_option == 0)
		getAustralianPlateMotionModelParameters<double, UINT32>(transformParameters);
	else
	{
		UINT32 plateIndex = DetermineTectonicPlate(stn_it->plate);
		setDefinedPlateMotionModelParameters<double, UINT32>(transformParameters,
			plate_motion_cartesians_.at(plateIndex).x_rotation,
			plate_motion_cartesians_.at(plateIndex).y_rotation,
			plate_motion_cartesians_.at(plateIndex).z_rotation);
	}

	timeElapsed = DetermineElapsedTime(datumFrom, datumTo, 
		transformParameters, transformation_type);
	ReduceParameters<double>(transformParameters.parameters_, 
		reduced_parameters, timeElapsed);
}
	
// This function is called only when:
//	- there are no direct parameters between the 'from' and 'to' frames.
//  - at least one of the input/output frames is static
// NOTE: Joining uses itrf2014 as the "stepping" frame
// If an exception is thrown, there's not much more we can do
void dna_reftran::JoinTransformationParameters(it_vstn_t& stn_it, double* reduced_parameters, 
	const CDnaDatum& datumFrom, const CDnaDatum& datumTo, transformation_parameter_set& transformParameters,
	transformationType transType, const matrix_2d& coordinates)
{
	
	transformation_parameter_set transP_a, transP_b;
	CDnaDatum datumStep(epsgCodeFromName<UINT32, std::string>(ITRF2014_s));
	
	// Set the reference epoch to the epoch of the data
	switch (transType)
	{
	case __static_to_static__:
		// this case doesn't need an epoch
		break;
	case __static_to_dynamic__:
	case __dynamic_to_dynamic__:
		datumStep.SetEpoch(datumTo.GetEpoch());
		break;
	case __dynamic_to_static__:
		datumStep.SetEpoch(datumFrom.GetEpoch());
		break;
	default:
		break;
	}

	double reduced_parameters_step[7];
	double timeElapsed_a(0.0), timeElapsed_b(0.0);
	std::string epoch_step;

	transformationType transformation_type;


	// datumFrom -> Step
	try
	{
		// reset transformation_type
		if (datumFrom.isStatic())
			transformation_type = __static_to_step__;
		else // (datumFrom.isDynamic())
			transformation_type = __dynamic_to_step__;
		
		epoch_step = datumStep.GetEpoch_s();

		// 1. Obtain parameters for the first step (to ITRF2014)
		// The time elapsed will be a function of:
		//		- epoch of datumFrom, and 
		//		- reference epoch of the transformation parameters
		ObtainHelmertParameters(datumFrom, datumStep, transP_a,
			timeElapsed_a, transformation_type);

		// 2. Reduce the first-step parameters to the appropriate unit and format
		if (datumFrom.isDynamic() || datumStep.isDynamic())
			ReduceParameters<double>(transP_a.parameters_, reduced_parameters, timeElapsed_a);
		else
			ReduceParameters<double>(transP_a.parameters_, reduced_parameters, timeElapsed_a, false);

#ifdef _MSDEBUG
		TRACE("Reduced parameters (1):\n");
		TRACE("%11.8f\n", reduced_parameters[0]);
		TRACE("%11.8f\n", reduced_parameters[1]);
		TRACE("%11.8f\n", reduced_parameters[2]);
		TRACE("%11.8g\n", reduced_parameters[3]);
		TRACE("%11.8g\n", reduced_parameters[4]);
		TRACE("%11.8g\n", reduced_parameters[5]);
		TRACE("%11.8g\n\n", reduced_parameters[6]);
#endif

	}
	catch (RefTranException& rft) {
		// No parameters exist
		switch (rft.exception_type())
		{
		case REFTRAN_TRANS_ON_PLATE_REQUIRED:
			ObtainPlateMotionParameters(stn_it, reduced_parameters, datumFrom, datumTo, transformParameters, timeElapsed_a);
			epoch_step = datumTo.GetEpoch_s();
			break;
		default:
			std::stringstream error_msg;
			error_msg << "Attempting to join parameters between " << 
				datumFrom_.GetName() << " and " << datumStep.GetName() << ":" << std::endl <<
				"    " << rft.what();
			throw RefTranException(error_msg.str());
		}
	}

	if (projectSettings_.g.verbose > 1 && data_type_ == stn_data)
	{
		matrix_2d coordinates_step(3, 1);
		Transform_7parameter<double>(coordinates, coordinates_step, reduced_parameters);

		*rft_file << std::setw(PAD3) << std::left << "JN" << 
			std::setw(STATION) << std::left << stn_it->stationName << 
			std::setw(REL) << std::right << ITRF2014_s <<
			std::setw(REL) << std::right << epoch_step <<
			std::setw(PAD3) << " " <<
			std::setw(PAD) << std::left << stn_it->plate << 
			std::setw(MEASR) << std::right << std::fixed << std::setprecision(4) << coordinates_step.get(0, 0) <<
			std::setw(MEASR) << std::right << std::fixed << std::setprecision(4) << coordinates_step.get(1, 0) <<
			std::setw(MEASR) << std::right << std::fixed << std::setprecision(4) << coordinates_step.get(2, 0);

		if (projectSettings_.g.verbose > 2)
		{
			*rft_file << std::setw(PACORR) << std::right << std::fixed << std::setprecision(4) << reduced_parameters[0] <<
				std::setw(PACORR) << std::right << std::fixed << std::setprecision(4) << reduced_parameters[1] <<
				std::setw(PACORR) << std::right << std::fixed << std::setprecision(4) << reduced_parameters[2] <<
				std::setw(PACORR) << std::right << std::scientific << std::setprecision(4) << reduced_parameters[3] <<
				std::setw(PACORR) << std::right << std::scientific << std::setprecision(4) << reduced_parameters[4] <<
				std::setw(PACORR) << std::right << std::scientific << std::setprecision(4) << reduced_parameters[5] <<
				std::setw(PACORR) << std::right << std::scientific << std::setprecision(4) << reduced_parameters[6] <<
				std::setw(PACORR) << std::right << std::fixed << std::setprecision(4) << timeElapsed_a;
		}
		*rft_file << std::endl;
	}

	// Step -> datumTo
	try
	{
		// reset transformation_type
		if (datumTo.isStatic())
			transformation_type = __step_to_static__;
		else // (datumFrom.isDynamic())
			transformation_type = __step_to_dynamic__;

		// 3. Obtain parameters for the second step (to "DatumTo_")
		//    timeElapsed_b is set, but not used.
		ObtainHelmertParameters(datumStep, datumTo, transP_b,
			timeElapsed_b, transformation_type);
		
		// 4. Reduce the second-step parameters to the appropriate unit and format
		if (datumStep.isDynamic() || datumTo.isDynamic())
			ReduceParameters<double>(transP_b.parameters_, reduced_parameters_step, timeElapsed_b);
		else
			ReduceParameters<double>(transP_b.parameters_, reduced_parameters_step, timeElapsed_b, false);

#ifdef _MSDEBUG
		TRACE("Reduced parameters (2):\n");
		TRACE("%11.8f\n", reduced_parameters_step[0]);
		TRACE("%11.8f\n", reduced_parameters_step[1]);
		TRACE("%11.8f\n", reduced_parameters_step[2]);
		TRACE("%11.8g\n", reduced_parameters_step[3]);
		TRACE("%11.8g\n", reduced_parameters_step[4]);
		TRACE("%11.8g\n", reduced_parameters_step[5]);
		TRACE("%11.8g\n\n", reduced_parameters_step[6]);
#endif

	}
	catch (RefTranException& rft) {
		// No parameters exist
		switch (rft.exception_type())
		{
		case REFTRAN_TRANS_ON_PLATE_REQUIRED:
			ObtainPlateMotionParameters(stn_it, reduced_parameters_step, datumFrom, datumTo, transformParameters, timeElapsed_b);
			break;
		default:
			std::stringstream error_msg;
			error_msg << "Attempting to join parameters between " <<
				datumStep.GetName() << " and " << datumTo.GetName() << ":" << std::endl <<
				"    " << rft.what();
			throw RefTranException(error_msg.str());
		}
	}

	// 6. Sum the reduced parameters
	reduced_parameters[0] += reduced_parameters_step[0];
	reduced_parameters[1] += reduced_parameters_step[1];
	reduced_parameters[2] += reduced_parameters_step[2];
	reduced_parameters[3] += reduced_parameters_step[3];
	reduced_parameters[4] += reduced_parameters_step[4];
	reduced_parameters[5] += reduced_parameters_step[5];
	reduced_parameters[6] += reduced_parameters_step[6];

#ifdef _MSDEBUG
	TRACE("Joined (1 + 2) parameters:\n");
	TRACE("%11.8f\n", reduced_parameters[0]);
	TRACE("%11.8f\n", reduced_parameters[1]);
	TRACE("%11.8f\n", reduced_parameters[2]);
	TRACE("%11.8g\n", reduced_parameters[3]);
	TRACE("%11.8g\n", reduced_parameters[4]);
	TRACE("%11.8g\n", reduced_parameters[5]);
	TRACE("%11.8g\n\n", reduced_parameters[6]);
#endif

	// Since this was the result of a 'join', force 
	// reinitialisation of datumFrom and datumTo codes
	transformParameters.from_to_ = uint32_uint32_pair(0, 0);
}
	

void dna_reftran::TransformEpochs_PlateMotionModel(it_vstn_t& stn_it, const matrix_2d& coordinates, matrix_2d& coordinates_mod,
	const CDnaDatum& datumFrom, const CDnaDatum& datumTo)
{
	// Propagate parameters using the appropriate Plate Motion Model.
	// If this attempt raises an exception, it will be caught by the 
	// calling method
	double reduced_parameters[7];
	double timeElapsed;
	transformation_parameter_set transformParameters;
	ObtainPlateMotionParameters(stn_it, reduced_parameters, datumFrom, datumTo, transformParameters, timeElapsed);

#ifdef _MSDEBUG
	TRACE("Final reduced parameters:\n");
	TRACE("%11.8f\n", reduced_parameters[0]);
	TRACE("%11.8f\n", reduced_parameters[1]);
	TRACE("%11.8f\n", reduced_parameters[2]);
	TRACE("%11.8g\n", reduced_parameters[3]);
	TRACE("%11.8g\n", reduced_parameters[4]);
	TRACE("%11.8g\n", reduced_parameters[5]);
	TRACE("%11.8g\n\n", reduced_parameters[6]);
#endif

	// Transform!
	Transform_7parameter<double>(coordinates, coordinates_mod, reduced_parameters);

	if (projectSettings_.g.verbose > 1 && data_type_ == stn_data)
	{
		*rft_file << std::setw(PAD3) << std::left << "PM" << 
			std::setw(STATION) << std::left << stn_it->stationName << 
			std::setw(REL) << std::right << datumTo.GetName() <<
			std::setw(REL) << std::right << datumTo.GetEpoch_s() <<
			std::setw(PAD3) << " " <<
			std::setw(PAD) << std::left << stn_it->plate << 
			std::setw(MEASR) << std::right << std::fixed << std::setprecision(4) << coordinates_mod.get(0, 0) <<
			std::setw(MEASR) << std::right << std::fixed << std::setprecision(4) << coordinates_mod.get(1, 0) <<
			std::setw(MEASR) << std::right << std::fixed << std::setprecision(4) << coordinates_mod.get(2, 0);

		if (projectSettings_.g.verbose > 2)
		{
			*rft_file << std::setw(PACORR) << std::right << std::fixed << std::setprecision(4) << reduced_parameters[0] <<
				std::setw(PACORR) << std::right << std::fixed << std::setprecision(4) << reduced_parameters[1] <<
				std::setw(PACORR) << std::right << std::fixed << std::setprecision(4) << reduced_parameters[2] <<
				std::setw(PACORR) << std::right << std::scientific << std::setprecision(4) << reduced_parameters[3] <<
				std::setw(PACORR) << std::right << std::scientific << std::setprecision(4) << reduced_parameters[4] <<
				std::setw(PACORR) << std::right << std::scientific << std::setprecision(4) << reduced_parameters[5] <<
				std::setw(PACORR) << std::right << std::scientific << std::setprecision(4) << reduced_parameters[6] <<
				std::setw(PACORR) << std::right << std::fixed << std::setprecision(4) << timeElapsed;
		}
		*rft_file << std::endl;
	}

#ifdef _MSDEBUG
	std::stringstream ss;
	ss << "coords, " << datumFrom.GetName() << " @ " << referenceEpoch<double>(datumFrom.GetEpoch());
	coordinates.trace(ss.str(), "%.4f ");
	ss.str("");
	ss << "coords_mod, " << datumTo.GetName() << " @ " << transformParameters.reference_epoch_;
	coordinates_mod.trace(ss.str(), "%.4f ");
#endif
}
	

void dna_reftran::TransformFrames_PlateMotionModel(it_vstn_t& stn_it, const matrix_2d& coordinates, matrix_2d& coordinates_mod,
	const CDnaDatum& datumFrom, const CDnaDatum& datumTo, transformation_parameter_set& transformParameters)
{
	// For this scenario, three steps are involved.  The sequence is very similar to
	// TransformFrames_Join, the exception being the PMM must be used to transform
	// between epochs once the coordinates have been transformed to the step frame (ITRF).
	// The steps in this scenario are:
	//	1. Transform datumFrom to ITRF2014 (using the epoch of the input dynamic frame)
	//	2. Apply PMM (transform from input epoch to output epoch on ITRF2014)
	//  3. Transform ITRF2014 to datumTo (using the epoch of the output dynamic frame)

	// Create the step datum and set the epoch to the epoch of the input data
	CDnaDatum datumStep1(epsgCodeFromName<UINT32, std::string>(ITRF2014_s), datumFrom.GetEpoch());
	CDnaDatum datumStep2(epsgCodeFromName<UINT32, std::string>(ITRF2014_s), datumTo.GetEpoch());

	matrix_2d coordinates_tmp(coordinates);

	//	1. Transform datumFrom to ITRF2014 (using the epoch of the input dynamic frame)
	if (datumFrom.GetEpsgCode_i() != datumStep1.GetEpsgCode_i())
	{
#ifdef _MSDEBUG
		TRACE("Step 1: Helmert transformation to ITRF2014 @ reference epoch\n");
#endif

		TransformFrames_WithoutPlateMotionModel(stn_it, coordinates, coordinates_mod, datumFrom, datumStep1,
			transformParameters, __dynamic_to_dynamic__);
		coordinates_tmp = coordinates_mod;
	}

#ifdef _MSDEBUG
	TRACE("Step 2: Plate motion model transformation on ITRF from epoch of input frame to epoch of output frame\n");
	std::stringstream ss;
	ss << "Transforming from " << datumStep1.GetName() << " @ " << datumStep1.GetEpoch() << std::endl;
	ss << "               to " << datumStep2.GetName() << " @ " << datumStep2.GetEpoch();
	TRACE("%s\n", ss.str().c_str());
#endif


	//	2. Apply PMM (transform from input epoch to output epoch on ITRF2014)
	// Create the step datum and set the epoch to the epoch of the output data
		
	TransformEpochs_PlateMotionModel(stn_it, coordinates_tmp, coordinates_mod, datumStep1, datumStep2);

	//  3. Transform ITRF2014 to datumTo (using the epoch of the output dynamic frame)
	if (datumStep2.GetEpsgCode_i() != datumTo.GetEpsgCode_i())
	{
#ifdef _MSDEBUG
		TRACE("Step 3: Helmert transformation from ITRF to output frame @ epoch\n");
#endif

		coordinates_tmp = coordinates_mod;
		TransformFrames_WithoutPlateMotionModel(stn_it, coordinates_tmp, coordinates_mod, datumStep2, datumTo,
			transformParameters, __dynamic_to_dynamic__);
	}
}
	
// Called by:
// - Transform(), when transforming between static-static, static-dynamic and dynamic-static
// - TransformDynamic(), when transforming between dynamic-dynamic when input/output epochs are the same
// Primary scenario - parameters exist between the two frames.
// Secondary scenario - no parameters exist, in which case a REFTRAN_DIRECT_PARAMS_UNAVAILABLE exception is
// thrown and caught, and TransformFrames_Join handles the transformation using ITRF2014 as a step.

void dna_reftran::TransformFrames_WithoutPlateMotionModel(it_vstn_t& stn_it, const matrix_2d& coordinates, matrix_2d& coordinates_mod,
	const CDnaDatum& datumFrom, const CDnaDatum& datumTo, transformation_parameter_set& transformParameters,
	transformationType transType)
{
	double timeElapsed;
	double reduced_parameters[7];

	try {
		// Primary scenario:   parameters exist between datumFrom and datumTo_
		// Secondary scenario: no parameters exist, in which case an exception
		//   of type REFTRAN_DIRECT_PARAMS_UNAVAILABLE will be thrown. In this case,
		//   the exception will be caught and the transformation handled in two steps
		ObtainHelmertParameters(datumFrom, datumTo, transformParameters,
			timeElapsed, transType);

#ifdef _MSDEBUG
		std::stringstream ss;
		ss << "Transforming from " << datumFrom.GetName() << " @ " << std::fixed << std::setprecision(4) << referenceEpoch<double>(datumFrom.GetEpoch()) << std::endl;
		ss << "               to " << datumTo.GetName() << " @ " << std::fixed << std::setprecision(4) << transformParameters.reference_epoch_;
		TRACE("%s\n", ss.str().c_str());

		//TRACE("Raw parameters:\n");
		//TRACE("%11.8f\n", transformParameters.parameters_[0]);
		//TRACE("%11.8f\n", transformParameters.parameters_[1]);
		//TRACE("%11.8f\n", transformParameters.parameters_[2]);
		//TRACE("%11.8f\n", transformParameters.parameters_[3]);
		//TRACE("%11.8f\n", transformParameters.parameters_[4]);
		//TRACE("%11.8f\n", transformParameters.parameters_[5]);
		//TRACE("%11.8f\n", transformParameters.parameters_[6]);
		//TRACE("%11.8f\n", transformParameters.parameters_[7]);
		//TRACE("%11.8f\n", transformParameters.parameters_[8]);
		//TRACE("%11.8f\n", transformParameters.parameters_[9]);
		//TRACE("%11.8f\n", transformParameters.parameters_[10]);
		//TRACE("%11.8f\n", transformParameters.parameters_[11]);
		//TRACE("%11.8f\n", transformParameters.parameters_[12]);
		//TRACE("%11.8f\n\n", transformParameters.parameters_[13]);
#endif

		// Reduce the parameters to the appropriate unit and format
		if (datumFrom.isDynamic() || datumTo.isDynamic())
			ReduceParameters<double>(transformParameters.parameters_, reduced_parameters, timeElapsed);
		else
			ReduceParameters<double>(transformParameters.parameters_, reduced_parameters, timeElapsed, false);

#ifdef _MSDEBUG
		//TRACE("Final reduced parameters:\n");
		//TRACE("%11.8f\n", reduced_parameters[0]);
		//TRACE("%11.8f\n", reduced_parameters[1]);
		//TRACE("%11.8f\n", reduced_parameters[2]);
		//TRACE("%11.8g\n", reduced_parameters[3]);
		//TRACE("%11.8g\n", reduced_parameters[4]);
		//TRACE("%11.8g\n", reduced_parameters[5]);
		//TRACE("%11.8g\n\n", reduced_parameters[6]);
#endif

		// Transform!
		Transform_7parameter<double>(coordinates, coordinates_mod, reduced_parameters);

		if (projectSettings_.g.verbose > 1 && data_type_ == stn_data)
		{
			*rft_file << std::setw(PAD3) << std::left << TransformationType<std::string, transformationType>(transType) << 
				std::setw(STATION) << std::left << stn_it->stationName << 
				std::setw(REL) << std::right << datumTo.GetName() <<
				std::setw(REL) << std::right << datumTo.GetEpoch_s() <<
				std::setw(PAD3) << " " <<
				std::setw(PAD) << std::left << stn_it->plate << 
				std::setw(MEASR) << std::right << std::fixed << std::setprecision(4) << coordinates_mod.get(0, 0) <<
				std::setw(MEASR) << std::right << std::fixed << std::setprecision(4) << coordinates_mod.get(1, 0) <<
				std::setw(MEASR) << std::right << std::fixed << std::setprecision(4) << coordinates_mod.get(2, 0);

			if (projectSettings_.g.verbose > 2)
			{
				*rft_file << std::setw(PACORR) << std::right << std::fixed << std::setprecision(4) << reduced_parameters[0] <<
					std::setw(PACORR) << std::right << std::fixed << std::setprecision(4) << reduced_parameters[1] <<
					std::setw(PACORR) << std::right << std::fixed << std::setprecision(4) << reduced_parameters[2] <<
					std::setw(PACORR) << std::right << std::scientific << std::setprecision(4) << reduced_parameters[3] <<
					std::setw(PACORR) << std::right << std::scientific << std::setprecision(4) << reduced_parameters[4] <<
					std::setw(PACORR) << std::right << std::scientific << std::setprecision(4) << reduced_parameters[5] <<
					std::setw(PACORR) << std::right << std::scientific << std::setprecision(4) << reduced_parameters[6] <<
					std::setw(PACORR) << std::right << std::fixed << std::setprecision(4) << timeElapsed;
			}
			*rft_file << std::endl;
		}
	

#ifdef _MSDEBUG
		ss.str("");
		ss << "coords, " << datumFrom.GetName() << " @ " << referenceEpoch<double>(datumFrom.GetEpoch());
		coordinates.trace(ss.str(), "%.4f ");
		ss.str("");
		ss << "coords_mod, " << datumTo.GetName() << " @ " << transformParameters.reference_epoch_;
		coordinates_mod.trace(ss.str(), "%.4f ");
#endif

	}
	catch (RefTranException& rft) {
		// No parameters exist
		switch (rft.exception_type())
		{
		case REFTRAN_DIRECT_PARAMS_UNAVAILABLE:
			TransformFrames_Join(stn_it, coordinates, coordinates_mod, datumFrom, datumTo, transformParameters, transType);
			break;
		case REFTRAN_WGS84_TRANS_UNSUPPORTED:
			throw RefTranException(rft.what(), REFTRAN_WGS84_TRANS_UNSUPPORTED);
			break;
		default:
			throw RefTranException(rft.what());
		}
	}
}


void dna_reftran::TransformDynamic(it_vstn_t& stn_it, const matrix_2d& coordinates, matrix_2d& coordinates_mod,
	const CDnaDatum& datumFrom, const CDnaDatum& datumTo, transformation_parameter_set& transformParameters,
	transformationType transType)
{
	epochSimilarity epoch_similarity;
	frameSimilarity frame_similarity;

	//   If REFTRAN_TRANS_ON_PLATE_REQUIRED is thrown, the exception is caught 
	//   and the plate motion model is applied.  The only problem is, the 
	//   plate motion model currently implemented within DynAdjust is the 
	//   Australian plate motion model. Hence, this cannot be used for
	//   points not on the Australian plate!!!

	// Are input and output epochs different?
	if (datumFrom.GetEpoch() == datumTo.GetEpoch())
		epoch_similarity = __epoch_epoch_same__;
	else
		epoch_similarity = __epoch_epoch_diff__;

	// Are input and output frames different?
	if (datumFrom.GetEpsgCode_i() == datumTo_.GetEpsgCode_i())
		frame_similarity = __frame_frame_same__;
	else
		frame_similarity = __frame_frame_diff__;

#ifdef _MSDEBUG
	std::stringstream ss;
	ss << "Transforming from " << datumFrom.GetName() << " @ " << datumFrom.GetEpoch() << std::endl;
	ss << "               to " << datumTo.GetName() << " @ " << datumTo.GetEpoch();
	TRACE("%s\n", ss.str().c_str());
#endif

	if (frame_similarity == __frame_frame_diff__ &&
		epoch_similarity == __epoch_epoch_same__)
	{
		TransformFrames_WithoutPlateMotionModel(stn_it, coordinates, coordinates_mod, datumFrom, datumTo,
			transformParameters, transType);
	}
	else
	{
		// if the input and output reference frames are equal, OR
		// if the input and output reference frames are different AND
		// the input and output epochs are different

		TransformFrames_PlateMotionModel(stn_it, coordinates, coordinates_mod, datumFrom, datumTo,
			transformParameters);
	}
}
	

void dna_reftran::Transform(it_vstn_t& stn_it, const matrix_2d& coordinates, matrix_2d& coordinates_mod,
	const CDnaDatum& datumFrom, transformation_parameter_set& transformParameters)
{
	transformationType transformation_type;	

	datumFrom_ = datumFrom;

	if (datumFrom.isStatic() && datumTo_.isStatic())
		transformation_type = __static_to_static__;
	else if (datumFrom.isStatic() && datumTo_.isDynamic())
		transformation_type = __static_to_dynamic__;
	else if (datumFrom.isDynamic() && datumTo_.isStatic())
		transformation_type = __dynamic_to_static__;
	else //if (datumFrom.isDynamic() && datumTo_.isDynamic())
		transformation_type = __dynamic_to_dynamic__;

	// Handle transformation based upon type
	switch (transformation_type)
	{
	// 1. Static to static
	//   a) if params exist, direct trans
	//   b) if no params exist, throw/catch, join on ITRF2014
	case __static_to_static__:
	
	// 2. Static to dynamic
	//   a) if params exist, direct trans
	//   b) if no params exist, throw/catch, join on ITRF2014
	case __static_to_dynamic__:

	// 3. Dynamic to static
	//   a) if params exist, direct trans
	//   b) if no params exist, throw/catch, join on ITRF2014
	case __dynamic_to_static__:

		// At least one frame is static
		TransformFrames_WithoutPlateMotionModel(stn_it, coordinates, coordinates_mod, datumFrom, datumTo_,
			transformParameters, transformation_type);
		break;

	// 4. Dynamic to dynamic
	//   a) if different frames, and epoch_in = epoch_out
	//     i) if params exist, direct trans
	//    ii) if no params exist, throw/catch join on ITRF2014
	//
	//   b) if different frames, and epoch_in != epoch_out 
	//     ** 3 steps involving ITRF2014 and PMM are req'd regardless 
	//        of whether params exist or not
	//     1. trans from-datum to ITRF2014
	//     2. apply PMM between epochs
	//     3. trans ITRF2014 to to-datum
	//   c) if same frame, and epoch_in != epoch_out 
	//     1. apply PMM between epochs
	case __dynamic_to_dynamic__:
		// Both frames are dynamic
		TransformDynamic(stn_it, coordinates, coordinates_mod, datumFrom, datumTo_,
			transformParameters, transformation_type);
		break;
	default:
		break;
	}
}
	
// At this point, a REFTRAN_DIRECT_PARAMS_UNAVAILABLE exception has been 
// thrown, caught and re-directed to here.
// Try joining two sets associated with ITRF2014
// If this attempt raises an exception, it will be caught by the 
// calling method
void dna_reftran::TransformFrames_Join(it_vstn_t& stn_it, const matrix_2d& coordinates, matrix_2d& coordinates_mod,
	const CDnaDatum& datumFrom, const CDnaDatum& datumTo, transformation_parameter_set& transformParameters,
	transformationType transType)
{
	double reduced_parameters[7];
	JoinTransformationParameters(stn_it, reduced_parameters, datumFrom, datumTo, transformParameters, transType, coordinates);

	// Transform!
	Transform_7parameter<double>(coordinates, coordinates_mod, reduced_parameters);

#ifdef _MSDEBUG
	coordinates.trace("coords", "%.4f ");
	coordinates_mod.trace("coords_mod", "%.4f ");
#endif
}
	

double dna_reftran::DetermineElapsedTime(const CDnaDatum& datumFrom, const CDnaDatum& datumTo, 
	transformation_parameter_set& transParams, transformationType transType)
{
	std::stringstream ss;
	double dTime(0.0);

	try
	{
		// Formula is always (Dt = t = t0) where the transformation is in the direction of
		// the published parameters, and where:
		//	  t  = epoch of the 'data' coordinate
		//    t0 = epoch of the reference epoch of the transformation parameters
		//  
		// If the reverse direction is required, then:
		//    t  = epoch of the 'to' coordinate+
		//    parameters are negated
		//
		// For example:
		// 1. From ITRF2008@2016.202 -> To GDA94
		//    - Parameters are Dawson Woods ITRF->GDA94 (ref. epoch t0 = 1994.0)
		//    - parameterSet.paramDirection_ is set to __paramForward__
		//    - elapsedTime is calculated as:
		//        Dt =  t(from) - t0
		//        Dt = 2016.202 - 1994.0
		//
		// 2. From GDA94 -> To ITRF2008@2016.202
		//    - Parameters are Dawson Woods ITRF->GDA94 (ref. epoch t0 = 1994.0)
		//    - parameterSet.paramDirection_ is set to __paramReverse__
		//    - values in parameterSet.parameters_ are multiplied by -1 by reverse() method
		//    - elapsedTime is calculated as:
		//        Dt =    t(to) - t0
		//        Dt = 2016.202 - 1994.0
		//
		// TOD - this is not right.  Since from and to epochs are different,
		//       a step is needed
		// 3. From ITRF2008@2016.202 -> To ITRF2005@2012.112
		//    - Parameters are IERS ITRF2008->ITRF2005 (ref. epoch t0 = 2000.0)
		//    - parameterSet.paramDirection_ is set to __paramForward__
		//    - elapsedTime is calculated as:
		//        Dt =    t(to) - t(from)
		//        Dt = 2016.202 - 2012.112
		//
		// TOD - this is not right.  Since from and to epochs are different,
		//       a step is needed
		// 4. From ITRF2005@2012.112 -> To ITRF2008@2016.202
		//    - Parameters are IERS ITRF2008->ITRF2005 (ref. epoch t0 = 2000.0)
		//    - parameterSet.paramDirection_ is set to __paramReverse__
		//    - values in parameterSet.parameters_ are multiplied by -1 by reverse() method
		//    - elapsedTime is calculated as:
		//        Dt =    t(to) - t(from)
		//        Dt = 2016.202 - 2012.112
		//

#ifdef _MSDEBUG
		ss << "From frame: " << datumFrom.GetName() << "  -> to frame: " << datumTo.GetName();
		TRACE("%s\n", ss.str().c_str());
#endif
		ss.str("");

		boost::gregorian::date dt;
		double dt0(transParams.reference_epoch_);
		
		switch (transType)
		{
		case __static_to_static__:
			// There is no influence from time on the parameters
			dTime = 0.;

#ifdef _MSDEBUG
			ss << "Static epochs - 0.0 (zero) time difference";
			TRACE("%s\n", ss.str().c_str());
#endif
			return dTime;
			break;

		case __static_to_dynamic__:
		case __dynamic_to_static__:
		
			if (transParams.paramDirection_ == __paramForward__)
				dt = datumFrom.GetEpoch();
			else 
				dt = datumTo.GetEpoch();
			break;

		case __dynamic_to_dynamic__:

			dt = datumFrom.GetEpoch();			
			break;

		case __plate_motion_model__:
			
			// Assume PMM parameters are always given in the positive direction
			dt = datumFrom.GetEpoch();
			dt0 = referenceEpoch<double>(datumTo.GetEpoch());
			break;


			// When a step frame is involved, either the input datum or output datum
			// will be dynamic.  In both cases, the epoch of the coordinates on the 
			// dynamic datum is needed. 

		case __static_to_step__:
			// datumFrom  = static frame
			// datumTo    = step frame (ITRF2014)
			// datumTo_   = dynamic frame*
		case __step_to_dynamic__:
			// datumFrom  = step frame (ITRF2014)
			// datumTo    = dynamic frame
			// datumTo_   = dynamic frame*
			dt = datumTo_.GetEpoch();
			break;
		
		case __dynamic_to_step__:
			// datumFrom_ = dynamic frame*
			// datumFrom  = dynamic frame
			// datumTo    = step frame (ITRF2014)
		case __step_to_static__:
			// datumFrom_ = dynamic frame*
			// datumFrom  = step frame (ITRF2014)
			// datumTo    = static frame
			dt = datumFrom_.GetEpoch();
		}

		dTime = elapsedTime<double>(dt, dt0);

#ifdef _MSDEBUG
		ss << "From epoch: " << std::fixed << std::setprecision(4) << referenceEpoch<double>(dt) << " -> to epoch: " << std::fixed << std::setprecision(4) << dt0 <<
			" = " << std::setprecision(4) << std::fixed << dTime;
		TRACE("%s\n", ss.str().c_str());
#endif

	}
	catch (std::out_of_range& e)
	{
		throw RefTranException(e.what());
	}
	catch (...)
	{
		std::stringstream ss;
		ss << "DetermineElapsedTime(): an error occurred whilst computing the elapsed time.";
		throw RefTranException(ss.str());
	}

	return dTime;
}
	

void dna_reftran::ObtainHelmertParameters(const CDnaDatum& datumFrom, const CDnaDatum& datumTo,
	transformation_parameter_set& transformParameters, double& timeElapsed, transformationType transType)
{
	// Does this transformation require a different set of parameters?
	if (transformParameters.from_to_.first != datumFrom.GetEpsgCode_i() ||
		transformParameters.from_to_.second != datumTo.GetEpsgCode_i())
	{		
		transformParameters.from_to_.first = datumFrom.GetEpsgCode_i();
		transformParameters.from_to_.second = datumTo.GetEpsgCode_i();

		// Primary scenario:   parameters exist between datumFrom and datumTo_
		// Secondary scenario: no parameters exist, in which case an exception
		//   of type REFTRAN_DIRECT_PARAMS_UNAVAILABLE will be thrown.
		determineHelmertParameters<UINT32>(transformParameters);
	}

	// TODO - here, timeElapsed is computed for every record.  In the interest of achieving
	// better performance, this function call could be modified such that DetermineElapsedTime
	// is called only when the epochs from the previous measurement are different.
	if (datumFrom.isDynamic() || datumTo.isDynamic())
		timeElapsed = DetermineElapsedTime(datumFrom, datumTo, transformParameters, transType);
	else
		timeElapsed = 0.;
}

void dna_reftran::TransformStationRecords(const std::string& newFrame, const std::string& newEpoch)
{
	it_vstn_t stn_it;
	CDnaDatum datumFrom;
	data_type_ = stn_data;

	UINT32 j(0);

	if (projectSettings_.g.verbose > 1)
	{
		j = (PAD3 * 2) + PAD + STATION + (REL * 2) + (MEASR * 3);
		*rft_file << std::endl << std::endl << "Station coordinate transformations" << std::endl <<
			"-------------------------------------------" << std::endl << std::endl;
		*rft_file << std::setw(PAD3) << std::left << "ID" << 
			std::setw(STATION) << std::left << "Station" << 
			std::setw(REL) << std::right << "Frame" <<
			std::setw(REL) << std::right << "Epoch" <<
			std::setw(PAD3) << " " <<
			std::setw(PAD) << std::left << "Plate" << 
			std::setw(MEASR) << std::right << "X" <<
			std::setw(MEASR) << std::right << "Y" <<
			std::setw(MEASR) << std::right << "Z";
		
		// Print reduced transformation parameters
		if (projectSettings_.g.verbose > 2)
		{
			*rft_file << std::setw(PACORR) << std::right << "dX" <<
				std::setw(PACORR) << std::right << "dY" <<
				std::setw(PACORR) << std::right << "dZ" <<
				std::setw(PACORR) << std::right << "Sc" <<
				std::setw(PACORR) << std::right << "rX" <<
				std::setw(PACORR) << std::right << "rY" <<
				std::setw(PACORR) << std::right << "rZ" <<
				std::setw(PACORR) << std::right << "dt";
			j += (8 * PACORR);
		}

		*rft_file << std::endl;

		UINT32 i;
		for (i=0; i<j; ++i)
			*rft_file << "-";
		*rft_file << std::endl;
	}

	// Create the transformation parameters to be used for the
	// entire set of station records.  If a station is in a different
	// frame, obtain new parameters
	transformation_parameter_set transformationParameters;
	
	transformationPerformed_ = false;
	m_stnsTransformed = m_stnsNotTransformed = 0;

#ifdef _MSDEBUG
	TRACE("\nTransforming stations...\n\n");
#endif

	try {
		// 1. Get the datum (and epoch) of the desired system
		datumTo_.SetDatumFromName(newFrame, newEpoch);
		
		// 2. For every station, get the datum, then transform
		//    TransformStation takes
		for (stn_it=bstBinaryRecords_.begin(); stn_it!=bstBinaryRecords_.end(); ++stn_it)
		{
			// a. Get datum of current station
			if (trimstr(std::string(stn_it->epoch)).empty())
				datumFrom.SetDatum(stn_it->epsgCode);
			else
				datumFrom.SetDatumFromEpsg(stn_it->epsgCode, stn_it->epoch);

			// b. test if a transformation is required
			if (datumFrom == datumTo_)
			{
				m_stnsNotTransformed++;
				continue;
			}

			// c. Transform
			TransformStation(stn_it, datumFrom, transformationParameters);

			// d. Update meta
			snprintf(stn_it->epsgCode, sizeof(stn_it->epsgCode), "%s", datumTo_.GetEpsgCode_s().c_str());
			snprintf(stn_it->epoch, sizeof(stn_it->epoch), "%s", datumTo_.GetEpoch_s().c_str());
			transformationPerformed_ = true;
			m_stnsTransformed++;
		}
	}
	catch (const std::runtime_error& e)
	{
		std::stringstream error_msg;
		error_msg << e.what() << std::endl <<
			"    - Station:          " << stn_it->stationName << std::endl <<
			"    - Frame and epoch:  " << datumFromEpsgString<std::string>(stn_it->epsgCode) << " @ " <<
			stn_it->epoch << std::endl;
		throw RefTranException(e.what());
	}
	catch (const RefTranException& e)
	{
		std::stringstream error_msg;
		error_msg << e.what() << std::endl <<
			"    - Station:          " << stn_it->stationName << std::endl <<
			"    - Frame and epoch:  " << datumFromEpsgString<std::string>(stn_it->epsgCode) << " @ " <<
			stn_it->epoch << std::endl;		
		throw RefTranException(e.what());
				
	}
}
	

void dna_reftran::TransformStation(it_vstn_t& stn_it, const CDnaDatum& datumFrom, 
		transformation_parameter_set& transformParameters)
{
	matrix_2d coordinates(3, 1), coordinates_mod(3, 1);

	// 1. Convert to cartesian.  Why?  Native form of coordinates in bst file 
	// is geographic
	GeoToCart<double>(
		stn_it->currentLatitude, stn_it->currentLongitude, stn_it->currentHeight, 
		coordinates.getelementref(0, 0), 
		coordinates.getelementref(1, 0),
		coordinates.getelementref(2, 0), 
		datumFrom.GetEllipsoidRef());

	if (projectSettings_.g.verbose > 1)
		*rft_file << std::setw(PAD3) << std::left << "FR" << 
			std::setw(STATION) << std::left << stn_it->stationName << 
			std::setw(REL) << std::right << datumFrom.GetName() <<
			std::setw(REL) << std::right << datumFrom.GetEpoch_s() <<
			std::setw(PAD3) << " " <<
			std::setw(PAD) << std::left << stn_it->plate << 
			std::setw(MEASR) << std::right << std::fixed << std::setprecision(4) << coordinates.get(0, 0) <<
			std::setw(MEASR) << std::right << std::fixed << std::setprecision(4) << coordinates.get(1, 0) <<
			std::setw(MEASR) << std::right << std::fixed << std::setprecision(4) << coordinates.get(2, 0) << std::endl;

	// 2. Transform!
	Transform(stn_it, coordinates, coordinates_mod, datumFrom, transformParameters);

	if (projectSettings_.g.verbose > 1)
		*rft_file << std::setw(PAD3) << std::left << "TO" << 
			std::setw(STATION) << std::left << stn_it->stationName << 
			std::setw(REL) << std::right << datumTo_.GetName() <<
			std::setw(REL) << std::right << datumTo_.GetEpoch_s() <<
			std::setw(PAD3) << " " <<
			std::setw(PAD) << std::left << stn_it->plate << 
			std::setw(MEASR) << std::right << std::fixed << std::setprecision(4) << coordinates_mod.get(0, 0) <<
			std::setw(MEASR) << std::right << std::fixed << std::setprecision(4) << coordinates_mod.get(1, 0) <<
			std::setw(MEASR) << std::right << std::fixed << std::setprecision(4) << coordinates_mod.get(2, 0) << std::endl;

	// 3. Convert back to geographic
	CartToGeo<double>(coordinates_mod.get(0, 0), 
		coordinates_mod.get(1, 0), coordinates_mod.get(2, 0),
		&stn_it->currentLatitude, &stn_it->currentLongitude, &stn_it->currentHeight, 
		datumTo_.GetEllipsoidRef());
}
	

void dna_reftran::TransformMeasurementRecords(const std::string& newFrame, const std::string& newEpoch)
{
	it_vmsr_t msr_it;
	CDnaDatum datumFrom;
	data_type_ = msr_data;
	
	// Create the transformation parameters to be used for the
	// entire set of measurement records.  If a measurement is in a different
	// frame, obtain new parameters
	transformation_parameter_set transformationParameters;
	
	transformationPerformed_ = false;
	m_msrsTransformed = m_msrsNotTransformed = 0;

#ifdef _MSDEBUG
	TRACE("\nTransforming measurements...\n\n");
#endif

	try {
		// 1. Get the datum (and epoch) of the desired system
		datumTo_.SetDatumFromName(newFrame, newEpoch);
	
		// 2. For every measurement, get the datum, determine parameters, then transform
		for (msr_it=bmsBinaryRecords_.begin(); msr_it!=bmsBinaryRecords_.end(); ++msr_it)
		{
			// a. ignore measurements not subject to geodetic datum
			switch (msr_it->measType)
			{
			// Local reference frame measurements not subject to
			// ellipsoid or reference frame/epoch
			case 'A':	// Horizontal angle
			case 'B':	// Geodetic azimuth
			case 'C':	// Chord dist
			case 'D':	// Directions
			case 'E':	// Ellipsoid arc
			case 'H':	// Orthometric height
			case 'I':	// Astronomic latitude
			case 'J':	// Astronomic longitude
			case 'K':	// Astronomic azimuth
			case 'L':	// Level difference
			case 'M':	// MSL arc
			case 'P':	// Geodetic latitude
			case 'Q':	// Geodetic longitude
			case 'R':	// Ellipsoidal height
			case 'S':	// Slope distance
			case 'V':	// Zenith distance
			case 'Z':	// Vertical angle
				continue;
			}

			if (msr_it->measStart != xMeas)
				continue;

			// b. Get datum of current measurement
			if (trimstr(std::string(msr_it->epoch)).empty())
				datumFrom.SetDatum(msr_it->epsgCode);
			else
				datumFrom.SetDatumFromEpsg(msr_it->epsgCode, msr_it->epoch);
			
			// c. test if a transformation is required
			if (datumFrom == datumTo_)
			{
				m_msrsNotTransformed++;
				continue;
			}

			// d. Transform!
			TransformMeasurement(msr_it, datumFrom, transformationParameters);
			
			// e. Update meta
			transformationPerformed_ = true;
		}
	}
	catch (const std::runtime_error& e)
	{
		std::stringstream error_msg;
		error_msg << e.what() << std::endl <<
			"    - Measurement type: " << measurement_name<char, std::string>(msr_it->measType) << std::endl <<
			"    - From:             " << bstBinaryRecords_.at(msr_it->station1).stationName << std::endl <<
			"    - To:               " << bstBinaryRecords_.at(msr_it->station2).stationName << std::endl <<
			"    - Frame and epoch:  " << datumFromEpsgString<std::string>(msr_it->epsgCode) << " @ " << 
			msr_it->epoch << std::endl;
		throw RefTranException(error_msg.str());
	}
	catch (const RefTranException& e) 
	{
		std::stringstream error_msg;
		error_msg << e.what() << std::endl <<
			"    - Measurement type: " << measurement_name<char, std::string>(msr_it->measType) << std::endl <<
			"    - From:             " << bstBinaryRecords_.at(msr_it->station1).stationName << std::endl <<
			"    - To:               " << bstBinaryRecords_.at(msr_it->station2).stationName << std::endl <<
			"    - Frame and epoch:  " << datumFromEpsgString<std::string>(msr_it->epsgCode) << " @ " <<
			msr_it->epoch << std::endl;

		throw RefTranException(error_msg.str());
	}
}

void dna_reftran::TransformMeasurement(it_vmsr_t& msr_it, const CDnaDatum& datumFrom, 
		transformation_parameter_set& transformParameters)
{
	// a. ignore measurements not subject to geodetic datum
	switch (msr_it->measType)
	{
	case 'B':	// Geodetic azimuth
	case 'C':	// Chord dist
	case 'E':	// Ellipsoid arc
	case 'I':	// Astronomic latitude
	case 'J':	// Astronomic longitude
	case 'K':	// Astronomic azimuth
	case 'M':	// MSL arc
	case 'P':	// Geodetic latitude
	case 'Q':	// Geodetic longitude
	case 'R':	// Ellipsoidal height
		// Some form of transformation if reference frame/epoch is different?
		break;
	case 'G':	// GPS baseline
	case 'X':	// GPS baseline cluster
		TransformMeasurement_GX(msr_it, datumFrom, transformParameters);
		break;
	case 'Y':	// GPS point cluster
		// 7 or 14 parameter transformation required
		TransformMeasurement_Y(msr_it, datumFrom, transformParameters);
		break;
	}
}


void dna_reftran::TransformMeasurement_GX(it_vmsr_t& msr_it, const CDnaDatum& datumFrom, 
		transformation_parameter_set& transformParameters)
{
	UINT32 cluster_bsl, bsl_count(msr_it->vectorCount1);
	UINT32 covariance_count;
	it_vstn_t stn1_it, stn2_it;

	matrix_2d coordinates1(3, 1), coordinates2(3, 1), 
		coordinates1_mod(3, 1), coordinates2_mod(3, 1);

	//// 1. Get upper triangular a-priori measurements variance matrix
	//matrix_2d vmat;
	//GetGPSVarianceMatrix<it_vmsr_t>(msr_it, vmat);
	
	for (cluster_bsl=0; cluster_bsl<bsl_count; ++cluster_bsl)
	{
		covariance_count = msr_it->vectorCount2;
		
		// Get stations
		stn1_it = bstBinaryRecords_.begin() + msr_it->station1;
		stn2_it = bstBinaryRecords_.begin() + msr_it->station2;

		// 1. Convert station 1 coordinates to cartesian.
		GeoToCart<double>(
			stn1_it->currentLatitude, stn1_it->currentLongitude, stn1_it->currentHeight, 
			coordinates1.getelementref(0, 0), 
			coordinates1.getelementref(1, 0),
			coordinates1.getelementref(2, 0), 
			datumFrom.GetEllipsoidRef());

		// 2. Compute station 2 coordinates from station 1 and GPS vector
		coordinates2 = coordinates1;
		
		// X
		coordinates2.elementadd(0, 0, msr_it->term1);
		msr_it++;

		// Y
		coordinates2.elementadd(1, 0, msr_it->term1);
		msr_it++;

		// Z
		coordinates2.elementadd(2, 0, msr_it->term1);
		
		// 3. Transform station 1
		Transform(stn1_it, coordinates1, coordinates1_mod, datumFrom, transformParameters);
				
		// 5. Transform station 2
		Transform(stn2_it, coordinates2, coordinates2_mod, datumFrom, transformParameters);
		
		// update transformation count
		m_msrsTransformed += 3;

		// Calculate and assign new 'transformed' baseline elements
		// Go back to X element
		msr_it -= 2;
		msr_it->term1 = coordinates2_mod.get(0, 0) - coordinates1_mod.get(0, 0);
		//TRACE("\nTransformed baseline\n");
		//TRACE("%.4f\n", msr_it->term1);
		snprintf(msr_it->epsgCode, sizeof(msr_it->epsgCode), "%s", datumTo_.GetEpsgCode_s().c_str());
		snprintf(msr_it->epoch, sizeof(msr_it->epoch), "%s", datumTo_.GetEpoch_s().c_str());
		msr_it++;
		
		// Y
		msr_it->term1 = coordinates2_mod.get(1, 0) - coordinates1_mod.get(1, 0);
		//TRACE("%.4f\n", msr_it->term1);
		snprintf(msr_it->epsgCode, sizeof(msr_it->epsgCode), "%s", datumTo_.GetEpsgCode_s().c_str());
		snprintf(msr_it->epoch, sizeof(msr_it->epoch), "%s", datumTo_.GetEpoch_s().c_str());
		msr_it++;
		// Z
		msr_it->term1 = coordinates2_mod.get(2, 0) - coordinates1_mod.get(2, 0);
		//TRACE("%.4f\n", msr_it->term1);
		snprintf(msr_it->epsgCode, sizeof(msr_it->epsgCode), "%s", datumTo_.GetEpsgCode_s().c_str());
		snprintf(msr_it->epoch, sizeof(msr_it->epoch), "%s", datumTo_.GetEpoch_s().c_str());

		// skip covariances until next baseline
		if (covariance_count > 0)
			msr_it += covariance_count * 3;

		if (covariance_count > 0)
			msr_it++;
	}
			
}

void dna_reftran::TransformMeasurement_Y(it_vmsr_t& msr_it, const CDnaDatum& datumFrom, 
		transformation_parameter_set& transformParameters)
{
	UINT32 cluster_pnt, pnt_count(msr_it->vectorCount1);
	UINT32 covariance_count;
	it_vstn_t stn_it;

	matrix_2d coordinates(3, 1), coordinates_mod(3, 1);
	
	_COORD_TYPE_ coordType(CDnaStation::GetCoordTypeC(msr_it->coordType));

	for (cluster_pnt=0; cluster_pnt<pnt_count; ++cluster_pnt)
	{
		covariance_count = msr_it->vectorCount2;
		
		// Get X
		coordinates.put(0, 0, msr_it->term1);
		msr_it++;

		// Get Y
		coordinates.put(1, 0, msr_it->term1);
		msr_it++;

		// Get Z
		coordinates.put(2, 0, msr_it->term1);

		// Go back to X element
		msr_it -= 2;
		
		// LLH?  Convert to Cartesian?
		if (coordType == LLH_type_i || coordType == LLh_type_i)
		{
			GeoToCart<double>(
				coordinates.get(0, 0), 
				coordinates.get(1, 0),
				coordinates.get(2, 0), 
				coordinates.getelementref(0, 0), 
				coordinates.getelementref(1, 0),
				coordinates.getelementref(2, 0), 
				datumFrom.GetEllipsoidRef());
		}
		
		stn_it = bstBinaryRecords_.begin() + msr_it->station1;

		// Transform
		Transform(stn_it, coordinates, coordinates_mod, datumFrom, transformParameters);

		// update transformation count
		m_msrsTransformed += 3;

		// LLH?  Convert to geographic?
		if (coordType == LLH_type_i || coordType == LLh_type_i)
		{
			CartToGeo<double>(
				coordinates_mod.get(0, 0), 
				coordinates_mod.get(1, 0),
				coordinates_mod.get(2, 0), 
				coordinates_mod.getelementref(0, 0), 
				coordinates_mod.getelementref(1, 0),
				coordinates_mod.getelementref(2, 0), 
				datumFrom.GetEllipsoidRef());
		}
		
		// Assign 'transformed' elements
		msr_it->term1 = coordinates_mod.get(0, 0);
		snprintf(msr_it->epsgCode, sizeof(msr_it->epsgCode), "%s", datumTo_.GetEpsgCode_s().c_str());
		snprintf(msr_it->epoch, sizeof(msr_it->epoch), "%s", datumTo_.GetEpoch_s().c_str());
		msr_it++;
		// Y
		msr_it->term1 = coordinates_mod.get(1, 0);
		snprintf(msr_it->epsgCode, sizeof(msr_it->epsgCode), "%s", datumTo_.GetEpsgCode_s().c_str());
		snprintf(msr_it->epoch, sizeof(msr_it->epoch), "%s", datumTo_.GetEpoch_s().c_str());
		msr_it++;
		// Z
		msr_it->term1 = coordinates_mod.get(2, 0);
		snprintf(msr_it->epsgCode, sizeof(msr_it->epsgCode), "%s", datumTo_.GetEpsgCode_s().c_str());
		snprintf(msr_it->epoch, sizeof(msr_it->epoch), "%s", datumTo_.GetEpoch_s().c_str());

		// skip covariances until next point		
		if (covariance_count > 0)
			msr_it += covariance_count * 3;

		if (covariance_count > 0)
			msr_it++;
	}
			
}
	

// First item in the file is a UINT32 value - the number of records in the file
// All records are of type UINT32
void dna_reftran::LoadDatabaseId()
{
	if (databaseIDsLoaded_)
		return;

	std::string dbid_filename = formPath<std::string>(projectSettings_.g.output_folder,
		projectSettings_.g.network_name, "dbid");

	std::stringstream ss;
	v_msr_db_map_.clear();

	std::ifstream dbid_file;
	try {
		// Create geoid file.  Throws runtime_error on failure.
		file_opener(dbid_file, dbid_filename,
			std::ios::in | std::ios::binary, binary);
	}
	catch (const std::runtime_error& e) {
		ss << e.what();
		throw std::runtime_error(ss.str());
	}

	UINT32 r, recordCount;
	msr_database_id_map rec;

	try {
		// get size and reserve vector size
		dbid_file.read(reinterpret_cast<char*>(&recordCount), sizeof(UINT32));
		v_msr_db_map_.reserve(recordCount);

		UINT16 val;

		for (r = 0; r < recordCount; r++)
		{
			// Read data
			dbid_file.read(reinterpret_cast<char*>(&rec.msr_id), sizeof(UINT32));
			dbid_file.read(reinterpret_cast<char*>(&rec.cluster_id), sizeof(UINT32));
			dbid_file.read(reinterpret_cast<char*>(&val), sizeof(UINT16));
			rec.is_msr_id_set = val_uint<bool, UINT16>(val);
			dbid_file.read(reinterpret_cast<char*>(&val), sizeof(UINT16));
			rec.is_cls_id_set = val_uint<bool, UINT16>(val);

			// push back
			v_msr_db_map_.push_back(rec);
		}

		dbid_file.close();
		databaseIDsLoaded_ = true;
		if (v_msr_db_map_.size() > 0)
			databaseIDsSet_ = true;
		else
			databaseIDsSet_ = false;
	}
	catch (const std::ifstream::failure& f) {
		ss << f.what();
		throw std::runtime_error(ss.str());
	}
}
	

void dna_reftran::SerialiseDNA(const std::string& stnfilename, const std::string& msrfilename, bool flagUnused)
{
	try {
		// Load Database IDs
		LoadDatabaseId();
	}
	catch (const std::runtime_error& e) {
		throw RefTranException(e.what());
	}

	CDnaProjection projection;
	try {
		// write binary stations data.  Throws runtime_error on failure.
		dna_io_dna dna;
		std::string comment(" transformed to ");
		comment.append(datumTo_.GetName());
		if (datumTo_.isDynamic())
			comment.append(", epoch ").append(datumTo_.GetEpoch_s());
		comment.append(".  Exported by reftran.");
		dna.set_dbid_ptr(&v_msr_db_map_);
		dna.write_dna_files(&bstBinaryRecords_, &bmsBinaryRecords_, 
			stnfilename, msrfilename, projectSettings_.g.network_name,
			datumTo_, projection, flagUnused, "Station coordinates" + comment, "GNSS measurements" + comment);
	}
	catch (const std::runtime_error& e) {
		throw RefTranException(e.what());
	}
}

void dna_reftran::SerialiseDynaML(const std::string& stnfilename, const std::string& msrfilename, bool flagUnused)
{
	try {
		// Load Database IDs
		LoadDatabaseId();
	}
	catch (const std::runtime_error& e) {
		throw RefTranException(e.what());
	}


	// Open DynaML Station file
	std::ofstream dynaml_stn_file;
	try {
		// Create DynaML station file.  Throws runtime_error on failure.
		file_opener(dynaml_stn_file, stnfilename);
		// write header
		dynaml_header(dynaml_stn_file, "Station File", datumTo_.GetName(), datumTo_.GetEpoch_s());
		dynaml_comment(dynaml_stn_file, "File type:    Station file");
		dynaml_comment(dynaml_stn_file, "Project name: " + projectSettings_.g.network_name);
		// Write header comment line
		std::string comment("Station coordinates transformed to ");
		comment.append(datumTo_.GetName());
		if (datumTo_.isDynamic())
			comment.append(", epoch ").append(datumTo_.GetEpoch_s());
		comment.append(".  Exported by reftran.");
		dynaml_comment(dynaml_stn_file, comment);
	}
	catch (const std::runtime_error& e) {
		throw RefTranException(e.what());
	}

	CDnaProjection projection;
	
	// Write DynaML Station file.  Throws runtime_error on failure.
	try {
		SerialiseDynaMLStn(&dynaml_stn_file, projection, flagUnused);
		dynaml_footer(dynaml_stn_file);
		dynaml_stn_file.close();
	}
	catch (const std::runtime_error& e) {
		throw RefTranException(e.what());
	}

	// Open DynaML Measurement file
	std::ofstream dynaml_msr_file;
	try {
		// Create DynaML measurement file.  Throws runtime_error on failure.
		file_opener(dynaml_msr_file, msrfilename);
		// write header
		dynaml_header(dynaml_msr_file, "Measurement File", datumTo_.GetName(), datumTo_.GetEpoch_s());
		dynaml_comment(dynaml_msr_file, "File type:    Measurement file");
		dynaml_comment(dynaml_msr_file, "Project name: " + projectSettings_.g.network_name);
		// Write header comment line
		std::string comment("GNSS measurements transformed to ");
		comment.append(datumTo_.GetName());
		if (datumTo_.isDynamic())
			comment.append(", epoch ").append(datumTo_.GetEpoch_s());
		comment.append(".  Exported by reftran.");
		dynaml_comment(dynaml_msr_file, comment);
	}
	catch (const std::runtime_error& e) {
		throw RefTranException(e.what());
	}
	
	// Write DynaML Measurement file.  Throws runtime_error on failure.
	try {
		SerialiseDynaMLMsr(&dynaml_msr_file);
		dynaml_footer(dynaml_msr_file);
		dynaml_msr_file.close();
	}
	catch (const std::runtime_error& e) {
		throw RefTranException(e.what());
	}	
}
	

void dna_reftran::SerialiseDynaML(const std::string& xmlfilename, bool flagUnused)
{
	try {
		// Load Database IDs
		LoadDatabaseId();
	}
	catch (const std::runtime_error& e) {
		throw RefTranException(e.what());
	}

	std::ofstream dynaml_xml_file;
	try {
		// Create DynaML station and measurement file.  Throws runtime_error on failure.
		file_opener(dynaml_xml_file, xmlfilename);
		// write header
		dynaml_header(dynaml_xml_file, "Combined File", datumTo_.GetName(), datumTo_.GetEpoch_s());
		// Write header comment line
		std::string comment("Station coordinates and measurements transformed to ");
		comment.append(datumTo_.GetName());
		if (datumTo_.isDynamic())
			comment.append(", epoch ").append(datumTo_.GetEpoch_s());
		comment.append(".  Exported by reftran.");
		dynaml_comment(dynaml_xml_file, comment);
	}
	catch (const std::runtime_error& e) {
		throw RefTranException(e.what());
	}

	CDnaProjection projection;

	try {
		// Write combined DynaML station and binary measurement data.  Throws runtime_error on failure.

		//1. Stations
		SerialiseDynaMLStn(&dynaml_xml_file, projection, flagUnused);

		// 2. Measurements
		SerialiseDynaMLMsr(&dynaml_xml_file);
	}
	catch (const std::runtime_error& e) {
		throw RefTranException(e.what());
	}

	dynaml_footer(dynaml_xml_file);
	dynaml_xml_file.close();


}

void dna_reftran::SerialiseDynaMLStn(std::ofstream* xml_file, CDnaProjection& projection, bool flagUnused/*=false*/)
{
	UINT32 epsgCode(LongFromString<UINT32>(bstBinaryRecords_.at(0).epsgCode));
	std::string datum(datumFromEpsgCode<std::string, UINT32>(epsgCode));
	std::string epoch(referenceepochFromEpsgCode<UINT32>(epsgCode));
	dnaStnPtr stnPtr(new CDnaStation(datum, epoch));
	
	it_vstn_t _it_stn;

	if (flagUnused) 
	{
		for (_it_stn=bstBinaryRecords_.begin(); _it_stn!=bstBinaryRecords_.end(); ++_it_stn)
		{
			if (_it_stn->unusedStation)
				continue;
			stnPtr->SetStationRec(*_it_stn);
			stnPtr->WriteDNAXMLStnCurrentEstimates(xml_file,
				datumTo_.GetEllipsoidRef(), &projection, dynaml);
		}
	}
	else
	{
		for (_it_stn=bstBinaryRecords_.begin(); _it_stn!=bstBinaryRecords_.end(); ++_it_stn)
		{
			stnPtr->SetStationRec(*_it_stn);
			stnPtr->WriteDNAXMLStnCurrentEstimates(xml_file,
				datumTo_.GetEllipsoidRef(), &projection, dynaml);
		}	
	}
}
	

void dna_reftran::SerialiseDynaMLMsr(std::ofstream* xml_file)
{
	dnaMsrPtr msrPtr;
	it_vmsr_t _it_msr;
	std::string comment("");
	size_t dbindex;
	it_vdbid_t _it_dbid;
	
	for (_it_msr=bmsBinaryRecords_.begin(); _it_msr!=bmsBinaryRecords_.end(); ++_it_msr)
	{
		ResetMeasurementPtr<char>(&msrPtr, _it_msr->measType);

		if (databaseIDsSet_)
		{
			dbindex = std::distance(bmsBinaryRecords_.begin(), _it_msr);
			_it_dbid = v_msr_db_map_.begin() + dbindex;
		}

		msrPtr->SetMeasurementRec(bstBinaryRecords_, _it_msr, _it_dbid);
		msrPtr->WriteDynaMLMsr(xml_file, comment);
	}
}

//bool dna_reftran::PrintTransformedStationCoordinatestoSNX()
//{
//	std::ofstream sinex_file;
//
//	try {
//		// Open output file stream.  Throws runtime_error on failure.
//		file_opener(sinex_file, projectSettings_.o._snx_file);
//	}
//	catch (const std::runtime_error& e) {
//		throw RefTranException(e.what());
//	}
//
//	//dna_io_snx snx;
//
//	sinex_file.close();
//
//	return true;
//}


}	// namespace referenceframe
}	// namespace dynadjust
