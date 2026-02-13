//============================================================================
// Name         : dnainterop.cpp
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
// Description  : DynAdjust Interoperability library
//============================================================================

#include <dynadjust/dnaimport/dnainterop.hpp>

//#include <include/io/DynaML-schema.hxx>

using namespace dynadjust::epsg;

MsrTally	g_map_tally;
MsrTally	g_parsemsr_tally;
StnTally	g_parsestn_tally;
UINT32		g_fileOrder;

//boost::random::mt19937 rng;
//boost::random::uniform_real_distribution<double> stdev(0.0, 3.0);
//boost::random::uniform_real_distribution<double> pertu(0.0, 1.0);

//int compare_chararray(const void *a, const void *b)
//{
//	const char* a1 = *(const char**)a;
//	const char* b1 = *(const char**)b;
//	return strcmp(a1,b1);
//}

namespace dynadjust {
namespace dynamlinterop {

std::mutex import_file_mutex;

dna_import::dna_import()
	: percentComplete_(-99.)
	, isProcessing_(false)
	, databaseIDsLoaded_(false)
	, databaseIDsSet_(false)
{
	ifsInputFILE_ = 0;
	g_fileOrder = 0;
	p_parsemsr_tally = &g_parsemsr_tally;
	p_parsestn_tally = &g_parsestn_tally;

#ifdef _MSC_VER
#if (_MSC_VER < 1900)
	{
		// this function is obsolete in MS VC++ 14.0, VS2015
		// Set scientific format to print two places for the exponent
		_set_output_format(_TWO_DIGIT_EXPONENT);
	}
#endif
#endif
	
	m_strProjectDefaultEpsg = DEFAULT_EPSG_S;
	m_strProjectDefaultEpoch = DEFAULT_EPOCH;

	stn_discontinuities_.clear();
}

dna_import::~dna_import()
{

}
	

double dna_import::GetProgress()
{
	// calculate on the fly
	if (!isProcessing_)
		return -1;

	// Obtain exclusive use of the input file pointer
	import_file_mutex.lock();

	try
	{
		if (ifsInputFILE_)
			percentComplete_ = fabs(ifsInputFILE_->tellg() * 100. / sifsFileSize_);
	}
	// Catch any type of error; do nothing.
	catch (...)
	{
		//if (ifsInputFILE_->eof())
		//	return percentComplete_;
		//if (ifsInputFILE_->rdstate() & std::ifstream::eofbit)
		//	return percentComplete_;
	}

	import_file_mutex.unlock();

	return percentComplete_;
}

void dna_import::DetermineBoundingBox()
{
	if (projectSettings_.i.bounding_box.empty())
		return;
	
	if (GetFields(const_cast<char*>(projectSettings_.i.bounding_box.c_str()), ',', false, "ffff", 
		&bbox_upperLat_, &bbox_lowerLon_, &bbox_lowerLat_, &bbox_upperLon_) < 4)
		return;
	
	bbox_upperLat_ = DmstoRad(bbox_upperLat_);
	bbox_upperLon_ = DmstoRad(bbox_upperLon_);
	bbox_lowerLat_ = DmstoRad(bbox_lowerLat_);
	bbox_lowerLon_ = DmstoRad(bbox_lowerLon_);
}
	
void dna_import::BuildExtractStationsList(const std::string& stnList, pvstring vstnList)
{
	// Extract stations from comma delimited string
	try {
		SplitDelimitedString<std::string>(stnList, std::string(","), vstnList);
	}
	catch (...) {
		std::stringstream ss;
		ss << "BuildExtractStationsList(): An error was encountered when parsing " << stnList << "." << std::endl;
		SignalExceptionParse(ss.str(), 0);
	}

	// Have discontinuities been applied?
	if (!projectSettings_.i.apply_discontinuities)
		return;

	if (stn_renamed_.empty())
		return;

	it_string_pair _it_discont_ren;

	_it_vstr _it_stn;
	std::string station_name;
	vstring renamed_stations;

	// loop through all station names in vstnList
	for (_it_stn = vstnList->begin();
		_it_stn != vstnList->end();
		_it_stn++)
	{
		station_name = *_it_stn;
		
		_it_discont_ren = stn_renamed_.begin();

		// Advance through _it_discont_ren for all occurrences of station_name
		while ((_it_discont_ren = lower_bound(_it_discont_ren, stn_renamed_.end(), 
			station_name, ComparePairFirst<std::string>())) != stn_renamed_.end())
		{
			// found a station that has been renamed.			
			// add the discontinuity name to the list
			renamed_stations.push_back(_it_discont_ren->second);
			_it_discont_ren++;
		}
	}

	// Add discontinuity sites and remove duplicates
	if (!renamed_stations.empty())
	{
		vstnList->insert(vstnList->end(), renamed_stations.begin(), renamed_stations.end());
		strip_duplicates(vstnList);
	}
}
	

void dna_import::InitialiseDatum(const std::string& reference_frame, const std::string epoch)
{
	try {
		// Take the default reference frame, set either by the user or
		// the datum_ constructor (GDA2020).  Epoch is that of the default
		// reference frame (indicated by "")
		datum_.SetDatumFromName(reference_frame, epoch);
	}
	catch (const std::runtime_error& e) {
		std::stringstream ss;
		ss << "InitialiseDatum(): An error occurred while initialising " << std::endl << "  the default reference frame.  Details:" << std::endl <<
			"  " << e.what() << std::endl;
		SignalExceptionParse(ss.str(), 0);
	}

	// Get epsg code and epoch from the 'parsed' datum and epoch
	m_strProjectDefaultEpsg = datum_.GetEpsgCode_s();
	m_strProjectDefaultEpoch = datum_.GetEpoch_s();

	if (datum_.GetEpoch() == timeImmemorial<boost::gregorian::date>())
		m_strProjectDefaultEpoch = "";

	// Update binary file meta
	// Note: the following rule applies each time a set of files is loaded via import:
	//	* This method (InitialiseDatum) is called (from dnaimportwrapper) before any files are loaded.
	//    By default, the bst & bms meta are initialised with the reference frame and reference epoch.
	//  * The datum and epoch within the first file (if present) is used to set the default project 
	//    datum. If a datum isn't provided in the first input file, e.g. SINEX file, the default datum
	//    (GDA2020) is used. As each subsequent file is loaded, the default frame and epoch are assumed.
	//  * After all files have been loaded, InitialiseDatum is called again to set the metadata.
	//  * import does not attempt to reconcile multiple default datums and epochs found across the input
	//    files. It will however produce a warning when something different to the default is discovered.
    snprintf(bst_meta_.epsgCode, sizeof(bst_meta_.epsgCode), "%s", m_strProjectDefaultEpsg.substr(0, STN_EPSG_WIDTH).c_str());
    snprintf(bms_meta_.epsgCode, sizeof(bms_meta_.epsgCode), "%s", m_strProjectDefaultEpsg.substr(0, STN_EPSG_WIDTH).c_str());
    snprintf(bst_meta_.epoch, sizeof(bst_meta_.epoch), "%s", m_strProjectDefaultEpoch.substr(0, STN_EPOCH_WIDTH).c_str());
    snprintf(bms_meta_.epoch, sizeof(bms_meta_.epoch), "%s", m_strProjectDefaultEpoch.substr(0, STN_EPOCH_WIDTH).c_str());
}
	

_PARSE_STATUS_ dna_import::ParseInputFile(const std::string& fileName, vdnaStnPtr* vStations, PUINT32 stnCount, 
							   vdnaMsrPtr* vMeasurements, PUINT32 msrCount, 
							   PUINT32 clusterID, input_file_meta_t* input_file_meta, bool firstFile,
							   std::string* success_msg, project_settings* p)
{
	projectSettings_ = *p;
	*success_msg = "";

	DetermineBoundingBox();

	percentComplete_ = -99.0;
	isProcessing_ = true;
	_filespecifiedreferenceframe = false;
	_filespecifiedepoch = false;
	std::stringstream ss;

	try 
	{
		// Obtain exclusive use of the input file pointer
		import_file_mutex.lock();

		if (ifsInputFILE_)
		{
			ifsInputFILE_->close();
			delete ifsInputFILE_;
		}

		ifsInputFILE_ = new std::ifstream;

		// Open and seek to end immediately after opening.
		file_opener(ifsInputFILE_, fileName, std::ios::in | std::ios::ate, ascii, true);
		
		// get file size and return to start
		sifsFileSize_ = (size_t)ifsInputFILE_->tellg();
		ifsInputFILE_->seekg(0, std::ios::beg);

		// release file pointer mutex
		import_file_mutex.unlock();
	}
	catch (const std::ios_base::failure& f) {	
		ss.str("");
		ss << "ParseInputFile(): An error was encountered when opening " << fileName << "." << std::endl << "  " << f.what() << std::endl << "  Check that the file exists and that the file is not already opened.";
		SignalExceptionParse(ss.str(), 0);
	}
	catch (...) {	
		ss.str("");
		ss << "ParseInputFile(): An error was encountered when opening " << fileName << "." << std::endl << "  Check that the file exists and that the file is not already opened.";
		SignalExceptionParse(ss.str(), 0);		
	}

	// Firstly, see what type of file this is, then decide what to do with it
	char first_chars[PRINT_LINE_LENGTH+1];

	std::string fileEpsg, fileEpoch;
	fileEpsg = "";
	
	try
	{
		// Obtain exclusive use of the input file pointer
		import_file_mutex.lock();

		ifsInputFILE_->get(first_chars, PRINT_LINE_LENGTH, '\n');
		ifsInputFILE_->seekg(0, std::ios::beg);				// put back to beginning

		// release file pointer mutex
		import_file_mutex.unlock();
	}
	catch (const std::ios_base::failure& f) {	
		ss.str("");
		ss << "ParseInputFile(): An error was encountered when reading " << fileName << "." << std::endl << "  " << f.what() << std::endl;
		SignalExceptionParse(ss.str(), 0);
	}
	
	/////////////////////////////////////////
	// Test file type and handle accordingly
	//
	// XML file
	if (strncmp(first_chars, "<?xml", 5) == 0)
	{
		// Set the file type
		input_file_meta->filetype = dynaml;
		m_ift = dynaml;
		
		// Parse the DynaML file
		ParseXML(fileName, vStations, stnCount, vMeasurements, msrCount, clusterID, fileEpsg, fileEpoch, firstFile, success_msg);

		if (fileEpsg.empty())
			fileEpsg = m_strProjectDefaultEpsg;
		
		// record the file's default reference frame
        snprintf(input_file_meta->epsgCode, sizeof(input_file_meta->epsgCode), "%s", fileEpsg.substr(0, STN_EPSG_WIDTH).c_str());
        snprintf(input_file_meta->epoch, sizeof(input_file_meta->epoch), "%s", fileEpoch.substr(0, STN_EPOCH_WIDTH).c_str());
	}
	// SNX
	else if (strncmp(first_chars, "%=SNX", 5) == 0)
	{
		// Set the file type
		input_file_meta->filetype = sinex;
		m_ift = sinex;
		m_idt = stn_msr_data;

		// Parse the SINEX file and capture the epoch
		ParseSNX(fileName, vStations, stnCount, vMeasurements, msrCount, clusterID);

		// Since SINEX files do not permit recording of reference frame within the file, set
		// the frame to the default reference frame
        snprintf(input_file_meta->epsgCode, sizeof(input_file_meta->epsgCode), "%s", m_strProjectDefaultEpsg.substr(0, STN_EPSG_WIDTH).c_str());
        snprintf(input_file_meta->epoch, sizeof(input_file_meta->epoch), "%s", datum_.GetEpoch_s().substr(0, STN_EPOCH_WIDTH).c_str());

		if (firstFile)
		{
			fileEpoch = input_file_meta->epoch;
			if (!fileEpoch.empty())
				m_strProjectDefaultEpoch = fileEpoch;
		}

		SignalComplete();
	}
	// STN or MSR
	else if (
		// use boost::algorithm::ifind_first, which is a case insensitive implementation of the find first algorithm. 
		strncmp(first_chars, "!#=DNA", 6) == 0 ||	// dna file?
		icontains(fileName, ".stn") ||			// dna station file
		icontains(fileName, ".msr"))				// dna measurement file
	{
		// Set the file type
		input_file_meta->filetype = dna;
		m_ift = dna;

		// Parse the DNA file
		ParseDNA(fileName, vStations, stnCount, vMeasurements, msrCount, clusterID, fileEpsg, fileEpoch, firstFile);
		
		if (fileEpsg.empty())
			fileEpsg = m_strProjectDefaultEpsg;

		// record the file's default reference frame
        snprintf(input_file_meta->epsgCode, sizeof(input_file_meta->epsgCode), "%s", fileEpsg.substr(0, STN_EPSG_WIDTH).c_str());
        snprintf(input_file_meta->epoch, sizeof(input_file_meta->epoch), "%s", fileEpoch.substr(0, STN_EPOCH_WIDTH).c_str());
		
		SignalComplete();
	}
	else
	{
		//throw XMLInteropException("ParseInputFile(): Could not deduce file type from extension or contents.", 0);		
		ss.str("");
		ss << "ParseInputFile(): " << leafStr<std::string>(fileName) << " is not a recognised station or" << std::endl <<
			"  measurement input file.";
		(*success_msg) = ss.str();
		//SignalExceptionParse(ss.str(), 0);
		parseStatus_ = PARSE_UNRECOGNISED_FILE;
		throw XMLInteropException(ss.str(), 0);
	}
	///////////////////////////////////

	// Apply discontinuities (if they exist) to each file except SINEX
	// SINEX files are automatically handled
	if (m_ift != sinex)
	{
		if (p->i.apply_discontinuities)
			ApplyDiscontinuities(vMeasurements);
	}
	
	// Populate metadata
    snprintf(input_file_meta->filename, sizeof(input_file_meta->filename), "%s", fileName.c_str());
	if (*stnCount > 0 && *msrCount > 0)
		input_file_meta->datatype = stn_msr_data;
	else if (*stnCount > 0)
		input_file_meta->datatype = stn_data;
	else if (*msrCount > 0)
		input_file_meta->datatype = msr_data;

	return parseStatus_;
}
	
void dna_import::ParseXML(const std::string& fileName, vdnaStnPtr* vStations, PUINT32 stnCount,
							   vdnaMsrPtr* vMeasurements, PUINT32 msrCount, PUINT32 clusterID,
							   std::string& fileEpsg, std::string& fileEpoch, bool firstFile, std::string* success_msg)
{
    // Lock before parsing the file
    import_file_mutex.lock();

	parseStatus_ = PARSE_SUCCESS;
	_filespecifiedreferenceframe = false;
	_filespecifiedepoch = false;

	// Check if DynaML.xsd exists in the current directory
	// This prevents the XML parser from hanging when the schema file is missing
	if (!std::filesystem::exists("DynaML.xsd"))
	{
		import_file_mutex.unlock();
		std::stringstream ss;
		ss << "ParseXML(): DynaML.xsd schema file not found in the current directory." << std::endl;
		ss << "  The XML parser requires this file to validate XML input files." << std::endl;
		ss << "  Please ensure DynaML.xsd is present in the working directory.";
		SignalExceptionParse(ss.str(), 0);
	}

	try
	{
		// Instantiate individual parsers.
		DnaXmlFormat_pimpl DnaXmlFormat_p(ifsInputFILE_,		// pass file stream to enable progress to be calculated
			clusterID,											// pass cluster ID so that a unique number can be retained across multiple files
			datum_.GetName(),									// pass the default reference frame
			datum_.GetEpoch_s(),								// pass the default epoch
			firstFile,											// is this the first file to be loaded?
			projectSettings_.i.user_supplied_frame==1,			// Has a reference frame been supplied?
			projectSettings_.i.user_supplied_epoch == 1,		// Has an epoch been supplied?
			projectSettings_.i.override_input_rfame==1);		// Should this reference frame override all others?
																		
		::DnaStation_pimpl DnaStation_p;
		::xml_schema::string_pimpl string_p;
		::StationCoord_pimpl StationCoord_p;
		::Height_pimpl Height_p;
		::DnaMeasurement_pimpl DnaMeasurement_p;
		::Directions_pimpl Directions_p;
		::GeoidModel_pimpl GeoidModel_p;
		::GPSBaseline_pimpl GPSBaseline_p;
		::GPSCovariance_pimpl GPSCovariance_p;
		::Clusterpoint_pimpl Clusterpoint_p;
		::PointCovariance_pimpl PointCovariance_p;
		::type_pimpl type_p;
		::referenceframe_pimpl referenceframe_p;
		::epoch_pimpl epoch_p;

		// Connect the parsers together.
		//
		DnaXmlFormat_p.parsers (DnaStation_p, DnaMeasurement_p, 
			// attributes
			type_p, referenceframe_p, epoch_p);

		DnaStation_p.parsers (string_p, string_p, string_p, StationCoord_p, string_p);

		StationCoord_p.parsers (string_p, string_p, string_p, Height_p, string_p, GeoidModel_p);

		DnaMeasurement_p.parsers (string_p, string_p, string_p, string_p, string_p, string_p, string_p, string_p,
				string_p, string_p, Directions_p, string_p, string_p, string_p, GPSBaseline_p, string_p, string_p,
				string_p, Clusterpoint_p, string_p, string_p, string_p, string_p,
				(projectSettings_.i.prefer_single_x_as_g == TRUE ? true : false));

		Directions_p.parsers (string_p, string_p, string_p, string_p, string_p);

		GPSBaseline_p.parsers (string_p, string_p, string_p, string_p, string_p, string_p, string_p, string_p,
				 string_p, string_p, GPSCovariance_p);

		GPSCovariance_p.parsers (string_p, string_p, string_p, string_p, string_p, string_p,
				string_p, string_p, string_p);

		Clusterpoint_p.parsers (string_p, string_p, string_p, string_p, string_p, string_p, string_p, string_p,
				string_p, string_p, PointCovariance_p);

		PointCovariance_p.parsers (string_p, string_p, string_p, string_p, string_p, string_p,
				string_p, string_p, string_p);

		// Parse the XML document.
		//
		::xml_schema::document doc_p (DnaXmlFormat_p, "DnaXmlFormat");

		DnaXmlFormat_p.pre();
		doc_p.parse (*ifsInputFILE_);
		DnaXmlFormat_p.post_DnaXmlFormat (vStations, vMeasurements);

        // unlock after parsing
        import_file_mutex.unlock();

		SignalComplete();
		*clusterID = DnaXmlFormat_p.CurrentClusterID();
		*stnCount = DnaXmlFormat_p.NumStationsRead();
		*msrCount = DnaXmlFormat_p.NumMeasurementsRead();
		*success_msg = DnaXmlFormat_p.DnaXmlParseMessage() + "\n";
		_filespecifiedreferenceframe = DnaXmlFormat_p.filespecifiedreferenceframe();
		_filespecifiedepoch = DnaXmlFormat_p.filespecifiedepoch();

		try
		{
			// 1. Capture file epsg
			// Get the reference frame from the XML file (i.e. referenceframe attribute in DnaXmlFormat element)
			fileEpsg = DnaXmlFormat_p.FileEpsg();
			// Is this attribute value an empty string?  As long as a default value 
			// is specified in DynaML.xsd, this value will never be empty, unless the user
			// has inadvertently set in the xml file, e.g.:
			//  <DnaXmlFormat referenceframe="" ... >
			if (fileEpsg.empty())
				// No, so get the epsg code from the default datum 
				fileEpsg = datum_.GetEpsgCode_s();

			fileEpoch = DnaXmlFormat_p.FileEpoch();
			// Is this attribute value an empty string?  As long as a default value 
			// is specified in DynaML.xsd, this value will never be empty, unless the user
			// has inadvertently set in the xml file, e.g.:
			//  <DnaXmlFormat epoch="" ... >
			// Or, is the file reference frame a static datum?
			if (fileEpoch.empty() || isEpsgDatumStatic(LongFromString<UINT32>(fileEpsg)))
			{
				// Get the epoch of the nominated epsg (whether default or from the file)
				fileEpoch = referenceepochFromEpsgString<std::string>(fileEpsg);
			}

			// Has the user supplied a reference frame on the command line (--reference-frame)?
			if (projectSettings_.i.user_supplied_frame)
			{
				// A frame has been supplied. Now check for epoch.

				// Has the user supplied an epoch on the command line (--epoch)?
				if (!projectSettings_.i.user_supplied_epoch)
				{
					// No epoch supplied. Take the epoch from the first file.
					if (firstFile)
					{
						projectSettings_.i.epoch = fileEpoch;
						projectSettings_.r.epoch = fileEpoch;
						m_strProjectDefaultEpoch = fileEpoch;
						datum_.SetEpoch(fileEpoch);
					}
				}
			}
			else
			{
				// No frame supplied. Set the datum from the first file. 
				// If the datum field is blank in the first file, the default will be used.
				if (firstFile)
				{
					projectSettings_.i.reference_frame = datumFromEpsgString<std::string>(fileEpsg);
					projectSettings_.r.reference_frame = projectSettings_.i.reference_frame;
					m_strProjectDefaultEpsg = fileEpsg;

					// If the user has not provided a reference frame on the command line, then
					// take the epoch from the first file. In this case, do not take the epoch
					// from the command line (dnaimportwrapper will throw in this scenario anyway)
					projectSettings_.i.epoch = fileEpoch;
					projectSettings_.r.epoch = fileEpoch;
					m_strProjectDefaultEpoch = fileEpoch;
				}
			}
		}
		catch (...)
		{
			std::stringstream ss;
			ss << "The default input file reference frame \"" << referenceframe_p.str() << "\" is not recognised.";
			SignalExceptionParse(static_cast<std::string>(ss.str()), 0);
		}

		if (!vStations->empty() && !vMeasurements->empty())
			m_idt = stn_msr_data;
		else if (!vStations->empty())
			m_idt = stn_data;
		else if (!vMeasurements->empty())
			m_idt = msr_data;

	}
	catch (const std::ios_base::failure& f)
	{
		if (ifsInputFILE_->eof())
		{
			// release file pointer mutex
			import_file_mutex.unlock();
			return;
		}
		if (ifsInputFILE_->rdstate() & std::ifstream::eofbit)
		{
			// release file pointer mutex
			import_file_mutex.unlock();
			return;
		}
		std::stringstream ss;
		ss << "ParseXML(): An std::ios_base failure was encountered while parsing " << fileName << "." << std::endl << "  " << f.what();
		SignalExceptionParse(static_cast<std::string>(ss.str()), 0);
	}
	catch (const std::system_error& e)
	{
		if (ifsInputFILE_->eof())
		{
			// release file pointer mutex
			import_file_mutex.unlock();
			return;
		}
		if (ifsInputFILE_->rdstate() & std::ifstream::eofbit)
		{
			// release file pointer mutex
			import_file_mutex.unlock();
			return;
		}
		std::stringstream ss;
		ss << "ParseXML(): An std::ios_base failure was encountered while parsing " << fileName << "." << std::endl << "  " << e.what();
		SignalExceptionParse(static_cast<std::string>(ss.str()), 0);
	}
	catch (const XMLInteropException& e) 
	{
		std::stringstream ss;
		ss << "ParseXML(): An exception was encountered while parsing " << fileName << "." << std::endl << "  " << e.what();
		SignalExceptionParse(static_cast<std::string>(ss.str()), 0);
	}
	catch (const ::xml_schema::parsing& e)
	{
		std::stringstream ss("");
		ss << e.what();

		::xsd::cxx::parser::diagnostics<char>::const_iterator _it;		
		for (_it=e.diagnostics().begin(); _it!=e.diagnostics().end(); _it++)
		{
			ss << std::endl;
			ss << "  - line " << _it->line();
			ss << ", column " <<  _it->column();
			ss << ", severity " <<  _it->severity() << std::endl;
			ss << "  - " << _it->message();
		}
		SignalExceptionParse(ss.str(), 0);
	}
	catch (const ::xml_schema::exception& e)
	{
		std::stringstream ss;
		ss << "ParseXML(): An xml_schema exception was encountered while parsing " << fileName << "." << std::endl << "  " << e.what();
		SignalExceptionParse(static_cast<std::string>(ss.str()), 0);
	}
	catch (...)
	{
		std::stringstream ss;
		ss << "ParseXML(): An unknown error was encountered while parsing " << fileName << "." << std::endl;
		SignalExceptionParse(ss.str(), 0);	
	}

	if (parseStatus_ != PARSE_SUCCESS)
	{
		std::stringstream ss("");
		ss.str("");
		ss << "- Warning: Parse success code = " << PARSE_SUCCESS << "." << std::endl;
		(*success_msg) += ss.str();
	}
}

void dna_import::ParseSNX(const std::string& fileName, vdnaStnPtr* vStations, PUINT32 stnCount, 
							   vdnaMsrPtr* vMeasurements, PUINT32 msrCount, PUINT32 clusterID)
{
	try {
        std::lock_guard<std::mutex> lock(import_file_mutex);

		// Load sinex file and capture epoch.  Throws runtime_error on failure.
		DnaIoSnx snx;
		snx.ParseSinex(&ifsInputFILE_, fileName, vStations, stnCount, vMeasurements, msrCount, clusterID,
			g_parsestn_tally, g_parsemsr_tally, g_fileOrder, 
			datum_, projectSettings_.i.apply_discontinuities==1, &stn_discontinuities_, m_discontsSortedbyName,
			m_lineNo, m_columnNo, parseStatus_);
	}
	catch (const std::runtime_error& e) {
		SignalExceptionParse(e.what(), 0);
	}
	_filespecifiedepoch = true;
}
	
// Parse discontinuities and create discontinuity tuple
void dna_import::ParseDiscontinuities(const std::string& fileName)
{
	std::ifstream discont_file;
	try {
		// Open discontinuity file.  Throws runtime_error on failure.
		file_opener(discont_file, fileName, 
			std::ios::in, ascii, true);
	}
	catch (const std::runtime_error& e) {
		SignalExceptionInterop(e.what(), 0, NULL);
	}

	try {
		// Load discontinuity file.  Throws runtime_error on failure.
		DnaIoSnx snx;
		projectSettings_.i.apply_discontinuities = true;
		snx.ParseDiscontinuityFile(&discont_file, fileName, 
			&stn_discontinuities_, m_discontsSortedbyName, 
			m_lineNo, m_columnNo, parseStatus_);
	}
	catch (const std::runtime_error& e) {
		SignalExceptionParse(e.what(), 0);
	}

	discont_file.close();
}
	
// Station names in SINEX files are automatically renamed to account for discontinuities and added
// to vStations.  For non SINEX files, any renamed stations are tracked in stn_renamed_.  After
// all files have been loaded, this function is called to (1) find any discontinuity stations not in the 
// station list and (2) add a duplicate for it, renamed to the respective discontinuity site name.
void dna_import::AddDiscontinuityStations(vdnaStnPtr* vStations)
{
	_it_vdnastnptr _it_stn(vStations->begin());
	it_string_pair stn_renames_it(stn_renamed_.begin());

	dnaStnPtr stn_ptr;
	std::string stationName;
	
	UINT32 i, station_count(static_cast<UINT32>(vStations->size()));
	UINT32 station_index(station_count);

	std::sort(vStations->begin(), vStations->end(), CompareStationName<dnaStnPtr>());
	
	for (i=0; i<station_count; ++i)
	{
		stationName = vStations->at(i)->GetName();
	
		// For each occurrence of stationName in stn_renamed_, create a clone and add it to 
		// vStations so that the measurement stations in non SINEX files that are
		// renamed by ApplyDiscontinuitiesMeasurements (via TrackDiscontinuitySite) are 
		// supported by corresponding stations.
		while ((stn_renames_it = binary_search_index_pair(stn_renames_it, stn_renamed_.end(), stationName)) != stn_renamed_.end())
		{
			// If a sinex file has been loaded and stationName exists in the sinex file, stn_renames_it->second may 
			// already exist in vStations.  In this case, the following code will add a duplicate station.
			// For this reason, duplicates will be removed at the end of this function.
			// Whilst this might seem inefficient, the alternative (i.e. search for stn_renames_it->second
			// in vStations before push_back) is far less efficient.  Also, if the first occurrence of this station
			// is the correct station loaded from a sinex file, the (cloned) duplicate will be removed on account of
			// the duplicate having a larger file order, via SetfileOrder(++station_index) below.
			
			stn_ptr.reset(vStations->at(i)->clone());
			stn_ptr->SetOriginalName();
			stn_ptr->SetName(stn_renames_it->second);
			stn_ptr->SetfileOrder(++station_index);
			vStations->push_back(stn_ptr);
			
			stn_renames_it++;
		}

		stn_renames_it = stn_renamed_.begin();
	}

	// Remove duplicates which may have been added through the above process
	std::sort(vStations->begin(), vStations->end(), CompareStationName<dnaStnPtr>());
	_it_vdnastnptr _it_stn_newend = unique(vStations->begin(), vStations->end(), EqualStationName<dnaStnPtr>());
	if (_it_stn_newend != vStations->end())
		vStations->resize(_it_stn_newend - vStations->begin());

}
	

void dna_import::ApplyDiscontinuities(vdnaMsrPtr* vMeasurements)
{
	if (stn_discontinuities_.empty())
		return;

	if (!m_discontsSortedbyName)
	{
		std::sort(stn_discontinuities_.begin(), stn_discontinuities_.end(),
			CompareSiteTuplesByName<discontinuity_tuple>());
		m_discontsSortedbyName = true;
	}

	// If the only input file is a SINEX file, there isn't any need to apply discontinuities 
	// to a station file, since those stations will be introduced via the process of 
	// handling a discontinuity file. 
	// However, if other station and measurement file types are loaded, there may be a case 
	// where station discontinuities exist but have not been accounted for at this point.  
	// For this purpose, a call to ApplyDiscontinuitiesStations will be made after all files 
	// have been imported.
	// If a station becomes renamed through the process of applying discontinuities to 
	// measurements, then it will be flagged as unused (since it won't be found in the 
	// measurements).

	if (!vMeasurements->empty())
		ApplyDiscontinuitiesMeasurements(vMeasurements);
}

void dna_import::TrackDiscontinuitySite(const std::string& site, const std::string& site_renamed)
{
	if (!binary_search(stn_renamed_.begin(), stn_renamed_.end(), site_renamed, ComparePairSecond<std::string>()))
	{
		stn_renamed_.push_back(string_string_pair(site, site_renamed));
		std::sort(stn_renamed_.begin(), stn_renamed_.end(), ComparePairFirst<std::string>());
	}
}
	

// This function renames stations in each measurement based on a discontinuity file
// NOTE: A fundamental prerequisite for the proper renaming of discontinuity sites
// in this function is a valid epoch for each measurement!
void dna_import::ApplyDiscontinuitiesMeasurements(vdnaMsrPtr* vMeasurements)
{
	_it_vdiscontinuity_tuple _it_discont(stn_discontinuities_.begin());

	std::string site_renamed, stn1, stn2, stn3;
	
	boost::gregorian::date site_date;

	std::vector<CDnaDirection>* vdirns;
	std::vector<CDnaGpsBaseline>* vgpsBsls;
	std::vector<CDnaGpsPoint>* vgpsPnts;

	_it_vdnamsrptr _it_msr(vMeasurements->begin());

	// Update station names
	for (_it_msr = vMeasurements->begin(); _it_msr != vMeasurements->end(); _it_msr++)
	{
		// 1. Handle nested type measurements (D, G, X, Y) separately
		switch (_it_msr->get()->GetTypeC())
		{
		case 'G':	// GPS Baseline (treat as single-baseline cluster)
		case 'X':	// GPS Baseline cluster
			vgpsBsls = _it_msr->get()->GetBaselines_ptr();
			ApplyDiscontinuitiesMeasurements_GX(vgpsBsls);
			continue;
		case 'Y':	// GPS point cluster
			vgpsPnts = _it_msr->get()->GetPoints_ptr();
			ApplyDiscontinuitiesMeasurements_Y(vgpsPnts);
			continue;
		}

		// Check if an epoch been provided with this measurement
		if (_it_msr->get()->GetEpoch().empty())
			continue;

		// Capture the epoch of the measurement
		site_date = dateFromString<boost::gregorian::date>(_it_msr->get()->GetEpoch());

		// 2. Handle 'first' station for every measurement type
		stn1 = _it_msr->get()->GetFirst();

		_it_discont = stn_discontinuities_.begin();
		if ((_it_discont = binary_search_discontinuity_site(
			_it_discont,
			stn_discontinuities_.end(),
			stn1)) != stn_discontinuities_.end())
		{
			// Site found in discontinuity file!
			// Handle renaming, and add this to the list of 'renamed stations'
			if (rename_discont_station(_it_discont, stn1, site_date, site_renamed))
			{
				_it_msr->get()->SetFirst(site_renamed);
				TrackDiscontinuitySite(stn1, site_renamed);
			}
		}

		// Finished with single station measurements
		switch (_it_msr->get()->GetTypeC())
		{
		case 'H':	// Orthometric height
		case 'R':	// Ellipsoidal height
		case 'I':	// Astronomic latitude
		case 'J':	// Astronomic longitude
		case 'P':	// Geodetic latitude
		case 'Q':	// Geodetic longitude
			continue;
		}

		// 3. Handle 'second' station for every measurement type
		stn2 = _it_msr->get()->GetTarget();

		_it_discont = stn_discontinuities_.begin();
		if ((_it_discont = binary_search_discontinuity_site(
			_it_discont,
			stn_discontinuities_.end(),
			stn2)) != stn_discontinuities_.end())
		{
			// Site found in discontinuity file!
			// Handle renaming, and add this to the list of 'renamed stations'
			if (rename_discont_station(_it_discont, stn2, site_date, site_renamed))
			{
				_it_msr->get()->SetTarget(site_renamed);
				TrackDiscontinuitySite(stn2, site_renamed);
			}
		}

		// 3. Handle D measurement directions
		switch (_it_msr->get()->GetTypeC())
		{
		case 'D':	// Direction set
					// handle the directions in the set
			vdirns = _it_msr->get()->GetDirections_ptr();
			ApplyDiscontinuitiesMeasurements_D(vdirns, site_date);
			continue;
			// Finished with single station measurements
		case 'B':	// Geodetic azimuth
		case 'C':	// Chord dist
		case 'E':	// Ellipsoid arc
		case 'K':	// Astronomic azimuth
		case 'L':	// Level difference
		case 'M':	// MSL arc
		case 'S':	// Slope distance
		case 'V':	// Zenith distance
		case 'Z':	// Vertical angle
			continue;
		}

		// 4. multiple station measurements
		stn3 = _it_msr->get()->GetTarget2();

		_it_discont = stn_discontinuities_.begin();
		if ((_it_discont = binary_search_discontinuity_site(
			_it_discont,
			stn_discontinuities_.end(),
			stn3)) != stn_discontinuities_.end())
		{
			// Site found in discontinuity file!
			// Handle renaming, and add this to the list of 'renamed stations'
			if (_it_discont->discontinuity_exists)
			{
				// Rename station name according to date if it is a discontinuity site
				if (rename_discont_station(_it_discont, stn3, site_date, site_renamed))
				{
					_it_msr->get()->SetTarget2(site_renamed);
					TrackDiscontinuitySite(stn3, site_renamed);
				}
			}
		}		
	}
}
	

void dna_import::ApplyDiscontinuitiesMeasurements_GX(std::vector<CDnaGpsBaseline>* vGpsBaselines)
{
	std::vector< CDnaGpsBaseline >::iterator _it_msr(vGpsBaselines->begin());
	_it_vdiscontinuity_tuple _it_discont(stn_discontinuities_.begin());
	std::string site_renamed, stn1, stn2, stn3;
	boost::gregorian::date site_date;

	for (_it_msr = vGpsBaselines->begin();
		_it_msr != vGpsBaselines->end();
		_it_msr++)
	{
		// Capture the start date of the site
		site_date = dateFromString<boost::gregorian::date>(_it_msr->GetEpoch());

		// Station 1
		stn1 = _it_msr->GetFirst();

		_it_discont = stn_discontinuities_.begin();
		if ((_it_discont = binary_search_discontinuity_site(
			_it_discont,
			stn_discontinuities_.end(),
			stn1)) != stn_discontinuities_.end())
		{
			// Site found in discontinuity file!
			// Handle renaming, and add this to the list of 'renamed stations'
			if (rename_discont_station(_it_discont, stn1, site_date, site_renamed))
			{
				_it_msr->SetFirst(site_renamed);
				TrackDiscontinuitySite(stn1, site_renamed);
			}
		}

		// Station 2
		stn2 = _it_msr->GetTarget();

		_it_discont = stn_discontinuities_.begin();
		if ((_it_discont = binary_search_discontinuity_site(
			_it_discont,
			stn_discontinuities_.end(),
			stn2)) != stn_discontinuities_.end())
		{
			// Site found in discontinuity file!
			// Handle renaming, and add this to the list of 'renamed stations'
			if (rename_discont_station(_it_discont, stn2, site_date, site_renamed))
			{
				_it_msr->SetTarget(site_renamed);
				TrackDiscontinuitySite(stn2, site_renamed);
			}
		}
	}
}


void dna_import::ApplyDiscontinuitiesMeasurements_Y(std::vector<CDnaGpsPoint>* vGpsPoints)
{
	std::vector< CDnaGpsPoint >::iterator _it_msr(vGpsPoints->begin());
	_it_vdiscontinuity_tuple _it_discont(stn_discontinuities_.begin());
	std::string site_renamed, stn1, stn2, stn3;
	boost::gregorian::date site_date;

	for (_it_msr = vGpsPoints->begin();
		_it_msr != vGpsPoints->end();
		_it_msr++)
	{
		// Capture the start date of the site
		site_date = dateFromString<boost::gregorian::date>(_it_msr->GetEpoch());

		// Station 1
		stn1 = _it_msr->GetFirst();

		_it_discont = stn_discontinuities_.begin();
		if ((_it_discont = binary_search_discontinuity_site(
			_it_discont,
			stn_discontinuities_.end(),
			stn1)) != stn_discontinuities_.end())
		{
			// Site found in discontinuity file!
			// Handle renaming, and add this to the list of 'renamed stations'
			if (rename_discont_station(_it_discont, stn1, site_date, site_renamed))
			{
				_it_msr->SetFirst(site_renamed);
				TrackDiscontinuitySite(stn1, site_renamed);
			}
		}
	}
}


void dna_import::ApplyDiscontinuitiesMeasurements_D(std::vector<CDnaDirection>* vDirections, const boost::gregorian::date& site_date)
{
	std::vector< CDnaDirection >::iterator _it_msr(vDirections->begin());
	_it_vdiscontinuity_tuple _it_discont(stn_discontinuities_.begin());
	std::string site_renamed, stn1, stn2, stn3;

	for (_it_msr = vDirections->begin();
		_it_msr != vDirections->end();
		_it_msr++)
	{
		// Station 1
		stn1 = _it_msr->GetFirst();

		_it_discont = stn_discontinuities_.begin();
		if ((_it_discont = binary_search_discontinuity_site(
			_it_discont,
			stn_discontinuities_.end(),
			stn1)) != stn_discontinuities_.end())
		{
			// Site found in discontinuity file!
			// Handle renaming, and add this to the list of 'renamed stations'
			if (rename_discont_station(_it_discont, stn1, site_date, site_renamed))
			{
				_it_msr->SetFirst(site_renamed);
				TrackDiscontinuitySite(stn1, site_renamed);
			}
		}

		// Station 2
		stn2 = _it_msr->GetTarget();

		_it_discont = stn_discontinuities_.begin();
		if ((_it_discont = binary_search_discontinuity_site(
			_it_discont,
			stn_discontinuities_.end(),
			stn2)) != stn_discontinuities_.end())
		{
			// Site found in discontinuity file!
			// Handle renaming, and add this to the list of 'renamed stations'
			if (rename_discont_station(_it_discont, stn2, site_date, site_renamed))
			{
				_it_msr->SetTarget(site_renamed);
				TrackDiscontinuitySite(stn2, site_renamed);
			}
		}
	}
}
	

void dna_import::ParseDNA(const std::string& fileName, vdnaStnPtr* vStations, PUINT32 stnCount, 
							   vdnaMsrPtr* vMeasurements, PUINT32 msrCount, PUINT32 clusterID, 
							   std::string& fileEpsg, std::string& fileEpoch, bool firstFile)
{
	parseStatus_ = PARSE_SUCCESS;

	(*stnCount) = 0;
	(*msrCount) = 0;
	g_parsestn_tally.initialise();
	g_parsemsr_tally.initialise();

	std::string stn_file_type(".stn"), msr_file_type(".msr");
	std::string version, geoversion;
	
	INPUT_DATA_TYPE idt;
	UINT32 count;
	size_t pos = 0;

	m_lineNo = 1;
	m_columnNo = 0;

	dna_io_dna dnaFile;

	// Get header information
	try {
		// Read DNA header
		
		// Obtain exclusive use of the input file pointer
		import_file_mutex.lock();
		// Read the dna file header, and set the
		// reference frame based on the header and user preferences
		dnaFile.read_dna_header(ifsInputFILE_, version, idt,			
			datum_,											// project datum
			fileEpsg, fileEpoch, geoversion, count);
		// release file pointer mutex
		import_file_mutex.unlock();
	}
	catch (const std::runtime_error& e) {
		import_file_mutex.unlock();
		parseStatus_ = PARSE_EXCEPTION_RAISED;
		throw XMLInteropException(e.what(), 0);
	}

	// Set to defaults if reference frame and epoch fields are empty
	if (fileEpsg.empty())
		fileEpsg = datum_.GetEpsgCode_s();	
	if (fileEpoch.empty())
		fileEpoch = datum_.GetEpoch_s();

	// Is the file reference frame a static datum?
	if (isEpsgDatumStatic(LongFromString<UINT32>(fileEpsg)))
	{
		// Get the reference epoch of the nominated epsg
		fileEpoch = referenceepochFromEpsgString<std::string>(fileEpsg);
	}

	// Has the user supplied a reference frame on the command line (--reference-frame)?
	if (projectSettings_.i.user_supplied_frame)
	{
		// A frame has been supplied. Now check for epoch.

		// Has the user supplied an epoch on the command line (--epoch)?
		if (!projectSettings_.i.user_supplied_epoch)
		{
			// No epoch supplied. Take the epoch from the first file.
			if (firstFile)
			{
				projectSettings_.i.epoch = fileEpoch;
				projectSettings_.r.epoch = fileEpoch;
				m_strProjectDefaultEpoch = fileEpoch;
				datum_.SetEpoch(fileEpoch);
			}
		}
	}
	else
	{
		// No frame supplied. Set the datum from the first file. 
		// If the datum field is blank in the first file, the default will be used. See
		// dnaFile.read_dna_header(..) for details.
		if (firstFile)
		{
			projectSettings_.i.reference_frame = datumFromEpsgString<std::string>(fileEpsg);
			projectSettings_.r.reference_frame = projectSettings_.i.reference_frame;
			m_strProjectDefaultEpsg = fileEpsg;

			// If the user has not provided a reference frame on the command line, then
			// take the epoch from the first file.  In this case, do not take the epoch
			// from the command line (dnaimportwrapper will throw in this scenario anyway)
			projectSettings_.i.epoch = fileEpoch;
			projectSettings_.r.epoch = fileEpoch;
			m_strProjectDefaultEpoch = fileEpoch;

			InitialiseDatum(projectSettings_.i.reference_frame, projectSettings_.i.epoch);
		}
	}

	_filespecifiedreferenceframe = dnaFile.filespecifiedReferenceFrame();
	_filespecifiedepoch = dnaFile.filespecifiedEpoch();
		
	// Station file
	if (idt == stn_data ||
		(pos = fileName.find(stn_file_type, 0)) != std::string::npos)
	{
		// set the column positions and widths (using version obtained from the header)
		dsl_ = dnaFile.dna_stn_positions();
		dsw_ = dnaFile.dna_stn_widths();
		
		try {
			// Does the user want to override the datum provided in the file with the project
			// datum?
			if (projectSettings_.i.override_input_rfame)
				ParseDNASTN(vStations, stnCount, datum_.GetEpsgCode_s(), datum_.GetEpoch_s());
			else
				ParseDNASTN(vStations, stnCount, fileEpsg, fileEpoch);
			m_idt = stn_data;
		}
		catch (const std::ios_base::failure& f) {
			if (ifsInputFILE_->eof())
			{
				// release file pointer mutex
				import_file_mutex.unlock();
				return;
			}
			std::stringstream ss;
			ss << "ParseDNA(): An std::ios_base failure was encountered when attempting to read stations file  " << fileName << "." << std::endl << "  " << f.what();
			throw XMLInteropException(ss.str(), 0);
		}
		catch (const XMLInteropException& f)
		{
			std::stringstream ss;
			ss << "ParseInputFile(): An exception was encountered while parsing " << fileName << "." << std::endl;
			ss << "  - line " << m_lineNo;
			ss << ", column " << m_columnNo << std::endl;
			ss << "  - " << f.what();
			SignalExceptionParse(static_cast<std::string>(ss.str()), 0);
		}
		catch (...) {
			if (ifsInputFILE_->eof())
			{
				// release file pointer mutex
				import_file_mutex.unlock();
				return;
			}
			std::stringstream ss;
			ss << "ParseDNA(): An error was encountered when attempting to read stations file  " << fileName << ".";
			throw XMLInteropException(ss.str(), 0);
		}
	}
	else if (idt == msr_data ||
		(pos = fileName.find(msr_file_type, 0)) != std::string::npos)
	{
		// Determine the file format version
		//ParseDNAVersion(msr_data);
		dml_ = dnaFile.dna_msr_positions();
		dmw_ = dnaFile.dna_msr_widths();
		
		try {
			if (projectSettings_.i.override_input_rfame)
				ParseDNAMSR(vMeasurements, msrCount, clusterID, datum_.GetEpsgCode_s(), datum_.GetEpoch_s());
			else
				ParseDNAMSR(vMeasurements, msrCount, clusterID, fileEpsg, fileEpoch);
			m_idt = msr_data;
		}
		catch (const std::ios_base::failure& f) {
			if (ifsInputFILE_->eof())
			{
				// release file pointer mutex
				import_file_mutex.unlock();
				return;
			}
			std::stringstream ss;
			ss << "ParseDNA(): An std::ios_base failure was encountered when attempting to read measurements file  " << fileName << "." << std::endl << "  " << f.what();
			throw XMLInteropException(ss.str(), 0);
		}
		catch (const XMLInteropException& f)
		{
			std::stringstream ss;
			ss << "ParseInputFile(): An exception was encountered while parsing " << fileName << "." << std::endl;
			ss << "  - line " << m_lineNo;
			ss << ", column " << m_columnNo << std::endl;
			ss << "  - " << f.what();
			SignalExceptionParse(static_cast<std::string>(ss.str()), 0);
		}
		catch (...) {
			if (ifsInputFILE_->eof())
			{
				// release file pointer mutex
				import_file_mutex.unlock();
				return;
			}
			std::stringstream ss;
			ss << "ParseDNA(): An error was encountered when attempting to read measurements file  " << fileName << "." << std::endl;
			ss << "  - line " << m_lineNo;
			ss << ", column " <<  m_columnNo << std::endl;
			throw XMLInteropException(ss.str(), 0);
		}
	}
}
	
void dna_import::ParseDNASTN(vdnaStnPtr* vStations, PUINT32 stnCount, const std::string& epsg, const std::string& epoch)
{
	std::string sBuf, tmp;

	dnaStnPtr stn_ptr;
	vStations->clear();

	//while (!ifsInputFILE_->eof())			// while EOF not found
	while (ifsInputFILE_)
	{
		// Obtain exclusive use of the input file pointer
		import_file_mutex.lock();

		if (ifsInputFILE_->eof())
		{
			// release file pointer mutex
			import_file_mutex.unlock();
			break;
		}

		m_lineNo++;
		
		try {
			getline((*ifsInputFILE_), sBuf);
		}
		catch (...) {
			if (ifsInputFILE_->eof())
			{
				// release file pointer mutex
				import_file_mutex.unlock();
				return;
			}
			std::stringstream ss;
			ss << "ParseDNASTN(): Could not read from the station file." << std::endl;
			m_columnNo = 0;
			throw XMLInteropException(ss.str(), m_lineNo);
		}

		// release file pointer mutex
		import_file_mutex.unlock();
		
		// blank or whitespace?
		if (trimstr(sBuf).empty())			
			continue;

		// Ignore lines with blank station name
		if (trimstr(sBuf.substr(dsl_.stn_name, dsw_.stn_name)).empty())			
			continue;
		
		// Ignore lines with comments
		if (sBuf.compare(0, 1, "*") == 0)
			continue;
		
		// initialise new station
		stn_ptr.reset(new CDnaStation(datumFromEpsgString<std::string>(epsg), epoch));

		stn_ptr->SetfileOrder(g_fileOrder++);

		// name
		try {
			tmp = trimstr(sBuf.substr(dsl_.stn_name, dsw_.stn_name));	
			stn_ptr->SetName(tmp);
		}
		catch (...) {
			std::stringstream ss;
			ss << "ParseDNASTN(): Could not extract station name from the record:  " << std::endl << "    " << sBuf << std::endl;
			m_columnNo = dsl_.stn_name+1;
			throw XMLInteropException(ss.str(), m_lineNo);
		}

		// constraints
		try {
			tmp = trimstr(sBuf.substr(dsl_.stn_const, dsw_.stn_const));	
			stn_ptr->SetConstraints(tmp);
			g_parsestn_tally.addstation(tmp);
		}
		catch (...) {
			std::stringstream ss;
			ss << "ParseDNASTN(): Could not extract station constraints from the record:  " << std::endl << "    " << sBuf << std::endl;
			m_columnNo = dsl_.stn_const+1;
			throw XMLInteropException(ss.str(), m_lineNo);
		}

		// coordinate type
		try {
			tmp = trimstr(sBuf.substr(dsl_.stn_type, dsw_.stn_type));
			stn_ptr->SetCoordType(tmp);
		}
		catch (...) {
			std::stringstream ss;
			ss << "ParseDNASTN(): Could not extract coordinate type from the record:  " << std::endl << "    " << sBuf << std::endl;
			m_columnNo = dsl_.stn_type+1;
			throw XMLInteropException(ss.str(), m_lineNo);
		}

		// coordinates
		try {
			tmp = trimstr(sBuf.substr(dsl_.stn_e_phi_x, dsw_.stn_e_phi_x));		// easting, latitude, X
			stn_ptr->SetXAxis(tmp);
		}
		catch (...) {
			std::stringstream ss;
			ss << "ParseDNASTN(): Could not extract station ";
			switch (stn_ptr.get()->GetMyCoordTypeC())
			{
			case XYZ_type_i:
				ss << "X";
				break;
			case UTM_type_i:
				ss << "easting";
				break;
			case LLH_type_i:
			case LLh_type_i:
				ss << "latitude";
				break;
			default:
				break;
			}
			ss << " value from the record:  " << std::endl << "    " << sBuf << std::endl;
			m_columnNo = dsl_.stn_e_phi_x+1;
			throw XMLInteropException(ss.str(), m_lineNo);
		}

		try {
			tmp = trimstr(sBuf.substr(dsl_.stn_n_lam_y, dsw_.stn_n_lam_y));		// northing, longitude, Y
			stn_ptr->SetYAxis(tmp);
		}
		catch (...) {
			std::stringstream ss;
			ss << "ParseDNASTN(): Could not extract station ";
			switch (stn_ptr->GetMyCoordTypeC())
			{
			case XYZ_type_i:
				ss << "Y";
				break;
			case UTM_type_i:
				ss << "northing";
				break;
			case LLH_type_i:
			case LLh_type_i:
				ss << "longitude";
				break;
			default:
				break;
			}
			ss << " value from the record:  " << std::endl << "    " << sBuf << std::endl;
			m_columnNo = dsl_.stn_n_lam_y+1;
			throw XMLInteropException(ss.str(), m_lineNo);
		}

		try {
			tmp = trimstr(sBuf.substr(dsl_.stn_ht_z, dsw_.stn_ht_z));		// orthometric height, Z
			if (stn_ptr->GetMyCoordTypeC() == XYZ_type_i)
				stn_ptr->SetZAxis(tmp);
			else
				stn_ptr->SetHeight(tmp);			
		}
		catch (...) {
			std::stringstream ss;
			ss << "ParseDNASTN(): Could not extract station ";
			switch (stn_ptr->GetMyCoordTypeC())
			{
			case XYZ_type_i:
				ss << "Z";
				break;
			case LLH_type_i:
			case LLh_type_i:
			case UTM_type_i:
				ss << "height";
				break;
			default:
				break;
			}
			ss << " value from the record:  " << std::endl << "    " << sBuf << std::endl;
			m_columnNo = dsl_.stn_ht_z+1;
			throw XMLInteropException(ss.str(), m_lineNo);
		}

		if (sBuf.length() > dsl_.stn_hemi_zo)
		{
			try {
				tmp = trimstr(sBuf.substr(dsl_.stn_hemi_zo, dsw_.stn_hemi_zo));		// hemisphere-zone
				stn_ptr->SetHemisphereZone(tmp);
			}
			catch (...) {
				switch (stn_ptr->GetMyCoordTypeC())
				{
				case LLh_type_i:
				case LLH_type_i:
				case ENU_type_i:
				case AED_type_i:
				case XYZ_type_i:
					break;
				case UTM_type_i:	// Hemisphere and zone is only essential for UTM types
					std::stringstream ss;
					ss << "ParseDNASTN(): Could not extract station hemisphere and zone from the record:  " << std::endl << "    " << sBuf << std::endl;
					m_columnNo = dsl_.stn_hemi_zo+1;
					throw XMLInteropException(ss.str(), m_lineNo);
				}
			}
		}

		if (sBuf.length() > dsl_.stn_desc)
		{
			try {
				tmp = trimstr(sBuf.substr(dsl_.stn_desc));		// description
				stn_ptr->SetDescription(tmp);
			}
			catch (...) {		// do nothing (description is not compulsory)
			}
		}

		(*stnCount)++;
		vStations->push_back(stn_ptr);
	}
	
}
	

void dna_import::ParseDNAMSR(pvdnaMsrPtr vMeasurements, PUINT32 msrCount, PUINT32 clusterID, const std::string& fileEpsg, const std::string& fileEpoch)
{
	std::string sBuf, tmp;

	dnaMsrPtr msr_ptr;
	vMeasurements->clear();

	char cType;
	bool ignoreMsr;
	bool measurementRead(false);

	m_msrComments.clear();

	//while (!ifsInputFILE_->eof())			// while EOF not found
	while (ifsInputFILE_)
	{
		// Obtain exclusive use of the input file pointer
		import_file_mutex.lock();

		if (ifsInputFILE_->eof())
		{
			// release file pointer mutex
			import_file_mutex.unlock();
			break;
		}

		m_lineNo++;
		
		try {
			getline((*ifsInputFILE_), sBuf);
		}
		catch (...) {
			if (ifsInputFILE_->eof())
			{
				// release file pointer mutex
				import_file_mutex.unlock();
				return;		
			}
			std::stringstream ss;
			ss << "ParseDNAMSR(): Could not read from the measurement file." << std::endl;
			m_columnNo = 0;
			throw XMLInteropException(ss.str(), m_lineNo);
		}
		
		// release file pointer mutex
		import_file_mutex.unlock();
		
		// blank or whitespace?
		if (trimstr(sBuf).empty())			
			continue;
		
		// one character (most likely '*') line
		if (trimstr(sBuf).length() < 2)
			continue;

		// no station value?
		if (trimstr(sBuf.substr(dml_.msr_inst, dmw_.msr_inst)).empty())			
			continue;

		// Capture comment (which may apply to several measurements)
		if (sBuf.compare(0, 1, "*") == 0)
		{
			// Is this a new comment?
			if (measurementRead)
				// The last record was a measurement, so this is a new comment
				m_msrComments.clear();

			measurementRead = false;

			if (!m_msrComments.empty())
				m_msrComments += "\n";
			m_msrComments += sBuf.substr(1);
			continue;
		}		
		
		// no station value?
		if (trimstr(sBuf.substr(dml_.msr_inst, dmw_.msr_inst)).empty())			
			continue;

		try {
			tmp = trimstr(sBuf.substr(dml_.msr_type, 1));
			cType = (tmp.c_str())[0];
			cType = static_cast<char>(toupper(cType));
		}
		catch (...) {
			std::stringstream ss;
			ss << "ParseDNAMSR(): Could not extract measurement type from the record:  " << std::endl << "    " << sBuf << std::endl;
			m_columnNo = dml_.msr_type+1;
			throw XMLInteropException(ss.str(), m_lineNo);
		}

		ignoreMsr = iequals("*", sBuf.substr(dml_.msr_ignore, dmw_.msr_ignore));

		switch (cType)
		{
		case 'A': // Horizontal angle
			g_parsemsr_tally.A++;
			msr_ptr.reset(new CDnaAngle);
			ParseDNAMSRAngular(sBuf, msr_ptr);
			(*msrCount) += 1;
			break;
		case 'B': // Geodetic azimuth
			g_parsemsr_tally.B++;
			msr_ptr.reset(new CDnaAzimuth);
			ParseDNAMSRAngular(sBuf, msr_ptr);
			(*msrCount) += 1;
			break;
		case 'C': // Chord dist
			g_parsemsr_tally.C++;
			msr_ptr.reset(new CDnaDistance);
			ParseDNAMSRLinear(sBuf, msr_ptr);
			(*msrCount) += 1;
			break;
		case 'D': // Direction set
			msr_ptr.reset(new CDnaDirectionSet(++(*clusterID)));
			if (ParseDNAMSRDirections(sBuf, msr_ptr, ignoreMsr) == 0)
			{
				measurementRead = true;
				continue;
			}
			(*msrCount) += static_cast<UINT32>(msr_ptr->GetDirections_ptr()->size());
			break;
		case 'E': // Ellipsoid arc
			g_parsemsr_tally.E++;
			msr_ptr.reset(new CDnaDistance);
			ParseDNAMSRLinear(sBuf, msr_ptr);
			(*msrCount) += 1;
			break;
		case 'G': // GPS Baseline (treat as single-baseline cluster)
		case 'X': // GPS Baseline cluster
			// Default to the fileEpsg and fileEpoch (see read_dna_header(..) in ParseDNA)
			msr_ptr.reset(new CDnaGpsBaselineCluster(++(*clusterID), datumFromEpsgString<std::string>(fileEpsg), fileEpoch));
			ParseDNAMSRGPSBaselines(sBuf, msr_ptr, ignoreMsr);
			(*msrCount) += static_cast<UINT32>(msr_ptr->GetBaselines_ptr()->size() * 3);
			break;
		case 'H': // Orthometric height
			g_parsemsr_tally.H++;
			msr_ptr.reset(new CDnaHeight);
			ParseDNAMSRLinear(sBuf, msr_ptr);
			(*msrCount) += 1;
			break;
		case 'I': // Astronomic latitude
			g_parsemsr_tally.I++;
			msr_ptr.reset(new CDnaCoordinate);
			ParseDNAMSRCoordinate(sBuf, msr_ptr);
			(*msrCount) += 1;
			break;
		case 'J': // Astronomic longitude
			g_parsemsr_tally.J++;
			msr_ptr.reset(new CDnaCoordinate);
			ParseDNAMSRCoordinate(sBuf, msr_ptr);
			(*msrCount) += 1;
			break;
		case 'K': // Astronomic azimuth
			g_parsemsr_tally.K++;
			msr_ptr.reset(new CDnaAzimuth);
			ParseDNAMSRAngular(sBuf, msr_ptr);
			(*msrCount) += 1;
			break;
		case 'L': // Level difference
			g_parsemsr_tally.L++;
			msr_ptr.reset(new CDnaHeightDifference);
			ParseDNAMSRLinear(sBuf, msr_ptr);
			(*msrCount) += 1;
			break;
		case 'M': // MSL arc
			g_parsemsr_tally.M++;
			msr_ptr.reset(new CDnaDistance);
			ParseDNAMSRLinear(sBuf, msr_ptr);
			(*msrCount) += 1;
			break;
		case 'P': // Geodetic latitude
			g_parsemsr_tally.P++;
			msr_ptr.reset(new CDnaCoordinate);
			ParseDNAMSRCoordinate(sBuf, msr_ptr);
			(*msrCount) += 1;
			break;
		case 'Q': // Geodetic longitude
			g_parsemsr_tally.Q++;
			msr_ptr.reset(new CDnaCoordinate);
			ParseDNAMSRCoordinate(sBuf, msr_ptr);
			(*msrCount) += 1;
			break;
		case 'R': // Ellipsoidal height
			g_parsemsr_tally.R++;
			msr_ptr.reset(new CDnaHeight);
			ParseDNAMSRLinear(sBuf, msr_ptr);
			(*msrCount) += 1;
			break;
		case 'S': // Slope distance
			g_parsemsr_tally.S++;
			msr_ptr.reset(new CDnaDistance);
			ParseDNAMSRLinear(sBuf, msr_ptr);
			(*msrCount) += 1;
			break;
		case 'V': // Zenith distance
			g_parsemsr_tally.V++;
			msr_ptr.reset(new CDnaDirection);
			ParseDNAMSRAngular(sBuf, msr_ptr);
			(*msrCount) += 1;
			break;
		case 'Y': // GPS point cluster
			// Default to the fileEpsg and fileEpoch (see read_dna_header(..) in ParseDNA)
			msr_ptr.reset(new CDnaGpsPointCluster(++(*clusterID), datumFromEpsgString<std::string>(fileEpsg), fileEpoch));
			ParseDNAMSRGPSPoints(sBuf, msr_ptr, ignoreMsr);
			(*msrCount) += static_cast<UINT32>(msr_ptr->GetPoints_ptr()->size() * 3);
			break;
		case 'Z': // Vertical angle
			g_parsemsr_tally.Z++;
			msr_ptr.reset(new CDnaDirection);
			ParseDNAMSRAngular(sBuf, msr_ptr);
			(*msrCount) += 1;
			break;
		default:
			std::stringstream ss;
			ss << "ParseDNAMSR(): Unknown measurement type:  " << cType << std::endl;
			throw XMLInteropException(ss.str(), m_lineNo); 
		}

		measurementRead = true;
		
		// set ignore flag
		msr_ptr->SetIgnore(ignoreMsr);

		// add the measurement to the list
		vMeasurements->push_back(msr_ptr);
	}
}
	
void dna_import::ParseDNAMSRLinear(const std::string& sBuf, dnaMsrPtr& msr_ptr)
{
	// Measurement type
	try {
		msr_ptr->SetType(trimstr(sBuf.substr(dml_.msr_type, dmw_.msr_type)));
	}
	catch (...) {
		SignalExceptionParseDNA("ParseDNAMSRLinear(): Could not extract measurement type from the record:  ",
			sBuf, dml_.msr_type);
	}

	// Instrument station
	msr_ptr->SetFirst(ParseInstrumentValue(sBuf, "ParseDNAMSRLinear"));

	// Target station
	if (msr_ptr->GetMsrStnCount() >= TWO_STATION)
		msr_ptr->SetTarget(ParseTargetValue(sBuf, "ParseDNAMSRLinear"));

	// If import is being used to simulate measurements, return as the
	// expected input is:
	//	 Measurement
	//	 From station
	//	 To station
	// So no need to read any further
	if (projectSettings_.i.simulate_measurements)
		return;

	// Value
	msr_ptr->SetValue(ParseLinearValue(sBuf, 
		measurement_name<char, std::string>(msr_ptr->GetTypeC()), "ParseDNAMSRLinear"));

	// Standard deviation
	msr_ptr->SetStdDev(ParseStdDevValue(sBuf, "ParseDNAMSRLinear"));

	// Epoch
	msr_ptr->SetEpoch(ParseEpochValue(sBuf, "ParseDNAMSRLinear"));

	// Capture msr_id and cluster_id (for database referencing)
	ParseDatabaseIds(sBuf, "ParseDNAMSRLinear", msr_ptr->GetTypeC());
	msr_ptr->SetDatabaseMap(m_msr_db_map);

	// instrument and target heights only make sense for 
	// slope distances, vertical angles and zenith distances
	switch (msr_ptr->GetTypeC())
	{
	case 'S':
		break;
	default:
		// don't need instrument height for all other measurement types
		return;
	}

	// Instrument and target heights
	if (sBuf.length() <= dml_.msr_inst_ht)
	{
		msr_ptr->SetInstrumentHeight("0.0");
		msr_ptr->SetTargetHeight("0.0");
		return;
	}
	msr_ptr->SetInstrumentHeight(ParseInstHeightValue(sBuf, "ParseDNAMSRLinear"));

	// Target height
	if (sBuf.length() <= dml_.msr_targ_ht)
	{
		msr_ptr->SetTargetHeight("0.0");
		return;
	}
	msr_ptr->SetTargetHeight(ParseTargHeightValue(sBuf, "ParseDNAMSRLinear"));
}
	

void dna_import::ParseDNAMSRCoordinate(const std::string& sBuf, dnaMsrPtr& msr_ptr)
{
	// Measurement type
	try {
		msr_ptr->SetType(trimstr(sBuf.substr(dml_.msr_type, dmw_.msr_type)));
	}
	catch (...) {
		SignalExceptionParseDNA("ParseDNAMSRCoordinate(): Could not extract measurement type from the record:  ",
			sBuf, dml_.msr_type);
	}

	// Instrument station
	msr_ptr->SetFirst(ParseInstrumentValue(sBuf, "ParseDNAMSRCoordinate"));
	
	// If import is being used to simulate measurements, return as the
	// expected input is:
	//	 Measurement
	//	 From station
	//	 To station
	// So no need to read any further
	if (projectSettings_.i.simulate_measurements)
		return;

	// Value
	msr_ptr->SetValue(ParseLinearValue(sBuf, 
		measurement_name<char, std::string>(msr_ptr->GetTypeC()), "ParseDNAMSRCoordinate"));

	// Standard deviation
	msr_ptr->SetStdDev(ParseStdDevValue(sBuf, "ParseDNAMSRCoordinate"));

	// Epoch
	msr_ptr->SetEpoch(ParseEpochValue(sBuf, "ParseDNAMSRLinear"));

	// Capture msr_id and cluster_id (for database referencing), then set
	// database id info
	ParseDatabaseIds(sBuf, "ParseDNAMSRCoordinate", msr_ptr->GetTypeC());
	msr_ptr->SetDatabaseMap(m_msr_db_map);
}
	

void dna_import::ParseDNAMSRGPSBaselines(std::string& sBuf, dnaMsrPtr& msr_ptr, bool ignoreMsr)
{
	CDnaGpsBaseline bslTmp;
	CDnaCovariance covTmp;

	bslTmp.SetReferenceFrame(msr_ptr->GetReferenceFrame());
	bslTmp.SetEpoch(msr_ptr->GetEpoch());

	// Measurement type
	std::string tmp;
	try {
		tmp = trimstr(sBuf.substr(dml_.msr_type, dmw_.msr_type));
	}
	catch (...) {
		SignalExceptionParseDNA("ParseDNAMSRGPSBaselines(): Could not extract measurement type from the record:  ",
			sBuf, dml_.msr_type);
	}
		
	msr_ptr->SetType(tmp);
	bslTmp.SetType(tmp);
	bslTmp.SetIgnore(ignoreMsr);
	covTmp.SetType(tmp);
	covTmp.SetIgnore(ignoreMsr);

	// So no need to read database ID
	if (!projectSettings_.i.simulate_measurements)
	{
		// Capture msr_id and cluster_id (for database referencing)
		ParseDatabaseIds(sBuf, "ParseDNAMSRGPSBaselines", msr_ptr->GetTypeC());
		msr_ptr->SetDatabaseMap(m_msr_db_map);
		bslTmp.SetDatabaseMap(m_msr_db_map);
		bslTmp.SetClusterID(msr_ptr->GetClusterID());
		covTmp.SetClusterID(msr_ptr->GetClusterID());
	}

	// Instrument station
	msr_ptr->SetFirst(ParseInstrumentValue(sBuf, "ParseDNAMSRGPSBaselines"));
	bslTmp.SetFirst(msr_ptr->GetFirst());

	// Target station
	bslTmp.SetTarget(ParseTargetValue(sBuf, "ParseDNAMSRGPSBaselines"));
	
	// Number of baselines
	UINT32 bslCount(1);
	if (iequals(msr_ptr->GetType(), "X"))
		msr_ptr->SetTotal(ParseMsrCountValue(sBuf, bslCount, "ParseDNAMSRGPSBaselines"));
	msr_ptr->SetRecordedTotal(bslCount);
	bslTmp.SetRecordedTotal(bslCount);

	msr_ptr->GetBaselines_ptr()->reserve(bslCount);
	if (iequals(msr_ptr->GetType(), "X"))
		g_parsemsr_tally.X += bslCount * 3;
	else
		g_parsemsr_tally.G += bslCount * 3;

	// V-scale
	msr_ptr->SetVscale(ParseScaleVValue(sBuf, "ParseDNAMSRGPSBaselines"));
	bslTmp.SetVscale(msr_ptr->GetVscale());

	// P-scale
	msr_ptr->SetPscale(ParseScalePValue(sBuf, "ParseDNAMSRGPSBaselines"));
	bslTmp.SetPscale(msr_ptr->GetPscale());

	// L-scale supplied?
	msr_ptr->SetLscale(ParseScaleLValue(sBuf, "ParseDNAMSRGPSBaselines"));
	bslTmp.SetLscale(msr_ptr->GetLscale());
	
	// H-scale supplied?
	msr_ptr->SetHscale(ParseScaleHValue(sBuf, "ParseDNAMSRGPSBaselines"));
	bslTmp.SetHscale(msr_ptr->GetHscale());

	// Process for setting reference frame is:
	// 1. Set the default reference frame in ParseDNAMSR when msr_ptr is initialised
	// 2. Initialise bslTmp using msr_ptr
	// 3. When a reference frame element is found here, change the reference frame
	//    using checks on override
	// same as:
	//	- DnaXmlFormat_pskel::_start_element_impl when DnaMeasurement is reached
	//	- DnaMeasurement_pimpl::ReferenceFrame when ReferenceFrame tag is found
	//	- GPSBaseline_pimpl::pre()
	if (projectSettings_.i.override_input_rfame)
	{
		// Set the cluster frame
		msr_ptr->SetReferenceFrame(projectSettings_.i.reference_frame);
		// Set the baseline frame
		bslTmp.SetReferenceFrame(projectSettings_.i.reference_frame);
	}
	else //if (!projectSettings_.i.override_input_rfame)
	{
		try {
			// Get the reference frame from the file
			tmp = ParseRefFrameValue(sBuf, "ParseDNAMSRGPSBaselines");
		
			if (!tmp.empty())
			{
				// Set the cluster frame
				msr_ptr->SetReferenceFrame(tmp);
				// Set the baseline frame
				bslTmp.SetReferenceFrame(tmp);
			}
		}
		catch (std::runtime_error& e) {
			std::stringstream ss;
			ss << "ParseDNAMSRGPSBaselines(): Error parsing reference frame:  " << std::endl << 
				"    " << e.what();
			SignalExceptionParseDNA(ss.str(), "", dml_.msr_gps_reframe);
		}

		try {
			// Get the epoch from the file
			tmp = ParseEpochValue(sBuf, "ParseDNAMSRGPSBaselines");

			if (tmp.empty())
			{
				if (isEpsgWGS84Ensemble(epsgCodeFromName<UINT32, std::string>(msr_ptr->GetReferenceFrame())))
				{
					// Set the cluster epoch
					msr_ptr->SetEpoch("");
					// Set the baseline epoch
					bslTmp.SetEpoch("");
				}
			}
			else
			{
				// Set the cluster epoch
				msr_ptr->SetEpoch(tmp);
				// Set the baseline epoch
				bslTmp.SetEpoch(tmp);
			}
		}
		catch (std::runtime_error& e) {
			std::stringstream ss;
			ss << "ParseDNAMSRGPSBaselines(): Error parsing epoch:  " << 
			 	std::endl << "    " << e.what();
			SignalExceptionParseDNA(ss.str(), "", dml_.msr_gps_epoch);
		}
	}
	
	bool first_run(true);
	UINT32 cov, covCount(bslCount - 1);

	for (UINT32 b(0); b<bslCount; ++b, --covCount)
	{
		if (!first_run)
		{
			m_lineNo++;
			
			// Obtain exclusive use of the input file pointer
			import_file_mutex.lock();
			getline((*ifsInputFILE_), sBuf);
			// release file pointer mutex
			import_file_mutex.unlock();

			// Instrument station
			msr_ptr->SetFirst(ParseInstrumentValue(sBuf, "ParseDNAMSRGPSBaselines"));
			bslTmp.SetFirst(msr_ptr->GetFirst());

			// Target station
			bslTmp.SetTarget(ParseTargetValue(sBuf, "ParseDNAMSRGPSBaselines"));

			// No need to read database ID for simulation
			if (!projectSettings_.i.simulate_measurements)
			{
				// Capture msr_id (for database referencing)
				ParseDatabaseIds(sBuf, "ParseDNAMSRGPSBaselines", msr_ptr->GetTypeC());
				bslTmp.SetDatabaseMap(m_msr_db_map);
			}
		}

		first_run = false;

		if (!projectSettings_.i.simulate_measurements)
		{
			m_lineNo++;
			// Obtain exclusive use of the input file pointer
			import_file_mutex.lock();
			getline((*ifsInputFILE_), sBuf);
			// release file pointer mutex
			import_file_mutex.unlock();

			bslTmp.SetX(ParseGPSMsrValue(sBuf, "X", "ParseDNAMSRGPSBaselines"));
			bslTmp.SetSigmaXX(ParseGPSVarValue(sBuf, "X", dml_.msr_gps_vcv_1, dmw_.msr_gps_vcv_1, "ParseDNAMSRGPSBaselines"));

			m_lineNo++;
			// Obtain exclusive use of the input file pointer
			import_file_mutex.lock();
			getline((*ifsInputFILE_), sBuf);
			// release file pointer mutex
			import_file_mutex.unlock();
	
			bslTmp.SetY(ParseGPSMsrValue(sBuf, "Y", "ParseDNAMSRGPSBaselines"));
			bslTmp.SetSigmaXY(ParseGPSVarValue(sBuf, "Y", dml_.msr_gps_vcv_1, dmw_.msr_gps_vcv_1, "ParseDNAMSRGPSBaselines"));
			bslTmp.SetSigmaYY(ParseGPSVarValue(sBuf, "Y", dml_.msr_gps_vcv_2, dmw_.msr_gps_vcv_2, "ParseDNAMSRGPSBaselines"));

			m_lineNo++;
			// Obtain exclusive use of the input file pointer
			import_file_mutex.lock();
			getline((*ifsInputFILE_), sBuf);
			// release file pointer mutex
			import_file_mutex.unlock();
	
			bslTmp.SetZ(ParseGPSMsrValue(sBuf, "Z", "ParseDNAMSRGPSBaselines"));
			bslTmp.SetSigmaXZ(ParseGPSVarValue(sBuf, "Z", dml_.msr_gps_vcv_1, dmw_.msr_gps_vcv_1, "ParseDNAMSRGPSBaselines"));
			bslTmp.SetSigmaYZ(ParseGPSVarValue(sBuf, "Z", dml_.msr_gps_vcv_2, dmw_.msr_gps_vcv_2, "ParseDNAMSRGPSBaselines"));
			bslTmp.SetSigmaZZ(ParseGPSVarValue(sBuf, "Z", dml_.msr_gps_vcv_3, dmw_.msr_gps_vcv_3, "ParseDNAMSRGPSBaselines"));
		}
		
		bslTmp.ResizeGpsCovariancesCount();
		bslTmp.ReserveGpsCovariancesCount(covCount);

		for (cov=0; cov<covCount; ++cov)
		{
			if (!projectSettings_.i.simulate_measurements)
				ParseDNAMSRCovariance(covTmp);
			bslTmp.AddGpsCovariance(&covTmp);
		}

		msr_ptr->AddGpsBaseline(((CDnaGpsBaseline*)&bslTmp));
	}
}
	
void dna_import::ParseDNAMSRGPSPoints(std::string& sBuf, dnaMsrPtr& msr_ptr, bool ignoreMsr)
{
	CDnaGpsPoint pntTmp;
	CDnaCovariance covTmp;

	pntTmp.SetReferenceFrame(msr_ptr->GetReferenceFrame());
	pntTmp.SetEpoch(msr_ptr->GetEpoch());

	// Measurement type
	std::string tmp;
	try {
		tmp = trimstr(sBuf.substr(dml_.msr_type, dmw_.msr_type));
	}
	catch (...) {
		SignalExceptionParseDNA("ParseDNAMSRGPSPoints(): Could not extract measurement type from the record:  ",
			sBuf, dml_.msr_type);
	}
	
	msr_ptr->SetType(tmp);
	pntTmp.SetType(tmp);
	pntTmp.SetIgnore(ignoreMsr);
	covTmp.SetType(tmp);
	covTmp.SetIgnore(ignoreMsr);

	// So no need to read database ID
	if (!projectSettings_.i.simulate_measurements)
	{
		// Capture msr_id and cluster_id (for database referencing)
		ParseDatabaseIds(sBuf, "ParseDNAMSRGPSPoints", msr_ptr->GetTypeC());
		msr_ptr->SetDatabaseMap(m_msr_db_map);
		pntTmp.SetDatabaseMap(m_msr_db_map);
		pntTmp.SetClusterID(msr_ptr->GetClusterID());
		covTmp.SetClusterID(msr_ptr->GetClusterID());
	}

	// Instrument station
	msr_ptr->SetFirst(ParseInstrumentValue(sBuf, "ParseDNAMSRGPSPoints"));
	pntTmp.SetFirst(msr_ptr->GetFirst());

	try {
		// Measurement type (i.e. LLH or XYZ)
		tmp = trimstr(sBuf.substr(dml_.msr_targ1, dmw_.msr_targ1));
	}
	catch (...) {
		SignalExceptionParseDNA("ParseDNAMSRGPSPoints(): Could not extract Y cluster coordinate type from the record:  ",
			sBuf, dml_.msr_targ1);
	}

	try {
		if (tmp.empty())
			msr_ptr->SetCoordType(XYZ_type);
		else
			msr_ptr->SetCoordType(tmp);
		pntTmp.SetCoordType(msr_ptr->GetCoordType());
	}
	catch (std::runtime_error& e) {
		SignalExceptionParseDNA("ParseDNAMSRGPSPoints(): " + std::string(e.what()),
			"", dml_.msr_targ1);
	}

	// Number of points
	UINT32 pntCount(1);
	msr_ptr->SetTotal(ParseMsrCountValue(sBuf, pntCount, "ParseDNAMSRGPSPoints"));
	msr_ptr->SetRecordedTotal(pntCount);
	pntTmp.SetRecordedTotal(pntCount);

	msr_ptr->GetPoints_ptr()->reserve(pntCount);
	g_parsemsr_tally.Y += pntCount * 3;

	// V-scale
	msr_ptr->SetVscale(ParseScaleVValue(sBuf, "ParseDNAMSRGPSPoints"));
	pntTmp.SetVscale(msr_ptr->GetVscale());

	// P-scale
	msr_ptr->SetPscale(ParseScalePValue(sBuf, "ParseDNAMSRGPSPoints"));
	pntTmp.SetPscale(msr_ptr->GetPscale());

	// L-scale supplied?
	msr_ptr->SetLscale(ParseScaleLValue(sBuf, "ParseDNAMSRGPSPoints"));
	pntTmp.SetLscale(msr_ptr->GetLscale());
	
	// H-scale supplied?
	msr_ptr->SetHscale(ParseScaleHValue(sBuf, "ParseDNAMSRGPSPoints"));
	pntTmp.SetHscale(msr_ptr->GetHscale());

	// Process for setting reference frame is:
	// 1. Set the default reference frame in ParseDNAMSR when msr_ptr is initialised
	// 2. Initialise bslTmp using msr_ptr
	// 3. When a reference frame element is found here, change the reference frame
	//    using checks on override
	// same as:
	//	- DnaXmlFormat_pskel::_start_element_impl when DnaMeasurement is reached
	//	- DnaMeasurement_pimpl::ReferenceFrame when ReferenceFrame tag is found
	//	- Clusterpoint_pimpl::pre()
	if (projectSettings_.i.override_input_rfame)
	{
		// Set the cluster frame
		msr_ptr->SetReferenceFrame(projectSettings_.i.reference_frame);
		// Set the point frame
		pntTmp.SetReferenceFrame(projectSettings_.i.reference_frame);
		pntTmp.SetEpoch(msr_ptr->GetEpoch());
	}
	else //if (!projectSettings_.i.override_input_rfame)
	{
		try {
			// Get the reference frame from the file
			tmp = ParseRefFrameValue(sBuf, "ParseDNAMSRGPSPoints");

			if (!tmp.empty())
			{
				// Set the cluster frame
				msr_ptr->SetReferenceFrame(tmp);
				// Set the point frame
				pntTmp.SetReferenceFrame(tmp);
			}
		}
		catch (std::runtime_error& e) {
			std::stringstream ss;
			ss << "ParseDNAMSRGPSPoints(): Error parsing reference frame:  " << std::endl << 
				"    " << e.what();
			SignalExceptionParseDNA(ss.str(), "", dml_.msr_gps_reframe);
		}

		try {
			// Get the epoch from the file
			tmp = ParseEpochValue(sBuf, "ParseDNAMSRGPSPoints");

			if (tmp.empty())
			{
				if (isEpsgWGS84Ensemble(epsgCodeFromName<UINT32, std::string>(msr_ptr->GetReferenceFrame())))
				{
					// Set the cluster epoch
					msr_ptr->SetEpoch("");
					// Set the baseline epoch
					pntTmp.SetEpoch("");
				}
			}
			else
			{
				// Set the cluster epoch
				msr_ptr->SetEpoch(tmp);
				// Set the point epoch
				pntTmp.SetEpoch(tmp);
			}
		}
		catch (std::runtime_error& e) {
			std::stringstream ss;
			ss << "ParseDNAMSRGPSPoints(): Error parsing epoch:  " << std::endl << 
				"    " << e.what();
			SignalExceptionParseDNA(ss.str(), "", dml_.msr_gps_epoch);
		}
	}
	
	bool first_run(true);
	UINT32 cov, covCount(pntCount - 1);

	for (UINT32 b(0); b<pntCount; ++b, --covCount)
	{
		if (!first_run)
		{
			m_lineNo++;
			// Obtain exclusive use of the input file pointer
			import_file_mutex.lock();
			getline((*ifsInputFILE_), sBuf);
			// release file pointer mutex
			import_file_mutex.unlock();

			// Instrument station
			msr_ptr->SetFirst(ParseInstrumentValue(sBuf, "ParseDNAMSRGPSPoints"));
			pntTmp.SetFirst(msr_ptr->GetFirst());

			// No need to read database ID for simulation
			if (!projectSettings_.i.simulate_measurements)
			{
				// Capture msr_id (for database referencing)
				ParseDatabaseIds(sBuf, "ParseDNAMSRGPSPoints", msr_ptr->GetTypeC());
				pntTmp.SetDatabaseMap(m_msr_db_map);
			}
		}

		first_run = false;

		if (!projectSettings_.i.simulate_measurements)
		{
			m_lineNo++;
			// Obtain exclusive use of the input file pointer
			import_file_mutex.lock();
			getline((*ifsInputFILE_), sBuf);
			// release file pointer mutex
			import_file_mutex.unlock();

			pntTmp.SetX(ParseGPSMsrValue(sBuf, "X", "ParseDNAMSRGPSPoints"));
			pntTmp.SetSigmaXX(ParseGPSVarValue(sBuf, "X", dml_.msr_gps_vcv_1, dmw_.msr_gps_vcv_1, "ParseDNAMSRGPSPoints"));

			m_lineNo++;
			// Obtain exclusive use of the input file pointer
			import_file_mutex.lock();
			getline((*ifsInputFILE_), sBuf);
			// release file pointer mutex
			import_file_mutex.unlock();

			pntTmp.SetY(ParseGPSMsrValue(sBuf, "Y", "ParseDNAMSRGPSPoints"));
			pntTmp.SetSigmaXY(ParseGPSVarValue(sBuf, "Y", dml_.msr_gps_vcv_1, dmw_.msr_gps_vcv_1, "ParseDNAMSRGPSPoints"));
			pntTmp.SetSigmaYY(ParseGPSVarValue(sBuf, "Y", dml_.msr_gps_vcv_2, dmw_.msr_gps_vcv_2, "ParseDNAMSRGPSPoints"));

			m_lineNo++;
			// Obtain exclusive use of the input file pointer
			import_file_mutex.lock();
			getline((*ifsInputFILE_), sBuf);
			// release file pointer mutex
			import_file_mutex.unlock();

			pntTmp.SetZ(ParseGPSMsrValue(sBuf, "Z", "ParseDNAMSRGPSPoints"));
			pntTmp.SetSigmaXZ(ParseGPSVarValue(sBuf, "Z", dml_.msr_gps_vcv_1, dmw_.msr_gps_vcv_1, "ParseDNAMSRGPSPoints"));
			pntTmp.SetSigmaYZ(ParseGPSVarValue(sBuf, "Z", dml_.msr_gps_vcv_2, dmw_.msr_gps_vcv_2, "ParseDNAMSRGPSPoints"));
			pntTmp.SetSigmaZZ(ParseGPSVarValue(sBuf, "Z", dml_.msr_gps_vcv_3, dmw_.msr_gps_vcv_3, "ParseDNAMSRGPSPoints"));
		}

		pntTmp.ResizeGpsCovariancesCount();
		pntTmp.ReserveGpsCovariancesCount(covCount);

		for (cov=0; cov<covCount; ++cov)
		{
			if (!projectSettings_.i.simulate_measurements)
				ParseDNAMSRCovariance(covTmp);
			pntTmp.AddPointCovariance(&covTmp);
		}

		msr_ptr->AddGpsPoint(((CDnaGpsPoint*)&pntTmp));
	}
}
	
void dna_import::ParseDNAMSRCovariance(CDnaCovariance& cov)
{
	std::string sBuf;
		
	m_lineNo++;
	// Obtain exclusive use of the input file pointer
	import_file_mutex.lock();
	getline((*ifsInputFILE_), sBuf);
	// release file pointer mutex
	import_file_mutex.unlock();

	// m11, m12, m13
	cov.SetM11(ParseGPSVarValue(sBuf, "co", dml_.msr_gps_vcv_1, dmw_.msr_gps_vcv_1, "ParseDNAMSRCovariance"));
	cov.SetM12(ParseGPSVarValue(sBuf, "co", dml_.msr_gps_vcv_2, dmw_.msr_gps_vcv_2, "ParseDNAMSRCovariance"));
	cov.SetM13(ParseGPSVarValue(sBuf, "co", dml_.msr_gps_vcv_3, dmw_.msr_gps_vcv_3, "ParseDNAMSRCovariance"));

	m_lineNo++;
	// Obtain exclusive use of the input file pointer
	import_file_mutex.lock();
	getline((*ifsInputFILE_), sBuf);
	// release file pointer mutex
	import_file_mutex.unlock();

	// m21, m22, m23
	cov.SetM21(ParseGPSVarValue(sBuf, "co", dml_.msr_gps_vcv_1, dmw_.msr_gps_vcv_1, "ParseDNAMSRCovariance"));
	cov.SetM22(ParseGPSVarValue(sBuf, "co", dml_.msr_gps_vcv_2, dmw_.msr_gps_vcv_2, "ParseDNAMSRCovariance"));
	cov.SetM23(ParseGPSVarValue(sBuf, "co", dml_.msr_gps_vcv_3, dmw_.msr_gps_vcv_3, "ParseDNAMSRCovariance"));

	m_lineNo++;
	// Obtain exclusive use of the input file pointer
	import_file_mutex.lock();
	getline((*ifsInputFILE_), sBuf);
	// release file pointer mutex
	import_file_mutex.unlock();

	// m31, m32, m33
	cov.SetM31(ParseGPSVarValue(sBuf, "co", dml_.msr_gps_vcv_1, dmw_.msr_gps_vcv_1, "ParseDNAMSRCovariance"));
	cov.SetM32(ParseGPSVarValue(sBuf, "co", dml_.msr_gps_vcv_2, dmw_.msr_gps_vcv_2, "ParseDNAMSRCovariance"));
	cov.SetM33(ParseGPSVarValue(sBuf, "co", dml_.msr_gps_vcv_3, dmw_.msr_gps_vcv_3, "ParseDNAMSRCovariance"));
}
	

void dna_import::ParseDatabaseIds(const std::string& sBuf, const std::string& calling_function, const char msrType)
{
	m_msr_db_map.is_msr_id_set = false;
	m_msr_db_map.is_cls_id_set = false;

	// No msr id?
	if (sBuf.length() <= dml_.msr_id_msr)
		return;

	// all measurement types: capture msr id 
	ParseDatabaseMsrId(sBuf, calling_function);

	// Initialise bms index.  This member will be set
	// in SerialiseBms() 
	//m_msr_db_map.bms_index = 0;

	// Initialise cluster id.  ParseDatabaseClusterId will 
	// update this member if a value is found
	m_msr_db_map.cluster_id = 0;
	
	switch (msrType)
	{
	// Cluster types: capture cluster ID
	case 'D':
	case 'G':
	case 'X':
	case 'Y':
		if (sBuf.length() > dml_.msr_id_cluster)
			ParseDatabaseClusterId(sBuf, calling_function);
	}
}
	

void dna_import::ParseDatabaseClusterId(const std::string& sBuf, const std::string& calling_function)
{
	std::string parsed_value;
	// Cluster ID
	try {
		parsed_value = trimstr(sBuf.substr(dml_.msr_id_cluster, dmw_.msr_id_cluster));
		if (!parsed_value.empty())
		{
			m_msr_db_map.cluster_id = val_uint<UINT32, std::string>(parsed_value);
			m_msr_db_map.is_cls_id_set = true;
		}
	}
	catch (std::runtime_error& e) {
		std::stringstream ss;
		ss << calling_function << "(): Could not extract database cluster id from the record:  " << 
		 	std::endl << "  " << e.what();
		SignalExceptionParseDNA(ss.str(), sBuf, dml_.msr_id_msr);
	}
}
	

void dna_import::ParseDatabaseMsrId(const std::string& sBuf, const std::string& calling_function)
{
	std::string parsed_value;
	// Measurement ID
	try {
		parsed_value = trimstr(sBuf.substr(dml_.msr_id_msr, dmw_.msr_id_msr));
		if (!parsed_value.empty())
		{
			m_msr_db_map.msr_id = val_uint<UINT32, std::string>(parsed_value);
			m_msr_db_map.is_msr_id_set = true;
		}
	}
	catch (std::runtime_error& e) {
		std::stringstream ss;
		ss << calling_function << "(): Could not extract database msr id from the record:  " << 
		 	std::endl << "  " << e.what();
		SignalExceptionParseDNA(ss.str(), sBuf, dml_.msr_id_msr);
	}
}
	

std::string dna_import::ParseAngularValue(const std::string& sBuf, const std::string& calling_function)
{
	std::string parsed_value, tmp;
	double d;
	UINT32 u;

	// degrees value
	try {
		parsed_value = trimstr(sBuf.substr(dml_.msr_ang_d, dmw_.msr_ang_d)) + ".";
	}
	catch (...) {
		SignalExceptionParseDNA(calling_function + "(): Could not extract degrees value from the record:  ",
			sBuf, dml_.msr_ang_d);
	}

	// minutes value
	try {
		tmp = trimstr(sBuf.substr(dml_.msr_ang_m, dmw_.msr_ang_m));
		u = LongFromString<UINT32>(tmp);
		if (u < 10)
			parsed_value.append("0");
		parsed_value.append(StringFromT<UINT32>(u));
	}
	catch (...) {
		SignalExceptionParseDNA(calling_function + "(): Could not extract minutes value from the record:  ",
			sBuf, dml_.msr_ang_m);
	}

	// seconds value
	size_t pos = 0;
	try {
		tmp = trimstr(sBuf.substr(dml_.msr_ang_s, dmw_.msr_ang_s));
		d = DoubleFromString<double>(tmp);
		if (d < 10 && tmp.at(0) != '0')
			parsed_value.append("0");
		if ((pos = tmp.find(".", pos)) != std::string::npos)
			parsed_value.append(tmp.substr(0, pos) + tmp.substr(pos+1));
		else
			parsed_value.append(tmp);
	}
	catch (...) {
		SignalExceptionParseDNA(calling_function + "(): Could not extract seconds value from the record:  ",
			sBuf, dml_.msr_ang_s);
	}
	return parsed_value;
}

std::string dna_import::ParseLinearValue(const std::string& sBuf, const std::string& msrName, const std::string& calling_function)
{
	try {
		return trimstr(sBuf.substr(dml_.msr_linear, dmw_.msr_linear));		// coordinate value
	}
	catch (...) {
		SignalExceptionParseDNA(calling_function + "(): Could not extract the " + msrName + " value from the record:  ",
			sBuf, dml_.msr_linear);
	}
	return "";
}

std::string dna_import::ParseInstrumentValue(const std::string& sBuf, const std::string& calling_function)
{

	try {
		// Capture string from the designated columns; throws on failure
		std::string stn = trimstr(sBuf.substr(dml_.msr_inst, dmw_.msr_inst));		// instrument station
	
		// No value supplied?
		if (stn.empty())
			SignalExceptionParseDNA(calling_function + "(): The following record doesn't contain an instrument station name:  ",
				sBuf, dml_.msr_inst);

		return stn;
	}
	catch (...) {
		SignalExceptionParseDNA(calling_function + "(): Could not extract first station name from the record:  ",
			sBuf, dml_.msr_inst);
	}
	return "";
}

std::string dna_import::ParseTargetValue(const std::string& sBuf, const std::string& calling_function)
{
	try {
		// Capture string from the designated columns; throws on failure
		std::string stn = trimstr(sBuf.substr(dml_.msr_targ1, dmw_.msr_targ1));		// first target station
	
		// No value supplied?
		if (stn.empty())
			SignalExceptionParseDNA(calling_function + "(): The following record doesn't contain a target station name:  ",
				sBuf, dml_.msr_targ1);

		return stn;		
	}
	catch (...) {
		SignalExceptionParseDNA(calling_function + "(): Could not extract second station name from the record:  ",
			sBuf, dml_.msr_targ1);
	}
	return "";
}

std::string dna_import::ParseTarget2Value(const std::string& sBuf, const std::string& calling_function)
{
	try {
		// Capture string from the designated columns; throws on failure
		std::string stn = trimstr(sBuf.substr(dml_.msr_targ2, dmw_.msr_targ2));		// second target station
	
		// No value supplied?
		if (stn.empty())
			SignalExceptionParseDNA(calling_function + "(): The following record doesn't contain a second target station name:  ",
				sBuf, dml_.msr_targ2);

		return stn;		
	}
	catch (...) {
		SignalExceptionParseDNA(calling_function + "(): Could not extract third station name from the record:  ",
			sBuf, dml_.msr_targ2);
	}
	return "";
}

std::string dna_import::ParseStdDevValue(const std::string& sBuf, const std::string& calling_function)
{
	std::string tmp;
	try {
		tmp = trimstr(sBuf.substr(dml_.msr_stddev, dmw_.msr_stddev));		// standard deviation
	}
	catch (...) {
		SignalExceptionParseDNA(calling_function + "(): Could not extract standard deviation from the record:  ",
			sBuf, dml_.msr_stddev);
	}

	if (DoubleFromString<double>(tmp) < PRECISION_1E25)
		SignalExceptionParseDNA(calling_function + "(): Invalid standard deviation (" + tmp + "). Values cannot be zero or negative:  ",
			sBuf, dml_.msr_stddev);
	return tmp;
}

std::string dna_import::ParseInstHeightValue(const std::string& sBuf, const std::string& calling_function)
{
	try {
		if (sBuf.length() > dml_.msr_targ_ht)
			return trimstr(sBuf.substr(dml_.msr_inst_ht, dmw_.msr_inst_ht));		// instrument height
		else
			return trimstr(sBuf.substr(dml_.msr_inst_ht));
	}
	catch (...) {
		SignalExceptionParseDNA(calling_function + "(): Could not extract instrument height from the record:  ",
			sBuf, dml_.msr_inst_ht);
	}
	return "";
}

std::string dna_import::ParseTargHeightValue(const std::string& sBuf, const std::string& calling_function)
{
	try {
		if (sBuf.length() > static_cast<std::string::size_type>(dml_.msr_targ_ht + 1 + dmw_.msr_targ_ht))
			return trimstr(sBuf.substr(dml_.msr_targ_ht, dmw_.msr_targ_ht));		// target height
		else
			return trimstr(sBuf.substr(dml_.msr_targ_ht));		
	}
	catch (...) {
		SignalExceptionParseDNA(calling_function + "(): Could not extract target height from the record:  ",
			sBuf, dml_.msr_targ_ht);
	}
	return "";
}

std::string dna_import::ParseMsrCountValue(const std::string& sBuf, UINT32& msrCount, const std::string& calling_function)
{
	try {
		std::string count(trimstr(sBuf.substr(dml_.msr_targ2, dmw_.msr_targ2)));		// number of measurements
		if (count.empty())
			SignalExceptionParseDNA(calling_function + "(): Could not extract number of measurements from the record:  ",
				sBuf, dml_.msr_targ2);
		msrCount = LongFromString<UINT32>(count);
		return count;
	}
	catch (...) {
		SignalExceptionParseDNA(calling_function + "(): Could not extract number of measurements from the record:  ",
			sBuf, dml_.msr_targ2);
	}
	return "";
}

std::string dna_import::ParseScaleVValue(const std::string& sBuf, const std::string& calling_function)
{
	if (sBuf.length() <= dml_.msr_gps_vscale)
		return "1";
	
	std::string scalar;
	try {
		if (sBuf.length() > dml_.msr_gps_pscale)
			scalar = trimstr(sBuf.substr(dml_.msr_gps_vscale, dmw_.msr_gps_vscale));		// v-scale
		else
			scalar = trimstr(sBuf.substr(dml_.msr_gps_vscale));					

		if (scalar.empty())
			return "1";
	}
	catch (...) {
		SignalExceptionParseDNA(calling_function + "(): Could not extract v-scale from the record:  ",
			sBuf, dml_.msr_gps_vscale);
	}

	return scalar;
}

std::string dna_import::ParseScalePValue(const std::string& sBuf, const std::string& calling_function)
{
	if (sBuf.length() <= dml_.msr_gps_pscale)
		return "1";
	
	std::string scalar;
	
	try {
		if (sBuf.length() > dml_.msr_gps_lscale)
			scalar = trimstr(sBuf.substr(dml_.msr_gps_pscale, dmw_.msr_gps_pscale));		// p-scale
		else
			scalar = trimstr(sBuf.substr(dml_.msr_gps_pscale));
	
		if (scalar.empty())
			return "1";
	}
	catch (...) {
		SignalExceptionParseDNA(calling_function + "(): Could not extract p-scale from the record:  ",
			sBuf, dml_.msr_gps_pscale);
	}

	return scalar;
}

std::string dna_import::ParseScaleLValue(const std::string& sBuf, const std::string& calling_function)
{
	if (sBuf.length() <= dml_.msr_gps_lscale)
		return "1";
	
	std::string scalar;

	try {
		if (sBuf.length() > dml_.msr_gps_hscale)
			scalar = trimstr(sBuf.substr(dml_.msr_gps_lscale, dmw_.msr_gps_lscale));		// l-scale
		else
			scalar = trimstr(sBuf.substr(dml_.msr_gps_lscale));					
	
		if (scalar.empty())
			return "1";
	}
	catch (...) {
		SignalExceptionParseDNA(calling_function + "(): Could not extract l-scale from the record:  ",
			sBuf, dml_.msr_gps_lscale);
	}

	return scalar;
}

std::string dna_import::ParseScaleHValue(const std::string& sBuf, const std::string& calling_function)
{
	if (sBuf.length() <= dml_.msr_gps_hscale)
		return "1";
	
	std::string scalar;

	try {
		if (sBuf.length() > (dml_.msr_gps_reframe))
			scalar = trimstr(sBuf.substr(dml_.msr_gps_hscale, dmw_.msr_gps_hscale));		// h-scale
		else
			scalar = trimstr(sBuf.substr(dml_.msr_gps_hscale));		
	
		if (scalar.empty())
			return "1";
	}
	catch (...) {
		SignalExceptionParseDNA(calling_function + "(): Could not extract h-scale from the record:  ",
			sBuf, dml_.msr_gps_hscale);
	}

	return scalar;
}

std::string dna_import::ParseRefFrameValue(const std::string& sBuf, const std::string& calling_function)
{
	// Override the reference frame with the user specified frame?
	if (projectSettings_.i.override_input_rfame)
		return projectSettings_.i.reference_frame;
	
	if (sBuf.length() <= dml_.msr_gps_reframe)
		return "";
	
	std::string frame;
	try {
		if (sBuf.length() > (dml_.msr_gps_epoch))
			frame = trimstr(sBuf.substr(dml_.msr_gps_reframe, dmw_.msr_gps_reframe));		// reference frame
		else
			frame = trimstr(sBuf.substr(dml_.msr_gps_reframe));		
	}
	catch (...) {
		SignalExceptionParseDNA(calling_function + "(): Could not extract reference frame from the record:  ",
			sBuf, dml_.msr_gps_reframe);
	}

	return frame;
}

std::string dna_import::ParseEpochValue(const std::string& sBuf, const std::string& calling_function)
{
	if (sBuf.length() <= dml_.msr_gps_epoch)
		return "";
	
	std::string epoch;
	try {
		if (sBuf.length() > static_cast<std::string::size_type>(dml_.msr_gps_epoch + dmw_.msr_gps_epoch))
			epoch = trimstr(sBuf.substr(dml_.msr_gps_epoch, dmw_.msr_gps_epoch));		// epoch
		else
			epoch = trimstr(sBuf.substr(dml_.msr_gps_epoch));
	}
	catch (...) {
		SignalExceptionParseDNA(calling_function + "(): Could not extract epoch from the record:  ",
			sBuf, dml_.msr_gps_epoch);
	}

	return epoch;
}

std::string dna_import::ParseGPSMsrValue(const std::string& sBuf, const std::string& element, const std::string& calling_function)
{
	try {
		return trimstr(sBuf.substr(dml_.msr_gps, dmw_.msr_gps));				// value
	}
	catch (...) {
		SignalExceptionParseDNA(calling_function + "(): Could not extract GNSS " + element + " measurement from the record:  ",
			sBuf, dml_.msr_gps);
	}
	return "";
}

std::string dna_import::ParseGPSVarValue(const std::string& sBuf, const std::string& element, const UINT32 location, const UINT32 width, const std::string& calling_function)
{
	try {
		return trimstr(sBuf.substr(location, width));		// variance
	}
	catch (...) {
		SignalExceptionParseDNA(calling_function + "(): Could not extract GNSS " + element + " variance from the record:  ",
			sBuf, location);
	}
	return "";
}

void dna_import::ParseDNAMSRAngular(const std::string& sBuf, dnaMsrPtr& msr_ptr)
{
	// Measurement type
	try {
		msr_ptr->SetType(trimstr(sBuf.substr(dml_.msr_type, dmw_.msr_type)));
	}
	catch (...) {
		SignalExceptionParseDNA("ParseDNAMSRAngular(): Could not extract measurement type from the record:  ",
			sBuf, dml_.msr_type);
	}
	
	// Instrument station
	msr_ptr->SetFirst(ParseInstrumentValue(sBuf, "ParseDNAMSRAngular"));
	
	// Target station
	if (msr_ptr->GetMsrStnCount() >= TWO_STATION)
		msr_ptr->SetTarget(ParseTargetValue(sBuf, "ParseDNAMSRAngular"));
	
	// Second target station
	if (msr_ptr->GetMsrStnCount() == THREE_STATION)
		msr_ptr->SetTarget2(ParseTarget2Value(sBuf, "ParseDNAMSRAngular"));
	
	// If import is being used to simulate measurements, return as the
	// expected input is:
	//	 Measurement
	//	 From station
	//	 To station
	// So no need to read any further
	if (projectSettings_.i.simulate_measurements)
		return;

	// Value
	msr_ptr->SetValue(ParseAngularValue(sBuf, "ParseDNAMSRAngular"));
	
	// Standard deviation
	msr_ptr->SetStdDev(ParseStdDevValue(sBuf, "ParseDNAMSRAngular"));

	// Epoch
	msr_ptr->SetEpoch(ParseEpochValue(sBuf, "ParseDNAMSRAngular"));

	// Capture msr_id and cluster_id (for database referencing), then set
	// database id info
	ParseDatabaseIds(sBuf, "ParseDNAMSRAngular", msr_ptr->GetTypeC());
	msr_ptr->SetDatabaseMap(m_msr_db_map);

	// instrument and target heights only make sense for 
	// slope distances, vertical angles and zenith distances
	switch (msr_ptr->GetTypeC())
	{
	case 'V':
	case 'Z':
		break;
	default:
		// don't need instrument height for all other measurement types
		return;
	}

	// Instrument and target heights
	if (sBuf.length() <= dml_.msr_inst_ht)
	{
		msr_ptr->SetInstrumentHeight("0.0");
		msr_ptr->SetTargetHeight("0.0");
		return;
	}
	msr_ptr->SetInstrumentHeight(ParseInstHeightValue(sBuf, "ParseDNAMSRAngular"));
	
	// Target height
	if (sBuf.length() <= dml_.msr_targ_ht)
	{
		msr_ptr->SetTargetHeight("0.0");
		return;
	}	
	msr_ptr->SetTargetHeight(ParseTargHeightValue(sBuf, "ParseDNAMSRAngular"));
}
	

UINT32 dna_import::ParseDNAMSRDirections(std::string& sBuf, dnaMsrPtr& msr_ptr, bool ignoreMsr)
{
	// Measurement type
	try {
		msr_ptr->SetType(trimstr(sBuf.substr(dml_.msr_type, dmw_.msr_type)));
	}
	catch (...) {
		SignalExceptionParseDNA("ParseDNAMSRDirections(): Could not extract measurement type from the record:  ",
			sBuf, dml_.msr_type);
	}

	// Instrument station
	msr_ptr->SetFirst(ParseInstrumentValue(sBuf, "ParseDNAMSRDirections"));
	
	// Target station
	if (msr_ptr->GetMsrStnCount() >= TWO_STATION)
		msr_ptr->SetTarget(ParseTargetValue(sBuf, "ParseDNAMSRDirections"));
	
	// Number of directions
	UINT32 dirnCount, dirnCountLessIgnored;
	msr_ptr->SetTotal(ParseMsrCountValue(sBuf, dirnCount, "ParseDNAMSRDirections"));
	msr_ptr->GetDirections_ptr()->reserve(dirnCount);
	
	CDnaDirection dirnTmp;
	dirnTmp.SetRecordedTotal(0);
	
	// If import is being used to simulate measurements, return as the
	// expected input is:
	//	 Measurement
	//	 From station
	//	 To station
	// So no need to read direction values or database ID
	if (!projectSettings_.i.simulate_measurements)
	{
		// Value
		msr_ptr->SetValue(ParseAngularValue(sBuf, "ParseDNAMSRDirections"));

		// Standard deviation
		msr_ptr->SetStdDev(ParseStdDevValue(sBuf, "ParseDNAMSRDirections"));

		// Epoch
		msr_ptr->SetEpoch(ParseEpochValue(sBuf, "ParseDNAMSRDirections"));

		// Capture msr_id and cluster_id (for database referencing)
		ParseDatabaseIds(sBuf, "ParseDNAMSRDirections", msr_ptr->GetTypeC());
		msr_ptr->SetDatabaseMap(m_msr_db_map);
		dirnTmp.SetClusterID(msr_ptr->GetClusterID());
	}

	bool subignoreMsr;
	dirnCountLessIgnored = dirnCount;

	if (dirnCount == 0)
		throw XMLInteropException("Zero direction count.", m_lineNo);

	for (UINT32 dirn=0; dirn<dirnCount; ++dirn)
	{
		m_lineNo++;
		// Obtain exclusive use of the input file pointer
		import_file_mutex.lock();
		getline((*ifsInputFILE_), sBuf);
		// release file pointer mutex
		import_file_mutex.unlock();

		// get ignore flag for sub direction and remove accordingly
		subignoreMsr = iequals("*", sBuf.substr(dml_.msr_ignore, dmw_.msr_ignore));

		if (subignoreMsr)
		{
			dirnCountLessIgnored--;
			if (dirnCountLessIgnored == 0)
			{
				// throw if the entire direction set is not ignored
				if (!ignoreMsr)
					throw XMLInteropException("There aren't any non-ignored directions in the set.", m_lineNo);
			}
		}

		dirnTmp.SetFirst(msr_ptr->GetFirst());
		dirnTmp.SetIgnore(subignoreMsr);
		g_parsemsr_tally.D++;

		// Second target station
		dirnTmp.SetTarget(ParseTarget2Value(sBuf, "ParseDNAMSRDirections"));

		// Again, no need to read direction values
		if (!projectSettings_.i.simulate_measurements)
		{
			// Value
			dirnTmp.SetValue(ParseAngularValue(sBuf, "ParseDNAMSRDirections"));

			// Standard deviation
			dirnTmp.SetStdDev(ParseStdDevValue(sBuf, "ParseDNAMSRDirections"));

			// Capture msr_id and cluster_id (for database referencing)
			ParseDatabaseIds(sBuf, "ParseDNAMSRDirections", msr_ptr->GetTypeC());
			dirnTmp.SetDatabaseMap(m_msr_db_map);
		}

		// add the direction
		msr_ptr->AddDirection(((CDnaMeasurement*)&dirnTmp));
	}

	msr_ptr->SetNonIgnoredDirns(dirnCountLessIgnored);

	return dirnCount;
}

void dna_import::RemoveIgnoredMeasurements(vdnaMsrPtr* vMeasurements, MsrTally* parsemsrTally)
{
	// Remove ignored measurements
	erase_if(vMeasurements, CompareIgnoredMsr<CDnaMeasurement>());
	// Rebuild measurement tally
	parsemsrTally->CreateTally(*vMeasurements);
}
	

void dna_import::IncludeMeasurementTypes(const std::string& includeMsrs, vdnaMsrPtr* vMeasurements, MsrTally* parsemsrTally)
{
	_it_vdnamsrptr _it_msr;
	std::string msrTypes = includeMsrs;
	// convert to upper case
	str_toupper<int>(msrTypes);

	CompareNonMeasType<CDnaMeasurement> meastypeCompareFunc(msrTypes);
	erase_if(vMeasurements, meastypeCompareFunc);

	MsrTally msrtallyTmp;

	for (UINT16 i=0; i<msrTypes.length(); ++i)
	{
		switch ((msrTypes.c_str())[i])
		{
		case 'A': // Horizontal angle
			msrtallyTmp.A = parsemsrTally->A;
			break;
		case 'B': // Geodetic azimuth
			msrtallyTmp.B = parsemsrTally->B;
			break;
		case 'C': // Chord dist
			msrtallyTmp.C = parsemsrTally->C;
			break;
		case 'D': // Direction set
			msrtallyTmp.D = parsemsrTally->D;
			break;
		case 'E': // Ellipsoid arc
			msrtallyTmp.E = parsemsrTally->E;
			break;
		case 'G': // GPS Baseline
			msrtallyTmp.G = parsemsrTally->G;
			break;
		case 'H': // Orthometric height
			msrtallyTmp.H = parsemsrTally->H;
			break;
		case 'I': // Astronomic latitude
			msrtallyTmp.I = parsemsrTally->I;
			break;
		case 'J': // Astronomic longitude
			msrtallyTmp.J = parsemsrTally->J;
			break;
		case 'K': // Astronomic azimuth
			msrtallyTmp.K = parsemsrTally->K;
			break;
		case 'L': // Level difference
			msrtallyTmp.L = parsemsrTally->L;
			break;
		case 'M': // MSL arc
			msrtallyTmp.M = parsemsrTally->M;
			break;
		case 'P': // Geodetic latitude
			msrtallyTmp.P = parsemsrTally->P;
			break;
		case 'Q': // Geodetic longitude
			msrtallyTmp.Q = parsemsrTally->Q;
			break;
		case 'R': // Ellipsoidal height
			msrtallyTmp.R = parsemsrTally->R;
			break;
		case 'S': // Slope distance
			msrtallyTmp.S = parsemsrTally->S;
			break;
		case 'V': // Zenith distance
			msrtallyTmp.V = parsemsrTally->V;
			break;
		case 'X': // GPS Baseline cluster
			msrtallyTmp.X = parsemsrTally->X;
			break;
		case 'Y': // GPS point cluster
			msrtallyTmp.Y = parsemsrTally->Y;
			break;
		case 'Z': // Vertical angle
			msrtallyTmp.Z = parsemsrTally->Z;
			break;
		}
	}

	*parsemsrTally = msrtallyTmp;
}

void dna_import::ExcludeMeasurementTypes(const std::string& excludeMsrs, vdnaMsrPtr* vMeasurements, MsrTally* parsemsrTally)
{
	_it_vdnamsrptr _it_msr;
	std::string msrTypes = excludeMsrs;
	// convert to upper case
	str_toupper<int>(msrTypes);

	CompareMeasType<CDnaMeasurement> meastypeCompareFunc(msrTypes);
	erase_if(vMeasurements, meastypeCompareFunc);

	for (UINT16 i=0; i<msrTypes.length(); ++i)
	{
		switch ((msrTypes.c_str())[i])
		{
		case 'A': // Horizontal angle
			parsemsrTally->A = 0;
			break;
		case 'B': // Geodetic azimuth
			parsemsrTally->B = 0;
			break;
		case 'C': // Chord dist
			parsemsrTally->C = 0;
			break;
		case 'D': // Direction set
			parsemsrTally->D = 0;
			break;
		case 'E': // Ellipsoid arc
			parsemsrTally->E = 0;
			break;
		case 'G': // GPS Baseline
			parsemsrTally->G = 0;
			break;
		case 'H': // Orthometric height
			parsemsrTally->H = 0;
			break;
		case 'I': // Astronomic latitude
			parsemsrTally->I = 0;
			break;
		case 'J': // Astronomic longitude
			parsemsrTally->J = 0;
			break;
		case 'K': // Astronomic azimuth
			parsemsrTally->K = 0;
			break;
		case 'L': // Level difference
			parsemsrTally->L = 0;
			break;
		case 'M': // MSL arc
			parsemsrTally->M = 0;
			break;
		case 'P': // Geodetic latitude
			parsemsrTally->P = 0;
			break;
		case 'Q': // Geodetic longitude
			parsemsrTally->Q = 0;
			break;
		case 'R': // Ellipsoidal height
			parsemsrTally->R = 0;
			break;
		case 'S': // Slope distance
			parsemsrTally->S = 0;
			break;
		case 'V': // Zenith distance
			parsemsrTally->V = 0;
			break;
		case 'X': // GPS Baseline cluster
			parsemsrTally->X = 0;
			break;
		case 'Y': // GPS point cluster
			parsemsrTally->Y = 0;
			break;
		case 'Z': // Vertical angle
			parsemsrTally->Z = 0;
			break;
		}
	}
}
	

void dna_import::ExcludeAllOutsideBoundingBox(vdnaStnPtr* vStations, vdnaMsrPtr* vMeasurements, 
		StnTally* parsestnTally, MsrTally* parsemsrTally, pvstring vExcludedStns, const project_settings& p,
		bool& splitXmsrs, bool& splitYmsrs)
{
	splitXmsrs = splitYmsrs = false;

	// backup station vector
	vdnaStnPtr bvStation = *vStations;

	// TestNotEqualStationName requires vStations to be sorted (for binary_search)
	std::sort(vStations->begin(), vStations->end());

	// Remove all stations outside the bounding box
	// vExcludedStns will contain the names of all the stations (if any) that were outside the bounding box.
	// This list will be used to strip the corresponding measurements.
	vExcludedStns->clear();
	FindStnsWithinBoundingBox<double> boundingBoxFunc(bbox_upperLat_, bbox_upperLon_, bbox_lowerLat_, bbox_lowerLon_, vExcludedStns);
	erase_if(vStations, boundingBoxFunc);

	// create vstring of used stations
	vstring vIncludedStns;
	vIncludedStns.reserve(vStations->size());
	_it_vdnastnptr _it_stn(vStations->begin());
	for_each(vStations->begin(), vStations->end(), 
		[&vIncludedStns] (dnaStnPtr& s) { 
			vIncludedStns.push_back(s.get()->GetName());
	});
	
	// FindMsrsConnectedToStns requires vExcludedStns to be sorted (for binary_search)
	std::sort(vExcludedStns->begin(), vExcludedStns->end());

	// OK, now measurements...
	if (p.i.include_transcending_msrs == 1)
		ExtractAssociatedMsrsConnectedToStns(vMeasurements, parsemsrTally, &vIncludedStns, vExcludedStns, p, splitXmsrs, splitYmsrs);
	else
		ExtractAssociatedMsrsBoundingBox(vMeasurements, parsemsrTally, &vIncludedStns, vExcludedStns, p, splitXmsrs, splitYmsrs);

	strip_duplicates(vIncludedStns);

	// restore station vector, and erase stations not in vIncludedStns
	*vStations = bvStation;
	vExcludedStns->clear();
	TestNotEqualStationName<dnaStnPtr, std::string> selectStnFunc(&vIncludedStns, vExcludedStns);
	erase_if(vStations, selectStnFunc);

	// Rebuild station tally
	parsestnTally->CreateTally(*vStations);
}
	

void dna_import::ExtractStnsAndAssociatedMsrs(const std::string& stnListInclude, const std::string& stnListExclude, vdnaStnPtr* vStations, vdnaMsrPtr* vMeasurements, 
		StnTally* parsestnTally, MsrTally* parsemsrTally, pvstring vExcludedStns, const project_settings& p,
		bool& splitXmsrs, bool& splitYmsrs)
{
	splitXmsrs = splitYmsrs = false;

	// backup station vector
	vdnaStnPtr bvStations = *vStations;

	const std::string *stnListIn;

	vstring vIncludedStns;
	pvstring pvStnsIn, pvStnsEx;
	
	// Has the user provided a list of stations to include?
	if (!stnListInclude.empty())
	{
		pvStnsIn = &vIncludedStns;
		pvStnsEx = vExcludedStns;		
		stnListIn = &stnListInclude; 
	}
	// Has the user provided a list of stations to exclude?
	else if (!stnListExclude.empty()) 
	{
		pvStnsIn = vExcludedStns;
		pvStnsEx = &vIncludedStns;
		stnListIn = &stnListExclude; 
	}
	else
	{
		// nothing to do
		return;
	}

	// Form vector of included stations from user-input string
	BuildExtractStationsList(*stnListIn, pvStnsIn);
	
	// TestNotEqualStationName requires stations to be sorted (for binary_search)
	strip_duplicates(pvStnsIn);		// Strip duplicates and sort
	std::sort(pvStnsIn->begin(), pvStnsIn->end());
		
	pvStnsEx->clear();		

	// Remove all stations not in vIncludedStns
	// vExcludedStns will contain the names of all the stations (if any) that were outside the bounding box.
	// This list will be used to strip the corresponding measurements.
	if (!stnListInclude.empty())
	{
		TestNotEqualStationName<dnaStnPtr, std::string> selectStnFunc(pvStnsIn, pvStnsEx);	
		erase_if(vStations, selectStnFunc);
	}
	else
	{
		TestEqualStationName<dnaStnPtr, std::string> selectStnFunc(pvStnsIn, pvStnsEx);	
		erase_if(vStations, selectStnFunc);
		pvStnsIn = &vIncludedStns;
		pvStnsEx = vExcludedStns;
		std::sort(pvStnsIn->begin(), pvStnsIn->end());
	}

	// FindMsrsConnectedToStns requires vExcludedStns to be sorted (for binary_search)
	std::sort(pvStnsEx->begin(), pvStnsEx->end());

	// OK, get all measurements connected associated with vIncludedStns, splitting clusters as necessary
	if (!stnListInclude.empty())
		ExtractAssociatedMsrsConnectedToStns(vMeasurements, parsemsrTally, pvStnsIn, pvStnsEx, p, splitXmsrs, splitYmsrs);
	else
		ExtractAssociatedMsrsBoundingBox(vMeasurements, parsemsrTally, pvStnsIn, pvStnsEx, p, splitXmsrs, splitYmsrs);

	strip_duplicates(pvStnsIn);

	// restore station vector, and erase stations not in pvStnsIn
	// This repeat step is required since pvStnsIn may contain more 
	// or less depending on how many stations the function
	// ExtractAssociatedMsrsConnectedToStns retrieves
	*vStations = bvStations;
	pvStnsEx->clear();
	TestNotEqualStationName<dnaStnPtr, std::string> selectStnFunc(pvStnsIn, pvStnsEx);	
	erase_if(vStations, selectStnFunc);

	// Rebuild station tally
	parsestnTally->CreateTally(*vStations);	
}


void dna_import::SplitClusterMsrs(vdnaMsrPtr& msrsConnectedToStns, 
	pvstring pvIncludedStns, pvstring pvExcludedStns, 
	vdnaMsrPtr* vMeasurements, bool& splitXmsrs, bool& splitYmsrs)
{
	// Backup (all) measurements straddling selection
	vdnaMsrPtr XYmsrsConnectedToStns(msrsConnectedToStns);
		
	// Get non-XY measurements
	CompareMeasType<CDnaMeasurement> meastypeCompareFunc_XY("XY");
	erase_if(msrsConnectedToStns, meastypeCompareFunc_XY);
		
	// Get XY measurements
	CompareNonMeasType<CDnaMeasurement> meastypeCompareFunc_nonXY("XY");
	erase_if(XYmsrsConnectedToStns, meastypeCompareFunc_nonXY);

	// Split the XY measurements, keeping only those measurements directly
	// connected to pvIncludedStns
	SplitClusterMsrsConnectedToStns(&XYmsrsConnectedToStns, pvIncludedStns, pvExcludedStns, splitXmsrs, splitYmsrs);

	// Add the split XY measurements
	vMeasurements->insert(vMeasurements->end(), 
		XYmsrsConnectedToStns.begin(), 
		XYmsrsConnectedToStns.end());

	// Add the non-XY measurements
	vMeasurements->insert(vMeasurements->end(), 
		msrsConnectedToStns.begin(), 
		msrsConnectedToStns.end());
}
	
void dna_import::SplitClusterMsrsConnectedToStns(vdnaMsrPtr* vClusterMsrs, pvstring pvIncludedStns, pvstring pvExcludedStns, bool& splitXmsrs, bool& splitYmsrs)
{
	_it_vdnamsrptr _it_msr;
	std::vector<CDnaGpsBaseline>* vgpsBsls;
	std::vector<CDnaGpsPoint>* vgpsPnts;
	std::vector<CDnaCovariance>* vgpsCovs;
	std::vector<CDnaGpsBaseline>::iterator _it_gps_bsl;
	std::vector<CDnaGpsPoint>::iterator _it_gps_pnt;
	std::vector<CDnaCovariance>::iterator _it_gps_cov;

	UINT32 msr(0), index(0), subindex(0), keepCount;
	vUINT32 vIndices;

	for (_it_msr=vClusterMsrs->begin(); _it_msr!=vClusterMsrs->end(); ++_it_msr, ++msr)
	{
		switch (_it_msr->get()->GetTypeC())
		{
		case 'X':
			vgpsBsls = _it_msr->get()->GetBaselines_ptr();
			vIndices.clear();
			index = 0;

			// Step 1. Set ignore flag for all baselines to true
			for_each(vgpsBsls->begin(), vgpsBsls->end(), 
				[this](CDnaGpsBaseline& bsl) {
					bsl.SetIgnore(true);
			});

			// Step 2. Set the baselines connected to the included stations to false
			for_each(vgpsBsls->begin(), vgpsBsls->end(), 
				[pvIncludedStns](CDnaGpsBaseline& bsl) {
					if (binary_search(pvIncludedStns->begin(), pvIncludedStns->end(), bsl.GetFirst()))
						bsl.SetIgnore(false);
					else if (binary_search(pvIncludedStns->begin(), pvIncludedStns->end(), bsl.GetTarget()))
						bsl.SetIgnore(false);
			});

			// Step 3. Create a vector of measurement indices which are to be removed
			for_each(vgpsBsls->begin(), vgpsBsls->end(), 
				[&vIndices, &index](CDnaGpsBaseline& bsl) {
					if (bsl.GetIgnore())
						vIndices.push_back(index);
					index++;
			});

			// Any baselines to remove?
			if (vIndices.empty())
				continue;

			splitXmsrs = true;

			// Step 4. Set the ignore flag for the unwanted covariances using the vector of
			// indices created above
			index = 0;
			for (_it_gps_bsl=vgpsBsls->begin(); _it_gps_bsl!=vgpsBsls->end(); ++_it_gps_bsl, ++index)
			{
				// is this cluster measurement ignored?  If so, continue as the entire record will be
				// removed later
				if (binary_search(vIndices.begin(), vIndices.end(), index))
					continue;
				subindex = index+1;
				vgpsCovs = _it_gps_bsl->GetCovariances_ptr();
				for (_it_gps_cov=vgpsCovs->begin(); _it_gps_cov!=vgpsCovs->end(); ++_it_gps_cov, ++subindex)
					if (binary_search(vIndices.begin(), vIndices.end(), subindex))
						_it_gps_cov->SetIgnore(true);
			}

			// Step 5. Remove the 'ignored' cluster measurements, and reset the measurement count
			erase_if(vgpsBsls, CompareIgnoreedClusterMeas<CDnaGpsBaseline>());
			
			// Step 6. Set the total count for the cluster and for each clustered measurement
			keepCount = static_cast<UINT32>(vgpsBsls->size());
			_it_msr->get()->SetTotal(keepCount);
			
			if (vgpsBsls->empty())
				continue;

			for_each(vgpsBsls->begin(), vgpsBsls->end(), 
				[&keepCount](CDnaGpsBaseline& bsl){
					bsl.SetTotal(keepCount);
			});

			// Step 7. Remove the 'ignored' cluster covariances, and reset the measurement count
			for (_it_gps_bsl=vgpsBsls->begin(); _it_gps_bsl!=vgpsBsls->end(); ++_it_gps_bsl, ++index)
			{
				vgpsCovs = _it_gps_bsl->GetCovariances_ptr();
				erase_if(vgpsCovs, CompareIgnoreedClusterMeas<CDnaCovariance>());
			}

			break;
		case 'Y':
			vgpsPnts = _it_msr->get()->GetPoints_ptr();
			vIndices.clear();
			index = 0;

			// Step 1. Create a vector of measurement indices which are to be removed, and set
			// the 'ignore' flag for the unwanted GPS point clustered measurements
			for_each(vgpsPnts->begin(), vgpsPnts->end(), 
				[&vIndices, &index, pvExcludedStns](CDnaGpsPoint& pnt) {
					if (binary_search(pvExcludedStns->begin(), pvExcludedStns->end(), pnt.GetFirst()))
					{
						pnt.SetIgnore(true);
						vIndices.push_back(index);
					}
					index++;
			});

			// Any points to remove?
			if (vIndices.empty())
				continue;
			
			splitYmsrs = true;

			// Step 2. Set the ignore flag for the unwanted covariances using the vector of
			// indices created above
			index = 0;
			for (_it_gps_pnt=vgpsPnts->begin(); _it_gps_pnt!=vgpsPnts->end(); ++_it_gps_pnt, ++index)
			{
				// is this cluster measurement ignored?  If so, continue as the entire record will be
				// removed later
				if (binary_search(vIndices.begin(), vIndices.end(), index))
					continue;
				subindex = index+1;
				vgpsCovs = _it_gps_pnt->GetCovariances_ptr();
				for (_it_gps_cov=vgpsCovs->begin(); _it_gps_cov!=vgpsCovs->end(); ++_it_gps_cov, ++subindex)
					if (binary_search(vIndices.begin(), vIndices.end(), subindex))
						_it_gps_cov->SetIgnore(true);
			}			

			// Step 3. Remove the 'ignored' cluster measurements, and reset the measurement count
			erase_if(vgpsPnts, CompareIgnoreedClusterMeas<CDnaGpsPoint>());
			
			// Step 4. Set the total count for the cluster and for each clustered measurement
			keepCount = static_cast<UINT32>(vgpsPnts->size());
			_it_msr->get()->SetTotal(keepCount);

			if (vgpsPnts->empty())
				continue;
			
			for_each(vgpsPnts->begin(), vgpsPnts->end(), 
				[&keepCount](CDnaGpsPoint& pnt){
					pnt.SetTotal(keepCount);
			});

			// Step 5. Remove the 'ignored' cluster covariances, and reset the measurement count
			for (_it_gps_pnt=vgpsPnts->begin(); _it_gps_pnt!=vgpsPnts->end(); ++_it_gps_pnt, ++index)
			{
				vgpsCovs = _it_gps_pnt->GetCovariances_ptr();
				erase_if(vgpsCovs, CompareIgnoreedClusterMeas<CDnaCovariance>());
			}

			break;
		}
	}

	// Step 6. Erase cluster with no measurements
	erase_if(vClusterMsrs, CompareEmptyClusterMeas<CDnaMeasurement>());
}
	
void dna_import::ExtractAssociatedMsrsConnectedToStns(vdnaMsrPtr* vMeasurements, MsrTally* parsemsrTally,
		pvstring pvIncludedStns, pvstring pvExcludedStns, const project_settings& p,
		bool& splitXmsrs, bool& splitYmsrs)
{
	// Create functor to pull measurements connected to used stations
	// FindMsrsConnectedToStns will find measurements that has any one station
	// in pvIncludedStns for first, second or third station
	FindMsrsConnectedToStns<pvstring> msrstoIncludedStnsFunc(pvIncludedStns);

	vdnaMsrPtr msrsConnectedToStns;

	// Copy all measurements conected to selected stations
	// This will pull all measurements tied directly or indirectly to stations in 
	// pvIncludedStns.  If a measurement tied to a used station is part of a 
	// cluster, then the whole cluster will be pulled in
	copy_if(vMeasurements->begin(), vMeasurements->end(), 
		back_inserter(msrsConnectedToStns), msrstoIncludedStnsFunc);

	// clear the original measurement list
	vMeasurements->clear();	
	
	// Split clusters?
	if (p.i.split_clusters == 1)
		// Split the GNSS clusters so that measurements not connected
		// to the included stations are removed
		SplitClusterMsrs(msrsConnectedToStns, pvIncludedStns, pvExcludedStns, 
			vMeasurements, splitXmsrs, splitYmsrs);
	else
		// Add the measurements wholly connected to the included stations
		vMeasurements->insert(vMeasurements->end(), 
			msrsConnectedToStns.begin(), msrsConnectedToStns.end());
	
	// Get the additional stations introduced by connected measurements
	ExtractAssociatedStns(vMeasurements, pvIncludedStns);
	
	// Rebuild measurement tally
	parsemsrTally->CreateTally(*vMeasurements);
}
	

void dna_import::ExtractAssociatedMsrsBoundingBox(vdnaMsrPtr* vMeasurements, MsrTally* parsemsrTally,
		pvstring pvIncludedStns, pvstring pvExcludedStns, const project_settings& p,
		bool& splitXmsrs, bool& splitYmsrs)
{

	if (p.i.split_clusters == 1)
	{
		// Create functor to pull measurements connected to used stations
		// FindMsrsConnectedToStns will find measurements that has any one station
		// in pvIncludedStns for first, second or third station
		FindMsrsConnectedToStns<pvstring> msrstoIncludedStnsFunc(pvIncludedStns);

		vdnaMsrPtr msrsConnectedToStns;

		// Copy all measurements conected to selected stations
		// This will pull all measurements tied directly or indirectly to stations in 
		// pvIncludedStns.  If a measurement tied to a used station is part of a 
		// cluster, then the whole cluster will be pulled in
		copy_if(vMeasurements->begin(), vMeasurements->end(), 
			back_inserter(msrsConnectedToStns), msrstoIncludedStnsFunc);

		// clear the original measurement list
		vMeasurements->clear();	

		// Split the GNSS clusters so that measurements not connected
		// to the included stations are removed
		SplitClusterMsrs(msrsConnectedToStns, pvIncludedStns, pvExcludedStns, 
			vMeasurements, splitXmsrs, splitYmsrs);
	}
	else
	{
		// If a measurement is connected to a station that is outside the box or not on the 
		// selected station list, then delete it
		FindMsrsConnectedToStns<pvstring> msrstoExcludedStnsFunc(pvExcludedStns);

		// Delete all measurements not wholly within the box
		erase_if(vMeasurements, msrstoExcludedStnsFunc);
	}

	// Get the reduced set of stations connected to the measurements
	ExtractAssociatedStns(vMeasurements, pvIncludedStns);

	// Rebuild measurement tally
	parsemsrTally->CreateTally(*vMeasurements);
}


void dna_import::ExtractAssociatedStns_GX(std::vector<CDnaGpsBaseline>* vGpsBaselines, pvstring pvUsedStns)
{
	std::string station;
	std::vector<CDnaGpsBaseline>::iterator _it_msr(vGpsBaselines->begin());
	_it_pair_vstring it_msr_stns;

	for (_it_msr=vGpsBaselines->begin(); _it_msr!=vGpsBaselines->end(); _it_msr++)
	{
		station = _it_msr->GetFirst();
		it_msr_stns = equal_range(pvUsedStns->begin(), pvUsedStns->end(), station);
		if (it_msr_stns.first == it_msr_stns.second)
		{
			pvUsedStns->push_back(station);
			std::sort(pvUsedStns->begin(), pvUsedStns->end());
		}

		station = _it_msr->GetTarget();
		it_msr_stns = equal_range(pvUsedStns->end(), pvUsedStns->end(), station);
		if (it_msr_stns.first == it_msr_stns.second)
		{
			pvUsedStns->push_back(station);
			std::sort(pvUsedStns->begin(), pvUsedStns->end());
		}
	}
}

void dna_import::ExtractAssociatedStns_Y(std::vector<CDnaGpsPoint>* vGpsPoints, pvstring pvUsedStns)
{
	std::string station;
	std::vector< CDnaGpsPoint >::iterator _it_msr(vGpsPoints->begin());
	_it_pair_vstring it_msr_stns;

	for (_it_msr=vGpsPoints->begin(); _it_msr!=vGpsPoints->end(); _it_msr++)
	{
		station = _it_msr->GetFirst();
		it_msr_stns = equal_range(pvUsedStns->begin(), pvUsedStns->end(), station);
		if (it_msr_stns.first == it_msr_stns.second)
		{
			pvUsedStns->push_back(station);
			std::sort(pvUsedStns->begin(), pvUsedStns->end());
		}
	}
}

void dna_import::ExtractAssociatedStns_D(std::vector<CDnaDirection>* vDirections, pvstring pvUsedStns)
{
	std::string station;
	std::vector<CDnaDirection>::iterator _it_msr(vDirections->begin());
	_it_pair_vstring it_msr_stns;

	for (_it_msr=vDirections->begin(); _it_msr!=vDirections->end(); _it_msr++)
	{
		station = _it_msr->GetFirst();
		it_msr_stns = equal_range(pvUsedStns->begin(), pvUsedStns->end(), station);
		if (it_msr_stns.first == it_msr_stns.second)
		{
			pvUsedStns->push_back(station);
			std::sort(pvUsedStns->begin(), pvUsedStns->end());
		}

		station = _it_msr->GetTarget();
		it_msr_stns = equal_range(pvUsedStns->begin(), pvUsedStns->end(), station);
		if (it_msr_stns.first == it_msr_stns.second)
		{
			pvUsedStns->push_back(station);
			std::sort(pvUsedStns->begin(), pvUsedStns->end());
		}
	}
}

void dna_import::ExtractAssociatedStns(vdnaMsrPtr* vMeasurements, pvstring pvUsedStns)
{
	_it_vdnamsrptr _it_msr;
	std::vector<CDnaDirection>* vdirns;
	std::vector<CDnaGpsBaseline>* vgpsBsls;
	std::vector<CDnaGpsPoint>* vgpsPnts;

	std::string station;

	_it_pair_vstring it_msr_stns;

	for (_it_msr=vMeasurements->begin(); _it_msr!=vMeasurements->end(); ++_it_msr)
	{
		// 1. Handle nested type measurements (D, G, X, Y) separately
		switch (_it_msr->get()->GetTypeC())
		{
		case 'G':	// GPS Baseline (treat as single-baseline cluster)
		case 'X':	// GPS Baseline cluster
			vgpsBsls = _it_msr->get()->GetBaselines_ptr();
			ExtractAssociatedStns_GX(vgpsBsls, pvUsedStns);
			continue;
		case 'Y':	// GPS point cluster
			vgpsPnts = _it_msr->get()->GetPoints_ptr();
			ExtractAssociatedStns_Y(vgpsPnts, pvUsedStns);
			continue;
		}

		if (_it_msr->get()->GetFirst().empty())
		{
			// size_t msr_no(std::distance(vMeasurements->begin(), _it_msr));
			throw XMLInteropException("Empty \"First\" station name.", 0);
		}

		// 2. All other measurements which have <First>
		station = _it_msr->get()->GetFirst();
		it_msr_stns = equal_range(pvUsedStns->begin(), pvUsedStns->end(), station);
		if (it_msr_stns.first == it_msr_stns.second)
		{
			pvUsedStns->push_back(station);
			std::sort(pvUsedStns->begin(), pvUsedStns->end());
		}
		
		// Finished with single station measurements
		switch (_it_msr->get()->GetTypeC())
		{
		case 'H':	// Orthometric height
		case 'R':	// Ellipsoidal height
		case 'I':	// Astronomic latitude
		case 'J':	// Astronomic longitude
		case 'P':	// Geodetic latitude
		case 'Q':	// Geodetic longitude
			continue;
		}

		if (_it_msr->get()->GetTarget().empty())
			throw XMLInteropException("Empty \"Second\" station name", 0);

		// 3. All measurements which have <Target> or <Second>
		station = _it_msr->get()->GetTarget();
		it_msr_stns = equal_range(pvUsedStns->begin(), pvUsedStns->end(), station);
		if (it_msr_stns.first == it_msr_stns.second)
		{
			pvUsedStns->push_back(station);
			std::sort(pvUsedStns->begin(), pvUsedStns->end());
		}

		// Dual station measurements...
		switch (_it_msr->get()->GetTypeC())
		{
		case 'D':	// Direction set
			// now map the directions in the set
			vdirns = _it_msr->get()->GetDirections_ptr();
			ExtractAssociatedStns_D(vdirns, pvUsedStns);
			continue;
		case 'B':	// Geodetic azimuth
		case 'K':	// Astronomic azimuth
		case 'C':	// Chord dist
		case 'E':	// Ellipsoid arc
		case 'M':	// MSL arc
		case 'S':	// Slope distance
		case 'L':	// Level difference
		case 'V':	// Zenith distance
		case 'Z':	// Vertical angle
			continue;
		}

		if (_it_msr->get()->GetTarget2().empty())
			throw XMLInteropException("Empty \"Third\" station name", 0);

		// 4. All measurements which have <Target2> or <Third>
		station = _it_msr->get()->GetTarget2();
		it_msr_stns = equal_range(pvUsedStns->begin(), pvUsedStns->end(), station);
		if (it_msr_stns.first == it_msr_stns.second)
		{
			pvUsedStns->push_back(station);
			std::sort(pvUsedStns->begin(), pvUsedStns->end());
		}
	}

}
	

void dna_import::LoadNetworkFiles(pvstn_t binaryStn, pvmsr_t binaryMsr, 
									const project_settings& projectSettings, bool loadSegmentFile)
{
	projectSettings_ = projectSettings;

	LoadBinaryFiles(binaryStn, binaryMsr);
	if (loadSegmentFile)
		LoadSegmentationFile(binaryMsr);
}

void dna_import::LoadBinaryFiles(pvstn_t binaryStn, pvmsr_t binaryMsr)
{
	g_parsestn_tally.initialise();
	g_parsemsr_tally.initialise();

	try {
		// Load binary stations data.  Throws runtime_error on failure.
		BstFile bst;
		bst.LoadFile(projectSettings_.i.bst_file, binaryStn, bst_meta_);

		// Load binary stations data.  Throws runtime_error on failure.
		BmsFile bms;
		bms.LoadFile(projectSettings_.i.bms_file, binaryMsr, bms_meta_);
	}
	catch (const std::runtime_error& e) {
		SignalExceptionInterop(e.what(), 0, NULL);
	}	
}

void dna_import::LoadSegmentationFile(pvmsr_t binaryMsr)
{
	UINT32 blockCount, blockThreshold, minInnerStns;

	try {
		// Load segmentation file.  Throws runtime_error on failure.
		SegFile seg;
		seg.LoadSegFile(projectSettings_.i.seg_file, 
			blockCount, blockThreshold, minInnerStns,
			v_ISL_, v_JSL_, v_CML_,
			true, binaryMsr,
			&v_measurementCount_, &v_unknownsCount_, &v_ContiguousNetList_, &v_parameterStationCount_);
	}
	catch (const std::runtime_error& e) {
		SignalExceptionInterop(e.what(), 0, NULL);
	}
}

// First item in the file is a UINT32 value - the number of records in the file
// All records are of type UINT32
void dna_import::LoadDatabaseId()
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

void dna_import::ImportStnsMsrsFromNetwork(vdnaStnPtr* vStations, vdnaMsrPtr* vMeasurements, const project_settings& p)
{
	vstn_t binaryStn;
	vmsr_t binaryMsr;

	LoadNetworkFiles(&binaryStn, &binaryMsr, p, true);

	try {
		// Load Database IDs
		LoadDatabaseId();
	}
	catch (const std::runtime_error& e) {
		throw XMLInteropException(e.what(), 0);
	}

	if ((p.i.import_network_number + 1) > v_ContiguousNetList_.size())
	{
		std::stringstream ss;
		ss << "The specified network ID " << p.i.import_network_number <<
			" does not exist in the list of " << std::endl <<
			"         identified contiguous Network IDs (0 to " << (v_ContiguousNetList_.size() - 1) << ").";
		throw XMLInteropException(ss.str(), 0);
	}

	it_vUINT32 _it_data, _it_netid;
	UINT32 segmentedBlock(0);

	dnaStnPtr stnPtr(new CDnaStation(datum_.GetName(), datum_.GetEpoch_s()));

	vStations->clear();

	for (_it_netid = v_ContiguousNetList_.begin(); _it_netid != v_ContiguousNetList_.end(); ++_it_netid)
	{
		if (*_it_netid != p.i.import_network_number)
		{
			segmentedBlock++;
			continue;
		}

		// get stations on ISL
		for (_it_data = v_ISL_.at(segmentedBlock).begin();
			_it_data != v_ISL_.at(segmentedBlock).end();
			++_it_data)
		{
			stnPtr->SetStationRec(binaryStn.at(*_it_data));
			vStations->push_back(stnPtr);
			g_parsestn_tally.addstation(stnPtr->GetConstraints());
			stnPtr.reset(new CDnaStation(datum_.GetName(), datum_.GetEpoch_s()));
		}

		// get stations on JSL
		for (_it_data = v_JSL_.at(segmentedBlock).begin();
			_it_data != v_JSL_.at(segmentedBlock).end();
			++_it_data)
		{
			stnPtr->SetStationRec(binaryStn.at(*_it_data));
			vStations->push_back(stnPtr);
			g_parsestn_tally.addstation(stnPtr->GetConstraints());
			stnPtr.reset(new CDnaStation(datum_.GetName(), datum_.GetEpoch_s()));
		}

		segmentedBlock++;
	}

	dnaMsrPtr msrPtr;
	msrPtr.reset();

	vMeasurements->clear();

	it_vmsr_t it_msr;
	size_t dbindex;
	it_vdbid_t it_dbid;
	UINT32 msr_no(0);

	segmentedBlock = 0;

	for (_it_netid = v_ContiguousNetList_.begin(); _it_netid != v_ContiguousNetList_.end(); ++_it_netid)
	{
		if (*_it_netid != p.i.import_network_number)
		{
			segmentedBlock++;
			continue;
		}

		// remove measurements that are in a block other than those 
		// belonging to network 0
		RemoveNonMeasurements(segmentedBlock, &binaryMsr);

		// get stations on CML
		for (_it_data = v_CML_.at(segmentedBlock).begin();
			_it_data != v_CML_.at(segmentedBlock).end();
			++_it_data)
		{
			it_msr = binaryMsr.begin() + *_it_data;
			ResetMeasurementPtr(&msrPtr, binaryMsr.at(*_it_data).measType);
			msr_no++;

			// build measurement tally - do determine whether measurements of a particular type
			// have been supplied
			switch (it_msr->measType)
			{
			case 'A': // Horizontal angle
				g_parsemsr_tally.A++;
				break;
			case 'B': // Geodetic azimuth
				g_parsemsr_tally.B++;
				break;
			case 'C': // Chord dist
				g_parsemsr_tally.C++;
				break;
			case 'D': // Direction set
				if (it_msr->measStart == xMeas)
					g_parsemsr_tally.D += it_msr->vectorCount1;
				break;
			case 'E': // Ellipsoid arc
				g_parsemsr_tally.E++;
				break;
			case 'G': // GPS Baseline
				g_parsemsr_tally.G += 3;
				break;
			case 'H': // Orthometric height
				g_parsemsr_tally.H++;
				break;
			case 'I': // Astronomic latitude
				g_parsemsr_tally.I++;
				break;
			case 'J': // Astronomic longitude
				g_parsemsr_tally.J++;
				break;
			case 'K': // Astronomic azimuth
				g_parsemsr_tally.K++;
				break;
			case 'L': // Level difference
				g_parsemsr_tally.L++;
				break;
			case 'M': // MSL arc
				g_parsemsr_tally.M++;
				break;
			case 'P': // Geodetic latitude
				g_parsemsr_tally.P++;
				break;
			case 'Q': // Geodetic longitude
				g_parsemsr_tally.Q++;
				break;
			case 'R': // Ellipsoidal height
				g_parsemsr_tally.R++;
				break;
			case 'S': // Slope distance
				g_parsemsr_tally.S++;
				break;
			case 'V': // Zenith distance
				g_parsemsr_tally.V++;
				break;
			case 'X': // GPS Baseline cluster
				if (it_msr->measStart == xMeas)
					g_parsemsr_tally.X += it_msr->vectorCount1 * 3;
				break;
			case 'Y': // GPS point cluster
				if (it_msr->measStart == xMeas)
					g_parsemsr_tally.Y += it_msr->vectorCount1 * 3;
				break;
			case 'Z': // Vertical angle
				g_parsemsr_tally.Z++;
				break;
			}

			if (databaseIDsSet_)
			{
				dbindex = std::distance(binaryMsr.begin(), it_msr);
				it_dbid = v_msr_db_map_.begin() + dbindex;
			}

			msrPtr->SetMeasurementRec(binaryStn, it_msr, it_dbid);
			msrPtr->ResolveSourceFile(bms_meta_.sourceFileMeta, bms_meta_.sourceFileCount);
			vMeasurements->push_back(msrPtr);
		}

		segmentedBlock++;
	}
}
	

void dna_import::ImportStnsMsrsFromBlock(vdnaStnPtr* vStations, vdnaMsrPtr* vMeasurements, const project_settings& p)
{
	vstn_t binaryStn;
	vmsr_t binaryMsr;

	LoadNetworkFiles(&binaryStn, &binaryMsr, p, true);

	try {
		// Load Database IDs
		LoadDatabaseId();
	}
	catch (const std::runtime_error& e) {
		throw XMLInteropException(e.what(), 0);
	}

	if (p.i.import_block_number > v_ISL_.size())
	{
		std::stringstream ss;
		ss << "The specified block number " << p.i.import_block_number << 
			" exceeds the total number of " << std::endl <<
			"         segmented blocks (" << v_ISL_.size() << ").";
		throw XMLInteropException(ss.str(), 0);
	}

	it_vUINT32 _it_data;

	dnaStnPtr stnPtr(new CDnaStation(datum_.GetName(), datum_.GetEpoch_s()));

	vStations->clear();

	// get stations on ISL
	for (_it_data=v_ISL_.at(p.i.import_block_number-1).begin();
		_it_data!=v_ISL_.at(p.i.import_block_number-1).end();
		++_it_data)
	{
		stnPtr->SetStationRec(binaryStn.at(*_it_data));
		vStations->push_back(stnPtr);
		g_parsestn_tally.addstation(stnPtr->GetConstraints());
		stnPtr.reset(new CDnaStation(datum_.GetName(), datum_.GetEpoch_s()));
	}

	// get stations on JSL
	for (_it_data=v_JSL_.at(p.i.import_block_number-1).begin();
		_it_data!=v_JSL_.at(p.i.import_block_number-1).end();
		++_it_data)
	{
		stnPtr->SetStationRec(binaryStn.at(*_it_data));
		vStations->push_back(stnPtr);
		g_parsestn_tally.addstation(stnPtr->GetConstraints());
		stnPtr.reset(new CDnaStation(datum_.GetName(), datum_.GetEpoch_s()));
	}

	dnaMsrPtr msrPtr;
	msrPtr.reset();

	vMeasurements->clear();

	it_vmsr_t it_msr;
	size_t dbindex;
	it_vdbid_t it_dbid;
	UINT32 msr_no(0);

	RemoveNonMeasurements(p.i.import_block_number-1, &binaryMsr);

	
	// get stations on CML
	for (_it_data=v_CML_.at(p.i.import_block_number-1).begin();
		_it_data!=v_CML_.at(p.i.import_block_number-1).end();
		++_it_data)
	{
		it_msr = binaryMsr.begin() + *_it_data;
		ResetMeasurementPtr(&msrPtr, binaryMsr.at(*_it_data).measType);
		msr_no++;

		// build measurement tally - do determine whether measurements of a particular type
		// have been supplied
		switch (it_msr->measType)
		{
		case 'A': // Horizontal angle
			g_parsemsr_tally.A++;
			break;
		case 'B': // Geodetic azimuth
			g_parsemsr_tally.B++;
			break;
		case 'C': // Chord dist
			g_parsemsr_tally.C++;
			break;
		case 'D': // Direction set
			if (it_msr->measStart == xMeas)
				g_parsemsr_tally.D += it_msr->vectorCount1;
			break;
		case 'E': // Ellipsoid arc
			g_parsemsr_tally.E++;
			break;
		case 'G': // GPS Baseline
			g_parsemsr_tally.G += 3;
			break;
		case 'H': // Orthometric height
			g_parsemsr_tally.H++;
			break;
		case 'I': // Astronomic latitude
			g_parsemsr_tally.I++;
			break;
		case 'J': // Astronomic longitude
			g_parsemsr_tally.J++;
			break;
		case 'K': // Astronomic azimuth
			g_parsemsr_tally.K++;
			break;
		case 'L': // Level difference
			g_parsemsr_tally.L++;
			break;
		case 'M': // MSL arc
			g_parsemsr_tally.M++;
			break;
		case 'P': // Geodetic latitude
			g_parsemsr_tally.P++;
			break;
		case 'Q': // Geodetic longitude
			g_parsemsr_tally.Q++;
			break;
		case 'R': // Ellipsoidal height
			g_parsemsr_tally.R++;
			break;
		case 'S': // Slope distance
			g_parsemsr_tally.S++;
			break;
		case 'V': // Zenith distance
			g_parsemsr_tally.V++;
			break;
		case 'X': // GPS Baseline cluster
			if (it_msr->measStart == xMeas)
				g_parsemsr_tally.X += it_msr->vectorCount1 * 3;
			break;
		case 'Y': // GPS point cluster
			if (it_msr->measStart == xMeas)
				g_parsemsr_tally.Y += it_msr->vectorCount1 * 3;
			break;
		case 'Z': // Vertical angle
			g_parsemsr_tally.Z++;
			break;
		}

		if (databaseIDsSet_)
		{
			dbindex = std::distance(binaryMsr.begin(), it_msr);
 			it_dbid = v_msr_db_map_.begin() + dbindex;
		}

		msrPtr->SetMeasurementRec(binaryStn, it_msr, it_dbid);
		msrPtr->ResolveSourceFile(bms_meta_.sourceFileMeta, bms_meta_.sourceFileCount);
		vMeasurements->push_back(msrPtr);
	}
}


void dna_import::RemoveNonMeasurements(const UINT32& block, pvmsr_t binaryMsr)
{
	if (v_CML_.at(block).size() < 2)
		return;
	CompareNonMeasStart<measurement_t, UINT32> measstartCompareFunc(binaryMsr, xMeas);
	std::sort(v_CML_.at(block).begin(), v_CML_.at(block).end(), measstartCompareFunc);
	erase_if(v_CML_.at(block), measstartCompareFunc);
	CompareMsrFileOrder<measurement_t, UINT32> fileorderCompareFunc(binaryMsr);
	std::sort(v_CML_.at(block).begin(), v_CML_.at(block).end(), fileorderCompareFunc);
}
	

void dna_import::SignalComplete()
{
	isProcessing_ = false;
	percentComplete_ = -99.0;
	
	// Obtain exclusive use of the input file pointer
	import_file_mutex.lock();
	
	try {
		if (ifsInputFILE_ != 0)
		{
			if (ifsInputFILE_->good())
				ifsInputFILE_->close();
			delete ifsInputFILE_;
		}
	}
	catch (const std::ifstream::failure& e) 
	{
		if (ifsInputFILE_->rdstate() & std::ifstream::eofbit)
		{
			ifsInputFILE_ = 0;
		}
		throw XMLInteropException(e.what(), 0);
	}
	
	ifsInputFILE_ = 0;

	// release file pointer mutex
	import_file_mutex.unlock();
}

void dna_import::SignalExceptionParse(std::string msg, int i)
{
	SignalComplete();
	size_t s = msg.find("is not allowed for content model");
	if (s != std::string::npos)
		msg += "\n  - check that the order of elements in the XML file matches the order of the XSD complex type elements.";
	parseStatus_ = PARSE_EXCEPTION_RAISED;
	throw XMLInteropException(msg, i);
}
	
void dna_import::SignalExceptionParseDNA(const std::string& message, const std::string& sBuf, const int& column_no)
{
	std::stringstream ss;
	ss << message << std::endl;
	if (!sBuf.empty())
		ss << "    " << sBuf << std::endl;
	m_columnNo = column_no + 1;
	throw XMLInteropException(ss.str(), m_lineNo);
}
	

// Name:				SignalExceptionInterop
// Purpose:				Closes all files (if file pointers are passed in) and throws XMLInteropException
// Called by:			Any
// Calls:				XMLInteropException()
void dna_import::SignalExceptionInterop(std::string msg, int i, const char *streamType, ...)
{
	if (streamType == NULL)
		throw XMLInteropException(msg, i);

	std::ofstream* outstream;
	std::ifstream* instream;

	va_list argptr; 
	va_start(argptr, streamType);

	while (*streamType != '\0')
	{
		//ifstream
		if (*streamType == 'i' )
		{
			instream = va_arg(argptr, std::ifstream*);
			instream->close();
		}
		//ofstream
		if (*streamType == 'o' )
		{
			outstream = va_arg(argptr, std::ofstream*);
			outstream->close();
		}
		streamType++;
	}
	
	va_end(argptr);
	
	throw XMLInteropException(msg, i);
}

void dna_import::SerialiseMSR(vdnaMsrPtr* vMeasurements, const std::string& msrfilename, 
	const project_settings& p, vifm_t* vinput_file_meta)
{
	// Measurements
	std::ofstream dna_msr_file;
	try {
		// Create DynAdjust MSR file. 
		file_opener(dna_msr_file, msrfilename);
	}
	catch (const std::runtime_error& e) {
		SignalExceptionInterop(e.what(), 0, NULL);
	}

	_it_vdnamsrptr _it_msr(vMeasurements->begin());

	try {
		// Write version line
		dna_header(dna_msr_file, "3.01", "MSR", datum_.GetName(), datum_.GetEpoch_s(), vMeasurements->size());

		// Capture source files
		std::string source_files(formatStnMsrFileSourceString<std::string>(vinput_file_meta, msr_data));

		// Write header comment lines about this file
		dna_comment(dna_msr_file, "File type:    Measurement file");
		dna_comment(dna_msr_file, "Project name: " + p.g.network_name);
		dna_comment(dna_msr_file, "Source files: " + source_files);

		// print measurements
		for (_it_msr = vMeasurements->begin(); _it_msr != vMeasurements->end(); _it_msr++)
			_it_msr->get()->WriteDNAMsr(&dna_msr_file, dmw_, dml_);

		dna_msr_file.close();

	}
	catch (const std::ifstream::failure& f) {
		SignalExceptionInterop(static_cast<std::string>(f.what()), 0, "o", &dna_msr_file);
	}
	catch (const XMLInteropException& e) {
		SignalExceptionInterop(static_cast<std::string>(e.what()), 0, "o", &dna_msr_file);
	}
}
	

void dna_import::SerialiseSTN(vdnaStnPtr* vStations, const std::string& stnfilename, 
	const project_settings& p, vifm_t* vinput_file_meta, bool flagUnused)
{
	// Stations
	std::ofstream dna_stn_file;
	try {
		// Create DynAdjust STN file. 
		file_opener(dna_stn_file, stnfilename);
	}
	catch (const std::runtime_error& e) {
		SignalExceptionInterop(e.what(), 0, NULL);
	}

	determineDNASTNFieldParameters<UINT16>("3.01", dsl_, dsw_);
	determineDNAMSRFieldParameters<UINT16>("3.01", dml_, dmw_);

	_it_vdnastnptr _it_stn(vStations->begin());
	CDnaProjection projection(UTM);

	try {
		// print stations
		// Has the user specified --flag-unused-stations, in which case, do not
		// print stations marked unused?
		UINT32 count(0);
		if (flagUnused)
		{
			for_each(vStations->begin(), vStations->end(),
				[&count](const dnaStnPtr& stn) {
					if (stn.get()->IsNotUnused())
						count++;
				});
		}
		else
			count = static_cast<UINT32>(vStations->size());

		// Write version line
		dna_header(dna_stn_file, "3.01", "STN", datum_.GetName(), datum_.GetEpoch_s(), count);

		// Capture source files
		std::string source_files(formatStnMsrFileSourceString<std::string>(vinput_file_meta, stn_data));

		// Write header comment lines about this file
		dna_comment(dna_stn_file, "File type:    Station file");
		dna_comment(dna_stn_file, "Project name: " + p.g.network_name);
		dna_comment(dna_stn_file, "Source files: " + source_files);

		// print stations
		// Has the user specified --flag-unused-stations, in wich case, do not
		// print stations marked unused?
		if (flagUnused)
		{
			for (_it_stn = vStations->begin(); _it_stn != vStations->end(); _it_stn++)
				if (_it_stn->get()->IsNotUnused())
					_it_stn->get()->WriteDNAXMLStnInitialEstimates(&dna_stn_file,
						datum_.GetEllipsoidRef(), &projection,
						dna, &dsw_);
		}
		else
		{
			// print all stations regardless of whether they are unused or not
			for (_it_stn = vStations->begin(); _it_stn != vStations->end(); _it_stn++)
				_it_stn->get()->WriteDNAXMLStnInitialEstimates(&dna_stn_file,
					datum_.GetEllipsoidRef(), &projection,
					dna, &dsw_);
		}

		dna_stn_file.close();

	}
	catch (const std::ifstream::failure& f) {
		SignalExceptionInterop(static_cast<std::string>(f.what()), 0, "o", &dna_stn_file);
	}
	catch (const XMLInteropException& e) {
		SignalExceptionInterop(static_cast<std::string>(e.what()), 0, "o", &dna_stn_file);
	}
}

void dna_import::SerialiseDNA(vdnaStnPtr* vStations, vdnaMsrPtr* vMeasurements,
	const std::string& stnfilename, const std::string& msrfilename, 
	const project_settings& p, vifm_t* vinput_file_meta, bool flagUnused)
{
	try
	{
		// Reset the default datum.
		datum_.SetDatumFromEpsg(m_strProjectDefaultEpsg, m_strProjectDefaultEpoch);
	}
	catch (const std::runtime_error& e) {
		SignalExceptionInterop(e.what(), 0, NULL);
	}

	// Write stations (if present)
	if (!vStations->empty())
		SerialiseSTN(vStations, stnfilename, p, vinput_file_meta, flagUnused);

	// Write measurements (if present)
	if (!vMeasurements->empty())
		SerialiseMSR(vMeasurements, msrfilename, p, vinput_file_meta);
}
	

void dna_import::InitialiseDynaMLFile(const project_settings& p, vifm_t* vinput_file_meta,
	const std::string& outfilename, std::ofstream* dynaml_file)
{
	// create combined filename
	std::string dynamlfilename(outfilename);
	size_t dp = std::string::npos;
	if ((dp = dynamlfilename.rfind(".")) == std::string::npos)
		dynamlfilename += ".xml";

	try
	{
		// Reset the default datum.
		datum_.SetDatumFromEpsg(m_strProjectDefaultEpsg, m_strProjectDefaultEpoch);
	}
	catch (const std::runtime_error& e) {
		SignalExceptionInterop(e.what(), 0, NULL);
	}

	try {
		// Create DynaML station file.
		file_opener(*dynaml_file, dynamlfilename);

		// write header
		dynaml_header(*dynaml_file, "Combined File", datum_.GetName(), datum_.GetEpoch_s());

		// Capture source files
		std::string source_files(formatStnMsrFileSourceString<std::string>(vinput_file_meta, stn_msr_data));

		// Write header comment lines about this file
		dynaml_comment(*dynaml_file, "File type:    Combined station and measurements file");
		dynaml_comment(*dynaml_file, "Project name: " + p.g.network_name);
		dynaml_comment(*dynaml_file, "Source files: " + source_files);		
	}
	catch (const std::runtime_error& e) {
		SignalExceptionInterop(e.what(), 0, NULL);
	}
}

// A useful, but unnecessary function	
//
//void dna_import::SerialiseDynaMLfromBinary(const std::string& outfilename, 
//	const project_settings& p, vifm_t* vinput_file_meta, bool flagUnused)
//{
//	// Load stations from binary file and serialise to XML
//	std::ifstream bst_file;
//	try {
//		// Load binary stations data.  Throws runtime_error on failure.
//		file_opener(bst_file, p.i.bst_file, std::ios::in | std::ios::binary, binary, true);
//	}
//	catch (const std::runtime_error& e) {
//		SignalExceptionInterop(e.what(), 0, NULL);
//	}
//
//	std::ofstream dynaml_file;
//	InitialiseDynaMLFile(p, vinput_file_meta, outfilename, &dynaml_file);
//
//	// Write the stations (from memory)
//	try {
//		SerialiseXmlStn(&bst_file, &dynaml_file);
//	}
//	catch (const std::ios_base::failure& f) {
//		SignalExceptionInterop(static_cast<std::string>(f.what()), 0, 
//			"io", &bst_file, &dynaml_file);	
//	}
//
//	std::ifstream bms_file;
//	try {
//		// Load binary measurements data.  Throws runtime_error on failure.
//		file_opener(bms_file, p.i.bms_file, std::ios::in | std::ios::binary, binary, true);
//	}
//	catch (const std::runtime_error& e) {
//		SignalExceptionInterop(e.what(), 0, NULL);
//	}
//
//	// Write the measurements (from binary file)
//	try {
//		SerialiseXmlMsr(&bst_file, &bms_file, &dynaml_file);
//	}
//	catch (const std::ifstream::failure& f) {
//		SignalExceptionInterop(static_cast<std::string>(f.what()), 0, 
//			"iio", &bst_file, &bms_file, &dynaml_file);
//	}
//	catch (const XMLInteropException& e)  {
//		SignalExceptionInterop(static_cast<std::string>(e.what()), 0, 
//			"iio", &bst_file, &bms_file, &dynaml_file);
//	}
//
//	bst_file.close();
//	bms_file.close();
//	dynaml_file << "</DnaXmlFormat>" << std::endl;
//	dynaml_file.close();
//}


void dna_import::SerialiseDynaMLfromMemory(vdnaStnPtr* vStations, vdnaMsrPtr* vMeasurements, 
	const std::string& outfilename, const project_settings& p, vifm_t* vinput_file_meta, bool flagUnused)
{
	// create combined filename
	std::string dynamlfilename(outfilename);
	size_t dp = std::string::npos;
	if ((dp = dynamlfilename.rfind(".")) == std::string::npos)
		dynamlfilename += ".xml";
	
	try
	{
		// Reset the default datum.
		datum_.SetDatumFromEpsg(m_strProjectDefaultEpsg, m_strProjectDefaultEpoch);
	}
	catch (const std::runtime_error& e) {
		SignalExceptionInterop(e.what(), 0, NULL);
	}

	std::ofstream dynaml_file;
	InitialiseDynaMLFile(p, vinput_file_meta, outfilename, &dynaml_file);

	CDnaProjection projection(UTM);
	std::string comment("");

	try {
		// Write the stations (from memory)
		if (flagUnused) 
		{
			for_each(vStations->begin(), vStations->end(), 
				[this, &projection, &dynaml_file] (dnaStnPtr stn) {
					if (stn->IsNotUnused())
						stn->WriteDNAXMLStnInitialEstimates(&dynaml_file,
							datum_.GetEllipsoidRef(), &projection, dynaml);
			}); 
		}
		else
		{
			for_each(vStations->begin(), vStations->end(), 
				[this, &projection, &dynaml_file] (dnaStnPtr stn) {
					stn->WriteDNAXMLStnInitialEstimates(&dynaml_file,
						datum_.GetEllipsoidRef(), &projection, dynaml);
			});
		}
		
		// Write the measurements (from memory)
		for_each(vMeasurements->begin(), vMeasurements->end(), 
			[this, &projection, &comment, &dynaml_file] (dnaMsrPtr msr) {
				msr->WriteDynaMLMsr(&dynaml_file, comment);
		});
	}
	catch (const std::ifstream::failure& f) {
		SignalExceptionInterop(static_cast<std::string>(f.what()), 0, 
			"o", &dynaml_file);
	}
	catch (const XMLInteropException& e)  {
		SignalExceptionInterop(static_cast<std::string>(e.what()), 0, 
			"o", &dynaml_file);
	}
	
	dynaml_file << "</DnaXmlFormat>" << std::endl;
	dynaml_file.close();
}
	

void dna_import::InitialiseDynaMLSepStationFile(const project_settings& p, vifm_t* vinput_file_meta,
	const std::string& stnfilename, std::ofstream* dynaml_stn_file)
{
	try
	{
		// Reset the default datum.
		datum_.SetDatumFromEpsg(m_strProjectDefaultEpsg, m_strProjectDefaultEpoch);
	}
	catch (const std::runtime_error& e) {
		SignalExceptionInterop(e.what(), 0, NULL);
	}

	try {
		// Create DynaML station file. 
		file_opener(*dynaml_stn_file, stnfilename);

		// write header
		dynaml_header(*dynaml_stn_file, "Station File", datum_.GetName(), datum_.GetEpoch_s());

		// Capture source files
		std::string source_files(formatStnMsrFileSourceString<std::string>(vinput_file_meta, stn_data));

		// Write header comment lines about this file
		dynaml_comment(*dynaml_stn_file, "File type:    Station file");
		dynaml_comment(*dynaml_stn_file, "Project name: " + p.g.network_name);
		dynaml_comment(*dynaml_stn_file, "Source files: " + source_files);
	}
	catch (const std::runtime_error& e) {
		SignalExceptionInterop(e.what(), 0, NULL);
	}
}

void dna_import::InitialiseDynaMLSepMeasurementFile(const project_settings& p, vifm_t* vinput_file_meta,
	const std::string& msrfilename, std::ofstream* dynaml_msr_file)
{
	try {
		// Create DynaML measurement file. 
		file_opener(*dynaml_msr_file, msrfilename);

		// write header
		dynaml_header(*dynaml_msr_file, "Measurement File", datum_.GetName(), datum_.GetEpoch_s());

		// Capture source files
		std::string source_files(formatStnMsrFileSourceString<std::string>(vinput_file_meta, msr_data));

		// Write header comment lines about this file
		dynaml_comment(*dynaml_msr_file, "File type:    Measurement file");
		dynaml_comment(*dynaml_msr_file, "Project name: " + p.g.network_name);
		dynaml_comment(*dynaml_msr_file, "Source files: " + source_files);
	}
	catch (const std::runtime_error& e) {
		SignalExceptionInterop(e.what(), 0, NULL);
	}
}

// A useful, but unnecessary function	
//
//void dna_import::SerialiseDynaMLSepfromBinary(const std::string& stnfilename, const std::string& msrfilename, 
//	const project_settings& p, vifm_t* vinput_file_meta, bool flagUnused)
//{
//	// Load stations from binary file and serialise to XML
//	std::ifstream bst_file;
//	try {
//		// Load binary stations data.  Throws runtime_error on failure.
//		file_opener(bst_file, p.i.bst_file, std::ios::in | std::ios::binary, binary, true);
//	}
//	catch (const std::runtime_error& e) {
//		SignalExceptionInterop(e.what(), 0, NULL);
//	}
//
//	// Create Station file and initialise header
//	std::ofstream dynaml_stn_file;
//	InitialiseDynaMLSepStationFile(p, vinput_file_meta, stnfilename, &dynaml_stn_file);
//
//	// Write the stations (from binary file)
//	try {		
//		SerialiseXmlStn(&bst_file, &dynaml_stn_file);
//		dynaml_stn_file << "</DnaXmlFormat>" << std::endl;
//		dynaml_stn_file.close();
//	}
//	catch (const std::ios_base::failure& f) {
//		SignalExceptionInterop(static_cast<std::string>(f.what()), 0, 
//			"io", &bst_file, &dynaml_stn_file);	
//	}
//
//	std::ifstream bms_file;
//	try {
//		// Load binary measurements data.  Throws runtime_error on failure.
//		file_opener(bms_file, p.i.bms_file, std::ios::in | std::ios::binary, binary, true);
//	}
//	catch (const std::runtime_error& e) {
//		SignalExceptionInterop(e.what(), 0, NULL);
//	}
//
//	// Create Measurement file and initialise header
//	std::ofstream dynaml_msr_file;
//	InitialiseDynaMLSepMeasurementFile(p, vinput_file_meta, msrfilename, &dynaml_msr_file);
//	
//	// Write the measurements (from binary file)
//	try {
//		SerialiseXmlMsr(&bst_file, &bms_file, &dynaml_msr_file);
//		dynaml_msr_file << "</DnaXmlFormat>" << std::endl;
//		dynaml_msr_file.close();
//	}
//	catch (const std::ifstream::failure& f) {
//		SignalExceptionInterop(static_cast<std::string>(f.what()), 0, 
//			"iio", &bst_file, &bms_file, &dynaml_msr_file);
//	}
//	catch (const XMLInteropException& e)  {
//		SignalExceptionInterop(static_cast<std::string>(e.what()), 0, 
//			"iio", &bst_file, &bms_file, &dynaml_msr_file);
//	}
//
//	bst_file.close();
//	bms_file.close();
//}
	

void dna_import::SerialiseDynaMLSepfromMemory(vdnaStnPtr* vStations, vdnaMsrPtr* vMeasurements, 
	const std::string& stnfilename, const std::string& msrfilename, 
	const project_settings& p, vifm_t* vinput_file_meta, bool flagUnused)
{	
	// Create Station file and initialise header
	std::ofstream dynaml_stn_file;
	InitialiseDynaMLSepStationFile(p, vinput_file_meta, stnfilename, &dynaml_stn_file);

	CDnaProjection projection(UTM);
	std::string comment("");
	
	// Write the stations (from memory)
	try {	
		if (flagUnused)
		{
			for_each(vStations->begin(), vStations->end(), 
				[this, &projection, &dynaml_stn_file] (dnaStnPtr stn) {
					if (stn->IsNotUnused())
						stn->WriteDNAXMLStnInitialEstimates(&dynaml_stn_file,
						datum_.GetEllipsoidRef(), &projection, dynaml);
			}); 
		}
		else
		{
			for_each(vStations->begin(), vStations->end(), 
				[this, &projection, &dynaml_stn_file] (dnaStnPtr stn) {
					stn->WriteDNAXMLStnInitialEstimates(&dynaml_stn_file,
						datum_.GetEllipsoidRef(), &projection, dynaml);
			});
		}

		dynaml_stn_file << "</DnaXmlFormat>" << std::endl;
		dynaml_stn_file.close();

	}
	catch (const std::ifstream::failure& f) {
		SignalExceptionInterop(static_cast<std::string>(f.what()), 0, "o", &dynaml_stn_file);
	}
	catch (const XMLInteropException& e)  {
		SignalExceptionInterop(static_cast<std::string>(e.what()), 0, "o", &dynaml_stn_file);
	}

	// Create Measurement file and initialise header
	std::ofstream dynaml_msr_file;
	InitialiseDynaMLSepMeasurementFile(p, vinput_file_meta, msrfilename, &dynaml_msr_file);

	// Write the measurements (from memory)
	try {	
		for_each(vMeasurements->begin(), vMeasurements->end(), 
			[this, &projection, &comment, &dynaml_msr_file] (dnaMsrPtr msr) {
				msr->WriteDynaMLMsr(&dynaml_msr_file, comment);
		});

		dynaml_msr_file << "</DnaXmlFormat>" << std::endl;
		dynaml_msr_file.close();

	}
	catch (const std::ifstream::failure& f) {
		SignalExceptionInterop(static_cast<std::string>(f.what()), 0, "o", &dynaml_msr_file);
	}
	catch (const XMLInteropException& e)  {
		SignalExceptionInterop(static_cast<std::string>(e.what()), 0, "o", &dynaml_msr_file);
	}
}


//void dna_import::SerialiseGeoidData(vdnaStnPtr* vStations, const std::string& geofilename)
//{
//	std::ofstream dynaml_geo_file;
//	try {
//		// Create geoid file.  Throws runtime_error on failure.
//		file_opener(dynaml_geo_file, geofilename);
//	}
//	catch (const std::runtime_error& e) {
//		SignalExceptionInterop(e.what(), 0, NULL);
//	}
//	
//	_it_vdnastnptr _it_stn;
//
//	try {
//		// Write header line
//		dna_comment(dynaml_geo_file, "DNA geoid file.");
//
//		// Print data to geoid file
//		for (_it_stn=vStations->begin(); _it_stn!=vStations->end(); _it_stn++)
//			_it_stn->get()->WriteGeoidfile(&dynaml_geo_file);
//		dynaml_geo_file.close();
//		
//	}
//	catch (const std::ifstream::failure& f) {
//		SignalExceptionInterop(static_cast<std::string>(f.what()), 0, "o", &dynaml_geo_file);
//	}
//	catch (const XMLInteropException& e)  {
//		SignalExceptionInterop(static_cast<std::string>(e.what()), 0, "o", &dynaml_geo_file);
//	}
//}
	

// A useful, but unnecessary function	
//
//void dna_import::SerialiseXmlStn(std::ifstream* ifs_stns, std::ofstream* ofs_dynaml)
//{
//	try
//	{
//		// Reset the default datum.
//		datum_.SetDatumFromEpsg(m_strProjectDefaultEpsg, m_strProjectDefaultEpoch);
//	}
//	catch (const std::runtime_error& e) {
//		SignalExceptionInterop(e.what(), 0, NULL);
//	}
//
//	station_t stationRecord;
//	dnaStnPtr stnPtr(new CDnaStation(datum_.GetName(), datum_.GetEpoch_s()));
//	stnPtr->SetCoordType(LLH_type);
//
//	CDnaProjection projection(UTM);
//
//	// get number of stations
//	UINT32 i, stnCount;
//	ifs_stns->read(reinterpret_cast<char *>(&stnCount), sizeof(UINT32)); 
//
//	for (i=0; i<stnCount; i++)
//	{
//		if (ifs_stns->eof() || !ifs_stns->good())
//			throw XMLInteropException("SerialiseXMLStn(): Errors were encountered when reading from the binary station file.", 0);
//		ifs_stns->read(reinterpret_cast<char *>(&stationRecord), sizeof(station_t)); 
//		stnPtr->SetStationRec(stationRecord);
//		stnPtr->WriteDNAXMLStnInitialEstimates(ofs_dynaml,
//			datum_.GetEllipsoidRef(), &projection,
//			dynaml);
//	}
//}
	

//void dna_import::SerialiseXmlMsr(std::ifstream* ifs_stns, std::ifstream* ifs_msrs, std::ofstream* ofs_dynaml)
//{
//	try
//	{
//		// Reset the default datum.
//		datum_.SetDatumFromEpsg(m_strProjectDefaultEpsg, m_strProjectDefaultEpoch);
//	}
//	catch (const std::runtime_error& e) {
//		SignalExceptionInterop(e.what(), 0, NULL);
//	}
//
//	measurement_t measRecord;
//	dnaMsrPtr msrPtr;
//	msrPtr.reset();
//
//	// get number of measurements
//	UINT32 msrCount;
//	std::string comment("");
//	ifs_msrs->read(reinterpret_cast<char *>(&msrCount), sizeof(UINT32)); 
//
//	for (UINT32 i=0; i<msrCount; i++)
//	{
//		if (ifs_msrs->eof() || !ifs_msrs->good())
//			throw XMLInteropException("SerialiseXMLMsr(): Errors were encountered when reading from the binary measurement file.", 0);
//		ifs_msrs->read(reinterpret_cast<char *>(&measRecord), sizeof(measurement_t));
//		ResetMeasurementPtr(&msrPtr, measRecord.measType);
//		i += msrPtr->SetMeasurementRec(ifs_stns, ifs_msrs, &measRecord); // increment by number of elements read from binary file
//		msrPtr->WriteDynaMLMsr(ofs_dynaml, comment);
//	}
//}
	

void dna_import::SerialiseBms(const std::string& bms_filename, vdnaMsrPtr* vMeasurements,
	vifm_t& vinput_file_meta)
{
	// Calculate number of measurement records
	m_binaryRecordCount = 0;
	//m_dbidRecordCount = 0;
	for_each(vMeasurements->begin(), vMeasurements->end(),
		[this](dnaMsrPtr msr) {
			// update database map
			//msr->SetDatabaseMap_bmsIndex(m_binaryRecordCount);
			m_binaryRecordCount += msr->CalcBinaryRecordCount();
			//m_dbidRecordCount += msr->CalcDbidRecordCount();
	});

	// the database ID vector is set to the same size as the 
	// binary measurement vector to enable efficient 1:1 lookup
	// of database IDs
	m_dbidRecordCount = m_binaryRecordCount;

	BmsFile bms;

	try {
        snprintf(bms_meta_.modifiedBy, sizeof(bms_meta_.modifiedBy), "%s", __BINARY_NAME__);
		bms_meta_.binCount = m_binaryRecordCount;
		bms_meta_.inputFileCount = bms.CreateMsrInputFileMeta(vinput_file_meta, &(bms_meta_.inputFileMeta));
		bms_meta_.reftran = false;
		bms_meta_.geoid = false;
		bms.WriteFile(bms_filename, vMeasurements, bms_meta_);
	}
	catch (const std::runtime_error& e) {
		SignalExceptionInterop(e.what(), 0, NULL);
	}
}
	

void dna_import::SerialiseBst(const std::string& bst_filename, vdnaStnPtr* vStations, 
	pvstring vUnusedStns, vifm_t& vinput_file_meta,
	bool flagUnused)
{	
	BstFile bst;
	
	try {
        snprintf(bms_meta_.modifiedBy, sizeof(bms_meta_.modifiedBy), "%s", __BINARY_NAME__);
		bst_meta_.binCount = static_cast<std::uint64_t>(vStations->size());
		bst_meta_.inputFileCount = bst.CreateStnInputFileMeta(vinput_file_meta, &(bst_meta_.inputFileMeta));
		bst_meta_.reduced = true;
		bst_meta_.reftran = false;
		bst_meta_.geoid = false;
		bst.WriteFile(bst_filename, vStations, vUnusedStns, bst_meta_, flagUnused);
	}
	catch (const std::runtime_error& e) {
		SignalExceptionInterop(e.what(), 0, NULL);
	}
}
	

// First item in the file is a UINT32 value - the number of records in the file
// All records are of type UINT32
void dna_import::SerialiseDatabaseId(const std::string& dbid_filename, pvdnaMsrPtr vMeasurements)
{
	std::ofstream dbid_file;
	try {
		// Create geoid file.  Throws runtime_error on failure.
		file_opener(dbid_file, dbid_filename,
			std::ios::out | std::ios::binary, binary);
	}
	catch (const std::runtime_error& e) {
		SignalExceptionInterop(e.what(), 0, NULL);
	}
	
	try {
		// Write header line
		dbid_file.write(reinterpret_cast<char *>(&m_dbidRecordCount), sizeof(UINT32));

		// Print data
		for_each (vMeasurements->begin(), vMeasurements->end(),
			[&dbid_file](const dnaMsrPtr& msr) {			
				msr->SerialiseDatabaseMap(&dbid_file);
		});

		dbid_file.close();
	}
	catch (const std::ifstream::failure& f) {
		SignalExceptionInterop(static_cast<std::string>(f.what()), 0, NULL);
	}
	catch (const XMLInteropException& e)  {
		SignalExceptionInterop(static_cast<std::string>(e.what()), 0, NULL);
	}
}
	

void dna_import::PrintMeasurementsToStations(std::string& m2s_file, MsrTally* parsemsrTally,
	std::string& bst_file, std::string& bms_file, std::string& aml_file, pvASLPtr vAssocStnList)
{
	AmlFile aml;
	vmsrtally v_stnmsrTally;
	v_aml_pair vAssocMsrList;

	vstn_t bstBinaryRecords;
	vmsr_t bmsBinaryRecords;
	
	try {
		// Load binary stations data.  Throws runtime_error on failure.
		BstFile bst;
		bst.LoadFile(bst_file, &bstBinaryRecords, bst_meta_);

		// Load binary measurements data.  Throws runtime_error on failure.
		BmsFile bms;
		bms.LoadFile(bms_file, &bmsBinaryRecords, bms_meta_);

		// Load aml file.  Throws runtime_error on failure.
		aml.load_aml_file(aml_file, &vAssocMsrList, &bmsBinaryRecords);
	}
	catch (const std::runtime_error& e) {
		SignalExceptionInterop(e.what(), 0, NULL);
	}

	// Create stns to msrs tally from AML and ASL
	aml.create_msr_to_stn_tally(vAssocStnList, vAssocMsrList, v_stnmsrTally, bmsBinaryRecords);
	
	vUINT32 vStationList(bstBinaryRecords.size());
	// initialise vector with 0,1,2,...,n-2,n-1,n
	initialiseIncrementingIntegerVector<UINT32>(vStationList, static_cast<UINT32>(bstBinaryRecords.size()));

	// Print measurement to station summary, sort stations as required
	switch (projectSettings_.o._sort_msr_to_stn)
	{
	case name_stn_sort_ui:
	{
		CompareStnNameOrder<station_t, UINT32> stnnameCompareFunc(&bstBinaryRecords);
		std::sort(vStationList.begin(), vStationList.end(), stnnameCompareFunc);
		break;
	}
	case meas_stn_sort_ui:
	{
		CompareMeasCount2<ASLPtr, UINT32> msrcountCompareFunc(vAssocStnList);
		std::sort(vStationList.begin(), vStationList.end(), msrcountCompareFunc);
		break;
	}
	case saem_stn_sort_ui:
	{
		CompareMeasCount2Desc<ASLPtr, UINT32> msrcountCompareFunc(vAssocStnList);
		std::sort(vStationList.begin(), vStationList.end(), msrcountCompareFunc);
		break;
	}
	case orig_stn_sort_ui:
	default:
	{
		CompareStnFileOrder<station_t, UINT32> stnorderCompareFunc(&bstBinaryRecords);
		std::sort(vStationList.begin(), vStationList.end(), stnorderCompareFunc);
	}
	break;
	}

	std::ofstream m2s_stream;

	try {
		// Create m2s file.  Throws runtime_error on failure.
		file_opener(m2s_stream, m2s_file);

		// Print formatted header
		print_file_header(m2s_stream, "DYNADJUST MEASUREMENT TO STATION OUTPUT FILE");

		m2s_stream << std::setw(PRINT_VAR_PAD) << std::left << "File name:" << safe_absolute_path(m2s_file) << std::endl << std::endl;

		m2s_stream << std::setw(PRINT_VAR_PAD) << std::left << "Associated measurement file: " << safe_absolute_path(aml_file) << std::endl;
		m2s_stream << std::setw(PRINT_VAR_PAD) << std::left << "Stations file:" << safe_absolute_path(bst_file) << std::endl;
		m2s_stream << std::setw(PRINT_VAR_PAD) << std::left << "Measurements file:" << safe_absolute_path(bms_file) << std::endl;

		// Print station count
		m2s_stream << std::endl;
		m2s_stream << std::setw(PRINT_VAR_PAD) << std::left << "No. stations:" << vStationList.size() << std::endl;

		m2s_stream << OUTPUTLINE << std::endl << std::endl;

		aml.write_msr_to_stn(m2s_stream, &bstBinaryRecords, &vStationList, v_stnmsrTally, parsemsrTally);

		m2s_stream.close();
	}
	catch (const std::runtime_error& e) {
		SignalExceptionInterop(e.what(), 0, NULL);
	}
}
	

// First item in the file is a UINT32 value - the number of records in the file
// All records are of type UINT32
void dna_import::SerialiseAml(const std::string& aml_filename, pvUINT32 vAML)
{
	try {
		// write the aml file.
		AmlFile aml;
		aml.write_aml_file(aml_filename, vAML);
	}
	catch (const std::runtime_error& e) {
		SignalExceptionInterop(e.what(), 0, NULL);
	}
}
	

void dna_import::SerialiseAmlTextFile(const std::string& bms_filename, const std::string& aml_filename, pvUINT32 vAML, pvASLPtr vAssocStnList, vdnaStnPtr* vStations)
{
	try {
		// write the aml file as raw text.
		AmlFile aml;
		aml.write_aml_file_txt(bms_filename, aml_filename + ".txt", vAML, vAssocStnList, vStations);
	}
	catch (const std::runtime_error& e) {
		SignalExceptionInterop(e.what(), 0, NULL);
	}
}
	

// First item in the file is a UINT32 value - the number of records in the file
// All records are of type ASLPtr
void dna_import::SerialiseAsl(const std::string& filename, pvASLPtr vAssocStnList)
{
	try {
		// write the asl file.
		AslFile asl(filename);
		asl.Write(*vAssocStnList);
	}
	catch (const std::runtime_error& e) {
		SignalExceptionInterop(e.what(), 0, NULL);
	}
}
	

// First item in the file is a UINT32 value - the number of records in the file
// All records are of type ASLPtr
void dna_import::SerialiseAslTextFile(const std::string& filename, pvASLPtr vAssocStnList, vdnaStnPtr* vStations)
{
	try {
		// write the ASCII version of the asl file.
		AslFile asl(filename + ".txt");
		asl.WriteText(*vAssocStnList, *vStations);
	}
	catch (const std::runtime_error& e) {
		SignalExceptionInterop(e.what(), 0, NULL);
	}
}
	

UINT32 dna_import::RemoveDuplicateStations(vdnaStnPtr* vStations, 
			pvstring vduplicateStations, pv_stringstring_doubledouble_pair vnearbyStations)
{
	size_t old_stn_count(vStations->size());
	// A prior sort on name is essential, since the criteria for removing duplicates 
	// is based upon two successive station entries in an ordered vector being equal

	_it_vdnastnptr _it_stn_newend;

	// sort vStations on name and if equal, sort on file order
	std::sort(vStations->begin(), vStations->end(), CompareStationName<dnaStnPtr>());
	
	vduplicateStations->clear();
	EqualStationNameSaveDuplicates<dnaStnPtr, std::string> duplicateStnCompareFunc(vduplicateStations);

	_it_stn_newend = unique(vStations->begin(), vStations->end(), duplicateStnCompareFunc);
	if (_it_stn_newend != vStations->end())
		vStations->resize(_it_stn_newend - vStations->begin());

	if (!projectSettings_.i.search_nearby_stn)
		return static_cast<UINT32>(old_stn_count - vStations->size());

	// search nearby stations using a radial search	
	//
	// begin by sorting on latitude
	std::sort(vStations->begin(), vStations->end(), CompareLatitude<dnaStnPtr>());

	vnearbyStations->clear();
	NearbyStation_LowAcc<dnaStnPtr, stringstring_doubledouble_pair, double> 
		nearbyStnLCompareFunc(projectSettings_.i.search_stn_radius, vnearbyStations);
	NearbyStation_HighAcc<dnaStnPtr, double, stringstring_doubledouble_pair, CDnaEllipsoid> 
		nearbyStnHCompareFunc(projectSettings_.i.search_stn_radius, vnearbyStations, 
			datum_.GetEllipsoid());

	// Find all occurrences of nearby stations, using a function appropriate for the search distance
	if (projectSettings_.i.search_stn_radius < 10.0)
		copy_if_all_occurrences(vStations->begin(), vStations->end(), nearbyStnLCompareFunc);
	else
		copy_if_all_occurrences(vStations->begin(), vStations->end(), nearbyStnHCompareFunc);
	
	// sort station pairs by name
	std::sort(vnearbyStations->begin(), vnearbyStations->end(), CompareStationPairs<stringstring_doubledouble_pair>());
	
	// sort by station name
	std::sort(vStations->begin(), vStations->end());

	return static_cast<UINT32>(vnearbyStations->size());
}
	

UINT32 dna_import::FindSimilarMeasurements(vdnaMsrPtr* vMeasurements, vdnaMsrPtr* vSimilarMeasurements)
{
	// sort measurements list by Type then by First station
	std::sort(vMeasurements->begin(), vMeasurements->end(), CompareMsr<CDnaMeasurement>());

	//bool similar;
	int similar_msrs_found(0);
	vSimilarMeasurements->clear();
	
	_it_vdnamsrptr _it_msr(vMeasurements->begin() + 1);
	_it_vdnamsrptr _it_msrprev(vMeasurements->begin());
	
	for (; _it_msr!=vMeasurements->end(); ++_it_msr, ++_it_msrprev)
	{
		// different type?
		if (_it_msr->get()->GetTypeC() != _it_msrprev->get()->GetTypeC())
			continue;		

		//similar = false;

		switch (_it_msr->get()->GetTypeC())
		{
		case 'A':	// Horizontal angles
			if (*(static_cast<const CDnaAngle*>(&(*_it_msr->get()))) ==
				*(static_cast<const CDnaAngle*>(&(*_it_msrprev->get()))))
				++similar_msrs_found;
			else
				continue;
			break;
		
		case 'C':	// Chord dist
		case 'E':	// Ellipsoid arc
		case 'M':	// MSL arc
		case 'S':	// Slope distance
			if (*(static_cast<const CDnaDistance*>(&(*_it_msr->get()))) == 
				*(static_cast<const CDnaDistance*>(&(*_it_msrprev->get()))))
				++similar_msrs_found;
			else
				continue;
			break;

		case 'D':	// Direction set
			if (*(static_cast<const CDnaDirectionSet*>(&(*_it_msr->get()))) ==
				*(static_cast<const CDnaDirectionSet*>(&(*_it_msrprev->get()))))
			{
				similar_msrs_found += static_cast<UINT32>(_it_msr->get()->GetDirections_ptr()->size());
				// Ignore the directions
				if (projectSettings_.i.ignore_similar_msr)
				{
					for_each(_it_msr->get()->GetDirections_ptr()->begin(),
						_it_msr->get()->GetDirections_ptr()->end(),
						[this] (CDnaDirection& dir) {
							dir.SetIgnore(true);
					});
				}
			}
			else
				continue;
			break;

		case 'B':	// Geodetic azimuth
		case 'K':	// Astronomic azimuth
			if (*(static_cast<const CDnaAzimuth*>(&(*_it_msr->get()))) ==
				*(static_cast<const CDnaAzimuth*>(&(*_it_msrprev->get()))))
				++similar_msrs_found;
			else
				continue;
			break;

		case 'I':	// Astronomic latitude
		case 'J':	// Astronomic longitude
		case 'P':	// Geodetic latitude
		case 'Q':	// Geodetic longitude
			if (*(static_cast<const CDnaCoordinate*>(&(*_it_msr->get()))) ==
				*(static_cast<const CDnaCoordinate*>(&(*_it_msrprev->get()))))
				++similar_msrs_found;
			else
				continue;
			break;

		case 'H':	// Orthometric height
		case 'R':	// Ellipsoidal height
			if (*(static_cast<const CDnaHeight*>(&(*_it_msr->get()))) ==
				*(static_cast<const CDnaHeight*>(&(*_it_msrprev->get()))))
				++similar_msrs_found;
			else
				continue;
			break;
				
		case 'L':	// Level difference
			if (*(static_cast<const CDnaHeightDifference*>(&(*_it_msr->get()))) ==
				*(static_cast<const CDnaHeightDifference*>(&(*_it_msrprev->get()))))
				++similar_msrs_found;
			else
				continue;
			break;
				
		case 'V':	// Zenith distance
		case 'Z':	// Vertical angle
			if (*(static_cast<const CDnaDirection*>(&(*_it_msr->get()))) ==
				*(static_cast<const CDnaDirection*>(&(*_it_msrprev->get()))))
				++similar_msrs_found;
			else
				continue;
			break;
				
		case 'G':	// GPS Baseline (treat as single-baseline cluster)
		case 'X':	// GPS Baseline cluster
			if (*(static_cast<const CDnaGpsBaselineCluster*>(&(*_it_msr->get()))) ==
				*(static_cast<const CDnaGpsBaselineCluster*>(&(*_it_msrprev->get()))))
			{
				similar_msrs_found += static_cast<UINT32>(_it_msr->get()->GetBaselines_ptr()->size() * 3);
				// Ignore the baselines
				if (projectSettings_.i.ignore_similar_msr)
				{
					for_each(_it_msr->get()->GetBaselines_ptr()->begin(),
						_it_msr->get()->GetBaselines_ptr()->end(),
						[this] (CDnaGpsBaseline& bsl) {
							bsl.SetIgnore(true);
					});
				}
			}
			else
				continue;
			break;
				
		case 'Y':	// GPS point cluster			
			if (*(static_cast<const CDnaGpsPointCluster*>(&(*_it_msr->get()))) ==
				*(static_cast<const CDnaGpsPointCluster*>(&(*_it_msrprev->get()))))
			{
				similar_msrs_found += static_cast<UINT32>(_it_msr->get()->GetPoints_ptr()->size() * 3);
				// Ignore the points
				if (projectSettings_.i.ignore_similar_msr)
				{
					for_each(_it_msr->get()->GetPoints_ptr()->begin(),
						_it_msr->get()->GetPoints_ptr()->end(),
						[this] (CDnaGpsPoint& pnt) {
							pnt.SetIgnore(true);
					});
				}
			}
			else
				continue;
			break;
		default:
			std::stringstream ss;
			ss << "FindSimilarMeasurements(): Unknown measurement type:  " << _it_msr->get()->GetTypeC() << std::endl;
			throw XMLInteropException(ss.str(), 0); 
		}
		
		vSimilarMeasurements->push_back(*_it_msr);
		if (projectSettings_.i.ignore_similar_msr)
			_it_msr->get()->SetIgnore(true);
	}

	return similar_msrs_found;
}
	

UINT32 dna_import::FindSimilarGXMeasurements(vdnaMsrPtr* vMeasurements, vdnaMsrPtr* vSimilarMeasurements)
{
	// copy all measurements
	vdnaMsrPtr vMeasurementsG(*vMeasurements);

	// Strip non-GX measurements
	std::string msrTypes("GX");
	CompareNonMeasType<CDnaMeasurement> meastypeCompareFuncGX(msrTypes);
	erase_if(vMeasurementsG, meastypeCompareFuncGX);

	msrTypes = "X";
	// Copy X measurements
	vdnaMsrPtr vMeasurementsX(vMeasurementsG);
	CompareNonMeasType<CDnaMeasurement> meastypeCompareFuncX(msrTypes);
	erase_if(vMeasurementsX, meastypeCompareFuncX);

	msrTypes = "G";
	meastypeCompareFuncGX.SetComparand(msrTypes);
	erase_if(vMeasurementsG, meastypeCompareFuncGX);

	// Don't think the following sorts are needed
	// sort measurements list by Type then by First station
	//sort(vMeasurementsG.begin(), vMeasurementsG.end(), CompareMsr<CDnaMeasurement>());
	//sort(vMeasurementsX.begin(), vMeasurementsX.end(), CompareMsr<CDnaMeasurement>());

	bool similarG(false);
	int similar_msrs_found(0);
	vSimilarMeasurements->clear();

	_it_vdnamsrptr _it_msrG(vMeasurementsG.begin());
	vstring stationsX;
	std::vector<CDnaGpsBaseline> *vgpsBslsX, *vgpsBslsG;
	std::string epochX, epochG;
	boost::gregorian::date dateObsX, dateObsG;
	std::stringstream strdiffDays;
	size_t diffDays;
	vUINT32 cluster_ids;

	// For each X Cluster, get the list of stations
	for (_it_vdnamsrptr _it_msrX(vMeasurementsX.begin()); _it_msrX!=vMeasurementsX.end(); ++_it_msrX)
	{
		vgpsBslsX = _it_msrX->get()->GetBaselines_ptr();
		GetGXMsrStations<std::string>(vgpsBslsX, stationsX);

		std::sort(stationsX.begin(), stationsX.end());

		similarG = false;

		// Check whether both stations in the G baseline exist
		// in the list of X Cluster stations
		for (_it_msrG=vMeasurementsG.begin(); _it_msrG!=vMeasurementsG.end(); ++_it_msrG)
		{
			vgpsBslsG = _it_msrG->get()->GetBaselines_ptr();
			
			// Is the first G baseline station in the cluster?
			if (!binary_search(stationsX.begin(), stationsX.end(), vgpsBslsG->at(0).GetFirst()))
				// No, then this isn't a "duplicate"
				continue;

			// Is the second G baseline station in the cluster?
			if (!binary_search(stationsX.begin(), stationsX.end(), vgpsBslsG->at(0).GetTarget()))
				// No, then this isn't a "duplicate"
				continue;

			// Okay, both stations are in the cluster... check if this is a "duplicate"
			// by examining epoch.  Here, the assumption is - if the observation epochs of the 
			// X and G measurements are the same, then the X and G measurements have come from
			// the same source data and are therefore duplicates
			epochX = _it_msrX->get()->GetEpoch();
			epochG = _it_msrG->get()->GetEpoch();

			dateObsX = dateFromString<boost::gregorian::date>(epochX);
			dateObsG = dateFromString<boost::gregorian::date>(epochG);

			boost::gregorian::days dateDifference = dateObsX - dateObsG;
			
			diffDays = abs(dateDifference.days());
			//strdiffDays << "Elapsed days: " << dateDifference;
			//TRACE("%s\n", strdiffDays.str().c_str());

			if (diffDays < 5)
			{
				++similar_msrs_found;
				similarG = true;
				vSimilarMeasurements->push_back(*_it_msrG);
				if (projectSettings_.i.ignore_similar_msr)
					cluster_ids.push_back(_it_msrG->get()->GetClusterID());
			}
		}

		// Were similar G baselines found for this cluster?
		if (similarG)
			vSimilarMeasurements->push_back(*_it_msrX);

	}

	// If required, ignore the measurements
	if (projectSettings_.i.ignore_similar_msr)
	{
		std::sort(cluster_ids.begin(), cluster_ids.end());
		for (_it_msrG=vMeasurements->begin(); _it_msrG != vMeasurements->end(); _it_msrG++)
			IgnoreGXMeasurements(_it_msrG->get(), cluster_ids.begin(), cluster_ids.end());
	}

	return similar_msrs_found;
}

void dna_import::SortandMapStations(vdnaStnPtr* vStations)
{
	size_t stnCount(vStations->size());
	vStnsMap_sortName_.clear();
	vStnsMap_sortName_.reserve(stnCount);

	UINT32 stnIndex(0);

	// sort on station name (by string, not int!!!)
	// Note that the sort order after this will be the final order
	std::sort(vStations->begin(), vStations->end());

	// Create the Station-Name / ID map
	string_uint32_pair stnID;
	for (stnIndex = 0; stnIndex < stnCount; stnIndex++)
	{
		stnID.first = vStations->at(stnIndex)->GetName();
		stnID.second = stnIndex;
		vStnsMap_sortName_.push_back(stnID);
		vStations->at(stnIndex)->SetnameOrder(stnIndex);
	}

	// sort on station name (i.e. first of the pair)
	std::sort(vStnsMap_sortName_.begin(), vStnsMap_sortName_.end(), StationNameIDCompareName());

	if (vStnsMap_sortName_.size() < stnCount)
		throw XMLInteropException("SortandMapStations(): Could not allocate sufficient memory for the Station map.", 0);
}


void dna_import::FullSortandMapStations(vdnaStnPtr* vStations, pv_string_uint32_pair vStnsMap_sortName)
{
	// Sort the stations list and develop the station map
	SortandMapStations(vStations);
	*vStnsMap_sortName = vStnsMap_sortName_;
}
	
void dna_import::SortStationsForExport(vdnaStnPtr* vStations)
{
	// Sort on original file order
	std::sort(vStations->begin(), vStations->end(), CompareStnFileOrder_CDnaStn<CDnaStation>());
}
		

void dna_import::ReduceStations(vdnaStnPtr* vStations, const CDnaProjection& projection)
{
	// reduce cartesian or projection coordinates to geographic (radians)
	for_each(vStations->begin(), vStations->end(),
		[this, &projection](dnaStnPtr stn){
			stn->ReduceStations_LLH(datum_.GetEllipsoidRef(), &projection);
	});

	// Try to reduce stations concurrently...
	//parallel_for_each(vStations->begin(), vStations->end(),
	//	[this, &projection](dnaStnPtr stn){
	//		stn->ReduceStations_LLH(datum_.GetEllipsoidRef(), &projection);
	//});
}	
	

void dna_import::RenameStations(vdnaStnPtr* vStations, vdnaMsrPtr* vMeasurements, project_settings* p)
{
	if (!std::filesystem::exists(p->i.stn_renamingfile))
	{
		std::string s("The station renaming file cannot be found:\n");
		s.append("    ").append(p->i.stn_renamingfile);
		throw XMLInteropException(s, 0);
	}

	// station renaming file
	//
	// <--- (20 chars) ---><--- (20 chars) ---><--- (20 chars) ---><--- (20 chars) --->
	// PREFERRED NAME      ALIAS               ALIAS               ALIAS               ...
	//
	// First column is the preferred name
	// Remaining columns are aliases (which may include the preferred name)
	// To increase performance, don't provide records where an alias equals the preferred name
	v_string_vstring_pair stationNames;
	dna_io_dna dna;
	dna.read_ren_file(p->i.stn_renamingfile, &stationNames);
	std::sort(stationNames.begin(), stationNames.end());

	// rename stations in stations vector
	for_each(vStations->begin(), vStations->end(), 
		[&stationNames](dnaStnPtr& stn) {
			for (_it_string_vstring_pair it = stationNames.begin();
				it != stationNames.end();
				++it)
			{
				// search the aliases for this station's name
				if (binary_search(it->second.begin(), it->second.end(), stn->GetName()))
				{
					// This name is one of the aliases, so replace it 
					// with the preferred name and break out (no need to continue)
					stn->SetName(it->first);
					break;
				}
			}
	});

	_it_vdnamsrptr _it_msr(vMeasurements->begin());
	std::vector<CDnaDirection>* vdirns;
	std::vector<CDnaGpsBaseline>* vgpsBsls;
	std::vector<CDnaGpsPoint>* vgpsPnts;

	// rename stations in each measurement
	for (_it_msr=vMeasurements->begin(); _it_msr != vMeasurements->end(); _it_msr++)
	{
		// 1. Handle nested type measurements (D, G, X, Y) separately
		switch (_it_msr->get()->GetTypeC())
		{
		case 'G':	// GPS Baseline (treat as single-baseline cluster)
		case 'X':	// GPS Baseline cluster
			vgpsBsls = _it_msr->get()->GetBaselines_ptr();
			RenameStationsBsl(vgpsBsls, stationNames);			
			continue;
		case 'Y':	// GPS point cluster
			vgpsPnts = _it_msr->get()->GetPoints_ptr();
			RenameStationsPnt(vgpsPnts, stationNames);
			continue;
		case 'D':	// GPS point cluster
			// Rename stations in first direction
			RenameStationsMsr(_it_msr->get(), stationNames.begin(), stationNames.end());
			
			// Rename stations in all other directions
			vdirns = _it_msr->get()->GetDirections_ptr();
			RenameStationsDir(vdirns, stationNames);
			continue;
		}

		RenameStationsMsr(_it_msr->get(), stationNames.begin(), stationNames.end());
	}
	
}	

void dna_import::RenameStationsBsl(std::vector<CDnaGpsBaseline>* vGpsBaselines, 
	v_string_vstring_pair& stnRenaming)
{
	for_each(vGpsBaselines->begin(), vGpsBaselines->end(),
		[this, &stnRenaming] (CDnaGpsBaseline& bsl) {
			RenameStationsMsr(&bsl, stnRenaming.begin(), stnRenaming.end());
	});
}
	

void dna_import::RenameStationsPnt(std::vector<CDnaGpsPoint>* vGpsPoints, 
	v_string_vstring_pair& stnRenaming)
{
	for_each(vGpsPoints->begin(), vGpsPoints->end(),
		[this, &stnRenaming] (CDnaGpsPoint& pnt) {
			RenameStationsMsr(&pnt, stnRenaming.begin(), stnRenaming.end());
	});
}
	

void dna_import::RenameStationsDir(std::vector<CDnaDirection>* vDirections, 
	v_string_vstring_pair& stnRenaming)
{
	for_each(vDirections->begin(), vDirections->end(),
		[this, &stnRenaming] (CDnaDirection& dir) {
			RenameStationsMsr(&dir, stnRenaming.begin(), stnRenaming.end());
	});
}

void dna_import::ApplyGNSSMsrScalar(std::vector<CDnaGpsBaseline>::iterator& _it_bsl, vscl_t& bslScalars)
{
	vscl_t::iterator _it_scalar;
	CompareScalarStations< scalar_t, std::string> scalarComparisonFunc(_it_bsl->GetFirst(), _it_bsl->GetTarget());

	if ((_it_scalar = find_if(bslScalars.begin(), bslScalars.end(), scalarComparisonFunc)) != bslScalars.end())
	{		
		_it_bsl->SetVscale(_it_bsl->GetVscale() * _it_scalar->v_scale);
		_it_bsl->SetPscale(_it_bsl->GetPscale() * _it_scalar->p_scale);
		_it_bsl->SetLscale(_it_bsl->GetLscale() * _it_scalar->l_scale);
		_it_bsl->SetHscale(_it_bsl->GetHscale() * _it_scalar->h_scale);
	}
}
	
	
void dna_import::EditGNSSMsrScalars(vdnaMsrPtr* vMeasurements, project_settings* p)
{
	// Set user-defined GNSS measurement scalars
	_it_vdnamsrptr _it_msr(vMeasurements->begin());

	// Scalar file only makes sense for individual GNSS measurements, not
	// GNSS baseline or point clusters!!!
	bool applyScalarFile(false);
	bool applyVScale(fabs(p->i.vscale - 1.0) > PRECISION_1E5);
	bool applyPScale(fabs(p->i.pscale - 1.0) > PRECISION_1E5);
	bool applyLScale(fabs(p->i.lscale - 1.0) > PRECISION_1E5);
	bool applyHScale(fabs(p->i.hscale - 1.0) > PRECISION_1E5);

	vscl_t bslScalars;

	if (!p->i.scalar_file.empty())
	{
		if (!std::filesystem::exists(p->i.scalar_file))
		{
			std::string s("The GNSS scalar file cannot be found:\n");
			s.append("    ").append(p->i.scalar_file);
			throw XMLInteropException(s, 0);
		}

		// load up scalars
		applyScalarFile = true;
		dna_io_scalar scalar;
		scalar.load_scalar_file(p->i.scalar_file, &bslScalars);
	}

	std::sort(bslScalars.begin(), bslScalars.end(), CompareScalars<scalar_t>());

	std::vector<CDnaGpsBaseline>* vgpsBsls;
	std::vector<CDnaGpsBaseline>::iterator _it_bsl;
	std::vector<CDnaGpsPoint>* vgpsPnts;
	std::vector<CDnaGpsPoint>::iterator _it_pnt;

	// set scalars first
	for (_it_msr=vMeasurements->begin(); _it_msr != vMeasurements->end(); _it_msr++)
	{		
		switch (_it_msr->get()->GetTypeC())
		{
		case 'G':	// GPS Baseline (treat as single-baseline cluster)
		case 'X':	// GPS Baseline cluster
	
			if (applyVScale)
				(static_cast<CDnaGpsBaselineCluster*>(&(*_it_msr->get()))
				)->SetVscale(_it_msr->get()->GetVscale() * p->i.vscale);
			if (applyPScale)
				(static_cast<CDnaGpsBaselineCluster*>(&(*_it_msr->get()))
				)->SetPscale(_it_msr->get()->GetPscale() * p->i.pscale);
			if (applyLScale)
				(static_cast<CDnaGpsBaselineCluster*>(&(*_it_msr->get()))
				)->SetLscale(_it_msr->get()->GetLscale() * p->i.lscale);
			if (applyHScale)
				(static_cast<CDnaGpsBaselineCluster*>(&(*_it_msr->get()))
				)->SetHscale(_it_msr->get()->GetVscale() * p->i.hscale);

			vgpsBsls = _it_msr->get()->GetBaselines_ptr();
			for (_it_bsl=vgpsBsls->begin(); _it_bsl!=vgpsBsls->end(); ++_it_bsl)
			{

				if (_it_msr->get()->GetTypeC() == 'G')
				{
					// Since the only way to associate scaling with measurements in DNA and
					// DynaML files is to match the from and to stations (unlike GeodeticXML
					// which allows for a unique ID to be used to identify each measurement),
					// the scalars file is only applied to baselines at this stage.
					if (applyScalarFile)
					{
						// Find the stations in the list contained in the scalar file.
						// If the first and second stations match, then apply.
						ApplyGNSSMsrScalar(_it_bsl, bslScalars);						
					}
				}
				if (applyVScale)
					_it_bsl->SetVscale(_it_bsl->GetVscale() * p->i.vscale);
				if (applyPScale)
					_it_bsl->SetPscale(_it_bsl->GetPscale() * p->i.pscale);
				if (applyLScale)
					_it_bsl->SetLscale(_it_bsl->GetLscale() * p->i.lscale);
				if (applyHScale)
					_it_bsl->SetHscale(_it_bsl->GetVscale() * p->i.hscale);
			}
			
			break;
		case 'Y':	// GPS point cluster
			if (applyVScale)
				(static_cast<CDnaGpsPointCluster*>(&(*_it_msr->get()))
				)->SetVscale(_it_msr->get()->GetVscale() * p->i.vscale);
			if (applyPScale)
				(static_cast<CDnaGpsPointCluster*>(&(*_it_msr->get()))
				)->SetPscale(_it_msr->get()->GetPscale() * p->i.pscale);
			if (applyLScale)
				(static_cast<CDnaGpsPointCluster*>(&(*_it_msr->get()))
				)->SetLscale(_it_msr->get()->GetLscale() * p->i.lscale);
			if (applyHScale)
				(static_cast<CDnaGpsPointCluster*>(&(*_it_msr->get()))
				)->SetHscale(_it_msr->get()->GetVscale() * p->i.hscale);

			vgpsPnts = _it_msr->get()->GetPoints_ptr();
			for (_it_pnt=vgpsPnts->begin(); _it_pnt!=vgpsPnts->end(); ++_it_pnt)
			{
				if (applyVScale)
					_it_pnt->SetVscale(_it_pnt->GetVscale() * p->i.vscale);
				if (applyPScale)
					_it_pnt->SetPscale(_it_pnt->GetPscale() * p->i.pscale);
				if (applyLScale)
					_it_pnt->SetLscale(_it_pnt->GetLscale() * p->i.lscale);
				if (applyHScale)
					_it_pnt->SetHscale(_it_pnt->GetHscale() * p->i.hscale);
			}
			
			break;
		}
	}
}
	

void dna_import::SerialiseMap(const std::string& stnmap_file)
{
	try {
		// write the aml file.
		dynadjust::iostreams::MapFile map;
		map.WriteFile(stnmap_file, &vStnsMap_sortName_);
	}
	catch (const std::runtime_error& e) {
		SignalExceptionInterop(e.what(), 0, NULL);
	}
}
	

void dna_import::SerialiseMapTextFile(const std::string& stnmap_file)
{
	try {
		// write the aml file as raw text.
		dynadjust::iostreams::MapFile map;
		map.WriteTextFile(stnmap_file + ".txt", &vStnsMap_sortName_);
	}
	catch (const std::runtime_error& e) {
		SignalExceptionInterop(e.what(), 0, NULL);
	}
}


void dna_import::SerialiseDiscontTextFile(const std::string& discont_file)
{
	std::ofstream discont_outfile;
	
	try {
		std::filesystem::path discontFile(discont_file);
		std::string outfileName = discontFile.filename().string();
		outfileName.append(".discont");
		// Open discontinuity output file.  Throws runtime_error on failure.
		file_opener(discont_outfile, outfileName,
			std::ios::out, ascii);
	}
	catch (const std::runtime_error& e) {
		SignalExceptionInterop(e.what(), 0, NULL);
	}

	try {
		discont_outfile <<
			std::setw(14) << std::left << "file_index" <<
			std::setw(14) << std::left << "solution_id" <<
			std::setw(STATION) << std::left << "site_name" <<
			std::setw(14) << "date_start" <<
			std::setw(14) << "date_end" <<
			std::setw(24) << "discontinuity_exists" << std::endl;

		// Print discontinuity file.  Throws runtime_error on failure.
		for_each(stn_discontinuities_.begin(), stn_discontinuities_.end(),
			[this, &discont_outfile](discontinuity_tuple& discont) {			// use lambda expression
			discont_outfile <<
				std::setw(14) << std::left << discont.file_index <<
				std::setw(14) << std::left << discont.solution_id <<
				std::setw(STATION) << std::left << discont.site_name <<
				std::setw(14) << stringFromDate(discont.date_start) <<
				std::setw(14) << stringFromDate(discont.date_end) <<
				std::setw(24) << discont.discontinuity_exists << std::endl;
		}
		);
	}
	catch (const std::runtime_error& e) {
		SignalExceptionParse(e.what(), 0);
	}

	discont_outfile.close();
}
	

void dna_import::SimulateMSR(vdnaStnPtr* vStations, vdnaMsrPtr* vMeasurements, 
	const std::string& msrfilename, 
	const project_settings& p)
{
	CDnaProjection projection(UTM);

	try {
		// reduce geographic or projection coordinates to cartesian
		for_each(vStations->begin(), vStations->end(),
			[this, &projection](dnaStnPtr stn){
				stn->ReduceStations_XYZ(datum_.GetEllipsoidRef(), &projection);
		});
	}
	catch (...) {
		std::stringstream ss;
		ss << "SimulateMSR(): An unexpected error occurred when reducing station coordinates.";
		SignalExceptionInterop(ss.str(), 0, NULL);
	}

	// Simulate measurements
	_it_vdnamsrptr _it_msr;

	try {
		for_each(vMeasurements->begin(), vMeasurements->end(),
			[this, &vStations](dnaMsrPtr _it_msr){
				_it_msr->SimulateMsr(vStations, datum_.GetEllipsoidRef());
		});
	}
	catch (...) {
		std::stringstream ss;
		ss << "SimulateMSR(): An unexpected error occurred when simulating measurements.";
		SignalExceptionInterop(ss.str(), 0, NULL);
	}


	try {
		// reduce cartesian coordinates back to geographic
		for_each(vStations->begin(), vStations->end(),
			[this, &projection](dnaStnPtr stn){
				stn->ReduceStations_LLH(datum_.GetEllipsoidRef(), &projection);
		});
	}
	catch (...) {
		std::stringstream ss;
		ss << "SimulateMSR(): An unexpected error occurred when reducing station coordinates.";
		SignalExceptionInterop(ss.str(), 0, NULL);
	}

	determineDNAMSRFieldParameters<UINT16>("3.01", dml_, dmw_);

	// Print simulated measurements
	std::ofstream dna_msr_file;
	try {
		// Create DynAdjust MSR file. 
		file_opener(dna_msr_file, msrfilename);
	}
	catch (const std::runtime_error& e) {
		SignalExceptionInterop(e.what(), 0, NULL);
	}

	try {	
		// Write version line
		dna_header(dna_msr_file, "3.01", "MSR", datum_.GetName(), datum_.GetEpoch_s(), vMeasurements->size());

		// Write header line
		dna_comment(dna_msr_file, p.g.network_name + " measurement file.");

		// print measurements
		for_each(vMeasurements->begin(), vMeasurements->end(),
			[this, &dna_msr_file](dnaMsrPtr _it_msr){
				_it_msr->WriteDNAMsr(&dna_msr_file, dmw_, dml_);
		});

		dna_msr_file.close();

	}
	catch (const std::ifstream::failure& f) {
		SignalExceptionInterop(static_cast<std::string>(f.what()), 0, "o", &dna_msr_file);
	}
	catch (const XMLInteropException& e)  {
		SignalExceptionInterop(static_cast<std::string>(e.what()), 0, "o", &dna_msr_file);
	}
}
	

void dna_import::IgnoreInsufficientMeasurements(vdnaStnPtr* vStations, vdnaMsrPtr* vMeasurements,
	pvstring vPoorlyConstrainedStns)
{
	// This function requires the formation of a station map and the ASL file, which 
	// is not normally formed until much later in the process.

	// 1. Develop station map
	SortandMapStations(vStations);

	vASLPtr vAssocStnList;		// Temporary ASL
	vUINT32 vAssocMsrList;		// Temporary AML

	vstring vunusedStations;
	vUINT32 vignoredMeasurements;
	UINT32 lMapCount;

	vPoorlyConstrainedStns->clear();
	vstring vInsufficientMsrStns;
	UINT32 stnIndex, stnmsr_indexAML;
	bool foundInsufficientMeasurement;

	vdnaMsrPtr vignMsr, vignMsrTotal;
	std::vector<CDnaDirection>* vdirns;

	vignMsrTotal.clear();

	do {

		// 2. Develop ASL
		MapMeasurementStations(vMeasurements, &vAssocStnList, &lMapCount,
			&vunusedStations, &vignoredMeasurements);

		// 3. Develop AML
		CompleteAssociationLists(vMeasurements, &vAssocStnList, &vAssocMsrList, cls_msr);

		// 4. Detect insufficient measurements
		vASLPtr::iterator _it_asl(vAssocStnList.begin()), _it_asl_begin(vAssocStnList.begin());
		
		foundInsufficientMeasurement = false;
		vInsufficientMsrStns.clear();

		// find stations for which there are insufficient measurements
		for (_it_asl = vAssocStnList.begin(); _it_asl != vAssocStnList.end(); _it_asl++)
		{
			// Is there only one measurement to this station?
			if (_it_asl->get()->GetAssocMsrCount() == 1)
			{
				stnIndex = static_cast<UINT32>(std::distance(_it_asl_begin, _it_asl));

				stnmsr_indexAML = vAssocMsrList.at(_it_asl->get()->GetAMLStnIndex());				
				switch (vMeasurements->at(stnmsr_indexAML)->GetTypeC())
				{
				case 'G':	// GPS Baseline (treat as single-baseline cluster)
				case 'X':	// GPS Baseline cluster
				case 'Y':	// GPS point cluster
					// Continue if a GNSS measurement is found, since GNSS is able to 
					// estimate station coordinates in 3D
					continue;
				case 'H':	// Orthometric height
				case 'R':	// Ellipsoidal height
				case 'V':	// Zenith distance
				case 'Z':	// Vertical angle
				case 'L':	// Level difference
					// Continue if:
					// - an ellipsoid/orthometric height or vertical angle/zenith distance is found, and either
					// - the station is constrained in 3D, or
					// - the station is constrained in 2D
					switch (vStations->at(stnIndex)->GetConstraintType())
					{
					case constrained_3D:		// CCC: 3D is constrained, no measurements are needed
					case free_1D:				// CCF: 2D is constrained, one measurement is needed
						// All good
						continue;
					case free_2D:				// FFC: 2D is free, two measurements are needed
					case free_3D:				// FFF: 3D is free, three measurements are needed
					case custom_constraint:		// Unsure: could be anything.
						break;
					}
					// Insufficient vertical measurements
					foundInsufficientMeasurement = true;
					vInsufficientMsrStns.push_back(vStnsMap_sortName_.at(stnIndex).first);
					break;
				default:	
					// All other 2D measurement types, which, alone, do not 
					// sufficiently constrain coordinates in 2D.  So, continue if:
					// - a latitude or longitude or ellipsoid/orthometric height or vertical angle/zenith distance is found, and
					// - the station is held constrained in two dimensions
					switch (vStations->at(stnIndex)->GetConstraintType())
					{
					case constrained_3D:		// CCC: 3D is constrained, no measurements are needed
						// All good.
						continue;
					case free_1D:				// CCF: 2D is constrained, one measurement is needed
					case free_2D:				// FFC: 2D is free, two measurements are needed
					case free_3D:				// FFF: 3D is free, three measurements are needed
					case custom_constraint:		// Unsure: could be anything.
						break;
					}
					// Insufficient horizontal measurements, because:
					// - CCF needs one measurement in vertical (H, R, V, Z, L)
					// - FFC needs two measurements
					// - FFF needs three measurements
					foundInsufficientMeasurement = true;					
					vInsufficientMsrStns.push_back(vStnsMap_sortName_.at(stnIndex).first);
					break;
				}				
			}
		}

		// Are there no stations with insufficient measurements?
		if (!foundInsufficientMeasurement)
			break;

		if (!projectSettings_.g.quiet)
			std::cout << ".";

		// add the newly found insufficient measurements
		vPoorlyConstrainedStns->insert(vPoorlyConstrainedStns->end(),
			vInsufficientMsrStns.begin(), vInsufficientMsrStns.end());

		vignMsr.clear();

		// find measurements connected to these stations and ignore them
		std::sort(vInsufficientMsrStns.begin(), vInsufficientMsrStns.end());
		for_each(vMeasurements->begin(), vMeasurements->end(),
			[&vInsufficientMsrStns, &vMeasurements, 
			&vdirns, &vignMsr](dnaMsrPtr msr) {

				bool insufficientMsrFound(false);

				switch (msr->GetTypeC())
				{
				// GNSS measurements do not qualify in this context.
				// three station measurements
				case 'A':	// Horizontal angle

					if (binary_search(vInsufficientMsrStns.begin(),
						vInsufficientMsrStns.end(), msr->GetTarget2()))
					{
						msr->SetInsufficient(true);
						msr->SetIgnore(true);
						vignMsr.push_back(msr);
						break;
					}
					[[fallthrough]];

				// two station measurements
				case 'B':	// Geodetic azimuth
				case 'C':	// Chord dist
				case 'E':	// Ellipsoid arc
				case 'K':	// Astronomic azimuth
				case 'L':	// Level difference
				case 'M':	// MSL arc
				case 'S':	// Slope distance
				case 'V':	// Zenith distance
				case 'Z':	// Vertical angle

					if (binary_search(vInsufficientMsrStns.begin(),
						vInsufficientMsrStns.end(), msr->GetTarget()))
					{
						msr->SetIgnore(true);
						msr->SetInsufficient(true);
						vignMsr.push_back(msr);
						break;
					}
					[[fallthrough]];

				// one station measurements
				case 'H':	// Orthometric height
				case 'I':	// Astronomic latitude
				case 'J':	// Astronomic longitude
				case 'P':	// Geodetic latitude
				case 'Q':	// Geodetic longitude
				case 'R':	// Ellipsoidal height

					if (binary_search(vInsufficientMsrStns.begin(),
						vInsufficientMsrStns.end(), msr->GetFirst()))
					{
						msr->SetInsufficient(true);
						msr->SetIgnore(true);
						vignMsr.push_back(msr);
						break;
					}
					break;

				case 'D':	// Direction set
				
					// test first station, which if found, ignore the entire measurement.
					if (binary_search(vInsufficientMsrStns.begin(),
						vInsufficientMsrStns.end(), msr->GetFirst()))
					{
						msr->SetInsufficient(true);
						msr->SetIgnore(true);
						insufficientMsrFound = true;
						break;
					}

					// test the target, which if found, ignore the entire measurement
					if (binary_search(vInsufficientMsrStns.begin(),
						vInsufficientMsrStns.end(), msr->GetTarget()))
					{
						msr->SetInsufficient(true);
						msr->SetIgnore(true);
						insufficientMsrFound = true;
						break;
					}

					vdirns = msr.get()->GetDirections_ptr();
					std::vector< CDnaDirection >::iterator _it_dir;
					for (_it_dir = vdirns->begin(); _it_dir != vdirns->end(); ++_it_dir)
					{
						// test the individual directions, which if found, remove the direction from the vector
						if (binary_search(vInsufficientMsrStns.begin(),
							vInsufficientMsrStns.end(), _it_dir->GetTarget()))
						{
							_it_dir->SetInsufficient(true);
							_it_dir->SetIgnore(true);
							insufficientMsrFound = true;
							break;
						}
					}

					if (insufficientMsrFound)
					{
						msr->SetInsufficient(true);
						msr->SetIgnore(true);

						vignMsr.push_back(msr);
					}

					break;
				}

			} // lambda
		);

		// Were any measurements found that contain stations in vInsufficientMsrStns? 
		if (vignMsr.size() > 0)
		{
			// exclude insufficient measurements
			erase_if(vMeasurements, CompareInsufficientMsr<CDnaMeasurement>());
			
			// capture the insufficient measurements
			vignMsrTotal.insert(vignMsrTotal.end(), vignMsr.begin(), vignMsr.end());
		}
		else
			foundInsufficientMeasurement = false;

	} while (foundInsufficientMeasurement);

	
	if (vPoorlyConstrainedStns->size() > 0)
		std::sort(vPoorlyConstrainedStns->begin(), vPoorlyConstrainedStns->end());

	// add the insufficient measurements back to the final vector.
	// Note, these were ignored as part of the process above
	if (vignMsrTotal.size() > 0)
		vMeasurements->insert(vMeasurements->end(), vignMsrTotal.begin(), vignMsrTotal.end());

}
	

void dna_import::MapMeasurementStations(vdnaMsrPtr* vMeasurements, pvASLPtr vAssocStnList, PUINT32 lMapCount, 
	pvstring vUnusedStns, pvUINT32 vIgnoredMsrs)
{
	const size_t mapsize = vStnsMap_sortName_.size();
	if (mapsize < 1)
		throw XMLInteropException("A station map has not been created. Run \"SortStations\" to create a station map.", 0);
	
	g_map_tally.initialise();

	std::ostringstream ss;

	std::vector<CDnaDirection>* vdirns;
	std::vector<CDnaGpsBaseline>* vgpsBsls;
	std::vector<CDnaGpsPoint>* vgpsPnts;

	// Association list initialisation and variables
	std::string at_station_name, to_station_name, to2_station_name;
	UINT32 at_station_index, to_station_index, to2_station_index;
	vAssocStnList->clear();
	vAssocStnList->resize(vStnsMap_sortName_.size());

	//size_t index = 0;
	*lMapCount = 0;
	_it_vdnamsrptr _it_msr(vMeasurements->begin());
	// set map values
	for (_it_msr=vMeasurements->begin(); _it_msr != vMeasurements->end(); _it_msr++)
	{
		//index = _it_msr - vMeasurements->begin();

		// The following is provided so as to cater for the opportunity to 
		// reintroduce an ignored measurement at a later stage, such as via a GUI.
		// The only requirement is that a check be made on stations which become
		// unused as a consequence of ignoring a measurement.
		if (_it_msr->get()->GetIgnore())
			vIgnoredMsrs->push_back(static_cast<UINT32>(std::distance(vMeasurements->begin(), _it_msr)));
		
		// 1. Handle nested type measurements (G, X, Y) separately
		switch (_it_msr->get()->GetTypeC())
		{
		case 'G':	// GPS Baseline (treat as single-baseline cluster)
		case 'X':	// GPS Baseline cluster
			vgpsBsls = _it_msr->get()->GetBaselines_ptr();
			_it_msr->get()->SetStn3Index(static_cast<UINT32>(vgpsBsls->size()));		// set number of GpsBaselines
			MapMeasurementStationsBsl(vgpsBsls, vAssocStnList, lMapCount);			
			continue;
		case 'Y':	// GPS point cluster
			vgpsPnts = _it_msr->get()->GetPoints_ptr();
			_it_msr->get()->SetStn3Index(static_cast<UINT32>(vgpsPnts->size()));		// set number of GpsPoints
			MapMeasurementStationsPnt(vgpsPnts, vAssocStnList, lMapCount);
			continue;
		}

		if (_it_msr->get()->GetFirst().empty())
		{
			//UINT32 msr_no(std::distance(vMeasurements->begin(), _it_msr));
			throw XMLInteropException("Empty \"First\" station name.", 0);
		}

		// 1. All other measurements which have <First>
		at_station_name = _it_msr->get()->GetFirst();		// GetFirst() invokes CDnaMeasurement member function
		
		it_stnmap_range = equal_range(vStnsMap_sortName_.begin(), vStnsMap_sortName_.end(), at_station_name, StationNameIDCompareName());

		if (it_stnmap_range.first == it_stnmap_range.second)
		{
			ss.str("");
			ss << at_station_name << " for a " << _it_msr->get()->GetTypeC() << " measurement";
			
			// Capture cluster ID
			switch (_it_msr->get()->GetTypeC())
			{
			case 'X':
			case 'Y':
			case 'D':
				ss << " (cluster ID " << _it_msr->get()->GetClusterID() << ")";
			}
			
			// Capture 'to' station
			switch (_it_msr->get()->GetTypeC())
			{
			case 'A':	// Horizontal angle
			case 'B':	// Geodetic azimuth
			case 'C':	// Chord dist
			case 'D':	// Direction set
			case 'E':	// Ellipsoid arc
			case 'K':	// Astronomic azimuth
			case 'L':	// Level difference
			case 'M':	// MSL arc
			case 'S':	// Slope distance
			case 'V':	// Zenith distance
			case 'Z':	// Vertical angle
				ss << " to " << _it_msr->get()->GetTarget();
			}

			ss << " was not found in the list of network stations." << std::endl <<
				"  Please ensure that " << at_station_name << 
				" is included in the input station file.";
			throw XMLInteropException(ss.str(), 0);
		}

		// Increment ASL associated measurement count (i.e. CAStationList.assocMsrs_)
		at_station_index = it_stnmap_range.first->second;
		
		_it_msr->get()->SetStn1Index(at_station_index);
		if (!vAssocStnList->at(at_station_index))
			vAssocStnList->at(at_station_index).reset(new CAStationList);

		if (_it_msr->get()->GetTypeC() != 'D')
			// measurement count for the 'at' station will be incremented in MapDirectionStations
			vAssocStnList->at(at_station_index).get()->IncrementMsrCount();

		// Increment single station measurement counters...
		switch (_it_msr->get()->GetTypeC())
		{
		case 'H':	// Orthometric height
			g_map_tally.H ++;
			(*lMapCount)++;
			continue;
		case 'R':	// Ellipsoidal height
			g_map_tally.R ++;
			(*lMapCount)++;
			continue;
		case 'I':	// Astronomic latitude
			g_map_tally.I ++;
			(*lMapCount)++;
			continue;
		case 'J':	// Astronomic longitude
			g_map_tally.J ++;
			(*lMapCount)++;
			continue;
		case 'P':	// Geodetic latitude
			g_map_tally.P ++;
			(*lMapCount)++;
			continue;
		case 'Q':	// Geodetic longitude
			g_map_tally.Q ++;
			(*lMapCount)++;
			continue;
		}

		if (_it_msr->get()->GetTarget().empty())
			throw XMLInteropException("Empty \"Second\" station name", 0);

		// 2. All measurements which have <Target> or <Second>
		to_station_name = _it_msr->get()->GetTarget();		// GetTarget() invokes specialised class (i.e. CDnaDistance, etc) member function

		it_stnmap_range = equal_range(vStnsMap_sortName_.begin(), vStnsMap_sortName_.end(), to_station_name, StationNameIDCompareName());

		if (it_stnmap_range.first == it_stnmap_range.second)
		{
			ss.str("");
			ss << to_station_name << " was not loaded from the stations map.  Please ensure that " << to_station_name << " is included in the list of network stations.";
			throw XMLInteropException(ss.str(), 0);
		}

		// Increment ASL associated measurement count (i.e. CAStationList.assocMsrs_)
		to_station_index = it_stnmap_range.first->second;
		
		_it_msr->get()->SetStn2Index(to_station_index);

		if (!vAssocStnList->at(to_station_index))
			vAssocStnList->at(to_station_index).reset(new CAStationList);
		vAssocStnList->at(to_station_index).get()->IncrementMsrCount();

		// Dual station measurements...
		switch (_it_msr->get()->GetTypeC())
		{
		case 'B':	// Geodetic azimuth
			g_map_tally.B ++;
			(*lMapCount)++;
			continue;
		case 'D':	// Direction set
			// now map the directions in the set
			vdirns = _it_msr->get()->GetDirections_ptr();
			_it_msr->get()->SetStn3Index(static_cast<UINT32>(vdirns->size()));
			MapMeasurementStationsDir(vdirns, vAssocStnList, lMapCount);
			continue;
		case 'K':	// Astronomic azimuth
			g_map_tally.K ++;
			(*lMapCount)++;
			continue;
		case 'C':	// Chord dist
			g_map_tally.C ++;
			(*lMapCount)++;
			continue;
		case 'E':	// Ellipsoid arc
			g_map_tally.E ++;
			(*lMapCount)++;
			continue;
		case 'M':	// MSL arc
			g_map_tally.M ++;
			(*lMapCount)++;
			continue;
		case 'S':	// Slope distance
			g_map_tally.S ++;
			(*lMapCount)++;
			continue;
		case 'L':	// Level difference
			g_map_tally.L ++;
			(*lMapCount)++;
			continue;
		case 'V':	// Zenith distance
			g_map_tally.V ++;
			(*lMapCount)++;
			continue;
		case 'Z':	// Vertical angle
			g_map_tally.Z ++;
			(*lMapCount)++;
			continue;
		}

		if (_it_msr->get()->GetTarget2().empty())
			throw XMLInteropException("Empty \"Third\" station name", 0);

		// 3. All measurements which have <Target2> or <Third>
		to2_station_name = _it_msr->get()->GetTarget2();	// GetTarget2() invokes specialised class (i.e. CDnaAngle, etc) member function
		
		it_stnmap_range = equal_range(vStnsMap_sortName_.begin(), vStnsMap_sortName_.end(), to2_station_name, StationNameIDCompareName());

		if (it_stnmap_range.first == it_stnmap_range.second)
		{
			ss.str("");
			ss << to2_station_name << " was not loaded from the stations map.  Please ensure that " << to2_station_name << " is included in the list of network stations.";
			throw XMLInteropException(ss.str(), 0);
		}

		// Increment ASL associated measurement count (i.e. CAStationList.assocMsrs_)
		to2_station_index = it_stnmap_range.first->second;
		
		_it_msr->get()->SetStn3Index(to2_station_index);

		if (!vAssocStnList->at(to2_station_index))
			vAssocStnList->at(to2_station_index).reset(new CAStationList);
		vAssocStnList->at(to2_station_index).get()->IncrementMsrCount();

		// at this stage, only 'A' type measurements have a Third
		switch (_it_msr->get()->GetTypeC())
		{
		case 'A':	// Horizontal angles
			g_map_tally.A ++;
			(*lMapCount)++;
			continue;
		}
	}

	vASLPtr::iterator _it_asl(vAssocStnList->begin()), _it_asl_begin(vAssocStnList->begin());
	UINT32 d;

	// find unused stations and store them in a vector, then reset 
	for (; _it_asl!=vAssocStnList->end(); _it_asl++)
	{
		if (!_it_asl->get())	// unused station
		{
			d = static_cast<UINT32>(std::distance(_it_asl_begin, _it_asl));
			vUnusedStns->push_back(vStnsMap_sortName_.at(d).first);
			_it_asl->reset(new CAStationList(false));
		}			
	}
}
	

UINT32 dna_import::ComputeMeasurementCount(vdnaMsrPtr* vMeasurements, const vUINT32& vIgnoredMsrs)
{
	MsrTally msr_tally;

	it_vUINT32_const it_ignmsr;
	_it_vdnamsrptr _it_msr;

	for (it_ignmsr=vIgnoredMsrs.begin(); it_ignmsr!=vIgnoredMsrs.end(); ++it_ignmsr)
	{
		_it_msr = vMeasurements->begin() + *it_ignmsr;

		// Increment single station measurement counters...
		switch (_it_msr->get()->GetTypeC())
		{
		case 'A': // Horizontal angle
			msr_tally.A++;
			break;
		case 'B': // Geodetic azimuth
			msr_tally.B++;
			break;
		case 'C': // Chord dist
			msr_tally.C++;
			break;
		case 'D': // Direction set
			msr_tally.D += static_cast<UINT32>(_it_msr->get()->GetDirections_ptr()->size());
			break;
		case 'E': // Ellipsoid arc
			msr_tally.E++;
			break;
		case 'G': // GPS Baseline (treat as single-baseline cluster)
			msr_tally.G += static_cast<UINT32>(_it_msr->get()->GetBaselines_ptr()->size() * 3);
			break;
		case 'X': // GPS Baseline cluster
			msr_tally.X += static_cast<UINT32>(_it_msr->get()->GetBaselines_ptr()->size() * 3);
			break;
		case 'H': // Orthometric height
			msr_tally.H++;
			break;
		case 'I': // Astronomic latitude
			msr_tally.I++;
			break;
		case 'J': // Astronomic longitude
			msr_tally.J++;
			break;
		case 'K': // Astronomic azimuth
			msr_tally.K++;
			break;
		case 'L': // Level difference
			msr_tally.L++;
			break;
		case 'M': // MSL arc
			msr_tally.M++;
			break;
		case 'P': // Geodetic latitude
			msr_tally.P++;
			break;
		case 'Q': // Geodetic longitude
			msr_tally.Q++;
			break;
		case 'R': // Ellipsoidal height
			msr_tally.R++;
			break;
		case 'S': // Slope distance
			msr_tally.S++;
			break;
		case 'V': // Zenith distance
			msr_tally.V++;
			break;
		case 'Y': // GPS point cluster
			msr_tally.Y += static_cast<UINT32>(_it_msr->get()->GetPoints_ptr()->size() * 3);
			break;
		case 'Z': // Vertical angle
			msr_tally.Z++;
			break;
		default:
			std::stringstream ss;
			ss << "ComputeMeasurementCount(): Unknown measurement type:  " << _it_msr->get()->GetTypeC() << std::endl;
			throw XMLInteropException(ss.str(), m_lineNo); 
		}
	}

	return msr_tally.TotalCount();
}

void dna_import::MapMeasurementStationsBsl(std::vector<CDnaGpsBaseline>* vGpsBaselines, pvASLPtr vAssocStnList, PUINT32 lMapCount)
{
	const size_t mapsize = vStnsMap_sortName_.size();
	if (mapsize < 1)
		throw XMLInteropException("A station map has not been created. Run \"SortStations\" to create a station map.", 0);

	std::ostringstream ss;
	//_it_string_uint32_pair _it_stnmap(vStnsMap_sortName_.end());
	std::vector< CDnaGpsBaseline >::iterator _it_msr(vGpsBaselines->begin());
	
	std::string station_name;
	UINT32 station_index; //, msrs_per_cluster_row;

#ifdef _MSDEBUG
	CAStationList* stnList;
#endif

	// Unique list of stations involved in this cluster
	vUINT32 msrStations;
	//msrStations.reserve(vGpsBaselines->size());

	for (_it_msr = vGpsBaselines->begin(); 
		_it_msr != vGpsBaselines->end(); 
		_it_msr++)
	{
		//msrs_per_cluster_row = 3 + static_cast<UINT32>(_it_msr->GetCovariances_ptr()->size() * 3);
		
		// <First> station
		station_name = _it_msr->GetFirst();
		
		it_stnmap_range = equal_range(vStnsMap_sortName_.begin(), vStnsMap_sortName_.end(), station_name, StationNameIDCompareName());

		if (it_stnmap_range.first == it_stnmap_range.second)
		{
			ss.str("");
			ss << station_name << " is referenced in measurement '" << _it_msr->GetType() << "' " <<
				_it_msr->GetFirst() << " - " << _it_msr->GetTarget() << std::endl;
			ss << "  but is not in the stations map.  " << std::endl;
			ss << "  Please ensure that " << station_name << " is included in the list of stations.";
			throw XMLInteropException(ss.str(), 0);
		}

		// Increment ASL associated measurement count (i.e. CAStationList.assocMsrs_)
		station_index = it_stnmap_range.first->second;
		
		_it_msr->SetStn1Index(station_index);
		if (!vAssocStnList->at(station_index))
			vAssocStnList->at(station_index).reset(new CAStationList);

		msrStations.push_back(station_index);

#ifdef _MSDEBUG
		stnList = vAssocStnList->at(station_index).get();
#endif

		// <Second> station
		station_name = _it_msr->GetTarget();

		it_stnmap_range = equal_range(vStnsMap_sortName_.begin(), vStnsMap_sortName_.end(), station_name, StationNameIDCompareName());

		if (it_stnmap_range.first == it_stnmap_range.second)
		{
			ss.str("");
			ss << station_name << " is referenced in measurement '" << _it_msr->GetType() << "' " <<
				_it_msr->GetFirst() << " - " << _it_msr->GetTarget() << std::endl;
			ss << "  but is not in the stations map.  " << std::endl;
			ss << "  Please ensure that " << station_name << " is included in the list of stations.";
			throw XMLInteropException(ss.str(), 0);
		}

		// Increment ASL associated measurement count (i.e. CAStationList.assocMsrs_)
		station_index = it_stnmap_range.first->second;
		
		_it_msr->SetStn2Index(station_index);
		if (!vAssocStnList->at(station_index))
			vAssocStnList->at(station_index).reset(new CAStationList);

		msrStations.push_back(station_index);

#ifdef _MSDEBUG
		stnList = vAssocStnList->at(station_index).get();
#endif

		(*lMapCount) += 3;		// one per vector component
		switch (_it_msr->GetTypeC())
		{
		case 'G':
			g_map_tally.G += 3;
			break;
		case 'X':
			g_map_tally.X += 3;
			break;
		}		
	}

	// Strip duplicates from msrStations, then increment station count for each of the stations tied to this cluster
	strip_duplicates(msrStations);
	it_vUINT32 _it_stn;
	for (_it_stn=msrStations.begin(); _it_stn!=msrStations.end(); ++_it_stn)
		vAssocStnList->at(*_it_stn).get()->IncrementMsrCount();

	if (vGpsBaselines->begin()->GetTypeC() != 'X')
		return;

	// set ID for each covariance term.  This is for the convenience of assigning covariance terms
	// to the right index in dnaAdjust
	std::vector< CDnaCovariance >::iterator _it_cov;
	std::vector< CDnaGpsBaseline >::iterator _it_msr2(vGpsBaselines->begin());

	for (_it_msr = vGpsBaselines->begin(); _it_msr != vGpsBaselines->end(); _it_msr++)
	{
		_it_msr2 = _it_msr;
		for (_it_cov = _it_msr->GetCovariances_ptr()->begin(); _it_cov != _it_msr->GetCovariances_ptr()->end(); _it_cov++)
		{
			_it_cov->SetStn1Index(_it_msr2->GetStn1Index());
			_it_cov->SetStn2Index((++_it_msr2)->GetStn2Index());
		}
	}	
}

void dna_import::MapMeasurementStationsPnt(std::vector<CDnaGpsPoint>* vGpsPoints, pvASLPtr vAssocStnList, PUINT32 lMapCount)
{
	const size_t mapsize = vStnsMap_sortName_.size();
	if (mapsize < 1)
		throw XMLInteropException("A station map has not been created. Run \"SortStations\" to create a station map.", 0);

	std::ostringstream ss;
	//_it_string_uint32_pair _it_stnmap;
	std::vector< CDnaGpsPoint >::iterator _it_msr;

	std::string station_name;
	UINT32 station_index; //, msrs_per_cluster_row;

#ifdef _MSDEBUG
	CAStationList* stnList;
#endif

	for (_it_msr=vGpsPoints->begin(); _it_msr != vGpsPoints->end(); _it_msr++)
	{
		// <First> station
		station_name = _it_msr->GetFirst();
		
		it_stnmap_range = equal_range(vStnsMap_sortName_.begin(), vStnsMap_sortName_.end(), station_name, StationNameIDCompareName());

		if (it_stnmap_range.first == it_stnmap_range.second)
		{
			ss.str("");
			ss << station_name << " is not in the stations map.  Please ensure that " << station_name << " is included in the list of stations.";
			throw XMLInteropException(ss.str(), 0);
		}

		// Increment ASL associated measurement count (i.e. CAStationList.assocMsrs_)
		station_index = it_stnmap_range.first->second;
		
		_it_msr->SetStn1Index(station_index);
		if (!vAssocStnList->at(station_index))
			vAssocStnList->at(station_index).reset(new CAStationList);
		vAssocStnList->at(station_index).get()->IncrementMsrCount();

#ifdef _MSDEBUG
		stnList = vAssocStnList->at(station_index).get();
#endif

		(*lMapCount) += 3;		// one per vector component
		g_map_tally.Y += 3;
	}

	// set ID for each covariance term.  This is for the convenience of assigning covariance terms
	// to the right index in dnaAdjust
	std::vector< CDnaCovariance >::iterator _it_cov;
	std::vector< CDnaGpsPoint >::iterator _it_msr2(vGpsPoints->begin());

	for (_it_msr = vGpsPoints->begin(); _it_msr != vGpsPoints->end(); _it_msr++)
	{
		_it_msr2 = _it_msr;
		for (_it_cov = _it_msr->GetCovariances_ptr()->begin(); _it_cov != _it_msr->GetCovariances_ptr()->end(); _it_cov++)
		{
			++_it_msr2;
			_it_cov->SetStn1Index(_it_msr2->GetStn1Index());
		}
	}
}


void dna_import::MapMeasurementStationsDir(std::vector<CDnaDirection>* vDirections, pvASLPtr vAssocStnList, PUINT32 lMapCount)
{
	const size_t mapsize = vStnsMap_sortName_.size();
	if (mapsize < 1)
		throw XMLInteropException("A station map has not been created. Run \"SortStations\" to create a station map.", 0);

	std::ostringstream ss;
	//_it_string_uint32_pair _it_stnmap;
	std::vector< CDnaDirection >::iterator _it_msr(vDirections->begin());

	std::string at_station_name, to_station_name;
	UINT32 at_station_index, to_station_index;

	// Unique list of stations involved in this cluster
	vUINT32 msrStations;
	msrStations.reserve(vDirections->size());

	for (_it_msr=vDirections->begin(); _it_msr!=vDirections->end(); _it_msr++)
	{
		// <First> station.  Strictly speaking, this value should not change from 
		// one direction to the next in the set
		at_station_name = _it_msr->GetFirst();
		
		it_stnmap_range = equal_range(vStnsMap_sortName_.begin(), vStnsMap_sortName_.end(), at_station_name, StationNameIDCompareName());

		if (it_stnmap_range.first == it_stnmap_range.second)
		{
			ss.str("");
			ss << at_station_name << " is not in the stations map.  Please ensure that " << at_station_name << " is included in the list of stations.";
			throw XMLInteropException(ss.str(), 0);
		}

		at_station_index = it_stnmap_range.first->second;
		
		// This is needed to initialise the at station index for every
		// target direction
		_it_msr->SetStn1Index(at_station_index);

		msrStations.push_back(at_station_index);

		// <Second> station
		to_station_name = _it_msr->GetTarget();
		
		it_stnmap_range = equal_range(vStnsMap_sortName_.begin(), vStnsMap_sortName_.end(), to_station_name, StationNameIDCompareName());

		if (it_stnmap_range.first == it_stnmap_range.second)
		{
			ss.str("");
			ss << to_station_name << " was not loaded from the stations map.  Please ensure that " << to_station_name << " is included in the list of network stations.";
			throw XMLInteropException(ss.str(), 0);
		}

		to_station_index = it_stnmap_range.first->second;
		
		_it_msr->SetStn2Index(to_station_index);
		if (!vAssocStnList->at(to_station_index))
			vAssocStnList->at(to_station_index).reset(new CAStationList);
		
		msrStations.push_back(to_station_index);
		
		if (_it_msr->NotIgnored())
		{
			(*lMapCount)++;
			g_map_tally.D++;
		}
	}

	// Strip duplicates from msrStations, then increment station count for each of the stations tied to this cluster
	strip_duplicates(msrStations);
	for (it_vUINT32 _it_stn=msrStations.begin(); _it_stn!=msrStations.end(); ++_it_stn)
		vAssocStnList->at(*_it_stn).get()->IncrementMsrCount();
}

void dna_import::CompleteAssociationLists(vdnaMsrPtr* vMeasurements, pvASLPtr vAssocStnList, pvUINT32 vAssocMsrList, const _AML_TYPE_ aml_type)
{
	CAStationList* currentASL;
	
	UINT32 msrstoAllStations = 0;
	UINT32 currentBmsFileIndex;
	UINT32 stnIndex, stnmsr_indexAML, msrCount;

	vASLPtr::iterator _it_asl(vAssocStnList->begin());

#ifdef _MSDEBUG
	CAStationList* stnList;
#endif

	// Complete ASL (first stage, i.e. IncrementMsrCount, was done in MapMeasurementStations)
	for (_it_asl=vAssocStnList->begin(); _it_asl!=vAssocStnList->end(); _it_asl++)
	{
		if (!_it_asl->get())		// unused station
			continue;
		if (_it_asl->get()->IsInvalid())		// unused station
			continue;
		
#ifdef _MSDEBUG
		stnList = _it_asl->get();
#endif
		_it_asl->get()->SetAMLStnIndex(msrstoAllStations);		// On first run, msrstoAllStations = 0!
		// increment total count (for dimension of AML)
		msrCount = _it_asl->get()->GetAssocMsrCount();
		msrstoAllStations += msrCount;

	}

	// Now create AML
	vAssocMsrList->clear();
	vAssocMsrList->resize(msrstoAllStations);

	std::vector<CDnaDirection>* vDirns;
	std::vector<CDnaGpsBaseline>* vGpsBsls;
	std::vector<CDnaGpsPoint>* vGpsPnts;

	// reset ASL #measurements to zero
	for (_it_asl=vAssocStnList->begin(); _it_asl!=vAssocStnList->end(); _it_asl++)
	{
		if (!_it_asl->get())					// unused station
			continue;
		if (_it_asl->get()->IsInvalid())		// unused station
			continue;
		_it_asl->get()->InitAssocMsrCount();

		// initialise all stations to 'invalid', then 
		// set 'valid' upon reading a non-ignored measurement
		_it_asl->get()->SetInvalid();
	}

	currentBmsFileIndex = 0;
	_it_vdnamsrptr _it_msr(vMeasurements->begin());
	for (_it_msr=vMeasurements->begin(); _it_msr != vMeasurements->end(); _it_msr++)
	{
		// 1. Handle nested type measurements (G, X, Y) separately
		switch (_it_msr->get()->GetTypeC())
		{
		case 'D':	// Direction set
			vDirns = _it_msr->get()->GetDirections_ptr();
			CompleteASLDirections(_it_msr, vDirns, vAssocStnList, vAssocMsrList, &currentBmsFileIndex, aml_type);
			continue;
		case 'G':	// GPS Baseline (treat as single-baseline cluster)
		case 'X':	// GPS Baseline cluster
			vGpsBsls = _it_msr->get()->GetBaselines_ptr();
			CompleteASLGpsBaselines(vGpsBsls, vAssocStnList, vAssocMsrList, &currentBmsFileIndex, aml_type);
			continue;
		case 'Y':	// GPS point cluster
			vGpsPnts = _it_msr->get()->GetPoints_ptr();
			CompleteASLGpsPoints(vGpsPnts, vAssocStnList, vAssocMsrList, &currentBmsFileIndex, aml_type);
			continue;
		}

		// calc AML index for <First>
		stnIndex = ((CDnaMeasurement*)_it_msr->get())->GetStn1Index();
		if (!vAssocStnList->at(stnIndex))
			throw XMLInteropException("CompleteAssociationLists(): An invalid index was found in the station map \
									while \ntrying to determine ASL index for the First station.", 0);

		currentASL = vAssocStnList->at(stnIndex).get();
		stnmsr_indexAML = currentASL->GetAMLStnIndex() + currentASL->GetAssocMsrCount();
		if (stnmsr_indexAML > msrstoAllStations)
			throw XMLInteropException("CompleteAssociationLists(): An error occurred while trying to determine ASL index.", 0);
		vAssocMsrList->at(stnmsr_indexAML) = currentBmsFileIndex;
		currentASL->IncrementMsrCount();

		if (_it_msr->get()->NotIgnored())
			currentASL->SetValid();

		// single station measurements
		switch (_it_msr->get()->GetTypeC())
		{
		case 'H':	// Orthometric height
		case 'R':	// Ellipsoidal height
		case 'I':	// Astronomic latitude
		case 'J':	// Astronomic longitude
		case 'P':	// Geodetic latitude
		case 'Q':	// Geodetic longitude
			currentBmsFileIndex ++;
			continue;
		}

		// calc AML index for <Second>
		stnIndex = ((CDnaMeasurement*)_it_msr->get())->GetStn2Index();
		if (!vAssocStnList->at(stnIndex))
			throw XMLInteropException("CompleteAssociationLists(): An invalid index was found in the station map \
									while \ntrying to determine ASL index for the Second station.", 0);
		
		currentASL = vAssocStnList->at(stnIndex).get();
		stnmsr_indexAML = currentASL->GetAMLStnIndex() + currentASL->GetAssocMsrCount();
		if (stnmsr_indexAML > msrstoAllStations)
			throw XMLInteropException("CompleteAssociationLists(): An error occurred while trying to determine ASL index.", 0);
		vAssocMsrList->at(stnmsr_indexAML) = currentBmsFileIndex;
		currentASL->IncrementMsrCount();

		if (_it_msr->get()->NotIgnored())
			currentASL->SetValid();

		// dual station measurements
		switch (_it_msr->get()->GetTypeC())
		{
		case 'B':	// Geodetic azimuth
		case 'K':	// Astronomic azimuth
		case 'C':	// Chord dist
		case 'E':	// Ellipsoid arc
		case 'M':	// MSL arc
		case 'S':	// Slope distance
		case 'L':	// Level difference
		case 'V':	// Zenith distance
		case 'Z':	// Vertical angle
			currentBmsFileIndex ++;
			continue;
		}

		// calc AML index for <Third>
		stnIndex = ((CDnaMeasurement*)_it_msr->get())->GetStn3Index();
		if (!vAssocStnList->at(stnIndex))
			throw XMLInteropException("CompleteAssociationLists(): An invalid index was found in the station map \
									while \ntrying to determine ASL index for the Third station.", 0);
		
		currentASL = vAssocStnList->at(stnIndex).get();
		stnmsr_indexAML = currentASL->GetAMLStnIndex() + currentASL->GetAssocMsrCount();
		if (stnmsr_indexAML > msrstoAllStations)
			throw XMLInteropException("CompleteAssociationLists(): An error occurred while trying to determine ASL index.", 0);
		vAssocMsrList->at(stnmsr_indexAML) = currentBmsFileIndex;
		currentASL->IncrementMsrCount();

		if (_it_msr->get()->NotIgnored())
			currentASL->SetValid();

		// triple station measurements
		switch (_it_msr->get()->GetTypeC())
		{
		case 'A':	// Horizontal angles
			currentBmsFileIndex ++;
			continue;
		}
	}
	return;
}
	

//// Find unused stations resulting from ignored measurements and add to the vector
//// NOTE - this function is not used
//void dna_import::FindUnusedStationsInIgnoredMeasurements(vdnaMsrPtr* vMeasurements, pvASLPtr vAssocStnList, pvUINT32 vAssocMsrList, pvstring vUnusedStns, pvUINT32 vIgnoredMsrs)
//{
//	vASLPtr::iterator _it_asl(vAssocStnList->begin()), _it_asl_begin(vAssocStnList->begin());
//
//	it_vUINT32 it_ignmsr, _it_stn_newend;
//
//	std::vector<CDnaDirection>* vDirns;
//	std::vector<CDnaDirection>::iterator it_Dirns;
//	std::vector<CDnaGpsBaseline>* vGpsBsls;
//	std::vector<CDnaGpsBaseline>::iterator it_GpsBsls;
//	std::vector<CDnaGpsPoint>* vGpsPnts;
//	std::vector<CDnaGpsPoint>::iterator it_GpsPnts;
//
//	vUINT32 vIgnoredMsrStations;
//	// Iterate through the ignored measurement list and build unique list of stations in ignored measurements
//	for (it_ignmsr=vIgnoredMsrs->begin(); it_ignmsr!=vIgnoredMsrs->end(); ++it_ignmsr)
//	{
//		switch (vMeasurements->at(*it_ignmsr)->GetTypeC())
//		{
//		case 'A':	// Horizontal angles
//			vIgnoredMsrStations.push_back(it_Dirns->GetStn3Index());			
//		case 'B':	// Geodetic azimuth
//		case 'K':	// Astronomic azimuth
//		case 'C':	// Chord dist
//		case 'E':	// Ellipsoid arc
//		case 'M':	// MSL arc
//		case 'S':	// Slope distance
//		case 'L':	// Level difference
//		case 'V':	// Zenith distance
//		case 'Z':	// Vertical angle
//			vIgnoredMsrStations.push_back(it_Dirns->GetStn2Index());
//		case 'H':	// Orthometric height
//		case 'R':	// Ellipsoidal height
//		case 'I':	// Astronomic latitude
//		case 'J':	// Astronomic longitude
//		case 'P':	// Geodetic latitude
//		case 'Q':	// Geodetic longitude
//			vIgnoredMsrStations.push_back(it_Dirns->GetStn1Index());
//			break;
//		case 'D':	// Direction set
//			vDirns = vMeasurements->at(*it_ignmsr)->GetDirections_ptr();
//			vIgnoredMsrStations.push_back(it_Dirns->GetStn1Index());
//			vIgnoredMsrStations.push_back(it_Dirns->GetStn2Index());
//			for (it_Dirns=vDirns->begin(); it_Dirns!=vDirns->end(); ++it_Dirns)
//				vIgnoredMsrStations.push_back(it_Dirns->GetStn2Index());
//			break;
//		case 'G':	// GPS Baseline (treat as single-baseline cluster)
//		case 'X':	// GPS Baseline cluster
//			vGpsBsls = vMeasurements->at(*it_ignmsr)->GetBaselines_ptr();
//			for (it_GpsBsls=vGpsBsls->begin(); it_GpsBsls!=vGpsBsls->end(); ++it_GpsBsls)
//			{
//				vIgnoredMsrStations.push_back(it_GpsBsls->GetStn1Index());
//				vIgnoredMsrStations.push_back(it_GpsBsls->GetStn2Index());
//			}
//			break;
//		case 'Y':	// GPS point cluster
//			vGpsPnts = vMeasurements->at(*it_ignmsr)->GetPoints_ptr();
//			for (it_GpsPnts=vGpsPnts->begin(); it_GpsPnts!=vGpsPnts->end(); ++it_GpsPnts)
//				vIgnoredMsrStations.push_back(it_GpsPnts->GetStn1Index());
//			break;
//		}
//	}
//
//	//sort(vIgnoredMsrStations.begin(), vIgnoredMsrStations.end());
//	//_it_stn_newend = unique(vIgnoredMsrStations.begin(), vIgnoredMsrStations.end());
//	//if (_it_stn_newend != vIgnoredMsrStations.end())
//	//	vIgnoredMsrStations.resize(_it_stn_newend - vIgnoredMsrStations.begin());
//	strip_duplicates(vIgnoredMsrStations);
//}
	

void dna_import::CompleteASLDirections(_it_vdnamsrptr _it_msr, std::vector<CDnaDirection>* vDirections, pvASLPtr vAssocStnList, 
										 pvUINT32 vAssocMsrList, PUINT32 currentBmsFileIndex, const _AML_TYPE_ aml_type)
{
	std::ostringstream ss;
	
	CAStationList* currentASL;

	std::string station_name;
	UINT32 stn_indexAML, stnmsr_indexAML;
	const UINT32 msrstoAllStations = static_cast<UINT32>(vAssocMsrList->size());
	const UINT32 bmsrIndex(*currentBmsFileIndex);

	std::vector< CDnaDirection >::iterator _it_dir(vDirections->begin());

	// Unique list of stations involved in this cluster
	vUINT32 msrStations;
	msrStations.reserve(vDirections->size() + 1);

	// calc AML index for <First>
	stn_indexAML = _it_msr->get()->GetStn1Index();
	if (!vAssocStnList->at(stn_indexAML))
		throw XMLInteropException("CompleteASLDirections(): An invalid index was found in the station map \
								while \ntrying to determine ASL index for the First station.", 0);

	currentASL = vAssocStnList->at(stn_indexAML).get();
	stnmsr_indexAML = currentASL->GetAMLStnIndex() + currentASL->GetAssocMsrCount();
	if (stnmsr_indexAML > msrstoAllStations)
		throw XMLInteropException("CompleteASLDirections(): An error occurred while trying to determine ASL index.", 0);
	vAssocMsrList->at(stnmsr_indexAML) = bmsrIndex;

	if (_it_msr->get()->NotIgnored())
		currentASL->SetValid();

	// Add to unique list of stations
	msrStations.push_back(stn_indexAML);

	// calc AML index for <Second>
	stn_indexAML = _it_msr->get()->GetStn2Index();
	if (!vAssocStnList->at(stn_indexAML))
		throw XMLInteropException("CompleteASLDirections(): An invalid index was found in the station map \
								while \ntrying to determine ASL index for the Second station.", 0);
		
	currentASL = vAssocStnList->at(stn_indexAML).get();
	stnmsr_indexAML = currentASL->GetAMLStnIndex() + currentASL->GetAssocMsrCount();
	if (stnmsr_indexAML > msrstoAllStations)
		throw XMLInteropException("CompleteASLDirections(): An error occurred while trying to determine ASL index.", 0);
	vAssocMsrList->at(stnmsr_indexAML) = bmsrIndex;

	if (_it_msr->get()->NotIgnored())
		currentASL->SetValid();

	// Add to unique list of stations
	msrStations.push_back(stn_indexAML);
	
	// Increment measurement count, same for aml_type=str_msr and aml_type=cls_msr:
	// - Increment binary measurement count for Inst->RO, then once per direction
	// - Increment by 1 since a cluster of directions is handled by a single measurement class instance
	(*currentBmsFileIndex)++;

	// Calc AML index for Station - Target directions
	for (_it_dir=vDirections->begin(); _it_dir != vDirections->end(); ++_it_dir)
	{
		// Calc AML index for <Second>
		stn_indexAML = _it_dir->GetStn2Index();
		if (!vAssocStnList->at(stn_indexAML))
			throw XMLInteropException("CompleteASLDirections(): An invalid index was found in the station map \
									while \ntrying to determine ASL index for the Second station.", 0);
		currentASL = vAssocStnList->at(stn_indexAML).get();
		stnmsr_indexAML = currentASL->GetAMLStnIndex() + currentASL->GetAssocMsrCount();
		if (stnmsr_indexAML > msrstoAllStations)
			throw XMLInteropException("CompleteASLDirections(): An error occurred while trying to determine \
									AML index for the Second station.", 0);
		
		// All ASL entries are set to invalid on initial run of 
		// CompleteAssociationLists, so this is vital so as to ensure a station 
		// is to be included.  But stations should only be included if the 
		// parent direction set measurement is not ignored
		if (_it_dir->NotIgnored() && _it_msr->get()->NotIgnored())
			currentASL->SetValid();

		// Set binary msr index for all stations in this cluster to the first element of the 
		// cluster in the binary file
		vAssocMsrList->at(stnmsr_indexAML) = bmsrIndex;
		
		// To station
		msrStations.push_back(stn_indexAML);

		// increment file index for next measurement
		switch (aml_type)
		{
		case cls_msr:		// CDnaMeasurement class
			// do nothing, handled above.
			break;
		case str_msr:		// msr_t struct
		default:
			(*currentBmsFileIndex)++;
			break;
		}
	}

	// Strip duplicates from msrStations, then increment station count once
	// for each station associated with this cluster
	strip_duplicates(msrStations);
	for (it_vUINT32 _it_stn=msrStations.begin(); _it_stn!=msrStations.end(); ++_it_stn)
		vAssocStnList->at(*_it_stn).get()->IncrementMsrCount();
}

void dna_import::CompleteASLGpsBaselines(std::vector<CDnaGpsBaseline>* vGpsBaselines, pvASLPtr vAssocStnList,
										   pvUINT32 vAssocMsrList, PUINT32 currentBmsFileIndex, const _AML_TYPE_ aml_type)
{
	CAStationList* currentASL;

	UINT32 stn_indexAML, stnmsr_indexAML;
	const UINT32 msrstoAllStations = static_cast<UINT32>(vAssocMsrList->size());
	const UINT32 bmsrIndex(*currentBmsFileIndex);
	
	std::vector< CDnaGpsBaseline >::iterator _it_msr(vGpsBaselines->begin());

	// Unique list of stations involved in this cluster
	vUINT32 msrStations;
	msrStations.reserve(vGpsBaselines->size());

	for (_it_msr=vGpsBaselines->begin(); _it_msr != vGpsBaselines->end(); _it_msr++)
	{
		// Calc AML index for <First>
		stn_indexAML = _it_msr->GetStn1Index();
		if (!vAssocStnList->at(stn_indexAML))
			throw XMLInteropException("CompleteASLGpsBaselines(): An invalid index was found in the station \
									map while \ntrying to determine ASL index for the First station.", 0);
		currentASL = vAssocStnList->at(stn_indexAML).get();
		stnmsr_indexAML = currentASL->GetAMLStnIndex() + currentASL->GetAssocMsrCount();
		if (stnmsr_indexAML > msrstoAllStations)
			throw XMLInteropException("CompleteASLGpsBaselines(): An error occurred while trying to determine the \
									AML index for the First station.", 0);

		if (_it_msr->NotIgnored())
			currentASL->SetValid();

		// Set binary msr index for all stations in this cluster to the first element of the 
		// cluster in the binary file
		vAssocMsrList->at(stnmsr_indexAML) = bmsrIndex;

		msrStations.push_back(stn_indexAML);

		// Calc AML index for <Second>
		stn_indexAML = _it_msr->GetStn2Index();
		if (!vAssocStnList->at(stn_indexAML))
			throw XMLInteropException("CompleteASLGpsBaselines(): An invalid index was found in the station map \
									while \ntrying to determine ASL index for the Second station.", 0);
		currentASL = vAssocStnList->at(stn_indexAML).get();
		stnmsr_indexAML = currentASL->GetAMLStnIndex() + currentASL->GetAssocMsrCount();
		if (stnmsr_indexAML > msrstoAllStations)
			throw XMLInteropException("CompleteASLGpsBaselines(): An error occurred while trying to determine \
									AML index for the Second station.", 0);
		
		if (_it_msr->NotIgnored())
			currentASL->SetValid();

		// Set binary msr index for all stations in this cluster to the first element of the 
		// cluster in the binary file
		vAssocMsrList->at(stnmsr_indexAML) = bmsrIndex;

		msrStations.push_back(stn_indexAML);

		// increment file index for next measurement
		switch (aml_type)
		{
		case cls_msr:		// CDnaMeasurement class
			// do nothing, handled later.
			break;
		case str_msr:		// msr_t struct
		default:
			// add:
			// - 3 for (X, X var); (Y, Y var, XY covar); (Z, Zvar, XZ, YZ covar), and
			// - 3 for each covariance
			(*currentBmsFileIndex) = (*currentBmsFileIndex) + 3 + static_cast<UINT32>(_it_msr->GetCovariances_ptr()->size() * 3);
			break;
		}
	}

	switch (aml_type)
	{
	case cls_msr:
		// increment by 1 since a cluster of GNSS baselines is handled by a single measurement class instance
		(*currentBmsFileIndex)++;
		break;
	case str_msr:
		break;
	}

	// Strip duplicates from msrStations, then increment station count for each of the stations tied to this cluster
	strip_duplicates(msrStations);
	for (it_vUINT32 _it_stn=msrStations.begin(); _it_stn!=msrStations.end(); ++_it_stn)
		vAssocStnList->at(*_it_stn).get()->IncrementMsrCount();
}

void dna_import::CompleteASLGpsPoints(std::vector<CDnaGpsPoint>* vGpsPoints, pvASLPtr vAssocStnList, 
										pvUINT32 vAssocMsrList, PUINT32 currentBmsFileIndex, const _AML_TYPE_ aml_type)
{
	CAStationList* currentASL;

	UINT32 stn_indexAML, stnmsr_indexAML;
	const UINT32 msrstoAllStations = static_cast<UINT32>(vAssocMsrList->size());
	const UINT32 bmsrIndex(*currentBmsFileIndex);
	
	std::vector< CDnaGpsPoint >::iterator _it_msr(vGpsPoints->begin());

	for (_it_msr=vGpsPoints->begin(); _it_msr != vGpsPoints->end(); _it_msr++)
	{		
		// Calc AML index for <First>
		stn_indexAML = _it_msr->GetStn1Index();
		if (!vAssocStnList->at(stn_indexAML))
			throw XMLInteropException("CompleteASLGpsPoints(): An invalid index was found in the station \
									map while \ntrying to determine ASL index for the First station.", 0);
		currentASL = vAssocStnList->at(stn_indexAML).get();
		stnmsr_indexAML = currentASL->GetAMLStnIndex() + currentASL->GetAssocMsrCount();
		if (stnmsr_indexAML > msrstoAllStations)
			throw XMLInteropException("CompleteASLGpsPoints(): An error occurred while trying to determine the \
									AML index for the First station.", 0);

		if (_it_msr->NotIgnored())
			currentASL->SetValid();

		// Set binary msr index for all stations in this cluster to the first element of the 
		// cluster in the binary file
		vAssocMsrList->at(stnmsr_indexAML) = bmsrIndex;

		// Increment the measurement count
		currentASL->IncrementMsrCount();

		// increment file index for next measurement
		switch (aml_type)
		{
		case cls_msr:		// CDnaMeasurement class
			// do nothing, handled later.
			break;
		case str_msr:		// msr_t struct
		default:
			// add:
			// - 3 for (X, X var); (Y, Y var, XY covar); (Z, Zvar, XZ, YZ covar), and
			// - 3 for each covariance
			(*currentBmsFileIndex) += 3 + static_cast<UINT32>(_it_msr->GetCovariances_ptr()->size() * 3);
			break;
		}
	}

	switch (aml_type)
	{
	case cls_msr:
		// increment by 1 since a cluster of GNSS points is handled by a single measurement class instance
		(*currentBmsFileIndex) ++;
		break;
	case str_msr:
		break;
	}
}

_PARSE_STATUS_ dna_import::LoadDNAGeoidFile(const std::string& fileName, vdnaStnPtr* vStations)
{
	std::ifstream geo_file;
	try {
		// Load geoid file.  Throws runtime_error on failure.
		file_opener(geo_file, fileName, 
			std::ios::in, ascii, true);
	}
	catch (const std::runtime_error& e) {
		SignalExceptionInterop(e.what(), 0, NULL);
	}

	std::string geoidRec, station, nValue, primeMeridian, primeVertical;
	double value;

	_it_vdnastnptr _it_stn;
	m_lineNo = 0;

	while (geo_file)
	{
		if (geo_file.eof())
			break;
		
		try {
			getline(geo_file, geoidRec);
		}
		catch (...) {
			if (geo_file.eof())
			{
				geo_file.close();
				return parseStatus_;
			}
			std::stringstream ss;
			ss << "LoadDNAGeoidFile(): Could not read from the geoid file." << std::endl;
			m_columnNo = 0;
			throw XMLInteropException(ss.str(), m_lineNo);
		}

		m_lineNo++;

		if (m_lineNo < 2)
			continue;

		// blank or whitespace?
		if (trimstr(geoidRec).empty())			
			continue;

		// Ignore lines with comments
		if (geoidRec.compare(0, 1, "*") == 0)
			continue;

		geoidRec = trimstr(geoidRec);

		// name
		try {
			station = trimstr(geoidRec.substr(0, 40));
		}
		catch (...) {
			std::stringstream ss;
			ss << "LoadDNAGeoidFile(): Could not extract station name from the record:  " << std::endl << "    " << geoidRec << std::endl;
			m_columnNo = 0;
			throw XMLInteropException(ss.str(), m_lineNo);
		}

		it_stnmap_range = equal_range(vStnsMap_sortName_.begin(), vStnsMap_sortName_.end(), 
			station, StationNameIDCompareName());

		if (it_stnmap_range.first == it_stnmap_range.second)
			continue;	// No station exists for this geoid information

		// N value
		try {
			nValue = trimstr(geoidRec.substr(40, 10));
		}
		catch (...) {
			std::stringstream ss;
			ss << "LoadDNAGeoidFile(): Could not extract the N-value from the record:  " << std::endl << "    " << geoidRec << std::endl;
			m_columnNo = 41;
			throw XMLInteropException(ss.str(), m_lineNo);
		}
		
		// Deflection in prime meridian
		try {
			primeMeridian = trimstr(geoidRec.substr(50, 19));
		}
		catch (...) {
			std::stringstream ss;
			ss << "LoadDNAGeoidFile(): Could not extract the deflection in prime meridian from the record:  " << std::endl << "    " << geoidRec << std::endl;
			m_columnNo = 51;
			throw XMLInteropException(ss.str(), m_lineNo);
		}

		// Deflection in prime vertical
		try {
			primeVertical = trimstr(geoidRec.substr(69));
		}
		catch (...) {
			std::stringstream ss;
			ss << "LoadDNAGeoidFile(): Could not extract the deflection in prime meridian from the record:  " << std::endl << "    " << geoidRec << std::endl;
			m_columnNo = 70;
			throw XMLInteropException(ss.str(), m_lineNo);
		}		

		try {
			_it_stn = vStations->begin() + it_stnmap_range.first->second;
			_it_stn->get()->SetgeoidSep(DoubleFromString<float>(nValue));
			
			// correct station height only if orthometric
			if (_it_stn->get()->GetMyHeightSystem() == ORTHOMETRIC_type_i)
			{
				value = DoubleFromString<double>(nValue);
				_it_stn->get()->SetcurrentHeight_d(_it_stn->get()->GetHeight() + value);
				//_it_stn->get()->SetHeightSystem(ELLIPSOIDAL_type_i);
			}

			// deflection values in the grid file are in seconds, so convert to radians
			value = DoubleFromString<double>(primeMeridian);
			_it_stn->get()->SetmeridianDef(SecondstoRadians(value));

			value = DoubleFromString<double>(primeVertical);
			_it_stn->get()->SetverticalDef(SecondstoRadians(value));
		}
		catch (...)
		{
			SignalExceptionInterop("Could not extract data for " + station + " from " + fileName + ", ", 0, "i", &geo_file);	
		}		
	}
	
	geo_file.close();

	return parseStatus_;
}

}	// namespace dynamlinterop
}	// namespace dynadjust
	

