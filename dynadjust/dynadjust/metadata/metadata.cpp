// metadata.cpp : Defines the entry point for the console application.
//
/// \cond
#include <iostream>
#include <fstream>
#include <algorithm>
#include <iomanip>
#include <sstream>
#include <string>
#include <time.h>
#include <mutex>

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

#include <filesystem>
/// \endcond

std::mutex cout_mutex;
std::mutex import_file_mutex;

#include <include/exception/dnaexception.hpp>

#include <include/config/dnaversion.hpp>
#include <include/config/dnaconsts.hpp>
#include <include/config/dnaconsts-interface.hpp>
#include <include/config/dnaoptions.hpp>
#include <include/config/dnaoptions-interface.hpp>

#include <include/functions/dnastringfuncs.hpp>
#include <include/functions/dnaiostreamfuncs.hpp>
#include <include/functions/dnafilepathfuncs.hpp>
#include <include/functions/dnastrmanipfuncs.hpp>
#include <include/functions/dnastrutils.hpp>
#include <include/functions/dnatimer.hpp>

#include <include/parameters/dnadatum.hpp>

using namespace dynadjust;
using namespace dynadjust::exception;
using namespace dynadjust::datum_parameters;

void processGNSSMeasurement(std::string& dnaMeasurementGNSS, std::string& referenceFrame, std::string& epoch, size_t& lineNo, size_t& measurementUpdateCount)
{
	std::string tmp(dnaMeasurementGNSS);

	// Validate reference frame and epoch.  Throws on error
	CDnaDatum datum;
	datum.SetDatumFromName(referenceFrame, epoch);

	std::string FirstTag = "<First>";
	std::string referenceframeTagOpen = "<ReferenceFrame>";
	std::string referenceframeTagClose = "</ReferenceFrame>";
	std::string epochTagOpen = "<Epoch>";
	std::string epochTagClose = "</Epoch>";
	
	std::string referenceFrameElement, epochElement;
	referenceFrameElement = referenceframeTagOpen + referenceFrame + referenceframeTagClose + "\n        ";
	epochElement = epochTagOpen + epoch + epochTagClose + "\n        ";

	size_t pos;

	bool updated(false);

	///////////////////////////////////////////////////////////
	// Look for the existence of <ReferenceFrame> in data
	// If none, add new <ReferenceFrame>
	if ((pos = tmp.find(referenceframeTagOpen)) == std::string::npos)
	{
		// Okay, no reference frame tag.  Add it before <First>
		if ((pos = tmp.find(FirstTag)) == std::string::npos)
		{
			std::stringstream ss;
			ss << "processGNSSMeasurement(): Could not find <First> element in:" << std::endl << 
				dnaMeasurementGNSS << std::endl;
			throw XMLInteropException(ss.str(), (int)lineNo);
		}

		tmp.insert(pos, referenceFrameElement);
		updated = true;
	}

	///////////////////////////////////////////////////////////
	// Look for the existence of <Epoch> in data
	// If none, add new <ReferenceFrame>
	if ((pos = tmp.find(epochTagOpen)) == std::string::npos)
	{
		// Okay, no reference frame tag.  Add it before <First>
		if ((pos = tmp.find(FirstTag)) == std::string::npos)
		{
			std::stringstream ss;
			ss << "processGNSSMeasurement(): Could not find <First> element in:" << std::endl << 
				dnaMeasurementGNSS << std::endl;
			throw XMLInteropException(ss.str(), (int)lineNo);
		}

		tmp.insert(pos, epochElement);
		updated = true;
	}

	dnaMeasurementGNSS = tmp;
	if (updated)
		measurementUpdateCount++;
}

void GetDefaultReferenceFrameEpoch(std::string& sBuf, std::string& referenceFrame, std::string& epoch)
{
	size_t pos;
	std::string referenceFrameAttribute = "referenceframe";
	std::string epochAttribute = "epoch";

	std::string tmp;

	if ((pos = sBuf.find(referenceFrameAttribute)) != std::string::npos)
	{
		// Okay, get the reference frame
		// Advance to =
		while (sBuf.at(pos) != '\"')
			pos++;
		pos++;

		tmp = sBuf.substr(pos);
		if ((pos = tmp.find("\"")) != std::string::npos)
			tmp = tmp.substr(0, pos);

		referenceFrame = tmp;
	}

	if ((pos = sBuf.find(epochAttribute)) != std::string::npos)
	{
		// Okay, get the epoch
		// Advance to =
		while (sBuf.at(pos) != '\"')
			pos++;
		pos++;

		tmp = sBuf.substr(pos);
		if ((pos = tmp.find("\"")) != std::string::npos)
			tmp = tmp.substr(0, pos);

		epoch = tmp;
	}
}
	

void processFile(std::ifstream* ifsFile, std::ofstream* ofsFile, std::string& defaultFrame, std::string& defaultEpoch, 
	size_t& measurementUpdateCount, size_t& unsureFrameCount, size_t& unsureEpochCount)
{
	std::string sBuf, tmp;

	std::string referenceFrame(GDA94_s);
	std::string epoch("01.01.1994");
	
	std::string openingDnaMeasurementTag = "<DnaMeasurement>";
	std::string closingDnaMeasurementTag = "</DnaMeasurement>";

	std::string dnaxmlFormatTag = "<DnaXmlFormat";

	std::string gType = "<Type>G</Type>";
	std::string xType = "<Type>X</Type>";
	std::string yType = "<Type>Y</Type>";

	const std::string xmlComment = "<!--";
	// <!--[BASELINE REFFRAME] ITRF2005-->
	const std::string metadataReferenceFrame = "[BASELINE REFFRAME]";
	// <!--[BASELINE EPOCH] 19.08.2010-->
	const std::string metadataEpoch = "[BASELINE EPOCH]";

	std::stringstream dnameasurementStream;
	bool newMeasurement(false);

	dnameasurementStream.str("");

	size_t lineNo = 0;
	size_t pos;
	
	CDnaDatum datum;

	while (ifsFile)
	{
		if (ifsFile->eof())
			break;

		lineNo++;

		try {
			getline((*ifsFile), sBuf);			
		}
		catch (...) {
			if (ifsFile->eof())
				return;
			
			std::stringstream ss;
			ss << "processFile(): Could not read from the station file." << std::endl;
			throw XMLInteropException(ss.str(), (int)lineNo);
		}

		// Reached a DnaMeasurement opening tag?  
		
		if (sBuf.find(dnaxmlFormatTag) != std::string::npos)
		{
			// Get the default reference frame and epoch
			GetDefaultReferenceFrameEpoch(sBuf, referenceFrame, epoch);

			// Validate reference frame and epoch.  Throws on error
			datum.SetDatumFromName(referenceFrame, epoch);

			// Set defaults
			defaultFrame = referenceFrame;
			defaultEpoch = epoch;
		}

		// Reached a DnaMeasurement opening tag?
		if (sBuf.find(openingDnaMeasurementTag) != std::string::npos)
		{	
			// yes, found an opening tag!
			newMeasurement = true;

			// Reset the stream
			dnameasurementStream.str("");

			// Add this string to the stringstream (for subsequent processing)
			dnameasurementStream << sBuf << std::endl;
		
		}
		// Reached a DnaMeasurement closing tag?
		else if (sBuf.find(closingDnaMeasurementTag) != std::string::npos)
		{	
			// yes, found a closing tag!
			dnameasurementStream << sBuf << std::endl;

			// capture the measurement
			std::string dnaMeasurementGNSS = dnameasurementStream.str();

			// Is this a GNSS measurement?  if not, dump and continue
			if (dnaMeasurementGNSS.find(gType) != std::string::npos)
			{
				// Process GNSS measurement
				processGNSSMeasurement(dnaMeasurementGNSS, referenceFrame, epoch, lineNo, measurementUpdateCount);
			}
			else if (dnaMeasurementGNSS.find(xType) != std::string::npos)
			{
				// Process GNSS measurement
				processGNSSMeasurement(dnaMeasurementGNSS, referenceFrame, epoch, lineNo, measurementUpdateCount);
			}
			else if (dnaMeasurementGNSS.find(yType) != std::string::npos)
			{
				// Process GNSS measurement
				processGNSSMeasurement(dnaMeasurementGNSS, referenceFrame, epoch, lineNo, measurementUpdateCount);
			}
			
			// Write dnaMeasurementGNSS to output file
			(*ofsFile) << dnaMeasurementGNSS;

			newMeasurement = false;
		}
		else
		{
			// At this point, the line is one of the following:
			// - A line before the first measurement (in which case, newMeasurement is false)
			// - An XML element within DnaMeasurement (in which case, newMeasurement is true)
			// - A comment within DnaMeasurement (in which case, newMeasurement is true)

			// Is this line within a DnaMeasurement
			if (newMeasurement)
			{
				// Yes, add this string to the stringstream (for subsequent processing)
				dnameasurementStream << sBuf << std::endl;
			}
			else
			{
				// A new measurement hasn't been reached yet, so this line is either:
				// - Before the first DnaMeasurement
				// - A comment preceding a DnaMeasurement

				// Dump the line to output file
				(*ofsFile) << sBuf << std::endl;
			}			

			// Now check for Reference frame and epoch comments
			tmp = trimstr(sBuf);

			// Is this line a comment?
			if (tmp.find(xmlComment) == std::string::npos)
			{
				// No, this is a measurement element.
				
				// Check for an existing ReferenceFrame tag
			}
			else
			{
				// At this point, tmp is a comment.

				// Is this line a reference frame metadata tag in the comment?
				if (tmp.find(metadataReferenceFrame) != std::string::npos)
				{
					if ((pos = tmp.find("]", 0)) != std::string::npos)
					{
						referenceFrame = tmp.substr(pos+1);
						if ((pos = referenceFrame.find("-->", 0)) != std::string::npos)
							referenceFrame = trimstr(referenceFrame.substr(0, pos));

						if (iequals(referenceFrame, "unsure"))
						{
							referenceFrame = defaultFrame;
							unsureFrameCount++;
						}

						if (iequals(referenceFrame, "WGS84"))
						{
							referenceFrame = "ITRF1989";
						}
					}
				}

				// Is this line a epoch metadata tag in the comment?
				else if (tmp.find(metadataEpoch) != std::string::npos)
				{
					if ((pos = tmp.find("]", 0)) != std::string::npos)
					{
						epoch = tmp.substr(pos+1);
						if ((pos = epoch.find("-->", 0)) != std::string::npos)
							epoch = trimstr(epoch.substr(0, pos));

						if (iequals(epoch, "unsure"))
						{
							epoch = "01.01.2005";
							unsureEpochCount++;
						}
					}
				}
			}
		}		
	}

}

int main(int argc, char* argv[])
{
	std::stringstream ss_err;

	std::string cmd_line_banner;
	fileproc_help_header(&cmd_line_banner);

	project_settings p;
	boost::program_options::variables_map vm;
	boost::program_options::positional_options_description positional_options;

	boost::program_options::options_description standard_options("+ " + std::string(ALL_MODULE_STDOPT), PROGRAM_OPTIONS_LINE_LENGTH);
	boost::program_options::options_description generic_options("+ " + std::string(ALL_MODULE_GENERIC), PROGRAM_OPTIONS_LINE_LENGTH);

	std::string cmd_line_usage("+ ");
	cmd_line_usage.append(__BINARY_NAME__).append(" usage:  ").append(__BINARY_NAME__).append(" ").append(" [options]");
	boost::program_options::options_description allowable_options(cmd_line_usage, PROGRAM_OPTIONS_LINE_LENGTH);

	try {
		standard_options.add_options()
			(IMPORT_FILE_F, boost::program_options::value< vstring >(&p.i.input_files), 
				"Station and measurement input file(s). Switch is not required.")
			;

		generic_options.add_options()
			(VERBOSE, boost::program_options::value<UINT16>(&p.g.verbose),
				std::string("Give detailed information about what ").append(__BINARY_NAME__).append(" is doing.\n  0: No information (default)\n  1: Helpful information\n  2: Extended information\n  3: Debug level information").c_str())
			(QUIET,
				std::string("Suppresses all explanation of what ").append(__BINARY_NAME__).append(" is doing unless an error occurs.").c_str())
			(VERSION_V, "Display the current program version.")
			(HELP_H, "Show this help message.")
			(HELP_MODULE_H, boost::program_options::value<std::string>(),
				"Provide help for a specific help category.")
			;

		allowable_options.add(standard_options).add(generic_options);

		// add "positional options" to handle command line tokens which have no option name
		positional_options.add(IMPORT_FILE, -1);

		boost::program_options::command_line_parser parser(argc, argv);
		store(parser.options(allowable_options).positional(positional_options).run(), vm);
		notify(vm);
	} 
	catch (const std::exception& e) {
		cout_mutex.lock();
		std::cout << "- Error: " << e.what() << std::endl;
		std::cout << cmd_line_banner << allowable_options << std::endl;
		cout_mutex.unlock();
		return EXIT_FAILURE;
	}

	if (argc < 2)
	{
		std::cout << std::endl << "- Nothing to do - no files or options provided. " << std::endl << std::endl;  
		std::cout << cmd_line_banner << allowable_options << std::endl;
		return EXIT_FAILURE;
	}

	if (vm.count(VERSION))
	{
		std::cout << cmd_line_banner << std::endl;
		return EXIT_SUCCESS;
	}

	if (vm.count(HELP))
	{
		std::cout << cmd_line_banner << allowable_options << std::endl;
		return EXIT_SUCCESS;
	}

	if (vm.count(HELP_MODULE))
	{
		std::cout << cmd_line_banner;
		std::string original_text = vm[HELP_MODULE].as<std::string>();
		std::string help_text = str_upper<std::string>(original_text);

		if (str_upper<std::string, char>(ALL_MODULE_STDOPT).find(help_text) != std::string::npos) {
			std::cout << standard_options << std::endl;
		}
		else if (str_upper<std::string, char>(ALL_MODULE_GENERIC).find(help_text) != std::string::npos) {
			std::cout << generic_options << std::endl;
		}
		else {
			std::cout << std::endl << "- Error: Help module '" <<
				original_text << "' is not in the list of options." << std::endl;
			return EXIT_FAILURE;
		}

		return EXIT_SUCCESS;
	}

	// Normalise files using input folder
	for_each(p.i.input_files.begin(), p.i.input_files.end(),
		[&p] (std::string& file) { 
			formPath<std::string>(p.g.input_folder, file);
	}
	);

	UINT32 stnCount(0), msrCount(0), clusterID(0);

	size_t i, nfiles(p.i.input_files.size());		// for each file...
	std::string input_file, output_file, ss;
	std::string defaultFrame, defaultEpoch;
	
	std::ostringstream ss_time, ss_msg;
	boost::posix_time::milliseconds elapsed_time(boost::posix_time::milliseconds(0));
	boost::posix_time::ptime pt;

	if (!p.g.quiet)
	{
		std::cout << std::endl << cmd_line_banner;
		std::cout << "+ Parsing: " << std::endl;
	}
	
	dynadjust::cpu_timer time;	// constructor of cpu_timer calls start()

	std::ifstream*	ifsDynaML_;
	size_t		sifsFileSize_, measurementUpdateCount(0), unsureFrameCount(0), unsureEpochCount(0);

	std::ofstream*	ofsDynaML_;

	ifsDynaML_ = 0;
	ofsDynaML_ = 0;

	size_t strlen_arg = 0;
	for_each(p.i.input_files.begin(), p.i.input_files.end(),
		[&strlen_arg](std::string& file) {
			if (leafStr<std::string>(file).length() > strlen_arg)
				strlen_arg = leafStr<std::string>(file).length();
	});
	strlen_arg += (6 + PROGRESS_PERCENT_04);

	for (i=0; i<nfiles; i++)
	{
		input_file = p.i.input_files.at(i);

		time.start();

		if (!std::filesystem::exists(input_file))
		{
			input_file = formPath<std::string>(p.g.input_folder, input_file);
			if (!std::filesystem::exists(input_file))
			{	
				std::cout << "- Error:  " << input_file << " does not exist" << std::endl;
				return EXIT_FAILURE;
			}
		}

		// Form output file path
		std::stringstream ss_outputfile("");
		ss_outputfile << std::filesystem::path(input_file).stem().generic_string();
		ss_outputfile << ".edit.xml";
		output_file = ss_outputfile.str();

		ss = leafStr<std::string>(p.i.input_files.at(i)) + "... ";
		if (!p.g.quiet)
			std::cout << "  " << std::setw(strlen_arg) << std::left << ss;
		
		// Obtain exclusive use of the input file pointer
		import_file_mutex.lock();
		
		try 
		{
			if (ifsDynaML_)
			{
				ifsDynaML_->close();
				delete ifsDynaML_;
			}

			ifsDynaML_ = new std::ifstream;

			// Open and seek to end immediately after opening.
			file_opener(ifsDynaML_, input_file, std::ios::in | std::ios::ate, ascii, true);
		
			// get file size and return to start
			sifsFileSize_ = (size_t)ifsDynaML_->tellg();
			ifsDynaML_->seekg(0, std::ios::beg);
		}
		catch (const std::ios_base::failure& f) {	
			std::cout << "- Error: An error was encountered when opening " << input_file << "." << std::endl;
			std::cout << "  Check that the file exists and that the file is not already opened." << std::endl << f.what();
			return EXIT_FAILURE;
		}
		catch (...) {	
			std::cout << "- Error: An error was encountered when opening " << input_file << "." << std::endl;
			std::cout << "  Check that the file exists and that the file is not already opened.";
			return EXIT_FAILURE;
		}

		try 
		{
			if (ofsDynaML_)
			{
				ofsDynaML_->close();
				delete ofsDynaML_;
			}

			ofsDynaML_ = new std::ofstream;

			// Open output file
			file_opener(ofsDynaML_, output_file, std::ios::out, ascii);

		}
		catch (const std::ios_base::failure& f) {	
			std::cout << "- Error: An error was encountered when opening " << output_file << "." << std::endl << f.what();
			return EXIT_FAILURE;
		}
		catch (...) {	
			std::cout << "- Error: An error was encountered when opening " << output_file << "." << std::endl;
			return EXIT_FAILURE;
		}

		try {
			measurementUpdateCount = 0;
			unsureFrameCount = 0;
			unsureEpochCount = 0;

			// Process the file
			processFile(ifsDynaML_, ofsDynaML_, defaultFrame, defaultEpoch, 
				measurementUpdateCount, unsureFrameCount, unsureEpochCount);
		}
		catch (const std::runtime_error& e) {
			std::cout << std::endl << "- Error: " << e.what() << std::endl;
			return EXIT_FAILURE;
		}catch (const XMLInteropException& e) {
			std::cout << std::endl << "- Error: " << e.what() << std::endl;
			return EXIT_FAILURE;
		}
		catch (...) 
		{
			std::cout << "+ Exception of unknown type!\n";
			return EXIT_FAILURE;
		} 

		// Finish up
		time.stop();
		elapsed_time = boost::posix_time::milliseconds(time.elapsed().wall/MILLI_TO_NANO);

		if (!p.g.quiet)
		{
			std::cout << "done." << std::endl << std::endl;
			std::cout << "+ Output file:                          " << output_file << std::endl;
			std::cout << "+ Default reference frame:              " << defaultFrame << std::endl;
			std::cout << "+ Default epoch:                        " << defaultEpoch << std::endl;
			std::cout << "+ No. DnaMeasurement records updated:   " << measurementUpdateCount << std::endl;
			std::cout << "+ No. comments with \"unsure\" frame:     " << unsureFrameCount << std::endl;
			std::cout << "+ No. comments with \"unsure\" epoch:     " << unsureEpochCount << std::endl;
			        
		}

		// Close input file
		try {
			if (ifsDynaML_)
			{
				if (ifsDynaML_->good())
					ifsDynaML_->close();
				delete ifsDynaML_;
			}
		}
		catch (const std::ifstream::failure& e)
		{
			if (ifsDynaML_->rdstate() & std::ifstream::eofbit)
			{
				ifsDynaML_ = 0;
			}
			std::cout << "- Error: An error was encountered when closing " << input_file << "." << std::endl << std::endl << e.what() << std::endl << "  Check that the file exists and that the file is not already opened.";
			return EXIT_FAILURE;
		}

		ifsDynaML_ = 0;

		// Close output file
		try {
			if (ofsDynaML_)
			{
				if (ofsDynaML_->good())
					ofsDynaML_->close();
				delete ofsDynaML_;
			}
		}
		catch (const std::ifstream::failure& e)
		{
			if (ofsDynaML_->rdstate() & std::ifstream::eofbit)
			{
				ofsDynaML_ = 0;
			}
			std::cout << "- Error: An error was encountered when closing " << output_file << "." << std::endl << std::endl << e.what() << std::endl;
			return EXIT_FAILURE;
		}

		ofsDynaML_ = 0;

		if (!p.g.quiet)
		{
			std::string time_message = formatedElapsedTime<std::string>(&elapsed_time, "+ File processing time:                 ");
			std::cout << time_message << std::endl << std::endl;
		}

		// release file pointer mutex
		import_file_mutex.unlock();

	}

	return EXIT_SUCCESS;

}




	
