//============================================================================
// Name         : dnaplot.cpp
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
// Description  : DynAdjust Plot library
//============================================================================

#include <dynadjust/dnaplot/dnaplot.hpp>

#include <include/functions/dnastrutils.hpp>
#include <include/functions/dnafilepathfuncs.hpp>
namespace dynadjust { 
namespace networkplot {

dna_plot::dna_plot()
	: datum_(DEFAULT_EPSG_U)
	, output_folder_("")
	, network_name_("")
	, stationCount_(0)
	, blockCount_(1)
{
	v_msr_file_.clear();
	v_isl_const_file_.clear();
	v_isl_pts_file_.clear();
	v_jsl_const_file_.clear();
	v_isl_lbl_file_.clear();
	v_stn_cor_file_.clear();
	v_tectonic_plate_file_.clear();

	v_parameterStationList_.clear();

	InitialiseAppsandSystemCommands();

#ifdef _MSC_VER
#if (_MSC_VER < 1900)
	{
		// this function is obsolete in MS VC++ 14.0, VS2015
		// Set scientific format to print two places for the exponent
		_set_output_format(_TWO_DIGIT_EXPONENT);
	}
#endif
#endif
}
	

dna_plot::~dna_plot()
{
	
}

void dna_plot::InitialiseAppsandSystemCommands()
{
	_APP_GMTSET_ = "gmt gmtset";
	_APP_PSCOAST_ = "gmt pscoast";
	_APP_PSCONVERT_ = "gmt psconvert";
	_APP_PSTEXT_ = "gmt pstext";
	_APP_PSVELO_ = "gmt psvelo";
	_APP_PSXY_ = "gmt psxy";
	_APP_PSLEGEND_ = "gmt pslegend";
	
	// common commands across Windows / Linux / Apple, or 
	// unique to one particular system
	_MAKEDIR_CMD_ = "mkdir ";
	_ECHO_CMD_ = "echo ";
	_CHMOD_CMD_ = "chmod +x ";
	_GMT_TMP_DIR_ = "GMT_TMPDIR";

	// system-specific variables
#if defined(_WIN32) || defined(__WIN32__)
	_LEGEND_ECHO_ = _ECHO_CMD_;
	_LEGEND_CMD_1_ = " > ";
	_LEGEND_CMD_2_ = " >> ";
	_COMMENT_PREFIX_ = ":: ";
	_CMD_EXT_ = ".bat";
	_CMD_HEADER_ = "@echo off";
	_PDF_AGGREGATE_ = "pdftk ";
	_DELETE_CMD_ = "del /Q /F ";
	_COPY_CMD_ = "copy /Y ";
	_MOVE_CMD_ = "move /Y ";
	_NULL_OUTPUT_ = " > NUL";
	_ENV_GMT_TMP_DIR_ = "%" + _GMT_TMP_DIR_ + "%";
	_MAKEENV_CMD_ = "set ";
	_RMDIR_CMD_ = "rmdir /Q /S ";

#elif defined(__linux) || defined(sun) || defined(__unix__) || defined(__APPLE__)
	_LEGEND_ECHO_ = "";
	_LEGEND_CMD_1_ = "";
	_LEGEND_CMD_2_ = "";
	_COMMENT_PREFIX_ = "# ";
	_CMD_EXT_ = ".sh";
	_CMD_HEADER_ = "#!/bin/bash";
	_PDF_AGGREGATE_ = "pdfunite ";
	_DELETE_CMD_ = "rm -rf ";
	_COPY_CMD_ = "cp ";
	_MOVE_CMD_ = "mv ";
	_NULL_OUTPUT_ = " > /dev/null";
	_MAKETEMP_CMD_ = "mktemp ";
	_ENV_GMT_TMP_DIR_ = "$" + _GMT_TMP_DIR_;
	_MAKEENV_CMD_ = "export ";
	_RMDIR_CMD_ = "rm -rf ";

#endif

}
	

void dna_plot::CreateSegmentationGraph(const plotGraphMode& graphMode)
{

	// Execute gnuplot in a separate process
	InvokeGnuplot();

	// clean up gnuplot command file and input data file
	CleanupGnuplotFiles(graphMode);
}

void dna_plot::CleanupGnuplotFiles(const plotGraphMode& graphMode)
{
	if (pprj_->p._keep_gen_files)
		return;
	
	std::stringstream ss;
	ss << _DELETE_CMD_;

	// remove gnuplot command file
	ss << pprj_->p._gnuplot_cmd_file << " ";

	// remove gnuplot input file
	switch (graphMode)
	{
	case StationsMode:
		ss << seg_stn_graph_file_;
		break;
	case MeasurementsMode:
		ss << seg_msr_graph_file_;
		break;
	default:
		break;
	}

	// delete
	std::string system_file_cmd = ss.str();
	std::system(system_file_cmd.c_str());
}

void dna_plot::CreategnuplotGraphEnvironment(project_settings* pprj, const plotGraphMode& graphMode)
{
	InitialiseAppsandSystemCommands();

	// Set up the environment
	pprj_ = pprj;

	if (!std::filesystem::exists(pprj_->g.output_folder))
	{
		std::stringstream ss("CreategnuplotGraphEnvironment(): Output path does not exist... \n\n    ");
		ss << pprj_->g.output_folder << ".";
		SignalExceptionPlot(ss.str(), 0, NULL);
	}

	output_folder_ = pprj_->g.output_folder;
	network_name_ = pprj_->g.network_name;

	/////////////////////////////////////////////////////////
	// create gnuplot command file and set gnuplot parameters 
	std::string gnuplot_cmd_filename("graph_" + network_name_);

	std::string gnuplot_pic_name;
	switch (graphMode)
	{
	case StationsMode:
		gnuplot_pic_name = (network_name_ + "_graph_stn");
		gnuplot_cmd_filename.append("_stns");
		break;
	case MeasurementsMode:
		gnuplot_pic_name = (network_name_ + "_graph_msr");
		gnuplot_cmd_filename.append("_msrs");
		break;
	}

	gnuplot_cmd_filename.append(_CMD_EXT_);
	std::string gnuplot_cmd_file(output_folder_ + FOLDER_SLASH + gnuplot_cmd_filename);

	pprj_->p._gnuplot_cmd_file = gnuplot_cmd_file;
	
	// create pdf filename
	pprj_->p._pdf_file_name = output_folder_ + FOLDER_SLASH + gnuplot_pic_name + ".pdf";
	
	switch (graphMode)
	{
	case StationsMode:
		PlotGnuplotDatFileStns();
		break;
	case MeasurementsMode:
		PlotGnuplotDatFileMsrs();
		break;
	}

	PrintGnuplotCommandFile(gnuplot_cmd_file, graphMode);

}

void dna_plot::InvokeGnuplot()
{
	// Invoke gnuplot using absolute path				
	std::string system_file_cmd = "gnuplot " + safe_absolute_path(pprj_->p._gnuplot_cmd_file);

	// set up a thread group to execute the gnuplot in parallel
	std::thread gnuplot_thread{dna_create_threaded_process(system_file_cmd)};
	
	// go!
	gnuplot_thread.join();
}


	

void dna_plot::PlotGnuplotDatFileStns()
{
	seg_stn_graph_file_ = output_folder_ + FOLDER_SLASH + network_name_ + "-stn.seg.data";
	
	std::ofstream seg_data;
	try {
		// Create gnuplot station segment data file.  Throws runtime_error on failure.
		file_opener(seg_data, seg_stn_graph_file_);
	}
	catch (const std::runtime_error& e) {
		SignalExceptionPlot(e.what(), 0, NULL);
	}

	std::stringstream ss(""), st("");
	ss << "\"Max block size (" << std::setprecision(0) << blockThreshold_ << ")\" ";
	st << "\"Min inner stns (" << minInnerStns_ << ")\" ";
	seg_data << std::setw(HEADER_18) << std::left << "Block" <<
		std::setw(HEADER_32) << std::left << ss.str() <<
		std::setw(HEADER_32) << std::left << st.str() <<
		std::setw(HEADER_25) << std::left << "\"Total block size\"" <<
		std::setw(HEADER_18) << std::left << "\"Inner stns\"" <<
		std::setw(HEADER_18) << std::left << "\"Junction stns\"" << std::endl;

	for (UINT32 block=0; block<blockCount_; ++block)
		seg_data << std::setw(HEADER_18) << std::left << block+1 << 
			std::setw(HEADER_32) << std::left << std::setprecision(0) << blockThreshold_ <<				// threshold
			std::setw(HEADER_32) << std::left << std::setprecision(0) << minInnerStns_ <<				// threshold
			std::setw(HEADER_25) << std::left << v_ISL_.at(block).size() + v_JSL_.at(block).size() <<	// block size
			std::setw(HEADER_18) << std::left << v_ISL_.at(block).size() <<		// inner station count
			std::setw(HEADER_18) << std::left << v_JSL_.at(block).size() << std::endl;	// junction station count

	seg_data.close();
}
	

void dna_plot::PlotGnuplotDatFileMsrs()
{
	seg_msr_graph_file_ = output_folder_ + FOLDER_SLASH + network_name_ + "-msr.seg.data";
	
	std::ofstream seg_data;
	try {
		// Create gnuplot measurement segment data file.  Throws runtime_error on failure.
		file_opener(seg_data, seg_msr_graph_file_);
	}
	catch (const std::runtime_error& e) {
		SignalExceptionPlot(e.what(), 0, NULL);
	}

	UINT16 c;

	// print header
	seg_data << std::setw(PRINT_VAR_PAD) << std::left << "Block";
	std::stringstream ss;
	ss << "\"Total measurements\"";
	seg_data << std::setw(PRINT_VAR_PAD) << std::left << ss.str();

	for (c=0; c<_combined_msr_list.size(); c++)
	{		
		if (parsemsrTally_.MeasurementCount(_combined_msr_list.at(c)) == 0)
			continue;
		ss.str("");
		ss << "\"" << measurement_name<char, std::string>(_combined_msr_list.at(c)) << " (" <<
			parsemsrTally_.MeasurementCount(_combined_msr_list.at(c)) << ")\"";
		seg_data << std::setw(PRINT_VAR_PAD) << std::left << ss.str();
	}
	seg_data << std::endl;

	// Tally up measurement types for each block
	ComputeMeasurementCount();

	// print measurement tally for each block
	UINT32 block;
	for (block=0; block<blockCount_; ++block)
	{
		seg_data << std::setw(PRINT_VAR_PAD) << std::left << block + 1;
		seg_data << std::setw(PRINT_VAR_PAD) << std::left << v_msr_tally_.at(block).TotalCount();

		for (c=0; c<_combined_msr_list.size(); c++)
		{
			// do any measurements of this type exist at all?
			if (parsemsrTally_.MeasurementCount(_combined_msr_list.at(c)) > 0)
			{
				// yes, so test if this block has such measurements...
				//if (v_msr_tally_.at(block).MeasurementCount(_combined_msr_list.at(c)) == 0)			// no 
				//	seg_data << std::setw(PRINT_VAR_PAD) << std::left << "-";
				//else
					seg_data << std::setw(PRINT_VAR_PAD) << std::left << v_msr_tally_.at(block).MeasurementCount(_combined_msr_list.at(c));
			}
		}
		seg_data << std::endl;
	}
	
	seg_data.close();
}
	

void dna_plot::ComputeMeasurementCount()
{
	v_msr_tally_.resize(blockCount_);

	UINT32 block;
	it_vUINT32 _it_msr;

	for (block=0; block<blockCount_; ++block)
	{
		for (_it_msr=v_CML_.at(block).begin(); _it_msr<v_CML_.at(block).end(); _it_msr++)
		{
			// Increment single station measurement counters...
			switch (bmsBinaryRecords_.at(*_it_msr).measType)
			{
			case 'A': // Horizontal angle
				v_msr_tally_.at(block).A++;
				break;
			case 'B': // Geodetic azimuth
				v_msr_tally_.at(block).B++;
				break;
			case 'C': // Chord dist
				v_msr_tally_.at(block).C++;
				break;
			case 'D': // Direction set
				if (bmsBinaryRecords_.at(*_it_msr).measStart == xMeas)
					v_msr_tally_.at(block).D += bmsBinaryRecords_.at(*_it_msr).vectorCount1;
				break;
			case 'E': // Ellipsoid arc
				v_msr_tally_.at(block).E++;
				break;
			case 'G': // GPS Baseline (treat as single-baseline cluster)
				v_msr_tally_.at(block).G ++;
				break;
			case 'X': // GPS Baseline cluster
				if (bmsBinaryRecords_.at(*_it_msr).measStart == xMeas)
				{
					if (bmsBinaryRecords_.at(*_it_msr).vectorCount1 == bmsBinaryRecords_.at(*_it_msr).vectorCount2 + 1)
						v_msr_tally_.at(block).X += (bmsBinaryRecords_.at(*_it_msr).vectorCount1 * 3);
				}					
				break;
			case 'H': // Orthometric height
				v_msr_tally_.at(block).H++;
				break;
			case 'I': // Astronomic latitude
				v_msr_tally_.at(block).I++;
				break;
			case 'J': // Astronomic longitude
				v_msr_tally_.at(block).J++;
				break;
			case 'K': // Astronomic azimuth
				v_msr_tally_.at(block).K++;
				break;
			case 'L': // Level difference
				v_msr_tally_.at(block).L++;
				break;
			case 'M': // MSL arc
				v_msr_tally_.at(block).M++;
				break;
			case 'P': // Geodetic latitude
				v_msr_tally_.at(block).P++;
				break;
			case 'Q': // Geodetic longitude
				v_msr_tally_.at(block).Q++;
				break;
			case 'R': // Ellipsoidal height
				v_msr_tally_.at(block).R++;
				break;
			case 'S': // Slope distance
				v_msr_tally_.at(block).S++;
				break;
			case 'V': // Zenith distance
				v_msr_tally_.at(block).V++;
				break;
			case 'Y': // GPS point cluster
				if (bmsBinaryRecords_.at(*_it_msr).measStart == xMeas)
				{
					if (bmsBinaryRecords_.at(*_it_msr).vectorCount1 == bmsBinaryRecords_.at(*_it_msr).vectorCount2 + 1)
						v_msr_tally_.at(block).Y += (bmsBinaryRecords_.at(*_it_msr).vectorCount1 * 3);
				}					
				break;
			case 'Z': // Vertical angle
				v_msr_tally_.at(block).Z++;
				break;
			default:
				std::stringstream ss;
				ss << "ComputeMeasurementCount(): Unknown measurement type:  " << bmsBinaryRecords_.at(*_it_msr).measType << std::endl;
				throw NetPlotException(ss.str(), 0);
			}
		}
	}
}


void dna_plot::PrintGnuplotCommandFile(const std::string& gnuplot_cmd_file, const plotGraphMode& graphMode)
{
	try {
		// Create gnuplot batch file.  Throws runtime_error on failure.
		file_opener(gnuplotbat_file_, gnuplot_cmd_file);
	}
	catch (const std::runtime_error& e) {
		SignalExceptionPlot(e.what(), 0, NULL);
	}

	//if (output_folder_ != ".")
	//	gnuplotbat_file_ << "cd '" << output_folder_ << "'" << std::endl << std::endl;

	//gnuplotbat_file_ << "set terminal postscript eps enhanced color solid colortext" << std::endl << std::endl;
	gnuplotbat_file_ << "set terminal pdf enhanced color solid linewidth 0.75" << std::endl << std::endl;
	// gnuplot requires single quotes for filenames
	gnuplotbat_file_ << "set output '" << pprj_->p._pdf_file_name << "'" << std::endl;

	// histogram style
	gnuplotbat_file_ << "set style fill transparent solid 0.4" << std::endl;
	gnuplotbat_file_ << "set boxwidth 0.5" << std::endl;

	UINT32 upperLimit(0), block(0);
	
	switch (graphMode)
	{
	case StationsMode:
		for (block=0; block<blockCount_; ++block)
			if (upperLimit < (v_ISL_.at(block).size() + v_JSL_.at(block).size()))
				upperLimit = static_cast<UINT32>(v_ISL_.at(block).size() + v_JSL_.at(block).size());
		
		upperLimit = std::max(blockThreshold_, upperLimit);

		gnuplotbat_file_ << "set title \"" << "Station segmentation summary for " << network_name_ << "\" font \"Calibri,16\" noenhanced" << std::endl << std::endl;
		gnuplotbat_file_ << "set key outside center bottom horizontal Left reverse enhanced autotitles samplen 2.5 font \"Calibri,8\"" << std::endl;
		gnuplotbat_file_ << "set key width -2 height 2.5" << std::endl << std::endl;

		break;
	case MeasurementsMode:
		for (block=0; block<blockCount_; ++block)
			if (upperLimit < v_msr_tally_.at(block).TotalCount())
				upperLimit = v_msr_tally_.at(block).TotalCount();
	
		gnuplotbat_file_ << "set title \"" << "Measurement segmentation summary for " << network_name_ << "\" font \"Calibri,20\" noenhanced" << std::endl << std::endl;
		//gnuplotbat_file_ << "set key outside right top vertical Left reverse enhanced autotitles columnhead box samplen 2.5 font \"Calibri,8\"" << std::endl;
		gnuplotbat_file_ << "set key outside center bottom horizontal Left reverse enhanced autotitles columnhead samplen 2.5 font \"Calibri,8\"" << std::endl;
		//gnuplotbat_file_ << "set key width -15 height 0" << std::endl << std::endl;
		gnuplotbat_file_ << "set key width -2 height 2.5" << std::endl << std::endl;

		gnuplotbat_file_ << "set style histogram rowstacked title offset character 0, 0, 0" << std::endl;
		gnuplotbat_file_ << "set style data histograms" << std::endl;
		gnuplotbat_file_ << "set datafile missing '-'" << std::endl;
		gnuplotbat_file_ << std::endl;
		break;
	}

	upperLimit = upperLimit + upperLimit / 10;

	gnuplotbat_file_ << "set format x '%.0f'" << std::endl;
	gnuplotbat_file_ << "set format y '%.0f'" << std::endl;
	gnuplotbat_file_ << "set yrange[0:" << upperLimit << "]" << std::endl;
	gnuplotbat_file_ << "set auto x" << std::endl;
	gnuplotbat_file_ << "set ytics scale 0.25 font \"Calibri,8\"" << std::endl;

	gnuplotbat_file_ << "set xtics scale 0.25 nomirror" << std::endl;

	UINT32 fontSize(8);
	if (blockCount_ > 5000)
	{
		fontSize = 5;
		gnuplotbat_file_ << "set xtics 0,500 font \"Calibri,6\"" << std::endl << std::endl;
	}
	else if (blockCount_ > 1000)
	{
		fontSize = 5;
		gnuplotbat_file_ << "set xtics 0,100 font \"Calibri,6\"" << std::endl << std::endl;
	}
	else if (blockCount_ > 500)
	{
		fontSize = 5;
		gnuplotbat_file_ << "set xtics 0,50 font \"Calibri,5\"" << std::endl << std::endl;
	}
	else if (blockCount_ > 100)
	{
		fontSize = 5;
		gnuplotbat_file_ << "set xtics 0,10 font \"Calibri,6\"" << std::endl << std::endl;
	}
	else if (blockCount_ > 50)
	{
		fontSize = 6;
		gnuplotbat_file_ << "set xtics 0,5 font \"Calibri,8\"" << std::endl << std::endl;
	}
	else
		gnuplotbat_file_ << "set xtics 0,1 font \"Calibri,8\"" << std::endl << std::endl;
	

	// x-axis label
	std::stringstream ss("");
	ss << "Segmented Network Blocks (Total " << std::fixed << std::setprecision(0) << blockCount_ << ")";
	gnuplotbat_file_ << "set xlabel '" << ss.str() << "' font \"Calibri,10\"" << std::endl;

	switch (graphMode)
	{
	case StationsMode:
		PrintGnuplotCommandFileStns(fontSize);
		break;
	case MeasurementsMode:
		PrintGnuplotCommandFileMsrs(fontSize);
		break;
	}

	gnuplotbat_file_.close();
}
	
void dna_plot::PrintGnuplotCommandFileStns(const UINT32& fontSize)
{
	std::stringstream ss("");
	ss << "Station Count (Total " << std::fixed << std::setprecision(0) << stationCount_ << ")";
	gnuplotbat_file_ << "set ylabel '" << ss.str() << "' font \"Calibri,10\"" << std::endl << std::endl;

	// All colours based on a palette:
	//   https://coolors.co/ffd275-235789-da5552-43aa8b-39a9db
	gnuplotbat_file_ << "set style line 1 lw 0.75 lt 1 pt 7 ps 0.25 lc rgb \"#35A7FF\"         # total block size" << std::endl;		// royalblue
	gnuplotbat_file_ << "set style line 2 lw 2 lt 5 pt 7 ps 0.25 lc rgb \"#43AA8B\"            # threshold" << std::endl;			// zomp (green)
	gnuplotbat_file_ << "set style line 3 lw 2 lt 5 pt 7 ps 0.25 lc rgb \"#FFD275\"            # minimum inner size" << std::endl;	// orange yellow crayola
	gnuplotbat_file_ << "set style line 4 lw 2 lt 1 pt 7 ps 0.25 lc rgb \"#235789\"            # inners" << std::endl;				// bdazzled blue
	gnuplotbat_file_ << "set style line 5 lw 2 lt 1 pt 7 ps 0.25 lc rgb \"#DA5552\"            # junctions" << std::endl << std::endl;	// indian red

	gnuplotbat_file_ << "plot '" << seg_stn_graph_file_ << "' using 1:4 with boxes ls 1 title columnheader(4), \\" << std::endl;
	gnuplotbat_file_ << "     '" << seg_stn_graph_file_ << "' using 1:4:(sprintf(\"%.0f\",$4)) with labels font \"Calibri," << 
		fontSize << "\" center offset 0,0.5 notitle, \\" << std::endl;
	gnuplotbat_file_ << "     '" << seg_stn_graph_file_ << "' using 1:2 with lines ls 2 title columnheader(2), \\" << std::endl;
	gnuplotbat_file_ << "     '" << seg_stn_graph_file_ << "' using 1:3 with lines ls 3 title columnheader(3), \\" << std::endl;
	gnuplotbat_file_ << "     '" << seg_stn_graph_file_ << "' using 1:5 with linespoints ls 4 title columnheader(5), \\" << std::endl;
	gnuplotbat_file_ << "     '" << seg_stn_graph_file_ << "' using 1:5:(sprintf(\"%.0f\",$5)) with labels tc ls 4 font \"Calibri," << 
		fontSize << "\" center offset 1,0 notitle, \\" << std::endl;
	gnuplotbat_file_ << "     '" << seg_stn_graph_file_ << "' using 1:6 with linespoints ls 5 title columnheader(6), \\" << std::endl;
	gnuplotbat_file_ << "     '" << seg_stn_graph_file_ << "' using 1:6:(sprintf(\"%.0f\",$6)) with labels tc ls 5 font \"Calibri," << 
		fontSize << "\" center offset 1,0 notitle" << std::endl << std::endl;
}
	

void dna_plot::PrintGnuplotCommandFileMsrs(const UINT32& fontSize)
{
	std::stringstream ss("");
	ss << "Measurement Count (Total " << std::fixed << std::setprecision(0) << measurementCount_ << ")";
	gnuplotbat_file_ << "set ylabel '" << ss.str() << "' font \"Calibri,10\"" << std::endl << std::endl;

	// All colours based on a palette:
	//   https://coolors.co/ffd275-235789-da5552-43aa8b-39a9db

	UINT32 line(1);
	ss.str("");
	ss << "\"#4169e1\"";
	gnuplotbat_file_ << "set style line " << line++ << " lw 1 lt 1 pt 7 ps 0.5 lc rgb " << std::left << std::setw(PRINT_VAR_PAD) << ss.str() << " # total block size" << std::endl;	// royalblue

	// print measurements for each block
	UINT32 c;
	std::string colour;
	it_pair_string _it_colour;

	std::sort(pprj_->p._msr_colours.begin(), pprj_->p._msr_colours.end(), ComparePairFirst<std::string>());

	for (c=0; c<_combined_msr_list.size(); c++)
	{
		if (parsemsrTally_.MeasurementCount(_combined_msr_list.at(c)) < 1)
			continue;
		colour = _combined_msr_list.at(c);
		_it_colour = equal_range(pprj_->p._msr_colours.begin(), pprj_->p._msr_colours.end(), 
			colour, ComparePairFirst<std::string>());

		if (_it_colour.first == _it_colour.second)
			colour = "light-gray";
		else
			colour = _it_colour.first->second;
		
		ss.str("");
		ss << "\"" << colour << "\"";
		gnuplotbat_file_ << "set style line " << line++ << " lw 1 lt 1 pt 7 ps 0.5 lc rgb " << std::left << std::setw(PRINT_VAR_PAD) << ss.str() << 
			" # \"" << measurement_name<char, std::string>(_combined_msr_list.at(c)) << "\"" << std::endl;
		
	}
	gnuplotbat_file_ << std::endl;
	
	//UINT32 block(0);
	//if (v_msr_tally_.at(block).MeasurementCount(_combined_msr_list.at(c)) == 0)
	//	continue;

	// plot total measurement count
	line = 3;
	UINT32 linestyle(line-1), msrs(0);
	gnuplotbat_file_ << "plot '";
	measurementCategories_ = 0;
	
	for (c=0; c<_combined_msr_list.size(); c++)
	{
		if (parsemsrTally_.MeasurementCount(_combined_msr_list.at(c)) < 1)
			continue;
		if (msrs++ > 0)
			gnuplotbat_file_ << ", \\" << std::endl << "     '";
		gnuplotbat_file_ << seg_msr_graph_file_ << "' using " << line++ << ":xtic(1) ls " << linestyle++;
		measurementCategories_++;
	}

	gnuplotbat_file_ << ", \\" << std::endl << "     '" << seg_msr_graph_file_ << 
		"' using 0:2:(sprintf(\"%d\",$2)) with labels font \"Calibri," << 
			fontSize << "\" center offset 0,0.5 notitle" << std::endl;
	
}


void dna_plot::InitialiseGMTParameters()
{
	// Set initial parameters
	if (!std::filesystem::exists(pprj_->g.output_folder))
	{
		std::stringstream ss("InitialiseGMTParameters(): Output path does not exist... \n\n    ");
		ss << pprj_->g.output_folder << ".";
		SignalExceptionPlot(ss.str(), 0, NULL);
	}

	InitialiseAppsandSystemCommands();

	output_folder_ = pprj_->g.output_folder;
	network_name_ = pprj_->g.network_name;

	v_isl_const_file_.clear();
	v_jsl_const_file_.clear();
	v_isl_pts_file_.clear();
	v_isl_lbl_file_.clear();
	v_stn_cor_file_.clear();
	v_tectonic_plate_file_.clear();
	v_stn_err_file_.clear();
	v_stn_apu_file_.clear();
	v_jsl_pts_file_.clear();
	v_jsl_lbl_file_.clear();

	default_paper_width_ = 59.4;
	default_paper_height_ = 42.0;

	plotConstraints_ = false;

	lowerDeg_ = 1000.0;
	leftDeg_ = 1000.0;
	upperDeg_ = -1000.0;
	rightDeg_ = -1000.0;

	default_limits_ = true;
	
	// check for valid plot limits
	if (!pprj_->p._bounding_box.empty())
		default_limits_ = false;
	else if (pprj_->p._plot_area_radius > 0. && 
		pprj_->p._plot_centre_latitude > -90.00000001 && 
		pprj_->p._plot_centre_latitude < 90.00000001 && 
		pprj_->p._plot_centre_longitude > -180.00000001 && 
		pprj_->p._plot_centre_longitude < 180.00000001)
			default_limits_ = false;
	else if (pprj_->p._plot_area_radius > 0. && 
		!pprj_->p._plot_station_centre.empty())
			default_limits_ = false;
}


void dna_plot::FinaliseGMTParameters()
{
	//
	// At this point, all limit values are in decimal degrees
	//

	if (rightDeg_ < leftDeg_)
	{
		std::stringstream ss;
		ss << "Right limit cannot be less than left limit." << std::endl;
		throw NetPlotException(ss.str(), 0);
	}

	if (upperDeg_ < lowerDeg_)
	{
		std::stringstream ss;
		ss << "Upper limit cannot be less than lower limit." << std::endl;
		throw NetPlotException(ss.str(), 0);
	}

	// calculate dimensions of data extents
	dWidth_ = (rightDeg_ - leftDeg_);
	dHeight_ = (upperDeg_ - lowerDeg_);
	
	// if data is a single point (or a very small area)
	if (dHeight_ < 0.00005)
		dHeight_ += seconds15;
	if (dWidth_ < 0.00005)
		dWidth_ += seconds15;	
	
	// capture smallest dimension
	dDimension_ = (std::min(dWidth_, dHeight_));
	
	// Determine a buffer to envelope the entire plot, set to
	// 10% of the width/height (whichever is smaller)
	dBuffer_ = (dDimension_ * (0.1));

	// calculate latitude at which to place the scale bar
	dScaleLat_ = (lowerDeg_ - dBuffer_ / 2.0);

	// calculate centre point of data
	centre_width_ = (rightDeg_ + leftDeg_) / 2.0;
	centre_height_ = (upperDeg_ + lowerDeg_) / 2.0;

	// Has the user specified certain limits, or will the data be used to 
	// define the limits
	if (default_limits_)
	{
		// OK, the spatial extent of the data sets the limits, so
		// add a buffer accordingly
		upperDeg_ += dBuffer_;
		lowerDeg_ -= dBuffer_;
		leftDeg_ -= dBuffer_;
		rightDeg_ += dBuffer_;

		if (dDimension_ > seconds60)
		{
			// round down to nearest 5"
			leftDeg_ = leftDeg_ - fmod(leftDeg_, seconds05);
			lowerDeg_ = lowerDeg_ - fmod(lowerDeg_, seconds05);
			rightDeg_ = rightDeg_ + (seconds05 - fmod(rightDeg_, seconds05));
			upperDeg_ = upperDeg_ + (seconds05 - fmod(upperDeg_, seconds05));
		}
	}
	else
	{
		// Limits have been calculated via user input, so simply calculate
		// calculate latitude at which the scale bar
		// will be placed
		centre_width_ = pprj_->p._plot_centre_longitude;
		centre_height_ = pprj_->p._plot_centre_latitude;

		dScaleLat_ = lowerDeg_ + dBuffer_ / 2.0;
	}

	// prevent limit values from exceeding 'normal' limits
	if (default_limits_)
	{
		if (lowerDeg_ < -90.)
			lowerDeg_ = -90.;
		if (lowerDeg_ > 90.)
			lowerDeg_ = 90.;
	}

	// set the position of the error ellipse legend
	uncertainty_legend_long_ = leftDeg_ + dBuffer_;
	uncertainty_legend_lat_ = dScaleLat_;

	if (pprj_->p._plot_correction_arrows)
	{
		// set the position of the correction arrow legend
		// put the corrections legend on the right hand side
		arrow_legend_long_ = rightDeg_ - dBuffer_ * 2.0;
		arrow_legend_lat_ = dScaleLat_;
	}

	switch (pprj_->p._projection)
	{
	case world:
	case orthographic:
	case robinson:
		default_paper_width_ = 27.7;
		default_paper_height_ = 20.1;
		break;
	default:
		default_paper_width_ = 59.4;
		default_paper_height_ = 42.0;
	}	

	// Portrait or Landscape?
	// A3 paper width (landscape) is 42cm, and (portrait) is 29.7cm.  Less 
	// nomenclature (~2 cm), this leaves an available width for plotting 
	// of 40cm L or 27.7cm P
	//
	title_block_height_ = 6;
	page_width_ = default_paper_width_;		// centimetres

	if (dWidth_ > dHeight_)
		avg_data_scale_ = dHeight_ / default_paper_height_;
	else
		avg_data_scale_ = dHeight_ / default_paper_width_;
	
	isLandscape = true;
	if (dWidth_ < (dHeight_ + (title_block_height_ * avg_data_scale_)) || !default_limits_)
	{
		page_width_ = default_paper_height_;			// centimetres
		isLandscape = false;
	}	

	pprj_->p._page_width = page_width_ + 2;

	// Scale (for A3 page width)
	if ((rightDeg_ - leftDeg_) > Degrees(PI) ||			// wider than 180 degrees?
		(upperDeg_ - lowerDeg_) > Degrees(PI_135))		// taller than 135 degrees?
	{
		// Set the appropriate map projection
		if(!pprj_->p._user_defined_projection)
		{
			if ((upperDeg_ - lowerDeg_) > Degrees(PI_135))		// taller than 135 degrees?
			{
				// Print the entire world on a flat sheet
				pprj_->p._projection = world;

				// Calculate ground width, which in this case is the circumference of the world
				dDimension_ = TWO_PI * datum_.GetEllipsoid().GetSemiMajor();		// set distance circumference of a circle
				
				// Centre world map according to data centre
				leftDeg_ = centre_width_ - 180.0;
				rightDeg_ = centre_width_ + 180.0;

				// Normalise limits according to 180 boundary
				if (leftDeg_ < -180.0)
					leftDeg_ += 180.0;
				if (rightDeg_ > 180.0)
					rightDeg_ -= 180.0;
				
				// Set latitude limits
				lowerDeg_ = -90;
				upperDeg_ = 90;
			}
			else
			{
				// Print a globe plot.  Limits are not required
				pprj_->p._projection = orthographic;
				
				// Calculate ground width between west and east limits (including buffer)
				dDimension_ = Radians(rightDeg_ - leftDeg_) / TWO_PI * datum_.GetEllipsoid().GetSemiMajor();
			}
		}
		else
			// Calculate ground width between west and east limits (including buffer)
			dDimension_ = Radians(rightDeg_ - leftDeg_) / TWO_PI * datum_.GetEllipsoid().GetSemiMajor();

		// Calculate scale
		scale_ = dDimension_ / page_width_ * 100;		// metres
	}
	else
	{
		double azimuth;

		// Calculate accurate map width from limits
		dDimension_ = RobbinsReverse<double>(		// calculate distance (in metres) of map width
			Radians(centre_height_), Radians(leftDeg_), Radians(centre_height_), Radians(rightDeg_),
			&azimuth, datum_.GetEllipsoidRef());

		// Calculate scale
		scale_ = dDimension_ / page_width_ * 100;		// metres

		if (!pprj_->p._user_defined_projection)
		{
			// default (large regions)
			pprj_->p._projection = mercator;

			// very high latitudes
			if (fabs(centre_height_) > 80.)
				pprj_->p._projection = orthographic;
			//
			// high latitudes
			else if (fabs(centre_height_) > 60.)
				pprj_->p._projection = albersConic;
			// sparse latitude coverage
			else if ((upperDeg_ - lowerDeg_) > Degrees(QUART_PI))	// taller than 45 degrees?
			{
				// tall, narrow plots
				if ((rightDeg_ - leftDeg_) < Degrees(QUART_PI))		// not more than 45 degrees wide?
					pprj_->p._projection = orthographic;
				// reasonably large area
				else
					pprj_->p._projection = mercator;
			}
			//
			// Smallish areas
			else if (scale_ < 750000 && scale_ > 5000)
				pprj_->p._projection = transverseMercator;
			//
			// wider than 20 degrees and taller than 20 degrees?
			else if ((rightDeg_ - leftDeg_) > Degrees(PI_20) && (upperDeg_ - lowerDeg_) > Degrees(PI_20))
				pprj_->p._projection = stereographicConformal;
			// 
			// wide plots
			else if ((rightDeg_ - leftDeg_) > Degrees(THIRD_PI))	// wider than 90 degrees?
				
				pprj_->p._projection = lambertEqualArea;
		}
	}

	// Update calling app's parameters
	pprj_->p._projection = pprj_->p._projection;
	pprj_->p._ground_width = dDimension_;

	// Normalise scale
	NormaliseScale(scale_);
	mapScale_ = scale_;

	// Calculate scale bar, then round
	scale_bar_width_ = dDimension_ / 3000.;		// convert to kilometres
	NormaliseScaleBar(scale_bar_width_);

	// Calculate best graticule width, then round to the best integer
	graticule_width_ = ((rightDeg_ - leftDeg_) / 4.);
	graticule_width_precision_ = 12;
	NormaliseGraticule(graticule_width_, graticule_width_precision_);

	scale_precision_ = 0;
	if (scale_bar_width_ < 1.0)
		scale_precision_ = 4;

	line_width_ = projectSettings_.p._msr_line_width;
	circle_radius_ = 0.2;
	circle_radius_2_ = circle_radius_;
	circle_line_width_ = (projectSettings_.p._msr_line_width * 1.5);

	//if (scale >= 900000)
	//{
	//	//line_width = 0.05;
	//	circle_radius = 0.05;
	//	circle_line_width = 0.05;
	//}

	// Determine coastline resolution
	coastResolution_ = "l";
	SelectCoastlineResolution(dDimension_, coastResolution_, &pprj_->p);
}


bool dna_plot::InitialiseandValidateStartingBlock(UINT32& block)
{
	block = 0;
	if (pprj_->p._plot_block_number > 0)
	{
		block = pprj_->p._plot_block_number - 1;
		return true;
	}
	return false;
}

void dna_plot::CreateExtraInputFiles()
{
	UINT32 block(0);
	bool oneBlockOnly(true);
		
	// one and only block?
	oneBlockOnly = InitialiseandValidateStartingBlock(block);

	// Now print stations labels (based on font size determined by PrintGMTPlotBatfile)
	if (pprj_->p._plot_station_labels)
	{
		if (plotBlocks_)
		{
			if (oneBlockOnly)
				PrintStationLabelsBlock(block);
			else
				for (block=0; block<blockCount_; ++block)
					PrintStationLabelsBlock(block);
		}
		else
			PrintStationLabels();
	}

	if (pprj_->p._plot_positional_uncertainty)
	{
		if (oneBlockOnly)
			PrintPositionalUncertainty(block);
		else
			for (block=0; block<blockCount_; ++block)
				PrintPositionalUncertainty(block);
	}

	if (pprj_->p._plot_error_ellipses)
	{
		if (oneBlockOnly)
			PrintErrorEllipses(block);
		else
			for (block=0; block<blockCount_; ++block)
				PrintErrorEllipses(block);
	}

	if (pprj_->p._plot_correction_arrows)
	{
		if (oneBlockOnly)
			PrintCorrectionArrows(block);
		else
			for (block=0; block<blockCount_; ++block)
				PrintCorrectionArrows(block);
	}

	if (pprj_->p._plot_plate_boundaries)
	{
		if (oneBlockOnly)
			PrintPlateBoundaries(block);
		else
			for (block=0; block<blockCount_; ++block)
				PrintPlateBoundaries(block);
	}
}

void dna_plot::CreateGMTInputFiles()
{	
	if (plotBlocks_)
	{
		UINT32 block(0);
		bool oneBlockOnly(true);
		
		// one and only block?
		if ((oneBlockOnly = InitialiseandValidateStartingBlock(block)))
			block = 0;
			
		for (; block<blockCount_; ++block)
		{
			PrintStationsDataFileBlock(block);
			if (!pprj_->p._omit_measurements)
				PrintMeasurementsDatFilesBlock(block);
			if (oneBlockOnly)
				break;
		}
	}
	else
	{
		v_msr_file_.resize(1);
		// PrintStationsDataFile calculates plot limits if user has 
		// specified centre station or lat/long, or bounding box
		PrintStationsDataFile();
		if (!pprj_->p._omit_measurements)
			PrintMeasurementsDatFiles();
	}
}

void dna_plot::InitialiseGMTFilenames()
{
	// create bat file name(s), one per block, as:
	//   create_<network_name>_block_#.[bat|sh]
	//
	// Simultaneous plot cmd files are created as:
	//   create_<network_name>_block_0.[bat|sh]
	//
	// Specific block (n) plot cmd files are created as:
	//   create_<network_name>_block_n.[bat|sh]
	v_gmt_cmd_filenames_.clear();
	v_gmt_pdf_filenames_.clear();
	std::string gmt_filename, gmt_cmd_basename("create_" + network_name_ + "_block_");
		
	UINT32 block;
	bool oneBlockOnly = InitialiseandValidateStartingBlock(block);

	for (; block<blockCount_; ++block)
	{
		// gmt command filename
		gmt_filename = gmt_cmd_basename;
		gmt_filename.append(StringFromT(block)).append(_CMD_EXT_);

		// Add folder path
		gmt_filename = pprj_->g.output_folder + FOLDER_SLASH + gmt_filename;

		// Create absolute path				
		gmt_filename = safe_absolute_path(gmt_filename);

		// Add to the list
		v_gmt_cmd_filenames_.push_back(gmt_filename);

		// generated pdf filename
		gmt_filename = network_name_;
		gmt_filename.append("_block_").append(StringFromT(block)).append(".pdf");
		v_gmt_pdf_filenames_.push_back(output_folder_ + FOLDER_SLASH + gmt_filename);

		// break out for specific block (n) plots
		if (oneBlockOnly)
			break;
	}

	pprj_->p._gmt_cmd_file = output_folder_ + FOLDER_SLASH + leafStr(v_gmt_cmd_filenames_.at(0));

	pprj_->p._pdf_file_name = output_folder_ + FOLDER_SLASH + network_name_;

	if (plotBlocks_)
		pprj_->p._pdf_file_name.append("-phased");
	else
		pprj_->p._pdf_file_name.append("-simult");

	if (plotBlocks_ && oneBlockOnly)
		pprj_->p._pdf_file_name.append("-block-").append(StringFromT(block));

	pprj_->p._pdf_file_name.append(".pdf");
	
	// simultaneous
	if (!plotBlocks_)
		v_gmt_pdf_filenames_.at(0) = pprj_->p._pdf_file_name;
}

void dna_plot::CreateGMTCommandFiles()
{
	UINT32 block(0);
	bool oneBlockOnly(true);
	
	if (plotBlocks_)
	{
		// one and only block?
		if ((oneBlockOnly = InitialiseandValidateStartingBlock(block)))
			block = 0;
	}

	for (; block<blockCount_; ++block)
	{
		try {
			// Create GMT batch file.  Throws runtime_error on failure.
			file_opener(gmtbat_file_, v_gmt_cmd_filenames_.at(block));
		}
		catch (const std::runtime_error& e) {
			SignalExceptionPlot(e.what(), 0, NULL);
		}

		// set header
		gmtbat_file_ << _CMD_HEADER_ << std::endl;
		
		// GMT bat file is printed last to reflect the options and dimensions as determined
		// by PrintStationsDataFile and PrintMeasurementsDatFiles
		CreateGMTCommandFile(block);

		// close the file
		gmtbat_file_.close();

		// change file permission to executable
#if defined(__linux) || defined(sun) || defined(__unix__) || defined(__APPLE__)		
		std::string system_file_cmd = _CHMOD_CMD_ + v_gmt_cmd_filenames_.at(block);
		std::system(system_file_cmd.c_str());
#endif

		// break out for specific block (n) plots
		if (oneBlockOnly)
			break;
	}
}


void dna_plot::CreateGMTCommandFile(const UINT32& block)
{
	it_pair_string _it_colour;
	string_string_pair this_msr;

	float annot_font_size_primary(9), label_font_size(11.);
	double symbol_offset(1.0);
	double label_offset(1.5);
	double error_ellipse_scale(uncertainty_legend_length_/largest_uncertainty_);
	UINT32 uncertainty_legend_precision(4);

	std::string psTempFile("tmp-");
	psTempFile.append(StringFromT(block)).append(".ps");
	std::string pdfTempFile("tmp-");
	pdfTempFile.append(StringFromT(block)).append(".pdf");

	std::string legendTempFile("map-block");
	legendTempFile.append(StringFromT(block)).append(".legend");

	// make temporary directory for gmt.conf gmt.history to prevent corruption of
	// gmt.conf and gmt.history during parallel processing
	gmtbat_file_ << std::endl << _COMMENT_PREFIX_ << "Create temporary folder for gmt.conf and gmt.history" << std::endl;
	
#if defined(_WIN32) || defined(__WIN32__)
	// set, create and provide access to temporary folder for gmt.conf file
	gmtbat_file_ << _MAKEENV_CMD_ << _GMT_TMP_DIR_ << "=%TEMP%\\gmt.block-" << block << std::endl;
	gmtbat_file_ << _MAKEDIR_CMD_ << _ENV_GMT_TMP_DIR_ << std::endl;
	//gmtbat_file_ << "icacls " << _ENV_GMT_TMP_DIR_ << " /grant Everyone:(f)" << std::endl;

#elif defined(__linux) || defined(sun) || defined(__unix__) || defined(__APPLE__)
	// temporary gmt.conf file location
	gmtbat_file_ << _MAKEENV_CMD_ << _GMT_TMP_DIR_ << "=$(" << _MAKETEMP_CMD_ << "-d ${TMPDIR:-/tmp}/gmt.XXXXXX)" << std::endl << std::endl;
#endif

	UINT32 colours(0), columns(5);
	
	if (plotBlocks_)
		gmtbat_file_ << std::endl << 
			_COMMENT_PREFIX_ << "GMT command file for segmented network block " << (block + 1) << std::endl;		
	else
		gmtbat_file_ << std::endl << 
			_COMMENT_PREFIX_ << "GMT command file for simultaneous network" << std::endl << std::endl;

	// write GMT parameters
	PrintGMTParameters();

	//if (isLandscape)
	//	gmtbat_file_ << _APP_GMTSET_ << " PS_PAGE_ORIENTATION landscape" << std::endl << std::endl;
	//else
	//	gmtbat_file_ << _APP_GMTSET_ << " PS_PAGE_ORIENTATION portrait" << std::endl << std::endl;
	
	gmtbat_file_ << _APP_GMTSET_ << " FONT_ANNOT " << std::fixed << std::setprecision(1) << annot_font_size_primary << "p" << std::endl;
	gmtbat_file_ << _APP_GMTSET_ << " FONT_ANNOT_PRIMARY " << std::fixed << std::setprecision(1) << annot_font_size_primary << "p" << std::endl;
	gmtbat_file_ << _APP_GMTSET_ << " FONT_LABEL " << std::fixed << std::setprecision(1) << label_font_size << "p" << std::endl;
	if (graticule_width_ < 60./3600.)
		gmtbat_file_ << _APP_GMTSET_ << " FORMAT_GEO_MAP ddd:mm:ss.x" << std::endl << std::endl;
	else
		gmtbat_file_ << _APP_GMTSET_ << " FORMAT_GEO_MAP ddd:mm" << std::endl << std::endl;

	page_width_ = pprj_->p._page_width;

	switch (pprj_->p._projection)
	{
	case world:
		// World map centered on the dateline
		// pscoast -Rg -JQ4.5i -B60f30g30 -Dc -A5000 -Gblack -P > tmp.ps
		
		// Override coastline resolution
		coastResolution_ = "i";
		pprj_->p._coasline_resolution = intermediate;
		circle_radius_ = 0.02;
		circle_line_width_ = 0.02;
		if (pprj_->p._label_font_size < 0.)
			pprj_->p._label_font_size = 5;

		gmtbat_file_ << _APP_PSCOAST_ << " -Rg" << 
			// Carree Cylindrical equidistant projection, which looks the nicest
			" -JQ" << page_width_  << "c -B60g30 -D" << coastResolution_ << 
			" -A10000 -W0.75p,16/169/243 -G245/245/245 -K > " << psTempFile << std::endl;
			//
			// Miller's Cylindrical projection, which is neither equal nor conformal. All meridians and parallels are straight lines.
			// " -JJ" << page_width  << "c -B60g30 -D" << coastResolution << " -A10000 -W0.75p,16/169/243 -G245/245/245 -P -K > " << psTempFile << std::endl;
			//
			// Cylindrical equal-area projection
			//" -JY" << page_width  << "c -B60g30 -D" << coastResolution << " -A10000 -W0.75p,16/169/243 -G245/245/245 -P -K > " << psTempFile << std::endl;
		break;
	//
	// Orthographic projection
	case orthographic:
			
		// Override coastline resolution
		coastResolution_ = "h";
		pprj_->p._coasline_resolution = high;
		circle_radius_ = 0.05;
		circle_line_width_ = 0.05;
		if (pprj_->p._label_font_size < 0.)
			pprj_->p._label_font_size = 5;

		gmtbat_file_ << _APP_PSCOAST_ << " -Rg -JG" << 
			std::fixed << std::setprecision(7) << centre_width_ << "/" <<		// longitude
			std::fixed << std::setprecision(7) << centre_height_ << "/" <<		// latitude
			std::fixed << std::setprecision(1) << page_width_ << "c -B30g15 -D" << coastResolution_ << 
			" -A10000 -W0.75p,16/169/243 -G245/245/245 -P -K > " << psTempFile << std::endl;
		break;
	//
	// Mercator projection
	case mercator:
		if (pprj_->p._label_font_size < 0.)
			pprj_->p._label_font_size = 6;

		gmtbat_file_ << _APP_PSCOAST_ << " -R" << 
		std::fixed << std::setprecision(7) << leftDeg_ << "/" <<
		std::fixed << std::setprecision(7) << lowerDeg_ << "/" <<
		std::fixed << std::setprecision(7) << rightDeg_ << "/" <<
		std::fixed << std::setprecision(7) << upperDeg_;

		// example: -Jm1.2e-2i
		gmtbat_file_ << "r -JM" << std::fixed << std::setprecision(1) << page_width_ << "c";
		gmtbat_file_ << " -D" << coastResolution_ << " -N2/0.25p -W0.75p,16/169/243 -G255/255/255 -S233/246/255 -Lf" <<
			std::fixed << std::setprecision(5) << centre_width_ << "/" << std::fixed << std::setprecision(5) << dScaleLat_ << "/" << 
			std::fixed << std::setprecision(5) << centre_height_ << "/" << std::fixed << std::setprecision(scale_precision_) << scale_bar_width_ << "k+lKilometres+jt ";
        if (graticule_width_ < 0.0001)
            gmtbat_file_ << "-Bafg";
        else
            gmtbat_file_ << "-B" <<
			std::fixed << std::setprecision(graticule_width_precision_) << graticule_width_ << "g" << std::fixed << std::setprecision(graticule_width_precision_) << graticule_width_ << "/" << 
			std::fixed << std::setprecision(graticule_width_precision_) << graticule_width_ << "g" << std::fixed << std::setprecision(graticule_width_precision_) << graticule_width_;

		if (!isLandscape)
			gmtbat_file_ << " -P";
		gmtbat_file_ << " -K > " << psTempFile << std::endl;
		
		break;
	//
	// Transverse Mercator projection
	case transverseMercator:
		if (pprj_->p._label_font_size < 0.)
			pprj_->p._label_font_size = 6;

		gmtbat_file_ << _APP_PSCOAST_ << " -R" << 
		std::fixed << std::setprecision(7) << leftDeg_ << "/" <<
		std::fixed << std::setprecision(7) << lowerDeg_ << "/" <<
		std::fixed << std::setprecision(7) << rightDeg_ << "/" <<
		std::fixed << std::setprecision(7) << upperDeg_;

		// example: -Jt139.9944444/-24.1486111/1:1000000
		//gmtbat_file_ << "r -Jt" << std::fixed << std::setprecision(7) << centre_width_ << "/" <<
		//	std::fixed << std::setprecision(7) << centre_height_ << "/1:" << std::fixed << std::setprecision(0) << scale;
		gmtbat_file_ << "r -JT" << std::fixed << std::setprecision(7) << centre_width_ << "/" <<
			std::fixed << std::setprecision(1) << page_width_ << "c";

		gmtbat_file_ << " -D" << coastResolution_ << " -N2/0.25p -W0.75p,16/169/243 -G255/255/255 -S233/246/255 -Lf" <<
			std::fixed << std::setprecision(5) << centre_width_ << "/" << std::fixed << std::setprecision(5) << dScaleLat_ << "/" << 
			std::fixed << std::setprecision(5) << centre_height_ << "/" << std::fixed << std::setprecision(scale_precision_) << scale_bar_width_ << "k+lKilometres+jt ";
		if (graticule_width_ < 0.0001)
            gmtbat_file_ << "-Bafg";
        else
            gmtbat_file_ << "-B" <<
			std::fixed << std::setprecision(graticule_width_precision_) << graticule_width_ << "g" << std::fixed << std::setprecision(graticule_width_precision_) << graticule_width_ << "/" << 
			std::fixed << std::setprecision(graticule_width_precision_) << graticule_width_ << "g" << std::fixed << std::setprecision(graticule_width_precision_) << graticule_width_;

		if (!isLandscape)
			gmtbat_file_ << " -P";
		gmtbat_file_ << " -K > " << psTempFile << std::endl;
		
		break;
	//
	// Albers conic equal-area projection
	case albersConic:
		if (pprj_->p._label_font_size < 0.)
			pprj_->p._label_font_size = 6;

		gmtbat_file_ << _APP_PSCOAST_ << " -R" << 
			std::fixed << std::setprecision(7) << leftDeg_ << "/" <<
			std::fixed << std::setprecision(7) << lowerDeg_ << "/" <<
			std::fixed << std::setprecision(7) << rightDeg_ << "/" <<
			std::fixed << std::setprecision(7) << upperDeg_;

		// example: -Jb136.5/-36/-18/-36/1:45000000
		//gmtbat_file_ << "r -Jb" << std::fixed << std::setprecision(7) << centre_width_ << "/" <<
		//	std::fixed << std::setprecision(7) << centre_height_ << "/" <<
		//	std::fixed << std::setprecision(7) << upperDeg_ - fabs(dHeight/3.) << "/" <<
		//	std::fixed << std::setprecision(7) << lowerDeg_ + fabs(dHeight/3.) << "/" <<				
		//	"1:" << std::fixed << std::setprecision(0) << scale;
		gmtbat_file_ << "r -JB" << std::fixed << std::setprecision(7) << centre_width_ << "/" <<
			std::fixed << std::setprecision(7) << centre_height_ << "/" <<
			std::fixed << std::setprecision(7) << upperDeg_ - fabs(dHeight_/3.) << "/" <<
			std::fixed << std::setprecision(7) << lowerDeg_ + fabs(dHeight_/3.) << "/" <<				
			std::fixed << std::setprecision(1) << page_width_ << "c";

		gmtbat_file_ << " -D" << coastResolution_ << " -N2/0.25p -W0.75p,16/169/243 -G255/255/255 -S233/246/255 -Lf" <<
			std::fixed << std::setprecision(5) << centre_width_ << "/" << std::fixed << std::setprecision(5) << dScaleLat_ << "/" << 
			std::fixed << std::setprecision(5) << centre_height_ << "/" << std::fixed << std::setprecision(scale_precision_) << scale_bar_width_ << "k+lKilometres+jt ";
		if (graticule_width_ < 0.0001)
            gmtbat_file_ << "-Bafg";
        else
            gmtbat_file_ << "-B" <<
			std::fixed << std::setprecision(graticule_width_precision_) << graticule_width_ << "g" << std::fixed << std::setprecision(graticule_width_precision_) << graticule_width_ << "/" << 
			std::fixed << std::setprecision(graticule_width_precision_) << graticule_width_ << "g" << std::fixed << std::setprecision(graticule_width_precision_) << graticule_width_;

		if (!isLandscape)
			gmtbat_file_ << " -P";
		gmtbat_file_ << " -K > " << psTempFile << std::endl;
		break;
	//
	// Lambert Azimuthal Equal-Area		
	case lambertEqualArea:		
		if (pprj_->p._label_font_size < 0.)
			pprj_->p._label_font_size = 6;

		gmtbat_file_ << _APP_PSCOAST_ << " -R" << 
		std::fixed << std::setprecision(7) << leftDeg_ << "/" <<
		std::fixed << std::setprecision(7) << lowerDeg_ << "/" <<
		std::fixed << std::setprecision(7) << rightDeg_ << "/" <<
		std::fixed << std::setprecision(7) << upperDeg_;

		// example: -JA30/-30/4.5i
		gmtbat_file_ << "r -JA" << std::fixed << std::setprecision(7) << centre_width_ << "/" <<
			std::fixed << std::setprecision(7) << centre_height_ << "/" <<
			std::fixed << std::setprecision(1) << page_width_ << "c";

		gmtbat_file_ << " -D" << coastResolution_ << " -N2/0.25p -W0.75p,16/169/243 -G255/255/255 -S233/246/255 -Lf" <<
			std::fixed << std::setprecision(5) << centre_width_ << "/" << std::fixed << std::setprecision(5) << dScaleLat_ << "/" << 
			std::fixed << std::setprecision(5) << centre_height_ << "/" << std::fixed << std::setprecision(scale_precision_) << scale_bar_width_ << "k+lKilometres+jt ";
        if (graticule_width_ < 0.0001)
            gmtbat_file_ << "-Bafg";
        else
            gmtbat_file_ << "-B" <<
			std::fixed << std::setprecision(graticule_width_precision_) << graticule_width_ << "g" << std::fixed << std::setprecision(graticule_width_precision_) << graticule_width_ << "/" << 
			std::fixed << std::setprecision(graticule_width_precision_) << graticule_width_ << "g" << std::fixed << std::setprecision(graticule_width_precision_) << graticule_width_;

		if (!isLandscape)
			gmtbat_file_ << " -P";
		gmtbat_file_ << " -K > " << psTempFile << std::endl;
		
		break;
	//
	// General stereographic map
	case stereographicConformal:
		if (pprj_->p._label_font_size < 0.)
			pprj_->p._label_font_size = 6;

		gmtbat_file_ << _APP_PSCOAST_ << " -R" << 
		std::fixed << std::setprecision(7) << leftDeg_ << "/" <<
		std::fixed << std::setprecision(7) << lowerDeg_ << "/" <<
		std::fixed << std::setprecision(7) << rightDeg_ << "/" <<
		std::fixed << std::setprecision(7) << upperDeg_;

		// example: -R100/-40/160/-10r -JS130/-30/4i
		gmtbat_file_ << "r -JS" << std::fixed << std::setprecision(7) << centre_width_ << "/" <<
			std::fixed << std::setprecision(7) << centre_height_ << "/" <<
			std::fixed << std::setprecision(1) << page_width_ << "c";

		gmtbat_file_ << " -D" << coastResolution_ << " -N2/0.25p -W0.75p,16/169/243 -G255/255/255 -S233/246/255 -Lf" <<
			std::fixed << std::setprecision(5) << centre_width_ << "/" << std::fixed << std::setprecision(5) << dScaleLat_ << "/" << 
			std::fixed << std::setprecision(5) << centre_height_ << "/" << std::fixed << std::setprecision(scale_precision_) << scale_bar_width_ << "k+lKilometres+jt ";
        if (graticule_width_ < 0.0001)
            gmtbat_file_ << "-Bafg";
        else
            gmtbat_file_ << "-B" <<
			std::fixed << std::setprecision(graticule_width_precision_) << graticule_width_ << "g" << std::fixed << std::setprecision(graticule_width_precision_) << graticule_width_ << "/" << 
			std::fixed << std::setprecision(graticule_width_precision_) << graticule_width_ << "g" << std::fixed << std::setprecision(graticule_width_precision_) << graticule_width_;

		if (!isLandscape)
			gmtbat_file_ << " -P";
		gmtbat_file_ << " -K > " << psTempFile << std::endl;
		
		break;
	// Robinson projection
	case robinson:
	default:

		// Robinson map centered on the dateline
		// pscoast -Rwest/east/south/north -JN7.5i -B60f30g30 -Dc -A5000 -Gblack -P > tmp.ps
		// pscoast -R(centre_width_-180)/(centre_width_+180)/-90/90

		// Override coastline resolution
		coastResolution_ = "i";
		pprj_->p._coasline_resolution = intermediate;
		circle_radius_ = 0.02;
		circle_line_width_ = 0.02;
		if (pprj_->p._label_font_size < 0.)
			pprj_->p._label_font_size = 5;

		gmtbat_file_ << _APP_PSCOAST_ << " -R" << 
			(centre_width_-180) << "/" << (centre_width_+180) << "/-90/90"
			" -JN" << page_width_  << "c -B60g30 -D" << coastResolution_ << 
			" -A10000 -W0.75p,16/169/243 -G245/245/245 -K > " << psTempFile << std::endl;
	}

	std::sort(pprj_->p._msr_colours.begin(), pprj_->p._msr_colours.end(), ComparePairFirst<std::string>());		

	if (pprj_->p._plot_plate_boundaries)
	{
		// print plate boundaries first
		gmtbat_file_ << _APP_PSXY_ << " -R -J \"" << v_tectonic_plate_file_.at(block) << 
			"\" -W0.75p,#DA5552 -O -K >> " << psTempFile << std::endl;
	}

	// Does the user want to print measurements?
	if (!pprj_->p._omit_measurements)
	{
		// print latitude, longitude and height measurements (large circles) before stations
		for (UINT32 i=0; i<v_msr_file_.at(block).size(); i++)
		{
			if (v_msr_file_.at(block).at(i).second != "H" &&
				v_msr_file_.at(block).at(i).second != "I" &&
				v_msr_file_.at(block).at(i).second != "J" &&
				v_msr_file_.at(block).at(i).second != "P" &&
				v_msr_file_.at(block).at(i).second != "Q" &&
				v_msr_file_.at(block).at(i).second != "R" &&
				v_msr_file_.at(block).at(i).second != "Y")
				continue;

			circle_radius_2_ = circle_radius_;

			if (v_msr_file_.at(block).at(i).second == "I" ||		// astronomic latitude
				v_msr_file_.at(block).at(i).second == "P")			//   geodetic latitude
				circle_radius_2_ = circle_radius_ * 2.75;
			else if (v_msr_file_.at(block).at(i).second == "J" ||	// astronomic longitude
				v_msr_file_.at(block).at(i).second == "Q")			//   geodetic longitude
				circle_radius_2_ = circle_radius_ * 2.4;
			else if (v_msr_file_.at(block).at(i).second == "H" ||	// orthometric height
				v_msr_file_.at(block).at(i).second == "R")			// ellipsoidal height
				circle_radius_2_ = circle_radius_ * 1.95;
			else if (v_msr_file_.at(block).at(i).second == "Y")		// GPS point cluster
				circle_radius_2_ = circle_radius_ * 3.0;

			_it_colour = equal_range(pprj_->p._msr_colours.begin(), pprj_->p._msr_colours.end(), 
				v_msr_file_.at(block).at(i).second, ComparePairFirst<std::string>());

			if (_it_colour.first == _it_colour.second)
			{
				gmtbat_file_ << _APP_PSXY_ << " -R -J -Skcircle/" << std::fixed << std::setprecision(2) << circle_radius_2_ << 
					" \"" << v_msr_file_.at(block).at(i).first << "\" -W" << std::setprecision(2) << circle_line_width_ << 
					"p,darkgray -Glightgray";
				gmtbat_file_ << " -O -K >> " << psTempFile << std::endl;
			}
			else
			{
				colours++;
				gmtbat_file_ << _APP_PSXY_ << " -R -J -Skcircle/" << std::fixed << std::setprecision(2) << circle_radius_2_ << 
					" \"" << v_msr_file_.at(block).at(i).first << "\" -W" << std::setprecision(2) << circle_line_width_ << 
					"p,";

				// Line colour
				if (v_msr_file_.at(block).at(i).second == "Y")						// GNSS Point Cluster
					gmtbat_file_ << "#235789";									// bdazzled blue
				else if (v_msr_file_.at(block).at(i).second == "I" ||				// astronomic latitude
					v_msr_file_.at(block).at(i).second == "P")						// geodetic latitude
					gmtbat_file_ << "#A393BF";									// glossy grape
				else if (v_msr_file_.at(block).at(i).second == "J" ||				// astronomic longitude
					v_msr_file_.at(block).at(i).second == "Q")						// geodetic longitude
					gmtbat_file_ << "#4C1C00";									// seal brown
				else if (v_msr_file_.at(block).at(i).second == "H" ||				// orthometric height
					v_msr_file_.at(block).at(i).second == "R")						// ellipsoidal height
					gmtbat_file_ << "#6622CC";									// french violet
				else
					gmtbat_file_ << _it_colour.first->second;

				// Fill colour
				gmtbat_file_ << " -G" << _it_colour.first->second  << " -O -K >> " << psTempFile << std::endl;
			}
		}
	}

	// print stations first to enable measurements to be seen
	gmtbat_file_ << _APP_PSXY_ << " -R -J -Skcircle/" << std::fixed << std::setprecision(2) << circle_radius_ << 
		" \"" << v_isl_pts_file_.at(block) << "\" -W" << std::setprecision(2) << circle_line_width_ << 
		"p,#235789 -Gwhite -O -K >> " << psTempFile << std::endl;

	if (plotConstraints_)
		// don't plot line, just fill
		gmtbat_file_ << _APP_PSXY_ << " -R -J -Skcircle/" << std::fixed << std::setprecision(2) << circle_radius_ << 
			" \"" << v_isl_const_file_.at(block) << "\" -G#235789 -O -K >> " << psTempFile << std::endl;

	if (plotBlocks_)
	{
		gmtbat_file_ << _APP_PSXY_ << " -R -J -Skcircle/" << std::fixed << std::setprecision(2) << circle_radius_ << 
			" \"" << v_jsl_pts_file_.at(block) << "\" -W" << std::setprecision(2) << circle_line_width_ * 2.0 << 
			"p,#DA5552 -Gwhite -O -K >> " << psTempFile << std::endl;

		if (plotConstraints_)
			// don't plot line, just fill
			gmtbat_file_ << _APP_PSXY_ << " -R -J -Skcircle/" << std::fixed << std::setprecision(2) << circle_radius_ << 
				" \"" << v_jsl_const_file_.at(block) << "\" -G#DA5552 -O -K >> " << psTempFile << std::endl;
	}

	// Does the user want to print measurements?
	if (!pprj_->p._omit_measurements)
	{
		// now print lines
		for (UINT32 i=0; i<v_msr_file_.at(block).size(); i++)
		{
			if (v_msr_file_.at(block).at(i).second == "H")		// already printed
				continue;
			if (v_msr_file_.at(block).at(i).second == "I")		// ''
				continue;
			if (v_msr_file_.at(block).at(i).second == "J")		// ''
				continue;
			if (v_msr_file_.at(block).at(i).second == "P")		// ''
				continue;
			if (v_msr_file_.at(block).at(i).second == "Q")		// ''
				continue;			
			if (v_msr_file_.at(block).at(i).second == "R")		// ''
				continue;
			if (v_msr_file_.at(block).at(i).second == "Y")		// ''
				continue;
				
			_it_colour = equal_range(pprj_->p._msr_colours.begin(), pprj_->p._msr_colours.end(), 
				v_msr_file_.at(block).at(i).second, ComparePairFirst<std::string>());

			if (_it_colour.first == _it_colour.second)
			{
				gmtbat_file_ << _APP_PSXY_ << " -R -J \"" << v_msr_file_.at(block).at(i).first << 
					"\" " << "-W" << std::setprecision(2) << line_width_ << 
					"p,lightgray";
				//if (v_msr_file_.at(block).at(i).second == "Y")		// not a vector measurement, so represent as dashed
				//	gmtbat_file_ << ",6_8:1p";
				gmtbat_file_ << " -O -K >> " << psTempFile << std::endl;
			}
			else
			{
				colours++;
				gmtbat_file_ << _APP_PSXY_ << " -R -J \"" << v_msr_file_.at(block).at(i).first << 
					"\" " << "-W" << std::setprecision(2) << line_width_ << 
					"p," << _it_colour.first->second;
				//if (v_msr_file_.at(block).at(i).second == "Y")		// not a vector measurement, so represent as dashed
				//	gmtbat_file_ << ",6_8:1p";
				gmtbat_file_ << " -O -K >> " << psTempFile << std::endl;
			}
		}
	}

	if (pprj_->p._plot_positional_uncertainty || 
		pprj_->p._plot_error_ellipses)
	{
		// calculate scale for uncertainty circles, to the end that 
		// largest_uncertainty_ is uncertainty_legend_length_ cm on the A3 page
		//if (largest_uncertainty_ >= 1.0)
		//	uncertainty_legend_precision = 0;
		/*else*/ if ((largest_uncertainty_ / pprj_->p._pu_ellipse_scale) >= 0.1)
			uncertainty_legend_precision = 1;
		else if ((largest_uncertainty_ / pprj_->p._pu_ellipse_scale) >= 0.01)
			uncertainty_legend_precision = 2;
		else if ((largest_uncertainty_ / pprj_->p._pu_ellipse_scale) >= 0.001)
			uncertainty_legend_precision = 3;
		else
			uncertainty_legend_precision = 4;
	}

	UINT32 precision(10);

	if (pprj_->p._plot_positional_uncertainty)
	{
		// print positional uncertainty
		gmtbat_file_ << _APP_PSVELO_ << " -R -J \"" << v_stn_apu_file_.at(block) << 
			"\" -Sr" << error_ellipse_scale << "/0.95/0c -L -W1.25p,#FFD275 -O -K >> " << psTempFile << std::endl;

		// Shift plot north to account for error ellipse legend
		if (pprj_->p._plot_error_ellipses)
			gmtbat_file_ << _APP_PSXY_ << " -R -J -T -Y" << 
				label_font_size * 1.5 << "p " << 
				" -O -K >> " << psTempFile << std::endl;
		
		// Add text for positional uncertainty legend 
		gmtbat_file_ << _ECHO_CMD_ << 
			std::setprecision(precision) << std::fixed << std::left << uncertainty_legend_long_ << " " << 
			std::setprecision(precision) << std::fixed << std::left << uncertainty_legend_lat_ << " " << 
			" 95% positional uncertainty " <<						// the label
			std::setprecision(uncertainty_legend_precision) <<			// ''
			std::fixed << largest_uncertainty_ / pprj_->p._pu_ellipse_scale << " radius \\(m\\)" <<					// radius
			" | ";													// Push to pstext
		gmtbat_file_ << 
			_APP_PSTEXT_ << " -R -J -F+f" << 
			std::fixed << std::setprecision(0) << 
			pprj_->p._label_font_size * 2.0 << "p,Helvetica" <<		// font size and face
			"=~" << pprj_->p._label_font_size/2.0 << "p,white" << 		// outline (or glow)
			"+jLM " <<														// justification
			"-Dj" << pprj_->p._label_font_size * 2.0 <<				// x shift
			"p/0p -O -K >> " << psTempFile << std::endl;					// y shift
	}
	
	if (pprj_->p._plot_error_ellipses)
	{
		// Terminate a sequence of GMT plotting commands without producing any plotting output
		if (pprj_->p._plot_positional_uncertainty)
			gmtbat_file_ << _APP_PSXY_ << " -R -J -T -Y-" << 
				label_font_size * 1.5 << "p " << 
				" -O -K >> " << psTempFile << std::endl;

		double ellipse_scale(1.);

		// The error ellipse legend is printed using the largest uncertainty (semi-major) found, unless
		// positional uncertainty is being printed, in which case the error ellipse is printed using
		// half the size.  This is purely to give the impression that error ellipses are 1 sigma (68%) and
		// positional uncertainty is 95%, or close to 2 sigma.
		if (pprj_->p._plot_positional_uncertainty)
			ellipse_scale = 2.;

		// print error ellipses in indian red
		gmtbat_file_ << _APP_PSVELO_ << " -R -J \"" << v_stn_err_file_.at(block) << 
			"\" -Sr" << error_ellipse_scale << "/0.95/0c -L -W0.75p,#DA5552 -O -K >> " << psTempFile << std::endl;

		// Add text for error ellipse legend 
		gmtbat_file_ << _ECHO_CMD_ <<
			std::setprecision(precision) << std::fixed << std::left << uncertainty_legend_long_ << " " << 
			std::setprecision(precision) << std::fixed << std::left << uncertainty_legend_lat_ << " " << 
			" 1 sigma error ellipse " <<							// the label
			std::setprecision(uncertainty_legend_precision) <<			// ''
			std::fixed << largest_uncertainty_ / ellipse_scale / pprj_->p._pu_ellipse_scale <<			// semi-major
			", " <<				
			std::fixed << largest_uncertainty_ / 3.0 / ellipse_scale / pprj_->p._pu_ellipse_scale <<	// semi-minor (make it third the height)
			" \\(m\\)" <<
			" | ";													// Push to pstext
		gmtbat_file_ << 
			_APP_PSTEXT_ << " -R -J -F+f" << 
			std::fixed << std::setprecision(0) << 
			pprj_->p._label_font_size * 2.0 << "p,Helvetica" <<				// font size and face
			"=~" << pprj_->p._label_font_size/2.0 << "p,white" << 		// outline (or glow)
			"+jLM " <<														// justification
			"-Dj" << pprj_->p._label_font_size * 2.0 <<				// x shift
			"p/0p -O -K >> " << psTempFile << std::endl;					// y shift
	}
	
	if (pprj_->p._plot_correction_arrows)
	{
		// calculate scale for arrows, to the end that 
		// average_correction_ is 3cm on the A3 page
		if (average_correction_ < PRECISION_1E4)
			average_correction_ = PRECISION_1E4;
		double correction_arrow_scale(arrow_legend_length_/average_correction_);
		UINT32 correction_legend_precision(4);

		//if (average_correction_ >= 1.0)
		//	correction_legend_precision = 0;
		/*else*/ if (average_correction_ >= 0.1)
			correction_legend_precision = 1;
		else if (average_correction_ >= 0.01)
			correction_legend_precision = 2;
		else if (average_correction_ >= 0.001)
			correction_legend_precision = 3;
		else
			correction_legend_precision = 4;


		// print text in black (no arrows)
		if (pprj_->p._plot_correction_labels)
			gmtbat_file_ << _APP_PSVELO_ << " -R -J \"" << v_stn_cor_file_.at(block) << 
				"\" -Se0.0001/0.95+f" << std::setprecision(0) << pprj_->p._label_font_size << 
				"p,Helvetica,black -L -A0.0001/0.0001/0.0001c -O -K >> " << psTempFile << std::endl;

		// print arrows (without black outline!)
		gmtbat_file_ << _APP_PSVELO_ << " -R -J \"" << v_stn_cor_file_.at(block) << 
			"\" -Se" << correction_arrow_scale << "/0.95/0c -L -A" << line_width_ / 2.0 * CM_TO_INCH << 
			"/0.5/0.1c -G#DA5552 -W0.0p,#DA5552 -O -K >> " << psTempFile << std::endl;

		// Add text below the corrections arrow legend 
		gmtbat_file_ << _ECHO_CMD_ <<
			std::setprecision(precision) << std::fixed << std::left << arrow_legend_long_ << " " << 
			std::setprecision(precision) << std::fixed << std::left << arrow_legend_lat_ << " " << 
			std::setprecision(correction_legend_precision) <<			// the label
			std::fixed << average_correction_ / pprj_->p._correction_scale <<	" \\(m\\)" <<				// ''
			" | ";													// Push to pstext
		gmtbat_file_ << 
			_APP_PSTEXT_ << " -R -J -F+f" << 
			std::fixed << std::setprecision(0) << 
			pprj_->p._label_font_size * 2.0 << "p,Helvetica" <<		// font size and face
			"=~" << pprj_->p._label_font_size/2.0 << "p,white" << 		// outline (or glow)
			"+jRM " <<														// justification
			"-Dj" << pprj_->p._label_font_size * 2.0 << 		 		// x shift
			"p/0p -O -K >> " << psTempFile << std::endl;					// y shift
		
		gmtbat_file_ << _APP_PSXY_ << " -R -J -T -Y-" << 
			label_font_size * 1.5 << "p " << 
			" -O -K >> " << psTempFile << std::endl;

		// Add text above for corrections legend 
		gmtbat_file_ << _ECHO_CMD_ <<
			std::setprecision(precision) << std::fixed << std::left << arrow_legend_long_ << " " << 
			std::setprecision(precision) << std::fixed << std::left << arrow_legend_lat_ << " " << 
			"Corrections scale" <<									// the label
			" | ";													// Push to pstext
		gmtbat_file_ << 
			_APP_PSTEXT_ << " -R -J -F+f" << 
			std::fixed << std::setprecision(0) << 
			pprj_->p._label_font_size * 2.0 << "p,Helvetica" <<		// font size and face
			"=~" << pprj_->p._label_font_size/2.0 << "p,white" << 		// outline (or glow)
			"+jCM " <<														// justification
			"-Dj-" << pprj_->p._label_font_size * 2.0 <<				// x shift
			"p/0p -O -K >> " << psTempFile << std::endl;					// y shift
		
		gmtbat_file_ << _APP_PSXY_ << " -R -J -T -Y" << 
			label_font_size * 1.5 << "p " << 
			" -O -K >> " << psTempFile << std::endl;
	}

	if (pprj_->p._plot_station_labels)
	{
		gmtbat_file_ << _APP_PSTEXT_ << " -R -J \"" << v_isl_lbl_file_.at(block) << "\"";
		
		// print shadow
		gmtbat_file_ << " -F+f" << 
			std::fixed << std::setprecision(0) << 
			pprj_->p._label_font_size << "p,Helvetica" <<				// font size and face
			"=~" << pprj_->p._label_font_size/3.0 << "p,white" << 		// outline (or glow)
			"+jLM " <<														// justification
			"-Dj" << pprj_->p._label_font_size + 						
				sqrt(pprj_->p._label_font_size) << 					// x shift
			"p/" << pprj_->p._label_font_size/2.0 <<					// y shift
			"p -O -K >> " << psTempFile << std::endl;

		if (plotBlocks_)
		{
			gmtbat_file_ << _APP_PSTEXT_ << " -R -J \"" << v_jsl_lbl_file_.at(block) << "\"";
			
			// print shadow
			gmtbat_file_ << " -F+f" << 
				std::fixed << std::setprecision(0) << 
				pprj_->p._label_font_size << "p,Helvetica" <<				// font size and face
				"=~" << pprj_->p._label_font_size/3.0 << "p,white" << 		// outline (or glow)
				"+jLM " <<														// justification
				"-Dj" << pprj_->p._label_font_size + 						
					sqrt(pprj_->p._label_font_size) << 					// x shift
				"p/" << pprj_->p._label_font_size/2.0 <<					// y shift
				"p -O -K >> " << psTempFile << std::endl;
		}
	}

	// Prepare to terminate a sequence of GMT plotting commands without producing 
	// any plotting output (depending on whether a title block is required)
	gmtbat_file_ << _APP_PSXY_ << " -R -J -T -O";

	std::string legendCommand1(_LEGEND_CMD_1_), legendCommand2(_LEGEND_CMD_2_);

#if defined(_WIN32) || defined(__WIN32__)
	// temporary legend files
	legendCommand1.append(legendTempFile);
	legendCommand2.append(legendTempFile);
#endif
	
	// No title block?
	if (pprj_->p._omit_title_block)
	{
		// Terminate a sequence of GMT plotting commands without producing any plotting output
		gmtbat_file_ << " >> " << psTempFile << std::endl << std::endl;
	}
	else
	{
		gmtbat_file_ << " -K >> " << psTempFile << std::endl << std::endl;
		
		// legend
		gmtbat_file_ << _APP_GMTSET_ << " FONT_TITLE 1" << std::endl;
		gmtbat_file_ << _APP_GMTSET_ << " FONT_ANNOT_PRIMARY " << label_font_size * 1.1 << "p";

		gmtbat_file_ << std::endl;

#if defined(__linux) || defined(sun) || defined(__unix__) || defined(__APPLE__)
		gmtbat_file_ << "cat > " << legendTempFile << " << END" << std::endl;
#endif

		bool isnameaNumber(is_number<std::string>(pprj_->p._title));

		gmtbat_file_ << _LEGEND_ECHO_ << "G 0.25" << legendCommand1 << std::endl;
		gmtbat_file_ << _LEGEND_ECHO_ << "N 1" << legendCommand2 << std::endl;
		gmtbat_file_ << _LEGEND_ECHO_ << "H " << std::fixed << std::setprecision(0) << label_font_size * 2.0 << "p,Helvetica-Bold ";
		if (isnameaNumber) 
			gmtbat_file_ << "'" << pprj_->p._title << "'";
		else
			gmtbat_file_ << pprj_->p._title;

		if (blockCount_ > 1)
			gmtbat_file_ << "  (Block " << block + 1 << " of " << blockCount_ << ") ";

		if (!pprj_->p._plot_station_centre.empty())
			gmtbat_file_ << " centred on " << pprj_->p._plot_station_centre;

		gmtbat_file_ << legendCommand2 << std::endl;
		gmtbat_file_ << _LEGEND_ECHO_ << "G 0.25" << legendCommand2 << std::endl;

		// Print stations legend
		gmtbat_file_ << _LEGEND_ECHO_ << "D 0 1p" << legendCommand2 << std::endl;		// horizontal line
		gmtbat_file_ << _LEGEND_ECHO_ << "G 0.25" << legendCommand2 << std::endl;		// space	

		UINT32 station_count(1);		// Simultaneous or ISL

		if (plotBlocks_)
			station_count++;			// JSL
		
		if (plotConstraints_)
		{
			station_count++;			// Simultaneous or ISL constraint stations
			if (plotBlocks_)
				station_count++;		// JSL constraint stations				
		}
		
		//if (pprj_->p._plot_correction_arrows)
		//	station_count++;			// Corrections to station coordinates

		gmtbat_file_ << _LEGEND_ECHO_ << "N " << station_count << legendCommand2 << std::endl;
		gmtbat_file_ << _LEGEND_ECHO_ << "V 0 1p" << legendCommand2 << std::endl;		// vertical line
		
		circle_radius_2_ = circle_radius_ * 1.75;

		// Simultaneous stations or Phased inner stations
		gmtbat_file_ << _LEGEND_ECHO_ <<
			"S " << std::fixed << std::setprecision(1) << symbol_offset << " c " << 
			circle_radius_2_ << " white " << 
			circle_line_width_ * 2 << "p,#235789 " <<
			std::fixed << std::setprecision(1) << label_offset;
		if (plotBlocks_)
			gmtbat_file_ << " Free Inner stations ";
		else 
			gmtbat_file_ << " Free Stations";
		gmtbat_file_ << legendCommand2 << std::endl;

		// Simultaneous stations or Phased inner stations
		if (plotBlocks_)
		{
			gmtbat_file_ << _LEGEND_ECHO_ <<
				"S " << std::fixed << std::setprecision(1) << symbol_offset << " c " << 
				circle_radius_2_ << " white " << 
				circle_line_width_ * 2 << "p,#DA5552 " <<
				std::fixed << std::setprecision(1) << label_offset <<
				" Free Junction stations" << legendCommand2 << std::endl;
		}

		if (plotConstraints_)
		{
			if (plotBlocks_)
			{
				gmtbat_file_ << _LEGEND_ECHO_ <<
					"S " << std::fixed << std::setprecision(1) << symbol_offset << " c " << 
					circle_radius_2_ << " #235789 " << 
					circle_line_width_ * 2 << "p,#235789 " <<
					std::fixed << std::setprecision(1) << label_offset <<
					" Constrained inner stations" << legendCommand2 << std::endl;
				gmtbat_file_ << _LEGEND_ECHO_ <<
					"S " << std::fixed << std::setprecision(1) << symbol_offset << " c " << 
					circle_radius_2_ << " #DA5552 " << 
					circle_line_width_ * 2 << "p,#DA5552 " <<
					std::fixed << std::setprecision(1) << label_offset <<
					" Constrained junction stations" << legendCommand2 << std::endl;
			}
			else
				gmtbat_file_ << _LEGEND_ECHO_ <<
					"S " << std::fixed << std::setprecision(1) << symbol_offset << " c " << 
					circle_radius_2_ << " #235789 " << 
					circle_line_width_ * 2 << "p,#235789 " <<
					std::fixed << std::setprecision(1) << label_offset <<
					" Constraint stations" << legendCommand2 << std::endl;
		}			
		
		title_block_height_ = 5;

		if (pprj_->p._plot_correction_arrows)
		{
		//	// Add corrections to end of stations legend
		//	gmtbat_file_ << _LEGEND_ECHO_ <<
		//		"S " << std::fixed << std::setprecision(1) << symbol_offset << 
		//		" v 1.5c/0.04/0.5/0.1 red " <<			// arrowlength/linewidth/arrowheadwidth/arrowheadlength
		//		line_width * 2 << "p,#DA5552 " 
		//		<< std::fixed << std::setprecision(1) << label_offset * 1.5 << " Corrections to stations" << legendCommand2 << std::endl;
		}
		// If corrections are not being printed, then print measurements legend
		else
		{
			gmtbat_file_ << _LEGEND_ECHO_ << "G 0.25" << legendCommand2 << std::endl;		// space
			gmtbat_file_ << _LEGEND_ECHO_ << "D 0 1p" << legendCommand2 << std::endl;		// horizontal line
			gmtbat_file_ << _LEGEND_ECHO_ << "G 0.25" << legendCommand2 << std::endl;		// space	

			if (v_msr_file_.at(block).size() > columns)
				gmtbat_file_ << _LEGEND_ECHO_ << "N " << columns << legendCommand2 << std::endl;
			else
				gmtbat_file_ << _LEGEND_ECHO_ << "N " << v_msr_file_.at(block).size() << legendCommand2 << std::endl;

			bool bOtherTypes(false);
			
			for (UINT32 i=0; i<v_msr_file_.at(block).size(); i++)
			{
				_it_colour = equal_range(pprj_->p._msr_colours.begin(), pprj_->p._msr_colours.end(), 
					v_msr_file_.at(block).at(i).second, ComparePairFirst<std::string>());

				if (_it_colour.first != _it_colour.second)
				{
					gmtbat_file_ << _LEGEND_ECHO_ << "S " << std::fixed << std::setprecision(1) << symbol_offset;

					circle_radius_2_ = circle_radius_ * 1.75;

					// print a circle or line?
					if (v_msr_file_.at(block).at(i).second == "I" ||			// astronomic latitude
						v_msr_file_.at(block).at(i).second == "P" ||			//   geodetic latitude
						v_msr_file_.at(block).at(i).second == "J" ||			// astronomic longitude
						v_msr_file_.at(block).at(i).second == "Q" ||			//   geodetic longitude
						v_msr_file_.at(block).at(i).second == "H" ||			// orthometric height
						v_msr_file_.at(block).at(i).second == "R" ||			// ellipsoidal height
						v_msr_file_.at(block).at(i).second == "Y")				// GPS point cluster
						gmtbat_file_ << " c " << circle_radius_2_ << " " <<
								_it_colour.first->second << " " << 
								circle_line_width_ * 2;								// circle
					else
						gmtbat_file_ << " - " << line_width_ * 3 << " " << 
						_it_colour.first->second <<	" " << line_width_ * 7.5;	// line

					gmtbat_file_ << "p,";

					if (v_msr_file_.at(block).at(i).second == "Y")
						gmtbat_file_ << "#235789";							// bdazzled blue
					else if (v_msr_file_.at(block).at(i).second == "I" ||		// astronomic latitude
						v_msr_file_.at(block).at(i).second == "P")				// geodetic latitude
						gmtbat_file_ << "#A393BF";							// glossy grape
					else if (v_msr_file_.at(block).at(i).second == "J" ||		// astronomic longitude
						v_msr_file_.at(block).at(i).second == "Q")				// geodetic longitude
						gmtbat_file_ << "#4C1C00";							// seal brown
					else if (v_msr_file_.at(block).at(i).second == "H" ||		// orthometric height
						v_msr_file_.at(block).at(i).second == "R")				// ellipsoidal height
						gmtbat_file_ << "#6622CC";							// french violet
					else
						gmtbat_file_ << _it_colour.first->second;

					gmtbat_file_ << " " << std::fixed << std::setprecision(1) << label_offset << " " << 
						measurement_name<char, std::string>(static_cast<char>(_it_colour.first->first.at(0))) << legendCommand2 << std::endl;
				}
				else
					bOtherTypes = true;
			}

			if (bOtherTypes)
				gmtbat_file_ << _LEGEND_ECHO_ <<
					"S " << std::fixed << std::setprecision(1) << symbol_offset << " - 0.5 lightgray " << line_width_ * 2.5 << "p,lightgray " 
					<< std::fixed << std::setprecision(1) << label_offset << "  All other types" << legendCommand2 << std::endl;			

			title_block_height_ += floor(((double)colours)/columns) * 0.25;
		
		}

		gmtbat_file_ << _LEGEND_ECHO_ << "G 0.25" << legendCommand2 << std::endl;		// space
		gmtbat_file_ << _LEGEND_ECHO_ << "D 0 1p" << legendCommand2 << std::endl;		// horizontal line
		gmtbat_file_ << _LEGEND_ECHO_ << "G 0.25" << legendCommand2 << std::endl;		// space
		gmtbat_file_ << _LEGEND_ECHO_ << "N 5" << legendCommand2 << std::endl;
		gmtbat_file_ << _LEGEND_ECHO_ << "S 0.01 c 0.01 white 1p,white 1 " << pprj_->p._title_block_subname << legendCommand2 << std::endl;
		gmtbat_file_ << _LEGEND_ECHO_ << "S 0.01 c 0.01 white 1p,white 0 " << pprj_->p._title_block_name << legendCommand2 << std::endl;
		gmtbat_file_ << _LEGEND_ECHO_ << "S 0.01 c 0.01 white 1p,white 3 " << reference_frame_ << legendCommand2 << std::endl;
		gmtbat_file_ << _LEGEND_ECHO_ << "S 0.01 c 0.01 white 1p,white 0 " << projectionTypes[pprj_->p._projection] << " projection" << legendCommand2 << std::endl;
		gmtbat_file_ << _LEGEND_ECHO_ << "S 0.01 c 0.01 white 1p,white 1 Scale 1:" << static_cast<UINT32>(scale_) << " (A3)" << legendCommand2 << std::endl;

#if defined(__linux) || defined(sun) || defined(__unix__) || defined(__APPLE__)
		gmtbat_file_ << "END" << std::endl << std::endl;
#endif

		gmtbat_file_ << _APP_PSLEGEND_ << " -R -J -DJTL+w" << std::fixed << std::setprecision(1) << page_width_ << 
			"c+jBL+l1.5+o0/1.5c -C0.3/0.3 -O -F+p+gwhite " <<
			legendTempFile << " -P >> " << psTempFile << std::endl;

	}


	////////////////////////////////////////////////
	// conversion to PDF 
	//
	
	gmtbat_file_ << std::endl << _COMMENT_PREFIX_ << "convert ps to pdf" << std::endl;

	size_t lastindex = v_gmt_pdf_filenames_.at(block).find_last_of("."); 
	std::string pdf_filename = v_gmt_pdf_filenames_.at(block).substr(0, lastindex); 	
	gmtbat_file_ << _APP_PSCONVERT_ << " -A0.2c+white -Tf -F\"" << pdf_filename << "\" " << psTempFile << std::endl << std::endl;

	if (pprj_->p._export_png)
	{
		gmtbat_file_ << std::endl << _COMMENT_PREFIX_ << "convert ps to png" << std::endl;
		gmtbat_file_ << _APP_PSCONVERT_ << " -A0.2c+white -Tg -F\"" << pdf_filename << "\" " << psTempFile << std::endl << std::endl;
	}

	////////////////////////////////////////////////
	
	////////////////////////////////////////////////
	// clean up legend files
	gmtbat_file_ << _COMMENT_PREFIX_ << "cleanup" << std::endl;
	gmtbat_file_ << _DELETE_CMD_ << psTempFile;

	if (!pprj_->p._omit_title_block)
		gmtbat_file_ << " " << legendTempFile << std::endl;
	else
		gmtbat_file_ << std::endl;

	////////////////////////////////////////////////

	// remove temporary directory for gmt.conf and gmt.history
	gmtbat_file_ << _RMDIR_CMD_ << _ENV_GMT_TMP_DIR_ << std::endl;

#if defined(__linux) || defined(sun) || defined(__unix__) || defined(__APPLE__)	
	gmtbat_file_ << "unset " << _GMT_TMP_DIR_ << std::endl;

#endif		

		
	gmtbat_file_ << std::endl;
		
	//////////////////////////////////////////////////////////////////////////////////////////
	// delete data files
	if (!pprj_->p._keep_gen_files)
	{
		gmtbat_file_ << _DELETE_CMD_;

		// stations
		gmtbat_file_ << 
			" \"" << v_isl_pts_file_.at(block) << "\"" <<
			" \"" << v_isl_const_file_.at(block) << "\"";

		if (plotBlocks_)
		{
			// junction stations5
			gmtbat_file_ << 
				" \"" << v_jsl_const_file_.at(block) << "\"" <<
				" \"" << v_jsl_pts_file_.at(block) << "\"";
		}

		// measurements
		for (UINT32 i=0; i<v_msr_file_.at(block).size(); i++)
			gmtbat_file_ << " \"" << v_msr_file_.at(block).at(i).first << "\"";
	
		// station labels
		if (pprj_->p._plot_station_labels)
		{
			gmtbat_file_ << " \"" << v_isl_lbl_file_.at(block) << "\"";
		
			if (plotBlocks_)
				gmtbat_file_ << " \"" << v_jsl_lbl_file_.at(block) << "\"";
		}

		if (pprj_->p._plot_positional_uncertainty)
			gmtbat_file_ << " \"" << v_stn_apu_file_.at(block) << "\"";

		if (pprj_->p._plot_error_ellipses)
			gmtbat_file_ << " \"" << v_stn_err_file_.at(block) << "\"";

		if (pprj_->p._plot_correction_arrows)
			gmtbat_file_ << " \"" << v_stn_cor_file_.at(block) << "\"";

		if (pprj_->p._plot_plate_boundaries)
			gmtbat_file_ << " \"" << v_tectonic_plate_file_.at(block) << "\"";

		gmtbat_file_ << std::endl;
	}
	else
		gmtbat_file_ << std::endl;
	//////////////////////////////////////////////////////////////////////////////////////////

	gmtbat_file_ << std::endl;

}

void dna_plot::CreateGMTPlotEnvironment(project_settings* pprj)
{
	// Set up the environment
	pprj_ = pprj;

	InitialiseGMTParameters();

	InitialiseGMTFilenames();
	CreateGMTInputFiles();

	FinaliseGMTParameters();

	// for 1..n blocks, create a shell script that generates a PDF for 
	// a single block.  dnaplotwrapper can then execute each shell
	// script in parallel.
	CreateGMTCommandFiles();
	
	// Now, based on font size determined by CreateGMTCommandFiles,
	// Print:
	//  - stations labels, 
	//  - positional uncertainty, 
	//  - error ellipses, and 
	//  - correction arrows
	CreateExtraInputFiles();
}

void dna_plot::InvokeGMT()
{
	// set up threads to execute the GMT scripts in parallel
	std::vector<std::thread> gmt_plot_threads;
	
	for (UINT32 plot=0; plot<v_gmt_cmd_filenames_.size(); ++plot)
	{
		gmt_plot_threads.emplace_back(dna_create_threaded_process(v_gmt_cmd_filenames_.at(plot)));
	}

	// go!
	for (auto& t : gmt_plot_threads) {
		t.join();
	}
}

// Aggregate individual PDFs created for each block (for
// phased adjustments only).
void dna_plot::AggregateGMTPDFs()
{
	if (!plotBlocks_)
		return;

	std::string system_file_cmd;
	std::stringstream ss;

	ss << _PDF_AGGREGATE_;

	for_each(v_gmt_pdf_filenames_.begin(), v_gmt_pdf_filenames_.end(),
		[&ss](const std::string& gmt_pdf_file) {
			ss << " " << gmt_pdf_file;
		});

#if defined(_WIN32) || defined(__WIN32__)
	// assumes pdftk is used.  add extra options
	ss << " cat output " << pprj_->p._pdf_file_name << " dont_ask";
#elif defined(__linux) || defined(sun) || defined(__unix__) || defined(__APPLE__)
	// assumes pdfunite is used (nothing to do here)
	ss << " " << pprj_->p._pdf_file_name;
#endif

	// aggregate!
	system_file_cmd = ss.str();

	std::system(system_file_cmd.c_str());
}


void dna_plot::CreateGMTPlot()
{	
	// Execute GMT in parallel
	InvokeGMT();

	// Combine all multi-page PDFs in one
	AggregateGMTPDFs();

	// clean up config and PDF files
	CleanupGMTFiles();
}

void dna_plot::CleanupGMTFiles()
{
	if (pprj_->p._keep_gen_files)
		return;
	
	std::stringstream ss;
	ss << _DELETE_CMD_;

	// shell scripts
	for_each(v_gmt_cmd_filenames_.begin(), v_gmt_cmd_filenames_.end(),
		[&ss](const std::string& gmt_cmd_file) {
			ss << " \"" << gmt_cmd_file << "\"";
		});

	if (plotBlocks_)
	{
		// PDF files for each block (except simultaneous)
		for_each(v_gmt_pdf_filenames_.begin(), v_gmt_pdf_filenames_.end(),
			[&ss](const std::string& gmt_pdf_file) {
				ss << " \"" << gmt_pdf_file << "\"";
			});
	}

	// delete
	std::string clean_up_gmt_config_files = ss.str();
	std::system(clean_up_gmt_config_files.c_str());

}
	

void dna_plot::DetermineBoundingBox()
{
	if (pprj_->p._bounding_box.empty())
		return;
	
	if (GetFields(const_cast<char*>(pprj_->p._bounding_box.c_str()), ',', false, "ffff", 
		&upperDeg_, &leftDeg_, &lowerDeg_, &rightDeg_ ) < 4)
		return;
	
	upperDeg_ = DmstoDeg(upperDeg_);
	rightDeg_ = DmstoDeg(rightDeg_);
	lowerDeg_ = DmstoDeg(lowerDeg_);
	leftDeg_ = DmstoDeg(leftDeg_);

	pprj_->p._plot_centre_longitude = (rightDeg_ + leftDeg_) / 2.0;
	pprj_->p._plot_centre_latitude = (upperDeg_ + lowerDeg_) / 2.0;

	default_limits_ = false;
}
	

void dna_plot::CalculateLimitsFromStation()
{
	it_pair_string_vUINT32 _it_stnmap;
	_it_stnmap = equal_range(stnsMap_.begin(), stnsMap_.end(), pprj_->p._plot_station_centre, StationNameIDCompareName());

	if (_it_stnmap.first == _it_stnmap.second)
	{
		std::stringstream ss;
		ss.str("");
		ss << pprj_->p._plot_station_centre << " is not in the list of network stations.";
		SignalExceptionPlot(ss.str(), 0, NULL);
	}

	pprj_->p._plot_centre_latitude = bstBinaryRecords_.at(_it_stnmap.first->second).currentLatitude;
	pprj_->p._plot_centre_longitude = bstBinaryRecords_.at(_it_stnmap.first->second).currentLongitude;
	
	CalculateLimitsFromPoint();
}
	

void dna_plot::CalculateLimitsFromPoint()
{	
	double temp;
	
	CDnaEllipsoid e(DEFAULT_EPSG_U);

	VincentyDirect<double>(Radians<double>(pprj_->p._plot_centre_latitude), Radians<double>(pprj_->p._plot_centre_longitude),
		Radians<double>(0.), pprj_->p._plot_area_radius, &upperDeg_, &temp, &e);
	upperDeg_ = Degrees(upperDeg_);
	
	VincentyDirect<double>(Radians<double>(pprj_->p._plot_centre_latitude), Radians<double>(pprj_->p._plot_centre_longitude),
		Radians<double>(90.), pprj_->p._plot_area_radius, &temp, &rightDeg_, &e);
	rightDeg_ = Degrees(rightDeg_);
	
	VincentyDirect<double>(Radians<double>(pprj_->p._plot_centre_latitude), Radians<double>(pprj_->p._plot_centre_longitude), 
		Radians<double>(180.), pprj_->p._plot_area_radius, &lowerDeg_, &temp, &e);
	lowerDeg_ = Degrees(lowerDeg_);
	
	VincentyDirect<double>(Radians<double>(pprj_->p._plot_centre_latitude), Radians<double>(pprj_->p._plot_centre_longitude),
		Radians<double>(270.), pprj_->p._plot_area_radius, &temp, &leftDeg_, &e);
	leftDeg_ = Degrees(leftDeg_);
}
	

void dna_plot::FormGMTDataFileNames(const UINT32& block)
{
	std::string firstPartISL = pprj_->g.output_folder + FOLDER_SLASH + network_name_;
	std::string firstPartSTN = pprj_->g.output_folder + FOLDER_SLASH + network_name_;

	firstPartSTN.append("_stn");

	if (plotBlocks_)
	{
		std::stringstream firstPart;
		std::string firstPartJSL;
		
		// add isl and jsl filenames for this block
		firstPart << firstPartSTN << "_block" << block + 1;
		firstPartISL = firstPart.str() + "_isl";
		firstPartJSL = firstPart.str() + "_jsl";
		firstPartSTN = firstPart.str();

		v_jsl_pts_file_.push_back(firstPartJSL + ".stn.d");
		v_jsl_const_file_.push_back(firstPartJSL + ".const.d");

		if (pprj_->p._plot_station_labels)
			v_jsl_lbl_file_.push_back(firstPartJSL + ".stn.lbl");
	}

	v_isl_pts_file_.push_back(firstPartISL + ".stn.d");
	v_isl_const_file_.push_back(firstPartISL + ".const.d");
	
	if (pprj_->p._plot_station_labels)
		v_isl_lbl_file_.push_back(firstPartISL + ".stn.lbl");

	if (pprj_->p._plot_positional_uncertainty)
		v_stn_apu_file_.push_back(firstPartSTN + ".apu.d");
	
	if (pprj_->p._plot_error_ellipses)
		v_stn_err_file_.push_back(firstPartSTN + ".err.d");
	
	if (pprj_->p._plot_correction_arrows)
		v_stn_cor_file_.push_back(firstPartSTN + ".cor.d");

	if (pprj_->p._plot_plate_boundaries)
		v_tectonic_plate_file_.push_back(firstPartSTN + ".plt.d");
}

void dna_plot::PrintStationsDataFileBlock(const UINT32& block)
{
	// Form file names for this block
	FormGMTDataFileNames(block);
	
	std::ofstream isl_pts, isl_const, jsl_pts, jsl_const;

	UINT32 block_index(block);
	if (plotBlocks_ && pprj_->p._plot_block_number > 0)
		block_index = 0;	// one and only block

	try {
		// Create stations data files for ISL/JSL.  Throws runtime_error on failure.
		file_opener(isl_pts, v_isl_pts_file_.at(block_index));
		file_opener(jsl_pts, v_jsl_pts_file_.at(block_index));
		// Create constraint stations data files for ISL/JSL.  Throws runtime_error on failure.
		file_opener(isl_const, v_isl_const_file_.at(block_index));
		file_opener(jsl_const, v_jsl_const_file_.at(block_index));
	}
	catch (const std::runtime_error& e) {
		SignalExceptionPlot(e.what(), 0, NULL);
	}

	// print isl coordinates
	it_vstn_t_const _it_stn;
	for (it_vUINT32 _it_isl=v_ISL_.at(block).begin();
		_it_isl<v_ISL_.at(block).end(); 
		++_it_isl)
	{
		// Advance iterator
		_it_stn = bstBinaryRecords_.begin() + *_it_isl;

		// Constrained?
		if (_it_stn->stationConst[0] == 'C' || 
			_it_stn->stationConst[1] == 'C' || 
			_it_stn->stationConst[2] == 'C')
		{
			// Constrained
			plotConstraints_ = true;
			PrintStationDataFile(isl_const, _it_stn);
		}
		else
			// Free
			PrintStationDataFile(isl_pts, _it_stn);
	}

	isl_pts.close();
	isl_const.close();

	// print jsl coordinates
	it_vUINT32 _it_jsl(v_JSL_.at(block).begin());
	for (it_vUINT32 _it_jsl=v_JSL_.at(block).begin();
		_it_jsl<v_JSL_.at(block).end(); 
		++_it_jsl)
	{
		// Advance iterator
		_it_stn = bstBinaryRecords_.begin() + *_it_jsl;

		// Constrained?
		if (_it_stn->stationConst[0] == 'C' || 
			_it_stn->stationConst[1] == 'C' || 
			_it_stn->stationConst[2] == 'C')
		{
			// Constrained
			plotConstraints_ = true;
			PrintStationDataFile(jsl_const, _it_stn);
		}
		else
			// Free
			PrintStationDataFile(jsl_pts, _it_stn);
	}

	jsl_pts.close();
	jsl_const.close();
}
	
void dna_plot::PrintStationsDataFile()
{
	// Form file names
	FormGMTDataFileNames();
	
	std::ofstream stn_pts, stn_const;
	try {
		// Create stations data file.  Throws runtime_error on failure.
		file_opener(stn_pts, v_isl_pts_file_.at(0));
		// Create constraint stations data file.  Throws runtime_error on failure.
		file_opener(stn_const, v_isl_const_file_.at(0));
	}
	catch (const std::runtime_error& e) {
		SignalExceptionPlot(e.what(), 0, NULL);
	}

	std::stringstream ss;

	for (it_vstn_t_const _it_stn(bstBinaryRecords_.begin());
		_it_stn<bstBinaryRecords_.end(); 
		++_it_stn)
	{
		// Unused Station - don't print.
		if (_it_stn->unusedStation == 1)
			continue;

		// Constrained?
		if (_it_stn->stationConst[0] == 'C' || 
			_it_stn->stationConst[1] == 'C' || 
			_it_stn->stationConst[2] == 'C')
		{
			// Constrained
			plotConstraints_ = true;
			PrintStationDataFile(stn_const, _it_stn);
		}
		else
			// Free
			PrintStationDataFile(stn_pts, _it_stn);
	}

	stn_pts.close();
	stn_const.close();

	if (!default_limits_)
	{
		if (!pprj_->p._bounding_box.empty())
			DetermineBoundingBox();
		else if (pprj_->p._plot_station_centre.empty())
			CalculateLimitsFromPoint();
		else
			CalculateLimitsFromStation();
	}
}
	

void dna_plot::PrintStationDataFile(std::ostream& os, it_vstn_t_const _it_stn)
{
	if (default_limits_) 
	{
        if (_it_stn->currentLatitude < lowerDeg_) 
			lowerDeg_ = _it_stn->currentLatitude;
        if (_it_stn->currentLatitude > upperDeg_)
            upperDeg_ = _it_stn->currentLatitude;
        if (_it_stn->currentLongitude < leftDeg_) 
			leftDeg_ = _it_stn->currentLongitude;
        if (_it_stn->currentLongitude > rightDeg_) 
			rightDeg_ = _it_stn->currentLongitude;
	}

	// print longitude and latitude
    os << std::setprecision(10) << std::fixed << 
		_it_stn->currentLongitude << "  " << _it_stn->currentLatitude << std::endl;
}

void dna_plot::PrintStationLabel(std::ostream& os, it_vstn_t_const _it_stn)
{
	if (_it_stn->unusedStation)
		return;

	// Print longitude, latitude and label
	os << std::setprecision(10) << std::fixed << 
		_it_stn->currentLongitude << "  " << // E/W
        _it_stn->currentLatitude << "   ";                                           // N/S

	// Print the label, default is station name (i.e. 200100350)
	std::string label(_it_stn->stationName);
	
	// plot alternate name?
	if (pprj_->p._plot_alt_name)
	{
		label = _it_stn->description;		// description (i.e. Acheron PM 35)
		if (trimstr(label).empty())
			label = "_";
	}
	os << label;			

	// Plot constraints?
	if (pprj_->p._plot_station_constraints)
		if (_it_stn->stationConst[0] == 'C' || 
			_it_stn->stationConst[1] == 'C' || 
			_it_stn->stationConst[2] == 'C')
				os << " (" << _it_stn->stationConst << ")";
	os << std::endl;
}

void dna_plot::PrintStationLabels()
{
	std::ofstream stn_lbl;
	try {
		// Create stations data file.  Throws runtime_error on failure.
		file_opener(stn_lbl, v_isl_lbl_file_.at(0));
	}
	catch (const std::runtime_error& e) {
		SignalExceptionPlot(e.what(), 0, NULL);
	}

	// initialise a station list and sort on longitude
	vUINT32 stnList(bstBinaryRecords_.size());
	it_vUINT32 _it_stn(stnList.begin());

	// initialise vector with 0,1,2,...,n-2,n-1,n
	initialiseIncrementingIntegerVector<UINT32>(stnList, static_cast<UINT32>(bstBinaryRecords_.size()));

	// sort stations on longitude so labels to the left are placed last
	CompareStnLongitude<station_t, UINT32> stnorderCompareFunc(&bstBinaryRecords_, false);
	std::sort(stnList.begin(), stnList.end(), stnorderCompareFunc);

	it_vstn_t _it_bstn;

	for (_it_stn=stnList.begin(); _it_stn<stnList.end(); _it_stn++)
	{
		_it_bstn = bstBinaryRecords_.begin() + *_it_stn;
		
		// print the label
		PrintStationLabel(stn_lbl, _it_bstn);
	}

	stn_lbl.close();
}
	

void dna_plot::PrintStationLabelsBlock(const UINT32& block)
{	
	UINT32 block_index(block);
	if (plotBlocks_ && pprj_->p._plot_block_number > 0)
		block_index = 0;
	
	std::ofstream isl_lbl;
	try {
		// Create stations data file.  Throws runtime_error on failure.
		file_opener(isl_lbl, v_isl_lbl_file_.at(block_index));
	}
	catch (const std::runtime_error& e) {
		SignalExceptionPlot(e.what(), 0, NULL);
	}

	// sort ISLs on longitude so labels to the left are placed last
	CompareStnLongitude<station_t, UINT32> stnorderCompareFunc(&bstBinaryRecords_, false);
	std::sort(v_ISL_.at(block).begin(), v_ISL_.at(block).end(), stnorderCompareFunc);

	it_vstn_t _it_bstn;

	for (it_vUINT32 _it_isl(v_ISL_.at(block).begin());
		_it_isl<v_ISL_.at(block).end(); 
		++_it_isl)
	{
		_it_bstn = bstBinaryRecords_.begin() + *_it_isl;
		
		// print the label
		PrintStationLabel(isl_lbl, _it_bstn);
	}

	isl_lbl.close();
	
	// return to former sort order
	std::sort(v_ISL_.at(block).begin(), v_ISL_.at(block).end());


	std::ofstream jsl_lbl;
	try {
		// Create stations data file.  Throws runtime_error on failure.
		file_opener(jsl_lbl, v_jsl_lbl_file_.at(block_index));
	}
	catch (const std::runtime_error& e) {
		SignalExceptionPlot(e.what(), 0, NULL);
	}

	// sort JSLs on longitude so labels to the left are placed last
	std::sort(v_JSL_.at(block).begin(), v_JSL_.at(block).end(), stnorderCompareFunc);

	for (it_vUINT32 _it_jsl(v_JSL_.at(block).begin());
		_it_jsl<v_JSL_.at(block).end(); 
		++_it_jsl)
	{
		_it_bstn = bstBinaryRecords_.begin() + *_it_jsl;
		
		// print the label
		PrintStationLabel(jsl_lbl, _it_bstn);	
	}
	
	jsl_lbl.close();

	// return to former sort order
	std::sort(v_JSL_.at(block).begin(), v_JSL_.at(block).end());
}
	

void dna_plot::PrintPositionalUncertainty(const UINT32& block)
{
	UINT32 block_index(block);
	if (plotBlocks_ && pprj_->p._plot_block_number > 0)
		block_index = 0;
	
	std::ofstream stn_apu;
	try {
		// Create apu file.  Throws runtime_error on failure.
		file_opener(stn_apu, v_stn_apu_file_.at(block_index));
	}
	catch (const std::runtime_error& e) {
		SignalExceptionPlot(e.what(), 0, NULL);
	}

	it_pair_string_vUINT32 _it_stnmap;
	it_vstn_t_const _it_stn(bstBinaryRecords_.begin());
	UINT32 precision(10);

	it_vstnPU_t _it_pu;

	// Print uncertainty legend, using average xx length.  See LoadPosUncertaintyFile().
	stn_apu << 
		std::setw(MSR) << std::setprecision(precision) << std::fixed << std::left << uncertainty_legend_long_ << "  " <<
		std::setw(MSR) << std::setprecision(precision) << std::fixed << std::left << uncertainty_legend_lat_ << "  " << 
		" 0 0 " <<
		std::setw(MSR) << std::setprecision(4) << std::scientific << std::right << largest_uncertainty_  << 		// semi-major
		std::setw(MSR) << std::setprecision(4) << std::scientific << std::right << largest_uncertainty_ <<			// semi-minor
		std::setw(MSR) << std::setprecision(2) << std::fixed << std::right << "45.0" << std::endl;						// orientation from e-axis
	
	for (_it_pu=v_stn_pu_.begin();
		_it_pu!=v_stn_pu_.end();
		++_it_pu)
	{
		_it_stnmap = equal_range(stnsMap_.begin(), stnsMap_.end(), _it_pu->_station, StationNameIDCompareName());

		// Not in the list? (Unlikely since adjust output will not contain stations not in the map, but test anyway)
		if (_it_stnmap.first == _it_stnmap.second)
			continue;
		
		// no need to plot constrained stations - the uncertainties will be zero!!!	
		// Horizontal only!!!
		if (bstBinaryRecords_.at(_it_stnmap.first->second).stationConst[0] == 'C' &&
			bstBinaryRecords_.at(_it_stnmap.first->second).stationConst[1] == 'C')
			continue;

		if (plotBlocks_)
			// If this station is not in this block, don't print
			if (!binary_search(v_parameterStationList_.at(block).begin(), v_parameterStationList_.at(block).end(),
				_it_stnmap.first->second))
				continue;

		stn_apu << 
			std::setw(MSR) << std::setprecision(precision) << std::fixed << std::left << _it_pu->_longitude << "  " <<
			std::setw(MSR) << std::setprecision(precision) << std::fixed << std::left << _it_pu->_latitude << "  " << 
			" 0 0 " <<
			std::setw(MSR) << std::setprecision(4) << std::scientific << std::right << _it_pu->_hzPosU * pprj_->p._pu_ellipse_scale <<
			std::setw(MSR) << std::setprecision(4) << std::scientific << std::right << _it_pu->_hzPosU * pprj_->p._pu_ellipse_scale <<
			std::setw(MSR) << std::setprecision(1) << std::fixed << std::right << "0." << std::endl;	
	}
	
	stn_apu.close();
}
	

void dna_plot::PrintErrorEllipses(const UINT32& block)
{	
	UINT32 block_index(block);
	if (plotBlocks_ && pprj_->p._plot_block_number > 0)
		block_index = 0;
	
	std::ofstream stn_err;
	try {
		// Create error ellipse file.  Throws runtime_error on failure.
		file_opener(stn_err, v_stn_err_file_.at(block_index));
	}
	catch (const std::runtime_error& e) {
		SignalExceptionPlot(e.what(), 0, NULL);
	}

	it_pair_string_vUINT32 _it_stnmap;
	it_vstn_t_const _it_stn(bstBinaryRecords_.begin());
	UINT32 precision(10);

	it_vstnPU_t _it_pu;

	double ellipse_scale(1.);

	// The error ellipse legend is printed using the largest uncertainty (semi-major) found, unless
	// positional uncertainty is being printed, in which case the error ellipse is printed using
	// half the size.  This is purely to give the impression that error ellipses are 1 sigma (68%) and
	// positional uncertainty is 95%, or close to 2 sigma.
	if (pprj_->p._plot_positional_uncertainty)
		ellipse_scale = 2.;

	// Print uncertainty legend, using average xx length.  See LoadPosUncertaintyFile().
	stn_err << 
		std::setw(MSR) << std::setprecision(precision) << std::fixed << std::left << uncertainty_legend_long_ << "  " <<
		std::setw(MSR) << std::setprecision(precision) << std::fixed << std::left << uncertainty_legend_lat_ << "  " << 
		" 0 0 " <<
		std::setw(MSR) << std::setprecision(4) << std::scientific << std::right << 
		largest_uncertainty_ / ellipse_scale << 					// semi-major
		std::setw(MSR) << std::setprecision(4) << std::scientific << std::right << 
		largest_uncertainty_ / 3.0 / ellipse_scale <<				// semi-minor
		std::setw(MSR) << std::setprecision(2) << std::fixed << std::right <<
		"45.0" << std::endl;												// orientation from e-axis
	
	double angle;

	for (_it_pu=v_stn_pu_.begin();
		_it_pu!=v_stn_pu_.end();
		++_it_pu)
	{
		_it_stnmap = equal_range(stnsMap_.begin(), stnsMap_.end(), _it_pu->_station, StationNameIDCompareName());

		// Not in the list? (Unlikely since adjust output will not contain stations not in the map, but test anyway)
		if (_it_stnmap.first == _it_stnmap.second)
			continue;
		
		// no need to plot constrained stations - the uncertainties will be zero!!!	
		// Horizontal only!!!
		if (bstBinaryRecords_.at(_it_stnmap.first->second).stationConst[0] == 'C' &&
			bstBinaryRecords_.at(_it_stnmap.first->second).stationConst[1] == 'C')
			continue;

		if (plotBlocks_)
			// If this station is not in this block, don't print
			if (!binary_search(v_parameterStationList_.at(block).begin(), v_parameterStationList_.at(block).end(),
				_it_stnmap.first->second))
				continue;

		// reduce error ellipse orientation in terms of
		// positive anti-clockwise direction from e-axis
		if (_it_pu->_orientation > 90.)
			angle = -(_it_pu->_orientation - 90.);
		else
			angle = DegtoDms(90. - DmstoDeg(_it_pu->_orientation));

		stn_err << 
			std::setw(MSR) << std::setprecision(precision) << std::fixed << std::left << _it_pu->_longitude << "  " <<
			std::setw(MSR) << std::setprecision(precision) << std::fixed << std::left << _it_pu->_latitude << "  " << 
			" 0 0 " <<
			std::setw(MSR) << std::setprecision(4) << std::scientific << std::right << _it_pu->_semimMajor * pprj_->p._pu_ellipse_scale <<
			std::setw(MSR) << std::setprecision(4) << std::scientific << std::right << _it_pu->_semimMinor * pprj_->p._pu_ellipse_scale <<
			std::setw(MSR) << std::setprecision(4) << std::fixed << std::right << angle << std::endl;	
	}
	
	stn_err.close();
}
	

void dna_plot::PrintCorrectionArrows(const UINT32& block)
{
	UINT32 block_index(block);
	if (plotBlocks_ && pprj_->p._plot_block_number > 0)
		block_index = 0;
	
	std::ofstream stn_cor;
	try {
		// Create corrections file.  Throws runtime_error on failure.
		file_opener(stn_cor, v_stn_cor_file_.at(block_index));
	}
	catch (const std::runtime_error& e) {
		SignalExceptionPlot(e.what(), 0, NULL);
	}

	it_pair_string_vUINT32 _it_stnmap;
	it_vstn_t_const _it_stn(bstBinaryRecords_.begin());
	UINT32 precision(10);

	it_vstnCor_t _it_cor;

	// Print correction legend, using average correction length.  See LoadCorrectionsFile().
	stn_cor << 
		std::setw(MSR) << std::setprecision(precision) << std::fixed << std::left << arrow_legend_long_ << "  " <<
		std::setw(MSR) << std::setprecision(precision) << std::fixed << std::left << arrow_legend_lat_ << "  " << 
		std::setw(MSR) << std::setprecision(4) << std::scientific << std::right << average_correction_ <<
		std::setw(MSR) << std::setprecision(2) << std::fixed << std::right << 0.0 <<
		" 0 0 0 " << std::endl;
	
	for (_it_cor=v_stn_corrs_.begin();
		_it_cor!=v_stn_corrs_.end();
		++_it_cor)
	{
		// Find the id of this station
		_it_stnmap = equal_range(stnsMap_.begin(), stnsMap_.end(), _it_cor->_station, StationNameIDCompareName());

		// Not in the list? (Unlikely since adjust will not output stations that are not 
		// in the map, but test anyway)
		if (_it_stnmap.first == _it_stnmap.second)
			continue;
		
		// no need to plot constrained stations - these won't have moved!	
		// Horizontal only!!!
		if (bstBinaryRecords_.at(_it_stnmap.first->second).stationConst[0] == 'C' &&
			bstBinaryRecords_.at(_it_stnmap.first->second).stationConst[1] == 'C')
			continue;

		if (plotBlocks_)
			// If this station is not in this block, don't print
			if (!binary_search(v_parameterStationList_.at(block).begin(), v_parameterStationList_.at(block).end(),
				_it_stnmap.first->second))
				continue;

		stn_cor << 
			std::setw(MSR) << std::setprecision(precision) << std::fixed << std::left << bstBinaryRecords_.at(_it_stnmap.first->second).currentLongitude << "  " <<
			std::setw(MSR) << std::setprecision(precision) << std::fixed << std::left << bstBinaryRecords_.at(_it_stnmap.first->second).currentLatitude << "  " << 
			std::setw(MSR) << std::setprecision(4) << std::scientific << std::right << _it_cor->_east * pprj_->p._correction_scale <<
			std::setw(MSR) << std::setprecision(4) << std::scientific << std::right << _it_cor->_north * pprj_->p._correction_scale <<
			" 0 0 0 ";

		if (pprj_->p._plot_correction_labels)
			stn_cor << 
				std::setw(HEIGHT) << std::setprecision(3) << std::fixed << std::right << _it_cor->_east << "e," <<
				std::setprecision(3) << std::fixed << _it_cor->_north << "n   ";
		stn_cor << std::endl;	
	}
	
	stn_cor.close();
}
	
void dna_plot::PrintPlateBoundaries(const UINT32& block)
{
	// Load tectonic plate boundaries file and export
	dna_io_tpb tpb;

	// Tectonic plate boundaries
	v_string_v_doubledouble_pair global_plates;				

	std::stringstream ss;
	ss << "PrintPlateBoundaries(): An error was encountered when loading tectonic plate information." << std::endl;

	try {
		// Load tectonic plate polygons.  Throws runtime_error on failure.
		tpb.load_tpb_file(pprj_->r.tpb_file, global_plates);
		std::sort(global_plates.begin(), global_plates.end());
	}
	catch (const std::runtime_error& e) {
		ss << e.what();
		throw std::runtime_error(ss.str());
	}

	UINT32 block_index(block);
	if (plotBlocks_ && pprj_->p._plot_block_number > 0)
		block_index = 0;
	
	std::ofstream tectonic_plate;
	try {
		// Create corrections file.  Throws runtime_error on failure.
		file_opener(tectonic_plate, v_tectonic_plate_file_.at(block_index));
	}
	catch (const std::runtime_error& e) {
		SignalExceptionPlot(e.what(), 0, NULL);
	}

	for_each(global_plates.begin(), global_plates.end(),
		[&tectonic_plate](const string_v_doubledouble_pair& plate){
			tectonic_plate << ">" << std::endl;

			for_each(plate.second.begin(), plate.second.end(),
			[&tectonic_plate](const doubledouble_pair& coordinate){
				tectonic_plate << std::fixed << std::setprecision(3) <<
					std::setw(10) << std::right << coordinate.first << "," <<
					std::setw(10) << std::right << coordinate.second << std::endl;
			});
		});	

	tectonic_plate.close();
}

// Plots the specified measurement types to individual files
// All other measurement types (_combined_msr_list) are printed to a single file
// This enables colours to be specified for individual measurement types
void dna_plot::PrintMeasurementsDatFilesBlock(const UINT32& block)
{
	v_msr_file_.at(block).clear();

	std::string type;
	std::stringstream ss, filename;

	std::ofstream msr_file_stream;
	
	for (UINT16 c=0; c<pprj_->p._separate_msrs.size(); c++)
	{
		// no measurement for this type?
		if (parsemsrTally_.MeasurementCount(pprj_->p._separate_msrs.at(c)) == 0)
			continue;

		ss.str("");
		ss << pprj_->p._separate_msrs.at(c);
		type = ss.str();
		filename.str("");
		filename << output_folder_ << FOLDER_SLASH << network_name_ << "_msr_block" << block + 1 << "_" << pprj_->p._separate_msrs.at(c) << ".d";
		v_msr_def_file_.push_back(filename.str());
		v_msr_file_.at(block).push_back(string_string_pair(v_msr_def_file_.back(), type));

		try {
			// Create msr data file.  Throws runtime_error on failure.
			file_opener(msr_file_stream, v_msr_def_file_.back());
		}
		catch (const std::runtime_error& e) {
			SignalExceptionPlot(e.what(), 0, NULL);
		}

		PrintMeasurementsDatFileBlock(block, pprj_->p._separate_msrs.at(c), &msr_file_stream);

		msr_file_stream.close();

		// remove this measurement type from the global list
		it_vchar_range = equal_range(_combined_msr_list.begin(), _combined_msr_list.end(), 
			pprj_->p._separate_msrs.at(c));

		if (it_vchar_range.first != it_vchar_range.second)
		{
			// found - remove it from the list.  The remaining elements will be printed to a single file
			_combined_msr_list.erase(it_vchar_range.first);
		}
	}

	std::vector<char>::const_iterator _it_type;
	bool printOtherTypes(false);
	for (_it_type=_combined_msr_list.begin(); _it_type!=_combined_msr_list.end(); ++_it_type)
	{
		if (parsemsrTally_.MeasurementCount(*_it_type) > 0)
		{
			printOtherTypes = true;
			break;
		}
	}

	// don't bother printing this file if there are no 'other' measurements
	if (!printOtherTypes)
		return;

	// now print remaining measurement types to a single file
	filename.str("");
	filename << output_folder_ << FOLDER_SLASH << network_name_ << "_msr_block" << block + 1 << "_";
	ss.str("");
	for (_it_type=_combined_msr_list.begin(); _it_type!=_combined_msr_list.end(); _it_type++)
		ss << *_it_type;
	filename << ss.str() << ".d";
	v_msr_all_file_.push_back(filename.str());
	v_msr_file_.at(block).insert(v_msr_file_.at(block).begin(), string_string_pair(v_msr_all_file_.back(), ss.str()));

	try {
		// Create msr data file.  Throws runtime_error on failure.
		file_opener(msr_file_stream, v_msr_all_file_.back());
	}
	catch (const std::runtime_error& e) {
		SignalExceptionPlot(e.what(), 0, NULL);
	}

	for (_it_type=_combined_msr_list.begin(); _it_type!=_combined_msr_list.end(); _it_type++)
		PrintMeasurementsDatFileBlock(block, *_it_type, &msr_file_stream);

	msr_file_stream.close();

}


// Plots the specified measurement types to individual files
// All other measurement types (_combined_msr_list) are printed to a single file
// This enables colours to be specified for individual measurement types
void dna_plot::PrintMeasurementsDatFiles()
{
	_combined_msr_list.clear();
	_combined_msr_list.push_back('A'); // Horizontal angle
	_combined_msr_list.push_back('B'); // Geodetic azimuth
	_combined_msr_list.push_back('C'); // Chord dist
	_combined_msr_list.push_back('D'); // Direction set
	_combined_msr_list.push_back('E'); // Ellipsoid arc
	_combined_msr_list.push_back('G'); // GPS Baseline (treat as single-baseline cluster)
	_combined_msr_list.push_back('H'); // Orthometric height
	_combined_msr_list.push_back('I'); // Astronomic latitude
	_combined_msr_list.push_back('J'); // Astronomic longitude
	_combined_msr_list.push_back('K'); // Astronomic azimuth
	_combined_msr_list.push_back('L'); // Level difference
	_combined_msr_list.push_back('M'); // MSL arc
	_combined_msr_list.push_back('P'); // Geodetic latitude
	_combined_msr_list.push_back('Q'); // Geodetic longitude
	_combined_msr_list.push_back('R'); // Ellipsoidal height
	_combined_msr_list.push_back('S'); // Slope distance
	_combined_msr_list.push_back('V'); // Zenith distance
	_combined_msr_list.push_back('X'); // GPS Baseline cluster
	_combined_msr_list.push_back('Y'); // GPS point cluster
	_combined_msr_list.push_back('Z'); // Vertical angle

	//vector<char>::iterator _it_type(_combined_msr_list.end());
	
	std::string msr_file_name;
	v_msr_file_.at(0).clear();

	YClusterStations_.clear();

	std::string type;
	std::stringstream ss;

	std::ofstream msr_file_stream;

	for (UINT16 c=0; c<pprj_->p._separate_msrs.size(); c++)
	{
		// no measurement for this type?
		if (parsemsrTally_.MeasurementCount(pprj_->p._separate_msrs.at(c)) == 0)
			continue;

		ss.str("");
		ss << pprj_->p._separate_msrs.at(c);
		type = ss.str();
		msr_file_name = output_folder_ + FOLDER_SLASH + network_name_ + "_msr_" + pprj_->p._separate_msrs.at(c) + ".d";
		v_msr_file_.at(0).push_back(string_string_pair(msr_file_name, type));

		try {
			// Create msr data file.  Throws runtime_error on failure.
			file_opener(msr_file_stream, msr_file_name);
		}
		catch (const std::runtime_error& e) {
			SignalExceptionPlot(e.what(), 0, NULL);
		}

		PrintMeasurementsDatFile(pprj_->p._separate_msrs.at(c), &msr_file_stream);

		msr_file_stream.close();

		// remove this measurement type from the global list
		it_vchar_range = equal_range(_combined_msr_list.begin(), _combined_msr_list.end(), 
			pprj_->p._separate_msrs.at(c));

		if (it_vchar_range.first != it_vchar_range.second)
		{
			// found - remove it from the list.  The remaining elements will be printed to a single file
			_combined_msr_list.erase(it_vchar_range.first);
		}
	}

	std::vector<char>::const_iterator _it_type;
	bool printOtherTypes(false);
	for (_it_type=_combined_msr_list.begin(); _it_type!=_combined_msr_list.end(); ++_it_type)
	{
		if (parsemsrTally_.MeasurementCount(*_it_type) > 0)
		{
			printOtherTypes = true;
			break;
		}
	}

	// don't bother printing this file if there are no 'other' measurements
	if (!printOtherTypes)
		return;

	// now print remaining measurement types to a single file
	msr_file_name = output_folder_ + FOLDER_SLASH + network_name_ + "_msr_";
	ss.str("");
	for (_it_type=_combined_msr_list.begin(); _it_type!=_combined_msr_list.end(); ++_it_type)
		ss << *_it_type;
	type = ss.str();
	msr_file_name.append(type).append(".d");
	v_msr_file_.at(0).insert(v_msr_file_.at(0).begin(), string_string_pair(msr_file_name, type));

	try {
		// Create msr data file.  Throws runtime_error on failure.
		file_opener(msr_file_stream, msr_file_name);
	}
	catch (const std::runtime_error& e) {
		SignalExceptionPlot(e.what(), 0, NULL);
	}

	for (_it_type=_combined_msr_list.begin(); _it_type!=_combined_msr_list.end(); _it_type++)
		PrintMeasurementsDatFile(*_it_type, &msr_file_stream);

	msr_file_stream.close();

}


bool dna_plot::WithinLimits(const double& latitude, const double& longitude)
{
	if (latitude >= lowerDeg_ && latitude <= upperDeg_ &&
		longitude >= leftDeg_ && longitude <= rightDeg_)
		return true;
	else
		return false;
}
	

void dna_plot::PrintMeasurementsDatFileBlock(const UINT32& block, char msrType, std::ofstream* msr_file_stream)
{
	it_vstn_t _it_stn(bstBinaryRecords_.begin());
	
	UINT32 precision(10);
	bool FirstisWithinLimits, SecondisWithinLimits, ThirdisWithinLimits;

	std::stringstream ss;
	
	it_vUINT32 _it_block_msr(v_CML_.at(block).begin());
	vUINT32 msrIndices;
	it_vmsr_t _it_msr;

	(*msr_file_stream) << ">" << std::endl;

	for (_it_block_msr=v_CML_.at(block).begin(); _it_block_msr<v_CML_.at(block).end(); _it_block_msr++)
	{
		_it_msr = bmsBinaryRecords_.begin() + (*_it_block_msr);

		if (msrType != _it_msr->measType)
			continue;

		if (_it_msr->ignore && !pprj_->p._plot_ignored_msrs)
			continue;

		ss.str("");
		ss.clear();
		FirstisWithinLimits = SecondisWithinLimits = ThirdisWithinLimits = false;
		switch (msrType)
		{
		case 'A':	// Horizontal angles
			// Print 3 - 1 - 2
			//
			// Third
			ss << std::setprecision(precision) << std::fixed << bstBinaryRecords_.at(_it_msr->station3).currentLongitude;
			ss << "  ";
			ss << std::setprecision(precision) << std::fixed << bstBinaryRecords_.at(_it_msr->station3).currentLatitude;
			ss << std::endl;

			// test if point is within limits of custom defined box
			if (!default_limits_)
				ThirdisWithinLimits = WithinLimits(bstBinaryRecords_.at(_it_msr->station3).currentLatitude, bstBinaryRecords_.at(_it_msr->station3).currentLongitude);

			[[fallthrough]];

		case 'G': // GPS Baseline
		case 'X': // GPS Baseline cluster
		
			// The following prevents GPS measurements from printing three times (i.e. once per X, Y, Z)
			if (_it_msr->measStart > xMeas)
				continue;

			[[fallthrough]];
			
		case 'B': // Geodetic azimuth
		case 'D': // Direction set
		case 'K': // Astronomic azimuth
		case 'C': // Chord dist
		case 'E': // Ellipsoid arc
		case 'M': // MSL arc
		case 'S': // Slope distance
		case 'L': // Level difference
		case 'V': // Zenith distance
		case 'Z': // Vertical angle
			// Get the measurement indices involved in this measurement
			GetMsrIndices<UINT32>(bmsBinaryRecords_, *_it_block_msr, msrIndices);

			for (it_vUINT32 msr=msrIndices.begin(); msr!=msrIndices.end(); ++msr)
			{
				if (bmsBinaryRecords_.at(*msr).measStart > xMeas)
					continue;

				if (!default_limits_)
				{
					FirstisWithinLimits = WithinLimits(
						bstBinaryRecords_.at(bmsBinaryRecords_.at(*msr).station1).currentLatitude, 
						bstBinaryRecords_.at(bmsBinaryRecords_.at(*msr).station1).currentLongitude);
					SecondisWithinLimits = WithinLimits(
						bstBinaryRecords_.at(bmsBinaryRecords_.at(*msr).station2).currentLatitude, 
						bstBinaryRecords_.at(bmsBinaryRecords_.at(*msr).station2).currentLongitude);

					if (!FirstisWithinLimits && !SecondisWithinLimits && !ThirdisWithinLimits)
						continue;
				}

				// Third.  If case 'A' was reached before case 'B', ss will be non-empty
				(*msr_file_stream) << ss.str();
			
				// First
				(*msr_file_stream) << std::setprecision(precision) << std::fixed << bstBinaryRecords_.at(bmsBinaryRecords_.at(*msr).station1).currentLongitude;
				(*msr_file_stream) << "  ";
				(*msr_file_stream) << std::setprecision(precision) << std::fixed << bstBinaryRecords_.at(bmsBinaryRecords_.at(*msr).station1).currentLatitude;
				(*msr_file_stream) << std::endl;
			
				// Second
				(*msr_file_stream) << std::setprecision(precision) << std::fixed << bstBinaryRecords_.at(bmsBinaryRecords_.at(*msr).station2).currentLongitude;
				(*msr_file_stream) << "  ";
				(*msr_file_stream) << std::setprecision(precision) << std::fixed << bstBinaryRecords_.at(bmsBinaryRecords_.at(*msr).station2).currentLatitude;
				(*msr_file_stream) << std::endl << ">" << std::endl;
			}
			break;
		
		// single station
		case 'H': // Orthometric height
		case 'I': // Astronomic latitude
		case 'J': // Astronomic longitude
		case 'P': // Geodetic latitude
		case 'Q': // Geodetic longitude
		case 'R': // Ellipsoidal height
			if (!default_limits_)
				if (!WithinLimits(
					bstBinaryRecords_.at(_it_msr->station1).currentLatitude, 
					bstBinaryRecords_.at(_it_msr->station1).currentLongitude))
					continue;

			// First
			(*msr_file_stream) << std::setprecision(precision) << std::fixed << bstBinaryRecords_.at(_it_msr->station1).currentLongitude;
			(*msr_file_stream) << "  ";
			(*msr_file_stream) << std::setprecision(precision) << std::fixed << bstBinaryRecords_.at(_it_msr->station1).currentLatitude;
			(*msr_file_stream) << std::endl << ">" << std::endl;
			
			break;

		case 'Y': // GPS point cluster
			
			// Get the measurement indices involved in this measurement
			GetMsrIndices<UINT32>(bmsBinaryRecords_, *_it_block_msr, msrIndices);

			for (it_vUINT32 msr=msrIndices.begin(); msr!=msrIndices.end(); ++msr)
			{
				if (bmsBinaryRecords_.at(*msr).measStart > xMeas)
					continue;

				if (!default_limits_)
					if (!WithinLimits(
						bstBinaryRecords_.at(bmsBinaryRecords_.at(*msr).station1).currentLatitude, 
						bstBinaryRecords_.at(bmsBinaryRecords_.at(*msr).station1).currentLongitude))
						continue;

				// First
				(*msr_file_stream) << std::setprecision(precision) << std::fixed << bstBinaryRecords_.at(bmsBinaryRecords_.at(*msr).station1).currentLongitude;
				(*msr_file_stream) << "  ";
				(*msr_file_stream) << std::setprecision(precision) << std::fixed << bstBinaryRecords_.at(bmsBinaryRecords_.at(*msr).station1).currentLatitude;
				(*msr_file_stream) << std::endl << ">" << std::endl;
			}
		}
	}

}
	

void dna_plot::PrintMeasurementsDatFile(char msrType, std::ofstream* msr_file_stream)
{
	it_vstn_t _it_stn(bstBinaryRecords_.begin());
	it_vmsr_t _it_msr(bmsBinaryRecords_.begin());
	
	UINT32 precision(10);
	bool FirstisWithinLimits, SecondisWithinLimits, ThirdisWithinLimits;

	std::stringstream ss;

	(*msr_file_stream) << ">" << std::endl;
	
	for (_it_msr=bmsBinaryRecords_.begin(); _it_msr!=bmsBinaryRecords_.end(); _it_msr++)
	{
		if (msrType != _it_msr->measType)
			continue;
		
		if (_it_msr->ignore && !pprj_->p._plot_ignored_msrs)
			continue;
		
		ss.str("");
		ss.clear();
		FirstisWithinLimits = SecondisWithinLimits = ThirdisWithinLimits = false;
		switch (msrType)
		{
		case 'A':	// Horizontal angles
			// Print 3 - 1 - 2
			//
			// Third
			ss << std::setprecision(precision) << std::fixed << bstBinaryRecords_.at(_it_msr->station3).currentLongitude;
			ss << "  ";
			ss << std::setprecision(precision) << std::fixed << bstBinaryRecords_.at(_it_msr->station3).currentLatitude;
			ss << std::endl;

			// test if point is within limits of custom defined box
			if (!default_limits_)
				ThirdisWithinLimits = WithinLimits(bstBinaryRecords_.at(_it_msr->station3).currentLatitude, bstBinaryRecords_.at(_it_msr->station3).currentLongitude);

			[[fallthrough]];

		case 'G': // GPS Baseline
		case 'X': // GPS Baseline cluster
		
			// The following prevents GPS measurements from printing three times (i.e. once per X, Y, Z)
			if (_it_msr->measStart > xMeas)
				continue;

			[[fallthrough]];

		case 'B': // Geodetic azimuth
		case 'D': // Direction set
		case 'K': // Astronomic azimuth
		case 'C': // Chord dist
		case 'E': // Ellipsoid arc
		case 'M': // MSL arc
		case 'S': // Slope distance
		case 'L': // Level difference
		case 'V': // Zenith distance
		case 'Z': // Vertical angle

			if (!default_limits_)
			{
				FirstisWithinLimits = WithinLimits(bstBinaryRecords_.at(_it_msr->station1).currentLatitude, bstBinaryRecords_.at(_it_msr->station1).currentLongitude);
				SecondisWithinLimits = WithinLimits(bstBinaryRecords_.at(_it_msr->station2).currentLatitude, bstBinaryRecords_.at(_it_msr->station2).currentLongitude);

				if (!FirstisWithinLimits && !SecondisWithinLimits && !ThirdisWithinLimits)
					continue;
			}

			// Third 
			// If this measurement is an angle, then ss will contain data for third station
			// Otherwise, ss will be blank
			(*msr_file_stream) << ss.str();
			
			// First
			(*msr_file_stream) << std::setprecision(precision) << std::fixed << bstBinaryRecords_.at(_it_msr->station1).currentLongitude;
			(*msr_file_stream) << "  ";
			(*msr_file_stream) << std::setprecision(precision) << std::fixed << bstBinaryRecords_.at(_it_msr->station1).currentLatitude;
			(*msr_file_stream) << std::endl;
			
			// Second
			(*msr_file_stream) << std::setprecision(precision) << std::fixed << bstBinaryRecords_.at(_it_msr->station2).currentLongitude;
			(*msr_file_stream) << "  ";
			(*msr_file_stream) << std::setprecision(precision) << std::fixed << bstBinaryRecords_.at(_it_msr->station2).currentLatitude;
			(*msr_file_stream) << std::endl << ">" << std::endl;
			break;
		
		// single station
		case 'H': // Orthometric height
		case 'I': // Astronomic latitude
		case 'J': // Astronomic longitude
		case 'P': // Geodetic latitude
		case 'Q': // Geodetic longitude
		case 'R': // Ellipsoidal height
			
			if (!default_limits_)
			{
				FirstisWithinLimits = WithinLimits(bstBinaryRecords_.at(_it_msr->station1).currentLatitude, bstBinaryRecords_.at(_it_msr->station1).currentLongitude);
				if (!FirstisWithinLimits)
					continue;
			}

			// First
			(*msr_file_stream) << std::setprecision(precision) << std::fixed << bstBinaryRecords_.at(_it_msr->station1).currentLongitude;
			(*msr_file_stream) << "  ";
			(*msr_file_stream) << std::setprecision(precision) << std::fixed << bstBinaryRecords_.at(_it_msr->station1).currentLatitude;
			(*msr_file_stream) << std::endl << ">" << std::endl;			
			break;

		case 'Y': // GPS point cluster
			
			// The following prevents GPS measurements from printing three times (i.e. once per X, Y, Z)
			if (_it_msr->measStart > xMeas)
				continue;

			//// capture first station in the cluster
			//if (_it_tmp == bmsBinaryRecords_.end())
			//{
			//	_it_tmp = _it_msr;
			//	continue;
			//}

			//// Is this a new measurement cluster?
			//if (_it_tmp->clusterID != _it_msr->clusterID)
			//{
			//	_it_tmp = _it_msr;
			//	continue;
			//}

			if (binary_search(YClusterStations_.begin(), YClusterStations_.end(), _it_msr->station1))
				// already have this one
				continue;
			
			YClusterStations_.push_back(_it_msr->station1);
			std::sort(YClusterStations_.begin(), YClusterStations_.end());

			// First
			(*msr_file_stream) << std::setprecision(precision) << std::fixed << bstBinaryRecords_.at(_it_msr->station1).currentLongitude;
			(*msr_file_stream) << "  ";
			(*msr_file_stream) << std::setprecision(precision) << std::fixed << bstBinaryRecords_.at(_it_msr->station1).currentLatitude;
			(*msr_file_stream) << std::endl << ">" << std::endl;

		}
	}

}
	

void dna_plot::PrintGMTParameters()
{
	it_string_pair _it_vpstr(pprj_->p._gmt_params.begin());

	try {
		for (; _it_vpstr!=pprj_->p._gmt_params.end(); _it_vpstr++)
			gmtbat_file_ << _APP_GMTSET_ << " " << _it_vpstr->first << " " << _it_vpstr->second << std::endl;
		gmtbat_file_ << std::endl;
	}
	catch (const std::ios_base::failure& f) {
		SignalExceptionPlot(static_cast<std::string>(f.what()), 0, "o", &gmtbat_file_);	
	}
	
}
	

void dna_plot::LoadNetworkFiles(const project_settings& projectSettings)
{
	projectSettings_ = projectSettings;

	_combined_msr_list.clear();
	MsrTally::FillMsrList(_combined_msr_list);

	blockCount_ = 1;
	plotBlocks_ = projectSettings.p._plot_phased_blocks;
	LoadBinaryFiles();
	LoadStationMap();
	if (plotBlocks_)
		LoadSegmentationFile();

	// Printing uncertainties?
	if (projectSettings.p._plot_error_ellipses || 
		projectSettings.p._plot_positional_uncertainty)
	{
		uncertainty_legend_length_ = 0.5;
		LoadPosUncertaintyFile();
	}
	
	// Printing corrections?
	if (projectSettings.p._plot_correction_arrows)
	{
		arrow_legend_length_ = 3.;
		if (projectSettings.p._compute_corrections)
			ComputeStationCorrections();
		else
			LoadCorrectionsFile();
	}

	if (plotBlocks_)
		// Parameter station list is needed for correction arrows
		BuildParameterStationList();
}


void dna_plot::LoadBinaryFiles()
{
	try {
		// Load binary stations data.  Throws runtime_error on failure.
		BstFile bst;
		stationCount_ = bst.LoadFile(projectSettings_.i.bst_file, &bstBinaryRecords_, bst_meta_);

		// Load binary stations data.  Throws runtime_error on failure.
		BmsFile bms;
		bms.LoadFile(projectSettings_.i.bms_file, &bmsBinaryRecords_, bms_meta_);
	}
	catch (const std::runtime_error& e) {
		SignalExceptionPlot(e.what(), 0, NULL);
	}

	reference_frame_ = projectSettings_.i.reference_frame;
	
	try {
		reference_frame_ = datumFromEpsgString<std::string>(bst_meta_.epsgCode);
	}
	catch (...) {
		// do nothing
	}
	
	for_each(
		bstBinaryRecords_.begin(), 
		bstBinaryRecords_.end(), 
		[this] (station_t& stn) {			// use lambda expression
			ReduceStationCoordinates(&stn);
		}
	);

	// Create measurement tally to determine the measurement
	// types that have been supplied
	for_each(
		bmsBinaryRecords_.begin(), 
		bmsBinaryRecords_.end(), 
		[this] (measurement_t& msr) {			// use lambda expression
			
			switch (msr.measType)
			{
			case 'A': // Horizontal angle
				parsemsrTally_.A++;
				break;
			case 'B': // Geodetic azimuth
				parsemsrTally_.B++;
				break;
			case 'C': // Chord dist
				parsemsrTally_.C++;
				break;
			case 'D': // Direction set
				if (msr.vectorCount1 > 0)
					parsemsrTally_.D += msr.vectorCount1 - 1;
				break;
			case 'E': // Ellipsoid arc
				parsemsrTally_.E++;
				break;
			case 'G': // GPS Baseline
				parsemsrTally_.G++;
				break;
			case 'H': // Orthometric height
				parsemsrTally_.H++;
				break;
			case 'I': // Astronomic latitude
				parsemsrTally_.I++;
				break;
			case 'J': // Astronomic longitude
				parsemsrTally_.J++;
				break;
			case 'K': // Astronomic azimuth
				parsemsrTally_.K++;
				break;
			case 'L': // Level difference
				parsemsrTally_.L++;
				break;
			case 'M': // MSL arc
				parsemsrTally_.M++;
				break;
			case 'P': // Geodetic latitude
				parsemsrTally_.P++;
				break;
			case 'Q': // Geodetic longitude
				parsemsrTally_.Q++;
				break;
			case 'R': // Ellipsoidal height
				parsemsrTally_.R++;
				break;
			case 'S': // Slope distance
				parsemsrTally_.S++;
				break;
			case 'V': // Zenith distance
				parsemsrTally_.V++;
				break;
			case 'X': // GPS Baseline cluster
				if (msr.measStart == xMeas)
				{
					if (msr.vectorCount1 == msr.vectorCount2 + 1)
						parsemsrTally_.X += (msr.vectorCount1 * 3);
				}
				break;
			case 'Y': // GPS point cluster
				if (msr.measStart == xMeas)
				{
					if (msr.vectorCount1 == msr.vectorCount2 + 1)
						parsemsrTally_.Y += (msr.vectorCount1 * 3);
				}
				break;
			case 'Z': // Vertical angle
				parsemsrTally_.Z++;
				break;
			}
		}
	);

	measurementCount_ = parsemsrTally_.TotalCount();	

}


void dna_plot::LoadStationMap()
{
	try {
		// Load station map.  Throws runtime_error on failure.
		dynadjust::iostreams::MapFile map;
		map.LoadFile(projectSettings_.i.map_file, &stnsMap_);
	}
	catch (const std::runtime_error& e) {
		SignalExceptionPlot(e.what(), 0, NULL);
	}
}
	

//void dna_plot::SortandMapStations()
//{
//	UINT32 stnCount(static_cast<UINT32>(bstBinaryRecords_.size()));
//	stnsMap_.clear();
//	stnsMap_.reserve(stnCount);
//
//	// Create the Station-Name / ID map
//	string_uint32_pair stnID;
//	for (UINT32 stnIndex=0; stnIndex<stnCount; stnIndex++)
//	{
//		stnID.first = bstBinaryRecords_.at(stnIndex).stationName;
//		stnID.second = stnIndex;
//		stnsMap_.push_back(stnID);
//	}
//
//	// sort on station name (i.e. first of the pair)
//	std::sort(stnsMap_.begin(), stnsMap_.end(), StationNameIDCompareName());
//
//	if (stnsMap_.size() < stnCount)
//		SignalExceptionPlot("SortandMapStations(): Could not allocate sufficient memory for the Station map.", 0, NULL);
//}
	

void dna_plot::LoadSegmentationFile()
{
	vUINT32 v_ContiguousNetList, v_parameterStationCount;
	try {
		// Load segmentation file.  Throws runtime_error on failure.
		SegFile seg;
		seg.LoadSegFile(projectSettings_.s.seg_file, 
			blockCount_, blockThreshold_, minInnerStns_,
			v_ISL_, v_JSL_, v_CML_,
			true, &bmsBinaryRecords_,
			&v_measurementCount_, &v_unknownsCount_, &v_ContiguousNetList,
			&v_parameterStationCount);
	}
	catch (const std::runtime_error& e) {
		SignalExceptionPlot(e.what(), 0, NULL);
	}

	if (projectSettings_.p._plot_block_number > blockCount_)
	{
		std::stringstream ss;
		ss << 
			"Specified block number (" << projectSettings_.p._plot_block_number << 
			") exceeds the total number of blocks (" <<
			blockCount_ << ")." << std::endl;
		SignalExceptionPlot(ss.str(), 0, NULL);	
	}

	v_msr_file_.resize(blockCount_);
}
	

// The positional uncertainty file is produced by adjust.
void dna_plot::LoadPosUncertaintyFile()
{
	// Load corrections file generated by adjust
	std::ifstream apu_file;
	try {
		// Load apu file.  Throws runtime_error on failure.
		file_opener(apu_file, projectSettings_.o._apu_file, std::ios::in, ascii, true);
	}
	catch (const std::runtime_error& e) {
		SignalExceptionPlot(e.what(), 0, NULL);
	}

	bool fullCovarianceMatrix(false);
	UINT32 block(0), stn(0), stn_cov, blockstnCount;
	std::string strLine;
	stationPosUncertainty_t posUnc;

	v_stn_pu_.resize(blockCount_);
	largest_uncertainty_ = 0.0;

	char line[PRINT_LINE_LENGTH];
	bool dataBlocks(false);
	it_pair_string_vUINT32 _it_stnmap;

	APU_UNITS_UI vcv_units;
	matrix_2d vcv_cart(3, 3), vcv_local(3, 3);

	try {
		
		// read header line
		apu_file.ignore(PRINT_LINE_LENGTH, '\n');		// ------------------------
		apu_file.ignore(PRINT_LINE_LENGTH, '\n');		//            DYNADJUST POSITIONAL UNCERTAINTY...
		apu_file.ignore(PRINT_LINE_LENGTH, '\n');		// 
		apu_file.ignore(PRINT_LINE_LENGTH, '\n');		// Version
		apu_file.ignore(PRINT_LINE_LENGTH, '\n');		// Build
		apu_file.ignore(PRINT_LINE_LENGTH, '\n');		// File created
		apu_file.ignore(PRINT_LINE_LENGTH, '\n');		// File name
		apu_file.ignore(PRINT_LINE_LENGTH, '\n');		// 
		apu_file.ignore(PRINT_LINE_LENGTH, '\n');		// PU Confidence interval
		apu_file.ignore(PRINT_LINE_LENGTH, '\n');		// Error ellipse axes
		apu_file.ignore(PRINT_LINE_LENGTH, '\n');		// Variances

		apu_file.getline(line, PRINT_LINE_LENGTH);		// Stations printed in blocks
		strLine = trimstr(std::string(line));		
		if (!iequals(strLine.substr(0, 16), "Stations printed"))
		{
			std::stringstream ss;
			ss << "LoadPosUncertaintyFile(): " << projectSettings_.o._apu_file << " is corrupt." << std::endl;
			ss << "  Expected to find 'Stations printed in blocks' field." << std::endl;
			SignalExceptionPlot(ss.str(), 0, "i", &apu_file);
		}
		
		if ((dataBlocks = yesno_uint<bool, std::string>(strLine.substr(PRINT_VAR_PAD))))
			v_stn_pu_.resize(bstBinaryRecords_.size());
		else
			v_stn_pu_.reserve(bstBinaryRecords_.size());

		apu_file.getline(line, PRINT_LINE_LENGTH);		// Variance matrix units
		strLine = trimstr(std::string(line));
		if (iequals(strLine.substr(PRINT_VAR_PAD, 3), "XYZ"))
			vcv_units = XYZ_apu_ui;
		else
			vcv_units = ENU_apu_ui;
		

		apu_file.getline(line, PRINT_LINE_LENGTH);		// Full covariance matrix
		strLine = trimstr(std::string(line));
		if (!iequals(strLine.substr(0, 15), "Full covariance"))
		{
			std::stringstream ss;
			ss << "LoadPosUncertaintyFile(): " << projectSettings_.o._apu_file << " is corrupt." << std::endl;
			ss << "  Expected to find 'Full covariance matrix' field." << std::endl;
			SignalExceptionPlot(ss.str(), 0, "i", &apu_file);
		}
		
		fullCovarianceMatrix = yesno_uint<bool, std::string>(strLine.substr(PRINT_VAR_PAD));

		// If covariances are printed, then the data will appear in blocks
		if (fullCovarianceMatrix)
			dataBlocks = true;
		
		apu_file.ignore(PRINT_LINE_LENGTH, '\n');		// ------------------------
		apu_file.ignore(PRINT_LINE_LENGTH, '\n');		// 
		apu_file.ignore(PRINT_LINE_LENGTH, '\n');		// Positional uncertainty of ...
		apu_file.ignore(PRINT_LINE_LENGTH, '\n');		// ------------------------
		
		for (block=0; block<blockCount_; ++block)
		{
			if (dataBlocks)
			{
				blockstnCount = static_cast<UINT32>(v_ISL_.at(block).size() + v_JSL_.at(block).size());

				apu_file.ignore(PRINT_LINE_LENGTH, '\n');		// 
				apu_file.getline(line, PRINT_LINE_LENGTH);		// Block #

				strLine = trimstr(std::string(line));
				
				if (!iequals(strLine.substr(0, 5), "Block"))
				{
					std::stringstream ss;
					ss << "LoadPosUncertaintyFile(): " << projectSettings_.o._apu_file << " is corrupt." << std::endl;
					ss << "  Expected to read Block data, but found " << strLine << std::endl;
					SignalExceptionPlot(ss.str(), 0, "i", &apu_file);
				}

				if (LongFromString<UINT32>(strLine.substr(6)) != block + 1)
				{
					std::stringstream ss;
					ss << "LoadPosUncertaintyFile(): " << projectSettings_.o._apu_file << " is corrupt." << std::endl;
					ss << "  Expected to read data for block " << strLine.substr(6) << ", but found " << strLine << std::endl;
					SignalExceptionPlot(ss.str(), 0, "i", &apu_file);
				}
			}
			else
			{
				blockstnCount = stationCount_;
				apu_file.ignore(PRINT_LINE_LENGTH, '\n');		// 
			}
			
			apu_file.ignore(PRINT_LINE_LENGTH, '\n');		// Station
			apu_file.ignore(PRINT_LINE_LENGTH, '\n');		// ------------------------
			
			for (stn=0; stn<blockstnCount;  ++stn)
			{
				apu_file.getline(line, PRINT_LINE_LENGTH);		// Now the data...

				strLine = trimstr(std::string(line));
				if (strLine.length() < STATION)
				{
					if (apu_file.eof())
						break;				
					continue;
				}

				posUnc._station =	  trimstr(strLine.substr(0, STATION));
				posUnc._latitude =	  DmstoDeg(DoubleFromString<double>(trimstr(strLine.substr(STATION+PAD2, LAT_EAST))));
				posUnc._longitude =   DmstoDeg(DoubleFromString<double>(trimstr(strLine.substr(STATION+PAD2+LAT_EAST, LON_NORTH))));
				posUnc._hzPosU =	  DoubleFromString<double>(trimstr(strLine.substr(STATION+PAD2+LAT_EAST+LON_NORTH, STAT)));
				posUnc._vtPosU =	  DoubleFromString<double>(trimstr(strLine.substr(STATION+PAD2+LAT_EAST+LON_NORTH+STAT, STAT)));
				posUnc._semimMajor =  DoubleFromString<double>(trimstr(strLine.substr(STATION+PAD2+LAT_EAST+LON_NORTH+STAT+STAT, PREC)));
				posUnc._semimMinor =  DoubleFromString<double>(trimstr(strLine.substr(STATION+PAD2+LAT_EAST+LON_NORTH+STAT+STAT+PREC, PREC)));
				posUnc._orientation = DoubleFromString<double>(trimstr(strLine.substr(STATION+PAD2+LAT_EAST+LON_NORTH+STAT+STAT+PREC+PREC, PREC)));
				posUnc._xx =		  DoubleFromString<double>(trimstr(strLine.substr(STATION+PAD2+LAT_EAST+LON_NORTH+STAT+STAT+PREC+PREC+PREC, MSR)));
				posUnc._xy =		  DoubleFromString<double>(trimstr(strLine.substr(STATION+PAD2+LAT_EAST+LON_NORTH+STAT+STAT+PREC+PREC+PREC+MSR, MSR)));
				posUnc._xz =		  DoubleFromString<double>(trimstr(strLine.substr(STATION+PAD2+LAT_EAST+LON_NORTH+STAT+STAT+PREC+PREC+PREC+MSR+MSR, MSR)));

				// get yy from next row of variance matrix
				apu_file.getline(line, PRINT_LINE_LENGTH);
				strLine = std::string(line);
				posUnc._yy = DoubleFromString<double>(trimstr(strLine.substr(STATION+PAD2+LAT_EAST+LON_NORTH+STAT+STAT+PREC+PREC+PREC+MSR, MSR)));
				posUnc._yz = DoubleFromString<double>(trimstr(strLine.substr(STATION+PAD2+LAT_EAST+LON_NORTH+STAT+STAT+PREC+PREC+PREC+MSR+MSR, MSR)));

				// Get reference of this station in bstBinaryRecords.
				if (dataBlocks || vcv_units == XYZ_apu_ui)
				{
					_it_stnmap = equal_range(stnsMap_.begin(), stnsMap_.end(), posUnc._station, StationNameIDCompareName());

					// Not in the list? (Unlikely since adjust output will not contain stations not in the map, but test anyway)
					if (_it_stnmap.first == _it_stnmap.second)
						continue;
				}

				// Calculate vcv in local reference frame
				if (vcv_units == XYZ_apu_ui)
				{
					// get yy from next row of variance matrix
					apu_file.getline(line, PRINT_LINE_LENGTH);
					strLine = trimstr(std::string(line));
					posUnc._zz = DoubleFromString<double>(trimstr(strLine));

					vcv_cart.put(0, 0, posUnc._xx);
					vcv_cart.put(0, 1, posUnc._xy);
					vcv_cart.put(0, 2, posUnc._xz);
					vcv_cart.put(1, 1, posUnc._yy);
					vcv_cart.put(1, 2, posUnc._yz);
					vcv_cart.put(2, 2, posUnc._zz);

					vcv_cart.filllower();

					// Calculate correlations in local reference frame
					PropagateVariances_LocalCart<double>(vcv_cart, vcv_local, 
						Radians(bstBinaryRecords_.at(_it_stnmap.first->second).currentLatitude), 
						Radians(bstBinaryRecords_.at(_it_stnmap.first->second).currentLongitude), false);

					posUnc._xx = vcv_local.get(0, 0);
					posUnc._xy = vcv_local.get(0, 1);
					posUnc._xz = vcv_local.get(0, 2);
					posUnc._yy = vcv_local.get(1, 1);
					posUnc._yz = vcv_local.get(1, 2);
					posUnc._zz = vcv_local.get(2, 2);
				}
				else
					// ENU - don't need up component - ignore last (height) row 
					// of variance matrix
					apu_file.ignore(PRINT_LINE_LENGTH, '\n');		// ------------------------

				// Data in blocks?  Set the vector element accordingly, otherwise add to the vector
				if (dataBlocks)
					v_stn_pu_.at(_it_stnmap.first->second) = posUnc;
				else
					v_stn_pu_.push_back(posUnc);

				if (fabs(posUnc._semimMajor) > largest_uncertainty_)
					largest_uncertainty_ = posUnc._semimMajor;

				if (!fullCovarianceMatrix)
					continue;
				
				// skip over covariances
				for (stn_cov=stn+1; stn_cov<blockstnCount;  ++stn_cov)
				{
					// ignore covariance block
					apu_file.ignore(PRINT_LINE_LENGTH, '\n');		// ------------------------
					apu_file.ignore(PRINT_LINE_LENGTH, '\n');		// ------------------------
					apu_file.ignore(PRINT_LINE_LENGTH, '\n');		// ------------------------
				}
			}
		}
	}
	catch (const std::ios_base::failure& f) {
		if (!apu_file.eof())
		{
			std::stringstream ss;
			ss << "LoadPosUncertaintyFile(): An error was encountered when reading " << projectSettings_.o._apu_file << "." << std::endl << "  " << f.what();
			SignalExceptionPlot(ss.str(), 0, "i", &apu_file);
		}
	}

	apu_file.close();
}

void dna_plot::ComputeStationCorrection(it_vstn_t_const _it_stn, stationCorrections_t& stnCor,
	matrix_2d& currentEstimates, matrix_2d& initialEstimates)
{
	// Compute XYZ for current estimates
	GeoToCart<double>(
		Radians(_it_stn->currentLatitude),
		Radians(_it_stn->currentLongitude),
		_it_stn->currentHeight,
		currentEstimates.getelementref(0, 0),
		currentEstimates.getelementref(1, 0),
		currentEstimates.getelementref(2, 0),
		datum_.GetEllipsoidRef());

	// Compute XYZ for initial estimates
	GeoToCart<double>(
		Radians(_it_stn->initialLatitude),
		Radians(_it_stn->initialLongitude),
		_it_stn->initialHeight + _it_stn->geoidSep,
		initialEstimates.getelementref(0, 0),
		initialEstimates.getelementref(1, 0),
		initialEstimates.getelementref(2, 0),
		datum_.GetEllipsoidRef());

	// compute vertical angle
	stnCor._vAngle = VerticalAngle<double>(
		initialEstimates.get(0, 0),		// X1
		initialEstimates.get(1, 0),		// Y1
		initialEstimates.get(2, 0),		// Z1
		currentEstimates.get(0, 0),		// X2
		currentEstimates.get(1, 0),		// Y2
		currentEstimates.get(2, 0),		// Z2
		Radians(_it_stn->initialLatitude),
		Radians(_it_stn->initialLongitude),
		Radians(_it_stn->currentLatitude),
		Radians(_it_stn->currentLongitude),
		0., 0.,
		&stnCor._east, &stnCor._north, &stnCor._up);

	if (fabs(stnCor._east) < PRECISION_1E5)
		if (fabs(stnCor._north) < PRECISION_1E5)
			if (fabs(stnCor._up) < PRECISION_1E5)
				stnCor._vAngle = 0.;
	
	// compute horizontal correction
	stnCor._hDistance = magnitude(stnCor._east, stnCor._north);
			
	// calculate azimuth
	stnCor._azimuth = Direction<double>(
		initialEstimates.get(0, 0),		// X1
		initialEstimates.get(1, 0),		// Y1
		initialEstimates.get(2, 0),		// Z1
		currentEstimates.get(0, 0),		// X2
		currentEstimates.get(1, 0),		// Y2
		currentEstimates.get(2, 0),		// Z2
		Radians(_it_stn->currentLatitude),
		Radians(_it_stn->currentLongitude),
		&stnCor._east, &stnCor._north);

	if (fabs(stnCor._east) < PRECISION_1E5)
		if (fabs(stnCor._north) < PRECISION_1E5)
			stnCor._azimuth = 0.;

	// calculate distances
	stnCor._sDistance = magnitude<double>(
		initialEstimates.get(0, 0),		// X1
		initialEstimates.get(1, 0),		// Y1
		initialEstimates.get(2, 0),		// Z1
		currentEstimates.get(0, 0),		// X2
		currentEstimates.get(1, 0),		// Y2
		currentEstimates.get(2, 0));	// Z2
}
	

void dna_plot::ComputeStationCorrections()
{
	// reserve enough memory for each station
	v_stn_corrs_.reserve(bstBinaryRecords_.size());

	stationCorrections_t stnCor;
	matrix_2d currentEstimates(3, 1), initialEstimates(3, 1);

	average_correction_ = 0.0;
	UINT32 correction_count(0);

	// compute corrections for every station in the network
	for (it_vstn_t_const _it_stn=bstBinaryRecords_.begin();
		_it_stn!=bstBinaryRecords_.end(); 
		++_it_stn)
	{
		// Constrained?  No correction!
		if (_it_stn->stationConst[0] == 'C' && 
			_it_stn->stationConst[1] == 'C' && 
			_it_stn->stationConst[2] == 'C')
		{
			// push a zero correction
			v_stn_corrs_.push_back(stationCorrections(_it_stn->stationName));
			continue;
		}

		stnCor._station = _it_stn->stationName;
		ComputeStationCorrection(_it_stn, stnCor,
			currentEstimates, initialEstimates);

		average_correction_ += stnCor._hDistance;
		correction_count++;

		v_stn_corrs_.push_back(stnCor);
	}

	// determine average correction length
	average_correction_ /= correction_count;
}
	
void dna_plot::BuildParameterStationList()
{
	if (!plotBlocks_)
		return;

	v_parameterStationList_.resize(blockCount_);

	vUINT32 vStns;
	it_vUINT32_const _it_stn;

	for (UINT32 block(0); block!=blockCount_; ++block)
	{
		// Form combined stations list
		v_parameterStationList_.at(block) = v_ISL_.at(block);
		v_parameterStationList_.at(block).insert(v_parameterStationList_.at(block).end(),
			v_JSL_.at(block).begin(), v_JSL_.at(block).end());
		std::sort(v_parameterStationList_.at(block).begin(), v_parameterStationList_.at(block).end());

	}
}

// The station corrections file is produced by adjust.  Note that corrections for *all* stations may not 
// be printed, since the user has the option to restrict the output of station corrections using
// horizontal and vertical thresholds.
// Also, the coorctions file may contain a list of stations printed in blocks, or in a contiguous list.
// These options create complexities!
//
// The best way to deal with station corrections is to ignore the cor file and to recalculate corrections
// from the original coordinates in the bst file.
void dna_plot::LoadCorrectionsFile()
{
	// Load corrections file generated by adjust
	std::ifstream cor_file;
	try {
		// Load corrections file.  Throws runtime_error on failure.
		file_opener(cor_file, projectSettings_.o._cor_file, std::ios::in, ascii, true);
	}
	catch (const std::runtime_error& e) {
		SignalExceptionPlot(e.what(), 0, NULL);
	}

	UINT32 correction_count(0);
	UINT32 block(0), stn(0), blockstnCount, corFileBlockCount(blockCount_);
	std::string strLine;
	stationCorrections_t stnCor;

	average_correction_ = 0.0;

	char line[PRINT_LINE_LENGTH];
	bool dataBlocks(false);
	it_pair_string_vUINT32 _it_stnmap;

	try {
		
		// read header line
		cor_file.ignore(PRINT_LINE_LENGTH, '\n');		// ------------------------
		cor_file.ignore(PRINT_LINE_LENGTH, '\n');		//            DYNADJUST CORRECT...
		cor_file.ignore(PRINT_LINE_LENGTH, '\n');		// 
		cor_file.ignore(PRINT_LINE_LENGTH, '\n');		// Version
		cor_file.ignore(PRINT_LINE_LENGTH, '\n');		// Build
		cor_file.ignore(PRINT_LINE_LENGTH, '\n');		// File created
		cor_file.ignore(PRINT_LINE_LENGTH, '\n');		// File name
		cor_file.ignore(PRINT_LINE_LENGTH, '\n');		// 
		cor_file.getline(line, PRINT_LINE_LENGTH);		// Stations printed in blocks
		
		// Get the yes/no string and convert to bool
		strLine = trimstr(std::string(line));
		if (!iequals(strLine.substr(0, 16), "Stations printed"))
		{
			std::stringstream ss;
			// TODO - make use of Boost current function
			// http://www.boost.org/doc/libs/1_58_0/boost/current_function.hpp
			// and if required, print the filename and line number using __FILE__ and __LINE__
			// https://stackoverflow.com/questions/15305310/predefined-macros-for-function-name-func
			ss << "LoadCorrectionsFile(): " << projectSettings_.o._cor_file << " is corrupt." << std::endl;
			ss << "  Expected to find 'Stations printed in blocks' field." << std::endl;
			SignalExceptionPlot(ss.str(), 0, "i", &cor_file);
		}

		if ((dataBlocks = yesno_uint<bool, std::string>(strLine.substr(PRINT_VAR_PAD))))		// Data printed in blocks?
			v_stn_corrs_.resize(bstBinaryRecords_.size());
		else
		{
			v_stn_corrs_.reserve(bstBinaryRecords_.size());
			corFileBlockCount = 1;
		}

		cor_file.ignore(PRINT_LINE_LENGTH, '\n');		// 
		cor_file.ignore(PRINT_LINE_LENGTH, '\n');		// Corrections to Stations
		cor_file.ignore(PRINT_LINE_LENGTH, '\n');		// ------------------------
		
		for (block=0; block<corFileBlockCount; ++block)
		{
			if (dataBlocks)
			{
				blockstnCount = static_cast<UINT32>(v_ISL_.at(block).size() + v_JSL_.at(block).size());

				cor_file.ignore(PRINT_LINE_LENGTH, '\n');		// 
				cor_file.getline(line, PRINT_LINE_LENGTH);		// Block #

				strLine = trimstr(std::string(line));
				
				if (!iequals(strLine.substr(0, 5), "Block"))
				{
					std::stringstream ss;
					ss << "LoadCorrectionsFile(): " << projectSettings_.o._cor_file << " is corrupt." << std::endl;
					ss << "  Expected to read Block data, but found " << strLine << std::endl;
					SignalExceptionPlot(ss.str(), 0, "i", &cor_file);
				}

				if (LongFromString<UINT32>(strLine.substr(6)) != block + 1)
				{
					std::stringstream ss;
					ss << "LoadCorrectionsFile(): " << projectSettings_.o._cor_file << " is corrupt." << std::endl;
					ss << "  Expected to read data for block " << strLine.substr(6) << ", but found " << strLine << std::endl;
					SignalExceptionPlot(ss.str(), 0, "i", &cor_file);
				}
			}
			else
			{
				blockstnCount = stationCount_;
				cor_file.ignore(PRINT_LINE_LENGTH, '\n');		// 
			}
			
			cor_file.ignore(PRINT_LINE_LENGTH, '\n');		// Station
			cor_file.ignore(PRINT_LINE_LENGTH, '\n');		// ------------------------

			for (stn=0; stn<blockstnCount;  ++stn)
			{
				cor_file.getline(line, PRINT_LINE_LENGTH);		// Now the data...

				strLine = trimstr(std::string(line));
				if (strLine.length() < STATION)
				{
					if (cor_file.eof())
						break;				
					continue;
				}

				//// print...
				//// station and constraint
				//os << std::setw(STATION) << std::left << bstBinaryRecords_.at(stn).stationName << std::setw(PAD2) << " ";
				//// data
				//os << std::setw(MSR) << std::right << FormatDmsString(RadtoDms(azimuth), 4, true, false) << 
				//	std::setw(MSR) << std::right << FormatDmsString(RadtoDms(vertical_angle), 4, true, false) << 
				//	std::setw(MSR) << std::setprecision(PRECISION_MTR_STN) << std::fixed << std::right << slope_distance << 
				//	std::setw(MSR) << std::setprecision(PRECISION_MTR_STN) << std::fixed << std::right << horiz_distance << 
				//	std::setw(HEIGHT) << std::setprecision(PRECISION_MTR_STN) << std::fixed << std::right << local_12e << 
				//	std::setw(HEIGHT) << std::setprecision(PRECISION_MTR_STN) << std::fixed << std::right << local_12n << 
				//	std::setw(HEIGHT) << std::setprecision(PRECISION_MTR_STN) << std::fixed << std::right << local_12up << std::endl;

				stnCor._station   = trimstr(strLine.substr(0, STATION));
				stnCor._azimuth   = ParseDmsString<double>(trimstr(strLine.substr(STATION+PAD2, MSR)), " ");
				stnCor._vAngle	  = ParseDmsString<double>(trimstr(strLine.substr(STATION+PAD2+MSR, MSR)), " ");
				stnCor._sDistance = DoubleFromString<double>(trimstr(strLine.substr(STATION+PAD2+MSR+MSR, MSR)));
				stnCor._hDistance = DoubleFromString<double>(trimstr(strLine.substr(STATION+PAD2+MSR+MSR+MSR, MSR)));
				stnCor._east      = DoubleFromString<double>(trimstr(strLine.substr(STATION+PAD2+MSR+MSR+MSR+MSR, HEIGHT)));
				stnCor._north     = DoubleFromString<double>(trimstr(strLine.substr(STATION+PAD2+MSR+MSR+MSR+MSR+HEIGHT, HEIGHT)));
				stnCor._up        = DoubleFromString<double>(trimstr(strLine.substr(STATION+PAD2+MSR+MSR+MSR+MSR+HEIGHT+HEIGHT, HEIGHT)));


				// Get reference of this station in bstBinaryRecords.
				if (dataBlocks)
				{
					_it_stnmap = equal_range(stnsMap_.begin(), stnsMap_.end(), stnCor._station, StationNameIDCompareName());

					// Not in the list? (Unlikely since adjust output will not contain stations not in the map, but test anyway)
					if (_it_stnmap.first == _it_stnmap.second)
						continue;
		
					v_stn_corrs_.at(_it_stnmap.first->second) = stnCor;
				}
				else
					v_stn_corrs_.push_back(stnCor);

				average_correction_ += stnCor._hDistance;
				correction_count++;
			}
		}
	}
	catch (const std::ios_base::failure& f) {
		if (!cor_file.eof())
		{
			std::stringstream ss;
			ss << "LoadCorrectionsFile(): An error was encountered when reading " << projectSettings_.o._cor_file << "." << std::endl << "  " << f.what();
			SignalExceptionPlot(ss.str(), 0, "i", &cor_file);
		}
	}

	// determine average correction length
	average_correction_ /= correction_count;

	cor_file.close();
}
	

// Name:				SignalExceptionPlot
// Purpose:				Closes all files (if file pointers are passed in) and throws runtime_error
// Called by:			Any
// Calls:				runtime_error()
void dna_plot::SignalExceptionPlot(const std::string& msg, const int& line_no, const char *streamType, ...)
{
	plotStatus_ = PLOT_EXCEPTION_RAISED;

	if (streamType == NULL)
		throw NetPlotException(msg, line_no);

	std::ofstream* ofsDynaML;
	std::ifstream* ifsbinaryFile;

	va_list argptr; 
	va_start(argptr, streamType);

	while (*streamType != '\0')
	{
		//ifstream
		if (*streamType == 'i' )
		{
			ifsbinaryFile = va_arg(argptr, std::ifstream*);
			ifsbinaryFile->close();
		}
		//ofstream
		if (*streamType == 'o' )
		{
			ofsDynaML = va_arg(argptr, std::ofstream*);
			ofsDynaML->close();
		}
		streamType++;
	}

	va_end(argptr);

	throw NetPlotException(msg, line_no);
}


//void dna_plot::PrintGMTPlotCoords(vdnaStnPtr* vStations, vdnaMsrPtr* vMeasurements, ostream& osStn, ostream& osStn2, ostream& osMsr)
//{
//	_it_vdnastnptr _it_stn(vStations->begin());
//	vdnaMsrPtr::iterator _it_msr(vMeasurements->begin());
//	
//	// Print station coords
//	for (; _it_stn!=vStations->end(); _it_stn++)
//		_it_stn->get()->coutStationData(osStn, osStn2, GMT_OUT);
//
//	//v_string_uint32_pair::iterator _it_stnmap(stnsMap_.begin());
//
//	UINT32 precision;
//	_COORD_TYPE_ ctType;
//
//	std::vector<CDnaDirection>* vdirns;
//	std::vector<CDnaGpsBaseline>* vgpsBsls;
//
//	for (; _it_msr!=vMeasurements->end(); _it_msr++)
//	{
//		switch (_it_msr->get()->GetTypeC())
//		{
//		case 'A':	// Horizontal angles
//			// Print 3 - 1 - 2
//			//
//			// Third
//			it_stnmap_range = equal_range(stnsMap_.begin(), stnsMap_.end(), 
//				_it_msr->get()->GetTarget2(), StationNameIDCompareName());
//
//			if (it_stnmap_range.first == it_stnmap_range.second)
//			{
//				std::cout << _it_msr->get()->GetTarget2() << " is not in the list of network stations.";
//				continue;
//			}
//			precision = 3;
//			if ((ctType = vStations->at(it_stnmap_range.first->second)->GetMyCoordTypeC()) == LLH_type_i)
//				precision = 10;
//
//			osMsr << std::setprecision(precision) << std::fixed << 
//				(ctType == LLH_type_i ? 
//					Degrees(vStations->at(it_stnmap_range.first->second)->GetYAxis()) :
//					vStations->at(it_stnmap_range.first->second)->GetYAxis());
//			osMsr << "  ";
//			osMsr << std::setprecision(precision) << std::fixed << 
//				(ctType == LLH_type_i ? 
//					Degrees(vStations->at(it_stnmap_range.first->second)->GetXAxis()) :
//					vStations->at(it_stnmap_range.first->second)->GetXAxis());
//			osMsr << std::endl;
//		case 'B': // Geodetic azimuth
//		case 'K': // Astronomic azimuth
//		case 'C': // Chord dist
//		case 'E': // Ellipsoid arc
//		case 'M': // MSL arc
//		case 'S': // Slope distance
//		case 'L': // Level difference
//		case 'V': // Zenith distance
//		case 'Z': // Vertical angle
//			// First
//			it_stnmap_range = equal_range(stnsMap_.begin(), stnsMap_.end(), 
//				_it_msr->get()->GetFirst(), StationNameIDCompareName());
//
//			if (it_stnmap_range.first == it_stnmap_range.second)
//			{
//				std::cout << _it_msr->get()->GetFirst() << " is not in the list of network stations.";
//				continue;
//			}
//			precision = 3;
//			if ((ctType = vStations->at(it_stnmap_range.first->second)->GetMyCoordTypeC()) == LLH_type_i)
//				precision = 10;
//
//			osMsr << std::setprecision(precision) << std::fixed << (ctType == LLH_type_i ? Degrees(vStations->at(it_stnmap_range.first->second)->GetYAxis()): vStations->at(it_stnmap_range.first->second)->GetYAxis());
//			osMsr << "  ";
//			osMsr << std::setprecision(precision) << std::fixed << (ctType == LLH_type_i ? Degrees(vStations->at(it_stnmap_range.first->second)->GetXAxis()): vStations->at(it_stnmap_range.first->second)->GetXAxis());
//			osMsr << std::endl;
//			
//			// Second
//			it_stnmap_range = equal_range(stnsMap_.begin(), stnsMap_.end(), 
//				_it_msr->get()->GetTarget(), StationNameIDCompareName());
//
//			if (it_stnmap_range.first == it_stnmap_range.second)
//			{
//				std::cout << _it_msr->get()->GetTarget() << " is not in the list of network stations.";
//				continue;
//			}
//			precision = 3;
//			if ((ctType = vStations->at(it_stnmap_range.first->second)->GetMyCoordTypeC()) == LLH_type_i)
//				precision = 10;
//
//			osMsr << std::setprecision(precision) << std::fixed << (ctType == LLH_type_i ? Degrees(vStations->at(it_stnmap_range.first->second)->GetYAxis()): vStations->at(it_stnmap_range.first->second)->GetYAxis());
//			osMsr << "  ";
//			osMsr << std::setprecision(precision) << std::fixed << (ctType == LLH_type_i ? Degrees(vStations->at(it_stnmap_range.first->second)->GetXAxis()): vStations->at(it_stnmap_range.first->second)->GetXAxis());
//			osMsr << std::endl << "#" << std::endl;
//			break;
//		case 'D': // Direction set
//			vdirns = _it_msr->get()->GetDirections_ptr();
//			PrintGMTPlotCoords_D(vStations, vdirns, osMsr);
//			break;
//		case 'G': // GPS Baseline
//		case 'X': // GPS Baseline cluster
//			vgpsBsls = _it_msr->get()->GetBaselines_ptr();
//			PrintGMTPlotCoords_GX(vStations, vgpsBsls, osMsr);
//			break;
//		}
//	}
//}
//
//
//void dna_plot::PrintGMTPlotCoords_D(vdnaStnPtr* vStations, std::vector<CDnaDirection>* vDirections, ostream& osMsr)
//{
//	std::vector<CDnaDirection>::iterator _it_dirn(vDirections->begin());
//	//v_string_uint32_pair::iterator _it_stnmap(stnsMap_.begin());
//
//	UINT32 precision;
//	_COORD_TYPE_ ctType;
//
//	for (; _it_dirn!=vDirections->end(); _it_dirn++)
//	{
//		// First
//		it_stnmap_range = equal_range(stnsMap_.begin(), stnsMap_.end(), 
//			_it_dirn->GetFirst(), StationNameIDCompareName());
//
//		if (it_stnmap_range.first == it_stnmap_range.second)
//		{
//			std::cout << _it_dirn->GetFirst() << " is not in the list of network stations.";
//			continue;
//		}
//		precision = 3;
//		if ((ctType = vStations->at(it_stnmap_range.first->second)->GetMyCoordTypeC()) == LLH_type_i)
//			precision = 10;
//
//		osMsr << std::setprecision(precision) << std::fixed << (ctType == LLH_type_i ? Degrees(vStations->at(it_stnmap_range.first->second)->GetYAxis()): vStations->at(it_stnmap_range.first->second)->GetYAxis());
//		osMsr << "  ";
//		osMsr << std::setprecision(precision) << std::fixed << (ctType == LLH_type_i ? Degrees(vStations->at(it_stnmap_range.first->second)->GetXAxis()): vStations->at(it_stnmap_range.first->second)->GetXAxis());
//		osMsr << std::endl;
//
//		// Second
//		it_stnmap_range = equal_range(stnsMap_.begin(), stnsMap_.end(), 
//			_it_dirn->GetTarget(), StationNameIDCompareName());
//
//		if (it_stnmap_range.first == it_stnmap_range.second)
//		{
//			std::cout << _it_dirn->GetTarget() << " is not in the list of network stations.";
//			continue;
//		}
//		precision = 3;
//		if ((ctType = vStations->at(it_stnmap_range.first->second)->GetMyCoordTypeC()) == LLH_type_i)
//			precision = 10;
//
//		osMsr << std::setprecision(precision) << std::fixed << (ctType == LLH_type_i ? Degrees(vStations->at(it_stnmap_range.first->second)->GetYAxis()): vStations->at(it_stnmap_range.first->second)->GetYAxis());
//		osMsr << "  ";
//		osMsr << std::setprecision(precision) << std::fixed << (ctType == LLH_type_i ? Degrees(vStations->at(it_stnmap_range.first->second)->GetXAxis()): vStations->at(it_stnmap_range.first->second)->GetXAxis());
//		osMsr << std::endl << "#" << std::endl;
//	}
//}
//
//
//void dna_plot::PrintGMTPlotCoords_GX(vdnaStnPtr* vStations, std::vector<CDnaGpsBaseline>* vGPSBaselines, ostream& osMsr)
//{
//	std::vector<CDnaGpsBaseline>::iterator _it_bsl(vGPSBaselines->begin());
//	//v_string_uint32_pair::iterator _it_stnmap(stnsMap_.begin());
//
//	UINT32 precision;
//	_COORD_TYPE_ ctType;
//
//	for (; _it_bsl!=vGPSBaselines->end(); _it_bsl++)
//	{
//		// First
//		it_stnmap_range = equal_range(stnsMap_.begin(), stnsMap_.end(), 
//			_it_bsl->GetFirst(), StationNameIDCompareName());
//
//		if (it_stnmap_range.first == it_stnmap_range.second)
//		{
//			std::cout << _it_bsl->GetFirst() << " is not in the list of network stations.";
//			continue;
//		}
//		precision = 3;
//		if ((ctType = vStations->at(it_stnmap_range.first->second)->GetMyCoordTypeC()) == LLH_type_i)
//			precision = 10;
//
//		osMsr << std::setprecision(precision) << std::fixed << (ctType == LLH_type_i ? Degrees(vStations->at(it_stnmap_range.first->second)->GetYAxis()): vStations->at(it_stnmap_range.first->second)->GetYAxis());
//		osMsr << "  ";
//		osMsr << std::setprecision(precision) << std::fixed << (ctType == LLH_type_i ? Degrees(vStations->at(it_stnmap_range.first->second)->GetXAxis()): vStations->at(it_stnmap_range.first->second)->GetXAxis());
//		osMsr << std::endl;
//
//		// Second
//		it_stnmap_range = equal_range(stnsMap_.begin(), stnsMap_.end(), 
//			_it_bsl->GetTarget(), StationNameIDCompareName());
//
//		if (it_stnmap_range.first == it_stnmap_range.second)
//		{
//			std::cout << _it_bsl->GetTarget() << " is not in the list of network stations.";
//			continue;
//		}
//		precision = 3;
//		if ((ctType = vStations->at(it_stnmap_range.first->second)->GetMyCoordTypeC()) == LLH_type_i)
//			precision = 10;
//
//		osMsr << std::setprecision(precision) << std::fixed << (ctType == LLH_type_i ? Degrees(vStations->at(it_stnmap_range.first->second)->GetYAxis()): vStations->at(it_stnmap_range.first->second)->GetYAxis());
//		osMsr << "  ";
//		osMsr << std::setprecision(precision) << std::fixed << (ctType == LLH_type_i ? Degrees(vStations->at(it_stnmap_range.first->second)->GetXAxis()): vStations->at(it_stnmap_range.first->second)->GetXAxis());
//		osMsr << std::endl << "#" << std::endl;
//	}
//}



void dna_plot::ReduceStationCoordinates(station_t* stationRecord)
{
	stationRecord->initialLatitude = Degrees(stationRecord->initialLatitude);
	stationRecord->initialLongitude = Degrees(stationRecord->initialLongitude);
    stationRecord->currentLatitude = Degrees(stationRecord->currentLatitude);
    stationRecord->currentLongitude = Degrees(stationRecord->currentLongitude);
}

void dna_plot::NormaliseScale(double& scale)
{
	if (scale > 5000000. && scale < 10000000.)
		scale = 5000000.;
	else if (scale > 1250000.)
		scale = 1500000.;
	else if (scale > 1000000.)
		scale = 1250000.;
	else if (scale > 750000.)
		scale = 1000000.;
	else if (scale > 500000.)
		scale = 750000.;
	else if (scale > 350000.)
		scale = 500000.;
	else if (scale > 250000.)
		scale = 350000.;
	else if (scale > 200000.)
		scale = 250000.;
	else if (scale > 150000.)
		scale = 200000.;
	else if (scale > 100000.)
		scale = 150000.;
	else if (scale > 75000.)
		scale = 100000.;
	else if (scale > 50000.)
		scale = 75000.;
	else if (scale > 35000.)
		scale = 50000.;
	else if (scale > 25000.)
		scale = 35000.;
	else if (scale > 10000.)
		scale = 25000.;
	else if (scale > 7500.)
		scale = 10000.;
	else if (scale > 5000.)
		scale = 7500.;
	else if (scale > 2500.)
		scale = 5000.;
	else if (scale > 1000.)
		scale = 2500.;
	else if (scale > 500.)
		scale = 1000.;
	else
		scale = 500.;
}

void dna_plot::NormaliseScaleBar(double& scale_bar_width)
{
	if (scale_bar_width > 1.)
		scale_bar_width = floor(scale_bar_width);

	if (scale_bar_width > 10000)
		scale_bar_width -= (int)scale_bar_width % 10000;
	else if (scale_bar_width > 1000)
		scale_bar_width -= (int)scale_bar_width % 1000;
	else if (scale_bar_width > 100)
		scale_bar_width -= (int)scale_bar_width % 100;
	else if (scale_bar_width > 10)
		scale_bar_width -= (int)scale_bar_width % 10;
	else if (scale_bar_width > 5)
		scale_bar_width -= (int)scale_bar_width % 5;
	else if (scale_bar_width > 1)
		scale_bar_width -= (int)scale_bar_width % 1;
	else if (scale_bar_width > 0.5)
		scale_bar_width = 0.5;
	else if (scale_bar_width > 0.25)
		scale_bar_width = 0.25;
	
	else if (scale_bar_width > 0.1)
		scale_bar_width = 0.1;
	else if (scale_bar_width > 0.075)
		scale_bar_width = 0.075;
	else if (scale_bar_width > 0.05)
		scale_bar_width = 0.05;
	else if (scale_bar_width > 0.025)
		scale_bar_width = 0.025;
	
	else if (scale_bar_width > 0.01)
		scale_bar_width = 0.01;
	else if (scale_bar_width > 0.0075)
		scale_bar_width = 0.0075;
	else if (scale_bar_width > 0.005)
		scale_bar_width = 0.005;
	else if (scale_bar_width > 0.0025)
		scale_bar_width = 0.0025;
	
	else if (scale_bar_width > 0.0001)
		scale_bar_width = 0.0001;
	else if (scale_bar_width > 0.00075)
		scale_bar_width = 0.00075;
	else if (scale_bar_width > 0.0005)
		scale_bar_width = 0.0005;
	else if (scale_bar_width > 0.00025)
		scale_bar_width = 0.00025;
	
}


void dna_plot::NormaliseGraticule(double& graticule_width_, UINT32& graticule_width_precision)
{
	double graticule = graticule_width_ * 3600.0;	// convert to seconds
	if (graticule > 1.)
		graticule = floor(graticule);
	else
		graticule_width_precision = 16;

	if (graticule > 7200)			// 2 degrees
		graticule -= (int)graticule % 7200;
	else if (graticule > 3600)		// 1 degrees
		graticule -= (int)graticule % 3600;
	else if (graticule > 1800)		// 30 minutes
		graticule -= (int)graticule % 1800;
	else if (graticule > 900)		// 15 minutes
		graticule -= (int)graticule % 900;
	else if (graticule > 600)		// 10 minutes
		graticule -= (int)graticule % 600;
	else if (graticule > 300)		// 5 minutes
		graticule -= (int)graticule % 300;
	else if (graticule > 60)		// 1 minute
		graticule -= (int)graticule % 60;
	else if (graticule > 30)		// 30 seconds
		graticule -= (int)graticule % 30;
	else if (graticule > 15)		// 15 seconds
		graticule -= (int)graticule % 15;
	else if (graticule > 10)		// 10 seconds
		graticule -= (int)graticule % 10;
	else if (graticule > 5)			// 5 seconds
		graticule -= (int)graticule % 5;
	else if (graticule > 1)			// 5 seconds
		graticule -= (int)graticule % 1;
	else if (graticule > 0.5)		// .5 seconds
		graticule = 0.5;
	else if (graticule > 0.25)		// .25 seconds
		graticule = 0.25;


	// convert to degrees
	graticule_width_ = graticule / 3600.;
}

void dna_plot::SelectCoastlineResolution(const double& dDimension, std::string& coastResolution, plot_settings* plotCriteria)
{
	plotCriteria->_coasline_resolution = low;
	coastResolution = "l";
	
	if (dDimension < 900000)
	{
		plotCriteria->_coasline_resolution = full;
		coastResolution = "f";
		return;
	}
	if (dDimension < 2000000)
	{
		coastResolution = "h";
		plotCriteria->_coasline_resolution = high;
		return;
	}
	if (dDimension < 9000000)
	{
		plotCriteria->_coasline_resolution = intermediate;
		coastResolution = "i";
	}	
}

}	// namespace mathcomp
}	// namespace dynadjust
