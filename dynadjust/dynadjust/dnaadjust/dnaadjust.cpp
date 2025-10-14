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
// Description  : DynAdjust Network Adjustment library
//============================================================================

#include <dynadjust/dnaadjust/dnaadjust.hpp>
#include <dynadjust/dnaadjust/dnaadjust-multi.cpp>

namespace dynadjust {
namespace networkadjust {

std::mutex combine_blockMutex;
std::mutex current_blockMutex, current_iterationMutex, maxCorrMutex;
std::mutex adj_file_mutex, xyz_file_mutex, dbg_file_mutex;

// multi thread adjustment variables
concurrent_queue<UINT32> combineAdjustmentQueue;
concurrent_queue<UINT32> prepareAdjustmentQueue;
std::exception_ptr fwd_error, rev_error, cmb_error, prep_error;

dna_adjust::dna_adjust()
	: isPreparing_(false)
	, isAdjusting_(false)
	, isCombining_(false)
	, forward_(true)
	, isFirstTimeAdjustment_(true)
	, isIterationComplete_(false)
	, isAdjustmentQuestionable_(false)
	, blockCount_(1)
	, currentBlock_(0)
	, total_time_(0)
	, adjustStatus_(ADJUST_SUCCESS)
	, currentIteration_(0)
	, datum_(DEFAULT_EPSG_U)
	, bmsr_count_(0)
	, bstn_count_(0)
	, asl_count_(0)
	, measurementParams_(0)
	, measurementCount_(0)
	, unknownParams_(0)
	, unknownsCount_(0)
	, chiSquared_(0.)
	, sigmaZero_(0.)
	, sigmaZeroSqRt_(0.)
	, chiSquaredUpperLimit_(0.)
	, chiSquaredLowerLimit_(0.)
	, globalPelzerReliability_(0.)
	, degreesofFreedom_(0)
	, passFail_(test_stat_pass)
	, maxCorr_(0.)
	, criticalValue_(1.68)
	, allStationsFixed_(false)
	, databaseIDsLoaded_(false)
	, isCancelled_(false)
{
	statusMessages_.clear();
	bstBinaryRecords_.clear();
	bmsBinaryRecords_.clear();
	vAssocMsrList_.clear();

	debug_file.clear();
	v_pseudoMeasCountFwd_.clear();
	v_measurementParams_.clear();
	v_measurementCount_.clear();
	v_measurementVarianceCount_.clear();
	v_unknownParams_.clear();
	v_unknownsCount_.clear();
	v_sigmaZero_.clear();
	v_chiSquaredUpperLimit_.clear();
	v_chiSquaredLowerLimit_.clear();
	v_passFail_.clear();
	
	v_originalStations_.clear();
	v_design_.clear();
	v_measMinusComp_.clear();

	v_AtVinv_.clear();
	v_normals_.clear();
	v_estimatedStations_.clear();
	v_rigorousStations_.clear();
	v_junctionVariances_.clear();
	v_junctionVariancesFwd_.clear();
	v_junctionEstimatesFwd_.clear();
	v_junctionEstimatesRev_.clear();
	v_rigorousVariances_.clear();
	v_precAdjMsrsFull_.clear();
	v_corrections_.clear();
	v_blockStationsMap_.clear();

	v_parameterStationCount_.clear();
	v_parameterStationList_.clear();

	// multi thread (for reverse and combine passes)
	v_designR_.clear();
	v_measMinusCompR_.clear();
	v_AtVinvR_.clear();
	v_normalsR_.clear();
	v_estimatedStationsR_.clear();
	v_junctionVariancesR_.clear();
	v_normalsRC_.clear();

#ifdef _MSC_VER
#if (_MSC_VER < 1900)
	{
		// this function is obsolete in MS VC++ 14.0, VS2015
		// Set scientific format to print two places for the exponent
		_set_output_format(_TWO_DIGIT_EXPONENT);
	}
#endif
#endif

	// v_correctionsR_ is only used in phased adjustment.
	// In single thread mode, it holds rigorous corrections 
	// for last block.  In multi-thread mode, it is used
	// for in the reverse thread.
	switch (projectSettings_.a.adjust_mode)
	{
	case PhasedMode:
		v_correctionsR_.clear();
	}

	// Initialize the printer
	printer_ = std::make_unique<DynAdjustPrinter>(*this);
}
	

dna_adjust::~dna_adjust()
{
	if (adjustStatus_ == ADJUST_EXCEPTION_RAISED)
	{
		try {
			debug_file.close();
		}
		catch (...) { }

		try {
			adj_file.close();
		}
		catch (...) { }

		try {
			xyz_file.close();
		}
		catch (...) { }
	}
}
	

UINT32 dna_adjust::CurrentIteration() const 
{ 
	std::lock_guard<std::mutex> lock(current_iterationMutex);
	return currentIteration_; 
}

UINT32& dna_adjust::incrementIteration() 
{ 
	std::lock_guard<std::mutex> lock(current_iterationMutex);
	return ++currentIteration_; 
}

void dna_adjust::initialiseIteration(const UINT32& iteration) 
{ 
	std::lock_guard<std::mutex> lock(current_iterationMutex);
	currentIteration_ = iteration; 
}

void dna_adjust::InitialiseAdjustment()
{
	adj_file << std::endl << "+ Initialising adjustment" << std::endl;
	
	// Confidence interval
	double conf(projectSettings_.a.confidence_interval * 0.01);
	boost::math::normal dist(0.0, 1.0);
	conf += (1.0 - conf) / 2.0;
	criticalValue_ = boost::math::quantile(dist, conf);

	potentialOutlierCount_ = 0;

	blockCount_ = 1;
	v_blockStationsMap_.resize(blockCount_);

	adjustStatus_ = ADJUST_SUCCESS;
	statusMessages_.clear();
	currentBlock_ = 0;
	initialiseIteration();
	
	bstBinaryRecords_.clear();
	bmsBinaryRecords_.clear();

	// output precision helpers
	PRECISION_SEC_STN = projectSettings_.o._precision_seconds_stn;
	PRECISION_MTR_STN = projectSettings_.o._precision_metres_stn;
	PRECISION_SEC_MSR = projectSettings_.o._precision_seconds_msr;
	PRECISION_MTR_MSR = projectSettings_.o._precision_metres_msr;

	// constraint values
	// cartesian reference frame
	_var_C = projectSettings_.a.fixed_std_dev * projectSettings_.a.fixed_std_dev;		// metres
	_var_F = projectSettings_.a.free_std_dev * projectSettings_.a.free_std_dev;			// metres
	
	// form inverse of variance for fixed stations
	_inv_var_cart_c.redim(3, 3);
	_inv_var_cart_c.put(0, 0, 1./_var_C);
	_inv_var_cart_c.put(1, 1, 1./_var_C);
	_inv_var_cart_c.put(2, 2, 1./_var_C);

	// form inverse of variance for free stations
	_inv_var_cart_f.redim(3, 3);
	_inv_var_cart_f.put(0, 0, 1./_var_F);
	_inv_var_cart_f.put(1, 1, 1./_var_F);
	_inv_var_cart_f.put(2, 2, 1./_var_F);

}

void dna_adjust::CreateMeasurementTally(const UINT32& block)
{
	v_msrTally_.at(block).CreateTally(bmsBinaryRecords_, v_CML_.at(block));
}
	


	

void dna_adjust::PrepareAdjustment(const project_settings& projectSettings)
{
	isPreparing_ = true;
	isAdjusting_ = true;
	isCombining_ = false;
	isFirstTimeAdjustment_ = true;
	projectSettings_ = projectSettings;

	if (projectSettings_.a.stage && 
		(projectSettings_.a.adjust_mode == SimultaneousMode || 
			projectSettings_.a.multi_thread))
		projectSettings_.a.stage = false;

	// Load the bst/bms meta and set the default 
	// reference frame (via binary station file)
	SetDefaultReferenceFrame();

	// Open output files for printing adjustment results
	OpenOutputFileStreams();
	printer_->PrintOutputFileHeaderInfo();
	
	// Initialise adjustment consts
	InitialiseAdjustment();
	
	// Load network files
	LoadNetworkFiles();

	// Load type b uncertainties, method handler, and the station map
	InitialiseTypeBUncertainties();

	isFirstTimeAdjustment_ = !bms_meta_.reduced;

	// Resizes matrix vectors using blockCount_.  
	// For phased adjustments, segmentation file is loaded
	ResizeMatrixVectors();

	if (!projectSettings_.a.report_mode)
	{
		std::string block_str(" block");
		if (blockCount_ > 1)
			block_str.append("s");
		block_str.append(")... ");

		adj_file << std::endl << "+ Preparing for adjustment";
		switch (projectSettings_.a.adjust_mode)
		{
		case Phased_Block_1Mode:
		case PhasedMode:
			adj_file << std::left << " (" << blockCount_<< block_str;
			if (projectSettings_.a.multi_thread && projectSettings_.g.verbose > 2)
				adj_file << std::endl;
			break;
		case SimultaneousMode:
			adj_file << std::left << "... ";
		}
		adj_file.flush();
	}
	
	try {
		
		// PrepareAdjustmentMultiThread works well but fails if there 
		// is insufficient memory.  Needs a little more work to catch
		// out of memory exception.

#ifndef _MSDEBUG
		// Only call PrepareAdjustmentMultiThread if:
		// - not in staged adjustment mode, and
		// - not running debugger
		if (!projectSettings_.a.stage)
			PrepareAdjustmentMultiThread();
		else
#endif // #ifndef _MSDEBUG
		{
			// Prepare initial matrices for least squares adjustment
			for (UINT32 block(0); block<blockCount_; ++block)
			{
				// update block no. for progress count
				SetcurrentBlock(block);

				CreateMeasurementTally(block);

				// Prepare adjustment for this block
				PrepareAdjustmentBlock(block);
			}
		}
	}
	catch (const NetMemoryException& e) {
		std::stringstream ss;
		ss << "PrepareAdjustment(): Process terminated while allocating memory for the " << std::endl << "  adjustment matrices. " <<
			"Details: " << std::endl << "  " << e.what() << std::endl;
		adj_file << "- Error" << std::endl << "  " << ss.str();
		SignalExceptionAdjustment(ss.str(), currentBlock_);
	}
	catch (const std::runtime_error& e) {
		std::stringstream ss;
		ss << "PrepareAdjustment(): Process terminated while preparing the " << std::endl << "  adjustment matrices. " <<
			"Details: " << std::endl << "  " << e.what() << std::endl;
		adj_file << "- Error" << std::endl << "  " << ss.str();
		SignalExceptionAdjustment(ss.str(), currentBlock_);
	}
	catch (...) {
		std::stringstream ss;
		ss << "PrepareAdjustment(): Process terminated while preparing the " << std::endl << "  adjustment matrices. " << std::endl;
		adj_file << "- Error" << std::endl << "  " << ss.str();
		SignalExceptionAdjustment(ss.str(), 0);
	}

	bms_meta_.reduced = true;

	if (projectSettings_.a.stage)
	{
		// Close binary file streams
		if (projectSettings_.a.recreate_stage_files)
			CloseStageFileStreams();

		// Set memory mapped file regions
		// NOTE - previously created files must exist and match the 
		// dimensions of the current matrix sizes			
		try {
			SetMapRegions();
		}
		catch (boost::interprocess::interprocess_exception& e){
			std::stringstream ss;
			ss << "PrepareMappedRegions() terminated while creating memory map" << std::endl;
			ss << "  regions from .mtx stage files. Details:\n  " << e.what() << std::endl << std::endl;
			ss << "  Please ensure the .mtx stage files from a previous staged" << std::endl;
			ss << "  adjustment exist, or re-run the adjustment using the" << std::endl;
			ss << "  --" << RECREATE_STAGE_FILES << " option." << std::endl;
			adj_file << std::endl << "- Error: " << ss.str() << std::endl;
			SignalExceptionAdjustment(ss.str(), 0);
		}
	}

	degreesofFreedom_ = measurementParams_ - unknownParams_;

	isPreparing_ = false;

	if (projectSettings_.a.multi_thread)
	{
		try {
			v_designR_ = v_design_;
			v_measMinusCompR_ = v_measMinusComp_;
			v_AtVinvR_ = v_AtVinv_;
			v_estimatedStationsR_ = v_estimatedStations_;
			v_junctionVariancesR_ = v_junctionVariances_;
		}
		catch (...)
		{
			std::stringstream ss;
			ss << "Failed to assign matrices for multi-thread adjustment." << std::endl;
			SignalExceptionAdjustment(ss.str(), 0);
		}
	}

	if (projectSettings_.a.adjust_mode == SimultaneousMode)
		BuildSimultaneousStnAppearance();

	BuildUniqueBlockStationMap();

	if (projectSettings_.o._adj_msr_final || 
		projectSettings_.o._adj_msr_iteration)
		// Create a vector of pairs (cml_id , block)
		// FormUniqueMsrList returns msr_block sorted according to block
		FormUniqueMsrList();
	
	if (projectSettings_.g.verbose > 3)
		debug_file << std::endl << std::endl;

	if (!projectSettings_.a.report_mode)
	{
		if (projectSettings_.a.multi_thread && projectSettings_.g.verbose > 2)
			adj_file << "+ Preparing for adjustment is";

		adj_file <<  " done.";

	}
}


void dna_adjust::UpdateBinaryFiles()
{
	snprintf(bst_meta_.modifiedBy, sizeof(bst_meta_.modifiedBy), "%s", __BINARY_NAME__);
	bst_meta_.reduced = true;

	try {
		// Write binary stations data.  Throws runtime_error on failure.
		BstFile bst;
		bst.WriteFile(projectSettings_.a.bst_file, &bstBinaryRecords_, bst_meta_);
	}
	catch (const std::runtime_error& e) {
		SignalExceptionAdjustment(e.what(), 0);
	}

	snprintf(bms_meta_.modifiedBy, sizeof(bms_meta_.modifiedBy), "%s", __BINARY_NAME__);
	bms_meta_.reduced = true;

	try {
		// Write binary measurements data.  Throws runtime_error on failure.
		BmsFile bms;
		bms.WriteFile(projectSettings_.a.bms_file, &bmsBinaryRecords_, bms_meta_);
	}
	catch (const std::runtime_error& e) {
		SignalExceptionAdjustment(e.what(), 0);
	}	
}
	

void dna_adjust::UpdateAdjustment(bool iterate)
{
	isPreparing_ = true;

	UINT32 block;

	bool lastBlock, updateGeographicCoordinates;
		
	try {
	
		// Prepare initial matrices for least squares adjustment
		for (block=0; block<blockCount_; ++block)
		{
			if (IsCancelled())
				break;

			switch (projectSettings_.a.adjust_mode)
			{
			case PhasedMode:
			case Phased_Block_1Mode:

				lastBlock = v_blockMeta_.at(block)._blockLast;
				updateGeographicCoordinates = 
					v_msrTally_.at(block).ContainsNonGPS() ||		// if this block contains measurements in the local reference frame, or
					lastBlock ||									// if this is the last block, or
					!iterate;										// if no further iterations are required

				if (updateGeographicCoordinates || lastBlock)
					// For staged adjustments, load block info
					if (projectSettings_.a.stage)
						DeserialiseBlockFromMappedFile(block, 1, sf_estimated_stns);

				// The last block in the forward pass is rigorous.  However, the coordinates were not updated 
				// since both phased and multi thread adjustments require original coordinates for the reverse
				// pass, both of which take place after the forward pass has completed.
				if (lastBlock)
				{
					// For staged adjustments, load block info
					if (projectSettings_.a.stage)
						DeserialiseBlockFromMappedFile(block, 4,
							sf_normals, sf_original_stns, 
							sf_rigorous_vars, sf_rigorous_stns);

					v_estimatedStations_.at(block) = v_rigorousStations_.at(block);
					v_originalStations_.at(block) = v_rigorousStations_.at(block);
					v_normals_.at(block) = v_rigorousVariances_.at(block);

					if (projectSettings_.a.stage)
						SerialiseBlockToMappedFile(block, 2,
							sf_estimated_stns, sf_original_stns);
				}

				// Update geographic coordinates if:
				//	- The network contains non-GPS measurements, so that partial
				//	  derivatives can be correctly formed
				//	- This is the last block in a contiguous network
				//	- There are no more iterations required
				if (updateGeographicCoordinates)
					UpdateGeographicCoordsPhased(block, &v_estimatedStations_.at(block));

				if (projectSettings_.a.stage)
				{
					// For staged adjustments, unload all matrix data
					UnloadBlock(block);
					continue;
				}

				break;
			case SimultaneousMode:
				// Update geographic coordinates
				if (v_msrTally_.at(block).ContainsNonGPS() ||		// if this block contains measurements in the local reference frame, or
					!iterate)										// if no further iterations are required
					UpdateGeographicCoords();
				break;
			}

			// Update measurements-computed vector using new estimates
			FillDesignNormalMeasurementsMatrices(false, block, false);

			// If no further iterations are required, then don't update the normals.
			// This is because the inverse of the normals is needed for computing
			// precision of adjusted measurements, whereas UpdateNormals() only
			// produces normals prior to inversion.
			// Hence, only design and msr-comp are required so as to compute stats
			if (!iterate)
				continue;
			
			switch (projectSettings_.a.adjust_mode)
			{
			case PhasedMode:
			case Phased_Block_1Mode:
				v_normals_.at(block).zero();
				UpdateNormals(block, false);
				
				if (projectSettings_.a.multi_thread)
				{
					v_estimatedStationsR_.at(block) = v_rigorousStations_.at(block);

					// Update measurements-computed vector for reverse thread using new estimates
					FillDesignNormalMeasurementsMatrices(false, block, true);
				}

				// Back up normals.  This copy contains the contributions from all
				// apriori measurement variances, excluding parameter station 
				// variances and junction station variances
				v_normalsR_.at(block) = v_normals_.at(block);				
				AddConstraintStationstoNormalsForward(block);

				break;
			case SimultaneousMode:
				if (v_msrTally_.at(0).ContainsNonGPS())
				{
					// update normals
					v_normals_.at(0).zero();
					UpdateNormals(0, false);
					AddConstraintStationstoNormalsSimultaneous(0);
				}
				break;
			}			
		}
	}
	catch (const NetMemoryException& e) {
		std::stringstream ss;
		ss << "UpdateAdjustment(): Process terminated while allocating memory for the " << std::endl << "  adjustment matrices. " <<
			"Details: " << std::endl << "  " << e.what() << std::endl;
		adj_file << "- Error:" << std::endl << "  " << ss.str();
		SignalExceptionAdjustment(ss.str(), block);
	}
	catch (const std::runtime_error& e) {
		std::stringstream ss;
		ss << "UpdateAdjustment(): Process terminated while updating the " << std::endl << "  adjustment matrices. " <<
			"Details: " << std::endl << "  " << e.what() << std::endl;
		adj_file << "- Error:" << std::endl << "  " << ss.str();
		SignalExceptionAdjustment(ss.str(), block);
	}		
	
	isPreparing_ = false;
}
		

	

void dna_adjust::PopulateEstimatedStationMatrix(const UINT32& block, UINT32& unknownParams)
{	
	// assign inner and junction station coordinates
	it_vUINT32_const _it_stn;

	double X, Y, Z;
	UINT32 j(0);

	// add inner and junction stations (all stations in phased mode are kept in v_parameterStationList_)
	for (_it_stn=v_parameterStationList_.at(block).begin();
		_it_stn!=v_parameterStationList_.at(block).end(); 
		++_it_stn)
	{
		GeoToCart<double>(
			bstBinaryRecords_.at(*_it_stn).currentLatitude,
			bstBinaryRecords_.at(*_it_stn).currentLongitude,
			bstBinaryRecords_.at(*_it_stn).currentHeight,
			&X, &Y, &Z, datum_.GetEllipsoidRef());

		// add junction stations
		switch (projectSettings_.a.adjust_mode)
		{
		case Phased_Block_1Mode:
		case PhasedMode:
			v_rigorousStations_.at(block).put(j, 0, X);
		}
		
		v_estimatedStations_.at(block).put(j, 0, X);
		v_originalStations_.at(block).put(j++, 0, X);
		
		// add junction stations
		switch (projectSettings_.a.adjust_mode)
		{
		case Phased_Block_1Mode:
		case PhasedMode:
			v_rigorousStations_.at(block).put(j, 0, Y);	
		}
		
		v_estimatedStations_.at(block).put(j, 0, Y);
		v_originalStations_.at(block).put(j++, 0, Y);
		
		// add junction stations
		switch (projectSettings_.a.adjust_mode)
		{
		case Phased_Block_1Mode:
		case PhasedMode:
			v_rigorousStations_.at(block).put(j, 0, Z);
		}
		
		v_estimatedStations_.at(block).put(j, 0, Z);
		v_originalStations_.at(block).put(j++, 0, Z);
		

		// assign appropriate variances for constraint stations
		if (bstBinaryRecords_.at((*_it_stn)).stationConst[0] == 'C')
			unknownParams--;
			
		if (bstBinaryRecords_.at((*_it_stn)).stationConst[1] == 'C')
			unknownParams--;
			
		if (bstBinaryRecords_.at((*_it_stn)).stationConst[2] == 'C')
			unknownParams--;
	}

	if (unknownParams_ == 0 && unknownsCount_ > 0)
		allStationsFixed_ = true;
}
	

void dna_adjust::PrepareStationandVarianceMatrices(const UINT32& block)
{	
	// Is this a staged adjustment for which the matrix data is to be loaded from
	// existing stage files created from a previous run?	
	if (projectSettings_.a.stage && !projectSettings_.a.recreate_stage_files && bms_meta_.reduced)
	{
		// Staged adjustments only.
		// Load matrix data from existing .mtx stage files 
		// created from a previous run.

		// redimension the estimated stations matrix
		v_estimatedStations_.at(block).setsize(v_unknownsCount_.at(block), 1);
	
		// redimension the original stations
		v_originalStations_.at(block).setsize(v_unknownsCount_.at(block), 1);
	}
	else
	{
		// Simultaneous, phased (multithreaded and staged) adjustments
		// For staged adjustments, create new .mtx stage files

		// redimension the estimated stations matrix
		v_estimatedStations_.at(block).redim(v_unknownsCount_.at(block), 1);
	
		// redimension the original stations
		v_originalStations_.at(block).redim(v_unknownsCount_.at(block), 1);
	}
	
	// Populate matrices with coordinates and
	// for phased adjustment, add junction stations
	switch (projectSettings_.a.adjust_mode)
	{
	case SimultaneousMode:
		PopulateEstimatedStationMatrix(block, unknownParams_);
		break;	
	case Phased_Block_1Mode:
	case PhasedMode:
		
		UINT32 j(static_cast<UINT32>(v_JSL_.at(block).size() * 3));

		// Is this a staged adjustment for which the matrix data is to be loaded from
		// existing stage files created from a previous run?	
		if (projectSettings_.a.stage && !projectSettings_.a.recreate_stage_files && bms_meta_.reduced)
		{
			// Staged adjustments only.
			// Load matrix data from existing .mtx stage files 
			// created from a previous run.

			// resize rigorous coordinate estimate array
			v_rigorousStations_.at(block).setsize(v_unknownsCount_.at(block), 1);
		
			// resize rigorous variances
			v_rigorousVariances_.at(block).setsize(v_unknownsCount_.at(block), v_unknownsCount_.at(block));

			// resize junction variances
			v_junctionVariances_.at(block).setsize(j, j);

			// resize junction station coordinate estimate arrays
			if (!v_blockMeta_.at(block)._blockLast && !v_blockMeta_.at(block)._blockIsolated)
			{
				v_junctionVariancesFwd_.at(block).setsize(j, j);
				v_junctionEstimatesFwd_.at(block).setsize(j, 1);		
				v_junctionEstimatesRev_.at(block+1).setsize(j, 1);
			}
		}
		else
		{		
			// Phased (multithreaded and staged) adjustments
			// For staged adjustments, create new .mtx stage files

			// resize rigorous coordinate estimate array
			v_rigorousStations_.at(block).redim(v_unknownsCount_.at(block), 1);
		
			// resize rigorous variances
			v_rigorousVariances_.at(block).redim(v_unknownsCount_.at(block), v_unknownsCount_.at(block));

			// resize junction variances
			v_junctionVariances_.at(block).redim(j, j);

			// resize junction station coordinate estimate arrays
			if (!v_blockMeta_.at(block)._blockLast && !v_blockMeta_.at(block)._blockIsolated)
			{
				v_junctionVariancesFwd_.at(block).redim(j, j);
				v_junctionEstimatesFwd_.at(block).redim(j, 1);
				v_junctionEstimatesRev_.at(block+1).redim(j, 1);
			}
		
			PopulateEstimatedStationMatrix(block, v_unknownParams_.at(block));
		}
		break;
	}
}


// Add parameter stations to measurements matrix.
void dna_adjust::PrepareDesignAndMsrMnsCmpMatrices(const UINT32& block)
{
	UINT32 pseudoMsrElemCount(0), pseudoMsrCount(0);

	switch (projectSettings_.a.adjust_mode)
	{
	case SimultaneousMode:

		// redimension the measurements matrix
		v_measMinusComp_.at(block).redim(
			v_measurementCount_.at(block), 1);
		
		// redimension the At*V-1 matrix
		v_AtVinv_.at(block).redim(
			v_unknownsCount_.at(block), 
			v_measurementCount_.at(block));

		break;
	case Phased_Block_1Mode:
	case PhasedMode:
		// Phased mode also includes sequential phased adjustment in
		// stage and multi-thread mode

		// Prepare all blocks except the first and any isolated blocks
		if (!v_blockMeta_.at(block)._blockIsolated)
		{
			// During forward and reverse (in isolation) adjustments, this block will include
			// JSLs as measurements from either v_JSL_.at(block-1) or v_JSL_.at(block).
			// During combination adjustments, this block will include JSLs as measurements 
			// from both.  So set pseudomsrJSLCount to be the sum of both and grow/shrink accordingly
			pseudoMsrCount = static_cast<UINT32>(v_JSL_.at(block).size());

			if (!v_blockMeta_.at(block)._blockFirst)
				pseudoMsrCount += static_cast<UINT32>(v_JSL_.at(block-1).size());

			pseudoMsrElemCount = pseudoMsrCount * 3;
		}

		// Is this a staged adjustment for which the matrix data is to be loaded from
		// existing stage files created from a previous run?
		if (projectSettings_.a.stage)
		{
			// Staged adjustments only.
			// Load matrix data from existing .mtx stage files 
			// created from a previous run.
			
			// redimension the measurements matrix
			if (projectSettings_.a.recreate_stage_files) {
				v_measMinusComp_.at(block).redim(
					v_measurementCount_.at(block) + pseudoMsrElemCount,	1);
			}	// real measurements + JSLs (fwd + reverse)
			else
				v_measMinusComp_.at(block).setsize(
					v_measurementCount_.at(block) + pseudoMsrElemCount,	1);

			// redimension the At*V-1 matrix
			v_AtVinv_.at(block).setsize(
				v_unknownsCount_.at(block), 
				v_measurementCount_.at(block) + pseudoMsrElemCount);	// real measurements + JSLs (fwd + reverse)
		}
		else
		{		
			// redimension the measurements matrix
			v_measMinusComp_.at(block).redim(
				v_measurementCount_.at(block) + pseudoMsrElemCount,	1);	// real measurements + JSLs (fwd + reverse)

			// redimension the At*V-1 matrix
			v_AtVinv_.at(block).redim(
				v_unknownsCount_.at(block), 
				v_measurementCount_.at(block) + pseudoMsrElemCount);	// real measurements + JSLs (fwd + reverse)
		}

		if (pseudoMsrElemCount > 0)
		{
			// Now reset these matrices back to actual measurement dimensions.
			// The forward, reverse and combination pass will use grow() to incorporate 
			// JSLs from previous and next blocks as required
			v_AtVinv_.at(block).shrink(0, pseudoMsrElemCount);
			v_measMinusComp_.at(block).shrink(pseudoMsrElemCount, 0);
		}

		// Store pseudo measurement count
		v_pseudoMeasCountFwd_.at(block) = pseudoMsrCount;
	}

	// Is this a phased adjustment in stage mode?
	if (projectSettings_.a.stage)
	{
		PrepareDesignAndMsrMnsCmpMatricesStage(block);
		return;
	}
	
	// Simultaneous, phased (multithreaded and block1) adjustments

	// Redim all matrices
	v_normals_.at(block).redim(v_unknownsCount_.at(block), v_unknownsCount_.at(block));
	
	v_design_.at(block).redim(v_measurementCount_.at(block), v_unknownsCount_.at(block));
	v_corrections_.at(block).redim(v_unknownsCount_.at(block), 1);
	v_precAdjMsrsFull_.at(block).redim(v_measurementVarianceCount_.at(block), 1);
	
	// All blocks in v_correctionsR_ are used for multi thread mode, but only the last is used
	// in single thread mode
	if (projectSettings_.a.adjust_mode == PhasedMode && 
		projectSettings_.a.adjust_mode != Phased_Block_1Mode)
	{
		if (v_blockMeta_.at(block)._blockLast) {
			v_correctionsR_.at(block).redim(v_unknownsCount_.at(block), 1);
		}
		else if (projectSettings_.a.multi_thread) {
			v_correctionsR_.at(block).redim(v_unknownsCount_.at(block), 1);
		}
	}
	

	// Now, form the design matrices
	FillDesignNormalMeasurementsMatrices(true, block, false);
	
}
	

void dna_adjust::PrepareDesignAndMsrMnsCmpMatricesStage(const UINT32& block)
{
	// Resize all matrices, which will be filled as and when needed
	v_normals_.at(block).setsize(v_unknownsCount_.at(block), v_unknownsCount_.at(block));
	v_design_.at(block).setsize(v_measurementCount_.at(block), v_unknownsCount_.at(block));
	
	// Does the matrix data need to be loaded from existing stage files
	// created from a previous run?
	if (projectSettings_.a.recreate_stage_files || !bms_meta_.reduced)
	{
		// Redim NormalsR
		v_normalsR_.at(block).redim(
			v_unknownsCount_.at(block), 
			v_unknownsCount_.at(block));

		// Redimension the corrections matrices
		v_corrections_.at(block).redim(v_unknownsCount_.at(block), 1);
		if (v_blockMeta_.at(block)._blockLast) {
			v_correctionsR_.at(block).redim(v_unknownsCount_.at(block), 1);
		}

		// Redimension precision of adjusted measurements
		v_precAdjMsrsFull_.at(block).redim(v_measurementVarianceCount_.at(block), 1);
	}
	else
	{
		// Resize the normals back-up matrix
		v_normalsR_.at(block).setsize(
			v_unknownsCount_.at(block), 
			v_unknownsCount_.at(block));

		// Resize the corrections matrices
		v_corrections_.at(block).setsize(v_unknownsCount_.at(block), 1);
		if (v_blockMeta_.at(block)._blockLast)
			v_correctionsR_.at(block).setsize(v_unknownsCount_.at(block), 1);

		// Resize precision of adjusted measurements
		v_precAdjMsrsFull_.at(block).setsize(v_measurementVarianceCount_.at(block), 1);
	}

	if (projectSettings_.a.recreate_stage_files || !bms_meta_.reduced)
		// Creating new memory mapped files?  Form the design matrices
		FillDesignNormalMeasurementsMatrices(true, block, false);
}
	

// Re-form At * V-1 for next block using estimated junction parameter station variances
// nextBlock = currentBlock+1
void dna_adjust::CarryStnEstimatesandVariancesForward(const UINT32& thisBlock, const UINT32& nextBlock)
{
	UINT32 jsl, jsl_cov;
	UINT32 jslvar, jslcovar;
	UINT32 est(0);

	it_vUINT32 _it_jsl(v_JSL_.at(thisBlock).begin()), _it_jsl_cov;

	// 1. Copy full covariance matrix and coordinate estimates of junctions to temporary
	for (_it_jsl=v_JSL_.at(thisBlock).begin(); _it_jsl!=v_JSL_.at(thisBlock).end(); ++_it_jsl)
	{
		// get index of this JSL
		jslvar = static_cast<UINT32>(std::distance(v_JSL_.at(thisBlock).begin(), _it_jsl) * 3);
		jsl = v_blockStationsMap_.at(thisBlock)[*_it_jsl] * 3;
		
		// copy variance elements for this JSL
		v_junctionVariances_.at(thisBlock).copyelements(jslvar, jslvar, 
			v_normals_.at(thisBlock), jsl, jsl, 3, 3);
		
		// copy junction estimates
		v_junctionEstimatesFwd_.at(thisBlock).copyelements(est, 0, 
			v_estimatedStations_.at(thisBlock), jsl, 0, 3, 1); 

		est += 3;

		// copy covariance elements for this JSL
		for (_it_jsl_cov=_it_jsl; _it_jsl_cov!=v_JSL_.at(thisBlock).end(); ++_it_jsl_cov)
		{
			// don't re-copy variances
			if (*_it_jsl_cov == *_it_jsl)
				continue;

			jslcovar = static_cast<UINT32>(std::distance(v_JSL_.at(thisBlock).begin(), _it_jsl_cov) * 3);
			jsl_cov = v_blockStationsMap_.at(thisBlock)[*_it_jsl_cov] * 3;
			v_junctionVariances_.at(thisBlock).copyelements(jslvar, jslcovar, 
				v_normals_.at(thisBlock), jsl, jsl_cov, 3, 3);
			v_junctionVariances_.at(thisBlock).copyelements(jslcovar, jslvar, 
				v_normals_.at(thisBlock), jsl_cov, jsl, 3, 3);
		}
	}

	if (projectSettings_.g.verbose > 5)
	{
		debug_file << "Variance matrix of junction station(s) ";
		for (_it_jsl=v_JSL_.at(thisBlock).begin(); _it_jsl!=v_JSL_.at(thisBlock).end(); ++_it_jsl)
			debug_file << bstBinaryRecords_.at(*_it_jsl).stationName << " ";
		debug_file << "carried forward: " << std::scientific << std::setprecision(16) << v_junctionVariances_.at(thisBlock);
	}	

	// 2. Perform inverse
	FormInverseVarianceMatrix(&(v_junctionVariances_.at(thisBlock)));

	// 3. Copy junction station variances for use in reverse combination adjustment
	v_junctionVariancesFwd_.at(thisBlock) = v_junctionVariances_.at(thisBlock);

	// 4. Grow msr-comp and AtVinv matrices for next block to include junction stations as measurements
	UINT32 pseudoMsrCount(static_cast<UINT32>(v_JSL_.at(thisBlock).size()));
	UINT32 pseudoMsrElemCount(pseudoMsrCount * 3);
	UINT32 paramCount(0), paramCount2, msrCountNext(v_measurementParams_.at(nextBlock));
	
	// grow matrices to accommodate JSL measurements
	v_measMinusComp_.at(nextBlock).grow(pseudoMsrElemCount, 0);
	v_AtVinv_.at(nextBlock).grow(0, pseudoMsrElemCount);
	
	// zero elements used in a prior iteration
	v_AtVinv_.at(nextBlock).zero(0, msrCountNext, v_AtVinv_.at(nextBlock).rows(), pseudoMsrElemCount);
	
	// 5. Add the temporary variance matrix to the normals of the next block, and form the 
	// msr-comp elements
	for (_it_jsl=v_JSL_.at(thisBlock).begin(); 
		_it_jsl!=v_JSL_.at(thisBlock).end(); 
		++_it_jsl, paramCount+=3)
	{
		// get index of this JSL in the next block
		jsl = v_blockStationsMap_.at(nextBlock)[*_it_jsl] * 3;		// next block
		jslvar = static_cast<UINT32>(std::distance(v_JSL_.at(thisBlock).begin(), _it_jsl) * 3);

		// add variance elements for this JSL to normals of the next block
		v_normals_.at(nextBlock).blockadd(jsl, jsl,
			v_junctionVariances_.at(thisBlock), jslvar, jslvar, 3, 3);

		// copy variances elements for this JSL to v_AtVinv_ of the next block
		v_AtVinv_.at(nextBlock).copyelements(jsl, v_measurementParams_.at(nextBlock) + paramCount,
			v_junctionVariances_.at(thisBlock), jslvar, jslvar, 3, 3);

		// msr-comp elements
		// Measured-computed for junction stations carried forward (from thisBlock-1)
		// add station X measurement
		v_measMinusComp_.at(nextBlock).put(msrCountNext++, 0, 
			(v_junctionEstimatesFwd_.at(thisBlock).get(jslvar, 0) -
			v_estimatedStations_.at(nextBlock).get(jsl, 0)));			// X (meas - computed)
		
		// add station Y measurement
		v_measMinusComp_.at(nextBlock).put(msrCountNext++, 0, 
			(v_junctionEstimatesFwd_.at(thisBlock).get(jslvar+1, 0) -
			v_estimatedStations_.at(nextBlock).get(jsl+1, 0)));		// Y (meas - computed)
		
		// add station Z measurement
		v_measMinusComp_.at(nextBlock).put(msrCountNext++, 0, 
			(v_junctionEstimatesFwd_.at(thisBlock).get(jslvar+2, 0) -
			v_estimatedStations_.at(nextBlock).get(jsl+2, 0)));		// Z (meas - computed)

		paramCount2 = 0;

		// get index of all covariances for this JSL in the next block
		for (_it_jsl_cov=_it_jsl; 
			_it_jsl_cov!=v_JSL_.at(thisBlock).end();
			++_it_jsl_cov, paramCount2+=3)
		{
			// don't re-copy variances
			if (*_it_jsl_cov == *_it_jsl)
				continue;

			jslcovar = static_cast<UINT32>(std::distance(v_JSL_.at(thisBlock).begin(), _it_jsl_cov) * 3);
			jsl_cov = v_blockStationsMap_.at(nextBlock)[*_it_jsl_cov] * 3;		// next block

			// copy covariance elements for this JSL to normals of the next block
			v_normals_.at(nextBlock).blockadd(jsl, jsl_cov,
				v_junctionVariances_.at(thisBlock), jslvar, jslcovar, 3, 3);
			v_normals_.at(nextBlock).blockadd(jsl_cov, jsl,
				v_junctionVariances_.at(thisBlock), jslcovar, jslvar, 3, 3);

			// copy covariance elements for this JSL to v_AtVinv_ of the next block
			v_AtVinv_.at(nextBlock).copyelements(jsl_cov, v_measurementParams_.at(nextBlock) + paramCount, 
				v_junctionVariances_.at(thisBlock), jslcovar, jslvar, 3, 3);

			v_AtVinv_.at(nextBlock).copyelements(jsl, v_measurementParams_.at(nextBlock) + paramCount + paramCount2, 
				v_junctionVariances_.at(thisBlock), jslvar, jslcovar, 3, 3);
		}
	}
}
	

// Re-form At * V-1 for the next block using estimated junction parameter station variances from this block (in reverse direction)
// nextBlock = thisBlock - 1
void dna_adjust::CarryStnEstimatesandVariancesReverse(const UINT32& nextBlock, const UINT32& thisBlock, bool MT_ReverseOrCombine)
{
	UINT32 jsl_var_order_this, jsl_var_order_next, jsl_covar_order_this, jsl_covar_order_next;
	UINT32 jsl_var_next, jsl_covar_next;
	UINT32 est(0);

	matrix_2d* junctionVariances(&v_junctionVariances_.at(nextBlock));
	matrix_2d* aposterioriVariances(&v_normals_.at(thisBlock));
	matrix_2d* estimatedStationsThis(&v_estimatedStations_.at(thisBlock));
	matrix_2d* estimatedStationsNext(&v_estimatedStations_.at(nextBlock));
	matrix_2d* normals(&v_normals_.at(nextBlock));
	matrix_2d* measMinusCompNext(&v_measMinusComp_.at(nextBlock));
	matrix_2d* AtVinvNext(&v_AtVinv_.at(nextBlock));

	if (MT_ReverseOrCombine)
	{
		junctionVariances = &v_junctionVariancesR_.at(nextBlock);
		aposterioriVariances = &v_normalsR_.at(thisBlock);
		estimatedStationsThis = &v_estimatedStationsR_.at(thisBlock);
		estimatedStationsNext = &v_estimatedStationsR_.at(nextBlock);
		normals = &v_normalsR_.at(nextBlock);
		measMinusCompNext = &v_measMinusCompR_.at(nextBlock);
		AtVinvNext = &v_AtVinvR_.at(nextBlock);
	}

	it_vUINT32 _it_jsl(v_JSL_.at(nextBlock).begin()), _it_jsl_cov;

	// 1. Copy full covariance matrix and coordinate estimates of junctions to temporary
	for (_it_jsl=v_JSL_.at(nextBlock).begin(); _it_jsl!=v_JSL_.at(nextBlock).end(); ++_it_jsl)
	{
		// get index of this JSL
		jsl_var_next = static_cast<UINT32>(std::distance(v_JSL_.at(nextBlock).begin(), _it_jsl) * 3);
		jsl_var_order_this = v_blockStationsMap_.at(thisBlock)[*_it_jsl] * 3;
		
		// copy junction variance elements
		junctionVariances->copyelements(jsl_var_next, jsl_var_next, 
			aposterioriVariances, jsl_var_order_this, jsl_var_order_this, 3, 3);
		
		// copy junction estimates
		v_junctionEstimatesRev_.at(thisBlock).copyelements(est, 0, 
			estimatedStationsThis, jsl_var_order_this, 0, 3, 1); 

		est += 3;

		// copy upper covariances for this JSL in the next block
		for (_it_jsl_cov=_it_jsl; _it_jsl_cov!=v_JSL_.at(nextBlock).end(); ++_it_jsl_cov)
		{
			// don't re-copy variances
			if (*_it_jsl_cov == *_it_jsl)
				continue;

			jsl_covar_next = static_cast<UINT32>(std::distance(v_JSL_.at(nextBlock).begin(), _it_jsl_cov) * 3);
			jsl_covar_order_this = v_blockStationsMap_.at(thisBlock)[*_it_jsl_cov] * 3;
			junctionVariances->copyelements(jsl_var_next, jsl_covar_next, 
				aposterioriVariances, jsl_var_order_this, jsl_covar_order_this, 3, 3);
			junctionVariances->copyelements(jsl_covar_next, jsl_var_next, 
				aposterioriVariances, jsl_covar_order_this, jsl_var_order_this, 3, 3);
		}
	}

	if (projectSettings_.g.verbose > 5)
	{
		debug_file << "Variance matrix of junction station(s) ";
		for (_it_jsl=v_JSL_.at(nextBlock).begin(); _it_jsl!=v_JSL_.at(nextBlock).end(); ++_it_jsl)
			debug_file << bstBinaryRecords_.at(*_it_jsl).stationName << " ";
		debug_file << "carried reverse: " << std::scientific << std::setprecision(16) << v_junctionVariances_.at(thisBlock);
	}	

	// 2. Perform inverse
	FormInverseVarianceMatrix(junctionVariances);

	// 3. Grow msr-comp and AtVinv matrices for next block to include junction stations as measurements
	UINT32 pseudoMsrCount(static_cast<UINT32>(v_JSL_.at(nextBlock).size()));
	UINT32 pseudoMsrElemCount(pseudoMsrCount * 3);
	UINT32 paramCount(0), paramCount2, msrCountNext(v_measurementParams_.at(nextBlock));
	
	// grow matrices to accommodate JSL measurements
	measMinusCompNext->grow(pseudoMsrElemCount, 0);
	AtVinvNext->grow(0, pseudoMsrElemCount);
	
	// zero elements used in the forward pass
	AtVinvNext->zero(0, msrCountNext, AtVinvNext->rows(), pseudoMsrElemCount);
	
	// 4. Add the temporary variance matrix to the normals of the next block, and form the 
	// msr-comp elements
	for (_it_jsl=v_JSL_.at(nextBlock).begin(); 
		_it_jsl!=v_JSL_.at(nextBlock).end(); 
		++_it_jsl, paramCount+=3)
	{
		// get index of this JSL in the next block
		jsl_var_next = static_cast<UINT32>(std::distance(v_JSL_.at(nextBlock).begin(), _it_jsl) * 3);
		jsl_var_order_this = v_blockStationsMap_.at(thisBlock)[*_it_jsl] * 3;		// previous block
		jsl_var_order_next = v_blockStationsMap_.at(nextBlock)[*_it_jsl] * 3;

		// add variance elements for this JSL to normals of the next block
		normals->blockadd(jsl_var_order_next, jsl_var_order_next,
			*junctionVariances, jsl_var_next, jsl_var_next, 3, 3);

		// copy variances elements for this JSL to v_AtVinv_ of the next block
		AtVinvNext->copyelements(jsl_var_order_next, v_measurementParams_.at(nextBlock) + paramCount,
			junctionVariances, jsl_var_next, jsl_var_next, 3, 3);

		// msr-comp elements
		// Measured-computed for junction stations carried forward (from thisBlock-1)
		// add station X measurement
		measMinusCompNext->put(msrCountNext++, 0, 
			(v_junctionEstimatesRev_.at(thisBlock).get(jsl_var_next, 0) -
			estimatedStationsNext->get(jsl_var_order_next, 0)));			// X (meas - computed)
		
		// add station Y measurement
		measMinusCompNext->put(msrCountNext++, 0, 
			(v_junctionEstimatesRev_.at(thisBlock).get(jsl_var_next+1, 0) -
			estimatedStationsNext->get(jsl_var_order_next+1, 0)));		// Y (meas - computed)
		
		// add station Z measurement
		measMinusCompNext->put(msrCountNext++, 0, 
			(v_junctionEstimatesRev_.at(thisBlock).get(jsl_var_next+2, 0) -
			estimatedStationsNext->get(jsl_var_order_next+2, 0)));		// Z (meas - computed)

		paramCount2 = 0;

		// get index of all covariances for this JSL in the next block
		for (_it_jsl_cov=_it_jsl; 
			_it_jsl_cov!=v_JSL_.at(nextBlock).end();
			++_it_jsl_cov, paramCount2+=3)
		{
			// don't re-copy variances
			if (*_it_jsl_cov == *_it_jsl)
				continue;

			jsl_covar_next = static_cast<UINT32>(std::distance(v_JSL_.at(nextBlock).begin(), _it_jsl_cov) * 3);
			jsl_covar_order_next = v_blockStationsMap_.at(nextBlock)[*_it_jsl_cov] * 3;

			// add covariance elements for this JSL to normals of the next block
			normals->blockadd(jsl_covar_order_next, jsl_var_order_next, 
				*junctionVariances, jsl_covar_next, jsl_var_next, 3, 3);

			normals->blockadd(jsl_var_order_next, jsl_covar_order_next,  
				*junctionVariances, jsl_var_next, jsl_covar_next, 3, 3);

			// copy covariance elements for this JSL to v_AtVinv_ of the next block
			AtVinvNext->copyelements(jsl_covar_order_next, v_measurementParams_.at(nextBlock) + paramCount,
				junctionVariances, jsl_covar_next, jsl_var_next, 3, 3);

			AtVinvNext->copyelements(jsl_var_order_next, v_measurementParams_.at(nextBlock) + paramCount + paramCount2,
				junctionVariances, jsl_var_next, jsl_covar_next, 3, 3);
		}
	}
}
	

// Update AtVinv based on new design matrix elements
void dna_adjust::UpdateAtVinv(pit_vmsr_t _it_msr, const UINT32& stn1, const UINT32& stn2, const UINT32& stn3, 
							UINT32& design_row, matrix_2d* design, matrix_2d* AtVinv, bool buildnewMatrices)
{
	// variance
	double variance(1./(*_it_msr)->term2);

	if (buildnewMatrices)
		if (projectSettings_.g.verbose > 5)
			debug_file << "V" << (*_it_msr)->measType << " " << 
				std::fixed << std::setprecision(16) << std::setw(26) << variance << std::endl;

	// Build  At * V-1
	AtVinv->put(stn1, design_row, variance * design->get(design_row, stn1));
	AtVinv->put(stn1+1, design_row, variance * design->get(design_row, stn1+1));
	AtVinv->put(stn1+2, design_row, variance * design->get(design_row, stn1+2));
	
	// Single station measurements
	switch ((*_it_msr)->measType)
	{
	case 'H':	// Orthometric height
	case 'I':	// Astronomic latitude
	case 'J':	// Astronomic longitude
	case 'P':	// Geodetic latitude
	case 'Q':	// Geodetic longitude
	case 'R':	// Ellipsoidal height    
		return;
	}

	// Two station measurements
	AtVinv->put(stn2, design_row, variance * design->get(design_row, stn2));
	AtVinv->put(stn2+1, design_row, variance * design->get(design_row, stn2+1));
	AtVinv->put(stn2+2, design_row, variance * design->get(design_row, stn2+2));

	switch ((*_it_msr)->measType)
	{
	case 'A':	// Horizontal angle
		AtVinv->put(stn3, design_row, variance * design->get(design_row, stn3));
		AtVinv->put(stn3+1, design_row, variance * design->get(design_row, stn3+1));
		AtVinv->put(stn3+2, design_row, variance * design->get(design_row, stn3+2));
	}
}
	

void dna_adjust::UpdateAtVinv_D(const UINT32& stn1, const UINT32& stn2, const UINT32& stn3, 
							const UINT32& angle, const UINT32& angle_count,
							UINT32& design_row, UINT32& design_row_begin,
							matrix_2d* Vinv, matrix_2d* design, matrix_2d* AtVinv)
{
	for (UINT32 j, i(0); i<3; ++i)													// for each coordinate element (x, y, z)
	{
		// add variances
		AtVinv->elementadd(stn1+i, design_row,											// station1
			design->get(design_row, stn1+i) * Vinv->get(angle, angle));
		AtVinv->elementadd(stn2+i, design_row,											// station2
			design->get(design_row, stn2+i) * Vinv->get(angle, angle));
		AtVinv->elementadd(stn3+i, design_row,											// station3
			design->get(design_row, stn3+i) * Vinv->get(angle, angle));
		
		// add covariances
		for (j=0; j<angle_count; ++j)
		{
			if (j == angle)
				continue;
			AtVinv->elementadd(stn1+i, design_row_begin+j,								// station1
				design->get(design_row, stn1+i) * Vinv->get(angle, j));
			AtVinv->elementadd(stn2+i, design_row_begin+j,								// station2
				design->get(design_row, stn2+i) * Vinv->get(angle, j));
			AtVinv->elementadd(stn3+i, design_row_begin+j,								// station3
				design->get(design_row, stn3+i) * Vinv->get(angle, j));
		}			
	}
}
	

// Re-form normals for next block using measurement variances and all parameter station variances
// Called by:
//		- AdjustPhasedReverseCombine()
//		- AdjustPhasedForward()
//		- PrepareFwdAdj (used by adjust_forward_thread)
void dna_adjust::UpdateNormals(const UINT32& block, bool MT_ReverseOrCombine)
{
	UINT32 stn1, stn2, stn3, design_row(0);
	
	it_vUINT32 _it_block_msr;
	it_vmsr_t _it_msr;

	matrix_2d* normals(&v_normals_.at(block));
	matrix_2d* design(&v_design_.at(block));
	matrix_2d* AtVinv(&v_AtVinv_.at(block));

	if (MT_ReverseOrCombine)
	{
		normals = &v_normalsR_.at(block);
		design = &v_designR_.at(block);
		AtVinv = &v_AtVinvR_.at(block);
	}	

	for (_it_block_msr=v_CML_.at(block).begin(); _it_block_msr!=v_CML_.at(block).end(); ++_it_block_msr)
	{
		if (InitialiseandValidateMsrPointer(_it_block_msr, _it_msr))
			continue;
	
		stn1 = GetBlkMatrixElemStn1(block, &_it_msr);

		switch (_it_msr->measType)
		{
		case 'H':	// Orthometric height
		case 'I':	// Astronomic latitude
		case 'J':	// Astronomic longitude
		case 'P':	// Geodetic latitude
		case 'Q':	// Geodetic longitude
		case 'R':	// Ellipsoidal height
			UpdateNormals_HIJPQR(stn1, design_row, normals, design, AtVinv);
			break;

		case 'A':	// Horizontal angle
			stn2 = GetBlkMatrixElemStn2(block, &_it_msr);
			stn3 = GetBlkMatrixElemStn3(block, &_it_msr);
			UpdateNormals_A(stn1, stn2, stn3, design_row, normals, design, AtVinv);
			break;

		case 'D':	// Direction set	
			// When a target direction is found, continue to next element.  
			if (_it_msr->vectorCount1 < 1)
				continue;
			UpdateNormals_D(block, _it_msr, design_row, normals, design, AtVinv);
			break;

		case 'B':	// Geodetic azimuth
		case 'C':	// Chord dist
		case 'E':	// Ellipsoid arc
		case 'K':	// Astronomic azimuth
		case 'L':	// Level difference
		case 'M':	// MSL arc
		case 'S':	// Slope distance
		case 'V':	// Zenith distance
		case 'Z':	// Vertical angle
			stn2 = GetBlkMatrixElemStn2(block, &_it_msr);
			UpdateNormals_BCEKLMSVZ(stn1, stn2, design_row, normals, design, AtVinv);
			break;

		case 'G':	// GPS Baseline
			// multiplies only the elements of AT * V-1 and A that contain this measurement
			stn2 = GetBlkMatrixElemStn2(block, &_it_msr);
			UpdateNormals_G(stn1, stn2, design_row, normals, AtVinv);
			break;
		
		case 'X':	// GPS Baseline cluster
			UpdateNormals_X(block, _it_msr, design_row, normals, design, AtVinv);
			break;
		
		case 'Y':	// GPS Point cluster
			UpdateNormals_Y(block, _it_msr, design_row, normals, AtVinv);
			break;		
		
		default:
			std::stringstream ss;
			ss << "UpdateNormals(): Unknown measurement type - '" << static_cast<std::string>(&(_it_msr->measType)) << 
				"'." << std::endl;
			SignalExceptionAdjustment(ss.str(), block);
		}		
	}	
}
	

void dna_adjust::AddMsrtoNormalsVar(const UINT32& design_row, const UINT32& stn,
										 matrix_2d* normals, matrix_2d* design, matrix_2d* AtVinv)
{
	// Add weighted measurement contributions to normal matrix
	for (UINT32 row, col(0); col<3; ++col)
	{
		for (row=0; row<3; ++row)
		{
			normals->elementadd(stn+row, stn+col,
				AtVinv->get(stn+row, design_row) * design->get(design_row, stn+col));
		}
	}
}
	

void dna_adjust::AddMsrtoNormalsCoVar2(const UINT32& design_row, const UINT32& stn1, const UINT32& stn2,
										 matrix_2d* normals, matrix_2d* design, matrix_2d* AtVinv)
{
	// Add covariance terms (station 1 and station 2) to normal matrix
	for (UINT32 row, col(0); col<3; ++col)
	{
		for (row=0; row<3; ++row)
		{
			// 1-2
			normals->elementadd(stn1+row, stn2+col,
				AtVinv->get(stn1+row, design_row) * design->get(design_row, stn2+col));

			normals->elementadd(stn2+row, stn1+col,
				AtVinv->get(stn2+row, design_row) * design->get(design_row, stn1+col));
		}
	}
}

void dna_adjust::AddMsrtoNormalsCoVar3(const UINT32& design_row, const UINT32& stn1, const UINT32& stn2, const UINT32& stn3,
										 matrix_2d* normals, matrix_2d* design, matrix_2d* AtVinv)
{
	// Add covariance terms (station 1, station 2, station 3) to normal matrix
	for (UINT32 row, col(0); col<3; ++col)
	{
		for (row=0; row<3; ++row)
		{
			// 1-2
			normals->elementadd(stn1+row, stn2+col,
				AtVinv->get(stn1+row, design_row) * design->get(design_row, stn2+col));

			normals->elementadd(stn2+row, stn1+col,
				AtVinv->get(stn2+row, design_row) * design->get(design_row, stn1+col));

			// 1-3
			normals->elementadd(stn1+row, stn3+col,
				AtVinv->get(stn1+row, design_row) * design->get(design_row, stn3+col));

			normals->elementadd(stn3+row, stn1+col,
				AtVinv->get(stn3+row, design_row) * design->get(design_row, stn1+col));
			
			// 2-3
			normals->elementadd(stn2+row, stn3+col,
				AtVinv->get(stn2+row, design_row) * design->get(design_row, stn3+col));

			normals->elementadd(stn3+row, stn2+col,
				AtVinv->get(stn3+row, design_row) * design->get(design_row, stn2+col));
		}
	}
}


void dna_adjust::UpdateNormals_A(const UINT32& stn1, const UINT32& stn2, const UINT32& stn3, UINT32& design_row,
										 matrix_2d* normals, matrix_2d* design, matrix_2d* AtVinv)
{
	// station 1
	AddMsrtoNormalsVar(design_row, stn1, normals, design, AtVinv);
	// station 2
	AddMsrtoNormalsVar(design_row, stn2, normals, design, AtVinv);
	// station 3
	AddMsrtoNormalsVar(design_row, stn3, normals, design, AtVinv);
	// covariance terms (station 1, station 2, station 3)
	AddMsrtoNormalsCoVar3(design_row, stn1, stn2, stn3, normals, design, AtVinv);
	design_row ++;

}
	

void dna_adjust::UpdateNormals_D(const UINT32& block, it_vmsr_t& _it_msr, UINT32& design_row,
										 matrix_2d* normals, matrix_2d* design, matrix_2d* AtVinv)
{
	UINT32 row, col, a, angle_count(_it_msr->vectorCount2 - 1);
	UINT32 skip(0), ignored(_it_msr->vectorCount1 - _it_msr->vectorCount2);

	std::vector<UINT32> stations;
	UINT32 stn1, stn2, stn3;

	vmsr_t angleRec;
	angleRec.push_back(*_it_msr);
	it_vmsr_t it_angle(angleRec.begin());

	_it_msr++;

	// Build  At * V-1 * A variances and form list of 
	// stations to assist formulation of covariances
	for (a=0; a<angle_count; ++a)																  // for each angle
	{
		// cater for ignored directions
		if (_it_msr->ignore)
		{
			while (skip < ignored)
			{
				skip++;
				_it_msr++;
				if (!_it_msr->ignore)
					break;
			}
		}

		it_angle->station3 = _it_msr->station2;

		stn1 = (v_blockStationsMap_.at(block)[it_angle->station1] * 3); 
		stn2 = (v_blockStationsMap_.at(block)[it_angle->station2] * 3);
		stn3 = (v_blockStationsMap_.at(block)[it_angle->station3] * 3);

		if (!binary_search(stations.begin(), stations.end(), it_angle->station1))
		{
			stations.push_back(it_angle->station1);
			std::sort(stations.begin(), stations.end());
		}
		if (!binary_search(stations.begin(), stations.end(), it_angle->station2))
		{
			stations.push_back(it_angle->station2);
			std::sort(stations.begin(), stations.end());
		}
		if (!binary_search(stations.begin(), stations.end(), it_angle->station3))
		{
			stations.push_back(it_angle->station3);
			std::sort(stations.begin(), stations.end());
		}

		// station 1
		for (col=0; col<3; ++col)
			for (row=0; row<3; ++row)
				normals->elementadd(stn1+row, stn1+col,
					AtVinv->get(stn1+row, design_row+a) * design->get(design_row+a, stn1+col));
		//
		// station 2
		for (col=0; col<3; ++col)
			for (row=0; row<3; ++row)
				normals->elementadd(stn2+row, stn2+col,
					AtVinv->get(stn2+row, design_row+a) * design->get(design_row+a, stn2+col));
		//
		// station 3
		for (col=0; col<3; ++col)
			for (row=0; row<3; ++row)
				normals->elementadd(stn3+row, stn3+col,
					AtVinv->get(stn3+row, design_row+a) * design->get(design_row+a, stn3+col));

		if (a+1 == angle_count)
			break;

		// prepare for next angle
		angleRec.push_back(*_it_msr);
		it_angle = angleRec.end() - 1;

		_it_msr++;	
	}

	it_vUINT32 it_stn, it_cov;
	it_angle = angleRec.begin();
	
	// Build  At * V-1 * A covariances
	for (a=0; a<angle_count; ++a)																  // for each angle
	{
		//
		// covariance terms (station 1, station 2, station 3)
		// requires std::vector<UINT32> stations which was built during formation of v_AtVinv_
		for (it_stn=stations.begin(); it_stn!=stations.end(); ++it_stn)
		{
			stn1 = (v_blockStationsMap_.at(block)[*it_stn] * 3);

			for (it_cov=it_stn; it_cov!=stations.end(); ++it_cov)
			{
				stn2 = (v_blockStationsMap_.at(block)[*it_cov] * 3);
				if (stn2 == stn1)
					continue;

				for (col=0; col<3; ++col)
				{
					for (row=0; row<3; ++row)
					{
						// 1-2
						normals->elementadd(stn1+row, stn2+col,
							AtVinv->get(stn1+row, design_row+a) * design->get(design_row+a, stn2+col));
						normals->elementadd(stn2+row, stn1+col,
							AtVinv->get(stn2+row, design_row+a) * design->get(design_row+a, stn1+col));
					}				
				}
			}
		}
		it_angle++;
	}

	design_row += angle_count;
}
	

// This function can be used for all two-station measurements
void dna_adjust::UpdateNormals_BCEKLMSVZ(const UINT32& stn1, const UINT32& stn2, UINT32& design_row,
										 matrix_2d* normals, matrix_2d* design, matrix_2d* AtVinv)
{
	// station 1
	AddMsrtoNormalsVar(design_row, stn1, normals, design, AtVinv);
	// station 2
	AddMsrtoNormalsVar(design_row, stn2, normals, design, AtVinv);
	// covariance terms (station 1 and station 2)
	AddMsrtoNormalsCoVar2(design_row, stn1, stn2, normals, design, AtVinv);
	design_row ++;
}
	

// This function can be used for all two-station measurements
void dna_adjust::UpdateNormals_HIJPQR(const UINT32& stn1, UINT32& design_row,
										 matrix_2d* normals, matrix_2d* design, matrix_2d* AtVinv)
{
	// station 1
	AddMsrtoNormalsVar(design_row, stn1, normals, design, AtVinv);
	design_row ++;
}
	

void dna_adjust::UpdateNormals_G(const UINT32& stn1, const UINT32& stn2, UINT32& design_row,
										 matrix_2d* normals, matrix_2d* AtVinv)
{
	UINT32 col, row;

	// station 2
	for (col=0; col<3; ++col)
		for (row=0; row<3; ++row)
			normals->elementadd(stn2+row, stn2+col,
				// AtVinv->get(stn2+row, design_row+col) * design->get(design_row+col, stn2+col));
				// No need to multiply by 1 as stn2 design element is always 1.  See UpdateDesignMeasMatrices_GX()
				AtVinv->get(stn2+row, design_row+col));

	// station 1
	for (col=0; col<3; ++col)
		for (row=0; row<3; ++row)
			normals->elementadd(stn1+row, stn1+col,
				// AtVinv->get(stn1+row, design_row+col) * design->get(design_row+col, stn1+col));
				// No need to multiply by -1 as stn1 design element is always -1.  See UpdateDesignMeasMatrices_GX()
				-AtVinv->get(stn1+row, design_row+col));

	// covariance terms (station 1 and station 2)
	for (col=0; col<3; ++col)
		for (row=0; row<3; ++row)
			normals->elementadd(stn1+row, stn2+col,
				// AtVinv->get(stn1+row, design_row+col) * design->get(design_row+col, stn2+col));
				// No need to multiply as stn2 is always 1.  See UpdateDesignMeasMatrices_GX()
				AtVinv->get(stn1+row, design_row+col));

	// covariance terms (station 2 and station 1)
	for (col=0; col<3; ++col)
		for (row=0; row<3; ++row)
			normals->elementadd(stn2+row, stn1+col,
				// AtVinv->get(stn2+row, design_row+col) * design->get(design_row+col, stn1+col));
				// No need to multiply by -1 as stn1 design element is always -1.  See UpdateDesignMeasMatrices_GX()
				-AtVinv->get(stn2+row, design_row+col));

	design_row += 3;
}
	

void dna_adjust::UpdateNormals_X(const UINT32& block, it_vmsr_t& _it_msr, UINT32& design_row,
										 matrix_2d* normals, matrix_2d* design, matrix_2d* AtVinv)
{
	UINT32 cluster_bsl, baseline_count(_it_msr->vectorCount1);
	UINT32 cluster_cov, covariance_count;
	UINT32 stn1, stn2, covr(0);

	vUINT32 baseline_stations;

	baseline_stations.push_back(GetBlkMatrixElemStn1(block, &_it_msr));
	baseline_stations.push_back(GetBlkMatrixElemStn2(block, &_it_msr));
	std::sort(baseline_stations.begin(), baseline_stations.end());

	matrix_2d tmp(3, 3);

	// Build  At * V-1 * A variances
	for (cluster_bsl=0; cluster_bsl<baseline_count; ++cluster_bsl)
	{
		stn1 = GetBlkMatrixElemStn1(block, &_it_msr);
		stn2 = GetBlkMatrixElemStn2(block, &_it_msr);
	
		// create a unique list of stations used in this cluster
		if (!binary_search(baseline_stations.begin(), baseline_stations.end(), stn1))
		{
			baseline_stations.push_back(stn1);
			std::sort(baseline_stations.begin(), baseline_stations.end());
		}

		if (!binary_search(baseline_stations.begin(), baseline_stations.end(), stn2))
		{
			baseline_stations.push_back(stn2);
			std::sort(baseline_stations.begin(), baseline_stations.end());
		}

		covr = cluster_bsl * 3;

		// add variances for these stations
		//normals->blockadd(stn1, stn1,							// Station 1.  This variance matrix
		//	AtVinv->submatrix(stn1,								// component must be multiplied by -1 to
		//		design_row+covr, 3, 3).scale(-1.), 0, 0, 3, 3);	// effect the design matrix elements
		
		// copy the relevant AtV-1 elements to temp and multiply 
		// by -1 to effect the design matrix elements
		AtVinv->submatrix(stn1, design_row+covr, &tmp, 3, 3);
		tmp.scale(-1.);

		// add variances for these stations
		normals->blockadd(stn1, stn1,							// Station 1.
			tmp, 0, 0, 3, 3);
		normals->blockadd(stn2, stn2,							// Station 2
			*AtVinv, stn2, design_row+covr, 3, 3);	

		covariance_count = _it_msr->vectorCount2;
		_it_msr += 3;			// move to covariances

		for (cluster_cov=0; cluster_cov<covariance_count; ++cluster_cov)
			_it_msr += 3;
	}

	matrix_2d tmp1(3, baseline_count*3), tmp2(baseline_count*3, 3);

	// Build  At * V-1 * A covariances
	for (cluster_bsl=0; cluster_bsl<baseline_stations.size(); ++cluster_bsl)
	{
		stn1 = baseline_stations.at(cluster_bsl);

		for (cluster_cov=0; cluster_cov<baseline_stations.size(); ++cluster_cov)
		{
			if (cluster_cov==cluster_bsl)
				continue;

			stn2 = baseline_stations.at(cluster_cov);

			//// add covariances for these stations
			//normals->blockadd(stn1, stn2,		
			//	AtVinv->submatrix(stn1,			
			//		design_row, 3, baseline_count*3).multiply(
			//			design->submatrix(design_row, stn2, 
			//				baseline_count*3, 3)), 0, 0, 3, 3);

			// The following is more code, but is more efficient as it
			// recycles matrix instances (as opposed to creating a
			// new instance for every call to submatrix and multiply)
			// Get AtVinv component
			AtVinv->submatrix(stn1, design_row, 
				&tmp1, 3, baseline_count*3);
			// Get design component
			design->submatrix(design_row, stn2, 
				&tmp2, baseline_count*3, 3);
			// multiply the components to form covariance
			// for stations stn1 and stn2
			//tmp.multiply(tmp1, tmp2);
			tmp.multiply(tmp1, "N", tmp2, "N");
			// Now add covariances to normals
			normals->blockadd(stn1, stn2,		
				tmp, 0, 0, 3, 3);
		}
	}

	design_row += baseline_count * 3;
}
	

void dna_adjust::UpdateNormals_Y(const UINT32& block, it_vmsr_t& _it_msr, UINT32& design_row,
										 matrix_2d* normals, matrix_2d* AtVinv)
{
	UINT32 cluster_pnt, point_count(_it_msr->vectorCount1);
	UINT32 cluster_cov, covariance_count;
	UINT32 stn1, stn2, cov_c;

	// Add to At * V-1 * A
	for (cluster_pnt=0; cluster_pnt<point_count; ++cluster_pnt)
	{		
		stn1 = GetBlkMatrixElemStn1(block, &_it_msr);
		covariance_count = _it_msr->vectorCount2;
		
		// add variance for this station
		normals->blockadd(stn1, stn1, *AtVinv, stn1, design_row, 3, 3);

		if (projectSettings_.g.verbose > 4)
			debug_file << "- Adding variances for " << bstBinaryRecords_.at(_it_msr->station1).stationName << " to v_normals_ (block " << block+1 << "), (" << stn1 << ", " << stn1 << "): " << 
				std::fixed << std::setprecision(16) << AtVinv->submatrix(stn1, design_row, 3, 3);

		if (covariance_count < 1)
		{
			design_row += 3;
			break;
		}

		_it_msr += 3;
		cov_c = 0;

		for (cluster_cov=0; cluster_cov<covariance_count; ++cluster_cov)	// number of baseline/point covariances
		{
			stn2 = GetBlkMatrixElemStn1(block, &_it_msr);		// this becomes stn1 in next cluster_pnt loop
			cov_c += 3;

			// add covariance between stn1 and this station
			normals->blockadd(stn1, stn2, *AtVinv, stn1, design_row + cov_c, 3, 3);
			normals->blockadd(stn2, stn1, *AtVinv, stn2, design_row, 3, 3);
			
			_it_msr += 3;
		}
		design_row += 3;
	}
}

void dna_adjust::BuildSimultaneousStnAppearance()
{
	it_vstn_appear _it_appear;

	v_paramStnAppearance_.resize(1);
	v_paramStnAppearance_.at(0).resize(v_blockStationsMap_.at(0).size());

	UINT32 stn(0);
	for_each(v_paramStnAppearance_.at(0).begin(), v_paramStnAppearance_.at(0).end(),
		[&stn](stn_appear& appear) {
			appear.set_id(stn++);
			appear.first_fwd();
	});
}
	
void dna_adjust::BuildUniqueBlockStationMap()
{
	v_blockStationsMapUnique_.reserve(bstBinaryRecords_.size());
	it_vvstn_appear _it_va(v_paramStnAppearance_.begin());
	it_vstn_appear _it_a;

	it_v_uint32_uint32_map _it_vm(v_blockStationsMap_.begin());
	it_uint32_uint32_map _it_m;
	
	// for each block
	for (_it_vm=v_blockStationsMap_.begin(), _it_va=v_paramStnAppearance_.begin(); 
		_it_vm!=v_blockStationsMap_.end(); 
		++_it_vm, ++_it_va)
	{
		// for each station
		for (_it_m=_it_vm->begin(), _it_a=_it_va->begin(); 
			_it_m!=_it_vm->end(); 
			++_it_m, ++_it_a)
		{
			if (_it_a->first_appearance_fwd)
			{
				v_blockStationsMapUnique_.push_back(
					u32u32_uint32_pair(
						uint32_uint32_pair(_it_m->first, _it_m->second),	// station, block index
						static_cast<UINT32>(std::distance(v_blockStationsMap_.begin(), _it_vm))));	// block
			}			
		}
	}
}

// Called by PrepareAdjustment.  Prepares the normals for a forward pass.
// The normals (prior to adding junction and parameter station variances) are
// copied to v_normalsR_.at(block) prior to this function call so that
// the reverse and combination adjustments can correctly populate the normals
// with the required junction and parameter station variances.
void dna_adjust::AddConstraintStationstoNormalsForward(const UINT32& block)
{
	it_vUINT32 _it_const;
	it_vstn_appear _it_appear;
	
	UINT32 stn;
	matrix_2d var_cart(3, 3);

	for (_it_const=v_parameterStationList_.at(block).begin(),
		_it_appear=v_paramStnAppearance_.at(block).begin(); 
		_it_const!=v_parameterStationList_.at(block).end(); 
		++_it_const, ++_it_appear)
	{

		// In phased adjustment mode, parameter station variances only 
		// need to be added once
		if (!_it_appear->first_appearance_fwd)
			continue;
		
		// compute _var_cart for this station
		FormConstraintStationVarianceMatrix(_it_const, var_cart);

		if (projectSettings_.g.verbose > 6)
			debug_file << "Parameter station " << bstBinaryRecords_.at((*_it_const)).stationName 
				<< " (" << bstBinaryRecords_.at((*_it_const)).stationConst << ") " 
				<< std::scientific << std::setprecision(16) << var_cart << std::endl;
		
		// Add the variance to the normals
		stn = v_blockStationsMap_.at(block)[(*_it_const)] * 3;
		v_normals_.at(block).blockadd(stn, stn, var_cart, 0, 0, 3, 3);
	}
}
	

// Called by PrepareAdjustmentReverse and prepares the normals for a reverse pass
// The normals (prior to adding junction and parameter station variances) are
// copied to v_normalsR_.at(block) prior to this function call so that
// the reverse and combination adjustments can correctly populate the normals
// with the required junction and parameter station variances.
void dna_adjust::AddConstraintStationstoNormalsReverse(const UINT32& block, bool MT_ReverseOrCombine)
{
	it_vUINT32 _it_const;
	it_vstn_appear _it_appear;
	
	UINT32 stn;
	matrix_2d var_cart(3, 3);
	matrix_2d* normals(&v_normals_.at(block));

	if (MT_ReverseOrCombine)
		normals = &v_normalsR_.at(block);	

	for (_it_const=v_parameterStationList_.at(block).begin(),
		_it_appear=v_paramStnAppearance_.at(block).begin(); 
		_it_const!=v_parameterStationList_.at(block).end(); 
		++_it_const, ++_it_appear)
	{
		// In phased adjustment mode, parameter station variances only 
		// need to be added once
		if (!_it_appear->first_appearance_rev)
			continue;		
		
		// compute _var_cart for this station
		FormConstraintStationVarianceMatrix(_it_const, var_cart);

		// Add the variance to the normals
		stn = v_blockStationsMap_.at(block)[(*_it_const)] * 3;
		normals->blockadd(stn, stn, var_cart, 0, 0, 3, 3);
	}
}
	

// Called by PrepareAdjustmentCombine to prepares the normals for a combination adjustment
// On entry, the normals v_normals_.at(block) will contain junction station variances from 
// forward and reverse passes, and will contain parameter station variances from the reverse
// pass.  The purpose of this function is to remove, if required, parameter station variances
// so that constraints are not applied twice.
void dna_adjust::AddConstraintStationstoNormalsCombine(const UINT32& block, bool MT_ReverseOrCombine)
{
	it_vUINT32 _it_const;
	it_vstn_appear _it_appear;
	
	UINT32 stn;
	matrix_2d var_cart(3, 3);
	matrix_2d* normals(&v_normals_.at(block));

	if (MT_ReverseOrCombine)
		normals = &v_normalsR_.at(block);	

	for (_it_const=v_parameterStationList_.at(block).begin(),
		_it_appear=v_paramStnAppearance_.at(block).begin(); 
		_it_const!=v_parameterStationList_.at(block).end(); 
		++_it_const, ++_it_appear)
	{

		// In combination adjustment, parameter station variances only need to be added 
		// if this station appears for the first time in both forward and reverse runs,
		// and subtracted if it appears for the first time in reverse but not in forward

		// If this station appears for the first time in a forward pass, then
		// constraints have not yet been added for this station in the forward pass.
		// So, continue to the next station (do nothing)
		if (_it_appear->first_appearance_fwd)
			continue;
		
		// At this point, this station has appeared in a previous block in a forward pass,
		// in which case, constraints would have been added and the effect carried through.
		// So, remove the effect.

		// compute _var_cart for this station
		FormConstraintStationVarianceMatrix(_it_const, var_cart);

		// Add the variance to the normals
		stn = v_blockStationsMap_.at(block)[(*_it_const)] * 3;

		if (!_it_appear->first_appearance_fwd)
			normals->blocksubtract(stn, stn, var_cart,
				0, 0, 3, 3);
	}
}
	

// Called by PrepareAdjustment and prepares the normals for a forward pass
// The normals (prior to adding junction and parameter station variances) are
// copied to v_normalsR_.at(block) prior to this function call so that
// the reverse and combination adjustments can correctly populate the normals
// with the required junction and parameter station variances.
void dna_adjust::AddConstraintStationstoNormalsSimultaneous(const UINT32& block)
{
	it_vUINT32 _it_const;
	
	UINT32 stn;
	matrix_2d var_cart(3, 3);

	for (_it_const=v_parameterStationList_.at(block).begin(); 
		_it_const!=v_parameterStationList_.at(block).end(); 
		++_it_const)
	{

		// compute _var_cart for this station
		FormConstraintStationVarianceMatrix(_it_const, var_cart);

		if (projectSettings_.g.verbose > 6)
            debug_file << "Parameter station " << bstBinaryRecords_.at((*_it_const)).stationName 
				<< " (" << bstBinaryRecords_.at((*_it_const)).stationConst << ") " 
                << std::scientific << std::setprecision(16) << var_cart << std::endl;
		
		// Add the variance to the normals
		stn = v_blockStationsMap_.at(block)[(*_it_const)] * 3;
		v_normals_.at(block).blockadd(stn, stn, var_cart, 0, 0, 3, 3);
	}

	if (projectSettings_.g.verbose > 6)
        debug_file << "Constrained normals " << std::scientific << std::setprecision(16) << v_normals_.at(block) << std::endl;
}
	

// used in phased adjustment to compute variances for all inner stations
void dna_adjust::FormConstraintStationVarianceMatrix(
	const it_vUINT32& _it_param_stn, matrix_2d& var_cart)
{
	vstn_t::const_iterator _it_stn(bstBinaryRecords_.begin() + (*_it_param_stn));
	
	// Constrain this station so that it will not move?
	if (_it_stn->stationConst[0] == 'C' &&
		_it_stn->stationConst[1] == 'C' &&
		_it_stn->stationConst[2] == 'C') 
	{
		// set to inverse of 'constrained' variance matrix
		var_cart = _inv_var_cart_c;
		return;
	}
	
	// free; no constraint.
	if (_it_stn->stationConst[0] == 'F' &&
		_it_stn->stationConst[1] == 'F' &&
		_it_stn->stationConst[2] == 'F')
	{
		// set to inverse of 'free' variance matrix
		var_cart = _inv_var_cart_f;
		return;
	}
	
	// If this point is reached, the station constraints will be a
	// combination of F and C.  In this case, form a variance matrix
	// in the local reference frame and propagate to the cartesian
	
	// CCF = constrained horizontally
	// FFC = constrained vertically
	
	// For CF* and FC* (irrespective of height), read on...
	//
	// If FCC is supplied for geographic coordinates, then:
	//	- latitude (north-south) is Free 
	//	- longitude (east-west) is Constrained
	//
	// If FCC is supplied for projection coordinates, then:
	//	- easting (east-west) is Free
	//	- northing (north-south) is Constrained

	matrix_2d var_local(3, 3);

	// first element
	if (_it_stn->stationConst[0] == 'F')
	{
        if (_it_stn->suppliedStationType == LLH_type_i ||
            _it_stn->suppliedStationType == LLh_type_i)
			var_local.put(1, 1, _var_F);	// Latitude:  'F'ree in the north-south (n/y) direction
		else
			var_local.put(0, 0, _var_F);	// Easting/X: 'F'ree in the east-west or X direction
	}
	else
	{
        if (_it_stn->suppliedStationType == LLH_type_i ||
            _it_stn->suppliedStationType == LLh_type_i)
			var_local.put(1, 1, _var_C);	// Latitude:  'C'onstrained in the north-south (n/y) direction
		else
			var_local.put(0, 0, _var_C);	// Easting/X: 'C'onstrained in the east-west or X direction
	}

	// second element
	if (_it_stn->stationConst[1] == 'F')
	{
        if (_it_stn->suppliedStationType == LLH_type_i ||
            _it_stn->suppliedStationType == LLh_type_i)
			var_local.put(0, 0, _var_F);
		else
			var_local.put(1, 1, _var_F);
	}
	else
	{
        if (_it_stn->suppliedStationType == LLH_type_i ||
            _it_stn->suppliedStationType == LLh_type_i)
			var_local.put(0, 0, _var_C);
		else
			var_local.put(1, 1, _var_C);
	}
	
	// ellipsoidal height - corresponds to up, whether coordinates are geographic, projection or cartesian
	if (_it_stn->stationConst[2] == 'F')
		var_local.put(2, 2, _var_F);
	else
		var_local.put(2, 2, _var_C);

	// If station is in cartesian, assume constraints are in cartesian reference frame
	if (_it_stn->suppliedStationType == XYZ_type_i)
		var_cart = var_local;
	else
		// Propagate variance in local reference frame to cartesian reference frame
		PropagateVariances_LocalCart<double>(var_local, var_cart, 
			_it_stn->currentLatitude, _it_stn->currentLongitude, 
			true);		// local to cart

	FormInverseVarianceMatrix(&var_cart);
}
	

_ADJUST_STATUS_ dna_adjust::AdjustNetwork()
{
	isAdjusting_ = true;
	isIterationComplete_ = false;
	isAdjustmentQuestionable_ = false;
	iterationCorrections_.clear_messages();

	if (projectSettings_.o._database_ids)
		LoadDatabaseId();

	switch (projectSettings_.a.adjust_mode)
	{
	case SimultaneousMode:
		if (!projectSettings_.a.report_mode)
			adj_file << std::endl << "+ Commencing simultaneous adjustment" << std::endl << std::endl;
		AdjustSimultaneous();
		return adjustStatus_;

	case PhasedMode:
		if (!projectSettings_.a.report_mode)
			adj_file << std::endl << "+ Commencing sequential phased adjustment";

		if (projectSettings_.a.multi_thread)
		{	
			if (!projectSettings_.a.report_mode)
			{
				adj_file << std::endl << "  Optimised for concurrent processing via multi-threading." << std::endl;
				adj_file << "  The active CPU supports the execution of " << std::thread::hardware_concurrency() << " concurrent threads." << std::endl;
				adj_file << std::endl;
			}
			AdjustPhasedMultiThread();
			return adjustStatus_;
		}

		if (!projectSettings_.a.report_mode)
			adj_file << std::endl << std::endl;
		AdjustPhased();

		return adjustStatus_;
	
	case Phased_Block_1Mode:
		if (!projectSettings_.a.report_mode)
		{
			adj_file << std::endl << "+ Commencing sequential phased adjustment, producing rigorous estimates" << std::endl;
			adj_file << "  for block 1 only" << std::endl << std::endl;
			adj_file << 
				"- Warning: Depending on the quality of the apriori station estimates, further" << std::endl <<
				"  iterations may be needed. --block1-phased mode should only be used once" << std::endl <<
				"  rigorous estimates have been produced for the entire network."<< std::endl << std::endl;
		}

		AdjustPhasedBlock1();
		return adjustStatus_;
	
	default:
		std::stringstream ss;
		ss << "AdjustNetwork(): Unknown adjustment type" << std::endl;
		SignalExceptionAdjustment(ss.str(), 0);
	}

	return adjustStatus_;
}
	
void dna_adjust::PrintAdjustedNetworkStations()
{
	networkadjust::DynAdjustPrinter printer(*this);
	printer.PrintAdjustedNetworkStations();
}
	
// First item in the file is a UINT32 value - the number of records in the file
// All records are of type UINT32
void dna_adjust::LoadDatabaseId()
{
	if (!projectSettings_.o._database_ids)
		return;

	if (databaseIDsLoaded_)
		return;

	std::string dbid_filename = formPath<std::string>(projectSettings_.g.output_folder,
		projectSettings_.g.network_name, "dbid");

	// When printing database ids, force printing adjusted measurements 
	// as a contiguous list in original sort order.  Why?  To simplify 
	// reading adj file when loading adjusted measurement info to
	// the database.
	projectSettings_.o._output_msr_blocks = 0;
	
	// Do not print computed measurements
	//projectSettings_.o._cmp_msr_iteration = 0;
	
	v_msr_db_map_.clear();
	
	std::ifstream dbid_file;
	try {
		// Create geoid file.  Throws runtime_error on failure.
		file_opener(dbid_file, dbid_filename,
			std::ios::in | std::ios::binary, binary);
	}
	catch (const std::runtime_error& e) {
		SignalExceptionAdjustment(e.what(), 0);
	}
	
	UINT32 r, recordCount;
	msr_database_id_map rec;
	
	try {
		// get size and reserve vector size
		dbid_file.read(reinterpret_cast<char *>(&recordCount), sizeof(UINT32));
		v_msr_db_map_.reserve(recordCount);
		
		UINT16 val;

		for (r = 0; r < recordCount; r++)
		{
			// Read data
			dbid_file.read(reinterpret_cast<char *>(&rec.msr_id), sizeof(UINT32));
			dbid_file.read(reinterpret_cast<char *>(&rec.cluster_id), sizeof(UINT32));
			dbid_file.read(reinterpret_cast<char *>(&val), sizeof(UINT16));
			rec.is_msr_id_set = val_uint<bool, UINT16>(val);
			dbid_file.read(reinterpret_cast<char *>(&val), sizeof(UINT16));
			rec.is_cls_id_set = val_uint<bool, UINT16>(val);
			
			// push back
			v_msr_db_map_.push_back(rec);
		}

		dbid_file.close();
		databaseIDsLoaded_ = true;
	}
	catch (const std::ifstream::failure& f) {
		SignalExceptionAdjustment(f.what(), 0);
	}
	catch (const XMLInteropException& e)  {
		SignalExceptionAdjustment(e.what(), 0);
	}
}

void dna_adjust::FormUniqueMsrList()
{
	// Determine actual CML size
	size_t cml_size(0);
	for_each(v_CML_.begin(), v_CML_.end(),
		[&cml_size](vUINT32& blockCML){
			cml_size += blockCML.size();
	});

	UINT32 block(0);
	it_vvUINT32 _it_cml;
	it_vUINT32 _it_block_msr;
	
	v_msr_block_.clear();
	v_msr_block_.reserve(cml_size);
	
	// Fill vector of msr-block pairs
	for (_it_cml=v_CML_.begin(); _it_cml!=v_CML_.end(); ++_it_cml, ++block)
		for (_it_block_msr=_it_cml->begin(); _it_block_msr!=_it_cml->end(); ++_it_block_msr)
			v_msr_block_.push_back(
				uint32_u32u32_pair(*_it_block_msr,		// msr index
					uint32_uint32_pair(block, 0)));		// block, precision adj msr row

	// Is the index to the first row of the precision of adjusted measurements vector needed?
	if (projectSettings_.o._adj_gnss_units == XYZ_adj_gnss_ui)
		return;

	it_vmsr_t _it_msr;
	_it_uint32_u32u32_pair _it_b_pam;
	UINT32 precadjmsr_row(0);
	block = UINT_MAX;

	for (_it_b_pam=v_msr_block_.begin(); _it_b_pam!=v_msr_block_.end(); ++_it_b_pam)
	{	
		// new block?  Initialise row counter to the first row
		if (block != _it_b_pam->second.first)
		{
			precadjmsr_row = 0;
			block = _it_b_pam->second.first;
		}
		
		// Set index
		_it_b_pam->second.second = precadjmsr_row;

		// Get the next measurement and determine how many rows until the next index
		_it_msr = bmsBinaryRecords_.begin() + (_it_b_pam->first);

		switch (_it_msr->measType)
		{
		// Single row measurements
		case 'A':
		case 'B':	// Geodetic azimuth
		case 'C':	// Chord dist
		case 'E':	// Ellipsoid arc
		case 'K':	// Astronomic azimuth
		case 'L':	// Level difference
		case 'M':	// MSL arc
		case 'S':	// Slope distance
		case 'V':	// Zenith distance
		case 'Z':	// Vertical angle
		case 'H':	// Orthometric height
		case 'I':	// Astronomic latitude
		case 'J':	// Astronomic longitude
		case 'P':	// Geodetic latitude
		case 'Q':	// Geodetic longitude
		case 'R':	// Ellipsoidal height
			precadjmsr_row++;
			break;
		// Single row per subtended angle (# directions less 1)
		case 'D':	// Direction set
			precadjmsr_row += (_it_msr->vectorCount2 - 1);
			break;
		// Three rows (6 elements for upper triangular)
		case 'G':	// GPS Baseline
			precadjmsr_row += 6;
			break;
		// Three rows (6 elements for upper triangular) per
		// cluster measurement
		case 'X':	// GPS Baseline cluster
		case 'Y':	// GPS Point cluster
			precadjmsr_row += _it_msr->vectorCount1 * 6;
			break;
		}
		
	}

}
	
void dna_adjust::CreateMsrToStnTally()
{
	// 1. Load the aml file. The aml file is needed to determine the measurements
	// connected to each station.
	try {
		// Load aml file.  Throws runtime_error on failure.
		AmlFile aml;
		aml.load_aml_file(projectSettings_.s.aml_file, &vAssocMsrList_, &bmsBinaryRecords_);
	}
	catch (const std::runtime_error& e) {
		SignalExceptionAdjustment(e.what(), 0);
	}

	v_stnmsrTally_.clear();
	v_stnmsrTally_.resize(vAssocStnList_.size());

	_it_vasl _it_asl;
	UINT32 i, m, amlIndex, msrCount;

	for (_it_asl=vAssocStnList_.begin(), i=0; _it_asl!=vAssocStnList_.end(); ++_it_asl, ++i)
	{
		msrCount = _it_asl->GetAssocMsrCount();

		for (m=0; m<msrCount; ++m)
		{
			amlIndex = _it_asl->GetAMLStnIndex() + m;

			if (bmsBinaryRecords_.at(vAssocMsrList_.at(amlIndex).bmsr_index).ignore)
				continue;

			v_stnmsrTally_.at(i).IncrementMsrType(bmsBinaryRecords_.at(vAssocMsrList_.at(amlIndex).bmsr_index).measType);
		}

		v_stnmsrTally_.at(i).TotalCount();
	}
}
	
void dna_adjust::CloseOutputFiles()
{
	adj_file.close();
	xyz_file.close();

	if (projectSettings_.g.verbose > 0)
		debug_file.close();
}
	

void dna_adjust::AdjustSimultaneous()
{
	adjustStatus_ = ADJUST_SUCCESS;
	initialiseIteration();
	
	std::ostringstream ss;
	std::string corr_msg;
	std::chrono::milliseconds elapsed_time(std::chrono::milliseconds(0));
	cpu_timer it_time, tot_time;

	bool iterate;

	for (UINT32 i=0; i<projectSettings_.a.max_iterations; ++i)
	{
		if (IsCancelled())
			break;

		isIterationComplete_ = false;

		blockLargeCorr_ = 0;
		largestCorr_ = 0.0;

		// Print the iteration # to adj file
		printer_->PrintIteration(incrementIteration());

		// Does the user want to print computed measurements on each iteration?
		if (projectSettings_.o._cmp_msr_iteration)
			printer_->PrintCompMeasurements(0, "a-priori");		

		ss.str("");
		
		it_time.start();

		// Least Squares Solution
		// Inverse is only required if:
		//	- This is the first iteration
		//	- The network contains non-GPS measurements, in which an updated 
		//	  normals matrix would be available based upon reformed partial
		//	  derivatives in the design matrix
		SolveTry(CurrentIteration() < 2 || v_msrTally_.at(0).ContainsNonGPS());

		// calculate and print total time
		printer_->PrintAdjustmentTime(it_time, iteration_time);
		
		// Add corrections to estimates
		v_estimatedStations_.at(0).add(v_corrections_.at(0));
		
		// Compute and print largest correction
		maxCorr_ = v_corrections_.at(0).compute_maximum_value();
		OutputLargestCorrection(corr_msg);

		// update data for messages
		iterationCorrections_.add_message(corr_msg);
		iterationQueue_.push_and_notify(CurrentIteration());	// currentIteration begins at 1, so not zero-indexed
		isIterationComplete_ = true;
		
		// continue iterating?
		iterate = !IsCancelled() && fabs(maxCorr_) > projectSettings_.a.iteration_threshold;
		if (!iterate)
			break;

		// Does the user want to print statistics on each iteration?
		if (projectSettings_.o._adj_stat_iteration)
		{
			// Compute network statistics
			ComputeStatisticsOnIteration();
			
			// Print statistics summary to adj file
			printer_->PrintStatistics(false);
		}

		// Does the user want to print adjusted measurements
		// on each iteration?
		if (projectSettings_.o._adj_msr_iteration)
			ComputeandPrintAdjMsrOnIteration();

		// Does the user want to print adjusted station coordinates
		// on each iteration?
		if (projectSettings_.o._adj_stn_iteration)
			// computes geographic coordinates if required
			printer_->PrintBlockStations(adj_file, 0, &v_estimatedStations_.at(0), &v_normals_.at(0), 
				false, !v_msrTally_.at(0).ContainsNonGPS(), !v_msrTally_.at(0).ContainsNonGPS(), true, false);

		// Update normals and measured-computed matrices for the next iteration.
		UpdateAdjustment(iterate);
	}

	ValidateandFinaliseAdjustment(tot_time);
}

void dna_adjust::ValidateandFinaliseAdjustment(cpu_timer& tot_time)
{
	isAdjusting_ = false;

	if (adjustStatus_ > ADJUST_TEST_FAILED)
		return;

	if (IsCancelled())
	{
		adjustStatus_ = ADJUST_CANCELLED;
		return;
	}

	if (CurrentIteration() == projectSettings_.a.max_iterations &&
		fabs(maxCorr_) > projectSettings_.a.iteration_threshold)
		adjustStatus_ = ADJUST_MAX_ITERATIONS_EXCEEDED;

	// Back up simultaneous rigorous variance estimates (for serialising 
	// to disk at SerialiseAdjustedVarianceMatrices() ), so that executing adjust
	// in report-results mode has access to the latest variance estimates
	switch (projectSettings_.a.adjust_mode)
	{
	case SimultaneousMode:
		v_rigorousVariances_.at(0) = v_normals_.at(0);
		break;
	case Phased_Block_1Mode:
	case PhasedMode:
	default:
		break;
	}

	// Print status
	printer_->PrintAdjustmentStatus();
	// Compute and print time taken to run adjustment
	printer_->PrintAdjustmentTime(tot_time, total_time);
}
	
void dna_adjust::PrintAdjustmentTime(cpu_timer& time, _TIMER_TYPE_ timerType)
{
	// Store total time if needed
	if (timerType != iteration_time) {
		auto elapsed = time.elapsed();
		total_time_ = std::chrono::milliseconds(elapsed.wall.count() / MILLI_TO_NANO);
	}
	
	printer_->PrintAdjustmentTime(time, static_cast<int>(timerType));
}

void dna_adjust::AdjustPhased()
{
	initialiseIteration();

	std::string corr_msg;
	std::ostringstream ss;
	UINT32 i;
	bool iterate(true);

	cpu_timer it_time, tot_time;
		
	// do until convergence criteria is met
	for (i=0; i<projectSettings_.a.max_iterations; ++i)
	{
		if (IsCancelled())
			break;

		isIterationComplete_ = false;

		blockLargeCorr_ = 0;
		largestCorr_ = 0.0;
		maxCorr_ = 0.0;

		// initialise potential outlier count and statistical
		// quantities
		potentialOutlierCount_ = 0;
		chiSquaredStage_ = 0.;
		measurementParams_ = 0;
	
		// Print the iteration # to adj file
		printer_->PrintIteration(incrementIteration());
		
		it_time.start();

		AdjustPhasedForward();
		if (IsCancelled())
			break;

		AdjustPhasedReverseCombine();
		if (IsCancelled())
			break;

		// calculate and print total time
		printer_->PrintAdjustmentTime(it_time, iteration_time);
		
		// Calculate and print largest adjustment correction and station ID
		OutputLargestCorrection(corr_msg);
		
		iterationCorrections_.add_message(corr_msg);
		iterationQueue_.push_and_notify(CurrentIteration());	// currentIteration begins at 1, so not zero-indexed
		isIterationComplete_ = true;

		// Continue iterating?
		iterate = !IsCancelled() && fabs(maxCorr_) > projectSettings_.a.iteration_threshold;
		if (!iterate)
			break;

		// Update normals and measured-computed matrices for the next iteration.
		// Similar to PrepareAdjustment, UpdateAdjustment prepares every block
		// in the network so that forward and reverse adjustments can commence
		// at the same time.
		UpdateAdjustment(iterate);	
		if (IsCancelled())
			break;

		// Does the user want to print statistics on each iteration?
		if (projectSettings_.o._adj_stat_iteration)
		{
			// Compute network statistics
			ComputeStatisticsOnIteration();
			
			// Print statistics summary to adj file
			printer_->PrintStatistics(false);
		}

		// Does the user want to print adjusted measurements
		// on each iteration?
		if (projectSettings_.o._adj_msr_iteration)
			ComputeandPrintAdjMsrOnIteration();
	}

	chiSquared_ = chiSquaredStage_;

	ValidateandFinaliseAdjustment(tot_time);
}
	

// In the phased block 1 approach, only one iteration is performed since
// the adjustment assumes that rigorous estimates have already been produced
void dna_adjust::AdjustPhasedBlock1()
{
	initialiseIteration(1);
	
	std::ostringstream ss;
	std::string corr_msg;

	adjustStatus_ = ADJUST_SUCCESS;
	
	maxCorr_ = 0.0;
	blockLargeCorr_ = 0;
	largestCorr_ = 0.0;
		
	adj_file << std::endl << std::endl;

	if (projectSettings_.g.verbose)
		debug_file << std::endl << std::endl;
		
	cpu_timer it_time, tot_time;
		
	// Only the reverse adjustment is needed to achieve rigorous
	// estimates for block 1
	AdjustPhasedReverse();

	blockLargeCorr_ = 0;
	largestCorr_ = v_corrections_.at(0).maxvalue();

	// calculate and print total time
	printer_->PrintAdjustmentTime(it_time, iteration_time);	
	
	// Compute and print largest correction for block 1 only
	maxCorr_ = v_corrections_.at(0).compute_maximum_value();
	OutputLargestCorrection(corr_msg);

	if (fabs(maxCorr_) > projectSettings_.a.iteration_threshold)
		adjustStatus_ = ADJUST_THRESHOLD_EXCEEDED;
		
	iterationCorrections_.add_message(corr_msg);
	iterationQueue_.push_and_notify(CurrentIteration());	// currentIteration begins at 1, so not zero-indexed

	ValidateandFinaliseAdjustment(tot_time);
}
	

// Used to rebuild normals for stage adjustments
void dna_adjust::RebuildNormals(const UINT32 block, adjustOperation direction, bool AddConstraintStationstoNormals, bool BackupNormals)
{
	// Update measurements-computed vector using new estimates
	FillDesignNormalMeasurementsMatrices(false, block, false);

	// Update normal equations and add parameter station variances
	v_normals_.at(block).zero();
	UpdateNormals(block, false);

	if (!AddConstraintStationstoNormals)
		return;

	// Add parameter stations
	switch (direction)
	{
	case __forward__:

		// Back up normals.  This copy contains the contributions from all
		// apriori measurement variances, excluding parameter station 
		// variances and junction station variances
		if (BackupNormals)
			v_normalsR_.at(block) = v_normals_.at(block);

		AddConstraintStationstoNormalsForward(block);
		break;

	case __reverse__:
		AddConstraintStationstoNormalsReverse(block, false);
		break;
	default:
		break;
	}
}
	
 
void dna_adjust::AdjustPhasedForward()
{	
	forward_ = true;

	UINT32 currentBlock(0);

	// For staged adjustments, load the first block from mapped memory file
	if (projectSettings_.a.stage)
	{
		DeserialiseBlockFromMappedFile(currentBlock);
		RebuildNormals(currentBlock, __forward__, true, true);
		if (IsCancelled())
			return;
	}

	_it_uint32_u32u32_pair begin, end;
	begin = v_msr_block_.begin();

	for (currentBlock=0; currentBlock<blockCount_; ++currentBlock)
	{		
		if (IsCancelled())
			break;

		SetcurrentBlock(currentBlock);

		// Does the user want to print computed measurements?
		if (projectSettings_.o._cmp_msr_iteration)
			printer_->PrintCompMeasurements(currentBlock, "a-priori");
		
		if (projectSettings_.o._adj_stn_iteration || projectSettings_.o._cmp_msr_iteration)
		{
			adj_file_mutex.lock();
			adj_file << std::endl << std::left << "Adjusting ";
			if (v_blockMeta_.at(currentBlock)._blockIsolated)
				adj_file << "(isolated) ";
			adj_file << "block " << currentBlock+1;
			if (v_blockMeta_.at(currentBlock)._blockLast || v_blockMeta_.at(currentBlock)._blockIsolated)
				adj_file << " (forward, rigorous)... ";
			else
				adj_file << " (forward, in isolation)... ";
			adj_file_mutex.unlock();
		}
		
		// At this point, whether first iteration or not, if currentBlock is the first block, 
		// the normals will have been initialised.  For all later blocks, the normals will contain
		// the contribution of junction station coordinates and variances from preceding blocks.
		// In either case, the block is ready for adjustment.  Junction station coordinates and 
		// variances are carried forward below

		// Least Squares Solution
		SolveTry(true, currentBlock);

		// Does the user want to print adjusted measurements
		// on each iteration?
		if (projectSettings_.o._adj_msr_iteration || 
			projectSettings_.o._adj_stn_iteration || 
			projectSettings_.o._cmp_msr_iteration)
			adj_file << " done." << std::endl;

		// Add corrections to estimates and update degrees of freedom.
		UpdateEstimatesForward(currentBlock);

		// Compute and print adjusted measurements for this 
		// block during an iteration
		if (projectSettings_.o._adj_msr_iteration)
		{
			end = begin + v_CML_.at(currentBlock).size();
			ComputeandPrintAdjMsrBlockOnIteration(currentBlock, v_uint32_u32u32_pair(begin, end), true);
			begin = end;
		}
		
		// Debug and diagnose (if required)
		debug_BlockInformation(currentBlock, (forward_ ? " (Forward)" : " (Reverse)"));
		
		// OK, now shrink matrices back to normal size
		ShrinkForwardMatrices(currentBlock);

		// Carry the estimated junction station coordinates and variances 
		// to the next block for applicable blocks only.  For staged
		// Deserialise the next block from mapped files and rebuild 
		// normals
		CarryForwardJunctions(currentBlock, currentBlock+1);

		// For staged adjustments, write to disk and unload matrix data
		if (projectSettings_.a.stage)
			// Don't offload last block since the reverse adjustment 
			// will need this block
			if (currentBlock < (blockCount_ - 1))
				OffloadBlockToMappedFile(currentBlock);
	}
}
	

void dna_adjust::PurgeMatricesFromDisk()
{
	// remove all the files
	for_each(
		v_stageFileStreams_.begin(), 
		v_stageFileStreams_.end(), 
		[](std::string str){			// use lambda expression
			
			// remove mapping
			boost::interprocess::file_mapping::remove(str.c_str());
			// If the file still exists, remove it!
			if (std::filesystem::exists(str))
			boost::interprocess::file_mapping::remove(str.c_str());
		}
	);	
}
	

void dna_adjust::PrepareAdjustmentBlock(const UINT32 block, const UINT32 thread_id)
{

	if (projectSettings_.a.multi_thread && projectSettings_.g.verbose > 2)
	{
		adj_file_mutex.lock();
		adj_file << thread_id << "> Preparing block " << block + 1 << "..." << std::endl;
		adj_file_mutex.unlock();
	}

	UINT32 i;
	it_vUINT32 _it_stn;

	if (projectSettings_.g.verbose > 3)
	{
		debug_file << "BLOCK " << block + 1 << " Station summary -------------------------------------------------" << std::endl;
		debug_file << "Stations: ";
		for (i=0; i<v_parameterStationList_.at(block).size(); ++i)
		{
			debug_file << v_parameterStationList_.at(block).at(i);
			if ((_it_stn = find(v_ISL_.at(block).begin(), v_ISL_.at(block).end(), v_parameterStationList_.at(block).at(i))) != v_ISL_.at(block).end())
				debug_file << "i ";
			else if ((_it_stn = find(v_JSL_.at(block).begin(), v_JSL_.at(block).end(), v_parameterStationList_.at(block).at(i))) != v_JSL_.at(block).end())
				debug_file << "j ";
			else 
				debug_file << " ";				
		}
		debug_file << " (";				
		for (i=0; i<v_parameterStationList_.at(block).size(); ++i)
		{
			if (i>0)
					debug_file << ", ";
			debug_file << bstBinaryRecords_.at(v_parameterStationList_.at(block).at(i)).stationName;				
		}
		debug_file << ")" << std::endl << std::endl;
	}

	// Form Estimated Stations matrices. For phased adjustments, matrices 
	// for junction estimates and variances carried forward and reverse 
	// are also formed.  Initial geographic coordinates are computed also.
	PrepareStationandVarianceMatrices(block);

	// Dimension matrices for design, msr-comp, AtVinv, normals, corrections and 
	// precision of adjusted msrs.  msr-comp and AtVinv have additional elements
	// for junction station measurements.
	// Form design matrix from observation equations using AddMsrtoNormalandMeasMatrices.
	// In this case, the more 'difficult' measurement types are reduced to simpler
	// types to allow direct expression of measurements in terms of X, Y, Z.
	// For ever iteration thereafter, RecomputeMeasurements... is used
			
	// For phased adjustment, all stations in ISL and JSL are treated as parameters 
	// to be estimated.  The inclusion of all stations is required since:
	//		"as a consequence of the segmentation algorithm, many blocks 
	//		are formed where there are insufficient connections to datum 
	//		for the section as a whole or, some stations coordinates within 
	//		the section cannot be estimated from the set of measurements 
	//		making up the section."
	//											(Leahy, PhD thesis, p. 3.21).	
	//
	// To begin with, the number of rows in the block matrices corresponds with the ISL.  To facilitate
	// the above, rows must be added to all blocks (not including the first) so as to
	// cater for junction station estimates and variances carried forward, reverse and during the combination stage.
	// Depending on whether an "in isolation" or "rigorous" adjustment is being performed, the matrix 
	// dimensions are "grown" or "shrunk" accordingly (without altering the containing data in buffer).

	PrepareDesignAndMsrMnsCmpMatrices(block);
	
	// Is this a staged adjustment for which the matrix data is to be loaded from
	// existing stage files created from a previous run?	
	if (projectSettings_.a.stage && !projectSettings_.a.recreate_stage_files && bms_meta_.reduced)
	{
		// Staged adjustments only.
		// Load matrix data from existing .mtx stage files 
		// created from a previous run.
		PrepareMappedRegions(block);
		return;
	}

	if (!projectSettings_.a.stage)
	{
		// Back up normals.  This copy contains the contributions from all apriori measurement
		// variances, excluding parameter station variances and junction station variances
		v_normalsR_.at(block) = v_normals_.at(block);

		if (projectSettings_.a.multi_thread) {
			v_normalsRC_.at(block).redim(
				v_normalsR_.at(block).rows(), v_normalsR_.at(block).columns());
		}

		switch (projectSettings_.a.adjust_mode)
		{
		case Phased_Block_1Mode:
		case PhasedMode:
			AddConstraintStationstoNormalsForward(block);
			break;
		case SimultaneousMode:
			AddConstraintStationstoNormalsSimultaneous(block);
			break;
		}
	}
	else
	{
		// Offload matrices to new .mtx stage files

		// Set memory mapped file region offsets for all 
		// matrices associated with this block
		SetRegionOffsets(block);
		// Purge from memory all matrices associated with this
		// block and write to disk (the memory mapped file)
		OffloadBlockToDisk(block);

		// A return here is needed to prevent any attempt to 
		// print v_design_, v_AtVinv_ or v_normals_ in verbose mode
		return;
	}

	if (projectSettings_.a.multi_thread && projectSettings_.g.verbose > 2)
	{
		adj_file_mutex.lock();
		adj_file << thread_id << "> Done." << std::endl;
		adj_file_mutex.unlock();
	}

	if (projectSettings_.g.verbose > 3)
	{
		debug_file << std::endl;
		debug_file << "Design " << std::scientific << std::setprecision(16) << v_design_.at(block) << std::endl;
		debug_file << "AtVinv " << std::scientific << std::setprecision(16) << v_AtVinv_.at(block) << std::endl;
		debug_file << "Normals " << std::scientific << std::setprecision(16) << v_normals_.at(block) << std::endl;
	}
}
	

void dna_adjust::ShrinkForwardMatrices(const UINT32 currentBlock)
{
	UINT32 pseudoMsrElemCount;
	if (!v_blockMeta_.at(currentBlock)._blockIsolated && !v_blockMeta_.at(currentBlock)._blockFirst)
	{
		pseudoMsrElemCount = static_cast<UINT32>(v_JSL_.at(currentBlock-1).size() * 3);

		v_measMinusComp_.at(currentBlock).shrink(pseudoMsrElemCount, 0);
		v_AtVinv_.at(currentBlock).shrink(0, pseudoMsrElemCount);
	}
}	
	

// Called by:
// - adjust_forward_thread::operator()()
void dna_adjust::UpdateEstimatesForward(const UINT32 currentBlock)
{
	// update station coordinates with lsq-estimated corrections
	v_estimatedStations_.at(currentBlock).add(v_corrections_.at(currentBlock));

	// compute degrees of freedom (this is only performed once here during forward pass)
	v_statSummary_.at(currentBlock)._degreesofFreedom = 
		v_measurementParams_.at(currentBlock) + 
		(v_pseudoMeasCountFwd_.at(currentBlock) * 3) - 
		v_unknownParams_.at(currentBlock);

	if (v_blockMeta_.at(currentBlock)._blockLast || v_blockMeta_.at(currentBlock)._blockIsolated)
	{
		// update max correction
		if (fabs(v_corrections_.at(currentBlock).compute_maximum_value()) > fabs(maxCorr_))
			SetmaxCorr(v_corrections_.at(currentBlock).maxvalue());
		
		// update largest correction
		if (fabs(v_corrections_.at(currentBlock).maxvalue()) > fabs(largestCorr_))
		{
			largestCorr_ = v_corrections_.at(currentBlock).maxvalue();
			blockLargeCorr_ = currentBlock;
		}

		// Now copy 'estimated' coordinates to 'rigorous' for comparison on the next iteration
		v_rigorousStations_.at(currentBlock) = v_estimatedStations_.at(currentBlock);
		v_rigorousVariances_.at(currentBlock) = v_normals_.at(currentBlock);

		// Temporarily hold corrections for last block.  These corrections will be restored
		// in UpdateEstimatesFinal()
		if (!projectSettings_.a.multi_thread)
			if (v_blockMeta_.at(currentBlock)._blockLast)
				v_correctionsR_.at(currentBlock) = v_corrections_.at(currentBlock);
	}
}
	

// multi-thread
// Add junction stations from prevBlock to thisBlock as measurements
// thisBlock = prevBlock + 1
// this function should not be called if from block is the first block
void dna_adjust::CarryForwardJunctions(const UINT32 thisBlock, const UINT32 nextBlock)
{
	// If required, print the stations for this block.  Geographic coordinates are updated
	// only if necessary for output
	if (projectSettings_.o._adj_stn_iteration)
	{
		adj_file_mutex.lock();

		if (projectSettings_.a.multi_thread)
		{	
			adj_file << std::endl << "1> Adjusted block " << thisBlock + 1;
			if (v_blockMeta_.at(thisBlock)._blockLast || v_blockMeta_.at(thisBlock)._blockIsolated)
				adj_file << " (forward, rigorous) " << std::endl;
			else
				adj_file << " (forward, in isolation) " << std::endl;
		}

		printer_->PrintBlockStations(adj_file, thisBlock, 
			&v_estimatedStations_.at(thisBlock), &v_normals_.at(thisBlock), 
			false, true, false, true, false);
		adj_file_mutex.unlock();
	}

	// Carry the estimated junction station coordinates and variances to the next block if
	// this block and the next block are both in the same network.  Thus, carry data except when:
	//	- thisBlock is a single block, or
	//	- nextBlock is a single block, or
	//	- thisBlock is the last block
	if (v_blockMeta_.at(thisBlock)._blockIsolated)
		return;
	if (v_blockMeta_.at(thisBlock)._blockLast)
		return;
	if (v_blockMeta_.at(nextBlock)._blockIsolated)
		return;
	
	// For staged adjustments, load block info for nextBlock
	if (projectSettings_.a.stage)
	{
		DeserialiseBlockFromMappedFile(nextBlock);
		RebuildNormals(nextBlock, __forward__, true, true);
	}

	CarryStnEstimatesandVariancesForward(thisBlock, nextBlock);
}
	

bool dna_adjust::PrepareAdjustmentReverse(const UINT32 currentBlock, bool MT_ReverseOrCombine)
{
	// if this is a single block, then there is no need to perform a reverse adjustment
	if (v_blockMeta_.at(currentBlock)._blockIsolated)
	{
		if (projectSettings_.o._adj_stn_iteration)
		{
			adj_file_mutex.lock();
			adj_file << "Skipping reverse adjustment of (isolated) block " << currentBlock+1 << ".  Rigorous coordinates produced in forward pass. " << std::endl << std::endl;
			adj_file_mutex.unlock();
		}

		if (projectSettings_.a.stage)
			OffloadBlockToMappedFile(currentBlock);

		return false;
	}

	// If currentBlock is not the last block, then return true so that 
	// AdjustPhasedReverseCombine can proceed.  In this instance, the
	// Normals were already formed in CarryReverseEstimates
	if (!v_blockMeta_.at(currentBlock)._blockLast)
		return true;

	// OK. currentBlock is the last block, and this is the commencement of
	// a reverse run, so reset coordinates
	matrix_2d* estimatedStations(&v_estimatedStations_.at(currentBlock));	

	if (MT_ReverseOrCombine)
	{
		estimatedStations = &v_estimatedStationsR_.at(currentBlock);
	}
	else
	{
		// Restore back up copy of normals for reverse adjustment in 
		// single thread mode
		v_normals_.at(currentBlock) = v_normalsR_.at(currentBlock);

	}

	// Reset coordinates to original values
	*estimatedStations = v_originalStations_.at(currentBlock);
	
	AddConstraintStationstoNormalsReverse(currentBlock, MT_ReverseOrCombine);

	return true;
}
	

// Back up normals prior to inversion for combination adjustment
void dna_adjust::BackupNormals(const UINT32 block, bool MT_ReverseOrCombine)
{
	// This procedure is only required if this block requires a combination
	// adjustment to produce rigorous estimates
	if (!CombineRequired(block))
		return;

	matrix_2d* normals(&v_normals_.at(block));
	matrix_2d* normalsCopy(&v_normalsR_.at(block));

	if (MT_ReverseOrCombine)
	{
		normals = &v_normalsR_.at(block);
		normalsCopy = &v_normalsRC_.at(block);
	}

	normalsCopy->copyelements(0, 0, normals, 0, 0,
		normals->rows(), normals->columns());
}
	

// Copies station estimates and variances from forward pass
// Closely matches functionality of CarryStnEstimatesandVariancesForward
// Used by:
//	- AdjustPhasedReverseCombine
// nextBlock = thisBlock-1
void dna_adjust::CarryStnEstimatesandVariancesCombine(
	const UINT32& nextBlock, const UINT32& thisBlock, 
	UINT32& pseudomsrJSLCount, bool MT_ReverseOrCombine)
{
	UINT32 pseudomsrJSLBegin(v_measurementParams_.at(thisBlock));

	//   intermediate                                   last
	//   _____________                              _____________
	//  |             |                            |             |
	//  | msr-comp    |                            | msr-comp    |
	//  |_____________|                            |_____________|  /__ pseudomsrJSLBegin
	//  |             |                            |             |  
	//  |  jrev-comp  |                            |  jfwd-comp  |
	//  |_____________|  /__ pseudomsrJSLBegin     |_____________|
	//  |             |  
	//  |  jfwd-comp  |
	//  |_____________|

	if (!v_blockMeta_.at(thisBlock)._blockLast)
		// PrepareAdjustmentCombine will have already weeded out
		// irrelevant blocks by testing _blockIsolated and _blockFirst 
		// So, add sufficient elements for junction station elements from
		// reverse pass.
		pseudomsrJSLBegin += static_cast<UINT32>(v_JSL_.at(thisBlock).size() * 3);	// junctions from reverse

	UINT32 pseudoMsrCount(static_cast<UINT32>(v_JSL_.at(nextBlock).size()));			// junctions from forward
	UINT32 pseudoMsrElemCount(pseudoMsrCount * 3);

	matrix_2d* AtVinv(&v_AtVinv_.at(thisBlock));
	matrix_2d* normals(&v_normals_.at(thisBlock));
	matrix_2d* measMinusComp(&v_measMinusComp_.at(thisBlock));
	matrix_2d* estimatedStations(&v_estimatedStations_.at(thisBlock));

	if (MT_ReverseOrCombine)
	{
		// restore backup copy of normals
		v_normalsR_.at(thisBlock) = v_normalsRC_.at(thisBlock);
		normals = &v_normalsR_.at(thisBlock);
	
		AtVinv = &v_AtVinvR_.at(thisBlock);
		measMinusComp = &v_measMinusCompR_.at(thisBlock);
		estimatedStations = &v_estimatedStationsR_.at(thisBlock);
	}
	else
	{

		// restore normals for single thread mode
		// v_normalsR_.at(thisBlock) is the set of normals with constraints
		// formed during the reverse adjustment
		v_normals_.at(thisBlock) = v_normalsR_.at(thisBlock);

	}

	AtVinv->grow(0, pseudoMsrElemCount);
	measMinusComp->grow(pseudoMsrElemCount, 0);

	pseudomsrJSLCount = measMinusComp->rows() - v_measurementParams_.at(thisBlock);

	UINT32 jsl, jsl_cov1, jslvar, jslcovar, paramCount, paramCount2;
	it_vUINT32 _it_jsl(v_JSL_.at(nextBlock).begin()), _it_jsl_cov;

	paramCount = 0;

	// Copy the junction variances from nextBlock (from forward) to thisBlock
	for (_it_jsl=v_JSL_.at(nextBlock).begin(); 
		_it_jsl!=v_JSL_.at(nextBlock).end(); 
		++_it_jsl, paramCount+=3)
	{
		// get index of this JSL in the next block
		jsl = v_blockStationsMap_.at(thisBlock)[*_it_jsl] * 3;
		jslvar = static_cast<UINT32>(std::distance(v_JSL_.at(nextBlock).begin(), _it_jsl) * 3);

		// add variances for this JSL to normals
		normals->blockadd(jsl, jsl, 
			v_junctionVariancesFwd_.at(nextBlock), jslvar, jslvar, 3, 3);

		// copy variances from v_junctionVariancesFwd_ to v_AtVinv_
		AtVinv->copyelements(jsl, pseudomsrJSLBegin + paramCount,
			v_junctionVariancesFwd_.at(nextBlock), jslvar, jslvar, 3, 3);

		// msr-comp elements
		// Measured-computed for junction stations carried forward (from thisBlock-1)
		// add station X measurement
		measMinusComp->put(pseudomsrJSLBegin + paramCount, 0, 
			(v_junctionEstimatesFwd_.at(nextBlock).get(jslvar, 0) -
			estimatedStations->get(jsl, 0)));			// X (meas - computed)
		
		// add station Y measurement
		measMinusComp->put(pseudomsrJSLBegin + paramCount+1, 0, 
			(v_junctionEstimatesFwd_.at(nextBlock).get(jslvar+1, 0) -
			estimatedStations->get(jsl+1, 0)));		// Y (meas - computed)
		
		// add station Z measurement
		measMinusComp->put(pseudomsrJSLBegin + paramCount+2, 0, 
			(v_junctionEstimatesFwd_.at(nextBlock).get(jslvar+2, 0) -
			estimatedStations->get(jsl+2, 0)));		// Z (meas - computed)

		paramCount2 = 0;

		// copy all covariances for this JSL in the next block
		for (_it_jsl_cov=_it_jsl; 
			 _it_jsl_cov!=v_JSL_.at(nextBlock).end(); 
			 ++_it_jsl_cov, paramCount2+=3)
		{
			// don't re-copy variances (i.e. diagonals)
			if (*_it_jsl_cov == *_it_jsl)
				continue;
			
			// get index of all covariances for this JSL in the next block
			jslcovar = static_cast<UINT32>(std::distance(v_JSL_.at(nextBlock).begin(), _it_jsl_cov) * 3);
			jsl_cov1 = v_blockStationsMap_.at(thisBlock)[*_it_jsl_cov] * 3;
			
			// add covariances for this JSL to normals
			normals->blockadd(jsl_cov1, jsl, v_junctionVariancesFwd_.at(nextBlock), jslcovar, jslvar, 3, 3);
			normals->blockadd(jsl, jsl_cov1, v_junctionVariancesFwd_.at(nextBlock), jslvar, jslcovar, 3, 3);

			// copy covariances from v_junctionVariancesFwd_ to v_AtVinv_
			AtVinv->copyelements(jsl_cov1, pseudomsrJSLBegin + paramCount, 
				v_junctionVariancesFwd_.at(nextBlock), jslcovar, jslvar, 3, 3);

			AtVinv->copyelements(jsl, pseudomsrJSLBegin + paramCount + paramCount2, 
				v_junctionVariancesFwd_.at(nextBlock), jslvar, jslcovar, 3, 3);
		}
	}

	if (projectSettings_.g.verbose > 0)
	{
		adj_file << std::endl << "Reverse block " << thisBlock;
		if (projectSettings_.a.stage)
			adj_file << " phased " << std::endl;
		else
			adj_file << " staged " << std::endl;
		adj_file << std::scientific << std::setprecision(16) << 
			"Design " << v_design_.at(thisBlock) << std::endl <<
			"AtVinv " << v_AtVinv_.at(thisBlock) << std::endl <<
			"Normals " << v_normals_.at(thisBlock) << std::endl;
	}
}
	

bool dna_adjust::PrepareAdjustmentCombine(const UINT32 currentBlock, UINT32& pseudomsrJSLCount, bool MT_ReverseOrCombine)
{
	// No combination is required for:
	// - first or last blocks of a contiguous network
	// - isolated blocks
	// - phased adjustment (block 1) only
	if (!CombineRequired(currentBlock))
		return false;
	
	matrix_2d* estimatedStations(&v_estimatedStations_.at(currentBlock));
	
	if (projectSettings_.a.multi_thread)
		estimatedStations = &v_estimatedStationsR_.at(currentBlock);

	// Now, update normals taking contribution from the junction station estimates and variances
	// of the preceding block estimated in the forward pass
	// 
	// Grow AtVinv and measMinusComp to include JSLs as 'pseudo measurements' from forward adjustment.
	//
	// For all intermediate blocks, CarryStnEstimatesandVariancesReverse is called immediately after
	// an adjustment in the reverse direction, in which measMinusComp is grown to include msr-comp
	// values using junction estimates obtained from currentBlock+1.  Then CarryStnEstimatesandVariancesCombine
	// is called, measMinusComp is grown again to include msr-comp values using junction estimates 
	// obtained from currentBlock+1.  Hence intermediate blocks will be formed as follows:
	//	 _____________
	//	|			  |
	//	|	msr-comp  |
	//	|_____________|
	//	|			  |
	//	|  jrev-comp  |	 msr-comp values using junction estimates obtained from currentBlock+1 reverse pass
	//	|_____________|
	//	|			  |
	//	|  jfwd-comp  |  msr-comp values using junction estimates obtained from currentBlock-1 forward pass
	//	|_____________|
	//
	// For the last block in a contiguous net, CarryStnEstimatesandVariancesCombine is called but measMinusComp 
	// is grown only to include msr-comp values using junction estimates obtained from currentBlock+1.  Hence,
	// the last block will be formed as follows:
	//	 _____________
	//	|			  |
	//	|	msr-comp  |
	//	|_____________|
	//	|			  |
	//	|  jfwd-comp  |  msr-comp values using junction estimates obtained from currentBlock-1 forward pass
	//	|_____________|
	//
	// The first block does not need a combination adjustment.

	// Reset coordinates to originals
	*estimatedStations = v_originalStations_.at(currentBlock);

	// Carry junction station estimates from currentBlock-1 (obtained during forward
	// adjustment) to this block.
	CarryStnEstimatesandVariancesCombine(currentBlock-1, currentBlock, 
		pseudomsrJSLCount, MT_ReverseOrCombine);

	AddConstraintStationstoNormalsCombine(currentBlock, MT_ReverseOrCombine);
	
	return true;

}
	

void dna_adjust::debug_SolutionInformation(const UINT32& currentBlock)
{
	if (projectSettings_.g.verbose > 0)
	{
		debug_file << "Block " << currentBlock + 1;
		if (projectSettings_.a.adjust_mode != SimultaneousMode)
			debug_file << (forward_ ? " (Forward)" : " (Reverse)");
		debug_file << std::endl;
		debug_file << "Pre-adjustment Estimates " << std::fixed << std::setprecision(16) << v_estimatedStations_.at(currentBlock);

		debug_file << "Block " << currentBlock + 1;
		if (projectSettings_.a.adjust_mode != SimultaneousMode)
			debug_file << (forward_ ? " (Forward)" : " (Reverse)");
		debug_file << std::endl;
		debug_file << "Design " << std::fixed << std::setprecision(16) << v_design_.at(currentBlock);

		debug_file << "Block " << currentBlock + 1;
		if (projectSettings_.a.adjust_mode != SimultaneousMode)
			debug_file << (forward_ ? " (Forward)" : " (Reverse)");
		debug_file << std::endl;
		debug_file << "Measurements " << std::fixed << std::setprecision(16) << v_measMinusComp_.at(currentBlock);

		debug_file << "Block " << currentBlock + 1;
		if (projectSettings_.a.adjust_mode != SimultaneousMode)
			debug_file << (forward_ ? " (Forward)" : " (Reverse)");
		debug_file << std::endl;
		debug_file << "At * V-inv " << std::fixed << std::setprecision(16) << v_AtVinv_.at(currentBlock) << std::endl;

		debug_file << "Block " << currentBlock + 1;
		if (projectSettings_.a.adjust_mode != SimultaneousMode)
			debug_file << (forward_ ? " (Forward)" : " (Reverse)");
		debug_file << std::endl;
		debug_file << "Normals " << std::fixed << std::setprecision(16) << v_normals_.at(currentBlock) << std::endl;
		debug_file.flush();
	}
}
	

void dna_adjust::debug_BlockInformation(const UINT32& currentBlock, const std::string& adjustment_method)
{
	if (projectSettings_.g.verbose > 3)
	{
		it_vUINT32 _it_stn;
		for (UINT32 i=0; i<v_parameterStationList_.at(currentBlock).size(); ++i)
		{
			debug_file << bstBinaryRecords_.at(v_parameterStationList_.at(currentBlock).at(i)).stationName;
			if ((_it_stn = find(v_ISL_.at(currentBlock).begin(), v_ISL_.at(currentBlock).end(), v_parameterStationList_.at(currentBlock).at(i))) != v_ISL_.at(currentBlock).end())
				debug_file << "i ";
			else if ((_it_stn = find(v_JSL_.at(currentBlock).begin(), v_JSL_.at(currentBlock).end(), v_parameterStationList_.at(currentBlock).at(i))) != v_JSL_.at(currentBlock).end())
				debug_file << "j ";
			else 
				debug_file << " ";				
		}
		debug_file << std::endl;

		debug_file << "Block " << currentBlock + 1 << " " << adjustment_method << std::endl;
		debug_file << "Adjustment Estimates " << std::fixed << std::setprecision(16) << v_estimatedStations_.at(currentBlock);
	}
}
	

// Single thread reverse and combination adjustment
void dna_adjust::AdjustPhasedReverseCombine()
{
	forward_ = false;
	isCombining_ = false;
	UINT32 currentBlock(blockCount_ - 1);
	UINT32 pseudomsrJSLCount(0);

	for (UINT32 block=0; block<blockCount_; ++block, --currentBlock)
	{
		if (IsCancelled())
			break;

		SetcurrentBlock(currentBlock);
	
		// If currentBlock is a single block, then there is no need to perform a 
		// reverse adjustment (i.e. continue);
		// Otherwise, if currentBlock is the last block, then this is the 
		// commencement of a reverse run.  In this case, recompute meas-minus-comp,
		// update normals
		if (!PrepareAdjustmentReverse(currentBlock, false))
			continue;		

		if (projectSettings_.o._adj_stn_iteration)
		{
			adj_file << std::endl << std::left << "Adjusting block " << currentBlock+1;
			if (!v_blockMeta_.at(currentBlock)._blockFirst)
				adj_file << " (reverse, in isolation)... ";
			else
				// all single blocks and first block of a contiguous list
				adj_file << " (reverse, rigorous)... ";
		}

		// At this point, whether first iteration or not, if currentBlock is the last block, 
		// the normals will have been initialised.  For all subsequent blocks, the normals will contain
		// the contribution of junction station coordinates and variances from preceding blocks.
		// In either case, the block is ready for adjustment.  Junction station coordinates and 
		// variances are carried in the reverse direction below

		if (projectSettings_.g.verbose > 3)
			debug_file << "In isolation" << std::endl;

		// Backup normals prior to inversion for re-use in combination
		// adjustment... only if a combination is required
		BackupNormals(currentBlock, false);

		// Least Squares Solution
		SolveTry(true, currentBlock);

		if (projectSettings_.o._adj_stn_iteration)
			adj_file << " done." << std::endl;

		// Add corrections to estimates.
		// If --output-adj-iter-stat argument is supplied or currentBlock
		// is the first block, then compute geographic coordinates, recompute 
		// meas-minus-computed vector, then print statistics.
		// Output of stats also prints largest correction.
		// If --output-adj-iter-stn argument is supplied, print station coords
		UpdateEstimatesReverse(currentBlock, false);

		// Debug and diagnose (if required)
		debug_BlockInformation(currentBlock, "(reverse, in isolation)");
		
		// Carry the estimated junction station coordinates and variances to the next block and
		// combine.  CarryReverseEstimates will return false if currentBlock is the first block
		// or an isolated block.  If estimates can be carried, CarryReverseEstimates updates 
		// geographic coords (if req'd), recomputes meas-minus-comp vector and updates normals for
		// the next block.
		if (CarryReverseJunctions(currentBlock, currentBlock-1, false))
		{			
			//////////////////////////////////////////////////////////////////////////
			// Combination Adjustment 
			//
			// Perform combination on all blocks except the first and last as they will
			// be rigorous from reverse and forward adjustments respectively
			//
			// First, carry the junction station estimates and variances of the next block
			// (obtained during the forward pass).  These were copied to v_junctionEstimatesFwd_
			// and v_junctionVariances_ in CarryStnEstimatesandVariancesForward(..) during
			// the forward pass. 
			if (PrepareAdjustmentCombine(currentBlock, pseudomsrJSLCount, false))
			{
				isCombining_ = true;
				if (projectSettings_.g.verbose > 3)
					debug_file << "Rigorous" << std::endl;
			
				if (projectSettings_.o._adj_stn_iteration)
					if (!v_blockMeta_.at(currentBlock)._blockFirst)
						adj_file << std::endl << std::left << "Adjusting block " << currentBlock+1 << " (reverse, rigorous)... ";

				// Least Squares Solution
				SolveTry(true, currentBlock);

				if (projectSettings_.o._adj_stn_iteration)
					adj_file << " done." << std::endl;

				// shrink msr-comp and AtVinv by pseudomsrJSLCount from forward pass
				UpdateEstimatesCombine(currentBlock, pseudomsrJSLCount);

				// Debug and diagnose (if required)
				debug_BlockInformation(currentBlock, "(reverse, rigorous)");

				isCombining_ = false;

			}	// if (PrepareAdjustmentCombine(currentBlock, pseudomsrJSLCount, false))  ...  combination adjustment

			// shrink msr-comp and AtVinv by JSL count from reverse pass

		} // if (CarryReverseJunctions(currentBlock, currentBlock-1, false))

		if (projectSettings_.a.stage)
		{
			ComputeAdjMsrBlockOnIteration(currentBlock);
			chiSquaredStage_ += chiSquared_;
			measurementParams_ += v_measurementParams_.at(currentBlock);
		}

		// Update the rigorous coordinates
		UpdateEstimatesFinal(currentBlock);

		// For staged adjustments, write to disk and unload matrix data
		if (projectSettings_.a.stage)
			OffloadBlockToMappedFile(currentBlock);

	}	// for (UINT32 block=0; block<blockCount_; ++block, --currentBlock)
}
	

// Single thread reverse and combination adjustment
void dna_adjust::AdjustPhasedReverse()
{
	forward_ = false;
	isCombining_ = false;
	UINT32 currentBlock(blockCount_ - 1);

	// For staged adjustments, load the last block from mapped memory file
	if (projectSettings_.a.stage)
		DeserialiseBlockFromMappedFile(currentBlock);

	for (UINT32 block=0; block<blockCount_; ++block, --currentBlock)
	{
		if (IsCancelled())
			break;

		SetcurrentBlock(currentBlock);

		// If currentBlock is a single block, then there is no need to perform a 
		// reverse adjustment (i.e. continue);
		// Otherwise, if currentBlock is the last block, then this is the 
		// commencement of a reverse run.  In this case, recompute meas-minus-comp,
		// update normals
		if (!PrepareAdjustmentReverse(currentBlock, false))
			continue;		

		if (projectSettings_.o._adj_stn_iteration)
		{
			adj_file << std::endl << std::left << "Adjusting block " << currentBlock+1;
			if (!v_blockMeta_.at(currentBlock)._blockFirst)
				adj_file << " (reverse, in isolation)... ";
			else
				// all single blocks and first block of a contiguous list
				adj_file << " (reverse, rigorous)... ";
		}

		// At this point, whether first iteration or not, if currentBlock is the last block, 
		// the normals will have been initialised.  For all subsequent blocks, the normals will contain
		// the contribution of junction station coordinates and variances from preceding blocks.
		// In either case, the block is ready for adjustment.  Junction station coordinates and 
		// variances are carried in the reverse direction below
		
		if (projectSettings_.g.verbose > 3)
			debug_file << "In isolation" << std::endl;

		// Least Squares Solution
		SolveTry(true, currentBlock);

		if (projectSettings_.o._adj_stn_iteration)
			adj_file << " done." << std::endl;
		
		// Add corrections to estimates.
		// If --output-adj-iter-stat argument is supplied or currentBlock
		// is the first block, then compute geographic coordinates, recompute 
		// meas-minus-computed vector, then print statistics.
		// Output of stats also prints largest correction.
		// If --output-adj-iter-stn argument is supplied, print station coords
		UpdateEstimatesReverse(currentBlock, false);

		// Debug (if required)
		debug_BlockInformation(currentBlock, "(reverse, in isolation)");

		// Carry the estimated junction station coordinates and variances to the next block and
		// combine.  CarryReverseEstimates will return false if currentBlock is the first block
		// or an isolated block.  If estimates can be carried, CarryReverseEstimates updates 
		// geographic coords (if req'd), recomputes meas-minus-comp vector and updates normals for
		// the next block.
		CarryReverseJunctions(currentBlock, currentBlock-1, false);

		// Update the rigorous coordinates
		UpdateEstimatesFinal(currentBlock);

		// For staged adjustments, write to disk and unload matrix data
		if (projectSettings_.a.stage)
			OffloadBlockToMappedFile(currentBlock);

	}	// for (UINT32 block=0; block<blockCount_; ++block, --currentBlock)
}
	

// This method will not be reached if currentBlock is a single (isolated) block, since 
// PrepareAdjustmentReverse performs this test.
void dna_adjust::UpdateEstimatesReverse(const UINT32 currentBlock, bool MT_ReverseOrCombine)
{
	matrix_2d* estimatedStations(&v_estimatedStations_.at(currentBlock));
	matrix_2d* corrections(&v_corrections_.at(currentBlock));
	matrix_2d* aposterioriVariances(&v_normals_.at(currentBlock));

	if (MT_ReverseOrCombine)
	{
		estimatedStations = &v_estimatedStationsR_.at(currentBlock);
		corrections = &v_correctionsR_.at(currentBlock);
		aposterioriVariances = &v_normalsR_.at(currentBlock);
	}	

	// update station coordinates with lsq-estimated corrections
	estimatedStations->add(*corrections);

	// Print the estimated stations for this block
	// At this point, all station estimates are "in isolation".  Rigorous to follow...
	if (projectSettings_.o._adj_stn_iteration)
	{
		adj_file_mutex.lock();

		if (projectSettings_.a.multi_thread)
		{
			adj_file << std::endl << "2> Adjusted block " << currentBlock + 1;
			if (v_blockMeta_.at(currentBlock)._blockFirst)
				adj_file << " (reverse, rigorous) " << std::endl;
			else
				adj_file << " (reverse, in isolation) " << std::endl;	
		}
		
		printer_->PrintBlockStations(adj_file, currentBlock, 
			estimatedStations, aposterioriVariances, 
			false, true, false, true, false);
		adj_file_mutex.unlock();
	}
	
}
	

void dna_adjust::UpdateEstimatesCombine(const UINT32 currentBlock, UINT32 pseudomsrJSLCount)
{
	matrix_2d* estimatedStations(&v_estimatedStations_.at(currentBlock));
	matrix_2d* corrections(&v_corrections_.at(currentBlock));
	matrix_2d* AtVinv(&v_AtVinv_.at(currentBlock));
	matrix_2d* measMinusComp(&v_measMinusComp_.at(currentBlock));

	if (projectSettings_.a.multi_thread)
	{
		estimatedStations = &v_estimatedStationsR_.at(currentBlock);
		corrections = &v_correctionsR_.at(currentBlock);
		AtVinv = &v_AtVinvR_.at(currentBlock);
		measMinusComp = &v_measMinusCompR_.at(currentBlock);
	}	

	// copy estimated parameter station coordinates
	estimatedStations->add(*corrections);

	// Now, combination adjustment is complete, so shrink to "exclude" JSLs from 
	// the next reverse adjustment.
	AtVinv->shrink(0, pseudomsrJSLCount);
	measMinusComp->shrink(pseudomsrJSLCount, 0);
}
	

// Updates rigorous estimates after a reverse/combine adjustment
void dna_adjust::UpdateEstimatesFinal(const UINT32 currentBlock)
{
	// Is this the last block?  If so, the rigorous estimates were updated during the 
	// forward pass - copy values from v_correctionsR_.
	if (v_blockMeta_.at(currentBlock)._blockLast)
	{
		if (projectSettings_.a.multi_thread)
			return;
		if (projectSettings_.a.adjust_mode == Phased_Block_1Mode)
			return;
		v_corrections_.at(currentBlock) = v_correctionsR_.at(currentBlock);
		return;
	}
	
	matrix_2d* estimatedStations(&v_estimatedStations_.at(currentBlock));
	matrix_2d* corrections(&v_corrections_.at(currentBlock));
	matrix_2d* AtVinv(&v_AtVinv_.at(currentBlock));
	matrix_2d* measMinusComp(&v_measMinusComp_.at(currentBlock));
	matrix_2d* aposterioriVariances(&v_normals_.at(currentBlock));

	if (projectSettings_.a.multi_thread)
	{
		estimatedStations = &v_estimatedStationsR_.at(currentBlock);
		corrections = &v_correctionsR_.at(currentBlock);
		AtVinv = &v_AtVinvR_.at(currentBlock);
		measMinusComp = &v_measMinusCompR_.at(currentBlock);
		aposterioriVariances = &v_normalsR_.at(currentBlock);
	}	

	if (v_blockMeta_.at(currentBlock)._blockFirst)
	{
		// If this point is reached, rigorous adjustment is complete for a contiguous network, so
		// shrink to "exclude" JSLs from the next iteration.
		AtVinv->shrink(0, static_cast<UINT32>(v_JSL_.at(currentBlock).size() * 3));
		measMinusComp->shrink(static_cast<UINT32>(v_JSL_.at(currentBlock).size() * 3), 0);
	}
	
	// update max correction
	if (fabs(corrections->compute_maximum_value()) > fabs(maxCorr_))
		SetmaxCorr(corrections->maxvalue());

	// update largest correction
	if (fabs(corrections->maxvalue()) > fabs(largestCorr_))
	{
		largestCorr_ = corrections->maxvalue();
		blockLargeCorr_ = currentBlock;
	}

	if (projectSettings_.a.multi_thread)
	{
		v_corrections_.at(currentBlock) = v_correctionsR_.at(currentBlock);
		v_estimatedStations_.at(currentBlock) = v_estimatedStationsR_.at(currentBlock);
		v_normals_.at(currentBlock) = v_normalsR_.at(currentBlock);
	}

	// Now copy 'estimated' coordinates to 'rigorous' for comparison on the next iteration
	v_rigorousStations_.at(currentBlock) = *estimatedStations;
	v_rigorousVariances_.at(currentBlock) = *aposterioriVariances;
	
	// update original coordinates
	v_originalStations_.at(currentBlock) = v_rigorousStations_.at(currentBlock);

	// If this is the first block, return since UpdateEstimateReverse prints station
	// coordinates
	if (v_blockMeta_.at(currentBlock)._blockFirst)
		return;

	// print the 'rigorous' stations for this block
	if (projectSettings_.o._adj_stn_iteration)
	{

		adj_file_mutex.lock();
		if (projectSettings_.a.multi_thread)
			adj_file << std::endl << "3> Adjusted block " << currentBlock + 1 << " (rigorous)" << std::endl;

		printer_->PrintBlockStations(adj_file, currentBlock, 
			&v_rigorousStations_.at(currentBlock), &v_rigorousVariances_.at(currentBlock), 
			false, true, false, true, false);		// update coordinates

		adj_file_mutex.unlock();
	}
}
	

// nextBlock = currentBlock-1
bool dna_adjust::CarryReverseJunctions(const UINT32 currentBlock, const UINT32 nextBlock, bool MT_ReverseOrCombine)
{
	// Carrying of junction parameters not required for single blocks, nor 
	// if this is the first block
	if (v_blockMeta_.at(currentBlock)._blockIsolated)
		return false;
	if (v_blockMeta_.at(currentBlock)._blockFirst)
		return false;

	// matrices from previous block
	matrix_2d* estimatedStationsNext(&v_estimatedStations_.at(nextBlock));
	
	if (MT_ReverseOrCombine)
		estimatedStationsNext = &v_estimatedStationsR_.at(nextBlock);
	else
	{
		// Restore back up copy of normals for reverse adjustment in 
		// single thread mode
		if (!projectSettings_.a.stage)
			v_normals_.at(nextBlock) = v_normalsR_.at(nextBlock);
	}

	// For staged adjustments, load blocks info for nextBlock
	if (projectSettings_.a.stage)
		DeserialiseBlockFromMappedFile(nextBlock);

	// Next, update estimates for next block to original estimates
	*estimatedStationsNext = v_originalStations_.at(nextBlock);
	
	// For staged adjustments, rebuild design and AtVinv
	if (projectSettings_.a.stage)
		RebuildNormals(nextBlock, __reverse__, false, false);
	
	// Now, carry the estimated junction station coordinates and variances 
	// to the next block (to be used in the next loop).
	// Remember - the junction station estimates and variances of the next block
	// (obtained during the forward pass) were copied during the forward pass
	// in CarryStnEstimatesandVariancesForward(..), so no need to back-up prior
	// to carrying in reverse.

	// copy estimated coordinates and variances for the junction stations from
	// this block to the next block (which is nextBlock)
	CarryStnEstimatesandVariancesReverse(nextBlock, currentBlock, MT_ReverseOrCombine);
	
	AddConstraintStationstoNormalsReverse(nextBlock, MT_ReverseOrCombine);

	return true;
}



// go through each of the measurements in the binary measurements file and formulate partial derivatives
void dna_adjust::FillDesignNormalMeasurementsMatrices(bool buildnewMatrices, const UINT32& block, bool MT_ReverseOrCombine)
{
	UINT32 design_row(0);
	
	it_vUINT32 _it_block_msr;
	it_vmsr_t _it_msr; 

	for (_it_block_msr=v_CML_.at(block).begin(); _it_block_msr!=v_CML_.at(block).end(); ++_it_block_msr)
	{
		if (InitialiseandValidateMsrPointer(_it_block_msr, _it_msr))
			continue;

		// When a target direction is found, continue to next element.  
		if (_it_msr->measType == 'D')
			if (_it_msr->vectorCount2 < 1)
				continue;

		// Build AtVinv, Normals and Meas minus Comp vectors
		UpdateDesignNormalMeasMatrices(&_it_msr, design_row, buildnewMatrices, block, MT_ReverseOrCombine);
	}
}

// Initialise msr pointer if new matrices need to be built
// If the data has already been reduced, then restore via the
// pre-adjustment value, otherwise back up the current value.
bool dna_adjust::InitialiseMeasurement(pit_vmsr_t _it_msr, bool buildnewMatrices)
{
	if (buildnewMatrices)
	{
		// Was this measurement reduced from previous
		// adjustment, which was serialised to disk?
		if (bms_meta_.reduced)
			// Restore measured value
			(*_it_msr)->term1 = (*_it_msr)->preAdjMeas;
		else
			// This is the first time adjust has run, so
			// back up the raw measurement
			(*_it_msr)->preAdjMeas = (*_it_msr)->term1;

		if (projectSettings_.a.stage)
			return true;
	}

	return false;
}
	

void dna_adjust::UpdateDesignNormalMeasMatrices(pit_vmsr_t _it_msr, UINT32& design_row, bool buildnewMatrices, const UINT32& block, bool MT_ReverseOrCombine)
{
	std::stringstream ss;

	matrix_2d* estimatedStations(&v_estimatedStations_.at(block));
	matrix_2d* design(&v_design_.at(block));
	matrix_2d* AtVinv(&v_AtVinv_.at(block));
	matrix_2d* measMinusComp(&v_measMinusComp_.at(block));
	matrix_2d* normals(&v_normals_.at(block));

	if (MT_ReverseOrCombine)
	{
		estimatedStations = &v_estimatedStationsR_.at(block);
		design = &v_designR_.at(block);
		AtVinv = &v_AtVinvR_.at(block);
		measMinusComp = &v_measMinusCompR_.at(block);
	}	

	switch ((*_it_msr)->measType)
	{
	case 'A':	// Horizontal angle
		UpdateDesignNormalMeasMatrices_A(_it_msr, design_row, block,
			measMinusComp, estimatedStations, normals, design, AtVinv, buildnewMatrices);
		break;
	case 'B':	// Geodetic azimuth
	case 'K':	// Astronomic azimuth
		// Note: UpdateDesignNormalMeasMatrices_BK reduces K measurements to the geodetic reference frame, 
		// after which UpdateDesignNormalMeasMatrices_BK treats K measurements as B measurements upon forming
		// design matrix elements.
		UpdateDesignNormalMeasMatrices_BK(_it_msr, design_row, block,
			measMinusComp, estimatedStations, normals, design, AtVinv, buildnewMatrices);	
		break;
	case 'C':	// Chord dist
		UpdateDesignNormalMeasMatrices_C(_it_msr, design_row, block,
			measMinusComp, estimatedStations, normals, design, AtVinv, buildnewMatrices);
		break;
	case 'D':	// Direction set	
		UpdateDesignNormalMeasMatrices_D(_it_msr, design_row, block,
			measMinusComp, estimatedStations, normals, design, AtVinv, buildnewMatrices);
		break;
	case 'E':	// Ellipsoid arc
		UpdateDesignNormalMeasMatrices_E(_it_msr, design_row, block,
			measMinusComp, estimatedStations, normals, design, AtVinv, buildnewMatrices);
		break;
	case 'G':	// GPS Baseline
		UpdateDesignNormalMeasMatrices_G(_it_msr, design_row, block,
			measMinusComp, estimatedStations, normals, design, AtVinv, buildnewMatrices);
		break;
	case 'H':	// Orthometric height
		// Note: UpdateDesignNormalMeasMatrices_H reduces term1 to ellipsoid height, after which
		// UpdateDesignNormalMeasMatrices_HR is used to form design elements.
		UpdateDesignNormalMeasMatrices_H(_it_msr, design_row, block,
			measMinusComp, estimatedStations, normals, design, AtVinv, buildnewMatrices);
		break;
	case 'I':	// Astronomic latitude
		// Note: UpdateDesignNormalMeasMatrices_I reduces term1 to geodetic latitude, after which
		// UpdateDesignNormalMeasMatrices_IP is used to form design elements.
		UpdateDesignNormalMeasMatrices_I(_it_msr, design_row, block,
			measMinusComp, estimatedStations, normals, design, AtVinv, buildnewMatrices);
		break;
	case 'J':	// Astronomic longitude
		// Note: UpdateDesignNormalMeasMatrices_J reduces term1 to geodetic longitude, after which
		// UpdateDesignNormalMeasMatrices_JQ is used to form design elements.
		UpdateDesignNormalMeasMatrices_J(_it_msr, design_row, block,
			measMinusComp, estimatedStations, normals, design, AtVinv, buildnewMatrices);
		break;		
	case 'L':	// Level difference
		UpdateDesignNormalMeasMatrices_L(_it_msr, design_row, block,
			measMinusComp, estimatedStations, normals, design, AtVinv, buildnewMatrices);
		break;		
	case 'M':	// MSL arc
		UpdateDesignNormalMeasMatrices_M(_it_msr, design_row, block,
			measMinusComp, estimatedStations, normals, design, AtVinv, buildnewMatrices);
		break;		
	case 'P':	// Geodetic latitude
		// Note: UpdateDesignNormalMeasMatrices_P archives the raw measurement, after which
		// UpdateDesignNormalMeasMatrices_IP is used to form design elements.
		UpdateDesignNormalMeasMatrices_P(_it_msr, design_row, block,
			measMinusComp, estimatedStations, normals, design, AtVinv, buildnewMatrices);
		break;		
	case 'Q':	// Geodetic longitude
		// Note: UpdateDesignNormalMeasMatrices_Q archives the raw measurement, after which
		// UpdateDesignNormalMeasMatrices_JQ is used to form design elements.
		UpdateDesignNormalMeasMatrices_Q(_it_msr, design_row, block,
			measMinusComp, estimatedStations, normals, design, AtVinv, buildnewMatrices);
		break;		
	case 'R':	// Ellipsoidal height
		// Note: UpdateDesignNormalMeasMatrices_R archives the raw measurement, after which
		// UpdateDesignNormalMeasMatrices_HR is used to form design elements.
		UpdateDesignNormalMeasMatrices_R(_it_msr, design_row, block,
			measMinusComp, estimatedStations, normals, design, AtVinv, buildnewMatrices);
		break;
	case 'S':	// Slope distance
		UpdateDesignNormalMeasMatrices_S(_it_msr, design_row, block,
			measMinusComp, estimatedStations, normals, design, AtVinv, buildnewMatrices);
		break;
	case 'V':	// Zenith distance
		UpdateDesignNormalMeasMatrices_V(_it_msr, design_row, block,
			measMinusComp, estimatedStations, normals, design, AtVinv, buildnewMatrices);
		break;		
	case 'X':	// GPS Baseline cluster
		UpdateDesignNormalMeasMatrices_X(_it_msr, design_row, block,
			measMinusComp, estimatedStations, design, AtVinv, buildnewMatrices);
		break;		
	case 'Y':	// GPS Point cluster
		UpdateDesignNormalMeasMatrices_Y(_it_msr, design_row, block,
			measMinusComp, estimatedStations, design, AtVinv, buildnewMatrices);
		break;		
	case 'Z':	// Vertical angle
		UpdateDesignNormalMeasMatrices_Z(_it_msr, design_row, block,
			measMinusComp, estimatedStations, normals, design, AtVinv, buildnewMatrices);
		break;
	default:
		ss << "UpdateDesignNormalMeasMatrices(): Unknown measurement type - '" <<
			(*_it_msr)->measType << "'." << std::endl;
		SignalExceptionAdjustment(ss.str(), block);
	}
}

// Forms a variance matrix representing correlated angles formed from
// a round of directions
void dna_adjust::LoadVarianceMatrix_D(it_vmsr_t _it_msr, matrix_2d* var_dirn, bool buildnewMatrices)
{
	// load the correlated angles variance matrix from the binary file
	if (bms_meta_.reduced || !buildnewMatrices)
	{
		try {
			// load as-is, assuming that an adjustment has been run
			// and the variances have been reduced (scaled)
			GetDirectionsVarianceMatrix(_it_msr, var_dirn);
			FormInverseVarianceMatrix(var_dirn, true);
			return;
		}
		catch (const std::runtime_error& e) {

			// Print error message to adj file and throw exception
			std::stringstream ss;
			printer_->PrintMsrVarianceMatrixException(_it_msr, e, ss, "LoadVarianceMatrix_D", 0);
			SignalExceptionAdjustment(ss.str(), 0);
		}
	}

	// OK, create the correlated angles variance matrix from the binary 
	// file

	it_vmsr_t _it_msr_first(_it_msr);
	
	// angle, variance and covariance set in the binary records as follows:
	// term1 = direction
	// term2 = variance (direction)
	// term3 = instrument height (leave as-is)
	// term4 = target height (leave as-is)
	// scale1 = derived angle corrected for deflection of the vertical
	// scale2 = variance (angle)
	// scale3 = covariance (angle)
	// vectorCount2 = number of non-ignored directions
	// preAdjMeas = original derived angle

	UINT32 a, angle_count(_it_msr->vectorCount2 - 1);		// number of directions excluding the RO
	
	matrix_2d A(angle_count, _it_msr->vectorCount2);
	matrix_2d AV(angle_count, _it_msr->vectorCount2);
	
	var_dirn->redim(angle_count, angle_count);
	var_dirn->zero();
	
	double previousVariance(_it_msr->term2);
	UINT32 skip(0), ignored(_it_msr->vectorCount1 - _it_msr->vectorCount2);

	if (projectSettings_.g.verbose > 5)
		debug_file << std::endl << "Std dev " << std::scientific << std::setprecision(16) << _it_msr->term2 << " (" << std::fixed << std::setprecision(2) << Seconds(sqrt(_it_msr->term2)) << " seconds)" << std::endl;
	
	_it_msr++;

	// propagate (uncorrelated) directions variances to angles variances
	try {
		for (a=0; a<angle_count; ++a)
		{
			// cater for ignored directions
			if (_it_msr->ignore)
			{
				while (skip < ignored)
				{
					skip++;
					_it_msr++;
					if (!_it_msr->ignore)
						break;
				}
			}

			// Fill design & variance matrices to propagate variances from directions to angles
			A.put(a, a, -1);
			A.put(a, a+1, 1);
			AV.put(a, a, previousVariance * -1);
			AV.put(a, a+1, _it_msr->term2);

			if (projectSettings_.g.verbose > 5)
				debug_file << "Std dev " << std::scientific << std::setprecision(16) << _it_msr->term2 << " (" << std::fixed << std::setprecision(2) << Seconds(sqrt(_it_msr->term2)) << " seconds)" << std::endl;

			if (a+1 == angle_count)
				break;

			previousVariance = _it_msr->term2;
			_it_msr++;
		}

		if (projectSettings_.g.verbose > 6)
		{
			debug_file << std::endl << "Directions variance matrix:" << std::endl;
			debug_file << "A " << A << std::endl;
			debug_file << "AV " << std::scientific << std::setprecision(16) << AV << std::endl;
		}

		// AVAT for correlated angles will not be a fully populated matrix, but a 
		// banded one whereby correlations exist only for angles adjacent to each other.  
		// So, propagate manually.
		UINT32 d;
		for (a=0; a<angle_count; ++a)
		{
			// variances
			for (d=0; d<2; ++d)
				var_dirn->elementadd(a, a, AV.get(a, a+d) * A.get(a, a+d));
			if (a+1 == angle_count)
				break;
			// covariances
			var_dirn->elementadd(a, a+1, AV.get(a, a+1) * A.get(a+1, a+1));
			var_dirn->elementadd(a+1, a, var_dirn->get(a, a+1));
		}

		// Now the variance matrix has been formed, set the variances and 
		// covariances to the binary file
		SetDirectionsVarianceMatrix(_it_msr_first, *var_dirn);

		if (projectSettings_.g.verbose > 5)
			debug_file << "V.dxyz " << std::scientific << std::setprecision(16) << std::setw(26) << *var_dirn;
		
		// Form inverse
		FormInverseVarianceMatrix(var_dirn);
	}
	catch (const std::runtime_error& e) {

		// Print error message to adj file and throw exception
		std::stringstream ss;
		printer_->PrintMsrVarianceMatrixException(_it_msr_first, e, ss, "LoadVarianceMatrix_D", 0);
		SignalExceptionAdjustment(ss.str(), 0);
	}

	if (projectSettings_.g.verbose > 5)
		debug_file << "Inv V.dxyz " << std::scientific << std::setprecision(16) << std::setw(26) << *var_dirn;
	
}

bool dna_adjust::FormInverseVarianceMatrixReduced(it_vmsr_t _it_msr, matrix_2d* var_cart,
	const std::string& method_name)
{
	// Has the variance matrix already been reduced (i.e. propagated and scaled)
	if (bms_meta_.reduced)
	{
		try {
			// load as-is, assuming that an adjustment has been run
			// and the variances have been reduced (propagated and scaled)
			FormInverseGPSVarianceMatrix(_it_msr, var_cart);
			return true;
		}
		catch (const std::runtime_error& e) {

			// Print error message to adj file and throw exception
			std::stringstream ss;
			printer_->PrintMsrVarianceMatrixException(_it_msr, e, ss, method_name, 0);
			SignalExceptionAdjustment(ss.str(), 0);
		}
	}
	return false;
}
	
// load the GPS variance matrix from the binary file
void dna_adjust::LoadVarianceMatrix_G(it_vmsr_t _it_msr, matrix_2d* var_cart)
{
	if (FormInverseVarianceMatrixReduced(_it_msr, var_cart, "LoadVarianceMatrix_G"))
		return;

	// OK, load the GPS variance matrix from the binary 
	// file, scaling where necessary
	
	it_vmsr_t _it_msr_first(_it_msr);
	UINT32 v(0);
	
	// Determine scaling
	double vScale, pScale, lScale, hScale;
	bool scaleMatrix, scalePartial;
	LoadVarianceScaling(_it_msr, vScale, pScale, lScale, hScale,
		scaleMatrix, scalePartial);

	var_cart->redim(3, 3);

	// vScale is performed on the fly, assuming that all VCVs for single baselines
	// will always be in the cartesian reference frame

	var_cart->put(0, v, scaleMatrix ? _it_msr->term2 * vScale : _it_msr->term2);	// XX

	_it_msr++;
	v++;

	var_cart->put(0, v, scaleMatrix ? _it_msr->term2 * vScale : _it_msr->term2);	// YX
	var_cart->put(1, v, scaleMatrix ? _it_msr->term3 * vScale : _it_msr->term3);	// YY
	
	_it_msr++;
	v++;

	var_cart->put(0, v, scaleMatrix ? _it_msr->term2 * vScale : _it_msr->term2);		// ZX
	var_cart->put(1, v, scaleMatrix ? _it_msr->term3 * vScale : _it_msr->term3);		// ZY
	var_cart->put(2, v, scaleMatrix ? _it_msr->term4 * vScale : _it_msr->term4);		// ZZ
	
	bool lowerisClear(true);
	if (scaleMatrix || scalePartial)
	{
		var_cart->filllower();
		lowerisClear = false;
	}

	try {

		// Scale variance matrix using phi, lambda and height scalars?
		if (scalePartial)
		{
			it_vstn_t_const stn1_it(bstBinaryRecords_.begin() + _it_msr->station1);
			it_vstn_t_const stn2_it(bstBinaryRecords_.begin() + _it_msr->station2);

			// Note - it makes little difference whether the start point, end point or mid point
			// is chosen to form the rotation matrix.  Tests on a 1987.6 Km baseline show
			// sub-millimetre difference in the propagated results!
			// So no need to average - just use the starting point

			ScaleGPSVCV<double>(*var_cart, var_cart,
				stn1_it->currentLatitude, stn1_it->currentLongitude, stn1_it->currentHeight,
				datum_.GetEllipsoidRef(), 
				pScale, lScale, hScale);
		}

		// Copy elements of scaled variance matrix to variances held in internal memory (_it_msr)
		// only if VCV has been scaled (by vScale, or by pScale+lScale+hScale).
		// This alleviates having to repeat any scaling for subsequent computations
		// involving variances (i.e. sigma zero, Pelzer's reliability, etc)
		if (scaleMatrix || scalePartial)
			SetGPSVarianceMatrix(_it_msr_first, *var_cart);
		
		if (projectSettings_.a.stage)
			return;
		
		// inverse
		FormInverseVarianceMatrix(var_cart, lowerisClear);
		
		if (boost::math::isnan(var_cart->get(0, 0)) || boost::math::isinf(var_cart->get(0, 0)))
		{
			std::stringstream ss;
			ss << "Invalid variance matrix:" << std::endl;
			ss << std::setprecision(6) << std::fixed << *var_cart;
			SignalExceptionAdjustment(ss.str(), 0);
		}
	}
	catch (const std::runtime_error& e) {

		// Print error message to adj file and throw exception
		std::stringstream ss;
		printer_->PrintMsrVarianceMatrixException(_it_msr_first, e, ss, "LoadVarianceMatrix_G", 0);
		SignalExceptionAdjustment(ss.str(), 0);
	}

	if (projectSettings_.g.verbose > 5)
		debug_file << std::endl << "Inv V.dxyz " << std::fixed << std::setprecision(16) << std::setw(26) << *var_cart;
	
}
	
// load the GPS variance matrix from the binary file
void dna_adjust::LoadVarianceMatrix_X(it_vmsr_t _it_msr, matrix_2d* var_cart)
{
	if (FormInverseVarianceMatrixReduced(_it_msr, var_cart, "LoadVarianceMatrix_X"))
		return;

	// OK, load the GPS variance matrix from the binary 
	// file, scaling where necessary
	
	it_vmsr_t _it_msr_first(_it_msr);

	UINT32 covr(0), covc(0);
	UINT32 cluster_bsl, baseline_count(_it_msr->vectorCount1);
	UINT32 cluster_cov, covariance_count;
	var_cart->redim(baseline_count * 3, baseline_count * 3);

	matrix_2d mpositions(baseline_count * 3, 1);

	it_vstn_t stn1_it;
	it_vstn_t stn2_it;

	// Determine scaling
	double vScale, pScale, lScale, hScale;
	bool scaleMatrix, scalePartial;
	LoadVarianceScaling(_it_msr, vScale, pScale, lScale, hScale,
		scaleMatrix, scalePartial);

	for (cluster_bsl=0; cluster_bsl<baseline_count; ++cluster_bsl)	// number of baselines/points
	{
		covr = cluster_bsl * 3;
		covariance_count = _it_msr->vectorCount2;

		stn1_it = bstBinaryRecords_.begin() + _it_msr->station1;
		stn2_it = bstBinaryRecords_.begin() + _it_msr->station2;
	
		if (scalePartial)
		{
			// Note - it makes little difference whether the start point, end point or mid point
			// is chosen to form the rotation matrix.  Tests on a 1987.6 Km baseline show
			// sub-millimetre difference in the propagated results!
			// So no need to average - just use the starting point
			mpositions.put(covr, 0, stn1_it->currentLatitude);
			mpositions.put(covr+1, 0, stn1_it->currentLongitude);
			mpositions.put(covr+2, 0, stn1_it->currentHeight);
		}

		// vScale is performed on the fly, assuming that all VCVs for single baselines
		// will always be in the cartesian reference frame

		var_cart->put(covr, covr, scaleMatrix ? _it_msr->term2 * vScale : _it_msr->term2);	// XX
		
		_it_msr++;

		var_cart->put(covr,   covr+1, scaleMatrix ? _it_msr->term2 * vScale : _it_msr->term2);	// YX
		var_cart->put(covr+1, covr+1, scaleMatrix ? _it_msr->term3 * vScale : _it_msr->term3);	// YY
		
		_it_msr++;

		var_cart->put(covr,   covr+2, scaleMatrix ? _it_msr->term2 * vScale : _it_msr->term2);	// ZX
		var_cart->put(covr+1, covr+2, scaleMatrix ? _it_msr->term3 * vScale : _it_msr->term3);	// ZY
		var_cart->put(covr+2, covr+2, scaleMatrix ? _it_msr->term4 * vScale : _it_msr->term4);	// ZZ

		if (covariance_count == 0)
			break;

		covc = cluster_bsl * 3 + 3;

		for (cluster_cov=0; cluster_cov<covariance_count; ++cluster_cov)	// number of baseline/point covariances
		{
			covr = cluster_bsl * 3;
			_it_msr++;		// X, lat
			var_cart->put(covr, covc, scaleMatrix ? _it_msr->term1 * vScale : _it_msr->term1);		// m11
			var_cart->put(covr, covc+1, scaleMatrix ? _it_msr->term2 * vScale : _it_msr->term2);	// m12
			var_cart->put(covr, covc+2, scaleMatrix ? _it_msr->term3 * vScale : _it_msr->term3);	// m13

			covr++;
			_it_msr++;		// Y, long
			var_cart->put(covr, covc, scaleMatrix ? _it_msr->term1 * vScale : _it_msr->term1);		// m21
			var_cart->put(covr, covc+1, scaleMatrix ? _it_msr->term2 * vScale : _it_msr->term2);	// m22
			var_cart->put(covr, covc+2, scaleMatrix ? _it_msr->term3 * vScale : _it_msr->term3);	// m23

			covr++;
			_it_msr++;		// Z, height
			var_cart->put(covr, covc, scaleMatrix ? _it_msr->term1 * vScale : _it_msr->term1);		// m31
			var_cart->put(covr, covc+1, scaleMatrix ? _it_msr->term2 * vScale : _it_msr->term2);	// m32
			var_cart->put(covr, covc+2, scaleMatrix ? _it_msr->term3 * vScale : _it_msr->term3);	// m33			

			covc+=3;
		}

		// next cluster measurement
		_it_msr++;
	}
	
	bool lowerisClear(true);

	try {

		if (scaleMatrix || scalePartial)
		{	
			if (scalePartial)
			{
				var_cart->filllower();
				lowerisClear = false;
				ScaleGPSVCV_Cluster<double>(*var_cart, var_cart, 
					mpositions, 
					datum_.GetEllipsoidRef(), 
					pScale, lScale, hScale);

			}	// if (scalePartial)

			// Copy elements of scaled variance matrix to variances held in internal memory (_it_msr)
			// This alleviates having to repeat any scaling for subsequent computations
			// involving variances (i.e. sigma zero, Pelzer's reliability, etc)
			SetGPSVarianceMatrix<it_vmsr_t>(_it_msr_first, *var_cart);

		}	// if (scaleMatrix || scalePartial)
		// inverse
		FormInverseVarianceMatrix(var_cart, lowerisClear);

		if (boost::math::isnan(var_cart->get(0, 0)) || boost::math::isinf(var_cart->get(0, 0)))
		{
			std::stringstream ss;
			ss << "Invalid variance matrix:" << std::endl;
			ss << std::setprecision(6) << std::fixed << *var_cart;
			SignalExceptionAdjustment(ss.str(), 0);
		}
	}
	catch (const std::runtime_error& e) {

		// Print error message to adj file and throw exception
		std::stringstream ss;
		printer_->PrintMsrVarianceMatrixException(_it_msr_first, e, ss, "LoadVarianceMatrix_X", baseline_count);
		SignalExceptionAdjustment(ss.str(), 0);
	}

	if (projectSettings_.g.verbose > 5)
		debug_file << std::endl << "Inv V.dxyz " << std::fixed << std::setprecision(16) << std::setw(26) << *var_cart;
	
}
	

void dna_adjust::LoadVarianceScaling(it_vmsr_t _it_msr, 
	double& vScale, double& pScale, double& lScale, double& hScale,
		bool& scaleMatrix, bool& scalePartial)
{
	vScale = _it_msr->scale4;
	
	// trap vscale equal to very small numbers or zero
	if (vScale < minVal(PRECISION_1E5, projectSettings_.a.fixed_std_dev))		
		vScale = 1.0;

	// non zero matrix scalar?
	scaleMatrix = (fabs(vScale - 1.0) > PRECISION_1E5);	
	
	// partial matrix scalars
	pScale = _it_msr->scale1;
	lScale = _it_msr->scale2;
	hScale = _it_msr->scale3;
	
	// trap vscale equal to very small numbers or zero
	if (pScale < minVal(PRECISION_1E5, projectSettings_.a.fixed_std_dev))
		pScale = 1.0;
	if (lScale < minVal(PRECISION_1E5, projectSettings_.a.fixed_std_dev))
		lScale = 1.0;	
	if (hScale < minVal(PRECISION_1E5, projectSettings_.a.fixed_std_dev))
		hScale = 1.0;

	if (fabs(pScale - 1.0) > PRECISION_1E5 || fabs(lScale - 1.0) > PRECISION_1E5 || fabs(hScale - 1.0) > PRECISION_1E5)
		scalePartial = true;
	else
		scalePartial = false;

	if (scalePartial && scaleMatrix)
	{
		// partial matrix scalars
		pScale *= vScale;
		lScale *= vScale;
		hScale *= vScale;
	}
}
	
// load the GPS variance matrix from the binary file
void dna_adjust::LoadVarianceMatrix_Y(it_vmsr_t _it_msr, matrix_2d* var_cart, const _COORD_TYPE_ coordType)
{
	if (FormInverseVarianceMatrixReduced(_it_msr, var_cart, "LoadVarianceMatrix_Y"))
		return;

	// OK, load the GPS variance matrix from the binary 
	// file, scaling where necessary
	
	it_vmsr_t _it_msr_first(_it_msr);

	//**********************************************************************
	// NOTE!!!
	//
	// Scaling needs careful consideration.  Three things are relevant
	// 1. Cluster may be in X,Y,Z or Lat,Long,Height
	// 2. Single v-scale may be supplied
	// 3. Phi, lambda and height scalars may be supplied
	//
	// if XYZ, scaling is as follows:
	//   if both p,l,h scalars and v-scale are supplied, then
	//		* rotate VCV to geographic,
	//      * apply scalars in single step (i.e. multiply v-scale and p,l,h scalars)
	//		* rotate VCV to cartesian
	//   
	//   if only v-scale supplied, then
	//		* scale VCV
	//
	// else if PLH, scaling is as follows:
	//   if p,l,h scalars are supplied then
	//		* apply scalars in single step (i.e. multiply v-scale and p,l,h scalars)
	//		* rotate VCV to cartesian
	//   if v-scale supplied then
	//		* scale VCV
	//
	// Hence, it makes no difference whether v-scale is performed
	// on a cartesian or geographic VCV.
	//
	//**********************************************************************	
	
	// Units for GPS point clusters are in:
	// - <Coords>LLH</Coords> -> radians^^2 and metres^2
	// - <Coords>XYZ</Coords> -> metres^2
	//
	// An example of LLH is as follows (from "uni_sqrmsr.xml")
	//
	// Worked example:  <SigmaXX>9.402E-09</SigmaXX>
	// 	= 0.000000009402 (radians^2)
	// 	= 0.000096963911 (radians)
	// 	= 0.005555622855 (deg.ddddddd)
	// 	= 0.002000024228 (deg.mmsssss)
	//	= 20.0           (seconds)


	UINT32 covr(0), covc(0);
	UINT32 cluster_pnt, point_count(_it_msr->vectorCount1);
	UINT32 cluster_cov, covariance_count;
	
	// Determine scaling
	double vScale, pScale, lScale, hScale;
	bool scaleMatrix, scalePartial;
	LoadVarianceScaling(_it_msr, vScale, pScale, lScale, hScale,
		scaleMatrix, scalePartial);

	var_cart->redim(point_count * 3, point_count * 3);		// performs zero on creation of new elements
	
	matrix_2d mpositions(point_count * 3, 1);
	it_vstn_t stn1_it;

	for (cluster_pnt=0; cluster_pnt<point_count; ++cluster_pnt)
	{
		covr = cluster_pnt * 3;
		covariance_count = _it_msr->vectorCount2;

		stn1_it = bstBinaryRecords_.begin() + _it_msr->station1;

		// store latitude, longitude and ellipsoid height
        if (scalePartial || coordType == LLH_type_i || 
			coordType == LLh_type_i)
		{
			mpositions.put(covr, 0, stn1_it->currentLatitude);
			mpositions.put(covr+1, 0, stn1_it->currentLongitude);
			mpositions.put(covr+2, 0, stn1_it->currentHeight);
		}

		var_cart->put(covr, covr, _it_msr->term2);		// XX
		_it_msr++;

		var_cart->put(covr, covr+1, _it_msr->term2);	// YX
		var_cart->put(covr+1, covr+1, _it_msr->term3);	// YY
		_it_msr++;	
		
		var_cart->put(covr,   covr+2, _it_msr->term2);	// ZX
		var_cart->put(covr+1, covr+2, _it_msr->term3);	// ZY
		var_cart->put(covr+2, covr+2, _it_msr->term4);	// ZZ

		if (covariance_count == 0)
			break;

		covc = cluster_pnt * 3 + 3;

		for (cluster_cov=0; cluster_cov<covariance_count; ++cluster_cov)	// number of baseline/point covariances
		{
			covr = cluster_pnt * 3;
			_it_msr++;		// X, lat
			var_cart->put(covr, covc, _it_msr->term1);		// m11
			var_cart->put(covr, covc+1, _it_msr->term2);	// m12
			var_cart->put(covr, covc+2, _it_msr->term3);	// m13

			covr++;
			_it_msr++;		// Y, long
			var_cart->put(covr, covc, _it_msr->term1);		// m11
			var_cart->put(covr, covc+1, _it_msr->term2);	// m12
			var_cart->put(covr, covc+2, _it_msr->term3);	// m13

			covr++;
			_it_msr++;		// Z, height
			var_cart->put(covr, covc, _it_msr->term1);		// m11
			var_cart->put(covr, covc+1, _it_msr->term2);	// m12
			var_cart->put(covr, covc+2, _it_msr->term3);	// m13			

			covc+=3;
		}

		_it_msr++;
	}

	bool lowerisClear(true);
	if (scaleMatrix || scalePartial || coordType == LLH_type_i || coordType == LLh_type_i)
	{
		var_cart->filllower();
		lowerisClear = false;
	}

	try {		

		// Scale variance matrix using phi, lambda and height scalars?
		if (scalePartial)
			// If v-scale has been supplied, partial scalars will have been 
			// multiplied by v-scale in LoadVarianceScaling.
			// ScaleGPSVCV_Cluster also converts var_cart to cartesian system
			ScaleGPSVCV_Cluster<double>(*var_cart, var_cart, 
				mpositions, 
				datum_.GetEllipsoidRef(),
				pScale, lScale, hScale, coordType);

		// Propagate variance matrix into cartesian reference frame?
        else if (coordType == LLH_type_i || coordType == LLh_type_i)
			PropagateVariances_GeoCart_Cluster<double>(*var_cart, var_cart,
				mpositions, 
				datum_.GetEllipsoidRef(), 
				true); 		// Geographic -> Cartesian

		if (scaleMatrix && !scalePartial)
			var_cart->scale(vScale);

		// Copy elements of scaled variance matrix to variances held in internal memory (_it_msr)
		// only if VCV has been scaled (by vScale, or by pScale+lScale+hScale) or if VCV was
		// in PLH reference frame.
		// This alleviates having to repeat any scaling for subsequent computations
		// involving variances (i.e. sigma zero, Pelzer's reliability, etc)
		if (scaleMatrix || scalePartial || coordType == LLH_type_i || coordType == LLh_type_i)
			SetGPSVarianceMatrix<it_vmsr_t>(_it_msr_first, *var_cart);
		
		// Perform inverse
		FormInverseVarianceMatrix(var_cart, lowerisClear);

		if (boost::math::isnan(var_cart->get(0, 0)) || boost::math::isinf(var_cart->get(0, 0)))
		{
			std::stringstream ss;
			ss << "Invalid variance matrix:" << std::endl;
			ss << std::setprecision(6) << std::fixed << *var_cart;
			SignalExceptionAdjustment(ss.str(), 0);
		}
	}
	catch (const std::runtime_error& e) {

		// Print error message to adj file and throw exception
		std::stringstream ss;
		printer_->PrintMsrVarianceMatrixException(_it_msr_first, e, ss, "LoadVarianceMatrix_Y", point_count);
		SignalExceptionAdjustment(ss.str(), 0);
	}

	if (projectSettings_.g.verbose > 5)
		debug_file << std::endl << "Inv V.dxyz " << std::fixed << std::setprecision(16) << std::setw(26) << *var_cart;
	
}
	

inline void dna_adjust::AddElementtoDesign(const UINT32& row, const UINT32& col, const double value, matrix_2d* design)
{
	design->put(row, col, value);
}
	

inline void dna_adjust::AddMsrtoDesign(const UINT32& design_row, const UINT32& stn,
				const double& dmdx, const double& dmdy, const double& dmdz, matrix_2d* design)
{
		// design matrix elements for dA/dX1, dA/dY1, dA/dZ1 for 1 station
	AddElementtoDesign(design_row, stn, dmdx, design);
	AddElementtoDesign(design_row, stn+1, dmdy, design);
	AddElementtoDesign(design_row, stn+2, dmdz, design);
}
	

inline void dna_adjust::AddMsrtoDesign_L(const UINT32& design_row, const UINT32& stn1, const UINT32& stn2, 
				const double dmdx1, const double dmdy1, const double dmdz1, 
				const double dmdx2, const double dmdy2, const double dmdz2, matrix_2d* design)
{
	// design matrix dA/dX1, dA/dY1, dA/dZ1
	AddMsrtoDesign(design_row, stn1, dmdx1, dmdy1, dmdz1, design);	
	// design matrix dA/dX2, dA/dY2, dA/dZ2
	AddMsrtoDesign(design_row, stn2, dmdx2, dmdy2, dmdz2, design);
}
	

inline void dna_adjust::AddMsrtoDesign_BCEKMSVZ(const UINT32& design_row, const UINT32& stn1, const UINT32& stn2, 
				const double dmdx, const double dmdy, const double dmdz, matrix_2d* design)
{
	// design matrix dA/dX1, dA/dY1, dA/dZ1
	AddMsrtoDesign(design_row, stn1, dmdx, dmdy, dmdz, design);	
	// design matrix dA/dX2, dA/dY2, dA/dZ2
	AddMsrtoDesign(design_row, stn2, -dmdx, -dmdy, -dmdz, design);
}
	

void dna_adjust::AddMsrtoMeasMinusComp(pit_vmsr_t _it_msr, const UINT32& design_row, const double comp_msr, 
	matrix_2d* measMinusComp, bool printBlock)
{
	double mmc((*_it_msr)->term1 - comp_msr);
	
	switch ((*_it_msr)->measType)
	{
	case 'A':
	case 'B':
	case 'D':
	case 'K':
		// check if residual is close to +/- 360 (make close to 0)
		if (mmc < -5.5)
			mmc += TWO_PI;
		else if (mmc > 5.5)
			mmc -= TWO_PI;
	}

	// measured - computed
	measMinusComp->put(design_row, 0, mmc);		

	if (projectSettings_.g.verbose > 5)
	{
		if (printBlock)
			debug_file << bstBinaryRecords_.at((*_it_msr)->station1).stationName <<
				" - " << bstBinaryRecords_.at((*_it_msr)->station2).stationName;
		debug_file << std::endl;
		debug_file << "d" << (*_it_msr)->measType << " " << std::fixed << std::setprecision(16) << std::setw(26) << 
			mmc <<
			" = " << std::setw(26) << (*_it_msr)->term1 << 
			" - " << std::setw(26) << comp_msr << "  ";
	}
}
	

void dna_adjust::UpdateDesignNormalMeasMatrices_A(pit_vmsr_t _it_msr, UINT32& design_row, const UINT32& block,
											  matrix_2d* measMinusComp, matrix_2d* estimatedStations, 
											  matrix_2d* normals, matrix_2d* design, matrix_2d* AtVinv, bool buildnewMatrices)
{
	UINT32 stn1(GetBlkMatrixElemStn1(block, _it_msr)); 
	UINT32 stn2(GetBlkMatrixElemStn2(block, _it_msr));
	UINT32 stn3(GetBlkMatrixElemStn3(block, _it_msr));

	it_vstn_t_const stn1_it(bstBinaryRecords_.begin() + (*_it_msr)->station1);

	double direction12, direction13, local_12e, local_12n, local_13e, local_13n;

	// compute angle 1 -> 2 -> 3 from estimated coordinates
	double comp_msr(HorizontalAngle(
		estimatedStations->get(stn1, 0),		// X1
		estimatedStations->get(stn1+1, 0),		// Y1
		estimatedStations->get(stn1+2, 0),		// Z1
		estimatedStations->get(stn2, 0), 		// X2
		estimatedStations->get(stn2+1, 0),		// Y2
		estimatedStations->get(stn2+2, 0),		// Z2
		estimatedStations->get(stn3, 0), 		// X3
		estimatedStations->get(stn3+1, 0),		// Y3
		estimatedStations->get(stn3+2, 0),		// Z3
		stn1_it->currentLatitude,
		stn1_it->currentLongitude,
		&direction12, &direction13, 
		&local_12e, &local_12n, &local_13e, &local_13n));
	
	// Initialise measurement.  No need to test if no further calculations are 
	// required (as in stage mode), as this is done later (below)
	InitialiseMeasurement(_it_msr, buildnewMatrices);

	if (buildnewMatrices)
	{
		// deflections available?
		if (fabs(stn1_it->verticalDef) > E4_SEC_DEFLECTION || fabs(stn1_it->meridianDef) > E4_SEC_DEFLECTION)
		{
			it_vstn_t stn2_it(bstBinaryRecords_.begin() + (*_it_msr)->station2);
			it_vstn_t stn3_it(bstBinaryRecords_.begin() + (*_it_msr)->station3);

			/////////////////////////////////////////////////////////////////////////////////
			// Angles (observed or derived from directions) must be corrected for deflection 
			// of the vertical via "Laplace correction".  This correction requires zenith 
			// distance (zenith12, zenith13) and geodetic azimuth (direction12, direction13), 
			// both of which must be computed from coordinates.

			////////////////////////////////////////////////////////////////////////////
			// Compute zenith distance 1 -> 2
			double zenith12(ZenithDistance<double>(
				estimatedStations->get(stn1, 0),			// X1
				estimatedStations->get(stn1+1, 0),			// Y1
				estimatedStations->get(stn1+2, 0),			// Z1
				estimatedStations->get(stn2, 0), 			// X2
				estimatedStations->get(stn2+1, 0),			// Y2
				estimatedStations->get(stn2+2, 0),			// Z2
				stn1_it->currentLatitude,
				stn1_it->currentLongitude,
				stn2_it->currentLatitude,
				stn2_it->currentLongitude,
				(*_it_msr)->term3,							// instrument height
				(*_it_msr)->term4));						// target height

			////////////////////////////////////////////////////////////////////////////
			// Compute zenith distance 1 -> 3
			double zenith13(ZenithDistance<double>(
				estimatedStations->get(stn1, 0),			// X1
				estimatedStations->get(stn1+1, 0),			// Y1
				estimatedStations->get(stn1+2, 0),			// Z1
				estimatedStations->get(stn3, 0), 			// X2
				estimatedStations->get(stn3+1, 0),			// Y2
				estimatedStations->get(stn3+2, 0),			// Z2
				stn1_it->currentLatitude,
				stn1_it->currentLongitude,
				stn3_it->currentLatitude,
				stn3_it->currentLongitude,
				(*_it_msr)->term3,							// instrument height
				(*_it_msr)->term4));						// target height

			// Laplace correction 1 -> 2 -> 3
			(*_it_msr)->preAdjCorr = HzAngleDeflectionCorrection<double>(		
				direction12,								// geodetic azimuth 1 -> 2
				zenith12,									// zenith distance 1 -> 2
				direction13,								// geodetic azimuth 1 -> 3
				zenith13,									// zenith distance 1 -> 3
				stn1_it->verticalDef,						// deflection in prime vertical
				stn1_it->meridianDef);						// deflection in prime meridian

			(*_it_msr)->term1 -= (*_it_msr)->preAdjCorr;	// apply deflection correction
		}
		else
			(*_it_msr)->preAdjCorr = 0.0;

		if (projectSettings_.a.stage)
			return;
	}

	// Update measured minus computed value
	AddMsrtoMeasMinusComp(_it_msr, design_row, comp_msr, measMinusComp);

	double cos_lat(cos(stn1_it->currentLatitude));
	double sin_lat(sin(stn1_it->currentLatitude));
	double cos_long(cos(stn1_it->currentLongitude));
	double sin_long(sin(stn1_it->currentLongitude));
	double sinlat_coslong(sin_lat*cos_long);
	double sinlat_sinlong(sin_lat*sin_long);
	double cos2_dir12_div_n122(cos(direction12) * cos(direction12) / (local_12n * local_12n));
	double cos2_dir13_div_n132(cos(direction13) * cos(direction13) / (local_13n * local_13n));

	// compute partial derivatives for normals
	// design matrix dA/dX1, dA/dY1, dA/dZ1
	AddElementtoDesign(design_row, stn1, 
		cos2_dir13_div_n132 * (local_13n * sin_long - local_13e * sinlat_coslong) - 
		cos2_dir12_div_n122 * (local_12n * sin_long - local_12e * sinlat_coslong),
		design);																		// X1
	AddElementtoDesign(design_row, stn1+1,  
		cos2_dir13_div_n132 * (-local_13n * cos_long - local_13e * sinlat_sinlong) - 
		cos2_dir12_div_n122 * (-local_12n * cos_long - local_12e * sinlat_sinlong),
		design);																		// Y1
	AddElementtoDesign(design_row, stn1+2, 
		cos2_dir13_div_n132 * local_13e * cos_lat - 
		cos2_dir12_div_n122 * local_12e * cos_lat, design);								// Z1

	// design matrix dA/dX2, dA/dY2, dA/dZ2
	AddElementtoDesign(design_row, stn2, 
		cos2_dir12_div_n122 * (local_12n * sin_long - local_12e * sinlat_coslong),
		design);																		// X2
	AddElementtoDesign(design_row, stn2+1, 
		cos2_dir12_div_n122 * (-local_12n * cos_long - local_12e * sinlat_sinlong),
		design);																		// Y2
	AddElementtoDesign(design_row, stn2+2, 
		cos2_dir12_div_n122 * local_12e * cos_lat, design);								// Z2

	// design matrix dA/dX3, dA/dY3, dA/dZ3
	AddElementtoDesign(design_row, stn3, 
		-cos2_dir13_div_n132 * (local_13n * sin_long - local_13e * sinlat_coslong),
		design);																		// X3
	AddElementtoDesign(design_row, stn3+1, 
		-cos2_dir13_div_n132 * (-local_13n * cos_long - local_13e * sinlat_sinlong),
		design);																		// Y3
	AddElementtoDesign(design_row, stn3+2, 
		-cos2_dir13_div_n132 * local_13e * cos_lat, design);							// Z3
	
	switch ((*_it_msr)->measType)
	{
	case 'D':
		// return to UpdateDesignNormalMeasMatrices_D(...)
		return;
	}

	UpdateAtVinv(_it_msr, stn1, stn2, stn3, design_row, design, AtVinv, buildnewMatrices);

	if (buildnewMatrices)
		// Add weighted measurement contributions to normal matrix
		UpdateNormals_A(stn1, stn2, stn3, design_row, normals, design, AtVinv);
	else
		design_row++;
}


void dna_adjust::UpdateDesignNormalMeasMatrices_BK(pit_vmsr_t _it_msr, UINT32& design_row, const UINT32& block,
											  matrix_2d* measMinusComp, matrix_2d* estimatedStations, 
											  matrix_2d* normals, matrix_2d* design, matrix_2d* AtVinv, bool buildnewMatrices)
{
	UINT32 stn1(GetBlkMatrixElemStn1(block, _it_msr)); 
	UINT32 stn2(GetBlkMatrixElemStn2(block, _it_msr));

	it_vstn_t_const stn1_it(bstBinaryRecords_.begin() + (*_it_msr)->station1);
	it_vstn_t_const stn2_it(bstBinaryRecords_.begin() + (*_it_msr)->station2);

	double local_12e, local_12n;

	// compute bearing from estimated coordinates
	double comp_msr(Direction(
		estimatedStations->get(stn1, 0),		// X1
		estimatedStations->get(stn1+1, 0),		// Y1
		estimatedStations->get(stn1+2, 0),		// Z1
		estimatedStations->get(stn2, 0), 		// X2
		estimatedStations->get(stn2+1, 0),		// Y2
		estimatedStations->get(stn2+2, 0),		// Z2
		stn1_it->currentLatitude,
		stn1_it->currentLongitude,
		&local_12e, &local_12n));
	
	// Initialise measurement.  No need to test if no further calculations are 
	// required (as in stage mode), as this is done later (below)
	InitialiseMeasurement(_it_msr, buildnewMatrices);
	
	if (buildnewMatrices)
	{
		// deflections available?
		if ((*_it_msr)->measType == 'K' &&						// Astro
			(fabs(stn1_it->verticalDef) > E4_SEC_DEFLECTION ||	// deflections available?
			 fabs(stn1_it->meridianDef) > E4_SEC_DEFLECTION))	
		{
			////////////////////////////////////////////////////////////////////////////
			// Astronomic azimuths must be corrected for deflection of the vertical via
			// "Laplace correction".  This correction requires zenith distance and geodetic 
			// azimuth (comp_msr), both of which must be computed from coordinates.
		
			////////////////////////////////////////////////////////////////////////////
			// Compute zenith distance
			double zenith(ZenithDistance<double>(
				estimatedStations->get(stn1, 0),					// X1
				estimatedStations->get(stn1+1, 0),					// Y1
				estimatedStations->get(stn1+2, 0),					// Z1
				estimatedStations->get(stn2, 0), 					// X2
				estimatedStations->get(stn2+1, 0),					// Y2
				estimatedStations->get(stn2+2, 0),					// Z2
				stn1_it->currentLatitude,
				stn1_it->currentLongitude,
				stn2_it->currentLatitude,
				stn2_it->currentLongitude,
				(*_it_msr)->term3,									// instrument height
				(*_it_msr)->term4));								// target height

			(*_it_msr)->preAdjCorr = LaplaceCorrection<double>(		// Laplace correction
				comp_msr,											// geodetic azimuth
				zenith,												// zenith distance
				stn1_it->verticalDef,								// deflection in prime vertical
				stn1_it->meridianDef,								// deflection in prime meridian
				stn1_it->currentLatitude);

			// scale4 must be used since term3 and term4 are station and target heights respectively
			(*_it_msr)->term1 -= (*_it_msr)->preAdjCorr;			// apply deflection correction
		}
		else
			(*_it_msr)->preAdjCorr = 0.0;

		if (projectSettings_.a.stage)
			return;
	}

	// Update measured minus computed value
	AddMsrtoMeasMinusComp(_it_msr, design_row, comp_msr, measMinusComp);

	////////////////////////////////////////////////////////////////////////////
	// compute partial derivatives for normals
	double cos_lat(cos(stn1_it->currentLatitude));
	double sin_lat(sin(stn1_it->currentLatitude));
	double cos_long(cos(stn1_it->currentLongitude));
	double sin_long(sin(stn1_it->currentLongitude));
	double sinlat_coslong(sin_lat*cos_long);
	double sinlat_sinlong(sin_lat*sin_long);
	double cos2_dir12_div_n122(cos(comp_msr) * cos(comp_msr) / (local_12n * local_12n));
	
	// Add partial derivatives dA/dX1, dA/dY1, dA/dZ1 to design matrix
	AddMsrtoDesign_BCEKMSVZ(design_row, stn1, stn2,
		cos2_dir12_div_n122 * (local_12n * sin_long - local_12e * sinlat_coslong),
		cos2_dir12_div_n122 * (-local_12n * cos_long - local_12e * sinlat_sinlong),
		cos2_dir12_div_n122 * local_12e * cos_lat,
		design);
	
	// Update AtVinv based on new design matrix elements
	UpdateAtVinv(_it_msr, stn1, stn2, 0, design_row, design, AtVinv, buildnewMatrices);

	if (buildnewMatrices)
		// Add weighted measurement contributions to normal matrix
		UpdateNormals_BCEKLMSVZ(stn1, stn2, design_row, normals, design, AtVinv);
	else
		design_row++;
}
	

void dna_adjust::UpdateDesignNormalMeasMatrices_C(pit_vmsr_t _it_msr, UINT32& design_row, const UINT32& block,
											  matrix_2d* measMinusComp, matrix_2d* estimatedStations, 
											  matrix_2d* normals, matrix_2d* design, matrix_2d* AtVinv, bool buildnewMatrices)
{
	// Initialise measurement and test if no further calculations are 
	// required (as in stage mode)
	if (InitialiseMeasurement(_it_msr, buildnewMatrices))
		return;

	// As a C measurement is a direct vector between two points in the cartesian frame,
	// there will be no pre adjustment correction
	(*_it_msr)->preAdjCorr = 0.;

	// Now call UpdateDesignNormalMeasMatrices_CEM
	UpdateDesignNormalMeasMatrices_CEM(_it_msr, design_row, block,
		measMinusComp, estimatedStations, normals, design, AtVinv, buildnewMatrices);
}
	

void dna_adjust::UpdateDesignNormalMeasMatrices_CEM(pit_vmsr_t _it_msr, UINT32& design_row, const UINT32& block,
											  matrix_2d* measMinusComp, matrix_2d* estimatedStations, 
											  matrix_2d* normals, matrix_2d* design, matrix_2d* AtVinv, bool buildnewMatrices)
{
	UINT32 stn1(GetBlkMatrixElemStn1(block, _it_msr)); 
	UINT32 stn2(GetBlkMatrixElemStn2(block, _it_msr));

	it_vstn_t_const stn1_it(bstBinaryRecords_.begin() + (*_it_msr)->station1);
	it_vstn_t_const stn2_it(bstBinaryRecords_.begin() + (*_it_msr)->station2);

	double dX, dY, dZ;

	// calculated distance
	double comp_msr(EllipsoidChordDistance<double>(
		estimatedStations->get(stn1, 0),
		estimatedStations->get(stn1+1, 0),
		estimatedStations->get(stn1+2, 0),
		estimatedStations->get(stn2, 0),
		estimatedStations->get(stn2+1, 0),
		estimatedStations->get(stn2+2, 0),
		stn1_it->currentLatitude,
		stn2_it->currentLatitude,
		stn1_it->currentHeight,
		stn2_it->currentHeight,
		&dX, &dY, &dZ, 
		datum_.GetEllipsoidRef()));
	
	// Update measured minus computed value
	AddMsrtoMeasMinusComp(_it_msr, design_row, comp_msr, measMinusComp);

	// Add partial derivatives dA/dX1, dA/dY1, dA/dZ1 to design matrix
	AddMsrtoDesign_BCEKMSVZ(design_row, stn1, stn2,
		-dX/comp_msr, -dY/comp_msr,	-dZ/comp_msr,
		design);

	// Update AtVinv based on new design matrix elements
	UpdateAtVinv(_it_msr, stn1, stn2, 0, design_row, design, AtVinv, buildnewMatrices);

	if (buildnewMatrices)
		// Add weighted measurement contributions to normal matrix
		UpdateNormals_BCEKLMSVZ(stn1, stn2, design_row, normals, design, AtVinv);
	else
		design_row++;
}
	

void dna_adjust::UpdateDesignNormalMeasMatrices_D(pit_vmsr_t _it_msr, UINT32& design_row, const UINT32& block,
											  matrix_2d* measMinusComp, matrix_2d* estimatedStations, 
											  matrix_2d* normals, matrix_2d* design, matrix_2d* AtVinv, bool buildnewMatrices)
{
	it_vmsr_t _it_msr_first(*_it_msr);
	UINT32 design_row_begin(design_row);

	UINT32 a, angle_count((*_it_msr)->vectorCount2 - 1);		// number of directions excluding the RO
	UINT32 skip(0), ignored((*_it_msr)->vectorCount1 - (*_it_msr)->vectorCount2);
	
	vmsr_t angleRec;
	angleRec.push_back(*(*_it_msr));
	it_vmsr_t it_angle(angleRec.begin());

	double previousDirection((*_it_msr)->term1);

	(*_it_msr)++;

	// set derived angle, variance and covariance to the binary records
	// term1 = measured direction
	// term2 = variance (direction)
	// term3 = instrument height (not used)
	// term4 = target height (not used)
	// scale1 = derived angle corrected for deflection of the vertical
	// scale2 = variance (angle)
	// scale3 = covariance (angle) - for the context of vmsr_t angleRec only, so as to 
	//          properly form the normals from covariances formed from directions SDs
	// vectorCount2 = number of non-ignored directions
	// preAdjMeas = original derived angle
	
	if (projectSettings_.g.verbose > 6)
		debug_file << "Reduced angles from raw directions: ";

	try
	{
		for (a=0; a<angle_count; ++a)
		{
			// cater for ignored directions
			if ((*_it_msr)->ignore)
			{
				while (skip < ignored)
				{					
					skip++;
					(*_it_msr)++;
					if (!(*_it_msr)->ignore)
						break;						
				}
			}

			it_angle->station3 = (*_it_msr)->station2;
			
			if (buildnewMatrices)
			{
				// Was this measurement reduced from previous
				// adjustment, which was serialised to disk?
				if (bms_meta_.reduced)
					// derived angle uncorrected for deflection of vertical
					it_angle->preAdjMeas = (*_it_msr)->preAdjMeas;
				else
				{
					// Derive the angle from two directions
					it_angle->term1 = (*_it_msr)->term1 - previousDirection;
					if (it_angle->term1 < 0)
						it_angle->term1 += TWO_PI;
					if (it_angle->term1 > TWO_PI)
						it_angle->term1 -= TWO_PI;

					if (projectSettings_.g.verbose > 6)
					{
						debug_file << std::endl << "angle " <<
							StringFromT(RadtoDms(it_angle->term1), 6) << " = " <<
							StringFromT(RadtoDms((*_it_msr)->term1), 6) << " - " <<
							StringFromT(RadtoDms(previousDirection), 6) << std::endl;
					}
				}
			}
			else
				// derived angle corrected for deflection of vertical
				it_angle->term1 = (*_it_msr)->scale1;
			
			// normals not needed
			UpdateDesignNormalMeasMatrices_A(&it_angle, design_row, block,
				measMinusComp, estimatedStations, 0, design, AtVinv, buildnewMatrices);

			if (buildnewMatrices)
			{
				// Update derived angle, corrected for deflection of vertical
				(*_it_msr)->scale1 = it_angle->term1;
				// Update derived angle, uncorrected for deflection of vertical
				(*_it_msr)->preAdjMeas = it_angle->preAdjMeas;
				// Update correction for deflection of the vertical
				(*_it_msr)->preAdjCorr = it_angle->preAdjCorr;
			}

			if (a + 1 == angle_count)
			{
				// check for any ignored directions
				if (_it_msr_first->vectorCount1 != _it_msr_first->vectorCount2)
					(*_it_msr) = _it_msr_first + _it_msr_first->vectorCount1;
				break;
			}

			if (buildnewMatrices)
				previousDirection = (*_it_msr)->term1;

			// prepare for next angle
			angleRec.push_back(*(*_it_msr));
			it_angle = angleRec.end() - 1;

			(*_it_msr)++;	
			design_row++;
		}
	}
	catch (...) {

		// Print error message to adj file and throw exception
		std::stringstream ss;
		ss << "UpdateDesignNormalMeasMatrices_D(): An error was encountered whilst" << std::endl <<
			"  updating the normal matrices" << std::endl;
		SignalExceptionAdjustment(ss.str(), block);
	}

	// Load apriori variance matrix, and assign to binary measurement
	matrix_2d var_dirn(angle_count, angle_count);
	LoadVarianceMatrix_D(_it_msr_first, &var_dirn, buildnewMatrices);

	if (buildnewMatrices && projectSettings_.a.stage)
		return;

	UINT32 stn1, stn2, stn3;
	it_angle = angleRec.begin();
	design_row = design_row_begin;

	// clear rows in AtVinv relating to this measurement
	AtVinv->zero(0, design_row, AtVinv->rows(), angle_count);

	if (!buildnewMatrices)
		debug_file << "block " << block + 1 << ", design_row " << design_row + 1 << ", directions var" << var_dirn << std::endl;

	// Update AtVinv based on new design matrix elements
	for (a=0; a<angle_count; ++a)																  // for each angle
	{
		stn1 = (v_blockStationsMap_.at(block)[it_angle->station1] * 3); 
		stn2 = (v_blockStationsMap_.at(block)[it_angle->station2] * 3);
		stn3 = (v_blockStationsMap_.at(block)[it_angle->station3] * 3);
		
		// Update AtVinv
		UpdateAtVinv_D(stn1, stn2, stn3, a, angle_count, 
			design_row, design_row_begin, 
			&var_dirn, design, AtVinv);

		design_row++;
		it_angle++;
	}

	if (buildnewMatrices)
		// Add weighted measurement contributions to normal matrix
		UpdateNormals_D(block, _it_msr_first, design_row_begin, normals, design, AtVinv);
}

void dna_adjust::UpdateDesignNormalMeasMatrices_E(pit_vmsr_t _it_msr, UINT32& design_row, const UINT32& block,
											  matrix_2d* measMinusComp, matrix_2d* estimatedStations, 
											  matrix_2d* normals, matrix_2d* design, matrix_2d* AtVinv, bool buildnewMatrices)
{
	// Initialise measurement and test if no further calculations are 
	// required (as in stage mode)
	if (InitialiseMeasurement(_it_msr, buildnewMatrices))
		return;

	UINT32 stn1(GetBlkMatrixElemStn1(block, _it_msr)); 
	UINT32 stn2(GetBlkMatrixElemStn2(block, _it_msr));

	it_vstn_t_const stn1_it(bstBinaryRecords_.begin() + (*_it_msr)->station1);
	it_vstn_t_const stn2_it(bstBinaryRecords_.begin() + (*_it_msr)->station2);

	// Note: Each measured Ellipsoid arc is reduced to an ellipsoid chord on first run.
	// Whist it would be more efficient to not reduce this distance on each iteration,
	// the ellipsoid chord is recomputed every time using the original measurement
	// to cater for (potentially) large changes in latitude and longitude, 
	
	// Reduce Ellipsoid arc to Ellipsoid chord
	(*_it_msr)->term1 = EllipsoidArctoEllipsoidChord<double>(
		(*_it_msr)->preAdjMeas,
		estimatedStations->get(stn1, 0), 
		estimatedStations->get(stn1+1, 0),
		estimatedStations->get(stn1+2, 0),
		estimatedStations->get(stn2, 0),
		estimatedStations->get(stn2+1, 0),
		estimatedStations->get(stn2+2, 0),
		stn1_it->currentLatitude,
		stn1_it->currentLongitude,
		stn2_it->currentLatitude,
		datum_.GetEllipsoidRef());

	(*_it_msr)->preAdjCorr = (*_it_msr)->term1 - (*_it_msr)->preAdjMeas;

	// Now that the ellipsoid arc has been reduced to a chord, call UpdateDesignNormalMeasMatrices_CEM
	UpdateDesignNormalMeasMatrices_CEM(_it_msr, design_row, block,
		measMinusComp, estimatedStations, normals, design, AtVinv, buildnewMatrices);
}

void dna_adjust::UpdateDesignMeasMatrices_GX(pit_vmsr_t _it_msr, UINT32& design_row,
											matrix_2d* measMinusComp, matrix_2d* estimatedStations, matrix_2d* design,
											const UINT32& stn1, const UINT32& stn2, bool buildnewMatrices)
{
	// If this method is called via PrepareAdjustment() and the adjustment 
	// mode is staged, then don't calculate measured-computed values or
	// fill design matrix.  This will be done during an adjustment
	// via AdjustPhasedForward().
	if (buildnewMatrices && projectSettings_.a.stage)
		return;

	// For all adjustment modes, when this method is called during an adjustment
	// to update the normals, only the the measured-computed values need to be 
	// updated.  This is because the (Jacobian) design matrix elements for GPS
	// are unity (1 or -1) and do not change as coordinates are updated (unlike 
	// the Jacobian elements formed for other measurements).  Hence, for all 
	// adjustment modes except staged, the design matrix is updated once via
	// PrepareAdjustment.  For staged adjustments, the design matrix is updated
	// on each iteration.

	// Add X elements to measured minus computed
	AddMsrtoMeasMinusComp(_it_msr, design_row, 
		(estimatedStations->get(stn2, 0) - estimatedStations->get(stn1, 0)), 
		measMinusComp);

	if (buildnewMatrices || projectSettings_.a.stage)
	{
		// Add X elements to design matrix
		AddElementtoDesign(design_row, stn1, -1., design);		// X1
		AddElementtoDesign(design_row, stn2, 1., design);		// X2
		(*_it_msr)->preAdjMeas = (*_it_msr)->term1;
	}

	// move to Y element
	design_row++;
	(*_it_msr)++;

	// Add Y elements to measured minus computed
	AddMsrtoMeasMinusComp(_it_msr, design_row, 
		(estimatedStations->get(stn2+1, 0) - estimatedStations->get(stn1+1, 0)), 
		measMinusComp, false);
			
	if (buildnewMatrices || projectSettings_.a.stage)
	{
		// Add Y elements to design matrix
		AddElementtoDesign(design_row, stn1+1, -1., design);	// Y1
		AddElementtoDesign(design_row, stn2+1, 1., design);		// Y2
		(*_it_msr)->preAdjMeas = (*_it_msr)->term1;
	}
		
	// move to Z element
	design_row++;
	(*_it_msr)++;
	
	// Add Z elements to measured minus computed
	AddMsrtoMeasMinusComp(_it_msr, design_row, 
		(estimatedStations->get(stn2+2, 0) - estimatedStations->get(stn1+2, 0)), 
		measMinusComp, false);
    if (projectSettings_.g.verbose > 5)
        debug_file << std::endl;
		
	if (buildnewMatrices || projectSettings_.a.stage)
	{
		// Add Z elements to design matrix
		AddElementtoDesign(design_row, stn1+2, -1., design);	// Z1
		AddElementtoDesign(design_row, stn2+2, 1., design);		// Z2
		(*_it_msr)->preAdjMeas = (*_it_msr)->term1;
	}
}

void dna_adjust::UpdateDesignNormalMeasMatrices_G(pit_vmsr_t _it_msr, UINT32& design_row, const UINT32& block,
											  matrix_2d* measMinusComp, matrix_2d* estimatedStations, 
											  matrix_2d* normals, matrix_2d* design, matrix_2d* AtVinv, bool buildnewMatrices)
{
	it_vmsr_t _it_msr_first(*_it_msr);
	UINT32 stn1(GetBlkMatrixElemStn1(block, _it_msr)); 
	UINT32 stn2(GetBlkMatrixElemStn2(block, _it_msr));
	UINT32 design_row_begin(design_row);
	
	UpdateDesignMeasMatrices_GX(_it_msr, design_row,
			measMinusComp, estimatedStations, design,
			stn1, stn2, buildnewMatrices);

	// For all adjustment modes except staged, when this method is called during 
	// an adjustment to update the normals, the weighted design matrix (AtVinv)
	// is not updated.  This is because the (Jacobian) design matrix elements for 
	// GPS are unity (1 or -1) and do not change as coordinates are updated (unlike 
	// the Jacobian elements formed for other measurements).  For all adjustment
	// modes except staged, the weighted design matrix and normals are updated once
	// via PrepareAdjustment.  For staged adjustments, the weighted design matrix 
	// and normals are updated on each iteration.


	// For first time run or staged adjustment mode
	if (buildnewMatrices || projectSettings_.a.stage)
	{
		// Load GPS Variance matrix.  For the first run, scaling is applied
		// and the variances are written to the binary measurements list.
		matrix_2d var_cart(3, 3);
		LoadVarianceMatrix_G(_it_msr_first, &var_cart);

		if (buildnewMatrices && projectSettings_.a.stage)
			return;
	
		// Build  At * V-1
		AtVinv->replace(stn1, design_row_begin, var_cart * -1);
		AtVinv->replace(stn2, design_row_begin, var_cart);

		if (!projectSettings_.a.stage)
			// Add weighted measurement contributions to normal matrix
			UpdateNormals_G(stn1, stn2, design_row_begin, normals, AtVinv);
	}
	
	design_row++;
}
	

void dna_adjust::UpdateDesignNormalMeasMatrices_M(pit_vmsr_t _it_msr, UINT32& design_row, const UINT32& block,
											  matrix_2d* measMinusComp, matrix_2d* estimatedStations, 
											  matrix_2d* normals, matrix_2d* design, matrix_2d* AtVinv, bool buildnewMatrices)
{
	// Initialise measurement and test if no further calculations are 
	// required (as in stage mode)
	if (InitialiseMeasurement(_it_msr, buildnewMatrices))
		return;

	it_vstn_t_const stn1_it(bstBinaryRecords_.begin() + (*_it_msr)->station1);
	it_vstn_t_const stn2_it(bstBinaryRecords_.begin() + (*_it_msr)->station2);

	// Note: Each measured MSL arc is reduced to an ellipsoid chord on the first run.
	// However, to cater for (potentially) large changes in latitude and longitude, recompute the
	// chord using the original measurement

	// Reduce MSL arc to ellipsoid chord
	(*_it_msr)->term1 = MSLArctoEllipsoidChord<double>(
		(*_it_msr)->preAdjMeas,									// use measured MSL arc
		stn1_it->currentLatitude, stn2_it->currentLatitude,
		stn1_it->geoidSep, stn2_it->geoidSep,
		datum_.GetEllipsoidRef());

	(*_it_msr)->preAdjCorr = (*_it_msr)->term1 - (*_it_msr)->preAdjMeas;

	// Now that the MSL arc has been reduced to a chord, call UpdateDesignNormalMeasMatrices_CEM
	UpdateDesignNormalMeasMatrices_CEM(_it_msr, design_row, block,
		measMinusComp, estimatedStations, normals, design, AtVinv, buildnewMatrices);
}

// Like zenith distances and vertical angles, the relationship between slope distances and the
// coordinates of p1 and p2 requires a little extra work to take into consideration instrument
// and target height.  For the measurements-minus-computed vector, the "computed" distance
// is the true distance between the instrument and target, and so must take into consideration 
// instrument and target heights.  However, the dX, dY, dZ components for the partial 
// derivatives represent the true geometric difference between the two stations (not 
// instrument and target).
void dna_adjust::UpdateDesignNormalMeasMatrices_S(pit_vmsr_t _it_msr, UINT32& design_row, const UINT32& block,
											  matrix_2d* measMinusComp, matrix_2d* estimatedStations, 
											  matrix_2d* normals, matrix_2d* design, matrix_2d* AtVinv, bool buildnewMatrices)
{
	// preAdjMeas is used to store original measured MSL arc distance
	if (buildnewMatrices)
	{
		// This is the first time adjust has run, so
		// back up the raw measurement
		(*_it_msr)->preAdjMeas = (*_it_msr)->term1;

		if (projectSettings_.a.stage)
			return;
	}

	UINT32 stn1(GetBlkMatrixElemStn1(block, _it_msr)); 
	UINT32 stn2(GetBlkMatrixElemStn2(block, _it_msr));

	it_vstn_t_const stn1_it(bstBinaryRecords_.begin() + (*_it_msr)->station1);

	// compute dX, dY, dZ for instrument height (ih) and target height (th)
	double dXih, dYih, dZih, dXth, dYth, dZth;
	CartesianElementsFromInstrumentHeight((*_it_msr)->term3,				// instrument height
		&dXih, &dYih, &dZih, 
		stn1_it->currentLatitude,
		stn1_it->currentLongitude);
	CartesianElementsFromInstrumentHeight((*_it_msr)->term4,				// target height
		&dXth, &dYth, &dZth,
		stn1_it->currentLatitude,
		stn1_it->currentLongitude);

	// compute distance between instrument and target, taking into consideration
	// instrument height, target height, and vector between stations
	double dX(estimatedStations->get(stn2, 0) - estimatedStations->get(stn1, 0) + dXth - dXih);
	double dY(estimatedStations->get(stn2+1, 0) - estimatedStations->get(stn1+1, 0) + dYth - dYih);
	double dZ(estimatedStations->get(stn2+2, 0) - estimatedStations->get(stn1+2, 0) + dZth - dZih);

	// calculated distance
	double comp_msr(magnitude(dX, dY, dZ));

	// Add measured minus computed value
	AddMsrtoMeasMinusComp(_it_msr, design_row, comp_msr, measMinusComp);

	// Add partial derivatives dA/dX1, dA/dY1, dA/dZ1 to design matrix
	AddMsrtoDesign_BCEKMSVZ(design_row, stn1, stn2,
		-dX/comp_msr, -dY/comp_msr,	-dZ/comp_msr,
		design);

	// Update AtVinv based on new design matrix elements
	UpdateAtVinv(_it_msr, stn1, stn2, 0, design_row, design, AtVinv, buildnewMatrices);

	if (buildnewMatrices)
		// Add weighted measurement contributions to normal matrix
		UpdateNormals_BCEKLMSVZ(stn1, stn2, design_row, normals, design, AtVinv);
	else
		design_row++;
}

// Zenith distance
//
// Like vertical angles, the relationship between zenith distances and the
// coordinates of p1 and p2 requires a little extra work to take into consideration instrument
// and target height.  For the measurements-minus-computed vector, the "computed" zenith distance
// is the true angle between the ellipsoid normal and instrument-target vector, and so must take 
// into consideration instrument and target heights.  However, the dX, dY, dZ components for the
// partial derivatives represent the true geometric difference between the two stations (not 
// instrument and target).
void dna_adjust::UpdateDesignNormalMeasMatrices_V(pit_vmsr_t _it_msr, UINT32& design_row, const UINT32& block,
											  matrix_2d* measMinusComp, matrix_2d* estimatedStations, 
											  matrix_2d* normals, matrix_2d* design, matrix_2d* AtVinv, bool buildnewMatrices)
{
	UINT32 stn1(GetBlkMatrixElemStn1(block, _it_msr)); 
	UINT32 stn2(GetBlkMatrixElemStn2(block, _it_msr));

	it_vstn_t_const stn1_it(bstBinaryRecords_.begin() + (*_it_msr)->station1);
	it_vstn_t_const stn2_it(bstBinaryRecords_.begin() + (*_it_msr)->station2);

	double local_12e, local_12n, local_12up;

	// Initialise measurement.  No need to test if no further calculations are 
	// required (as in stage mode), as this is done later (below)
	InitialiseMeasurement(_it_msr, buildnewMatrices);

	if (buildnewMatrices)
	{
		// deflections available?
		if (fabs(stn1_it->verticalDef) > E4_SEC_DEFLECTION || fabs(stn1_it->meridianDef) > E4_SEC_DEFLECTION)
		{
			////////////////////////////////////////////////////////////////////////////
			// Correct for deflection of the vertical
			// 1. compute bearing from estimated coordinates
			double azimuth(Direction(
				estimatedStations->get(stn1, 0),		// X1
				estimatedStations->get(stn1+1, 0),		// Y1
				estimatedStations->get(stn1+2, 0),		// Z1
				estimatedStations->get(stn2, 0), 		// X2
				estimatedStations->get(stn2+1, 0),		// Y2
				estimatedStations->get(stn2+2, 0),		// Z2
				stn1_it->currentLatitude,
				stn1_it->currentLongitude));
		
			// 2. Compute correction
			(*_it_msr)->preAdjCorr = ZenithDeflectionCorrection<double>(		// Correction to zenith distance for deflection of vertical
				azimuth,														// geodetic azimuth
				stn1_it->verticalDef,											// deflection in prime vertical
				stn1_it->meridianDef);											// deflection in prime meridian

			(*_it_msr)->term1 += (*_it_msr)->preAdjCorr;						// apply deflection correction
			////////////////////////////////////////////////////////////////////////////
		}
		else
			(*_it_msr)->preAdjCorr = 0.0;

		if (projectSettings_.a.stage)
			return;
	}

	// compute zenith distance from estimated coordinates
	double comp_msr(ZenithDistance(
		estimatedStations->get(stn1, 0),					// X1
		estimatedStations->get(stn1+1, 0),					// Y1
		estimatedStations->get(stn1+2, 0),					// Z1
		estimatedStations->get(stn2, 0), 					// X2
		estimatedStations->get(stn2+1, 0),					// Y2
		estimatedStations->get(stn2+2, 0),					// Z2
		stn1_it->currentLatitude,
		stn1_it->currentLongitude,
		stn2_it->currentLatitude,
		stn2_it->currentLongitude,
		(*_it_msr)->term3,									// instrument height
		(*_it_msr)->term4,									// target height
		&local_12e,											// local_12e, ..12n, ..12up represent
		&local_12n,											// the geometric difference between
		&local_12up));										// station1 and station2
	
	// Update measured minus computed value
	AddMsrtoMeasMinusComp(_it_msr, design_row, comp_msr, measMinusComp);

	// compute partial derivatives for normals
	double e2n2(local_12e*local_12e + local_12n*local_12n);
	double sqrt_e2n2(sqrt(e2n2));
	double se2n2_up2(sqrt_e2n2/(local_12up*local_12up));
	double up_se2n2(local_12up*sqrt_e2n2);
	double cos2v(cos(comp_msr)*cos(comp_msr));
	double cos_lat(cos(stn1_it->currentLatitude));
	double sin_lat(sin(stn1_it->currentLatitude));
	double cos_long(cos(stn1_it->currentLongitude));
	double sin_long(sin(stn1_it->currentLongitude));

	// Add partial derivatives dA/dX1, dA/dY1, dA/dZ1 to design matrix
	AddMsrtoDesign_BCEKMSVZ(design_row, stn1, stn2,
		cos2v * (((local_12e * sin_long + local_12n * sin_lat * cos_long) / up_se2n2) + cos_lat * cos_long * se2n2_up2),
		cos2v * (((-local_12e * cos_long + local_12n * sin_lat * sin_long) / up_se2n2) + cos_lat * sin_long * se2n2_up2),	
		cos2v * ((-local_12n * cos_lat / up_se2n2) + sin_lat * se2n2_up2),
		design);

	// Update AtVinv based on new design matrix elements
	UpdateAtVinv(_it_msr, stn1, stn2, 0, design_row, design, AtVinv, buildnewMatrices);

	if (buildnewMatrices)
		// Add weighted measurement contributions to normal matrix
		UpdateNormals_BCEKLMSVZ(stn1, stn2, design_row, normals, design, AtVinv);
	else
		design_row++;
}


// Vertical angles
//
// Like zenith distances, the relationship between slope distances and the
// coordinates of p1 and p2 requires a little extra work to take into consideration instrument
// and target height.  For the measurements-minus-computed vector, the "computed" distance
// is the true distance between the instrument and target, and so must take into consideration 
// instrument and target heights.  However, the dX, dY, dZ components for the partial 
// derivatives represent the true geometric difference between the two stations (not 
// instrument and target).
void dna_adjust::UpdateDesignNormalMeasMatrices_Z(pit_vmsr_t _it_msr, UINT32& design_row, const UINT32& block,
											  matrix_2d* measMinusComp, matrix_2d* estimatedStations, 
											  matrix_2d* normals, matrix_2d* design, matrix_2d* AtVinv, bool buildnewMatrices)
{
	UINT32 stn1(GetBlkMatrixElemStn1(block, _it_msr)); 
	UINT32 stn2(GetBlkMatrixElemStn2(block, _it_msr));

	it_vstn_t_const stn1_it(bstBinaryRecords_.begin() + (*_it_msr)->station1);
	it_vstn_t_const stn2_it(bstBinaryRecords_.begin() + (*_it_msr)->station2);

	double local_12e, local_12n, local_12up;

	// Initialise measurement.  No need to test if no further calculations are 
	// required (as in stage mode), as this is done later (below)
	InitialiseMeasurement(_it_msr, buildnewMatrices);

	if (buildnewMatrices)
	{
		// deflections available?
		if (fabs(stn1_it->verticalDef) > E4_SEC_DEFLECTION || fabs(stn1_it->meridianDef) > E4_SEC_DEFLECTION)
		{
			////////////////////////////////////////////////////////////////////////////
			// Correct for deflection of the vertical
			// 1. compute bearing from estimated coordinates
			double azimuth(Direction(
				estimatedStations->get(stn1, 0),		// X1
				estimatedStations->get(stn1+1, 0),		// Y1
				estimatedStations->get(stn1+2, 0),		// Z1
				estimatedStations->get(stn2, 0), 		// X2
				estimatedStations->get(stn2+1, 0),		// Y2
				estimatedStations->get(stn2+2, 0),		// Z2
				stn1_it->currentLatitude,
				stn1_it->currentLongitude));
		
			// 2. Compute correction
			(*_it_msr)->preAdjCorr = ZenithDeflectionCorrection<double>(		// Correction to vertical angle for deflection of vertical
				azimuth,														// geodetic azimuth
				stn1_it->verticalDef,											// deflection in prime vertical
				stn1_it->meridianDef);											// deflection in prime meridian
			
			(*_it_msr)->term1 -= (*_it_msr)->preAdjCorr;						// apply deflection correction
			////////////////////////////////////////////////////////////////////////////
		}
		else
			(*_it_msr)->preAdjCorr = 0.0;

		if (projectSettings_.a.stage)
			return;
	}

	// compute vertical angle from estimated coordinates
	double comp_msr(VerticalAngle(
		estimatedStations->get(stn1, 0),					// X1
		estimatedStations->get(stn1+1, 0),					// Y1
		estimatedStations->get(stn1+2, 0),					// Z1
		estimatedStations->get(stn2, 0), 					// X2
		estimatedStations->get(stn2+1, 0),					// Y2
		estimatedStations->get(stn2+2, 0),					// Z2
		stn1_it->currentLatitude,
		stn1_it->currentLongitude,
		stn2_it->currentLatitude,
		stn2_it->currentLongitude,
		(*_it_msr)->term3,									// instrument height
		(*_it_msr)->term4,									// target height
		&local_12e,											// local_12e, ..12n, ..12up represent
		&local_12n,											// the geometric difference between
		&local_12up));										// station1 and station2

	// Update measured minus computed value
	AddMsrtoMeasMinusComp(_it_msr, design_row, comp_msr, measMinusComp);

	// compute partial derivative terms
	double e2n2(local_12e*local_12e + local_12n*local_12n);
	double sqrt_e2n2(sqrt(e2n2));
	double se2n2_d_e2n2(sqrt_e2n2/e2n2);
	double up_d_se2n2_e2n2(local_12up/(sqrt_e2n2*e2n2));
	double cos2v(cos(comp_msr)*cos(comp_msr));
	double cos_lat(cos(stn1_it->currentLatitude));
	double sin_lat(sin(stn1_it->currentLatitude));
	double cos_long(cos(stn1_it->currentLongitude));
	double sin_long(sin(stn1_it->currentLongitude));

	// Add partial derivatives dA/dX1, dA/dY1, dA/dZ1 to design matrix
	AddMsrtoDesign_BCEKMSVZ(design_row, stn1, stn2,
		cos2v * ((-cos_lat * cos_long * se2n2_d_e2n2) - ((local_12e * sin_long + local_12n * sin_lat * cos_long) * up_d_se2n2_e2n2)),
		cos2v * ((-cos_lat * sin_long * se2n2_d_e2n2) + ((local_12e * cos_long - local_12n * sin_lat * sin_long) * up_d_se2n2_e2n2)),
		cos2v * ((-sin_lat * se2n2_d_e2n2) + (local_12n * cos_lat * up_d_se2n2_e2n2)),
		design);
	
	// Update AtVinv based on new design matrix elements
	UpdateAtVinv(_it_msr, stn1, stn2, 0, design_row, design, AtVinv, buildnewMatrices);

	if (buildnewMatrices)
		// Add weighted measurement contributions to normal matrix
		UpdateNormals_BCEKLMSVZ(stn1, stn2, design_row, normals, design, AtVinv);
	else
		design_row++;
}
	

// Orthometric height difference.
// This function requires geoid-ellipsoid separation in order to reduce 
// to ellipsoidal height difference counterpart.  
// Hence, run geoid with -f, -s and -n options
void dna_adjust::UpdateDesignNormalMeasMatrices_L(pit_vmsr_t _it_msr, UINT32& design_row, const UINT32& block,
											  matrix_2d* measMinusComp, matrix_2d* estimatedStations, 
											  matrix_2d* normals, matrix_2d* design, matrix_2d* AtVinv, bool buildnewMatrices)
{
	UINT32 stn1(GetBlkMatrixElemStn1(block, _it_msr)); 
	UINT32 stn2(GetBlkMatrixElemStn2(block, _it_msr));

	it_vstn_t_const stn1_it(bstBinaryRecords_.begin() + (*_it_msr)->station1);
	it_vstn_t_const stn2_it(bstBinaryRecords_.begin() + (*_it_msr)->station2);

	// Zn is the z coordinate element of the point on the z-axis 
	// which intersects with the the normal at the given Latitude
	double h1, h2, nu1, nu2, Zn1, Zn2;

	// calculated diff height
	double comp_msr(EllipsoidHeightDifference<double>(
		estimatedStations->get(stn1, 0), 
		estimatedStations->get(stn1+1, 0),
		estimatedStations->get(stn1+2, 0),
		estimatedStations->get(stn2, 0),
		estimatedStations->get(stn2+1, 0),
		estimatedStations->get(stn2+2, 0),
		stn1_it->currentLatitude,
		stn2_it->currentLatitude,
		&h1, &h2, &nu1, &nu2, &Zn1, &Zn2, 
		datum_.GetEllipsoidRef()));

	// Initialise measurement.  No need to test if no further calculations are 
	// required (as in stage mode), as this is done later (below)
	InitialiseMeasurement(_it_msr, buildnewMatrices);

	if (buildnewMatrices)
	{
		// N value available?
		if (fabs(stn1_it->geoidSep) > PRECISION_1E4 ||
			fabs(stn2_it->geoidSep) > PRECISION_1E4)
		{
			// Reduce to ellipsoidal difference
			(*_it_msr)->preAdjCorr = stn2_it->geoidSep - stn1_it->geoidSep;
			(*_it_msr)->term1 += (*_it_msr)->preAdjCorr;						
		}			

		if (projectSettings_.a.stage)
			return;
	}

	// Update measured minus computed value
	AddMsrtoMeasMinusComp(_it_msr, design_row, comp_msr, measMinusComp);

	// Add partial derivatives dA/dX1, dA/dY1, dA/dZ1, dA/dX2, dA/dY2, dA/dZ2 to design matrix
	AddMsrtoDesign_L(design_row, stn1, stn2,
		-estimatedStations->get(stn1, 0)/(nu1+h1), 
		-estimatedStations->get(stn1+1, 0)/(nu1+h1),
		-(estimatedStations->get(stn1+2, 0)+Zn1)/(nu1+h1),
		estimatedStations->get(stn2, 0)/(nu2+h2),
		estimatedStations->get(stn2+1, 0)/(nu2+h2),
		(estimatedStations->get(stn2+2, 0)+Zn2)/(nu2+h2),
		design);

	// Update AtVinv based on new design matrix elements
	UpdateAtVinv(_it_msr, stn1, stn2, 0, design_row, design, AtVinv, buildnewMatrices);

	if (buildnewMatrices)
		// Add weighted measurement contributions to normal matrix
		UpdateNormals_BCEKLMSVZ(stn1, stn2, design_row, normals, design, AtVinv);
	else
		design_row++;
}

void dna_adjust::UpdateDesignNormalMeasMatrices_I(pit_vmsr_t _it_msr, UINT32& design_row, const UINT32& block,
											  matrix_2d* measMinusComp, matrix_2d* estimatedStations, 
											  matrix_2d* normals, matrix_2d* design, matrix_2d* AtVinv, bool buildnewMatrices)
{
	// Initialise measurement.  No need to test if no further calculations are 
	// required (as in stage mode), as this is done later (below)
	InitialiseMeasurement(_it_msr, buildnewMatrices);

	if (buildnewMatrices)
	{
		it_vstn_t_const stn1_it(bstBinaryRecords_.begin() + (*_it_msr)->station1);
	
		// deflections available?
		if (fabs(stn1_it->meridianDef) > E4_SEC_DEFLECTION)
		{
			(*_it_msr)->preAdjCorr = stn1_it->meridianDef;						// deflection in the prime meridian
			(*_it_msr)->term1 -= (*_it_msr)->preAdjCorr;						// apply deflection correction
		}
		else
			(*_it_msr)->preAdjCorr = 0.0;

		if (projectSettings_.a.stage)
			return;
	}    

	UpdateDesignNormalMeasMatrices_IP(_it_msr, design_row, block,
		measMinusComp, estimatedStations, normals, design, AtVinv, buildnewMatrices);
}


void dna_adjust::UpdateDesignNormalMeasMatrices_J(pit_vmsr_t _it_msr, UINT32& design_row, const UINT32& block,
											  matrix_2d* measMinusComp, matrix_2d* estimatedStations, 
											  matrix_2d* normals, matrix_2d* design, matrix_2d* AtVinv, bool buildnewMatrices)
{
	// Initialise measurement.  No need to test if no further calculations are 
	// required (as in stage mode), as this is done later (below)
	InitialiseMeasurement(_it_msr, buildnewMatrices);

	if (buildnewMatrices)
	{
		it_vstn_t_const stn1_it(bstBinaryRecords_.begin() + (*_it_msr)->station1);
	
		// deflections available?
		if (fabs(stn1_it->verticalDef) > E4_SEC_DEFLECTION)
		{
			(*_it_msr)->preAdjCorr =										// deflection in the prime vertical
				stn1_it->verticalDef / cos(stn1_it->currentLatitude);		// sec(a) = 1/cos(a)
			(*_it_msr)->term1 -= (*_it_msr)->preAdjCorr;					// apply deflection correction
		}
		else
			(*_it_msr)->preAdjCorr = 0.0;

		if (projectSettings_.a.stage)
			return;
	}    

	UpdateDesignNormalMeasMatrices_JQ(_it_msr, design_row, block,
		measMinusComp, estimatedStations, normals, design, AtVinv, buildnewMatrices);
}


void dna_adjust::UpdateDesignNormalMeasMatrices_P(pit_vmsr_t _it_msr, UINT32& design_row, const UINT32& block,
											  matrix_2d* measMinusComp, matrix_2d* estimatedStations, 
											  matrix_2d* normals, matrix_2d* design, matrix_2d* AtVinv, bool buildnewMatrices)
{
	// Initialise measurement and test if no further calculations are 
	// required (as in stage mode)
	if (InitialiseMeasurement(_it_msr, buildnewMatrices))
		return;

	UpdateDesignNormalMeasMatrices_IP(_it_msr, design_row, block,
		measMinusComp, estimatedStations, normals, design, AtVinv, buildnewMatrices);
}
	

void dna_adjust::UpdateDesignNormalMeasMatrices_IP(pit_vmsr_t _it_msr, UINT32& design_row, const UINT32& block,
											  matrix_2d* measMinusComp, matrix_2d* estimatedStations, 
											  matrix_2d* normals, matrix_2d* design, matrix_2d* AtVinv, bool buildnewMatrices)
{
	UINT32 stn1(GetBlkMatrixElemStn1(block, _it_msr)); 
	
	// Due to the lack of an explicit relationship between latitude and cartesian elements,
	// use mechanical differentiation.
	double latitude;
	// partial derivative for x
	double partialDerivative = PartialD_Latitude_F(
		estimatedStations->get(stn1, 0),		// X1
		estimatedStations->get(stn1+1, 0),		// Y1
		estimatedStations->get(stn1+2, 0),		// Z1
		x_element,								// compute for X
		&latitude, 
		datum_.GetEllipsoidRef());
	AddElementtoDesign(design_row, stn1, partialDerivative, 
		design);								// X

	// Update measured minus computed value
	AddMsrtoMeasMinusComp(_it_msr, design_row, latitude, measMinusComp);

	// partial derivative for Y
	partialDerivative = PartialD_Latitude(
		estimatedStations->get(stn1, 0),		// X1
		estimatedStations->get(stn1+1, 0),		// Y1
		estimatedStations->get(stn1+2, 0),		// Z1
		y_element,								// compute for Y
		latitude, 
		datum_.GetEllipsoidRef());
	AddElementtoDesign(design_row, stn1+1, partialDerivative, 
		design);								// Y

	// partial derivative for Z
	partialDerivative = PartialD_Latitude(
		estimatedStations->get(stn1, 0),		// X1
		estimatedStations->get(stn1+1, 0),		// Y1
		estimatedStations->get(stn1+2, 0),		// Z1
		z_element,								// compute for Z
		latitude, 
		datum_.GetEllipsoidRef());
	AddElementtoDesign(design_row, stn1+2, partialDerivative, 
		design);								// Z

	// Update AtVinv based on new design matrix elements
	UpdateAtVinv(_it_msr, stn1, 0, 0, design_row, design, AtVinv, buildnewMatrices);
	
	if (buildnewMatrices)
		// Add weighted measurement contributions to normal matrix
		UpdateNormals_HIJPQR(stn1, design_row, normals, design, AtVinv);
	else
		design_row++;
}
	

void dna_adjust::UpdateDesignNormalMeasMatrices_Q(pit_vmsr_t _it_msr, UINT32& design_row, const UINT32& block,
											  matrix_2d* measMinusComp, matrix_2d* estimatedStations, 
											  matrix_2d* normals, matrix_2d* design, matrix_2d* AtVinv, bool buildnewMatrices)
{
	// Initialise measurement and test if no further calculations are 
	// required (as in stage mode)
	if (InitialiseMeasurement(_it_msr, buildnewMatrices))
		return;

	UpdateDesignNormalMeasMatrices_JQ(_it_msr, design_row, block,
		measMinusComp, estimatedStations, normals, design, AtVinv, buildnewMatrices);
}
	

void dna_adjust::UpdateDesignNormalMeasMatrices_JQ(pit_vmsr_t _it_msr, UINT32& design_row, const UINT32& block,
											  matrix_2d* measMinusComp, matrix_2d* estimatedStations, 
											  matrix_2d* normals, matrix_2d* design, matrix_2d* AtVinv, bool buildnewMatrices)
{
	UINT32 stn1(GetBlkMatrixElemStn1(block, _it_msr)); 
	
	it_vstn_t_const stn1_it(bstBinaryRecords_.begin() + (*_it_msr)->station1);
	
	// Add measured minus computed value
	AddMsrtoMeasMinusComp(_it_msr, design_row, stn1_it->currentLongitude, measMinusComp);

	// xy / (x^2 + y^2)^1.5
	double term1(estimatedStations->get(stn1, 0) * estimatedStations->get(stn1+1, 0) / 
		pow(estimatedStations->get(stn1, 0) * estimatedStations->get(stn1, 0) +
		estimatedStations->get(stn1+1, 0) * estimatedStations->get(stn1+1, 0), 1.5));

	// design matrix dlambda/dX1
	AddElementtoDesign(design_row, stn1, term1 * -1. / cos(stn1_it->currentLongitude), design);	// X1

	// design matrix dlambda/dY1
	AddElementtoDesign(design_row, stn1+1, term1 / sin(stn1_it->currentLongitude), design);		// Y1

	// design matrix dlambda/dZ1
	// the differentiation dlambda/dZ1 is undefined as longitude is not related
	// to the Z-coordinate. So use zero.
	AddElementtoDesign(design_row, stn1+2, 0., design);											// Z1
	
	// Update AtVinv based on new design matrix elements
	UpdateAtVinv(_it_msr, stn1, 0, 0, design_row, design, AtVinv, buildnewMatrices);
	
	if (buildnewMatrices)
		// Add weighted measurement contributions to normal matrix
		UpdateNormals_HIJPQR(stn1, design_row, normals, design, AtVinv);
	else
		design_row++;
}
	

void dna_adjust::UpdateDesignNormalMeasMatrices_H(pit_vmsr_t _it_msr, UINT32& design_row, const UINT32& block,
											  matrix_2d* measMinusComp, matrix_2d* estimatedStations, 
											  matrix_2d* normals, matrix_2d* design, matrix_2d* AtVinv, bool buildnewMatrices)
{
	// Initialise measurement.  No need to test if no further calculations are 
	// required (as in stage mode), as this is done later (below)
	InitialiseMeasurement(_it_msr, buildnewMatrices);

	if (buildnewMatrices)
	{
		it_vstn_t_const stn1_it(bstBinaryRecords_.begin() + (*_it_msr)->station1);
	
		// N value available?
		if (fabs(stn1_it->geoidSep) > PRECISION_1E4)
		{
			// reduce to ellipsoidal height
			(*_it_msr)->preAdjCorr = stn1_it->geoidSep;
			(*_it_msr)->term1 += (*_it_msr)->preAdjCorr;
		}

		if (projectSettings_.a.stage)
			return;
	}

	UpdateDesignNormalMeasMatrices_HR(_it_msr, design_row, block,
		measMinusComp, estimatedStations, normals, design, AtVinv, buildnewMatrices);
}


void dna_adjust::UpdateDesignNormalMeasMatrices_HR(pit_vmsr_t _it_msr, UINT32& design_row, const UINT32& block,
											  matrix_2d* measMinusComp, matrix_2d* estimatedStations, 
											  matrix_2d* normals, matrix_2d* design, matrix_2d* AtVinv, bool buildnewMatrices)
{
	UINT32 stn1(GetBlkMatrixElemStn1(block, _it_msr));

	it_vstn_t_const stn1_it(bstBinaryRecords_.begin() + (*_it_msr)->station1);
	
	// Zn is the z coordinate element of the point on the z-axis 
	// which intersects with the the normal at the given Latitude
	double nu1, Zn1;

	// calculated height
	double comp_msr(EllipsoidHeight<double>(
		estimatedStations->get(stn1, 0), 
		estimatedStations->get(stn1+1, 0),
		estimatedStations->get(stn1+2, 0),
		stn1_it->currentLatitude,
		&nu1, &Zn1, 
		datum_.GetEllipsoidRef()));
	
	// Add measured minus computed value
	AddMsrtoMeasMinusComp(_it_msr, design_row, comp_msr, measMinusComp);

	// design matrix dR/dX1
	AddElementtoDesign(design_row, stn1, estimatedStations->get(stn1, 0)/(nu1+comp_msr), design);			// X1

	// design matrix dR/dY1
	AddElementtoDesign(design_row, stn1+1, estimatedStations->get(stn1+1, 0)/(nu1+comp_msr), design);		// Y1

	// design matrix dR/dZ1
	AddElementtoDesign(design_row, stn1+2, (estimatedStations->get(stn1+2, 0)+Zn1)/(nu1+comp_msr), design);	// Z1
	
	// Update AtVinv based on new design matrix elements
	UpdateAtVinv(_it_msr, stn1, 0, 0, design_row, design, AtVinv, buildnewMatrices);
	
	if (buildnewMatrices)
		// Add weighted measurement contributions to normal matrix
		UpdateNormals_HIJPQR(stn1, design_row, normals, design, AtVinv);
	else
		design_row++;
}
	

void dna_adjust::UpdateDesignNormalMeasMatrices_R(pit_vmsr_t _it_msr, UINT32& design_row, const UINT32& block,
											  matrix_2d* measMinusComp, matrix_2d* estimatedStations, 
											  matrix_2d* normals, matrix_2d* design, matrix_2d* AtVinv, bool buildnewMatrices)
{
	// Initialise measurement and test if no further calculations are 
	// required (as in stage mode)
	if (InitialiseMeasurement(_it_msr, buildnewMatrices))
		return;

	UpdateDesignNormalMeasMatrices_HR(_it_msr, design_row, block,
		measMinusComp, estimatedStations, normals, design, AtVinv, buildnewMatrices);
}
	

void dna_adjust::UpdateDesignNormalMeasMatrices_X(pit_vmsr_t _it_msr, UINT32& design_row, const UINT32& block,
											  matrix_2d* measMinusComp, matrix_2d* estimatedStations, 
											  matrix_2d* design, matrix_2d* AtVinv, bool buildnewMatrices)
{
	it_vmsr_t _it_msr_first(*_it_msr);

	UINT32 stn1, stn2, covr(0), covc(0);
	
	UINT32 cluster_bsl, baseline_count((*_it_msr)->vectorCount1);
	UINT32 covariance_count;
	UINT32 design_row_begin(design_row);
		
	for (cluster_bsl=0; cluster_bsl<baseline_count; ++cluster_bsl)	// number of baselines/points
	{
		stn1 = GetBlkMatrixElemStn1(block, _it_msr); 
		stn2 = GetBlkMatrixElemStn2(block, _it_msr);

		UpdateDesignMeasMatrices_GX(_it_msr, design_row,
			measMinusComp, estimatedStations, design,
			stn1, stn2, buildnewMatrices);

		// UpdateDesignMeasMatrices_GX increments design_row by 2 internally (for X->Y and Y->Z)
		// We need one more increment to account for the Z component
		design_row++;
		
		covariance_count = (*_it_msr)->vectorCount2;
				
		// skip covariances until next point
		(*_it_msr) += covariance_count * 3;
		
		if (covariance_count > 0)
			(*_it_msr)++;

		if (projectSettings_.g.verbose > 5)
			debug_file << std::endl;
	}

	matrix_2d var_cart(baseline_count * 3, baseline_count * 3);

	if (buildnewMatrices || projectSettings_.a.stage)
	{
		// Load apriori variance matrix, and assign to binary measurement
		// If required, apply scalars
		LoadVarianceMatrix_X(_it_msr_first, &var_cart);

		if (buildnewMatrices && projectSettings_.a.stage)
			return;
	}

	// If this method is called via PrepareAdjustment() and the adjustment 
	// mode is staged, then don't update the AtVinv matrix.  This will be
	// done during an adjustment via AdjustPhasedForward().
	if (!buildnewMatrices && !projectSettings_.a.stage)
		return;

	UINT32 cluster_cov;

	it_vmsr_t _it_msr_temp(_it_msr_first);
	vUINT32 baseline_stations;

	baseline_stations.push_back(GetBlkMatrixElemStn1(block, &_it_msr_temp));
	baseline_stations.push_back(GetBlkMatrixElemStn2(block, &_it_msr_temp));
	std::sort(baseline_stations.begin(), baseline_stations.end());

	covc = baseline_count * 3;
	matrix_2d tmp(3, covc);

	// Build  At * V-1
	for (cluster_bsl=0; cluster_bsl<baseline_count; ++cluster_bsl)
	{
		stn1 = GetBlkMatrixElemStn1(block, &_it_msr_temp);
		stn2 = GetBlkMatrixElemStn2(block, &_it_msr_temp);
	
		// create a unique list of stations used in this cluster
		if (!binary_search(baseline_stations.begin(), baseline_stations.end(), stn1))
		{
			baseline_stations.push_back(stn1);
			std::sort(baseline_stations.begin(), baseline_stations.end());
		}

		if (!binary_search(baseline_stations.begin(), baseline_stations.end(), stn2))
		{
			baseline_stations.push_back(stn2);
			std::sort(baseline_stations.begin(), baseline_stations.end());
		}

		covr = cluster_bsl * 3;
		//covc = baseline_count * 3;

		// add variances/covariances for this baseline to At * V-1
		//AtVinv->blockadd(stn1, design_row_begin, 
		//	var_cart.submatrix(covr, 
		//	0, 3, covc).scale(-1.), 0, 0, 3, covc);						// add entire row of this vcv
		
		// copy the relevant var_cart elements to temp and multiply 
		// by -1 to effect the design matrix elements
		var_cart.submatrix(covr, 0, &tmp, 3, covc);
		tmp.scale(-1.);
		// add variances/covariances for this baseline to At * V-1
		AtVinv->blockadd(stn1, design_row_begin, 
			tmp, 0, 0, 3, covc);										// add entire row of this vcv		
		AtVinv->blockadd(stn2, design_row_begin, 
			var_cart, covr, 0, 3, covc);								// add entire row of this vcv
		
		covariance_count = _it_msr_temp->vectorCount2;
		_it_msr_temp += 3;			// move to covariances
		
		for (cluster_cov=0; cluster_cov<covariance_count; ++cluster_cov)
			_it_msr_temp += 3;
	}

	if (projectSettings_.a.stage)
		return;
	
	_it_msr_temp = _it_msr_first;
	
	matrix_2d tmp0(3, 3);

	// Build  At * V-1 * A variances
	for (cluster_bsl=0; cluster_bsl<baseline_count; ++cluster_bsl)
	{
		stn1 = GetBlkMatrixElemStn1(block, &_it_msr_temp);
		stn2 = GetBlkMatrixElemStn2(block, &_it_msr_temp);
	
		covr = cluster_bsl * 3;

		//// add variances for these stations
		//v_normals_.at(block).blockadd(stn1, stn1,						// Station 1.  This variance matrix
		//	AtVinv->submatrix(stn1,										// component must be multiplied by -1 to
		//		design_row_begin+covr, 3, 3).scale(-1.), 0, 0, 3, 3);	// effect the design matrix elements
		
		// copy the relevant var_cart elements to temp and multiply 
		// by -1 to effect the design matrix elements
		AtVinv->submatrix(stn1, design_row_begin+covr, &tmp0, 3, 3);
		tmp0.scale(-1.);

		// add variances for these stations
		v_normals_.at(block).blockadd(stn1, stn1,						// Station 1.
			tmp0, 0, 0, 3, 3);
		v_normals_.at(block).blockadd(stn2, stn2,						// Station 2
			v_AtVinv_.at(block), stn2, design_row_begin+covr, 3, 3);	

		covariance_count = _it_msr_temp->vectorCount2;
		_it_msr_temp += 3;			// move to covariances

		// add covariances for these stations
		for (cluster_cov=0; cluster_cov<covariance_count; ++cluster_cov)
			_it_msr_temp += 3;
	}

	matrix_2d tmp1(3, baseline_count*3), tmp2(baseline_count*3, 3);
	tmp2.zero();

	// Build  At * V-1 * A covariances
	for (cluster_bsl=0; cluster_bsl<baseline_stations.size(); ++cluster_bsl)
	{
		stn1 = baseline_stations.at(cluster_bsl);

		for (cluster_cov=0; cluster_cov<baseline_stations.size(); ++cluster_cov)
		{
			if (cluster_cov==cluster_bsl)
				continue;

			stn2 = baseline_stations.at(cluster_cov);

			//// add covariances for these stations
			//v_normals_.at(block).blockadd(stn1, stn2,		
			//	AtVinv->submatrix(stn1,			
			//		design_row_begin, 3, baseline_count*3).multiply(
			//			design->submatrix(design_row_begin, stn2, 
			//				baseline_count*3, 3)), 0, 0, 3, 3);

			// The following is more code, but is more efficient as it
			// recycles matrix instances (as opposed to creating a
			// new instance for every call to submatrix and multiply)
			// Get AtVinv component
			AtVinv->submatrix(stn1, design_row_begin, 
				&tmp1, 3, baseline_count*3);
			// Get design component
			design->submatrix(design_row_begin, stn2, 
				&tmp2, baseline_count*3, 3);
			// multiply the components to form covariance
			// for stations stn1 and stn2
			//tmp0.multiply(tmp1, tmp2);
			tmp0.multiply(tmp1, "N", tmp2, "N");
			// Now add covariances to normals
			v_normals_.at(block).blockadd(stn1, stn2,		
				tmp0, 0, 0, 3, 3);
		}
	}
}


void dna_adjust::UpdateDesignNormalMeasMatrices_Y(pit_vmsr_t _it_msr, UINT32& design_row, const UINT32& block,
											  matrix_2d* measMinusComp, matrix_2d* estimatedStations, 
											  matrix_2d* design, matrix_2d* AtVinv, bool buildnewMatrices)
{
	it_vmsr_t _it_msr_first(*_it_msr);
	it_vmsr_t tmp_msr;

	UINT32 stn1, covr(0), covc(0);
	
	UINT32 cluster_pnt, point_count((*_it_msr)->vectorCount1);
	UINT32 covariance_count;
	UINT32 design_row_begin(design_row);
	
	// the following is required in the case when Y cluster coordinates are LLH
	double latitude, longitude, ellipsoidHeight, x, y, z;
		
	it_vstn_t stn1_it;

	_COORD_TYPE_ coordType(CDnaStation::GetCoordTypeC((*_it_msr)->coordType));

	for (cluster_pnt=0; cluster_pnt<point_count; ++cluster_pnt)
	{
		covr = cluster_pnt * 3;
		covariance_count = (*_it_msr)->vectorCount2;

		stn1_it = bstBinaryRecords_.begin() + (*_it_msr)->station1;
	
		tmp_msr = *_it_msr;

		stn1 = GetBlkMatrixElemStn1(block, _it_msr);

		// Convert to cartesian reference frame?
        	if (coordType == LLH_type_i || coordType == LLh_type_i)
		{
			// This section (and others in this function where coordType == LLH_type_i) 
			// will only be performed once because this measurement will be converted 
			// to XYZ and coordType will be set to XYZ_type_i

			// Get latitude value
			latitude = tmp_msr->term1;
			// retain original value
			tmp_msr->preAdjMeas = tmp_msr->term1;
			tmp_msr++;

			// Get longitude value
			longitude = tmp_msr->term1; 
			// retain original value
			tmp_msr->preAdjMeas = tmp_msr->term1;
			tmp_msr++;

			// Get height
			// NOTE: heights for Y clusters are orthometric.
			// Retain original value
			ellipsoidHeight = tmp_msr->term1;
			tmp_msr->preAdjMeas = tmp_msr->term1;

			// Reduce to ellipsoid height
			if (fabs(stn1_it->geoidSep) > PRECISION_1E4)
			{
				tmp_msr->preAdjCorr = stn1_it->geoidSep;
				ellipsoidHeight += tmp_msr->preAdjCorr;				
			}

			// convert
			GeoToCart<double>(latitude, longitude, ellipsoidHeight, 
				&x, &y, &z, 
				datum_.GetEllipsoidRef());

			// Update bms record
			(*_it_msr)->term1 = x;
			snprintf((*_it_msr)->coordType, sizeof((*_it_msr)->coordType), "%s", XYZ_type);

			// retain original reference frame
            		(*_it_msr)->station3 = coordType;
		}

		// If this method is called via PrepareAdjustment() and the adjustment 
		// mode is staged, then don't calculate measured-computed values or
		// fill design matrix.  This will be done during an adjustment
		// via AdjustPhasedForward().
		if (!(buildnewMatrices && projectSettings_.a.stage))
		{
			// For all adjustment modes, when this method is called during an adjustment
			// to update the normals, only the the measured-computed values need to be 
			// updated.  This is because the (Jacobian) design matrix elements for GPS
			// are unity (1 or -1) and do not change as coordinates are updated (unlike 
			// the Jacobian elements formed for other measurements).  Hence, for all 
			// adjustment modes except staged, the design matrix is updated once via
			// PrepareAdjustment.  For staged adjustments, the design matrix is updated
			// on each iteration.

			// Add X element to measured minus computed
			AddMsrtoMeasMinusComp(_it_msr, design_row, 
				estimatedStations->get(stn1, 0), 
				measMinusComp);
		
			if (buildnewMatrices || projectSettings_.a.stage)
				// Add X elements to design matrix
				AddElementtoDesign(design_row, stn1, 1., design);		// X
			
			design_row++;
		}

		if (buildnewMatrices)
			if (coordType != LLH_type_i && (*_it_msr)->station3 != LLH_type_i &&
                coordType != LLh_type_i && (*_it_msr)->station3 != LLh_type_i)
				(*_it_msr)->preAdjMeas = (*_it_msr)->term1;

		// move to Y element
		(*_it_msr)++;

		// measurements matrix Y
        	if (coordType == LLH_type_i || coordType == LLh_type_i)
		{
			// Update bms record
			(*_it_msr)->term1 = y;
			snprintf((*_it_msr)->coordType, sizeof((*_it_msr)->coordType), "%s", XYZ_type);

			// retain original reference frame
            		(*_it_msr)->station3 = coordType;
		}

		// If this method is called via PrepareAdjustment() and the adjustment 
		// mode is staged, then don't calculate measured-computed values or
		// fill design matrix.  This will be done during an adjustment
		// via AdjustPhasedForward().
		if (!(buildnewMatrices && projectSettings_.a.stage))
		{
			// For all adjustment modes, when this method is called during an adjustment
			// to update the normals, only the the measured-computed values need to be 
			// updated.  This is because the (Jacobian) design matrix elements for GPS
			// are unity (1 or -1) and do not change as coordinates are updated (unlike 
			// the Jacobian elements formed for other measurements).  Hence, for all 
			// adjustment modes except staged, the design matrix is updated once via
			// PrepareAdjustment.  For staged adjustments, the design matrix is updated
			// on each iteration.

			// Add Y element to measured minus computed
			AddMsrtoMeasMinusComp(_it_msr, design_row, 
				estimatedStations->get(stn1+1, 0),
				measMinusComp, false);
		
			if (buildnewMatrices || projectSettings_.a.stage)
				// Add Y elements to design matrix
				AddElementtoDesign(design_row, stn1+1, 1., design);		// Y

			design_row++;
		}

		if (buildnewMatrices)
			if ((coordType != LLH_type_i && (*_it_msr)->station3 != LLH_type_i) &&
				(coordType != LLh_type_i && (*_it_msr)->station3 != LLh_type_i))
				(*_it_msr)->preAdjMeas = (*_it_msr)->term1;

		// move to Z element
		(*_it_msr)++;

		// measurements matrix Z
        	if (coordType == LLH_type_i || coordType == LLh_type_i)
		{
			// Update bms record
			(*_it_msr)->term1 = z;
			snprintf((*_it_msr)->coordType, sizeof((*_it_msr)->coordType), "%s", XYZ_type);

			// retain original reference frame
			(*_it_msr)->station3 = coordType;
		}

		// If this method is called via PrepareAdjustment() and the adjustment 
		// mode is staged, then don't calculate measured-computed values or
		// fill design matrix.  This will be done during an adjustment
		// via AdjustPhasedForward().
		if (!(buildnewMatrices && projectSettings_.a.stage))
		{
			// For all adjustment modes, when this method is called during an adjustment
			// to update the normals, only the the measured-computed values need to be 
			// updated.  This is because the (Jacobian) design matrix elements for GPS
			// are unity (1 or -1) and do not change as coordinates are updated (unlike 
			// the Jacobian elements formed for other measurements).  Hence, for all 
			// adjustment modes except staged, the design matrix is updated once via
			// PrepareAdjustment.  For staged adjustments, the design matrix is updated
			// on each iteration.

			// Add Z element to measured minus computed
			AddMsrtoMeasMinusComp(_it_msr, design_row, 
				estimatedStations->get(stn1+2, 0),
				measMinusComp, false);
            if (projectSettings_.g.verbose > 5)
                debug_file << std::endl;
		
			if (buildnewMatrices || projectSettings_.a.stage)
				// Add Z elements to design matrix
				AddElementtoDesign(design_row, stn1+2, 1., design);		// Z

			design_row++;
		}

		if (buildnewMatrices)
			if ((coordType != LLH_type_i && (*_it_msr)->station3 != LLH_type_i) &&
				(coordType != LLh_type_i && (*_it_msr)->station3 != LLh_type_i))
				(*_it_msr)->preAdjMeas = (*_it_msr)->term1;
		
		covariance_count = (*_it_msr)->vectorCount2;
				
		// skip covariances until next point
		(*_it_msr) += covariance_count * 3;
		
		if (covariance_count > 0)
			(*_it_msr)++;
	}

	matrix_2d var_cart(point_count * 3, point_count * 3);

	if (buildnewMatrices || projectSettings_.a.stage)
	{
		// Load apriori variance matrix, and assign to binary measurement
		// If required, propagate to cartesian reference frame and apply 
		// scalars
		LoadVarianceMatrix_Y(_it_msr_first, &var_cart, coordType);
	
		// If preparing for a stage adjustment, return
		// Normals will be built for each block as needed
		if (buildnewMatrices && projectSettings_.a.stage)
			return;
	}
	
	// If this method is called via PrepareAdjustment() and the adjustment 
	// mode is staged, then don't update the AtVinv matrix.  This will be
	// done during an adjustment via AdjustPhasedForward().
	if (!buildnewMatrices && !projectSettings_.a.stage)
		return;

	UINT32 cluster_cov;

	UINT32 stn2, cov_c;
	it_vmsr_t _it_msr_temp(_it_msr_first);

	// Build  At * V-1
	for (cluster_pnt=0; cluster_pnt<point_count; ++cluster_pnt)
	{		
		stn1 = GetBlkMatrixElemStn1(block, &_it_msr_temp);
		covariance_count = _it_msr_temp->vectorCount2;
		
		covr = cluster_pnt * 3;
		covc = cluster_pnt * 3;

		// add variance for this station
		AtVinv->copyelements(stn1, design_row_begin, var_cart, covr, covc, 3, 3);

		if (covariance_count < 1)
			break;

		_it_msr_temp += 3;	// advance to first cluster point (if clusters exist, otherwise, next point)
		cov_c = 0;

		for (cluster_cov=0; cluster_cov<covariance_count; ++cluster_cov)		// number of baseline/point covariances
		{
			stn2 = GetBlkMatrixElemStn1(block, &_it_msr_temp);	// this becomes stn1 in next cluster_pnt loop

			covc += 3;
			cov_c += 3;

			// add covariance between stn1 and this station
			AtVinv->copyelements(stn1, design_row_begin + cov_c, var_cart, covr, covc, 3, 3);
			AtVinv->copyelements(stn2, design_row_begin, var_cart, covc, covr, 3, 3);

			_it_msr_temp += 3;
		}

		design_row_begin += 3;
	}

	if (projectSettings_.a.stage)
		return;
	
	_it_msr_temp = _it_msr_first;

	// Add to At * V-1 * A
	for (cluster_pnt=0; cluster_pnt<point_count; ++cluster_pnt)
	{		
		stn1 = GetBlkMatrixElemStn1(block, &_it_msr_temp);
		covariance_count = _it_msr_temp->vectorCount2;
		
		covr = cluster_pnt * 3;
		covc = cluster_pnt * 3;

		// add variance for this station
		v_normals_.at(block).blockadd(stn1, stn1, var_cart, covr, covc, 3, 3);

		if (covariance_count < 1)
			break;

		_it_msr_temp += 3;

		for (cluster_cov=0; cluster_cov<covariance_count; ++cluster_cov)	// number of baseline/point covariances
		{
			stn2 = GetBlkMatrixElemStn1(block, &_it_msr_temp);		// this becomes stn1 in next cluster_pnt loop
			covc += 3;

			if (stn1 < stn2)
			{
				// add covariance between stn1 and this station
				v_normals_.at(block).blockadd(stn1, stn2, var_cart, covr, covc, 3, 3);
				v_normals_.at(block).blockadd(stn2, stn1, var_cart, covc, covr, 3, 3);
			}
			else
			{
				// Transpose
				// add covariance between stn1 and this station
				v_normals_.at(block).blockTadd(stn1, stn2, var_cart, covr, covc, 3, 3);
				v_normals_.at(block).blockTadd(stn2, stn1, var_cart, covc, covr, 3, 3);
			}
			_it_msr_temp += 3;
		}
	}
}
	

void dna_adjust::SolveTry(bool COMPUTE_INVERSE, const UINT32& block)
{
	// Least Squares Solution
	try {            
		Solve(COMPUTE_INVERSE, block);
	}
	catch (const std::runtime_error& e) {

		// debug matrices if required
		debug_SolutionInformation(block);

		// Could not invert matrix.  Fire an exception
		SignalExceptionAdjustment(e.what(), block);
	}
}
	

void dna_adjust::Solve(bool COMPUTE_INVERSE, const UINT32& block)
{
	// debug matrices if required
	debug_SolutionInformation(block);

	if (COMPUTE_INVERSE)
	{
		// When non-GPS measurements exist, partial derivatives will vary upon
		// each iteration due to changes in the latest estimates, and so
		// normals will be reformed.  Hence, compute inverse.

		//////////////////
		// The procedure for constraining free coordinates causes the elements 
		// in the normal matrix to vary greatly in magnitude.  In some cases,
		// this opens the possibility of a loss of accuracy in the 
		// matrix inversion process as the algorithm will be repeatedly 
		// computing the difference between numbers of vastly different 
		// magnitude. Inaccuracy in the inversion can be mitigated by scaling 
		// the normal matrix before inversion and subsequently reversing the effect.
		//

		// 1. Create scalar matrix
		matrix_2d *S = nullptr, *SN = nullptr;

		if (projectSettings_.a.scale_normals_to_unity)
		{
			S = new matrix_2d(v_normals_.at(block).rows(), v_normals_.at(block).rows());
			SN = new matrix_2d(v_normals_.at(block).rows(), v_normals_.at(block).rows());
			for (UINT32 i(0); i<v_normals_.at(block).rows(); ++i)
				S->put(i, i, sqrt(v_normals_.at(block).get(i, i)));
			// 2. Scale Normals to reduce the diagonal elements of Normals to unity
			//SN->multiply(*S, v_normals_.at(block));
			SN->multiply(*S, "N", v_normals_.at(block), "N");
			
			//v_normals_.at(block).multiply(*SN, *S);
			v_normals_.at(block).multiply(*SN, "N", *S, "N");
		}
		//////////////////
	
		// Clear upper triangle to match expected structure for cholesky_inverse
		v_normals_.at(block).clearupper();

		// Calculate Inverse of AT * V-1 * A
		// Explicitly set LOWER_IS_CLEARED to false (data is in lower triangle)
		try {
			FormInverseVarianceMatrix(&(v_normals_.at(block)), false);
		}
	catch (const exception::MatrixInversionException& e) {
		std::stringstream enhanced_msg;
		enhanced_msg << "Solve(): Matrix inversion failed:\n" << e.what();

		UINT32 dpotrf_row = UINT_MAX;
		if (e.dpotrf_info() > 0) {
			dpotrf_row = e.dpotrf_info() - 1;
		}

		UINT32 gershgorin_row = UINT_MAX;
		if (e.gershgorin_critical_row() >= 0) {
			gershgorin_row = static_cast<UINT32>(e.gershgorin_critical_row());
		}

		auto report_station = [&](UINT32 row, const char* source) {
			UINT32 station_pos = row / 3;
			UINT32 param_idx = row % 3;
			const char* param_names[] = {"Latitude", "Longitude", "Height"};

			if (station_pos < v_parameterStationList_.at(block).size()) {
				UINT32 global_stn_idx = v_parameterStationList_.at(block).at(station_pos);
				const station_t& stn = bstBinaryRecords_.at(global_stn_idx);

				enhanced_msg << "\n\nProblematic Station (" << source << "):\n";
				enhanced_msg << "  Station name: " << stn.stationName << "\n";
				enhanced_msg << "  Constraint:   " << stn.stationConst << "\n";
				enhanced_msg << "  Coordinates:  " << std::fixed << std::setprecision(9);
				enhanced_msg << "Lat=" << stn.currentLatitude << ", ";
				enhanced_msg << "Lon=" << stn.currentLongitude << ", ";
				enhanced_msg << "Hgt=" << std::setprecision(4) << stn.currentHeight << "\n";
				enhanced_msg << "  Critical row: " << row << " (station " << station_pos
							 << ", parameter: " << param_names[param_idx] << ")\n";
			}
		};

		if (dpotrf_row != UINT_MAX) {
			report_station(dpotrf_row, "dpotrf failure point");
		}

		if (gershgorin_row != UINT_MAX && gershgorin_row != dpotrf_row) {
			report_station(gershgorin_row, "Gershgorin analysis");
		}

		std::set<UINT32> problematic_stations;
		if (e.has_eigendata()) {
			for (const auto& evd : e.eigenvectors()) {
				for (const auto& comp : evd.components) {
					UINT32 station_pos = comp.row / 3;
					if (station_pos < v_parameterStationList_.at(block).size()) {
						UINT32 global_stn_idx = v_parameterStationList_.at(block).at(station_pos);
						problematic_stations.insert(global_stn_idx);
					}
				}
			}
		}

		if (!problematic_stations.empty()) {
			std::vector<const measurement_t*> related_measurements;
			std::set<UINT32> added_measurements;

			for (UINT32 station_idx : problematic_stations) {
				if (station_idx >= vAssocStnList_.size())
					continue;

				const CAStationList& asl = vAssocStnList_.at(station_idx);
				UINT32 msr_count = asl.GetAssocMsrCount();
				UINT32 aml_start = asl.GetAMLStnIndex();

				for (UINT32 m = 0; m < msr_count; ++m) {
					UINT32 aml_index = aml_start + m;
					if (aml_index >= vAssocMsrList_.size())
						continue;

					UINT32 bmsr_index = vAssocMsrList_.at(aml_index).bmsr_index;
					if (bmsr_index >= bmsBinaryRecords_.size())
						continue;

					if (added_measurements.count(bmsr_index) > 0)
						continue;

					const measurement_t& msr = bmsBinaryRecords_.at(bmsr_index);
					if (!msr.ignore) {
						related_measurements.push_back(&msr);
						added_measurements.insert(bmsr_index);
					}
				}
			}

			if (!related_measurements.empty()) {
				enhanced_msg << "\n\nMeasurements involving problematic stations ("
							 << related_measurements.size() << " total):\n";
				enhanced_msg << "  " << std::left
							 << std::setw(6) << "Type"
							 << std::setw(8) << "FileOrd"
							 << std::setw(25) << "Station1"
							 << std::setw(25) << "Station2"
							 << std::setw(25) << "Station3"
							 << std::right
							 << std::setw(12) << "Value"
							 << std::setw(12) << "StdDev"
							 << std::setw(12) << "Residual"
							 << std::setw(10) << "NStat"
							 << std::setw(8) << "Cluster" << "\n";
				enhanced_msg << "  " << std::string(145, '-') << "\n";

				for (size_t i = 0; i < related_measurements.size(); ++i) {
					const measurement_t* msr = related_measurements[i];
					const char* stn1_name = bstBinaryRecords_.at(msr->station1).stationName;
					const char* stn2_name = (msr->measurementStations >= 2) ?
						bstBinaryRecords_.at(msr->station2).stationName : "-";
					const char* stn3_name = (msr->measurementStations >= 3) ?
						bstBinaryRecords_.at(msr->station3).stationName : "-";

					double stddev = sqrt(msr->term2);

					enhanced_msg << "  " << std::left
								 << std::setw(6) << msr->measType
								 << std::setw(8) << msr->fileOrder
								 << std::setw(25) << stn1_name
								 << std::setw(25) << stn2_name
								 << std::setw(25) << stn3_name
								 << std::right << std::fixed << std::setprecision(4)
								 << std::setw(12) << msr->term1
								 << std::setw(12) << stddev
								 << std::setw(12) << msr->measCorr
								 << std::setw(10) << std::setprecision(2) << msr->NStat
								 << std::setw(8) << msr->clusterID << "\n";
				}
			}
		}

		SignalExceptionAdjustment(enhanced_msg.str(), 0);
	}
	catch (const std::runtime_error& e) {
		std::stringstream ss;
		ss << "Solve(): Matrix inversion failed:\n" << e.what();
		SignalExceptionAdjustment(ss.str(), 0);
	}

		// Check for a failed inverse solution
		if (boost::math::isnan(v_normals_.at(block).get(0, 0)) ||
			boost::math::isinf(v_normals_.at(block).get(0, 0)))
		{
			std::stringstream ss;
			ss << "Solve(): Invalid variance matrix:" << std::endl;
			ss << std::setprecision(6) << std::fixed << v_normals_.at(block);
			SignalExceptionAdjustment(ss.str(), 0);
		}

		//////////////////
		// 2. Compute inverse of N (via S * (SNS)-1 * S)
		if (projectSettings_.a.scale_normals_to_unity)
		{
			//SN->multiply(*S, v_normals_.at(block));
			SN->multiply(*S, "N", v_normals_.at(block), "N");

			//v_normals_.at(block).multiply(*SN, *S);
			v_normals_.at(block).multiply(*SN, "N", *S, "N");

			delete S;
			delete SN;
		}
		//////////////////
	}
	
	if (projectSettings_.g.verbose > 0)
	{
		debug_file << "Block " << block + 1;
		if (projectSettings_.a.adjust_mode != SimultaneousMode)
			debug_file << (forward_ ? " (Forward)" : " (Reverse)");
		debug_file << std::endl;
		debug_file << "Precisions " << std::fixed << std::setprecision(16) << v_normals_.at(block) << std::endl;
	}
	
	// compute weighted "measured minus computed"
	matrix_2d At_Vinv_m(v_design_.at(block).columns(), 1);
	At_Vinv_m.multiply(v_AtVinv_.at(block), "N", v_measMinusComp_.at(block), "N");
	
	// Solve corrections from normal equations
	v_corrections_.at(block).redim(v_design_.at(block).columns(), 1);
	v_corrections_.at(block).multiply(v_normals_.at(block), "N", At_Vinv_m, "N");

	if (projectSettings_.g.verbose > 0)
	{

		if (projectSettings_.a.multi_thread)
			dbg_file_mutex.lock();

		it_vUINT32 _it_stn;
		UINT32 i;
		
		switch (projectSettings_.a.adjust_mode)
		{
		case Phased_Block_1Mode:
		case PhasedMode:
			for (i=0; i<v_parameterStationList_.at(block).size(); ++i)
			{
				debug_file << bstBinaryRecords_.at(v_parameterStationList_.at(block).at(i)).stationName;
				if ((_it_stn = find(v_ISL_.at(block).begin(), v_ISL_.at(block).end(), v_parameterStationList_.at(block).at(i))) != v_ISL_.at(block).end())
					debug_file << "i ";
				else if ((_it_stn = find(v_JSL_.at(block).begin(), v_JSL_.at(block).end(), v_parameterStationList_.at(block).at(i))) != v_JSL_.at(block).end())
						debug_file << "j ";
				else 
					debug_file << " ";				
			}
			break;
		case SimultaneousMode:
			for (i=0; i<v_ISL_.at(block).size(); ++i)
				debug_file << bstBinaryRecords_.at(v_ISL_.at(block).at(i)).stationName << " ";				
			break;
		}

		debug_file << std::endl;
		debug_file << "Block " << block + 1;
		if (projectSettings_.a.adjust_mode != SimultaneousMode)
			debug_file << (forward_ ? " (Forward)" : " (Reverse)");
		debug_file << std::endl;
		debug_file << "Weighted measurements " << std::fixed << std::setprecision(16) << At_Vinv_m << std::endl;

		debug_file << "Block " << block + 1;
		if (projectSettings_.a.adjust_mode != SimultaneousMode)
			debug_file << (forward_ ? " (Forward)" : " (Reverse)");
		debug_file << std::endl;
		debug_file << "Corrections " << std::fixed << std::setprecision(16) << v_corrections_.at(block) << std::endl;
		debug_file.flush();

		if (projectSettings_.a.multi_thread)
			dbg_file_mutex.unlock();

	}
}
	

void dna_adjust::DeSerialiseAdjustedVarianceMatrices()
{
	// No need to facilitate serialising if network adjustment is in stage,
	// as this will already be taken care of
	if (projectSettings_.a.stage)
		return;

	// Create file streams and memory map regions
	OpenStageFileStreams(2, sf_rigorous_vars, sf_prec_adj_msrs);
	ReserveBlockMapRegions(2, sf_rigorous_vars, sf_prec_adj_msrs);

	for (UINT32 block=0; block<blockCount_; ++block)
	{
		// Set Region offsets
		SetRegionOffsets(block, 2, sf_rigorous_vars, sf_prec_adj_msrs);

		// Prepare mapping regions
		PrepareMemoryMapRegions(block, 2, sf_rigorous_vars, sf_prec_adj_msrs);
	}

	
	try {
		// Set memory mapped file regions
		// NOTE - previously created files must exist and match the 
		// dimensions of the current matrix sizes			
		SetMapRegions(2, sf_rigorous_vars, sf_prec_adj_msrs);
		
	}
	catch (boost::interprocess::interprocess_exception& e){
		std::stringstream ss;
		ss << "DeSerialiseAdjustedVarianceMatrices() terminated while creating memory map" << std::endl;
		ss << "  regions from mtx file. Details:\n  " << e.what() << std::endl;
		adj_file << std::endl << "- Error: " << ss.str() << std::endl;
		SignalExceptionAdjustment(ss.str(), 0);
	}

	// Now load up matrix from .mtx file
	for (UINT32 block=0; block<blockCount_; ++block)
	{
		// Read this block
		DeserialiseBlockFromMappedFile(block, 2,
			sf_rigorous_vars, sf_prec_adj_msrs);
	}

	// Close rigorous variances stream
	f_rigorousVariances_.close();
	f_precAdjMsrs_.close();
}
	

void dna_adjust::SerialiseAdjustedVarianceMatrices()
{
	// No need to facilitate serialising if network adjustment is in stage,
	// as this will already be taken care of
	if (projectSettings_.a.stage)
		return;

	// Create file streams and memory map regions
	projectSettings_.a.recreate_stage_files = 1;
	OpenStageFileStreams(2, sf_rigorous_vars, sf_prec_adj_msrs);
	ReserveBlockMapRegions(2, sf_rigorous_vars, sf_prec_adj_msrs);

	// Now offload matrices to new .mtx stage files

	UINT32 block;
	
	for (block=0; block<blockCount_; ++block)
	{
		// Set memory mapped file region offsets
		SetRegionOffsets(block, 2, sf_rigorous_vars, sf_prec_adj_msrs);

		// Write to disk (the memory mapped file)
		f_rigorousVariances_ << v_rigorousVariances_.at(block);			// ...rva.mtx
		f_precAdjMsrs_ << v_precAdjMsrsFull_.at(block);					// ...pam.mtx
	}

	// Close rigorous variances stream
	f_rigorousVariances_.close();
	f_precAdjMsrs_.close();
}
	

void dna_adjust::GenerateStatistics()
{
	try {
		// Backup the latest estimates and, if this is not a staged 
		// adjustment, update normals and measured-computed matrices.
		UpdateAdjustment(false);

		// Compute whole-of-network statistics.
		ComputeStatistics();

		// Print statistics summary to adj file
		switch (projectSettings_.a.adjust_mode)
		{
		case Phased_Block_1Mode:
			printer_->PrintStatistics(false);
			break;
		//case SimultaneousMode:
		//case PhasedMode:
		default:
			printer_->PrintStatistics();
			break;
		}
	}
	catch (const std::out_of_range& e) {
		std::cerr << "ERROR in GenerateStatistics: out_of_range exception - " << e.what() << std::endl;
		throw;
	}
	catch (const std::exception& e) {
		std::cerr << "ERROR in GenerateStatistics: " << e.what() << std::endl;
		throw;
	}

	isAdjustmentQuestionable_ = (
		// non zero means something is amiss
		adjustStatus_ ||
		// a much larger than expected outcome
		sigmaZero_ > 10.0 * chiSquaredUpperLimit_ ||
		// didn't converge
		fabs(maxCorr_) > projectSettings_.a.iteration_threshold);
}
	

void dna_adjust::ComputeChiSquareSimultaneous()
{
	measurementCount_ = v_measurementCount_.at(0);
	measurementParams_ = v_measurementParams_.at(0);

	// Compute adjusted measurement statistics
	ComputeChiSquare(0);
}
	

void dna_adjust::ComputeGlobalNetStat()
{
	degreesofFreedom_ = measurementParams_ - unknownParams_;
	sigmaZero_ = sigmaZeroSqRt_ = 0.;
	if (degreesofFreedom_ != 0)
	{
		sigmaZero_ = chiSquared_ / degreesofFreedom_;
		sigmaZeroSqRt_ = sqrt(sigmaZero_);
	}
}
	

void dna_adjust::ComputeTestStat(const double& dof, double& chiUpper, double& chiLower, double& sigmaZero, UINT32& passFail)
{
	// Calculate limits using boost libraries
	// Confidence is halved (i.e. * 0.5) because it is a two-sided (upper and lower limit) test.
	double confidence = (100. - projectSettings_.a.confidence_interval) * 0.01;
	double conf = confidence * 0.5;

	try {
		// chi_squared throws when dof == 0
		boost::math::chi_squared dist(dof);
		
		chiUpper = boost::math::quantile(boost::math::complement(dist, conf)) / degreesofFreedom_;
		chiLower = boost::math::quantile(dist, conf) / degreesofFreedom_;

		switch (projectSettings_.a.adjust_mode)
		{
		case PhasedMode:
		case Phased_Block_1Mode:
			sigmaZero = chiSquared_ / dof;
			break;
		}

		// A pass is when sigma-zero is less than the chi-square upper limit.  Strictly speaking,
		// a value lower than the lower limit is not a fail.
		if (sigmaZero < chiLower)
			passFail = test_stat_warning;
		else if (sigmaZero > chiUpper)
			passFail = test_stat_fail;
		else
			passFail = test_stat_pass;
	}
	catch (const std::exception& e)
	{ 
		// This is not a critical error, so don't throw an exception.  Rather,
		// just print message to output files

		if (projectSettings_.g.verbose > 0)
			debug_file << std::endl << "ComputeTestStat():\n     " << e.what() << std::endl;

		if (projectSettings_.g.verbose > 0)
			if (dof == 0)
				adj_file << std::endl << "Cannot perform chi-square test with zero degrees of freedom." << std::endl << std::endl; 
	
		passFail = test_stat_fail;
	}
}
	

void dna_adjust::ComputeGlobalTestStat()
{
	ComputeTestStat(degreesofFreedom_, chiSquaredUpperLimit_, chiSquaredLowerLimit_, 
		sigmaZero_, passFail_);
}
	

void dna_adjust::ComputeBlockTestStat(const UINT32& block)
{
	UINT32 passFail(test_stat_pass);

	ComputeTestStat(
		v_statSummary_.at(block)._degreesofFreedom, 
		v_chiSquaredUpperLimit_.at(block), 
		v_chiSquaredLowerLimit_.at(block), 
		v_sigmaZero_.at(block), 
		passFail);

	// std::vector<bool> is a specialisation, and as such the value 
	// returned from at() is not a bool, but a wrapper around a value
	v_passFail_.at(block) = passFail;
}
	

void dna_adjust::ComputeAdjustedMsrPrecisions()
{
	// Compute adjusted measurement precisions
	UINT32 block;
	switch (projectSettings_.a.adjust_mode)
	{
	case PhasedMode:
		for (block=0; block<blockCount_; ++block)
		{
			// For staged adjustments, load block info
			if (projectSettings_.a.stage)
			{
				DeserialiseBlockFromMappedFile(block, 7,
					sf_normals, sf_rigorous_vars, 
					sf_design, sf_atvinv, sf_estimated_stns,
					sf_meas_minus_comp, sf_prec_adj_msrs);
				v_normals_.at(block) = v_rigorousVariances_.at(block);

				FillDesignNormalMeasurementsMatrices(false, block, false);
			}

			// Compute adjusted measurement precisions (v_precAdjMsrsFull_)
			// from design and rigorous station variances 
			ComputePrecisionAdjMsrs(block);

			// Update measurement records and Pelzer's Global reliability
			UpdateMsrRecords(block);

			// For staged adjustments, unload all matrix data
			if (projectSettings_.a.stage)
				UnloadBlock(block);
		}
		break;
	case SimultaneousMode:
	case Phased_Block_1Mode:			// only block 1 is rigorous
		// Compute adjusted measurement precisions (v_precAdjMsrsFull_)
		// from design and rigorous station variances 
		ComputePrecisionAdjMsrs();
		
		// Update measurement records and Pelzer's Global reliability
		UpdateMsrRecords();
		break;
	}
}
	
void dna_adjust::UpdateMsrTstatistic_D(it_vmsr_t& _it_msr)
{
	UINT32 a, angle_count(_it_msr->vectorCount2 - 1);
	UINT32 skip(0), ignored(_it_msr->vectorCount1 - _it_msr->vectorCount2);

	// move to first direction record which contains the derived angles
	_it_msr++;

	for (a=0; a<angle_count; ++a)		// for each angle
	{
		// cater for ignored directions
		if (_it_msr->ignore)
		{
			while (skip < ignored)
			{
				skip++;
				_it_msr++;
				if (!_it_msr->ignore)
					break;
			}
		}

		if (fabs(sigmaZeroSqRt_ - 0.0) < PRECISION_1E10)
			_it_msr->TStat = 0.0;
		else
			_it_msr->TStat = _it_msr->NStat / sigmaZeroSqRt_;
		_it_msr++;
	}

}

void dna_adjust::UpdateMsrTstatistic_GXY(it_vmsr_t& _it_msr)
{
	UINT32 cluster_msr, cluster_count(_it_msr->vectorCount1);
	UINT32 covariance_count;

	for (cluster_msr=0; cluster_msr<cluster_count; ++cluster_msr)
	{
		covariance_count = _it_msr->vectorCount2;

		// X measurement
		if (fabs(sigmaZeroSqRt_ - 0.0) < PRECISION_1E10)
			_it_msr->TStat = 0.0;
		else
			_it_msr->TStat = _it_msr->NStat / sigmaZeroSqRt_;
		_it_msr++;

		// Y measurement
		if (fabs(sigmaZeroSqRt_ - 0.0) < PRECISION_1E10)
			_it_msr->TStat = 0.0;
		else
			_it_msr->TStat = _it_msr->NStat / sigmaZeroSqRt_;
		_it_msr++;
		

		// Z measurement
		if (fabs(sigmaZeroSqRt_ - 0.0) < PRECISION_1E10)
			_it_msr->TStat = 0.0;
		else
			_it_msr->TStat = _it_msr->NStat / sigmaZeroSqRt_;

		// skip covariances until next point
		_it_msr += covariance_count * 3;

		if (covariance_count > 0)
			_it_msr++;
	}
}
	

void dna_adjust::UpdateMsrTstatistic(const UINT32& block)
{
	it_vUINT32 _it_block_msr;
	it_vmsr_t _it_msr;

	for (_it_block_msr=v_CML_.at(block).begin(); 
		_it_block_msr!=v_CML_.at(block).end(); 
		++_it_block_msr)
	{
		if (InitialiseandValidateMsrPointer(_it_block_msr, _it_msr))
			continue;

		switch (_it_msr->measType)
		{
		case 'D':	// Direction set
			// When a target direction is found, 
			// continue to next element.  
			if (_it_msr->vectorCount1 < 1)
				continue;
			UpdateMsrTstatistic_D(_it_msr);
			continue;
		case 'G':	// GPS Baseline  (treat as single-baseline cluster)
		case 'X':	// GPS Baseline cluster
		case 'Y':	// GPS Point cluster
			UpdateMsrTstatistic_GXY(_it_msr);
			continue;
		}

		if (fabs(sigmaZeroSqRt_ - 0.0) < PRECISION_1E10)
			_it_msr->TStat = 0.0;
		else
			_it_msr->TStat = _it_msr->NStat / sigmaZeroSqRt_;

	}
}
	
void dna_adjust::ComputeTstatistics()
{
	// Now that rigorous sigma zero has been computed, compute Student's T statistic
	UINT32 block;
	switch (projectSettings_.a.adjust_mode)
	{
	case PhasedMode:
		for (block=0; block<blockCount_; ++block)
		{
			// Update measurement t statistic
			UpdateMsrTstatistic(block);

		}
		break;
	case SimultaneousMode:
	case Phased_Block_1Mode:			// only block 1 is rigorous
		// Update measurement t statistic
		UpdateMsrTstatistic();
		break;
	}
}

void dna_adjust::ComputeStatistics()
{
	// initialise potential outlier count 
	potentialOutlierCount_ = 0;

	// Compute adjusted measurement precisions.  If this is a staged adjustment,
	// update normals and measured-computed matrices.
	ComputeAdjustedMsrPrecisions();

	// Compute whole-of-network statistics.
	ComputeChiSquareNetwork();
	
	// Compute DOF and sigma zero
	ComputeGlobalNetStat();

	// Update Student's T statistic
	if (projectSettings_.o._adj_msr_tstat)
		ComputeTstatistics();
	
	// Compute Pelzer's reliability for whole-of-network
	ComputeGlobalPelzer();

	switch (projectSettings_.a.adjust_mode)
	{
	case Phased_Block_1Mode:
		return;
	}
	
	// Perform global test
	ComputeGlobalTestStat();

}
	

void dna_adjust::ComputeandPrintAdjMsrOnIteration()
{
	// initialise potential outlier count 
	potentialOutlierCount_ = 0;
	bool printHeader(true);

	_it_uint32_u32u32_pair begin, end;
	begin = v_msr_block_.begin();

	// Compute and print adjusted measurements for all 
	// blocks at the end of an iteration
	for (UINT32 block(0); block<blockCount_; ++block)
	{
		if (projectSettings_.a.stage)
			DeserialiseBlockFromMappedFile(block, 4, sf_normals,
				sf_meas_minus_comp, sf_design, sf_prec_adj_msrs);

		// send subvector of measurements from this block
		end = begin + v_CML_.at(block).size();
		ComputeandPrintAdjMsrBlockOnIteration(block, v_uint32_u32u32_pair(begin, end), printHeader);
		begin = end;
		printHeader = false;

		if (projectSettings_.a.stage)
			UnloadBlock(block, 4, sf_normals,
				sf_meas_minus_comp, sf_design, sf_prec_adj_msrs);
	}
}

void dna_adjust::ComputeAdjMsrBlockOnIteration(const UINT32& block)
{
	// Compute adjusted measurement precisions (v_precAdjMsrsFull_)
	// from design and rigorous station variances 
	ComputePrecisionAdjMsrs(block);

	// Update measurement records and Pelzer's Global reliability
	UpdateMsrRecords(block);

	if (projectSettings_.a.stage)
		ComputeChiSquarePhased(block);		// This initialises chiSquared_ on each call

	if (projectSettings_.o._adj_msr_iteration ||
			projectSettings_.o._adj_stn_iteration ||
			projectSettings_.o._cmp_msr_iteration)
	{
		ComputeChiSquarePhased(block);		// This initialises chiSquared_ on each call

		switch (projectSettings_.a.adjust_mode)
		{
		case PhasedMode:
		case Phased_Block_1Mode:
			
			ComputeBlockTestStat(block);

			isAdjustmentQuestionable_ = (
				v_passFail_.at(block) != test_stat_pass ||
				// a much larger than expected outcome
				v_sigmaZero_.at(block) > 3.0 * v_chiSquaredUpperLimit_.at(block));
			break;
		default:
			ComputeGlobalNetStat();
			ComputeGlobalTestStat();

			isAdjustmentQuestionable_ = (
				// a much larger than expected outcome
				sigmaZero_ > 3.0 * chiSquaredUpperLimit_);
			break;
		}		
		
		
	}
}
	

void dna_adjust::ComputeandPrintAdjMsrBlockOnIteration(const UINT32& block, v_uint32_u32u32_pair msr_block, bool printHeader)
{
	// Compute adjusted measurements
	ComputeAdjMsrBlockOnIteration(block);
	
	// Print adjusted measurements
	printer_->PrintAdjMeasurements(msr_block, printHeader);
	
}
	

void dna_adjust::ComputeStatisticsOnIteration()
{
	// Compute whole-of-network statistics.
	ComputeChiSquareNetwork();
	
	// Compute DOF and sigma zero
	ComputeGlobalNetStat();

	// Perform global test
	ComputeGlobalTestStat();
}
	

void dna_adjust::ComputeChiSquare(const UINT32& block)
{
	UINT32 measurement_index(0);

	it_vUINT32 _it_block_msr;
	it_vmsr_t _it_msr;

	// Initialise chi-square
	chiSquared_ = 0.;	

	for (_it_block_msr=v_CML_.at(block).begin(); 
		_it_block_msr!=v_CML_.at(block).end(); 
		++_it_block_msr)
	{
		if (InitialiseandValidateMsrPointer(_it_block_msr, _it_msr))
			continue;

		switch (_it_msr->measType)
		{
		case 'A':	// Horizontal angle
		case 'B':	// Geodetic azimuth
		case 'C':	// Chord dist
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
			ComputeChiSquare_ABCEHIJKLMPQRSVZ(_it_msr, measurement_index, &v_measMinusComp_.at(block));
			break;
		
		case 'D':
			// When a target direction is found, continue to next element.  
			if (_it_msr->vectorCount1 < 1)
				continue;
			ComputeChiSquare_D(_it_msr, measurement_index, &v_measMinusComp_.at(block));
			break;
		
		case 'G':	// GPS Baseline
			ComputeChiSquare_G(_it_msr, measurement_index, &v_measMinusComp_.at(block));
			break;
		
		case 'X':
		case 'Y':
			ComputeChiSquare_XY(_it_msr, measurement_index, &v_measMinusComp_.at(block));
			break;
		}
	}
}
	

void dna_adjust::ComputeChiSquareNetwork()
{
	UINT32 block;
	double chiSquared(0.);
	
	// Compute measurement statistics
	switch (projectSettings_.a.adjust_mode)
	{
	case PhasedMode:
		measurementParams_ = 0;
		for (block=0; block<blockCount_; ++block)
		{
			// For staged adjustments, load block info
			if (projectSettings_.a.stage)
				DeserialiseBlockFromMappedFile(block, 1, sf_meas_minus_comp);

			ComputeChiSquarePhased(block);		// This initialises chiSquared_ on each call, 
			chiSquared += chiSquared_;			// so increment a local and update chiSquared_ later
			measurementParams_ += v_measurementParams_.at(block);

			// For staged adjustments, unload all matrix data
			if (projectSettings_.a.stage)
				UnloadBlock(block);
		}

		// update global
		chiSquared_ = chiSquared;
		break;
	case Phased_Block_1Mode:					// only block 1 is rigorous
	case SimultaneousMode:
		ComputeChiSquareSimultaneous();
		break;
	}
}
	

void dna_adjust::ComputeChiSquarePhased(const UINT32& block)
{
	// Compute adjusted measurement statistics
	ComputeChiSquare(block);
}
	

void dna_adjust::OutputLargestCorrection(std::string& formatted_msg)
{
	UINT32 stnIndex(0), x_coordElement(0);

	// load corrections matrix for block blockLargeCorr_
	if (projectSettings_.a.stage)
		DeserialiseBlockFromMappedFile(blockLargeCorr_, 1, 
			sf_corrections);

	// print maximum correction
	adj_file << std::setw(PRINT_VAR_PAD) << std::left << "Maximum station correction";

	// which station?
	stnIndex = (UINT32)floor(v_corrections_.at(blockLargeCorr_).maxvalueRow() / 3.);
	stnIndex = v_parameterStationList_.at(blockLargeCorr_).at(stnIndex);
	x_coordElement = (UINT32)(v_corrections_.at(blockLargeCorr_).maxvalueRow() % 3);

	switch (projectSettings_.a.adjust_mode)
	{
	case Phased_Block_1Mode:
	case PhasedMode:
		adj_file << "Block " << blockLargeCorr_+1 << 
			", station " << bstBinaryRecords_.at(stnIndex).stationName;
		break;
	case SimultaneousMode:
	default:
		adj_file << "Station " << bstBinaryRecords_.at(stnIndex).stationName;
		break;
	}

	adj_file << std::endl << std::setw(PRINT_VAR_PAD) << " ";

	// calculate error vector in local reference frame
	matrix_2d cart(3, 1), local(3, 1);
	switch (x_coordElement)
	{
	case 0:
		x_coordElement = v_corrections_.at(blockLargeCorr_).maxvalueRow();
		break;
	case 1:
		x_coordElement = v_corrections_.at(blockLargeCorr_).maxvalueRow() - 1;
		break;
	case 2:
		x_coordElement = v_corrections_.at(blockLargeCorr_).maxvalueRow() - 2;
		break;
	}

	cart.put(0, 0, v_corrections_.at(blockLargeCorr_).get(x_coordElement, 0));
	cart.put(1, 0, v_corrections_.at(blockLargeCorr_).get(x_coordElement+1, 0));
	cart.put(2, 0, v_corrections_.at(blockLargeCorr_).get(x_coordElement+2, 0));

	// OK, finished with the corrections, so unload if staged
	if (projectSettings_.a.stage)
		UnloadBlock(blockLargeCorr_, 1, sf_corrections);

	Rotate_CartLocal<double>(cart, &local,
		bstBinaryRecords_.at(stnIndex).currentLatitude,
		bstBinaryRecords_.at(stnIndex).currentLongitude);	

	double largestCorr = local.compute_maximum_value();
	
	if (fabs(largestCorr) > 0.000999)
		adj_file << std::fixed << std::setprecision(3) << 
			local.get(0, 0) << ", " << local.get(1, 0) << ", " << local.get(2, 0);
	else if (fabs(largestCorr) > 0.00009)
		adj_file << std::fixed << std::setprecision(4) << 
			local.get(0, 0) << ", " << local.get(1, 0) << ", " << local.get(2, 0);
	else
		adj_file << std::scientific << std::setprecision(1) << 
			local.get(0, 0) << ", " << local.get(1, 0) << ", " << local.get(2, 0);
	

	adj_file << " (e, n, up)" << std::endl << std::endl;

	std::stringstream ss;
	ss << std::fixed << std::setprecision(PRECISION_MTR_STN) << local.get(local.maxvalueRow(), 0) << " ";

	switch (local.maxvalueRow())
	{
	case 0:
		ss << std::setw(2) << std::left << "e";
		break;
	case 1:
		ss << std::setw(2) << std::left << "n";
		break;
	case 2:
		ss << std::setw(2) << std::left << "up";
		break;
	}

	formatted_msg = ss.str();
}
	
void dna_adjust::ComputePrecisionAdjMsrs(const UINT32& block /*= 0*/)
{
	if (projectSettings_.a.report_mode)
		return;

	// A*V-1*At, where:
	//   - A is design matrix
	//   - V is the inverse of the normals (i.e. precision of estimates)
	v_precAdjMsrsFull_.at(block).zero();

	UINT32 design_row(0);
	UINT32 precadjmsr_row(0);
	
	it_vUINT32 _it_block_msr;
	it_vmsr_t _it_msr;

	matrix_2d *design(&v_design_.at(block)), *aposterioriVariances(&v_normals_.at(block));

	// Measurements can only ever appear once in the whole CML.  That is, no one measurement will be found
	// in two or more blocks.  Therefore, unlike precisions of adjusted stations (which may appear in one
	// or more blocks), precisions of adjusted measurements are unique.
	for (_it_block_msr=v_CML_.at(block).begin(); _it_block_msr!=v_CML_.at(block).end(); ++_it_block_msr)
	{
		if (InitialiseandValidateMsrPointer(_it_block_msr, _it_msr))
			continue;

		// Build  At * V-1 (diagonals only as full covariances are not required)
		switch (_it_msr->measType)
		{
		case 'A':	// Horizontal angle
			ComputePrecisionAdjMsrs_A(block, 
				GetBlkMatrixElemStn1(block, &_it_msr), 
				GetBlkMatrixElemStn2(block, &_it_msr), 
				GetBlkMatrixElemStn3(block, &_it_msr),
				design, aposterioriVariances,
				design_row, precadjmsr_row);
			break;
		case 'D':	// Direction set
			// When a target direction is found, continue to next element.  
			if (_it_msr->vectorCount1 < 1)
				continue;
			ComputePrecisionAdjMsrs_D(block, _it_msr, 
				design, aposterioriVariances,
				design_row, precadjmsr_row);
			break;
		// Single station measurements
		case 'H':	// Orthometric height
		case 'I':	// Astronomic latitude
		case 'J':	// Astronomic longitude
		case 'P':	// Geodetic latitude
		case 'Q':	// Geodetic longitude
		case 'R':	// Ellipsoidal height
			ComputePrecisionAdjMsrs_HIJPQR(block,
				GetBlkMatrixElemStn1(block, &_it_msr),				
				design, aposterioriVariances,
				design_row, precadjmsr_row);
			break;
		// Two station measurements
		case 'B':	// Geodetic azimuth
		case 'C':	// Chord dist
		case 'E':	// Ellipsoid arc
		case 'K':	// Astronomic azimuth
		case 'L':	// Level difference
		case 'M':	// MSL arc
		case 'S':	// Slope distance
		case 'V':	// Zenith distance
		case 'Z':	// Vertical angle
			ComputePrecisionAdjMsrs_BCEKLMSVZ(block, 
				GetBlkMatrixElemStn1(block, &_it_msr), 
				GetBlkMatrixElemStn2(block, &_it_msr), 
				design, aposterioriVariances,
				design_row, precadjmsr_row);
			break;
		case 'G':	// GPS Baseline
		case 'X':	// GPS Baseline cluster
			ComputePrecisionAdjMsrs_GX(block, _it_msr, 
				aposterioriVariances, design_row, precadjmsr_row);
			break;
		case 'Y':	// GPS Point cluster
			ComputePrecisionAdjMsrs_Y(block, _it_msr, 
				aposterioriVariances, design_row, precadjmsr_row);
			break;		
		default:
			std::stringstream ss;
			ss << "ComputePrecisionAdjMsrs(): Unknown measurement type - '" << static_cast<std::string>(&(_it_msr->measType)) <<
				"'." << std::endl;
			SignalExceptionAdjustment(ss.str(), block);
		}
	}
}
	

void dna_adjust::ComputePrecisionAdjMsrs_A(const UINT32& block, const UINT32& stn1, const UINT32& stn2, const UINT32& stn3, 
											  matrix_2d* design, matrix_2d* aposterioriVariances, 
											  UINT32& design_row, UINT32& precadjmsr_row)
{
	// Horizontal angle
	double part_1[9] = {0.,0.,0.,0.,0.,0.,0.,0.,0.};
	UINT32 stations[3] = {stn1, stn2, stn3};
	UINT32 station_count(3);
	UINT32 station_i[3] = {0, 3, 6};
	UINT32 var, elem(0);
	UINT32 j, i, s;

	for (s=0; s<station_count; ++s)	// for every station
	{
		for (i=0; i<3; ++i)		// X, Y, Z
		{
			for (j=0; j<station_count; ++j)			// for every correlated station
			{	
				var = stations[j];
				part_1[elem] += design->get(design_row, var) * aposterioriVariances->get(var, stations[s]+i);
				part_1[elem] += design->get(design_row, var+1) * aposterioriVariances->get(var+1, stations[s]+i);
				part_1[elem] += design->get(design_row, var+2) * aposterioriVariances->get(var+2, stations[s]+i);
			}
			elem++;
		}		
	}

	for (s=0; s<station_count; ++s)		// for every station
		for (i=0; i<3; ++i)				// X, Y, Z
			v_precAdjMsrsFull_.at(block).elementadd(precadjmsr_row, 0, part_1[station_i[s]+i] * design->get(design_row, stations[s]+i));
	
	design_row++;
	precadjmsr_row++;
}
	

void dna_adjust::ComputePrecisionAdjMsrs_D(const UINT32& block, it_vmsr_t& _it_msr, 
											  matrix_2d* design, matrix_2d* aposterioriVariances, 
											  UINT32& design_row, UINT32& precadjmsr_row)
{
	UINT32 stn1, stn2, stn3;
	UINT32 a, angle_count(_it_msr->vectorCount2 - 1);		// number of directions excluding the RO
	UINT32 skip(0), ignored(_it_msr->vectorCount1 - _it_msr->vectorCount2);

	for (a=0; a<angle_count; ++a)
	{
		// On the first time this loop is entered, the stn1 and stn2 will be instrument and RO
		// Then, all following directions will be instrument and target.
		stn1 = GetBlkMatrixElemStn1(block, &_it_msr);
		stn2 = GetBlkMatrixElemStn2(block, &_it_msr);
		_it_msr++;

		// cater for ignored directions
		if (_it_msr->ignore)
		{
			while (skip < ignored)
			{
				skip++;
				_it_msr++;
				if (!_it_msr->ignore)
					break;
			}
		}

		stn3 = GetBlkMatrixElemStn2(block, &_it_msr);

		ComputePrecisionAdjMsrs_A(block, stn1, stn2, stn3, 
			design, aposterioriVariances, design_row, precadjmsr_row);
	}
}
	

void dna_adjust::ComputePrecisionAdjMsrs_BCEKLMSVZ(const UINT32& block, const UINT32& stn1, const UINT32& stn2, 
											  matrix_2d* design, matrix_2d* aposterioriVariances, 
											  UINT32& design_row, UINT32& precadjmsr_row)
{
	// Two station measurement
	double part_1[6] = {0.,0.,0.,0.,0.,0.};
	UINT32 stations[2] = {stn1, stn2};
	UINT32 station_count(2);
	UINT32 station_i[2] = {0, 3};
	UINT32 var, elem(0);
	UINT32 j, i, s;

	for (s=0; s<station_count; ++s)		// for every station
	{
		for (i=0; i<3; ++i)				// X, Y, Z
		{
			for (j=0; j<station_count; ++j)			// for every correlated station
			{	
				var = stations[j];
				part_1[elem] += design->get(design_row, var) * aposterioriVariances->get(var, stations[s]+i);
				part_1[elem] += design->get(design_row, var+1) * aposterioriVariances->get(var+1, stations[s]+i);
				part_1[elem] += design->get(design_row, var+2) * aposterioriVariances->get(var+2, stations[s]+i);
			}
			elem++;
		}		
	}

	for (s=0; s<station_count; ++s)	// for every station
		for (i=0; i<3; ++i)		// X, Y, Z
			v_precAdjMsrsFull_.at(block).elementadd(precadjmsr_row, 0, part_1[station_i[s]+i] * design->get(design_row, stations[s]+i));
	
	design_row++;
	precadjmsr_row++;
}
	

void dna_adjust::ComputePrecisionAdjMsrs_HIJPQR(const UINT32& block, const UINT32& stn1, 
											  matrix_2d* design, matrix_2d* aposterioriVariances, 
											  UINT32& design_row, UINT32& precadjmsr_row)
{
	// Single station measurement
	double part_1[3] = {0.,0.,0.};
	UINT32 i, elem;

	for (elem=0,i=0; i<3; ++i)		// X, Y, Z
	{	
		part_1[elem] += design->get(design_row, stn1) * aposterioriVariances->get(stn1, stn1+i);
		part_1[elem] += design->get(design_row, stn1+1) * aposterioriVariances->get(stn1+1, stn1+i);
		part_1[elem] += design->get(design_row, stn1+2) * aposterioriVariances->get(stn1+2, stn1+i);
		elem++;
	}

	for (i=0; i<3; ++i)				// X, Y, Z
		v_precAdjMsrsFull_.at(block).elementadd(precadjmsr_row, 0, part_1[i] * design->get(design_row, stn1+i));
	
	design_row++;
	precadjmsr_row++;
}
		

void dna_adjust::ComputePrecisionAdjMsrs_GX(const UINT32& block, it_vmsr_t& _it_msr, 
											  matrix_2d* aposterioriVariances, 
											  UINT32& design_row, UINT32& precadjmsr_row)
{
	UINT32 cluster_bsl, baseline_count(_it_msr->vectorCount1);
	UINT32 stn1, stn2;
	UINT32 i, j;
	matrix_2d precision_bsl(3, 3);
	precision_bsl.zero();

	for (cluster_bsl=0; cluster_bsl<baseline_count; ++cluster_bsl)
	{
		stn1 = GetBlkMatrixElemStn1(block, &_it_msr);
		stn2 = GetBlkMatrixElemStn2(block, &_it_msr);

		Precision_Adjusted_GNSS_bsl<double>(*aposterioriVariances,
			stn1, stn2, &precision_bsl, false);

		for (i=0; i<3; ++i)
			for (j=i; j<3; ++j, ++precadjmsr_row)
				v_precAdjMsrsFull_.at(block).put(precadjmsr_row, 0, precision_bsl.get(i, j));
		
		design_row += 3;
		_it_msr += 3;
	}
}
	

void dna_adjust::ComputePrecisionAdjMsrs_Y(const UINT32& block, it_vmsr_t& _it_msr, 
											  matrix_2d* aposterioriVariances, 
											  UINT32& design_row, UINT32& precadjmsr_row)
{
	UINT32 cluster_pnt, point_count(_it_msr->vectorCount1);
	UINT32 stn1, i, j;

	for (cluster_pnt=0; cluster_pnt<point_count; ++cluster_pnt)
	{
		stn1 = GetBlkMatrixElemStn1(block, &_it_msr);

		for (i=0; i<3; ++i)
		{
			for (j=i; j<3; ++j, ++precadjmsr_row)
			{
				v_precAdjMsrsFull_.at(block).put(precadjmsr_row, 0, 
					aposterioriVariances->get(stn1+i, stn1+j));
			}
		}
		
		design_row +=3;
		_it_msr += 3;
	}
}
	

bool dna_adjust::InitialiseandValidateMsrPointer(const it_vUINT32& _it_block_msr, it_vmsr_t& _it_msr)
{
	_it_msr = bmsBinaryRecords_.begin() + (*_it_block_msr);

	// Is this measurement ignored?
	if (_it_msr->ignore)
		return true;

	// Only pointers to the starting element of each measurement are kept. 
	// Hence, Y and Z and variance-covariance components are regarded as
	// non-measurement elements.  When one of these is found, 
	// continue to next element.  See declaration of measRecord.measStart 
	// for more info.
	if (_it_msr->measStart != xMeas)
		return true;

	return false;
}

// store adjusted measurements and corrections
void dna_adjust::UpdateMsrRecords(const UINT32& block)
{	
	UINT32 msr_row(0), precadjmsr_row(0);

	it_vUINT32 _it_block_msr;
	it_vmsr_t _it_msr;

	for (_it_block_msr=v_CML_.at(block).begin(); 
		_it_block_msr!=v_CML_.at(block).end();	
		++_it_block_msr)
	{
		if (InitialiseandValidateMsrPointer(_it_block_msr, _it_msr))
			continue;

		switch (_it_msr->measType)
		{
		case 'D':	// Direction set
			// When a target direction is found, 
			// continue to next element.  
			if (_it_msr->vectorCount1 < 1)
				continue;
			UpdateMsrRecords_D(block, _it_msr, msr_row, precadjmsr_row);
			continue;
		case 'G':	// GPS Baseline  (treat as single-baseline cluster)
		case 'X':	// GPS Baseline cluster
		case 'Y':	// GPS Point cluster
			UpdateMsrRecords_GXY(block, _it_msr, msr_row, precadjmsr_row);
			continue;
		}

		// Update measurement record with adjusted measurement values
		UpdateMsrRecord(block, _it_msr, msr_row++, precadjmsr_row++, _it_msr->term2);

	}
}
	

// store adjusted measurements and corrections
void dna_adjust::UpdateMsrRecords_D(const UINT32& block, it_vmsr_t& _it_msr, UINT32& msr_row, UINT32& precadjmsr_row)
{
	UINT32 a, angle_count(_it_msr->vectorCount2 - 1);
	UINT32 skip(0), ignored(_it_msr->vectorCount1 - _it_msr->vectorCount2);

	// move to first direction record which contains the derived angles
	_it_msr++;
	
	for (a=0; a<angle_count; ++a)		// for each angle
	{
		// cater for ignored directions
		if (_it_msr->ignore)
		{
			while (skip < ignored)
			{
				skip++;
				_it_msr++;
				if (!_it_msr->ignore)
					break;
			}
		}

		UpdateMsrRecord(block, _it_msr, msr_row, precadjmsr_row, _it_msr->scale2);
		_it_msr++;
		msr_row++;
		precadjmsr_row++;
	}
}
	

// store adjusted measurements and corrections
void dna_adjust::UpdateMsrRecords_GXY(const UINT32& block, it_vmsr_t& _it_msr, UINT32& msr_row, UINT32& precadjmsr_row)
{
	UINT32 cluster_msr, cluster_count(_it_msr->vectorCount1);
	UINT32 covariance_count;
	
	for (cluster_msr=0; cluster_msr<cluster_count; ++cluster_msr)
	{
		covariance_count = _it_msr->vectorCount2;

		// X measurement
		UpdateMsrRecord(block, _it_msr, msr_row, precadjmsr_row, _it_msr->term2);
		_it_msr++;
		msr_row++;
		precadjmsr_row += 3;

		// Y measurement
		UpdateMsrRecord(block, _it_msr, msr_row, precadjmsr_row, _it_msr->term3);
		_it_msr++;
		msr_row++;
		precadjmsr_row += 2;

		// Z measurement
		UpdateMsrRecord(block, _it_msr, msr_row, precadjmsr_row, _it_msr->term4);
		msr_row++;
		precadjmsr_row++;

		// skip covariances until next point
		_it_msr += covariance_count * 3;
		
		if (covariance_count > 0)
			_it_msr++;
	}
}
	

void dna_adjust::UpdateMsrRecord(const UINT32& block, it_vmsr_t& _it_msr, 
	const UINT32& msr_row, const UINT32& precadjmsr_row, const double& measPrec)
{
	// set adjusted measurement correction
	_it_msr->measCorr = -v_measMinusComp_.at(block).get(msr_row, 0);
	
	// apply adjusted measurement correction
	if (_it_msr->measType == 'D')
	{
		_it_msr->measAdj = _it_msr->scale1 + _it_msr->measCorr;
		if (_it_msr->measAdj > TWO_PI)
			_it_msr->measAdj -= TWO_PI;
	}
	else
		_it_msr->measAdj = _it_msr->term1 + _it_msr->measCorr;

	UINT32 stn1, stn2;
	it_vstn_t_const stn1_it, stn2_it;
	matrix_2d* estimatedStations(&v_estimatedStations_.at(block));

	// Recompute measurements using the original types
	switch (_it_msr->measType)
	{
	case 'E':
		stn1 = GetBlkMatrixElemStn1(block, &_it_msr); 
		stn2 = GetBlkMatrixElemStn2(block, &_it_msr);
		stn1_it = bstBinaryRecords_.begin() + _it_msr->station1;
		stn2_it = bstBinaryRecords_.begin() + _it_msr->station2;
	
		// reduce adjusted measurement to ellipsoid arc
		_it_msr->measAdj = EllipsoidChordtoEllipsoidArc<double>(
			_it_msr->measAdj,
			estimatedStations->get(stn1, 0),
			estimatedStations->get(stn1+1, 0),
			estimatedStations->get(stn1+2, 0),
			estimatedStations->get(stn2, 0),
			estimatedStations->get(stn2+1, 0),
			estimatedStations->get(stn2+2, 0),
			stn1_it->currentLatitude,
			stn1_it->currentLongitude,
			stn2_it->currentLatitude,
			datum_.GetEllipsoidRef());
		break;
	
	case 'M':
		stn1 = GetBlkMatrixElemStn1(block, &_it_msr); 
		stn2 = GetBlkMatrixElemStn2(block, &_it_msr);
		stn1_it = bstBinaryRecords_.begin() + _it_msr->station1;
		stn2_it = bstBinaryRecords_.begin() + _it_msr->station2;
	
		// reduce adjusted measurement to MSL arc
		_it_msr->measAdj = EllipsoidChordtoMSLArc<double>(
			_it_msr->measAdj,
			stn1_it->currentLatitude,
			stn2_it->currentLatitude,
			stn1_it->geoidSep,
			stn2_it->geoidSep,
			datum_.GetEllipsoidRef());
		break;

	// Apply N value
	case 'H':
	case 'L':
		_it_msr->measAdj -= _it_msr->preAdjCorr;
		break;

	// Apply correction for deflections of the vertical
	case 'A':
	case 'D':
	case 'I':
	case 'J':
	case 'K':
	case 'Z':
		// add deflection of the vertical to adjusted value
		_it_msr->measAdj += _it_msr->preAdjCorr;
		break;
	case 'V':
		// apply deflection of the vertical to the adjusted
		// measurement to get the commensurate measurement value
		_it_msr->measAdj -= _it_msr->preAdjCorr;
		break;
	}

	_it_msr->measAdjPrec = v_precAdjMsrsFull_.at(block).get(precadjmsr_row, 0);		// precision of adjusted measurement
	_it_msr->residualPrec = measPrec - _it_msr->measAdjPrec;						// residual precision

	// Often, residualPrec contains a very small negative 
	// value, e.g. -1.5 e-35, effectively zero!
	// To circumvent errors arising from taking the square root 
	// of a negative number, first compute the absolute value.
	// While this is theoretically incorrect, the impact is of no practical
	// consequence. P. Collier
	if (_it_msr->residualPrec < 0.0)
		_it_msr->residualPrec = fabs(_it_msr->residualPrec);

	// Update Pelzer reliability, N-Stat, T-Stat
	UpdateMsrRecordStats(_it_msr, measPrec);

	if (fabs(_it_msr->NStat) > criticalValue_)
		potentialOutlierCount_++;
}

// Compute and update Pelzer's reliability value and N statistic
// Student's T statistic is computed once global sigma-zero has been computed
void dna_adjust::UpdateMsrRecordStats(it_vmsr_t& _it_msr, const double& measPrec)
{
	_it_msr->PelzerRel = sqrt(measPrec) / sqrt(_it_msr->residualPrec);
	if (_it_msr->PelzerRel < 0. || _it_msr->PelzerRel > STABLE_LIMIT)
		_it_msr->PelzerRel = UNRELIABLE;

	_it_msr->NStat = _it_msr->measCorr / sqrt(_it_msr->residualPrec);
}
	

// Compute Pelzer's global reliability
void dna_adjust::ComputeGlobalPelzer()
{
	UINT32 block, numMsr(0);
	double sum(0.);

	it_vUINT32 _it_block_msr;
	it_vmsr_t _it_msr;

	for (block=0; block<blockCount_; ++block)
	{
		for (_it_block_msr=v_CML_.at(block).begin(); 
			_it_block_msr!=v_CML_.at(block).end();	
			++_it_block_msr)
		{
			if (InitialiseandValidateMsrPointer(_it_block_msr, _it_msr))
				continue;

			switch (_it_msr->measType)
			{
			case 'D':	// Direction set
				// When a target direction is found, continue to next element.  
				if (_it_msr->vectorCount1 < 1)
					continue;
				ComputeGlobalPelzer_D(_it_msr, numMsr, sum);
				continue;

			case 'G':	// GPS Baseline (treat as single-baseline cluster)
			case 'X':	// GPS Baseline cluster
			case 'Y':	// GPS Point cluster
				ComputeGlobalPelzer_GXY(_it_msr, numMsr, sum);
				continue;
			}
			
			// All measurement types
			if (_it_msr->PelzerRel > 0. && _it_msr->PelzerRel < STABLE_LIMIT)
			{
				sum += (_it_msr->PelzerRel * _it_msr->PelzerRel - 1.);
				numMsr++;
			}
			else
				_it_msr->PelzerRel = UNRELIABLE;
		}

		switch (projectSettings_.a.adjust_mode)
		{
		case Phased_Block_1Mode:
			break;
		}
	}

	// Ok, now compute Pelzer's global reliability
	if (numMsr > 0)
		globalPelzerReliability_ = sqrt(sum / numMsr);
	else
		globalPelzerReliability_  = UNRELIABLE;

}

// Compute Pelzer's global reliability
void dna_adjust::ComputeGlobalPelzer_D(it_vmsr_t& _it_msr, UINT32& numMsr, double& sum)
{
	UINT32 a, angle_count(_it_msr->vectorCount2 - 1);
	UINT32 skip(0), ignored(_it_msr->vectorCount1 - _it_msr->vectorCount2);

	// move to first direction record which contains the derived angles
	_it_msr++;
	
	for (a=0; a<angle_count; ++a)		// for each angle
	{
		// cater for ignored directions
		if (_it_msr->ignore)
		{
			while (skip < ignored)
			{
				skip++;
				_it_msr++;
				if (!_it_msr->ignore)
					break;
			}
		}

		if (_it_msr->PelzerRel > 0. && _it_msr->PelzerRel < UNRELIABLE)
		{
			sum += (_it_msr->PelzerRel * _it_msr->PelzerRel - 1.);
			numMsr++;
		}
		else
			_it_msr->PelzerRel = UNRELIABLE;
		_it_msr++;
	}
}
		

// Compute Pelzer's global reliability
void dna_adjust::ComputeGlobalPelzer_GXY(it_vmsr_t& _it_msr, UINT32& numMsr, double& sum)
{
	UINT32 cluster_msr, cluster_count(_it_msr->vectorCount1);
	UINT32 covariance_count, i;

	for (cluster_msr=0; cluster_msr<cluster_count; ++cluster_msr)
	{
		covariance_count = _it_msr->vectorCount2;

		// X, Y, Z measurement
		for (i=0; i<3; ++i)
		{
			if (_it_msr->PelzerRel > 0. && _it_msr->PelzerRel < UNRELIABLE)
			{
				sum += (_it_msr->PelzerRel * _it_msr->PelzerRel - 1.);
				numMsr++;
			}
			else
				_it_msr->PelzerRel = UNRELIABLE;

			if (i < 2)
				_it_msr++;
		}

		// skip covariances until next point
		_it_msr += covariance_count * 3;
		
		if (covariance_count > 0)
			_it_msr++;
	}

}
		

void dna_adjust::ComputeChiSquare_ABCEHIJKLMPQRSVZ(const it_vmsr_t& _it_msr, UINT32& measurement_index, matrix_2d* measMinusComp)
{
	chiSquared_ +=  
		measMinusComp->get(measurement_index, 0) * 
		measMinusComp->get(measurement_index, 0) / _it_msr->term2;

	measurement_index++;
}
	

void dna_adjust::ComputeChiSquare_D(it_vmsr_t& _it_msr, UINT32& measurement_index, matrix_2d* measMinusComp)
{
	UINT32 a, angle_count(_it_msr->vectorCount2 - 1);
	UINT32 skip(0), ignored(_it_msr->vectorCount1 - _it_msr->vectorCount2);

	// move to first direction record which contains the derived angles
	_it_msr++;
	
	for (a=0; a<angle_count; ++a)		// for each angle
	{
		// cater for ignored directions
		if (_it_msr->ignore)
		{
			while (skip < ignored)
			{
				skip++;
				_it_msr++;
				if (!_it_msr->ignore)
					break;
			}
		}

		chiSquared_ +=
			measMinusComp->get(measurement_index, 0) * 
			measMinusComp->get(measurement_index, 0) / _it_msr->scale2;		//variance (angle)

		measurement_index++;
		_it_msr++;
	}
}
	

void dna_adjust::FormInverseVarianceMatrix(matrix_2d* vmat, bool LOWER_IS_CLEARED)
{
	if (vmat->rows() == 1)
	{
		vmat->put(0, 0, 1./vmat->get(0, 0));
		return;
	}
	
	// As of version 3.2.0, force all inversions to use MKL.  This change
	// is enforced for two reasons:
	// 1. Sweep, Gaussian inverse and numerical recipes cholesky
	//    all require matrix data to be stored in row wise fashion, upper 
	//    triangle only, whereas contiguous matrix class stores matrix data 
	//    in column wise fashion, lower triangle.
	// 2. Sweep, Gaussian never really offered a stable solution.
	// The following switch is kept in case future development warrants
	// re-offering the option to choose (in which case, those matrix inversion
	// functions must be revised to cater for column wise, lower triangular
	// order.

	// TODO: All functions which load variance matrix from binary files
	// store the data in upper triangular form.  This could be changed to
	// lower triangular form, thus alleviating the need to pass 
	// LOWER_IS_CLEARED.  That is, force all operations to use a lower 
	// triangular matrix.  Not sure if this would create an efficiency or not.
	switch (projectSettings_.a.inverse_method_msr)
	{
//	case Gaussian:
//		vmat->gaussianinverse();
//		break;
//	case Sweep:
//		vmat->sweepinverse();
//		break;
	case Cholesky_mkl:
	default:
		// Inversion using Intel MKL
		vmat->cholesky_inverse(LOWER_IS_CLEARED);
		break;
	// choleskyinverse broke once the storage order of the matrix buffer was
	// changed from row-wise to column wise.
//	case Cholesky:
//	default:
		// Inversion using Numerical Recipes
//		vmat->choleskyinverse(LOWER_IS_CLEARED);
//		break;
	}
}
	

void dna_adjust::FormInverseGPSVarianceMatrix(const it_vmsr_t& _it_msr, matrix_2d* vmat)
{
	// 1. Get upper triangular a-priori measurements variance matrix
	GetGPSVarianceMatrix<it_vmsr_t>(_it_msr, vmat);

	// 2. Inverse
	FormInverseVarianceMatrix(vmat, true);
}
	

void dna_adjust::ComputeChiSquare_G(const it_vmsr_t& _it_msr, UINT32& measurement_index, matrix_2d* measMinusComp)
{
	matrix_2d V(3, 3);		// a-priori measurements variance matrix
	
	// 1. Form inverse variance matrix for GPS measurement
	FormInverseGPSVarianceMatrix(_it_msr, &V);

	// 2. Form residuals matrix
	UINT32 row, col;
	double cs(0.);
	for (row=0; row<3; ++row)
		for (col=0; col<3; ++col)
			cs += V.get(row, col) * 
				measMinusComp->get(measurement_index + row, 0) * 
				measMinusComp->get(measurement_index + col, 0);

	chiSquared_ += cs;
	measurement_index += 3;
}
	

void dna_adjust::ComputeChiSquare_XY(const it_vmsr_t& _it_msr, UINT32& measurement_index, matrix_2d* measMinusComp)
{
	// compute At * Vm-1 * A
	UINT32 element_count(_it_msr->vectorCount1);
	UINT32 variance_dim(element_count * 3);
	
	matrix_2d V(variance_dim, variance_dim);		// a-priori measurements variance matrix
		
	// 1. Form inverse variance matrix for GPS measurement
	FormInverseGPSVarianceMatrix(_it_msr, &V);

	// 2. Form residuals matrix
	matrix_2d r(measMinusComp->submatrix(measurement_index, 0, variance_dim, 1));
	
	matrix_2d rt(r.transpose());
	matrix_2d rt_Vinv(1, variance_dim);
	matrix_2d rt_Vinv_r(1, 1);
	
	//rt_Vinv.multiply(rt, V);
	rt_Vinv.multiply(r, "T", V, "N");

	//rt_Vinv_r.multiply(rt_Vinv, r);
	rt_Vinv_r.multiply(rt_Vinv, "N", r, "N");

	chiSquared_ += rt_Vinv_r.get(0, 0);
	measurement_index += variance_dim;
}
	

void dna_adjust::OpenOutputFileStreams()
{
	try {
		// Create adj file.  Throws runtime_error on failure.
		file_opener(adj_file, projectSettings_.o._adj_file);
	}
	catch (const std::runtime_error& e) {
		SignalExceptionAdjustment(e.what(), 0);
	}

	try {
		// Create xyz file.  Throws runtime_error on failure.
		file_opener(xyz_file, projectSettings_.o._xyz_file);
	}
	catch (const std::runtime_error& e) {
		SignalExceptionAdjustment(e.what(), 0);
	}

	if (projectSettings_.g.verbose > 0)
	{
		std::string debug_output_file = 
			projectSettings_.g.output_folder + FOLDER_SLASH + projectSettings_.g.network_name;

		switch (projectSettings_.a.adjust_mode)
		{
		case Phased_Block_1Mode:
		case PhasedMode:

			debug_output_file += ".phased";

			if (projectSettings_.a.adjust_mode == Phased_Block_1Mode)
				debug_output_file += "-block1";
			else if (projectSettings_.a.stage)
				debug_output_file += "-stage";
			else if (projectSettings_.a.multi_thread)
				debug_output_file += "-mt";
			break;
		case SimultaneousMode:
			debug_output_file += ".simult";
			break;
		}
		
		debug_output_file += ".dbg";
		
		try {
			// Create debug file.  Throws runtime_error on failure.
			file_opener(debug_file, debug_output_file);
		}
		catch (const std::runtime_error& e) {
			SignalExceptionAdjustment(e.what(), 0);
		}
	}

	if (projectSettings_.a.stage)
	{
		OpenStageFileStreams();
		ReserveBlockMapRegions();
	}
}
	
void dna_adjust::SortStationsbyFileOrder(vUINT32& v_blockStations)
{
	CompareStnFileOrder<station_t, UINT32> stnorderCompareFunc(&bstBinaryRecords_);
	std::sort(v_blockStations.begin(), v_blockStations.end(), stnorderCompareFunc);
}

void dna_adjust::SortStationsbyID(vUINT32& v_blockStations)
{
	std::sort(v_blockStations.begin(), v_blockStations.end());
}
	
void dna_adjust::PrintCorStations(std::ostream &cor_file, const UINT32& block)
{
	vUINT32 v_blockStations(v_parameterStationList_.at(block));
	
	// if required, sort stations according to original station file order
	if (projectSettings_.o._sort_stn_file_order)
		SortStationsbyFileOrder(v_blockStations);
	
	// Use the already initialized printer
	printer_->PrintStationCorrelations(cor_file, block);
	
	// Continue with existing detailed implementation
	switch (projectSettings_.a.adjust_mode)
	{
	case PhasedMode:
	case Phased_Block_1Mode:		// only the first block is rigorous
		cor_file << "Block " << block + 1 << std::endl;
		break;
	}

	cor_file << std::setw(STATION) << std::left << "Station" <<
		std::setw(PAD2) << " " << 
		std::right << std::setw(MSR) << "Azimuth" << 
		std::right << std::setw(MSR) << "V. Angle" << 
		std::right << std::setw(MSR) << "S. Distance" << 
		std::right << std::setw(MSR) << "H. Distance" << 
		std::right << std::setw(HEIGHT) << "east" << 
		std::right << std::setw(HEIGHT) << "north" << 
		std::right << std::setw(HEIGHT) << "up" << std::endl;

	UINT32 i, j = STATION+PAD2+MSR+MSR+MSR+MSR+HEIGHT+HEIGHT+HEIGHT;

	for (i=0; i<j; ++i)
		cor_file << "-";
	cor_file << std::endl;

	UINT32 mat_idx, stn;

	pv_mat_2d estimates(&v_rigorousStations_);
	if (projectSettings_.a.adjust_mode == SimultaneousMode)
		estimates = &v_estimatedStations_;

	// Print stations according to the user-defined sort order
	for (i=0; i<v_blockStationsMap_.at(block).size(); ++i)
	{
		stn = v_blockStations.at(i);
		mat_idx = v_blockStationsMap_.at(block)[stn] * 3;

		printer_->PrintCorStation(cor_file, block, stn, mat_idx, 
			&estimates->at(block));
	}

	cor_file << std::endl;

	// return sort order to alpha-numeric
	if (projectSettings_.o._sort_stn_file_order)
		SortStationsbyID(v_blockStations);

}

void dna_adjust::UpdateGeographicCoordsPhased(const UINT32& block, matrix_2d* estimatedStations)
{
	UINT32 i(0), j(0);
	vUINT32 v_blockStations(v_parameterStationList_.at(block));
	it_vstn_appear _it_appear(v_paramStnAppearance_.at(block).begin());

	for (j=0, i=0; i<estimatedStations->rows(); j++, i+=3, ++_it_appear)
	{
		// The same station may appear in several blocks.  So, only
		// update (once) when this is the first time this station 
		// appears
		if (!_it_appear->first_appearance_fwd)
			continue;

		CartToGeo<double>(estimatedStations->get(i, 0), estimatedStations->get(i+1, 0), estimatedStations->get(i+2, 0),
			&(bstBinaryRecords_.at(v_blockStations.at(j)).currentLatitude), 
			&(bstBinaryRecords_.at(v_blockStations.at(j)).currentLongitude), 
			&(bstBinaryRecords_.at(v_blockStations.at(j)).currentHeight), 
			datum_.GetEllipsoidRef());
	}
}
	

void dna_adjust::UpdateGeographicCoords()
{
	it_vUINT32 _it_stn;
	UINT32 stn(0);
	
	// all stations in simultaneous mode are kept in ISL
	for (_it_stn=v_ISL_.at(0).begin(); _it_stn!=v_ISL_.at(0).end(); ++_it_stn, stn+=3)
	{
		CartToGeo<double>(v_estimatedStations_.at(0).get(stn, 0), v_estimatedStations_.at(0).get(stn+1, 0), v_estimatedStations_.at(0).get(stn+2, 0),
			&(bstBinaryRecords_.at(*_it_stn).currentLatitude), 
			&(bstBinaryRecords_.at(*_it_stn).currentLongitude), 
			&(bstBinaryRecords_.at(*_it_stn).currentHeight), 
			datum_.GetEllipsoidRef());
	}
}

void dna_adjust::UpdateIgnoredMeasurements_A(pit_vmsr_t _it_msr, bool storeOriginalMeasurement)
{
	// initialise measurement (on the first adjustment only!)
	if (storeOriginalMeasurement)
		(*_it_msr)->preAdjMeas = (*_it_msr)->term1;

	// Use v_blockStationsMapUnique_ to get the correct block for each 
	// station in the measurement
	_it_pair_u32u32_uint32 it_stnmap_range;
	_it_u32u32_uint32_pair _it_bsmu;

	// Get estimated station coordinates matrix and index of the 
	// station within the matrix for station 1
	it_stnmap_range = equal_range(v_blockStationsMapUnique_.begin(), v_blockStationsMapUnique_.end(),
		(*_it_msr)->station1, CompareBlockStationMapUnique_Station<UINT32, u32u32_uint32_pair>());
	if (it_stnmap_range.first == it_stnmap_range.second)
	{
		//TRACE("%d not found\n", (*_it_msr)->station1);
		std::stringstream ss;
		ss << "UpdateIgnoredMeasurements(): Station " << (*_it_msr)->station1 <<
			" was not found in the station map." << std::endl;
		SignalExceptionAdjustment(ss.str(), 0);
	}

	// Get the index of the station in the map
	_it_bsmu = it_stnmap_range.first;

	matrix_2d* estimatedStations_stn1(&v_estimatedStations_.at(_it_bsmu->second));
	UINT32 stn1(GetBlkMatrixElemStn1(_it_bsmu->second, _it_msr));
	
	// Get estimated station coordinates matrix and index of the 
	// station within the matrix for station 2
	it_stnmap_range = equal_range(v_blockStationsMapUnique_.begin(), v_blockStationsMapUnique_.end(),
		(*_it_msr)->station2, CompareBlockStationMapUnique_Station<UINT32, u32u32_uint32_pair>());
	if (it_stnmap_range.first == it_stnmap_range.second)
	{
		//TRACE("%d not found\n", (*_it_msr)->station2);
		std::stringstream ss;
		ss << "UpdateIgnoredMeasurements(): Station " << (*_it_msr)->station1 <<
			" was not found in the station map." << std::endl;
		SignalExceptionAdjustment(ss.str(), 0);
	}
	// Get the index of the station in the map
	_it_bsmu = it_stnmap_range.first;

	matrix_2d* estimatedStations_stn2(&v_estimatedStations_.at(_it_bsmu->second));
	UINT32 stn2(GetBlkMatrixElemStn2(_it_bsmu->second, _it_msr));
	
	// Get estimated station coordinates matrix and index of the 
	// station within the matrix for station 3
	it_stnmap_range = equal_range(v_blockStationsMapUnique_.begin(), v_blockStationsMapUnique_.end(),
		(*_it_msr)->station3, CompareBlockStationMapUnique_Station<UINT32, u32u32_uint32_pair>());
	if (it_stnmap_range.first == it_stnmap_range.second)
	{
		//TRACE("%d not found\n", (*_it_msr)->station3);
		std::stringstream ss;
		ss << "UpdateIgnoredMeasurements(): Station " << (*_it_msr)->station3 <<
			" was not found in the station map." << std::endl;
		SignalExceptionAdjustment(ss.str(), 0);
	}

	// Get the index of the station in the map
	_it_bsmu = it_stnmap_range.first;

	matrix_2d* estimatedStations_stn3(&v_estimatedStations_.at(_it_bsmu->second));
	UINT32 stn3(GetBlkMatrixElemStn3(_it_bsmu->second, _it_msr));

	it_vstn_t_const stn1_it(bstBinaryRecords_.begin() + (*_it_msr)->station1);

	double direction12, direction13, local_12e, local_12n, local_13e, local_13n;

	// compute angle 1 -> 2 -> 3 from estimated coordinates
	(*_it_msr)->measAdj = (HorizontalAngle(
		estimatedStations_stn1->get(stn1, 0),		// X1
		estimatedStations_stn1->get(stn1 + 1, 0),		// Y1
		estimatedStations_stn1->get(stn1 + 2, 0),		// Z1
		estimatedStations_stn2->get(stn2, 0), 		// X2
		estimatedStations_stn2->get(stn2 + 1, 0),		// Y2
		estimatedStations_stn2->get(stn2 + 2, 0),		// Z2
		estimatedStations_stn3->get(stn3, 0), 		// X3
		estimatedStations_stn3->get(stn3 + 1, 0),		// Y3
		estimatedStations_stn3->get(stn3 + 2, 0),		// Z3
		stn1_it->currentLatitude,
		stn1_it->currentLongitude,
		&direction12, &direction13,
		&local_12e, &local_12n, &local_13e, &local_13n));

	// deflections available?
	if (fabs(stn1_it->verticalDef) > E4_SEC_DEFLECTION || fabs(stn1_it->meridianDef) > E4_SEC_DEFLECTION)
	{
		it_vstn_t stn2_it(bstBinaryRecords_.begin() + (*_it_msr)->station2);
		it_vstn_t stn3_it(bstBinaryRecords_.begin() + (*_it_msr)->station3);

		/////////////////////////////////////////////////////////////////////////////////
		// Angles (observed or derived from directions) must be corrected for deflection 
		// of the vertical via "Laplace correction".  This correction requires zenith 
		// distance (zenith12, zenith13) and geodetic azimuth (direction12, direction13), 
		// both of which must be computed from coordinates.

		////////////////////////////////////////////////////////////////////////////
		// Compute zenith distance 1 -> 2
		double zenith12(ZenithDistance<double>(
			estimatedStations_stn1->get(stn1, 0),			// X1
			estimatedStations_stn1->get(stn1 + 1, 0),			// Y1
			estimatedStations_stn1->get(stn1 + 2, 0),			// Z1
			estimatedStations_stn2->get(stn2, 0), 			// X2
			estimatedStations_stn2->get(stn2 + 1, 0),			// Y2
			estimatedStations_stn2->get(stn2 + 2, 0),			// Z2
			stn1_it->currentLatitude,
			stn1_it->currentLongitude,
			stn2_it->currentLatitude,
			stn2_it->currentLongitude,
			(*_it_msr)->term3,							// instrument height
			(*_it_msr)->term4));						// target height

		////////////////////////////////////////////////////////////////////////////
		// Compute zenith distance 1 -> 3
		double zenith13(ZenithDistance<double>(
			estimatedStations_stn1->get(stn1, 0),			// X1
			estimatedStations_stn1->get(stn1 + 1, 0),			// Y1
			estimatedStations_stn1->get(stn1 + 2, 0),			// Z1
			estimatedStations_stn3->get(stn3, 0), 			// X2
			estimatedStations_stn3->get(stn3 + 1, 0),			// Y2
			estimatedStations_stn3->get(stn3 + 2, 0),			// Z2
			stn1_it->currentLatitude,
			stn1_it->currentLongitude,
			stn3_it->currentLatitude,
			stn3_it->currentLongitude,
			(*_it_msr)->term3,							// instrument height
			(*_it_msr)->term4));						// target height

		// Laplace correction 1 -> 2 -> 3
		(*_it_msr)->preAdjCorr = HzAngleDeflectionCorrection<double>(
			direction12,								// geodetic azimuth 1 -> 2
			zenith12,									// zenith distance 1 -> 2
			direction13,								// geodetic azimuth 1 -> 3
			zenith13,									// zenith distance 1 -> 3
			stn1_it->verticalDef,						// deflection in prime vertical
			stn1_it->meridianDef);						// deflection in prime meridian
	}
	else
		(*_it_msr)->preAdjCorr = 0.;

	// compute adjustment correction
	(*_it_msr)->measAdj += (*_it_msr)->preAdjCorr;
	// compute adjustment correction
	(*_it_msr)->measCorr = (*_it_msr)->measAdj - (*_it_msr)->preAdjMeas;
}

void dna_adjust::UpdateIgnoredMeasurements_B(pit_vmsr_t _it_msr, bool storeOriginalMeasurement)
{
	// initialise measurement (on the first adjustment only!)
	if (storeOriginalMeasurement)
		(*_it_msr)->preAdjMeas = (*_it_msr)->term1;

	// Compute the geodetic azimuth
	UpdateIgnoredMeasurements_BK(_it_msr);

	// compute adjustment correction
	(*_it_msr)->measCorr = (*_it_msr)->measAdj - (*_it_msr)->preAdjMeas;
}

void dna_adjust::UpdateIgnoredMeasurements_BK(pit_vmsr_t _it_msr)
{	
	// Use v_blockStationsMapUnique_ to get the correct block for each 
	// station in the measurement
	_it_u32u32_uint32_pair _it_bsmu;

	// Get estimated station coordinates matrix and index of the 
	// station within the matrix for station 1
	_it_bsmu = v_blockStationsMapUnique_.begin() + (*_it_msr)->station1;
	matrix_2d* estimatedStations_stn1(&v_estimatedStations_.at(_it_bsmu->second));
	UINT32 stn1(GetBlkMatrixElemStn1(_it_bsmu->second, _it_msr));
	
	// Get estimated station coordinates matrix and index of the 
	// station within the matrix for station 2
	_it_bsmu = v_blockStationsMapUnique_.begin() + (*_it_msr)->station2;
	matrix_2d* estimatedStations_stn2(&v_estimatedStations_.at(_it_bsmu->second));
	UINT32 stn2(GetBlkMatrixElemStn2(_it_bsmu->second, _it_msr));

	it_vstn_t_const stn1_it(bstBinaryRecords_.begin() + (*_it_msr)->station1);
	it_vstn_t_const stn2_it(bstBinaryRecords_.begin() + (*_it_msr)->station2);

	double local_12e, local_12n;

	// compute bearing from estimated coordinates
	(*_it_msr)->measAdj = (Direction(
		estimatedStations_stn1->get(stn1, 0),		// X1
		estimatedStations_stn1->get(stn1 + 1, 0),	// Y1
		estimatedStations_stn1->get(stn1 + 2, 0),	// Z1
		estimatedStations_stn2->get(stn2, 0), 		// X2
		estimatedStations_stn2->get(stn2 + 1, 0),	// Y2
		estimatedStations_stn2->get(stn2 + 2, 0),	// Z2
		stn1_it->currentLatitude,
		stn1_it->currentLongitude,
		&local_12e, &local_12n));
}
	

void dna_adjust::UpdateIgnoredMeasurements_C(pit_vmsr_t _it_msr, bool storeOriginalMeasurement)
{
	// initialise measurement (on the first adjustment only!)
	if (storeOriginalMeasurement)
		(*_it_msr)->preAdjMeas = (*_it_msr)->term1;

	// As a C measurement is a direct vector between two points in the cartesian frame,
	// there will be no pre adjustment correction
	(*_it_msr)->preAdjCorr = 0.;

	// Compute the chord distance
	UpdateIgnoredMeasurements_CEM(_it_msr);

	// compute adjustment correction
	(*_it_msr)->measCorr = (*_it_msr)->measAdj - (*_it_msr)->preAdjMeas;
}
	

void dna_adjust::UpdateIgnoredMeasurements_CEM(pit_vmsr_t _it_msr)
{	

	// Use v_blockStationsMapUnique_ to get the correct block for each 
	// station in the measurement
	_it_u32u32_uint32_pair _it_bsmu;

	// Get estimated station coordinates matrix and index of the 
	// station within the matrix for station 1
	_it_bsmu = v_blockStationsMapUnique_.begin() + (*_it_msr)->station1;
	matrix_2d* estimatedStations_stn1(&v_estimatedStations_.at(_it_bsmu->second));
	UINT32 stn1(GetBlkMatrixElemStn1(_it_bsmu->second, _it_msr));
	
	// Get estimated station coordinates matrix and index of the 
	// station within the matrix for station 2
	_it_bsmu = v_blockStationsMapUnique_.begin() + (*_it_msr)->station2;
	matrix_2d* estimatedStations_stn2(&v_estimatedStations_.at(_it_bsmu->second));
	UINT32 stn2(GetBlkMatrixElemStn2(_it_bsmu->second, _it_msr));

	it_vstn_t_const stn1_it(bstBinaryRecords_.begin() + (*_it_msr)->station1);
	it_vstn_t_const stn2_it(bstBinaryRecords_.begin() + (*_it_msr)->station2);

	double dX, dY, dZ;

	// calculate chord distance
	(*_it_msr)->measAdj = (EllipsoidChordDistance<double>(
		estimatedStations_stn1->get(stn1, 0),
		estimatedStations_stn1->get(stn1 + 1, 0),
		estimatedStations_stn1->get(stn1 + 2, 0),
		estimatedStations_stn2->get(stn2, 0),
		estimatedStations_stn2->get(stn2 + 1, 0),
		estimatedStations_stn2->get(stn2 + 2, 0),
		stn1_it->currentLatitude,
		stn2_it->currentLatitude,
		stn1_it->currentHeight,
		stn2_it->currentHeight,
		&dX, &dY, &dZ,
		datum_.GetEllipsoidRef()));
}


void dna_adjust::UpdateIgnoredMeasurements_D(pit_vmsr_t _it_msr, bool storeOriginalMeasurement)
{
	it_vmsr_t _it_msr_first(*_it_msr);
	UINT32 a, angle_count((*_it_msr)->vectorCount1 - 1);		// number of directions excluding the RO
	
	vmsr_t angleRec;
	angleRec.push_back(*(*_it_msr));
	it_vmsr_t it_angle(angleRec.begin());

	double previousDirection((*_it_msr)->term1);
	double previousVariance((*_it_msr)->term2);

	(*_it_msr)++;

	// set derived angle, variance and covariance to the binary records
	// term1 = measured direction
	// term2 = variance (direction)
	// term3 = instrument height (not used)
	// term4 = target height (not used)
	// scale1 = derived angle corrected for deflection of the vertical
	// scale2 = variance (angle)
	// scale3 = covariance (angle) - for the context of vmsr_t angleRec only, so as to 
	//          properly form the normals from covariances formed from directions SDs
	// preAdjMeas = original derived angle

	if (projectSettings_.g.verbose > 6)
		debug_file << "Reduced angles from raw directions: ";

	try
	{
		for (a=0; a<angle_count; ++a)
		{
			it_angle->station3 = (*_it_msr)->station2;

			if (storeOriginalMeasurement)
			{				
				// Derive the angle from two directions
				it_angle->term1 = (*_it_msr)->term1 - previousDirection;
				if (it_angle->term1 < 0)
					it_angle->term1 += TWO_PI;
				if (it_angle->term1 > TWO_PI)
					it_angle->term1 -= TWO_PI;
				// Derive the variance of the angle
				(*_it_msr)->scale2 = previousVariance + (*_it_msr)->term2;
			}
			else
				it_angle->preAdjMeas = (*_it_msr)->scale1;
			
			UpdateIgnoredMeasurements_A(&it_angle, storeOriginalMeasurement);
						
			// apply correction for deflections in the vertical
			(*_it_msr)->scale1 = it_angle->preAdjMeas;
			// compute adjustment correction
			(*_it_msr)->measCorr = it_angle->measCorr;
			// Update correction for deflection of the vertical
			(*_it_msr)->preAdjCorr = it_angle->preAdjCorr;

			if (a + 1 == angle_count)
				break;

			if (storeOriginalMeasurement)
			{
				previousDirection = (*_it_msr)->term1;
				previousVariance = (*_it_msr)->term2;
			}

			// prepare for next angle
			angleRec.push_back(*(*_it_msr));
			it_angle = angleRec.end() - 1;

			(*_it_msr)++;
		}

	}
	catch (...) {

		// Print error message to adj file and throw exception
		std::stringstream ss;
		ss << "UpdateIgnoredMeasurements_D(): An error was encountered whilst" << std::endl <<
			"  calculating ignored measurement details" << std::endl;
		SignalExceptionAdjustment(ss.str(), 0);
	}
	
}
	

void dna_adjust::UpdateIgnoredMeasurements_E(pit_vmsr_t _it_msr, bool storeOriginalMeasurement)
{	
	// initialise measurement (on the first adjustment only!)
	if (storeOriginalMeasurement)
		(*_it_msr)->preAdjMeas = (*_it_msr)->term1;

	// Compute the chord distance
	UpdateIgnoredMeasurements_CEM(_it_msr);

	// Use v_blockStationsMapUnique_ to get the correct block for each 
	// station in the measurement
	_it_u32u32_uint32_pair _it_bsmu;

	// Get estimated station coordinates matrix and index of the 
	// station within the matrix for station 1
	_it_bsmu = v_blockStationsMapUnique_.begin() + (*_it_msr)->station1;
	matrix_2d* estimatedStations_stn1(&v_estimatedStations_.at(_it_bsmu->second));
	UINT32 stn1(GetBlkMatrixElemStn1(_it_bsmu->second, _it_msr));
	
	// Get estimated station coordinates matrix and index of the 
	// station within the matrix for station 2
	_it_bsmu = v_blockStationsMapUnique_.begin() + (*_it_msr)->station2;
	matrix_2d* estimatedStations_stn2(&v_estimatedStations_.at(_it_bsmu->second));
	UINT32 stn2(GetBlkMatrixElemStn2(_it_bsmu->second, _it_msr));

	it_vstn_t_const stn1_it(bstBinaryRecords_.begin() + (*_it_msr)->station1);
	it_vstn_t_const stn2_it(bstBinaryRecords_.begin() + (*_it_msr)->station2);

	// Compute Ellipsoid arc from Ellipsoid chord
	double ellipsoid_arc = EllipsoidChordtoEllipsoidArc<double>(
		(*_it_msr)->measAdj,		// use Ellipsoid chord computed by UpdateIgnoredMeasurements_CEM
		estimatedStations_stn1->get(stn1, 0),
		estimatedStations_stn1->get(stn1 + 1, 0),
		estimatedStations_stn1->get(stn1 + 2, 0),
		estimatedStations_stn2->get(stn2, 0),
		estimatedStations_stn2->get(stn2 + 1, 0),
		estimatedStations_stn2->get(stn2 + 2, 0),
		stn1_it->currentLatitude,
		stn1_it->currentLongitude,
		stn2_it->currentLatitude,
		datum_.GetEllipsoidRef());

	// compute correction from arc to chord
	(*_it_msr)->preAdjCorr = ellipsoid_arc - (*_it_msr)->measAdj;
	// update adjusted measurement
	(*_it_msr)->measAdj = ellipsoid_arc;
	// compute adjustment correction
	(*_it_msr)->measCorr = (*_it_msr)->measAdj - (*_it_msr)->preAdjMeas;
}
	

void dna_adjust::UpdateIgnoredMeasurements_G(pit_vmsr_t _it_msr, bool storeOriginalMeasurement)
{
	UpdateIgnoredMeasurements_GX(_it_msr, storeOriginalMeasurement);
}

void dna_adjust::UpdateIgnoredMeasurements_GX(pit_vmsr_t _it_msr, bool storeOriginalMeasurement)
{
	// Use v_blockStationsMapUnique_ to get the correct block for each 
	// station in the measurement
	_it_u32u32_uint32_pair _it_bsmu;

	// Get estimated station coordinates matrix and index of the 
	// station within the matrix for station 1
	_it_bsmu = v_blockStationsMapUnique_.begin() + (*_it_msr)->station1;
	matrix_2d* estimatedStations_stn1(&v_estimatedStations_.at(_it_bsmu->second));
	UINT32 stn1(GetBlkMatrixElemStn1(_it_bsmu->second, _it_msr));

	// Get estimated station coordinates matrix and index of the 
	// station within the matrix for station 2
	_it_bsmu = v_blockStationsMapUnique_.begin() + (*_it_msr)->station2;
	matrix_2d* estimatedStations_stn2(&v_estimatedStations_.at(_it_bsmu->second));
	UINT32 stn2(GetBlkMatrixElemStn2(_it_bsmu->second, _it_msr));

	// initialise measurement (on the first adjustment only!)
	if (storeOriginalMeasurement)
		// X element
		(*_it_msr)->preAdjMeas = (*_it_msr)->term1;
	// update adjusted measurement
	(*_it_msr)->measAdj = (estimatedStations_stn2->get(stn2, 0) - estimatedStations_stn1->get(stn1, 0));
	// compute adjustment correction
	(*_it_msr)->measCorr = (*_it_msr)->measAdj - (*_it_msr)->preAdjMeas;

	// move to Y element
	(*_it_msr)++;
	if (storeOriginalMeasurement)
		(*_it_msr)->preAdjMeas = (*_it_msr)->term1;
	// update adjusted measurement
	(*_it_msr)->measAdj = (estimatedStations_stn2->get(stn2+1, 0) - estimatedStations_stn1->get(stn1+1, 0));
	// compute adjustment correction
	(*_it_msr)->measCorr = (*_it_msr)->measAdj - (*_it_msr)->preAdjMeas;

	// move to Z element
	(*_it_msr)++;
	if (storeOriginalMeasurement)
		(*_it_msr)->preAdjMeas = (*_it_msr)->term1;
	// update adjusted measurement
	(*_it_msr)->measAdj = (estimatedStations_stn2->get(stn2+2, 0) - estimatedStations_stn1->get(stn1+2, 0));
	// compute adjustment correction
	(*_it_msr)->measCorr = (*_it_msr)->measAdj - (*_it_msr)->preAdjMeas;
}
	

void dna_adjust::UpdateIgnoredMeasurements_H(pit_vmsr_t _it_msr, bool storeOriginalMeasurement)
{
	// initialise measurement (on the first adjustment only!)
	if (storeOriginalMeasurement)
		(*_it_msr)->preAdjMeas = (*_it_msr)->term1;

	// Compute the ellipsoid height
	UpdateIgnoredMeasurements_HR(_it_msr);

	it_vstn_t_const stn1_it(bstBinaryRecords_.begin() + (*_it_msr)->station1);

	// N value available?
	if (fabs(stn1_it->geoidSep) > PRECISION_1E4)
	{
		// get ellipsoid - geoid separation
		(*_it_msr)->preAdjCorr = stn1_it->geoidSep;
	}
	else
		(*_it_msr)->preAdjCorr = 0.;

	(*_it_msr)->measAdj -= (*_it_msr)->preAdjCorr;
	// compute adjustment correction
	(*_it_msr)->measCorr = (*_it_msr)->measAdj - (*_it_msr)->preAdjMeas;
}
	

void dna_adjust::UpdateIgnoredMeasurements_HR(pit_vmsr_t _it_msr)
{
	// Use v_blockStationsMapUnique_ to get the correct block for each 
	// station in the measurement
	_it_u32u32_uint32_pair _it_bsmu;

	// Get estimated station coordinates matrix and index of the 
	// station within the matrix for station 1
	_it_bsmu = v_blockStationsMapUnique_.begin() + (*_it_msr)->station1;
	matrix_2d* estimatedStations(&v_estimatedStations_.at(_it_bsmu->second));
	UINT32 stn1(GetBlkMatrixElemStn1(_it_bsmu->second, _it_msr));

	it_vstn_t_const stn1_it(bstBinaryRecords_.begin() + (*_it_msr)->station1);

	// Zn is the z coordinate element of the point on the z-axis 
	// which intersects with the the normal at the given Latitude
	double nu1, Zn1;

	// compute the ellipsoid height height
	(*_it_msr)->measAdj = (EllipsoidHeight<double>(
		estimatedStations->get(stn1, 0),
		estimatedStations->get(stn1 + 1, 0),
		estimatedStations->get(stn1 + 2, 0),
		stn1_it->currentLatitude,
		&nu1, &Zn1,
		datum_.GetEllipsoidRef()));
}


void dna_adjust::UpdateIgnoredMeasurements_I(pit_vmsr_t _it_msr, bool storeOriginalMeasurement)
{
	// initialise measurement (on the first adjustment only!)
	if (storeOriginalMeasurement)
		(*_it_msr)->preAdjMeas = (*_it_msr)->term1;

	// calculate geodetic latitude
	UpdateIgnoredMeasurements_IP(_it_msr);

	it_vstn_t_const stn1_it(bstBinaryRecords_.begin() + (*_it_msr)->station1);

	// deflections available?
	if (fabs(stn1_it->meridianDef) > E4_SEC_DEFLECTION)
		// deflection in the prime meridian
		(*_it_msr)->preAdjCorr = stn1_it->meridianDef;						
	else
		(*_it_msr)->preAdjCorr = 0.;

	// apply deflection correction
	(*_it_msr)->measAdj += (*_it_msr)->preAdjCorr;						
	// compute adjustment correction
	(*_it_msr)->measCorr = (*_it_msr)->measAdj - (*_it_msr)->preAdjMeas;
}
	

void dna_adjust::UpdateIgnoredMeasurements_IP(pit_vmsr_t _it_msr)
{
	// Use v_blockStationsMapUnique_ to get the correct block for each 
	// station in the measurement
	_it_u32u32_uint32_pair _it_bsmu;

	// Get estimated station coordinates matrix and index of the 
	// station within the matrix for station 1
	_it_bsmu = v_blockStationsMapUnique_.begin() + (*_it_msr)->station1;
	matrix_2d* estimatedStations(&v_estimatedStations_.at(_it_bsmu->second));
	UINT32 stn1(GetBlkMatrixElemStn1(_it_bsmu->second, _it_msr));

	// compute the geodetic latitude
	(*_it_msr)->measAdj = (CartToLat<double>(
		estimatedStations->get(stn1, 0),			// X1
		estimatedStations->get(stn1 + 1, 0),		// Y1
		estimatedStations->get(stn1 + 2, 0),		// Z1
		datum_.GetEllipsoidRef()));
}


void dna_adjust::UpdateIgnoredMeasurements_J(pit_vmsr_t _it_msr, bool storeOriginalMeasurement)
{
	// initialise measurement (on the first adjustment only!)
	if (storeOriginalMeasurement)
		(*_it_msr)->preAdjMeas = (*_it_msr)->term1;

	// calculate geodetic longitude
	UpdateIgnoredMeasurements_JQ(_it_msr);

	it_vstn_t_const stn1_it(bstBinaryRecords_.begin() + (*_it_msr)->station1);

	// deflections available?
	if (fabs(stn1_it->verticalDef) > E4_SEC_DEFLECTION)
		// deflection in the prime vertical
		(*_it_msr)->preAdjCorr =
		stn1_it->verticalDef / cos(stn1_it->currentLatitude);		// sec(a) = 1/cos(a)
	else
		(*_it_msr)->preAdjCorr = 0.;
	
	// apply deflection correction
	(*_it_msr)->measAdj += (*_it_msr)->preAdjCorr;
	// compute adjustment correction
	(*_it_msr)->measCorr = (*_it_msr)->measAdj - (*_it_msr)->preAdjMeas;
	
}
	

void dna_adjust::UpdateIgnoredMeasurements_JQ(pit_vmsr_t _it_msr)
{
	// Use v_blockStationsMapUnique_ to get the correct block for each 
	// station in the measurement
	_it_u32u32_uint32_pair _it_bsmu;

	// Get estimated station coordinates matrix and index of the 
	// station within the matrix for station 1
	_it_bsmu = v_blockStationsMapUnique_.begin() + (*_it_msr)->station1;
	matrix_2d* estimatedStations(&v_estimatedStations_.at(_it_bsmu->second));
	UINT32 stn1(GetBlkMatrixElemStn1(_it_bsmu->second, _it_msr));
	
	double latitude, longitude, ellipsoidHeight;

	// compute the geodetic longitude
	CartToGeo<double>(
		estimatedStations->get(stn1, 0),			// X1
		estimatedStations->get(stn1 + 1, 0),		// Y1
		estimatedStations->get(stn1 + 2, 0),		// Z1
		&latitude,
		&longitude,
		&ellipsoidHeight,
		datum_.GetEllipsoidRef());

	(*_it_msr)->measAdj = longitude;
}

void dna_adjust::UpdateIgnoredMeasurements_K(pit_vmsr_t _it_msr, bool storeOriginalMeasurement)
{
	// initialise measurement (on the first adjustment only!)
	if (storeOriginalMeasurement)
		(*_it_msr)->preAdjMeas = (*_it_msr)->term1;

	// Compute the geodetic azimuth
	UpdateIgnoredMeasurements_BK(_it_msr);

	// Use v_blockStationsMapUnique_ to get the correct block for each 
	// station in the measurement
	_it_u32u32_uint32_pair _it_bsmu;

	// Get estimated station coordinates matrix and index of the 
	// station within the matrix for station 1
	_it_bsmu = v_blockStationsMapUnique_.begin() + (*_it_msr)->station1;
	matrix_2d* estimatedStations_stn1(&v_estimatedStations_.at(_it_bsmu->second));
	UINT32 stn1(GetBlkMatrixElemStn1(_it_bsmu->second, _it_msr));
	
	// Get estimated station coordinates matrix and index of the 
	// station within the matrix for station 2
	_it_bsmu = v_blockStationsMapUnique_.begin() + (*_it_msr)->station2;
	matrix_2d* estimatedStations_stn2(&v_estimatedStations_.at(_it_bsmu->second));
	UINT32 stn2(GetBlkMatrixElemStn2(_it_bsmu->second, _it_msr));
	
	it_vstn_t_const stn1_it(bstBinaryRecords_.begin() + (*_it_msr)->station1);
	it_vstn_t_const stn2_it(bstBinaryRecords_.begin() + (*_it_msr)->station2);

	// deflections available?
	if ((*_it_msr)->measType == 'K' &&						// Astro
		(fabs(stn1_it->verticalDef) > E4_SEC_DEFLECTION ||	// deflections available?
			fabs(stn1_it->meridianDef) > E4_SEC_DEFLECTION))
	{
		////////////////////////////////////////////////////////////////////////////
		// Astronomic azimuths must be corrected for deflection of the vertical via
		// "Laplace correction".  This correction requires zenith distance and geodetic 
		// azimuth (comp_msr), both of which must be computed from coordinates.

		////////////////////////////////////////////////////////////////////////////
		// Compute zenith distance
		double zenith(ZenithDistance<double>(
			estimatedStations_stn1->get(stn1, 0),					// X1
			estimatedStations_stn1->get(stn1 + 1, 0),				// Y1
			estimatedStations_stn1->get(stn1 + 2, 0),				// Z1
			estimatedStations_stn2->get(stn2, 0), 					// X2
			estimatedStations_stn2->get(stn2 + 1, 0),				// Y2
			estimatedStations_stn2->get(stn2 + 2, 0),				// Z2
			stn1_it->currentLatitude,
			stn1_it->currentLongitude,
			stn2_it->currentLatitude,
			stn2_it->currentLongitude,
			(*_it_msr)->term3,									// instrument height
			(*_it_msr)->term4));								// target height

		// Compute pre adjustment correction
		(*_it_msr)->preAdjCorr = LaplaceCorrection<double>(		// Laplace correction
			(*_it_msr)->measAdj,								// geodetic azimuth
			zenith,												// zenith distance
			stn1_it->verticalDef,								// deflection in prime vertical
			stn1_it->meridianDef,								// deflection in prime meridian
			stn1_it->currentLatitude);
	}
	else
		(*_it_msr)->preAdjCorr = 0.;

	// apply correction for deflections in the vertical
	(*_it_msr)->measAdj += (*_it_msr)->preAdjCorr;
	// compute adjustment correction
	(*_it_msr)->measCorr = (*_it_msr)->measAdj - (*_it_msr)->preAdjMeas;
}

void dna_adjust::UpdateIgnoredMeasurements_L(pit_vmsr_t _it_msr, bool storeOriginalMeasurement)
{
	// initialise measurement (on the first adjustment only!)
	if (storeOriginalMeasurement)
		(*_it_msr)->preAdjMeas = (*_it_msr)->term1;

	// Use v_blockStationsMapUnique_ to get the correct block for each 
	// station in the measurement
	_it_u32u32_uint32_pair _it_bsmu;

	// Get estimated station coordinates matrix and index of the 
	// station within the matrix for station 1
	_it_bsmu = v_blockStationsMapUnique_.begin() + (*_it_msr)->station1;
	matrix_2d* estimatedStations_stn1(&v_estimatedStations_.at(_it_bsmu->second));
	UINT32 stn1(GetBlkMatrixElemStn1(_it_bsmu->second, _it_msr));
	
	// Get estimated station coordinates matrix and index of the 
	// station within the matrix for station 2
	_it_bsmu = v_blockStationsMapUnique_.begin() + (*_it_msr)->station2;
	matrix_2d* estimatedStations_stn2(&v_estimatedStations_.at(_it_bsmu->second));
	UINT32 stn2(GetBlkMatrixElemStn2(_it_bsmu->second, _it_msr));

	it_vstn_t_const stn1_it(bstBinaryRecords_.begin() + (*_it_msr)->station1);
	it_vstn_t_const stn2_it(bstBinaryRecords_.begin() + (*_it_msr)->station2);

	// Zn is the z coordinate element of the point on the z-axis 
	// which intersects with the the normal at the given Latitude
	double h1, h2, nu1, nu2, Zn1, Zn2;

	// calculated diff height
	(*_it_msr)->measAdj = (EllipsoidHeightDifference<double>(
		estimatedStations_stn1->get(stn1, 0),
		estimatedStations_stn1->get(stn1 + 1, 0),
		estimatedStations_stn1->get(stn1 + 2, 0),
		estimatedStations_stn2->get(stn2, 0),
		estimatedStations_stn2->get(stn2 + 1, 0),
		estimatedStations_stn2->get(stn2 + 2, 0),
		stn1_it->currentLatitude,
		stn2_it->currentLatitude,
		&h1, &h2, &nu1, &nu2, &Zn1, &Zn2,
		datum_.GetEllipsoidRef()));

	// N value available?
	if (fabs(stn1_it->geoidSep) > PRECISION_1E4 ||
		fabs(stn2_it->geoidSep) > PRECISION_1E4)
		// Compute ellipsoid-geoid separation correction
		(*_it_msr)->preAdjCorr = stn2_it->geoidSep - stn1_it->geoidSep;
	else
		(*_it_msr)->preAdjCorr = 0.;

	// apply the ellipsoid-geoid separation correction
	(*_it_msr)->measAdj -= (*_it_msr)->preAdjCorr;
	// compute adjustment correction
	(*_it_msr)->measCorr = (*_it_msr)->measAdj - (*_it_msr)->preAdjMeas;
}
	

void dna_adjust::UpdateIgnoredMeasurements_M(pit_vmsr_t _it_msr, bool storeOriginalMeasurement)
{
	// initialise measurement (on the first adjustment only!)
	if (storeOriginalMeasurement)
		(*_it_msr)->preAdjMeas = (*_it_msr)->term1;
	
	// Compute the chord distance
	UpdateIgnoredMeasurements_CEM(_it_msr);
		
	it_vstn_t_const stn1_it(bstBinaryRecords_.begin() + (*_it_msr)->station1);
	it_vstn_t_const stn2_it(bstBinaryRecords_.begin() + (*_it_msr)->station2);

	// Compute MSL arc from Ellipsoid chord
	double msl_arc = EllipsoidChordtoMSLArc<double>(
		(*_it_msr)->measAdj,			// use Ellipsoid chord computed by UpdateIgnoredMeasurements_CEM
		stn1_it->currentLatitude, stn2_it->currentLatitude,
		stn1_it->geoidSep, stn2_it->geoidSep,
		datum_.GetEllipsoidRef());

	// compute correction from arc to chord
	(*_it_msr)->preAdjCorr = msl_arc - (*_it_msr)->measAdj;
	// update adjusted measurement
	(*_it_msr)->measAdj = msl_arc;
	// compute adjustment correction
	(*_it_msr)->measCorr = (*_it_msr)->measAdj - (*_it_msr)->preAdjMeas;
}
	

void dna_adjust::UpdateIgnoredMeasurements_P(pit_vmsr_t _it_msr, bool storeOriginalMeasurement)
{
	// initialise measurement (on the first adjustment only!)
	if (storeOriginalMeasurement)
		(*_it_msr)->preAdjMeas = (*_it_msr)->term1;

	UpdateIgnoredMeasurements_IP(_it_msr);

	// As a P measurement is the result of a direct conversion of cartesian coordinates,
	// there will be no pre adjustment correction
	(*_it_msr)->preAdjCorr = 0.;
	// compute adjustment correction
	(*_it_msr)->measCorr = (*_it_msr)->measAdj - (*_it_msr)->preAdjMeas;
}
	

void dna_adjust::UpdateIgnoredMeasurements_Q(pit_vmsr_t _it_msr, bool storeOriginalMeasurement)
{
	// initialise measurement (on the first adjustment only!)
	if (storeOriginalMeasurement)
		(*_it_msr)->preAdjMeas = (*_it_msr)->term1;

	UpdateIgnoredMeasurements_JQ(_it_msr);

	// As a Q measurement is the result of a direct conversion of cartesian coordinates,
	// there will be no pre adjustment correction
	(*_it_msr)->preAdjCorr = 0.;
	// compute adjustment correction
	(*_it_msr)->measCorr = (*_it_msr)->measAdj - (*_it_msr)->preAdjMeas;
}
	

void dna_adjust::UpdateIgnoredMeasurements_R(pit_vmsr_t _it_msr, bool storeOriginalMeasurement)
{
	// initialise measurement (on the first adjustment only!)
	if (storeOriginalMeasurement)
		(*_it_msr)->preAdjMeas = (*_it_msr)->term1;

	UpdateIgnoredMeasurements_HR(_it_msr);

	// As a R measurement is the result of a direct conversion of cartesian coordinates,
	// there will be no pre adjustment correction
	(*_it_msr)->preAdjCorr = 0.;
	// compute adjustment correction
	(*_it_msr)->measCorr = (*_it_msr)->measAdj - (*_it_msr)->preAdjMeas;
}
	

void dna_adjust::UpdateIgnoredMeasurements_S(pit_vmsr_t _it_msr, bool storeOriginalMeasurement)
{
	// initialise measurement (on the first adjustment only!)
	if (storeOriginalMeasurement)
		(*_it_msr)->preAdjMeas = (*_it_msr)->term1;

	// Use v_blockStationsMapUnique_ to get the correct block for each 
	// station in the measurement
	_it_pair_u32u32_uint32 it_stnmap_range;
	_it_u32u32_uint32_pair _it_bsmu;

	// Get estimated station coordinates matrix and index of the 
	// station within the matrix for station 1
	it_stnmap_range = equal_range(v_blockStationsMapUnique_.begin(), v_blockStationsMapUnique_.end(), 
		(*_it_msr)->station1, CompareBlockStationMapUnique_Station<UINT32, u32u32_uint32_pair>());
	if (it_stnmap_range.first == it_stnmap_range.second)
	{
		//TRACE("%d not found\n", (*_it_msr)->station1);
		std::stringstream ss;
		ss << "UpdateIgnoredMeasurements(): Station " << (*_it_msr)->station1  << 
			" was not found in the station map." << std::endl;
		SignalExceptionAdjustment(ss.str(), 0);
	}	

	// Get the index of the station in the map
	_it_bsmu = it_stnmap_range.first;

	matrix_2d* estimatedStations_stn1(&v_estimatedStations_.at(_it_bsmu->second));
	UINT32 stn1(GetBlkMatrixElemStn1(_it_bsmu->second, _it_msr));
	
	// Get estimated station coordinates matrix and index of the 
	// station within the matrix for station 2
	it_stnmap_range = equal_range(v_blockStationsMapUnique_.begin(), v_blockStationsMapUnique_.end(),
		(*_it_msr)->station2, CompareBlockStationMapUnique_Station<UINT32, u32u32_uint32_pair>());
	if (it_stnmap_range.first == it_stnmap_range.second)
	{
		//TRACE("%d not found\n", (*_it_msr)->station2);
		std::stringstream ss;
		ss << "UpdateIgnoredMeasurements(): Station " << (*_it_msr)->station2 <<
			" was not found in the station map." << std::endl;
		SignalExceptionAdjustment(ss.str(), 0);
	}
	
	// Get the index of the station in the map
	_it_bsmu = it_stnmap_range.first;

	matrix_2d* estimatedStations_stn2(&v_estimatedStations_.at(_it_bsmu->second));
	UINT32 stn2(GetBlkMatrixElemStn2(_it_bsmu->second, _it_msr));
	
	it_vstn_t_const stn1_it(bstBinaryRecords_.begin() + (*_it_msr)->station1);

	// compute dX, dY, dZ for instrument height (ih) and target height (th)
	double dXih, dYih, dZih, dXth, dYth, dZth;
	CartesianElementsFromInstrumentHeight((*_it_msr)->term3,				// instrument height
		&dXih, &dYih, &dZih,
		stn1_it->currentLatitude,
		stn1_it->currentLongitude);
	CartesianElementsFromInstrumentHeight((*_it_msr)->term4,				// target height
		&dXth, &dYth, &dZth,
		stn1_it->currentLatitude,
		stn1_it->currentLongitude);

	// compute distance between instrument and target, taking into consideration
	// instrument height, target height, and vector between stations
	double dX(estimatedStations_stn2->get(stn2, 0)     - estimatedStations_stn1->get(stn1, 0) + dXth - dXih);
	double dY(estimatedStations_stn2->get(stn2 + 1, 0) - estimatedStations_stn1->get(stn1 + 1, 0) + dYth - dYih);
	double dZ(estimatedStations_stn2->get(stn2 + 2, 0) - estimatedStations_stn1->get(stn1 + 2, 0) + dZth - dZih);

	// calculated distance
	(*_it_msr)->measAdj = (magnitude(dX, dY, dZ));
	// As an S measurement is the result of a direct computation from cartesian coordinates,
	// there will be no pre adjustment correction
	(*_it_msr)->preAdjCorr = 0.;
	// compute adjustment correction
	(*_it_msr)->measCorr = (*_it_msr)->measAdj - (*_it_msr)->preAdjMeas;
}
	
// Zenith distance
//
void dna_adjust::UpdateIgnoredMeasurements_V(pit_vmsr_t _it_msr, bool storeOriginalMeasurement)
{
	// initialise measurement (on the first adjustment only!)
	if (storeOriginalMeasurement)
		(*_it_msr)->preAdjMeas = (*_it_msr)->term1;

	// Use v_blockStationsMapUnique_ to get the correct block for each 
	// station in the measurement
	_it_u32u32_uint32_pair _it_bsmu;

	// Get estimated station coordinates matrix and index of the 
	// station within the matrix for station 1
	_it_bsmu = v_blockStationsMapUnique_.begin() + (*_it_msr)->station1;
	matrix_2d* estimatedStations_stn1(&v_estimatedStations_.at(_it_bsmu->second));
	UINT32 stn1(GetBlkMatrixElemStn1(_it_bsmu->second, _it_msr));

	// Get estimated station coordinates matrix and index of the 
	// station within the matrix for station 2
	_it_bsmu = v_blockStationsMapUnique_.begin() + (*_it_msr)->station2;
	matrix_2d* estimatedStations_stn2(&v_estimatedStations_.at(_it_bsmu->second));
	UINT32 stn2(GetBlkMatrixElemStn2(_it_bsmu->second, _it_msr));

	it_vstn_t_const stn1_it(bstBinaryRecords_.begin() + (*_it_msr)->station1);
	it_vstn_t_const stn2_it(bstBinaryRecords_.begin() + (*_it_msr)->station2);

	double local_12e, local_12n, local_12up;

	// compute zenith distance from estimated coordinates
	(*_it_msr)->measAdj = (ZenithDistance(
		estimatedStations_stn1->get(stn1, 0),					// X1
		estimatedStations_stn1->get(stn1 + 1, 0),					// Y1
		estimatedStations_stn1->get(stn1 + 2, 0),					// Z1
		estimatedStations_stn2->get(stn2, 0), 					// X2
		estimatedStations_stn2->get(stn2 + 1, 0),					// Y2
		estimatedStations_stn2->get(stn2 + 2, 0),					// Z2
		stn1_it->currentLatitude,
		stn1_it->currentLongitude,
		stn2_it->currentLatitude,
		stn2_it->currentLongitude,
		(*_it_msr)->term3,									// instrument height
		(*_it_msr)->term4,									// target height
		&local_12e,											// local_12e, ..12n, ..12up represent
		&local_12n,											// the geometric difference between
		&local_12up));										// station1 and station2

	// deflections available?
	if (fabs(stn1_it->verticalDef) > E4_SEC_DEFLECTION || fabs(stn1_it->meridianDef) > E4_SEC_DEFLECTION)
	{
		////////////////////////////////////////////////////////////////////////////
		// Correct for deflections in the vertical
		// 1. compute bearing from estimated coordinates
		double azimuth(Direction(
			estimatedStations_stn1->get(stn1, 0),		// X1
			estimatedStations_stn1->get(stn1 + 1, 0),		// Y1
			estimatedStations_stn1->get(stn1 + 2, 0),		// Z1
			estimatedStations_stn2->get(stn2, 0), 		// X2
			estimatedStations_stn2->get(stn2 + 1, 0),		// Y2
			estimatedStations_stn2->get(stn2 + 2, 0),		// Z2
			stn1_it->currentLatitude,
			stn1_it->currentLongitude));

		// 2. Compute correction
		(*_it_msr)->preAdjCorr = ZenithDeflectionCorrection<double>(		// Correction to vertical angle for deflection of vertical
			azimuth,														// geodetic azimuth
			stn1_it->verticalDef,											// deflection in prime vertical
			stn1_it->meridianDef);											// deflection in prime meridian
		////////////////////////////////////////////////////////////////////////////
	}
	else
		(*_it_msr)->preAdjCorr = 0.;

	// apply correction for deflections in the vertical
	(*_it_msr)->measAdj -= (*_it_msr)->preAdjCorr;
	// compute adjustment correction
	(*_it_msr)->measCorr = (*_it_msr)->measAdj - (*_it_msr)->preAdjMeas;
}
	

void dna_adjust::UpdateIgnoredMeasurements_X(pit_vmsr_t _it_msr, bool storeOriginalMeasurement)
{
	UINT32 cluster_bsl, baseline_count((*_it_msr)->vectorCount1);
	UINT32 covariance_count;

	for (cluster_bsl = 0; cluster_bsl < baseline_count; ++cluster_bsl)	// number of baselines/points
	{
		UpdateIgnoredMeasurements_GX(_it_msr, storeOriginalMeasurement);

		covariance_count = (*_it_msr)->vectorCount2;

		// skip covariances until next point
		(*_it_msr) += covariance_count * 3;

		if (covariance_count > 0)
			(*_it_msr)++;
	}
}
	

void dna_adjust::UpdateIgnoredMeasurements_Y(pit_vmsr_t _it_msr, bool storeOriginalMeasurement)
{
	// Use v_blockStationsMapUnique_ to get the correct block for each 
	// station in the measurement
	_it_u32u32_uint32_pair _it_bsmu;

	// Get estimated station coordinates matrix and index of the 
	// station within the matrix for station 1.  Since Y clusters are
	// correlated measurements which are never split, all stations will
	// be in the same block as station 1.
	_it_bsmu = v_blockStationsMapUnique_.begin() + (*_it_msr)->station1;
	matrix_2d* estimatedStations(&v_estimatedStations_.at(_it_bsmu->second));

	UINT32 stn1;
	UINT32 cluster_pnt, point_count((*_it_msr)->vectorCount1);
	UINT32 covariance_count;
	it_vstn_t stn1_it;
	double latitude, longitude, height, x, y, z;

	_COORD_TYPE_ coordType(CDnaStation::GetCoordTypeC((*_it_msr)->coordType));

	for (cluster_pnt = 0; cluster_pnt < point_count; ++cluster_pnt)
	{
		covariance_count = (*_it_msr)->vectorCount2;

		stn1_it = bstBinaryRecords_.begin() + (*_it_msr)->station1;

		stn1 = GetBlkMatrixElemStn1(_it_bsmu->second, _it_msr);

		// Get latest cartesian coordinates
		x = estimatedStations->get(stn1, 0);
		y = estimatedStations->get(stn1 + 1, 0);
		z = estimatedStations->get(stn1 + 2, 0);
		
		// initialise measurement (on the first adjustment only!)
		if (storeOriginalMeasurement)
			(*_it_msr)->preAdjMeas = (*_it_msr)->term1;

		// Convert to geographic coordinates?
        if (coordType == LLH_type_i || coordType == LLh_type_i)
		{
			CartToGeo<double>(x, y, z, &latitude, &longitude, &height, datum_.GetEllipsoidRef());
			(*_it_msr)->measAdj = latitude;
		}
		else
			(*_it_msr)->measAdj = x;
			
		// move to Y element
		(*_it_msr)++;
		// initialise measurement (on the first adjustment only!)
		if (storeOriginalMeasurement)
			(*_it_msr)->preAdjMeas = (*_it_msr)->term1;
		
		if (coordType == LLH_type_i || coordType == LLh_type_i)
			(*_it_msr)->measAdj = longitude;
		else
			(*_it_msr)->measAdj = y;

		// move to Z element
		(*_it_msr)++;
		// initialise measurement (on the first adjustment only!)
		if (storeOriginalMeasurement)
			(*_it_msr)->preAdjMeas = (*_it_msr)->term1;
			
		if (coordType == LLH_type_i || coordType == LLh_type_i)
		{
			// Reduce to ellipsoid height?
			if (fabs(stn1_it->geoidSep) > PRECISION_1E4)
			{
				(*_it_msr)->preAdjCorr = stn1_it->geoidSep;
				(*_it_msr)->term1 += (*_it_msr)->preAdjCorr;
			}

			(*_it_msr)->measAdj = height;
		}
		else
			(*_it_msr)->measAdj = z;

		covariance_count = (*_it_msr)->vectorCount2;

		// skip covariances until next point
		(*_it_msr) += covariance_count * 3;

		if (covariance_count > 0)
			(*_it_msr)++;
	}
}
	

// Vertical angle
//
void dna_adjust::UpdateIgnoredMeasurements_Z(pit_vmsr_t _it_msr, bool storeOriginalMeasurement)
{
	// initialise measurement (on the first adjustment only!)
	if (storeOriginalMeasurement)
		(*_it_msr)->preAdjMeas = (*_it_msr)->term1;

	// Use v_blockStationsMapUnique_ to get the correct block for each 
	// station in the measurement
	_it_u32u32_uint32_pair _it_bsmu;

	// Get estimated station coordinates matrix and index of the 
	// station within the matrix for station 1
	_it_bsmu = v_blockStationsMapUnique_.begin() + (*_it_msr)->station1;
	matrix_2d* estimatedStations_stn1(&v_estimatedStations_.at(_it_bsmu->second));
	UINT32 stn1(GetBlkMatrixElemStn1(_it_bsmu->second, _it_msr));
	
	// Get estimated station coordinates matrix and index of the 
	// station within the matrix for station 2
	_it_bsmu = v_blockStationsMapUnique_.begin() + (*_it_msr)->station2;
	matrix_2d* estimatedStations_stn2(&v_estimatedStations_.at(_it_bsmu->second));
	UINT32 stn2(GetBlkMatrixElemStn2(_it_bsmu->second, _it_msr));

	it_vstn_t_const stn1_it(bstBinaryRecords_.begin() + (*_it_msr)->station1);
	it_vstn_t_const stn2_it(bstBinaryRecords_.begin() + (*_it_msr)->station2);

	double local_12e, local_12n, local_12up;

	// compute vertical angle from estimated coordinates
	(*_it_msr)->measAdj = (VerticalAngle(
		estimatedStations_stn1->get(stn1, 0),					// X1
		estimatedStations_stn1->get(stn1 + 1, 0),					// Y1
		estimatedStations_stn1->get(stn1 + 2, 0),					// Z1
		estimatedStations_stn2->get(stn2, 0), 					// X2
		estimatedStations_stn2->get(stn2 + 1, 0),					// Y2
		estimatedStations_stn2->get(stn2 + 2, 0),					// Z2
		stn1_it->currentLatitude,
		stn1_it->currentLongitude,
		stn2_it->currentLatitude,
		stn2_it->currentLongitude,
		(*_it_msr)->term3,									// instrument height
		(*_it_msr)->term4,									// target height
		&local_12e,											// local_12e, ..12n, ..12up represent
		&local_12n,											// the geometric difference between
		&local_12up));										// station1 and station2

	// deflections available?
	if (fabs(stn1_it->verticalDef) > E4_SEC_DEFLECTION || fabs(stn1_it->meridianDef) > E4_SEC_DEFLECTION)
	{
		////////////////////////////////////////////////////////////////////////////
		// Correct for deflection of the vertical
		// 1. compute bearing from estimated coordinates
		double azimuth(Direction(
			estimatedStations_stn1->get(stn1, 0),		// X1
			estimatedStations_stn1->get(stn1 + 1, 0),		// Y1
			estimatedStations_stn1->get(stn1 + 2, 0),		// Z1
			estimatedStations_stn2->get(stn2, 0), 		// X2
			estimatedStations_stn2->get(stn2 + 1, 0),		// Y2
			estimatedStations_stn2->get(stn2 + 2, 0),		// Z2
			stn1_it->currentLatitude,
			stn1_it->currentLongitude));

		// 2. Compute correction
		(*_it_msr)->preAdjCorr = ZenithDeflectionCorrection<double>(		// Correction to vertical angle for deflection of vertical
			azimuth,														// geodetic azimuth
			stn1_it->verticalDef,											// deflection in prime vertical
			stn1_it->meridianDef);											// deflection in prime meridian
		////////////////////////////////////////////////////////////////////////////
	}
	else
		(*_it_msr)->preAdjCorr = 0.;

	// apply deflection correction
	(*_it_msr)->measAdj += (*_it_msr)->preAdjCorr;						
	// compute adjustment correction
	(*_it_msr)->measCorr = (*_it_msr)->measAdj - (*_it_msr)->preAdjMeas;
}
	

void dna_adjust::UpdateIgnoredMeasurements(pit_vmsr_t _it_msr, bool storeOriginalMeasurement)
{
	std::stringstream ss;

	// Since the stations connected by an ignored measurement may appear
	// in different blocks, it is essential to get the correct block 
	// number for each station.  For this, use v_blockStationsMapUnique_ which
	// is sorted on station id.
	std::sort(v_blockStationsMapUnique_.begin(), v_blockStationsMapUnique_.end());

	switch ((*_it_msr)->measType)
	{
	case 'A':	// Horizontal angle
		UpdateIgnoredMeasurements_A(_it_msr, storeOriginalMeasurement);
		break;
	case 'B':	// Geodetic azimuth
		UpdateIgnoredMeasurements_B(_it_msr, storeOriginalMeasurement);
		break;
	case 'C':	// Chord dist
		UpdateIgnoredMeasurements_C(_it_msr, storeOriginalMeasurement);
		break;
	case 'D':	// Direction set
		UpdateIgnoredMeasurements_D(_it_msr, storeOriginalMeasurement);
		break;
	case 'E':	// Ellipsoid arc
		UpdateIgnoredMeasurements_E(_it_msr, storeOriginalMeasurement);
		break;
	case 'G':	// GPS Baseline
		UpdateIgnoredMeasurements_G(_it_msr, storeOriginalMeasurement);
		break;
	case 'H':	// Orthometric height
		UpdateIgnoredMeasurements_H(_it_msr, storeOriginalMeasurement);
		break;
	case 'I':	// Astronomic latitude
		UpdateIgnoredMeasurements_I(_it_msr, storeOriginalMeasurement);
		break;
	case 'J':	// Astronomic longitude
		UpdateIgnoredMeasurements_J(_it_msr, storeOriginalMeasurement);
		break;
	case 'K':	// Astronomic azimuth
		UpdateIgnoredMeasurements_K(_it_msr, storeOriginalMeasurement);
		break;
	case 'L':	// Level difference
		UpdateIgnoredMeasurements_L(_it_msr, storeOriginalMeasurement);
		break;
	case 'M':	// MSL arc
		UpdateIgnoredMeasurements_M(_it_msr, storeOriginalMeasurement);
		break;
	case 'P':	// Geodetic latitude
		UpdateIgnoredMeasurements_P(_it_msr, storeOriginalMeasurement);
		break;
	case 'Q':	// Geodetic longitude
		UpdateIgnoredMeasurements_Q(_it_msr, storeOriginalMeasurement);
		break;
	case 'R':	// Ellipsoidal height
		UpdateIgnoredMeasurements_R(_it_msr, storeOriginalMeasurement);
		break;
	case 'S':	// Slope distance
		UpdateIgnoredMeasurements_S(_it_msr, storeOriginalMeasurement);
		break;
	case 'V':	// Zenith distance
		UpdateIgnoredMeasurements_V(_it_msr, storeOriginalMeasurement);
		break;
	case 'X':	// GPS Baseline cluster
		UpdateIgnoredMeasurements_X(_it_msr, storeOriginalMeasurement);
		break;
	case 'Y':	// GPS Point cluster
		UpdateIgnoredMeasurements_Y(_it_msr, storeOriginalMeasurement);
		break;
	case 'Z':	// Vertical angle
		UpdateIgnoredMeasurements_Z(_it_msr, storeOriginalMeasurement);
		break;
	default:
		ss << "UpdateIgnoredMeasurements(): Unknown measurement type - '" <<
			(*_it_msr)->measType << "'." << std::endl;
		SignalExceptionAdjustment(ss.str(), 0);
	}
}
	
void dna_adjust::ReduceYLLHMeasurementsforPrinting(vmsr_t& y_msr, matrix_2d& mpositions, printMeasurementsMode print_mode)
{
	it_vmsr_t _it_y_msr(y_msr.begin());
	UINT32 covr, cluster_msr, cluster_count(_it_y_msr->vectorCount1), covariance_count;
	it_vstn_t stn1_it;
	double latitude, longitude, height, x, y, z;

	for (cluster_msr=0; cluster_msr<cluster_count; ++cluster_msr)
	{
		covr = cluster_msr * 3;
		covariance_count = _it_y_msr->vectorCount2;

		stn1_it = bstBinaryRecords_.begin() + _it_y_msr->station1;

		// transform computed values to geographic
		switch (print_mode)
		{
		case computedMsrs:
			latitude = stn1_it->currentLatitude;
            longitude = stn1_it->currentLongitude;
            height = stn1_it->currentHeight;
			break;
		case adjustedMsrs:
		default:
			x = _it_y_msr->measAdj;
			_it_y_msr++;
			y = _it_y_msr->measAdj;
			_it_y_msr++;
			z = _it_y_msr->measAdj;

			// Convert to geographic
            CartToGeo<double>(x, y, z, &latitude, &longitude, &height, datum_.GetEllipsoidRef());
            _it_y_msr -= 2;
			break;
		}

		// Reduce ellipsoidal height to orthometric height
		if (fabs(stn1_it->geoidSep) > PRECISION_1E4)
			height -= stn1_it->geoidSep;		

		// Assign computed values		
		_it_y_msr->measAdj = latitude;
		_it_y_msr->measCorr = latitude - _it_y_msr->preAdjMeas;
		mpositions.put(covr, 0, latitude);
		_it_y_msr++;
		_it_y_msr->measAdj = longitude;
		_it_y_msr->measCorr = longitude - _it_y_msr->preAdjMeas;
		mpositions.put(covr+1, 0, longitude);
		_it_y_msr++;
		_it_y_msr->measAdj = height;
		_it_y_msr->measCorr = height - _it_y_msr->preAdjMeas;
		mpositions.put(covr+2, 0, height);

		// skip covariances until next point
		_it_y_msr += covariance_count * 3;

		if (covariance_count > 0)
			_it_y_msr++;
	}
}


void dna_adjust::SignalExceptionAdjustment(const std::string& msg, const UINT32 block_no)
{
	adjustStatus_ = ADJUST_EXCEPTION_RAISED;

	isPreparing_ = false;
	isCombining_ = false;
	isAdjusting_ = false;

	std::string error_msg(msg);

	switch (projectSettings_.a.adjust_mode)
	{
	case Phased_Block_1Mode:
	case PhasedMode:
		std::stringstream ss;
		ss << msg << std::endl << "  Phased adjustment terminated whilst processing block " << block_no + 1 << std::endl;
		error_msg = ss.str();
	}
	
	throw std::runtime_error(error_msg);
}

void dna_adjust::SetDefaultReferenceFrame()
{
	// Capture the datum recorded within binary file meta.
	// Note: One of two applications will have been the last to update the
	//		 reference frame in the bst/bms_meta and in every station record:
	//		 1. dnaimport
	//		 2. dnareftran
	//
	//	  *	When dnaimport runs, the default reference frame and its epoch is 
	//		set within within the bst/bms_meta, which may or may not be set by 
	//		the user, and is intended to be used to define the datum (ellipsoid 
	//		and epochal coordinate set) to be used throughout all computations.
	//	  *	For input files which do not specify datum, the datum in bst/bms_meta
	//		will be the same as datum defined for all stations and measurement
	//		records.  For input files which do specify a datum for each station
	//		and measurement, the datum in the bst/bms_meta may be different to what
	//		is set for each station and measurement.  In either case, dnaadjust
	//		takes what is set in the bst/bms_meta.  The user is responsible for
	//		knowing the data and its reference frame.
	//	  *	If dnareftran has been run, every station and measurement record will
	//		have been transformed to the one frame and epoch. These reference frame
	//		records will be consistent with the bst/bms_meta.
	//    * If the user has not specified an epoch with WGS84, then the epoch will
	//      be set to "time immemorial" being 01.01.1900. Any occurrence of this will
	//      cause output to be set to an empty string ("").
	try {
		// Load binary stations data.  Throws runtime_error on failure.
		BstFile bst;
		bst.LoadFileMeta(projectSettings_.a.bst_file, bst_meta_);
		datum_.SetDatumFromEpsg(bst_meta_.epsgCode, bst_meta_.epoch);
	}
	catch (const std::runtime_error& e) {
		SignalExceptionAdjustment(e.what(), 0);
	}
}

void dna_adjust::LoadNetworkFiles()
{
    adj_file << "+ Loading network files" << std::endl;
    
    try {
        NetworkDataLoader loader(projectSettings_);
        
        UINT32 measurementVarianceCount = 0;
        
        bool success = false;
        if (projectSettings_.a.adjust_mode == SimultaneousMode) {
            success = loader.LoadForSimultaneous(
                &bstBinaryRecords_,
                bst_meta_,
                &vAssocStnList_,
                &bmsBinaryRecords_,
                bms_meta_,
                v_ISL_,
                &v_blockStationsMap_,
                &v_CML_,
                bstn_count_,
                asl_count_,
                bmsr_count_,
                unknownParams_,
                unknownsCount_,
                measurementParams_,
                measurementCount_,
                measurementVarianceCount,
                &blockCount_,
                &v_JSL_,
                &v_unknownsCount_,
                &v_measurementCount_,
                &v_measurementVarianceCount_,
                &v_measurementParams_,
                &v_ContiguousNetList_,
                &v_blockMeta_,
                &v_parameterStationList_,
                &v_paramStnAppearance_,
                &v_junctionVariances_,
                &v_junctionVariancesFwd_);
        } else {
            // For phased mode, we need to initialize blockCount_ to 1 initially
            // The actual block count will be loaded from the segmentation file later
            blockCount_ = 1;
            success = loader.LoadForPhased(
                &bstBinaryRecords_,
                bst_meta_,
                &vAssocStnList_,
                &bmsBinaryRecords_,
                bms_meta_,
                v_ISL_,
                &v_blockStationsMap_,
                &v_CML_,
                bstn_count_,
                asl_count_,
                bmsr_count_,
                unknownParams_,
                unknownsCount_,
                measurementParams_,
                measurementCount_,
                measurementVarianceCount);
        }
           
        if (!success) {
            SignalExceptionAdjustment("LoadNetworkFiles(): Failed to load network files", 0);
        }
        
        // Check if all measurements are ignored (no measurements to process)
        if (measurementCount_ == 0 || (v_CML_.size() > 0 && v_CML_.at(0).empty())) {
            SignalExceptionAdjustment("LoadNetworkFiles(): No valid measurements to process. All measurements may be ignored.", 0);
        }
        
        // Ensure v_blockStationsMap_ matches expected state
        if (v_blockStationsMap_.empty()) {
            v_blockStationsMap_.resize(1);
        }
    }
    catch (const std::exception& e) {
        SignalExceptionAdjustment(std::string("LoadNetworkFiles(): ") + e.what(), 0);
    }
}


void dna_adjust::AddDiscontinuitySites(vstring& constraintStns)
{
	// sort on original name
	std::sort(bstBinaryRecords_.begin(), bstBinaryRecords_.end(), CompareStnOriginalName<station_t, std::string>());
	
	std::pair<it_vstn_t, it_vstn_t> it_stn_range;
	std::string constraint;

	_it_vstr const_it;

	vstring constraintDiscontStns;

	// for every constraint (note every second value will be the constraint code)
	for (const_it=constraintStns.begin(); const_it!=constraintStns.end(); ++const_it) 
	{
		// find this constraint station in the station map
		it_stn_range = equal_range(bstBinaryRecords_.begin(), bstBinaryRecords_.end(), (*const_it), CompareStnOriginalName<station_t, std::string>());
		if (it_stn_range.first == it_stn_range.second)
		{
			// Not in the map.  Get the next one
			continue;
		}

		// this constraint station is one of the discontinuity sites.
		// Add the "discontinuity name"
		constraintDiscontStns.push_back(it_stn_range.first->stationName);
		// Add the constraint
		if (++const_it == constraintStns.end())
			break;
		constraintDiscontStns.push_back(*const_it);
	}

	// add the constraint discontinuity sites
	constraintStns.insert(constraintStns.end(), constraintDiscontStns.begin(), constraintDiscontStns.end());

	// restore original sort order
	std::sort(bstBinaryRecords_.begin(), bstBinaryRecords_.end(), CompareStnName<station_t>());
}

// Take user-supplied type B uncertainties and compute relevant elements for
// application to the output results
void dna_adjust::InitialiseTypeBUncertainties()
{
	if (!projectSettings_.o._apply_type_b_global && !projectSettings_.o._apply_type_b_file)
		return;

	dna_io_tbu tbu;

	// Apply global type B uncertainties to all sites
	if (projectSettings_.o._apply_type_b_global)
	{
		// Load the typeb argument. The typeb argument contains type b uncertainties to be 
		// applied (by default) to all sies
		try {
			
			tbu.load_tbu_argument(projectSettings_.a.type_b_global, typeBUncertaintyGlobal_);
		}
		catch (const std::runtime_error& e) {
			SignalExceptionAdjustment(e.what(), 0);
		}
		catch (...) {
			return;
		}		
	}

	// Apply (or overwrite) type B uncertainties to specific sites
	if (projectSettings_.o._apply_type_b_file)
	{
		// load station map
		v_string_uint32_pair vStnsMap;
		if (projectSettings_.a.map_file.empty())
			projectSettings_.a.map_file = projectSettings_.g.input_folder + FOLDER_SLASH + projectSettings_.g.network_name + ".map";
		LoadStationMap(&vStnsMap, projectSettings_.a.map_file);

		// Load the typeb file. The typeb file contains site-specific type b uncertainties
		try {

			tbu.load_tbu_file(projectSettings_.a.type_b_file, v_typeBUncertaintiesLocal_, vStnsMap);
			std::sort(v_typeBUncertaintiesLocal_.begin(), v_typeBUncertaintiesLocal_.end());

			// reduce to the cartesian reference frame
			tbu.reduce_uncertainties_local(v_typeBUncertaintiesLocal_, bstBinaryRecords_);
		}
		catch (const std::ifstream::failure& f) {
			SignalExceptionAdjustment(f.what(), 0);
		}
		catch (const std::runtime_error& e) {
			SignalExceptionAdjustment(e.what(), 0);
		}
	}

	// Now, set the method of application
	UINT32 stn(0);
	it_pair_type_b_uncertainty it_tbu;

	// Create a vector of elements (with a size equal to the station count) that manages how
	// type b uncertainties are to be applied.
	v_typeBUncertaintyMethod_.resize(bstBinaryRecords_.size());
	// Create a map between station ID and the type b uncertainties vector
	v_stationTypeBMap_.resize(bstBinaryRecords_.size());

	for (stn = 0; stn < bstBinaryRecords_.size(); ++stn)
	{
		v_typeBUncertaintyMethod_.at(stn).station_id = stn;
		v_stationTypeBMap_.at(stn).first = stn;
		v_stationTypeBMap_.at(stn).second = 0;

		// set method to be global by default
		if (projectSettings_.o._apply_type_b_global)
		{
			v_typeBUncertaintyMethod_.at(stn).apply = true;
			v_typeBUncertaintyMethod_.at(stn).method = type_b_global;
		}

		// overwrite global with site specific type b
		if (projectSettings_.o._apply_type_b_file)
		{
			it_tbu = equal_range(v_typeBUncertaintiesLocal_.begin(), v_typeBUncertaintiesLocal_.end(), 
				stn, CompareTypeBStationID<UINT32, type_b_uncertainty>());
			if (it_tbu.first != it_tbu.second)
			{
				v_typeBUncertaintyMethod_.at(stn).apply = true;
				v_typeBUncertaintyMethod_.at(stn).method = type_b_local;

				// map the local uncertainty index
				v_stationTypeBMap_.at(stn).second = static_cast<UINT32>(std::distance(v_typeBUncertaintiesLocal_.begin(), it_tbu.first));
			}
		}

	}
}

	

void dna_adjust::LoadStationMap(pv_string_uint32_pair stnsMap, const std::string& stnmap_file)
{
	try {
		// Load station map.  Throws runtime_error on failure.
		dynadjust::iostreams::MapFile map;
		map.LoadFile(stnmap_file, stnsMap);
	}
	catch (const std::runtime_error& e) {
		SignalExceptionAdjustment(e.what(), 0);
	}
}
	

void dna_adjust::ResizeMatrixVectors()
{
	adj_file << "+ Allocating memory" << std::endl;
	
	switch (projectSettings_.a.adjust_mode)
	{
	case SimultaneousMode:		
		v_parameterStationCount_.resize(blockCount_);
		v_parameterStationList_.resize(blockCount_);
		v_parameterStationList_.at(0) = v_ISL_.at(0);
		break;
	case PhasedMode:
	case Phased_Block_1Mode:
		// Load Segmentation data.  LoadPhasedBlocks() calls LoadSegmentationFile() which
		// and assigns blockCount_
		LoadPhasedBlocks();

		v_msrTally_.resize(blockCount_);

		v_statSummary_.resize(blockCount_);
		v_pseudoMeasCountFwd_.resize(blockCount_);
		v_sigmaZero_.resize(blockCount_);
		v_chiSquaredUpperLimit_.resize(blockCount_);
		v_chiSquaredLowerLimit_.resize(blockCount_);
		v_passFail_.resize(blockCount_);

		v_correctionsR_.resize(blockCount_);
		break;
	default:
		// thrown in LoadNetworkFiles()
		break;
	}

	v_msrTally_.resize(blockCount_);

	v_design_.resize(blockCount_);
	v_measMinusComp_.resize(blockCount_);
	v_AtVinv_.resize(blockCount_);
	v_normals_.resize(blockCount_);
	v_normalsR_.resize(blockCount_);		// used to back up normals prior to adding parameter stations and junctions
	v_precAdjMsrsFull_.resize(blockCount_);

	v_rigorousVariances_.resize(blockCount_);
	v_estimatedStations_.resize(blockCount_);
	v_originalStations_.resize(blockCount_);
	v_corrections_.resize(blockCount_);

	if (projectSettings_.a.multi_thread)
		v_normalsRC_.resize(blockCount_);

	for (UINT32 block(0); block<blockCount_; ++block)
	{
		// sparse matrix
		v_design_.at(block).matrixType(mtx_sparse);

		// lower triangular
		v_normals_.at(block).matrixType(mtx_lower);
		v_normalsR_.at(block).matrixType(mtx_lower);
		v_rigorousVariances_.at(block).matrixType(mtx_lower);

		if (projectSettings_.a.adjust_mode == PhasedMode)
		{
			v_junctionVariances_.at(block).matrixType(mtx_lower);
			v_junctionVariancesFwd_.at(block).matrixType(mtx_lower);
		}
	}
	
}
	
void dna_adjust::LoadPhasedBlocks()
{
	LoadSegmentationFile();

	if (v_ISL_.size() != v_JSL_.size())
	{
		std::stringstream ss;
		ss << "LoadPhasedBlocks(): An unrecoverable error was encountered when loading the phased adjustment blocks." << std::endl;
		SignalExceptionAdjustment(ss.str(), 0);
	}
	if (v_JSL_.size() != v_CML_.size())
	{
		std::stringstream ss;
		ss << "LoadPhasedBlocks(): An unrecoverable error was encountered when loading the phased adjustment blocks." << std::endl;
		SignalExceptionAdjustment(ss.str(), 0);
	}
	
}


void dna_adjust::LoadSegmentationMetrics()
{
	v_blockMeta_.resize(blockCount_);
	v_blockStationsMap_.resize(blockCount_);
	v_parameterStationList_.resize(blockCount_);
	v_paramStnAppearance_.resize(blockCount_);
	v_junctionVariances_.resize(blockCount_);
	v_junctionVariancesFwd_.resize(blockCount_);
	v_junctionEstimatesFwd_.resize(blockCount_);
	v_junctionEstimatesRev_.resize(blockCount_);
	v_rigorousStations_.resize(blockCount_);

	v_measurementVarianceCount_ = v_measurementCount_;
	
	vUINT32 parameterStations;
	UINT32 b, netID(999999), stn, c;

	try {
		// read block data
		for (b=0; b<blockCount_; ++b)
		{
			// set Meta
			// If this block has a different ID, then it must be the first.  
			if (netID != v_ContiguousNetList_.at(b))
				v_blockMeta_.at(b)._blockFirst = true;
			else
				v_blockMeta_.at(b)._blockFirst = false;

			netID = v_ContiguousNetList_.at(b);

			// if this block's ID is different to the next, then it must be the last
			if (b < blockCount_ - 1)
			{
				// the reason why the code in this for loop can't be placed in the preceding loop
				// is because element b+1 is required
				if (v_ContiguousNetList_.at(b) != v_ContiguousNetList_.at(b+1))		
					v_blockMeta_.at(b)._blockLast = true;
				else
					v_blockMeta_.at(b)._blockLast = false;
			}
			else // if (b == blockCount_ - 1)		// final block in whole segmentation list?
				v_blockMeta_.at(b)._blockLast = true;

			// now, test whether this block is a single block
			if (v_blockMeta_.at(b)._blockFirst && v_blockMeta_.at(b)._blockLast)
				v_blockMeta_.at(b)._blockIsolated = true;
			// or, test whether this block is an intermediate block
			else if (!v_blockMeta_.at(b)._blockFirst && !v_blockMeta_.at(b)._blockLast)
				v_blockMeta_.at(b)._blockIntermediate = true;

			// add station map to enable ISL/JSL stations to be referenced according to block
			v_parameterStationList_.at(b) = v_ISL_.at(b);
			v_parameterStationList_.at(b).insert(v_parameterStationList_.at(b).end(),
				v_JSL_.at(b).begin(), v_JSL_.at(b).end());
			std::sort(v_parameterStationList_.at(b).begin(), v_parameterStationList_.at(b).end());

			v_paramStnAppearance_.at(b).resize(v_parameterStationList_.at(b).size());

			// Add all stations to parameterStations
			parameterStations.insert(parameterStations.end(), v_parameterStationList_.at(b).begin(), v_parameterStationList_.at(b).end());

			stn = 0;

			// check memory availability for block station map
			if (v_blockStationsMap_.at(b).max_size() <= v_parameterStationCount_.at(b))
			{
				std::stringstream ss;
				ss << "LoadSegmentationMetrics(): Could not allocate sufficient memory for blockStationsMap" << std::endl;
				SignalExceptionAdjustment(ss.str(), 0);
			}

			// fill block station map
			for (c=0; c<v_parameterStationCount_.at(b); ++c)
				v_blockStationsMap_.at(b)[v_parameterStationList_.at(b).at(c)] = stn++;	
		}

	}
	catch (...) {
		std::stringstream ss;
		ss << "LoadSegmentationMetrics(): An error was encountered when attempting to calculate " << std::endl <<
			"  the segmentation metrics" << std::endl;
		SignalExceptionAdjustment(ss.str(), 0);
	}

	v_unknownParams_ = v_unknownsCount_;
	v_measurementParams_ = v_measurementCount_;	

	// update v_measurementVarianceCount_ if adjusted GNSS measurement
	// variances are to be propagated into another system. Initially, 
	// v_measurementVarianceCount_ is copied from v_measurementCount_ 
	// which contains a count of the rows in the design matrix.  This
	// count represents standard deviations for non-GNSS measurements
	// and the diagonals of GNSS measurement variances. The
	// following code increments the count by 3 for GNSS measurements 
	// so that the count represents an upper triangular matrix,
	// being 6 (xx, xy, xz, yy, yz, zz).
	
	it_vUINT32 _it_block_msr;
	it_vmsr_t _it_msr;

	for (b=0; b<blockCount_; ++b)
	{
		for (_it_block_msr=v_CML_.at(b).begin(); 
			_it_block_msr!=v_CML_.at(b).end(); 
			++_it_block_msr)
		{
			_it_msr = bmsBinaryRecords_.begin() + (*_it_block_msr);

			// Is this measurement ignored?
			if (_it_msr->ignore)
				continue;

			// When a GNSS measurement is found, add three 
			// elements for X,Y,Z covariance
			switch (_it_msr->measType)
			{
			case 'G':
			case 'X':
			case 'Y':
				
				// Increment by 3 to take count of diagonals
				// to 6, representing upper triangular matrix
				v_measurementVarianceCount_.at(b) += _it_msr->vectorCount1 * 3;

				//v_measurementVarianceCount_.at(b) += 3;
				break;
			default:
				continue;
			}
		}
	}


	// Remove duplicates from the parameterStations list and
	// compute the total unknowns for the entire network	
	strip_duplicates(parameterStations);

	// Update "unknowns' count for whole of adjustment
	for_each(parameterStations.begin(), parameterStations.end(),
		[this](const UINT32& stn) {

			if (bstBinaryRecords_.at(stn).stationConst[0] == 'C')
				unknownParams_--;
		
			if (bstBinaryRecords_.at(stn).stationConst[1] == 'C')
				unknownParams_--;
		
			if (bstBinaryRecords_.at(stn).stationConst[2] == 'C')
				unknownParams_--;		
	});

	if (unknownParams_ == 0 && unknownsCount_ > 0)
		allStationsFixed_ = true;

	try {
		// Create station appearance list.  Throws runtime_error on failure.
		SegFile seg;
		seg.CreateStnAppearanceList(v_paramStnAppearance_, 
			v_parameterStationList_, vAssocStnList_);
	}
	catch (const std::runtime_error& e) {
		SignalExceptionAdjustment(e.what(), 0);
	}

	if (projectSettings_.g.verbose > 1)
	{
		debug_file << std::endl << "Station appearance list" << std::endl;

		it_vUINT32_const _it_const;
		it_vstn_appear _it_appear;

		for (UINT32 block(0); block!=blockCount_; ++block)
		{
			debug_file << "Block " << block << std::endl;
			debug_file << std::setw(HEADER_20) << "Station" << std::setw(HEADER_20) << "Forward" << std::setw(HEADER_20) << "Reverse" << std::setw(HEADER_20) << " " << std::endl;

			for (_it_const=v_parameterStationList_.at(block).begin(),
				_it_appear=v_paramStnAppearance_.at(block).begin();
				_it_appear!=v_paramStnAppearance_.at(block).end();
				++_it_const, ++_it_appear)
			{
				// station id
				debug_file << std::setw(HEADER_20) << *_it_const;
				// forward
				if (_it_appear->first_appearance_fwd)
					debug_file << std::setw(HEADER_20) << "true";
				else
					debug_file << std::setw(HEADER_20) << "false";
				// reverse
				if (_it_appear->first_appearance_rev)
					debug_file << std::setw(HEADER_20) << "true" << std::endl;
				else
					debug_file << std::setw(HEADER_20) << "false" << std::endl;
			}
			debug_file << std::endl;
		}

		debug_file << std::endl;
	}

}

void dna_adjust::LoadSegmentationFileParameters(const std::string& seg_filename)
{
	UINT32 blockThreshold, minInnerStns;

	try {
		// Load segmentation file.  Throws runtime_error on failure.
		SegFile seg;
		seg.LoadSegFileHeaderF(seg_filename, 
			blockCount_, blockThreshold, minInnerStns);
	}
	catch (const std::runtime_error& e) {
		SignalExceptionAdjustment(e.what(), 0);
	}
}


void dna_adjust::LoadSegmentationFile()
{
	UINT32 blockThreshold, minInnerStns;

	try {
		// Load segmentation file.  Throws runtime_error on failure.
		SegFile seg;
		seg.LoadSegFile(projectSettings_.a.seg_file, 
			blockCount_, blockThreshold, minInnerStns,
			v_ISL_, v_JSL_, v_CML_,
			true, &bmsBinaryRecords_,
			&v_measurementCount_, &v_unknownsCount_, &v_ContiguousNetList_,
			&v_parameterStationCount_);
	}
	catch (const std::runtime_error& e) {
		SignalExceptionAdjustment(e.what(), 0);
	}

	// load segmentation metrics into appropriate vectors
	LoadSegmentationMetrics();
}
	

void dna_adjust::RemoveDuplicateStations(vUINT32& vStations)
{
	if (vStations.size() < 2)
		return;
	
	strip_duplicates(vStations);
}
	

// Used only for simultaneous adjustment
// Required since the binary stations list may contain stations with no measurements
// Not required for phased adjustment since these stations are removed during the
// segmentation algorithm.
	

void dna_adjust::SortMeasurementsbyType(v_uint32_u32u32_pair& msr_block)
{
	if (msr_block.size() < 2)
		return;
	CompareMeasType_PairFirst<measurement_t, UINT32> meastypeCompareFunc(&bmsBinaryRecords_);
	std::sort(msr_block.begin(), msr_block.end(), meastypeCompareFunc);
}


void dna_adjust::SortMeasurementsbyFromStn(v_uint32_u32u32_pair& msr_block)
{
	if (msr_block.size() < 2)
		return;
	CompareMeasFromStn_PairFirst<measurement_t, UINT32> measfromstnCompareFunc(&bmsBinaryRecords_);
	std::sort(msr_block.begin(), msr_block.end(), measfromstnCompareFunc);
}


void dna_adjust::SortMeasurementsbyToStn(v_uint32_u32u32_pair& msr_block)
{
	if (msr_block.size() < 2)
		return;
	CompareMeasToStn_PairFirst<measurement_t, UINT32> meastostnCompareFunc(&bmsBinaryRecords_);
	std::sort(msr_block.begin(), msr_block.end(), meastostnCompareFunc);
}


void dna_adjust::SortMeasurementsbyValue(v_uint32_u32u32_pair& msr_block)
{
	if (msr_block.size() < 2)
		return;
	CompareMeasValue_PairFirst<measurement_t, UINT32> measvalueCompareFunc(&bmsBinaryRecords_);
	std::sort(msr_block.begin(), msr_block.end(), measvalueCompareFunc);
}


void dna_adjust::SortMeasurementsbyResidual(v_uint32_u32u32_pair& msr_block)
{
	if (msr_block.size() < 2)
		return;
	CompareMeasResidual_PairFirst<measurement_t, UINT32> measresidualCompareFunc(&bmsBinaryRecords_);
	std::sort(msr_block.begin(), msr_block.end(), measresidualCompareFunc);
}


void dna_adjust::SortMeasurementsbyAdjSD(v_uint32_u32u32_pair& msr_block)
{
	if (msr_block.size() < 2)
		return;
	CompareMeasAdjSD_PairFirst<measurement_t, UINT32> measadjsdCompareFunc(&bmsBinaryRecords_);
	std::sort(msr_block.begin(), msr_block.end(), measadjsdCompareFunc);
}


void dna_adjust::SortMeasurementsbyNstat(v_uint32_u32u32_pair& msr_block)
{
	if (msr_block.size() < 2)
		return;
	CompareMeasNstat_PairFirst<measurement_t, UINT32> measnstatCompareFunc(&bmsBinaryRecords_);
	std::sort(msr_block.begin(), msr_block.end(), measnstatCompareFunc);
}
	



}	// namespace networkadjust
}	// namespace dynadjust
