//============================================================================
// Name         : dnaadjust-multi.cpp
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
// Description  : DynAdjust Network Adjustment (multi-thread) library
//============================================================================

#include <dynadjust/dnaadjust/dnaadjust.hpp>

namespace dynadjust { namespace networkadjust {

concurrent_block_adjustment<UINT32> concurrentAdjustments;

extern concurrent_queue<UINT32> combineAdjustmentQueue;
extern concurrent_queue<UINT32> prepareAdjustmentQueue;
extern std::mutex dbg_file_mutex;
extern std::exception_ptr fwd_error;
extern std::exception_ptr rev_error;
extern std::exception_ptr cmb_error;
extern std::exception_ptr prep_error;

bool combineAdjustmentExceptionThrown(const std::vector<std::exception_ptr>& cmb_errors_)
{
	// got a forward or reverse exception, or 
	// have any of the other combination adjustments failed?
	if (fwd_error || rev_error)
	{
		if (combineAdjustmentQueue.is_empty())
			combineAdjustmentQueue.push_and_notify(99999999);
		return true;
	}

	UINT32 err_size(static_cast<UINT32>(cmb_errors_.size()));
	for (UINT32 err=0; err<err_size; ++err)
	{
		if (cmb_errors_.at(err))
		{
			if (combineAdjustmentQueue.is_empty())
				combineAdjustmentQueue.push_and_notify(99999999);
			return true;
		}
	}
	return false;
}

bool prepareAdjustmentExceptionThrown(const std::vector<std::exception_ptr>& prep_errors_)
{
	// have any of the other combination adjustments failed?
	UINT32 err_size(static_cast<UINT32>(prep_errors_.size()));
	for (UINT32 err=0; err<err_size; ++err)
	{
		if (prep_errors_.at(err))
		{
			if (prepareAdjustmentQueue.is_empty())
				prepareAdjustmentQueue.push_and_notify(99999999);
			return true;
		}
	}
	return false;
}

// Multi thread phased adjustment
// Notes for general understanding of the multi thread operation:
//	- Thread 1 undertakes a forward adjustment
//	- Thread 2 undertakes a reverse adjustment
//	- Thread 3 undertakes a combine adjustment
//	- Threads 1 and 2 are executed simultaneously without any due regard for each other.  The 
//	  coordinates and uncertainties produced from the forward and reverse runs are managed in
//	  separate matrices, hence the two threads may run without conflict.
//    After each block is adjusted, the corresponding forward and reverse mutexes for each
//	  block are unlocked.
//	- Immediately after block n has been adjusted in both forward and reverse directions, 
//	  Thread 3 undertakes combination adjustment of forward/reverse blocks.  This is handled
//	  by waiting for the two (forward and reverse) mutexes for each block.  Hence, it is
//	  only when both Threads 1 and 2 have triggered 'unlock' on the respective forward/reverse
//	  mutexes that Thread 3 commence.
//
void dna_adjust::AdjustPhasedMultiThread()
{
	initialiseIteration();

	std::string corr_msg;
	std::ostringstream ss;
	UINT32 i;
	bool iterate(true);

	std::chrono::milliseconds iteration_time(std::chrono::milliseconds(0));
	cpu_timer it_time, tot_time;

#if defined(__ICC) || defined(__INTEL_COMPILER)		// Intel compiler
	std::shared_ptr<std::thread> f, r, c;
	std::vector< std::shared_ptr<std::thread> > mt_adjust_threads_sp;
#else
	std::vector<std::thread> mt_adjust_threads;
#endif	

	concurrentAdjustments.resize_runs(blockCount_);

	// do until convergence criteria is met
	for (i=0; i<projectSettings_.a.max_iterations; ++i)
	{
		if (IsCancelled())
			break;

		SetcurrentBlock(0);
		blockLargeCorr_ = 0;
		largestCorr_ = 0.0;
		maxCorr_ = 0.0;

		///////////////////////////////////
		// Print the iteration # to adj file.
		// protected write to adj file (not needed here since write to
		// adj file at this stage is via single thread
		printer_->PrintIteration(incrementIteration());
		///////////////////////////////////

		concurrentAdjustments.reset_adjustment_runs();		
		concurrentAdjustments.combine_run_stopped();

		combineAdjustmentQueue.reset_blocks_coming();

		// Forward and reverse threads sequentially adjust all blocks in
		// forward and reverse directions.
#if defined(__ICC) || defined(__INTEL_COMPILER)		// Intel compiler
		f.reset(new thread(adjust_forward_thread(this, std::ref(fwd_error))));
		mt_adjust_threads_sp.push_back(f);
		r.reset(new thread(adjust_reverse_thread(this, std::ref(rev_error))));
		mt_adjust_threads_sp.push_back(r);
#else
		mt_adjust_threads.push_back(std::thread(adjust_forward_thread(this, std::ref(fwd_error))));
		mt_adjust_threads.push_back(std::thread(adjust_reverse_thread(this, std::ref(rev_error))));
#endif			
		
		// Combination thread fires up as many threads as there are cores
		// on the current PC to concurrently combine blocks
#if defined(__ICC) || defined(__INTEL_COMPILER)		// Intel compiler
		if (CombinationThreadRequired())
		{
			c.reset(new thread(adjust_combine_thread(this, std::ref(cmb_error))));
			mt_adjust_threads_sp.push_back(c);
		}
#else
		if (CombinationThreadRequired())
			mt_adjust_threads.push_back(std::thread(adjust_combine_thread(this, std::ref(cmb_error))));
#endif	

		// Start the clock
		it_time.start();

		// Start the forward, reverse and combine threads, which commences the
		// network adjustment
#if defined(__ICC) || defined(__INTEL_COMPILER)		// Intel compiler
		for_each(mt_adjust_threads_sp.begin(), mt_adjust_threads_sp.end(), std::mem_fn(&thread::join));
#else
		for_each(mt_adjust_threads.begin(), mt_adjust_threads.end(), std::mem_fn(&std::thread::join));
#endif
		// This point is reached when the threads have finished
		iteration_time = std::chrono::duration_cast<std::chrono::milliseconds>(it_time.elapsed().wall);

		//delete mt_adjust_threads;
#if defined(__ICC) || defined(__INTEL_COMPILER)		// Intel compiler
		mt_adjust_threads_sp.clear();
#else
		mt_adjust_threads.clear();
#endif
		
		concurrentAdjustments.combine_run_stopped();
				
		// Were any exceptions thrown?  If so, re-throw and let
		// test stub handle the exception
		if (fwd_error)
			// exception thrown in the forward pass
			std::rethrow_exception(fwd_error);
		else if (rev_error)
			// exception thrown in the reverse pass
			std::rethrow_exception(rev_error);
		else if (cmb_error)
			// exception thrown in one of the combination adjustments
			std::rethrow_exception(cmb_error);

		if (IsCancelled())
			break;

		ss.str("");
		if (iteration_time >= std::chrono::seconds(1)) {
			auto seconds = std::chrono::duration_cast<std::chrono::seconds>(iteration_time);
			ss << seconds.count() << "s";
		} else {
			ss << iteration_time.count() << "ms";
		}

		///////////////////////////////////
		// protected write to adj file (not needed here since write to
		// adj file at this stage is via single thread
		adj_file << std::setw(PRINT_VAR_PAD) << std::left << "Elapsed time" << ss.str() << std::endl;
		OutputLargestCorrection(corr_msg);
		UpdateAdaptiveRelaxation();
		///////////////////////////////////

		if (projectSettings_.g.verbose)
			debug_file << concurrentAdjustments.print_adjusted_blocks();

		iterationCorrections_.add_message(corr_msg);
		iterationQueue_.push_and_notify(CurrentIteration());				// currentIteration begins at 1, so not zero-indexed

		// continue iterating?
		iterate = fabs(maxCorr_) > projectSettings_.a.iteration_threshold;

		// Update normals and measured-computed matrices for the next iteration.
		// Similar to PrepareAdjustment, UpdateAdjustment prepares every block
		// in the network so that forward and reverse adjustments can commence
		// at the same time.
		UpdateAdjustment(iterate);

		if (!iterate)
			break;
	}

	isAdjusting_ = false;

	if (adjustStatus_ > ADJUST_TEST_FAILED)
		return;

	// This is similar to the check in ValidateandFinaliseAdjustment API for
	// phased single-thread mode.
	if (IsCancelled())
	{
		adjustStatus_ = ADJUST_CANCELLED;
		return;
	}

	if (blowupDetected_)
		adjustStatus_ = ADJUST_MAX_ITERATIONS_EXCEEDED;
	else if (CurrentIteration() == projectSettings_.a.max_iterations)
		adjustStatus_ = ADJUST_MAX_ITERATIONS_EXCEEDED;

	// Print status
	printer_->PrintAdjustmentStatus();
	// Compute and print time taken to run adjustment
	PrintAdjustmentTime(tot_time, total_time);
}


bool dna_adjust::CombinationThreadRequired()
{
	for (UINT32 block(0); block<blockCount_; ++block) 
	{
		// Combination adjustment is required if there
		// is at least one intermediate block.
		if (v_blockMeta_.at(block)._blockIntermediate)
			return true;
	}
	return false;
}


void dna_adjust::UpdateEstimatesFinalNoCombine()
{
	for (UINT32 block(0); block<blockCount_; ++block) 
		UpdateEstimatesFinal(block);
}


void dna_adjust::SolveMTTry(bool COMPUTE_INVERSE, const UINT32& block)
{
	// Least Squares Solution
	try {            
		SolveMT(COMPUTE_INVERSE, block);
	}
	catch (const std::runtime_error& e) {

		// could not invert matrix.  throw error

		if (projectSettings_.g.verbose)
		{
			debug_file << "Block " << block + 1 << " (reverse, multi-thread)" << std::endl;
			debug_file << "Pre-adjustment Estimates" << std::fixed << std::setprecision(16) << v_estimatedStationsR_.at(block);
			debug_file << "Design" << std::fixed << std::setprecision(16) << v_designR_.at(block);
			debug_file << "Measurements" << std::fixed << std::setprecision(16) << v_measMinusCompR_.at(block);
			debug_file << "At * V-inv" << std::fixed << std::setprecision(16) << v_AtVinvR_.at(block);
			debug_file << "Normals " << std::fixed << std::setprecision(16) << v_normalsR_.at(block) << std::endl;
			debug_file.flush();
		}

		SignalExceptionAdjustment(e.what(), block);
	}
}

void dna_adjust::SolveMT(bool COMPUTE_INVERSE, const UINT32& block)
{
	if (COMPUTE_INVERSE)
	{
		// Compute inverse of normals (aposteriori variance matrix)
		// (AT * V-1 * A)-1
		FormInverseVarianceMatrix(&(v_normalsR_.at(block)), false);
	}

	// compute weighted "measured minus computed"
	matrix_2d At_Vinv_m(v_designR_.at(block).columns(), 1);
	At_Vinv_m.multiply(v_AtVinvR_.at(block), "N", v_measMinusCompR_.at(block), "N");

	// Solve corrections from normal equations
	v_correctionsR_.at(block).redim(v_designR_.at(block).columns(), 1);
	v_correctionsR_.at(block).multiply(v_normalsR_.at(block), "N", At_Vinv_m, "N");

	// debug output?
	if (projectSettings_.g.verbose > 3)
	{
		dbg_file_mutex.lock();

		UINT32 i;
		switch (projectSettings_.a.adjust_mode)
		{
		case PhasedMode:
			for (i=0; i<v_parameterStationList_.at(block).size(); ++i)
			{
				debug_file << bstBinaryRecords_.at(v_parameterStationList_.at(block).at(i)).stationName;
				if (find(v_ISL_.at(block).begin(), v_ISL_.at(block).end(), v_parameterStationList_.at(block).at(i)) != v_ISL_.at(block).end())
					debug_file << "i ";
				else if (find(v_JSL_.at(block).begin(), v_JSL_.at(block).end(), v_parameterStationList_.at(block).at(i)) != v_JSL_.at(block).end())
						debug_file << "j ";
				else 
					debug_file << " ";				
			}
			break;
		}

		debug_file << std::endl;
		
		debug_file << "Block " << block + 1 << " (Reverse)";
		debug_file << std::endl;
		debug_file << "Pre-adjustment Estimates" << std::fixed << std::setprecision(16) << v_estimatedStationsR_.at(block);
		debug_file << "Block " << block + 1 << " (Reverse)";
		debug_file << std::endl;
		debug_file << "Design" << std::fixed << std::setprecision(16) << v_designR_.at(block);
		debug_file << "Block " << block + 1 << " (Reverse)";
		debug_file << std::endl;
		debug_file << "Measurements" << std::fixed << std::setprecision(16) << v_measMinusCompR_.at(block);
		debug_file << "Block " << block + 1 << " (Reverse)";
		debug_file << std::endl;
		debug_file << "At * V-inv" << std::fixed << std::setprecision(16) << v_AtVinvR_.at(block) << std::endl;
		debug_file << "Block " << block + 1 << " (Reverse)";
		debug_file << std::endl;
		debug_file << "Normals " << std::fixed << std::setprecision(16) << v_normalsR_.at(block) << std::endl;
		debug_file << "Block " << block + 1 << " (Reverse)";
		debug_file << std::endl;
		debug_file << "Precisions" << std::fixed << std::setprecision(16) << v_normalsR_.at(block) << std::endl;
		debug_file << "Block " << block + 1 << " (Reverse)";
		debug_file << std::endl;
		debug_file << "Corrections" << std::fixed << std::setprecision(16) << v_correctionsR_.at(block) << std::endl;
		debug_file.flush();

		dbg_file_mutex.unlock();

	}
}

		
void adjust_forward_thread::operator()()
{
	std::stringstream ss;
	UINT32 currentBlock;

	concurrentAdjustments.forward_run_started();
	
	try {
		// Set exception
		error_ = std::exception_ptr();

		for (currentBlock=0; currentBlock<main_adj_->blockCount_; ++currentBlock)
		{
			if (main_adj_->IsCancelled())
				return;

			//ss.str("");
			//ss << "1> Adjusted block " << currentBlock + 1 << " (forward, in isolation)... " << std::endl;
			//main_adj_->ThreadSafeWritetoAdjFile(ss.str());

			// Check if an exception was thrown in the reverse or combine threads
			if (rev_error || cmb_error)
				return;
			
			// At this point, whether first iteration or not, if currentBlock is the first block, 
			// the normals will have been initialised.  For all later blocks, the normals will contain
			// the contribution of junction station coordinates and variances from preceding blocks.
			// In either case, the block is ready for adjustment.  Junction station coordinates and 
			// variances are carried forward below

			// Least Squares Solution
			main_adj_->SolveTry(true, currentBlock);

			// Since Solve() can be computationally intensive, re-check
			// if an exception was thrown in the reverse or combine threads
			if (rev_error || cmb_error)
				return;
		
			main_adj_->UpdateEstimatesForward(currentBlock);
		
			/////////////////////////////////////
			//// protected write to adj file
			//ss.str("");
			//ss << "1> Forward adjustment (in isolation) of block " << currentBlock + 1 << " complete." << std::endl;

			//adj_file_mutex.lock();
			//main_adj_->ThreadSafeWritetoAdjFile(ss.str());
			//adj_file_mutex.unlock();
			/////////////////////////////////////

			// OK, now shrink matrices back to normal size
			main_adj_->ShrinkForwardMatrices(currentBlock);

			// This step is needed to carry coordinates and variance estimates for junctions only
			// for the forward run.  It is not needed for the combination stage
			main_adj_->CarryForwardJunctions(currentBlock, currentBlock+1);
	
			// Does this block need combining?
			if (!main_adj_->CombineRequired(currentBlock))
				continue;
				
			// Ok finished adjustment, set adjusted 
			// state of this block in the forward list
			concurrentAdjustments.set_forward_block_adjusted(currentBlock);

			// Has the reverseThread dealt with this block? 
			if (concurrentAdjustments.reverse_block_adjusted(currentBlock))
				// add this block # to the queue and notify the waiting 
				// combine thread
				combineAdjustmentQueue.push_and_notify(currentBlock);
		
			///////////////////////////////////
			// protected write to adj file
			//ss.str("");
			//ss << "1> Forward adjustment of block " << currentBlock << " complete." << std::endl;

			//main_adj_->ThreadSafeWritetoAdjFile(ss.str());		
			///////////////////////////////////
		}

		concurrentAdjustments.forward_run_finished();

		if (concurrentAdjustments.finished_all_runs())
			combineAdjustmentQueue.queue_exhausted();

		// notify all threads waiting on combineAdjustmentQueue
		combineAdjustmentQueue.notify_all();

		// Reset exception
		error_ = std::exception_ptr();

	} 
	catch (...) {
		error_ = std::current_exception();
		if (combineAdjustmentQueue.is_empty())
			combineAdjustmentQueue.push_and_notify(99999999);
		//combineAdjustmentQueue.notify_all();
		return;
	}

	// notify all threads waiting on combineAdjustmentQueue
	//combineAdjustmentQueue.notify_all();

	//std::this_thread::sleep_for(std::chrono::milliseconds(10));

	// notify all threads waiting on combineAdjustmentQueue
	//combineAdjustmentQueue.notify_all();	
}
	

void adjust_reverse_thread::operator()()
{
	std::stringstream ss;
	UINT32 currentBlock(main_adj_->blockCount_-1), block;

	concurrentAdjustments.reverse_run_started();
	
	try {
		// Set exception
		error_ = std::exception_ptr();

		for (block=0; block<main_adj_->blockCount_; ++block, --currentBlock)
		{
			if (main_adj_->IsCancelled())
				return;

			// Check if an exception was thrown in the reverse or combine threads
			if (fwd_error || cmb_error)
				return;
		
			// if this is a single block, then there is no need to perform a reverse adjustment
			if (!main_adj_->PrepareAdjustmentReverse(currentBlock, true))
				continue;
		
			/////////////////////////////////////
			//// protected write to adj file
			//ss.str("");
			//ss << "2> Adjusted block " << currentBlock + 1 << " (reverse, in isolation)... " << std::endl;
			//main_adj_->ThreadSafeWritetoAdjFile(ss.str());
			/////////////////////////////////////

			// Backup normals prior to inversion for re-use in combination
			// adjustment... only if a combination is required
			main_adj_->BackupNormals(currentBlock, true);

			// At this point, whether first iteration or not, if currentBlock is the last block, 
			// the normals will have been initialised.  For all subsequent blocks, the normals will contain
			// the contribution of junction station coordinates and variances from preceding blocks.
			// In either case, the block is ready for adjustment.  Junction station coordinates and 
			// variances are carried in the reverse direction below

			// Least Squares Solution
			main_adj_->SolveMTTry(true, currentBlock);

			// Since SolveMT() can be computationally intensive, re-check
			// if an exception was thrown in the forward or combine threads
			if (fwd_error || cmb_error)
				return;
		
			main_adj_->UpdateEstimatesReverse(currentBlock, true);

			/////////////////////////////////////
			//// protected write to adj file
			//ss.str("");
			//ss << "2> Reverse adjustment (in isolation) of block " << currentBlock + 1 << " complete." << std::endl;

			//adj_file_mutex.lock();
			//main_adj_->ThreadSafeWritetoAdjFile(ss.str());
			//adj_file_mutex.unlock();
			/////////////////////////////////////

			// Now, carry the estimated junction station coordinates and variances 
			// to the next block (to be used in the next loop), except when currentBlock
			// is the first block and currentBlock-1 is an isolated block.
			//
			// Remember - the junction station estimates and variances of the next block
			// (obtained during the forward pass) were copied during the forward pass
			// in CarryStnEstimatesandVariancesForward(..), so no need to re-copy.
			main_adj_->CarryReverseJunctions(currentBlock, currentBlock-1, true);
	
			// Does this block need combining?
			if (!main_adj_->CombineRequired(currentBlock))
			{
				if (main_adj_->FirstBlock(currentBlock))
					main_adj_->UpdateEstimatesFinal(currentBlock);
			
				// This block does not need combination - continue
				continue;
			}

			// Ok finished adjustment, set adjusted 
			// state of this block in the reverse list
			concurrentAdjustments.set_reverse_block_adjusted(currentBlock);			

			// Has the forwardThread dealt with this block? 
			if (concurrentAdjustments.forward_block_adjusted(currentBlock))
				// add this block # to the queue and notify the waiting 
				// combine thread
				combineAdjustmentQueue.push_and_notify(currentBlock);
		
			
			///////////////////////////////////
			// protected write to adj file
			//ss.str("");
			//ss << "2> Reverse adjustment of block " << currentBlock << " complete." << std::endl;

			//main_adj_->ThreadSafeWritetoAdjFile(ss.str());
			///////////////////////////////////
		}

		concurrentAdjustments.reverse_run_finished();

		if (concurrentAdjustments.finished_all_runs())
			combineAdjustmentQueue.queue_exhausted();

		// Reset exception
		error_ = std::exception_ptr();

	} 
	catch (...) {
		error_ = std::current_exception();
		if (combineAdjustmentQueue.is_empty())
			combineAdjustmentQueue.push_and_notify(99999999);
		return;
	}
}
	

void adjust_process_combine_thread::operator()()
{
	std::stringstream ss;
	UINT32 pseudomsrJSLCount;
	UINT32 currentBlock;

	try {
		while (true)		
		{
			if (main_adj_->IsCancelled())
				return;

			if (combineAdjustmentQueue.is_queue_exhausted())
				break;
			
			// Have any other adjustments failed?
			if (combineAdjustmentExceptionThrown(cmb_errors_))
				return;

			// Wait here until blocks have been placed on the queue
			if (!combineAdjustmentQueue.front_and_pop(currentBlock))
			{
				std::this_thread::sleep_for(std::chrono::milliseconds(2));
				continue;
			}

			// Have any other adjustments failed?
			if (combineAdjustmentExceptionThrown(cmb_errors_))
				return;

			main_adj_->isCombining_ = true;
			main_adj_->SetcurrentBlock(currentBlock);

			if (main_adj_->PrepareAdjustmentCombine(currentBlock, pseudomsrJSLCount, true))
			{
				// Least Squares Solution
				main_adj_->SolveMTTry(true, currentBlock);
				main_adj_->UpdateEstimatesCombine(currentBlock, pseudomsrJSLCount);
				main_adj_->UpdateEstimatesFinal(currentBlock);
			}
		}

		error_ = std::exception_ptr();
	}
	catch (...) {
		error_ = std::current_exception();
		return;
	}
}
	

void adjust_combine_thread::operator()()
{
	main_adj_->isCombining_ = false;
	
	// Wait here until notification has been received
	combineAdjustmentQueue.wait_if_queue_is_empty();

	// Has either a forward or reverse adjustment failed?
	if (fwd_error || rev_error)
		return;

	// OK, at this point, there is a block # in the queue - let's combine!

	// Get number of cores available
	UINT32 thread_id, cores(std::thread::hardware_concurrency());

	// Set up exception pointers
	std::vector<std::exception_ptr> cmb_errors;
	cmb_errors.resize(cores);

	// Create the thread pool
#if defined(__ICC) || defined(__INTEL_COMPILER)		// Intel compiler
	std::shared_ptr<std::thread> c;
	std::vector< std::shared_ptr<std::thread> > mt_combine_threads_sp;
#else
	std::vector<std::thread> mt_combine_threads;
#endif	

	for (thread_id=0; thread_id<cores; ++thread_id)
	{

#if defined(__ICC) || defined(__INTEL_COMPILER)			// Intel compiler
		c.reset(new thread(
			adjust_process_combine_thread(
				main_adj_,								// pointer to main adjustment object
				thread_id,								// this new thread's ID (1..n, where n = # cores)
				std::ref(cmb_errors.at(thread_id)),	// reference to this thread's exception_ptr
				std::ref(cmb_errors))));				// reference to all other threads' exception_ptrs

		mt_combine_threads_sp.push_back(c);
#else
		mt_combine_threads.push_back(std::thread(
			adjust_process_combine_thread(
				main_adj_,								// pointer to main adjustment object
				thread_id,								// this new thread's ID (1..n, where n = # cores)
				std::ref(cmb_errors.at(thread_id)),	// reference to this thread's exception_ptr
				std::ref(cmb_errors))));				// reference to all other threads' exception_ptrs
#endif			
		
	}

#if defined(__ICC) || defined(__INTEL_COMPILER)		// Intel compiler
	for_each(mt_combine_threads_sp.begin(), mt_combine_threads_sp.end(), std::mem_fn(&thread::join));
#else
	for_each(mt_combine_threads.begin(), mt_combine_threads.end(), std::mem_fn(&std::thread::join));
#endif	

	// cmb_error is set during combineThreadProcessor
	for (UINT32 err=0; err<cores; ++err)
	{
		if (cmb_errors.at(err))
		{
			error_ = cmb_errors.at(err);
			return;
		}
	}	
}



// Prepare adjustment

void adjust_process_prepare_thread::operator()()
{
	std::stringstream ss;
	UINT32 currentBlock;

	try {
		while (true)		
		{
			if (prepareAdjustmentQueue.is_queue_exhausted())
				break;
			
			// Have any other adjustments failed?
			if (prepareAdjustmentExceptionThrown(prep_errors_))
				return;

			// Wait here until blocks have been placed on the queue
			if (!prepareAdjustmentQueue.front_and_pop(currentBlock))
			{
				std::this_thread::sleep_for(std::chrono::milliseconds(2));
				continue;
			}

			// Have any other adjustments failed?
			if (prepareAdjustmentExceptionThrown(prep_errors_))
				return;

			main_adj_->SetcurrentBlock(currentBlock);
			main_adj_->CreateMeasurementTally(currentBlock);
			main_adj_->PrepareAdjustmentBlock(currentBlock, thread_id_);

			if (prepareAdjustmentQueue.is_empty())
				prepareAdjustmentQueue.queue_exhausted();
		}

		error_ = std::exception_ptr();
	}
	catch (...) {
		error_ = std::current_exception();
	}
}
	

void adjust_prepare_thread::operator()()
{	
	// load up the queue with the block ids
	prepareAdjustmentQueue.push_data(createIncrementingIntegerVector<UINT32>(main_adj_->blockCount()));

	// OK, at this point, there is a block # in the queue - let's combine!
	// Get number of cores available
	UINT32 thread_id, cores(std::thread::hardware_concurrency());

	// Set up exception pointers
	std::vector<std::exception_ptr> prep_errors;
	prep_errors.resize(cores);

	// Create the thread pool
#if defined(__ICC) || defined(__INTEL_COMPILER)		// Intel compiler
	std::shared_ptr<std::thread> p;
	std::vector< std::shared_ptr<std::thread> > mt_prepare_threads_sp;
#else
	std::vector<std::thread> mt_prepare_threads;
#endif	
	

	for (thread_id=0; thread_id<cores; ++thread_id)
	{
#if defined(__ICC) || defined(__INTEL_COMPILER)		// Intel compiler
		p.reset(new thread(
			adjust_process_prepare_thread(
				main_adj_,				// pointer to main adjustment object
				thread_id,				// this new thread's ID (1..n, where n = # cores)
				std::ref(prep_errors.at(thread_id)),	// reference to this thread's exception_ptr
				std::ref(prep_errors))));		// reference to all other threads' exception_ptrs

		mt_prepare_threads_sp.push_back(p);
#else
		mt_prepare_threads.push_back(std::thread(
			adjust_process_prepare_thread(
				main_adj_,				// pointer to main adjustment object
				thread_id,				// this new thread's ID (1..n, where n = # cores)
				std::ref(prep_errors.at(thread_id)),	// reference to this thread's exception_ptr
				std::ref(prep_errors))));		// reference to all other threads' exception_ptrs
#endif	
	}

#if defined(__ICC) || defined(__INTEL_COMPILER)		// Intel compiler
	for_each(mt_prepare_threads_sp.begin(), mt_prepare_threads_sp.end(), std::mem_fn(&thread::join));
#else
	for_each(mt_prepare_threads.begin(), mt_prepare_threads.end(), std::mem_fn(&std::thread::join));
#endif	

	// prep_error is set during prepareThreadProcessor
	for (UINT32 err=0; err<cores; ++err)
	{
		if (prep_errors.at(err))
		{
			error_ = prep_errors.at(err);
			return;
		}
	}
}


void dna_adjust::PrepareAdjustmentMultiThread()
{
	prepareAdjustmentQueue.reset_blocks_coming();
	
	std::thread mt_prepare_thread(adjust_prepare_thread(this, std::ref(prep_error)));

	// Start the forward, reverse and combine threads, which commences the
	// network adjustment
	mt_prepare_thread.join();

	// Were any exceptions thrown?  If so, re-throw and let
	// test stub handle the exception
	if (prep_error)
		// exception thrown in one of the combination adjustments
		std::rethrow_exception(prep_error);
}



}	// namespace networkadjust
}	// namespace dynadjust

