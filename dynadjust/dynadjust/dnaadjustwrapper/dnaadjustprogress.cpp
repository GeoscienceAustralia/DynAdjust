//============================================================================
// Name         : dnaadjustprogress.cpp
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
// Description  : DynAdjust Adjustment executable interface and progress
//============================================================================

#include <dynadjust/dnaadjustwrapper/dnaadjustprogress.hpp>
#include <dynadjust/dnaadjust/dnaadjust.hpp>

/// \cond
#include <time.h>
#include <algorithm>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <thread>
#include <vector>
/// \endcond

using namespace dynadjust::networkadjust;
using namespace dynadjust::exception;

bool running;
std::mutex cout_mutex;

dna_adjust_thread::dna_adjust_thread(dna_adjust* dnaAdj, project_settings* p,
		_ADJUST_STATUS_* adjustStatus)
		: _dnaAdj(dnaAdj), _p(p), _adjustStatus(adjustStatus) 
{
}

void dna_adjust_thread::operator()()
{
	running = true;
		
	if (!prepareAdjustment())
	{
		running = false;
		return;
	}
	if (!processAdjustment())
	{
		running = false;
		return;
	}

	// Do addition things here...

	running = false;
}

bool dna_adjust_thread::prepareAdjustment()
{
	try {
		*_adjustStatus = ADJUST_EXCEPTION_RAISED;
		_dnaAdj->PrepareAdjustment(*_p);
		*_adjustStatus = ADJUST_SUCCESS;
		return true;
	}
	catch (const NetAdjustException& e) {
		handlePrepareAdjustError(e.what());
		_dnaAdj->SetExceptionRaised();
	}
	catch (const NetMemoryException& e) {
		handlePrepareAdjustError(e.what());
		_dnaAdj->SetExceptionRaised();
	}
	catch (const std::runtime_error& e) {
		handlePrepareAdjustError(e.what());
		_dnaAdj->SetExceptionRaised();
	}
	catch (const std::bad_alloc& e) {
		handlePrepareAdjustError(e.what());
		_dnaAdj->SetExceptionRaised();
	}
	catch (const std::out_of_range& e) {
		handlePrepareAdjustError(e.what());
		_dnaAdj->SetExceptionRaised();
	}
	catch (const std::exception& e) {
		std::string err = std::string("Standard exception: ") + e.what();
		handlePrepareAdjustError(err);
		_dnaAdj->SetExceptionRaised();
	}
	catch (...) {
		std::string err("Undefined error.\n  This could be the result of insufficient memory or poorly configured data.");
		handlePrepareAdjustError(err);
		_dnaAdj->SetExceptionRaised();
	}
	*_adjustStatus = ADJUST_EXCEPTION_RAISED;
	return false;
}

bool dna_adjust_thread::processAdjustment()
{
	try {
		*_adjustStatus = _dnaAdj->AdjustNetwork();
		*_adjustStatus = ADJUST_SUCCESS;
		return true;
	}
	catch (const NetAdjustException& e) {
		handleProcessAdjustError(e.what());
		_dnaAdj->SetExceptionRaised();
	}
	catch (const NetMemoryException& e) {
		handleProcessAdjustError(e.what());
		_dnaAdj->SetExceptionRaised();
	}
	catch (const std::runtime_error& e) {
		handleProcessAdjustError(e.what());
		_dnaAdj->SetExceptionRaised();
	}
	catch (...) {
		std::string err("- Error: Undefined error.\n  This could be the result of insufficient memory or poorly configured data.");
		handleProcessAdjustError(err);
		_dnaAdj->SetExceptionRaised();
	}
	*_adjustStatus = ADJUST_EXCEPTION_RAISED;
	return false;
}

void dna_adjust_thread::handlePrepareAdjustError(const std::string& error_msg)
{
	running = false;
	std::stringstream ss;
	ss << /*PROGRESS_BACKSPACE_12 << */setw(PROGRESS_ADJ_BLOCK_12) << " " << std::endl
		<< "- Error: " << std::endl << "  " << error_msg << std::endl;
	coutMessage(ss.str());
}

void dna_adjust_thread::handleProcessAdjustError(const std::string& error_msg)
{
	running = false;
	std::this_thread::sleep_for(std::chrono::milliseconds(50));
	printErrorMsg(error_msg);
}

void dna_adjust_thread::printErrorMsg(const std::string& error_msg)
{
	std::stringstream ss, sst;
	
	switch (_p->a.adjust_mode)
	{
	case SimultaneousMode:
		// "  Iteration " = 12
		// setw(2)        =  2
		// -------------------
		//                  14
		sst << "  Iteration " << std::right << std::setw(2) << std::fixed << std::setprecision(0) << _dnaAdj->CurrentIteration();
		ss << PROGRESS_BACKSPACE_14 << std::setw(PROGRESS_ADJ_BLOCK_14) << std::right << sst.str();
		break;
	case Phased_Block_1Mode:
	case PhasedMode:
		// "  Iteration " = 12
		// setw(2)        =  2
		// ", block"      =  7
		// setw(6)        =  6
		// "(forward) "   = 10
		// -------------------
		//                  37
		sst << "  Iteration " << 
			std::right << std::setw(2) << std::fixed << std::setprecision(0) << _dnaAdj->CurrentIteration() << ", block" << std::right << std::setw(6) << std::fixed << std::setprecision(0) << 
			_dnaAdj->CurrentBlock() + 1 << (_dnaAdj->processingForward() ? " (forward)" : " (reverse)");
		ss << PROGRESS_BACKSPACE_37 << std::setw(37) << std::right << sst.str();
	}

	ss << std::endl << "- Error: Cannot compute the least squares estimates.  Reason:\n  " << error_msg << std::endl;
	coutMessage(ss);
}

void dna_adjust_thread::coutMessage(std::stringstream& message)
{
	coutMessage(message.str());
}

void dna_adjust_thread::coutMessage(const std::string& message)
{
	cout_mutex.lock();
	std::cout.flush();
	std::cout << message;
	std::cout.flush();
	cout_mutex.unlock();
}








dna_adjust_progress_thread::dna_adjust_progress_thread(dna_adjust* dnaAdj, project_settings* p)
		: _dnaAdj(dnaAdj), _p(p) {}

void dna_adjust_progress_thread::prepareAdjustment()
{
	if (_p->g.quiet)
		// Nothing to be displayed
		return;
	
	if (_p->a.report_mode)
		// Nothing to be displayed
		return;

	std::stringstream ss;
	ss << "+ Preparing for ";

	switch (_p->a.adjust_mode)
	{
	case Phased_Block_1Mode:
	case PhasedMode:
		ss << "a " << _dnaAdj->blockCount() << " block ";
		break;
	}
			
	ss << "adjustment... ";		
	coutMessage(ss.str());
	ss.str("");
	
	UINT32 block, currentBlock(0);
	std::stringstream sst;
	bool first_time(true);
	
	switch (_p->a.adjust_mode)
	{
	case SimultaneousMode:
		while (running && _dnaAdj->IsPreparing())
		{
			std::this_thread::sleep_for(std::chrono::milliseconds(40));
		}
		
		if (_dnaAdj->ExceptionRaised())
			return;

		ss << " done." << std::endl;
		coutMessage(ss.str());
		
		break;
	case Phased_Block_1Mode:
	case PhasedMode:
		while (running && _dnaAdj->IsPreparing())
		{
			block = _dnaAdj->CurrentBlock();

			if (block != currentBlock)
			{
				ss.str("");
				ss << " block " << std::left << std::setw(5) << std::fixed << std::setprecision(0) << _dnaAdj->CurrentBlock() + 1;

				sst.str("");
				if (first_time)
				{
					sst << std::setw(PROGRESS_ADJ_BLOCK_12) << std::left << " ";
					first_time = false;
				}
				sst << PROGRESS_BACKSPACE_12 << std::setw(PROGRESS_ADJ_BLOCK_12) << std::left << ss.str();
				coutMessage(sst.str());
				currentBlock = block;
			}
			std::this_thread::sleep_for(std::chrono::milliseconds(40));
		}
		
		// Check if an exception occurred before printing "done"
		if (_dnaAdj->ExceptionRaised())
			return;
		
		ss.str("");
		ss << " done.";
		sst.str("");
		if (!first_time)
			sst << PROGRESS_BACKSPACE_12;
		sst << std::setw(PROGRESS_ADJ_BLOCK_12) << std::left << ss.str() << std::endl;
		coutMessage(sst.str());
		
	}
}

void dna_adjust_progress_thread::processAdjustment()
{
	if (_p->g.quiet)
		// Nothing to be displayed
		return;
	
	if (_p->a.max_iterations > 0)
		coutMessage(std::string("+ Adjusting network...\n"));

	UINT32 block, currentBlock(0);
	UINT32 currentIteration(0);
	std::stringstream ss, sst;
	bool first_time(true);

	switch (_p->a.adjust_mode)
	{
	case SimultaneousMode:
		
		// Check if an exception was raised during preparation
		if (_dnaAdj->ExceptionRaised())
			return;
			
		while (_dnaAdj->IsAdjusting())
		{
			// use message bank
			while (_dnaAdj->NewMessagesAvailable())
			{
				if (!_dnaAdj->GetMessageIteration(currentIteration))
					break;
				ss.str("");
				ss << "  Iteration " << std::right << std::setw(2) << std::fixed << std::setprecision(0) << currentIteration;
				ss << ", max station corr: " << std::right << std::setw(PROGRESS_ADJ_BLOCK_12) <<
					_dnaAdj->GetMaxCorrection(currentIteration) << std::endl;
					
				coutMessage(ss.str());
			}
			std::this_thread::sleep_for(std::chrono::milliseconds(80));
		}
		
		break;
	case Phased_Block_1Mode:
	case PhasedMode:
		
		// Check if an exception was raised during preparation
		if (_dnaAdj->ExceptionRaised())
			return;
		
		while (_dnaAdj->IsAdjusting())
		{
			block = _dnaAdj->CurrentBlock();
			
			// use message bank
			while (_dnaAdj->NewMessagesAvailable())
			{
				if (!_dnaAdj->GetMessageIteration(currentIteration))
					break;
						
				ss.str("");
				ss << "  Iteration " << std::right << std::setw(2) << std::fixed << std::setprecision(0) << currentIteration;
				ss << ", max station corr: " << std::right << std::setw(PROGRESS_ADJ_BLOCK_12) << _dnaAdj->GetMaxCorrection(currentIteration) << std::endl;
				
				sst.str("");
				if (first_time)
					sst << std::setw(PROGRESS_ADJ_BLOCK_28) << std::left << " ";
				sst << PROGRESS_BACKSPACE_28 << std::setw(PROGRESS_ADJ_BLOCK_28) << std::left << ss.str();
				coutMessage(sst.str());
				first_time = true;
			}

			std::this_thread::sleep_for(std::chrono::milliseconds(40));
					
			// print new block to screen when adjusting only
			if (block != currentBlock && _dnaAdj->IsAdjusting())
			{						
				ss.str("");
				ss << "  Iteration " << std::right << std::setw(2) << std::fixed << std::setprecision(0) << _dnaAdj->CurrentIteration();

				if (_p->a.multi_thread && !_dnaAdj->processingCombine())
					ss << std::left << std::setw(13) << ", adjusting...";
				else
					ss << ", block " << std::left << std::setw(6) << std::fixed << std::setprecision(0) << _dnaAdj->CurrentBlock() + 1;
						
				sst.str("");
				if (first_time)
				{
					sst << std::setw(PROGRESS_ADJ_BLOCK_28) << std::left << " ";
					first_time = false;
				}
				
				sst << PROGRESS_BACKSPACE_28 << std::setw(PROGRESS_ADJ_BLOCK_28) << std::left << ss.str();
				coutMessage(sst.str());

				currentBlock = block;
			}
		}
	}
}

void dna_adjust_progress_thread::operator()()
{
	if (_p->g.quiet)
		// Nothing to be displayed
		return;

	// cout messages related to preparing adjustment
	// prepareAdjustment() will continue while the following is true:
	//	running
	//	_dnaAdj->IsPreparing()
	//	!_dnaAdj->ExceptionRaised()
	prepareAdjustment();
	
	if (_dnaAdj->ExceptionRaised())
		return;
	
	processAdjustment();	
}

void dna_adjust_progress_thread::coutMessage(const std::string& message)
{
	cout_mutex.lock();
	std::cout.flush();
	std::cout << message;
	std::cout.flush();
	cout_mutex.unlock();
}

