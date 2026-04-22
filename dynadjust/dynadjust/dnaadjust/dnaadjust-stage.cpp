//============================================================================
// Name         : dnaadjust-stage.cpp
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
// Description  : DynAdjust Network Adjustment (stage) library
//============================================================================

#include <dynadjust/dnaadjust/dnaadjust.hpp>

namespace dynadjust { 
namespace networkadjust {

// Prepare mapped regions for staged adjustments that read from
// existing mapped files.
void dna_adjust::PrepareMappedRegions(const UINT32& block)
{
	if (projectSettings_.a.adjust_mode == SimultaneousMode)
		return;
	if (projectSettings_.a.multi_thread)
		return;

	// Set memory map region offsets for all 
	// matrices associated with this block
	SetRegionOffsets(block);

	// Set memory mapped file regions
	// NOTE - previously created files must exist and match the 
	// dimensions of the current matrix sizes			
	try {
		PrepareMemoryMapRegions(block);
	}
	catch (boost::interprocess::interprocess_exception& e){
		std::stringstream ss;
		ss << "PrepareMappedRegions() terminated while creating memory map" << std::endl;
		ss << "  regions from .mtx stage files. Details:\n  " << e.what() << std::endl << std::endl;
		ss << "  Please ensure the .mtx stage files from a previous staged" << std::endl;
		ss << "  adjustment exist, or re-run the adjustment using the" << std::endl;
		ss << "  --" << RECREATE_STAGE_FILES << " option." << std::endl;
		adj_file << std::endl << "- Error: " << ss.str() << std::endl;
		throw std::runtime_error(ss.str());
	}
}
	

void dna_adjust::SetRegionOffsets(const UINT32& block, const int file_count, ...)
{
	if (file_count == 0)
	{
		// deserialise all.  That is, call this function again, but with 
		// arguments for all files
		SetRegionOffsets(block, 12, 
			sf_normals_r, sf_meas_minus_comp,
			sf_estimated_stns, sf_original_stns, sf_rigorous_stns,
			sf_junction_vars, sf_junction_vars_f, sf_junction_ests_f, sf_junction_ests_r,
			sf_rigorous_vars, sf_prec_adj_msrs, sf_corrections);
		return;
	}

	va_list vlist;
	va_start(vlist, file_count);

	for (UINT16 file(0); file<file_count; ++file)
	{
		switch (va_arg(vlist, int))
		{
		case sf_normals_r:
			normalsR_map_.addblockMapRegion(block_map_t(v_normalsR_.at(block).get_size()));
			SetRegionOffset(normalsR_map_);
			break;
		case sf_meas_minus_comp:
			measMinusComp_map_.addblockMapRegion(block_map_t(v_measMinusComp_.at(block).get_size()));
			SetRegionOffset(measMinusComp_map_);
			break;
		case sf_estimated_stns:
			estimatedStations_map_.addblockMapRegion(block_map_t(v_estimatedStations_.at(block).get_size()));
			SetRegionOffset(estimatedStations_map_);
			break;
		case sf_original_stns:
			originalStations_map_.addblockMapRegion(block_map_t(v_originalStations_.at(block).get_size()));
			SetRegionOffset(originalStations_map_);
			break;
		case sf_rigorous_stns:
			rigorousStations_map_.addblockMapRegion(block_map_t(v_rigorousStations_.at(block).get_size()));
			SetRegionOffset(rigorousStations_map_);
			break;
		case sf_junction_vars:
			junctionVariances_map_.addblockMapRegion(block_map_t(v_junctionVariances_.at(block).get_size()));
			SetRegionOffset(junctionVariances_map_);
			break;
		case sf_junction_vars_f:
			junctionVariancesFwd_map_.addblockMapRegion(block_map_t(v_junctionVariancesFwd_.at(block).get_size()));
			SetRegionOffset(junctionVariancesFwd_map_);
			break;
		case sf_junction_ests_f:
			junctionEstimatesFwd_map_.addblockMapRegion(block_map_t(v_junctionEstimatesFwd_.at(block).get_size()));
			SetRegionOffset(junctionEstimatesFwd_map_);
			break;
		case sf_junction_ests_r:
			junctionEstimatesRev_map_.addblockMapRegion(block_map_t(v_junctionEstimatesRev_.at(block).get_size()));
			SetRegionOffset(junctionEstimatesRev_map_);
			break;
		case sf_rigorous_vars:
			rigorousVariances_map_.addblockMapRegion(block_map_t(v_rigorousVariances_.at(block).get_size()));
			SetRegionOffset(rigorousVariances_map_);
			break;
		case sf_prec_adj_msrs:
			precAdjMsrs_map_.addblockMapRegion(block_map_t(v_precAdjMsrsFull_.at(block).get_size()));
			SetRegionOffset(precAdjMsrs_map_);
			break;
		case sf_corrections:
			corrections_map_.addblockMapRegion(block_map_t(v_corrections_.at(block).get_size()));
			SetRegionOffset(corrections_map_);
			break;
		}
	}
	va_end(vlist);
}
	

void dna_adjust::SetRegionOffset(vmat_file_map& file_map)
{
	// Get iterator to the second last mapped region
	_it_block_map _it_bmap_prev(file_map.vblockMapRegions_.end() - 1);
	if (_it_bmap_prev == file_map.vblockMapRegions_.begin())
		return;
	_it_bmap_prev--;
	// Get iterator to the last mapped region
	_it_block_map _it_bmap(file_map.vblockMapRegions_.end() - 1);
	// Set region offset for the last mapped region
	_it_bmap->SetRegionOffset(_it_bmap_prev->GetCumulativeRegionOffset());	
}
	

void dna_adjust::DeserialiseBlockFromMappedFile(const UINT32& block, const int file_count, ...)
{
	if (file_count == 0)
	{
		// deserialise all.  That is, call this function again, but with 
		// arguments for all files
		DeserialiseBlockFromMappedFile(block, 15, 
			sf_normals, sf_normals_r, sf_atvinv, sf_design, sf_meas_minus_comp,
			sf_estimated_stns, sf_original_stns, sf_rigorous_stns,
			sf_junction_vars, sf_junction_vars_f, sf_junction_ests_f, sf_junction_ests_r,
			sf_rigorous_vars, sf_prec_adj_msrs, sf_corrections);
		return;
	}

	va_list vlist;
	va_start(vlist, file_count);
	
	void* addr;

	for (UINT16 file(0); file<file_count; ++file)
	{
		switch (va_arg(vlist, int))
		{
		case sf_normals:
			v_normals_.at(block).allocate();
			break;
		case sf_normals_r:
			addr = normalsR_map_.GetBlockRegionAddr(block);
			v_normalsR_.at(block).AttachMappedFileRegion(addr);
			break;
		case sf_atvinv:
			v_AtVinv_.at(block).allocate();
			break;
		case sf_design:
			v_design_.at(block).allocate();
			break;
		case sf_meas_minus_comp:
			addr = measMinusComp_map_.GetBlockRegionAddr(block);
			v_measMinusComp_.at(block).AttachMappedFileRegion(addr);
			break;
		case sf_estimated_stns:
			addr = estimatedStations_map_.GetBlockRegionAddr(block);
			v_estimatedStations_.at(block).AttachMappedFileRegion(addr);
			break;
		case sf_original_stns:
			addr = originalStations_map_.GetBlockRegionAddr(block);
			v_originalStations_.at(block).AttachMappedFileRegion(addr);
			break;
		case sf_rigorous_stns:
			addr = rigorousStations_map_.GetBlockRegionAddr(block);
			v_rigorousStations_.at(block).AttachMappedFileRegion(addr);
			break;
		case sf_junction_vars:
			addr = junctionVariances_map_.GetBlockRegionAddr(block);
			v_junctionVariances_.at(block).AttachMappedFileRegion(addr);
			break;
		case sf_junction_vars_f:
			addr = junctionVariancesFwd_map_.GetBlockRegionAddr(block);
			v_junctionVariancesFwd_.at(block).AttachMappedFileRegion(addr);
			break;
		case sf_junction_ests_f:
			addr = junctionEstimatesFwd_map_.GetBlockRegionAddr(block);
			v_junctionEstimatesFwd_.at(block).AttachMappedFileRegion(addr);
			break;
		case sf_junction_ests_r:
			addr = junctionEstimatesRev_map_.GetBlockRegionAddr(block);
			v_junctionEstimatesRev_.at(block).AttachMappedFileRegion(addr);
			break;
		case sf_rigorous_vars:
			addr = rigorousVariances_map_.GetBlockRegionAddr(block);
			v_rigorousVariances_.at(block).AttachMappedFileRegion(addr);
			break;
		case sf_prec_adj_msrs:
			addr = precAdjMsrs_map_.GetBlockRegionAddr(block);
			v_precAdjMsrsFull_.at(block).AttachMappedFileRegion(addr);
			break;
		case sf_corrections:
			addr = corrections_map_.GetBlockRegionAddr(block);
			v_corrections_.at(block).AttachMappedFileRegion(addr);

			if (v_blockMeta_.at(block)._blockLast)
				v_correctionsR_.at(block).allocate();

			break;
		}
	}
	va_end(vlist);
}

void dna_adjust::SerialiseBlockToMappedFile(const UINT32& block, const int file_count, ...)
{
	if (file_count == 0)
	{
		// serialise all.  That is, call this function again, but with 
		// arguments for all files
		SerialiseBlockToMappedFile(block, 12,
			sf_normals_r, sf_meas_minus_comp,
			sf_estimated_stns, sf_original_stns, sf_rigorous_stns,
			sf_junction_vars, sf_junction_vars_f, sf_junction_ests_f, sf_junction_ests_r,
			sf_rigorous_vars, sf_prec_adj_msrs, sf_corrections);
		return;
	}

	va_list vlist;
	va_start(vlist, file_count);
	
	void* addr;

	for (UINT16 file(0); file<file_count; ++file)
	{
		switch (va_arg(vlist, int))
		{
		case sf_normals_r:
			addr = normalsR_map_.GetBlockRegionAddr(block);
			v_normalsR_.at(block).WriteMappedFileRegion(addr);
			break;
		case sf_meas_minus_comp:
			addr = measMinusComp_map_.GetBlockRegionAddr(block);
			v_measMinusComp_.at(block).WriteMappedFileRegion(addr);
			break;
		case sf_estimated_stns:
			addr = estimatedStations_map_.GetBlockRegionAddr(block);
			v_estimatedStations_.at(block).WriteMappedFileRegion(addr);
			break;
		case sf_original_stns:
			addr = originalStations_map_.GetBlockRegionAddr(block);
			v_originalStations_.at(block).WriteMappedFileRegion(addr);
			break;
		case sf_rigorous_stns:
			addr = rigorousStations_map_.GetBlockRegionAddr(block);
			v_rigorousStations_.at(block).WriteMappedFileRegion(addr);
			break;
		case sf_junction_vars:
			addr = junctionVariances_map_.GetBlockRegionAddr(block);
			v_junctionVariances_.at(block).WriteMappedFileRegion(addr);
			break;
		case sf_junction_vars_f:
			addr = junctionVariancesFwd_map_.GetBlockRegionAddr(block);
			v_junctionVariancesFwd_.at(block).WriteMappedFileRegion(addr);
			break;
		case sf_junction_ests_f:
			addr = junctionEstimatesFwd_map_.GetBlockRegionAddr(block);
			v_junctionEstimatesFwd_.at(block).WriteMappedFileRegion(addr);
			break;
		case sf_junction_ests_r:
			addr = junctionEstimatesRev_map_.GetBlockRegionAddr(block);
			v_junctionEstimatesRev_.at(block).WriteMappedFileRegion(addr);
			break;
		case sf_rigorous_vars:
			addr = rigorousVariances_map_.GetBlockRegionAddr(block);
			v_rigorousVariances_.at(block).WriteMappedFileRegion(addr);
			break;
		case sf_prec_adj_msrs:
			addr = precAdjMsrs_map_.GetBlockRegionAddr(block);
			v_precAdjMsrsFull_.at(block).WriteMappedFileRegion(addr);
			break;
		case sf_corrections:
			addr = corrections_map_.GetBlockRegionAddr(block);
			v_corrections_.at(block).WriteMappedFileRegion(addr);
			break;
		}
	}
	va_end(vlist);
}
	

void dna_adjust::ReserveBlockMapRegions(const int file_count, ...)
{
	if (file_count == 0)
	{
		// serialise all.  That is, call this function again, but with 
		// arguments for all files
		ReserveBlockMapRegions(12, 
			sf_normals_r, sf_meas_minus_comp,
			sf_estimated_stns, sf_original_stns, sf_rigorous_stns,
			sf_junction_vars, sf_junction_vars_f, sf_junction_ests_f, sf_junction_ests_r,
			sf_rigorous_vars, sf_prec_adj_msrs, sf_corrections);
		return;
	}

	va_list vlist;
	va_start(vlist, file_count);

	for (UINT16 file(0); file<file_count; ++file)
	{
		switch (va_arg(vlist, int))
		{
		case sf_normals_r:
			normalsR_map_.reserveblockMapRegions(blockCount_);
			break;
		case sf_meas_minus_comp:
			measMinusComp_map_.reserveblockMapRegions(blockCount_);
			break;
		case sf_estimated_stns:
			estimatedStations_map_.reserveblockMapRegions(blockCount_);
			break;
		case sf_original_stns:
			originalStations_map_.reserveblockMapRegions(blockCount_);
			break;
		case sf_rigorous_stns:
			rigorousStations_map_.reserveblockMapRegions(blockCount_);
			break;
		case sf_junction_vars:
			junctionVariances_map_.reserveblockMapRegions(blockCount_);
			break;
		case sf_junction_vars_f:
			junctionVariancesFwd_map_.reserveblockMapRegions(blockCount_);
			break;
		case sf_junction_ests_f:
			junctionEstimatesFwd_map_.reserveblockMapRegions(blockCount_);
			break;
		case sf_junction_ests_r:
			junctionEstimatesRev_map_.reserveblockMapRegions(blockCount_);
			break;
		case sf_rigorous_vars:
			rigorousVariances_map_.reserveblockMapRegions(blockCount_);
			break;
		case sf_prec_adj_msrs:
			precAdjMsrs_map_.reserveblockMapRegions(blockCount_);
			break;
		case sf_corrections:
			corrections_map_.reserveblockMapRegions(blockCount_);
			break;
		}
	}
	va_end(vlist);
}
	

void dna_adjust::OpenStageFileStreams(const int file_count, ...)
{
	if (file_count == 0)
	{
		// serialise all.  That is, call this function again, but with 
		// arguments for all files
		OpenStageFileStreams(12, 
			sf_normals_r, sf_meas_minus_comp,
			sf_estimated_stns, sf_original_stns, sf_rigorous_stns,
			sf_junction_vars, sf_junction_vars_f, sf_junction_ests_f, sf_junction_ests_r,
			sf_rigorous_vars, sf_prec_adj_msrs, sf_corrections);
		return;
	}

	std::stringstream ss;
	if (projectSettings_.a.stage_path.empty())
		ss << projectSettings_.g.output_folder;
	else
		ss << projectSettings_.a.stage_path;
	ss << FOLDER_SLASH << projectSettings_.g.network_name << "-";
	std::string filePath(ss.str());

	v_stageFileStreams_.clear();
	std::string fullfilePath;

	va_list vlist;
	va_start(vlist, file_count);
	
	try {		

		for (UINT16 file(0); file<file_count; ++file)
		{
			switch (va_arg(vlist, int))
			{
			case sf_normals_r:
		
				// Normals - reverse
				fullfilePath = filePath + "neqr.mtx";
				v_stageFileStreams_.push_back(fullfilePath);
				if (projectSettings_.a.recreate_stage_files)
					file_opener(f_normalsR_, fullfilePath, std::ios::out | std::ios::binary, binary);
				normalsR_map_.setnewFilePath(fullfilePath, projectSettings_.a.purge_stage_files);
			
				break;
			case sf_meas_minus_comp:

				// Measured minus computed
				fullfilePath = filePath + "mmc.mtx";
				v_stageFileStreams_.push_back(fullfilePath);
				if (projectSettings_.a.recreate_stage_files)
					file_opener(f_measMinusComp_, fullfilePath, std::ios::out | std::ios::binary, binary);
				measMinusComp_map_.setnewFilePath(fullfilePath, projectSettings_.a.purge_stage_files);
		
				break;
			case sf_estimated_stns:

				// Estimated stations
				fullfilePath = filePath + "est.mtx";
				v_stageFileStreams_.push_back(fullfilePath);
				if (projectSettings_.a.recreate_stage_files)
					file_opener(f_estimatedStations_, fullfilePath, std::ios::out | std::ios::binary, binary);
				estimatedStations_map_.setnewFilePath(fullfilePath, projectSettings_.a.purge_stage_files);

				break;
			case sf_original_stns:

				// Original stations
				fullfilePath = filePath + "ost.mtx";
				v_stageFileStreams_.push_back(fullfilePath);
				if (projectSettings_.a.recreate_stage_files)
					file_opener(f_originalStations_, fullfilePath, std::ios::out | std::ios::binary, binary);
				originalStations_map_.setnewFilePath(fullfilePath, projectSettings_.a.purge_stage_files);

				break;
			case sf_rigorous_stns:

				// Rigorous stations
				fullfilePath = filePath + "rst.mtx";
				v_stageFileStreams_.push_back(fullfilePath);
				if (projectSettings_.a.recreate_stage_files)
					file_opener(f_rigorousStations_, fullfilePath, std::ios::out | std::ios::binary, binary);
				rigorousStations_map_.setnewFilePath(fullfilePath, projectSettings_.a.purge_stage_files);

				break;
			case sf_junction_vars:

				// Junction variances
				fullfilePath = filePath + "jva.mtx";
				v_stageFileStreams_.push_back(fullfilePath);
				if (projectSettings_.a.recreate_stage_files)
					file_opener(f_junctionVariances_, fullfilePath, std::ios::out | std::ios::binary, binary);
				junctionVariances_map_.setnewFilePath(fullfilePath, projectSettings_.a.purge_stage_files);

				break;
			case sf_junction_vars_f:

				// Junction variances forward
				fullfilePath = filePath + "jvf.mtx";
				v_stageFileStreams_.push_back(fullfilePath);
				if (projectSettings_.a.recreate_stage_files)
					file_opener(f_junctionVariancesFwd_, fullfilePath, std::ios::out | std::ios::binary, binary);
				junctionVariancesFwd_map_.setnewFilePath(fullfilePath, projectSettings_.a.purge_stage_files);

				break;
			case sf_junction_ests_f:

				// Junction estimates forward
				fullfilePath = filePath + "jef.mtx";
				v_stageFileStreams_.push_back(fullfilePath);
				if (projectSettings_.a.recreate_stage_files)
					file_opener(f_junctionEstimatesFwd_, fullfilePath, std::ios::out | std::ios::binary, binary);
				junctionEstimatesFwd_map_.setnewFilePath(fullfilePath, projectSettings_.a.purge_stage_files);

				break;
			case sf_junction_ests_r:

				// Junction estimates reverse
				fullfilePath = filePath + "jer.mtx";
				v_stageFileStreams_.push_back(fullfilePath);
				if (projectSettings_.a.recreate_stage_files)
					file_opener(f_junctionEstimatesRev_, fullfilePath, std::ios::out | std::ios::binary, binary);
				junctionEstimatesRev_map_.setnewFilePath(fullfilePath, projectSettings_.a.purge_stage_files);

				break;
			case sf_rigorous_vars:

				// Rigorous variances
				fullfilePath = filePath + "rva.mtx";
				v_stageFileStreams_.push_back(fullfilePath);
				if (projectSettings_.a.recreate_stage_files)
					file_opener(f_rigorousVariances_, fullfilePath, std::ios::out | std::ios::binary, binary);
				rigorousVariances_map_.setnewFilePath(fullfilePath, projectSettings_.a.purge_stage_files);

				break;
			case sf_prec_adj_msrs:

				// Precision adjusted measurements
				fullfilePath = filePath + "pam.mtx";
				v_stageFileStreams_.push_back(fullfilePath);
				if (projectSettings_.a.recreate_stage_files)
					file_opener(f_precAdjMsrs_, fullfilePath, std::ios::out | std::ios::binary, binary);
				precAdjMsrs_map_.setnewFilePath(fullfilePath, projectSettings_.a.purge_stage_files);

				break;
			case sf_corrections:

				// Corrections
				fullfilePath = filePath + "cor.mtx";
				v_stageFileStreams_.push_back(fullfilePath);
				if (projectSettings_.a.recreate_stage_files)
					file_opener(f_corrections_, fullfilePath, std::ios::out | std::ios::binary, binary);
				corrections_map_.setnewFilePath(fullfilePath, projectSettings_.a.purge_stage_files);

				break;
			}
		}
		va_end(vlist);
	}
	catch (const std::runtime_error& e) {
		va_end(vlist);
		SignalExceptionAdjustment(e.what(), 0);
	}
}
	
void dna_adjust::SetMapRegions(const int file_count, ...)
{
	if (file_count == 0)
	{
		// serialise all.  That is, call this function again, but with 
		// arguments for all files
		SetMapRegions(12, 
			sf_normals_r, sf_meas_minus_comp,
			sf_estimated_stns, sf_original_stns, sf_rigorous_stns,
			sf_junction_vars, sf_junction_vars_f, sf_junction_ests_f, sf_junction_ests_r,
			sf_rigorous_vars, sf_prec_adj_msrs, sf_corrections);
		return;
	}

	va_list vlist;
	va_start(vlist, file_count);	

	for (UINT16 file(0); file<file_count; ++file)
	{
		switch (va_arg(vlist, int))
		{
		case sf_normals_r:

			// Create file map, and associate the regions with the normals blocks in the file
			normalsR_map_.CreateFileMap();		// File path set in OpenStageFileStreams
			for_each(
				normalsR_map_.vblockMapRegions_.begin(), 
				normalsR_map_.vblockMapRegions_.end(), 
				[this] (block_map_t& block_map) {			// use lambda expression
					block_map.MapRegion(normalsR_map_.getFileMapPtr());
				}
			);

			break;
		case sf_meas_minus_comp:

			// Create file map, and associate the regions with the normals blocks in the file
			measMinusComp_map_.CreateFileMap();		// File path set in OpenStageFileStreams
			for_each(
				measMinusComp_map_.vblockMapRegions_.begin(), 
				measMinusComp_map_.vblockMapRegions_.end(), 
				[this] (block_map_t& block_map) {			// use lambda expression
					block_map.MapRegion(measMinusComp_map_.getFileMapPtr());
				}
			);

			break;
		case sf_estimated_stns:

			// Create file map, and associate the regions with the normals blocks in the file
			estimatedStations_map_.CreateFileMap();		// File path set in OpenStageFileStreams
			for_each(
				estimatedStations_map_.vblockMapRegions_.begin(), 
				estimatedStations_map_.vblockMapRegions_.end(), 
				[this] (block_map_t& block_map) {			// use lambda expression
					block_map.MapRegion(estimatedStations_map_.getFileMapPtr());
				}
			);

			break;
		case sf_original_stns:

			// Create file map, and associate the regions with the normals blocks in the file
			originalStations_map_.CreateFileMap();		// File path set in OpenStageFileStreams
			for_each(
				originalStations_map_.vblockMapRegions_.begin(), 
				originalStations_map_.vblockMapRegions_.end(), 
				[this] (block_map_t& block_map) {			// use lambda expression
					block_map.MapRegion(originalStations_map_.getFileMapPtr());
				}
			);

			break;
		case sf_rigorous_stns:

			// Create file map, and associate the regions with the normals blocks in the file
			rigorousStations_map_.CreateFileMap();		// File path set in OpenStageFileStreams
			for_each(
				rigorousStations_map_.vblockMapRegions_.begin(), 
				rigorousStations_map_.vblockMapRegions_.end(), 
				[this] (block_map_t& block_map) {			// use lambda expression
					block_map.MapRegion(rigorousStations_map_.getFileMapPtr());
				}
			);

			break;
		case sf_junction_vars:

			// Create file map, and associate the regions with the normals blocks in the file
			junctionVariances_map_.CreateFileMap();		// File path set in OpenStageFileStreams
			for_each(
				junctionVariances_map_.vblockMapRegions_.begin(), 
				junctionVariances_map_.vblockMapRegions_.end(), 
				[this] (block_map_t& block_map) {			// use lambda expression
					block_map.MapRegion(junctionVariances_map_.getFileMapPtr());
				}
			);

			break;
		case sf_junction_vars_f:

			// Create file map, and associate the regions with the normals blocks in the file
			junctionVariancesFwd_map_.CreateFileMap();		// File path set in OpenStageFileStreams
			for_each(
				junctionVariancesFwd_map_.vblockMapRegions_.begin(), 
				junctionVariancesFwd_map_.vblockMapRegions_.end(), 
				[this] (block_map_t& block_map) {			// use lambda expression
					block_map.MapRegion(junctionVariancesFwd_map_.getFileMapPtr());
				}
			);

			break;
		case sf_junction_ests_f:

			// Create file map, and associate the regions with the normals blocks in the file
			junctionEstimatesFwd_map_.CreateFileMap();		// File path set in OpenStageFileStreams
			for_each(
				junctionEstimatesFwd_map_.vblockMapRegions_.begin(), 
				junctionEstimatesFwd_map_.vblockMapRegions_.end(), 
				[this] (block_map_t& block_map) {			// use lambda expression
					block_map.MapRegion(junctionEstimatesFwd_map_.getFileMapPtr());
				}
			);

			break;
		case sf_junction_ests_r:

			// Create file map, and associate the regions with the normals blocks in the file
			junctionEstimatesRev_map_.CreateFileMap();		// File path set in OpenStageFileStreams
			for_each(
				junctionEstimatesRev_map_.vblockMapRegions_.begin(), 
				junctionEstimatesRev_map_.vblockMapRegions_.end(), 
				[this] (block_map_t& block_map) {			// use lambda expression
					block_map.MapRegion(junctionEstimatesRev_map_.getFileMapPtr());
				}
			);

			break;
		case sf_rigorous_vars:

			// Create file map, and associate the regions with the normals blocks in the file
			rigorousVariances_map_.CreateFileMap();		// File path set in OpenStageFileStreams
			for_each(
				rigorousVariances_map_.vblockMapRegions_.begin(), 
				rigorousVariances_map_.vblockMapRegions_.end(), 
				[this] (block_map_t& block_map) {			// use lambda expression
					block_map.MapRegion(rigorousVariances_map_.getFileMapPtr());
				}
			);

			break;
		case sf_prec_adj_msrs:

			// Create file map, and associate the regions with the normals blocks in the file
			precAdjMsrs_map_.CreateFileMap();		// File path set in OpenStageFileStreams
			for_each(
				precAdjMsrs_map_.vblockMapRegions_.begin(), 
				precAdjMsrs_map_.vblockMapRegions_.end(), 
				[this] (block_map_t& block_map) {			// use lambda expression
					block_map.MapRegion(precAdjMsrs_map_.getFileMapPtr());
				}
			);

			break;
		case sf_corrections:

			// Create file map, and associate the regions with the normals blocks in the file
			corrections_map_.CreateFileMap();		// File path set in OpenStageFileStreams
			for_each(
				corrections_map_.vblockMapRegions_.begin(), 
				corrections_map_.vblockMapRegions_.end(), 
				[this] (block_map_t& block_map) {			// use lambda expression
					block_map.MapRegion(corrections_map_.getFileMapPtr());
				}
			);

			break;
		}
	}
	va_end(vlist);
}
	

void dna_adjust::PrepareMemoryMapRegions(const UINT32& block, const int file_count, ...)
{
	if (file_count == 0)
	{
		// serialise all.  That is, call this function again, but with 
		// arguments for all files
		PrepareMemoryMapRegions(block, 12, 
			sf_normals_r, sf_meas_minus_comp,
			sf_estimated_stns, sf_original_stns, sf_rigorous_stns,
			sf_junction_vars, sf_junction_vars_f, sf_junction_ests_f, sf_junction_ests_r,
			sf_rigorous_vars, sf_prec_adj_msrs, sf_corrections);
		return;
	}

	va_list vlist;
	va_start(vlist, file_count);	

	for (UINT16 file(0); file<file_count; ++file)
	{
		switch (va_arg(vlist, int))
		{
		case sf_normals_r:

			// Create file map, and associate the regions with the normals blocks in the file
			normalsR_map_.CreateFileMap();				// File path set in OpenStageFileStreams
			normalsR_map_.MapRegion(block);

			break;
		case sf_meas_minus_comp:

			// Create file map, and associate the regions with the normals blocks in the file
			measMinusComp_map_.CreateFileMap();			// File path set in OpenStageFileStreams
			measMinusComp_map_.MapRegion(block);
		
			break;
		case sf_estimated_stns:

			// Create file map, and associate the regions with the normals blocks in the file
			estimatedStations_map_.CreateFileMap();		// File path set in OpenStageFileStreams
			estimatedStations_map_.MapRegion(block);
		
			break;
		case sf_original_stns:

			// Create file map, and associate the regions with the normals blocks in the file
			originalStations_map_.CreateFileMap();		// File path set in OpenStageFileStreams
			originalStations_map_.MapRegion(block);
		
			break;
		case sf_rigorous_stns:

			// Create file map, and associate the regions with the normals blocks in the file
			rigorousStations_map_.CreateFileMap();		// File path set in OpenStageFileStreams
			rigorousStations_map_.MapRegion(block);
		
			break;
		case sf_junction_vars:

			// Create file map, and associate the regions with the normals blocks in the file
			junctionVariances_map_.CreateFileMap();		// File path set in OpenStageFileStreams
			junctionVariances_map_.MapRegion(block);
		
			break;
		case sf_junction_vars_f:

			// Create file map, and associate the regions with the normals blocks in the file
			junctionVariancesFwd_map_.CreateFileMap();	// File path set in OpenStageFileStreams
			junctionVariancesFwd_map_.MapRegion(block);
		
			break;
		case sf_junction_ests_f:

			// Create file map, and associate the regions with the normals blocks in the file
			junctionEstimatesFwd_map_.CreateFileMap();	// File path set in OpenStageFileStreams
			junctionEstimatesFwd_map_.MapRegion(block);
		
			break;
		case sf_junction_ests_r:

			// Create file map, and associate the regions with the normals blocks in the file
			junctionEstimatesRev_map_.CreateFileMap();	// File path set in OpenStageFileStreams
			junctionEstimatesRev_map_.MapRegion(block);
		
			break;
		case sf_rigorous_vars:

			// Create file map, and associate the regions with the normals blocks in the file
			rigorousVariances_map_.CreateFileMap();		// File path set in OpenStageFileStreams
			rigorousVariances_map_.MapRegion(block);
		
			break;
		case sf_prec_adj_msrs:

			// Create file map, and associate the regions with the normals blocks in the file
			precAdjMsrs_map_.CreateFileMap();			// File path set in OpenStageFileStreams
			precAdjMsrs_map_.MapRegion(block);
		
			break;
		case sf_corrections:

			// Create file map, and associate the regions with the normals blocks in the file
			corrections_map_.CreateFileMap();			// File path set in OpenStageFileStreams
			corrections_map_.MapRegion(block);

			break;
		}
	}
	va_end(vlist);
}
	

void dna_adjust::OffloadBlockToDisk(const UINT32& block)
{
	// Write block matrix data to disk
	SerialiseBlockToDisk(block);

	// Unload block matrix data from memory
	UnloadBlock(block);
}
	

void dna_adjust::OffloadBlockToMappedFile(const UINT32& block)
{
	// Write block matrix data to disk
	SerialiseBlockToMappedFile(block);

	// Unload block matrix data from memory
	UnloadBlock(block);

	// Advise kernel that mapped pages for this block are no longer needed,
	// freeing page cache for upcoming blocks
	AdviseBlockDontNeed(block);
}


void dna_adjust::AdviseBlockDontNeed(const UINT32& block)
{
	using advice = boost::interprocess::mapped_region::advice_types;
	normalsR_map_.AdviseRegion(block, advice::advice_dontneed);
	measMinusComp_map_.AdviseRegion(block, advice::advice_dontneed);
	estimatedStations_map_.AdviseRegion(block, advice::advice_dontneed);
	originalStations_map_.AdviseRegion(block, advice::advice_dontneed);
	rigorousStations_map_.AdviseRegion(block, advice::advice_dontneed);
	junctionVariances_map_.AdviseRegion(block, advice::advice_dontneed);
	junctionVariancesFwd_map_.AdviseRegion(block, advice::advice_dontneed);
	junctionEstimatesFwd_map_.AdviseRegion(block, advice::advice_dontneed);
	junctionEstimatesRev_map_.AdviseRegion(block, advice::advice_dontneed);
	rigorousVariances_map_.AdviseRegion(block, advice::advice_dontneed);
	precAdjMsrs_map_.AdviseRegion(block, advice::advice_dontneed);
	corrections_map_.AdviseRegion(block, advice::advice_dontneed);
}


void dna_adjust::AdviseBlockWillNeed(const UINT32& block)
{
	using advice = boost::interprocess::mapped_region::advice_types;
	normalsR_map_.AdviseRegion(block, advice::advice_willneed);
	measMinusComp_map_.AdviseRegion(block, advice::advice_willneed);
	estimatedStations_map_.AdviseRegion(block, advice::advice_willneed);
	originalStations_map_.AdviseRegion(block, advice::advice_willneed);
	rigorousStations_map_.AdviseRegion(block, advice::advice_willneed);
	junctionVariances_map_.AdviseRegion(block, advice::advice_willneed);
	junctionVariancesFwd_map_.AdviseRegion(block, advice::advice_willneed);
	junctionEstimatesFwd_map_.AdviseRegion(block, advice::advice_willneed);
	junctionEstimatesRev_map_.AdviseRegion(block, advice::advice_willneed);
	rigorousVariances_map_.AdviseRegion(block, advice::advice_willneed);
	precAdjMsrs_map_.AdviseRegion(block, advice::advice_willneed);
	corrections_map_.AdviseRegion(block, advice::advice_willneed);
}
	

void dna_adjust::SerialiseBlockToDisk(const UINT32& block)
{
	// Write block matrix data to disk
	try {
		// Normals - reverse
		f_normalsR_ << v_normalsR_.at(block);							// ...neqr.mtx		
		// Measured minus computed
		f_measMinusComp_ << v_measMinusComp_.at(block);					// ...mmc.mtx		
		// Estimated stations
		f_estimatedStations_ << v_estimatedStations_.at(block);			// ...est.mtx
		// original stations
		f_originalStations_ << v_originalStations_.at(block);			// ...ost.mtx
		// Rigorous stations
		f_rigorousStations_ << v_rigorousStations_.at(block);			// ...rst.mtx
		// Junction variances
		f_junctionVariances_ << v_junctionVariances_.at(block);			// ...jva.mtx
		// Junction variances forward
		f_junctionVariancesFwd_ << v_junctionVariancesFwd_.at(block);	// ...jvf.mtx
		// Junction estimates forward
		f_junctionEstimatesFwd_ << v_junctionEstimatesFwd_.at(block);	// ...jef.mtx
		// Junction estimates reverse
		f_junctionEstimatesRev_ << v_junctionEstimatesRev_.at(block);	// ...jer.mtx
		// Rigorous variances
		f_rigorousVariances_ << v_rigorousVariances_.at(block);			// ...rva.mtx
		// Precision adjusted measurements
		f_precAdjMsrs_ << v_precAdjMsrsFull_.at(block);					// ...pam.mtx
		// Corrections
		f_corrections_ << v_corrections_.at(block);						// ...cor.mtx
	}
	catch (...) {
		std::stringstream ss;
		ss << "SerialiseBlockToDisk(): An error was encountered when writing matrix data to disk." << std::endl;
		SignalExceptionAdjustment(ss.str(), 0);
	}	
}
	

void dna_adjust::UnloadBlock(const UINT32& block, const int file_count, ...)
{
	if (file_count == 0)
	{
		// deserialise all
		UnloadBlock(block, 15,
			sf_normals, sf_normals_r, sf_atvinv, sf_design, sf_meas_minus_comp,
			sf_estimated_stns, sf_original_stns, sf_rigorous_stns,
			sf_junction_vars, sf_junction_vars_f, sf_junction_ests_f, sf_junction_ests_r,
			sf_rigorous_vars, sf_prec_adj_msrs, sf_corrections);
		return;
	}

	va_list vlist;
	va_start(vlist, file_count);

	// Unload block matrix data from memory by freeing buffers.
	// Use deallocate() instead of explicit destructor calls to keep
	// the matrix objects in a valid state for later reuse.
	for (UINT16 file(0); file<file_count; ++file)
	{
		switch (va_arg(vlist, int))
		{
		case sf_normals:
			v_normals_.at(block).deallocate();
			break;
		case sf_normals_r:
			v_normalsR_.at(block).deallocate();
			break;
		case sf_atvinv:
			v_AtVinv_.at(block).deallocate();
			break;
		case sf_design:
			v_design_.at(block).deallocate();
			break;
		case sf_meas_minus_comp:
			v_measMinusComp_.at(block).deallocate();
			break;
		case sf_estimated_stns:
			v_estimatedStations_.at(block).deallocate();
			break;
		case sf_original_stns:
			v_originalStations_.at(block).deallocate();
			break;
		case sf_rigorous_stns:
			v_rigorousStations_.at(block).deallocate();
			break;
		case sf_junction_vars:
			v_junctionVariances_.at(block).deallocate();
			break;
		case sf_junction_vars_f:
			v_junctionVariancesFwd_.at(block).deallocate();
			break;
		case sf_junction_ests_f:
			v_junctionEstimatesFwd_.at(block).deallocate();
			break;
		case sf_junction_ests_r:
			v_junctionEstimatesRev_.at(block).deallocate();
			break;
		case sf_rigorous_vars:
			v_rigorousVariances_.at(block).deallocate();
			break;
		case sf_prec_adj_msrs:
			v_precAdjMsrsFull_.at(block).deallocate();
			break;
		case sf_corrections:
			v_corrections_.at(block).deallocate();

			if (v_blockMeta_.at(block)._blockLast)
				v_correctionsR_.at(block).deallocate();
			break;
		}
	}
	va_end(vlist);
}
	

void dna_adjust::CloseStageFileStreams()
{
	// Normals -- reverse
	f_normalsR_.close();

	// Measured minus computed
	f_measMinusComp_.close();

	// Estimated stations
	f_estimatedStations_.close();

	// Original stations
	f_originalStations_.close();

	// Rigorous stations
	f_rigorousStations_.close();

	// Junction variances
	f_junctionVariances_.close();

	// Junction variances forward
	f_junctionVariancesFwd_.close();

	// Junction estimates forward
	f_junctionEstimatesFwd_.close();

	// Junction estimates reverse
	f_junctionEstimatesRev_.close();

	// Rigorous variances
	f_rigorousVariances_.close();

	// Precision adjusted measurements
	f_precAdjMsrs_.close();

	// Corrections
	f_corrections_.close();

}


}	// namespace networkadjust
}	// namespace dynadjust

