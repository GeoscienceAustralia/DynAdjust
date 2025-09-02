//============================================================================
// Name         : dnafile_mapping.hpp
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
// Description  : DynAdjust Memory Mapped File Library
//============================================================================

#ifndef DNAFILEMAPPING_H_
#define DNAFILEMAPPING_H_

#if defined(_MSC_VER)
#if defined(LIST_INCLUDES_ON_BUILD) 
#pragma message("  " __FILE__) 
#endif
#endif

/// \cond
#include <memory>
#include <filesystem>

#include <boost/interprocess/file_mapping.hpp>
#include <boost/interprocess/mapped_region.hpp>
/// \endcond

#include <include/config/dnatypes.hpp>

namespace dynadjust { 
namespace memory {

typedef std::shared_ptr<boost::interprocess::mapped_region> MapRegPtr;
typedef std::shared_ptr<boost::interprocess::file_mapping> FileMapPtr;

class block_map_t
{
public:
	block_map_t(const size_t& size);
	block_map_t(const block_map_t &p);

private:
	// Disallow use of compiler generated functions. 
	block_map_t();
	// virtual inline block_map_t* clone() const { return new block_map_t(*this); }
	block_map_t& operator=(const block_map_t& rhs);
	bool operator==(const block_map_t& rhs) const;

public:

	inline void* GetRegionAddr() { return region_ptr_->get_address(); }
	inline size_t GetDataSize() const { return data_size_; }
	inline size_t GetRegionOffset() const { return region_offset_; }
	inline size_t GetCumulativeRegionOffset() const { return region_offset_ + data_size_; }
	
	inline void SetDataSize(const size_t& size) { data_size_ = size; }
	inline void SetRegionOffset(const size_t& size) { region_offset_ = size; }
	
	void MapRegion(FileMapPtr file_map_ptr);

	size_t			data_size_;		// Size of this matrix.  
	size_t			region_offset_;		// Offset from the beginning of the region
	MapRegPtr		region_ptr_;		// shared pointer to the region
};

typedef std::vector<block_map_t> vblock_map;
typedef vblock_map::iterator _it_block_map;

// class to hold addresses and sizes for all matrices 
// in a vector of segmented blocks
class vmat_file_map {
public:
	vmat_file_map();
	vmat_file_map(const std::string& filePath, bool remove_mapped_file);
	~vmat_file_map();

	void reserveblockMapRegions(const UINT32& size);
	void addblockMapRegion(const block_map_t& map);

	void setnewFilePath(const std::string& filePath, bool remove_mapped_file);
	void CreateFileMap();
	void MapRegion(const UINT32 block);

	inline FileMapPtr getFileMapPtr() const { return file_map_ptr_; }
	inline void* GetBlockRegionAddr(const UINT32 block) const { 
		return vblockMapRegions_.at(block).region_ptr_->get_address(); 
	}

	vmat_file_map(const vmat_file_map&);				// prevent copying
	vmat_file_map& operator=(const vmat_file_map&);		//   ''      ''

	vblock_map		vblockMapRegions_;
	std::string			filePath_;
	FileMapPtr		file_map_ptr_;
	bool			remove_mapped_file_;
};

typedef std::vector<vmat_file_map> vvmatf_map;
typedef vvmatf_map::iterator _it_vvmatf_map;


}	// namespace memory 
}	// namespace dynadjust 

#endif	// DNAFILEMAPPING_H_

