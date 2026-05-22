//============================================================================
// Name         : dnafile_mapping.cpp
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

#include <include/memory/dnafile_mapping.hpp>

#ifdef __linux__
#include <cstdint>
#include <sys/mman.h>
#include <unistd.h>
#endif

namespace dynadjust {
namespace memory {

//block_map_t::block_map_t()
//		: data_size_(0), region_offset_(0) 
//{
//
//}
	

block_map_t::block_map_t(const size_t& size)
		: data_size_(size), region_offset_(0) 
{

}
	

block_map_t::block_map_t(const block_map_t &p) 
		: data_size_(p.data_size_)
		, region_offset_(p.region_offset_)
		, region_ptr_(p.region_ptr_) 
{

}
	

//block_map_t& block_map_t::operator=(const block_map_t& rhs) 
//{
//	if (this == &rhs)
//		return *this;
//	data_size_ = rhs.data_size_;
//	region_offset_ = rhs.region_offset_;
//	region_ptr_ = rhs.region_ptr_;
//	return *this;
//}
	

//bool block_map_t::operator==(const block_map_t& rhs) const 
//{
//	return (
//		data_size_ == rhs.data_size_ &&
//		region_offset_ == rhs.region_offset_ &&
//		region_ptr_ == rhs.region_ptr_
//		);
//}
	

void block_map_t::MapRegion(FileMapPtr file_map_ptr) {
	region_ptr_.reset(
		new boost::interprocess::mapped_region(
			*file_map_ptr,
			boost::interprocess::read_write,
			region_offset_,
			data_size_
			)
		);
}


void block_map_t::AdviseSequential() {
	if (region_ptr_)
		region_ptr_->advise(boost::interprocess::mapped_region::advice_sequential);
}


void block_map_t::AdviseDontNeed() {
	if (region_ptr_)
		region_ptr_->advise(boost::interprocess::mapped_region::advice_dontneed);
}


void block_map_t::AdviseWillNeed() {
	if (region_ptr_)
		region_ptr_->advise(boost::interprocess::mapped_region::advice_willneed);
}


// class to hold addresses and sizes for all matrices 
// in a vector of segmented blocks
vmat_file_map::vmat_file_map()
{

}
	

vmat_file_map::vmat_file_map(const std::string& filePath, bool remove_mapped_file) 
	: filePath_(filePath), remove_mapped_file_(remove_mapped_file) 
{ 
		file_map_ptr_.reset(new boost::interprocess::file_mapping(filePath_.c_str(), boost::interprocess::read_write));
}
	

vmat_file_map::~vmat_file_map() 
{
	if (remove_mapped_file_)
		if (std::filesystem::exists(filePath_))
			boost::interprocess::file_mapping::remove(filePath_.c_str());
}
	

void vmat_file_map::reserveblockMapRegions(const UINT32& size)
{
	vblockMapRegions_.reserve(size); 
}
	

void vmat_file_map::addblockMapRegion(const block_map_t& map) 
{
	vblockMapRegions_.push_back(map); 
}


void vmat_file_map::setnewFilePath(const std::string& filePath, bool remove_mapped_file) 
{
	filePath_ = filePath; 
	remove_mapped_file_ = remove_mapped_file;
}
	

void vmat_file_map::CreateFileMap() 
{
	file_map_ptr_.reset(new boost::interprocess::file_mapping(filePath_.c_str(), boost::interprocess::read_write));
}
	

void vmat_file_map::MapRegion(const UINT32 block)
{
	vblockMapRegions_.at(block).MapRegion(file_map_ptr_);
}


void vmat_file_map::AdviseRegion(const UINT32 block, boost::interprocess::mapped_region::advice_types advice)
{
	auto& region = vblockMapRegions_.at(block);
#ifdef __linux__
	// Boost's mapped_region::advise calls madvise() without
	// page-aligning the address.  Stage regions are packed
	// end-to-end so only the first region's address is page-
	// aligned by chance; madvise on the rest fails silently
	// with EINVAL.  Round start up and end down to whole
	// pages, then call madvise directly.  For DONTNEED also
	// flush dirty pages asynchronously so the kernel can
	// evict them promptly.
	if (region.region_ptr_ && region.data_size_ > 0) {
		const std::size_t page_size =
			static_cast<std::size_t>(sysconf(_SC_PAGESIZE));
		auto base = reinterpret_cast<std::uintptr_t>(
			region.region_ptr_->get_address());
		std::uintptr_t aligned_start =
			(base + page_size - 1) & ~(page_size - 1);
		std::uintptr_t aligned_end =
			(base + region.data_size_) & ~(page_size - 1);
		if (aligned_end > aligned_start) {
			void* addr = reinterpret_cast<void*>(aligned_start);
			std::size_t len = aligned_end - aligned_start;
			int madv = MADV_NORMAL;
			switch (advice) {
			case boost::interprocess::mapped_region::advice_sequential:
				madv = MADV_SEQUENTIAL; break;
			case boost::interprocess::mapped_region::advice_willneed:
				madv = MADV_WILLNEED; break;
			case boost::interprocess::mapped_region::advice_dontneed:
				msync(addr, len, MS_ASYNC);
				madv = MADV_DONTNEED; break;
			default: break;
			}
			madvise(addr, len, madv);
			return;
		}
	}
#endif
	region.region_ptr_->advise(advice);
}
	

}	// namespace memory 
}	// namespace dynadjust 


