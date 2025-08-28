//============================================================================
// Name         : dnafilepathfuncs.hpp
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
// Description  : File Path Manipulation Functions
//============================================================================

#ifndef DNAFILEPATHFUNCS_H_
#define DNAFILEPATHFUNCS_H_

#if defined(_MSC_VER)
	#if defined(LIST_INCLUDES_ON_BUILD) 
		#pragma message("  " __FILE__) 
	#endif
#endif

/// \cond
#include <stdio.h>
#include <stdarg.h>
#include <string>
#include <string_view>

#include <algorithm>
#include <functional>

#include <filesystem>
/// \endcond

#include <include/config/dnaconsts.hpp>

template <typename T>
T formPath(const T& folder, const T& file, const T& ext)
{
	return T(folder + FOLDER_SLASH + file + "." + ext);
}

template <typename T>
T formPath(const T& folder, const T file)
{
	return T(folder + FOLDER_SLASH + file);
}

template <typename T>
T leafStr(const T& filePath)
{
	return std::filesystem::path(filePath).filename().string();
}

// Safe wrapper for std::filesystem::absolute that handles empty paths
// On Linux, std::filesystem::absolute("") throws an exception
// This function returns empty string for empty input
// Accepts std::string_view for flexibility with both std::string and string_view
inline std::string safe_absolute_path(std::string_view path)
{
	if (path.empty())
		return "";
	return std::filesystem::absolute(std::string(path)).string();
}

#endif //DNAFILEPATHFUNCS_H_
