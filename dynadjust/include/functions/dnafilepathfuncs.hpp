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

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string>

#include <algorithm>
#include <functional>

#include <boost/filesystem.hpp>
#include <boost/iostreams/detail/absolute_path.hpp>

#include <include/config/dnaconsts.hpp>

using namespace std;
using namespace boost;
using namespace boost::filesystem;

using namespace std;

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
	return path(filePath).filename().string();
}

#endif //DNAFILEPATHFUNCS_H_