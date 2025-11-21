//============================================================================
// Name         : dnaprocessfuncs.cpp
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
// Description  : Common process and fork functions
//============================================================================

#include <include/functions/dnaprocessfuncs.hpp>

#include <cstdlib>

bool run_command(const string& executable_path, const UINT16& quiet)
{	
	// use boost's platform independent code to invoke a process
	// see https://www.boost.org/doc/libs/develop/doc/html/process.html
	//
	// An exit code of zero means the process was successful, 
	// while one different than zero indicates an error.

	
	// For windows batch files, add cmd to execute the batch file.
#if defined(_WIN32) || defined(__WIN32__)

	STARTUPINFO startup;
	PROCESS_INFORMATION process;

	memset(&startup, 0, sizeof(STARTUPINFO));
	memset(&process, 0, sizeof(PROCESS_INFORMATION));

	startup.cb = sizeof(STARTUPINFO);

	DWORD return_code(EXIT_SUCCESS);
	if (CreateProcess(0, (LPSTR)executable_path.c_str(), 0, 0, FALSE,
		0, 0, 0, &startup, &process))
	{
		WaitForSingleObject(process.hProcess, INFINITE);
		// Capture the return code
		GetExitCodeProcess(process.hProcess, &return_code);
		CloseHandle(process.hThread);
		CloseHandle(process.hProcess);

		return (return_code == EXIT_SUCCESS);
	}
	return EXIT_SUCCESS;

#elif defined(__linux) || defined(sun) || defined(__unix__) || defined(__APPLE__)

	int return_value(0);

	if (quiet) {
		// Redirect stdout and stderr to /dev/null
		std::string cmd = executable_path + " > /dev/null 2>&1";
		return_value = std::system(cmd.c_str());
	} else {
		return_value = std::system(executable_path.c_str());
	}

	return (return_value == EXIT_SUCCESS);

#endif

	return EXIT_SUCCESS;
}


