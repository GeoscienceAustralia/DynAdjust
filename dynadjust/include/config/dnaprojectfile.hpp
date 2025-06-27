//============================================================================
// Name         : dnaprojectfile.hpp
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
// Description  : DynAdjust Project class
//============================================================================

#ifndef DNAPROJECTFILE_H_
#define DNAPROJECTFILE_H_

#if defined(_MSC_VER)
	#if defined(LIST_INCLUDES_ON_BUILD) 
		#pragma message("  " __FILE__) 
	#endif
#endif

/// \cond
#include <string>
#include <sstream>
#include <iomanip>

#include <boost/timer/timer.hpp>
/// \endcond

#include <include/config/dnatypes.hpp>
#include <include/config/dnaoptions-interface.hpp>
#include <include/config/dnaoptions.hpp>
#include <include/functions/dnastrmanipfuncs.hpp>

namespace dynadjust {

class CDnaProjectFile
{
public:

	CDnaProjectFile(void);
	CDnaProjectFile(const std::string& projectFile, const UINT16& verifyContext);
	
	virtual inline ~CDnaProjectFile(void) {}
	//virtual inline CDnaProjectFile* clone() const { 
	//	return new CDnaProjectFile(*this); 
	//}

private:
	// Disallow
	CDnaProjectFile(const project_settings& project);
	CDnaProjectFile(const CDnaProjectFile& newProject);
	CDnaProjectFile& operator=(const CDnaProjectFile& rhs);
	bool operator==(const CDnaProjectFile& rhs) const;
	
public:
	inline project_settings GetSettings() const { return settings_; }
	//inline CDnaProjectFile& operator[](int iIndex) { return this[iIndex]; }

	void LoadProjectFile(const std::string& projectFile);
	void LoadProjectFile();
	void PrintProjectFile();

	template <typename T, typename U>
	void PrintRecord(std::ostream& os, const T& variable, const U& value);
	template <typename T>
	void PrintVariable(std::ostream& os, const T& variable);
	template <typename T>
	void PrintValue(std::ostream& os, const T& value);

	inline void UpdateSettingsGeneral(const general_settings& g) { 
		settings_.g = g;
	}
	inline void UpdateSettingsImport(const project_settings& p) { 
		UpdateSettingsGeneral(p.g);
		settings_.i = p.i;
	}
	inline void UpdateSettingsGeoid(const project_settings& p) { 
		UpdateSettingsGeneral(p.g);
		settings_.n = p.n; 
	}
	inline void UpdateSettingsReftran(const project_settings& p) { 
		UpdateSettingsGeneral(p.g);
		settings_.r = p.r;
	}
	inline void UpdateSettingsAdjust(const project_settings& p) { 
		UpdateSettingsGeneral(p.g);
		settings_.a = p.a;
	}
	inline void UpdateSettingsSegment(const project_settings& p) { 
		UpdateSettingsGeneral(p.g);
		settings_.s = p.s;
	}
	inline void UpdateSettingsPlot(const project_settings& p) { 
		UpdateSettingsGeneral(p.g);
		settings_.p = p.p; 
	}
	inline void UpdateSettingsOutput(const project_settings& p) { 
		UpdateSettingsGeneral(p.g);
		settings_.o = p.o; 
	}

	template <typename T>
	void AddOptionValue(std::ostream& os, const char* const option, const T& value);
	template <typename T>
	void AddDefaultValue(std::ostream& os, const T& value);

	std::string FormCommandLineOptionsStringImport();
	std::string FormCommandLineOptionsStringGeneral();

private:
	void InitialiseGeneralSettings();
	void InitialiseImportSettings();
	void InitialiseReftranSettings();
	void InitialiseGeoidSettings();
	void InitialiseSegmentSettings();
	void InitialiseAdjustSettings();

	void LoadSettingGeneral(const settingMode mSetting, const std::string& var, std::string& val);
	void LoadSettingImport(const settingMode mSetting, const std::string& var, std::string& val);
	void LoadSettingReftran(const std::string& var, std::string& val);
	void LoadSettingGeoid(const std::string& var, std::string& val);
	void LoadSettingSegment(const std::string& var, std::string& val);
	void LoadSettingAdjust(const std::string& var, std::string& val);
	void LoadSettingOutput(const std::string& var, std::string& val);
	//void LoadSettingPlot(const std::string& var, std::string& val);
	//void LoadSettingDisplay(const std::string& var, std::string& val);
	
	project_settings settings_;
};


}	// namespace dynadjust

#endif /* DNADATUMPROJECTION_H_ */
