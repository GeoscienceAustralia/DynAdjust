//============================================================================
// Name         : dnaadjust.hpp
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

#ifndef DNAADJUST_H_
#define DNAADJUST_H_

#if defined(_MSC_VER)
#if defined(LIST_INCLUDES_ON_BUILD)
#pragma message("  " __FILE__)
#endif
#endif

/// \cond
#include <cstdarg>
#include <exception>
#include <fstream>
#include <iostream>
#include <math.h>
#include <queue>
#include <stdexcept>
#include <string>
#include <utility>
#include <memory>
#include <vector>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <chrono>
#include <functional>
#include <atomic>
#include <filesystem>

#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/normal.hpp>

/// \endcond

#include <include/config/dnaversion.hpp>
#include <include/config/dnaconsts.hpp>
#include <include/config/dnaexports.hpp>
#include <include/config/dnaoptions-interface.hpp>
#include <include/config/dnatypes-gui.hpp>
#include <include/config/dnatypes.hpp>
#include <include/exception/dnaexception.hpp>

#include <include/io/adj_file.hpp>
#include <include/io/aml_file.hpp>
#include <include/io/asl_file.hpp>
#include <include/io/bms_file.hpp>
#include <include/io/bst_file.hpp>
#include <include/io/map_file.hpp>
#include <include/io/seg_file.hpp>
#include <include/io/dnaiosnx.hpp>
#include <include/io/dnaiotbu.hpp>

#include <include/functions/dnafilepathfuncs.hpp>
#include <include/functions/dnaintegermanipfuncs.hpp>
#include <include/functions/dnaiostreamfuncs.hpp>
#include <include/functions/dnastringfuncs.hpp>
#include <include/functions/dnatemplatematrixfuncs.hpp>
#include <include/functions/dnatemplatestnmsrfuncs.hpp>
#include <include/functions/dnatimer.hpp>

#include <include/math/dnamatrix_contiguous.hpp>
#include <include/memory/dnafile_mapping.hpp>
#include <include/parameters/dnadatum.hpp>
#include <include/parameters/dnaepsg.hpp>
#include <include/parameters/dnaprojection.hpp>
#include <include/thread/dnathreading.hpp>

#include <dynadjust/dnaadjust/network_data_loader.hpp>
#include <dynadjust/dnaadjust/dnaadjust_printer.hpp>

using namespace dynadjust::datum_parameters;
using namespace dynadjust::measurements;
using namespace dynadjust::math;
using namespace dynadjust::exception;
using namespace dynadjust::memory;
using namespace dynadjust::iostreams;

namespace dynadjust {
namespace networkadjust {

extern std::mutex maxCorrMutex;

using dynadjust::cpu_timer;

// forward declaration of dna_adjust
class dna_adjust;
class DynAdjustPrinter;

class adjust_prepare_thread {
  public:
    adjust_prepare_thread(dna_adjust* dnaAdj, std::exception_ptr& error)
        : main_adj_(dnaAdj), error_(error) {}
    void operator()();

  private:
    dna_adjust* main_adj_;
    std::exception_ptr& error_;

    // Prevent assignment operator
    adjust_prepare_thread& operator=(const adjust_prepare_thread& rhs);
};

class adjust_process_prepare_thread {
  public:
    adjust_process_prepare_thread(dna_adjust* dnaAdj, const UINT32& id,
                                  std::exception_ptr& error,
                                  std::vector<std::exception_ptr>& prep_errors)
        : main_adj_(dnaAdj), thread_id_(id), error_(error),
          prep_errors_(prep_errors) {}
    void operator()();

  private:
    dna_adjust* main_adj_;
    UINT32 thread_id_;
    std::exception_ptr& error_;
    std::vector<std::exception_ptr>& prep_errors_;

    // Prevent assignment operator
    adjust_process_prepare_thread&
    operator=(const adjust_process_prepare_thread& rhs);
};

class adjust_forward_thread {
  public:
    adjust_forward_thread(dna_adjust* dnaAdj, std::exception_ptr& error)
        : main_adj_(dnaAdj), error_(error) {}
    void operator()();

  private:
    dna_adjust* main_adj_;
    std::exception_ptr& error_;

    // Prevent assignment operator
    adjust_forward_thread& operator=(const adjust_forward_thread& rhs);
};

class adjust_reverse_thread {
  public:
    adjust_reverse_thread(dna_adjust* dnaAdj, std::exception_ptr& error)
        : main_adj_(dnaAdj), error_(error) {}
    void operator()();

  private:
    dna_adjust* main_adj_;
    std::exception_ptr& error_;

    // Prevent assignment operator
    adjust_reverse_thread& operator=(const adjust_reverse_thread& rhs);
};

class adjust_combine_thread {
  public:
    adjust_combine_thread(dna_adjust* dnaAdj, std::exception_ptr& error)
        : main_adj_(dnaAdj), error_(error) {}
    void operator()();

  private:
    dna_adjust* main_adj_;
    std::exception_ptr& error_;

    // Prevent assignment operator
    adjust_combine_thread& operator=(const adjust_combine_thread& rhs);
};

class adjust_process_combine_thread {
  public:
    adjust_process_combine_thread(dna_adjust* dnaAdj, const UINT32& id,
                                  std::exception_ptr& error,
                                  std::vector<std::exception_ptr>& cmb_errors)
        : main_adj_(dnaAdj), thread_id_(id), error_(error),
          cmb_errors_(cmb_errors) {}
    void operator()();

  private:
    dna_adjust* main_adj_;
    UINT32 thread_id_;
    std::exception_ptr& error_;
    std::vector<std::exception_ptr>& cmb_errors_;

    // Prevent assignment operator
    adjust_process_combine_thread&
    operator=(const adjust_process_combine_thread& rhs);
};

struct robust_statistics_t {
    double medianSigmaZero;
    struct outlier_type_count {
        char measType;
        UINT32 totalCount;
        UINT32 outlierCount;
    };
    std::vector<outlier_type_count> outlierBreakdown;
    double trimmedChiSquared;
    UINT32 trimmedMeasurementCount;
    int trimmedDOF;
    double trimmedSigmaZero;
    double trimmedChiUpperLimit;
    double trimmedChiLowerLimit;
    UINT32 trimmedPassFail;
};

// This class is exported from the dnaAdjust.dll
#ifdef _MSC_VER
class DNAADJUST_API dna_adjust {
#else
class dna_adjust {
#endif
    friend class DynAdjustPrinter;
    
  public:
    dna_adjust();
    virtual ~dna_adjust();

  private:
    // From Williams, A. 2018, C++ Concurrency in Action, Manning Publications.
    //
    // The standard idiom for preventing copies of a class used to be to declare
    // the copy constructor and copy assignment operator private and then not
    // provide an implementation. This would cause a compile error if any code
    // outside the class in question tried to copy an instance and a link-time
    // error (due to lack of an implementation) if any of the class' member
    // functions or friends tried to copy an instance. With C++11, the committee
    // realized that this was a common idiom but also realized that it's a bit
    // of a hack. The committee therefore provided a more general mechanism that
    // can be applied in other cases too: you can declare a function as deleted.
    //
    // Disallow use of compiler generated functions. For backward compatibility
    // with older C++ compilers, keep the "idiom" of private copy constructor
    // and assignment operator, rather than use delete.
    dna_adjust(const dna_adjust&);
    dna_adjust& operator=(const dna_adjust&);

  public:
    void ShrinkForwardMatrices(const UINT32 currentBlock);
    void
    CarryForwardJunctions(const UINT32 currentBlock, const UINT32 nextBlock);
    bool
    CarryReverseJunctions(const UINT32 currentBlock, const UINT32 nextBlock,
                          bool MT_ReverseOrCombine);
    void PrepareAdjustmentMultiThread();
    void PrepareAdjustmentBlock(const UINT32 block, const UINT32 thread_id = 0);
    bool PrepareAdjustmentReverse(const UINT32 block, bool MT_ReverseOrCombine);
    bool PrepareAdjustmentCombine(const UINT32 block, UINT32& pseudomsrJSLCount,
                                  bool MT_ReverseOrCombine);
    void BackupNormals(const UINT32 block, bool MT_ReverseOrCombine);
    void UpdateEstimatesForward(const UINT32 block);
    void UpdateEstimatesReverse(const UINT32 block, bool MT_ReverseOrCombine);
    void UpdateEstimatesCombine(const UINT32 block, UINT32 prevJSLCount);
    void UpdateEstimatesFinal(const UINT32 currentBlock);
    void UpdateEstimatesFinalNoCombine();

    void GenerateStatistics();
    void PrepareAdjustment(const project_settings& adjustmentSettings);

    inline void CancelAdjustment() { isCancelled_.store(true); }
    inline bool IsCancelled() const { return isCancelled_.load(); };

    inline _ADJUST_STATUS_ GetStatus() const { return adjustStatus_; }
    inline bool IsPreparing() { return isPreparing_; }
    inline bool IsAdjusting() { return isAdjusting_; }
    inline bool ExceptionRaised() {
        return adjustStatus_ == ADJUST_EXCEPTION_RAISED;
    }
    inline void SetExceptionRaised() {
        adjustStatus_ = ADJUST_EXCEPTION_RAISED;
        isPreparing_ = false;
        isCombining_ = false;
        isAdjusting_ = false;
    }
    inline DynAdjustPrinter* GetPrinter() { return printer_.get(); }

    // Printing methods moved to DynAdjustPrinter - use GetPrinter()->methodName() instead
    void PrintAdjustedNetworkStations();
    // void PrintPositionalUncertainty();
    // void PrintNetworkStationCorrections();
    // void PrintAdjustedNetworkMeasurements();
    // void PrintMeasurementsToStation();

    // Printing methods moved to DynAdjustPrinter - use GetPrinter()->methodName() instead
    // bool PrintEstimatedStationCoordinatestoSNX(std::string& sinex_file);
    // void PrintEstimatedStationCoordinatestoDNAXML(const std::string& stnFile,
    //                                               INPUT_FILE_TYPE t,
    //                                               bool flagUnused = false);
    // void PrintEstimatedStationCoordinatestoDNAXML_Y(const std::string& msrFile,
    //                                                 INPUT_FILE_TYPE t);

    void CloseOutputFiles();
    void UpdateBinaryFiles();

    UINT32 CurrentIteration() const;
    UINT32& incrementIteration();
    void initialiseIteration(const UINT32& iteration = 0);

    inline UINT32 CurrentBlock() const { return currentBlock_; };

    inline void SetcurrentBlock(const UINT32 b) { currentBlock_ = b; };

    inline void SetmaxCorr(const double c) {

        std::lock_guard<std::mutex> lock(maxCorrMutex);
        maxCorr_ = c;
    };

    inline UINT32 blockCount() const { return blockCount_; }
    inline std::chrono::milliseconds adjustTime() const {
        return total_time_;
    }
    inline bool processingForward() { return forward_; }
    inline bool processingCombine() { return isCombining_; }
    inline int GetDegreesOfFreedom() const { return degreesofFreedom_; }
    inline UINT32 GetMeasurementCount() const { return measurementParams_; }
    inline UINT32 GetUnknownsCount() const { return unknownParams_; }
    inline double GetChiSquared() const { return chiSquared_; }
    inline double GetSigmaZero() const { return sigmaZero_; }
    inline UINT32 GetPotentialOutlierCount() const {
        return potentialOutlierCount_;
    }
    inline double GetChiSquaredUpperLimit() const {
        return chiSquaredUpperLimit_;
    }
    inline double GetChiSquaredLowerLimit() const {
        return chiSquaredLowerLimit_;
    }
    inline double GetGlobalPelzerRel() const {
        return globalPelzerReliability_;
    }
    inline UINT32 GetTestResult() const { return passFail_; }
    inline bool GetAllFixed() const { return allStationsFixed_; }
    inline const robust_statistics_t& GetRobustStats() const { return robustStats_; }


    void LoadSegmentationFileParameters(const std::string& seg_filename);

    void SerialiseAdjustedVarianceMatrices();
    void DeSerialiseAdjustedVarianceMatrices();

    /////////////////////////////////
    // Messages
    inline bool NewMessagesAvailable() { return iterationQueue_.not_empty(); }
    inline bool GetMessageIteration(UINT32& iteration) {
        return iterationQueue_.front_and_pop(iteration);
    }

    // iteration is indexed from 1.  That is, not zero-indexed.
    inline std::string GetMaxCorrection(const UINT32& iteration) const {
        if (iteration == 0)
            return iterationCorrections_.get_message(iteration); // safe guard
        return iterationCorrections_.get_message(iteration - 1);
    };
    /////////////////////////////

    concurrent_ofstream<std::string> concurrent_adj_ofstream;

    inline void ThreadSafeWritetoAdjFile(const std::string& s) {
        concurrent_adj_ofstream.wrtie(adj_file, s);
    }
    inline void WriteStatstoAdjFile(const std::string& s) {
        if (projectSettings_.o._adj_stat_iteration)
            ThreadSafeWritetoAdjFile(s);
    }
    inline void WriteStntoAdjFile(const std::string& s) {
        if (projectSettings_.o._adj_stn_iteration)
            ThreadSafeWritetoAdjFile(s);
    }
    inline void WriteMsrtoAdjFile(const std::string& s) {
        if (projectSettings_.o._adj_msr_iteration)
            ThreadSafeWritetoAdjFile(s);
    }

    inline void ThreadSafeWritetoDbgFile(const std::string& s) {
        concurrent_adj_ofstream.wrtie(debug_file, s);
    }

    _ADJUST_STATUS_ AdjustNetwork();

    bool isPreparing_;
    bool isAdjusting_;
    bool isCombining_;
    bool forward_;
    bool isFirstTimeAdjustment_;
    bool isIterationComplete_;
    bool isAdjustmentQuestionable_;

    UINT32 blockCount_;
    UINT32 currentBlock_;

    std::chrono::milliseconds total_time_;
    _ADJUST_STATUS_ adjustStatus_;
    vstring statusMessages_;
    UINT32 currentIteration_;
    std::vector<blockMeta_t> v_blockMeta_;

    inline bool FirstBlock(const UINT32& block) const {
        return v_blockMeta_.at(block)._blockFirst;
    }

    inline bool LastBlock(const UINT32& block) const {
        return v_blockMeta_.at(block)._blockLast;
    }

    void
    SignalExceptionAdjustment(const std::string& msg, const UINT32 block_no);

    void CreateMeasurementTally(const UINT32& block = 0);

    ///////////////////////////////////////////////////////////////////////
    // Used for reverse and combine adjustments in multi thread environment
    //
    // Wrappers including try/catch statements
    void SolveTry(bool COMPUTE_INVERSE, const UINT32& block = 0);
    void SolveMTTry(bool COMPUTE_INVERSE, const UINT32& block = 0);

    void Solve(bool COMPUTE_INVERSE, const UINT32& block = 0);
    void SolveMT(bool COMPUTE_INVERSE, const UINT32& block);

    inline bool CombineRequired(const UINT32& block) const {
        if (v_blockMeta_.at(block)._blockLast)
            return false;
        if (v_blockMeta_.at(block)._blockIsolated)
            return false;
        if (v_blockMeta_.at(block)._blockFirst)
            return false;
        return true;
    }

  private:
    bool InitialiseandValidateMsrPointer(const it_vUINT32& _it_block_msr,
                                         it_vmsr_t& _it_msr);
    bool InitialiseMeasurement(pit_vmsr_t _it_msr, bool buildnewMatrices);

    void
    RebuildNormals(const UINT32 block, adjustOperation direction,
                   bool AddConstraintStationstoNormals, bool BackupNormals);
    void UpdateAdjustment(bool iterate);
    void ValidateandFinaliseAdjustment(cpu_timer& tot_time);
    void PrintAdjustmentTime(cpu_timer& time, _TIMER_TYPE_);
    void InitialiseAdjustment();
    void SetDefaultReferenceFrame();
    void LoadNetworkFiles();
    void CreateMsrToStnTally();
    // void CreateStationAppearanceList(const vUINT32& parameterStations);

    // Simultaneous adjustment
    void AdjustSimultaneous();

    // Phased adjustment
    void AdjustPhased();
    void AdjustPhasedForward();
    void AdjustPhasedReverseCombine();
    bool CombinationThreadRequired();

    // Phased adjustment using multiple cores
    void AdjustPhasedMultiThread();

    // Phased adjustment producing rigorous
    // coordinates for block 1 only
    void AdjustPhasedBlock1();
    void AdjustPhasedReverse();

    // Adjustment helps
    void AddDiscontinuitySites(vstring& constraintStns);
    void LoadStationMap(pv_string_uint32_pair stnsMap,
                        const std::string& stnmap_file);
    void ResizeMatrixVectors();
    void LoadPhasedBlocks();
    void LoadSegmentationFile();
    void LoadSegmentationMetrics();
    void RemoveDuplicateStations(vUINT32& vStns);
    void InitialiseTypeBUncertainties();

    // Adjusted measurement sorting
    void SortMeasurementsbyType(v_uint32_u32u32_pair& msr_block);
    void SortMeasurementsbyFromStn(v_uint32_u32u32_pair& msr_block);
    void SortMeasurementsbyToStn(v_uint32_u32u32_pair& msr_block);
    void SortMeasurementsbyValue(v_uint32_u32u32_pair& msr_block);
    void SortMeasurementsbyResidual(v_uint32_u32u32_pair& msr_block);
    void SortMeasurementsbyAdjSD(v_uint32_u32u32_pair& msr_block);
    void SortMeasurementsbyNstat(v_uint32_u32u32_pair& msr_block);

    // Adjustment station sorting
    void SortStationsbyFileOrder(vUINT32& v_blockStations);
    void SortStationsbyID(vUINT32& v_blockStations);

    // Stage helps
    void
    SetRegionOffsets(const UINT32& block, const int file_count = 0, ...);
    void SetRegionOffset(vmat_file_map& file_map);
    void SetMapRegions(const int file_count = 0, ...);
    void PrepareMemoryMapRegions(const UINT32& block,
                                 const int file_count = 0, ...);
    void OffloadBlockToDisk(const UINT32& block);
    void OffloadBlockToMappedFile(const UINT32& block);
    void SerialiseBlockToDisk(const UINT32& block);
    void SerialiseBlockToMappedFile(const UINT32& block,
                                    const int file_count = 0, ...);
    void DeserialiseBlockFromMappedFile(const UINT32& block,
                                        const int count = 0, ...);
    void UnloadBlock(const UINT32& block, const int file_count = 0, ...);
    void PurgeMatricesFromDisk();

    // Helpers
    void AddMsrtoMeasMinusComp(pit_vmsr_t _it_msr, const UINT32& design_row,
                               const double comp_msr, matrix_2d* measMinusComp,
                               bool printBlock = true);
    void AddMsrtoNormalsVar(const UINT32& design_row, const UINT32& stn,
                            matrix_2d* normals, matrix_2d* design,
                            matrix_2d* AtVinv);
    void AddMsrtoNormalsCoVar2(const UINT32& design_row, const UINT32& stn1,
                               const UINT32& stn2, matrix_2d* normals,
                               matrix_2d* design, matrix_2d* AtVinv);
    void AddMsrtoNormalsCoVar3(const UINT32& design_row, const UINT32& stn1,
                               const UINT32& stn2, const UINT32& stn3,
                               matrix_2d* normals, matrix_2d* design,
                               matrix_2d* AtVinv);

    inline void AddMsrtoDesign(const UINT32& design_row, const UINT32& stn,
                               const double& dmdx, const double& dmdy,
                               const double& dmdz, matrix_2d* design);
    inline void
    AddMsrtoDesign_L(const UINT32& design_row, const UINT32& stn1,
                     const UINT32& stn2, const double dmdx1, const double dmdy1,
                     const double dmdz1, const double dmdx2, const double dmdy2,
                     const double dmdz2, matrix_2d* design);
    inline void AddMsrtoDesign_BCEKMSVZ(const UINT32& design_row,
                                        const UINT32& stn1, const UINT32& stn2,
                                        const double dmdx, const double dmdy,
                                        const double dmdz, matrix_2d* design);
    inline void AddElementtoDesign(const UINT32& row, const UINT32& col,
                                   const double value, matrix_2d* design);

    void
    UpdateDesignNormalMeasMatrices(pit_vmsr_t _it_msr, UINT32& design_row,
                                   bool buildnewMatrices, const UINT32& block,
                                   bool MT_ReverseOrCombine);

    void
    UpdateDesignNormalMeasMatrices_A(pit_vmsr_t _it_msr, UINT32& design_row,
                                     const UINT32& block,
                                     matrix_2d* measMinusComp,
                                     matrix_2d* estimatedStations,
                                     matrix_2d* normals, matrix_2d* design,
                                     matrix_2d* AtVinv, bool buildnewMatrices);
    void
    UpdateDesignNormalMeasMatrices_BK(pit_vmsr_t _it_msr, UINT32& design_row,
                                      const UINT32& block,
                                      matrix_2d* measMinusComp,
                                      matrix_2d* estimatedStations,
                                      matrix_2d* normals, matrix_2d* design,
                                      matrix_2d* AtVinv, bool buildnewMatrices);
    void
    UpdateDesignNormalMeasMatrices_C(pit_vmsr_t _it_msr, UINT32& design_row,
                                     const UINT32& block,
                                     matrix_2d* measMinusComp,
                                     matrix_2d* estimatedStations,
                                     matrix_2d* normals, matrix_2d* design,
                                     matrix_2d* AtVinv, bool buildnewMatrices);
    void UpdateDesignNormalMeasMatrices_CEM(
        pit_vmsr_t _it_msr, UINT32& design_row, const UINT32& block,
        matrix_2d* measMinusComp, matrix_2d* estimatedStations,
        matrix_2d* normals, matrix_2d* design, matrix_2d* AtVinv,
        bool buildnewMatrices);
    void
    UpdateDesignNormalMeasMatrices_D(pit_vmsr_t _it_msr, UINT32& design_row,
                                     const UINT32& block,
                                     matrix_2d* measMinusComp,
                                     matrix_2d* estimatedStations,
                                     matrix_2d* normals, matrix_2d* design,
                                     matrix_2d* AtVinv, bool buildnewMatrices);
    void
    UpdateDesignNormalMeasMatrices_E(pit_vmsr_t _it_msr, UINT32& design_row,
                                     const UINT32& block,
                                     matrix_2d* measMinusComp,
                                     matrix_2d* estimatedStations,
                                     matrix_2d* normals, matrix_2d* design,
                                     matrix_2d* AtVinv, bool buildnewMatrices);
    void UpdateDesignMeasMatrices_GX(pit_vmsr_t _it_msr, UINT32& design_row,
                                     matrix_2d* measMinusComp,
                                     matrix_2d* estimatedStations,
                                     matrix_2d* design, const UINT32& stn1,
                                     const UINT32& stn2, bool buildnewMatrices);
    void
    UpdateDesignNormalMeasMatrices_G(pit_vmsr_t _it_msr, UINT32& design_row,
                                     const UINT32& block,
                                     matrix_2d* measMinusComp,
                                     matrix_2d* estimatedStations,
                                     matrix_2d* normals, matrix_2d* design,
                                     matrix_2d* AtVinv, bool buildnewMatrices);
    void
    UpdateDesignNormalMeasMatrices_H(pit_vmsr_t _it_msr, UINT32& design_row,
                                     const UINT32& block,
                                     matrix_2d* measMinusComp,
                                     matrix_2d* estimatedStations,
                                     matrix_2d* normals, matrix_2d* design,
                                     matrix_2d* AtVinv, bool buildnewMatrices);
    void
    UpdateDesignNormalMeasMatrices_HR(pit_vmsr_t _it_msr, UINT32& design_row,
                                      const UINT32& block,
                                      matrix_2d* measMinusComp,
                                      matrix_2d* estimatedStations,
                                      matrix_2d* normals, matrix_2d* design,
                                      matrix_2d* AtVinv, bool buildnewMatrices);
    void
    UpdateDesignNormalMeasMatrices_I(pit_vmsr_t _it_msr, UINT32& design_row,
                                     const UINT32& block,
                                     matrix_2d* measMinusComp,
                                     matrix_2d* estimatedStations,
                                     matrix_2d* normals, matrix_2d* design,
                                     matrix_2d* AtVinv, bool buildnewMatrices);
    void
    UpdateDesignNormalMeasMatrices_IP(pit_vmsr_t _it_msr, UINT32& design_row,
                                      const UINT32& block,
                                      matrix_2d* measMinusComp,
                                      matrix_2d* estimatedStations,
                                      matrix_2d* normals, matrix_2d* design,
                                      matrix_2d* AtVinv, bool buildnewMatrices);
    void
    UpdateDesignNormalMeasMatrices_J(pit_vmsr_t _it_msr, UINT32& design_row,
                                     const UINT32& block,
                                     matrix_2d* measMinusComp,
                                     matrix_2d* estimatedStations,
                                     matrix_2d* normals, matrix_2d* design,
                                     matrix_2d* AtVinv, bool buildnewMatrices);
    void
    UpdateDesignNormalMeasMatrices_JQ(pit_vmsr_t _it_msr, UINT32& design_row,
                                      const UINT32& block,
                                      matrix_2d* measMinusComp,
                                      matrix_2d* estimatedStations,
                                      matrix_2d* normals, matrix_2d* design,
                                      matrix_2d* AtVinv, bool buildnewMatrices);
    void
    UpdateDesignNormalMeasMatrices_L(pit_vmsr_t _it_msr, UINT32& design_row,
                                     const UINT32& block,
                                     matrix_2d* measMinusComp,
                                     matrix_2d* estimatedStations,
                                     matrix_2d* normals, matrix_2d* design,
                                     matrix_2d* AtVinv, bool buildnewMatrices);
    void
    UpdateDesignNormalMeasMatrices_M(pit_vmsr_t _it_msr, UINT32& design_row,
                                     const UINT32& block,
                                     matrix_2d* measMinusComp,
                                     matrix_2d* estimatedStations,
                                     matrix_2d* normals, matrix_2d* design,
                                     matrix_2d* AtVinv, bool buildnewMatrices);
    void
    UpdateDesignNormalMeasMatrices_P(pit_vmsr_t _it_msr, UINT32& design_row,
                                     const UINT32& block,
                                     matrix_2d* measMinusComp,
                                     matrix_2d* estimatedStations,
                                     matrix_2d* normals, matrix_2d* design,
                                     matrix_2d* AtVinv, bool buildnewMatrices);
    void
    UpdateDesignNormalMeasMatrices_Q(pit_vmsr_t _it_msr, UINT32& design_row,
                                     const UINT32& block,
                                     matrix_2d* measMinusComp,
                                     matrix_2d* estimatedStations,
                                     matrix_2d* normals, matrix_2d* design,
                                     matrix_2d* AtVinv, bool buildnewMatrices);
    void
    UpdateDesignNormalMeasMatrices_R(pit_vmsr_t _it_msr, UINT32& design_row,
                                     const UINT32& block,
                                     matrix_2d* measMinusComp,
                                     matrix_2d* estimatedStations,
                                     matrix_2d* normals, matrix_2d* design,
                                     matrix_2d* AtVinv, bool buildnewMatrices);
    void
    UpdateDesignNormalMeasMatrices_S(pit_vmsr_t _it_msr, UINT32& design_row,
                                     const UINT32& block,
                                     matrix_2d* measMinusComp,
                                     matrix_2d* estimatedStations,
                                     matrix_2d* normals, matrix_2d* design,
                                     matrix_2d* AtVinv, bool buildnewMatrices);
    void
    UpdateDesignNormalMeasMatrices_V(pit_vmsr_t _it_msr, UINT32& design_row,
                                     const UINT32& block,
                                     matrix_2d* measMinusComp,
                                     matrix_2d* estimatedStations,
                                     matrix_2d* normals, matrix_2d* design,
                                     matrix_2d* AtVinv, bool buildnewMatrices);
    void UpdateDesignNormalMeasMatrices_X(
        pit_vmsr_t _it_msr, UINT32& design_row, const UINT32& block,
        matrix_2d* measMinusComp, matrix_2d* estimatedStations,
        matrix_2d* design, matrix_2d* AtVinv, bool buildnewMatrices);
    void UpdateDesignNormalMeasMatrices_Y(
        pit_vmsr_t _it_msr, UINT32& design_row, const UINT32& block,
        matrix_2d* measMinusComp, matrix_2d* estimatedStations,
        matrix_2d* design, matrix_2d* AtVinv, bool buildnewMatrices);
    void
    UpdateDesignNormalMeasMatrices_Z(pit_vmsr_t _it_msr, UINT32& design_row,
                                     const UINT32& block,
                                     matrix_2d* measMinusComp,
                                     matrix_2d* estimatedStations,
                                     matrix_2d* normals, matrix_2d* design,
                                     matrix_2d* AtVinv, bool buildnewMatrices);

    void UpdateIgnoredMeasurements(pit_vmsr_t _it_msr,
                                   bool storeOriginalMeasurement);
    void UpdateIgnoredMeasurements_A(pit_vmsr_t _it_msr,
                                     bool storeOriginalMeasurement);
    void UpdateIgnoredMeasurements_B(pit_vmsr_t _it_msr,
                                     bool storeOriginalMeasurement);
    void UpdateIgnoredMeasurements_BK(pit_vmsr_t _it_msr);
    void UpdateIgnoredMeasurements_C(pit_vmsr_t _it_msr,
                                     bool storeOriginalMeasurement);
    void UpdateIgnoredMeasurements_CEM(pit_vmsr_t _it_msr);
    void UpdateIgnoredMeasurements_D(pit_vmsr_t _it_msr,
                                     bool storeOriginalMeasurement);
    void UpdateIgnoredMeasurements_E(pit_vmsr_t _it_msr,
                                     bool storeOriginalMeasurement);
    void UpdateIgnoredMeasurements_G(pit_vmsr_t _it_msr,
                                     bool storeOriginalMeasurement);
    void UpdateIgnoredMeasurements_GX(pit_vmsr_t _it_msr,
                                      bool storeOriginalMeasurement);
    void UpdateIgnoredMeasurements_H(pit_vmsr_t _it_msr,
                                     bool storeOriginalMeasurement);
    void UpdateIgnoredMeasurements_HR(pit_vmsr_t _it_msr);
    void UpdateIgnoredMeasurements_I(pit_vmsr_t _it_msr,
                                     bool storeOriginalMeasurement);
    void UpdateIgnoredMeasurements_IP(pit_vmsr_t _it_msr);
    void UpdateIgnoredMeasurements_J(pit_vmsr_t _it_msr,
                                     bool storeOriginalMeasurement);
    void UpdateIgnoredMeasurements_JQ(pit_vmsr_t _it_msr);
    void UpdateIgnoredMeasurements_K(pit_vmsr_t _it_msr,
                                     bool storeOriginalMeasurement);
    void UpdateIgnoredMeasurements_L(pit_vmsr_t _it_msr,
                                     bool storeOriginalMeasurement);
    void UpdateIgnoredMeasurements_M(pit_vmsr_t _it_msr,
                                     bool storeOriginalMeasurement);
    void UpdateIgnoredMeasurements_P(pit_vmsr_t _it_msr,
                                     bool storeOriginalMeasurement);
    void UpdateIgnoredMeasurements_Q(pit_vmsr_t _it_msr,
                                     bool storeOriginalMeasurement);
    void UpdateIgnoredMeasurements_R(pit_vmsr_t _it_msr,
                                     bool storeOriginalMeasurement);
    void UpdateIgnoredMeasurements_S(pit_vmsr_t _it_msr,
                                     bool storeOriginalMeasurement);
    void UpdateIgnoredMeasurements_V(pit_vmsr_t _it_msr,
                                     bool storeOriginalMeasurement);
    void UpdateIgnoredMeasurements_X(pit_vmsr_t _it_msr,
                                     bool storeOriginalMeasurement);
    void UpdateIgnoredMeasurements_Y(pit_vmsr_t _it_msr,
                                     bool storeOriginalMeasurement);
    void UpdateIgnoredMeasurements_Z(pit_vmsr_t _it_msr,
                                     bool storeOriginalMeasurement);

    void LoadVarianceScaling(it_vmsr_t _it_msr, double& vScale, double& pScale,
                             double& lScale, double& hScale, bool& scaleMatrix,
                             bool& scalePartial);
    void LoadVarianceMatrix_D(it_vmsr_t _it_msr, matrix_2d* var_dirn,
                              bool buildnewMatrices);
    void LoadVarianceMatrix_G(it_vmsr_t _it_msr, matrix_2d* var_cart);
    void LoadVarianceMatrix_X(it_vmsr_t _it_msr, matrix_2d* var_cart);
    void LoadVarianceMatrix_Y(it_vmsr_t _it_msr, matrix_2d* var_cart,
                              const _COORD_TYPE_ coordType);

    // Printing methods moved to DynAdjustPrinter - use printer_->methodName() instead
    // void PrintMsrVarianceMatrixException(const it_vmsr_t& _it_msr,
    //                                      const std::runtime_error& e,
    //                                      std::stringstream& ss,
    //                                      const std::string& calling_function,
    //                                      const UINT32 msr_count = 0);

    // Estimated stations matrix
    void PrepareStationandVarianceMatrices(const UINT32& block);
    void PrepareMappedRegions(const UINT32& block);
    void
    PopulateEstimatedStationMatrix(const UINT32& block, UINT32& unknownParams);
    void CarryStnEstimatesandVariancesForward(const UINT32& thisBlock,
                                              const UINT32& nextBlock);
    void CarryStnEstimatesandVariancesReverse(const UINT32& nextBlock,
                                              const UINT32& thisBlock,
                                              bool MT_ReverseOrCombine);
    void CarryStnEstimatesandVariancesCombine(const UINT32& nextBlock,
                                              const UINT32& thisBlock,
                                              UINT32& pseudomsrJSLCount,
                                              bool MT_ReverseOrCombine);

    // Update AtVinv based on new design matrix elements
    void
    UpdateAtVinv(pit_vmsr_t _it_msr, const UINT32& stn1, const UINT32& stn2,
                 const UINT32& stn3, UINT32& design_row, matrix_2d* design,
                 matrix_2d* AtVinv, bool buildnewMatrices = true);
    // Direction (cluster measurement)
    void
    UpdateAtVinv_D(const UINT32& stn1, const UINT32& stn2, const UINT32& stn3,
                   const UINT32& angle, const UINT32& angle_count,
                   UINT32& design_row, UINT32& design_row_begin,
                   matrix_2d* Vinv, matrix_2d* design, matrix_2d* AtVinv);

    // Update Normals based on new design matrix elements (i.e. when non-GPS
    // msrs are involved)
    void UpdateNormals(const UINT32& block, bool MT_ReverseOrCombine);
    // three station measurements
    void
    UpdateNormals_A(const UINT32& stn1, const UINT32& stn2, const UINT32& stn3,
                    UINT32& design_row, matrix_2d* normals, matrix_2d* design,
                    matrix_2d* AtVinv);
    // two station measurements
    void UpdateNormals_BCEKLMSVZ(const UINT32& stn1, const UINT32& stn2,
                                 UINT32& design_row, matrix_2d* normals,
                                 matrix_2d* design, matrix_2d* AtVinv);
    // single station measurements
    void UpdateNormals_HIJPQR(const UINT32& stn1, UINT32& design_row,
                              matrix_2d* normals, matrix_2d* design,
                              matrix_2d* AtVinv);
    // Direction (cluster measurement)
    void
    UpdateNormals_D(const UINT32& block, it_vmsr_t& _it_msr, UINT32& design_row,
                    matrix_2d* normals, matrix_2d* design, matrix_2d* AtVinv);
    // GPS specific
    void
    UpdateNormals_G(const UINT32& stn1, const UINT32& stn2, UINT32& design_row,
                    matrix_2d* normals, matrix_2d* AtVinv);
    void
    UpdateNormals_X(const UINT32& block, it_vmsr_t& _it_msr, UINT32& design_row,
                    matrix_2d* normals, matrix_2d* design, matrix_2d* AtVinv);
    void
    UpdateNormals_Y(const UINT32& block, it_vmsr_t& _it_msr, UINT32& design_row,
                    matrix_2d* normals, matrix_2d* AtVinv);

    void OutputLargestCorrection(std::string& formatted_msg);

    // Design/Normals matrix
    void PrepareDesignAndMsrMnsCmpMatrices(const UINT32& block);
    void PrepareDesignAndMsrMnsCmpMatricesStage(const UINT32& block);
    void FillDesignNormalMeasurementsMatrices(bool buildnewMatrices,
                                              const UINT32& block,
                                              bool MT_ReverseOrCombine);

    // void RecomputeMeasurementsCommonJunctions(const UINT32& nextBlock, const
    // UINT32& thisBlock, const UINT32& prevBlock);
    void AddConstraintStationstoNormalsForward(const UINT32& block);
    void AddConstraintStationstoNormalsReverse(const UINT32& block,
                                               bool MT_ReverseOrCombine);
    void AddConstraintStationstoNormalsCombine(const UINT32& block,
                                               bool MT_ReverseOrCombine);
    void AddConstraintStationstoNormalsSimultaneous(const UINT32& block);
    void FormConstraintStationVarianceMatrix(const it_vUINT32& _it_param_stn,
                                             matrix_2d& var_cart);

    void UpdateNormalsR(const UINT32& block);
    void UpdateNormalsC(const UINT32& block);
    void CarryStnEstimatesandVariancesReverseR(const UINT32& nextBlock,
                                               const UINT32& thisBlock);

    static void fillSinexExample();

    void ComputeStatistics();
    void ComputeStatisticsOnIteration();

    void ComputeTstatistics();
    void UpdateMsrTstatistic(const UINT32& block = 0);
    void UpdateMsrTstatistic_D(it_vmsr_t& _it_msr);
    void UpdateMsrTstatistic_GXY(it_vmsr_t& _it_msr);

    void ComputeandPrintAdjMsrOnIteration();
    void ComputeandPrintAdjMsrBlockOnIteration(const UINT32& block,
                                               v_uint32_u32u32_pair msr_block,
                                               bool printHeader);
    void ComputeAdjMsrBlockOnIteration(const UINT32& block);

    void ComputeAdjustedMsrPrecisions();

    void ComputeChiSquareNetwork();
    void ComputeChiSquare(const UINT32& block);
    void ComputeChiSquareSimultaneous();
    void ComputeChiSquarePhased(const UINT32& block);

    void ComputeTestStat(const double& dof, double& chiUpper, double& chiLower,
                         double& sigmaZero, UINT32& passFail);
    void ComputeBlockTestStat(const UINT32& block);
    void ComputeGlobalTestStat();
    void ComputeGlobalNetStat();

    void ComputeRobustStatistics(robust_statistics_t& stats);
    void ComputeMedianSigmaZero(robust_statistics_t& stats);
    void ComputeOutlierBreakdown(robust_statistics_t& stats);
    void ComputeTrimmedChiSquare(robust_statistics_t& stats);

    void ComputePrecisionAdjMsrs(const UINT32& block = 0);
    void ComputePrecisionAdjMsrs_A(const UINT32& block, const UINT32& stn1,
                                   const UINT32& stn2, const UINT32& stn3,
                                   matrix_2d* design,
                                   matrix_2d* aposterioriVariances,
                                   UINT32& design_row, UINT32& precadjmsr_row);
    void ComputePrecisionAdjMsrs_D(const UINT32& block, it_vmsr_t& _it_msr,
                                   matrix_2d* design,
                                   matrix_2d* aposterioriVariances,
                                   UINT32& design_row, UINT32& precadjmsr_row);
    void
    ComputePrecisionAdjMsrs_BCEKLMSVZ(const UINT32& block, const UINT32& stn1,
                                      const UINT32& stn2, matrix_2d* design,
                                      matrix_2d* aposterioriVariances,
                                      UINT32& design_row,
                                      UINT32& precadjmsr_row);
    void
    ComputePrecisionAdjMsrs_HIJPQR(const UINT32& block, const UINT32& stn1,
                                   matrix_2d* design,
                                   matrix_2d* aposterioriVariances,
                                   UINT32& design_row, UINT32& precadjmsr_row);
    void ComputePrecisionAdjMsrs_GX(const UINT32& block, it_vmsr_t& _it_msr,
                                    matrix_2d* aposterioriVariances,
                                    UINT32& design_row, UINT32& precadjmsr_row);
    void ComputePrecisionAdjMsrs_Y(const UINT32& block, it_vmsr_t& _it_msr,
                                   matrix_2d* aposterioriVariances,
                                   UINT32& design_row, UINT32& precadjmsr_row);

    void UpdateMsrRecords(const UINT32& block = 0);
    void UpdateMsrRecord(const UINT32& block, it_vmsr_t& _it_msr,
                         const UINT32& msr_row, const UINT32& precadjmsr_row,
                         const double& measPrec);
    void UpdateMsrRecords_D(const UINT32& block, it_vmsr_t& _it_msr,
                            UINT32& msr_row, UINT32& precadjmsr_row);
    void UpdateMsrRecords_GXY(const UINT32& block, it_vmsr_t& _it_msr,
                              UINT32& msr_row, UINT32& precadjmsr_row);
    void UpdateMsrRecordStats(it_vmsr_t& _it_msr, const double& measPrec);

    void ComputeGlobalPelzer();
    void ComputeGlobalPelzer_D(it_vmsr_t& _it_msr, UINT32& numMsr, double& sum);
    void
    ComputeGlobalPelzer_GXY(it_vmsr_t& _it_msr, UINT32& numMsr, double& sum);

    void ComputeChiSquare_ABCEHIJKLMPQRSVZ(const it_vmsr_t& _it_msr,
                                           UINT32& measurement_index,
                                           matrix_2d* measMinusComp);
    void ComputeChiSquare_D(it_vmsr_t& _it_msr, UINT32& measurement_index,
                            matrix_2d* measMinusComp);
    void ComputeChiSquare_G(const it_vmsr_t& _it_msr, UINT32& measurement_index,
                            matrix_2d* measMinusComp);
    void
    ComputeChiSquare_XY(const it_vmsr_t& _it_msr, UINT32& measurement_index,
                        matrix_2d* measMinusComp);

    void
    FormInverseVarianceMatrix(matrix_2d* vmat, bool LOWER_IS_CLEARED = false);
    void
    FormInverseGPSVarianceMatrix(const it_vmsr_t& _it_msr, matrix_2d* vmat);
    bool
    FormInverseVarianceMatrixReduced(it_vmsr_t _it_msr, matrix_2d* var_cart,
                                     const std::string& method_name);

    // Printing methods moved to DynAdjustPrinter - use printer_->methodName() instead
    // void PrintStatistics(bool printPelzer = true);

    // Output files
    void OpenOutputFileStreams();
    // Printing methods moved to DynAdjustPrinter - use printer_->methodName() instead
    // void PrintOutputFileHeaderInfo();
    // Printing methods moved to DynAdjustPrinter - use printer_->methodName() instead
    // void
    // PrintCompMeasurements(const UINT32& block, const std::string& msg = "");
    // void
    // PrintCompMeasurementsAngular(const char cardinal, const double& computed,
    //                              const double& correction,
    //                              const it_vmsr_t& _it_msr);
    // void
    // PrintCompMeasurementsLinear(const char cardinal, const double& computed,
    //                             const double& correction,
    //                             const it_vmsr_t& _it_msr);
    // void PrintCompMeasurements_A(const UINT32& block, it_vmsr_t& _it_msr,
    //                              UINT32& design_row,
    //                              printMeasurementsMode printMode);
    // void PrintCompMeasurements_CELMS(it_vmsr_t& _it_msr, UINT32& design_row,
    //                                  printMeasurementsMode printMode);
    // void PrintCompMeasurements_D(it_vmsr_t& _it_msr, UINT32& design_row,
    //                              bool printIgnored = false);
    // void PrintCompMeasurements_HR(const UINT32& block, it_vmsr_t& _it_msr,
    //                               UINT32& design_row,
    //                               printMeasurementsMode printMode);
    // void PrintCompMeasurements_IJPQ(const UINT32& block, it_vmsr_t& _it_msr,
    //                                 UINT32& design_row,
    //                                 printMeasurementsMode printMode);
    // void PrintCompMeasurements_BKVZ(const UINT32& block, it_vmsr_t& _it_msr,
    //                                 UINT32& design_row,
    //                                 printMeasurementsMode printMode);
    // void PrintCompMeasurements_GXY(const UINT32& block, it_vmsr_t& _it_msr,
    //                                UINT32& design_row,
    //                                printMeasurementsMode printMode);
    // void PrintCompMeasurements_YLLH(it_vmsr_t& _it_msr, UINT32& design_row);

    // Printing methods moved to DynAdjustPrinter - use printer_->methodName() instead
    // void
    // PrintMeasurementsAngular(const char cardinal, const double& measurement,
    //                          const double& correction, const it_vmsr_t& _it_msr,
    //                          bool printAdjMsr = true);
    // void
    // PrintMeasurementsLinear(const char cardinal, const double& measurement,
    //                         const double& correction, const it_vmsr_t& _it_msr,
    //                         bool printAdjMsr = true);
    // void
    // PrintMeasurementCorrection(const char cardinal, const it_vmsr_t& _it_msr);
    // void PrintMeasurementDatabaseID(const it_vmsr_t& _it_msr,
    //                                 bool initialise_dbindex = true);

    void FormUniqueMsrList();
    // Printing methods moved to DynAdjustPrinter - use printer_->methodName() instead
    // void PrintAdjMeasurementsHeader(bool printHeader,
    //                                 const std::string& table_heading,
    //                                 printMeasurementsMode printMode,
    //                                 UINT32 block, bool printBlocks = false);
    // void PrintAdjMeasurements(v_uint32_u32u32_pair msr_block, bool printHeader);
    // void
    // PrintAdjMeasurementsAngular(const char cardinal, const it_vmsr_t& _it_msr,
    //                             bool initialise_dbindex = true);
    // void
    // PrintAdjMeasurementsLinear(const char cardinal, const it_vmsr_t& _it_msr,
    //                            bool initialise_dbindex = true);
    // void PrintAdjGNSSAlternateUnits(it_vmsr_t& _it_msr,
    //                                 const uint32_uint32_pair& b_pam);
    // void
    // PrintAdjMeasurementStatistics(const char cardinal, const it_vmsr_t& _it_msr,
    //                               bool initialise_dbindex = true);

    // void PrintIgnoredAdjMeasurements(bool printHeader);

    // Printing methods moved to DynAdjustPrinter - use printer_->methodName() instead
    // void PrintAdjMeasurements_A(it_vmsr_t& _it_msr);
    // void PrintAdjMeasurements_CELMS(it_vmsr_t& _it_msr);
    // void PrintAdjMeasurements_D(it_vmsr_t& _it_msr);
    // void PrintAdjMeasurements_HR(it_vmsr_t& _it_msr);
    // void PrintAdjMeasurements_IJPQ(it_vmsr_t& _it_msr);
    // void PrintAdjMeasurements_BKVZ(it_vmsr_t& _it_msr);
    // void PrintAdjMeasurements_GXY(it_vmsr_t& _it_msr,
    //                               const uint32_uint32_pair& b_pam);
    // void PrintAdjMeasurements_YLLH(it_vmsr_t& _it_msr);
    void ReduceYLLHMeasurementsforPrinting(vmsr_t& y_msr, matrix_2d& mpositions,
                                           printMeasurementsMode print_mode);


    void PrintCorStations(std::ostream& cor_file, const UINT32& block);
    // void
    // PrintCorStation(std::ostream& os, const UINT32& block, const UINT32& stn,
    //                 const UINT32& mat_idx, const matrix_2d* stationEstimates);

    // Printing methods moved to DynAdjustPrinter - use printer_->methodName() instead
    // void
    // PrintPosUncertainty(std::ostream& os, /*ostream* csv,*/ const UINT32& block,
    //                     const UINT32& stn, const UINT32& mat_idx,
    //                     const matrix_2d* stationVariances,
    //                     const UINT32& map_idx = 0,
    //                     const vUINT32* blockStations = NULL);
    // void PrintPosUncertainties(std::ostream& os, const UINT32& block,
    //                            const matrix_2d* stationVariances);
    // void PrintPosUncertaintiesUniqueList(std::ostream& os,
    //                                      const v_mat_2d* stationVariances);
    // void PrintPosUncertaintiesHeader(std::ostream& os);

    void UpdateGeographicCoordsPhased(const UINT32& block,
                                      matrix_2d* estimatedStations);
    void UpdateGeographicCoords();

    inline UINT32
    GetBlkMatrixElemStn1(const UINT32& block, const pit_vmsr_t _it_msr) {
        return v_blockStationsMap_.at(block)[(*_it_msr)->station1] * 3;
    }
    inline UINT32
    GetBlkMatrixElemStn2(const UINT32& block, const pit_vmsr_t _it_msr) {
        return v_blockStationsMap_.at(block)[(*_it_msr)->station2] * 3;
    }
    inline UINT32
    GetBlkMatrixElemStn3(const UINT32& block, const pit_vmsr_t _it_msr) {
        return v_blockStationsMap_.at(block)[(*_it_msr)->station3] * 3;
    }

    void debug_BlockInformation(const UINT32& currentBlock,
                                const std::string& adjustment_method);
    void debug_SolutionInformation(const UINT32& currentBlock);

    CDnaDatum datum_;
    CDnaProjection projection_;
    double _var_C, _var_F /*, _var_Cg, _var_Fg*/;
    matrix_2d _inv_var_cart_c, _inv_var_cart_f;

    project_settings projectSettings_;

    binary_file_meta_t bst_meta_;
    binary_file_meta_t bms_meta_;

    vstn_t bstBinaryRecords_;
    vmsr_t bmsBinaryRecords_;
    vASL vAssocStnList_;
    v_aml_pair vAssocMsrList_;

    bool blowupDetected_;
    double currentRelaxation_;
    double prevMaxCorr_;
    UINT32 oscillationCount_;
    UINT32 monotonicDecreaseCount_;
    void UpdateAdaptiveRelaxation();

    UINT32 bmsr_count_;
    UINT32 bstn_count_;
    UINT32 asl_count_;

    std::ofstream debug_file;
    std::ofstream adj_file;
    std::ofstream xyz_file;

    std::unique_ptr<DynAdjustPrinter> printer_;

    UINT32 blockLargeCorr_;
    double largestCorr_;

    UINT16 PRECISION_SEC_MSR, PRECISION_SEC_STN;
    UINT16 PRECISION_MTR_MSR, PRECISION_MTR_STN;

    // ----------------------------------------------
    // Adjustment parameters
    //
    // Entire network
    UINT32 measurementParams_; // number of raw measurements (less ignored
                               // measurements)
    UINT32 measurementCount_;  // number of raw measurements and constrained
                               // stations (used to resize matrices)
    UINT32 unknownParams_;     // number of free stations
    UINT32 unknownsCount_;     // number of all stations (constrained and free)
    double chiSquared_, chiSquaredStage_;
    double sigmaZero_;
    double sigmaZeroSqRt_;
    double chiSquaredUpperLimit_;
    double chiSquaredLowerLimit_;
    double globalPelzerReliability_;
    int degreesofFreedom_;
    UINT32 passFail_;
    double maxCorr_;
    double criticalValue_;
    UINT32 potentialOutlierCount_;
    bool allStationsFixed_;
    robust_statistics_t robustStats_;

    message_bank<std::string> iterationCorrections_;

    // For each block
    vUINT32 v_ContiguousNetList_; // vector of contiguous network IDs
                                  // (corresponding to each block)
    vsummary_t v_statSummary_;
    vUINT32
        v_pseudoMeasCountFwd_; // number of pseudo measurements in forward pass
    vUINT32 v_measurementParams_; // number of raw measurements (less ignored
                                  // measurements)
    vUINT32 v_measurementCount_;  // number of raw measurements and constrained
                                  // stations (used to resize matrices)
    vUINT32 v_measurementVarianceCount_; // number of raw measurements and
                                         // constrained stations (used to resize
                                         // matrices)
    vUINT32 v_unknownParams_;            // number of free stations
    vUINT32 v_unknownsCount_; // number of all station elements (X, Y and Z -
                              // constrained and free)
    vdouble v_sigmaZero_;
    vdouble v_chiSquaredUpperLimit_;
    vdouble v_chiSquaredLowerLimit_;
    vUINT32 v_passFail_;
    // ----------------------------------------------

    v_uint32_uint32_map v_blockStationsMap_;
    v_u32u32_uint32_pair
        v_blockStationsMapUnique_; // [ [station, block index] , [block] ]
    void BuildUniqueBlockStationMap();
    void BuildSimultaneousStnAppearance();

    vvUINT32 v_ISL_; // Inner stations
    vvUINT32 v_JSL_; // Junction stations
    vvUINT32 v_CML_; // Measurements.  Each index refers to:
                     //  - Non-ignored measurements
                     //  - The fist measurement in a cluster

    v_uint32_u32u32_pair v_msr_block_; // map of measurements and block number,
                                       // prec adj msr matrix index

    vmsrtally v_msrTally_;    // total measurements tally
    vmsrtally v_stnmsrTally_; // Measurements to station tally

    vv_stn_appear v_paramStnAppearance_; // The appearance of stations in blocks
    vUINT32 v_parameterStationCount_;
    vvUINT32 v_parameterStationList_; // Inner and Junction stations, sorted.
                                      // See LoadSegmentationFile()

    type_b_uncertainty typeBUncertaintyGlobal_; // Type B uncertainties to be
                                                // applied to all stations
    v_type_b_uncertainty
        v_typeBUncertaintiesLocal_; // Type B uncertainties to be applied to
                                    // selected stations
    v_type_b_method v_typeBUncertaintyMethod_; // How Type B uncertainties are
                                               // to be handled for each station
    v_uint32_uint32_pair
        v_stationTypeBMap_; // Maps station ID to the index of the local type B
                            // uncertainties vector

    // ----------------------------------------------
    // Adjustment matrices for phased adjustment
    // In the case where MULTI_THREAD is defined, these
    // matrices are used for the forward thread
    v_mat_2d v_normals_;           // vector of ((At * V-1) * A) matrices
    v_mat_2d v_normalsR_;          // vector of ((At * V-1) * A) matrices
    v_mat_2d v_AtVinv_;            // vector of (At * V-1) matrices
    v_mat_2d v_design_;            // vector of design matrices
    v_mat_2d v_rigorousVariances_; // Precisions of rigorous coordinates

    v_mat_2d v_measMinusComp_;     // vector of measurement matrices
    v_mat_2d v_estimatedStations_; // Coordinate estimates for each block after
                                   // each block adjustment (in isolation)
    v_mat_2d v_originalStations_;  // vector of initial station coordinate
                                   // estimates matrices
    v_mat_2d v_rigorousStations_;  // Coordinate estimates for each block after
                                   // rigorous phased adjustment
    v_mat_2d v_junctionVariances_; // used to carry junction variances between
                                   // all successive blocks
    v_mat_2d
        v_junctionVariancesFwd_; // retains junction variances from forward pass
    v_mat_2d v_junctionEstimatesFwd_; // retains junctions estimates from
                                      // forward pass
    v_mat_2d v_junctionEstimatesRev_; // retains junctions estimates from
                                      // reverse pass
    v_mat_2d v_precAdjMsrsFull_;      // vector of (A * Vx * At) matrices (Vx is
                                      // aposteriori Variance)
    v_mat_2d v_corrections_;          // vector of residuals matrices
    v_mat_2d v_correctionsR_;         // vector of residuals matrices

    // ----------------------------------------------
    // Adjustment functions and variables for staged adjustment

    // Stage adjustment file streams
    void OpenStageFileStreams(const int file_count = 0, ...);
    void ReserveBlockMapRegions(const int file_count = 0, ...);
    void CloseStageFileStreams();

    vmat_file_map normals_map_;
    vmat_file_map normalsR_map_;
    vmat_file_map AtVinv_map_;
    vmat_file_map design_map_;
    vmat_file_map measMinusComp_map_;
    vmat_file_map estimatedStations_map_;
    vmat_file_map originalStations_map_;
    vmat_file_map rigorousStations_map_;
    vmat_file_map junctionVariances_map_;
    vmat_file_map junctionVariancesFwd_map_;
    vmat_file_map junctionEstimatesFwd_map_;
    vmat_file_map junctionEstimatesRev_map_;
    vmat_file_map rigorousVariances_map_;
    vmat_file_map precAdjMsrs_map_;
    vmat_file_map corrections_map_;

    vstring v_stageFileStreams_;
    std::fstream f_normals_;
    std::fstream f_normalsR_;
    std::fstream f_AtVinv_;
    std::fstream f_design_;
    std::fstream f_measMinusComp_;
    std::fstream f_estimatedStations_;
    std::fstream f_originalStations_;
    std::fstream f_rigorousStations_;
    std::fstream f_junctionVariances_;
    std::fstream f_junctionVariancesFwd_;
    std::fstream f_junctionEstimatesFwd_;
    std::fstream f_junctionEstimatesRev_;
    std::fstream f_rigorousVariances_;
    std::fstream f_precAdjMsrs_;
    std::fstream f_corrections_;

    // queue to handle notification of messages for each iteration
    concurrent_queue<UINT32> iterationQueue_;

    // ----------------------------------------------
    // Adjustment matrices for multi-threaded phased adjustment
    // These matrices are used for reverse and combine threads
    v_mat_2d v_designR_;   // vector of design matrices
    v_mat_2d v_AtVinvR_;   // vector of (At * V-1) matrices
    v_mat_2d v_normalsRC_; // vector of ((At * V-1) * A) matrices

    v_mat_2d v_measMinusCompR_;     // vector of measurement matrices
    v_mat_2d v_estimatedStationsR_; // Coordinate estimates for each block after
                                    // each block adjustment (in isolation)
    v_mat_2d v_junctionVariancesR_; // used to carry junction variances between
                                    // all successive blocks

    vUINT32
        v_blockStationsR; // Stations in the current block (used for printing);

    // Database management
    v_msr_database_id_map v_msr_db_map_;
    it_vdbid_t _it_dbid;
    bool databaseIDsLoaded_;

    void LoadDatabaseId();

    // flag to tell if users have cancelled running dynajust
    std::atomic<bool> isCancelled_;
};

} // namespace networkadjust
} // namespace dynadjust

#endif /* DNAADJUST_H_ */
