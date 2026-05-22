//============================================================================
// Name         : dnaadjust_json_printer.hpp
// Author       : Dale Roberts <dale.o.roberts@gmail.com>
// Copyright    : Copyright 2026 Geoscience Australia
//
//                Licensed under the Apache License, Version 2.0 (the "License");
//                you may not use this file except in compliance with the License.
//                You may obtain a copy of the License at
//
//                http://www.apache.org/licenses/LICENSE-2.0
//
//                Unless required by applicable law or agreed to in writing, software
//                distributed under the License is distributed on an "AS IS" BASIS,
//                WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//                See the License for the specific language governing permissions and
//                limitations under the License.
//
// Description  : JSONL side-channel printer emitted alongside the default
//                DynAdjust text reports.  Subclass of DynAdjustPrinter that
//                overrides the virtual hooks to write one JSON object per
//                line for each adjusted record, using the same key names as
//                the JSONL import parser so output can round-trip.
//============================================================================

#ifndef DNAADJUST_JSON_PRINTER_H_
#define DNAADJUST_JSON_PRINTER_H_

#if defined(_MSC_VER)
#if defined(LIST_INCLUDES_ON_BUILD)
#pragma message("  " __FILE__)
#endif
#endif

#include <fstream>
#include <mutex>

#include <include/config/dnaexports.hpp>

#include <dynadjust/dnaadjust/dnaadjust_printer.hpp>

namespace dynadjust {
namespace networkadjust {

class DNAADJUST_API DynAdjustJsonPrinter : public DynAdjustPrinter {
 public:
  explicit DynAdjustJsonPrinter(dna_adjust& adjust_instance);
  ~DynAdjustJsonPrinter() override;

  // Open each JSONL output stream whose text sibling is enabled for this
  // run and write the per-file header line.  Must be called after the
  // printer is constructed but before any Print* method is invoked.
  void OpenStreams();

  // Flush and close any open JSONL output streams.  Called from the
  // destructor so callers do not need to invoke it explicitly.
  void CloseStreams();

  // Pairs an output stream with its serialising mutex.  One instance per
  // report file.  Public so helpers in the .cpp translation unit (which
  // is where nlohmann::json is in scope) can hold references — the
  // members themselves stay private.  Concurrent adjustment threads
  // writing to different sibling files don't serialise on a single lock,
  // and the critical section is only the stream write (the JSON dump
  // happens outside).  Non-copyable and non-movable because std::mutex is.
  struct JsonStream {
    std::ofstream stream;
    std::mutex mu;

    JsonStream() = default;
    JsonStream(const JsonStream&) = delete;
    JsonStream& operator=(const JsonStream&) = delete;
    JsonStream(JsonStream&&) = delete;
    JsonStream& operator=(JsonStream&&) = delete;
  };

 protected:
  void OnAdjustedStation(ReportKind kind,
                         const it_vstn_t& stn_it,
                         const matrix_2d* estimates,
                         const matrix_2d* variances,
                         const UINT32& mat_idx,
                         const AdjustedStationContext& ctx) override;
  void OnAdjustedMeasurement(const it_vmsr_t& it_msr) override;
  void OnPositionalUncertainty(const it_vstn_t& stn_it,
                               const matrix_2d* variances,
                               const UINT32& mat_idx,
                               const PositionalUncertaintyContext& ctx) override;
  void OnStationCorrection(const UINT32& block,
                           const it_vstn_t& stn_it,
                           const matrix_2d* estimates,
                           const UINT32& mat_idx,
                           const StationCorrectionContext& ctx) override;
  void OnM2SRecord(const it_vstn_t& stn_it,
                   MsrTally& stn_msr_tally) override;
  void OnStatistics() override;

 private:
  JsonStream adj_;
  JsonStream xyz_;
  JsonStream apu_;
  JsonStream cor_;
  JsonStream m2s_;
};

}  // namespace networkadjust
}  // namespace dynadjust

#endif  // DNAADJUST_JSON_PRINTER_H_
