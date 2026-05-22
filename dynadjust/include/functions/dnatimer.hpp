//============================================================================
// Name         : dnatimer.hpp
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
// Description  : DynAdjust CPU timer replacement for boost::timer
//============================================================================

#ifndef DNATIMER_H_
#define DNATIMER_H_

#include <chrono>
#include <iomanip>
#include <sstream>
#include <string>

namespace dynadjust {

// High-precision timer class to replace boost::timer::cpu_timer
class cpu_timer {
   public:
    struct cpu_times {
        std::chrono::nanoseconds wall;
        std::chrono::nanoseconds user;
        std::chrono::nanoseconds system;
    };

    cpu_timer() { start(); }

    void start() { start_time_ = std::chrono::high_resolution_clock::now(); }

    void resume() { start(); }

    void stop() {
        // For compatibility, but no-op since we calculate elapsed on demand
    }

    cpu_times elapsed() const {
        auto end_time = std::chrono::high_resolution_clock::now();
        auto wall_duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time_);
        return {wall_duration, wall_duration, wall_duration};  // For simplicity, user and system = wall
    }

    std::string format(int places = 6) const {
        auto times = elapsed();
        double wall_seconds = times.wall.count() / 1e9;

        std::ostringstream oss;
        oss << std::fixed << std::setprecision(places);
        oss << wall_seconds << "s wall";
        return oss.str();
    }

   private:
    std::chrono::high_resolution_clock::time_point start_time_;
};

// Format an elapsed duration (in seconds) for display.
//   < 1 second   -> "X.XXXms"
//   1..60 second -> "X.XXXs"
//   >= 60 second -> "hh:mm:ss"
inline std::string FormatElapsedTime(double seconds) {
    std::ostringstream oss;
    if (seconds >= 60.0) {
        long total_seconds = static_cast<long>(seconds);
        long hours = total_seconds / 3600;
        long minutes = (total_seconds % 3600) / 60;
        long secs = total_seconds % 60;
        oss << std::setfill('0') << std::setw(2) << hours << ':' << std::setw(2) << minutes << ':' << std::setw(2)
            << secs;
    } else if (seconds >= 1.0) {
        oss << std::fixed << std::setprecision(3) << seconds << 's';
    } else {
        oss << std::fixed << std::setprecision(3) << (seconds * 1000.0) << "ms";
    }
    return oss.str();
}

}  // namespace dynadjust

#endif  // DNATIMER_H_
