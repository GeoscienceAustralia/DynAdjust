//============================================================================
// Name         : threading_init.hpp
// Author       : Dale Roberts <dale.o.roberts@gmail.com>
// Contributors : 
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
// Description  : DynAdjust implementation
//============================================================================

#pragma once
/// \cond
#include <algorithm>
#include <cstdlib>
/// \endcond

#define STRINGIFY_HELPER(x) #x
#define STRINGIFY(x) STRINGIFY_HELPER(x)

/// \cond
#if defined(_OPENMP)
#include <omp.h>
#endif

#if defined(USE_MKL) || defined(__MKL__)
#include <mkl.h>

#elif defined(OPENBLAS_VERSION) || defined(__OPENBLAS_CONFIG_H) || defined(OPENBLAS_LOPT_H)

#ifndef USE_OPENBLAS
#define USE_OPENBLAS
#endif

#include <cblas.h>
#include <openblas_config.h>

extern "C" { void openblas_set_num_threads(int); }

#elif defined(__APPLE__)
#include <Accelerate/Accelerate.h>
#endif
/// \endcond

inline int init_linear_algebra_threads(int requested_threads = 0) {
    int n = requested_threads;

#if defined(_OPENMP)
    if (n <= 0) {
        if (const char* env = std::getenv("OMP_NUM_THREADS"); env && *env) {
            int v = std::atoi(env);
            if (v > 0) n = v;
        }
    }
    if (n <= 0) n = std::max(1, omp_get_max_threads());
#elif defined(__APPLE__)
    if (n == 1) {
        BLASSetThreading(BLAS_THREADING_SINGLE_THREADED);
    } else {
        BLASSetThreading(BLAS_THREADING_MULTI_THREADED);
    }
#else
    if (n <= 0) n = 1;
#endif

#if !defined(USE_MKL)
    if (n <= 0) {
        if (const char* env = std::getenv("OPENBLAS_NUM_THREADS"); env && *env) {
            int v = std::atoi(env);
            if (v > 0) n = v;
        }
    }
#endif

#if defined(_OPENMP)
    omp_set_dynamic(0);
#if _OPENMP >= 201511
    omp_set_max_active_levels(1);
#else
    omp_set_nested(0);
#endif
    omp_set_num_threads(n);
#endif

#ifdef USE_MKL
    mkl_set_dynamic(0);
    mkl_set_num_threads(n);
#if defined(MKL_THREAD_LOCAL)
    mkl_set_num_threads_local(n);
#endif
#elif defined(USE_OPENBLAS)
    openblas_set_num_threads(n);
#endif

    return n;
}
