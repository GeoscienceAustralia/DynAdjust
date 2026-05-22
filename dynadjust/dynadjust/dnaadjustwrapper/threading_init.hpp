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

#elif !defined(__APPLE__)

#include <cblas.h>

#if defined(OPENBLAS_VERSION) || defined(__OPENBLAS_CONFIG_H) || defined(OPENBLAS_LOPT_H)

#ifndef USE_OPENBLAS
#define USE_OPENBLAS
#endif

#elif defined(__has_include)
#if __has_include(<openblas_config.h>)
#ifndef USE_OPENBLAS
#define USE_OPENBLAS
#endif
#include <openblas_config.h>
#endif
#endif

#ifdef USE_OPENBLAS
extern "C" { void openblas_set_num_threads(int); }
#endif

#elif defined(__APPLE__)
#include <Accelerate/Accelerate.h>
#endif
/// \endcond

inline int positive_env_int(const char* name) {
    const char* env = std::getenv(name);
    if (!env || !*env) return 0;
    int value = std::atoi(env);
    return value > 0 ? value : 0;
}

inline int init_linear_algebra_threads(int requested_threads = 0) {
    int n = requested_threads;

    if (n <= 0) {
        n = positive_env_int("OMP_NUM_THREADS");
#if defined(USE_MKL) || defined(__MKL__)
        if (n <= 0) n = positive_env_int("MKL_NUM_THREADS");
#elif defined(USE_OPENBLAS)
        if (n <= 0) n = positive_env_int("OPENBLAS_NUM_THREADS");
#elif defined(__APPLE__)
        if (n <= 0) n = positive_env_int("VECLIB_MAXIMUM_THREADS");
#endif
    }

    if (n <= 0) return 0;

#if defined(_OPENMP)
    omp_set_dynamic(0);
#if _OPENMP >= 201511
    omp_set_max_active_levels(1);
#else
    omp_set_nested(0);
#endif
    omp_set_num_threads(n);
#endif

#if defined(__APPLE__)
    if (n == 1) {
        BLASSetThreading(BLAS_THREADING_SINGLE_THREADED);
    } else {
        BLASSetThreading(BLAS_THREADING_MULTI_THREADED);
    }
#elif defined(USE_MKL) || defined(__MKL__)
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
