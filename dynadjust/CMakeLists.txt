cmake_minimum_required(VERSION 3.15)
project(dynadjust VERSION 1.2.7)

option(USE_MKL "Use Intel MKL (Linux and Windows)" OFF)
option(BUILD_TESTING "Enable testing" OFF)
option(BUILD_STATIC "Build static binaries in the 'static' directory" OFF)
option(USE_OPENMP "Use OpenMP for parallel processing" OFF)
option(ILP64 "Use ILP64 interface (64-bit matrix indexing) for BLAS and LAPACK" OFF)

# Build optimization options
option(USE_PRECOMPILED_HEADERS "Use precompiled headers to speed up compilation" ON)
option(USE_UNITY_BUILD "Use unity builds to speed up compilation" OFF)

# Enable precompiled headers
set(CMAKE_DISABLE_PRECOMPILE_HEADERS OFF)

# Enable unity builds if requested
if(USE_UNITY_BUILD)
    set(CMAKE_UNITY_BUILD ON)
    set(CMAKE_UNITY_BUILD_BATCH_SIZE 16)
    message(STATUS "Unity builds enabled (batch size: ${CMAKE_UNITY_BUILD_BATCH_SIZE})")
endif()

# Fix Boost 1.70+ issue with CMake 3.14+
if(POLICY CMP0167)
    cmake_policy(SET CMP0167 NEW)
endif()

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

message("")
message("=====================================================================")
message("Configuring DynAdjust ${PROJECT_VERSION} ${CMAKE_BUILD_TYPE} Build")
message("=====================================================================")
message("")

if (WIN32)
    message(STATUS "Toolchain file: [${CMAKE_TOOLCHAIN_FILE}]")
    message(STATUS "ENV VCPKG_ROOT: [$ENV{VCPKG_ROOT}]")
endif()

if(BUILD_TESTING)
    message(STATUS "Testing enabled (BUILD_TESTING=ON)")
    include(CTest)
    enable_testing()
    set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
endif()

# Set output directories for all platforms to use ../bin/
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/../bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/../bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/../bin)

# For multi-config generators (Visual Studio, Xcode)
foreach(OUTPUTCONFIG ${CMAKE_CONFIGURATION_TYPES})
    string(TOUPPER ${OUTPUTCONFIG} OUTPUTCONFIG_UPPER)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_${OUTPUTCONFIG_UPPER} ${CMAKE_SOURCE_DIR}/../bin)
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY_${OUTPUTCONFIG_UPPER} ${CMAKE_SOURCE_DIR}/../bin)
    set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_${OUTPUTCONFIG_UPPER} ${CMAKE_SOURCE_DIR}/../bin)
endforeach()

if (BUILD_STATIC)
    message(STATUS "Static building enabled (BUILD_STATIC=ON)")

    # Build internal libraries as static instead of shared
    set(DNA_LIBRARY_TYPE STATIC)

    # For maximum portability, statically link libstdc++ and libgcc
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static-libstdc++ -static-libgcc")
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -static-libstdc++ -static-libgcc")

    # Include static build optimizations
    include(${CMAKE_CURRENT_LIST_DIR}/cmake/StaticBuildOptimizations.cmake)
    
    # Enable link-time optimization for static builds
    if(CMAKE_BUILD_TYPE STREQUAL "Release" OR CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
        include(CheckIPOSupported)
        check_ipo_supported(RESULT IPO_SUPPORTED OUTPUT IPO_OUTPUT)
        if(IPO_SUPPORTED)
            message(STATUS "Link-time optimization (LTO/IPO) is supported and enabled for static builds")
            set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
            
            # Add compiler-specific LTO flags
            if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
                add_compile_options(-flto=auto -ffunction-sections -fdata-sections)
                if(NOT APPLE)
                    add_link_options(-flto=auto -Wl,--gc-sections)
                else()
                    add_link_options(-flto=auto -Wl,-dead_strip)
                endif()
            elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
                add_compile_options(-flto=thin -ffunction-sections -fdata-sections)
                if(NOT APPLE)
                    add_link_options(-flto=thin -Wl,--gc-sections)
                else()
                    add_link_options(-flto=thin -Wl,-dead_strip)
                endif()
            elseif(MSVC)
                add_compile_options(/GL)
                add_link_options(/LTCG)
            endif()
        else()
            message(STATUS "Link-time optimization (LTO/IPO) is not supported: ${IPO_OUTPUT}")
        endif()
    else()
        message(STATUS "Link-time optimization disabled for non-Release builds")
    endif()
    
    # Function to print static linking information for each target
    function(print_static_linking_info TARGET_NAME PROGRAM_NAME)
        message(STATUS "")
        message(STATUS "Static linking for ${PROGRAM_NAME}:")
        message(STATUS "  Target: ${TARGET_NAME}")
        message(STATUS "  Libraries being linked:")
        foreach(lib ${STATIC_LIBS_LIST})
            get_filename_component(LIB_NAME ${lib} NAME)
            message(STATUS "    - ${LIB_NAME}")
        endforeach()
    endfunction()
    
    if (WIN32)
        message(STATUS "Static building on Windows requires vcpkg")
        set(VCPKG_CRT_LINKAGE dynamic)
        set(VCPKG_LIBRARY_LINKAGE static)
        set(VCPKG_TARGET_ARCHITECTURE x64-windows)
    endif()
else()
    message(STATUS "Static building disabled (BUILD_STATIC=OFF)")
    set(DNA_LIBRARY_TYPE SHARED)
endif()

if(UNIX)
    set(DNA_PROGRAM_PREFIX "dna")
endif()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}")

# ----------------------------------------------------------------------------
# C++ standard
# ----------------------------------------------------------------------------

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
message(STATUS "C++ standard: ${CMAKE_CXX_STANDARD}")

# ----------------------------------------------------------------------------
# Compiler flags
# ----------------------------------------------------------------------------

# Suppress Boost warnings
add_compile_definitions(
    BOOST_ALLOW_DEPRECATED_HEADERS
    BOOST_BIND_GLOBAL_PLACEHOLDERS
    BOOST_DISABLE_PRAGMA_MESSAGE
)

if(BUILD_TESTING)
    # Suppress all warnings during testing
    if(MSVC)
        add_compile_options(/W0)
    else()
        add_compile_options(-w)
    endif()
endif()

if(APPLE)
    add_compile_options(-Wno-unguarded-availability-new)
    add_compile_options(-Wno-deprecated-declarations)
    add_compile_options(-Wno-nullability-completeness)
endif()

if (WIN32)
    add_compile_options(-DNOMINMAX)
endif()

# ----------------------------------------------------------------------------
# XercesC
# ----------------------------------------------------------------------------

find_package(XercesC REQUIRED)
message(STATUS "Found XercesC (${XercesC_INCLUDE_DIRS} ${XercesC_LIBRARIES})")

# ----------------------------------------------------------------------------
# XSD
# ----------------------------------------------------------------------------

if(WIN32)
    set(LIBXSD_VERSION "4.2.0")
    set(LIBXSD_DIR "${CMAKE_BINARY_DIR}/libxsd-${LIBXSD_VERSION}")
    set(LIBXSD_TAR "${CMAKE_BINARY_DIR}/libxsd-${LIBXSD_VERSION}.tar.gz")
    
    if(NOT EXISTS ${LIBXSD_DIR})
        message(STATUS "Downloading libxsd ${LIBXSD_VERSION} for Xerces 3.2+ compatibility...")
        
        execute_process(
            COMMAND powershell -Command "
                try {
                    Invoke-WebRequest -Uri 'http://www.codesynthesis.com/download/xsd/4.2/libxsd-${LIBXSD_VERSION}.tar.gz' -OutFile '${LIBXSD_TAR}' -UseBasicParsing
                    exit 0
                } catch {
                    exit 1
                }"
            RESULT_VARIABLE DOWNLOAD_RESULT
            OUTPUT_QUIET
            ERROR_QUIET
        )
        
        if(EXISTS ${LIBXSD_TAR})
            message(STATUS "Extracting libxsd...")
            execute_process(
                COMMAND ${CMAKE_COMMAND} -E tar xzf ${LIBXSD_TAR}
                WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
                RESULT_VARIABLE EXTRACT_RESULT
            )
            if(EXTRACT_RESULT)
                message(FATAL_ERROR "Failed to extract libxsd. Please delete ${LIBXSD_TAR} and re-run CMake.")
            endif()
        else()
            message(FATAL_ERROR 
                "\nFailed to download libxsd automatically.\n"
                "Please manually:\n"
                "1. Download: http://www.codesynthesis.com/download/xsd/4.2/libxsd-${LIBXSD_VERSION}.tar.gz\n"
                "2. Extract to: ${CMAKE_BINARY_DIR}\n"
                "3. Re-run CMake\n"
                "\nOr run this PowerShell command:\n"
                "Invoke-WebRequest -Uri 'http://www.codesynthesis.com/download/xsd/4.2/libxsd-${LIBXSD_VERSION}.tar.gz' -OutFile '${LIBXSD_TAR}'")
        endif()
    endif()
    
    # Set XSD include directory to the downloaded libxsd
    set(XSD_INCLUDE_DIRS "${LIBXSD_DIR}")
    set(XSD_INCLUDE_DIR "${LIBXSD_DIR}")
    set(XSD_ROOT "${LIBXSD_DIR}")
    
    message(STATUS "Using libxsd ${LIBXSD_VERSION} from: ${XSD_INCLUDE_DIR}")

else()
    # For non-Windows, check if XSD is available or download libxsd
    find_path(XSD_INCLUDE_DIR xsd/cxx/parser/elements.hxx
        PATHS /usr/include /usr/local/include
    )
    
    if(XSD_INCLUDE_DIR)
        set(XSD_INCLUDE_DIRS ${XSD_INCLUDE_DIR})
        message(STATUS "Found XSD include: ${XSD_INCLUDE_DIRS}")
    endif()
endif()

add_compile_definitions(XSD_CXX11) # Tell XSD to use C++11 and later

# Only call find_package(XSD) if we haven't already set up libxsd
if(NOT XSD_INCLUDE_DIR)
    find_package(XSD REQUIRED)
    message(STATUS "Found XSD (${XSD_INCLUDE_DIR})")
else()
    message(STATUS "Using pre-configured XSD: ${XSD_INCLUDE_DIR}")
endif()

include_directories(${XSD_INCLUDE_DIR} ${XSD_INCLUDE_DIR}/xsd)

# ----------------------------------------------------------------------------
# Boost
# ----------------------------------------------------------------------------

# Force CMake to use FindBoost module instead of BoostConfig
set(Boost_NO_BOOST_CMAKE ON)
set(Boost_NO_CMAKE ON)
find_package(Boost REQUIRED
    COMPONENTS program_options filesystem timer chrono thread)
message(STATUS "Found Boost version: ${Boost_VERSION}")

# ----------------------------------------------------------------------------
# Set up libraries
# ----------------------------------------------------------------------------

set(DNA_LIBRARIES
    ${Boost_LIBRARIES}
    ${XercesC_LIBRARIES}
)

# ----------------------------------------------------------------------------
# OpenMP
# ----------------------------------------------------------------------------

if(USE_OPENMP)
    find_package(OpenMP)

    if(OpenMP_CXX_FOUND)
        # Check for libomp on macOS
        if(APPLE)
            if(NOT "${OpenMP_CXX_LIBRARIES}" MATCHES ".*libomp.*")
                message(WARNING "libomp not found (install with `brew install libomp`)")
                message(STATUS "OpenMP disabled")
            else()
                message(STATUS "Enabled OpenMP on macOS using libomp")
                add_compile_options(-Xpreprocessor -fopenmp)
                add_link_options(-L/opt/homebrew/lib)
                list(APPEND DNA_LIBRARIES ${OpenMP_CXX_LIBRARIES})
            endif()
        elseif(MSVC)
            message(STATUS "Enabled OpenMP for MSVC")
            add_compile_options(${OpenMP_CXX_FLAGS})
            list(APPEND DNA_LIBRARIES ${OpenMP_CXX_LIBRARIES})
        else()
            # Linux and other Unix systems
            message(STATUS "Enabled OpenMP")
            add_compile_options(${OpenMP_CXX_FLAGS})
            add_link_options(${OpenMP_CXX_FLAGS})
            list(APPEND DNA_LIBRARIES ${OpenMP_CXX_LIBRARIES})
        endif()
    else()
        message(WARNING "OpenMP not found. Building without parallel processing support.")
    endif()
else()
    message(STATUS "OpenMP disabled (enable with -DUSE_OPENMP=ON)")
endif()

# ----------------------------------------------------------------------------
# Math libraries
# ----------------------------------------------------------------------------

# ILP64 interface for BLAS and LAPACK
if(ILP64)
    add_compile_options("-DUSE_ILP64")
    message(STATUS "Using ILP64 interface for BLAS and LAPACK")
else()
    message(STATUS "Using LP64 interface for BLAS and LAPACK")
endif()

if(APPLE)
    message(STATUS "Looking for BLAS and LAPACK on MacOS")

    find_library(ACCELERATE_LIB Accelerate)
    if(NOT ACCELERATE_LIB)
        message(FATAL_ERROR "Cannot find Apple Accelerate framework.")
    else()
        message(STATUS "Found Accelerate framework for fast linear algebra")
        add_compile_options("-DACCELERATE_NEW_LAPACK")
    endif()
    list(APPEND DNA_LIBRARIES ${ACCELERATE_LIB})

else()
    # Windows and Linux/Unix handling
    if(WIN32)
        message(STATUS "Looking for BLAS and LAPACK on Windows")
    else()
        message(STATUS "Looking for BLAS and LAPACK on Linux/Unix")
        if(BUILD_STATIC)
            set(BLA_STATIC ON)
        endif()
    endif()

    if(USE_MKL)
        # Use Intel MKL's native CMake config
        add_compile_options("-DUSE_MKL")

        # Explicitly control the MKL interface at compile time
        if(NOT ILP64)
            # Undefine MKL_ILP64 to prevent OneAPI environment from forcing ILP64
            # Use CMAKE_CXX_FLAGS to ensure this is applied early
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -UMKL_ILP64")
            set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -UMKL_ILP64")
            message(STATUS "LP64 interface selected: added -UMKL_ILP64 to compiler flags")
        else()
            message(STATUS "ILP64 interface selected")
        endif()

        # Use FindBLAS to locate Intel MKL with explicit interface selection
        # Note: MKLConfig.cmake in OneAPI 2025.x has issues with interface selection,
        # so we use FindBLAS which allows explicit LP64/ILP64 selection via BLA_VENDOR
        message(STATUS "Using FindBLAS with BLA_VENDOR set to Intel MKL")

        if(ILP64)
            set(BLA_VENDOR "Intel10_64ilp")
        else()
            set(BLA_VENDOR "Intel10_64lp")
        endif()

        find_package(BLAS MODULE REQUIRED)
        find_package(LAPACK MODULE)

        # Verify we actually got MKL
        if(NOT "${BLAS_LIBRARIES}" MATCHES "mkl")
            message(FATAL_ERROR "USE_MKL=ON but Intel MKL was not found. Found: ${BLAS_LIBRARIES}")
        endif()

        message(STATUS "Found Intel MKL via FindBLAS: ${BLAS_LIBRARIES}")
        list(APPEND DNA_LIBRARIES ${BLAS_LIBRARIES})
        if(LAPACK_FOUND)
            list(APPEND DNA_LIBRARIES ${LAPACK_LIBRARIES})
        endif()

        # MKL provides mkl_cblas.h, not cblas.h - no need to search for cblas.h
        
    else()
        # Use standard BLAS and LAPACK search
        find_package(BLAS MODULE REQUIRED)
        find_package(LAPACK MODULE)
        
        message(STATUS "Found BLAS: ${BLAS_LIBRARIES}")
        if(LAPACK_FOUND)
            message(STATUS "Found LAPACK: ${LAPACK_LIBRARIES}")
        endif()
        
        # Add libraries and linker flags
        if(WIN32)
            add_link_options(${BLAS_LINKER_FLAGS})
            if(LAPACK_FOUND)
                add_link_options(${LAPACK_LINKER_FLAGS})
            endif()
        endif()
        
        list(APPEND DNA_LIBRARIES ${BLAS_LIBRARIES})
        if(LAPACK_FOUND)
            list(APPEND DNA_LIBRARIES ${LAPACK_LIBRARIES})
        endif()
        
        # When statically linking OpenBLAS, we need to link against gfortran
        if(BUILD_STATIC AND NOT APPLE AND NOT WIN32)
            string(TOLOWER "${BLAS_LIBRARIES}" BLAS_LIBRARIES_LOWER_CHECK)
            if(BLAS_LIBRARIES_LOWER_CHECK MATCHES "openblas")
                # Find gfortran library
                find_library(GFORTRAN_LIBRARY
                    NAMES gfortran
                    PATHS
                        /usr/lib/gcc/x86_64-linux-gnu/13
                        /usr/lib/gcc/x86_64-linux-gnu/12
                        /usr/lib/gcc/x86_64-linux-gnu/11
                        /usr/lib/gcc/x86_64-linux-gnu/10
                        /usr/lib/gcc/x86_64-linux-gnu/9
                        /usr/lib/gcc/x86_64-linux-gnu/8
                        /usr/lib/gcc/x86_64-linux-gnu/7
                        /usr/lib/gcc/x86_64-linux-gnu
                        /usr/lib/x86_64-linux-gnu
                        /usr/lib64
                        /usr/lib
                        /usr/local/lib
                    PATH_SUFFIXES gcc
                )
                
                if(GFORTRAN_LIBRARY)
                    list(APPEND DNA_LIBRARIES ${GFORTRAN_LIBRARY})
                    message(STATUS "Found gfortran for static OpenBLAS linking: ${GFORTRAN_LIBRARY}")
                else()
                    message(WARNING "gfortran library not found - static OpenBLAS linking may fail")
                    # Try adding it anyway in case it's in the default search path
                    list(APPEND DNA_LIBRARIES gfortran)
                endif()
            endif()
        endif()
        
        # Check if MKL was found by FindBLAS
        string(TOLOWER "${BLAS_LIBRARIES}" BLAS_LIBRARIES_LOWER)
        if(BLAS_LIBRARIES_LOWER MATCHES "mkl")
            message(FATAL_ERROR "Intel MKL detected via FindBLAS but USE_MKL=OFF. If you want to use MKL then set USE_MKL=ON.")
        else()
            # Try to find cblas.h header for non-MKL BLAS implementations (OpenBLAS, etc.)
            set(BLAS_INCLUDE_SEARCH_PATHS
                ${CMAKE_PREFIX_PATH}/include
                ${CMAKE_PREFIX_PATH}/include/openblas
            )
            
            if(WIN32)
                list(APPEND BLAS_INCLUDE_SEARCH_PATHS
                    ${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}/include
                    ${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}/include/openblas
                    C:/vcpkg/installed/x64-windows/include
                    C:/vcpkg/installed/x64-windows/include/openblas
                )
            else()
                list(APPEND BLAS_INCLUDE_SEARCH_PATHS
                    /usr/include
                    /usr/local/include
                    /usr/include/openblas
                    /usr/local/include/openblas
                    # Architecture-specific paths for Debian/Ubuntu
                    /usr/include/x86_64-linux-gnu
                    /usr/include/aarch64-linux-gnu
                    /usr/include/i386-linux-gnu
                )
            endif()
            
            find_path(BLAS_INCLUDE_DIRS 
                NAMES cblas.h
                PATHS ${BLAS_INCLUDE_SEARCH_PATHS}
                REQUIRED
            )
            
            if(BLAS_INCLUDE_DIRS)
                message(STATUS "Found cblas.h at: ${BLAS_INCLUDE_DIRS}")
                include_directories(${BLAS_INCLUDE_DIRS})
            endif()
        endif()
    endif()

endif()

message(STATUS "DynAdjust library dependencies: ${DNA_LIBRARIES}")

include_directories(
    ${CMAKE_SOURCE_DIR}
    ${Boost_INCLUDE_DIRS}
    ${XercesC_INCLUDE_DIRS}
    ${XSD_INCLUDE_DIRS}
)

# ----------------------------------------------------------------------------
# Static build  (i.e., no dependencies on shared libraries)
# ----------------------------------------------------------------------------

if(BUILD_STATIC)
    message("")
    message("=====================================================================")
    message("Static build configuration")
    message("=====================================================================")
    message("")

    add_compile_definitions(
        BUILD_STATIC
        XERCES_STATIC_LIBRARY
    )

    if(MSVC)
        set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Release>:Release>")
        if(NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
            add_compile_options(/O2)
        endif()
    else()
        if(CMAKE_BUILD_TYPE STREQUAL "Release" OR CMAKE_BUILD_TYPE STREQUAL "")
            add_compile_options(-O3)
        elseif(CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
            add_compile_options(-O2)
        endif()
    endif()

    if(BLAS_LIBRARIES)
        message(STATUS "BLAS_LIBRARIES=${BLAS_LIBRARIES}")
    endif()
    if(LAPACK_LIBRARIES)
        message(STATUS "LAPACK_LIBRARIES=${LAPACK_LIBRARIES}")
    endif()

    if (UNIX)
        # Find the static library equivalents
        set(STATIC_LIBS_LIST "")
        set(DNA_LIB_LIST ${DNA_LIBRARIES})
        list(REMOVE_DUPLICATES DNA_LIB_LIST)
        foreach(lib ${DNA_LIB_LIST})
            get_filename_component(libdir ${lib} DIRECTORY)
            get_filename_component(libname ${lib} NAME_WE)
            string(REGEX REPLACE "-l" "" libname ${libname})
            string(REGEX REPLACE "^Boost::" "boost_" libname ${libname})
            string(REGEX REPLACE "^XercesC::XercesC" "xerces-c" libname ${libname})
            string(REGEX REPLACE "^lib" "" libname ${libname})
            if(libname)
                if(NOT (libname STREQUAL "m" OR 
                        libname STREQUAL "dl" OR 
                        libname MATCHES "^-W*" OR 
                        libname STREQUAL "Accelerate"))
                    message(CHECK_START "Finding static '${lib}'")
                    find_library(STATIC_LIB_${libname}
                        NAMES lib${libname}.a lib${libname}_nonshared.a
                        PATHS ${libdir})
                    if (${STATIC_LIB_${libname}} MATCHES "NOTFOUND")
                        message(CHECK_FAIL "Not Found: 'lib${libname}.a'")
                    else()
                        message(CHECK_PASS "Found")
                        list(APPEND STATIC_LIBS_LIST ${STATIC_LIB_${libname}})
                    endif()
                endif()
            endif()
        endforeach()

        message("")
        message(STATUS "STATIC_LIBS_LIST=${STATIC_LIBS_LIST}")

        # Find static versions of xerces-c dependencies (curl and ICU)
        find_library(CURL_STATIC_LIB NAMES libcurl.a HINTS
            ENV LIBRARY_PATH /g/data/ms9/.envbin/lib /usr/lib64 /usr/lib)
        find_library(ICUUC_STATIC_LIB NAMES libicuuc.a HINTS
            ENV LIBRARY_PATH /g/data/ms9/.envbin/lib /usr/lib64 /usr/lib)
        find_library(ICUI18N_STATIC_LIB NAMES libicui18n.a HINTS
            ENV LIBRARY_PATH /g/data/ms9/.envbin/lib /usr/lib64 /usr/lib)
        find_library(ICUDATA_STATIC_LIB NAMES libicudata.a HINTS
            ENV LIBRARY_PATH /g/data/ms9/.envbin/lib /usr/lib64 /usr/lib)
        find_library(SSL_STATIC_LIB NAMES libssl.a HINTS
            ENV LIBRARY_PATH /g/data/ms9/.envbin/lib /usr/lib64 /usr/lib)
        find_library(CRYPTO_STATIC_LIB NAMES libcrypto.a HINTS
            ENV LIBRARY_PATH /g/data/ms9/.envbin/lib /usr/lib64 /usr/lib)
        find_library(Z_STATIC_LIB NAMES libz.a HINTS
            ENV LIBRARY_PATH /g/data/ms9/.envbin/lib /usr/lib64 /usr/lib)

        # Build the list of xerces-c dependencies, preferring static
        set(XERCES_DEPS "")

        # ICU libraries (order matters: i18n depends on uc, both depend on data)
        if(ICUI18N_STATIC_LIB AND ICUUC_STATIC_LIB AND ICUDATA_STATIC_LIB)
            message(STATUS "Found static ICU: ${ICUUC_STATIC_LIB}")
            list(APPEND XERCES_DEPS ${ICUI18N_STATIC_LIB} ${ICUUC_STATIC_LIB} ${ICUDATA_STATIC_LIB})
        else()
            message(STATUS "Static ICU not found, using dynamic")
            list(APPEND XERCES_DEPS -licui18n -licuuc -licudata)
        endif()

        # curl and its dependencies
        if(CURL_STATIC_LIB)
            message(STATUS "Found static curl: ${CURL_STATIC_LIB}")
            list(APPEND XERCES_DEPS ${CURL_STATIC_LIB})
            # curl's dependencies
            if(SSL_STATIC_LIB AND CRYPTO_STATIC_LIB)
                message(STATUS "Found static OpenSSL")
                list(APPEND XERCES_DEPS ${SSL_STATIC_LIB} ${CRYPTO_STATIC_LIB})
            else()
                list(APPEND XERCES_DEPS -lssl -lcrypto)
            endif()
            if(Z_STATIC_LIB)
                message(STATUS "Found static zlib: ${Z_STATIC_LIB}")
                list(APPEND XERCES_DEPS ${Z_STATIC_LIB})
            else()
                list(APPEND XERCES_DEPS -lz)
            endif()
            # Other curl deps that are usually only available as shared
            list(APPEND XERCES_DEPS -lnghttp2 -lidn2 -lssh -lpsl -lbrotlidec -lbrotlicommon)
        else()
            message(STATUS "Static curl not found, using dynamic")
            list(APPEND XERCES_DEPS -lcurl)
        endif()

        # Set DNA_LIBRARIES with all dependencies
        # System libs at the end
        set(DNA_LIBRARIES ${STATIC_LIBS_LIST} ${XERCES_DEPS} -lpthread -lrt -lm -ldl)

        # Print sizes of static libraries
        message(STATUS "")
        message(STATUS "Static library sizes:")
        message(STATUS "---------------------")
        foreach(lib ${STATIC_LIBS_LIST})
            if(EXISTS ${lib})
                file(SIZE ${lib} LIB_SIZE)
                math(EXPR LIB_SIZE_MB "${LIB_SIZE} / 1048576")
                math(EXPR LIB_SIZE_KB "(${LIB_SIZE} % 1048576) / 1024")
                get_filename_component(LIB_NAME ${lib} NAME)
                message(STATUS "  ${LIB_NAME}: ${LIB_SIZE_MB}.${LIB_SIZE_KB} MB (${LIB_SIZE} bytes)")
            else()
                message(STATUS "  ${lib}: File not found")
            endif()
        endforeach()
        
        # Calculate total size
        set(TOTAL_SIZE 0)
        foreach(lib ${STATIC_LIBS_LIST})
            if(EXISTS ${lib})
                file(SIZE ${lib} LIB_SIZE)
                math(EXPR TOTAL_SIZE "${TOTAL_SIZE} + ${LIB_SIZE}")
            endif()
        endforeach()
        math(EXPR TOTAL_SIZE_MB "${TOTAL_SIZE} / 1048576")
        message(STATUS "---------------------")
        message(STATUS "Total static libraries size: ${TOTAL_SIZE_MB} MB")
        message(STATUS "")
    endif()

    file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/static)

    if(APPLE)
        message(STATUS "Using macOS-specific settings for static binaries")
        # On macOS, we'll create dynamically linked binaries with maximum compatibility
        set(CMAKE_MACOSX_RPATH ON)

        # Use a more recent deployment target to avoid LAPACK availability issues
        # as *new* Accelerate framework available from macOS 13.3+
        set(CMAKE_OSX_DEPLOYMENT_TARGET "15.0" CACHE STRING "Minimum macOS version" FORCE)
        add_compile_options("-mmacosx-version-min=15.0")
    endif()

    message("")
    message("---------------------------------------------------------------------")
    message("")
endif()

# ----------------------------------------------------------------------------
# Subdirectories
# ----------------------------------------------------------------------------


message(STATUS "Configuring import...")
add_subdirectory(dynadjust/dnaimportwrapper)

message(STATUS "Configuring reftran...")
add_subdirectory(dynadjust/dnareftranwrapper)

message(STATUS "Configuring segment...")
add_subdirectory(dynadjust/dnasegmentwrapper)

message(STATUS "Configuring geoid...")
add_subdirectory(dynadjust/dnageoidwrapper)

message(STATUS "Configuring plot...")
add_subdirectory(dynadjust/dnaplotwrapper)

message(STATUS "Configuring dynadjust...")
add_subdirectory(dynadjust/dynadjust)

message(STATUS "Configuring adjust...")
add_subdirectory(dynadjust/dnaadjustwrapper)

message(STATUS "Configuring diff...")
# add_subdirectory(dynadjust/dnadiff)  # Not present in version 1.2.7

# ----------------------------------------------------------------------------
# Static build info
# ----------------------------------------------------------------------------


if(BUILD_STATIC)
    message("")
    message(" **")
    message(" * Static binaries will be available in `static` subdirectory *")
    message("                                                             **")
    message("")
endif()


# ----------------------------------------------------------------------------
# Unified Installation
# ----------------------------------------------------------------------------

# Set default install prefix if not specified
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    # Install to repo root bin directory (DynAdjust/bin)
    set(CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}/.." CACHE PATH "Installation Directory" FORCE)
endif()

message(STATUS "Install prefix: ${CMAKE_INSTALL_PREFIX}")

# Install executables based on build type
if(BUILD_STATIC)
    # Install static executables (they already output with correct names)
    # Note: static targets not yet implemented, commented out
    # install(TARGETS
    #     dnaimportwrapper_static dnareftranwrapper_static dnasegmentwrapper_static
    #     dnageoidwrapper_static dnaplotwrapper_static dnaadjustwrapper_static
    #     dynadjust_static
    #     RUNTIME DESTINATION bin)

    # For now, install regular executables with static linking
    install(TARGETS
        dnaimportwrapper dnareftranwrapper dnasegmentwrapper
        dnageoidwrapper dnaplotwrapper dnaadjustwrapper dynadjust
        RUNTIME DESTINATION bin)
else()
    # Install wrapper executables
    install(TARGETS
        dnaimportwrapper dnareftranwrapper dnasegmentwrapper
        dnageoidwrapper dnaplotwrapper dnaadjustwrapper dynadjust
        RUNTIME DESTINATION bin)

    # Install libraries (only for non-static builds)
    install(TARGETS
        dnaimport dnareftran dnasegment dnageoid dnaplot dnaadjust
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION bin
        ARCHIVE DESTINATION bin)
endif()

# ----------------------------------------------------------------------------
# Testing
# ----------------------------------------------------------------------------

if (BUILD_TESTING)
    message("=====================================================================")
    message("Testing enabled (BUILD_TESTING=ON)")
    message("=====================================================================")
    message("")
    message(STATUS "Setting up tests...")


    # test execution of dynadjust binaries
    # 1. gnss network (simultaneous, phased-block 1, phased, staged)
    add_test (NAME import-gnss-network COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> -n gnss gnss-network.stn gnss-network.msr --export-dna --export-xml --export-asl --export-aml --export-map --output-msr-to-stn --test-integrity -r GDA2020 --flag-unused-stations) 
    add_test (NAME geoid-gnss-network COMMAND $<TARGET_FILE:${DNAGEOID_TARGET}> gnss -g gnss-network-geoid.gsb --export-dna-geo)
    add_test (NAME reftran-gnss-network COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss --export-dna --export-xml)
    add_test (NAME adjust-gnss-network COMMAND $<TARGET_FILE:${DNAADJUST_TARGET}> gnss --output-adj-msr --scale-normals-to-unity)

    # Check results with dnadiff
    add_test(NAME test-gnss-network
        COMMAND $<TARGET_FILE:${DNADIFF_TARGET}> gnss.simult.adj gnss.simult.adj.expected --skip-headers 52 -t 0.001)

    add_test (NAME import-gnss-network-similar COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> -n gnss_similar gnss-network.stn gnss-network.msr -r itrf2008 --override-input-ref-frame --search-similar-gnss-msr --quiet) 
    add_test (NAME import-gnss-network-exclude COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> -n gnss_excl gnss-network.stn gnss-network.msr --exclude-stns-assoc-msrs "BEEC,261000380,356000780,349800490" --split) 
    add_test (NAME import-gnss-network-scaled COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> -n gnss_scaled gnss-network.stn gnss-network.msr --baseline-scalar-file gnss-network.scalars -r GDA94) 
    add_test (NAME import-gnss-network-block1 COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> -n gnss_b1 gnss-network.stn gnss-network.msr --export-dna --export-xml --single-xml-file --export-asl --export-aml --export-map --output-msr-to-stn --test-integrity -r GDA94 --flag-unused-stations) 
    add_test (NAME segment-gnss-network-block1-a COMMAND $<TARGET_FILE:${DNASEGMENT_TARGET}> gnss_b1 --min 90 --max 90 --starting-stns 220700210)
    add_test (NAME segment-gnss-network-block1-b COMMAND $<TARGET_FILE:${DNASEGMENT_TARGET}> gnss_b1 --min 110 --max 110 --net-file)
    add_test (NAME segment-gnss-network-block1-c COMMAND $<TARGET_FILE:${DNASEGMENT_TARGET}> gnss_b1 --min 80 --max 200 --net-file --starting-stns 211300940 --verbose 3)
    add_test (NAME adjust-gnss-network-block1 COMMAND $<TARGET_FILE:${DNAADJUST_TARGET}> gnss_b1 --block1-phased --output-adj-msr --output-adj-gnss-units 1)
    add_test (NAME adjust-gnss-network-phased COMMAND $<TARGET_FILE:${DNAADJUST_TARGET}> gnss_b1 --phased --output-adj-msr --output-adj-gnss-units 2 --output-iter-adj-msr)
    add_test (NAME adjust-gnss-network-verb COMMAND $<TARGET_FILE:${DNAADJUST_TARGET}> gnss_b1 --verbose 5)
    add_test (NAME adjust-gnss-network-staged-a COMMAND $<TARGET_FILE:${DNAADJUST_TARGET}> gnss_b1 --staged-adjustment --create-stage-files --output-adj-msr --output-adj-gnss-units 2 --sort-adj-msr-field 4)
    add_test (NAME adjust-gnss-network-staged-b COMMAND $<TARGET_FILE:${DNAADJUST_TARGET}> gnss_b1 --staged-adjustment --purge-stage-files --output-adj-msr --output-adj-gnss-units 3 --output-stn-blocks --output-msr-blocks --sort-adj-msr-field 5)

    # 2. urban network (phased-sequential)
    add_test (NAME import-urban-network COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> -n urban urban-network.stn urban-network.msr --flag-unused-stations)
    add_test (NAME geoid-urban-network COMMAND $<TARGET_FILE:${DNAGEOID_TARGET}> urban -g urban-network-geoid.gsb --export-dna-geo)
    add_test (NAME segment-urban-network COMMAND $<TARGET_FILE:${DNASEGMENT_TARGET}> urban --min 50 --max 150 --test-integrity)
    add_test (NAME adjust-urban-network-verbose COMMAND $<TARGET_FILE:${DNAADJUST_TARGET}> urban --verbose 3)
    add_test (NAME adjust-urban-network COMMAND $<TARGET_FILE:${DNAADJUST_TARGET}> urban --output-adj-msr --phased --stn-corrections --export-sinex-file --export-xml-stn-file --export-dna-stn-file --output-pos-uncertainty --export-dna-msr --export-xml-msr)
    add_test (NAME plot-urban-network-01 COMMAND $<TARGET_FILE:${DNAPLOT_TARGET}> urban --phased --label-sta --correction-arrows --label-corr --compute-corrections --scale-arrows 10.5 --error-ellipse --positional-uncertainty --scale-ellipse-c 10.5)
    add_test (NAME plot-urban-network-02 COMMAND $<TARGET_FILE:${DNAPLOT_TARGET}> urban --phased --label-sta --label-constraints --correction-arrows --label-corr --compute-corrections --scale-arrows 10.5 --error-ellipse --positional-uncertainty --scale-ellipse-c 10.5 --block-number 2 --alternate-name)
    add_test (NAME plot-urban-seg-stn COMMAND $<TARGET_FILE:${DNAPLOT_TARGET}> urban --supress-pdf --graph-stn)
    add_test (NAME plot-urban-seg-msr COMMAND $<TARGET_FILE:${DNAPLOT_TARGET}> urban --supress-pdf --graph-msr)

    # 3. urban network (transform to GDA2020, phased-concurrent)
    add_test (NAME import-urban-network-thread COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> -n urban_mt urban-network.stn urban-network.msr)
    add_test (NAME reftran-urban-network-thread COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> urban_mt -r gda2020)
    add_test (NAME geoid-urban-network-thread COMMAND $<TARGET_FILE:${DNAGEOID_TARGET}> urban_mt -g urban-network-geoid.gsb)
    add_test (NAME segment-urban-network-thread COMMAND $<TARGET_FILE:${DNASEGMENT_TARGET}> urban_mt --min 50 --max 85 --verb 3)
    add_test (NAME adjust-urban-network-thread-01 COMMAND $<TARGET_FILE:${DNAADJUST_TARGET}> urban_mt --output-adj-msr --multi  --free-stn-sd 4.0 --fixed-stn-sd 0.000001 --max-iterations 20 --output-tstat-adj-msr --sort-adj-msr-field 2 --sort-stn-orig-order --stn-coord-types PLHhENz --angular-stn-type 1 --angular-msr-type 1 --precision-stn-linear 3 --precision-msr-linear 3 --precision-stn-angular 4 --precision-msr-angular 4 --output-pos-uncertainty --output-all-covariances --output-corrections-file)
    add_test (NAME plot-urban-network-thread COMMAND $<TARGET_FILE:${DNAPLOT_TARGET}> urban_mt --phased --label-sta --label-font 16 --msr-line-w 0.5 --map-projection 3)
    add_test (NAME adjust-urban-network-thread-02 COMMAND $<TARGET_FILE:${DNAADJUST_TARGET}> urban_mt --output-adj-msr --multi)

    # 4. urban network (phased-staged)
    add_test (NAME import-urban-network-stage COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> -n urban_st urban-network.stn urban-network.msr --flag-unused-stations)
    add_test (NAME geoid-urban-network-stage COMMAND $<TARGET_FILE:${DNAGEOID_TARGET}> urban_st -g urban-network-geoid.gsb --export-dna-geo)
    add_test (NAME segment-urban-network-stage COMMAND $<TARGET_FILE:${DNASEGMENT_TARGET}> urban_st --min 90 --max 90)
    add_test (NAME adjust-urban-network-stage COMMAND $<TARGET_FILE:${DNAADJUST_TARGET}> urban_st --phased --staged-adjustment --create-stage-files --output-adj-msr --export-sinex-file --output-pos-uncertainty --export-xml-stn-file --export-xml-msr-file --export-dna-stn-file --export-dna-msr --output-iter-adj-stn --output-iter-adj-stat --output-iter-adj-msr --output-iter-cmp-msr --stn-corrections --output-corrections-file)

    # test all frame labels
    add_test (NAME imp-frame-misc-01 COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> -n impframe-01 urban-network.stn urban-network.msr -r itrf1988 -e 03.12.1988)
    add_test (NAME ref-frame-misc-01 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> impframe-01 --verb 6 --plate-model-option 1 -b PB2002_plates.dig -m PB2002_poles.dat)
    add_test (NAME imp-frame-misc-02 COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> -n impframe-02 urban-network.stn urban-network.msr -r itrf1988 -e 03.12.1988)
    add_test (NAME ref-frame-misc-02 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> impframe-02 --verb 6 --plate-model-option 1 -b PB2002_plates.dig -m PB2002_poles.dat)
    add_test (NAME imp-frame-misc-03 COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> -n impframe-03 urban-network.stn urban-network.msr -r itrf1989 -e 01.09.1989)
    add_test (NAME ref-frame-misc-03 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> impframe-03 --verb 6 --plate-model-option 1 -b PB2002_plates.dig -m PB2002_poles.dat)
    add_test (NAME imp-frame-misc-04 COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> -n impframe-04 urban-network.stn urban-network.msr -r itrf1990 -e 01.10.1990)
    add_test (NAME ref-frame-misc-04 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> impframe-04 --verb 6 --plate-model-option 1 -b PB2002_plates.dig -m PB2002_poles.dat)
    add_test (NAME imp-frame-misc-05 COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> -n impframe-05 urban-network.stn urban-network.msr -r itrf1991 -e 01.01.1991)
    add_test (NAME ref-frame-misc-05 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> impframe-05 --verb 6 --plate-model-option 1 -b PB2002_plates.dig -m PB2002_poles.dat)
    add_test (NAME imp-frame-misc-06 COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> -n impframe-06 urban-network.stn urban-network.msr -r itrf1992 -e 01.02.1992)
    add_test (NAME ref-frame-misc-06 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> impframe-06 --verb 6 --plate-model-option 1 -b PB2002_plates.dig -m PB2002_poles.dat)
    add_test (NAME imp-frame-misc-07 COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> -n impframe-07 urban-network.stn urban-network.msr -r itrf1993 -e 01.03.1993)
    add_test (NAME ref-frame-misc-07 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> impframe-07 --verb 6 --plate-model-option 1 -b PB2002_plates.dig -m PB2002_poles.dat)
    add_test (NAME imp-frame-misc-08 COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> -n impframe-08 urban-network.stn urban-network.msr -r itrf1994 -e 01.04.1994)
    add_test (NAME ref-frame-misc-08 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> impframe-08 --verb 6 --plate-model-option 1 -b PB2002_plates.dig -m PB2002_poles.dat)
    add_test (NAME imp-frame-misc-09 COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> -n impframe-09 urban-network.stn urban-network.msr -r itrf1996 -e 01.06.1996)
    add_test (NAME ref-frame-misc-09 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> impframe-09 --verb 6 --plate-model-option 1 -b PB2002_plates.dig -m PB2002_poles.dat)
    add_test (NAME imp-frame-misc-10 COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> -n impframe-10 urban-network.stn urban-network.msr -r itrf1997 -e 01.07.1997)
    add_test (NAME ref-frame-misc-10 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> impframe-10 --verb 6 --plate-model-option 1 -b PB2002_plates.dig -m PB2002_poles.dat)
    add_test (NAME imp-frame-misc-11 COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> -n impframe-11 urban-network.stn urban-network.msr -r itrf2000 -e 22.2.2000)
    add_test (NAME ref-frame-misc-11 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> impframe-11 --verb 6 --plate-model-option 1 -b PB2002_plates.dig -m PB2002_poles.dat)
    add_test (NAME imp-frame-misc-12 COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> -n impframe-12 urban-network.stn urban-network.msr -r itrf2005 -e 5.5.2005)
    add_test (NAME ref-frame-misc-12 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> impframe-12 --verb 6 --plate-model-option 1 -b PB2002_plates.dig -m PB2002_poles.dat)
    add_test (NAME imp-frame-misc-13 COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> -n impframe-13 urban-network.stn urban-network.msr -r itrf2008 -e 8.8.2008)
    add_test (NAME ref-frame-misc-13 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> impframe-13 --verb 6 --plate-model-option 1 -b PB2002_plates.dig -m PB2002_poles.dat)
    add_test (NAME imp-frame-misc-14 COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> -n impframe-14 urban-network.stn urban-network.msr -r itrf2014 -e 2020)
    add_test (NAME ref-frame-misc-14 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> impframe-14 --verb 6 --plate-model-option 1 -b PB2002_plates.dig -m PB2002_poles.dat)
    add_test (NAME imp-frame-misc-15 COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> -n impframe-15 urban-network.stn urban-network.msr -r itrf2020 -e today)
    add_test (NAME ref-frame-misc-15 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> impframe-15 --verb 6 --plate-model-option 1 -b PB2002_plates.dig -m PB2002_poles.dat)
    add_test (NAME imp-frame-misc-16 COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> -n impframe-16 urban-network.stn urban-network.msr -r gda94)
    add_test (NAME ref-frame-misc-16 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> impframe-16 --verb 6 --plate-model-option 1 -b PB2002_plates.dig -m PB2002_poles.dat)
    add_test (NAME imp-frame-misc-17 COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> -n impframe-17 urban-network.stn urban-network.msr -r gda2020)
    add_test (NAME ref-frame-misc-17 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> impframe-17 --verb 6 --plate-model-option 1 -b PB2002_plates.dig -m PB2002_poles.dat)
    add_test (NAME imp-frame-misc-18 COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> -n impframe urban-network.stn urban-network.msr -r wgs84)
    add_test (NAME imp-frame-misc-19 COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> -n impframe urban-network.stn urban-network.msr -r "wgs 84")
    add_test (NAME imp-frame-misc-20 COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> -n impframe urban-network.stn urban-network.msr -r "wgs84 (ensemble)")
    add_test (NAME imp-frame-misc-21 COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> -n impframe urban-network.stn urban-network.msr -r "wgs84 (G730)")
    add_test (NAME imp-frame-misc-22 COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> -n impframe urban-network.stn urban-network.msr -r "wgs84 (G873)")
    add_test (NAME imp-frame-misc-23 COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> -n impframe urban-network.stn urban-network.msr -r "wgs84 (G1150)")
    add_test (NAME imp-frame-misc-24 COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> -n impframe urban-network.stn urban-network.msr -r "wgs84 (G1674)")
    add_test (NAME imp-frame-misc-25 COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> -n impframe urban-network.stn urban-network.msr -r "wgs84 (G1762)")
    add_test (NAME imp-frame-misc-26 COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> -n impframe urban-network.stn urban-network.msr -r "wgs84 (G2139)")
    add_test (NAME imp-frame-misc-27 COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> -n impframe urban-network.stn urban-network.msr -r "wgs84 (transit)")

    # 5. dynadjust project mode    
    add_test (NAME dynadjust-name-01 COMMAND $<TARGET_FILE:${DYNADJUST_TARGET}> -n sample --import --geoid --reftran --segment --adjust)
    set_tests_properties(dynadjust-name-01 PROPERTIES ENVIRONMENT "PATH=${CMAKE_BINARY_DIR}/bin:$ENV{PATH}")

    # 6. urban data manipulation
    add_test (NAME imp-urb00-dna-data COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> -n urban_dna urban-network.stn urban-network.msr --export-xml --export-asl --export-aml --export-map --output-msr-to-stn)
    add_test (NAME imp-urb01-xml-data COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> -n urban_xml ./urban_dnastn.xml ./urban_dnamsr.xml --include-msr-types gxyemsabdpqr --exclude-msr-types vz)
    add_test (NAME imp-urb02-xml-data COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> -n urban_xml ./urban_dnastn.xml ./urban_dnamsr.xml --search-nearby-stn --nearby-stn-buffer 0.4)
    add_test (NAME imp-urb03-xml-data COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> -n urban_xml ./urban_dnastn.xml ./urban_dnamsr.xml --search-similar-msr --ignore-similar-msr --flag-unused-stations)
    add_test (NAME imp-urb04-xml-data COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> -n urban_xml ./urban_dnastn.xml ./urban_dnamsr.xml --v-scale 2.5 --p-scale 5.0 --l-scale 5.0 --h-scale 10.0)
    add_test (NAME imp-urb05-xml-data COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> -n urban_xml ./urban_dnastn.xml ./urban_dnamsr.xml --include-stns-assoc-msrs "4000,1,1004,1042,1046,2,2215,4000,5,9004,1002,1003,1034,1039,1050,13,2101,2102,2105,2106,2109,2118,2119,2122,2123,2124,2125,2206,2214")
    add_test (NAME imp-urb06-xml-data COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> -n urban_xml ./urban_dnastn.xml ./urban_dnamsr.xml --export-dna --export-asl --export-aml --export-map --output-msr-to-stn --flag-unused-stations)
    add_test (NAME seg-urb00-xml-data COMMAND $<TARGET_FILE:${DNASEGMENT_TARGET}> urban_xml --min 50 --max 50)
    add_test (NAME imp-urb07-xml-data COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> -n urban_xml --import-block-stn-msr 2)
    add_test (NAME imp-urb08-xml-data COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> -n urban_xml ./urban_dnastn.xml ./urban_dnamsr.xml --geo-file urban.geo)
    add_test (NAME plt-urb01-xml-data COMMAND $<TARGET_FILE:${DNAPLOT_TARGET}> urban_xml --centre-station 1004 --area-radius 5000 --plot-msr-types gxyasvz --msr-line-width 0.5)
    add_test (NAME plt-urb02-xml-data COMMAND $<TARGET_FILE:${DNAPLOT_TARGET}> urban_xml --centre-latitude -37.480489013 --centre-longitude 144.573579814 --area-radius 3000 --omit-measurements)
    add_test (NAME imp-psu01-dna-data COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> -n pseudo pseudo.stn --flag-unused-stations)
    add_test (NAME geo-psu00-dna-data COMMAND $<TARGET_FILE:${DNAGEOID_TARGET}> pseudo -g gnss-network-geoid.gsb --export-dna-geo-file)
    add_test (NAME imp-psu02-dna-data COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> -n pseudo pseudo.stn pseudo-driver-file.msr --simulate --geo-file pseudo.geo --export-xml)
    add_test (NAME imp-psu03-xml-data COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> -n pseudo pseudostn.xml pseudomsr.xml --export-asl --export-aml --export-map)
    add_test (NAME geo-psu01-xml-data COMMAND $<TARGET_FILE:${DNAGEOID_TARGET}> pseudo -g gnss-network-geoid.gsb --convert --interpolation-method 0)
    add_test (NAME adj-psu00-xml-data COMMAND $<TARGET_FILE:${DNAADJUST_TARGET}> pseudo --output-adj-msr --free-stn-sd 5.0 --fixed-stn-sd 0.000005 --max-iterations 15 --output-tstat-adj-msr --sort-adj-msr-field 7 --sort-stn-orig-order --stn-coord-types XYZPLHhENz --angular-stn-type 1 --angular-msr-type 1 --precision-stn-linear 3 --precision-msr-linear 3 --precision-stn-angular 4 --precision-msr-angular 4 --output-pos-uncertainty --output-all-covariances --output-corrections-file)
    add_test (NAME imp-gnss0-xml-data COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> -n gnss_xml gnssstn.xml gnssmsr.xml --bounding-box "-36.3000,145.4500,-36.4500,146.4500" --get-msrs-transcending-box --split-gnss-cluster-msrs --prefer-single-x-as-g)
    add_test (NAME imp-gnss1-xml-data COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> -n gnss_xml gnssstn.xml gnssmsr.xml --include-stns-assoc-msrs "BEEC,BNLA,EURA,HOTH,MNSF,MYRT,320500750,385900240" --get-msrs-transcending-box --split-gnss-cluster-msrs --prefer-single-x-as-g -r itrf2014 --override-input-ref-frame)

    # 7. geoid data manipulation
    add_test (NAME geo-interpolate-01 COMMAND $<TARGET_FILE:${DNAGEOID_TARGET}> -g gnss-network-geoid.gsb --interactive --latitude -36.3348253617 --longitude 145.5741006771)
    add_test (NAME geo-interpolate-02 COMMAND $<TARGET_FILE:${DNAGEOID_TARGET}> -g gnss-network-geoid.gsb --interactive --latitude -36.2704890597 --longitude 146.4809871537 --interpolation-method 0)
    add_test (NAME geo-interpolate-03 COMMAND $<TARGET_FILE:${DNAGEOID_TARGET}> -g gnss-network-geoid.gsb --summary)
    add_test (NAME geo-interpolate-04 COMMAND $<TARGET_FILE:${DNAGEOID_TARGET}> --version)
    add_test (NAME geo-interpolate-05 COMMAND $<TARGET_FILE:${DNAGEOID_TARGET}> -g gnss-network-geoid.gsb -t geoid-points-dat.txt)
    add_test (NAME geo-interpolate-06 COMMAND $<TARGET_FILE:${DNAGEOID_TARGET}> -g gnss-network-geoid.gsb -t geoid-points-csv.csv)    
    add_test (NAME geo-interpolate-07 COMMAND $<TARGET_FILE:${DNAGEOID_TARGET}> -g gnss-network-geoid.gsb -t geoid-points-dat.txt --decimal-degrees --direction 0 --export-dna-geo-file)
    add_test (NAME geo-interpolate-08 COMMAND $<TARGET_FILE:${DNAGEOID_TARGET}> -g gnss-network-geoid.gsb -t geoid-points-csv.csv --decimal-degrees --direction 0 --export-dna-geo-file)
    add_test (NAME geo-ntv2-create-01 COMMAND $<TARGET_FILE:${DNAGEOID_TARGET}> -d ausgeoid09_gda94_v1.01_clip_1x1.dat -c -g ausgeoid_clip_default.gsb)
    add_test (NAME geo-ntv2-create-02 COMMAND $<TARGET_FILE:${DNAGEOID_TARGET}> -d ausgeoid09_gda94_v1.01_clip_1x1.dat -c -g ausgeoid_clip_1.0.1.0.gsb --grid-shift radians --grid-version 1.0.1.0 --system-fr ___GDA94 --system-to ___AHD71 --sub-grid-n 1D-grid --creation 21.04.2021 --update 22.04.2021)
    add_test (NAME geo-ntv2-export-01 COMMAND $<TARGET_FILE:${DNAGEOID_TARGET}> -g gnss-network-geoid.gsb --export-ntv2-asc)
    add_test (NAME geo-ntv2-export-02 COMMAND $<TARGET_FILE:${DNAGEOID_TARGET}> -g gnss-network-geoid.gsb --export-ntv2-asc --grid-shift-type radians)
    add_test (NAME geo-ntv2-export-03 COMMAND $<TARGET_FILE:${DNAGEOID_TARGET}> -g ./gnss-network-geoid.gsb.asc --export-ntv2-gsb)
    add_test (NAME geo-ntv2-export-04 COMMAND $<TARGET_FILE:${DNAGEOID_TARGET}> -g ./gnss-network-geoid.gsb.asc --export-ntv2-gsb --grid-shift-type radians)
    add_test (NAME imp-outside COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> -n outside urban-network.stn urban-network.msr)
    add_test (NAME geo-outside-01 COMMAND $<TARGET_FILE:${DNAGEOID_TARGET}> outside -g gnss-network-geoid.gsb --convert --verbose 1)
    add_test (NAME geo-outside-02 COMMAND $<TARGET_FILE:${DNAGEOID_TARGET}> outside -g ./gnss-network-geoid.gsb.asc --convert --verbose 2)
    add_test (NAME geo-inside-01 COMMAND $<TARGET_FILE:${DNAGEOID_TARGET}> gnss -g ./gnss-network-geoid.gsb.asc --convert --export-dna-geo)
    add_test (NAME geo-inside-02 COMMAND $<TARGET_FILE:${DNAGEOID_TARGET}> -g ./gnss-network-geoid.gsb.asc -t geoid-points-dat.txt)
    add_test (NAME geo-inside-03 COMMAND $<TARGET_FILE:${DNAGEOID_TARGET}> -g ./gnss-network-geoid.gsb.asc -t geoid-points-dat.txt --interpolation-method 0)

    # Check results with dnadiff
    add_test (NAME test-urban-phased-network COMMAND $<TARGET_FILE:${DNADIFF_TARGET}> urban.phased.adj urban.phased.adj.expected --skip-headers 62 -t 0.001)
    add_test (NAME test-urban-thread-network COMMAND $<TARGET_FILE:${DNADIFF_TARGET}> urban_mt.phased-mt.adj urban_mt.phased-mt.adj.expected -t 0.01 -v)

    # 8. gnss reference frame transformations
    add_test (NAME ref-gnss01-network COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2014 -e 01.01.2020)
    add_test (NAME ref-gnss02-network COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1988 -e 03.12.1995)
    add_test (NAME ref-gnss03-network COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1989)
    add_test (NAME ref-gnss04-network COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1990)
    add_test (NAME ref-gnss05-network COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1991)
    add_test (NAME ref-gnss06-network COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2005 -e 5.5.2005 --export-xml --single-xml-file --verb 3)
    add_test (NAME ref-gnss07-network COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1992 --verb 3)
    add_test (NAME ref-gnss08-network COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r gda94 --export-xml --single-xml-file --verb 3)
    add_test (NAME ref-gnss09-network COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1988)
    add_test (NAME ref-gnss10-network COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1993)
    add_test (NAME ref-gnss11-network COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1994)
    add_test (NAME ref-gnss12-network COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1996)
    add_test (NAME ref-gnss13-network COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1997)
    add_test (NAME ref-gnss14-network COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2020)

    # gda2020 
    add_test (NAME gda2020-param-01 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r gda2020)
    add_test (NAME gda2020-param-02 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r gda94)
    add_test (NAME gda2020-param-03 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r gda2020)
    add_test (NAME gda2020-param-04 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1988)
    add_test (NAME gda2020-param-05 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r gda2020)
    add_test (NAME gda2020-param-06 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1989)
    add_test (NAME gda2020-param-07 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r gda2020)
    add_test (NAME gda2020-param-08 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1990)
    add_test (NAME gda2020-param-09 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r gda2020)
    add_test (NAME gda2020-param-10 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1991)
    add_test (NAME gda2020-param-11 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r gda2020)
    add_test (NAME gda2020-param-12 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1992)
    add_test (NAME gda2020-param-13 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r gda2020)
    add_test (NAME gda2020-param-14 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1993)
    add_test (NAME gda2020-param-15 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r gda2020)
    add_test (NAME gda2020-param-16 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1994)
    add_test (NAME gda2020-param-17 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r gda2020)
    add_test (NAME gda2020-param-18 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1996)
    add_test (NAME gda2020-param-19 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r gda2020)
    add_test (NAME gda2020-param-20 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1997)
    add_test (NAME gda2020-param-21 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r gda2020)
    add_test (NAME gda2020-param-22 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2000)
    add_test (NAME gda2020-param-23 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r gda2020)
    add_test (NAME gda2020-param-24 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2005)
    add_test (NAME gda2020-param-25 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r gda2020)
    add_test (NAME gda2020-param-26 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2008)
    add_test (NAME gda2020-param-27 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r gda2020)
    add_test (NAME gda2020-param-28 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2020)
    # gda94 
    add_test (NAME gda94-param-01 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r gda94)
    add_test (NAME gda94-param-02 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r gda2020)
    add_test (NAME gda94-param-03 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r gda94)
    add_test (NAME gda94-param-04 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1988)
    add_test (NAME gda94-param-05 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r gda94)
    add_test (NAME gda94-param-06 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1989)
    add_test (NAME gda94-param-07 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r gda94)
    add_test (NAME gda94-param-08 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1990)
    add_test (NAME gda94-param-09 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r gda94)
    add_test (NAME gda94-param-10 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1991)
    add_test (NAME gda94-param-11 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r gda94)
    add_test (NAME gda94-param-12 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1992)
    add_test (NAME gda94-param-13 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r gda94)
    add_test (NAME gda94-param-14 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1993)
    add_test (NAME gda94-param-15 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r gda94)
    add_test (NAME gda94-param-16 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1994)
    add_test (NAME gda94-param-17 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r gda94)
    add_test (NAME gda94-param-18 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1996)
    add_test (NAME gda94-param-19 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r gda94)
    add_test (NAME gda94-param-20 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1997)
    add_test (NAME gda94-param-21 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r gda94)
    add_test (NAME gda94-param-22 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2000)
    add_test (NAME gda94-param-23 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r gda94)
    add_test (NAME gda94-param-24 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2005)
    add_test (NAME gda94-param-25 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r gda94)
    add_test (NAME gda94-param-26 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2008)
    add_test (NAME gda94-param-27 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r gda94)
    add_test (NAME gda94-param-28 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2020)
    # itrf1988 
    add_test (NAME itrf1988-param-01 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1988)
    add_test (NAME itrf1988-param-02 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r gda94)
    add_test (NAME itrf1988-param-03 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1988)
    add_test (NAME itrf1988-param-04 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r gda2020)
    add_test (NAME itrf1988-param-05 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1988)
    add_test (NAME itrf1988-param-06 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1989)
    add_test (NAME itrf1988-param-07 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1988)
    add_test (NAME itrf1988-param-08 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1990)
    add_test (NAME itrf1988-param-09 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1988)
    add_test (NAME itrf1988-param-10 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1991)
    add_test (NAME itrf1988-param-11 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1988)
    add_test (NAME itrf1988-param-12 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1992)
    add_test (NAME itrf1988-param-13 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1988)
    add_test (NAME itrf1988-param-14 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1993)
    add_test (NAME itrf1988-param-15 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1988)
    add_test (NAME itrf1988-param-16 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1994)
    add_test (NAME itrf1988-param-17 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1988)
    add_test (NAME itrf1988-param-18 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1996)
    add_test (NAME itrf1988-param-19 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1988)
    add_test (NAME itrf1988-param-20 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1997)
    add_test (NAME itrf1988-param-21 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1988)
    add_test (NAME itrf1988-param-22 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2000)
    add_test (NAME itrf1988-param-23 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1988)
    add_test (NAME itrf1988-param-24 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2005)
    add_test (NAME itrf1988-param-25 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1988)
    add_test (NAME itrf1988-param-26 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2008)
    add_test (NAME itrf1988-param-27 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1988)
    add_test (NAME itrf1988-param-28 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2020)
    # itrf1989 
    add_test (NAME itrf1989-param-01 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1989)
    add_test (NAME itrf1989-param-02 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r gda94)
    add_test (NAME itrf1989-param-03 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1989)
    add_test (NAME itrf1989-param-04 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r gda2020)
    add_test (NAME itrf1989-param-05 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1989)
    add_test (NAME itrf1989-param-06 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1988)
    add_test (NAME itrf1989-param-07 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1989)
    add_test (NAME itrf1989-param-08 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1990)
    add_test (NAME itrf1989-param-09 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1989)
    add_test (NAME itrf1989-param-10 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1991)
    add_test (NAME itrf1989-param-11 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1989)
    add_test (NAME itrf1989-param-12 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1992)
    add_test (NAME itrf1989-param-13 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1989)
    add_test (NAME itrf1989-param-14 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1993)
    add_test (NAME itrf1989-param-15 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1989)
    add_test (NAME itrf1989-param-16 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1994)
    add_test (NAME itrf1989-param-17 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1989)
    add_test (NAME itrf1989-param-18 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1996)
    add_test (NAME itrf1989-param-19 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1989)
    add_test (NAME itrf1989-param-20 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1997)
    add_test (NAME itrf1989-param-21 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1989)
    add_test (NAME itrf1989-param-22 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2000)
    add_test (NAME itrf1989-param-23 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1989)
    add_test (NAME itrf1989-param-24 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2005)
    add_test (NAME itrf1989-param-25 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1989)
    add_test (NAME itrf1989-param-26 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2008)
    add_test (NAME itrf1989-param-27 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1989)
    add_test (NAME itrf1989-param-28 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2020)
    # itrf1990 
    add_test (NAME itrf1990-param-01 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1990)
    add_test (NAME itrf1990-param-02 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r gda94)
    add_test (NAME itrf1990-param-03 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1990)
    add_test (NAME itrf1990-param-04 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r gda2020)
    add_test (NAME itrf1990-param-05 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1990)
    add_test (NAME itrf1990-param-06 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1988)
    add_test (NAME itrf1990-param-07 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1990)
    add_test (NAME itrf1990-param-08 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1989)
    add_test (NAME itrf1990-param-09 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1990)
    add_test (NAME itrf1990-param-10 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1991)
    add_test (NAME itrf1990-param-11 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1990)
    add_test (NAME itrf1990-param-12 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1992)
    add_test (NAME itrf1990-param-13 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1990)
    add_test (NAME itrf1990-param-14 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1993)
    add_test (NAME itrf1990-param-15 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1990)
    add_test (NAME itrf1990-param-16 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1994)
    add_test (NAME itrf1990-param-17 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1990)
    add_test (NAME itrf1990-param-18 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1996)
    add_test (NAME itrf1990-param-19 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1990)
    add_test (NAME itrf1990-param-20 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1997)
    add_test (NAME itrf1990-param-21 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1990)
    add_test (NAME itrf1990-param-22 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2000)
    add_test (NAME itrf1990-param-23 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1990)
    add_test (NAME itrf1990-param-24 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2005)
    add_test (NAME itrf1990-param-25 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1990)
    add_test (NAME itrf1990-param-26 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2008)
    add_test (NAME itrf1990-param-27 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1990)
    add_test (NAME itrf1990-param-28 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2020)
    # itrf1991 
    add_test (NAME itrf1991-param-01 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1991)
    add_test (NAME itrf1991-param-02 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r gda94)
    add_test (NAME itrf1991-param-03 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1991)
    add_test (NAME itrf1991-param-04 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r gda2020)
    add_test (NAME itrf1991-param-05 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1991)
    add_test (NAME itrf1991-param-06 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1988)
    add_test (NAME itrf1991-param-07 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1991)
    add_test (NAME itrf1991-param-08 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1989)
    add_test (NAME itrf1991-param-09 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1991)
    add_test (NAME itrf1991-param-10 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1990)
    add_test (NAME itrf1991-param-11 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1991)
    add_test (NAME itrf1991-param-12 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1992)
    add_test (NAME itrf1991-param-13 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1991)
    add_test (NAME itrf1991-param-14 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1993)
    add_test (NAME itrf1991-param-15 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1991)
    add_test (NAME itrf1991-param-16 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1994)
    add_test (NAME itrf1991-param-17 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1991)
    add_test (NAME itrf1991-param-18 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1996)
    add_test (NAME itrf1991-param-19 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1991)
    add_test (NAME itrf1991-param-20 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1997)
    add_test (NAME itrf1991-param-21 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1991)
    add_test (NAME itrf1991-param-22 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2000)
    add_test (NAME itrf1991-param-23 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1991)
    add_test (NAME itrf1991-param-24 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2005)
    add_test (NAME itrf1991-param-25 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1991)
    add_test (NAME itrf1991-param-26 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2008)
    add_test (NAME itrf1991-param-27 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1991)
    add_test (NAME itrf1991-param-28 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2020)
    # itrf1992 
    add_test (NAME itrf1992-param-01 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1992)
    add_test (NAME itrf1992-param-02 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r gda94)
    add_test (NAME itrf1992-param-03 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1992)
    add_test (NAME itrf1992-param-04 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r gda2020)
    add_test (NAME itrf1992-param-05 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1992)
    add_test (NAME itrf1992-param-06 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1988)
    add_test (NAME itrf1992-param-07 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1992)
    add_test (NAME itrf1992-param-08 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1989)
    add_test (NAME itrf1992-param-09 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1992)
    add_test (NAME itrf1992-param-10 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1990)
    add_test (NAME itrf1992-param-11 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1992)
    add_test (NAME itrf1992-param-12 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1991)
    add_test (NAME itrf1992-param-13 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1992)
    add_test (NAME itrf1992-param-14 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1993)
    add_test (NAME itrf1992-param-15 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1992)
    add_test (NAME itrf1992-param-16 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1994)
    add_test (NAME itrf1992-param-17 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1992)
    add_test (NAME itrf1992-param-18 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1996)
    add_test (NAME itrf1992-param-19 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1992)
    add_test (NAME itrf1992-param-20 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1997)
    add_test (NAME itrf1992-param-21 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1992)
    add_test (NAME itrf1992-param-22 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2000)
    add_test (NAME itrf1992-param-23 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1992)
    add_test (NAME itrf1992-param-24 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2005)
    add_test (NAME itrf1992-param-25 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1992)
    add_test (NAME itrf1992-param-26 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2008)
    add_test (NAME itrf1992-param-27 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1992)
    add_test (NAME itrf1992-param-28 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2020)
    # itrf1993 
    add_test (NAME itrf1993-param-01 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1993)
    add_test (NAME itrf1993-param-02 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r gda94)
    add_test (NAME itrf1993-param-03 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1993)
    add_test (NAME itrf1993-param-04 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r gda2020)
    add_test (NAME itrf1993-param-05 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1993)
    add_test (NAME itrf1993-param-06 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1988)
    add_test (NAME itrf1993-param-07 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1993)
    add_test (NAME itrf1993-param-08 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1989)
    add_test (NAME itrf1993-param-09 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1993)
    add_test (NAME itrf1993-param-10 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1990)
    add_test (NAME itrf1993-param-11 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1993)
    add_test (NAME itrf1993-param-12 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1991)
    add_test (NAME itrf1993-param-13 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1993)
    add_test (NAME itrf1993-param-14 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1992)
    add_test (NAME itrf1993-param-15 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1993)
    add_test (NAME itrf1993-param-16 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1994)
    add_test (NAME itrf1993-param-17 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1993)
    add_test (NAME itrf1993-param-18 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1996)
    add_test (NAME itrf1993-param-19 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1993)
    add_test (NAME itrf1993-param-20 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1997)
    add_test (NAME itrf1993-param-21 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1993)
    add_test (NAME itrf1993-param-22 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2000)
    add_test (NAME itrf1993-param-23 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1993)
    add_test (NAME itrf1993-param-24 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2005)
    add_test (NAME itrf1993-param-25 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1993)
    add_test (NAME itrf1993-param-26 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2008)
    add_test (NAME itrf1993-param-27 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1993)
    add_test (NAME itrf1993-param-28 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2020)
    # itrf1994 
    add_test (NAME itrf1994-param-01 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1994)
    add_test (NAME itrf1994-param-02 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r gda94)
    add_test (NAME itrf1994-param-03 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1994)
    add_test (NAME itrf1994-param-04 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r gda2020)
    add_test (NAME itrf1994-param-05 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1994)
    add_test (NAME itrf1994-param-06 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1988)
    add_test (NAME itrf1994-param-07 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1994)
    add_test (NAME itrf1994-param-08 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1989)
    add_test (NAME itrf1994-param-09 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1994)
    add_test (NAME itrf1994-param-10 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1990)
    add_test (NAME itrf1994-param-11 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1994)
    add_test (NAME itrf1994-param-12 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1991)
    add_test (NAME itrf1994-param-13 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1994)
    add_test (NAME itrf1994-param-14 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1992)
    add_test (NAME itrf1994-param-15 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1994)
    add_test (NAME itrf1994-param-16 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1993)
    add_test (NAME itrf1994-param-17 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1994)
    add_test (NAME itrf1994-param-18 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1996)
    add_test (NAME itrf1994-param-19 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1994)
    add_test (NAME itrf1994-param-20 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1997)
    add_test (NAME itrf1994-param-21 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1994)
    add_test (NAME itrf1994-param-22 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2000)
    add_test (NAME itrf1994-param-23 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1994)
    add_test (NAME itrf1994-param-24 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2005)
    add_test (NAME itrf1994-param-25 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1994)
    add_test (NAME itrf1994-param-26 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2008)
    add_test (NAME itrf1994-param-27 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1994)
    add_test (NAME itrf1994-param-28 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2020)
    # itrf1996 
    add_test (NAME itrf1996-param-01 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1996)
    add_test (NAME itrf1996-param-02 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r gda94)
    add_test (NAME itrf1996-param-03 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1996)
    add_test (NAME itrf1996-param-04 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r gda2020)
    add_test (NAME itrf1996-param-05 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1996)
    add_test (NAME itrf1996-param-06 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1988)
    add_test (NAME itrf1996-param-07 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1996)
    add_test (NAME itrf1996-param-08 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1989)
    add_test (NAME itrf1996-param-09 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1996)
    add_test (NAME itrf1996-param-10 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1990)
    add_test (NAME itrf1996-param-11 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1996)
    add_test (NAME itrf1996-param-12 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1991)
    add_test (NAME itrf1996-param-13 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1996)
    add_test (NAME itrf1996-param-14 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1992)
    add_test (NAME itrf1996-param-15 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1996)
    add_test (NAME itrf1996-param-16 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1993)
    add_test (NAME itrf1996-param-17 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1996)
    add_test (NAME itrf1996-param-18 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1994)
    add_test (NAME itrf1996-param-19 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1996)
    add_test (NAME itrf1996-param-20 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1997)
    add_test (NAME itrf1996-param-21 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1996)
    add_test (NAME itrf1996-param-22 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2000)
    add_test (NAME itrf1996-param-23 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1996)
    add_test (NAME itrf1996-param-24 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2005)
    add_test (NAME itrf1996-param-25 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1996)
    add_test (NAME itrf1996-param-26 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2008)
    add_test (NAME itrf1996-param-27 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1996)
    add_test (NAME itrf1996-param-28 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2020)
    # itrf1997 
    add_test (NAME itrf1997-param-01 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1997)
    add_test (NAME itrf1997-param-02 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r gda94)
    add_test (NAME itrf1997-param-03 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1997)
    add_test (NAME itrf1997-param-04 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r gda2020)
    add_test (NAME itrf1997-param-05 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1997)
    add_test (NAME itrf1997-param-06 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1988)
    add_test (NAME itrf1997-param-07 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1997)
    add_test (NAME itrf1997-param-08 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1989)
    add_test (NAME itrf1997-param-09 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1997)
    add_test (NAME itrf1997-param-10 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1990)
    add_test (NAME itrf1997-param-11 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1997)
    add_test (NAME itrf1997-param-12 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1991)
    add_test (NAME itrf1997-param-13 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1997)
    add_test (NAME itrf1997-param-14 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1992)
    add_test (NAME itrf1997-param-15 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1997)
    add_test (NAME itrf1997-param-16 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1993)
    add_test (NAME itrf1997-param-17 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1997)
    add_test (NAME itrf1997-param-18 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1994)
    add_test (NAME itrf1997-param-19 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1997)
    add_test (NAME itrf1997-param-20 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1996)
    add_test (NAME itrf1997-param-21 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1997)
    add_test (NAME itrf1997-param-22 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2000)
    add_test (NAME itrf1997-param-23 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1997)
    add_test (NAME itrf1997-param-24 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2005)
    add_test (NAME itrf1997-param-25 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1997)
    add_test (NAME itrf1997-param-26 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2008)
    add_test (NAME itrf1997-param-27 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1997)
    add_test (NAME itrf1997-param-28 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2020)
    # itrf2000 
    add_test (NAME itrf2000-param-01 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2000)
    add_test (NAME itrf2000-param-02 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r gda94)
    add_test (NAME itrf2000-param-03 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2000)
    add_test (NAME itrf2000-param-04 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r gda2020)
    add_test (NAME itrf2000-param-05 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2000)
    add_test (NAME itrf2000-param-06 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1988)
    add_test (NAME itrf2000-param-07 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2000)
    add_test (NAME itrf2000-param-08 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1989)
    add_test (NAME itrf2000-param-09 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2000)
    add_test (NAME itrf2000-param-10 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1990)
    add_test (NAME itrf2000-param-11 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2000)
    add_test (NAME itrf2000-param-12 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1991)
    add_test (NAME itrf2000-param-13 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2000)
    add_test (NAME itrf2000-param-14 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1992)
    add_test (NAME itrf2000-param-15 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2000)
    add_test (NAME itrf2000-param-16 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1993)
    add_test (NAME itrf2000-param-17 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2000)
    add_test (NAME itrf2000-param-18 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1994)
    add_test (NAME itrf2000-param-19 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2000)
    add_test (NAME itrf2000-param-20 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1996)
    add_test (NAME itrf2000-param-21 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2000)
    add_test (NAME itrf2000-param-22 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1997)
    add_test (NAME itrf2000-param-23 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2000)
    add_test (NAME itrf2000-param-24 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2005)
    add_test (NAME itrf2000-param-25 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2000)
    add_test (NAME itrf2000-param-26 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2008)
    add_test (NAME itrf2000-param-27 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2000)
    add_test (NAME itrf2000-param-28 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2020)
    # itrf2005 
    add_test (NAME itrf2005-param-01 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2005)
    add_test (NAME itrf2005-param-02 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r gda94)
    add_test (NAME itrf2005-param-03 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2005)
    add_test (NAME itrf2005-param-04 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r gda2020)
    add_test (NAME itrf2005-param-05 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2005)
    add_test (NAME itrf2005-param-06 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1988)
    add_test (NAME itrf2005-param-07 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2005)
    add_test (NAME itrf2005-param-08 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1989)
    add_test (NAME itrf2005-param-09 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2005)
    add_test (NAME itrf2005-param-10 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1990)
    add_test (NAME itrf2005-param-11 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2005)
    add_test (NAME itrf2005-param-12 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1991)
    add_test (NAME itrf2005-param-13 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2005)
    add_test (NAME itrf2005-param-14 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1992)
    add_test (NAME itrf2005-param-15 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2005)
    add_test (NAME itrf2005-param-16 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1993)
    add_test (NAME itrf2005-param-17 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2005)
    add_test (NAME itrf2005-param-18 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1994)
    add_test (NAME itrf2005-param-19 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2005)
    add_test (NAME itrf2005-param-20 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1996)
    add_test (NAME itrf2005-param-21 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2005)
    add_test (NAME itrf2005-param-22 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1997)
    add_test (NAME itrf2005-param-23 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2005)
    add_test (NAME itrf2005-param-24 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2000)
    add_test (NAME itrf2005-param-25 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2005)
    add_test (NAME itrf2005-param-26 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2008)
    add_test (NAME itrf2005-param-27 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2005)
    add_test (NAME itrf2005-param-28 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2020)
    # itrf2008 
    add_test (NAME itrf2008-param-01 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2008)
    add_test (NAME itrf2008-param-02 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r gda94)
    add_test (NAME itrf2008-param-03 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2008)
    add_test (NAME itrf2008-param-04 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r gda2020)
    add_test (NAME itrf2008-param-05 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2008)
    add_test (NAME itrf2008-param-06 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1988)
    add_test (NAME itrf2008-param-07 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2008)
    add_test (NAME itrf2008-param-08 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1989)
    add_test (NAME itrf2008-param-09 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2008)
    add_test (NAME itrf2008-param-10 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1990)
    add_test (NAME itrf2008-param-11 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2008)
    add_test (NAME itrf2008-param-12 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1991)
    add_test (NAME itrf2008-param-13 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2008)
    add_test (NAME itrf2008-param-14 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1992)
    add_test (NAME itrf2008-param-15 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2008)
    add_test (NAME itrf2008-param-16 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1993)
    add_test (NAME itrf2008-param-17 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2008)
    add_test (NAME itrf2008-param-18 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1994)
    add_test (NAME itrf2008-param-19 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2008)
    add_test (NAME itrf2008-param-20 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1996)
    add_test (NAME itrf2008-param-21 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2008)
    add_test (NAME itrf2008-param-22 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1997)
    add_test (NAME itrf2008-param-23 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2008)
    add_test (NAME itrf2008-param-24 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2000)
    add_test (NAME itrf2008-param-25 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2008)
    add_test (NAME itrf2008-param-26 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2005)
    add_test (NAME itrf2008-param-27 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2008)
    add_test (NAME itrf2008-param-28 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2020)
    # itrf2020 
    add_test (NAME itrf2020-param-01 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r gda94)
    add_test (NAME itrf2020-param-02 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2020)
    add_test (NAME itrf2020-param-03 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r gda2020)
    add_test (NAME itrf2020-param-04 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2020)
    add_test (NAME itrf2020-param-05 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1988)
    add_test (NAME itrf2020-param-06 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2020)
    add_test (NAME itrf2020-param-07 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1989)
    add_test (NAME itrf2020-param-08 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2020)
    add_test (NAME itrf2020-param-09 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1990)
    add_test (NAME itrf2020-param-10 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2020)
    add_test (NAME itrf2020-param-11 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1991)
    add_test (NAME itrf2020-param-12 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2020)
    add_test (NAME itrf2020-param-13 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1992)
    add_test (NAME itrf2020-param-14 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2020)
    add_test (NAME itrf2020-param-15 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1993)
    add_test (NAME itrf2020-param-16 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2020)
    add_test (NAME itrf2020-param-17 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1994)
    add_test (NAME itrf2020-param-18 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2020)
    add_test (NAME itrf2020-param-19 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1996)
    add_test (NAME itrf2020-param-20 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2020)
    add_test (NAME itrf2020-param-21 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf1997)
    add_test (NAME itrf2020-param-22 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2020)
    add_test (NAME itrf2020-param-23 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2000)
    add_test (NAME itrf2020-param-24 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2020)
    add_test (NAME itrf2020-param-25 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2005)
    add_test (NAME itrf2020-param-26 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2020)
    add_test (NAME itrf2020-param-27 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> gnss -r itrf2008)

    # 9. gnss reference frame transformations with plate motion model
    add_test (NAME ref-itrf-pmm-01 COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> -n itrf05_itrf1997_a -r itrf2005 TEST_ITRF05.SNX --flag-unused-stations)
    add_test (NAME ref-itrf-pmm-02 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> itrf05_itrf1997_a -r itrf1997 -e 25.10.2017 --export-dna --verb 2 --plate-model-option 1 -b PB2002_plates.dig -m PB2002_poles.dat)
    add_test (NAME ref-itrf-pmm-03 COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> -n itrf05_itrf1997_b -r itrf2005 TEST_ITRF05.SNX --flag-unused-stations)
    add_test (NAME ref-itrf-pmm-04 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> itrf05_itrf1997_b -r itrf1997 -e 01.01.2020 --export-dna --verb 3 --plate-model-option 1 -b PB2002_plates.dig -m PB2002_poles.dat)
    add_test (NAME ref-itrf-pmm-05 COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> -n apr -r itrf2008 apr.ITRF2008.04.06.2020.stn apr.ITRF2008.04.06.2020.msr --flag-unused-stations)
    add_test (NAME ref-itrf-pmm-06 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> apr -r itrf2008 -e 01.01.2021 --export-dna --verb 2 --plate-model-option 1 -b MORVEL56_plates.dig -m NNR-MORVEL56_poles.dat)

    #add_test (NAME ref-itrf-pmm-07 COMMAND $<TARGET_FILE:${DNADIFF_TARGET}> apr.ITRF2008.01.01.2021.stn apr.ITRF2008.01.01.2021.stn.expected --skip-headers 7 -t 0.0001)

    # 10. station discontinuities
    add_test (NAME imp-discont-01 COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> -n discont -r itrf2005 TEST_ITRF05.SNX --discontinuity-file disconts20201205.snx --export-discont-file --flag-unused-stations)
    add_test (NAME imp-discont-02 COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> -n discont -r itrf2005 TEST_ITRF05.SNX --discontinuity-file disconts20201205.snx --split-gnss --include-stn "ALIC" --export-dna --flag-unused-stations)
    add_test (NAME imp-discont-03 COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> -n discont -r itrf2014 --override-input-ref-frame TEST_ITRF05.SNX sa-gnss-stn.xml sa-gnss-msr.xml --prefer-single-x-as-g --flag-unused --discontinuity-file disconts20201205.snx --remove-ignored-msr --split-gnss --include-stns-assoc-msrs "ALIC" --export-dna --export-xml)
    add_test (NAME adj-discont-01 COMMAND $<TARGET_FILE:${DNAADJUST_TARGET}> discont --constraint "ALIC,CCC")
    add_test (NAME imp-discont-04 COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> -n discont -r itrf2005 TEST_ITRF05.SNX --export-xml)
    add_test (NAME imp-discont-05 COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> -n discont -r itrf2000 --override-input-ref-frame ./discontstn.xml ./discontmsr.xml --discontinuity-file disconts20201205.snx)

    # 11. measurements to stations sorting
    add_test (NAME imp-m2s-sort-01 COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> -n gnss gnss-network.stn gnss-network.msr --output-msr-to-stn --sort-msr-to-stn-field 1 -r GDA94) 

    # 12. help modules and version
    add_test (NAME imp-help-00 COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> -h)
    add_test (NAME imp-help-01 COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> --help-module standard)
    add_test (NAME imp-help-02 COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> --help-module reference)
    add_test (NAME imp-help-03 COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> --help-module screening)
    add_test (NAME imp-help-04 COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> --help-module gnss)
    add_test (NAME imp-help-05 COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> --help-module simulation)
    add_test (NAME imp-help-06 COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> --help-module output)
    add_test (NAME imp-help-07 COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> --help-module export)
    add_test (NAME imp-help-08 COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> --help-module generic)
    add_test (NAME imp-version-01 COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> --version)

    add_test (NAME geo-help-00 COMMAND $<TARGET_FILE:${DNAGEOID_TARGET}> -h)
    add_test (NAME geo-help-01 COMMAND $<TARGET_FILE:${DNAGEOID_TARGET}> --help-module standard)
    add_test (NAME geo-help-02 COMMAND $<TARGET_FILE:${DNAGEOID_TARGET}> --help-module interpolation)
    add_test (NAME geo-help-03 COMMAND $<TARGET_FILE:${DNAGEOID_TARGET}> --help-module ntv2)
    add_test (NAME geo-help-04 COMMAND $<TARGET_FILE:${DNAGEOID_TARGET}> --help-module interactive)
    add_test (NAME geo-help-05 COMMAND $<TARGET_FILE:${DNAGEOID_TARGET}> --help-module file)
    add_test (NAME geo-help-06 COMMAND $<TARGET_FILE:${DNAGEOID_TARGET}> --help-module export)
    add_test (NAME geo-help-07 COMMAND $<TARGET_FILE:${DNAGEOID_TARGET}> --help-module generic)
    add_test (NAME geo-version-01 COMMAND $<TARGET_FILE:${DNAGEOID_TARGET}> --version)

    add_test (NAME seg-help-00 COMMAND $<TARGET_FILE:${DNASEGMENT_TARGET}> -h)
    add_test (NAME seg-help-01 COMMAND $<TARGET_FILE:${DNASEGMENT_TARGET}> --help-module standard)
    add_test (NAME seg-help-02 COMMAND $<TARGET_FILE:${DNASEGMENT_TARGET}> --help-module configuration)
    add_test (NAME seg-help-03 COMMAND $<TARGET_FILE:${DNASEGMENT_TARGET}> --help-module generic)
    add_test (NAME seg-version-01 COMMAND $<TARGET_FILE:${DNASEGMENT_TARGET}> --version)

    add_test (NAME ref-help-00 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> -h)
    add_test (NAME ref-help-01 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> --help-module standard)
    add_test (NAME ref-help-02 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> --help-module transformation)
    add_test (NAME ref-help-03 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> --help-module export)
    add_test (NAME ref-help-04 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> --help-module generic)
    add_test (NAME ref-version-01 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> --version)

    add_test (NAME adj-help-00 COMMAND $<TARGET_FILE:${DNAADJUST_TARGET}> -h)
    add_test (NAME adj-help-01 COMMAND $<TARGET_FILE:${DNAADJUST_TARGET}> --help-module standard)
    add_test (NAME adj-help-02 COMMAND $<TARGET_FILE:${DNAADJUST_TARGET}> --help-module adjustment)
    add_test (NAME adj-help-03 COMMAND $<TARGET_FILE:${DNAADJUST_TARGET}> --help-module phased)
    add_test (NAME adj-help-04 COMMAND $<TARGET_FILE:${DNAADJUST_TARGET}> --help-module configuration)
    add_test (NAME adj-help-05 COMMAND $<TARGET_FILE:${DNAADJUST_TARGET}> --help-module stage)
    add_test (NAME adj-help-06 COMMAND $<TARGET_FILE:${DNAADJUST_TARGET}> --help-module output)
    add_test (NAME adj-help-07 COMMAND $<TARGET_FILE:${DNAADJUST_TARGET}> --help-module export)
    add_test (NAME adj-help-08 COMMAND $<TARGET_FILE:${DNAADJUST_TARGET}> --help-module generic)
    add_test (NAME adj-version-01 COMMAND $<TARGET_FILE:${DNAADJUST_TARGET}> --version)

    add_test (NAME plt-help-00 COMMAND $<TARGET_FILE:${DNAPLOT_TARGET}> -h)
    add_test (NAME plt-help-01 COMMAND $<TARGET_FILE:${DNAPLOT_TARGET}> --help-module standard)
    add_test (NAME plt-help-02 COMMAND $<TARGET_FILE:${DNAPLOT_TARGET}> --help-module data)
    add_test (NAME plt-help-03 COMMAND $<TARGET_FILE:${DNAPLOT_TARGET}> --help-module mapping)
    add_test (NAME plt-help-04 COMMAND $<TARGET_FILE:${DNAPLOT_TARGET}> --help-module title)
    add_test (NAME plt-help-05 COMMAND $<TARGET_FILE:${DNAPLOT_TARGET}> --help-module generic)
    add_test (NAME plt-version-01 COMMAND $<TARGET_FILE:${DNAPLOT_TARGET}> --version)

    add_test (NAME dyn-help-00 COMMAND $<TARGET_FILE:${DYNADJUST_TARGET}> -h)
    add_test (NAME dyn-help-01 COMMAND $<TARGET_FILE:${DYNADJUST_TARGET}> --help-module standard)
    add_test (NAME dyn-help-02 COMMAND $<TARGET_FILE:${DYNADJUST_TARGET}> --help-module generic)
    add_test (NAME dyn-version-01 COMMAND $<TARGET_FILE:${DYNADJUST_TARGET}> --version)

    # 13. Handle database id, output measurement to stations, output adjustment information on each iteration, scale normals to unity, report results, type-b
    add_test (NAME import-dbid COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> -n dbid dsg.stn dsg.msr dsg-ignored.msr) 
    add_test (NAME adjust-dbid-01 COMMAND $<TARGET_FILE:${DNAADJUST_TARGET}> dbid --output-adj-msr --output-database-ids --output-ignored-msrs --sort-adj-msr-field 6 --output-msr-to-stn --sort-msr-to-stn-field 1 --output-iter-adj-stn --output-iter-adj-stat --output-iter-adj-msr --output-iter-cmp-msr)
    add_test (NAME adjust-dbid-02 COMMAND $<TARGET_FILE:${DNAADJUST_TARGET}> dbid --scale-normals-to-unity)
    add_test (NAME adjust-dbid-03 COMMAND $<TARGET_FILE:${DNAADJUST_TARGET}> dbid --report-results --output-adj-msr --output-pos-uncertainty --output-apu-vcv-units 1 --constraints "236300210,CCC")
    add_test (NAME adjust-dbid-04 COMMAND $<TARGET_FILE:${DNAADJUST_TARGET}> dbid --comments "This is a comment that is quite lengthy in content and is relatively meaningless.  Feel free to delete this comment.")
    add_test (NAME adjust-dbid-05 COMMAND $<TARGET_FILE:${DNAADJUST_TARGET}> -p dbid.dnaproj)
    add_test (NAME adjust-dbid-06 COMMAND $<TARGET_FILE:${DNAADJUST_TARGET}> dbid --type-b-sd-glob "0.01,0.01,0.035" --type-b-sd-file dsg.typeb --output-pos)

    # 13. Renaming, no export of flagged stations, no network name, non-contiguous networks, import block, nearby stations
    add_test (NAME import-misc-01 COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> -n misc dsg.stn dsg.msr dsg-ignored.msr --stn-renaming-file dsg.renaming --export-xml)
    add_test (NAME import-misc-02 COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> -n misc dsg.stn dsg-retained.msr dsg-retained.msr --search-similar-msr --export-dna --export-xml --single-xml-file --output-msr-to-stn)
    add_test (NAME import-misc-03 COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> -n misc dsg.stn dsg.msr dsg-retained.msr dsg-retained.msr dsg-ignored.msr --search-similar-msr --ignore-similar-msr --flag-unused-stations)
    add_test (NAME import-misc-04 COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> -n misc dsg.stn dsg.msr dsg-retained.msr dsg-ignored.msr --export-xml)
    add_test (NAME segment-misc-plt COMMAND $<TARGET_FILE:${DNASEGMENT_TARGET}> misc --min 500 --max 500)
    add_test (NAME plot-misc-omit COMMAND $<TARGET_FILE:${DNAPLOT_TARGET}> misc --omit-measurements --omit-title-block)
    add_test (NAME plot-misc-all COMMAND $<TARGET_FILE:${DNAPLOT_TARGET}> misc --plot-msr-types "ABCDEGHIJKLMNPQRSVXYZ")
    add_test (NAME plot-misc-gph-stn COMMAND $<TARGET_FILE:${DNAPLOT_TARGET}> misc --supress-pdf --graph-stn)
    add_test (NAME plot-misc-gph-msr COMMAND $<TARGET_FILE:${DNAPLOT_TARGET}> misc --supress-pdf --graph-msr)
    add_test (NAME import-misc-05 COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> -n misc ./misc.xml sa-gnss-stn.xml --flag-unused-stations --export-dna --export-xml)
    add_test (NAME reftran-misc-01 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> misc -r itrf1991 --export-xml --verb 3)
    add_test (NAME reftran-misc-02 COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> -p misc.dnaproj)
    add_test (NAME segment-misc-01 COMMAND $<TARGET_FILE:${DNASEGMENT_TARGET}> misc --min 500 --max 500 --starting-stns 230900140)
    add_test (NAME segment-misc-02 COMMAND $<TARGET_FILE:${DNASEGMENT_TARGET}> misc --min 50 --max 100 --starting-stns 230900140 --seg-file ./misc.seg)
    add_test (NAME import-misc-06 COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> -n misc --import-block-stn-msr 2)
    add_test (NAME import-misc-07 COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> -p misc.dnaproj)
    add_test (NAME import-noname-01 COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> dsg.stn dsg.msr)
    add_test (NAME import-noname-02 COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> dsg.stn dsg.msr)
    add_test (NAME import-noncontiguous COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> -n noncontig dsg.stn dsg.msr skye-tutorial.stn skye-tutorial.msr -r itrf2000 --override-input-ref-frame)
    add_test (NAME segment-noncontiguous-01 COMMAND $<TARGET_FILE:${DNASEGMENT_TARGET}> noncontig --min 3 --max 5 --search-level 1 --test-integrity --verbose 3)
    add_test (NAME segment-noncontiguous-02 COMMAND $<TARGET_FILE:${DNASEGMENT_TARGET}> noncontig --min 3 --max 5 --contiguous-blocks 0 --search-level 1 --test-integrity  --verbose 3)
    add_test (NAME segment-noncontiguous-03 COMMAND $<TARGET_FILE:${DNASEGMENT_TARGET}> -p noncontig.dnaproj)
    add_test (NAME import-block-01 COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> -n misc ./miscstn.xml ./miscmsr.xml)
    add_test (NAME segment-block COMMAND $<TARGET_FILE:${DNASEGMENT_TARGET}> misc --min 2 --max 3)
    add_test (NAME import-block-02 COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> -n misc dsg.stn dsg.msr --seg-file misc.seg --import-block 2)
    add_test (NAME import-include-msrs COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> -n misc ./miscstn.xml ./miscmsr.xml --include-msr-types "ABCDEGHIJKLMNPQRSVXYZ")
    add_test (NAME import-exclude-msrs COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> -n null ./miscstn.xml ./miscmsr.xml --exclude-msr-types "ABCDEGHIJKLMNPQRSVXYZ" --flag --export-dna)
    add_test (NAME import-stn-msr-connect COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> -n connect ./miscstn.xml ./miscmsr.xml --include-stns-assoc-msrs "365300060,360100260,269100210")
    add_test (NAME import-all COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> -n all 
        ./miscstn.xml 
        ./miscmsr.xml
        apr.ITRF2008.04.06.2020.stn 
        apr.ITRF2008.04.06.2020.msr
        gnss-network.stn 
        gnss-network.msr
        sa-gnss-stn.xml 
        sa-gnss-msr.xml
        skye-tutorial.stn 
        skye-tutorial.msr 
        TEST_ITRF05.SNX 
        urban-network.stn 
        urban-network.msr)
    add_test (NAME plot-all-01 COMMAND $<TARGET_FILE:${DNAPLOT_TARGET}> all --map-projection 1)
    add_test (NAME plot-all-02 COMMAND $<TARGET_FILE:${DNAPLOT_TARGET}> all --map-projection 2)
    add_test (NAME plot-all-03 COMMAND $<TARGET_FILE:${DNAPLOT_TARGET}> all --map-projection 3)
    add_test (NAME plot-all-04 COMMAND $<TARGET_FILE:${DNAPLOT_TARGET}> all --map-projection 4)
    add_test (NAME plot-all-05 COMMAND $<TARGET_FILE:${DNAPLOT_TARGET}> all --map-projection 5)
    add_test (NAME plot-all-06 COMMAND $<TARGET_FILE:${DNAPLOT_TARGET}> all --map-projection 6)
    add_test (NAME plot-all-07 COMMAND $<TARGET_FILE:${DNAPLOT_TARGET}> all --map-projection 7)
    add_test (NAME plot-all-08 COMMAND $<TARGET_FILE:${DNAPLOT_TARGET}> all --bounding "-36.3000,-200.0,-36.4500,200.0")
    add_test (NAME plot-all-09 COMMAND $<TARGET_FILE:${DNAPLOT_TARGET}> all --bounding "100.3000,130.0,-100.4500,150.0")
    add_test (NAME plot-all-10 COMMAND $<TARGET_FILE:${DNAPLOT_TARGET}> all --bounding "80.3000,130.0,-80.4500,150.0")
    add_test (NAME plot-all-11 COMMAND $<TARGET_FILE:${DNAPLOT_TARGET}> all --bounding "60.3000,100.0,-60.4500,180.0")
    add_test (NAME plot-all-12 COMMAND $<TARGET_FILE:${DNAPLOT_TARGET}> all --centre-latitude "-36.0101" --centre-longitude "145.0202")
    add_test (NAME plot-all-13 COMMAND $<TARGET_FILE:${DNAPLOT_TARGET}> all --centre-latitude "-66.0101" --centre-longitude "145.0202")
    add_test (NAME plot-all-14 COMMAND $<TARGET_FILE:${DNAPLOT_TARGET}> all --centre-latitude "-86.0101" --centre-longitude "145.0202")
    add_test (NAME plot-all-15 COMMAND $<TARGET_FILE:${DNAPLOT_TARGET}> all --title "A plot generated by cmake testing." --org-unit "Geoscience Australia" --org-subunit "National Geodesy")
    add_test (NAME plot-all-16 COMMAND $<TARGET_FILE:${DNAPLOT_TARGET}> all --export-png --omit-title)
    add_test (NAME plot-all-17 COMMAND $<TARGET_FILE:${DNAPLOT_TARGET}> all --export-png --plate MORVEL56_plates.dig --map 8)
    add_test (NAME import-misc-08 COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> -n contig 
        gnss-network.stn 
        gnss-network.msr
        sa-gnss-stn.xml 
        sa-gnss-msr.xml
        skye-tutorial.stn 
        skye-tutorial.msr)
    add_test (NAME segment-misc-08 COMMAND $<TARGET_FILE:${DNASEGMENT_TARGET}> contig --contiguous-blocks 0)
    add_test (NAME import-misc-08a COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> -n contig --import-block-stn-msr 1 --export-dna)
    add_test (NAME import-misc-08b COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> -n contig --import-contiguous-stn-msr 0 --export-dna)

    # 14. Exceptions
    ### IMPORT ###
    # no options
    add_test (NAME import-no-option COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}>)
    # invalid options
    add_test (NAME import-invalid-option COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> -x)
    # invalid input file
    add_test (NAME import-invalid-file COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> -n invalid invalidfile.txt)
    # the --search-nearby-stn option in import will return a failure if nearby stations are found.  In this case, check it is non-zero (i.e. it will fail!)
    add_test (NAME import-nearby COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> -n urban urban-network.stn urban-network.msr --search-nearby-stn --nearby-stn-buffer 10)
    # files don't exist
    add_test (NAME import-no-files COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> -n bigdata bigdata.stn bigdata.msr)
    # reference frame doesn't exist
    add_test (NAME import-no-frame COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> -n sample -r itrf1975 dsg.stn dsg.msr)
    # no stations and measurements (generated above)
    add_test (NAME import-no-data COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> -n null ./null.stn ./null.msr)
    # no DynaML.xsd - remove file in platform-independent manner (using CMake command) by setting up a test dependency
    add_test (NAME import-no-xsd-file-remove COMMAND ${CMAKE_COMMAND} -E remove DynaML.xsd)
    add_test (NAME import-no-xsd-file COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> -n gnss gnss-networkstn.xml gnss-networkmsr.xml)
    set_tests_properties(import-no-xsd-file PROPERTIES DEPENDS import-no-xsd-file-remove TIMEOUT 5)
    # Restore DynaML.xsd after testing its absence
    add_test (NAME import-no-xsd-file-restore COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_SOURCE_DIR}/../sampleData/DynaML.xsd DynaML.xsd)
    set_tests_properties(import-no-xsd-file-restore PROPERTIES DEPENDS import-no-xsd-file)
    # no geoid file
    add_test (NAME import-no-geoid COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> -n geoid --geo-file ./nofile.geo)
    set_tests_properties(import-no-geoid PROPERTIES DEPENDS import-no-xsd-file-restore)
    # no project file
    add_test (NAME import-no-project COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> -p ./nofile.dnaproj)
    # measurement stns not in station file
    add_test (NAME import-no-stns COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> -n nostn urban-network.stn dsg.msr)
    # no help module
    add_test (NAME import-no-help COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> --help-module cplusplus)
    # excessive block and network numbers
    add_test (NAME import-misc-block-size COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> -n contig --import-block-stn-msr 5000 --export-dna)
    add_test (NAME import-misc-network-size COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> -n contig --import-contiguous-stn-msr 5000 --export-dna)

    set_tests_properties(import-no-option import-invalid-option import-invalid-file import-nearby import-no-files import-no-frame
        import-no-data import-no-xsd-file import-no-geoid import-no-project import-no-stns import-no-help
        import-misc-block-size import-misc-network-size import-misc-07
        PROPERTIES WILL_FAIL TRUE)

    ### REFTRAN ###
    # no options
    add_test (NAME reftran-no-option COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}>)
    # invalid options
    add_test (NAME reftran-invalid-option COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> -x)
    # no project file
    add_test (NAME reftran-no-project COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> -p ./nofile.dnaproj)
    # no name
    add_test (NAME reftran-no-name COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> -r "gda94")
    # no binary files
    add_test (NAME reftran-no-bst-bms COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> extremelylargenetwork -r "itrf1990")
    # reference frame doesn't exist
    add_test (NAME reftran-bad-frame COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> all -r itrf1975)
    # epoch doesn't exist
    add_test (NAME reftran-bad-epoch COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> all -r "itrf2000" -e "yyyy")
    # no plate boundary file
    add_test (NAME reftran-no-discont-file COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> all -r "ITRF2008" --plate-boundary-file nofile.txt)
    # no euler pole file
    add_test (NAME reftran-no-euler-file COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> all -r "ITRF2008" --plate-pole-file nofile.txt)
    # no plate model option argument
    add_test (NAME reftran-no-model-option COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> all -r "ITRF2014" -e "31.12.9999" --plate-model-option)
    # no help module
    add_test (NAME reftran-no-help COMMAND $<TARGET_FILE:${DNAREFTRAN_TARGET}> --help-module cplusplus)


    set_tests_properties(
        reftran-no-option reftran-invalid-option reftran-no-project reftran-no-name reftran-no-bst-bms 
        reftran-bad-frame reftran-bad-epoch reftran-no-discont-file reftran-no-euler-file
        reftran-no-model-option reftran-no-help
        PROPERTIES WILL_FAIL TRUE)

    ### GEOID ###
    # no options
    add_test (NAME geoid-no-option COMMAND $<TARGET_FILE:${DNAGEOID_TARGET}>)
    # invalid options
    add_test (NAME geoid-invalid-option COMMAND $<TARGET_FILE:${DNAGEOID_TARGET}> -x)
    # no project file
    add_test (NAME geoid-no-project COMMAND $<TARGET_FILE:${DNAGEOID_TARGET}> -p ./nofile.dnaproj)
    # no standard options
    add_test (NAME geoid-no-std-options COMMAND $<TARGET_FILE:${DNAGEOID_TARGET}> -g gnss-network-geoid.gsb)
    # no interpolation input file
    add_test (NAME geoid-no-input-file COMMAND $<TARGET_FILE:${DNAGEOID_TARGET}> largenet)
    # no standard options
    add_test (NAME geoid-no-bst COMMAND $<TARGET_FILE:${DNAGEOID_TARGET}> -g gnss-network-geoid.gsb largenet)
    # no help module
    add_test (NAME geoid-no-help COMMAND $<TARGET_FILE:${DNAGEOID_TARGET}> --help-module cplusplus)
    # no geoid file
    add_test (NAME geoid-no-gsb COMMAND $<TARGET_FILE:${DNAGEOID_TARGET}> gnss -g nogridfile.gsb)

    set_tests_properties(
        geoid-no-option geoid-invalid-option geoid-no-project geoid-no-std-options geoid-no-input-file geoid-no-bst geoid-no-help geoid-no-gsb
        PROPERTIES WILL_FAIL TRUE)

    ### SEGMENT ###
    # no options
    add_test (NAME segment-no-option COMMAND $<TARGET_FILE:${DNASEGMENT_TARGET}>)
    # invalid options
    add_test (NAME segment-invalid-option COMMAND $<TARGET_FILE:${DNASEGMENT_TARGET}> -x)
    # no project file
    add_test (NAME segment-no-project COMMAND $<TARGET_FILE:${DNASEGMENT_TARGET}> -p ./nofile.dnaproj)
    # no binary files
    add_test (NAME segment-no-bst-bms COMMAND $<TARGET_FILE:${DNASEGMENT_TARGET}> extremelylargenetwork)
    # no net file
    add_test (NAME segment-no-net-01 COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> -n nonetfile dsg.stn dsg.msr)
    add_test (NAME segment-no-net-02 COMMAND $<TARGET_FILE:${DNASEGMENT_TARGET}> nonetfile --net-file)
    # no option arguments
    add_test (NAME segment-no-opt-arg-01 COMMAND $<TARGET_FILE:${DNASEGMENT_TARGET}> nonetfile --starting-stns)
    add_test (NAME segment-no-opt-arg-02 COMMAND $<TARGET_FILE:${DNASEGMENT_TARGET}> nonetfile --min)
    add_test (NAME segment-no-opt-arg-03 COMMAND $<TARGET_FILE:${DNASEGMENT_TARGET}> nonetfile --max)
    add_test (NAME segment-no-opt-arg-04 COMMAND $<TARGET_FILE:${DNASEGMENT_TARGET}> nonetfile --contiguous)
    add_test (NAME segment-no-opt-arg-05 COMMAND $<TARGET_FILE:${DNASEGMENT_TARGET}> nonetfile --search)
    # all measurements ignored
    add_test (NAME segment-all-ignored-01 COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> -n allignored stn_allignored.xml msr_allignored.xml)
    add_test (NAME segment-all-ignored-02 COMMAND $<TARGET_FILE:${DNASEGMENT_TARGET}> allignored)
    # no help module
    add_test (NAME segment-no-help COMMAND $<TARGET_FILE:${DNASEGMENT_TARGET}> --help-module cplusplus)

    set_tests_properties(
        segment-no-option segment-invalid-option segment-no-project segment-no-bst-bms segment-no-net-02 segment-no-help
        segment-no-opt-arg-01 segment-no-opt-arg-02 segment-no-opt-arg-03 segment-no-opt-arg-04 segment-no-opt-arg-05
        segment-all-ignored-02
        PROPERTIES WILL_FAIL TRUE)

    ### ADJUST ###
    # no options
    add_test (NAME adjust-no-option COMMAND $<TARGET_FILE:${DNAADJUST_TARGET}>)
    # invalid options
    add_test (NAME adjust-invalid-option COMMAND $<TARGET_FILE:${DNAADJUST_TARGET}> -x)
    # no project file
    add_test (NAME adjust-no-project COMMAND $<TARGET_FILE:${DNAADJUST_TARGET}> -p ./nofile.dnaproj)
    # no network name
    add_test (NAME adjust-no-name COMMAND $<TARGET_FILE:${DNAADJUST_TARGET}> --phased)
    # no help module
    add_test (NAME adjust-no-help COMMAND $<TARGET_FILE:${DNAADJUST_TARGET}> --help-module cplusplus)
    # no binary files
    add_test (NAME adjust-no-bst-bms COMMAND $<TARGET_FILE:${DNAADJUST_TARGET}> largenet)
    # no segment file
    add_test (NAME import-no-seg COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> -n noseg dsg.stn dsg.msr)
    add_test (NAME adjust-no-seg COMMAND $<TARGET_FILE:${DNAADJUST_TARGET}> noseg --phased)
    # no memory mapped files from previous adjustment
    add_test (NAME import-no-memmap COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> -n memmap dsg.stn dsg.msr)
    add_test (NAME segment-no-memmap COMMAND $<TARGET_FILE:${DNASEGMENT_TARGET}> memmap --min 2 --max 3)
    add_test (NAME adjust-no-memmap-01 COMMAND $<TARGET_FILE:${DNAADJUST_TARGET}> memmap --staged --create --purge)
    add_test (NAME adjust-no-memmap-02 COMMAND $<TARGET_FILE:${DNAADJUST_TARGET}> memmap --staged)
    add_test (NAME adjust-no-stn-exists COMMAND $<TARGET_FILE:${DNAADJUST_TARGET}> memmap --constraints "no-name,CCC")
    add_test (NAME adjust-no-constraint COMMAND $<TARGET_FILE:${DNAADJUST_TARGET}> memmap --constraints "236300210,AAA")
    # all measurements ignored
    add_test (NAME adjust-all-ignored-01 COMMAND $<TARGET_FILE:${DNAADJUST_TARGET}> allignored)

    set_tests_properties(
        adjust-no-option adjust-invalid-option adjust-no-project adjust-no-name adjust-no-help adjust-no-bst-bms adjust-no-seg 
        adjust-no-memmap-02 adjust-no-stn-exists adjust-no-constraint adjust-all-ignored-01
        PROPERTIES WILL_FAIL TRUE)

    ### PLOT ###
    # no options
    add_test (NAME plot-no-option COMMAND $<TARGET_FILE:${DNAPLOT_TARGET}>)
    # invalid options
    add_test (NAME plot-invalid-option COMMAND $<TARGET_FILE:${DNAPLOT_TARGET}> -x)
    # no project file
    add_test (NAME plot-no-project COMMAND $<TARGET_FILE:${DNAPLOT_TARGET}> -p ./nofile.dnaproj)
    # no standard options
    add_test (NAME plot-no-std-options COMMAND $<TARGET_FILE:${DNAPLOT_TARGET}> -g gnss-network-plot.gsb)
    # no such network
    add_test (NAME plot-no-network COMMAND $<TARGET_FILE:${DNAPLOT_TARGET}> largenet)
    # no network-name
    add_test (NAME plot-no-name COMMAND $<TARGET_FILE:${DNAPLOT_TARGET}> --centre-station "PM1")
    # no help module
    add_test (NAME plot-no-help COMMAND $<TARGET_FILE:${DNAPLOT_TARGET}> --help-module cplusplus)
    # bounding box error
    add_test (NAME plot-box-ew COMMAND $<TARGET_FILE:${DNAPLOT_TARGET}> all --bounding "-36.3000,146.4500,-36.4500,145.4500")
    add_test (NAME plot-box-ns COMMAND $<TARGET_FILE:${DNAPLOT_TARGET}> all --bounding "-36.4500,145.4500,-36.3000,146.4500")
    # no seg file, no apu file, no cor file, no plate
    add_test (NAME import-no-seg-apu-cor COMMAND $<TARGET_FILE:${DNAIMPORT_TARGET}> -n nosegapucor dsg.stn dsg.msr)
    add_test (NAME plot-no-seg COMMAND $<TARGET_FILE:${DNAPLOT_TARGET}> nosegapucor --phased)
    add_test (NAME plot-no-cor COMMAND $<TARGET_FILE:${DNAPLOT_TARGET}> nosegapucor --correction-arrows)
    add_test (NAME plot-no-apu COMMAND $<TARGET_FILE:${DNAPLOT_TARGET}> nosegapucor --error-ellipses)
    add_test (NAME plot-no-plate COMMAND $<TARGET_FILE:${DNAPLOT_TARGET}> nosegapucor --plate ./noplate.file)

    set_tests_properties(
        plot-no-option plot-invalid-option plot-no-project plot-no-std-options plot-no-network plot-no-name plot-no-help
        plot-box-ew plot-box-ns plot-no-seg plot-no-cor plot-no-apu plot-no-plate
        PROPERTIES WILL_FAIL TRUE)

    ### DYNADJUST ###
    # no options
    add_test (NAME dyna-no-option COMMAND $<TARGET_FILE:${DYNADJUST_TARGET}>)
    # invalid option
    add_test (NAME dyna-invalid-option COMMAND $<TARGET_FILE:${DYNADJUST_TARGET}> -x)
    # no help module
    add_test (NAME dyna-no-help COMMAND $<TARGET_FILE:${DYNADJUST_TARGET}> --help-module cplusplus)
    # no project file and network name
    add_test (NAME dyna-no-proj-no-name COMMAND $<TARGET_FILE:${DYNADJUST_TARGET}> --import)
    # no programs
    add_test (NAME dyna-no-programs COMMAND $<TARGET_FILE:${DYNADJUST_TARGET}>  -p sample.dnaproj)
    # no project file
    add_test (NAME dyna-no-project COMMAND $<TARGET_FILE:${DYNADJUST_TARGET}> -p ./nofile.dnaproj)
    # name doesn't match project file
    add_test (NAME dyna-name-proj-mismatch COMMAND $<TARGET_FILE:${DYNADJUST_TARGET}> -p ./nofile.dnaproj -n noname)

    set_tests_properties(
        dyna-no-option dyna-invalid-option dyna-no-help dyna-no-proj-no-name dyna-no-programs
        dyna-no-project dyna-name-proj-mismatch
        PROPERTIES WILL_FAIL TRUE)


    # set execution dependencies (the execution of tests must be sequential)
    set_tests_properties(
        import-gnss-network geoid-gnss-network adjust-gnss-network
        PROPERTIES RUN_SERIAL TRUE)
    set_tests_properties(
        import-urban-network geoid-urban-network segment-urban-network adjust-urban-network
        PROPERTIES RUN_SERIAL TRUE)
    set_tests_properties(
        import-urban-network-thread reftran-urban-network-thread geoid-urban-network-thread segment-urban-network-thread adjust-urban-network-thread-01
        PROPERTIES RUN_SERIAL TRUE)

    set_tests_properties(test-gnss-network PROPERTIES DEPENDS adjust-gnss-network)

    set_tests_properties(ref-itrf-pmm-06 PROPERTIES DEPENDS ref-itrf-pmm-05)
    #set_tests_properties(ref-itrf-pmm-07 PROPERTIES DEPENDS ref-itrf-pmm-06)

    # Force all tests to run serially (prevent parallel execution)
    get_property(all_tests DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY TESTS)
    if(all_tests)
        set_tests_properties(${all_tests} PROPERTIES RUN_SERIAL TRUE)
        message(STATUS "All tests set to run serially (RUN_SERIAL=TRUE)")
    endif()

    # Setup testing directory

    set(TESTING_DIR ${CMAKE_SOURCE_DIR}/../testing)

    if(EXISTS ${TESTING_DIR})
        file(REMOVE_RECURSE ${TESTING_DIR})
    endif()

    file(MAKE_DIRECTORY ${TESTING_DIR})
    file(COPY ${CMAKE_SOURCE_DIR}/../sampleData/ DESTINATION ${TESTING_DIR})

    get_property(all_tests DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY TESTS)
    foreach(test ${all_tests})
        set_tests_properties(${test} PROPERTIES WORKING_DIRECTORY ${TESTING_DIR})
    endforeach()

    message (STATUS "Finished setting up tests")
    message ("")
    message ("")
    message (" **")
    message (" * To run tests, type `ctest`")
    message ("   To run a specific test, type `ctest -R <testname>`")
    message ("   To run all tests in parallel, type `ctest -j <number>`")
    message ("   Note: Some tests depend on others and might fail in parallel.  *")
    message ("                                                                 **")
    message ("")
    message("---------------------------------------------------------------------")
endif ()

message ("")
message ("")

