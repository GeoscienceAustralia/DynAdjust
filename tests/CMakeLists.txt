cmake_minimum_required(VERSION 3.16)
project(DynAdjustTests LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set output directory for binaries
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

# Include directories
include_directories(../dynadjust)
include_directories(../dynadjust/include)

# Find required packages
find_package(Boost REQUIRED COMPONENTS filesystem system timer date_time)

# Platform-specific libraries
if(APPLE)
    find_library(ACCELERATE_FRAMEWORK Accelerate)
    set(PLATFORM_LIBS ${ACCELERATE_FRAMEWORK})
elseif(WIN32)
    # Windows-specific libraries can be added here
    set(PLATFORM_LIBS "")
else()
    # Linux-specific libraries (e.g., BLAS/LAPACK) can be added here
    set(PLATFORM_LIBS "")
endif()

# C++ filesystem library for C++17
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "9.0")
        list(APPEND PLATFORM_LIBS stdc++fs)
    endif()
elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    # Only needed on Linux with Clang (not on macOS/Apple Clang)
    if(UNIX AND NOT APPLE)
        list(APPEND PLATFORM_LIBS stdc++fs)
    endif()
endif()

# Common source files for I/O tests
set(IO_COMMON_SOURCES
    ../dynadjust/include/io/dynadjust_file.cpp
    ../dynadjust/include/measurement_types/dnastation.cpp
)

# Test 1: Matrix test
add_executable(test_matrix
    test_matrix.cpp
    ../dynadjust/include/math/dnamatrix_contiguous.cpp
    ../dynadjust/include/ide/trace.cpp
)

target_link_libraries(test_matrix
    ${PLATFORM_LIBS}
)

target_compile_definitions(test_matrix PRIVATE
    __BINARY_NAME__="test_matrix"
    __BINARY_DESC__="Unit tests for matrix operations"
)

# Test 2: BST File Loader test
add_executable(test_bst_file
    test_bst_file.cpp
    ../dynadjust/include/io/bst_file.cpp
    ${IO_COMMON_SOURCES}
)

target_link_libraries(test_bst_file
    ${PLATFORM_LIBS}
    Boost::filesystem
    Boost::system
)

target_compile_definitions(test_bst_file PRIVATE
    __BINARY_NAME__="test_bst_file"
    __BINARY_DESC__="Unit tests for BstFile class"
)

# Test 3: ASL File test
add_executable(test_asl_file
    test_asl_file.cpp
    ../dynadjust/include/io/asl_file.cpp
    ${IO_COMMON_SOURCES}
)

target_link_libraries(test_asl_file
    ${PLATFORM_LIBS}
    Boost::filesystem
    Boost::system
)

target_compile_definitions(test_asl_file PRIVATE
    __BINARY_NAME__="test_asl_file"
    __BINARY_DESC__="Unit tests for AslFile class"
)

# Test 4: AML File Loader test
add_executable(test_aml_file_loader
    test_aml_file_loader.cpp
    ../dynadjust/include/io/aml_file.cpp
    ../dynadjust/include/io/bms_file.cpp
    ../dynadjust/include/measurement_types/dnamsrtally.cpp
    ${IO_COMMON_SOURCES}
)

target_link_libraries(test_aml_file_loader
    ${PLATFORM_LIBS}
    Boost::filesystem
    Boost::system
)

target_compile_definitions(test_aml_file_loader PRIVATE
    __BINARY_NAME__="test_aml_file_loader"
    __BINARY_DESC__="Unit tests for AML file loader"
)

# Test 5: BMS File Loader test
add_executable(test_bms_file
    test_bms_file.cpp
    ../dynadjust/include/io/bms_file.cpp
    ${IO_COMMON_SOURCES}
)

target_link_libraries(test_bms_file
    ${PLATFORM_LIBS}
    Boost::filesystem
    Boost::system
)

target_compile_definitions(test_bms_file PRIVATE
    __BINARY_NAME__="test_bms_file"
    __BINARY_DESC__="Unit tests for BmsFile class"
)

# Test 6: Network Data Loader test
add_executable(test_network_data_loader
    test_network_data_loader.cpp
    ../dynadjust/dynadjust/dnaadjust/network_data_loader.cpp
    ../dynadjust/dynadjust/dnaadjust/measurement_processor.cpp
    ../dynadjust/include/io/asl_file.cpp
    ../dynadjust/include/io/bst_file.cpp
    ../dynadjust/include/io/bms_file.cpp
    ../dynadjust/include/io/map_file.cpp
    ${IO_COMMON_SOURCES}
)

target_link_libraries(test_network_data_loader
    ${PLATFORM_LIBS}
    Boost::filesystem
    Boost::system
)

target_compile_definitions(test_network_data_loader PRIVATE
    __BINARY_NAME__="test_network_data_loader"
    __BINARY_DESC__="Unit tests for NetworkDataLoader class"
)


# Test 7: Measurement Processor test
add_executable(test_measurement_processor
    test_measurement_processor.cpp
    ../dynadjust/dynadjust/dnaadjust/measurement_processor.cpp
    ${IO_COMMON_SOURCES}
)

target_link_libraries(test_measurement_processor
    ${PLATFORM_LIBS}
    Boost::filesystem
    Boost::system
)

target_compile_definitions(test_measurement_processor PRIVATE
    __BINARY_NAME__="test_measurement_processor"
    __BINARY_DESC__="Unit tests for MeasurementProcessor class"
)

# Test 8: DynAdjust Printer test
add_executable(test_dnaadjust_printer
    test_dnaadjust_printer_simple.cpp
    ../dynadjust/include/measurement_types/dnameasurement.cpp
    ${IO_COMMON_SOURCES}
)

target_link_libraries(test_dnaadjust_printer
    ${PLATFORM_LIBS}
    Boost::filesystem
    Boost::system
    Boost::timer
    Boost::chrono
)

target_compile_definitions(test_dnaadjust_printer PRIVATE
    __BINARY_NAME__="test_dnaadjust_printer"
    __BINARY_DESC__="Unit tests for DynAdjustPrinter class"
)

# Enable testing
enable_testing()

# Add tests
add_test(NAME MatrixTest COMMAND test_matrix)
add_test(NAME BstFileTest COMMAND test_bst_file)
add_test(NAME AslFileTest COMMAND test_asl_file)
add_test(NAME AmlFileLoaderTest COMMAND test_aml_file_loader)
add_test(NAME BmsFileTest COMMAND test_bms_file)
add_test(NAME NetworkDataLoaderTest COMMAND test_network_data_loader)
add_test(NAME MeasurementProcessorTest COMMAND test_measurement_processor)
add_test(NAME DynAdjustPrinterTest COMMAND test_dnaadjust_printer)

# Custom target to run all tests
add_custom_target(run_tests
    COMMAND ${CMAKE_CTEST_COMMAND} --verbose
    DEPENDS test_matrix test_bst_file_loader test_asl_file_loader test_aml_file_loader test_bms_file_loader test_network_data_loader test_measurement_processor test_dnaadjust_printer
    COMMENT "Running all tests"
)

# Custom target equivalent to 'make all'
add_custom_target(tests_all
    DEPENDS test_matrix test_bst_file_loader test_asl_file_loader test_aml_file_loader test_bms_file_loader test_network_data_loader test_measurement_processor test_dnaadjust_printer
    COMMENT "Building all tests"
)
