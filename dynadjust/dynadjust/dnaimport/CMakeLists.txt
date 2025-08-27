# <dnaimport/...> build rules
project (dnaimport)

add_definitions(-DBUILD_DNAIMPORT_DLL -DBUILD_EXPORTDNATYPES)

include_directories (${PROJECT_SOURCE_DIR})

add_library (${PROJECT_NAME} SHARED
             ${CMAKE_SOURCE_DIR}/include/io/aml_file.cpp
             ${CMAKE_SOURCE_DIR}/include/io/asl_file.cpp
             ${CMAKE_SOURCE_DIR}/include/io/dynadjust_file.cpp
             ${CMAKE_SOURCE_DIR}/include/io/bms_file.cpp
             ${CMAKE_SOURCE_DIR}/include/io/bst_file.cpp
             ${CMAKE_SOURCE_DIR}/include/io/dnaiodna.cpp
             ${CMAKE_SOURCE_DIR}/include/io/map_file.cpp
             ${CMAKE_SOURCE_DIR}/include/io/dnaioscalar.cpp
             ${CMAKE_SOURCE_DIR}/include/io/seg_file.cpp
             ${CMAKE_SOURCE_DIR}/include/io/dnaiosnxread.cpp
             ${CMAKE_SOURCE_DIR}/include/math/dnamatrix_contiguous.cpp
             ${CMAKE_SOURCE_DIR}/include/ide/trace.cpp
             ${CMAKE_SOURCE_DIR}/include/functions/dnastringfuncs.cpp
             ${CMAKE_SOURCE_DIR}/include/parameters/dnadatum.cpp
             ${CMAKE_SOURCE_DIR}/include/parameters/dnaellipsoid.cpp
             ${CMAKE_SOURCE_DIR}/include/parameters/dnaprojection.cpp
             ${CMAKE_SOURCE_DIR}/include/measurement_types/dnaangle.cpp
             ${CMAKE_SOURCE_DIR}/include/measurement_types/dnacoordinate.cpp
             ${CMAKE_SOURCE_DIR}/include/measurement_types/dnadirection.cpp
             ${CMAKE_SOURCE_DIR}/include/measurement_types/dnadirectionset.cpp
             ${CMAKE_SOURCE_DIR}/include/measurement_types/dnadistance.cpp
             ${CMAKE_SOURCE_DIR}/include/measurement_types/dnagpsbaseline.cpp
             ${CMAKE_SOURCE_DIR}/include/measurement_types/dnagpspoint.cpp
             ${CMAKE_SOURCE_DIR}/include/measurement_types/dnaheight.cpp
             ${CMAKE_SOURCE_DIR}/include/measurement_types/dnaheightdifference.cpp
             ${CMAKE_SOURCE_DIR}/include/measurement_types/dnameasurement.cpp
             ${CMAKE_SOURCE_DIR}/include/measurement_types/dnastation.cpp
             ${CMAKE_SOURCE_DIR}/include/measurement_types/dnamsrtally.cpp
             ${CMAKE_SOURCE_DIR}/include/measurement_types/dnastntally.cpp
             dnaparser_pskel.cxx
             dnaparser_pimpl.cxx
             dnainterop.cpp
             ${CMAKE_SOURCE_DIR}/dynadjust.rc)

target_precompile_headers (${PROJECT_NAME} PUBLIC precompile.h)

target_link_libraries (${PROJECT_NAME} ${DNA_LIBRARIES})
