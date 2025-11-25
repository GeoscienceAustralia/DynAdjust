# <dnaimport/...> build rules
project (dnaimport)

add_definitions(-DBUILD_DNAIMPORT_DLL -DBUILD_EXPORTDNATYPES -DMKL_ILP64)

include_directories (${PROJECT_SOURCE_DIR})

add_library (${PROJECT_NAME} ${DNA_LIBRARY_TYPE}
             ${CMAKE_SOURCE_DIR}/include/io/dnaioaml.cpp
             ${CMAKE_SOURCE_DIR}/include/io/dnaioasl.cpp
             ${CMAKE_SOURCE_DIR}/include/io/dnaiobase.cpp
             ${CMAKE_SOURCE_DIR}/include/io/dnaiobms.cpp
             ${CMAKE_SOURCE_DIR}/include/io/dnaiobst.cpp
             ${CMAKE_SOURCE_DIR}/include/io/dnaiodna.cpp
             ${CMAKE_SOURCE_DIR}/include/io/dnaiomap.cpp
             ${CMAKE_SOURCE_DIR}/include/io/dnaioscalar.cpp
             ${CMAKE_SOURCE_DIR}/include/io/dnaioseg.cpp
             ${CMAKE_SOURCE_DIR}/include/io/dnaiosnxread.cpp
             ${CMAKE_SOURCE_DIR}/include/math/dnamatrix_contiguous.cpp
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

target_link_libraries (${PROJECT_NAME} ${DNA_LIBRARIES})

install (TARGETS ${PROJECT_NAME} DESTINATION lib)
