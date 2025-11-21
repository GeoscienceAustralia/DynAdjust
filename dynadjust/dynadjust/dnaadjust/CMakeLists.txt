# <dnaadjust/...> build rules
project (dnaadjust)

add_definitions(-DBUILD_DNAADJUST_DLL -DBUILD_EXPORTDNATYPES -DMULTI_THREAD_ADJUST -DMKL_ILP64)

include_directories (${PROJECT_SOURCE_DIR})

add_library (${PROJECT_NAME} SHARED
             ${CMAKE_SOURCE_DIR}/include/io/dnaioadj.cpp
             ${CMAKE_SOURCE_DIR}/include/io/dnaioaml.cpp
             ${CMAKE_SOURCE_DIR}/include/io/dnaioasl.cpp
             ${CMAKE_SOURCE_DIR}/include/io/dnaiobase.cpp
             ${CMAKE_SOURCE_DIR}/include/io/dnaiobms.cpp
             ${CMAKE_SOURCE_DIR}/include/io/dnaiobst.cpp
             ${CMAKE_SOURCE_DIR}/include/io/dnaiomap.cpp
             ${CMAKE_SOURCE_DIR}/include/io/dnaioseg.cpp
             ${CMAKE_SOURCE_DIR}/include/io/dnaiotbu.cpp
             ${CMAKE_SOURCE_DIR}/include/io/dnaiosnxwrite.cpp
             ${CMAKE_SOURCE_DIR}/include/parameters/dnadatum.cpp
             ${CMAKE_SOURCE_DIR}/include/parameters/dnaellipsoid.cpp
             ${CMAKE_SOURCE_DIR}/include/parameters/dnaprojection.cpp
             ${CMAKE_SOURCE_DIR}/include/functions/dnastringfuncs.cpp
             ${CMAKE_SOURCE_DIR}/include/math/dnamatrix_contiguous.cpp
             ${CMAKE_SOURCE_DIR}/include/measurement_types/dnagpspoint.cpp
             ${CMAKE_SOURCE_DIR}/include/measurement_types/dnameasurement.cpp
             ${CMAKE_SOURCE_DIR}/include/measurement_types/dnastation.cpp
             ${CMAKE_SOURCE_DIR}/include/measurement_types/dnamsrtally.cpp
             ${CMAKE_SOURCE_DIR}/include/memory/dnafile_mapping.cpp
             dnaadjust-stage.cpp
             dnaadjust.cpp
             ${CMAKE_SOURCE_DIR}/dynadjust.rc)

target_link_libraries (${PROJECT_NAME} ${DNA_LIBRARIES})

install (TARGETS ${PROJECT_NAME} DESTINATION lib)
