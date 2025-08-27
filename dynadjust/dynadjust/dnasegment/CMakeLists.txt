# <dnasegment/...> build rules
project (dnasegment)

add_definitions(-DBUILD_DNASEGMENT_DLL -DBUILD_EXPORTDNATYPES)

include_directories (${PROJECT_SOURCE_DIR})

add_library (${PROJECT_NAME} SHARED
             ${CMAKE_SOURCE_DIR}/include/io/aml_file.cpp
             ${CMAKE_SOURCE_DIR}/include/io/asl_file.cpp
             ${CMAKE_SOURCE_DIR}/include/io/dynadjust_file.cpp
             ${CMAKE_SOURCE_DIR}/include/io/bms_file.cpp
             ${CMAKE_SOURCE_DIR}/include/io/bst_file.cpp
             ${CMAKE_SOURCE_DIR}/include/io/map_file.cpp
             ${CMAKE_SOURCE_DIR}/include/io/seg_file.cpp
             ${CMAKE_SOURCE_DIR}/include/measurement_types/dnastation.cpp
             ${CMAKE_SOURCE_DIR}/include/measurement_types/dnamsrtally.cpp
             ${CMAKE_SOURCE_DIR}/include/parameters/dnaellipsoid.cpp
             dnasegment.cpp
             ${CMAKE_SOURCE_DIR}/dynadjust.rc)

target_precompile_headers (${PROJECT_NAME} PUBLIC precompile.h)

target_link_libraries (${PROJECT_NAME} ${DNA_LIBRARIES})

