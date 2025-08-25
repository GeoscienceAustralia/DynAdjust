# <dnaplot/...> build rules
project (dnaplot)

#add_definitions(-DBUILD_DNAPLOT_DLL -DBUILD_EXPORTDNATYPES -DMKL_ILP64 -fopenmp)
add_definitions(-DBUILD_DNAPLOT_DLL -DBUILD_EXPORTDNATYPES)

include_directories (${PROJECT_SOURCE_DIR})

add_library (${PROJECT_NAME} SHARED
             ${CMAKE_SOURCE_DIR}/include/io/bms_file.cpp
             ${CMAKE_SOURCE_DIR}/include/io/bst_file.cpp
             ${CMAKE_SOURCE_DIR}/include/io/map_file.cpp
             ${CMAKE_SOURCE_DIR}/include/io/dnaiotpb.cpp
             ${CMAKE_SOURCE_DIR}/include/io/seg_file.cpp
             ${CMAKE_SOURCE_DIR}/include/parameters/dnadatum.cpp
             ${CMAKE_SOURCE_DIR}/include/parameters/dnaellipsoid.cpp
             ${CMAKE_SOURCE_DIR}/include/parameters/dnaprojection.cpp
             ${CMAKE_SOURCE_DIR}/include/functions/dnaprocessfuncs.cpp
             ${CMAKE_SOURCE_DIR}/include/functions/dnastringfuncs.cpp
             ${CMAKE_SOURCE_DIR}/include/math/dnamatrix_contiguous.cpp
             ${CMAKE_SOURCE_DIR}/include/ide/trace.cpp
             ${CMAKE_SOURCE_DIR}/include/measurement_types/dnamsrtally.cpp
             ${CMAKE_SOURCE_DIR}/include/measurement_types/dnastation.cpp
             ${CMAKE_SOURCE_DIR}/include/io/dynadjust_file.cpp
             dnaplot.cpp
             ${CMAKE_SOURCE_DIR}/dynadjust.rc)

target_precompile_headers (${PROJECT_NAME} PUBLIC precompile.h)

target_link_libraries (${PROJECT_NAME} ${DNA_LIBRARIES})

install (TARGETS ${PROJECT_NAME} DESTINATION lib)

# Install to the x64 directory
set_target_properties (${PROJECT_NAME} PROPERTIES 
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/x64")    
    
install (TARGETS ${PROJECT_NAME} DESTINATION x64)
