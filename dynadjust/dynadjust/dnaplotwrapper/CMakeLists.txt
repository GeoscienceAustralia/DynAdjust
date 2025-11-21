# <dnaplotwrapper/...> build rules
project (dnaplotwrapper)

# build dnaplot library first
add_subdirectory (${CMAKE_SOURCE_DIR}/dynadjust/dnaplot ${CMAKE_BINARY_DIR}/dynadjust/dnaplot)

add_definitions (-DBUILD_PLOT_EXE -DMKL_ILP64)

include_directories (${PROJECT_SOURCE_DIR})

add_executable (${PROJECT_NAME} 
			${CMAKE_SOURCE_DIR}/include/io/dnaiobase.cpp
			${CMAKE_SOURCE_DIR}/include/io/dnaiobst.cpp
			${CMAKE_SOURCE_DIR}/include/io/dnaiobms.cpp
			${CMAKE_SOURCE_DIR}/include/functions/dnastringfuncs.cpp
			${CMAKE_SOURCE_DIR}/include/functions/dnaprocessfuncs.cpp
			${CMAKE_SOURCE_DIR}/include/parameters/dnaellipsoid.cpp
			${CMAKE_SOURCE_DIR}/include/parameters/dnaprojection.cpp
			dnaplotwrapper.cpp
			${CMAKE_SOURCE_DIR}/dynadjust.rc)

target_link_libraries (${PROJECT_NAME} dnaplot ${DNA_LIBRARIES})

set_target_properties (${PROJECT_NAME} PROPERTIES OUTPUT_NAME "${DNA_PROGRAM_PREFIX}plot")

install (TARGETS ${PROJECT_NAME} DESTINATION bin)
