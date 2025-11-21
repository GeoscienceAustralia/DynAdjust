# <dnaadjustwrapper/...> build rules
project (dnaadjustwrapper)

# build dnaadjust library first
add_subdirectory (${CMAKE_SOURCE_DIR}/dynadjust/dnaadjust ${CMAKE_BINARY_DIR}/dynadjust/dnaadjust)

add_definitions (-DBUILD_ADJUST_EXE -DMULTI_THREAD_ADJUST -DMKL_ILP64)

include_directories (${PROJECT_SOURCE_DIR})

add_executable (${PROJECT_NAME}
                ${CMAKE_SOURCE_DIR}/include/config/dnaprojectfile.cpp
                ${CMAKE_SOURCE_DIR}/include/io/dnaiobase.cpp
                ${CMAKE_SOURCE_DIR}/include/io/dnaiobms.cpp
                ${CMAKE_SOURCE_DIR}/include/io/dnaiobst.cpp
                ${CMAKE_SOURCE_DIR}/include/parameters/dnadatum.cpp
                ${CMAKE_SOURCE_DIR}/include/functions/dnastringfuncs.cpp
                ${CMAKE_SOURCE_DIR}/include/parameters/dnaellipsoid.cpp
                ${CMAKE_SOURCE_DIR}/include/parameters/dnaprojection.cpp
                dnaadjustprogress.cpp
                dnaadjustwrapper.cpp
                ${CMAKE_SOURCE_DIR}/dynadjust.rc)

target_link_libraries (${PROJECT_NAME} dnaadjust ${DNA_LIBRARIES})

set_target_properties (${PROJECT_NAME} PROPERTIES OUTPUT_NAME "${DNA_PROGRAM_PREFIX}adjust")

install (TARGETS ${PROJECT_NAME} DESTINATION bin)
