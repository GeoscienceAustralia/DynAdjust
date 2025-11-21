# <dnaplot/...> build rules
project (dnaplot)

add_definitions(-DBUILD_DNAPLOT_DLL -DBUILD_EXPORTDNATYPES -DMKL_ILP64)

include_directories (${PROJECT_SOURCE_DIR})

add_library (${PROJECT_NAME} SHARED
             ${CMAKE_SOURCE_DIR}/include/io/dnaiobms.cpp
             ${CMAKE_SOURCE_DIR}/include/io/dnaiobst.cpp
             ${CMAKE_SOURCE_DIR}/include/io/dnaiomap.cpp
             ${CMAKE_SOURCE_DIR}/include/io/dnaiotpb.cpp
             ${CMAKE_SOURCE_DIR}/include/io/dnaioseg.cpp
             ${CMAKE_SOURCE_DIR}/include/parameters/dnadatum.cpp
             ${CMAKE_SOURCE_DIR}/include/parameters/dnaellipsoid.cpp
             ${CMAKE_SOURCE_DIR}/include/parameters/dnaprojection.cpp
             ${CMAKE_SOURCE_DIR}/include/functions/dnastringfuncs.cpp
             ${CMAKE_SOURCE_DIR}/include/functions/dnaprocessfuncs.cpp
             ${CMAKE_SOURCE_DIR}/include/math/dnamatrix_contiguous.cpp
             ${CMAKE_SOURCE_DIR}/include/measurement_types/dnamsrtally.cpp
             ${CMAKE_SOURCE_DIR}/include/measurement_types/dnastation.cpp
             ${CMAKE_SOURCE_DIR}/include/io/dnaiobase.cpp
             dnaplot.cpp
             ${CMAKE_SOURCE_DIR}/dynadjust.rc)

target_link_libraries (${PROJECT_NAME} ${DNA_LIBRARIES})

install (TARGETS ${PROJECT_NAME} DESTINATION lib)
