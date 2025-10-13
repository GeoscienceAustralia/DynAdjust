# <dnaadjustwrapper/...> build rules
project (dnaadjustwrapper)

add_subdirectory (${CMAKE_SOURCE_DIR}/dynadjust/dnaadjust ${CMAKE_BINARY_DIR}/dynadjust/dnaadjust)

add_definitions (-DBUILD_ADJUST_EXE)

include_directories (${PROJECT_SOURCE_DIR})

if(NOT BUILD_STATIC)
  add_executable (${PROJECT_NAME}
                  ${CMAKE_SOURCE_DIR}/include/config/dnaprojectfile.cpp
                  ${CMAKE_SOURCE_DIR}/include/functions/dnastringfuncs.cpp
                  dnaadjustprogress.cpp
                  dnaadjustwrapper.cpp
                  ${CMAKE_SOURCE_DIR}/dynadjust.rc)

  target_precompile_headers (${PROJECT_NAME} PUBLIC precompile.h)

  target_link_libraries (${PROJECT_NAME} dnaadjust ${DNA_LIBRARIES})

  # Link C++ filesystem library if needed
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "9.0")
      target_link_libraries (${PROJECT_NAME} stdc++fs)
    endif()
  endif()

  set_target_properties (${PROJECT_NAME} PROPERTIES OUTPUT_NAME "${DNA_PROGRAM_PREFIX}adjust")
endif()

if(BUILD_STATIC)

  # Define source files for static build
  set(TARGET_SOURCES
    ${CMAKE_SOURCE_DIR}/include/config/dnaprojectfile.cpp
    ${CMAKE_SOURCE_DIR}/include/functions/dnastringfuncs.cpp
    dnaadjustprogress.cpp
    dnaadjustwrapper.cpp
    ${CMAKE_SOURCE_DIR}/dynadjust.rc)
  
  # Create the static executable with _static suffix
  set(STATIC_TARGET_NAME "${PROJECT_NAME}_static")
  add_executable(${STATIC_TARGET_NAME} ${TARGET_SOURCES})
  
  # Add dnaadjust library sources directly
  get_target_property(THE_SOURCES dnaadjust SOURCES)
  get_target_property(THE_SOURCE_DIR dnaadjust SOURCE_DIR)
  
  # Add these sources directly to the static executable
  foreach(src ${THE_SOURCES})
    # Check if it's a relative path
    if(NOT IS_ABSOLUTE ${src})
      set(src "${THE_SOURCE_DIR}/${src}")
    endif()
    message(STATUS "Adding source ${src}")
    # Add the source to the static target
    target_sources(${STATIC_TARGET_NAME} PRIVATE ${src})
  endforeach()

  # When building static, exclude the library from the build
  set_target_properties(dnaadjust PROPERTIES EXCLUDE_FROM_ALL 1 EXCLUDE_FROM_DEFAULT_BUILD 1)
  
  print_static_linking_info(${STATIC_TARGET_NAME} "${DNA_PROGRAM_PREFIX}adjust")
  
  # Link with static libraries
  target_link_libraries(${STATIC_TARGET_NAME} PRIVATE ${STATIC_LIBS_LIST})

  # Apply static build optimizations
  optimize_static_target(${STATIC_TARGET_NAME})
  
  if(UNIX)
    if(APPLE)
      # On macOS, link against system libraries for curl
      #find_package(CURL REQUIRED)

      # On macOS, true static linking is challenging but we can improve it
      # Ensure we're not using @rpath for the DynAdjust libraries
      set_target_properties(${STATIC_TARGET_NAME} PROPERTIES
        BUILD_WITH_INSTALL_RPATH FALSE
        INSTALL_RPATH ""
        SKIP_BUILD_RPATH TRUE)
      
      # Request static linking to the extent possible on macOS
      target_link_options(${STATIC_TARGET_NAME} PRIVATE
        "-Wl,-search_paths_first"
        "-Wl,-dead_strip"
        "-Wl,-all_load"
        "-Wl,-dead_strip_dylibs")
     
      target_link_libraries(${STATIC_TARGET_NAME} PRIVATE
        #CURL::libcurl
        "-framework Accelerate"
        "-framework CoreFoundation")
    else()
      target_link_options(${STATIC_TARGET_NAME} PRIVATE
        -static)

    endif()
  elseif(WIN32)
    # On Windows with MSVC
    if(MSVC)
      set_property(TARGET ${STATIC_TARGET_NAME} PROPERTY MSVC_RUNTIME_LIBRARY "MultiThreaded")
    endif()
  endif()
  
  # Set the output name to be the same as the dynamic version
  set_target_properties(${STATIC_TARGET_NAME}
    PROPERTIES OUTPUT_NAME "${DNA_PROGRAM_PREFIX}adjust")

endif()

if(NOT BUILD_STATIC)
  set(DNAADJUST_TARGET ${PROJECT_NAME} PARENT_SCOPE)
else()
  set(DNAADJUST_TARGET ${STATIC_TARGET_NAME} PARENT_SCOPE)
endif()
