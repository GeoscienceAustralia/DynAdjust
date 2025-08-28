# <dynadjust/...> build rules
project (dynadjust)

add_definitions (-DBUILD_DYNADJUST_EXE)

include_directories (${PROJECT_SOURCE_DIR})

add_executable (${PROJECT_NAME} 
                ${CMAKE_SOURCE_DIR}/include/config/dnaprojectfile.cpp
                ${CMAKE_SOURCE_DIR}/include/functions/dnaprocessfuncs.cpp
                ${CMAKE_SOURCE_DIR}/include/functions/dnastringfuncs.cpp
                dynadjust.cpp
                ${CMAKE_SOURCE_DIR}/dynadjust.rc)

target_precompile_headers (${PROJECT_NAME} PUBLIC precompile.h)

target_link_libraries (${PROJECT_NAME} ${DNA_LIBRARIES})

set_target_properties (${PROJECT_NAME} PROPERTIES OUTPUT_NAME "dynadjust")



if(BUILD_STATIC)

  # Create a static version with the same source files
  get_target_property(TARGET_SOURCES ${PROJECT_NAME} SOURCES)
  
  # Create the static executable with _static suffix
  set(STATIC_TARGET_NAME "${PROJECT_NAME}_static")
  add_executable(${STATIC_TARGET_NAME} ${TARGET_SOURCES})
  
  # Add dnadynadjust library sources directly
  get_target_property(DNAREFTRAN_SOURCES dynadjust SOURCES)
  get_target_property(DNAREFTRAN_SOURCE_DIR dynadjust SOURCE_DIR)
  
  # Add these sources directly to the static executable
  foreach(src ${DNAREFTRAN_SOURCES})
    # Check if it's a relative path
    if(NOT IS_ABSOLUTE ${src})
      set(src "${DNAREFTRAN_SOURCE_DIR}/${src}")
    endif()
    message(STATUS "Adding source ${src}")
    # Add the source to the static target
    target_sources(${STATIC_TARGET_NAME} PRIVATE ${src})
  endforeach()

  # Speed things up by not bothering with normal build
  set_target_properties(${PROJECT_NAME} PROPERTIES EXCLUDE_FROM_ALL 1 EXCLUDE_FROM_DEFAULT_BUILD 1)
  
  # Link with static libraries
  target_link_libraries(${STATIC_TARGET_NAME} PRIVATE ${STATIC_LIBS_LIST})
  
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
        "-framework CoreFoundation")
    else()
      # On Linux/Unix, we can use full static linking
      target_link_options(${STATIC_TARGET_NAME} PRIVATE
        -Wl,-Bstatic
        -fno-builtin
        -static-libgcc
        -static-libstdc++
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
    PROPERTIES OUTPUT_NAME "dynadjust")
  
  # Also set runtime output directory to the static subdirectory
  set_target_properties(${STATIC_TARGET_NAME} PROPERTIES 
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/static")
    
endif()
