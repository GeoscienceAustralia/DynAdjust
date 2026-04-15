# <dynadjust/...> build rules
project (dynadjust)

add_definitions (-DBUILD_DYNADJUST_EXE)

include_directories (${PROJECT_SOURCE_DIR})

if(NOT BUILD_STATIC)
  add_executable (${PROJECT_NAME} 
                  ${CMAKE_SOURCE_DIR}/include/config/dnaprojectfile.cpp
                  ${CMAKE_SOURCE_DIR}/include/functions/dnaprocessfuncs.cpp
                  ${CMAKE_SOURCE_DIR}/include/functions/dnastringfuncs.cpp
                  dynadjust.cpp
                  ${CMAKE_SOURCE_DIR}/dynadjust.rc)

  target_precompile_headers (${PROJECT_NAME} PUBLIC precompile.h)

  target_link_libraries (${PROJECT_NAME} ${DNA_LIBRARIES})

  set_target_properties (${PROJECT_NAME} PROPERTIES OUTPUT_NAME "dynadjust")
endif()

if(BUILD_STATIC)

  # Define source files for static build
  set(TARGET_SOURCES
    ${CMAKE_SOURCE_DIR}/include/config/dnaprojectfile.cpp
    ${CMAKE_SOURCE_DIR}/include/functions/dnaprocessfuncs.cpp
    ${CMAKE_SOURCE_DIR}/include/functions/dnastringfuncs.cpp
    dynadjust.cpp
    ${CMAKE_SOURCE_DIR}/dynadjust.rc)
  
  # Create the static executable with _static suffix
  set(STATIC_TARGET_NAME "${PROJECT_NAME}_static")
  add_executable(${STATIC_TARGET_NAME} ${TARGET_SOURCES})
  
  # Note: dynadjust is just an executable wrapper, not a library
  # No additional sources to add

  # When building static, exclude the library from the build
  
  print_static_linking_info(${STATIC_TARGET_NAME} "${DNA_PROGRAM_PREFIX}dynadjust")
  
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
        "-framework CoreFoundation")
    else()
      target_link_options(${STATIC_TARGET_NAME} PRIVATE
        -static-libgcc
        -static-libstdc++)
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

endif()

# Export the target name for use in tests
if(NOT BUILD_STATIC)
  set(DYNADJUST_TARGET ${PROJECT_NAME} PARENT_SCOPE)
else()
  set(DYNADJUST_TARGET ${STATIC_TARGET_NAME} PARENT_SCOPE)
endif()
