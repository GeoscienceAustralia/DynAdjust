# ./dynadjust/dynadjust/dnadiff/CMakeLists.txt
CMAKE_MINIMUM_REQUIRED(VERSION 3.16)
SET(CMAKE_CXX_STANDARD 17)

PROJECT (dnadiff CXX)

SET (CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/../../../../cmake" ${CMAKE_MODULE_PATH})

if(NOT BUILD_STATIC)
  add_executable (${PROJECT_NAME}
                  dnadiff.cpp)

  target_link_libraries (${PROJECT_NAME})

  # Link C++ filesystem library if needed
  if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "9.0")
      target_link_libraries (${PROJECT_NAME} stdc++fs)
    endif()
  elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    # Only needed on Linux with Clang (not on macOS/Apple Clang)
    if(UNIX AND NOT APPLE)
      target_link_libraries (${PROJECT_NAME} stdc++fs)
    endif()
  endif()

else()
  SET (STATIC_TARGET_NAME "${PROJECT_NAME}_static")

  add_executable(${STATIC_TARGET_NAME}
                dnadiff.cpp)

  # Always output as 'dnadiff' (not 'dnadiff_static')
  SET_TARGET_PROPERTIES(${STATIC_TARGET_NAME} PROPERTIES OUTPUT_NAME dnadiff)

  # Link C++ filesystem library if needed
  if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "9.0")
      target_link_libraries (${STATIC_TARGET_NAME} stdc++fs)
    endif()
  elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    # Only needed on Linux with Clang (not on macOS/Apple Clang)
    if(UNIX AND NOT APPLE)
      target_link_libraries (${STATIC_TARGET_NAME} stdc++fs)
    endif()
  endif()

  # Link statically
  if(UNIX)
    if(NOT APPLE)
      SET_TARGET_PROPERTIES(${STATIC_TARGET_NAME} PROPERTIES LINK_FLAGS "-static")
      target_link_libraries (${STATIC_TARGET_NAME} pthread)
    endif()
  endif()

endif()

if(NOT BUILD_STATIC)
  install (TARGETS ${PROJECT_NAME} DESTINATION bin)
else()
  install (TARGETS ${STATIC_TARGET_NAME} DESTINATION bin)
endif()

# Also export the target name for use in tests
if(NOT BUILD_STATIC)
  set(DNADIFF_TARGET ${PROJECT_NAME} PARENT_SCOPE)
else()
  set(DNADIFF_TARGET ${STATIC_TARGET_NAME} PARENT_SCOPE)
endif()