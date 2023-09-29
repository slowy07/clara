cmake_minimum_required(VERSION 3.0.0)
project(
  clara
  VERSION 1.0
  LANGUAGES CXX)

include_directories(include)

set(SOURCE_FILES testing/minimal.cpp)

message(STATUS "Detecting Eigen3")
set(EIGEN3_INCLUDE_DIR "" CACHE PATH "Path to Eigen3")
if(NOT ${EIGEN3_INCLUDE_DIR} STREQUAL "")
  if(IS_DIRECTORY ${EIGEN3_INCLUDE_DIR})
    message(STATUS "Detecting Eigen3 - done (in ${EIGEN3_INCLUDE_DIR})")
    include_directories(SYSTEM "${EIGEN3_INCLUDE_DIR}")
  else()
    message(FATAL_ERROR "Invalid path to Eigen3 installation")
  endif()
else()
  find_package(Eigen3 3.0 QUIET NO_MODULE)
  if(NOT TARGET Eigen3::Eigen)
    message(FATAL_ERROR "eigen 3 not detected please point EIGEN3_INCLUDE_DIR")
  endif()
  message(STATUS "Detecting Eigen3 - done")
endif()

if(CYGWIN)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++11")
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
endif()

set(ADDITIONAL_FLAGS "-pedantic -Wall -Wextra -Weffc++")
set(CMAKE_INCLUDE_SYSTEM_FLAG_CXX "-isystem")

if(NOT ${CMAKE_CXX_COMPILER_ID} MATCHES "Clang"
   AND NOT ${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
  message(
    FATAL_ERROR "this CMakeLists.txt work only with GNU gcc/Clang/AppleCLang")
endif()

if(${CMAKE_CXX_COMPILER_ID} STREQUAL "AppleClang")
  add_definitions(-DNO_THREAD_LOCAL_)
  message(
    WARNING
      "Detected compiler ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}. thread_local not supported"
  )
endif()

option(WITH_OPENMP "OpenMP support" ON)
if(${WITH_OPENMP})
  if(${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU"
     OR (${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang"
         AND (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER "3.8"
              OR CMAKE_CXX_COMPILER_VERSION VERSION_EQUAL "3.8")))
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
    add_definitions(-DWITH_OPENMP)
  else()
    message(
      WARNING
        "Detected compiler : ${CMAKE_CXX_COMPILER_VERSION}, support for OpenMP is enabled only for the gnu gcc compiler or the clang version 3.8 higher"
    )
  endif()
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ADDITIONAL_FLAGS}")

if(${CMAKE_CXX_FLAGS} MATCHES "GNU")
  if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fno-weak")
  endif()
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Og -D_GLIBCXX_DEBUG")
endif()

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG")
set(CMAKE_CXX_FLAGS_MINSIZEREL "${CMAKE_CXX_FLAGS_MINSIZEREL} -DEIGEN_NO_DEBUG")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO
    "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -DEIGEN_NO_DEBUG")

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE
      Release
      CACHE
        STRING
        "choose tyhe type build, oprtions are: None Debug Release MinSizeRel RelWithDebInfo"
        FORCE)
endif()

add_executable(clara ${SOURCE_FILES})

if(TARGET Eigen3::Eigen)
  target_link_libraries(clara Eigen3::Eigen)
endif()

if($WITH_OPENMP$
   AND ${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang"
   AND CMAKE_CXX_COMPILER_VERSION VERSION_GREATER "3.7")
  target_link_libraries(clara omp)
endif()
