cmake_minimum_required(VERSION 3.0.0)
project(
  clara
  VERSION 1.0
  LANGUAGES CXX)

# C++11
if(CYGWIN)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++11")
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
endif()

set(ADDITIONAL_FLAGS "-pedantic -Wall -Wextra -Weffc++")
set(CMAKE_INCLUDE_SYSTEM_FLAG_CXX "-isystem")

# CMake suppport only for GNUC gcc and Clang and apple clang
if(NOT ${CMAKE_CXX_COMPILER_ID} MATCHES "Clang"
   AND NOT ${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
  message(FATAL_ERROR "CMakeLists.txt only work with GNU gcc/Clang/AppleClang")
endif()

# clara header
include_directories(include)
set(CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} include)

# source file
set(SOURCE_FILES testing/_test.cpp)

include_directories(SYSTEM "$ENV{HOME}/eigen3")

# openmp
option(WITH_OPENMP "OpenMP support" ON)
if(${WITH_OPENMP})
  # inject definition (as #define) in the source file
  add_definitions(-DWITH_OPENMP_)
  # gnuc
  if(${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")
    # clang but not appleClang
  elseif(${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
    if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER "3.7")
      set(CLANG_LIBOMP "opt/local/lib/libomp/")
      set(CLANG_LIBOMP_INCLUDE "/opt/local/include/libomp/")
      include_directories(SYSTEM "${CLANG_LIBOMP_INCLUDE}")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp=libomp")
      set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -L${CLANG_LIBOMP}")
    endif()
  else()
    message(
      WARNING
        "detected compiler: ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION} A support for openMP is enable only the GNU gcc compiler or Clang >= 3.7"
    )
  endif()
endif()

# enable all warnings
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ADDITIONAL_FLAGS}")

# disable support for thread_local storage duration specifier when using apple
# clang as libc++ doesn't yet support
if(${CMAKE_CXX_COMPILER_ID} STREQUAL "AppleClang")
  # inject definition (as #define) in the source files
  add_definitions(-DNO_THREAD_LOCAL_)
  message(
    WARNING
      "detected compiler: ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION} thread_local not support"
  )
endif()

# gnu gcc additional debug settings
if(${CMAKE_CXX_COMPILER_ID} MATCHES "GNU")
  if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fno-weak")
  endif()
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Og -D_GLIBCXX_DEBUG")
endif()

# configuration
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG")
set(CMAKE_CXX_FLAGS_MINSIZEREL "${CMAKE_CXX_FLAGS_MINSIZEREL} -DEIGEN_NO_DEBUG")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -DEIGEN_NO_DEBUG")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO
    "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -DEIGEN_NO_DEBUG")

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE
      Release
      CACHE
        STRING
        "choose the type of build, options are: None Debug Release miniSizeRel RelWithDebugInfo"
        FORCE)
endif()

add_executable(clara ${SOURCE_FILES})

if($WITH_OPENMP$
   AND ${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang"
   AND CMAKE_CXX_COMPILER_VERSION VERSION_GREATER "3.7")
  target_link_libraries(clara omp)
endif()
