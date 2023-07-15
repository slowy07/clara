cmake_minimum_required(VERSION 3.0.0)
project(clara_testing)

if(CYGWIN)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++11")
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
endif()

set(ADDITIONAL_FLAGS "-pedantic -Wall -Wextra -Weffc++")
set(CMAKE_INCLUDE_SYSTEM_FLAG_CXX "-isystem")

if(NOT ${CMAKE_CXX_COMPILER_ID} MATCHES "Clang"
   AND NOT ${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
  message(FATAL_ERROR "cmakelists works only with gnu gcc/clang/appleClang")
endif()

include_directories(../include)
# eigen
include_directories(SYSTEM "$ENV{HOME}/eigen")

option(WITH_OPENMP "OpenMP support" ON)
if(${WITH_OPENMP})
  add_definitions(-DWITH_OPENMP_)
  if(${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")
  elseif(${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
    if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER "3.7")
      set(CLANG_LIBOMP "/opt/local/lib/libomp/")
      set(CLANG_LIBOMP_INCLUDE "/opt/local/include/libomp/")
      include_directories(SYSTEM "${CLANG_LIBOMP_INCLUDE}")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp=libomp")
      set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -L${CLANG_LIBOMP}")
    endif()
  else()
    message(
      WARNING
        "compiler: ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}. support for OpenMP is enabled only for GNU gcc compiler Clang >= 3.7"
    )
  endif()
endif()

# google test
add_subdirectory(googletest)

# unitesting
add_subdirectory(tests)

# enable all warnings
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ADDITIONAL_FLAGS}")

if(${CMAKE_CXX_COMPILER_ID} STREQUAL "AppleClang")
  add_definitions(-DNO_THREAD_LOCAL_)
  message(
    WARNING
      "Detect compiler: ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION} thread_local not supported"
  )
endif()

# gnu gcc aditional debug setting
if(${CMAKE_CXX_COMPILER_ID} MATCHES "GNU")
  if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fno-weak")
  endif()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_DEBUG} -Og -D_GLIBCXX_DEBUG")
endif()

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG")
set(CMAKE_CXX_FLAGS_MINSIZEREL "${CMAKE_CXX_FLAGS_MINSIZEREL} -DEIGEN_NO_DEBUG")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO
    "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -DEIGEN_NO_DEBUG")

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE
      Release
      CACHE STRING "Choose the type of build, options are: \
         None Debug Release MinSizeRel RelWithDebInfo." FORCE)
endif()

if($WITH_OPENMP$
   AND ${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang"
   AND CMAKE_CXX_COMPILER_VERSION VERSION_GREATER "3.7")
  target_link_libraries(clara_testing omp)
endif()
