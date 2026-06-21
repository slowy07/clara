cmake_minimum_required(VERSION 3.10...3.31)
project(
  clara
  VERSION 1.0
  LANGUAGES CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# ── Sources (all files for LSP compile_commands) ───────────────────────────
file(GLOB SOURCE_FILES testing/*.cpp)
list(REMOVE_ITEM SOURCE_FILES "${CMAKE_SOURCE_DIR}/testing/_test.cpp")

# ── C++ Standard ───────────────────────────────────────────────────────────
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if(CYGWIN)
  set(CMAKE_CXX_EXTENSIONS ON)
endif()

# ── Compiler check ─────────────────────────────────────────────────────────
if(NOT CMAKE_CXX_COMPILER_ID MATCHES "(Clang|GNU|AppleClang)")
  message(FATAL_ERROR "this project works only with GNU gcc/Clang/AppleClang")
endif()

# ── Warnings ───────────────────────────────────────────────────────────────
add_library(clara_warnings INTERFACE)
target_compile_options(clara_warnings INTERFACE
  -pedantic -Wall -Wextra -Weffc++)

# ── Eigen3 ─────────────────────────────────────────────────────────────────
message(STATUS "Detecting Eigen3")
set(EIGEN3_INCLUDE_DIR "" CACHE PATH "Path to Eigen3")
if(EIGEN3_INCLUDE_DIR)
  if(NOT IS_DIRECTORY "${EIGEN3_INCLUDE_DIR}")
    message(FATAL_ERROR "Invalid path to Eigen3 installation")
  endif()
  add_library(clara_eigen INTERFACE)
  target_include_directories(clara_eigen SYSTEM INTERFACE "${EIGEN3_INCLUDE_DIR}")
  message(STATUS "Detecting Eigen3 - done (in ${EIGEN3_INCLUDE_DIR})")
else()
  find_package(Eigen3 3.0 QUIET NO_MODULE)
  if(NOT TARGET Eigen3::Eigen)
    message(FATAL_ERROR "Eigen3 not detected -- set EIGEN3_INCLUDE_DIR or install")
  endif()
  message(STATUS "Detecting Eigen3 - done")
endif()

# ── OpenMP ─────────────────────────────────────────────────────────────────
option(WITH_OPENMP "OpenMP support" ON)
if(WITH_OPENMP)
  if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU"
     OR (CMAKE_CXX_COMPILER_ID STREQUAL "Clang"
         AND CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL "3.8"))
    find_package(OpenMP QUIET)
    if(OpenMP_CXX_FOUND)
      add_library(clara_openmp INTERFACE)
      target_compile_definitions(clara_openmp INTERFACE WITH_OPENMP)
      target_link_libraries(clara_openmp INTERFACE OpenMP::OpenMP_CXX)
    else()
      message(WARNING "OpenMP requested but not found -- continuing without it")
    endif()
  else()
    message(WARNING "OpenMP not supported for ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
  endif()
endif()

# ── AppleClang workaround ──────────────────────────────────────────────────
if(CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
  add_definitions(-DNO_THREAD_LOCAL_)
  message(WARNING "AppleClang: thread_local not supported")
endif()

# ── Library of common dependencies ────────────────────────────────────────
add_library(clara_deps INTERFACE)
target_include_directories(clara_deps BEFORE INTERFACE include)

if(TARGET clara_eigen)
  target_link_libraries(clara_deps INTERFACE clara_eigen)
elseif(TARGET Eigen3::Eigen)
  target_link_libraries(clara_deps INTERFACE Eigen3::Eigen)
endif()

if(TARGET clara_openmp)
  target_link_libraries(clara_deps INTERFACE clara_openmp)
endif()

target_link_libraries(clara_deps INTERFACE clara_warnings)

# ── Executables (one per test file) ───────────────────────────────────────
foreach(FILE ${SOURCE_FILES})
  get_filename_component(NAME ${FILE} NAME_WE)
  add_executable(${NAME} ${FILE})
  target_link_libraries(${NAME} PRIVATE clara_deps)
endforeach()

# ── Debug flags ────────────────────────────────────────────────────────────
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU" AND CMAKE_SYSTEM_NAME MATCHES "Darwin")
  string(APPEND CMAKE_CXX_FLAGS_DEBUG " -fno-weak")
endif()
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  string(APPEND CMAKE_CXX_FLAGS_DEBUG " -Og -D_GLIBCXX_DEBUG")
endif()
string(APPEND CMAKE_CXX_FLAGS_DEBUG " -DDEBUG")
string(APPEND CMAKE_CXX_FLAGS_MINSIZEREL " -DEIGEN_NO_DEBUG")
string(APPEND CMAKE_CXX_FLAGS_RELWITHDEBINFO " -DEIGEN_NO_DEBUG")

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE
      Release
      CACHE STRING "build type: None Debug Release MinSizeRel RelWithDebInfo" FORCE)
endif()
