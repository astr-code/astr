cmake_minimum_required(VERSION 3.0.2)                        
cmake_policy(SET CMP0048 NEW)

project(ASTR LANGUAGES Fortran)

set(AUTHOR "Jian Fang")
set(AUTHOR_DETAILS "jian.fang@stfc.ac.uk")
set(DESCRIPTION "Building ASTR using cmake")

message(STATUS "building ${PROJECT_NAME}")

include(GNUInstallDirs)
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR})
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "${PROJECT_BINARY_DIR}/opt" CACHE PATH "..." FORCE)
endif()

# Add support for CMAKE_DEPENDENT_OPTION
INCLUDE(CMakeDependentOption)
INCLUDE(CMakeParseArguments)

# Find the modules included with ASTR
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")

message(STATUS "CMAKE_MODULE_PATH:${CMAKE_MODULE_PATH}")

# make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE)
endif (NOT CMAKE_BUILD_TYPE)

if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "PGI" OR
   ${CMAKE_Fortran_COMPILER_ID} STREQUAL "NVHPC")
 
  set(CMAKE_Fortran_FLAGS "-cpp")
  set(CMAKE_Fortran_FLAGS_DEBUG "-g -O0 -traceback -Mbounds -Mchkptr -Ktrap=fp")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -fast")

elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
  set(CMAKE_Fortran_FLAGS "-cpp -std=f2018")
  set(CMAKE_Fortran_FLAGS_DEBUG "-g -Og -Wall -Wpedantic -Werror -Wimplicit-interface -Wimplicit-procedure -Wno-unused-dummy-argument")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -ffast-math")

  message(STATUS "CMAKE_Fortran_FLAGS: ${CMAKE_Fortran_FLAGS}")
endif()

if (CMAKE_BUILD_TYPE MATCHES "DEBUG")
  add_definitions("-DDEBG")
endif (CMAKE_BUILD_TYPE MATCHES "DEBUG")

set(ASTR_ROOT ${CMAKE_CURRENT_SOURCE_DIR})
message(STATUS "ASTR directory: ${ASTR_ROOT}")
set(CMAKE_INSTALL_PREFIX ${ASTR_ROOT})

# Create the ASTR executable
add_subdirectory(tgv_solver_3d)
