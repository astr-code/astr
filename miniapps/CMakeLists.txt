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

set(Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER_ID} )
message(STATUS "COMP ID ${Fortran_COMPILER_NAME}")
message(STATUS "Fortran compiler name ${Fortran_COMPILER_NAME}")
message(STATUS "Fortran compiler version ${CMAKE_Fortran_COMPILER_VERSION}")

FIND_PACKAGE(OpenMP_Fortran)
set(OMP_COMPILE_FLAGS "${OpenMP_Fortran_FLAGS}")
set(OMP_LINK_FLAGS "${OpenMP_Fortran_FLAGS}")
message(STATUS "OMP_COMPILE_FLAGS:" ${OMP_COMPILE_FLAGS})
message(STATUS "OMP_LINK_FLAGS:" ${OMP_LINK_FLAGS})

if (Fortran_COMPILER_NAME MATCHES "GNU")
  # gfortran
  message(STATUS "Setting gfortran flags")
  set(CMAKE_Fortran_FLAGS "-cpp -funroll-loops -floop-optimize -g -fcray-pointer -fbacktrace -ffree-line-length-none")
  if (CMAKE_Fortran_COMPILER_VERSION GREATER_EQUAL "10")
    message(STATUS "Set New Fortran basic flags")
    set(CMAKE_Fortran_FLAGS "-cpp -funroll-loops -floop-optimize -g -fcray-pointer -fbacktrace -ffree-line-length-none -fallow-argument-mismatch ${OMP_COMPILE_FLAGS}" )
  endif (CMAKE_Fortran_COMPILER_VERSION GREATER_EQUAL "10")
  set(CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3 ${OMP_LINK_FLAGS}")
  set(CMAKE_Fortran_FLAGS_DEBUG   "-DDEBG -O0 -g ${OMP_LINK_FLAGS}")
elseif (Fortran_COMPILER_NAME MATCHES "Intel")
  message(STATUS "Setting ifort flags")
  set(CMAKE_Fortran_FLAGS "-fpp -xHost -heaparrays -safe-cray-ptr -g -traceback")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O3 -ipo")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -DDEBG")
  #set(CMAKE_Fortran_FLAGS "-cpp xSSE4.2 -axAVX,CORE-AVX-I,CORE-AVX2 -ipo -fp-model fast=2 -mcmodel=large -safe-cray-ptr")
elseif (Fortran_COMPILER_NAME MATCHES "NAG")
  message(STATUS "Setting nagfor flags")
  set(CMAKE_Fortran_FLAGS "-fpp")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O3")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
elseif (Fortran_COMPILER_NAME MATCHES "Cray")
  message(STATUS "Setting cray fortran flags")
  set(CMAKE_Fortran_FLAGS "-eF -g -N 1023")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O3")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
elseif (Fortran_COMPILER_NAME MATCHES "PGI")
  message(STATUS "Setting PGI fortran flags")
  set(CMAKE_Fortran_FLAGS "-cpp -acc -Mfree -Kieee -Minfo=accel -g")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O3")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -DDEBG")
elseif (Fortran_COMPILER_NAME MATCHES "Fujitsu")
  message(STATUS "Setting Fujitsu fortran flags")
  set (CMAKE_Fortran_FLAGS "-Cpp")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O3")
  set (CMAKE_Fortran_FLAGS_DEBUG "-O0")
else (Fortran_COMPILER_NAME MATCHES "GNU")
  message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
  message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
  message ("No optimized Fortran compiler flags are known, we just try -O2...")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
endif (Fortran_COMPILER_NAME MATCHES "GNU")

if (CMAKE_BUILD_TYPE MATCHES "DEBUG")
  add_definitions("-DDEBG")
endif (CMAKE_BUILD_TYPE MATCHES "DEBUG")

set(ASTR_ROOT ${CMAKE_CURRENT_SOURCE_DIR})
message(STATUS "ASTR directory: ${ASTR_ROOT}")
set(CMAKE_INSTALL_PREFIX ${ASTR_ROOT})

# Create the ASTR executable
add_subdirectory(tgv_solver_3d)

