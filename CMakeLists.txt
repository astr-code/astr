cmake_minimum_required(VERSION 3.16)

project(ASTR
  VERSION 2.6.0
  DESCRIPTION "High-order finite-difference solver for compressible turbulence"
  LANGUAGES Fortran)

include(GNUInstallDirs)
include(CTest)

option(CHEMISTRY "Build ASTR with chemistry" OFF)
option(ASTR_MPI_OVERSUBSCRIBE "Pass --oversubscribe to MPI test runs" ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the build type" FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release RelWithDebInfo MinSizeRel)
endif()

set(CMAKE_Fortran_MODULE_DIRECTORY "${PROJECT_BINARY_DIR}/mod")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}")

message(STATUS "Building ${PROJECT_NAME} ${PROJECT_VERSION}")
message(STATUS "Fortran compiler: ${CMAKE_Fortran_COMPILER_ID} ${CMAKE_Fortran_COMPILER_VERSION}")

add_library(astr_build_options INTERFACE)

if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  target_compile_options(astr_build_options INTERFACE
    -cpp
    -funroll-loops
    -floop-optimize
    -fcray-pointer
    -fbacktrace
    -ffree-line-length-none
    $<$<CONFIG:Debug>:-DDEBG;-O0;-g>
    $<$<CONFIG:Release>:-funroll-all-loops;-fno-f2c;-O3>)
  if(CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 10)
    target_compile_options(astr_build_options INTERFACE -fallow-argument-mismatch)
  endif()
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel|IntelLLVM")
  target_compile_options(astr_build_options INTERFACE
    -fpp
    -heaparrays
    -safe-cray-ptr
    -traceback
    $<$<CONFIG:Debug>:-DDEBG;-O0;-g>
    $<$<CONFIG:Release>:-O3>)
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "NAG")
  target_compile_options(astr_build_options INTERFACE
    -fpp
    $<$<CONFIG:Debug>:-O0;-g>
    $<$<CONFIG:Release>:-O3>)
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Cray")
  target_compile_options(astr_build_options INTERFACE
    -eF
    -N
    1023
    $<$<CONFIG:Debug>:-O0;-g>
    $<$<CONFIG:Release>:-O3>)
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "PGI|NVHPC")
  target_compile_options(astr_build_options INTERFACE
    -cpp
    -Mfree
    -Kieee
    $<$<CONFIG:Debug>:-DDEBG;-O0>
    $<$<CONFIG:Release>:-O3>)
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Fujitsu")
  target_compile_options(astr_build_options INTERFACE
    -Cpp
    $<$<CONFIG:Debug>:-O0>
    $<$<CONFIG:Release>:-O3>)
else()
  target_compile_options(astr_build_options INTERFACE
    $<$<CONFIG:Debug>:-O0;-g>
    $<$<CONFIG:Release>:-O2>)
endif()

set(HDF5_PREFER_PARALLEL TRUE)
find_package(HDF5 REQUIRED COMPONENTS Fortran HL)

find_package(MPI REQUIRED COMPONENTS Fortran)

add_library(astr_dependencies INTERFACE)

target_compile_definitions(astr_dependencies INTERFACE HDF5)
target_include_directories(astr_dependencies INTERFACE ${HDF5_INCLUDE_DIRS})
target_link_libraries(astr_dependencies INTERFACE
  ${HDF5_LIBRARIES}
  ${HDF5_HL_LIBRARIES}
  ${HDF5_Fortran_LIBRARIES}
  ${HDF5_Fortran_HL_LIBRARIES}
  MPI::MPI_Fortran)

if(CHEMISTRY)
  if(NOT CANTERA_DIR)
    message(FATAL_ERROR "CHEMISTRY=ON requires -DCANTERA_DIR=<path>")
  endif()
  if(NOT EXISTS "${CANTERA_DIR}")
    message(FATAL_ERROR "CANTERA_DIR does not exist: ${CANTERA_DIR}")
  endif()

  target_compile_definitions(astr_dependencies INTERFACE COMB)
  target_include_directories(astr_dependencies INTERFACE "${CANTERA_DIR}/include/cantera")
  target_link_directories(astr_dependencies INTERFACE "${CANTERA_DIR}/lib")
  target_link_libraries(astr_dependencies INTERFACE cantera_fortran cantera stdc++ pthread)
endif()

set(ASTR_ROOT "${CMAKE_CURRENT_SOURCE_DIR}")
set(test_dir "${PROJECT_BINARY_DIR}/Test")

add_subdirectory(src)
add_subdirectory(examples)

if(BUILD_TESTING)
  set(ASTR_MPI_OVERSUBSCRIBE_FLAG)
  if(ASTR_MPI_OVERSUBSCRIBE)
    set(ASTR_MPI_OVERSUBSCRIBE_FLAG --oversubscribe)
  endif()

  function(astr_add_mpi_test test_name ranks workdir input_file label)
    add_test(
      NAME "${test_name}"
      COMMAND
        ${MPIEXEC_EXECUTABLE}
        ${MPIEXEC_NUMPROC_FLAG} ${ranks}
        ${ASTR_MPI_OVERSUBSCRIBE_FLAG}
        $<TARGET_FILE:astr>
        run
        "${input_file}"
      WORKING_DIRECTORY "${workdir}")
    set_tests_properties("${test_name}" PROPERTIES LABELS "${label}")
  endfunction()

  astr_add_mpi_test(SOD 2 "${test_dir}/sod" "./datin/input.sod" nondim)
  astr_add_mpi_test(Shuosher 2 "${test_dir}/Shuosher" "./datin/input.shuosher" nondim)
  astr_add_mpi_test(VT 8 "${test_dir}/VT" "./datin/input.2dvort" nondim)
  astr_add_mpi_test(mixinglayer 8 "${test_dir}/mixinglayer" "./datin/input.2d" nondim)
  astr_add_mpi_test(WT2D 8 "${test_dir}/WT2D" "./datin/input2d.windtunn" nondim)
  astr_add_mpi_test(TGV 8 "${test_dir}/TGV" "./datin/input.tgv" nondim)
  astr_add_mpi_test(RTI 8 "${test_dir}/RTI" "./datin/input.rti" nondim)
  astr_add_mpi_test(Channel 8 "${test_dir}/Channel" "./datin/input.chl" nondim)
  astr_add_mpi_test(SWLBI 8 "${test_dir}/SWLBI" "./datin/input.2d" nondim)
  astr_add_mpi_test(onedflame 1 "${test_dir}/onedflame" "./datin/input.onedflame" combustion)
  astr_add_mpi_test(airreactor 1 "${test_dir}/airreactor" "./datin/input.airreactor" combustion)
  astr_add_mpi_test(TGVflame 8 "${test_dir}/TGVflame" "./datin/input.tgvflame" combustion)
endif()
