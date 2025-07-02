cmake_minimum_required(VERSION 3.10)
project(astr_min Fortran)

# Set the Fortran standard
set(CMAKE_Fortran_STANDARD 90)

# Add the source files
set(SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/../src/derivative.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/../src/filter.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/../src/commfunc.F90 
    ${CMAKE_CURRENT_SOURCE_DIR}/../src/commtype.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/../src/commvar.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/../src/constdef.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/../src/singleton.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/tgv_solver_3d/tgvsolver.F90

)

# Define the executable
add_executable(astr.min ${SOURCES})