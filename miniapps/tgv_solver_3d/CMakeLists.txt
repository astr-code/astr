message(STATUS "CMAKE_SOURCE_DIR: ${CMAKE_SOURCE_DIR}")
message(STATUS "PROJECT_BINARY_DIR: ${PROJECT_BINARY_DIR}")

add_executable(astr_omp
                    com/constdef.f90
                    com/xdmwrite.f90
                    com/tecio.f90
                    com/comvardef.f90
                    com/utility.f90
                    com/io.f90
                    omp/fluidynamcs.f90
                    omp/field_init.f90
                    omp/bc.f90
                    omp/numerics.f90
                    omp/solver.f90
                    omp/statistics.f90
                    omp/tgvsolver_omp.f90)


# add_executable(astr_oacc
#                     com/constdef.f90
#                     com/xdmwrite.f90
#                     com/tecio.f90
#                     com/comvardef.f90
#                     com/utility.f90
#                     com/io.f90
#                     oacc/fluidynamcs.f90
#                     oacc/field_init.f90
#                     oacc/bc.f90
#                     oacc/numerics.f90
#                     oacc/solver.f90
#                     oacc/statistics.f90
#                     oacc/tgvsolver_oacc.f90)
# 
# 
add_executable(astr_cpu
                    com/constdef.f90
                    com/xdmwrite.f90
                    com/tecio.f90
                    com/comvardef.f90
                    com/utility.f90
                    com/io.f90
                    omp/fluidynamcs.f90
                    omp/field_init.f90
                    omp/bc.f90
                    omp/numerics.f90
                    omp/solver.f90
                    omp/statistics.f90
                    omp/tgvsolver_omp.f90)

target_link_libraries(astr_cpu)

find_package(OpenMP REQUIRED)
target_link_libraries(astr_omp PRIVATE OpenMP::OpenMP_Fortran)

# find_package(OpenACC REQUIRED)
# target_link_libraries(astr_oacc PRIVATE OpenACC::OpenACC_Fortran)