message(STATUS "CMAKE_SOURCE_DIR: ${CMAKE_SOURCE_DIR}")
message(STATUS "PROJECT_BINARY_DIR: ${PROJECT_BINARY_DIR}")

add_executable(astr_omp
                    constdef.f90
                    xdmwrite.f90
                    tecio.f90
                    comvardef.f90
                    utility.f90
                    io.f90
                    fluidynamcs.f90
                    field_init.f90
                    bc.f90
                    numerics.f90
                    solver.f90
                    statistics.f90
                    tgvsolver_omp.f90)
target_link_libraries(astr_omp)

install(TARGETS astr_omp
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
