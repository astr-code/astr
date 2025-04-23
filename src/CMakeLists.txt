message(STATUS "CMAKE_SOURCE_DIR: ${CMAKE_SOURCE_DIR}")
message(STATUS "PROJECT_BINARY_DIR: ${PROJECT_BINARY_DIR}")

# message(STATUS "the local path: ${PROJECT_BINARY_DIR}")

set(UDF_SRC ${CMAKE_SOURCE_DIR}/user_define_module)
file(GLOB UDF_SOURCE_FILES ${UDF_SRC}/*.F90)

message(STATUS "UDF path: ${UDF_SRC}")
    
add_executable(astr 
                    astr.F90
                    ${UDF_SOURCE_FILES}
                    bc.F90
                    cmdefne.F90
                    commarray.F90
                    comsolver.F90
                    commcal.F90
                    commfunc.F90
                    commtype.F90
                    commvar.F90
                    constdef.F90
                    fdnn.F90
                    fludyna.F90
                    geom.F90
                    gridgeneration.F90
                    hdf5io.F90
                    ibmethod.F90
                    initialisation.F90
                    interp.F90
                    io.F90
                    mainloop.F90
                    methodmoment.F90
                    models.F90
                    parallel.F90
                    pp.F90
                    readwrite.F90
                    riemann.F90
                    singleton.F90
                    solver.F90
                    statistic.F90
                    stlaio.F90
                    strings.F90
                    tecio.F90
                    test.F90
                    thermchem.F90
                    utility.F90
                    vtkio.F90)

target_link_libraries(astr)
if (MPI_FOUND)
  target_link_libraries(astr PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)
if (CHEMISTRY)
  message(STATUS "CANTERA library: ${CTRDIR}/lib")
  target_link_libraries(astr PRIVATE -L${CTRDIR}/lib -lcantera_fortran -lcantera -lstdc++ -pthread)
endif (CHEMISTRY)

install(TARGETS astr
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
