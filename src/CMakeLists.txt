add_executable(astr
                    astr.F90
                    bc.F90
                    cmdefne.F90
                    commarray.F90
                    comsolver.F90
                    commcal.F90
                    commfunc.F90
                    commtype.F90
                    commvar.F90
                    constdef.F90
                    derivative.F90
                    fdnn.F90
                    filter.F90
                    fludyna.F90
                    flux.F90
                    geom.F90
                    gridgeneration.F90
                    hdf5io.F90
                    ibmethod.F90
                    initialisation.F90
                    interp.F90
                    mainloop.F90
                    models.F90
                    mpiio.F90
                    parallel.F90
                    pp.F90
                    readwrite.F90
                    riemann.F90
                    singleton.F90
                    solver.F90
                    sponge_layer.F90
                    statistic.F90
                    stlaio.F90
                    strings.F90
                    tecio.F90
                    test.F90
                    thermchem.F90
                    ../user_define_module/userdefine.F90
                    utility.F90
                    vtkio.F90)

target_link_libraries(astr
  PRIVATE
    astr_build_options
    astr_dependencies)

install(TARGETS astr
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
