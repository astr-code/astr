install(DIRECTORY DESTINATION ${CMAKE_INSTALL_PREFIX}/examples/Taylor_Green_Vortex_MoM)
install(DIRECTORY datin DESTINATION ${CMAKE_INSTALL_PREFIX}/examples/Taylor_Green_Vortex_MoM)

# If testing active add test for TGV case involving moment equations
if (${BUILD_TESTING})
  set(case_dir "${test_dir}/TGVMoM")
  file(MAKE_DIRECTORY ${case_dir})
  file(COPY datin DESTINATION ${case_dir})
endif()
