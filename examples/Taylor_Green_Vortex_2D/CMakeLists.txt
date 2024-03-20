install(DIRECTORY DESTINATION ${CMAKE_INSTALL_PREFIX}/examples/Taylor_Green_Vortex)
install(DIRECTORY datin DESTINATION ${CMAKE_INSTALL_PREFIX}/examples/Taylor_Green_Vortex)

# If testing active add test for TGV case
if (${BUILD_TESTING})
  set(case_dir "${test_dir}/TGV2D")
  file(MAKE_DIRECTORY ${case_dir})
  file(COPY datin DESTINATION ${case_dir})
endif()
