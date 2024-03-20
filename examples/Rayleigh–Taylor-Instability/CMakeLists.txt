install(DIRECTORY DESTINATION ${CMAKE_INSTALL_PREFIX}/examples/Rayleigh–Taylor-Instability)
install(DIRECTORY datin DESTINATION ${CMAKE_INSTALL_PREFIX}/examples/Rayleigh–Taylor-Instability)

# If testing active add test for RTI case
if (${BUILD_TESTING})
  set(case_dir "${test_dir}/RTI")
  file(MAKE_DIRECTORY ${case_dir})
  file(COPY datin DESTINATION ${case_dir})
endif()
