install(DIRECTORY DESTINATION ${CMAKE_INSTALL_PREFIX}/examples/Riemann2D)
install(DIRECTORY datin DESTINATION ${CMAKE_INSTALL_PREFIX}/examples/Riemann2D)

# If testing active add test for Riemann2D case
if (${BUILD_TESTING})
  set(case_dir "${test_dir}/Riemann2D")
  file(MAKE_DIRECTORY ${case_dir})
  file(COPY datin DESTINATION ${case_dir})
endif()
