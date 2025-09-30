install(DIRECTORY DESTINATION ${CMAKE_INSTALL_PREFIX}/examples/air_reactor)
install(DIRECTORY datin DESTINATION ${CMAKE_INSTALL_PREFIX}/examples/air_reactor)

# If testing active add test for air_reactor case
if (${BUILD_TESTING})
  set(case_dir "${test_dir}/airreactor")
  file(MAKE_DIRECTORY ${case_dir})
  file(COPY datin DESTINATION ${case_dir})
endif()
