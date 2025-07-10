install(DIRECTORY DESTINATION ${CMAKE_INSTALL_PREFIX}/examples/sod)
install(DIRECTORY datin DESTINATION ${CMAKE_INSTALL_PREFIX}/examples/sod)

# If testing active add test for sod case
if (${BUILD_TESTING})
  set(case_dir "${test_dir}/sod")
  file(MAKE_DIRECTORY ${case_dir})
  file(COPY datin DESTINATION ${case_dir})
endif()
