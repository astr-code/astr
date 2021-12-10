install(DIRECTORY DESTINATION ${CMAKE_INSTALL_PREFIX}/examples/Shuosher)
install(DIRECTORY datin DESTINATION ${CMAKE_INSTALL_PREFIX}/examples/Shuosher)

# If testing active add test for Shuosher case
if (${BUILD_TESTING})
  set(case_dir "${test_dir}/Shuosher")
  file(MAKE_DIRECTORY ${case_dir})
  file(COPY datin DESTINATION ${case_dir})
endif()
