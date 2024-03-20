install(DIRECTORY DESTINATION ${CMAKE_INSTALL_PREFIX}/examples/onedflame)
install(DIRECTORY datin DESTINATION ${CMAKE_INSTALL_PREFIX}/examples/onedflame)

# If testing active add test for Shuosher case
if (${BUILD_TESTING})
  set(case_dir "${test_dir}/onedflame")
  file(MAKE_DIRECTORY ${case_dir})
  file(COPY datin DESTINATION ${case_dir})
endif()
