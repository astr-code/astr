install(DIRECTORY DESTINATION ${CMAKE_INSTALL_PREFIX}/examples/Channel)
install(DIRECTORY datin DESTINATION ${CMAKE_INSTALL_PREFIX}/examples/Channel)

# If testing active add test for Channel case
if (${BUILD_TESTING})
  set(case_dir "${test_dir}/Channel")
  file(MAKE_DIRECTORY ${case_dir})
  file(COPY datin DESTINATION ${case_dir})
endif()
