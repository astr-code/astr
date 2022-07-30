install(DIRECTORY DESTINATION ${CMAKE_INSTALL_PREFIX}/examples/Wind_Tunnel_2D)
install(DIRECTORY datin DESTINATION ${CMAKE_INSTALL_PREFIX}/examples/Wind_Tunnel_2D)

# If testing active add test for Shuosher case
if (${BUILD_TESTING})
  set(case_dir "${test_dir}/WT2D")
  file(MAKE_DIRECTORY ${case_dir})
  file(COPY datin DESTINATION ${case_dir})
endif()
