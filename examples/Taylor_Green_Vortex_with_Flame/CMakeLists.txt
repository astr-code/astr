install(DIRECTORY DESTINATION ${CMAKE_INSTALL_PREFIX}/examples/Taylor_Green_Vortex_with_Flame)
install(DIRECTORY datin DESTINATION ${CMAKE_INSTALL_PREFIX}/examples/Taylor_Green_Vortex_with_Flame)

# If testing active add test for TGV case
if (${BUILD_TESTING})
  set(case_dir "${test_dir}/TGVflame")
  file(MAKE_DIRECTORY ${case_dir})
  file(COPY datin DESTINATION ${case_dir})
endif()
