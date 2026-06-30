install(DIRECTORY DESTINATION examples/Taylor_Green_Vortex)
install(DIRECTORY datin DESTINATION examples/Taylor_Green_Vortex)

# If testing active add test for TGV case
if (${BUILD_TESTING})
  set(case_dir "${test_dir}/TGV")
  file(MAKE_DIRECTORY ${case_dir})
  file(COPY datin DESTINATION ${case_dir})
endif()
