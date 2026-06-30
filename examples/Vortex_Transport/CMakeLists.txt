install(DIRECTORY DESTINATION examples/Vortex_Transport)
install(DIRECTORY datin DESTINATION examples/Vortex_Transport)

# If testing active add test for vortex trasport case
if (${BUILD_TESTING})
  set(case_dir "${test_dir}/VT")
  file(MAKE_DIRECTORY ${case_dir})
  file(COPY datin DESTINATION ${case_dir})
endif()
