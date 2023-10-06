install(DIRECTORY DESTINATION ${CMAKE_INSTALL_PREFIX}/examples/MixingLayer)
install(DIRECTORY datin DESTINATION ${CMAKE_INSTALL_PREFIX}/examples/MixingLayer)

# If testing active add test for mixing layer case
if (${BUILD_TESTING})
  set(case_dir "${test_dir}/mixinglayer")
  file(MAKE_DIRECTORY ${case_dir})
  file(COPY datin DESTINATION ${case_dir})
endif()