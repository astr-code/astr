install(DIRECTORY DESTINATION ${CMAKE_INSTALL_PREFIX}/examples/Lid-Driven-Cavity-Rarefied)
install(DIRECTORY datin DESTINATION ${CMAKE_INSTALL_PREFIX}/examples/Lid-Driven-Cavity-Rarefied)

# If testing active add test for cavity flow involving moment equations
if (${BUILD_TESTING})
  set(case_dir "${test_dir}/CavMoM")
  file(MAKE_DIRECTORY ${case_dir})
  file(COPY datin DESTINATION ${case_dir})
endif()