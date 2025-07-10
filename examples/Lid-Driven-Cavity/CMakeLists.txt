install(DIRECTORY DESTINATION ${CMAKE_INSTALL_PREFIX}/examples/Lid-Driven-Cavity)
install(DIRECTORY datin DESTINATION ${CMAKE_INSTALL_PREFIX}/examples/Lid-Driven-Cavity)

# If testing active add test for Lid-Driven-Cavity case
if (${BUILD_TESTING})
  set(case_dir "${test_dir}/Lid-Driven-Cavity")
  file(MAKE_DIRECTORY ${case_dir})
  file(COPY datin DESTINATION ${case_dir})
endif()
