message(STATUS "CMAKE_SOURCE_DIR: ${CMAKE_SOURCE_DIR}")
message(STATUS "PROJECT_BINARY_DIR: ${PROJECT_BINARY_DIR}")

# Set the root directory for the source files
set(SRC_DIR ${CMAKE_SOURCE_DIR}/src)

# Collect all Fortran source files recursively
file(GLOB_RECURSE Fortran_SOURCES
    ${SRC_DIR}/*.f90
    ${SRC_DIR}/*.F90
    ${SRC_DIR}/*.f
    ${SRC_DIR}/*.F 
)

add_executable(astr ${Fortran_SOURCES})

target_link_libraries(astr)
if (MPI_FOUND)
  target_link_libraries(astr PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

if (CHEMISTRY)
  message(STATUS "CANTERA library: ${CTRDIR}/lib")
  target_link_libraries(astr PRIVATE -L${CTRDIR}/lib -lcantera_fortran -lcantera -lstdc++ -pthread)
endif (CHEMISTRY)

install(TARGETS astr
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
