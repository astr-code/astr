message(STATUS "CMAKE_SOURCE_DIR: ${CMAKE_SOURCE_DIR}")
message(STATUS "PROJECT_BINARY_DIR: ${PROJECT_BINARY_DIR}")

# Automatically gather all Fortran sources under src/
file(GLOB PASTR_SOURCES
    "*.F90"
)

# Ensure modules compile first if needed
list(SORT PASTR_SOURCES)  # optional: ensures *.mod files first

# Create the executable
add_executable(pastr ${PASTR_SOURCES})

# Include directories
target_include_directories(pastr PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

# Link optional libraries
if (MPI_FOUND)
  target_link_libraries(pastr PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)
if (PASTR_OpenMP_LIBS)
    target_link_libraries(pastr PRIVATE ${PASTR_OpenMP_LIBS})
endif()
# Installation
install(TARGETS pastr
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
