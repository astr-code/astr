cmake_minimum_required(VERSION 3.12)
project(PASTR LANGUAGES Fortran)

#-----------------------------------------
#  Build type
#-----------------------------------------
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING
        "Choose the type of build: Debug or Release" FORCE)
endif()

#-----------------------------------------
#  Fortran compiler flags
#-----------------------------------------
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(CMAKE_Fortran_FLAGS "-cpp -ffree-line-length-none -fbacktrace -g")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -funroll-loops")
    set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -DDEBUG")
endif()

#-----------------------------------------
#  MPI support
#-----------------------------------------
find_package(MPI REQUIRED)
if(MPI_FOUND)
    message(STATUS "MPI found: ${MPI_Fortran_COMPILER}")
    enable_language(Fortran)
    set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
endif()

#-----------------------------------------
#  OpenMP support
#-----------------------------------------
find_package(OpenMP)
if(OpenMP_Fortran_FOUND)
    message(STATUS "OpenMP found")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
endif()

#-----------------------------------------
#  HDF5 support (with MPI)
#-----------------------------------------
set(HDF5_PREFER_PARALLEL TRUE)
find_package(HDF5 REQUIRED COMPONENTS Fortran HL)

message(STATUS "HDF5_Fortran_COMPILER found: ${HDF5_Fortran_COMPILER_EXECUTABLE}")
message(STATUS "HDF5_DIR: ${HDF5_DIR}")

# Warning if Include are not found => can be fixed with more recent cmake version
if (HDF5_FOUND)
    message(STATUS "HDF5 FOUND: ${HDF5_FOUND}")
    include_directories(${HDF5_INCLUDE_DIRS})
    message(STATUS "HDF5 INCL: ${HDF5_INCLUDE_DIRS}")
    link_libraries(${HDF5_LIBRARIES})
    message(STATUS "HDF5_LIB: ${HDF5_LIBRARIES}")
    link_libraries(${HDF5_HL_LIBRARIES})
    message(STATUS "HDF5_HL_LIB: ${HDF5_HL_LIBRARIES}")
    link_libraries(${HDF5_Fortran_LIBRARIES})
    message(STATUS "HDF5_Fortran_LIBRARIES: ${HDF5_Fortran_LIBRARIES}")
    link_libraries(${HDF5_Fortran_HL_LIBRARIES})
    message(STATUS "HDF5_Fortran_HL_LIBRARIES: ${HDF5_Fortran_HL_LIBRARIES}")
    add_definitions("-DHDF5")
else (HDF5_FOUND)
    message(STATUS "NO HDF5 include have been found. The executable won't be targeted with HDF5 include")
    message(STATUS "Code will compile but performaces can be compromised")
    message(STATUS "Using a CMake vers > 3.10 should solve the problem")
    message(STATUS "Alternatively use ccmake to manually set the include if available")
endif (HDF5_FOUND)

#-----------------------------------------
#  Output directories
#-----------------------------------------
include(GNUInstallDirs)
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR})

#-----------------------------------------
#  Add source
#-----------------------------------------
add_subdirectory(src)

#-----------------------------------------
#  Installation
#-----------------------------------------
install(TARGETS pastr
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
