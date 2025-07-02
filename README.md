ASTR Code

Version 2.5

ASTR is a high-order finite-difference flow solver designed for high-fidelity simulation of compressible turbulence. It supports multi-physics extensions including combustion and is optimized for modern high-performance computing systems.

Dependencies
Before building ASTR, ensure the following dependencies are installed:

Fortran 90 compiler (e.g., gfortran, ifort)
MPI (e.g., OpenMPI, MPICH)
HDF5 with Fortran bindings
(Optional) Cantera (for combustion simulations)

Download the Source Code
Clone the official repository:
git clone git@github.com:astr-code/astr.git

Compilation and Installation

Option 1: Using make
A simple build using GNU Make:

cd astr
make
The executable will be located at:
./bin/astr

Option 2: Using CMake (Recommended)
CMake provides a safer and more flexible build environment. To compile and install using CMake:

Create a case directory:

mkdir test_case

cmake path_to_the_source
cmake --build .
cmake --install .
ctest -L nondim
The binary will be installed under:

test_case/bin/astr

Enabling Combustion Module
ASTR supports detailed chemical kinetics via Cantera. To enable this feature:

Install Cantera (Fortran interface required):

python scons/scripts/scons.py build python_package=none \
    FORTRAN=<your fortran compiler> f90_interface=y \
    prefix=<installation_dir> boost_inc_dir=<boost_include_dir>

python scons/scripts/scons.py test
python scons/scripts/scons.py install

Configure ASTR with Cantera support:
cmake -DCHEMISTRY=TRUE -DCANTERA_DIR=path_to_cantera path_to_the_source
cmake --build .
cmake --install .
ctest -L combustion

Running a Simulation
To execute a simulation:
mpirun -np 8 ./astr run ./datin/input_file

Mini Apps
ASTR includes lightweight mini-applications for testing and development.

To compile and run mini-apps:

mkdir test_mini_apps
cd test_mini_apps
cmake path_to_the_source/miniapps/
cmake --build .
./astr.min

Directory Structure Overview

astr/
├── src/                # Core solver source code
├── script/             # Utility scripts (e.g. case creation)
├── miniapps/           # Testing and development tools
├── bin/                # Compiled binaries
└── examples/           # Sample cases

Contact
For questions, bug reports, or contributions, please open an issue on GitHub or contact the development team.

Let me know if you want to include badges (e.g., build status), citation info, or extended documentation links.