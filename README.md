# ASTR code 
Version 2.0 

ASTR code is a high-order finite-difference flow solver for compressible turbulence research.

# Download, Installation and Compilation
Required dependencies: Fortran 90, MPI, HDF5

## Download the astr code:
#git clone git@gitlab.com:astr-code/astr.git

git clone git@github.com:astr-code/astr.git

## Compilation and installation:
ASTR 2.0 supports both gnu make and cmake.
For the use of gnu make, do:
make 
in the directory containing src folder, and the executable will be found as ./bin/astr 

The cmake gives a more complete and safe way of compiling and installing the code.
By default ASTR solves equations under non-dimensional form, unless the chemstry is included.

create a case folder, e.g.
mkdir test_case
cd test_case
cmake path_to_the_source
make
make install
ctest -L nondim

The code will be installed in test_case and excutable can be found at test_case/bin/astr

If you want to use the chemstry function, you need first to install cantera. After download and unpack the cantera, you can use the following script to install:
python scons/scripts/scons.py build python_package=none FORTRAN=<your fortran compiler> f90_interface=y prefix=<the directory of centera to install> boost_inc_dir=<to boost directory>

python scons/scripts/scons.py test

python scons/scripts/scons.py install

you may need to make and test ASTR with chemstry with following cmake commands:
export CANTERA_DIR=<the directory of centera to install>
cmake
make 
make install
ctest -L combustion


## Run the solver:

Once the excutable is built, a typical simulation can be run as,
mpirun -np 8 ./astr run ./datin/input_file






