download cantera-2.5.1.tgz
tar -xvf cantera-2.5.1.tgz
cd cantera-2.5.1
module load python3
python3 scons/scripts/scons.py build python_package=none FORTRAN=mpif90 f90_interface=y prefix=/home/abr01399/opt/cantera-2.5.1/ boost_inc_dir=/usr/include/boost/
python scons/scripts/scons.py install
