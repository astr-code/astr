wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.12/hdf5-1.12.0/src/hdf5-1.12.0.tar.bz2
tar -xvf hdf5-1.12.0.tar.bz2
cd hdf5-1.12.0
./configure --prefix=$HD5_INSTALL_DIR --enable-fortran --disable-shared --enable-parallel CC=mpicc FC=mpif90
make 
make install

echo 'export HDF5_INSTALL_PATH=$HD5_INSTALL_DIR' >> .bashrc
echo 'export PATH=$PATH:$HD5_INSTALL_DIR/bin/' >> .bashrc