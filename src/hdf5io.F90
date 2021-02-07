!+---------------------------------------------------------------------+
!| This module contains HDF5 IO subroutines.                           |
!| ==============                                                      |
!| CHANGE RECORD                                                       |
!| -------------                                                       |
!| 02-Jun-2020  | Created by J. Fang STFC Daresbury Laboratory         |
!+---------------------------------------------------------------------+
module hdf5io
  !
#ifdef HDF5
  use hdf5
  use h5lt
#endif
  !
  use parallel, only : lio,mpirank,ig0,jg0,kg0
  !
  implicit none
  !
  Interface h5read
    !
    module procedure h5ra3d_r8
    !
  end Interface h5read
  !
#ifdef HDF5
  integer(hid_t) :: h5file_id
#endif
  !+------------------------------------------+
  ! h5file_id : file id to operate h5 file io |
  !+------------------------------------------+
  !
  contains
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to open the h5file interface and assign   |
  !| h5file_id. For write each new file, this will be called first, but|
  !| once it is called, the file will be overwriten.                   |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 03-Jun-2020 | Created by J. Fang STFC Daresbury Laboratory        |
  !+-------------------------------------------------------------------+
  subroutine h5io_init(filename,mode)
    !
    use mpi, only: mpi_comm_world,mpi_info_null
    use parallel, only: mpirank
    !
    ! arguments
    character(len=*),intent(in) :: filename
    character(len=*),intent(in) :: mode
    ! h5file_id is returned
    !
#ifdef HDF5
    ! local data
    integer :: h5error
    integer(hid_t) :: plist_id
    !
    call h5open_f(h5error)
    if(h5error.ne.0)  stop ' !! error in h5io_init call h5open_f'
    !
    ! create access property list and set mpi i/o
    call h5pcreate_f(h5p_file_access_f,plist_id,h5error)
    if(h5error.ne.0)  stop ' !! error in h5io_init call h5pcreate_f'
    !
    call h5pset_fapl_mpio_f(plist_id,mpi_comm_world,mpi_info_null,     &
                                                                h5error)
    if(h5error.ne.0)  stop ' !! error in h5io_init call h5pset_fapl_mpio_f'
    !
    if(mode=='writ') then
      call h5fcreate_f(filename,h5f_acc_trunc_f,h5file_id,             &
                                            h5error,access_prp=plist_id)
      if(h5error.ne.0)  stop ' !! error in h5io_init call h5fcreate_f'
    elseif(mode=='read') then
      call h5fopen_f(filename,h5f_acc_rdwr_f,h5file_id,                &
                                            h5error,access_prp=plist_id)
      if(h5error.ne.0)  stop ' !! error in h5io_init call h5fopen_f'
    else
        stop ' !! mode not defined @ h5io_init'
    endif
    !
    call h5pclose_f(plist_id,h5error)
    if(h5error.ne.0)  stop ' !! error in h5io_init call h5pclose_f'
    !
    if(mpirank==0) print*,' ** open h5 file: ',filename
    !
#endif
    !
  end subroutine h5io_init
  !+-------------------------------------------------------------------+
  !| This end of the subroutine h5io_init.                             |
  !+-------------------------------------------------------------------+
  !!
  !+-------------------------------------------------------------------+
  !| This subroutine is used to close hdf5 interface after finish      |
  !| input/output a hdf5 file.                                         |
  !| the only data needed is h5file_id                                 |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 03-Jun-2020 | Created by J. Fang STFC Daresbury Laboratory        |
  !+-------------------------------------------------------------------+
  subroutine h5io_end
    !
    ! local data
    integer :: h5error
    !
#ifdef HDF5
    call h5fclose_f(h5file_id,h5error)
    if(h5error.ne.0)  stop ' !! error in h5io_end call h5fclose_f'
    !
    call h5close_f(h5error)
    if(h5error.ne.0)  stop ' !! error in h5io_end call h5close_f'
    !
#endif
    !
  end subroutine h5io_end
  !+-------------------------------------------------------------------+
  !| This end of the subroutine h5io_end.                              |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to read 3D real8 araay using hdf5         |
  !| the only data needed is h5file_id                                 |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 03-Jun-2020 | Created by J. Fang STFC Daresbury Laboratory        |
  !+-------------------------------------------------------------------+
  subroutine h5ra3d_r8(varname,var)
    !
    use mpi, only: mpi_comm_world,mpi_info_null
    !
    ! arguments
    character(LEN=*),intent(in) :: varname
    real(8),intent(inout) :: var(:,:,:)
    !
#ifdef HDF5
    ! local data
    integer :: jrk
    integer :: dim(3)
    integer(hsize_t), dimension(3) :: offset
    integer :: h5error
    !
    integer(hid_t) :: dset_id,filespace,memspace,plist_id
    integer(hsize_t) :: dimt(3)
    !
    dim(1)=size(var,1)
    dim(2)=size(var,2)
    dim(3)=size(var,3)
    !
    dimt=dim
    offset=(/ig0,jg0,kg0/)
    !
    ! read the data
    !
    call h5dopen_f(h5file_id,varname,dset_id,h5error)
    call h5screate_simple_f(3,dimt,memspace,h5error)
    call h5dget_space_f(dset_id,filespace,h5error)
    call h5sselect_hyperslab_f(filespace,h5s_select_set_f,offset,      &
                                                           dimt,h5error)
    ! Create property list for collective dataset read
    call h5pcreate_f(h5p_dataset_xfer_f,plist_id,h5error)
    call h5pset_dxpl_mpio_f(plist_id,h5fd_mpio_collective_f,h5error)
    !
    ! Read 3D array
    call h5dread_f(dset_id,h5t_native_double,var,dimt,h5error,         &
                   mem_space_id=memspace,file_space_id=filespace,      &
                                                      xfer_prp=plist_id)
    !Close dataspace
    call h5sclose_f(filespace,h5error)
    !
    call h5sclose_f(memspace,h5error)
    !
    call h5dclose_f(dset_id,h5error)
    !
    call h5pclose_f(plist_id,h5error)
    !
    if(lio) print*,' >> ',varname
    !
#endif
    !
  end subroutine h5ra3d_r8
  !+-------------------------------------------------------------------+
  !| This end of the function h5ra3d_r8.                               |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This function is used to get the dimension of the hdf5 array.     |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 30-JuL-2020 | Coped from ASTR Post by J. Fang STFC Daresbury Lab. |
  !+-------------------------------------------------------------------+
  function h5getdim3d(varname,filenma) result(dims)
    !
    character(len=*),intent(in) :: varname,filenma
    integer :: dims(3)
    !
#ifdef HDF5
    ! local data
    integer(hid_t)  :: file, space, dset
    integer(hsize_t) :: ndims(3)
    integer         :: h5error ! error flag
    integer(hid_t) :: dims_h5(3)
    !
    !
    call h5open_f(h5error)
    call h5fopen_f(filenma,h5f_acc_rdonly_f,file,h5error)
    call h5dopen_f (file,varname,dset, h5error)
    call h5dget_space_f(dset, space, h5error)
    call h5sget_simple_extent_dims_f(space,dims_h5,ndims,h5error)
    !
    dims=dims_h5
    !
    call h5dclose_f(dset , h5error)
    call h5sclose_f(space, h5error)
    call h5fclose_f(file , h5error)
    call h5close_f(h5error)
    !
#endif
    !
  end function h5getdim3d
  !+-------------------------------------------------------------------+
  !| This end of the function h5getdim3d.                              |
  !+-------------------------------------------------------------------+
  !
end module hdf5io