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
  use parallel, only : lio,irk,jrk,krk,mpirank,ig0,jg0,kg0,            &
                       mpi_comm_world,mpi_info_null,mpistop
  !
  implicit none
  !
  Interface h5read
    !
    module procedure h5ra2d_r8
    module procedure h5ra3d_r8
    module procedure h5ra3d_i4
    module procedure h5r_int4
    module procedure h5r_real8
    module procedure h5r_real8_1d
    !
  end Interface h5read
  !
  Interface h5write
    !
    module procedure h5wa2d_r8
    module procedure h5wa3d_r8
    module procedure h5w_int4
    module procedure h5w_real8
    !
  end Interface h5write
  !
  Interface h5sread
    !
    module procedure h5_readarray1d
    module procedure h5_read1rl8
    module procedure h5_read1int
    module procedure h5_readarray3d
    !
  end Interface h5sread
  !
  Interface h5srite
    !
    module procedure h5_write1int
    module procedure h5_write1rl8
    module procedure h5_writearray1dint
    module procedure h5_writearray1d
    module procedure h5_writearray3d
    !
  end Interface h5srite
  !
  Interface H5ReadSubset
    !
    module procedure h5_read2dfrom3d
    module procedure h5_read1dfrom3d
    module procedure h5_read1dfrom2d
    !
  end Interface H5ReadSubset
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
  subroutine h5io_init(filename,mode,comm)
    !
    use parallel, only: mpirank
    !
    ! arguments
    character(len=*),intent(in) :: filename
    character(len=*),intent(in) :: mode
    integer,intent(in),optional :: comm
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
    if(present(comm)) then
      call h5pset_fapl_mpio_f(plist_id,comm,mpi_info_null,h5error)
    else
      call h5pset_fapl_mpio_f(plist_id,mpi_comm_world,mpi_info_null,   &
                                                                h5error)
    endif
    if(h5error.ne.0)  stop ' !! error in h5io_init call h5pset_fapl_mpio_f'
    !
    if(mode=='write') then
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
  !| This subroutine is a interface of read 3D data array using hdf5   |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 05-04-2022 | Created by J. Fang STFC Daresbury Laboratory         |
  !+-------------------------------------------------------------------+
  subroutine h5ra3d_r8(varname,var,mode)
    !
    ! arguments
    character(LEN=*),intent(in) :: varname
    real(8),intent(inout) :: var(:,:,:)
    character(len=1),intent(in) :: mode
    !
    if(mode=='h') then
      call h5ra3d_r8_struct(varname,var)
    elseif(mode=='s') then
      call h5ra3d_r8_stream(varname,var)
    else
      stop ' !! error @  h5ra3d_r8'
    endif
    !
  end subroutine h5ra3d_r8
  !+-------------------------------------------------------------------+
  !| This end of the subroutine h5ra3d_r8.                             |
  !+-------------------------------------------------------------------+
  !!
  !+-------------------------------------------------------------------+
  !| This subroutine is used to read 3D real8 araay using hdf5         |
  !| the only data needed is h5file_id                                 |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 03-Jun-2020 | Created by J. Fang STFC Daresbury Laboratory        |
  !| 30-Jun-2021 | Add h5ra3d_i4 to read 3-d integer array by J. Fang  |
  !+-------------------------------------------------------------------+
  subroutine h5ra3d_r8_struct(varname,var)
    !
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
  end subroutine h5ra3d_r8_struct
  !
  subroutine h5ra3d_r8_stream(varname,var)
    !
    use parallel,only: lio,mpirankmax,ptabupd
    !
    ! arguments
    character(LEN=*),intent(in) :: varname
    real(8),intent(inout) :: var(:,:,:)
    !
#ifdef HDF5
    ! local data
    real(8),allocatable :: var_1d(:)
    integer :: dima,jrank,i,j,k,n
    integer :: data_size(0:mpirankmax)
    integer :: h5error
    integer,save :: dim(1)
    !
    integer(hid_t) :: dset_id,filespace,memspace,plist_id
    integer(hsize_t),dimension(1),save :: offset,dimt,dimat
    integer,save :: size_sav(3)=0
    !
    ! read the data
    if( size_sav(1)==size(var,1) .and. size_sav(2)==size(var,2) .and. & 
        size_sav(3)==size(var,3) ) then
      continue
    else
      size_sav(1)=size(var,1)
      size_sav(2)=size(var,2)
      size_sav(3)=size(var,3)
      !
      dim(1)=size_sav(1)*size_sav(2)*size_sav(3)
      !
      call ptabupd(var=dim(1),table=data_size)
      !
      dima=0
      do jrank=0,mpirankmax
        dima=dima+data_size(jrank)
      enddo
      !
      offset(1)=0
      do jrank=0,mpirank-1
        offset(1)=offset(1)+data_size(jrank)
      enddo
      !
      dimat=dima
      !
      dimt=dim
    endif
    !
    ! print*,mpirank,'|',offset(1),offset(1)+dimt
    !
    allocate(var_1d(1:dim(1)))
    !
    call h5dopen_f(h5file_id,varname,dset_id,h5error)
    call h5screate_simple_f(1,dimt,memspace,h5error)
    call h5dget_space_f(dset_id,filespace,h5error)
    call h5sselect_hyperslab_f(filespace,h5s_select_set_f,offset,      &
                                                           dimt,h5error)
    ! Create property list for collective dataset read
    call h5pcreate_f(h5p_dataset_xfer_f,plist_id,h5error)
    call h5pset_dxpl_mpio_f(plist_id,h5fd_mpio_collective_f,h5error)
    !
    ! Read 1-d array
    call h5dread_f(dset_id,h5t_native_double,var_1d,dimt,h5error,      &
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
    n=0
    do k=1,size_sav(3)
    do j=1,size_sav(2)
    do i=1,size_sav(1)
      n=n+1
      var(i,j,k)=var_1d(n)
    enddo
    enddo
    enddo
    !
    deallocate(var_1d)
    !
    if(lio) print*,' >> ',varname
    !
#endif
    !
  end subroutine h5ra3d_r8_stream
  !
  subroutine h5ra3d_i4(varname,var)
    !
    !
    ! arguments
    character(LEN=*),intent(in) :: varname
    integer,intent(inout) :: var(:,:,:)
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
    call h5dread_f(dset_id,h5t_native_integer,var,dimt,h5error,        &
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
  end subroutine h5ra3d_i4
  !
  subroutine h5r_int4(varname,var)
    !
    ! arguments
    character(LEN=*),intent(in) :: varname
    integer,intent(out) :: var
    !
#ifdef HDF5
    ! local data
    integer :: nvar(1)
    integer :: h5error
    integer(hsize_t) :: dimt(1)=(/1/)
    !
    ! read the data
    !
    call h5ltread_dataset_f(h5file_id,varname,h5t_native_integer,nvar, &
                                                           dimt,h5error)
    var=nvar(1)
    if(h5error.ne.0)  stop ' !! error in h5r_int4 call h5ltread_dataset_f'
    !
    if(lio) print*,' >> ',varname,' ',var
    !
#endif
    !
  end subroutine h5r_int4
  !
  subroutine h5r_real8(varname,var)
    !
    ! arguments
    character(LEN=*),intent(in) :: varname
    real(8),intent(out) :: var
    !
#ifdef HDF5
    ! local data
    real(8) :: rvar(1)
    integer :: h5error
    integer(hsize_t) :: dimt(1)=(/1/)
    !
    ! read the data
    !
    call h5ltread_dataset_f(h5file_id,varname,h5t_native_double,rvar,  &
                                                           dimt,h5error)
    var=rvar(1)
    if(h5error.ne.0)  stop ' !! error in h5r_int4 call h5ltread_dataset_f'
    !
    if(lio) print*,' >> ',varname,' ',var
    !
#endif
    !
  end subroutine h5r_real8
  !
  subroutine h5r_real8_1d(varname,var,dim)
    !
    ! arguments
    character(LEN=*),intent(in) :: varname
    real(8),intent(out) :: var(:)
    integer,intent(in) :: dim
    !
#ifdef HDF5
    ! local data
    real(8) :: rvar(dim)
    integer :: h5error
    integer(hsize_t) :: dimt(1)
    !
    ! read the data
    dimt(1) =dim
    !
    call h5ltread_dataset_f(h5file_id,varname,h5t_native_double,rvar,  &
                                                           dimt,h5error)
    var=rvar
    if(h5error.ne.0)  stop ' !! error in h5r_int4 call h5ltread_dataset_f'
    !
    ! if(lio) print*,' >> ',varname,' ',var
    !
#endif
    !
  end subroutine h5r_real8_1d
  !
  subroutine h5ra2d_r8(varname,var,dir,display)
    !
    !
    ! arguments
    character(LEN=*),intent(in) :: varname
    real(8),intent(inout) :: var(:,:)
    character(len=1),intent(in) :: dir
    logical,intent(in),optional :: display
    !
#ifdef HDF5
    ! local data
    integer :: jrk
    integer :: dim(2)
    integer(hsize_t), dimension(2) :: offset
    integer :: h5error
    logical :: lexplicit
    !
    integer(hid_t) :: dset_id,filespace,memspace,plist_id
    integer(hsize_t) :: dimt(2)
    !
    if (present(display)) then
       lexplicit = display
    else
       lexplicit = .true.
    end if
    !
    dim(1)=size(var,1)
    dim(2)=size(var,2)
    !
    dimt=dim
    !
    if(dir=='i') then
      offset=(/jg0,kg0/)
    elseif(dir=='j') then
      offset=(/ig0,kg0/)
    elseif(dir=='k') then
      offset=(/ig0,jg0/)
    else
      stop ' !! error in dir @ h5wa2d_r8'
    endif
    !
    ! read the data
    !
    call h5dopen_f(h5file_id,varname,dset_id,h5error)
    call h5screate_simple_f(2,dimt,memspace,h5error)
    call h5dget_space_f(dset_id,filespace,h5error)
    call h5sselect_hyperslab_f(filespace,h5s_select_set_f,offset,      &
                                                           dimt,h5error)
    ! Create property list for collective dataset read
    call h5pcreate_f(h5p_dataset_xfer_f,plist_id,h5error)
    call h5pset_dxpl_mpio_f(plist_id,h5fd_mpio_collective_f,h5error)
    !
    ! Read 2D array
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
    if(lio .and. lexplicit) print*,' >> ',varname
    !
#endif
    !
  end subroutine h5ra2d_r8

  subroutine h5read_istrip_r8(varname,var)
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
    offset=(/0,jg0,kg0/)
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
  end subroutine h5read_istrip_r8
  !+-------------------------------------------------------------------+
  !| This end of the function h5read_istrip_r8.                        |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to write a 1D array with hdf5 interface.  |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 02-Jun-2020 | Created by J. Fang STFC Daresbury Laboratory        |
  !+-------------------------------------------------------------------+
  subroutine h5w_int4(varname,var)
    !
    ! arguments
    character(LEN=*),intent(in) :: varname
    integer,intent(in) :: var
    !
#ifdef HDF5
    ! local data
    integer :: nvar(1)
    integer :: h5error
    integer(hsize_t) :: dimt(1)=(/1/)
    !
    ! writing the data
    !
    nvar=var
    call h5ltmake_dataset_f(h5file_id,varname,1,dimt,                  &
                                        h5t_native_integer,nvar,h5error)
    if(h5error.ne.0)  stop ' !! error in h5w_int4 call h5ltmake_dataset_f'
    !
    if(lio) print*,' << ',varname
    !
#endif
    !
  end subroutine h5w_int4
  !
  subroutine h5w_real8(varname,var)
    !
    ! arguments
    character(LEN=*),intent(in) :: varname
    real(8),intent(in) :: var
    !
#ifdef HDF5
    ! local data
    real(8) :: rvar(1)
    integer :: h5error
    integer(hsize_t) :: dimt(1)=(/1/)
    !
    rvar=var
    call h5ltmake_dataset_f(h5file_id,varname,1,dimt,                  &
                                        h5t_native_double,rvar,h5error)
    if(h5error.ne.0)  stop ' !! error in h5w_real8 call h5ltmake_dataset_f'
    !
    if(mpirank==0) print*,' << ',varname
    !
#endif
    !
  end subroutine h5w_real8
  !
  subroutine h5wa3d_r8(varname,var,mode)
    !
    ! arguments
    character(LEN=*),intent(in) :: varname
    real(8),intent(in) :: var(:,:,:)
    character(len=1),intent(in) :: mode
    !
#ifdef HDF5
    if(mode=='h') then
      call h5wa3d_r8_struct(varname,var)
    elseif(mode=='s') then
      call h5wa3d_r8_stream(varname,var)
    else
      stop ' !! error @  h5wa3d_r8'
    endif
    !
#endif
    !
  end subroutine h5wa3d_r8
  !
  subroutine h5wa3d_r8_struct(varname,var)
    !
    use commvar, only: ia,ja,ka
    use parallel,only: lio,ig0,jg0,kg0
    !
    ! arguments
    character(LEN=*),intent(in) :: varname
    real(8),intent(in) :: var(:,:,:)
    !
#ifdef HDF5
    ! local data
    integer :: dim(3),dima
    integer(hsize_t), dimension(3) :: offset
    integer :: h5error
    !
    integer(hid_t) :: dset_id,filespace,memspace,plist_id
    integer(hsize_t) :: dimt(3),dimat(3)
    !
    dim(1)=size(var,1)
    dim(2)=size(var,2)
    dim(3)=size(var,3)
    !
    dimat=(/ia+1,ja+1,ka+1/)
    !
    dimt=dim
    offset=(/ig0,jg0,kg0/)
    !
    call h5screate_simple_f(3,dimat,filespace,h5error)
    call h5dcreate_f(h5file_id,varname,H5T_NATIVE_DOUBLE,filespace,    &
                                                       dset_id,h5error)
    call h5screate_simple_f(3,dimt,memspace,h5error)
    call h5sclose_f(filespace,h5error)
    call h5dget_space_f(dset_id,filespace,h5error)
    call h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,offset,     &
                                                          dimt,h5error)
    call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,h5error) 
    !
    call h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,h5error)
    call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,var,dimt,h5error,       &
                    file_space_id=filespace,mem_space_id=memspace,    &
                                                     xfer_prp=plist_id)
    call h5sclose_f(filespace,h5error)
    call h5sclose_f(memspace,h5error)
    call h5dclose_f(dset_id,h5error)
    call h5pclose_f(plist_id,h5error)
    !
    if(lio) print*,' << ',varname
    !
#endif
    !
  end subroutine h5wa3d_r8_struct
  !
  subroutine h5wa3d_r8_stream(varname,var)
    !
    use parallel,only: lio,mpirankmax,ptabupd
    !
    ! arguments
    character(LEN=*),intent(in) :: varname
    real(8),intent(in) :: var(:,:,:)
    !
#ifdef HDF5
    ! local data
    real(8),allocatable :: var_1d(:)
    integer :: dima,jrank,i,j,k,n
    integer :: data_size(0:mpirankmax)
    integer :: h5error
    integer,save :: dim(1)
    !
    integer(hid_t) :: dset_id,filespace,memspace,plist_id
    integer(hsize_t),dimension(1),save :: offset
    integer(hsize_t),save :: dimt(1),dimat(1)
    integer,save :: size_sav(3)=0
    !
    if( size_sav(1)==size(var,1) .and. size_sav(2)==size(var,2) .and. & 
        size_sav(3)==size(var,3) ) then
      continue
    else
      !
      size_sav(1)=size(var,1)
      size_sav(2)=size(var,2)
      size_sav(3)=size(var,3)
      !
      dim(1)=size_sav(1)*size_sav(2)*size_sav(3)
      !
      call ptabupd(var=dim(1),table=data_size)
      !
      dima=0
      do jrank=0,mpirankmax
        dima=dima+data_size(jrank)
      enddo
      !
      offset(1)=0
      do jrank=0,mpirank-1
        offset(1)=offset(1)+data_size(jrank)
      enddo
      !
      dimat=dima
      !
      dimt=dim
    endif
    !
    allocate(var_1d(1:dim(1)))
    n=0
    do k=1,size_sav(3)
    do j=1,size_sav(2)
    do i=1,size_sav(1)
      n=n+1
      var_1d(n)=var(i,j,k)
    enddo
    enddo
    enddo
    !
    call h5screate_simple_f(1,dimat,filespace,h5error)
    call h5dcreate_f(h5file_id,varname,H5T_NATIVE_DOUBLE,filespace,    &
                                                       dset_id,h5error)
    call h5screate_simple_f(1,dimt,memspace,h5error)
    call h5sclose_f(filespace,h5error)
    call h5dget_space_f(dset_id,filespace,h5error)
    call h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,offset,       &
                                                          dimt,h5error)
    call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,h5error) 
    !
    call h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,h5error)
    call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,var_1d,dimt,h5error,    &
                    file_space_id=filespace,mem_space_id=memspace,    &
                                                     xfer_prp=plist_id)
    call h5sclose_f(filespace,h5error)
    call h5sclose_f(memspace,h5error)
    call h5dclose_f(dset_id,h5error)
    call h5pclose_f(plist_id,h5error)
    !
    deallocate(var_1d)
    !
    if(lio) print*,' << ',varname
    !
#endif
    !
  end subroutine h5wa3d_r8_stream
  !
  subroutine h5wa2d_r8(varname,var,dir)
    !
    use commvar, only: ia,ja,ka
    use parallel,only: lio,ig0,jg0,kg0
    !
    ! arguments
    character(LEN=*),intent(in) :: varname
    real(8),intent(in) :: var(:,:)
    character(len=1),intent(in) :: dir
    !
#ifdef HDF5
    ! local data
    integer :: dim(2),dima
    integer(hsize_t), dimension(2) :: offset
    integer :: h5error
    !
    integer(hid_t) :: dset_id,filespace,memspace,plist_id
    integer(hsize_t) :: dimt(2),dimat(2)
    !
    dim(1)=size(var,1)
    dim(2)=size(var,2)
    !
    if(dir=='i') then
      dimat=(/ja+1,ka+1/)
      !
      offset=(/jg0,kg0/)
    elseif(dir=='j') then
      dimat=(/ia+1,ka+1/)
      !
      offset=(/ig0,kg0/)
    elseif(dir=='k') then
      dimat=(/ia+1,ja+1/)
      !
      offset=(/ig0,jg0/)
    else
      stop ' !! error in dir @ h5wa2d_r8'
    endif
    !
    dimt=dim
    !
    call h5screate_simple_f(2,dimat,filespace,h5error)
    call h5dcreate_f(h5file_id,varname,H5T_NATIVE_DOUBLE,filespace,    &
                                                       dset_id,h5error)
    call h5screate_simple_f(2,dimt,memspace,h5error)
    call h5sclose_f(filespace,h5error)
    call h5dget_space_f(dset_id,filespace,h5error)
    call h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,offset,     &
                                                          dimt,h5error)
    call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,h5error) 
    !
    call h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,h5error)
    call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,var,dimt,h5error,       &
                    file_space_id=filespace,mem_space_id=memspace,    &
                                                     xfer_prp=plist_id)
    call h5sclose_f(filespace,h5error)
    call h5sclose_f(memspace,h5error)
    call h5dclose_f(dset_id,h5error)
    call h5pclose_f(plist_id,h5error)
    !
    if(dir=='i' .and. jrk==0 .and. krk==0) print*,' << ',varname
    if(dir=='j' .and. irk==0 .and. krk==0) print*,' << ',varname
    if(dir=='k' .and. irk==0 .and. jrk==0) print*,' << ',varname
    !
#endif
    !
  end subroutine h5wa2d_r8
  !+-------------------------------------------------------------------+
  !| This end of the function h5wa3d_r8.                               |
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
    integer(hsize_t) :: dims_h5(3)
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
  !+-------------------------------------------------------------------+
  !| This subroutine is used to write a 3D array to a hdf-5 file in    |
  !| the sequential mode.                                              |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 28-05-2021  | Created by J. Fang @ Warrington                     |
  !+-------------------------------------------------------------------+
  subroutine h5_write1int(var,varname,filename,explicit,newfile)
    !
    integer,intent(in) :: var
    character(len=*),intent(in) :: varname,filename
    logical,intent(in), optional:: explicit,newfile
    logical :: lexplicit,lnew
    logical :: lfilalive
    !
    integer :: v(1)
    !
#ifdef HDF5
    integer(hid_t) :: file_id
    ! file identifier
    integer(hid_t) :: dset_id1
    ! dataset identifier
    integer :: h5error ! error flag
    integer(hsize_t) :: dimt(1)
    !
    if (present(explicit)) then
       lexplicit = explicit
    else
       lexplicit = .true.
    end if
    if (present(newfile)) then
       lnew      = newfile
    else
       lnew = .false. 
    end if
    !
    v=var
    !
    call h5open_f(h5error)
    !
    if(lnew) then
      call h5fcreate_f(filename,h5f_acc_trunc_f,file_id,h5error)
    else
      call h5fopen_f(filename,h5f_acc_rdwr_f,file_id,h5error)
    endif
    !
    if(h5error.ne.0)  stop ' !! error in h5_write1int 1'
    !
    dimt=(/1/)
    !
    ! write the dataset.
    call h5ltmake_dataset_f(file_id,varname,1,dimt,h5t_native_integer,v,h5error)
    if(h5error .ne. 0) stop 'h5 write error h5_write1int'
    if(lexplicit) print*,' << ',varname,' to ',filename
    !
    call h5fclose_f(file_id, h5error)
    if(h5error.ne.0)  stop ' !! error in h5_write1int 2'
    !
    ! close fortran interface.
    call h5close_f(h5error)
    if(h5error.ne.0)  stop ' !! error in h5_write1int 3'
    !
#endif
    !
  end subroutine h5_write1int
  !
  subroutine h5_writearray1dint(var,varname,filename,explicit,newfile)
    !
    integer,intent(in) :: var(:)
    character(len=*),intent(in) :: varname,filename
    logical,intent(in), optional:: explicit,newfile
    logical :: lexplicit,lnew
    logical :: lfilalive
    integer :: dim1
    !
#ifdef HDF5
    integer(hid_t) :: file_id
    ! file identifier
    integer(hid_t) :: dset_id1
    ! dataset identifier
    integer :: h5error ! error flag
    integer(hsize_t) :: dimt(1)
    !
    if (present(explicit)) then
       lexplicit = explicit
    else
       lexplicit = .true.
    end if
    if (present(newfile)) then
       lnew      = newfile
    else
       lnew = .false. 
    end if
    !
    call h5open_f(h5error)
    !
    if(lnew) then
      call h5fcreate_f(filename,h5f_acc_trunc_f,file_id,h5error)
    else
      call h5fopen_f(filename,h5f_acc_rdwr_f,file_id,h5error)
    endif
    !
    if(h5error.ne.0)  stop ' !! error in h5_writearray1dint 1'
    !
    dim1=size(var)
    dimt=(/dim1+1/)
    !
    ! write the dataset.
    call h5ltmake_dataset_f(file_id,varname,1,dimt,h5t_native_integer,var,h5error)
    if(h5error .ne. 0) stop 'h5 write error h5_writearray1dint'
    if(lexplicit) print*,' << ',varname,' to ',filename
    !
    call h5fclose_f(file_id, h5error)
    !
    if(h5error.ne.0)  stop ' !! error in h5_writearray1dint 2'
    ! close fortran interface.
    call h5close_f(h5error)
    !
    if(h5error.ne.0)  stop ' !! error in h5_writearray1dint 3'
    !
#endif 
    !
  end subroutine h5_writearray1dint
  !
  subroutine h5_write1rl8(var,varname,filename,explicit,newfile)
    !
    real(8),intent(in) :: var
    character(len=*),intent(in) :: varname,filename
    logical,intent(in), optional:: explicit,newfile
    logical :: lexplicit,lnew
    logical :: lfilalive
    !
    real(8) :: v(1)
    !
#ifdef HDF5
    !
    integer(hid_t) :: file_id
    ! file identifier
    integer(hid_t) :: dset_id1
    ! dataset identifier
    integer :: h5error ! error flag
    integer(hsize_t) :: dimt(1)
    !
    if (present(explicit)) then
       lexplicit = explicit
    else
       lexplicit = .true.
    end if
    if (present(newfile)) then
       lnew      = newfile
    else
       lnew = .false. 
    end if
    !
    v=var
    !
    call h5open_f(h5error)
    !
    if(lnew) then
      call h5fcreate_f(filename,h5f_acc_trunc_f,file_id,h5error)
    else
      call h5fopen_f(filename,h5f_acc_rdwr_f,file_id,h5error)
    endif
    !
    if(h5error.ne.0)  stop ' !! error in h5_write1rl8 1'
    !
    dimt=(/1/)
    !
    ! write the dataset.
    call h5ltmake_dataset_f(file_id,varname,1,dimt,h5t_native_double,v,h5error)
    if(h5error .ne. 0) stop 'h5 write error h5_write1rl8'
    if(lexplicit) print*,' << ',varname,' to ',filename
    !
    call h5fclose_f(file_id, h5error)
    if(h5error.ne.0)  stop ' !! error in h5_write1rl8 2'
    !
    ! close fortran interface.
    call h5close_f(h5error)
    if(h5error.ne.0)  stop ' !! error in h5_write1rl8 3'
    !
#endif 
    !
  end subroutine h5_write1rl8
  !
  subroutine h5_writearray1d(var,varname,filename,explicit,newfile)
    !
    real(8),intent(in) :: var(:)
    character(len=*),intent(in) :: varname,filename
    logical,intent(in), optional:: explicit,newfile
    logical :: lexplicit,lnew
    logical :: lfilalive
    integer :: dim1
    !
#ifdef HDF5
    !
    integer(hid_t) :: file_id
    ! file identifier
    integer(hid_t) :: dset_id1
    ! dataset identifier
    integer :: h5error ! error flag
    integer(hsize_t) :: dimt(1)
    !
    if (present(explicit)) then
       lexplicit = explicit
    else
       lexplicit = .true.
    end if
    if (present(newfile)) then
       lnew      = newfile
    else
       lnew = .false. 
    end if
    !
    call h5open_f(h5error)
    !
    if(lnew) then
      call h5fcreate_f(filename,h5f_acc_trunc_f,file_id,h5error)
    else
      call h5fopen_f(filename,h5f_acc_rdwr_f,file_id,h5error)
    endif
    !
    if(h5error.ne.0)  stop ' !! error in h5_writearray1d 1'
    !
    dim1=size(var)
    dimt=(/dim1/)
    !
    ! write the dataset.
    call h5ltmake_dataset_f(file_id,varname,1,dimt,h5t_native_double,var,h5error)
    if(h5error .ne. 0) stop 'h5 write error h5_writearray1d'
    if(lexplicit) print*,' << ',varname,' to ',filename
    !
    call h5fclose_f(file_id, h5error)
    if(h5error.ne.0)  stop ' !! error in h5_writearray1d 2'
    !
    ! close fortran interface.
    call h5close_f(h5error)
    if(h5error.ne.0)  stop ' !! error in h5_writearray1d 3'
    !
#endif 
    !
  end subroutine h5_writearray1d
  !
  subroutine h5_writearray3d(var,varname,filename,explicit,newfile)
    !
    real(8),intent(in) :: var(:,:,:)
    character(len=*),intent(in) :: varname,filename
    logical,intent(in), optional:: explicit,newfile
    !
    integer :: dim1,dim2,dim3
    logical :: lexplicit,lnew
    logical :: lfilalive
    !
#ifdef HDF5
    !
    integer(hid_t) :: file_id
    ! file identifier
    integer(hid_t) :: dset_id1
    ! dataset identifier
    integer :: h5error ! error flag
    integer(hsize_t) :: dimt(3)
    !
    if (present(explicit)) then
       lexplicit = explicit
    else
       lexplicit = .true.
    end if
    if (present(newfile)) then
       lnew      = newfile
    else
       lnew = .false. 
    end if
    !
    dim1=size(var,1)
    dim2=size(var,2)
    dim3=size(var,3)
    !
    call h5open_f(h5error)
    !
    if(lnew) then
      call h5fcreate_f(filename,h5f_acc_trunc_f,file_id,h5error)
    else
      call h5fopen_f(filename,h5f_acc_rdwr_f,file_id,h5error)
    endif
    !
    if(h5error.ne.0)  stop ' !! error in h5_writearray3d 1'
    !
    dimt=(/dim1,dim2,dim3/)
    !
    ! write the dataset.
    call h5ltmake_dataset_f(file_id,varname,3,dimt,h5t_native_double,var,h5error)
    if(h5error .ne. 0) stop 'h5 write error h5_writearray3d'
    if(lexplicit) print*,' << ',varname,' to ',filename
    if(h5error.ne.0)  stop ' !! error in h5_writearray3d 2'
    !
    call h5fclose_f(file_id, h5error)
    if(h5error.ne.0)  stop ' !! error in h5_writearray3d 3'
    !
    ! close fortran interface.
    call h5close_f(h5error)
    if(h5error.ne.0)  stop ' !! error in h5_writearray3d 4'
    !
#endif 
    !
  end subroutine h5_writearray3d
  !+-------------------------------------------------------------------+
  !| This end of the subroutine h5_writearray3d.                       |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to remove a dataset from a hdf5 file.     |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 31-03-2022  | Created by J. Fang @ Warrington                     |
  !+-------------------------------------------------------------------+
  subroutine h5delete(varname)
    !
    character(len=*),intent(in) :: varname
    !
    ! local data
    integer :: h5error
    !
#ifdef HDF5
    call h5ldelete_f(h5file_id,varname,h5error)
    !
    print*,' ** dataset: ',varname,'removed'
    !
#endif 
    !
  end subroutine h5delete
  !+-------------------------------------------------------------------+
  !| This end of the subroutine h5delete.                              |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to read 1-D array via hdf5 interface.     |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 31-03-2022  | Created by J. Fang @ Warrington                     |
  !+-------------------------------------------------------------------+
  subroutine h5_readarray1d(varname,var,dim,filename,explicit)
    !
    real(8),intent(out) :: var(:)
    integer,intent(in) :: dim
    character(len=*),intent(in) :: varname,filename
    logical,intent(in), optional:: explicit
    logical :: lexplicit
    !
#ifdef HDF5
    !
    integer(hid_t) :: file_id
    ! file identifier
    integer(hid_t) :: dset_id1
    ! dataset identifier
    integer :: h5error ! error flag
    integer(hsize_t) :: dimt(1)
    !
    if (present(explicit)) then
       lexplicit = explicit
    else
       lexplicit = .true.
    end if
    !
    call h5open_f(h5error)
    !
    call h5fopen_f(filename,h5f_acc_rdwr_f,file_id,h5error)

    ! open an existing dataset.
    call h5dopen_f(file_id,varname,dset_id1,h5error)
    !
    dimt=(/dim/)
    !
    ! read the dataset.
    call h5dread_f(dset_id1,h5t_native_double,var,dimt,h5error)

    if(h5error.ne.0)  stop ' !! error in h5_readarray1d 1'
    !
    ! close the dataset
    call h5dclose_f(dset_id1, h5error)
    if(h5error.ne.0)  stop ' !! error in h5_readarray1d 2'
    ! close the file.
    call h5fclose_f(file_id, h5error)
    if(h5error.ne.0)  stop ' !! error in h5_readarray1d 3'
    !
    ! close fortran interface.
    call h5close_f(h5error)
    if(h5error.ne.0)  stop ' !! error in h5_readarray1d 4'
    !
    if(lexplicit)  print*,' >> ',varname,' from ',filename,' ... done'
    !
#endif 
    !
  end subroutine h5_readarray1d
  !
  subroutine h5_read1int(var,varname,filename,explicit)
    !
    integer,intent(out) :: var
    character(len=*),intent(in) :: varname,filename
    logical,intent(in), optional:: explicit
    logical :: lexplicit
    !
#ifdef HDF5
    !
    integer(hid_t) :: file_id
    ! file identifier
    integer(hid_t) :: dset_id1
    ! dataset identifier
    integer :: v(1)
    integer :: h5error ! error flag
    integer(hsize_t) :: dimt(1)
    !
    if (present(explicit)) then
       lexplicit = explicit
    else
       lexplicit = .true.
    end if
    !
    dimt=(/1/)
    !
    call h5open_f(h5error)
    print*,' ** open hdf5 interface'
    !
    call h5fopen_f(filename,h5f_acc_rdwr_f,file_id,h5error)
    !
    call h5ltread_dataset_f(file_id,varname,h5t_native_integer,v,dimt,h5error)
    !
    call h5fclose_f(file_id,h5error)
    !
    if(h5error.ne.0)  stop ' !! error in h5_readarray1dint 1'
    !
    ! close fortran interface.
    call h5close_f(h5error)
    !
    var=v(1)
    if(h5error.ne.0)  stop ' !! error in h5_readarray1dint 2'
    !
    if(lexplicit)  print*,' >> ',varname,' from ',filename,' ... done'
    !
#endif 
    !
  end subroutine h5_read1int
  !
  subroutine h5_read1rl8(var,varname,filename,explicit)
    !
    real(8),intent(out) :: var
    character(len=*),intent(in) :: varname,filename
    logical,intent(in), optional:: explicit
    logical :: lexplicit
    !
#ifdef HDF5
    !
    integer(hid_t) :: file_id
    ! file identifier
    integer(hid_t) :: dset_id1
    ! dataset identifier
    real(8) :: v(1)
    integer :: h5error ! error flag
    integer(hsize_t) :: dimt(1)
    !
    if (present(explicit)) then
       lexplicit = explicit
    else
       lexplicit = .true.
    end if
    !
    dimt=(/1/)
    !
    call h5open_f(h5error)
    if(lexplicit)  print*,' ** open hdf5 interface'
    !
    call h5fopen_f(filename,h5f_acc_rdwr_f,file_id,h5error)
    !
    call h5ltread_dataset_f(file_id,varname,h5t_native_double,v,dimt,h5error)
    !
    call h5fclose_f(file_id,h5error)
    !
    if(h5error.ne.0)  stop ' !! error in h5_readarray1dint 1'
    !
    ! close fortran interface.
    call h5close_f(h5error)
    !
    var=v(1)
    if(h5error.ne.0)  stop ' !! error in h5_readarray1dint 2'
    !
    if(lexplicit)  print*,' >> ',varname,' from ',filename,' ... done'
    !
#endif 
    !
  end subroutine h5_read1rl8
  !+-------------------------------------------------------------------+
  !| This end of the subroutine h5_readarray1d.                        |
  !+-------------------------------------------------------------------+
  !
  subroutine h5_readarray3d(varin,vname,dim1,dim2,dim3,fname,explicit)
    !
    integer :: dim1,dim2,dim3
    real(8) :: varin(0:dim1,0:dim2,0:dim3)
    character(len=*),intent(in) :: vname,fname
    logical,intent(in), optional:: explicit
    logical :: lexplicit
    !
#ifdef HDF5
    !
    integer(hid_t) :: file_id
    ! file identifier
    integer(hid_t) :: dset_id1
    ! dataset identifier
    integer :: h5error ! error flag
    integer(hsize_t) :: dimt(3)
    !
    if (present(explicit)) then
       lexplicit = explicit
    else
       lexplicit = .true.
    end if
    !
    call h5open_f(h5error)
    if(h5error.ne.0)  stop ' !! error in h5_readarray3d 1'
    !
    call h5fopen_f(fname,h5f_acc_rdwr_f,file_id,h5error)
    if(h5error.ne.0)  stop ' !! error in h5_readarray3d 2'

    ! open an existing dataset.
    call h5dopen_f(file_id,vname,dset_id1,h5error)
    if(h5error.ne.0)  stop ' !! error in h5_readarray3d 3'
    !
    dimt=(/dim1+1,dim2+1,dim3+1/)
    !
    ! read the dataset.
    call h5dread_f(dset_id1,h5t_native_double,varin,dimt,h5error)
    if(h5error.ne.0)  stop ' !! error in h5_readarray3d 4'
    !
    ! close the dataset
    call h5dclose_f(dset_id1, h5error)
    if(h5error.ne.0)  stop ' !! error in h5_readarray3d 5'
    ! close the file.
    call h5fclose_f(file_id, h5error)
    if(h5error.ne.0)  stop ' !! error in h5_readarray3d 6'
    !
    ! close fortran interface.
    call h5close_f(h5error)
    if(h5error.ne.0)  stop ' !! error in h5_readarray3d 7'
    !
    if(lexplicit) print*,' >> ',vname,' from ',fname,' ... done'
    !
#endif 
    !
  end subroutine h5_readarray3d
  !
  !+-------------------------------------------------------------------+
  !| this subroutine is to read a 2D slice from 3D array               |
  !+-------------------------------------------------------------------+
  subroutine h5_read2dfrom3d(vread,dim1,dim2,dim3,vname,fname,         &
                                                   islice,jslice,kslice)
    !
    ! argument2
    integer,intent(in) :: dim1,dim2,dim3
    real(8),allocatable,intent(out) :: vread(:,:)
    character(len=*),intent(in) :: vname,fname
    integer,intent(in),optional :: islice,jslice,kslice
    !
#ifdef HDF5
    !
    ! local data
    integer :: error                          ! error flag
    integer(hid_t) :: file_id                 ! file identifier 
    integer(hid_t) :: dataspace               ! dataspace identifier 
    integer(hid_t) :: memspace                ! memspace identifier 
    integer(hid_t) :: dset_id                 ! dataset identifier 
    integer :: rank                           ! dataset rank ( in file )
    integer(hsize_t) :: offset(1:3)           ! hyperslab offset
    integer(hsize_t) :: dimsm(1:3)            ! dataset dimensions
    !
    ! Initialize FORTRAN interface.
    call h5open_f(error) 
    !
    ! Open the file.
    call h5fopen_f(fname,h5f_acc_rdwr_f,file_id,error)
    !
    ! Open the dataset.
    call h5dopen_f(file_id,vname,dset_id,error)
    !
    if(present(islice)) then
      offset=(/islice,0,0/)
      dimsm=(/1,dim2+1,dim3+1/)
      !
      allocate(vread(0:dim2,0:dim3))
    elseif(present(jslice)) then
      offset=(/0,jslice,0/)
      dimsm=(/dim1+1,1,dim3+1/)
      !
      allocate(vread(0:dim1,0:dim3))
    elseif(present(kslice)) then
      offset=(/0,0,kslice/)
      dimsm=(/dim1+1,dim2+1,1/)
      !
      allocate(vread(0:dim1,0:dim2))
    else
      stop ' !! slice set error !!'
    endif
    !
    rank=3 
    ! Create memory dataspace.
    call h5screate_simple_f(rank,dimsm,memspace,error)
    !
    ! Get dataset's dataspace identifier and select subset.
    call h5dget_space_f(dset_id,dataspace,error)
    call h5sselect_hyperslab_f(dataspace,h5s_select_set_f,offset,dimsm,error) 
    !
    call h5dread_f(dset_id,h5t_native_double,vread,dimsm,error,        &
                   mem_space_id=memspace,file_space_id=dataspace)
    !
    call h5sclose_f(dataspace, error)
    call h5sclose_f(memspace, error)
    call h5dclose_f(dset_id, error)
    call h5fclose_f(file_id, error)
    !
    ! Close FORTRAN interface.
    !
    call h5close_f(error)
    !
    if(error==0) then
      print*,' >> ',vname,' from ',fname,' ... done'
    else
      stop ' !! error in h5_read2dfrom3d'
    endif
    !
#endif 
    !
  end subroutine h5_read2dfrom3d
  !
  subroutine h5_read1dfrom3d(vread,dim1,dim2,dim3,vname,fname,         &
                                                   islice,jslice,kslice)
    !
    ! argument2
    integer,intent(in) :: dim1,dim2,dim3
    real(8),allocatable,intent(out) :: vread(:)
    character(len=*),intent(in) :: vname,fname
    integer,intent(in),optional :: islice,jslice,kslice
    !
#ifdef HDF5
    ! local data
    integer :: error                          ! error flag
    integer(hid_t) :: file_id                 ! file identifier 
    integer(hid_t) :: dataspace               ! dataspace identifier 
    integer(hid_t) :: memspace                ! memspace identifier 
    integer(hid_t) :: dset_id                 ! dataset identifier 
    integer :: rank                           ! dataset rank ( in file )
    integer(hsize_t) :: offset(1:3)           ! hyperslab offset
    integer(hsize_t) :: dimsm(1:3)            ! dataset dimensions
    !
    ! Initialize FORTRAN interface.
    call h5open_f(error) 
    !
    ! Open the file.
    call h5fopen_f(fname,h5f_acc_rdwr_f,file_id,error)
    !
    ! Open the dataset.
    call h5dopen_f(file_id,vname,dset_id,error)
    !
    if(present(islice) .and. present(jslice)) then
      offset=(/islice+1,jslice+1,0/)
      dimsm=(/1,1,dim3+1/)
      !
      allocate(vread(0:dim3))
    elseif(present(islice) .and. present(kslice)) then
      offset=(/islice+1,0,kslice+1/)
      dimsm=(/1,dim2+1,1/)
      !
      allocate(vread(0:dim2))
    elseif(present(jslice) .and. present(kslice)) then
      offset=(/0,jslice+1,kslice+1/)
      dimsm=(/dim1+1,1,1/)
      !
      allocate(vread(0:dim1))
    else
      stop ' !! slice set error !!'
    endif
    !
    rank=3 
    ! Create memory dataspace.
    call h5screate_simple_f(rank,dimsm,memspace,error)
    !
    ! Get dataset's dataspace identifier and select subset.
    call h5dget_space_f(dset_id,dataspace,error)
    call h5sselect_hyperslab_f(dataspace,h5s_select_set_f,offset,dimsm,error) 
    !
    call h5dread_f(dset_id,h5t_native_double,vread,dimsm,error,        &
                   mem_space_id=memspace,file_space_id=dataspace)
    !
    call h5sclose_f(dataspace, error)
    call h5sclose_f(memspace, error)
    call h5dclose_f(dset_id, error)
    call h5fclose_f(file_id, error)
    !
    ! Close FORTRAN interface.
    !
    call h5close_f(error)
    !
    if(error==0) then
      print*,' >> ',vname,' from ',fname,' ... done'
    else
      stop ' !! error in h5_read1dfrom3d'
    endif
    !
#endif 
    !
  end subroutine h5_read1dfrom3d

  subroutine h5_read1dfrom2d(vread,dim1,dim2,vname,fname,islice,jslice)
    !
    ! argument2
    integer,intent(in) :: dim1,dim2
    real(8),allocatable,intent(out) :: vread(:)
    character(len=*),intent(in) :: vname,fname
    integer,intent(in),optional :: islice,jslice
    !
#ifdef HDF5
    ! local data
    integer :: error                          ! error flag
    integer(hid_t) :: file_id                 ! file identifier 
    integer(hid_t) :: dataspace               ! dataspace identifier 
    integer(hid_t) :: memspace                ! memspace identifier 
    integer(hid_t) :: dset_id                 ! dataset identifier 
    integer :: rank                           ! dataset rank ( in file )
    integer(hsize_t) :: offset(1:2)           ! hyperslab offset
    integer(hsize_t) :: dimsm(1:2)            ! dataset dimensions
    !
    ! Initialize FORTRAN interface.
    call h5open_f(error) 
    !
    ! Open the file.
    call h5fopen_f(fname,h5f_acc_rdwr_f,file_id,error)
    !
    ! Open the dataset.
    call h5dopen_f(file_id,vname,dset_id,error)
    !
    if(present(islice)) then
      offset=(/0,jslice+1/)
      dimsm=(/1,dim2+1/)
      !
      allocate(vread(0:dim2))
    elseif(present(jslice)) then
      offset=(/islice+1,0/)
      dimsm=(/dim1+1,1/)
      !
      allocate(vread(0:dim1))
    else
      stop ' !! slice set error !!'
    endif
    !
    rank=2 
    ! Create memory dataspace.
    call h5screate_simple_f(rank,dimsm,memspace,error)
    !
    ! Get dataset's dataspace identifier and select subset.
    call h5dget_space_f(dset_id,dataspace,error)
    call h5sselect_hyperslab_f(dataspace,h5s_select_set_f,offset,dimsm,error) 
    !
    call h5dread_f(dset_id,h5t_native_double,vread,dimsm,error,        &
                   mem_space_id=memspace,file_space_id=dataspace)
    !
    call h5sclose_f(dataspace, error)
    call h5sclose_f(memspace, error)
    call h5dclose_f(dset_id, error)
    call h5fclose_f(file_id, error)
    !
    ! Close FORTRAN interface.
    !
    call h5close_f(error)
    !
    if(error==0) then
      print*,' >> ',vname,' from ',fname,' ... done'
    else
      stop ' !! error in h5_read1dfrom2d'
    endif
    !
#endif 
    !
  end subroutine h5_read1dfrom2d
  !+-------------------------------------------------------------------+
  !! the end of the subroutine h5_read2dfrom3d                         |
  !+-------------------------------------------------------------------+
  !!
end module hdf5io