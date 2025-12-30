!+---------------------------------------------------------------------+
!| This module is to write field via h5 interface.                     |
!+---------------------------------------------------------------------+
!| Writen by Jian Fang, 2020-04-16                                     |
!+---------------------------------------------------------------------+
module pastr_h5io
  !
  use hdf5
  use h5lt
  !
  implicit none
  !
  Interface H5ReadArray
    !
    module procedure h5_read1int
    module procedure h5_readarray1dint
    module procedure h5_read1rl8
    module procedure h5_readarray1d
    module procedure h5_readarray2d
    module procedure h5_readarray3d
    module procedure h5_read1rl4
    module procedure h5_readarray1d_r4
    module procedure h5_readarray2d_r4
    module procedure h5_readarray3d_r4
    !
  end Interface H5ReadArray
  !
  Interface H5ReadSubset
    !
    module procedure h5_read2dfrom3d
    module procedure h5_read2dfrom3d_r4
    module procedure h5_read1dfrom3d
    module procedure h5_read1dfrom3d_r4
    module procedure h5_read1dfrom2d
    !
  end Interface H5ReadSubset
  !
  Interface H5WriteArray
    !
    module procedure h5_writestring
    module procedure h5_write1int
    module procedure h5_write1rl8
    module procedure h5_write1rl4
    module procedure h5_writearray1dint
    module procedure h5_writearray1d
    module procedure h5_writearray1d_r4
    module procedure h5_writearray2d
    module procedure h5_writearray2d_r4
    module procedure h5_writearray3d
    module procedure h5_writearray3d_r4
    module procedure h5_writearray3d_int
    !
  end Interface H5WriteArray
  !
  contains
  !
  !+-------------------------------------------------------------------+
  !|this subroutine is used to read a array file with hdf5
  !+-------------------------------------------------------------------+
  subroutine h5_read1int(varin,vname,fname,explicit)
    !
    integer,intent(out) :: varin
    character(len=*),intent(in) :: vname,fname
    logical,intent(in), optional:: explicit
    logical :: lexplicit
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
    call h5fopen_f(fname,h5f_acc_rdwr_f,file_id,h5error)
    !
    call h5ltread_dataset_f(file_id,vname,h5t_native_integer,v,dimt,h5error)
    !
    call h5fclose_f(file_id,h5error)
    !
    if(h5error.ne.0)  stop ' !! error in h5_readarray1dint 1'
    !
    ! close fortran interface.
    call h5close_f(h5error)
    !
    varin=v(1)
    if(h5error.ne.0)  stop ' !! error in h5_readarray1dint 2'
    !
    if(lexplicit)  print*,' >> ',vname,' from ',fname,' ... done'
    !
  end subroutine h5_read1int
  !
  subroutine h5_read1rl8(varin,vname,fname,explicit)
    !
    real(8),intent(out) :: varin
    character(len=*),intent(in) :: vname,fname
    logical,intent(in), optional:: explicit
    logical :: lexplicit
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
    call h5fopen_f(fname,h5f_acc_rdwr_f,file_id,h5error)
    !
    call h5ltread_dataset_f(file_id,vname,h5t_native_double,v,dimt,h5error)
    !
    call h5fclose_f(file_id,h5error)
    !
    if(h5error.ne.0)  stop ' !! error in h5_readarray1dint 1'
    !
    ! close fortran interface.
    call h5close_f(h5error)
    !
    varin=v(1)
    if(h5error.ne.0)  stop ' !! error in h5_readarray1dint 2'
    !
    if(lexplicit)  print*,' >> ',vname,' from ',fname,' ... done'
    !
  end subroutine h5_read1rl8
  

  subroutine h5_read1rl4(varin,vname,fname,explicit)
    !
    real(4),intent(out) :: varin
    character(len=*),intent(in) :: vname,fname
    logical,intent(in), optional:: explicit
    logical :: lexplicit
    !
    integer(hid_t) :: file_id
    ! file identifier
    integer(hid_t) :: dset_id1
    ! dataset identifier
    real(4) :: v(1)
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
    call h5fopen_f(fname,h5f_acc_rdwr_f,file_id,h5error)
    !
    call h5ltread_dataset_f(file_id,vname,h5t_native_real,v,dimt,h5error)
    !
    call h5fclose_f(file_id,h5error)
    !
    if(h5error.ne.0)  stop ' !! error in h5_read1rl4 1'
    !
    ! close fortran interface.
    call h5close_f(h5error)
    !
    varin=v(1)
    if(h5error.ne.0)  stop ' !! error in h5_read1rl4 2'
    !
    if(lexplicit)  print*,' >> ',vname,' from ',fname,' ... done'
    !
  end subroutine h5_read1rl4

  subroutine h5_readarray1dint(varin,dim1,vname,fname,explicit)
    !
    integer :: dim1
    integer :: varin(0:dim1)
    character(len=*),intent(in) :: vname,fname
    logical,intent(in), optional:: explicit
    logical :: lexplicit
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
    dimt=(/dim1+1/)
    !
    call h5open_f(h5error)
    print*,' ** open hdf5 interface'
    !
    call h5fopen_f(fname,h5f_acc_rdwr_f,file_id,h5error)
    !
    call h5ltread_dataset_f(file_id,vname,h5t_native_integer,          &
                                                     varin,dimt,h5error)


    !
    call h5fclose_f(file_id,h5error)
    !
    if(h5error.ne.0)  stop ' !! error in h5_readarray1dint 1'
    !
    ! close fortran interface.
    call h5close_f(h5error)
    !
    if(h5error.ne.0)  stop ' !! error in h5_readarray1dint 2'
    !
    if(lexplicit)  print*,' >> ',vname,' from ',fname,' ... done'
    !
  end subroutine h5_readarray1dint
  !
  subroutine h5_readarray1d(varin,dim1,vname,fname,explicit)
    !
    integer :: dim1
    real(8) :: varin(0:dim1)
    character(len=*),intent(in) :: vname,fname
    logical,intent(in), optional:: explicit
    logical :: lexplicit
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
    call h5fopen_f(fname,h5f_acc_rdwr_f,file_id,h5error)

    ! open an existing dataset.
    call h5dopen_f(file_id,vname,dset_id1,h5error)
    !
    dimt=(/dim1+1/)
    !
    ! read the dataset.
    call h5dread_f(dset_id1,h5t_native_double,varin,dimt,h5error)

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
    if(lexplicit)  print*,' >> ',vname,' from ',fname,' ... done'
    !
  end subroutine h5_readarray1d
  !
  subroutine h5_readarray2d(varin,dim1,dim2,vname,fname,explicit)
    !
    integer :: dim1,dim2
    real(8) :: varin(0:dim1,0:dim2)
    character(len=*),intent(in) :: vname,fname
    logical,intent(in), optional:: explicit
    logical :: lexplicit
    !
    integer(hid_t) :: file_id
    ! file identifier
    integer(hid_t) :: dset_id1
    ! dataset identifier
    integer :: h5error ! error flag
    integer(hsize_t) :: dimt(2)
    !
    if (present(explicit)) then
       lexplicit = explicit
    else
       lexplicit = .true.
    end if
    !
    call h5open_f(h5error)
    if(h5error.ne.0)  stop ' !! error in h5_readarray2d 1'
    !
    call h5fopen_f(fname,h5f_acc_rdwr_f,file_id,h5error)
    if(h5error.ne.0)  stop ' !! error in h5_readarray2d 2'

    ! open an existing dataset.
    call h5dopen_f(file_id,vname,dset_id1,h5error)
    if(h5error.ne.0)  stop ' !! error in h5_readarray2d 3'
    !
    dimt=(/dim1+1,dim2+1/)
    !
    ! read the dataset.
    call h5dread_f(dset_id1,h5t_native_double,varin,dimt,h5error)
    if(h5error.ne.0)  stop ' !! error in h5_readarray2d 4'
    !
    ! close the dataset
    call h5dclose_f(dset_id1, h5error)
    if(h5error.ne.0)  stop ' !! error in h5_readarray2d 5'
    ! close the file.
    call h5fclose_f(file_id, h5error)
    if(h5error.ne.0)  stop ' !! error in h5_readarray2d 6'
    !
    ! close fortran interface.
    call h5close_f(h5error)
    if(h5error.ne.0)  stop ' !! error in h5_readarray2d 7'
    !
    if(h5error==0) then
      if(lexplicit) print*,' >> ',vname,' from ',fname,' ... done'
      return
    else
      stop ' !! error in h5_writearray2d'
    endif
    !
    !
  end subroutine h5_readarray2d
  !
  subroutine h5_readarray3d(varin,dim1,dim2,dim3,vname,fname,explicit)
    !
    integer :: dim1,dim2,dim3
    real(8) :: varin(0:dim1,0:dim2,0:dim3)
    character(len=*),intent(in) :: vname,fname
    logical,intent(in), optional:: explicit
    logical :: lexplicit
    !
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
  end subroutine h5_readarray3d

  subroutine h5_readarray1d_r4(varin,dim1,vname,fname,explicit)
    !
    integer :: dim1
    real(4) :: varin(0:dim1)
    character(len=*),intent(in) :: vname,fname
    logical,intent(in), optional:: explicit
    logical :: lexplicit
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
    call h5fopen_f(fname,h5f_acc_rdwr_f,file_id,h5error)

    ! open an existing dataset.
    call h5dopen_f(file_id,vname,dset_id1,h5error)
    !
    dimt=(/dim1+1/)
    !
    ! read the dataset.
    call h5dread_f(dset_id1,h5t_native_real,varin,dimt,h5error)

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
    if(lexplicit)  print*,' >> ',vname,' from ',fname,' ... done'
    !
  end subroutine h5_readarray1d_r4
  !
  subroutine h5_readarray2d_r4(varin,dim1,dim2,vname,fname,explicit)
    !
    integer :: dim1,dim2
    real(4) :: varin(0:dim1,0:dim2)
    character(len=*),intent(in) :: vname,fname
    logical,intent(in), optional:: explicit
    logical :: lexplicit
    !
    integer(hid_t) :: file_id
    ! file identifier
    integer(hid_t) :: dset_id1
    ! dataset identifier
    integer :: h5error ! error flag
    integer(hsize_t) :: dimt(2)
    !
    if (present(explicit)) then
       lexplicit = explicit
    else
       lexplicit = .true.
    end if
    !
    call h5open_f(h5error)
    if(h5error.ne.0)  stop ' !! error in h5_readarray2d 1'
    !
    call h5fopen_f(fname,h5f_acc_rdwr_f,file_id,h5error)
    if(h5error.ne.0)  stop ' !! error in h5_readarray2d 2'

    ! open an existing dataset.
    call h5dopen_f(file_id,vname,dset_id1,h5error)
    if(h5error.ne.0)  stop ' !! error in h5_readarray2d 3'
    !
    dimt=(/dim1+1,dim2+1/)
    !
    ! read the dataset.
    call h5dread_f(dset_id1,h5t_native_real,varin,dimt,h5error)
    if(h5error.ne.0)  stop ' !! error in h5_readarray2d 4'
    !
    ! close the dataset
    call h5dclose_f(dset_id1, h5error)
    if(h5error.ne.0)  stop ' !! error in h5_readarray2d 5'
    ! close the file.
    call h5fclose_f(file_id, h5error)
    if(h5error.ne.0)  stop ' !! error in h5_readarray2d 6'
    !
    ! close fortran interface.
    call h5close_f(h5error)
    if(h5error.ne.0)  stop ' !! error in h5_readarray2d 7'
    !
    if(h5error==0) then
      if(lexplicit) print*,' >> ',vname,' from ',fname,' ... done'
      return
    else
      stop ' !! error in h5_writearray2d'
    endif
    !
    !
  end subroutine h5_readarray2d_r4
  !
  subroutine h5_readarray3d_r4(varin,dim1,dim2,dim3,vname,fname,explicit)
    !
    integer :: dim1,dim2,dim3
    real(4) :: varin(0:dim1,0:dim2,0:dim3)
    character(len=*),intent(in) :: vname,fname
    logical,intent(in), optional:: explicit
    logical :: lexplicit
    !
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
    call h5dread_f(dset_id1,h5t_native_real,varin,dimt,h5error)
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
  end subroutine h5_readarray3d_r4
  !+-------------------------------------------------------------------+
  ! the end of the subroutine h5_readarray3d
  !+-------------------------------------------------------------------+
  !!
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
  end subroutine h5_read2dfrom3d
  !
  subroutine h5_read2dfrom3d_r4(vread,dim1,dim2,dim3,vname,fname,         &
                                                   islice,jslice,kslice)
    !
    ! argument2
    integer,intent(in) :: dim1,dim2,dim3
    real(4),allocatable,intent(out) :: vread(:,:)
    character(len=*),intent(in) :: vname,fname
    integer,intent(in),optional :: islice,jslice,kslice
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
    call h5dread_f(dset_id,h5t_native_real,vread,dimsm,error,        &
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
  end subroutine h5_read2dfrom3d_r4

  subroutine h5_read2dfrom3dint(vread,dim1,dim2,dim3,vname,fname,      &
                                                   islice,jslice,kslice)
    !
    ! argument2
    integer,intent(in) :: dim1,dim2,dim3
    integer,allocatable,intent(out) :: vread(:,:)
    character(len=*),intent(in) :: vname,fname
    integer,intent(in),optional :: islice,jslice,kslice
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
      offset=(/islice+1,0,0/)
      dimsm=(/1,dim2+1,dim3+1/)
      !
      allocate(vread(0:dim2,0:dim3))
    elseif(present(jslice)) then
      offset=(/0,jslice+1,0/)
      dimsm=(/dim1+1,1,dim3+1/)
      !
      allocate(vread(0:dim1,0:dim3))
    elseif(present(kslice)) then
      offset=(/0,0,kslice+1/)
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
    call h5dread_f(dset_id,h5t_native_integer,vread,dimsm,error,        &
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
  end subroutine h5_read2dfrom3dint
  !+-------------------------------------------------------------------+
  !! the end of the subroutine h5_read2dfrom3d                         |
  !+-------------------------------------------------------------------+
  !!
  !+-------------------------------------------------------------------+
  !| this subroutine is to read a 1D slice from 3D array               |
  !+-------------------------------------------------------------------+
  subroutine h5_read1dfrom3d(vread,dim1,dim2,dim3,vname,fname,         &
                                                   islice,jslice,kslice)
    !
    ! argument2
    integer,intent(in) :: dim1,dim2,dim3
    real(8),allocatable,intent(out) :: vread(:)
    character(len=*),intent(in) :: vname,fname
    integer,intent(in),optional :: islice,jslice,kslice
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
  end subroutine h5_read1dfrom3d

  subroutine h5_read1dfrom3d_r4(vread,dim1,dim2,dim3,vname,fname,         &
                                                   islice,jslice,kslice)
    !
    ! argument2
    integer,intent(in) :: dim1,dim2,dim3
    real(4),allocatable,intent(out) :: vread(:)
    character(len=*),intent(in) :: vname,fname
    integer,intent(in),optional :: islice,jslice,kslice
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
    call h5dread_f(dset_id,h5t_native_real,vread,dimsm,error,        &
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
  end subroutine h5_read1dfrom3d_r4

  subroutine h5_read1dfrom2d(vread,dim1,dim2,vname,fname,islice,jslice)
    !
    ! argument2
    integer,intent(in) :: dim1,dim2
    real(8),allocatable,intent(out) :: vread(:)
    character(len=*),intent(in) :: vname,fname
    integer,intent(in),optional :: islice,jslice
    !
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
  end subroutine h5_read1dfrom2d
  !+-------------------------------------------------------------------+
  !! the end of the subroutine h5_read1dfrom3d                         |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  ! this subroutine is used to read a 1d array file with hdf5          !
  !+-------------------------------------------------------------------+
  !
  subroutine h5_writestring(varin,vname,fname,explicit,ierr)
    !
    character(len=*),intent(in) :: varin,vname,fname
    logical,intent(in), optional:: explicit
    integer,intent(out), optional:: ierr
    logical :: lexplicit
    logical :: lfilalive,dexists
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
    inquire(file=fname, exist=lfilalive)
    if(lfilalive) then
      call h5fopen_f(fname,H5F_ACC_RDWR_F,file_id,h5error)
      !
      call h5lexists_f(file_id,vname,dexists,h5error)
      !
      if(dexists) call h5ldelete_f(file_id,vname,h5error)
    else
      call h5fcreate_f(fname,H5F_ACC_TRUNC_F,file_id,h5error)
    end if
    !
    if(h5error.ne.0)  stop ' !! error in h5_writestring 1'
    !
    dimt=(/1/)
    !
    ! write the dataset.
    call h5ltmake_dataset_string_f(file_id,vname,varin,h5error)
    if(h5error .ne. 0) stop 'h5 write error h5_writestring'
    if(lexplicit) print*,' << ',vname,' to ',fname
    if(present(ierr)) ierr=h5error
    !
    call h5fclose_f(file_id, h5error)
    if(h5error.ne.0)  stop ' !! error in h5_writestring 2'
    !
    ! close fortran interface.
    call h5close_f(h5error)
    if(h5error.ne.0)  stop ' !! error in h5_write1int 3'
    !
  end subroutine h5_writestring
  !
  subroutine h5_write1int(varin,vname,fname,explicit,ierr)
    !
    integer,intent(in) :: varin
    character(len=*),intent(in) :: vname,fname
    logical,intent(in), optional:: explicit
    integer,intent(out), optional:: ierr
    logical :: lexplicit
    logical :: lfilalive,dexists
    !
    integer :: v(1)
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
    v=varin
    !
    call h5open_f(h5error)
    !
    inquire(file=fname, exist=lfilalive)
    if(lfilalive) then
      call h5fopen_f(fname,H5F_ACC_RDWR_F,file_id,h5error)
      !
      call h5lexists_f(file_id,vname,dexists,h5error)
      !
      if(dexists) call h5ldelete_f(file_id,vname,h5error)
    else
      call h5fcreate_f(fname,H5F_ACC_TRUNC_F,file_id,h5error)
    end if
    !
    if(h5error.ne.0)  stop ' !! error in h5_write1int 1'
    !
    dimt=(/1/)
    !
    ! write the dataset.
    call h5ltmake_dataset_f(file_id,vname,1,dimt,h5t_native_integer,v,h5error)
    if(h5error .ne. 0) stop 'h5 write error h5_write1int'
    if(lexplicit) print*,' << ',vname,' to ',fname
    if(present(ierr)) ierr=h5error
    !
    call h5fclose_f(file_id, h5error)
    if(h5error.ne.0)  stop ' !! error in h5_write1int 2'
    !
    ! close fortran interface.
    call h5close_f(h5error)
    if(h5error.ne.0)  stop ' !! error in h5_write1int 3'
    !
  end subroutine h5_write1int
  !
  subroutine h5_writearray1dint(varin,dim1,vname,fname,explicit,ierr)
    !
    integer,intent(in) :: dim1
    integer,intent(in) :: varin(0:dim1)
    character(len=*),intent(in) :: vname,fname
    logical,intent(in), optional:: explicit
    integer,intent(out), optional:: ierr
    logical :: lexplicit
    logical :: lfilalive
    !
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
    inquire(file=fname, exist=lfilalive)
    if(lfilalive) then
      call h5fopen_f(fname,H5F_ACC_RDWR_F,file_id,h5error)
    else
      call h5fcreate_f(fname,H5F_ACC_TRUNC_F,file_id,h5error)
    end if
    !
    if(h5error.ne.0)  stop ' !! error in h5_writearray1dint 1'
    !
    dimt=(/dim1+1/)
    !
    ! write the dataset.
    call h5ltmake_dataset_f(file_id,vname,1,dimt,h5t_native_integer,varin,h5error)
    if(h5error .ne. 0) stop 'h5 write error h5_writearray1dint'
    if(lexplicit) print*,' << ',vname,' to ',fname
    if(present(ierr)) ierr=h5error
    !
    call h5fclose_f(file_id, h5error)
    !
    if(h5error.ne.0)  stop ' !! error in h5_writearray1dint 2'
    ! close fortran interface.
    call h5close_f(h5error)
    !
    if(h5error.ne.0)  stop ' !! error in h5_writearray1dint 3'
    !
  end subroutine h5_writearray1dint
  !
  subroutine h5_write1rl4(varin,vname,fname,explicit,ierr)
    !
    real(4),intent(in) :: varin
    character(len=*),intent(in) :: vname,fname
    logical,intent(in), optional:: explicit
    integer,intent(out), optional:: ierr
    logical :: lexplicit
    logical :: lfilalive
    !
    real(4) :: v(1)
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
    v=varin
    !
    call h5open_f(h5error)
    !
    inquire(file=fname, exist=lfilalive)
    if(lfilalive) then
      call h5fopen_f(fname,H5F_ACC_RDWR_F,file_id,h5error)
    else
      call h5fcreate_f(fname,H5F_ACC_TRUNC_F,file_id,h5error)
    end if
    !
    if(h5error.ne.0)  stop ' !! error in h5_write1rl4 1'
    !
    dimt=(/1/)
    !
    ! write the dataset.
    call h5ltmake_dataset_f(file_id,vname,1,dimt,h5t_native_real,v,h5error)
    if(h5error .ne. 0) stop 'h5 write error h5_write1rl4'
    if(lexplicit) print*,' << ',vname,' to ',fname
    if(present(ierr)) ierr=h5error
    !
    call h5fclose_f(file_id, h5error)
    if(h5error.ne.0)  stop ' !! error in h5_write1rl4 2'
    !
    ! close fortran interface.
    call h5close_f(h5error)
    if(h5error.ne.0)  stop ' !! error in h5_write1rl4 3'
    !
  end subroutine h5_write1rl4

  subroutine h5_write1rl8(varin,vname,fname,explicit,ierr)
    !
    real(8),intent(in) :: varin
    character(len=*),intent(in) :: vname,fname
    logical,intent(in), optional:: explicit
    integer,intent(out), optional:: ierr
    logical :: lexplicit
    logical :: lfilalive
    !
    real(8) :: v(1)
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
    v=varin
    !
    call h5open_f(h5error)
    !
    inquire(file=fname, exist=lfilalive)
    if(lfilalive) then
      call h5fopen_f(fname,H5F_ACC_RDWR_F,file_id,h5error)
    else
      call h5fcreate_f(fname,H5F_ACC_TRUNC_F,file_id,h5error)
    end if
    !
    if(h5error.ne.0)  stop ' !! error in h5_write1rl8 1'
    !
    dimt=(/1/)
    !
    ! write the dataset.
    call h5ltmake_dataset_f(file_id,vname,1,dimt,h5t_native_double,v,h5error)
    if(h5error .ne. 0) stop 'h5 write error h5_write1rl8'
    if(lexplicit) print*,' << ',vname,' to ',fname
    if(present(ierr)) ierr=h5error
    !
    call h5fclose_f(file_id, h5error)
    if(h5error.ne.0)  stop ' !! error in h5_write1rl8 2'
    !
    ! close fortran interface.
    call h5close_f(h5error)
    if(h5error.ne.0)  stop ' !! error in h5_write1rl8 3'
    !
  end subroutine h5_write1rl8
  !
  subroutine h5_writearray1d(varin,dim1,vname,fname,explicit,ierr)
    !
    integer,intent(in) :: dim1
    real(8),intent(in) :: varin(0:dim1)
    character(len=*),intent(in) :: vname,fname
    logical,intent(in), optional:: explicit
    integer,intent(out), optional:: ierr
    logical :: lexplicit
    logical :: lfilalive
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
    !
    call h5open_f(h5error)
    !
    inquire(file=fname, exist=lfilalive)
    if(lfilalive) then
      call h5fopen_f(fname,H5F_ACC_RDWR_F,file_id,h5error)
    else
      call h5fcreate_f(fname,H5F_ACC_TRUNC_F,file_id,h5error)
    end if
    !
    if(h5error.ne.0)  stop ' !! error in h5_writearray1d 1'
    !
    dimt=(/dim1+1/)
    !
    ! write the dataset.
    call h5ltmake_dataset_f(file_id,vname,1,dimt,h5t_native_double,varin,h5error)
    if(h5error .ne. 0) stop 'h5 write error h5_writearray1d'
    if(lexplicit) print*,' << ',vname,' to ',fname
    if(present(ierr)) ierr=h5error
    !
    call h5fclose_f(file_id, h5error)
    if(h5error.ne.0)  stop ' !! error in h5_writearray1d 2'
    !
    ! close fortran interface.
    call h5close_f(h5error)
    if(h5error.ne.0)  stop ' !! error in h5_writearray1d 3'
    !
  end subroutine h5_writearray1d


  subroutine h5_writearray1d_r4(varin,dim1,vname,fname,explicit,ierr)
    !
    integer,intent(in) :: dim1
    real(4),intent(in) :: varin(0:dim1)
    character(len=*),intent(in) :: vname,fname
    logical,intent(in), optional:: explicit
    integer,intent(out), optional:: ierr
    logical :: lexplicit
    logical :: lfilalive
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
    !
    call h5open_f(h5error)
    !
    inquire(file=fname, exist=lfilalive)
    if(lfilalive) then
      call h5fopen_f(fname,H5F_ACC_RDWR_F,file_id,h5error)
    else
      call h5fcreate_f(fname,H5F_ACC_TRUNC_F,file_id,h5error)
    end if
    !
    if(h5error.ne.0)  stop ' !! error in h5_writearray1d 1'
    !
    dimt=(/dim1+1/)
    !
    ! write the dataset.
    call h5ltmake_dataset_f(file_id,vname,1,dimt,h5t_native_real,varin,h5error)
    if(h5error .ne. 0) stop 'h5 write error h5_writearray1d'
    if(lexplicit) print*,' << ',vname,' to ',fname
    if(present(ierr)) ierr=h5error
    !
    call h5fclose_f(file_id, h5error)
    if(h5error.ne.0)  stop ' !! error in h5_writearray1d 2'
    !
    ! close fortran interface.
    call h5close_f(h5error)
    if(h5error.ne.0)  stop ' !! error in h5_writearray1d 3'
    !
  end subroutine h5_writearray1d_r4
  !+-------------------------------------------------------------------+
  !|the end of the subroutine h5_readarray1d                           |
  !+-------------------------------------------------------------------+
  !!
  !+-------------------------------------------------------------------+
  ! this subroutine is used to read a 2d array file with hdf5
  !+-------------------------------------------------------------------+
  subroutine h5_writearray2d(varin,dim1,dim2,vname,fname,explicit,ierr)
    !
    integer,intent(in) :: dim1,dim2
    real(8),intent(in) :: varin(0:dim1,0:dim2)
    character(len=*),intent(in) :: vname,fname
    logical,intent(in), optional:: explicit
    integer,intent(out), optional:: ierr
    logical :: lexplicit
    logical :: lfilalive
    !
    integer(hid_t) :: file_id
    ! file identifier
    integer(hid_t) :: dset_id1
    ! dataset identifier
    integer :: h5error ! error flag
    integer(hsize_t) :: dimt(2)
    !
    if (present(explicit)) then
       lexplicit = explicit
    else
       lexplicit = .true.
    end if
    !
    !
    call h5open_f(h5error)
    !
    if(h5error.ne.0)  stop ' !! error in h5_writearray2d 1'
    !
    inquire(file=fname, exist=lfilalive)
    if(lfilalive) then
      call h5fopen_f(fname,H5F_ACC_RDWR_F,file_id,h5error)
    else
      call h5fcreate_f(fname,H5F_ACC_TRUNC_F,file_id,h5error)
    end if
    !
    if(h5error.ne.0)  stop ' !! error in h5_writearray2d 2'
    !
    dimt=(/dim1+1,dim2+1/)
    !
    ! write the dataset.
    call h5ltmake_dataset_f(file_id,vname,2,dimt,h5t_native_double,varin,h5error)
    if(h5error .ne. 0) stop 'h5 write error h5_writearray2d'
    if(lexplicit) print*,' << ',vname,' to ',fname
    if(present(ierr)) ierr=h5error
    !
    call h5fclose_f(file_id, h5error)
    if(h5error.ne.0)  stop ' !! error in h5_writearray2d 3'
    !
    ! close fortran interface.
    call h5close_f(h5error)
    if(h5error.ne.0)  stop ' !! error in h5_writearray2d 4'
    !
  end subroutine h5_writearray2d

  subroutine h5_writearray2d_r4(varin,dim1,dim2,vname,fname,explicit,ierr)
    !
    integer,intent(in) :: dim1,dim2
    real(4),intent(in) :: varin(0:dim1,0:dim2)
    character(len=*),intent(in) :: vname,fname
    logical,intent(in), optional:: explicit
    integer,intent(out), optional:: ierr
    logical :: lexplicit
    logical :: lfilalive
    !
    integer(hid_t) :: file_id
    ! file identifier
    integer(hid_t) :: dset_id1
    ! dataset identifier
    integer :: h5error ! error flag
    integer(hsize_t) :: dimt(2)
    !
    if (present(explicit)) then
       lexplicit = explicit
    else
       lexplicit = .true.
    end if
    !
    !
    call h5open_f(h5error)
    !
    if(h5error.ne.0)  stop ' !! error in h5_writearray2d 1'
    !
    inquire(file=fname, exist=lfilalive)
    if(lfilalive) then
      call h5fopen_f(fname,H5F_ACC_RDWR_F,file_id,h5error)
    else
      call h5fcreate_f(fname,H5F_ACC_TRUNC_F,file_id,h5error)
    end if
    !
    if(h5error.ne.0)  stop ' !! error in h5_writearray2d 2'
    !
    dimt=(/dim1+1,dim2+1/)
    !
    ! write the dataset.
    call h5ltmake_dataset_f(file_id,vname,2,dimt,h5t_native_real,varin,h5error)
    if(h5error .ne. 0) stop 'h5 write error h5_writearray2d'
    if(lexplicit) print*,' << ',vname,' to ',fname
    if(present(ierr)) ierr=h5error
    !
    call h5fclose_f(file_id, h5error)
    if(h5error.ne.0)  stop ' !! error in h5_writearray2d 3'
    !
    ! close fortran interface.
    call h5close_f(h5error)
    if(h5error.ne.0)  stop ' !! error in h5_writearray2d 4'
    !
  end subroutine h5_writearray2d_r4
  !+-------------------------------------------------------------------+
  !|the end of the subroutine h5_readarray2d                           |
  !+-------------------------------------------------------------------+
  !!
  subroutine h5_writesubset(var,dim1,dim2,dim3,offset,varname,filename,initialise)
    !
    real(8),intent(in) :: var(:,:,:)
    integer,intent(in) :: dim1,dim2,dim3
    integer,intent(in) :: offset(3)
    character(len=*),intent(in) :: varname,filename
    logical,intent(in),optional :: initialise
    !
    ! local data
    logical :: linit,lfilalive
    integer :: h5error              ! error flag
    integer(hid_t) :: file_id       ! file identifier 
    integer(hid_t) :: dset_id       ! dataset identifier
    integer(hid_t) :: dataspace     ! dataspace identifier 
    integer(hid_t) :: memspace      ! memspace identifier  
    !
    integer(hsize_t) :: dimt(3)
    integer(hsize_t) :: dim_subset(3)
    integer(hsize_t) :: offseti8(3)
    !
    real(8),allocatable :: inidata(:,:,:)
    !
    if(present(initialise)) then
      linit=initialise
    else
      linit=.false.
    endif
    !
    dimt=(/dim1+1,dim2+1,dim3+1/)
    !
    ! Initialize FORTRAN interface. 
    call h5open_f(h5error) 
    !
    if(h5error.ne.0)  stop ' !! error in h5_writesubset: Initialize FORTRAN interface'
    !
    if(linit) then
      ! Create a new file using default properties.
      inquire(file=filename,exist=lfilalive)
      if(lfilalive) then
        call h5fopen_f(filename,H5F_ACC_RDWR_F,file_id,h5error)
        print*,' ** creat file:',filename
      else
        call h5fcreate_f(filename,h5f_acc_trunc_f,file_id,h5error)
        print*,' ** open file:',filename
      end if
      !
      if(h5error.ne.0)  stop ' !! error in h5_writesubset: h5fcreate_f'
      ! Create the data space for the  dataset. 
      call h5screate_simple_f(3,dimt,dataspace,h5error)
      if(h5error.ne.0)  stop ' !! error in h5_writesubset: h5screate_simple_f'
      ! Create the dataset with default properties.
      call h5dcreate_f(file_id,varname,h5t_native_double,dataspace,dset_id,h5error)
      if(h5error.ne.0)  stop ' !! error in h5_writesubset: h5dcreate_f'
      !
      ! Write the initial data set
      allocate(inidata(0:dim1,0:dim2,0:dim3))
      inidata=0.d0
      !
      call h5dwrite_f(dset_id,h5t_native_double,inidata,dimt,h5error)
      deallocate(inidata)
      !
      if(h5error.ne.0)  stop ' !! error in h5_writesubset: h5dwrite_f 1'
      !
      call h5sclose_f(dataspace,h5error)
      if(h5error.ne.0)  stop ' !! error in h5_writesubset: h5sclose_f'
      call h5dclose_f(dset_id,h5error)
      if(h5error.ne.0)  stop ' !! error in h5_writesubset: h5dclose_f'
      call h5fclose_f(file_id,h5error)
      if(h5error.ne.0)  stop ' !! error in h5_writesubset: h5fclose_f'
      !
      print*,' ** initialise h5 data: ',varname
      !
    endif
    !
    call h5fopen_f(filename,h5f_acc_rdwr_f,file_id,h5error)
    if(h5error.ne.0)  stop ' !! error in h5_writesubset: h5fopen_f'
    !
    ! Open the  dataset.
    CALL h5dopen_f(file_id,varname,dset_id,h5error)
    !
    print*,' ** write data: ',varname
    !
    if(h5error.ne.0)  stop ' !! error in h5_writesubset: Open the  dataset.'
    !
    ! Get dataset's dataspace identifier and select subset.
    call h5dget_space_f(dset_id,dataspace,h5error)
    !
    if(h5error.ne.0)  stop ' !! error in h5_writesubset: Get datasets dataspace identifier and select subset.'
    !
    dim_subset(1)=size(var,1)
    dim_subset(2)=size(var,2)
    dim_subset(3)=size(var,3)
    !
    offseti8=offset
    !
    call h5sselect_hyperslab_f(dataspace,h5s_select_set_f,offseti8,dim_subset,h5error) 
    !
    if(h5error.ne.0)  stop ' !! error in h5_writesubset: h5sselect_hyperslab_f.'
    !
    ! Create memory dataspace.
    CALL h5screate_simple_f(3,dim_subset,memspace,h5error)
    !
    if(h5error.ne.0)  stop ' !! error in h5_writesubset: Create memory dataspace'
    !
    ! Write subset to dataset 
    ! print*,offseti8,'|',dim_subset
    call h5dwrite_f(dset_id,h5t_native_double,var,dim_subset,h5error,memspace,dataspace)
    !
    if(h5error.ne.0)  stop ' !! error in h5_writesubset: h5dwrite_f 2'
    !
    ! Close everything opened.
    !
    call h5sclose_f(dataspace,h5error)
    if(h5error.ne.0)  stop ' !! error in h5_writesubset: h5sclose_f'
    call h5sclose_f(memspace,h5error)
    if(h5error.ne.0)  stop ' !! error in h5_writesubset: h5sclose_f'
    call h5dclose_f(dset_id,h5error)
    if(h5error.ne.0)  stop ' !! error in h5_writesubset: h5dclose_f'
    call h5fclose_f(file_id,h5error)
    if(h5error.ne.0)  stop ' !! error in h5_writesubset: h5fclose_f'
    !
    ! Close FORTRAN interface.
    call h5close_f(h5error)
    if(h5error.ne.0)  stop ' !! error in h5_writesubset: Close FORTRAN interface'
    !
    write(*,'(A,3(I5,A,I5,A))')'  <<',                                 &
                      offseti8(1),'-',offseti8(1)+dim_subset(1)-1,',', &
                      offseti8(2),'-',offseti8(2)+dim_subset(2)-1,',', &
                      offseti8(3),'-',offseti8(3)+dim_subset(3)-1,'|'
    !
  end subroutine h5_writesubset
  !!
  !+-------------------------------------------------------------------+
  ! this subroutine is used to read a 3d array file with hdf5
  !+-------------------------------------------------------------------+

  subroutine h5_write_blocks(blocks,fname)
    
    use pastr_commtype, only : tblock

    type(tblock),intent(in),target :: blocks(:)
    character(len=*),intent(in) :: fname

    logical :: lfilalive
    integer :: i,n
    type(tblock),pointer :: b
    integer(hid_t) :: file_id,group_id
    ! file identifier
    integer(hid_t) :: dataspace_id,dset_id
    ! dataset identifier
    integer :: h5error ! error flag
    integer(hsize_t) :: dimt(3)
    
    call h5open_f(h5error)
    !
    inquire(file=fname, exist=lfilalive)
    if(lfilalive) then
      call h5fopen_f(fname,H5F_ACC_RDWR_F,file_id,h5error)
    else
      call h5fcreate_f(fname,H5F_ACC_TRUNC_F,file_id,h5error)
    end if
    if(h5error.ne.0)  stop ' !! error in h5_write_blocks 1'
    
    do i=1,size(blocks)

      b=>blocks(i)

      call h5gcreate_f(file_id, b%name, group_id, h5error)

      dimt=(/b%im+1,b%jm+1,b%km+1/)

      do n=1,b%nvar
        CALL h5screate_simple_f(3, dimt, dataspace_id, h5error)
        call h5dcreate_f(group_id, trim(b%varname(n)), h5t_native_double, &
                         dataspace_id, dset_id, h5error)
        call h5dwrite_f(dset_id, h5t_native_double, b%var(:,:,:,n), dimt, h5error)
        CALL h5sclose_f(dataspace_id, h5error)
        call h5dclose_f(dset_id, h5error)
      enddo
      call h5gclose_f(group_id, h5error)

    enddo
    call h5fclose_f(file_id, h5error)

    call h5close_f(h5error)

  end subroutine h5_write_blocks

  subroutine h5_writearray3d(varin,dim1,dim2,dim3,vname,fname,explicit,ierr)
    !
    integer,intent(in) :: dim1,dim2,dim3
    real(8),intent(in) :: varin(0:dim1,0:dim2,0:dim3)
    character(len=*),intent(in) :: vname,fname
    logical,intent(in), optional:: explicit
    integer,intent(out), optional:: ierr
    logical :: lexplicit
    logical :: lfilalive
    !
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
    !
    inquire(file=fname, exist=lfilalive)
    if(lfilalive) then
      call h5fopen_f(fname,H5F_ACC_RDWR_F,file_id,h5error)
    else
      call h5fcreate_f(fname,H5F_ACC_TRUNC_F,file_id,h5error)
    end if
    if(h5error.ne.0)  stop ' !! error in h5_writearray3d 1'
    !
    dimt=(/dim1+1,dim2+1,dim3+1/)
    !
    ! write the dataset.
    call h5ltmake_dataset_f(file_id,vname,3,dimt,h5t_native_double,varin,h5error)
    if(h5error .ne. 0) stop 'h5 write error h5_writearray3d'
    if(lexplicit) print*,' << ',vname,' to ',fname
    if(present(ierr)) ierr=h5error
    if(h5error.ne.0)  stop ' !! error in h5_writearray3d 2'
    !
    call h5fclose_f(file_id, h5error)
    if(h5error.ne.0)  stop ' !! error in h5_writearray3d 3'
    !
    ! close fortran interface.
    call h5close_f(h5error)
    if(h5error.ne.0)  stop ' !! error in h5_writearray3d 4'
    !
  end subroutine h5_writearray3d
  !
  subroutine h5_writearray3d_r4(varin,dim1,dim2,dim3,vname,fname,explicit,ierr)
    !
    integer,intent(in) :: dim1,dim2,dim3
    real(4),intent(in) :: varin(0:dim1,0:dim2,0:dim3)
    character(len=*),intent(in) :: vname,fname
    logical,intent(in), optional:: explicit
    integer,intent(out), optional:: ierr
    logical :: lexplicit
    logical :: lfilalive
    !
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
    !
    inquire(file=fname, exist=lfilalive)
    if(lfilalive) then
      call h5fopen_f(fname,H5F_ACC_RDWR_F,file_id,h5error)
    else
      call h5fcreate_f(fname,H5F_ACC_TRUNC_F,file_id,h5error)
    end if
    if(h5error.ne.0)  stop ' !! error in h5_writearray3d 1'
    !
    dimt=(/dim1+1,dim2+1,dim3+1/)
    !
    ! write the dataset.
    call h5ltmake_dataset_f(file_id,vname,3,dimt,h5t_native_real,varin,h5error)
    if(h5error .ne. 0) stop 'h5 write error h5_writearray3d'
    if(lexplicit) print*,' << ',vname,' to ',fname
    if(present(ierr)) ierr=h5error
    if(h5error.ne.0)  stop ' !! error in h5_writearray3d 2'
    !
    call h5fclose_f(file_id, h5error)
    if(h5error.ne.0)  stop ' !! error in h5_writearray3d 3'
    !
    ! close fortran interface.
    call h5close_f(h5error)
    if(h5error.ne.0)  stop ' !! error in h5_writearray3d 4'
    !
  end subroutine h5_writearray3d_r4
  !
  subroutine h5_writearray3d_int(varin,dim1,dim2,dim3,vname,fname,explicit,ierr)
    !
    integer,intent(in) :: dim1,dim2,dim3
    integer,intent(in) :: varin(0:dim1,0:dim2,0:dim3)
    character(len=*),intent(in) :: vname,fname
    logical,intent(in), optional:: explicit
    integer,intent(out), optional:: ierr
    logical :: lexplicit
    logical :: lfilalive
    !
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
    !
    inquire(file=fname, exist=lfilalive)
    if(lfilalive) then
      call h5fopen_f(fname,H5F_ACC_RDWR_F,file_id,h5error)
    else
      call h5fcreate_f(fname,H5F_ACC_TRUNC_F,file_id,h5error)
    end if
    if(h5error.ne.0)  stop ' !! error in h5_writearray3d 1'
    !
    dimt=(/dim1+1,dim2+1,dim3+1/)
    !
    ! write the dataset.
    call h5ltmake_dataset_f(file_id,vname,3,dimt,h5t_native_integer,varin,h5error)
    if(h5error .ne. 0) stop 'h5 write error h5_writearray3d'
    if(lexplicit) print*,' << ',vname,' to ',fname
    if(present(ierr)) ierr=h5error
    if(h5error.ne.0)  stop ' !! error in h5_writearray3d 2'
    !
    call h5fclose_f(file_id, h5error)
    if(h5error.ne.0)  stop ' !! error in h5_writearray3d 3'
    !
    ! close fortran interface.
    call h5close_f(h5error)
    if(h5error.ne.0)  stop ' !! error in h5_writearray3d 4'
    !
  end subroutine h5_writearray3d_int
  !+-------------------------------------------------------------------+
  !|the end of the subroutine h5_readarray3d                           |
  !+-------------------------------------------------------------------+
  !
  function h5_getdimensio(varname,filenma) result(dims)
    !
    integer(hsize_t) :: dims(3)
    character(len=*),intent(in) :: varname,filenma
    !
    ! local data
    integer(hid_t)  :: file, space, dset
    integer         :: h5error ! error flag
    integer(hsize_t) :: ndims(3)
    !
    !
    call h5open_f(h5error) 
    call h5fopen_f(filenma,h5f_acc_rdonly_f,file,h5error)
    call h5dopen_f (file,varname,dset, h5error)
    call h5dget_space_f(dset, space, h5error)
    call h5sget_simple_extent_dims_f(space,dims,ndims,h5error)
    !
    call h5dclose_f(dset , h5error)
    call h5sclose_f(space, h5error)
    call h5fclose_f(file , h5error)
    call h5close_f(h5error)
    !
  end function h5_getdimensio
  !
  !+-------------------------------------------------------------------+
  !| This function is used to get the dimension of the hdf5 array.     |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 30-JuL-2020 | Coped from ASTR Post by J. Fang STFC Daresbury Lab. |
  !+-------------------------------------------------------------------+
  function h5getdim1d(varname,filenma) result(dims)
    !
    character(len=*),intent(in) :: varname,filenma
    integer :: dims
    !
    ! local data
    integer(hid_t)  :: file, space, dset
    integer(hsize_t) :: ndims(1)
    integer         :: h5error ! error flag
    integer(hsize_t) :: dims_h5(1)
    !
    !
    call h5open_f(h5error)
    call h5fopen_f(filenma,h5f_acc_rdonly_f,file,h5error)
    call h5dopen_f (file,varname,dset, h5error)
    call h5dget_space_f(dset, space, h5error)
    call h5sget_simple_extent_dims_f(space,dims_h5,ndims,h5error)
    !
    dims=dims_h5(1)
    !
    call h5dclose_f(dset , h5error)
    call h5sclose_f(space, h5error)
    call h5fclose_f(file , h5error)
    call h5close_f(h5error)
    !
  end function h5getdim1d
  !
  !+-------------------------------------------------------------------+
  !| This end of the function h5getdim3d.                              |
  !+-------------------------------------------------------------------+
  subroutine h5rcw(filename,dataname,datavalue)
    !
    character(len=*),intent(in) :: filename,dataname
    integer,intent(in) :: datavalue
    !
    integer :: oldvalue
    !
    call H5ReadArray(oldvalue,dataname,filename)
    print*,' ** the value of data',dataname,'is : ',oldvalue
    !
    call h5_write1int(datavalue,dataname,filename,explicit=.true.)
    !
    call H5ReadArray(oldvalue,dataname,filename)
    print*,' ** new value is : ',oldvalue
    !
  end subroutine h5rcw
  !
end module pastr_h5io
!+---------------------------------------------------------------------+
!| The end of the module pastr_h5io.                                   |
!+---------------------------------------------------------------------+