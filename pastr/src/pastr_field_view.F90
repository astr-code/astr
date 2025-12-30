module pastr_field_view

    use iso_fortran_env, only: wp => real64

    implicit none


    private

    public :: write_xy_slice,write_xz_slice,write_3d_field

contains

    subroutine write_xy_slice(filein,fileout,nfirst,nlast,slice,format)

      use pastr_io, only : read_grid,parse_command_line
      use pastr_multiblock_type, only: block_define
      use pastr_commvar, only : im,jm,km
      use pastr_h5io
      use pastr_tecio
      use pastr_xdmf
      use pastr_data_convert
      use pastr_commtype

      character(len=*),intent(in), optional :: filein,fileout
      integer,intent(in), optional :: nfirst,nlast
      character(len=*),intent(in) :: format
      integer,intent(in) :: slice

      real(wp),allocatable :: x(:,:),y(:,:),z(:,:)
      integer :: num_var_in,num_var_out,nsta,nend,n,i,j,is,ie,js,je
      character(len=6),allocatable :: nam_var_in(:),nam_var_out(:)
      character(len=128) :: file2read,file2write
      character(len=5) :: fname
      character(len=6) :: src_of_data
      real(wp),allocatable :: r8dat(:,:)
      real(wp),allocatable,dimension(:,:,:) :: dat_var_in,dat_var_out
      type(tblock),allocatable,target :: pblocks(:)
      type(tblock),pointer :: b
      integer :: nblocks
      logical :: multi_block

      call parse_command_line( string=src_of_data )

      if(len(trim(src_of_data))<1) src_of_data='outdat'

      print*,' ** source of data: ',src_of_data

      allocate(x(0:im,0:jm),y(0:im,0:jm),z(0:im,0:jm))

      call read_grid(x=x,y=y,z=z,kslice=slice)

      call var_define(num_var_in,nam_var_in,num_var_out,nam_var_out)

      call range_define(is,ie,js,je)

      allocate(dat_var_in(0:im,0:jm,1:num_var_in),dat_var_out(0:im,0:jm,1:num_var_out))

      if(present(filein) .and. present(fileout)) then
        nsta=-1
        nend=-1
      elseif(present(nfirst) .and. present(nlast)) then
        nsta=nfirst
        nend=nlast
      endif

      allocate(r8dat(0:im,0:jm))

      do n=nsta,nend

        if(n<0) then
          file2read =filein
          file2write=fileout
        else
          if(src_of_data=='outdat') then
            write(fname,'(i4.4)') n
            file2read ='outdat/flowfield'//trim(fname)//'.h5'
          elseif(src_of_data=='kslice') then
            write(fname,'(i5.5)') n
            file2read ='kslice/kslice'//trim(fname)//'.h5'
          endif
          
          file2write='visuxy'//trim(fname)
          if(trim(format)=='plt') file2write=trim(file2write)//'.plt'

        endif

        do i=1,num_var_in
          if(src_of_data=='outdat') then
            call h5_read2dfrom3d(r8dat,im,jm,km,trim(nam_var_in(i)),trim(file2read),kslice=slice)
          elseif(src_of_data=='kslice') then
            call H5ReadArray(r8dat,im,jm,trim(nam_var_in(i)),trim(file2read))
          endif
          dat_var_in(:,:,i)=real(r8dat)
        enddo

        dat_var_out=dat_out_cal(dat_var_in,nam_var_in,nam_var_out,x,y)

        call block_define(multi_block,nblocks,pblocks)

        do i=1,nblocks
          b=>pblocks(i)
          b%nvar=num_var_out
          call b%init_data()
          b%x(0:b%im,0:b%jm,0)=x(b%ilo:b%ihi,b%jlo:b%jhi)
          b%y(0:b%im,0:b%jm,0)=y(b%ilo:b%ihi,b%jlo:b%jhi)
          b%var(0:b%im,0:b%jm,0,1:num_var_out)=dat_var_out(b%ilo:b%ihi,b%jlo:b%jhi,1:num_var_out )
          b%varname=nam_var_out
        enddo

        if(multi_block) then
            call xdmfwriter(dir='snapshot/',filename=trim(file2write),blocks=pblocks )
        else
          if(trim(format)=='plt') then
            call writetecbin2dlist_xy('snapshot/'//trim(file2write),x(is:ie,js:je), &
                                                                    y(is:ie,js:je), &
                                                                    dat_var_out(is:ie,js:je,:),nam_var_out)
          elseif(trim(format)=='xdmf') then
            call xdmfwriter(dir='snapshot/',filename=trim(file2write),x=x(is:ie,0),y=y(0,js:je), &
                                               var=dat_var_out(is:ie,js:je,:),varname=nam_var_out )
          endif
        endif

      enddo

    end subroutine write_xy_slice

    subroutine write_xz_slice(filein,fileout,nfirst,nlast,slice,format)

      use pastr_io, only : read_grid,parse_command_line
      use pastr_commvar, only : im,jm,km
      use pastr_h5io
      use pastr_tecio
      use pastr_xdmf
      use pastr_data_convert

      character(len=*),intent(in), optional :: filein,fileout
      integer,intent(in), optional :: nfirst,nlast
      character(len=*),intent(in) :: format
      integer,intent(in) :: slice

      real(wp),allocatable :: x(:,:),y(:,:),z(:,:)
      integer :: num_var_in,num_var_out,nsta,nend,n,i,j,is,ie,js,je,ks,ke
      character(len=6),allocatable :: nam_var_in(:),nam_var_out(:)
      character(len=128) :: file2read,file2write
      character(len=5) :: fname
      character(len=6) :: src_of_data
      real(wp),allocatable :: r8dat(:,:)
      real(wp),allocatable,dimension(:,:,:) :: dat_var_in,dat_var_out

      call parse_command_line( string=src_of_data )

      if(len(trim(src_of_data))<1) src_of_data='outdat'

      print*,' ** source of data: ',src_of_data

      allocate(x(0:im,0:km),y(0:im,0:km),z(0:im,0:km))

      call read_grid(x=x,y=y,z=z,jslice=slice)

      call var_define(num_var_in,nam_var_in,num_var_out,nam_var_out)

      call range_define(is,ie,js,je,ks,ke)

      allocate(dat_var_in(0:im,0:km,1:num_var_in),dat_var_out(0:im,0:km,1:num_var_out))

      if(present(filein) .and. present(fileout)) then
        nsta=-1
        nend=-1
      elseif(present(nfirst) .and. present(nlast)) then
        nsta=nfirst
        nend=nlast
      endif

      allocate(r8dat(0:im,0:km))

      do n=nsta,nend

        if(n<0) then
          file2read =filein
          file2write=fileout
        else
          if(src_of_data=='outdat') then
            write(fname,'(i4.4)') n
            file2read ='outdat/flowfield'//trim(fname)//'.h5'
          elseif(src_of_data=='slice') then
            write(fname,'(i5.5)') n
            file2read ='slice/slice'//trim(fname)//'.h5'
          endif
          
          file2write='visuxz'//trim(fname)
          if(trim(format)=='plt') file2write=trim(file2write)//'.plt'

        endif

        do i=1,num_var_in
          if(src_of_data=='outdat') then
            call h5_read2dfrom3d(r8dat,im,jm,km,trim(nam_var_in(i)),trim(file2read),jslice=slice)
          elseif(src_of_data=='jslice') then
            call H5ReadArray(r8dat,im,km,trim(nam_var_in(i)),trim(file2read))
          endif
          dat_var_in(:,:,i)=real(r8dat)
        enddo

        dat_var_out=dat_out_cal(dat_var_in,nam_var_in,nam_var_out,x,y)

        if(trim(format)=='plt') then
          call writetecbin2dlist_xy('snapshot/'//trim(file2write),x(is:ie,ks:ke), &
                                                                  z(is:ie,ks:ke), &
                                                        dat_var_out(is:ie,ks:ke,:),nam_var_out)
        elseif(trim(format)=='xdmf') then
          call xdmfwriter(dir='snapshot/',filename=trim(file2write),x=x(is:ie,0),y=z(0,ks:ke), &
                                             var=dat_var_out(is:ie,ks:ke,:),varname=nam_var_out )
        endif

      enddo

    end subroutine write_xz_slice

    subroutine write_3d_field(filein,fileout,nfirst,nlast,format)

      use pastr_io, only : read_grid,parse_command_line
      use pastr_commvar, only : im,jm,km
      use pastr_h5io
      use pastr_tecio
      use pastr_xdmf
      use pastr_data_convert

      character(len=*),intent(in), optional :: filein,fileout
      integer,intent(in), optional :: nfirst,nlast
      character(len=*),intent(in) :: format

      real(wp),allocatable :: x(:,:,:),y(:,:,:),z(:,:,:)
      integer :: num_var_in,num_var_out,nsta,nend,n,m,i,j,k,is,ie,js,je,ks,ke
      character(len=6),allocatable :: nam_var_in(:),nam_var_out(:)
      character(len=128) :: file2read,file2write
      character(len=5) :: fname
      character(len=6) :: src_of_data
      real(wp),allocatable :: r8dat(:,:,:)
      real(wp),allocatable,dimension(:,:,:,:) :: dat_var_in,dat_var_out

      src_of_data='outdat'

      allocate(x(0:im,0:jm,0:km),y(0:im,0:jm,0:km),z(0:im,0:jm,0:km))

      call read_grid(x=x,y=y,z=z)

      call var_define(num_var_in,nam_var_in,num_var_out,nam_var_out)

      call range_define(is,ie,js,je,ks,ke)

      allocate(dat_var_in(0:im,0:jm,0:km,1:num_var_in),dat_var_out(0:im,0:jm,0:km,1:num_var_out))

      if(present(filein) .and. present(fileout)) then
        nsta=-1
        nend=-1
      elseif(present(nfirst) .and. present(nlast)) then
        nsta=nfirst
        nend=nlast
      endif

      allocate(r8dat(0:im,0:jm,0:km))

      do n=nsta,nend

        if(n<0) then
          file2read =filein
          file2write=fileout
        else

          write(fname,'(i4.4)') n
          file2read ='outdat/flowfield'//trim(fname)//'.h5'
          
          file2write='visu3d'//trim(fname)
          if(trim(format)=='plt') file2write=trim(file2write)//'.plt'

        endif

        do m=1,num_var_in

          call H5ReadArray(r8dat,im,jm,km,trim(nam_var_in(m)),trim(file2read))

          dat_var_in(:,:,:,m)=r8dat

        enddo

        dat_var_out=dat_out_cal(dat_var_in,nam_var_in,nam_var_out,x,y,z)
        
        if(trim(format)=='plt') then
          call writetecbin3dlist('snapshot/'//trim(file2write),x(is:ie,js:je,ks:ke), &
                                                               y(is:ie,js:je,ks:ke), &
                                                               z(is:ie,js:je,ks:ke), &
                                                     dat_var_out(is:ie,js:je,ks:ke,:),nam_var_out)
        elseif(trim(format)=='xdmf') then
          call xdmfwriter(dir='snapshot/',filename=trim(file2write),x=x(is:ie,0,0),y=y(0,js:je,0),z=z(0,0,ks:ke), &
                                               var=dat_var_out(is:ie,js:je,ks:ke,:),varname=nam_var_out )
        endif

      enddo

    end subroutine write_3d_field

    subroutine range_define(is,ie,js,je,ks,ke)

      use pastr_commvar, only : im,jm,km

      integer,intent(out),optional :: is,ie,js,je,ks,ke

      logical :: lfex

      inquire(file='randef.txt',exist=lfex)
      if(lfex) then
        open(14,file='randef.txt')
        if( present(is) .and. present(ie) .and. &
            present(js) .and. present(je) .and. &
            present(ks) .and. present(ke)) then
          read(14,*)is,ie,js,je,ks,ke
          write(*,'(6(A,I0),A)')'  **        output range: (',is,'-',ie,')x(',js,'-',je,')x(',ks,'-',ke,')'
        elseif(present(is) .and. present(ie) .and. &
               present(js) .and. present(je)) then
          read(14,*)is,ie,js,je
          write(*,'(4(A,I0),A)')'  **        output range: (',is,'-',ie,')x(',js,'-',je,')'
        elseif(present(is) .and. present(ie) ) then
          read(14,*)is,ie
          write(*,'(2(A,I0),A)')'  **        output range: (',is,'-',ie,')'
        endif
        close(14)
        print*,' >> randef.txt'
      else
        if( present(is) )is=0
        if( present(ie) )ie=im
        if( present(js) )js=0
        if( present(je) )je=jm
        if( present(ks) )ks=0
        if( present(ke) )ke=km
      endif

    end subroutine range_define

    subroutine var_define(num_var_in,nam_var_in,num_var_out,nam_var_out)

      integer,intent(out) :: num_var_in,num_var_out
      character(len=6),allocatable,intent(out) :: nam_var_in(:),nam_var_out(:)

      integer :: i
      logical :: lfex

      inquire(file='vardef.txt',exist=lfex)
      if(lfex) then
        open(14,file='vardef.txt')
        read(14,*)num_var_in,num_var_out
        close(14)
        num_var_out=num_var_out
      else
        num_var_in=6
        num_var_out=6
      endif
      print*,' **   number of inputs :',num_var_in
      print*,' **   number of outputs:',num_var_out

      allocate(nam_var_in(num_var_in),nam_var_out(num_var_out))

      if(lfex) then
        open(14,file='vardef.txt')
        read(14,*)
        do i=1,num_var_in
          read(14,*)nam_var_in(i)
        enddo
        read(14,*)
        do i=1,num_var_out
          read(14,*)nam_var_out(i)
        enddo
        close(14)
        print*,' >> vardef.txt'
      else
        nam_var_in(1)='ro'
        nam_var_in(2)='u1'
        nam_var_in(3)='u2'
        nam_var_in(4)='u3'
        nam_var_in(5)= 'p'
        nam_var_in(6)= 't'

        nam_var_out=nam_var_in
      endif

      print*,' **     variables read : ',nam_var_in
      print*,' **     variables write: ',nam_var_out
    
    end subroutine var_define

end module pastr_field_view