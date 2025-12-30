module pastr_io

    use iso_fortran_env, only: wp => real64

    implicit none

    private

    public :: parse_command_line,read_grid
    public :: global_n_arguments,read_monitor_data,write_monitor_data
    public :: read_stats

    integer :: global_n_arguments=0

    Interface read_grid
      module procedure read_grid_2d
      module procedure read_grid_3d
    end Interface read_grid
    !
contains

    !------------------------------------------------------------------
    !> Automatically parse command-line arguments into strings and integers
    !!  - Arguments that can be read as integers go to `numbers`
    !!  - All other arguments go to `strings`
    !------------------------------------------------------------------
    subroutine parse_command_line(string, inumber,rnumber)
        
        use pastr_utility

        character(len=*), intent(inout), optional :: string
        integer, intent(inout), optional :: inumber
        real(wp), intent(inout), optional :: rnumber

        integer :: n_args, i, n_count, s_count
        character(len=255) :: arg
        integer :: num
        logical :: proper_assign

        proper_assign=.false.

        global_n_arguments=global_n_arguments+1

        call get_cmd_argument(global_n_arguments,arg)
        
        if(present(string)) string=''
        if(present(inumber)) inumber=-1
        if(present(rnumber)) rnumber=0.d0

        if(isnum(arg)==0) then
          ! non-numeric
          if(present(string)) then
            string=trim(arg)
            proper_assign=.true.
          endif

        elseif(isnum(arg)==1) then
          ! integer
          if(present(inumber)) then
            read(arg,*)inumber
            proper_assign=.true.
          endif

        elseif(isnum(arg)>1) then
          ! real number
          if(present(rnumber)) then
            read(arg,*)rnumber
            proper_assign=.true.
          endif
        else
          stop ' !!  error1 @ parse_command_line !!'
        endif

        if(.not.  proper_assign) global_n_arguments=global_n_arguments-1
        

    end subroutine parse_command_line

    subroutine read_grid_2d(x,y,z,islice,jslice,kslice)
      
      use pastr_commvar, only: gridfile,im,jm,km
      use pastr_h5io

      real(wp),intent(inout),allocatable,optional :: x(:,:),y(:,:),z(:,:)
      integer,intent(in),optional :: islice,jslice,kslice

      if(present(islice)) then
        if(present(x)) call h5_read2dfrom3d(x,im,jm,km,'x',trim(gridfile),islice=islice)
        if(present(y)) call h5_read2dfrom3d(y,im,jm,km,'y',trim(gridfile),islice=islice)
        if(present(z)) call h5_read2dfrom3d(z,im,jm,km,'z',trim(gridfile),islice=islice)
      elseif(present(jslice)) then
        if(present(x)) call h5_read2dfrom3d(x,im,jm,km,'x',trim(gridfile),jslice=jslice)
        if(present(y)) call h5_read2dfrom3d(y,im,jm,km,'y',trim(gridfile),jslice=jslice)
        if(present(z)) call h5_read2dfrom3d(z,im,jm,km,'z',trim(gridfile),jslice=jslice)
      elseif(present(kslice)) then
        if(present(x)) call h5_read2dfrom3d(x,im,jm,km,'x',trim(gridfile),kslice=kslice)
        if(present(y)) call h5_read2dfrom3d(y,im,jm,km,'y',trim(gridfile),kslice=kslice)
        if(present(z)) call h5_read2dfrom3d(z,im,jm,km,'z',trim(gridfile),kslice=kslice)
      endif

    end subroutine read_grid_2d

    subroutine read_grid_3d(x,y,z)
      
      use pastr_commvar, only: gridfile,im,jm,km,lx,ly,lz
      use pastr_h5io

      real(wp),intent(inout),allocatable :: x(:,:,:),y(:,:,:),z(:,:,:)

      call H5ReadArray(x,im,jm,km,'x',trim(gridfile))
      call H5ReadArray(y,im,jm,km,'y',trim(gridfile))
      call H5ReadArray(z,im,jm,km,'z',trim(gridfile))

      lx=x(im,0,0)-x(0,0,0)
      ly=y(0,jm,0)-y(0,0,0)
      lz=z(0,0,km)-z(0,0,0)

    end subroutine read_grid_3d

    subroutine read_monitor_data(num_first_file,num_last_file,n_data,n_col,mon_data)

      use pastr_commvar, only: montype

      integer,intent(in) :: num_first_file,num_last_file,n_data,n_col
      type(montype),allocatable,intent(out) :: mon_data(:)
      integer :: recl_size,stat,nfiles
      
      character(len=4) :: mfname
      character(len=16) :: varname(n_data-1)
      integer :: i,j,n
      
      open(12,file='mondef.txt')
      do i=1,n_data-1
        read(12,*)varname(i)
      enddo
      close(12)
      print*,' >> mondef.txt'
      print*,' ** variables recorded in monitor:',(trim(varname(i)),';',i=1,n_data-1)

      recl_size=8*(n_data+1)
      nfiles=num_last_file-num_first_file+1
      allocate(mon_data(nfiles))

      do n=1,nfiles

        mon_data(n)%npoints=n_col
        mon_data(n)%nvariables=n_data-1
        call mon_data(n)%init()
        
        mon_data(n)%varname=varname

      enddo

      do j=num_first_file,num_last_file

        write(mfname,'(I4.4)')j
        
        open(12,file='monitor/monitor'//mfname//'.dat',access='direct',recl=recl_size,action='read')
        do i=1,n_col
          read(12,rec=i,iostat=stat)mon_data(j)%nstep(i),mon_data(j)%time(i),mon_data(j)%data(:,i)
        enddo
        close(12)
        print*,' >> monitor/monitor'//mfname//'.dat'

      enddo

    end subroutine read_monitor_data

    subroutine write_monitor_data(mon_data)
      
      use pastr_commvar, only: montype

      type(montype),intent(in) :: mon_data(:)

      integer :: msize,i,j,k
      character(len=4) :: mfname
      character(len=120) :: txtformat,datformat

      msize=size(mon_data)

      write(txtformat,'(A,I0,A)')'(',mon_data(1)%nvariables+1,'(1X,A20))'
      write(datformat,'(A,I0,A)')'(',mon_data(1)%nvariables+1,'(1X,E20.13E2))'


      do j=1,msize
        write(mfname,'(I4.4)')j

        open(18,file='monitor'//mfname//'.txt')
        write(18,txtformat)'time',(trim(mon_data(j)%varname(k)),k=1,mon_data(j)%nvariables)
        do i=1,mon_data(j)%npoints
          write(18,datformat)mon_data(j)%time(i),mon_data(j)%data(:,i)
        enddo
        close(18)
        print*,' << monitor'//mfname//'.txt'

      enddo

    end subroutine write_monitor_data

    subroutine read_stats(var,varname)

      use pastr_commvar, only: im,jm,km
      use pastr_h5io

      real(wp) :: var(0:im,0:jm,0:km)
      character(len=*),intent(in) :: varname

      character(len=18) :: fname
      integer,save :: nsamples=0

      fname='outdat/meanflow.h5'

      if(nsamples==0) then
        call H5ReadArray(nsamples,'nsamples',fname)
        print*,' ** nsamples: ',nsamples
      endif

      call H5ReadArray(var,im,jm,km,varname,fname)
      var = var/dble(nsamples)

    end subroutine read_stats

end module pastr_io
