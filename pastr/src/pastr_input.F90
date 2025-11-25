module pastr_input

    use iso_fortran_env, only: wp => real64

    implicit none

    private

    public ::  read_astr_input,read_moniter_input

contains

    !----------------------------------------------------------------------
    !> Read ASTR input file and fill astr_input_type
    !!  Format assumed:
    !!     key = value
    !!  Example:
    !!     nx = 200
    !!     ny = 150
    !!     dt = 0.0005
    !!  Missing keys use defaults
    !----------------------------------------------------------------------
    subroutine read_astr_input()

        use pastr_commvar
        use pastr_io, only: parse_command_line

        character(len=512) :: input_file,line, key, val
        integer :: ios, iunit
        logical :: eof_flag
        integer :: eq_pos
        character(len=6) :: iden
        logical :: file_exist

        ! Initialize with defaults

        call parse_command_line(string=iden)
        if(iden=='-input') then
          call parse_command_line(string=input_file)
        endif

        inquire(file=trim(input_file), exist=file_exist)
        if (.not. file_exist) then
            write(*,*) "Error: ASTR input file '", trim(input_file), "' does not exist!"
            stop 1
        end if

        ! Open file
        open(newunit=iunit, file=trim(input_file), status='old', action='read', iostat=ios)
        read(iunit,'(///////)')
        read(iunit,*)im,jm,km
        read(iunit,"(/)")
        read(iunit,*)lihomo,ljhomo,lkhomo
        read(iunit,"(/)")
        read(iunit,*)nondimen
        read(iunit,'(///////)')
        if(nondimen) then
          read(iunit,*)ref_t,reynolds,mach
        else
          read(iunit,*)ref_t,ref_vel,ref_len,ref_den
        endif
        read(iunit,'(///////////////////////////)')
        read(iunit,'(A)')gridfile
        close(iunit)
        write(*,*)' >> ',trim(input_file)

        write(*,*)' ==========================informout=========================='
        write(*,*)'                    im                  jm                  km'
        write(*,"(3X,3(I20))")im,jm,km
        write(*,*)' -------------------------------------------------------------'
        write(*,*)'              Reynolds                Mach               ref_t'
        write(*,"(3X,3(F20.5))")Reynolds,mach,ref_t
        write(*,*)' -------------------------------------------------------------'
        write(*,"(3X,A40,A20)")'grid: ',trim(gridfile)
        write(*,*)' ==========================informout=========================='

    end subroutine read_astr_input

    subroutine read_moniter_input(imon,jmon,kmon)

      use pastr_commvar, only: im,jm,km
      use pastr_io, only: read_grid
      use pastr_tecio, only: tecbin
      
      integer,allocatable :: imon(:),jmon(:),kmon(:)
      real(wp),allocatable :: x(:,:),y(:,:)

      integer :: ios,nmonitor,i,n

      open(12,file='datin/monitor.dat')
      read(12,'()')
      ios=0
      nmonitor=0
      do while(ios==0)
        read(12,*,iostat=ios)n
        !
        if(ios==0) nmonitor=nmonitor+1
        !
      enddo
      close(12)
      print*,' ** number of monitor points',nmonitor

      allocate(imon(nmonitor),jmon(nmonitor),kmon(nmonitor))
      open(12,file='datin/monitor.dat')
      read(12,'()')
      ios=0
      do i=1,nmonitor
        read(12,*)imon(i),jmon(i),kmon(i)
      enddo
      close(12)

      allocate(x(0:im,0:jm),y(0:im,0:jm))
      call read_grid(x=x,y=y,kslice=kmon(1))
      do i=1,nmonitor
        print*,' ** montir',i,':',x(imon(i),jmon(i)),y(imon(i),jmon(i))
      enddo

      call tecbin('tecgrid.plt',x,'x',y,'y')
      
      open(18,file='monitor_point.dat')
      write(18,*)'TITLE     = "scalar field for tecplot"'
      write(18,*)'VARIABLES = "x" "y" "n"'
      write(18,*)'ZONE T="ZONE 001"'
      write(18,*)'I=',nmonitor,', J=1, K=1, ZONETYPE=Ordered'
      write(18,*)'DATAPACKING=POINT'
      do i=1,nmonitor
        write(18,*)x(imon(i),jmon(i)),y(imon(i),jmon(i)),i
      enddo
      close(18)
      print*,' << monitor_point.dat'

      deallocate(x,y)

    end subroutine read_moniter_input

end module pastr_input
