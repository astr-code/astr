module pastr_commtype

    use iso_fortran_env, only: wp => real64
    use pastr_constdef

    implicit none

    type :: montype
      integer :: npoints,nvariables
      character(len=16),allocatable :: varname(:)
      integer,allocatable :: nstep(:)
      real(wp),allocatable :: time(:),data(:,:)
      contains
      procedure :: init => alloc_monitor
    end type montype

    type :: tblock
      character(len=6) :: name
      integer :: im,jm,km,nvar
      integer :: ilo,ihi,jlo,jhi,klo,khi
      character(len=16),allocatable :: varname(:)
      real(wp), allocatable, dimension(:,:,:) :: x,y,z
      real(wp), allocatable, dimension(:,:,:,:) :: var
      contains
      procedure :: init_data => alloc_block_data
    end type tblock

contains

    subroutine alloc_monitor(amonitor)

      class(montype),target :: amonitor

      allocate(amonitor%varname(amonitor%nvariables))
      allocate(amonitor%time(amonitor%npoints))
      allocate(amonitor%nstep(amonitor%npoints))
      allocate(amonitor%data(amonitor%nvariables,amonitor%npoints))

    end subroutine alloc_monitor

    subroutine alloc_block_data(ablock)

      class(tblock),target :: ablock

      allocate(ablock%varname(ablock%nvar))
      allocate(ablock%x(0:ablock%im,0:ablock%jm,0:ablock%km))
      allocate(ablock%y(0:ablock%im,0:ablock%jm,0:ablock%km))
      allocate(ablock%var(0:ablock%im,0:ablock%jm,0:ablock%km,1:ablock%nvar))

    end subroutine alloc_block_data

end module pastr_commtype