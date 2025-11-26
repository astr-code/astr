module pastr_commvar

    use iso_fortran_env, only: wp => real64

    implicit none

    integer :: im, jm, km

    logical :: nondimen=.true.
    logical :: lihomo, ljhomo, lkhomo

    character(len=128) :: gridfile

    real(wp) :: ref_t, reynolds, mach, ref_vel, ref_len, ref_den
    real(wp) :: prandtl, gamma, rgas, ref_miu, tempconst, tempconst1, const1, &
                const2, const3, const4, const5, const6, const7, const8, roinf, &
                uinf, vinf, tinf, pinf
    real(wp) :: deltat    

    real(wp), allocatable, dimension(:, :, :) :: x, y, z, ro, u1, u2, u3, p, t
    real(wp), allocatable, dimension(:, :, :) :: ro_m, u1_m, u2_m, u3_m, p_m, t_m
    real(wp), allocatable, dimension(:, :)   :: ro_zm, u1_zm, u2_zm, u3_zm, &
                                               p_zm, t_zm
    real(wp), allocatable, dimension(:)     :: ro_xzm, u1_xzm, u2_xzm, u3_xzm, &
                                              p_xzm, t_xzm
    real(wp), allocatable, dimension(:, :, :, :) :: du1_m, du2_m, du3_m, dt_m
    real(wp), allocatable, dimension(:, :, :) :: du1_zm, du2_zm, du3_zm, dt_zm
    real(wp), allocatable, dimension(:, :) :: du1_xzm, du2_xzm, du3_xzm, dt_xzm

    type :: montype
      integer :: npoints,nvariables
      character(len=16),allocatable :: varname(:)
      integer,allocatable :: nstep(:)
      real(wp),allocatable :: time(:),data(:,:)
      contains
      procedure :: init => alloc_monitor
    end type montype

    type :: conntype

      integer :: nblock_local,npatch_local
      integer :: nblock_remot,npatch_remot

    end type conntype

    type :: pachtype
      
      integer :: imin,imax,jmin,jmax,kmin,kmax

      character(len=1) :: ntyp ! internal (i) or bc (n)
      
      type(conntype),allocatable :: connector(:)

      character(len=6),allocatable :: varname(:)

      real(wp),allocatable :: x(:,:,:,:),var(:,:,:,:)


    end type pachtype

    type :: bloktype

      integer :: im,jm,km
      integer :: nvar,nhalo,npatch

      character(len=6) :: name

      character(len=6),allocatable :: varname(:)

      real(wp),allocatable :: x(:,:,:,:),var(:,:,:,:)

      type(pachtype),allocatable :: patch(:)

      contains

      procedure :: init => alloc_block

    end type bloktype


contains

    subroutine alloc_monitor(amonitor)
    !
    class(montype),target :: amonitor

    allocate(amonitor%varname(amonitor%nvariables))
    allocate(amonitor%time(amonitor%npoints))
    allocate(amonitor%nstep(amonitor%npoints))
    allocate(amonitor%data(amonitor%nvariables,amonitor%npoints))

  end subroutine alloc_monitor

    subroutine alloc_block(ablock)

    class(bloktype),target :: ablock

    ablock%nhalo=4

    allocate(ablock%varname(ablock%nvar))
    allocate(ablock%x( -ablock%nhalo:ablock%im+ablock%nhalo, &
                       -ablock%nhalo:ablock%jm+ablock%nhalo, &
                       -ablock%nhalo:ablock%km+ablock%nhalo,3) )
    allocate(ablock%var(-ablock%nhalo:ablock%im+ablock%nhalo, &
                        -ablock%nhalo:ablock%jm+ablock%nhalo, &
                        -ablock%nhalo:ablock%km+ablock%nhalo, ablock%nvar))

    allocate(ablock%patch(ablock%npatch))

    print*,' ** block ',ablock%name,' initiated.'

  end subroutine alloc_block

end module pastr_commvar