module pastr_thermo_phys
    
    use iso_fortran_env, only: wp => real64
    
    !!  Basic thermophysical property functions:
    !!     - speed of sound
    !!     - temperature from ideal gas EOS
    !!     - viscosity (Sutherland)
    !!     - thermal conductivity
    !!     - specific heats cp, cv
    !!
    implicit none
    private

    public :: sos,upluscal,viscosity_suth
    ! public :: viscosity_suth
    ! public :: thermal_conductivity
    ! public :: cp, cv
    ! public :: temperature_from_rho_p

    ! ===== Default air properties =====
    real(8), parameter :: rgas    = 287.1d0         ! J/kg/K
    real(8), parameter :: gamma_g = 1.4d0

    ! Sutherland constants for air
    real(8), parameter :: mu0     = 1.716e-5        ! reference viscosity [Pa·s]
    real(8), parameter :: T0_suth = 273.15d0        ! reference temperature [K]
    real(8), parameter :: Suth    = 110.4d0         ! Sutherland temperature [K]

    ! Prandtl number (approx constant)
    real(8), parameter :: Pr      = 0.72d0

contains
    !==============================================================
    ! SPEED OF SOUND (ideal gas)
    !    a = sqrt(γ R T)
    !==============================================================

    real(wp) function sos(tmp,spc)
      !
      use pastr_commvar,   only: nondimen,mach,gamma
      !
      real(wp),intent(in) :: tmp
      real(wp),intent(in),optional :: spc(:)
      !
      ! local data
      real(wp) :: cpcmix,gamrgc
      !
      if(nondimen) then
        sos=sqrt(tmp)/mach
      else
        sos=sqrt(gamma*rgas*tmp)
      endif
      
      return
      !
    end function sos

!==============================================================
! VISCOSITY — Sutherland Law
!    μ(T) = μ0 * (T/T0)**3/2 * (T0 + S) / (T + S)
!==============================================================
function viscosity_suth(T,nondim,ref_temperature) result(mu)

    real(wp), intent(in) :: T
    logical,intent(in),optional :: nondim
    real(wp), intent(in),optional :: ref_temperature
    real(wp) :: mu

    logical :: lnond
    real(wp) :: tnondim,tempconst,tempconst1

    if(present(nondim)) then
      lnond=nondim
    else
      lnond=.true.
    endif

    if(lnond) then

      if(.not. present(ref_temperature)) stop ' reference temperature requred in viscosity_suth'

      tempconst=Suth/ref_temperature
      tempconst1=1._wp+tempconst
      mu=T*sqrt(T)*tempconst1/(T+tempconst)

    else

      tnondim=T/T0_suth
      mu = mu0 * tnondim*sqrt(tnondim) * (T0_suth + Suth) / (T + Suth)

    endif

    return

end function viscosity_suth

!==============================================================
! SPECIFIC HEATS FOR IDEAL AIR (constant cp, cv)
!==============================================================
pure function cp() result(cpval)
    real(8) :: cpval
    cpval = gamma_g * Rgas / (gamma_g - 1.0d0)
end function cp

pure function cv() result(cvval)
    real(8) :: cvval
    cvval = Rgas / (gamma_g - 1.0d0)
end function cv

!==============================================================
! TEMPERATURE FROM RHO AND PRESSURE (ideal gas)
!    p = rho * R * T
!==============================================================
pure function temperature_from_rho_p(rho, p) result(T)
    real(8), intent(in) :: rho, p
    real(8) :: T
    T = p / (rho * Rgas)
end function temperature_from_rho_p

!==============================================================
! THERMAL CONDUCTIVITY
!    k = μ cp / Pr
!==============================================================
function thermal_conductivity(T) result(k)
    real(8), intent(in) :: T
    real(8) :: k

    real(8) :: muT

    muT = viscosity_suth(T)
    k   = muT * cp() / Pr
end function thermal_conductivity

!+-------------------------------------------------------------------+
  !| this subroutine is to calculate uplus-yplse                       |
  !+-------------------------------------------------------------------+
  subroutine upluscal(uplus,yplus,u,y,ro,t,utau)
    !
    use pastr_gradients,only: grad_y
    !
    ! arguments
    real(wp),intent(out),allocatable :: uplus(:),yplus(:)
    real(wp),intent(in) :: u(0:),y(0:),ro(0:),t(0:)
    real(wp),intent(out),optional :: utau
    !
    ! local data
    integer :: dim,j,j1
    real(wp),allocatable :: dudy(:),uvd(:),yin(:)
    real(wp) :: miu,utaw_local,var1,var2,tawx
    
    dim=size(u)-1

    allocate(dudy(0:dim),uvd(0:dim),yin(0:dim))

    yin=y-y(0)

    miu=viscosity_suth(t(0))

    dudy=grad_y(u,yin)

    tawx=abs(miu*dudy(0))

    utaw_local=sqrt(tawx/ro(0))

    if(present(utau)) utau=utaw_local

    allocate(uplus(0:dim),yplus(0:dim))

    uvd(0)=0._wp
    do j=1,dim
      uvd(j)=0._wp
      do j1=1,j
        var1=sqrt(0.5_wp*(ro(j1)+ro(j1-1))/ro(0))
        var2=u(j1)-u(j1-1)
        uvd(j)=uvd(j)+var1*var2
      end do
    end do
    !
    do j=0,dim
      yplus(j)=ro(0)*utaw_local*yin(j)/miu
      uplus(j)=uvd(j)/utaw_local
    end do
    !
    if(allocated(dudy)) deallocate(dudy)
    
    deallocate(uvd)
    !
    write(*,'(A)')'  ------ first point from wall ------'
    write(*,'(A)')'                y1+               u1+'
    write(*,'(1x,F18.7,F18.7)')yplus(1),uplus(1)
    write(*,'(A)')'  -----------------------------------'
    print*,' ** tawx=',tawx
    print*,' ** utaw=',utaw_local
    print*,' ** uplus-yplus calculated.'
    !
  end subroutine upluscal
  !+-------------------------------------------------------------------+
  ! End of subroutine upluscal.                                        |
  !+-------------------------------------------------------------------+

end module pastr_thermo_phys
