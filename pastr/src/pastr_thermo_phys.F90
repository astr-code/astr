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

    public :: sos
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
pure function viscosity_suth(T) result(mu)
    real(8), intent(in) :: T
    real(8) :: mu

    mu = mu0 * (T / T0_suth)**1.5d0 * (T0_suth + Suth) / (T + Suth)
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
pure function thermal_conductivity(T) result(k)
    real(8), intent(in) :: T
    real(8) :: k

    real(8) :: muT

    muT = viscosity_suth(T)
    k   = muT * cp() / Pr
end function thermal_conductivity

end module pastr_thermo_phys
