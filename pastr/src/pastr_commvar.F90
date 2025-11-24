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

end module pastr_commvar