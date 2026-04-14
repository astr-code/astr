module commarray

   implicit none

   real(8) :: dx, dy, dz

end module commarray

module commvar

   implicit none

   integer :: im, jm, km, hm
   parameter(hm=5)

end module commvar

program compare_outputs
   use commvar
   use commarray
   use vtkio
   implicit none
   real(8), allocatable, dimension(:, :, :, :) :: vel_1
   real(8), allocatable, dimension(:, :, :) :: rho_1, prs_1, tmp_1
   real(8), allocatable, dimension(:, :, :, :) :: vel_2
   real(8), allocatable, dimension(:, :, :) :: rho_2, prs_2, tmp_2

   real(8), allocatable, dimension(:, :, :, :) :: vel_err
   real(8), allocatable, dimension(:, :, :) :: rho_err, prs_err, tmp_err
   real(8) :: max_rho, min_rho
   real(8) :: max_prs, min_prs
   real(8) :: max_tmp, min_tmp
   real(8) :: max_vel, min_vel

   integer :: lallo, i, j, k
   im = 128
   jm = 128
   km = 128

   allocate (vel_1(0:im, 0:jm, 0:km, 1:3), stat=lallo)
   allocate (rho_1(0:im, 0:jm, 0:km), stat=lallo)
   allocate (prs_1(0:im, 0:jm, 0:km), stat=lallo)
   allocate (tmp_1(0:im, 0:jm, 0:km), stat=lallo)

   allocate (vel_2(0:im, 0:jm, 0:km, 1:3), stat=lallo)
   allocate (rho_2(0:im, 0:jm, 0:km), stat=lallo)
   allocate (prs_2(0:im, 0:jm, 0:km), stat=lallo)
   allocate (tmp_2(0:im, 0:jm, 0:km), stat=lallo)

   call read_vtk_binary("flow_000010_cpu.vtk", rho_1, vel_1, prs_1, tmp_1)
   call read_vtk_binary("flow_000010_gpu.vtk", rho_2, vel_2, prs_2, tmp_2)

   allocate (vel_err(0:im, 0:jm, 0:km, 1:3), stat=lallo)
   allocate (rho_err(0:im, 0:jm, 0:km), stat=lallo)
   allocate (prs_err(0:im, 0:jm, 0:km), stat=lallo)
   allocate (tmp_err(0:im, 0:jm, 0:km), stat=lallo)


   vel_err = vel_1 - vel_2
   rho_err = rho_1 - rho_2
   prs_err = prs_1 - prs_2
   tmp_err = tmp_1 - tmp_2

   max_rho = maxval(rho_err)
   min_rho = minval(rho_err)

   max_prs = maxval(prs_err)
   min_prs = minval(prs_err)

   max_tmp = maxval(tmp_err)
   min_tmp = minval(tmp_err)

   max_vel = maxval(vel_err)
   min_vel = minval(vel_err)

   print *, 'rho: min = ', min_rho, ' max = ', max_rho
   print *, 'prs: min = ', min_prs, ' max = ', max_prs
   print *, 'tmp: min = ', min_tmp, ' max = ', max_tmp
   print *, 'vel: min = ', min_vel, ' max = ', max_vel
   
   call write_vtk_ascii("vel_err.vtk", rho_err, vel_err, prs_err, tmp_err)
end program compare_outputs
