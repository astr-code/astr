!+---------------------------------------------------------------------+
!| This module contains subroutines and functions from the x3d2 solver |
!| The source:  https://github.com/xcompact3d/x3d2                     |
!| ==============                                                      |
!| CHANGE RECORD                                                       |
!| -------------                                                       |
!| 14-04-2025: Created by J. Fang @ Liverpool, UK                      |
!+---------------------------------------------------------------------+
module x3d2

  use iso_fortran_env, only: stderr => error_unit
  use commvar, only: dp

  implicit none

  integer :: dims_block(3)
  integer, parameter :: SZ=16,n_halo=4
  integer, parameter :: BC_PERIODIC = 0, BC_NEUMANN = 1, BC_DIRICHLET = 2,BC_HALO=4

  integer, parameter :: RDR_X2Y = 12, RDR_X2Z = 13, RDR_Y2X = 21, &
                        RDR_Y2Z = 23, RDR_Z2X = 31, RDR_Z2Y = 32, &
                        RDR_C2X = 41, RDR_C2Y = 42, RDR_C2Z = 43, &
                        RDR_X2C = 14, RDR_Y2C = 24, RDR_Z2C = 34
  integer, parameter :: DIR_X = 1, DIR_Y = 2, DIR_Z = 3, DIR_C = 4

  integer, protected :: &
    rdr_map(4, 4) = reshape([0, RDR_Y2X, RDR_Z2X, RDR_C2X, &
                             RDR_X2Y, 0, RDR_Z2Y, RDR_C2Y, &
                             RDR_X2Z, RDR_Y2Z, 0, RDR_C2Z, &
                             RDR_X2C, RDR_Y2C, RDR_Z2C, 0], shape=[4, 4])

  

  interface get_index_reordering
    procedure get_index_reordering_rdr, get_index_reordering_dirs
  end interface

  type :: tdsops_t
      !! Tridiagonal Solver Operators class.
      !!
      !! Operator arrays are preprocessed in this class based on the arguments
      !! provided. dist_fw and dist_bw are used in the first phase of the
      !! distributed tridiagonal solver algorithm. dist_sa and dist_sc are used
      !! in the final substitution phase. See the kernels_dist.f90 files in the
      !! relevant backend folders.
      !! coeff arrays define the specific rules of building the RHS
      !! corresponding to the tridiagonal system to be solved, and used only in
      !! the first phase of the distributed algorithm when building the RHS.
      !! If a boundary condition is defined then coeffs_s and coeffs_e differ
      !! from coeffs array and define the RHS rule for the first and last 4
      !! entries in the tridiagonal system (n_halo = 4).
      !!
      !! This class does not know about the current rank or its relative
      !! location among other ranks. All the operator arrays here are used when
      !! executing a distributed tridiagonal solver phase one or two.
    real(dp), allocatable, dimension(:) :: dist_fw, dist_bw, & !! fw/bw phase
                                           dist_sa, dist_sc, & !! back subs.
                                           dist_af !! the auxiliary factors
    real(dp), allocatable, dimension(:) :: thom_f, thom_s, thom_w, thom_p
    real(dp), allocatable :: stretch(:), stretch_correct(:)
    real(dp), allocatable :: coeffs(:), coeffs_s(:, :), coeffs_e(:, :)
    real(dp) :: alpha, a, b, c = 0._dp, d = 0._dp !! Compact scheme coeffs
    logical :: periodic
    integer :: n_tds !! Tridiagonal system size
    integer :: n_rhs !! Right-hand-side builder size
    integer :: move = 0 !! move between vertices and cell centres
    integer :: n_halo !! number of halo points
  contains
    procedure :: deriv_1st
    procedure :: preprocess_dist, preprocess_thom
  end type tdsops_t

  interface tdsops_t
    module procedure tdsops_init
  end interface tdsops_t

  type(tdsops_t),target :: tdsops_i,tdsops_j,tdsops_k

  contains

  subroutine x3d2init(idim,jdim,kdim)

    use parallel, only : mpirank,mpistop

    integer,intent(in) :: idim,jdim,kdim

    dims_block(1)=idim
    dims_block(2)=jdim
    dims_block(3)=kdim

    tdsops_i = tdsops_init(idim, operation='first-deriv', scheme='compact6', bc_start=BC_PERIODIC, bc_end=BC_PERIODIC)

    tdsops_j = tdsops_init(jdim, operation='first-deriv', scheme='compact6', bc_start=BC_PERIODIC, bc_end=BC_PERIODIC)

    tdsops_k = tdsops_init(kdim, operation='first-deriv', scheme='compact6', bc_start=BC_PERIODIC, bc_end=BC_PERIODIC)

    if(mpirank==0) print*,' ** x3d2 solver initialised' 

  end subroutine x3d2init

  function ddf_x3d2(u,nvar,dir) result(du)
    
    use commvar, only: im,jm,km,hm
    use parallel,only: sendrecv_fields,mpileft,mpiright,mpistop,datasync,mpidown,mpiup,mpiback,mpifront

    integer,intent(in) :: nvar,dir
    
    real(dp), dimension(:,:,:,:), intent(in) :: u
    real(dp), dimension(0:im,0:jm,0:km, 1:nvar) :: du

    ! local data
    integer :: i,j,k,n,n_groups,lb1,lb2,lb3,ub1,ub2,ub3,narray,pprev,pnext
    type(tdsops_t) :: tdsops
    real(dp),dimension(:,:,:,:),allocatable :: ubuff
    real(dp),dimension(:,:,:),allocatable :: u_x,du_x,u_temp
    real(dp),dimension(:,:,:),allocatable :: u_send_s,u_send_e,u_recv_s, &
                                             u_recv_e,send_s,send_e,     &
                                             recv_s,recv_e

    if(dir==1) then
      lb1=-hm; ub1=im+hm
      lb2=0;   ub2=jm
      lb3=0;   ub3=km

      n_groups=jm*km/SZ
      narray=im
      pprev=mpileft
      pnext=mpiright
      tdsops=tdsops_i
    elseif(dir==2) then
      lb1=0;   ub1=im
      lb2=-hm; ub2=jm+hm
      lb3=0;   ub3=km

      n_groups=im*km/SZ
      narray=jm
      pprev=mpidown
      pnext=mpiup
      tdsops=tdsops_j
    elseif(dir==3) then
      lb1=0;   ub1=im
      lb2=0;   ub2=jm
      lb3=-hm; ub3=km+hm

      n_groups=im*jm/SZ
      narray=km
      pprev=mpiback
      pnext=mpifront
      tdsops=tdsops_k
    else
      stop ' !! dir error @ function ddf_x3d2'
    endif

    allocate(ubuff(lb1:ub1,lb2:ub2,lb3:ub3,1:nvar))

    allocate(u_x(SZ,narray,n_groups),du_x(SZ,narray,n_groups),u_temp(1:im,1:jm,1:km))

    allocate (u_send_s(SZ, n_halo, n_groups))
    allocate (u_send_e(SZ, n_halo, n_groups))
    allocate (u_recv_s(SZ, n_halo, n_groups))
    allocate (u_recv_e(SZ, n_halo, n_groups))

    allocate (send_s(SZ, 1, n_groups), send_e(SZ, 1, n_groups))
    allocate (recv_s(SZ, 1, n_groups), recv_e(SZ, 1, n_groups))

    ubuff=u

    do n=1,nvar

      u_temp=ubuff(0:im-1,0:jm-1,0:km-1,n)

      call reorder_omp(u_x, u_temp, dir_from=DIR_C, dir_to=dir)

      do k = 1, n_groups
      do j = 1, n_halo
      do i = 1, SZ
        u_send_s(i, j, k) = u_x(i, j, k)
        u_send_e(i, j, k) = u_x(i, narray - n_halo + j, k)
      end do
      end do
      end do

      call sendrecv_fields(u_recv_s, u_recv_e, u_send_s, u_send_e, &
                           SZ*n_halo*n_groups, pprev, pnext)


      call exec_dist_tds_compact(du_x, u_x, u_recv_s, u_recv_e, &
                                 send_s, send_e, recv_s, recv_e, &
                                 tdsops,  pprev, pnext, n_groups)

      call reorder_omp(u_temp, du_x, dir_from=dir, dir_to=DIR_C)

      du(0:im-1,0:jm-1,0:km-1,n)=u_temp

    enddo

    call datasync(du,nvar)

    return

  end function ddf_x3d2

  subroutine reorder_omp(u_,u,dir_from,dir_to)

    implicit none

    real(dp), intent(inout) :: u_(:,:,:)
    real(dp), intent(in) :: u(:,:,:)
    integer, intent(in) :: dir_from, dir_to
    integer, dimension(3) :: dims
    integer :: i, j, k
    integer :: out_i, out_j, out_k

    dims = shape(u)

    ! print*,dir_from, dir_to, direction

    !$omp parallel do private(out_i, out_j, out_k) collapse(2)
    do k = 1, dims(3)
      do j = 1, dims(2)
        do i = 1, dims(1)
          call get_index_reordering(out_i, out_j, out_k, i, j, k, dir_from, dir_to, dims_block)
          ! print*,i,j,k,out_i, out_j, out_k
          u_(out_i, out_j, out_k) = u(i, j, k)
        end do
      end do
    end do
    !$omp end parallel do

  end subroutine reorder_omp

  pure subroutine get_index_ijk(i, j, k, dir_i, dir_j, dir_k, dir, &
                                nx_padded, ny_padded, nz_padded)
      !! Get cartesian index from application storage directional one
    integer, intent(out) :: i, j, k                   ! cartesian indices
    integer, intent(in) :: dir_i, dir_j, dir_k        ! application storage indices
    integer, intent(in) :: dir                        ! direction of the applicatino storage indices
    integer, intent(in) :: nx_padded, ny_padded, nz_padded ! dimensions of the block

    select case (dir)
    case (DIR_X)
      i = dir_j
      j = mod(dir_k - 1, ny_padded/SZ)*SZ + dir_i
      k = 1 + (dir_k - 1)/(ny_padded/SZ)
    case (DIR_Y)
      i = mod(dir_k - 1, nx_padded/SZ)*SZ + dir_i
      j = dir_j
      k = 1 + (dir_k - 1)/(nx_padded/SZ)
    case (DIR_Z)
      i = mod(dir_k - 1, nx_padded/SZ)*SZ + dir_i
      j = 1 + (dir_k - 1)/(nx_padded/SZ)
      k = dir_j
    case (DIR_C)
      i = dir_i
      j = dir_j
      k = dir_k
    end select

  end subroutine get_index_ijk

  pure subroutine get_index_dir(dir_i, dir_j, dir_k, i, j, k, dir, &
                                nx_padded, ny_padded, nz_padded)
      !! Get application storage directional index from cartesian index
    integer, intent(out) :: dir_i, dir_j, dir_k        ! application storage indices
    integer, intent(in) :: i, j, k                     ! cartesian indices
    integer, intent(in) :: dir                        ! direction of the application storage indices
    integer, intent(in) :: nx_padded, ny_padded, nz_padded ! dimensions of the block

    select case (dir)
    case (DIR_X)
      dir_i = mod(j - 1, SZ) + 1
      dir_j = i
      dir_k = (ny_padded/SZ)*(k - 1) + 1 + (j - 1)/SZ
    case (DIR_Y)
      dir_i = mod(i - 1, SZ) + 1
      dir_j = j
      dir_k = (nx_padded/SZ)*(k - 1) + 1 + (i - 1)/SZ
    case (DIR_Z)
      dir_i = mod(i - 1, SZ) + 1
      dir_j = k
      dir_k = (nx_padded/SZ)*(j - 1) + 1 + (i - 1)/SZ
    case (DIR_C)
      dir_i = i
      dir_j = j
      dir_k = k
    end select

  end subroutine get_index_dir

  pure subroutine get_dirs_from_rdr(dir_from, dir_to, rdr_dir)
    integer, intent(out) :: dir_from, dir_to
    integer, intent(in) :: rdr_dir
    integer, dimension(2) :: dirs

    dirs = findloc(rdr_map, rdr_dir)
    dir_from = dirs(1)
    dir_to = dirs(2)

  end subroutine

  pure subroutine get_index_reordering_dirs(out_i, out_j, out_k, in_i, in_j, in_k, dir_from, dir_to, dims)
      !! Converts a set of application storage directional index to an other direction.
      !! The two directions are defined by the reorder_dir variable, RDR_X2Y will go from storage in X to Y etc.
    integer, intent(out) :: out_i, out_j, out_k         ! new indices in the application storage
    integer, intent(in) :: in_i, in_j, in_k             ! original indices
    integer, intent(in) :: dir_from, dir_to
    integer, intent(in), dimension(3) :: dims
    integer :: i, j, k        ! Intermediary cartesian indices

    call get_index_ijk(i, j, k, in_i, in_j, in_k, dir_from, &
                       dims(1), dims(2), dims(3))
    call get_index_dir(out_i, out_j, out_k, i, j, k, dir_to, &
                       dims(1), dims(2), dims(3))

  end subroutine get_index_reordering_dirs

  pure subroutine get_index_reordering_rdr(out_i, out_j, out_k, &
                                           in_i, in_j, in_k, reorder_dir,dims)
    integer, intent(out) :: out_i, out_j, out_k         ! new indices in the application storage
    integer, intent(in) :: in_i, in_j, in_k             ! original indices
    integer, intent(in) :: reorder_dir
    integer, intent(in), dimension(3) :: dims
    integer :: dir_from, dir_to

    call get_dirs_from_rdr(dir_from, dir_to, reorder_dir)
    call get_index_reordering(out_i, out_j, out_k, in_i, in_j, in_k, &
                              dir_from, dir_to, dims)

  end subroutine get_index_reordering_rdr

  subroutine exec_dist_tds_compact( &
    du, u, u_recv_s, u_recv_e, du_send_s, du_send_e, du_recv_s, du_recv_e, &
    tdsops, pprev, pnext, n_groups)
    

    use parallel,only: sendrecv_fields,mpileft,mpiright,mpistop



    ! du = d(u)
    real(dp), dimension(:, :, :), intent(out) :: du
    real(dp), dimension(:, :, :), intent(in) :: u, u_recv_s, u_recv_e

    ! The ones below are intent(out) just so that we can write data in them,
    ! not because we actually need the data they store later where this
    ! subroutine is called. We absolutely don't care about the data they pass back
    real(dp), dimension(:, :, :), intent(out) :: &
      du_send_s, du_send_e, du_recv_s, du_recv_e

    type(tdsops_t), intent(in) :: tdsops
    integer, intent(in) :: pprev, pnext
    integer, intent(in) :: n_groups

    integer :: n_data
    integer :: k

    n_data = SZ*n_groups

    !$omp parallel do
    do k = 1, n_groups
      call der_univ_dist( &
        du(:, :, k), du_send_s(:, :, k), du_send_e(:, :, k), &
        u(:, :, k), u_recv_s(:, :, k), u_recv_e(:, :, k), &
        tdsops%n_tds, tdsops%n_rhs, &
        tdsops%coeffs_s, tdsops%coeffs_e, tdsops%coeffs, &
        tdsops%dist_fw, tdsops%dist_bw, tdsops%dist_af &
        )
    end do
    !$omp end parallel do

    ! halo exchange for 2x2 systems
    call sendrecv_fields(du_recv_s, du_recv_e, du_send_s, du_send_e, &
                         n_data, pprev, pnext)

    !$omp parallel do
    do k = 1, n_groups
      call der_univ_subs( &
        du(:, :, k), du_recv_s(:, :, k), du_recv_e(:, :, k), &
        tdsops%n_tds, tdsops%dist_sa, tdsops%dist_sc, tdsops%stretch &
        )
    end do
    !$omp end parallel do

  end subroutine exec_dist_tds_compact


  subroutine der_univ_dist( &
    du, send_u_s, send_u_e, u, u_s, u_e, &
    n_tds, n_rhs, coeffs_s, coeffs_e, coeffs, ffr, fbc, faf &
    )
    implicit none

    ! Arguments
    real(dp), intent(out), dimension(:, :) :: du, send_u_s, send_u_e
    real(dp), intent(in), dimension(:, :) :: u, u_s, u_e
    integer, intent(in) :: n_tds, n_rhs
    real(dp), intent(in), dimension(:, :) :: coeffs_s, coeffs_e ! start/end
    real(dp), intent(in), dimension(:) :: coeffs
    real(dp), intent(in), dimension(:) :: ffr, fbc, faf

    ! Local variables
    integer :: i, j

    real(dp) :: c_m4, c_m3, c_m2, c_m1, c_j, c_p1, c_p2, c_p3, c_p4, &
                alpha, last_r

    ! store bulk coeffs in the registers
    c_m4 = coeffs(1); c_m3 = coeffs(2); c_m2 = coeffs(3); c_m1 = coeffs(4)
    c_j = coeffs(5)
    c_p1 = coeffs(6); c_p2 = coeffs(7); c_p3 = coeffs(8); c_p4 = coeffs(9)
    last_r = ffr(1)

    !$omp simd
    do i = 1, SZ
      du(i, 1) = coeffs_s(1, 1)*u_s(i, 1) &
                 + coeffs_s(2, 1)*u_s(i, 2) &
                 + coeffs_s(3, 1)*u_s(i, 3) &
                 + coeffs_s(4, 1)*u_s(i, 4) &
                 + coeffs_s(5, 1)*u(i, 1) &
                 + coeffs_s(6, 1)*u(i, 2) &
                 + coeffs_s(7, 1)*u(i, 3) &
                 + coeffs_s(8, 1)*u(i, 4) &
                 + coeffs_s(9, 1)*u(i, 5)
      du(i, 1) = du(i, 1)*faf(1)
      du(i, 2) = coeffs_s(1, 2)*u_s(i, 2) &
                 + coeffs_s(2, 2)*u_s(i, 3) &
                 + coeffs_s(3, 2)*u_s(i, 4) &
                 + coeffs_s(4, 2)*u(i, 1) &
                 + coeffs_s(5, 2)*u(i, 2) &
                 + coeffs_s(6, 2)*u(i, 3) &
                 + coeffs_s(7, 2)*u(i, 4) &
                 + coeffs_s(8, 2)*u(i, 5) &
                 + coeffs_s(9, 2)*u(i, 6)
      du(i, 2) = du(i, 2)*faf(2)
      du(i, 3) = coeffs_s(1, 3)*u_s(i, 3) &
                 + coeffs_s(2, 3)*u_s(i, 4) &
                 + coeffs_s(3, 3)*u(i, 1) &
                 + coeffs_s(4, 3)*u(i, 2) &
                 + coeffs_s(5, 3)*u(i, 3) &
                 + coeffs_s(6, 3)*u(i, 4) &
                 + coeffs_s(7, 3)*u(i, 5) &
                 + coeffs_s(8, 3)*u(i, 6) &
                 + coeffs_s(9, 3)*u(i, 7)
      du(i, 3) = ffr(3)*(du(i, 3) - faf(3)*du(i, 2))
      du(i, 4) = coeffs_s(1, 4)*u_s(i, 4) &
                 + coeffs_s(2, 4)*u(i, 1) &
                 + coeffs_s(3, 4)*u(i, 2) &
                 + coeffs_s(4, 4)*u(i, 3) &
                 + coeffs_s(5, 4)*u(i, 4) &
                 + coeffs_s(6, 4)*u(i, 5) &
                 + coeffs_s(7, 4)*u(i, 6) &
                 + coeffs_s(8, 4)*u(i, 7) &
                 + coeffs_s(9, 4)*u(i, 8)
      du(i, 4) = ffr(4)*(du(i, 4) - faf(4)*du(i, 3))
    end do
    !$omp end simd

    ! alpha is always the same in the bulk region for us
    alpha = faf(5)

    do j = 5, n_rhs - 4
      !$omp simd
      do i = 1, SZ
        du(i, j) = c_m4*u(i, j - 4) + c_m3*u(i, j - 3) &
                   + c_m2*u(i, j - 2) + c_m1*u(i, j - 1) &
                   + c_j*u(i, j) &
                   + c_p1*u(i, j + 1) + c_p2*u(i, j + 2) &
                   + c_p3*u(i, j + 3) + c_p4*u(i, j + 4)
        du(i, j) = ffr(j)*(du(i, j) - alpha*du(i, j - 1))
      end do
      !$omp end simd
    end do

    !$omp simd
    do i = 1, SZ
      j = n_rhs - 3
      du(i, j) = coeffs_e(1, 1)*u(i, j - 4) &
                 + coeffs_e(2, 1)*u(i, j - 3) &
                 + coeffs_e(3, 1)*u(i, j - 2) &
                 + coeffs_e(4, 1)*u(i, j - 1) &
                 + coeffs_e(5, 1)*u(i, j) &
                 + coeffs_e(6, 1)*u(i, j + 1) &
                 + coeffs_e(7, 1)*u(i, j + 2) &
                 + coeffs_e(8, 1)*u(i, j + 3) &
                 + coeffs_e(9, 1)*u_e(i, 1)
      du(i, j) = ffr(j)*(du(i, j) - faf(j)*du(i, j - 1))
      j = n_rhs - 2
      du(i, j) = coeffs_e(1, 2)*u(i, j - 4) &
                 + coeffs_e(2, 2)*u(i, j - 3) &
                 + coeffs_e(3, 2)*u(i, j - 2) &
                 + coeffs_e(4, 2)*u(i, j - 1) &
                 + coeffs_e(5, 2)*u(i, j) &
                 + coeffs_e(6, 2)*u(i, j + 1) &
                 + coeffs_e(7, 2)*u(i, j + 2) &
                 + coeffs_e(8, 2)*u_e(i, 1) &
                 + coeffs_e(9, 2)*u_e(i, 2)
      du(i, j) = ffr(j)*(du(i, j) - faf(j)*du(i, j - 1))
      j = n_rhs - 1
      du(i, j) = coeffs_e(1, 3)*u(i, j - 4) &
                 + coeffs_e(2, 3)*u(i, j - 3) &
                 + coeffs_e(3, 3)*u(i, j - 2) &
                 + coeffs_e(4, 3)*u(i, j - 1) &
                 + coeffs_e(5, 3)*u(i, j) &
                 + coeffs_e(6, 3)*u(i, j + 1) &
                 + coeffs_e(7, 3)*u_e(i, 1) &
                 + coeffs_e(8, 3)*u_e(i, 2) &
                 + coeffs_e(9, 3)*u_e(i, 3)
      du(i, j) = ffr(j)*(du(i, j) - faf(j)*du(i, j - 1))
      j = n_rhs
      du(i, j) = coeffs_e(1, 4)*u(i, j - 4) &
                 + coeffs_e(2, 4)*u(i, j - 3) &
                 + coeffs_e(3, 4)*u(i, j - 2) &
                 + coeffs_e(4, 4)*u(i, j - 1) &
                 + coeffs_e(5, 4)*u(i, j) &
                 + coeffs_e(6, 4)*u_e(i, 1) &
                 + coeffs_e(7, 4)*u_e(i, 2) &
                 + coeffs_e(8, 4)*u_e(i, 3) &
                 + coeffs_e(9, 4)*u_e(i, 4)
      du(i, j) = ffr(j)*(du(i, j) - faf(j)*du(i, j - 1))
    end do
    !$omp end simd

    !$omp simd
    do i = 1, SZ
      send_u_e(i, 1) = du(i, n_tds)
    end do
    !$omp end simd

    ! Backward pass of the hybrid algorithm
    do j = n_tds - 2, 2, -1
      !$omp simd
      do i = 1, SZ
        du(i, j) = du(i, j) - fbc(j)*du(i, j + 1)
      end do
      !$omp end simd
    end do
    !$omp simd
    do i = 1, SZ
      du(i, 1) = last_r*(du(i, 1) - fbc(1)*du(i, 2))
      send_u_s(i, 1) = du(i, 1)
    end do
    !$omp end simd

  end subroutine der_univ_dist

  subroutine der_univ_subs(du, recv_u_s, recv_u_e, n, dist_sa, dist_sc, strch)
    implicit none

    ! Arguments
    real(dp), intent(out), dimension(:, :) :: du
    real(dp), intent(in), dimension(:, :) :: recv_u_s, recv_u_e
    real(dp), intent(in), dimension(:) :: dist_sa, dist_sc, strch
    integer, intent(in) :: n

    ! Local variables
    integer :: i, j!, b
    real(dp) :: ur, bl, recp
    real(dp), dimension(SZ) :: du_s, du_e

    !$omp simd
    do i = 1, SZ
      ! A small trick we do here is valid for symmetric Toeplitz matrices.
      ! In our case our matrices satisfy this criteria in the (5:n-4) region
      ! and as long as a rank has around at least 20 entries the assumptions
      ! we make here are perfectly valid.

      ! bl is the bottom left entry in the 2x2 matrix
      ! ur is the upper right entry in the 2x2 matrix

      ! Start
      ! At the start we have the 'bl', and assume 'ur'
      bl = dist_sa(1)
      ur = dist_sa(1)
      recp = 1._dp/(1._dp - ur*bl)
      du_s(i) = recp*(du(i, 1) - bl*recv_u_s(i, 1))

      ! End
      ! At the end we have the 'ur', and assume 'bl'
      bl = dist_sc(n)
      ur = dist_sc(n)
      recp = 1._dp/(1._dp - ur*bl)
      du_e(i) = recp*(du(i, n) - ur*recv_u_e(i, 1))
    end do
    !$omp end simd

    !$omp simd
    do i = 1, SZ
      du(i, 1) = du_s(i)*strch(1)
    end do
    !$omp end simd
    do j = 2, n - 1
      !$omp simd
      do i = 1, SZ
        du(i, j) = (du(i, j) - dist_sa(j)*du_s(i) - dist_sc(j)*du_e(i)) &
                   *strch(j)
      end do
      !$omp end simd
    end do
    !$omp simd
    do i = 1, SZ
      du(i, n) = du_e(i)*strch(n)
    end do
    !$omp end simd

  end subroutine der_univ_subs

  function tdsops_init( &
    n_tds,  operation, scheme, bc_start, bc_end, &
    stretch, stretch_correct, n_halo, from_to, sym, c_nu, nu0_nu &
    ) result(tdsops)
    !! Constructor function for the tdsops_t class.
    !!
    !! 'n_tds', 'delta', 'operation', 'scheme', 'bc_start', and 'bc_end' are
    !! necessary arguments. The remaining arguments are optional.
    !!
    !! 'stretch' is for obtaining the correct derivations in a stretched mesh
    !! 'stretch_correct' is for correcting the second derivative with the first
    !!
    !! 'from_to' is necessary for interpolation and staggared derivative, and
    !! it can be 'v2p' or 'p2v'.
    !! If the specific region the instance is operating is not a boundary
    !! region, then 'bc_start' and 'bc_end' are BC_HALO.
    !!
    !! 'sym' is relevant when the BC is free-slip. If sym is .true. then it
    !! means the field we operate on is assumed to be an even function
    !! (symmetric, cos type) accross the boundary. If it is .false. it means
    !! the field is assumed to be an odd function (anti-symmetric, sin type).
    !!
    !! 'c_nu', 'nu0_nu' are relevant when operation is second order
    !! derivative and scheme is compact6-hyperviscous.

    type(tdsops_t) :: tdsops !! return value of the function

    integer, intent(in) :: n_tds !! Tridiagonal system size
    character(*), intent(in) :: operation, scheme
    integer, intent(in) :: bc_start, bc_end !! Boundary Cond.
    real(dp), optional, intent(in) :: stretch(:) !! Stretching coefficients
    real(dp), optional, intent(in) :: stretch_correct(:) !! Stretch correction
    integer, optional, intent(in) :: n_halo !! Number of halo cells
    character(*), optional, intent(in) :: from_to !! 'v2p' or 'p2v'
    logical, optional, intent(in) :: sym !! (==npaire), only for Neumann BCs
    real(dp), optional, intent(in) :: c_nu, nu0_nu !! params for hypervisc.

    integer :: n, n_stencil

    tdsops%n_tds = n_tds

    ! we need special treatment in the right-hand-side build stage for
    ! the very last point in the domain if output length is smaller than
    ! the input length
    if (present(from_to)) then
      if ((bc_end == BC_NEUMANN .or. bc_end == BC_DIRICHLET) &
          .and. from_to == 'v2p') then
        tdsops%n_rhs = n_tds + 1
      else
        tdsops%n_rhs = n_tds
      end if
    else
      tdsops%n_rhs = n_tds
    end if

    if (present(n_halo)) then
      tdsops%n_halo = n_halo
      if (n_halo /= 4) then
        write (stderr, '("Warning: n_halo is set to ", i2, "be careful! &
                          &The default is 4 and there are quite a few &
                          &places where things are hardcoded assuming &
                          &n_halo is 4.")') n_halo
      end if
    else
      tdsops%n_halo = 4
    end if

    ! n_rhs >= n_tds, n is used when its better to allocate a larger size
    n = tdsops%n_rhs

    ! preprocessed coefficient arrays for the distributed algorithm
    allocate (tdsops%dist_fw(n), tdsops%dist_bw(n))
    allocate (tdsops%dist_sa(n), tdsops%dist_sc(n))
    allocate (tdsops%dist_af(n))

    ! preprocessed coefficient arrays for the Thomas algorithm
    allocate (tdsops%thom_f(n), tdsops%thom_s(n))
    allocate (tdsops%thom_w(n), tdsops%thom_p(n))

    ! RHS coefficient arrays
    n_stencil = 2*tdsops%n_halo + 1
    allocate (tdsops%coeffs(n_stencil))
    allocate (tdsops%coeffs_s(n_stencil, tdsops%n_halo))
    allocate (tdsops%coeffs_e(n_stencil, tdsops%n_halo))

    allocate (tdsops%stretch(n_tds))
    if (present(stretch)) then
      tdsops%stretch(:) = stretch(:)
    else
      tdsops%stretch(:) = 1._dp
    end if

    allocate (tdsops%stretch_correct(n_tds))
    if (present(stretch_correct)) then
      tdsops%stretch_correct(:) = stretch_correct(:)
    else
      tdsops%stretch_correct(:) = 0._dp
    end if

    tdsops%periodic = bc_start == BC_PERIODIC .and. bc_end == BC_PERIODIC

    if (operation == 'first-deriv') then
      call tdsops%deriv_1st( scheme, bc_start, bc_end, sym)
    else
      error stop 'operation is not defined'
    end if

    select case (from_to)
    case ('v2p')
      tdsops%move = 1
    case ('p2v')
      tdsops%move = -1
    case default
      tdsops%move = 0
    end select

    if (tdsops%dist_sa(n_tds) > 1d-16) then
      print *, 'There are ', n_tds, 'points in a subdomain, it may be too few!'
      print *, 'The entry distributed solver disregards in "' &
        //operation//'" operation is:', tdsops%dist_sa(n_tds)
      print *, 'It may result in numerical errors with the distributed solver!'
    end if

  end function tdsops_init

  subroutine deriv_1st(self,  scheme, bc_start, bc_end, sym)

    class(tdsops_t), intent(inout) :: self
    character(*), intent(in) :: scheme
    integer, intent(in) :: bc_start, bc_end
    logical, optional, intent(in) :: sym

    real(dp), allocatable :: dist_b(:)
    real(dp) :: alpha, afi, bfi
    integer :: i, n, n_halo
    logical :: symmetry

    if (self%n_halo < 2) error stop 'First derivative require n_halo >= 2'

    if (present(sym)) then
      symmetry = sym
    else
      symmetry = .false.
    end if

    ! alpha is alfa

    select case (scheme)
    case ('compact6')
      alpha = 1._dp/3._dp
      afi = 7._dp/9._dp
      bfi = 1._dp/36._dp
    case default
      error stop 'scheme is not defined'
    end select

    self%alpha = alpha
    self%a = afi; self%b = bfi

    self%coeffs(:) = [0._dp, 0._dp, -bfi, -afi, &
                      0._dp, &
                      afi, bfi, 0._dp, 0._dp]

    do i = 1, self%n_halo
      self%coeffs_s(:, i) = self%coeffs(:)
      self%coeffs_e(:, i) = self%coeffs(:)
    end do

    self%dist_sa(:) = alpha; self%dist_sc(:) = alpha

    n = self%n_tds
    n_halo = self%n_halo

    allocate (dist_b(self%n_rhs))
    dist_b(:) = 1._dp

    select case (bc_start)
    case (BC_NEUMANN)
      if (symmetry) then
        ! sym == .true.; d(uu)/dx, dv/dx, dw/dx
        !                d(vv)/dy, du/dy, dw/dy
        !                d(ww)/dz, du/dz, dv/dz
        self%dist_sa(1) = 0._dp
        self%dist_sc(1) = 0._dp
        self%coeffs_s(:, 1) = [0._dp, 0._dp, 0._dp, 0._dp, &
                               0._dp, &
                               0._dp, 0._dp, 0._dp, 0._dp]
        self%coeffs_s(:, 2) = [0._dp, 0._dp, 0._dp, -afi, &
                               -bfi, &
                               afi, bfi, 0._dp, 0._dp]
      else
        ! sym == .false.; d(uv)/dx, d(uw)/dx, du/dx
        !                 d(vu)/dy, d(vw)/dy, dv/dy
        !                 d(wu)/dz, d(wv)/dz, dw/dz
        self%dist_sa(1) = 0._dp
        self%dist_sc(1) = 2*alpha
        self%coeffs_s(:, 1) = [0._dp, 0._dp, 0._dp, 0._dp, &
                               0._dp, &
                               2*afi, 2*bfi, 0._dp, 0._dp]
        self%coeffs_s(:, 2) = [0._dp, 0._dp, 0._dp, -afi, &
                               bfi, &
                               afi, bfi, 0._dp, 0._dp]
      end if
    case (BC_DIRICHLET)
      ! first line
      self%dist_sa(1) = 0._dp
      self%dist_sc(1) = 2._dp
      self%coeffs_s(:, 1) = [0._dp, 0._dp, 0._dp, 0._dp, &
                             -2.5_dp, &
                             2._dp, 0.5_dp, 0._dp, 0._dp]
      self%coeffs_s(:, 1) = self%coeffs_s(:, 1)
      ! second line
      self%dist_sa(2) = 0.25_dp
      self%dist_sc(2) = 0.25_dp
      self%coeffs_s(:, 2) = [0._dp, 0._dp, 0._dp, -0.75_dp, &
                             0._dp, &
                             0.75_dp, 0._dp, 0._dp, 0._dp]
      self%coeffs_s(:, 2) = self%coeffs_s(:, 2)
    end select

    select case (bc_end)
    case (BC_NEUMANN)
      if (symmetry) then
        ! sym == .true.; d(uu)/dx, dv/dx, dw/dx
        !                d(vv)/dy, du/dy, dw/dy
        !                d(ww)/dz, du/dz, dv/dz
        self%dist_sa(n) = 0._dp
        self%dist_sc(n) = 0._dp
        self%coeffs_e(:, n_halo) = [0._dp, 0._dp, 0._dp, 0._dp, &
                                    0._dp, &
                                    0._dp, 0._dp, 0._dp, 0._dp]
        self%coeffs_e(:, n_halo - 1) = [0._dp, 0._dp, -bfi, -afi, &
                                        bfi, &
                                        afi, 0._dp, 0._dp, 0._dp]
      else
        ! sym == .false.; d(uv)/dx, d(uw)/dx, du/dx
        !                 d(vu)/dy, d(vw)/dy, dv/dy
        !                 d(wu)/dz, d(wv)/dz, dw/dz
        self%dist_sa(n) = 2*alpha
        self%dist_sc(n) = 0._dp
        self%coeffs_e(:, n_halo) = [0._dp, 0._dp, -2*bfi, -2*afi, &
                                    0._dp, &
                                    0._dp, 0._dp, 0._dp, 0._dp]
        self%coeffs_e(:, n_halo - 1) = [0._dp, 0._dp, -bfi, -afi, &
                                        -bfi, &
                                        afi, 0._dp, 0._dp, 0._dp]
      end if
    case (BC_DIRICHLET)
      ! last line
      self%dist_sa(n) = 2._dp
      self%dist_sc(n) = 0._dp
      self%coeffs_e(:, n_halo) = [0._dp, 0._dp, -0.5_dp, -2._dp, &
                                  2.5_dp, &
                                  0._dp, 0._dp, 0._dp, 0._dp]
      self%coeffs_e(:, n_halo) = self%coeffs_e(:, n_halo)
      ! second last line
      self%dist_sa(n - 1) = 0.25_dp
      self%dist_sc(n - 1) = 0.25_dp
      self%coeffs_e(:, n_halo - 1) = [0._dp, 0._dp, 0._dp, -0.75_dp, &
                                      0._dp, &
                                      0.75_dp, 0._dp, 0._dp, 0._dp]
      self%coeffs_e(:, n_halo - 1) = self%coeffs_e(:, n_halo - 1)
    end select

    call self%preprocess_thom(dist_b)
    call self%preprocess_dist(dist_b)

  end subroutine deriv_1st


  subroutine preprocess_dist(self, dist_b)

    class(tdsops_t), intent(inout) :: self

    real(dp), dimension(:), intent(in) :: dist_b

    integer :: i

    ! Ref DOI: 10.1109/MCSE.2021.3130544
    ! Algorithm 3 in page 4
    ! First two lines first
    do i = 1, 2
      self%dist_sa(i) = self%dist_sa(i)/dist_b(i)
      self%dist_sc(i) = self%dist_sc(i)/dist_b(i)
      self%dist_bw(i) = self%dist_sc(i)
      self%dist_af(i) = 1._dp/dist_b(i)
    end do

    ! Then the remaining in the forward pass
    do i = 3, self%n_tds
      ! Algorithm 3 in ref obtains 'r' coeffs on the fly in line 7.
      ! As we have to solve many RHSs with the same tridiagonal system,
      ! it is better to do a preprocessing first.
      ! So lets store 'r' coeff in dist_fw array.
      self%dist_fw(i) = 1._dp/(dist_b(i) &
                               - self%dist_sa(i)*self%dist_sc(i - 1))
      ! dist_af is 'a_i' in line 7 of Algorithm 3 in ref.
      self%dist_af(i) = self%dist_sa(i)
      ! We store a_i^* and c_i^* in dist_sa and dist_sc because
      ! we need them later in the substitution phase.
      self%dist_sa(i) = -self%dist_fw(i)*self%dist_sa(i) &
                        *self%dist_sa(i - 1)
      self%dist_sc(i) = self%dist_fw(i)*self%dist_sc(i)
    end do

    ! backward pass starting in line 12 of Algorithm 3.
    do i = self%n_tds - 2, 2, -1
      self%dist_sa(i) = self%dist_sa(i) &
                        - self%dist_sc(i)*self%dist_sa(i + 1)
      self%dist_bw(i) = self%dist_sc(i)
      self%dist_sc(i) = -self%dist_sc(i)*self%dist_sc(i + 1)
    end do

    ! Line 17 and 18 are tricky
    ! First we have a new 'r', we need it.
    ! And for 'r' we need c_0^*...
    ! Now examine closely, c_0^* is set in line 4 and never overwritten!
    ! So we can use dist_sc(1) as is in place of c_0^*.
    ! We need to store this new 'r' somewhere ...
    ! dist_fw(1) is never used, so store this extra 'r' factor here instead
    self%dist_fw(1) = 1._dp/(1._dp - self%dist_sc(1)*self%dist_sa(2))

    ! Finally Line 19 and 20 in Algorithm 3 in ref.
    self%dist_sa(1) = self%dist_fw(1)*self%dist_sa(1)
    self%dist_sc(1) = -self%dist_fw(1)*self%dist_sc(1)*self%dist_sc(2)

  end subroutine preprocess_dist

  subroutine preprocess_thom(self, b)

    class(tdsops_t), intent(inout) :: self
    real(dp), dimension(:), intent(in) :: b

    integer :: i, n

    n = self%n_tds

    self%thom_w = b
    self%thom_f = self%dist_sc
    if (self%periodic) then
      self%thom_w(1) = 2._dp
      self%thom_w(n) = 1._dp + self%alpha*self%alpha
    end if

    self%thom_s(1) = 0._dp
    do i = 2, n
      self%thom_s(i) = self%dist_sa(i)/self%thom_w(i - 1)
      self%thom_w(i) = self%thom_w(i) - self%thom_f(i - 1)*self%thom_s(i)
    end do
    do i = 1, n
      self%thom_w(i) = 1._dp/self%thom_w(i)
    end do

    self%thom_p = [-1._dp, (0._dp, i=2, n - 1), self%alpha]
    do i = 2, n
      self%thom_p(i) = self%thom_p(i) - self%thom_p(i - 1)*self%thom_s(i)
    end do
    self%thom_p(n) = self%thom_p(n)*self%thom_w(n)
    do i = n - 1, 1, -1
      self%thom_p(i) = self%thom_w(i)*(self%thom_p(i) &
                                       - self%thom_f(i)*self%thom_p(i + 1))
    end do

  end subroutine preprocess_thom

end module x3d2
!+---------------------------------------------------------------------+
!| The end of the module x3d2.                                         |
!+---------------------------------------------------------------------+