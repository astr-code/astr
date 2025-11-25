module pastr_flowvis
    !!
    use iso_fortran_env, only: wp => real64
    !!  Common flow-visualisation diagnostics:
    !!     - Vorticity (ω)
    !!     - Q-criterion
    !!     - Schlieren magnitude (|∇ρ|)
    !!     - Density gradient magnitude
    !!
    !!  Assumptions:
    !!     - All arrays are real(8), shape (nx, ny, nz)
    !!     - Uniform grid spacing (dx, dy, dz)
    !!
    implicit none
    private

    public :: schlieren_xy
    ! public :: compute_vorticity, compute_Q_criterion
    ! public :: compute_schlieren, compute_grad_rho_mag
    ! public :: central_diff_x, central_diff_y, central_diff_z

contains

  !==============================================================
  ! CENTRAL DERIVATIVES (PERIODIC or ZERO-GRADIENT ghost layer)
  !==============================================================
  ! function central_diff_x(f, dx) result(df_dx)
  !     real(8), intent(in) :: f(:,:,:)
  !     real(8), intent(in) :: dx
  !     real(8) :: df_dx(size(f,1), size(f,2), size(f,3))
  !     integer :: i, j, k, nx, ny, nz
  
  !     nx = size(f,1); ny = size(f,2); nz = size(f,3)
  
  !     do k=1,nz
  !         do j=1,ny
  !             do i=2,nx-1
  !                 df_dx(i,j,k) = (f(i+1,j,k) - f(i-1,j,k)) / (2.0d0*dx)
  !             end do
  !             df_dx(1,j,k)    = df_dx(2,j,k)
  !             df_dx(nx,j,k)   = df_dx(nx-1,j,k)
  !         end do
  !     end do
  ! end function central_diff_x
  
  ! function central_diff_y(f, dy) result(df_dy)
  !     real(8), intent(in) :: f(:,:,:)
  !     real(8), intent(in) :: dy
  !     real(8) :: df_dy(size(f,1), size(f,2), size(f,3))
  !     integer :: i, j, k, nx, ny, nz
  
  !     nx = size(f,1); ny = size(f,2); nz = size(f,3)
  
  !     do k=1,nz
  !         do i=1,nx
  !             do j=2,ny-1
  !                 df_dy(i,j,k) = (f(i,j+1,k) - f(i,j-1,k)) / (2.0d0*dy)
  !             end do
  !             df_dy(i,1,k)  = df_dy(i,2,k)
  !             df_dy(i,ny,k) = df_dy(i,ny-1,k)
  !         end do
  !     end do
  ! end function central_diff_y
  
  ! function central_diff_z(f, dz) result(df_dz)
  !     real(8), intent(in) :: f(:,:,:)
  !     real(8), intent(in) :: dz
  !     real(8) :: df_dz(size(f,1), size(f,2), size(f,3))
  !     integer :: i, j, k, nx, ny, nz
  
  !     nx = size(f,1); ny = size(f,2); nz = size(f,3)
  
  !     do j=1,ny
  !         do i=1,nx
  !             do k=2,nz-1
  !                 df_dz(i,j,k) = (f(i,j,k+1) - f(i,j,k-1)) / (2.0d0*dz)
  !             end do
  !             df_dz(i,j,1)  = df_dz(i,j,2)
  !             df_dz(i,j,nz) = df_dz(i,j,nz-1)
  !         end do
  !     end do
  ! end function central_diff_z
  
  ! !==============================================================
  ! ! COMPUTE VORTICITY VECTOR
  ! !==============================================================
  ! subroutine compute_vorticity(u, v, w, dx, dy, dz, wx, wy, wz)
  !     !!
  !     !! Vorticity = ∇ × U
  !     !!   wx = dw/dy - dv/dz
  !     !!   wy = du/dz - dw/dx
  !     !!   wz = dv/dx - du/dy
  !     !!
  !     real(8), intent(in)  :: u(:,:,:), v(:,:,:), w(:,:,:)
  !     real(8), intent(in)  :: dx, dy, dz
  !     real(8), intent(out) :: wx(:,:,:), wy(:,:,:), wz(:,:,:)
  
  !     real(8), allocatable :: du_dx(:,:,:), du_dy(:,:,:), du_dz(:,:,:)
  !     real(8), allocatable :: dv_dx(:,:,:), dv_dy(:,:,:), dv_dz(:,:,:)
  !     real(8), allocatable :: dw_dx(:,:,:), dw_dy(:,:,:), dw_dz(:,:,:)
  
  !     allocate(du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, dw_dx, dw_dy, dw_dz)
  
  !     du_dx = central_diff_x(u, dx)
  !     du_dy = central_diff_y(u, dy)
  !     du_dz = central_diff_z(u, dz)
  
  !     dv_dx = central_diff_x(v, dx)
  !     dv_dy = central_diff_y(v, dy)
  !     dv_dz = central_diff_z(v, dz)
  
  !     dw_dx = central_diff_x(w, dx)
  !     dw_dy = central_diff_y(w, dy)
  !     dw_dz = central_diff_z(w, dz)
  
  !     wx = dw_dy - dv_dz
  !     wy = du_dz - dw_dx
  !     wz = dv_dx - du_dy
  
  !     deallocate(du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, dw_dx, dw_dy, dw_dz)
  ! end subroutine compute_vorticity
  
  ! !==============================================================
  ! ! Q-CRITERION
  ! ! Q = 0.5 (||Ω||^2 - ||S||^2)
  ! !==============================================================
  ! subroutine compute_Q_criterion(u, v, w, dx, dy, dz, Q)
  !     real(8), intent(in) :: u(:,:,:), v(:,:,:), w(:,:,:)
  !     real(8), intent(in) :: dx, dy, dz
  !     real(8), intent(out) :: Q(:,:,:)
  
  !     real(8), allocatable :: du_dx, du_dy, du_dz
  !     real(8), allocatable :: dv_dx, dv_dy, dv_dz
  !     real(8), allocatable :: dw_dx, dw_dy, dw_dz
  
  !     real(8) :: Sxx, Syy, Szz, Sxy, Sxz, Syz
  !     real(8) :: Oxy, Oxz, Oyz
  !     integer :: i,j,k, nx,ny,nz
  
  !     allocate(du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, dw_dx, dw_dy, dw_dz)
  
  !     du_dx = central_diff_x(u, dx)
  !     du_dy = central_diff_y(u, dy)
  !     du_dz = central_diff_z(u, dz)
  
  !     dv_dx = central_diff_x(v, dx)
  !     dv_dy = central_diff_y(v, dy)
  !     dv_dz = central_diff_z(v, dz)
  
  !     dw_dx = central_diff_x(w, dx)
  !     dw_dy = central_diff_y(w, dy)
  !     dw_dz = central_diff_z(w, dz)
  
  !     nx=size(u,1); ny=size(u,2); nz=size(u,3)
  
  !     do k=1,nz
  !         do j=1,ny
  !             do i=1,nx
  
  !                 ! Strain rate tensor (S)
  !                 Sxx = du_dx(i,j,k)
  !                 Syy = dv_dy(i,j,k)
  !                 Szz = dw_dz(i,j,k)
  
  !                 Sxy = 0.5d0*(du_dy(i,j,k) + dv_dx(i,j,k))
  !                 Sxz = 0.5d0*(du_dz(i,j,k) + dw_dx(i,j,k))
  !                 Syz = 0.5d0*(dv_dz(i,j,k) + dw_dy(i,j,k))
  
  !                 ! Rotation tensor (Ω = antisymmetric)
  !                 Oxy = 0.5d0*(du_dy(i,j,k) - dv_dx(i,j,k))
  !                 Oxz = 0.5d0*(du_dz(i,j,k) - dw_dx(i,j,k))
  !                 Oyz = 0.5d0*(dv_dz(i,j,k) - dw_dy(i,j,k))
  
  !                 Q(i,j,k) = 0.5d0 * ( &
  !                     (Oxy**2 + Oxz**2 + Oyz**2)*2.0d0   &
  !                   - (Sxx**2 + Syy**2 + Szz**2        &
  !                     + 2.0d0*(Sxy**2 + Sxz**2 + Syz**2)) )
  
  !             end do
  !         end do
  !     end do
  
  !     deallocate(du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, dw_dx, dw_dy, dw_dz)
  ! end subroutine compute_Q_criterion
  
  !==============================================================
  ! SCHLIEREN: magnitude of density gradient
  !==============================================================
  function schlieren_xy(rho, x, y) result(schl)

      use pastr_commvar, only : im,jm
      use pastr_gradients, only : grad_xy

      real(wp), intent(in) :: rho(0:im,0:jm)
      real(wp), intent(in),optional :: x(0:im,0:jm), y(0:im,0:jm)
      real(wp) :: schl(0:im,0:jm)
  
      real(wp) :: dro(1:2,0:im,0:jm)
      real(wp) ::var1
      real(wp) ::rnsmin=0._wp,rnsmax=0._wp,c1=1._wp,c2=10._wp

      integer :: i,j
      logical :: lfex

      save rnsmin,rnsmax
  
      if(present(x) .and. present(y)) then
        dro=grad_xy(rho,x,y)
      else
        dro=grad_xy(rho)
      endif

      schl(:,:)=sqrt(dro(1,:,:)**2+dro(2,:,:)**2)

    
      if(rnsmax<=1.d-10 .and. rnsmin<1.d-10) then
        
        inquire(file='rnsdef.txt',exist=lfex)
        if(lfex) then
          open(12,file='rnsdef.txt')
          read(12,*)rnsmin,rnsmax
          close(12)
          print*,'rnsmin',rnsmin,'rnsmax',rnsmax
        else
          rnsmax=maxval(schl)
          rnsmin=minval(schl)
        endif
      endif
      print*,'rnsmin',rnsmin,'rnsmax',rnsmax

      schl=c1*exp(-c2*(schl-rnsmin)/(rnsmax-rnsmin))
  
  end function schlieren_xy
  
  ! !==============================================================
  ! ! GRADIENT MAGNITUDE FOR GENERAL SCALAR
  ! !==============================================================
  ! subroutine compute_grad_rho_mag(f, dx, dy, dz, gradmag)
  !     real(8), intent(in) :: f(:,:,:)
  !     real(8), intent(in) :: dx, dy, dz
  !     real(8), intent(out):: gradmag(:,:,:)
  
  !     real(8), allocatable :: fx, fy, fz
  
  !     allocate(fx, fy, fz)
  
  !     fx = central_diff_x(f, dx)
  !     fy = central_diff_y(f, dy)
  !     fz = central_diff_z(f, dz)
  
  !     gradmag = sqrt(fx**2 + fy**2 + fz**2)
  
  !     deallocate(fx, fy, fz)
  ! end subroutine compute_grad_rho_mag

end module pastr_flowvis
