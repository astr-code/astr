module comvardef
  !
  implicit none
  !
  integer :: nthread,mthread,ncore
  
  integer :: im=64,jm=64,km=64,hm=5,numq=5
  !
  integer :: nstep
  real(8) :: time,deltat
  !
  real(8) :: gamma,mach,reynolds,prandtl,ref_t
  real(8) :: const1,const2,const3,const4,const5,const6,const7
  real(8) :: dx,dy,dz
  !
  real :: ctime(12)
  !
  real(8),allocatable,dimension(:,:,:) :: rho,prs,tmp
  real(8),allocatable,dimension(:,:,:,:) :: x,q,qrhs,vel,dtmp
  real(8),allocatable,dimension(:,:,:,:,:) :: dvel
  !
  contains
  !
  subroutine alloarray
    !
    allocate( x(0:im,0:jm,0:km,1:3))
    allocate(   q(-hm:im+hm,-hm:jm+hm,-hm:km+hm,1:numq))
    allocate(qrhs(0:im,0:jm,0:km,1:numq))
    !
    allocate( rho(-hm:im+hm,-hm:jm+hm,-hm:km+hm))
    allocate( vel(-hm:im+hm,-hm:jm+hm,-hm:km+hm,1:3))
    allocate( prs(-hm:im+hm,-hm:jm+hm,-hm:km+hm))
    allocate( tmp(-hm:im+hm,-hm:jm+hm,-hm:km+hm))
    allocate(dvel(0:im,0:jm,0:km,1:3,1:3))
    allocate(dtmp(0:im,0:jm,0:km,1:3))
    !
    print*,' ** common array allocated'
    !
  end subroutine alloarray
  !
end module comvardef