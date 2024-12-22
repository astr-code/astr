module comvardef
  
  implicit none
  
  integer :: nthread,mthread,ncore
  
  integer :: im,jm,km,hm,numq,ndims
  
  integer :: nstep,maxstep,file_number
  real(8) :: time,deltat
  
  real(8) :: gamma,mach,reynolds,prandtl,ref_t,pinf
  real(8) :: const1,const2,const3,const4,const5,const6,const7
  real(8) :: dx,dy,dz
  
  real :: ctime(12)

  character(len=64) :: flowtype
  
  real(8),allocatable,dimension(:,:,:) :: rho,prs,tmp
  real(8),allocatable,dimension(:,:,:,:) :: x,q,qrhs,vel,dtmp
  real(8),allocatable,dimension(:,:,:,:,:) :: dvel
  
  contains
  
  subroutine alloarray
    
    if(im>0 .and. jm>0 .and. km>0) then
      ndims=3
    elseif(im>0 .and. jm>0 .and. km==0) then
      ndims=2
    elseif(im>0 .and. jm==0 .and. km==0) then
      ndims=1
    elseif(im==0 .and. jm==0 .and. km==0) then
      ndims=0
    endif

    const1=1.d0/(gamma*(gamma-1.d0)*mach**2)
    const2=gamma*mach**2
    const3=(gamma-1.d0)/3.d0*prandtl*(mach**2)
    const4=(gamma-1.d0)*mach**2*reynolds*prandtl
    const5=(gamma-1.d0)*mach**2
    const6=1.d0/(gamma-1.d0)
    const7=(gamma-1.d0)*mach**2*Reynolds*prandtl
    
    pinf=1.d0/const2

    allocate( x(0:im,0:jm,0:km,1:3))
    allocate(   q(-hm:im+hm,-hm:jm+hm,-hm:km+hm,1:numq))
    allocate(qrhs(0:im,0:jm,0:km,1:numq))

    allocate( rho(-hm:im+hm,-hm:jm+hm,-hm:km+hm))
    allocate( vel(-hm:im+hm,-hm:jm+hm,-hm:km+hm,1:3))
    allocate( prs(-hm:im+hm,-hm:jm+hm,-hm:km+hm))
    allocate( tmp(-hm:im+hm,-hm:jm+hm,-hm:km+hm))
    allocate(dvel(0:im,0:jm,0:km,1:3,1:3))
    allocate(dtmp(0:im,0:jm,0:km,1:3))

    print*,' ** common array allocated'

  end subroutine alloarray
  
end module comvardef