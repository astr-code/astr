module comvardef
  
  use constdef
  
  implicit none
  
  integer :: nthread,mthread,ncore
  
  integer :: im,jm,km,hm,numq,ndims
  
  integer :: nstep,maxstep,file_number
  real(rtype) :: time,deltat
  
  real(rtype) :: gamma,mach,reynolds,prandtl,ref_t,pinf
  real(rtype) :: const1,const2,const3,const4,const5,const6,const7
  real(rtype) :: dx,dy,dz
  
  real :: ctime(12)

  character(len=64) :: flowtype
  
  real(rtype),allocatable,dimension(:,:,:) :: rho,prs,tmp
  real(rtype),allocatable,dimension(:,:,:,:) :: x,q,qrhs,vel,dtmp
  real(rtype),allocatable,dimension(:,:,:,:,:) :: dvel
  
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

    const1=1._rtype/(gamma*(gamma-1._rtype)*mach**2)
    const2=gamma*mach**2
    const3=(gamma-1._rtype)/3._rtype*prandtl*(mach**2)
    const4=(gamma-1._rtype)*mach**2*reynolds*prandtl
    const5=(gamma-1._rtype)*mach**2
    const6=1._rtype/(gamma-1._rtype)
    const7=(gamma-1._rtype)*mach**2*Reynolds*prandtl
    
    pinf=1._rtype/const2

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