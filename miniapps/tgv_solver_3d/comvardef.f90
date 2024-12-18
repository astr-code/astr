module comvardef
  !
  implicit none
  !
  integer :: im,jm,km,hm,numq
  integer :: ndim,ih,jh,kh

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
    if(im==0 .and. jm==0 .and. km==0) then
      ndim=0
      ih=0
      jh=0
      kh=0
    elseif(im>0 .and. jm==0 .and. km==0) then
      ndim=1
      ih=hm
      jh=0
      kh=0
    elseif(im>0 .and. jm>0 .and. km==0) then
      ndim=2
      ih=hm
      jh=hm
      kh=0
    elseif(im>0 .and. jm>0 .and. km>0) then
      ndim=3
      ih=hm
      jh=hm
      kh=hm
    else
      stop ' error ndim set error @ alloarray'
    endif
    write(*,'(A,I0)')'  ** dimension problem  :',ndim
    ! 
    allocate( x(0:im,0:jm,0:km,1:ndim))
    allocate( q(-ih:im+im,-jm:jm+jm,-kh:km+kh,1:numq))
    allocate(qrhs(0:im,0:jm,0:km,1:numq))
    !
    allocate( rho(-ih:im+im,-jm:jm+jm,-kh:km+kh))
    allocate( vel(-ih:im+im,-jm:jm+jm,-kh:km+kh,1:ndim))
    allocate( prs(-ih:im+im,-jm:jm+jm,-kh:km+kh))
    allocate( tmp(-ih:im+im,-jm:jm+jm,-kh:km+kh))
    allocate(dvel(0:im,0:jm,0:km,1:ndim,1:ndim))
    allocate(dtmp(0:im,0:jm,0:km,1:ndim))
    !
    print*,' ** common array allocated'
    !
  end subroutine alloarray
  !
end module comvardef