!+---------------------------------------------------------------------+
!| This module is to define common array.                              |
!+---------------------------------------------------------------------+
!| CHANGE RECORD                                                       |
!| -------------                                                       |
!| 06-02-2021  | Created by J. Fang                                    |
!+---------------------------------------------------------------------+
module commarray
  !
  implicit none
  !
  real(8),allocatable,dimension(:,:,:,:) :: x,q,qrhs,vel,spc,dtmp
  real(8),allocatable,dimension(:,:,:) :: jacob,celvol,rho,prs,tmp
  real(8),allocatable,dimension(:,:,:,:,:) :: dxi,dvel,dspc
  real(8),allocatable,dimension(:,:,:) :: rom,u1m,u2m,u3m,pm,tm,u11,   &
                                          u22,u33,u12,u13,u23,tt,pp
  !+---------------------+---------------------------------------------+
  !|                   x | coordinates.                                |
  !|               jacob | geometrical Jacobian.                       |
  !|                 dxi | geometrical transform matrix                |
  !|                   q | indepedent conservation variables.          |
  !|                qrhs | R.H.S of equations.                         |
  !|                 rho | density.                                    |
  !|                 prs | pressure.                                   |
  !|                 tmp | temperature.                                |
  !|                 vel | velocity.                                   |
  !|                 spc | species.                                    |
  !+---------------------+---------------------------------------------+
  !
  contains
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to allocate common array.                 |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 07-02-2021  | Created by J. Fang STFC Daresbury Laboratory        |
  !+-------------------------------------------------------------------+
  subroutine allocommarray
    !
    use commvar, only : im,jm,km,hm,numq,num_species
    !
    ! local data
    integer :: lallo
    !
    allocate( x(-hm:im+hm,-hm:jm+hm,-hm:km+hm,1:3),stat=lallo)
    if(lallo.ne.0) stop ' !! error at allocating x'
    !
    allocate( jacob(-hm:im+hm,-hm:jm+hm,-hm:km+hm),stat=lallo)
    if(lallo.ne.0) stop ' !! error at allocating jacob'
    !
    allocate( celvol(1:im,1:jm,1:km),stat=lallo)
    if(lallo.ne.0) stop ' !! error at allocating celvol'
    !
    allocate( dxi(-hm:im+hm,-hm:jm+hm,-hm:km+hm,1:3,1:3),stat=lallo)
    if(lallo.ne.0) stop ' !! error at allocating dxi'
    !
    allocate( q(-hm:im+hm,-hm:jm+hm,-hm:km+hm,1:numq),stat=lallo)
    if(lallo.ne.0) stop ' !! error at allocating q'
    !
    allocate( rho(-hm:im+hm,-hm:jm+hm,-hm:km+hm),stat=lallo)
    if(lallo.ne.0) stop ' !! error at allocating rho'
    !
    allocate( prs(-hm:im+hm,-hm:jm+hm,-hm:km+hm),stat=lallo)
    if(lallo.ne.0) stop ' !! error at allocating prs'
    !
    allocate( tmp(-hm:im+hm,-hm:jm+hm,-hm:km+hm),stat=lallo)
    if(lallo.ne.0) stop ' !! error at allocating tmp'
    !
    allocate( vel(-hm:im+hm,-hm:jm+hm,-hm:km+hm,1:3),stat=lallo)
    if(lallo.ne.0) stop ' !! error at allocating vel'
    !
    allocate( spc(-hm:im+hm,-hm:jm+hm,-hm:km+hm,1:num_species),stat=lallo)
    if(lallo.ne.0) stop ' !! error at allocating spc'
    !
    allocate(qrhs(0:im,0:jm,0:km,1:numq),stat=lallo)
    if(lallo.ne.0) stop ' !! error at allocating qrhs'
    !
    allocate(dvel(0:im,0:jm,0:km,1:3,1:3),stat=lallo)
    if(lallo.ne.0) stop ' !! error at allocating dvel'
    !
    allocate(dspc(0:im,0:jm,0:km,1:num_species,1:3),stat=lallo)
    if(lallo.ne.0) stop ' !! error at allocating dvel'
    !
    allocate(dtmp(0:im,0:jm,0:km,1:3),stat=lallo)
    if(lallo.ne.0) stop ' !! error at allocating dvel'
    !
  end subroutine allocommarray
  !+-------------------------------------------------------------------+
  !| The end of the subroutine allocommarray.                          |
  !+-------------------------------------------------------------------+
  !
end module commarray
!+---------------------------------------------------------------------+
!| The end of the module commarray.                                    |
!+---------------------------------------------------------------------+