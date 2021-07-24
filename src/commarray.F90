!+---------------------------------------------------------------------+
!| This module is to define common array.                              |
!+---------------------------------------------------------------------+
!| CHANGE RECORD                                                       |
!| -------------                                                       |
!| 06-02-2021  | Created by J. Fang                                    |
!+---------------------------------------------------------------------+
module commarray
  !
  use commtype, only : nodcel
  !
  implicit none
  !
  real(8),allocatable,dimension(:,:,:,:) :: x,q,qrhs,vel,spc,dtmp
  real(8),allocatable,dimension(:,:,:) :: jacob,rho,prs,tmp
  real(8),allocatable,dimension(:,:,:,:,:) :: dxi,dvel,dspc
  real(8),allocatable,dimension(:,:) :: lspg_imin,lspg_imax,           &
                                        lspg_jmin,lspg_jmax,           &
                                        lspg_kmin,lspg_kmax
  integer,allocatable,dimension(:,:,:) :: nodestat
  logical,allocatable,dimension(:,:,:) :: lsolid
  type(nodcel),allocatable,dimension(:,:,:) :: cell
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
  !|              lspg_* | length of sponge layer                      |
  !|            nodestat | node status: fluid or solid                 |
  !+---------------------+---------------------------------------------+
  real(8),allocatable :: acctest_ref(:)
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
    use commvar, only : im,jm,km,hm,numq,num_species,ndims
    !
    ! local data
    integer :: lallo
    !
    allocate( x(-hm:im+hm,-hm:jm+hm,-hm:km+hm,1:3),stat=lallo)
    if(lallo.ne.0) stop ' !! error at allocating x'
    !
    allocate( nodestat(-hm:im+hm,-hm:jm+hm,-hm:km+hm),stat=lallo)
    if(lallo.ne.0) stop ' !! error at allocating nodestat'
    !
    allocate( jacob(-hm:im+hm,-hm:jm+hm,-hm:km+hm),stat=lallo)
    if(lallo.ne.0) stop ' !! error at allocating jacob'
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
    allocate(lsolid(-hm:im+hm,-hm:jm+hm,-hm:km+hm),stat=lallo)
    if(lallo.ne.0) stop ' !! error at allocating lsolid'
    !
    if(ndims==1) then
      allocate( cell(1:im,0:0,0:0),stat=lallo)
    elseif(ndims==2) then
      allocate( cell(1:im,1:jm,0:0),stat=lallo)
    elseif(ndims==3) then
      allocate( cell(1:im,1:jm,1:km),stat=lallo)
    endif
    if(lallo.ne.0) stop ' !! error at allocating cell'
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