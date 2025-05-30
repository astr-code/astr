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
  real(8),allocatable,dimension(:,:,:,:) :: x,q,qrhs,vel,spc,dtmp,     &
                                            dgrid,vor
  real(8),allocatable,dimension(:,:,:) :: jacob,rho,prs,tmp
  real(8),allocatable,dimension(:,:,:,:,:) :: dxi,dvel,dspc
  real(8),allocatable,dimension(:,:,:) :: bnorm_i0,bnorm_im,bnorm_j0,  &
                                          bnorm_jm,bnorm_k0,bnorm_km
  real(8),allocatable,dimension(:,:,:) :: dis2wall
  integer,allocatable,dimension(:,:,:) :: nodestat
  logical,allocatable,dimension(:,:,:) :: lsolid,lshock,crinod
  type(nodcel),allocatable,dimension(:,:,:) :: cell
  !
  real(8),allocatable,dimension(:,:,:) :: sponge_damp_coef_i0,sponge_damp_coef_im, &
                                          sponge_damp_coef_j0,sponge_damp_coef_jm, &
                                          sponge_damp_coef_k0,sponge_damp_coef_km, &
                                          sponge_damp_coef
  !
  real(8),allocatable,dimension(:,:,:) :: tke,omg,miut,res12,ssf
  real(8),allocatable,dimension(:,:,:,:) :: dtke,domg
  !
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
  !|            lenspg_* | length of sponge layer                      |
  !|              xspg_* | starting x,y,z of sponge layer              |
  !|            nodestat | node status: fluid or solid                 |
  !|             bnorm_* | the vector of boundary normal direction,    |
  !|                     | towards the inside of the domain            |
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
    use commvar, only : im,jm,km,hm,numq,num_species,ndims,turbmode
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
    allocate(vor(0:im,0:jm,0:km,1:3),stat=lallo)
    if(lallo.ne.0) stop ' !! error at allocating vor'
    !
    allocate(lsolid(-hm:im+hm,-hm:jm+hm,-hm:km+hm),stat=lallo)
    if(lallo.ne.0) stop ' !! error at allocating lsolid'
    !
    allocate(crinod(0:im,0:jm,0:km),stat=lallo)
    if(lallo.ne.0) stop ' !! error at allocating crinod'
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