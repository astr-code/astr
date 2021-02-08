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
  real(8),allocatable,dimension(:,:,:,:) :: x,q
  real(8),allocatable,dimension(:,:,:) :: jacob
  real(8),allocatable,dimension(:,:,:,:,:) :: dxi
  !+---------------------+---------------------------------------------+
  !|                   x | coordinates.                                |
  !|               jacob | geometrical jacobian.                       |
  !|                 dxi | geometrical transform matrix                |
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
    use commvar, only : im,jm,km,hm,numq
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
    allocate( dxi(-hm:im+hm,-hm:jm+hm,-hm:km+hm,1:3,1:3),stat=lallo)
    if(lallo.ne.0) stop ' !! error at allocating dxi'
    !
    allocate( q(-hm:im+hm,-hm:jm+hm,-hm:km+hm,1:numq),stat=lallo)
    if(lallo.ne.0) stop ' !! error at allocating q'
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