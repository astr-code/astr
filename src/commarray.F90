!+---------------------------------------------------------------------+
!| This module is to define common array.                          |
!+---------------------------------------------------------------------+
!| CHANGE RECORD                                                       |
!| -------------                                                       |
!| 06-02-2021  | Created by J. Fang                                    |
!+---------------------------------------------------------------------+
module commarray
  !
  implicit none
  !
  real(8),allocatable,dimension(:,:,:,:) :: x
  !+---------------------+---------------------------------------------+
  !|                   x | coordinates.                                |
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
    use commvar, only : im,jm,km,nh
    !
    ! local data
    integer :: lallo
    !
    allocate( x(-nh:im+nh,-nh:jm+nh,-nh:km+nh,1:3),stat=lallo)
    !
    if(lallo.ne.0) stop ' !! error at allocating x'
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