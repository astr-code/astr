!+---------------------------------------------------------------------+
! this module contains subroutines and functions related to tridiagonal| 
! solver within a parallel envionment                                  |
!+---------------------------------------------------------------------+
!| CHANGE RECORD                                                       |
!| -------------                                                       |
!| 30-12-2024  | Created by J. Fang                                    |
!+---------------------------------------------------------------------+
module tridiagonal
  
  implicit none

  contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! This subroutine is used to initial for solving the tridiagonal 
  ! martix with two layer of boundary scheme:
  ! A*x=b
  !   |1,c11,.........................|
  !   |c22,1,c12,.....................|
  !   |..c23,1,c13,...................|
  !   |...............................|
  ! A=|...,c2i,1,c1i..................|
  !   |...............................|
  !   |.........,c(2,m-2),1,c(1,m-2)..|
  !   |...........,c(2,m-1),1,c(1,m-1)|
  !   |.........................,c2m,1|
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Writen by Fang Jian, 2008-11-04.
  ! Redesigned by Fang Jian for more general cases, 2024-12-30.
  ! ref: https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine ptds_ini_com(c)
    
    real(8),intent(inout) :: c(:,:)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! c: output array
    ! dim: input dat
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! local data
    integer :: i,dim

    dim=size(c,2)

    c(3,1)=c(1,1)
    do i=2,dim-1
      c(3,i)=c(1,i)/(1.d0-c(2,i)*c(3,i-1))
    enddo
    
    c(4,1)=1.d0
    do i=2,dim
      c(4,i)=1.d0/(1.d0-c(2,i)*c(3,i-1))
    enddo

    return

  end subroutine ptds_ini_com

  function ptds_cal_com(b,c) result(x)
    
    real(8),intent(in) :: c(:,:),b(:)
    real(8) :: x(size(b))
    
    ! local data
    integer :: i,dim
    real(8),allocatable :: c4(:)

    dim=size(b)

    allocate(c4(1:dim)) 
    
    c4(1)=b(1)
    do i=2,dim
      c4(i)=(b(i)-c(2,i)*c4(i-1))*c(4,i)
    enddo

    x(dim)=c4(dim)
    do i=dim-1,1,-1
      x(i)=c4(i)-c(3,i)*x(i+1)
    enddo

    deallocate(c4)

    return

  end function ptds_cal_com



  ! function solve_tridiag(a,c,d,n) result(x)
  !   ! https://en.wikibooks.org/wiki/Algorithm_Implementation/Linear_Algebra/Tridiagonal_matrix_algorithm
  !   !	 a - sub-diagonal (means it is the diagonal below the main diagonal)
  !   !	 b - the main diagonal
  !   !	 c - sup-diagonal (means it is the diagonal above the main diagonal)
  !   !	 d - right part
  !   !	 x - the answer
  !   !	 n - number of equations

  !   integer,intent(in) :: n
  !   real(8),dimension(n),intent(in) :: a,c,d
  !   real(8),dimension(n),intent(out) :: x
  !   real(8),dimension(n) :: cp,dp
  !   real(8) :: m
  !   integer i

  !   ! initialize c-prime and d-prime
  !   cp(1) = c(1)/1.d0
  !   dp(1) = d(1)/1.d0
  !   ! solve for vectors c-prime and d-prime
  !   do i = 2,n
  !     m = 1.d0-cp(i-1)*a(i)
  !     cp(i) = c(i)/m
  !     dp(i) = (d(i)-dp(i-1)*a(i))/m
  !   end do

  !   ! initialize x
  !   x(n) = dp(n)
  !   ! solve for x from the vectors c-prime and d-prime
  !   do i = n-1, 1, -1
  !     x(i) = dp(i)-cp(i)*x(i+1)
  !   end do

  ! end function solve_tridiag

end module tridiagonal
