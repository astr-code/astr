!+---------------------------------------------------------------------+
!| This module contains subroutine of common functions.                |
!+---------------------------------------------------------------------+
!| CHANGE RECORD                                                       |
!| -------------                                                       |
!| 08-02-2021  | Created by J. Fang                                    |
!+---------------------------------------------------------------------+
module commfunc
  !
  use commvar, only : hm,kcutoff,ltimrpt
  use parallel, only: mpirank,mpistop,lio,lreport,ptime
  use constdef
  use utility,  only : timereporter

  implicit none
  !
  interface extrapolate
    module procedure extrapolate_1o
    module procedure extrapolate_2o
    module procedure extrapolate_3o
    module procedure extrapolate_4o
  end interface
  !
  interface deriv
    module procedure deriv_1o
    module procedure deriv_2o
    module procedure deriv_3o
    module procedure deriv_4o
  end interface
  !
  complex(8) :: ci=(0.d0,1.d0)
  !
  contains
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! This function is used to calculate the area of a quadrilateral.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  real(8) function arquad(x1,x2,x3,x4)
    !
    real(8),intent(in) :: x1(3),x2(3),x3(3),x4(3)
    !
    real(8) :: a,b,c,d,e
    !
    a=sqrt((x1(1)-x2(1))*(x1(1)-x2(1))+(x1(2)-x2(2))*(x1(2)-x2(2)))
    b=sqrt((x2(1)-x3(1))*(x2(1)-x3(1))+(x2(2)-x3(2))*(x2(2)-x3(2)))
    c=sqrt((x3(1)-x4(1))*(x3(1)-x4(1))+(x3(2)-x4(2))*(x3(2)-x4(2)))
    d=sqrt((x4(1)-x1(1))*(x4(1)-x1(1))+(x4(2)-x1(2))*(x4(2)-x1(2)))
    !
    e=sqrt((x1(1)-x3(1))*(x1(1)-x3(1))+(x1(2)-x3(2))*(x1(2)-x3(2)))
    !
    arquad=artria(a,b,e)+artria(c,d,e)
    !
  end function arquad
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! End of the function arquad.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! This function is used to calculate the area of a triangle.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  real(8) function artria(a,b,c)
    !
    real(8),intent(in) :: a,b,c
    !
    real(8) :: s
    !
    s=0.5d0*(a+b+c)
    artria=sqrt(s*(s-a)*(s-b)*(s-c))
    !
    return
    !
  end function artria
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! End of the function artria.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! This function is used to calculate the volume of a cell.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  real(8) function volhex(x1,x2,x3,x4,x5,x6,x7,x8 )
    !
    real(8),intent(in) :: x1(3),x2(3),x3(3),x4(3),x5(3),x6(3),x7(3),x8(3)
    !
    ! local data
    real(8) :: martix(3,3),centroid(3),pyramid(6,3,3)
    integer :: ii
    !
    centroid(1)=0.125d0*(x1(1)+x2(1)+x3(1)+x4(1)+x5(1)+x6(1)+x7(1)+x8(1))
    centroid(2)=0.125d0*(x1(2)+x2(2)+x3(2)+x4(2)+x5(2)+x6(2)+x7(2)+x8(2))
    centroid(3)=0.125d0*(x1(3)+x2(3)+x3(3)+x4(3)+x5(3)+x6(3)+x7(3)+x8(3))
    !
    pyramid(1,1,1)=centroid(1)-x1(1)
    pyramid(1,1,2)=centroid(2)-x1(2)
    pyramid(1,1,3)=centroid(3)-x1(3)
    pyramid(1,2,1)=x3(1)-x1(1)
    pyramid(1,2,2)=x3(2)-x1(2)
    pyramid(1,2,3)=x3(3)-x1(3)
    pyramid(1,3,1)=x2(1)-x4(1)
    pyramid(1,3,2)=x2(2)-x4(2)
    pyramid(1,3,3)=x2(3)-x4(3)
    !
    pyramid(2,1,1)=centroid(1)-x5(1)
    pyramid(2,1,2)=centroid(2)-x5(2)
    pyramid(2,1,3)=centroid(3)-x5(3)
    pyramid(2,2,1)=x7(1)-x5(1)
    pyramid(2,2,2)=x7(2)-x5(2)
    pyramid(2,2,3)=x7(3)-x5(3)
    pyramid(2,3,1)=x8(1)-x6(1)
    pyramid(2,3,2)=x8(2)-x6(2)
    pyramid(2,3,3)=x8(3)-x6(3)
 
    pyramid(3,1,1)=centroid(1)-x1(1)
    pyramid(3,1,2)=centroid(2)-x1(2)
    pyramid(3,1,3)=centroid(3)-x1(3)
    pyramid(3,2,1)=x5(1)-x4(1)
    pyramid(3,2,2)=x5(2)-x4(2)
    pyramid(3,2,3)=x5(3)-x4(3)
    pyramid(3,3,1)=x8(1)-x1(1)
    pyramid(3,3,2)=x8(2)-x1(2)
    pyramid(3,3,3)=x8(3)-x1(3) 
 
    pyramid(4,1,1)=centroid(1)-x2(1)
    pyramid(4,1,2)=centroid(2)-x2(2)
    pyramid(4,1,3)=centroid(3)-x2(3)
    pyramid(4,2,1)=x6(1)-x3(1)
    pyramid(4,2,2)=x6(2)-x3(2)
    pyramid(4,2,3)=x6(3)-x3(3)
    pyramid(4,3,1)=x2(1)-x7(1)
    pyramid(4,3,2)=x2(2)-x7(2)
    pyramid(4,3,3)=x2(3)-x7(3) 
 
    pyramid(5,1,1)=centroid(1)-x3(1)
    pyramid(5,1,2)=centroid(2)-x3(2)
    pyramid(5,1,3)=centroid(3)-x3(3)
    pyramid(5,2,1)=x8(1)-x3(1)
    pyramid(5,2,2)=x8(2)-x3(2)
    pyramid(5,2,3)=x8(3)-x3(3)
    pyramid(5,3,1)=x7(1)-x4(1)
    pyramid(5,3,2)=x7(2)-x4(2)
    pyramid(5,3,3)=x7(3)-x4(3) 
 
    pyramid(6,1,1)=centroid(1)-x1(1)
    pyramid(6,1,2)=centroid(2)-x1(2)
    pyramid(6,1,3)=centroid(3)-x1(3)
    pyramid(6,2,1)=x5(1)-x2(1)
    pyramid(6,2,2)=x5(2)-x2(2)
    pyramid(6,2,3)=x5(3)-x2(3)
    pyramid(6,3,1)=x1(1)-x6(1)
    pyramid(6,3,2)=x1(2)-x6(2)
    pyramid(6,3,3)=x1(3)-x6(3) 
    !
    volhex=0.d0
    do ii=1,6                            
      martix(1,1)=pyramid(ii,1,1)
      martix(1,2)=pyramid(ii,1,2)
      martix(1,3)=pyramid(ii,1,3)
      martix(2,1)=pyramid(ii,2,1)
      martix(2,2)=pyramid(ii,2,2)
      martix(2,3)=pyramid(ii,2,3)
      martix(3,1)=pyramid(ii,3,1)
      martix(3,2)=pyramid(ii,3,2)
      martix(3,3)=pyramid(ii,3,3)
      !
      volhex=volhex+determinant33(martix)*num1d6
    enddo  
    !
    return
    !
  end function volhex
  !
  pure real(8) function determinant33(a)
    real(8),intent(in) :: a(3,3)
    !
    determinant33=a(1,1)*a(2,2)*a(3,3)+a(1,2)*a(2,3)*a(3,1)+           &
                  a(1,3)*a(2,1)*a(3,2)-a(1,3)*a(2,2)*a(3,1)-           &
                  a(1,1)*a(2,3)*a(3,2)-a(1,2)*a(2,1)*a(3,3)
    return
  end function determinant33
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! End of the subroutine volhex.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! this subroutine is to extraoplation according to the gradient
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function gradextrp(qbou,q1st,dq) result(q)
    !
    real(8) :: q(4)
    real(8),intent(in) :: qbou,q1st
    real(8),intent(in),optional :: dq
    !
    real(8) :: dqr
    !
    if(present(dq)) then
      dqr=dq
    else
      dqr=qbou-q1st
    endif
    !
    q(1)=q1st+2.d0*dqr
    q(2)=-2.d0*q1st -3.d0*qbou +6.d0*q(1) -6.d0*dqr
    q(3)= 3.d0*q1st+10.d0*qbou-18.d0*q(1) +6.d0*q(2)+12.d0*dqr
    q(4)=-4.d0*q1st-num65d3*qbou+40.d0*q(1)-20.d0*q(2)+num20d3*q(3)-20.d0*dqr
    !
  end function gradextrp
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! End of the subroutine gradextrp.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !+-------------------------------------------------------------------+
  !| This function is to calculate local derivative at point 0.        |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 25-Feb-2021: Created by J. Fang @ STFC Daresbury Laboratory       |
  !+-------------------------------------------------------------------+
  function deriv_1o(v0,v1) result(dv)
    !
    real(8) :: dv
    real(8),intent(in) :: v0,v1
    !
    dv=-v0+v1
    !
    return
    !
  end function deriv_1o
  !
  function deriv_2o(v0,v1,v2) result(dv)
    !
    real(8) :: dv
    real(8),intent(in) :: v0,v1,v2
    !
    dv=-1.5d0*v0+2.d0*v1-0.5d0*v2
    !
    return
    !
  end function deriv_2o
  !
  function deriv_3o(v0,v1,v2,v3) result(dv)
    !
    real(8) :: dv
    real(8),intent(in) :: v0,v1,v2,v3
    !
    dv=-num11d6*v0+3.d0*v1-1.5d0*v2+num1d3*v3
    !
    return
    !
  end function deriv_3o
  !
  function deriv_4o(v0,v1,v2,v3,v4) result(dv)
    !
    real(8) :: dv
    real(8),intent(in) :: v0,v1,v2,v3,v4
    !
    dv=-num25d12*v0+4.d0*v1-3.d0*v2+num4d3*v3-0.25d0*v4
    !
    return
    !
  end function deriv_4o
  !+-------------------------------------------------------------------+
  !| The end of the function deriv.                                    |
  !+-------------------------------------------------------------------+
  !!
  !+-------------------------------------------------------------------+
  !| This function is to extrapolate variables according to gradient.  |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 24-Feb-2021: Created by J. Fang @ STFC Daresbury Laboratory       |
  !+-------------------------------------------------------------------+
  function extrapolate_1o(v1,dv) result(v0)
    !
    ! arguments
    real(8) :: v0
    real(8),intent(in) :: v1,dv
    !
    v0=v1-dv
    !
  end function extrapolate_1o
  !
  function extrapolate_2o(v1,v2,dv) result(v0)
    !
    ! arguments
    real(8) :: v0
    real(8),intent(in) :: v1,v2,dv
    !
    v0=num1d3*(4.d0*v1-v2-2.d0*dv)
    !
  end function extrapolate_2o
  !
  function extrapolate_3o(v1,v2,v3,dv) result(v0)
    !
    ! arguments
    real(8) :: v0
    real(8),intent(in) :: v1,v2,v3,dv
    !
    v0=num1d11*(18.d0*v1-9.d0*v2+2.d0*v3-6.d0*dv)
    !
  end function extrapolate_3o
  !
  function extrapolate_4o(v1,v2,v3,v4,dv) result(v0)
    !
    ! arguments
    real(8) :: v0
    real(8),intent(in) :: v1,v2,v3,v4,dv
    !
    v0=0.04d0*(48.d0*v1-36.d0*v2+16.d0*v3-3.d0*v4-12.d0*dv)
    !
  end function extrapolate_4o
  !+-------------------------------------------------------------------+
  !| The end of the function extrapolate.                              |
  !+-------------------------------------------------------------------+
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! This function is used for cubic interpolation.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  pure function cubic(y1,y2,y3,y4,x1,delta,xx)
    !
    real(8) :: cubic
    real(8),intent(in) :: y1,y2,y3,y4,x1,delta,xx
    real(8) :: a,b,c,d,varx
    !
    a=(y4-3.d0*y3+3.d0*y2-y1)/(6.d0*delta**3)
    b=((y3-2.d0*y2+y1)-6.d0*a*delta**3)/(2.d0*delta**2)
    c=(y2-y1-a*delta**3-b*delta**2)/delta
    d=y1
    !
    varx=xx-x1
    !
    cubic=a*varx**3+b*varx**2+c*varx+d
    !
    return
    !
  end function cubic
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! End of the function cubic.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! This function is the the inverse function of the function tanh(x).
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  pure real(8) function argtanh(var)
    !
    real(8),intent(in) :: var
    !
    argtanh=0.5d0*dlog((1.d0+var)/(1.d0-var))
    !
    return
    !
  end function argtanh
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! End of the function argtanh.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !+-------------------------------------------------------------------+
  !| This function is to do the cross product for a 3-D vector.        | 
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 04-07-2021  | Created by J. Fang @ Warrington                     |
  !+-------------------------------------------------------------------+
  pure function cross_product(a,b)
    !
    real(8) :: cross_product(3)
    real(8),dimension(3),intent(in) :: a(3), b(3)
  
    cross_product(1) = a(2) * b(3) - a(3) * b(2)
    cross_product(2) = a(3) * b(1) - a(1) * b(3)
    cross_product(3) = a(1) * b(2) - a(2) * b(1)
    !
    return
    !
  end function cross_product
  !+-------------------------------------------------------------------+
  !| The end of the function cross_product.                            |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This function is to calculate the distance between 2 points.      | 
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 04-07-2021  | Created by J. Fang @ Warrington                     |
  !+-------------------------------------------------------------------+
  pure real(8) function dis2point(p1,p2)
    !
    ! arguments
    real(8),intent(in) :: p1(:),p2(:)
    !
    ! local data
    integer :: n,nd
    !
    nd=size(p1)
    !
    dis2point=0.d0
    do n=1,nd
      dis2point=dis2point+(p1(n)-p2(n))**2
    enddo
    !
    dis2point=sqrt(dis2point)
    !
    return
    !
  end function dis2point
  !
  pure real(8) function dis2point2(p1,p2)
    !
    ! arguments
    real(8),intent(in) :: p1(:),p2(:)
    !
    ! local data
    integer :: n,nd
    !
    nd=size(p1)
    !
    dis2point2=0.d0
    do n=1,nd
      dis2point2=dis2point2+(p1(n)-p2(n))**2
    enddo
    !
    return
    !
  end function dis2point2
  !+-------------------------------------------------------------------+
  !| The end of the subroutine dis2point.                              |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This function is to calculate the area of a triangle using        |
  !| Heron's formula.                                                  | 
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 04-07-2021  | Created by J. Fang @ Warrington                     |
  !+-------------------------------------------------------------------+
  pure real(8) function areatriangle(p1,p2,p3)
    !
    ! arguments
    real(8),intent(in) :: p1(3),p2(3),p3(3)
    !
    ! local data
    real(8) :: a,b,c,var1
    !
    a=dis2point(p1,p2)
    b=dis2point(p2,p3)
    c=dis2point(p1,p3)
    !
    var1=0.5d0*(a+b+c)
    !
    areatriangle=sqrt(var1*(var1-a)*(var1-b)*(var1-c))
    !
    return
    !
  end function areatriangle
  !+-------------------------------------------------------------------+
  !| The end of the subroutine areatriangle.                           |
  !+-------------------------------------------------------------------+
  !!
  !+-------------------------------------------------------------------+
  !| This subroutine is to determin if two points the same.            | 
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 04-07-2021  | Created by J. Fang @ Warrington                     |
  !+-------------------------------------------------------------------+
  pure logical function samepoint(p1,p2)
    !
    ! arguments
    real(8),intent(in) :: p1(3),p2(3)
    !
    ! local data
    real(8) :: epsilon
    !
    epsilon=1.d-16
    !
    if(abs(p1(1)-p2(1))<=epsilon .and. abs(p1(2)-p2(2))<=epsilon .and. &
       abs(p1(3)-p2(3))<=epsilon) then
      samepoint=.true.
    else
      samepoint=.false.
    endif
    !
    return
    !
  end function samepoint
  !+-------------------------------------------------------------------+
  !| The end of the subroutine areatriangle.                           |
  !+-------------------------------------------------------------------+
  !!
  !+-------------------------------------------------------------------+
  !| ref: http://fortranwiki.org/fortran/show/Matrix+inversion
  !+-------------------------------------------------------------------+
  !| LAPACK is great when you wish to invert huge N×N matrices, 
  !| but it can be really slow for inverting smaller 2×2, 3×3, and 
  !| 4×4 matrices. For my use case, where I need to invert billions of 
  !| 2×2 and 4×4 matrices instead of a few large N×N matrices, I got a 
  !| 30% speedup of my program replacing the LAPACK calls by direct 
  !| calculations of the matrix inversions. I have attached the code 
  !| that I’ve used for the 2×2, 3×3, and 4×4 cases below. 
  !| The 2×2 version is quite easy to derive analytically. 
  !| The 3×3 and 4×4 versions are based on the subroutines M33INV and
  !|  M44INV by David G. Simpson; I just converted them from subroutines
  !| to pure functions
  !+-------------------------------------------------------------------+
  pure function matinv2(A) result(B)
    !! Performs a direct calculation of the inverse of a 2×2 matrix.
    real(8), intent(in) :: A(2,2)   !! Matrix
    real(8)             :: B(2,2)   !! Inverse matrix
    real(8)             :: detinv

    ! Calculate the inverse determinant of the matrix
    detinv = 1.d0/(A(1,1)*A(2,2) - A(1,2)*A(2,1))

    ! Calculate the inverse of the matrix
    B(1,1) = +detinv * A(2,2)
    B(2,1) = -detinv * A(2,1)
    B(1,2) = -detinv * A(1,2)
    B(2,2) = +detinv * A(1,1)
    !
    return
    !
  end function matinv2
  !!
  pure function matinv3(A) result(B)
    !! Performs a direct calculation of the inverse of a 3×3 matrix.
    real(8), intent(in) :: A(3,3)   !! Matrix
    real(8)             :: B(3,3)   !! Inverse matrix
    real(8)             :: detinv

    ! Calculate the inverse determinant of the matrix
    detinv = 1.d0/(A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2)&
                 - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1)&
                 + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1))

    ! Calculate the inverse of the matrix
    B(1,1) = +detinv * (A(2,2)*A(3,3) - A(2,3)*A(3,2))
    B(2,1) = -detinv * (A(2,1)*A(3,3) - A(2,3)*A(3,1))
    B(3,1) = +detinv * (A(2,1)*A(3,2) - A(2,2)*A(3,1))
    B(1,2) = -detinv * (A(1,2)*A(3,3) - A(1,3)*A(3,2))
    B(2,2) = +detinv * (A(1,1)*A(3,3) - A(1,3)*A(3,1))
    B(3,2) = -detinv * (A(1,1)*A(3,2) - A(1,2)*A(3,1))
    B(1,3) = +detinv * (A(1,2)*A(2,3) - A(1,3)*A(2,2))
    B(2,3) = -detinv * (A(1,1)*A(2,3) - A(1,3)*A(2,1))
    B(3,3) = +detinv * (A(1,1)*A(2,2) - A(1,2)*A(2,1))
    !
    return
    !
  end function matinv3
  !!
  pure function matinv4(A) result(B)
    !! Performs a direct calculation of the inverse of a 4×4 matrix.
    real(8), intent(in) :: A(4,4)   !! Matrix
    real(8)             :: B(4,4)   !! Inverse matrix
    real(8)             :: detinv

    ! Calculate the inverse determinant of the matrix
    detinv = &
     1.d0/(A(1,1)*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2)))&
         - A(1,2)*(A(2,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))&
         + A(1,3)*(A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))&
         - A(1,4)*(A(2,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+A(2,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(2,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1))))

    ! Calculate the inverse of the matrix
    B(1,1) = detinv*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2)))
    B(2,1) = detinv*(A(2,1)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))+A(2,3)*(A(3,1)*A(4,4)-A(3,4)*A(4,1))+A(2,4)*(A(3,3)*A(4,1)-A(3,1)*A(4,3)))
    B(3,1) = detinv*(A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))
    B(4,1) = detinv*(A(2,1)*(A(3,3)*A(4,2)-A(3,2)*A(4,3))+A(2,2)*(A(3,1)*A(4,3)-A(3,3)*A(4,1))+A(2,3)*(A(3,2)*A(4,1)-A(3,1)*A(4,2)))
    B(1,2) = detinv*(A(1,2)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))+A(1,3)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(1,4)*(A(3,3)*A(4,2)-A(3,2)*A(4,3)))
    B(2,2) = detinv*(A(1,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(1,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(1,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))
    B(3,2) = detinv*(A(1,1)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(1,2)*(A(3,1)*A(4,4)-A(3,4)*A(4,1))+A(1,4)*(A(3,2)*A(4,1)-A(3,1)*A(4,2)))
    B(4,2) = detinv*(A(1,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+A(1,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(1,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))
    B(1,3) = detinv*(A(1,2)*(A(2,3)*A(4,4)-A(2,4)*A(4,3))+A(1,3)*(A(2,4)*A(4,2)-A(2,2)*A(4,4))+A(1,4)*(A(2,2)*A(4,3)-A(2,3)*A(4,2)))
    B(2,3) = detinv*(A(1,1)*(A(2,4)*A(4,3)-A(2,3)*A(4,4))+A(1,3)*(A(2,1)*A(4,4)-A(2,4)*A(4,1))+A(1,4)*(A(2,3)*A(4,1)-A(2,1)*A(4,3)))
    B(3,3) = detinv*(A(1,1)*(A(2,2)*A(4,4)-A(2,4)*A(4,2))+A(1,2)*(A(2,4)*A(4,1)-A(2,1)*A(4,4))+A(1,4)*(A(2,1)*A(4,2)-A(2,2)*A(4,1)))
    B(4,3) = detinv*(A(1,1)*(A(2,3)*A(4,2)-A(2,2)*A(4,3))+A(1,2)*(A(2,1)*A(4,3)-A(2,3)*A(4,1))+A(1,3)*(A(2,2)*A(4,1)-A(2,1)*A(4,2)))
    B(1,4) = detinv*(A(1,2)*(A(2,4)*A(3,3)-A(2,3)*A(3,4))+A(1,3)*(A(2,2)*A(3,4)-A(2,4)*A(3,2))+A(1,4)*(A(2,3)*A(3,2)-A(2,2)*A(3,3)))
    B(2,4) = detinv*(A(1,1)*(A(2,3)*A(3,4)-A(2,4)*A(3,3))+A(1,3)*(A(2,4)*A(3,1)-A(2,1)*A(3,4))+A(1,4)*(A(2,1)*A(3,3)-A(2,3)*A(3,1)))
    B(3,4) = detinv*(A(1,1)*(A(2,4)*A(3,2)-A(2,2)*A(3,4))+A(1,2)*(A(2,1)*A(3,4)-A(2,4)*A(3,1))+A(1,4)*(A(2,2)*A(3,1)-A(2,1)*A(3,2)))
    B(4,4) = detinv*(A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))+A(1,2)*(A(2,3)*A(3,1)-A(2,1)*A(3,3))+A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1)))
    !
    return
    !
  end function matinv4
  !+-------------------------------------------------------------------+
  !| The end of the function matinv.                                   |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  ! Inverse matrix
  ! Method: Invert matrix by Gauss method
  !-----------------------------------------------------------
  ! input ...
  ! a(n,n) - array of coefficients for matrix A
  ! n      - dimension
  ! output ...
  ! aa(n,n) - inverse matrix of A
  ! comments ...
  !+-------------------------------------------------------------------+
  ! http://computer-programming-forum.com/49-fortran/55ed2ad2a9e4cef7.htm
  !+-------------------------------------------------------------------+
  function matinv (a,n) result(aa)  
    !
    real(8),intent(in) :: a(n,n)
    integer,intent(in) :: n
    real(8) :: aa(n,n)
    !
    ! - - - local variables - - -
    real(8) :: b(n,n), c, d, temp(n)
    integer :: i, j, k, m, imax(1), ipvt(n)
    ! - - - - - - - - - - - - - -
    b = a
    ipvt = (/ (i, i = 1, n) /)
    do k = 1,n
       imax = maxloc(abs(b(k:n,k)))
       m = k-1+imax(1)
       if (m /= k) then
          ipvt( (/m,k/) ) = ipvt( (/k,m/) )
          b((/m,k/),:) = b((/k,m/),:)
       end if
       d = 1/b(k,k)
       temp = b(:,k)
       do j = 1, n
          c = b(k,j)*d
          b(:,j) = b(:,j)-temp*c
          b(k,j) = c
       end do
       b(:,k) = temp*(-d)
       b(k,k) = d
    end do
    aa(:,ipvt) = b
    !
    return
    !
  end function matinv
  !+-------------------------------------------------------------------+
  !| The end of the function matinv.                                   |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to check if a matrix is an Identity matrix.    | 
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 18-08-2021  | Created by J. Fang @ Warrington                     |
  !+-------------------------------------------------------------------+
  pure logical function isidenmar(a,epsilon) result(iden)
    !
    ! arguments
    real(8),intent(in) :: a(:,:)
    real(8),intent(in) :: epsilon
    !
    ! local data
    integer :: n,i,j
    !
    n=size(a,1)
    !
    iden=.true.
    !
    do j=1,n
      if(abs(a(j,j)-1.d0)>epsilon) then
        iden=.false.
        return
      endif
    enddo
    !
    do j=1,n
    do i=1,n
      if(i==j) cycle
      if(abs(a(i,j))>epsilon) then
        iden=.false.
        return
      endif
    enddo
    enddo
    !
    return
    !
  end function isidenmar
  !+-------------------------------------------------------------------+
  !| The end of the function isidenmar.                                |
  !+-------------------------------------------------------------------+
  
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! this fuction is the 3-variable median() function.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    pure function median(var1,var2,var3)
      !
      real(8),intent(in) :: var1,var2,var3
      real(8) :: median
      !
      median= var1+minmod2(var2-var1,var3-var1)
      !
      return
      !
    end function median
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! End of the function median.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! this fuction is the 2-variable minmod() function.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    pure function minmod2(var1,var2)
      
      real(8),intent(in) :: var1,var2
      real(8) :: minmod2
      !
      if(var1>0.d0 .and. var2>0.d0) then
        minmod2=min(dabs(var1),dabs(var2))
      elseif(var1<0.d0 .and. var2<0.d0) then
        minmod2=-1.d0*min(dabs(var1),dabs(var2))
      else
        minmod2=0.d0
      end if
      !
      return
      !
    end function minmod2
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! End of the function minmod2.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! this fuction is the 4-variable minmod() function.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    pure function minmod4(var1,var2,var3,var4)
      !
      real(8),intent(in) :: var1,var2,var3,var4
      real(8) :: minmod4
      !
      if(var1>0.d0 .and. var2>0.d0 .and. var3>0.d0 .and. var4>0.d0) then
        minmod4=min(dabs(var1),dabs(var2),dabs(var3),dabs(var4))
      elseif(var1<0.d0 .and. var2<0.d0 .and. var3<0.d0 .and. var4<0.d0)&
                                                                    then
        minmod4=-1.d0*min(dabs(var1),dabs(var2),dabs(var3),dabs(var4))
      else
        minmod4=0.d0
      end if
      !
      return
      !
    end function minmod4
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! End of the function minmod4.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !+-------------------------------------------------------------------+
    !| This subroutine is a preprocessor of the Thomas algorithm.        |
    !+-------------------------------------------------------------------+
    !|ref: https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm    |
    !+-------------------------------------------------------------------+
    !| CHANGE RECORD                                                     |
    !| -------------                                                     |
    !| 27-05-2025  | Rewrite by J. Fang @ Liverpool.                     |
    !+-------------------------------------------------------------------+
    pure function tridiagonal_thomas_proprocess(a,c) result(ac)
        
        real(8),intent(in) :: a(:),c(:)

        real(8) :: ac(3,size(a))

        integer :: i

        ac(1,1)=c(1)

        do i=2,size(a)

            ac(1,i)=c(i)/(1.d0-a(i)*ac(1,i-1))

            ac(2,i)=1.d0/(1.d0-a(i)*ac(1,i-1))

            ac(3,i)=a(i)/(1.d0-a(i)*ac(1,i-1))

        enddo

        return

    end function tridiagonal_thomas_proprocess
    !+-------------------------------------------------------------------+
    !| The end of the subroutine tridiagonal_thomas_proprocess.          |
    !+-------------------------------------------------------------------+

    !+-------------------------------------------------------------------+
    !| This is a standard tridiagonal matrix algorithm, also known as the|
    !| Thomas algorithm to solve tridiagonal systems of equations.       |
    !+-------------------------------------------------------------------+
    !|ref: https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm    |
    !+-------------------------------------------------------------------+
    !| CHANGE RECORD                                                     |
    !| -------------                                                     |
    !| 25-05-2025  | Rewrite by J. Fang @ Liverpool.                     |
    !| 26-05-2025  | Separate the solver part by J. Fang @ Liverpool.    |
    !+-------------------------------------------------------------------+
    function tridiagonal_thomas_solver(ac,d) result(x)
        
        ! arguments
        real(8),intent(in) :: ac(:,:)
        real(8) :: d(:)
        real(8) :: x(size(d))

        ! local data 
        integer :: i,n

        n=size(d)

        do i=2,n
            d(i)=d(i)*ac(2,i)-d(i-1)*ac(3,i)
        enddo

        x(n)=d(n)
        do i=n-1,1,-1
            x(i)=d(i)-ac(1,i)*x(i+1)
        enddo

        return

    end function tridiagonal_thomas_solver
    !+-------------------------------------------------------------------+
    !| The end of the function tridiagonal_thomas_solver.                |
    !+-------------------------------------------------------------------+

    !+-------------------------------------------------------------------+
    !| This is a standard tridiagonal matrix algorithm, also known as the|
    !| Thomas algorithm to solve tridiagonal systems of equations.       |
    !+-------------------------------------------------------------------+
    !|ref: https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm    |
    !+-------------------------------------------------------------------+
    !| CHANGE RECORD                                                     |
    !| -------------                                                     |
    !| 25-05-2025  | Rewrite by J. Fang @ Liverpool.                     |
    !+-------------------------------------------------------------------+
    function tridiagonal_thomas(a,c,d) result(x)
        
        ! arguments
        real(8),intent(in) :: a(:),c(:),d(:)
        real(8) :: x(size(d))

        ! local data 
        integer :: i,n
        real(8),allocatable :: cc(:),dd(:)

        n=size(d)

        allocate(cc(n),dd(n))

        cc(1)=c(1)
        dd(1)=d(1)
        do i=2,n

            cc(i)=c(i)/(1.d0-a(i)*cc(i-1))

            dd(i)=(d(i)-a(i)*dd(i-1))/(1.d0-a(i)*cc(i-1))

        enddo

        x(n)=dd(n)
        do i=n-1,1,-1
            x(i)=dd(i)-cc(i)*x(i+1)
        enddo

        deallocate(cc,dd)

        return

    end function tridiagonal_thomas
    !+-------------------------------------------------------------------+
    !| The end of the function tridiagonal_thomas.                       |
    !+-------------------------------------------------------------------+

end module commfunc
!+---------------------------------------------------------------------+
!| The end of the module commfunc.                                     |
!+---------------------------------------------------------------------+