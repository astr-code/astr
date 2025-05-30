!+---------------------------------------------------------------------+
!| This module contains subroutines to do interpolation.               |
!| ==============                                                      |
!| CHANGE RECORD                                                       |
!| -------------                                                       |
!| 22-Jul-2022  | Created by J. Fang @ Warrington                      |
!+---------------------------------------------------------------------+
module interp
  !
  implicit none
  !
  interface interlinear
    module procedure linear1d_s
    module procedure linear1d_a1
    module procedure linear1d_a2
    module procedure linear1d_arrayin
  end interface interlinear
  !
  contains
  !
  !+-------------------------------------------------------------------+
  !| This function is a linear interpolation function.                 |
  !+-------------------------------------------------------------------+
  function linear1d_s(xx1,xx2,yy1,yy2,xx) result(yy)
    !
    real(8),intent(in) :: xx1,xx2,yy1,yy2,xx
    real(8) :: yy
    !
    real(8) :: var1
    !
    var1=(yy2-yy1)/(xx2-xx1)
    yy=var1*(xx-xx1)+yy1
    !
    return
    !
  end function linear1d_s
  !
  function linear1d_a1(xx1,xx2,yy1,yy2,xx) result(yy)
    !
    real(8),intent(in) :: xx1,xx2,xx
    real(8),intent(in) ::  yy1(:),yy2(:)
    real(8) :: yy(1:size(yy1))
    !
    real(8) :: var1
    !
    var1=(xx-xx1)/(xx2-xx1)
    yy=(yy2-yy1)*var1+yy1
    !
    return
    !
  end function linear1d_a1
  !
  function linear1d_a2(xx1,xx2,yy1,yy2,xx) result(yy)
    !
    real(8),intent(in) :: xx1,xx2,xx
    real(8),intent(in) ::  yy1(:,:),yy2(:,:)
    real(8) :: yy(1:size(yy1,1),1:size(yy1,2))
    !
    real(8) :: var1
    !
    var1=(xx-xx1)/(xx2-xx1)
    yy=(yy2-yy1)*var1+yy1
    !
    return
    !
  end function linear1d_a2
  !
  function linear1d_arrayin(x1,y1,xx) result(yy)
    !
    real(8),intent(in) :: x1(:),y1(:),xx
    real(8) :: yy
    !
    integer :: dim,i
    !
    dim=size(x1)
    !
    if(xx<x1(1)) then
      ! yy=linear1d_s(x1(1),x1(2),y1(1),y1(2),xx)
      yy=2.d0*y1(1)-y1(2)
    elseif(xx>=x1(dim)) then
      ! yy=linear1d_s(x1(dim-1),x1(dim),y1(dim-1),y1(dim),xx)
      yy=2.d0*y1(dim)-y1(dim-1)
    else
      do i=2,dim
        if(xx>=x1(i-1) .and. xx<x1(i)) then
          yy=linear1d_s(x1(i-1),x1(i),y1(i-1),y1(i),xx)
        endif
      enddo
    endif
    !
    return
    !
  end function linear1d_arrayin
  !+-------------------------------------------------------------------+
  !| The end of the function linear1d                                  |
  !+-------------------------------------------------------------------+
  !
  !
end module interp
!+---------------------------------------------------------------------+
!| The end of the module interp.                                       |
!+---------------------------------------------------------------------+