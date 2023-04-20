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
  !+-------------------------------------------------------------------+
  !| The end of the function linear1d                                  |
  !+-------------------------------------------------------------------+
  !
  subroutine regularlinearinterp1d1v(x1,v11,x2,v21)
    !
    !
    real(8),intent(in) :: x1(:),v11(:),x2(:)
    real(8),intent(out) :: v21(:)
    !
    integer :: dim1,dim2
    integer :: i,j,k,i1,j1,k1,nper
    !
    dim1=size(x1)
    dim2=size(x2)
    !
    nper=0
    do i=1,dim2
      !
      do i1=2,dim1
        !
        if( x2(i)>=x1(i1-1) .and. x2(i)<=x1(i1) ) then
          !
          v21(i)=linear1d_s(x1(i1-1),x1(i1),v11(i1-1),v11(i1),x2(i))
          exit
          !
        elseif(x2(i)>=x1(dim1)) then
          v21(i)=v11(dim1)
        elseif(x2(i)<=x1(1)) then
          v21(i)=v11(1)
        end if
        !
      end do
      !
      nper=nper+1
      !
      write(*,'(1A1,A23,I4,A2,$)')char(13),                            &
                          '  ** Interpolating ... ',100*nper/(dim2+1),' %'
    end do
    !
    print*
    !
  end subroutine regularlinearinterp1d1v
  !
end module interp
!+---------------------------------------------------------------------+
!| The end of the module interp.                                       |
!+---------------------------------------------------------------------+