module pastr_interpolation

    use iso_fortran_env, only: wp => real64

    implicit none

    interface interpolat
      module procedure linear1d_s
      module procedure linear1d_search_array
    end interface interpolat

contains

    pure function linear1d_s(xx1,xx2,yy1,yy2,xx) result(yy)
      real(wp),intent(in) :: xx1,xx2,yy1,yy2,xx
      real(wp) :: yy
      real(wp) :: var1
      var1=(yy2-yy1)/(xx2-xx1)
      yy=var1*(xx-xx1)+yy1
      return
    end function linear1d_s

    pure function linear1d_search_array(x,y,xx) result(yy)

      real(wp),intent(in) :: x(:),y(:),xx
      real(wp) :: yy

      integer :: i,m

      m=size(x)

      if(xx<=x(1)) then
        yy=linear1d_s(x(1),x(2),y(1),y(2),xx)
      elseif(xx>=x(m)) then
        yy=linear1d_s(x(m-1),x(m),y(m-1),y(m),xx)
      else
        do i=2,m
          if(xx>=x(i-1) .and. xx<=x(i)) then
            yy=linear1d_s(x(i-1),x(i),y(i-1),y(i),xx)
            exit
          endif
        enddo
      endif

      return

   end function linear1d_search_array

end module pastr_interpolation