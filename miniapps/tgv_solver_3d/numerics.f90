
module numerics
  !
  use constdef, only : num1d60
  implicit none
  !
  contains
  !
  subroutine diff6ec(vin,dim,n,vout,comptime)
    !
    integer,intent(in) :: dim,n
    real(8),intent(in) :: vin(-n:dim+n)
    real(8) :: vout(0:dim)
    !
    real,intent(inout),optional :: comptime
    
    ! local data
    integer :: i
    real :: tstart,tfinish

    if(present(comptime)) then
        call cpu_time(tstart)
    endif
    !
    do i=0,dim
      vout(i)  =0.75d0 *(vin(i+1)-vin(i-1))- &
                0.15d0 *(vin(i+2)-vin(i-2))+ &
               num1d60 *(vin(i+3)-vin(i-3))
    enddo
    !
    if(present(comptime)) then
      call cpu_time(tfinish)
      !
      comptime=comptime+tfinish-tstart
    endif
    !
  end subroutine diff6ec
    !
  subroutine filter10ec(vin,dim,n,vout)
    !
    integer,intent(in) :: dim,n
    real(8),intent(in) :: vin(-n:dim+n)
    real(8) :: vout(0:dim)
    !
    integer :: i
    !
    do i=0,dim
      !
      vout(i)=   0.376953125d0*(vin(i)  +vin(i))     &
               + 0.205078125d0*(vin(i-1)+vin(i+1))   &
               -   0.1171875d0*(vin(i-2)+vin(i+2))   &
               +0.0439453125d0*(vin(i-3)+vin(i+3))   &
               - 0.009765625d0*(vin(i-4)+vin(i+4))   &
               +0.0009765625d0*(vin(i-5)+vin(i+5))
    enddo
    !
  end subroutine filter10ec
  !
end module numerics