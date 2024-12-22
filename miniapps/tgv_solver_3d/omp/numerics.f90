
module numerics
  !
  use constdef, only : num1d60
  implicit none

  real(8),allocatable :: cci(:,:),ccj(:,:),cck(:,:)
  real(8),allocatable :: cfi(:,:),cfj(:,:),cfk(:,:)
  real(8),allocatable :: alfa_fdm(:),coef10i(:)
  real(8) :: compact_filer_coefficient=0.49

  contains

  subroutine fdm_solver_init

    use comvardef, only : im,jm,km

    ! alfa_fdm=coeffcompac('cc6')

    ! call ptds_ini(cci,alfa_fdm,im)
    ! call ptds_ini(ccj,alfa_fdm,jm)
    ! call ptds_ini(cck,alfa_fdm,km)

    call qtds_solver_init(cci,'cc6',im)
    call qtds_solver_init(ccj,'cc6',jm)
    call qtds_solver_init(cck,'cc6',km)

  end subroutine fdm_solver_init
  
  ! this subroutine is a 1d finite-difference solver
  subroutine fdm_solver_1d(f,df,dir)

    use comvardef, only : hm,im,jm,km
    
    real(8),intent(in) :: f(:,:)
    character(len=1),intent(in) :: dir
    real(8),intent(out) :: df(:,:)

    integer :: dim,nclo,n

    dim =size(f,1)-2*hm-1
    nclo=size(f,2)

    do n=1,nclo
      ! call diff6ec(f(:,n),dim,hm,df(:,n))
      call diff6cc(f(:,n),dim,df(:,n),dir)
    enddo

    return

  end subroutine fdm_solver_1d

  subroutine diff6cc(vin,dim,vout,dir)

    use comvardef, only : hm
    use constdef

    integer,intent(in) :: dim
    real(8),intent(in) :: vin(-hm:dim+hm)
    character(len=1),intent(in) :: dir
    real(8) :: vout(0:dim),b(1:dim)

    ! local data
    integer :: i

    ! b(0)=0.75d0* (vin(1)-vin(-1)) -              &
    !      0.15d0* (vin(2)-vin(-2))+               &
    !      num1d60*(vin(3)-vin(-3))
    do i=1,dim
      b(i)=num7d9* (vin(i+1)-vin(i-1))+          &
           num1d36*(vin(i+2)-vin(i-2))
    end do
    ! b(dim)=0.75d0* (vin(dim+1)-vin(dim-1))-      &
    !        0.15d0* (vin(dim+2)-vin(dim-2))+      &
    !       num1d60* (vin(dim+3)-vin(dim-3))

    if(dir=='i') then
      call qtds_solver(b,vout(1:dim),cci,dim)
    elseif(dir=='j') then
      call qtds_solver(b,vout(1:dim),ccj,dim)
    elseif(dir=='k') then
      call qtds_solver(b,vout(1:dim),cck,dim)
    endif

    vout(0)=vout(dim)

  end subroutine diff6cc

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

  subroutine filter_init

    use comvardef, only : im,jm,km

    allocate(coef10i(0:5))

    coef10i(0)=(193.d0 +126.d0*compact_filer_coefficient)/512.d0
    coef10i(1)=(105.d0 +302.d0*compact_filer_coefficient)/512.d0
    coef10i(2)=(-15.d0 + 30.d0*compact_filer_coefficient)/128.d0
    coef10i(3)=( 45.d0 - 90.d0*compact_filer_coefficient)/1024.d0
    coef10i(4)=( -5.d0 + 10.d0*compact_filer_coefficient)/512.d0
    coef10i(5)=(  1.d0 -  2.d0*compact_filer_coefficient)/1024.d0

    call qtds_solver_init(cfi,'cf10',im)
    call qtds_solver_init(cfj,'cf10',jm)
    call qtds_solver_init(cfk,'cf10',km)

  end subroutine filter_init
  
  subroutine low_pass_filter(f,ff,dir,dim)

    use comvardef, only : hm,im,jm,km
    
    integer,intent(in) :: dim
    real(8),intent(in) :: f(:)
    character(len=1),intent(in) :: dir
    real(8),intent(out) :: ff(:)

    integer :: nclo,n

    ! call filter10ec(f,dim,hm,ff)
    call filter10cc(f,dim,hm,ff,dir)

  end subroutine low_pass_filter

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

  subroutine filter10cc(vin,dim,n,vout,dir)
    !
    integer,intent(in) :: dim,n
    character(len=1),intent(in) :: dir
    real(8),intent(in) :: vin(-n:dim+n)
    real(8),intent(inout) :: vout(0:dim)
    
    real(8) :: b(1:dim)
    real(8) :: var0,var1,var2,var3,var4,var5

    integer :: i
    !
    do i=1,dim
      !
      var0=vin(i)  +vin(i)
      var1=vin(i+1)+vin(i-1)
      var2=vin(i+2)+vin(i-2)
      var3=vin(i+3)+vin(i-3)
      var4=vin(i+4)+vin(i-4)
      var5=vin(i+5)+vin(i-5)
      !
      b(i)=coef10i(0)*var0+coef10i(1)*var1+coef10i(2)*var2+          &
           coef10i(3)*var3+coef10i(4)*var4+coef10i(5)*var5
      !
    end do

    if(dir=='i') then
      call qtds_solver(b,vout(1:dim),cfi,dim)
    elseif(dir=='j') then
      call qtds_solver(b,vout(1:dim),cfj,dim)
    elseif(dir=='k') then
      call qtds_solver(b,vout(1:dim),cfk,dim)
    endif

    vout(0)=vout(dim)

  end subroutine filter10cc
  !
  subroutine qtds_init(scheme,cc,m)

    character(len=*),intent(in) :: scheme
    integer,intent(in) :: m
    real(8),intent(out) :: cc(1:m,1:6)

    integer :: i

    call matrix_coef(cc(:,1),cc(:,2),scheme)
    !
    cc(1,3)=1.d0
    do i=1,m-2
      cc(i,4)=cc(i,1)/cc(i,3)
      cc(i+1,3)=1.d0-cc(i+1,2)*cc(i,4)
    enddo
    !
    cc(1,5)=cc(1,2)/cc(1,3)
    do i=2,m-2
      cc(i,5)=-cc(i,2)*cc(i-1,5)/cc(i,3)
    enddo
    cc(m-1,5)=(cc(m-1,1)-cc(m-1,2)*cc(m-2,5))/cc(m-1,3)
    !
    cc(1,6)=cc(m,1)
    do i=2,m-2
      cc(i,6)=-cc(i-1,6)*cc(i-1,4)
    enddo
    cc(m-1,6)=cc(m,2)-cc(m-2,6)*cc(m-2,4)
    !
    cc(m,6)=1.d0
    do i=1,m-1
      cc(m,6)=cc(m,6)-cc(i,6)*cc(i,5)
    enddo

  end subroutine qtds_init
  
  function qtds(b,cc) result(x)
    !
    ! arguments
    real(8),intent(in) :: b(:),cc(:,:)
    real(8) :: x(size(b))
    !
    ! local data
    integer :: m,i
    real(8),allocatable :: y(:)
    !
    m=size(b)
    !
    allocate(y(1:m))

    y(1)=b(1)/cc(1,3)
    do i=2,m-1
      y(i)=(b(i)-cc(i,2)*y(i-1))/cc(i,3)
    enddo
    !
    y(m)=b(m)
    do i=1,m-1
      y(m)=y(m)-cc(i,6)*y(i)
    enddo
    y(m)=y(m)/cc(m,6)
    !
    x(m)=y(m)
    x(m-1)=y(m-1)-cc(m-1,5)*x(m)
    do i=m-2,1,-1
      x(i)=y(i)-cc(i,4)*x(i+1)-cc(i,5)*x(m)
    enddo
    
    deallocate(y)

    return
    !
  end function qtds

  !+-------------------------------------------------------------------+
  !| This subroutine is to set A matrix.                               |
  !+-------------------------------------------------------------------+
  subroutine matrix_coef(c,d,scheme)
    !
    use constdef
    !
    character(len=*),intent(in) :: scheme
    real(8),intent(out) :: c(:),d(:)
    !
    !
    if(scheme=='cu5') then
      !
      d(:)=0.5d0
      c(:)=num1d6
      !
    elseif(scheme=='cc6') then
      !
      d(:)=num1d3
      c(:)=num1d3
      !
    else
      stop ' !! scheme not defined @ qtds_a_init'
    endif
    !
    return
    !
  end subroutine matrix_coef
  !+-------------------------------------------------------------------+
  !| The end of the subroutine qtds  system                            |
  !+-------------------------------------------------------------------+
  !
  subroutine ptds_ini(cc,af,dim)
    !
    integer,intent(in) :: dim
    real(8),intent(in) :: af(1:3)
    real(8),allocatable,intent(out) :: cc(:,:)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! af(1),af2,af: input dat
    ! cc: output array
    ! dim: input dat
    ! l: temporary variable
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! local data
    integer :: l
    !
    allocate(cc(1:2,0:dim))
    !
    cc(1,0)=1.d0
    cc(2,0)=0.d0
    !
    cc(1,1)=1.d0-af(3)*cc(2,0)
    cc(2,1)=af(3)/cc(1,1)
    !
    cc(1,2)=1.d0-af(3)*cc(2,1)
    !
    do l=2,dim-3
      cc(2,l)=af(3)/cc(1,l)
      cc(1,l+1)=1.d0-af(3)*cc(2,l)
      ! if(mpirank==0) print*,l,cc(2,l),cc(1,l+1)
    end do
    cc(2,dim-2)=af(3)/cc(1,dim-2)
    !
    cc(1,dim-1)=1.d0-af(3)*cc(2,dim-2)
    cc(2,dim-1)=af(3)/cc(1,dim-1)
    !
    cc(1,dim)=1.d0

    cc(1,:)=1.d0/cc(1,:)
    !
    return
    !
  end subroutine ptds_ini

  function ptds_cal(bd,af,cc,dim) result(xd)
    !
    integer,intent(in) :: dim
    real(8),intent(in) :: af(3),bd(0:dim),cc(1:2,0:dim)
    real(8) :: xd(0:dim)
    
    ! local data
    integer :: l
    real(8),allocatable,dimension(:) :: yd

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! af(3): input dat
    ! bd: input array
    ! xd: output array
    ! cc: input array
    ! dim: input dat
    ! l, yd: temporary variable
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    allocate ( yd(0:dim) )

    ! inner block
    yd(0)=bd(0)*cc(1,0)
    do l=1,dim-1
      yd(l)=(bd(l)-af(3)*yd(l-1))*cc(1,l)
    end do
    yd(dim)=bd(dim)*cc(1,dim)

    xd(dim)=yd(dim)
    do l=dim-1,0,-1
      xd(l)=yd(l)-cc(2,l)*xd(l+1)
    end do

    deallocate( yd )

    return
    !
  end function ptds_cal

  !+-------------------------------------------------------------------+
  !| This function is to return the coefficients of a compact scheme.  |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 08-02-2021  | Created by J. Fang @ Warrington.                    |
  !+-------------------------------------------------------------------+
  function coeffcompac(scheme) result(alfa)
    
    use constdef
    
    character(len=*),intent(in) :: scheme
    real(8),allocatable :: alfa(:)
    !
    if(scheme=='cc6') then
      allocate(alfa(3))
      alfa(3)=num1d3
      alfa(2)=0.25d0
      alfa(1)=1.d0
    else
      stop ' !! scheme not defined @ coef_diffcompac'
    endif
    !
    return
    !
  end function coeffcompac
  !+-------------------------------------------------------------------+
  !| The end of the function coeffcompac.                              |
  !+-------------------------------------------------------------------+
  !
  ! https://www.reddit.com/r/fortran/comments/l5doqt/optimizing_solver_for_almost_tridiagonal_matrix/
  ! http://inis.jinr.ru/sl/Simulation/Hirsch,_Numerical_Computation_of_Internal&External_Flows,1994/Hirsch,_Numerical_Computation_of_Internal&External_Flows,v1,1994/conclusion.pdf
  subroutine qtds_solver_init(cc,scheme,n)

    use constdef

    character(len=*),intent(in) :: scheme
    integer,intent(in) :: n

    real(8),intent(out),allocatable :: cc(:,:)
    real(8) :: fac

    integer :: i,ii

    allocate(cc(n,4))

    if(scheme=='cc6') then
      cc(:,3)=num1d3
      cc(:,4)=num1d3
    endif

    if(scheme=='cf10') then
      cc(:,3)=compact_filer_coefficient
      cc(:,4)=compact_filer_coefficient
    endif


    cc(1,1) = -cc(1,3)
    cc(1,2) = 1.d0
    do i = 2, n - 1, 1
        ii = i - 1
        fac = cc(i,3) / cc(ii,2)
        cc(i,2) = 1.d0 - (fac * cc(ii,4))
        cc(i,1) = -fac * cc(ii,1)
    end do
    cc(n,2) = 1.d0

  end subroutine qtds_solver_init

  subroutine qtds_solver(b,x,cc,n)

    integer, intent(in) :: n
    real(8), intent(in) :: b(:)    ! b vector
    real(8), intent(inout) :: x(:) ! output
    real(8), intent(in) :: cc(:,:)   ! prestored matrix

    real(8), dimension(n) :: dl      ! lower-diagonal
    real(8), dimension(n) :: du      ! upper-diagonal
    real(8), dimension(n) :: maind,w ! work array

    integer :: i, ii
    real(8) :: fac

    w(1:n)    =cc(1:n,1)
    maind(1:n)=cc(1:n,2)
    dl(1:n)   =cc(1:n,3)
    du(1:n)   =cc(1:n,4)

    x(1) = b(1)
    do i = 2, n - 1, 1
        ii = i - 1
        fac = dl(i) / maind(ii)
        x(i) = b(i) - (fac * x(ii))
    end do
    x(n) = b(n)

    ii = n - 1
    x(ii) = x(ii) / maind(ii)
    w(ii) = (w(ii) - du(ii)) / maind(ii)

    do i = n - 2, 1, -1
        ii = i + 1
        x(i) = (x(i) - du(i) * x(ii)) / maind(i)
        w(i) = (w(i) - du(i) * w(ii)) / maind(i)
    end do

    i = n
    ii = n - 1
    fac = maind(i) + (du(i) * w(1)) + (dl(i) * w(ii))
    x(i) = (x(i) - ((du(i) * x(1)) + (dl(i) * x(ii)))) / fac

    fac = x(n)
    do i = 1, n - 1, 1
        x(i) = x(i) + (w(i) * fac)
    end do

  end subroutine qtds_solver

end module numerics