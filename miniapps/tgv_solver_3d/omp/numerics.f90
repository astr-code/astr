module numerics

  use constdef

  implicit none

  real(rtype),allocatable :: cci(:,:),ccj(:,:),cck(:,:)
  real(rtype),allocatable :: cdi(:,:),cdj(:,:),cdk(:,:)
  real(rtype),allocatable :: cfi(:,:),cfj(:,:),cfk(:,:)
  real(rtype),allocatable :: alfa_fdm(:),coef10i(:)
  real(rtype) :: compact_filer_coefficient=0.49_rtype

  contains

  subroutine fdm_solver_init

    use comvardef, only : im,jm,km

    call qtds_solver_init(cci,'cflux4',im)
    call qtds_solver_init(ccj,'cflux4',jm)
    call qtds_solver_init(cck,'cflux4',km)

    call qtds_solver_init(cdi,'cc6',im)
    call qtds_solver_init(cdj,'cc6',jm)
    call qtds_solver_init(cdk,'cc6',km)

  end subroutine fdm_solver_init
  
  ! this subroutine is a 1d finite-difference solver
  subroutine fdm_solver_1d(f,df,dir)

    use comvardef, only : hm
    
    real(rtype),intent(in) :: f(:,:)
    character(len=1),intent(in) :: dir
    real(rtype),intent(out) :: df(:,:)

    integer :: dim,nclo,n

    dim =size(f,1)-2*hm-1
    nclo=size(f,2)


    do n=1,nclo
      ! call diff2ec(f(:,n),dim,hm,df(:,n))
      ! call diff6ec(f(:,n),dim,hm,df(:,n))
      call diff6cc(f(:,n),dim,df(:,n),dir)
      ! call diff4ec(f(:,n),dim,hm,df(:,n))
    enddo

    return

  end subroutine fdm_solver_1d

  subroutine diff6cc(vin,dim,vout,dir)

    use comvardef, only : hm
    use constdef

    integer,intent(in) :: dim
    real(rtype),intent(in) :: vin(-hm:dim+hm)
    character(len=1),intent(in) :: dir
    real(rtype) :: vout(0:dim),b(1:dim)

    ! local data
    integer :: i

    ! b(0)=0.75_rtype* (vin(1)-vin(-1)) -              &
    !      0.15_rtype* (vin(2)-vin(-2))+               &
    !      num1d60*(vin(3)-vin(-3))
    do i=1,dim
      b(i)=num7d9* (vin(i+1)-vin(i-1))+          &
           num1d36*(vin(i+2)-vin(i-2))
    end do
    ! b(dim)=0.75_rtype* (vin(dim+1)-vin(dim-1))-      &
    !        0.15_rtype* (vin(dim+2)-vin(dim-2))+      &
    !       num1d60* (vin(dim+3)-vin(dim-3))

    if(dir=='i') then
      call qtds_solver(b,vout(1:dim),cdi,dim)
    elseif(dir=='j') then
      call qtds_solver(b,vout(1:dim),cdj,dim)
    elseif(dir=='k') then
      call qtds_solver(b,vout(1:dim),cdk,dim)
    endif

    vout(0)=vout(dim)

  end subroutine diff6cc

  subroutine diff6ec(vin,dim,n,vout,comptime)
    !
    integer,intent(in) :: dim,n
    real(rtype),intent(in) :: vin(-n:dim+n)
    real(rtype) :: vout(0:dim)
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
      vout(i)  =0.75_rtype *(vin(i+1)-vin(i-1))- &
                0.15_rtype *(vin(i+2)-vin(i-2))+ &
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

  subroutine diff2ec(vin,dim,n,vout)

    integer,intent(in) :: dim,n
    real(rtype),intent(in) :: vin(-n:dim+n)
    real(rtype) :: vout(0:dim)
    
    ! local data
    integer :: i


    do i=0,dim
      vout(i)  =0.5_rtype *(vin(i+1)-vin(i-1))
    enddo

  end subroutine diff2ec

  subroutine diff4ec(vin,dim,n,vout)

    integer,intent(in) :: dim,n
    real(rtype),intent(in) :: vin(-n:dim+n)
    real(rtype) :: vout(0:dim)
    
    ! local data
    integer :: i


    do i=0,dim
      vout(i)  =num2d3 *(vin(i+1)-vin(i-1)) - &
                num1d12*(vin(i+2)-vin(i-2))
    enddo

  end subroutine diff4ec

  subroutine filter_init

    use comvardef, only : im,jm,km

    allocate(coef10i(0:5))

    coef10i(0)=(193._rtype +126._rtype*compact_filer_coefficient)/512._rtype
    coef10i(1)=(105._rtype +302._rtype*compact_filer_coefficient)/512._rtype
    coef10i(2)=(-15._rtype + 30._rtype*compact_filer_coefficient)/128._rtype
    coef10i(3)=( 45._rtype - 90._rtype*compact_filer_coefficient)/1024._rtype
    coef10i(4)=( -5._rtype + 10._rtype*compact_filer_coefficient)/512._rtype
    coef10i(5)=(  1._rtype -  2._rtype*compact_filer_coefficient)/1024._rtype

    call qtds_solver_init(cfi,'cf10',im)
    call qtds_solver_init(cfj,'cf10',jm)
    call qtds_solver_init(cfk,'cf10',km)

  end subroutine filter_init
  
  subroutine low_pass_filter(f,ff,dir,dim)

    use comvardef, only : hm
    
    integer,intent(in) :: dim
    real(rtype),intent(in) :: f(:)
    character(len=1),intent(in) :: dir
    real(rtype),intent(out) :: ff(:)

    ! call filter10ec(f,dim,hm,ff)
    call filter10cc(f,dim,hm,ff,dir)

  end subroutine low_pass_filter

  subroutine filter10ec(vin,dim,n,vout)
    !
    integer,intent(in) :: dim,n
    real(rtype),intent(in) :: vin(-n:dim+n)
    real(rtype) :: vout(0:dim)
    !
    integer :: i
    !
    do i=0,dim
      !
      vout(i)=   0.376953125_rtype*(vin(i)  +vin(i))     &
               + 0.205078125_rtype*(vin(i-1)+vin(i+1))   &
               -   0.1171875_rtype*(vin(i-2)+vin(i+2))   &
               +0.0439453125_rtype*(vin(i-3)+vin(i+3))   &
               - 0.009765625_rtype*(vin(i-4)+vin(i+4))   &
               +0.0009765625_rtype*(vin(i-5)+vin(i+5))
    enddo
    !
  end subroutine filter10ec

  subroutine filter10cc(vin,dim,n,vout,dir)
    !
    integer,intent(in) :: dim,n
    character(len=1),intent(in) :: dir
    real(rtype),intent(in) :: vin(-n:dim+n)
    real(rtype),intent(inout) :: vout(0:dim)
    
    real(rtype) :: b(1:dim)
    real(rtype) :: var0,var1,var2,var3,var4,var5

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
    real(rtype),intent(out) :: cc(1:m,1:6)

    integer :: i

    call matrix_coef(cc(:,1),cc(:,2),scheme)
    !
    cc(1,3)=1._rtype
    do i=1,m-2
      cc(i,4)=cc(i,1)/cc(i,3)
      cc(i+1,3)=1._rtype-cc(i+1,2)*cc(i,4)
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
    cc(m,6)=1._rtype
    do i=1,m-1
      cc(m,6)=cc(m,6)-cc(i,6)*cc(i,5)
    enddo

  end subroutine qtds_init
  
  function qtds(b,cc) result(x)
    !
    ! arguments
    real(rtype),intent(in) :: b(:),cc(:,:)
    real(rtype) :: x(size(b))
    !
    ! local data
    integer :: m,i
    real(rtype),allocatable :: y(:)
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
    real(rtype),intent(out) :: c(:),d(:)
    !
    !
    if(scheme=='cu5') then
      !
      d(:)=0.5_rtype
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
    real(rtype),intent(in) :: af(1:3)
    real(rtype),allocatable,intent(out) :: cc(:,:)
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
    cc(1,0)=1._rtype
    cc(2,0)=0._rtype
    !
    cc(1,1)=1._rtype-af(3)*cc(2,0)
    cc(2,1)=af(3)/cc(1,1)
    !
    cc(1,2)=1._rtype-af(3)*cc(2,1)
    !
    do l=2,dim-3
      cc(2,l)=af(3)/cc(1,l)
      cc(1,l+1)=1._rtype-af(3)*cc(2,l)
      ! if(mpirank==0) print*,l,cc(2,l),cc(1,l+1)
    end do
    cc(2,dim-2)=af(3)/cc(1,dim-2)
    !
    cc(1,dim-1)=1._rtype-af(3)*cc(2,dim-2)
    cc(2,dim-1)=af(3)/cc(1,dim-1)
    !
    cc(1,dim)=1._rtype

    cc(1,:)=1._rtype/cc(1,:)
    !
    return
    !
  end subroutine ptds_ini

  function ptds_cal(bd,af,cc,dim) result(xd)
    !
    integer,intent(in) :: dim
    real(rtype),intent(in) :: af(3),bd(0:dim),cc(1:2,0:dim)
    real(rtype) :: xd(0:dim)
    
    ! local data
    integer :: l
    real(rtype),allocatable,dimension(:) :: yd

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
    real(rtype),allocatable :: alfa(:)
    !
    if(scheme=='cc6') then
      allocate(alfa(3))
      alfa(3)=num1d3
      alfa(2)=0.25_rtype
      alfa(1)=1._rtype
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

    real(rtype),intent(out),allocatable :: cc(:,:)
    real(rtype) :: fac

    integer :: i,ii

    allocate(cc(n,4))

    if(scheme=='cc6') then
      cc(:,3)=num1d3
      cc(:,4)=num1d3
    endif

    if(scheme=='cflux4') then
      cc(:,3)=0.25_rtype
      cc(:,4)=0.25_rtype
    endif

    if(scheme=='cf10') then
      cc(:,3)=compact_filer_coefficient
      cc(:,4)=compact_filer_coefficient
    endif


    cc(1,1) = -cc(1,3)
    cc(1,2) = 1._rtype
    do i = 2, n - 1, 1
        ii = i - 1
        fac = cc(i,3) / cc(ii,2)
        cc(i,2) = 1._rtype - (fac * cc(ii,4))
        cc(i,1) = -fac * cc(ii,1)
    end do
    cc(n,2) = 1._rtype

  end subroutine qtds_solver_init

  subroutine qtds_solver(b,x,cc,n)

    integer, intent(in) :: n
    real(rtype), intent(in) :: b(:)    ! b vector
    real(rtype), intent(inout) :: x(:) ! output
    real(rtype), intent(in) :: cc(:,:)   ! prestored matrix

    real(rtype), dimension(n) :: dl      ! lower-diagonal
    real(rtype), dimension(n) :: du      ! upper-diagonal
    real(rtype), dimension(n) :: maind,w ! work array

    integer :: i, ii
    real(rtype) :: fac

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