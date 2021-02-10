!+---------------------------------------------------------------------+
!| This module contains subroutine of common functions.                |
!+---------------------------------------------------------------------+
!| CHANGE RECORD                                                       |
!| -------------                                                       |
!| 08-02-2021  | Created by J. Fang                                    |
!+---------------------------------------------------------------------+
module commfunc
  !
  use commvar, only : hm
  use parallel, only: mpirank,mpistop
  use constdef
  !
  implicit none
  !
  real(8),allocatable :: coef2i(:),coef4i(:),coef6i(:),coef8i(:),     &
                         coef10i(:),coefb(:,:)
  !
  contains
  !
  !+-------------------------------------------------------------------+
  !| This function is to return the finite-difference value of a input |
  !| function.                                                         |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 08-02-2021  | Created by J. Fang @ Warrington.                    |
  !+-------------------------------------------------------------------+
  function ddfc(f,stype,ntype,dim,af,cc)
    !
    ! arguments
    character(len=4),intent(in) :: stype
    integer,intent(in) :: ntype,dim
    real(8),intent(in) :: f(-hm:dim+hm)
    real(8),intent(in),optional :: af(3),cc(1:2,0:dim)
    real(8) :: ddfc(0:dim)
    !
    ! local data
    integer :: nscheme
    real(8) :: b(0:dim)
    !
    ! print*,mpirank,'|',dim,im
    !
    read(stype(1:3),*) nscheme
    !
    if(dim==0) then
      ddfc=f(1)-f(0)
    else
      !
      b   =ptds_rhs(f,dim,nscheme,ntype)
      ddfc=ptds_cal(b,af,cc,dim,ntype)
      !
    endif
    !
    return
    !
  end function ddfc
  !+-------------------------------------------------------------------+
  !| The end of the function ddf.                                      |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This function is to return the coefficients of a compact scheme.  |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 08-02-2021  | Created by J. Fang @ Warrington.                    |
  !+-------------------------------------------------------------------+
  function coeffcompac(scheme) result(alfa)
    !
    integer,intent(in) :: scheme
    real(8) :: alfa(1:3)
    !
    if(scheme==642) then
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
  !+-------------------------------------------------------------------+
  !| This function is to apply compact low-pass filter to a input array|
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 10-02-2021  | Created by J. Fang @ Warrington.                    |
  !+-------------------------------------------------------------------+
  real(8) function filter8exp(a)
    !
    real(8) :: a(-4:4)
    !
    real(8) :: coef(0:4)
    integer :: n
    !
    coef(0)= 93.d0
    coef(1)= 56.d0
    coef(2)=-28.d0
    coef(3)=  8.d0
    coef(4)= -1.d0
    !
    filter8exp=0.d0
    do n=0,4
      filter8exp=filter8exp+coef(n)/256.d0*(a(n)+a(-n))
    enddo 
    !
  end function filter8exp
  !
  function spafilter10(f,ntype,dim,af,fc) result(ff)
    !
    ! arguments
    integer,intent(in) :: ntype,dim
    real(8),intent(in) :: f(-hm:dim+hm)
    real(8),intent(in) :: af,fc(1:2,-hm:dim+hm)
    real(8) :: ff(0:dim)
    !
    ! local data
    real(8),allocatable :: b(:)
    integer :: i
    real(8) :: aff(3)
    !
    aff(1)=0.d0
    aff(2)=af
    aff(3)=af
    ! b =pfilterrhs(f,dim,ntype)
    ! ff=ptdsfilter_cal(b,af,fc,dim,ntype)
    b=PFilterRHS2(f,dim,ntype)
    ff=ptds_cal(b,aff,fc,dim,ntype)
    !
    ! do i=0,dim
    !     ff(i)=filter8exp(f(i-4:i+4))
    ! enddo
    !
    return
    !
  end function spafilter10
  !+-------------------------------------------------------------------+
  !| The end of the function spafilter10.                              |
  !+-------------------------------------------------------------------+
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! This subroutine is used to assign the right-hand side of the spatial
  ! filter in parallel environment, where the boundary are not filtered.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Writen by Fang Jian, 2010-09-07.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function pfilterrhs(var,dim,ntype) result(b)
    !
    ! arguments
    real(8),intent(in) :: var(-hm:dim+hm)
    integer,intent(in) :: dim,ntype
    real(8) :: b(-hm:dim+hm)
    !
    ! local data
    integer :: i,m
    real(8) :: var0,var1,var2,var3,var4,var5
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! dim: dimensions in the direction.
    ! var: the variable of the filetering.
    ! var* : temporary variable.
    ! b: return variable.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    if(ntype==1) then
      !
      ! boundary filter
      do i=0,3
        b(i)=0.d0
        do m=0,8
          b(i)=b(i)+coefb(i,m)*var(m)
        enddo
      enddo
      !
      ! 8th-order filter
      i=4
      var0=var(i)+var(i)
      var1=var(i+1)+var(i-1)
      var2=var(i+2)+var(i-2)
      var3=var(i+3)+var(i-3)
      var4=var(i+4)+var(i-4)
      !
      b(i)=coef8i(0)*var0+coef8i(1)*var1+coef8i(2)*var2+               &
           coef8i(3)*var3+coef8i(4)*var4
      !
      ! inner block
      do i=5,dim
        !
        var0=var(i)+var(i)
        var1=var(i+1)+var(i-1)
        var2=var(i+2)+var(i-2)
        var3=var(i+3)+var(i-3)
        var4=var(i+4)+var(i-4)
        var5=var(i+5)+var(i-5)
        !
        b(i)=coef10i(0)*var0+coef10i(1)*var1+coef10i(2)*var2+          &
             coef10i(3)*var3+coef10i(4)*var4+coef10i(5)*var5
        !
      end do
      !
      ! 8th-order filter
      i=dim+1
      var0=var(i)+var(i)
      var1=var(i+1)+var(i-1)
      var2=var(i+2)+var(i-2)
      var3=var(i+3)+var(i-3)
      var4=var(i+4)+var(i-4)
      !
      b(i)=coef8i(0)*var0+coef8i(1)*var1+coef8i(2)*var2+               &
           coef8i(3)*var3+coef8i(4)*var4
      !
      ! boundary filter
      do i=dim+2,dim+hm
        b(i)=0.d0
        do m=0,8
          b(i)=b(i)+coefb(dim+hm-i,m)*var(dim+hm-m)
        enddo
      enddo
      !
    elseif(ntype==2) then
      !
      ! boundary filter
      do i=-hm,-2
        b(i)=0.d0
        do m=0,8
          b(i)=b(i)+coefb(i+hm,m)*var(-hm+m)
        enddo
      enddo
      !
      ! 8th-order filter
      i=-1
      var0=var(i)+var(i)
      var1=var(i+1)+var(i-1)
      var2=var(i+2)+var(i-2)
      var3=var(i+3)+var(i-3)
      var4=var(i+4)+var(i-4)
      !
      b(i)=coef8i(0)*var0+coef8i(1)*var1+coef8i(2)*var2+               &
           coef8i(3)*var3+coef8i(4)*var4

      ! inner block
      do i=0,dim-5
        !
        var0=var(i)+var(i)
        var1=var(i+1)+var(i-1)
        var2=var(i+2)+var(i-2)
        var3=var(i+3)+var(i-3)
        var4=var(i+4)+var(i-4)
        var5=var(i+5)+var(i-5)
        !
        b(i)=coef10i(0)*var0+coef10i(1)*var1+coef10i(2)*var2+          &
             coef10i(3)*var3+coef10i(4)*var4+coef10i(5)*var5
        !
      end do
      !
      ! 8th-order filter
      i=dim-4
      var0=var(i)+var(i)
      var1=var(i+1)+var(i-1)
      var2=var(i+2)+var(i-2)
      var3=var(i+3)+var(i-3)
      var4=var(i+4)+var(i-4)
      !
      b(i)=coef8i(0)*var0+coef8i(1)*var1+coef8i(2)*var2+               &
           coef8i(3)*var3+coef8i(4)*var4
      !
      ! boundary filter
      do i=dim-3,dim
        b(i)=0.d0
        do m=0,8
          b(i)=b(i)+coefb(dim-i,m)*var(dim-m)
        enddo
      enddo
      !
    elseif(ntype==3) then
      !
      ! boundary filter
      do i=-hm,-2
        b(i)=0.d0
        do m=0,8
          b(i)=b(i)+coefb(i+hm,m)*var(-hm+m)
        enddo
      enddo
      !
      ! 8th-order filter
      i=-1
      var0=var(i)+var(i)
      var1=var(i+1)+var(i-1)
      var2=var(i+2)+var(i-2)
      var3=var(i+3)+var(i-3)
      var4=var(i+4)+var(i-4)
      !
      b(i)=coef8i(0)*var0+coef8i(1)*var1+coef8i(2)*var2+               &
           coef8i(3)*var3+coef8i(4)*var4

      ! inner block
      do i=0,dim
        !
        var0=var(i)+var(i)
        var1=var(i+1)+var(i-1)
        var2=var(i+2)+var(i-2)
        var3=var(i+3)+var(i-3)
        var4=var(i+4)+var(i-4)
        var5=var(i+5)+var(i-5)
        !
        b(i)=coef10i(0)*var0+coef10i(1)*var1+coef10i(2)*var2+          &
             coef10i(3)*var3+coef10i(4)*var4+coef10i(5)*var5
        !
      end do
      !
      ! 8th-order filter
      i=dim+1
      var0=var(i)+var(i)
      var1=var(i+1)+var(i-1)
      var2=var(i+2)+var(i-2)
      var3=var(i+3)+var(i-3)
      var4=var(i+4)+var(i-4)
      !
      b(i)=coef8i(0)*var0+coef8i(1)*var1+coef8i(2)*var2+               &
           coef8i(3)*var3+coef8i(4)*var4
      !
      ! boundary filter
      do i=dim+2,dim+hm
        b(i)=0.d0
        do m=0,8
          b(i)=b(i)+coefb(dim+hm-i,m)*var(dim+hm-m)
        enddo
      enddo
      !
    else
      print*,' !! ntype error @ function pfilterrhs!'
    end if
    !
    !
  end function pfilterrhs
  !
  function PFilterRHS2(var,dim,ntype) result(b)
    !
    
    integer :: dim,ntype
    integer :: l
    real(8) :: var(-hm:dim+hm),b(0:dim)
    real(8) :: var0,var1,var2,var3,var4,var5
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! dim: dimensions in the direction.
    ! var: the variable of the filetering.
    ! var* : temporary variable.
    ! b: return variable.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    if(ntype==3) then
      !
      b(0)=    0.376953125d0*(var(0)  +var(0))   &
             + 0.205078125d0*(var(-1) +var(1))   &
             -   0.1171875d0*(var(-2) +var(2))   &
             +0.0439453125d0*(var(-3) +var(3))   &
             - 0.009765625d0*(var(-4) +var(4))   &
             +0.0009765625d0*(var(-5) +var(5))
    
      ! inner block
      do l=1,dim-1
        !
        var0=var(l)+var(l)
        var1=var(l+1)+var(l-1)
        var2=var(l+2)+var(l-2)
        var3=var(l+3)+var(l-3)
        var4=var(l+4)+var(l-4)
        var5=var(l+5)+var(l-5)
        !
        b(l)=coef10i(0)*var0+coef10i(1)*var1+coef10i(2)*var2+          &
             coef10i(3)*var3+coef10i(4)*var4+coef10i(5)*var5
        !
      end do
      !
      b(dim)=  0.376953125d0*(var(dim)  +var(dim))     &
             + 0.205078125d0*(var(dim-1)+var(dim+1))   &
             -   0.1171875d0*(var(dim-2)+var(dim+2))   &
             +0.0439453125d0*(var(dim-3)+var(dim+3))   &
             - 0.009765625d0*(var(dim-4)+var(dim+4))   &
             +0.0009765625d0*(var(dim-5)+var(dim+5))
    
    else
      print*,' !! Error in subroutine PFilterRHS!'
    end if
    !
    !
  end function PFilterRHS2
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! End of the function pfilterrhs
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! This subroutine is used to generate the coefficient for the filter,
  ! and bounday filter is included.
  ! The boundary filter is lower order (0-6-6-6-8...) using one side
  ! filter.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Writen by Fang Jian, 2008-11-05.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine genfilt10coef(alfa)
    !
    real(8),intent(in) :: alfa
    !
    allocate(coef2i(0:1),coef4i(0:2),coef6i(0:3),coef8i(0:4),coef10i(0:5))
    allocate(coefb(0:4,0:8))
    !
    ! inernal coefficent
    ! 2nd-order
    coef2i(0)=( 1.d0 + 2.d0*alfa)   /4.d0
    coef2i(1)=( 1.d0 + 2.d0*alfa)   /4.d0
    !
    ! 4th-order
    coef4i(0)=( 5.d0 + 6.d0*alfa)   /16.d0
    coef4i(1)=( 1.d0 + 2.d0*alfa)   /4.d0
    coef4i(2)=(-1.d0 + 2.d0*alfa)   /16.d0
    !
    ! 6th-order
    coef6i(0)=(11.d0 +10.d0*alfa)   /32.d0
    coef6i(1)=(15.d0 +34.d0*alfa)   /64.d0
    coef6i(2)=(-3.d0 + 6.d0*alfa)   /32.d0
    coef6i(3)=( 1.d0 - 2.d0*alfa)   /64.d0
    !
    ! 8th-order
    coef8i(0)=(93.d0 +70.d0*alfa)   /256.d0
    coef8i(1)=( 7.d0 +18.d0*alfa)   /32.d0
    coef8i(2)=(-7.d0 +14.d0*alfa)   /64.d0
    coef8i(3)=( 1.d0 - 2.d0*alfa)   /32.d0
    coef8i(4)=(-1.d0 + 2.d0*alfa)   /256.d0
    !
    ! 10th-order
    coef10i(0)=(193.d0 +126.d0*alfa)/512.d0
    coef10i(1)=(105.d0 +302.d0*alfa)/512.d0
    coef10i(2)=(-15.d0 + 30.d0*alfa)/128.d0
    coef10i(3)=( 45.d0 - 90.d0*alfa)/1024.d0
    coef10i(4)=( -5.d0 + 10.d0*alfa)/512.d0
    coef10i(5)=(  1.d0 -  2.d0*alfa)/1024.d0
    !
    ! boundary coef10ifficent
    !
    ! 8th-order for the point3
    coefb(3,0)=(  1.d0 - 2.d0*alfa)  /256.d0
    coefb(3,1)=( -1.d0 + 2.d0*alfa)  /32.d0
    coefb(3,2)=(  7.d0 +50.d0*alfa)  /64.d0
    coefb(3,3)=( 25.d0 +14.d0*alfa)  /32.d0
    coefb(3,4)=( 35.d0 +58.d0*alfa)  /128.d0
    coefb(3,5)=( -7.d0 +14.d0*alfa)  /32.d0
    coefb(3,6)=(  7.d0 -14.d0*alfa)  /64.d0
    coefb(3,7)=( -1.d0 + 2.d0*alfa)  /32.d0
    coefb(3,8)=(  1.d0 - 2.d0*alfa)  /256.d0
    !
    ! 6-order asymmetry scheme for point 32
    coefb(2,0)=( -1.d0 + 2.d0*alfa)  /64.d0
    coefb(2,1)=(  3.d0 +26.d0*alfa)  /32.d0
    coefb(2,2)=( 49.d0 +30.d0*alfa)  /64.d0
    coefb(2,3)=(  5.d0 + 6.d0*alfa)  /16.d0
    coefb(2,4)=(-15.d0 +30.d0*alfa)  /64.d0
    coefb(2,5)=(  3.d0 - 6.d0*alfa)  /32.d0
    coefb(2,6)=( -1.d0 + 2.d0*alfa)  /64.d0
    coefb(2,7)=0.d0
    coefb(2,8)=0.d0
    !
    ! 6-order asymmetry scheme for point 1
    coefb(1,0)=( 1.d0 +62.d0*alfa)  /64.d0
    coefb(1,1)=(29.d0 + 6.d0*alfa)  /32.d0
    coefb(1,2)=(15.d0 +34.d0*alfa)  /64.d0
    coefb(1,3)=(-5.d0 +10.d0*alfa)  /16.d0
    coefb(1,4)=(15.d0 -30.d0*alfa)  /64.d0
    coefb(1,5)=(-3.d0 + 6.d0*alfa)  /32.d0
    coefb(1,6)=( 1.d0 - 2.d0*alfa)  /64.d0
    coefb(1,7)=0.d0
    coefb(1,8)=0.d0
    !
    ! 6-order asymmetry scheme for point 0
    coefb(0,0)=( 63.d0 + 1.d0*alfa)  /64.d0
    coefb(0,1)=(  3.d0 +29.d0*alfa)  /32.d0
    coefb(0,2)=(-15.d0 +15.d0*alfa)  /64.d0
    coefb(0,3)=(  5.d0 - 5.d0*alfa)  /16.d0
    coefb(0,4)=(-15.d0 +15.d0*alfa)  /64.d0
    coefb(0,5)=(  3.d0 - 3.d0*alfa)  /32.d0
    coefb(0,6)=( -1.d0 + 1.d0*alfa)  /64.d0
    coefb(0,7)=0.d0
    coefb(0,8)=0.d0
    !
  end subroutine genfilt10coef
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! End of the Subroutine genfilt10coef.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! This subroutine is used to initial for solving the tridiagonal 
  ! martix with two layer of boundary scheme:
  ! A*x=b
  !   |1,af1,.............|
  !   |af2,1,af2,.... ....|
  !   |..af,1,af,.........|
  !   |...................|
  ! A=|...,af,1,af,.......|
  !   |...................|
  !   |.........,af,1,af..|
  !   |.........,af2,1,af2|
  !   |.............,af1,1|
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Writen by Fang Jian, 2008-11-04.
  ! Ref: ÊýÖµ·ÖÎö ±±º½³ö°æ P41
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine ptds_ini(cc,af,scheme,dim,ntype)
    !
    integer,intent(in) :: dim,ntype,scheme
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
    if(dim==0) then
      cc=0.d0
      return
    endif
    !
    if(ntype==1) then
      ! the block with boundary at i==0
      !
      cc(1,0)=1.d0
      cc(2,0)=af(1)/cc(1,0)
      !
      cc(1,1)=1.d0-af(2)*cc(2,0)
      cc(2,1)=af(2)/cc(1,1)
      !
      cc(1,2)=1.d0-af(3)*cc(2,1)
      !
      do l=2,dim-3
        cc(2,l)=af(3)/cc(1,l)
        cc(1,l+1)=1.d0-af(3)*cc(2,l)
      end do
      cc(2,dim-2)=af(3)/cc(1,dim-2)
      !
      cc(1,dim-1)=1.d0-af(3)*cc(2,dim-2)
      cc(2,dim-1)=af(3)/cc(1,dim-1)
      !
      cc(1,dim)=1.d0
      !
    elseif(ntype==2) then
      ! the block with boundary at i==im
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
      end do
      cc(2,dim-2)=af(3)/cc(1,dim-2)
      !
      cc(1,dim-1)=1.d0-af(2)*cc(2,dim-2)
      cc(2,dim-1)=af(2)/cc(1,dim-1)
      !
      cc(1,dim)=1.d0-af(1)*cc(2,dim-1)
    elseif(ntype==3) then
      ! inner block
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
    else
      print*, ' !! error in subroutine ptds_ini !'
      stop
    end if
    !
    cc(1,:)=1.d0/cc(1,:)
    !
    return
    !
  end subroutine ptds_ini
  !
  subroutine ptdsfilter_ini(cc,af,dim,ntype)
    !
    integer,intent(in) :: dim,ntype
    real(8),intent(in) :: af
    real(8),allocatable,intent(out) :: cc(:,:)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! af: input dat
    ! cc: output array
    ! dim: input dat
    ! l: temporary variable
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! local data
    integer :: l
    real(8) :: af_1,af_2,af_3
    !
    allocate(cc(1:2,0:dim))
    !
    if(dim==0) then
      cc=0.d0
      return
    endif
    !
    af_1=0.d0
    af_2=af
    af_3=af
    !
    if(ntype==1) then
      ! the block with boundary at i==0
      !
      cc(1,0)=1.d0
      cc(2,0)=af_1/cc(1,0)
      !
      cc(1,1)=1.d0-af_2*cc(2,0)
      cc(2,1)=af_2/cc(1,1)
      !
      cc(1,2)=1.d0-af_3*cc(2,1)
      !
      do l=2,dim-3
        cc(2,l)=af_3/cc(1,l)
        cc(1,l+1)=1.d0-af_3*cc(2,l)
      end do
      cc(2,dim-2)=af_3/cc(1,dim-2)
      !
      cc(1,dim-1)=1.d0-af_3*cc(2,dim-2)
      cc(2,dim-1)=af_3/cc(1,dim-1)
      !
      cc(1,dim)=1.d0
      !
    elseif(ntype==2) then
      ! the block with boundary at i==im
      !
      cc(1,0)=1.d0
      cc(2,0)=0.d0
      !
      cc(1,1)=1.d0-af_3*cc(2,0)
      cc(2,1)=af_3/cc(1,1)
      !
      cc(1,2)=1.d0-af_3*cc(2,1)
      !
      do l=2,dim-3
        cc(2,l)=af_3/cc(1,l)
        cc(1,l+1)=1.d0-af_3*cc(2,l)
      end do
      cc(2,dim-2)=af_3/cc(1,dim-2)
      !
      cc(1,dim-1)=1.d0-af_2*cc(2,dim-2)
      cc(2,dim-1)=af_2/cc(1,dim-1)
      !
      cc(1,dim)=1.d0-af_1*cc(2,dim-1)
    elseif(ntype==3) then
      ! inner block
      !
      cc(1,0)=1.d0
      cc(2,0)=0.d0
      !
      cc(1,1)=1.d0-af_3*cc(2,0)
      cc(2,1)=af_3/cc(1,1)
      !
      cc(1,2)=1.d0-af_3*cc(2,1)
      !
      do l=2,dim-3
        cc(2,l)=af_3/cc(1,l)
        cc(1,l+1)=1.d0-af_3*cc(2,l)
        ! if(mpirank==0) print*,l,cc(2,l),cc(1,l+1)
      end do
      cc(2,dim-2)=af_3/cc(1,dim-2)
      !
      cc(1,dim-1)=1.d0-af_3*cc(2,dim-2)
      cc(2,dim-1)=af_3/cc(1,dim-1)
      !
      cc(1,dim)=1.d0
    else
      print*, ' !! error in subroutine ptds_ini !'
      stop
    end if
    !
    cc(1,:)=1.d0/cc(1,:)
    !
    return
    !
  end subroutine ptdsfilter_ini
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! end of the subroutine ptds_ini.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! 
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! This subroutine is used to calculate the RHS of tridiagonal martix 
  ! for A*x=b.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function ptds_rhs(vin,dim,ns,ntype) result(vout)
    !
    integer,intent(in) :: dim,ns,ntype
    real(8),intent(in) :: vin(-hm:dim+hm)
    real(8) :: vout(0:dim)
    !
    ! local data
    integer :: l
    real(8) :: var1,var2,var3
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! ns=601: 6-o compact central scheme with
    ! boundary scheme: 2-4-6.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    if(ntype==1) then
      ! the block with boundary at i==0
      if(ns==642) then
        ! ns==642: 2-4-6-6-6-...-6-6-6-4-2
        vout(0)=2.d0*  (-vin(0)+vin(1))
        vout(1)=0.75d0*( vin(2)-vin(0))
        !
      end if
      !
      ! first order deritive
      do l=2,dim-1
        vout(l)=num7d9* (vin(l+1)-vin(l-1))+                           &
                num1d36*(vin(l+2)-vin(l-2))
      end do
      vout(dim)=0.75d0 *(vin(dim+1)-vin(dim-1))-                       &
                0.15d0 *(vin(dim+2)-vin(dim-2))+                       &
                num1d60*(vin(dim+3)-vin(dim-3))
      !
    elseif(ntype==2) then
      ! the block with boundary at i==im
      ! first order deritive
      vout(0)=0.75d0* (vin(1)-vin(-1)) -                               &
              0.15d0* (vin(2)-vin(-2)) +                               &
              num1d60*(vin(3)-vin(-3))
      do l=1,dim-2
        vout(l)=num7d9* (vin(l+1)-vin(l-1))+                           &
                num1d36*(vin(l+2)-vin(l-2))
      end do
      !
      if(ns==642) then
        ! ns==642: 2-4-6-6-6-...-6-6-6-4-2
        vout(dim-1)=0.75d0*( vin(dim)  -vin(dim-2))
        vout(dim)  =2.d0*  (-vin(dim-1)+vin(dim))
      end if
      !
    elseif(ntype==3) then
      ! inner block
      vout(0)=0.75d0* (vin(1)-vin(-1)) -                               &
              0.15d0* (vin(2)-vin(-2))+                                &
              num1d60*(vin(3)-vin(-3))
      do l=1,dim-1
        vout(l)=num7d9* (vin(l+1)-vin(l-1))+                           &
                num1d36*(vin(l+2)-vin(l-2))
      end do
        vout(dim)=0.75d0* (vin(dim+1)-vin(dim-1))-                     &
                  0.15d0* (vin(dim+2)-vin(dim-2))+                     &
                  num1d60*(vin(dim+3)-vin(dim-3))
     
    else
      print*, ' !! error in subroutine ptds_rhs !'
      stop
    end if
    !
    return
    !
  end function ptds_rhs
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! End of the function ptds_rhs
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! This subroutine is used to solve the tridiagonal martix in the 
  ! programme with two layer of boundary scheme:
  ! A*x=b
  !   |1,af1,.............|
  !   |af2,1,af2,.... ....|
  !   |..af,1,af,.........|
  !   |...................|
  ! A=|...,af,1,af,.......|
  !   |...................|
  !   |.........,af,1,af..|
  !   |.........,af2,1,af2|
  !   |.............,af1,1|
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Writen by Fang Jian, 2008-11-04.
  ! Ref: ÊýÖµ·ÖÎö ±±º½³ö°æ P41
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function ptds_cal(bd,af,cc,dim,ntype) result(xd)
    !
    integer,intent(in) :: dim,ntype
    real(8),intent(in) :: af(3),bd(0:dim),cc(1:2,0:dim)
    real(8) :: xd(0:dim)
    !
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
    !
    if(ntype==1) then
      ! the block with boundary at i==0
      !
      yd(0)=bd(0)*cc(1,0)
      yd(1)=(bd(1)-af(2)*yd(0))*cc(1,1)
      do l=2,dim-2
        yd(l)=(bd(l)-af(3)*yd(l-1))*cc(1,l)
      end do
      yd(dim-1)=(bd(dim-1)-af(3)*yd(dim-2))*cc(1,dim-1)
      yd(dim)=(bd(dim)-0.d0*yd(dim-1))*cc(1,dim)
    elseif(ntype==2) then
      ! the block with boundary at i==im
      !
      yd(0)=bd(0)*cc(1,0)
      yd(1)=(bd(1)-af(3)*yd(0))*cc(1,1)
      do l=2,dim-2
        yd(l)=(bd(l)-af(3)*yd(l-1))*cc(1,l)
      end do
      yd(dim-1)=(bd(dim-1)-af(2)*yd(dim-2))*cc(1,dim-1)
      yd(dim)=(bd(dim)-af(1)*yd(dim-1))*cc(1,dim)
    elseif(ntype==3) then
      ! inner block
      yd(0)=bd(0)*cc(1,0)
      yd(1)=(bd(1)-af(3)*yd(0))*cc(1,1)
      do l=2,dim-2
        yd(l)=(bd(l)-af(3)*yd(l-1))*cc(1,l)
      end do
      yd(dim-1)=(bd(dim-1)-af(3)*yd(dim-2))*cc(1,dim-1)
      yd(dim)=(bd(dim)-0.d0*yd(dim-1))*cc(1,dim)
    else
      print*, ' !! error in subroutine ptds_cal !'
      stop
    end if
    !
    xd(dim)=yd(dim)
    do l=dim-1,0,-1
      xd(l)=yd(l)-cc(2,l)*xd(l+1)
    end do
    !
    deallocate( yd )
    !
    return
    !
  end function ptds_cal
  !
  function ptdsfilter_cal(bd,af,cc,dim,ntype) result(xd)
    !
    integer,intent(in) :: dim,ntype
    real(8),intent(in) :: af,bd(-hm:dim+hm),cc(1:2,-hm:dim+hm)
    real(8) :: xd(0:dim)
    !
    ! local data
    integer :: i
    real(8),allocatable,dimension(:) :: yd,xyd
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! af: input dat
    ! bd: input array
    ! xd: output array
    ! cc: input array
    ! dim: input dat
    ! i, yd: temporary variable
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    allocate ( yd(-hm:dim+hm),xyd(-hm:dim+hm) )
    !
    if(ntype==1) then
      ! the block with boundary at i==0
      yd(0)=bd(0)*cc(1,0)
      do i=1,dim+hm
        yd(i)=(bd(i)-af*yd(i-1))*cc(1,i)
      end do
      !
      xyd(dim+hm)=yd(dim+hm)
      do i=dim+hm-1,0,-1
        xyd(i)=yd(i)-cc(2,i)*xyd(i+1)
      end do
      !
    elseif(ntype==2) then
      ! the block with boundary at i==im
      yd(-hm)=bd(-hm)*cc(1,-hm)
      do i=-hm+1,dim
        yd(i)=(bd(i)-af*yd(i-1))*cc(1,i)
      end do
      !
      xyd(dim)=yd(dim)
      do i=dim-1,-hm,-1
        xyd(i)=yd(i)-cc(2,i)*xyd(i+1)
      end do
      !
    elseif(ntype==3) then
      ! inner block
      yd(-hm)=bd(-hm)*cc(1,-hm)
      do i=-hm+1,dim+hm
        yd(i)=(bd(i)-af*yd(i-1))*cc(1,i)
      end do
      !
      xyd(dim+hm)=yd(dim+hm)
      do i=dim+hm-1,-hm,-1
        xyd(i)=yd(i)-cc(2,i)*xyd(i+1)
      end do
      !
    else
      print*, ' !! error in subroutine ptds_cal !'
      stop
    end if
    !
    xd(0:dim)=xyd(0:dim)
    !
    deallocate( yd,xyd )
    !
    return
    !
  end function ptdsfilter_cal
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! End of the function ptds_cal.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
end module commfunc
!+---------------------------------------------------------------------+
!| The end of the module commfunc.                                     |
!+---------------------------------------------------------------------+