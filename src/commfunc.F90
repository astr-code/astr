!+---------------------------------------------------------------------+
!| This module contains subroutine of common functions.                |
!+---------------------------------------------------------------------+
!| CHANGE RECORD                                                       |
!| -------------                                                       |
!| 08-02-2021  | Created by J. Fang                                    |
!+---------------------------------------------------------------------+
module commfunc
  !
  use commvar, only : hm,kcutoff
  use parallel, only: mpirank,mpistop,lio
  use constdef
  !
  implicit none
  !
  ! include 'fftw3.f'
  !
  interface ddfc
    module procedure ddfc_basic
    ! module procedure ddfc_2d
    ! module procedure ddfc_3d
    ! module procedure ddfc_4d
  end interface
  !
  interface spafilter10
    module procedure spafilter10_basic
    ! module procedure spafilter10_2d
  end interface
  !
  real(8),allocatable :: coef2i(:),coef4i(:),coef6i(:),coef8i(:),     &
                         coef10i(:),coefb(:,:)
  complex(8) :: ci=(0.d0,1.d0)
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
  function ddfc_basic(f,stype,ntype,dim,af,cc,lfft) result(ddfc)
    !
    use singleton
    !
    ! arguments
    character(len=4),intent(in) :: stype
    integer,intent(in) :: ntype,dim
    real(8),intent(in) :: f(-hm:dim+hm)
    real(8),intent(in),optional :: af(3),cc(1:2,0:dim)
    logical,intent(in),optional :: lfft
    real(8) :: ddfc(0:dim)
    !
    ! local data
    integer :: nscheme,k
    real(8) :: b(0:dim)
    logical :: ffton
    ! integer(8),save :: plan_f,plan_b
    !
    complex(8),allocatable :: cf(:)
    real(8),allocatable :: kama(:)
    !
    if(present(lfft)) then
      ffton=lfft
    else
      ffton=.false.
    endif
    !
    ! print*,mpirank,'|',dim,im
    !
    read(stype(1:3),*) nscheme
    !
    if(dim==0) then
      ddfc=f(1)-f(0)
    else
      !
      if(ffton) then
        !
        ! double complex in, out
        ! dimension in(N), out(N)
        ! integer*8 plan

        ! call dfftw_plan_dft_1d(plan,N,in,out,FFTW_FORWARD,FFTW_ESTIMATE)
        ! call dfftw_execute_dft(plan, in, out)
        ! call dfftw_destroy_plan(plan)
        !
        allocate(cf(1:dim),kama(1:dim))
        !
        ! call dfftw_plan_dft_1d(plan_f,dim,cf,cf,fftw_forward,fftw_estimate)
        ! call dfftw_plan_dft_1d(plan_b,dim,cf,cf,fftw_backward,fftw_estimate)
        !
        cf=dcmplx(f(1:dim))
        !
        cf=fft(cf,inv=.false.)
        ! call dfftw_execute_dft(plan_f, cf, cf)
        !
        do k=1,dim/2
          kama(k)=dble(k-1)
          !
          ! cutoff filter
          ! if(k>=kcutoff) kama(k)=0.d0
          !
        enddo
        kama(dim/2+1)=0.d0
        do k=dim/2+2,dim
          kama(k)=dble(k-dim-1)
          !
          ! cutoff filter
          ! if(dim+2-k>=kcutoff) kama(k)=0.d0
          !
        enddo
        !
        kama=kama*2.d0*pi/dble(dim)
        !
        cf=cf*kama*ci
        !
        cf=fft(cf,inv=.true.)
        ! call dfftw_execute_dft(plan_b, cf, cf)
        !
        ddfc(1:dim)=dble(cf(1:dim))
        ddfc(0)=ddfc(dim)
        !
        deallocate(cf,kama)
        !
        ! call dfftw_destroy_plan(plan_f)
        ! call dfftw_destroy_plan(plan_b)
        !
      else
        !
        if(stype(4:4)=='c') then
          b   =ptds_rhs(f,dim,nscheme,ntype)
          ddfc=ptds_cal(b,af,cc,dim,ntype)
        elseif(stype(4:4)=='e') then
          !
          if(nscheme/100==6) then
            ddfc=diff6ec(f,dim,nscheme,ntype)
          else
            print*,' !! nscheme',nscheme
            stop ' !! scheme not defined @ ddfc_basic'
          endif
        else
          !
          print*,' !! stype',stype(4:4)
          stop ' !! error at ddfc_basic' 
          !
        endif
        !
      endif
      !
    endif
    !
    return
    !
  end function ddfc_basic
  !
  function ddfc_2d(f,stype,ntype,dim,af,cc,ncolm,lfft) result(ddfc)
    !
    use singleton
    !
    ! arguments
    character(len=4),intent(in) :: stype
    integer,intent(in) :: ntype,dim,ncolm
    real(8),intent(in) :: f(-hm:dim+hm,1:ncolm)
    real(8),intent(in),optional :: af(3),cc(1:2,0:dim)
    logical,intent(in),optional :: lfft
    real(8) :: ddfc(0:dim,1:ncolm)
    !
    ! local data
    integer :: nscheme,n,k
    logical :: ffton
    real(8) :: b(0:dim,1:ncolm)
    complex(8),allocatable :: cf(:,:)
    real(8),allocatable :: kama(:)
    !
    if(present(lfft)) then
      ffton=lfft
    else
      ffton=.false.
    endif
    !
    ! print*,mpirank,'|',dim,im
    !
    read(stype(1:3),*) nscheme
    !
    if(dim==0) then
      do n=1,ncolm
        ddfc(:,n)=f(1,n)-f(0,n)
      enddo
    else
      !
      if(ffton) then
        !
        allocate(cf(1:dim,1:ncolm),kama(1:dim))
        !
        do k=1,dim/2
          kama(k)=dble(k-1)
          ! cutoff filt
          ! if(k>=kcutoff) kama(k)=0.d0
        enddo
        kama(dim/2+1)=0.d0
        do k=dim/2+2,dim
          kama(k)=dble(k-dim-1)
          ! cutoff filter
          ! if(dim+2-k>=kcutoff) kama(k)=0.d0
        enddo
        kama=kama*2.d0*pi/dble(dim)
        !
        cf=dcmplx(f(1:dim,1:ncolm))
        !
        do n=1,ncolm
          !
          cf(:,n)=fft(cf(:,n),inv=.false.)
          !
          cf(:,n)=cf(:,n)*kama*ci
          !
          cf(:,n)=fft(cf(:,n),inv=.true.)
          !
          ddfc(1:dim,n)=dble(cf(1:dim,n))
          ddfc(0,n)=ddfc(dim,n)
          !
        enddo
        !
        deallocate(cf,kama)
        !
      else
        b   =ptds2d_rhs(f,dim,nscheme,ntype,ncolm)
        ddfc=ptds2d_cal(b,af,cc,dim,ntype,ncolm)
      endif
      !
    endif
    !
    return
    !
  end function ddfc_2d
  ! !
  ! function ddfc_3d(f,stype,ntype,dim,af,cc,ncolm1,ncolm2) result(ddfc)
  !   !
  !   ! arguments
  !   character(len=4),intent(in) :: stype
  !   integer,intent(in) :: ntype,dim,ncolm1,ncolm2
  !   real(8),intent(in) :: f(-hm:dim+hm,1:ncolm1,1:ncolm2)
  !   real(8),intent(in),optional :: af(3),cc(1:2,0:dim)
  !   real(8) :: ddfc(0:dim,1:ncolm1,1:ncolm2)
  !   !
  !   ! local data
  !   integer :: nscheme,n1,n2
  !   real(8) :: b(0:dim,1:ncolm1,1:ncolm2)
  !   !
  !   ! print*,mpirank,'|',dim,im
  !   !
  !   read(stype(1:3),*) nscheme
  !   !
  !   if(dim==0) then
  !     do n2=1,ncolm2
  !     do n1=1,ncolm1
  !       ddfc(:,n1,n2)=f(1,n1,n2)-f(0,n1,n2)
  !     enddo
  !     enddo
  !   else
  !     !
  !     b   =ptds3d_rhs(f,dim,nscheme,ntype,ncolm1,ncolm2)
  !     ddfc=ptds3d_cal(b,af,cc,dim,ntype,ncolm1,ncolm2)
  !     !
  !   endif
  !   !
  !   return
  !   !
  ! end function ddfc_3d
  ! !
  ! function ddfc_4d(f,stype,ntype,dim,af,cc,ncolm1,ncolm2,ncolm3) result(ddfc)
  !   !
  !   ! arguments
  !   character(len=4),intent(in) :: stype
  !   integer,intent(in) :: ntype,dim,ncolm1,ncolm2,ncolm3
  !   real(8),intent(in) :: f(-hm:dim+hm,1:ncolm1,1:ncolm2,1:ncolm3)
  !   real(8),intent(in),optional :: af(3),cc(1:2,0:dim)
  !   real(8) :: ddfc(0:dim,1:ncolm1,1:ncolm2,1:ncolm3)
  !   !
  !   ! local data
  !   integer :: nscheme,n1,n2,n3
  !   real(8) :: b(0:dim,1:ncolm1,1:ncolm2,1:ncolm3)
  !   !
  !   ! print*,mpirank,'|',dim,im
  !   !
  !   read(stype(1:3),*) nscheme
  !   !
  !   if(dim==0) then
  !     do n3=1,ncolm3
  !     do n2=1,ncolm2
  !     do n1=1,ncolm1
  !       ddfc(:,n1,n2,n3)=f(1,n1,n2,n3)-f(0,n1,n2,n3)
  !     enddo
  !     enddo
  !     enddo
  !   else
  !     !
  !     b   =ptds4d_rhs(f,dim,nscheme,ntype,ncolm1,ncolm2,ncolm3)
  !     ddfc=ptds4d_cal(b,af,cc,dim,ntype,ncolm1,ncolm2,ncolm3)
  !     !
  !   endif
  !   !
  !   return
  !   !
  ! end function ddfc_4d
  !+-------------------------------------------------------------------+
  !| The end of the function ddf.                                      |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This function is to return the finit difference value of 6th-order|
  !| central scheme, inculding boundary closure.                       |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 19-02-2021  | Created by J. Fang @ Warrington.                    |
  !+-------------------------------------------------------------------+
  function diff6ec(vin,dim,ns,ntype) result(vout)
    !
    integer,intent(in) :: dim,ns,ntype
    real(8),intent(in) :: vin(-hm:dim+hm)
    real(8) :: vout(0:dim)
    !
    ! local data
    integer :: i
    !
    if(ntype==1) then
      !
      if(ns==642) then
        ! ns==642: 2-4-6-6-6-...-6-6-6-4-2
        vout(0)=-0.5d0*vin(2)+2.d0*vin(1)-1.5d0*vin(0)
        vout(1)=0.5d0*(vin(2)-vin(0))
        vout(2)=num2d3*(vin(3)-vin(1))-num1d12*(vin(4)-vin(0))
      else
        print*,' !! ns=',ns
        stop ' error 1 @ diff6c'
      endif
      !
      do i=3,dim
        vout(i)  =0.75d0 *(vin(i+1)-vin(i-1))-                         &
                  0.15d0 *(vin(i+2)-vin(i-2))+                         &
                  num1d60*(vin(i+3)-vin(i-3))
      enddo
      !
    elseif(ntype==2) then
      !
      do i=0,dim-3
        vout(i)  =0.75d0 *(vin(i+1)-vin(i-1))-                         &
                  0.15d0 *(vin(i+2)-vin(i-2))+                         &
                  num1d60*(vin(i+3)-vin(i-3))
      enddo
      !
      if(ns==642) then
        ! ns==642: 2-4-6-6-6-...-6-6-6-4-2
        vout(dim-2) =num2d3*(vin(dim-1)-vin(dim-3))-                   &
                    num1d12*(vin(dim)  -vin(dim-4))
        vout(dim-1)=0.5d0*(vin(dim)-vin(dim-2))
        vout(dim)  =0.5d0*vin(dim-2)-2.d0*vin(dim-1)+1.5d0*vin(dim)
      else
        print*,' !! ns=',ns
        stop ' error 2 @ diff6c'
      endif
      !
    elseif(ntype==3) then
      do i=0,dim
        vout(i)  =0.75d0 *(vin(i+1)-vin(i-1))-                         &
                  0.15d0 *(vin(i+2)-vin(i-2))+                         &
                  num1d60*(vin(i+3)-vin(i-3))
      enddo
    else
      print*,' !! ntype=',ntype
      stop ' !! errpr 3 @ diff6c'
    endif
    !
  end function diff6ec
  !+-------------------------------------------------------------------+
  !| The end of the diffcen ddf.                                       |
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
    elseif(scheme==644) then
      alfa(3)=num1d3
      alfa(2)=0.25d0
      alfa(1)=0.d0
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
  function spafilter10_basic(f,ntype,dim,af,fc,lfft) result(ff)
    !
    use singleton
    !
    ! arguments
    integer,intent(in) :: ntype,dim
    real(8),intent(in) :: f(-hm:dim+hm)
    real(8),intent(in) :: af,fc(1:2,-hm:dim+hm)
    real(8) :: ff(0:dim)
    logical,intent(in),optional :: lfft
    !
    ! local data
    real(8) :: b(0:dim)
    integer :: k
    real(8) :: aff(3)
    logical :: ffton
    complex(8),allocatable :: cf(:)
    real(8),allocatable :: kama(:)
    !
    if(present(lfft)) then
      ffton=lfft
    else
      ffton=.false.
    endif
    !
    if(ffton) then
      !
      allocate(cf(1:dim),kama(1:dim))
      !
      cf=dcmplx(f(1:dim))
      !
      cf=fft(cf,inv=.false.)
      !
      cf(dim/2+1)=0.d0 ! remove grid-grid oscillation
      ! filter the wavenumber higher than kcutoff
      do k=kcutoff,dim/2
        cf(k)=0.d0
        cf(dim+2-k)=0.d0
      enddo
      !
      ! do k=1,dim
      !   if(lio) print*,k,real(cf(k)),AIMAG(cf(k))
      ! enddo
      !
      !
      cf=fft(cf,inv=.true.)
      !
      ff(1:dim)=dble(cf(1:dim))
      ff(0)=ff(dim)
      !
      deallocate(cf,kama)
        !
    else
      aff(1)=0.d0
      aff(2)=af
      aff(3)=af
      ! b =pfilterrhs(f,dim,ntype)
      ! ff=ptdsfilter_cal(b,af,fc,dim,ntype)
      b =PFilterRHS2(f,dim,ntype)
      ff=ptds_cal(b,aff,fc,dim,ntype)
    endif
    !
    ! do i=0,dim
    !     ff(i)=filter8exp(f(i-4:i+4))
    ! enddo
    !
    return
    !
  end function spafilter10_basic
  !
  function spafilter10_2d(f,ntype,dim,af,fc,ncolm) result(ff)
    !
    ! arguments
    integer,intent(in) :: ntype,dim,ncolm
    real(8),intent(in) :: f(-hm:dim+hm,1:ncolm)
    real(8),intent(in) :: af,fc(1:2,-hm:dim+hm)
    real(8) :: ff(0:dim,1:ncolm)
    !
    ! local data
    real(8) :: b(0:dim,1:ncolm)
    integer :: i
    real(8) :: aff(3)
    !
    aff(1)=0.d0
    aff(2)=af
    aff(3)=af
    ! b =pfilterrhs(f,dim,ntype)
    ! ff=ptdsfilter_cal(b,af,fc,dim,ntype)
    b =PFilterRHS2_2d(f,dim,ntype,ncolm)
    ff=ptds2d_cal(b,aff,fc,dim,ntype,ncolm)
    !
    ! do i=0,dim
    !     ff(i)=filter8exp(f(i-4:i+4))
    ! enddo
    !
    return
    !
  end function spafilter10_2d
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
  function pfilterrhs2(var,dim,ntype) result(b)
    !
    
    integer,intent(in) :: dim,ntype 
    real(8),intent(in) :: var(-hm:dim+hm)
    real(8) :: b(0:dim)
    !
    ! local data
    integer :: l
    real(8) :: var0,var1,var2,var3,var4,var5
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! dim: dimensions in the direction.
    ! var: the variable of the filetering.
    ! var* : temporary variable.
    ! b: return variable.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    if(ntype==1) then
      ! the block with boundary at 0
      b(0)=var(0)
      !
      b(1)=coefb(1,0)*var(0)+coefb(1,1)*var(1)+coefb(1,2)*var(2)+      &
           coefb(1,3)*var(3)+coefb(1,4)*var(4)+coefb(1,5)*var(5)+      &
           coefb(1,6)*var(6)
      !
      b(2)=coefb(2,0)*var(0)+coefb(2,1)*var(1)+coefb(2,2)*var(2)+      &
           coefb(2,3)*var(3)+coefb(2,4)*var(4)+coefb(2,5)*var(5)+      &
           coefb(2,6)*var(6)
      !
      var0=var(3)+var(3)
      var1=var(4)+var(2)
      var2=var(5)+var(1)
      var3=var(6)+var(0)
      b(3)=coef6i(0)*var0+coef6i(1)*var1+coef6i(2)*var2+coef6i(3)*var3
      !
      var0=var(4)+var(4)
      var1=var(5)+var(3)
      var2=var(6)+var(2)
      var3=var(7)+var(1)
      var4=var(8)+var(0)
      b(4)=coef8i(0)*var0+coef8i(1)*var1+coef8i(2)*var2 +              &
           coef8i(3)*var3+coef8i(4)*var4
      !
      ! inner block
      do l=5,dim-1
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
      if(hm>=5) then
        b(dim)=  0.376953125d0*(var(dim)  +var(dim))     &
               + 0.205078125d0*(var(dim-1)+var(dim+1))   &
               -   0.1171875d0*(var(dim-2)+var(dim+2))   &
               +0.0439453125d0*(var(dim-3)+var(dim+3))   &
               - 0.009765625d0*(var(dim-4)+var(dim+4))   &
               +0.0009765625d0*(var(dim-5)+var(dim+5))
      elseif(hm>=4) then
        b(dim)=var(dim)-(0.2734375d0*var(dim)                          &
                     -0.21875d0*(var(dim-1)+var(dim+1))                &
                     +0.109375d0*(var(dim-2)+var(dim+2))               &
                     -3.125d-2*(var(dim-3)+var(dim+3))                 &
                   +3.90625d-3*(var(dim-4)+var(dim+4)))*1.d0
      else
        print*,' !! halo cell not sufficient !!'
        stop ' error01 @ pfilterrhs2'
      endif
      !
    elseif(ntype==2) then
      ! the block with boundary at m
      if(hm>=5) then
        b(0)=    0.376953125d0*(var(0)  +var(0))   &
               + 0.205078125d0*(var(-1) +var(1))   &
               -   0.1171875d0*(var(-2) +var(2))   &
               +0.0439453125d0*(var(-3) +var(3))   &
               - 0.009765625d0*(var(-4) +var(4))   &
               +0.0009765625d0*(var(-5) +var(5))
      elseif(hm>=4) then
        b(0)=var(0)-(0.2734375d0*var(0)-0.21875d0*(var(-1)+var(1))     &
                    +0.109375d0*(var(-2)+var(2))                       &
                    -3.125d-2*(var(-3)+var(3))                         &
                    +3.90625d-3*(var(-4)+var(4)))*1.d0
      else
        print*,' !! halo cell not sufficient !!'
        stop ' error02 @ pfilterrhs2'
      endif
      !
      ! inner block
      do l=1,dim-5
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
      var0=var(dim-4)+var(dim-4)
      var1=var(dim-3)+var(dim-5)
      var2=var(dim-2)+var(dim-6)
      var3=var(dim-1)+var(dim-7)
      var4=var(dim)  +var(dim-8)
      b(dim-4)=coef8i(0)*var0+coef8i(1)*var1+coef8i(2)*var2 +          &
               coef8i(3)*var3+coef8i(4)*var4
      !
      var0=var(dim-3)+var(dim-3)
      var1=var(dim-2)+var(dim-4)
      var2=var(dim-1)+var(dim-5)
      var3=var(dim)  +var(dim-6)
      b(dim-3)=coef6i(0)*var0+coef6i(1)*var1+                          &
               coef6i(2)*var2+coef6i(3)*var3
      !
      b(dim-2)=coefb(2,0)*var(dim)  +coefb(2,1)*var(dim-1)+            &
               coefb(2,2)*var(dim-2)+coefb(2,3)*var(dim-3)+            &
               coefb(2,4)*var(dim-4)+coefb(2,5)*var(dim-5)+            &
               coefb(2,6)*var(dim-6)
      !
      b(dim-1)=coefb(1,0)*var(dim)  +coefb(1,1)*var(dim-1)+            &
               coefb(1,2)*var(dim-2)+coefb(1,3)*var(dim-3)+            &
               coefb(1,4)*var(dim-4)+coefb(1,5)*var(dim-5)+            &
               coefb(1,6)*var(dim-6)
      !
      b(dim)=var(dim)
      !
    elseif(ntype==3) then
      ! inner block
      !
      if(hm>=5) then
        b(0)=    0.376953125d0*(var(0)  +var(0))   &
               + 0.205078125d0*(var(-1) +var(1))   &
               -   0.1171875d0*(var(-2) +var(2))   &
               +0.0439453125d0*(var(-3) +var(3))   &
               - 0.009765625d0*(var(-4) +var(4))   &
               +0.0009765625d0*(var(-5) +var(5))
      elseif(hm>=4) then
        b(0)=var(0)-(0.2734375d0*var(0)-0.21875d0*(var(-1)+var(1))     &
                    +0.109375d0*(var(-2)+var(2))                       &
                    -3.125d-2*(var(-3)+var(3))                         &
                    +3.90625d-3*(var(-4)+var(4)))*1.d0
      else
        print*,' !! halo cell not sufficient !!'
        stop ' error03 @ pfilterrhs2'
      endif
      !
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
      if(hm>=5) then
        b(dim)=  0.376953125d0*(var(dim)  +var(dim))     &
               + 0.205078125d0*(var(dim-1)+var(dim+1))   &
               -   0.1171875d0*(var(dim-2)+var(dim+2))   &
               +0.0439453125d0*(var(dim-3)+var(dim+3))   &
               - 0.009765625d0*(var(dim-4)+var(dim+4))   &
               +0.0009765625d0*(var(dim-5)+var(dim+5))
      elseif(hm>=4) then
        b(dim)=var(dim)-(0.2734375d0*var(dim)                          &
                     -0.21875d0*(var(dim-1)+var(dim+1))                &
                     +0.109375d0*(var(dim-2)+var(dim+2))               &
                     -3.125d-2*(var(dim-3)+var(dim+3))                 &
                   +3.90625d-3*(var(dim-4)+var(dim+4)))*1.d0
      else
        print*,' !! halo cell not sufficient !!'
        stop ' error04 @ pfilterrhs2'
      endif
      !
    
    else
      print*,' !! error in subroutine pfilterrhs2!'
    end if
    !
    !
  end function pfilterrhs2
  !
  function PFilterRHS2_2d(var,dim,ntype,ncolm) result(b)
    !
    
    integer,intent(in) :: dim,ntype ,ncolm
    real(8),intent(in) :: var(-hm:dim+hm,1:ncolm)
    real(8) :: b(0:dim,1:ncolm)
    !
    ! local data
    integer :: l,n
    real(8) :: var0,var1,var2,var3,var4,var5
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! dim: dimensions in the direction.
    ! var: the variable of the filetering.
    ! var* : temporary variable.
    ! b: return variable.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    if(ntype==1) then
      ! the block with boundary at 0
      do n=1,ncolm
        !
        b(0,n)=var(0,n)
        !
        b(1,n)=coefb(1,0)*var(0,n)+coefb(1,1)*var(1,n)+                  &
               coefb(1,2)*var(2,n)+coefb(1,3)*var(3,n)+                  &
               coefb(1,4)*var(4,n)+coefb(1,5)*var(5,n)+                  &
               coefb(1,6)*var(6,n)
        !
        b(2,n)=coefb(2,0)*var(0,n)+coefb(2,1)*var(1,n)+                  &
               coefb(2,2)*var(2,n)+coefb(2,3)*var(3,n)+                  &
               coefb(2,4)*var(4,n)+coefb(2,5)*var(5,n)+                  &
               coefb(2,6)*var(6,n)
        !
        b(3,n)=coef6i(0)*(var(3,n)+var(3,n))+                            &
               coef6i(1)*(var(4,n)+var(2,n))+                            &
               coef6i(2)*(var(5,n)+var(1,n))+                            &
               coef6i(3)*(var(6,n)+var(0,n))
        !
        b(4,n)=coef8i(0)*(var(4,n)+var(4,n)) +                           &
               coef8i(1)*(var(5,n)+var(3,n)) +                           &
               coef8i(2)*(var(6,n)+var(2,n)) +                           &
               coef8i(3)*(var(7,n)+var(1,n)) +                           &
               coef8i(4)*(var(8,n)+var(0,n))
        !
        ! inner block
        do l=5,dim-1
          !
          b(l,n)=coef10i(0)*(var(l,n)  +var(l,n))  +                     &
                 coef10i(1)*(var(l+1,n)+var(l-1,n))+                     &
                 coef10i(2)*(var(l+2,n)+var(l-2,n))+                     &
                 coef10i(3)*(var(l+3,n)+var(l-3,n))+                     &
                 coef10i(4)*(var(l+4,n)+var(l-4,n))+                     &
                 coef10i(5)*(var(l+5,n)+var(l-5,n)) 
          !
        end do
        !
        if(hm>=5) then
          b(dim,n)=  0.376953125d0*(var(dim,n)  +var(dim,n))             &
                   + 0.205078125d0*(var(dim-1,n)+var(dim+1,n))           &
                   -   0.1171875d0*(var(dim-2,n)+var(dim+2,n))           &
                   +0.0439453125d0*(var(dim-3,n)+var(dim+3,n))           &
                   - 0.009765625d0*(var(dim-4,n)+var(dim+4,n))           &
                   +0.0009765625d0*(var(dim-5,n)+var(dim+5,n))
        elseif(hm>=4) then
          b(dim,n)=var(dim,n)-(0.2734375d0*var(dim,n)                    &
                       -0.21875d0*(var(dim-1,n)+var(dim+1,n))            &
                       +0.109375d0*(var(dim-2,n)+var(dim+2,n))           &
                       -3.125d-2*(var(dim-3,n)+var(dim+3,n))             &
                     +3.90625d-3*(var(dim-4,n)+var(dim+4,n)))*1.d0
        else
          print*,' !! halo cell not sufficient !!'
          stop ' error01 @ PFilterRHS2_2d'
        endif
        !
      enddo
      !
    elseif(ntype==2) then
      ! the block with boundary at m
      do n=1,ncolm
        !
        if(hm>=5) then
          b(0,n)=    0.376953125d0*(var(0,n)  +var(0,n))                 &
                   + 0.205078125d0*(var(-1,n) +var(1,n))                 &
                   -   0.1171875d0*(var(-2,n) +var(2,n))                 &
                   +0.0439453125d0*(var(-3,n) +var(3,n))                 &
                   - 0.009765625d0*(var(-4,n) +var(4,n))                 &
                   +0.0009765625d0*(var(-5,n) +var(5,n))              
        elseif(hm>=4) then
          b(0,n)=var(0,n)-(0.2734375d0*var(0,n)-                         &
                        0.21875d0*(var(-1,n)+var(1,n))                   &
                      +0.109375d0*(var(-2,n)+var(2,n))                   &
                        -3.125d-2*(var(-3,n)+var(3,n))                   &
                      +3.90625d-3*(var(-4,n)+var(4,n)))*1.d0
        else
          print*,' !! halo cell not sufficient !!'
          stop ' error02 @ PFilterRHS2_2d'
        endif
        !
        ! inner block
        do l=1,dim-5
          !
          b(l,n)=coef10i(0)*(var(l,n)  +var(l,n))  +                     &
                 coef10i(1)*(var(l+1,n)+var(l-1,n))+                     &
                 coef10i(2)*(var(l+2,n)+var(l-2,n))+                     &
                 coef10i(3)*(var(l+3,n)+var(l-3,n))+                     &
                 coef10i(4)*(var(l+4,n)+var(l-4,n))+                     &
                 coef10i(5)*(var(l+5,n)+var(l-5,n)) 
          !
        end do
        !
        b(dim-4,n)=coef8i(0)*(var(dim-4,n)+var(dim-4,n)) +               &
                   coef8i(1)*(var(dim-3,n)+var(dim-5,n)) +               &
                   coef8i(2)*(var(dim-2,n)+var(dim-6,n)) +               &
                   coef8i(3)*(var(dim-1,n)+var(dim-7,n)) +               &
                   coef8i(4)*(var(dim,n)  +var(dim-8,n))
        !
        b(dim-3,n)=coef6i(0)*(var(dim-3,n)+var(dim-3,n))+                &
                   coef6i(1)*(var(dim-2,n)+var(dim-4,n))+                &
                   coef6i(2)*(var(dim-1,n)+var(dim-5,n))+                &
                   coef6i(3)*(var(dim,n)  +var(dim-6,n))
        !
        b(dim-2,n)=coefb(2,0)*var(dim,n)  +coefb(2,1)*var(dim-1,n)+      &
                   coefb(2,2)*var(dim-2,n)+coefb(2,3)*var(dim-3,n)+      &
                   coefb(2,4)*var(dim-4,n)+coefb(2,5)*var(dim-5,n)+      &
                   coefb(2,6)*var(dim-6,n)
        !  
        b(dim-1,n)=coefb(1,0)*var(dim,n)  +coefb(1,1)*var(dim-1,n)+      &
                   coefb(1,2)*var(dim-2,n)+coefb(1,3)*var(dim-3,n)+      &
                   coefb(1,4)*var(dim-4,n)+coefb(1,5)*var(dim-5,n)+      &
                   coefb(1,6)*var(dim-6,n)
        !
        b(dim,n)=var(dim,n)
        !
      enddo
      !
    elseif(ntype==3) then
      ! inner block
      !
      do n=1,ncolm
        !
        if(hm>=5) then
          b(0,n)=    0.376953125d0*(var(0,n)  +var(0,n))                 &
                   + 0.205078125d0*(var(-1,n) +var(1,n))                 &
                   -   0.1171875d0*(var(-2,n) +var(2,n))                 &
                   +0.0439453125d0*(var(-3,n) +var(3,n))                 &
                   - 0.009765625d0*(var(-4,n) +var(4,n))                 &
                   +0.0009765625d0*(var(-5,n) +var(5,n))              
        elseif(hm>=4) then
          b(0,n)=var(0,n)-(0.2734375d0*var(0,n)-                         &
                        0.21875d0*(var(-1,n)+var(1,n))                   &
                      +0.109375d0*(var(-2,n)+var(2,n))                   &
                        -3.125d-2*(var(-3,n)+var(3,n))                   &
                      +3.90625d-3*(var(-4,n)+var(4,n)))*1.d0
        else
          print*,' !! halo cell not sufficient !!'
          stop ' error03 @ PFilterRHS2_2d'
        endif
        !
        ! inner block
        do l=1,dim-1
          !
          b(l,n)=coef10i(0)*(var(l,n)  +var(l,n))  +                     &
                 coef10i(1)*(var(l+1,n)+var(l-1,n))+                     &
                 coef10i(2)*(var(l+2,n)+var(l-2,n))+                     &
                 coef10i(3)*(var(l+3,n)+var(l-3,n))+                     &
                 coef10i(4)*(var(l+4,n)+var(l-4,n))+                     &
                 coef10i(5)*(var(l+5,n)+var(l-5,n)) 
          !
        end do
        !
        if(hm>=5) then
          b(dim,n)=  0.376953125d0*(var(dim,n)  +var(dim,n))             &
                   + 0.205078125d0*(var(dim-1,n)+var(dim+1,n))           &
                   -   0.1171875d0*(var(dim-2,n)+var(dim+2,n))           &
                   +0.0439453125d0*(var(dim-3,n)+var(dim+3,n))           &
                   - 0.009765625d0*(var(dim-4,n)+var(dim+4,n))           &
                   +0.0009765625d0*(var(dim-5,n)+var(dim+5,n))
        elseif(hm>=4) then
          b(dim,n)=var(dim,n)-(0.2734375d0*var(dim,n)                    &
                       -0.21875d0*(var(dim-1,n)+var(dim+1,n))            &
                       +0.109375d0*(var(dim-2,n)+var(dim+2,n))           &
                       -3.125d-2*(var(dim-3,n)+var(dim+3,n))             &
                     +3.90625d-3*(var(dim-4,n)+var(dim+4,n)))*1.d0
        else
          print*,' !! halo cell not sufficient !!'
          stop ' error04 @ PFilterRHS2_2d'
        endif
        !
      enddo
    else
      print*,' !! error in subroutine PFilterRHS2_2d!'
    end if
    !
    !
  end function PFilterRHS2_2d
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
      elseif(ns==644) then
        ! ns==644: 4-4-6-6-6-...-6-6-6-4-4
        vout(0)=num2d3*( vin(1)-vin(-1)) - num1d12*( vin(2)-vin(-2))
        ! vout(0)=-1.5d0*vin(0)+2.d0*vin(1)-0.5d0*vin(2)
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
      elseif(ns==644) then
        ! ns==644: 4-4-6-6-6-...-6-6-6-4-4
        vout(dim-1)=0.75d0*( vin(dim)  -vin(dim-2))
        vout(dim)=num2d3*( vin(dim+1)-vin(dim-1)) -                    &
                  num1d12*( vin(dim+2)-vin(dim-2))
        
        ! vout(dim)=1.5d0*vin(dim)-2.d0*vin(dim-1)+0.5d0*vin(dim-2)
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
  !
  function ptds2d_rhs(vin,dim,ns,ntype,ncolm) result(vout)
    !
    integer,intent(in) :: dim,ns,ntype,ncolm
    real(8),intent(in) :: vin(-hm:dim+hm,1:ncolm)
    real(8) :: vout(0:dim,1:ncolm)
    !
    ! local data
    integer :: l,n
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
        do n=1,ncolm
          vout(0,n)=2.d0*  (-vin(0,n)+vin(1,n))
          vout(1,n)=0.75d0*( vin(2,n)-vin(0,n))
        enddo
        !
      elseif(ns==644) then
        do n=1,ncolm
          ! ns==644: 4-4-6-6-6-...-6-6-6-4-4
          ! vout(0,n)=-1.5d0*vin(0,n)+2.d0*vin(1,n)-0.5d0*vin(2,n)
          vout(0,n)=num2d3*( vin(1,n)-vin(-1,n)) -                     &
                   num1d12*( vin(2,n)-vin(-2,n))
          vout(1,n)=0.75d0*( vin(2,n)-vin(0,n))
        enddo
      end if
      !
      do n=1,ncolm
        ! first order deritive
        do l=2,dim-1
          vout(l,n)=num7d9* (vin(l+1,n)-vin(l-1,n))+                   &
                    num1d36*(vin(l+2,n)-vin(l-2,n))
        end do
        vout(dim,n)=0.75d0 *(vin(dim+1,n)-vin(dim-1,n))-               &
                    0.15d0 *(vin(dim+2,n)-vin(dim-2,n))+               &
                    num1d60*(vin(dim+3,n)-vin(dim-3,n))
      enddo
      !
    elseif(ntype==2) then
      ! the block with boundary at i==im
      ! first order deritive
      do n=1,ncolm
        vout(0,n)=0.75d0* (vin(1,n)-vin(-1,n)) -                      &
                  0.15d0* (vin(2,n)-vin(-2,n)) +                      &
                  num1d60*(vin(3,n)-vin(-3,n))
        do l=1,dim-2
          vout(l,n)=num7d9* (vin(l+1,n)-vin(l-1,n))+                  &
                    num1d36*(vin(l+2,n)-vin(l-2,n))
        end do
      enddo
      !
      if(ns==642) then
        ! ns==642: 2-4-6-6-6-...-6-6-6-4-2
        do n=1,ncolm
          vout(dim-1,n)=0.75d0*( vin(dim,n)  -vin(dim-2,n))
          vout(dim,n)  =2.d0*  (-vin(dim-1,n)+vin(dim,n))
        end do
      elseif(ns==644) then
        ! ns==644: 4-4-6-6-6-...-6-6-6-4-4
        do n=1,ncolm
          vout(dim-1,n)=0.75d0*( vin(dim,n)  -vin(dim-2,n))
            vout(dim,n)=num2d3*( vin(dim+1,n)-vin(dim-1,n)) -          &
                       num1d12*( vin(dim+2,n)-vin(dim-2,n))
          ! vout(dim,n)=1.5d0*vin(dim,n)-2.d0*vin(dim-1,n)+0.5d0*vin(dim-2,n)
        end do
      end if
      !
    elseif(ntype==3) then
      ! inner block
      do n=1,ncolm
        vout(0,n)=0.75d0* (vin(1,n)-vin(-1,n)) -                       &
                  0.15d0* (vin(2,n)-vin(-2,n))+                        &
                  num1d60*(vin(3,n)-vin(-3,n))
        do l=1,dim-1
          vout(l,n)=num7d9* (vin(l+1,n)-vin(l-1,n))+                   &
                    num1d36*(vin(l+2,n)-vin(l-2,n))
        end do
          vout(dim,n)=0.75d0* (vin(dim+1,n)-vin(dim-1,n))-            &
                      0.15d0* (vin(dim+2,n)-vin(dim-2,n))+            &
                      num1d60*(vin(dim+3,n)-vin(dim-3,n))
      enddo
     
    else
      print*, ' !! error in subroutine ptds_rhs !'
      stop
    end if
    !
    return
    !
  end function ptds2d_rhs
  !
  ! function ptds3d_rhs(vin,dim,ns,ntype,ncolm1,ncolm2) result(vout)
  !   !
  !   integer,intent(in) :: dim,ns,ntype,ncolm1,ncolm2
  !   real(8),intent(in) :: vin(-hm:dim+hm,1:ncolm1,1:ncolm2)
  !   real(8) :: vout(0:dim,1:ncolm1,1:ncolm2)
  !   !
  !   ! local data
  !   integer :: l,n1,n2
  !   real(8) :: var1,var2,var3
  !   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !   ! ns=601: 6-o compact central scheme with
  !   ! boundary scheme: 2-4-6.
  !   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !   !
  !   if(ntype==1) then
  !     ! the block with boundary at i==0
  !     if(ns==642) then
  !       ! ns==642: 2-4-6-6-6-...-6-6-6-4-2
  !       do n2=1,ncolm2
  !       do n1=1,ncolm1
  !         vout(0,n1,n2)=2.d0*  (-vin(0,n1,n2)+vin(1,n1,n2))
  !         vout(1,n1,n2)=0.75d0*( vin(2,n1,n2)-vin(0,n1,n2))
  !       enddo
  !       enddo
  !       !
  !     end if
  !     !
  !     do n2=1,ncolm2
  !     do n1=1,ncolm1
  !       ! first order deritive
  !       do l=2,dim-1
  !         vout(l,n1,n2)=num7d9* (vin(l+1,n1,n2)-vin(l-1,n1,n2))+      &
  !                   num1d36*(vin(l+2,n1,n2)-vin(l-2,n1,n2))
  !       end do
  !       vout(dim,n1,n2)=0.75d0 *(vin(dim+1,n1,n2)-vin(dim-1,n1,n2))-   &
  !                   0.15d0 *(vin(dim+2,n1,n2)-vin(dim-2,n1,n2))+       &
  !                   num1d60*(vin(dim+3,n1,n2)-vin(dim-3,n1,n2))
  !     enddo
  !     enddo
  !     !
  !   elseif(ntype==2) then
  !     ! the block with boundary at i==im
  !     ! first order deritive
  !     do n2=1,ncolm2
  !     do n1=1,ncolm1
  !       vout(0,n1,n2)=0.75d0* (vin(1,n1,n2)-vin(-1,n1,n2)) -           &
  !                 0.15d0* (vin(2,n1,n2)-vin(-2,n1,n2)) +               &
  !                 num1d60*(vin(3,n1,n2)-vin(-3,n1,n2))
  !       do l=1,dim-2
  !         vout(l,n1,n2)=num7d9* (vin(l+1,n1,n2)-vin(l-1,n1,n2))+       &
  !                   num1d36*(vin(l+2,n1,n2)-vin(l-2,n1,n2))
  !       end do
  !     enddo
  !     enddo
  !     !
  !     if(ns==642) then
  !       ! ns==642: 2-4-6-6-6-...-6-6-6-4-2
  !       !
  !       do n2=1,ncolm2
  !       do n1=1,ncolm1
  !         vout(dim-1,n1,n2)=0.75d0*( vin(dim,n1,n2)  -vin(dim-2,n1,n2))
  !         vout(dim,n1,n2)  =2.d0*  (-vin(dim-1,n1,n2)+vin(dim,n1,n2))
  !       end do
  !       end do
  !     end if
  !     !
  !   elseif(ntype==3) then
  !     ! inner block
  !     do n2=1,ncolm2
  !     do n1=1,ncolm1
  !       vout(0,n1,n2)=0.75d0* (vin(1,n1,n2)-vin(-1,n1,n2)) -           &
  !                 0.15d0* (vin(2,n1,n2)-vin(-2,n1,n2))+                &
  !                 num1d60*(vin(3,n1,n2)-vin(-3,n1,n2))
  !       do l=1,dim-1
  !         vout(l,n1,n2)=num7d9* (vin(l+1,n1,n2)-vin(l-1,n1,n2))+       &
  !                   num1d36*(vin(l+2,n1,n2)-vin(l-2,n1,n2))
  !       end do
  !         vout(dim,n1,n2)=0.75d0* (vin(dim+1,n1,n2)-vin(dim-1,n1,n2))- &
  !                     0.15d0* (vin(dim+2,n1,n2)-vin(dim-2,n1,n2))+     &
  !                     num1d60*(vin(dim+3,n1,n2)-vin(dim-3,n1,n2))
  !     enddo
  !     enddo
     
  !   else
  !     print*, ' !! error in subroutine ptds_rhs !'
  !     stop
  !   end if
  !   !
  !   return
  !   !
  ! end function ptds3d_rhs
  ! !
  ! function ptds4d_rhs(vin,dim,ns,ntype,ncolm1,ncolm2,ncolm3) result(vout)
  !   !
  !   integer,intent(in) :: dim,ns,ntype,ncolm1,ncolm2,ncolm3
  !   real(8),intent(in) :: vin(-hm:dim+hm,1:ncolm1,1:ncolm2,1:ncolm3)
  !   real(8) :: vout(0:dim,1:ncolm1,1:ncolm2,1:ncolm3)
  !   !
  !   ! local data
  !   integer :: l,n1,n2,n3
  !   real(8) :: var1,var2,var3
  !   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !   ! ns=601: 6-o compact central scheme with
  !   ! boundary scheme: 2-4-6.
  !   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !   !
  !   if(ntype==1) then
  !     ! the block with boundary at i==0
  !     if(ns==642) then
  !       ! ns==642: 2-4-6-6-6-...-6-6-6-4-2
  !       do n3=1,ncolm3
  !       do n2=1,ncolm2
  !       do n1=1,ncolm1
  !         vout(0,n1,n2,n3)=2.d0*  (-vin(0,n1,n2,n3)+vin(1,n1,n2,n3))
  !         vout(1,n1,n2,n3)=0.75d0*( vin(2,n1,n2,n3)-vin(0,n1,n2,n3))
  !       enddo
  !       enddo
  !       enddo
  !       !
  !     end if
  !     !
  !     do n3=1,ncolm3
  !     do n2=1,ncolm2
  !     do n1=1,ncolm1
  !       ! first order deritive
  !       do l=2,dim-1
  !         vout(l,n1,n2,n3)=num7d9* (vin(l+1,n1,n2,n3)-vin(l-1,n1,n2,n3))+      &
  !                   num1d36*(vin(l+2,n1,n2,n3)-vin(l-2,n1,n2,n3))
  !       end do
  !       vout(dim,n1,n2,n3)=0.75d0 *(vin(dim+1,n1,n2,n3)-vin(dim-1,n1,n2,n3))-   &
  !                   0.15d0 *(vin(dim+2,n1,n2,n3)-vin(dim-2,n1,n2,n3))+       &
  !                   num1d60*(vin(dim+3,n1,n2,n3)-vin(dim-3,n1,n2,n3))
  !     enddo
  !     enddo
  !     enddo
  !     !
  !   elseif(ntype==2) then
  !     ! the block with boundary at i==im
  !     ! first order deritive
  !     do n3=1,ncolm3
  !     do n2=1,ncolm2
  !     do n1=1,ncolm1
  !       vout(0,n1,n2,n3)=0.75d0* (vin(1,n1,n2,n3)-vin(-1,n1,n2,n3)) -           &
  !                 0.15d0* (vin(2,n1,n2,n3)-vin(-2,n1,n2,n3)) +               &
  !                 num1d60*(vin(3,n1,n2,n3)-vin(-3,n1,n2,n3))
  !       do l=1,dim-2
  !         vout(l,n1,n2,n3)=num7d9* (vin(l+1,n1,n2,n3)-vin(l-1,n1,n2,n3))+       &
  !                   num1d36*(vin(l+2,n1,n2,n3)-vin(l-2,n1,n2,n3))
  !       end do
  !     enddo
  !     enddo
  !     enddo
  !     !
  !     if(ns==642) then
  !       ! ns==642: 2-4-6-6-6-...-6-6-6-4-2
  !       !
  !       do n3=1,ncolm3
  !       do n2=1,ncolm2
  !       do n1=1,ncolm1
  !         vout(dim-1,n1,n2,n3)=0.75d0*( vin(dim,n1,n2,n3)  -vin(dim-2,n1,n2,n3))
  !         vout(dim,n1,n2,n3)  =2.d0*  (-vin(dim-1,n1,n2,n3)+vin(dim,n1,n2,n3))
  !       end do
  !       end do
  !       end do
  !       !
  !     end if
  !     !
  !   elseif(ntype==3) then
  !     ! inner block
  !     do n3=1,ncolm3
  !     do n2=1,ncolm2
  !     do n1=1,ncolm1
  !       vout(0,n1,n2,n3)=0.75d0* (vin(1,n1,n2,n3)-vin(-1,n1,n2,n3)) -           &
  !                 0.15d0* (vin(2,n1,n2,n3)-vin(-2,n1,n2,n3))+                &
  !                 num1d60*(vin(3,n1,n2,n3)-vin(-3,n1,n2,n3))
  !       do l=1,dim-1
  !         vout(l,n1,n2,n3)=num7d9* (vin(l+1,n1,n2,n3)-vin(l-1,n1,n2,n3))+       &
  !                   num1d36*(vin(l+2,n1,n2,n3)-vin(l-2,n1,n2,n3))
  !       end do
  !         vout(dim,n1,n2,n3)=0.75d0* (vin(dim+1,n1,n2,n3)-vin(dim-1,n1,n2,n3))- &
  !                     0.15d0* (vin(dim+2,n1,n2,n3)-vin(dim-2,n1,n2,n3))+     &
  !                     num1d60*(vin(dim+3,n1,n2,n3)-vin(dim-3,n1,n2,n3))
  !     enddo
  !     enddo
  !     enddo
     
  !   else
  !     print*, ' !! error in subroutine ptds_rhs !'
  !     stop
  !   end if
  !   !
  !   return
  !   !
  ! end function ptds4d_rhs
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
  !!
  function ptds2d_cal(bd,af,cc,dim,ntype,ncolm) result(xd)
    !
    integer,intent(in) :: dim,ntype,ncolm
    real(8),intent(in) :: af(3),bd(0:dim,1:ncolm),cc(1:2,0:dim)
    real(8) :: xd(0:dim,1:ncolm)
    !
    ! local data
    integer :: l,n
    real(8),allocatable :: yd(:,:)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! af(3): input dat
    ! bd: input array
    ! xd: output array
    ! cc: input array
    ! dim: input dat
    ! l, yd: temporary variable
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    allocate ( yd(0:dim,1:ncolm) )
    !
    if(ntype==1) then
      ! the block with boundary at i==0
      !
      do n=1,ncolm
        yd(0,n)=bd(0,n)*cc(1,0)
        yd(1,n)=(bd(1,n)-af(2)*yd(0,n))*cc(1,1)
        do l=2,dim-2
          yd(l,n)=(bd(l,n)-af(3)*yd(l-1,n))*cc(1,l)
        end do
        yd(dim-1,n)=(bd(dim-1,n)-af(3)*yd(dim-2,n))*cc(1,dim-1)
        yd(dim,n)=(bd(dim,n)-0.d0*yd(dim-1,n))*cc(1,dim)
      enddo
    elseif(ntype==2) then
      ! the block with boundary at i==im
      !
      do n=1,ncolm
        yd(0,n)=bd(0,n)*cc(1,0)
        yd(1,n)=(bd(1,n)-af(3)*yd(0,n))*cc(1,1)
        do l=2,dim-2
          yd(l,n)=(bd(l,n)-af(3)*yd(l-1,n))*cc(1,l)
        end do
        yd(dim-1,n)=(bd(dim-1,n)-af(2)*yd(dim-2,n))*cc(1,dim-1)
        yd(dim,n)=(bd(dim,n)-af(1)*yd(dim-1,n))*cc(1,dim)
      enddo
    elseif(ntype==3) then
      ! inner block
      do n=1,ncolm
        yd(0,n)=bd(0,n)*cc(1,0)
        yd(1,n)=(bd(1,n)-af(3)*yd(0,n))*cc(1,1)
        do l=2,dim-2
          yd(l,n)=(bd(l,n)-af(3)*yd(l-1,n))*cc(1,l)
        end do
        yd(dim-1,n)=(bd(dim-1,n)-af(3)*yd(dim-2,n))*cc(1,dim-1)
        yd(dim,n)=(bd(dim,n)-0.d0*yd(dim-1,n))*cc(1,dim)
      enddo
    else
      print*, ' !! error in subroutine ptds_cal !'
      stop
    end if
    !
    xd(dim,:)=yd(dim,:)
    do l=dim-1,0,-1
      xd(l,:)=yd(l,:)-cc(2,l)*xd(l+1,:)
    end do
    !
    deallocate( yd )
    !
    return
    !
  end function ptds2d_cal
  !
  ! function ptds3d_cal(bd,af,cc,dim,ntype,ncolm1,ncolm2) result(xd)
  !   !
  !   integer,intent(in) :: dim,ntype,ncolm1,ncolm2
  !   real(8),intent(in) :: af(3),bd(0:dim,1:ncolm1,1:ncolm2),cc(1:2,0:dim)
  !   real(8) :: xd(0:dim,1:ncolm1,1:ncolm2)
  !   !
  !   ! local data
  !   integer :: l,n1,n2
  !   real(8),allocatable :: yd(:,:,:)
  !   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !   ! af(3): input dat
  !   ! bd: input array
  !   ! xd: output array
  !   ! cc: input array
  !   ! dim: input dat
  !   ! l, yd: temporary variable
  !   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !   allocate ( yd(0:dim,1:ncolm1,1:ncolm2) )
  !   !
  !   if(ntype==1) then
  !     ! the block with boundary at i==0
  !     !
  !     do n2=1,ncolm2
  !     do n1=1,ncolm1
  !       yd(0,n1,n2)=bd(0,n1,n2)*cc(1,0)
  !       yd(1,n1,n2)=(bd(1,n1,n2)-af(2)*yd(0,n1,n2))*cc(1,1)
  !       do l=2,dim-2
  !         yd(l,n1,n2)=(bd(l,n1,n2)-af(3)*yd(l-1,n1,n2))*cc(1,l)
  !       end do
  !       yd(dim-1,n1,n2)=(bd(dim-1,n1,n2)-af(3)*yd(dim-2,n1,n2))*cc(1,dim-1)
  !       yd(dim,n1,n2)=(bd(dim,n1,n2)-0.d0*yd(dim-1,n1,n2))*cc(1,dim)
  !     enddo
  !     enddo
  !   elseif(ntype==2) then
  !     ! the block with boundary at i==im
  !     !
  !     do n2=1,ncolm2
  !     do n1=1,ncolm1
  !       yd(0,n1,n2)=bd(0,n1,n2)*cc(1,0)
  !       yd(1,n1,n2)=(bd(1,n1,n2)-af(3)*yd(0,n1,n2))*cc(1,1)
  !       do l=2,dim-2
  !         yd(l,n1,n2)=(bd(l,n1,n2)-af(3)*yd(l-1,n1,n2))*cc(1,l)
  !       end do
  !       yd(dim-1,n1,n2)=(bd(dim-1,n1,n2)-af(2)*yd(dim-2,n1,n2))*cc(1,dim-1)
  !       yd(dim,n1,n2)=(bd(dim,n1,n2)-af(1)*yd(dim-1,n1,n2))*cc(1,dim)
  !     enddo
  !     enddo
  !   elseif(ntype==3) then
  !     ! inner block
  !     do n2=1,ncolm2
  !     do n1=1,ncolm1
  !       yd(0,n1,n2)=bd(0,n1,n2)*cc(1,0)
  !       yd(1,n1,n2)=(bd(1,n1,n2)-af(3)*yd(0,n1,n2))*cc(1,1)
  !       do l=2,dim-2
  !         yd(l,n1,n2)=(bd(l,n1,n2)-af(3)*yd(l-1,n1,n2))*cc(1,l)
  !       end do
  !       yd(dim-1,n1,n2)=(bd(dim-1,n1,n2)-af(3)*yd(dim-2,n1,n2))*cc(1,dim-1)
  !       yd(dim,n1,n2)=(bd(dim,n1,n2)-0.d0*yd(dim-1,n1,n2))*cc(1,dim)
  !     enddo
  !     enddo
  !   else
  !     print*, ' !! error in subroutine ptds_cal !'
  !     stop
  !   end if
  !   !
  !   xd(dim,:,:)=yd(dim,:,:)
  !   do l=dim-1,0,-1
  !     xd(l,:,:)=yd(l,:,:)-cc(2,l)*xd(l+1,:,:)
  !   end do
  !   !
  !   deallocate( yd )
  !   !
  !   return
  !   !
  ! end function ptds3d_cal
  ! !
  ! function ptds4d_cal(bd,af,cc,dim,ntype,ncolm1,ncolm2,ncolm3) result(xd)
  !   !
  !   integer,intent(in) :: dim,ntype,ncolm1,ncolm2,ncolm3
  !   real(8),intent(in) :: af(3),bd(0:dim,1:ncolm1,1:ncolm2,1:ncolm3),cc(1:2,0:dim)
  !   real(8) :: xd(0:dim,1:ncolm1,1:ncolm2,1:ncolm3)
  !   !
  !   ! local data
  !   integer :: l,n1,n2,n3
  !   real(8),allocatable :: yd(:,:,:,:)
  !   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !   ! af(3): input dat
  !   ! bd: input array
  !   ! xd: output array
  !   ! cc: input array
  !   ! dim: input dat
  !   ! l, yd: temporary variable
  !   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !   allocate ( yd(0:dim,1:ncolm1,1:ncolm2,1:ncolm3) )
  !   !
  !   if(ntype==1) then
  !     ! the block with boundary at i==0
  !     !
  !     do n3=1,ncolm3
  !     do n2=1,ncolm2
  !     do n1=1,ncolm1
  !       yd(0,n1,n2,n3)=bd(0,n1,n2,n3)*cc(1,0)
  !       yd(1,n1,n2,n3)=(bd(1,n1,n2,n3)-af(2)*yd(0,n1,n2,n3))*cc(1,1)
  !       do l=2,dim-2
  !         yd(l,n1,n2,n3)=(bd(l,n1,n2,n3)-af(3)*yd(l-1,n1,n2,n3))*cc(1,l)
  !       end do
  !       yd(dim-1,n1,n2,n3)=(bd(dim-1,n1,n2,n3)-af(3)*yd(dim-2,n1,n2,n3))*cc(1,dim-1)
  !       yd(dim,n1,n2,n3)=(bd(dim,n1,n2,n3)-0.d0*yd(dim-1,n1,n2,n3))*cc(1,dim)
  !     enddo
  !     enddo
  !     enddo
  !   elseif(ntype==2) then
  !     ! the block with boundary at i==im
  !     !
  !     do n3=1,ncolm3
  !     do n2=1,ncolm2
  !     do n1=1,ncolm1
  !       yd(0,n1,n2,n3)=bd(0,n1,n2,n3)*cc(1,0)
  !       yd(1,n1,n2,n3)=(bd(1,n1,n2,n3)-af(3)*yd(0,n1,n2,n3))*cc(1,1)
  !       do l=2,dim-2
  !         yd(l,n1,n2,n3)=(bd(l,n1,n2,n3)-af(3)*yd(l-1,n1,n2,n3))*cc(1,l)
  !       end do
  !       yd(dim-1,n1,n2,n3)=(bd(dim-1,n1,n2,n3)-af(2)*yd(dim-2,n1,n2,n3))*cc(1,dim-1)
  !       yd(dim,n1,n2,n3)=(bd(dim,n1,n2,n3)-af(1)*yd(dim-1,n1,n2,n3))*cc(1,dim)
  !     enddo
  !     enddo
  !     enddo
  !   elseif(ntype==3) then
  !     ! inner block
  !     do n3=1,ncolm3
  !     do n2=1,ncolm2
  !     do n1=1,ncolm1
  !       yd(0,n1,n2,n3)=bd(0,n1,n2,n3)*cc(1,0)
  !       yd(1,n1,n2,n3)=(bd(1,n1,n2,n3)-af(3)*yd(0,n1,n2,n3))*cc(1,1)
  !       do l=2,dim-2
  !         yd(l,n1,n2,n3)=(bd(l,n1,n2,n3)-af(3)*yd(l-1,n1,n2,n3))*cc(1,l)
  !       end do
  !       yd(dim-1,n1,n2,n3)=(bd(dim-1,n1,n2,n3)-af(3)*yd(dim-2,n1,n2,n3))*cc(1,dim-1)
  !       yd(dim,n1,n2,n3)=(bd(dim,n1,n2,n3)-0.d0*yd(dim-1,n1,n2,n3))*cc(1,dim)
  !     enddo
  !     enddo
  !     enddo
  !   else
  !     print*, ' !! error in subroutine ptds_cal !'
  !     stop
  !   end if
  !   !
  !   xd(dim,:,:,:)=yd(dim,:,:,:)
  !   do l=dim-1,0,-1
  !     xd(l,:,:,:)=yd(l,:,:,:)-cc(2,l)*xd(l+1,:,:,:)
  !   end do
  !   !
  !   deallocate( yd )
  !   !
  !   return
  !   !
  ! end function ptds4d_cal
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
end module commfunc
!+---------------------------------------------------------------------+
!| The end of the module commfunc.                                     |
!+---------------------------------------------------------------------+