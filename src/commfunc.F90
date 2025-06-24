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
  !
  interface ddfc
    module procedure ddfc_basic
    module procedure ddfc_2d_lastcolum
    module procedure ddfc_3d
    module procedure ddfc_4d
  end interface
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
  !+-------------------------------------------------------------------+
  !| This function is to return the reconstructed value of a input     |
  !| function.                                                         |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 08-02-2021  | Created by J. Fang @ Warrington.                    |
  !+-------------------------------------------------------------------+
  function recons(f,stype,ntype,dim,af,cc,windir) result(fc)
    !
    ! arguments
    character(len=4),intent(in) :: stype
    integer,intent(in) :: ntype,dim
    real(8),intent(in) :: f(-hm:dim+hm)
    real(8),intent(in) :: af(1:5),cc(1:2,-2:dim+1)
    character(len=1) :: windir
    real(8) :: fc(-1:dim)
    !
    ! local data
    integer :: nscheme,k
    real(8) :: b(-2:dim+1)
    ! integer(8),save :: plan_f,plan_b
    !
    complex(8),allocatable :: cf(:)
    real(8),allocatable :: kama(:)
    !
    ! print*,mpirank,'|',dim,im
    !
    read(stype(1:3),*) nscheme
    !
    if(dim==0) then
      fc=f(0)
    else
      b =ptds_recon_rhs(f,dim,nscheme,ntype,windir)
      fc=ptds_recon_cal(b,af,cc,dim,ntype,windir)
    endif
    !
    return
    !
  end function recons
  !+-------------------------------------------------------------------+
  !| The end of the function recons.                                   |
  !+-------------------------------------------------------------------+
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
          b   =ptds_rhs(f,dim,nscheme,ntype,timerept=.true.)
          ddfc=ptds_cal(b,af,cc,dim,ntype,timerept=.true.)
        elseif(stype(4:4)=='e') then
          !
          ! if(nscheme/100==8) then
          !   ddfc=diff8ec(f,dim,nscheme,ntype)
          ! elseif(nscheme/100==6) then
          !   ddfc=diff6ec(f,dim,nscheme,ntype)
          ! elseif(nscheme/100==4) then
          !   ddfc=diff4ec(f,dim,nscheme,ntype)
          ! elseif(nscheme/100==2) then
          !   ddfc=diff2ec(f,dim,nscheme,ntype)
          ! else
          !   print*,' !! nscheme',nscheme
          !   stop ' !! scheme not defined @ ddfc_basic'
          ! endif
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
  function ddfc_2d_lastcolum(f,stype,ntype,dim,af,cc) result(ddfc)
    !
    ! arguments
    character(len=4),intent(in) :: stype
    integer,intent(in) :: ntype,dim
    real(8),intent(in) :: f(:,-hm:)
    real(8),intent(in),optional :: af(3),cc(1:2,0:dim)
    real(8) :: ddfc(ubound(f,1),0:dim)
    !
    ! local data
    integer :: nscheme,n1,n2,ncolm1
    real(8),allocatable :: b(:,:)
    !
    ncolm1=ubound(f,1)
    !
    read(stype(1:3),*) nscheme
    !
    if(dim==0) then
      do n1=1,ncolm1
        ddfc(n1,:)=f(n1,1)-f(n1,0)
      enddo
    else
      !
      allocate(b(1:ncolm1,0:dim))
      !
      b   =ptds2d_rhs_lastcolum(f,dim,nscheme,ntype,ncolm1,timerept=.true.)
      ddfc=ptds2d_cal_lastcolum(b,af,cc,dim,ntype,ncolm1,timerept=.true.)
      !
      deallocate(b)
      !
    endif
    !
    return
    !
  end function ddfc_2d_lastcolum
  !
  function ddfc_2d_firstcolum(f,stype,ntype,dim,af,cc) result(ddfc)
    !
    ! arguments
    character(len=4),intent(in) :: stype
    integer,intent(in) :: ntype,dim
    real(8),intent(in) :: f(-hm:,:)
    real(8),intent(in),optional :: af(3),cc(1:2,0:dim)
    real(8) :: ddfc(0:dim,ubound(f,2))
    !
    ! local data
    integer :: nscheme,n1,n2,ncolm1
    real(8),allocatable :: b(:,:)
    !
    ncolm1=ubound(f,2)
    !
    read(stype(1:3),*) nscheme
    !
    if(dim==0) then
      do n1=1,ncolm1
        ddfc(:,n1)=f(1,n1)-f(0,n1)
      enddo
    else
      !
      allocate(b(0:dim,1:ncolm1))
      !
      b   =ptds2d_rhs_firstcolum(f,dim,nscheme,ntype,ncolm1,timerept=.true.)
      ddfc=ptds2d_cal_firstcolum(b,af,cc,dim,ntype,ncolm1,timerept=.true.)
      !
      deallocate(b)
      !
    endif
    !
    return
    !
  end function ddfc_2d_firstcolum
  !
  function ddfc_3d(f,stype,ntype,af,cc) result(ddfc)
    !
    ! arguments
    character(len=4),intent(in) :: stype
    integer,intent(in) :: ntype
    real(8),intent(in) :: f(:,:,-hm:)
    real(8),intent(in),optional :: af(3),cc(1:2,0:ubound(f,3)-hm)
    real(8) :: ddfc(ubound(f,1),ubound(f,2),0:ubound(f,3)-hm)
    !
    ! local data
    integer :: nscheme,n1,n2,dim,ncolm1,ncolm2
    real(8),allocatable :: b(:,:,:)
    !
    ncolm1=ubound(f,1)
    ncolm2=ubound(f,2)
    dim   =ubound(f,3)-hm
    !
    read(stype(1:3),*) nscheme
    !
    if(dim==0) then
      do n2=1,ncolm2
      do n1=1,ncolm1
        ddfc(n1,n2,:)=f(n1,n2,1)-f(n1,n2,0)
      enddo
      enddo
    else
      !
      allocate(b(1:ncolm1,1:ncolm2,0:dim))
      !
      b   =ptds3d_rhs(f,dim,nscheme,ntype,ncolm1,ncolm2,timerept=.true.)
      ddfc=ptds3d_cal(b,af,cc,dim,ntype,ncolm1,ncolm2,timerept=.true.)
      !
      deallocate(b)
      !
    endif
    !
    return
    !
  end function ddfc_3d
  ! !
  function ddfc_4d(f,stype,ntype,dim,af,cc) result(ddfc)
    !
    ! arguments
    character(len=4),intent(in) :: stype
    integer,intent(in) :: ntype,dim
    real(8),intent(in) :: f(:,:,:,-hm:)
    real(8),intent(in),optional :: af(3),cc(1:2,0:dim)
    real(8) :: ddfc(ubound(f,1),ubound(f,2),ubound(f,3),0:dim)
    !
    ! local data
    integer :: nscheme,n1,n2,n3,ncolm1,ncolm2,ncolm3
    real(8),allocatable :: b(:,:,:,:)
    !
    ncolm1=ubound(f,1)
    ncolm2=ubound(f,2)
    ncolm3=ubound(f,3)
    !
    read(stype(1:3),*) nscheme
    !
    if(dim==0) then
      do n3=1,ncolm3
      do n2=1,ncolm2
      do n1=1,ncolm1
        ddfc(n1,n2,n3,:)=f(n1,n2,n3,1)-f(n1,n2,n3,0)
      enddo
      enddo
      enddo
    else
      !
      allocate(b(1:ncolm1,1:ncolm2,1:ncolm3,0:dim))
      !
      b   =ptds4d_rhs(f,dim,nscheme,ntype,ncolm1,ncolm2,ncolm3,timerept=.true.)
      ddfc=ptds4d_cal(b,af,cc,dim,ntype,ncolm1,ncolm2,ncolm3,timerept=.true.)
      !
      deallocate(b)
      !
    endif
    !
    return
    !
  end function ddfc_4d
  !+-------------------------------------------------------------------+
  !| This function is to return the finit difference value of 6th-order|
  !| central scheme, inculding boundary closure.                       |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 19-02-2021  | Created by J. Fang @ Warrington.                    |
  !+-------------------------------------------------------------------+
  function diff8ec(vin,dim,ns,ntype) result(vout)
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
        vout(i)  =(672.d0*(vin(i+1)-vin(i-1))-                         &
                   168.d0*(vin(i+2)-vin(i-2))+                         &
                    32.d0*(vin(i+3)-vin(i-3))-                         &
                     3.d0*(vin(i+4)-vin(i-4)))*num1d840
      enddo
    elseif(ntype==4) then
      do i=0,dim
        vout(i)  =(672.d0*(vin(i+1)-vin(i-1))-                         &
                   168.d0*(vin(i+2)-vin(i-2))+                         &
                    32.d0*(vin(i+3)-vin(i-3))-                         &
                     3.d0*(vin(i+4)-vin(i-4)))*num1d840
      enddo
    else
      print*,' !! ntype=',ntype
      stop ' !! errpr 3 @ diff6c'
    endif
    !
  end function diff8ec
  !
  !
  function diff4ec(vin,dim,ns,ntype) result(vout)
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
      if(ns==422) then
        ! ns==422: 2-4-4-...-4-4-2
        vout(0)=-0.5d0*vin(2)+2.d0*vin(1)-1.5d0*vin(0)
        vout(1)=0.5d0*(vin(2)-vin(0))
      else
        print*,' !! ns=',ns
        stop ' error 1 @ diff6c'
      endif
      !
      do i=2,dim
        vout(i)  =num2d3*(vin(i+1)-vin(i-1))-                         &
                  num1d12*(vin(i+2)-vin(i-2))
      enddo
      !
    elseif(ntype==2) then
      !
      do i=0,dim-2
        vout(i)  =num2d3*(vin(i+1)-vin(i-1))-                         &
                  num1d12*(vin(i+2)-vin(i-2))
      enddo
      !
      if(ns==422) then
        ! ns==422: 2-2-4-...-4-2-2
        vout(dim-1)=0.5d0*(vin(dim)-vin(dim-2))
        vout(dim)  =0.5d0*vin(dim-2)-2.d0*vin(dim-1)+1.5d0*vin(dim)
      else
        print*,' !! ns=',ns
        stop ' error 2 @ diff6c'
      endif
      !
    elseif(ntype==3) then
      do i=0,dim
        vout(i)  =num2d3*(vin(i+1)-vin(i-1))-                         &
                  num1d12*(vin(i+2)-vin(i-2))
      enddo
    else
      print*,' !! ntype=',ntype
      stop ' !! errpr 3 @ diff6c'
    endif
    !
  end function diff4ec
  !
  function diff2ec(vin,dim,ns,ntype) result(vout)
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
      ! vout(0)=-0.5d0*vin(2)+2.d0*vin(1)-1.5d0*vin(0)
      vout(0)=vin(1)-vin(0)
      !
      do i=1,dim
        vout(i)  =0.5d0*(vin(i+1)-vin(i-1))
      enddo
      !
    elseif(ntype==2) then
      !
      do i=0,dim-1
        vout(i)  =0.5d0*(vin(i+1)-vin(i-1))
      enddo
      !
      ! vout(dim)  =0.5d0*vin(dim-2)-2.d0*vin(dim-1)+1.5d0*vin(dim)
      vout(dim)  =-vin(dim-1)+vin(dim)
      !
    elseif(ntype==3) then
      do i=0,dim
        vout(i)  =0.5d0*(vin(i+1)-vin(i-1))
      enddo
    else
      print*,' !! ntype=',ntype
      stop ' !! errpr 3 @ diff6c'
    endif
    !
  end function diff2ec
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
    use commvar,  only: bfacmpld
    !
    integer,intent(in) :: scheme
    real(8),allocatable :: alfa(:)
    !
    if(scheme==642) then
      allocate(alfa(3))
      alfa(3)=num1d3
      alfa(2)=0.25d0
      alfa(1)=1.d0
    elseif(scheme==644) then
      allocate(alfa(3))
      alfa(3)=num1d3
      alfa(2)=0.25d0
      alfa(1)=0.d0
    elseif(scheme==553) then
      allocate(alfa(5))
      alfa(1)=0.5d0
      alfa(2)=0.5d0 -1.d0/6.d0*bfacmpld
      alfa(3)=num1d6+1.d0/6.d0*bfacmpld
      alfa(4)=0.5d0 -1.d0/6.d0*bfacmpld
      alfa(5)=num1d6+1.d0/6.d0*bfacmpld
    elseif(scheme==753) then
      allocate(alfa(5))
      alfa(1)=0.5d0
      alfa(2)=0.5d0 -1.d0/6.d0*bfacmpld
      alfa(3)=num1d6+1.d0/6.d0*bfacmpld
      alfa(4)=0.5d0 -1.d0/8.d0*bfacmpld
      alfa(5)=0.25d0+1.d0/8.d0*bfacmpld
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
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! This subroutine is used to initial for solving the tridiagonal 
  ! martix with two layer of boundary scheme:, considering asymmetry 
  ! upwind scheme
  ! A*x=b
  !   |1,af1,...............|
  !   |af2,1,af3,.... ......|
  !   |..af4,1,af5,.........|
  !   |.....................|
  ! A=|...,af4,1,af5,.......|
  !   |.....................|
  !   |.........,af4,1,af5..|
  !   |...........,af2,1,af3|
  !   |...............,af1,1|
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Writen by Fang Jian, 2008-11-04.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine ptds_aym_ini(cc,af,dim,ntype,windir)
    !
    integer,intent(in) :: dim,ntype
    real(8),intent(in) :: af(1:5)
    real(8),allocatable,intent(out) :: cc(:,:)
    character(len=1),intent(in) :: windir
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
    allocate(cc(1:2,-2:dim+1))
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
      if(windir=='+') then
        cc(1,1)=1.d0-af(2)*cc(2,0)
        cc(2,1)=af(3)/cc(1,1)
        do l=2,dim+1
          cc(1,l)=1.d0-af(4)*cc(2,l-1)
          cc(2,l)=af(5)/cc(1,l)
        end do
      elseif(windir=='-') then
        cc(1,1)=1.d0-af(3)*cc(2,0)
        cc(2,1)=af(2)/cc(1,1)
        do l=2,dim+1
          cc(1,l)=1.d0-af(5)*cc(2,l-1)
          cc(2,l)=af(4)/cc(1,l)
        end do
      endif
      !
      cc(1,dim+1)=1.d0
      !
    elseif(ntype==2) then
      ! the block with boundary at i==im
      !
      cc(1,-2)=1.d0
      cc(2,-2)=0.d0
      !
      if(windir=='+') then
        do l=-1,dim-3
          cc(1,l)=1.d0-af(4)*cc(2,l-1)
          cc(2,l)=af(5)/cc(1,l)
        end do
        l=dim-2
        cc(1,l)=1.d0-af(2)*cc(2,l-1)
        cc(2,l)=af(3)/cc(1,l)
        l=dim-1
        cc(1,l)=1.d0-af(1)*cc(2,l-1)
      elseif(windir=='-') then
        do l=-1,dim+1
          cc(1,l)=1.d0-af(5)*cc(2,l-1)
          cc(2,l)=af(4)/cc(1,l)
        end do
        l=dim-2
        cc(1,l)=1.d0-af(3)*cc(2,l-1)
        cc(2,l)=af(2)/cc(1,l)
        l=dim-1
        cc(1,l)=1.d0-af(1)*cc(2,l-1)
      endif
      !
    elseif(ntype==3) then
      ! inner block
      !
      cc(1,-2)=1.d0
      cc(2,-2)=0.d0
      !
      if(windir=='+') then
        do l=-1,dim+1
          cc(1,l)=1.d0-af(4)*cc(2,l-1)
          cc(2,l)=af(5)/cc(1,l)
        end do
      elseif(windir=='-') then
        do l=-1,dim+1
          cc(1,l)=1.d0-af(5)*cc(2,l-1)
          cc(2,l)=af(4)/cc(1,l)
        end do
      endif
      !
      cc(1,dim+1)=1.d0
      !
    else
      print*, ' !! error in subroutine ptds_aym_ini !'
      stop
    end if
    !
    cc(1,:)=1.d0/cc(1,:)
    !
    return
    !
  end subroutine ptds_aym_ini
  !
  function ptds_recon_rhs(vin,dim,ns,ntype,windir) result(vout)
    !
    use commvar,  only: bfacmpld
    !
    integer,intent(in) :: dim,ns,ntype
    real(8),intent(in) :: vin(-hm:dim+hm)
    character(len=1),intent(in) :: windir
    real(8) :: vout(-2:dim+1)
    !
    ! local data
    integer :: l
    !
    if(ntype==1) then
      ! the block with boundary at i==0
      !
      if(ns==553) then
        ! ns==553 3-5-5-...-5-5-3
        l=0
        vout(0)=0.25d0*vin(l)+1.25d0*vin(l+1)
        if(windir=='+') then
          do l=1,dim
            vout(l)=(num1d18 -num1d36*bfacmpld)*vin(l-1) + &
                    (num19d18-num9d36*bfacmpld)*vin(l)   + &
                    (num5d9  +num9d36*bfacmpld)*vin(l+1) + &
                              num1d36*bfacmpld *vin(l+2)
          end do
          l=dim+1
          vout(l)=(2.d0*vin(l-2)-13.d0*vin(l-1)+47.d0*vin(l)+          &
                   27.d0*vin(l+1)-3.d0*vin(l+2))/60.d0
        elseif(windir=='-') then
          do l=1,dim
            vout(l)=(num1d18 -num1d36*bfacmpld)*vin(l+2) + &
                    (num19d18-num9d36*bfacmpld)*vin(l+1) + &
                    (num5d9  +num9d36*bfacmpld)*vin(l)   + &
                              num1d36*bfacmpld *vin(l-1)
          end do
          l=dim+1
          vout(l)=(2.d0*vin(l+3)-13.d0*vin(l+2)+47.d0*vin(l+1)+        &
                   27.d0*vin(l)  -3.d0*vin(l-1))/60.d0
        endif
        !
      elseif(ns==753) then
        ! ns==553 3-5-7-...-7-5-3
        !
        ! 3rd-order compact upwind
        l=0
        vout(0)=0.25d0*vin(l)+1.25d0*vin(l+1)
        !
        if(windir=='+') then
          !
          ! 5th-order compact upwind near boundary
          l=1
          vout(l)=(num1d18 -num1d36*bfacmpld)*vin(l-1) + &
                  (num19d18-num9d36*bfacmpld)*vin(l)   + &
                  (num5d9  +num9d36*bfacmpld)*vin(l+1) + &
                            num1d36*bfacmpld *vin(l+2)
          !
          ! 7th-order compact upwind
          do l=2,dim
            vout(l)=(( -1.d0 + 0.5d0*bfacmpld)*vin(l-2) +   &
                     ( 19.d0 - 7.5d0*bfacmpld)*vin(l-1) +   &
                     (239.d0 - 40.d0*bfacmpld)*vin(l)   +   &
                     (159.d0 + 40.d0*bfacmpld)*vin(l+1) +   &
                     (  4.d0 + 7.5d0*bfacmpld)*vin(l+2)     &
                             - 0.5d0*bfacmpld *vin(l+3))/240.d0
          end do
          !
          ! 7th-order explicit upwind
          l=dim+1
          vout(l)=(   -3.d0*vin(l-3)+  &
                      25.d0*vin(l-2)-  &
                     101.d0*vin(l-1)+  &
                     319.d0*vin(l)+    &
                     214.d0*vin(l+1)-  &
                      38.d0*vin(l+2)+  &
                       4.d0*vin(l+3) )/420.d0
          !
        elseif(windir=='-') then
          !
          ! 5th-order compact upwind near boundary
          l=1
          vout(l)=(num1d18 -num1d36*bfacmpld)*vin(l+2) + &
                  (num19d18-num9d36*bfacmpld)*vin(l+1) + &
                  (num5d9  +num9d36*bfacmpld)*vin(l)   + &
                            num1d36*bfacmpld *vin(l-1)
          !
          ! 7th-order compact upwind
          do l=2,dim
            vout(l)=(( -1.d0 + 0.5d0*bfacmpld)*vin(l+3) +   &
                     ( 19.d0 - 7.5d0*bfacmpld)*vin(l+2) +   &
                     (239.d0 - 40.d0*bfacmpld)*vin(l+1) +   &
                     (159.d0 + 40.d0*bfacmpld)*vin(l)   +   &
                     (  4.d0 + 7.5d0*bfacmpld)*vin(l-1)     &
                             - 0.5d0*bfacmpld *vin(l-2))/240.d0
          end do
          !
          ! 7th-order explicit upwind for interface
          l=dim+1
          vout(l)=(   -3.d0*vin(l+4)+  &
                      25.d0*vin(l+3)-  &
                     101.d0*vin(l+2)+  &
                     319.d0*vin(l+1)+  &
                     214.d0*vin(l)-    &
                      38.d0*vin(l-1)+  &
                       4.d0*vin(l-2) )/420.d0
        endif
        !
      end if
      !
    elseif(ntype==2) then
      ! the block with boundary at i==im
      !
      if(ns==553) then
        ! ns==553 3-5-5-...-5-5-3
        !
        if(windir=='+') then
          l=-2
          vout(l)=(2.d0*vin(l-2)-13.d0*vin(l-1)+47.d0*vin(l)+          &
                   27.d0*vin(l+1)-3.d0*vin(l+2))/60.d0
          do l=-1,dim-2
            vout(l)=(num1d18 -num1d36*bfacmpld)*vin(l-1) + &
                    (num19d18-num9d36*bfacmpld)*vin(l)   + &
                    (num5d9  +num9d36*bfacmpld)*vin(l+1) + &
                              num1d36*bfacmpld *vin(l+2)
          end do
        elseif(windir=='-') then
          l=-2
          vout(l)=(2.d0*vin(l+3)-13.d0*vin(l+2)+47.d0*vin(l+1)+        &
                   27.d0*vin(l)  -3.d0*vin(l-1))/60.d0
          do l=-1,dim-2
            vout(l)=(num1d18 -num1d36*bfacmpld)*vin(l+2) + &
                    (num19d18-num9d36*bfacmpld)*vin(l+1) + &
                    (num5d9  +num9d36*bfacmpld)*vin(l)   + &
                              num1d36*bfacmpld *vin(l-1)
          end do
        endif
        l=dim-1
        vout(l)=0.25d0*vin(l+1)+1.25d0*vin(l)
        !
      elseif(ns==753) then
        ! ns==753 3-5-7-...-7-5-3
        !
        if(windir=='+') then
          !
          l=-2
          vout(l)=(   -3.d0*vin(l-3)+  &
                      25.d0*vin(l-2)-  &
                     101.d0*vin(l-1)+  &
                     319.d0*vin(l)+    &
                     214.d0*vin(l+1)-  &
                      38.d0*vin(l+2)+  &
                       4.d0*vin(l+3) )/420.d0
          ! 7th-order compact upwind
          do l=-1,dim-3
            vout(l)=(( -1.d0 + 0.5d0*bfacmpld)*vin(l-2) +   &
                     ( 19.d0 - 7.5d0*bfacmpld)*vin(l-1) +   &
                     (239.d0 - 40.d0*bfacmpld)*vin(l)   +   &
                     (159.d0 + 40.d0*bfacmpld)*vin(l+1) +   &
                     (  4.d0 + 7.5d0*bfacmpld)*vin(l+2)     &
                             - 0.5d0*bfacmpld *vin(l+3))/240.d0
          end do
          !
          ! 5th-order compact upwind near the boundary
          l=dim-2
          vout(l)=(num1d18 -num1d36*bfacmpld)*vin(l-1) + &
                  (num19d18-num9d36*bfacmpld)*vin(l)   + &
                  (num5d9  +num9d36*bfacmpld)*vin(l+1) + &
                            num1d36*bfacmpld *vin(l+2)
          !
        elseif(windir=='-') then
          !
          ! 7th-order explicit upwind for interface
          l=-2
          vout(l)=(   -3.d0*vin(l+4)+  &
                      25.d0*vin(l+3)-  &
                     101.d0*vin(l+2)+  &
                     319.d0*vin(l+1)+  &
                     214.d0*vin(l)-    &
                      38.d0*vin(l-1)+  &
                       4.d0*vin(l-2) )/420.d0
          ! 7th-order compact upwind
          do l=-1,dim-3
            vout(l)=(( -1.d0 + 0.5d0*bfacmpld)*vin(l+3) +   &
                     ( 19.d0 - 7.5d0*bfacmpld)*vin(l+2) +   &
                     (239.d0 - 40.d0*bfacmpld)*vin(l+1) +   &
                     (159.d0 + 40.d0*bfacmpld)*vin(l)   +   &
                     (  4.d0 + 7.5d0*bfacmpld)*vin(l-1)     &
                             - 0.5d0*bfacmpld *vin(l-2))/240.d0
          end do
          !
          ! 5th-order compact upwind near the boundary
          l=dim-2
          vout(l)=(num1d18 -num1d36*bfacmpld)*vin(l+2) + &
                  (num19d18-num9d36*bfacmpld)*vin(l+1) + &
                  (num5d9  +num9d36*bfacmpld)*vin(l)   + &
                            num1d36*bfacmpld *vin(l-1)
          !
        endif
        !
        ! 3rd-order compact upwind at the boundary
        l=dim-1
        vout(l)=0.25d0*vin(l+1)+1.25d0*vin(l)
        !
      else
        print*,' ** ntype: ',ntype
        print*, ' !! error 2 in subroutine ptds_recon_rhs !'
        stop
      endif
      !
    elseif(ntype==3) then
      !
      ! inner block
      if(ns/100==5) then
        !
        if(windir=='+') then
          !
          l=-2
          vout(l)=(2.d0*vin(l-2)-13.d0*vin(l-1)+47.d0*vin(l)+          &
                   27.d0*vin(l+1)-3.d0*vin(l+2))/60.d0
          do l=-1,dim
            vout(l)=(num1d18 -num1d36*bfacmpld)*vin(l-1) + &
                    (num19d18-num9d36*bfacmpld)*vin(l)   + &
                    (num5d9  +num9d36*bfacmpld)*vin(l+1) + &
                              num1d36*bfacmpld *vin(l+2)
          end do
          l=dim+1
          vout(l)=(2.d0*vin(l-2)-13.d0*vin(l-1)+47.d0*vin(l)+          &
                   27.d0*vin(l+1)-3.d0*vin(l+2))/60.d0
        elseif(windir=='-') then
          l=-2
          vout(l)=(2.d0*vin(l+3)-13.d0*vin(l+2)+47.d0*vin(l+1)+        &
                   27.d0*vin(l)  -3.d0*vin(l-1))/60.d0
          do l=-1,dim
            vout(l)=(num1d18 -num1d36*bfacmpld)*vin(l+2) + &
                    (num19d18-num9d36*bfacmpld)*vin(l+1) + &
                    (num5d9  +num9d36*bfacmpld)*vin(l)   + &
                              num1d36*bfacmpld *vin(l-1)
          end do
          l=dim+1
          vout(l)=(2.d0*vin(l+3)-13.d0*vin(l+2)+47.d0*vin(l+1)+        &
                   27.d0*vin(l)  -3.d0*vin(l-1))/60.d0
        endif
        !
      elseif(ns/100==7) then
        !
        if(windir=='+') then
          !
          ! 7th-order explicit upwind for interface
          l=-2
          vout(l)=(   -3.d0*vin(l-3)+  &
                      25.d0*vin(l-2)-  &
                     101.d0*vin(l-1)+  &
                     319.d0*vin(l)+    &
                     214.d0*vin(l+1)-  &
                      38.d0*vin(l+2)+  &
                       4.d0*vin(l+3) )/420.d0
          ! 
          ! 7th-order compact scheme
          do l=-1,dim
            vout(l)=(( -1.d0 + 0.5d0*bfacmpld)*vin(l-2) +   &
                     ( 19.d0 - 7.5d0*bfacmpld)*vin(l-1) +   &
                     (239.d0 - 40.d0*bfacmpld)*vin(l)   +   &
                     (159.d0 + 40.d0*bfacmpld)*vin(l+1) +   &
                     (  4.d0 + 7.5d0*bfacmpld)*vin(l+2)     &
                             - 0.5d0*bfacmpld *vin(l+3))/240.d0
          end do
          !
          ! 7th-order explicit upwind for interface
          l=dim+1
          vout(l)=(   -3.d0*vin(l-3)+  &
                      25.d0*vin(l-2)-  &
                     101.d0*vin(l-1)+  &
                     319.d0*vin(l)+    &
                     214.d0*vin(l+1)-  &
                      38.d0*vin(l+2)+  &
                       4.d0*vin(l+3) )/420.d0
          !
        elseif(windir=='-') then
          !
          ! 7th-order explicit upwind for interface
          l=-2
          vout(l)=(   -3.d0*vin(l+4)+  &
                      25.d0*vin(l+3)-  &
                     101.d0*vin(l+2)+  &
                     319.d0*vin(l+1)+  &
                     214.d0*vin(l)-    &
                      38.d0*vin(l-1)+  &
                       4.d0*vin(l-2) )/420.d0
          ! 7th-order compact scheme
          do l=-1,dim
            vout(l)=(( -1.d0 + 0.5d0*bfacmpld)*vin(l+3) +   &
                     ( 19.d0 - 7.5d0*bfacmpld)*vin(l+2) +   &
                     (239.d0 - 40.d0*bfacmpld)*vin(l+1) +   &
                     (159.d0 + 40.d0*bfacmpld)*vin(l)   +   &
                     (  4.d0 + 7.5d0*bfacmpld)*vin(l-1)     &
                             - 0.5d0*bfacmpld *vin(l-2))/240.d0
          end do
          ! 7th-order explicit upwind for interface
          l=dim+1
          vout(l)=(   -3.d0*vin(l+4)+  &
                      25.d0*vin(l+3)-  &
                     101.d0*vin(l+2)+  &
                     319.d0*vin(l+1)+  &
                     214.d0*vin(l)-    &
                      38.d0*vin(l-1)+  &
                       4.d0*vin(l-2) )/420.d0
        endif
        !
      endif
      !
    else
      print*,' ** ntype: ',ntype
      print*, ' !! error in subroutine ptds_recon_rhs !'
      stop
    end if
  end function ptds_recon_rhs
  !
  function ptds_recon_cal(bd,af,cc,dim,ntype,windir) result(xd)
    !
    integer,intent(in) :: dim,ntype
    real(8),intent(in) :: af(5),bd(-2:dim+1),cc(1:2,-2:dim+1)
    character(len=1),intent(in) :: windir
    real(8) :: xd(-1:dim)
    !
    ! local data
    integer :: l
    real(8),allocatable,dimension(:) :: yd,md
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! af(3): input dat
    ! bd: input array
    ! xd: output array
    ! cc: input array
    ! dim: input dat
    ! l, yd: temporary variable
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    allocate ( yd(-2:dim+1),md(-2:dim+1) )
    !
    if(ntype==1) then
      ! the block with boundary at i==0
      !
      yd(0)=bd(0)*cc(1,0)
      if(windir=='+') then
        yd(1)=(bd(1)-af(2)*yd(0))*cc(1,1)
        do l=2,dim
          yd(l)=(bd(l)-af(4)*yd(l-1))*cc(1,l)
        end do
      elseif(windir=='-') then
        yd(1)=(bd(1)-af(3)*yd(0))*cc(1,1)
        do l=2,dim
          yd(l)=(bd(l)-af(5)*yd(l-1))*cc(1,l)
        end do
      endif
      yd(dim+1)=bd(dim+1)*cc(1,dim+1)
      !
      md(dim+1)=yd(dim+1)
      do l=dim,-2,-1
        md(l)=yd(l)-cc(2,l)*md(l+1)
      end do
      !
      xd(0:dim)=md(0:dim)
      !
    elseif(ntype==2) then
      ! the block with boundary at i==im
      !
      yd(-2)=bd(-2)*cc(1,-2)
      if(windir=='+') then
        do l=-1,dim-3
          yd(l)=(bd(l)-af(4)*yd(l-1))*cc(1,l)
        end do
        l=dim-2
        yd(l)=(bd(l)-af(2)*yd(l-1))*cc(1,l)
      elseif(windir=='-') then
        do l=-1,dim-3
          yd(l)=(bd(l)-af(5)*yd(l-1))*cc(1,l)
        end do
        l=dim-2
        yd(l)=(bd(l)-af(3)*yd(l-1))*cc(1,l)
      endif
      l=dim-1
      yd(l)=(bd(l)-af(1)*yd(l-1))*cc(1,l)
      !
      md(dim-1)=yd(dim-1)
      do l=dim-2,-2,-1
        md(l)=yd(l)-cc(2,l)*md(l+1)
      end do
      !
      xd(-1:dim-1)=md(-1:dim-1)
      !
    elseif(ntype==3) then
      !
      ! inner block
      !
      yd(-2)=bd(-2)*cc(1,-2)
      if(windir=='+') then
        do l=-1,dim
          yd(l)=(bd(l)-af(4)*yd(l-1))*cc(1,l)
        end do
      elseif(windir=='-') then
        do l=-1,dim
          yd(l)=(bd(l)-af(5)*yd(l-1))*cc(1,l)
        end do
      endif
      yd(dim+1)=bd(dim+1)*cc(1,dim+1)
      !
      md(dim+1)=yd(dim+1)
      do l=dim,-2,-1
        md(l)=yd(l)-cc(2,l)*md(l+1)
      end do
      !
      xd(-1:dim)=md(-1:dim)
      !
    else
      print*,' ** ntype: ',ntype
      print*, ' !! error in subroutine ptds_recon_cal !'
      stop
    end if
    !
    deallocate( yd )
    !
    return
    !
  end function ptds_recon_cal
  !+-------------------------------------------------------------------+
  !| The end of the function ptds_recon_.                              |
  !+-------------------------------------------------------------------+
  !
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
  subroutine ptds_ini(cc,af,dim,ntype)
    !
    integer,intent(in) :: dim,ntype
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
    elseif(ntype==4) then
      ! the block with boundary at i==0 .and. i==im

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
      cc(1,dim-1)=1.d0-af(2)*cc(2,dim-2)
      cc(2,dim-1)=af(2)/cc(1,dim-1)
      !
      cc(1,dim)=1.d0-af(1)*cc(2,dim-1)
      !
    else
      print*,'ntype=',ntype
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
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! This subroutine is used to calculate the RHS of tridiagonal martix 
  ! for A*x=b.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function ptds_rhs(vin,dim,ns,ntype,timerept) result(vout)
    !
    integer,intent(in) :: dim,ns,ntype
    real(8),intent(in) :: vin(-hm:dim+hm)
    real(8) :: vout(0:dim)
    logical,intent(in),optional :: timerept
    !
    ! local data
    integer :: l
    real(8) :: var1,var2,var3
    !
    real(8) :: time_beg
    real(8),save :: subtime=0.d0
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! ns=601: 6-o compact central scheme with
    ! boundary scheme: 2-4-6.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    if(present(timerept) .and. timerept) time_beg=ptime() 
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
        vout(0)=0.5d0*( vin(1)-vin(-1))
        ! vout(0)=num2d3*( vin(1)-vin(-1)) - num1d12*( vin(2)-vin(-2))
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
        vout(dim)=0.5d0 *( vin(dim+1)-vin(dim-1))
        ! vout(dim)=num2d3 *( vin(dim+1)-vin(dim-1)) -                   &
        !           num1d12*( vin(dim+2)-vin(dim-2))
        
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
     
    elseif(ntype==4) then
      ! the block with boundary at i=0 and i=im
      !
      if(ns==642) then
        ! ns==642: 2-4-6-6-6-...-6-6-6-4-2
        vout(0)=2.d0*  (-vin(0)+vin(1))
        vout(1)=0.75d0*( vin(2)-vin(0))
        !
      elseif(ns==644) then
        ! ns==644: 4-4-6-6-6-...-6-6-6-4-4
        vout(0)=0.5d0*( vin(1)-vin(-1))
        ! vout(0)=num2d3*( vin(1)-vin(-1)) - num1d12*( vin(2)-vin(-2))
        ! vout(0)=-1.5d0*vin(0)+2.d0*vin(1)-0.5d0*vin(2)
        vout(1)=0.75d0*( vin(2)-vin(0))
        !
      end if
      !
      ! first order deritive
      do l=2,dim-2
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
        vout(dim)=0.5d0 *( vin(dim+1)-vin(dim-1))
        ! vout(dim)=num2d3 *( vin(dim+1)-vin(dim-1)) -                   &
        !           num1d12*( vin(dim+2)-vin(dim-2))
        
        ! vout(dim)=1.5d0*vin(dim)-2.d0*vin(dim-1)+0.5d0*vin(dim-2)
      end if
      !
    else
      print*,ntype
      print*, ' !! ntype error in subroutine ptds_rhs !'
      stop
    end if
    !
    if(present(timerept) .and. timerept) then
      !
      subtime=subtime+ptime()-time_beg
      !
      if(lio .and. lreport .and. ltimrpt) call timereporter(routine='ptds_rhs', &
                                             timecost=subtime, &
                                              message='rhs of compact scheme')
    endif
    !
    return
    !
  end function ptds_rhs
  !
  function ptds2d_rhs_lastcolum(vin,dim,ns,ntype,ncolm1,timerept) result(vout)
    !
    integer,intent(in) :: dim,ns,ntype,ncolm1
    real(8),intent(in) :: vin(1:ncolm1,-hm:dim+hm)
    real(8) :: vout(1:ncolm1,0:dim)
    logical,intent(in),optional :: timerept
    !
    ! local data
    integer :: l
    real(8) :: var1,var2,var3
    !
    real(8) :: time_beg
    real(8),save :: subtime=0.d0
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! ns=601: 6-o compact central scheme with
    ! boundary scheme: 2-4-6.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    if(present(timerept) .and. timerept) time_beg=ptime() 
    !
    if(ntype==1) then
      ! the block with boundary at i==0
      if(ns==642) then
        ! ns==642: 2-4-6-6-6-...-6-6-6-4-2
        vout(:,0)=2.d0*  (-vin(:,0)+vin(:,1))
        vout(:,1)=0.75d0*( vin(:,2)-vin(:,0))
        !
      elseif(ns==644) then
        ! ns==644: 4-4-6-6-6-...-6-6-6-4-4
        vout(:,0)= 0.5d0*( vin(:,1)-vin(:,-1))
        vout(:,1)=0.75d0*( vin(:,2)-vin(:,0))
        !
      end if
      !
      ! first order deritive
      do l=2,dim-1
        vout(:,l)=num7d9* (vin(:,l+1)-vin(:,l-1))+               &
                 num1d36*(vin(:,l+2)-vin(:,l-2))
      end do
      vout(:,dim)=0.75d0 *(vin(:,dim+1)-vin(:,dim-1))-           &
                  0.15d0 *(vin(:,dim+2)-vin(:,dim-2))+           &
                  num1d60*(vin(:,dim+3)-vin(:,dim-3))
      !
    elseif(ntype==2) then
      ! the block with boundary at i==im
      ! first order deritive
      vout(:,0)=0.75d0* (vin(:,1)-vin(:,-1)) -                  &
                0.15d0* (vin(:,2)-vin(:,-2)) +                  &
                num1d60*(vin(:,3)-vin(:,-3))
      do l=1,dim-2
        vout(:,l)=num7d9* (vin(:,l+1)-vin(:,l-1))+              &
                  num1d36*(vin(:,l+2)-vin(:,l-2))
      end do
      !
      if(ns==642) then
        ! ns==642: 2-4-6-6-6-...-6-6-6-4-2
        vout(:,dim-1)=0.75d0*( vin(:,dim)  -vin(:,dim-2))
        vout(:,dim)  =2.d0*  (-vin(:,dim-1)+vin(:,dim))
      elseif(ns==644) then
        ! ns==644: 4-4-6-6-6-...-6-6-6-4-4
        vout(:,dim-1)=0.75d0*( vin(:,dim)  -vin(:,dim-2))
        vout(:,dim)=  0.5d0 *( vin(:,dim+1)-vin(:,dim-1))
      end if
      !
    elseif(ntype==3) then
      ! inner block
      vout(:,0)=0.75d0* (vin(:,1)-vin(:,-1)) -                 &
                0.15d0* (vin(:,2)-vin(:,-2))+                  &
                num1d60*(vin(:,3)-vin(:,-3))
      do l=1,dim-1
        vout(:,l)=num7d9* (vin(:,l+1)-vin(:,l-1))+             &
                  num1d36*(vin(:,l+2)-vin(:,l-2))
      end do
        vout(:,dim)=0.75d0* (vin(:,dim+1)-vin(:,dim-1))-       &
                    0.15d0* (vin(:,dim+2)-vin(:,dim-2))+       &
                    num1d60*(vin(:,dim+3)-vin(:,dim-3))
     
    elseif(ntype==4) then
      ! the block with boundary at i=0 and i=im
      !
      if(ns==642) then
        ! ns==642: 2-4-6-6-6-...-6-6-6-4-2
        vout(:,0)=2.d0*  (-vin(:,0)+vin(:,1))
        vout(:,1)=0.75d0*( vin(:,2)-vin(:,0))
        !
      elseif(ns==644) then
        ! ns==644: 4-4-6-6-6-...-6-6-6-4-4
        vout(:,0)=0.5d0*( vin(:,1)-vin(:,-1))
        vout(:,1)=0.75d0*( vin(:,2)-vin(:,0))
        !
      end if
      !
      ! first order deritive
      do l=2,dim-2
        vout(:,l)=num7d9* (vin(:,l+1)-vin(:,l-1))+               &
                  num1d36*(vin(:,l+2)-vin(:,l-2))
      end do
      !
      if(ns==642) then
        ! ns==642: 2-4-6-6-6-...-6-6-6-4-2
        vout(:,dim-1)=0.75d0*( vin(:,dim)  -vin(:,dim-2))
        vout(:,dim)  =2.d0*  (-vin(:,dim-1)+vin(:,dim))
      elseif(ns==644) then
        ! ns==644: 4-4-6-6-6-...-6-6-6-4-4
        vout(:,dim-1)=0.75d0*( vin(:,dim)  -vin(:,dim-2))
        vout(:,dim)  =0.5d0 *( vin(:,dim+1)-vin(:,dim-1))
      end if
      !
    else
      print*,ntype
      print*, ' !! ntype error in subroutine ptds_rhs !'
      stop
    end if
    !
    if(present(timerept) .and. timerept) then
      !
      subtime=subtime+ptime()-time_beg
      !
      if(lio .and. lreport .and. ltimrpt) call timereporter(routine='ptds2d_rhs', &
                                             timecost=subtime,      &
                                              message='rhs of compact scheme')
    endif
    !
    return
    !
  end function ptds2d_rhs_lastcolum
  !
  function ptds2d_rhs_firstcolum(vin,dim,ns,ntype,ncolm1,timerept) result(vout)
    !
    integer,intent(in) :: dim,ns,ntype,ncolm1
    real(8),intent(in) :: vin(-hm:dim+hm,1:ncolm1)
    real(8) :: vout(0:dim,1:ncolm1)
    logical,intent(in),optional :: timerept
    !
    ! local data
    integer :: l
    real(8) :: var1,var2,var3
    !
    real(8) :: time_beg
    real(8),save :: subtime=0.d0
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! ns=601: 6-o compact central scheme with
    ! boundary scheme: 2-4-6.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    if(present(timerept) .and. timerept) time_beg=ptime() 
    !
    if(ntype==1) then
      ! the block with boundary at i==0
      if(ns==642) then
        ! ns==642: 2-4-6-6-6-...-6-6-6-4-2
        vout(0,:)=2.d0*  (-vin(0,:)+vin(1,:))
        vout(1,:)=0.75d0*( vin(2,:)-vin(0,:))
        !
      elseif(ns==644) then
        ! ns==644: 4-4-6-6-6-...-6-6-6-4-4
        vout(0,:)= 0.5d0*( vin(1,:)-vin(-1,:))
        vout(1,:)=0.75d0*( vin(2,:)-vin(0,:))
        !
      end if
      !
      ! first order deritive
      do l=2,dim-1
        vout(l,:)=num7d9* (vin(l+1,:)-vin(l-1,:))+               &
                 num1d36* (vin(l+2,:)-vin(l-2,:))
      end do
      vout(dim,:)=0.75d0 *(vin(dim+1,:)-vin(dim-1,:))-           &
                  0.15d0 *(vin(dim+2,:)-vin(dim-2,:))+           &
                  num1d60*(vin(dim+3,:)-vin(dim-3,:))
      !
    elseif(ntype==2) then
      ! the block with boundary at i==im
      ! first order deritive
      vout(0,:)=0.75d0* (vin(1,:)-vin(-1,:)) -                  &
                0.15d0* (vin(2,:)-vin(-2,:)) +                  &
                num1d60*(vin(3,:)-vin(-3,:))
      do l=1,dim-2
        vout(l,:)=num7d9* (vin(l+1,:)-vin(l-1,:))+              &
                  num1d36*(vin(l+2,:)-vin(l-2,:))
      end do
      !
      if(ns==642) then
        ! ns==642: 2-4-6-6-6-...-6-6-6-4-2
        vout(dim-1,:)=0.75d0*( vin(dim,:)  -vin(dim-2,:))
        vout(dim,:)  =2.d0*  (-vin(dim-1,:)+vin(dim,:))
      elseif(ns==644) then
        ! ns==644: 4-4-6-6-6-...-6-6-6-4-4
        vout(dim-1,:)=0.75d0*( vin(dim,:)  -vin(dim-2,:))
        vout(dim,:)=  0.5d0 *( vin(dim+1,:)-vin(dim-1,:))
      end if
      !
    elseif(ntype==3) then
      ! inner block
      vout(0,:)=0.75d0* (vin(1,:)-vin(-1,:)) -                 &
                0.15d0* (vin(2,:)-vin(-2,:))+                  &
                num1d60*(vin(3,:)-vin(-3,:))
      do l=1,dim-1
        vout(l,:)=num7d9* (vin(l+1,:)-vin(l-1,:))+             &
                  num1d36*(vin(l+2,:)-vin(l-2,:))
      end do
        vout(dim,:)=0.75d0* (vin(dim+1,:)-vin(dim-1,:))-       &
                    0.15d0* (vin(dim+2,:)-vin(dim-2,:))+       &
                    num1d60*(vin(dim+3,:)-vin(dim-3,:))
     
    elseif(ntype==4) then
      ! the block with boundary at i=0 and i=im
      !
      if(ns==642) then
        ! ns==642: 2-4-6-6-6-...-6-6-6-4-2
        vout(0,:)=2.d0*  (-vin(0,:)+vin(1,:))
        vout(1,:)=0.75d0*( vin(2,:)-vin(0,:))
        !
      elseif(ns==644) then
        ! ns==644: 4-4-6-6-6-...-6-6-6-4-4
        vout(0,:)=0.5d0*( vin(1,:)-vin(-1,:))
        vout(1,:)=0.75d0*( vin(2,:)-vin(0,:))
        !
      end if
      !
      ! first order deritive
      do l=2,dim-2
        vout(l,:)=num7d9* (vin(l+1,:)-vin(l-1,:))+               &
                  num1d36*(vin(l+2,:)-vin(l-2,:))
      end do
      !
      if(ns==642) then
        ! ns==642: 2-4-6-6-6-...-6-6-6-4-2
        vout(dim-1,:)=0.75d0*( vin(dim,:)  -vin(dim-2,:))
        vout(dim,:)  =2.d0*  (-vin(dim-1,:)+vin(dim,:))
      elseif(ns==644) then
        ! ns==644: 4-4-6-6-6-...-6-6-6-4-4
        vout(dim-1,:)=0.75d0*( vin(dim,:)  -vin(dim-2,:))
        vout(dim,:)  =0.5d0 *( vin(dim+1,:)-vin(dim-1,:))
      end if
      !
    else
      print*,ntype
      print*, ' !! ntype error in subroutine ptds_rhs !'
      stop
    end if
    !
    if(present(timerept) .and. timerept) then
      !
      subtime=subtime+ptime()-time_beg
      !
      if(lio .and. lreport .and. ltimrpt) call timereporter(routine='ptds2d_rhs', &
                                             timecost=subtime,      &
                                              message='rhs of compact scheme')
    endif
    !
    return
    !
  end function ptds2d_rhs_firstcolum
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! This subroutine is used to calculate the RHS of tridiagonal martix 
  ! for A*x=b.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function ptds3d_rhs(vin,dim,ns,ntype,ncolm1,ncolm2,timerept) result(vout)
    !
    integer,intent(in) :: dim,ns,ntype,ncolm1,ncolm2
    real(8),intent(in) :: vin(1:ncolm1,1:ncolm2,-hm:dim+hm)
    real(8) :: vout(1:ncolm1,1:ncolm2,0:dim)
    logical,intent(in),optional :: timerept
    !
    ! local data
    integer :: l
    real(8) :: var1,var2,var3
    !
    real(8) :: time_beg
    real(8),save :: subtime=0.d0
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! ns=601: 6-o compact central scheme with
    ! boundary scheme: 2-4-6.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    if(present(timerept) .and. timerept) time_beg=ptime() 
    !
    if(ntype==1) then
      ! the block with boundary at i==0
      if(ns==642) then
        ! ns==642: 2-4-6-6-6-...-6-6-6-4-2
        vout(:,:,0)=2.d0*  (-vin(:,:,0)+vin(:,:,1))
        vout(:,:,1)=0.75d0*( vin(:,:,2)-vin(:,:,0))
        !
      elseif(ns==644) then
        ! ns==644: 4-4-6-6-6-...-6-6-6-4-4
        vout(:,:,0)= 0.5d0*( vin(:,:,1)-vin(:,:,-1))
        vout(:,:,1)=0.75d0*( vin(:,:,2)-vin(:,:,0))
        !
      end if
      !
      ! first order deritive
      do l=2,dim-1
        vout(:,:,l)=num7d9* (vin(:,:,l+1)-vin(:,:,l-1))+               &
                    num1d36*(vin(:,:,l+2)-vin(:,:,l-2))
      end do
      vout(:,:,dim)=0.75d0 *(vin(:,:,dim+1)-vin(:,:,dim-1))-           &
                    0.15d0 *(vin(:,:,dim+2)-vin(:,:,dim-2))+           &
                    num1d60*(vin(:,:,dim+3)-vin(:,:,dim-3))
      !
    elseif(ntype==2) then
      ! the block with boundary at i==im
      ! first order deritive
      vout(:,:,0)=0.75d0* (vin(:,:,1)-vin(:,:,-1)) -                   &
                  0.15d0* (vin(:,:,2)-vin(:,:,-2)) +                   &
                  num1d60*(vin(:,:,3)-vin(:,:,-3))
      do l=1,dim-2
        vout(:,:,l)=num7d9* (vin(:,:,l+1)-vin(:,:,l-1))+                           &
                    num1d36*(vin(:,:,l+2)-vin(:,:,l-2))
      end do
      !
      if(ns==642) then
        ! ns==642: 2-4-6-6-6-...-6-6-6-4-2
        vout(:,:,dim-1)=0.75d0*( vin(:,:,dim)  -vin(:,:,dim-2))
        vout(:,:,dim)  =2.d0*  (-vin(:,:,dim-1)+vin(:,:,dim))
      elseif(ns==644) then
        ! ns==644: 4-4-6-6-6-...-6-6-6-4-4
        vout(:,:,dim-1)=0.75d0*( vin(:,:,dim)  -vin(:,:,dim-2))
        vout(:,:,dim)=  0.5d0 *( vin(:,:,dim+1)-vin(:,:,dim-1))
      end if
      !
    elseif(ntype==3) then
      ! inner block
      vout(:,:,0)=0.75d0* (vin(:,:,1)-vin(:,:,-1)) -                   &
                  0.15d0* (vin(:,:,2)-vin(:,:,-2))+                    &
                  num1d60*(vin(:,:,3)-vin(:,:,-3))
      do l=1,dim-1
        vout(:,:,l)=num7d9* (vin(:,:,l+1)-vin(:,:,l-1))+                           &
                num1d36*(vin(:,:,l+2)-vin(:,:,l-2))
      end do
        vout(:,:,dim)=0.75d0* (vin(:,:,dim+1)-vin(:,:,dim-1))-         &
                      0.15d0* (vin(:,:,dim+2)-vin(:,:,dim-2))+         &
                      num1d60*(vin(:,:,dim+3)-vin(:,:,dim-3))
     
    elseif(ntype==4) then
      ! the block with boundary at i=0 and i=im
      !
      if(ns==642) then
        ! ns==642: 2-4-6-6-6-...-6-6-6-4-2
        vout(:,:,0)=2.d0*  (-vin(:,:,0)+vin(:,:,1))
        vout(:,:,1)=0.75d0*( vin(:,:,2)-vin(:,:,0))
        !
      elseif(ns==644) then
        ! ns==644: 4-4-6-6-6-...-6-6-6-4-4
        vout(:,:,0)=0.5d0*( vin(:,:,1)-vin(:,:,-1))
        vout(:,:,1)=0.75d0*( vin(:,:,2)-vin(:,:,0))
        !
      end if
      !
      ! first order deritive
      do l=2,dim-2
        vout(:,:,l)=num7d9* (vin(:,:,l+1)-vin(:,:,l-1))+               &
                    num1d36*(vin(:,:,l+2)-vin(:,:,l-2))
      end do
      !
      if(ns==642) then
        ! ns==642: 2-4-6-6-6-...-6-6-6-4-2
        vout(:,:,dim-1)=0.75d0*( vin(:,:,dim)  -vin(:,:,dim-2))
        vout(:,:,dim)  =2.d0*  (-vin(:,:,dim-1)+vin(:,:,dim))
      elseif(ns==644) then
        ! ns==644: 4-4-6-6-6-...-6-6-6-4-4
        vout(:,:,dim-1)=0.75d0*( vin(:,:,dim)  -vin(:,:,dim-2))
        vout(:,:,dim)  =0.5d0 *( vin(:,:,dim+1)-vin(:,:,dim-1))
      end if
      !
    else
      print*,ntype
      print*, ' !! ntype error in subroutine ptds_rhs !'
      stop
    end if
    !
    if(present(timerept) .and. timerept) then
      !
      subtime=subtime+ptime()-time_beg
      !
      if(lio .and. lreport .and. ltimrpt) call timereporter(routine='ptds3d_rhs', &
                                             timecost=subtime,      &
                                              message='rhs of compact scheme')
    endif
    !
    return
    !
  end function ptds3d_rhs
  !
  function ptds4d_rhs(vin,dim,ns,ntype,ncolm1,ncolm2,ncolm3,timerept) result(vout)
    !
    integer,intent(in) :: dim,ns,ntype,ncolm1,ncolm2,ncolm3
    real(8),intent(in) :: vin(1:ncolm1,1:ncolm2,1:ncolm3,-hm:dim+hm)
    real(8) :: vout(1:ncolm1,1:ncolm2,1:ncolm3,0:dim)
    logical,intent(in),optional :: timerept
    !
    ! local data
    integer :: l
    real(8) :: var1,var2,var3
    !
    real(8) :: time_beg
    real(8),save :: subtime=0.d0
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! ns=601: 6-o compact central scheme with
    ! boundary scheme: 2-4-6.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    if(present(timerept) .and. timerept) time_beg=ptime() 
    !
    if(ntype==1) then
      ! the block with boundary at i==0
      if(ns==642) then
        ! ns==642: 2-4-6-6-6-...-6-6-6-4-2
        vout(:,:,:,0)=2.d0*  (-vin(:,:,:,0)+vin(:,:,:,1))
        vout(:,:,:,1)=0.75d0*( vin(:,:,:,2)-vin(:,:,:,0))
        !
      elseif(ns==644) then
        ! ns==644: 4-4-6-6-6-...-6-6-6-4-4
        vout(:,:,:,0)= 0.5d0*( vin(:,:,:,1)-vin(:,:,:,-1))
        vout(:,:,:,1)=0.75d0*( vin(:,:,:,2)-vin(:,:,:,0))
        !
      end if
      !
      ! first order deritive
      do l=2,dim-1
        vout(:,:,:,l)=num7d9* (vin(:,:,:,l+1)-vin(:,:,:,l-1))+         &
                    num1d36*(vin(:,:,:,l+2)-vin(:,:,:,l-2))
      end do
      vout(:,:,:,dim)=0.75d0 *(vin(:,:,:,dim+1)-vin(:,:,:,dim-1))-     &
                    0.15d0 *(vin(:,:,:,dim+2)-vin(:,:,:,dim-2))+       &
                    num1d60*(vin(:,:,:,dim+3)-vin(:,:,:,dim-3))
      !
    elseif(ntype==2) then
      ! the block with boundary at i==im
      ! first order deritive
      vout(:,:,:,0)=0.75d0* (vin(:,:,:,1)-vin(:,:,:,-1)) -             &
                  0.15d0* (vin(:,:,:,2)-vin(:,:,:,-2)) +               &
                  num1d60*(vin(:,:,:,3)-vin(:,:,:,-3))
      do l=1,dim-2
        vout(:,:,:,l)=num7d9* (vin(:,:,:,l+1)-vin(:,:,:,l-1))+                           &
                    num1d36*(vin(:,:,:,l+2)-vin(:,:,:,l-2))
      end do
      !
      if(ns==642) then
        ! ns==642: 2-4-6-6-6-...-6-6-6-4-2
        vout(:,:,:,dim-1)=0.75d0*( vin(:,:,:,dim)  -vin(:,:,:,dim-2))
        vout(:,:,:,dim)  =2.d0*  (-vin(:,:,:,dim-1)+vin(:,:,:,dim))
      elseif(ns==644) then
        ! ns==644: 4-4-6-6-6-...-6-6-6-4-4
        vout(:,:,:,dim-1)=0.75d0*( vin(:,:,:,dim)  -vin(:,:,:,dim-2))
        vout(:,:,:,dim)=  0.5d0 *( vin(:,:,:,dim+1)-vin(:,:,:,dim-1))
      end if
      !
    elseif(ntype==3) then
      ! inner block
      vout(:,:,:,0)=0.75d0* (vin(:,:,:,1)-vin(:,:,:,-1)) -             &
                  0.15d0* (vin(:,:,:,2)-vin(:,:,:,-2))+                &
                  num1d60*(vin(:,:,:,3)-vin(:,:,:,-3))
      do l=1,dim-1
        vout(:,:,:,l)=num7d9* (vin(:,:,:,l+1)-vin(:,:,:,l-1))+                           &
                num1d36*(vin(:,:,:,l+2)-vin(:,:,:,l-2))
      end do
        vout(:,:,:,dim)=0.75d0* (vin(:,:,:,dim+1)-vin(:,:,:,dim-1))-   &
                      0.15d0* (vin(:,:,:,dim+2)-vin(:,:,:,dim-2))+     &
                      num1d60*(vin(:,:,:,dim+3)-vin(:,:,:,dim-3))
     
    elseif(ntype==4) then
      ! the block with boundary at i=0 and i=im
      !
      if(ns==642) then
        ! ns==642: 2-4-6-6-6-...-6-6-6-4-2
        vout(:,:,:,0)=2.d0*  (-vin(:,:,:,0)+vin(:,:,:,1))
        vout(:,:,:,1)=0.75d0*( vin(:,:,:,2)-vin(:,:,:,0))
        !
      elseif(ns==644) then
        ! ns==644: 4-4-6-6-6-...-6-6-6-4-4
        vout(:,:,:,0)=0.5d0*( vin(:,:,:,1)-vin(:,:,:,-1))
        vout(:,:,:,1)=0.75d0*( vin(:,:,:,2)-vin(:,:,:,0))
        !
      end if
      !
      ! first order deritive
      do l=2,dim-2
        vout(:,:,:,l)=num7d9* (vin(:,:,:,l+1)-vin(:,:,:,l-1))+         &
                    num1d36*(vin(:,:,:,l+2)-vin(:,:,:,l-2))
      end do
      !
      if(ns==642) then
        ! ns==642: 2-4-6-6-6-...-6-6-6-4-2
        vout(:,:,:,dim-1)=0.75d0*( vin(:,:,:,dim)  -vin(:,:,:,dim-2))
        vout(:,:,:,dim)  =2.d0*  (-vin(:,:,:,dim-1)+vin(:,:,:,dim))
      elseif(ns==644) then
        ! ns==644: 4-4-6-6-6-...-6-6-6-4-4
        vout(:,:,:,dim-1)=0.75d0*( vin(:,:,:,dim)  -vin(:,:,:,dim-2))
        vout(:,:,:,dim)  =0.5d0 *( vin(:,:,:,dim+1)-vin(:,:,:,dim-1))
      end if
      !
    else
      print*,ntype
      print*, ' !! ntype error in subroutine ptds_rhs !'
      stop
    end if
    !
    if(present(timerept) .and. timerept) then
      !
      subtime=subtime+ptime()-time_beg
      !
      if(lio .and. lreport .and. ltimrpt) call timereporter(routine='ptds4d_rhs', &
                                             timecost=subtime,      &
                                              message='rhs of compact scheme')
    endif
    !
    return
    !
  end function ptds4d_rhs


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
  function ptds_cal(bd,af,cc,dim,ntype,timerept) result(xd)
    !
    integer,intent(in) :: dim,ntype
    real(8),intent(in) :: af(3),bd(0:dim),cc(1:2,0:dim)
    real(8) :: xd(0:dim)
    !
    logical,intent(in),optional :: timerept
    !
    ! local data
    integer :: l
    real(8),allocatable,dimension(:) :: yd
    real(8) :: time_beg
    real(8),save :: subtime=0.d0
    !
    if(present(timerept) .and. timerept) time_beg=ptime() 
    !
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
      do l=2,dim-1
        yd(l)=(bd(l)-af(3)*yd(l-1))*cc(1,l)
      end do
      yd(dim)=bd(dim)*cc(1,dim)
    elseif(ntype==2) then
      ! the block with boundary at i==im
      !
      yd(0)=bd(0)*cc(1,0)
      do l=1,dim-2
        yd(l)=(bd(l)-af(3)*yd(l-1))*cc(1,l)
      end do
      yd(dim-1)=(bd(dim-1)-af(2)*yd(dim-2))*cc(1,dim-1)
      yd(dim)=(bd(dim)-af(1)*yd(dim-1))*cc(1,dim)
      !
    elseif(ntype==3) then
      ! inner block
      yd(0)=bd(0)*cc(1,0)
      do l=1,dim-1
        yd(l)=(bd(l)-af(3)*yd(l-1))*cc(1,l)
      end do
      yd(dim)=bd(dim)*cc(1,dim)
    elseif(ntype==4) then
      ! the block with boundary at i=0 and i=im
      !
      yd(0)=bd(0)*cc(1,0)
      yd(1)=(bd(1)-af(2)*yd(0))*cc(1,1)
      do l=2,dim-2
        yd(l)=(bd(l)-af(3)*yd(l-1))*cc(1,l)
      end do
      yd(dim-1)=(bd(dim-1)-af(2)*yd(dim-2))*cc(1,dim-1)
      yd(dim)=(bd(dim)-af(1)*yd(dim-1))*cc(1,dim)
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
    if(present(timerept) .and. timerept) then
      !
      subtime=subtime+ptime()-time_beg
      !
      if(lio .and. lreport .and. ltimrpt) call timereporter(routine='ptds_cal', &
                                              timecost=subtime, &
                                              message='solve compact tridiagonal system')
    endif
    !
    return
    !
  end function ptds_cal
  !!
  function ptds2d_cal_lastcolum(bd,af,cc,dim,ntype,ncolm1,timerept) result(xd)
    !
    integer,intent(in) :: dim,ntype,ncolm1
    real(8),intent(in) :: af(3),bd(1:ncolm1,0:dim),cc(1:2,0:dim)
    real(8) :: xd(1:ncolm1,0:dim)
    !
    logical,intent(in),optional :: timerept
    !
    ! local data
    integer :: l
    real(8),allocatable :: yd(:,:)
    real(8) :: time_beg
    real(8),save :: subtime=0.d0
    !
    if(present(timerept) .and. timerept) time_beg=ptime() 
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! af(3): input dat
    ! bd: input array
    ! xd: output array
    ! cc: input array
    ! dim: input dat
    ! l, yd: temporary variable
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    allocate ( yd(1:ncolm1,0:dim) )
    !
    if(ntype==1) then
      ! the block with boundary at i==0
      !
      yd(:,0)=bd(:,0)*cc(1,0)
      yd(:,1)=(bd(:,1)-af(2)*yd(:,0))*cc(1,1)
      do l=2,dim-1
        yd(:,l)=(bd(:,l)-af(3)*yd(:,l-1))*cc(1,l)
      end do
      yd(:,dim)=bd(:,dim)*cc(1,dim)
    elseif(ntype==2) then
      ! the block with boundary at i==im
      !
      yd(:,0)=bd(:,0)*cc(1,0)
      do l=1,dim-2
        yd(:,l)=(bd(:,l)-af(3)*yd(:,l-1))*cc(1,l)
      end do
      yd(:,dim-1)=(bd(:,dim-1)-af(2)*yd(:,dim-2))*cc(1,dim-1)
      yd(:,dim)=(bd(:,dim)-af(1)*yd(:,dim-1))*cc(1,dim)
      !
    elseif(ntype==3) then
      ! inner block
      yd(:,0)=bd(:,0)*cc(1,0)
      do l=1,dim-1
        yd(:,l)=(bd(:,l)-af(3)*yd(:,l-1))*cc(1,l)
      end do
      yd(:,dim)=bd(:,dim)*cc(1,dim)
    elseif(ntype==4) then
      ! the block with boundary at i=0 and i=im
      !
      yd(:,0)=bd(:,0)*cc(1,0)
      yd(:,1)=(bd(:,1)-af(2)*yd(:,0))*cc(1,1)
      do l=2,dim-2
        yd(:,l)=(bd(:,l)-af(3)*yd(:,l-1))*cc(1,l)
      end do
      yd(:,dim-1)=(bd(:,dim-1)-af(2)*yd(:,dim-2))*cc(1,dim-1)
      yd(:,dim)=(bd(:,dim)-af(1)*yd(:,dim-1))*cc(1,dim)
    else
      print*, ' !! error in subroutine ptds_cal !'
      stop
    end if
    !
    xd(:,dim)=yd(:,dim)
    do l=dim-1,0,-1
      xd(:,l)=yd(:,l)-cc(2,l)*xd(:,l+1)
    end do
    !
    deallocate( yd )
    !
    if(present(timerept) .and. timerept) then
      !
      subtime=subtime+ptime()-time_beg
      !
      if(lio .and. lreport .and. ltimrpt) call timereporter(routine='ptds2d_cal', &
                                              timecost=subtime,     &
                                              message='solve compact tridiagonal system')
    endif
    !
    return
    !
  end function ptds2d_cal_lastcolum
  !
  function ptds2d_cal_firstcolum(bd,af,cc,dim,ntype,ncolm1,timerept) result(xd)
    !
    integer,intent(in) :: dim,ntype,ncolm1
    real(8),intent(in) :: af(3),cc(1:2,0:dim)
    real(8),intent(in) :: bd(0:dim,1:ncolm1)
    real(8) :: xd(0:dim,1:ncolm1)
    !
    logical,intent(in),optional :: timerept
    !
    ! local data
    integer :: l
    real(8),allocatable :: yd(:,:)
    real(8) :: time_beg
    real(8),save :: subtime=0.d0
    !
    if(present(timerept) .and. timerept) time_beg=ptime() 
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! af(3): input dat
    ! bd: input array
    ! xd: output array
    ! cc: input array
    ! dim: input dat
    ! l, yd: temporary variable
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    allocate ( yd(0:dim,1:ncolm1) )
    !
    if(ntype==1) then
      ! the block with boundary at i==0
      !
      yd(0,:)=bd(0,:)*cc(1,0)
      yd(1,:)=(bd(1,:)-af(2)*yd(0,:))*cc(1,1)
      do l=2,dim-1
        yd(l,:)=(bd(l,:)-af(3)*yd(l-1,:))*cc(1,l)
      end do
      yd(dim,:)=bd(dim,:)*cc(1,dim)
    elseif(ntype==2) then
      ! the block with boundary at i==im
      !
      yd(0,:)=bd(0,:)*cc(1,0)
      do l=1,dim-2
        yd(l,:)=(bd(l,:)-af(3)*yd(l-1,:))*cc(1,l)
      end do
      yd(dim-1,:)=(bd(dim-1,:)-af(2)*yd(dim-2,:))*cc(1,dim-1)
      yd(dim,:)=(bd(dim,:)-af(1)*yd(dim-1,:))*cc(1,dim)
      !
    elseif(ntype==3) then
      ! inner block
      yd(0,:)=bd(0,:)*cc(1,0)
      do l=1,dim-1
        yd(l,:)=(bd(l,:)-af(3)*yd(l-1,:))*cc(1,l)
      end do
      yd(dim,:)=bd(dim,:)*cc(1,dim)
    elseif(ntype==4) then
      ! the block with boundary at i=0 and i=im
      !
      yd(0,:)=bd(0,:)*cc(1,0)
      yd(1,:)=(bd(1,:)-af(2)*yd(0,:))*cc(1,1)
      do l=2,dim-2
        yd(l,:)=(bd(l,:)-af(3)*yd(l-1,:))*cc(1,l)
      end do
      yd(dim-1,:)=(bd(dim-1,:)-af(2)*yd(dim-2,:))*cc(1,dim-1)
      yd(dim,:)=(bd(dim,:)-af(1)*yd(dim-1,:))*cc(1,dim)
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
    if(present(timerept) .and. timerept) then
      !
      subtime=subtime+ptime()-time_beg
      !
      if(lio .and. lreport .and. ltimrpt) call timereporter(routine='ptds2d_cal', &
                                              timecost=subtime,     &
                                              message='solve compact tridiagonal system')
    endif
    !
    return
    !
  end function ptds2d_cal_firstcolum
  !
  function ptds3d_cal(bd,af,cc,dim,ntype,ncolm1,ncolm2,timerept) result(xd)
    !
    integer,intent(in) :: dim,ntype,ncolm1,ncolm2
    real(8),intent(in) :: af(3),bd(1:ncolm1,1:ncolm2,0:dim),cc(1:2,0:dim)
    real(8) :: xd(1:ncolm1,1:ncolm2,0:dim)
    !
    logical,intent(in),optional :: timerept
    !
    ! local data
    integer :: l
    real(8),allocatable :: yd(:,:,:)
    real(8) :: time_beg
    real(8),save :: subtime=0.d0
    !
    if(present(timerept) .and. timerept) time_beg=ptime() 
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! af(3): input dat
    ! bd: input array
    ! xd: output array
    ! cc: input array
    ! dim: input dat
    ! l, yd: temporary variable
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    allocate ( yd(1:ncolm1,1:ncolm2,0:dim) )
    !
    if(ntype==1) then
      ! the block with boundary at i==0
      !
      yd(:,:,0)=bd(:,:,0)*cc(1,0)
      yd(:,:,1)=(bd(:,:,1)-af(2)*yd(:,:,0))*cc(1,1)
      do l=2,dim-1
        yd(:,:,l)=(bd(:,:,l)-af(3)*yd(:,:,l-1))*cc(1,l)
      end do
      yd(:,:,dim)=bd(:,:,dim)*cc(1,dim)
    elseif(ntype==2) then
      ! the block with boundary at i==im
      !
      yd(:,:,0)=bd(:,:,0)*cc(1,0)
      do l=1,dim-2
        yd(:,:,l)=(bd(:,:,l)-af(3)*yd(:,:,l-1))*cc(1,l)
      end do
      yd(:,:,dim-1)=(bd(:,:,dim-1)-af(2)*yd(:,:,dim-2))*cc(1,dim-1)
      yd(:,:,dim)=(bd(:,:,dim)-af(1)*yd(:,:,dim-1))*cc(1,dim)
      !
    elseif(ntype==3) then
      ! inner block
      yd(:,:,0)=bd(:,:,0)*cc(1,0)
      do l=1,dim-1
        yd(:,:,l)=(bd(:,:,l)-af(3)*yd(:,:,l-1))*cc(1,l)
      end do
      yd(:,:,dim)=bd(:,:,dim)*cc(1,dim)
    elseif(ntype==4) then
      ! the block with boundary at i=0 and i=im
      !
      yd(:,:,0)=bd(:,:,0)*cc(1,0)
      yd(:,:,1)=(bd(:,:,1)-af(2)*yd(:,:,0))*cc(1,1)
      do l=2,dim-2
        yd(:,:,l)=(bd(:,:,l)-af(3)*yd(:,:,l-1))*cc(1,l)
      end do
      yd(:,:,dim-1)=(bd(:,:,dim-1)-af(2)*yd(:,:,dim-2))*cc(1,dim-1)
      yd(:,:,dim)=(bd(:,:,dim)-af(1)*yd(:,:,dim-1))*cc(1,dim)
    else
      print*, ' !! error in subroutine ptds_cal !'
      stop
    end if
    !
    xd(:,:,dim)=yd(:,:,dim)
    do l=dim-1,0,-1
      xd(:,:,l)=yd(:,:,l)-cc(2,l)*xd(:,:,l+1)
    end do
    !
    deallocate( yd )
    !
    if(present(timerept) .and. timerept) then
      !
      subtime=subtime+ptime()-time_beg
      !
      if(lio .and. lreport .and. ltimrpt) call timereporter(routine='ptds3d_cal', &
                                              timecost=subtime,     &
                                              message='solve compact tridiagonal system')
    endif
    !
    return
    !
  end function ptds3d_cal
  !!
  function ptds4d_cal(bd,af,cc,dim,ntype,ncolm1,ncolm2,ncolm3,timerept) result(xd)
    !
    integer,intent(in) :: dim,ntype,ncolm1,ncolm2,ncolm3
    real(8),intent(in) :: af(3),cc(1:2,0:dim)
    real(8),intent(in) :: bd(1:ncolm1,1:ncolm2,1:ncolm3,0:dim)
    real(8) :: xd(1:ncolm1,1:ncolm2,1:ncolm3,0:dim)
    !
    logical,intent(in),optional :: timerept
    !
    ! local data
    integer :: l
    real(8),allocatable :: yd(:,:,:,:)
    real(8) :: time_beg
    real(8),save :: subtime=0.d0
    !
    if(present(timerept) .and. timerept) time_beg=ptime() 
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! af(3): input dat
    ! bd: input array
    ! xd: output array
    ! cc: input array
    ! dim: input dat
    ! l, yd: temporary variable
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    allocate ( yd(1:ncolm1,1:ncolm2,1:ncolm3,0:dim) )
    !
    if(ntype==1) then
      ! the block with boundary at i==0
      !
      yd(:,:,:,0)=bd(:,:,:,0)*cc(1,0)
      yd(:,:,:,1)=(bd(:,:,:,1)-af(2)*yd(:,:,:,0))*cc(1,1)
      do l=2,dim-1
        yd(:,:,:,l)=(bd(:,:,:,l)-af(3)*yd(:,:,:,l-1))*cc(1,l)
      end do
      yd(:,:,:,dim)=bd(:,:,:,dim)*cc(1,dim)
    elseif(ntype==2) then
      ! the block with boundary at i==im
      !
      yd(:,:,:,0)=bd(:,:,:,0)*cc(1,0)
      do l=1,dim-2
        yd(:,:,:,l)=(bd(:,:,:,l)-af(3)*yd(:,:,:,l-1))*cc(1,l)
      end do
      yd(:,:,:,dim-1)=(bd(:,:,:,dim-1)-af(2)*yd(:,:,:,dim-2))*cc(1,dim-1)
      yd(:,:,:,dim)=(bd(:,:,:,dim)-af(1)*yd(:,:,:,dim-1))*cc(1,dim)
      !
    elseif(ntype==3) then
      ! inner block
      yd(:,:,:,0)=bd(:,:,:,0)*cc(1,0)
      do l=1,dim-1
        yd(:,:,:,l)=(bd(:,:,:,l)-af(3)*yd(:,:,:,l-1))*cc(1,l)
      end do
      yd(:,:,:,dim)=bd(:,:,:,dim)*cc(1,dim)
    elseif(ntype==4) then
      ! the block with boundary at i=0 and i=im
      !
      yd(:,:,:,0)=bd(:,:,:,0)*cc(1,0)
      yd(:,:,:,1)=(bd(:,:,:,1)-af(2)*yd(:,:,:,0))*cc(1,1)
      do l=2,dim-2
        yd(:,:,:,l)=(bd(:,:,:,l)-af(3)*yd(:,:,:,l-1))*cc(1,l)
      end do
      yd(:,:,:,dim-1)=(bd(:,:,:,dim-1)-af(2)*yd(:,:,:,dim-2))*cc(1,dim-1)
      yd(:,:,:,dim)=(bd(:,:,:,dim)-af(1)*yd(:,:,:,dim-1))*cc(1,dim)
    else
      print*, ' !! error in subroutine ptds_cal !'
      stop
    end if
    !
    xd(:,:,:,dim)=yd(:,:,:,dim)
    do l=dim-1,0,-1
      xd(:,:,:,l)=yd(:,:,:,l)-cc(2,l)*xd(:,:,:,l+1)
    end do
    !
    deallocate( yd )
    !
    if(present(timerept) .and. timerept) then
      !
      subtime=subtime+ptime()-time_beg
      !
      if(lio .and. lreport .and. ltimrpt) call timereporter(routine='ptds4d_cal', &
                                              timecost=subtime,     &
                                              message='solve compact tridiagonal system')
    endif
    !
    return
    !
    return
    !
  end function ptds4d_cal
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
    subroutine tridiagonal_thomas_proprocess(a,c,ac)
        
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

    end subroutine tridiagonal_thomas_proprocess
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