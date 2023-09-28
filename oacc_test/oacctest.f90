module comvardef
  !
  implicit none
  !
  integer :: im=256,jm=256,km=256,hm=5
  real(8),parameter :: pi=4.d0*atan(1.0_8),                            &
                       num1d35 =1.d0/35.d0,  num1d3  =1.d0/3.d0,       &
                       num2d3  =2.d0/3.d0,   num1d24 =1.d0/24.d0,      &
                       num4d3  =4.d0/3.d0,   num1d6  =1.d0/6.d0,       &
                       num1d12 =1.d0/12.d0,  num7d12 =7.d0/12.d0,      &
                       num7d9  =7.d0/9.d0,   num1d36 =1.d0/36.d0,      &
                       num1d60 =1.d0/60.d0,  num65d3 =65.d0/3.d0,      &
                       num20d3 =20.d0/3.d0,  num1d11 =1.d0/11.d0,      &
                       num25d12=25.d0/12.d0, num11d6 =11.d0/6.d0,      &
                       num1d840=1.d0/840.d0, num13d60=13.d0/60.d0,     &
                       num1d30 =1.d0/30.d0,  num47d60=47.d0/60.d0,     &
                       num5d6  =5.d0/6.d0,   num1d18 =1.d0/18.d0,      &
                       num19d18=19.d0/18.d0, num5d9  =5.d0/9.d0,       &
                       num9d36 =9.d0/36.d0
  !
  contains
  !
end module comvardef
!
module test
  !
  use comvardef
  !
  implicit none
  !
  contains
  !
  subroutine arrayoperate
    !
    integer :: ncolm,i,j,k,n
    real(8) :: var1,start,finish
    real(8),allocatable,dimension(:,:,:,:) :: vtest
    !
    ncolm=12
    !
    allocate(vtest(-hm:im+hm,-hm:jm+hm,-hm:km+hm,1:ncolm))
    ! allocate(dvtes(0:im,0:jm,0:km,1:ncolm,1:3))
    !
    ! a random number assign
    call cpu_time(start)
    !$acc parallel loop collapse(4)
    do n=1,ncolm
    do k=0,km
    do j=0,jm
    do i=0,im
      !
      ! call random_number(var1)
      vtest(i,j,k,n)=sin(dble(i)*2.d0*pi)
      !
    enddo
    enddo
    enddo
    enddo
    call cpu_time(finish)
    print*,'i,j,k,n - number assign', (finish - start)
    !
    call cpu_time(start)
    do n=1,ncolm
    !$acc parallel loop collapse(3)
    do k=0,km
    do j=0,jm
    do i=0,im
      !
      ! call random_number(var1)
      vtest(i,j,k,n)=sin(dble(i)*2.d0*pi)
      !
    enddo
    enddo
    enddo
    enddo
    call cpu_time(finish)
    print*,'i,j,k - number assign', (finish - start)
    !
    call cpu_time(start)
    do n=1,ncolm
    do k=0,km
    !$acc parallel loop collapse(2)
    do j=0,jm
    do i=0,im
      !
      ! call random_number(var1)
      vtest(i,j,k,n)=sin(dble(i)*2.d0*pi)
      !
    enddo
    enddo
    enddo
    enddo
    call cpu_time(finish)
    print*,'i,j - number assign', (finish - start)
    !
    call cpu_time(start)
    do n=1,ncolm
    do k=0,km
    do j=0,jm
    !$acc parallel loop collapse(1)
    do i=0,im
      !
      ! call random_number(var1)
      vtest(i,j,k,n)=sin(dble(i)*2.d0*pi)
      !
    enddo
    enddo
    enddo
    enddo
    call cpu_time(finish)
    print*,'i - number assign', (finish - start)
    !
  end subroutine arrayoperate
  !
  subroutine gradientcal
    !
    integer :: ncolm,i,j,k,n,npdci,nscheme
    real(8) :: var1,time_beg,start,finish
    real(8),allocatable,dimension(:,:,:,:) :: vtest
    real(8),allocatable,dimension(:,:,:,:,:) :: dvtes
    real(8),allocatable :: f1(:,:),df1(:,:)
    real(8),allocatable :: alfa(:),ci(:,:)
    !
    ncolm=12
    !
    allocate(vtest(-hm:im+hm,-hm:jm+hm,-hm:km+hm,1:ncolm))
    allocate(dvtes(0:im,0:jm,0:km,1:ncolm,1:3))
    !
    nscheme=642
    npdci=3
    allocate(alfa(3))
    alfa(3)=num1d3
    alfa(2)=0.25d0
    alfa(1)=1.d0
    !
    call ptds_ini(ci,alfa,im,npdci)
    !
    call cpu_time(time_beg)
    !
    call cpu_time(start)
    !
    do n=1,ncolm
    do k=0,km
    do j=0,jm
    do i=0,im
      !
      ! call random_number(var1)
      vtest(i,j,k,n)=sin(dble(i)/dble(im)*2.d0*pi)
      !
    enddo
    enddo
    enddo
    enddo
    !
    call cpu_time(finish)
    !
    print*,' gradient - i dir: initial data generation: ', (finish - start)
    !
    call cpu_time(start)
    !
    allocate(f1(-hm:im+hm,ncolm),df1(0:im,ncolm))
    do k=0,km
    do j=0,jm
      !
      do n=1,ncolm
        f1(:,n) =vtest(:,j,k,n)
        df1(:,n)=ddfc_basic(f1(:,n),nscheme,npdci,im,alfa,ci)
      enddo
      !
    enddo
    enddo
    !
    call cpu_time(finish)
    print*,' gradient - i dir: deritive calculation:   ', (finish - start)
    !
    call cpu_time(start)
    !
    do k=0,km
    do j=0,jm
      !
      do n=1,ncolm
        dvtes(:,j,k,n,1)=df1(:,n)
      enddo
      !
    enddo
    enddo
    deallocate(f1,df1)
    !
    call cpu_time(finish)
    !
    print*,' gradient - i dir: gradient assignment:   ', (finish - start)
    !
    print*,' gradient - i dir: total time cost:       ', (finish - time_beg)
    !
  end subroutine gradientcal
  !
  function ddfc_basic(f,nscheme,ntype,dim,af,cc) result(ddfc)
    !
    ! arguments
    integer,intent(in) :: nscheme,ntype,dim
    real(8),intent(in) :: f(-hm:dim+hm)
    real(8),intent(in) :: af(3),cc(1:2,0:dim)
    real(8) :: ddfc(0:dim)
    !
    ! local data
    integer :: k
    real(8) :: b(0:dim)
    !
    b   =ptds_rhs(f,dim,nscheme,ntype,timerept=.true.)
    ddfc=ptds_cal(b,af,cc,dim,ntype,timerept=.true.)
    !
    return
    !
  end function ddfc_basic
  !
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
    return
    !
  end function ptds_rhs
  !
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
    return
    !
  end function ptds_cal
  !
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
end module test
!+------------------------------------------------------------------+
!| this program is a mini-app to test and optimise ASTR within      |
!| multicore and GPU environment                                    |
!| ==============                                                   |
!| CHANGE RECORD                                                    |
!| -------------                                                    |
!| 27-09-2023: Created by J. Fang @ STFC Daresbury Laboratory       |
!+------------------------------------------------------------------+
program oaccapp
  !
  use test
  !
  implicit none
  !
  ! call arrayoperate
  !
  call gradientcal
  !
end program oaccapp
!+------------------------------------------------------------------+
!| The end of the program oaccapp                                   |
!+------------------------------------------------------------------+