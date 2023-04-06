!+---------------------------------------------------------------------+
!| This module contains subroutines for testing and debug purpose      |
!+---------------------------------------------------------------------+
!| CHANGE RECORD                                                       |
!| -------------                                                       |
!| 05-04-2023  | Created by J. Fang                                    |
!+---------------------------------------------------------------------+
module test
  !
  use constdef
  !
  implicit none
  !
  contains
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is the entrance of the test                       |  
  !+-------------------------------------------------------------------+
  subroutine codetest
    !
    use commvar,   only : testmode
    use parallel,  only : mpistop,bcast
    use cmdefne,   only : readkeyboad
    !
    call readkeyboad(testmode)
    !
    call bcast(testmode)

    !
    if(testmode=='grad') then 
      call gradtest
    elseif(testmode=='enst') then 
      call enstest
    elseif(testmode=='filt') then 
      call filtertest
    else
      return
    endif
    !
    call mpistop
    !
  end subroutine codetest
  !+-------------------------------------------------------------------+
  !| The end of the subroutine codetest                                |  
  !+-------------------------------------------------------------------+
  !
  subroutine enstest

    use commvar,   only : im,jm,km,ia,ja,ka,roinf,uinf,hm,is,ie,js,je,ks,ke
    use commarray, only : x,vel,dvel,rho
    use solver,    only : gradcal
    use statistic, only : enstophycal
    use parallel,  only : qswap,psum,irk,jrk,krk,mpirank

    integer :: i,j,k
    real(8) :: enstrophy,omega(3),omegam

    
    call qswap

    call gradcal

    
    ! enstrophy=0.d0
    ! do k=1,km
    ! do j=1,jm
    ! do i=1,im
    !   ! dx=
    !   omega(1)=dvel(i,j,k,3,2)-dvel(i,j,k,2,3)
    !   omega(2)=dvel(i,j,k,1,3)-dvel(i,j,k,3,1)
    !   omega(3)=dvel(i,j,k,2,1)-dvel(i,j,k,1,2)
    !   omegam=omega(1)*omega(1)+omega(2)*omega(2)+omega(3)*omega(3)
    !   !
    !   enstrophy=enstrophy+rho(i,j,k)*omegam
    !   !
    !   if(mpirank==0 .and. j==1 .and. k==1) then
    !     print*,x(i,j,k,1),dvel(i,j,k,2,1),dvel(i,j,k,1,2)
    !   endif
    !   !
    ! enddo
    ! enddo
    ! enddo
    ! enstrophy=0.5d0*psum(enstrophy)/real(ia*ja*ka,8)
    ! enstrophy=enstrophy/(roinf*(uinf)**2)
    !
    enstrophy=enstophycal()
    print*,' ** enstrophy=',enstrophy
    !
    if(mpirank==0) then
      open(18,file='testout/profilex.dat')
      do i=0,im
        write(18,*)x(i,0,0,1),vel(i,0,0,1),dvel(i,0,0,1,1)
      enddo
      close(18)
      print*,' << profilex.dat'
    endif
    !
  end subroutine enstest
  !
  subroutine gradtest
    !
    use commvar,   only : im,jm,km,npdci,npdcj,npdck,conschm,          &
                          alfa_filter,numq,is,ie
    use commarray, only : x,q,dxi
    use commfunc,  only : ddfc,recons,spafilter10,spafilter6exp
    use bc,        only : boucon
    use parallel,  only : dataswap,mpirankname
    use solver,    only : alfa_con,cci,ccj,cck
    !
    ! local data
    integer :: i,j,k,n
    real(8) :: dx
    real(8),allocatable :: vtest(:,:,:)
    real(8),allocatable :: dq(:),qhp(:),qhm(:),dqref(:)
    !
    ! print*,x(:,0,0,1)
    !
    ! testing ddx
    do k=0,km
    do j=0,jm
    do i=0,im
      !
      q(i,j,k,1)=sin(x(i,j,k,1))
      !
    enddo
    enddo
    enddo
    !
    call dataswap(q)
    !
    allocate(dq(0:im))
    !
    dq(:)=ddfc(q(:,0,0,1),conschm,npdci,im,alfa_con,cci)*dxi(:,0,0,1,1)
    !
    open(18,file='testout/ddx_j=0_k=0.dat')
    write(18,'(3(1X,A15))')'x','dqdx','cosx'
    write(18,'(3(1X,E15.7E3))')(x(i,0,0,1),dq(i),cos(x(i,0,0,1)),i=0,im)
    close(18)
    print*,' << testout/ddx_j=0_k=0.dat'
    !
    deallocate(dq)
    !
    ! testing ddy
    do k=0,km
    do j=0,jm
    do i=0,im
      !
      q(i,j,k,1)=sin(0.5d0*pi*x(i,j,k,2))
      !
    enddo
    enddo
    enddo
    !
    call dataswap(q)
    !
    allocate(dq(0:jm),dqref(0:jm))
    !
    dq(:)=ddfc(q(0,:,0,1),conschm,npdcj,jm,alfa_con,ccj)*dxi(0,0:jm,0,2,2)
    !
    i=0; k=0
    do j=0,jm
      dqref(j)=0.5d0*pi*cos(0.5d0*pi*x(i,j,k,2))
    end do
    !
    open(18,file='testout/ddy_i=0_k=0.dat')
    write(18,'(4(1X,A15))')'y','q','dqdy','dq_ref'
    write(18,'(4(1X,E15.7E3))')(x(0,j,0,2),q(0,j,0,1),dq(j),dqref(j),j=0,jm)
    close(18)
    print*,' << testout/ddy_i=0_k=0.dat'
    !
    deallocate(dq,dqref)
    !
    ! testing ddz
    do k=0,km
    do j=0,jm
    do i=0,im
      !
      q(i,j,k,1)=sin(2.d0*x(i,j,k,3))
      !
    enddo
    enddo
    enddo
    !
    call dataswap(q)
    !
    allocate(dq(0:km))
    !
    dq(:)=ddfc(q(0,0,:,1),conschm,npdck,km,alfa_con,cck)*dxi(0,0,0:km,3,3)
    !
    open(18,file='testout/ddz_i=0_j=0.dat')
    write(18,'(3(1X,A15))')'y','dqdy','cosy'
    write(18,'(3(1X,E15.7E3))')(x(0,0,k,3),dq(k),cos(x(0,0,k,3)),k=0,km)
    close(18)
    print*,' << testout/ddz_i=0_j=0.dat'
    !
    deallocate(dq)
    !
    !
    ! call boucon
    !
    ! allocate(dq(0:im,1:2),qhp(is-1:ie),qhm(is-1:ie))
    ! !
    ! qhp(:)=recons(q(:,0,0,1),conschm,npdci,im,alfa_con,uci,windir='+')
    ! !
    ! qhm(:)=recons(q(:,0,0,1),conschm,npdci,im,alfa_con,bci,windir='-')
    ! !
    ! dx=x(1,0,0,1)-x(0,0,0,1)
    ! !
    ! do i=is,ie
    !   dq(i,1)=(qhp(i)-qhp(i-1))/dx
    ! enddo
    ! !
    ! do i=is,ie
    !   dq(i,2)=(qhm(i)-qhm(i-1))/dx
    ! enddo
    !
    ! dq=ddfc(q(:,1,1,1),conschm,npdci,im,alfa_con,cci)*dxi(:,1,1,1,1)
    
    ! dq=spafilter10(q(:,jm,0,2),npdci,im,alfa_filter,fci)
    !
    ! allocate(dq(0:jm,1:1))
    ! ! dq(:,1)=ddfc(q(0,:,0,2),conschm,npdcj,jm,alfa_con,ccj)*dxi(0,0:jm,0,2,2)
    ! !
    ! ! dq(:,1)=spafilter10(q(1,:,0,3),npdcj,jm,alfa_filter,fcj)
    ! dq(:,1)=spafilter6exp(q(1,:,0,3),npdcj,jm)
    !
    ! allocate(dq(0:km,1:numq))
    ! do n=1,numq
    !   ! dq(:,n)=ddfc(q(1,1,:,n),conschm,npdck,km,alfa_con,cck,lfft=.true.)
    !   dq(:,n)=spafilter10(q(1,1,:,n),npdck,km,alfa_filter,fck,lfft=lfftk)
    ! enddo
    !
    ! do n=1,numq
    !   dq(:,n)=dq(:,n)*dxi(1,1,0:km,3,3)
    ! enddo
    ! !
    ! if(mpirank==0) then
    !   do k=0,km
    !     print*,k,dxi(1,1,k,3,3)
    !   enddo
    ! endif
    !
    ! dq=ddfc(q(1,1,:,1),conschm,npdck,km,alfa_con,cck)/(x(1,1,1,3)-x(1,1,0,3))
    !
    ! open(18,file='testout/profile'//mpirankname//'.dat')
    ! do j=0,jm
    !   write(18,'(3(1X,E15.7E3))')x(1,j,0,2),q(1,j,0,3),dq(j,1)
    ! enddo
    ! close(18)
    ! print*,' << testout/profile',mpirankname,'.dat'
    ! !
    ! deallocate(dq)
    !
  end subroutine gradtest
  !
  subroutine filtertest
    !
    use commvar,   only : im,jm,km,npdci,npdcj,npdck,conschm,          &
                          alfa_filter,numq,is,ie
    use commarray, only : x,q,dxi
    use commfunc,  only : ddfc,recons,spafilter10,spafilter6exp
    use bc,        only : boucon
    use parallel,  only : dataswap,mpirankname
    use solver,    only : alfa_con,fci,fcj,fck
    !
    ! local data
    integer :: i,j,k,n
    real(8) :: dx
    real(8),allocatable :: vtest(:,:,:)
    real(8),allocatable :: dq(:),qhp(:),qhm(:),dqref(:)
    !
    ! print*,x(:,0,0,1)
    !
    ! testing ddx
    do k=0,km
    do j=0,jm
    do i=0,im
      !
      q(i,j,k,1)=sin(10.d0*x(i,j,k,1))
      !
    enddo
    enddo
    enddo
    !
    call dataswap(q)
    !
    allocate(dq(0:im))
    !
    dq(:)=spafilter10(q(:,0,0,1),npdci,im,alfa_filter,fci)
    !
    open(18,file='testout/filterx.dat')
    write(18,'(3(1X,A15))')'x','q','q~'
    write(18,'(3(1X,E15.7E3))')(x(i,0,0,1),q(i,0,0,1),dq(i),i=0,im)
    close(18)
    print*,' << testout/filterx.dat'
    !
    deallocate(dq)
    !
    ! testing ddy
    do k=0,km
    do j=0,jm
    do i=0,im
      !
      q(i,j,k,1)=sin(6.5d0*pi*x(i,j,k,2))
      !
    enddo
    enddo
    enddo
    !
    call dataswap(q)
    !
    allocate(dq(0:jm))
    !
    dq(:)=spafilter10(q(0,:,0,1),npdcj,jm,alfa_filter,fcj)
    !
    open(18,file='testout/filtery.dat')
    write(18,'(3(1X,A15))')'y','q','q~'
    write(18,'(3(1X,E15.7E3))')(x(0,j,0,2),q(0,j,0,1),dq(j),j=0,jm)
    close(18)
    print*,' << testout/filtery.dat'
    !
    deallocate(dq)
    !
    ! testing ddz
    do k=0,km
    do j=0,jm
    do i=0,im
      !
      q(i,j,k,1)=sin(6.d0*x(i,j,k,3))
      !
    enddo
    enddo
    enddo
    !
    call dataswap(q)
    !
    allocate(dq(0:km))
    !
    dq(:)=spafilter10(q(0,0,:,1),npdck,km,alfa_filter,fck)
    !
    open(18,file='testout/filterz.dat')
    write(18,'(3(1X,A15))')'y','q','q~'
    write(18,'(3(1X,E15.7E3))')(x(0,0,k,3),q(0,0,k,1),dq(k),k=0,km)
    close(18)
    print*,' << testout/filterz.dat'
    !
    deallocate(dq)
    !
  end subroutine filtertest
  !
    !
end module test
!+---------------------------------------------------------------------+
!| The end of the module test.                                         |
!+---------------------------------------------------------------------+