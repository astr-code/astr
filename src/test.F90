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
  use parallel, only: mpistop
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
    elseif(testmode=='accu') then 
      call accuracytest
    elseif(testmode=='enst') then 
      call enstest
    elseif(testmode=='filt') then 
      call filtertest
    elseif(testmode=='flux') then 
      call fluxtest
    elseif(testmode=='bc') then 
      call testbc
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
  subroutine testbc
    !
    use bc,        only : inflowintx
    use parallel,  only : qswap,psum,irk,jrk,krk,mpirank,mpirankname
    use tecio
    !
    if(irk==0) then
      !
      call inflowintx
      !
    endif
    !
  end subroutine testbc
  !
  subroutine enstest

    use commvar,   only : im,jm,km,ia,ja,ka,roinf,uinf,hm,is,ie,js,je,ks,ke
    use commarray, only : x,vel,dvel,rho
    use comsolver, only : gradcal
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
  subroutine accuracytest
    !
    use commvar,   only : im,jm,km,npdci,npdcj,npdck,conschm,          &
                          alfa_filter,numq,is,ie,ia
    use commarray, only : x,q,dxi
    use bc,        only : boucon
    use parallel,  only : dataswap,mpirankname,psum,pmax
    use comsolver, only : alfa_con,cci,ccj,cck
    use derivative, only : fds,fds_compact_i
    !
    ! local data
    integer :: i,j,k,n
    real(8) :: dx,error1,error2,errorinf
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
      q(i,j,k,1)=sin(4.d0*x(i,j,k,1))
      !
    enddo
    enddo
    enddo
    !
    allocate(dqref(0:im))
    do i=0,im
      dqref(i)=4.d0*cos(4.d0*x(i,0,0,1))
    enddo
    !
    call dataswap(q)
    !
    allocate(dq(0:im))
    !
    dq(:)=fds%central(fds_compact_i,f=q(:,0,0,1),dim=im)*dxi(:,0,0,1,1)
    !
    error1=0.d0
    error2=0.d0
    errorinf=0.d0
    do i=1,im
      error1=error1+abs(dq(i)-dqref(i))
      error2=error2+(dq(i)-dqref(i))**2
      errorinf=max(errorinf,abs(dq(i)-dqref(i)))
    enddo
    !
    error1=psum(error1)/dble(ia)
    error2=sqrt(psum(error1)/dble(ia))
    errorinf=pmax(errorinf)
    !
    print*,' ** total number of nodes:',ia
    print*,' **   L1 error:',error1
    print*,' **   L2 error:',error2
    print*,' ** Linf error:',errorinf
    !
  end subroutine accuracytest
  !
  subroutine filtertest2

    use commvar,   only : im,jm,km
    use commarray, only : x,q
    use comsolver, only: filterq
    use parallel,  only : dataswap,mpirankname,jrkm,jrk

    integer :: i,j,k,n

    do k=0,km
    do j=0,jm
    do i=0,im
      !
      q(i,j,k,1)=sin(10.d0*x(i,j,k,1))
      !
    enddo
    enddo
    enddo

    call filterq

    j=jm/2
    k=jm/2
    open(18,file='testout/filterqx'//mpirankname//'.dat')
    write(18,'(3(1X,A15))')'x','q','q~'
    write(18,'(3(1X,E15.7E3))')(x(i,j,k,1),sin(10.d0*x(i,j,k,1)),q(i,j,k,1),i=0,im)
    close(18)
    print*,' << testout/filterqx.dat'
    !
    ! open(18,file='testout/filterqy'//mpirankname//'.dat')
    ! write(18,'(2(1X,A15))')'y','q','q~'
    ! write(18,'(2(1X,E15.7E3))')(x(im/2,j,km/2,2),q(im/2,j,km/2,1),j=0,jm)
    ! close(18)
    ! print*,' << testout/filterqy.dat'
    ! !
    ! open(18,file='testout/filterqz'//mpirankname//'.dat')
    ! write(18,'(2(1X,A15))')'y','q','q~'
    ! write(18,'(2(1X,E15.7E3))')(x(im/2,jm/2,k,3),q(im/2,jm/2,k,1),k=0,km)
    ! close(18)
    ! print*,' << testout/filterqz.dat'

  end subroutine filtertest2

  subroutine filtertest
    !
    use commvar,   only : im,jm,km,npdci,npdcj,npdck,conschm,          &
                          alfa_filter,numq,is,ie,ia,ja,ka
    use commarray, only : x,q,dxi
    use bc,        only : boucon
    use parallel,  only : dataswap,mpirankname,jrkm,jrk,mpirank,psum,ptime
    use comsolver, only : alfa_con,fci,fcj,fck
    use filter,    only : compact_filter,filter_i,filter_j,filter_k,filter_ii,filter_jj,filter_kk
    !
    ! local data
    integer :: i,j,k,n
    real(8) :: dx,error,time_beg,subtime
    real(8),allocatable :: vtest(:,:,:)
    real(8),allocatable :: dq(:),qhp(:),qhm(:),dqref(:)
    !
    ! print*,x(:,0,0,1)
    !
    time_beg=ptime() 

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
    ! dq(:)=spafilter10(q(:,0,0,1),npdci,im,alfa_filter,fci)
    dq(:)=compact_filter(afilter=filter_i,f=q(:,0,0,1),dim=im)

    error=0.d0
    do i=1,im
      error=error+(dq(i)-q(i,0,0,1))**2
    enddo

    error=psum(error)/(ia)

    if(mpirank==0) then
      print*,' ** filtering a smooth function in x direction,         error:',error
    endif
    ! open(18,file='testout/filterx'//mpirankname//'.dat')
    ! write(18,'(3(1X,A15))')'x','q','q~'
    ! write(18,'(3(1X,E15.7E3))')(x(i,0,0,1),q(i,0,0,1),dq(i),i=0,im)
    ! close(18)
    ! print*,' << testout/filterx.dat'
    ! !
    deallocate(dq)

    ! testing ddy
    do k=0,km
    do j=0,jm
      !
      i=0
      q(i,j,k,1)=sin(0.5d0*pi*x(i,j,k,2))+0.1d0*sin((dble(j)+0.5d0)*pi)
      !
      i=1
      q(i,j,k,1)=sin(0.5d0*pi*x(i,j,k,2))
    enddo
    enddo
    if(jrk==0) q(:, 0,:,1)=0.d0
    if(jrk==jrkm) q(:,jm,:,1)=0.d0
    !
    call dataswap(q)
    !
    allocate(dq(0:jm))
    !
    ! dq(:)=compact_filter(afilter=filter_jj,f=q(0,:,0,1),ntype=npdcj,dim=jm,note='boundary_no_filter')
    dq(:)=compact_filter(afilter=filter_j,f=q(0,:,0,1),dim=jm)
    
    error=0.d0
    do j=1,jm
      error=error+(dq(j)-q(1,j,0,1))**2
    enddo

    error=psum(error)/(ja)

    if(mpirank==0) then
      print*,' ** filtering a smooth+wiggled function in y direction, error:',error
    endif

    ! open(18,file='testout/filtery'//mpirankname//'.dat')
    ! write(18,'(3(1X,A15))')'y','q','q~'
    ! write(18,'(3(1X,E15.7E3))')(x(0,j,0,2),q(0,j,0,1),dq(j),j=0,jm)
    ! close(18)
    ! print*,' << testout/filtery.dat'
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
    dq(:)=compact_filter(afilter=filter_k,f=q(0,0,:,1),dim=km)
    ! dq(:)=compact_filter(afilter=filter_kk,f=q(0,0,:,1),ntype=npdck,dim=km,note='boundary_no_filter')
    !
    error=0.d0
    do k=1,km
      error=error+(dq(k)-q(0,0,k,1))**2
    enddo

    error=psum(error)/(ka)

    if(mpirank==0) then
      print*,' ** filtering a smooth function in z direction,         error:',error
    endif

    subtime=subtime+ptime()-time_beg
    !
    if(mpirank==0) print*,' ** time cost:',subtime
    ! open(18,file='testout/filterz'//mpirankname//'.dat')
    ! write(18,'(3(1X,A15))')'y','q','q~'
    ! write(18,'(3(1X,E15.7E3))')(x(0,0,k,3),q(0,0,k,1),dq(k),k=0,km)
    ! close(18)
    ! print*,' << testout/filterz.dat'
    !
    deallocate(dq)
    !
  end subroutine filtertest
  
  subroutine gradtest
    !
    use commvar,   only : im,jm,km,ia,ja,ka,npdci,npdcj,npdck,conschm,          &
                          alfa_filter,numq,is,ie,hm,difschm
    use commarray, only : x,q,dxi
    use comsolver, only : alfa_con,alfa_dif,cci,ccj,cck,dci,dcj,dck
    use bc,        only : boucon
    use parallel,  only : dataswap,mpirank,mpirankname,ptime,psum,mpistop
    use derivative, only : fds,fds_compact_i,fds_compact_j,fds_compact_k
    use tecio
    !
    ! local data
    integer :: i,j,k,n,s,asize,ncolm,counter
    real(8) :: dx,var1
    real(8),allocatable :: vtest(:,:,:,:),dvtes(:,:,:,:,:),dvref(:,:,:,:,:)
    real(8),allocatable :: f1(:,:),df1(:,:)
    real(8),allocatable :: ff(:,:,:,:),dff(:,:,:,:)
    real(8),allocatable :: error(:,:)
    !
    real(8) :: time_beg
    real(8),save :: subtime=0.d0
    !
    ! call x3d2init(im,jm,km)
    ! print*,x(:,0,0,1)
    !
    ! testing ddx
    ncolm=1
    allocate(vtest(-hm:im+hm,-hm:jm+hm,-hm:km+hm,1:ncolm))
    allocate(dvtes(0:im,0:jm,0:km,1:ncolm,1:3),dvref(0:im,0:jm,0:km,1:ncolm,1:3))
    allocate(error(ncolm,3))

    dvtes=0.d0
    !
    time_beg=ptime() 
    counter=0
    !
    do k=0,km
    do j=0,jm
    do i=0,im
      !
      ! do n=1,ncolm
      !   call random_number(var1)
      !   vtest(i,j,k,n)=var1
      ! enddo
      ! vtest(i,j,k,1)=sin(x(i,j,k,1))*sin(0.5d0*pi*x(i,j,k,2))
      vtest(i,j,k,1)=sin(x(i,j,k,1))*sin(0.5d0*pi*x(i,j,k,2))*cos(2.d0*x(i,j,k,3))
      ! !
      dvref(i,j,k,1,1)= cos(x(i,j,k,1))*sin(0.5d0*pi*x(i,j,k,2))*cos(2.d0*x(i,j,k,3))
      dvref(i,j,k,1,2)= 0.5d0*pi*sin(x(i,j,k,1))*cos(0.5d0*pi*x(i,j,k,2))*cos(2.d0*x(i,j,k,3))
      dvref(i,j,k,1,3)=-2.d0*sin(x(i,j,k,1))*sin(0.5d0*pi*x(i,j,k,2))*sin(2.d0*x(i,j,k,3))

    enddo
    enddo
    enddo

    call dataswap(vtest)
    !
    allocate(dff(0:im,0:jm,0:km,ncolm))

    allocate(ff(-hm:im+hm,0:jm,0:km,ncolm))

    ff(:,0:jm,0:km,1:ncolm)=vtest(:,0:jm,0:km,1:ncolm)
    
    do n=1,1
    do k=0,km
    do j=0,jm
        dff(:,j,k,n)=fds%central(fds_compact_i,f=ff(:,j,k,n),dim=im)
    enddo
    enddo
    enddo


    do k=0,km
    do j=0,jm
    do i=0,im
      do n=1,ncolm
        dvtes(i,j,k,n,1)=dff(i,j,k,n)*dxi(i,j,k,1,1)
        dvtes(i,j,k,n,2)=dff(i,j,k,n)*dxi(i,j,k,1,2)
        dvtes(i,j,k,n,3)=dff(i,j,k,n)*dxi(i,j,k,1,3)
      enddo
    enddo
    enddo
    enddo

    deallocate(ff)

    allocate(ff(0:im,-hm:jm+hm,0:km,ncolm))

    ff(0:im,:,0:km,1:ncolm)=vtest(0:im,:,0:km,1:ncolm)
    
    do n=1,1
    do k=0,km
    do i=0,im
        dff(i,:,k,n)=fds%central(fds_compact_j,f=ff(i,:,k,n),dim=jm)
    enddo
    enddo
    enddo

    do k=0,km
    do j=0,jm
    do i=0,im
      do n=1,ncolm
        dvtes(i,j,k,n,1)=dvtes(i,j,k,n,1)+dff(i,j,k,n)*dxi(i,j,k,2,1)
        dvtes(i,j,k,n,2)=dvtes(i,j,k,n,2)+dff(i,j,k,n)*dxi(i,j,k,2,2)
        dvtes(i,j,k,n,3)=dvtes(i,j,k,n,3)+dff(i,j,k,n)*dxi(i,j,k,2,3)
      enddo
    enddo
    enddo
    enddo

    deallocate(ff)

    ! i=im
    ! k=km/4
    ! open(18,file='testout/test'//mpirankname//'.dat')
    ! write(18,'(3(1X,A15))')'x','f','df'
    ! write(18,'(3(1X,E15.7E3))')(x(i,j,k,2),ff(i,j,k,1),dff(i,j,k,1),j=0,jm)
    ! close(18)
    ! print*,' << testout/test'//mpirankname//'.dat'

    allocate(ff(0:im,0:jm,-hm:km+hm,ncolm))

    ff(0:im,0:jm,:,1:ncolm)=vtest(0:im,0:jm,:,1:ncolm)

    do n=1,1
    do j=0,jm
    do i=0,im
        ! dff(i,j,:,n)=ddfc(ff(i,j,:,n),conschm,npdck,km,alfa_con,cck)
        dff(i,j,:,n)=fds%central(fds_compact_k,f=ff(i,j,:,n),dim=km)
    enddo
    enddo
    enddo

    do k=0,km
    do j=0,jm
    do i=0,im
      do n=1,ncolm
        dvtes(i,j,k,n,1)=dvtes(i,j,k,n,1)+dff(i,j,k,n)*dxi(i,j,k,3,1)
        dvtes(i,j,k,n,2)=dvtes(i,j,k,n,2)+dff(i,j,k,n)*dxi(i,j,k,3,2)
        dvtes(i,j,k,n,3)=dvtes(i,j,k,n,3)+dff(i,j,k,n)*dxi(i,j,k,3,3)
      enddo
    enddo
    enddo
    enddo
    deallocate(ff)

    error=0.d0
    do k=1,km
    do j=1,jm
    do i=1,im
      do n=1,ncolm
        error(n,1)=error(n,1)+(dvtes(i,j,k,n,1)-dvref(i,j,k,n,1))**2
        error(n,2)=error(n,2)+(dvtes(i,j,k,n,2)-dvref(i,j,k,n,2))**2
        error(n,3)=error(n,3)+(dvtes(i,j,k,n,3)-dvref(i,j,k,n,3))**2
      enddo
    enddo
    enddo
    enddo

    error=psum(error)
    ! error=sqrt(error)

    if(mpirank==0) then
      error=sqrt(error/(ia*ja*ka))

      print*,' ** error in x gradient',error(:,1)
      print*,' ** error in y gradient',error(:,2)
      print*,' ** error in z gradient',error(:,3)
    endif

    subtime=subtime+ptime()-time_beg
    !
    if(mpirank==0) print*,' ** time cost:',subtime


    ! call tecbin('testout/tectest'//mpirankname//'.plt',            &
    !                                   x(0:im,0:jm,0:km,1),'x',     &
    !                                   x(0:im,0:jm,0:km,2),'y',     &
    !                                   x(0:im,0:jm,0:km,3),'z',     &
    !                               vtest(0:im,0:jm,0:km,1),'v1',    &
    !                               dvtes(0:im,0:jm,0:km,1,1),'dv1dy',    &
    !                               dvref(0:im,0:jm,0:km,1,2),'dv1dy_ref' )
  end subroutine gradtest

  subroutine fluxtest
    !
    use commvar,   only : im,jm,km,ia,ja,ka,npdci,npdcj,npdck,conschm,          &
                          alfa_filter,numq,is,ie,hm,difschm,is,ie,js,je,ks,ke
    use commarray, only : x,q,dxi
    use comsolver, only : alfa_con,alfa_dif,cci,ccj,cck,dci,dcj,dck
    use bc,        only : boucon
    use parallel,  only : dataswap,mpirank,mpirankname,ptime,psum,mpistop
    use derivative,only : fds
    use flux,      only : flux_compact,flux_uw_i,flux_dw_i,flux_uw_j,flux_dw_j, &
                          flux_uw_k,flux_dw_k

    use tecio
    !
    ! local data
    integer :: i,j,k,n,s,asize,ncolm,counter
    real(8) :: dx,var1
    real(8),allocatable :: vtest(:,:,:,:),dvtes(:,:,:,:,:),dvref(:,:,:,:,:)
    real(8),allocatable :: f1(:,:),df1(:,:)
    real(8),allocatable :: ff(:,:,:,:),dff(:,:,:,:),fh(:,:,:,:)
    real(8),allocatable :: error(:,:)
    !
    real(8) :: time_beg
    real(8),save :: subtime=0.d0
    !
    ! call x3d2init(im,jm,km)
    ! print*,x(:,0,0,1)
    !
    ! testing ddx
    ncolm=1
    allocate(vtest(-hm:im+hm,-hm:jm+hm,-hm:km+hm,1:ncolm))
    allocate(dvtes(0:im,0:jm,0:km,1:ncolm,1:3),dvref(0:im,0:jm,0:km,1:ncolm,1:3))
    allocate(error(ncolm,3))

    dvtes=0.d0
    !
    time_beg=ptime() 
    counter=0
    !
    do k=0,km
    do j=0,jm
    do i=0,im
      !
      ! do n=1,ncolm
      !   call random_number(var1)
      !   vtest(i,j,k,n)=var1
      ! enddo
      ! vtest(i,j,k,1)=sin(x(i,j,k,1))*sin(0.5d0*pi*x(i,j,k,2))
      vtest(i,j,k,1)=sin(x(i,j,k,1))*sin(0.5d0*pi*x(i,j,k,2))*cos(2.d0*x(i,j,k,3))
      ! !
      dvref(i,j,k,1,1)= cos(x(i,j,k,1))*sin(0.5d0*pi*x(i,j,k,2))*cos(2.d0*x(i,j,k,3))
      dvref(i,j,k,1,2)= 0.5d0*pi*sin(x(i,j,k,1))*cos(0.5d0*pi*x(i,j,k,2))*cos(2.d0*x(i,j,k,3))
      dvref(i,j,k,1,3)=-2.d0*sin(x(i,j,k,1))*sin(0.5d0*pi*x(i,j,k,2))*sin(2.d0*x(i,j,k,3))

    enddo
    enddo
    enddo

    call dataswap(vtest)
    !
    allocate(dff(0:im,0:jm,0:km,ncolm),fh(-1:im,-1:jm,-1:km,ncolm))

    allocate(ff(-hm:im+hm,0:jm,0:km,ncolm))

    ff(:,0:jm,0:km,1:ncolm)=vtest(:,0:jm,0:km,1:ncolm)
    
    do n=1,1
    do k=0,km
    do j=0,jm
        fh(:,j,k,n)=flux_compact(asolver=flux_uw_i,f=ff(:,j,k,n),dim=im)
        do i=0,im
            dff(i,j,k,1)=fh(i,j,k,1)-fh(i-1,j,k,1)
        enddo
    enddo
    enddo
    enddo

    ! j=jm/2
    ! k=km/4
    ! open(18,file='testout/fh_i'//mpirankname//'.dat')
    ! write(18,'(5(1X,A15))')'x','f','df','xh','fh'
    ! do i=0,im
    !   write(18,'(5(1X,E15.7E3))')x(i,j,k,1),ff(i,j,k,1),0.5d0*(x(i,j,k,1)+x(i+1,j,k,1)),fh(i,j,k,1)
    ! enddo
    ! close(18)
    ! print*,' << testout/fh_i'//mpirankname//'.dat'

    ! open(18,file='testout/dfh_i'//mpirankname//'.dat')
    ! write(18,'(3(1X,A15))')'x','df_ref','df'
    ! do i=0,im
    !   write(18,'(3(1X,E15.7E3))')x(i,j,k,1),dvref(i,j,k,1,1),dff(i,j,k,1)*dxi(i,j,k,1,1)
    ! enddo
    ! close(18)
    ! print*,' << testout/dfh_i'//mpirankname//'.dat'

    do k=0,km
    do j=0,jm
    do i=0,im
      do n=1,ncolm
        dvtes(i,j,k,n,1)=dff(i,j,k,n)*dxi(i,j,k,1,1)
        dvtes(i,j,k,n,2)=dff(i,j,k,n)*dxi(i,j,k,1,2)
        dvtes(i,j,k,n,3)=dff(i,j,k,n)*dxi(i,j,k,1,3)
      enddo
    enddo
    enddo
    enddo

    deallocate(ff)

    allocate(ff(0:im,-hm:jm+hm,0:km,ncolm))

    ff(0:im,:,0:km,1:ncolm)=vtest(0:im,:,0:km,1:ncolm)
    

    do n=1,1
    do k=0,km
    do i=0,im
        fh(i,:,k,n)=flux_compact(asolver=flux_uw_j,f=ff(i,:,k,n),dim=jm)
        do j=0,jm
            dff(i,j,k,1)=fh(i,j,k,1)-fh(i,j-1,k,1)
        enddo
    enddo
    enddo
    enddo

    i=im/4
    k=0
    ! open(18,file='testout/fh_j'//mpirankname//'.dat')
    ! write(18,'(4(1X,A15))')'x','f','xh','fh'
    ! do j=0,jm
    !   write(18,'(4(1X,E15.7E3))')x(i,j,k,2),ff(i,j,k,1),0.5d0*(x(i,j,k,2)+x(i,j+1,k,2)),fh(i,j,k,1)
    ! enddo
    ! close(18)
    ! print*,' << testout/fh_j'//mpirankname//'.dat'
    open(18,file='testout/dfh_j'//mpirankname//'.dat')
    write(18,'(3(1X,A15))')'x','df_ref','df'
    do j=js,je
      write(18,'(3(1X,E15.7E3))')x(i,j,k,2),dvref(i,j,k,1,2),dff(i,j,k,1)*dxi(i,j,k,2,2)
    enddo
    close(18)
    print*,' << testout/dfh_j'//mpirankname//'.dat'

    do k=0,km
    do j=0,jm
    do i=0,im
      do n=1,ncolm
        dvtes(i,j,k,n,1)=dvtes(i,j,k,n,1)+dff(i,j,k,n)*dxi(i,j,k,2,1)
        dvtes(i,j,k,n,2)=dvtes(i,j,k,n,2)+dff(i,j,k,n)*dxi(i,j,k,2,2)
        dvtes(i,j,k,n,3)=dvtes(i,j,k,n,3)+dff(i,j,k,n)*dxi(i,j,k,2,3)
      enddo
    enddo
    enddo
    enddo


    deallocate(ff)


    allocate(ff(0:im,0:jm,-hm:km+hm,ncolm))

    ff(0:im,0:jm,:,1:ncolm)=vtest(0:im,0:jm,:,1:ncolm)

    do n=1,1
    do j=0,jm
    do i=0,im
        ! dff(i,j,:,n)=ddfc(ff(i,j,:,n),conschm,npdck,km,alfa_con,cck)
        fh(i,j,:,n)=flux_compact(asolver=flux_uw_k,f=ff(i,j,:,n),dim=km)
        do k=0,km
            dff(i,j,k,1)=fh(i,j,k,1)-fh(i,j,k-1,1)
        enddo
    enddo
    enddo
    enddo

    ! i=im/2
    ! j=jm/2
    ! open(18,file='testout/fh_k'//mpirankname//'.dat')
    ! write(18,'(4(1X,A15))')'x','f','xh','fh'
    ! do k=0,km
    !   write(18,'(4(1X,E15.7E3))')x(i,j,k,3),ff(i,j,k,1),0.5d0*(x(i,j,k,3)+x(i,j,k+1,3)),fh(i,j,k,1)
    ! enddo
    ! close(18)
    ! print*,' << testout/fh_k'//mpirankname//'.dat'

    ! open(18,file='testout/dfh_k'//mpirankname//'.dat')
    ! write(18,'(3(1X,A15))')'x','df_ref','df'
    ! do k=0,km
    !   write(18,'(3(1X,E15.7E3))')x(i,j,k,3),dvref(i,j,k,1,3),dff(i,j,k,1)*dxi(i,j,k,3,3)
    ! enddo
    ! close(18)
    ! print*,' << testout/dfh_k'//mpirankname//'.dat'

    do k=0,km
    do j=0,jm
    do i=0,im
      do n=1,ncolm
        dvtes(i,j,k,n,1)=dvtes(i,j,k,n,1)+dff(i,j,k,n)*dxi(i,j,k,3,1)
        dvtes(i,j,k,n,2)=dvtes(i,j,k,n,2)+dff(i,j,k,n)*dxi(i,j,k,3,2)
        dvtes(i,j,k,n,3)=dvtes(i,j,k,n,3)+dff(i,j,k,n)*dxi(i,j,k,3,3)
      enddo
    enddo
    enddo
    enddo

    deallocate(ff)

    error=0.d0
    do k=ks,ke
    do j=js,je
    do i=is,ie
      do n=1,ncolm
        error(n,1)=error(n,1)+(dvtes(i,j,k,n,1)-dvref(i,j,k,n,1))**2
        error(n,2)=error(n,2)+(dvtes(i,j,k,n,2)-dvref(i,j,k,n,2))**2
        error(n,3)=error(n,3)+(dvtes(i,j,k,n,3)-dvref(i,j,k,n,3))**2
      enddo
    enddo
    enddo
    enddo

    error=psum(error)
    ! error=sqrt(error)

    if(mpirank==0) then
      error=sqrt(error/(ia*ja*ka))

      print*,' ** error in x gradient',error(:,1)
      print*,' ** error in y gradient',error(:,2)
      print*,' ** error in z gradient',error(:,3)
    endif

    subtime=subtime+ptime()-time_beg
    !
    if(mpirank==0) print*,' ** time cost:',subtime

    ! call tecbin('testout/tectest'//mpirankname//'.plt',            &
    !                                   x(0:im,0:jm,0:km,1),'x',     &
    !                                   x(0:im,0:jm,0:km,2),'y',     &
    !                                   x(0:im,0:jm,0:km,3),'z',     &
    !                               vtest(0:im,0:jm,0:km,1),'v1',    &
    !                               dvtes(0:im,0:jm,0:km,1,2),'dv1dy',    &
    !                               dvref(0:im,0:jm,0:km,1,2),'dv1dy_ref' )
  end subroutine fluxtest

end module test
!+---------------------------------------------------------------------+
!| The end of the module test.                                         |
!+---------------------------------------------------------------------+