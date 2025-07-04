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
  use parallel, only: mpistop,mpirank
  !
  implicit none
  !
  contains
  !!+-------------------------------------------------------------------+
  !| This subroutine is the entrance of the test                        |
  !+-------------------------------------------------------------------+
  subroutine codetest
  
    use commvar,   only : testmode
    use parallel,  only : mpistop, bcast
    use cmdefne,   only : readkeyboad

    !-------------------------------------------------------------------
    ! Read test mode and broadcast it
    !-------------------------------------------------------------------
    call readkeyboad(testmode)

    call bcast(testmode)
  
    !-------------------------------------------------------------------
    ! Execute selected test
    !-------------------------------------------------------------------
    select case (trim(testmode))
    
    case ('grad')
      call gradtest
  
    case ('flux')
      call fluxtest
  
    case ('filt')
      call filtertest

    case ('accu')
      call accuracytest
  
    case ('enst')
      call enstest

    case ('bc')
      call testbc
  
    case default

      if(mpirank==0) then
        write(*,*) ' +------------------------------------------------------------+'
        write(*,*) ' | Available test options                                     |'
        write(*,*) ' +------------------------------------------------------------+'
        write(*,*) ' | grad    - Test gradient calculation                        |'
        write(*,*) ' | flux    - Test flux schemes                                |'
        write(*,*) ' | filt    - Test filter routines                             |'
        write(*,*) ' | accu    - Test numerical accuracy                          |'
        write(*,*) ' | enst    - Test enstrophy evaluation                        |'
        write(*,*) ' | bc      - Test boundary condition treatment                |'
        write(*,*) ' +------------------------------------------------------------+'
      endif
    
    end select
  
    call mpistop
  
  end subroutine codetest
  !+-------------------------------------------------------------------+
  !| End of subroutine codetest                                        |
  !+-------------------------------------------------------------------+

  subroutine test_init

  end subroutine test_init

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

  !+-------------------------------------------------------------------+
  !| Subroutine: filtertest                                            |
  !|                                                                   |
  !| Purpose:                                                          |
  !|   This subroutine verifies the performance of 1D compact filters  |
  !|   in the x, y, and z directions using known smooth or perturbed   |
  !|   test functions. It calculates the L2 error between filtered and |
  !|   original fields to assess filter accuracy and boundary behavior.|
  !|                                                                   |
  !| Functionality:                                                    |
  !|   - Initializes a grid and analytical test fields                 |
  !|   - Applies compact filters in x/y/z directions                   |
  !|   - Computes L2 error and reports performance                     |
  !|                                                                   |
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 28-05-2021  | Created by J. Fang @ Warrington                     |
  !| 03-07-2025  | Remastered by J. Fang @ Imech, Beijing, with ChatGPT|
  !+-------------------------------------------------------------------+
  subroutine filtertest
  
    use commvar,   only : im, jm, km, npdci, npdcj, npdck, ia, ja, ka, hm, &
                          alfa_filter, lihomo, ljhomo, lkhomo
    use commarray, only : x
    use parallel,  only : mpisizedis, parapp, parallelini, dataswap,      &
                          mpirankname, mpirank, psum, ptime, jrk, jrkm
    use solver,    only : refcal
    use gridgeneration, only : gridcube
    use filter,    only : compact_filter_initiate, filter_coefficient_cal, &
                          filter_i, filter_j, filter_k, compact_filter
  
    implicit none
  
    ! Local variables
    integer :: i, j, k, n
    real(8) :: dx, error, time_beg, subtime
    real(8), allocatable :: fv(:,:,:)
    real(8), allocatable :: ffv(:)
  
    !---------------------------------------------------------------
    ! Initialization
    !---------------------------------------------------------------
    ia = 128
    ja = 128
    ka = 128
  
    lihomo = .true.
    ljhomo = .false.
    lkhomo = .true.
  
    alfa_filter = 0.49d0
  
    call mpisizedis
    call parapp
    call parallelini
    call refcal
  
    call filter_coefficient_cal(alfa=alfa_filter, beter_halo=1.11d0, beter_bouond=0.98d0)
    call compact_filter_initiate(afilter=filter_i, ntype=npdci, dim=im)
    call compact_filter_initiate(afilter=filter_j, ntype=npdcj, dim=jm)
    call compact_filter_initiate(afilter=filter_k, ntype=npdck, dim=km)
  
    allocate(x(-hm:im+hm, -hm:jm+hm, -hm:km+hm, 1:3))
    call gridcube(2.d0*pi, 2.d0, 2.d0*pi)
    allocate(fv(-hm:im+hm, -hm:jm+hm, -hm:km+hm))
  
    !---------------------------------------------------------------
    ! Test in x-direction
    !---------------------------------------------------------------
    do k = 0, km
      do j = 0, jm
        do i = 0, im
          fv(i,j,k) = sin(10.d0 * x(i,j,k,1))
        end do
      end do
    end do
  
    time_beg = ptime()
  
    call dataswap(fv)
  
    allocate(ffv(0:im))
    ffv(:) = compact_filter(afilter=filter_i, f=fv(:,0,0), dim=im)
  
    error = 0.d0
    do i = 1, im
      error = error + (ffv(i) - fv(i,0,0))**2
    end do
    error = psum(error) / ia
  
    if (mpirank == 0) then
      write(*,*) ' +------------------------------------------------------------+'
      write(*,*) ' | errors caused by compact_filter                            |'
      write(*,*) ' +-------------------------------+----------------------------+'
      write(*,*) ' | x (homo+smooth)               |', error,' |'
    end if
  
    deallocate(ffv)
  
    !---------------------------------------------------------------
    ! Test in y-direction with wiggle
    !---------------------------------------------------------------
    do k = 0, km
      do j = 0, jm
        fv(0,j,k) = sin(0.5d0 * pi * x(0,j,k,2)) + 0.1d0 * sin((dble(j)+0.5d0)*pi)
        fv(1,j,k) = sin(0.5d0 * pi * x(1,j,k,2))
      end do
    end do
  
    if (jrk == 0)    fv(:, 0, :) = 0.d0
    if (jrk == jrkm) fv(:, jm, :) = 0.d0
  
    call dataswap(fv)
  
    allocate(ffv(0:jm))
    ffv(:) = compact_filter(afilter=filter_j, f=fv(0,:,0), dim=jm)
  
    error = 0.d0
    do j = 1, jm
      error = error + (ffv(j) - fv(1,j,0))**2
    end do
    error = psum(error) / ja
  
    if (mpirank == 0) then
      write(*,*) ' | y (non-homo + wiggles)        |', error,' |'
    end if
  
    deallocate(ffv)
  
    !---------------------------------------------------------------
    ! Test in z-direction
    !---------------------------------------------------------------
    do k = 0, km
      do j = 0, jm
        do i = 0, im
          fv(i,j,k) = sin(10.d0 * x(i,j,k,3))
        end do
      end do
    end do
  
    call dataswap(fv)
  
    allocate(ffv(0:km))
    ffv(:) = compact_filter(afilter=filter_k, f=fv(0,0,:), dim=km)
  
    error = 0.d0
    do k = 1, km
      error = error + (ffv(k) - fv(0,0,k))**2
    end do
    error = psum(error) / ka
  
    if (mpirank == 0) then
      write(*,*) ' | z (homo+smooth)               |', error,' |'
      write(*,'(A,F11.6,A)') '  | Time cost                     |', ptime() - time_beg,' s               |'
      write(*,*) ' +-------------------------------+----------------------------+'
    end if
  
    deallocate(ffv)
    deallocate(fv)
  
  end subroutine filtertest
  !+-------------------------------------------------------------------+
  !| End of subroutine filtertest                                      |
  !+-------------------------------------------------------------------+

  !+-------------------------------------------------------------------+
  !| Subroutine: gradtest                                              |
  !|                                                                   |
  !| Purpose:                                                          |
  !|   This subroutine tests the accuracy of spatial gradient          |
  !|   computations using high-order compact finite difference schemes.|
  !|   It computes the numerical derivatives of a known trigonometric  |
  !|   field and compares them to analytical gradients in x, y, z.     |
  !|   The L2 error of the derivatives is reported for validation.     |
  !|                                                                   |
  !| Functionality:                                                    |
  !|   - Initializes a 3D periodic grid                                |
  !|   - Constructs an analytic test field                             |
  !|   - Computes reference (exact) derivatives                        |
  !|   - Applies finite difference schemes in x/y/z directions         |
  !|   - Computes and outputs L2-norm errors of numerical gradients    |
  !|                                                                   |
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 28-05-2021  | Created by J. Fang @ Warrington                     |
  !| 03-07-2025  | Remastered by J. Fang @ Imech, CAS, Beijing         |
  !+-------------------------------------------------------------------+
  subroutine gradtest
  
    use commvar,   only : im,jm,km,hm,ia,ja,ka,conschm,difschm,        &
                          lihomo,ljhomo,lkhomo
    use commarray, only : x
    use comsolver, only : solvrinit
    use parallel,  only : mpiinitial,mpistop,mpisizedis,lio,           &
                          parallelini,parapp,dataswap,mpirank,         &
                          mpirankname,ptime,psum
    use solver,    only : refcal
    use derivative,only : fds,fds_compact_i,fds_compact_j,fds_compact_k
    use gridgeneration, only : gridcube
    use tecio
  
    implicit none
  
    !-------------------------------------------------------------------
    ! Local variables
    !-------------------------------------------------------------------
    integer :: i, j, k, n, s, asize, ncolm, counter
    real(8) :: dx, dy, dz, var1
    real(8), allocatable :: vtest(:,:,:,:), dvtes(:,:,:,:,:), dvref(:,:,:,:,:)
    real(8), allocatable :: f1(:,:), df1(:,:)
    real(8), allocatable :: ff(:,:,:,:), dff(:,:,:,:)
    real(8), allocatable :: error(:,:)
    real(8) :: time_beg
    real(8), save :: subtime = 0.d0
  
    !-------------------------------------------------------------------
    ! Initialization
    !-------------------------------------------------------------------
    difschm  = '643c'
    conschm  = '643c'
  
    ia = 128
    ja = 128
    ka = 128
  
    lihomo = .true.
    ljhomo = .false.
    lkhomo = .true.
  
    call mpisizedis
    call parapp
    call parallelini
    call refcal
    call solvrinit
  
    allocate(x(-hm:im+hm, -hm:jm+hm, -hm:km+hm, 1:3))
    call gridcube(2.d0*pi, 2.d0, 2.d0*pi)
  
    dx = x(1,0,0,1) - x(0,0,0,1)
    dy = x(0,1,0,2) - x(0,0,0,2)
    dz = x(0,0,1,3) - x(0,0,0,3)
  
    ncolm = 1
    allocate(vtest(-hm:im+hm, -hm:jm+hm, -hm:km+hm, ncolm))
    allocate(dvtes(0:im,0:jm,0:km,ncolm,3), dvref(0:im,0:jm,0:km,ncolm,3))
    allocate(error(ncolm, 3))
    dvtes = 0.d0
  
    !-------------------------------------------------------------------
    ! Construct test field and reference derivatives
    !-------------------------------------------------------------------
    time_beg = ptime()
    counter = 0
  
    do k = 0, km
      do j = 0, jm
        do i = 0, im
          vtest(i,j,k,1)   = sin(x(i,j,k,1)) * sin(0.5d0*pi*x(i,j,k,2)) * sin(x(i,j,k,3))
          dvref(i,j,k,1,1) = cos(x(i,j,k,1)) * sin(0.5d0*pi*x(i,j,k,2)) * sin(x(i,j,k,3))
          dvref(i,j,k,1,2) = 0.5d0*pi*sin(x(i,j,k,1)) * cos(0.5d0*pi*x(i,j,k,2)) * sin(x(i,j,k,3))
          dvref(i,j,k,1,3) = sin(x(i,j,k,1)) * sin(0.5d0*pi*x(i,j,k,2)) * cos(x(i,j,k,3))
        end do
      end do
    end do
  
    call dataswap(vtest)
  
    !-------------------------------------------------------------------
    ! x-direction gradient
    !-------------------------------------------------------------------
    allocate(dff(0:im,0:jm,0:km,ncolm))
    allocate(ff(-hm:im+hm,0:jm,0:km,ncolm))
    ff(:,0:jm,0:km,:) = vtest(:,0:jm,0:km,:)
  
    do n = 1, ncolm
      do k = 0, km
        do j = 0, jm
          dff(:,j,k,n) = fds%central(fds_compact_i, f=ff(:,j,k,n), dim=im)
        end do
      end do
    end do
  
    do k = 0, km
      do j = 0, jm
        do i = 0, im
          do n = 1, ncolm
            dvtes(i,j,k,n,1) = dff(i,j,k,n) / dx
          end do
        end do
      end do
    end do
  
    deallocate(ff)
  
    !-------------------------------------------------------------------
    ! y-direction gradient
    !-------------------------------------------------------------------
    allocate(ff(0:im,-hm:jm+hm,0:km,ncolm))
    ff(0:im,:,0:km,:) = vtest(0:im,:,0:km,:)
  
    do n = 1, ncolm
      do k = 0, km
        do i = 0, im
          dff(i,:,k,n) = fds%central(fds_compact_j, f=ff(i,:,k,n), dim=jm)
        end do
      end do
    end do
  
    do k = 0, km
      do j = 0, jm
        do i = 0, im
          do n = 1, ncolm
            dvtes(i,j,k,n,2) = dff(i,j,k,n) / dy
          end do
        end do
      end do
    end do
  
    deallocate(ff)
  
    i=im/4
    k=km/4
    open(18,file='testout/grad_y'//mpirankname//'.dat')
    write(18,'(3(1X,A15))')'y','f','df'
    write(18,'(3(1X,E15.7E3))')(x(i,j,k,2),dvref(i,j,k,1,2),dvtes(i,j,k,1,2),j=0,jm)
    close(18)
    print*,' << testout/grad_y'//mpirankname//'.dat'
  
    !-------------------------------------------------------------------
    ! z-direction gradient
    !-------------------------------------------------------------------
    allocate(ff(0:im,0:jm,-hm:km+hm,ncolm))
    ff(0:im,0:jm,:,:) = vtest(0:im,0:jm,:,:)
  
    do n = 1, ncolm
      do j = 0, jm
        do i = 0, im
          dff(i,j,:,n) = fds%central(fds_compact_k, f=ff(i,j,:,n), dim=km)
        end do
      end do
    end do
  
    do k = 0, km
      do j = 0, jm
        do i = 0, im
          do n = 1, ncolm
            dvtes(i,j,k,n,3) = dff(i,j,k,n) / dz
          end do
        end do
      end do
    end do
  
    deallocate(ff)
  
    !-------------------------------------------------------------------
    ! Compute L2 error
    !-------------------------------------------------------------------
    error = 0.d0
    do k = 1, km
      do j = 1, jm
        do i = 1, im
          do n = 1, ncolm
            error(n,1) = error(n,1) + (dvtes(i,j,k,n,1) - dvref(i,j,k,n,1))**2
            error(n,2) = error(n,2) + (dvtes(i,j,k,n,2) - dvref(i,j,k,n,2))**2
            error(n,3) = error(n,3) + (dvtes(i,j,k,n,3) - dvref(i,j,k,n,3))**2
          end do
        end do
      end do
    end do
  
    error = psum(error)
  
    if (mpirank == 0) then
      error = sqrt(error / (ia * ja * ka))

      write(*,*) ' +------------------------------------------------------------+'
      write(*,*) ' | errors of gradient calculation by fds%central              |'
      write(*,*) ' +-------------------------------+----------------------------+'
      write(*,*) ' | Error in x gradient           |', error(:,1),' |'
      write(*,*) ' | Error in y gradient (non-homo)|', error(:,2),' |'
      write(*,*) ' | Error in z gradient           |', error(:,3),' |'
      write(*,'(A,F11.6,A)') '  | Time cost                     |', ptime() - time_beg,' s               |'
      write(*,*) ' +-------------------------------+----------------------------+'

    end if
  
    ! Uncomment for Tecplot output
    ! call tecbin('testout/tectest'//mpirankname//'.plt',                 &
    !             x(0:im,0:jm,0:km,1), 'x',                              &
    !             x(0:im,0:jm,0:km,2), 'y',                              &
    !             x(0:im,0:jm,0:km,3), 'z',                              &
    !             vtest(0:im,0:jm,0:km,1), 'v1',                         &
    !             dvtes(0:im,0:jm,0:km,1,2), 'dv1dy',                    &
    !             dvref(0:im,0:jm,0:km,1,2), 'dv1dy_ref')
  
  end subroutine gradtest
  !+-------------------------------------------------------------------+
  !| End of subroutine gradtest                                        |
  !+-------------------------------------------------------------------+

  !+-------------------------------------------------------------------+
  !| Subroutine: fluxtest                                              |
  !|                                                                   |
  !| Purpose:                                                          |
  !|   This subroutine tests the accuracy of flux-difference-based     |
  !|   derivative approximations using high-order compact flux schemes.|
  !|                                                                   |
  !| Functionality:                                                    |
  !|   - Constructs a 3D periodic analytical field                     |
  !|   - Computes reference (analytical) derivatives in x, y, z        |
  !|   - Applies upwind flux schemes in each direction to evaluate     |
  !|     numerical derivatives using flux differences                  |
  !|   - Computes L2 error between numerical and analytical results    |
  !|                                                                   |
  !| Author: Fang Jian                                                 |
  !| Date  : 2021                                                      |
  !+-------------------------------------------------------------------+
  subroutine fluxtest
  
    use commvar,   only : im, jm, km, hm, ia, ja, ka, conschm, difschm, &
                          lihomo, ljhomo, lkhomo, bfacmpld
    use commarray, only : x
    use comsolver, only : solvrinit
    use parallel,  only : mpiinitial, mpistop, mpisizedis, lio,         &
                          parallelini, parapp, dataswap, mpirank,       &
                          mpirankname, ptime, psum
    use solver,    only : refcal
    use flux,      only : flux_compact, flux_uw_i, flux_dw_i,          &
                          flux_uw_j, flux_dw_j, flux_uw_k, flux_dw_k
    use gridgeneration, only : gridcube
    use tecio
  
    implicit none
  
    ! Local variables
    integer :: i, j, k, n, ncolm, counter
    real(8) :: dx, dy, dz, time_beg
    real(8), save :: subtime = 0.d0
    real(8), allocatable :: vtest(:,:,:,:), dvtes(:,:,:,:,:), dvref(:,:,:,:,:)
    real(8), allocatable :: ff(:,:,:,:), dff(:,:,:,:), fh(:,:,:,:)
    real(8), allocatable :: error(:,:)
  
    !---------------------------------------------------------------
    ! Initialization
    !---------------------------------------------------------------
    conschm = '543c'
    difschm = '643c'
  
    ia = 128
    ja = 128
    ka = 128
  
    lihomo = .true.
    ljhomo = .false.
    lkhomo = .true.

    bfacmpld=1.d0
  
    call mpisizedis
    call parapp
    call parallelini
    call refcal
    call solvrinit
  
    allocate(x(-hm:im+hm, -hm:jm+hm, -hm:km+hm, 1:3))
    call gridcube(2.d0*pi, 2.d0, 2.d0*pi)
  
    dx = x(1,0,0,1) - x(0,0,0,1)
    dy = x(0,1,0,2) - x(0,0,0,2)
    dz = x(0,0,1,3) - x(0,0,0,3)
  
    ncolm = 1
    allocate(vtest(-hm:im+hm,-hm:jm+hm,-hm:km+hm,ncolm))
    allocate(dvtes(0:im,0:jm,0:km,ncolm,3), dvref(0:im,0:jm,0:km,ncolm,3))
    allocate(error(ncolm,3))
    dvtes = 0.d0
  
    time_beg = ptime()
    counter = 0
  
    !---------------------------------------------------------------
    ! Construct test function and analytical gradient
    !---------------------------------------------------------------
    do k = 0, km
      do j = 0, jm
        do i = 0, im
          vtest(i,j,k,1)   = sin(x(i,j,k,1)) * sin(0.5d0*pi*x(i,j,k,2)) * sin(x(i,j,k,3))
          dvref(i,j,k,1,1) = cos(x(i,j,k,1)) * sin(0.5d0*pi*x(i,j,k,2)) * sin(x(i,j,k,3))
          dvref(i,j,k,1,2) = 0.5d0*pi*sin(x(i,j,k,1)) * cos(0.5d0*pi*x(i,j,k,2)) * sin(x(i,j,k,3))
          dvref(i,j,k,1,3) = sin(x(i,j,k,1)) * sin(0.5d0*pi*x(i,j,k,2)) * cos(x(i,j,k,3))
        end do
      end do
    end do
  
    call dataswap(vtest)
  
    !---------------------------------------------------------------
    ! X-direction flux computation
    !---------------------------------------------------------------
    allocate(dff(0:im,0:jm,0:km,ncolm), fh(-1:im,-1:jm,-1:km,ncolm))
    allocate(ff(-hm:im+hm,0:jm,0:km,ncolm))
    ff(:,0:jm,0:km,:) = vtest(:,0:jm,0:km,:)
  
    do n = 1, 1
      do k = 0, km
        do j = 0, jm
          fh(:,j,k,n) = flux_compact(asolver=flux_uw_i, f=ff(:,j,k,n), dim=im)
          do i = 0, im
            dff(i,j,k,n) = fh(i,j,k,n) - fh(i-1,j,k,n)
          end do
        end do
      end do
    end do
  
    do k = 0, km
      do j = 0, jm
        do i = 0, im
          dvtes(i,j,k,1,1) = dff(i,j,k,1) / dx
        end do
      end do
    end do
    deallocate(ff)
  
    !---------------------------------------------------------------
    ! Y-direction flux computation
    !---------------------------------------------------------------
    allocate(ff(0:im,-hm:jm+hm,0:km,ncolm))
    ff(0:im,:,0:km,:) = vtest(0:im,:,0:km,:)
  
    do n = 1, 1
      do k = 0, km
        do i = 0, im
          fh(i,:,k,n) = flux_compact(asolver=flux_uw_j, f=ff(i,:,k,n), dim=jm)
          do j = 0, jm
            dff(i,j,k,n) = fh(i,j,k,n) - fh(i,j-1,k,n)
          end do
        end do
      end do
    end do
  
    do k = 0, km
      do j = 0, jm
        do i = 0, im
          dvtes(i,j,k,1,2) = dff(i,j,k,1) / dy
        end do
      end do
    end do
    deallocate(ff)

    i=im/4
    k=km/4
    open(18,file='testout/flux_y'//mpirankname//'.dat')
    write(18,'(3(1X,A15))')'y','f','df'
    write(18,'(3(1X,E15.7E3))')(x(i,j,k,2),dvref(i,j,k,1,2),dvtes(i,j,k,1,2),j=0,jm)
    close(18)
    print*,' << testout/flux_y'//mpirankname//'.dat'
  
    !---------------------------------------------------------------
    ! Z-direction flux computation
    !---------------------------------------------------------------
    allocate(ff(0:im,0:jm,-hm:km+hm,ncolm))
    ff(0:im,0:jm,:,:) = vtest(0:im,0:jm,:,:)
  
    do n = 1, 1
      do j = 0, jm
        do i = 0, im
          fh(i,j,:,n) = flux_compact(asolver=flux_uw_k, f=ff(i,j,:,n), dim=km)
          do k = 0, km
            dff(i,j,k,n) = fh(i,j,k,n) - fh(i,j,k-1,n)
          end do
        end do
      end do
    end do
  
    do k = 0, km
      do j = 0, jm
        do i = 0, im
          dvtes(i,j,k,1,3) = dff(i,j,k,1) / dz
        end do
      end do
    end do
    deallocate(ff)
  
    !---------------------------------------------------------------
    ! Compute L2 error against analytical gradient
    !---------------------------------------------------------------
    error = 0.d0
    do k = 1, km
      do j = 1, jm
        do i = 1, im
          error(1,1) = error(1,1) + (dvtes(i,j,k,1,1) - dvref(i,j,k,1,1))**2
          error(1,2) = error(1,2) + (dvtes(i,j,k,1,2) - dvref(i,j,k,1,2))**2
          error(1,3) = error(1,3) + (dvtes(i,j,k,1,3) - dvref(i,j,k,1,3))**2
        end do
      end do
    end do
  
    error = psum(error)
  
    if (mpirank == 0) then

      error = sqrt(error / (ia * ja * ka))

      write(*,*) ' +------------------------------------------------------------+'
      write(*,*) ' | errors of gradient calculation by flux_compact             |'
      write(*,*) ' +-------------------------------+----------------------------+'
      write(*,*) ' | Error in x gradient           |', error(:,1),' |'
      write(*,*) ' | Error in y gradient (non-homo)|', error(:,2),' |'
      write(*,*) ' | Error in z gradient           |', error(:,3),' |'
      write(*,'(A,F11.6,A)') '  | Time cost                     |', ptime() - time_beg,' s               |'
      write(*,*) ' +-------------------------------+----------------------------+'

    end if
  
    subtime = subtime + ptime() - time_beg
  
  end subroutine fluxtest
  !+-------------------------------------------------------------------+
  !| End of subroutine fluxtest                                        |
  !+-------------------------------------------------------------------+


end module test
!+---------------------------------------------------------------------+
!| The end of the module test.                                         |
!+---------------------------------------------------------------------+