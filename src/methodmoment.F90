!+---------------------------------------------------------------------+
!| This module contains subroutines related to method of moment.       |
!| ==============                                                      |
!| CHANGE RECORD                                                       |
!| -------------                                                       |
!| 15-08-2023: Created by J. Fang @ STFC Daresbury Laboratory          |
!+---------------------------------------------------------------------+
module methodmoment
  !
  use constdef
  use commvar,  only : im,jm,km,hm,ctime,lreport,ltimrpt
  use parallel, only : ptime,lio,mpirank,mpistop
  use utility,  only : timereporter
  !
  implicit none
  !
  ! real(8),parameter ::  Am=1.d0,AR1=1.d0,Aq=1.d0,Asigma=1.d0,          &
  !                       Adelta1=1.d0,Aphi=1.d0 ! BGK model
  real(8),parameter ::  Asigma=1.d0,Aq=num2d3,Am=1.5d0,AR1=num7d6,    &
                        AR2 = num2d3,                                 & 
                        Adelta1=num2d3,Adelta2=num2d3,Aphi1=2.097d0,  &
                        Apsi1=1.698d9,Aomega1=1.d0 ! Maxwell Molecule model
  !
  !
  ! real(8),allocatable,dimension(:,:,:) :: mijk,1,mijk,2,mijk,3,mijk,4,mijk,5,    &
  !                                         mijk,6,mijk,7
  !
  real(8),allocatable,dimension(:,:,:,:) :: dqx,dqy,dqz
  real(8),allocatable,dimension(:,:,:,:) :: dsigmaxx,dsigmaxy,       &
                          dsigmaxz,dsigmayy,dsigmayz,dsigmazz,       &
                          dprs,ddelta
  real(8),allocatable,dimension(:,:,:,:,:) :: dRij, domega, dmijk

  real(8),allocatable,dimension(:,:,:,:) :: Eijkl
  real(8),allocatable,dimension(:,:,:,:) :: psi
  real(8),allocatable,dimension(:,:,:,:) :: omega
  real(8),allocatable,dimension(:,:,:) :: delta
  real(8),allocatable,dimension(:,:,:,:) :: mijk,rij
  !
  real(8),allocatable,dimension(:,:,:,:) :: q_mom,qrhs_mom
  !
  interface fvar2qmom
     module procedure sigma_flux2qmom_3da
     module procedure sigma_flux2qmom_1da
     module procedure sigma_flux2qmom_sca
  end interface
  !
  interface qmom2fvar
     module procedure qmom2sigma_flux_3da
     module procedure qmom2sigma_flux_1da
     module procedure qmom2sigma_flux_sca
  end interface
  !
  contains
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to allocate moment equations.                  |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 10-04-2024: Created by J. Fang @ STFC Daresbury Laboratory        |
  !+-------------------------------------------------------------------+
  subroutine allo_moment
    !
    use commvar, only: im,jm,km,hm,num_xtrmom
    !
    allocate(q_mom(-hm:im+hm,-hm:jm+hm,-hm:km+hm,1:num_xtrmom),source=0.d0)
    allocate(qrhs_mom(0:im,0:jm,0:km,1:num_xtrmom),source=0.d0)
    allocate( dprs(0:im,0:jm,0:km,1:3),source=0.d0)
    allocate(  mijk(-hm:im+hm,-hm:jm+hm,-hm:km+hm,1:7),source=0.d0)
    allocate(   rij(-hm:im+hm,-hm:jm+hm,-hm:km+hm,1:5),source=0.d0)
    !
    allocate( dsigmaxx(0:im,0:jm,0:km,1:3),   &
              dsigmaxy(0:im,0:jm,0:km,1:3),   &
              dsigmaxz(0:im,0:jm,0:km,1:3),   &
              dsigmayy(0:im,0:jm,0:km,1:3),   &
              dsigmayz(0:im,0:jm,0:km,1:3),   &
              dsigmazz(0:im,0:jm,0:km,1:3),source=0.d0    )
    allocate(delta(-hm:im+hm,-hm:jm+hm,-hm:km+hm),source=0.d0)
    allocate( ddelta(0:im,0:jm,0:km,1:3),source=0.d0 )
    !
    allocate( dqx(0:im,0:jm,0:km,1:3),source=0.d0 )
    allocate( dqy(0:im,0:jm,0:km,1:3),source=0.d0 )
    allocate( dqz(0:im,0:jm,0:km,1:3),source=0.d0 )
    !
    allocate(Eijkl(-hm:im+hm,-hm:jm+hm,-hm:km+hm,1:9),source=0.d0)
    allocate(psi(-hm:im+hm,-hm:jm+hm,-hm:km+hm,1:7),source=0.d0)
    allocate(omega(-hm:im+hm,-hm:jm+hm,-hm:km+hm,1:7),source=0.d0)
    allocate(dmijk(0:im,0:jm,0:km,1:3,1:7),source=0.d0)
    allocate(dRij(0:im,0:jm,0:km,1:3,1:6),source=0.d0)
    allocate(domega(0:im,0:jm,0:km,1:3,1:3),source=0.d0)
  end subroutine allo_moment
  !+-------------------------------------------------------------------+
  !| The end of the subroutine allo_moment.                            |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine advances the moment equation in time using 3-step |
  !| 3rd-rder Rungle-Kutta scheme.                                     |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 27-Nov-2018: Created by J. Fang @ STFC Daresbury Laboratory       |
  !| 15-Aug-2023: Move to this module, J. Fang @ STFC                  |
  !+-------------------------------------------------------------------+
  subroutine rk3mom(timestep,timerept)
    !
    use commvar,  only : num_xtrmom,deltat,lfilter,feqchkpt,lavg,feqavg, &
                         limmbou,turbmode,feqslice,feqwsequ,lwslic,      &
                         flowtype,num_species,maxstep,ltimrpt
    use commarray,only : x,jacob
    use bc,       only : bctype
    use comsolver,only : filterq,spongefilter,filter2e

    implicit none

    !
    ! argument
    real(8),intent(in) :: timestep
    logical,intent(in),optional :: timerept
    !
    ! local data
    logical,save :: firstcall = .true.
    real(8),save :: rkcoe(3,3)
    integer :: nrk,i,j,k,m,n
    real(8) :: time_beg,time_beg_rhs,time_beg_sta,time_beg_io
    real(8),allocatable,save :: qsave(:,:,:,:)
    integer :: dt_ratio,jdnn,idnn
    real(8) :: hrr,time_beg_2
    real(8),save :: subtime=0.d0
    !
    time_beg=ptime()
    !
    if(firstcall) then
      !
      rkcoe(1,1)=1.d0
      rkcoe(2,1)=0.d0
      rkcoe(3,1)=1.d0
      !
      rkcoe(1,2)=0.75d0
      rkcoe(2,2)=0.25d0
      rkcoe(3,2)=0.25d0
      !
      rkcoe(1,3)=num1d3
      rkcoe(2,3)=num2d3
      rkcoe(3,3)=num2d3
      !
      allocate(qsave(0:im,0:jm,0:km,1:num_xtrmom),source=0.d0)
      !
      firstcall=.false.
      !
    endif
    !
    do nrk=1,3
      !
      call qmomswap(timerept=ltimrpt)
      !
      !----------------------------
      !-- wall boundary condition -
      !----------------------------
      do n = 1, 6
         if(bctype(n)==413) then
           call MOM_wall_boundary(n)
         end if
      end do
      !
      call rhsmomcal(timerept=ltimrpt)
      !
      time_beg_2=ptime()
      !
      if(nrk==1) then
        !
        do m=1,num_xtrmom
          qsave(0:im,0:jm,0:km,m)=q_mom(0:im,0:jm,0:km,m)
        enddo
        !
      endif
      !
      do m=1,num_xtrmom
        !
        q_mom(0:im,0:jm,0:km,m)=rkcoe(1,nrk)*qsave(0:im,0:jm,0:km,m)+    &
                                rkcoe(2,nrk)*q_mom(0:im,0:jm,0:km,m)    +    &
                                rkcoe(3,nrk)*qrhs_mom(0:im,0:jm,0:km,m)/jacob(0:im,0:jm,0:km)*timestep
        !
      enddo
      !
      ctime(14)=ctime(14)+ptime()-time_beg_2
      !
      if(lfilter) call filterq(q_mom,num_xtrmom,timerept=ltimrpt)
      !
      time_beg_2=ptime()
      !
      call updatemoment(0,im,0,jm,0,km)
      !
      ctime(15)=ctime(15)+ptime()-time_beg_2
      !
    enddo !rk
    !
    if(present(timerept) .and. timerept) then
      !
      subtime=subtime+ptime()-time_beg
      !
    endif
    !
    ! call mpistop
    !
    return
    !
  end subroutine rk3mom
  !+-------------------------------------------------------------------+
  !| The end of the subroutine rk3mom.                                 |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to calculate the convectional terms of    |
  !| moment transport equations.                                       |
  !+-------------------------------------------------------------------+
  subroutine convection_mom(timerept)
    !
    use commvar,  only: im,jm,km,hm,conschm,npdci,npdcj,npdck,         &
                        is,ie,js,je,ks,ke,num_xtrmom,lfftk,ndims
    use commarray,only: vel,dxi,jacob
    use commfunc, only: ddfc
    use comsolver,only : cci,ccj,cck,alfa_con

    implicit none
    !
    ! arguments
    logical,intent(in),optional :: timerept
    !
    ! local data
    integer :: i,j,k,jmom
    real(8),allocatable :: fcs(:,:),dfcs(:,:),uu(:)
    !
    real(8) :: time_beg
    real(8),save :: subtime=0.d0
    !
    if(present(timerept) .and. timerept) time_beg=ptime() 
    !

    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! calculating along i direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    allocate(fcs(-hm:im+hm,1:num_xtrmom),dfcs(0:im,1:num_xtrmom),uu(-hm:im+hm))
    ! !
    fcs=0.d0; dfcs=0.d0
    do k=ks,ke
    do j=js,je
      !
      uu(:)=dxi(:,j,k,1,1)*vel(:,j,k,1)+dxi(:,j,k,1,2)*vel(:,j,k,2) +  &
            dxi(:,j,k,1,3)*vel(:,j,k,3)
      !
      do jmom=1,num_xtrmom
        fcs(:,jmom)=jacob(:,j,k)*q_mom(:,j,k,jmom)*uu
      enddo
      !
      do jmom=1,num_xtrmom
        dfcs(:,jmom)=ddfc(fcs(:,jmom),conschm,npdci,im,alfa_con,cci)
      enddo
      !
      qrhs_mom(is:ie,j,k,:)=qrhs_mom(is:ie,j,k,:)+dfcs(is:ie,:) 
      !
    enddo
    enddo
    deallocate(fcs,dfcs,uu)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! end calculating along i direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! calculating along j direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    allocate(fcs(-hm:jm+hm,1:num_xtrmom),dfcs(0:jm,1:num_xtrmom),uu(-hm:jm+hm))
    fcs=0.d0; dfcs=0.d0
    do k=ks,ke
    do i=is,ie
      !
      uu(:)=dxi(i,:,k,2,1)*vel(i,:,k,1)+dxi(i,:,k,2,2)*vel(i,:,k,2) +  &
            dxi(i,:,k,2,3)*vel(i,:,k,3)
      !
      do jmom=1,num_xtrmom
        fcs(:,jmom)=jacob(i,:,k)*q_mom(i,:,k,jmom)*uu
      enddo
      !
      do jmom=1,num_xtrmom
        dfcs(:,jmom)=ddfc(fcs(:,jmom),conschm,npdcj,jm,alfa_con,ccj)
      enddo
      !
      qrhs_mom(i,js:je,k,:)=qrhs_mom(i,js:je,k,:)+dfcs(js:je,:)
      !
      !
    enddo
    enddo
    deallocate(fcs,dfcs,uu)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! end calculating along j direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    if(ndims==3) then
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! calculating along j direction
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      allocate(fcs(-hm:km+hm,1:num_xtrmom),dfcs(0:km,1:num_xtrmom),uu(-hm:km+hm))
      fcs=0.d0; dfcs=0.d0
      do j=js,je
      do i=is,ie
        !
        uu(:)=dxi(i,j,:,3,1)*vel(i,j,:,1)+dxi(i,j,:,3,2)*vel(i,j,:,2) +  &
              dxi(i,j,:,3,3)*vel(i,j,:,3)
        !
        do jmom=1,num_xtrmom
          fcs(:,jmom)=jacob(i,j,:)*q_mom(i,j,:,jmom)*uu
        enddo
        !
        do jmom=1,num_xtrmom
          dfcs(:,jmom)=ddfc(fcs(:,jmom),conschm,npdck,km,alfa_con,cck,lfft=lfftk)
        enddo
        !
        qrhs_mom(i,j,ks:ke,:)=qrhs_mom(i,j,ks:ke,:)+dfcs(ks:ke,:)
        !
      enddo
      enddo
      deallocate(fcs,dfcs,uu)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! end calculating along j direction
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
    endif
    !
    if(present(timerept) .and. timerept) then
      !
      subtime=subtime+ptime()-time_beg
      !
      if(lio .and. lreport) call timereporter(routine='convection_mom', &
                                             timecost=subtime, &
                                              message='diffusion term with central scheme')
    endif
    !
    return
    !
  end subroutine convection_mom
  !+-------------------------------------------------------------------+
  !| The of the subroutine convection_mom                              |
  !+-------------------------------------------------------------------+
  !!
  !+-------------------------------------------------------------------+
  !| This subroutine is to calculate the diffusion term from the r13   |
  !| moment equations.                                                 |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 01-11-2022: Created by J. Fang @ STFC Daresbury Laboratory        |
  !| 27-01-2023: Add heat flux term by J. Fang @ STFC DL               |
  !+-------------------------------------------------------------------+
  subroutine diffu_sigma_q
    !
    use constdef
    use commvar,  only : im,jm,km,gamma,Reynolds,mach,difschm,         &
                         npdci,npdcj,npdck,num_species,num_modequ,     &
                         is,ie,js,je,ks,ke,moment,lfftk,ndims
    use commarray,only : tmp,prs,sigma,qflux,dxi,jacob
    use commfunc, only : ddfc
    use tecio
    use comsolver,only: dci,dcj,dck,alfa_dif
    !
    ! local data
    !
    implicit none

    real(8),allocatable,dimension(:,:,:) :: tsgmxx,tsgmxy,tsgmxz,      &
                                            tsgmyy,tsgmyz,tsgmzz
    real(8),allocatable,dimension(:,:,:,:) :: dtsgmxx,dtsgmxy,dtsgmxz, &
                                              dtsgmyy,dtsgmyz,dtsgmzz
    real(8),allocatable :: df(:,:),ff(:,:)
    integer :: i,j,k,n
    real(8) :: co1,co2
    logical,save :: firstcall=.true.
    !
    allocate(ff(-hm:im+hm,1:8),df(0:im,1:8))
    do k=0,km
    do j=0,jm
      !
      ! sigmaxx
      ff(:,1)=( mijk(:,j,k,1)*dxi(:,j,k,1,1) +  &
                mijk(:,j,k,2)*dxi(:,j,k,1,2) +  &
                mijk(:,j,k,3)*dxi(:,j,k,1,3) )*jacob(:,j,k)
      ! sigmaxy
      ff(:,2)=( mijk(:,j,k,2)*dxi(:,j,k,1,1) +  &
                mijk(:,j,k,4)*dxi(:,j,k,1,2) +  &
                mijk(:,j,k,7)*dxi(:,j,k,1,3) )*jacob(:,j,k)
      ! sigmaxz
      ff(:,3)=( mijk(:,j,k,3)*dxi(:,j,k,1,1) +  &
                mijk(:,j,k,7)*dxi(:,j,k,1,2) +  &
              (-mijk(:,j,k,1)-mijk(:,j,k,4))*dxi(:,j,k,1,3) )*jacob(:,j,k)
      ! sigmayy
      ff(:,4)=( mijk(:,j,k,4)*dxi(:,j,k,1,1) +  &
                mijk(:,j,k,5)*dxi(:,j,k,1,2) +  &
                mijk(:,j,k,6)*dxi(:,j,k,1,3) )*jacob(:,j,k)
      ! sigmayz
      ff(:,5)=( mijk(:,j,k,7)*dxi(:,j,k,1,1) +  &
                mijk(:,j,k,6)*dxi(:,j,k,1,2) +  &
              (-mijk(:,j,k,2)-mijk(:,j,k,5))*dxi(:,j,k,1,3) )*jacob(:,j,k)
      !
      ! qx
      ff(:,6)=( Rij(:,j,k,1)*dxi(:,j,k,1,1) +  &
                Rij(:,j,k,2)*dxi(:,j,k,1,2) +  &
                Rij(:,j,k,3)*dxi(:,j,k,1,3) )*0.5d0*jacob(:,j,k)
      ! qy
      ff(:,7)=( Rij(:,j,k,2)*dxi(:,j,k,1,1) +  &
                Rij(:,j,k,4)*dxi(:,j,k,1,2) +  &
                Rij(:,j,k,5)*dxi(:,j,k,1,3) )*0.5d0*jacob(:,j,k)
      ! qz
      ff(:,8)=( Rij(:,j,k,3)*dxi(:,j,k,1,1) +  &
                Rij(:,j,k,5)*dxi(:,j,k,1,2) +  &
                (-Rij(:,j,k,1)-Rij(:,j,k,4))*dxi(:,j,k,1,3) )*0.5d0*jacob(:,j,k)
      !
      do n=1,8
        df(:,n)=ddfc(ff(:,n),difschm,npdci,im,alfa_dif,dci)
      enddo
      !
      qrhs_mom(is:ie,j,k,1)=qrhs_mom(is:ie,j,k,1)-df(is:ie,1)
      qrhs_mom(is:ie,j,k,2)=qrhs_mom(is:ie,j,k,2)-df(is:ie,2)
      qrhs_mom(is:ie,j,k,3)=qrhs_mom(is:ie,j,k,3)-df(is:ie,3)
      qrhs_mom(is:ie,j,k,4)=qrhs_mom(is:ie,j,k,4)-df(is:ie,4)
      qrhs_mom(is:ie,j,k,5)=qrhs_mom(is:ie,j,k,5)-df(is:ie,5)
      !
      qrhs_mom(is:ie,j,k,6)=qrhs_mom(is:ie,j,k,6)-df(is:ie,6)
      qrhs_mom(is:ie,j,k,7)=qrhs_mom(is:ie,j,k,7)-df(is:ie,7)
      qrhs_mom(is:ie,j,k,8)=qrhs_mom(is:ie,j,k,8)-df(is:ie,8)
      !
    enddo
    enddo
    deallocate(ff,df)
    !
    allocate(ff(-hm:jm+hm,1:8),df(0:jm,1:8))
    do k=0,km
    do i=0,im
      !
      ! sigmaxx
      ff(:,1)=( mijk(i,:,k,1)*dxi(i,:,k,2,1) +  &
                mijk(i,:,k,2)*dxi(i,:,k,2,2) +  &
                mijk(i,:,k,3)*dxi(i,:,k,2,3) )*jacob(i,:,k)
      ! sigmaxy
      ff(:,2)=( mijk(i,:,k,2)*dxi(i,:,k,2,1) +  &
                mijk(i,:,k,4)*dxi(i,:,k,2,2) +  &
                mijk(i,:,k,7)*dxi(i,:,k,2,3) )*jacob(i,:,k)
      ! sigmaxz
      ff(:,3)=( mijk(i,:,k,3)*dxi(i,:,k,2,1) +  &
                mijk(i,:,k,7)*dxi(i,:,k,2,2) +  &
              (-mijk(i,:,k,1)-mijk(i,:,k,4))*dxi(i,:,k,2,3) )*jacob(i,:,k)
      ! sigmayy
      ff(:,4)=( mijk(i,:,k,4)*dxi(i,:,k,2,1) +  &
                mijk(i,:,k,5)*dxi(i,:,k,2,2) +  &
                mijk(i,:,k,6)*dxi(i,:,k,2,3) )*jacob(i,:,k)
      ! sigmayz
      ff(:,5)=( mijk(i,:,k,7)*dxi(i,:,k,2,1) +  &
                mijk(i,:,k,6)*dxi(i,:,k,2,2) +  &
              (-mijk(i,:,k,2)-mijk(i,:,k,5))*dxi(i,:,k,2,3) )*jacob(i,:,k)
      !
      ! qx
      ff(:,6)=( Rij(i,:,k,1)*dxi(i,:,k,2,1) +  &
                Rij(i,:,k,2)*dxi(i,:,k,2,2) +  &
                Rij(i,:,k,3)*dxi(i,:,k,2,3) )*0.5d0*jacob(i,:,k)
      ! qy
      ff(:,7)=( Rij(i,:,k,2)*dxi(i,:,k,2,1) +  &
                Rij(i,:,k,4)*dxi(i,:,k,2,2) +  &
                Rij(i,:,k,5)*dxi(i,:,k,2,3) )*0.5d0*jacob(i,:,k)
      ! qz
      ff(:,8)=( Rij(i,:,k,3)*dxi(i,:,k,2,1) +  &
                Rij(i,:,k,5)*dxi(i,:,k,2,2) +  &
              (-Rij(i,:,k,1)-Rij(i,:,k,4))*dxi(i,:,k,2,3) )*0.5d0*jacob(i,:,k)
      !
      do n=1,8
        df(:,n)=ddfc(ff(:,n),difschm,npdcj,jm,alfa_dif,dcj)
      enddo
      !
      qrhs_mom(i,js:je,k,1)=qrhs_mom(i,js:je,k,1)-df(js:je,1)
      qrhs_mom(i,js:je,k,2)=qrhs_mom(i,js:je,k,2)-df(js:je,2)
      qrhs_mom(i,js:je,k,3)=qrhs_mom(i,js:je,k,3)-df(js:je,3)
      qrhs_mom(i,js:je,k,4)=qrhs_mom(i,js:je,k,4)-df(js:je,4)
      qrhs_mom(i,js:je,k,5)=qrhs_mom(i,js:je,k,5)-df(js:je,5)
      !
      qrhs_mom(i,js:je,k,6)=qrhs_mom(i,js:je,k,6)-df(js:je,6)
      qrhs_mom(i,js:je,k,7)=qrhs_mom(i,js:je,k,7)-df(js:je,7)
      qrhs_mom(i,js:je,k,8)=qrhs_mom(i,js:je,k,8)-df(js:je,8)
      !
    enddo
    enddo
    deallocate(ff,df)
    !
    if(ndims==3) then
      !
      allocate(ff(-hm:km+hm,1:8),df(0:km,1:8))
      do j=0,jm
      do i=0,im
        !
        ! sigmaxx
        ff(:,1)=( mijk(i,j,:,1)*dxi(i,j,:,3,1) +  &
                  mijk(i,j,:,2)*dxi(i,j,:,3,2) +  &
                  mijk(i,j,:,3)*dxi(i,j,:,3,3) )*jacob(i,j,:)
        ! sigmaxy
        ff(:,2)=( mijk(i,j,:,2)*dxi(i,j,:,3,1) +  &
                  mijk(i,j,:,4)*dxi(i,j,:,3,2) +  &
                  mijk(i,j,:,7)*dxi(i,j,:,3,3) )*jacob(i,j,:)
        ! sigmaxz
        ff(:,3)=( mijk(i,j,:,3)*dxi(i,j,:,3,1) +  &
                  mijk(i,j,:,7)*dxi(i,j,:,3,2) +  &
                (-mijk(i,j,:,1)-mijk(i,j,:,4))*dxi(i,j,:,3,3) )*jacob(i,j,:)
        ! sigmayy
        ff(:,4)=( mijk(i,j,:,4)*dxi(i,j,:,3,1) +  &
                  mijk(i,j,:,5)*dxi(i,j,:,3,2) +  &
                  mijk(i,j,:,6)*dxi(i,j,:,3,3) )*jacob(i,j,:)
        ! sigmayz
        ff(:,5)=( mijk(i,j,:,7)*dxi(i,j,:,3,1) +  &
                  mijk(i,j,:,6)*dxi(i,j,:,3,2) +  &
                (-mijk(i,j,:,2)-mijk(i,j,:,5))*dxi(i,j,:,3,3) )*jacob(i,j,:)
        !
        ! qx
        ff(:,6)=( Rij(i,j,:,1)*dxi(i,j,:,3,1) +  &
                  Rij(i,j,:,2)*dxi(i,j,:,3,2) +  &
                  Rij(i,j,:,3)*dxi(i,j,:,3,3) )*0.5d0*jacob(i,j,:)
        ! qy
        ff(:,7)=( Rij(i,j,:,2)*dxi(i,j,:,3,1) +  &
                  Rij(i,j,:,4)*dxi(i,j,:,3,2) +  &
                  Rij(i,j,:,5)*dxi(i,j,:,3,3) )*0.5d0*jacob(i,j,:)
        ! qz
        ff(:,8)=( Rij(i,j,:,3)*dxi(i,j,:,3,1) +  &
                  Rij(i,j,:,5)*dxi(i,j,:,3,2) +  &
                (-Rij(i,j,:,1)-Rij(i,j,:,4))*dxi(i,j,:,3,3) )*0.5d0*jacob(i,j,:)
        !
        do n=1,8
          df(:,n)=ddfc(ff(:,n),difschm,npdck,km,alfa_dif,dck,lfft=lfftk)
        enddo
        !
        qrhs_mom(i,j,ks:ke,1)=qrhs_mom(i,j,ks:ke,1)-df(ks:ke,1)
        qrhs_mom(i,j,ks:ke,2)=qrhs_mom(i,j,ks:ke,2)-df(ks:ke,2)
        qrhs_mom(i,j,ks:ke,3)=qrhs_mom(i,j,ks:ke,3)-df(ks:ke,3)
        qrhs_mom(i,j,ks:ke,4)=qrhs_mom(i,j,ks:ke,4)-df(ks:ke,4)
        qrhs_mom(i,j,ks:ke,5)=qrhs_mom(i,j,ks:ke,5)-df(ks:ke,5)
        !
        qrhs_mom(i,j,ks:ke,6)=qrhs_mom(i,j,ks:ke,6)-df(ks:ke,6)
        qrhs_mom(i,j,ks:ke,7)=qrhs_mom(i,j,ks:ke,7)-df(ks:ke,7)
        qrhs_mom(i,j,ks:ke,8)=qrhs_mom(i,j,ks:ke,8)-df(ks:ke,8)
        !
      enddo
      enddo
      deallocate(ff,df)
      !
    endif

  end subroutine diffu_sigma_q
  !+-------------------------------------------------------------------+
  !| The end of the subroutine diffu_sigma_q.                          |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to solve the diffusion term of the r26 eq.     |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 09-05-2023: Add heat flux term by J. Fang @ STFC DL               |
  !+-------------------------------------------------------------------+
  subroutine diffu_mijk
    !
    use commvar,  only : im,jm,km,hm,Reynolds,npdci,npdcj,npdck,     &
                         is,ie,js,je,ks,ke,moment,difschm,lfftk,ndims
    use commarray,only : tmp,rho,dxi,jacob,miu
    use commfunc, only : ddfc
    use comsolver,only: dci,dcj,dck,alfa_dif,grad
    use parallel, only : dataswap
    !
    real(8),allocatable :: df(:,:),ff(:,:)
    !
    integer :: i,j,k,n
    !
    if (.not. allocated(dmijk)) then
       allocate(dmijk(0:im,0:jm,0:km,1:3,1:7))
       dmijk = 0.0d0
    end if
    !
    dmijk(:,:,:,:,1)=grad(mijk(:,:,:,1))
    dmijk(:,:,:,:,2)=grad(mijk(:,:,:,2))
    dmijk(:,:,:,:,3)=grad(mijk(:,:,:,3))
    dmijk(:,:,:,:,4)=grad(mijk(:,:,:,4))
    dmijk(:,:,:,:,5)=grad(mijk(:,:,:,5))
    dmijk(:,:,:,:,6)=grad(mijk(:,:,:,6))
    dmijk(:,:,:,:,7)=grad(mijk(:,:,:,7))
    !
    if (.not. allocated(Eijkl)) then
       allocate(Eijkl(-hm:im+hm,-hm:jm+hm,-hm:km+hm,1:9))
       Eijkl = 0.0d0
    end if
    !
    !E_{xxxx}
    Eijkl(0:im,0:jm,0:km,1)= num16d7 * dmijk(:,:,:,1,1)- num12d7 * (dmijk(:,:,:,2,2)+dmijk(:,:,:,3,3))
     !
    !E_{xxxy}
    Eijkl(0:im,0:jm,0:km,2)= num15d7 * dmijk(:,:,:,1,2) + dmijk(:,:,:,2,1)  &
                           - num6d7 * (dmijk(:,:,:,2,4) + dmijk(:,:,:,3,7) )
    !
    !E_{xxxz}
    Eijkl(0:im,0:jm,0:km,3)= num15d7 * dmijk(:,:,:,1,3) - num6d7 * dmijk(:,:,:,2,7) + &
                             num13d7 * dmijk(:,:,:,3,1) + num6d7 * dmijk(:,:,:,3,4) 
    !
    !E_{xxyy}
    Eijkl(0:im,0:jm,0:km,4)= num12d7 * (dmijk(:,:,:,1,4)+dmijk(:,:,:,2,2)) - &
                              num2d7 * (dmijk(:,:,:,1,1)+dmijk(:,:,:,2,5) + &
                                        dmijk(:,:,:,3,6)+dmijk(:,:,:,3,3))
    !
    !E_{xxyz}
    Eijkl(0:im,0:jm,0:km,5)=dmijk(:,:,:,2,3) + num9d7 * dmijk(:,:,:,3,2) + &
                            num12d7 * dmijk(:,:,:,1,7) - num2d7 * (dmijk(:,:,:,2,6)-dmijk(:,:,:,3,5))
    !
    !E_{xyyy}
    Eijkl(0:im,0:jm,0:km,6)= num15d7 * dmijk(:,:,:,2,4)+dmijk(:,:,:,1,5) - &
                              num6d7 * (dmijk(:,:,:,1,2)+dmijk(:,:,:,3,7))
    !
    !E_{xyyz}
    Eijkl(0:im,0:jm,0:km,7)= num12d7 * dmijk(:,:,:,2,7) + num9d7 * dmijk(:,:,:,3,4) + &
                             dmijk(:,:,:,1,6) - num2d7 * (dmijk(:,:,:,1,3)-dmijk(:,:,:,3,1))
    !
    !E_{yyyy}
    Eijkl(0:im,0:jm,0:km,8)= num16d7 * dmijk(:,:,:,2,5) - num12d7 * (dmijk(:,:,:,1,4)+dmijk(:,:,:,3,6))
    !
    !E_{yyyz}
    Eijkl(0:im,0:jm,0:km,9)= num15d7 * dmijk(:,:,:,2,6)-dmijk(:,:,:,3,3) -dmijk(:,:,:,3,6) &
                           -  num6d7 * dmijk(:,:,:,1,7)+ num6d7 * dmijk(:,:,:,3,2)+ num6d7 * dmijk(:,:,:,3,5)
    !
    do n=1,9
      Eijkl(0:im,0:jm,0:km,n)=-miu(0:im,0:jm,0:km)/rho(0:im,0:jm,0:km)/Aphi1*Eijkl(0:im,0:jm,0:km,n)
    enddo
    !
    call dataswap(Eijkl)
    !
    allocate(ff(-hm:im+hm,1:7),df(0:im,1:7))
    do k=0,km
    do j=0,jm
      !
      ! mxxx
      ff(:,1)=( Eijkl(:,j,k,1)*dxi(:,j,k,1,1) +  &
                Eijkl(:,j,k,2)*dxi(:,j,k,1,2) +  &
                Eijkl(:,j,k,3)*dxi(:,j,k,1,3) )*jacob(:,j,k)
      ! mxxy
      ff(:,2)=( Eijkl(:,j,k,2)*dxi(:,j,k,1,1) +  &
                Eijkl(:,j,k,4)*dxi(:,j,k,1,2) +  &
                Eijkl(:,j,k,5)*dxi(:,j,k,1,3) )*jacob(:,j,k)
      ! mxxz
      ff(:,3)=( Eijkl(:,j,k,3)*dxi(:,j,k,1,1) +  &
                Eijkl(:,j,k,5)*dxi(:,j,k,1,2) +  &
              (-Eijkl(:,j,k,1)-Eijkl(:,j,k,4))*dxi(:,j,k,1,3) )*jacob(:,j,k)
      ! mxyy
      ff(:,4)=( Eijkl(:,j,k,4)*dxi(:,j,k,1,1) +  &
                Eijkl(:,j,k,6)*dxi(:,j,k,1,2) +  &
                Eijkl(:,j,k,7)*dxi(:,j,k,1,3) )*jacob(:,j,k)
      ! myyy
      ff(:,5)=( Eijkl(:,j,k,6)*dxi(:,j,k,1,1) +  &
                Eijkl(:,j,k,8)*dxi(:,j,k,1,2) +  &
                Eijkl(:,j,k,9)*dxi(:,j,k,1,3) )*jacob(:,j,k)
      ! myyz
      ff(:,6)=( Eijkl(:,j,k,7)*dxi(:,j,k,1,1) +  &
                Eijkl(:,j,k,9)*dxi(:,j,k,1,2) +  &
              (-Eijkl(:,j,k,4)-Eijkl(:,j,k,8))*dxi(:,j,k,1,3) )*jacob(:,j,k)
      ! mxyz
      ff(:,7)=( Eijkl(:,j,k,5)*dxi(:,j,k,1,1) +  &
                Eijkl(:,j,k,7)*dxi(:,j,k,1,2) +  &
              (-Eijkl(:,j,k,2)-Eijkl(:,j,k,6))*dxi(:,j,k,1,3) )*jacob(:,j,k)
      !
      !
      do n=1,7
        df(:,n)=ddfc(ff(:,n),difschm,npdci,im,alfa_dif,dci)
      enddo
      !
      do n=1,7
        qrhs_mom(is:ie,j,k,8+n)=qrhs_mom(is:ie,j,k,8+n)-df(is:ie,n)
      enddo
      !
    enddo
    enddo
    deallocate(ff,df)
    !
    allocate(ff(-hm:jm+hm,1:7),df(0:jm,1:7))
    do k=0,km
    do i=0,im
      !
      ! mxxx
      ff(:,1)=( Eijkl(i,:,k,1)*dxi(i,:,k,2,1) +  &
                Eijkl(i,:,k,2)*dxi(i,:,k,2,2) +  &
                Eijkl(i,:,k,3)*dxi(i,:,k,2,3) )*jacob(i,:,k)
      ! mxxy
      ff(:,2)=( Eijkl(i,:,k,2)*dxi(i,:,k,2,1) +  &
                Eijkl(i,:,k,4)*dxi(i,:,k,2,2) +  &
                Eijkl(i,:,k,5)*dxi(i,:,k,2,3) )*jacob(i,:,k)
      ! mxxz
      ff(:,3)=( Eijkl(i,:,k,3)*dxi(i,:,k,2,1) +  &
                Eijkl(i,:,k,5)*dxi(i,:,k,2,2) +  &
              (-Eijkl(i,:,k,1)-Eijkl(i,:,k,4))*dxi(i,:,k,2,3) )*jacob(i,:,k)
      ! mxyy
      ff(:,4)=( Eijkl(i,:,k,4)*dxi(i,:,k,2,1) +  &
                Eijkl(i,:,k,6)*dxi(i,:,k,2,2) +  &
                Eijkl(i,:,k,7)*dxi(i,:,k,2,3) )*jacob(i,:,k)
      ! myyy
      ff(:,5)=( Eijkl(i,:,k,6)*dxi(i,:,k,2,1) +  &
                Eijkl(i,:,k,8)*dxi(i,:,k,2,2) +  &
                Eijkl(i,:,k,9)*dxi(i,:,k,2,3) )*jacob(i,:,k)
      ! myyz
      ff(:,6)=( Eijkl(i,:,k,7)*dxi(i,:,k,2,1) +  &
                Eijkl(i,:,k,9)*dxi(i,:,k,2,2) +  &
              (-Eijkl(i,:,k,4)-Eijkl(i,:,k,8))*dxi(i,:,k,2,3) )*jacob(i,:,k)
      ! mxyz
      ff(:,7)=( Eijkl(i,:,k,5)*dxi(i,:,k,2,1) +  &
                Eijkl(i,:,k,7)*dxi(i,:,k,2,2) +  &
              (-Eijkl(i,:,k,2)-Eijkl(i,:,k,6))*dxi(i,:,k,2,3) )*jacob(i,:,k)
      !
      do n=1,7
        df(:,n)=ddfc(ff(:,n),difschm,npdcj,jm,alfa_dif,dcj)
      enddo
      !
      do n=1,7
        qrhs_mom(i,js:je,k,8+n)=qrhs_mom(i,js:je,k,8+n)-df(js:je,n)
      enddo
      !
    enddo
    enddo
    deallocate(ff,df)
    !
    if(ndims==3) then
      !
      allocate(ff(-hm:km+hm,1:7),df(0:km,1:7))
      do j=0,jm
      do i=0,im
        !
        ! mxxx
        ff(:,1)=( Eijkl(i,j,:,1)*dxi(i,j,:,3,1) +  &
                  Eijkl(i,j,:,2)*dxi(i,j,:,3,2) +  &
                  Eijkl(i,j,:,3)*dxi(i,j,:,3,3) )*jacob(i,j,:)
        ! mxxy
        ff(:,2)=( Eijkl(i,j,:,2)*dxi(i,j,:,3,1) +  &
                  Eijkl(i,j,:,4)*dxi(i,j,:,3,2) +  &
                  Eijkl(i,j,:,5)*dxi(i,j,:,3,3) )*jacob(i,j,:)
        ! mxxz
        ff(:,3)=( Eijkl(i,j,:,3)*dxi(i,j,:,3,1) +  &
                  Eijkl(i,j,:,5)*dxi(i,j,:,3,2) +  &
                (-Eijkl(i,j,:,1)-Eijkl(i,j,:,4))*dxi(i,j,:,3,3) )*jacob(i,j,:)
        ! mxyy
        ff(:,4)=( Eijkl(i,j,:,4)*dxi(i,j,:,3,1) +  &
                  Eijkl(i,j,:,6)*dxi(i,j,:,3,2) +  &
                  Eijkl(i,j,:,7)*dxi(i,j,:,3,3) )*jacob(i,j,:)
        ! myyy
        ff(:,5)=( Eijkl(i,j,:,6)*dxi(i,j,:,3,1) +  &
                  Eijkl(i,j,:,8)*dxi(i,j,:,3,2) +  &
                  Eijkl(i,j,:,9)*dxi(i,j,:,3,3) )*jacob(i,j,:)
        ! myyz
        ff(:,6)=( Eijkl(i,j,:,7)*dxi(i,j,:,3,1) +  &
                  Eijkl(i,j,:,9)*dxi(i,j,:,3,2) +  &
                (-Eijkl(i,j,:,4)-Eijkl(i,j,:,8))*dxi(i,j,:,3,3) )*jacob(i,j,:)
        ! mxyz
        ff(:,7)=( Eijkl(i,j,:,5)*dxi(i,j,:,3,1) +  &
                  Eijkl(i,j,:,7)*dxi(i,j,:,3,2) +  &
                (-Eijkl(i,j,:,2)-Eijkl(i,j,:,6))*dxi(i,j,:,3,3) )*jacob(i,j,:)
        !
        do n=1,7
          df(:,n)=ddfc(ff(:,n),difschm,npdck,km,alfa_dif,dck,lfft=lfftk)
        enddo
        !
        do n=1,7
          qrhs_mom(i,j,ks:ke,8+n)=qrhs_mom(i,j,ks:ke,8+n)-df(ks:ke,n)
        enddo
        !
      enddo
      enddo
      deallocate(ff,df)
      !
    endif
    !
  end subroutine diffu_mijk
  !+-------------------------------------------------------------------+
  !| The end of the subroutine diffu_mijk.                             |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to solve the diffusion term of the rij.        |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 09-05-2023: Add heat flux term by J. Fang @ STFC DL               |
  !+-------------------------------------------------------------------+
  subroutine diffu_rij
    !
    use commvar,   only: npdci,npdcj,npdck,is,ie,js,je,ks,ke,difschm,lfftk,ndims
    use comsolver, only: grad,dci,dcj,dck,alfa_dif
    use commarray, only: miu,rho,dxi,jacob
    use commfunc,  only: ddfc
    use parallel, only : dataswap
    !
    real(8),allocatable :: coef(:,:,:)
    real(8),allocatable :: df(:,:),ff(:,:)
    integer :: i,j,k,n
    !
    if (.not. allocated(dRij)) then
       allocate(dRij(0:im,0:jm,0:km,1:3,1:6))
       dRij = 0.0d0
    end if
    !
    dRij(:,:,:,:,1)=grad(Rij(:,:,:,1))  ! Rxx 
    dRij(:,:,:,:,2)=grad(Rij(:,:,:,2))  ! Rxy 
    dRij(:,:,:,:,3)=grad(Rij(:,:,:,3))  ! Rxz
    dRij(:,:,:,:,4)=grad(Rij(:,:,:,4))  ! Ryy 
    dRij(:,:,:,:,5)=grad(Rij(:,:,:,5))  ! Ryz 
    dRij(:,:,:,:,6)=-dRij(:,:,:,:,1)-dRij(:,:,:,:,4) ! Rzz
    !
    if (.not. allocated(psi))  then
       allocate(psi(-hm:im+hm,-hm:jm+hm,-hm:km+hm,1:7))
       psi = 0.0d0
    end if
    allocate(coef(0:im,0:jm,0:km))
    !
    coef=miu(0:im,0:jm,0:km)/rho(0:im,0:jm,0:km)/Apsi1
    !
    ! psi_xxx
    psi(0:im,0:jm,0:km,1)=-num27d35*coef*( 3.d0*dRij(0:im,0:jm,0:km,1,1) - & ! dRxx/dx
                                          2.d0*(dRij(0:im,0:jm,0:km,2,2) + & ! dRxy/dy
                                                dRij(0:im,0:jm,0:km,3,3)) )  ! dRxz/dz
    ! psi_xxy
    psi(0:im,0:jm,0:km,2)=-num27d21*coef*(  dRij(0:im,0:jm,0:km,2,1) + & ! dRxx/dy 
                                      1.6d0*dRij(0:im,0:jm,0:km,1,2) - & ! dRxy/dx
                                     0.4d0*(dRij(0:im,0:jm,0:km,2,4) + & ! dRyy/dy
                                            dRij(0:im,0:jm,0:km,3,5)) )  ! dRyz/dz
    ! psi_xxz
    psi(0:im,0:jm,0:km,3)=-num27d21*coef*(  dRij(0:im,0:jm,0:km,3,1) + & ! dRxx/dz 
                                      1.6d0*dRij(0:im,0:jm,0:km,1,3) - & ! dRxz/dx
                                     0.4d0*(dRij(0:im,0:jm,0:km,2,5) + & ! dRyz/dy
                                            dRij(0:im,0:jm,0:km,3,6)) )  ! dRzz/dz
    ! psi_xyy
    psi(0:im,0:jm,0:km,4)=-num27d21*coef*(  dRij(0:im,0:jm,0:km,1,4) + & ! dRyy/dx 
                                      1.6d0*dRij(0:im,0:jm,0:km,2,2) - & ! dRxy/dy
                                     0.4d0*(dRij(0:im,0:jm,0:km,1,1) + & ! dRxx/dx
                                            dRij(0:im,0:jm,0:km,3,3)) )  ! dRxz/dz
    ! psi_yyy
    psi(0:im,0:jm,0:km,5)=-num27d35*coef*( 3.d0*dRij(0:im,0:jm,0:km,2,4) - & ! dRyy/dy
                                          2.d0*(dRij(0:im,0:jm,0:km,1,2) + & ! dRxy/dx
                                                dRij(0:im,0:jm,0:km,3,5)) )  ! dRyz/dz
    ! psi_yyz
    psi(0:im,0:jm,0:km,6)=-num27d21*coef*(  dRij(0:im,0:jm,0:km,3,4) + & ! dRyy/dz 
                                      1.6d0*dRij(0:im,0:jm,0:km,2,5) - & ! dRyz/dy
                                     0.4d0*(dRij(0:im,0:jm,0:km,1,3) + & ! dRxz/dx
                                            dRij(0:im,0:jm,0:km,3,6)) )  ! dRzz/dz
    ! psi_xyz
    psi(0:im,0:jm,0:km,7)=-num27d21*coef*(  dRij(0:im,0:jm,0:km,3,2) + & ! dRxy/dz 
                                            dRij(0:im,0:jm,0:km,2,3) + & ! dRxz/dy
                                            dRij(0:im,0:jm,0:km,1,5) )   ! dRyz/dx
    !
    call dataswap(psi)
    !
    allocate(ff(-hm:im+hm,1:5),df(0:im,1:5))
    do k=0,km
    do j=0,jm
      !
      ! Rxx
      ff(:,1)=( psi(:,j,k,1)*dxi(:,j,k,1,1) +  &
                psi(:,j,k,2)*dxi(:,j,k,1,2) +  &
                psi(:,j,k,3)*dxi(:,j,k,1,3) )*jacob(:,j,k)
      ! Rxy
      ff(:,2)=( psi(:,j,k,2)*dxi(:,j,k,1,1) +  &
                psi(:,j,k,4)*dxi(:,j,k,1,2) +  &
                psi(:,j,k,7)*dxi(:,j,k,1,3) )*jacob(:,j,k)
      ! Rxz
      ff(:,3)=( psi(:,j,k,3)*dxi(:,j,k,1,1) +  &
                psi(:,j,k,7)*dxi(:,j,k,1,2) +  &
              (-psi(:,j,k,1)-psi(:,j,k,4))*dxi(:,j,k,1,3) )*jacob(:,j,k)
      ! Ryy
      ff(:,4)=( psi(:,j,k,4)*dxi(:,j,k,1,1) +  &
                psi(:,j,k,5)*dxi(:,j,k,1,2) +  &
                psi(:,j,k,6)*dxi(:,j,k,1,3) )*jacob(:,j,k)
      ! Ryz
      ff(:,5)=( psi(:,j,k,7)*dxi(:,j,k,1,1) +  &
                psi(:,j,k,6)*dxi(:,j,k,1,2) +  &
              (-psi(:,j,k,2)-psi(:,j,k,5))*dxi(:,j,k,1,3) )*jacob(:,j,k)
      !
      do n=1,5
        df(:,n)=ddfc(ff(:,n),difschm,npdci,im,alfa_dif,dci)
      enddo
      !
      do n=1,5
        qrhs_mom(is:ie,j,k,15+n)=qrhs_mom(is:ie,j,k,15+n)-df(is:ie,n)
      enddo
      !
    enddo
    enddo
    deallocate(ff,df)
    !
    allocate(ff(-hm:jm+hm,1:5),df(0:jm,1:5))
    do k=0,km
    do i=0,im
      !
      ! Rxx
      ff(:,1)=( psi(i,:,k,1)*dxi(i,:,k,2,1) +  &
                psi(i,:,k,2)*dxi(i,:,k,2,2) +  &
                psi(i,:,k,3)*dxi(i,:,k,2,3) )*jacob(i,:,k)
      ! Rxy
      ff(:,2)=( psi(i,:,k,2)*dxi(i,:,k,2,1) +  &
                psi(i,:,k,4)*dxi(i,:,k,2,2) +  &
                psi(i,:,k,7)*dxi(i,:,k,2,3) )*jacob(i,:,k)
      ! Rxz
      ff(:,3)=( psi(i,:,k,3)*dxi(i,:,k,2,1) +  &
                psi(i,:,k,7)*dxi(i,:,k,2,2) +  &
              (-psi(i,:,k,1)-psi(i,:,k,4))*dxi(i,:,k,2,3) )*jacob(i,:,k)
      ! Ryy
      ff(:,4)=( psi(i,:,k,4)*dxi(i,:,k,2,1) +  &
                psi(i,:,k,5)*dxi(i,:,k,2,2) +  &
                psi(i,:,k,6)*dxi(i,:,k,2,3) )*jacob(i,:,k)
      ! Ryz
      ff(:,5)=( psi(i,:,k,7)*dxi(i,:,k,2,1) +  &
                psi(i,:,k,6)*dxi(i,:,k,2,2) +  &
              (-psi(i,:,k,2)-psi(i,:,k,5))*dxi(i,:,k,2,3) )*jacob(i,:,k)
      !
      do n=1,5
        df(:,n)=ddfc(ff(:,n),difschm,npdcj,jm,alfa_dif,dcj)
      enddo
      !
      do n=1,5
        qrhs_mom(i,js:je,k,15+n)=qrhs_mom(i,js:je,k,15+n)-df(js:je,n)
      enddo
      !
    enddo
    enddo
    deallocate(ff,df)
    !
    if(ndims==3) then
      !
      allocate(ff(-hm:km+hm,1:5),df(0:km,1:5))
      do j=0,jm
      do i=0,im
        !
        ! Rxx
        ff(:,1)=( psi(i,j,:,1)*dxi(i,j,:,3,1) +  &
                  psi(i,j,:,2)*dxi(i,j,:,3,2) +  &
                  psi(i,j,:,3)*dxi(i,j,:,3,3) )*jacob(i,j,:)
        ! Rxy
        ff(:,2)=( psi(i,j,:,2)*dxi(i,j,:,3,1) +  &
                  psi(i,j,:,4)*dxi(i,j,:,3,2) +  &
                  psi(i,j,:,7)*dxi(i,j,:,3,3) )*jacob(i,j,:)
        ! Rxz
        ff(:,3)=( psi(i,j,:,3)*dxi(i,j,:,3,1) +  &
                  psi(i,j,:,7)*dxi(i,j,:,3,2) +  &
                (-psi(i,j,:,1)-psi(i,j,:,4))*dxi(i,j,:,3,3) )*jacob(i,j,:)
        ! Ryy
        ff(:,4)=( psi(i,j,:,4)*dxi(i,j,:,3,1) +  &
                  psi(i,j,:,5)*dxi(i,j,:,3,2) +  &
                  psi(i,j,:,6)*dxi(i,j,:,3,3) )*jacob(i,j,:)
        ! Ryz
        ff(:,5)=( psi(i,j,:,7)*dxi(i,j,:,3,1) +  &
                  psi(i,j,:,6)*dxi(i,j,:,3,2) +  &
                (-psi(i,j,:,2)-psi(i,j,:,5))*dxi(i,j,:,3,3) )*jacob(i,j,:)
        !
        do n=1,5
          df(:,n)=ddfc(ff(:,n),difschm,npdck,km,alfa_dif,dck,lfft=lfftk)
        enddo
        !
        do n=1,5
          qrhs_mom(i,j,ks:ke,15+n)=qrhs_mom(i,j,ks:ke,15+n)-df(ks:ke,n)
        enddo
        !
      enddo
      enddo
      deallocate(ff,df)
      !
    endif
    !
    deallocate(coef)
    !
  end subroutine diffu_rij
  !+-------------------------------------------------------------------+
  !| The end of the subroutine diffu_rij.                              |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to solve the diffusion term of delta equation. |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 16-08-2023: Add heat flux term by J. Fang @ STFC DL               |
  !+-------------------------------------------------------------------+
  subroutine diffu_delta
    !
    use commvar,   only: npdci,npdcj,npdck,is,ie,js,je,ks,ke,difschm,lfftk,ndims
    use commarray, only: miu,rho,dxi,jacob
    use comsolver, only: dci,dcj,dck,alfa_dif
    use commfunc,  only: ddfc
    use parallel,  only : dataswap
    !
    real(8) :: var1
    integer :: i,j,k
    real(8),allocatable :: df(:),ff(:)
    !
    if (.not. allocated(omega)) then
       allocate(omega(-hm:im+hm,-hm:jm+hm,-hm:km+hm,1:3))
       omega = 0.0d0
    end if
    !
    do k=0,km
    do j=0,jm
    do i=0,im
      var1=miu(i,j,k)/rho(i,j,k)/Aomega1
      !
      omega(i,j,k,1)=-var1*(num7d3*ddelta(i,j,k,1) + &
                      4.d0*(dRij(i,j,k,1,1)+dRij(i,j,k,2,2)+dRij(i,j,k,3,3)))
      !
      omega(i,j,k,2)=-var1*(num7d3*ddelta(i,j,k,2) + &
                      4.d0*(dRij(i,j,k,1,2)+dRij(i,j,k,2,4)+dRij(i,j,k,3,5)))
      !
      omega(i,j,k,3)=-var1*(num7d3*ddelta(i,j,k,3) + &
                      4.d0*(dRij(i,j,k,1,3)+dRij(i,j,k,2,5)+dRij(i,j,k,3,6)))
    enddo
    enddo
    enddo
    !
    call dataswap(omega)
    !
    allocate(ff(-hm:im+hm),df(0:im))
    do k=0,km
    do j=0,jm
      !
      ff=( omega(:,j,k,1)*dxi(:,j,k,1,1) +  &
           omega(:,j,k,2)*dxi(:,j,k,1,2) +  &
           omega(:,j,k,3)*dxi(:,j,k,1,3) )*jacob(:,j,k)
      !
      df=ddfc(ff,difschm,npdci,im,alfa_dif,dci)
      !
      qrhs_mom(is:ie,j,k,21)=qrhs_mom(is:ie,j,k,21)-df(is:ie)
      !
    enddo
    enddo
    deallocate(ff,df)
    !
    allocate(ff(-hm:jm+hm),df(0:jm))
    do k=0,km
    do i=0,im
      !
      ff=( omega(i,:,k,1)*dxi(i,:,k,2,1) +  &
           omega(i,:,k,2)*dxi(i,:,k,2,2) +  &
           omega(i,:,k,3)*dxi(i,:,k,2,3) )*jacob(i,:,k)
      !
      df=ddfc(ff,difschm,npdcj,jm,alfa_dif,dcj)
      !
      qrhs_mom(i,js:je,k,21)=qrhs_mom(i,js:je,k,21)-df(js:je)
      !
    enddo
    enddo
    deallocate(ff,df)
    !
    if(ndims==3) then
      !
      allocate(ff(-hm:km+hm),df(0:km))
      do j=0,jm
      do i=0,im
        !
        ! Rxx
        ff=( omega(i,j,:,1)*dxi(i,j,:,3,1) +  &
             omega(i,j,:,2)*dxi(i,j,:,3,2) +  &
             omega(i,j,:,3)*dxi(i,j,:,3,3) )*jacob(i,j,:)
        !
        df=ddfc(ff,difschm,npdck,km,alfa_dif,dck,lfft=lfftk)
        !
        qrhs_mom(i,j,ks:ke,21)=qrhs_mom(i,j,ks:ke,21)-df(ks:ke)
        !
      enddo
      enddo
      deallocate(ff,df)
      !
    endif
    !
    !
  end subroutine diffu_delta
  !+-------------------------------------------------------------------+
  !| The end of the subroutine diffu_delta.                            |
  !+-------------------------------------------------------------------+

  !+-------------------------------------------------------------------+
  !| The following subroutines are to calculate the source term to the |
  !| r26 moment equations.                                             |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 18-04-2024: Created by J. Fang @ STFC Daresbury Laboratory        |
  !+-------------------------------------------------------------------+
  subroutine src_mijk_A
    !
    use constdef
    use commvar,   only : im,jm,km,hm,Reynolds,const2
    use commarray, only : tmp,sigma,miu,prs,jacob
    use comsolver,only : grad

    real(8),allocatable,dimension(:,:,:,:) :: Tsigma
    real(8),allocatable,dimension(:,:,:,:,:) :: dTsigma
    real(8) :: vasrc,Fijk(7)
    integer :: i,j,k,n
    !
    allocate(Tsigma(-hm:im+hm,-hm:jm+hm,-hm:km+hm,1:5), &
             dTsigma(0:im,0:jm,0:km,1:3,1:6) )
    !
    Tsigma(:,:,:,1)=tmp(:,:,:)*sigma(:,:,:,1)
    Tsigma(:,:,:,2)=tmp(:,:,:)*sigma(:,:,:,2)
    Tsigma(:,:,:,3)=tmp(:,:,:)*sigma(:,:,:,3)
    Tsigma(:,:,:,4)=tmp(:,:,:)*sigma(:,:,:,4)
    Tsigma(:,:,:,5)=tmp(:,:,:)*sigma(:,:,:,5)
    !
    dTsigma(:,:,:,:,1)=grad(Tsigma(:,:,:,1))
    dTsigma(:,:,:,:,2)=grad(Tsigma(:,:,:,2))
    dTsigma(:,:,:,:,3)=grad(Tsigma(:,:,:,3))
    dTsigma(:,:,:,:,4)=grad(Tsigma(:,:,:,4))
    dTsigma(:,:,:,:,5)=grad(Tsigma(:,:,:,5))
    dTsigma(:,:,:,:,6)=-dTsigma(:,:,:,:,1)-dTsigma(:,:,:,:,4)
    !
    do k=0,km
    do j=0,jm
    do i=0,im
      !
      Fijk(1)=  0.6d0 * dTsigma(i,j,k,1,1) -  &
                0.4d0 * (dTsigma(i,j,k,2,2)+dTsigma(i,j,k,3,3))                                                                                        
      Fijk(2)= num1d3 * (dTsigma(i,j,k,2,1) + 1.6d0 * dTsigma(i,j,k,1,2) - &
                         0.4d0 * (dTsigma(i,j,k,2,4)+dTsigma(i,j,k,3,5)) )                                                 
      Fijk(3)= num1d3 * (dTsigma(i,j,k,3,1) + 1.6d0 * dTsigma(i,j,k,1,3) - &
                         0.4d0 * (dTsigma(i,j,k,2,5)+dTsigma(i,j,k,3,6)) )                                                 
      Fijk(4)= num1d3 * (dTsigma(i,j,k,1,4)+  1.6d0 * dTsigma(i,j,k,2,2) - &
                         0.4d0 * (dTsigma(i,j,k,1,1)+dTsigma(i,j,k,3,3)) )
      Fijk(5)=  0.6d0 * dTsigma(i,j,k,2,4) - &
                0.4d0 * (dTsigma(i,j,k,1,2)+dTsigma(i,j,k,3,5))
      Fijk(6)= num1d3 * (dTsigma(i,j,k,3,4) + &
                         1.6d0 * dTsigma(i,j,k,2,5) -  &
                         0.4d0 * (dTsigma(i,j,k,1,3)+dTsigma(i,j,k,3,6)) )
      Fijk(7)= num1d3 * (dTsigma(i,j,k,3,2)+dTsigma(i,j,k,2,3)+dTsigma(i,j,k,1,5))
      !
      do n=1,7
        vasrc= -Am*prs(i,j,k)/miu(i,j,k)*mijk(i,j,k,n)-3.d0/const2*Fijk(n)
        !
        qrhs_mom(i,j,k,8+n)=qrhs_mom(i,j,k,8+n) + vasrc*jacob(i,j,k)
      enddo
      !
    enddo
    enddo
    enddo
    !
    !
  end subroutine src_mijk_A
  !
  subroutine src_delta
    !
    use commvar,   only : const2
    use commarray, only : rho,tmp,prs,miu,jacob,sigma,qflux,dvel,dtmp

    real(8) :: var1,var2,var3
    real(8) :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz, &
               dpdx,dpdy,dpdz,dtdx,dtdy,dtdz,qx,qy,qz,       &
               sigmaxx,sigmaxy,sigmaxz,sigmayy,sigmayz,sigmazz
    real(8) :: yx,yy,yz,eksxx,eksxy,eksxz,eksyy,eksyz,ekszz
    integer :: i,j,k
    !
    do k=0,km
    do j=0,jm
    do i=0,im
      !
      var1=-Adelta1*prs(i,j,k)*delta(i,j,k)/miu(i,j,k)
      var2=-8.d0/const2*tmp(i,j,k)*(dqx(i,j,k,1)+dqy(i,j,k,2)+dqz(i,j,k,3))
      !
      qx=qflux(i,j,k,1)
      qy=qflux(i,j,k,2)
      qz=qflux(i,j,k,3)
      !
      dudx=dvel(i,j,k,1,1)
      dudy=dvel(i,j,k,1,2)
      dudz=dvel(i,j,k,1,3)
      dvdx=dvel(i,j,k,2,1)
      dvdy=dvel(i,j,k,2,2)
      dvdz=dvel(i,j,k,2,3)
      dwdx=dvel(i,j,k,3,1)
      dwdy=dvel(i,j,k,3,2)
      dwdz=dvel(i,j,k,3,3)
      !
      dpdx=dprs(i,j,k,1)
      dpdy=dprs(i,j,k,2)
      dpdz=dprs(i,j,k,3)
      !
      dtdx=dtmp(i,j,k,1)
      dtdy=dtmp(i,j,k,2)
      dtdz=dtmp(i,j,k,3)
      !
      sigmaxx=sigma(i,j,k,1)
      sigmaxy=sigma(i,j,k,2)
      sigmaxz=sigma(i,j,k,3)
      sigmayy=sigma(i,j,k,4)
      sigmayz=sigma(i,j,k,5)
      sigmazz=-sigmaxx-sigmayy
      !
      yx=8.d0*dpdx/prs(i,j,k)-28.d0*dtdx/tmp(i,j,k)
      yy=8.d0*dpdy/prs(i,j,k)-28.d0*dtdy/tmp(i,j,k)
      yz=8.d0*dpdz/prs(i,j,k)-28.d0*dtdz/tmp(i,j,k)
      !
      eksxx=2.d0/const2*tmp(i,j,k)*sigmaxx+rij(i,j,k,1)
      eksxy=2.d0/const2*tmp(i,j,k)*sigmaxy+rij(i,j,k,2)
      eksxz=2.d0/const2*tmp(i,j,k)*sigmaxz+rij(i,j,k,3)
      eksyy=2.d0/const2*tmp(i,j,k)*sigmayy+rij(i,j,k,4)
      eksyz=2.d0/const2*tmp(i,j,k)*sigmayz+rij(i,j,k,5)
      ekszz=2.d0/const2*tmp(i,j,k)*sigmazz-rij(i,j,k,1)-rij(i,j,k,4)
      !
      var3=-num4d3*delta(i,j,k)*(dudx+dvdy+dwdz) + tmp(i,j,k)/const2*(qx*yx+qy*yy+qz*yz) - &
            4.d0*(eksxx*dudx+eksyy*dvdy+ekszz*dwdz+eksxy*(dudy+dvdx)+eksxz*(dudz+dwdx)+eksyz*(dvdz+dwdy)) + &
            8.d0/rho(i,j,k)*( qx*(dsigmaxx(i,j,k,1)+dsigmaxy(i,j,k,2)+dsigmaxz(i,j,k,3)) + &
                              qy*(dsigmaxy(i,j,k,1)+dsigmayy(i,j,k,2)+dsigmayz(i,j,k,3)) + &
                              qz*(dsigmaxz(i,j,k,1)+dsigmayz(i,j,k,2)+dsigmazz(i,j,k,3)) ) - &
            Adelta2*prs(i,j,k)/(rho(i,j,k)*miu(i,j,k))*( sigmaxx**2 + sigmayy**2 + sigmazz**2 + &
                                                  2.d0*( sigmaxy**2 + sigmaxz**2 + sigmayz**2) )

      !
      qrhs_mom(i,j,k,21)=qrhs_mom(i,j,k,21)+(var1+var2+var3)*jacob(i,j,k)
      !
    enddo
    enddo
    enddo
    !
  end subroutine src_delta
  !
  subroutine src_rij_A
    !
    use commvar,   only : const2
    use commarray, only : tmp,prs,dvel,miu,jacob

    integer :: i,j,k,n
    real(8) :: dqxdx,dqxdy,dqxdz,dqydx,dqydy,      &
               dqydz,dqzdx,dqzdy,dqzdz
    real(8) :: var1,var2,srcterm(5)
    !
    do k=0,km
    do j=0,jm
    do i=0,im
      !
      dqxdx=dqx(i,j,k,1)
      dqxdy=dqx(i,j,k,2)
      dqxdz=dqx(i,j,k,3)
      dqydx=dqy(i,j,k,1)
      dqydy=dqy(i,j,k,2)
      dqydz=dqy(i,j,k,3)
      dqzdx=dqz(i,j,k,1)
      dqzdy=dqz(i,j,k,2)
      dqzdz=dqz(i,j,k,3)
      !
      var1=AR1*prs(i,j,k)/miu(i,j,k)
      var2=tmp(i,j,k)/const2
      !
      srcterm(1)=-var1*rij(i,j,k,1)-num28d15*var2*(2.d0*dqxdx-dqydy-dqzdz)
      srcterm(2)=-var1*rij(i,j,k,2)- num14d5*var2*(dqydx+dqxdy)
      srcterm(3)=-var1*rij(i,j,k,3)- num14d5*var2*(dqzdx+dqxdz)
      srcterm(4)=-var1*rij(i,j,k,4)-num28d15*var2*(2.d0*dqydy-dqxdx-dqzdz)
      srcterm(5)=-var1*rij(i,j,k,5)- num14d5*var2*(dqzdy+dqydz)
      !
      do n=1,5
        !
        qrhs_mom(i,j,k,15+n)=qrhs_mom(i,j,k,15+n) + srcterm(n)*jacob(i,j,k)
        !
      enddo
      !
    enddo
    enddo
    enddo
    !
  end subroutine src_rij_A
  !
  subroutine src_rij_B
    !   ---------------------------
    !
    use constdef
    use commvar,   only : im,jm,km,const2
    use commarray, only : rho,dvel,sigma,qflux,tmp,prs,dtmp,jacob,miu
    use comsolver, only : grad
    !
    integer :: i,j,k,n
    real(8) :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz,             &
               sxx,sxy,sxz,syy,syz,szz, dqxdx,dqxdy,dqxdz,dqydx,dqydy,   &
               dqydz,dqzdx,dqzdy,dqzdz,dtdx,dtdy,dtdz,co1,               &
               dpx,dpy,dpz,qx,qy,qz,dpdx,dpdy,dpdz,rrho
    real(8) :: sigmaxx,sigmaxy,sigmaxz,sigmayy,sigmayz,sigmazz
    real(8) :: mxxx,mxxy,mxyy,myyy,mxxz,mxyz,myzz,mzzz,myyz,mxzz
    real(8) :: rxx, ryy, rxy, rxz, ryz, rzz
    real(8) :: a,b,c,d,e,srcterm(5)
    real(8) :: t12, t13, t23, t12a, t13a, t23a, term1, term2, term3, term4, term5, term6

    real(8) :: axx, ayy, axy, axz, ayz, bxx, byy, bxy, bxz, byz, cxx, cyy, cxy, cxz, cyz, &
               dxx, dyy, dxy, dxz, dyz, exx, eyy, exy, exz, eyz, fxx, fyy, fxy, fxz, fyz, &
               gxx, gyy, gxy, gxz, gyz, hxx, hyy, hxy, hxz, hyz, kxx, kyy, kxy, kxz, kyz

    real(8) :: qxxxx, qyyyy, qxxxy, qxxxz, qxxyy, qxyyy, qxxyz, qxyyz, qyyyz, &
               qxzzz, qyzzz, qxxzz, qyyzz, qxyzz, qzzzz

    real(8) :: dmxxxdx, dmxxydy, dmxxzdz, dmxyydx, dmyyydy, dmyyzdz, dmxxydx,  & 
               dmxyydy, dmxyzdz, dmxxzdx, dmxyzdy, dmxzzdz, dmxyzdx, dmyyzdy,  &
               dmyzzdz 

    real(8) :: doxdx, doxdy, doxdz, doydx, doydy, doydz, dozdx, dozdy, dozdz 

    if(.not. allocated(domega)) allocate(domega(0:im,0:jm,0:km,1:3,1:3))
    domega(:,:,:,:,1)=grad(omega(:,:,:,1))  
    domega(:,:,:,:,2)=grad(omega(:,:,:,2))  
    domega(:,:,:,:,3)=grad(omega(:,:,:,3))  
    !
    srcterm=0.d0
    !
    do k=0,km
    do j=0,jm
    do i=0,im

      dudx = dvel(i,j,k,1,1)
      dudy = dvel(i,j,k,1,2)
      dudz = dvel(i,j,k,1,3)
      dvdx = dvel(i,j,k,2,1)
      dvdy = dvel(i,j,k,2,2)
      dvdz = dvel(i,j,k,2,3)
      dwdx = dvel(i,j,k,3,1)
      dwdy = dvel(i,j,k,3,2)
      dwdz = dvel(i,j,k,3,3)
      !
      dtdx=dtmp(i,j,k,1)
      dtdy=dtmp(i,j,k,2)
      dtdz=dtmp(i,j,k,3)      
      !
      sigmaxx = sigma(i,j,k,1)
      sigmaxy = sigma(i,j,k,2)
      sigmaxz = sigma(i,j,k,3)
      sigmayy = sigma(i,j,k,4)
      sigmayz = sigma(i,j,k,5)
      sigmazz = -sigmaxx-sigmayy
      !
      qx = qflux(i,j,k,1)
      qy = qflux(i,j,k,2)
      qz = qflux(i,j,k,3)
      !
      dqxdx = dqx(i,j,k,1)
      dqxdy = dqx(i,j,k,2)
      dqxdz = dqx(i,j,k,3)
      dqydx = dqy(i,j,k,1)
      dqydy = dqy(i,j,k,2) 
      dqydz = dqy(i,j,k,3)
      dqzdx = dqz(i,j,k,1)
      dqzdy = dqz(i,j,k,2)
      dqzdz = dqz(i,j,k,3)
      !
      dpdx = dprs(i,j,k,1)
      dpdy = dprs(i,j,k,2)
      dpdz = dprs(i,j,k,3)
      !
      mxxx = mijk(i,j,k,1)
      mxxy = mijk(i,j,k,2)
      mxxz = mijk(i,j,k,3)
      mxyy = mijk(i,j,k,4)
      myyy = mijk(i,j,k,5)
      myyz = mijk(i,j,k,6)
      mxyz = mijk(i,j,k,7)     
      mxzz = - mxxx - mxyy
      myzz = - mxxy - myyy
      mzzz = - mxxz - myyz      
      
      !
      rxx = rij(i,j,k,1)
      rxy = rij(i,j,k,2)
      rxz = rij(i,j,k,3)
      ryy = rij(i,j,k,4)
      ryz = rij(i,j,k,5)
      rzz = -rxx-ryy

      a = dudx + dvdy + dwdz
      b = tmp(i,j,k)/const2

      c = num8d3*b
      
      axx = a*(c*sigmaxx - num2d7*rxx)
      ayy = a*(c*sigmayy - num2d7*ryy)
      axy = a*(c*sigmaxy - num2d7*rxy)
      axz = a*(c*sigmaxz - num2d7*rxz)
      ayz = a*(c*sigmayz - num2d7*ryz)
      
      c = 7.0d0*b
      sxx = -num4d7*(c*sigmaxx + rxx)
      syy = -num4d7*(c*sigmayy + ryy)
      sxy = -num4d7*(c*sigmaxy + rxy)
      sxz = -num4d7*(c*sigmaxz + rxz)
      syz = -num4d7*(c*sigmayz + ryz)
      szz = -sxx - syy
      
      t12 = dudy + dvdx
      t13 = dudz + dwdx
      t23 = dvdz + dwdy

      t12a = dudx + dvdy
      t13a = dudx + dwdz
      t23a = dvdy + dwdz
      
      bxx = (2.0d0*rxx*dudx + rxy*(2.0d0*dudy - dvdx)  & 
          + rxz*(2.0d0*dudz - dwdx)                    &
          - (ryy*dvdy + rzz*dwdz + ryz*t23) )/3.0d0         
          
      byy = (2.0d0*ryy*dvdy + rxy*(2.0d0*dvdx - dudy)  &
          + ryz*(2.0d0*dvdz - dwdy)                    &
          - (rxx*dudx + rzz*dwdz + rxz*t13) )/3.0d0 
          
      bxy = rxy*t12a + rxx*dvdx + ryy*dudy + rxz*dvdz + ryz*dudz
          
      bxz = rxz*t13a + rxx*dwdx + rzz*dudz + rxy*dwdy + ryz*dudy      
      
      byz = ryz*t23a + ryy*dwdy + rzz*dvdz + rxy*dwdx + rxz*dvdx
          
      bxx = (4.0d0*sxx*dudx + sxy*t12  + sxz*t13         &
          - 2.0d0*(syy*dvdy + szz*dwdz + syz*t23))/3.0d0 &
          - 2.0d0*bxx

      byy = (4.0d0*syy*dvdy + sxy*t12  + syz*t23         &
          - 2.0d0*(sxx*dudx + szz*dwdz + sxz*t13))/3.0d0 &
          - 2.0d0*byy

      bxy = sxy*t12a + 0.5d0*( (sxx + syy)*t12 + sxz*t23 + syz*t13) &
          - bxy

      bxz = sxz*t13a + 0.5d0*( (sxx + szz)*t13 + sxy*t23 + syz*t12) &
          - bxz

      byz = syz*t23a + 0.5d0*( (syy + szz)*t23 + sxy*t13 + sxz*t12) &
          - byz

     
      d = tmp(i,j,k)/prs(i,j,k)
    
      a = (d*dpdx - 2.0d0*dtdx)
      b = (d*dpdy - 2.0d0*dtdy)
      c = (d*dpdz - 2.0d0*dtdz)
      
      d =  num28d15/const2

      dxx = d*(2.0d0*qx*a - qy*b - qz*c)
      dyy = d*(2.0d0*qy*b - qx*a - qz*c)
      dxy = 1.5d0*d*(qx*b + qy*a)
      dxz = 1.5d0*d*(qx*c + qz*a)
      dyz = 1.5d0*d*(qy*c + qz*b)                  

      a = dsigmaxx(i,j,k,1) + dsigmaxy(i,j,k,2) + dsigmaxz(i,j,k,3)
      b = dsigmaxy(i,j,k,1) + dsigmayy(i,j,k,2) + dsigmayz(i,j,k,3)
      c = dsigmaxz(i,j,k,1) + dsigmayz(i,j,k,2) + dsigmazz(i,j,k,3)


      d =  num28d15/rho(i,j,k)
       
      cxx = d*(2.0d0*qx*a - qy*b - qz*c)
      cyy = d*(2.0d0*qy*b - qx*a - qz*c)
      cxy = 1.5d0*d*(qx*b + qy*a)
      cxz = 1.5d0*d*(qx*c + qz*a)
      cyz = 1.5d0*d*(qy*c + qz*b)
      
      d = 2.0d0/rho(i,j,k)
      e = 9.0d0/const2
      
      a = d*(a + dpdx) - e*dtdx   
      b = d*(b + dpdy) - e*dtdy    
      c = d*(c + dpdz) - e*dtdz 

      exx = mxxx*a + mxxy*b + mxxz*c
      eyy = mxyy*a + myyy*b + myyz*c
      exy = mxxy*a + mxyy*b + mxyz*c
      exz = mxxz*a + mxyz*b + mxzz*c
      eyz = mxyz*a + myyz*b + myzz*c
         
      dmxxxdx = dmijk(i,j,k,1,1)
      dmxxydy = dmijk(i,j,k,2,2)
      dmxxzdz = dmijk(i,j,k,3,3)
      
      dmxyydx = dmijk(i,j,k,1,4)
      dmyyydy = dmijk(i,j,k,2,5)
      dmyyzdz = dmijk(i,j,k,3,6)
      
      dmxxydx = dmijk(i,j,k,1,2)
      dmxyydy = dmijk(i,j,k,2,4)
      dmxyzdz = dmijk(i,j,k,3,7)
      
      dmxxzdx = dmijk(i,j,k,1,3)
      dmxyzdy = dmijk(i,j,k,2,7)
      dmxzzdz = -dmijk(i,j,k,3,1) - dmijk(i,j,k,3,4)
      
      dmxyzdx = dmijk(i,j,k,1,7)
      dmyyzdy = dmijk(i,j,k,2,6)
      dmyzzdz = -dmijk(i,j,k,3,2) - dmijk(i,j,k,3,5)
      
      d = 2.0d0*Tmp(i,j,k)/const2
      
      exx = exx - d*(dmxxxdx + dmxxydy + dmxxzdz)
      eyy = eyy - d*(dmxyydx + dmyyydy + dmyyzdz)
      exy = exy - d*(dmxxydx + dmxyydy + dmxyzdz)
      exz = exz - d*(dmxxzdx + dmxyzdy + dmxzzdz)
      eyz = eyz - d*(dmxyzdx + dmyyzdy + dmyzzdz)
      
      a = (sigmaxx*dudx + sigmayy*dvdy + sigmazz*dwdz       &
         + sigmaxy*t12 + sigmaxz*t13  + sigmayz*t23         &
         + dqxdx + dqydy + dqzdz)*num14d3/rho(i,j,k)

      fxx = sigmaxx*a ; fyy = sigmayy*a ; fxy = sigmaxy*a
      fxz = sigmaxz*a ; fyz = sigmayz*a
   
      if (.false.) then
      qxxxx = Eijkl(i,j,k,1)
      qyyyy = Eijkl(i,j,k,8)
      qxxxy = Eijkl(i,j,k,2)
      qxxxz = Eijkl(i,j,k,3)
      qxxyy = Eijkl(i,j,k,4)
      qxyyy = Eijkl(i,j,k,6)
      qxxyz = Eijkl(i,j,k,5)
      qxyyz = Eijkl(i,j,k,7)
      qyyyz = Eijkl(i,j,k,9)

      qxzzz = - qxxxz - qxyyz
      qyzzz = - qxxyz - qyyyz
      qxxzz = - qxxxx - qxxyy
      qyyzz = - qxxyy - qyyyy
      qxyzz = - qxxxy - qxyyy
      qzzzz = - qxxzz - qyyzz
                
      gxx = -2.0d0*(qxxxx*dudx + qxxyy*dvdy + qxxzz*dwdz    &
          + qxxxy*t12 + qxxxz*t13 + qxxyz*t23)

      gyy = -2.0d0*(qxxyy*dudx + qyyyy*dvdy + qyyzz*dwdz    &
          + qxyyy*t12 + qxyyz*t13 + qyyyz*t23)

      gxy = -2.0d0*(qxxxy*dudx + qxyyy*dvdy + qxyzz*dwdz    &
          + qxxyy*t12 + qxxyz*t13 + qxyyz*t23)

      gxz = -2.0d0*(qxxxz*dudx + qxyyz*dvdy + qxzzz*dwdz    &
          + qxxyz*t12 + qxxzz*t13 + qxyzz*t23)

      gyz = -2.0d0*(qxxyz*dudx + qyyyz*dvdy + qyzzz*dwdz    &
          + qxyyz*t12 + qxyzz*t13 + qyyzz*t23)
           

      doxdx = domega(i,j,k,1,1)
      doxdy = domega(i,j,k,2,1)
      doxdz = domega(i,j,k,3,1)
      doydx = domega(i,j,k,1,2)
      doydy = domega(i,j,k,2,2)
      doydz = domega(i,j,k,3,2)
      dozdx = domega(i,j,k,1,3)
      dozdy = domega(i,j,k,2,3)
      dozdz = domega(i,j,k,3,3)


      a = num7d3*delta(i,j,k)
      hxx = -num2d15*(a*(2.0d0*dudx - dvdy - dwdz) + 2.0d0*doxdx - doydy - dozdz) 
      hyy = -num2d15*(a*(2.0d0*dvdy - dudx - dwdz) + 2.0d0*doydy - doxdx - dozdz)  
      hxy = -0.2d0*(a*t12 + doydx + doxdy)
      hxz = -0.2d0*(a*t13 + dozdx + doxdz)
      hyz = -0.2d0*(a*t23 + dozdy + doydz)
      
      term1 = sigmaxx*sigmaxx ; term2 = sigmayy*sigmayy ; term3 = sigmazz*sigmazz
      term4 = sigmaxy*sigmaxy ; term5 = sigmaxz*sigmaxz ; term6 = sigmayz*sigmayz 

      a = -AR2*prs(i,j,k)/(rho(i,j,k)*miu(i,j,k))

      kxx = a*(2.0d0*(term1-term6)+term4+term5-term2-term3)/3.0d0
      kyy = a*(2.0d0*(term2-term5)+term4+term6-term1-term3)/3.0d0
      kxy = a*( (sigmaxx + sigmayy)*sigmaxy + sigmaxz*sigmayz)
      kxz = a*( (sigmaxx + sigmazz)*sigmaxz + sigmaxy*sigmayz)
      kyz = a*( (sigmayy + sigmazz)*sigmayz + sigmaxy*sigmaxz)
      end if

! e
!     Rxx
      srcterm(1) = axx + bxx + cxx + dxx + exx + fxx !+ gxx + hxx + kxx
      
!     Rxy
      srcterm(2) = axy + bxy + cxy + dxy + exy + fxy !+ gxy + hxy + kxy

!     Rxz
      srcterm(3) = axz + bxz + cxz + dxz + exz + fxz !+ gxz + hxz + kxz
      
!     Ryy
      srcterm(4) = ayy + byy + cyy + dyy + eyy + fyy !+ gyy + hyy + kyy

!     Ryz
      srcterm(5) = ayz + byz + cyz + dyz + eyz + fyz !+ gyz + hyz + kyz
                    
      do n=1,5
        !
        qrhs_mom(i,j,k,15+n)=qrhs_mom(i,j,k,15+n) + srcterm(n)*jacob(i,j,k)
        !
      enddo

 
      end do ; end do; end do
      return
      !
  end subroutine src_rij_B

  subroutine src_mijk_B
    !
    use constdef
    use commvar,   only : im,jm,km
    use commarray, only : rho,dvel,sigma,qflux,jacob
    use comsolver,only : grad
    !
    integer :: i,j,k,n
    real(8) :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz,       &
               sxx,sxy,sxz,syy,syz,dqxdx,dqxdy,dqxdz,dqydx,dqydy,      &
               dqydz,dqzdx,dqzdy,dqzdz,dtdx,dtdy,dtdz,                 &
               sigmaxx,sigmaxy,sigmaxz,sigmayy,sigmayz,sigmazz,co1,    &
               dpx,dpy,dpz,qx,qy,qz,dpdx,dpdy,dpdz,rrho
    real(8) :: mxxx,mxxy,mxyy,myyy,mxxz,mxyz,myzz,mzzz,myyz,mxzz
    real(8) :: a,b,c,d,e,srcterm(7)

    real(8) :: terma, termb, termc, t12, t13, t23
    !
    srcterm=0.d0
    !
    do k=0,km
    do j=0,jm
    do i=0,im
      !
      rrho=1.d0/rho(i,j,k)
      !
      dudx=dvel(i,j,k,1,1)
      dudy=dvel(i,j,k,1,2)
      dudz=dvel(i,j,k,1,3)
      dvdx=dvel(i,j,k,2,1)
      dvdy=dvel(i,j,k,2,2)
      dvdz=dvel(i,j,k,2,3)
      dwdx=dvel(i,j,k,3,1)
      dwdy=dvel(i,j,k,3,2)
      dwdz=dvel(i,j,k,3,3)
      !
      sigmaxx=sigma(i,j,k,1)
      sigmaxy=sigma(i,j,k,2)
      sigmaxz=sigma(i,j,k,3)
      sigmayy=sigma(i,j,k,4)
      sigmayz=sigma(i,j,k,5)
      sigmazz=-sigmaxx-sigmayy
      !
      qx=qflux(i,j,k,1)
      qy=qflux(i,j,k,2)
      qz=qflux(i,j,k,3)
      !
      dpdx=dprs(i,j,k,1)
      dpdy=dprs(i,j,k,2)
      dpdz=dprs(i,j,k,3)
      !
      mxxx=mijk(i,j,k,1)
      mxxy=mijk(i,j,k,2)
      mxxz=mijk(i,j,k,3)
      mxyy=mijk(i,j,k,4)
      myyy=mijk(i,j,k,5) 
      myyz=mijk(i,j,k,6)
      mxyz=mijk(i,j,k,7)
      mxzz=-mxxx-mxyy
      myzz=-mxxy-myyy
      mzzz=-mxxz-myyz


      terma = (dsigmaxx(i,j,k,1)+dsigmaxy(i,j,k,2)+dsigmaxz(i,j,k,3)) + dpdx 
      termb = (dsigmaxy(i,j,k,1)+dsigmayy(i,j,k,2)+dsigmayz(i,j,k,3)) + dpdy
      termc = (dsigmaxz(i,j,k,1)+dsigmayz(i,j,k,2)+dsigmazz(i,j,k,3)) + dpdz

      t12 = dudy+dvdx
      t13 = dudz+dwdx
      t23 = dwdy+dvdz
      !
      a = num9d35*dRij(i,j,k,1,1) - num6d35*(dRij(i,j,k,2,2) + dRij(i,j,k,3,3))
      b = 1.8d0*sigmaxx*terma - 1.2d0*(sigmaxy*termb + sigmaxz*termc)
      d = 0.48d0*(qx*(2.d0*dudx - t23) - qy*t12 - qz*t13)
      e = 1.8d0*mxxx*dudx + mxxy*(1.8d0*dudy - 1.2d0*dvdx)                                  & 
        + mxxz*(1.8d0*dudz - 1.2d0*dwdx) - 1.2d0*(mxyy*dvdy + mxyz*t23 + mxzz*dwdz)
      !
      srcterm(1) = - a + rrho*b - d - e
      !
      a = num1d7*dRij(i,j,k,2,1) + num8d35*dRij(i,j,k,1,2)                                  &
        - num2d35*(dRij(i,j,k,2,4) + dRij(i,j,k,3,5))
      b = (sigmaxx  - 0.4d0*sigmayy)*termb + 1.6d0*sigmaxy*terma - 0.4d0*sigmayz*termc
      d = 0.64d0*qx*t12 + 0.16d0*(qy*(4.d0*dudx - 3.d0*dvdy - dwdz) - qz*t23)
      e = mxxx*dvdx +  (1.6d0*dudx+dvdy)*mxxy + mxxz*dvdz + 0.4d0*(mxyy*(4.d0*dudy - dvdx)  & 
        + mxyz*(4.d0*dudz - dwdx) - myyy*dvdy - myyz*t23 - myzz*dwdz)
      !
      srcterm(2) = - a + rrho*b - d - e
      !
      a = num1d7*dRij(i,j,k,3,1) + num8d35*dRij(i,j,k,1,3)                                  &
        - num2d35*(dRij(i,j,k,2,5) + dRij(i,j,k,3,6))  
      b = (sigmaxx - 0.4d0*sigmazz)*termc  + 1.6d0*sigmaxz*terma - 0.4d0*sigmayz*termb  
      d = 0.64d0*qx*t13 - 0.16d0*qy*t23 + 0.16d0*qz*(4.d0*dudx-dvdy-3.d0*dwdz)
      e = mxxx*dwdx + mxxz*(1.6d0*dudx + dwdz) + mxxy*dwdy + 0.4d0*(mxyz*(4.d0*dudy - dvdx)  &
        + mxzz*(4.d0*dudz - dwdx) - myyz*dvdy - myzz*t23 - mzzz*dwdz)

      srcterm(3) = - a + rrho*b - d - e
      !
      a = num1d7*dRij(i,j,k,1,4) + num8d35*dRij(i,j,k,2,2)                                  &
        - num2d35*(dRij(i,j,k,1,1) + dRij(i,j,k,3,3)) 
      b = (sigmayy - 0.4d0*sigmaxx)*terma + 1.6d0*sigmaxy*termb - 0.4d0*sigmaxz*termc 
      d = 0.16d0*qx*(4.d0*dvdy - 3.d0*dudx - dwdz)+ 0.64d0*qy*t12 - 0.16d0*qz*(dwdx - dudz)
      e = mxyy*(1.6d0*dvdy + dudx) + myyy*dudy + myyz*dudz                                  &
        + 0.4d0*(mxxy*(4.d0*dvdx - dudy) + mxyz*(4.d0*dvdz - dwdy) - mxxx*dudx              &
        - mxxz*t13 - mxzz*dwdz) 

      srcterm(4) = - a + rrho*b - d - e
      !
      a = num9d35*dRij(i,j,k,2,4) - num6d35*(dRij(i,j,k,1,2) + dRij(i,j,k,3,5))
      b = - 1.2d0*sigmaxy*terma + 1.8d0*sigmayy*termb - 1.2d0*sigmayz*termc
      d = 0.48d0*(- qx*t12 + qy*(2.d0*dvdy - dudx - dwdz) - qz*t23 )
      e = 1.8d0*myyy*dvdy + mxyy*(1.8d0*dvdx - 1.2d0*dudy) + myyz*(1.8d0*dvdz - 1.2d0*dwdy) &
        - 1.2d0*(mxxy*dudx + mxyz*t13 + myzz*dwdy)

      srcterm(5) = - a + rrho*b - d - e
      !
      a = num1d7*dRij(i,j,k,3,4) + num8d35*dRij(i,j,k,2,5)                                  &
        - num2d35*(dRij(i,j,k,1,3) + dRij(i,j,k,3,6)) 
      b = - 0.4d0*sigmaxz*terma + (sigmayy - 0.4d0*sigmazz)*termc + 1.6d0*sigmayz*termb 
      d = - 0.16d0*qx*t13 + 0.64d0*qy*t23 + 0.16d0*qz*(4.d0*dvdy - dudx - 3.d0*dwdz)
      e = myyy*dwdy + myyz*(1.6d0*dvdy + dwdz) + mxyy*dwdx + 0.4d0*(mxyz*(4.d0*dvdx - dudy)  &
        - mxxz*dudx - mxzz *t13 + myzz*(4.d0*dvdz - dwdy) - mzzz*dwdz)

      srcterm(6) = - a + rrho*b - d - e
      !
      a = num1d7*(dRij(i,j,k,1,5) + dRij(i,j,k,2,3) + dRij(i,j,k,3,2))
      b = sigmayz*terma + sigmaxy*termc + sigmaxz*termb
      d = 0.4d0*(qx*t23 + qy*t13 + qz * t12)
      e = mxxy*dwdx + mxxz*dvdx + mxyz*(dudx+dvdy+dwdz) + myzz*dudz + myyz*dudy + mxzz*dvdz  &
        + mxyy*dwdy

      srcterm(7) = - a + rrho*b - d - e
      !
      do n=1,7
        !
        qrhs_mom(i,j,k,8+n) = qrhs_mom(i,j,k,8+n) + srcterm(n)*jacob(i,j,k)
        !
      enddo
      !
    enddo
    enddo
    enddo
    !
  end subroutine src_mijk_B
  !+-------------------------------------------------------------------+
  !| The end of the subroutines src_***.                               |
  !+-------------------------------------------------------------------+

  !+-------------------------------------------------------------------+
  !| This subroutine is to calculate the source term from the r13      |
  !| moment equations.                                                 |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 01-11-2022: Created by J. Fang @ STFC Daresbury Laboratory        |
  !+-------------------------------------------------------------------+
  subroutine src_sigma_q
    !
    use commvar,   only : im,jm,km,reynolds,Mach,gamma,num_species,    &
                          num_modequ,const2,const5
    use commarray, only : rho,miu,sigma,qflux,prs,tmp,jacob,dvel,dtmp
    use comsolver, only : grad
    use parallel,  only : dataswap

    implicit none
    !
    ! local data
    real(8) :: miup,dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz,       &
               sxx,sxy,sxz,syy,syz,dqxdx,dqxdy,dqxdz,dqydx,dqydy,      &
               dqydz,dqzdx,dqzdy,dqzdz,dtdx,dtdy,dtdz,                 &
               sigmaxx,sigmaxy,sigmaxz,sigmayy,sigmayz,sigmazz,co1,    &
               pxx,pyy,pzz,pxx0,pyy0,pzz0,qx,qy,qz,mxzz,myzz,mzzz, pres
    !
    real(8) :: srcterm(8),srcmax,srcmin
    integer :: i,j,k
    logical,save :: firstcall=.true.
    !
    call dataswap(prs)
    !
    dprs=grad(prs)
    !
    srcmax=0.d0
    srcmin=1.d10
    do k=0,km
    do j=0,jm
    do i=0,im
      !
      miup=miu(i,j,k)
      !
      dudx=dvel(i,j,k,1,1)
      dudy=dvel(i,j,k,1,2)
      dudz=dvel(i,j,k,1,3)
      dvdx=dvel(i,j,k,2,1)
      dvdy=dvel(i,j,k,2,2)
      dvdz=dvel(i,j,k,2,3)
      dwdx=dvel(i,j,k,3,1)
      dwdy=dvel(i,j,k,3,2)
      dwdz=dvel(i,j,k,3,3)
      !
      dqxdx=dqx(i,j,k,1)
      dqxdy=dqx(i,j,k,2)
      dqxdz=dqx(i,j,k,3)
      dqydx=dqy(i,j,k,1)
      dqydy=dqy(i,j,k,2)
      dqydz=dqy(i,j,k,3)
      dqzdx=dqz(i,j,k,1)
      dqzdy=dqz(i,j,k,2)
      dqzdz=dqz(i,j,k,3)
      !
      sigmaxx=sigma(i,j,k,1)
      sigmaxy=sigma(i,j,k,2)
      sigmaxz=sigma(i,j,k,3)
      sigmayy=sigma(i,j,k,4)
      sigmayz=sigma(i,j,k,5)
      sigmazz=-sigmaxx-sigmayy
      !
      dtdx=dtmp(i,j,k,1)
      dtdy=dtmp(i,j,k,2)
      dtdz=dtmp(i,j,k,3)
      !
      pxx0=dsigmaxx(i,j,k,1)+dsigmaxy(i,j,k,2)+dsigmaxz(i,j,k,3)
      pyy0=dsigmaxy(i,j,k,1)+dsigmayy(i,j,k,2)+dsigmayz(i,j,k,3)
      pzz0=dsigmaxz(i,j,k,1)+dsigmayz(i,j,k,2)+dsigmazz(i,j,k,3)
      !
      pxx=dprs(i,j,k,1)+pxx0 
      pyy=dprs(i,j,k,2)+pyy0
      pzz=dprs(i,j,k,3)+pzz0 

      !
      qx=qflux(i,j,k,1)
      qy=qflux(i,j,k,2)
      qz=qflux(i,j,k,3)
      !
      mxzz=-mijk(i,j,k,1)-mijk(i,j,k,4)
      myzz=-mijk(i,j,k,2)-mijk(i,j,k,5)
      mzzz=-mijk(i,j,k,3)-mijk(i,j,k,6)

      pres = prs(i,j,k)
      !
      ! sigmaxx
      srcterm(1)= -pres*(sigmaxx/miup + num2d3*(2.d0*dudx-dvdy-dwdz))  &
                  -num4d15*(2.d0*dqxdx-dqydy-dqzdz)                    &
                  -num2d3*( sigmaxx*(2.d0*dudx+dwdz)    &
                           +sigmayy*(dwdz-dvdy)         &
                           +sigmaxy*(2.d0*dudy-dvdx)    &
                           +sigmaxz*(2.d0*dudz-dwdx)    &
                           -sigmayz*(dwdy+dvdz))
      ! sigmaxy
      srcterm(2)= -pres*(sigmaxy/miup + (dvdx+dudy) )                  &
                  -0.4d0*(dqydx+dqxdy)                                 &
                  -( sigmaxx*dvdx + sigmayy*dudy + sigmaxy*(dudx+dvdy) &
                    +sigmaxz*dvdz + sigmayz*dudz )
      ! sigmaxz
      srcterm(3)= -pres*(sigmaxz/miup + (dwdx+dudz) )                  &
                  -0.4d0*(dqzdx+dqxdz)                                 &
                  -( sigmaxx*(dwdx-dudz) - sigmayy*dudz   &
                    +sigmaxz*(dudx+dwdz) + sigmaxy*dwdy + sigmayz*dudy )
      ! sigmayy
      srcterm(4)= -pres*(sigmayy/miup + num2d3*(2.d0*dvdy-dudx-dwdz))  &
                  -num4d15*(2.d0*dqydy-dqxdx-dqzdz)                    &
                  -num2d3*( sigmaxx*(dwdz-dudx)         &
                           +sigmayy*(2.d0*dvdy+dwdz)    &
                           +sigmaxy*(2.d0*dvdx-dudy)    &
                           -sigmaxz*(dudz+dwdx)         &
                           +sigmayz*(2.d0*dvdz-dwdy))
      ! sigmayz
      srcterm(5)= -pres*(sigmayz/miup + (dwdy+dvdz) )                  &
                  -0.4d0*(dqydz+dqzdy)                                 &
                  -(-sigmaxx*dvdz + sigmayy*(dwdy-dvdz)   &
                    +sigmaxy*dwdx + sigmaxz*dvdx + sigmayz*(dvdy+dwdz) )
      !
      ! qx
      srcterm(6)= -pres*(Aq*qx/miup + dtmp(i,j,k,1)/const5)          &
                  -( 3.5d0*(sigmaxx*dtdx+sigmaxy*dtdy+sigmaxz*dtdz)  &
                    + tmp(i,j,k)*pxx0 )/const2                       &
                  +(Pxx*sigmaxx+Pyy*sigmaxy+Pzz*sigmaxz)/rho(i,j,k)  &
                  -0.4d0*( qx*(5.5d0*dudx+dvdy+dwdz)                 &
                         + qy*(3.5d0*dudy+dvdx) + qz*(3.5d0*dudz+dwdx) ) &
                  -( mijk(i,j,k,1)*dudx+mijk(i,j,k,4)*dvdy+mxzz*dwdz     &
                   + mijk(i,j,k,2)*(dudy+dvdx)+mijk(i,j,k,3)*(dudz+dwdx) &
                   + mijk(i,j,k,7)*(dwdy+dvdz) ) - num1d6*ddelta(i,j,k,1)

      ! qy
      srcterm(7)= -pres*(Aq*qy/miup + dtmp(i,j,k,2)/const5)           &
                  -( 3.5d0*(sigmaxy*dtdx+sigmayy*dtdy+sigmayz*dtdz)   &
                    + tmp(i,j,k)*pyy0 )/const2                        &
                  +(Pxx*sigmaxy+Pyy*sigmayy+Pzz*sigmayz)/rho(i,j,k)   &
                  -0.4d0*( qx*(3.5d0*dvdx+dudy)                       &
                         + qy*(dudx+5.5d0*dvdy+dwdz) + qz*(3.5d0*dvdz+dwdy) ) &
                  -( mijk(i,j,k,2)*dudx + mijk(i,j,k,5)*dvdy + myzz*dwdz    &
                   + mijk(i,j,k,4)*(dudy+dvdx) + mijk(i,j,k,7)*(dudz+dwdx)  &
                   + mijk(i,j,k,6)*(dwdy+dvdz) ) - num1d6*ddelta(i,j,k,2)

      ! qz
      srcterm(8)= -pres*(Aq*qz/miup + dtmp(i,j,k,3)/const5)           &
                  -( 3.5d0*(sigmaxz*dtdx+sigmayz*dtdy+sigmazz*dtdz)   &
                    + tmp(i,j,k)*pzz0 )/const2                        &
                  +(Pxx*sigmaxz+Pyy*sigmayz+Pzz*sigmazz)/rho(i,j,k)   &
                  -0.4d0*( qx*(3.5d0*dwdx+dudz)                       &
                         + qy*(3.5d0*dwdy+dvdz) + qz*(dudx+dvdy+5.5d0*dwdz) ) &
                  -( mijk(i,j,k,3)*dudx + mijk(i,j,k,6)*dvdy + mzzz*dwdz  &
                   + mijk(i,j,k,7)*(dudy+dvdx) + mxzz*(dudz+dwdx)         &
                   + myzz*(dwdy+dvdz) ) - num1d6*ddelta(i,j,k,3)

      !
      qrhs_mom(i,j,k,1)=qrhs_mom(i,j,k,1)+srcterm(1)*jacob(i,j,k)
      qrhs_mom(i,j,k,2)=qrhs_mom(i,j,k,2)+srcterm(2)*jacob(i,j,k)
      qrhs_mom(i,j,k,3)=qrhs_mom(i,j,k,3)+srcterm(3)*jacob(i,j,k)
      qrhs_mom(i,j,k,4)=qrhs_mom(i,j,k,4)+srcterm(4)*jacob(i,j,k)
      qrhs_mom(i,j,k,5)=qrhs_mom(i,j,k,5)+srcterm(5)*jacob(i,j,k)
      !
      qrhs_mom(i,j,k,6)=qrhs_mom(i,j,k,6)+srcterm(6)*jacob(i,j,k)
      qrhs_mom(i,j,k,7)=qrhs_mom(i,j,k,7)+srcterm(7)*jacob(i,j,k)
      qrhs_mom(i,j,k,8)=qrhs_mom(i,j,k,8)+srcterm(8)*jacob(i,j,k)
      !
      ! srcmax=max(srcmax,maxval(srcterm))
      ! srcmin=min(srcmin,minval(srcterm))
    enddo
    enddo 
    enddo
    !
    ! srcmax=pmax(srcmax)
    ! srcmin=pmin(srcmin)
    ! !
    ! if(mpirank==0) print*,' **  r13 max src=',srcmax,'min src=',srcmin
    !
  end subroutine src_sigma_q
  !+-------------------------------------------------------------------+
  !| The end of the subroutine src_sigma_q.                            |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to calculate the RHS of the moment        |
  !| transport equations.                                              |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 10-04-2024 | Copied from gitlab by J. Fang @ Warringon.           |
  !+-------------------------------------------------------------------+
  subroutine rhsmomcal(timerept)
    !
    use commvar,  only : flowtype,conschm,diffterm,recon_schem,limmbou,moment
    use comsolver,only : gradcal
    !
    ! arguments
    logical,intent(in),optional :: timerept
    !
    ! local data

    logical,save :: firstcall=.true.
    integer :: j, n
    integer :: nconv
    !
    real(8) :: time_beg
    real(8),save :: subtime=0.d0
    !
    !+-------------------------+
    !| Begining part           |
    !+-------------------------+
    if(present(timerept) .and. timerept) time_beg=ptime()
    !
    call gradcal()
    !
    qrhs_mom=0.d0
    !
    read(conschm(1:1),*) nconv
    !
    if(mod(nconv,2)==0) then
      !
      call convection_mom(timerept=ltimrpt)
      !
    else
      !
      stop ' !! scheme not enabled for moment equations @ rhsmomcal.'
      !
    endif
    !
    qrhs_mom=-qrhs_mom
    !
    ! call heatflux
    ! call dqcal
    !
    if(moment=='r13') then
      !
      call mijkcal
      !
      call rijcal
      !
      call deltacal
      !
    endif
    !
    call ddeltacal
    !
    if(moment=='r13' .or. moment=='r26') then
      !
      call diffu_sigma_q
      
      call src_sigma_q
      !
    endif
    !
    if(moment=='r26') then
      !
      call diffu_mijk
      !
      call diffu_rij
      !
      call diffu_delta
      !
      call src_mijk_A
      !
      call src_mijk_B
      !
      call src_rij_A
      !
      call src_rij_B
      !
      call src_delta
      !
    endif
    !+-------------------------+
    !| Ending part             |
    !+-------------------------+
    if(present(timerept) .and. timerept) then
      !
      subtime=subtime+ptime()-time_beg
      !
      if(lio .and. lreport) call timereporter(routine='rhsmomcal', &
                                              timecost=subtime,    &
                                              message='RHS MoM term')
    endif
    !
    return
    !
  end subroutine rhsmomcal
  !+-------------------------------------------------------------------+
  !| The of the subroutine rhsmomcal                                   |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutin is to calculate the mijk term in the r13 moment    |
  !+-------------------------------------------------------------------+
  !  checked on 22/08/23                                               |
  !+-------------------------------------------------------------------+
  subroutine mijkcal

    use constdef
    use commvar,  only : im,jm,km,gamma,Reynolds,mach,difschm,         &
                         npdci,npdcj,npdck,num_species,num_modequ,     &
                         is,ie,js,je,ks,ke,const2
    use commarray,only : tmp,prs,sigma,qflux,miu
    use comsolver,only : grad
    use parallel, only : dataswap
    use tecio

    implicit none
    !
    ! local data
    !
    integer :: i,j,k
    real(8) :: co1,co2,Fxxx,Fxxy,Fxxz,Fxyy,Fyyy,Fyyz,Fxyz
    !
    logical,save :: firstcall=.true.
    !
    call dataswap(sigma)
    !
    dsigmaxx=grad(sigma(:,:,:,1))
    dsigmaxy=grad(sigma(:,:,:,2))
    dsigmaxz=grad(sigma(:,:,:,3))
    dsigmayy=grad(sigma(:,:,:,4))
    dsigmayz=grad(sigma(:,:,:,5))
    dsigmazz= - dsigmaxx - dsigmayy
    !
    !
    do k=0,km
    do j=0,jm
    do i=0,im
      !
      !
      co1=-3.d0*miu(i,j,k)*tmp(i,j,k)/(const2*Am*prs(i,j,k))
      !
      Fxxx=  0.6d0 * (dsigmaxx(i,j,k,1))-  0.4d0 *(dsigmaxy(i,j,k,2)+dsigmaxz(i,j,k,3))
      Fxxy= num1d3 * (dsigmaxx(i,j,k,2) +  1.6d0 * dsigmaxy(i,j,k,1)-  0.4d0 * (dsigmayy(i,j,k,2)+dsigmayz(i,j,k,3)))
      Fxxz= num1d3 * (dsigmaxx(i,j,k,3) +  1.6d0 * dsigmaxz(i,j,k,1)-  0.4d0 * (dsigmayz(i,j,k,2)+dsigmazz(i,j,k,3)))
      Fxyy= num1d3 * (dsigmayy(i,j,k,1) +  1.6d0 * dsigmaxy(i,j,k,2)-  0.4d0 * (dsigmaxx(i,j,k,1)+dsigmaxz(i,j,k,3)))
      Fyyy=  0.6d0 * (dsigmayy(i,j,k,2))-  0.4d0 *(dsigmaxy(i,j,k,1)+dsigmayz(i,j,k,3))
      Fyyz= num1d3 * (dsigmayy(i,j,k,3) +  1.6d0 *(dsigmayz(i,j,k,2)) -0.4d0 * (dsigmaxz(i,j,k,1)+dsigmazz(i,j,k,3)))
      Fxyz= num1d3 * (dsigmaxy(i,j,k,3)+dsigmaxz(i,j,k,2)+dsigmayz(i,j,k,1))
      !
      mijk(i,j,k,1)=Fxxx*co1
      mijk(i,j,k,2)=Fxxy*co1
      mijk(i,j,k,3)=Fxxz*co1
      mijk(i,j,k,4)=Fxyy*co1
      mijk(i,j,k,5)=Fyyy*co1
      mijk(i,j,k,6)=Fyyz*co1
      mijk(i,j,k,7)=Fxyz*co1
      !
    enddo
    enddo
    enddo
    !
    call dataswap(mijk,timerept=ltimrpt)
    !
  end subroutine mijkcal
  !+-------------------------------------------------------------------+
  !| The of the subroutine mijkcal                                     |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutin is to calculate the rij  term in the r13 moment    |
  !+-------------------------------------------------------------------+
  !  checked on 22/08/23                                               |
  !+-------------------------------------------------------------------+
  subroutine rijcal 
    !
    use commvar,   only : im,jm,km,reynolds,gamma,mach,const2
    use commarray, only : tmp,prs,miu
    use parallel,  only : dataswap

    implicit none
    !
    real(8) :: dqxdx,dqxdy,dqxdz,dqydx,dqydy,dqydz,dqzdx,dqzdy,dqzdz
    real(8) :: co2
    integer :: i,j,k
    !
    call dqfluxcal
    !
    do k=0,km
    do j=0,jm
    do i=0,im
      !
      !
      co2=-14.d0*miu(i,j,k)*tmp(i,j,k)/(15.d0*AR1*prs(i,j,k)*const2)
      !
      dqxdx=dqx(i,j,k,1)
      dqxdy=dqx(i,j,k,2)
      dqxdz=dqx(i,j,k,3)
      dqydx=dqy(i,j,k,1)
      dqydy=dqy(i,j,k,2)
      dqydz=dqy(i,j,k,3)
      dqzdx=dqz(i,j,k,1)
      dqzdy=dqz(i,j,k,2)
      dqzdz=dqz(i,j,k,3)
      !
      rij(i,j,k,1)=2.d0*co2*(2.d0*dqxdx-dqydy-dqzdz)  ! rxx
      rij(i,j,k,2)=3.d0*co2*(dqydx+dqxdy)             ! rxy    
      rij(i,j,k,3)=3.d0*co2*(dqzdx+dqxdz)             ! rxz
      rij(i,j,k,4)=2.d0*co2*(2.d0*dqydy-dqxdx-dqzdz)  ! ryy   
      rij(i,j,k,5)=3.d0*co2*(dqydz+dqzdy)             ! ryz
      !
    enddo
    enddo
    enddo
    !
    call dataswap(rij,timerept=ltimrpt)
    !
  end subroutine rijcal
  !+-------------------------------------------------------------------+
  !| The of the subroutine rijcal                                      |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutin is to calculate the delta term in the r13 moment   |
  !+-------------------------------------------------------------------+
  !  checked on 22/08/23                                               |
  !+-------------------------------------------------------------------+
  subroutine deltacal
    !
    use commvar,   only : const2
    use commarray, only : miu,tmp,prs

    implicit none
    !
    real(8) :: co1
    !
    co1=-8.d0/Adelta1/const2
    delta(0:im,0:jm,0:km)=co1*miu(0:im,0:jm,0:km)*tmp(0:im,0:jm,0:km)/prs(0:im,0:jm,0:km)* &
                           (dqx(0:im,0:jm,0:km,1)+dqy(0:im,0:jm,0:km,2)+dqy(0:im,0:jm,0:km,3))
    !
  end subroutine deltacal
  !+-------------------------------------------------------------------+
  !| The of the subroutine deltacal                                    |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is calculate gradient of delta                    |
  !+-------------------------------------------------------------------+
  subroutine ddeltacal
    !
    use comsolver,only : grad
    use parallel,  only : dataswap
    !
    call dataswap(delta)
    !
    ddelta=grad(delta)
    !
  end subroutine ddeltacal
  !+-------------------------------------------------------------------+
  !| The of the subroutine ddeltacal                                   |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is calculate gradient of qflux                    |
  !+-------------------------------------------------------------------+
  subroutine dqfluxcal
    !
    use commarray, only : qflux, vel, dvel, tmp, dtmp
    use comsolver, only : grad
    use parallel,  only : dataswap

    implicit none
    !
    call dataswap(qflux)
    !
    dqx=grad(qflux(:,:,:,1))
    dqy=grad(qflux(:,:,:,2))
    dqz=grad(qflux(:,:,:,3))

    ! dvel(0:im,0:jm,0:km,1,:)=grad(vel(:,:,:,1))
    ! dvel(0:im,0:jm,0:km,2,:)=grad(vel(:,:,:,2))
    ! dvel(0:im,0:jm,0:km,3,:)=grad(vel(:,:,:,3))

    ! call dataswap(tmp)

    ! dtmp = grad(tmp)
    !
  end subroutine dqfluxcal
  !+-------------------------------------------------------------------+
  !| The of the subroutine dqfluxcal                                   |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to swap qmom and update the moment        |
  !| variables. The variables at interfaces are also synchronized.     |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 08-Feb-2021 | Created by J. Fang @ Warringon.                     |
  !+-------------------------------------------------------------------+
  subroutine qmomswap(timerept)
    !
    use mpi
    use commvar,  only: num_xtrmom,ia,ja,ka,lreport
    use parallel, only: isize,jsize,ksize,lio,lihomo,ljhomo,lkhomo,    &
                        mpileft,mpiright,mpidown,mpiup,mpifront,mpiback, &
                        mpitag,status
    !
    ! argument
    logical,intent(in),optional :: timerept
    !
    ! local data
    integer :: ncou
    integer :: ierr,j,k
    real(8),allocatable,dimension(:,:,:,:) :: sbuf1,sbuf2,rbuf1,rbuf2
    real(8) :: time_beg
    real(8),save :: subtime=0.d0
    !
    if(present(timerept) .and. timerept) time_beg=ptime() 
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! buf1: send buffer
    ! buf2: send buffer
    ! buf1: redevice buffer
    ! buf2: redevice buffer
    !!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Message pass in i direction.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if(isize==1) then
    else
      ncou=(jm+1)*(km+1)*num_xtrmom*(hm+1)
      !
      allocate( sbuf1(0:hm, 0:jm,0:km,1:num_xtrmom),                        &
                sbuf2(0:hm, 0:jm,0:km,1:num_xtrmom),                        &
                rbuf1(0:hm, 0:jm,0:km,1:num_xtrmom),                        &
                rbuf2(-hm:0,0:jm,0:km,1:num_xtrmom) )
      !
      if(mpileft .ne. MPI_PROC_NULL) then
        ! pack the left send buffer
        sbuf1(0:hm,0:jm,0:km,:)=q_mom(0:hm,0:jm,0:km,:)
      endif
      if(mpiright .ne. MPI_PROC_NULL) then
        ! pack the right send buffer
        sbuf2(0:hm,0:jm,0:km,:)=q_mom(im-hm:im,0:jm,0:km,:)
      endif
      !
      ! Message passing
      call mpi_sendrecv(sbuf1,ncou,mpi_real8,mpileft, mpitag,            &
                        rbuf1,ncou,mpi_real8,mpiright,mpitag,            &
                                               mpi_comm_world,status,ierr)
      mpitag=mpitag+1
      call mpi_sendrecv(sbuf2,ncou,mpi_real8,mpiright,mpitag,            &
                        rbuf2,ncou,mpi_real8,mpileft, mpitag,            &
                                               mpi_comm_world,status,ierr)
      mpitag=mpitag+1
      !
      if(mpiright .ne. MPI_PROC_NULL) then
        !
        ! unpack the received the packet from right
        q_mom(im+1:im+hm,0:jm,0:km,:)=rbuf1(1:hm,0:jm,0:km,:)
        !
        q_mom(im,0:jm,0:km,:)=0.5d0*( q_mom(im,0:jm,0:km,:) +            &
                                  rbuf1(0,0:jm,0:km,:) )
        !
        call updatemoment(im,im+hm,0,jm,0,km)
        ! 
      end if
      !
      if(mpileft .ne. MPI_PROC_NULL) then
        !
        ! unpack the received the packet from left
        q_mom(-hm:-1,0:jm,0:km,:)=rbuf2(-hm:-1,0:jm,0:km,:)
        !
        q_mom(0,0:jm,0:km,:)=0.5d0*( q_mom(0,0:jm,0:km,:) +              &
                                 rbuf2(0,0:jm,0:km,:) )
        !
        call updatemoment(-hm,0,0,jm,0,km)
        !
      end if
      !
      deallocate( sbuf1,sbuf2,rbuf1,rbuf2 )
      !
    endif
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Finish message pass in i direction.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    if(jsize==1 .and. ljhomo) then
      !
      if(jm==0) then
        do j=-hm,hm
          q_mom(0:im,j,0:km,:)=q_mom(0:im,0,0:km,:)
        enddo
      else
        q_mom(0:im, -hm:-1   ,0:km,:)=q_mom(0:im,jm-hm:jm-1,0:km,:)
        q_mom(0:im,jm+1:jm+hm,0:km,:)=q_mom(0:im,    1:hm,  0:km,:)
        !
        q_mom(0:im, 0,0:km,:)=0.5d0*(q_mom(0:im,0,0:km,:)+q_mom(0:im,jm,0:km,:))
        q_mom(0:im,jm,0:km,:)=q_mom(0:im,0,0:km,:)
      endif
      !
      call updatemoment(0,im,-hm,0,0,km)
      call updatemoment(0,im,jm,jm+hm,0,km)
      !
    else
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Message pass in j direction.
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ncou=(im+1)*(km+1)*num_xtrmom*(hm+1)
      !
      allocate( sbuf1(0:im,0:hm, 0:km,1:num_xtrmom),                        &
                sbuf2(0:im,0:hm, 0:km,1:num_xtrmom),                        &
                rbuf1(0:im,0:hm, 0:km,1:num_xtrmom),                        &
                rbuf2(0:im,-hm:0,0:km,1:num_xtrmom) )
      !
      if(mpidown .ne. MPI_PROC_NULL) then
        ! pack the upper send buffer
        sbuf1(0:im,0:hm,0:km,:)=q_mom(0:im,0:hm,0:km,:)
      endif
      if(mpiup .ne. MPI_PROC_NULL) then
        ! pack the down send buffer
        sbuf2(0:im,0:hm,0:km,:)=q_mom(0:im,jm-hm:jm,0:km,:)
      end if
      !
      ! Message passing
      call mpi_sendrecv(sbuf1,ncou,mpi_real8,mpidown,mpitag,             &
                        rbuf1,ncou,mpi_real8,mpiup,  mpitag,             &
                                               mpi_comm_world,status,ierr)
      mpitag=mpitag+1
      call mpi_sendrecv(sbuf2,ncou,mpi_real8,mpiup,  mpitag,             &
                        rbuf2,ncou,mpi_real8,mpidown,mpitag,             &
                                               mpi_comm_world,status,ierr)
      mpitag=mpitag+1
      !
      if(mpiup .ne. MPI_PROC_NULL) then
        ! unpack the received the packet from up
        q_mom(0:im,jm+1:jm+hm,0:km,:)=rbuf1(0:im,1:hm,0:km,:)
        !
        q_mom(0:im,jm,0:km,:)=0.5d0*( q_mom(0:im,jm,0:km,:) +            &
                                  rbuf1(0:im, 0,0:km,:) )
        !
        call updatemoment(0,im,jm,jm+hm,0,km)
        !
      endif
      !
      if(mpidown .ne. MPI_PROC_NULL) then
        ! unpack the received the packet from down
        q_mom(0:im,-hm:-1,0:km,:)=rbuf2(0:im,-hm:-1,0:km,:) 
        !
        q_mom(0:im,0,0:km,:)=0.5d0*( q_mom(0:im, 0,0:km,:) +            &
                                 rbuf2(0:im, 0,0:km,:) )
        !
        call updatemoment(0,im,-hm,0,0,km)
        !
      end if
      !
      deallocate( sbuf1,sbuf2,rbuf1,rbuf2 )
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Finish message pass in j direction.
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    endif
    !
    if(ksize==1 .and. lkhomo) then
      !
      if(ka==0) then
        do k=-hm,hm
          q_mom(0:im,0:jm,k,:)=q_mom(0:im,0:jm,0,:)
        enddo
      else
        q_mom(0:im,0:jm, -hm:-1   ,:)=q_mom(0:im,0:jm,km-hm:km-1,:)
        q_mom(0:im,0:jm,km+1:km+hm,:)=q_mom(0:im,0:jm,    1:hm,  :)
        !
        q_mom(0:im,0:jm,0,:)=0.5d0*(q_mom(0:im,0:jm,0,:)+q_mom(0:im,0:jm,km,:))
        q_mom(0:im,0:jm,km,:)=q_mom(0:im,0:jm,0,:)
      endif
      !
      call updatemoment(0,im,0,jm,-hm,0)
      call updatemoment(0,im,0,jm,km,km+hm)
      !
    else
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Message pass in k direction.
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ncou=(im+1)*(jm+1)*num_xtrmom*(hm+1)
      !
      allocate( sbuf1(0:im,0:jm, 0:hm,1:num_xtrmom),                      &
                sbuf2(0:im,0:jm, 0:hm,1:num_xtrmom),                      &
                rbuf1(0:im,0:jm, 0:hm,1:num_xtrmom),                      &
                rbuf2(0:im,0:jm,-hm:0,1:num_xtrmom) )
      !
      if(mpiback .ne. MPI_PROC_NULL) then
        ! pack the back send buffer
        sbuf1(0:im,0:jm,0:hm,:)=q_mom(0:im,0:jm,0:hm,:)
      endif
      if(mpifront .ne. MPI_PROC_NULL) then
        ! pack the front send buffer
        sbuf2(0:im,0:jm,0:hm,:)=q_mom(0:im,0:jm,km-hm:km,:)
      endif
      !
      ! Message passing
      call mpi_sendrecv(sbuf1,ncou,mpi_real8,mpiback, mpitag,          &
                        rbuf1,ncou,mpi_real8,mpifront,mpitag,          &
                                             mpi_comm_world,status,ierr)
      mpitag=mpitag+1
      call mpi_sendrecv(sbuf2,ncou,mpi_real8,mpifront,mpitag,          &
                        rbuf2,ncou,mpi_real8,mpiback, mpitag,          &
                                             mpi_comm_world,status,ierr)
      mpitag=mpitag+1
      !
      if(mpifront .ne. MPI_PROC_NULL) then
        !
        ! unpack the received the packet from front
        q_mom(0:im,0:jm,km+1:km+hm,:)=rbuf1(0:im,0:jm,1:hm,:)
        !
        q_mom(0:im,0:jm,km,:)=0.5d0*( q_mom(0:im,0:jm,km,:) +              &
                                  rbuf1(0:im,0:jm, 0,:) )
        !
        call updatemoment(0,im,0,jm,km,km+hm)
        !
      end if
      !
      if(mpiback .ne. MPI_PROC_NULL) then
        !
        ! unpack the received the packet back
        q_mom(0:im,0:jm,-hm:-1,:)=rbuf2(0:im,0:jm,-hm:-1,:)
        !
        q_mom(0:im,0:jm,0,:)=0.5d0*( q_mom(0:im,0:jm,0,:) +                &
                                 rbuf2(0:im,0:jm,0,:)  )
        !
        call updatemoment(0,im,0,jm,-hm,0)
        !
      end if
      !
      deallocate( sbuf1,sbuf2,rbuf1,rbuf2 )
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Finish message pass in k direction.
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    endif
    !
    if(mpitag>10000) mpitag=100
    ! reset mpitag
    !
    if(present(timerept) .and. timerept) then
      !
      subtime=subtime+ptime()-time_beg
      !
      if(lio .and. lreport) call timereporter(routine='qswap', &
                                              timecost=subtime )
    endif
    !
    return
    !
  end subroutine qmomswap
  !+-------------------------------------------------------------------+
  !| The end of the subroutine qmomswap.                               |
  !+-------------------------------------------------------------------+
  !!
  !+-------------------------------------------------------------------+
  !| This subroutine updated q_mom                                     |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 23-11-2023: Created by J. Fang @ Warrington.                      |
  !+-------------------------------------------------------------------+
  subroutine updateqmom(ibeg,iend,jbeg,jend,kbeg,kend)
    !
    use commarray,only : sigma,qflux
    use commvar,  only : moment
    !
    ! arguments
    integer,intent(in),optional :: ibeg,iend,jbeg,jend,kbeg,kend
    ! local data
    integer :: i,j,k
    integer :: iss,iee,jss,jee,kss,kee
    !
    iss=-hm
    iee=im+hm
    jss=-hm
    jee=jm+hm
    kss=-hm
    kee=km+hm
    if(present(ibeg)) iss=ibeg
    if(present(iend)) iee=iend
    if(present(jbeg)) jss=jbeg
    if(present(jend)) jee=jend
    if(present(kbeg)) kss=kbeg
    if(present(kend)) kee=kend
    !
    if(moment=='r13') then
      call fvar2qmom(q=q_mom(iss:iee,jss:jee,kss:kee,:),                        &
                                      stress=sigma(iss:iee,jss:jee,kss:kee,:),  &
                                    heatflux=qflux(iss:iee,jss:jee,kss:kee,:) )
    elseif(moment=='r26') then
      call fvar2qmom(q=q_mom(iss:iee,jss:jee,kss:kee,:),                        &
                                      stress= sigma(iss:iee,jss:jee,kss:kee,:), &
                                    heatflux= qflux(iss:iee,jss:jee,kss:kee,:), &
                                        mijk=  mijk(iss:iee,jss:jee,kss:kee,:), &
                                         rij=   rij(iss:iee,jss:jee,kss:kee,:), &
                                       delta= delta(iss:iee,jss:jee,kss:kee) )
    endif
    !
  end subroutine updateqmom
  !+-------------------------------------------------------------------+
  !| The end of the subroutine updateqmom.                             |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine updated moment from q_mom                         |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 23-11-2023: Created by J. Fang @ Warrington.                      |
  !+-------------------------------------------------------------------+
  subroutine updatemoment(ibeg,iend,jbeg,jend,kbeg,kend)
    !
    use commarray,only : sigma,qflux
    use commvar,  only : moment
    !
    ! arguments
    integer,intent(in),optional :: ibeg,iend,jbeg,jend,kbeg,kend
    ! local data
    integer :: i,j,k
    integer :: iss,iee,jss,jee,kss,kee
    !
    iss=-hm
    iee=im+hm
    jss=-hm
    jee=jm+hm
    kss=-hm
    kee=km+hm
    if(present(ibeg)) iss=ibeg
    if(present(iend)) iee=iend
    if(present(jbeg)) jss=jbeg
    if(present(jend)) jee=jend
    if(present(kbeg)) kss=kbeg
    if(present(kend)) kee=kend
    !
    if(moment=='r13') then
      call qmom2fvar(q=q_mom(iss:iee,jss:jee,kss:kee,:),                        &
                                      stress=sigma(iss:iee,jss:jee,kss:kee,:),  &
                                    heatflux=qflux(iss:iee,jss:jee,kss:kee,:) )
    elseif(moment=='r26') then
      call qmom2fvar(q=q_mom(iss:iee,jss:jee,kss:kee,:),                        &
                                      stress= sigma(iss:iee,jss:jee,kss:kee,:), &
                                    heatflux= qflux(iss:iee,jss:jee,kss:kee,:), &
                                        mijk=  mijk(iss:iee,jss:jee,kss:kee,:), &
                                         rij=   rij(iss:iee,jss:jee,kss:kee,:), &
                                       delta= delta(iss:iee,jss:jee,kss:kee) )
    endif
    !
  end subroutine updatemoment
  !+-------------------------------------------------------------------+
  !| The end of the subroutine updatemoment.                           |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine calcualtes stress and heat flux from q.           |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 23-11-2022: Created by J. Fang @ Warrington.                      |
  !+-------------------------------------------------------------------+
  subroutine qmom2sigma_flux_3da(q,stress,heatflux,mijk,rij,delta)
    !
    real(8),intent(in) :: q(:,:,:,:)
    real(8),intent(out),optional :: stress(:,:,:,:),heatflux(:,:,:,:), &
                                    mijk(:,:,:,:),rij(:,:,:,:),        &
                                    delta(:,:,:)
    !
    if(present(stress)) then
      stress(:,:,:,1)=q(:,:,:,1)
      stress(:,:,:,2)=q(:,:,:,2)
      stress(:,:,:,3)=q(:,:,:,3)
      stress(:,:,:,4)=q(:,:,:,4)
      stress(:,:,:,5)=q(:,:,:,5)
      stress(:,:,:,6)=-stress(:,:,:,1)-stress(:,:,:,4)
    endif
    !
    if(present(heatflux)) then
      heatflux(:,:,:,1)=q(:,:,:,6)
      heatflux(:,:,:,2)=q(:,:,:,7)
      heatflux(:,:,:,3)=q(:,:,:,8)
    endif
    !
    if(present(mijk)) then
      mijk(:,:,:,1)=q(:,:,:,9)    ! mxxx
      mijk(:,:,:,2)=q(:,:,:,10)   ! mxxy
      mijk(:,:,:,3)=q(:,:,:,11)   ! mxxz
      mijk(:,:,:,4)=q(:,:,:,12)   ! mxyy
      mijk(:,:,:,5)=q(:,:,:,13)   ! myyy 
      mijk(:,:,:,6)=q(:,:,:,14)   ! myyz
      mijk(:,:,:,7)=q(:,:,:,15)   ! mxyz 
    endif
    !
    if(present(rij)) then
      rij(:,:,:,1)=q(:,:,:,16)   ! Rxx
      rij(:,:,:,2)=q(:,:,:,17)   ! Rxy
      rij(:,:,:,3)=q(:,:,:,18)   ! Rxz
      rij(:,:,:,4)=q(:,:,:,19)   ! Ryy
      rij(:,:,:,5)=q(:,:,:,20)   ! Ryz 
    endif
    !
    if(present(delta)) then
      delta(:,:,:)=q(:,:,:,21)
    endif
    !
  end subroutine qmom2sigma_flux_3da
  !
  subroutine qmom2sigma_flux_1da(q,stress,heatflux,mijk,rij,delta)
    !
    real(8),intent(in) :: q(:,:)
    real(8),intent(out),optional :: stress(:,:),heatflux(:,:),mijk(:,:),rij(:,:),delta(:)
    !
    if(present(stress)) then
      stress(:,1)=q(:,1)
      stress(:,2)=q(:,2)
      stress(:,3)=q(:,3)
      stress(:,4)=q(:,4)
      stress(:,5)=q(:,5)
      stress(:,6)=-stress(:,1)-stress(:,4)
    endif
    !
    if(present(heatflux)) then
      heatflux(:,1)=q(:,6)
      heatflux(:,2)=q(:,7)
      heatflux(:,3)=q(:,8)
    endif
    !
    if(present(mijk)) then
      mijk(:,1)=q(:,9)    ! mxxx
      mijk(:,2)=q(:,10)   ! mxxy
      mijk(:,3)=q(:,11)   ! mxxz
      mijk(:,4)=q(:,12)   ! mxyy
      mijk(:,5)=q(:,13)   ! myyy 
      mijk(:,6)=q(:,14)   ! myyz 
      mijk(:,7)=q(:,15)   ! mxyz
    endif
    !
    if(present(rij)) then
      rij(:,1)=q(:,16)   ! Rxx
      rij(:,2)=q(:,17)   ! Rxy
      rij(:,3)=q(:,18)   ! Rxz
      rij(:,4)=q(:,19)   ! Ryy
      rij(:,5)=q(:,20)   ! Ryz 
    endif
    !
    if(present(delta)) then
      delta(:)=q(:,21)
    endif
    !
  end subroutine qmom2sigma_flux_1da
  !!
  subroutine qmom2sigma_flux_sca(q,stress,heatflux,mijk,rij,delta)
    !
    real(8),intent(in) :: q(:)
    real(8),intent(out),optional :: stress(:),heatflux(:),mijk(:),rij(:),delta
    !
    if(present(stress)) then
      stress(1)=q(1)
      stress(2)=q(2)
      stress(3)=q(3)
      stress(4)=q(4)
      stress(5)=q(5)
      stress(6)=-stress(1)-stress(4)
    endif
    !
    if(present(heatflux)) then
      heatflux(1)=q(6)
      heatflux(2)=q(7)
      heatflux(3)=q(8)
    endif
    !
    if(present(mijk)) then
      mijk(1)=q(9)    ! mxxx
      mijk(2)=q(10)   ! mxxy
      mijk(3)=q(11)   ! mxxz
      mijk(4)=q(12)   ! mxyy
      mijk(5)=q(13)   ! myyy
      mijk(6)=q(14)   ! myyz
      mijk(7)=q(15)   ! mxyz
    endif
    !
    if(present(rij)) then
      rij(1)=q(16)   ! Rxx
      rij(2)=q(17)   ! Rxy
      rij(3)=q(18)   ! Rxz
      rij(4)=q(19)   ! Ryy
      rij(5)=q(20)   ! Ryz 
    endif
    !
    if(present(delta)) then
      delta=q(21)
    endif
    !
  end subroutine qmom2sigma_flux_sca
  !

  subroutine sigma_flux2qmom_sca(q,stress,heatflux,mijk,rij,delta)
    !
    real(8),intent(in),optional :: stress(:),heatflux(:),mijk(:),rij(:),delta
    real(8),intent(out) :: q(:)
    !
    if(present(stress)) then
      q(1)=stress(1)
      q(2)=stress(2)
      q(3)=stress(3)
      q(4)=stress(4)
      q(5)=stress(5)
    endif
    !
    if(present(heatflux)) then
      q(6)=heatflux(1)
      q(7)=heatflux(2)
      q(8)=heatflux(3)
    endif
    !
    if(present(mijk)) then
      q(9)  = mijk(1)   ! mxxx
      q(10) = mijk(2)   ! mxxy
      q(11) = mijk(3)   ! mxxz
      q(12) = mijk(4)   ! mxyy
      q(13) = mijk(5)   ! myyy 
      q(14) = mijk(6)   ! myyz
      q(15) = mijk(7)   ! mxyz 
    endif
    !
    if(present(rij)) then
      q(16)=rij(1)   ! Rxx
      q(17)=rij(2)   ! Rxy
      q(18)=rij(3)   ! Rxz
      q(19)=rij(4)   ! Ryy
      q(20)=rij(5)   ! Ryz 
    endif
    !
    if(present(delta)) then
      q(21)=delta
    endif
    !
  end subroutine sigma_flux2qmom_sca
  !
  subroutine sigma_flux2qmom_1da(q,stress,heatflux,mijk,rij,delta)
    !
    real(8),intent(in),optional :: stress(:,:),heatflux(:,:),mijk(:,:),rij(:,:),delta(:)
    real(8),intent(out) :: q(:,:)
    !
    if(present(stress)) then
      !
      q(:,1)=stress(:,1)
      q(:,2)=stress(:,2)
      q(:,3)=stress(:,3)
      q(:,4)=stress(:,4)
      q(:,5)=stress(:,5)
    endif
    !
    if(present(heatflux)) then
      !
      q(:,6)=heatflux(:,1)
      q(:,7)=heatflux(:,2)
      q(:,8)=heatflux(:,3)
    endif
    !
    if(present(mijk)) then
      q(:,9)  = mijk(:,1)   ! mxxx
      q(:,10) = mijk(:,2)   ! mxxy
      q(:,11) = mijk(:,3)   ! mxxz
      q(:,12) = mijk(:,4)   ! mxyy
      q(:,13) = mijk(:,5)   ! myyy 
      q(:,14) = mijk(:,6)   ! myyz
      q(:,15) = mijk(:,7)   ! mxyz 
    endif
    !
    if(present(rij)) then
      q(:,16)=rij(:,1)   ! Rxx
      q(:,17)=rij(:,2)   ! Rxy
      q(:,18)=rij(:,3)   ! Rxz
      q(:,19)=rij(:,4)   ! Ryy
      q(:,20)=rij(:,5)   ! Ryz 
    endif
    !
    if(present(delta)) then
      q(:,21)=delta(:)
    endif
    !
  end subroutine sigma_flux2qmom_1da
  !
  subroutine sigma_flux2qmom_3da(q,stress,heatflux,mijk,rij,delta)
    !
    real(8),intent(in),optional :: stress(:,:,:,:),heatflux(:,:,:,:),mijk(:,:,:,:),  &
                                   rij(:,:,:,:),delta(:,:,:)
    real(8),intent(out) :: q(:,:,:,:)
    !
    if(present(stress)) then
      q(:,:,:,1)=stress(:,:,:,1)
      q(:,:,:,2)=stress(:,:,:,2)
      q(:,:,:,3)=stress(:,:,:,3)
      q(:,:,:,4)=stress(:,:,:,4)
      q(:,:,:,5)=stress(:,:,:,5)
    endif
    !
    if(present(heatflux)) then
      q(:,:,:,6)=heatflux(:,:,:,1)
      q(:,:,:,7)=heatflux(:,:,:,2)
      q(:,:,:,8)=heatflux(:,:,:,3)
    endif
    !
    if(present(mijk)) then
      q(:,:,:,9)  = mijk(:,:,:,1)   ! mxxx
      q(:,:,:,10) = mijk(:,:,:,2)   ! mxxy
      q(:,:,:,11) = mijk(:,:,:,3)   ! mxxz
      q(:,:,:,12) = mijk(:,:,:,4)   ! mxyy
      q(:,:,:,13) = mijk(:,:,:,5)   ! myyy 
      q(:,:,:,14) = mijk(:,:,:,6)   ! myyz
      q(:,:,:,15) = mijk(:,:,:,7)   ! mxyz 
    endif
    !
    if(present(rij)) then
      q(:,:,:,16)=rij(:,:,:,1)   ! Rxx
      q(:,:,:,17)=rij(:,:,:,2)   ! Rxy
      q(:,:,:,18)=rij(:,:,:,3)   ! Rxz
      q(:,:,:,19)=rij(:,:,:,4)   ! Ryy
      q(:,:,:,20)=rij(:,:,:,5)   ! Ryz 
    endif
    !
    if(present(delta)) then
      q(:,:,:,21)=delta(:,:,:)
    endif
    !
  end subroutine sigma_flux2qmom_3da
  !+-------------------------------------------------------------------+
  !| The end of the subroutine fvar2qom.                               |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to call b.c. for N-S with moment.              |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| xx-xx-xxxx: Created by X. Gu.                                     |
  !+-------------------------------------------------------------------+
  subroutine MOM_wall_boundary(ndir)
    !==================================

    use parallel,only: lio,mpistop,mpirank,mpirankname,irk,jrk,krk,      &
                        irkm,jrkm,krkm,pmax,ptime
    use commvar, only: hm,im,jm,km,uinf,vinf,winf,pinf,roinf,tinf,ndims, &
                     num_species,flowtype,gamma,numq,npdci,npdcj,      &
                     npdck,is,ie,js,je,ks,ke,xmin,xmax,ymin,ymax,      &
                     zmin,zmax,time,num_modequ,ltimrpt, moment

    integer,intent(in) :: ndir
    !
    ! local data
    integer :: i,j,k,l,ip, jp, kp
    real(8) n1, n2, n3, uwall, vwall, wwall, twall, alpha
    !
    !
    if (ndir==1) then
      ! left wall
      !
      if (irk==0) then
        !
        i = 0 ; ip = i + 1
        !
        !--------------------------------
        !-- wall normal vector ------
        !---------------------------
        n1 = 1.d0; n2=0.d0;  n3=0.d0
        !----------------------------
        !-- wall velocity --
        !---------------------------
        uwall = 0.0d0; vwall = 0.0d0 ; wwall = 0.0d0
        !--------------------------------------------
        !-- wall temperature --
        !---------------------
        twall = 1.0d0   ! wall temperature is equal to the ref temperature
        !-------------------------------------
        !-- Accomodation coff
        !------------------------
        alpha = 1.0d0   ! diffusive wall
        if (moment=='r13') then
          do k=0,km
          do j=0,jm
            jp = j ; kp = k
            call R13wbc(n1, n2, n3, i, j, k, ip, jp, kp, &
                        uwall, vwall, wwall, twall, alpha)
          end do
          end do
        elseif (moment=='r26') then
          do k=0,km
          do j=0,jm
            jp = j ; kp = k
            call R26wbc(n1, n2, n3, i, j, k,ip, jp, kp, &
                        uwall, vwall, wwall, twall, alpha)
          end do
          end do
        end if
        !
      end if
      !
    elseif (ndir==2) then
      ! right wall
      !
      if (irk==irkm) then
        !
        i=im ; ip = i - 1
        !
        !--------------------------------
        !-- wall normal vector ------
        !---------------------------
        n1 = -1.d0; n2=0.d0;  n3=0.d0
        !----------------------------
        !-- wall velocity --
        !---------------------------
        uwall = 0.0d0; vwall = 0.0d0 ; wwall = 0.0d0
        !--------------------------------------------
        !-- wall temperature --
        !---------------------
        twall = 1.0d0   ! wall temperature is equal to the ref temperature
        !-------------------------------------
        !-- Accomodation coff
        !------------------------
        alpha = 1.0d0   ! diffusive wall
        if (moment=='r13') then
          do k=0,km
          do j=0,jm
            jp = j ; kp = k
            call R13wbc(n1, n2, n3, i, j, k, ip, jp, kp, &
                        uwall, vwall, wwall, twall, alpha)
          end do
          end do
        elseif (moment=='r26') then
          do k=0,km
          do j=0,jm
            jp = j ; kp = k
            call R26wbc(n1, n2, n3, i, j, k,ip, jp, kp, &
                        uwall, vwall, wwall, twall, alpha)
          end do
          end do
        end if
        !
      end if
      !
    elseif (ndir==3) then
      ! bottom wall
      !
      if (jrk==0) then
        !
        j = 0 ; jp = j + 1
        !
        !--------------------------------
        !-- wall normal vector ------
        !---------------------------
        n1 = 0.d0; n2=1.d0;  n3=0.d0
        !----------------------------
        !-- wall velocity --
        !---------------------------
        uwall = 0.0d0; vwall = 0.0d0 ; wwall = 0.0d0
        !--------------------------------------------
        !-- wall temperature --
        !---------------------
        twall = 1.0d0   ! wall temperature is equal to the ref temperature
        !-------------------------------------
        !-- Accomodation coff
        !------------------------
        alpha = 1.0d0   ! diffusive wall
        if (moment=='r13') then
          do k=0,km
          do i=0,im
            ip = i ; kp = k
            call R13wbc(n1, n2, n3, i, j, k, ip, jp, kp, &
                        uwall, vwall, wwall, twall, alpha)
          end do
          end do
        elseif (moment=='r26') then
          do k=0,km
          do i=0,im
            ip = i ; kp = k
            call R26wbc(n1, n2, n3, i, j, k,ip, jp, kp, &
                        uwall, vwall, wwall, twall, alpha)
          end do
          end do
        end if
        !
      end if
      !
    else if(ndir==4) then
      !
      if (jrk==jrkm) then
        ! upper wall
        !
        j=jm ; jp = j - 1
        !
        !--------------------------------
        !-- wall normal vector ------
        !---------------------------
        n1 = 0.d0; n2=-1.d0;  n3=0.d0
        !----------------------------
        !-- wall velocity --
        !---------------------------
        uwall =  0.0d0; vwall = 0.0d0 ; wwall = 0.0d0
        !--------------------------------------------
        !-- wall temperature --
        !---------------------
        twall = 1.0d0   ! wall temperature is equal to the ref temperature
        !-------------------------------------
        !-- Accomodation coff
        !------------------------
        alpha = 1.0d0   ! diffusive wall
        !
        if (moment=='r13') then
          do k=0,km
          do i=0,im
            ip = i ; kp = k
            call R13wbc(n1, n2, n3, i, j, k, ip, jp, kp, &
                        uwall, vwall, wwall, twall,alpha)
          end do
          end do
        elseif(moment=='r26') then
          do k=0,km
          do i=0,im
            ip = i ; kp = k
            call R26wbc(n1, n2, n3, i, j, k, ip, jp, kp, &
                        uwall, vwall, wwall,twall,alpha)
           end do
           end do
        end if
        !
      end if
      !
    else
      stop ' !! ndir not defined @ MOM_wall_boundary'
    end if

    return

  end subroutine MOM_wall_boundary 
  !+-------------------------------------------------------------------+
  !| The end of the subroutine MOM_wall_boundary.                      |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to b.c. for N-S with R13 moment.               |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| xx-xx-xxxx: Created by X. Gu.                                     |
  !+-------------------------------------------------------------------+
  subroutine R13wbc(n1, n2, n3, iw, jw, kw, ip, jp, kp,  &
                         uwall, vwall, wwall, twall, alpha) 
    !===================================================================
    use commvar,   only : deltat,const2
    use commarray, only : q,rho, vel, sigma, qflux,tmp, prs, miu, x
    use fludyna,   only : fvar2q,thermal

    integer, intent(in) :: iw, jw, kw, ip, jp, kp
    real(8), intent(in) :: n1, n2, n3, uwall, vwall, wwall, twall, alpha
    !
    real(8) tempw, presw, u2w, v2w, w2w, uvw, uww, vww, qxw, qyw, qzw
    real(8) m111, m112, m113, m122, m223, m123, m222, m133, m233, m333
    real(8) r11, r22, r12, r13, r23, r33, deltap

    real(8) t1, t2, t3 , s1, s2 , s3
    real(8) n11, n22, n33, n12, n13, n23
    real(8) t11, t22, t33, t12, t13, t23
    real(8) s11, s22, s33, s12, s13, s23
    
    real(8) sig_nn, sig_tt, sig_ss, sig_nt, sig_ns, sig_ts
    real(8) r_nn, r_tt, r_ss, r_ts, r_nt, r_ns, qn, qt, qs
    real(8) m_ttt, m_sss, m_tts, m_tss, m_nts, m_ntt, m_nss, m_nnt, m_nns, m_nnn
    
    real(8) upb, vpb, wpb, t_jump, slip_t, slip_s, slipnorm, G, co1, co2
    real(8) slip2, term1, term2, tr1, tempp, rhop, presp, deltatp

    real(8) dx , dy , dz

    real(8) sig_nn_new, sig_tt_new, sig_ss_new, sig_ts_new

    deltatp = deltat 
    

    G = (2.0d0 - alpha)/alpha*sqrt(0.5d0*pi) 
    
    n11 = n1*n1;  n22 = n2*n2; n33 = n3*n3; n12 = n1*n2; n13 = n1*n3; n23 = n2*n3


    dx = x(iw, jw, kw, 1) - x(ip, jp, kp, 1)
    dy = x(iw, jw, kw, 2) - x(ip, jp, kp, 2)
    dz = x(iw, jw, kw, 3) - x(ip, jp, kp, 3)


    tempw = tmp(iw, jw, kw)
    presw  = prs(iw, jw, kw)
   
    u2w = sigma(iw, jw, kw, 1) 
    uvw = sigma(iw, jw, kw, 2) 
    uww = sigma(iw, jw, kw, 3) 
    v2w = sigma(iw, jw, kw, 4) 
    vww = sigma(iw, jw, kw, 5) 
    w2w = - u2w - v2w
   
    qxw = qflux(iw, jw, kw, 1) 
    qyw = qflux(iw, jw, kw, 2) 
    qzw = qflux(iw, jw, kw, 3) 

   
    m111 = mijk(iw,jw,kw,1)
    m112 = mijk(iw,jw,kw,2)
    m113 = mijk(iw,jw,kw,3)
    m122 = mijk(iw,jw,kw,4)
    m222 = mijk(iw,jw,kw,5)
    m223 = mijk(iw,jw,kw,6)
    m123 = mijk(iw,jw,kw,7)
    m133 = - m111 - m122 
    m233 = - m112 - m222
    m333 = - m113 - m223 

    r11 = rij(iw,jw,kw,1)
    r12 = rij(iw,jw,kw,2)
    r13 = rij(iw,jw,kw,3)
    r22 = rij(iw,jw,kw,4)
    r23 = rij(iw,jw,kw,5)
    r33 = - r11 - r22

    deltap = delta(iw, jw, kw)
   
    co2 = tempw/const2
    co1 = G*sqrt(co2)
    tr1 = twall/tempw

    upb = vel(iw, jw, kw,1) - uwall
    vpb = vel(iw, jw, kw,2) - vwall
    wpb = vel(iw, jw, kw,3) - wwall

    call slip_orientation(upb, vpb, wpb, n1,n2,n3,t1, t2, t3, s1, s2, s3,  &
                t11, t22, t33, t12, t13, t23, s11, s22, s33, s12, s13, s23)

    slip_t = (upb*t1 + vpb*t2 + wpb*t3)/sqrt(co2)
    slip_s = (upb*s1 + vpb*s2 + wpb*s3)/sqrt(co2)

    slip2 = slip_t*slip_t + slip_s*slip_s

    call trans_mijk_to_wall(n1, n2, n3, t1, t2, t3, s1, s2, s3,             &
                            n11, n22, n33, n12, n13, n23, t11, t22, t33,    &
                            t12, t13, t23, s11, s22, s33, s12, s13, s23,    &
                m111, m112, m113, m122, m223, m123, m222, m133, m233, m333, &
                m_ttt, m_sss, m_tts, m_tss, m_nts, m_ntt, m_nss, m_nnt, m_nns, m_nnn)

    call trans_sigma_to_wall(n1, n2, n3, t1, t2, t3, s1, s2, s3,            &
                             n11, n22, n33, n12, n13, n23, t11, t22, t33,   &
                             t12, t13, t23, s11, s22, s33, s12, s13, s23,   &
                             u2w, v2w, w2w, uvw, uww, vww,                  &
                             sig_nn, sig_tt, sig_ss, sig_nt, sig_ns, sig_ts)


    call trans_rij_to_wall(n1, n2, n3, t1, t2, t3, s1, s2, s3,             &
                           n11, n22, n33, n12, n13, n23, t11, t22, t33,    &
                           t12, t13, t23, s11, s22, s33, s12, s13, s23,    &
                           r11, r22, r12, r13, r23, r33,                   &
                           r_nn, r_tt, r_ss, r_ts, r_nt, r_ns)

    !=======================================
    !==          normal heat flux        ==
    !=======================================
    qn = n1*qxw + n2*qyw + n3*qzw
    !=======================================
    !             velocity slip
    !=======================================

    presw = presw + 0.5d0*sig_nn                   &
          - (30.0d0*r_nn + 7.0d0*deltap)/(840.0d0*co2)

    upb = n1*(u2w - term2 - 2.0d0*n23*vww)         &
        + (1.0d0 - 2.0d0*n11)*(n2*uvw + n3*uww) 
    
    vpb = n2*(v2w - term2 - 2.0d0*n13*uww)         &
        + (1.0d0 - 2.0d0*n22)*(n1*uvw + n3*vww)     

    wpb = n3*(w2w - term2 - 2.0d0*n12*uvw)         &
         + (1.0d0 - 2.0d0*n33)*(n1*uww + n2*vww) 
    
    upb = -upb*co1 - 0.2d0*(qxw - qn*n1)  
    vpb = -vpb*co1 - 0.2d0*(qyw - qn*n2)  
    wpb = -wpb*co1 - 0.2d0*(qzw - qn*n3) 

    !==========================================
    call slip_orientation(upb, vpb, wpb, n1,n2,n3,t1, t2, t3, s1, s2, s3,  &
                t11, t22, t33, t12, t13, t23, s11, s22, s33, s12, s13, s23)

    upb = (upb - 0.5d0*(t1*m_nnt + s1*m_nns))/presw
    vpb = (vpb - 0.5d0*(t2*m_nnt + s2*m_nns))/presw
    wpb = (wpb - 0.5d0*(t3*m_nnt + s3*m_nns))/presw
        
    upb = uwall + upb
    vpb = vwall + vpb
    wpb = wwall + wpb

    vel(iw, jw, kw, 1) = vel(iw,jw,kw,1) + (upb-vel(iw,jw,kw,1))*deltatp 
    vel(iw, jw, kw, 2) = vel(iw,jw,kw,2) + (vpb-vel(iw,jw,kw,2))*deltatp
    vel(iw, jw, kw, 3) = vel(iw,jw,kw,3) + (wpb-vel(iw,jw,kw,3))*deltatp
   
    !========================
    !==   temperature jump ==
    !========================
    t_jump = -0.5d0*co1*qn - (deltap/30.d0 + 5.0d0/56.0d0*r_nn) 
    t_jump = t_jump*const2/presw + 0.25d0*(slip2*const2 - tempw*sig_nn/presw)

    tempp = twall + t_jump 
    
    tmp(iw, jw, kw) = tmp(iw, jw, kw) + (tempp - tmp(iw, jw, kw) )*deltatp  
    !==================
    !== pressure ======
    !==================
    prs(iw, jw, kw) = prs(iw, jw, kw) 

    !==================
    !== density =======
    !==================
    rho(iw,jw,kw) = thermal(pressure=prs(iw,jw,kw),temperature=tmp(iw,jw,kw))

    call fvar2q(      q=  q(iw,jw,kw,:),                            &
                   density=rho(iw,jw,kw),                              &
                  velocity=vel(iw,jw,kw,:),                            &
                  temperature=tmp(iw,jw,kw) )
    !==============
    !== stresses == 
    !==============
    sig_tt_new = (-co1*(m_ntt + 0.4d0*qn) + r_ss/14.0d0 - deltap/30.0d0)/co2  & 
           + (tr1 - 1.0d0 + slip_t**2)*presw  
          
    sig_nn_new =  (-co1*(0.5d0*m_nnn + 0.6d0*qn) - (r_nn/7.0d0 + deltap/30.0d0) )/co2  &
           + (tr1 - 1.0d0)*presw
 
    sig_ss_new = - sig_nn_new - sig_tt_new

    sig_ts_new = (-co1*m_nts - r_ts/14.d0)/co2 + slip_t*slip_s*presw

    !====================
    ! == heat fluxes =====
    !====================
    term1 = (tr1 + slip2/6.0d0 )*presw/0.6d0*sqrt(co2)
    qt = - co1*(r_nt/co2 + 7.0d0*sig_nt)/3.6d0  - m_nnt/0.9d0 - slip_t*term1
    qs = - co1*(r_ns/co2 + 7.0d0*sig_ns)/3.6d0  - m_nns/0.9d0 - slip_s*term1

    u2w = sigma(ip,jp,kp,1) + dx*dsigmaxx(ip,jp,kp,1) + dy*dsigmaxx(ip,jp,kp,2) &
                            + dz*dsigmaxx(ip,jp,kp,3)
    uvw = sigma(ip,jp,kp,2) + dx*dsigmaxy(ip,jp,kp,1) + dy*dsigmaxy(ip,jp,kp,2) &
                            + dz*dsigmaxy(ip,jp,kp,3)
    uww = sigma(ip,jp,kp,3) + dx*dsigmaxz(ip,jp,kp,1) + dy*dsigmaxz(ip,jp,kp,2) &
                            + dz*dsigmaxz(ip,jp,kp,3)
    v2w = sigma(ip,jp,kp,4) + dx*dsigmayy(ip,jp,kp,1) + dy*dsigmayy(ip,jp,kp,2) &
                            + dz*dsigmayy(ip,jp,kp,3)
    vww = sigma(ip,jp,kp,5) + dx*dsigmayz(ip,jp,kp,1) + dy*dsigmayz(ip,jp,kp,2) &
                            + dz*dsigmayz(ip,jp,kp,3)

    w2w = - u2w - v2w

    call trans_sigma_to_wall(n1, n2, n3, t1, t2, t3, s1, s2, s3, &
                             n11, n22, n33, n12, n13, n23, t11, t22, t33, &
                             t12, t13, t23, s11, s22, s33, s12, s13, s23, &
                             u2w, v2w, w2w, uvw, uww, vww,                 &
                             sig_nn, sig_tt, sig_ss, sig_nt, sig_ns, sig_ts)


    sig_tt = sig_tt_new
    sig_nn = sig_nn_new
    sig_ss = sig_ss_new
    sig_ts = sig_ts_new

    call trans_sigma_to_cartisian(n1, n2, n3, t1, t2, t3, s1, s2, s3,             &
                                  n11, n22, n33, n12, n13, n23, t11, t22, t33,    &
                                  t12, t13, t23, s11, s22, s33, s12, s13, s23,    &
                                  u2w, v2w, w2w, uvw, uww, vww,                   &
                                  sig_nn, sig_tt, sig_ss, sig_nt, sig_ns, sig_ts)
                                  
    sigma(iw,jw,kw,1) = sigma(iw,jw,kw,1) + (u2w-sigma(iw,jw,kw,1))*deltatp 
    sigma(iw,jw,kw,2) = sigma(iw,jw,kw,2) + (uvw-sigma(iw,jw,kw,2))*deltatp
    sigma(iw,jw,kw,3) = sigma(iw,jw,kw,3) + (uww-sigma(iw,jw,kw,3))*deltatp
    sigma(iw,jw,kw,4) = sigma(iw,jw,kw,4) + (v2w-sigma(iw,jw,kw,4))*deltatp
    sigma(iw,jw,kw,5) = sigma(iw,jw,kw,5) + (vww-sigma(iw,jw,kw,5))*deltatp 
    !==================================

    qxw = qflux(ip,jp,kp,1) + dx*dqx(ip,jp,kp,1) + dy*dqx(ip,jp,kp,2) &
                            + dz*dqx(ip,jp,kp,3)
    qyw = qflux(ip,jp,kp,2) + dx*dqy(ip,jp,kp,1) + dy*dqy(ip,jp,kp,2) &
                            + dz*dqy(ip,jp,kp,3)
    qzw = qflux(ip,jp,kp,3) + dx*dqz(ip,jp,kp,1) + dy*dqz(ip,jp,kp,2) &
                               + dz*dqz(ip,jp,kp, 3)

    qn = n1*qxw + n2*qyw + n3*qzw
      
    qxw = n1*qn + t1*qt + s1*qs
    qyw = n2*qn + t2*qt + s2*qs
    qzw = n3*qn + t3*qt + s3*qs

    qflux(iw,jw,kw,1) = qflux(iw,jw,kw,1) + (qxw - qflux(iw,jw,kw,1))*deltatp
    qflux(iw,jw,kw,2) = qflux(iw,jw,kw,2) + (qyw - qflux(iw,jw,kw,2))*deltatp
    qflux(iw,jw,kw,3) = qflux(iw,jw,kw,3) + (qzw - qflux(iw,jw,kw,3))*deltatp 

    call fvar2qmom(       q=q_mom(iw, jw, kw, :),               &
                      stress=sigma(iw, jw, kw,:),               &
                    heatflux=qflux(iw, jw, kw,:)                )

   
    return

  end  subroutine R13wbc
  !+-------------------------------------------------------------------+
  !| The end of the subroutine R13wbc.                                 |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to b.c. for N-S with R26 moment.               |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| xx-xx-xxxx: Created by X. Gu.                                     |
  !+-------------------------------------------------------------------+
  subroutine R26wbc(n1, n2, n3, iw, jw, kw, ip, jp, kp,  &
                  uwall, vwall, wwall, twall, alpha) 

    !==================================================================
       
    use commvar,   only : deltat,const2
    use commarray, only : q, rho, vel, sigma, qflux, tmp, prs, miu, spc, x
    use fludyna,   only : fvar2q, thermal


    integer, intent(in) :: iw, jw, kw, ip, jp, kp
    real(8), intent(in) :: n1, n2, n3, uwall, vwall, wwall, twall, alpha
    
    real(8) tempw, presw, u2w, v2w, w2w, uvw, uww, vww, qxw, qyw, qzw
    real(8) m111, m112, m113, m122, m223, m123, m222, m133, m233, m333
    real(8) r11, r22, r12, r13, r23, r33, deltap
    real(8) q1111, q1112, q1122, q1222, q2222, q1113,  &
            q1123, q1223, q2223, q1333, q2333, q1133,  &
            q2233, q1233, q3333
            
    real(8) psi111, psi222, psi112, psi122, psi123,   &
            psi113, psi223, psi133, psi233, psi333

    real(8) omega_x, omega_y, omega_z

    real(8) t1, t2, t3 , s1, s2 , s3
    real(8) n11, n22, n33, n12, n13, n23
    real(8) t11, t22, t33, t12, t13, t23
    real(8) s11, s22, s33, s12, s13, s23
    
    real(8) sig_nn, sig_tt, sig_ss, sig_nt, sig_ns, sig_ts
    real(8) r_nn, r_tt, r_ss, r_ts, r_nt, r_ns, qn, qt, qs
    real(8) m_ttt, m_sss, m_tts, m_tss, m_nts, m_ntt, m_nss, m_nnt, m_nns, m_nnn
    real(8) q_tttt, q_ttts, q_ttss, q_tsss, q_ssss,   &
            q_nttt, q_ntts, q_ntss, q_nsss, q_nnnt,   &
            q_nnns, q_nntt, q_nnss, q_nnts, q_nnnn

    real(8) psi_ttt, psi_sss, psi_tts, psi_tss, psi_nts,   &
            psi_ntt, psi_nss, psi_nnt, psi_nns, psi_nnn
    
    real(8) omega_n, omega_t, omega_s
    
    real(8) upb, vpb, wpb, t_jump, slip_t, slip_s, G, co1, co2, co3
    real(8) slip2, term1, term2, tr1, tempp, rhop, presp, deltatp, deltatt

    real(8) dx, dy , dz

    real(8) qt_new, qs_new
    real(8) sig_nn_new, sig_tt_new, sig_ss_new, sig_ts_new

    real(8) r_nn_new, r_tt_new, r_ss_new, r_ts_new

    real(8) m_ttt_new, m_sss_new, m_nnt_new, m_nns_new, m_nts_new,  &
            m_tts_new, m_tss_new

    deltatp = deltat*0.01d0 
    deltatt = deltat*0.01d0
    
    G = (2.0d0 - alpha)/alpha*sqrt(0.5d0*pi) 
    
    n11 = n1*n1;  n22 = n2*n2; n33 = n3*n3; n12 = n1*n2; n13 = n1*n3; n23 = n2*n3

    upb = vel(iw, jw, kw, 1) - uwall
    vpb = vel(iw, jw, kw, 2) - vwall
    wpb = vel(iw, jw, kw, 3) - wwall
 
    dx = x(iw, jw, kw, 1) - x(ip, jp, kp, 1)
    dy = x(iw, jw, kw, 2) - x(ip, jp, kp, 2)
    dz = x(iw, jw, kw, 3) - x(ip, jp, kp, 3)

    tempw = tmp(iw, jw, kw)
    presw = prs(iw, jw, kw)
   
    u2w = sigma(iw, jw, kw, 1) 
    uvw = sigma(iw, jw, kw, 2) 
    uww = sigma(iw, jw, kw, 3) 
    v2w = sigma(iw, jw, kw, 4) 
    vww = sigma(iw, jw, kw, 5) 
    w2w = - u2w - v2w

    qxw = qflux(iw, jw, kw, 1) 
    qyw = qflux(iw, jw, kw, 2) 
    qzw = qflux(iw, jw, kw, 3)

    m111 = mijk(iw,jw,kw,1) 
    m112 = mijk(iw,jw,kw,2) 
    m113 = mijk(iw,jw,kw,3) 
    m122 = mijk(iw,jw,kw,4) 
    m222 = mijk(iw,jw,kw,5) 
    m223 = mijk(iw,jw,kw,6) 
    m123 = mijk(iw,jw,kw,7) 
    m133 = - m111 - m122 
    m233 = - m112 - m222
    m333 = - m113 - m223 

    r11 = rij(iw,jw,kw,1) 
    r12 = rij(iw,jw,kw,2) 
    r13 = rij(iw,jw,kw,3) 
    r22 = rij(iw,jw,kw,4) 
    r23 = rij(iw,jw,kw,5) 
    r33 = - r11 - r22
    
    q1111 = Eijkl(iw,jw,kw,1)
    q2222 = Eijkl(iw,jw,kw,8)
    q1112 = Eijkl(iw,jw,kw,2)
    q1113 = Eijkl(iw,jw,kw,3)
    q1122 = Eijkl(iw,jw,kw,4)
    q1222 = Eijkl(iw,jw,kw,6)
    q1123 = Eijkl(iw,jw,kw,5)
    q1223 = Eijkl(iw,jw,kw,7)
    q2223 = Eijkl(iw,jw,kw,9)

    q1333 = - q1113 - q1223
    q2333 = - q1123 - q2223
    q1133 = - q1111 - q1122
    q2233 = - q1122 - q2222
    q1233 = - q1112 - q1222
    q3333 = - q1133 - q2233
    
    psi111 = psi(iw,jw,kw,1)
    psi112 = psi(iw,jw,kw,2)
    psi113 = psi(iw,jw,kw,3)
    psi122 = psi(iw,jw,kw,4)
    psi222 = psi(iw,jw,kw,5)
    psi223 = psi(iw,jw,kw,6)
    psi123 = psi(iw,jw,kw,7)
    psi133 = - psi111 - psi122 
    psi233 = - psi112 - psi222
    psi333 = - psi113 - psi223 

    omega_x = omega(iw,jw,kw,1)
    omega_y = omega(iw,jw,kw,2)
    omega_z = omega(iw,jw,kw,3)

    deltap = delta(iw, jw, kw)
   

    call slip_orientation(upb, vpb, wpb, n1,n2,n3,t1, t2, t3, s1, s2, s3,  &
                t11, t22, t33, t12, t13, t23, s11, s22, s33, s12, s13, s23)

    co2 = tempw/const2
    co1 = G*sqrt(co2)
    tr1 = twall/tempw
     
    slip_t = (upb*t1 + vpb*t2 + wpb*t3)/sqrt(co2)
    slip_s = (upb*s1 + vpb*s2 + wpb*s3)/sqrt(co2)
    slip2 = slip_t*slip_t + slip_s*slip_s

    call trans_sigma_to_wall(n1, n2, n3, t1, t2, t3, s1, s2, s3,           &
                             n11, n22, n33, n12, n13, n23, t11, t22, t33,  &
                             t12, t13, t23, s11, s22, s33, s12, s13, s23,  &
                             u2w, v2w, w2w, uvw, uww, vww,                 &
                             sig_nn, sig_tt, sig_ss, sig_nt, sig_ns, sig_ts)

    call trans_rij_to_wall(n1, n2, n3, t1, t2, t3, s1, s2, s3,           &
                           n11, n22, n33, n12, n13, n23, t11, t22, t33,  &
                           t12, t13, t23, s11, s22, s33, s12, s13, s23,  &
                           r11, r22, r12, r13, r23, r33,                 &
                           r_nn, r_tt, r_ss, r_ts, r_nt, r_ns)

    call trans_mijk_to_wall(n1, n2, n3, t1, t2, t3, s1, s2, s3,          &
                            n11, n22, n33, n12, n13, n23, t11, t22, t33, &
                            t12, t13, t23, s11, s22, s33, s12, s13, s23, &
                m111, m112, m113, m122, m223, m123, m222, m133, m233, m333, &
                m_ttt, m_sss, m_tts, m_tss, m_nts, m_ntt, m_nss, m_nnt, m_nns, m_nnn)

    call trans_phi_to_wall(n1, n2, n3, t1, t2, t3, s1, s2, s3,           &
                         n11, n22, n33, n12, n13, n23, t11, t22, t33,    &
                         t12, t13, t23, s11, s22, s33, s12, s13, s23,    &
                         q1111, q1112, q1122, q1222, q2222, q1113,       &
                         q1123, q1223, q2223, q1333, q2333, q1133,       &
                         q2233, q1233, q3333,                            &
                         q_tttt, q_ttts, q_ttss, q_tsss, q_ssss,         &
                         q_nttt, q_ntts, q_ntss, q_nsss, q_nnnt,         &
                         q_nnns, q_nntt, q_nnss, q_nnts, q_nnnn)

    call trans_psi_to_wall(n1, n2, n3, t1, t2, t3, s1, s2, s3,           &
                         n11, n22, n33, n12, n13, n23, t11, t22, t33,    &
                         t12, t13, t23, s11, s22, s33, s12, s13, s23,    &
                         psi111, psi112, psi113, psi122, psi223, psi123, &
                         psi222, psi133, psi233, psi333,                 &
                         psi_ttt, psi_sss, psi_tts, psi_tss, psi_nts,    &
                         psi_ntt, psi_nss, psi_nnt, psi_nns, psi_nnn)

    omega_t = t1*omega_x + t2*omega_y + t3*omega_z
    omega_s = s1*omega_x + s2*omega_y + s3*omega_z
    omega_n = n1*omega_x + n2*omega_y + n3*omega_z

    qt = t1*qxw + t2*qyw + t3*qzw
    qs = s1*qxw + s2*qyw + s3*qzw
    qn = n1*qxw + n2*qyw + n3*qzw
       
    !=======================================
    !             velocity slip
    !=======================================

    presw = presw + 0.5d0*sig_nn                              &
           - ( (30.0d0*r_nn + 7.0d0*deltap)/840.d0 + q_nnnn/24.0d0)/co2

    upb = n1*(u2w - term2 - 2.0d0*n23*vww)         &
        + (1.0d0 - 2.0d0*n11)*(n2*uvw + n3*uww) 
    
    vpb = n2*(v2w - term2 - 2.0d0*n13*uww)         &
        + (1.0d0 - 2.0d0*n22)*(n1*uvw + n3*vww)     

    wpb = n3*(w2w - term2 - 2.0d0*n12*uvw)         &
         + (1.0d0 - 2.0d0*n33)*(n1*uww + n2*vww) 
    
    co3 = 9.0d0/2520.0d0/co2
    upb = -upb*co1 - 0.2d0*(qxw - qn*n1) + co3*(omega_x - omega_n*n1)
    vpb = -vpb*co1 - 0.2d0*(qyw - qn*n2) + co3*(omega_y - omega_n*n2)
    wpb = -wpb*co1 - 0.2d0*(qzw - qn*n3) + co3*(omega_z - omega_n*n3)

    !==========================================

    call slip_orientation(upb, vpb, wpb, n1,n2,n3,t1, t2, t3, s1, s2, s3,  &
                t11, t22, t33, t12, t13, t23, s11, s22, s33, s12, s13, s23)

    co3 = 7.0d0/252.0d0/co2
    upb = (upb - 0.5d0*(t1*m_nnt + s1*m_nns) + co3*(t1*psi_nnt + s1*psi_nns))/presw
    vpb = (vpb - 0.5d0*(t2*m_nnt + s2*m_nns) + co3*(t2*psi_nnt + s2*psi_nns))/presw
    wpb = (wpb - 0.5d0*(t3*m_nnt + s3*m_nns) + co3*(t3*psi_nnt + s3*psi_nns))/presw

    upb = uwall + upb
    vpb = vwall + vpb
    wpb = wwall + wpb

    vel(iw, jw, kw, 1) = vel(iw,jw,kw,1) + (upb-vel(iw,jw,kw,1))*deltatp 
    vel(iw, jw, kw, 2) = vel(iw,jw,kw,2) + (vpb-vel(iw,jw,kw,2))*deltatp
    vel(iw, jw, kw, 3) = vel(iw,jw,kw,3) + (wpb-vel(iw,jw,kw,3))*deltatp

    !========================
    !==   temperature jump ==
    !========================
    t_jump = -0.5d0*co1*qn - (deltap/30.d0 + 5.0d0/56.0d0*r_nn) + q_nnnn/24.0d0
               
    t_jump = t_jump*const2/presw + 0.25d0*(slip2*const2 - tempw*sig_nn/presw)

    tempp = twall + t_jump 
    
    tmp(iw, jw, kw) = tmp(iw, jw, kw) + (tempp - tmp(iw, jw, kw) )*deltatt  

    !==================
    !== pressure ======
    !==================
    prs(iw,jw,kw) = prs(ip,jp,kp)

    !==================
    !== density =======
    !==================
    rho(iw,jw,kw) =thermal(pressure=prs(iw,jw,kw),temperature=tmp(iw,jw,kw))

    call fvar2q(      q=  q(iw,jw,kw,:),         &
                   density=rho(iw,jw,kw),        &
                  velocity=vel(iw,jw,kw,:),      &
                  temperature=tmp(iw,jw,kw)        ) 

    !==============
    !== stresses == 
    !==============
    sig_tt_new = (-co1*(m_ntt + 0.4d0*qn) + r_ss/14.0d0 - deltap/30.0d0 -0.5d0*q_nntt)/co2  & 
           + (tr1 - 1.0d0 + slip_t**2)*presw  
          
    sig_nn_new =  (-co1*(0.5d0*m_nnn + 0.6d0*qn)               &
           - (r_nn/7.0d0 + deltap/30.0d0 + q_nnnn/6.d0) )/co2  &
           + (tr1 - 1.0d0)*presw
 
    sig_ss_new = - sig_nn_new - sig_tt_new

    sig_ts_new = (-co1*m_nts - r_ts/14.d0 - 0.5d0*q_nnts)/co2 + slip_t*slip_s*presw

    !====================
    ! == heat fluxes =====
    !====================
    term1 = (tr1 + slip2/6.0d0 )*presw/0.6d0*sqrt(co2)
    qt_new = - co1*(r_nt/co2 + 7.0d0*sig_nt)/3.6d0  - m_nnt/0.9d0  &
           - slip_t*term1  - (psi_nnt/16.2d0 + omega_t/56.0d0)/co2
    qs_new = - co1*(r_ns/co2 + 7.0d0*sig_ns)/3.6d0  - m_nns/0.9d0  &
           - slip_s*term1  - (psi_nns/16.2d0 + omega_s/56.0d0)/co2

    !===============
    !==  mijk ======
    !===============
    term1 = (3.0d0*tr1 + slip2)*presw*sqrt(co2)

    m_ttt_new = -co1*(3.0d0*sig_nt + (q_nttt + 3.0d0/7.0d0*r_nt)/co2) &
      - slip_t*term1 - 1.8d0*qt - 1.5d0*m_nnt                      &
      -(psi_nnt/12.0d0 + psi_ttt/18.0 + 9.0d0*omega_t/280.0d0)/co2

    m_sss_new = -co1*(3.0d0*sig_ns + (q_nsss + 3.0d0/7.0d0*r_ns)/co2) &
      - slip_s*term1 - 1.8d0*qs - 1.5d0*m_nns                      &
      -(psi_nns/12.0d0 + psi_sss/18.0 + 9.0d0*omega_s/280.0d0)/co2

    m_nnt_new = -co1*(sig_nt + (q_nnnt/3.0d0 + r_nt/7.0d0)/co2)       &
      - slip_t*tr1*presw*sqrt(co2)/1.5d0 - 0.4d0*qt                &
      -(omega_t/140.0d0 + psi_nnt/18.0d0)/co2

    m_nns_new = -co1*(sig_ns + (q_nnns/3.0d0 + r_ns/7.0d0)/co2)       &
      - tr1*slip_s*presw*sqrt(co2)/1.5d0 - 0.4d0*qs               &
      -(omega_s/140.0d0 + psi_nns/18.0d0)/co2

    m_nts_new = -0.5d0*co1*(sig_ts + (q_nnts + r_ts/7.0d0)/co2) - psi_nts/18.0d0/co2

    m_tts_new = -(m_sss_new + m_nns_new)
    m_tss_new = -(m_ttt_new + m_nnt_new)

    !===============
    !==  Rij ======
    !===============
    r_tt_new = -co1*(28.0d0/15.0d0*qn + 14.0d0/3.0d0*m_ntt    &
         + (omega_n/15.0d0 +14.0d0/27.0d0*psi_ntt)/co2)   &
         + (7.0d0*presw *(slip_t**4 + 6.0d0*tr1*slip_t**2  &
         + 3.0d0*tr1**2 - 3.0d0)/9.0d0 - 14.0d0/3.0d0*sig_tt)*co2  &
         - R_nn/3.0d0 - 14.0d0/45.0d0*deltap               &
         - 7.0d0/9.0d0*(q_tttt + 3.0d0*q_nntt)


    r_nn_new = -co1*(21.0d0/8.0d0*qn + 35.0d0/16.0d0*m_nnn               &
         + (35.0d0/144.0d0*psi_nnn + 3.0d0/32.0d0*omega_n)/co2 )     &
         + (7.0d0/4.0d0*presw*(tr1*tr1 - 1.0d0) - 3.5d0*sig_nn)*co2 &
         - 7.0d0*deltap/30.0d0 - 7.0d0/6.0d0*q_nnnn


    r_ts_new = - co1*(21.0d0*m_nts + 7.0d0/3.0d0*psi_nts/co2)/4.0d0      &
        + 7.0d0/3.0d0*(slip_t*slip_s*(2.0d0*tr1 + 0.25*slip2)*presw  &
        - 2.0d0*sig_ns)*co2 - 35.0d0/12.0d0*q_nnts


    r_ss_new = - r_nn_new - r_tt_new

    !============
    !== Delta ==
    !===========

    deltap = -5.0d0/16.0d0*co1*(28.0d0*qn + Omega_n/co2)               &
    + 1.25d0*co2*( (6.0d0*(tr1*tr1 - 1.0d0)                                 &
    + slip2*(3.0d0*tr1*tr1 + 0.25d0*slip2) )*presw - 3.0d0*sig_nn)  &
    - 1.875d0*r_nn + 35.0d0/48.0d0*q_nnnn

    delta(iw, jw, kw) = delta(iw, jw, kw) + (deltap - delta(iw, jw, kw))*deltatt

    !=====================
    !== extrapolation ==
    !=====================


    u2w = sigma(ip, jp, kp, 1) + dx*dsigmaxx(ip, jp, kp, 1) &
                               + dy*dsigmaxx(ip, jp, kp, 2) &
                               + dz*dsigmaxx(ip, jp, kp, 3)
    uvw = sigma(ip, jp, kp, 2) + dx*dsigmaxy(ip, jp, kp, 1) &
                               + dy*dsigmaxy(ip, jp, kp, 2) &
                               + dz*dsigmaxy(ip, jp, kp, 3)
    uww = sigma(ip, jp, kp, 3) + dx*dsigmaxz(ip, jp, kp, 1) &
                               + dy*dsigmaxz(ip, jp, kp, 2) &
                               + dz*dsigmaxz(ip, jp, kp, 3)
    v2w = sigma(ip, jp, kp, 4) + dx*dsigmayy(ip, jp, kp, 1) &
                               + dy*dsigmayy(ip, jp, kp, 2) &
                               + dz*dsigmayy(ip, jp, kp, 3)
    vww = sigma(ip, jp, kp, 5) + dx*dsigmayz(ip, jp, kp, 1) &
                               + dy*dsigmayz(ip, jp, kp, 2) &
                               + dz*dsigmayz(ip, jp, kp, 3)

    w2w = - u2w - v2w

    call trans_sigma_to_wall(n1, n2, n3, t1, t2, t3, s1, s2, s3,           &
                             n11, n22, n33, n12, n13, n23, t11, t22, t33,  &
                             t12, t13, t23, s11, s22, s33, s12, s13, s23,  &
                             u2w, v2w, w2w, uvw, uww, vww,                 &
                             sig_nn, sig_tt, sig_ss, sig_nt, sig_ns, sig_ts)


    sig_tt = sig_tt_new 
    sig_nn = sig_nn_new 
    sig_ss = - (sig_nn + sig_tt)
    sig_ts = sig_ts_new 

    call trans_sigma_to_cartisian(n1, n2, n3, t1, t2, t3, s1, s2, s3,          &
                                  n11, n22, n33, n12, n13, n23, t11, t22, t33, &
                                  t12, t13, t23, s11, s22, s33, s12, s13, s23, &
                                  u2w, v2w, w2w, uvw, uww, vww,                &
                                  sig_nn, sig_tt, sig_ss, sig_nt, sig_ns, sig_ts)
                                  
    sigma(iw,jw,kw,1) = sigma(iw,jw,kw,1) + (u2w-sigma(iw,jw,kw,1))*deltatp 
    sigma(iw,jw,kw,2) = sigma(iw,jw,kw,2) + (uvw-sigma(iw,jw,kw,2))*deltatp
    sigma(iw,jw,kw,3) = sigma(iw,jw,kw,3) + (uww-sigma(iw,jw,kw,3))*deltatp
    sigma(iw,jw,kw,4) = sigma(iw,jw,kw,4) + (v2w-sigma(iw,jw,kw,4))*deltatp
    sigma(iw,jw,kw,5) = sigma(iw,jw,kw,5) + (vww-sigma(iw,jw,kw,5))*deltatp 
    !==================================

    qxw = qflux(ip, jp, kp, 1) + dx*dqx(ip, jp, kp, 1) &
                               + dy*dqx(ip, jp, kp, 2) &
                               + dz*dqx(ip, jp, kp, 3)

    qyw = qflux(ip, jp, kp, 2) + dx*dqy(ip, jp, kp, 1) &
                               + dy*dqy(ip, jp, kp, 2) &
                               + dz*dqy(ip, jp, kp, 3)

    qzw = qflux(ip, jp, kp, 3) + dx*dqz(ip, jp, kp, 1) &
                               + dy*dqz(ip, jp, kp, 2) &
                               + dz*dqz(ip, jp, kp, 3)

    qn = n1*qxw + n2*qyw + n3*qzw
    qt = qt_new
    qs = qs_new

    qxw = n1*qn + t1*qt + s1*qs
    qyw = n2*qn + t2*qt + s2*qs
    qzw = n3*qn + t3*qt + s3*qs

    qflux(iw,jw,kw,1) = qflux(iw,jw,kw,1) + (qxw - qflux(iw,jw,kw,1))*deltatt
    qflux(iw,jw,kw,2) = qflux(iw,jw,kw,2) + (qyw - qflux(iw,jw,kw,2))*deltatt
    qflux(iw,jw,kw,3) = qflux(iw,jw,kw,3) + (qzw - qflux(iw,jw,kw,3))*deltatt 
    
    !===============
    !==  mijk ======
    !===============
    m111 = mijk(ip,jp,kp,1) + dx*dmijk(ip,jp,kp,1,1)   &
                            + dy*dmijk(ip,jp,kp,2,1)   &
                            + dz*dmijk(ip,jp,kp,3,1)
    m112 = mijk(ip,jp,kp,2) + dx*dmijk(ip,jp,kp,1,2)   &
                            + dy*dmijk(ip,jp,kp,2,2)   &
                            + dz*dmijk(ip,jp,kp,3,2)
    m113 = mijk(ip,jp,kp,3) + dx*dmijk(ip,jp,kp,1,3)   &
                            + dy*dmijk(ip,jp,kp,2,3)   &
                            + dz*dmijk(ip,jp,kp,3,3)
    m122 = mijk(ip,jp,kp,4) + dx*dmijk(ip,jp,kp,1,4)   &
                            + dy*dmijk(ip,jp,kp,2,4)   &
                            + dz*dmijk(ip,jp,kp,3,4)
    m222 = mijk(ip,jp,kp,5) + dx*dmijk(ip,jp,kp,1,5)   &
                            + dy*dmijk(ip,jp,kp,2,5)   &
                            + dz*dmijk(ip,jp,kp,3,5)
    m223 = mijk(ip,jp,kp,6) + dx*dmijk(ip,jp,kp,1,6)   &
                            + dy*dmijk(ip,jp,kp,2,6)   &
                            + dz*dmijk(ip,jp,kp,3,6)
    m123 = mijk(ip,jp,kp,7) + dx*dmijk(ip,jp,kp,1,7)   &
                            + dy*dmijk(ip,jp,kp,2,7)   &
                            + dz*dmijk(ip,jp,kp,3,7)

    m133 = - m111 - m122
    m233 = - m112 - m222
    m333 = - m113 - m223

    call trans_mijk_to_wall(n1, n2, n3, t1, t2, t3, s1, s2, s3,          &
                            n11, n22, n33, n12, n13, n23, t11, t22, t33, &
                            t12, t13, t23, s11, s22, s33, s12, s13, s23, &
                m111, m112, m113, m122, m223, m123, m222, m133, m233, m333, &
                m_ttt, m_sss, m_tts, m_tss, m_nts, m_ntt, m_nss, m_nnt, m_nns, m_nnn)

    m_ttt = m_ttt_new 
    m_sss = m_sss_new 
    m_nnt = m_nnt_new 
    m_nns = m_nns_new 
    m_nts = m_nts_new 
    m_tts = -(m_sss + m_nns)
    m_tss = -(m_ttt + m_nnt)

    call trans_mijk_to_cartisian(n1, n2, n3, t1, t2, t3, s1, s2, s3,        &
                            n11, n22, n33, n12, n13, n23, t11, t22, t33,    &
                            t12, t13, t23, s11, s22, s33, s12, s13, s23,    &
                m111, m112, m113, m122, m223, m123, m222, m133, m233, m333, &
       m_ttt, m_sss, m_tts, m_tss, m_nts, m_ntt, m_nss, m_nnt, m_nns, m_nnn)
       
    mijk(iw,jw,kw,1) = mijk(iw,jw,kw,1) + (m111 - mijk(iw,jw,kw,1))*deltatp
    mijk(iw,jw,kw,2) = mijk(iw,jw,kw,2) + (m112 - mijk(iw,jw,kw,2))*deltatp
    mijk(iw,jw,kw,3) = mijk(iw,jw,kw,3) + (m113 - mijk(iw,jw,kw,3))*deltatp
    mijk(iw,jw,kw,4) = mijk(iw,jw,kw,4) + (m122 - mijk(iw,jw,kw,4))*deltatp
    mijk(iw,jw,kw,5) = mijk(iw,jw,kw,5) + (m222 - mijk(iw,jw,kw,5))*deltatp
    mijk(iw,jw,kw,6) = mijk(iw,jw,kw,6) + (m223 - mijk(iw,jw,kw,6))*deltatp
    mijk(iw,jw,kw,7) = mijk(iw,jw,kw,7) + (m123 - mijk(iw,jw,kw,7))*deltatp    
       

    r11 = rij(ip,jp,kp,1) + dx * dRij(ip,jp,kp,1,1)  &
                          + dy * dRij(ip,jp,kp,2,1)  &
                          + dz * dRij(ip,jp,kp,3,1)
    r12 = rij(ip,jp,kp,2) + dx * dRij(ip,jp,kp,1,2)  &
                          + dy * dRij(ip,jp,kp,2,2)  &
                          + dz * dRij(ip,jp,kp,3,2)

    r13 = rij(ip,jp,kp,3) + dx * dRij(ip,jp,kp,1,3)  &
                          + dy * dRij(ip,jp,kp,2,3)  &
                          + dz * dRij(ip,jp,kp,3,3)

    r22 = rij(ip,jp,kp,4) + dx * dRij(ip,jp,kp,1,4)  &
                          + dy * dRij(ip,jp,kp,2,4)  &
                          + dz * dRij(ip,jp,kp,3,4)

    r23 = rij(ip,jp,kp,5) + dx * dRij(ip,jp,kp,1,5)  &
                          + dy * dRij(ip,jp,kp,2,5)  &
                          + dz * dRij(ip,jp,kp,3,5)

    r33 = - r11 - r22


    call trans_rij_to_wall(n1, n2, n3, t1, t2, t3, s1, s2, s3,           &
                           n11, n22, n33, n12, n13, n23, t11, t22, t33,  &
                           t12, t13, t23, s11, s22, s33, s12, s13, s23,  &
                           r11, r22, r12, r13, r23, r33,                 &
                           r_nn, r_tt, r_ss, r_ts, r_nt, r_ns)


    r_tt = r_tt_new  
    r_nn = r_nn_new 
    r_ts = r_ts_new 
    r_ss = r_ss_new 
    
    call trans_rij_to_cartisian(n1, n2, n3, t1, t2, t3, s1, s2, s3,    &
                             n11, n22, n33, n12, n13, n23, t11, t22, t33,    &
                             t12, t13, t23, s11, s22, s33, s12, s13, s23,    &
                             r11, r22, r12, r13, r23, r33,                   &
                             r_nn, r_tt, r_ss, r_ts, r_nt, r_ns)
                             
    rij(iw,jw,kw,1) = rij(iw,jw,kw,1) + (r11 - rij(iw,jw,kw,1))*deltatt
    rij(iw,jw,kw,2) = rij(iw,jw,kw,2) + (r12 - rij(iw,jw,kw,2))*deltatt
    rij(iw,jw,kw,3) = rij(iw,jw,kw,3) + (r13 - rij(iw,jw,kw,3))*deltatt
    rij(iw,jw,kw,4) = rij(iw,jw,kw,4) + (r22 - rij(iw,jw,kw,4))*deltatt
    rij(iw,jw,kw,5) = rij(iw,jw,kw,5) + (r23 - rij(iw,jw,kw,5))*deltatt      
    
    
    call fvar2qmom(       q=q_mom(iw, jw, kw,:),                &
                      stress=sigma(iw, jw, kw,:),               &
                    heatflux=qflux(iw, jw, kw,:),               &
                        mijk=mijk(iw, jw, kw,:),                &
                         rij=rij(iw, jw, kw,:),                 &
                     delta=delta(iw, jw, kw)                    )

   
    return

  end  subroutine R26wbc
  !+-------------------------------------------------------------------+
  !| The end of the subroutine R26wbc.                                 |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to calculate slip orientation                  |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| xx-xx-xxxx: Created by X. Gu.                                     |
  !+-------------------------------------------------------------------+
  subroutine slip_orientation(upb, vpb, wpb, n1,n2,n3,t1, t2, t3, s1, s2, s3,  &
                  t11, t22, t33, t12, t13, t23, s11, s22, s33, s12, s13, s23)
    !==============================================================================
    implicit none
    real(8), intent(in) :: upb, vpb, wpb, n1, n2, n3
    real(8), intent(out) :: t1, t2, t3, s1, s2, s3, t11, t22, t33, t12, t13, t23, &
                            s11, s22, s33, s12, s13, s23 
    real(8) slipnorm

    slipnorm = upb*upb + vpb*vpb + wpb*wpb

    if (slipnorm == 0.0d0) then
       if (abs(n1) == 1.0d0) then
          t1 = 0.0d0
          t2 = 1.0d0
          t3 = 0.0d0
       else if (abs(n2) == 1.0d0) then
          t1 = 1.0d0
          t2 = 0.0d0
          t3 = 0.0d0
       else if (abs(n3) == 1.0d0) then
          t1 = 1.0d0
          t2 = 0.0d0
          t3 = 0.0d0
       else
          t1 = 0.2d0
          t2 = 0.2d0
          t3 = 1.0d0 - sqrt(0.2d0**2 * 2.0d0)
       end if
    else
       slipnorm = 1.0d0/sqrt(slipnorm)
       t1 = upb*slipnorm
       t2 = vpb*slipnorm
       t3 = wpb*slipnorm
    end if

    s1 = n2*t3 - n3*t2
    s2 = n3*t1 - n1*t3
    s3 = n1*t2 - n2*t1

    t11 = t1*t1 ; t22 = t2*t2 ; t33 = t3*t3
    t12 = t1*t2 ; t13 = t1*t3 ; t23 = t2*t3

    s11 = s1*s1 ; s22 = s2*s2 ; s33 = s3*s3
    s12 = s1*s2 ; s13 = s1*s3 ; s23 = s2*s3
  
    return

  end subroutine slip_orientation
  !+-------------------------------------------------------------------+
  !| The end of the subroutine slip_orientation.                       |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| These subroutines are to transfer moments to wall.                |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| xx-xx-xxxx: Created by X. Gu.                                     |
  !+-------------------------------------------------------------------+
  !+-------------------------------------------------------------------+
  !| The subroutine trans_sigma_to_wall.                               |
  !+-------------------------------------------------------------------+
  subroutine trans_sigma_to_wall(n1, n2, n3, t1, t2, t3, s1, s2, s3,       &
                                n11, n22, n33, n12, n13, n23, t11, t22, t33,    &
                                t12, t13, t23, s11, s22, s33, s12, s13, s23,    &
                                u2w, v2w, w2w, uvw, uww, vww,                   &
                                sig_nn, sig_tt, sig_ss, sig_nt, sig_ns, sig_ts)
    !============================================================================    
    real(8), intent(in)  :: n1, n2, n3, t1, t2, t3, s1, s2, s3,   &
                            n11, n22, n33, n12, n13, n23,  &
                            t11, t22, t33, t12, t13, t23,  &
                            s11, s22, s33, s12, s13, s23
    real(8), intent(in)  :: u2w, v2w, w2w, uvw, uww, vww
    real(8), intent(out) :: sig_nn, sig_tt, sig_ss, sig_nt, sig_ns, sig_ts
       
    sig_nn = u2w*n11 + v2w*n22 + w2w*n33 + 2.0d0*(uvw*n12 + uww*n13 + vww*n23)

    sig_tt = u2w*t11 + v2w*t22 + w2w*t33 + 2.0d0*(uvw*t12 + uww*t13 + vww*t23)

    sig_ss = -sig_nn - sig_tt

    sig_nt = u2w*n1*t1 + v2w*n2*t2 + w2w*n3*t3 + (n2*t1 + n1*t2)*uvw  &
           + (n3*t1 + n1*t3)*uww + (n3*t2 + n2*t3)*vww

    sig_ns = u2w*n1*s1 + v2w*n2*s2 + w2w*n3*s3 + (n2*s1 + n1*s2)*uvw  &
           + (n3*s1 + n1*s3)*uww + (n3*s2 + n2*s3)*vww

    sig_ts = u2w*t1*s1 + v2w*t2*s2 + w2w*t3*s3 + (t2*s1 + t1*s2)*uvw  &
           + (t3*s1 + t1*s3)*uww  + (t3*s2 + t2*s3)*vww 
    return

  end subroutine trans_sigma_to_wall
  !+-------------------------------------------------------------------+
  !| The end of the subroutine trans_sigma_to_wall.                    |
  !+-------------------------------------------------------------------+
  !!
  !+-------------------------------------------------------------------+
  !| The subroutine trans_sigma_to_cartisian.                          |
  !+-------------------------------------------------------------------+
  subroutine trans_sigma_to_cartisian(n1, n2, n3, t1, t2, t3, s1, s2, s3,       &
                                     n11, n22, n33, n12, n13, n23, t11, t22, t33,    &
                                     t12, t13, t23, s11, s22, s33, s12, s13, s23,    &
                                     u2w, v2w, w2w, uvw, uww, vww,                   &
                                     sig_nn, sig_tt, sig_ss, sig_nt, sig_ns, sig_ts)

    real(8), intent(in)  :: n1, n2, n3, t1, t2, t3, s1, s2, s3,   &
                            n11, n22, n33, n12, n13, n23,  &
                            t11, t22, t33, t12, t13, t23,  &
                            s11, s22, s33, s12, s13, s23

    real(8), intent(in) :: sig_nn, sig_tt, sig_ss, sig_nt, sig_ns, sig_ts
    real(8), intent(out)  :: u2w, v2w, w2w, uvw, uww, vww

    
    u2w = n11*sig_nn + t11*sig_tt + s11*sig_ss + 2.0d0*(n1*s1*sig_ns + n1*t1*sig_nt + s1*t1*sig_ts)

    v2w = n22*sig_nn + t22*sig_tt + s22*sig_ss + 2.0d0*(n2*s2*sig_ns + n2*t2*sig_nt + t2*s2*sig_ts)

    uvw = n12*sig_nn + t12*sig_tt + s12*sig_ss + (s1*n2+s2*n1)*sig_ns + (n1*t2+t1*n2)*sig_nt       &
        + (t1*s2+s1*t2)*sig_ts

    uww = n13*sig_nn + t13*sig_tt + s13*sig_ss + (n1*s3+n3*s1)*sig_ns + (n1*t3+n3*t1)*sig_nt       &
        + (t1*s3+t3*s1)*sig_ts

    vww = n23*sig_nn + t23*sig_tt + s23*sig_ss  + (n2*s3+n3*s2)*sig_ns + (n2*t3+t2*n3)*sig_nt      &
        + (t3*s2+t2*s3)*sig_ts

    w2w = - u2w - v2w

    return

  end  subroutine trans_sigma_to_cartisian 
  !+-------------------------------------------------------------------+
  !| The end of the subroutine trans_sigma_to_cartisian.               |
  !+-------------------------------------------------------------------+
  !!
  !+-------------------------------------------------------------------+
  !| The subroutine trans_mijk_to_wall                                 |
  !+-------------------------------------------------------------------+
  subroutine trans_mijk_to_wall(n1, n2, n3, t1, t2, t3, s1, s2, s3,             &
                                     n11, n22, n33, n12, n13, n23, t11, t22, t33,    &
                                     t12, t13, t23, s11, s22, s33, s12, s13, s23,    &
                            m111, m112, m113, m122, m223, m123, m222, m133, m233, m333, &
                     m_ttt, m_sss, m_tts, m_tss, m_nts, m_ntt, m_nss, m_nnt, m_nns, m_nnn)
    !============================

    real(8), intent(in)  :: n1, n2, n3, t1, t2, t3, s1, s2, s3,   &
                            n11, n22, n33, n12, n13, n23,  &
                            t11, t22, t33, t12, t13, t23,  &
                            s11, s22, s33, s12, s13, s23
    real(8), intent(in) ::  m111, m112, m113, m122, m223, m123, m222, m133, m233, m333
    real(8), intent(out) :: m_ttt, m_sss, m_tts, m_tss, m_nts, m_ntt, m_nss, m_nnt, m_nns, m_nnn

    m_ttt = m111*t11*t1 + m222*t22*t2 + m333*t33*t3   &
          + 3.0d0*( (m112*t2 + m113*t3)*t11            &
          +         (m122*t1 + m223*t3)*t22            &
          +         (m133*t1 + m233*t2)*t33            &
          +           2.0d0*m123*t3*t12  ) 

    m_nnn = m111*n11*n1 + m222*n22*n2 + m333*n33*n3   &
          + 3.0d0*( (m112*n2 + m113*n3)*n11            &
          +         (m122*n1 + m223*n3)*n22            &
          +         (m133*n1 + m233*n2)*n33            &
          +         2.0d0*m123*n3*n12 )

    m_ntt = m111*n1*t11 + m222*n2*t22 + m333*n3*t33   &
          + (n2*t11 + 2.0d0*n1*t12)*m112                &
          + (n3*t11 + 2.0d0*n1*t13)*m113                &
          + (n1*t22 + 2.0d0*n2*t12)*m122                &
          + (n3*t22 + 2.0d0*n2*t23)*m223                &
          + (n1*t33 + 2.0d0*n3*t13)*m133                &
          + (n2*t33 + 2.0d0*n3*t23)*m233                &
          + 2.0d0*(n2*t13 + n3*t12 + n1*t23)*m123 
 
    m_nnt = m111*n11*t1 + m222*n22*t2 + m333*n33*t3   &
          + (n11*t2 + 2.0d0*n12*t1)*m112                &
          + (n11*t3 + 2.0d0*n13*t1)*m113                &
          + (n22*t1 + 2.0d0*n12*t2)*m122                &
          + (n33*t1 + 2.0d0*n13*t3)*m133                &
          + (n22*t3 + 2.0d0*n23*t2)*m223                &
          + (n33*t2 + 2.0d0*n23*t3)*m233                &
          + 2.0d0*(n12*t3 + n13*t2 + n23*t1)*m123

    m_nts = m111*n1*t1*s1 + m222*n2*t2*s2 + m333*n3*t3*s3 &
          + (n2*t1*s1 + n1*t2*s1 + n1*t1*s2)*m112           &
          + (n3*t1*s1 + n1*t1*s3 + n1*t3*s1)*m113           &
          + (n2*t2*s1 + n2*t1*s2 + n1*t2*s2)*m122           &
          + (n2*t3*s1 + n3*t1*s2 + n2*t1*s3                  &
          +  n1*t3*s2 + n1*t2*s3 + n3*t2*s1)*m123           &
          + (n3*t3*s1 + n3*t1*s3 + n1*t3*s3)*m133           & 
          + (n3*t2*s2 + n2*t2*s3 + n2*t3*s2)*m223           &
          + (n2*t3*s3 + n3*t2*s3 + n3*t3*s2)*m233 

    m_tts = m111*t11*s1 + m222*t22*s2 + m333*t33*s3    &
          + (t11*s2 + 2.0d0*t12*s1)*m112                 &
          + (t11*s3 + 2.0d0*t13*s1)*m113                 &
          + (t22*s1 + 2.0d0*t12*s2)*m122                 &
          + (t33*s1 + 2.0d0*t13*s3)*m133                 &
          + (t22*s3 + 2.0d0*t23*s2)*m223                 &
          + (t33*s2 + 2.0d0*t23*s3)*m233                 &
          + 2.0d0*(t12*s3 + t13*s2 + t23*s1)*m123            

    m_nns = m111*n11*s1 + m222*n22*s2 + m333*n33*s3    &
          + (n11*s2 + 2.0d0*n12*s1)*m112                 &
          + (n11*s3 + 2.0d0*n13*s1)*m113                 &
          + (n22*s1 + 2.0d0*n12*s2)*m122                 &
          + (n33*s1 + 2.0d0*n13*s3)*m133                 &
          + (n22*s3 + 2.0d0*n23*s2)*m223                 &
          + (n33*s2 + 2.0d0*n23*s3)*m233                 &
          + 2.0d0*(n12*s3 + n13*s2 + n23*s1)*m123

    m_sss = - m_nns - m_tts  
    m_tss = - m_nnt - m_ttt
    m_nss = - m_ntt - m_nnn
 
    return

  end subroutine trans_mijk_to_wall  
  !+-------------------------------------------------------------------+
  !| The end of the subroutine trans_mijk_to_wall.                     |
  !+-------------------------------------------------------------------+
  !!
  !+-------------------------------------------------------------------+
  !| The subroutine trans_mijk_to_cartisian.                           |
  !+-------------------------------------------------------------------+
  subroutine trans_mijk_to_cartisian(n1, n2, n3, t1, t2, t3, s1, s2, s3,             &
                                     n11, n22, n33, n12, n13, n23, t11, t22, t33,    &
                                     t12, t13, t23, s11, s22, s33, s12, s13, s23,    &
                      m111, m112, m113, m122, m223, m123, m222, m133, m233, m333,     &
                     m_ttt, m_sss, m_tts, m_tss, m_nts, m_ntt, m_nss, m_nnt, m_nns, m_nnn)
    
    !======================================================================

    real(8), intent(in)  :: n1, n2, n3, t1, t2, t3, s1, s2, s3,   &
                            n11, n22, n33, n12, n13, n23,  &
                            t11, t22, t33, t12, t13, t23,  &
                            s11, s22, s33, s12, s13, s23

    real(8), intent(in) :: m_ttt, m_sss, m_tts, m_tss, m_nts, m_ntt, m_nss, m_nnt, m_nns, m_nnn
    real(8), intent(out) ::  m111, m112, m113, m122, m223, m123, m222, m133, m233, m333


    m111 = t11*t1*m_ttt + s11*s1*m_sss + n11*n1*m_nnn  &
          + 3.0d0*(t11*(s1*m_tts + n1*m_ntt)            &
                 + s11*(t1*m_tss + n1*m_nss)            &
                 + n11*(t1*m_nnt + s1*m_nns)            &
                 +  2.0d0*n1*s1*t1*m_nts)    

    m222 = t22*t2*m_ttt + s22*s2*m_sss +  n22*n2*m_nnn &
          + 3.0d0*(t22*(s2*m_tts + n2*m_ntt)            &
                 + s22*(t2*m_tss + n2*m_nss)            &
                 + n22*(t2*m_nnt + s2*m_nns)            &
                 + 2.0d0*n2*s2*t2*m_nts)   

    m112 = t11*t2*m_ttt + s11*s2*m_sss + n11*n2*m_nnn  &
          + (t11*s2 + 2.0d0*t12*s1)*m_tts               &
          + (t11*n2 + 2.0d0*t12*n1)*m_ntt               &
          + (s11*t2 + 2.0d0*t1*s12)*m_tss               &
          + (2.0d0*n12*t1 + n11*t2)*m_nnt               &                 
          + (2.0d0*s12*n1 + s11*n2)*m_nss               &
          + (2.0d0*n12*s1 + n11*s2)*m_nns               &
          + 2.0d0*(s1*t1*n2 + n1*t1*s2 + n1*s1*t2)*m_nts 

    m122 = t22*t1*m_ttt + s22*s1*m_sss + n22*n1*m_nnn   &
          + (2.0d0*t12*s2 + t22*s1)*m_tts               &
          + (t22*n1 + 2.0d0*t12*n2)*m_ntt               &
          + (2.0d0*s12*t2 + s22*t1)*m_tss               &
          + (n22*t1 + 2.0d0*n12*t2)*m_nnt               &
          + (s22*n1 + 2.0d0*s12*n2)*m_nss               &
          + (n22*s1 + 2.0d0*n12*s2)*m_nns               &
          + 2.0d0*(s2*t2*n1 + n2*t2*s1 + n2*s2*t1)*m_nts                 
 
    m123 = t12*t3*m_ttt + s12*s3*m_sss + n12*n3*m_nnn  &                  
          + (t12*s3 + t13*s2 + t23*s1)*m_tts           &
          + (t23*n1 + t13*n2 + t12*n3)*m_ntt           &
          + (s13*t2 + s12*t3 + s23*t1)*m_tss           &
          + (n23*t1 + n12*t3 + n13*t2)*m_nnt           &                 
          + (s13*n2 + s12*n3 + s23*n1)*m_nss           &
          + (n23*s1 + n12*s3 + n13*s2)*m_nns           &
          + (t1*s2*n3 + t1*s3*n2 + t2*s1*n3            &
           + t2*s3*n1 + t3*s2*n1 + t3*s1*n2)*m_nts  

    m113 = t11*t3*m_ttt + s11*s3*m_sss + n11*n3*m_nnn &
          + (t11*s3 + 2.0d0*t13*s1)*m_tts  &
          + (t11*n3 + 2.0d0*t13*n1)*m_ntt  &
          + (2.0d0*s13*t1 + s11*t3)*m_tss  &
          + (2.0d0*n13*t1 + n11*t3)*m_nnt  &
          + (2.0d0*s13*n1 + s11*n3)*m_nss  &
          + (2.0d0*n13*s1 + n11*s3)*m_nns  &
          + 2.0d0*(s1*t1*n3 + n1*t1*s3 + n1*s1*t3)*m_nts                 
                       
    m223 = t22*t3*m_ttt + s22*s3*m_sss + n22*n3*m_nnn &
          + (2.0d0*t23*s2 + t22*s3)*m_tts  &
          + (t22*n3 + 2.0d0*t23*n2)*m_ntt  &
          + (s22*t3 + 2.0d0*s23*t2)*m_tss  &
          + (2.0d0*n23*t2 + n22*t3)*m_nnt  &
          + (s22*n3 + 2.0d0*s23*n2)*m_nss  &
          + (n22*s3 + 2.0d0*n23*s2)*m_nns  &
          + 2.0d0*(s2*n2*t3 + t2*n2*s3 + t2*s2*n3)*m_nts 

    m133 = - m111 - m122
    m233 = - m112 - m222
    m333 = - m113 - m223

    return

  end subroutine trans_mijk_to_cartisian
  !+-------------------------------------------------------------------+
  !| The end of the subroutine trans_mijk_to_cartisian.                |
  !+-------------------------------------------------------------------+
  !!
  !+-------------------------------------------------------------------+
  !| The subroutine trans_rij_to_wall                                  |
  !+-------------------------------------------------------------------+
  subroutine trans_rij_to_wall(n1, n2, n3, t1, t2, t3, s1, s2, s3,         &
                           n11, n22, n33, n12, n13, n23, t11, t22, t33,    &
                           t12, t13, t23, s11, s22, s33, s12, s13, s23,    &
                           r11, r22, r12, r13, r23, r33,                   &
                          r_nn, r_tt, r_ss, r_ts, r_nt, r_ns)

    real(8), intent(in)  :: n1, n2, n3, t1, t2, t3, s1, s2, s3,   &
                            n11, n22, n33, n12, n13, n23,  &
                            t11, t22, t33, t12, t13, t23,  &
                            s11, s22, s33, s12, s13, s23
       
    real(8), intent(in)  :: r11, r22, r12, r13, r23, r33
    real(8), intent(out) :: r_nn, r_tt, r_ss, r_ts, r_nt, r_ns
    
    
    r_tt = t11*r11 + t22*r22 + t33*r33        & 
          + 2.0d0*(r12*t12 + r13*t13 + r23*t23)

    r_ss = s11*r11 + s22*r22 + s33*r33        & 
          + 2.0d0*(r12*s12 + r13*s13 + r23*s23)

    r_ts = r11*t1*s1 + r22*t2*s2 + r33*t3*s3  &
          + r12*(t1*s2 + t2*s1)               &   
          + r13*(t1*s3 + t3*s1)               & 
          + r23*(t2*s3 + t3*s2)  

    r_nt = r11*n1*t1 + r22*n2*t2 + r33*n3*t3  &
          + r12*(n1*t2 + n2*t1)               &
          + r13*(n1*t3 + n3*t1)               &
          + r23*(n2*t3 + n3*t2)

    r_ns = r11*n1*s1 + r22*n2*s2 + r33*n3*s3  &
          + r12*(n1*s2 + n2*s1)               &
          + r13*(n1*s3 + n3*s1)               &
          + r23*(n2*s3 + n3*s2)
    r_nn = - r_tt - r_ss

    return

  end subroutine trans_rij_to_wall
  !+-------------------------------------------------------------------+
  !| The end of the subroutine trans_rij_to_wall.                      |
  !+-------------------------------------------------------------------+
  !!
  !+-------------------------------------------------------------------+
  !| The subroutine trans_rij_to_cartisian.                            |
  !+-------------------------------------------------------------------+
  subroutine trans_rij_to_cartisian(n1, n2, n3, t1, t2, t3, s1, s2, s3,    &
                           n11, n22, n33, n12, n13, n23, t11, t22, t33,    &
                           t12, t13, t23, s11, s22, s33, s12, s13, s23,    &
                           r11, r22, r12, r13, r23, r33,                   &
                           r_nn, r_tt, r_ss, r_ts, r_nt, r_ns)


    real(8), intent(in)  :: n1, n2, n3, t1, t2, t3, s1, s2, s3,   &
                            n11, n22, n33, n12, n13, n23,  &
                            t11, t22, t33, t12, t13, t23,  &
                            s11, s22, s33, s12, s13, s23
    

    real(8), intent(in)  :: r_nn, r_tt, r_ss, r_ts, r_nt, r_ns
    real(8), intent(out) :: r11, r22, r12, r13, r23, r33

    r11 = n11*r_nn + t11*r_tt + s11*r_ss   &
       + 2.0d0*(n1*s1*r_ns + n1*t1*r_nt    &
       + s1*t1*r_ts )

    r22 = n22*r_nn + t22*r_tt + s22*r_ss   &
       + 2.0d0*(n2*s2*r_ns + n2*t2*r_nt    &
       + t2*s2*r_ts)

    r12 = n12*r_nn + t12*r_tt + s12*r_ss   &
       + (s1*n2 + s2*n1)*r_ns              &
       + (n1*t2 + t1*n2)*r_nt              &
       + (t1*s2 + s1*t2)*r_ts

    r13 = n13*r_nn + t13*r_tt + s13*r_ss   &
       + (n1*s3 + n3*s1)*r_ns              &
       + (n1*t3 + n3*t1)*r_nt              &
       + (t1*s3 + t3*s1)*r_ts

    r23 = n23*r_nn + t23*r_tt + s23*r_ss   &
       + (n2*s3 + n3*s2)*r_ns              &
       + (n2*t3 + t2*n3)*r_nt              &
       + (t3*s2 + t2*s3)*r_ts

    r33 = - r11 - r22

    return

  end  subroutine trans_rij_to_cartisian
  !+-------------------------------------------------------------------+
  !| The end of the subroutine trans_rij_to_cartisian.                 |
  !+-------------------------------------------------------------------+
    
  !+-------------------------------------------------------------------+
  !| The subroutine trans_psi_to_wall.                                 |
  !+-------------------------------------------------------------------+
    subroutine trans_psi_to_wall(n1, n2, n3, t1, t2, t3, s1, s2, s3,     &
                         n11, n22, n33, n12, n13, n23, t11, t22, t33,    &
                         t12, t13, t23, s11, s22, s33, s12, s13, s23,    &
                         psi111, psi112, psi113, psi122, psi223, psi123, &
                         psi222, psi133, psi233, psi333,                 &
                         psi_ttt, psi_sss, psi_tts, psi_tss, psi_nts,    &
                         psi_ntt, psi_nss, psi_nnt, psi_nns, psi_nnn)


    real(8), intent(in)  :: n1, n2, n3, t1, t2, t3, s1, s2, s3,   &
                            n11, n22, n33, n12, n13, n23,  &
                            t11, t22, t33, t12, t13, t23,  &
                            s11, s22, s33, s12, s13, s23
                               
    real(8), intent(in) ::  psi111, psi222, psi112, psi122, psi123,   &
                            psi113, psi223, psi133, psi233, psi333
    
    real(8), intent(out) :: psi_ttt, psi_sss, psi_tts, psi_tss, psi_nts,   &
                            psi_ntt, psi_nss, psi_nnt, psi_nns, psi_nnn


    psi_ttt = psi111*t11*t1 + psi222*t22*t2 + psi333*t33*t3   &
          + 3.0d0*( (psi112*t2 + psi113*t3)*t11            &
          +         (psi122*t1 + psi223*t3)*t22            &
          +         (psi133*t1 + psi233*t2)*t33            &
          +           2.0d0*psi123*t3*t12  ) 

    psi_nnn = psi111*n11*n1 + psi222*n22*n2 + psi333*n33*n3   &
          + 3.0d0*( (psi112*n2 + psi113*n3)*n11            &
          +         (psi122*n1 + psi223*n3)*n22            &
          +         (psi133*n1 + psi233*n2)*n33            &
          +         2.0d0*psi123*n3*n12 )

    psi_ntt = psi111*n1*t11 + psi222*n2*t22 + psi333*n3*t33   &
          + (n2*t11 + 2.0d0*n1*t12)*psi112                &
          + (n3*t11 + 2.0d0*n1*t13)*psi113                &
          + (n1*t22 + 2.0d0*n2*t12)*psi122                &
          + (n3*t22 + 2.0d0*n2*t23)*psi223                &
          + (n1*t33 + 2.0d0*n3*t13)*psi133                &
          + (n2*t33 + 2.0d0*n3*t23)*psi233                &
          + 2.0d0*(n2*t13 + n3*t12 + n1*t23)*psi123 
 
    psi_nnt = psi111*n11*t1 + psi222*n22*t2 + psi333*n33*t3   &
          + (n11*t2 + 2.0d0*n12*t1)*psi112                &
          + (n11*t3 + 2.0d0*n13*t1)*psi113                &
          + (n22*t1 + 2.0d0*n12*t2)*psi122                &
          + (n33*t1 + 2.0d0*n13*t3)*psi133                &
          + (n22*t3 + 2.0d0*n23*t2)*psi223                &
          + (n33*t2 + 2.0d0*n23*t3)*psi233                &
          + 2.0d0*(n12*t3 + n13*t2 + n23*t1)*psi123

    psi_nts = psi111*n1*t1*s1 + psi222*n2*t2*s2 + psi333*n3*t3*s3 &
          + (n2*t1*s1 + n1*t2*s1 + n1*t1*s2)*psi112           &
          + (n3*t1*s1 + n1*t1*s3 + n1*t3*s1)*psi113           &
          + (n2*t2*s1 + n2*t1*s2 + n1*t2*s2)*psi122           &
          + (n2*t3*s1 + n3*t1*s2 + n2*t1*s3                  &
          +  n1*t3*s2 + n1*t2*s3 + n3*t2*s1)*psi123           &
          + (n3*t3*s1 + n3*t1*s3 + n1*t3*s3)*psi133           & 
          + (n3*t2*s2 + n2*t2*s3 + n2*t3*s2)*psi223           &
          + (n2*t3*s3 + n3*t2*s3 + n3*t3*s2)*psi233 

    psi_tts = psi111*t11*s1 + psi222*t22*s2 + psi333*t33*s3    &
          + (t11*s2 + 2.0d0*t12*s1)*psi112                 &
          + (t11*s3 + 2.0d0*t13*s1)*psi113                 &
          + (t22*s1 + 2.0d0*t12*s2)*psi122                 &
          + (t33*s1 + 2.0d0*t13*s3)*psi133                 &
          + (t22*s3 + 2.0d0*t23*s2)*psi223                 &
          + (t33*s2 + 2.0d0*t23*s3)*psi233                 &
          + 2.0d0*(t12*s3 + t13*s2 + t23*s1)*psi123            

    psi_nns = psi111*n11*s1 + psi222*n22*s2 + psi333*n33*s3    &
          + (n11*s2 + 2.0d0*n12*s1)*psi112                 &
          + (n11*s3 + 2.0d0*n13*s1)*psi113                 &
          + (n22*s1 + 2.0d0*n12*s2)*psi122                 &
          + (n33*s1 + 2.0d0*n13*s3)*psi133                 &
          + (n22*s3 + 2.0d0*n23*s2)*psi223                 &
          + (n33*s2 + 2.0d0*n23*s3)*psi233                 &
          + 2.0d0*(n12*s3 + n13*s2 + n23*s1)*psi123

    psi_sss = - psi_nns - psi_tts  
    psi_tss = - psi_nnt - psi_ttt
    psi_nss = - psi_ntt - psi_nnn
 
    return

  end subroutine trans_psi_to_wall
  !+-------------------------------------------------------------------+
  !| The end of the subroutine trans_psi_to_wall.                      |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| The subroutine trans_phi_to_wall.                                 |
  !+-------------------------------------------------------------------+
  subroutine trans_phi_to_wall(n1, n2, n3, t1, t2, t3, s1, s2, s3,     &
                       n11, n22, n33, n12, n13, n23, t11, t22, t33,    &
                       t12, t13, t23, s11, s22, s33, s12, s13, s23,    &
                       q1111, q1112, q1122, q1222, q2222, q1113,       &
                       q1123, q1223, q2223, q1333, q2333, q1133,       & 
                       q2233, q1233, q3333,                            &
                       q_tttt, q_ttts, q_ttss, q_tsss, q_ssss,         &
                       q_nttt, q_ntts, q_ntss, q_nsss, q_nnnt,         &
                       q_nnns, q_nntt, q_nnss, q_nnts, q_nnnn)

    
    real(8), intent(in)  :: n1, n2, n3, t1, t2, t3, s1, s2, s3,   &
                            n11, n22, n33, n12, n13, n23,  &
                            t11, t22, t33, t12, t13, t23,  &
                            s11, s22, s33, s12, s13, s23

    real(8), intent(in) :: q1111, q1112, q1122, q1222, q2222, q1113,  &
                           q1123, q1223, q2223, q1333, q2333, q1133,  & 
                           q2233, q1233, q3333

    real(8), intent(out) :: q_tttt, q_ttts, q_ttss, q_tsss, q_ssss,   &
                            q_nttt, q_ntts, q_ntss, q_nsss, q_nnnt,   &
                            q_nnns, q_nntt, q_nnss, q_nnts, q_nnnn



    q_nnnn = (      q1111*n11 +  4.0d0*q1112*n12             &
             +  4.0d0*q1113*n13 +  6.0d0*q1122*n22           &
             + 12.0d0*q1123*n23 +  6.0d0*q1133*n33)*n11      &
             + (4.0d0*q1222*n12 + 12.0d0*q1223*n13           &
             +        q2222*n22 +  4.0d0*q2223*n23)*n22      &
             + (4.0d0*q1333*n13 +  6.0d0*q2233*n22           &
             +  4.0d0*q2333*n23 + 12.0d0*q1233*n12           &
             +        q3333*n33)*n33

    q_tttt = (      q1111*t11 +  4.0d0*q1112*t12             &
             +  4.0d0*q1113*t13 +  6.0d0*q1122*t22           &
             + 12.0d0*q1123*t23 +  6.0d0*q1133*t33)*t11      &
             + (4.0d0*q1222*t12 + 12.0d0*q1223*t13           &
             +        q2222*t22 +  4.0d0*q2223*t23)*t22      &
             + (6.0d0*q2233*t22 +  4.0d0*q2333*t23           &
             + 12.0d0*q1233*t12 +  4.0d0*q1333*t13           &
             +        q3333*t33)*t33

    q_ssss = (      q1111*s11 +  4.0d0*q1112*s12             &
             +  4.0d0*q1113*s13 +  6.0d0*q1122*s22           &
             + 12.0d0*q1123*s23 +  6.0d0*q1133*s33)*s11      &
             + (4.0d0*q1222*s12 + 12.0d0*q1223*s13           &
             +        q2222*s22 +  4.0d0*q2223*s23)*s22      &
             + (6.0d0*q2233*s22 +  4.0d0*q2333*s23           &
             + 12.0d0*q1233*s12 +  4.0d0*q1333*s13           &
             +        q3333*s33)*s33

    q_nnnt = q1111*n11*n1*t1 + q2222*n22*n2*t2 + q3333*n33*n3*t3    &
             + n11*(n1*t2 + 3.0d0*n2*t1)*q1112                        &
             + n22*(n2*t1 + 3.0d0*n1*t2)*q1222                        &
             + n22*(n2*t3 + 3.0d0*n3*t2)*q2223                        &
             + n11*(n1*t3 + 3.0d0*n3*t1)*q1113                        &
             + n33*(n3*t1 + 3.0d0*n1*t3)*q1333                        &
             + n33*(n3*t2 + 3.0d0*n2*t3)*q2333                        &
             + 3.0d0*((n2*t2*n11 + n22*n1*t1)*q1122                   &
             +        (n3*t3*n11 + n33*n1*t1)*q1133                   &
             +        (n2*t2*n33 + n22*n3*t3)*q2233                   &
             +        (2.0d0*n3*t1*n12 + (n3*t2 + n2*t3)*n11)*q1123   &
             +        (2.0d0*n3*t2*n12 + (n1*t3 + n3*t1)*n22)*q1223   &
             +        (2.0d0*n3*t3*n12 + (n2*t1 + n1*t2)*n33)*q1233 )

    q_nntt = q1111*n11*t11 + q2222*n22*t22 + q3333*n33*t33                 &
             + (4.0d0*n13*t13 + n11*t33 + n33*t11)*q1133                       &
             + (4.0d0*n12*t12 + n11*t22 + n22*t11)*q1122                       &
             + (4.0d0*n23*t23 + n22*t33 + n33*t22)*q2233                       &
             + 2.0d0*( (n11*t12 + n12*t11)*q1112                               &
             +         (n33*t13 + n13*t33)*q1333                               &
             +         (n22*t12 + n12*t22)*q1222                               &
             +         (n11*t13 + n13*t11)*q1113                               &
             +         (n23*t33 + n33*t23)*q2333                               &
             +         (n23*t22 + n22*t23)*q2223                               &
             +         (2.0d0*(n13*t12 + n12*t13) + n23*t11 + n11*t23)*q1123   &
             +         (2.0d0*(n12*t23 + n23*t12) + n13*t22 + n22*t13)*q1223   &
             +         (2.0d0*(n23*t13 + n13*t23) + n33*t12 + n12*t33)*q1233 )


    q_nttt = q1111*n1*t11*t1 + q2222*n2*t22*t2 + q3333*n3*t33*t3     &
             + (3.0d0*n1*t2 + n2*t1)*t11*q1112                         &
             + (3.0d0*n1*t3 + n3*t1)*t11*q1113                         &
             + (3.0d0*n2*t1 + n1*t2)*t22*q1222                         &
             + (3.0d0*n2*t3 + n3*t2)*t22*q2223                         &
             + (3.0d0*n3*t1 + n1*t3)*t33*q1333                         &
             + (3.0d0*n3*t2 + n2*t3)*t33*q2333                         &
             + 3.0d0*( (n2*t2*t11 + n1*t1*t22)*q1122                   &
             +         (n3*t3*t11 + n1*t1*t33)*q1133                   &
             +         (n2*t2*t33 + n3*t3*t22)*q2233                   &
             +        ((n2*t3 + n3*t2)*t11 + 2.0d0*n1*t3*t12)*q1123    &
             +        ((n3*t1 + n1*t3)*t22 + 2.0d0*n2*t3*t12)*q1223    &
             +        ((n2*t1 + n1*t2)*t33 + 2.0d0*n3*t3*t12)*q1233 )

    q_nsss = q1111*n1*s11*s1 + q2222*n2*s22*s2 +q3333*n3*s33*s3      &
             + (n2*s1 + 3.0d0*n1*s2)*s11*q1112                         &
             + (n2*s3 + 3.0d0*n3*s2)*s33*q2333                         &
             + (n1*s2 + 3.0d0*n2*s1)*s22*q1222                         &
             + (n3*s1 + 3.0d0*n1*s3)*s11*q1113                         &
             + (n3*s2 + 3.0d0*n2*s3)*s22*q2223                         &
             + (n1*s3 + 3.0d0*n3*s1)*s33*q1333                         &
             + 3.0d0*( (n3*s3*s22 + n2*s2*s33)*q2233                   &
             +         (n2*s2*s11 + n1*s1*s22)*q1122                   &
             +         (n3*s3*s11 + n1*s1*s33)*q1133                   &
             +        ((n3*s2 + n2*s3)*s11 + 2.0d0*n1*s3*s12)*q1123    &
             +        ((n1*s3 + n3*s1)*s22 + 2.0d0*n2*s3*s12)*q1223    &
             +        ((n1*s2 + n2*s1)*s33 + 2.0d0*n3*s3*s12)*q1233 )


    q_nnts = q1111*n11*s1*t1 + q2222*n22*s2*t2 + q3333*n33*s3*t3        &
             + (2.0d0*n12*s1*t1 + n11*(s2*t1 + s1*t2))*q1112              &
             + (2.0d0*n12*s2*t2 + n22*(s2*t1 + s1*t2))*q1222              &
             + (2.0d0*n13*s1*t1 + n11*(s3*t1 + s1*t3))*q1113              &
             + (2.0d0*n23*s2*t2 + n22*(s3*t2 + s2*t3))*q2223              &
             + (2.0d0*n23*s3*t3 + n33*(s2*t3 + s3*t2))*q2333              &
             + (2.0d0*n13*s3*t3 + n33*(s1*t3 + s3*t1))*q1333              &
             + (2.0d0*n12*(s2*t1 + s1*t2) + n11*s2*t2 + n22*s1*t1)*q1122  &
             + (2.0d0*n13*(s3*t1 + s1*t3) + n33*s1*t1 + n11*s3*t3)*q1133  &
             + (2.0d0*n23*(s3*t2 + s2*t3) + n22*s3*t3 + n33*s2*t2)*q2233  &
             + (2.0d0*(n12*(s3*t1 + s1*t3) + n13*(s2*t1 + s1*t2)           &
                               + n23*s1*t1) + n11*(s2*t3 + s3*t2))*q1123  &
             + (2.0d0*(n12*(s3*t2 + s2*t3) + n23*(s2*t1 + s1*t2)           &
                               + n13*s2*t2) + n22*(s1*t3 + s3*t1))*q1223  &
             + (2.0d0*(n23*(s1*t3 + s3*t1) + n13*(s2*t3 + s3*t2)           &
                               + n12*s3*t3) + n33*(s1*t2 + s2*t1))*q1233


    q_ntts = q1111*n1*t11*s1 + q2222*n2*t22*s2 + q3333*n3*t33*s3      &
             + (t11*(n1*s2 + n2*s1) + 2.0d0*n1*s1*t12)*q1112              &
             + (t22*(n2*s1 + n1*s2) + 2.0d0*n2*s2*t12)*q1222              &
             + (t11*(n1*s3 + n3*s1) + 2.0d0*n1*s1*t13)*q1113              &
             + (t22*(n3*s2 + n2*s3) + 2.0d0*n2*s2*t23)*q2223              &
             + (t33*(n3*s1 + n1*s3) + 2.0d0*n3*s3*t13)*q1333              &
             + (t33*(n3*s2 + n2*s3) + 2.0d0*n3*s3*t23)*q2333              &
             + (t22*n1*s1 + t11*n2*s2 + 2.0d0*t12*(n2*s1 + n1*s2))*q1122  &
             + (t11*n3*s3 + t33*n1*s1 + 2.0d0*t13*(n1*s3 + n3*s1))*q1133  &
             + (t22*n3*s3 + t33*n2*s2 + 2.0d0*t23*(n2*s3 + n3*s2))*q2233  &
             + (2.0d0*((n2*s1 + n1*s2)*t13 + (n3*s1 + n1*s3)*t12           &
                               + n1*s1*t23) + t11*(n3*s2 + n2*s3))*q1123  &
             + (2.0d0*((n1*s2 + n2*s1)*t23 + (n3*s2 + n2*s3)*t12           &
                               + n2*s2*t13) + t22*(n1*s3 + n3*s1))*q1223  &
             + (2.0d0*((n1*s3 + n3*s1)*t23 + (n3*s2 + n2*s3)*t13           &
                               + n3*s3*t12) + t33*(n2*s1 + n1*s2))*q1233


    q_ntss = q1111*n1*t1*s11 + q2222*n2*t2*s22 +q3333*n3*t3*s33       &
             + ((n1*t2 + n2*t1)*s11 + 2.0d0*n1*t1*s12)*q1112              &
             + ((n2*t1 + n1*t2)*s22 + 2.0d0*n2*t2*s12)*q1222              &
             + ((n1*t3 + n3*t1)*s11 + 2.0d0*n1*t1*s13)*q1113              & 
             + ((n3*t2 + n2*t3)*s22 + 2.0d0*n2*t2*s23)*q2223              &
             + ((n3*t1 + n1*t3)*s33 + 2.0d0*n3*t3*s13)*q1333              &
             + ((n2*t3 + n3*t2)*s33 + 2.0d0*n3*t3*s23)*q2333              &
             + (2.0d0*(n1*t2 + n2*t1)*s12 + n1*t1*s22 + n2*t2*s11)*q1122  &
             + (2.0d0*(n3*t1 + n1*t3)*s13 + n3*t3*s11 + n1*t1*s33)*q1133  &
             + (2.0d0*(n3*t2 + n2*t3)*s23 + n3*t3*s22 + n2*t2*s33)*q2233  &
             + (2.0d0*((n2*t1 + n1*t2)*s13 + (n1*t3 + n3*t1)*s12           &
                               + n1*t1*s23) + (n3*t2 + n2*t3)*s11)*q1123  &
             + (2.0d0*((n2*t1 + n1*t2)*s23 + (n3*t2 + n2*t3)*s12           &
                               + n2*t2*s13) + (n3*t1 + n1*t3)*s22)*q1223  &
             + (2.0d0*((n1*t3 + n3*t1)*s23 + (n3*t2 + n2*t3)*s13           &
                               + n3*t3*s12) + (n1*t2 + n2*t1)*s33)*q1233

    q_tsss = q1111*s11*s1*t1 + q2222*s22*s2*t2 + q3333*s33*s3*t3    &
             + s11*(s1*t2 + 3.0d0*s2*t1)*q1112                        &
             + s22*(s2*t1 + 3.0d0*s1*t2)*q1222                        &
             + s22*(s2*t3 + 3.0d0*s3*t2)*q2223                        &
             + s11*(s1*t3 + 3.0d0*s3*t1)*q1113                        &
             + s33*(s3*t1 + 3.0d0*s1*t3)*q1333                        &
             + s33*(s3*t2 + 3.0d0*s2*t3)*q2333                        &
             + 3.0d0*((s2*t2*s11 + s22*s1*t1)*q1122                   &
             +        (s3*t3*s11 + s33*s1*t1)*q1133                   &
             +        (s2*t2*s33 + s22*s3*t3)*q2233                   &
             +        (2.0d0*s3*t1*s12 + (s3*t2 + s2*t3)*s11)*q1123   &
             +        (2.0d0*s3*t2*s12 + (s1*t3 + s3*t1)*s22)*q1223   &
             +        (2.0d0*s3*t3*s12 + (s2*t1 + s1*t2)*s33)*q1233 )

    q_nnss = - q_nntt - q_nnnn
    q_ttss = - q_tttt - q_nntt
    q_ttts = - q_tsss - q_nnts
    q_nnns = - q_ntts - q_nsss
    
    return

  end subroutine trans_phi_to_wall
  !+-------------------------------------------------------------------+
  !| The end of the subroutine trans_phi_to_wall.                      |
  !+-------------------------------------------------------------------+

end module methodmoment
!+---------------------------------------------------------------------+
!| The end of the module  methodmoment                                 |
!+---------------------------------------------------------------------+