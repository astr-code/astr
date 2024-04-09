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
  !| This subroutine is to initilsie moment equations.                 |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 10-04-2024: Created by J. Fang @ STFC Daresbury Laboratory        |
  !+-------------------------------------------------------------------+
  subroutine init_moment
    !
    use commvar, only: im,jm,km,hm,num_xtrmom
    !
    print*,allocated(q_mom)
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
  end subroutine init_moment
  !+-------------------------------------------------------------------+
  !| The end of the subroutine init_moment.                            |
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
    use bc,       only : boucon,immbody
    use comsolver,only: filterq,spongefilter,filter2e

    implicit none

    !
    ! argument
    real(8),intent(in) :: timestep
    logical,intent(in),optional :: timerept
    !
    ! local data
    logical,save :: firstcall = .true.
    real(8),save :: rkcoe(3,3)
    integer :: nrk,i,j,k,m
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
    use mpi,      only: MPI_PROC_NULL,mpi_sendrecv,mpi_real8,mpi_comm_world
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
end module methodmoment
!+---------------------------------------------------------------------+
!| The end of the module  methodmoment                                 |
!+---------------------------------------------------------------------+