!+---------------------------------------------------------------------+
!| This module contains subroutines of calculating statistics.         |
!+---------------------------------------------------------------------+
!| CHANGE RECORD                                                       |
!| -------------                                                       |
!| 12-02-2021  | Created by J. Fang                                    |
!+---------------------------------------------------------------------+
module statistic
  !
  use constdef
  use commvar,  only : im,jm,km,ia,ja,ka,ndims,xmin,xmax,ymin,ymax,    &
                       zmin,zmax,nstep,deltat,force,numq,num_species,  &
                       lreport,ltimrpt
  use parallel, only : mpirank,lio,psum,pmax,pmin,bcast,mpistop,       &
                       irk,jrk,irkm,jrkm,ptime
  use stlaio,  only: get_unit
  use utility, only: timereporter
  !
  implicit none
  !
  type :: tsta
    real(8) :: urms,machrms,macht,kolmoglength,kolmogvelocity,        &
               kolmogtime,taylorlength,retaylor
    real(8) :: rhoe
  end type tsta
  !
  integer :: nsamples,nstep_sbeg
  logical :: lmeanallocated=.false.
  logical :: liosta=.false.
  real(8) :: time_sbeg
  real(8) :: enstophy,kenergy,fbcx,massflux,massflux_target,wrms,      &
             wallheatflux,dissipation,nominal_thickness,xflame,vflame, &
             poutrt
  real(8) :: maxT,overall_qdot,v_H2O,v_HO2
  real(8) :: vel_incom,prs_incom,rho_incom
  real(8) :: umax,rhomax,tmpmax,qdotmax
  real(8),allocatable :: max_q(:),min_q(:)
  !
  real(8),allocatable,dimension(:) :: ro_xzm,u1_xzm,u2_xzm,t_xzm,p_xzm,&
                                      eng_xzm,tke_xzm,ee_xzm
  real(8),allocatable,dimension(:,:) :: ro_zm,u1_zm,u2_zm,t_zm,p_zm
  real(8),allocatable,dimension(:,:,:) :: rom,u1m,u2m,u3m,pm,tm,       &
                                          u11,u22,u33,u12,u13,u23,pp,  &
                                          tt,tu1,tu2,tu3,              &
                                          u111,u222,u333,u112,u113,    &
                                          u122,u133,u223,u233,u123,    &
                                          u1rem,u2rem,u3rem,pu1,pu2,   &
                                          pu3,sgmam11,sgmam22,sgmam33, &
                                          sgmam12,sgmam13,sgmam23,     &
                                          disspa,predil,               &
                                          visdif1,visdif2,visdif3
  !
  type(tsta) :: hitsta
  !
  contains
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to allocate mean flow arraies.            |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 19-02-2021  | Created by J. Fang STFC Daresbury Laboratory        |
  !+-------------------------------------------------------------------+
  subroutine allomeanflow
    !
    ! local data
    integer :: lallo
    !
    allocate(rom(0:im,0:jm,0:km),stat=lallo)
    if(lallo.ne.0) stop ' !! error at allocating rom'
    allocate(u1m(0:im,0:jm,0:km),stat=lallo)
    if(lallo.ne.0) stop ' !! error at allocating u1m'
    allocate(u2m(0:im,0:jm,0:km),stat=lallo)
    if(lallo.ne.0) stop ' !! error at allocating u2m'
    allocate(u3m(0:im,0:jm,0:km),stat=lallo)
    if(lallo.ne.0) stop ' !! error at allocating u3m'
    allocate(pm(0:im,0:jm,0:km),stat=lallo)
    if(lallo.ne.0) stop ' !! error at allocating pm'
    allocate(tm(0:im,0:jm,0:km),stat=lallo)
    if(lallo.ne.0) stop ' !! error at allocating tm'
    !
    allocate(u11(0:im,0:jm,0:km),stat=lallo)
    if(lallo.ne.0) stop ' !! error at allocating u11'
    allocate(u22(0:im,0:jm,0:km),stat=lallo)
    if(lallo.ne.0) stop ' !! error at allocating u22'
    allocate(u33(0:im,0:jm,0:km),stat=lallo)
    if(lallo.ne.0) stop ' !! error at allocating u33'
    allocate(u12(0:im,0:jm,0:km),stat=lallo)
    if(lallo.ne.0) stop ' !! error at allocating u12'
    allocate(u13(0:im,0:jm,0:km),stat=lallo)
    if(lallo.ne.0) stop ' !! error at allocating u13'
    allocate(u23(0:im,0:jm,0:km),stat=lallo)
    if(lallo.ne.0) stop ' !! error at allocating u23'
    allocate( pp(0:im,0:jm,0:km),stat=lallo)
    if(lallo.ne.0) stop ' !! error at allocating pp'
    allocate( tt(0:im,0:jm,0:km),stat=lallo)
    if(lallo.ne.0) stop ' !! error at allocating tt'
    allocate( tu1(0:im,0:jm,0:km),stat=lallo)
    if(lallo.ne.0) stop ' !! error at allocating tu1'
    allocate( tu2(0:im,0:jm,0:km),stat=lallo)
    if(lallo.ne.0) stop ' !! error at allocating tu2'
    allocate( tu3(0:im,0:jm,0:km),stat=lallo)
    if(lallo.ne.0) stop ' !! error at allocating tu3'
    !
    allocate( u111(0:im,0:jm,0:km),stat=lallo)
    allocate( u222(0:im,0:jm,0:km),stat=lallo)
    allocate( u333(0:im,0:jm,0:km),stat=lallo)
    allocate( u112(0:im,0:jm,0:km),stat=lallo)
    allocate( u113(0:im,0:jm,0:km),stat=lallo)
    allocate( u122(0:im,0:jm,0:km),stat=lallo)
    allocate( u133(0:im,0:jm,0:km),stat=lallo)
    allocate( u223(0:im,0:jm,0:km),stat=lallo)
    allocate( u233(0:im,0:jm,0:km),stat=lallo)
    allocate( u123(0:im,0:jm,0:km),stat=lallo)
    !
    allocate( disspa(0:im,0:jm,0:km),stat=lallo)
    allocate( predil(0:im,0:jm,0:km),stat=lallo)
    allocate( pu1(0:im,0:jm,0:km),stat=lallo)
    allocate( pu2(0:im,0:jm,0:km),stat=lallo)
    allocate( pu3(0:im,0:jm,0:km),stat=lallo)
    allocate( u1rem(0:im,0:jm,0:km),stat=lallo)
    allocate( u2rem(0:im,0:jm,0:km),stat=lallo)
    allocate( u3rem(0:im,0:jm,0:km),stat=lallo)
    allocate( visdif1(0:im,0:jm,0:km),stat=lallo)
    allocate( visdif2(0:im,0:jm,0:km),stat=lallo)
    allocate( visdif3(0:im,0:jm,0:km),stat=lallo)
    allocate( sgmam11(0:im,0:jm,0:km),stat=lallo)
    allocate( sgmam22(0:im,0:jm,0:km),stat=lallo)
    allocate( sgmam33(0:im,0:jm,0:km),stat=lallo)
    allocate( sgmam12(0:im,0:jm,0:km),stat=lallo)
    allocate( sgmam13(0:im,0:jm,0:km),stat=lallo)
    allocate( sgmam23(0:im,0:jm,0:km),stat=lallo)
    !
    lmeanallocated=.true.
    !
  end subroutine allomeanflow
  !+-------------------------------------------------------------------+
  !| The end of the subroutine allomeanflow.                           |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to calculate mean flow variables.         |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 19-02-2021  | Created by J. Fang STFC Daresbury Laboratory        |
  !+-------------------------------------------------------------------+
  subroutine meanflowcal(timerept)
    !
    use commvar,   only : reynolds,nstep,time,nondimen
    use commarray, only : rho,vel,prs,tmp,dvel
    use fludyna,   only : miucal
    !
    ! arguments
    logical,intent(in),optional :: timerept
    !
    ! local data
    real(8) :: rho_t,prs_t,tmp_t,u,v,w
    integer :: i,j,k
    real(8) :: d11,d12,d13,d21,d22,d23,d31,d32,d33,s11,s12,s13,s22,    &
               s23,s33,skk,miu,div
    !
    real(8) :: time_beg
    real(8),save :: subtime=0.d0
    !
    if(present(timerept) .and. timerept) time_beg=ptime() 
    !
    if(nsamples==0) then
      !
      if(.not. lmeanallocated) call allomeanflow
      !
      rom=0.d0
      u1m=0.d0
      u2m=0.d0
      u3m=0.d0
      pm =0.d0
      tm =0.d0
      !
      u11=0.d0
      u22=0.d0
      u33=0.d0
      u12=0.d0
      u13=0.d0
      u23=0.d0
      !
      pp =0.d0
      tt =0.d0
      tu1=0.d0
      tu2=0.d0
      tu3=0.d0
      !
      u111=0.d0
      u222=0.d0 
      u333=0.d0
      u112=0.d0
      u113=0.d0
      u122=0.d0
      u133=0.d0
      u223=0.d0
      u233=0.d0
      u123=0.d0
      !
      u1rem=0.d0
      u2rem=0.d0
      u3rem=0.d0
      !
      pu1=0.d0
      pu2=0.d0
      pu3=0.d0
      !
      sgmam11=0.d0
      sgmam12=0.d0
      sgmam13=0.d0
      sgmam22=0.d0
      sgmam13=0.d0
      sgmam33=0.d0
      !
      disspa=0.d0
      predil=0.d0
      !
      visdif1=0.d0
      visdif2=0.d0
      visdif3=0.d0
      !
      if(lio) print*,' ** meanflow initilised, nsamples=',nsamples,'nstep=',nstep
      !
      nstep_sbeg=nstep
      time_sbeg=time
    endif
    !
    do k=0,km
    do j=0,jm
    do i=0,im
      !
      rho_t=rho(i,j,k)
      prs_t=prs(i,j,k)
      tmp_t=tmp(i,j,k)
      u=vel(i,j,k,1)
      v=vel(i,j,k,2)
      w=vel(i,j,k,3)
      !
      rom(i,j,k)=rom(i,j,k) + rho_t
      u1m(i,j,k)=u1m(i,j,k) + rho_t*u
      u2m(i,j,k)=u2m(i,j,k) + rho_t*v
      u3m(i,j,k)=u3m(i,j,k) + rho_t*w
      pm(i,j,k) =pm(i,j,k)  + prs(i,j,k)
      tm(i,j,k) =tm(i,j,k)  + rho_t*tmp(i,j,k)
      !
      u11(i,j,k)=u11(i,j,k) + rho_t*u*u
      u22(i,j,k)=u22(i,j,k) + rho_t*v*v
      u33(i,j,k)=u33(i,j,k) + rho_t*w*w
      u12(i,j,k)=u12(i,j,k) + rho_t*u*v
      u13(i,j,k)=u13(i,j,k) + rho_t*u*w
      u23(i,j,k)=u23(i,j,k) + rho_t*v*w
      pp(i,j,k) =pp(i,j,k)  + prs_t*prs_t
      !
      tt(i,j,k) =tt(i,j,k)  + rho_t*tmp_t*tmp_t
      tu1(i,j,k)=tu1(i,j,k) + rho_t*u*tmp_t
      tu2(i,j,k)=tu2(i,j,k) + rho_t*v*tmp_t
      tu3(i,j,k)=tu3(i,j,k) + rho_t*w*tmp_t
      !
      u111(i,j,k)=u111(i,j,k)+rho_t*u**3
      u222(i,j,k)=u222(i,j,k)+rho_t*v**3
      u333(i,j,k)=u333(i,j,k)+rho_t*w**3
      u112(i,j,k)=u112(i,j,k)+rho_t*u*u*v
      u113(i,j,k)=u113(i,j,k)+rho_t*u*u*w
      u122(i,j,k)=u122(i,j,k)+rho_t*u*v*v
      u133(i,j,k)=u133(i,j,k)+rho_t*u*w*w
      u223(i,j,k)=u223(i,j,k)+rho_t*v*v*w
      u233(i,j,k)=u233(i,j,k)+rho_t*v*w*w
      u123(i,j,k)=u123(i,j,k)+rho_t*u*v*w
      !
      u1rem(i,j,k)=u1rem(i,j,k)+u
      u2rem(i,j,k)=u2rem(i,j,k)+v
      u3rem(i,j,k)=u3rem(i,j,k)+w
      !
      pu1(i,j,k)=pu1(i,j,k)+prs_t*u
      pu2(i,j,k)=pu2(i,j,k)+prs_t*v
      pu3(i,j,k)=pu3(i,j,k)+prs_t*w
      !
      d11=dvel(i,j,k,1,1); d12=dvel(i,j,k,1,2); d13=dvel(i,j,k,1,3)
      d21=dvel(i,j,k,2,1); d22=dvel(i,j,k,2,2); d23=dvel(i,j,k,2,3)
      d31=dvel(i,j,k,3,1); d32=dvel(i,j,k,3,2); d33=dvel(i,j,k,3,3)
      !
      if(nondimen) then
        miu=miucal(tmp(i,j,k))/reynolds
      else
        miu=miucal(tmp(i,j,k))
      endif
      !
      skk=num1d3*(d11+d22+d33)
      !
      s11=2.d0*miu*(d11-skk)
      s22=2.d0*miu*(d22-skk)
      s33=2.d0*miu*(d33-skk)
      !
      s12=miu*(d12+d21)
      s13=miu*(d13+d31)
      s23=miu*(d23+d32)
      !
      sgmam11(i,j,k)=sgmam11(i,j,k)+s11
      sgmam22(i,j,k)=sgmam22(i,j,k)+s22
      sgmam33(i,j,k)=sgmam33(i,j,k)+s33
      sgmam12(i,j,k)=sgmam12(i,j,k)+s12
      sgmam13(i,j,k)=sgmam13(i,j,k)+s13
      sgmam23(i,j,k)=sgmam23(i,j,k)+s23
      !
      visdif1(i,j,k)=visdif1(i,j,k)+s11*u+s12*v+s13*w
      visdif2(i,j,k)=visdif2(i,j,k)+s12*u+s22*v+s23*w
      visdif3(i,j,k)=visdif3(i,j,k)+s13*u+s23*v+s33*w
      !
      disspa(i,j,k)=disspa(i,j,k)+s11*d11+s12*d12+s13*d13+           &
                                  s12*d21+s22*d22+s23*d23+           &
                                  s13*d31+s23*d32+s33*d33
      !
      predil(i,j,k)=predil(i,j,k)+prs_t*(d11+d22+d33)
      !
    enddo
    enddo
    enddo
    !
    nsamples=nsamples+1
    !
    liosta=.true.
    !
    if(lio) print*,' ** meanflow calculated, nsamples=',nsamples,'nstep=',nstep
    !
    if(present(timerept) .and. timerept) then
      !
      subtime=subtime+ptime()-time_beg
      !
      if(lio .and. lreport .and. ltimrpt) call timereporter(routine='meanflowcal', &
                                              timecost=subtime, &
                                              message='data collection for statistics')
      !
    endif
    !
    return
    !
  end subroutine meanflowcal
  !+-------------------------------------------------------------------+
  !| The end of the subroutine meanflowcal.                            |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to calculate spatical averaged variables. |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 15-05-2022  | Created by J. Fang STFC Daresbury Laboratory        |
  !+-------------------------------------------------------------------+
  subroutine meanflowxzm
    !
    use commvar,  only : ndims
    use commarray,only : vel,prs,rho,tmp,x,q
    use parallel, only : mpi_ikgroup
    !
    ! local data
    integer :: i,j,k
    integer,save :: ks,ke
    logical,save :: linit=.true.
    real(8),save :: varsamp
    real(8) :: eng,ro,u1,u2,u3
    !
    if(linit) then
      allocate(ro_xzm(0:jm),u1_xzm(0:jm),u2_xzm(0:jm),t_xzm(0:jm),     &
               p_xzm(0:jm),eng_xzm(0:jm),tke_xzm(0:jm),ee_xzm(0:jm))
      !
      if(ndims==2) then
        ks=0
        ke=0
        varsamp=dble(ia)
      elseif(ndims==3) then
        ks=1
        ke=km
        varsamp=dble(ia*ka)
      endif
      !
      linit=.false.
      !
    endif
    !
    ro_xzm=0.d0
    u1_xzm=0.d0
    u2_xzm=0.d0
     t_xzm=0.d0
     p_xzm=0.d0
   eng_xzm=0.d0
    !
    tke_xzm=0.d0
    ee_xzm=0.d0
    !
    do k=ks,ke
    do i=1,im
      !
      do j=0,jm
        !
        ro=rho(i,j,k)
        eng=q(i,j,k,5)/ro
        !
        ro_xzm(j)=ro_xzm(j)+ro
         p_xzm(j)= p_xzm(j)+prs(i,j,k)
        !
        u1_xzm(j)=u1_xzm(j)+ro*vel(i,j,k,1)
        u2_xzm(j)=u2_xzm(j)+ro*vel(i,j,k,2)
         t_xzm(j)= t_xzm(j)+ro*tmp(i,j,k)
       eng_xzm(j)=eng_xzm(j)+q(i,j,k,5)
        !
        tke_xzm(j)=tke_xzm(j)+ro*(vel(i,j,k,1)**2+vel(i,j,k,2)**2+     &
                                  vel(i,j,k,3)**2)
        ee_xzm(j)=ee_xzm(j)+ro*eng**2
        !
      enddo
      !
    enddo
    enddo
    !
    p_xzm =psum( p_xzm,comm=mpi_ikgroup(jrk))
    ro_xzm=psum(ro_xzm,comm=mpi_ikgroup(jrk))
    u1_xzm=psum(u1_xzm,comm=mpi_ikgroup(jrk))
    u2_xzm=psum(u2_xzm,comm=mpi_ikgroup(jrk))
     t_xzm=psum( t_xzm,comm=mpi_ikgroup(jrk))
    eng_xzm=psum(eng_xzm,comm=mpi_ikgroup(jrk))
    !
    tke_xzm=psum(tke_xzm,comm=mpi_ikgroup(jrk))
    ee_xzm=psum(ee_xzm,comm=mpi_ikgroup(jrk))
    !
    u1_xzm=u1_xzm/ro_xzm
    u2_xzm=u2_xzm/ro_xzm
    t_xzm = t_xzm/ro_xzm
    eng_xzm=eng_xzm/ro_xzm
    !
    tke_xzm=tke_xzm/ro_xzm
    ee_xzm=ee_xzm/ro_xzm
    !
    ro_xzm=ro_xzm/varsamp
     p_xzm= p_xzm/varsamp
    !
    tke_xzm=sqrt(tke_xzm-(u1_xzm**2+u2_xzm**2))
    ee_xzm=sqrt(ee_xzm-eng_xzm**2)
    !
  end subroutine meanflowxzm
  !+-------------------------------------------------------------------+
  !| The end of the subroutine meanflowxzm.                            |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to calculate turbulence statistics.            |
  !+-------------------------------------------------------------------+
  !| ref: Samtaney, et al., Physics of Fluids, 2001,13,1415-1430.      |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 17-08-2023  | Created by J. Fang STFC Daresbury Laboratory        |
  !+-------------------------------------------------------------------+
  subroutine turbstats
    !
    use commvar,  only : reynolds,lrestart,mach,nondimen
    use commarray,only : vel,rho,tmp,dvel,q
    use fludyna,  only : miucal,sos
    use utility,  only : listinit,listwrite
    !
    type(tsta) :: stas
    !
    integer :: i,j,k,ns
    real(8) :: s11,s12,s13,s22,s23,s33,div,miu,dissa
    real(8) :: rsamples,miudrho,dudx2,csavg,v2,cs,ufluc
    !
    logical :: fex
    logical,save :: linit=.true.
    integer,save :: hand_fs
    !
    if(linit) then
      !
      if(lio) then
        call listinit(filename='turbstats.dat',handle=hand_fs, &
                          firstline='nstep time urms taylorlength kolmoglength Retaylor machrms macht roE')
      endif
      !
      linit=.false.
      !
    endif
    !
    rsamples=dble(ia*ja*ka)
    !
    stas%urms=0.d0
    stas%machrms=0.d0
    !
    miudrho=0.d0
    dudx2=0.d0
    csavg=0.d0
    dissa=0.d0
    do k=1,km
    do j=1,jm
    do i=1,im
      !
      if(nondimen) then
        miu=miucal(tmp(i,j,k))/Reynolds
      else
        miu=miucal(tmp(i,j,k))
      endif
      !
      s11=dvel(i,j,k,1,1)
      s12=0.5d0*(dvel(i,j,k,1,2)+dvel(i,j,k,2,1))
      s13=0.5d0*(dvel(i,j,k,1,3)+dvel(i,j,k,3,1))
      s22=dvel(i,j,k,2,2)
      s23=0.5d0*(dvel(i,j,k,2,3)+dvel(i,j,k,3,2))
      s33=dvel(i,j,k,3,3)
      !
      ! div=s11+s22+s33
      !
      v2=vel(i,j,k,1)**2+vel(i,j,k,2)**2+vel(i,j,k,3)**2
      !
      cs=sos(tmp(i,j,k))
      !
      stas%urms=stas%urms+v2
      !
      dudx2=dudx2+s11**2+s22**2+s33**2
      !
      miudrho=miudrho+miu/rho(i,j,k)
      !
      csavg=csavg+cs
      !
      stas%machrms=stas%machrms+v2/(cs*cs)
      !
      stas%rhoe=stas%rhoe+q(i,j,k,5)
      !
      dissa=dissa+2.d0*miu*(s11**2+s22**2+s33**2+2.d0*(s12**2+s13**2+s23**2)-num1d3*div**2)
      !
    enddo
    enddo
    enddo
    stas%urms  = sqrt(psum(stas%urms)/rsamples)
    dudx2      = num1d3*psum(dudx2)/rsamples
    miudrho    = psum(miudrho)/rsamples
    csavg      = psum(csavg)/rsamples
    dissa      = psum(dissa)/rsamples
    !
    stas%machrms=sqrt(psum(stas%machrms)/rsamples)
    !
    stas%rhoe=psum(stas%rhoe)/rsamples
    !
    ufluc=stas%urms/sqrt(3.d0)
    !
    stas%macht         = stas%urms/csavg
    stas%taylorlength  = ufluc/sqrt(dudx2)
    stas%retaylor      = ufluc*stas%taylorlength/miudrho
    stas%kolmoglength  = sqrt(sqrt(miudrho**3/dissa))
    ! stas%kolmogvelocity= sqrt(sqrt(dissipation*miudrho))
    ! stas%kolmogtime    = sqrt(miudrho/dissipation)
    !
    if(lio) call listwrite(hand_fs,stas%urms,stas%taylorlength,       &
                stas%kolmoglength,stas%Retaylor,stas%machrms,stas%macht,stas%rhoe)
    !
  end subroutine turbstats
  !+-------------------------------------------------------------------+
  !| The end of the subroutine turbstats.                              |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to calculate and output instantous status.|
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 12-02-2021  | Created by J. Fang STFC Daresbury Laboratory        |
  !+-------------------------------------------------------------------+
  subroutine statcal(timerept)
    !
    use commvar,  only : flowtype,voldom
    use commarray,only : vel,prs,rho,tmp,spc,cell,x
    use interp,   only : interlinear
#ifdef COMB
    use thermchem,only : heatrate,spcindex
#endif
    !
    ! arguments
    logical,intent(in),optional :: timerept
    !
    ! local data
    integer :: i,j,k
    real(8) :: time_beg,qdot,var1,var2
    real(8),save :: subtime=0.d0
    !
    if(present(timerept) .and. timerept) time_beg=ptime() 
    !
    if(trim(flowtype)=='tgv' .or. trim(flowtype)=='hit') then
      enstophy=enstophycal()
      kenergy =kenergycal()
      dissipation=diss_rate_cal()
      !
      if(trim(flowtype)=='hit') then
        !
        call turbstats
        !
      endif
      !
    elseif(trim(flowtype)=='bl' .or. trim(flowtype)=='swbli') then
      fbcx=fbcxbl()
      massflux=massfluxcal()
      wrms=wrmscal()
      wallheatflux=whfbl()
      !
      ! can be used for reacting case
      ! maxT=maxval(tmp(0:im,0:jm,0:km))
      ! maxT=pmax(maxT)
      ! !
      ! overall_qdot=0.d0
      ! v_H2O=0.d0
      ! v_HO2=0.d0
      ! !
      ! do i=1,im
      !   do j=1,jm
      !     do k=1,km
      !       overall_qdot=overall_qdot+heatrate(rho(i,j,k),tmp(i,j,k),spc(i,j,k,:))*cell(i,j,k)%vol
      !       v_H2O =v_H2O+spc(i,j,k,spcindex('H2O'))*rho(i,j,k)*cell(i,j,k)%vol
      !       v_HO2 =v_HO2+spc(i,j,k,spcindex('HO2'))*rho(i,j,k)*cell(i,j,k)%vol
      !     enddo 
      !   enddo 
      ! enddo   
      ! !
      ! overall_qdot=psum(overall_qdot)
      ! v_H2O =psum(v_H2O) /voldom
      ! v_HO2 =psum(v_HO2) /voldom
      !
    elseif(trim(flowtype)=='tbl') then
      !
      call meanflowxzm
      !
      fbcx=fbcxbl()
      massflux=massfluxcal()
      wrms=wrmscal()
      wallheatflux=whfbl()
      nominal_thickness=blthickness('nominal')
      !
    elseif(trim(flowtype)=='channel') then
      !
      fbcx=fbcxchan()
      !
      massflux=massfluxchan()
      !
      wrms=wrmscal()
      !
      if(nstep==0) then
        !
        massflux_target=massflux
        force=0.d0
        !
        if(lio) print*,' ** target mass flux= ',massflux_target
        !
      endif
      !
      force(1)=chanfoce(force(1),massflux,fbcx,massflux_target)
      !
      call bcast(force)
      !
    elseif(trim(flowtype)=='cylinder') then
      !
      vel_incom=0.d0
      prs_incom=0.d0
      rho_incom=0.d0
      !
      if(irk==0 .and. jrk==jrkm) then
        vel_incom=vel(0,jm,0,1)
        prs_incom=prs(0,jm,0)
        rho_incom=rho(0,jm,0)
      endif
      !
      vel_incom=psum(vel_incom)
      prs_incom=psum(prs_incom)
      rho_incom=psum(rho_incom)
      !
    elseif(flowtype=='tgvflame'.or. flowtype=='1dflame'.or.flowtype=='0dreactor' &
        .or.flowtype=='h2supersonic') then
      
#ifdef COMB
      tmpmax=maxval(tmp(0:im,0:jm,0:km))
      tmpmax=pmax(tmpmax)
      rhomax=maxval(rho(0:im,0:jm,0:km))
      rhomax=pmax(rhomax)
      umax=maxval(vel(0:im,0:jm,0:km,1))
      umax=pmax(umax)
      !
      enstophy=enstophycal()
      !
      var1=0.d0
      var2=0.d0
      !
      qdotmax=-1.d30
      !
      do i=0,im
        do j=0,jm
          do k=0,km
            !
            qdot=heatrate(rho(i,j,k),tmp(i,j,k),spc(i,j,k,:))
            if(qdot>qdotmax) then
              qdotmax=qdot
            endif
            !
          enddo 
        enddo 
      enddo  
      qdotmax=pmax(qdotmax)
      !
#endif
      !
    endif
    !
    call maxmincal
    !
    if(present(timerept) .and. timerept) then
      !
      subtime=subtime+ptime()-time_beg
      !
      if(lio .and. lreport .and. ltimrpt) call timereporter(routine='statcal', &
                                              timecost=subtime, &
                                              message='on-fly statistics')
      !
    endif
    !
  end subroutine statcal
  !+-------------------------------------------------------------------+
  !| The end of the subroutine statcal.                                |
  !+-------------------------------------------------------------------+
  !!
  !+-------------------------------------------------------------------+
  !| This subroutine is to return min and max variables.               |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 13-02-2021  | Created by J. Fang STFC Daresbury Laboratory        |
  !+-------------------------------------------------------------------+
  subroutine maxmincal
    !
    use commvar,   only : im,jm,km
    use commarray, only : q
    !
    ! local data
    integer :: n
    logical,save :: linit=.true.
    !
    if(linit) then
      allocate(max_q(1:numq),min_q(1:numq))
      linit=.false.
    endif
    !
    do n=1,numq
      !
      max_q(n)=maxval(q(0:im,0:jm,0:km,n))
      min_q(n)=minval(q(0:im,0:jm,0:km,n))
      !
      max_q(n)=pmax(max_q(n))
      min_q(n)=pmin(min_q(n))
      !
    enddo
    !
    return
    !
  end subroutine maxmincal
  !+-------------------------------------------------------------------+
  !| The end of the subroutine maxmincal.                              |
  !+-------------------------------------------------------------------+
  !!
  !+-------------------------------------------------------------------+
  !| This subroutine is used to print state on screen and write to a   |
  !| file.                                                             |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 12-02-2021  | Created by J. Fang STFC Daresbury Laboratory        |
  !+-------------------------------------------------------------------+
  subroutine statout(time_verybegin)
    !
    use commvar, only : flowtype,nstep,time,maxstep,lrestart,feqlist,uinf
    use parallel,only : ptime
    use utility,  only : listinit,listwrite
    !
    ! arguments
    real(8),intent(in),optional :: time_verybegin
    !
    ! local data
    logical,save :: linit=.true.
    logical :: fex
    integer :: nprthead=200
    integer :: i,ferr,ns,walltime
    integer,save :: hand_fs
    real(8) :: l_0
    character(len=1024) :: fstitle
    ! 
    if(lio) then
      !
      if(linit) then
        !
        if(trim(flowtype)=='tgv' .or. trim(flowtype)=='hit') then
          fstitle='nstep time kenergy enstophy dissipation'
        elseif(trim(flowtype)=='channel') then
          fstitle='nstep time massflux fbcx forcex wrms'
        elseif(trim(flowtype)=='bl' .or. trim(flowtype)=='swbli') then
          fstitle='nstep time massflux fbcx wallheatflux wrms'
          !can be used for reacting case
          ! write(hand_fs,"(A7,1X,A13,4(1X,A20))")'nstep time',      &
          !                  'maxT overall-qdot H2O_aver HO2_aver'
        elseif(trim(flowtype)=='tbl') then
          fstitle='nstep time massflux fbcx blthickness wrms'
        elseif(flowtype=='tgvflame' .or. flowtype=='1dflame' .or. flowtype=='0dreactor' .or. &
               flowtype=='h2supersonic') then
          fstitle='nstep time maxT maxU maxHRR xflame vflame pout'
        else
          fstitle='nstep time maxq1 maxq2 maxq3 maxq4 maxq5'
        endif
        !
        call listinit(filename='flowstate.dat',handle=hand_fs, &
                      firstline=trim(fstitle))
        !
        linit=.false.
        !
      endif
      !
      if(trim(flowtype)=='tgv' .or. trim(flowtype)=='hit') then
        call listwrite(hand_fs,kenergy,enstophy,dissipation)
      elseif(trim(flowtype)=='channel') then
        call listwrite(hand_fs,massflux,fbcx,force(1),wrms)
      elseif(trim(flowtype)=='bl' .or. trim(flowtype)=='tbl' .or. trim(flowtype)=='swbli') then
        call listwrite(hand_fs,massflux,fbcx,wallheatflux,wrms)
      elseif(trim(flowtype)=='cylinder') then
        call listwrite(hand_fs,vel_incom,prs_incom,rho_incom)
      elseif(flowtype=='tgvflame' .or. flowtype=='1dflame' .or. flowtype=='0dreactor' .or. &
             flowtype=='h2supersonic') then
        call listwrite(hand_fs,tmpmax,umax,qdotmax,xflame,vflame,poutrt)
      else
        ! general flowstate
        call listwrite(hand_fs,max_q(1),max_q(2),max_q(3),max_q(4),max_q(5))
      endif
      !
      if(mod(nstep,feqlist)==0) then
        !
        if(mod(nstep,nprthead*feqlist)==0) then
          write(*,"(2X,A7,6(1X,A13))")'nstep','time','q1max','q2max','q3max','q4max','q5max'
          !
          write(*,'(2X,92A)')('-',i=1,91)
        endif
        !
        write(*,"(2X,I7,1X,E13.6E2,5(1X,E13.6E2))")nstep,time,(max_q(i),i=1,5)
        !
      endif
      !
      if(nstep==maxstep) then
        close(hand_fs)
        print*,' << flowstate.dat'
      endif
      !
    endif
    !
  end subroutine statout
  !+-------------------------------------------------------------------+
  !| The end of the subroutine statprint.                              |
  !+-------------------------------------------------------------------+
  !!
  !+-------------------------------------------------------------------+
  !| This function is to return spatial averaged enstrophy.            |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 12-02-2021  | Created by J. Fang STFC Daresbury Laboratory        |
  !+-------------------------------------------------------------------+
  function enstophycal() result(vout)
    !
    use commvar,   only : im,jm,km,ia,ja,ka,roinf,uinf
    use commarray, only : vel,cell,rho,dvel
    use commfunc,  only : volhex
    !
    real(8) :: vout
    !
    ! local data
    integer :: i,j,k
    real(8) :: omega(3),omegam
    real(8) :: l_0,u_0
    !
    vout=0.d0
    !
    if(ndims==2) then
      k=0
      do j=1,jm
      do i=1,im
        !
        omega(3)=dvel(i,j,k,2,1)-dvel(i,j,k,1,2)
        omegam=omega(3)*omega(3)
        !
        vout=vout+rho(i,j,k)*omegam
        !
        !
      enddo
      enddo
      vout=0.5d0*psum(vout)/real(ia*ja,8)
    elseif(ndims==3) then
      do k=1,km
      do j=1,jm
      do i=1,im
        ! dx=
        omega(1)=dvel(i,j,k,3,2)-dvel(i,j,k,2,3)
        omega(2)=dvel(i,j,k,1,3)-dvel(i,j,k,3,1)
        omega(3)=dvel(i,j,k,2,1)-dvel(i,j,k,1,2)
        omegam=omega(1)*omega(1)+omega(2)*omega(2)+omega(3)*omega(3)
        !
        vout=vout+rho(i,j,k)*omegam
        !
        !
      enddo
      enddo
      enddo
      !
      vout=0.5d0*psum(vout)/real(ia*ja*ka,8)
    endif
    !
    l_0=xmax/(2.d0*pi)
    !
    vout=vout/(roinf*(uinf/l_0)**2)
    !
    return
    !
  end function enstophycal
  !+-------------------------------------------------------------------+
  !| The end of the subroutine enstophycal.                            |
  !+-------------------------------------------------------------------+
  !!
  !+-------------------------------------------------------------------+
  !| This function is to return spatial averaged kinetic energy.       |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 12-02-2021  | Created by J. Fang STFC Daresbury Laboratory        |
  !+-------------------------------------------------------------------+
  function kenergycal() result(vout)
    !
    use commvar,   only : im,jm,km,ia,ja,ka,roinf,uinf
    use commarray, only : vel,cell,rho,dvel
    use commfunc,  only : volhex
    !
    real(8) :: vout
    !
    ! local data
    integer :: i,j,k
    real(8) :: var1
    !
    vout=0.d0
    if(ndims==2) then
      k=0
      do j=1,jm
      do i=1,im
        ! 
        var1=vel(i,j,k,1)**2+vel(i,j,k,2)**2+vel(i,j,k,3)**2
        !
        vout=vout+rho(i,j,k)*var1
      enddo
      enddo
      !
      vout=0.5d0*psum(vout)/real(ia*ja,8)
    elseif(ndims==3) then
      do k=1,km
      do j=1,jm
      do i=1,im
        ! 
        var1=vel(i,j,k,1)**2+vel(i,j,k,2)**2+vel(i,j,k,3)**2
        !
        vout=vout+rho(i,j,k)*var1
      enddo
      enddo
      enddo
      !
      vout=0.5d0*psum(vout)/real(ia*ja*ka,8)
    endif
    !
    vout=vout/(roinf*uinf*uinf)
    !
    return
    !
  end function kenergycal
  !+-------------------------------------------------------------------+
  !| The end of the subroutine kenergycal.                             |
  !+-------------------------------------------------------------------+
  !!
  !+-------------------------------------------------------------------+
  !| This function is to return spatial averaged dissipation rate.     |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 12-02-2021  | Created by J. Fang STFC Daresbury Laboratory        |
  !+-------------------------------------------------------------------+
  function diss_rate_cal() result(vout)
    !
    use commvar,  only : im,jm,km,ia,ja,ka,reynolds,nondimen
    use commarray,only : tmp,dvel
    use fludyna,  only : miucal
    !
    real(8) :: vout
    !
    ! local data
    integer :: i,j,k
    real(8) :: var1,miu
    real(8) :: du11,du12,du13,du21,du22,du23,du31,du32,du33,           &
               s11,s12,s13,s22,s23,s33,div
    !
    vout=0.d0
    !
    do k=1,km
    do j=1,jm
    do i=1,im
      ! 
      if(nondimen) then
        miu=miucal(tmp(i,j,k))/reynolds
      else
        miu=miucal(tmp(i,j,k))
      endif
      !
      du11=dvel(i,j,k,1,1); du12=dvel(i,j,k,1,2); du13=dvel(i,j,k,1,3)
      du21=dvel(i,j,k,2,1); du22=dvel(i,j,k,2,2); du23=dvel(i,j,k,2,3)
      du31=dvel(i,j,k,3,1); du32=dvel(i,j,k,3,2); du33=dvel(i,j,k,3,3)
      !
      s11=du11; s12=0.5d0*(du12+du21); s13=0.5d0*(du13+du31)
                s22=du22;              s23=0.5d0*(du23+du32)
                                       s33=du33
      !
      div=s11+s22+s33
      !
      var1=2.d0*miu*(s11**2+s22**2+s33**2+2.d0*(s12**2+s13**2+s23**2)- &
                                                          num1d3*div**2)
      !
      vout=vout+var1
      !
    enddo
    enddo
    enddo
    !
    vout=psum(vout)/dble(ia*ja*ka)
    !
    return
    !
  end function diss_rate_cal
  !+-------------------------------------------------------------------+
  !| The end of the subroutine diss_rate_cal.                          |
  !+-------------------------------------------------------------------+
  !+-------------------------------------------------------------------+
  !| This function is to return rms spanwise velocity fluctuation.     |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 12-02-2021  | Created by J. Fang STFC Daresbury Laboratory        |
  !+-------------------------------------------------------------------+
  function wrmscal() result(vout)
    !
    use commvar,   only : im,jm,km,ia,ja,ka
    use commarray, only : vel
    !
    real(8) :: vout
    !
    ! local data
    integer :: i,j,k,k1,k2
    real(8) :: norm
    !
    if(ndims==2) then
      k1=0
      k2=0
      norm=real(ia*ja,8)
    elseif(ndims==3) then
      k1=1
      k2=km
      norm=real(ia*ja*ka,8)
    else
      print*,' !! ndims=',ndims
      stop ' !! error @ fbcxchan !!'
    endif
    !
    vout=0.d0
    do k=k1,k2
    do j=1,jm
    do i=1,im
      vout=vout+vel(i,j,k,3)*vel(i,j,k,3)
    enddo
    enddo
    enddo
    !
    vout=sqrt(psum(vout)/norm)
    !
    return
    !
  end function wrmscal
  !+-------------------------------------------------------------------+
  !| The end of the subroutine wrmscal.                                |
  !+-------------------------------------------------------------------+
  !!
  !+-------------------------------------------------------------------+
  !| This function is to return wall skin-friction of bl flow.         |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 29-09-2021  | Created by J. Fang STFC Daresbury Laboratory        |
  !+-------------------------------------------------------------------+
  real(8) function fbcxbl()
    !
    use commvar,  only : reynolds,nondimen
    use commarray,only : tmp,dvel,x
    use parallel, only : jrk,jrkm
    use fludyna,  only : miucal
    use commfunc, only : arquad
    !
    ! local data
    integer :: i,j,k
    real(8) :: miu,norm,area
    real(8) :: tau(0:im,0:km)
    !
    fbcxbl=0.d0
    norm=0.d0
    !
    if(jrk==0) then
      !
      j=0
      !
      do k=0,km
      do i=0,im
        !
        if(nondimen) then
          miu=miucal(tmp(i,j,k))/reynolds
        else
          miu=miucal(tmp(i,j,k))
        endif
        !
        tau(i,k)=miu*dvel(i,j,k,1,2)
      enddo
      enddo
      !
      if(ndims==2) then
        !
        k=0
        do i=1,im
          area=abs(x(i,j,k,1)-x(i-1,j,k,1))
          fbcxbl=fbcxbl+0.5d0*area*(tau(i,k)+tau(i-1,k))
          norm=norm+area
        enddo
        !
      elseif(ndims==3) then
        do k=1,km
        do i=1,im
          area=arquad(x(i-1,j-1,k,:),x(i,j-1,k,:),x(i,j,k,:),x(i-1,j,k,:))
          fbcxbl=fbcxbl+0.25d0*area*(tau(i,k)+tau(i-1,k)+tau(i,k-1)+tau(i-1,k-1))
          norm=norm+area
        enddo
        enddo
      else
        print*,' !! ndims=',ndims
        stop ' !! error @ fbcxbl !!'
      endif
      !
    endif
    !
    fbcxbl=psum(fbcxbl)/psum(norm)
    !
    return
    !
  end function fbcxbl
  !+-------------------------------------------------------------------+
  !| The end of the subroutine fbcxbl.                                 |
  !+-------------------------------------------------------------------+
  !!
  !+-------------------------------------------------------------------+
  !| This function is to return wall skin-friction of bl flow.         |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 29-09-2021  | Created by J. Fang STFC Daresbury Laboratory        |
  !+-------------------------------------------------------------------+
  real(8) function blthickness(which_thickness) result(thickness)
    !
    use commvar,  only : reynolds
    use commarray,only : tmp,dvel,x
    use parallel, only : jrk,jrkm
    use fludyna,  only : miucal
    use commfunc, only : arquad
    use interp,   only : interlinear
    !
    character(len=*),intent(in) :: which_thickness
    !
    ! local data
    integer :: i,j,k
    real(8) :: miu,norm,area
    real(8) :: tau(0:im,0:km)
    !
    if(which_thickness=='nominal') then
      !
      thickness=0.d0
      do j=1,jm
        !
        if(u1_xzm(j-1)<0.99d0 .and. u1_xzm(j)>=0.99d0) then
          thickness=interlinear(u1_xzm(j-1),u1_xzm(j),                 &
                                x(0,j-1,0,2),x(0,j,0,2),0.99d0)
        endif
        !
      enddo
      !
    else
      stop ' !! which_thickness not defined @ blthickness'
    endif
    !
    thickness=pmax(thickness)
    !
    ! if(thickness>1.05d0) then
    !   if(lio) print*,' ** target BL thickness achieved:',thickness
    !   call mpistop
    ! endif
    !
    return
    !
  end function blthickness
  !+-------------------------------------------------------------------+
  !| The end of the subroutine blthickness.                            |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This function is to return wall wallheatflux of bl flow.          |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 29-09-2021  | Created by J. Fang STFC Daresbury Laboratory        |
  !+-------------------------------------------------------------------+
  real(8) function whfbl()
    !
    use commvar,  only : reynolds,const5,prandtl,nondimen,cp
    use commarray,only : tmp,dtmp,x
    use parallel, only : jrk,jrkm
    use fludyna,  only : miucal
    use commfunc, only : arquad
    !
    ! local data
    integer :: i,j,k
    real(8) :: miu,hcc,norm,area
    real(8) :: qw(0:im,0:km)
    !
    whfbl=0.d0
    norm=0.d0
    !
    if(jrk==0) then
      !
      j=0
      !
      do k=0,km
      do i=0,im
        if(nondimen) then
          miu=miucal(tmp(i,j,k))/reynolds
          hcc=(miu/prandtl)/const5
        else
          miu=miucal(tmp(i,j,k))
          hcc=cp*miu/prandtl
        endif
        !
        qw(i,k)=hcc*dtmp(i,j,k,2)
      enddo
      enddo
      !
      if(ndims==2) then
        !
        k=0
        do i=1,im
          area=abs(x(i,j,k,1)-x(i-1,j,k,1))
          whfbl=whfbl+0.5d0*area*(qw(i,k)+qw(i-1,k))
          norm=norm+area
        enddo
        !
      elseif(ndims==3) then
        do k=1,km
        do i=1,im
          area=arquad(x(i-1,j-1,k,:),x(i,j-1,k,:),x(i,j,k,:),x(i-1,j,k,:))
          whfbl=whfbl+0.25d0*area*(qw(i,k)+qw(i-1,k)+qw(i,k-1)+qw(i-1,k-1))
          norm=norm+area
        enddo
        enddo
      else
        print*,' !! ndims=',ndims
        stop ' !! error @ whfbl !!'
      endif
      !
    endif
    !
    whfbl=psum(whfbl)/psum(norm)
    !
    return
    !
  end function whfbl
  !+-------------------------------------------------------------------+
  !| The end of the subroutine whfbl.                                  |
  !+-------------------------------------------------------------------+
  !!
  !+-------------------------------------------------------------------+
  !| This function is to return wall skin-friction of channel.         |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 13-02-2021  | Created by J. Fang STFC Daresbury Laboratory        |
  !+-------------------------------------------------------------------+
  real(8) function fbcxchan()
    !
    use commvar,  only : reynolds,nondimen
    use commarray,only : tmp,dvel
    use parallel, only : jrk,jrkm
    use fludyna,  only : miucal
    !
    ! local data
    integer :: i,j,k,k1,k2
    real(8) :: miu,norm
    !
    fbcxchan=0.d0
    !
    if(ndims==2) then
      k1=0
      k2=0
      norm=real(ia,8)
    elseif(ndims==3) then
      k1=1
      k2=km
      norm=real(ia*ka,8)
    else
      print*,' !! ndims=',ndims
      stop ' !! error @ fbcxchan !!'
    endif
    !
    if(jrk==0) then
      j=0
      do k=k1,k2
      do i=1,im
        !
        if(nondimen) then
          miu=miucal(tmp(i,j,k))/reynolds
        else
          miu=miucal(tmp(i,j,k))
        endif
        !
        fbcxchan=fbcxchan+miu*dvel(i,j,k,1,2)
      enddo
      enddo
      !
    endif
    !
    if(jrk==jrkm) then
      !
      j=jm
      do k=k1,k2
      do i=1,im
        !
        if(nondimen) then
          miu=miucal(tmp(i,j,k))/reynolds
        else
          miu=miucal(tmp(i,j,k))
        endif
        !
        fbcxchan=fbcxchan-miu*dvel(i,j,k,1,2)
      enddo
      enddo
    endif
    !
    fbcxchan=psum(fbcxchan)/norm
    !
    return
    !
  end function fbcxchan
  !+-------------------------------------------------------------------+
  !| The end of the subroutine fbcxchan.                               |
  !+-------------------------------------------------------------------+
  !!
  !+-------------------------------------------------------------------+
  !| This function is to return mass flux of flow.                     |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 29-09-2021  | Created by J. Fang STFC Daresbury Laboratory        |
  !+-------------------------------------------------------------------+
  real(8) function massfluxcal()
    !
    use commvar,  only : voldom
    use commarray,only : x,q,rho,cell
    !
    ! local data
    integer :: i,j,k
    real(8) :: dy,var1,norm
    !
    massfluxcal=0.d0
    !
    if(ndims==2) then
      !
      k=0
      do j=1,jm
      do i=1,im
        var1=0.25d0*(q(i-1,j-1,k,2)+q(i,j-1,k,2)+q(i-1,j,k,2)+q(i,j,k,2))
        !
        massfluxcal=massfluxcal+var1*cell(i,j,k)%vol
        !
      enddo
      enddo
      !
    elseif(ndims==3) then
      !
      do k=1,km
      do j=1,jm
      do i=1,im
        var1=0.125d0*( q(i-1,j-1,k,2)  +q(i,j-1,k,2)  +q(i-1,j,k,2)  +q(i,j,k,2) + &
                       q(i-1,j-1,k-1,2)+q(i,j-1,k-1,2)+q(i-1,j,k-1,2)+q(i,j,k-1,2) )
        !
        massfluxcal=massfluxcal+var1*cell(i,j,k)%vol
        !
      enddo
      enddo
      enddo
      !
    else
      print*,' !! ndims=',ndims
      stop ' !! error @ massfluxcal !!'
    endif
    !
    massfluxcal=psum(massfluxcal)/voldom
    !
    return
    !
  end function massfluxcal
  !+-------------------------------------------------------------------+
  !| The end of the subroutine massfluxcal.                            |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This function is to return mass flux of channel flow.             |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 13-02-2021  | Created by J. Fang STFC Daresbury Laboratory        |
  !+-------------------------------------------------------------------+
  real(8) function massfluxchan()
    !
    use commarray,only : x,q,rho
    !
    ! local data
    integer :: i,j,k,k1,k2
    real(8) :: dy,var1,norm
    !
    massfluxchan=0.d0
    !
    if(ndims==2) then
      k1=0
      k2=0
      norm=real(ia,8)
    elseif(ndims==3) then
      k1=1
      k2=km
      norm=real(ia*ka,8)
    else
      print*,' !! ndims=',ndims
      stop ' !! error @ massfluxchan !!'
    endif
    !
    do k=k1,k2
    do j=1,jm
    do i=1,im
      dy=x(i,j,k,2)-x(i,j-1,k,2)
      var1=0.5d0*(q(i,j,k,2)+q(i,j-1,k,2))
      !
      massfluxchan=massfluxchan+var1*dy
      !
    enddo
    enddo
    enddo
    !
    massfluxchan=psum(massfluxchan)/norm
    !
    return
    !
  end function massfluxchan
  !+-------------------------------------------------------------------+
  !| The end of the subroutine massfluxchan.                           |
  !+-------------------------------------------------------------------+
  !!
  !+-------------------------------------------------------------------+
  !| This subroutine is used to calcualte the body force in the channel|
  !| to main the target mass flux.                                     |
  !+-------------------------------------------------------------------+
  !| ref: Lenormand, E., Sagaut, P.  and Phuoc, L.T., Large eddy       |
  !|      simulation of subsonic and supersonic channel flow at        |
  !|      moderate Reynolds number. international journal for numerical|
  !|       methods in fluids, 2000, 32: 369-406.                       |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 22-Mar-2019  | Created by J. Fang STFC Daresbury Laboratory       |
  !+-------------------------------------------------------------------+
  real(8) function chanfoce(force,massflux,friction,massfluxtarget)
    !
    use commvar, only: time,roinf,uinf,ref_len,ref_tim,ref_miu,nondimen
    !
    ! arguments
    real(8),intent(in) :: force,massflux,friction,massfluxtarget
    !
    ! local data
    real(8) :: gn,qn1,ly
    logical,save :: linit=.true.
    !
    ly=(ymax-ymin)
    !
    if(nstep==0) then
      chanfoce=friction/ly
    else
      gn=(ly*force-friction)
      qn1=massflux+deltat*gn
      chanfoce=force-(2.d0*(qn1-massfluxtarget)-                       &
               0.2d0*(massflux-massfluxtarget))/ly*(uinf/ref_len)
    endif
    !
    ! if(nondimen) then
    !   if(lio) write(*,"(I7,6(1X,E20.13E2))")                           &
    !                      nstep,time,deltat*gn 
    !                      !massflux,qn1,force,chanfoce,friction
    ! else
    !   if(lio) write(*,"(I7,6(1X,E20.13E2))")                           &
    !           nstep,time/ref_tim,                                      &
    !           deltat*gn/(roinf*uinf*ref_len)
    !           ! massflux/(ref_den*uinf*ref_len),                      &
    !           ! qn1/(ref_den*uinf*ref_len),                           &
    !           ! force/(ref_den*uinf*uinf/ref_len),                 &
    !           ! chanfoce/(ref_den*uinf*uinf/ref_len),              &
    !           ! friction/(ref_miu*uinf/ref_len)
    ! endif
    ! !
    ! if(nstep==100) call mpistop
    !
  end function chanfoce
  !+-------------------------------------------------------------------+
  !| The end of the subroutine chanfoce.                               |
  !+-------------------------------------------------------------------+
  !!
end module statistic
!+---------------------------------------------------------------------+
!| The end of the module statistic.                                    |
!+---------------------------------------------------------------------+
