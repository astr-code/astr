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
                       zmin,zmax,nstep,deltat,force,numq,num_species
  use parallel, only : mpirank,lio,psum,pmax,pmin,bcast,mpistop,       &
                       irk,jrk,jrkm,ptime
  use stlaio,  only: get_unit
  !
  implicit none
  !
  integer :: nsamples,nstep_sbeg
  logical :: lmeanallocated=.false.
  logical :: liosta=.false.
  real(8) :: time_sbeg
  real(8) :: enstophy,kenergy,fbcx,massflux,massflux_target,wrms,      &
             wallheatflux,dissipation
  real(8) :: vel_incom,prs_incom,rho_incom
  real(8) :: umax,rhomax,tmpmax,qdotmax
  real(8),allocatable :: max_q(:),min_q(:)
  !
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
  subroutine meanflowcal(subtime)
    !
    use commvar,   only : reynolds,nstep,time
    use commarray, only : rho,vel,prs,tmp,dvel
    use fludyna,   only : miucal
    !
    ! arguments
    real(8),intent(inout),optional :: subtime
    !
    ! local data
    real(8) :: rho_t,prs_t,tmp_t,u,v,w
    integer :: i,j,k
    real(8) :: d11,d12,d13,d21,d22,d23,d31,d32,d33,s11,s12,s13,s22,    &
               s23,s33,skk,miu,div
    !
    real(8) :: time_beg
    !
    if(present(subtime)) time_beg=ptime() 
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
      miu=miucal(tmp_t)/reynolds
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
    if(present(subtime)) subtime=subtime+ptime()-time_beg
    !
    return
    !
  end subroutine meanflowcal
  !+-------------------------------------------------------------------+
  !| The end of the subroutine meanflowcal.                            |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to calculate and output instantous status.|
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 12-02-2021  | Created by J. Fang STFC Daresbury Laboratory        |
  !+-------------------------------------------------------------------+
  subroutine statcal(subtime)
    !
    use commvar,  only : flowtype
    use commarray,only : vel,prs,rho,tmp,spc
    use thermchem,only : heatrate
    !
    ! arguments
    real(8),intent(inout),optional :: subtime
    !
    ! local data
    integer :: i,j,k
    real(8) :: time_beg,qdot
    !
    if(present(subtime)) time_beg=ptime() 
    !
    if(trim(flowtype)=='tgv') then
      enstophy=enstophycal()
      kenergy =kenergycal()
      dissipation=diss_rate_cal()
    elseif(trim(flowtype)=='bl' .or. trim(flowtype)=='swbli') then
      fbcx=fbcxbl()
      massflux=massfluxcal()
      wrms=wrmscal()
      wallheatflux=whfbl()
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
    elseif(flowtype=='1dflame'.or.flowtype=='0dreactor' &
        .or.flowtype=='h2supersonic') then
      tmpmax=maxval(tmp(0:im,0:jm,0:km))
      tmpmax=pmax(tmpmax)
      rhomax=maxval(rho(0:im,0:jm,0:km))
      rhomax=pmax(rhomax)
      umax=maxval(vel(0:im,0:jm,0:km,1))
      umax=pmax(umax)
      !
      qdotmax=-1.d30
      do i=0,im
        do j=0,jm
          do k=0,km
            qdot=heatrate(rho(i,j,k),tmp(i,j,k),spc(i,j,k,:))
            if(qdot>qdotmax)qdotmax=qdot
          enddo 
        enddo 
      enddo  
      qdotmax=pmax(qdotmax)
      !
    endif
    !
    call maxmincal
    !
    if(present(subtime)) subtime=subtime+ptime()-time_beg
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
    ! 
    if(lio) then
      !
      if(linit) then
        !
        inquire(file='flowstate.dat',exist=fex)
        hand_fs=get_unit()
        open(hand_fs,file='flowstate.dat')
        !
        if(lrestart .and. fex) then
          ! resume flowstate.dat
          ns=0
          read(hand_fs,*)
          do while(ns<nstep)
            !
            read(hand_fs,*,iostat=ferr)ns
            !
            if(ferr< 0) then
              print*,' ** nstep=',ns,nstep
              print*,' ** end of flowstate.dat is reached.'
              exit
            endif
            !
          enddo
          !
          backspace(hand_fs)
          write(*,'(A,I0)')'   ** resume flowstate.dat at step: ',ns
          !
        else
          ! a new flowstat file
          !
          if(trim(flowtype)=='tgv') then
            write(hand_fs,"(A7,1X,A13,3(1X,A20))")'nstep','time',      &
                                      'kenergy','enstophy','dissipation'
          elseif(trim(flowtype)=='channel') then
            write(hand_fs,"(A7,1X,A13,4(1X,A20))")'nstep','time',      &
                                       'massflux','fbcx','forcex','wrms'
          elseif(trim(flowtype)=='bl' .or. trim(flowtype)=='swbli') then
            write(hand_fs,"(A7,1X,A13,4(1X,A20))")'nstep','time',      &
                                 'massflux','fbcx','wallheatflux','wrms'
          elseif(flowtype=='1dflame' .or. flowtype=='0dreactor'  &
                  .or.flowtype=='h2supersonic') then
            write(hand_fs,"(2(A10,1X),5(A12,1X))") &
              'nstep','clock','time','maxT','maxU','maxRho','maxHRR'
          else
            write(hand_fs,"(A7,1X,A13,5(1X,A20))")'nstep','time',      &
                                 'q1max','q2max','q3max','q4max','q5max'
          endif
          write(*,'(A)')'  ** create new flowstate.dat'
          !
        endif
        !
        linit=.false.
      endif
      !
      if(trim(flowtype)=='tgv') then
        write(hand_fs,"(I7,1X,E13.6E2,3(1X,E20.13E2))")nstep,time,kenergy,enstophy,dissipation
      elseif(trim(flowtype)=='channel') then
        write(hand_fs,"(I7,1X,E13.6E2,4(1X,E20.13E2))")nstep,time,massflux,fbcx,force(1),wrms
      elseif(trim(flowtype)=='bl' .or. trim(flowtype)=='swbli') then
        write(hand_fs,"(I7,1X,E13.6E2,4(1X,E20.13E2))")nstep,time,massflux,fbcx,wallheatflux,wrms
      elseif(trim(flowtype)=='cylinder') then
        write(hand_fs,"(I7,1X,E13.6E2,3(1X,E20.13E2))")nstep,time,vel_incom,prs_incom,rho_incom
      elseif(flowtype=='1dflame' .or. flowtype=='0dreactor' .or. flowtype=='h2supersonic') then
        walltime=int(ptime()-time_verybegin)
        write(hand_fs,'(2(I10,1X),5(E12.5E2,1X))') nstep,walltime,time,tmpmax,umax,rhomax,qdotmax
      else
        ! general flowstate
        write(hand_fs,"(I7,1X,E13.6E2,5(1X,E20.13E2))")nstep,time,(max_q(i),i=1,5)
      endif
      !
      if(mod(nstep,feqlist)==0) then
        !
        if(mod(nstep,nprthead*feqlist)==0) then
          if(trim(flowtype)=='tgv') then
            write(*,"(2X,A7,3(1X,A13))")'nstep','time','enstophy','kenergy'
          elseif(trim(flowtype)=='channel') then
            write(*,"(2X,A7,5(1X,A13))")'nstep','time','massflux','fbcx','forcex','wrms'
          elseif(trim(flowtype)=='bl' .or. trim(flowtype)=='swbli') then
            write(*,"(2X,A7,5(1X,A13))")'nstep','time','massflux','fbcx','wallheatflux','wrms'
          elseif(trim(flowtype)=='cylinder') then
            write(*,"(2X,A7,4(1X,A13))")'nstep','time','u_inf','p_inf','ro_inf'
          elseif(flowtype=='1dflame' .or. flowtype=='0dreactor'  &
                  .or.flowtype=='h2supersonic') then
            write(*,"(2(A10,1X),5(A12,1X))") 'nstep','clock','time','maxT', &
              'maxU','maxRho','maxHRR'
          else
            write(*,"(2X,A7,6(1X,A13))")'nstep','time','q1max','q2max','q3max','q4max','q5max'
          endif
          write(*,'(2X,78A)')('-',i=1,77)
        endif
        !
        if(trim(flowtype)=='tgv') then
          l_0=xmax/(2.d0*pi)
          write(*,"(2X,I7,1X,F13.7,2(1X,E13.6E2))")nstep,time/(l_0/uinf),enstophy,kenergy
        elseif(trim(flowtype)=='channel') then
          write(*,"(2X,I7,1X,F13.7,4(1X,E13.6E2))")nstep,time,massflux,fbcx,force(1),wrms
        elseif(trim(flowtype)=='bl' .or. trim(flowtype)=='swbli') then
          write(*,"(2X,I7,1X,F13.7,4(1X,E13.6E2))")nstep,time,massflux,fbcx,wallheatflux,wrms
        elseif(trim(flowtype)=='cylinder') then
          write(*,"(2X,I7,1X,F13.7,3(1X,E13.6E2))")nstep,time,vel_incom,prs_incom,rho_incom
        elseif(flowtype=='1dflame' .or. flowtype=='0dreactor'  &
                .or.flowtype=='h2supersonic') then
          write(*,'(2(I10,1X),5(E12.5E2,1X))') nstep,walltime,time,tmpmax,umax,rhomax,qdotmax
        else
          write(*,"(2X,I7,1X,F13.7,5(1X,E13.6E2))")nstep,time,(max_q(i),i=1,5)
        endif
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
    use commvar,  only : im,jm,km,ia,ja,ka,reynolds
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
      miu=miucal(tmp(i,j,k))/reynolds
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
    use commvar,  only : reynolds
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
        miu=miucal(tmp(i,j,k))/reynolds
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
  !| This function is to return wall wallheatflux of bl flow.          |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 29-09-2021  | Created by J. Fang STFC Daresbury Laboratory        |
  !+-------------------------------------------------------------------+
  real(8) function whfbl()
    !
    use commvar,  only : reynolds,const5,prandtl
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
        miu=miucal(tmp(i,j,k))/reynolds
        hcc=(miu/prandtl)/const5
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
    use commvar,  only : reynolds
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
        miu=miucal(tmp(i,j,k))/reynolds
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
        miu=miucal(tmp(i,j,k))/reynolds
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
               0.2d0*(massflux-massfluxtarget))/ly
    endif
    !
    ! if(lio) write(*,"(I7,5(1X,E20.13E2))")                             &
    !                           nstep,massflux,qn1,force,chanfoce,friction
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