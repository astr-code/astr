!+---------------------------------------------------------------------+
!| This module contains some solver related subroutines.               |
!+---------------------------------------------------------------------+
!| CHANGE RECORD                                                       |
!| -------------                                                       |
!| 08-02-2021  | Created by J. Fang                                    |
!+---------------------------------------------------------------------+
module solver
  !
  use constdef
  use parallel, only : mpirankname,mpistop,mpirank,lio,dataswap,       &
                       datasync,ptime,irk,jrk,krk,irkm,jrkm,krkm
  use commvar,  only : ndims,ks,ke,hm,lfftk,ctime,nondimen,lreport,    &
                       ltimrpt
  use utility,  only : timereporter
  !
  implicit none
  !
  contains
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to calculate some constant parameters          |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 08-02-2021  | Created by J. Fang @ Warrington                     |
  !+-------------------------------------------------------------------+
  subroutine refcal
    !
    use commvar, only : numq,num_species,prandtl,gamma,rgas,cp,cv,     &
                        ia,ja,ka,uinf,vinf,winf,roinf,pinf,tinf,const1,&
                        const2,const3,const4,const5,const6,const7,     &
                        tempconst,tempconst1,reynolds,mach,num_modequ, &
                        turbmode,spcinf,nondimen,ref_tem,ref_vel,      &
                        ref_len,ref_den,ref_miu,ref_tim
    use thermchem, only: spcindex
    use fludyna,   only: thermal,sos,miucal
    use userdefine,only: udf_setflowenv
    use parallel,  only: mpisize
    !
    ! local data
    character(len=8) :: mpimaxname
    !
    if(trim(turbmode)=='k-omega') then
      num_modequ=2
    else
      num_modequ=0
    endif
    !
    numq=5+num_species+num_modequ
    !
    if(ia>0 .and. ja>0 .and. ka>0) then
      ndims=3
      !
      if(ia<hm .or. ja<hm .or. ka<hm) then
        print*,' !! input dimension smaller than the halo cells'
      endif
      !
    elseif(ka==0 .and. ja==0 .and. ia==0) then
      ndims=0
    elseif(ka==0 .and. ja==0) then
      ndims=1
      !
      if(ia<hm) then
        print*,' !! input dimension smaller than the halo cells'
      endif
      !
    elseif(ka==0) then
      ndims=2
      !
      if(ia<hm .or. ja<hm) then
        print*,' !! input dimension smaller than the halo cells'
      endif
      !
    else
      print*,ndims
      stop ' !! ndims error @ refcal'
    endif
    !
    prandtl=0.72d0
    !
#ifdef COMB
    rgas=287.1d0
    uinf=1.d0
    vinf=0.d0
    winf=0.d0
    pinf=1.01325d5
    tinf=300.d0

    allocate(spcinf(num_species))
    
    spcinf(:)=0.d0
    spcinf(spcindex('O2'))=0.233d0
    spcinf(spcindex('N2'))=1.d0-sum(spcinf(:))

    roinf=thermal(temperature=tinf,pressure=pinf,species=spcinf)

    ref_tim=1.d0
#else
    !
    gamma=1.4d0
    !
    if(nondimen) then 
      !
      const1=1.d0/(gamma*(gamma-1.d0)*mach**2)
      const2=gamma*mach**2
      const3=(gamma-1.d0)/3.d0*prandtl*(mach**2)
      const4=(gamma-1.d0)*mach**2*reynolds*prandtl
      const5=(gamma-1.d0)*mach**2
      const6=1.d0/(gamma-1.d0)
      const7=(gamma-1.d0)*mach**2*Reynolds*prandtl
      !
      roinf=1.d0
      uinf=1.d0
      vinf=0.d0
      winf=0.d0
      tinf=1.d0
      !
      pinf=roinf*tinf/const2
      !
      tempconst=110.3d0/ref_tem
      tempconst1=1.d0+tempconst
      !
      ref_vel=1.d0
      ref_len=1.d0
      ref_tim=ref_len/ref_vel
      !
    else 
      !
      ! rgas=287.1d0
      !
      rgas=376.177d0
      !
      cp  =gamma/(gamma-1.d0)*rgas
      cv  = rgas/(gamma-1.d0)
      !
      uinf =ref_vel
      vinf =0.d0
      winf =0.d0
      tinf =ref_tem
      roinf=ref_den
      pinf =thermal(temperature=tinf,density=roinf)
      !
      ref_miu=miucal(ref_tem)
      ref_tim=ref_len/ref_vel
      !
      mach    =ref_vel/sos(ref_tem)
      reynolds=ref_den*ref_vel*ref_len/ref_miu
      !
      if(num_species>1) then
        allocate(spcinf(num_species))
        spcinf=0.d0
      endif
      !
      const1=1.d0/(gamma*(gamma-1.d0)*mach**2)
      const2=gamma*mach**2
      const3=(gamma-1.d0)/3.d0*prandtl*(mach**2)
      const4=(gamma-1.d0)*mach**2*reynolds*prandtl
      const5=(gamma-1.d0)*mach**2
      const6=1.d0/(gamma-1.d0)
      const7=(gamma-1.d0)*mach**2*Reynolds*prandtl
      !
    endif 
    !
#endif
    !
    call udf_setflowenv
    !
    if(lio .and. ltimrpt) then
      write(mpimaxname,'(i8.8)')mpisize
      call timereporter(message=mpimaxname)
    endif
    !
  end subroutine refcal
  !+-------------------------------------------------------------------+
  !| The end of the subroutine refcal.                                 |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to calulate the rhs of N-S equations.          |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 09-02-2021  | Created by J. Fang @ Warrington                     |
  !+-------------------------------------------------------------------+
  subroutine rhscal(timerept)
    !
    use commarray, only : qrhs,x,q
    use commvar,   only : flowtype,conschm,diffterm,im,jm,             &
                          recon_schem,limmbou,lchardecomp,lihomo
    use commcal,   only : ShockSolid,ducrossensor
    use comsolver, only : gradcal
    use userdefine,only : udf_src
    use tecio
    !
    ! arguments
    logical,intent(in),optional :: timerept
    logical,save :: firstcall=.true.
    !
    ! local data
    integer :: j
    integer :: nconv
    !
    real(8) :: time_beg
    real(8),save :: subtime=0.d0
    !
    if(present(timerept) .and. timerept) time_beg=ptime() 
    !
    if(flowtype(1:2)/='0d') then
      !
      read(conschm(1:1),*) nconv
      !
      if(firstcall) then
        if(limmbou) call ShockSolid
        firstcall=.false.
      endif
      !
      call gradcal(timerept=ltimrpt)
      !
      if(mod(nconv,2)==0) then
        call convrsdcal6(timerept=ltimrpt)
      else
        !
        if(conschm(4:4)=='e') then
          !
          if(recon_schem==5 .or. lchardecomp) call ducrossensor(timerept=ltimrpt)
          !
          call convrsduwd(timerept=ltimrpt)
          !
        elseif(conschm(4:4)=='c') then
          !
          if(lchardecomp) call ducrossensor(timerept=ltimrpt)
          !
          call convrsdcmp(timerept=ltimrpt)
          !
        else
          stop ' !! error @ conschm'
        endif
        !
      endif
      !
    endif !flowtype
    !
    !
    qrhs=-qrhs
    !
    if(diffterm) call diffrsdcal6(timerept=ltimrpt)
    !
    if(trim(flowtype)=='channel') then 
      if(lihomo) call src_chan
    elseif(trim(flowtype)=='rti') then 
      call src_rti
    elseif(trim(flowtype)=='tbl') then 
      call src_tbl
    endif
    !
    call udf_src
    !
#ifdef COMB
    call srccomb(timerept=ltimrpt)
#endif
    !
    if(present(timerept) .and. timerept) then
      !
      subtime=subtime+ptime()-time_beg
      !
      if(lio .and. lreport .and. ltimrpt) call timereporter(routine='rhscal', &
                                              timecost=subtime, &
                                              message='RHS term')
    endif
    !
    ! call tecbin('testout/tecqrhs'//mpirankname//'.plt',            &
    !                                   x(0:im,0:jm,0,1),'x',        &
    !                                   x(0:im,0:jm,0,2),'y',        &
    !                                qrhs(0:im,0:jm,0,1),'qrhs1',    &
    !                                qrhs(0:im,0:jm,0,2),'qrhs2',    &
    !                                qrhs(0:im,0:jm,0,3),'qrhs3',    &
    !                                qrhs(0:im,0:jm,0,4),'qrhs4',    &
    !                                qrhs(0:im,0:jm,0,5),'qrhs5',    &
    !                                qrhs(0:im,0:jm,0,6),'qrhs6')
    ! call mpistop
    !
    return
    !
  end subroutine rhscal
  !+-------------------------------------------------------------------+
  !| The end of the subroutine rhscal.                                 |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine add a source term to the rsd of the equation to   |
  !| drive channel flow.                                               |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 13-02-2021: Created by J. Fang @ STFC Daresbury Laboratory        |
  !+-------------------------------------------------------------------+
  subroutine src_chan
    !
    use commvar,  only : force,im,jm,km
    use parallel, only : psum 
    use commarray,only : q,qrhs,x,jacob
    !
    ! local data
    integer :: i,j,k,k1,k2
    !
    real(8) :: u1bulk,u2bulk,u3bulk,robulk,dy
    !
    if(ndims==2) then
      k1=0
      k2=0
    elseif(ndims==3) then
      k1=1
      k2=km
    else
      print*,' !! ndims=',ndims
      stop ' !! error @ massfluxchan !!'
    endif
    !
    u1bulk=0.d0
    u2bulk=0.d0
    u3bulk=0.d0
    robulk=0.d0
    !
    do k=k1,k2
    do j=1,jm
    do i=1,im
      !
      dy=x(i,j,k,2)-x(i,j-1,k,2)
      robulk=robulk+0.5d0*(q(i,j-1,k,1)+q(i,j,k,1))*dy
      u1bulk=u1bulk+0.5d0*(q(i,j-1,k,2)+q(i,j,k,2))*dy
      u2bulk=u2bulk+0.5d0*(q(i,j-1,k,3)+q(i,j,k,3))*dy
      u3bulk=u3bulk+0.5d0*(q(i,j-1,k,4)+q(i,j,k,4))*dy
      !
    end do
    end do
    end do
    !
    robulk=psum(robulk)
    u1bulk=psum(u1bulk)/robulk
    u2bulk=psum(u2bulk)/robulk
    u3bulk=psum(u3bulk)/robulk
    !
    do k=0,km
    do j=0,jm
    do i=0,im
      qrhs(i,j,k,2)=qrhs(i,j,k,2)+force(1)*jacob(i,j,k)
      qrhs(i,j,k,3)=qrhs(i,j,k,3)+force(2)*jacob(i,j,k)
      qrhs(i,j,k,4)=qrhs(i,j,k,4)+force(3)*jacob(i,j,k)
      qrhs(i,j,k,5)=qrhs(i,j,k,5)+( force(1)*u1bulk+force(2)*u2bulk+   &
                                    force(3)*u3bulk )*jacob(i,j,k)
    end do
    end do
    end do
    !
  end subroutine src_chan
  !+-------------------------------------------------------------------+
  !| The end of the subroutine src_chan.                               |
  !+-------------------------------------------------------------------+
  !!
  !+-------------------------------------------------------------------+
  !| This subroutine add gravity effect to the Rayleigh–Taylor         |
  !| instability                                                       |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 09-05-2022: Created by J. Fang @ STFC Daresbury Laboratory        |
  !+-------------------------------------------------------------------+
  subroutine src_rti
    !
    use commvar,  only : im,jm,km
    use commarray,only : rho,vel,qrhs,jacob
    !
    ! local data
    integer :: i,j,k
    !
    do k=0,km
    do j=0,jm
    do i=0,im
      qrhs(i,j,k,3)=qrhs(i,j,k,3)+rho(i,j,k)*jacob(i,j,k)
      !
      qrhs(i,j,k,5)=qrhs(i,j,k,5)+rho(i,j,k)*vel(i,j,k,2)*jacob(i,j,k)
    end do
    end do
    end do
    !
  end subroutine src_rti
  !+-------------------------------------------------------------------+
  !| The end of the subroutine src_rti.                                |
  !+-------------------------------------------------------------------+
  !!
  !+-------------------------------------------------------------------+
  !| This subroutine add the slow growth terms to the temporal BL flow |
  !+-------------------------------------------------------------------+
  !| ref: V. Topalian, T. A. Oliver, R. Ulerich, and R. D. Moser,      |
  !|      Temporal slow-growth formulation for direct numerical        |
  !|      simulation of compressible wall-bounded flows,               |
  !|      Physical Review Fluids 2,(2017).                             |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 22-07-2022: Created by J. Fang @ STFC Daresbury Laboratory        |
  !+-------------------------------------------------------------------+
  subroutine src_tbl
    ! 
    use commvar,  only : im,jm,km
    use commarray,only : rho,vel,qrhs,jacob,x,q
    use comsolver,only : grad,gradcal
    use statistic,only : ro_xzm,u1_xzm,u2_xzm,eng_xzm,tke_xzm,ee_xzm,  &
                         nominal_thickness
    !
    ! local data
    integer :: i,j,k
    real(8) :: arytmp(-hm:im+hm,-hm:jm+hm,-hm:km+hm,1:6)
    real(8) :: drho(0:im,0:jm,0:km,1:3),du1(0:im,0:jm,0:km,1:3),       &
               du2(0:im,0:jm,0:km,1:3),deng(0:im,0:jm,0:km,1:3),       &
               dtke(0:im,0:jm,0:km,1:3),dee(0:im,0:jm,0:km,1:3)
    real(8) :: bl_growth_rate,s_rho,s_u1,s_u2,s_u3,s_eng,uf,vf,wf,ef
    !
    bl_growth_rate=0.001d0
    ! bl_growth_rate=0.044672d0
    !
    return
    !
    ! if(nominal_thickness<0.9d0) return
    !
    if(.not. allocated(ro_xzm)) return
    !
    do k=0,km
    do j=0,jm
    do i=0,im
      arytmp(i,j,k,1)=ro_xzm(j)
      arytmp(i,j,k,2)=u1_xzm(j)
      arytmp(i,j,k,3)=u2_xzm(j)
      arytmp(i,j,k,4)=eng_xzm(j)
      !
      arytmp(i,j,k,5)=tke_xzm(j)
      arytmp(i,j,k,6)=ee_xzm(j)
    end do
    end do
    end do
    !
    call dataswap(arytmp,timerept=ltimrpt)
    !
    drho=grad(arytmp(:,:,:,1))
    du1 =grad(arytmp(:,:,:,2))
    du2 =grad(arytmp(:,:,:,3))
    deng=grad(arytmp(:,:,:,4))
    !
    dtke=grad(arytmp(:,:,:,5))
    dee =grad(arytmp(:,:,:,6))
    !
    do k=0,km
    do j=0,jm
    do i=0,im
      !
      uf=vel(i,j,k,1)-u1_xzm(j)
      vf=vel(i,j,k,2)-u2_xzm(j)
      wf=vel(i,j,k,3)
      ef=q(i,j,k,5)/rho(i,j,k)-eng_xzm(j)
      !
      s_rho= x(i,j,k,2)*bl_growth_rate*drho(i,j,k,2)/ro_xzm(j)
      s_u1 = x(i,j,k,2)*bl_growth_rate* (du1(i,j,k,2) + uf/tke_xzm(j)*dtke(i,j,k,2))
      s_u2 = x(i,j,k,2)*bl_growth_rate* (du2(i,j,k,2) + vf/tke_xzm(j)*dtke(i,j,k,2))
      s_u3 = x(i,j,k,2)*bl_growth_rate* (               wf/tke_xzm(j)*dtke(i,j,k,2))
      s_eng= x(i,j,k,2)*bl_growth_rate* (deng(i,j,k,2)+ ef/ee_xzm(j) *dee(i,j,k,2))
      !
      ! s_rho= x(i,j,k,2)*bl_growth_rate*drho(i,j,k,2)/ro_xzm(j)
      ! s_u1 = x(i,j,k,2)*bl_growth_rate* du1(i,j,k,2)
      ! s_u2 = x(i,j,k,2)*bl_growth_rate* du2(i,j,k,2)
      ! s_u3 = 0.d0
      ! s_eng= x(i,j,k,2)*bl_growth_rate* deng(i,j,k,2)
      !
      if(isnan(s_rho)) s_rho=0.d0
      if(isnan(s_u1 )) s_u1 =0.d0
      if(isnan(s_u2 )) s_u2 =0.d0
      if(isnan(s_u3 )) s_u3 =0.d0
      if(isnan(s_eng)) s_eng=0.d0
      !
      qrhs(i,j,k,1)=qrhs(i,j,k,1)+s_rho*jacob(i,j,k)
      qrhs(i,j,k,2)=qrhs(i,j,k,2)+rho(i,j,k)*jacob(i,j,k)*(s_u1+vel(i,j,k,1)*s_rho)
      qrhs(i,j,k,3)=qrhs(i,j,k,3)+rho(i,j,k)*jacob(i,j,k)*(s_u2+vel(i,j,k,2)*s_rho)
      qrhs(i,j,k,4)=qrhs(i,j,k,4)+rho(i,j,k)*jacob(i,j,k)*(s_u3+vel(i,j,k,3)*s_rho)
      qrhs(i,j,k,5)=qrhs(i,j,k,5)+(rho(i,j,k)*s_eng+q(i,j,k,5)*s_rho)*jacob(i,j,k)
      !
    end do
    end do
    end do
    !
  end subroutine src_tbl
  !+-------------------------------------------------------------------+
  !| The end of the subroutine src_tbl.                                |
  !+-------------------------------------------------------------------+
  !+-------------------------------------------------------------------+
  !| This subroutine add a source term to the rsd of the equation to   |
  !| drive channel flow.                                               |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 13-02-2021: Created by J. Fang @ STFC Daresbury Laboratory        |
  !+-------------------------------------------------------------------+
#ifdef COMB
  subroutine srccomb(timerept)
    !
    use commvar,  only : im,jm,km,numq,num_species,odetype,lcomb
    use commarray,only : qrhs,rho,tmp,spc,jacob
    use thermchem,only : chemrate,wirate
    use parallel, only : ptime 
    !
    ! arguments
    logical,intent(in),optional :: timerept
    !
    ! local data
    integer :: i,j,k
    real(8) :: time_beg
    real(8),save :: subtime=0.d0
    !
    if(present(timerept) .and. timerept) time_beg=ptime() 
    !
    if(odetype(1:2)=='rk' .and. lcomb) then
      do k=0,km
      do j=0,jm
      do i=0,im
        call chemrate(den=rho(i,j,k),tmp=tmp(i,j,k),spc=spc(i,j,k,:))
        qrhs(i,j,k,6:numq)=qrhs(i,j,k,6:numq)+wirate(:)*jacob(i,j,k)
      enddo
      enddo 
      enddo
    endif
    !
    if(present(timerept) .and. timerept) then
      !
      subtime=subtime+ptime()-time_beg
      !
      if(lio .and. lreport .and. ltimrpt) call timereporter(routine='srccomb', &
                                             timecost=subtime, &
                                              message='SRC term for combustion')
    endif
    !
  end subroutine srccomb
#endif 
  !+-------------------------------------------------------------------+
  !| this subroutine is to solve the convectional term with upwind     |
  !| biased schemes using Steger-Warming flux splitting scheme.        |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 22-03-2021: Created by J. Fang @ Warrington.                      |
  !| add the round scheme: by Xi Deng.                                 |
  !+-------------------------------------------------------------------+
  subroutine convrsduwd(timerept)
    !
    use commvar,  only: im,jm,km,hm,numq,num_species,num_modequ,       &
                        npdci,npdcj,npdck,is,ie,js,je,ks,ke,gamma,     &
                        recon_schem,lchardecomp,conschm,bfacmpld,      &
                        nondimen
    use commarray,only: q,vel,rho,prs,tmp,spc,dxi,jacob,qrhs,lsolid,   &
                        lshock,crinod
    use commfunc, only: recons_exp
    use riemann,  only: flux_steger_warming
    use fludyna,  only: sos
#ifdef COMB
    use thermchem,only: gammarmix
#endif
    !
    ! arguments
    logical,intent(in),optional :: timerept
    !
    ! local data
    integer :: i,j,k,iss,iee,jss,jee,kss,kee,nwd
    integer :: m,n,jvar
    real(8) :: eps,gm2,var1,var2
    !
    real(8) :: lmda(5),lmdap(5),lmdam(5),gpd(3),REV(5,5),LEV(5,5),     &
               Pmult(5,5),Flcp(1:numq,1:8),Flcm(1:numq,1:8),Fhc(numq)
    !
    real(8), allocatable, dimension(:,:) :: fswp,fswm,Fh
    !
    real(8) :: time_beg
    real(8),save :: subtime=0.d0
    !
    logical :: lsh,lso,sson,hdiss,lvar,ldebug
    !
    if(present(timerept) .and. timerept) time_beg=ptime() 
    !
    sson=allocated(lshock)
    !
    eps=0.04d0
    ! gm2=0.5d0/gamma
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! calculating along i direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    allocate( fswp(-hm:im+hm,1:numq),fswm(-hm:im+hm,1:numq))
    allocate( Fh(-1:im,1:numq) )
    !
    if(npdci==1) then
      iss=0
      iee=im+hm
    elseif(npdci==2) then
      iss=-hm
      iee=im
    elseif(npdci==3) then
      iss=-hm
      iee=im+hm
    else
      stop ' !! error 1 @ subroutine convrsduwd'
    endif
    !
    do k=ks,ke
    do j=js,je
      !
      ! flux split by using Steger-Warming method
      !
      call flux_steger_warming( fplus= fswp(iss:iee,:),         &
                                fmius= fswm(iss:iee,:),         &
                                  rho=  rho(iss:iee,j,k),       &
                                  vel=  vel(iss:iee,j,k,:),     &
                                  prs=  prs(iss:iee,j,k),       &
                                  tmp=  tmp(iss:iee,j,k),       &
                                  spc=  spc(iss:iee,j,k,:),     &
                                    q=    q(iss:iee,j,k,:),     &
                                  dxi=  dxi(iss:iee,j,k,1,:),   &
                                jacob=jacob(iss:iee,j,k))
      !
      ! End of flux split by using Steger-Warming method
      !
      ! Calculating Matrix of Left and Right eigenvectors using Roe 
      ! Average and flux split as well as i+1/2 construction.
      do i=is-1,ie
        !
        ! if(i==is-1) then
        !   lso=lsolid(i+1,j,k)
        ! elseif(i==ie) then
        !   lso=lsolid(i,j,k)
        ! else
        !   lso=lsolid(i,j,k) .or. lsolid(i+1,j,k)
        ! end if
        !
        ! get value of gamma of this point

#ifdef COMB
        gamma = gammarmix(tmp(i,j,k),spc(i,j,k,:))
#endif
        !
        lso=.false.
        !
        if(sson) then
          !
          if(i<0) then
            lsh=lshock(i+1,j,k)
          elseif(i+1>im) then
            lsh=lshock(i,j,k)
          else
            lsh=lshock(i,j,k) .or. lshock(i+1,j,k)
          endif
          !
        else
          !
          lsh=.true.
          !
        endif
        !
        if(i<0) then
          hdiss=crinod(i+1,j,k)
        elseif(i+1>im) then
          hdiss=crinod(i,j,k)
        else
          hdiss=crinod(i,j,k) .or. crinod(i+1,j,k)
        endif
        !
        ! hdiss=.true.
        !
        if(lchardecomp .and. lsh) then
        ! if(lchardecomp) then
          !
          ! if(mpirank==0 .and. i==11 .and. j==0 .and. k==0) then
          !   ldebug=.true.
          ! else
          !   ldebug=.false.
          ! endif
          !
          call chardecomp(tmp(i,j,k),rho(i,j,k),prs(i,j,k),q(i,j,k,5),    &
                          vel(i,j,k,:),  spc(i,j,k,:),dxi(i,j,k,1,:),     &
                          tmp(i+1,j,k),rho(i+1,j,k), prs(i+1,j,k),q(i+1,j,k,5),  &
                          vel(i+1,j,k,:),spc(i+1,j,k,:),dxi(i+1,j,k,1,:),         &
                          REV,LEV)
          !
          ! if(mpirank==0) call check_LR55_unit(LEV,REV,i,j,k)
          ! if(irk==0) then
          !   Pmult=MatMul(LEV,REV)
          !   print*,'---------------------------------------------------------'
          !   write(*,"(5(F7.4))")Pmult(1,:)
          !   write(*,"(5(F7.4))")Pmult(2,:)
          !   write(*,"(5(F7.4))")Pmult(3,:)
          !   write(*,"(5(F7.4))")Pmult(4,:)
          !   write(*,"(5(F7.4))")Pmult(5,:)
          ! end if
          !
          ! Project to characteristic space using local eigenvector
          do m=1,5
            !
            do n=1,8
              !
              ! plus flux
              nwd=iwind8(i,n,iss,iee,'+')
              !
              Flcp(m,n)=LEV(m,1)*Fswp(nwd,1)+                        &
                        LEV(m,2)*Fswp(nwd,2)+                        &
                        LEV(m,3)*Fswp(nwd,3)+                        &
                        LEV(m,4)*Fswp(nwd,4)+                        &
                        LEV(m,5)*Fswp(nwd,5)
              !
              ! minus flux
              nwd=iwind8(i,n,iss,iee,'-')
              Flcm(m,n)=LEV(m,1)*Fswm(nwd,1)+                        &
                        LEV(m,2)*Fswm(nwd,2)+                        &
                        LEV(m,3)*Fswm(nwd,3)+                        &
                        LEV(m,4)*Fswm(nwd,4)+                        &
                        LEV(m,5)*Fswm(nwd,5)
            end do
            !
          end do
          !
          if(numq>5) then
            !
            do n=1,8
              !
              nwd=iwind8(i,n,iss,iee,'+')
              Flcp(6:numq,n)=Fswp(nwd,6:numq)
              !
              nwd=iwind8(i,n,iss,iee,'-')
              Flcm(6:numq,n)=Fswm(nwd,6:numq)
              !
            end do
            !
          endif
          ! End of characteristic decomposition.
          !
        else
          ! No characteristic decomposition
          !
          do n=1,8
            ! plus flux
            nwd=iwind8(i,n,iss,iee,'+')
            Flcp(1:numq,n)=Fswp(nwd,1:numq)
            !
            ! minus flux
            nwd=iwind8(i,n,iss,iee,'-')
            Flcm(1:numq,n)=Fswm(nwd,1:numq)
          end do
          !
        endif
        !
        ! Calculating values at i+1/2 using shock-capturing scheme.
        !
        do m=1,numq
          !
          if(hdiss) then
            !
            var1=Flcp(m,4)
            var2=Flcm(m,4)
            !
          else
            !
            var1=recons_exp(   f=Flcp(m,:), inode=i,       &
                             dim=im,        ntype=npdci,   &
                             reschem=recon_schem,shock=lsh,solid=lso)
            var2=recons_exp(   f=Flcm(m,:), inode=i,       &
                             dim=im,        ntype=npdci,   &
                             reschem=recon_schem,shock=lsh,solid=lso)
            !
          endif
          !
          Fhc(m)=var1+var2
          !
        enddo
        !
        if(lchardecomp .and. lsh) then
        ! if(lchardecomp) then
          do m=1,5
            Fh(i,m)=REV(m,1)*Fhc(1)+REV(m,2)*Fhc(2)+REV(m,3)*Fhc(3)+   &
                    REV(m,4)*Fhc(4)+REV(m,5)*Fhc(5) 
          end do
          !
          if(numq>5) then
            Fh(i,6:numq)=Fhc(6:numq)
          endif
          !
        else
          Fh(i,1:numq)=Fhc(1:numq)
        endif
        !
      enddo
      
      do i=is,ie
        do m=1,numq
          qrhs(i,j,k,m)=qrhs(i,j,k,m)+Fh(i,m)-Fh(i-1,m)
        enddo
        !
      enddo
      !
    enddo
    enddo
    deallocate( Fswp,Fswm,Fh )
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! end of calculation at i direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! calculating along j direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    if(npdcj==1) then
      jss=0
      jee=jm+hm
    elseif(npdcj==2) then
      jss=-hm
      jee=jm
    elseif(npdcj==3) then
      jss=-hm
      jee=jm+hm
    else
      stop ' !! error 2 @ subroutjne convrsduwd'
    endif
    !
    allocate( fswp(-hm:jm+hm,1:numq),fswm(-hm:jm+hm,1:numq))
    allocate( Fh(-1:jm,1:numq) )
    !
    do k=ks,ke
    do i=is,ie
      !
      ! flux split by using Steger-Warming method
      call flux_steger_warming( fplus= fswp(jss:jee,:),         &
                                fmius= fswm(jss:jee,:),         &
                                  rho=  rho(i,jss:jee,k),       &
                                  vel=  vel(i,jss:jee,k,:),     &
                                  prs=  prs(i,jss:jee,k),       &
                                  tmp=  tmp(i,jss:jee,k),       &
                                  spc=  spc(i,jss:jee,k,:),     &
                                    q=    q(i,jss:jee,k,:),     &
                                  dxi=  dxi(i,jss:jee,k,2,:),   &
                                jacob=jacob(i,jss:jee,k))
      !
      ! Calculating Matrix of Left and Right eigenvectors using Roe 
      ! Average and flux split as well as i+1/2 construction.
      do j=js-1,je
        !
        ! if(j==js-1) then
        !   lso=lsolid(i,j+1,k)
        ! elseif(j==je) then
        !   lso=lsolid(i,j,k)
        ! else
        !   lso=lsolid(i,j,k) .or. lsolid(i,j+1,k)
        ! end if
        ! 
        ! get value of gamma of this point
#ifdef COMB
          gamma = gammarmix(tmp(i,j,k),spc(i,j,k,:))
#endif
        !
        lso=.false.
        !
        if(sson) then
          !
          if(j<0) then
            lsh=lshock(i,j+1,k)
          elseif(j+1>jm) then
            lsh=lshock(i,j,k)
          else
            lsh=lshock(i,j,k) .or. lshock(i,j+1,k)
          endif
          !
        else
          lsh=.true.
        endif
        !
        if(j<0) then
          hdiss=crinod(i,j+1,k)
        elseif(j+1>jm) then
          hdiss=crinod(i,j,k)
        else
          hdiss=crinod(i,j,k) .or. crinod(i,j+1,k)
        endif
        !
        ! hdiss=.true.
        !
        ! lchardecomp=.false.
        if(lchardecomp .and. lsh) then
        ! if(lchardecomp) then
          !
          call chardecomp(tmp(i,j,k),rho(i,j,k),prs(i,j,k),q(i,j,k,5),   &
                          vel(i,j,k,:),  spc(i,j,k,:),dxi(i,j,k,2,:),    &
                          tmp(i,j+1,k),rho(i,j+1,k),  prs(i,j+1,k),q(i,j+1,k,5),  &
                          vel(i,j+1,k,:),spc(i,j+1,k,:),dxi(i,j+1,k,2,:),         &
                          REV,LEV)
          !
          ! Pmult=MatMul(LEV,REV)
          !
          ! call check_mat55_unit(Pmult,lvar)
          ! if(irk==0) then
          !   Pmult=MatMul(LEV,REV)
          !   print*,'---------------------------------------------------------'
          !   write(*,"(5(F7.4))")Pmult(1,:)
          !   write(*,"(5(F7.4))")Pmult(2,:)
          !   write(*,"(5(F7.4))")Pmult(3,:)
          !   write(*,"(5(F7.4))")Pmult(4,:)
          !   write(*,"(5(F7.4))")Pmult(5,:)
          ! end if
          !
          ! Project to characteristic space using local eigenvector
          do m=1,5
            !
            do n=1,8
              !
              ! plus flux
              nwd=iwind8(j,n,jss,jee,'+')
              Flcp(m,n)=LEV(m,1)*Fswp(nwd,1)+                        &
                        LEV(m,2)*Fswp(nwd,2)+                        &
                        LEV(m,3)*Fswp(nwd,3)+                        &
                        LEV(m,4)*Fswp(nwd,4)+                        &
                        LEV(m,5)*Fswp(nwd,5)
              !
              ! minus flux
              nwd=iwind8(j,n,jss,jee,'-')
              Flcm(m,n)=LEV(m,1)*Fswm(nwd,1)+                        &
                        LEV(m,2)*Fswm(nwd,2)+                        &
                        LEV(m,3)*Fswm(nwd,3)+                        &
                        LEV(m,4)*Fswm(nwd,4)+                        &
                        LEV(m,5)*Fswm(nwd,5)
            end do
            !
          end do
          !
          if(numq>5) then
            do n=1,8
              nwd=iwind8(j,n,jss,jee,'+')
              Flcp(6:numq,n)=Fswp(nwd,6:numq)
              !
              nwd=iwind8(j,n,jss,jee,'-')
              Flcm(6:numq,n)=Fswm(nwd,6:numq)
            end do
          endif
          !
          ! End of characteristic decomposition.
          !
        else
          ! No characteristic decomposition
          !
          do n=1,8
            ! plus flux
            nwd=iwind8(j,n,jss,jee,'+')
            Flcp(1:numq,n)=Fswp(nwd,1:numq)
            !
            ! minus flux
            nwd=iwind8(j,n,jss,jee,'-')
            Flcm(1:numq,n)=Fswm(nwd,1:numq)
          end do
          !
        endif
        !
        ! Calculating values at i+1/2 using shock-capturing scheme.
        do m=1,numq
          !
          if(hdiss) then
            !
            var1=Flcp(m,4)
            var2=Flcm(m,4)
            !
          else
            !
            var1=recons_exp(   f=Flcp(m,:), inode=j,       &
                             dim=jm,        ntype=npdcj,   &
                             reschem=recon_schem,shock=lsh,solid=lso)
            var2=recons_exp(   f=Flcm(m,:), inode=j,       &
                             dim=jm,        ntype=npdcj,   &
                             reschem=recon_schem,shock=lsh,solid=lso)
            !
          endif
          !
          Fhc(m)=var1+var2
          !
        enddo
        !
        if(lchardecomp .and. lsh) then
        ! if(lchardecomp) then
          do m=1,5
            Fh(j,m)=REV(m,1)*Fhc(1)+REV(m,2)*Fhc(2)+REV(m,3)*Fhc(3)+   &
                    REV(m,4)*Fhc(4)+REV(m,5)*Fhc(5) 
          end do
          !
          if(numq>5) then
            Fh(j,6:numq)=Fhc(6:numq)
          endif
        else
          Fh(j,1:numq)=Fhc(1:numq)
        endif
        !
      enddo
      
      do j=js,je
        do m=1,numq
          qrhs(i,j,k,m)=qrhs(i,j,k,m)+Fh(j,m)-Fh(j-1,m)
        enddo
      enddo
      !
    enddo
    enddo
    deallocate( Fswp,Fswm,Fh )
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! end of calculation at j direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    if(ndims==3) then
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! calculating along k direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    if(npdck==1) then
      kss=0
      kee=km+hm
    elseif(npdck==2) then
      kss=-hm
      kee=km
    elseif(npdck==3) then
      kss=-hm
      kee=km+hm
    else
      stop ' !! error 2 @ subroutkne convrsduwd'
    endif
    !
    allocate( fswp(-hm:km+hm,1:numq),fswm(-hm:km+hm,1:numq))
    allocate( Fh(-1:km,1:numq) )
    !
    do j=js,je
    do i=is,ie
      !
      ! flux split by using Steger-Warming method
      call flux_steger_warming( fplus= fswp(kss:kee,:),         &
                                fmius= fswm(kss:kee,:),         &
                                  rho=  rho(i,j,kss:kee),       &
                                  vel=  vel(i,j,kss:kee,:),     &
                                  prs=  prs(i,j,kss:kee),       &
                                  tmp=  tmp(i,j,kss:kee),       &
                                  spc=  spc(i,j,kss:kee,:),     &
                                    q=    q(i,j,kss:kee,:),     &
                                  dxi=  dxi(i,j,kss:kee,3,:),   &
                                jacob=jacob(i,j,kss:kee))
      ! End of flux split by using Steger-Warming method
      !
      !
      ! Calculating Matrix of Left and Right eigenvectors using Roe 
      ! Average and flux split as well as i+1/2 construction.
      do k=ks-1,ke
        !
        ! if(k==ks-1) then
        !   lso=lsolid(i,j,k+1)
        ! elseif(k==ke) then
        !   lso=lsolid(i,j,k)
        ! else
        !   lso=lsolid(i,j,k) .or. lsolid(i,j,k+1)
        ! end if
        !
        ! get value of gamma of this point
        
#ifdef COMB
        gamma = gammarmix(tmp(i,j,k),spc(i,j,k,:))
#endif
        ! 
        lso=.false.
        !
        if(sson) then
          !
          if(k<0) then
            lsh=lshock(i,j,k+1)
          elseif(k+1>km) then
            lsh=lshock(i,j,k)
          else
            lsh=lshock(i,j,k) .or. lshock(i,j,k+1)
          endif
          !
        else
          lsh=.true.
        endif
        !
        if(k<0) then
          hdiss=crinod(i,j,k+1)
        elseif(k+1>km) then
          hdiss=crinod(i,j,k)
        else
          hdiss=crinod(i,j,k) .or. crinod(i,j,k+1)
        endif
        !
        ! hdiss=.true.
        !
        ! if(lchardecomp .and. lsh) then
        ! if(lchardecomp) then
        if(.false.) then
          !
          call chardecomp(tmp(i,j,k),rho(i,j,k), prs(i,j,k), q(i,j,k,5),      &
                          vel(i,j,k,:),  spc(i,j,k,:),dxi(i,j,k,3,:),         &
                          tmp(i,j,k+1),rho(i,j,k+1),  prs(i,j,k+1),q(i,j,k+1,5), &
                          vel(i,j,k+1,:),spc(i,j,k+1,:),dxi(i,j,k+1,3,:), &
                          REV,LEV)
          !
          ! if(mpirank==0) call check_LR55_unit(LEV,REV,i,j,k)
          !
          ! Project to characteristic space using local eigenvector
          do m=1,5
            !
            do n=1,8
              !
              ! plus flux
              nwd=iwind8(k,n,kss,kee,'+')
              Flcp(m,n)=LEV(m,1)*Fswp(nwd,1)+                        &
                        LEV(m,2)*Fswp(nwd,2)+                        &
                        LEV(m,3)*Fswp(nwd,3)+                        &
                        LEV(m,4)*Fswp(nwd,4)+                        &
                        LEV(m,5)*Fswp(nwd,5)
              !
              ! minus flux
              nwd=iwind8(k,n,kss,kee,'-')
              Flcm(m,n)=LEV(m,1)*Fswm(nwd,1)+                        &
                        LEV(m,2)*Fswm(nwd,2)+                        &
                        LEV(m,3)*Fswm(nwd,3)+                        &
                        LEV(m,4)*Fswm(nwd,4)+                        &
                        LEV(m,5)*Fswm(nwd,5)
            end do
            !
          end do
          !
          if(numq>5) then
            do n=1,8
              nwd=iwind8(k,n,kss,kee,'+')
              Flcp(6:numq,n)=Fswp(nwd,6:numq)
              !
              nwd=iwind8(k,n,kss,kee,'-')
              Flcm(6:numq,n)=Fswm(nwd,6:numq)
            end do
          endif
          !
          ! End of characteristic decomposition.
          !
        else
          ! No characteristic decomposition
          !
          do n=1,8
            ! plus flux
            nwd=iwind8(k,n,kss,kee,'+')
            Flcp(1:numq,n)=Fswp(nwd,1:numq)
            !
            ! minus flux
            nwd=iwind8(k,n,kss,kee,'-')
            Flcm(1:numq,n)=Fswm(nwd,1:numq)
          end do
          !
        endif
        !
        ! Calculating values at i+1/2 using shock-capturing scheme.
        do m=1,numq
          !
          if(hdiss) then
            !
            var1=Flcp(m,4)
            var2=Flcm(m,4)
            !
          else
            !
            var1=recons_exp(   f=Flcp(m,:), inode=k,       &
                             dim=km,        ntype=npdck,   &
                             reschem=recon_schem,shock=lsh,solid=lso)
            var2=recons_exp(   f=Flcm(m,:), inode=k,       &
                             dim=km,        ntype=npdck,   &
                             reschem=recon_schem,shock=lsh,solid=lso)
            !
          endif
          !
          Fhc(m)=var1+var2
          !
        enddo
        !
        ! if(lchardecomp .and. lsh) then
        ! if(lchardecomp) then
        if(.false.) then
          do m=1,5
            Fh(k,m)=REV(m,1)*Fhc(1)+REV(m,2)*Fhc(2)+REV(m,3)*Fhc(3)+   &
                    REV(m,4)*Fhc(4)+REV(m,5)*Fhc(5) 
          end do
          !
          if(numq>5) then
            Fh(k,6:numq)=Fhc(6:numq)
          endif
        else
          Fh(k,1:numq)=Fhc(1:numq)
        endif
        !
      enddo
      
      do k=ks,ke
        do m=1,numq
          qrhs(i,j,k,m)=qrhs(i,j,k,m)+Fh(k,m)-Fh(k-1,m)
        enddo
      enddo
      !
    enddo
    enddo
    deallocate( Fswp,Fswm,Fh )
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! end of calculation at k direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    endif
    !
    if(present(timerept) .and. timerept) then
      !
      subtime=subtime+ptime()-time_beg
      !
      if(lio .and. lreport .and. ltimrpt) call timereporter(routine='convrsduwd', &
                                             timecost=subtime, &
                                              message='convection term using explicit upwind scheme')
    endif
    !  
    return
    !
  end subroutine convrsduwd
  !+-------------------------------------------------------------------+
  !| The end of the subroutine convrsduwd.                             |
  !+-------------------------------------------------------------------+
  !!
  !+-------------------------------------------------------------------+
  !| this function is to retune the i in a scheme's stencile.          |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 25-06-2022: Finished by J. Fang @ Warrington.                     |
  !+-------------------------------------------------------------------+
  pure integer function iwind8(i,n,imin,imax,dir)
    !
    integer,intent(in) :: i,n,imin,imax
    character(len=*),intent(in) :: dir
    !
    if(dir=='+') then
      !
      iwind8=i+n-4
      !
      if(iwind8<imin) iwind8=imin
      if(iwind8>imax) iwind8=imax
      !
    elseif(dir=='-') then
      !
      iwind8=i+5-n
      !
      if(iwind8<imin) iwind8=imin
      if(iwind8>imax) iwind8=imax
      !
    endif
    !
  end function iwind8
  !
  pure integer function iwind6(i,n,imin,imax,dir)
    !
    integer,intent(in) :: i,n,imin,imax
    character(len=*),intent(in) :: dir
    !
    if(dir=='+') then
      !
      iwind6=i+n-3
      !
      if(iwind6<imin) iwind6=imin
      if(iwind6>imax) iwind6=imax
      !
    elseif(dir=='-') then
      !
      iwind6=i+4-n
      !
      if(iwind6<imin) iwind6=imin
      if(iwind6>imax) iwind6=imax
      !
    endif
    !
  end function iwind6
  !+-------------------------------------------------------------------+
  !| The end of the function iwind.                                    |
  !+-------------------------------------------------------------------+
  !!
  !+-------------------------------------------------------------------+
  !| this subroutine is to solve the convectional term with compact    |
  !| MP schemes using Steger-Warming flux splitting scheme.            |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 22-03-2021: Created by J. Fang @ Warrington.                      |
  !| 15-02-2022: Finished by J. Fang @ Warrington.                     |
  !+-------------------------------------------------------------------+
  subroutine convrsdcmp(timerept)
    !
    use commvar,  only: im,jm,km,hm,numq,                              &
                        npdci,npdcj,npdck,is,ie,js,je,ks,ke,gamma,     &
                        recon_schem,lchardecomp,conschm
    use commarray,only: q,vel,rho,prs,tmp,spc,dxi,jacob,qrhs,lshock,crinod
    use commfunc, only: recons,mplimiter
    use comsolver, only : alfa_con,uci,ucj,uck,bci,bcj,bck
    use riemann,  only: flux_steger_warming
    !
    ! arguments
    logical,intent(in),optional :: timerept
    !
    ! local data
    integer :: i,j,k,iss,iee,jss,jee,kss,kee,nwd
    integer :: m,n,jvar
    real(8) :: eps,gm2,var1,var2
    !
    real(8) :: lmda(5),lmdap(5),lmdam(5),gpd(3),REV(5,5),LEV(5,5),     &
               Flcp(1:numq,1:5),Flcm(1:numq,1:5),Fhc(1:numq),          &
               fhcpc(1:numq),fhcmc(1:numq)
    !
    real(8), allocatable, dimension(:,:) :: fswp,fswm,fhcp,fhcm,Fh
    !
    real(8) :: time_beg
    real(8),save :: subtime=0.d0
    !
    logical :: lsh,sson,hdiss
    !
    if(present(timerept) .and. timerept) time_beg=ptime() 
    !
    sson=allocated(lshock)
    !
    eps=0.04d0
    gm2=0.5d0/gamma
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! calculating along i direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    allocate( fswp(-hm:im+hm,1:numq),fswm(-hm:im+hm,1:numq),           &
              fhcp(-1:im,1:numq),    fhcm(-1:im,1:numq) )
    allocate( Fh(-1:im,1:numq) )
    !
    if(npdci==1) then
      iss=0
      iee=im+hm
    elseif(npdci==2) then
      iss=-hm
      iee=im
    elseif(npdci==3) then
      iss=-hm
      iee=im+hm
    else
      stop ' !! error 1 @ subroutine convrsdcmp'
    endif
    !
    do k=ks,ke
    do j=js,je
      !
      ! flux split by using Steger-Warming method
      call flux_steger_warming( fplus= fswp(iss:iee,:),         &
                                fmius= fswm(iss:iee,:),         &
                                  rho=  rho(iss:iee,j,k),       &
                                  vel=  vel(iss:iee,j,k,:),     &
                                  prs=  prs(iss:iee,j,k),       &
                                  tmp=  tmp(iss:iee,j,k),       &
                                  spc=  spc(iss:iee,j,k,:),     &
                                    q=    q(iss:iee,j,k,:),     &
                                  dxi=  dxi(iss:iee,j,k,1,:),   &
                                jacob=jacob(iss:iee,j,k))
      ! End of flux split by using Steger-Warming method
      !
      ! calculating interface flux using compact upwind scheme ay i+1/2
      do jvar=1,numq
        fhcp(:,jvar)=recons(fswp(:,jvar),conschm,npdci,im,alfa_con,uci,windir='+')
        fhcm(:,jvar)=recons(fswm(:,jvar),conschm,npdci,im,alfa_con,bci,windir='-')
      enddo
      !
      ! Calculating Matrix of Left and Right eigenvectors using Roe 
      ! Average and flux split as well as i+1/2 construction.
      do i=is-1,ie
        !
        if(lchardecomp) then
          !
          call chardecomp(tmp(i,j,k),rho(i,j,k),prs(i,j,k),q(i,j,k,5),      &
                          vel(i,j,k,:),  spc(i,j,k,:),dxi(i,j,k,1,:),       &
                          tmp(i+1,j,k),rho(i+1,j,k),  prs(i+1,j,k),q(i+1,j,k,5),& 
                          vel(i+1,j,k,:),spc(i+1,j,k,:),dxi(i+1,j,k,1,:), &
                          REV,LEV)
          ! Project to characteristic space using local eigenvector
          do m=1,5
            !
            do n=1,5
              ! plus flux
              nwd=iwind6(i,n,iss,iee,'+')
              !
              flcp(m,n)=LEV(m,1)*Fswp(nwd,1)+                        &
                        LEV(m,2)*Fswp(nwd,2)+                        &
                        LEV(m,3)*Fswp(nwd,3)+                        &
                        LEV(m,4)*Fswp(nwd,4)+                        &
                        LEV(m,5)*Fswp(nwd,5)
              !
              ! minus flux
              nwd=iwind6(i,n,iss,iee,'-')
              flcm(m,n)=LEV(m,1)*Fswm(nwd,1)+                        &
                        LEV(m,2)*Fswm(nwd,2)+                        &
                        LEV(m,3)*Fswm(nwd,3)+                        &
                        LEV(m,4)*Fswm(nwd,4)+                        &
                        LEV(m,5)*Fswm(nwd,5)
            end do
            !
            fhcpc(m)=LEV(m,1)*fhcp(i,1)+                        &
                     LEV(m,2)*fhcp(i,2)+                        &
                     LEV(m,3)*fhcp(i,3)+                        &
                     LEV(m,4)*fhcp(i,4)+                        &
                     LEV(m,5)*fhcp(i,5)
            !
            fhcmc(m)=LEV(m,1)*fhcm(i,1)+                        &
                     LEV(m,2)*fhcm(i,2)+                        &
                     LEV(m,3)*fhcm(i,3)+                        &
                     LEV(m,4)*fhcm(i,4)+                        &
                     LEV(m,5)*fhcm(i,5)
            !
          end do
          !
          if(numq>5) then
            !
            do n=1,5
              ! plus flux
              nwd=iwind6(i,n,iss,iee,'+')
              flcp(6:numq,n)=Fswp(nwd,6:numq)
              !
              ! minus flux
              nwd=iwind6(i,n,iss,iee,'-')
              flcm(6:numq,n)=Fswm(nwd,6:numq)
            end do
            !
            fhcpc(6:numq)=fhcp(i,6:numq)
            fhcmc(6:numq)=fhcm(i,6:numq)
            !
          endif
          ! End of characteristic decomposition.
          !
        else
          ! No characteristic decomposition
          !
          do n=1,5
            ! plus flux
            nwd=iwind6(i,n,iss,iee,'+')
            flcp(1:numq,n)=Fswp(nwd,1:numq)
            !
            ! minus flux
            nwd=iwind6(i,n,iss,iee,'-')
            flcm(1:numq,n)=Fswm(nwd,1:numq)
          end do
          !
          fhcpc(1:numq)=fhcp(i,1:numq)
          fhcmc(1:numq)=fhcm(i,1:numq)
          !
        endif
        !
        ! Calculating values at i+1/2 using shock-capturing scheme.
        !
        if(sson) then
          !
          if(i<0) then
            lsh=lshock(i+1,j,k)
          elseif(i+1>im) then
            lsh=lshock(i,j,k)
          else
            lsh=lshock(i,j,k) .or. lshock(i+1,j,k)
          endif
          !
        else
          !
          lsh=.true.
          !
        endif
        !
        if(i<0) then
          hdiss=crinod(i+1,j,k)
        elseif(i+1>im) then
          hdiss=crinod(i,j,k)
        else
          hdiss=crinod(i,j,k) .or. crinod(i+1,j,k)
        endif
        !
        do m=1,numq
          !
          if(hdiss) then
            !
            var1=flcp(m,4)
            var2=flcm(m,4)
            !
          else
            !
            var1=mplimiter(flcp(m,1:5),fhcpc(m),shock=lsh,inode=i,   &
                                                   dim=im,ntype=npdci)
            var2=mplimiter(flcm(m,1:5),fhcmc(m),shock=lsh,inode=i,   &
                                                   dim=im,ntype=npdci)
            !
          endif
          !
          Fhc(m)=var1+var2
          !
        enddo
        !
        if(lchardecomp) then
          do m=1,5
            Fh(i,m)=REV(m,1)*Fhc(1)+REV(m,2)*Fhc(2)+REV(m,3)*Fhc(3)+   &
                    REV(m,4)*Fhc(4)+REV(m,5)*Fhc(5) 
          end do
          !
          if(numq>5) then
            Fh(i,6:numq)=Fhc(6:numq)
          endif
          !
        else
          Fh(i,1:numq)=Fhc(1:numq)
        endif
        !
      enddo
      
      do i=is,ie
        do m=1,numq
          qrhs(i,j,k,m)=qrhs(i,j,k,m)+Fh(i,m)-Fh(i-1,m)
        enddo
      enddo
      !
    enddo
    enddo
    deallocate( fswp,fswm,fhcp,fhcm,Fh )
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! end of calculation at i direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    if(ndims==1) return
    !
    if(npdcj==1) then
      jss=0
      jee=jm+hm
    elseif(npdcj==2) then
      jss=-hm
      jee=jm
    elseif(npdcj==3) then
      jss=-hm
      jee=jm+hm
    else
      stop ' !! error 2 @ subroutjne convrsdcmp'
    endif
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! ! calculating along j direction
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    allocate( fswp(-hm:jm+hm,1:numq),fswm(-hm:jm+hm,1:numq),           &
              fhcp(-1:jm,1:numq),    fhcm(-1:jm,1:numq) )
    allocate( Fh(-1:jm,1:numq) )
    !
    do k=ks,ke
    do i=is,ie
      !
      ! flux split by using Steger-Warming method
      call flux_steger_warming( fplus= fswp(jss:jee,:),         &
                                fmius= fswm(jss:jee,:),         &
                                  rho=  rho(i,jss:jee,k),       &
                                  vel=  vel(i,jss:jee,k,:),     &
                                  prs=  prs(i,jss:jee,k),       &
                                  tmp=  tmp(i,jss:jee,k),       &
                                  spc=  spc(i,jss:jee,k,:),     &
                                    q=    q(i,jss:jee,k,:),     &
                                  dxi=  dxi(i,jss:jee,k,2,:),   &
                                jacob=jacob(i,jss:jee,k))
      ! End of flux split by using Steger-Warming method
      !
      ! calculating interface flux using compact upwind scheme ay j+1/2
      do jvar=1,numq
        fhcp(:,jvar)=recons(fswp(:,jvar),conschm,npdcj,jm,alfa_con,ucj,windir='+')
        fhcm(:,jvar)=recons(fswm(:,jvar),conschm,npdcj,jm,alfa_con,bcj,windir='-')
      enddo
      !
      ! Calculating Matrix of Left and Right eigenvectors using Roe 
      ! Average and flux split as well as i+1/2 construction.
      do j=js-1,je
        !
        if(lchardecomp) then
          !
          call chardecomp(tmp(i,j,k),rho(i,j,k),prs(i,j,k),q(i,j,k,5),      &
                          vel(i,j,k,:),  spc(i,j,k,:),dxi(i,j,k,2,:),       &
                          tmp(i,j+1,k),rho(i,j+1,k),  prs(i,j+1,k),q(i,j+1,k,5), &
                          vel(i,j+1,k,:),spc(i,j+1,k,:),dxi(i,j+1,k,2,:),        &
                          REV,LEV)
          !
          ! Project to characteristic space using local eigenvector
          do m=1,5
            !
            do n=1,5
              ! plus flux 
              nwd=iwind6(j,n,jss,jee,'+')
              flcp(m,n)=LEV(m,1)*Fswp(nwd,1)+                        &
                        LEV(m,2)*Fswp(nwd,2)+                        &
                        LEV(m,3)*Fswp(nwd,3)+                        &
                        LEV(m,4)*Fswp(nwd,4)+                        &
                        LEV(m,5)*Fswp(nwd,5)
              !
              ! minus flux
              nwd=iwind6(j,n,jss,jee,'-')
              flcm(m,n)=LEV(m,1)*Fswm(nwd,1)+                        &
                        LEV(m,2)*Fswm(nwd,2)+                        &
                        LEV(m,3)*Fswm(nwd,3)+                        &
                        LEV(m,4)*Fswm(nwd,4)+                        &
                        LEV(m,5)*Fswm(nwd,5)
            end do
            !
            fhcpc(m)=LEV(m,1)*fhcp(j,1)+                        &
                     LEV(m,2)*fhcp(j,2)+                        &
                     LEV(m,3)*fhcp(j,3)+                        &
                     LEV(m,4)*fhcp(j,4)+                        &
                     LEV(m,5)*fhcp(j,5)
            !
            fhcmc(m)=LEV(m,1)*fhcm(j,1)+                        &
                     LEV(m,2)*fhcm(j,2)+                        &
                     LEV(m,3)*fhcm(j,3)+                        &
                     LEV(m,4)*fhcm(j,4)+                        &
                     LEV(m,5)*fhcm(j,5)
            !
          end do
          !
          if(numq>5) then
            !
            do n=1,5
              ! plus flux 
              nwd=iwind6(j,n,jss,jee,'+')
              flcp(6:numq,n)=Fswp(nwd,6:numq)
              !
              ! minus flux
              nwd=iwind6(j,n,jss,jee,'-')
              flcm(6:numq,n)=Fswm(nwd,6:numq)
            enddo
            !
            fhcpc(6:numq)=fhcp(j,6:numq)
            fhcmc(6:numq)=fhcm(j,6:numq)
            !
          endif
          ! End of characteristic decomposition.
          !
        else
          ! No characteristic decomposition
          !
          do n=1,5
            ! plus flux
            nwd=iwind6(j,n,jss,jee,'+')
            flcp(1:numq,n)=Fswp(nwd,1:numq)
            !
            ! minus flux
            nwd=iwind6(j,n,jss,jee,'-')
            flcm(1:numq,n)=Fswm(nwd,1:numq)
          end do
          !
          fhcpc(1:numq)=fhcp(j,1:numq)
          fhcmc(1:numq)=fhcm(j,1:numq)
          !
        endif
        !
        ! Calculating values at j+1/2 using shock-capturing scheme.
        !
        if(sson) then
          !
          if(j<0) then
            lsh=lshock(i,j+1,k)
          elseif(j+1>jm) then
            lsh=lshock(i,j,k)
          else
            lsh=lshock(i,j,k) .or. lshock(i,j+1,k)
          endif
          !
        else
          lsh=.true.
        endif
        !
        if(j<0) then
          hdiss=crinod(i,j+1,k)
        elseif(j+1>jm) then
          hdiss=crinod(i,j,k)
        else
          hdiss=crinod(i,j,k) .or. crinod(i,j+1,k)
        endif
        !
        do m=1,5
          !
          if(hdiss) then
            !
            var1=flcp(m,4)
            var2=flcm(m,4)
            !
          else
            !
            var1=mplimiter(flcp(m,1:5),fhcpc(m),shock=lsh,inode=j,   &
                                                   dim=jm,ntype=npdcj)
            var2=mplimiter(flcm(m,1:5),fhcmc(m),shock=lsh,inode=j,   &
                                                   dim=jm,ntype=npdcj)
            !
          endif
          !
          Fhc(m)=var1+var2
          !
        enddo
        !
        if(lchardecomp) then
          do m=1,5
            Fh(j,m)=REV(m,1)*Fhc(1)+REV(m,2)*Fhc(2)+REV(m,3)*Fhc(3)+   &
                    REV(m,4)*Fhc(4)+REV(m,5)*Fhc(5) 
          end do
          !
          if(numq>5) then
            Fh(j,6:numq)=Fhc(6:numq)
          endif
          !
        else
          Fh(j,1:numq)=Fhc(1:numq)
        endif
        !
        !
      enddo
      !
      do j=js,je
        do m=1,numq
          qrhs(i,j,k,m)=qrhs(i,j,k,m)+Fh(j,m)-Fh(j-1,m)
        enddo
      enddo
      !
    enddo
    enddo
    deallocate( fswp,fswm,fhcp,fhcm,Fh )
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! ! end of calculation at j direction
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    if(ndims==2) return
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! ! calculating along k direction
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    if(npdck==1) then
      kss=0
      kee=km+hm
    elseif(npdck==2) then
      kss=-hm
      kee=km
    elseif(npdck==3) then
      kss=-hm
      kee=km+hm
    else
      stop ' !! error 2 @ subroutine convrsdcmp'
    endif
    !
    allocate( fswp(-hm:km+hm,1:numq),fswm(-hm:km+hm,1:numq),           &
              fhcp(-1:km,1:numq),    fhcm(-1:km,1:numq) )
    allocate( Fh(-1:km,1:numq) )
    !
    do j=js,je
    do i=is,ie
      !
      ! flux split by using Steger-Warming method
      call flux_steger_warming( fplus= fswp(kss:kee,:),         &
                                fmius= fswm(kss:kee,:),         &
                                  rho=  rho(i,j,kss:kee),       &
                                  vel=  vel(i,j,kss:kee,:),     &
                                  prs=  prs(i,j,kss:kee),       &
                                  tmp=  tmp(i,j,kss:kee),       &
                                  spc=  spc(i,j,kss:kee,:),     &
                                    q=    q(i,j,kss:kee,:),     &
                                  dxi=  dxi(i,j,kss:kee,3,:),   &
                                jacob=jacob(i,j,kss:kee))
      ! End of flux split by using Steger-Warming method
      !
      ! calculating interface flux using compact upwind scheme ay k+1/2
      do jvar=1,numq
        fhcp(:,jvar)=recons(fswp(:,jvar),conschm,npdck,km,alfa_con,uck,windir='+')
        fhcm(:,jvar)=recons(fswm(:,jvar),conschm,npdck,km,alfa_con,bck,windir='-')
      enddo
      !
      ! Calculating Matrix of Left and Right eigenvectors using Roe 
      ! Average and flux split as well as k+1/2 construction.
      do k=ks-1,ke
        !
        if(lchardecomp) then
          !
          call chardecomp(tmp(i,j,k),rho(i,j,k),    prs(i,j,k),  q(i,j,k,5),      &
                          vel(i,j,k,:),  spc(i,j,k,:),dxi(i,j,k,3,:),             &
                          tmp(i,j,k+1),rho(i,j,k+1),  prs(i,j,k+1),q(i,j,k+1,5),  &
                          vel(i,j,k+1,:),spc(i,j,k+1,:),dxi(i,j,k+1,3,:),REV,LEV)
          !
          ! Project to characteristic space using local eigenvector
          do m=1,5
            !
            do n=1,5
              ! plus flux
              nwd=iwind6(k,n,kss,kee,'+')
              flcp(m,n)=LEV(m,1)*Fswp(nwd,1)+                        &
                        LEV(m,2)*Fswp(nwd,2)+                        &
                        LEV(m,3)*Fswp(nwd,3)+                        &
                        LEV(m,4)*Fswp(nwd,4)+                        &
                        LEV(m,5)*Fswp(nwd,5)
              !
              ! minus flux
              nwd=iwind6(k,n,kss,kee,'-')
              flcm(m,n)=LEV(m,1)*Fswm(nwd,1)+                        &
                        LEV(m,2)*Fswm(nwd,2)+                        &
                        LEV(m,3)*Fswm(nwd,3)+                        &
                        LEV(m,4)*Fswm(nwd,4)+                        &
                        LEV(m,5)*Fswm(nwd,5)
            end do
            !
            fhcpc(m)=LEV(m,1)*fhcp(k,1)+                        &
                     LEV(m,2)*fhcp(k,2)+                        &
                     LEV(m,3)*fhcp(k,3)+                        &
                     LEV(m,4)*fhcp(k,4)+                        &
                     LEV(m,5)*fhcp(k,5)
            !
            fhcmc(m)=LEV(m,1)*fhcm(k,1)+                        &
                     LEV(m,2)*fhcm(k,2)+                        &
                     LEV(m,3)*fhcm(k,3)+                        &
                     LEV(m,4)*fhcm(k,4)+                        &
                     LEV(m,5)*fhcm(k,5)
            !
          end do
          !
          if(numq>5) then
            !
            do n=1,5
              ! plus flux
              nwd=iwind6(k,n,kss,kee,'+')
              flcp(6:numq,n)=Fswp(nwd,6:numq)
              !
              ! minus flux
              nwd=iwind6(k,n,kss,kee,'-')
              flcm(6:numq,n)=Fswm(nwd,6:numq)
            enddo
            !
            fhcpc(6:numq)=fhcp(k,6:numq)
            fhcmc(6:numq)=fhcm(k,6:numq)
            !
          endif
          ! End of characteristic decomposition.
          !
        else
          ! No characteristic decomposition
          !
          do n=1,5
            ! plus flux
            nwd=iwind6(k,n,kss,kee,'+')
            flcp(1:numq,n)=Fswp(nwd,1:numq)
            !
            ! minus flux
            nwd=iwind6(k,n,kss,kee,'-')
            flcm(1:numq,n)=Fswm(nwd,1:numq)
          end do
          !
          fhcpc(1:numq)=fhcp(k,1:numq)
          fhcmc(1:numq)=fhcm(k,1:numq)
          !
        endif
        !
        ! Calculating values at k+1/2 using shock-capturing scheme.
        !
        if(sson) then
          !
          if(k<0) then
            lsh=lshock(i,j,k+1)
          elseif(k+1>km) then
            lsh=lshock(i,j,k)
          else
            lsh=lshock(i,j,k) .or. lshock(i,j,k+1)
          endif
          !
        else
          lsh=.true.
        endif
        !
        if(k<0) then
          hdiss=crinod(i,j,k+1)
        elseif(k+1>km) then
          hdiss=crinod(i,j,k)
        else
          hdiss=crinod(i,j,k) .or. crinod(i,j,k+1)
        endif
        !
        do m=1,5
          !
          if(hdiss) then
            !
            var1=flcp(m,4)
            var2=flcm(m,4)
            !
          else
            !
            var1=mplimiter(flcp(m,1:5),fhcpc(m),shock=lsh,inode=k,   &
                                                   dim=km,ntype=npdck)
            var2=mplimiter(flcm(m,1:5),fhcmc(m),shock=lsh,inode=k,   &
                                                   dim=km,ntype=npdck)
            !
          endif
          !
          Fhc(m)=var1+var2
          !
        enddo
        !
        if(lchardecomp) then
          !
          do m=1,5
            Fh(k,m)=REV(m,1)*Fhc(1)+REV(m,2)*Fhc(2)+REV(m,3)*Fhc(3)+   &
                    REV(m,4)*Fhc(4)+REV(m,5)*Fhc(5) 
          end do
          !
          if(numq>5) then
            Fh(k,6:numq)=Fhc(6:numq)
          endif
          !
        else
          Fh(k,1:numq)=Fhc(1:numq)
        endif
        !
      enddo
      !
      do k=ks,ke
        do m=1,numq
          qrhs(i,j,k,m)=qrhs(i,j,k,m)+Fh(k,m)-Fh(k-1,m)
        enddo
      enddo
      !
    enddo
    enddo
    !
    deallocate( fswp,fswm,fhcp,fhcm,Fh )
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! ! end of calculation at k direction
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    if(present(timerept) .and. timerept) then
      !
      subtime=subtime+ptime()-time_beg
      !
      if(lio .and. lreport .and. ltimrpt) call timereporter(routine='convrsdcmp', &
                                             timecost=subtime, &
                                              message='convection term using compact upwind scheme')
    endif
    !
    return
    !
  end subroutine convrsdcmp
  !+-------------------------------------------------------------------+
  !| The end of the subroutine convrsdcmp.                             |
  !+-------------------------------------------------------------------+
  !!
  !+-------------------------------------------------------------------+
  !| this subroutine is to return left or right matrix of              |
  !| characteristc vectpr                                              |
  !+-------------------------------------------------------------------+
  !| ref1: J. Wang, S. Pan, X. Y. Hu, and N. A. Adams, Partial         |
  !| characteristic decomposition for multi-species Euler equations,   |
  !| Comput. Fluids 181, 364, (2019).                                  |
  !| ref2: R. P. Fedkiw, B. Merriman, and S. Osher, High accuracy      |
  !| numerical methods for thermally perfect gas flows with chemistry, |
  !| Journal of Computational Physics 132, 175-190, (1997).            |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 22-03-2021: Created by J. Fang @ Warrington.                      |
  !| 28-02-2023: Modified to adapt multi-species equations.            |
  !+-------------------------------------------------------------------+
  subroutine chardecomp(tmp_l,ro_l,p_l,E_l,vel_l,spc_l,ddi_l,             &
                        tmp_r,ro_r,p_r,E_r,vel_r,spc_r,ddi_r,REV,LEV,debug)
    !
    use commvar, only: gamma,num_species
#ifdef COMB
    use thermchem,only:aceval,temperature_calc,gammarmix
#endif
    !
    real(8),intent(in) :: tmp_l,ro_l,p_l,E_l,vel_l(3),ddi_l(3),       &
                          tmp_r,ro_r,p_r,E_r,vel_r(3),ddi_r(3)
    real(8),intent(in) :: spc_l(:),spc_r(:)
    real(8),intent(out) :: REV(5,5),LEV(5,5)
    !
    logical,optional,intent(in) :: debug
    !
    ! local data
    real(8),parameter :: rero=1.d-12
    real(8) :: WRoe,WRoe1,u1Roe,u2Roe,u3Roe,HL,HR,HRoe,CssRoe,        &
               KRoe,ugp,rcs,var1,var2,var3,var4,rgp,gpd(3),           &
               ee_L,ee_R,eeRoe,tmpRoe,roRoe,b1,b2
    real(8) :: spcRoe(num_species),gamavg,gamL,gamR,CssL,CssR,KL,KR
    logical :: ldebug
    !
    if(present(debug)) then
      ldebug=debug
    else
      ldebug=.false.
    endif
    !
    roRoe=sqrt(ro_l*ro_R)
    !
    WRoe=sqrt(ro_l)/(sqrt(ro_l)+sqrt(ro_r))
    WRoe1=1.d0-WRoe
    u1Roe=WRoe*vel_l(1)+WRoe1*vel_r(1)
    u2Roe=WRoe*vel_l(2)+WRoe1*vel_r(2)
    u3Roe=WRoe*vel_l(3)+WRoe1*vel_r(3)
    !
    KRoe=0.5d0*(u1Roe*u1Roe+u2Roe*u2Roe+u3Roe*u3Roe)
    KL=0.5d0*(vel_l(1)*vel_l(1)+vel_l(2)*vel_l(2)+vel_l(3)*vel_l(3))
    KR=0.5d0*(vel_r(1)*vel_r(1)+vel_r(2)*vel_r(2)+vel_r(3)*vel_r(3))
    !
#ifdef COMB
    gamL = gammarmix(tmp_l,spc_l(:))
    gamR = gammarmix(tmp_r,spc_r(:))
    gamavg=0.5d0*(gamL+gamR)
    !
    call aceval(tmp_l,spc_l(:),CssL)
    call aceval(tmp_R,spc_R(:),CssR)

    HL=CssL*CssL/(gamL-1.d0)+KL
    HR=CssR*CssR/(gamR-1.d0)+KR
    HRoe=WRoe*HL+WRoe1*HR

    CssRoe=sqrt((gamavg-1.d0)*(HRoe-KRoe))
#else
    HL=(E_l+p_l)/ro_l
    HR=(E_r+p_r)/ro_r
    HRoe=WRoe*HL+WRoe1*HR
    CssRoe=sqrt((gamma-1.d0)*(HRoe-KRoe))
#endif
    !
    rcs=1.d0/CssRoe
    !
    ! print*,' ** CssRoe:',CssRoe
    !
    var1=0.5d0*(ddi_l(1)+ddi_r(1))
    var2=0.5d0*(ddi_l(2)+ddi_r(2))
    var3=0.5d0*(ddi_l(3)+ddi_r(3))
    var4=1.d0/sqrt(var1*var1+var2*var2+var3*var3)
    !
    gpd(1)=var1*var4
    gpd(2)=var2*var4
    gpd(3)=var3*var4
    !
    ugp=u1Roe*gpd(1)+u2Roe*gpd(2)+u3Roe*gpd(3)
    !
    ! Calculating Left and Right eigenvectors
    b1=(gamma-1.d0)/(CssRoe*CssRoe)
    ! b2=b1*KRoe
    b2=1.d0+2.d0*b1*KRoe-b1*HRoe
    !
    LEV(1,1)= 0.5d0*(b2+ugp*rcs)
    LEV(1,2)=-0.5d0*(b1*u1Roe+gpd(1)*rcs)
    LEV(1,3)=-0.5d0*(b1*u2Roe+gpd(2)*rcs)
    LEV(1,4)=-0.5d0*(b1*u3Roe+gpd(3)*rcs)
    LEV(1,5)= 0.5d0*b1
    !
    REV(1,1)=1.d0
    REV(2,1)=u1Roe-CssRoe*gpd(1)
    REV(3,1)=u2Roe-CssRoe*gpd(2)
    REV(4,1)=u3Roe-CssRoe*gpd(3)
    REV(5,1)=HRoe-ugp*CssRoe
    !
    LEV(2,1)=1.d0-b2
    LEV(2,2)=b1*u1Roe
    LEV(2,3)=b1*u2Roe
    LEV(2,4)=b1*u3Roe
    LEV(2,5)=-b1
    !
    REV(1,2)=1.d0
    REV(2,2)=u1Roe
    REV(3,2)=u2Roe
    REV(4,2)=u3Roe
    REV(5,2)=HRoe-1.d0/b1
    !
    if(abs(var1)>rero) then
      rgp=1.d0/gpd(1)
      !
      LEV(3,1)=(ugp*gpd(2)-u2Roe)*rgp
      LEV(3,2)=-gpd(2)
      LEV(3,3)=(1.d0-gpd(2)*gpd(2))*rgp
      LEV(3,4)=-gpd(2)*gpd(3)*rgp
      LEV(3,5)=0.d0
      !
      LEV(4,1)=(ugp*gpd(3)-u3Roe)*rgp
      LEV(4,2)=-gpd(3)
      LEV(4,3)=-gpd(2)*gpd(3)*rgp
      LEV(4,4)=(1.d0-gpd(3)*gpd(3))*rgp
      LEV(4,5)=0.d0
      !
      REV(1,3)=0.d0
      REV(2,3)=-gpd(2)
      REV(3,3)= gpd(1)
      REV(4,3)=0.d0
      REV(5,3)=u2Roe*gpd(1)-u1Roe*gpd(2)
      !
      REV(1,4)=0.d0
      REV(2,4)=-gpd(3)
      REV(3,4)= 0.d0
      REV(4,4)= gpd(1)
      REV(5,4)=u3Roe*gpd(1)-u1Roe*gpd(3)
      !
    elseif(abs(var2)>rero) then
      rgp=1.d0/gpd(2)
      !
      LEV(3,1)=(ugp*gpd(1)-u1Roe)*rgp
      LEV(3,2)=(1.d0-gpd(1)*gpd(1))*rgp
      LEV(3,3)=-gpd(1)
      LEV(3,4)=-gpd(1)*gpd(3)*rgp
      LEV(3,5)=0.d0
      !
      LEV(4,1)=(ugp*gpd(3)-u3Roe)*rgp
      LEV(4,2)=-gpd(1)*gpd(3)*rgp
      LEV(4,3)=-gpd(3)
      LEV(4,4)=(1.d0-gpd(3)*gpd(3))*rgp
      LEV(4,5)=0.d0
      !
      REV(1,3)=0.d0
      REV(2,3)= gpd(2)
      REV(3,3)=-gpd(1)
      REV(4,3)=0.d0
      REV(5,3)=u1Roe*gpd(2)-u2Roe*gpd(1)
      !
      REV(1,4)=0.d0
      REV(2,4)=0.d0
      REV(3,4)=-gpd(3)
      REV(4,4)= gpd(2)
      REV(5,4)=u3Roe*gpd(2)-u2Roe*gpd(3)
      !
    elseif(abs(var3)>rero) then
      rgp=1.d0/gpd(3)
      !
      LEV(3,1)=(ugp*gpd(1)-u1Roe)*rgp
      LEV(3,2)=(1.d0-gpd(1)*gpd(1))*rgp
      LEV(3,3)=-gpd(1)*gpd(2)*rgp
      LEV(3,4)=-gpd(1)
      LEV(3,5)=0.d0
      !
      LEV(4,1)=(ugp*gpd(2)-u2Roe)*rgp
      LEV(4,2)=-gpd(1)*gpd(2)*rgp
      LEV(4,3)=(1.d0-gpd(2)*gpd(2))*rgp
      LEV(4,4)=-gpd(2)
      LEV(4,5)=0.d0
      !
      REV(1,3)=0.d0
      REV(2,3)= gpd(3)
      REV(3,3)=0.d0
      REV(4,3)=-gpd(1)
      REV(5,3)=u1Roe*gpd(3)-u3Roe*gpd(1)
      !
      REV(1,4)=0.d0
      REV(2,4)=0.d0
      REV(3,4)= gpd(3)
      REV(4,4)=-gpd(2)
      REV(5,4)=u2Roe*gpd(3)-u3Roe*gpd(2)
      !
    else
      stop ' !! ERROR 1 @ chardecomp'
    endif
    !
    LEV(5,1)= 0.5d0*(b2-ugp*rcs)           
    LEV(5,2)=-0.5d0*(b1*u1Roe-gpd(1)*rcs)  
    LEV(5,3)=-0.5d0*(b1*u2Roe-gpd(2)*rcs)  
    LEV(5,4)=-0.5d0*(b1*u3Roe-gpd(3)*rcs)  
    LEV(5,5)= 0.5d0*b1 
    !
    REV(1,5)= 1.d0
    REV(2,5)=u1Roe+CssRoe*gpd(1)
    REV(3,5)=u2Roe+CssRoe*gpd(2)
    REV(4,5)=u3Roe+CssRoe*gpd(3)
    REV(5,5)=HRoe+ugp*CssRoe
    !
    return
    !
  end subroutine chardecomp
  !+-------------------------------------------------------------------+
  !| The end of the subroutine chardecomp.                             |
  !+-------------------------------------------------------------------+
  !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! This subroutine is used to calculate the convectional residual
  ! terms with compact six-order central scheme.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Writen by Fang Jian, 2009-06-03.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine convrsdcal6(timerept)
    !
    use commvar,  only: im,jm,km,hm,numq,num_species,num_modequ,       &
                        conschm,npdci,npdcj,npdck,is,ie,js,je,ks,ke
    use commarray,only: q,vel,rho,prs,tmp,spc,dxi,jacob,qrhs
    use derivative,only : fds
    use comsolver, only : alfa_con,cci,ccj,cck
    !
    ! arguments
    logical,intent(in),optional :: timerept
    !
    ! local data
    integer :: i,j,k,jspc,jmod,n
    real(8),allocatable :: fcs(:,:),dfcs(:,:),uu(:)
    !
    real(8) :: time_beg
    real(8),save :: subtime=0.d0
    !
    if(present(timerept) .and. timerept) time_beg=ptime() 
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! calculating along i direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    allocate(fcs(-hm:im+hm,1:numq),dfcs(0:im,1:numq),uu(-hm:im+hm))
    do k=ks,ke
    do j=js,je
      !
      uu(:)=dxi(:,j,k,1,1)*vel(:,j,k,1)+dxi(:,j,k,1,2)*vel(:,j,k,2) +  &
            dxi(:,j,k,1,3)*vel(:,j,k,3)
      fcs(:,1)=jacob(:,j,k)*  q(:,j,k,1)*uu
      fcs(:,2)=jacob(:,j,k)*( q(:,j,k,2)*uu+dxi(:,j,k,1,1)*prs(:,j,k) )
      fcs(:,3)=jacob(:,j,k)*( q(:,j,k,3)*uu+dxi(:,j,k,1,2)*prs(:,j,k) )
      fcs(:,4)=jacob(:,j,k)*( q(:,j,k,4)*uu+dxi(:,j,k,1,3)*prs(:,j,k) )
      fcs(:,5)=jacob(:,j,k)*( q(:,j,k,5)+prs(:,j,k) )*uu
      !
      ! do i=0,im
      !   print*,prs(i,j,k)
      ! enddo
      !
      if(num_species>0) then
        n=5
        do jspc=1,num_species
          fcs(:,n+jspc)=jacob(:,j,k)*q(:,j,k,n+jspc)*uu
        enddo
      endif
      !
      if(num_modequ>0) then
        n=5+num_species
        do jmod=1,num_modequ
          fcs(:,n+jmod)=jacob(:,j,k)*q(:,j,k,n+jmod)*uu
        enddo
      endif
      !
      do n=1,numq
        dfcs(:,n)=fds%central(fcs(:,n),ntype=npdci,dim=im,dir=1)
      enddo
      !
      qrhs(is:ie,j,k,:)=qrhs(is:ie,j,k,:)+dfcs(is:ie,:) 
      !
    enddo
    enddo
    deallocate(fcs,dfcs,uu)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! end calculating along i direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    if(ndims>=2) then
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! calculating along j direction
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      allocate(fcs(-hm:jm+hm,1:numq),dfcs(0:jm,1:numq),uu(-hm:jm+hm))
      do k=ks,ke
      do i=is,ie
        !
        uu(:)=dxi(i,:,k,2,1)*vel(i,:,k,1)+dxi(i,:,k,2,2)*vel(i,:,k,2) +  &
              dxi(i,:,k,2,3)*vel(i,:,k,3)
        fcs(:,1)=jacob(i,:,k)*  q(i,:,k,1)*uu
        fcs(:,2)=jacob(i,:,k)*( q(i,:,k,2)*uu+dxi(i,:,k,2,1)*prs(i,:,k) )
        fcs(:,3)=jacob(i,:,k)*( q(i,:,k,3)*uu+dxi(i,:,k,2,2)*prs(i,:,k) )
        fcs(:,4)=jacob(i,:,k)*( q(i,:,k,4)*uu+dxi(i,:,k,2,3)*prs(i,:,k) )
        fcs(:,5)=jacob(i,:,k)*( q(i,:,k,5)+prs(i,:,k) )*uu
        !
        if(num_species>0) then
          n=5
          do jspc=1,num_species
            fcs(:,n+jspc)=jacob(i,:,k)*q(i,:,k,n+jspc)*uu
          enddo
        endif
        !
        if(num_modequ>0) then
          n=5+num_species
          do jmod=1,num_modequ
            fcs(:,n+jmod)=jacob(i,:,k)*q(i,:,k,n+jmod)*uu
          enddo
        endif
        !
        do n=1,numq
          dfcs(:,n)=fds%central(fcs(:,n),ntype=npdcj,dim=jm,dir=2)
        enddo
        !
        qrhs(i,js:je,k,:)=qrhs(i,js:je,k,:)+dfcs(js:je,:)
        !
        !
      enddo
      enddo
      deallocate(fcs,dfcs,uu)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! end calculating along j direction
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    endif
    !
    if(ndims==3) then
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! calculating along j direction
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      allocate(fcs(-hm:km+hm,1:numq),dfcs(0:km,1:numq),uu(-hm:km+hm))
      do j=js,je
      do i=is,ie
        !
        uu(:)=dxi(i,j,:,3,1)*vel(i,j,:,1)+dxi(i,j,:,3,2)*vel(i,j,:,2) +  &
              dxi(i,j,:,3,3)*vel(i,j,:,3)
        fcs(:,1)=jacob(i,j,:)*  q(i,j,:,1)*uu
        fcs(:,2)=jacob(i,j,:)*( q(i,j,:,2)*uu+dxi(i,j,:,3,1)*prs(i,j,:) )
        fcs(:,3)=jacob(i,j,:)*( q(i,j,:,3)*uu+dxi(i,j,:,3,2)*prs(i,j,:) )
        fcs(:,4)=jacob(i,j,:)*( q(i,j,:,4)*uu+dxi(i,j,:,3,3)*prs(i,j,:) )
        fcs(:,5)=jacob(i,j,:)*( q(i,j,:,5)+prs(i,j,:) )*uu
        !
        if(num_species>0) then
          n=5
          do jspc=1,num_species
            fcs(:,n+jspc)=jacob(i,j,:)*q(i,j,:,n+jspc)*uu
          enddo
        endif
        !
        if(num_modequ>0) then
          n=5+num_species
          do jmod=1,num_modequ
            fcs(:,n+jmod)=jacob(i,j,:)*q(i,j,:,n+jmod)*uu
          enddo
        endif
        !
        do n=1,numq
          dfcs(:,n)=fds%central(fcs(:,n),ntype=npdck,dim=km,dir=3)
        enddo
        !
        qrhs(i,j,ks:ke,:)=qrhs(i,j,ks:ke,:)+dfcs(ks:ke,:)
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
      if(lio .and. lreport .and. ltimrpt) call timereporter(routine='convrsdcal6', &
                                             timecost=subtime, &
                                              message='diffusion term with central scheme')
    endif
    !
    return
    !
  end subroutine convrsdcal6
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! End of the subroutine ConvRsdCal6.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! This subroutine is used to calculate the diffusion term with 6-order
  ! Compact Central scheme.
  !   sixth-order Compact Central scheme.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Writen by Fang Jian, 2009-06-09.
  ! Add scalar transport equation by Fang Jian, 2022-01-12.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine diffrsdcal6(timerept)
    !
    use commvar,   only : im,jm,km,numq,npdci,npdcj,npdck,difschm,     &
                          conschm,ndims,num_species,num_modequ,        &
                          reynolds,prandtl,const5,is,ie,js,je,ks,ke,   &
                          turbmode,nondimen,schmidt,nstep,deltat,      &
                          cp,flowtype
    use commarray, only : vel,tmp,spc,dvel,dtmp,dspc,dxi,x,jacob,qrhs, &
                          rho,vor,omg,tke,miut,dtke,domg,res12
    use derivative, only : fds
    use comsolver, only : alfa_dif,dci,dcj,dck
    use fludyna,   only : miucal
    use models,    only : komega,src_komega
    use tecio
    use parallel,  only : yflux_sendrecv
#ifdef COMB
    use thermchem, only : tranmod,tranco,enthpy,convertxiyi,wmolar
#endif
    !
    ! arguments
    logical,intent(in),optional :: timerept
    !
    ! local data
    integer :: i,j,k,n,ncolm,jspc,idir
    real(8),allocatable :: df(:,:),ff(:,:)
    real(8),allocatable,dimension(:,:,:,:),  save :: sigma,qflux,dkflux,doflux
    real(8),allocatable,dimension(:,:,:,:,:),save :: yflux
    real(8) :: miu,miu2,miu3,miu4,hcc,s11,s12,s13,s22,s23,s33,skk
    real(8) :: d11,d12,d13,d21,d22,d23,d31,d32,d33,miueddy,var1,var2
    real(8) :: tau11,tau12,tau13,tau22,tau23,tau33
    real(8) :: detk
    real(8),allocatable :: dispec(:,:)
    real(8) :: corrdiff,hispec(num_species),xi(num_species),cpe,kama, &
               gradyi(num_species),sum1,sum2,mw
    real(8),allocatable :: dfu(:)
    logical,save :: firstcall=.true.
    !
    real(8) :: time_beg
    real(8),save :: subtime=0.d0
    !
    if(present(timerept) .and. timerept) time_beg=ptime() 
    !
    if(firstcall) then
      allocate( sigma(-hm:im+hm,-hm:jm+hm,-hm:km+hm,1:6),                &
                qflux(-hm:im+hm,-hm:jm+hm,-hm:km+hm,1:3) )
    endif
    !
    sigma =0.d0
    qflux =0.d0
    !
    if(num_species>0) then
      !
      if(firstcall) then
        allocate( yflux(-hm:im+hm,-hm:jm+hm,-hm:km+hm,1:num_species,1:3) )
      endif
      !
#ifndef  COMB
      allocate(dfu(1:num_species))
#endif
      yflux=0.d0
      !
    endif
    !
    if(trim(turbmode)=='k-omega') then
      !
      if(firstcall) then
        allocate( dkflux(-hm:im+hm,-hm:jm+hm,-hm:km+hm,1:3),             &
                  doflux(-hm:im+hm,-hm:jm+hm,-hm:km+hm,1:3) )
      endif
      !
      dkflux=0.d0
      doflux=0.d0
      !
    endif
    !
    if(firstcall) firstcall=.false.
    !
    if(trim(turbmode)=='k-omega') call src_komega
    !
    do k=0,km
    do j=0,jm
    do i=0,im
      !
#ifdef COMB

     call enthpy(tmp(i,j,k),hispec(:))
     call convertxiyi(spc(i,j,k,:),xi(:),'Y2X')
     mw=sum(wmolar(:)*xi(:))
     !
     select case(tranmod)
       case('multi')
         if(.not.allocated(dispec))allocate(dispec(num_species,num_species))
         call tranco(den=rho(i,j,k),tmp=tmp(i,j,k),cp=cpe,mu=miu,lam=kama, &
                     spc=spc(i,j,k,:),rhodij=dispec(:,:))
     case default
       if(.not.allocated(dispec)) allocate(dispec(num_species,1))
       call tranco(den=rho(i,j,k),tmp=tmp(i,j,k),cp=cpe,mu=miu,lam=kama, &
                   spc=spc(i,j,k,:),rhodi=dispec(:,1))
       ! print*,' ** miu=',miu,'cp=',cpe,'kamma=',kama,'rhodi=',dispec(:,1)
     end select

#else
      if(nondimen) then 
        miu=miucal(tmp(i,j,k))/reynolds
      else
        miu=miucal(tmp(i,j,k))
      endif 
#endif
      !
      s11=dvel(i,j,k,1,1)
      s12=0.5d0*(dvel(i,j,k,1,2)+dvel(i,j,k,2,1))
      s13=0.5d0*(dvel(i,j,k,1,3)+dvel(i,j,k,3,1))
      s22=dvel(i,j,k,2,2)
      s23=0.5d0*(dvel(i,j,k,2,3)+dvel(i,j,k,3,2))
      s33=dvel(i,j,k,3,3)
      !
      skk=num1d3*(s11+s22+s33)
      !
      vor(i,j,k,1)=dvel(i,j,k,3,2)-dvel(i,j,k,2,3)
      vor(i,j,k,2)=dvel(i,j,k,1,3)-dvel(i,j,k,3,1)
      vor(i,j,k,3)=dvel(i,j,k,2,1)-dvel(i,j,k,1,2)
      !
      tau11=0.d0
      tau12=0.d0
      tau13=0.d0
      tau22=0.d0
      tau23=0.d0
      tau33=0.d0
      !
      if(trim(turbmode)=='k-omega') then
        !
        miu2=2.d0*(miu+miut(i,j,k))
        !
        hcc=(miu/prandtl+miut(i,j,k)/komega%prt)/const5
        !
        detk=num2d3*rho(i,j,k)*tke(i,j,k)
        !
        miu3=miu+komega%sigma_k(i,j,k)    *miut(i,j,k)
        miu4=miu+komega%sigma_omega(i,j,k)*miut(i,j,k)
        !
        dkflux(i,j,k,:)=miu3*dtke(i,j,k,:)
        doflux(i,j,k,:)=miu4*domg(i,j,k,:)
      elseif(trim(turbmode)=='udf1') then
        !
        ! miu2=2.d0*(miu+miut(i,j,k))
        ! hcc=(miu/prandtl+miut(i,j,k)/0.9d0)/const5
        !
        miu2=2.d0*miu
        hcc=(miu/prandtl)/const5
        !
        tau11=0.d0
        tau12=-res12(i,j,k)*rho(i,j,k)
        tau13=0.d0
        tau22=0.d0
        tau23=0.d0
        tau33=0.d0
        !
        detk=0.d0
        !
      elseif(trim(turbmode)=='none') then
        miu2=2.d0*miu
        !
#ifdef COMB
        hcc=kama
#else
        if(nondimen) then 
          hcc=(miu/prandtl)/const5
        else
          hcc=cp*miu/prandtl
        endif 
#endif
        !
        detk=0.d0
        !
      endif
      !
      sigma(i,j,k,1)=miu2*(s11-skk)-detk + tau11 !s11   
      sigma(i,j,k,2)=miu2* s12           + tau12 !s12  
      sigma(i,j,k,3)=miu2* s13           + tau13 !s13   
      sigma(i,j,k,4)=miu2*(s22-skk)-detk + tau22 !s22   
      sigma(i,j,k,5)=miu2* s23           + tau23 !s23  
      sigma(i,j,k,6)=miu2*(s33-skk)-detk + tau33 !s33  
      !
      qflux(i,j,k,1)=hcc*dtmp(i,j,k,1)+sigma(i,j,k,1)*vel(i,j,k,1) +   &
                                       sigma(i,j,k,2)*vel(i,j,k,2) +   &
                                       sigma(i,j,k,3)*vel(i,j,k,3)
      qflux(i,j,k,2)=hcc*dtmp(i,j,k,2)+sigma(i,j,k,2)*vel(i,j,k,1) +   &
                                       sigma(i,j,k,4)*vel(i,j,k,2) +   &
                                       sigma(i,j,k,5)*vel(i,j,k,3)
      qflux(i,j,k,3)=hcc*dtmp(i,j,k,3)+sigma(i,j,k,3)*vel(i,j,k,1) +   &
                                       sigma(i,j,k,5)*vel(i,j,k,2) +   &
                                       sigma(i,j,k,6)*vel(i,j,k,3)
      !                                      
      if(num_species>0) then
        !
#ifdef COMB
        !
        do idir=1,3
          !
          gradyi(:)=dspc(i,j,k,:,idir)
          !
          sum1=mw*sum(gradyi(:)/wmolar(:))
          !
          select case(tranmod)
            !
            case('multi')
              !
              do jspc=1,num_species
                sum2=sum(dispec(jspc,:)*(gradyi(:)-spc(i,j,k,:)*sum1))
                !species diffusive flux
                yflux(i,j,k,jspc,idir)=-1.d0*wmolar(jspc)/mw*sum2
                !energy flux due to species diffusion
                qflux(i,j,k,idir)=qflux(i,j,k,idir)-yflux(i,j,k,jspc,idir)*hispec(jspc)
              enddo
              !
            case default
              !
              sum2=sum(dispec(:,1)*(gradyi(:)-spc(i,j,k,:)*sum1))
              !
              do jspc=1,num_species
                !Corretion diffusion velocity for continuity
                corrdiff=sum2*spc(i,j,k,jspc)
                !species diffusive flux
                yflux(i,j,k,jspc,idir)=dispec(jspc,1)*(gradyi(jspc)-(spc(i,j,k,jspc)*sum1)) &
                                        -corrdiff
                !energy flux due to species diffusion
                qflux(i,j,k,idir)=qflux(i,j,k,idir)+yflux(i,j,k,jspc,idir)*hispec(jspc)
              enddo
              !
          end select
          !
        enddo 
#else
        if(nondimen) then 
          dfu(1:num_species)=miu/schmidt(1:num_species)
        else
          dfu(1:num_species)=schmidt(1:num_species)
        endif
        !
        do idir=1,3
          yflux(i,j,k,:,idir)=dfu(:)*dspc(i,j,k,:,idir)
        enddo
#endif
        !
      endif !num_species>0 
      !
    enddo !k
    enddo !j
    enddo !i
    !
    call dataswap(sigma,timerept=ltimrpt)
    !
    call dataswap(qflux,timerept=ltimrpt)
    !
    if(num_species>0) then
      call dataswap(yflux,timerept=ltimrpt)
      ! call yflux_sendrecv(yflux,timerept=ltimrpt)
    endif
    !
    if(trim(turbmode)=='k-omega') then
      call dataswap(dkflux,timerept=ltimrpt)
      !
      call dataswap(doflux,timerept=ltimrpt)
    endif
    !
    ! Calculating along i direction.
    !
    ncolm=5+num_species+num_modequ
    !
    allocate(ff(-hm:im+hm,2:ncolm),df(0:im,2:ncolm))
    do k=0,km
    do j=0,jm
      !
      ff(:,2)=( sigma(:,j,k,1)*dxi(:,j,k,1,1) +                        &
                sigma(:,j,k,2)*dxi(:,j,k,1,2) +                        &
                sigma(:,j,k,3)*dxi(:,j,k,1,3) )*jacob(:,j,k)
      ff(:,3)=( sigma(:,j,k,2)*dxi(:,j,k,1,1) +                        &
                sigma(:,j,k,4)*dxi(:,j,k,1,2) +                        &
                sigma(:,j,k,5)*dxi(:,j,k,1,3) )*jacob(:,j,k)
      ff(:,4)=( sigma(:,j,k,3)*dxi(:,j,k,1,1) +                        &
                sigma(:,j,k,5)*dxi(:,j,k,1,2) +                        &
                sigma(:,j,k,6)*dxi(:,j,k,1,3) )*jacob(:,j,k)
      ff(:,5)=( qflux(:,j,k,1)*dxi(:,j,k,1,1) +                        &
                qflux(:,j,k,2)*dxi(:,j,k,1,2) +                        &
                qflux(:,j,k,3)*dxi(:,j,k,1,3) )*jacob(:,j,k)
      !
      if(num_species>0) then
        do jspc=1,num_species
          ff(:,5+jspc)=( yflux(:,j,k,jspc,1)*dxi(:,j,k,1,1) +               &
                         yflux(:,j,k,jspc,2)*dxi(:,j,k,1,2) +               &
                         yflux(:,j,k,jspc,3)*dxi(:,j,k,1,3) )*jacob(:,j,k)
        enddo
      endif
      !
      if(trim(turbmode)=='k-omega') then
        n=5+num_species
        !
        ff(:,1+n)=( dkflux(:,j,k,1)*dxi(:,j,k,1,1) +                   &
                    dkflux(:,j,k,2)*dxi(:,j,k,1,2) +                   &
                    dkflux(:,j,k,3)*dxi(:,j,k,1,3) )*jacob(:,j,k)
        ff(:,2+n)=( doflux(:,j,k,1)*dxi(:,j,k,1,1) +                   &
                    doflux(:,j,k,2)*dxi(:,j,k,1,2) +                   &
                    doflux(:,j,k,3)*dxi(:,j,k,1,3) )*jacob(:,j,k)
      endif
      !
      !+------------------------------+
      !|    calculate derivative      |
      !+------------------------------+
      do n=2,ncolm
        df(:,n)=fds%central(f=ff(:,n),ntype=npdci,dim=im,dir=1)
      enddo
      !
      !+------------------------------+
      !| end of calculate derivative  |
      !+------------------------------+
      !
      qrhs(is:ie,j,k,2)=qrhs(is:ie,j,k,2)+df(is:ie,2)
      qrhs(is:ie,j,k,3)=qrhs(is:ie,j,k,3)+df(is:ie,3)
      qrhs(is:ie,j,k,4)=qrhs(is:ie,j,k,4)+df(is:ie,4)
      qrhs(is:ie,j,k,5)=qrhs(is:ie,j,k,5)+df(is:ie,5)
      ! freeze energy flux for startup
      ! if(flowtype=='tgvflame'.and.nstep*deltat<5.d-5)qrhs(is:ie,j,k,5)=0.d0
      
      if(num_species>0) then
        do jspc=6,5+num_species
          qrhs(is:ie,j,k,jspc)=qrhs(is:ie,j,k,jspc)+df(is:ie,jspc)
        enddo
      endif
      !
      if(trim(turbmode)=='k-omega') then
        n=5+num_species
        !
        qrhs(is:ie,j,k,1+n)=qrhs(is:ie,j,k,1+n)+df(is:ie,1+n)
        qrhs(is:ie,j,k,2+n)=qrhs(is:ie,j,k,2+n)+df(is:ie,2+n)
      endif
      !
    enddo
    enddo
    !
    deallocate(ff,df)
    !!!!!!!!!!!!!!!!!!!!!!!!!!
    ! End calculating along i
    !!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    if(ndims>=2) then
      ! Calculating along j direction.
      !
      allocate(ff(-hm:jm+hm,2:ncolm),df(0:jm,2:ncolm))
      do k=0,km
      do i=0,im
        !
        ff(:,2)=( sigma(i,:,k,1)*dxi(i,:,k,2,1) +                        &
                  sigma(i,:,k,2)*dxi(i,:,k,2,2) +                        &
                  sigma(i,:,k,3)*dxi(i,:,k,2,3) )*jacob(i,:,k)
        ff(:,3)=( sigma(i,:,k,2)*dxi(i,:,k,2,1) +                        &
                  sigma(i,:,k,4)*dxi(i,:,k,2,2) +                        &
                  sigma(i,:,k,5)*dxi(i,:,k,2,3) )*jacob(i,:,k)
        ff(:,4)=( sigma(i,:,k,3)*dxi(i,:,k,2,1) +                        &
                  sigma(i,:,k,5)*dxi(i,:,k,2,2) +                        &
                  sigma(i,:,k,6)*dxi(i,:,k,2,3) )*jacob(i,:,k)
        ff(:,5)=( qflux(i,:,k,1)*dxi(i,:,k,2,1) +                        &
                  qflux(i,:,k,2)*dxi(i,:,k,2,2) +                        &
                  qflux(i,:,k,3)*dxi(i,:,k,2,3) )*jacob(i,:,k)
        !
        if(num_species>0) then
          do jspc=1,num_species
            ff(:,5+jspc)=( yflux(i,:,k,jspc,1)*dxi(i,:,k,2,1) +          &
                           yflux(i,:,k,jspc,2)*dxi(i,:,k,2,2) +          &
                           yflux(i,:,k,jspc,3)*dxi(i,:,k,2,3) )*jacob(i,:,k)
          enddo
        endif
        !
        if(trim(turbmode)=='k-omega') then
          n=5+num_species
          ff(:,1+n)=( dkflux(i,:,k,1)*dxi(i,:,k,2,1) +                   &
                      dkflux(i,:,k,2)*dxi(i,:,k,2,2) +                   &
                      dkflux(i,:,k,3)*dxi(i,:,k,2,3) )*jacob(i,:,k)
          ff(:,2+n)=( doflux(i,:,k,1)*dxi(i,:,k,2,1) +                   &
                      doflux(i,:,k,2)*dxi(i,:,k,2,2) +                   &
                      doflux(i,:,k,3)*dxi(i,:,k,2,3) )*jacob(i,:,k)
        endif
        !
        !+------------------------------+
        !|    calculate derivative      |
        !+------------------------------+
        do n=2,ncolm
          df(:,n)=fds%central(f=ff(:,n),ntype=npdcj,dim=jm,dir=2)
        enddo
        !+------------------------------+
        !| end of calculate derivative  |
        !+------------------------------+
        !
        qrhs(i,js:je,k,2)=qrhs(i,js:je,k,2)+df(js:je,2)
        qrhs(i,js:je,k,3)=qrhs(i,js:je,k,3)+df(js:je,3)
        qrhs(i,js:je,k,4)=qrhs(i,js:je,k,4)+df(js:je,4)
        qrhs(i,js:je,k,5)=qrhs(i,js:je,k,5)+df(js:je,5)
        ! freeze energy flux for startup
        ! if(flowtype=='tgvflame'.and.nstep*deltat<5.d-5)qrhs(i,js:je,k,5)=0.d0
        
        if(num_species>0) then
          do jspc=6,5+num_species
            qrhs(i,js:je,k,jspc)=qrhs(i,js:je,k,jspc)+df(js:je,jspc)
          enddo
        endif
        !
        if(trim(turbmode)=='k-omega') then
          n=5+num_species
          !
          qrhs(i,js:je,k,1+n)=qrhs(i,js:je,k,1+n)+df(js:je,1+n)
          qrhs(i,js:je,k,2+n)=qrhs(i,js:je,k,2+n)+df(js:je,2+n)
        endif
        !
      enddo
      enddo
      !
      deallocate(ff,df)
      !!!!!!!!!!!!!!!!!!!!!!!!!!
      ! End calculating along j
      !!!!!!!!!!!!!!!!!!!!!!!!!!
      !
    endif
    !
    if(ndims==3) then
      ! Calculating along k direction.
      !
      allocate(ff(-hm:km+hm,2:ncolm),df(0:km,2:ncolm))
      do j=0,jm
      do i=0,im
        !
        ff(:,2)=( sigma(i,j,:,1)*dxi(i,j,:,3,1) +                      &
                  sigma(i,j,:,2)*dxi(i,j,:,3,2) +                      &
                  sigma(i,j,:,3)*dxi(i,j,:,3,3) )*jacob(i,j,:)
        ff(:,3)=( sigma(i,j,:,2)*dxi(i,j,:,3,1) +                      &
                  sigma(i,j,:,4)*dxi(i,j,:,3,2) +                      &
                  sigma(i,j,:,5)*dxi(i,j,:,3,3) )*jacob(i,j,:)
        ff(:,4)=( sigma(i,j,:,3)*dxi(i,j,:,3,1) +                      &
                  sigma(i,j,:,5)*dxi(i,j,:,3,2) +                      &
                  sigma(i,j,:,6)*dxi(i,j,:,3,3) )*jacob(i,j,:)
        ff(:,5)=( qflux(i,j,:,1)*dxi(i,j,:,3,1) +                      &
                  qflux(i,j,:,2)*dxi(i,j,:,3,2) +                      &
                  qflux(i,j,:,3)*dxi(i,j,:,3,3) )*jacob(i,j,:)
        !
        if(num_species>0) then
          do jspc=1,num_species
            ff(:,5+jspc)=( yflux(i,j,:,jspc,1)*dxi(i,j,:,3,1) +        &
                           yflux(i,j,:,jspc,2)*dxi(i,j,:,3,2) +        &
                           yflux(i,j,:,jspc,3)*dxi(i,j,:,3,3) )*jacob(i,j,:)
          enddo
        endif
        !
        if(trim(turbmode)=='k-omega') then
          n=5+num_species
          !
          ff(:,1+n)=( dkflux(i,j,:,1)*dxi(i,j,:,3,1) +                 &
                      dkflux(i,j,:,2)*dxi(i,j,:,3,2) +                 &
                      dkflux(i,j,:,3)*dxi(i,j,:,3,3) )*jacob(i,j,:)
          ff(:,2+n)=( doflux(i,j,:,1)*dxi(i,j,:,3,1) +                 &
                      doflux(i,j,:,2)*dxi(i,j,:,3,2) +                 &
                      doflux(i,j,:,3)*dxi(i,j,:,3,3) )*jacob(i,j,:)
        endif
        !
        !+------------------------------+
        !|    calculate derivative      |
        !+------------------------------+
        do n=2,ncolm
          df(:,n)=fds%central(f=ff(:,n),ntype=npdck,dim=km,dir=3)
        enddo
        !+------------------------------+
        !| end of calculate derivative  |
        !+------------------------------+
        !
        qrhs(i,j,ks:ke,2)=qrhs(i,j,ks:ke,2)+df(ks:ke,2)
        qrhs(i,j,ks:ke,3)=qrhs(i,j,ks:ke,3)+df(ks:ke,3)
        qrhs(i,j,ks:ke,4)=qrhs(i,j,ks:ke,4)+df(ks:ke,4)
        qrhs(i,j,ks:ke,5)=qrhs(i,j,ks:ke,5)+df(ks:ke,5)
        ! freeze energy flux for startup
        ! if(flowtype=='tgvflame'.and.nstep*deltat<5.d-5)qrhs(i,j,ks:ke,5)=0.d0
        !
        if(num_species>0) then
          do jspc=6,5+num_species
            qrhs(i,j,ks:ke,jspc)=qrhs(i,j,ks:ke,jspc)+df(ks:ke,jspc)
          enddo
        endif
        !
        if(trim(turbmode)=='k-omega') then
          n=5+num_species
          !
          qrhs(i,j,ks:ke,1+n)=qrhs(i,j,ks:ke,1+n)+df(ks:ke,1+n)
          qrhs(i,j,ks:ke,2+n)=qrhs(i,j,ks:ke,2+n)+df(ks:ke,2+n)
        endif
        !
      enddo
      enddo
      !
      deallocate(ff,df)
      !!!!!!!!!!!!!!!!!!!!!!!!!!
      ! End calculating along j
      !!!!!!!!!!!!!!!!!!!!!!!!!!
    endif
    !
    ! if(allocated(sigma))  deallocate(sigma)
    ! if(allocated(qflux))  deallocate(qflux)
    ! if(allocated(dkflux)) deallocate(dkflux)
    ! if(allocated(doflux)) deallocate(doflux)
    ! if(allocated(yflux))  deallocate(yflux)
    if(allocated(dispec)) deallocate(dispec)
    if(allocated(dfu))    deallocate(dfu)
    !
    if(present(timerept) .and. timerept) then
      !
      subtime=subtime+ptime()-time_beg
      !
      if(lio .and. lreport .and. ltimrpt) call timereporter(routine='diffrsdcal6', &
                                             timecost=subtime,      &
                                              message='diffusion term')
    endif
    !
    return
    !
  end subroutine diffrsdcal6
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! End of the subroutine diffrsdcal6.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  subroutine check_LR55_unit(Lvec,Rvec,ii,jj,kk)
    !
    real(8),intent(in) :: Lvec(5,5),Rvec(5,5)
    integer,intent(in) :: ii,jj,kk
    !
    real(8) :: matrix(5,5)
    real(8) :: epslion
    integer :: i,j
    logical :: normal
    !
    epslion=1.d-8
    !
    normal=.true.
    !
    matrix=MatMul(Lvec,Rvec)
    !
    do j=1,5
    do i=1,5
     if( i==j ) then
       if(abs(matrix(i,j)-1.d0)<epslion) then
        continue
       else
         ! print*,' !! WARNING of UNIT MARTIX'
         ! print*, ii,jj,kk,matrix(i,j),Lvec(i,j),Rvec(i,j)
         normal=.false.
       endif
     else
       if(abs(matrix(i,j))<epslion) then
        continue
       else
         ! print*,' !! WARNING of UNIT MARTIX'
         ! print*, ii,jj,kk,matrix(i,j),Lvec(i,j),Rvec(i,j)
         normal=.false.
       endif
     endif
    enddo
    enddo
    !
    if(.not. normal) then
      print*,' local characteristic decomposition error at',ii,jj,kk
      print*,' Left vector'
      print*,'---------------------------------------------------------'
      write(*,"(5(F7.4))")Lvec(1,:)
      write(*,"(5(F7.4))")Lvec(2,:)
      write(*,"(5(F7.4))")Lvec(3,:)
      write(*,"(5(F7.4))")Lvec(4,:)
      write(*,"(5(F7.4))")Lvec(5,:)
      print*,' righ vector'
      print*,'---------------------------------------------------------'
      write(*,"(5(F7.4))")Rvec(1,:)
      write(*,"(5(F7.4))")Rvec(2,:)
      write(*,"(5(F7.4))")Rvec(3,:)
      write(*,"(5(F7.4))")Rvec(4,:)
      write(*,"(5(F7.4))")Rvec(5,:)
      !
      stop
    endif
    ! if(normal) then
    !   print*,'---------------------------------------------------------'
    !   write(*,"(5(F7.4))")matrix(1,:)
    !   write(*,"(5(F7.4))")matrix(2,:)
    !   write(*,"(5(F7.4))")matrix(3,:)
    !   write(*,"(5(F7.4))")matrix(4,:)
    !   write(*,"(5(F7.4))")matrix(5,:)
    ! endif
    !
  end subroutine check_LR55_unit
  !!
end module solver
!+---------------------------------------------------------------------+
!| The end of the module solver.                                       |
!+---------------------------------------------------------------------+
