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
  use commvar,  only : ndims,ks,ke,hm,hm,lfftk,ctime,nondimen
  !
  implicit none
  !
  real(8),allocatable :: alfa_con(:),alfa_dif(:)
  real(8), allocatable, dimension(:,:) :: cci,ccj,cck,dci,dcj,dck,     &
                                          fci,fcj,fck,uci,ucj,uck,     &
                                          bci,bcj,bck
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
    use commvar, only : numq,num_species,prandtl,gamma,rgas,ia,ja,ka,  &
                        uinf,vinf,winf,roinf,pinf,tinf,const1,const2,  &
                        const3,const4,const5,const6,const7,tempconst,  &
                        tempconst1,reynolds,ref_t,mach,                &
                        num_modequ,turbmode,spcinf,nondimen
    use thermchem, only: spcindex
    use fludyna, only: thermal
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
    if(nondimen) then 
      prandtl=0.72d0
      gamma=1.4d0
      rgas=287.1d0
      !
      const1=1.d0/(gamma*(gamma-1.d0)*mach**2)
      const2=gamma*mach**2
      const3=(gamma-1.d0)/3.d0*prandtl*(mach**2)
      const4=(gamma-1.d0)*mach**2*reynolds*prandtl
      const5=(gamma-1.d0)*mach**2
      const6=1.d0/(gamma-1.d0)
      const7=(gamma-1.d0)*mach**2*Reynolds*prandtl
      !
      uinf=1.d0
      vinf=0.d0
      winf=0.d0
      tinf=1.d0
      roinf=1.d0
      !
      pinf=roinf*tinf/const2
      !
      tempconst=110.3d0/ref_t
      tempconst1=1.d0+tempconst
      !
    else 
      !
      prandtl=0.72d0
      uinf=1.d0
      vinf=0.d0
      winf=0.d0
      pinf=1.01325d5
      tinf=300.d0
      allocate(spcinf(num_species))
      if(num_species==1) then
        spcinf(1)=1.d0
      else
        spcinf(:)=0.d0
        spcinf(spcindex('O2'))=0.233d0
        spcinf(spcindex('N2'))=1.d0-sum(spcinf(:))
      endif
      roinf=thermal(temperature=tinf,pressure=pinf,species=spcinf)
      !
    endif 
    !
  end subroutine refcal
  !+-------------------------------------------------------------------+
  !| The end of the subroutine refcal.                                 |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to initialise solver.                          |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 08-02-2021  | Created by J. Fang @ Warrington                     |
  !+-------------------------------------------------------------------+
  subroutine solvrinit
    !
    use commvar, only : im,jm,km,numq,npdci,npdcj,npdck,               &
                        conschm,difschm,lfilter,alfa_filter,hm,turbmode
    use commfunc,only : coeffcompac,ptds_ini,ptdsfilter_ini,           &
                        ptds_aym_ini,genfilt10coef
    use models,  only : init_komegasst
    !
    ! local data
    integer :: nscheme,i
    !
    ! convectional term
    if(conschm(4:4)=='c') then
      ! a compact scheme is used
      !
      read(conschm(1:3),*) nscheme
      !
      alfa_con=coeffcompac(nscheme)
      !
      if(mod(nscheme/100,2)==0) then
        ! symmetrical central scheme
        call ptds_ini(cci,alfa_con,im,npdci)
        call ptds_ini(ccj,alfa_con,jm,npdcj)
        call ptds_ini(cck,alfa_con,km,npdck)
      else
        ! asymmetrical reconstruction upwind scheme
        call ptds_aym_ini(uci,alfa_con,im,npdci,windir='+')
        call ptds_aym_ini(ucj,alfa_con,jm,npdcj,windir='+')
        call ptds_aym_ini(uck,alfa_con,km,npdck,windir='+')
        !
        call ptds_aym_ini(bci,alfa_con,im,npdci,windir='-')
        call ptds_aym_ini(bcj,alfa_con,jm,npdcj,windir='-')
        call ptds_aym_ini(bck,alfa_con,km,npdck,windir='-')
      endif
      !
    endif
    !
    ! diffusional term
    if(difschm(4:4)=='c') then
      ! a compact scheme is used
      !
      read(difschm(1:3),*) nscheme
      !
      alfa_dif=coeffcompac(nscheme)
      !
      call ptds_ini(dci,alfa_dif,im,npdci)
      call ptds_ini(dcj,alfa_dif,jm,npdcj)
      call ptds_ini(dck,alfa_dif,km,npdck)
      !
    endif
    !
    if(lfilter) then
      !
      call ptdsfilter_ini(fci,alfa_filter,im,npdci)
      call ptdsfilter_ini(fcj,alfa_filter,jm,npdcj)
      call ptdsfilter_ini(fck,alfa_filter,km,npdck)
      !
    endif
    !
    call genfilt10coef(alfa_filter)
    !
    if(trim(turbmode)=='k-omega') then
      call init_komegasst
    endif
    !
    if(lio) print*,' ** numerical solver initilised.'
    !
  end subroutine solvrinit
  !+-------------------------------------------------------------------+
  !| The end of the subroutine solvrinit.                              |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to calulate the rhs of N-S equations.          |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 09-02-2021  | Created by J. Fang @ Warrington                     |
  !+-------------------------------------------------------------------+
  subroutine rhscal(subtime)
    !
    use commarray, only : qrhs,x,q
    use commvar,   only : flowtype,conschm,diffterm,im,jm,             &
                          recon_schem,limmbou,lchardecomp
    use commcal,   only : ShockSolid,ducrossensor
    use tecio
    !
    ! arguments
    real(8),intent(inout),optional :: subtime
    logical,save :: firstcall=.true.
    !
    ! local data
    integer :: j
    integer :: nconv
    !
    real(8) :: time_beg
    !
    if(present(subtime)) time_beg=ptime() 
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
      call gradcal(subtime=ctime(10))
      !
      if(mod(nconv,2)==0) then
        call convrsdcal6(subtime=ctime(9))
      else
        !
        if(conschm(4:4)=='e') then
          !
          if(recon_schem==5 .or. lchardecomp) call ducrossensor(ctime(12))
          !
          call convrsduwd(subtime=ctime(9))
          !
        elseif(conschm(4:4)=='c') then
          !
          if(lchardecomp) call ducrossensor(ctime(12))
          !
          call convrsdcmp(subtime=ctime(9))
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
    if(diffterm) call diffrsdcal6(subtime=ctime(10))
    !
    if(trim(flowtype)=='channel') then 
      call src_chan
    elseif(trim(flowtype)=='rti') then 
      call src_rti
    endif
    !
#ifdef COMB
    call srccomb(subtime=ctime(13))
#endif
    !
    if(present(subtime)) subtime=subtime+ptime()-time_beg
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
  !| This subroutine add a source term to the rsd of the equation to   |
  !| drive channel flow.                                               |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 13-02-2021: Created by J. Fang @ STFC Daresbury Laboratory        |
  !+-------------------------------------------------------------------+
#ifdef COMB
  subroutine srccomb(subtime)
    !
    use commvar,  only : im,jm,km,numq,num_species,odetype,lcomb
    use commarray,only : qrhs,rho,tmp,spc,jacob
    use thermchem,only : chemrate,wirate
    use parallel, only : ptime 
    !
    ! arguments
    real(8),intent(inout),optional :: subtime
    !
    ! local data
    integer :: i,j,k
    real(8) :: time_beg
    !
    if(present(subtime)) time_beg=ptime() 
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
    if(present(subtime)) subtime=subtime+ptime()-time_beg
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
  subroutine convrsduwd(subtime)
    !
    use commvar,  only: im,jm,km,hm,numq,num_species,num_modequ,       &
                        npdci,npdcj,npdck,is,ie,js,je,ks,ke,gamma,     &
                        recon_schem,lchardecomp,conschm,bfacmpld
    use commarray,only: q,vel,rho,prs,tmp,spc,dxi,jacob,qrhs,lsolid,   &
                        lshock,crinod
    use commfunc, only: recons_exp
    use riemann,  only: flux_steger_warming
    use fludyna,  only: sos
    use thermchem,only: aceval,gammarmix
    !
    ! arguments
    real(8),intent(inout),optional :: subtime
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
    !
    logical :: lsh,lso,sson,hdiss
    !
    if(present(subtime)) time_beg=ptime() 
    !
    sson=allocated(lshock)
    !
    eps=0.04d0
    gm2=0.5d0/gamma
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
        if(lchardecomp) then
          !
          call chardecomp(rho(i,j,k),    prs(i,j,k),  q(i,j,k,5),      &
                          vel(i,j,k,:),  dxi(i,j,k,1,:),               &
                          rho(i+1,j,k),  prs(i+1,j,k),q(i+1,j,k,5),    &
                          vel(i+1,j,k,:),dxi(i+1,j,k,1,:),REV,LEV)
          !
          ! if(irk==0) then
          !   print*,'---------------------------------------------------------'
          !   Pmult=MatMul(LEV,REV)
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
        if(lchardecomp) then
          !
          call chardecomp(rho(i,j,k),    prs(i,j,k),  q(i,j,k,5),      &
                          vel(i,j,k,:),  dxi(i,j,k,2,:),               &
                          rho(i,j+1,k),  prs(i,j+1,k),q(i,j+1,k,5),    &
                          vel(i,j+1,k,:),dxi(i,j+1,k,2,:),REV,LEV)
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
        if(lchardecomp) then
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
        if(lchardecomp) then
          !
          call chardecomp(rho(i,j,k),    prs(i,j,k),  q(i,j,k,5),      &
                          vel(i,j,k,:),  dxi(i,j,k,3,:),               &
                          rho(i,j,k+1),  prs(i,j,k+1),q(i,j,k+1,5),    &
                          vel(i,j,k+1,:),dxi(i,j,k+1,3,:),REV,LEV)
          !
          ! Project to characteristic space using local eigenvector
          do m=1,5
            !
            do n=1,8
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
            var1=recons_exp(   f=Flcp(m,:), inode=k,       &
                             dim=km,        ntype=npdck,   &
                             reschem=recon_schem,shock=lsh,solid=lso)
            var2=recons_exp(   f=Flcm(m,:), inode=k,       &
                             dim=km,        ntype=npdck,   &
                             reschem=recon_schem,shock=lsh,solid=lso)
          endif
          !
          Fhc(m)=var1+var2
          !
        enddo
        !
        if(lchardecomp) then
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
    if(present(subtime)) subtime=subtime+ptime()-time_beg
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
  integer function iwind8(i,n,imin,imax,dir)
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
    else
      print*,dir
      stop ' !! ERROR dir@ function iwind8'
    endif
    !
  end function iwind8
  !+-------------------------------------------------------------------+
  !| The end of the function iwind8.                                   |
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
  subroutine convrsdcmp(subtime)
    !
    use commvar,  only: im,jm,km,hm,numq,                              &
                        npdci,npdcj,npdck,is,ie,js,je,ks,ke,gamma,     &
                        recon_schem,lchardecomp,conschm
    use commarray,only: q,vel,rho,prs,tmp,spc,dxi,jacob,qrhs,lshock,crinod
    use commfunc, only: recons,mplimiter
    use riemann,  only: flux_steger_warming
    !
    ! arguments
    real(8),intent(inout),optional :: subtime
    !
    ! local data
    integer :: i,j,k,iss,iee,jss,jee,kss,kee
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
    !
    logical :: lsh,sson,hdiss
    !
    if(present(subtime)) time_beg=ptime() 
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
          call chardecomp(rho(i,j,k),    prs(i,j,k),  q(i,j,k,5),      &
                          vel(i,j,k,:),  dxi(i,j,k,1,:),               &
                          rho(i+1,j,k),  prs(i+1,j,k),q(i+1,j,k,5),    &
                          vel(i+1,j,k,:),dxi(i+1,j,k,1,:),REV,LEV)
          ! Project to characteristic space using local eigenvector
          do m=1,5
            !
            do n=1,5
              ! plus flux
              flcp(m,n)=LEV(m,1)*Fswp(i+n-3,1)+                        &
                        LEV(m,2)*Fswp(i+n-3,2)+                        &
                        LEV(m,3)*Fswp(i+n-3,3)+                        &
                        LEV(m,4)*Fswp(i+n-3,4)+                        &
                        LEV(m,5)*Fswp(i+n-3,5)
              !
              ! minus flux
              flcm(m,n)=LEV(m,1)*Fswm(i+4-n,1)+                        &
                        LEV(m,2)*Fswm(i+4-n,2)+                        &
                        LEV(m,3)*Fswm(i+4-n,3)+                        &
                        LEV(m,4)*Fswm(i+4-n,4)+                        &
                        LEV(m,5)*Fswm(i+4-n,5)
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
            do n=1,5
              flcp(6:numq,n)=Fswp(i+n-3,6:numq)
              flcm(6:numq,n)=Fswm(i+5-n,6:numq)
            end do
            !
            fhcpc(6:numq)=fhcp(i,6:numq)
            fhcmc(6:numq)=fhcm(i,6:numq)
          endif
          ! End of characteristic decomposition.
          !
        else
          ! No characteristic decomposition
          !
          do n=1,5
            ! plus flux
            flcp(1:numq,n)=Fswp(i+n-3,1:numq)
            !
            ! minus flux
            flcm(1:numq,n)=Fswm(i+4-n,1:numq)
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
          call chardecomp(rho(i,j,k),    prs(i,j,k),  q(i,j,k,5),      &
                          vel(i,j,k,:),  dxi(i,j,k,2,:),               &
                          rho(i,j+1,k),  prs(i,j+1,k),q(i,j+1,k,5),    &
                          vel(i,j+1,k,:),dxi(i,j+1,k,2,:),REV,LEV)
          !
          ! Project to characteristic space using local eigenvector
          do m=1,5
            !
            do n=1,5
              ! plus flux
              flcp(m,n)=LEV(m,1)*Fswp(j+n-3,1)+                        &
                        LEV(m,2)*Fswp(j+n-3,2)+                        &
                        LEV(m,3)*Fswp(j+n-3,3)+                        &
                        LEV(m,4)*Fswp(j+n-3,4)+                        &
                        LEV(m,5)*Fswp(j+n-3,5)
              !
              ! minus flux
              flcm(m,n)=LEV(m,1)*Fswm(j+4-n,1)+                        &
                        LEV(m,2)*Fswm(j+4-n,2)+                        &
                        LEV(m,3)*Fswm(j+4-n,3)+                        &
                        LEV(m,4)*Fswm(j+4-n,4)+                        &
                        LEV(m,5)*Fswm(j+4-n,5)
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
              flcp(6:numq,n)=Fswp(j+n-3,6:numq)
              flcm(6:numq,n)=Fswm(j+4-n,6:numq)
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
            flcp(1:numq,n)=Fswp(j+n-3,1:numq)
            !
            ! minus flux
            flcm(1:numq,n)=Fswm(j+4-n,1:numq)
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
          call chardecomp(rho(i,j,k),    prs(i,j,k),  q(i,j,k,5),      &
                          vel(i,j,k,:),  dxi(i,j,k,3,:),               &
                          rho(i,j,k+1),  prs(i,j,k+1),q(i,j,k+1,5),    &
                          vel(i,j,k+1,:),dxi(i,j,k+1,3,:),REV,LEV)
          !
          ! Project to characteristic space using local eigenvector
          do m=1,5
            !
            do n=1,5
              ! plus flux
              flcp(m,n)=LEV(m,1)*Fswp(k+n-3,1)+                        &
                        LEV(m,2)*Fswp(k+n-3,2)+                        &
                        LEV(m,3)*Fswp(k+n-3,3)+                        &
                        LEV(m,4)*Fswp(k+n-3,4)+                        &
                        LEV(m,5)*Fswp(k+n-3,5)
              !
              ! minus flux
              flcm(m,n)=LEV(m,1)*Fswm(k+4-n,1)+                        &
                        LEV(m,2)*Fswm(k+4-n,2)+                        &
                        LEV(m,3)*Fswm(k+4-n,3)+                        &
                        LEV(m,4)*Fswm(k+4-n,4)+                        &
                        LEV(m,5)*Fswm(k+4-n,5)
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
              flcp(6:numq,n)=Fswp(k+n-3,6:numq)
              flcm(6:numq,n)=Fswm(k+4-n,6:numq)
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
            flcp(1:numq,n)=Fswp(k+n-3,1:numq)
            !
            ! minus flux
            flcm(1:numq,n)=Fswm(k+4-n,1:numq)
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
    if(present(subtime)) subtime=subtime+ptime()-time_beg
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
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 22-03-2021: Created by J. Fang @ Warrington.                      |
  !+-------------------------------------------------------------------+
  subroutine chardecomp(ro_l,p_l,E_l,vel_l,ddi_l,                      &
                        ro_r,p_r,E_r,vel_r,ddi_r,REV,LEV)
    !
    use commvar, only: gamma
    !
    real(8),intent(in) :: ro_l,p_l,E_l,vel_l(3),ddi_l(3),              &
                          ro_r,p_r,E_r,vel_r(3),ddi_r(3)
    real(8),intent(out) :: REV(5,5),LEV(5,5)
    !
    ! local data
    real(8) :: WRoe,WRoe1,u1Roe,u2Roe,u3Roe,HL,HR,HRoe,CssRoe,        &
               KRoe,ugp,rcs,var1,var2,var3,var4,rgp,gpd(3),b1,b2
    !
    WRoe=sqrt(ro_l)/(sqrt(ro_l)+sqrt(ro_r))
    WRoe1=1.d0-WRoe
    u1Roe=WRoe*vel_l(1)+WRoe1*vel_r(1)
    u2Roe=WRoe*vel_l(2)+WRoe1*vel_r(2)
    u3Roe=WRoe*vel_l(3)+WRoe1*vel_r(3)
    !
    HL=(E_l+p_l)/ro_l
    HR=(E_r+p_r)/ro_r
    HRoe=WRoe*HL+WRoe1*HR
    KRoe=0.5d0*(u1Roe*u1Roe+u2Roe*u2Roe+u3Roe*u3Roe)
    CssRoe=sqrt((gamma-1.d0)*(HRoe-KRoe))
    rcs=1.d0/CssRoe
    !
    var1=0.5d0*(ddi_l(1)+ddi_r(1))
    var2=0.5d0*(ddi_l(2)+ddi_r(2))
    var3=0.5d0*(ddi_l(3)+ddi_r(3))
    if(abs(var1)<1d-30) var1=1d-30
    var4=sqrt(var1*var1+var2*var2+var3*var3)
    gpd(1)=var1/var4
    gpd(2)=var2/var4
    gpd(3)=var3/var4
    ugp=u1Roe*gpd(1)+u2Roe*gpd(2)+u3Roe*gpd(3)
    rgp=1.d0/gpd(1)
    !
    ! Calculating Left and Right eigenvectors
    REV(1,1)=1.d0
    REV(1,2)=0.d0
    REV(1,3)=0.d0
    REV(1,4)=1.d0
    REV(1,5)=1.d0
    !                          
    REV(2,1)=u1Roe-CssRoe*gpd(1)
    REV(2,2)=-gpd(2)
    REV(2,3)=-gpd(3)
    REV(2,4)=u1Roe
    REV(2,5)=u1Roe+CssRoe*gpd(1)
    !
    REV(3,1)=u2Roe-CssRoe*gpd(2)
    REV(3,2)=gpd(1)
    REV(3,3)=0.d0
    REV(3,4)=u2Roe
    REV(3,5)=u2Roe+CssRoe*gpd(2)
    !
    REV(4,1)=u3Roe-CssRoe*gpd(3)
    REV(4,2)=0.d0
    REV(4,3)=gpd(1)
    REV(4,4)=u3Roe
    REV(4,5)=u3Roe+CssRoe*gpd(3)
    !
    REV(5,1)=HRoe-ugp*CssRoe
    REV(5,2)=u2Roe*gpd(1)-u1Roe*gpd(2)
    REV(5,3)=u3Roe*gpd(1)-u1Roe*gpd(3)
    REV(5,4)=KRoe
    REV(5,5)=HRoe+ugp*CssRoe
    !
    b1=(gamma-1.d0)/(CssRoe*CssRoe)
    b2=b1*KRoe
    !
    LEV(1,1)= 0.5d0*(b2+ugp*rcs)
    LEV(1,2)=-0.5d0*(b1*u1Roe+gpd(1)*rcs)
    LEV(1,3)=-0.5d0*(b1*u2Roe+gpd(2)*rcs)
    LEV(1,4)=-0.5d0*(b1*u3Roe+gpd(3)*rcs)
    LEV(1,5)= 0.5d0*b1
    !
    LEV(2,1)=u1Roe*gpd(2)-u2Roe*(1.d0-gpd(2)**2)*rgp+u3Roe*gpd(2)*gpd(3)*rgp
    LEV(2,2)=-gpd(2)
    LEV(2,3)=(1.d0-gpd(2)**2)*rgp
    LEV(2,4)=-gpd(2)*gpd(3)*rgp
    LEV(2,5)=0.d0
    !
    LEV(3,1)=u1Roe*gpd(3)+u2Roe*gpd(2)*gpd(3)*rgp-u3Roe*(1.d0-gpd(3)**2)*rgp
    LEV(3,2)=-gpd(3)
    LEV(3,3)=-gpd(2)*gpd(3)*rgp
    LEV(3,4)=(1.d0-gpd(3)**2)*rgp
    LEV(3,5)=0.d0
    !
    LEV(4,1)=1.d0-b2
    LEV(4,2)=b1*u1Roe
    LEV(4,3)=b1*u2Roe
    LEV(4,4)=b1*u3Roe
    LEV(4,5)=-b1
    
    LEV(5,1)= 0.5d0*(b2-ugp*rcs)           
    LEV(5,2)=-0.5d0*(b1*u1Roe-gpd(1)*rcs)  
    LEV(5,3)=-0.5d0*(b1*u2Roe-gpd(2)*rcs)  
    LEV(5,4)=-0.5d0*(b1*u3Roe-gpd(3)*rcs)  
    LEV(5,5)= 0.5d0*b1 
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
  subroutine convrsdcal6(subtime)
    !
    use commvar,  only: im,jm,km,hm,numq,num_species,num_modequ,       &
                        conschm,npdci,npdcj,npdck,is,ie,js,je,ks,ke
    use commarray,only: q,vel,rho,prs,tmp,spc,dxi,jacob,qrhs
    use commfunc, only: ddfc
    !
    ! arguments
    real(8),intent(inout),optional :: subtime
    !
    ! local data
    integer :: i,j,k,jspc,jmod,n
    real(8),allocatable :: fcs(:,:),dfcs(:,:),uu(:)
    !
    real(8) :: time_beg
    !
    if(present(subtime)) time_beg=ptime() 
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
        dfcs(:,n)=ddfc(fcs(:,n),conschm,npdci,im,alfa_con,cci)
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
        dfcs(:,n)=ddfc(fcs(:,n),conschm,npdcj,jm,alfa_con,ccj)
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
          dfcs(:,n)=ddfc(fcs(:,n),conschm,npdck,km,alfa_con,cck,lfft=lfftk)
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
    if(present(subtime)) subtime=subtime+ptime()-time_beg
    !
    return
    !
  end subroutine convrsdcal6
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! End of the subroutine ConvRsdCal6.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !+-------------------------------------------------------------------+
  !| This subroutine is to calculate gradients of flow variables.      |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 08-10-2021  | Moved from  diffrsdcal6 by J. Fang @ Warrington     |
  !+-------------------------------------------------------------------+
  subroutine gradcal(subtime)
    !
    use commvar,   only : im,jm,km,npdci,npdcj,npdck,difschm,ndims,    &
                          num_species,num_modequ,is,ie,js,je,ks,ke,    &
                          turbmode
    use commarray, only : vel,tmp,spc,dvel,dtmp,dspc,dxi,omg,tke,dtke, &
                          domg
    use commfunc,  only : ddfc
    !
    ! arguments
    real(8),intent(inout),optional :: subtime
    !
    ! local data
    integer :: i,j,k,n,ncolm
    real(8),allocatable :: df(:,:),ff(:,:)
    !
    real(8) :: time_beg
    !
    if(present(subtime)) time_beg=ptime() 
    !
    dvel=0.d0
    dtmp=0.d0
    dspc=0.d0
    !
    if(trim(turbmode)=='k-omega') then
      dtke=0.d0
      domg=0.d0
    endif
    !
    ncolm=4+num_species+num_modequ
    !
    ! calculate velocity and temperature gradient
    !
    allocate(ff(-hm:im+hm,ncolm),df(0:im,ncolm))
    !
    do k=0,km
    do j=0,jm
      !
      ff(:,1)=vel(:,j,k,1)
      ff(:,2)=vel(:,j,k,2)
      ff(:,3)=vel(:,j,k,3)
      ff(:,4)=tmp(:,j,k)
      !
      if(num_species>0) then
        do n=1,num_species
          ff(:,4+n)=spc(:,j,k,n)
        enddo
      endif
      !
      if(trim(turbmode)=='k-omega') then
        n=4+num_species
        !
        ff(:,n+1)=tke(:,j,k)
        ff(:,n+2)=omg(:,j,k)
      endif
      !
      do n=1,ncolm
        df(:,n)=ddfc(ff(:,n),difschm,npdci,im,alfa_dif,dci)
      enddo
      !
      dvel(:,j,k,1,1)=dvel(:,j,k,1,1)+df(:,1)*dxi(0:im,j,k,1,1)
      dvel(:,j,k,1,2)=dvel(:,j,k,1,2)+df(:,1)*dxi(0:im,j,k,1,2)
      dvel(:,j,k,1,3)=dvel(:,j,k,1,3)+df(:,1)*dxi(0:im,j,k,1,3)
      !
      dvel(:,j,k,2,1)=dvel(:,j,k,2,1)+df(:,2)*dxi(0:im,j,k,1,1)
      dvel(:,j,k,2,2)=dvel(:,j,k,2,2)+df(:,2)*dxi(0:im,j,k,1,2)
      dvel(:,j,k,2,3)=dvel(:,j,k,2,3)+df(:,2)*dxi(0:im,j,k,1,3)
      !
      dvel(:,j,k,3,1)=dvel(:,j,k,3,1)+df(:,3)*dxi(0:im,j,k,1,1)
      dvel(:,j,k,3,2)=dvel(:,j,k,3,2)+df(:,3)*dxi(0:im,j,k,1,2)
      dvel(:,j,k,3,3)=dvel(:,j,k,3,3)+df(:,3)*dxi(0:im,j,k,1,3)
      !
      dtmp(:,j,k,1)=dtmp(:,j,k,1)+df(:,4)*dxi(0:im,j,k,1,1)
      dtmp(:,j,k,2)=dtmp(:,j,k,2)+df(:,4)*dxi(0:im,j,k,1,2)
      dtmp(:,j,k,3)=dtmp(:,j,k,3)+df(:,4)*dxi(0:im,j,k,1,3)
      !
      if(num_species>0) then
        do n=1,num_species
          dspc(:,j,k,n,1)=dspc(:,j,k,n,1)+df(:,4+n)*dxi(0:im,j,k,1,1)
          dspc(:,j,k,n,2)=dspc(:,j,k,n,2)+df(:,4+n)*dxi(0:im,j,k,1,2)
          dspc(:,j,k,n,3)=dspc(:,j,k,n,3)+df(:,4+n)*dxi(0:im,j,k,1,3)
        enddo
      endif
      !
      if(trim(turbmode)=='k-omega') then
        n=4+num_species
        !
        dtke(:,j,k,1)=dtke(:,j,k,1)+df(:,1+n)*dxi(0:im,j,k,1,1)
        dtke(:,j,k,2)=dtke(:,j,k,2)+df(:,1+n)*dxi(0:im,j,k,1,2)
        dtke(:,j,k,3)=dtke(:,j,k,3)+df(:,1+n)*dxi(0:im,j,k,1,3)
        !
        domg(:,j,k,1)=domg(:,j,k,1)+df(:,2+n)*dxi(0:im,j,k,1,1)
        domg(:,j,k,2)=domg(:,j,k,2)+df(:,2+n)*dxi(0:im,j,k,1,2)
        domg(:,j,k,3)=domg(:,j,k,3)+df(:,2+n)*dxi(0:im,j,k,1,3)
      endif
      !
    enddo
    enddo
    !
    deallocate(ff,df)
    !
    allocate(ff(-hm:jm+hm,ncolm),df(0:jm,ncolm))
    do k=0,km
    do i=0,im
      !
      ff(:,1)=vel(i,:,k,1)
      ff(:,2)=vel(i,:,k,2)
      ff(:,3)=vel(i,:,k,3)
      ff(:,4)=tmp(i,:,k)
      !
      if(num_species>0) then
        do n=1,num_species
          ff(:,4+n)=spc(i,:,k,n)
        enddo
      endif
      !
      if(trim(turbmode)=='k-omega') then
        n=4+num_species
        !
        ff(:,n+1)=tke(i,:,k)
        ff(:,n+2)=omg(i,:,k)
      endif
      !
      do n=1,ncolm
        df(:,n)=ddfc(ff(:,n),difschm,npdcj,jm,alfa_dif,dcj)
      enddo
      !
      dvel(i,:,k,1,1)=dvel(i,:,k,1,1)+df(:,1)*dxi(i,0:jm,k,2,1)
      dvel(i,:,k,1,2)=dvel(i,:,k,1,2)+df(:,1)*dxi(i,0:jm,k,2,2)
      dvel(i,:,k,1,3)=dvel(i,:,k,1,3)+df(:,1)*dxi(i,0:jm,k,2,3)
      !
      dvel(i,:,k,2,1)=dvel(i,:,k,2,1)+df(:,2)*dxi(i,0:jm,k,2,1)
      dvel(i,:,k,2,2)=dvel(i,:,k,2,2)+df(:,2)*dxi(i,0:jm,k,2,2)
      dvel(i,:,k,2,3)=dvel(i,:,k,2,3)+df(:,2)*dxi(i,0:jm,k,2,3)
      !
      dvel(i,:,k,3,1)=dvel(i,:,k,3,1)+df(:,3)*dxi(i,0:jm,k,2,1)
      dvel(i,:,k,3,2)=dvel(i,:,k,3,2)+df(:,3)*dxi(i,0:jm,k,2,2)
      dvel(i,:,k,3,3)=dvel(i,:,k,3,3)+df(:,3)*dxi(i,0:jm,k,2,3)
      !
      dtmp(i,:,k,1)=dtmp(i,:,k,1)+df(:,4)*dxi(i,0:jm,k,2,1)
      dtmp(i,:,k,2)=dtmp(i,:,k,2)+df(:,4)*dxi(i,0:jm,k,2,2)
      dtmp(i,:,k,3)=dtmp(i,:,k,3)+df(:,4)*dxi(i,0:jm,k,2,3)
      !
      if(num_species>0) then
        do n=1,num_species
          dspc(i,:,k,n,1)=dspc(i,:,k,n,1)+df(:,4+n)*dxi(i,0:jm,k,2,1)
          dspc(i,:,k,n,2)=dspc(i,:,k,n,2)+df(:,4+n)*dxi(i,0:jm,k,2,2)
          dspc(i,:,k,n,3)=dspc(i,:,k,n,3)+df(:,4+n)*dxi(i,0:jm,k,2,3)
        enddo
      endif
      !
      if(trim(turbmode)=='k-omega') then
        n=4+num_species
        !
        dtke(i,:,k,1)=dtke(i,:,k,1)+df(:,1+n)*dxi(i,0:jm,k,2,1)
        dtke(i,:,k,2)=dtke(i,:,k,2)+df(:,1+n)*dxi(i,0:jm,k,2,2)
        dtke(i,:,k,3)=dtke(i,:,k,3)+df(:,1+n)*dxi(i,0:jm,k,2,3)
        !
        domg(i,:,k,1)=domg(i,:,k,1)+df(:,2+n)*dxi(i,0:jm,k,2,1)
        domg(i,:,k,2)=domg(i,:,k,2)+df(:,2+n)*dxi(i,0:jm,k,2,2)
        domg(i,:,k,3)=domg(i,:,k,3)+df(:,2+n)*dxi(i,0:jm,k,2,3)
        !
        !
      endif
      !
    enddo
    enddo
    deallocate(ff,df)
    !
    if(ndims==3) then
      allocate(ff(-hm:km+hm,ncolm),df(0:km,ncolm))
      do j=0,jm
      do i=0,im
        !
        ff(:,1)=vel(i,j,:,1)
        ff(:,2)=vel(i,j,:,2)
        ff(:,3)=vel(i,j,:,3)
        ff(:,4)=tmp(i,j,:)
        !
        if(num_species>0) then
          do n=1,num_species
            ff(:,4+n)=spc(i,j,:,n)
          enddo
        endif
        !
        if(trim(turbmode)=='k-omega') then
          n=4+num_species
          !
          ff(:,n+1)=tke(i,j,:)
          ff(:,n+2)=omg(i,j,:)
        endif
        !
        do n=1,ncolm
          df(:,n)=ddfc(ff(:,n),difschm,npdck,km,alfa_dif,dck,lfft=lfftk)
        enddo
        !
        dvel(i,j,:,1,1)=dvel(i,j,:,1,1)+df(:,1)*dxi(i,j,0:km,3,1)
        dvel(i,j,:,1,2)=dvel(i,j,:,1,2)+df(:,1)*dxi(i,j,0:km,3,2)
        dvel(i,j,:,1,3)=dvel(i,j,:,1,3)+df(:,1)*dxi(i,j,0:km,3,3)
        !
        dvel(i,j,:,2,1)=dvel(i,j,:,2,1)+df(:,2)*dxi(i,j,0:km,3,1)
        dvel(i,j,:,2,2)=dvel(i,j,:,2,2)+df(:,2)*dxi(i,j,0:km,3,2)
        dvel(i,j,:,2,3)=dvel(i,j,:,2,3)+df(:,2)*dxi(i,j,0:km,3,3)
        !
        dvel(i,j,:,3,1)=dvel(i,j,:,3,1)+df(:,3)*dxi(i,j,0:km,3,1)
        dvel(i,j,:,3,2)=dvel(i,j,:,3,2)+df(:,3)*dxi(i,j,0:km,3,2)
        dvel(i,j,:,3,3)=dvel(i,j,:,3,3)+df(:,3)*dxi(i,j,0:km,3,3)
        !
        dtmp(i,j,:,1)=dtmp(i,j,:,1)+df(:,4)*dxi(i,j,0:km,3,1)
        dtmp(i,j,:,2)=dtmp(i,j,:,2)+df(:,4)*dxi(i,j,0:km,3,2)
        dtmp(i,j,:,3)=dtmp(i,j,:,3)+df(:,4)*dxi(i,j,0:km,3,3)
        !
        if(num_species>0) then
          do n=1,num_species
            dspc(i,j,:,n,1)=dspc(i,j,:,n,1)+df(:,4+n)*dxi(i,j,0:km,3,1)
            dspc(i,j,:,n,2)=dspc(i,j,:,n,2)+df(:,4+n)*dxi(i,j,0:km,3,2)
            dspc(i,j,:,n,3)=dspc(i,j,:,n,3)+df(:,4+n)*dxi(i,j,0:km,3,3)
          enddo
        endif
        !
        if(trim(turbmode)=='k-omega') then
          n=4+num_species
          !
          dtke(i,j,:,1)=dtke(i,j,:,1)+df(:,1+n)*dxi(i,j,0:km,3,1)
          dtke(i,j,:,2)=dtke(i,j,:,2)+df(:,1+n)*dxi(i,j,0:km,3,2)
          dtke(i,j,:,3)=dtke(i,j,:,3)+df(:,1+n)*dxi(i,j,0:km,3,3)
          !
          domg(i,j,:,1)=domg(i,j,:,1)+df(:,2+n)*dxi(i,j,0:km,3,1)
          domg(i,j,:,2)=domg(i,j,:,2)+df(:,2+n)*dxi(i,j,0:km,3,2)
          domg(i,j,:,3)=domg(i,j,:,3)+df(:,2+n)*dxi(i,j,0:km,3,3)
        endif
        !
      enddo
      enddo
      deallocate(ff,df)
    endif
    !
    if(present(subtime)) subtime=subtime+ptime()-time_beg
    !
    return
    !
  end subroutine gradcal
  !+-------------------------------------------------------------------+
  !| The end of the subroutine gradcal.                                |
  !+-------------------------------------------------------------------+
  !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! This subroutine is used to calculate the diffusion term with 6-order
  ! Compact Central scheme.
  !   sixth-order Compact Central scheme.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Writen by Fang Jian, 2009-06-09.
  ! Add scalar transport equation by Fang Jian, 2022-01-12.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine diffrsdcal6(subtime)
    !
    use commvar,   only : im,jm,km,numq,npdci,npdcj,npdck,difschm,     &
                          conschm,ndims,num_species,num_modequ,        &
                          reynolds,prandtl,const5,is,ie,js,je,ks,ke,   &
                          turbmode,nondimen,schmidt,nstep,deltat,flowtype
    use commarray, only : vel,tmp,spc,dvel,dtmp,dspc,dxi,x,jacob,qrhs, &
                          rho,vor,omg,tke,miut,dtke,domg,res12
    use commfunc,  only : ddfc
    use fludyna,   only : miucal
    use models,    only : komega,src_komega
    use tecio
    use parallel,  only : yflux_sendrecv
#ifdef COMB
    use thermchem, only : tranmod,tranco,enthpy,convertxiyi,wmolar
#endif
    !
    ! arguments
    real(8),intent(inout),optional :: subtime
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
    !
    if(present(subtime)) time_beg=ptime() 
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
      if(nondimen) allocate(dfu(1:num_species))
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
      if(nondimen) then 
        miu=miucal(tmp(i,j,k))/reynolds
      else
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
#endif
        !
      endif 
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
        miu2=2.d0*(miu+miut(i,j,k))
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
        if(nondimen) then 
          hcc=(miu/prandtl)/const5
        else
          hcc=kama
        endif 
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
        if(nondimen) then 
          !
          dfu(1:num_species)=miu/schmidt(1:num_species)
          !
          do idir=1,3
            yflux(i,j,k,:,idir)=dfu(:)*dspc(i,j,k,:,idir)
          enddo
          !
        else
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
#endif
          !
        endif !nondimen
        !
      endif !num_species>0 
      !
    enddo !k
    enddo !j
    enddo !i
    !
    if(present(subtime)) subtime=subtime+ptime()-time_beg
    !
    call dataswap(sigma,subtime=ctime(7))
    !
    call dataswap(qflux,subtime=ctime(7))
    !
    if(num_species>0) then
      call dataswap(yflux,subtime=ctime(7))
      ! call yflux_sendrecv(yflux,subtime=ctime(7))
    endif
    !
    if(trim(turbmode)=='k-omega') then
      call dataswap(dkflux,subtime=ctime(7))
      !
      call dataswap(doflux,subtime=ctime(7))
    endif
    !
    if(present(subtime)) time_beg=ptime() 
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
        df(:,n)=ddfc(ff(:,n),difschm,npdci,im,alfa_dif,dci)
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
        df(:,n)=ddfc(ff(:,n),difschm,npdcj,jm,alfa_dif,dcj)
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
          df(:,n)=ddfc(ff(:,n),difschm,npdck,km,alfa_dif,dck,lfft=lfftk)
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
    if(present(subtime)) subtime=subtime+ptime()-time_beg
    !
    return
    !
  end subroutine diffrsdcal6
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! End of the subroutine diffrsdcal6.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  subroutine filter2e(phi)
    !
    use commvar,  only : is,ie,je,js,je,ks,ke,im,jm,km
    !
    real(8),intent(inout) :: phi(-hm:im+hm,-hm:jm+hm,-hm:km+hm) 
    !
    integer :: i,j,k
    real(8),allocatable :: phtemp(:,:,:)
    !
    call dataswap(phi)
    !
    allocate(phtemp(is:ie,js:je,ks:ke))
    do k=ks,ke
    do i=is,ie
      !
      do j=js,je
        phtemp(i,j,k)=0.01d0*(0.25d0*(phi(i,j-1,k)+phi(i,j+1,k))+0.5d0*phi(i,j,k))  + &
                      0.99d0*phi(i,j,k)
      enddo
      !
    enddo
    enddo
    !
    phi(is:ie,js:je,ks:ke)=phtemp(is:ie,js:je,ks:ke)
    !
    return
    !
  end subroutine filter2e
  !
  subroutine filter4e(phi)
    !
    use commvar,  only : is,ie,je,js,je,ks,ke,im,jm,km
    !
    real(8),intent(inout) :: phi(-hm:im+hm,-hm:jm+hm,-hm:km+hm) 
    !
    integer :: i,j,k
    real(8),allocatable :: phtemp(:,:,:)
    !
    call dataswap(phi)
    !
    allocate(phtemp(is:ie,js:je,ks:ke))
    do k=ks,ke
    do i=is,ie
      !
      do j=js+1,je-1
        phtemp(i,j,k)= 0.0001d0*(0.625d0*phi(i,j,k)               + &
                               0.25d0*(phi(i,j-1,k)+phi(i,j+1,k)) - &
                              0.06250*(phi(i,j-2,k)+phi(i,j+3,k)))+ &
                       0.9999d0*phi(i,j,k)
  
      enddo
      !
    enddo
    enddo
    !
    phi(is:ie,js:je,ks:ke)=phtemp(is:ie,js:je,ks:ke)
    !
    return
    !
  end subroutine filter4e
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! This subroutine is used for spatial filter the conservative variable
  ! for stabilizing the computation.
  ! 2-order Explicit filter is incorporated.
  ! Ref1: Datta V. Gaitonde and Miguel R. Visbal, AIAA JOURNAL Vol.38,
  !      No.11, November 2000. 
  ! Ref2: Xavier Gloerfelt and Philippe Lafon, Computers & Fluids, 2008,
  !       37: 388-401.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine spongefilter
    !
    use commvar,  only : spg_imin,spg_imax,spg_jmin,spg_jmax,          &
                         spg_kmin,spg_kmax,im,jm,km,ia,ja,ka,          &
                         lisponge,ljsponge,lksponge,is,je,js,je,ks,ke, &
                         numq
    use commarray,only: lspg_imin,lspg_imax,lspg_jmin,lspg_jmax,       &
                        lspg_kmin,lspg_kmax,x,q
    use commfunc, only : spafilter10
    !
    real(8),parameter :: dampfac=0.05d0
    !
    integer :: i,j,k,n
    real(8) :: dis,var1
    real(8),allocatable :: qtemp(:,:,:,:)
    !
    if(lisponge) then
      !
      call dataswap(q,direction=1,subtime=ctime(7))
      !
      if(spg_imax>0) then
        !
        allocate(qtemp(im-spg_imax:im-1,js:je,ks:ke,1:numq))
        !
        do k=ks,ke
        do j=js,je
          !
          dis=0.d0
          do i=im-spg_imax,im-1
            !
            dis=dis+ sqrt( (x(i+1,j,k,1)-x(i,j,k,1))**2+               &
                           (x(i+1,j,k,2)-x(i,j,k,2))**2+               &
                           (x(i+1,j,k,3)-x(i,j,k,3))**2                )
            var1=dampfac*(dis/lspg_imax(j,k))**2
            !
            do n=1,numq
              qtemp(i,j,k,n)=(1.d0-var1)*q(i,j,k,n)+        &
                             num1d6*var1*(q(i+1,j,k,n)+     &
                                          q(i-1,j,k,n)+     &
                                          q(i,j+1,k,n)+     &
                                          q(i,j-1,k,n)+     &
                                          q(i,j,k+1,n)+     &
                                          q(i,j,k-1,n)      )
            enddo
            !
          enddo
          !
        enddo
        enddo
        !
        do k=ks,ke
        do j=js,je
          !
          do i=im-spg_imax,im-1
            !
            do n=1,numq
              q(i,j,k,n)=qtemp(i,j,k,n)
            enddo
            !
          enddo
          !
        enddo
        enddo
        !
        deallocate(qtemp)
        !
      endif
      !
    endif
    !
  end subroutine spongefilter
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! End of the subroutine spongefilter.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! This subroutine is used for spatial filter the conservative variable
  ! for stabilizing the computation.
  ! 10-order filter is incorporated.
  ! for boundary filter: the high-order one side filter is used.
  ! the 0-6-6-6-8-10.............-10-8-6-6-6-0. boundary order is
  ! dopted.
  ! Ref: Datta V. Gaitonde and Miguel R. Visbal, AIAA JOURNAL Vol.38,
  !      No.11, November 2000.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Writen by Fang Jian, 2008-11-03.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine filterq(subtime)
    !
    use commvar,  only : im,jm,km,numq,npdci,npdcj,npdck,              &
                         alfa_filter,ndims,is,ie,js,je,ks,ke,turbmode
    use commarray,only : q
    use commfunc, only : spafilter10,spafilter6exp
    !
    ! arguments
    real(8),intent(inout),optional :: subtime
    !
    ! local data
    integer :: i,j,k,n,m
    real(8),allocatable :: phi(:,:),fph(:,:)
    !
    real(8) :: time_beg
    !
    if(present(subtime)) time_beg=ptime() 
    !
    ! filtering in i direction
    call dataswap(q,direction=1,subtime=ctime(7))
    !
    allocate(phi(-hm:im+hm,1:numq),fph(0:im,1:numq))
    !
    do k=0,km
    do j=0,jm
      !
      phi(:,:)=q(:,j,k,:)
      !
      do n=1,numq
        fph(:,n)=spafilter10(phi(:,n),npdci,im,alfa_filter,fci)
        ! fph(:,n)=spafilter6exp(phi(:,n),npdci,im)
      enddo
      !
      q(0:im,j,k,:)=fph(0:im,:)
      !
      ! if(npdci==1) then
      !   q(2:im,j,k,:)=fph(2:im,:)
      ! elseif(npdci==2) then
      !   q(0:im-2,j,k,:)=fph(0:im-2,:)
      ! elseif(npdci==3) then
      !   q(0:im,j,k,:)=fph(0:im,:)
      ! endif
      !
    end do
    end do
    !
    deallocate(phi,fph)
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! end filter in i direction.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! filtering in j direction
    call dataswap(q,direction=2,subtime=ctime(7))
    !
    allocate(phi(-hm:jm+hm,1:numq),fph(0:jm,1:numq))
    !
    do k=0,km
    do i=0,im
      !
      phi(:,:)=q(i,:,k,:)
      !
      do n=1,numq
        fph(:,n)=spafilter10(phi(:,n),npdcj,jm,alfa_filter,fcj)
        ! fph(:,n)=spafilter6exp(phi(:,n),npdcj,jm)
      enddo
      !
      q(i,0:jm,k,:)=fph(0:jm,:)
      !
      ! if(npdcj==1) then
      !   q(i,2:jm,k,:)=fph(2:jm,:)
      ! elseif(npdcj==2) then
      !   q(i,0:jm-2,k,:)=fph(0:jm-2,:)
      ! elseif(npdcj==3) then
      !   q(i,0:jm,k,:)=fph(0:jm,:)
      ! endif
      !
      !
    end do
    end do
    !
    deallocate(phi,fph)
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! end filter in j direction.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    if(ndims==3) then
      !
      call dataswap(q,direction=3,subtime=ctime(7))
      !
      !
      allocate(phi(-hm:km+hm,1:numq),fph(0:km,1:numq))
      !
      ! filtering in k direction
      do j=0,jm
      do i=0,im
        !
        phi(:,:)=q(i,j,:,:)
        !
        do n=1,numq
          fph(:,n)=spafilter10(phi(:,n),npdck,km,alfa_filter,fck,lfft=lfftk)
        enddo
        !
        q(i,j,0:km,:)=fph
        !
      end do
      end do
      !
      deallocate(phi,fph)
      !
    end if
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! end filter in k direction.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !
    ! call filter2e(q(:,:,:,1))
    ! call filter2e(q(:,:,:,2))
    ! call filter2e(q(:,:,:,3))
    ! call filter2e(q(:,:,:,4))
    ! call filter2e(q(:,:,:,5))
    ! call filter2e(q(:,:,:,6))
    if(trim(turbmode)=='k-omega') then
      call filter2e(q(:,:,:,7))
    endif
    !
    if(present(subtime)) subtime=subtime+ptime()-time_beg
    !
    return
    !
  end subroutine filterq
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! End of the subroutine filterq.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!

  !!
  subroutine gradtest
    !
    use commvar,   only : im,jm,km,npdci,npdcj,npdck,conschm,          &
                          alfa_filter,numq,is,ie
    use commarray, only : x,q,dxi
    use commfunc,  only : ddfc,recons,spafilter10,spafilter6exp
    use bc,       only : boucon
    !
    ! local data
    integer :: i,j,k,n
    real(8) :: dx
    real(8),allocatable :: dq(:,:),qhp(:),qhm(:)
    !
    ! do k=0,km
    ! do j=0,jm
    ! do i=0,im
    !   ! q(i,j,k,1)=sin(4.d0*x(1,1,k,3)) !+0.1d0*sin(pi*j+0.5d0*pi)
    !   do n=1,numq
    !     q(i,j,k,n)=sin(pi*x(i,j,k,2))
    !     ! sin(10.d0*(x(1,1,k,3)-0.5*pi))*exp(-4.d0*(x(1,1,k,3)-0.5*pi)**2)+0.1d0*sin(pi*k+0.5d0*pi)
    !   enddo
    !   !
    ! enddo
    ! enddo
    ! enddo
    !
    call dataswap(q,subtime=ctime(7))
    ! !
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
    allocate(dq(0:jm,1:1))
    ! dq(:,1)=ddfc(q(0,:,0,2),conschm,npdcj,jm,alfa_con,ccj)*dxi(0,0:jm,0,2,2)
    !
    ! dq(:,1)=spafilter10(q(1,:,0,3),npdcj,jm,alfa_filter,fcj)
    dq(:,1)=spafilter6exp(q(1,:,0,3),npdcj,jm)
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
    open(18,file='testout/profile'//mpirankname//'.dat')
    do j=0,jm
      write(18,'(3(1X,E15.7E3))')x(1,j,0,2),q(1,j,0,3),dq(j,1)
    enddo
    close(18)
    print*,' << testout/profile',mpirankname,'.dat'
    !
    deallocate(dq)
    !
  end subroutine gradtest
  !
end module solver
!+---------------------------------------------------------------------+
!| The end of the module commarray.                                    |
!+---------------------------------------------------------------------+