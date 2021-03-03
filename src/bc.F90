!+---------------------------------------------------------------------+
!| This module contains subroutines of applying boundary conditions.   |
!| ==============                                                      |
!| CHANGE RECORD                                                       |
!| -------------                                                       |
!| 13-02-2021: Created by J. Fang @ Warrington                         |
!+---------------------------------------------------------------------+
module bc
  !
  use constdef
  use parallel,only: lio,mpistop,mpirank,mpirankname,irk,jrk,krk,      &
                     irkm,jrkm,krkm
  use commvar, only: hm,im,jm,km,uinf,vinf,winf,pinf,roinf,tinf,ndims, &
                     num_species,flowtype,gamma,numq
  use tecio
  !
  implicit none
  !
  real(8) :: pout
  real(8),allocatable :: rho_in(:,:),vel_in(:,:,:),tmp_in(:,:),        &
                         prs_in(:,:),spc_in(:,:,:)
  !
  contains
  !
  ! Table of values for parameter (ibc(ilat))
  ! ilat = 1,6
  !
  !     ____________           
  !    /|         /|           
  !   / |  4     / |           
  !  /  |    5  /  |  j        
  ! /__________/   |  |        
  ! | 1 |______|_2 |  |____ i  
  ! |  / 6     |  /  /         
  ! | /    3   | /  /          
  ! |/         |/  k           
  ! /----------/
  !+-------------------------------------------------------------------+
  !| This subroutine is to allocate inflow array.                      |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 24-02-2021: Created by J. Fang @ Warrington                       |
  !+-------------------------------------------------------------------+
  subroutine alloinflow(ndir)
    !
    ! arguments
    integer,intent(in) :: ndir
    !
    allocate( rho_in(0:jm,0:km),vel_in(0:jm,0:km,1:3),               &
              tmp_in(0:jm,0:km),prs_in(0:jm,0:km),                   &
              spc_in(0:jm,0:km,1:num_species) )
    !
  end subroutine alloinflow
  !+-------------------------------------------------------------------+
  !| The end of the subroutine alloinflow.                             |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to apply bounday conditions to geometrical     |
  !| variables.                                                        |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 14-02-2021: Created by J. Fang @ Warrington                       |
  !+-------------------------------------------------------------------+
  subroutine geombc
    !
    use commvar, only : bctype
    use commarray, only : jacob,dxi
    use commfunc,  only : gradextrp
    !
    ! local data
    integer :: i,j,k,l
    real(8) :: fex(4,10)
    !
    if(jrk==0) then
      !
      j=0
      !
      if(bctype(3)==41) then
        do k=0,km
        do i=0,im
          !
          fex(:,1) =gradextrp(qbou=dxi(i,0,k,1,1),q1st=dxi(i,1,k,1,1))
          fex(:,2) =gradextrp(qbou=dxi(i,0,k,1,2),q1st=dxi(i,1,k,1,2))
          fex(:,3) =gradextrp(qbou=dxi(i,0,k,1,3),q1st=dxi(i,1,k,1,3))
          fex(:,4) =gradextrp(qbou=dxi(i,0,k,2,1),q1st=dxi(i,1,k,2,1))
          fex(:,5) =gradextrp(qbou=dxi(i,0,k,2,2),q1st=dxi(i,1,k,2,2))
          fex(:,6) =gradextrp(qbou=dxi(i,0,k,2,3),q1st=dxi(i,1,k,2,3))
          fex(:,7) =gradextrp(qbou=dxi(i,0,k,3,1),q1st=dxi(i,1,k,3,1))
          fex(:,8) =gradextrp(qbou=dxi(i,0,k,3,2),q1st=dxi(i,1,k,3,2))
          fex(:,9) =gradextrp(qbou=dxi(i,0,k,3,3),q1st=dxi(i,1,k,3,3))
          fex(:,10)=gradextrp(qbou=jacob(i,0,k),  q1st=jacob(i,1,k))
          !
          do l=1,4
            dxi(i,-l,k,1,1)=fex(l,1) 
            dxi(i,-l,k,1,2)=fex(l,2) 
            dxi(i,-l,k,1,3)=fex(l,3) 
            dxi(i,-l,k,2,1)=fex(l,4) 
            dxi(i,-l,k,2,2)=fex(l,5) 
            dxi(i,-l,k,2,3)=fex(l,6) 
            dxi(i,-l,k,3,1)=fex(l,7) 
            dxi(i,-l,k,3,2)=fex(l,8) 
            dxi(i,-l,k,3,3)=fex(l,9) 
            jacob(i,-l,k)  =fex(l,10)
          enddo
          !
        enddo
        enddo
      endif
      !
    endif
    !
    if(jrk==jrkm) then
      !
      j=0
      !
      if(bctype(4)==41) then
        do k=0,km
        do i=0,im
          !
          dxi(i,jm+1:jm+4,k,1,1)=gradextrp(qbou=dxi(i,jm,k,1,1),q1st=dxi(i,jm-1,k,1,1))
          dxi(i,jm+1:jm+4,k,1,2)=gradextrp(qbou=dxi(i,jm,k,1,2),q1st=dxi(i,jm-1,k,1,2))
          dxi(i,jm+1:jm+4,k,1,3)=gradextrp(qbou=dxi(i,jm,k,1,3),q1st=dxi(i,jm-1,k,1,3))
          dxi(i,jm+1:jm+4,k,2,1)=gradextrp(qbou=dxi(i,jm,k,2,1),q1st=dxi(i,jm-1,k,2,1))
          dxi(i,jm+1:jm+4,k,2,2)=gradextrp(qbou=dxi(i,jm,k,2,2),q1st=dxi(i,jm-1,k,2,2))
          dxi(i,jm+1:jm+4,k,2,3)=gradextrp(qbou=dxi(i,jm,k,2,3),q1st=dxi(i,jm-1,k,2,3))
          dxi(i,jm+1:jm+4,k,3,1)=gradextrp(qbou=dxi(i,jm,k,3,1),q1st=dxi(i,jm-1,k,3,1))
          dxi(i,jm+1:jm+4,k,3,2)=gradextrp(qbou=dxi(i,jm,k,3,2),q1st=dxi(i,jm-1,k,3,2))
          dxi(i,jm+1:jm+4,k,3,3)=gradextrp(qbou=dxi(i,jm,k,3,3),q1st=dxi(i,jm-1,k,3,3))
          jacob(1,jm+1:jm+4,k)  =gradextrp(qbou=jacob(i,jm,k),  q1st=jacob(i,jm-1,k))
          !
        enddo
        enddo
      endif
      !
    endif
    !
  end subroutine geombc
  !
  subroutine xyzbc
    !
    use commvar,  only : bctype
    use commarray,only : x
    use commfunc, only : gradextrp
    !
    ! local data
    integer :: i,j,k,l
    real(8) :: fex(4,3)
    !
    if(jrk==0) then
      !
      j=0
      !
      if(bctype(3)==41) then
        do k=0,km
        do i=0,im
          !
          fex(:,1) =gradextrp(qbou=x(i,0,k,1),q1st=x(i,1,k,1))
          fex(:,2) =gradextrp(qbou=x(i,0,k,2),q1st=x(i,1,k,2))
          fex(:,3) =gradextrp(qbou=x(i,0,k,3),q1st=x(i,1,k,3))
          !
          do l=1,4
            x(i,-l,k,1)=fex(l,1) 
            x(i,-l,k,2)=fex(l,2) 
            x(i,-l,k,3)=fex(l,3) 
          enddo
          !
        enddo
        enddo
      endif
      !
    endif
    !
    if(jrk==jrkm) then
      !
      j=0
      !
      if(bctype(4)==41) then
        do k=0,km
        do i=0,im
          !
          x(i,jm+1:jm+4,k,1)=gradextrp(qbou=x(i,jm,k,1),q1st=x(i,jm-1,k,1))
          x(i,jm+1:jm+4,k,2)=gradextrp(qbou=x(i,jm,k,2),q1st=x(i,jm-1,k,2))
          x(i,jm+1:jm+4,k,3)=gradextrp(qbou=x(i,jm,k,3),q1st=x(i,jm-1,k,3))
          !
        enddo
        enddo
      endif
      !
    endif
    !
  end subroutine xyzbc
  !+-------------------------------------------------------------------+
  !| This subroutine is to apply bounday conditions.                   |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 13-02-2021: Created by J. Fang @ Warrington                       |
  !+-------------------------------------------------------------------+
  subroutine boucon
    !
    use commvar, only : bctype,twall
    !
    ! local data
    integer :: n
    !
    do n=1,6
      !
      if(bctype(n)==41) then
        call noslip(n,twall(n))
      endif
      !
      if(bctype(n)==11) then
        call inflow(n)
      endif
      !
      if(bctype(n)==21) then
        call outflow(n)
      endif
      !
      if(bctype(n)==22) then
        call outflow_nscbc(n)
      endif
      !
    enddo
    !
  end subroutine boucon
  !+-------------------------------------------------------------------+
  !| The end of the subroutine boucon.                                 |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to apply inflow bc.                            |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 24-02-2021: Created by J. Fang @ Warrington                       |
  !+-------------------------------------------------------------------+
  subroutine inflow(ndir)
    !
    use commarray, only : prs,vel,tmp,rho,spc,q
    use fludyna,   only : thermal,fvar2q,q2fvar,sos
    use commfunc,  only : extrapolate
    !
    ! arguments
    integer,intent(in) :: ndir
    !
    ! local data
    integer :: i,j,k,l,jspec
    real(8) :: css,csse,ub,pe,roe,ue
    !
    logical,save :: lfirstcal=.true.
    !
    if(ndir==1 .and. irk==0) then
      !
      if(lfirstcal) then
        !
        call alloinflow(ndir)
        !
        if(trim(flowtype)=='jet') then
          call jetinflow
        else
          call freestreaminflow
        endif
        !
        lfirstcal=.false.
        !
      endif
      !
      i=0
      do k=0,km
      do j=0,jm
        !
        css=sos(tmp(i,j,k))
        ub =vel(i,j,k,1)
        !
        if(ub>=css) then
          ! supersonic inlet
          !
          vel(i,j,k,:)=vel_in(j,k,:)
          tmp(i,j,k)  =tmp_in(j,k)
          prs(i,j,k)  =prs_in(j,k)
          rho(i,j,k)  =rho_in(j,k)
          !
          spc(i,j,k,:)=spc_in(j,k,:)
          !
        elseif(ub<css .and. ub>=0.d0) then
          ! subsonic inlet
          ue  =extrapolate(vel(i+1,j,k,1),vel(i+2,j,k,1),dv=0.d0)
          pe  =extrapolate(prs(i+1,j,k),prs(i+2,j,k),dv=0.d0)
          roe =extrapolate(rho(i+1,j,k),rho(i+2,j,k),dv=0.d0)
          csse=extrapolate(sos(tmp(i+1,j,k)),sos(tmp(i+2,j,k)),dv=0.d0)
          !
          vel(i,j,k,1)=0.5d0*(prs_in(j,k)-pe)/(roe*csse)+0.5d0*(vel_in(j,k,1)+ue)
          vel(i,j,k,2)=vel_in(j,k,2)
          vel(i,j,k,3)=vel_in(j,k,3)
          prs(i,j,k)  =0.5d0*(prs_in(j,k)+pe)+0.5d0*roe*csse*(vel_in(j,k,1)-ue)
          rho(i,j,k)  =rho_in(j,k)*(prs(i,j,k)/prs_in(j,k))**(1.d0/gamma)
          !
          tmp(i,j,k)  =thermal(pressure=prs(i,j,k),density=rho(i,j,k))
          !
          spc(i,j,k,:)=spc_in(j,k,:)
          !
          ! print*,vel_in(j,k,1),rho(i,j,k)
          !
        else
          stop ' !! velocity at inflow error !! @ inflow'
        endif
        !
        call fvar2q(      q=  q(i,j,k,:),   density=rho(i,j,k),        &
                   velocity=vel(i,j,k,:),  pressure=prs(i,j,k),        &
                    species=spc(i,j,k,:)                               )
        !
      enddo
      enddo
      !
    endif
    !
    ! call mpistop
    !
  end subroutine inflow
  !+-------------------------------------------------------------------+
  !| The end of the subroutine inflow.                                 |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to apply outflow bc.                           |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 25-02-2021: Created by J. Fang @ Warrington                       |
  !+-------------------------------------------------------------------+
  subroutine outflow(ndir)
    !
    use commarray, only : prs,vel,tmp,rho,spc,q
    use fludyna,   only : thermal,fvar2q,q2fvar,sos
    use commfunc,  only : extrapolate
    !
    ! arguments
    integer,intent(in) :: ndir
    !
    ! local data
    integer :: i,j,k,l,jspec
    real(8) :: css,csse,ub,pe,roe,ue,ve,we,spce(1:num_species)
    !
    logical,save :: lfirstcal=.true.
    !
    if(ndir==2 .and. irk==irkm) then
      !
      i=im
      do k=0,km
      do j=0,jm
        !
        css=sos(tmp(i,j,k))
        ub =vel(i,j,k,1)
        !
        ue  =extrapolate(vel(i-1,j,k,1),vel(i-2,j,k,1),dv=0.d0)
        ve  =extrapolate(vel(i-1,j,k,2),vel(i-2,j,k,2),dv=0.d0)
        we  =extrapolate(vel(i-1,j,k,3),vel(i-2,j,k,3),dv=0.d0)
        pe  =extrapolate(prs(i-1,j,k),  prs(i-2,j,k),dv=0.d0)
        roe =extrapolate(rho(i-1,j,k),  rho(i-2,j,k),dv=0.d0)
        csse=extrapolate(sos(tmp(i-1,j,k)),sos(tmp(i-2,j,k)),dv=0.d0)
        !
        do jspec=1,num_species
          spce(jspec)=extrapolate(spc(i-1,j,k,jspec),                  &
                                  spc(i-2,j,k,jspec),dv=0.d0)
        enddo
        !
        if(ub>=css) then
          ! supersonic inlet
          !
          vel(i,j,k,1)=ue 
          prs(i,j,k)  =pe
          rho(i,j,k)  =roe
          !
        elseif(ub<css .and. ub>=0.d0) then
          ! subsonic inlet
          vel(i,j,k,1)=-0.5d0*(pinf-pe)/(rho(i,j,k)*css)+0.5d0*(uinf+ue)
          prs(i,j,k)  = 0.5d0*(pinf+pe)+0.5d0*rho(i,j,k)*css*(ue-uinf)
          rho(i,j,k)  = roe*(prs(i,j,k)/pe)**(1.d0/gamma)
          !
        else
          stop ' !! velocity at outflow error !! @ outflow'
        endif
        !
        vel(i,j,k,2)=ve 
        vel(i,j,k,3)=we 
        tmp(i,j,k)  =thermal(pressure=prs(i,j,k),density=rho(i,j,k))
        spc(i,j,k,:)=spce(:)
        !
        call fvar2q(      q=  q(i,j,k,:),   density=rho(i,j,k),        &
                   velocity=vel(i,j,k,:),  pressure=prs(i,j,k),        &
                    species=spc(i,j,k,:)                               )
        !
      enddo
      enddo
      !
    endif
    !
  end subroutine outflow
  !+-------------------------------------------------------------------+
  !| The end of the subroutine outflow.                                |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to apply outflow bc using nscbc.               |
  !+-------------------------------------------------------------------+
  !| ref: Jae Wook Kim, AIAA JOURNAL Vol. 38, No. 11, November 2000    |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 25-02-2021: Created by J. Fang @ Warrington                       |
  !+-------------------------------------------------------------------+
  subroutine outflow_nscbc(ndir)
    !
    use commvar,   only : xmin,xmax
    use commarray, only : prs,vel,tmp,rho,spc,q,qrhs,dxi,jacob
    use fludyna,   only : thermal,fvar2q,q2fvar,sos
    use commfunc,  only : deriv
    !
    ! arguments
    integer,intent(in) :: ndir
    !
    ! local data
    integer :: i,j,k,l,jspec,ii,n
    real(8) :: pinv(5,5),pnor(5,5),Pmult(5,5),E(5),F(5),G(5),Rest(5),  &
               jcbi(3),Ecs(0:2,1:5),dEcs(5),LODi1(5),LODi(5)
    real(8) :: uu,css,gmachmax2,kinout
    !
    if(ndir==2 .and. irk==irkm) then
      !
      i=im
      !
      do k=0,km
      do j=0,jm
        !
        pnor=pmatrix(rho(i,j,k),vel(i,j,k,1),vel(i,j,k,2),            &
                     vel(i,j,k,3),tmp(i,j,k),dxi(i,j,k,1,:),inv=.false.)
        pinv=pmatrix(rho(i,j,k),vel(i,j,k,1),vel(i,j,k,2),            &
                     vel(i,j,k,3),tmp(i,j,k),dxi(i,j,k,1,:),inv=.true.)
        !
        ! Pmult=MatMul(Pinv,pnor)
        ! if(jrk==0) then
        !   print*,'---------------------------------------------------------'
        !   write(*,"(5(F7.4))")Pmult(1,1),Pmult(1,2),Pmult(1,3),Pmult(1,4),Pmult(1,5)
        !   write(*,"(5(F7.4))")Pmult(2,1),Pmult(2,2),Pmult(2,3),Pmult(2,4),Pmult(2,5)
        !   write(*,"(5(F7.4))")Pmult(3,1),Pmult(3,2),Pmult(3,3),Pmult(3,4),Pmult(3,5)
        !   write(*,"(5(F7.4))")Pmult(4,1),Pmult(4,2),Pmult(4,3),Pmult(4,4),Pmult(4,5)
        !   write(*,"(5(F7.4))")Pmult(5,1),Pmult(5,2),Pmult(5,3),Pmult(5,4),Pmult(5,5)
        ! end if
        !
        do ii=0,2
          uu=dxi(i-ii,j,k,1,1)*vel(i-ii,j,k,1) +                       &
             dxi(i-ii,j,k,1,2)*vel(i-ii,j,k,2) +                       &
             dxi(i-ii,j,k,1,3)*vel(i-ii,j,k,3)
          !
          Ecs(ii,1)=jacob(i-ii,j,k)*  q(i-ii,j,k,1)*uu
          Ecs(ii,2)=jacob(i-ii,j,k)*( q(i-ii,j,k,2)*uu+dxi(i-ii,j,k,1,1)*prs(i-ii,j,k) )
          Ecs(ii,3)=jacob(i-ii,j,k)*( q(i-ii,j,k,3)*uu+dxi(i-ii,j,k,1,2)*prs(i-ii,j,k) )
          Ecs(ii,4)=jacob(i-ii,j,k)*( q(i-ii,j,k,4)*uu+dxi(i-ii,j,k,1,3)*prs(i-ii,j,k) )
          Ecs(ii,5)=jacob(i-ii,j,k)*( q(i-ii,j,k,5)+prs(i-ii,j,k) )*uu
        enddo
        !
        do n=1,5
          dEcs(n)=-deriv( Ecs(0,n),Ecs(1,n),Ecs(2,n) )
        enddo
        !
        E(1)= q(i,j,k,2)
        E(2)= q(i,j,k,2)*vel(i,j,k,1)+prs(i,j,k)
        E(3)= q(i,j,k,3)*vel(i,j,k,1)
        E(4)= q(i,j,k,4)*vel(i,j,k,1)
        E(5)=(q(i,j,k,5)+prs(i,j,k))*vel(i,j,k,1)
        !
        F(1)= q(i,j,k,3)
        F(2)= q(i,j,k,2)*vel(i,j,k,2)
        F(3)= q(i,j,k,3)*vel(i,j,k,2)+prs(i,j,k)
        F(4)= q(i,j,k,4)*vel(i,j,k,2)
        F(5)=(q(i,j,k,5)+prs(i,j,k))*vel(i,j,k,2)
        !
        G(1)= q(i,j,k,4)
        G(2)= q(i,j,k,2)*vel(i,j,k,3)
        G(3)= q(i,j,k,3)*vel(i,j,k,3)
        G(4)= q(i,j,k,4)*vel(i,j,k,3)+prs(i,j,k)
        G(5)=(q(i,j,k,5)+prs(i,j,k))*vel(i,j,k,3)
        !
        jcbi(1)=-deriv( dxi(i,j,k,1,1)  *jacob(i,j,k),                 &
                        dxi(i-1,j,k,1,1)*jacob(i-1,j,k),               &
                        dxi(i-2,j,k,1,1)*jacob(i-2,j,k) )
        jcbi(2)=-deriv( dxi(i,j,k,1,2)  *jacob(i,j,k),                 &
                        dxi(i-1,j,k,1,2)*jacob(i-1,j,k),               &
                        dxi(i-2,j,k,1,2)*jacob(i-2,j,k) )
        jcbi(3)=-deriv( dxi(i,j,k,1,3)  *jacob(i,j,k),                 &
                        dxi(i-1,j,k,1,3)*jacob(i-1,j,k),               &
                        dxi(i-2,j,k,1,3)*jacob(i-2,j,k) )
        !
        Rest(1)=E(1)*Jcbi(1)+F(1)*Jcbi(2)+G(1)*Jcbi(3)
        Rest(2)=E(2)*Jcbi(1)+F(2)*Jcbi(2)+G(2)*Jcbi(3)
        Rest(3)=E(3)*Jcbi(1)+F(3)*Jcbi(2)+G(3)*Jcbi(3)
        Rest(4)=E(4)*Jcbi(1)+F(4)*Jcbi(2)+G(4)*Jcbi(3)
        Rest(5)=E(5)*Jcbi(1)+F(5)*Jcbi(2)+G(5)*Jcbi(3)
        !
        LODi1(1)=dEcs(1)-Rest(1)
        LODi1(2)=dEcs(2)-Rest(2)
        LODi1(3)=dEcs(3)-Rest(3)
        LODi1(4)=dEcs(4)-Rest(4)
        LODi1(5)=dEcs(5)-Rest(5)
        !
        LODi=MatMul(Pinv,LODi1)/jacob(i,j,k)
        !
        css=sos(tmp(i,j,k))
        !
        kinout=0.25d0*(1.d0-gmachmax2)*css/(xmax-xmin)
        LODi(5)=kinout*(prs(i,j,k)-pinf)/rho(i,j,k)/css
      enddo
      enddo
      !
    endif
    !
  end subroutine outflow_nscbc
  !
  function pmatrix(rho,u,v,w,t,ddi,inv)
    !
    use commvar, only : gamma
    use fludyna, only : sos
    !
    ! arguments
    real(8) :: pmatrix(5,5)
    logical,intent(in) :: inv
    real(8),intent(in) :: rho,u,v,w,t,ddi(3)
    !
    ! local data
    real(8) :: ke,ma2,css,cs2,b0(3),lxi(3),vlxi(3),cplus(3),cminu(3),  &
               b(3)
    real(8) :: var1,gamm1,rhi,h,cvlxi
    !
    gamm1=gamma-1.d0
    rhi  =1.d0/rho
    !
    ke=(u**2+v**2+w**2)
    css=sos(t)
    cs2=css*css
    ma2=ke/cs2
    !
    var1=1.d0/sqrt(ddi(1)*ddi(1)+ddi(2)*ddi(2)+ddi(3)*ddi(3))
    lxi(1)=ddi(1)*var1
    lxi(2)=ddi(2)*var1
    lxi(3)=ddi(3)*var1
    !
    vlxi(1)=lxi(3)*v-lxi(2)*w
    vlxi(2)=lxi(1)*w-lxi(3)*u
    vlxi(3)=lxi(2)*u-lxi(1)*v
    !
    pmatrix=0.d0
    !
    if(inv) then
      !
      b0(1)=(1.d0-0.5d0*gamm1*ma2)*lxi(1)-vlxi(1)*rhi
      b0(2)=(1.d0-0.5d0*gamm1*ma2)*lxi(2)-vlxi(2)*rhi
      b0(3)=(1.d0-0.5d0*gamm1*ma2)*lxi(3)-vlxi(3)*rhi
      !
      cplus(1)=(lxi(1)-gamm1/css*u)*rhi
      cplus(2)=(lxi(2)-gamm1/css*v)*rhi
      cplus(3)=(lxi(3)-gamm1/css*w)*rhi
      !
      cminu(1)=(-lxi(1)-gamm1/css*u)*rhi
      cminu(2)=(-lxi(2)-gamm1/css*v)*rhi
      cminu(3)=(-lxi(3)-gamm1/css*w)*rhi
      !
      pmatrix(1,1)= b0(1)
      pmatrix(1,2)= gamm1*u/cs2*lxi(1)
      pmatrix(1,3)= gamm1*v/cs2*lxi(1)+lxi(3)*rhi
      pmatrix(1,4)= gamm1*w/cs2*lxi(1)-lxi(2)*rhi
      pmatrix(1,5)=-gamm1/cs2*lxi(1)
      !
      pmatrix(2,1)= b0(2)
      pmatrix(2,2)= gamm1*u/cs2*lxi(2)-lxi(3)*rhi
      pmatrix(2,3)= gamm1*v/cs2*lxi(2)
      pmatrix(2,4)= gamm1*w/cs2*lxi(2)+lxi(1)*rhi
      pmatrix(2,5)=-gamm1/cs2*lxi(2)
      !
      pmatrix(3,1)= b0(3)
      pmatrix(3,2)= gamm1*u/cs2*lxi(3)+lxi(2)*rhi
      pmatrix(3,3)= gamm1*v/cs2*lxi(3)-lxi(1)*rhi
      pmatrix(3,4)= gamm1*w/cs2*lxi(3)
      pmatrix(3,5)=-gamm1/cs2*lxi(3)
      !
      pmatrix(4,1)= css*rhi*(0.5d0*gamm1*ma2-(u*lxi(1)+v*lxi(2)+w*lxi(3))/css)
      pmatrix(4,2)= cplus(1)
      pmatrix(4,3)= cplus(2)
      pmatrix(4,4)= cplus(3)
      pmatrix(4,5)= gamm1/css*rhi
      !
      pmatrix(5,1)= css*rhi*(0.5d0*gamm1*ma2+(u*lxi(1)+v*lxi(2)+w*lxi(3))/css)
      pmatrix(5,2)= cminu(1)
      pmatrix(5,3)= cminu(2)
      pmatrix(5,4)= cminu(3)
      pmatrix(5,5)= gamm1/css*rhi
      !
    else
      !
      pmatrix(1,1)=lxi(1)
      pmatrix(1,2)=lxi(2)
      pmatrix(1,3)=lxi(3)
      pmatrix(1,4)=0.5d0*rho/css
      pmatrix(1,5)=0.5d0*rho/css
      !
      pmatrix(2,1)=lxi(1)*u
      pmatrix(2,2)=lxi(2)*u-lxi(3)*rho
      pmatrix(2,3)=lxi(3)*u+lxi(2)*rho
      pmatrix(2,4)=0.5d0*rho/css*(u+lxi(1)*css)
      pmatrix(2,5)=0.5d0*rho/css*(u-lxi(1)*css)
      !
      pmatrix(3,1)=lxi(1)*v+lxi(3)*rho
      pmatrix(3,2)=lxi(2)*v
      pmatrix(3,3)=lxi(3)*v-lxi(1)*rho
      pmatrix(3,4)=0.5d0*rho/css*(v+lxi(2)*css)
      pmatrix(3,5)=0.5d0*rho/css*(v-lxi(2)*css)
      !
      pmatrix(4,1)=lxi(1)*w-lxi(2)*rho
      pmatrix(4,2)=lxi(2)*w+lxi(1)*rho
      pmatrix(4,3)=lxi(3)*w
      pmatrix(4,4)=0.5d0*rho/css*(w+lxi(3)*css)
      pmatrix(4,5)=0.5d0*rho/css*(w-lxi(3)*css)
      !
      b(1)=0.5d0*ke*lxi(1)+rho*vlxi(1)
      b(2)=0.5d0*ke*lxi(2)+rho*vlxi(2)
      b(3)=0.5d0*ke*lxi(3)+rho*vlxi(3)
      !
      h=0.5d0*ke+css*css/gamm1
      cvlxi=css*(lxi(1)*u+lxi(2)*v+lxi(3)*w)
      !
      pmatrix(5,1)=b(1)
      pmatrix(5,2)=b(2)
      pmatrix(5,3)=b(3)
      pmatrix(5,4)=0.5d0*rho/css*(h+cvlxi)
      pmatrix(5,5)=0.5d0*rho/css*(h-cvlxi)
      !
    endif
    !
  end function pmatrix
  !
  ! subroutine outflow_nscbc(ndir)
  !   !
  !   use commarray, only : prs,vel,tmp,rho,spc,q,qrhs,dxi,jacob
  !   use fludyna,   only : thermal,fvar2q,q2fvar,sos
  !   use commfunc,  only : deriv
  !   !
  !   ! arguments
  !   integer,intent(in) :: ndir
  !   !
  !   ! local data
  !   integer :: i,j,k,l,jspec
  !   real(8) :: css,ri,uu,vv,ww,qq,gm
  !   real(8),allocatable :: qinf(:),vel_inf(:),spc_inf(:)
  !   real(8) :: jac(5,5),jacinv(5,5),el(5,5),er(5,5),ev(5),             &
  !              dwdxi(5),dwdxo(5),dwcdxi(5),dwcdxo(5),dwcdx(5),df(5)
  !   !
  !   if(ndir==2 .and. irk==irkm) then
  !     !
  !     i=im
  !     !
  !     allocate(qinf(numq),vel_inf(3),spc_inf(num_species))
  !     !
  !     vel_inf(1)=uinf
  !     vel_inf(2)=vinf
  !     vel_inf(3)=winf
  !     spc_inf(1)=0.d0
  !     !
  !     call fvar2q(      q=  qinf, density=roinf, velocity=vel_inf,     &
  !                                 pressure=pinf,  species=spc_inf      )
  !     !
  !     do k=0,km
  !     do j=0,jm
  !       !
  !       css=sos(tmp(i,j,k))
  !       !
  !       ri  =  1.d0/rho(i,j,k)
  !       uu  =  vel(i,j,k,1)
  !       vv  =  vel(i,j,k,2)
  !       ww  =  vel(i,j,k,3)
  !       qq  =  0.5d0 * (uu*uu  + vv*vv + ww*ww)
  !       !
  !       !   Jacobian of conservative/primitive transformation          
  !       !   (Eqn. (A.5) of Lodato et al, JCP 2008)
  !       jac(1,1)  =  1.d0
  !       jac(1,2)  =  0.d0
  !       jac(1,3)  =  0.d0
  !       jac(1,4)  =  0.d0
  !       jac(1,5)  =  0.d0
  !       jac(2,1)  =  uu
  !       jac(2,2)  =  rho(i,j,k)
  !       jac(2,3)  =  0.d0
  !       jac(2,4)  =  0.d0
  !       jac(2,5)  =  0.d0
  !       jac(3,1)  =  vv
  !       jac(3,2)  =  0.d0
  !       jac(3,3)  =  rho(i,j,k)
  !       jac(3,4)  =  0.d0
  !       jac(3,5)  =  0.d0
  !       jac(4,1)  =  ww
  !       jac(4,2)  =  0.d0
  !       jac(4,3)  =  0.d0
  !       jac(4,4)  =  rho(i,j,k)
  !       jac(4,5)  =  0.d0
  !       jac(5,1)  =  qq
  !       jac(5,2)  =  q(i,j,k,2)
  !       jac(5,3)  =  q(i,j,k,3)
  !       jac(5,4)  =  q(i,j,k,4)
  !       jac(5,5)  =  1.d0/(gamma-1.d0)
  !       !
  !       ! Jacobian of inverse conservative/primitive transformation 
  !       ! (Eqn. (A.5) of Lodato et al, JCP 2008)
  !       jacinv(1,1) =  1.d0
  !       jacinv(1,2) =  0.d0
  !       jacinv(1,3) =  0.d0
  !       jacinv(1,4) =  0.d0
  !       jacinv(1,5) =  0.d0
  !       jacinv(2,1) = -uu*ri
  !       jacinv(2,2) =  ri
  !       jacinv(2,3) =  0.d0
  !       jacinv(2,4) =  0.d0
  !       jacinv(2,5) =  0.d0
  !       jacinv(3,1) = -vv*ri
  !       jacinv(3,2) =  0.d0
  !       jacinv(3,3) =  ri
  !       jacinv(3,4) =  0.d0
  !       jacinv(3,5) =  0.d0
  !       jacinv(4,1) = -ww*ri
  !       jacinv(4,2) =  0.d0
  !       jacinv(4,3) =  0.d0
  !       jacinv(4,4) =  ri
  !       jacinv(4,5) =  0.d0
  !       jacinv(5,1) =  (gamma-1.d0)*qq
  !       jacinv(5,2) = -(gamma-1.d0)*uu
  !       jacinv(5,3) = -(gamma-1.d0)*vv
  !       jacinv(5,4) = -(gamma-1.d0)*ww
  !       jacinv(5,5) =  (gamma-1.d0)
  !       !
  !       ! left eigenvectors matrix (Eqn.d0 (A.12d0) of 
  !       ! Lodato et al, JCP 2008)
  !       el(1,1) =  0.d0
  !       el(1,2) = -rho(i,j,k)*css
  !       el(1,3) =  0.d0
  !       el(1,4) =  0.d0
  !       el(1,5) =  1.d0
  !       el(2,1) =  css*css
  !       el(2,2) =  0.d0
  !       el(2,3) =  0.d0
  !       el(2,4) =  0.d0
  !       el(2,5) = -1.d0
  !       el(3,1) =  0.d0
  !       el(3,2) =  0.d0
  !       el(3,3) =  1.d0
  !       el(3,4) =  0.d0
  !       el(3,5) =  0.d0
  !       el(4,1) =  0.d0
  !       el(4,2) =  0.d0
  !       el(4,3) =  0.d0
  !       el(4,4) =  1.d0
  !       el(4,5) =  0.d0
  !       el(5,1) =  0.d0
  !       el(5,2) =  rho(i,j,k)*css
  !       el(5,3) =  0.d0
  !       el(5,4) =  0.d0
  !       el(5,5) =  1.d0
  !       !
  !       ! left eigenvectors matrix (Eqn.d0 (A.11d0) of Lodato et al, JCP 2008)
  !       er(1,1) =  0.5d0/css/css
  !       er(2,1) = -0.5d0*ri/css
  !       er(3,1) =  0.d0
  !       er(4,1) =  0.d0
  !       er(5,1) =  0.5d0
  !       er(1,2) =  1.d0/css/css
  !       er(2,2) =  0.d0
  !       er(3,2) =  0.d0
  !       er(4,2) =  0.d0
  !       er(5,2) =  0.d0
  !       er(1,3) =  0.d0
  !       er(2,3) =  0.d0
  !       er(3,3) =  1.d0
  !       er(4,3) =  0.d0
  !       er(5,3) =  0.d0
  !       er(1,4) =  0.d0
  !       er(2,4) =  0.d0
  !       er(3,4) =  0.d0
  !       er(4,4) =  1.d0
  !       er(5,4) =  0.d0
  !       er(1,5) =  0.5d0/css/css
  !       er(2,5) =  0.5d0*ri/css
  !       er(3,5) =  0.d0
  !       er(4,5) =  0.d0
  !       er(5,5) =  0.5d0
  !       !
  !       !  Eigenvalues 
  !       ev(1)    =  uu-css
  !       ev(2)    =  uu
  !       ev(3)    =  uu
  !       ev(4)    =  uu
  !       ev(5)    =  uu+css
  !       !
  !       ! Derivatives of conservative variables
  !       ! Inner derivatives
  !       dwdxi(1) = -deriv(q(i,j,k,1),q(i-1,j,k,1),q(i-2,j,k,1))
  !       dwdxi(2) = -deriv(q(i,j,k,2),q(i-1,j,k,2),q(i-2,j,k,2))
  !       dwdxi(3) = -deriv(q(i,j,k,3),q(i-1,j,k,3),q(i-2,j,k,3))
  !       dwdxi(4) = -deriv(q(i,j,k,4),q(i-1,j,k,4),q(i-2,j,k,4))
  !       dwdxi(5) = -deriv(q(i,j,k,5),q(i-1,j,k,5),q(i-2,j,k,5))
  !       !
  !       ! Outer derivatives
  !       dwdxo(1) = -deriv(qinf(1),q(i,j,k,1))
  !       dwdxo(2) = -deriv(qinf(2),q(i,j,k,2))
  !       dwdxo(3) = -deriv(qinf(3),q(i,j,k,3))
  !       dwdxo(4) = -deriv(qinf(4),q(i,j,k,4))
  !       dwdxo(5) = -deriv(qinf(5),q(i,j,k,5))
  !       !
  !       !   Derivatives of characteristic variables
  !       dwcdxi = 0.d0
  !       dwcdxo = 0.d0
  !       !m=1
  !       !mm=1
  !       dwcdxi(1) = dwcdxi(1) +el(1,1)*jacinv(1,1)*dwdxi(1) !  inner
  !       dwcdxi(1) = dwcdxi(1) +el(1,1)*jacinv(1,2)*dwdxi(2) !  inner
  !       dwcdxi(1) = dwcdxi(1) +el(1,1)*jacinv(1,3)*dwdxi(3) !  inner
  !       dwcdxi(1) = dwcdxi(1) +el(1,1)*jacinv(1,4)*dwdxi(4) !  inner
  !       dwcdxi(1) = dwcdxi(1) +el(1,1)*jacinv(1,5)*dwdxi(5) !  inner
  !       !mm=2
  !       dwcdxi(1) =dwcdxi(1) +el(1,2)*jacinv(2,1)*dwdxi(1) !  inner
  !       dwcdxi(1) =dwcdxi(1) +el(1,2)*jacinv(2,2)*dwdxi(2) !  inner
  !       dwcdxi(1) =dwcdxi(1) +el(1,2)*jacinv(2,3)*dwdxi(3) !  inner
  !       dwcdxi(1) =dwcdxi(1) +el(1,2)*jacinv(2,4)*dwdxi(4) !  inner
  !       dwcdxi(1) =dwcdxi(1) +el(1,2)*jacinv(2,5)*dwdxi(5) !  inner
  !       !mm=3
  !       dwcdxi(1) =dwcdxi(1) +el(1,3)*jacinv(3,1)*dwdxi(1) !  inner
  !       dwcdxi(1) =dwcdxi(1) +el(1,3)*jacinv(3,2)*dwdxi(2) !  inner
  !       dwcdxi(1) =dwcdxi(1) +el(1,3)*jacinv(3,3)*dwdxi(3) !  inner
  !       dwcdxi(1) =dwcdxi(1) +el(1,3)*jacinv(3,4)*dwdxi(4) !  inner
  !       dwcdxi(1) =dwcdxi(1) +el(1,3)*jacinv(3,5)*dwdxi(5) !  inner
  !       !mm=4
  !       dwcdxi(1) =dwcdxi(1) +el(1,4)*jacinv(4,1)*dwdxi(1) !  inner
  !       dwcdxi(1) =dwcdxi(1) +el(1,4)*jacinv(4,2)*dwdxi(2) !  inner
  !       dwcdxi(1) =dwcdxi(1) +el(1,4)*jacinv(4,3)*dwdxi(3) !  inner
  !       dwcdxi(1) =dwcdxi(1) +el(1,4)*jacinv(4,4)*dwdxi(4) !  inner
  !       dwcdxi(1) =dwcdxi(1) +el(1,4)*jacinv(4,5)*dwdxi(5) !  inner
  !       !mm=5
  !       dwcdxi(1) =dwcdxi(1) +el(1,5)*jacinv(5,1)*dwdxi(1) !  inner
  !       dwcdxi(1) =dwcdxi(1) +el(1,5)*jacinv(5,2)*dwdxi(2) !  inner
  !       dwcdxi(1) =dwcdxi(1) +el(1,5)*jacinv(5,3)*dwdxi(3) !  inner
  !       dwcdxi(1) =dwcdxi(1) +el(1,5)*jacinv(5,4)*dwdxi(4) !  inner
  !       dwcdxi(1) =dwcdxi(1) +el(1,5)*jacinv(5,5)*dwdxi(5) !  inner
  !       !
  !       dwcdxo(1) =dwcdxo(1) +el(1,5)*jacinv(5,1)*dwdxo(1) !  outer
  !       dwcdxo(1) =dwcdxo(1) +el(1,5)*jacinv(5,2)*dwdxo(2) !  outer
  !       dwcdxo(1) =dwcdxo(1) +el(1,5)*jacinv(5,3)*dwdxo(3) !  outer
  !       dwcdxo(1) =dwcdxo(1) +el(1,5)*jacinv(5,4)*dwdxo(4) !  outer
  !       dwcdxo(1) =dwcdxo(1) +el(1,5)*jacinv(5,5)*dwdxo(5) !  outer
  !       !
  !       !m=2
  !       !
  !       !mm=1
  !       dwcdxi(2) =dwcdxi(2) +el(2,1)*jacinv(1,1)*dwdxi(1) !  inner
  !       dwcdxi(2) =dwcdxi(2) +el(2,1)*jacinv(1,2)*dwdxi(2) !  inner
  !       dwcdxi(2) =dwcdxi(2) +el(2,1)*jacinv(1,3)*dwdxi(3) !  inner
  !       dwcdxi(2) =dwcdxi(2) +el(2,1)*jacinv(1,4)*dwdxi(4) !  inner
  !       dwcdxi(2) =dwcdxi(2) +el(2,1)*jacinv(1,5)*dwdxi(5) !  inner
  !       !
  !       dwcdxi(2) =dwcdxi(2) +el(2,2)*jacinv(2,1)*dwdxi(1) !  inner
  !       dwcdxi(2) =dwcdxi(2) +el(2,2)*jacinv(2,2)*dwdxi(2) !  inner
  !       dwcdxi(2) =dwcdxi(2) +el(2,2)*jacinv(2,3)*dwdxi(3) !  inner
  !       dwcdxi(2) =dwcdxi(2) +el(2,2)*jacinv(2,4)*dwdxi(4) !  inner
  !       dwcdxi(2) =dwcdxi(2) +el(2,2)*jacinv(2,5)*dwdxi(5) !  inner
  !       !
  !       dwcdxi(2) =dwcdxi(2) +el(2,3)*jacinv(3,1)*dwdxi(1) !  inner
  !       dwcdxi(2) =dwcdxi(2) +el(2,3)*jacinv(3,2)*dwdxi(2) !  inner
  !       dwcdxi(2) =dwcdxi(2) +el(2,3)*jacinv(3,3)*dwdxi(3) !  inner
  !       dwcdxi(2) =dwcdxi(2) +el(2,3)*jacinv(3,4)*dwdxi(4) !  inner
  !       dwcdxi(2) =dwcdxi(2) +el(2,3)*jacinv(3,5)*dwdxi(5) !  inner
  !       !
  !       dwcdxi(2) =dwcdxi(2) +el(2,4)*jacinv(4,1)*dwdxi(1) !  inner
  !       dwcdxi(2) =dwcdxi(2) +el(2,4)*jacinv(4,2)*dwdxi(2) !  inner
  !       dwcdxi(2) =dwcdxi(2) +el(2,4)*jacinv(4,3)*dwdxi(3) !  inner
  !       dwcdxi(2) =dwcdxi(2) +el(2,4)*jacinv(4,4)*dwdxi(4) !  inner
  !       dwcdxi(2) =dwcdxi(2) +el(2,4)*jacinv(4,5)*dwdxi(5) !  inner
  !       !
  !       dwcdxi(2) =dwcdxi(2) +el(2,5)*jacinv(5,1)*dwdxi(1) !  inner
  !       dwcdxi(2) =dwcdxi(2) +el(2,5)*jacinv(5,2)*dwdxi(2) !  inner
  !       dwcdxi(2) =dwcdxi(2) +el(2,5)*jacinv(5,3)*dwdxi(3) !  inner
  !       dwcdxi(2) =dwcdxi(2) +el(2,5)*jacinv(5,4)*dwdxi(4) !  inner
  !       dwcdxi(2) =dwcdxi(2) +el(2,5)*jacinv(5,5)*dwdxi(5) !  inner
  !       !
  !       dwcdxo(2) =dwcdxo(2) +el(2,5)*jacinv(5,1)*dwdxo(1) !  outer
  !       dwcdxo(2) =dwcdxo(2) +el(2,5)*jacinv(5,2)*dwdxo(2) !  outer
  !       dwcdxo(2) =dwcdxo(2) +el(2,5)*jacinv(5,3)*dwdxo(3) !  outer
  !       dwcdxo(2) =dwcdxo(2) +el(2,5)*jacinv(5,4)*dwdxo(4) !  outer
  !       dwcdxo(2) =dwcdxo(2) +el(2,5)*jacinv(5,5)*dwdxo(5) !  outer
  !       !
  !       !m=3
  !       !
  !       !mm=1
  !       dwcdxi(3) =dwcdxi(3) +el(3,1)*jacinv(1,1)*dwdxi(1) !  inner
  !       dwcdxi(3) =dwcdxi(3) +el(3,1)*jacinv(1,2)*dwdxi(2) !  inner
  !       dwcdxi(3) =dwcdxi(3) +el(3,1)*jacinv(1,3)*dwdxi(3) !  inner
  !       dwcdxi(3) =dwcdxi(3) +el(3,1)*jacinv(1,4)*dwdxi(4) !  inner
  !       dwcdxi(3) =dwcdxi(3) +el(3,1)*jacinv(1,5)*dwdxi(5) !  inner
  !       !
  !       dwcdxi(3) =dwcdxi(3) +el(3,2)*jacinv(2,1)*dwdxi(1) !  inner
  !       dwcdxi(3) =dwcdxi(3) +el(3,2)*jacinv(2,2)*dwdxi(2) !  inner
  !       dwcdxi(3) =dwcdxi(3) +el(3,2)*jacinv(2,3)*dwdxi(3) !  inner
  !       dwcdxi(3) =dwcdxi(3) +el(3,2)*jacinv(2,4)*dwdxi(4) !  inner
  !       dwcdxi(3) =dwcdxi(3) +el(3,2)*jacinv(2,5)*dwdxi(5) !  inner
  !       !
  !       dwcdxi(3) =dwcdxi(3) +el(3,3)*jacinv(3,1)*dwdxi(1) !  inner
  !       dwcdxi(3) =dwcdxi(3) +el(3,3)*jacinv(3,2)*dwdxi(2) !  inner
  !       dwcdxi(3) =dwcdxi(3) +el(3,3)*jacinv(3,3)*dwdxi(3) !  inner
  !       dwcdxi(3) =dwcdxi(3) +el(3,3)*jacinv(3,4)*dwdxi(4) !  inner
  !       dwcdxi(3) =dwcdxi(3) +el(3,3)*jacinv(3,5)*dwdxi(5) !  inner
  !       !
  !       dwcdxi(3) =dwcdxi(3) +el(3,4)*jacinv(4,1)*dwdxi(1) !  inner
  !       dwcdxi(3) =dwcdxi(3) +el(3,4)*jacinv(4,2)*dwdxi(2) !  inner
  !       dwcdxi(3) =dwcdxi(3) +el(3,4)*jacinv(4,3)*dwdxi(3) !  inner
  !       dwcdxi(3) =dwcdxi(3) +el(3,4)*jacinv(4,4)*dwdxi(4) !  inner
  !       dwcdxi(3) =dwcdxi(3) +el(3,4)*jacinv(4,5)*dwdxi(5) !  inner
  !       !
  !       dwcdxi(3) =dwcdxi(3) +el(3,5)*jacinv(5,1)*dwdxi(1) !  inner
  !       dwcdxi(3) =dwcdxi(3) +el(3,5)*jacinv(5,2)*dwdxi(2) !  inner
  !       dwcdxi(3) =dwcdxi(3) +el(3,5)*jacinv(5,3)*dwdxi(3) !  inner
  !       dwcdxi(3) =dwcdxi(3) +el(3,5)*jacinv(5,4)*dwdxi(4) !  inner
  !       dwcdxi(3) =dwcdxi(3) +el(3,5)*jacinv(5,5)*dwdxi(5) !  inner
  !       !
  !       dwcdxo(3) =dwcdxo(3) +el(3,5)*jacinv(5,1)*dwdxo(1) !  outer
  !       dwcdxo(3) =dwcdxo(3) +el(3,5)*jacinv(5,2)*dwdxo(2) !  outer
  !       dwcdxo(3) =dwcdxo(3) +el(3,5)*jacinv(5,3)*dwdxo(3) !  outer
  !       dwcdxo(3) =dwcdxo(3) +el(3,5)*jacinv(5,4)*dwdxo(4) !  outer
  !       dwcdxo(3) =dwcdxo(3) +el(3,5)*jacinv(5,5)*dwdxo(5) !  outer
  !       !
  !       !m=4
  !       !
  !       !mm=1
  !       dwcdxi(4) =dwcdxi(4) +el(4,1)*jacinv(1,1)*dwdxi(1) !  inner
  !       dwcdxi(4) =dwcdxi(4) +el(4,1)*jacinv(1,2)*dwdxi(2) !  inner
  !       dwcdxi(4) =dwcdxi(4) +el(4,1)*jacinv(1,3)*dwdxi(3) !  inner
  !       dwcdxi(4) =dwcdxi(4) +el(4,1)*jacinv(1,4)*dwdxi(4) !  inner
  !       dwcdxi(4) =dwcdxi(4) +el(4,1)*jacinv(1,5)*dwdxi(5) !  inner
  !       !
  !       dwcdxi(4) =dwcdxi(4) +el(4,2)*jacinv(2,1)*dwdxi(1) !  inner
  !       dwcdxi(4) =dwcdxi(4) +el(4,2)*jacinv(2,2)*dwdxi(2) !  inner
  !       dwcdxi(4) =dwcdxi(4) +el(4,2)*jacinv(2,3)*dwdxi(3) !  inner
  !       dwcdxi(4) =dwcdxi(4) +el(4,2)*jacinv(2,4)*dwdxi(4) !  inner
  !       dwcdxi(4) =dwcdxi(4) +el(4,2)*jacinv(2,5)*dwdxi(5) !  inner
  !       !
  !       dwcdxi(4) =dwcdxi(4) +el(4,3)*jacinv(3,1)*dwdxi(1) !  inner
  !       dwcdxi(4) =dwcdxi(4) +el(4,3)*jacinv(3,2)*dwdxi(2) !  inner
  !       dwcdxi(4) =dwcdxi(4) +el(4,3)*jacinv(3,3)*dwdxi(3) !  inner
  !       dwcdxi(4) =dwcdxi(4) +el(4,3)*jacinv(3,4)*dwdxi(4) !  inner
  !       dwcdxi(4) =dwcdxi(4) +el(4,3)*jacinv(3,5)*dwdxi(5) !  inner
  !       !
  !       dwcdxi(4) =dwcdxi(4) +el(4,4)*jacinv(4,1)*dwdxi(1) !  inner
  !       dwcdxi(4) =dwcdxi(4) +el(4,4)*jacinv(4,2)*dwdxi(2) !  inner
  !       dwcdxi(4) =dwcdxi(4) +el(4,4)*jacinv(4,3)*dwdxi(3) !  inner
  !       dwcdxi(4) =dwcdxi(4) +el(4,4)*jacinv(4,4)*dwdxi(4) !  inner
  !       dwcdxi(4) =dwcdxi(4) +el(4,4)*jacinv(4,5)*dwdxi(5) !  inner
  !       !
  !       dwcdxi(4) =dwcdxi(4) +el(4,5)*jacinv(5,1)*dwdxi(1) !  inner
  !       dwcdxi(4) =dwcdxi(4) +el(4,5)*jacinv(5,2)*dwdxi(2) !  inner
  !       dwcdxi(4) =dwcdxi(4) +el(4,5)*jacinv(5,3)*dwdxi(3) !  inner
  !       dwcdxi(4) =dwcdxi(4) +el(4,5)*jacinv(5,4)*dwdxi(4) !  inner
  !       dwcdxi(4) =dwcdxi(4) +el(4,5)*jacinv(5,5)*dwdxi(5) !  inner
  !       !
  !       dwcdxo(4) =dwcdxo(4) +el(4,5)*jacinv(5,1)*dwdxo(1) !  outer
  !       dwcdxo(4) =dwcdxo(4) +el(4,5)*jacinv(5,2)*dwdxo(2) !  outer
  !       dwcdxo(4) =dwcdxo(4) +el(4,5)*jacinv(5,3)*dwdxo(3) !  outer
  !       dwcdxo(4) =dwcdxo(4) +el(4,5)*jacinv(5,4)*dwdxo(4) !  outer
  !       dwcdxo(4) =dwcdxo(4) +el(4,5)*jacinv(5,5)*dwdxo(5) !  outer
  !       !
  !       !m=5
  !       !
  !       !mm=1
  !       dwcdxi(5) =dwcdxi(5) +el(5,1)*jacinv(1,1)*dwdxi(1) !  inner
  !       dwcdxi(5) =dwcdxi(5) +el(5,1)*jacinv(1,2)*dwdxi(2) !  inner
  !       dwcdxi(5) =dwcdxi(5) +el(5,1)*jacinv(1,3)*dwdxi(3) !  inner
  !       dwcdxi(5) =dwcdxi(5) +el(5,1)*jacinv(1,4)*dwdxi(4) !  inner
  !       dwcdxi(5) =dwcdxi(5) +el(5,1)*jacinv(1,5)*dwdxi(5) !  inner
  !       !
  !       dwcdxi(5) =dwcdxi(5) +el(5,2)*jacinv(2,1)*dwdxi(1) !  inner
  !       dwcdxi(5) =dwcdxi(5) +el(5,2)*jacinv(2,2)*dwdxi(2) !  inner
  !       dwcdxi(5) =dwcdxi(5) +el(5,2)*jacinv(2,3)*dwdxi(3) !  inner
  !       dwcdxi(5) =dwcdxi(5) +el(5,2)*jacinv(2,4)*dwdxi(4) !  inner
  !       dwcdxi(5) =dwcdxi(5) +el(5,2)*jacinv(2,5)*dwdxi(5) !  inner
  !       !
  !       dwcdxi(5) =dwcdxi(5) +el(5,3)*jacinv(3,1)*dwdxi(1) !  inner
  !       dwcdxi(5) =dwcdxi(5) +el(5,3)*jacinv(3,2)*dwdxi(2) !  inner
  !       dwcdxi(5) =dwcdxi(5) +el(5,3)*jacinv(3,3)*dwdxi(3) !  inner
  !       dwcdxi(5) =dwcdxi(5) +el(5,3)*jacinv(3,4)*dwdxi(4) !  inner
  !       dwcdxi(5) =dwcdxi(5) +el(5,3)*jacinv(3,5)*dwdxi(5) !  inner
  !       !
  !       dwcdxi(5) =dwcdxi(5) +el(5,4)*jacinv(4,1)*dwdxi(1) !  inner
  !       dwcdxi(5) =dwcdxi(5) +el(5,4)*jacinv(4,2)*dwdxi(2) !  inner
  !       dwcdxi(5) =dwcdxi(5) +el(5,4)*jacinv(4,3)*dwdxi(3) !  inner
  !       dwcdxi(5) =dwcdxi(5) +el(5,4)*jacinv(4,4)*dwdxi(4) !  inner
  !       dwcdxi(5) =dwcdxi(5) +el(5,4)*jacinv(4,5)*dwdxi(5) !  inner
  !       !
  !       dwcdxi(5) =dwcdxi(5) +el(5,5)*jacinv(5,1)*dwdxi(1) !  inner
  !       dwcdxi(5) =dwcdxi(5) +el(5,5)*jacinv(5,2)*dwdxi(2) !  inner
  !       dwcdxi(5) =dwcdxi(5) +el(5,5)*jacinv(5,3)*dwdxi(3) !  inner
  !       dwcdxi(5) =dwcdxi(5) +el(5,5)*jacinv(5,4)*dwdxi(4) !  inner
  !       dwcdxi(5) =dwcdxi(5) +el(5,5)*jacinv(5,5)*dwdxi(5) !  inner
  !       !
  !       dwcdxo(5) =dwcdxo(5) +el(5,5)*jacinv(5,1)*dwdxo(1) !  outer
  !       dwcdxo(5) =dwcdxo(5) +el(5,5)*jacinv(5,2)*dwdxo(2) !  outer
  !       dwcdxo(5) =dwcdxo(5) +el(5,5)*jacinv(5,3)*dwdxo(3) !  outer
  !       dwcdxo(5) =dwcdxo(5) +el(5,5)*jacinv(5,4)*dwdxo(4) !  outer
  !       dwcdxo(5) =dwcdxo(5) +el(5,5)*jacinv(5,5)*dwdxo(5) !  outer
  !       !
  !       !   Enforce LODI relations
  !       !m=1
  !       if (ev(1)>0.d0) then
  !         ! Waves pointing out of the domain
  !         dwcdx(1) = dwcdxi(1)
  !       else
  !         ! Waves entering the domain
  !         dwcdx(1) = dwcdxo(1) ! n.r with or without relaxation
  !       endif
  !       !m=2
  !       if (ev(2)>0.d0) then
  !         ! Waves pointing out of the domain
  !         dwcdx(2) = dwcdxi(2)
  !       else
  !         ! Waves entering the domain
  !         dwcdx(2) = dwcdxo(2) ! n.r. with or without relaxation
  !       endif
  !       !m=3
  !       if (ev(3)>0.d0) then
  !         ! Waves pointing out of the domain
  !         dwcdx(3) = dwcdxi(3)
  !       else
  !         ! Waves entering the domain
  !         dwcdx(3) = dwcdxo(3) ! n.r with or without relaxation
  !       endif
  !       !m=4
  !       if (ev(4)>0.d0) then
  !         ! Waves pointing out of the domain
  !         dwcdx(4) = dwcdxi(4)
  !       else
  !         ! Waves entering the domain
  !         dwcdx(4) = dwcdxo(4) ! n.r with or without relaxation
  !       endif
  !       !m=5
  !       if (ev(5)>0.d0) then
  !         ! Waves pointing out of the domain
  !         dwcdx(5) = dwcdxi(5)
  !       else
  !         ! Waves entering the domain
  !         dwcdx(5) = dwcdxo(5) ! n.r with or without relaxation
  !       endif
  !       !
  !       !   Amplitude of characteristic waves
  !       dwcdx = dwcdx * ev
  !       !
  !       !  Return to conservative variables 
  !       df=0.d0
  !       !
  !       !mm=1
  !       df(1) = df(1) + jac(1,1)*er(1,1)*dwcdx(1)
  !       df(1) = df(1) + jac(1,1)*er(1,2)*dwcdx(2)
  !       df(1) = df(1) + jac(1,1)*er(1,3)*dwcdx(3)
  !       df(1) = df(1) + jac(1,1)*er(1,4)*dwcdx(4)
  !       df(1) = df(1) + jac(1,1)*er(1,5)*dwcdx(5)
  !       df(1) = df(1) + jac(1,2)*er(2,1)*dwcdx(1)
  !       df(1) = df(1) + jac(1,2)*er(2,2)*dwcdx(2)
  !       df(1) = df(1) + jac(1,2)*er(2,3)*dwcdx(3)
  !       df(1) = df(1) + jac(1,2)*er(2,4)*dwcdx(4)
  !       df(1) = df(1) + jac(1,2)*er(2,5)*dwcdx(5)
  !       df(1) = df(1) + jac(1,3)*er(3,1)*dwcdx(1)
  !       df(1) = df(1) + jac(1,3)*er(3,2)*dwcdx(2)
  !       df(1) = df(1) + jac(1,3)*er(3,3)*dwcdx(3)
  !       df(1) = df(1) + jac(1,3)*er(3,4)*dwcdx(4)
  !       df(1) = df(1) + jac(1,3)*er(3,5)*dwcdx(5)
  !       df(1) = df(1) + jac(1,4)*er(4,1)*dwcdx(1)
  !       df(1) = df(1) + jac(1,4)*er(4,2)*dwcdx(2)
  !       df(1) = df(1) + jac(1,4)*er(4,3)*dwcdx(3)
  !       df(1) = df(1) + jac(1,4)*er(4,4)*dwcdx(4)
  !       df(1) = df(1) + jac(1,4)*er(4,5)*dwcdx(5)
  !       df(1) = df(1) + jac(1,5)*er(5,1)*dwcdx(1)
  !       df(1) = df(1) + jac(1,5)*er(5,2)*dwcdx(2)
  !       df(1) = df(1) + jac(1,5)*er(5,3)*dwcdx(3)
  !       df(1) = df(1) + jac(1,5)*er(5,4)*dwcdx(4)
  !       df(1) = df(1) + jac(1,5)*er(5,5)*dwcdx(5)
  !       !
  !       df(2) = df(2) + jac(2,1)*er(1,1)*dwcdx(1)
  !       df(2) = df(2) + jac(2,1)*er(1,2)*dwcdx(2)
  !       df(2) = df(2) + jac(2,1)*er(1,3)*dwcdx(3)
  !       df(2) = df(2) + jac(2,1)*er(1,4)*dwcdx(4)
  !       df(2) = df(2) + jac(2,1)*er(1,5)*dwcdx(5)
  !       df(2) = df(2) + jac(2,2)*er(2,1)*dwcdx(1)
  !       df(2) = df(2) + jac(2,2)*er(2,2)*dwcdx(2)
  !       df(2) = df(2) + jac(2,2)*er(2,3)*dwcdx(3)
  !       df(2) = df(2) + jac(2,2)*er(2,4)*dwcdx(4)
  !       df(2) = df(2) + jac(2,2)*er(2,5)*dwcdx(5)
  !       df(2) = df(2) + jac(2,3)*er(3,1)*dwcdx(1)
  !       df(2) = df(2) + jac(2,3)*er(3,2)*dwcdx(2)
  !       df(2) = df(2) + jac(2,3)*er(3,3)*dwcdx(3)
  !       df(2) = df(2) + jac(2,3)*er(3,4)*dwcdx(4)
  !       df(2) = df(2) + jac(2,3)*er(3,5)*dwcdx(5)
  !       df(2) = df(2) + jac(2,4)*er(4,1)*dwcdx(1)
  !       df(2) = df(2) + jac(2,4)*er(4,2)*dwcdx(2)
  !       df(2) = df(2) + jac(2,4)*er(4,3)*dwcdx(3)
  !       df(2) = df(2) + jac(2,4)*er(4,4)*dwcdx(4)
  !       df(2) = df(2) + jac(2,4)*er(4,5)*dwcdx(5)
  !       df(2) = df(2) + jac(2,5)*er(5,1)*dwcdx(1)
  !       df(2) = df(2) + jac(2,5)*er(5,2)*dwcdx(2)
  !       df(2) = df(2) + jac(2,5)*er(5,3)*dwcdx(3)
  !       df(2) = df(2) + jac(2,5)*er(5,4)*dwcdx(4)
  !       df(2) = df(2) + jac(2,5)*er(5,5)*dwcdx(5)

  !       df(3) = df(3) + jac(3,1)*er(1,1)*dwcdx(1)
  !       df(3) = df(3) + jac(3,1)*er(1,2)*dwcdx(2)
  !       df(3) = df(3) + jac(3,1)*er(1,3)*dwcdx(3)
  !       df(3) = df(3) + jac(3,1)*er(1,4)*dwcdx(4)
  !       df(3) = df(3) + jac(3,1)*er(1,5)*dwcdx(5)
  !       df(3) = df(3) + jac(3,2)*er(2,1)*dwcdx(1)
  !       df(3) = df(3) + jac(3,2)*er(2,2)*dwcdx(2)
  !       df(3) = df(3) + jac(3,2)*er(2,3)*dwcdx(3)
  !       df(3) = df(3) + jac(3,2)*er(2,4)*dwcdx(4)
  !       df(3) = df(3) + jac(3,2)*er(2,5)*dwcdx(5)
  !       df(3) = df(3) + jac(3,3)*er(3,1)*dwcdx(1)
  !       df(3) = df(3) + jac(3,3)*er(3,2)*dwcdx(2)
  !       df(3) = df(3) + jac(3,3)*er(3,3)*dwcdx(3)
  !       df(3) = df(3) + jac(3,3)*er(3,4)*dwcdx(4)
  !       df(3) = df(3) + jac(3,3)*er(3,5)*dwcdx(5)
  !       df(3) = df(3) + jac(3,4)*er(4,1)*dwcdx(1)
  !       df(3) = df(3) + jac(3,4)*er(4,2)*dwcdx(2)
  !       df(3) = df(3) + jac(3,4)*er(4,3)*dwcdx(3)
  !       df(3) = df(3) + jac(3,4)*er(4,4)*dwcdx(4)
  !       df(3) = df(3) + jac(3,4)*er(4,5)*dwcdx(5)
  !       df(3) = df(3) + jac(3,5)*er(5,1)*dwcdx(1)
  !       df(3) = df(3) + jac(3,5)*er(5,2)*dwcdx(2)
  !       df(3) = df(3) + jac(3,5)*er(5,3)*dwcdx(3)
  !       df(3) = df(3) + jac(3,5)*er(5,4)*dwcdx(4)
  !       df(3) = df(3) + jac(3,5)*er(5,5)*dwcdx(5)
  !       !
  !       df(4) = df(4) + jac(4,1)*er(1,1)*dwcdx(1)
  !       df(4) = df(4) + jac(4,1)*er(1,2)*dwcdx(2)
  !       df(4) = df(4) + jac(4,1)*er(1,3)*dwcdx(3)
  !       df(4) = df(4) + jac(4,1)*er(1,4)*dwcdx(4)
  !       df(4) = df(4) + jac(4,1)*er(1,5)*dwcdx(5)
  !       df(4) = df(4) + jac(4,2)*er(2,1)*dwcdx(1)
  !       df(4) = df(4) + jac(4,2)*er(2,2)*dwcdx(2)
  !       df(4) = df(4) + jac(4,2)*er(2,3)*dwcdx(3)
  !       df(4) = df(4) + jac(4,2)*er(2,4)*dwcdx(4)
  !       df(4) = df(4) + jac(4,2)*er(2,5)*dwcdx(5)
  !       df(4) = df(4) + jac(4,3)*er(3,1)*dwcdx(1)
  !       df(4) = df(4) + jac(4,3)*er(3,2)*dwcdx(2)
  !       df(4) = df(4) + jac(4,3)*er(3,3)*dwcdx(3)
  !       df(4) = df(4) + jac(4,3)*er(3,4)*dwcdx(4)
  !       df(4) = df(4) + jac(4,3)*er(3,5)*dwcdx(5)
  !       df(4) = df(4) + jac(4,4)*er(4,1)*dwcdx(1)
  !       df(4) = df(4) + jac(4,4)*er(4,2)*dwcdx(2)
  !       df(4) = df(4) + jac(4,4)*er(4,3)*dwcdx(3)
  !       df(4) = df(4) + jac(4,4)*er(4,4)*dwcdx(4)
  !       df(4) = df(4) + jac(4,4)*er(4,5)*dwcdx(5)
  !       df(4) = df(4) + jac(4,5)*er(5,1)*dwcdx(1)
  !       df(4) = df(4) + jac(4,5)*er(5,2)*dwcdx(2)
  !       df(4) = df(4) + jac(4,5)*er(5,3)*dwcdx(3)
  !       df(4) = df(4) + jac(4,5)*er(5,4)*dwcdx(4)
  !       df(4) = df(4) + jac(4,5)*er(5,5)*dwcdx(5)
  !       !
  !       df(5) = df(5) + jac(5,1)*er(1,1)*dwcdx(1)
  !       df(5) = df(5) + jac(5,1)*er(1,2)*dwcdx(2)
  !       df(5) = df(5) + jac(5,1)*er(1,3)*dwcdx(3)
  !       df(5) = df(5) + jac(5,1)*er(1,4)*dwcdx(4)
  !       df(5) = df(5) + jac(5,1)*er(1,5)*dwcdx(5)
  !       df(5) = df(5) + jac(5,2)*er(2,1)*dwcdx(1)
  !       df(5) = df(5) + jac(5,2)*er(2,2)*dwcdx(2)
  !       df(5) = df(5) + jac(5,2)*er(2,3)*dwcdx(3)
  !       df(5) = df(5) + jac(5,2)*er(2,4)*dwcdx(4)
  !       df(5) = df(5) + jac(5,2)*er(2,5)*dwcdx(5)
  !       df(5) = df(5) + jac(5,3)*er(3,1)*dwcdx(1)
  !       df(5) = df(5) + jac(5,3)*er(3,2)*dwcdx(2)
  !       df(5) = df(5) + jac(5,3)*er(3,3)*dwcdx(3)
  !       df(5) = df(5) + jac(5,3)*er(3,4)*dwcdx(4)
  !       df(5) = df(5) + jac(5,3)*er(3,5)*dwcdx(5)
  !       df(5) = df(5) + jac(5,4)*er(4,1)*dwcdx(1)
  !       df(5) = df(5) + jac(5,4)*er(4,2)*dwcdx(2)
  !       df(5) = df(5) + jac(5,4)*er(4,3)*dwcdx(3)
  !       df(5) = df(5) + jac(5,4)*er(4,4)*dwcdx(4)
  !       df(5) = df(5) + jac(5,5)*er(5,1)*dwcdx(1)
  !       df(5) = df(5) + jac(5,5)*er(5,2)*dwcdx(2)
  !       df(5) = df(5) + jac(5,5)*er(5,3)*dwcdx(3)
  !       df(5) = df(5) + jac(5,5)*er(5,4)*dwcdx(4)
  !       df(5) = df(5) + jac(5,5)*er(5,5)*dwcdx(5)
  !       !
  !       qrhs(i,j,k,1) = df(1)*dxi(i,j,k,1,1)/jacob(i,j,k)
  !       qrhs(i,j,k,2) = df(2)*dxi(i,j,k,1,1)/jacob(i,j,k)
  !       qrhs(i,j,k,3) = df(3)*dxi(i,j,k,1,1)/jacob(i,j,k)
  !       qrhs(i,j,k,4) = df(4)*dxi(i,j,k,1,1)/jacob(i,j,k)
  !       qrhs(i,j,k,5) = df(5)*dxi(i,j,k,1,1)/jacob(i,j,k)
  !       !
  !     enddo
  !     enddo
  !     !
  !   endif
  !   !
  ! end subroutine outflow_nscbc
  !+-------------------------------------------------------------------+
  !| The end of the subroutine outflow_nscbc.                          |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to apply nonslip bc.                           |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 13-02-2021: Created by J. Fang @ Warrington                       |
  !+-------------------------------------------------------------------+
  subroutine noslip(ndir,tw)
    !
    use commarray, only : prs,vel,tmp,rho,spc,q
    use fludyna,   only : thermal,fvar2q,q2fvar
    !
    ! arguments
    integer,intent(in) :: ndir
    real(8),intent(in) :: tw
    !
    ! local data
    integer :: i,j,k,l,jspec
    real(8) :: pe
    !
    if(ndir==3) then
      !
      if(jrk==0) then
        !
        j=0
        do k=0,km
        do i=0,im
          !
          pe=num1d3*(4.d0*prs(i,1,k)-prs(i,2,k))
          !
          vel(i,j,k,1)=0.d0
          vel(i,j,k,2)=0.d0
          vel(i,j,k,3)=0.d0
          prs(i,j,k)  =pe
          tmp(i,j,k)  =tw
          !
          do jspec=1,num_species
            spc(i,j,k,jspec)=num1d3*(4.d0*spc(i,1,k,jspec)-spc(i,2,k,jspec))
          enddo
          !
          rho(i,j,k)  =thermal(pressure=prs(i,j,k),temperature=tmp(i,j,k))
          !
          call fvar2q(      q=  q(i,j,k,:),                            &
                      density=rho(i,j,k),                              &
                     velocity=vel(i,j,k,:),                            &
                     pressure=prs(i,j,k),                              &
                      species=spc(i,j,k,:)                             )
          !
          do l=1,hm
            q(i,j-l,k,1)= q(i,j+l,k,1) ! rho   is even
            q(i,j-l,k,2)=-q(i,j+l,k,2) ! rho*u is odd 
            q(i,j-l,k,3)=-q(i,j+l,k,3) ! rho*v is odd 
            q(i,j-l,k,4)=-q(i,j+l,k,4) ! rho*w is odd
            q(i,j-l,k,5)= q(i,j+l,k,5) ! rho*E is even
            !
            do jspec=1,num_species
              q(i,j-l,k,5+jspec)= q(i,j+l,k,5+jspec) ! rho*Yj is even
            enddo
            !
            call q2fvar(q=q(i,j-l,k,:),density=rho(i,j-l,k),           &
                                      velocity=vel(i,j-l,k,:),         &
                                      pressure=prs(i,j-l,k),           &
                                   temperature=tmp(i,j-l,k),           &
                                       species=spc(i,j-l,k,:))
            !
          enddo
          !
        enddo
        enddo
        !
      endif
      !
    elseif(ndir==4) then
      !
      if(jrk==jrkm) then
        !
        j=jm
        do k=0,km
        do i=0,im
          pe=num1d3*(4.d0*prs(i,j-1,k)-prs(i,j-2,k))
          !
          vel(i,j,k,1)=0.d0
          vel(i,j,k,2)=0.d0
          vel(i,j,k,3)=0.d0
          prs(i,j,k)  =pe
          tmp(i,j,k)  =tw
          !
          if(num_species>0) then
            spc(i,j,k,1)=num1d3*(4.d0*spc(i,j-1,k,1)-spc(i,j-2,k,1))
          endif
          !
          rho(i,j,k)  =thermal(pressure=prs(i,j,k),temperature=tmp(i,j,k))
          !
          call fvar2q(      q=  q(i,j,k,:),                            &
                      density=rho(i,j,k),                              &
                     velocity=vel(i,j,k,:),                            &
                     pressure=prs(i,j,k),                              &
                      species=spc(i,j,k,:)                             )
          !
          do l=1,hm
            q(i,j+l,k,1)= q(i,j-l,k,1) ! rho   is even
            q(i,j+l,k,2)=-q(i,j-l,k,2) ! rho*u is odd 
            q(i,j+l,k,3)=-q(i,j-l,k,3) ! rho*v is odd 
            q(i,j+l,k,4)=-q(i,j-l,k,4) ! rho*w is odd
            q(i,j+l,k,5)= q(i,j-l,k,5) ! rho*E is even
            !
            do jspec=1,num_species
              q(i,j+l,k,5+jspec)= q(i,j-l,k,5+jspec) ! rho*Yj is even
            enddo
            !
            call q2fvar(q=q(i,j+l,k,:),density=rho(i,j+l,k),           &
                                      velocity=vel(i,j+l,k,:),         &
                                      pressure=prs(i,j+l,k),           &
                                   temperature=tmp(i,j+l,k),           &
                                       species=spc(i,j+l,k,:))
            !
          enddo
          !
        enddo
        enddo
        !
      endif
      !
    else
      stop ' !! ndir not defined @ noslip'
    endif
    !
  end subroutine noslip
  !+-------------------------------------------------------------------+
  !| The end of the subroutine noslip.                                 |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to obtain the inlet flow for jet.              |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 24-02-2021: Created by J. Fang @ Warrington                       |
  !+-------------------------------------------------------------------+
  subroutine jetinflow
    !
    use fludyna, only : jetvel,thermal
    use commarray, only : x
    !
    ! local data
    integer :: i,j,k
    !
    i=0
    do k=0,km
    do j=0,jm
      !
      rho_in(j,k)  =roinf
      vel_in(j,k,:)=jetvel(x(i,j,k,2))
      tmp_in(j,k)  =tinf
      prs_in(j,k)  =thermal(density=rho_in(j,k),temperature=tmp_in(j,k))
      !
      if(vel_in(j,k,1)>(1.d0+1.d-10)*uinf) then
        spc_in(j,k,1)=1.d0
      else
        spc_in(j,k,1)=0.d0
      endif
      !
      spc_in(j,k,2)=1.d0-spc_in(j,k,1)
      !
    enddo
    enddo
    !
  end subroutine jetinflow
  !+-------------------------------------------------------------------+
  !| The end of the subroutine jetinflow.                              |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to obtain the inlet of freestream flow.        |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 25-02-2021: Created by J. Fang @ Warrington                       |
  !+-------------------------------------------------------------------+
  subroutine freestreaminflow
    !
    use fludyna, only : jetvel,thermal
    use commarray, only : x
    !
    ! local data
    integer :: i,j,k,jspec
    !
    i=0
    do k=0,km
    do j=0,jm
      !
      rho_in(j,k)  =roinf
      vel_in(j,k,1)=uinf
      vel_in(j,k,2)=vinf
      vel_in(j,k,3)=winf
      tmp_in(j,k)  =tinf
      prs_in(j,k)  =thermal(density=rho_in(j,k),temperature=tmp_in(j,k))
      !
      do jspec=1,num_species
        spc_in(j,k,jspec)=0.d0
      enddo
      !
    enddo
    enddo
    !
  end subroutine freestreaminflow
  !+-------------------------------------------------------------------+
  !| The end of the subroutine jetinflow.                              |
  !+-------------------------------------------------------------------+
  !
  !
end module bc
!+---------------------------------------------------------------------+
!| The end of the module bc.                                           |
!+---------------------------------------------------------------------+