!+---------------------------------------------------------------------+
!| This module contains some common calculater.                        |
!+---------------------------------------------------------------------+
!| CHANGE RECORD                                                       |
!| -------------                                                       |
!| 19-03-2021  | Created by J. Fang                                    |
!+---------------------------------------------------------------------+
module commcal
  !
  use parallel, only: mpirank,mpistop,mpirankname
  use utility,  only: timereporter
  !
  implicit none
  !
  contains
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to calculate CFL number and the           |
  !| corresponding time step.                                          |
  !+-------------------------------------------------------------------+
  !| ref: Adams, N. A., Shariff, K. 1996, J COMPUT PHYS. 127,27-51.    |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 21-03-2021: Created by J. Fang @ STFC Daresbury Laboratory        |
  !+-------------------------------------------------------------------+
  subroutine cflcal(deltat)
    !
    use commvar,  only: im,jm,km,nondimen
    use commarray,only: vel,rho,prs,tmp,dxi,jacob,spc
    use fludyna,  only: sos
    use parallel, only: pmax
    !
    ! arguments
    real(8),intent(in) :: deltat
    !
    ! local data
    real(8) :: cfl,ubar,vbar,wbar,css,csi,csj,csk
    integer :: i,j,k
    real(8) :: deltai,deltaj,deltak
    !
    deltai=0.d0
    deltaj=0.d0
    deltak=0.d0
    do k=0,km
    do j=0,jm
    do i=0,im
      ubar=dxi(i,j,k,1,1)*vel(i,j,k,1)+dxi(i,j,k,1,2)*vel(i,j,k,2)+    &
           dxi(i,j,k,1,3)*vel(i,j,k,3)
      vbar=dxi(i,j,k,2,1)*vel(i,j,k,1)+dxi(i,j,k,2,2)*vel(i,j,k,2)+    &
           dxi(i,j,k,2,3)*vel(i,j,k,3)
      wbar=dxi(i,j,k,3,1)*vel(i,j,k,1)+dxi(i,j,k,3,2)*vel(i,j,k,2)+    &
           dxi(i,j,k,3,3)*vel(i,j,k,3)
      !
      css=sos(tmp(i,j,k),spc(i,j,k,:))
      !
      csi=css*sqrt(dxi(i,j,k,1,1)**2+dxi(i,j,k,1,2)**2+dxi(i,j,k,1,3)**2)
      csj=css*sqrt(dxi(i,j,k,2,1)**2+dxi(i,j,k,2,2)**2+dxi(i,j,k,2,3)**2)
      csk=css*sqrt(dxi(i,j,k,3,1)**2+dxi(i,j,k,3,2)**2+dxi(i,j,k,3,3)**2)
      !
      deltai=max(deltai,ubar,ubar-csi,ubar+csi)
      deltaj=max(deltaj,vbar,vbar-csj,vbar+csj)
      deltak=max(deltak,wbar,wbar-csk,wbar+csk)
      !
    enddo
    enddo
    enddo
    !
    deltai=pmax(deltai)
    deltaj=pmax(deltaj)
    deltak=pmax(deltak)
    !
    cfl=deltat*(deltai+deltaj+deltak)
    !
    if(mpirank==0) then
      write(*,"(A38)")'  =========== CFL Condition==========='
      write(*,"(A24,1x,E13.5)")'     current time step: ',deltat
      write(*,"(A24,1x,F13.7)")'           current CFL: ',cfl
      write(*,"(A24,1x,E13.5)")'   time step for CFL=1: ',deltat/cfl
      write(*,"(A38)")'  ===================================='
    end if
    !
  end subroutine cflcal
  !+-------------------------------------------------------------------+
  !| The end of the subroutine cflcal.                                 |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to search monitor points.                 |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 13-07-2021: Created by J. Fang @ STFC Daresbury Laboratory        |
  !+-------------------------------------------------------------------+
  subroutine monitorsearch(ijk,xyz,ijkmon,nmon)
    !
    use commvar,   only : im,jm,km
    use commarray, only : x
    use parallel,  only : ig0,jg0,kg0,irk,jrk,krk
    !
    ! arguments 
    integer,intent(in) :: ijk(:,:)
    real(8),intent(in) :: xyz(:,:)
    integer,allocatable,intent(out) :: ijkmon(:,:)
    integer,intent(out) :: nmon
    !
    ! local data
    integer :: i,j,k,nsize,n
    integer :: is,js,ks,ig,jg,kg
    integer,allocatable :: ijktemp(:,:)
    !
    nsize=size(ijk,1)
    allocate(ijktemp(nsize,4))
    !
    if(irk==0) then
      is=0
    else
      is=1
    endif
    !
    if(jrk==0) then
      js=0
    else
      js=1
    endif
    !
    if(krk==0) then
      ks=0
    else
      ks=1
    endif
    !
    nmon=0
    !
    do k=ks,km
    do j=js,jm
    do i=is,im
      !
      ig=ig0+i
      jg=jg0+j
      kg=kg0+k
      !
      do n=1,nsize
        !
        if(ijk(n,1)==-1) then
          !
          ! xyz(n,1)
          ! xyz(n,2)
          ! xyz(n,3)
          !
        else
          !
          if(ig==ijk(n,1) .and. jg==ijk(n,2) .and. kg==ijk(n,3)) then
            !
            nmon=nmon+1
            !
            ijktemp(nmon,1)=i
            ijktemp(nmon,2)=j
            ijktemp(nmon,3)=k
            ijktemp(nmon,4)=ijk(n,4) 
            !
          endif
          !
        endif
        !
      enddo
      !
    enddo
    enddo
    enddo
    !
    allocate(ijkmon(nmon,4))
    ijkmon(1:nmon,1:4)=ijktemp(1:nmon,1:4)
    !
    deallocate(ijktemp)
    !
  end subroutine monitorsearch
  !+-------------------------------------------------------------------+
  !| The end of the subroutine monitorsearch.                          |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to be a shock senor.                      |
  !+-------------------------------------------------------------------+
  !! Ref1: F. Ducros, J. Comp. Phys. 1999, 152:517-549.                |
  !! Ref2: S. C. Lo, Int. J. Num. Meth. Fluid., 2009.                  |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 08-10-2021: Created by J. Fang @ STFC Daresbury Laboratory        |
  !+-------------------------------------------------------------------+
  subroutine ducrossensor(timerept)
    !
    use commvar,  only : im,jm,km,is,ie,js,je,ks,ke,ia,ja,ka,hm,      &
                         npdci,npdcj,npdck,shkcrt,lreport,ltimrpt
    use commarray,only : ssf,lshock,dvel,prs
    use parallel, only : dataswap,pmin,pmax,psum,lio,ptime
    !
    logical,intent(in),optional :: timerept
    !
    ! local data
    logical,save :: firstcall=.true.
    real(8) :: div2,vort,vortx,vorty,vortz,dpdi,dpdj,dpdk
    real(8) :: ssfmin,ssfmax,ssfavg,norm,ssfmax_local
    integer :: i,j,k,i1,j1,k1,ii,jj,kk,ip1,jp1,kp1,im1,jm1,km1,nshknod,nijka
    !
    real(8) :: time_beg
    real(8),save :: subtime=0.d0
    !
    if(present(timerept) .and. timerept) time_beg=ptime() 
    !
    if(firstcall) then
      !
      allocate(ssf(-hm:im+hm,-hm:jm+hm,-hm:km+hm))
      allocate(lshock(0:im,0:jm,0:km))
      !
      firstcall=.false.
      !
    endif
    !
    ssfmin= 1.d10
    ssfmax=-1.d10
    ssfavg=0.d0
    norm=0.d0
    do k=0,km
    do j=0,jm
    do i=0,im
      !
      div2=(dvel(i,j,k,1,1)+dvel(i,j,k,2,2)+dvel(i,j,k,3,3))**2
      !
      vortx=dvel(i,j,k,3,2)-dvel(i,j,k,2,3)
      vorty=dvel(i,j,k,1,3)-dvel(i,j,k,3,1)
      vortz=dvel(i,j,k,2,1)-dvel(i,j,k,1,2)
      vort=vortx*vortx+vorty*vorty+vortz*vortz
      !
      ip1=i+1; im1=i-1
      jp1=j+1; jm1=j-1
      kp1=k+1; km1=k-1
      !
      if(npdci==1 .and. im1<0)  im1=0
      if(npdcj==1 .and. jm1<0)  jm1=0
      if(npdck==1 .and. km1<0)  km1=0
      if(npdci==2 .and. ip1>im) ip1=im
      if(npdcj==2 .and. jp1>jm) jp1=jm
      if(npdck==2 .and. kp1>km) kp1=km
      !
      dpdi= abs(prs(ip1,j,k)-2.d0*prs(i,j,k)+prs(im1,j,k)) /  &
               (prs(ip1,j,k)+2.d0*prs(i,j,k)+prs(im1,j,k))
      dpdj= abs(prs(i,jp1,k)-2.d0*prs(i,j,k)+prs(i,jm1,k)) /  &
               (prs(i,jp1,k)+2.d0*prs(i,j,k)+prs(i,jm1,k))
      dpdk= abs(prs(i,j,kp1)-2.d0*prs(i,j,k)+prs(i,j,km1)) /  &
               (prs(i,j,kp1)+2.d0*prs(i,j,k)+prs(i,j,km1))
      !
      ssf(i,j,k)=div2/(div2+vort+1.d-30) * max(dpdi,dpdj,dpdk)
      !
      ssfmin=min(ssfmin,ssf(i,j,k))
      ssfmax=max(ssfmax,ssf(i,j,k))
      !
      ssfavg=ssfavg+ssf(i,j,k)
      norm=norm+1.d0
      !
    enddo
    enddo
    enddo
    !
    call dataswap(ssf,timerept=ltimrpt)
    !
    nshknod=0
    !
    do k=0,km
    do j=0,jm
    do i=0,im
      !
      ssfmax_local=0.d0
      do i1=-hm+1,hm
        !
        ii=i+i1
        !
        if(npdci==1 .and. ii<0)  ii=0
        if(npdci==2 .and. ii>im) ii=im
        !
        ssfmax_local=max(ssfmax_local,ssf(ii,j,k))
        !
      enddo
      do j1=-hm+1,hm
        !
        jj=j+j1
        !
        if(npdcj==1 .and. jj<0)  jj=0
        if(npdcj==2 .and. jj>jm) jj=jm
        !
        ssfmax_local=max(ssfmax_local,ssf(i,jj,k))
        !
      enddo
      do k1=-hm+1,hm
        !
        kk=k+k1
        !
        if(npdck==1 .and. kk<0)  kk=0
        if(npdck==2 .and. kk>km) kk=km
        !
        ssfmax_local=max(ssfmax_local,ssf(i,j,kk))
        !
      enddo
      !
      if(ssfmax_local>shkcrt) then
        lshock(i,j,k)=.true.
        nshknod=nshknod+1
      else
        lshock(i,j,k)=.false.
      endif
      !
    enddo
    enddo
    enddo
    !
    if(lreport) then
      !
      nshknod=psum(nshknod)
      !
      ssfmin=pmin(ssfmin)
      ssfmax=pmax(ssfmax)
      ssfavg=psum(ssfavg)/psum(norm)
      !
      if(lio) then
        nijka=(ia+1)*(ja+1)*(ka+1)
        print*,' ------------- shock sensor -------------'
        write(*,"(4x,A,11X,F12.5)")'      max ssf: ',ssfmax
        write(*,"(4x,A,11X,F12.5)")'      min ssf: ',ssfmin
        write(*,"(4x,A,11X,F12.5)")'      avg ssf: ',ssfavg
        write(*,"(4x,2(A,I0),(A,F10.5,A))")'  shock nodes:  ',       &
                 nshknod,'/',nijka,' = ',dble(100*nshknod)/dble(nijka),' %'
        print*,' ----------------------------------------'
      endif
      !
    endif
    !
    if(present(timerept) .and. timerept) then
      !
      subtime=subtime+ptime()-time_beg
      !
      if(lio .and. lreport) call timereporter(routine='ducrossensor', &
                                             timecost=subtime, &
                                              message='shock sensor')
    endif
    !
    return
    !
  end subroutine ducrossensor
  !+-------------------------------------------------------------------+
  !| The end of the subroutine ducrossensor.                           |
  !+-------------------------------------------------------------------+
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! This subroutine is used to be a shock senor to control the nonlinear
  ! shock-capturing schemes with Jameson's approach.
  ! Ref1: F. Ducros, J. Comp. Phys. 1999, 152:517-549.
  ! Ref2: S. C. Lo, Int. J. Num. Meth. Fluid., 2009.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Writen by Fang Jian, 2010-04-20.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine ShockSolid
    !
    use parallel, only : npdci,npdcj,npdck,dataswap
    use commvar,  only : im,jm,km,hm,ltimrpt
    use commarray,only : nodestat,lsolid,x
    use tecio
    !
    ! local data
    integer :: i,j,k,m,i1,j1,k1
    logical,allocatable :: lss(:,:,:)
    real(8),allocatable :: rss(:,:,:)
    !
    ! set the shock area 
    allocate(lss(-hm:im+hm,-hm:jm+hm,-hm:km+hm))
    !
    lss=.false.
    do k=0,km
    do j=0,jm
    do i=0,im
      !
      if(nodestat(i,j,k)>0 .or. nodestat(i,j,k)==-1) then
        lss(i,j,k)=.true.
      endif
      !
    end do
    end do
    end do
    !
    call dataswap(lss,timerept=ltimrpt)
    !
    lsolid=lss
    !
    ! expand the solid area 
    lsolid=.false.
    do k=-hm,km+hm
    do j=-hm,jm+hm
    do i=-hm,im+hm
      !
      if(lss(i,j,k)) then
        !
        do m=-4,4
          !
          i1=i+m
          !
          if(i1>im+hm) i1=im+hm
          if(i1<-hm)   i1=-hm
          !
          lsolid(i1,j,k)=.true.
          !
          j1=j+m
          !
          if(j1>jm+hm) j1=jm+hm
          if(j1<-hm)   j1=-hm
          !
          lsolid(i,j1,k)=.true.
          !
          k1=k+m
          !
          if(k1>km+hm) k1=km+hm
          if(k1<-hm)   k1=-hm
          !
          lsolid(i,j,k1)=.true.
          !
        enddo
        !
      endif
      !
    end do
    end do
    end do
    !
    call dataswap(lsolid,timerept=ltimrpt)
    !
    deallocate(lss)
    !
    ! allocate(rss(0:im,0:jm,0:km))
    ! !
    ! do k=0,km
    ! do j=0,jm
    ! do i=0,im
    !   !
    !   if(lsolid(i,j,k)) then
    !     rss(i,j,k)=1.d0
    !   else
    !     rss(i,j,k)=0.d0
    !   endif
    !   !
    ! end do
    ! end do
    ! end do
    ! !
    ! call tecbin('testout/tec_solid_sensor'//mpirankname//'.plt',      &
    !                                   x(0:im,0:jm,0:km,1),'x',    &
    !                                   x(0:im,0:jm,0:km,2),'y',    &
    !                                   x(0:im,0:jm,0:km,3),'z',    &
    !                                                   rss,'ss' )
    ! !
    ! deallocate(rss)
    !
    ! call mpistop
    !
    return
    !
  end subroutine ShockSolid
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! End of the subroutine ShockSolid
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !+-------------------------------------------------------------------+
  !| This function is to determin is a i,j,k is in the domain          |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 23-07-2021  | Created by J. Fang @ Warrington                     |
  !+-------------------------------------------------------------------+
  pure logical function ijkin(i,j,k)
    !
    use commvar, only : im,jm,km
    !
    integer,intent(in) :: i,j,k
    !
    if(i<0 .or. j<0 .or. k<0 .or. i>im .or. j>jm .or. k>km) then
      ijkin=.false.
    else
      ijkin=.true.
    endif
    !
    return
    !
  end function ijkin
  !+-------------------------------------------------------------------+
  !| The end of the function ijkin.                                    |
  !+-------------------------------------------------------------------+
  !!
  !+-------------------------------------------------------------------+
  !| This function is to determin is a i,j,k is in the domain          |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 24-08-2021  | Created by J. Fang @ Warrington                     |
  !+-------------------------------------------------------------------+
  pure logical function ijkcellin(i,j,k)
    !
    use commvar, only : im,jm,km,ndims
    !
    integer,intent(in) :: i,j,k
    !
    if(i>=1 .and. i<=im .and. j>=1 .and. j<=jm) then
      !
      if(ndims==3) then
        if(k>=1 .and. k<=km) then
          ijkcellin=.true.
        else
          ijkcellin=.false.
        endif
      elseif(ndims==2) then
        ijkcellin=.true.
      endif
    else
      ijkcellin=.false.
    endif
    !
    return
    !
  end function ijkcellin
  !+-------------------------------------------------------------------+
  !| The end of the function ijkcellin.                                |
  !+-------------------------------------------------------------------+

  !
end module commcal
!+---------------------------------------------------------------------+
!| The end of the module commcal.                                      |
!+---------------------------------------------------------------------+