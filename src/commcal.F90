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
    use commvar,  only: im,jm,km
    use commarray,only: vel,rho,prs,tmp,dxi,jacob
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
      css=sos(tmp(i,j,k))
      csi=css*sqrt(dxi(i,j,k,1,1)**2+dxi(i,j,k,1,2)**2+dxi(i,j,k,1,3)**2)
      csj=css*sqrt(dxi(i,j,k,2,1)**2+dxi(i,j,k,2,2)**2+dxi(i,j,k,2,3)**2)
      csk=css*sqrt(dxi(i,j,k,3,1)**2+dxi(i,j,k,3,2)**2+dxi(i,j,k,3,3)**2)
      !
      deltai=max(ubar,ubar-csi,ubar+csi)
      deltaj=max(vbar,vbar-csj,vbar+csj)
      deltak=max(wbar,wbar-csk,wbar+csk)
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
      write(*,"(A24,1x,F13.7)")'     current time step: ',deltat
      write(*,"(A24,1x,F13.7)")'           current CFL: ',cfl
      write(*,"(A24,1x,F13.7)")'   time step for CFL=1: ',deltat/cfl
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
    use commvar,  only : im,jm,km,hm
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
    call dataswap(lss)
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
    call dataswap(lsolid)
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

  !
end module commcal
!+---------------------------------------------------------------------+
!| The end of the module commcal.                                      |
!+---------------------------------------------------------------------+