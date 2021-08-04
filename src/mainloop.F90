!+---------------------------------------------------------------------+
!| This module contains subroutines in hte main computational loop.    |
!+---------------------------------------------------------------------+
!| CHANGE RECORD                                                       |
!| -------------                                                       |
!| 09-02-2021  | Created by J. Fang                                    |
!+---------------------------------------------------------------------+
module mainloop
  !
  use constdef
  use parallel, only: lio,mpistop,mpirank,qswap,mpirankname,pmax,      &
                      ptime,irk,jrk,irkm,jrkm
  use commvar,  only: im,jm,km,ia,ja,ka
  use tecio
  use stlaio,  only: get_unit
  !
  implicit none
  !
  integer :: loop_counter=0
  !
  contains
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to advance the solution.                       |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 09-02-2021  | Created by J. Fang @ Warrington                     |
  !+-------------------------------------------------------------------+
  subroutine steploop
    !
    use commvar,  only: maxstep,time,deltat,nstep,nwrite,ctime,        &
                        nlstep,rkscheme
    use readwrite,only: readcont,timerept
    use commcal,  only: cflcal
    !
    ! local data
    real(8) :: time_start,time_beg
    !
    time_start=ptime()
    !
    do while(nstep<=maxstep)
      !
      time_beg=ptime()
      !
      if(rkscheme=='rk3') then
        call rk3
      elseif(rkscheme=='rk4') then
        call rk4
      endif
      !
      ctime(2)=ctime(2)+ptime()-time_beg
      !
      if(loop_counter==nwrite .or. loop_counter==0) then
        !
        call readcont
        !
        call cflcal(deltat)
        !
        call timerept
        !
        loop_counter=0
        !
      endif
      !
      loop_counter=loop_counter+1
      nstep=nstep+1
      time=time+deltat
      !
    enddo
    !
    ctime(1)=ptime()-time_start
    !
    ! call errest
    !
  end subroutine steploop
  !+-------------------------------------------------------------------+
  !| The end of the subroutine steploop.                               |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to calculate error to estimate order of        |
  !| accuracy.                                                         |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 09-02-2021  | Created by J. Fang @ Warrington                     |
  !+-------------------------------------------------------------------+
  subroutine errest
    !
    use commvar, only : flowtype,xmin,xmax,ymin,ymax,zmin,zmax,nstep,  &
                        time
    use commarray,only: x,vel,rho,prs,spc,tmp,q,acctest_ref
    use parallel, only: psum,pmax,mpirankname
    !
    real(8) :: l1error,l2error,lineror
    !
    integer :: i,j,k,fh
    real(8) :: xc,yc,zc,rvor,radi2,var1
    !
    if(trim(flowtype)=='accutest') then
      !
      ! error calculation
      l1error=0.d0
      l2error=0.d0
      lineror=0.d0
      !
      do i=0,im
        !
        var1=(acctest_ref(i)-spc(i,0,0,1))
        !
        l1error=l1error+abs(var1)
        l2error=l2error+var1**2
        lineror=max(lineror,abs(var1))
        !
        ! if(lio) print*,i,acctest_ref(i),spc(i,0,0,1)
        !
      enddo
      !
      ! print*,l1error,dble(ia)
      l1error=psum(l1error)/dble(ia)
      l2error=sqrt(psum(l2error)/dble(ia))
      lineror=pmax(lineror)
      !
      if(lio) then
        print*,' ** nstep= ',nstep,'time= ',time
        write(*,'(2X,A7,3(1X,A13))')'mx','L1_error','L2_error','L∞_error'
        write(*,'(2X,I7,3(1X,E13.6E2))')ia,l1error,l2error,lineror
      endif
      !
      fh=get_unit()
      open(fh,file='prof'//mpirankname//'.dat')
      do i=0,im
        write(fh,*)x(i,0,0,1),acctest_ref(i),spc(i,0,0,1)
      enddo
      close(fh)
      print*,' << prof',mpirankname,'.dat ... done.'
      !
    endif
    !
  end subroutine errest
  !+-------------------------------------------------------------------+
  !| The end of the subroutine errest.                                 |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine advances the field solution in time using 3-step  |
  !| 3rd-rder Rungle-Kutta scheme.                                     |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 27-Nov-2018: Created by J. Fang @ STFC Daresbury Laboratory       |
  !+-------------------------------------------------------------------+
  subroutine rk3
    !
    use commvar,  only : im,jm,km,numq,deltat,lfilter,nstep,nwrite,    &
                         ctime,hm,lavg,navg,nstep
    use commarray,only : x,q,qrhs,rho,vel,prs,tmp,spc,jacob
    use fludyna,  only : q2fvar
    use solver,   only : rhscal,filterq,spongefilter
    use statistic,only : statcal,statout,meanflowcal,nsamples,liosta
    use readwrite,only : output
    use bc,       only : boucon
    !
    ! logical data
    logical,save :: firstcall = .true.
    real(8),save :: rkcoe(3,3)
    integer :: irk,i,j,k,m
    real(8) :: time_beg,time_beg_rhs,time_beg_sta,time_beg_io
    real(8),allocatable :: qsave(:,:,:,:)
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
      firstcall=.false.
      !
    endif
    !
    allocate(qsave(0:im,0:jm,0:km,1:numq))
    !
    do irk=1,3
      !
      qrhs=0.d0
      !
      call qswap(ctime(7))
      !
      call boucon
      !
      call rhscal(ctime(4))
      !
      if(irk==1) then
        !
        do m=1,numq
          qsave(0:im,0:jm,0:km,m)=q(0:im,0:jm,0:km,m)*jacob(0:im,0:jm,0:km)
        enddo
        !
        if(loop_counter.ne.0) then
          !
          call statcal(ctime(5))
          !
          call statout
          !
          if(lavg) then
            if(mod(nstep,navg)==0) call meanflowcal(ctime(5))
          else
            nsamples=0
            liosta=.false.
          endif
          !
        endif
        !
        !
        if(loop_counter==nwrite .or. loop_counter==0) then
          !
          call output(ctime(6))
          !
        endif
        !
      endif
      !
      do m=1,numq
        q(0:im,0:jm,0:km,m)=rkcoe(1,irk)*qsave(0:im,0:jm,0:km,m)+      &
                            rkcoe(2,irk)*q(0:im,0:jm,0:km,m)*          &
                                     jacob(0:im,0:jm,0:km)+            &
                            rkcoe(3,irk)*qrhs(0:im,0:jm,0:km,m)*deltat
        !
        q(0:im,0:jm,0:km,m)=q(0:im,0:jm,0:km,m)/jacob(0:im,0:jm,0:km)
      enddo
      !
      if(lfilter) call filterq(ctime(8))
      !
      call spongefilter
      !
      call updatefvar
      !
      call crashcheck
      !
      ! call tecbin('testout/tecinit'//mpirankname//'.plt',              &
      !                                    x(0:im,0:jm,-hm:km+hm,1),'x', &
      !                                    x(0:im,0:jm,-hm:km+hm,2),'y', &
      !                                    x(0:im,0:jm,-hm:km+hm,3),'z', &
      !                                  rho(0:im,0:jm,-hm:km+hm),'ro',  &
      !                                  vel(0:im,0:jm,-hm:km+hm,1),'u', &
      !                                  vel(0:im,0:jm,-hm:km+hm,2),'v', &
      !                                  prs(0:im,0:jm,-hm:km+hm),'p',   &
      !                                  spc(0:im,0:jm,-hm:km+hm,1),'T' )
      ! call mpistop
      !
    enddo
    !
    deallocate(qsave)
    !
    ctime(3)=ctime(3)+ptime()-time_beg
    !
    return
    !
  end subroutine rk3
  !+-------------------------------------------------------------------+
  !| The end of the subroutine rk3.                                    |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to update flow variables from q.               |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 04-Aug-2018: Created by J. Fang @ STFC Daresbury Laboratory       |
  !+-------------------------------------------------------------------+
  subroutine updatefvar
    !
    use commarray,only : q,rho,vel,prs,tmp,spc,tke,omg
    use commvar,  only : num_species,num_modequ,turbmode
    use fludyna,  only : q2fvar
    !
    if(trim(turbmode)=='k-omega') then
      !
      call q2fvar(q=q(0:im,0:jm,0:km,:),                               &
                                     density=rho(0:im,0:jm,0:km),      &
                                    velocity=vel(0:im,0:jm,0:km,:),    &
                                    pressure=prs(0:im,0:jm,0:km),      &
                                 temperature=tmp(0:im,0:jm,0:km),      &
                                     species=spc(0:im,0:jm,0:km,:),    &
                                         tke=tke(0:im,0:jm,0:km),      &
                                       omega=omg(0:im,0:jm,0:km) )
      !
    elseif(trim(turbmode)=='none') then
      !
      call q2fvar(q=q(0:im,0:jm,0:km,:),                               &
                                     density=rho(0:im,0:jm,0:km),      &
                                    velocity=vel(0:im,0:jm,0:km,:),    &
                                    pressure=prs(0:im,0:jm,0:km),      &
                                 temperature=tmp(0:im,0:jm,0:km),      &
                                     species=spc(0:im,0:jm,0:km,:) )
      !
    else
      print*,' !! ERROR @ updatefvar'
      stop
    endif
    !
  end subroutine updatefvar
  !+-------------------------------------------------------------------+
  !| The end of the subroutine updatefvar.                             |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine advances the field solution in time using 4-step  |
  !| 4th-rder Rungle-Kutta scheme.                                     |
  !+-------------------------------------------------------------------+
  !| ref: A. Dipankar, T.K. Sengupta, Symmetrized compact scheme for   |
  !| receptivity study of 2D transitional channel flow, Journal of     |
  !| Computational Physics 215 (2006) 245–273.                         |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 27-Nov-2018: Created by J. Fang @ STFC Daresbury Laboratory       |
  !+-------------------------------------------------------------------+
  subroutine rk4
    !
    use commvar,  only : im,jm,km,numq,deltat,lfilter,nstep,nwrite,    &
                         ctime,hm,lavg,navg,nstep,limmbou
    use commarray,only : x,q,qrhs,rho,vel,prs,tmp,spc,jacob
    use fludyna,  only : q2fvar
    use solver,   only : rhscal,filterq,spongefilter
    use statistic,only : statcal,statout,meanflowcal,liosta,nsamples
    use readwrite,only : output,writemon
    use bc,       only : boucon,immbody
    !
    !
    ! logical data
    logical,save :: firstcall = .true.
    real(8),save :: rkcoe(2,4)
    integer :: nrk,i,j,k,m
    real(8) :: time_beg,time_beg_rhs,time_beg_sta,time_beg_io
    real(8),allocatable :: qsave(:,:,:,:),rhsav(:,:,:,:)
    !
    time_beg=ptime() 
    !
    if(firstcall) then
      !
      rkcoe(1,1)=0.5d0
      rkcoe(1,2)=0.5d0
      rkcoe(1,3)=1.d0
      rkcoe(1,4)=num1d6
      !
      rkcoe(2,1)=1.d0
      rkcoe(2,2)=2.d0
      rkcoe(2,3)=2.d0
      rkcoe(2,4)=1.d0
      !
      firstcall=.false.
      !
    endif
    !
    allocate(qsave(0:im,0:jm,0:km,1:numq),rhsav(0:im,0:jm,0:km,1:numq))
    !
    do nrk=1,4
      !
      qrhs=0.d0
      !
      if(limmbou) call immbody(ctime(11))
      !
      call qswap(ctime(7))
      !
      call boucon
      !
      call rhscal(ctime(4))
      !
      if(nrk==1) then
        !
        do m=1,numq
          qsave(0:im,0:jm,0:km,m)=q(0:im,0:jm,0:km,m)*jacob(0:im,0:jm,0:km)
          rhsav(0:im,0:jm,0:km,m)=0.d0
        enddo
        !
        if(nstep==0 .or. loop_counter.ne.0) then
          !
          call statcal(ctime(5))
          !
          call statout
          !
          call writemon
          !
          if(lavg) then
            if(mod(nstep,navg)==0) call meanflowcal(ctime(5))
          else
            nsamples=0
            liosta=.false.
          endif
          !
        else
          ! call statcal(ctime(5))
          ! call statout
        endif
        !
        if(loop_counter==nwrite) then
          !
          call output(ctime(6))
          !
        endif
        !
      endif
      !
      if(nrk<=3) then
        do m=1,numq
          q(0:im,0:jm,0:km,m)=qsave(0:im,0:jm,0:km,m)+                 &
                              rkcoe(1,nrk)*deltat*qrhs(0:im,0:jm,0:km,m)
          !
          q(0:im,0:jm,0:km,m)=q(0:im,0:jm,0:km,m)/jacob(0:im,0:jm,0:km)
          !
          rhsav(0:im,0:jm,0:km,m)=rhsav(0:im,0:jm,0:km,m)+             &
                                  rkcoe(2,nrk)*qrhs(0:im,0:jm,0:km,m)
        enddo
      else
        do m=1,numq
          q(0:im,0:jm,0:km,m)=qsave(0:im,0:jm,0:km,m)+                 &
                                  rkcoe(1,nrk)*deltat*(                &
                                 qrhs(0:im,0:jm,0:km,m)+               &
                                 rhsav(0:im,0:jm,0:km,m) )
          !
          q(0:im,0:jm,0:km,m)=q(0:im,0:jm,0:km,m)/jacob(0:im,0:jm,0:km)
          !
        enddo
      endif
      !
      if(lfilter) call filterq(ctime(8))
      !
      call spongefilter
      !
      call q2fvar(q=q(0:im,0:jm,0:km,:),                               &
                                     density=rho(0:im,0:jm,0:km),      &
                                    velocity=vel(0:im,0:jm,0:km,:),    &
                                    pressure=prs(0:im,0:jm,0:km),      &
                                 temperature=tmp(0:im,0:jm,0:km),      &
                                     species=spc(0:im,0:jm,0:km,:)     )
      !
      call crashcheck
      !
    enddo
    !
    deallocate(qsave,rhsav)
    !
    ctime(3)=ctime(3)+ptime()-time_beg
    !
    return
    !
  end subroutine rk4
  !+-------------------------------------------------------------------+
  !| The end of the subroutine rk4.                                    |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to check if the computational is crashed.      |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 27-Nov-2018: Created by J. Fang @ Warrington                      |
  !+-------------------------------------------------------------------+
  subroutine crashcheck
    !
    use commvar,   only : hm
    use commarray, only : q,rho,tmp,prs,x
    use parallel, only : por
    !
    ! local data
    integer :: i,j,k,l
    logical :: ltocrash
    !
    ltocrash=.false.
    !
    do k=0,km
    do j=0,jm
    do i=0,im
      !
      if(rho(i,j,k)>=0.d0) then
        continue
      else
        print*,' !! non-positive density identified !!'
        write(*,'(2(A,I0),2(2X,I0),A,3(1X,E13.6E2))')'   ** mpirank= ',&
                           mpirank,' i,j,k= ',i,j,k,' x,y,z=',x(i,j,k,:)
        write(*,'(A,5(1X,E13.6E2))')'   ** q= ',q(i,j,k,:)
        ltocrash=.true.
      endif
      !
      if(tmp(i,j,k)>=0.d0) then
        continue
      else
        print*,' !! non-positive temperature identified !!'
        write(*,'(2(A,I0),2(2X,I0),A,3(1X,E13.6E2))')'   ** mpirank= ',&
                           mpirank,' i,j,k= ',i,j,k,' x,y,z=',x(i,j,k,:)
        write(*,'(A,5(1X,E13.6E2))')'   ** q= ',q(i,j,k,:)
        do l=-hm,hm
          print*,l,x(i,j+l,k,2),q(i,j+l,k,5),tmp(i,j+l,k)
        enddo
        ltocrash=.true.
      endif
      !
      if(prs(i,j,k)>=0.d0) then
        continue
      else
        print*,' !! non-positive pressure identified !!'
        write(*,'(2(A,I0),2(2X,I0),A,3(1X,E13.6E2))')'   ** mpirank= ',&
                           mpirank,' i,j,k= ',i,j,k,' x,y,z=',x(i,j,k,:)
        write(*,'(A,5(1X,E13.6E2))')'   ** q= ',q(i,j,k,:)
        ltocrash=.true.
      endif
      !
      if(q(i,j,k,5)>=0.d0) then
        continue
      else
        print*,' !! non-positive energy identified !!'
        write(*,'(2(A,I0),2(2X,I0),A,3(1X,E13.6E2))')'   ** mpirank= ',&
                           mpirank,' i,j,k= ',i,j,k,' x,y,z=',x(i,j,k,:)
        write(*,'(A,5(1X,E13.6E2))')'   ** q= ',q(i,j,k,:)
        ltocrash=.true.
      endif
      !
    enddo
    enddo
    enddo
    !
    ltocrash=por(ltocrash)
    !
    if(ltocrash) then
      if(lio) print*,' !! COMPUTATION CRASHED !!'
      call mpistop
    endif
    !
    return
    !
  end subroutine crashcheck
  !+-------------------------------------------------------------------+
  !| The end of the subroutine crashcheck.                             |
  !+-------------------------------------------------------------------+
  !
end module mainloop
!+---------------------------------------------------------------------+
!| The end of the module readwrite.                                    |
!+---------------------------------------------------------------------+