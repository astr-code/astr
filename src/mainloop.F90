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
  use parallel, only: lio,mpistop,mpirank,qswap,dataswap,mpirankname,  &
                      pmax,ptime
  use commvar,  only: im,jm,km
  use tecio
  !
  implicit none
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
    use commvar,  only: maxstep,time,deltat,nstep,nwrite,ctime,nlstep
    use readwrite,only: readcont,timerept
    !
    ! local data
    integer :: counter
    real(8) :: var1
    !
    var1=ptime()
    !
    counter=nwrite
    do while(nstep<=maxstep)
      !
      call rk3(counter)
      !
      if(counter==nwrite) then
        !
        call readcont
        !
        ctime(2)=ptime()-var1
        !
        call timerept
        !
        counter=0
        !
      endif
      !
      counter=counter+1
      nstep=nstep+1
      time=time+deltat
      !
    enddo
    !
    ctime(1)=ptime()-var1

    !
  end subroutine steploop
  !+-------------------------------------------------------------------+
  !| The end of the subroutine steploop.                               |
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
  subroutine rk3(counter)
    !
    use commvar,  only : im,jm,km,numq,deltat,lfilter,nstep,nwrite,    &
                         ctime
    use commarray,only : x,q,qrhs,rho,vel,prs,tmp,spc,jacob
    use fludyna,  only : q2fvar
    use solver,   only : rhscal,filterq
    use statistic,only : statcal,statout
    use readwrite,only : output
    !
    integer,intent(in) :: counter
    !
    ! logical data
    logical,save :: firstcall = .true.
    real(8),save :: rkcoe(3,3)
    integer :: irk,m
    real(8) :: time_beg,time_beg_rhs,time_beg_sta,time_beg_io
    real(8),allocatable :: qsave(:,:,:,:)
    !
#ifdef cputime
    time_beg=ptime()
#endif
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
    do m=1,numq
      qsave(0:im,0:jm,0:km,m)=q(0:im,0:jm,0:km,m)*jacob(0:im,0:jm,0:km)
    enddo
    !
    do irk=1,3
      !
      call qswap
      !
#ifdef cputime
      time_beg_rhs=ptime()
#endif
      call rhscal
#ifdef cputime
      ctime(4)=ctime(4)+ptime()-time_beg_rhs
#endif
      !
      if(irk==1) then
        !
#ifdef cputime
        time_beg_sta=ptime()
#endif
        !
        call statcal
        !
        call statout
        !
#ifdef cputime
        ctime(5)=ctime(5)+ptime()-time_beg_sta
#endif
        !
        if(counter==nwrite) then
#ifdef cputime
          time_beg_io=ptime()
#endif
          !
          call output
          !
#ifdef cputime
          ctime(6)=ctime(6)+ptime()-time_beg_io
#endif
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
      if(lfilter) call filterq
      !
      call q2fvar(q=q(0:im,0:jm,0:km,:),                               &
                                     density=rho(0:im,0:jm,0:km),      &
                                    velocity=vel(0:im,0:jm,0:km,:),    &
                                    pressure=prs(0:im,0:jm,0:km),      &
                                 temperature=tmp(0:im,0:jm,0:km),      &
                                     species=spc(0:im,0:jm,0:km,:)     )
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
end module mainloop
!+---------------------------------------------------------------------+
!| The end of the module readwrite.                                    |
!+---------------------------------------------------------------------+