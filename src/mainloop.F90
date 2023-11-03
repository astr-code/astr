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
  use commvar,  only: im,jm,km,ia,ja,ka,ctime,nstep,lcracon,lreport,   &
                      ltimrpt,rkstep
  use commarray,only: crinod
  use tecio
  use stlaio,   only: get_unit
  use utility,  only: timereporter
  !
  implicit none
  !
  integer :: loop_counter=0
  integer :: nxtbakup,feqbakup=5
  integer :: fhand_err
  integer :: nstep0
  real(8) :: time_start
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
    use commvar,  only: maxstep,time,deltat,feqchkpt,feqwsequ,feqlist, &
                        rkscheme,nsrpt,flowtype,limmbou
    use readwrite,only: readcont,timerept,nxtchkpt,nxtwsequ
    use commcal,  only: cflcal
    use ibmethod, only: ibforce
    use userdefine,only: udf_eom_set
    !
    ! local data
    real(8) :: time_beg,time_next_step
    integer :: hours,minus,secod
    logical,save :: firstcall = .true.
    integer,dimension(8) :: value
    real(8) :: time_total,time_dowhile
    !
    time_start=ptime()
    time_total  =0.d0
    time_dowhile=0.d0
    !
    nstep0=nstep
    !
    if(firstcall) then
      !
      crinod=.false.
      firstcall=.false.
      !
    endif
    !
    nxtchkpt=nstep+feqchkpt
    nxtwsequ=nstep+feqwsequ
    nxtbakup=nstep+feqbakup
    !
    nsrpt   =nstep
    !
    fhand_err=get_unit()
    open(fhand_err,file='errnode.log')
    !
    if(lio) call timereporter(routine='steploop',timecost=time_dowhile,  &
                              message='init file')
    !
    do while(nstep<=maxstep)
      !
      time_beg=ptime()
      !
      call crashcheck
      !
      if(rkscheme=='rk3') then
        call rk3(timerept=ltimrpt)
      elseif(rkscheme=='rk4') then
        call rk4
      endif
      !
      time_dowhile=time_dowhile+ptime()-time_beg
      !
      if(loop_counter==feqchkpt .or. loop_counter==0) then
        !
        call readcont
        !
        nxtchkpt=nstep+feqchkpt
        !
        call cflcal(deltat)
        !
        if(lio) then
          !
          write(*,'(A,I0)',advance='no')'  ** next checkpoint at step : ',nxtchkpt
          !
          if(loop_counter==0) then
            write(*,*)''
          else
            !
            time_next_step=(ctime(2)-ctime(6))*dble(feqchkpt)/dble(nstep-nsrpt)
            hours=int(time_next_step/3600.d0)
            minus=int((time_next_step-3600.d0*dble(hours))/60.d0)
            secod=time_next_step-3600.d0*dble(hours)-60.d0*dble(minus)
            !
            write(*,'(3(A,I0),A)',advance='no')', after ',hours,'h:',  &
                                               minus,'m:',secod,'s'
            call date_and_time(VALUES=value)
            !
            secod=secod+value(7)
            if(secod>=60) then
              secod=secod-60
              minus=minus+1
            endif
            minus=minus+value(6)
            if(minus>=60) then
              minus=minus-60
              hours=hours+1
            endif
            hours=hours+value(5)
            if(hours>=24) then
              hours=hours-24
            endif
            !
            write(*,'(3(A,I0))')', estimated at ',hours,':',minus,':',secod
            !
            time_next_step=ctime(2)*dble(maxstep-nstep)/dble(nstep-nsrpt)
            hours=int(time_next_step/3600.d0)
            minus=int((time_next_step-3600.d0*dble(hours))/60.d0)
            secod=time_next_step-3600.d0*dble(hours)-60.d0*dble(minus)
            !
            write(*,'(3(A,I0),A)',advance='no')'  ** job ends after ', &
                                        hours,'h:',minus,'m:',secod,'s'
            call date_and_time(VALUES=value)
            !
            secod=secod+value(7)
            if(secod>=60) then
              secod=secod-60
              minus=minus+1
            endif
            minus=minus+value(6)
            if(minus>=60) then
              minus=minus-60
              hours=hours+1
            endif
            hours=hours+value(5)
            if(hours>=24) then
              hours=hours-24
            endif
            !
            write(*,'(3(A,I0))')', estimated at ',hours,':',minus,':',secod
            !
            nsrpt=nstep
            !
            call timereporter(routine='steploop',timecost=time_dowhile,  &
                              message='main loop')
          endif
          !
        endif
        !
        ! call timerept
        !
        ! if(limmbou) call ibforce
        !
        loop_counter=0
        !
      endif
      !
      call udf_eom_set
      !
      loop_counter=loop_counter+1
      nstep=nstep+1
      time=time+deltat
      !
    enddo
    !
    ! if(limmbou) call timerept
    !
    ! call ibforce
    !
    time_total=ptime()-time_start
    !
    if(lio) call timereporter(timecost=time_total,mode='final')
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
    use commvar, only : flowtype,xmin,xmax,ymin,ymax,zmin,zmax,time
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
  subroutine rk3(timerept)
    !
    use commvar,  only : im,jm,km,numq,deltat,lfilter,feqchkpt,hm,     &
                         lavg,feqavg,nstep,limmbou,turbmode,feqslice,  &
                         feqwsequ,lwslic,lreport,flowtype,     &
                         ndims,num_species,maxstep
    use commarray,only : x,q,qrhs,rho,vel,prs,tmp,spc,jacob
    use fludyna,  only : updatefvar
    use comsolver,only : filterq,spongefilter,filter2e
    use solver,   only : rhscal
    use bc,       only : boucon,immbody
#ifdef COMB
    use thermchem,only : imp_euler_ode,heatrate
    use fdnn
    use commvar,  only : odetype
#endif 
    !
    ! argument
    logical,intent(in),optional :: timerept
    !
    ! local data
    logical,save :: firstcall = .true.
    real(8),save :: rkcoe(3,3)
    integer :: i,j,k,m,n
    real(8) :: time_beg,time_beg_rhs,time_beg_sta,time_beg_io
    real(8),allocatable,save :: qsave(:,:,:,:)
    integer :: dt_ratio,jdnn,idnn
    real(8) :: hrr,time_beg_2
    real(8),save :: subtime=0.d0
    !
    time_beg=ptime()
    !
#ifdef COMB
    if(odetype=='dnn') then 
      if(nstep==nstep0) then
        call initialization
        call initialize_locell(im,jm,km)
      endif 
      dnnidx(:,:,:)=1
      jdnn=0
      dt_ratio=delta_t/deltat
    endif 
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
      allocate(qsave(0:im,0:jm,0:km,1:numq))
      !
      firstcall=.false.
      !
    endif
    !
    do rkstep=1,3
      !
      if( (loop_counter==feqchkpt .or. loop_counter==0) .and. rkstep==1 ) then
        lreport=.true.
      else
        lreport=.false.
      endif
      !
      qrhs=0.d0
      !
      if(limmbou) call immbody(timerept=ltimrpt)
      !
      if(flowtype(1:2)/='0d') call qswap(timerept=ltimrpt)
      !
      if(flowtype(1:2)/='0d') call boucon
      !
      if(flowtype(1:2)/='0d') call qswap(timerept=ltimrpt)
      !
      call rhscal(timerept=ltimrpt)
      !
      !for debug
      ! do k=0,km
      ! do j=0,jm
      ! do i=0,im
      !   !
      !   do n=1,numq
      !     if(isnan(qrhs(i,j,k,n))) then
      !       print*,'!! have NaN in qrhs !! the ',n,'one is wrong, and rk = ',rkstep
      !       stop
      !     endif
      !   enddo
      !   !
      ! enddo
      ! enddo
      ! enddo
      !
      if(flowtype(1:2)=='0d') jacob=1.d0
      !
      time_beg_2=ptime()
      !
      if(rkstep==1) then
        !
        do m=1,numq
          qsave(0:im,0:jm,0:km,m)=q(0:im,0:jm,0:km,m)*jacob(0:im,0:jm,0:km)
        enddo
        !
        call rkfirst
        !
      endif
      !
      do m=1,numq
        !
        q(0:im,0:jm,0:km,m)=rkcoe(1,rkstep)*qsave(0:im,0:jm,0:km,m)+      &
                            rkcoe(2,rkstep)*q(0:im,0:jm,0:km,m)*          &
                                     jacob(0:im,0:jm,0:km)+            &
                            rkcoe(3,rkstep)*qrhs(0:im,0:jm,0:km,m)*deltat
        !
        q(0:im,0:jm,0:km,m)=q(0:im,0:jm,0:km,m)/jacob(0:im,0:jm,0:km)
        !
      enddo
      !
      ctime(14)=ctime(14)+ptime()-time_beg_2
      !
      if(lfilter) then
        call filterq(timerept=ltimrpt)
      endif
      !
      if(flowtype(1:2)/='0d') call spongefilter
      !
      time_beg_2=ptime()
      !
      call updatefvar
      !
      ! for debug
      ! do k=0,km
      ! do j=0,jm
      ! do i=0,im
      !   !
      !   print*,qrhs(i,j,k,:)
      !   !
      ! enddo
      ! enddo
      ! enddo
      !
      ! call mpistop
      !
      ctime(15)=ctime(15)+ptime()-time_beg_2
      !
      if(lcracon) call crashfix(ctime(16))
      !
#ifdef COMB
      if(rkstep==3) then
        !
        time_beg_2=ptime()
        !
        do i=0,im 
        do j=0,jm
        do k=0,km
              !
          if(odetype=='ime' .or. odetype=='rk3') then
            !
            if(odetype=='ime') &
              call imp_euler_ode(rho(i,j,k),tmp(i,j,k),spc(i,j,k,:),deltat)
            !
            q(i,j,k,6:numq)=max(0.d0,spc(i,j,k,:)*rho(i,j,k))
            q(i,j,k,6:numq)=q(i,j,k,1)*q(i,j,k,6:numq)/sum(q(i,j,k,6:numq))
            !
          elseif(odetype=='dnn' .and. mod(nstep,dt_ratio)==0 .and. nstep>0) then
            !
            ! if(locell(jcell)%tmp<1000.d0 .or. locell(jcell)%tmp>2200.d0) then 
            !
            hrr=heatrate(rho(i,j,k),tmp(i,j,k),spc(i,j,k,:))
            if(hrr<1.d8) then 
              !
              dnnidx(i,j,k)=0
              !
              do idnn=1,dt_ratio
                !
                call imp_euler_ode(rho(i,j,k),tmp(i,j,k),spc(i,j,k,:),deltat)
                !
              enddo 
              !
            else 
              !
              jdnn=jdnn+1
              inputholder(1,jdnn)=tmp(i,j,k)
              inputholder(2,jdnn)=1.0
              inputholder(3:num_species+1,jdnn)=spc(i,j,k,1:num_species-1)
              !
            endif
            !
          else
            !
            continue
            !
          endif !odetype
          !
        enddo !i
        enddo !j
        enddo !k
        !
        if(odetype=='dnn' .and. nstep>0 .and. mod(nstep,dt_ratio)==0) then
          !
          allocate(input(n_layer(1),jdnn),output(n_layer(1),jdnn))
          !
          input(:,:)=inputholder(:,1:jdnn)
          if(jdnn/=0) output=netOneStep(input,epoch,n_layer,jdnn)
          !
          idnn=0
          do i=0,im
          do j=0,jm
          do k=0,km
            !
            if(dnnidx(i,j,k)/=0) then 
              idnn=idnn+1
              spc(i,j,k,1:num_species-1)=output(3:num_species+1,idnn)
              spc(i,j,k,num_species)=0.d0
            endif 
            !
            q(i,j,k,6:numq)=max(0.d0,spc(i,j,k,:)*rho(i,j,k))
            q(i,j,k,6:numq)=q(i,j,k,1)*q(i,j,k,6:numq)/sum(q(i,j,k,6:numq))
            !
          enddo 
          enddo
          enddo
          !
        endif 
        !
        ctime(13)=ctime(13)+ptime()-time_beg_2
        !
        time_beg_2=ptime()
        !
        call updatefvar
        !
        ctime(15)=ctime(15)+ptime()-time_beg_2
        !
      endif !odetype
#endif
    !
    enddo !rk
    !
#ifdef COMB
    if(odetype=='dnn') then 
      if(nstep==maxstep) call finalize()
      if(allocated(input))deallocate(input,output)
    endif
    !
#endif
    !
    if(present(timerept) .and. timerept) then
      !
      subtime=subtime+ptime()-time_beg
      !
      if(lio .and. loop_counter==feqchkpt) call timereporter(routine='rk3',   &
                                                            timecost=subtime)
    endif
    !
    return
    !
  end subroutine rk3
  !+-------------------------------------------------------------------+
  !| The end of the subroutine rk3.                                    |
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
    use commvar,  only : im,jm,km,numq,deltat,lfilter,feqchkpt,hm,     &
                         lavg,feqavg,nstep,limmbou,turbmode,feqslice,  &
                         feqwsequ,lwslic
    use commarray,only : x,q,qrhs,rho,vel,prs,tmp,spc,jacob
    use fludyna,  only : updatefvar
    use comsolver,only : filterq,spongefilter
    use solver,   only : rhscal
    use bc,       only : boucon,immbody
    !
    !
    ! local data
    logical,save :: firstcall = .true.
    real(8),save :: rkcoe(2,4)
    integer :: i,j,k,m
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
    do rkstep=1,4
      !
      qrhs=0.d0
      !
      if(limmbou) call immbody(timerept=ltimrpt)
      !
      call qswap(timerept=ltimrpt)
      !
      call boucon
      !
      call qswap(timerept=ltimrpt)
      !
      call rhscal(timerept=ltimrpt)
      !
      if(rkstep==1) then
        !
        do m=1,numq
          qsave(0:im,0:jm,0:km,m)=q(0:im,0:jm,0:km,m)*jacob(0:im,0:jm,0:km)
          rhsav(0:im,0:jm,0:km,m)=0.d0
        enddo
        !
        call rkfirst
        !
      endif
      !
      if(rkstep<=3) then
        do m=1,numq
          q(0:im,0:jm,0:km,m)=qsave(0:im,0:jm,0:km,m)+                 &
                              rkcoe(1,rkstep)*deltat*qrhs(0:im,0:jm,0:km,m)
          !
          q(0:im,0:jm,0:km,m)=q(0:im,0:jm,0:km,m)/jacob(0:im,0:jm,0:km)
          !
          rhsav(0:im,0:jm,0:km,m)=rhsav(0:im,0:jm,0:km,m)+             &
                                  rkcoe(2,rkstep)*qrhs(0:im,0:jm,0:km,m)
        enddo
      else
        do m=1,numq
          q(0:im,0:jm,0:km,m)=qsave(0:im,0:jm,0:km,m)+                 &
                                  rkcoe(1,rkstep)*deltat*(                &
                                 qrhs(0:im,0:jm,0:km,m)+               &
                                 rhsav(0:im,0:jm,0:km,m) )
          !
          q(0:im,0:jm,0:km,m)=q(0:im,0:jm,0:km,m)/jacob(0:im,0:jm,0:km)
          !
        enddo
      endif
      !
      if(lfilter) call filterq(timerept=ltimrpt)
      !
      call spongefilter
      !
      call updatefvar
      !
      ! call crashcheck
      if(lcracon) call crashfix
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
  !| This subroutine is to conducte operation in the fist step of rk   |
  !| temporal advance.                                                 |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 28-Dec-2021: Created by J. Fang @ Warrington                      |
  !+-------------------------------------------------------------------+
  subroutine rkfirst
    !
    use commvar,  only : lavg,lwslic,lwsequ,feqavg,feqchkpt,feqwsequ,  &
                         feqslice,iomode
    use statistic,only : statcal,statout,meanflowcal,liosta,nsamples
    use readwrite,only : writechkpt,writemon,writeslice,writeflfed,    &
                         nxtchkpt,nxtwsequ
    use userdefine,only: udf_stalist,udf_write
    !
    ! local data
    integer,save :: nxtavg
    logical,save :: firstcall = .true.
    !
    if(firstcall) then
      nxtavg=nstep+feqavg
      firstcall = .false.
    endif
    !
    call statcal(timerept=ltimrpt)
    !
    call statout(time_start)
    !
    call udf_stalist
    !
    if(nstep==0 .or. loop_counter.ne.0) then
      ! the first step after reading ehecking out doesn't need to do this
      !
      call writemon
      !
      if(lavg) then
        if(nstep==nxtavg) then
          call meanflowcal(timerept=ltimrpt)
          !
          nxtavg=nstep+feqavg
        endif
      else
        nsamples=0
        liosta=.false.
        nxtavg=nstep+feqavg
      endif
      !
      if(lwslic .and. mod(nstep,feqslice)==0) then
        call writeslice(ctime(23))
      endif
      !
    endif
    !
    ! time to write checkpoint
    if(iomode == 'n') then
      !
      continue
      !
    else
      !
      if(nstep==nxtchkpt) then
        !
        ! the checkpoint and flowfield may be writen in the same time
        call writechkpt(nxtwsequ,timerept=ltimrpt)
        !
        call udf_write
        !
      endif
      !
      if(lwsequ .and. nstep==nxtwsequ) then
        !
        call writeflfed(timerept=ltimrpt)
        !
      endif
      !
    endif
    !
    return
    !
  end subroutine rkfirst
  !+-------------------------------------------------------------------+
  !| The end of the subroutine rkfirst.                                |
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
    use commvar,  only: hm
    use commarray,only: q,rho,tmp,prs,x,nodestat
    use parallel, only: por
    use fludyna,  only: updateq
    use readwrite,only: readcheckpoint
    use statistic,only: nsamples
    !
    ! local data
    integer :: i,j,k,l,fh,ii,jj,kk
    logical :: ltocrash
    !
    ltocrash=.false.
    !
    fh=get_unit()
    !
    ! open(fh,file='badpoint'//mpirankname//'.dat')
    do k=0,km
    do j=0,jm
    do i=0,im
      !
      if(nodestat(i,j,k)<=0.d0) then
        !
        ! only check fluid points
        if(q(i,j,k,1)>=0.d0) then
          continue
        else
          !
          crinod(i,j,k)=.true.
          !
          write(fhand_err,'(A,I0)')'error node cant wiped at nstep=',nstep
          write(fhand_err,'(2(A,I0),2(2X,I0),A,3(1X,E13.6E2))')'mpirank= ',mpirank, &
                                               ' i,j,k= ',i,j,k,' x,y,z=',x(i,j,k,:)
          write(fhand_err,'(A,I0)')'nodestat: ',nodestat(i,j,k)
          write(fhand_err,'(A,5(1X,E13.6E2))')'q= ',q(i,j,k,1:5)
          !
          do kk=-1,1
          do jj=-1,1
          do ii=-1,1
            if(ii==0 .and. jj==0 .and. kk==0) then
              continue
            else
              write(fhand_err,'(3(1X,I0),5(1X,E13.6E2),1X,I0)')ii,jj,kk,  &
                             q(i+ii,j+jj,k+kk,1:5),nodestat(i+ii,j+jj,k+kk)
            endif
          enddo
          enddo
          enddo
          !
          ltocrash=.true.
          !
          ! write(fh,'(3(1X,E13.6E2))')x(i,j,k,:)
        endif
        !
      endif
      !
    enddo
    enddo
    enddo
    ! close(fh)
    ! print*,' << badpoint',mpirankname,'.dat'
    !
    ltocrash=por(ltocrash)
    !
    if(ltocrash) then
      !
      if(lcracon) then
        if(lio) print*,' !! COMPUTATION CRASHED !!'
        if(lio) print*,' !! FETCH AN BAKUP FLOW FIELD !!'
        !
        call databakup('recovery')
        ! call readcheckpoint(folder='bakup')
        !
        loop_counter=0
      else
        if(lio) print*,' !! COMPUTATION CRASHED, JOB STOPS !!'
        call mpistop
      endif
      !
    elseif(lcracon) then
      !
      ! if not crash, backup data
      !
      if(nstep==nxtbakup) then
        !
        call databakup('backup')
        !
        nxtbakup=nstep+feqbakup
        !
        feqbakup=feqbakup*2
        !
        feqbakup=min(500,feqbakup) ! the max frequency of back is set to 500
        !
      endif
      !
    endif
    !
    return
    !
  end subroutine crashcheck
  !+-------------------------------------------------------------------+
  !| The end of the subroutine crashcheck.                             |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to backup or recovery data                     |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 06-04-2022: Created by J. Fang @ Warrington                       |
  !+-------------------------------------------------------------------+
  subroutine databakup(mode)
    !
    use commvar,  only:  filenumb,fnumslic,time,flowtype,force,numq
    use bc,       only : ninflowslice
    use statistic,only : massflux,massflux_target
    use commarray,only : q 
    use fludyna,  only : updatefvar
    !
    ! argument
    character(len=*),intent(in) :: mode
    !
    ! lcoal data
    type :: databack
      integer :: nstep,filenumb,fnumslic,ninflowslice
      real(8) :: time,massflux,massflux_target,force(3)
      real(8),allocatable :: q(:,:,:,:)
    end type databack
    !
    type(databack),save :: dat_a,dat_b
    character(len=1),save :: datpnt='o'
    !
    if(mode=='backup') then
      !
      if(datpnt=='o') datpnt='a'
      !
      if(datpnt=='a') then
        !
        dat_a%nstep        =nstep       
        dat_a%filenumb     =filenumb    
        dat_a%fnumslic     =fnumslic    
        dat_a%ninflowslice =ninflowslice
        dat_a%time         =time        
        if(flowtype=='channel') then
          dat_a%massflux       =massflux
          dat_a%massflux_target=massflux_target
          dat_a%force          =force
        endif
        !
        if(.not. allocated(dat_a%q)) then
          allocate(dat_a%q(0:im,0:jm,0:km,1:numq))
          !
          dat_a%q(0:im,0:jm,0:km,1:numq)=q(0:im,0:jm,0:km,1:numq)
        endif
        !
        if(lio) write(*,'(A,I0)')'  ** data backed at nstep= ',nstep
        !
        datpnt='b'
        !
      elseif(datpnt=='b') then
        !
        dat_b%nstep        =nstep       
        dat_b%filenumb     =filenumb    
        dat_b%fnumslic     =fnumslic    
        dat_b%ninflowslice =ninflowslice
        dat_b%time         =time        
        if(flowtype=='channel') then
          dat_b%massflux       =massflux
          dat_b%massflux_target=massflux_target
          dat_b%force          =force
        endif
        !
        if(.not. allocated(dat_b%q)) then
          allocate(dat_b%q(0:im,0:jm,0:km,1:numq))
          !
          dat_b%q(0:im,0:jm,0:km,1:numq)=q(0:im,0:jm,0:km,1:numq)
        endif
        !
        if(lio) write(*,'(A,I0)')'  ** data backed at nstep= ',nstep
        !
        datpnt='a'
        !
      else
        print*,' !! datpnt: ',datpnt
        stop ' !! error 1 of datpnt @ databakup'
      endif
      !
    elseif(mode=='recovery') then
      !
      if(datpnt=='o') then
        print*,' !! not backup data avaliable !!'
        stop
      elseif(datpnt=='a') then
        !
        nstep       =dat_a%nstep       
        filenumb    =dat_a%filenumb    
        fnumslic    =dat_a%fnumslic    
        ninflowslice=dat_a%ninflowslice
        time        =dat_a%time        
        if(flowtype=='channel') then
          massflux       =dat_a%massflux       
          massflux_target=dat_a%massflux_target
          force          =dat_a%force          
        endif
        !
        q(0:im,0:jm,0:km,1:numq)=dat_a%q(0:im,0:jm,0:km,1:numq)
        !
        if(lio) write(*,'(A,I0,A)')'  ** data recovered to nstep= ',nstep,' from dat_a'
        !
        datpnt='b'
        !
      elseif(datpnt=='b') then
        !
        nstep       =dat_b%nstep       
        filenumb    =dat_b%filenumb    
        fnumslic    =dat_b%fnumslic    
        ninflowslice=dat_b%ninflowslice
        time        =dat_b%time        
        if(flowtype=='channel') then
          massflux       =dat_b%massflux       
          massflux_target=dat_b%massflux_target
          force          =dat_b%force          
        endif
        !
        q(0:im,0:jm,0:km,1:numq)=dat_b%q(0:im,0:jm,0:km,1:numq)
        !
        if(lio) write(*,'(A,I0,A)')'  ** data recovered to nstep= ',nstep,' from dat_b'
        !
        datpnt='a'
        !
      else
        print*,' !! datpnt: ',datpnt
        stop ' !! error 2 of datpnt @ databakup'
      endif
      !
      call updatefvar
      !
    else
      stop ' !! mode error @ databakup'
    endif
    !
  end subroutine databakup
  !+-------------------------------------------------------------------+
  !| The end of the subroutine databakup.                              |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to wipe the point where the result is not good.|
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 04-10-2020: Created by J. Fang @ Warrington                       |
  !+-------------------------------------------------------------------+
  subroutine crashfix(subtime)
    !
    use commvar,   only : numq,lreport,nondimen,spcinf
    use commarray, only : q,rho,tmp,vel,prs,spc,x,nodestat
    use parallel,  only : por,ig0,jg0,kg0,psum,ptime
    use fludyna,   only : q2fvar,thermal
    !
    ! arguments
    real(8),intent(inout),optional :: subtime
    !
    ! local data
    integer :: i,j,k,l,fh,ii,jj,kk
    real(8) :: qavg(numq)
    integer :: norm,counter
    !
    logical,save :: firstcall = .true.
    real(8),save :: eps_rho,eps_prs,eps_tmp,time_beg
    integer,save :: step_normal = 0
    !
    if(present(subtime)) time_beg=ptime()
    !
    if(firstcall) then
      eps_rho=1.d-5
      eps_tmp=1.d-5
      !
      if(nondimen) then 
        eps_prs=thermal(density=eps_rho,temperature=eps_tmp)
      else 
        eps_prs=thermal(density=eps_rho,temperature=eps_tmp,species=spcinf)
      endif 
    endif
    !
    counter=0
    !
    do k=0,km
    do j=0,jm
    do i=0,im
      !
      if(nodestat(i,j,k)<=0.d0) then
        ! only check fluid points
        !
        if(rho(i,j,k)>=eps_rho .and. prs(i,j,k)>=eps_prs .and.         &
           tmp(i,j,k)>=eps_tmp) then
          continue
        else
          !
          crinod(i,j,k)=.true.
          !
          ! print*,' !! non-positive density/energy identified !!'
          ! write(*,'(2(A,I0),2(2X,I0),A,3(1X,E13.6E2))')'   ** mpirank= ',&
          !                    mpirank,' i,j,k= ',i,j,k,' x,y,z=',x(i,j,k,:)
          ! write(*,'(A,I0)')'   ** nodestat: ',nodestat(i,j,k)
          ! write(*,'(A,5(1X,E13.6E2))')'   ** q= ',q(i,j,k,1:5)
          ! !
          qavg=0.d0
          norm=0
          do kk=-1,1
          do jj=-1,1
          do ii=-1,1
            !
            if(ii==0 .and. jj==0 .and. kk==0) then
              continue
            elseif( ig0+i+ii<0 .or. ig0+i+ii>ia .or.  &
                    jg0+j+jj<0 .or. jg0+j+jj>ja ) then
              continue
            else
              if(rho(i+ii,j+jj,k+kk)>=0.d0 .and. prs(i+ii,j+jj,k+kk)>=0.d0 &
                  .and. tmp(i+ii,j+jj,k+kk)>=0.d0) then
                qavg(:)=qavg(:)+q(i+ii,j+jj,k+kk,:)
                norm=norm+1
              endif
            endif
            !
          enddo
          enddo
          enddo
          !
          if(norm>=1) then
            !
            write(fhand_err,'(A,I0)')' error nodes identified and wiped at nstep=',nstep
            write(fhand_err,'(2(A,I0),2(2X,I0),A,3(1X,E13.6E2))')'mpirank= ', mpirank,    &
                                                    ' i,j,k= ',i,j,k,' x,y,z=',x(i,j,k,:)
            write(fhand_err,'(3(A,E13.6E2))')'rho=',rho(i,j,k),' prs=',prs(i,j,k),        &
                                                               ' tmp=',tmp(i,j,k)
            write(fhand_err,'(2(A,5(1X,E13.6E2)))')'q= ',q(i,j,k,1:5),' -> ',qavg/dble(norm)
            !
            q(i,j,k,:)=qavg/dble(norm)
            !
            call q2fvar(q=q(i,j,k,:),                       &
                                       density=rho(i,j,k),  &
                                      velocity=vel(i,j,k,:),&
                                      pressure=prs(i,j,k),  &
                                   temperature=tmp(i,j,k),  &
                                       species=spc(i,j,k,:) )
            !
            counter=counter+1
            !
          endif
          !
        endif
        !
      endif
      !
    enddo
    enddo
    enddo
    !
    counter=psum(counter)
    !
    if(counter==0) then
      step_normal=step_normal+1
    else
      step_normal=0
    endif
    !
    if(lio .and. counter>1) then
      write(*,'(A,I0,A)')'  !! ',counter,' error nodes were wiped.'
    endif
    !
    if(lreport) then
      !
      if(lio) then
        if(step_normal>0) then
          write(*,'(A,I0,A)')'  ** the comput. been normal for: ',     &
                                                    step_normal,' steps'
        endif
      endif
      !
      counter=0
      do k=0,km
      do j=0,jm
      do i=0,im
        !
        if(crinod(i,j,k)) then
          counter=counter+1
        endif
        !
      enddo
      enddo
      enddo
      !
      counter=psum(counter)
      !
      if(lio) then
        write(*,'(A,I0,A)')'  ** ',counter,' critical nodes'
      endif
      !
    endif
    !
    if(present(subtime)) subtime=subtime+ptime()-time_beg
    !
    return
    !
  end subroutine crashfix
  !+-------------------------------------------------------------------+
  !| The end of the subroutine crashfix.                               |
  !+-------------------------------------------------------------------+
  !
end module mainloop
!+---------------------------------------------------------------------+
!| The end of the module readwrite.                                    |
!+---------------------------------------------------------------------+
