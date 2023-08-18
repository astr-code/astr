!+---------------------------------------------------------------------+
!| This module contains utility subroutines                            |
!+---------------------------------------------------------------------+
!| CHANGE RECORD                                                       |
!| -------------                                                       |
!| 10-08-2022  | Created by J. Fang                                    |
!+---------------------------------------------------------------------+
module utility
  !
  use stlaio,  only: get_unit
  !
  implicit none
  !
  contains
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to report time cost by each subroutine.   |
  !+-------------------------------------------------------------------+
  !| note: should only be called from one rank, usually the root       |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 12-02-2021  | Created by J. Fang @ Warrington                     |
  !| 10-08-2022  | Moved to this module, and called by subroutines     |
  !|             | by J. Fang @ Warrington.                            |
  !+-------------------------------------------------------------------+
  subroutine timereporter(routine,message,timecost,mode)
    !
    use commvar, only : nstep,maxstep,ctime,flowtype,conschm,          &
                        difschm,rkscheme,ia,ja,ka,preptime,nsrpt
    ! arguments
    character(len=*),intent(in),optional :: routine,message,mode
    real(8),intent(inout),optional :: timecost
    !
    ! local data
    !
    type :: trep
      character(len=16) :: rout,mode
      character(len=64) :: mesg
      character(len=1) :: cate
      integer :: order
      real(8) :: time
    end type trep
    !
    logical :: lexist
    integer :: i,ios,n,varorder
    integer,save :: counter=0
    integer,save :: hand_rp,repsp
    logical,save :: linit=.true.
    character(len=16) :: realmode,charinput
    character(len=64) :: messinput,messoutp
    character(len=20),save :: rptfname
    real(8),save :: total_time=1.d-10
    real(8) :: percent,datainput,commtime,iotime,xtratime,vartime
    !
    integer,parameter :: nmax=100
    !
    type(trep),save :: recorder(nmax)
    !
    if(linit) then
      !
      rptfname='time_report.'//message
      !
      inquire(file=rptfname, exist=lexist)
      !
      if(lexist) call system('mv -v '//rptfname//' '//rptfname//'.bak')
      !
      call system('echo "----------------------------------------------------------------" '//rptfname)
      call system('echo "CPU infomation" >> '//rptfname)
      call system('echo "----------------------------------------------------------------" >> '//rptfname)
      call system('lscpu | grep "Model name" >> '//rptfname)
      call system('lscpu | grep "CPU MHz" >> '//rptfname)
      call system('lscpu | grep "Socket(s):" >> '//rptfname)
      call system('lscpu | grep "Core(s) per socket:" >> '//rptfname)
      call system('lscpu | grep "Thread(s) per core:" >> '//rptfname)
      call system('lscpu | grep "cache" >> '//rptfname)
      call system('echo "----------------------------------------------------------------" >> '//rptfname)
      !
      hand_rp=get_unit()
      !
      open(hand_rp,file=rptfname,position="append")
      write(hand_rp,'(A)')'  statistic of computing time'
      write(hand_rp,'(A,A)')'     flowtype: ',trim(flowtype)
      write(hand_rp,'(A,A)')'  conv scheme: ',trim(conschm)
      write(hand_rp,'(A,A)')'  diff scheme: ',trim(difschm)
      write(hand_rp,'(A,A)')'    rk scheme: ',trim(rkscheme)
      write(hand_rp,'(4(A,I0))')'    grid size: ',ia,' x ',ja,' x ', &
                                         ka,' = ',(ia+1)*(ja+1)*(ka+1)
      ! write(hand_rp,'(2X,62A)')('-',i=1,62)
      !
      close(hand_rp)
      print*,' << ',rptfname
      !
      repsp=nstep
      !
      linit=.false.
      !
      return
      !
    endif
    !
    if(repsp==nstep) return
    !
    if(present(mode)) then
      realmode=mode
    else
      realmode='general'
    endif
    !
    if(trim(realmode)=='final') then
      hand_rp=get_unit()
      open (hand_rp,file=rptfname,position="append")
      write(hand_rp,'(2X,62A)')('-',i=1,62)
      write(hand_rp,'(2X,2(A20))')'total nsteps','computational time'
      write(hand_rp,'(2X,3(A20))')('----------',i=1,2)
      write(hand_rp,'(2X,I20,E20.6E2,10X,A)')nstep-1,timecost
      close(hand_rp)
      print*,' << ',rptfname
    else
      !
      do n=1,counter
        !
        if(present(routine)) then
         !
         if(trim(recorder(n)%rout)==routine) then
            !
            if(trim(recorder(n)%mesg)==message) then
                !
                ! it is the report from the same subroutine
                !
                recorder(n)%time=recorder(n)%time+timecost
                !
                exit
                !
            endif
            !
         endif
         !
        endif
        !
      enddo
      !
      if(n==counter+1) then
        ! normal exit from the previous step
        !
        counter=counter+1
        !
        if(present(routine))  recorder(counter)%rout=routine
        if(present(mode))     recorder(counter)%mode=mode
        if(present(message))  recorder(counter)%mesg=message
        !
        recorder(counter)%time=timecost
        !
      endif
      !
      if(counter==nmax) then
        print*,' !! WARNING MAX counter reached @ timereporter'
      endif
      !
      if(trim(routine)=='steploop') then
        ! this is the last subroutine reporting
        !
        total_time=timecost
        !
        total_time=max(total_time,1.d-10)
        !
        hand_rp=get_unit()
        open (hand_rp,file=rptfname,position="append")
        write(hand_rp,'(2X,62A)')('-',i=1,62)
        write(hand_rp,'(2X,A20,I7,A3,I7)')'    nsteps  : ',repsp,' - ',nstep
        ! write(hand_rp,'(7X,55A)')('-',i=1,55)
        ! write(hand_rp,'(2X,A16,A14,A11,A20)')'subroutine','time cost','%','note'
        !
        commtime=0.d0
        xtratime=0.d0
        iotime=0.d0
        vartime=0.d0
        varorder=1
        do n=1,counter
          !
          if(trim(recorder(n)%rout)=='qswap'             .or. &
             trim(recorder(n)%rout)=='array3d_sendrecv'  .or. & 
             trim(recorder(n)%rout)=='array4d_sendrecv'  .or. &  
             trim(recorder(n)%rout)=='updatable_rel2d_'  .or. & 
             trim(recorder(n)%rout)=='updatable_rel_a2'  .or. & 
             trim(recorder(n)%rout)=='updatable_rel2d_'  ) then
           !
           recorder(n)%cate='m'
           !
           commtime=commtime+recorder(n)%time
           !
          elseif(trim(recorder(n)%rout)=='writechkpt') then
           !
           recorder(n)%cate='o'
           !
           iotime=iotime+recorder(n)%time
           !
          else
           !
           recorder(n)%cate='x'
           !
           xtratime=xtratime+recorder(n)%time
           !
          endif
          !
        enddo
        !
        ! output the message part
        if(commtime>1.d-10) then
          !
          percent=commtime/total_time*100.d0
          write(hand_rp,'(5X,57A)')('-',i=1,57)
          write(hand_rp,'(2X,A16,E14.5E2,3X,F7.2,A)')'Comm Time',commtime,percent,'%'
          write(hand_rp,'(2X,A16,A14,A11,A18)')('--------',i=1,4)
          !
          do n=1,counter
            !
            percent=recorder(n)%time/total_time*100.d0
            !
            if(recorder(n)%cate=='m') then
              !
              write(hand_rp,'(2X,A16,E14.5E2,3X,F7.2,A,10X,A)')trim(recorder(n)%rout), &
                                     recorder(n)%time,percent,'%',trim(recorder(n)%mesg)
            endif
            !
          enddo
          !
        endif
        !
        if(iotime>1.d-10) then
          !
          ! output the other part
          percent=iotime/total_time*100.d0
          write(hand_rp,'(5X,57A)')('-',i=1,57)
          write(hand_rp,'(2X,A16,E14.5E2,3X,F7.2,A)')'IO Time',xtratime,percent,'%'
          write(hand_rp,'(2X,A16,A14,A11,A18)')('--------',i=1,4)
          !
          do n=1,counter
            !
            percent=recorder(n)%time/total_time*100.d0
            !
            if(recorder(n)%cate=='o') then
              !
              write(hand_rp,'(2X,A16,E14.5E2,3X,F7.2,A,10X,A)')trim(recorder(n)%rout), &
                                     recorder(n)%time,percent,'%',trim(recorder(n)%mesg)
            endif
            !
          enddo
          !
        endif
        !
        if(xtratime>1.d-10) then
          !
          ! output the other part
          percent=100.d0-(iotime+commtime)/total_time*100.d0
          write(hand_rp,'(5X,57A)')('-',i=1,57)
          write(hand_rp,'(2X,A16,E14.5E2,3X,F7.2,A)')'Other Part',xtratime,percent,'%'
          write(hand_rp,'(2X,A16,A14,A11,A18)')('--------',i=1,4)
          !
          do n=1,counter
            !
            percent=recorder(n)%time/total_time*100.d0
            !
            if(recorder(n)%cate=='x') then
              !
              write(hand_rp,'(2X,A16,E14.5E2,3X,F7.2,A,10X,A)')trim(recorder(n)%rout), &
                                     recorder(n)%time,percent,'%',trim(recorder(n)%mesg)
            endif
            !
          enddo
          !
        endif
        !
        close(hand_rp)
        !
        print*,' << ',rptfname
        !
        total_time=timecost
        !
        counter=0
        !
        repsp=nstep
        !
        do n=1,nmax
          recorder(n)%rout=''
          recorder(n)%mode=''
          recorder(n)%mesg=''
          recorder(n)%time=0.d0
        enddo
        !
      endif
      !
    endif
    !
    timecost=0.d0
    !
  end subroutine timereporter
  !+-------------------------------------------------------------------+
  !| The end of the subroutine timerept.                               |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to init a text file, either create a new file  |
  !| to resume an old file.                                            |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 17-Aug-2023: Created by J. Fang @ Appleton                        |
  !+-------------------------------------------------------------------+
  subroutine listinit(filename,handle,firstline)
    !
    use commvar,   only: lrestart,nstep
    use strings,   only: split
    !
    character(len=*),intent(in) :: filename
    integer,intent(out) :: handle
    character(len=*),intent(in),optional :: firstline
    !
    character(len=16),allocatable :: args(:)
    logical :: fex
    integer :: nargs,ns,ferr,n
    character(len=120) :: txtformat
    !
    inquire(file=filename,exist=fex)
    handle=get_unit()
    !
    open(handle,file=filename)
    !
    if(lrestart .and. fex) then
      ! resume a file
      ns=0
      read(handle,*)
      ! first line is alway a text 
      !
      do while(ns<nstep)
        !
        read(handle,*,iostat=ferr)ns
        !
        if(ferr< 0) then
          print*,' ** ns,nstep=',ns,nstep
          print*,' ** end of file is reached.'
          exit
        endif
        !
      enddo
      !
      backspace(handle)
      write(*,'(A,I0)')'   ** resume'//filename//'at step: ',ns
      !
    else
      ! create a file
      call split(firstline,args,' ')
      !
      nargs=size(args)
      !
      write(txtformat,'(A,I0,A)')'(',nargs,'(1X,A20))'

      write(handle,txtformat)(trim(args(n)),n=1,nargs)
      !
    endif
    !
  end subroutine listinit
  !+-------------------------------------------------------------------+
  !| The end of the subroutine listinit.                               |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to write listing data.                         |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 17-Aug-2023: Created by J. Fang @ Appleton                        |
  !+-------------------------------------------------------------------+
  subroutine listwrite(handle,var1,var2,var3,var4,var5,var6,var7,var8)
    !
    use commvar, only: nstep,time
    !
    integer,intent(in) :: handle
    real(8),intent(in),optional :: var1,var2,var3,var4,var5,var6,var7,var8
    !
    if(present(var8)) then
      write(handle,"(1X,I20,9(1X,E20.13E2))")nstep,time,var1,var2,    &
                                           var3,var4,var5,var6,var7,var8
    elseif(present(var7)) then
      write(handle,"(1X,I20,8(1X,E20.13E2))")nstep,time,var1,var2,    &
                                               var3,var4,var5,var6,var7
    elseif(present(var6)) then
      write(handle,"(1X,I20,7(1X,E20.13E2))")nstep,time,var1,var2,    &
                                               var3,var4,var5,var6
    elseif(present(var5)) then
      write(handle,"(1X,I20,6(1X,E20.13E2))")nstep,time,var1,var2,    &
                                               var3,var4,var5
    elseif(present(var4)) then
      write(handle,"(1X,I20,5(1X,E20.13E2))")nstep,time,var1,var2,    &
                                               var3,var4
    elseif(present(var3)) then
      write(handle,"(1X,I20,4(1X,E20.13E2))")nstep,time,var1,var2,    &
                                               var3
    elseif(present(var2)) then
      write(handle,"(1X,I20,3(1X,E20.13E2))")nstep,time,var1,var2
    elseif(present(var1)) then
      write(handle,"(1X,I20,2(1X,E20.13E2))")nstep,time,var1
    else
      stop ' !! error @ listwrite'
    endif
    !
  end subroutine listwrite
  !+-------------------------------------------------------------------+
  !| The end of the subroutine listwrite.                              |
  !+-------------------------------------------------------------------+
  !
end module utility
!+---------------------------------------------------------------------+
!| The end of the module utility                                       |
!+---------------------------------------------------------------------+