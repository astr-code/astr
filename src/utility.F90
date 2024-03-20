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
  !| Progress indicators library.                                      |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 31-01-2024  | copied by J. Fang via:                              |
  !| https://github.com/macie/fortran-libs                             |
  !|  Maciej Żok, 2010 MIT License                                     |
  !+-------------------------------------------------------------------+
  subroutine progress_bar(iteration,maximum,info2show,barlength)
     !
     ! Prints progress bar.
     !
     ! Args: 
     !     iteration - iteration number
     !     maximum - total iterations
     !     barlength - length of the bar
     !   
     ! use iso_fortran_env
     integer,intent(in) :: iteration,maximum
     character(len=*),intent(in),optional :: info2show
     integer,intent(in),optional :: barlength
     integer :: counter,nlength
     integer :: done
     real(4) :: perc
     !
     if(present(barlength)) then
         nlength=barlength
     else
         nlength=10
     endif
     !
     perc = 100.0*real(iteration)/real(maximum)
     done = floor(perc/(100.0/real(nlength)))  ! mark length
     !
     write(6,'(1A1,A,A)',advance='no')char(13),info2show,'['
     if (done .LE. 0) then
         do counter = 1, nlength
             write(6,'(1A1,A)',advance='no')'='
         end do
     else if ((done .GT. 0) .and. (done .LT. nlength)) then
         do counter = 1, done
             write(6,'(1A1,A)',advance='no')'>'
         end do
         do counter = done+1, nlength
             write(6,'(1A1,A)',advance='no')'='
         end do 
     else
         do counter = 1, nlength
             write(6,'(1A1,A)',advance='no')'>'
         end do
     end if
     write(6,'(A,F5.1,A)',advance='no')'] ',perc,'%'
     !
     if(iteration==maximum) write(6,*)
     !
  end subroutine progress_bar
  !+-------------------------------------------------------------------+
  !| The end of the subroutine progress_bar.                           |
  !+-------------------------------------------------------------------+
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
      !
      print*,nargs,txtformat
      !
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
  subroutine listwrite(handle,var1,var2,var3,var4,var5,var6,var7,var8,var9,var10,var11,var12)
    !
    use commvar, only: nstep,time,ref_tim
    !
    integer,intent(in) :: handle
    real(8),intent(in),optional :: var1,var2,var3,var4,var5,var6,      &
                                   var7,var8,var9,var10,var11,var12

    !
    if(present(var12)) then
      write(handle,"(1X,I20,13(1X,E20.13E2))")nstep,time/ref_tim,       &
           var1,var2,var3,var4,var5,var6,var7,var8,var9,var10,var11,var12
    elseif(present(var11)) then
      write(handle,"(1X,I20,12(1X,E20.13E2))")nstep,time/ref_tim,       &
                 var1,var2,var3,var4,var5,var6,var7,var8,var9,var10,var11
    elseif(present(var10)) then
      write(handle,"(1X,I20,11(1X,E20.13E2))")nstep,time/ref_tim,       &
                       var1,var2,var3,var4,var5,var6,var7,var8,var9,var10
    elseif(present(var9)) then
      write(handle,"(1X,I20,10(1X,E20.13E2))")nstep,time/ref_tim,       &
                       var1,var2,var3,var4,var5,var6,var7,var8,var9
    elseif(present(var8)) then
      write(handle,"(1X,I20,9(1X,E20.13E2))")nstep,time/ref_tim,       &
                       var1,var2,var3,var4,var5,var6,var7,var8
    elseif(present(var7)) then
      write(handle,"(1X,I20,8(1X,E20.13E2))")nstep,time/ref_tim,       &
                       var1,var2,var3,var4,var5,var6,var7
    elseif(present(var6)) then
      write(handle,"(1X,I20,7(1X,E20.13E2))")nstep,time/ref_tim,       &
                       var1,var2,var3,var4,var5,var6
    elseif(present(var5)) then
      write(handle,"(1X,I20,6(1X,E20.13E2))")nstep,time/ref_tim,       &
                       var1,var2,var3,var4,var5
    elseif(present(var4)) then
      write(handle,"(1X,I20,5(1X,E20.13E2))")nstep,time/ref_tim,       &
                       var1,var2,var3,var4
    elseif(present(var3)) then
      write(handle,"(1X,I20,4(1X,E20.13E2))")nstep,time/ref_tim,       &
                       var1,var2,var3
    elseif(present(var2)) then
      write(handle,"(1X,I20,3(1X,E20.13E2))")nstep,time/ref_tim,       &
                       var1,var2
    elseif(present(var1)) then
      write(handle,"(1X,I20,2(1X,E20.13E2))")nstep,time/ref_tim,       &
                       var1
    else
      stop ' !! error @ listwrite'
    endif
    !
  end subroutine listwrite
  !+-------------------------------------------------------------------+
  !| The end of the subroutine listwrite.                              |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This function Verifies that a character string represents a       |
  !|  numerical value                                                  |
  !+-------------------------------------------------------------------+
  ! ref: http://fcode.cn/code_gen-115-1.html
  !+-------------------------------------------------------------------+
  Integer Function IsNum(zval)
    ! 确定字符是否是数值类型：
    ! 0-非数值的字符串
    ! 1-整数(integer)
    ! 2-小数(fixed point real)
    ! 3-指数类型实数(exponent type real)
    ! 4-双精度实数指数形式(exponent type double)
    Character (Len=*), Intent (In) :: zval
    !
    Integer :: num, nmts, nexp, kmts, ifexp, ichr
    !
    Integer, Parameter :: kint = 1 ! integer
    Integer, Parameter :: kfix = 2 ! fixed point real
    Integer, Parameter :: kexp = 3 ! exponent type real
    Integer, Parameter :: kdbl = 4 ! exponent type double
    !
    ! initialise
    num = 0  ! 数字的格式，最后传递给ISNUM返回
    nmts = 0 ! 整数或浮点数的数字个数
    nexp = 0 ! 指数形式的数字个数
    kmts = 0 ! 有+-号为1，否则为0
    ifexp = 0! 似乎没用
    ! loop over characters
    ichr = 0
    !
    Do
    
      If (ichr>=len(zval)) Then
    
        ! last check
    
        If (nmts==0) Exit
    
        If (num>=kexp .And. nexp==0) Exit
    
        isnum = num
    
        Return
    
      End If
    
      ichr = ichr + 1
    
      Select Case (zval(ichr:ichr))
    
        ! process blanks
    
      Case (' ')
    
        Continue
    
        ! process digits
    
      Case ('0', '1', '2', '3', '4', '5', '6', '7', '8', '9')
    
        If (num==0) num = kint
    
        If (num<kexp) Then
    
          nmts = nmts + 1
    
          ! 整数或浮点数+1
    
        Else
    
          nexp = nexp + 1
    
          ! 指数形式+1
    
        End If
    
        ! process signs
    
      Case ('+', '-')
    
        If (num==0) Then
    
          If (kmts>0) Exit
    
          ! 出现2个符号，非数字
    
          kmts = 1
    
          num = kint
    
        Else
    
          If (num<kexp) Exit
    
          If (ifexp>0) Exit
    
          ifexp = 1
    
        End If
    
        ! process decimal point
    
      Case ('.')
    
        If (num/=kint .And. ichr/=1) Exit
    
        ! 前面不是整数，小数点也不是第一个字符，则非数字
    
        num = kfix
    
        ! process exponent
    
      Case ('e', 'E')
    
        If (num>=kexp) Exit
    
        If (nmts==0) Exit
    
        num = kexp
      Case ('d', 'D')
    
        If (num>=kexp) Exit
    
        If (nmts==0) Exit
    
        num = kdbl
    
        ! any other character means the string is non-numeric
    
      Case Default
    
        Exit
    
      End Select
    
    End Do
    
    ! if this point is reached, the string is non-numeric
    
    isnum = 0
    
    Return
    
  End Function IsNum
  !+-------------------------------------------------------------------+
  !| The end of the Function IsNum.                                    |
  !+-------------------------------------------------------------------+
  !
  function rnorm_box_muller(mode) result(variates) 
    !
    use constdef
    ! coded formulas from https://en.wikipedia.org/wiki/Box%E2%80%93Muller_transform
    !
    character(len=*),intent(in) :: mode
    !
    ! return two uncorrelated standard normal variates
    integer,allocatable :: seed(:)
    integer :: rantime(8)
    !
    integer :: n
    real(8) :: variates(2)
    real(8) :: u(2), factor, arg
    !
    logical,save :: firstcall=.true.
    !
    if(mode=='sync') then
      ! all processor generate same random numbers
      if(firstcall) then
        !
        call random_seed(size = n)
        allocate(seed(n))
        call date_and_time(values=rantime) 
        !  use date and minutes for synthetisation 
        !    1      2    3    4     5      6      7    8
        !-----------------------------------------------
        ! 2023     10    5   60    23      6     28  962
        !-----------------------------------------------
        ! year  month date   ??  hour minute second msec
        !
        seed=0
        seed(1:6)=rantime(1:6)
        !
        call random_seed(put=seed)
        !
        deallocate(seed)
        !
        firstcall=.false.
      endif
      !
    endif
    !
    do
       call random_number(u)
       if (u(1) > 0.d0) exit
    end do
    factor = sqrt(-2 * log(u(1)))
    arg = 2.d0*pi*u(2)
    variates = factor * [cos(arg),sin(arg)]
    !
  end function rnorm_box_muller
  !
end module utility
!+---------------------------------------------------------------------+
!| The end of the module utility                                       |
!+---------------------------------------------------------------------+