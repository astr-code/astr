!+---------------------------------------------------------------------+
!| This module contains subroutine of getting command via interactive  |
!+---------------------------------------------------------------------+
!| ==============                                                      |
!| CHANGE RECORD                                                       |
!| -------------                                                       |
!|  28-05-2021  | Created by J. Fang @ Warrington                      |
!+---------------------------------------------------------------------+
module cmdefne
  !
  implicit none
  !
  contains
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to input command from keyboard.           |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 11-02-2021  | Created by J. Fang @ Warrington                     |
  !| 28-05-2021  | Moved to the module command by J. Fang @ Warrington |
  !+-------------------------------------------------------------------+
  subroutine readkeyboad(inputfile,cmd,casename)
    !
    character(len=*),intent(out),optional :: inputfile,cmd,casename
    !
    ! local data
    integer :: ierr,cli_len,nkey,nlen,arg_count
    character(len=128) :: keyin
    !
    nkey=0
    cli_len=1
    !
    if(present(inputfile)) inputfile='.' ! default value
    !
    do while(cli_len>0) 
      !
      nkey=nkey+1
      call get_command_argument(nkey,keyin,cli_len,ierr)
      !
      if(trim(keyin)=='-input') then
        nkey=nkey+1
        call get_command_argument(nkey,keyin,cli_len,ierr)
        inputfile=trim(keyin)
      endif
      !
      if(trim(keyin)=='-cmd') then
        nkey=nkey+1
        call get_command_argument(nkey,keyin,cli_len,ierr)
        cmd=trim(keyin)
        !
        if(trim(cmd)=='init') then
          nkey=nkey+1
          call get_command_argument(nkey,keyin,cli_len,ierr)
          casename=trim(keyin)
        endif
        !
      endif
      !
    enddo
    !
  end subroutine readkeyboad
  !+-------------------------------------------------------------------+
  !| The end of the subroutine readkeyboad.                            |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to get command from the main program.     |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 28-05-2021  | Created by J. Fang @ Warrington                     |
  !+-------------------------------------------------------------------+
  subroutine getcmd(cmd,casename)
    !
    use parallel, only: mpistop,mpirank,bcast
    !
    character(len=*),intent(out) :: cmd
    character(len=*),intent(out),optional :: casename
    !
    cmd='run' 
    ! default value
    !
    if(mpirank==0) then
      !
      call readkeyboad(cmd=cmd,casename=casename)
      !
      if(cmd=='list' .or. cmd=='help') then
        call listcmd
      endif
      !
      print*,' ** input command: ',cmd,casename
      !
    endif
    !
    call bcast(cmd)
    !
  end subroutine getcmd
  !+-------------------------------------------------------------------+
  !| The end of the subroutine getcmd.                                 |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to list all command that is predefined.        |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 28-05-2021  | Created by J. Fang @ Warrington                     |
  !+-------------------------------------------------------------------+
  subroutine listcmd
    !
    use parallel, only: mpirank
    !
    if(mpirank==0) then
      !
      write(*,*)' +------------------------------------------------------------+'
      write(*,*)' |                          command line                      |'
      write(*,*)' +----------------+-------------------------------------------+'
      write(*,*)' |        command |                                  function |'
      write(*,*)' +----------------+-------------------------------------------+'
      write(*,*)' |      -cmd list |               to list all functionalities |'
      write(*,*)' +----------------+-------------------------------------------+'
      write(*,*)' |      -cmd help |                          the same as list |'
      write(*,*)' +----------------+-------------------------------------------+'
      write(*,*)' |       -cmd run |                      to run a computation |'
      write(*,*)' +----------------+-------------------------------------------+'
      write(*,*)' |      -cmd init |             to generation an example case |'
      write(*,*)' +----------------+-------------------------------------------+'
      write(*,*)' |     -input *** |                  to assign the input file |'
      write(*,*)' +----------------+-------------------------------------------+'
      !
    endif
    !
  end subroutine listcmd
  !+-------------------------------------------------------------------+
  !| The end of the subroutine listcmd.                                |
  !+-------------------------------------------------------------------+
  !
end module cmdefne
!+---------------------------------------------------------------------+
!| The end of the module  cmdefne                                      |
!+---------------------------------------------------------------------+