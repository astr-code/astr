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
  subroutine readkeyboad(keyin)
    !
    character(len=*),intent(out),optional :: keyin
    !
    ! local data
    integer :: ierr,cli_len,nlen,arg_count
    integer,save :: nkey=0
    !
    nkey=nkey+1
    call get_command_argument(nkey,keyin,cli_len,ierr)
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
      call readkeyboad(cmd)
      !
      if(cmd=='list' .or. cmd=='help') then
        call listcmd
      endif
      !
      print*,' ** input command: ',cmd
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
      write(*,*)' |    list / help |               to list all functionalities |'
      write(*,*)' +----------------+-------------------------------------------+'
      write(*,*)' |            run |                      to run a computation |'
      write(*,*)' +----------------+-------------------------------------------+'
      write(*,*)' |             pp |                          pre/post-process |'
      write(*,*)' |                | init              generate a example case |'
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