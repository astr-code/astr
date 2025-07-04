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
  !!+-------------------------------------------------------------------+
  !| Subroutine: listcmd                                               |
  !|                                                                   |
  !| Purpose:                                                          |
  !|   This subroutine prints out a list of all predefined             |
  !|   command-line options available in the ASTR code.                |
  !|   Only the root MPI rank (rank 0) performs the output.            |
  !|                                                                   |
  !| Commands supported:                                               |
  !|   - list / help : list all available functionalities              |
  !|   - run         : run a computation using a given input file      |
  !|                   (usage: mpirun -np 8 ./astr run datin/input)    |
  !|   - test        : run internal code tests                         |
  !|   - pp          : enter pre/post-processing mode                  |
  !|                                                                   |
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 28-05-2021  | Created by J. Fang @ Warrington                     |
  !| 03-07-2025  | Prettify by J. Fang @ IMech, Beijing, using ChatGPT |
  !+-------------------------------------------------------------------+
  subroutine listcmd
  
    use parallel, only : mpirank
  
    if (mpirank == 0) then
      write(*,*) ' +------------------------------------------------------------+'
      write(*,*) ' |                        Command Line Help                   |'
      write(*,*) ' +----------------+-------------------------------------------+'
      write(*,*) ' |     Command    | Description                               |'
      write(*,*) ' +----------------+-------------------------------------------+'
      write(*,*) ' | list / help    | List all functionalities                  |'
      write(*,*) ' | run            | Run a computation with an input file      |'
      write(*,*) ' |     usage: mpirun -np 8 ./astr run datin/input             |'
      write(*,*) ' | test           | Run code test routines                    |'
      write(*,*) ' |     usage: mpirun -np 8 ./astr test grad                   |'
      write(*,*) ' | pp             | Pre/Post-processing                       |'
      write(*,*) ' +----------------+-------------------------------------------+'
    end if
  
  end subroutine listcmd
  !+-------------------------------------------------------------------+
  !| End of subroutine listcmd                                         |
  !+-------------------------------------------------------------------+

end module cmdefne
!+---------------------------------------------------------------------+
!| The end of the module  cmdefne                                      |
!+---------------------------------------------------------------------+