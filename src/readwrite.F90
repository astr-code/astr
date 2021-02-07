!+---------------------------------------------------------------------+
!| This module contains subroutines of reading and writing files.      |
!| ==============                                                      |
!| CHANGE RECORD                                                       |
!| -------------                                                       |
!| 06-Oct-2018  | Created by J. Fang @ Warrington                      |
!+---------------------------------------------------------------------+
module readwrite
  !
  use parallel,only : mpirank,mpirankname
  use tecio
  !
  contains
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to print welcome infomation.              |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 08-Oct-2018  | Created by J. Fang @ Warrington                    |
  !+-------------------------------------------------------------------+
  subroutine statement
    !
    if(mpirank==0) then
      !
      write(*,*)
      write(*,*)'                ___   _____________            ___ '
      write(*,*)'               / _ | / __/_  __/ _ \   ____   / _\\'
      write(*,*)'              / __ |_\ \  / / / , _/  /___/  / , _/'
      write(*,*)'             /_/ |_/___/ /_/ /_/|_|         /_/|_| '
      write(*,*)
      write(*,*)
  
      print*,' +----------------------- Statement --------------------------+'
      print*,' |                                                            |'
      print*,' |                   Developed by Jian Fang                   |'
      print*,' |             <ASTR-R> Copyright Resvered <2020>             |'
      print*,' |                 ASTR code for reacting flow                |'
      print*,' |                                                            |'
      print*,' +------------------------------------------------------------+'
      print*,' |                                                            |'
      print*,' | Copyright 2020 Jian Fang                                   |'
      print*,' |                                                            |'
      print*,' | Licensed under the Apache License, Version 2.0             |'
      print*,' | you may not use this file except in compliance with the    |'
      print*,' | license. You may obtain a copy of the License at           |'
      print*,' |                                                            |'
      print*,' |        http://www.apache.org/licenses/LICENSE-2.0          |'
      print*,' |                                                            |'
      print*,' | Unless required by applicable law or agreed to in writing, |'
      print*,' | software distributed under the License is distributed on an|'
      print*,' | "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY     |'
      print*,' | KIND, either express or implied. See the License for the   |'
      print*,' | specific language governing permissions and limitations    |'
      print*,' | under the License.                                         |'
      print*,' |                                                            |'
      print*,' +--------------------- End of Statement ---------------------+'
    !
    endif
  end subroutine statement
  !+-------------------------------------------------------------------+
  !| The end of the subroutine statement.                              |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to print the state of computation         |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 07-02-2021  | Created by J. Fang @ Warrington                     |
  !+-------------------------------------------------------------------+
  subroutine infodisp
    !
    use parallel, only : lio,isize,jsize,ksize
    use commvar, only : ia,ja,ka
    !
    if(lio) then
      write(*,'(2X,62A)')('-',i=1,62)
      write(*,'(2X,A)')'                   *** computation Setup ***'
      write(*,'(2X,62A)')('-',i=1,62)
      write(*,'(2X,A)')'                     ***    MPI size    ***'
      write(*,"(4(1x,A15))")'isize','jsize','ksize','size'
      write(*,"(4(1x,I15))")isize,jsize,ksize,isize*jsize*ksize
      write(*,'(2X,62A)')('-',i=1,62)
      write(*,'(2X,A)')'                        *** Grid Size ***'
      write(*,"(1x,3(A21))")'im','jm','km'
      write(*,"(1x,3(I21))")ia,ja,ka
      write(*,'(2X,62A)')('-',i=1,62)
    endif
    !
  end subroutine infodisp
  !+-------------------------------------------------------------------+
  !| The end of the subroutine infodisp.                               |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to read input file of defining and        |
  !|  controlling a simulation                                         |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 06-02-2021  | Created by J. Fang @ Warrington                     |
  !+-------------------------------------------------------------------+
  subroutine readinput
    !
    use commvar, only : ia,ja,ka,lihomo,ljhomo,lkhomo
    !
    ! local data
    character(len=64) :: inputfile
    !
    inputfile='datin/input.dat'
    !
    if(mpirank==0) then
      !
      open(11,file=trim(inputfile),action='read')
      read(11,'(///)')
      read(11,*)ia,ja,ka
      write(*,'(3(a,i0))')'  ** Dimension: ',ia,' x ',ja,' x ',ka
      read(11,'(/)')
      read(11,*)lihomo,ljhomo,lkhomo
      write(*,'(A)',advance='no')'  ** homogeneous direction: '
      if(lihomo) write(*,'(A)',advance='no')'i,'
      if(ljhomo) write(*,'(A)',advance='no')' j,'
      if(lkhomo) write(*,'(A)')' k'
      close(11)
      print*,' >> ',trim(inputfile),' ... done'
      !
    endif
    !
  end subroutine readinput
  !+-------------------------------------------------------------------+
  !| The end of the subroutine readinput.                              |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to print the state of computation         |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 07-02-2021  | Created by J. Fang @ Warrington                     |
  !+-------------------------------------------------------------------+
  subroutine readgrid
    !
    use commvar,   only : im,jm,km
    use commarray, only : x
    use hdf5io
    !
    call h5io_init(filename='datin/grid.h5',mode='read')
    call h5read(varname='x',var=x(0:im,0:jm,0:km,1))
    call h5read(varname='y',var=x(0:im,0:jm,0:km,2))
    call h5read(varname='z',var=x(0:im,0:jm,0:km,3))
    call h5io_end
    !
    call tecbin('tecgrid'//mpirankname//'.plt',x(0:im,0:jm,0:km,1),'x',&
                                               x(0:im,0:jm,0:km,2),'y',&
                                               x(0:im,0:jm,0:km,3),'z' )
    !
  end subroutine readgrid
  !+-------------------------------------------------------------------+
  !| The end of the subroutine readinput.                              |
  !+-------------------------------------------------------------------+
  !
  !
end module readwrite
!+---------------------------------------------------------------------+
!| The end of the module readwrite.                                    |
!+---------------------------------------------------------------------+
