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
    use parallel,only : lio,isize,jsize,ksize,mpistop
    use commvar, only : ia,ja,ka,hm,numq,conschm,difschm,nondimen,     &
                        diffterm,ref_t,reynolds,mach,num_species,      &
                        flowtype,ndims
    !
    ! local data
    character(len=42) :: typedefine
    !
    if(lio) then
      !
      select case(trim(flowtype))
      case('2dvort')
        typedefine='                 2D inviscid vortical flow'
      case default
        print*,trim(flowtype)
        stop ' !! flowtype not defined @ infodisp'
      end select
      !
      write(*,'(2X,62A)')('-',i=1,62)
      write(*,'(2X,A)')'                   *** computation Setup ***'
      write(*,'(2X,62A)')('-',i=1,62)
      write(*,'(2X,A)')'                     ***    flow type   ***'
      write(*,'(2X,A59,I3)')' dimension: ',ndims
      write(*,'(14X,A,A,A42)')trim(flowtype),': ',typedefine
      write(*,'(2X,62A)')('-',i=1,62)
      write(*,'(2X,A)')'                     ***    MPI size    ***'
      write(*,"(4(1x,A15))")'isize','jsize','ksize','size'
      write(*,"(4(1x,I15))")isize,jsize,ksize,isize*jsize*ksize
      write(*,'(2X,62A)')('-',i=1,62)
      write(*,'(2X,A)')'                        *** Grid Size ***'
      write(*,"(4(1x,A15))")'im','jm','km','halo'
      write(*,"(4(1x,I15))")ia,ja,ka,hm
      write(*,'(2X,62A)')('-',i=1,62)
      write(*,'(2X,A)')'                         *** Equations ***'
      !
      if(nondimen) then
        write(*,'(2X,A62)')' non-dimensional'
      else
        write(*,'(2X,A62)')' metric equations'
      endif
      !
      if(diffterm) then
        write(*,'(2X,A62)')' N-S equations are solved'
      else
        write(*,'(2X,A62)')' Euler equations are solved'
      endif
      !
      write(*,'(2X,A59,I3)')' number of independent variables: ',numq
      write(*,'(2X,A59,I3)')' number of species: ',num_species
      !
      write(*,'(2X,62A)')('-',i=1,62)
      write(*,'(2X,A)')'                      *** Flow Parameters ***'
      if(nondimen) then
         write(*,'(4X,3(A20))')'ref_t','Reynolds','Mach'
         write(*,"(4X,3(F20.6))")ref_t,reynolds,mach
      endif
      write(*,'(2X,62A)')('-',i=1,62)
      write(*,'(2X,A)')'                          *** sceheme ***'
      write(*,'(19X,A)',advance='no')'  convection terms: '
      if(conschm(4:4)=='c') then
        write(*,'(A)',advance='no')' compact '
        write(*,'(11A)')conschm(3:3),'-',conschm(2:2),'-',             &
                        conschm(1:1),'......',conschm(1:1),'-',        &
                        conschm(2:2),'-',conschm(3:3)
      else
        stop ' !! error: conschm not defined !!'
      endif
      write(*,'(19X,A)',advance='no')'   diffusion terms: '
      if(difschm(4:4)=='c') then
        write(*,'(A)',advance='no')' compact '
        write(*,'(11A)')difschm(3:3),'-',difschm(2:2),'-',             &
                        difschm(1:1),'......',difschm(1:1),'-',        &
                        difschm(2:2),'-',difschm(3:3)
      else
        stop ' !! error: difschm not defined !!'
      endif
      write(*,'(2X,62A)')('-',i=1,62)
      !
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
    use commvar, only : ia,ja,ka,lihomo,ljhomo,lkhomo,conschm,difschm, &
                        nondimen,diffterm,ref_t,reynolds,mach,         &
                        num_species,flowtype
    use parallel,only : bcast
    !
    ! local data
    character(len=64) :: inputfile
    !
    inputfile='datin/input.dat'
    !
    if(mpirank==0) then
      !
      open(11,file=trim(inputfile),action='read')
      read(11,'(////)')
      read(11,*)flowtype
      read(11,'(/)')
      read(11,*)ia,ja,ka
      read(11,'(/)')
      read(11,*)lihomo,ljhomo,lkhomo
      write(*,'(A)',advance='no')'  ** homogeneous direction: '
      if(lihomo) write(*,'(A)',advance='no')'i,'
      if(ljhomo) write(*,'(A)',advance='no')' j,'
      if(lkhomo) write(*,'(A)')' k'
      read(11,'(/)')
      read(11,*)nondimen,diffterm
      !
      if(nondimen) then
        read(11,'(/)')
        read(11,*)ref_t,reynolds,mach
      else
        read(11,'(///)')
      endif
      !
      read(11,'(/)')
      read(11,*)conschm,difschm
      read(11,'(/)')
      read(11,*)num_species
      close(11)
      print*,' >> ',trim(inputfile),' ... done'
      !
    endif
    !
    call bcast(ia)
    call bcast(ja)
    call bcast(ka)
    !
    call bcast(lihomo)
    call bcast(ljhomo)
    call bcast(lkhomo)
    !
    call bcast(flowtype)
    call bcast(conschm)
    call bcast(difschm)
    !
    call bcast(nondimen)
    call bcast(diffterm)
    call bcast(ref_t)
    call bcast(reynolds)
    call bcast(mach)
    !
    call bcast(num_species)
    !
  end subroutine readinput
  !+-------------------------------------------------------------------+
  !| The end of the subroutine readinput.                              |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to read a file that is sued to control    |
  !| computation.                                                      |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 09-02-2021  | Created by J. Fang @ Warrington                     |
  !+-------------------------------------------------------------------+
  subroutine readcont
    !
    use commvar, only: maxstep,nwrite,deltat
    use parallel,only : bcast
    !
    ! local data
    character(len=64) :: inputfile
    !
    inputfile='datin/contr.dat'
    !
    if(mpirank==0) then
      !
      open(11,file=trim(inputfile),action='read')
      read(11,'(////)')
      read(11,*)maxstep,nwrite
      read(11,'(/)')
      read(11,*)deltat
      close(11)
      print*,' >> ',trim(inputfile),' ... done'
      !
    endif
    !
    call bcast(maxstep)
    call bcast(nwrite)
    call bcast(deltat)
    !
  end subroutine readcont
  !+-------------------------------------------------------------------+
  !| The end of the subroutine readcont.                               |
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
