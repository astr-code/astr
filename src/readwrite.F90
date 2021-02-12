!+---------------------------------------------------------------------+
!| This module contains subroutines of reading and writing files.      |
!| ==============                                                      |
!| CHANGE RECORD                                                       |
!| -------------                                                       |
!| 06-Oct-2018  | Created by J. Fang @ Warrington                      |
!+---------------------------------------------------------------------+
module readwrite
  !
  use parallel,only : mpirank,mpirankname,mpistop,lio
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
  !| This subroutine is used to initialise files and folders.          |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 12-02-2021  | Created by J. Fang @ Warrington                     |
  !+-------------------------------------------------------------------+
  subroutine fileini
    !
    use commvar, only : hand_fs,hand_rp
    !
    if(lio) then
      call system('mkdir testout/')
      call system('mkdir outdat/')
      !
      hand_fs=13
      hand_rp=14
      !
    endif
    !
  end subroutine fileini
  !+-------------------------------------------------------------------+
  !| The end of the subroutine fileini.                                |
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
                        flowtype,ndims,lfilter,alfa_filter
    !
    ! local data
    character(len=42) :: typedefine
    !
    if(lio) then
      !
      select case(trim(flowtype))
      case('2dvort')
        typedefine='                 2D inviscid vortical flow'
      case('channel')
        typedefine='                              channel flow'
      case('tgv')
        typedefine='                  Taylor-Green Vortex flow'
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
      !
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
      if(lfilter) then
        write(*,'(2X,A43)',advance='no')' low-pass filter is used, '
        write(*,'(A13,F6.3)')' coefficient:',alfa_filter
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
                        num_species,flowtype,lfilter,alfa_filter,      &
                        lreadgrid,gridfile
    use parallel,only : bcast
    !
    ! local data
    character(len=64) :: inputfile
    !
    if(mpirank==0) then
      !
      call readkeyboad(inputfile=inputfile)
      !
      if(trim(inputfile)=='.') inputfile='datin/input.dat'
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
      read(11,*)nondimen,diffterm,lfilter,lreadgrid
      read(11,'(/)')
      read(11,*)alfa_filter
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
      if(lreadgrid) then
        read(11,'(/)')
        read(11,'(A)')gridfile
      endif
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
    call bcast(lreadgrid)
    !
    call bcast(flowtype)
    call bcast(conschm)
    call bcast(difschm)
    call bcast(gridfile)
    !
    call bcast(nondimen)
    call bcast(diffterm)
    call bcast(lfilter)
    !
    call bcast(alfa_filter)
    !
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
    use commvar, only: maxstep,nwrite,deltat,nlstep
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
      read(11,*)maxstep,nwrite,nlstep
      read(11,'(/)')
      read(11,*)deltat
      close(11)
      print*,' >> ',trim(inputfile),' ... done'
      !
    endif
    !
    call bcast(maxstep)
    call bcast(nwrite)
    call bcast(nlstep)
    call bcast(deltat)
    !
  end subroutine readcont
  !+-------------------------------------------------------------------+
  !| The end of the subroutine readcont.                               |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to read grid file.                        |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 07-02-2021  | Created by J. Fang @ Warrington                     |
  !+-------------------------------------------------------------------+
  subroutine readgrid
    !
    use commvar,   only : im,jm,km,gridfile
    use commarray, only : x
    use hdf5io
    !
    call h5io_init(filename=trim(gridfile),mode='read')
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
  !+-------------------------------------------------------------------+
  !| This subroutine is used to write a grid file for postprocess.     |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 07-02-2021  | Created by J. Fang @ Warrington                     |
  !+-------------------------------------------------------------------+
  subroutine writegrid
    !
    use commvar,   only : im,jm,km
    use commarray, only : x
    use hdf5io
    !
    call h5io_init('./grid.h5',mode='write')
    call h5write(varname='x',var=x(0:im,0:jm,0:km,1))
    call h5write(varname='y',var=x(0:im,0:jm,0:km,2))
    call h5write(varname='z',var=x(0:im,0:jm,0:km,3))
    call h5io_end
    !
  end subroutine writegrid
  !+-------------------------------------------------------------------+
  !| The end of the subroutine readinput.                              |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to input command from keyboard.           |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 11-02-2021  | Created by J. Fang @ Warrington                     |
  !+-------------------------------------------------------------------+
  subroutine readkeyboad(inputfile)
    !
    character(len=*),intent(out),optional :: inputfile
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
    enddo
    !
  end subroutine readkeyboad
  !+-------------------------------------------------------------------+
  !| The end of the subroutine readkeyboad.                            |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to dump flow field data                   |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 12-02-2021  | Created by J. Fang @ Warrington                     |
  !+-------------------------------------------------------------------+
  subroutine output
    !
    use commvar, only: time,nstep,filenumb,num_species
    use commarray,only : x,rho,vel,prs,tmp,spc,qrhs
    use hdf5io
    !
    ! local data
    character(len=4) :: stepname
    character(len=2) :: spname
    integer :: jsp
    !
    write(stepname,'(i4.4)')filenumb
    !

    
    call h5io_init('outdat/flowfield'//stepname//'.h5',mode='write')
    call h5write(varname='nstep',var=nstep)
    call h5write(varname='time',var=time)
    call h5write(varname='ro',var=rho(0:im,0:jm,0:km))
    call h5write(varname='u', var=vel(0:im,0:jm,0:km,1))
    call h5write(varname='v', var=vel(0:im,0:jm,0:km,2))
    call h5write(varname='w', var=vel(0:im,0:jm,0:km,3))
    call h5write(varname='p', var=prs(0:im,0:jm,0:km))
    call h5write(varname='t', var=tmp(0:im,0:jm,0:km))
    do jsp=1,num_species
       write(spname,'(i2.2)')jsp
      call h5write(varname='sp'//spname,var=spc(0:im,0:jm,0:km,jsp))
    enddo
    call h5io_end
    !
    ! call tecbin('testout/tecfield'//stepname//mpirankname//'.plt',     &
    !                                           x(0:im,0:jm,0:km,1),'x', &
    !                                           x(0:im,0:jm,0:km,2),'y', &
    !                                           x(0:im,0:jm,0:km,3),'z', &
    !                                         rho(0:im,0:jm,0:km),'ro',  &
    !                                         vel(0:im,0:jm,0:km,1),'u', &
    !                                         vel(0:im,0:jm,0:km,2),'v', &
    !                                         prs(0:im,0:jm,0:km),'p',   &
    !                                         spc(0:im,0:jm,0:km,1),'Y1' )
    filenumb=filenumb+1
    !
  end subroutine output
  !+-------------------------------------------------------------------+
  !| The end of the subroutine output.                                 |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to write report file.                     |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 12-02-2021  | Created by J. Fang @ Warrington                     |
  !+-------------------------------------------------------------------+
  subroutine timerept
    !
    use commvar, only : hand_rp,nstep,maxstep,ctime
    !
    ! local data
    logical,save :: linit=.true.
    logical :: lexist
    integer :: i
    !
    if(lio) then
      !
      if(linit) then
        !
        inquire(file='report.txt', exist=lexist)
        !
        if(lexist) call system('mv -v report.txt report.bak')
        !
        open(hand_rp,file='report.txt')
        linit=.false.
      endif
      !
      write(hand_rp,'(2X,62A)')('-',i=1,62)
      write(hand_rp,'(2X,A,I7)')'time report at nstep ',nstep
      write(hand_rp,'(2X,62A)')('-',i=1,62)
      write(hand_rp,'(2X,A,E13.6E2,A,F6.2,A)')'total time cost : ',    &
                            ctime(2),' - ',100.d0*ctime(2)/ctime(2),' %'
      !
#ifdef cputime
      write(hand_rp,'(2X,A,E13.6E2,A,F6.2,A)')'  - rk          : ',    &
                            ctime(3),' - ',100.d0*ctime(3)/ctime(2),' %'
      write(hand_rp,'(2X,A,E13.6E2,A,F6.2,A)')'    - rhs       : ',    &
                            ctime(4),' - ',100.d0*ctime(4)/ctime(2),' %'
      write(hand_rp,'(2X,A,E13.6E2,A,F6.2,A)')'      - Convc   : ',    &
                            ctime(9),' - ',100.d0*ctime(9)/ctime(2),' %'
      write(hand_rp,'(2X,A,E13.6E2,A,F6.2,A)')'      - Diffu   : ',    &
                          ctime(10),' - ',100.d0*ctime(10)/ctime(2),' %'
      write(hand_rp,'(2X,A,E13.6E2,A,F6.2,A)')'    - filter    : ',    &
                            ctime(8),' - ',100.d0*ctime(9)/ctime(2),' %'
      write(hand_rp,'(2X,A,E13.6E2,A,F6.2,A)')'    - io        : ',    &
                            ctime(6),' - ',100.d0*ctime(6)/ctime(2),' %'
      write(hand_rp,'(2X,A,E13.6E2,A,F6.2,A)')'    - sta       : ',    &
                            ctime(5),' - ',100.d0*ctime(5)/ctime(2),' %'
      write(hand_rp,'(2X,A,E13.6E2,A,F6.2,A)')'  - com         : ',    &
                           ctime(7), ' - ',100.d0*ctime(7)/ctime(2),' %'
#endif
      !
      flush(hand_rp)
      !
      if(nstep==maxstep) then
        close(hand_rp)
        print*,' << report.txt'
      endif
      !
    endif
    !
  end subroutine timerept
  !+-------------------------------------------------------------------+
  !| The end of the subroutine timerept.                               |
  !+-------------------------------------------------------------------+
  !
  !
end module readwrite
!+---------------------------------------------------------------------+
!| The end of the module readwrite.                                    |
!+---------------------------------------------------------------------+
