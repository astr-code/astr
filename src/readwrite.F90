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
  implicit none
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
      write(*,*)
      print*,' +------------------------ Statement -------------------------+'
      print*,' |                                                            |'
      print*,' |              Developed by Jian Fang since 2008             |'
      print*,' |                   Re-mastered at 02-2021                   |'
      print*,' |               <ASTR> Copyright Resvered <2008>             |'
      print*,' |         Advanced Simulatior for Turbulence Research        |'
      print*,' |                                                            |'
      print*,' +------------------------------------------------------------+'
      print*,' |                                                            |'
      print*,' | Copyright 2008 Jian Fang                                   |'
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
      write(*,*)
      !
      write(*,*)'                        ___   _____________  '
      write(*,*)'                       / _ | / __/_  __/ _ \ '
      write(*,*)'                      / __ |_\ \  / / / , _/ '
      write(*,*)'                     /_/ |_/___/ /_/ /_/|_|  '
      write(*,*)
      write(*,*)
      !
      print*,' ** computation start ... '
      write(*,*)
      write(*,*)
    endif
    !
    return
    !
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
    use parallel,only : mpirank,mpirankmax,isize,jsize,ksize,mpistop
    use commvar, only : ia,ja,ka,hm,numq,conschm,difschm,nondimen,     &
                        diffterm,ref_t,reynolds,mach,num_species,      &
                        flowtype,ndims,lfilter,alfa_filter,bctype,     &
                        twall,lfftk,kcutoff,ninit,rkscheme
    !
    ! local data
    character(len=42) :: typedefine
    integer :: n,i
    character(len=4) :: bcdir(1:6)
    !
    bcdir(1)='imin'; bcdir(2)='imax'
    bcdir(3)='jmin'; bcdir(4)='jmax'
    bcdir(5)='kmin'; bcdir(6)='kmax'
    !
    if(mpirank==mpirankmax) then
      !
      select case(trim(flowtype))
      case('2dvort')
        typedefine='                 2D inviscid vortical flow'
      case('channel')
        typedefine='                              channel flow'
      case('tgv')
        typedefine='                  Taylor-Green Vortex flow'
      case('jet')
        typedefine='                                  Jet flow'
      case('accutest')
        typedefine='                             accuracy test'
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
      write(*,'(2X,A18,A2,A42)')trim(flowtype),': ',typedefine
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
      elseif(conschm(4:4)=='e') then
        write(*,'(A)',advance='no')' explicit '
      else
        stop ' !! error: conschm not defined !!'
      endif
      write(*,'(11A)')conschm(3:3),'-',conschm(2:2),'-',               &
                      conschm(1:1),'......',conschm(1:1),'-',          &
                      conschm(2:2),'-',conschm(3:3)
      write(*,'(19X,A)',advance='no')'   diffusion terms: '
      if(difschm(4:4)=='c') then
        write(*,'(A)',advance='no')' compact '
      elseif(difschm(4:4)=='e') then
        write(*,'(A)',advance='no')' explicit '
      else
        stop ' !! error: difschm not defined !!'
      endif
      write(*,'(11A)')difschm(3:3),'-',difschm(2:2),'-',               &
                      difschm(1:1),'......',difschm(1:1),'-',          &
                      difschm(2:2),'-',difschm(3:3)
      if(lfilter) then
        write(*,'(2X,A43)',advance='no')' low-pass filter is used, '
        write(*,'(A12,F7.3)')' coefficient:',alfa_filter
      endif
      if(lfftk) then
        write(*,'(2X,A41)',advance='no')' FFT used at the k direction,'
        write(*,'(A,I4,A,I4)')'  cutoff k:',kcutoff,' /',ka/2
        !
        if(kcutoff>=ka/2) then
          write(*,'(2X,A62)')'   !! Warning  cutoff wavenumber too big !!'
        endif
        !
      endif
      !
      write(*,'(14X,A)',advance='no')'  temporal scheme: '
      if(rkscheme=='rk3') then
        write(*,'(A)')'  3-stage 3rd-order Runge–Kutta'
      elseif(rkscheme=='rk4') then
        write(*,'(A)')'  4-stage 4th-order Runge–Kutta'
      else
        stop ' !! error: rk scheme not defined !!'
      endif
      !
      write(*,'(2X,62A)')('-',i=1,62)
      write(*,'(2X,A)')'                    *** Boundary Conditions ***'
      do n=1,6
        !
        if(bctype(n)==1) then
          write(*,'(45X,I0,2(A))')bctype(n),' periodic at: ',bcdir(n)
        elseif(bctype(n)==41) then
          write(*,'(24X,I0,3(A),F6.3)')bctype(n),' isothermal wall at ',bcdir(n),&
                                                     ' Twall= ',twall(n)
        elseif(bctype(n)==11) then
          write(*,'(45X,I0,2(A))')bctype(n),' inflow  at: ',bcdir(n)
        elseif(bctype(n)==21) then
          write(*,'(45X,I0,2(A))')bctype(n),' outflow at: ',bcdir(n)
        elseif(bctype(n)==22) then
          write(*,'(38X,I0,2(A))')bctype(n),' nscbcc outflow at: ',bcdir(n)
        else
          print*,n,bctype(n)
          stop ' !! BC not defined !!'
        endif
        !
      enddo
      write(*,'(2X,62A)')('-',i=1,62)
      !
      if(ninit==2) then
        write(*,'(2X,A62)')' initialise by rading a flowini2d file'
      elseif(ninit==3) then
        write(*,'(2X,A62)')' initialise by rading a flowini3d file'
      else
        write(*,'(2X,A62)')' initialise by a user defined way'
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
                        lreadgrid,lfftk,gridfile,bctype,twall,kcutoff, &
                        ninit,rkscheme
    use parallel,only : bcast
    !
    ! local data
    character(len=64) :: inputfile
    integer :: n
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
      read(11,*)nondimen,diffterm,lfilter,lreadgrid,lfftk
      read(11,'(/)')
      read(11,*)alfa_filter,kcutoff
      !
      if(nondimen) then
        read(11,'(/)')
        read(11,*)ref_t,reynolds,mach
      else
        read(11,'(///)')
      endif
      !
      read(11,'(/)')
      read(11,*)conschm,difschm,rkscheme
      read(11,'(/)')
      read(11,*)num_species
      read(11,'(/)')
      do n=1,6
        read(11,*)bctype(n)
        if(bctype(n)==41) then
          backspace(11)
          read(11,*)bctype(n),twall(n)
        endif
      enddo
      read(11,'(/)')
      read(11,*)ninit
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
    call bcast(lfftk)
    !
    call bcast(flowtype)
    call bcast(conschm)
    call bcast(difschm)
    call bcast(gridfile)
    !
    call bcast(nondimen)
    call bcast(diffterm)
    call bcast(rkscheme)
    call bcast(lfilter)
    !
    call bcast(alfa_filter)
    call bcast(kcutoff)
    !
    call bcast(ref_t)
    call bcast(reynolds)
    call bcast(mach)
    !
    call bcast(num_species)
    !
    call bcast(bctype)
    call bcast(twall)
    call bcast(ninit)
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
    use commvar, only: maxstep,nwrite,deltat,nlstep,lwrite,lavg,navg
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
      read(11,*)lwrite,lavg
      read(11,'(/)')
      read(11,*)maxstep,nwrite,nlstep,navg
      read(11,'(/)')
      read(11,*)deltat
      close(11)
      print*,' >> ',trim(inputfile),' ... done'
      !
    endif
    !
    call bcast(lwrite)
    call bcast(lavg)
    call bcast(maxstep)
    call bcast(nwrite)
    call bcast(nlstep)
    call bcast(navg)
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
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! This subroutine is used to read a initial flow filed.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine readflowini3d
    !
    use commvar,   only : im,jm,km,num_species
    use commarray, only : rho,vel,prs,tmp,spc
    use hdf5io
    !
    ! local data
    !
    integer :: jsp
    character(len=2) :: spname
    !
    call h5io_init(filename='datin/flowini3d.h5',mode='read')
    !
    call h5read(varname='ro',var=rho(0:im,0:jm,0:km))
    call h5read(varname='u1', var=vel(0:im,0:jm,0:km,1))
    call h5read(varname='u2', var=vel(0:im,0:jm,0:km,2))
    call h5read(varname='u3', var=vel(0:im,0:jm,0:km,3))
    call h5read(varname='p', var=prs(0:im,0:jm,0:km))
    call h5read(varname='t', var=tmp(0:im,0:jm,0:km))
    do jsp=1,num_species
       write(spname,'(i2.2)')jsp
      call h5read(varname='sp'//spname,var=spc(0:im,0:jm,0:km,jsp))
    enddo
    !
    call h5io_end
    !
  end subroutine readflowini3d
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! End of subroutine readflowini3d.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
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
    use commvar, only: time,nstep,filenumb,num_species,im,jm,km,       &
                       lwrite,lavg,nsamples
    use commarray,only : x,rho,vel,prs,tmp,spc,qrhs,rom,u1m,u2m,u3m,   &
                         pm,tm,u11,u22,u33,u12,u13,u23,tt,pp
    !
    use hdf5io
    !
    ! local data
    character(len=4) :: stepname
    character(len=2) :: spname
    character(len=64) :: outfilename
    integer :: jsp
    !
    if(lwrite) then
      write(stepname,'(i4.4)')filenumb
      filenumb=filenumb+1
      !
      outfilename='outdat/flowfield'//stepname//'.h5'
    else
      outfilename='outdat/flowfield.h5'
    endif
    !
    call h5io_init(trim(outfilename),mode='write')
    call h5write(varname='nstep',var=nstep)
    call h5write(varname='time',var=time)
    call h5write(varname='ro',var=rho(0:im,0:jm,0:km))
    call h5write(varname='u1', var=vel(0:im,0:jm,0:km,1))
    call h5write(varname='u2', var=vel(0:im,0:jm,0:km,2))
    call h5write(varname='u3', var=vel(0:im,0:jm,0:km,3))
    call h5write(varname='p', var=prs(0:im,0:jm,0:km))
    call h5write(varname='t', var=tmp(0:im,0:jm,0:km))
    do jsp=1,num_species
       write(spname,'(i2.2)')jsp
      call h5write(varname='sp'//spname,var=spc(0:im,0:jm,0:km,jsp))
    enddo
    call h5io_end
    !
    if(lavg .and. nsamples>0) then
      !
      call h5io_init('outdat/meanflow.h5',mode='write')
      call h5write(varname='nstep',var=nstep)
      call h5write(varname='nsamples',var=nsamples)
      call h5write(varname='rom',var=rom(0:im,0:jm,0:km))
      call h5write(varname='u1m',var=u1m(0:im,0:jm,0:km))
      call h5write(varname='u2m',var=u2m(0:im,0:jm,0:km))
      call h5write(varname='u3m',var=u3m(0:im,0:jm,0:km))
      call h5write(varname='pm ',var=pm (0:im,0:jm,0:km))
      call h5write(varname='tm ',var=tm (0:im,0:jm,0:km))
      call h5write(varname='u11',var=u11(0:im,0:jm,0:km))
      call h5write(varname='u22',var=u22(0:im,0:jm,0:km))
      call h5write(varname='u33',var=u33(0:im,0:jm,0:km))
      call h5write(varname='u12',var=u12(0:im,0:jm,0:km))
      call h5write(varname='u13',var=u13(0:im,0:jm,0:km))
      call h5write(varname='u23',var=u23(0:im,0:jm,0:km))
      call h5write(varname='tt ',var=tt (0:im,0:jm,0:km))
      call h5write(varname='pp ',var=pp (0:im,0:jm,0:km))
      call h5io_end
      !
    endif
    !
    call tecbin('testout/tecfield'//mpirankname//'.plt',               &
                                              x(0:im,0:jm,0:km,1),'x', &
                                              x(0:im,0:jm,0:km,2),'y', &
                                              x(0:im,0:jm,0:km,3),'z', &
                                            rho(0:im,0:jm,0:km),'ro',  &
                                            vel(0:im,0:jm,0:km,1),'u', &
                                            vel(0:im,0:jm,0:km,2),'v', &
                                            prs(0:im,0:jm,0:km),'p',   &
                                            spc(0:im,0:jm,0:km,1),'Y1' )
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
