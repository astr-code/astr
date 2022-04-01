!+---------------------------------------------------------------------+
!| This module contains subroutines of reading and writing files.      |
!| ==============                                                      |
!| CHANGE RECORD                                                       |
!| -------------                                                       |
!| 06-Oct-2018  | Created by J. Fang @ Warrington                      |
!+---------------------------------------------------------------------+
module readwrite
  !
  use constdef
  use parallel,only : mpirank,mpirankname,mpistop,lio,irk,jrkm,jrk,    &
                      ptime,bcast,mpirankmax,ig0,jg0,kg0
  use commvar, only : ndims,im,jm,km
  use tecio
  use stlaio,  only: get_unit
  !
  implicit none
  !
  integer :: nxtchkpt,nxtwsequ
  !+---------------------+---------------------------------------------+
  !|            nxtchkpt | the next nstep to checkpoint.               |
  !|            nxtwsequ | the next nstep to write flow field.         |
  !+---------------------+---------------------------------------------+
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
    if(lio) then
      !
      call system('mkdir testout/')
      call system('mkdir outdat/')
      call system('mkdir bakup/')
      !
      call system('mkdir monitor/')
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
    use parallel,only : mpirank,isize,jsize,ksize,mpistop
    use commvar, only : ia,ja,ka,hm,numq,conschm,difschm,nondimen,     &
                        diffterm,ref_t,reynolds,mach,num_species,      &
                        flowtype,ndims,lfilter,alfa_filter,            &
                        lfftk,kcutoff,ninit,rkscheme,                  &
                        spg_imin,spg_imax,spg_jmin,spg_jmax,           &
                        spg_kmin,spg_kmax,lchardecomp,recon_schem,     &
                        lrestart,limmbou,solidfile,bfacmpld,           &
                        turbmode,iomode,schmidt
    use bc,      only : bctype,twall,xslip,turbinf,xrhjump,angshk
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
      case('cylinder')
        typedefine='                      flow past a cylinder'
      case('mixlayer')
        typedefine='                              mixing layer'
      case('shuosher')
        typedefine='                         shu-osher problem'
      case('bl')
        typedefine='                       boundary layer flow'
      case('swbli')
        typedefine='     shock-wave/boundary layer interaction'
      case('windtunn')
        typedefine='                   a numerical wind tunnel'
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
      if(num_species>0) then
      write(*,'(1X,A,12(1X,F8.6))')' schmidt number: ',(schmidt(i),i=1,num_species)
      endif
      !
      write(*,'(2X,62A)')('-',i=1,62)
      write(*,'(2X,A)')'                      *** Flow Parameters ***'
      if(nondimen) then
         write(*,'(4X,3(A20))')'ref_t','Reynolds','Mach'
         write(*,"(4X,3(F20.6))")ref_t,reynolds,mach
      endif
      write(*,'(2X,62A)')('-',i=1,62)
      write(*,'(2X,A)')'                          *** sceheme ***'
      write(*,'(11X,A)',advance='no')'  convection terms: '
      if(conschm(4:4)=='c') then
        write(*,'(A)',advance='no')' compact '
      elseif(conschm(4:4)=='e') then
        write(*,'(A)',advance='no')'explicit '
      else
        stop ' !! error: conschm not defined, only c & e allowed !!'
      endif
      !
      read(conschm(1:1),*) n
      !
      if(mod(n,2)==0) then
        write(*,'(A)',advance='no')'central '
      else
        write(*,'(A)',advance='no')' upwind '
      endif

      !
      write(*,'(11A)')conschm(3:3),'-',conschm(2:2),'-',               &
                      conschm(1:1),'......',conschm(1:1),'-',          &
                      conschm(2:2),'-',conschm(3:3)
      if(mod(n,2).ne.0) then
        if(lchardecomp) then
          write(*,'(25X,39A)')'local characteristic decomposition used'
        else
          write(*,'(28X,39A)')'    reconstruction in physical space'
        endif
        !
        write(*,'(15X,A)',advance='no')' reconstruction method : '
        !
        if(conschm(4:4)=='e') then
          select case(recon_schem)
          case(-1)
            write(*,'(24A)')'  1st-order construction'
          case(0)
            write(*,'(24A)')'     linear construction'
          case(1)
            write(*,'(24A)')'       WENO construction'
          case(2)
            write(*,'(24A)')'     WENO-Z construction'
          case(3)
            write(*,'(24A)')'         MP construction'
          case(4)
            write(*,'(24A)')'   WENO-SYM construction'
          case(5)
            write(*,'(24A)')'      MP-LD construction'
            write(*,'(A56,F8.3)')'            b factor = ',bfacmpld
          case(6)
            write(*,'(24A)')'      round construction'
          case default
            stop ' !! reconstruction scheme not defined !!' 
          end select
        elseif(conschm(4:4)=='c') then
            write(*,'(24A)')'      MP-LD construction'
            write(*,'(A56,F8.3)')'            b factor = ',bfacmpld
        endif
        !
      endif
      !
      write(*,'(11X,A)',advance='no')'   diffusion terms: '
      if(difschm(4:4)=='c') then
        write(*,'(A)',advance='no')' compact central '
      elseif(difschm(4:4)=='e') then
        write(*,'(A)',advance='no')'explicit central '
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
          write(*,'(37X,I0,2(A))')bctype(n),' isothermal wall at: ',bcdir(n)
          write(*,'(33X,A,F12.6)')' wall temperature: ',twall(n)
        elseif(bctype(n)==411) then
          write(*,'(23X,I0,2(A))')bctype(n),' slip-nonslip isothermal wall at: ',bcdir(n)
          write(*,'(31X,A,F12.6)')' slip-nonslip point: ',xslip
          write(*,'(33X,A,F12.6)')' wall temperature: ',twall(n)
        elseif(bctype(n)==421) then
          write(*,'(23X,I0,2(A))')bctype(n),'  slip-nonslip adiabatic wall at: ',bcdir(n)
          write(*,'(31X,A,F12.6)')' slip-nonslip point: ',xslip
        elseif(bctype(n)==51) then
          write(*,'(44X,I0,2(A))')bctype(n),' farfield at: ',bcdir(n)
          !
          if(trim(flowtype)=='swbli') then
            write(*,'(30X,A,F12.6)')'      angle of shock: ',angshk
            write(*,'(30X,A,F12.6)')' R-H jump point at x: ',xrhjump
          endif
          !
        elseif(bctype(n)==52) then
          write(*,'(38X,I0,2(A))')bctype(n),' nscbc farfield at: ',bcdir(n)
        elseif(bctype(n)==11) then
          write(*,'(45X,I0,2(A))')bctype(n),' inflow  at: ',bcdir(n)
          write(*,'(18X,(A))',advance='no')'inflow turbulence : '
          !
          if(turbinf=='none') then
            write(*,'(A)')'      no inflow turbulence'
          elseif(turbinf=='intp') then
            write(*,'(A)')' interpolation of database'
          else
            print*,' !! turbinf: ',turbinf
            stop ' !! ERROR in defining turbinf @ bc 11 !!'
          endif
          !
        elseif(bctype(n)==12) then
          write(*,'(38X,I0,2(A))')bctype(n),' nscbc inflow  at: ',bcdir(n)
        elseif(bctype(n)==13) then
          write(*,'(44X,I0,2(A))')bctype(n),' fixed bc at: ',bcdir(n)
        elseif(bctype(n)==21) then
          write(*,'(45X,I0,2(A))')bctype(n),' outflow at: ',bcdir(n)
        elseif(bctype(n)==22) then
          write(*,'(39X,I0,2(A))')bctype(n),' nscbc outflow at: ',bcdir(n)
        elseif(bctype(n)==23) then
          write(*,'(36X,I0,2(A))')bctype(n),' gc-nscbc outflow at: ',bcdir(n)
        else
          print*,n,bctype(n)
          stop ' !! BC not defined !!'
        endif
        !
      enddo
      write(*,'(2X,62A)')('-',i=1,62)
      !
      if(limmbou) then
        write(*,'(2X,A)')'                    *** Solid Body Immersed ***'
        write(*,'(2X,A,A)')'solid body file: ',trim(solidfile)
        write(*,'(2X,62A)')('-',i=1,62)
      endif
      !

      write(*,'(2X,A)')'                      *** Turbulence model ***'
      !
      if(trim(turbmode)=='none') then
        write(*,'(35X,A)')'          no turbulence model'
      elseif(trim(turbmode)=='k-omega') then
        write(*,'(35X,A)')' the Menter SST k-omega model'
      elseif(trim(turbmode)=='udf1') then
        write(*,'(35X,A)')'         user defined model 1'
      else
        print*,' !! ERROR in defining turbulence model: ',turbmode
        stop
      endif
      !
      write(*,'(2X,62A)')('-',i=1,62)
      !
      write(*,'(2X,A)')'                          *** i/o model ***'
      if(iomode=='s') then
        write(*,'(35X,A)')'                stream output'
      elseif(iomode=='h') then
        write(*,'(35X,A)')'            structured output'
      else
        print*,' !! ERROR in defining iomode: ',iomode
        stop
      endif
      !
      write(*,'(2X,62A)')('-',i=1,62)
      !
      write(*,'(2X,A)')'                 *** Initialising method model ***'
      if(lrestart) then
        write(*,'(2X,A62)')' restart a previous computation'
      else
        if(ninit==1) then
          write(*,'(2X,A62)')'          initialise by using a profile'
        elseif(ninit==2) then
          write(*,'(2X,A62)')' initialise by reading a flowini2d file'
        elseif(ninit==3) then
          write(*,'(2X,A62)')' initialise by reading a flowini3d file'
        else
          write(*,'(2X,A62)')' initialise by a user defined way'
        endif
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
                        lreadgrid,lfftk,gridfile,kcutoff,              &
                        ninit,rkscheme,spg_imin,spg_imax,spg_jmin,     &
                        spg_jmax,spg_kmin,spg_kmax,lchardecomp,        &
                        recon_schem,lrestart,limmbou,solidfile,        &
                        bfacmpld,shkcrt,turbmode,iomode,schmidt
    use parallel,only : bcast
    use cmdefne, only : readkeyboad
    use bc,      only : bctype,twall,xslip,turbinf,xrhjump,angshk
    !
    ! local data
    character(len=64) :: inputfile
    character(len=5) :: char
    integer :: n,fh,i
    !
    if(mpirank==0) then
      !
      call readkeyboad(inputfile)
      !
      if(trim(inputfile)=='.') inputfile='datin/input.dat'
      !
      fh=get_unit()
      !
      open(fh,file=trim(inputfile),action='read')
      read(fh,'(////)')
      read(fh,*)flowtype
      read(fh,'(/)')
      read(fh,*)ia,ja,ka
      read(fh,'(/)')
      read(fh,*)lihomo,ljhomo,lkhomo
      write(*,'(A)',advance='no')'  ** homogeneous direction: '
      if(lihomo) write(*,'(A)',advance='no')'i,'
      if(ljhomo) write(*,'(A)',advance='no')' j,'
      if(lkhomo) write(*,'(A)')' k'
      read(fh,'(/)')
      read(fh,*)nondimen,diffterm,lfilter,lreadgrid,lfftk,limmbou
      read(fh,'(/)')
      read(fh,*)lrestart
      read(fh,'(/)')
      read(fh,*)alfa_filter,kcutoff
      !
      if(nondimen) then
        read(fh,'(/)')
        read(fh,*)ref_t,reynolds,mach
      else
        read(fh,'(///)')
      endif
      !
      read(fh,'(/)')
      read(fh,*)conschm,difschm,rkscheme
      read(fh,'(/)')
      read(fh,*)recon_schem,lchardecomp,bfacmpld,shkcrt
      read(fh,'(/)')
      read(fh,*)num_species
      if(num_species>0) then
        allocate(schmidt(1:num_species))
        backspace(fh)
        read(fh,*)num_species,(schmidt(i),i=1,num_species)
      endif
      read(fh,'(/)')
      read(fh,*)turbmode,iomode
      read(fh,'(/)')
      do n=1,6
        read(fh,*)bctype(n)
        if(bctype(n)==41) then
          backspace(fh)
          read(fh,*)bctype(n),twall(n)
        endif
        if(bctype(n)==411) then
          backspace(fh)
          read(fh,*)bctype(n),xslip,twall(n)
        endif
        if(bctype(n)==421) then
          backspace(fh)
          read(fh,*)bctype(n),xslip
        endif
        if(bctype(n)==11) then
          backspace(fh)
          read(fh,*)bctype(n),turbinf
        endif
        if(bctype(n)==51) then
          if(trim(flowtype)=='swbli') then
            backspace(fh)
            read(fh,*)bctype(n),xrhjump,angshk
          endif
        endif
      enddo
      read(fh,'(/)')
      read(fh,*)ninit
      read(fh,'(/)')
      read(fh,*)spg_imin,spg_imax,spg_jmin,spg_jmax,spg_kmin,spg_kmax
      read(fh,'(/)')
      read(fh,'(A)')gridfile
      if(limmbou) then
        read(fh,'(/)')
        read(fh,'(A)')solidfile
      endif
      close(fh)
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
    call bcast(limmbou)
    call bcast(lchardecomp)
    !
    call bcast(lrestart)
    !
    call bcast(flowtype)
    call bcast(turbmode)
    call bcast(iomode)
    !
    call bcast(conschm)
    call bcast(difschm)
    call bcast(gridfile)
    call bcast(solidfile)
    !
    call bcast(nondimen)
    call bcast(diffterm)
    call bcast(rkscheme)
    call bcast(lfilter)
    !
    call bcast(alfa_filter)
    call bcast(bfacmpld)
    call bcast(shkcrt)
    call bcast(kcutoff)
    !
    call bcast(ref_t)
    call bcast(reynolds)
    call bcast(mach)
    !
    call bcast(recon_schem)
    call bcast(num_species)
    call bcast(schmidt,num_species)
    !
    call bcast(bctype)
    call bcast(twall)
    call bcast(xslip)
    call bcast(ninit)
    !
    call bcast(spg_imin)
    call bcast(spg_imax)
    call bcast(spg_jmin)
    call bcast(spg_jmax)
    call bcast(spg_kmin)
    call bcast(spg_kmax)
    !
    call bcast(turbinf)
    call bcast(xrhjump)
    call bcast(angshk)
    !
    call readslic
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
    use commvar, only: deltat,lwsequ,lwslic,lavg,maxstep,feqchkpt,     &
                       feqwsequ,feqslice,feqlist,feqavg
    use parallel,only: bcast
    !
    ! local data
    character(len=64) :: inputfile
    integer :: fh
    !
    inputfile='datin/controller'
    !
    if(mpirank==0) then
      !
      fh=get_unit()
      !
      open(fh,file=trim(inputfile),action='read')
      read(fh,'(////)')
      read(fh,*)lwsequ,lwslic,lavg
      read(fh,'(/)')
      read(fh,*)maxstep,feqchkpt,feqwsequ,feqslice,feqlist,feqavg
      read(fh,'(/)')
      read(fh,*)deltat
      close(fh)
      print*,' >> ',trim(inputfile),' ... done'
      !
    endif
    !
    call bcast(lwsequ)
    call bcast(lwslic)
    call bcast(lavg)
    call bcast(maxstep)
    call bcast(feqchkpt)
    call bcast(feqwsequ)
    call bcast(feqslice)
    call bcast(feqlist)
    call bcast(feqavg)
    call bcast(deltat)
    !
  end subroutine readcont
  !+-------------------------------------------------------------------+
  !| The end of the subroutine readcont.                               |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to read a file that defines monitors.     |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 13-07-2021  | Created by J. Fang @ Warrington                     |
  !+-------------------------------------------------------------------+
  subroutine readslic
    !
    use commvar, only : islice,jslice,kslice,ia,ja,ka
    use parallel,only : bcast
    !
    ! local data
    character(len=64) :: inputfile
    logical :: lexist
    integer :: fh
    !
    if(mpirank==0) then
      !
      inputfile='datin/slice.dat'
      inquire(file=trim(inputfile), exist=lexist)
      if(lexist) then
        !
        fh=get_unit()
        open(fh,file=trim(inputfile),action='read')
        read(fh,*)
        read(fh,*)islice,jslice,kslice
        close(fh)
        print*,' >> ',trim(inputfile),' ... done'
        !
      else
        islice=-1
        jslice=-1
        kslice=-1
      endif
      !
      if(islice>=0 .and. islice<=ia) then
        write(*,'(A,I0,A)')'  ** silce cut at i: ',islice,' is to save'
      endif
      if(jslice>=0 .and. jslice<=ja) then
        write(*,'(A,I0,A)')'  ** silce cut at j: ',jslice,' is to save'
      endif
      if(kslice>=0 .and. kslice<=ka) then
        write(*,'(A,I0,A)')'  ** silce cut at k: ',kslice,' is to save'
      endif
      !
    endif
    !
    call bcast(islice)
    call bcast(jslice)
    call bcast(kslice)
    !
  end subroutine readslic
  !
  subroutine readmonc
    !
    use commvar, only : nmonitor,imon,islice,jslice,kslice,ia,ja,ka
    use parallel,only : bcast
    use commcal, only : monitorsearch
    !
    ! local data
    character(len=64) :: inputfile
    character(len=32) :: c1,c2,c3
    logical :: lexist
    integer :: ios,bignum,ncou,i,n,fh
    integer,allocatable :: ijk(:,:)
    real(8),allocatable :: xyz(:,:)
    !
    bignum=1024
    allocate(ijk(1:bignum,4),xyz(1:bignum,4))
    !
    ijk=-1
    !
    !
    if(mpirank==0) then
      !
      inputfile='datin/monitor.dat'
      inquire(file=trim(inputfile), exist=lexist)
      !
      ncou=0
      !
      if(lexist) then
        !
        fh=get_unit()
        !
        open(fh,file=trim(inputfile),action='read')
        read(fh,*,iostat=ios)
        if(ios.ne.0) then
          nmonitor=0
        endif
        write(*,'(2X,62A)')('-',i=1,62)
        write(*,'(2X,A)')'                     *** monitor points ***'
        do while(ios==0)
          !
          read(fh,*,iostat=ios)c1,c2,c3
          if(ios==0) then
            !
            ncou=ncou+1
            !
            if(isnum(c1)==1) then
              !
              read(c1,*)ijk(ncou,1)
              read(c2,*)ijk(ncou,2)
              read(c3,*)ijk(ncou,3)
              ijk(ncou,4)=ncou
              !
              write(*,'(24X,I3,A,3(I12))')ijk(ncou,4),':',ijk(ncou,1:3)
              !
            elseif(isnum(c1)==2 .or. isnum(c1)==3 .or. isnum(c1)==4) then
              !
              read(c1,*)xyz(ncou,1)
              read(c2,*)xyz(ncou,2)
              read(c3,*)xyz(ncou,3)
              xyz(ncou,4)=dble(ncou)
              !
              write(*,'(24X,I3,A,3(F12.6))')ncou,':',xyz(ncou,:)
            else
              print*,'ncou=',ncou
              stop ' !! ERROR 1 @ readmon'
            endif
            !
          else
            exit
          endif
          !
        enddo
        write(*,'(2X,62A)')('-',i=1,62)
        close(fh)
        print*,' >> ',trim(inputfile),' ... done'
        !
      endif
      !
      inputfile='datin/slice.dat'
      inquire(file=trim(inputfile), exist=lexist)
      if(lexist) then
        !
        fh=get_unit()
        open(fh,file=trim(inputfile),action='read')
        read(fh,*)
        read(fh,*)islice,jslice,kslice
        close(fh)
        print*,' >> ',trim(inputfile),' ... done'
        !
      else
        islice=-1
        jslice=-1
        kslice=-1
      endif
      !
      if(islice>=0 .and. islice<=ia) then
        write(*,'(A,I0,A)')'  ** silce cut at i: ',islice,' is to save'
      endif
      if(jslice>=0 .and. jslice<=ja) then
        write(*,'(A,I0,A)')'  ** silce cut at j: ',jslice,' is to save'
      endif
      if(kslice>=0 .and. kslice<=ka) then
        write(*,'(A,I0,A)')'  ** silce cut at k: ',kslice,' is to save'
      endif
      !
    endif
    !
    call bcast(ncou)
    call bcast(islice)
    call bcast(jslice)
    call bcast(kslice)
    !
    if(ncou>0) then
      !
      call bcast(ijk(1:ncou,:))
      call bcast(xyz(1:ncou,:))
      !
      call monitorsearch(ijk(1:ncou,:),xyz(1:ncou,:),imon,nmonitor)
      !
      ! if(nmonitor>0) then
      !   do n=1,nmonitor
      !     print*,mpirank,'|',n,':',imon(n,:)
      !   enddo
      ! endif
      !
    else
      nmonitor=0
    endif
    !
    ! call mpistop
    !
  end subroutine readmonc
  !+-------------------------------------------------------------------+
  !| The end of the subroutine readmonc.                               |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to output monitor files.                  |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 07-02-2021  | Created by J. Fang @ Warrington                     |
  !+-------------------------------------------------------------------+
  subroutine writemon
    !
    use commvar, only: nmonitor,imon,nstep,time,pinf,deltat
    use commarray, only : x,rho,vel,prs,tmp,dvel
    !
    ! local data
    integer :: n,i,j,k,ios,ns
    integer,allocatable,save :: fh(:),record(:)
    logical,save :: firstcall = .true.
    logical :: lexist
    character(len=4) :: moname
    character(len=64) :: filename
    character(len=32) :: c1
    real(8) :: rtime
    !
    if(nmonitor>0) then
      !
      if(firstcall) then
        !
        allocate(fh(nmonitor),record(1:nmonitor))
        !
        record=0
        !
        do n=1,nmonitor
          !
          i=imon(n,1)
          j=imon(n,2)
          k=imon(n,3)
          !
          write(moname,'(i4.4)')imon(n,4)
          filename='monitor/monitor'//moname//'.dat'
          !
          fh(n)=get_unit()
          !
          inquire(file=trim(filename), exist=lexist)
          open(fh(n),file=trim(filename),access='direct',recl=8*17)
          !
          if(nstep==0 .or. (.not.lexist)) then
            ! create new monitor files
            ! write(fh(n),'(A,3(1X,I0))')'#',imon(n,1)+ig0,              &
            !                                imon(n,2)+jg0,imon(n,3)+kg0
            ! write(fh(n),'(A,3(1X,E15.7E3))')'#',x(i,j,k,1:3)
            ! write(fh(n),"(A7,1X,A13,15(1X,A20))")'nstep','time',        &
            !             'u','v','w','ro','p','t','dudx','dudy','dudz',  &
            !                 'dvdx','dvdy','dvdz','dwdx','dwdy','dwdz'
            ! write(*,'(3(A),I0)')'   ** create monitor file ',          &
            !                       trim(filename),' using handle: ',fh(n)
          else
            !
            ! ios=0
            ! do while(ios==0)
            !   read(fh(n),*,iostat=ios)c1
            !   !
            !   if(IsNum(c1)==1) then
            !     read(c1,*)ns
            !     !
            !     if(ns==nstep) exit
            !     !
            !   endif
            !   !
            ! enddo
            ! if(ios<0) backspace(fh(n))
            !
            ns=0
            do while(ns<nstep-1)
              record(n)=record(n)+1
              !
              read(fh(n),rec=record(n))ns,rtime
              ! if(n==3) print*,mpirank,'|',record(n),ns,rtime
            enddo
            !
            write(*,*)'   ** resume monitor file ',          &
                                  trim(filename),' at nstep: ',ns
            !
          endif
          !
        enddo
        !
        firstcall=.false.
        !
      endif
      !
      do n=1,nmonitor
        i=imon(n,1)
        j=imon(n,2)
        k=imon(n,3)
        ! write(fh(n),"(I7,1X,E13.6E2,15(1X,E20.13E2))")nstep,time,     &
        !     vel(i,j,k,1:3),rho(i,j,k),prs(i,j,k)/pinf,tmp(i,j,k),     &
        !     dvel(i,j,k,1,:),dvel(i,j,k,2,:),dvel(i,j,k,3,:)
        record(n)=record(n)+1
        write(fh(n),rec=record(n))nstep,time,vel(i,j,k,:),rho(i,j,k),prs(i,j,k), &
                               tmp(i,j,k),dvel(i,j,k,:,:)
        ! write(*,*)nstep,time,vel(i,j,k,:),rho(i,j,k),prs(i,j,k), &
        !                        tmp(i,j,k),dvel(i,j,k,:,:)
      enddo
      !
    else
      return
    endif
    !
    return
    !
  end subroutine writemon
  !+-------------------------------------------------------------------+
  !| The end of the subroutine writemon.                               |
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
    return
    !
  end subroutine readgrid
  !+-------------------------------------------------------------------+
  !| The end of the subroutine readinput.                              |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to read STL file that defines a solid body!
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 01-07-2021  | Created by J. Fang @ Warrington                     |
  !+-------------------------------------------------------------------+
  subroutine readsolid(inputfile)
    !
    use commtype,  only : solid,triangle
    use commvar,   only : immbody,nsolid
    use stlaio
    !
    ! arguments
    character(len=*),intent(in) :: inputfile
    !
    ! local data
    logical :: lstl
    integer :: solid_num,node_num,face_num,text_num,jso
    !
    print*,' >> solid body: '
    !
    nsolid=num_solid(trim(inputfile))
    !
    allocate(immbody(1:nsolid))
    !
    write(*,'(A,I0)')'    ** number of solids : ',nsolid
    !
    do jso=1,nsolid
      immbody(jso)%num_face=num_face_in_solid(trim(inputfile),jso)
      call immbody(jso)%alloface()
      write(*,'(2(A,I0))')'    ** number of faces for solid ',jso, &
                                        ' : ',immbody(jso)%num_face
    enddo
    !
    call stla_read(trim(inputfile),immbody)
    !
    ! print*,' face num for solde 1: ',num_face_in_solid(trim(solidfile),1)
    ! print*,' face num for solde 2: ',num_face_in_solid(trim(solidfile),2)
    ! call stla_size(trim(solidfile),nsolid,node_num,face_num,text_num)
    ! call stla_size_print(trim(solidfile),nsolid,node_num,face_num,text_num)
    ! print*,'      file: ',trim(solidfile)
    ! print*,' solid num: ',solid_num
    ! print*,'  node num: ',node_num
    ! print*,'  face num: ',face_num
    ! print*,'  text num: ',text_num
    !
    return
    !
  end subroutine readsolid
  !+-------------------------------------------------------------------+
  !| The end of the subroutine readinput.                              |
  !+-------------------------------------------------------------------+
  !!
  !+-------------------------------------------------------------------+
  !| This subroutine is used to read a profile                         !
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 13-08-2021  | Created by J. Fang @ Warrington                     |
  !+-------------------------------------------------------------------+
  subroutine readprofile(filename,dir,var1,var2,var3,var4,var5,skipline)
    !
    use commvar,   only : im,jm,km,ia,ja,ka
    use parallel,  only : ptabupd,pscatter
    !
    ! arguments
    character(len=*),intent(in) :: filename,dir
    integer,intent(in),optional :: skipline
    real(8),optional,intent(out) ::  var1(0:),var2(0:),var3(0:),       &
                                     var4(0:),var5(0:)
    !
    ! local data
    integer :: fh,ios,i,j,n,offset,nsk,nvar
    real(8),allocatable,dimension(:,:) :: data,scat
    integer :: table(0:mpirankmax),displs(0:mpirankmax)
    character(len=32) :: c1
    !
    if(present(skipline)) then
      nsk=skipline
    else
      nsk=0
    endif
    !
    fh=get_unit()
    !
    if( present(var5) ) then
      nvar=5
    elseif( present(var4)  ) then
      nvar=4
    elseif( present(var3) ) then
      nvar=3
    elseif( present(var2) ) then
      nvar=2
    elseif( present(var1) ) then
      nvar=1
    else
      stop ' !! ERROR 1  define nvar @  readprofile '
    endif
    !
    allocate(data(1:nvar,0:ja))
    !
    if(dir=='j') then
      !
      if(mpirank==0) then
        !
        open(fh,file=filename,action='read',form='formatted')
        do n=1,nsk
          read(fh,*)
        enddo
        !
        ios=0
        do while(ios==0)
          !
          read(fh,*,iostat=ios)c1
          !
          if(isnum(c1)==2 .or. isnum(c1)==3 .or. isnum(c1)==4) then
            ! the first line with a real number
            backspace(fh)
            !
            exit
            !
          endif
          !
        enddo
        !
        ! start from real number
        do j=0,ja
          read(fh,*)(data(i,j),i=1,nvar)
        enddo
        !
        close(fh)
        write(*,'(3A,I0)')'  >> ',filename,' using file handle : ',fh
        !
      endif
      !
      call ptabupd(var=jm+1,table=table)
      !
      call ptabupd(var=jg0*nvar,table=displs)
      !
      call pscatter(data,scat,table,offset=displs,rank=0)
      !
      ! print*,mpirank,'|',displs(mpirank)
      ! do j=0,jm
      !   print*,mpirank,'|',j+jg0,scat(1,j+1)
      ! enddo
      !
      if( present(var5) )  var5(0:jm)=scat(5,1:1+jm)
      if( present(var4) )  var4(0:jm)=scat(4,1:1+jm)
      if( present(var3) )  var3(0:jm)=scat(3,1:1+jm)
      if( present(var2) )  var2(0:jm)=scat(2,1:1+jm)
      if( present(var1) )  var1(0:jm)=scat(1,1:1+jm)
      !
      deallocate(data,scat)
      !
    else
      !
      print*,' !! dir not set yet @ readprofile !!'
      stop
      !
    endif
    !
    return
    !
  end subroutine readprofile
  !+-------------------------------------------------------------------+
  !| The end of the subroutine readprofile.                            |
  !+-------------------------------------------------------------------+
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! This subroutine is used to read a initial flow filed.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine readflowini3d
    !
    use commvar,   only : im,jm,km,num_species,time,roinf,tinf,pinf
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
    ! call h5read(varname='time',var=time)
    call h5read(varname='ro',var=rho(0:im,0:jm,0:km))
    call h5read(varname='u1', var=vel(0:im,0:jm,0:km,1))
    call h5read(varname='u2', var=vel(0:im,0:jm,0:km,2))
    call h5read(varname='u3', var=vel(0:im,0:jm,0:km,3))
    call h5read(varname='p', var=prs(0:im,0:jm,0:km))
    call h5read(varname='t', var=tmp(0:im,0:jm,0:km))
    !
    ! rho=roinf
    ! prs=pinf
    ! tmp=tinf
    ! do jsp=1,num_species
    !    write(spname,'(i2.2)')jsp
    !   call h5read(varname='sp'//spname,var=spc(0:im,0:jm,0:km,jsp))
    ! enddo
    !
    call h5io_end
    !
  end subroutine readflowini3d
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! End of subroutine readflowini3d.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! This subroutine is used to read a initial flow filed.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine readflowini2d
    !
    use commvar,   only : im,jm,km,num_species
    use commarray, only : rho,vel,prs,tmp,spc
    use fludyna,   only : thermal
    use hdf5io
    !
    ! local data
    !
    integer :: jsp,i,j,k
    character(len=2) :: spname
    !
    call h5io_init(filename='datin/flowini2d.h5',mode='read')
    !
    call h5read(varname='ro', var=rho(0:im,0:jm,0),  dir='k')
    call h5read(varname='u1', var=vel(0:im,0:jm,0,1),dir='k')
    call h5read(varname='u2', var=vel(0:im,0:jm,0,2),dir='k')
    call h5read(varname='t',  var=tmp(0:im,0:jm,0),  dir='k')
    ! do jsp=1,num_species
    !    write(spname,'(i2.2)')jsp
    !   call h5read(varname='sp'//spname,var=spc(0:im,0:jm,0:km,jsp))
    ! enddo
    !
    call h5io_end
    !
    do k=0,km
    do j=0,jm
    do i=0,im
      !
      rho(i,j,k)  =rho(i,j,0)  
      vel(i,j,k,1)=vel(i,j,0,1)
      vel(i,j,k,2)=vel(i,j,0,2)
      vel(i,j,k,3)=0.d0
      tmp(i,j,k)  =tmp(i,j,0)
      prs(i,j,k)  =thermal(temperature=tmp(i,j,k),density=rho(i,j,k)) 
      !
      do jsp=1,num_species
        spc(i,j,k,jsp)=0.d0 !spc(i,j,0,jsp)
      enddo
      !
    enddo
    enddo
    enddo
    !
  end subroutine readflowini2d
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! End of subroutine readflowini2d.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !+-------------------------------------------------------------------+
  !| This subroutine is used to read a checkpoint file for restarting. |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 24-06-2021  | Created by J. Fang @ Warrington                     |
  !+-------------------------------------------------------------------+
  subroutine readcheckpoint(folder)
    !
    use commvar, only: nstep,filenumb,fnumslic,time,flowtype,          &
                       num_species,im,jm,km,force,numq,turbmode,lwsequ
    use commarray, only : rho,vel,prs,tmp,spc,q,tke,omg,miut
    use statistic, only : massflux,massflux_target,nsamples
    use hdf5io
    use bc,        only : ninflowslice
    !
    ! arguments
    character(len=*),intent(in) :: folder
    !
    ! local data
    integer :: nstep_1,jsp
    character(len=2) :: spname,qname
    character(len=128) :: infilename
    character(len=4) :: stepname
    !
    call h5io_init(filename=folder//'/auxiliary.h5',mode='read')
    call h5read(varname='nstep',var=nstep)
    call h5read(varname='filenumb',var=filenumb)
    call h5read(varname='fnumslic',var=fnumslic)
    call h5read(varname='ninflowslice',var=ninflowslice)
    if(flowtype=='channel') then
      call h5read(varname='massflux',var=massflux)
      call h5read(varname='massflux_target',var=massflux_target)
      call h5read(varname='force',var=force,dim=3)
    endif
    call h5read(varname='nsamples',var=nsamples)
    call h5io_end
    !
    if(lwsequ) then
      write(stepname,'(i4.4)')filenumb
      infilename='outdat/flowfield'//stepname//'.h5'
    else
      infilename=folder//'/flowfield.h5'
    endif
    !
    call h5io_init(filename=trim(infilename),mode='read')
    call h5read(varname='nstep',var=nstep_1)
    !
    if(nstep_1==nstep) then
      call h5read(varname='time',var=time)
      call h5read(varname='ro',  var=rho(0:im,0:jm,0:km))
      call h5read(varname='u1',  var=vel(0:im,0:jm,0:km,1))
      call h5read(varname='u2',  var=vel(0:im,0:jm,0:km,2))
      call h5read(varname='u3',  var=vel(0:im,0:jm,0:km,3))
      call h5read(varname='p',   var=prs(0:im,0:jm,0:km))
      call h5read(varname='t',   var=tmp(0:im,0:jm,0:km))
      do jsp=1,num_species
         write(spname,'(i2.2)')jsp
        call h5read(varname='sp'//spname,var=spc(0:im,0:jm,0:km,jsp))
      enddo
      !
      if(trim(turbmode)=='k-omega') then
        call h5read(varname='k',     var=tke(0:im,0:jm,0:km))
        call h5read(varname='omega', var=omg(0:im,0:jm,0:km))
        call h5read(varname='miut',  var=miut(0:im,0:jm,0:km))
      elseif(trim(turbmode)=='udf1') then
        call h5read(varname='miut',  var=miut(0:im,0:jm,0:km))
      endif
      !
      ! call h5io_init('checkpoint/qdata.h5',mode='read')
      ! do jsp=1,numq
      !   write(qname,'(i2.2)')jsp
      !   call h5read(varname='q'//qname,var=q(0:im,0:jm,0:km,jsp))
      ! enddo
      !
      call h5io_end
      !
    else
      if(lio)  print*,' !! flowfield.h5 NOT consistent with auxiliary.h5'
      if(lio)  print*,' nstep =',nstep,' in auxiliary.h5 '
      if(lio)  print*,' nstep =',nstep_1,' in flowfield.h5 '
      call mpistop
    endif
    !
    ! infilename='checkpoint/flowfield'//mpirankname
    ! open(16,file=trim(infilename),form='unformatted')
    ! read(16)nstep,time
    ! read(16)rho,vel,prs,tmp
    ! read(16)spc
    ! close(16)
    ! if(lio) print*,' >> ',trim(infilename)
    ! !
    ! infilename='checkpoint/qdata'//mpirankname
    ! open(16,file=trim(infilename),form='unformatted')
    ! read(16)nstep
    ! read(16)q
    ! close(16)
    ! if(lio) print*,' >> ',trim(infilename)
    ! !
    ! infilename='checkpoint/auxiliary'//mpirankname
    ! open(16,file=trim(infilename),form='unformatted')
    ! read(16)nstep,filenumb
    ! read(16)massflux,massflux_target
    ! read(16)force
    ! close(16)
    ! if(lio) print*,' >> ',trim(infilename)
    !
    if(lio)  print*,' ** checkpoint file read. '
    !
  end subroutine readcheckpoint
  !+-------------------------------------------------------------------+
  !| The end of the subroutine readcheckpoint.                         |
  !+-------------------------------------------------------------------+
  !!
  !+-------------------------------------------------------------------+
  !| This subroutine is used to read mean flow statistics.             |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 30-06-2021  | Created by J. Fang @ Warrington                     |
  !+-------------------------------------------------------------------+
  subroutine readmeanflow
    !
    use commvar, only: nstep,im,jm,km
    use statistic,only : nsamples,liosta,nstep_sbeg,time_sbeg,         &
                         rom,u1m,u2m,u3m,pm,tm,                        &
                         u11,u22,u33,u12,u13,u23,pp,tt,tu1,tu2,tu3,    &
                         u111,u222,u333,u112,u113,                     &
                         u122,u133,u223,u233,u123,                     &
                         u1rem,u2rem,u3rem,pu1,pu2,pu3,                &
                         sgmam11,sgmam22,sgmam33,                      &
                         sgmam12,sgmam13,sgmam23,                      &
                         disspa,predil,visdif1,visdif2,visdif3,        &
                         allomeanflow,lmeanallocated
    use hdf5io
    !
    ! local data
    integer :: nstep_1,nsamples_1
    !
    call allomeanflow
    lmeanallocated=.true.
    !
    call h5io_init('outdat/meanflow.h5',mode='read')
    call h5read(varname='nstep',     var=nstep_1)
    if(nstep_1 .ne. nstep) then
      print*,' !! meanflow.h5 is NOT consistent with flowfield.h5'
      print*,' !! nstep:',nstep,'nstep_1:',nstep_1
      stop 
    endif
    call h5read(varname='nsamples',  var=nsamples_1)
    if(nsamples_1 .ne. nsamples) then
      print*,' !! meanflow.h5 is NOT consistent with auxiliary.h5'
      print*,' !! nsamples:',nsamples,'nsamples_1:',nsamples_1
      stop 
    endif
    call h5read(varname='nstep_sbeg',var=nstep_sbeg)
    call h5read(varname='time_sbeg', var=time_sbeg)
    call h5read(varname='rom',var=rom(0:im,0:jm,0:km))
    call h5read(varname='u1m',var=u1m(0:im,0:jm,0:km))
    call h5read(varname='u2m',var=u2m(0:im,0:jm,0:km))
    call h5read(varname='u3m',var=u3m(0:im,0:jm,0:km))
    call h5read(varname='pm ',var=pm (0:im,0:jm,0:km))
    call h5read(varname='tm ',var=tm (0:im,0:jm,0:km))
    call h5io_end
    !
    call h5io_init('outdat/2ndsta.h5',mode='read')
    call h5read(varname='nstep',var=nstep_1)
    if(nstep_1 .ne. nstep) then
      print*,' !! 2ndsta.h5 is NOT consistent with flowfield.h5'
      stop 
    endif
    call h5read(varname='u11',var=u11(0:im,0:jm,0:km))
    call h5read(varname='u22',var=u22(0:im,0:jm,0:km))
    call h5read(varname='u33',var=u33(0:im,0:jm,0:km))
    call h5read(varname='u12',var=u12(0:im,0:jm,0:km))
    call h5read(varname='u13',var=u13(0:im,0:jm,0:km))
    call h5read(varname='u23',var=u23(0:im,0:jm,0:km))
    !
    call h5read(varname='pp', var=pp(0:im,0:jm,0:km))
    call h5read(varname='tt', var=tt(0:im,0:jm,0:km))
    call h5read(varname='tu1',var=tu1(0:im,0:jm,0:km))
    call h5read(varname='tu2',var=tu2(0:im,0:jm,0:km))
    call h5read(varname='tu3',var=tu2(0:im,0:jm,0:km))
    call h5io_end
    !
    call h5io_init('outdat/3rdsta.h5',mode='read')
    call h5read(varname='nstep',     var=nstep_1)
    if(nstep_1 .ne. nstep) then
      print*,' !! 3rdsta.h5 is NOT consistent with flowfield.h5'
      stop 
    endif
    call h5read(varname='u111',var=u111(0:im,0:jm,0:km))
    call h5read(varname='u222',var=u222(0:im,0:jm,0:km))
    call h5read(varname='u333',var=u333(0:im,0:jm,0:km))
    call h5read(varname='u112',var=u112(0:im,0:jm,0:km))
    call h5read(varname='u113',var=u113(0:im,0:jm,0:km))
    call h5read(varname='u122',var=u122(0:im,0:jm,0:km))
    call h5read(varname='u133',var=u133(0:im,0:jm,0:km))
    call h5read(varname='u223',var=u223(0:im,0:jm,0:km))
    call h5read(varname='u233',var=u233(0:im,0:jm,0:km))
    call h5read(varname='u123',var=u123(0:im,0:jm,0:km))
    call h5io_end
    !
    call h5io_init('outdat/budget.h5',mode='read')
    call h5read(varname='nstep',   var=nstep_1)
    if(nstep_1 .ne. nstep) then
      print*,' !! budget.h5 is NOT consistent with flowfield.h5'
      stop 
    endif
    call h5read(varname='disspa', var=disspa(0:im,0:jm,0:km))
    call h5read(varname='predil', var=predil(0:im,0:jm,0:km))
    call h5read(varname='pu1',    var=pu1(0:im,0:jm,0:km))
    call h5read(varname='pu2',    var=pu2(0:im,0:jm,0:km))
    call h5read(varname='pu3',    var=pu3(0:im,0:jm,0:km))
    call h5read(varname='u1rem',  var=u1rem(0:im,0:jm,0:km))
    call h5read(varname='u2rem',  var=u2rem(0:im,0:jm,0:km))
    call h5read(varname='u3rem',  var=u3rem(0:im,0:jm,0:km))
    call h5read(varname='visdif1',var=visdif1(0:im,0:jm,0:km))
    call h5read(varname='visdif2',var=visdif2(0:im,0:jm,0:km))
    call h5read(varname='visdif3',var=visdif3(0:im,0:jm,0:km))
    call h5read(varname='sgmam11',var=sgmam11(0:im,0:jm,0:km))
    call h5read(varname='sgmam22',var=sgmam22(0:im,0:jm,0:km))
    call h5read(varname='sgmam33',var=sgmam33(0:im,0:jm,0:km))
    call h5read(varname='sgmam12',var=sgmam12(0:im,0:jm,0:km))
    call h5read(varname='sgmam13',var=sgmam13(0:im,0:jm,0:km))
    call h5read(varname='sgmam23',var=sgmam23(0:im,0:jm,0:km))
    call h5io_end
    !
  end subroutine readmeanflow
  !+-------------------------------------------------------------------+
  !| The end of the subroutine readmeanflow.                           |
  !+-------------------------------------------------------------------+
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
  !| This subroutine is used to bakup files.                           |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 29-09-2021  | Created by J. Fang @ Warrington                     |
  !+-------------------------------------------------------------------+
  subroutine bakupfile(file)
    !
    ! arguments
    character(len=*),intent(in),optional :: file
    !
    ! local data
    logical :: lexist
    !
    inquire(file=file, exist=lexist)
    !
    if(lexist) call system('mv -v '//file//' bakup/')
    !
  end subroutine bakupfile
  !+-------------------------------------------------------------------+
  !| The end of the subroutine bakupfile.                              |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to write sequence of flow field.          |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 28-12-2021  | Created by J. Fang @ Warrington                     |
  !+-------------------------------------------------------------------+
  subroutine writeflfed(subtime)
    !
    use commvar, only: time,nstep,filenumb,fnumslic,num_species,im,jm, &
                       km,lwsequ,turbmode,feqwsequ,iomode
    use commarray,only : x,rho,vel,prs,tmp,spc,q,ssf,lshock,crinod
    use models,   only : tke,omg,miut
    use hdf5io
    !
    ! arguments
    real(8),intent(inout),optional :: subtime
    !
    ! local data
    integer :: i,j,k,jsp
    character(len=4) :: stepname
    character(len=64) :: outfilename
    character(len=2) :: spname
    real(8),allocatable :: rshock(:,:,:),rcrinod(:,:,:)
    !
    if(lwsequ .and. nstep==nxtwsequ) then
      !
      filenumb=filenumb+1
      !
      write(stepname,'(i4.4)')filenumb
      !
      outfilename='outdat/flowfield'//stepname//'.'//iomode//'5'
    else
      stepname=''
      outfilename='outdat/flowfield.'//iomode//'5'
    endif
    !
    ! outfilename='outdat/flowfield.h5'
    call h5io_init(trim(outfilename),mode='write')
    if(iomode=='s') then
      call h5write(varname='ro',var=rho(0:im,0:jm,0:km),char='s')
      call h5write(varname='u1',var=vel(0:im,0:jm,0:km,1),char='s')
      call h5write(varname='u2',var=vel(0:im,0:jm,0:km,2),char='s')
      call h5write(varname='u3',var=vel(0:im,0:jm,0:km,3),char='s')
      call h5write(varname='p', var=prs(0:im,0:jm,0:km),char='s')
      call h5write(varname='t', var=tmp(0:im,0:jm,0:km),char='s')
    elseif(iomode=='h') then
      call h5write(varname='ro',var=rho(0:im,0:jm,0:km))
      call h5write(varname='u1',var=vel(0:im,0:jm,0:km,1))
      call h5write(varname='u2',var=vel(0:im,0:jm,0:km,2))
      call h5write(varname='u3',var=vel(0:im,0:jm,0:km,3))
      call h5write(varname='p', var=prs(0:im,0:jm,0:km))
      call h5write(varname='t', var=tmp(0:im,0:jm,0:km))
      if(num_species>0) then
        do jsp=1,num_species
           write(spname,'(i2.2)')jsp
          call h5write(varname='sp'//spname,var=spc(0:im,0:jm,0:km,jsp))
        enddo
      endif
    else
      print*,' !! iomode error: ',iomode
      stop
    endif
    !
    ! if(allocated(ssf)) then
    !   call h5write(varname='ssf', var=ssf(0:im,0:jm,0:km))
    ! endif
    ! if(allocated(lshock)) then
    !   allocate(rshock(0:im,0:jm,0:km))
    !   allocate(rcrinod(0:im,0:jm,0:km))
    !   do k=0,km
    !   do j=0,jm
    !   do i=0,im
    !     !
    !     if(lshock(i,j,k)) then
    !       rshock(i,j,k)=1.d0
    !     else
    !       rshock(i,j,k)=0.d0
    !     endif
    !     !
    !     if(crinod(i,j,k)) then
    !       rcrinod(i,j,k)=1.d0
    !     else
    !       rcrinod(i,j,k)=0.d0
    !     endif
    !     !
    !   enddo
    !   enddo
    !   enddo
    !   !
    !   call h5write(varname='lshk', var=rshock(0:im,0:jm,0:km))
    !   call h5write(varname='crit', var=rcrinod(0:im,0:jm,0:km))
    ! endif
    !
    if(trim(turbmode)=='k-omega') then
      call h5write(varname='k',     var=tke(0:im,0:jm,0:km))
      call h5write(varname='omega', var=omg(0:im,0:jm,0:km))
      call h5write(varname='miut',  var=miut(0:im,0:jm,0:km))
    elseif(trim(turbmode)=='udf1') then
      call h5write(varname='miut',  var=miut(0:im,0:jm,0:km))
    endif
    !
    call h5io_end
    !
    if(lio) then
      call h5srite(varname='nstep',var=nstep,filename=trim(outfilename))
      call h5srite(varname='time',var=time,filename=trim(outfilename))
    endif
    !
    if(ndims==1) then
      !
      open(18,file='outdat/profile'//trim(stepname)//mpirankname//'.dat')
      write(18,"(5(1X,A15))")'x','ro','u','p','t'
      write(18,"(5(1X,E15.7E3))")(x(i,0,0,1),rho(i,0,0),vel(i,0,0,1),  &
                                  prs(i,0,0),tmp(i,0,0),i=0,im)
      close(18)
      print*,' << outdat/profile',trim(stepname),mpirankname,'.dat'
      !
    endif
    !
    nxtwsequ=nstep+feqwsequ
    !
  end subroutine writeflfed
  !+-------------------------------------------------------------------+
  !| The end of the subroutine writeflfed.                             |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to write mean flow statistics.            |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 28-12-2021  | Created by J. Fang @ Warrington                     |
  !+-------------------------------------------------------------------+
  subroutine writemeanflow
    !
    use commvar,  only : nstep
    use statistic,only : nsamples,nstep_sbeg,time_sbeg,                &
                         rom,u1m,u2m,u3m,pm,tm,                        &
                         u11,u22,u33,u12,u13,u23,pp,tt,tu1,tu2,tu3,    &
                         u111,u222,u333,u112,u113,                     &
                         u122,u133,u223,u233,u123,                     &
                         u1rem,u2rem,u3rem,pu1,pu2,pu3,                &
                         sgmam11,sgmam22,sgmam33,                      &
                         sgmam12,sgmam13,sgmam23,                      &
                         disspa,predil,visdif1,visdif2,visdif3
    use hdf5io
    !
    if(lio) call bakupfile(file='outdat/meanflow.h5')
    !
    call h5io_init('outdat/meanflow.h5',mode='write')
    call h5write(varname='rom',var=rom(0:im,0:jm,0:km))
    call h5write(varname='u1m',var=u1m(0:im,0:jm,0:km))
    call h5write(varname='u2m',var=u2m(0:im,0:jm,0:km))
    call h5write(varname='u3m',var=u3m(0:im,0:jm,0:km))
    call h5write(varname='pm ',var=pm (0:im,0:jm,0:km))
    call h5write(varname='tm ',var=tm (0:im,0:jm,0:km))
    call h5io_end
    !
    if(lio) call bakupfile(file='outdat/2ndsta.h5')
    call h5io_init('outdat/2ndsta.h5',mode='write')
    call h5write(varname='u11',var=u11(0:im,0:jm,0:km))
    call h5write(varname='u22',var=u22(0:im,0:jm,0:km))
    call h5write(varname='u33',var=u33(0:im,0:jm,0:km))
    call h5write(varname='u12',var=u12(0:im,0:jm,0:km))
    call h5write(varname='u13',var=u13(0:im,0:jm,0:km))
    call h5write(varname='u23',var=u23(0:im,0:jm,0:km))
    !
    call h5write(varname='pp', var=pp(0:im,0:jm,0:km))
    call h5write(varname='tt', var=tt(0:im,0:jm,0:km))
    call h5write(varname='tu1',var=tu1(0:im,0:jm,0:km))
    call h5write(varname='tu2',var=tu2(0:im,0:jm,0:km))
    call h5write(varname='tu3',var=tu2(0:im,0:jm,0:km))
    call h5io_end
    !
    if(lio) call bakupfile(file='outdat/3rdsta.h5')
    call h5io_init('outdat/3rdsta.h5',mode='write')
    call h5write(varname='u111',var=u111(0:im,0:jm,0:km))
    call h5write(varname='u222',var=u222(0:im,0:jm,0:km))
    call h5write(varname='u333',var=u333(0:im,0:jm,0:km))
    call h5write(varname='u112',var=u112(0:im,0:jm,0:km))
    call h5write(varname='u113',var=u113(0:im,0:jm,0:km))
    call h5write(varname='u122',var=u122(0:im,0:jm,0:km))
    call h5write(varname='u133',var=u133(0:im,0:jm,0:km))
    call h5write(varname='u223',var=u223(0:im,0:jm,0:km))
    call h5write(varname='u233',var=u233(0:im,0:jm,0:km))
    call h5write(varname='u123',var=u123(0:im,0:jm,0:km))
    call h5io_end
    !
    if(lio) call bakupfile(file='outdat/budget.h5')
    call h5io_init('outdat/budget.h5',mode='write')
    call h5write(varname='disspa', var=disspa(0:im,0:jm,0:km))
    call h5write(varname='predil', var=predil(0:im,0:jm,0:km))
    call h5write(varname='pu1',    var=pu1(0:im,0:jm,0:km))
    call h5write(varname='pu2',    var=pu2(0:im,0:jm,0:km))
    call h5write(varname='pu3',    var=pu3(0:im,0:jm,0:km))
    call h5write(varname='u1rem',  var=u1rem(0:im,0:jm,0:km))
    call h5write(varname='u2rem',  var=u2rem(0:im,0:jm,0:km))
    call h5write(varname='u3rem',  var=u3rem(0:im,0:jm,0:km))
    call h5write(varname='visdif1',var=visdif1(0:im,0:jm,0:km))
    call h5write(varname='visdif2',var=visdif2(0:im,0:jm,0:km))
    call h5write(varname='visdif3',var=visdif3(0:im,0:jm,0:km))
    call h5write(varname='sgmam11',var=sgmam11(0:im,0:jm,0:km))
    call h5write(varname='sgmam22',var=sgmam22(0:im,0:jm,0:km))
    call h5write(varname='sgmam33',var=sgmam33(0:im,0:jm,0:km))
    call h5write(varname='sgmam12',var=sgmam12(0:im,0:jm,0:km))
    call h5write(varname='sgmam13',var=sgmam13(0:im,0:jm,0:km))
    call h5write(varname='sgmam23',var=sgmam23(0:im,0:jm,0:km))
    call h5io_end
    !
    if(lio) then
      call h5srite(varname='nstep',   var=nstep,   filename='outdat/meanflow.h5')
      call h5srite(varname='nsamples',var=nsamples,filename='outdat/meanflow.h5')
      call h5srite(varname='nstep_sbeg',var=nstep_sbeg,filename='outdat/meanflow.h5')
      call h5srite(varname='time_sbeg', var=time_sbeg, filename='outdat/meanflow.h5')
      call h5srite(varname='nstep',   var=nstep,   filename='outdat/2ndsta.h5')
      call h5srite(varname='nsamples',var=nsamples,filename='outdat/2ndsta.h5')
      call h5srite(varname='nstep_sbeg',var=nstep_sbeg,filename='outdat/2ndsta.h5')
      call h5srite(varname='time_sbeg', var=time_sbeg, filename='outdat/2ndsta.h5')
      call h5srite(varname='nstep',   var=nstep,   filename='outdat/3rdsta.h5')
      call h5srite(varname='nsamples',var=nsamples,filename='outdat/3rdsta.h5')
      call h5srite(varname='nstep_sbeg',var=nstep_sbeg,filename='outdat/3rdsta.h5')
      call h5srite(varname='time_sbeg', var=time_sbeg, filename='outdat/3rdsta.h5')
      call h5srite(varname='nstep',   var=nstep,   filename='outdat/budget.h5')
      call h5srite(varname='nsamples',var=nsamples,filename='outdat/budget.h5')
      call h5srite(varname='nstep_sbeg',var=nstep_sbeg,filename='outdat/budget.h5')
      call h5srite(varname='time_sbeg', var=time_sbeg, filename='outdat/budget.h5')
    endif
    !
  end subroutine writemeanflow
  !+-------------------------------------------------------------------+
  !| The end of the subroutine writemeanflow.                          |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to dump all computation data              |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 12-02-2021  | Created by J. Fang @ Warrington                     |
  !+-------------------------------------------------------------------+
  subroutine writechkpt(stepsequ,subtime)
    !
    use commvar, only: time,nstep,filenumb,fnumslic,num_species,im,jm, &
                       km,lwsequ,lavg,force,numq,imbroot,limmbou,      &
                       turbmode,feqchkpt
    use commarray,only : x,rho,vel,prs,tmp,spc,q,ssf,lshock,crinod
    use models,   only : tke,omg,miut
    use statistic,only : nsamples,liosta,massflux,massflux_target
    use bc,       only : ninflowslice
    !
    use hdf5io
    !
    !
    ! arguments
    integer,intent(in) :: stepsequ
    real(8),intent(inout),optional :: subtime
    !
    ! local data
    real(8),allocatable :: rshock(:,:,:),rcrinod(:,:,:)
    character(len=2) :: spname,qname
    integer :: jsp,i,fh,ios,j,k
    logical :: lfopen
    real(8) :: time_beg
    !
    if(present(subtime)) time_beg=ptime() 
    !
    ! if(lwsequ) then
    !   write(stepname,'(i4.4)')filenumb
    !   filenumb=filenumb+1
    !   !
    !   outfilename='outdat/flowfield'//stepname//'.h5'
    ! else
    !   stepname=''
    !   outfilename='outdat/flowfield.h5'
    ! endif
    !
    ! outfilename='outdat/flowfield'//mpirankname
    ! open(21,file=trim(outfilename),form='unformatted')
    ! write(21)nstep,time
    ! write(21)rho,vel,prs,tmp
    ! write(21)spc
    ! close(21)
    ! if(lio) print*,' << ',trim(outfilename)
    ! !
    ! outfilename='outdat/qdata'//mpirankname
    ! open(21,file=trim(outfilename),form='unformatted')
    ! write(21)nstep
    ! write(21)q
    ! close(21)
    ! if(lio) print*,' << ',trim(outfilename)
    ! !
    ! outfilename='outdat/auxiliary'//mpirankname
    ! open(21,file=trim(outfilename),form='unformatted')
    ! write(21)nstep,filenumb
    ! write(21)massflux,massflux_target
    ! write(21)force
    ! close(21)
    ! if(lio) print*,' << ',trim(outfilename)
    !
    if(lio) call bakupfile('outdat/auxiliary.h5')
    if(lio) call bakupfile('outdat/flowfield.h5')
    !
    call writeflfed()
    !
    ! call h5io_init('outdat/qdata.h5',mode='write')
    ! do jsp=1,numq
    !    write(qname,'(i2.2)')jsp
    !   call h5write(varname='q'//qname,var=q(0:im,0:jm,0:km,jsp))
    ! enddo
    ! call h5io_end
    !
    if(lio) then
      !
      call h5srite(varname='nstep',var=nstep,                          &
                          filename='outdat/auxiliary.h5',newfile=.true.)
      call h5srite(varname='filenumb',var=filenumb,                    &
                                         filename='outdat/auxiliary.h5')
      call h5srite(varname='fnumslic',var=fnumslic,                    &
                                         filename='outdat/auxiliary.h5')
      call h5srite(varname='ninflowslice',var=ninflowslice,            &
                                         filename='outdat/auxiliary.h5')
      call h5srite(varname='massflux',var=massflux,                    &
                                         filename='outdat/auxiliary.h5')
      call h5srite(varname='massflux_target',var=massflux_target,      &
                                         filename='outdat/auxiliary.h5')
      call h5srite(varname='force',var=force,                          &
                                         filename='outdat/auxiliary.h5')
      call h5srite(varname='nsamples',var=nsamples,                    &
                                         filename='outdat/auxiliary.h5')
    endif
    !
    if(limmbou .and. mpirank==imbroot) then
      ! call imboundarydata 
      !
      ! call writecylimmbou
      !
    endif
    !
    ! if(irk==0 .and. jrk==jrkm) then
    !     print*,'------------------------------------'
    !     print*,x(0,jm,0,1),prs(0,jm,0)
    !     print*,'------------------------------------'
    ! endif
    !
    if(liosta) then
      !
      call writemeanflow
      !
    endif
    !
    ! flush all openned files
    fh=20
    lfopen=.true.
    !
    do while(lfopen)
      inquire(unit=fh, opened=lfopen,iostat=ios)
      !
      if(lfopen .and. ios==0) then
        flush(fh)
        fh=fh+1
      elseif(ios .ne. 0) then
        print*,' !! error with opening file unit:',fh
      endif
    enddo
    !
    nxtchkpt=nstep+feqchkpt
    !
    if(present(subtime)) subtime=subtime+ptime()-time_beg 
    !
    return
    !
  end subroutine writechkpt
  !+-------------------------------------------------------------------+
  !| The end of the subroutine writechkpt.                             |
  !+-------------------------------------------------------------------+
  !!
  !+-------------------------------------------------------------------+
  !| This subroutine is used to write slice cut.                       |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 10-10-2021  | Created by J. Fang @ Warrington                     |
  !+-------------------------------------------------------------------+
  subroutine writeslice
    !
    use parallel, only: irk_islice,jrk_jslice,krk_kslice,mpi_islice,   &
                        mpi_jslice,mpi_kslice
    use commvar,  only: fnumslic,nstep,time,islice,jslice,kslice,im,jm,km
    use commarray,only: rho,vel,prs,tmp,dvel
    use hdf5io
    !
    ! local data
    character(len=5) :: stepname
    character(len=14) :: filename
    integer :: i,j,k
    !
    if(irk==irk_islice) then
      !
      i=islice-ig0
      !
      if(i<0 .or. i>im) then
        print*,' !! islice location error !! i=',i
        stop
      endif
      !
      call system('mkdir islice/')
      !
      write(stepname,'(i5.5)')fnumslic
      !
      filename='islice'//stepname//'.h5'
      !
      call h5io_init(filename='islice/'//filename,mode='write',        &
                                                        comm=mpi_islice)
      call h5write(varname='ro',var=rho(i,0:jm,0:km),  dir='i')
      call h5write(varname='u1',var=vel(i,0:jm,0:km,1),dir='i')
      call h5write(varname='u2',var=vel(i,0:jm,0:km,2),dir='i')
      call h5write(varname='u3',var=vel(i,0:jm,0:km,3),dir='i')
      call h5write(varname='p', var=prs(i,0:jm,0:km)  ,dir='i')
      call h5write(varname='t', var=tmp(i,0:jm,0:km)  ,dir='i')
      !
      call h5write(varname='dudx',var=dvel(i,0:jm,0:km,1,1),dir='i')
      call h5write(varname='dudy',var=dvel(i,0:jm,0:km,1,2),dir='i')
      call h5write(varname='dudz',var=dvel(i,0:jm,0:km,1,3),dir='i')
      call h5write(varname='dvdx',var=dvel(i,0:jm,0:km,2,1),dir='i')
      call h5write(varname='dvdy',var=dvel(i,0:jm,0:km,2,2),dir='i')
      call h5write(varname='dvdz',var=dvel(i,0:jm,0:km,2,3),dir='i')
      call h5write(varname='dwdx',var=dvel(i,0:jm,0:km,3,1),dir='i')
      call h5write(varname='dwdy',var=dvel(i,0:jm,0:km,3,2),dir='i')
      call h5write(varname='dwdz',var=dvel(i,0:jm,0:km,3,3),dir='i')
      !
      call h5io_end
      !
      if(jrk==0 .and. krk==0) then
        call h5srite(filename='islice/'//filename,varname='nstep',var=nstep)
        call h5srite(filename='islice/'//filename,varname='time', var=time)
      endif
      !
    endif
    !
    if(jrk==jrk_jslice) then
      !
      j=jslice-jg0
      !
      if(j<0 .or. j>jm) then
        print*,' !! jslice location error !! j=',j
        stop
      endif
      !
      call system('mkdir jslice/')
      !
      write(stepname,'(i5.5)')fnumslic
      !
      filename='jslice'//stepname//'.h5'
      !
      call h5io_init(filename='jslice/'//filename,mode='write',        &
                                                        comm=mpi_jslice)
      call h5write(varname='ro',var=rho(0:im,j,0:km),  dir='j')
      call h5write(varname='u1',var=vel(0:im,j,0:km,1),dir='j')
      call h5write(varname='u2',var=vel(0:im,j,0:km,2),dir='j')
      call h5write(varname='u3',var=vel(0:im,j,0:km,3),dir='j')
      call h5write(varname='p', var=prs(0:im,j,0:km)  ,dir='j')
      call h5write(varname='t', var=tmp(0:im,j,0:km)  ,dir='j')
      !
      call h5write(varname='dudx',var=dvel(0:im,j,0:km,1,1),dir='j')
      call h5write(varname='dudy',var=dvel(0:im,j,0:km,1,2),dir='j')
      call h5write(varname='dudz',var=dvel(0:im,j,0:km,1,3),dir='j')
      call h5write(varname='dvdx',var=dvel(0:im,j,0:km,2,1),dir='j')
      call h5write(varname='dvdy',var=dvel(0:im,j,0:km,2,2),dir='j')
      call h5write(varname='dvdz',var=dvel(0:im,j,0:km,2,3),dir='j')
      call h5write(varname='dwdx',var=dvel(0:im,j,0:km,3,1),dir='j')
      call h5write(varname='dwdy',var=dvel(0:im,j,0:km,3,2),dir='j')
      call h5write(varname='dwdz',var=dvel(0:im,j,0:km,3,3),dir='j')
      !
      call h5io_end
      !
      if(irk==0 .and. krk==0) then
        call h5srite(filename='jslice/'//filename,varname='nstep',var=nstep)
        call h5srite(filename='jslice/'//filename,varname='time', var=time)
      endif
      !
    endif
    !
    if(krk==krk_kslice) then
      !
      k=kslice-kg0
      !
      if(k<0 .or. k>km) then
        print*,' !! kslice location error !! k=',k
        stop
      endif
      !
      call system('mkdir kslice/')
      !
      write(stepname,'(i5.5)')fnumslic
      !
      filename='kslice'//stepname//'.h5'
      !
      call h5io_init(filename='kslice/'//filename,mode='write',        &
                                                        comm=mpi_kslice)
      call h5write(varname='ro',var=rho(0:im,0:jm,k),  dir='k')
      call h5write(varname='u1',var=vel(0:im,0:jm,k,1),dir='k')
      call h5write(varname='u2',var=vel(0:im,0:jm,k,2),dir='k')
      call h5write(varname='u3',var=vel(0:im,0:jm,k,3),dir='k')
      call h5write(varname='p', var=prs(0:im,0:jm,k)  ,dir='k')
      call h5write(varname='t', var=tmp(0:im,0:jm,k)  ,dir='k')
      !
      call h5write(varname='dudx',var=dvel(0:im,0:jm,k,1,1),dir='k')
      call h5write(varname='dudy',var=dvel(0:im,0:jm,k,1,2),dir='k')
      call h5write(varname='dudz',var=dvel(0:im,0:jm,k,1,3),dir='k')
      call h5write(varname='dvdx',var=dvel(0:im,0:jm,k,2,1),dir='k')
      call h5write(varname='dvdy',var=dvel(0:im,0:jm,k,2,2),dir='k')
      call h5write(varname='dvdz',var=dvel(0:im,0:jm,k,2,3),dir='k')
      call h5write(varname='dwdx',var=dvel(0:im,0:jm,k,3,1),dir='k')
      call h5write(varname='dwdy',var=dvel(0:im,0:jm,k,3,2),dir='k')
      call h5write(varname='dwdz',var=dvel(0:im,0:jm,k,3,3),dir='k')
      !
      call h5io_end
      !
      if(irk==0 .and. jrk==0) then
        call h5srite(filename='kslice/'//filename,varname='nstep',var=nstep)
        call h5srite(filename='kslice/'//filename,varname='time', var=time)
      endif
      !
    endif
    !
    fnumslic=fnumslic+1
    !
  end subroutine writeslice
  !+-------------------------------------------------------------------+
  !| The end of the subroutine output.                                 |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to write data at the immersed boundary.   |
  !| for cylinder configuration.                                       |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 20-07-2021  | Created by J. Fang @ Warrington                     |
  !+-------------------------------------------------------------------+
  subroutine writecylimmbou
    !
    use commtype,  only : sboun
    use commvar,   only : immbond,num_species,pinf
    use fludyna,   only : fvar2q,q2fvar
    !
    ! local data
    integer :: fh,jb,kb,nb
    real(8) :: xc,yc,radi,varmin
    real(8) :: vel_bou(3),prs_bou,rho_bou,tmp_bou,spc_bou(num_species)
    real(8),allocatable,save :: alfa(:)
    integer,allocatable,save :: ialf(:)
    type(sboun),pointer :: pb
    !
    logical,save :: firstcall=.true.
    !
    nb=size(immbond)
    !
    if(firstcall) then
      !
      xc=5.d0
      yc=5.d0
      ! radi=0.5d0
      !
      allocate(alfa(nb),ialf(nb))
      !
      do jb=1,nb
        !
        pb=>immbond(jb)
        !
        radi=sqrt((pb%x(1)-xc)**2+(pb%x(2)-yc)**2)
        !
        if(pb%x(2)-yc>=0.d0) then
          alfa(jb)=pi-acos((pb%x(1)-xc)/radi)
        else
          alfa(jb)=pi+acos((pb%x(1)-xc)/radi)
        endif
        !
      enddo
      !
      ! sort the array according to the order of alfa
      ialf=0
      do kb=1,nb
        !
        varmin=1.d10
        do jb=1,nb
          !
          if(any(ialf(1:kb)==jb)) cycle
          !
          if(alfa(jb)<varmin) then
            varmin=alfa(jb)
            ialf(kb)=jb
          endif
          !
        enddo
        !
      enddo
      !
      firstcall=.false.
      !
    endif
    !
    ! do jb=1,nb
    !   write(*,*)jb,ialf(jb),alfa(ialf(jb))
    !   print*,'------------------------------------------'
    ! enddo
    !
    fh=get_unit()
    open(fh,file='outdat/immbou.dat')
    write(fh,'(2(1X,A20))')'alfa','p'
    !
    do jb=1,nb
      !
      kb=ialf(jb)
      !
      pb=>immbond(kb)
      !
      call q2fvar(      q=  pb%qimag,                 &
                    density=rho_bou,                  &
                   velocity=vel_bou(:),               &
                   pressure=prs_bou,                  &
                temperature=tmp_bou,                  &
                    species=spc_bou                   )
      !
      write(fh,'(2(1X,E20.13E2))')alfa(kb),2.d0*(prs_bou-pinf)
      !
    enddo
    !
    close(fh)
    print*,' << outdat/immbou.dat'
    !
    ! call mpistop
    !
    ! fh=get_unit()
    !
  end subroutine writecylimmbou
  !+-------------------------------------------------------------------+
  !| The end of the subroutine writecylimmbou.                         |
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
    use commvar, only : nstep,maxstep,ctime,flowtype,conschm,          &
                        difschm,rkscheme,ia,ja,ka,preptime,nsrpt
    use parallel,only : mpirankmax
    !
    ! local data
    logical,save :: linit=.true.
    logical :: lexist
    integer :: i
    integer,save :: hand_rp,repsp
    !
    if(lio) then
      !
      if(linit) then
        !
        preptime=ptime()-preptime
        !
        inquire(file='report.txt', exist=lexist)
        !
        if(lexist) call system('mv -v report.txt report.bak')
        !
        call system('echo "----------------------------------------------------------------" > report.txt')
        call system('echo "CPU infomation" >> report.txt')
        call system('echo "----------------------------------------------------------------" >> report.txt')
        call system('lscpu | grep "Model name" >> report.txt')
        call system('lscpu | grep "CPU MHz" >> report.txt')
        call system('lscpu | grep "Socket(s):" >> report.txt')
        call system('lscpu | grep "Core(s) per socket:" >> report.txt')
        call system('lscpu | grep "Thread(s) per core:" >> report.txt')
        call system('lscpu | grep "cache" >> report.txt')
        call system('echo "----------------------------------------------------------------" >> report.txt')
        !
        hand_rp=get_unit()
        !
        open(hand_rp,file='report.txt',position="append")
        write(hand_rp,'(A)')'  statistic of computing time'
        write(hand_rp,'(A,A)')'     flowtype: ',trim(flowtype)
        write(hand_rp,'(A,A)')'  conv scheme: ',trim(conschm)
        write(hand_rp,'(A,A)')'  diff scheme: ',trim(difschm)
        write(hand_rp,'(A,A)')'    rk scheme: ',trim(rkscheme)
        write(hand_rp,'(4(A,I0))')'    grid size: ',ia,' x ',ja,' x ', &
                                           ka,' = ',(ia+1)*(ja+1)*(ka+1)
        write(hand_rp,'(1(A,I0))')'     mpi size: ',mpirankmax+1
        write(hand_rp,'(2X,A,E13.6E2)')'time cost for preparation : ',preptime
        !
        repsp=0
        !
        linit=.false.
        !
      endif
      !
      write(hand_rp,'(2X,62A)')('-',i=1,62)
      write(hand_rp,'(2X,2(A,I7))')'time report from nstep : ',repsp,  &
                                                ' to nstep : ',nstep
      write(hand_rp,'(2X,62A)')('-',i=1,62)
      write(hand_rp,'(2X,A,E13.6E2,A,F6.2,A)')'total time cost : ',    &
                            ctime(2),' - ',100.d0*ctime(2)/ctime(2),' %'
      !
      write(hand_rp,'(2X,A,E13.6E2,A,F6.2,A)')'  - rk          : ',    &
                            ctime(3),' - ',100.d0*ctime(3)/ctime(2),' %'
      write(hand_rp,'(2X,A,E13.6E2,A,F6.2,A)')'    - rhs       : ',    &
                            ctime(4),' - ',100.d0*ctime(4)/ctime(2),' %'
      write(hand_rp,'(2X,A,E13.6E2,A,F6.2,A)')'      - Convc   : ',    &
                            ctime(9),' - ',100.d0*ctime(9)/ctime(2),' %'
      write(hand_rp,'(2X,A,E13.6E2,A,F6.2,A)')'        - ssf   : ',    &
                          ctime(12),' - ',100.d0*ctime(12)/ctime(2),' %'
      write(hand_rp,'(2X,A,E13.6E2,A,F6.2,A)')'      - Diffu   : ',    &
                          ctime(10),' - ',100.d0*ctime(10)/ctime(2),' %'
      write(hand_rp,'(2X,A,E13.6E2,A,F6.2,A)')'    - filter    : ',    &
                            ctime(8),' - ',100.d0*ctime(8)/ctime(2),' %'
      write(hand_rp,'(2X,A,E13.6E2,A,F6.2,A)')'    - io        : ',    &
                            ctime(6),' - ',100.d0*ctime(6)/ctime(2),' %'
      write(hand_rp,'(2X,A,E13.6E2,A,F6.2,A)')'    - sta       : ',    &
                            ctime(5),' - ',100.d0*ctime(5)/ctime(2),' %'
      write(hand_rp,'(2X,A,E13.6E2,A,F6.2,A)')'    - bc        : ',    &
                          ctime(11),' - ',100.d0*ctime(11)/ctime(2),' %'
      write(hand_rp,'(2X,A,E13.6E2,A,F6.2,A)')'  - com         : ',    &
                           ctime(7), ' - ',100.d0*ctime(7)/ctime(2),' %'
      !
      flush(hand_rp)
      !
      if(nstep>=maxstep) then
        close(hand_rp)
        print*,' << report.txt'
      endif
      !
    endif
    !
    repsp=nstep
    ctime=0.d0
    !
  end subroutine timerept
  !+-------------------------------------------------------------------+
  !| The end of the subroutine timerept.                               |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to write immersed boundary structure           |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 18-Aug-2021: Created by J. Fang @ Appleton                        |
  !+-------------------------------------------------------------------+
  subroutine write_sboun(abound,item2w)
    !
    use commtype, only : sboun
    use commvar,  only : ndims
    use commarray,only : x,cell
    use commcal,  only : ijkin,ijkcellin
    !
    ! arguments
    type(sboun),intent(in) :: abound(:)
    character(len=*),intent(in) :: item2w
    !
    ! local data
    integer :: counter,n,n1,m,i,j,k,fh
    real(8),allocatable :: xcell(:,:)
    !
    select case(item2w)
    case('bound')
      !
      counter=0
      do n=1,size(abound)
        if(abound(n)%localin) counter=counter+1
      enddo
      !
      if(counter>0) then
        fh=get_unit()
        open(fh,file='tecbound'//mpirankname//'.dat')
        write(fh,'(A)')' VARIABLES = "x" "y" "z" '
        write(fh,'(A)')'ZONE T="ZONE 001"'
        write(fh,'(A,I0,A)')'I=',counter,', J=1, K=1, ZONETYPE=Ordered'
        write(fh,'(A)')'DATAPACKING=POINT'
        do n=1,size(abound)
          if(abound(n)%localin) then
            write(fh,'(3(1X,E20.13E2))')abound(n)%x(:)
          endif
        enddo
        close(fh)
        print*,' << tecbound,',mpirankname,'.dat'
      endif
      !
    case('ghost')
      !
      counter=0
      do n=1,size(abound)
        i=abound(n)%igh(1)-ig0
        j=abound(n)%igh(2)-jg0
        k=abound(n)%igh(3)-kg0
        if(ijkin(i,j,k)) counter=counter+1
      enddo
      !
      if(counter>0) then
        !
        fh=get_unit()
        open(fh,file='tecghost'//mpirankname//'.dat')
        write(fh,'(A)')' VARIABLES = "x" "y" "z" '
        write(fh,'(A)')'ZONE T="ZONE 001"'
        write(fh,'(A,I0,A)')'I=',counter,', J=1, K=1, ZONETYPE=Ordered'
        write(fh,'(A)')'DATAPACKING=POINT'
        do n=1,size(abound)
          i=abound(n)%igh(1)-ig0
          j=abound(n)%igh(2)-jg0
          k=abound(n)%igh(3)-kg0
          if(ijkin(i,j,k)) then
            write(fh,'(3(1X,E20.13E2))')x(i,j,k,:)
          endif
        enddo
        close(fh)
        print*,' << tecghost',mpirankname,'.dat'
        !
      endif
      !
    case('image')
      !
      counter=0
      do n=1,size(abound)
        if(abound(n)%localin) counter=counter+1
      enddo
      !
      if(mpirank==0) then
        fh=get_unit()
        open(fh,file='tecimage'//mpirankname//'.dat')
        write(fh,'(A)')' VARIABLES = "x" "y" "z" '
        write(fh,'(A)')'ZONE T="ZONE 001"'
        write(fh,'(A,I0,A)')'I=',size(abound),', J=1, K=1, ZONETYPE=Ordered'
        write(fh,'(A)')'DATAPACKING=POINT'
        do n=1,size(abound)
          write(fh,'(3(1X,E20.13E2))')abound(n)%ximag(:)
        enddo
        close(fh)
        print*,' << tecimage',mpirankname,'.dat'
      endif
      !
    case('icell')
      !
      allocate(xcell(1:2**ndims,1:3))
      !
      counter=0
      do n=1,size(abound)
        !
        i=abound(n)%icell(1)-ig0
        j=abound(n)%icell(2)-jg0
        k=abound(n)%icell(3)-kg0
        !
        if(ijkcellin(i,j,k)) then
          counter=counter+1
        endif
        !
      enddo
      !
      if(counter>0) then
        !
        fh=get_unit()
        open(fh,file='tecicell'//mpirankname//'.dat')
        !
        if(ndims==2) then
          write(fh,'(A)')'TITLE = "FE-Volume QUADRILATERAL Data"'
          write(fh,'(a)')'variables = "x", "y", "z"'
          write(fh,'((A,I0,A))')'ZONE T="',mpirank,'"'
          write(fh,'(2(A,I0),A)')'DATAPACKING=POINT, NODES= ',counter*4, &
                                 ', ELEMENTS= ',counter,', ZONETYPE=FEQUADRILATERAL'
        elseif(ndims==3) then
          write(fh,'(A)')'TITLE = "FE-Volume BRICK Data"'
          write(fh,'(a)')'variables = "x", "y", "z"'
          write(fh,'((A,I0,A))')'ZONE T="',mpirank,'"'
          write(fh,'(2(A,I0),A)')'DATAPACKING=POINT, NODES= ',counter*8, &
                                 ', ELEMENTS= ',counter,', ZONETYPE=FEBRICK'
        endif
        !
        do n=1,size(abound)
          !
          i=abound(n)%icell(1)-ig0
          j=abound(n)%icell(2)-jg0
          k=abound(n)%icell(3)-kg0
          !
          if(ijkcellin(i,j,k)) then
            !
            write(fh,'(3(1X,E20.13E2))')cell(i,j,k)%x(1,:)
            write(fh,'(3(1X,E20.13E2))')cell(i,j,k)%x(2,:)
            write(fh,'(3(1X,E20.13E2))')cell(i,j,k)%x(3,:)
            write(fh,'(3(1X,E20.13E2))')cell(i,j,k)%x(4,:)
            !
            if(ndims==3 ) then
              write(fh,'(3(1X,E20.13E2))')cell(i,j,k)%x(5,:)
              write(fh,'(3(1X,E20.13E2))')cell(i,j,k)%x(6,:)
              write(fh,'(3(1X,E20.13E2))')cell(i,j,k)%x(7,:)
              write(fh,'(3(1X,E20.13E2))')cell(i,j,k)%x(8,:)
            endif
            !
          endif
          !
        enddo
        write(fh,'(8(1X,I0))')(m,m=1,counter*2**ndims)
        close(fh)
        print*,' << tecicell',mpirankname,'.dat'
        !
      endif
      !
      deallocate(xcell)
      !
    case default
      stop ' ERROR 1, item2w not defined @ write_sboun'
    end select
    !
    !
  end subroutine write_sboun
  !+-------------------------------------------------------------------+
  !| The end of the subroutine write_sboun.                            |
  !+-------------------------------------------------------------------+
  !!
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
end module readwrite
!+---------------------------------------------------------------------+
!| The end of the module readwrite.                                    |
!+---------------------------------------------------------------------+
