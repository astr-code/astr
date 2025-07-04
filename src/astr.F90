!=======================================================================
! DNS Solver for Compressible Flow - ASTR
!=======================================================================
! Re-mastered on: 2021-02-06
! Author        : Fang Jian
! Contact       : fangjian19@gmail.com
!=======================================================================
program astr

  !---------------------------------------------------------------------
  ! Module usage
  !---------------------------------------------------------------------
  use parallel,      only: mpiinitial,mpistop,mpisizedis,lio,          &
                           parallelini,parapp
  use readwrite,     only: statement,readinput,fileini,infodisp
  use commarray,     only: allocommarray
  use solver,        only: refcal
  use initialisation,only: spongelayerini,flowinit
  use mainloop,      only: steploop
  use gridgeneration,only: gridgen
  use cmdefne,       only: getcmd,listcmd
  use pp,            only: ppentrance
  use geom,          only: geomcal
  use ibmethod,      only: ibprocess
  use test,          only: codetest
  use comsolver,     only: solvrinit

  implicit none

  !---------------------------------------------------------------------
  ! Local variables
  !---------------------------------------------------------------------
  character(len=16) :: cmd

  !---------------------------------------------------------------------
  ! MPI Initialization and Command Processing
  !---------------------------------------------------------------------
  call mpiinitial

  call statement

  call listcmd

  call getcmd(cmd)

  !---------------------------------------------------------------------
  ! Select operation based on command
  !---------------------------------------------------------------------
  if (trim(cmd) == 'pp') then

    ! Pre/Post-processing
    call ppentrance

  elseif (trim(cmd) == 'test') then

    ! code tests
    call codetest

  elseif (trim(cmd) == 'run') then

    ! Main simulation run
    call readinput

    call mpisizedis

    call parapp

    call parallelini

    call refcal

    call fileini

    call infodisp

    call allocommarray

    call ibprocess

    call gridgen

    call solvrinit

    call geomcal

    call spongelayerini

    call flowinit

    call steploop

    call mpistop

  else

    if (lio) print *, ' ** All jobs (nothing) done. **'
    stop

  end if

end program astr
!=======================================================================
! End of program ASTR
!=======================================================================
