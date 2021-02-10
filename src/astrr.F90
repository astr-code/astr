!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This Program is a DNS solver for reacting flow 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! re-mastered at 2021-02-06, by Fang Jian 
! fangjian19@gmail.com
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program astrr
  !
  use parallel
  use readwrite
  use commarray
  use solver
  use initialisation
  use mainloop
  !
  implicit none
  !
  call mpiinitial
  !
  call statement
  !
  call readinput
  !
  call mpisizedis
  !
  call parapp
  !
  call parallelini
  !
  call refcal
  !
  call infodisp
  !
  call allocommarray
  !
  call readgrid
  !
  call geomcal
  !
  call solvrinit
  !
  call flowinit
  !
  call steploop
  !
  call mpistop
  !
end program astrr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! End of the program main
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
