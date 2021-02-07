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
  call infodisp
  !
  call allocommarray
  !
  call ReadGrid
  !
  if(mpirank==0) print*,' ** The job is done!'
  call mpistop
  !
end program astrr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! End of the program main
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
