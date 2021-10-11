!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This Program is a DNS solver for compressible flow 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! re-mastered at 2021-02-06, by Fang Jian 
! fangjian19@gmail.com
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program astr
  !
  use parallel
  use readwrite
  use commarray
  use solver
  use initialisation
  use mainloop
  use gridgeneration
  use cmdefne
  use pp
  use geom
  !
  implicit none
  !
  character(len=16) :: cmd
  !
  call mpiinitial
  !
  call statement
  !
  call getcmd(cmd)
  !
  call listcmd
  !
  if(trim(cmd)=='pp') then
    !
    ! generate an example channel flow case
    !
    call ppentrance
    !
  elseif(trim(cmd)=='run') then
    !
    ! computational loop
    !
    call readinput
    !
    call readmonc
    !
    call mpisizedis
    !
    call parapp
    !
    call parallelini
    !
    call refcal
    !
    call fileini
    !
    call infodisp
    !
    call allocommarray
    !
    call gridgen
    !
    call geomcal
    !
    call spongelayerini
    !
    call solvrinit
    !
    call flowinit
    !
    call steploop
    ! call gradcal
    !
    call mpistop
    !
  else
    if(mpirank==0) print*,' ** all jobs done. **'
    stop
  endif
  !
end program astr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! End of the program main
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
