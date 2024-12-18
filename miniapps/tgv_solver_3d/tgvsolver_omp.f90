program astr_mini_omp
  
  use utility
  use comvardef
  use io
  use field_init
  use solver

  implicit none

  real :: tstart,tfinish

  call statement

  call cpu_time(tstart)

  call readinput

  call casereport

  call alloarray

  call flowinit('tgv')

  call mainloop

  call cpu_time(tfinish)
  
  print*,' ** time cost from begining to end: ',tfinish-tstart

  ctime(1)=tfinish-tstart
  
  call timereport

  stop

end program astr_mini_omp