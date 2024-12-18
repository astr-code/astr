program astr_mini_omp
  
  use utility
  use comvardef
  use io
  use field_init
  use solver
  !$ use omp_lib

  implicit none

  real :: tstart,tfinish

  call statement

  nthread=16
  !$ write(*,'(A,I0,A)')'  ** The program is working in OpenMP with ',nthread,' threads'
  !$ call OMP_SET_DYNAMIC(.false.)
  !$ if (OMP_GET_DYNAMIC()) print *,' **!! Warning: dynamic adjustment of threads has been set!'
  !$ call OMP_SET_NUM_THREADS(nthread)
  !$ ncore=OMP_get_num_procs()
  !$ if(nthread>ncore) print *,' **!! Warning: the thread No. is larger than processers No.!'
  !
  !$ mthread=omp_get_num_threads()
  !$ write(*,'(A,I0,A)')'  ** The real number of threads ',mthread
  !
  tstart=time_in_second()

  call readinput

  call casereport

  call alloarray

  call flowinit('tgv')

  call mainloop

  tfinish=time_in_second()
  
  print*,' ** time cost from begining to end: ',tfinish-tstart

  ctime(1)=tfinish-tstart
  
  call timereport

  stop

end program astr_mini_omp