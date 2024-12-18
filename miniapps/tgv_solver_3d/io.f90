module io
   
   implicit none

   contains

   subroutine readinput

     use utility
     use comvardef

     character(len=128) :: input_file_name
     integer :: file_handle

     namelist /basicparam/ im,jm,km,hm,numq

     call readkeyboad(input_file_name)

     write(*,'(A)',advance='no')'  >> '//trim(input_file_name)//'...'
     
     ! Read parameters
     open(newunit=file_handle, file=trim(input_file_name))

     ! These are the 'essential' parameters
     read(file_handle, nml=basicparam); rewind(file_handle)

     close(file_handle)
     print*,' done. '

   end subroutine readinput

   subroutine casereport
     
     use utility
     use comvardef
     !$ use omp_lib

     !$ print*,' !! OpenMP Activated !!'

     !$ write(*,'(A,I0)')'  ** number of omp cores:',omp_get_num_procs()
     write(*,'(3(A,I0))')  '  ** number of equations:',numq
     write(*,'(3(A,I0))')  '  ** dimesion           :',im,'x',jm,'x',km
     write(*,'(3(A,I0))')  '  ** halo level         :',hm

   end subroutine casereport

end module io