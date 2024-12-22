module io
   
   implicit none

   contains

   subroutine readinput

     use utility
     use comvardef

     character(len=128) :: input_file_name
     integer :: file_handle

     namelist /basicparam/ im,jm,km,hm,numq,ref_t,gamma,mach,reynolds,prandtl,deltat,maxstep

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

     !$ write(*,'(A,I0)')'  ** number of omp threads:',omp_get_num_threads()
        write(*,'(3(A,I0))')  '  ** number of equations:',numq
        write(*,'(3(A,I0))')  '  ** dimesion           :',im,'x',jm,'x',km
        write(*,'(3(A,I0))')  '  ** halo level         :',hm
        write(*,'(3(A,I0))')  '  ** maxstep            :',maxstep
     write(*,'(3(A,F10.5))')  '  ** reynolds           :',reynolds
     write(*,'(3(A,F10.5))')  '  ** mach               :',mach
     write(*,'(3(A,F10.5))')  '  ** prandtl            :',prandtl
     write(*,'(3(A,F10.5))')  '  ** gamma              :',gamma
     write(*,'(3(A,F10.5))')  '  ** deltat             :',deltat

   end subroutine casereport


end module io