module io
   
   implicit none

   contains

   subroutine readinput

     use utility
     use comvardef, only: im,jm,km,hm,numq,ref_t,gamma,mach,reynolds, &
                          prandtl,deltat,maxstep,flowtype,ldiffusion, &
                          lfilter,fdmform,len_sten_max

     character(len=128) :: input_file_name
     integer :: file_handle

     namelist /basicparam/ im,jm,km,hm,numq,ref_t,gamma,mach,reynolds, &
                           prandtl,deltat,maxstep,flowtype,ldiffusion, &
                           lfilter,fdmform,len_sten_max

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

     write(*,'(3(A,I0))')     '  ** number of equations:',numq
     write(*,'(3(A,I0))')     '  ** dimesion           :',im,'x',jm,'x',km
     write(*,'(3(A,I0))')     '  ** halo level         :',hm
     write(*,'(3(A,I0))')     '  ** maxstep            :',maxstep
     if(ldiffusion) then
        write(*,'((A))')     '  ** diffusion term  activated'
        write(*,'(3(A,F10.5))')  '  ** reynolds           :',reynolds
     endif
     if(lfilter) then
        write(*,'((A))')     '  ** filter activated'
     endif
     write(*,'(3(A,I0))')     '  ** length of stencil  :',len_sten_max
     write(*,'(2(A))')        '  ** finite-difference  :',fdmform
     write(*,'(3(A,F10.5))')  '  ** mach               :',mach
     write(*,'(3(A,F10.5))')  '  ** prandtl            :',prandtl
     write(*,'(3(A,F10.5))')  '  ** gamma              :',gamma
     write(*,'(3(A,F10.5))')  '  ** deltat             :',deltat
     write(*,'(3(A,A))')      '  ** flowtype           :',trim(flowtype)

   end subroutine casereport

   subroutine write_field

    use comvardef, only: x,vel,dvel,prs,file_number,im,jm,km
    use tecio, only: tecbin

    character(len=4) :: file_name

    integer :: i
    real(8) :: omega

    file_number=file_number+1

    write(file_name,'(i4.4)')file_number

    call tecbin(filename='flow'//file_name//'.plt',var1=x(0:im,0:jm,0:km,1),  var1name='x', &
                                       var2=x(0:im,0:jm,0:km,2),  var2name='y', &
                                       var3=x(0:im,0:jm,0:km,3),  var3name='z', &
                                       var4=vel(0:im,0:jm,0:km,1),var4name='u', &
                                       var5=vel(0:im,0:jm,0:km,2),var5name='v', &
                                       var6=prs(0:im,0:jm,0:km),  var6name='p')

    open(18,file='xprofile'//file_name//'.dat')
    do i=0,im
      omega=dvel(i,jm/2,0,2,1)-dvel(i,jm/2,0,1,2)
      write(18,*)x(i,jm/2,0,1),vel(i,jm/2,0,1),vel(i,jm/2,0,2),prs(i,jm/2,0),omega
    enddo
    close(18)
    print*,' << '//'xprofile'//file_name//'.dat'

   end subroutine write_field


end module io