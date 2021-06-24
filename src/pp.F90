!+---------------------------------------------------------------------+
!| This module contains subroutines for pre and post-process           |
!+---------------------------------------------------------------------+
!| ==============                                                      |
!| CHANGE RECORD                                                       |
!| -------------                                                       |
!|  28-05-2021  | Created by J. Fang @ Warrington                      |
!+---------------------------------------------------------------------+
module pp
  !
  implicit none
  !
  contains
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to generate an example case.              |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 28-05-2021  | Created by J. Fang @ Warrington                     |
  !+-------------------------------------------------------------------+
  subroutine examplegen(folder)
    !
    use parallel, only: mpirank
    !
    character(len=*),intent(in) :: folder
    !
    if(mpirank==0) then
      !
      print*,' ** Generating an example case.'
      !
      call system('mkdir '//trim(folder))
      call system('mkdir '//trim(folder)//'/datin')
      !
      open(16,file=trim(folder)//'/datin/input.chl',form='formatted')
      write(16,'(A)')'########################################################################'
      write(16,'(A)')'#                     input file of ASTR code                          #'
      write(16,'(A)')'########################################################################'
      write(16,*)
      write(16,'(A)')'# flowtype                                              : The type of flow problem'
      write(16,'(A)')'channel'
      write(16,*)
      write(16,'(A)')'# im,jm,km                                              : The size of grid.'
      write(16,'(A)')'128,128,128'
      write(16,*)
      write(16,'(A)')'# lihomo,ljhomo,lkhomo                                  : The homogeneous directions'
      write(16,'(A)')'t,f,t'
      write(16,*)
      write(16,'(A)')'# nondimen,diffterm,lfilter,lreadgrid,lfftz             : Parameters'
      write(16,'(A)')'t,t,t,t,f'
      write(16,*)
      write(16,'(A)')'# lrestar                                               : start mode'
      write(16,'(A)')'f'
      write(16,*)
      write(16,'(A)')'# alfa_filter, kcutoff                                  : Filter parameters'
      write(16,'(A)')'0.49d0, 48'
      write(16,*)
      write(16,'(A)')'# ref_t,reynolds,mach                                   : Reference variables'
      write(16,'(A)')'273.15d0,  3000.d0,  0.5d0 '
      write(16,*)
      write(16,'(A)')'# conschm,difschm,rkscheme                              : Numerical scheme'
      write(16,'(A)')'642c, 642c, rk4 '
      write(16,*)
      write(16,'(A)')'# recon_schem, lchardecomp                              : Parameters for upwind-biased scheme'
      write(16,'(A)')'0, f '
      write(16,*)
      write(16,'(A)')'# num_species : number of species'
      write(16,'(A)')'1'
      write(16,*)
      write(16,'(A)')'# bctype                                                : Boundary condition definition '
      write(16,'(A)')'1'
      write(16,'(A)')'1'
      write(16,'(A)')'41, 1.d0'
      write(16,'(A)')'41, 1.d0'
      write(16,'(A)')'1'
      write(16,'(A)')'1'
      write(16,*)
      write(16,'(A)')'# ninit                                                 : Initial method'
      write(16,'(A)')'1'
      write(16,*)
      write(16,'(A)')'# spg_imin,spg_imax,spg_jmin,spg_jmax,spg_kmin,spg_kmax : Sponge layer range'
      write(16,'(A)')'0, 0, 0, 0, 0, 0'
      write(16,*)
      write(16,'(A)')'# gridfile                                              : grid'
      write(16,'(A)')'./datin/grid.chl'
      write(16,*)
      write(16,*)
      write(16,*)
      write(16,*)'########################################################################'
      write(16,*)'# bctype                                                               #'
      write(16,*)'#   1 : periodic bc,     nothing will be done.                         #'
      write(16,*)'#  41 : isothermal wall, wall temperature input.                       #'
      write(16,*)'########################################################################'
      close(16)
      print*,' << ',trim(folder),'/datin/input.chl'
      !
      call gridchannel(128,128,128,trim(folder))
      !
      open(16,file=trim(folder)//'/datin/contr.dat',form='formatted')
      write(16,'(A)')'############################################################'
      write(16,'(A)')'#               control file for ASTRR code                #'
      write(16,'(A)')'############################################################'
      write(16,*)
      write(16,'(A)')'# lwrite,lavg                : switch for IO and on-fly sta'
      write(16,'(A)')'f,f'
      write(16,*)
      write(16,'(A)')'# maxstep,nwrite,nlstep,navg : parameters'
      write(16,'(A)')'1000,100,1,20'
      write(16,*)
      write(16,'(A)')'# deltat                    : time step'
      write(16,'(A)')'1.d-3'
      close(16)
      print*,' << ',trim(folder),'/datin/contr.dat'
      !
      call system('cp -v ./bin/astr ./'//trim(folder))
      !
      print*,' ** An example case is generated.'
      print*,' ** you can now run a simulation of channel flow by : '
      print*,' mpirun -np 8 ./astr run -input datin/input.chl'
      !
    endif
    !
  end subroutine examplegen
  !+-------------------------------------------------------------------+
  !| The end of the subroutine examplegen.                             |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to generate an example grid.              |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 28-05-2021  | Created by J. Fang @ Warrington                     |
  !+-------------------------------------------------------------------+
  subroutine gridchannel(im,jm,km,folder)
    !
    use tecio
    use hdf5io
    use constdef
    use commfunc, only : argtanh
    !
    integer,intent(in) :: im,jm,km
    character(len=*),intent(in) :: folder
    !
    real(8),allocatable,dimension(:,:,:) :: x,y,z
    integer :: npa,n,i,j,k,jmm
    real(8) :: varc,var1,var2,Retau,dx
    !
    allocate(x(0:im,0:jm,0:km),y(0:im,0:jm,0:km),z(0:im,0:jm,0:km))
    !
    Retau=185.d0
    jmm=jm/2
    !
    varc=1.075d0
    do j=0,jm
      !
      var1=argtanh(1.d0/varc)
      var2=2.d0*(j)/(jm)*1.d0-1.d0
      !
      y(0,j,0)=1.d0*(1.d0+varc*dtanh(var1*var2))
      !
    end do
    !
    do k=0,km
    do j=0,jm
    do i=0,im
      x(i,j,k)=2.d0*pi/im*i
      y(i,j,k)=y(0,j,0)
      !
      if(km==0) then
        z(i,j,k)=0.d0
      else
        z(i,j,k)=pi/km*k
      endif
    end do
    end do
    end do
    print*,' ** y1+=',y(0,1,0)*Retau
    print*,' ** ym+=',(y(0,jmm,0)-y(0,jmm-1,0))*Retau
    print*,' ** dx+=',(x(1,0,0)-x(0,0,0))*Retau,(x(1,0,0)-x(0,0,0))
    print*,' ** dz+=',(z(0,0,1)-z(0,0,0))*Retau,(z(0,0,1)-z(0,0,0))
    print*,' ** lx+=',(x(im,0,0)-x(0,0,0))*Retau
    print*,' ** lz+=',(z(0,0,km)-z(0,0,0))*Retau
    !
    open(18,file='dy.dat')
    write(18,*)(j*1.d0,y(0,j,0),y(0,j,0)-y(0,j-1,0),j=1,jm)
    close(18)
    print*,' << dy.dat'
    !
    ! call writetecbin('Testout/tecgrid.plt',x,'x',y,'y',z,'z')
    !
    call h5srite(x,'x',folder//'/datin/grid.chl')
    call h5srite(y,'y',folder//'/datin/grid.chl')
    call h5srite(z,'z',folder//'/datin/grid.chl')
    !
    !
    deallocate(x,y,z)
    !
  end subroutine gridchannel
  !+-------------------------------------------------------------------+
  !| The end of the subroutine gridchannel.                            |
  !+-------------------------------------------------------------------+
  !
end module pp
!+---------------------------------------------------------------------+
!| The end of the module pp                                            |
!+---------------------------------------------------------------------+