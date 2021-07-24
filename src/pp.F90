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
  use constdef
  use stlaio,  only: get_unit
  !
  implicit none
  !
  contains
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is the entrance of post/pre-process.              |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 05-07-2021  | Created by J. Fang @ Warrington                     |
  !+-------------------------------------------------------------------+
  subroutine ppentrance
    !
    use cmdefne
    !
    ! local data
    character(len=64) :: cmd,casefolder
    !
    call readkeyboad(cmd=cmd,name1=casefolder)
    !
    if(trim(cmd)=='init') then
      call examplegen(trim(casefolder))
    elseif(trim(cmd)=='solid') then
      call solidpp(trim(casefolder))
    endif
    ! generate an example channel flow case
    ! 
    !
  end subroutine ppentrance
  !+-------------------------------------------------------------------+
  !| The end of the subroutine preprocess.                             |
  !+-------------------------------------------------------------------+
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
    integer :: fh
    !
    if(mpirank==0) then
      !
      print*,' ** Generating an example case.'
      !
      call system('mkdir '//trim(folder))
      call system('mkdir '//trim(folder)//'/datin')
      !
      fh=get_unit()
      !
      open(fh,file=trim(folder)//'/datin/input.chl',form='formatted')
      write(fh,'(A)')'########################################################################'
      write(fh,'(A)')'#                     input file of ASTR code                          #'
      write(fh,'(A)')'########################################################################'
      write(fh,*)
      write(fh,'(A)')'# flowtype                                              : The type of flow problem'
      write(fh,'(A)')'channel'
      write(fh,*)
      write(fh,'(A)')'# im,jm,km                                              : The size of grid.'
      write(fh,'(A)')'128,128,128'
      write(fh,*)
      write(fh,'(A)')'# lihomo,ljhomo,lkhomo                                  : The homogeneous directions'
      write(fh,'(A)')'t,f,t'
      write(fh,*)
      write(fh,'(A)')'# nondimen,diffterm,lfilter,lreadgrid,lfftz,limmbou     : Parameters'
      write(fh,'(A)')'t,t,t,t,f,f'
      write(fh,*)
      write(fh,'(A)')'# lrestar                                               : start mode'
      write(fh,'(A)')'f'
      write(fh,*)
      write(fh,'(A)')'# alfa_filter, kcutoff                                  : Filter parameters'
      write(fh,'(A)')'0.49d0, 48'
      write(fh,*)
      write(fh,'(A)')'# ref_t,reynolds,mach                                   : Reference variables'
      write(fh,'(A)')'273.15d0,  3000.d0,  0.5d0 '
      write(fh,*)
      write(fh,'(A)')'# conschm,difschm,rkscheme                              : Numerical scheme'
      write(fh,'(A)')'642c, 642c, rk4 '
      write(fh,*)
      write(fh,'(A)')'# recon_schem, lchardecomp                              : Parameters for upwind-biased scheme'
      write(fh,'(A)')'0, f '
      write(fh,*)
      write(fh,'(A)')'# num_species : number of species'
      write(fh,'(A)')'1'
      write(fh,*)
      write(fh,'(A)')'# bctype                                                : Boundary condition definition '
      write(fh,'(A)')'1'
      write(fh,'(A)')'1'
      write(fh,'(A)')'41, 1.d0'
      write(fh,'(A)')'41, 1.d0'
      write(fh,'(A)')'1'
      write(fh,'(A)')'1'
      write(fh,*)
      write(fh,'(A)')'# ninit                                                 : Initial method'
      write(fh,'(A)')'1'
      write(fh,*)
      write(fh,'(A)')'# spg_imin,spg_imax,spg_jmin,spg_jmax,spg_kmin,spg_kmax : Sponge layer range'
      write(fh,'(A)')'0, 0, 0, 0, 0, 0'
      write(fh,*)
      write(fh,'(A)')'# gridfile                                              : grid'
      write(fh,'(A)')'./datin/grid.chl'
      write(fh,*)
      write(fh,*)
      write(fh,*)
      write(fh,*)'########################################################################'
      write(fh,*)'# bctype                                                               #'
      write(fh,*)'#   1 : periodic bc,     nothing will be done.                         #'
      write(fh,*)'#  41 : isothermal wall, wall temperature input.                       #'
      write(fh,*)'########################################################################'
      close(fh)
      print*,' << ',trim(folder),'/datin/input.chl'
      !
      call gridchannel(128,128,128,trim(folder))
      !
      open(fh,file=trim(folder)//'/datin/contr.dat',form='formatted')
      write(fh,'(A)')'############################################################'
      write(fh,'(A)')'#               control file for ASTRR code                #'
      write(fh,'(A)')'############################################################'
      write(fh,*)
      write(fh,'(A)')'# lwrite,lavg                : switch for IO and on-fly sta'
      write(fh,'(A)')'f,f'
      write(fh,*)
      write(fh,'(A)')'# maxstep,nwrite,nlstep,navg : parameters'
      write(fh,'(A)')'1000,100,1,20'
      write(fh,*)
      write(fh,'(A)')'# deltat                    : time step'
      write(fh,'(A)')'1.d-3'
      close(fh)
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
    integer :: npa,n,i,j,k,jmm,fh
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
    fh=get_unit()
    open(fh,file='dy.dat')
    write(fh,*)(j*1.d0,y(0,j,0),y(0,j,0)-y(0,j-1,0),j=1,jm)
    close(fh)
    print*,' << dy.dat'
    !
    ! call writetecbin('Testout/tecgrid.plt',x,'x',y,'y',z,'z')
    !
    call h5srite(x,'x',folder//'/datin/grid.chl',newfile=.true.)
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
  !+-------------------------------------------------------------------+
  !| This subroutine is to process solid file to improve its accuracy. |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 05-07-2021  | Created by J. Fang @ Warrington                     |
  !+-------------------------------------------------------------------+
  subroutine solidpp(inputfile)
    !
    use commtype, only : solid,triangle
    use readwrite,only : readsolid
    use commvar,  only : nsolid,immbody
    use geom,     only : solidrange,solidresc,solidshif,solidimpro
    use tecio,     only : tecsolid
    use stlaio,    only : stla_write
    !
    ! arguments
    character(len=*),intent(in) :: inputfile
    !
    ! local data
    integer :: js
    !
    ! call readsolid(inputfile)
    call solidgen_circle
    ! call solidgen_cub
    !
    do js=1,nsolid
      call solidrange(immbody(js))
    enddo
    !
    ! call solidresc(immbody(1),0.025d0)
    !
    call solidshif(immbody(1),x=5.d0-immbody(1)%xcen(1),  &
                              y=5.d0-immbody(1)%xcen(2),  &
                              z=0.d0-immbody(1)%xcen(3))
    !
    call solidrange(immbody(1))
    !
    ! do js=1,nsolid
    !   call solidimpro(immbody(js))
    ! enddo
    !
    call tecsolid('tecsolid.plt',immbody)
    !
    call stla_write('solid_in.stl',immbody)
    !
  end subroutine solidpp
  !+-------------------------------------------------------------------+
  !| The end of the subroutine solidpp.                                |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to generate a solid usign stl file format.     |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 05-07-2021  | Created by J. Fang @ Warrington                     |
  !+-------------------------------------------------------------------+
  subroutine solidgen_cub
    !
    use commtype, only : solid,triangle
    use readwrite,only : readsolid
    use commvar,  only : nsolid,immbody
    use geom,     only : solidrange,solidresc,solidshif,solidimpro
    use tecio,     only : tecsolid
    use stlaio,    only : stla_write
    use commfunc,  only : cross_product
    !
    ! local data
    integer :: i,j,k,jf,km,jm,im,nface
    real(8) :: dthe,dphi,theter1,theter2,phi1,phi2,radi, &
               x1(3),x2(3),x3(3),x4(3),norm1(3),var1
    real(8) :: epsilon
    type(triangle),allocatable :: tempface(:)
    !
    epsilon=1.d-12
    !
    print*,' ** generating solid'
    !
    nsolid=1
    allocate(immbody(nsolid))
    immbody(1)%name='cube'
    !
    allocate(tempface(12))
    !
    x1(1)=0.d0
    x1(2)=0.d0
    x1(3)=0.d0
    !
    x2(1)=1.d0
    x2(2)=0.d0
    x2(3)=0.d0
    !
    x3(1)=1.d0
    x3(2)=1.d0
    x3(3)=0.d0
    !
    x4(1)=0.d0
    x4(2)=1.d0
    x4(3)=0.d0
    !
    nface=nface+1
    tempface(nface)%a=x1
    tempface(nface)%b=x2
    tempface(nface)%c=x3
    tempface(nface)%normdir=(/0.d0,0.d0,-1.d0/)
    !
    nface=nface+1
    tempface(nface)%a=x3
    tempface(nface)%b=x4
    tempface(nface)%c=x1
    tempface(nface)%normdir=(/0.d0,0.d0,-1.d0/)
    !
    x1(1)=1.d0
    x1(2)=0.d0
    x1(3)=0.d0
    !
    x2(1)=1.d0
    x2(2)=0.d0
    x2(3)=1.d0
    !
    x3(1)=1.d0
    x3(2)=1.d0
    x3(3)=1.d0
    !
    x4(1)=1.d0
    x4(2)=1.d0
    x4(3)=0.d0
    !
    nface=nface+1
    tempface(nface)%a=x1
    tempface(nface)%b=x2
    tempface(nface)%c=x3
    tempface(nface)%normdir=(/1.d0,0.d0,0.d0/)
    !
    nface=nface+1
    tempface(nface)%a=x3
    tempface(nface)%b=x4
    tempface(nface)%c=x1
    tempface(nface)%normdir=(/1.d0,0.d0,0.d0/)
    !
    x1(1)=0.d0
    x1(2)=0.d0
    x1(3)=0.d0
    !
    x2(1)=0.d0
    x2(2)=0.d0
    x2(3)=1.d0
    !
    x3(1)=0.d0
    x3(2)=1.d0
    x3(3)=1.d0
    !
    x4(1)=0.d0
    x4(2)=1.d0
    x4(3)=0.d0
    !
    nface=nface+1
    tempface(nface)%a=x1
    tempface(nface)%b=x2
    tempface(nface)%c=x3
    tempface(nface)%normdir=(/-1.d0,0.d0,0.d0/)
    !
    nface=nface+1
    tempface(nface)%a=x3
    tempface(nface)%b=x4
    tempface(nface)%c=x1
    tempface(nface)%normdir=(/-1.d0,0.d0,0.d0/)
    !
    x1(1)=0.d0
    x1(2)=0.d0
    x1(3)=1.d0
    !
    x2(1)=1.d0
    x2(2)=0.d0
    x2(3)=1.d0
    !
    x3(1)=1.d0
    x3(2)=1.d0
    x3(3)=1.d0
    !
    x4(1)=0.d0
    x4(2)=1.d0
    x4(3)=1.d0
    !
    nface=nface+1
    tempface(nface)%a=x1
    tempface(nface)%b=x2
    tempface(nface)%c=x3
    tempface(nface)%normdir=(/0.d0,0.d0,1.d0/)
    !
    nface=nface+1
    tempface(nface)%a=x3
    tempface(nface)%b=x4
    tempface(nface)%c=x1
    tempface(nface)%normdir=(/0.d0,0.d0,1.d0/)
    !
    x1(1)=0.d0
    x1(2)=0.d0
    x1(3)=0.d0
    !
    x2(1)=1.d0
    x2(2)=0.d0
    x2(3)=0.d0
    !
    x3(1)=1.d0
    x3(2)=0.d0
    x3(3)=1.d0
    !
    x4(1)=0.d0
    x4(2)=0.d0
    x4(3)=1.d0
    !
    nface=nface+1
    tempface(nface)%a=x1
    tempface(nface)%b=x2
    tempface(nface)%c=x3
    tempface(nface)%normdir=(/0.d0,-1.d0,0.d0/)
    !
    nface=nface+1
    tempface(nface)%a=x3
    tempface(nface)%b=x4
    tempface(nface)%c=x1
    tempface(nface)%normdir=(/0.d0,-1.d0,0.d0/)
    !
    x1(1)=0.d0
    x1(2)=1.d0
    x1(3)=0.d0
    !
    x2(1)=1.d0
    x2(2)=1.d0
    x2(3)=0.d0
    !
    x3(1)=1.d0
    x3(2)=1.d0
    x3(3)=1.d0
    !
    x4(1)=0.d0
    x4(2)=1.d0
    x4(3)=1.d0
    !
    nface=nface+1
    tempface(nface)%a=x1
    tempface(nface)%b=x2
    tempface(nface)%c=x3
    tempface(nface)%normdir=(/0.d0,1.d0,0.d0/)
    !
    nface=nface+1
    tempface(nface)%a=x3
    tempface(nface)%b=x4
    tempface(nface)%c=x1
    
    tempface(nface)%normdir=(/0.d0,1.d0,0.d0/)
    !
    immbody(1)%num_face=nface
    call immbody(1)%alloface()
    immbody(1)%face(1:nface)=tempface(1:nface)
    !
    call tecsolid('tecsolid.plt',immbody)
    !
    print*,' ** solid generated'
    !
  end subroutine solidgen_cub
  !+-------------------------------------------------------------------+
  !| The end of the subroutine solidgen_cub.                           |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to generate a solid usign stl file format.     |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 05-07-2021  | Created by J. Fang @ Warrington                     |
  !+-------------------------------------------------------------------+
  subroutine solidgen_sphere
    !
    use commtype, only : solid,triangle
    use readwrite,only : readsolid
    use commvar,  only : nsolid,immbody
    use geom,     only : solidrange,solidresc,solidshif,solidimpro
    use tecio,     only : tecsolid
    use stlaio,    only : stla_write
    use commfunc,  only : cross_product
    !
    ! local data
    integer :: i,j,k,jf,km,jm,im,nface
    real(8) :: dthe,dphi,theter1,theter2,phi1,phi2,radi, &
               x1(3),x2(3),x3(3),x4(3),norm1(3),var1
    real(8) :: epsilon
    type(triangle),allocatable :: tempface(:)
    !
    epsilon=1.d-12
    !
    print*,' ** generating solid'
    !
    nsolid=1
    allocate(immbody(nsolid))
    immbody(1)%name='sphere'
    !
    km=16
    jm=32
    !
    allocate(tempface(km*jm*2))
    !
    dthe=2.d0*pi/km
    dphi=pi/jm
    radi=0.5d0
    !
    nface=0
    !
    do k=0,km-1
    do j=0,jm-1
      !
      theter1=dthe*k
      theter2=dthe*(k+1)
      phi1=dphi*j
      phi2=dphi*(j+1)
      !
      x1(1)=radi*sin(phi1)*cos(theter1)
      x1(2)=radi*sin(phi1)*sin(theter1)
      x1(3)=radi*cos(phi1)
      !
      x2(1)=radi*sin(phi2)*cos(theter1)
      x2(2)=radi*sin(phi2)*sin(theter1)
      x2(3)=radi*cos(phi2)
      !
      x3(1)=radi*sin(phi2)*cos(theter2)
      x3(2)=radi*sin(phi2)*sin(theter2)
      x3(3)=radi*cos(phi2)
      !
      x4(1)=radi*sin(phi1)*cos(theter2)
      x4(2)=radi*sin(phi1)*sin(theter2)
      x4(3)=radi*cos(phi1)
      !
      if(abs(sin(phi1))<epsilon) then
        !
        nface=nface+1
        tempface(nface)%a=x1
        tempface(nface)%b=x2
        tempface(nface)%c=x3
        !
      elseif(abs(sin(phi2))<epsilon) then
        !
        nface=nface+1
        tempface(nface)%a=x3
        tempface(nface)%b=x4
        tempface(nface)%c=x1
        !
      else
        !
        nface=nface+1
        tempface(nface)%a=x1
        tempface(nface)%b=x2
        tempface(nface)%c=x3
        !
        nface=nface+1
        tempface(nface)%a=x3
        tempface(nface)%b=x4
        tempface(nface)%c=x1
        !
      endif
      !
    enddo
    enddo
    !
    immbody(1)%num_face=nface
    call immbody(1)%alloface()
    immbody(1)%face(1:nface)=tempface(1:nface)
    !
    do jf=1,immbody(1)%num_face
      !
      norm1=cross_product(immbody(1)%face(jf)%b-immbody(1)%face(jf)%a,         &
                          immbody(1)%face(jf)%c-immbody(1)%face(jf)%a )
      !
      var1=sqrt(norm1(1)**2+norm1(2)**2+norm1(3)**2)
      !
      if(abs(var1)<epsilon) then
        !
        norm1=cross_product(immbody(1)%face(jf)%a-immbody(1)%face(jf)%b,         &
                            immbody(1)%face(jf)%c-immbody(1)%face(jf)%b )
        !
        var1=sqrt(norm1(1)**2+norm1(2)**2+norm1(3)**2)
        !
      endif
      !
      if(abs(var1)<epsilon) then
        !
        norm1=cross_product(immbody(1)%face(jf)%a-immbody(1)%face(jf)%c,         &
                            immbody(1)%face(jf)%b-immbody(1)%face(jf)%c )
        !
        var1=sqrt(norm1(1)**2+norm1(2)**2+norm1(3)**2)
        !
      endif
      !
      immbody(1)%face(jf)%normdir=norm1/var1
      !
      if(isnan(immbody(1)%face(jf)%normdir(1))) then
        print*,jf
        print*,immbody(1)%face(jf)%a
        print*,immbody(1)%face(jf)%b
        print*,immbody(1)%face(jf)%c
        print*,norm1,'|',var1
        stop
      endif
      !
    enddo
    !
    call tecsolid('tecsolid.plt',immbody)
    !
    print*,' ** solid generated'
    !
  end subroutine solidgen_sphere
  !+-------------------------------------------------------------------+
  !| The end of the subroutine solidgen_sphere.                        |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to generate a solid usign stl file format.     |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 06-07-2021  | Created by J. Fang @ Warrington                     |
  !+-------------------------------------------------------------------+
  subroutine solidgen_circle
    !
    use commtype, only : solid,triangle
    use readwrite,only : readsolid
    use commvar,  only : nsolid,immbody
    use geom,     only : solidrange,solidresc,solidshif,solidimpro
    use tecio,     only : tecsolid
    use stlaio,    only : stla_write
    use commfunc,  only : cross_product
    !
    ! local data
    integer :: i,j,k,jf,km,jm,im,nface
    real(8) :: dthe,dphi,theter1,theter2,phi1,phi2,radi, &
               x1(3),x2(3),x3(3),x4(3),x5(3),x6(3),norm1(3),var1
    real(8) :: epsilon
    type(triangle),allocatable :: tempface(:)
    !
    epsilon=1.d-12
    !
    print*,' ** generating solid'
    !
    nsolid=1
    allocate(immbody(nsolid))
    immbody(1)%name='circle'
    !
    km=360
    !
    allocate(tempface(km*4))
    !
    dthe=2.d0*pi/km
    radi=0.5d0
    !
    nface=0
    !
    do k=0,km-1
      !
      theter1=dthe*k
      theter2=dthe*(k+1)
      !
      x1(1)=radi*cos(theter1)
      x1(2)=radi*sin(theter1)
      x1(3)=-0.01d0
      !
      x2(1)=radi*cos(theter2)
      x2(2)=radi*sin(theter2)
      x2(3)=-0.01d0
      !
      x3(1)=radi*cos(theter2)
      x3(2)=radi*sin(theter2)
      x3(3)=0.01d0
      !
      x4(1)=radi*cos(theter1)
      x4(2)=radi*sin(theter1)
      x4(3)=0.01d0
      !
      x5(1)=0.d0
      x5(2)=0.d0
      x5(3)=-0.01d0
      !
      x6(1)=0.d0
      x6(2)=0.d0
      x6(3)=0.01d0
      !
      nface=nface+1
      tempface(nface)%a=x1
      tempface(nface)%b=x2
      tempface(nface)%c=x3
      !
      nface=nface+1
      tempface(nface)%a=x3
      tempface(nface)%b=x4
      tempface(nface)%c=x1
      !
      ! nface=nface+1
      ! tempface(nface)%a=x1
      ! tempface(nface)%b=x2
      ! tempface(nface)%c=x5
      ! !
      ! nface=nface+1
      ! tempface(nface)%a=x3
      ! tempface(nface)%b=x4
      ! tempface(nface)%c=x6
      !
    enddo
    !
    immbody(1)%num_face=nface
    call immbody(1)%alloface()
    immbody(1)%face(1:nface)=tempface(1:nface)
    !
    do jf=1,immbody(1)%num_face
      !
      norm1=cross_product(immbody(1)%face(jf)%b-immbody(1)%face(jf)%a,         &
                          immbody(1)%face(jf)%c-immbody(1)%face(jf)%a )
      !
      var1=sqrt(norm1(1)**2+norm1(2)**2+norm1(3)**2)
      !
      if(abs(var1)<epsilon) then
        !
        norm1=cross_product(immbody(1)%face(jf)%a-immbody(1)%face(jf)%b,         &
                            immbody(1)%face(jf)%c-immbody(1)%face(jf)%b )
        !
        var1=sqrt(norm1(1)**2+norm1(2)**2+norm1(3)**2)
        !
      endif
      !
      if(abs(var1)<epsilon) then
        !
        norm1=cross_product(immbody(1)%face(jf)%a-immbody(1)%face(jf)%c,         &
                            immbody(1)%face(jf)%b-immbody(1)%face(jf)%c )
        !
        var1=sqrt(norm1(1)**2+norm1(2)**2+norm1(3)**2)
        !
      endif
      !
      immbody(1)%face(jf)%normdir=norm1/var1
      !
      if(isnan(immbody(1)%face(jf)%normdir(1))) then
        print*,jf
        print*,immbody(1)%face(jf)%a
        print*,immbody(1)%face(jf)%b
        print*,immbody(1)%face(jf)%c
        print*,norm1,'|',var1
        stop
      endif
      !
    enddo
    !
    call tecsolid('tecsolid.plt',immbody)
    !
    print*,' ** solid generated'
    !
  end subroutine solidgen_circle
  !+-------------------------------------------------------------------+
  !| The end of the subroutine solidgen_circle.                        |
  !+-------------------------------------------------------------------+
  !
end module pp
!+---------------------------------------------------------------------+
!| The end of the module pp                                            |
!+---------------------------------------------------------------------+