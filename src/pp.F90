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
    character(len=64) :: cmd,casefolder,inputfile,outputfile
    !
    call readkeyboad(cmd)
    print*,' ** pp command: ',cmd
    !
    if(trim(cmd)=='init') then
      !
      call readkeyboad(casefolder)
      call examplegen(trim(casefolder))
      ! generate an example channel flow case
      !
    elseif(trim(cmd)=='solid') then
      call solidpp
    elseif(trim(cmd)=='datacon') then
      !
      call readkeyboad(inputfile)
      !
      if(trim(inputfile)=='all') then
        call stream2struc('outdat/flowfield')
        call stream2struc('outdat/meanflow')
        call stream2struc('outdat/2ndsta')
        call stream2struc('outdat/3rdsta')
        call stream2struc('outdat/budget')
      else
        call stream2struc(trim(inputfile))
      endif
      !
    elseif(trim(cmd)=='parinfo') then
      !
      call parallelifogen
      !
    else
      stop ' !! pp command not defined. @ ppentrance'
    endif
    ! 
    !
  end subroutine ppentrance
  !+-------------------------------------------------------------------+
  !| The end of the subroutine preprocess.                             |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to generate parallel.info file.                |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 04-04-2022  | Created by J. Fang @ Warrington                     |
  !+-------------------------------------------------------------------+
  subroutine parallelifogen
    !
    use readwrite, only: readinput
    use parallel,  only: mpisizedis,parapp,mpisize,mpirankmax
    !
    call readinput
    !
    mpisize=128*128
    mpirankmax=mpisize-1
    !
    call mpisizedis
    !
    call parapp
    !
  end subroutine parallelifogen
  !+-------------------------------------------------------------------+
  !| The end of the subroutine parallelifogen.                         |
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
      write(fh,'(A)')'# recon_schem, lchardecomp,bfacmpld                     : Parameters for upwind-biased scheme'
      write(fh,'(A)')'0, f, 0.3d0 '
      write(fh,*)
      write(fh,'(A)')'# num_species                                           : number of species'
      write(fh,'(A)')'1'
      write(fh,*)
      write(fh,'(A)')'# turbmode                                              : turbulence model'
      write(fh,'(A)')'none'
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
      call gridchannel(320,128,128,trim(folder))
      !
      open(fh,file=trim(folder)//'/datin/controller',form='formatted')
      write(fh,'(A)')'############################################################'
      write(fh,'(A)')'#               control file for ASTRR code                #'
      write(fh,'(A)')'############################################################'
      write(fh,*)
      write(fh,'(A)')'# lwsequ,lwslic,lavg'
      write(fh,'(A)')'f,f'
      write(fh,*)
      write(fh,'(A)')'# maxstep,nwrite,ninst,nlstep,navg'
      write(fh,'(A)')'1000,100,20,1,20'
      write(fh,*)
      write(fh,'(A)')'# deltat '
      write(fh,'(A)')'1.d-3'
      write(fh,'(A)')'+----------------------------------------------------------+'
      write(fh,'(A)')'| This file will be read each time after checkpoint        |'
      write(fh,'(A)')'+--------------+-------------------------------------------+'
      write(fh,'(A)')'|       lwrite | to write a sequence of flowfield files.   |'
      write(fh,'(A)')'|       lwslic | to write a sequence of slice cut files.   |'
      write(fh,'(A)')'|         lavg | to conduct on-fly stastistics.            |'
      write(fh,'(A)')'+--------------+-------------------------------------------+'
      write(fh,'(A)')'|      maxstep | max step to run.                          |'
      write(fh,'(A)')'|       nwrite | frequency of dumping checkpoint.          |'
      write(fh,'(A)')'|        ninst | frequency of writing slice.               |'
      write(fh,'(A)')'|       nlstep | frequency of listing computing state.     |'
      write(fh,'(A)')'|         navg | frequency of calculating statistics.      |'
      write(fh,'(A)')'+--------------+-------------------------------------------+'
      write(fh,'(A)')'|       deltat | time step.                                |'
      write(fh,'(A)')'+--------------+-------------------------------------------+'
      close(fh)
      print*,' << ',trim(folder),'/datin/controller'
      !
      call system('cp -v ./bin/astr ./'//trim(folder))
      !
      print*,' ** An example case is generated.'
      print*,' ** you can now run a simulation of channel flow by : '
      print*,' mpirun -np 8 ./astr run datin/input.chl'
      !
    endif
    !
  end subroutine examplegen
  !+-------------------------------------------------------------------+
  !| The end of the subroutine examplegen.                             |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to convert a streammed data to a          |
  !| structured data.                                                  |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 31-03-2022  | Created by J. Fang @ Warrington                     |
  !+-------------------------------------------------------------------+
  subroutine stream2struc(fname)
    !
    use hdf5io
    !
    character(len=*),intent(in) :: fname
    !
    integer :: fh,isize,jsize,ksize,rankmax,jrank,vai1,irp,jrp,krp,    &
               asizemax,ima,jma,kma,i,j,k,n,nstep
    real(8) :: time
    integer,allocatable :: im(:),jm(:),km(:),i0(:),j0(:),k0(:),        &
                           offset(:),asize(:)
    !
    real(8),allocatable,dimension(:) :: data_1d
    real(8),allocatable,dimension(:,:,:) :: data_3d
    !
    fh=get_unit()
    open(fh,file='datin/parallel.info',form='formatted',action='read')
    read(fh,*)
    read(fh,*)isize,jsize,ksize
    rankmax=isize*jsize*ksize-1
    print*,' ** rankmax=',rankmax
    allocate(im(0:rankmax),jm(0:rankmax),km(0:rankmax),                &
             i0(0:rankmax),j0(0:rankmax),k0(0:rankmax),                &
             offset(0:rankmax),asize(0:rankmax))
    read(fh,*)
    do jrank=0,rankmax
      read(fh,*)vai1,irp,jrp,krp,im(jrank),jm(jrank),km(jrank),        &
                i0(jrank),j0(jrank),k0(jrank)
    enddo
    close(fh)
    print*,' >> datin/parallel.info'
    !
    do jrank=0,rankmax
      asize(jrank)=(im(jrank)+1)*(jm(jrank)+1)*(km(jrank)+1)
    enddo
    !
    offset(0)=0
    do jrank=1,rankmax
      offset(jrank)=offset(jrank-1)+asize(jrank-1)
    enddo
    !
    asizemax=offset(rankmax)+asize(rankmax)
    !
    ima=i0(rankmax)+im(rankmax)
    jma=j0(rankmax)+jm(rankmax)
    kma=k0(rankmax)+km(rankmax)
    !
    print*,' **     stream data size:',asizemax
    print*,' ** structured data size:',ima,jma,kma
    ! !
    allocate( data_1d(1:asizemax),data_3d(0:ima,0:jma,0:kma) )
    !
    call h5sread(varname='nstep',var=nstep,filename=fname//'.s5')
    call h5srite(var=nstep,varname='nstep',filename=fname//'.h5',newfile=.true.)
    !
    if(fname=='outdat/flowfield') then
      !
      call h5sread(varname='time',var=time,filename=fname//'.s5')
      call h5srite(var=time,varname='time',filename=fname//'.h5')
      !
      call data3d_rcw(varname='ro')
      call data3d_rcw(varname='u1')
      call data3d_rcw(varname='u2')
      call data3d_rcw(varname='u3')
      call data3d_rcw(varname='p')
      call data3d_rcw(varname='t')
      !
    elseif(fname=='outdat/meanflow') then
      !
      call h5sread(varname='nsamples',var=time,filename=fname//'.s5')
      call h5srite(var=time,varname='nsamples',filename=fname//'.h5')
      !
      call h5sread(varname='nstep_sbeg',var=time,filename=fname//'.s5')
      call h5srite(var=time,varname='nstep_sbeg',filename=fname//'.h5')
      !
      call h5sread(varname='time_sbeg',var=time,filename=fname//'.s5')
      call h5srite(var=time,varname='time_sbeg',filename=fname//'.h5')
      !
      call data3d_rcw(varname='rom')
      call data3d_rcw(varname='u1m')
      call data3d_rcw(varname='u2m')
      call data3d_rcw(varname='u3m')
      call data3d_rcw(varname='pm')
      call data3d_rcw(varname='tm')
      !
    elseif(fname=='outdat/2ndsta') then
      !
      call h5sread(varname='nsamples',var=time,filename=fname//'.s5')
      call h5srite(var=time,varname='nsamples',filename=fname//'.h5')
      !
      call h5sread(varname='nstep_sbeg',var=time,filename=fname//'.s5')
      call h5srite(var=time,varname='nstep_sbeg',filename=fname//'.h5')
      !
      call h5sread(varname='time_sbeg',var=time,filename=fname//'.s5')
      call h5srite(var=time,varname='time_sbeg',filename=fname//'.h5')
      !
      call data3d_rcw(varname= 'pp')
      call data3d_rcw(varname= 'tt')
      call data3d_rcw(varname='tu1')
      call data3d_rcw(varname='tu2')
      call data3d_rcw(varname='tu3')
      call data3d_rcw(varname='u11')
      call data3d_rcw(varname='u12')
      call data3d_rcw(varname='u13')
      call data3d_rcw(varname='u22')
      call data3d_rcw(varname='u23')
      call data3d_rcw(varname='u33')
      !
    elseif(fname=='outdat/3rdsta') then
      !
      call h5sread(varname='nsamples',var=time,filename=fname//'.s5')
      call h5srite(var=time,varname='nsamples',filename=fname//'.h5')
      !
      call h5sread(varname='nstep_sbeg',var=time,filename=fname//'.s5')
      call h5srite(var=time,varname='nstep_sbeg',filename=fname//'.h5')
      !
      call h5sread(varname='time_sbeg',var=time,filename=fname//'.s5')
      call h5srite(var=time,varname='time_sbeg',filename=fname//'.h5')
      !
      call data3d_rcw(varname='u111')
      call data3d_rcw(varname='u112')
      call data3d_rcw(varname='u113')
      call data3d_rcw(varname='u122')
      call data3d_rcw(varname='u123')
      call data3d_rcw(varname='u133')
      call data3d_rcw(varname='u222')
      call data3d_rcw(varname='u223')
      call data3d_rcw(varname='u233')
      call data3d_rcw(varname='u333')
      !
    elseif(fname=='outdat/budget') then
      !
      call h5sread(varname='nsamples',var=time,filename=fname//'.s5')
      call h5srite(var=time,varname='nsamples',filename=fname//'.h5')
      !
      call h5sread(varname='nstep_sbeg',var=time,filename=fname//'.s5')
      call h5srite(var=time,varname='nstep_sbeg',filename=fname//'.h5')
      !
      call h5sread(varname='time_sbeg',var=time,filename=fname//'.s5')
      call h5srite(var=time,varname='time_sbeg',filename=fname//'.h5')
      !
      call data3d_rcw(varname= 'predil')
      call data3d_rcw(varname=    'pu1')
      call data3d_rcw(varname=    'pu2')
      call data3d_rcw(varname=    'pu3')
      call data3d_rcw(varname='sgmam11')
      call data3d_rcw(varname='sgmam12')
      call data3d_rcw(varname='sgmam13')
      call data3d_rcw(varname='sgmam22')
      call data3d_rcw(varname='sgmam23')
      call data3d_rcw(varname='sgmam33')
      call data3d_rcw(varname=  'u1rem')
      call data3d_rcw(varname=  'u2rem')
      call data3d_rcw(varname=  'u3rem')
      call data3d_rcw(varname='visdif1')
      call data3d_rcw(varname='visdif2')
      call data3d_rcw(varname='visdif3')
      !
    else
      stop ' !! file name error @ stream2struc'
    endif
    !
    contains
    !
    subroutine data3d_rcw(varname)
      !
      character(len=*),intent(in) :: varname
      !
      call h5sread(varname=varname,var= data_1d,dim=asizemax,filename=fname//'.s5')
      !
      write(*,'(A)',advance='no')'  ** data converting ... '
      do jrank=0,rankmax
        !
        n=offset(jrank)
        !
        do k=k0(jrank),k0(jrank)+km(jrank)
        do j=j0(jrank),j0(jrank)+jm(jrank)
        do i=i0(jrank),i0(jrank)+im(jrank)
          !
          n=n+1
          !
          data_3d(i,j,k)=data_1d(n)
          !
        enddo
        enddo
        enddo
        !
        write(*,'(1A1,A,I4,A,$)')char(13),'  ** data converting ... ', &
                                                100*jrank/rankmax,'  % '
        !
      enddo
      write(*,*)''
      !
      call h5srite(var=data_3d,varname=varname,filename=fname//'.h5')
      !
    end subroutine data3d_rcw
    !
  end subroutine stream2struc
  !+-------------------------------------------------------------------+
  !| The end of the subroutine stream2struc.                           |
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
      x(i,j,k)=4.d0*pi/im*i
      y(i,j,k)=y(0,j,0)
      !
      if(km==0) then
        z(i,j,k)=0.d0
      else
        z(i,j,k)=num4d3*pi/km*k
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
  subroutine solidpp
    !
    use commtype, only : solid,triangle
    use readwrite,only : readsolid
    use commvar,  only : nsolid,immbody
    use geom,     only : solidrange,solidresc,solidshif,solidimpro,solidrota
    use tecio,    only : tecsolid
    use stlaio,   only : stla_write
    use cmdefne,  only : readkeyboad
    !
    ! arguments
    !
    ! local data
    integer :: js
    character(len=64) :: inputfile
    character(len=4) :: cmd
    real(8) :: resc_fact
    real(8) :: rot_vec(3),rot_theta,shift_cor(3)
    !
    call readkeyboad(cmd)
    !
    print*,cmd
    !
    if(cmd=='sgen') then
      call solidgen_mvg
    elseif(cmd=='proc') then
      call readkeyboad(inputfile)
      !
      call readsolid(trim(inputfile))
      !
      ! resc_fact=0.015d0
      resc_fact=0.06667d0
      rot_vec=(/0.d0,1.d0,0.d0/)
      rot_theta=-90.d0
      shift_cor=(/5.d0,5.d0,0.d0/)
    else
      !
      print*,' !! ERROR 1, cmd not defined !!'
      print*,' ** cmd=',cmd
      print*,' ** inputfile',inputfile
      !
      stop
      !
    endif
    !
    ! 
    ! call solidgen_circle
    ! call solidgen_cub
    ! call solidgen_triagnle
    ! call solidgen_airfoil
    !
    ! do js=1,nsolid
    !   call solidrange(immbody(js))
    !   !
    !   call solidresc(immbody(js),resc_fact)
    !   call solidrota(immbody(js),rot_theta,rot_vec)
    !   call solidshif(immbody(js),x=shift_cor(1)-immbody(js)%xcen(1),  &
    !                              y=shift_cor(2)-immbody(js)%xcen(2),  &
    !                              z=shift_cor(3)-immbody(js)%xcen(3))
    !   !
    ! enddo
    !
    !
    !
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
  subroutine solidgen_airfoil
    !
    use commtype, only : solid,triangle
    use readwrite,only : readsolid
    use commvar,  only : nsolid,immbody
    use geom,     only : solidrange,solidresc,solidshif,solidimpro
    use tecio,     only : tecsolid
    use stlaio,    only : stla_write
    use commfunc,  only : cross_product,argtanh
    !
    ! local data
    integer :: i,j,k,jf,km,jm,im,nface
    real(8) :: dthe,dphi,theter1,theter2,phi1,phi2,radi, &
               x1(3),x2(3),x3(3),x4(3),x5(3),x6(3),norm1(3),var1,var2, &
               xc,yc
    real(8) :: epsilon,varc,theter
    integer :: map
    character(len=4) :: nacaname
    real(8),allocatable :: xap(:),yap(:),xin(:)
    type(triangle),allocatable :: tempface(:)
    !
    epsilon=1.d-12
    !
    print*,' ** generating solid'
    !
    map=1024
    nacaname='4412'
    allocate(xin(0:map),xap(0:map),yap(0:map))
    !
    varc=1.02d0
    do i=0,map/2
      !
      var1=argtanh(1.d0/varc)
      var2=2.d0*(i)/(map)*1.d0-1.d0
      !
      xin(i)=1.d0*(1.d0+varc*dtanh(var1*var2))
      call naca4digit(xin(i),xap(i),yap(i),nacaname,'upper')
      !
    end do
    !
    do i=map/2+1,map
      xin(i)=xin(map-i)
      call naca4digit(xin(i),xap(i),yap(i),nacaname,'lower')
    enddo
    !
    theter=-10.d0/180.d0*pi
    do i=0,map
      var1=xap(i)*cos(theter)-yap(i)*sin(theter)
      var2=xap(i)*sin(theter)+yap(i)*cos(theter)
      xap(i)=var1
      yap(i)=var2
    enddo
    !
    open(18,file='naca'//nacaname//'.dat')
    do i=0,map
      write(18,*)xap(i),yap(i)
    enddo
    close(18)
    print*,' << naca',nacaname,'.dat'
    !
    ! open(12,file='naca4412.dat')
    ! read(12,*)map
    ! allocate(xap(map),yap(map))
    ! do i=1,map
    !   read(12,*)xap(i),yap(i)
    ! enddo
    ! close(12)
    ! print*,' << naca4412.dat'
    !
    nface=0
    !
    xc=0.d0
    yc=0.d0
    do i=1,map-1
      xc=xc+xap(i)
      yc=yc+yap(i)
    enddo
    !
    xc=xc/dble(map-1)
    yc=yc/dble(map-1)
    !
    nsolid=1
    allocate(immbody(nsolid))
    immbody(1)%name='cube'
    !
    allocate(tempface(map*16))
    !
    do i=0,map-1
      x1(1)=xap(i)
      x1(2)=yap(i)
      x1(3)=0.d0
      !
      x2(1)=xap(i+1)
      x2(2)=yap(i+1)
      x2(3)=0.d0
      !
      x3(1)=xap(i+1)
      x3(2)=yap(i+1)
      x3(3)=1.d0
      !
      x4(1)=xap(i)
      x4(2)=yap(i)
      x4(3)=1.d0
      !
      x5(1)=xc
      x5(2)=yc
      x5(3)=0.d0
      !
      x6(1)=xc
      x6(2)=yc
      x6(3)=1.d0
      !
      nface=nface+1
      tempface(nface)%a=x1
      tempface(nface)%b=x2
      tempface(nface)%c=x3
      tempface(nface)%normdir=cross_product(x3-x1,x2-x1)
      !
      nface=nface+1
      tempface(nface)%a=x3
      tempface(nface)%b=x4
      tempface(nface)%c=x1
      tempface(nface)%normdir=cross_product(x1-x3,x4-x3)
      !
      nface=nface+1
      tempface(nface)%a=x1
      tempface(nface)%b=x2
      tempface(nface)%c=x5
      tempface(nface)%normdir=(/0.d0,0.d0,-1.d0/)
      !
      nface=nface+1
      tempface(nface)%a=x4
      tempface(nface)%b=x3
      tempface(nface)%c=x6
      tempface(nface)%normdir=(/0.d0,0.d0,1.d0/)
      !
    enddo
    !
    immbody(1)%num_face=nface
    call immbody(1)%alloface()
    immbody(1)%face(1:nface)=tempface(1:nface)
    !
    call tecsolid('tecsolid.plt',immbody)
    !
    print*,' ** solid generated'
    !
  end subroutine solidgen_airfoil
  !+-------------------------------------------------------------------+
  !| The end of the subroutine solidgen_airfoil.                       |
  !+-------------------------------------------------------------------+
  !
  
  !+-------------------------------------------------------------------+
  !| This subroutine is to generate a mvg solid usign stl file format. |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 22-10-2021  | Created by J. Fang @ Warrington                     |
  !+-------------------------------------------------------------------+
  subroutine solidgen_mvg
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
               x1(3),x2(3),x3(3),x4(3),norm1(3),var1,var2
    real(8) :: epsilon,delta,x0,lx,lz,h,sz,chord
    type(triangle),allocatable :: tempface(:)
    !
    epsilon=1.d-12
    delta=1.d0
    x0=10.d0
    lx=2.523597432d0
    lz=2.25d0
    chord=sqrt(lx**2+0.25d0*lz**2)
    sz=0.25d0
    h =0.3838d0
    nface=0
    !
    print*,' ** generating solid'
    !
    nsolid=2
    allocate(immbody(nsolid))
    immbody(1)%name='mvg1'
    !
    allocate(tempface(12))
    !
    x1(1)=x0
    x1(2)=0.d0
    x1(3)=0.5d0*sz
    !
    x2(1)=x0
    x2(2)=0.d0
    x2(3)=0.5d0*sz+lz
    !
    x3(1)=x0+lx
    x3(2)=0.d0
    x3(3)=0.5d0*sz+0.5d0*lz
    !
    x4(1)=x0+lx
    x4(2)=0.d0+h
    x4(3)=0.5d0*sz+0.5d0*lz
    !
    nface=nface+1
    tempface(nface)%a=x1
    tempface(nface)%b=x2
    tempface(nface)%c=x3
    tempface(nface)%normdir=(/0.d0,-1.d0,0.d0/)
    !
    var1=0.5d0*lz/chord
    var2=      lx/chord
    nface=nface+1
    tempface(nface)%a=x1
    tempface(nface)%b=x3
    tempface(nface)%c=x4
    tempface(nface)%normdir=(/var1,0.d0,-var2/)
    !
    var1=0.5d0*lz/chord
    var2=      lx/chord
    nface=nface+1
    tempface(nface)%a=x2
    tempface(nface)%b=x3
    tempface(nface)%c=x4
    tempface(nface)%normdir=(/var1,0.d0,var2/)
    !
    var1=h /sqrt(lx**2+h**2)
    var2=lx/sqrt(lx**2+h**2)
    nface=nface+1
    tempface(nface)%a=x1
    tempface(nface)%b=x2
    tempface(nface)%c=x4
    tempface(nface)%normdir=(/-var1,var2,0.d0/)
    !
    immbody(1)%num_face=nface
    call immbody(1)%alloface()
    immbody(1)%face(1:nface)=tempface(1:nface)
    !
    nface=0
    !
    x1(3)=x1(3)+sz+lz
    x2(3)=x2(3)+sz+lz
    x3(3)=x3(3)+sz+lz
    x4(3)=x4(3)+sz+lz
    !
    nface=nface+1
    tempface(nface)%a=x1
    tempface(nface)%b=x2
    tempface(nface)%c=x3
    tempface(nface)%normdir=(/0.d0,-1.d0,0.d0/)
    !
    var1=0.5d0*lz/chord
    var2=      lx/chord
    nface=nface+1
    tempface(nface)%a=x1
    tempface(nface)%b=x3
    tempface(nface)%c=x4
    tempface(nface)%normdir=(/var1,0.d0,-var2/)
    !
    var1=0.5d0*lz/chord
    var2=      lx/chord
    nface=nface+1
    tempface(nface)%a=x2
    tempface(nface)%b=x3
    tempface(nface)%c=x4
    tempface(nface)%normdir=(/var1,0.d0, var2/)
    !
    var1=h /sqrt(lx**2+h**2)
    var2=lx/sqrt(lx**2+h**2)
    nface=nface+1
    tempface(nface)%a=x1
    tempface(nface)%b=x2
    tempface(nface)%c=x4
    tempface(nface)%normdir=(/-var1,var2,0.d0/)
    !
    immbody(2)%num_face=nface
    call immbody(2)%alloface()
    immbody(2)%face(1:nface)=tempface(1:nface)
    !
    call tecsolid('tecsolid.plt',immbody)
    !
    print*,' ** solid generated'
    !
  end subroutine solidgen_mvg
  !+-------------------------------------------------------------------+
  !| This subroutine is to generate a solid usign stl file format.     |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 05-07-2021  | Created by J. Fang @ Warrington                     |
  !+-------------------------------------------------------------------+
  subroutine solidgen_triagnle
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
    nface=0
    !
    x1(1)=0.d0
    x1(2)=0.d0
    x1(3)=0.d0
    !
    x2(1)=1.d0
    x2(2)=-0.5d0
    x2(3)=0.d0
    !
    x3(1)=1.d0
    x3(2)=0.5d0
    x3(3)=0.d0
    !
    nface=nface+1
    tempface(nface)%a=x1
    tempface(nface)%b=x2
    tempface(nface)%c=x3
    tempface(nface)%normdir=(/0.d0,0.d0,-1.d0/)
    !
    x1(1)=0.d0
    x1(2)=0.d0
    x1(3)=1.d0
    !
    x2(1)=1.d0
    x2(2)=-0.5d0
    x2(3)=1.d0
    !
    x3(1)=1.d0
    x3(2)=0.5d0
    x3(3)=1.d0
    !
    nface=nface+1
    tempface(nface)%a=x1
    tempface(nface)%b=x2
    tempface(nface)%c=x3
    tempface(nface)%normdir=(/0.d0,0.d0, 1.d0/)
    !
    x1(1)=1.d0
    x1(2)=-0.5d0
    x1(3)=0.d0
    !
    x2(1)=1.d0
    x2(2)=-0.5d0
    x2(3)=1.d0
    !
    x3(1)=1.d0
    x3(2)=0.5d0
    x3(3)=1.d0
    !
    x4(1)=1.d0
    x4(2)=0.5d0
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
    x3(1)=1.d0
    x3(2)=0.5d0
    x3(3)=1.d0
    !
    x4(1)=1.d0
    x4(2)=0.5d0
    x4(3)=0.d0
    !
    nface=nface+1
    tempface(nface)%a=x1
    tempface(nface)%b=x2
    tempface(nface)%c=x3
    tempface(nface)%normdir=(/-1.d0/sqrt(5.d0),2.d0/sqrt(5.d0),0.d0/)
    !
    nface=nface+1
    tempface(nface)%a=x3
    tempface(nface)%b=x4
    tempface(nface)%c=x1
    tempface(nface)%normdir=(/-1.d0/sqrt(5.d0),2.d0/sqrt(5.d0),0.d0/)
    !
    x1(1)=0.d0
    x1(2)=0.d0
    x1(3)=0.d0
    !
    x2(1)=0.d0
    x2(2)=0.d0
    x2(3)=1.d0
    !
    x3(1)=1.d0
    x3(2)=-0.5d0
    x3(3)=1.d0
    !
    x4(1)=1.d0
    x4(2)=-0.5d0
    x4(3)=0.d0
    !
    !
    nface=nface+1
    tempface(nface)%a=x1
    tempface(nface)%b=x2
    tempface(nface)%c=x3
    tempface(nface)%normdir=(/-1.d0/sqrt(5.d0),-2.d0/sqrt(5.d0),0.d0/)
    !
    nface=nface+1
    tempface(nface)%a=x3
    tempface(nface)%b=x4
    tempface(nface)%c=x1
    tempface(nface)%normdir=(/-1.d0/sqrt(5.d0),-2.d0/sqrt(5.d0),0.d0/)
    !
    !
    immbody(1)%num_face=nface
    call immbody(1)%alloface()
    immbody(1)%face(1:nface)=tempface(1:nface)
    !
    call tecsolid('tecsolid.plt',immbody)
    !
    print*,' ** solid generated'
    !
  end subroutine solidgen_triagnle
  !+-------------------------------------------------------------------+
  !| The end of the subroutine solidgen_triagnle.                      |
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
    nface=0
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
  !| This subroutine is to generate a sphere STL file with many        |
  !| equilateral triangle.                                             |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 05-07-2021  | Created by J. Fang @ Warrington                     |
  !+-------------------------------------------------------------------+
  subroutine solidgen_sphere_tri
    !
    use commtype, only : solid,triangle
    use readwrite,only : readsolid
    use commvar,  only : nsolid,immbody
    use geom,     only : solidrange,solidresc,solidshif,solidimpro
    use tecio,     only : tecsolid
    use stlaio,    only : stla_write
    use commfunc,  only : cross_product,dis2point
    !
    ! local data
    integer :: i,j,k,jf,km,jm,im,nface,nface2,jface,ntrimax,m,n,l
    real(8) :: dthe,dphi,theter1,theter2,phi1,phi2,radi,   &
               x1(3),x2(3),x3(3),x4(3),norm1(3),norm2(3),var1,var2, &
               ab(3),ac(3),bc(3),abc(3),xico(12,3),disab,disac,disbc
    integer :: lmnsav(120,3)
    real(8) :: epsilon
    type(triangle),allocatable :: tempface(:),tempface2(:)
    !
    epsilon=1.d-12
    ntrimax=10000
    !
    print*,' ** generating solid'
    !
    nsolid=1
    allocate(immbody(nsolid))
    immbody(1)%name='sphere'
    !
    allocate(tempface(ntrimax),tempface2(2*ntrimax))
    !
    ! first build a icosahedron
    !
    var1=0.5d0*(sqrt(5.d0)+1.d0)
    !
    xico(1,1)=0.d0; xico(2,1)= 0.d0; xico(3,1)= 0.d0; xico(4,1)= 0.d0
    xico(1,2)=1.d0; xico(2,2)= 1.d0; xico(3,2)=-1.d0; xico(4,2)=-1.d0
    xico(1,3)=var1; xico(2,3)=-var1; xico(3,3)=-var1; xico(4,3)= var1
    !
    xico(5,1)=1.d0; xico(6,1)= 1.d0; xico(7,1)=-1.d0; xico(8,1)=-1.d0
    xico(5,2)=var1; xico(6,2)=-var1; xico(7,2)=-var1; xico(8,2)= var1
    xico(5,3)=0.d0; xico(6,3)= 0.d0; xico(7,3)= 0.d0; xico(8,3)= 0.d0
    !
    xico(9,1)=var1; xico(10,1)=-var1; xico(11,1)=-var1; xico(12,1)= var1
    xico(9,2)=0.d0; xico(10,2)= 0.d0; xico(11,2)= 0.d0; xico(12,2)= 0.d0
    xico(9,3)=1.d0; xico(10,3)= 1.d0; xico(11,3)=-1.d0; xico(12,3)=-1.d0
    !
    radi=sqrt(xico(1,1)**2+xico(1,2)**2+xico(1,3)**2)
    print*,' ** radi= ',radi
    xico=xico/radi
    var2=2.d0/radi
    !
    nface=0
    lmnsav=0.d0
    do l=1,12
    do m=1,12
    do n=1,12
      !
      disab=abs(dis2point(xico(m,:),xico(n,:))-var2)
      disac=abs(dis2point(xico(l,:),xico(n,:))-var2)
      disbc=abs(dis2point(xico(m,:),xico(l,:))-var2)
      !
      if(disab<epsilon .and. disac<epsilon .and. disbc<epsilon) then
        !
        do i=1,nface
          !
          if( (lmnsav(i,1)==l .and. lmnsav(i,2)==m .and. lmnsav(i,3)==n) .or. &
              (lmnsav(i,1)==l .and. lmnsav(i,2)==n .and. lmnsav(i,3)==m) .or. &
              (lmnsav(i,1)==m .and. lmnsav(i,2)==l .and. lmnsav(i,3)==n) .or. &
              (lmnsav(i,1)==m .and. lmnsav(i,2)==n .and. lmnsav(i,3)==l) .or. &
              (lmnsav(i,1)==n .and. lmnsav(i,2)==l .and. lmnsav(i,3)==m) .or. &
              (lmnsav(i,1)==n .and. lmnsav(i,2)==m .and. lmnsav(i,3)==l) ) then
            exit
          endif
        enddo
        !
        if(i==nface+1) then
          !
          nface=nface+1
          !
          print*,' ** nface= ',nface
          !
          lmnsav(nface,1)=l
          lmnsav(nface,2)=m
          lmnsav(nface,3)=n
          !
          tempface(nface)%a=xico(l,:)
          tempface(nface)%b=xico(m,:)
          tempface(nface)%c=xico(n,:)
          !
          norm1=cross_product(tempface(nface)%b-tempface(nface)%a,         &
                              tempface(nface)%c-tempface(nface)%a )
          norm2=tempface(nface)%a
          !
          if(dot_product(norm1,norm2)<0.d0) then
            !
            tempface(nface)%a=xico(l,:)
            tempface(nface)%b=xico(n,:)
            tempface(nface)%c=xico(m,:)
            !
            norm1=cross_product(tempface(nface)%b-tempface(nface)%a,         &
                                tempface(nface)%c-tempface(nface)%a )
            norm2=tempface(nface)%a
            !
          endif
          !
          !
        endif
        !
      endif
      !
    enddo
    enddo
    enddo
    !
    do while(nface<ntrimax)
      !
      nface2=0
      !
      do jface=1,nface
        !
        ab=0.5d0*(tempface(jface)%a+tempface(jface)%b)
        ac=0.5d0*(tempface(jface)%a+tempface(jface)%c)
        bc=0.5d0*(tempface(jface)%b+tempface(jface)%c)
        abc=num1d3*(tempface(jface)%a+tempface(jface)%b+tempface(jface)%c)
        !
        var1=sqrt(ab(1)**2+ab(2)**2+ab(3)**2)
        ab=ab/var1
        var1=sqrt(ac(1)**2+ac(2)**2+ac(3)**2)
        ac=ac/var1
        var1=sqrt(bc(1)**2+bc(2)**2+bc(3)**2)
        bc=bc/var1
        var1=sqrt(abc(1)**2+abc(2)**2+abc(3)**2)
        abc=abc/var1
        !
        nface2=nface2+1
        tempface2(nface2)%a=tempface(jface)%a
        tempface2(nface2)%b=ab
        tempface2(nface2)%c=ac
        !
        nface2=nface2+1
        tempface2(nface2)%a=tempface(jface)%b
        tempface2(nface2)%b=ab
        tempface2(nface2)%c=bc
        !
        nface2=nface2+1
        tempface2(nface2)%a=tempface(jface)%c
        tempface2(nface2)%b=ac
        tempface2(nface2)%c=bc
        !
        nface2=nface2+1
        tempface2(nface2)%a=ab
        tempface2(nface2)%b=ac
        tempface2(nface2)%c=bc
        !
        if(nface2>ntrimax) exit
        !
      enddo
      !
      if(nface2>ntrimax) exit
      !
      nface=nface2
      tempface=tempface2
      !
      print*,' ** nface=',nface
      !
    enddo
    !
    !
    immbody(1)%num_face=nface
    call immbody(1)%alloface()
    immbody(1)%face(1:nface)=tempface(1:nface)
    !
    do jf=1,immbody(1)%num_face
      !
      norm1=cross_product(immbody(1)%face(jf)%a-immbody(1)%face(jf)%b,         &
                          immbody(1)%face(jf)%c-immbody(1)%face(jf)%b )
      norm2=immbody(1)%face(jf)%a
      !
      if(dot_product(norm1,norm2)<0.d0) then
        x1=immbody(1)%face(jf)%c
        immbody(1)%face(jf)%c=immbody(1)%face(jf)%a
        immbody(1)%face(jf)%a=x1
        !
        norm1=cross_product(immbody(1)%face(jf)%a-immbody(1)%face(jf)%b,       &
                            immbody(1)%face(jf)%c-immbody(1)%face(jf)%b )
      endif
      !
      var1=sqrt(norm1(1)**2+norm1(2)**2+norm1(3)**2)
      !
      if(abs(var1)<epsilon) then
        !
        norm1=cross_product(immbody(1)%face(jf)%C-immbody(1)%face(jf)%a,         &
                            immbody(1)%face(jf)%b-immbody(1)%face(jf)%a )
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
  end subroutine solidgen_sphere_tri
  !+-------------------------------------------------------------------+
  !| The end of the subroutine solidgen_sphere_tri.                    |
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
  !+-------------------------------------------------------------------+
  !| This function is to generate a  4-digit NACA airfoil.             |
  !+-------------------------------------------------------------------+
  !| ref: https://en.wikipedia.org/wiki/NACA_airfoil
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 28-Jul-2021: Created by J. Fang @ Appleton                        |
  !+-------------------------------------------------------------------+
  subroutine naca4digit(xin,x,y,name,surface)
    !
    ! arguments
    real(8),intent(in) :: xin
    real(8),intent(out) :: x,y
    character(len=4),intent(in) :: name
    character(len=5),intent(in) :: surface
    !
    ! local data 
    real(8) :: yt,t,yc,theter,mr,pr
    integer :: m,p
    !
    read(name(1:1),*)m
    read(name(2:2),*)p
    read(name(3:4),*)t
    !
    t=t/100.d0
    !
    yt=5.d0*t*(0.2969d0*sqrt(xin)-0.1260d0*xin-0.3516d0*xin*xin+       &
               0.2843d0*xin**3-0.1036d0*xin**4)
    !
    if(p==0 .and. m==0) then
      ! NACA-00XX, a symmetrical 4-digit NACA airfoil
      !
      if(surface=='upper') then
        x=xin
        y=yt
      elseif(surface=='lower') then
        x=xin
        y=-yt
      else
        print*,' surface=',surface
        stop '!! ERROR 1 @ naca4digit'
      endif
      !
    else
      !
      mr=0.01d0*dble(m)
      pr=0.1d0*dble(p)
      !
      if(xin>=0.d0 .and. xin<=pr) then
        yc=mr/pr/pr*(2.d0*pr*xin-xin*xin)
        theter=atan(2.d0*mr/pr/pr*(pr-xin))
      elseif(xin>=pr .and. xin<=1.d0) then
        yc=mr/(1-pr)**2*(1.d0-2.d0*pr+2.d0*pr*xin-xin**2)
        theter=atan(2.d0*mr/(1-pr)**2*(pr-xin))
      else
        stop '!! ERROR 2 @ naca4digit'
      endif
      !
      if(surface=='upper') then
        x=xin-yt*sin(theter)
        y= yc+yt*cos(theter)
      elseif(surface=='lower') then
        x=xin+yt*sin(theter)
        y= yc-yt*cos(theter)
      else
        print*,' surface=',surface
        stop '!! ERROR 3 @ naca4digit'
      endif
      !
    endif
    !
  end subroutine naca4digit
  !+-------------------------------------------------------------------+
  !| The end of the subroutine naca4digit.                             |
  !+-------------------------------------------------------------------+
  !
end module pp
!+---------------------------------------------------------------------+
!| The end of the module pp                                            |
!+---------------------------------------------------------------------+