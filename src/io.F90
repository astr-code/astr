!+---------------------------------------------------------------------+
!| This module contains subroutines for io purpose.                    |
!| ==============                                                      |
!| CHANGE RECORD                                                       |
!| -------------                                                       |
!| 18-June-2023  | Created by J. Fang @ Warrington                     |
!+---------------------------------------------------------------------+
module io
  !
  use commvar, only: ndims,im,jm,km,iomode,lreport,ltimrpt
  use parallel,only: mpirank,mpirankname,ptime
  use stlaio,  only: get_unit
  use utility, only: timereporter
  !
  implicit none
  !
  contains
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to read a checkpoint file for restarting. |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 24-06-2021  | Created by J. Fang @ Warrington                     |
  !+-------------------------------------------------------------------+
  subroutine readcheckpoint(folder,mode)
    !
    use commvar, only: nstep,filenumb,fnumslic,time,flowtype,           &
                       num_species,im,jm,km,force,numq,turbmode,lwsequ, &
                       nondimen,moment,iomode
    use commarray, only : rho,vel,prs,tmp,spc,q,tke,omg,miut,sigma,qflux
    use methodmoment, only: mijk, Rij
    use statistic, only : massflux,massflux_target,nsamples
    use hdf5io
    use bc,        only : ninflowslice
#ifdef COMB
    use thermchem, only: spcindex
#endif
    !
    ! arguments
    character(len=*),intent(in) :: folder
    character(len=1),intent(in),optional :: mode
    !
    ! local data
    integer :: nstep_1,jsp
    character(len=1) :: modeio
    character(len=2) :: qname
    character(len=3) :: spname
    character(len=128) :: infilename
    character(len=4) :: stepname
    logical :: lexist
    !
    if(present(mode)) then
      modeio=mode
    else
      modeio=iomode
    endif
    !
    call h5io_init(filename=folder//'/auxiliary.h5',mode='read')
    call h5read(varname='nstep',var=nstep)
    call h5read(varname='filenumb',var=filenumb)
    call h5read(varname='fnumslic',var=fnumslic)
    call h5read(varname='ninflowslice',var=ninflowslice)
    if(flowtype=='channel') then
      call h5read(varname='massflux',var=massflux)
      call h5read(varname='massflux_target',var=massflux_target)
      call h5read(varname='force',var=force,dim=3)
    endif
    call h5read(varname='nsamples',var=nsamples)
    call h5io_end
    !
    write(stepname,'(i4.4)')filenumb
    infilename='outdat/flowfield'//stepname//'.'//modeio//'5'
    !
    ! if the file is not found, just go to the default flow field file
    inquire(file=infilename, exist=lexist)
    if(.not. lexist) then
      infilename=folder//'/flowfield.'//modeio//'5'
    endif
    !
    call h5io_init(filename=trim(infilename),mode='read')
    call h5read(varname='nstep',var=nstep_1)
    !
    if(nstep_1==nstep) then
      call h5read(varname='time',var=time)
      !
      call h5read(varname='ro',  var=rho(0:im,0:jm,0:km)  ,mode=modeio)
      call h5read(varname='u1',  var=vel(0:im,0:jm,0:km,1),mode=modeio)
      call h5read(varname='u2',  var=vel(0:im,0:jm,0:km,2),mode=modeio)
      call h5read(varname='u3',  var=vel(0:im,0:jm,0:km,3),mode=modeio)
      call h5read(varname='p',   var=prs(0:im,0:jm,0:km)  ,mode=modeio)
      call h5read(varname='t',   var=tmp(0:im,0:jm,0:km)  ,mode=modeio)
      do jsp=1,num_species
         write(spname,'(i3.3)')jsp
        call h5read(varname='sp'//spname,var=spc(0:im,0:jm,0:km,jsp),mode=modeio)
      enddo
      !
      if(trim(turbmode)=='k-omega') then
        call h5read(varname='k',     var=tke(0:im,0:jm,0:km),mode=modeio)
        call h5read(varname='omega', var=omg(0:im,0:jm,0:km),mode=modeio)
        call h5read(varname='miut',  var=miut(0:im,0:jm,0:km),mode=modeio)
      elseif(trim(turbmode)=='udf1') then
        call h5read(varname='miut',  var=miut(0:im,0:jm,0:km),mode=modeio)
      endif
      !
      if(moment=='r13' .or. moment == 'r26') then
        !
        call h5read(varname='sigmaxx',var=sigma(0:im,0:jm,0:km,1),mode=modeio)
        call h5read(varname='sigmaxy',var=sigma(0:im,0:jm,0:km,2),mode=modeio)
        call h5read(varname='sigmaxz',var=sigma(0:im,0:jm,0:km,3),mode=modeio)
        call h5read(varname='sigmayy',var=sigma(0:im,0:jm,0:km,4),mode=modeio)
        call h5read(varname='sigmayz',var=sigma(0:im,0:jm,0:km,5),mode=modeio)
        call h5read(varname='sigmazz',var=sigma(0:im,0:jm,0:km,6),mode=modeio)
          !
        call h5read(varname='qx',var=qflux(0:im,0:jm,0:km,1),mode=modeio)
        call h5read(varname='qy',var=qflux(0:im,0:jm,0:km,2),mode=modeio)
        call h5read(varname='qz',var=qflux(0:im,0:jm,0:km,3),mode=modeio)
        !
        call h5read(varname='mxxx',var=mijk(0:im,0:jm,0:km,1),mode=modeio)
        call h5read(varname='mxxy',var=mijk(0:im,0:jm,0:km,2),mode=modeio)
        call h5read(varname='mxxz',var=mijk(0:im,0:jm,0:km,3),mode=modeio)
        call h5read(varname='mxyy',var=mijk(0:im,0:jm,0:km,4),mode=modeio)
        call h5read(varname='myyy',var=mijk(0:im,0:jm,0:km,5),mode=modeio)
        call h5read(varname='myyz',var=mijk(0:im,0:jm,0:km,6),mode=modeio)
        call h5read(varname='mxyz',var=mijk(0:im,0:jm,0:km,7),mode=modeio)

        call h5read(varname='Rxx',var=Rij(0:im,0:jm,0:km,1),mode=modeio)
        call h5read(varname='Rxy',var=Rij(0:im,0:jm,0:km,2),mode=modeio)
        call h5read(varname='Rxz',var=Rij(0:im,0:jm,0:km,3),mode=modeio)
        call h5read(varname='Ryy',var=Rij(0:im,0:jm,0:km,4),mode=modeio)
        call h5read(varname='Ryz',var=Rij(0:im,0:jm,0:km,5),mode=modeio)
        !
      endif
      ! call h5io_init('checkpoint/qdata.'//modeio//'5',mode='read')
      ! do jsp=1,numq
      !   write(qname,'(i2.2)')jsp
      !   call h5read(varname='q'//qname,var=q(0:im,0:jm,0:km,jsp))
      ! enddo
      !
      call h5io_end
      !
    else
      if(lio)  print*,' !! flowfield.'//modeio//'5 NOT consistent with auxiliary.'//modeio//'5'
      if(lio)  print*,' nstep =',nstep,' in auxiliary.h5 '
      if(lio)  print*,' nstep =',nstep_1,' in flowfield.'//modeio//'5 '
      call mpistop
    endif
    !
    ! infilename='checkpoint/flowfield'//mpirankname
    ! open(16,file=trim(infilename),form='unformatted')
    ! read(16)nstep,time
    ! read(16)rho,vel,prs,tmp
    ! read(16)spc
    ! close(16)
    ! if(lio) print*,' >> ',trim(infilename)
    ! !
    ! infilename='checkpoint/qdata'//mpirankname
    ! open(16,file=trim(infilename),form='unformatted')
    ! read(16)nstep
    ! read(16)q
    ! close(16)
    ! if(lio) print*,' >> ',trim(infilename)
    ! !
    ! infilename='checkpoint/auxiliary'//mpirankname
    ! open(16,file=trim(infilename),form='unformatted')
    ! read(16)nstep,filenumb
    ! read(16)massflux,massflux_target
    ! read(16)force
    ! close(16)
    ! if(lio) print*,' >> ',trim(infilename)
    !
    ! use for the initialization for the reacting supersonic backstep case
!     if(.not.nondimen) then
#ifdef COMB
!     ! phi = 0.4 
!     spc(0:im,0:jm,0:km,:)=0.d0
!     spc(0:im,0:jm,0:km,spcindex('O2'))=0.2302d0
!     spc(0:im,0:jm,0:km,spcindex('H2'))=0.0116d0
!     spc(0:im,0:jm,0:km,spcindex('N2'))=0.7582d0
!     print*,' ** phi=0.4 '
      ! phi=0.2
      ! spc(0:im,0:jm,0:km,:)=0.d0
      ! spc(0:im,0:jm,0:km,spcindex('O2'))=0.23154d0
      ! spc(0:im,0:jm,0:km,spcindex('H2'))=0.00583d0
      ! spc(0:im,0:jm,0:km,spcindex('N2'))=0.76263d0
      ! if(lio) print*,' ** phi=0.2 '
      ! phi=0.3
      ! spc(0:im,0:jm,0:km,:)=0.d0
      ! spc(0:im,0:jm,0:km,spcindex('O2'))=0.23087d0
      ! spc(0:im,0:jm,0:km,spcindex('H2'))=0.00873d0
      ! spc(0:im,0:jm,0:km,spcindex('N2'))=0.76040d0
      ! if(lio) print*,' ** phi=0.3 '
#endif
!     endif
    !
    if(lio)  print*,' ** checkpoint file read. '
    !
  end subroutine readcheckpoint
  !+-------------------------------------------------------------------+
  !| The end of the subroutine readcheckpoint.                         |
  !+-------------------------------------------------------------------+
  !!
  !+-------------------------------------------------------------------+
  !| This subroutine is used to write sequence of flow field.          |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 28-12-2021  | Created by J. Fang @ Warrington                     |
  !+-------------------------------------------------------------------+
  subroutine writeflfed(timerept)
    !
    use commvar,  only: time,nstep,filenumb,fnumslic,num_species,im,jm, &
                       km,lwsequ,turbmode,feqwsequ,force,ymin,ymax,moment,iomode
    use readwrite,only: nxtwsequ
    use commarray,only : x,rho,vel,prs,tmp,spc,q,ssf,lshock,crinod,sigma,qflux
    use methodmoment,only: mijk,Rij
    use models,   only : tke,omg,miut
    use statistic,only : nsamples,liosta,massflux,massflux_target
    use bc,       only : ninflowslice
    use hdf5io
#ifdef COMB
    use thermchem,only : heatrate
#endif
    !
    ! arguments
    logical,intent(in),optional :: timerept
    !
    ! local data
    integer :: i,j,k,jsp
    character(len=4) :: stepname
    character(len=64) :: outfilename,outauxiname
    character(len=64),save :: savfilenmae='first'
    character(len=3) :: spname
    real(8),allocatable :: rshock(:,:,:),rcrinod(:,:,:),hrr(:)
    !
    real(8) :: ypos
    logical :: lwprofile
    real(8) :: time_beg
    real(8),save :: subtime=0.d0
    !
    if(present(timerept) .and. timerept) time_beg=ptime()
    !
    if(lwsequ .and. nstep==nxtwsequ) then
      !
      if(nstep/=0) filenumb=filenumb+1
      !
      write(stepname,'(i4.4)')filenumb
      !
      outfilename='outdat/flowfield'//stepname//'.'//iomode//'5'
      outauxiname='outdat/auxiliary'//stepname//'.'//iomode//'5'
      !
    else
      stepname=''
      outfilename='outdat/flowfield.'//iomode//'5'
      outauxiname='outdat/auxiliary.'//iomode//'5'
    endif
    !
    ! outfilename='outdat/flowfield.h5'
    call h5io_init(trim(outfilename),mode='write')
    !
    call h5write(varname='ro',var=rho(0:im,0:jm,0:km),  mode=iomode)
    !
    call h5write(varname='u1',var=vel(0:im,0:jm,0:km,1),mode=iomode)
    call h5write(varname='u2',var=vel(0:im,0:jm,0:km,2),mode=iomode)
    call h5write(varname='u3',var=vel(0:im,0:jm,0:km,3),mode=iomode)
    call h5write(varname='p', var=prs(0:im,0:jm,0:km),  mode=iomode)
    call h5write(varname='t', var=tmp(0:im,0:jm,0:km),  mode=iomode)
    if(num_species>0) then
      do jsp=1,num_species
         write(spname,'(i3.3)')jsp
        call h5write(varname='sp'//spname,var=spc(0:im,0:jm,0:km,jsp),mode=iomode)
      enddo
    endif
    !
    if(allocated(ssf)) then
      call h5write(varname='ssf', var=ssf(0:im,0:jm,0:km),mode=iomode)
    endif
    if(allocated(lshock)) then
      allocate(rshock(0:im,0:jm,0:km))
      allocate(rcrinod(0:im,0:jm,0:km))
      do k=0,km
      do j=0,jm
      do i=0,im
        !
        if(lshock(i,j,k)) then
          rshock(i,j,k)=1.d0
        else
          rshock(i,j,k)=0.d0
        endif
        !
        if(crinod(i,j,k)) then
          rcrinod(i,j,k)=1.d0
        else
          rcrinod(i,j,k)=0.d0
        endif
        !
      enddo
      enddo
      enddo
      !
      ! call h5write(varname='lshk', var=rshock(0:im,0:jm,0:km),mode=iomode)
      ! call h5write(varname='crit', var=rcrinod(0:im,0:jm,0:km),mode=iomode)
    endif
    !
    if(trim(turbmode)=='k-omega') then
      call h5write(varname='k',     var=tke(0:im,0:jm,0:km),mode=iomode)
      call h5write(varname='omega', var=omg(0:im,0:jm,0:km),mode=iomode)
      call h5write(varname='miut',  var=miut(0:im,0:jm,0:km),mode=iomode)
    elseif(trim(turbmode)=='udf1') then
      call h5write(varname='miut',  var=miut(0:im,0:jm,0:km),mode=iomode)
    endif
    !
    if(moment=='r13' .or. moment == 'r26') then
      !
      call h5write(varname='sigmaxx',var=sigma(0:im,0:jm,0:km,1),mode=iomode)
      call h5write(varname='sigmaxy',var=sigma(0:im,0:jm,0:km,2),mode=iomode)
      call h5write(varname='sigmaxz',var=sigma(0:im,0:jm,0:km,3),mode=iomode)
      call h5write(varname='sigmayy',var=sigma(0:im,0:jm,0:km,4),mode=iomode)
      call h5write(varname='sigmayz',var=sigma(0:im,0:jm,0:km,5),mode=iomode)
      call h5write(varname='sigmazz',var=sigma(0:im,0:jm,0:km,6),mode=iomode)
      !
      call h5write(varname='qx',var=qflux(0:im,0:jm,0:km,1),mode=iomode)
      call h5write(varname='qy',var=qflux(0:im,0:jm,0:km,2),mode=iomode)
      call h5write(varname='qz',var=qflux(0:im,0:jm,0:km,3),mode=iomode)
      !
      call h5write(varname='mxxx',var=mijk(0:im,0:jm,0:km,1),mode=iomode)
      call h5write(varname='mxxy',var=mijk(0:im,0:jm,0:km,2),mode=iomode)
      call h5write(varname='mxxz',var=mijk(0:im,0:jm,0:km,3),mode=iomode)
      call h5write(varname='mxyy',var=mijk(0:im,0:jm,0:km,4),mode=iomode)
      call h5write(varname='myyy',var=mijk(0:im,0:jm,0:km,5),mode=iomode)
      call h5write(varname='myyz',var=mijk(0:im,0:jm,0:km,6),mode=iomode)
      call h5write(varname='mxyz',var=mijk(0:im,0:jm,0:km,7),mode=iomode)

      call h5write(varname='Rxx',var=Rij(0:im,0:jm,0:km,1),mode=iomode)
      call h5write(varname='Rxy',var=Rij(0:im,0:jm,0:km,2),mode=iomode)
      call h5write(varname='Rxz',var=Rij(0:im,0:jm,0:km,3),mode=iomode)
      call h5write(varname='Ryy',var=Rij(0:im,0:jm,0:km,4),mode=iomode)
      call h5write(varname='Ryz',var=Rij(0:im,0:jm,0:km,5),mode=iomode)
      !
    endif
    !
    call h5io_end
    !
    if(lio) then
      !
      call h5srite(varname='nstep',var=nstep,filename=trim(outfilename))
      call h5srite(varname='time',var=time,filename=trim(outfilename))
      !

      call h5srite(varname='nstep',var=nstep,                          &
                          filename=trim(outauxiname),newfile=.true.)
      call h5srite(varname='filenumb',var=filenumb,                    &
                      filename=trim(outauxiname))
      call h5srite(varname='fnumslic',var=fnumslic,                    &
                                         filename=trim(outauxiname))
      call h5srite(varname='ninflowslice',var=ninflowslice,            &
                                         filename=trim(outauxiname))
      call h5srite(varname='massflux',var=massflux,                    &
                                         filename=trim(outauxiname))
      call h5srite(varname='massflux_target',var=massflux_target,      &
                                         filename=trim(outauxiname))
      call h5srite(varname='force',var=force,                          &
                                         filename=trim(outauxiname))
      call h5srite(varname='nsamples',var=nsamples,                    &
                                         filename=trim(outauxiname))
    endif
    !
    ! if(trim(savfilenmae)=='first' .or. savfilenmae==outfilename) then
    !   call xdmfwriter(flowh5file=trim(outfilename),dataname='ro',timesec=time,mode='newvisu')
    ! else
    !   call xdmfwriter(flowh5file=trim(outfilename),dataname='ro',timesec=time,mode='animati')
    ! endif
    ! !
    ! call xdmfwriter(flowh5file=trim(outfilename),dataname='u1',mode='data')
    ! call xdmfwriter(flowh5file=trim(outfilename),dataname='u2',mode='data')
    ! call xdmfwriter(flowh5file=trim(outfilename),dataname='u3',mode='data')
    ! call xdmfwriter(flowh5file=trim(outfilename),dataname='p',mode='data')
    ! call xdmfwriter(flowh5file=trim(outfilename),dataname='t',mode='data')
    !
    if(ndims==1) then
      !
      ! open(18,file='outdat/profile'//trim(stepname)//mpirankname//'.dat')
      ! write(18,"(5(1X,A15))")'x','ro','u','p','t'
      ! write(18,"(5(1X,E15.7E3))")(x(i,0,0,1),rho(i,0,0),vel(i,0,0,1),  &
      !                             prs(i,0,0),tmp(i,0,0),i=0,im)
      ! close(18)
      ! print*,' << outdat/profile',trim(stepname),mpirankname,'.dat'
      !
    elseif(ndims==2) then
      !
      ! call tecbin('testout/tecfield'//mpirankname//'.plt',           &
      !                                   x(0:im,0:jm,0,1),'x',        &
      !                                   x(0:im,0:jm,0,2),'y',        &
      !                                   rho(0:im,0:jm,0),'ro',       &
      !                                 vel(0:im,0:jm,0,1),'u',        &
      !                                 vel(0:im,0:jm,0,2),'v',        &
      !                                   tmp(0:im,0:jm,0),'T',        &
      !                                   prs(0:im,0:jm,0),'p',        &
      !                               rcrinod(0:im,0:jm,0),'c' )
      !
    endif
    !
    !
    savfilenmae=outfilename
    !
    nxtwsequ=nstep+feqwsequ
    !
    if(present(timerept) .and. timerept) then
      !
      subtime=subtime+ptime()-time_beg
      !
      if(lio .and. lreport) call timereporter(routine='writeflfed', &
                                              timecost=subtime, &
                                              message='write flow data')
      !
    endif
    !
  end subroutine writeflfed
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to dump all computation data              |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 12-02-2021  | Created by J. Fang @ Warrington                     |
  !+-------------------------------------------------------------------+
  subroutine writechkpt(stepsequ,timerept)
    !
    use commvar, only: time,nstep,filenumb,fnumslic,num_species,im,jm, &
                       km,lwsequ,lavg,force,numq,imbroot,limmbou,      &
                       turbmode,feqchkpt,iomode
    use commarray,only : x,rho,vel,prs,tmp,spc,q,ssf,lshock,crinod
    use models,   only : tke,omg,miut
    use statistic,only : liosta!,nsamples,massflux,massflux_target
    ! use bc,       only : ninflowslice
    !
    use readwrite,only : nxtchkpt,bakupfile
    use hdf5io
    !
    !
    ! arguments
    integer,intent(in) :: stepsequ
    logical,intent(in),optional :: timerept
    !
    ! local data
    real(8),allocatable :: rshock(:,:,:),rcrinod(:,:,:)
    character(len=2) :: qname
    character(len=3) :: spname
    integer :: jsp,i,fh,ios,j,k
    logical :: lfopen
    real(8) :: time_beg
    real(8),save :: subtime=0.d0
    !
    if(present(timerept) .and. timerept) time_beg=ptime() 
    !
    ! if(lwsequ) then
    !   write(stepname,'(i4.4)')filenumb
    !   filenumb=filenumb+1
    !   !
    !   outfilename='outdat/flowfield'//stepname//'.h5'
    ! else
    !   stepname=''
    !   outfilename='outdat/flowfield.h5'
    ! endif
    !
    ! outfilename='outdat/flowfield'//mpirankname
    ! open(21,file=trim(outfilename),form='unformatted')
    ! write(21)nstep,time
    ! write(21)rho,vel,prs,tmp
    ! write(21)spc
    ! close(21)
    ! if(lio) print*,' << ',trim(outfilename)
    ! !
    ! outfilename='outdat/qdata'//mpirankname
    ! open(21,file=trim(outfilename),form='unformatted')
    ! write(21)nstep
    ! write(21)q
    ! close(21)
    ! if(lio) print*,' << ',trim(outfilename)
    ! !
    ! outfilename='outdat/auxiliary'//mpirankname
    ! open(21,file=trim(outfilename),form='unformatted')
    ! write(21)nstep,filenumb
    ! write(21)massflux,massflux_target
    ! write(21)force
    ! close(21)
    ! if(lio) print*,' << ',trim(outfilename)
    !
    if(lio.and.(.not.lwsequ)) call bakupfile('outdat/auxiliary.h5')
    if(lio.and.(.not.lwsequ)) call bakupfile('outdat/flowfield.'//iomode//'5')
    !
    call writeflfed()
    ! call writeflfed_2d()
    !
    ! call h5io_init('outdat/qdata.h5',mode='write')
    ! do jsp=1,numq
    !    write(qname,'(i2.2)')jsp
    !   call h5write(varname='q'//qname,var=q(0:im,0:jm,0:km,jsp))
    ! enddo
    ! call h5io_end
    !
    ! if(lio) then
    !   !
    !   call h5srite(varname='nstep',var=nstep,                          &
    !                       filename='outdat/auxiliary.h5',newfile=.true.)
    !   call h5srite(varname='filenumb',var=filenumb,                    &
    !                                      filename='outdat/auxiliary.h5')
    !   call h5srite(varname='fnumslic',var=fnumslic,                    &
    !                                      filename='outdat/auxiliary.h5')
    !   call h5srite(varname='ninflowslice',var=ninflowslice,            &
    !                                      filename='outdat/auxiliary.h5')
    !   call h5srite(varname='massflux',var=massflux,                    &
    !                                      filename='outdat/auxiliary.h5')
    !   call h5srite(varname='massflux_target',var=massflux_target,      &
    !                                      filename='outdat/auxiliary.h5')
    !   call h5srite(varname='force',var=force,                          &
    !                                      filename='outdat/auxiliary.h5')
    !   call h5srite(varname='nsamples',var=nsamples,                    &
    !                                      filename='outdat/auxiliary.h5')
    ! endif
    !
    if(limmbou .and. mpirank==imbroot) then
      ! call imboundarydata 
      !
      ! call writecylimmbou
      !
    endif
    !
    ! if(irk==0 .and. jrk==jrkm) then
    !     print*,'------------------------------------'
    !     print*,x(0,jm,0,1),prs(0,jm,0)
    !     print*,'------------------------------------'
    ! endif
    !
    if(liosta) then
      !
      call writemeanflow
      !
    endif
    !
    ! flush all openned files
    fh=20
    lfopen=.true.
    !
    do while(lfopen)
      inquire(unit=fh, opened=lfopen,iostat=ios)
      !
      if(lfopen .and. ios==0) then
        flush(fh)
        fh=fh+1
      elseif(ios .ne. 0) then
        print*,' !! error with opening file unit:',fh
      endif
    enddo
    !
    nxtchkpt=nstep+feqchkpt
    !
    if(present(timerept) .and. timerept) then
      !
      subtime=subtime+ptime()-time_beg
      !
      if(lio .and. lreport) call timereporter(routine='writechkpt', &
                                              timecost=subtime, &
                                              message='write checkpoint')
      !
    endif
    !
    return
    !
  end subroutine writechkpt
  !+-------------------------------------------------------------------+
  !| The end of the subroutine writechkpt.                             |
  !+-------------------------------------------------------------------+
  !!
  !+-------------------------------------------------------------------+
  !| This subroutine is used to write mean flow statistics.            |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 28-12-2021  | Created by J. Fang @ Warrington                     |
  !+-------------------------------------------------------------------+
  subroutine writemeanflow
    !
    use commvar,  only : nstep,iomode
    use statistic,only : nsamples,nstep_sbeg,time_sbeg,                &
                         rom,u1m,u2m,u3m,pm,tm,                        &
                         u11,u22,u33,u12,u13,u23,pp,tt,tu1,tu2,tu3,    &
                         u111,u222,u333,u112,u113,                     &
                         u122,u133,u223,u233,u123,                     &
                         u1rem,u2rem,u3rem,pu1,pu2,pu3,                &
                         sgmam11,sgmam22,sgmam33,                      &
                         sgmam12,sgmam13,sgmam23,                      &
                         disspa,predil,visdif1,visdif2,visdif3
    use readwrite,only : bakupfile
    use hdf5io
    !
    if(lio) call bakupfile('outdat/meanflow.'//iomode//'5')
    !
    call h5io_init('outdat/meanflow.'//iomode//'5',mode='write')
    call h5write(varname='rom',var=rom(0:im,0:jm,0:km),mode=iomode)
    call h5write(varname='u1m',var=u1m(0:im,0:jm,0:km),mode=iomode)
    call h5write(varname='u2m',var=u2m(0:im,0:jm,0:km),mode=iomode)
    call h5write(varname='u3m',var=u3m(0:im,0:jm,0:km),mode=iomode)
    call h5write(varname='pm ',var=pm (0:im,0:jm,0:km),mode=iomode)
    call h5write(varname='tm ',var=tm (0:im,0:jm,0:km),mode=iomode)
    call h5io_end
    !
    if(lio) call bakupfile('outdat/2ndsta.'//iomode//'5')
    call h5io_init('outdat/2ndsta.'//iomode//'5',mode='write')
    call h5write(varname='u11',var=u11(0:im,0:jm,0:km),mode=iomode)
    call h5write(varname='u22',var=u22(0:im,0:jm,0:km),mode=iomode)
    call h5write(varname='u33',var=u33(0:im,0:jm,0:km),mode=iomode)
    call h5write(varname='u12',var=u12(0:im,0:jm,0:km),mode=iomode)
    call h5write(varname='u13',var=u13(0:im,0:jm,0:km),mode=iomode)
    call h5write(varname='u23',var=u23(0:im,0:jm,0:km),mode=iomode)
    !
    call h5write(varname='pp', var=pp(0:im,0:jm,0:km) ,mode=iomode)
    call h5write(varname='tt', var=tt(0:im,0:jm,0:km) ,mode=iomode)
    call h5write(varname='tu1',var=tu1(0:im,0:jm,0:km),mode=iomode)
    call h5write(varname='tu2',var=tu2(0:im,0:jm,0:km),mode=iomode)
    call h5write(varname='tu3',var=tu2(0:im,0:jm,0:km),mode=iomode)
    call h5io_end
    !
    if(lio) call bakupfile('outdat/3rdsta.'//iomode//'5')
    call h5io_init('outdat/3rdsta.'//iomode//'5',mode='write')
    call h5write(varname='u111',var=u111(0:im,0:jm,0:km),mode=iomode)
    call h5write(varname='u222',var=u222(0:im,0:jm,0:km),mode=iomode)
    call h5write(varname='u333',var=u333(0:im,0:jm,0:km),mode=iomode)
    call h5write(varname='u112',var=u112(0:im,0:jm,0:km),mode=iomode)
    call h5write(varname='u113',var=u113(0:im,0:jm,0:km),mode=iomode)
    call h5write(varname='u122',var=u122(0:im,0:jm,0:km),mode=iomode)
    call h5write(varname='u133',var=u133(0:im,0:jm,0:km),mode=iomode)
    call h5write(varname='u223',var=u223(0:im,0:jm,0:km),mode=iomode)
    call h5write(varname='u233',var=u233(0:im,0:jm,0:km),mode=iomode)
    call h5write(varname='u123',var=u123(0:im,0:jm,0:km),mode=iomode)
    call h5io_end
    !
    if(lio) call bakupfile('outdat/budget.'//iomode//'5')
    call h5io_init('outdat/budget.'//iomode//'5',mode='write')
    call h5write(varname='disspa', var=disspa(0:im,0:jm,0:km),mode=iomode)
    call h5write(varname='predil', var=predil(0:im,0:jm,0:km),mode=iomode)
    call h5write(varname='pu1',    var=pu1(0:im,0:jm,0:km)   ,mode=iomode)
    call h5write(varname='pu2',    var=pu2(0:im,0:jm,0:km)   ,mode=iomode)
    call h5write(varname='pu3',    var=pu3(0:im,0:jm,0:km)   ,mode=iomode)
    call h5write(varname='u1rem',  var=u1rem(0:im,0:jm,0:km) ,mode=iomode)
    call h5write(varname='u2rem',  var=u2rem(0:im,0:jm,0:km) ,mode=iomode)
    call h5write(varname='u3rem',  var=u3rem(0:im,0:jm,0:km) ,mode=iomode)
    call h5write(varname='visdif1',var=visdif1(0:im,0:jm,0:km),mode=iomode)
    call h5write(varname='visdif2',var=visdif2(0:im,0:jm,0:km),mode=iomode)
    call h5write(varname='visdif3',var=visdif3(0:im,0:jm,0:km),mode=iomode)
    call h5write(varname='sgmam11',var=sgmam11(0:im,0:jm,0:km),mode=iomode)
    call h5write(varname='sgmam22',var=sgmam22(0:im,0:jm,0:km),mode=iomode)
    call h5write(varname='sgmam33',var=sgmam33(0:im,0:jm,0:km),mode=iomode)
    call h5write(varname='sgmam12',var=sgmam12(0:im,0:jm,0:km),mode=iomode)
    call h5write(varname='sgmam13',var=sgmam13(0:im,0:jm,0:km),mode=iomode)
    call h5write(varname='sgmam23',var=sgmam23(0:im,0:jm,0:km),mode=iomode)
    call h5io_end
    !
    if(lio) then
      call h5srite(varname='nstep',     var=nstep,     filename='outdat/meanflow.'//iomode//'5')
      call h5srite(varname='nsamples',  var=nsamples,  filename='outdat/meanflow.'//iomode//'5')
      call h5srite(varname='nstep_sbeg',var=nstep_sbeg,filename='outdat/meanflow.'//iomode//'5')
      call h5srite(varname='time_sbeg', var=time_sbeg, filename='outdat/meanflow.'//iomode//'5')
      call h5srite(varname='nstep',     var=nstep,     filename=  'outdat/2ndsta.'//iomode//'5')
      call h5srite(varname='nsamples',  var=nsamples,  filename=  'outdat/2ndsta.'//iomode//'5')
      call h5srite(varname='nstep_sbeg',var=nstep_sbeg,filename=  'outdat/2ndsta.'//iomode//'5')
      call h5srite(varname='time_sbeg', var=time_sbeg, filename=  'outdat/2ndsta.'//iomode//'5')
      call h5srite(varname='nstep',     var=nstep,     filename=  'outdat/3rdsta.'//iomode//'5')
      call h5srite(varname='nsamples',  var=nsamples,  filename=  'outdat/3rdsta.'//iomode//'5')
      call h5srite(varname='nstep_sbeg',var=nstep_sbeg,filename=  'outdat/3rdsta.'//iomode//'5')
      call h5srite(varname='time_sbeg', var=time_sbeg, filename=  'outdat/3rdsta.'//iomode//'5')
      call h5srite(varname='nstep',     var=nstep,     filename=  'outdat/budget.'//iomode//'5')
      call h5srite(varname='nsamples',  var=nsamples,  filename=  'outdat/budget.'//iomode//'5')
      call h5srite(varname='nstep_sbeg',var=nstep_sbeg,filename=  'outdat/budget.'//iomode//'5')
      call h5srite(varname='time_sbeg', var=time_sbeg, filename=  'outdat/budget.'//iomode//'5')
    endif
    !
  end subroutine writemeanflow
  !+-------------------------------------------------------------------+
  !| The end of the subroutine writemeanflow.                          |
  !+-------------------------------------------------------------------+
  !
end module io
!+---------------------------------------------------------------------+
!| The end of the module io.                                           |
!+---------------------------------------------------------------------+