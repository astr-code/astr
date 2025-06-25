!+---------------------------------------------------------------------+
!| This module contains subroutines related to the method of moment.   |
!| ==============                                                      |
!| CHANGE RECORD                                                       |
!| -------------                                                       |
!| 15-08-2023: Created by J. Fang @ STFC Daresbury Laboratory          |
!+---------------------------------------------------------------------+
module comsolver
  !
  use constdef
  use commvar,  only : im,jm,km,hm,deltat,nstep,ndims,lreport,ctime,  &
                       ltimrpt,lfftk
  use parallel, only : lio,ptime,mpirankname,mpistop,mpirank,lio,     &
                       dataswap,irk,jrk,krk,irkm,jrkm,krkm
  use utility,  only : timereporter
  !
  implicit none
  !
  real(8),allocatable :: alfa_con(:),alfa_dif(:)
  real(8), allocatable, dimension(:,:) :: cci,ccj,cck,dci,dcj,dck,     &
                                          fci,fcj,fck,uci,ucj,uck,     &
                                          bci,bcj,bck
  contains
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to initialise solver.                          |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 08-02-2021  | Created by J. Fang @ Warrington                     |
  !+-------------------------------------------------------------------+
  subroutine solvrinit
    !
    use commvar,    only : numq,npdci,npdcj,npdck,               &
                           conschm,difschm,lfilter,alfa_filter,hm,turbmode
    use models,     only : init_komegasst
    use filter,     only : compact_filter_initiate,filter_coefficient_cal, &
                           filter_i,filter_j,filter_k,filter_ii,filter_jj, &
                           filter_kk
    use derivative, only : fd_scheme_initiate,fds_compact_i,fds_compact_j, &
                           fds_compact_k,explicit_central,compact_central,fds
    use flux,       only : compact_flux_initiate,ptds_aym_ini,coeffcompac, &
                           flux_uw_i,flux_dw_i,flux_uw_j,flux_dw_j,flux_uw_k,flux_dw_k

    ! local data
    integer :: nscheme,i

    ! convectional and diffusional term

    call fd_scheme_initiate(asolver=fds_compact_i,scheme=difschm,ntype=npdci,dim=im,dir=1)
    call fd_scheme_initiate(asolver=fds_compact_j,scheme=difschm,ntype=npdcj,dim=jm,dir=2)
    call fd_scheme_initiate(asolver=fds_compact_k,scheme=difschm,ntype=npdck,dim=km,dir=3)

    if(difschm(4:4)=='e') then
      ! a explicit scheme is used

      allocate(explicit_central :: fds)

    elseif(difschm(4:4)=='c') then
      ! a compact scheme is used

      allocate(compact_central :: fds)

    else

      print*,' !! error !! the definition of a scheme must end with c/e, e.g., 642c'

      stop

    endif

    ! Extract the scheme number from the first 3 characters of conschm
    read(conschm(1:3), *) nscheme
    
    ! Check if nscheme is a multiple of 200 (i.e., using central scheme)
    if (mod(nscheme/100, 2) == 0) then
        ! Central scheme is used for the convection term
    
        ! Check for consistency between convection and diffusion schemes
        if (trim(conschm) /= trim(difschm)) then

            print *, ' !! WARNING !! The convection and diffusion terms are both central schemes,'
            print *, ' !! WARNING !! but they are not consistently defined.'
            print *, ' !! WARNING !! It is recommended to use the same scheme for both terms.'
            print *, ' !! WARNING !! Otherwise, unexpected behavior may occur.'

            print*,  ' ** conschm:',conschm
            print*,  ' ** difschm:',difschm

            stop
        end if
    else
      ! upwind-baised schemes
      if(conschm(4:4)=='c') then

        call compact_flux_initiate(asolver=flux_uw_i,scheme=nscheme,ntype=npdci,dim=im,wind='+')
        call compact_flux_initiate(asolver=flux_dw_i,scheme=nscheme,ntype=npdci,dim=im,wind='-')
        call compact_flux_initiate(asolver=flux_uw_j,scheme=nscheme,ntype=npdcj,dim=jm,wind='+')
        call compact_flux_initiate(asolver=flux_dw_j,scheme=nscheme,ntype=npdcj,dim=jm,wind='-')
        call compact_flux_initiate(asolver=flux_uw_k,scheme=nscheme,ntype=npdck,dim=km,wind='+')
        call compact_flux_initiate(asolver=flux_dw_k,scheme=nscheme,ntype=npdck,dim=km,wind='-')

        ! alfa_con=coeffcompac(nscheme)

        ! call ptds_aym_ini(uci,alfa_con,im,npdci,windir='+')
        ! call ptds_aym_ini(ucj,alfa_con,jm,npdcj,windir='+')
        ! call ptds_aym_ini(uck,alfa_con,km,npdck,windir='+')
        ! !
        ! call ptds_aym_ini(bci,alfa_con,im,npdci,windir='-')
        ! call ptds_aym_ini(bcj,alfa_con,jm,npdcj,windir='-')
        ! call ptds_aym_ini(bck,alfa_con,km,npdck,windir='-')

      endif

    end if

    !
    if(lfilter) then
      call compact_filter_initiate(afilter=filter_i,ntype=npdci,dim=im)
      call compact_filter_initiate(afilter=filter_j,ntype=npdcj,dim=jm)
      call compact_filter_initiate(afilter=filter_k,ntype=npdck,dim=km)
      ! call compact_filter_initiate(afilter=filter_ii,ntype=npdci,dim=im,note='boundary_no_filter')
      ! call compact_filter_initiate(afilter=filter_jj,ntype=npdcj,dim=jm,note='boundary_no_filter')
      ! call compact_filter_initiate(afilter=filter_kk,ntype=npdck,dim=km,note='boundary_no_filter')
    endif
    !
    call filter_coefficient_cal(alfa_filter)
    !
    if(trim(turbmode)=='k-omega') then
      call init_komegasst
    endif
    !
    if(lio) print*,' ** numerical solver initilised.'
    !
  end subroutine solvrinit
  !+-------------------------------------------------------------------+
  !| The end of the subroutine solvrinit.                              |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is a general gradient calculater                  |
  !|   input: scalar                                                   |
  !|   output: the gradient of the input scalar                        |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 22-07-2022  | Moved from  diffrsdcal6 by J. Fang @ Warrington     |
  !+-------------------------------------------------------------------+
  function grad(var) result(dvar)
    !
    use commvar,   only : im,jm,km,npdci,npdcj,npdck,difschm,ndims
    use commarray, only : dxi
    use derivative, only : fds,fds_compact_i,fds_compact_j,fds_compact_k
    !
    ! arguments
    real(8),intent(in) :: var(-hm:im+hm,-hm:jm+hm,-hm:km+hm)
    real(8) :: dvar(0:im,0:jm,0:km,1:3)
    !
    ! local data
    integer :: i,j,k
    real(8),allocatable :: df(:),ff(:)
    !
    allocate(ff(-hm:im+hm),df(0:im))
    !
    dvar=0.d0
    !
    do k=0,km
    do j=0,jm
      !
      ff(:)=var(:,j,k)
      !
      df(:)=fds%central(fds_compact_i,f=ff(:),dim=im)
      !
      dvar(:,j,k,1)=dvar(:,j,k,1)+df(:)*dxi(0:im,j,k,1,1)
      dvar(:,j,k,2)=dvar(:,j,k,2)+df(:)*dxi(0:im,j,k,1,2)
      dvar(:,j,k,3)=dvar(:,j,k,3)+df(:)*dxi(0:im,j,k,1,3)
      !
    enddo
    enddo
    !
    deallocate(ff,df)
    !
    allocate(ff(-hm:jm+hm),df(0:jm))
    do k=0,km
    do i=0,im
      !
      ff(:)=var(i,:,k)
      !
      df(:)=fds%central(fds_compact_j,f=ff(:),dim=jm)
      !
      dvar(i,:,k,1)=dvar(i,:,k,1)+df(:)*dxi(i,0:jm,k,2,1)
      dvar(i,:,k,2)=dvar(i,:,k,2)+df(:)*dxi(i,0:jm,k,2,2)
      dvar(i,:,k,3)=dvar(i,:,k,3)+df(:)*dxi(i,0:jm,k,2,3)
      !
    enddo
    enddo
    deallocate(ff,df)
    !
    if(ndims==3) then
      !
      allocate(ff(-hm:km+hm),df(0:km))
      do j=0,jm
      do i=0,im
        !
        ff(:)=var(i,j,:)
        !
        df(:)=fds%central(fds_compact_k,f=ff(:),dim=km)
        !
        dvar(i,j,:,1)=dvar(i,j,:,1)+df(:)*dxi(i,j,0:km,3,1)
        dvar(i,j,:,2)=dvar(i,j,:,2)+df(:)*dxi(i,j,0:km,3,2)
        dvar(i,j,:,3)=dvar(i,j,:,3)+df(:)*dxi(i,j,0:km,3,3)
        !
      enddo
      enddo
      deallocate(ff,df)
      !
    endif
    !
    return
    !
  end function grad
  !+-------------------------------------------------------------------+
  !| The end of the subroutine grad.                                   |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to calculate gradients of flow variables.      |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 08-10-2021  | Moved from  diffrsdcal6 by J. Fang @ Warrington     |
  !+-------------------------------------------------------------------+
  subroutine gradcal(timerept)
    !
    use commvar,   only : im,jm,km,npdci,npdcj,npdck,difschm,ndims,    &
                          num_species,num_modequ,is,ie,js,je,ks,ke,    &
                          turbmode
    use commarray, only : vel,tmp,spc,dvel,dtmp,dspc,dxi,omg,tke,dtke, &
                          domg
    use derivative, only : fds,fds_compact_i,fds_compact_j,fds_compact_k
    !
    ! arguments
    logical,intent(in),optional :: timerept
    !
    ! local data
    integer :: i,j,k,n,ncolm
    real(8),allocatable :: df(:,:),ff(:,:)
    !
    real(8) :: time_beg
    real(8),save :: subtime=0.d0
    !
    if(present(timerept) .and. timerept) time_beg=ptime() 
    !
    dvel=0.d0
    dtmp=0.d0
    dspc=0.d0
    !
    if(trim(turbmode)=='k-omega') then
      dtke=0.d0
      domg=0.d0
    endif
    !
    ncolm=4+num_species+num_modequ
    !
    ! calculate velocity and temperature gradient
    !
    allocate(ff(-hm:im+hm,ncolm),df(0:im,ncolm))
    !
    do k=0,km
    do j=0,jm
      !
      ff(:,1)=vel(:,j,k,1)
      ff(:,2)=vel(:,j,k,2)
      ff(:,3)=vel(:,j,k,3)
      ff(:,4)=tmp(:,j,k)
      !
      if(num_species>0) then
        do n=1,num_species
          ff(:,4+n)=spc(:,j,k,n)
        enddo
      endif
      !
      if(trim(turbmode)=='k-omega') then
        n=4+num_species
        !
        ff(:,n+1)=tke(:,j,k)
        ff(:,n+2)=omg(:,j,k)
      endif
      !
      do n=1,ncolm
        df(:,n)=fds%central(fds_compact_i,f=ff(:,n),dim=im)
      enddo
      !
      dvel(:,j,k,1,1)=dvel(:,j,k,1,1)+df(:,1)*dxi(0:im,j,k,1,1)
      dvel(:,j,k,1,2)=dvel(:,j,k,1,2)+df(:,1)*dxi(0:im,j,k,1,2)
      dvel(:,j,k,1,3)=dvel(:,j,k,1,3)+df(:,1)*dxi(0:im,j,k,1,3)
      !
      dvel(:,j,k,2,1)=dvel(:,j,k,2,1)+df(:,2)*dxi(0:im,j,k,1,1)
      dvel(:,j,k,2,2)=dvel(:,j,k,2,2)+df(:,2)*dxi(0:im,j,k,1,2)
      dvel(:,j,k,2,3)=dvel(:,j,k,2,3)+df(:,2)*dxi(0:im,j,k,1,3)
      !
      dvel(:,j,k,3,1)=dvel(:,j,k,3,1)+df(:,3)*dxi(0:im,j,k,1,1)
      dvel(:,j,k,3,2)=dvel(:,j,k,3,2)+df(:,3)*dxi(0:im,j,k,1,2)
      dvel(:,j,k,3,3)=dvel(:,j,k,3,3)+df(:,3)*dxi(0:im,j,k,1,3)
      !
      dtmp(:,j,k,1)=dtmp(:,j,k,1)+df(:,4)*dxi(0:im,j,k,1,1)
      dtmp(:,j,k,2)=dtmp(:,j,k,2)+df(:,4)*dxi(0:im,j,k,1,2)
      dtmp(:,j,k,3)=dtmp(:,j,k,3)+df(:,4)*dxi(0:im,j,k,1,3)
      !
      if(num_species>0) then
        do n=1,num_species
          dspc(:,j,k,n,1)=dspc(:,j,k,n,1)+df(:,4+n)*dxi(0:im,j,k,1,1)
          dspc(:,j,k,n,2)=dspc(:,j,k,n,2)+df(:,4+n)*dxi(0:im,j,k,1,2)
          dspc(:,j,k,n,3)=dspc(:,j,k,n,3)+df(:,4+n)*dxi(0:im,j,k,1,3)
        enddo
      endif
      !
      if(trim(turbmode)=='k-omega') then
        n=4+num_species
        !
        dtke(:,j,k,1)=dtke(:,j,k,1)+df(:,1+n)*dxi(0:im,j,k,1,1)
        dtke(:,j,k,2)=dtke(:,j,k,2)+df(:,1+n)*dxi(0:im,j,k,1,2)
        dtke(:,j,k,3)=dtke(:,j,k,3)+df(:,1+n)*dxi(0:im,j,k,1,3)
        !
        domg(:,j,k,1)=domg(:,j,k,1)+df(:,2+n)*dxi(0:im,j,k,1,1)
        domg(:,j,k,2)=domg(:,j,k,2)+df(:,2+n)*dxi(0:im,j,k,1,2)
        domg(:,j,k,3)=domg(:,j,k,3)+df(:,2+n)*dxi(0:im,j,k,1,3)
      endif
      !
    enddo
    enddo
    !
    deallocate(ff,df)
    !
    if(ndims>=2) then
      !
      allocate(ff(-hm:jm+hm,ncolm),df(0:jm,ncolm))
      do k=0,km
      do i=0,im
        !
        ff(:,1)=vel(i,:,k,1)
        ff(:,2)=vel(i,:,k,2)
        ff(:,3)=vel(i,:,k,3)
        ff(:,4)=tmp(i,:,k)
        !
        if(num_species>0) then
          do n=1,num_species
            ff(:,4+n)=spc(i,:,k,n)
          enddo
        endif
        !
        if(trim(turbmode)=='k-omega') then
          n=4+num_species
          !
          ff(:,n+1)=tke(i,:,k)
          ff(:,n+2)=omg(i,:,k)
        endif
        !
        do n=1,ncolm
          df(:,n)=fds%central(fds_compact_j,f=ff(:,n),dim=jm)
        enddo
        !
        dvel(i,:,k,1,1)=dvel(i,:,k,1,1)+df(:,1)*dxi(i,0:jm,k,2,1)
        dvel(i,:,k,1,2)=dvel(i,:,k,1,2)+df(:,1)*dxi(i,0:jm,k,2,2)
        dvel(i,:,k,1,3)=dvel(i,:,k,1,3)+df(:,1)*dxi(i,0:jm,k,2,3)
        !
        dvel(i,:,k,2,1)=dvel(i,:,k,2,1)+df(:,2)*dxi(i,0:jm,k,2,1)
        dvel(i,:,k,2,2)=dvel(i,:,k,2,2)+df(:,2)*dxi(i,0:jm,k,2,2)
        dvel(i,:,k,2,3)=dvel(i,:,k,2,3)+df(:,2)*dxi(i,0:jm,k,2,3)
        !
        dvel(i,:,k,3,1)=dvel(i,:,k,3,1)+df(:,3)*dxi(i,0:jm,k,2,1)
        dvel(i,:,k,3,2)=dvel(i,:,k,3,2)+df(:,3)*dxi(i,0:jm,k,2,2)
        dvel(i,:,k,3,3)=dvel(i,:,k,3,3)+df(:,3)*dxi(i,0:jm,k,2,3)
        !
        dtmp(i,:,k,1)=dtmp(i,:,k,1)+df(:,4)*dxi(i,0:jm,k,2,1)
        dtmp(i,:,k,2)=dtmp(i,:,k,2)+df(:,4)*dxi(i,0:jm,k,2,2)
        dtmp(i,:,k,3)=dtmp(i,:,k,3)+df(:,4)*dxi(i,0:jm,k,2,3)
        !
        if(num_species>0) then
          do n=1,num_species
            dspc(i,:,k,n,1)=dspc(i,:,k,n,1)+df(:,4+n)*dxi(i,0:jm,k,2,1)
            dspc(i,:,k,n,2)=dspc(i,:,k,n,2)+df(:,4+n)*dxi(i,0:jm,k,2,2)
            dspc(i,:,k,n,3)=dspc(i,:,k,n,3)+df(:,4+n)*dxi(i,0:jm,k,2,3)
          enddo
        endif
        !
        if(trim(turbmode)=='k-omega') then
          n=4+num_species
          !
          dtke(i,:,k,1)=dtke(i,:,k,1)+df(:,1+n)*dxi(i,0:jm,k,2,1)
          dtke(i,:,k,2)=dtke(i,:,k,2)+df(:,1+n)*dxi(i,0:jm,k,2,2)
          dtke(i,:,k,3)=dtke(i,:,k,3)+df(:,1+n)*dxi(i,0:jm,k,2,3)
          !
          domg(i,:,k,1)=domg(i,:,k,1)+df(:,2+n)*dxi(i,0:jm,k,2,1)
          domg(i,:,k,2)=domg(i,:,k,2)+df(:,2+n)*dxi(i,0:jm,k,2,2)
          domg(i,:,k,3)=domg(i,:,k,3)+df(:,2+n)*dxi(i,0:jm,k,2,3)
          !
          !
        endif
        !
      enddo
      enddo
      deallocate(ff,df)
      !
    endif
    !
    if(ndims==3) then
      allocate(ff(-hm:km+hm,ncolm),df(0:km,ncolm))
      do j=0,jm
      do i=0,im
        !
        ff(:,1)=vel(i,j,:,1)
        ff(:,2)=vel(i,j,:,2)
        ff(:,3)=vel(i,j,:,3)
        ff(:,4)=tmp(i,j,:)
        !
        if(num_species>0) then
          do n=1,num_species
            ff(:,4+n)=spc(i,j,:,n)
          enddo
        endif
        !
        if(trim(turbmode)=='k-omega') then
          n=4+num_species
          !
          ff(:,n+1)=tke(i,j,:)
          ff(:,n+2)=omg(i,j,:)
        endif
        !
        do n=1,ncolm
          df(:,n)=fds%central(fds_compact_k,f=ff(:,n),dim=km)
        enddo
        !
        dvel(i,j,:,1,1)=dvel(i,j,:,1,1)+df(:,1)*dxi(i,j,0:km,3,1)
        dvel(i,j,:,1,2)=dvel(i,j,:,1,2)+df(:,1)*dxi(i,j,0:km,3,2)
        dvel(i,j,:,1,3)=dvel(i,j,:,1,3)+df(:,1)*dxi(i,j,0:km,3,3)
        !
        dvel(i,j,:,2,1)=dvel(i,j,:,2,1)+df(:,2)*dxi(i,j,0:km,3,1)
        dvel(i,j,:,2,2)=dvel(i,j,:,2,2)+df(:,2)*dxi(i,j,0:km,3,2)
        dvel(i,j,:,2,3)=dvel(i,j,:,2,3)+df(:,2)*dxi(i,j,0:km,3,3)
        !
        dvel(i,j,:,3,1)=dvel(i,j,:,3,1)+df(:,3)*dxi(i,j,0:km,3,1)
        dvel(i,j,:,3,2)=dvel(i,j,:,3,2)+df(:,3)*dxi(i,j,0:km,3,2)
        dvel(i,j,:,3,3)=dvel(i,j,:,3,3)+df(:,3)*dxi(i,j,0:km,3,3)
        !
        dtmp(i,j,:,1)=dtmp(i,j,:,1)+df(:,4)*dxi(i,j,0:km,3,1)
        dtmp(i,j,:,2)=dtmp(i,j,:,2)+df(:,4)*dxi(i,j,0:km,3,2)
        dtmp(i,j,:,3)=dtmp(i,j,:,3)+df(:,4)*dxi(i,j,0:km,3,3)
        !
        if(num_species>0) then
          do n=1,num_species
            dspc(i,j,:,n,1)=dspc(i,j,:,n,1)+df(:,4+n)*dxi(i,j,0:km,3,1)
            dspc(i,j,:,n,2)=dspc(i,j,:,n,2)+df(:,4+n)*dxi(i,j,0:km,3,2)
            dspc(i,j,:,n,3)=dspc(i,j,:,n,3)+df(:,4+n)*dxi(i,j,0:km,3,3)
          enddo
        endif
        !
        if(trim(turbmode)=='k-omega') then
          n=4+num_species
          !
          dtke(i,j,:,1)=dtke(i,j,:,1)+df(:,1+n)*dxi(i,j,0:km,3,1)
          dtke(i,j,:,2)=dtke(i,j,:,2)+df(:,1+n)*dxi(i,j,0:km,3,2)
          dtke(i,j,:,3)=dtke(i,j,:,3)+df(:,1+n)*dxi(i,j,0:km,3,3)
          !
          domg(i,j,:,1)=domg(i,j,:,1)+df(:,2+n)*dxi(i,j,0:km,3,1)
          domg(i,j,:,2)=domg(i,j,:,2)+df(:,2+n)*dxi(i,j,0:km,3,2)
          domg(i,j,:,3)=domg(i,j,:,3)+df(:,2+n)*dxi(i,j,0:km,3,3)
        endif
        !
      enddo
      enddo
      deallocate(ff,df)
    endif
    !
    if(present(timerept) .and. timerept) then
      !
      subtime=subtime+ptime()-time_beg
      !
      if(lio .and. lreport .and. ltimrpt) call timereporter(routine='gradcal', &
                                             timecost=subtime,  &
                                              message='calculation of gradients')
    endif
    !
    return
    !
  end subroutine gradcal
  !+-------------------------------------------------------------------+
  !| The end of the subroutine gradcal.                                |
  !+-------------------------------------------------------------------+
  !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! This subroutine is used for spatial filter the conservative variable
  ! for stabilizing the computation.
  ! 10-order filter is incorporated.
  ! for boundary filter: the high-order one side filter is used.
  ! the 0-6-6-6-8-10.............-10-8-6-6-6-0. boundary order is
  ! dopted.
  ! Ref: Datta V. Gaitonde and Miguel R. Visbal, AIAA JOURNAL Vol.38,
  !      No.11, November 2000.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Writen by Fang Jian, 2008-11-03.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine filterq(timerept)
    !
    use commvar,  only : im,jm,km,numq,npdci,npdcj,npdck,              &
                         alfa_filter,ndims,is,ie,js,je,ks,ke,turbmode
    use commarray,only : q
    use filter,   only : compact_filter,filter_i,filter_j,filter_k,  &
                        filter_ii,filter_jj,filter_kk
    !
    ! arguments
    logical,intent(in),optional :: timerept
    !
    ! local data
    integer :: i,j,k,n,m
    real(8),allocatable :: phi(:,:),fph(:,:)
    !
    real(8) :: time_beg
    real(8),save :: subtime=0.d0
    !
    if(present(timerept) .and. timerept) time_beg=ptime() 
    !
    ! filtering in i direction
    call dataswap(q,direction=1,timerept=ltimrpt)
    !
    allocate(phi(-hm:im+hm,1:numq),fph(0:im,1:numq))
    !
    do k=0,km
    do j=0,jm
      !
      phi(:,:)=q(:,j,k,:)
      !
      do n=1,numq
        fph(:,n)=compact_filter(afilter=filter_i,f=phi(:,n),dim=im)
      enddo
      !
      q(0:im,j,k,:)=fph(0:im,:)
      !
    end do
    end do
    !
    deallocate(phi,fph)
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! end filter in i direction.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    if(ndims>=2) then
      !
      ! filtering in j direction
      call dataswap(q,direction=2,timerept=ltimrpt)
      !
      allocate(phi(-hm:jm+hm,1:numq),fph(0:jm,1:numq))
      !
      do k=0,km
      do i=0,im
        !
        phi(:,:)=q(i,:,k,:)
        !
        do n=1,numq
          fph(:,n)=compact_filter(afilter=filter_j,f=phi(:,n),dim=jm)
        enddo
        !
        q(i,0:jm,k,:)=fph(0:jm,:)
        !
      end do
      end do
      !
      deallocate(phi,fph)
      !
    endif
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! end filter in j direction.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    if(ndims==3) then
      !
      call dataswap(q,direction=3,timerept=ltimrpt)
      !
      !
      allocate(phi(-hm:km+hm,1:numq),fph(0:km,1:numq))
      !
      ! filtering in k direction
      do j=0,jm
      do i=0,im
        !
        phi(:,:)=q(i,j,:,:)
        !
        do n=1,numq
          fph(:,n)=compact_filter(afilter=filter_k,f=phi(:,n),dim=km)
        enddo
        !
        q(i,j,0:km,:)=fph
        !
      end do
      end do
      !
      deallocate(phi,fph)
      !
    end if
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! end filter in k direction.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    if(trim(turbmode)=='k-omega') then
      call filter2e(q(:,:,:,7))
    endif
    !
    if(present(timerept) .and. timerept) then
      !
      subtime=subtime+ptime()-time_beg
      !
      if(lio .and. lreport .and. ltimrpt) call timereporter(routine='filterq', &
                                             timecost=subtime, &
                                              message='low-pass filter')
    endif
    !
    return
    !
  end subroutine filterq
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! End of the subroutine filterq.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  subroutine filter2e(phi)
    !
    use commvar,  only : is,ie,je,js,je,ks,ke,im,jm,km
    !
    real(8),intent(inout) :: phi(-hm:im+hm,-hm:jm+hm,-hm:km+hm) 
    !
    integer :: i,j,k
    real(8),allocatable :: phtemp(:,:,:)
    !
    call dataswap(phi,timerept=ltimrpt)
    !
    allocate(phtemp(is:ie,js:je,ks:ke))
    do k=ks,ke
    do i=is,ie
      !
      do j=js,je
        phtemp(i,j,k)=0.01d0*(0.25d0*(phi(i,j-1,k)+phi(i,j+1,k))+0.5d0*phi(i,j,k))  + &
                      0.99d0*phi(i,j,k)
      enddo
      !
    enddo
    enddo
    !
    phi(is:ie,js:je,ks:ke)=phtemp(is:ie,js:je,ks:ke)
    !
    return
    !
  end subroutine filter2e
  !
  subroutine filter4e(phi)
    !
    use commvar,  only : is,ie,je,js,je,ks,ke,im,jm,km
    !
    real(8),intent(inout) :: phi(-hm:im+hm,-hm:jm+hm,-hm:km+hm) 
    !
    integer :: i,j,k
    real(8),allocatable :: phtemp(:,:,:)
    !
    call dataswap(phi,timerept=ltimrpt)
    !
    allocate(phtemp(is:ie,js:je,ks:ke))
    do k=ks,ke
    do i=is,ie
      !
      do j=js+1,je-1
        phtemp(i,j,k)= 0.0001d0*(0.625d0*phi(i,j,k)               + &
                               0.25d0*(phi(i,j-1,k)+phi(i,j+1,k)) - &
                              0.06250*(phi(i,j-2,k)+phi(i,j+3,k)))+ &
                       0.9999d0*phi(i,j,k)
  
      enddo
      !
    enddo
    enddo
    !
    phi(is:ie,js:je,ks:ke)=phtemp(is:ie,js:je,ks:ke)
    !
    return
    !
  end subroutine filter4e
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! This subroutine is used for spatial filter the conservative variable
  ! for stabilizing the computation.
  ! 2-order Explicit filter is incorporated.
  ! Ref1: Datta V. Gaitonde and Miguel R. Visbal, AIAA JOURNAL Vol.38,
  !      No.11, November 2000. 
  ! Ref2: Xavier Gloerfelt and Philippe Lafon, Computers & Fluids, 2008,
  !       37: 388-401.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine spongefilter
    
    use commvar,  only : spg_def

    if(spg_def=='layer') then
      call spongefilter_layer
    elseif(spg_def=='circl') then
      call spongefilter_global
    endif

  end subroutine spongefilter
  !
  subroutine spongefilter_layer
    !
    use commvar,  only : is,ie,js,je,ks,ke,numq,                       &
                         spg_i0,spg_im,spg_j0,spg_jm,spg_k0,spg_km,    &
                         lspg_i0,lspg_im,lspg_j0,lspg_jm,lspg_k0,      &
                         lspg_km,spg_i0_beg,spg_i0_end,spg_im_beg,     &
                         spg_im_end,spg_jm_beg,spg_jm_end
                         
    use commarray,only: sponge_damp_coef_i0,sponge_damp_coef_im, &
                        sponge_damp_coef_j0,sponge_damp_coef_jm, &
                        sponge_damp_coef_k0,sponge_damp_coef_km,x,q
    !
    integer :: i,j,k,n
    real(8) :: var1
    real(8),allocatable :: qtemp(:,:,:,:)
    !
    ! sponger layer attached at the im end.
    if(lspg_i0) then
      !
      call dataswap(q,direction=1,timerept=ltimrpt)
      !
      if(spg_i0_beg>=0) then
        !
        print*,mpirank,'|',spg_i0_beg,spg_i0_end
        !
        allocate(qtemp(spg_i0_beg:spg_i0_end,js:je,ks:ke,1:numq))
        !
        do k=ks,ke
        do j=js,je
          !
          do i=spg_i0_beg,spg_i0_end
            !
            var1=sponge_damp_coef_i0(i,j,k)
            !
            do n=1,numq
              qtemp(i,j,k,n)=(1.d0-var1)*q(i,j,k,n)+        &
                             num1d6*var1*(q(i+1,j,k,n)+     &
                                          q(i-1,j,k,n)+     &
                                          q(i,j+1,k,n)+     &
                                          q(i,j-1,k,n)+     &
                                          q(i,j,k+1,n)+     &
                                          q(i,j,k-1,n)      )
            enddo
            !
            ! if(j==0) print*,i,sponge_damp_coef_i0(i,j,k)
            !
          enddo
          !
        enddo
        enddo
        !
        do k=ks,ke
        do j=js,je
          !
          do i=spg_i0_beg,spg_i0_end
            !
            do n=1,numq
              q(i,j,k,n)=qtemp(i,j,k,n)
            enddo
            !
          enddo
          !
        enddo
        enddo
        !
        deallocate(qtemp)
        !
      endif
      !
    endif
    !
    ! sponger layer attached at the im end.
    if(lspg_im) then
      !
      call dataswap(q,direction=1,timerept=ltimrpt)
      !
      if(spg_im_beg>=0) then
        !
        allocate(qtemp(spg_im_beg:spg_im_end,js:je,ks:ke,1:numq))
        !
        do k=ks,ke
        do j=js,je
          !
          do i=spg_im_beg,spg_im_end
            !
            var1=sponge_damp_coef_im(i,j,k)
            !
            do n=1,numq
              qtemp(i,j,k,n)=(1.d0-var1)*q(i,j,k,n)+        &
                             num1d6*var1*(q(i+1,j,k,n)+     &
                                          q(i-1,j,k,n)+     &
                                          q(i,j+1,k,n)+     &
                                          q(i,j-1,k,n)+     &
                                          q(i,j,k+1,n)+     &
                                          q(i,j,k-1,n)      )
            enddo
            !
          enddo
          !
        enddo
        enddo
        !
        do k=ks,ke
        do j=js,je
          !
          do i=spg_im_beg,spg_im_end
            !
            do n=1,numq
              q(i,j,k,n)=qtemp(i,j,k,n)
            enddo
            !
          enddo
          !
        enddo
        enddo
        !
        deallocate(qtemp)
        !
      endif
      !
    endif
    !
    ! sponger layer attached at the jm end.
    if(lspg_jm) then
      !
      call dataswap(q,direction=1,timerept=ltimrpt)
      !
      if(spg_jm_beg>=0) then
        !
        allocate(qtemp(is:ie,spg_jm_beg:spg_jm_end,ks:ke,1:numq))
        !
        do k=ks,ke
        do i=is,ie
          !
          do j=spg_jm_beg,spg_jm_end
            !
            var1=sponge_damp_coef_jm(i,j,k)
            !
            do n=1,numq
              qtemp(i,j,k,n)=(1.d0-var1)*q(i,j,k,n)+        &
                             num1d6*var1*(q(i+1,j,k,n)+     &
                                          q(i-1,j,k,n)+     &
                                          q(i,j+1,k,n)+     &
                                          q(i,j-1,k,n)+     &
                                          q(i,j,k+1,n)+     &
                                          q(i,j,k-1,n)      )
            enddo
            !
          enddo
          !
        enddo
        enddo
        !
        do k=ks,ke
        do i=is,ie
          !
          do j=spg_jm_beg,spg_jm_end
            !
            do n=1,numq
              q(i,j,k,n)=qtemp(i,j,k,n)
            enddo
            !
          enddo
          !
        enddo
        enddo
        !
        deallocate(qtemp)
        !
      endif
      !
    endif
    !
  end subroutine spongefilter_layer

  subroutine spongefilter_global
    !
    use commvar,  only : is,ie,js,je,ks,ke,numq,lsponge,lsponge_loc
    use commarray,only: sponge_damp_coef,x,q
    !
    integer :: i,j,k,n
    real(8) :: var1
    real(8),allocatable :: qtemp(:,:,:,:)
    !
    ! sponger layer attached at the im end.

    if(lsponge) call dataswap(q)

    if(lsponge_loc) then

      allocate(qtemp(is:ie,js:je,ks:ke,1:numq))
  
      do k=ks,ke
      do j=js,je
      do i=is,ie
        var1=sponge_damp_coef(i,j,k)

        do n=1,numq
          qtemp(i,j,k,n)=(1.d0-var1)*q(i,j,k,n)+        &
                         num1d6*var1*(q(i+1,j,k,n)+     &
                                      q(i-1,j,k,n)+     &
                                      q(i,j+1,k,n)+     &
                                      q(i,j-1,k,n)+     &
                                      q(i,j,k+1,n)+     &
                                      q(i,j,k-1,n)      )
        enddo

      enddo
      enddo
      enddo

      q(is:ie,js:je,ks:ke,1:numq)=qtemp(is:ie,js:je,ks:ke,1:numq)

      deallocate(qtemp)

    endif

  end subroutine spongefilter_global
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! End of the subroutine spongefilter.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  subroutine check_mat55_unit(matrix,normal)
    !
    real(8),intent(in) :: matrix(5,5)
    logical,intent(out) :: normal
    !
    real(8) :: epslion
    integer :: i,j
    !
    epslion=1.d-8
    !
    normal=.true.
    !
    do j=1,5
    do i=1,5
     if( i==j ) then
       if(abs(matrix(i,j)-1.d0)<epslion) then
        continue
       else
         print*,' !! WARNING of UNIT MARTIX'
         print*, i,j,matrix(i,j)
         normal=.false.
       endif
     else
       if(abs(matrix(i,j))<epslion) then
        continue
       else
         print*,' !! WARNING of UNIT MARTIX'
         print*, i,j,matrix(i,j)
         normal=.false.
       endif
     endif
    enddo
    enddo
    !
  end subroutine check_mat55_unit
  !
end module comsolver
!+---------------------------------------------------------------------+
!| The end of the module comsolver.                                    |
!+---------------------------------------------------------------------+
