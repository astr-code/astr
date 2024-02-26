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
    use commvar, only : numq,npdci,npdcj,npdck,               &
                        conschm,difschm,lfilter,alfa_filter,hm,turbmode
    use commfunc,only : coeffcompac,ptds_ini,ptdsfilter_ini,           &
                        ptds_aym_ini,genfilt10coef
    use models,  only : init_komegasst
    !
    ! local data
    integer :: nscheme,i
    !
    ! convectional term
    if(conschm(4:4)=='c') then
      ! a compact scheme is used
      !
      read(conschm(1:3),*) nscheme
      !
      alfa_con=coeffcompac(nscheme)
      !
      if(mod(nscheme/100,2)==0) then
        ! symmetrical central scheme
        call ptds_ini(cci,alfa_con,im,npdci)
        call ptds_ini(ccj,alfa_con,jm,npdcj)
        call ptds_ini(cck,alfa_con,km,npdck)
      else
        ! asymmetrical reconstruction upwind scheme
        call ptds_aym_ini(uci,alfa_con,im,npdci,windir='+')
        call ptds_aym_ini(ucj,alfa_con,jm,npdcj,windir='+')
        call ptds_aym_ini(uck,alfa_con,km,npdck,windir='+')
        !
        call ptds_aym_ini(bci,alfa_con,im,npdci,windir='-')
        call ptds_aym_ini(bcj,alfa_con,jm,npdcj,windir='-')
        call ptds_aym_ini(bck,alfa_con,km,npdck,windir='-')
      endif
      !
    endif
    !
    ! diffusional term
    if(difschm(4:4)=='c') then
      ! a compact scheme is used
      !
      read(difschm(1:3),*) nscheme
      !
      alfa_dif=coeffcompac(nscheme)
      !
      call ptds_ini(dci,alfa_dif,im,npdci)
      call ptds_ini(dcj,alfa_dif,jm,npdcj)
      call ptds_ini(dck,alfa_dif,km,npdck)
      !
    endif
    !
    if(lfilter) then
      !
      call ptdsfilter_ini(fci,alfa_filter,im,npdci)
      call ptdsfilter_ini(fcj,alfa_filter,jm,npdcj)
      call ptdsfilter_ini(fck,alfa_filter,km,npdck)
      !
    endif
    !
    call genfilt10coef(alfa_filter)
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
    use commfunc,  only : ddfc
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
      df(:)=ddfc(ff(:),difschm,npdci,im,alfa_dif,dci)
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
      df(:)=ddfc(ff(:),difschm,npdcj,jm,alfa_dif,dcj)
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
        df(:)=ddfc(ff(:),difschm,npdck,km,alfa_dif,dck,lfft=lfftk)
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
    use commfunc,  only : ddfc
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
        df(:,n)=ddfc(ff(:,n),difschm,npdci,im,alfa_dif,dci)
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
          df(:,n)=ddfc(ff(:,n),difschm,npdcj,jm,alfa_dif,dcj)
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
          df(:,n)=ddfc(ff(:,n),difschm,npdck,km,alfa_dif,dck,lfft=lfftk)
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
      if(lio .and. lreport) call timereporter(routine='gradcal', &
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
    use commfunc, only : spafilter10,spafilter6exp
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
        fph(:,n)=spafilter10(phi(:,n),npdci,im,alfa_filter,fci)
        ! fph(:,n)=spafilter6exp(phi(:,n),npdci,im)
      enddo
      !
      q(0:im,j,k,:)=fph(0:im,:)
      !
      ! if(npdci==1) then
      !   q(2:im,j,k,:)=fph(2:im,:)
      ! elseif(npdci==2) then
      !   q(0:im-2,j,k,:)=fph(0:im-2,:)
      ! elseif(npdci==3) then
      !   q(0:im,j,k,:)=fph(0:im,:)
      ! endif
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
          fph(:,n)=spafilter10(phi(:,n),npdcj,jm,alfa_filter,fcj)
          ! fph(:,n)=spafilter6exp(phi(:,n),npdcj,jm)
        enddo
        !
        q(i,0:jm,k,:)=fph(0:jm,:)
        !
        ! if(npdcj==1) then
        !   q(i,2:jm,k,:)=fph(2:jm,:)
        ! elseif(npdcj==2) then
        !   q(i,0:jm-2,k,:)=fph(0:jm-2,:)
        ! elseif(npdcj==3) then
        !   q(i,0:jm,k,:)=fph(0:jm,:)
        ! endif
        !
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
          fph(:,n)=spafilter10(phi(:,n),npdck,km,alfa_filter,fck,lfft=lfftk)
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
    !
    ! call filter2e(q(:,:,:,1))
    ! call filter2e(q(:,:,:,2))
    ! call filter2e(q(:,:,:,3))
    ! call filter2e(q(:,:,:,4))
    ! call filter2e(q(:,:,:,5))
    ! call filter2e(q(:,:,:,6))
    if(trim(turbmode)=='k-omega') then
      call filter2e(q(:,:,:,7))
    endif
    !
    if(present(timerept) .and. timerept) then
      !
      subtime=subtime+ptime()-time_beg
      !
      if(lio .and. lreport) call timereporter(routine='filterq', &
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
    !
    use commvar,  only : is,ie,js,je,ks,ke,numq,                       &
                         spg_i0,spg_im,spg_j0,spg_jm,spg_k0,spg_km,    &
                         lspg_i0,lspg_im,lspg_j0,lspg_jm,lspg_k0,      &
                         lspg_km,spg_im_beg,spg_im_end,spg_jm_beg,     &
                         spg_jm_end
                         
    use commarray,only: lenspg_i0,lenspg_im,lenspg_j0,lenspg_jm,       &
                        lenspg_k0,lenspg_km,xspg_i0,xspg_im,xspg_j0,   &
                        xspg_jm,xspg_k0,xspg_km,x,q
    use commfunc, only : spafilter10
    !
    real(8),parameter :: dampfac=0.05d0
    !
    integer :: i,j,k,n
    real(8) :: dis,var1
    real(8),allocatable :: qtemp(:,:,:,:)
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
          dis=0.d0
          do i=spg_im_beg,spg_im_end
            !
            dis=  (x(i,j,k,1)-xspg_im(j,k,1))**2+               &
                  (x(i,j,k,2)-xspg_im(j,k,2))**2+               &
                  (x(i,j,k,3)-xspg_im(j,k,3))**2                
            !
            var1=dampfac*dis/lenspg_im(j,k)
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
          dis=0.d0
          do j=spg_jm_beg,spg_jm_end
            !
            dis=  (x(i,j,k,1)-xspg_jm(i,k,1))**2+               &
                  (x(i,j,k,2)-xspg_jm(i,k,2))**2+               &
                  (x(i,j,k,3)-xspg_jm(i,k,3))**2                
            !
            var1=dampfac*dis/lenspg_jm(i,k)
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
  end subroutine spongefilter
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