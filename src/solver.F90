!+---------------------------------------------------------------------+
!| This module contains some solver related subroutines.               |
!+---------------------------------------------------------------------+
!| CHANGE RECORD                                                       |
!| -------------                                                       |
!| 08-02-2021  | Created by J. Fang                                    |
!+---------------------------------------------------------------------+
module solver
  !
  use constdef
  use parallel, only : mpirankname,mpistop,mpirank,lio,dataswap
  use commvar,  only : ndims,ks,ke,hm
  !
  implicit none
  !
  real(8) :: alfa_con(3),alfa_dif(3)
  real(8), allocatable, dimension(:,:) :: cci,ccj,cck,dci,dcj,dck,     &
                                          fci,fcj,fck
  !
  contains
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to calculate some constant parameters          |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 08-02-2021  | Created by J. Fang @ Warrington                     |
  !+-------------------------------------------------------------------+
  subroutine refcal
    !
    use commvar, only : numq,num_species,prandtl,gamma,rgas,ia,ja,ka,  &
                        uinf,vinf,winf,roinf,pinf,tinf,const1,const2,  &
                        const3,const4,const5,const6,const7,tempconst,  &
                        tempconst1,reynolds,ref_t,mach
    !
    numq=5+num_species
    !
    if(ia>0 .and. ja>0 .and. ka>0) then
      ndims=3
      !
      if(ia<hm .or. ja<hm .or. ka<hm) then
        print*,' !! input dimension smaller than the halo cells'
      endif
      !
    elseif(ka==0 .and. ja==0 .and. ia==0) then
      ndims=0
    elseif(ka==0 .and. ja==0) then
      ndims=1
      !
      if(ia<hm) then
        print*,' !! input dimension smaller than the halo cells'
      endif
      !
    elseif(ka==0) then
      ndims=2
      !
      if(ia<hm .or. ja<hm) then
        print*,' !! input dimension smaller than the halo cells'
      endif
      !
    else
      print*,ndims
      stop ' !! ndims error @ refcal'
    endif
    !
    prandtl=0.72d0
    gamma=1.4d0
    rgas=287.1d0
    !
    const1=1.d0/(gamma*(gamma-1.d0)*mach**2)
    const2=gamma*mach**2
    const3=(gamma-1.d0)/3.d0*prandtl*(mach**2)
    const4=(gamma-1.d0)*mach**2*reynolds*prandtl
    const5=(gamma-1.d0)*mach**2
    const6=1.d0/(gamma-1.d0)
    const7=(gamma-1.d0)*mach**2*Reynolds*prandtl
    !
    uinf=1.d0
    vinf=0.d0
    winf=0.d0
    tinf=1.d0
    roinf=1.d0
    !
    pinf=roinf*tinf/const2
    !
    tempconst=110.3d0/ref_t
    tempconst1=1.d0+tempconst
    !
  end subroutine refcal
  !+-------------------------------------------------------------------+
  !| The end of the subroutine refcal.                                 |
  !+-------------------------------------------------------------------+
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
    use commvar, only : im,jm,km,numq,npdci,npdcj,npdck,               &
                        conschm,difschm,lfilter,alfa_filter,hm
    use commfunc,only : coeffcompac,ptds_ini,ptdsfilter_ini,           &
                        genfilt10coef
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
      call ptds_ini(cci,alfa_con,nscheme,im,npdci)
      call ptds_ini(ccj,alfa_con,nscheme,jm,npdcj)
      call ptds_ini(cck,alfa_con,nscheme,km,npdck)
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
      call ptds_ini(dci,alfa_dif,nscheme,im,npdci)
      call ptds_ini(dcj,alfa_dif,nscheme,jm,npdcj)
      call ptds_ini(dck,alfa_dif,nscheme,km,npdck)
      !
    endif
    !
    if(lfilter) then
      !
      call ptdsfilter_ini(fci,alfa_filter,im,npdci)
      call ptdsfilter_ini(fcj,alfa_filter,jm,npdcj)
      call ptdsfilter_ini(fck,alfa_filter,km,npdck)
      !
      call genfilt10coef(alfa_filter)
      !
    endif
    !
    if(lio) print*,' ** numerical solver initilised.'
    !
  end subroutine solvrinit
  !+-------------------------------------------------------------------+
  !| The end of the subroutine refcal.                                 |
  !+-------------------------------------------------------------------+
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! This subtoutine is used to calculate geometrical transform matrix
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Writen by Fang Jian, 2008-09-22.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine geomcal
    !
    use commvar,   only : ia,ja,ka,im,jm,km,hm,npdci,npdcj,npdck,hm,   &
                          xmax,xmin,ymax,ymin,zmax,zmin,voldom
    use commarray, only : x,jacob,dxi,celvol
    use parallel,  only : gridsendrecv,ksize,psum,pmax,pmin
    use commfunc,  only : coeffcompac,ptds_ini,ddfc,volhex,arquad
    use tecio
    !
    ! local data
    character(len=4) :: cscheme
    integer :: nscheme
    integer :: i,j,k,m
    real(8) :: alfa(3)
    real(8), allocatable, dimension(:,:) :: gci,gcj,gck
    real(8), allocatable :: dx(:,:,:,:,:)
    real(8),allocatable :: phi(:),can(:,:,:,:)
    real(8) :: can1av,can2av,can3av,can1var,can2var,can3var
    !
    allocate( dx(-hm:im+hm,-hm:jm+hm,-hm:km+hm,1:3,1:3) )
    !
    call gridsendrecv
    !
    cscheme='642c'
    !
    read(cscheme(1:3),*) nscheme
    !
    if(cscheme(4:4)=='c') then
      ! a compact scheme is used
      !
      alfa=coeffcompac(nscheme)
      !
      call ptds_ini(gci,alfa,nscheme,im,npdci)
      call ptds_ini(gcj,alfa,nscheme,jm,npdcj)
      if(ksize>1) call ptds_ini(gck,alfa,nscheme,km,npdck)
      !
    endif
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Calculating d<x,y,z>/d<i,j,k>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do k=0,km
    do j=0,jm
      do m=1,3
        dx(0:im,j,k,m,1)=ddfc(x(:,j,k,m),cscheme,npdci,im,alfa,gci)
      enddo
    enddo
    enddo
    !
    do k=0,km
    do i=0,im
      !
      do m=1,3
        dx(i,0:jm,k,m,2)=ddfc(x(i,:,k,m),cscheme,npdcj,jm,alfa,gcj)
      enddo
      !
    enddo
    enddo
    !
    do j=0,jm
    do i=0,im
      !
      do m=1,3
        dx(i,j,0:km,m,3)=ddfc(x(i,j,:,m),cscheme,npdck,km,alfa,gck)
      enddo
      !
    enddo
    enddo
    !
    if(ndims==3) then
      voldom=0.d0
      do k=1,km
      do j=1,jm
      do i=1,im
        celvol(i,j,k)=volhex( x(i-1,j-1,k-1,:), x(i,j-1,k-1,:),        &
                              x(i,j-1,k,:),     x(i-1,j-1,k,:),        &
                              x(i-1,j,k-1,:),   x(i,j,k-1,:)  ,        &
                              x(i,j,k,:),       x(i-1,j,k,:) )
        voldom=voldom+celvol(i,j,k)
      enddo
      enddo
      enddo
      voldom=psum(voldom)
    elseif(ndims==3) then
      voldom=0.d0
      do k=1,km
      do j=1,jm
      do i=1,im
        celvol(i,j,k)=arquad( x(i-1,j-1,k,:),   x(i,j-1,k,:),          &
                              x(i,j,k,:),       x(i-1,j,k,:) )
        voldom=voldom+celvol(i,j,k)
      enddo
      enddo
      enddo
      voldom=psum(voldom)
    endif
    if(lio) print*,' ** total volume of the domain is: ',voldom
    !
    if(lio) print*,' ** dxyz/dijk calculated.'
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! end of cal d<x,y,z>/d<i,j,k>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    call dataswap(dx)
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Calculating geometrical 
    ! Jacobian
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if(ndims==3) then
      jacob(0:im,0:jm,0:km)=dx(0:im,0:jm,0:km,1,1)*                    &
                            dx(0:im,0:jm,0:km,2,2)*                    &
                            dx(0:im,0:jm,0:km,3,3)                     &
                          + dx(0:im,0:jm,0:km,1,2)*                    &
                            dx(0:im,0:jm,0:km,2,3)*                    &
                            dx(0:im,0:jm,0:km,3,1)                     &
                          + dx(0:im,0:jm,0:km,1,3)*                    &
                            dx(0:im,0:jm,0:km,2,1)*                    &
                            dx(0:im,0:jm,0:km,3,2)                     &
                          - dx(0:im,0:jm,0:km,1,3)*                    &
                            dx(0:im,0:jm,0:km,2,2)*                    &
                            dx(0:im,0:jm,0:km,3,1)                     &
                          - dx(0:im,0:jm,0:km,1,2)*                    &
                            dx(0:im,0:jm,0:km,2,1)*                    &
                            dx(0:im,0:jm,0:km,3,3)                     &
                          - dx(0:im,0:jm,0:km,1,1)*                    &
                            dx(0:im,0:jm,0:km,2,3)*                    &
                            dx(0:im,0:jm,0:km,3,2)
    elseif(ndims==2) then  
      jacob(0:im,0:jm,0:km)=dx(0:im,0:jm,0:km,1,1)*                    &
                            dx(0:im,0:jm,0:km,2,2)                     &
                          - dx(0:im,0:jm,0:km,1,2)*                    &
                            dx(0:im,0:jm,0:km,2,1)
    else
      stop ' !! ndimes not defined at jacob calculation.'
    endif
    !
    call dataswap(jacob)
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! End of Calculating 
    ! geometrical Jacobian
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Calculating d<i,j,k>/d<x,y,z>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    dxi=0.d0
    !
    if(ndims==3) then
      !
      allocate( phi(-hm:im+hm)  )
      do k=0,km
      do j=0,jm
        !
        phi(:)=0.5d0*(dx(-hm:im+hm,j,k,2,3)*x(-hm:im+hm,j,k,3)-           &
                      dx(-hm:im+hm,j,k,3,3)*x(-hm:im+hm,j,k,2))
        dxi(0:im,j,k,2,1)=dxi(0:im,j,k,2,1)+ddfc(phi,cscheme,npdci,im,alfa,gci)
        !
        phi(:)=0.5d0*(dx(-hm:im+hm,j,k,3,3)*x(-hm:im+hm,j,k,1)-           &
                      dx(-hm:im+hm,j,k,1,3)*x(-hm:im+hm,j,k,3))
        dxi(0:im,j,k,2,2)=dxi(0:im,j,k,2,2)+ddfc(phi,cscheme,npdci,im,alfa,gci)
        !
        phi(:)=0.5d0*(dx(-hm:im+hm,j,k,1,3)*x(-hm:im+hm,j,k,2)-           &
                      dx(-hm:im+hm,j,k,2,3)*x(-hm:im+hm,j,k,1))
        dxi(0:im,j,k,2,3)=dxi(0:im,j,k,2,3)+ddfc(phi,cscheme,npdci,im,alfa,gci)
        !
        phi(:)=0.5d0*(dx(-hm:im+hm,j,k,3,2)*x(-hm:im+hm,j,k,2)-           &
                      dx(-hm:im+hm,j,k,2,2)*x(-hm:im+hm,j,k,3))
        dxi(0:im,j,k,3,1)=dxi(0:im,j,k,3,1)+ddfc(phi,cscheme,npdci,im,alfa,gci)
        !
        phi(:)=0.5d0*(dx(-hm:im+hm,j,k,1,2)*x(-hm:im+hm,j,k,3)-           &
                      dx(-hm:im+hm,j,k,3,2)*x(-hm:im+hm,j,k,1))
        dxi(0:im,j,k,3,2)=dxi(0:im,j,k,3,2)+ddfc(phi,cscheme,npdci,im,alfa,gci)
        !
        phi(:)=0.5d0*(dx(-hm:im+hm,j,k,2,2)*x(-hm:im+hm,j,k,1)-           &
                      dx(-hm:im+hm,j,k,1,2)*x(-hm:im+hm,j,k,2))
        dxi(0:im,j,k,3,3)=dxi(0:im,j,k,3,3)+ddfc(phi,cscheme,npdci,im,alfa,gci)
        !
      enddo
      enddo
      deallocate( phi )
      !
      allocate( phi(-hm:jm+hm)  )
      do k=0,km
      do i=0,im
        !
        phi(:)=0.5d0*(dx(i,-hm:jm+hm,k,3,3)*x(i,-hm:jm+hm,k,2)-           &
                      dx(i,-hm:jm+hm,k,2,3)*x(i,-hm:jm+hm,k,3))
        dxi(i,0:jm,k,1,1)=dxi(i,0:jm,k,1,1)+ddfc(phi,cscheme,npdcj,jm,alfa,gcj)
        !
        phi(:)=0.5d0*(dx(i,-hm:jm+hm,k,1,3)*x(i,-hm:jm+hm,k,3)-           &
                      dx(i,-hm:jm+hm,k,3,3)*x(i,-hm:jm+hm,k,1))
        dxi(i,0:jm,k,1,2)=dxi(i,0:jm,k,1,2)+ddfc(phi,cscheme,npdcj,jm,alfa,gcj)
        !
        phi(:)=0.5d0*(dx(i,-hm:jm+hm,k,2,3)*x(i,-hm:jm+hm,k,1)-           &
                      dx(i,-hm:jm+hm,k,1,3)*x(i,-hm:jm+hm,k,2))
        dxi(i,0:jm,k,1,3)=dxi(i,0:jm,k,1,3)+ddfc(phi,cscheme,npdcj,jm,alfa,gcj)
        !
        phi(:)=0.5d0*(dx(i,-hm:jm+hm,k,2,1)*x(i,-hm:jm+hm,k,3)-           &
                      dx(i,-hm:jm+hm,k,3,1)*x(i,-hm:jm+hm,k,2))
        dxi(i,0:jm,k,3,1)=dxi(i,0:jm,k,3,1)+ddfc(phi,cscheme,npdcj,jm,alfa,gcj)
        !
        phi(:)=0.5d0*(dx(i,-hm:jm+hm,k,3,1)*x(i,-hm:jm+hm,k,1)-           &
                      dx(i,-hm:jm+hm,k,1,1)*x(i,-hm:jm+hm,k,3))
        dxi(i,0:jm,k,3,2)=dxi(i,0:jm,k,3,2)+ddfc(phi,cscheme,npdcj,jm,alfa,gcj)
        !
        phi(:)=0.5d0*(dx(i,-hm:jm+hm,k,1,1)*x(i,-hm:jm+hm,k,2)-           &
                      dx(i,-hm:jm+hm,k,2,1)*x(i,-hm:jm+hm,k,1))
        dxi(i,0:jm,k,3,3)=dxi(i,0:jm,k,3,3)+ddfc(phi,cscheme,npdcj,jm,alfa,gcj)
      enddo
      enddo
      deallocate( phi )
      !
      allocate( phi(-hm:km+hm)  )
      do j=0,jm
      do i=0,im
        !
        phi(:)=0.5d0*(dx(i,j,-hm:km+hm,2,2)*x(i,j,-hm:km+hm,3)-           &
                      dx(i,j,-hm:km+hm,3,2)*x(i,j,-hm:km+hm,2))
        dxi(i,j,0:km,1,1)=dxi(i,j,0:km,1,1)+ddfc(phi,cscheme,npdck,km,alfa,gck)
        !
        phi(:)=0.5d0*(dx(i,j,-hm:km+hm,3,2)*x(i,j,-hm:km+hm,1)-           &
                      dx(i,j,-hm:km+hm,1,2)*x(i,j,-hm:km+hm,3))
        dxi(i,j,0:km,1,2)=dxi(i,j,0:km,1,2)+ddfc(phi,cscheme,npdck,km,alfa,gck)
        !
        phi(:)=0.5d0*(dx(i,j,-hm:km+hm,1,2)*x(i,j,-hm:km+hm,2)-           &
                      dx(i,j,-hm:km+hm,2,2)*x(i,j,-hm:km+hm,1))
        dxi(i,j,0:km,1,3)=dxi(i,j,0:km,1,3)+ddfc(phi,cscheme,npdck,km,alfa,gck)
        !
        phi(:)=0.5d0*(dx(i,j,-hm:km+hm,3,1)*x(i,j,-hm:km+hm,2)-           &
                      dx(i,j,-hm:km+hm,2,1)*x(i,j,-hm:km+hm,3))
        dxi(i,j,0:km,2,1)=dxi(i,j,0:km,2,1)+ddfc(phi,cscheme,npdck,km,alfa,gck)
        !
        phi(:)=0.5d0*(dx(i,j,-hm:km+hm,1,1)*x(i,j,-hm:km+hm,3)-           &
                      dx(i,j,-hm:km+hm,3,1)*x(i,j,-hm:km+hm,1))
        dxi(i,j,0:km,2,2)=dxi(i,j,0:km,2,2)+ddfc(phi,cscheme,npdck,km,alfa,gck)
        !
        phi(:)=0.5d0*(dx(i,j,-hm:km+hm,2,1)*x(i,j,-hm:km+hm,1)-           &
                      dx(i,j,-hm:km+hm,1,1)*x(i,j,-hm:km+hm,2))
        dxi(i,j,0:km,2,3)=dxi(i,j,0:km,2,3)+ddfc(phi,cscheme,npdck,km,alfa,gck)
        !
      enddo
      enddo
      deallocate( phi )
      !
    elseif(ndims==2) then
      !
      dxi(0:im,0:jm,0:km,1,1)= dx(0:im,0:jm,0:km,2,2)  
      dxi(0:im,0:jm,0:km,1,2)=-dx(0:im,0:jm,0:km,1,2)  
      dxi(0:im,0:jm,0:km,2,1)=-dx(0:im,0:jm,0:km,2,1)  
      dxi(0:im,0:jm,0:km,2,2)= dx(0:im,0:jm,0:km,1,1) 
      !
      dxi(0:im,0:jm,0:km,1,3)=0.d0
      dxi(0:im,0:jm,0:km,2,3)=0.d0
      dxi(0:im,0:jm,0:km,3,:)=0.d0
      !
    endif
    
    !
    call dataswap(dxi)
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Calculating geometrical
    ! metric identities
    !!!!!!!!!!!!!!!!!!!!!!!!!!!
    allocate(can(0:im,0:jm,0:km,1:3))
    !
    can=0.d0
    !
    do k=0,km
    do j=0,jm
      can(:,j,k,1)=can(:,j,k,1)+ddfc(dxi(:,j,k,1,1),cscheme,npdci,im,alfa,gci)
      can(:,j,k,2)=can(:,j,k,2)+ddfc(dxi(:,j,k,1,2),cscheme,npdci,im,alfa,gci)
      can(:,j,k,3)=can(:,j,k,3)+ddfc(dxi(:,j,k,1,3),cscheme,npdci,im,alfa,gci)
    enddo
    enddo
    !
    do k=0,km
    do i=0,im
      can(i,:,k,1)=can(i,:,k,1)+ddfc(dxi(i,:,k,2,1),cscheme,npdcj,jm,alfa,gcj)
      can(i,:,k,2)=can(i,:,k,2)+ddfc(dxi(i,:,k,2,2),cscheme,npdcj,jm,alfa,gcj)
      can(i,:,k,3)=can(i,:,k,3)+ddfc(dxi(i,:,k,2,3),cscheme,npdcj,jm,alfa,gcj)
    enddo
    enddo
    !
    do j=0,jm
    do i=0,im
      can(i,j,:,1)=can(i,j,:,1)+ddfc(dxi(i,j,:,3,1),cscheme,npdck,km,alfa,gck)
      can(i,j,:,2)=can(i,j,:,2)+ddfc(dxi(i,j,:,3,2),cscheme,npdck,km,alfa,gck)
      can(i,j,:,3)=can(i,j,:,3)+ddfc(dxi(i,j,:,3,3),cscheme,npdck,km,alfa,gck)
    enddo
    enddo
    !
    can1av=0.d0
    can2av=0.d0
    can3av=0.d0
    if(ndims==3) then
      do k=1,km
      do j=1,jm
      do i=1,im
        can1av=can1av+can(i,j,k,1)
        can2av=can2av+can(i,j,k,2)
        can3av=can3av+can(i,j,k,3)
        !
        if(isnan(can(i,j,k,1))) print*,mpirank,'|-1',i,j,k,can(i,j,k,:)
        if(isnan(can(i,j,k,2))) print*,mpirank,'|-2',i,j,k,can(i,j,k,:)
        if(isnan(can(i,j,k,3))) print*,mpirank,'|-3',i,j,k,can(i,j,k,:)
        !
      end do
      end do
      end do
      can1av=psum(can1av)/real(ia*ja*ka,8)
      can2av=psum(can2av)/real(ia*ja*ka,8)
      can3av=psum(can3av)/real(ia*ja*ka,8)
    elseif(ndims==2) then
      do j=1,jm
      do i=1,im
        can1av=can1av+can(i,j,k,1)
        can2av=can2av+can(i,j,k,2)
        can3av=can3av+can(i,j,k,3)
        !
        if(isnan(can(i,j,k,1))) print*,mpirank,'|-1',i,j,k,can(i,j,k,:)
        if(isnan(can(i,j,k,2))) print*,mpirank,'|-2',i,j,k,can(i,j,k,:)
        if(isnan(can(i,j,k,3))) print*,mpirank,'|-3',i,j,k,can(i,j,k,:)
        !
      end do
      end do
      can1av=psum(can1av)/real(ia*ja,8)
      can2av=psum(can2av)/real(ia*ja,8)
      can3av=psum(can3av)/real(ia*ja,8)
    endif
    !
    !
    can1var=0.d0
    can2var=0.d0
    can3var=0.d0
    !
    if(ndims==3) then
      do k=1,km
      do j=1,jm
      do i=1,im
        can1var=can1var+(can(i,j,k,1)-can1av)**2
        can2var=can2var+(can(i,j,k,2)-can2av)**2
        can3var=can3var+(can(i,j,k,3)-can3av)**2
      end do
      end do
      end do
      can1var=psum(can1var)/real(ia*ja*ka,8)
      can2var=psum(can2var)/real(ia*ja*ka,8)
      can3var=psum(can3var)/real(ia*ja*ka,8)
    elseif(ndims==2) then
      do j=1,jm
      do i=1,im
        can1var=can1var+(can(i,j,k,1)-can1av)**2
        can2var=can2var+(can(i,j,k,2)-can2av)**2
        can3var=can3var+(can(i,j,k,3)-can3av)**2
      end do
      end do
      can1var=psum(can1var)/real(ia*ja,8)
      can2var=psum(can2var)/real(ia*ja,8)
      can3var=psum(can3var)/real(ia*ja,8)
    endif
    !
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! End ofCalculating 
    ! geometrical metric 
    ! identities
    !!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! normlizing 
    ! d<i,j,k>/d<x,y,z> and
    ! geometrical vector.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    do j=1,3
    do i=1,3
      dxi(0:im,0:jm,0:km,i,j)=dxi(0:im,0:jm,0:km,i,j)/jacob(0:im,0:jm,0:km)
    enddo
    enddo
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! End of normlizing
    ! d<i,j,k>/d<x,y,z> and
    ! geometrical vector.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    call dataswap(dxi)
    !
    xmax=maxval(x(0:im,0:jm,0:km,1))
    ymax=maxval(x(0:im,0:jm,0:km,2))
    zmax=maxval(x(0:im,0:jm,0:km,3))
    !
    xmin=minval(x(0:im,0:jm,0:km,1))
    ymin=minval(x(0:im,0:jm,0:km,2))
    zmin=minval(x(0:im,0:jm,0:km,3))
    !
    xmax=pmax(xmax)
    ymax=pmax(ymax)
    zmax=pmax(zmax)
    !
    xmin=pmin(xmin)
    ymin=pmin(ymin)
    zmin=pmin(zmin)
    !
    !
    if(lio) then
      !
      write(*,'(2X,62A)')('-',i=1,62)
      write(*,'(2X,A)')'                    *** Grids Information *** '
      write(*,'(2X,62A)')('-',i=1,62)
      write(*,'(3X,62A)')'       xmin      xmax      ymin      ymax      zmin      zmax'
      write(*,"(4X,6(F10.6))")xmin,xmax,ymin,ymax,zmin,zmax
      write(*,'(2X,62A)')('-',i=1,62)
      write(*,'(2X,A)')'                   *** Averaged of Identity ***'
      write(*,"(1X,3(1X,E20.7E3))")can1av,can2av,can3av
      write(*,'(2X,62A)')('-',i=1,62)
      write(*,'(2X,A)')'                   *** Variance of Identity ***'
      write(*,"(1X,3(1X,E20.7E3))")can1var,can2var,can3var
      write(*,'(2X,62A)')('-',i=1,62)
      !
      if(can1av>1d-15 .or. can2av>1d-15 .or. can3av>1d-15) then
        write(*,*)' !! Warning: Averaged Grids Identity is too large'
        write(*,'(2X,62A)')('-',i=1,62)
      end if
      if(can1var>1d-15 .or. can2var>1d-15 .or. can3var>1d-15) then
        write(*,*)' !! Warning: Variance of Grids Identity is too large'
        write(*,'(2X,62A)')('-',i=1,62)
      end if
      !
      print*,' ** geometrical parameters calculated'
    end if
    !
    deallocate(dx,gci,gcj,can)
    !
    if(allocated(gck)) deallocate(gck)
    !
    ! if(mpirank==0) then
    !   do k=0,km
    !     print*,x(1,1,k,3),dx(1,1,k,3,3),x(1,1,k,3)-x(1,1,k-1,3)
    !   enddo
    ! endif
    ! call tecbin('testout/tecgrid'//mpirankname//'.plt',                &
    !                                   x(0:im,0:jm,0:km,1),'x',    &
    !                                   x(0:im,0:jm,0:km,2),'y',    &
    !                                   x(0:im,0:jm,0:km,3),'z',    &
    !                                 dxi(0:im,0:jm,0:km,1,1),'dxdi',    &
    !                                 dxi(0:im,0:jm,0:km,2,2),'dydj')
    ! ! 
    ! call mpistop
    !
  end subroutine geomcal
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! End of the subroutine GeomCal.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to calulate the rhs of N-S equations.          |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 09-02-2021  | Created by J. Fang @ Warrington                     |
  !+-------------------------------------------------------------------+
  subroutine rhscal
    !
    use commarray, only : qrhs
    !
    qrhs=0.d0
    !
    call convrsdcal6
    !
    qrhs=-qrhs
    !
    call diffrsdcal6
    !
  end subroutine rhscal
  !+-------------------------------------------------------------------+
  !| The end of the subroutine rhscal.                                 |
  !+-------------------------------------------------------------------+
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! This subroutine is used to calculate the convectional residual
  ! terms with compact six-order central scheme.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Writen by Fang Jian, 2009-06-03.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine convrsdcal6
    !
    use commvar,  only: im,jm,km,hm,numq,num_species,conschm,          &
                        npdci,npdcj,npdck
    use commarray,only: q,vel,rho,prs,tmp,spc,dxi,jacob,qrhs
    use commfunc, only: ddfc
    !
    ! local data
    integer :: i,j,k,jspc,m
    real(8),allocatable :: fcs(:,:),dfcs(:,:),uu(:) 
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! calculating along i direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    allocate(fcs(-hm:im+hm,1:numq),dfcs(0:im,1:numq),uu(-hm:im+hm))
    do k=0,km
    do j=0,jm
      !
      uu(:)=dxi(:,j,k,1,1)*vel(:,j,k,1)+dxi(:,j,k,1,2)*vel(:,j,k,2) +  &
            dxi(:,j,k,1,3)*vel(:,j,k,3)
      fcs(:,1)=jacob(:,j,k)*  q(:,j,k,1)*uu
      fcs(:,2)=jacob(:,j,k)*( q(:,j,k,2)*uu+dxi(:,j,k,1,1)*prs(:,j,k) )
      fcs(:,3)=jacob(:,j,k)*( q(:,j,k,3)*uu+dxi(:,j,k,1,2)*prs(:,j,k) )
      fcs(:,4)=jacob(:,j,k)*( q(:,j,k,4)*uu+dxi(:,j,k,1,3)*prs(:,j,k) )
      fcs(:,5)=jacob(:,j,k)*( q(:,j,k,5)+prs(:,j,k) )*uu
      !
      if(num_species>0) then
        do jspc=1,num_species
          fcs(:,5+jspc)=jacob(:,j,k)*q(:,j,k,5+jspc)*uu
        enddo
      endif
      !
      do m=1,numq
        dfcs(:,m)=ddfc(fcs(:,m),conschm,npdci,im,alfa_con,cci)
        !
        qrhs(:,j,k,m)=qrhs(:,j,k,m)+dfcs(:,m)
      enddo
      !
      ! print*,maxval(dfcs(:,2)),maxval(dfcs(:,3))
      !
      !
    enddo
    enddo
    deallocate(fcs,dfcs,uu)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! end calculating along i direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! calculating along j direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    allocate(fcs(-hm:jm+hm,1:numq),dfcs(0:jm,1:numq),uu(-hm:jm+hm))
    do k=0,km
    do i=0,im
      !
      uu(:)=dxi(i,:,k,2,1)*vel(i,:,k,1)+dxi(i,:,k,2,2)*vel(i,:,k,2) +  &
            dxi(i,:,k,2,3)*vel(i,:,k,3)
      fcs(:,1)=jacob(i,:,k)*  q(i,:,k,1)*uu
      fcs(:,2)=jacob(i,:,k)*( q(i,:,k,2)*uu+dxi(i,:,k,2,1)*prs(i,:,k) )
      fcs(:,3)=jacob(i,:,k)*( q(i,:,k,3)*uu+dxi(i,:,k,2,2)*prs(i,:,k) )
      fcs(:,4)=jacob(i,:,k)*( q(i,:,k,4)*uu+dxi(i,:,k,2,3)*prs(i,:,k) )
      fcs(:,5)=jacob(i,:,k)*( q(i,:,k,5)+prs(i,:,k) )*uu
      !
      if(num_species>0) then
        do jspc=1,num_species
          fcs(:,5+jspc)=jacob(i,:,k)*q(i,:,k,5+jspc)*uu
        enddo
      endif
      !
      do m=1,numq
        dfcs(:,m)=ddfc(fcs(:,m),conschm,npdcj,jm,alfa_con,ccj)
        !
        qrhs(i,:,k,m)=qrhs(i,:,k,m)+dfcs(:,m)
      enddo
      !
    enddo
    enddo
    deallocate(fcs,dfcs,uu)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! end calculating along j direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    if(ndims==3) then
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! calculating along j direction
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      allocate(fcs(-hm:km+hm,1:numq),dfcs(0:km,1:numq),uu(-hm:km+hm))
      do j=0,jm
      do i=0,im
        !
        uu(:)=dxi(i,j,:,3,1)*vel(i,j,:,1)+dxi(i,j,:,3,2)*vel(i,j,:,2) +  &
              dxi(i,j,:,3,3)*vel(i,j,:,3)
        fcs(:,1)=jacob(i,j,:)*  q(i,j,:,1)*uu
        fcs(:,2)=jacob(i,j,:)*( q(i,j,:,2)*uu+dxi(i,j,:,3,1)*prs(i,j,:) )
        fcs(:,3)=jacob(i,j,:)*( q(i,j,:,3)*uu+dxi(i,j,:,3,2)*prs(i,j,:) )
        fcs(:,4)=jacob(i,j,:)*( q(i,j,:,4)*uu+dxi(i,j,:,3,3)*prs(i,j,:) )
        fcs(:,5)=jacob(i,j,:)*( q(i,j,:,5)+prs(i,j,:) )*uu
        !
        if(num_species>0) then
          do jspc=1,num_species
            fcs(:,5+jspc)=jacob(i,j,:)*q(i,j,:,5+jspc)*uu
          enddo
        endif
        !
        do m=1,numq
          dfcs(:,m)=ddfc(fcs(:,m),conschm,npdck,km,alfa_con,cck)
          !
          qrhs(i,j,:,m)=qrhs(i,j,:,m)+dfcs(:,m)
        enddo
        !
      enddo
      enddo
      deallocate(fcs,dfcs,uu)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! end calculating along j direction
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
    endif
    !
  end subroutine convrsdcal6
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! End of the subroutine ConvRsdCal6.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! This subroutine is used to calculate the diffusion term with 6-order
  ! Compact Central scheme.
  !   sixth-order Compact Central scheme.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Writen by Fang Jian, 2009-06-09.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine diffrsdcal6
    !
    use commvar,   only : im,jm,km,numq,npdci,npdcj,npdck,difschm,     &
                          conschm,ndims,num_species,reynolds,prandtl,  &
                          const5
    use commarray, only : vel,tmp,spc,dvel,dtmp,dspc,dxi,x,jacob,qrhs
    use commfunc,  only : ddfc
    use fludyna,   only : miucal
    use tecio
    !
    ! local data
    integer :: i,j,k,n
    real(8),allocatable :: df(:),ff(:)
    real(8),allocatable,dimension(:,:,:,:) :: sigma,qflux
    real(8) :: miu,miu2,hcc,s11,s12,s13,s22,s23,s33,skk
    !
    dvel=0.d0
    dtmp=0.d0
    dspc=0.d0
    !
    ! calculate velocity and temperature gradient
    !
    allocate(df(0:im))
    !
    do k=0,km
    do j=0,jm
      !
      do n=1,3
        df=ddfc(vel(:,j,k,n),difschm,npdci,im,alfa_dif,dci)
        dvel(:,j,k,n,1)=dvel(:,j,k,n,1)+df*dxi(0:im,j,k,1,1)
        dvel(:,j,k,n,2)=dvel(:,j,k,n,2)+df*dxi(0:im,j,k,1,2)
        dvel(:,j,k,n,3)=dvel(:,j,k,n,3)+df*dxi(0:im,j,k,1,3)
        !
      enddo
      !
      df=ddfc(tmp(:,j,k),difschm,npdci,im,alfa_dif,dci)
      dtmp(:,j,k,1)=dtmp(:,j,k,1)+df*dxi(0:im,j,k,1,1)
      dtmp(:,j,k,2)=dtmp(:,j,k,2)+df*dxi(0:im,j,k,1,2)
      dtmp(:,j,k,3)=dtmp(:,j,k,3)+df*dxi(0:im,j,k,1,3)
      !
      if(num_species>1) then
        do n=1,num_species
          df=ddfc(spc(:,j,k,n),difschm,npdci,im,alfa_dif,dci)
          dspc(:,j,k,n,1)=dspc(:,j,k,n,1)+df*dxi(0:im,j,k,1,1)
          dspc(:,j,k,n,2)=dspc(:,j,k,n,2)+df*dxi(0:im,j,k,1,2)
          dspc(:,j,k,n,3)=dspc(:,j,k,n,3)+df*dxi(0:im,j,k,1,3)
        enddo
      endif
      !
    enddo
    enddo
    !
    deallocate(df)
    !
    allocate(df(0:jm))
    do k=0,km
    do i=0,im
      !
      do n=1,3
        df=ddfc(vel(i,:,k,n),difschm,npdcj,jm,alfa_dif,dcj)
        dvel(i,:,k,n,1)=dvel(i,:,k,n,1)+df*dxi(i,0:jm,k,2,1)
        dvel(i,:,k,n,2)=dvel(i,:,k,n,2)+df*dxi(i,0:jm,k,2,2)
        dvel(i,:,k,n,3)=dvel(i,:,k,n,3)+df*dxi(i,0:jm,k,2,3)
      enddo
      !
      df=ddfc(tmp(i,:,k),difschm,npdcj,jm,alfa_dif,dcj)
      dtmp(i,:,k,1)=dtmp(i,:,k,1)+df*dxi(i,0:jm,k,2,1)
      dtmp(i,:,k,2)=dtmp(i,:,k,2)+df*dxi(i,0:jm,k,2,2)
      dtmp(i,:,k,3)=dtmp(i,:,k,3)+df*dxi(i,0:jm,k,2,3)
      !
      if(num_species>1) then
        do n=1,num_species
          df=ddfc(spc(i,:,k,n),difschm,npdcj,jm,alfa_dif,dcj)
          dspc(i,:,k,n,1)=dspc(i,:,k,n,1)+df*dxi(i,0:jm,k,2,1)
          dspc(i,:,k,n,2)=dspc(i,:,k,n,2)+df*dxi(i,0:jm,k,2,2)
          dspc(i,:,k,n,3)=dspc(i,:,k,n,3)+df*dxi(i,0:jm,k,2,3)
        enddo
      endif
      !
    enddo
    enddo
    deallocate(df)
    !
    if(ndims==3) then
      allocate(df(0:km))
      do j=0,jm
      do i=0,im
        !
        do n=1,3
          df=ddfc(vel(i,j,:,n),difschm,npdck,km,alfa_dif,dck)
          dvel(i,j,:,n,1)=dvel(i,j,:,n,1)+df*dxi(i,j,0:km,3,1)
          dvel(i,j,:,n,2)=dvel(i,j,:,n,2)+df*dxi(i,j,0:km,3,2)
          dvel(i,j,:,n,3)=dvel(i,j,:,n,3)+df*dxi(i,j,0:km,3,3)
        enddo
        !
        df=ddfc(tmp(i,j,:),difschm,npdck,km,alfa_dif,dck)
        dtmp(i,j,:,1)=dtmp(i,j,:,1)+df*dxi(i,j,0:km,3,1)
        dtmp(i,j,:,2)=dtmp(i,j,:,2)+df*dxi(i,j,0:km,3,2)
        dtmp(i,j,:,3)=dtmp(i,j,:,3)+df*dxi(i,j,0:km,3,3)
        !
        if(num_species>1) then
          do n=1,num_species
            df=ddfc(spc(i,j,:,n),difschm,npdck,km,alfa_dif,dck)
            dspc(i,j,:,n,1)=dspc(i,j,:,n,1)+df*dxi(i,j,0:km,3,1)
            dspc(i,j,:,n,2)=dspc(i,j,:,n,2)+df*dxi(i,j,0:km,3,2)
            dspc(i,j,:,n,3)=dspc(i,j,:,n,3)+df*dxi(i,j,0:km,3,3)
          enddo
        endif
        !
      enddo
      enddo
      deallocate(df)
    endif
    !
    allocate( sigma(-hm:im+hm,-hm:jm+hm,-hm:km+hm,1:6),                &
              qflux(-hm:im+hm,-hm:jm+hm,-hm:km+hm,1:3) )
    ! 
    do k=0,km
    do j=0,jm
    do i=0,im
      !
      miu=miucal(tmp(i,j,k))
      !
      s11=dvel(i,j,k,1,1)
      s12=0.5d0*(dvel(i,j,k,1,2)+dvel(i,j,k,2,1))
      s13=0.5d0*(dvel(i,j,k,1,3)+dvel(i,j,k,3,1))
      s22=dvel(i,j,k,2,2)
      s23=0.5d0*(dvel(i,j,k,2,3)+dvel(i,j,k,3,2))
      s33=dvel(i,j,k,3,3)
      !
      skk=num1d3*(s11+s22+s33)
      !
      miu2=2.d0*miu/reynolds
      hcc=(miu/reynolds/prandtl)/const5
      !
      sigma(i,j,k,1)=miu2*(s11-skk)  !s11   
      sigma(i,j,k,2)=miu2*s12        !s12  
      sigma(i,j,k,3)=miu2*s13        !s13   
      sigma(i,j,k,4)=miu2*(s22-skk)  !s22   
      sigma(i,j,k,5)=miu2*s23        !s23  
      sigma(i,j,k,6)=miu2*(s33-skk)  !s33  
      !
      qflux(i,j,k,1)=hcc*dtmp(i,j,k,1)+sigma(i,j,k,1)*vel(i,j,k,1) +   &
                                       sigma(i,j,k,2)*vel(i,j,k,2) +   &
                                       sigma(i,j,k,3)*vel(i,j,k,3)
      qflux(i,j,k,2)=hcc*dtmp(i,j,k,2)+sigma(i,j,k,2)*vel(i,j,k,1) +   &
                                       sigma(i,j,k,4)*vel(i,j,k,2) +   &
                                       sigma(i,j,k,5)*vel(i,j,k,3)
      qflux(i,j,k,3)=hcc*dtmp(i,j,k,3)+sigma(i,j,k,3)*vel(i,j,k,1) +   &
                                       sigma(i,j,k,5)*vel(i,j,k,2) +   &
                                       sigma(i,j,k,6)*vel(i,j,k,3)
    enddo
    enddo
    enddo
    !
    call dataswap(sigma)
    !
    call dataswap(qflux)
    !
    ! Calculating along i direction.
    !
    allocate(ff(-hm:im+hm),df(0:im))
    do k=0,km
    do j=0,jm
      !
      ff(:)=( sigma(:,j,k,1)*dxi(:,j,k,1,1) +                          &
              sigma(:,j,k,2)*dxi(:,j,k,1,2) +                          &
              sigma(:,j,k,3)*dxi(:,j,k,1,3) )*jacob(:,j,k)
      !
      df=ddfc(ff,difschm,npdci,im,alfa_dif,dci)
      !
      qrhs(:,j,k,2)=qrhs(:,j,k,2)+df
      !
      ff(:)=( sigma(:,j,k,2)*dxi(:,j,k,1,1) +                          &
              sigma(:,j,k,4)*dxi(:,j,k,1,2) +                          &
              sigma(:,j,k,5)*dxi(:,j,k,1,3) )*jacob(:,j,k)
      !
      df=ddfc(ff,difschm,npdci,im,alfa_dif,dci)
      !
      qrhs(:,j,k,3)=qrhs(:,j,k,3)+df
      !
      ff(:)=( sigma(:,j,k,3)*dxi(:,j,k,1,1) +                          &
              sigma(:,j,k,5)*dxi(:,j,k,1,2) +                          &
              sigma(:,j,k,6)*dxi(:,j,k,1,3) )*jacob(:,j,k)
      !
      df=ddfc(ff,difschm,npdci,im,alfa_dif,dci)
      !
      qrhs(:,j,k,4)=qrhs(:,j,k,4)+df
      !
      ff(:)=( qflux(:,j,k,1)*dxi(:,j,k,1,1) +                          &
              qflux(:,j,k,2)*dxi(:,j,k,1,2) +                          &
              qflux(:,j,k,3)*dxi(:,j,k,1,3) )*jacob(:,j,k)
      !
      df=ddfc(ff,difschm,npdci,im,alfa_dif,dci)
      !
      qrhs(:,j,k,5)=qrhs(:,j,k,5)+df
      !
      !
    enddo
    enddo
    !
    deallocate(ff,df)
    !!!!!!!!!!!!!!!!!!!!!!!!!!
    ! End calculating along i
    !!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! Calculating along j direction.
    !
    allocate(ff(-hm:jm+hm),df(0:jm))
    do k=0,km
    do i=0,im
      !
      ff(:)=( sigma(i,:,k,1)*dxi(i,:,k,2,1) +                          &
              sigma(i,:,k,2)*dxi(i,:,k,2,2) +                          &
              sigma(i,:,k,3)*dxi(i,:,k,2,3) )*jacob(i,:,k)
      !
      df=ddfc(ff,difschm,npdcj,jm,alfa_dif,dcj)
      !
      qrhs(i,:,k,2)=qrhs(i,:,k,2)+df
      !
      ff(:)=( sigma(i,:,k,2)*dxi(i,:,k,2,1) +                          &
              sigma(i,:,k,4)*dxi(i,:,k,2,2) +                          &
              sigma(i,:,k,5)*dxi(i,:,k,2,3) )*jacob(i,:,k)
      !
      df=ddfc(ff,difschm,npdcj,jm,alfa_dif,dcj)
      !
      qrhs(i,:,k,3)=qrhs(i,:,k,3)+df
      !
      ff(:)=( sigma(i,:,k,3)*dxi(i,:,k,2,1) +                          &
              sigma(i,:,k,5)*dxi(i,:,k,2,2) +                          &
              sigma(i,:,k,6)*dxi(i,:,k,2,3) )*jacob(i,:,k)
      !
      df=ddfc(ff,difschm,npdcj,jm,alfa_dif,dcj)
      !
      qrhs(i,:,k,4)=qrhs(i,:,k,4)+df
      !
      ff(:)=( qflux(i,:,k,1)*dxi(i,:,k,2,1) +                          &
              qflux(i,:,k,2)*dxi(i,:,k,2,2) +                          &
              qflux(i,:,k,3)*dxi(i,:,k,2,3) )*jacob(i,:,k)
      !
      df=ddfc(ff,difschm,npdcj,jm,alfa_dif,dcj)
      !
      qrhs(i,:,k,5)=qrhs(i,:,k,5)+df
      !
      !
    enddo
    enddo
    !
    deallocate(ff,df)
    !!!!!!!!!!!!!!!!!!!!!!!!!!
    ! End calculating along j
    !!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    if(ndims==3) then
      ! Calculating along j direction.
      !
      allocate(ff(-hm:km+hm),df(0:km))
      do j=0,jm
      do i=0,im
        !
        ff(:)=( sigma(i,j,:,1)*dxi(i,j,:,3,1) +                        &
                sigma(i,j,:,2)*dxi(i,j,:,3,2) +                        &
                sigma(i,j,:,3)*dxi(i,j,:,3,3) )*jacob(i,j,:)
        !
        df=ddfc(ff,difschm,npdck,km,alfa_dif,dck)
        !
        qrhs(i,j,:,2)=qrhs(i,j,:,2)+df
        !
        ff(:)=( sigma(i,j,:,2)*dxi(i,j,:,3,1) +                        &
                sigma(i,j,:,4)*dxi(i,j,:,3,2) +                        &
                sigma(i,j,:,5)*dxi(i,j,:,3,3) )*jacob(i,j,:)
        !
        df=ddfc(ff,difschm,npdck,km,alfa_dif,dck)
        !
        qrhs(i,j,:,3)=qrhs(i,j,:,3)+df
        !
        ff(:)=( sigma(i,j,:,3)*dxi(i,j,:,3,1) +                        &
                sigma(i,j,:,5)*dxi(i,j,:,3,2) +                        &
                sigma(i,j,:,6)*dxi(i,j,:,3,3) )*jacob(i,j,:)
        !
        df=ddfc(ff,difschm,npdck,km,alfa_dif,dck)
        !
        qrhs(i,j,:,4)=qrhs(i,j,:,4)+df
        !
        ff(:)=( qflux(i,j,:,1)*dxi(i,j,:,3,1) +                        &
                qflux(i,j,:,2)*dxi(i,j,:,3,2) +                        &
                qflux(i,j,:,3)*dxi(i,j,:,3,3) )*jacob(i,j,:)
        !
        df=ddfc(ff,difschm,npdck,km,alfa_dif,dck)
        !
        qrhs(i,j,:,5)=qrhs(i,j,:,5)+df
        !
        !
      enddo
      enddo
      !
      deallocate(ff,df)
      !!!!!!!!!!!!!!!!!!!!!!!!!!
      ! End calculating along j
      !!!!!!!!!!!!!!!!!!!!!!!!!!
    endif
    !
    ! call tecbin('testout/tecgrad'//mpirankname//'.plt',                &
    !                                   x(0:im,0:jm,0:km,1),'x',         &
    !                                   x(0:im,0:jm,0:km,2),'y',         &
    !                                 vel(0:im,0:jm,0:km,1),'u',         &
    !                                 vel(0:im,0:jm,0:km,2),'v',         &
    !                                dvel(0:im,0:jm,0:km,1,1),'dudx',    &
    !                                dvel(0:im,0:jm,0:km,1,2),'dudy',    &
    !                                dvel(0:im,0:jm,0:km,2,1),'dvdx',    &
    !                                dvel(0:im,0:jm,0:km,2,2),'dvdy' )
    !
    deallocate(sigma,qflux)
    !
  end subroutine diffrsdcal6
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! End of the subroutine diffrsdcal6.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
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
  subroutine filterq
    !
    use commvar,  only : im,jm,km,numq,npdci,npdcj,npdck,alfa_filter,ndims
    use commarray,only : q
    use commfunc, only : spafilter10
    !
    ! local data
    integer :: i,j,k,m
    !
    ! filtering in i direction
    call dataswap(q,direction=1)
    !
    do k=0,km
    do j=0,jm
      !
      do m=1,numq
        q(0:im,j,k,m)=spafilter10(q(:,j,k,m),npdci,im,alfa_filter,fci)
      enddo
      !
    end do
    end do
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! end filter in i direction.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! filtering in j direction
    call dataswap(q,direction=2)
    !
    do k=0,km
    do i=0,im
      !
      do m=1,numq
        q(i,0:jm,k,m)=spafilter10(q(i,:,k,m),npdcj,jm,alfa_filter,fcj)
      enddo
      !
    end do
    end do
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! end filter in j direction.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    if(ndims==3) then
      !
      call dataswap(q,direction=3)
      !
      ! filtering in k direction
      do j=0,jm
      do i=0,im
        !
        do m=1,numq
          q(i,j,0:km,m)=spafilter10(q(i,j,:,m),npdck,km,alfa_filter,fck)
        enddo
        !
      end do
      end do
      !
    end if
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! end filter in k direction.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
  end subroutine filterq
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! End of the subroutine filterq.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!

  !!
  subroutine gradcal
    !
    use commvar,   only : im,jm,km,npdci,npdcj,npdck,conschm,alfa_filter
    use commarray, only : x,q,dxi
    use commfunc,  only : ddfc,spafilter10
    !
    ! local data
    integer :: i,j,k
    real(8),allocatable :: dq(:)
    !
    ! do k=0,km
    ! do j=0,jm
    ! do i=0,im
    !   q(i,j,k,1)=sin(0.1d0*pi*x(i,j,k,2))+0.1d0*sin(pi*j+0.5d0*pi)
    !   !
    ! enddo
    ! enddo
    ! enddo
    !
    call dataswap(q)
    !
    allocate(dq(0:im))
    !
    ! dq=ddfc(q(:,1,1,1),conschm,npdci,im,alfa_con,cci)*dxi(:,1,1,1,1)
    
    dq=spafilter10(q(:,jm,0,2),npdci,im,alfa_filter,fci)
    !
    ! allocate(dq(0:jm))
    ! dq=ddfc(q(1,:,1,1),conschm,npdcj,jm,alfa_con,ccj)*dxi(1,0:jm,1,2,2)
    ! dq=spafilter10(q(1,:,0,1),npdcj,jm,alfa_filter,fcj)
    !
    ! allocate(dq(0:km))
    ! dq=ddfc(q(1,1,:,1),conschm,npdck,km,alfa_con,cck)*dxi(1,1,0:km,3,3)
    !
    open(18,file='testout/profile'//mpirankname//'.dat')
    do i=0,im
      write(18,*)x(i,jm,0,1),q(i,jm,0,2),dq(i)
    enddo
    close(18)
    print*,' << testout/profile',mpirankname,'.dat'
    !
  end subroutine gradcal
  !
end module solver
!+---------------------------------------------------------------------+
!| The end of the module commarray.                                    |
!+---------------------------------------------------------------------+