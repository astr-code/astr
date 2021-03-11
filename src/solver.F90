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
  use parallel, only : mpirankname,mpistop,mpirank,lio,dataswap,ptime
  use commvar,  only : ndims,ks,ke,hm,ctime,hm,lfftk
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
                          xmax,xmin,ymax,ymin,zmax,zmin,voldom,difschm
    use commarray, only : x,jacob,dxi,celvol
    use parallel,  only : gridsendrecv,jsize,ksize,psum,pmax,pmin
    use commfunc,  only : coeffcompac,ptds_ini,ddfc,volhex,arquad
    use tecio
    use bc,       only : geombc,xyzbc
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
    call xyzbc
    !
    ! cscheme='442e'
    cscheme=difschm
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
      call ptds_ini(gck,alfa,nscheme,km,npdck)
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
      if(lfftk) then
        !
        do m=1,2
          dx(i,j,0:km,m,3)=ddfc(x(i,j,:,m),cscheme,npdck,km,alfa,gck,lfft=lfftk)
        enddo
        dx(i,j,0:km,3,3)=x(i,j,1,3)-x(i,j,0,3)
      else
        do m=1,3
          dx(i,j,0:km,m,3)=ddfc(x(i,j,:,m),cscheme,npdck,km,alfa,gck)
        enddo
      endif
      !
    enddo
    enddo
    !
    if(ndims==1) then
      !
      j=0
      k=0
      !
      voldom=0.d0
      do i=1,im
        celvol(i,j,k)=abs(x(i,j,k,1)-x(i-1,j,k,1))
        voldom=voldom+celvol(i,j,k)
      enddo
      voldom=psum(voldom)
      !
    elseif(ndims==2) then
      k=0
      !
      voldom=0.d0
      do j=1,jm
      do i=1,im
        celvol(i,j,k)=arquad( x(i-1,j-1,k,:),   x(i,j-1,k,:),          &
                              x(i,j,k,:),       x(i-1,j,k,:) )
        voldom=voldom+celvol(i,j,k)
      enddo
      enddo
      voldom=psum(voldom)
    elseif(ndims==3) then
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
    endif
    !
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
    elseif(ndims==1) then  
      jacob(0:im,0:jm,0:km)=dx(0:im,0:jm,0:km,1,1)*                    &
                            dx(0:im,0:jm,0:km,2,2)
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
        if(lfftk) then
          dxi(i,j,0:km,1,1)=dxi(i,j,0:km,1,1)+phi(1)-phi(0)
        else
          dxi(i,j,0:km,1,1)=dxi(i,j,0:km,1,1)+ddfc(phi,cscheme,npdck,km,alfa,gck,lfft=lfftk)
        endif
        !
        phi(:)=0.5d0*(dx(i,j,-hm:km+hm,3,2)*x(i,j,-hm:km+hm,1)-           &
                      dx(i,j,-hm:km+hm,1,2)*x(i,j,-hm:km+hm,3))
        if(lfftk) then
          dxi(i,j,0:km,1,2)=dxi(i,j,0:km,1,2)+phi(1)-phi(0)
        else
          dxi(i,j,0:km,1,2)=dxi(i,j,0:km,1,2)+ddfc(phi,cscheme,npdck,km,alfa,gck,lfft=lfftk)
        endif
        !
        phi(:)=0.5d0*(dx(i,j,-hm:km+hm,1,2)*x(i,j,-hm:km+hm,2)-           &
                      dx(i,j,-hm:km+hm,2,2)*x(i,j,-hm:km+hm,1))
        if(lfftk) then
          dxi(i,j,0:km,1,3)=dxi(i,j,0:km,1,3)+phi(1)-phi(0)
        else
          dxi(i,j,0:km,1,3)=dxi(i,j,0:km,1,3)+ddfc(phi,cscheme,npdck,km,alfa,gck,lfft=lfftk)
        endif
        !
        phi(:)=0.5d0*(dx(i,j,-hm:km+hm,3,1)*x(i,j,-hm:km+hm,2)-           &
                      dx(i,j,-hm:km+hm,2,1)*x(i,j,-hm:km+hm,3))
        if(lfftk) then
          dxi(i,j,0:km,2,1)=dxi(i,j,0:km,2,1)+phi(1)-phi(0)
        else
          dxi(i,j,0:km,2,1)=dxi(i,j,0:km,2,1)+ddfc(phi,cscheme,npdck,km,alfa,gck,lfft=lfftk)
        endif
        !
        phi(:)=0.5d0*(dx(i,j,-hm:km+hm,1,1)*x(i,j,-hm:km+hm,3)-           &
                      dx(i,j,-hm:km+hm,3,1)*x(i,j,-hm:km+hm,1))
        if(lfftk) then
          dxi(i,j,0:km,2,2)=dxi(i,j,0:km,2,2)+phi(1)-phi(0)
        else
          dxi(i,j,0:km,2,2)=dxi(i,j,0:km,2,2)+ddfc(phi,cscheme,npdck,km,alfa,gck,lfft=lfftk)
        endif
        !
        phi(:)=0.5d0*(dx(i,j,-hm:km+hm,2,1)*x(i,j,-hm:km+hm,1)-           &
                      dx(i,j,-hm:km+hm,1,1)*x(i,j,-hm:km+hm,2))
        if(lfftk) then
          dxi(i,j,0:km,2,3)=dxi(i,j,0:km,2,3)+phi(1)-phi(0)
        else
          dxi(i,j,0:km,2,3)=dxi(i,j,0:km,2,3)+ddfc(phi,cscheme,npdck,km,alfa,gck,lfft=lfftk)
        endif
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
    elseif(ndims==1) then
      !
      dxi(0:im,0:jm,0:km,1,1)= dx(0:im,0:jm,0:km,2,2)    
      dxi(0:im,0:jm,0:km,2,2)= dx(0:im,0:jm,0:km,1,1) 
      !
      dxi(0:im,0:jm,0:km,1,2)=0.d0
      dxi(0:im,0:jm,0:km,2,1)=0.d0
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
      can(i,j,:,1)=can(i,j,:,1)+ddfc(dxi(i,j,:,3,1),cscheme,npdck,km,alfa,gck,lfft=lfftk)
      can(i,j,:,2)=can(i,j,:,2)+ddfc(dxi(i,j,:,3,2),cscheme,npdck,km,alfa,gck,lfft=lfftk)
      can(i,j,:,3)=can(i,j,:,3)+ddfc(dxi(i,j,:,3,3),cscheme,npdck,km,alfa,gck,lfft=lfftk)
      !
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
      k=0
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
    call geombc
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
      write(*,"(4X,6(F10.3))")xmin,xmax,ymin,ymax,zmin,zmax
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
    ! if(mpirank==0) then
    !   do k=0,km
    !     print*,x(1,1,k,3),dx(1,1,k,3,3),x(1,1,k,3)-x(1,1,k-1,3)
    !   enddo
    ! ! endif
    ! call tecbin('testout/tecgrid'//mpirankname//'.plt',                &
    !                                   x(0:im,0:jm,0:km,1),'x',    &
    !                                   x(0:im,0:jm,0:km,2),'y',    &
    !                                   x(0:im,0:jm,0:km,3),'z',    &
    !                                 dxi(0:im,0:jm,0:km,1,1),'dxdi',    &
    !                                 dxi(0:im,0:jm,0:km,2,2),'dydj',    &
    !                                 dxi(0:im,0:jm,0:km,3,3),'dzdk')
    ! ! ! 
    deallocate(dx,can)
    !
    if(allocated(gci)) deallocate(gci)
    if(allocated(gcj)) deallocate(gcj)
    if(allocated(gck)) deallocate(gck)
    !
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
    use commarray, only : qrhs,x,q
    use commvar,   only : flowtype
    !
    integer :: j
#ifdef cputime
    real(8) :: time_beg
    !
    time_beg=ptime() 
#endif
    !
    call convrsdcal6
    ! call ConvRsdCal6_legc
    !
#ifdef cputime
    ctime(9)=ctime(9)+ptime()-time_beg
#endif
    !
    qrhs=-qrhs
    !
#ifdef cputime
    time_beg=ptime() 
#endif
    !
    call diffrsdcal6
    !
    if(trim(flowtype)=='channel') then 
      call srcchan
    endif
    !
#ifdef cputime
    ctime(10)=ctime(10)+ptime()-time_beg
#endif
    !
  end subroutine rhscal
  !+-------------------------------------------------------------------+
  !| The end of the subroutine rhscal.                                 |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine add a source term to the rsd of the equation to   |
  !| drive channel flow.                                               |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 13-02-2021: Created by J. Fang @ STFC Daresbury Laboratory        |
  !+-------------------------------------------------------------------+
  subroutine srcchan
    !
    use commvar,  only : force,im,jm,km
    use parallel, only : psum 
    use commarray,only : q,qrhs,x,jacob
    !
    ! local data
    integer :: i,j,k,k1,k2
    !
    real(8) :: u1bulk,u2bulk,u3bulk,robulk,dy
    !
    if(ndims==2) then
      k1=0
      k2=0
    elseif(ndims==3) then
      k1=1
      k2=km
    else
      print*,' !! ndims=',ndims
      stop ' !! error @ massfluxchan !!'
    endif
    !
    u1bulk=0.d0
    u2bulk=0.d0
    u3bulk=0.d0
    robulk=0.d0
    !
    do k=k1,k2
    do j=1,jm
    do i=1,im
      !
      dy=x(i,j,k,2)-x(i,j-1,k,2)
      robulk=robulk+0.5d0*(q(i,j-1,k,1)+q(i,j,k,1))*dy
      u1bulk=u1bulk+0.5d0*(q(i,j-1,k,2)+q(i,j,k,2))*dy
      u2bulk=u2bulk+0.5d0*(q(i,j-1,k,3)+q(i,j,k,3))*dy
      u3bulk=u3bulk+0.5d0*(q(i,j-1,k,4)+q(i,j,k,4))*dy
      !
    end do
    end do
    end do
    !
    robulk=psum(robulk)
    u1bulk=psum(u1bulk)/robulk
    u2bulk=psum(u2bulk)/robulk
    u3bulk=psum(u3bulk)/robulk
    !
    do k=0,km
    do j=0,jm
    do i=0,im
      qrhs(i,j,k,2)=qrhs(i,j,k,2)+force(1)*jacob(i,j,k)
      qrhs(i,j,k,3)=qrhs(i,j,k,3)+force(2)*jacob(i,j,k)
      qrhs(i,j,k,4)=qrhs(i,j,k,4)+force(3)*jacob(i,j,k)
      qrhs(i,j,k,5)=qrhs(i,j,k,5)+( force(1)*u1bulk+force(2)*u2bulk+   &
                                    force(3)*u3bulk )*jacob(i,j,k)
    end do
    end do
    end do
    !
  end subroutine srcchan
  !+-------------------------------------------------------------------+
  !| The end of the subroutine source1.                                |
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
                        npdci,npdcj,npdck,is,ie,js,je,ks,ke
    use commarray,only: q,vel,rho,prs,tmp,spc,dxi,jacob,qrhs
    use commfunc, only: ddfc
    !
    ! local data
    integer :: i,j,k,jspc,n
    real(8),allocatable :: fcs(:,:),dfcs(:,:),uu(:) 
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! calculating along i direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    allocate(fcs(-hm:im+hm,1:numq),dfcs(0:im,1:numq),uu(-hm:im+hm))
    do k=ks,ke
    do j=js,je
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
      do n=1,numq
        dfcs(:,n)=ddfc(fcs(:,n),conschm,npdci,im,alfa_con,cci)
      enddo
      !
      qrhs(is:ie,j,k,:)=qrhs(is:ie,j,k,:)+dfcs(is:ie,:) 
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
    do k=ks,ke
    do i=is,ie
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
      do n=1,numq
        dfcs(:,n)=ddfc(fcs(:,n),conschm,npdcj,jm,alfa_con,ccj)
      enddo
      !
      qrhs(i,js:je,k,:)=qrhs(i,js:je,k,:)+dfcs(js:je,:)
      !
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
      do j=js,je
      do i=is,ie
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
        do n=1,numq
          dfcs(:,n)=ddfc(fcs(:,n),conschm,npdck,km,alfa_con,cck,lfft=lfftk)
        enddo
        !
        qrhs(i,j,ks:ke,:)=qrhs(i,j,ks:ke,:)+dfcs(ks:ke,:)
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
                          const5,is,ie,js,je,ks,ke
    use commarray, only : vel,tmp,spc,dvel,dtmp,dspc,dxi,x,jacob,qrhs
    use commfunc,  only : ddfc
    use fludyna,   only : miucal
    use tecio
    !
    ! local data
    integer :: i,j,k,n,ncolm
    real(8),allocatable :: df(:,:),ff(:,:)
    real(8),allocatable,dimension(:,:,:,:) :: sigma,qflux
    real(8) :: miu,miu2,hcc,s11,s12,s13,s22,s23,s33,skk
    !
    dvel=0.d0
    dtmp=0.d0
    dspc=0.d0
    !
    ! calculate velocity and temperature gradient
    !
    ncolm=4+num_species
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
      if(num_species>1) then
        do n=1,num_species
          ff(:,4+n)=spc(:,j,k,n)
        enddo
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
      if(num_species>1) then
        do n=1,num_species
          dspc(:,j,k,n,1)=dspc(:,j,k,n,1)+df(:,4+n)*dxi(0:im,j,k,1,1)
          dspc(:,j,k,n,2)=dspc(:,j,k,n,2)+df(:,4+n)*dxi(0:im,j,k,1,2)
          dspc(:,j,k,n,3)=dspc(:,j,k,n,3)+df(:,4+n)*dxi(0:im,j,k,1,3)
        enddo
      endif
      !
    enddo
    enddo
    !
    deallocate(ff,df)
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
      if(num_species>1) then
        do n=1,num_species
          ff(:,4+n)=spc(i,:,k,n)
        enddo
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
      if(num_species>1) then
        do n=1,num_species
          dspc(i,:,k,n,1)=dspc(i,:,k,n,1)+df(:,4+n)*dxi(i,0:jm,k,2,1)
          dspc(i,:,k,n,2)=dspc(i,:,k,n,2)+df(:,4+n)*dxi(i,0:jm,k,2,2)
          dspc(i,:,k,n,3)=dspc(i,:,k,n,3)+df(:,4+n)*dxi(i,0:jm,k,2,3)
        enddo
      endif
      !
    enddo
    enddo
    deallocate(ff,df)
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
        if(num_species>1) then
          do n=1,num_species
            ff(:,4+n)=spc(i,j,:,n)
          enddo
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
        if(num_species>1) then
          do n=1,num_species
            dspc(i,j,:,n,1)=dspc(i,j,:,n,1)+df(:,4+n)*dxi(i,j,0:km,3,1)
            dspc(i,j,:,n,2)=dspc(i,j,:,n,2)+df(:,4+n)*dxi(i,j,0:km,3,2)
            dspc(i,j,:,n,3)=dspc(i,j,:,n,3)+df(:,4+n)*dxi(i,j,0:km,3,3)
          enddo
        endif
        !
      enddo
      enddo
      deallocate(ff,df)
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
    ncolm=4
    !
    allocate(ff(-hm:im+hm,1:ncolm),df(0:im,1:ncolm))
    do k=ks,ke
    do j=js,je
      !
      ff(:,1)=( sigma(:,j,k,1)*dxi(:,j,k,1,1) +                        &
                sigma(:,j,k,2)*dxi(:,j,k,1,2) +                        &
                sigma(:,j,k,3)*dxi(:,j,k,1,3) )*jacob(:,j,k)
      ff(:,2)=( sigma(:,j,k,2)*dxi(:,j,k,1,1) +                        &
                sigma(:,j,k,4)*dxi(:,j,k,1,2) +                        &
                sigma(:,j,k,5)*dxi(:,j,k,1,3) )*jacob(:,j,k)
      ff(:,3)=( sigma(:,j,k,3)*dxi(:,j,k,1,1) +                        &
                sigma(:,j,k,5)*dxi(:,j,k,1,2) +                        &
                sigma(:,j,k,6)*dxi(:,j,k,1,3) )*jacob(:,j,k)
      ff(:,4)=( qflux(:,j,k,1)*dxi(:,j,k,1,1) +                        &
                qflux(:,j,k,2)*dxi(:,j,k,1,2) +                        &
                qflux(:,j,k,3)*dxi(:,j,k,1,3) )*jacob(:,j,k)
      !
      do n=1,ncolm
        df(:,n)=ddfc(ff(:,n),difschm,npdci,im,alfa_dif,dci)
      enddo
      !
      qrhs(is:ie,j,k,2)=qrhs(is:ie,j,k,2)+df(is:ie,1)
      qrhs(is:ie,j,k,3)=qrhs(is:ie,j,k,3)+df(is:ie,2)
      qrhs(is:ie,j,k,4)=qrhs(is:ie,j,k,4)+df(is:ie,3)
      qrhs(is:ie,j,k,5)=qrhs(is:ie,j,k,5)+df(is:ie,4)
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
    allocate(ff(-hm:jm+hm,1:ncolm),df(0:jm,1:ncolm))
    do k=ks,ke
    do i=is,ie
      !
      ff(:,1)=( sigma(i,:,k,1)*dxi(i,:,k,2,1) +                        &
                sigma(i,:,k,2)*dxi(i,:,k,2,2) +                        &
                sigma(i,:,k,3)*dxi(i,:,k,2,3) )*jacob(i,:,k)
      ff(:,2)=( sigma(i,:,k,2)*dxi(i,:,k,2,1) +                        &
                sigma(i,:,k,4)*dxi(i,:,k,2,2) +                        &
                sigma(i,:,k,5)*dxi(i,:,k,2,3) )*jacob(i,:,k)
      ff(:,3)=( sigma(i,:,k,3)*dxi(i,:,k,2,1) +                        &
                sigma(i,:,k,5)*dxi(i,:,k,2,2) +                        &
                sigma(i,:,k,6)*dxi(i,:,k,2,3) )*jacob(i,:,k)
      ff(:,4)=( qflux(i,:,k,1)*dxi(i,:,k,2,1) +                        &
                qflux(i,:,k,2)*dxi(i,:,k,2,2) +                        &
                qflux(i,:,k,3)*dxi(i,:,k,2,3) )*jacob(i,:,k)
      !
      do n=1,ncolm
        df(:,n)=ddfc(ff(:,n),difschm,npdcj,jm,alfa_dif,dcj)
      enddo
      !
      qrhs(i,js:je,k,2)=qrhs(i,js:je,k,2)+df(js:je,1)
      qrhs(i,js:je,k,3)=qrhs(i,js:je,k,3)+df(js:je,2)
      qrhs(i,js:je,k,4)=qrhs(i,js:je,k,4)+df(js:je,3)
      qrhs(i,js:je,k,5)=qrhs(i,js:je,k,5)+df(js:je,4)
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
      allocate(ff(-hm:km+hm,1:ncolm),df(0:km,1:ncolm))
      do j=js,je
      do i=is,ie
        !
        ff(:,1)=( sigma(i,j,:,1)*dxi(i,j,:,3,1) +                      &
                  sigma(i,j,:,2)*dxi(i,j,:,3,2) +                      &
                  sigma(i,j,:,3)*dxi(i,j,:,3,3) )*jacob(i,j,:)
        ff(:,2)=( sigma(i,j,:,2)*dxi(i,j,:,3,1) +                      &
                  sigma(i,j,:,4)*dxi(i,j,:,3,2) +                      &
                  sigma(i,j,:,5)*dxi(i,j,:,3,3) )*jacob(i,j,:)
        ff(:,3)=( sigma(i,j,:,3)*dxi(i,j,:,3,1) +                      &
                  sigma(i,j,:,5)*dxi(i,j,:,3,2) +                      &
                  sigma(i,j,:,6)*dxi(i,j,:,3,3) )*jacob(i,j,:)
        ff(:,4)=( qflux(i,j,:,1)*dxi(i,j,:,3,1) +                      &
                  qflux(i,j,:,2)*dxi(i,j,:,3,2) +                      &
                  qflux(i,j,:,3)*dxi(i,j,:,3,3) )*jacob(i,j,:)
        !
        do n=1,ncolm
          df(:,n)=ddfc(ff(:,n),difschm,npdck,km,alfa_dif,dck,lfft=lfftk)
        enddo
        !
        qrhs(i,j,ks:ke,2)=qrhs(i,j,ks:ke,2)+df(ks:ke,1)
        qrhs(i,j,ks:ke,3)=qrhs(i,j,ks:ke,3)+df(ks:ke,2)
        qrhs(i,j,ks:ke,4)=qrhs(i,j,ks:ke,4)+df(ks:ke,3)
        qrhs(i,j,ks:ke,5)=qrhs(i,j,ks:ke,5)+df(ks:ke,4)
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
  ! 2-order Explicit filter is incorporated.
  ! Ref1: Datta V. Gaitonde and Miguel R. Visbal, AIAA JOURNAL Vol.38,
  !      No.11, November 2000. 
  ! Ref2: Xavier Gloerfelt and Philippe Lafon, Computers & Fluids, 2008,
  !       37: 388-401.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine spongefilter
    !
    use commvar,  only : spg_imin,spg_imax,spg_jmin,spg_jmax,          &
                         spg_kmin,spg_kmax,im,jm,km,ia,ja,ka,          &
                         lisponge,ljsponge,lksponge,is,je,js,je,ks,ke, &
                         numq
    use commarray,only: lspg_imin,lspg_imax,lspg_jmin,lspg_jmax,       &
                        lspg_kmin,lspg_kmax,x,q
    use commfunc, only : spafilter10
    !
    real(8),parameter :: dampfac=0.05d0
    !
    integer :: i,j,k,n
    real(8) :: dis,var1
    real(8),allocatable :: qtemp(:,:,:,:)
    !
    if(lisponge) then
      !
      call dataswap(q,direction=1)
      !
      !
      if(spg_imax>0) then
        !
        allocate(qtemp(im-spg_imax:im-1,js:je,ks:ke,1:numq))
        !
        do k=ks,ke
        do j=js,je
          !
          dis=0.d0
          do i=im-spg_imax,im-1
            !
            dis=dis+ sqrt( (x(i+1,j,k,1)-x(i,j,k,1))**2+               &
                           (x(i+1,j,k,2)-x(i,j,k,2))**2+               &
                           (x(i+1,j,k,3)-x(i,j,k,3))**2                )
            var1=dampfac*(dis/lspg_imax(j,k))**2
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
          do i=im-spg_imax,im-1
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
    use commvar,  only : im,jm,km,numq,npdci,npdcj,npdck,              &
                         alfa_filter,ndims
    use commarray,only : q
    use commfunc, only : spafilter10
    !
    ! local data
    integer :: i,j,k,n,m
    real(8),allocatable :: phi(:,:),fph(:,:)
    !
#ifdef cputime
    real(8) :: time_beg
    !
    time_beg=ptime()
#endif
    !
    ! filtering in i direction
    call dataswap(q,direction=1)
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
      enddo
      !
      q(0:im,j,k,:)=fph
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
    ! filtering in j direction
    call dataswap(q,direction=2)
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
      enddo
      !
      q(i,0:jm,k,:)=fph
      !
    end do
    end do
    !
    deallocate(phi,fph)
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! end filter in j direction.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    if(ndims==3) then
      !
      call dataswap(q,direction=3)
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
#ifdef cputime
    !
    ctime(8)=ctime(8)+ptime()-time_beg
    !
#endif
    !
    return
    !
  end subroutine filterq
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! End of the subroutine filterq.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!

  !!
  subroutine gradcal
    !
    use commvar,   only : im,jm,km,npdci,npdcj,npdck,conschm,          &
                          alfa_filter,numq
    use commarray, only : x,q,dxi
    use commfunc,  only : ddfc,spafilter10
    use bc,       only : boucon
    !
    ! local data
    integer :: i,j,k,n
    real(8),allocatable :: dq(:,:)
    !
    do k=0,km
    do j=0,jm
    do i=0,im
      ! q(i,j,k,1)=sin(4.d0*x(1,1,k,3)) !+0.1d0*sin(pi*j+0.5d0*pi)
      do n=1,numq
        q(i,j,k,n)=sin(0.5d0*pi*x(i,j,k,2))
        ! sin(10.d0*(x(1,1,k,3)-0.5*pi))*exp(-4.d0*(x(1,1,k,3)-0.5*pi)**2)+0.1d0*sin(pi*k+0.5d0*pi)
      enddo
      !
    enddo
    enddo
    enddo
    !
    call dataswap(q,debug=.true.)
    !
    call boucon
    !
    ! allocate(dq(0:im))
    !
    ! dq=ddfc(q(:,1,1,1),conschm,npdci,im,alfa_con,cci)*dxi(:,1,1,1,1)
    
    ! dq=spafilter10(q(:,jm,0,2),npdci,im,alfa_filter,fci)
    !
    allocate(dq(0:jm,1:1))
    dq(:,1)=ddfc(q(0,:,0,2),conschm,npdcj,jm,alfa_con,ccj)*dxi(0,0:jm,0,2,2)
    !
    ! dq=spafilter10(q(1,:,0,1),npdcj,jm,alfa_filter,fcj)
    !
    ! allocate(dq(0:km,1:numq))
    ! do n=1,numq
    !   ! dq(:,n)=ddfc(q(1,1,:,n),conschm,npdck,km,alfa_con,cck,lfft=.true.)
    !   dq(:,n)=spafilter10(q(1,1,:,n),npdck,km,alfa_filter,fck,lfft=lfftk)
    ! enddo
    !
    ! do n=1,numq
    !   dq(:,n)=dq(:,n)*dxi(1,1,0:km,3,3)
    ! enddo
    ! !
    ! if(mpirank==0) then
    !   do k=0,km
    !     print*,k,dxi(1,1,k,3,3)
    !   enddo
    ! endif
    !
    ! dq=ddfc(q(1,1,:,1),conschm,npdck,km,alfa_con,cck)/(x(1,1,1,3)-x(1,1,0,3))
    !
    open(18,file='testout/profile'//mpirankname//'.dat')
    do j=0,jm
      write(18,*)x(0,j,0,2),q(0,j,0,2),dq(j,1)
    enddo
    close(18)
    print*,' << testout/profile',mpirankname,'.dat'
    !
    deallocate(dq)
    !
  end subroutine gradcal
  !
end module solver
!+---------------------------------------------------------------------+
!| The end of the module commarray.                                    |
!+---------------------------------------------------------------------+