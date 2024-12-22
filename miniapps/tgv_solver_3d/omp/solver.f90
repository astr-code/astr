module solver
   
  use constdef
  use utility, only: time_in_second

  implicit none

  contains
  
  subroutine mainloop
  
    use comvardef, only: time,nstep,deltat,ctime,maxstep
    use utility, only: progress_bar
    use numerics, only: fdm_solver_init,filter_init
    use io, only: write_field

    call fdm_solver_init
    call filter_init

    ! call solver_test
    ! call filter_test

    ! stop
  
    do while(nstep<maxstep)
      
      call rk3(ctime(2))
      
      nstep=nstep+1
      time =time + deltat

      if(mod(nstep,100)==0) then
        call write_field
      endif

      ! call progress_bar(nstep,maxstep,'     ',50)
      
    enddo
  
  end subroutine mainloop

  subroutine solver_test
    
    use comvardef, only: hm,im,x,dx
    use numerics, only: fdm_solver_1d

    real(rtype) :: f(-hm:im+hm,1),df(0:im,1)
    real(rtype) :: dfref,error

    integer :: i

    do i=0,im
      f(i,1)=cos(2._rtype*x(i,0,0,1))
    enddo
    f(-hm:-1,1)=f(im-hm:im-1,1)
    f(im+1:im+hm,1)=f(1:hm,1)

    call fdm_solver_1d(f,df,'i')

    df=df/dx

    open(12,file='data.dat')
    do i=0,im
      write(12,*)i,x(i,0,0,1),df(i,1)
    enddo
    close(12)
    print*,' << data.dat'


    error=0._rtype
    do i=1,im
      dfref=-2._rtype*sin(2._rtype*x(i,0,0,1))
      error=error+(dfref-df(i,1))**2
    enddo
    error=error/real(im,rtype)
    print*,' ** error: ',error

  end subroutine solver_test

  subroutine filter_test
    
    use comvardef, only: hm,im,x
    use numerics,  only: low_pass_filter

    real(rtype) :: f(-hm:im+hm),ff(0:im)
    real(rtype) :: error

    integer :: i

    do i=0,im
      f(i)=cos(2._rtype*x(i,0,0,1))+0.1_rtype*sin(48._rtype*x(i,0,0,1))
    enddo
    f(-hm:-1)=f(im-hm:im-1)
    f(im+1:im+hm)=f(1:hm)

    call low_pass_filter(f,ff,'i',im)

    open(12,file='data_filter.dat')
    do i=0,im
      write(12,*)i,x(i,0,0,1),f(i),ff(i)
    enddo
    close(12)
    print*,'<< data_filter.dat'

    error=0._rtype
    do i=1,im
      error=error+(f(i)-ff(i))**2
    enddo
    error=error/real(im,rtype)
    print*,' ** error: ',error

  end subroutine filter_test

  subroutine rk3(comptime)
  
    use comvardef, only: im,jm,km,numq,q,qrhs,deltat,rho,vel,tmp,prs,nstep,ctime
    use fluids, only: q2fvar
    use bc
    use statistics, only: stacal
  
    real,intent(inout),optional :: comptime

    ! local data
    logical,save :: firstcall = .true.
    real(rtype),save :: rkcoe(3,3)
    real(rtype),allocatable,save :: qsave(:,:,:,:)
    integer :: rkstep,i,j,k,m
  
    real :: tstart,tfinish

    if(present(comptime)) then
      tstart=time_in_second()
    endif
  
    if(firstcall) then
      
      rkcoe(1,1)=1._rtype
      rkcoe(2,1)=0._rtype
      rkcoe(3,1)=1._rtype
      
      rkcoe(1,2)=0.75_rtype
      rkcoe(2,2)=0.25_rtype
      rkcoe(3,2)=0.25_rtype
      
      rkcoe(1,3)=num1d3
      rkcoe(2,3)=num2d3
      rkcoe(3,3)=num2d3
      
      allocate(qsave(0:im,0:jm,0:km,1:numq))
      
      firstcall=.false.
      
    endif
  
    do rkstep=1,3
      
      qrhs=0._rtype
      
      call boundarycondition
      
      call rhscal(ctime(7))
      
      if(rkstep==1) then
        !
        if(mod(nstep,10)==0) call stacal
        !
        !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k,m)
        !$OMP DO
        do k=0,km
        do j=0,jm
        do i=0,im
        do m=1,numq
          qsave(i,j,k,m)=q(i,j,k,m)
        enddo
        enddo
        enddo
        enddo
        !$OMP END DO
        !$OMP END PARALLEL
        !
      endif
      
      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k,m)
      !$OMP DO
      do k=0,km
      do j=0,jm
      do i=0,im
      do m=1,numq
        !
        q(i,j,k,m)=rkcoe(1,rkstep)*qsave(i,j,k,m)+   &
                   rkcoe(2,rkstep)*q(i,j,k,m)    +   &
                   rkcoe(3,rkstep)*qrhs(i,j,k,m)*deltat
      enddo
      enddo
      enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      
      call filterq(ctime(6))
      
      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k)
      !$OMP DO
      do k=0,km
      do j=0,jm
      do i=0,im
        call q2fvar( q=q(i,j,k,:),   density=rho(i,j,k),   &
                                    velocity=vel(i,j,k,:), &
                                    pressure=prs(i,j,k),   &
                                 temperature=tmp(i,j,k) )
      enddo
      enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      
      
    enddo
  
    if(present(comptime)) then
      tfinish=time_in_second()
      
      comptime=comptime+tfinish-tstart
    endif
  
  end subroutine rk3

  subroutine rhscal(comptime)
    !
    use comvardef, only: qrhs,ctime
    !
    real,intent(inout),optional :: comptime
    
    real :: tstart,tfinish

    if(present(comptime)) then
      tstart=time_in_second()
    endif
    !
    call gradcal(ctime(3))
    !
    call convection(ctime(4))
    !
    qrhs=-qrhs
    !
    call diffusion(ctime(5))
    !
    if(present(comptime)) then
      tfinish=time_in_second()
      !
      comptime=comptime+tfinish-tstart
    endif

  end subroutine rhscal

  subroutine gradcal(comptime)
    !
    use comvardef,only: im,jm,km,hm,vel,dvel,tmp,dtmp,dx,dy,dz,ndims
    use numerics, only: fdm_solver_1d
    !
    integer :: i,j,k
    !
    real(rtype),allocatable,dimension(:,:) :: f,df
    real,intent(inout),optional :: comptime
    !
    real :: tstart,tfinish
    !
    !$ save f,df
    !$OMP THREADPRIVATE(f,df)
    !
    !
    if(present(comptime)) then
      tstart=time_in_second()
    endif
    !
    !
    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k)

    allocate(f(-hm:im+hm,4),df(0:im,4))

    !$OMP DO
    do k=0,km
    do j=0,jm
      !
      do i=-hm,im+hm
        f(i,1)=vel(i,j,k,1)
        f(i,2)=vel(i,j,k,2)
        f(i,3)=vel(i,j,k,3)
        f(i,4)=tmp(i,j,k)
      enddo

      call fdm_solver_1d(f,df,'i')

      do i=0,im
        !
        dvel(i,j,k,1,1)=df(i,1)/dx
        dvel(i,j,k,2,1)=df(i,2)/dx
        dvel(i,j,k,3,1)=df(i,3)/dx
        dtmp(i,j,k,1)  =df(i,4)/dx
        !
      enddo
      !
    enddo
    enddo
    !$OMP END DO

    deallocate(f,df)
    !
    !call progress_bar(1,3,'  ** temperature and velocity gradient ',10)
    !
    allocate(f(-hm:jm+hm,4),df(0:jm,4))

    !$OMP DO
    do k=0,km
    do i=0,im
      !
      do j=-hm,jm+hm
        f(j,1)=vel(i,j,k,1)
        f(j,2)=vel(i,j,k,2)
        f(j,3)=vel(i,j,k,3)
        f(j,4)=tmp(i,j,k)
      enddo
      !
      call fdm_solver_1d(f,df,'j')
      !
      do j=0,jm
        !
        dvel(i,j,k,1,2)=df(j,1)/dy
        dvel(i,j,k,2,2)=df(j,2)/dy
        dvel(i,j,k,3,2)=df(j,3)/dy
        dtmp(i,j,k,2)  =df(j,4)/dy
        !
      enddo
      !
      !
    enddo
    enddo
    !$OMP END DO

    deallocate(f,df)
    !
    !call progress_bar(2,3,'  ** temperature and velocity gradient ',10)
    !
    if(ndims==3) then

      allocate(f(-hm:km+hm,4),df(0:km,4))
  
      !$OMP DO
      do j=0,jm
      do i=0,im
        !
        do k=-hm,km+hm
          f(k,1)=vel(i,j,k,1)
          f(k,2)=vel(i,j,k,2)
          f(k,3)=vel(i,j,k,3)
          f(k,4)=tmp(i,j,k)
        enddo
        !
        call fdm_solver_1d(f,df,'k')
        !
        do k=0,km
          !
          dvel(i,j,k,1,3)=df(k,1)/dz
          dvel(i,j,k,2,3)=df(k,2)/dz
          dvel(i,j,k,3,3)=df(k,3)/dz
          dtmp(i,j,k,3)  =df(k,4)/dz
          !
        enddo
        !
        !
      enddo
      enddo
      !$OMP END DO

      deallocate(f,df)

    endif

    !$OMP END PARALLEL

    if(present(comptime)) then
      tfinish=time_in_second()
      !
      comptime=comptime+tfinish-tstart
    endif
    !
    !call progress_bar(3,3,'  ** temperature and velocity gradient ',10)
    !
  end subroutine gradcal
  !
  subroutine convection(comptime)
    !
    use comvardef, only: im,jm,km,hm,numq,prs,vel,q,qrhs,dx,dy,dz,ndims
    use fluids,  only: var2q
    use numerics,  only: fdm_solver_1d
    !
    real,intent(inout),optional :: comptime
    !
    real(rtype),allocatable,dimension(:,:) :: fcs,dfcs
    integer :: i,j,k,n
    !
    real :: tstart,tfinish

    !$ save fcs,dfcs
    !$OMP THREADPRIVATE(fcs,dfcs)
    
    if(present(comptime)) then
      tstart=time_in_second()
    endif
    !
    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k)

    !call progress_bar(0,3,'  ** convection terms ',10)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! calculating along i direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    allocate(fcs(-hm:im+hm,1:numq),dfcs(0:im,1:numq))
    !
    !$OMP DO
    do k=0,km
    do j=0,jm
      !
      do i=-hm,im+hm
        fcs(i,1)=q(i,j,k,1)*vel(i,j,k,1)
        fcs(i,2)=q(i,j,k,2)*vel(i,j,k,1) + prs(i,j,k)
        fcs(i,3)=q(i,j,k,3)*vel(i,j,k,1)
        fcs(i,4)=q(i,j,k,4)*vel(i,j,k,1)
        fcs(i,5)=(q(i,j,k,5)+prs(i,j,k))*vel(i,j,k,1)
      enddo
      !
      call fdm_solver_1d(fcs,dfcs,'i')

      do n=1,numq
        !
        do i=0,im
          !
          qrhs(i,j,k,n)=qrhs(i,j,k,n)+dfcs(i,n)/dx
          !
        enddo
        !
      enddo
      !
    enddo
    enddo
    !$OMP END DO
    !
    deallocate(fcs,dfcs)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! end calculating along i direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !call progress_bar(1,3,'  ** convection terms ',10)
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! calculating along j direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    allocate(fcs(-hm:jm+hm,1:numq),dfcs(0:jm,1:numq))
    !
    !$OMP DO
    do k=0,km
    do i=0,im
      !
      do j=-hm,jm+hm
        fcs(j,1)=q(i,j,k,1)*vel(i,j,k,2)
        fcs(j,2)=q(i,j,k,2)*vel(i,j,k,2)
        fcs(j,3)=q(i,j,k,3)*vel(i,j,k,2) + prs(i,j,k)
        fcs(j,4)=q(i,j,k,4)*vel(i,j,k,2)
        fcs(j,5)=(q(i,j,k,5)+prs(i,j,k))*vel(i,j,k,2)
      enddo
      
      call fdm_solver_1d(fcs,dfcs,'j')

      do n=1,numq

        do j=0,jm
          !
          qrhs(i,j,k,n)=qrhs(i,j,k,n)+dfcs(j,n)/dy
          !
        enddo
        !
      enddo

    enddo
    enddo
    !$OMP END DO
    !
    deallocate(fcs,dfcs)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! end calculating along j direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !call progress_bar(2,3,'  ** convection terms ',10)
    !
    if(ndims==3) then
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! calculating along k direction
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      allocate(fcs(-hm:jm+hm,1:numq),dfcs(0:jm,1:numq))
      !
      !$OMP DO
      do j=0,jm
      do i=0,im
        !
        do k=-hm,im+hm
          fcs(k,1)=q(i,j,k,1)*vel(i,j,k,3)
          fcs(k,2)=q(i,j,k,2)*vel(i,j,k,3)
          fcs(k,3)=q(i,j,k,3)*vel(i,j,k,3)
          fcs(k,4)=q(i,j,k,4)*vel(i,j,k,3) + prs(i,j,k)
          fcs(k,5)=(q(i,j,k,5)+prs(i,j,k))*vel(i,j,k,3)
        enddo
        
        call fdm_solver_1d(fcs,dfcs,'k')
  
        do n=1,numq
  
          do k=0,km
            !
            qrhs(i,j,k,n)=qrhs(i,j,k,n)+dfcs(k,n)/dz
            !
          enddo
          !
        enddo
        !
      enddo
      enddo
      !$OMP END DO
      !
      deallocate(fcs,dfcs)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! end calculating along k direction
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    endif
    !
    !$OMP END PARALLEL

    if(present(comptime)) then
      tfinish=time_in_second()
      !
      comptime=comptime+tfinish-tstart
    endif
    !call progress_bar(3,3,'  ** convection terms ',10)
    !
  end subroutine convection
  !
  subroutine diffusion(comptime)
    !
    use constdef, only: num1d3
    use comvardef, only: reynolds,prandtl,const5,dx,dy,dz,ndims
    use comvardef, only: im,jm,km,hm,vel,dvel,dtmp,qrhs,tmp
    use fluids,  only: miucal
    use numerics,  only: fdm_solver_1d
    use bc, only: bchomovec
    !
    real,intent(inout),optional :: comptime
    !
    real(rtype),allocatable,dimension(:,:,:,:),save :: sigma,qflux
    real(rtype),allocatable :: f(:,:),df(:,:)
    !
    real :: tstart,tfinish
    !
    integer :: i,j,k
    real(rtype) :: s11,s12,s13,s22,s23,s33,skk,miu,miu2,hcc
    logical,save :: firstcall=.true.
    !
    !$ save f,df
    !$OMP THREADPRIVATE(f,df)


    if(present(comptime)) then
      tstart=time_in_second()
    endif
    !
    if(firstcall) then
      allocate( sigma(-hm:im+hm,-hm:jm+hm,-hm:km+hm,1:6),              &
                qflux(-hm:im+hm,-hm:jm+hm,-hm:km+hm,1:3) )
      !
      firstcall=.false.
    endif
    !
    !call progress_bar(0,4,'  ** diffusion terms ',10)
    !
    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k,s11,s12,s13,s22,s23,s33,skk,miu,miu2,hcc)
    !$OMP DO
    do k=0,km
    do j=0,jm
    do i=0,im
      !
      miu=miucal(tmp(i,j,k))/reynolds
      hcc=(miu/prandtl)/const5
      !
      miu2=2._rtype*miu
      !
      s11=dvel(i,j,k,1,1)
      s12=0.5_rtype*(dvel(i,j,k,1,2)+dvel(i,j,k,2,1))
      s13=0.5_rtype*(dvel(i,j,k,1,3)+dvel(i,j,k,3,1))
      s22=dvel(i,j,k,2,2)
      s23=0.5_rtype*(dvel(i,j,k,2,3)+dvel(i,j,k,3,2))
      s33=dvel(i,j,k,3,3)
      !
      skk=num1d3*(s11+s22+s33)
      !
      sigma(i,j,k,1)=miu2*(s11-skk) !s11   
      sigma(i,j,k,2)=miu2* s12      !s12  
      sigma(i,j,k,3)=miu2* s13      !s13   
      sigma(i,j,k,4)=miu2*(s22-skk) !s22   
      sigma(i,j,k,5)=miu2* s23      !s23  
      sigma(i,j,k,6)=miu2*(s33-skk) !s33 
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
      ! 
    enddo
    enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    !
    call bchomovec(sigma)
    !
    call bchomovec(qflux)
    !
    !call progress_bar(1,4,'  ** diffusion terms ',10)
    !
    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k)

    allocate(f(-hm:im+hm,1:4),df(0:im,1:4))
    !
    !$OMP DO
    do k=0,km
    do j=0,jm
      !
      do i=-hm,im+hm
        f(i,1)=sigma(i,j,k,1)
        f(i,2)=sigma(i,j,k,2)
        f(i,3)=sigma(i,j,k,3)
        f(i,4)=qflux(i,j,k,1)
      enddo
      !
      call fdm_solver_1d(f,df,'i')
      !
      do i=0,im
        qrhs(i,j,k,2)=qrhs(i,j,k,2)+df(i,1)/dx
        qrhs(i,j,k,3)=qrhs(i,j,k,3)+df(i,2)/dx
        qrhs(i,j,k,4)=qrhs(i,j,k,4)+df(i,3)/dx
        qrhs(i,j,k,5)=qrhs(i,j,k,5)+df(i,4)/dx
      enddo
      !
    enddo
    enddo
    !$OMP END DO
    !
    deallocate(f,df)
    !
    !call progress_bar(2,4,'  ** diffusion terms ',10)
    !
    allocate(f(-hm:jm+hm,1:4),df(0:jm,1:4))
    !
    !$OMP DO
    do k=0,km
    do i=0,im
      !
      do j=-hm,jm+hm
        f(j,1)=sigma(i,j,k,2)
        f(j,2)=sigma(i,j,k,4)
        f(j,3)=sigma(i,j,k,5)
        f(j,4)=qflux(i,j,k,2)
      enddo
      !
      call fdm_solver_1d(f,df,'j')
      !
      do j=0,jm
        qrhs(i,j,k,2)=qrhs(i,j,k,2)+df(j,1)/dy
        qrhs(i,j,k,3)=qrhs(i,j,k,3)+df(j,2)/dy
        qrhs(i,j,k,4)=qrhs(i,j,k,4)+df(j,3)/dy
        qrhs(i,j,k,5)=qrhs(i,j,k,5)+df(j,4)/dy
      enddo
      !
    enddo
    enddo
    !$OMP END DO
    !
    deallocate(f,df)
    !
    
    if(ndims==3) then

      allocate(f(-hm:km+hm,1:4),df(0:km,1:4))
      !
      !$OMP DO
      do j=0,jm
      do i=0,im
        !
        do k=-hm,km+hm
          f(k,1)=sigma(i,j,k,3)
          f(k,2)=sigma(i,j,k,5)
          f(k,3)=sigma(i,j,k,6)
          f(k,4)=qflux(i,j,k,3)
        enddo
        !
        call fdm_solver_1d(f,df,'k')
        !
        do k=0,km
          qrhs(i,j,k,2)=qrhs(i,j,k,2)+df(k,1)/dz
          qrhs(i,j,k,3)=qrhs(i,j,k,3)+df(k,2)/dz
          qrhs(i,j,k,4)=qrhs(i,j,k,4)+df(k,3)/dz
          qrhs(i,j,k,5)=qrhs(i,j,k,5)+df(k,4)/dz
        enddo
        !
      enddo
      enddo
      !$OMP END DO
      !
      deallocate(f,df)

    endif

    !$OMP END PARALLEL
    !
    if(present(comptime)) then
      tfinish=time_in_second()
      !
      comptime=comptime+tfinish-tstart
    endif
    !
    !call progress_bar(1,4,'  ** diffusion terms ',10)
    !
  end subroutine diffusion
  !
  subroutine filterq(comptime)
    !
    use comvardef, only: im,jm,km,hm,numq,q,ndims
    use numerics,  only: low_pass_filter
    use bc, only: bchomovec
    !
    real,intent(inout),optional :: comptime
    
    real :: tstart,tfinish

    real(rtype),allocatable :: phi(:),fph(:)
    !
    integer :: i,j,k,n
    !
    !$ save phi,fph
    !$OMP THREADPRIVATE(phi,fph)

    if(present(comptime)) then
      tstart=time_in_second()
    endif
    
    call bchomovec(q,dir='i')

    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k)

    allocate(phi(-hm:im+hm),fph(0:im))

    !$OMP DO
    do k=0,km
    do j=0,jm
      !
      do n=1,numq
        !
        do i=-hm,im+hm
          phi(i)=q(i,j,k,n)
        enddo
        !
        call low_pass_filter(phi,fph,'i',im)
        !
        !
        do i=0,im
          q(i,j,k,n)=fph(i)
        enddo

      enddo
      !
    enddo
    enddo
    !$OMP END DO

    deallocate(phi,fph)

    !$OMP END PARALLEL

    call bchomovec(q,dir='j')

    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k)

    allocate(phi(-hm:jm+hm),fph(0:jm))

    !$OMP DO
    do k=0,km
    do i=0,im
      !
      do n=1,numq
        !
        do j=-hm,jm+hm
          phi(j)=q(i,j,k,n)
        enddo
        !
        call low_pass_filter(phi,fph,'j',jm)
        !
        do j=0,jm
          q(i,j,k,n)=fph(j)
        enddo
        !
      enddo
      !
    enddo
    enddo
    !$OMP END DO

    deallocate(phi,fph)

    !$OMP END PARALLEL
    
    if(ndims==3) then

      call bchomovec(q,dir='k')

      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k)

      allocate(phi(-hm:km+hm),fph(0:km))
  
      !$OMP DO
      do j=0,jm
      do i=0,im
        !
        do n=1,numq
          !
          do k=-hm,km+hm
            phi(k)=q(i,j,k,n)
          enddo
          !
          call low_pass_filter(phi,fph,'k',km)
          !
          do k=0,km
            q(i,j,k,n)=fph(k)
          enddo
          !
        enddo
        !
      enddo
      enddo
      !$OMP END DO
  
      deallocate(phi,fph)

      !$OMP END PARALLEL

    endif

    if(present(comptime)) then
      tfinish=time_in_second()
      !
      comptime=comptime+tfinish-tstart
    endif
    !
    !call progress_bar(3,3,'  ** spatial filter ',10)
    !
  end subroutine filterq

end module solver