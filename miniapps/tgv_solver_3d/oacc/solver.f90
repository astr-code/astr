module solver
   
  use constdef
  use utility, only: time_in_second

  implicit none

  contains
  
  subroutine mainloop
  
    use comvardef, only: time,nstep,deltat,ctime,maxstep
    use utility, only: progress_bar
  
    do while(nstep<maxstep)
      
      call rk3(ctime(2))
      
      nstep=nstep+1
      time =time + deltat

      call progress_bar(nstep,maxstep,'     ',50)
      
    enddo
  
  end subroutine mainloop

  subroutine rk3(comptime)
  
    use comvardef, only: im,jm,km,numq,q,qrhs,deltat,rho,vel,tmp,prs,nstep,ctime
    use fluidynamcs, only: q2fvar
    use bc
    use statistics, only: stacal
  
    real,intent(inout),optional :: comptime

    ! local data
    logical,save :: firstcall = .true.
    real(8),save :: rkcoe(3,3)
    real(8),allocatable,save :: qsave(:,:,:,:)
    integer :: rkstep,i,j,k,m
  
    real :: tstart,tfinish

    if(present(comptime)) then
      tstart=time_in_second()
    endif
  
    if(firstcall) then
      
      rkcoe(1,1)=1.d0
      rkcoe(2,1)=0.d0
      rkcoe(3,1)=1.d0
      
      rkcoe(1,2)=0.75d0
      rkcoe(2,2)=0.25d0
      rkcoe(3,2)=0.25d0
      
      rkcoe(1,3)=num1d3
      rkcoe(2,3)=num2d3
      rkcoe(3,3)=num2d3
      
      allocate(qsave(0:im,0:jm,0:km,1:numq))
      
      firstcall=.false.
      
    endif
  
    do rkstep=1,3
      
      qrhs=0.d0
      
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
    !$acc routine(diff6ec) vector

    use comvardef,only: im,jm,km,hm,vel,dvel,tmp,dtmp,dx,dy,dz,ctime
    use numerics, only: diff6ec

    integer :: i,j,k,n

    real(8) :: fi(-hm:im+hm,4),dfi(0:im,4), &
               fj(-hm:im+hm,4),dfj(0:im,4), &
               fk(-hm:im+hm,4),dfk(0:im,4)
    real,intent(inout),optional :: comptime

    real :: tstart,tfinish
    
    if(present(comptime)) then
      tstart=time_in_second()
    endif

    !$acc data copyin(vel,tmp), copyout(dvel,dtmp)

    !$acc parallel loop collapse(2) private(i, fi, dfi)
    do k=0,km
    do j=0,jm
      !
      !$acc loop
      do i=-hm,im+hm
        fi(i,1)=vel(i,j,k,1)
        fi(i,2)=vel(i,j,k,2)
        fi(i,3)=vel(i,j,k,3)
        fi(i,4)=tmp(i,j,k)
      enddo
      !
      call diff6ec(fi(:,1),im,hm,dfi(:,1))
      call diff6ec(fi(:,2),im,hm,dfi(:,2))
      call diff6ec(fi(:,3),im,hm,dfi(:,3))
      call diff6ec(fi(:,4),im,hm,dfi(:,4))
      !
      !$acc loop
      do i=0,im
        !
        dvel(i,j,k,1,1)=dfi(i,1)/dx
        dvel(i,j,k,2,1)=dfi(i,2)/dx
        dvel(i,j,k,3,1)=dfi(i,3)/dx
        dtmp(i,j,k,1)  =dfi(i,4)/dx
        !
      enddo
      !
    enddo
    enddo

    !$acc parallel loop collapse(2) private(j, fj, dfj)
    do k=0,km
    do i=0,im
      !
      !$acc loop
      do j=-hm,jm+hm
        fj(j,1)=vel(i,j,k,1)
        fj(j,2)=vel(i,j,k,2)
        fj(j,3)=vel(i,j,k,3)
        fj(j,4)=tmp(i,j,k)
      enddo
      !
      call diff6ec(fj(:,1),jm,hm,dfj(:,1))
      call diff6ec(fj(:,2),jm,hm,dfj(:,2))
      call diff6ec(fj(:,3),jm,hm,dfj(:,3))
      call diff6ec(fj(:,4),jm,hm,dfj(:,4))
      !
      !$acc loop
      do j=0,jm
        !
        dvel(i,j,k,1,2)=dfj(j,1)/dy
        dvel(i,j,k,2,2)=dfj(j,2)/dy
        dvel(i,j,k,3,2)=dfj(j,3)/dy
        dtmp(i,j,k,2)  =dfj(j,4)/dy
        !
      enddo
      !
      !
    enddo
    enddo

    !$acc parallel loop collapse(2) private(k, fk, dfk)
    do j=0,jm
    do i=0,im
      !
      !$acc loop
      do k=-hm,km+hm
        fk(k,1)=vel(i,j,k,1)
        fk(k,2)=vel(i,j,k,2)
        fk(k,3)=vel(i,j,k,3)
        fk(k,4)=tmp(i,j,k)
      enddo
      !
      call diff6ec(fk(:,1),km,hm,dfk(:,1))
      call diff6ec(fk(:,2),km,hm,dfk(:,2))
      call diff6ec(fk(:,3),km,hm,dfk(:,3))
      call diff6ec(fk(:,4),km,hm,dfk(:,4))
      !
      !$acc loop
      do k=0,km
        !
        dvel(i,j,k,1,3)=dfk(k,1)/dz
        dvel(i,j,k,2,3)=dfk(k,2)/dz
        dvel(i,j,k,3,3)=dfk(k,3)/dz
        dtmp(i,j,k,3)  =dfk(k,4)/dz
        !
      enddo

    enddo
    enddo

    !$acc end data

    if(present(comptime)) then
      tfinish=time_in_second()
      !
      comptime=comptime+tfinish-tstart
    endif

  end subroutine gradcal
  !
  subroutine convection(comptime)
    !$acc routine(diff6ec) vector
    !
    use comvardef, only: im,jm,km,hm,numq,rho,prs,tmp,vel,q,qrhs,dx,dy,dz,ctime
    use fluidynamcs,  only: var2q
    use numerics,  only: diff6ec
    !
    real,intent(inout),optional :: comptime
    !
    real(8) :: fi(-hm:im+hm,5),dfi(0:im,5), &
               fj(-hm:im+hm,5),dfj(0:im,5), &
               fk(-hm:im+hm,5),dfk(0:im,5)

    integer :: i,j,k,n
    !
    real :: tstart,tfinish
    
    if(present(comptime)) then
      tstart=time_in_second()
    endif
    
    !$acc data copyin(q,vel,prs), copyout(qrhs)

    !call progress_bar(0,3,'  ** convection terms ',10)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! calculating along i direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !$acc parallel loop collapse(2) private(i, fi, dfi)
    do k=0,km
    do j=0,jm
      !
      !$acc loop
      do i=-hm,im+hm
        fi(i,1)=q(i,j,k,1)*vel(i,j,k,1)
        fi(i,2)=q(i,j,k,2)*vel(i,j,k,1) + prs(i,j,k)
        fi(i,3)=q(i,j,k,3)*vel(i,j,k,1)
        fi(i,4)=q(i,j,k,4)*vel(i,j,k,1)
        fi(i,5)=(q(i,j,k,5)+prs(i,j,k))*vel(i,j,k,1)
      enddo
      !
      do n=1,numq
        !
        call diff6ec(fi(:,n),im,hm,dfi(:,n))
        !
        !$acc loop
        do i=0,im
          !
          qrhs(i,j,k,n)=qrhs(i,j,k,n)+dfi(i,n)/dx
          !
        enddo
        !
      enddo
      !
    enddo
    enddo
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! end calculating along i direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !call progress_bar(1,3,'  ** convection terms ',10)
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! calculating along j direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !$acc parallel loop collapse(2) private(j, fj, dfj)
    do k=0,km
    do i=0,im
      !
      !$acc loop
      do j=-hm,jm+hm
        fj(j,1)=q(i,j,k,1)*vel(i,j,k,2)
        fj(j,2)=q(i,j,k,2)*vel(i,j,k,2)
        fj(j,3)=q(i,j,k,3)*vel(i,j,k,2) + prs(i,j,k)
        fj(j,4)=q(i,j,k,4)*vel(i,j,k,2)
        fj(j,5)=(q(i,j,k,5)+prs(i,j,k))*vel(i,j,k,2)
      enddo
      !
      do n=1,numq
        !
        call diff6ec(fj(:,n),jm,hm,dfj(:,n))
        !
        !$acc loop
        do j=0,jm
          !
          qrhs(i,j,k,n)=qrhs(i,j,k,n)+dfj(j,n)/dy
          !
        enddo
        !
      enddo
      !
    enddo
    enddo
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! end calculating along j direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !call progress_bar(2,3,'  ** convection terms ',10)
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! calculating along k direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !$acc parallel loop collapse(2) private(k, fk, dfk)
    do j=0,jm
    do i=0,im
      !
      !$acc loop
      do k=-hm,im+hm
        fk(k,1)=q(i,j,k,1)*vel(i,j,k,3)
        fk(k,2)=q(i,j,k,2)*vel(i,j,k,3)
        fk(k,3)=q(i,j,k,3)*vel(i,j,k,3)
        fk(k,4)=q(i,j,k,4)*vel(i,j,k,3) + prs(i,j,k)
        fk(k,5)=(q(i,j,k,5)+prs(i,j,k))*vel(i,j,k,3)
      enddo
      !
      do n=1,numq
        !
        call diff6ec(fk(:,n),km,hm,dfk(:,n))
        !
        !$acc loop
        do k=0,km
          !
          qrhs(i,j,k,n)=qrhs(i,j,k,n)+dfk(k,n)/dz
          !
        enddo
        !
      enddo
      !
    enddo
    enddo
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! end calculating along k direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !$acc end data

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
    use comvardef, only: reynolds,prandtl,const5,dx,dy,dz
    use comvardef, only: im,jm,km,hm,vel,dvel,dtmp,qrhs,tmp,ctime
    use fluidynamcs,  only: miucal
    use numerics,  only: diff6ec
    use bc, only: bchomovec
    !
    real,intent(inout),optional :: comptime
    !
    real(8),allocatable,dimension(:,:,:,:),save :: sigma,qflux
    real(8),allocatable :: f(:,:),df(:,:)
    !
    real :: tstart,tfinish
    !
    integer :: i,j,k
    real(8) :: s11,s12,s13,s22,s23,s33,skk,miu,miu2,hcc
    logical,save :: firstcall=.true.
    
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
    do k=0,km
    do j=0,jm
    do i=0,im
      !
      miu=miucal(tmp(i,j,k))/reynolds
      hcc=(miu/prandtl)/const5
      !
      miu2=2.d0*miu
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
    !
    call bchomovec(sigma)
    !
    call bchomovec(qflux)
    !
    !call progress_bar(1,4,'  ** diffusion terms ',10)
    !
    allocate(f(-hm:im+hm,1:4),df(0:im,1:4))
    !
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
      call diff6ec(f(:,1),im,hm,df(:,1))
      call diff6ec(f(:,2),im,hm,df(:,2))
      call diff6ec(f(:,3),im,hm,df(:,3))
      call diff6ec(f(:,4),im,hm,df(:,4))
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
    !
    deallocate(f,df)
    !
    !call progress_bar(2,4,'  ** diffusion terms ',10)
    !
    allocate(f(-hm:jm+hm,1:4),df(0:jm,1:4))
    !
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
      call diff6ec(f(:,1),jm,hm,df(:,1))
      call diff6ec(f(:,2),jm,hm,df(:,2))
      call diff6ec(f(:,3),jm,hm,df(:,3))
      call diff6ec(f(:,4),jm,hm,df(:,4))
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
    !
    deallocate(f,df)
    !
    !call progress_bar(3,4,'  ** diffusion terms ',10)
    !
    allocate(f(-hm:km+hm,1:4),df(0:km,1:4))
    !
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
      call diff6ec(f(:,1),km,hm,df(:,1))
      call diff6ec(f(:,2),km,hm,df(:,2))
      call diff6ec(f(:,3),km,hm,df(:,3))
      call diff6ec(f(:,4),km,hm,df(:,4))
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
    !
    deallocate(f,df)

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
    use comvardef, only: im,jm,km,hm,numq,q
    use numerics,  only: filter10ec
    !
    real,intent(inout),optional :: comptime
    
    real :: tstart,tfinish

    real(8),allocatable :: phi(:),fph(:)
    !
    integer :: i,j,k,n
    !
    !$ save phi,fph
    !$OMP THREADPRIVATE(phi,fph)

    if(present(comptime)) then
      tstart=time_in_second()
    endif
    !
    !call progress_bar(0,3,'  ** spatial filter ',10)
    !
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
        call filter10ec(phi,im,hm,fph)
        !
        do i=0,im
          q(i,j,k,n)=fph(i)
        enddo
        !
      enddo
      !
    enddo
    enddo
    !$OMP END DO

    deallocate(phi,fph)
    !
    !call progress_bar(1,3,'  ** spatial filter ',10)
    !
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
        call filter10ec(phi,jm,hm,fph)
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
    !
    !call progress_bar(2,3,'  ** spatial filter ',10)
    !
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
        call filter10ec(phi,km,hm,fph)
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
    !
    !$OMP END PARALLEL

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