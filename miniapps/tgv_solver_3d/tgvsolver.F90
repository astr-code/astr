module commarray

  implicit none

  real(8) :: dx,dy,dz

  real(8),allocatable,dimension(:,:,:,:,:) :: dvel
  real(8),allocatable,dimension(:,:,:,:) :: q,qrhs,vel,dtmp
  real(8),allocatable,dimension(:,:,:) :: rho,prs,tmp

  contains

  subroutine allocommarray

    use commvar, only : im,jm,km,hm,numq

    ! local data
    integer :: lallo

    allocate( q(-hm:im+hm,-hm:jm+hm,-hm:km+hm,1:numq),stat=lallo)

    allocate( rho(-hm:im+hm,-hm:jm+hm,-hm:km+hm),stat=lallo)

    allocate( prs(-hm:im+hm,-hm:jm+hm,-hm:km+hm),stat=lallo)

    allocate( tmp(-hm:im+hm,-hm:jm+hm,-hm:km+hm),stat=lallo)

    allocate( vel(-hm:im+hm,-hm:jm+hm,-hm:km+hm,1:3),stat=lallo)

    allocate(qrhs(0:im,0:jm,0:km,1:numq),stat=lallo)

    allocate(dvel(0:im,0:jm,0:km,1:3,1:3),stat=lallo)

    allocate(dtmp(0:im,0:jm,0:km,1:3),stat=lallo)

    print*,' ** global arrays allocated'

  end subroutine allocommarray

end module commarray

module fludyna

    implicit none

    contains

    subroutine fvar2q(q,density,velocity,pressure,temperature)

        use commvar, only: const1,const6
        !
        real(8),intent(in) :: density,velocity(:)
        real(8),intent(in),optional :: pressure,temperature
        real(8),intent(out) :: q(:)
        !
        ! local data
        real(8) :: var1
        !
        q(1)=density
        q(2)=density*velocity(1)
        q(3)=density*velocity(2)
        q(4)=density*velocity(3)
        
        var1=0.5d0*sum(velocity(:)*velocity(:))
        if(present(temperature)) then
          q(5)=density*(temperature*const1+var1)
        elseif(present(pressure)) then
          q(5)=pressure*const6+density*var1
        endif

  end subroutine fvar2q

  function thermal_scar(density,pressure,temperature,species) result(vout)

      use commvar,only : const2,rgas
      !
      ! arguments
      real(8) :: vout
      real(8),intent(in) ,optional :: density,pressure,temperature,species(:)
      !
      real(8) :: rloc

      if(present(density) .and. present(temperature)) then
        vout=density*temperature/const2
      elseif(present(density) .and. present(pressure)) then
        vout=pressure/density*const2
      elseif(present(temperature) .and. present(pressure)) then
        vout=pressure/temperature*const2
      endif
  end function thermal_scar

  subroutine q2fvar(q,density,velocity,pressure,temperature)
    !
    use commvar, only: const6
    !
    real(8),intent(in) :: q(:)
    real(8),intent(out) :: density
    real(8),intent(out),optional :: velocity(:),pressure,temperature
    !
    ! local data
    integer :: jspec
    real(8) :: var1
    !
    density   =q(1)
    !
    if(present(velocity) .or. present(pressure) .or. present(temperature)) then
      velocity(1)=q(2)/density
      velocity(2)=q(3)/density
      velocity(3)=q(4)/density
    endif

    if(present(pressure) .or. present(temperature)) then
      pressure  =( q(5)-0.5d0*density*(velocity(1)**2+velocity(2)**2+    &
                                       velocity(3)**2) )/const6
    endif
    !
    if(present(temperature)) then
      temperature=thermal_scar(pressure=pressure,density=density)
    endif

  end subroutine q2fvar

  pure real(8) function miucal(temper)
    !
    use commvar, only :  ref_tem
    !
    real(8),intent(in) :: temper
    ! temper represent temperature, dimensionless
    ! below calculate miucal using sutherland's law
    !
    real(8) :: tempconst,tempconst1
    ! 
    tempconst=110.4d0/ref_tem
    tempconst1=1.d0+tempconst
    !
    miucal=temper*sqrt(temper)*tempconst1/(temper+tempconst)
    !
    return
    !
  end function miucal

end module fludyna

module initilise
    
    implicit none

    integer :: maxstep

    contains

    subroutine program_init

        use commvar, only : im,jm,km,hm,numq,difschm,alfa_filter

        numq=5

        im=128
        jm=128
        km=128

        difschm='643c'

        alfa_filter=0.49d0

        maxstep=10

        write(*,'(3(A,I0))')'  ** dimension set: ',im,' x ',jm,' x ',km

    end subroutine program_init

    subroutine flowfield_init

        use constdef
        use commvar, only: im,jm,km,gamma,mach,reynolds,prandtl,ref_tem, &
                           const1,const2,const3,const4,const5,const6,const7,time,deltat,nstep
        use commarray, only: rho,prs,tmp,q,vel,dx,dy,dz,qrhs
        use fludyna,   only: thermal_scar,fvar2q
        !
        integer :: i,j,k
        real(8) :: var1,pinf,x1,x2,x3
        !
        ref_tem=273.15d0
        gamma=1.4d0
        mach =0.1d0
        reynolds=1600.d0
        prandtl=0.72d0
        !
        const1=1.d0/(gamma*(gamma-1.d0)*mach**2)
        const2=gamma*mach**2
        const3=(gamma-1.d0)/3.d0*prandtl*(mach**2)
        const4=(gamma-1.d0)*mach**2*reynolds*prandtl
        const5=(gamma-1.d0)*mach**2
        const6=1.d0/(gamma-1.d0)
        const7=(gamma-1.d0)*mach**2*Reynolds*prandtl
        !
        pinf=1.d0/const2
        !
        dx=2.d0*pi/dble(im)
        dy=2.d0*pi/dble(jm)
        dz=2.d0*pi/dble(km)
        !
        var1=0.d0
        do k=0,km
        do j=0,jm
        do i=0,im
          !
          x1  =2.d0*pi/dble(im)*dble(i)
          x2  =2.d0*pi/dble(jm)*dble(j)
          x3  =2.d0*pi/dble(km)*dble(k)
          !
          rho(i,j,k)  =1.d0
          vel(i,j,k,1)= sin(x1)*cos(x2)*cos(x3)
          vel(i,j,k,2)=-cos(x1)*sin(x2)*cos(x3)
          vel(i,j,k,3)=0.d0
          prs(i,j,k)  =pinf+1.d0/16.d0*(cos(2.d0*x1)+cos(2.d0*x2))*(cos(2.d0*x3)+2.d0)
          !
          tmp(i,j,k)  =thermal_scar(density=rho(i,j,k),pressure=prs(i,j,k))
          !
          call fvar2q(q=q(i,j,k,:),density=rho(i,j,k),velocity=vel(i,j,k,:),pressure=prs(i,j,k))
          !
        enddo
        enddo
        enddo

        qrhs=0.d0

        nstep=0
        time =0.d0
        deltat=1.d-3

        print*,' ** flow field initilised'

    end subroutine flowfield_init

    subroutine solver_init

        use commvar, only : difschm,im,jm,km,alfa_filter
        use derivative, only : fd_scheme_initiate,fds_compact_i,fds_compact_j, &
                               fds_compact_k,explicit_central,compact_central,fds
        use filter,     only : compact_filter_initiate,filter_coefficient_cal, &
                               filter_i,filter_j,filter_k,filter_ii,filter_jj, &
                               filter_kk

        call fd_scheme_initiate(asolver=fds_compact_i,scheme=difschm,ntype=3,dim=im,dir=1)
        call fd_scheme_initiate(asolver=fds_compact_j,scheme=difschm,ntype=3,dim=jm,dir=2)
        call fd_scheme_initiate(asolver=fds_compact_k,scheme=difschm,ntype=3,dim=km,dir=3)

        allocate(compact_central :: fds)

        print*,' ** finite-difference solver initilised'

        call filter_coefficient_cal(alfa=alfa_filter,beter_halo=1.11d0,beter_bouond=1.09d0)
  
        call compact_filter_initiate(afilter=filter_i,ntype=3,dim=im)
        call compact_filter_initiate(afilter=filter_j,ntype=3,dim=jm)
        call compact_filter_initiate(afilter=filter_k,ntype=3,dim=km)

        print*,' ** low-pass filter initilised'

    end subroutine solver_init

end module initilise

module bc

    implicit none

    contains

    subroutine bchomo
        !
        use commvar, only: im,jm,km,hm
        use commarray, only: rho,prs,tmp,vel,q
        use fludyna,  only: fvar2q
        !
        integer :: i,j,k
        !
        ! b.c. at the i direction
        do k=0,km
        do j=0,jm
          !
          do i=-hm,-1
            rho(i,j,k)  =rho(im+i,j,k)
            vel(i,j,k,:)=vel(im+i,j,k,:)
            tmp(i,j,k)  =tmp(im+i,j,k)
            prs(i,j,k)  =prs(im+i,j,k)
            !
            call fvar2q(q=q(i,j,k,:),density=rho(i,j,k),   &
                                    velocity=vel(i,j,k,:), &
                                    pressure=prs(i,j,k))
          enddo
          !
          do i=im+1,im+hm
            rho(i,j,k)  =rho(i-im,j,k)
            vel(i,j,k,:)=vel(i-im,j,k,:)
            tmp(i,j,k)  =tmp(i-im,j,k)
            prs(i,j,k)  =prs(i-im,j,k)
            !
            call fvar2q(q=q(i,j,k,:),density=rho(i,j,k),   &
                                    velocity=vel(i,j,k,:), &
                                    pressure=prs(i,j,k))
          enddo
          !
        enddo
        enddo
        ! end of applying b.c. along i direction
        !
        ! b.c. at the j direction
        do k=0,km
        do i=0,im
          !
          do j=-hm,-1
            rho(i,j,k)  =rho(i,jm+j,k)
            vel(i,j,k,:)=vel(i,jm+j,k,:)
            tmp(i,j,k)  =tmp(i,jm+j,k)
            prs(i,j,k)  =prs(i,jm+j,k)
            !
            call fvar2q(q=q(i,j,k,:),density=rho(i,j,k),   &
                                    velocity=vel(i,j,k,:), &
                                    pressure=prs(i,j,k))
          enddo
          !
          do j=jm+1,jm+hm
            rho(i,j,k)  =rho(i,j-jm,k)
            vel(i,j,k,:)=vel(i,j-jm,k,:)
            tmp(i,j,k)  =tmp(i,j-jm,k)
            prs(i,j,k)  =prs(i,j-jm,k)
            !
            call fvar2q(q=q(i,j,k,:),density=rho(i,j,k),   &
                                    velocity=vel(i,j,k,:), &
                                    pressure=prs(i,j,k))
          enddo
          !
        enddo
        enddo
        ! end of applying b.c. along j direction
        !
        ! b.c. at the k direction
        do j=0,jm
        do i=0,im
          !
          do k=-hm,-1
            rho(i,j,k)  =rho(i,j,km+k)
            vel(i,j,k,:)=vel(i,j,km+k,:)
            tmp(i,j,k)  =tmp(i,j,km+k)
            prs(i,j,k)  =prs(i,j,km+k)
            !
            call fvar2q(q=q(i,j,k,:),density=rho(i,j,k),   &
                                   velocity=vel(i,j,k,:), &
                                   pressure=prs(i,j,k))
          enddo
          !
          do k=km+1,km+hm
            rho(i,j,k)  =rho(i,j,k-km)
            vel(i,j,k,:)=vel(i,j,k-km,:)
            tmp(i,j,k)  =tmp(i,j,k-km)
            prs(i,j,k)  =prs(i,j,k-km)
            !
            call fvar2q(q=q(i,j,k,:),density=rho(i,j,k),   &
                                    velocity=vel(i,j,k,:), &
                                    pressure=prs(i,j,k))
          enddo
          !
        enddo
        enddo
        ! end of applying b.c. along k direction
        !
    end subroutine bchomo
    !
    subroutine bchomovec(var,dir)
      !
      use commvar, only: im,jm,km,hm
      !
      real(8),intent(inout) :: var(-hm:,-hm:,-hm:,1:)
      integer,intent(in),optional :: dir
      !
      integer :: i,j,k,n,ndir
      
      if(present(dir)) then
        ndir=dir
      else
        ndir=0
      endif
      
      if(ndir==1 .or. ndir==0) then
        ! b.c. at the i direction
        do k=0,km
        do j=0,jm
          !
          do i=-hm,-1
            var(i,j,k,:)  =var(im+i,j,k,:)
          enddo
          !
          do i=im+1,im+hm
            var(i,j,k,:)=var(i-im,j,k,:)
          enddo
          !
        enddo
        enddo
        ! end of applying b.c. along i direction
      endif
      !
      if(ndir==2 .or. ndir==0) then
        ! b.c. at the j direction
        do k=0,km
        do i=0,im
          !
          do j=-hm,-1
            var(i,j,k,:)=var(i,jm+j,k,:)
          enddo
          !
          do j=jm+1,jm+hm
            var(i,j,k,:)=var(i,j-jm,k,:)
          enddo
          !
        enddo
        enddo
        ! end of applying b.c. along j direction
      endif
      !
      if(ndir==3 .or. ndir==0) then
        ! b.c. at the k direction
        do j=0,jm
        do i=0,im
          !
          do k=-hm,-1
            var(i,j,k,:)=var(i,j,km+k,:)
          enddo
          !
          do k=km+1,km+hm
            var(i,j,k,:)=var(i,j,k-km,:)
          enddo
          !
        enddo
        enddo
        ! end of applying b.c. along k direction
      endif
      !
    end subroutine bchomovec

end module

module solver
    
    implicit none
    
    contains

    subroutine gradcal

        use commvar,    only : im,jm,km,hm,difschm,ndims
        use commarray,  only : vel,tmp,dvel,dtmp,dx,dy,dz
        use derivative, only : fds,fds_compact_i,fds_compact_j,fds_compact_k

        ! local data
        integer :: i,j,k,n,ncolm
        real(8),allocatable :: df(:,:),ff(:,:)

        real(8) :: time_beg

        ncolm=4

        ! calculate velocity and temperature gradient

        allocate(ff(-hm:im+hm,ncolm),df(0:im,ncolm))

        do k=0,km
        do j=0,jm

          ff(:,1)=vel(:,j,k,1)
          ff(:,2)=vel(:,j,k,2)
          ff(:,3)=vel(:,j,k,3)
          ff(:,4)=tmp(:,j,k)

          do n=1,ncolm
            df(:,n)=fds%central(fds_compact_i,f=ff(:,n),dim=im)
          enddo

          dvel(:,j,k,1,1)=df(:,1)/dx
          dvel(:,j,k,2,1)=df(:,2)/dx
          dvel(:,j,k,3,1)=df(:,3)/dx

          dtmp(:,j,k,1)=df(:,4)/dx

        enddo
        enddo

        deallocate(ff,df)

        allocate(ff(-hm:jm+hm,ncolm),df(0:jm,ncolm))
        do k=0,km
        do i=0,im

          ff(:,1)=vel(i,:,k,1)
          ff(:,2)=vel(i,:,k,2)
          ff(:,3)=vel(i,:,k,3)
          ff(:,4)=tmp(i,:,k)

          do n=1,ncolm
            df(:,n)=fds%central(fds_compact_j,f=ff(:,n),dim=jm)
          enddo

          dvel(i,:,k,1,2)=df(:,1)/dy
          dvel(i,:,k,2,2)=df(:,2)/dy
          dvel(i,:,k,3,2)=df(:,3)/dy

          dtmp(i,:,k,2)=df(:,4)/dy

        enddo
        enddo
        deallocate(ff,df)
    
        allocate(ff(-hm:km+hm,ncolm),df(0:km,ncolm))
        do j=0,jm
        do i=0,im

          ff(:,1)=vel(i,j,:,1)
          ff(:,2)=vel(i,j,:,2)
          ff(:,3)=vel(i,j,:,3)
          ff(:,4)=tmp(i,j,:)

          do n=1,ncolm
            df(:,n)=fds%central(fds_compact_k,f=ff(:,n),dim=km)
          enddo

          dvel(i,j,:,1,3)=df(:,1)/dz
          dvel(i,j,:,2,3)=df(:,2)/dz
          dvel(i,j,:,3,3)=df(:,3)/dz

          dtmp(i,j,:,3)=df(:,4)/dz

        enddo
        enddo
        deallocate(ff,df)

    end subroutine gradcal

    subroutine convection

        use commvar,  only: im,jm,km,hm,numq
        use commarray,only: q,vel,rho,prs,tmp,qrhs,dx,dy,dz
        use derivative,only : fds,fds_compact_i,fds_compact_j,fds_compact_k
        !
        !
        ! local data
        integer :: i,j,k,n
        real(8),allocatable :: fcs(:,:),dfcs(:,:),uu(:)
        !
        real(8),save :: subtime=0.d0
        !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! calculating along i direction
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        allocate(fcs(-hm:im+hm,1:numq),dfcs(0:im,1:numq),uu(-hm:im+hm))
        do k=0,km
        do j=0,jm
          !
          uu(:)=vel(:,j,k,1)
          fcs(:,1)= q(:,j,k,1)*uu
          fcs(:,2)= q(:,j,k,2)*uu+prs(:,j,k) 
          fcs(:,3)= q(:,j,k,3)*uu 
          fcs(:,4)= q(:,j,k,4)*uu 
          fcs(:,5)= (q(:,j,k,5)+prs(:,j,k) )*uu
    
          do n=1,numq
            dfcs(:,n)=fds%central(fds_compact_i,fcs(:,n),dim=im)
          enddo

          qrhs(0:im,j,k,:)=qrhs(0:im,j,k,:)+dfcs(0:im,:)/dx
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
          uu(:)=vel(i,:,k,2)
          fcs(:,1)=q(i,:,k,1)*uu
          fcs(:,2)=q(i,:,k,2)*uu 
          fcs(:,3)=q(i,:,k,3)*uu+prs(i,:,k) 
          fcs(:,4)=q(i,:,k,4)*uu 
          fcs(:,5)=( q(i,:,k,5)+prs(i,:,k) )*uu
    
          do n=1,numq
            dfcs(:,n)=fds%central(fds_compact_j,fcs(:,n),dim=jm)
          enddo

          qrhs(i,0:jm,k,:)=qrhs(i,0:jm,k,:)+dfcs(0:jm,:)/dy
    
        enddo
        enddo
        deallocate(fcs,dfcs,uu)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! end calculating along j direction
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! calculating along k direction
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        allocate(fcs(-hm:km+hm,1:numq),dfcs(0:km,1:numq),uu(-hm:km+hm))
        do j=0,jm
        do i=0,im
          !
          uu(:)=vel(i,j,:,3)
          fcs(:,1)= q(i,j,:,1)*uu
          fcs(:,2)= q(i,j,:,2)*uu 
          fcs(:,3)= q(i,j,:,3)*uu 
          fcs(:,4)= q(i,j,:,4)*uu+prs(i,j,:) 
          fcs(:,5)= ( q(i,j,:,5)+prs(i,j,:) )*uu

          do n=1,numq
            dfcs(:,n)=fds%central(fds_compact_k,fcs(:,n),dim=km)
          enddo

          qrhs(i,j,0:km,:)=qrhs(i,j,0:km,:)+dfcs(0:km,:)/dz

        enddo
        enddo
        deallocate(fcs,dfcs,uu)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! end calculating along k direction
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !
    end subroutine convection

    subroutine diffusion

        use constdef
        use commvar,  only: im,jm,km,hm,reynolds,prandtl,const5
        use commarray,only: vel,dvel,dtmp,qrhs,tmp,dx,dy,dz
        use fludyna,  only: miucal
        use derivative, only : fds,fds_compact_i,fds_compact_j,fds_compact_k
        use bc, only: bchomovec

        real(8),allocatable,dimension(:,:,:,:),save :: sigma,qflux
        real(8),allocatable :: ff(:,:),df(:,:)

        integer :: i,j,k,n,ncolm
        real(8) :: s11,s12,s13,s22,s23,s33,skk,miu,miu2,hcc

        logical,save :: firstcall=.true.

        if(firstcall) then
          allocate( sigma(-hm:im+hm,-hm:jm+hm,-hm:km+hm,1:6),              &
                    qflux(-hm:im+hm,-hm:jm+hm,-hm:km+hm,1:3) )
          !
          firstcall=.false.
        endif

        do k=0,km
        do j=0,jm
        do i=0,im

            s11=dvel(i,j,k,1,1)
            s12=0.5d0*(dvel(i,j,k,1,2)+dvel(i,j,k,2,1))
            s13=0.5d0*(dvel(i,j,k,1,3)+dvel(i,j,k,3,1))
            s22=dvel(i,j,k,2,2)
            s23=0.5d0*(dvel(i,j,k,2,3)+dvel(i,j,k,3,2))
            s33=dvel(i,j,k,3,3)

            skk=num1d3*(s11+s22+s33)

            miu=miucal(tmp(i,j,k))
            miu2=2.d0*miu
            hcc=(miu/prandtl)/const5

            sigma(i,j,k,1)=miu2*(s11-skk)
            sigma(i,j,k,2)=miu2* s12          
            sigma(i,j,k,3)=miu2* s13          
            sigma(i,j,k,4)=miu2*(s22-skk)
            sigma(i,j,k,5)=miu2* s23          
            sigma(i,j,k,6)=miu2*(s33-skk)
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

        call bchomovec(sigma)

        call bchomovec(qflux)

        ncolm=5

        allocate(ff(-hm:im+hm,2:ncolm),df(0:im,2:ncolm))

        do k=0,km
        do j=0,jm
          !
          ff(:,2)= sigma(:,j,k,1)
          ff(:,3)= sigma(:,j,k,2)
          ff(:,4)= sigma(:,j,k,3)
          ff(:,5)= qflux(:,j,k,1)
          !
          !+------------------------------+
          !|    calculate derivative      |
          !+------------------------------+
          do n=2,ncolm
            df(:,n)=fds%central(fds_compact_i,f=ff(:,n),dim=im)
          enddo
          !
          !+------------------------------+
          !| end of calculate derivative  |
          !+------------------------------+
          !
          qrhs(0:im,j,k,2)=qrhs(0:im,j,k,2)+df(0:im,2)/dx
          qrhs(0:im,j,k,3)=qrhs(0:im,j,k,3)+df(0:im,3)/dx
          qrhs(0:im,j,k,4)=qrhs(0:im,j,k,4)+df(0:im,4)/dx
          qrhs(0:im,j,k,5)=qrhs(0:im,j,k,5)+df(0:im,5)/dx

        enddo
        enddo
        !
        deallocate(ff,df)


        allocate(ff(-hm:jm+hm,2:ncolm),df(0:jm,2:ncolm))
        do k=0,km
        do i=0,im
          !
          ff(:,2)= sigma(i,:,k,2)
          ff(:,3)= sigma(i,:,k,4)
          ff(:,4)= sigma(i,:,k,5)
          ff(:,5)= qflux(i,:,k,2)

          !+------------------------------+
          !|    calculate derivative      |
          !+------------------------------+
          do n=2,ncolm
            df(:,n)=fds%central(fds_compact_j,f=ff(:,n),dim=jm)
          enddo
          !+------------------------------+
          !| end of calculate derivative  |
          !+------------------------------+
          !
          qrhs(i,0:jm,k,2)=qrhs(i,0:jm,k,2)+df(0:jm,2)/dy
          qrhs(i,0:jm,k,3)=qrhs(i,0:jm,k,3)+df(0:jm,3)/dy
          qrhs(i,0:jm,k,4)=qrhs(i,0:jm,k,4)+df(0:jm,4)/dy
          qrhs(i,0:jm,k,5)=qrhs(i,0:jm,k,5)+df(0:jm,5)/dy
          !
        enddo
        enddo
        !
        deallocate(ff,df)

        allocate(ff(-hm:km+hm,2:ncolm),df(0:km,2:ncolm))
        !
        do j=0,jm
        do i=0,im
          !
          do k=-hm,km+hm
            ff(k,2)=sigma(i,j,k,3)
            ff(k,3)=sigma(i,j,k,5)
            ff(k,4)=sigma(i,j,k,6)
            ff(k,5)=qflux(i,j,k,3)
          enddo

          do n=2,ncolm
            df(:,n)=fds%central(fds_compact_k,f=ff(:,n),dim=km)
          enddo

          qrhs(i,j,0:km,2)=qrhs(i,j,0:km,2)+df(0:km,1)/dz
          qrhs(i,j,0:km,3)=qrhs(i,j,0:km,3)+df(0:km,2)/dz
          qrhs(i,j,0:km,4)=qrhs(i,j,0:km,4)+df(0:km,3)/dz
          qrhs(i,j,0:km,5)=qrhs(i,j,0:km,5)+df(0:km,4)/dz
          !
        enddo
        enddo
        !
        deallocate(ff,df)

    end subroutine diffusion

    subroutine filterq
    !
    use commvar,  only : im,jm,km,hm,numq
    use commarray,only : q
    use filter,   only : compact_filter,filter_i,filter_j,filter_k
    use bc, only : bchomovec
    !
    ! local data
    integer :: i,j,k,n,m
    real(8),allocatable :: phi(:,:),fph(:,:)
    !
    real(8) :: time_beg
    real(8),save :: subtime=0.d0
    !
    ! filtering in i direction
    call bchomovec(var=q,dir=1)
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
    ! filtering in j direction
    call bchomovec(var=q,dir=2)
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
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! end filter in j direction.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! filtering in k direction
    call bchomovec(var=q,dir=3)

    allocate(phi(-hm:km+hm,1:numq),fph(0:km,1:numq))
    
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
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! end filter in k direction.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    return
    !
  end subroutine filterq

end module solver

module mainloop

  implicit none

  contains

  subroutine stacal
    !
    use commvar, only: im,jm,km,nstep,time
    use commarray, only: rho,vel,dvel
    !
    integer :: i,j,k
    real(8) :: var1,omega(3),tke,rhom,enst
    !
    logical,save :: firstcall = .true.
    integer,save :: filehand
    !
    if(firstcall) then
      filehand=13
      open(filehand,file='state.dat')
      firstcall=.false.
    endif
    !
    rhom=0.d0
    tke =0.d0
    enst=0.d0
    do k=1,km
    do j=1,jm
    do i=1,im
      ! 
      rhom=rhom+rho(i,j,k)
      !
      var1=vel(i,j,k,1)**2+vel(i,j,k,2)**2+vel(i,j,k,3)**2
      !
      tke=tke+rho(i,j,k)*var1
      !
      omega(1)=dvel(i,j,k,3,2)-dvel(i,j,k,2,3)
      omega(2)=dvel(i,j,k,1,3)-dvel(i,j,k,3,1)
      omega(3)=dvel(i,j,k,2,1)-dvel(i,j,k,1,2)
      !
      enst=enst+rho(i,j,k)*(omega(1)**2+omega(2)**2+omega(3)**2)
      !
    enddo
    enddo
    enddo
    rhom=rhom/dble(im*jm*km)
    tke =0.5d0*tke/dble(im*jm*km)
    enst=0.5d0*enst/dble(im*jm*km)
    !
    ! print*,' ** rho=',rhom,', tke=',tke,', enstrophy=',enst
    !
    write(filehand,*)nstep,time,tke,enst
    !
    if(mod(nstep,10)==0) flush(filehand)
    !
  end subroutine stacal

  subroutine steploop

    use commvar, only: time,nstep,deltat,ctime
    use initilise, only: maxstep
    
    do while(nstep<maxstep)

      call rk3

      nstep=nstep+1
      time =time + deltat

      print*,nstep,time

    enddo

  end subroutine steploop

  subroutine rhscal(comptime)
    !
    use commvar, only: ctime
    use commarray, only: qrhs
    use solver,   only: gradcal,convection,diffusion
    !
    real,intent(inout),optional :: comptime
    
    real :: tstart,tfinish

    if(present(comptime)) then
        call cpu_time(tstart)
    endif
    !
    call gradcal
    !
    call convection
    !
    qrhs=-qrhs
    !
    call diffusion
    !
    if(present(comptime)) then
      call cpu_time(tfinish)
      !
      comptime=comptime+tfinish-tstart
    endif

  end subroutine rhscal

  subroutine rk3(comptime)
    
    use constdef
    use commvar, only: im,jm,km,numq,deltat,nstep,ctime
    use commarray,only: q,qrhs,rho,vel,tmp,prs
    use fludyna,  only: q2fvar
    use bc, only: bchomo
    use solver, only: filterq
    !
    real,intent(inout),optional :: comptime
    
    ! local data
    logical,save :: firstcall = .true.
    real(8),save :: rkcoe(3,3)
    real(8),allocatable,save :: qsave(:,:,:,:)
    integer :: rkstep,i,j,k,m
    !
    real :: tstart,tfinish

    if(present(comptime)) then
        call cpu_time(tstart)
    endif
    !
    if(firstcall) then
      !
      rkcoe(1,1)=1.d0
      rkcoe(2,1)=0.d0
      rkcoe(3,1)=1.d0
      !
      rkcoe(1,2)=0.75d0
      rkcoe(2,2)=0.25d0
      rkcoe(3,2)=0.25d0
      !
      rkcoe(1,3)=num1d3
      rkcoe(2,3)=num2d3
      rkcoe(3,3)=num2d3
      !
      allocate(qsave(0:im,0:jm,0:km,1:numq))
      !
      firstcall=.false.
      !
    endif
    !
    do rkstep=1,3
      !
      qrhs=0.d0

      call bchomo

      call rhscal
      !
      if(rkstep==1) then
        !
        call stacal
        !
        do m=1,numq
          qsave(0:im,0:jm,0:km,m)=q(0:im,0:jm,0:km,m)
        enddo
        !
      endif
      !
      do m=1,numq
        !
        q(0:im,0:jm,0:km,m)=rkcoe(1,rkstep)*qsave(0:im,0:jm,0:km,m)+   &
                            rkcoe(2,rkstep)*q(0:im,0:jm,0:km,m)    +   &
                            rkcoe(3,rkstep)*qrhs(0:im,0:jm,0:km,m)*deltat
        !
      enddo
      !
      call filterq
      !
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
      !
    enddo
    !
    if(present(comptime)) then
      call cpu_time(tfinish)
      !
      comptime=comptime+tfinish-tstart
    endif
    !
  end subroutine rk3

end module mainloop

program boxsolver
    
    use initilise, only: program_init,flowfield_init,solver_init
    use commarray, only: allocommarray
    use mainloop,  only: steploop

    implicit none
    !
    real :: tstart,tfinish
    !
    call cpu_time(tstart)
    
    call program_init

    call solver_init

    call allocommarray

    call flowfield_init

    call steploop
    
    call cpu_time(tfinish)
    !
    ! ctime(1)=tfinish-tstart
    !
    ! call timereport
    !
    print*,' ** the job is done **'
    !
end program boxsolver

