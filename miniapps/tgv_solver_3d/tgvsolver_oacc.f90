module comvardef
  !
  implicit none
  !
  integer :: im=128,jm=128,km=128,hm=5,numq=5
  ! integer :: im=64,jm=64,km=64,hm=5,numq=5

  real(8),parameter :: pi=4.d0*atan(1.0_8),                            &
                       num1d35 =1.d0/35.d0,  num1d3  =1.d0/3.d0,       &
                       num2d3  =2.d0/3.d0,   num1d24 =1.d0/24.d0,      &
                       num4d3  =4.d0/3.d0,   num1d6  =1.d0/6.d0,       &
                       num1d12 =1.d0/12.d0,  num7d12 =7.d0/12.d0,      &
                       num7d9  =7.d0/9.d0,   num1d36 =1.d0/36.d0,      &
                       num1d60 =1.d0/60.d0,  num65d3 =65.d0/3.d0,      &
                       num20d3 =20.d0/3.d0,  num1d11 =1.d0/11.d0,      &
                       num25d12=25.d0/12.d0, num11d6 =11.d0/6.d0,      &
                       num1d840=1.d0/840.d0, num13d60=13.d0/60.d0,     &
                       num1d30 =1.d0/30.d0,  num47d60=47.d0/60.d0,     &
                       num5d6  =5.d0/6.d0,   num1d18 =1.d0/18.d0,      &
                       num19d18=19.d0/18.d0, num5d9  =5.d0/9.d0,       &
                       num9d36 =9.d0/36.d0
  !
  integer :: nstep
  real(8) :: time,deltat
  !
  real(8) :: gamma,mach,reynolds,prandtl,ref_t
  real(8) :: const1,const2,const3,const4,const5,const6,const7
  real(8) :: dx,dy,dz
  !
  real :: ctime(12)
  !
  logical, parameter :: performval1 = .true.
  ! logical, parameter :: performval1 = .false.
  !
  real(8),allocatable,dimension(:,:,:) :: rho,prs,tmp
  real(8),allocatable,dimension(:,:,:,:) :: x,q,qrhs,qsave,vel,dtmp,sigma,qflux
  real(8),allocatable,dimension(:,:,:,:,:) :: dvel
  !
  contains
  !
  subroutine alloarray
    !
    allocate( x(0:im,0:jm,0:km,1:3))
    allocate( q(-hm:im+hm,-hm:jm+hm,-hm:km+hm,1:numq),qsave(0:im,0:jm,0:km,1:numq))
    allocate(qrhs(0:im,0:jm,0:km,1:numq))
    !
    allocate( rho(-hm:im+hm,-hm:jm+hm,-hm:km+hm))
    allocate( vel(-hm:im+hm,-hm:jm+hm,-hm:km+hm,1:3))
    allocate( prs(-hm:im+hm,-hm:jm+hm,-hm:km+hm))
    allocate( tmp(-hm:im+hm,-hm:jm+hm,-hm:km+hm))
    allocate(dvel(0:im,0:jm,0:km,1:3,1:3))
    allocate(dtmp(0:im,0:jm,0:km,1:3))
    allocate( sigma(-hm:im+hm,-hm:jm+hm,-hm:km+hm,1:6),              &
              qflux(-hm:im+hm,-hm:jm+hm,-hm:km+hm,1:3) )
    !
    ! print*,' ** common array allocated'
    !
  end subroutine alloarray
  !
end module comvardef
!
module utlity
  !
  implicit none
  !
  contains
  !
  !+-------------------------------------------------------------------+
  !| Progress indicators library.                                      |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 31-01-2024  | copied by J. Fang via:                              |
  !| https://github.com/macie/fortran-libs                             |
  !|  Maciej Żok, 2010 MIT License                                     |
  !+-------------------------------------------------------------------+
  subroutine progress_bar(iteration,maximum,info2show,barlength)
     !
     ! Prints progress bar.
     !
     ! Args: 
     !     iteration - iteration number
     !     maximum - total iterations
     !     barlength - length of the bar
     !   
     ! use iso_fortran_env
     integer,intent(in) :: iteration,maximum
     character(len=*),intent(in),optional :: info2show
     integer,intent(in),optional :: barlength
     integer :: counter,nlength
     integer :: done
     real(4) :: perc
     !
     if(present(barlength)) then
         nlength=barlength
     else
         nlength=10
     endif
     !
     perc = 100.0*real(iteration)/real(maximum)
     done = floor(perc/(100.0/real(nlength)))  ! mark length
     !
     write(6,'(1A1,A,A)',advance='no')char(13),info2show,'['
     if (done .LE. 0) then
         do counter = 1, nlength
             write(6,'(1A1,A)',advance='no')'='
         end do
     else if ((done .GT. 0) .and. (done .LT. nlength)) then
         do counter = 1, done
             write(6,'(1A1,A)',advance='no')'>'
         end do
         do counter = done+1, nlength
             write(6,'(1A1,A)',advance='no')'='
         end do 
     else
         do counter = 1, nlength
             write(6,'(1A1,A)',advance='no')'>'
         end do
     end if
     write(6,'(A,F5.1,A)',advance='no')'] ',perc,'%'
     !
     if(iteration==maximum) write(6,*)
     !
  end subroutine progress_bar
  !+-------------------------------------------------------------------+
  !| The end of the subroutine progress_bar.                           |
  !+-------------------------------------------------------------------+
  !
end module utlity
!
module dataoper
  !
  implicit none
  !
  contains
  !
  function thermal_scar(density,pressure,temperature) result(vout)
    !
    use comvardef,only : const2
    !
    ! arguments
    real(8) :: vout
    real(8),intent(in) ,optional :: density,pressure,temperature
    !
    real(8) :: rloc
    !
    if(present(density) .and. present(temperature)) then
      vout=density*temperature/const2
    elseif(present(density) .and. present(pressure)) then
      vout=pressure/density*const2
    elseif(present(temperature) .and. present(pressure)) then
      vout=pressure/temperature*const2
    else
      stop ' !! unable to get thermal variable  @ thermal_scar !!'
    endif
    !
  end function thermal_scar
  !
  function var2q(density,velocity,pressure,temperature) result(q)
    !
    use comvardef, only: const1,const6,numq
    !
    real(8) :: q(1:numq)
    real(8),intent(in) :: density,velocity(3)
    real(8),intent(in),optional :: pressure,temperature
    !
    ! local data
    integer :: jspec,j
    real(8) :: var1,var2
    !
    q(1)=density
    q(2)=density*velocity(1)
    q(3)=density*velocity(2)
    q(4)=density*velocity(3)
    !
    var1=0.5d0*sum(velocity(:)*velocity(:))
    if(present(temperature)) then
      q(5)=density*(temperature*const1+var1)
    elseif(present(pressure)) then
      q(5)=pressure*const6+density*var1
    endif
    !
  end function var2q
  !
  subroutine q2fvar(q,density,velocity,pressure,temperature)
    !
    use comvardef, only: const6
    !
    real(8),intent(in) :: q(:)
    real(8),intent(out) :: density,velocity(:),pressure,temperature
    !
    density   =q(1)
    !
    velocity(1)=q(2)/density
    velocity(2)=q(3)/density
    velocity(3)=q(4)/density
    !
    pressure  =( q(5)-0.5d0*density*(velocity(1)**2+velocity(2)**2+  &
                                       velocity(3)**2) )/const6
    !
    temperature=thermal_scar(pressure=pressure,density=density)
    !
  end subroutine q2fvar
  !
  pure real(8) function miucal(temper)
    !
    use comvardef, only :  ref_t
    !
    real(8),intent(in) :: temper
    ! temper represent temperature, dimensionless
    ! below calculate miucal using sutherland's law
    !
    real(8) :: tempconst,tempconst1
    ! 
    tempconst=110.4d0/ref_t
    tempconst1=1.d0+tempconst
    !
    miucal=temper*sqrt(temper)*tempconst1/(temper+tempconst)
    !
    return
    !
  end function miucal
  !
  subroutine tgvinit
    !
    use comvardef, only: im,jm,km,pi,gamma,mach,reynolds,prandtl,dx,dy,dz, &
                         const1,const2,const3,const4,const5,const6,    &
                         const7,time,nstep,deltat,ref_t,               &
                         x,qrhs,q,vel,rho,prs,tmp,ctime
    use utlity, only: progress_bar
    !
    integer :: i,j,k
    real(8) :: var1,pinf
    !
    ref_t=273.15d0
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
      x(i,j,k,1)  =2.d0*pi/dble(im)*dble(i)
      x(i,j,k,2)  =2.d0*pi/dble(jm)*dble(j)
      x(i,j,k,3)  =2.d0*pi/dble(km)*dble(k)
      !
      rho(i,j,k)  =1.d0
      vel(i,j,k,1)= sin(x(i,j,k,1))*cos(x(i,j,k,2))*cos(x(i,j,k,3))
      vel(i,j,k,2)=-cos(x(i,j,k,1))*sin(x(i,j,k,2))*cos(x(i,j,k,3))
      vel(i,j,k,3)=0.d0
      prs(i,j,k)  =pinf+1.d0/16.d0*(cos(2.d0*x(i,j,k,1))+cos(2.d0*x(i,j,k,2)))*(cos(2.d0*x(i,j,k,3))+2.d0)
      !
      tmp(i,j,k)  =thermal_scar(density=rho(i,j,k),pressure=prs(i,j,k))
      !
      q(i,j,k,:)=var2q(density=rho(i,j,k),velocity=vel(i,j,k,:),       &
                       pressure=prs(i,j,k))
      !
    enddo
    enddo
      !call progress_bar(k,km,'  ** initilising data ',20)
    enddo
    !
    qrhs=0.d0
    !
    nstep=0
    time =0.d0
    deltat=1.d-3
    !
    ctime=0.0
    !
    ! print*,' ** data initilised'
    !
  end subroutine tgvinit
  !
end module dataoper
!
module numerics
  !
  use comvardef, only : num1d60
  !
  implicit none
  !
  contains
  !
  subroutine diff6ec(vin,dim,n,vout)
    !$acc routine seq
    !
    integer,intent(in) :: dim,n
    real(8),intent(in) :: vin(-n:dim+n)
    real(8) :: vout(0:dim)
    !
    ! local data
    integer :: i

    do i=0,dim
      vout(i)  =0.75d0 *(vin(i+1)-vin(i-1))- &
                0.15d0 *(vin(i+2)-vin(i-2))+ &
            1.d0/60.d0 *(vin(i+3)-vin(i-3))
    enddo
    !
  end subroutine diff6ec
  !
  subroutine diff6ec_vec(vin,dim,n,vout)
    !$acc routine vector
    !
    integer,intent(in) :: dim,n
    real(8),intent(in) :: vin(-n:dim+n)
    real(8) :: vout(0:dim)
    !
    ! local data
    integer :: i
    !
    !$acc parallel loop vector
    do i=0,dim
      vout(i)  =0.75d0 *(vin(i+1)-vin(i-1))- &
                0.15d0 *(vin(i+2)-vin(i-2))+ &
               num1d60 *(vin(i+3)-vin(i-3))
    enddo
    !
  end subroutine diff6ec_vec
  !
  subroutine filter10ec(vin,dim,n,vout)
    !
    !$acc routine seq
    !
    integer,intent(in) :: dim,n
    real(8),intent(in) :: vin(-n:dim+n)
    real(8) :: vout(0:dim)
    !
    integer :: i
    !
    do i=0,dim
      !
      vout(i)=   0.376953125d0*(vin(i)  +vin(i))     &
               + 0.205078125d0*(vin(i-1)+vin(i+1))   &
               -   0.1171875d0*(vin(i-2)+vin(i+2))   &
               +0.0439453125d0*(vin(i-3)+vin(i+3))   &
               - 0.009765625d0*(vin(i-4)+vin(i+4))   &
               +0.0009765625d0*(vin(i-5)+vin(i+5))
    enddo
    !
  end subroutine filter10ec
  !
end module numerics
  !
module solver
  !
  use utlity,   only: progress_bar
  !
  implicit none
  !
  contains
  !
  subroutine bchomo
    !
    use comvardef, only: im,jm,km,hm,rho,prs,tmp,vel,q
    use dataoper,  only: var2q
    !
    integer :: i,j,k
    !
    ! b.c. at the i direction
    !$acc parallel loop gang collapse(2) present(rho,prs,tmp,vel,q)
    do k=0,km
    do j=0,jm
      !
      do i=-hm,-1
        rho(i,j,k)  =rho(im+i,j,k)
        vel(i,j,k,:)=vel(im+i,j,k,:)
        tmp(i,j,k)  =tmp(im+i,j,k)
        prs(i,j,k)  =prs(im+i,j,k)
        !
        q(i,j,k,:)  =var2q(density=rho(i,j,k),   &
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
        q(i,j,k,:)  =var2q(density=rho(i,j,k),   &
                          velocity=vel(i,j,k,:), &
                          pressure=prs(i,j,k))
      enddo
      !
    enddo
    enddo
    ! end of applying b.c. along i direction
    !
    ! b.c. at the j direction
    !$acc parallel loop gang collapse(2) present(rho,prs,tmp,vel,q)
    do k=0,km
    do i=0,im
      !
      do j=-hm,-1
        rho(i,j,k)  =rho(i,jm+j,k)
        vel(i,j,k,:)=vel(i,jm+j,k,:)
        tmp(i,j,k)  =tmp(i,jm+j,k)
        prs(i,j,k)  =prs(i,jm+j,k)
        !
        q(i,j,k,:)  =var2q(density=rho(i,j,k),   &
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
        q(i,j,k,:)  =var2q(density=rho(i,j,k),   &
                          velocity=vel(i,j,k,:), &
                          pressure=prs(i,j,k))
      enddo
      !
    enddo
    enddo
    ! end of applying b.c. along j direction
    !
    ! b.c. at the k direction
    !$acc parallel loop gang collapse(2) present(rho,prs,tmp,vel,q)
    do j=0,jm
    do i=0,im
      !
      do k=-hm,-1
        rho(i,j,k)  =rho(i,j,km+k)
        vel(i,j,k,:)=vel(i,j,km+k,:)
        tmp(i,j,k)  =tmp(i,j,km+k)
        prs(i,j,k)  =prs(i,j,km+k)
        !
        q(i,j,k,:)  =var2q(density=rho(i,j,k),   &
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
        q(i,j,k,:)  =var2q(density=rho(i,j,k),   &
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
  subroutine bchomovec(var)
    !
    use comvardef, only: im,jm,km,hm
    use dataoper,  only: var2q
    !
    real(8),intent(inout) :: var(-hm:,-hm:,-hm:,1:)
    !
    integer :: i,j,k,n,nd4
    !
    nd4=size(var,4)
    !
    ! b.c. at the i direction
    !$acc parallel loop gang collapse(2) present(var)
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
    !
    ! b.c. at the j direction
    !$acc parallel loop gang collapse(2) present(var)
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
    !
    ! b.c. at the k direction
    !$acc parallel loop gang collapse(2) present(var)
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
    !
  end subroutine bchomovec
  !
  subroutine rhscal(comptime)
    !
    use comvardef, only: im,jm,km,numq,vel,tmp,dvel,dtmp,qrhs,rho,prs,q,ctime
    !
    real,intent(inout),optional :: comptime
    integer :: i,j,k,m
    
    real :: tstart,tfinish

    if(present(comptime)) then
        call cpu_time(tstart)
    endif
    !
    call gradcal(ctime(3))
    !
    call convection(ctime(4))
    !
    call diffusion(ctime(5))
    !
    if(present(comptime)) then
      call cpu_time(tfinish)
      !
      comptime=comptime+tfinish-tstart
    endif
    !
  end subroutine rhscal
  !
  subroutine gradcal(comptime)
    !
    use comvardef,only: im,jm,km,hm,vel,dvel,tmp,dtmp,x,dx,dy,dz,ctime
    use numerics, only: diff6ec,diff6ec_vec
    !
    real,intent(inout),optional :: comptime
    !
    integer :: i,j,k,n
    !
    real(8) :: fi(-hm:im+hm,4),fj(-hm:jm+hm,4),fk(-hm:km+hm,4), &
               dfi(0:im,4),dfj(0:jm,4),dfk(0:km,4)
    !
    ! real(8) :: du(0:im,0:jm,0:km,1:3),er1,er2,er3,e2r1,e2r2,e2r3,emr1,emr2,emr3
    !
    real :: tstart,tfinish
    !
    if(present(comptime)) then
        call cpu_time(tstart)
    endif
    !
    !$acc data create(fi,fj,fk,dfi,dfj,dfk)
    !
    !call progress_bar(0,3,'  ** temperature and velocity gradient ',10)
    !
    ! !$acc parallel loop gang collapse(3) present(x,du) 
    ! do k=0,km
    ! do j=0,jm
    ! do i=0,im
    !   du(i,j,k,1)= cos(x(i,j,k,1))*cos(x(i,j,k,2))*cos(x(i,j,k,3))
    !   du(i,j,k,2)=-sin(x(i,j,k,1))*sin(x(i,j,k,2))*cos(x(i,j,k,3))
    !   du(i,j,k,3)=-sin(x(i,j,k,1))*cos(x(i,j,k,2))*sin(x(i,j,k,3))
    ! enddo
    ! enddo
    ! enddo
    !
    !$acc parallel loop gang collapse(2) present(vel,tmp,dvel,dtmp) private(fi,dfi) 
    !
    do k=0,km
    do j=0,jm
      !
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
    !
    !call progress_bar(1,3,'  ** temperature and velocity gradient ',10)
    !
    !
    !$acc parallel loop gang collapse(2) present(vel,tmp,dvel,dtmp) private(fj,dfj)
    !
    do k=0,km
    do i=0,im
      !
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
    !
    !call progress_bar(2,3,'  ** temperature and velocity gradient ',10)
    !
    !$acc parallel loop gang collapse(2) present(vel,tmp,dvel,dtmp) private(fk,dfk)
    !
    do j=0,jm
    do i=0,im
      !
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
      do k=0,km
        !
        dvel(i,j,k,1,3)=dfk(k,1)/dz
        dvel(i,j,k,2,3)=dfk(k,2)/dz
        dvel(i,j,k,3,3)=dfk(k,3)/dz
        dtmp(i,j,k,3)  =dfk(k,4)/dz
        !
      enddo
      !
    enddo
    enddo
    !
    !$acc end data
    !
    if(present(comptime)) then
      call cpu_time(tfinish)
      !
      comptime=comptime+tfinish-tstart
    endif
    !
    !call progress_bar(3,3,'  ** temperature and velocity gradient ',10)
    !
    ! e2r1=0.d0
    ! e2r2=0.d0
    ! e2r3=0.d0
    ! emr1=0.d0
    ! emr2=0.d0
    ! emr3=0.d0
    ! do k=0,km
    ! do j=0,jm
    ! do i=0,im
    !   er1=du(i,j,k,1)-dvel(i,j,k,1,1)
    !   er2=du(i,j,k,2)-dvel(i,j,k,1,2)
    !   er3=du(i,j,k,3)-dvel(i,j,k,1,3)
    !   !
    !   e2r1=e2r1+er1**2
    !   e2r2=e2r2+er2**2
    !   e2r3=e2r3+er3**2
    !   !
    !   emr1=max(emr1,abs(er1))
    !   emr2=max(emr2,abs(er2))
    !   emr3=max(emr3,abs(er3))
    ! enddo
    ! enddo
    ! enddo
    ! !
    ! print*,' L2 error  :',e2r1,e2r2,e2r3
    ! print*,' Lmax error:',emr1,emr2,emr3
    !
  end subroutine gradcal
  !

  subroutine convection(comptime)
    !
    use comvardef, only: im,jm,km,hm,prs,vel,q,qrhs,dx,dy,dz,ctime
    use dataoper,  only: var2q
    use numerics,  only: diff6ec
    !
    real,intent(inout),optional :: comptime
    !
    real(8) :: fi(-hm:im+hm,5),dfi(0:im,5),fj(-hm:jm+hm,5),dfj(0:jm,5), &
               fk(-hm:km+hm,5),dfk(0:km,5)
    integer :: i,j,k,n
    !
    real :: tstart,tfinish
    !
    !$acc data create(fi,fj,fk,dfi,dfj,dfk)
    !
    if(present(comptime)) then
        call cpu_time(tstart)
    endif
    !
    !call progress_bar(0,3,'  ** convection terms ',10)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! calculating along i direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !$acc parallel loop gang collapse(2) present(q,vel,prs,qrhs) private(fi,dfi)
    do k=0,km
    do j=0,jm
      !
      do i=-hm,im+hm
        fi(i,1)=q(i,j,k,1)*vel(i,j,k,1)
        fi(i,2)=q(i,j,k,2)*vel(i,j,k,1) + prs(i,j,k)
        fi(i,3)=q(i,j,k,3)*vel(i,j,k,1)
        fi(i,4)=q(i,j,k,4)*vel(i,j,k,1)
        fi(i,5)=(q(i,j,k,5)+prs(i,j,k))*vel(i,j,k,1)
      enddo
      !
      do n=1,5
        !
        call diff6ec(fi(:,n),im,hm,dfi(:,n))
        !
        do i=0,im
          !
          qrhs(i,j,k,n)=qrhs(i,j,k,n)-dfi(i,n)/dx
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
    !$acc parallel loop gang collapse(2) present(q,vel,prs,qrhs) private(fj,dfj)
    do k=0,km
    do i=0,im
      !
      do j=-hm,jm+hm
        fj(j,1)=q(i,j,k,1)*vel(i,j,k,2)
        fj(j,2)=q(i,j,k,2)*vel(i,j,k,2)
        fj(j,3)=q(i,j,k,3)*vel(i,j,k,2) + prs(i,j,k)
        fj(j,4)=q(i,j,k,4)*vel(i,j,k,2)
        fj(j,5)=(q(i,j,k,5)+prs(i,j,k))*vel(i,j,k,2)
      enddo
      !
      do n=1,5
        !
        call diff6ec(fj(:,n),jm,hm,dfj(:,n))
        !
        do j=0,jm
          !
          qrhs(i,j,k,n)=qrhs(i,j,k,n)-dfj(j,n)/dy
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
    !$acc parallel loop gang collapse(2) present(q,vel,prs,qrhs) private(fk,dfk)
    do j=0,jm
    do i=0,im
      !
      do k=-hm,im+hm
        fk(k,1)=q(i,j,k,1)*vel(i,j,k,3)
        fk(k,2)=q(i,j,k,2)*vel(i,j,k,3)
        fk(k,3)=q(i,j,k,3)*vel(i,j,k,3)
        fk(k,4)=q(i,j,k,4)*vel(i,j,k,3) + prs(i,j,k)
        fk(k,5)=(q(i,j,k,5)+prs(i,j,k))*vel(i,j,k,3)
      enddo
      !
      do n=1,5
        !
        call diff6ec(fk(:,n),km,hm,dfk(:,n))
        !
        do k=0,km
          !
          qrhs(i,j,k,n)=qrhs(i,j,k,n)-dfk(k,n)/dz
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
    !
    if(present(comptime)) then
      call cpu_time(tfinish)
      !
      comptime=comptime+tfinish-tstart
    endif
    !call progress_bar(3,3,'  ** convection terms ',10)
    !
  end subroutine convection
  !
  subroutine diffusion(comptime)
    !
    use comvardef, only: num1d3,reynolds,prandtl,const5,dx,dy,dz
    use comvardef, only: im,jm,km,hm,vel,dvel,dtmp,qrhs,tmp,ctime,sigma,qflux
    use dataoper,  only: miucal
    use numerics,  only: diff6ec
    !
    real,intent(inout),optional :: comptime
    !
    real(8) :: fi(-hm:im+hm,4),dfi(0:im,4),fj(-hm:jm+hm,4),dfj(0:jm,4), &
               fk(-hm:km+hm,4),dfk(0:km,4)
    !
    real :: tstart,tfinish
    !
    integer :: i,j,k,n
    real(8) :: s11,s12,s13,s22,s23,s33,skk,miu,miu2,hcc
    logical,save :: firstcall=.true.
    !
    if(present(comptime)) then
        call cpu_time(tstart)
    endif
    !
    !$acc data create(fi,fj,fk,dfi,dfj,dfk)
    !
    !call progress_bar(0,4,'  ** diffusion terms ',10)
    !
    !$acc parallel loop collapse(3) present(dvel,vel,dtmp,tmp,sigma,qflux) &
    !$acc                           private(miu,miu2,hcc,s11,s12,s13,s22,s23,s33,skk)
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
    !$acc parallel loop gang collapse(2) present(sigma,qflux,qrhs) private(fi,dfi)
    do k=0,km
    do j=0,jm
      !
      do i=-hm,im+hm
        fi(i,1)=sigma(i,j,k,1)
        fi(i,2)=sigma(i,j,k,2)
        fi(i,3)=sigma(i,j,k,3)
        fi(i,4)=qflux(i,j,k,1)
      enddo
      !
      call diff6ec(fi(:,1),im,hm,dfi(:,1))
      call diff6ec(fi(:,2),im,hm,dfi(:,2))
      call diff6ec(fi(:,3),im,hm,dfi(:,3))
      call diff6ec(fi(:,4),im,hm,dfi(:,4))
      !
      do i=0,im
        qrhs(i,j,k,2)=qrhs(i,j,k,2)+dfi(i,1)/dx
        qrhs(i,j,k,3)=qrhs(i,j,k,3)+dfi(i,2)/dx
        qrhs(i,j,k,4)=qrhs(i,j,k,4)+dfi(i,3)/dx
        qrhs(i,j,k,5)=qrhs(i,j,k,5)+dfi(i,4)/dx
      enddo
      !
    enddo
    enddo
    !
    !call progress_bar(2,4,'  ** diffusion terms ',10)
    !
    !$acc parallel loop gang collapse(2) present(sigma,qflux,qrhs) private(fj,dfj)
    do k=0,km
    do i=0,im
      !
      do j=-hm,jm+hm
        fj(j,1)=sigma(i,j,k,2)
        fj(j,2)=sigma(i,j,k,4)
        fj(j,3)=sigma(i,j,k,5)
        fj(j,4)=qflux(i,j,k,2)
      enddo
      !
      call diff6ec(fj(:,1),jm,hm,dfj(:,1))
      call diff6ec(fj(:,2),jm,hm,dfj(:,2))
      call diff6ec(fj(:,3),jm,hm,dfj(:,3))
      call diff6ec(fj(:,4),jm,hm,dfj(:,4))
      !
      do j=0,jm
        qrhs(i,j,k,2)=qrhs(i,j,k,2)+dfj(j,1)/dy
        qrhs(i,j,k,3)=qrhs(i,j,k,3)+dfj(j,2)/dy
        qrhs(i,j,k,4)=qrhs(i,j,k,4)+dfj(j,3)/dy
        qrhs(i,j,k,5)=qrhs(i,j,k,5)+dfj(j,4)/dy
      enddo
      !
    enddo
    enddo
    !
    !call progress_bar(3,4,'  ** diffusion terms ',10)
    !
    !$acc parallel loop gang collapse(2) present(sigma,qflux,qrhs) private(fk,dfk)
    do j=0,jm
    do i=0,im
      !
      do k=-hm,km+hm
        fk(k,1)=sigma(i,j,k,3)
        fk(k,2)=sigma(i,j,k,5)
        fk(k,3)=sigma(i,j,k,6)
        fk(k,4)=qflux(i,j,k,3)
      enddo
      !
      call diff6ec(fk(:,1),km,hm,dfk(:,1))
      call diff6ec(fk(:,2),km,hm,dfk(:,2))
      call diff6ec(fk(:,3),km,hm,dfk(:,3))
      call diff6ec(fk(:,4),km,hm,dfk(:,4))
      !
      do k=0,km
        qrhs(i,j,k,2)=qrhs(i,j,k,2)+dfk(k,1)/dz
        qrhs(i,j,k,3)=qrhs(i,j,k,3)+dfk(k,2)/dz
        qrhs(i,j,k,4)=qrhs(i,j,k,4)+dfk(k,3)/dz
        qrhs(i,j,k,5)=qrhs(i,j,k,5)+dfk(k,4)/dz
      enddo
      !
    enddo
    enddo
    !
    !$acc end data
    !
    if(present(comptime)) then
      call cpu_time(tfinish)
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

    real(8) :: phi(-hm:im+hm),fphi(0:im),phj(-hm:jm+hm),fphj(0:jm), &
               phk(-hm:km+hm),fphk(0:km)
    !
    integer :: i,j,k,n
    !
    if(present(comptime)) then
        call cpu_time(tstart)
    endif
    !
    !$acc data create(phi,fphi,phj,fphj,phk,fphk) 
    !
    !call progress_bar(0,3,'  ** spatial filter ',10)
    !
    !
    !$acc parallel loop gang collapse(2) present(q) private(phi,fphi)
    !
    do k=0,km
    do j=0,jm
      !
      do n=1,numq
        !
        do i=-hm,im+hm
          phi(i)=q(i,j,k,n)
        enddo
        !
        call filter10ec(phi,im,hm,fphi)
        !
        do i=0,im
          q(i,j,k,n)=fphi(i)
        enddo
        !
      enddo
      !
    enddo
    enddo
    !
    !call progress_bar(1,3,'  ** spatial filter ',10)
    !
    !$acc parallel loop gang collapse(2) present(q) private(phj,fphj) 
    do k=0,km
    do i=0,im
      !
      do n=1,numq
        !
        do j=-hm,jm+hm
          phj(j)=q(i,j,k,n)
        enddo
        !
        call filter10ec(phj,jm,hm,fphj)
        !
        do j=0,jm
          q(i,j,k,n)=fphj(j)
        enddo
        !
      enddo
      !
    enddo
    enddo
    !
    !call progress_bar(2,3,'  ** spatial filter ',10)
    !
    !$acc parallel loop gang collapse(2) present(q) private(phk,fphk) 
    do j=0,jm
    do i=0,im
      !
      do n=1,numq
        !
        do k=-hm,km+hm
          phk(k)=q(i,j,k,n)
        enddo
        !
        call filter10ec(phk,km,hm,fphk)
        !
        do k=0,km
          q(i,j,k,n)=fphk(k)
        enddo
        !
      enddo
      !
    enddo
    enddo
    !
    !$acc end data
    !
    if(present(comptime)) then
      call cpu_time(tfinish)
      !
      comptime=comptime+tfinish-tstart
    endif
    !
    !call progress_bar(3,3,'  ** spatial filter ',10)
    !
  end subroutine filterq
  !
  subroutine stacal
    !
    use comvardef, only: im,jm,km,rho,vel,dvel,nstep,time
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
    !
    !$acc parallel loop collapse(3) reduction(+:tke,rhom,enst) private(var1)
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
    !
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
  !
  subroutine validation(comptime)
    use comvardef, only: im,jm,km,hm,nstep, &
                         dvel, dtmp, qrhs, q, rho, vel, prs, tmp
    ! 
    real,intent(inout),optional :: comptime
    real :: tstart,tfinish
    ! 
    if(present(comptime)) then
      call cpu_time(tstart)
    endif
    ! 
    ! gradcal
    write(99, *) nstep, "dvel", "(1)", dvel(:,0,0,1,1)
    write(99, *) nstep, "dvel", "(2)", dvel(im,0,:,2,2)
    write(99, *) nstep, "dtmp", "(1)", dtmp(:,0,0,1)
    write(99, *) nstep, "dtmp", "(2)", dtmp(im,:,0,2)
    ! 
    ! convection & diffusion
    write(98, *) nstep, "qrhs", "(i)", qrhs(:,0,0,1)
    write(98, *) nstep, "qrhs", "(i)", qrhs(0,:,im,4)
    ! 
    ! filterq
    write(96, *) nstep, "q", "(i)", q(:,0,im,1)
    write(96, *) nstep, "q", "(i)", q(im,:,km,1)
    ! 
    ! q2fvar, bchomo
    write(94, *) nstep, "rho", rho(:,:,km)
    write(94, *) nstep, "vel", vel(im,0,:,2)
    write(94, *) nstep, "prs", prs(:,:,km)
    write(94, *) nstep, "tmp", tmp(:,jm,:)
    ! 
    if(present(comptime)) then
      call cpu_time(tfinish)
      !
      comptime=comptime+tfinish-tstart
    endif
  end subroutine validation
  !
  subroutine rk3(comptime)
    !
    use comvardef, only: im,jm,km,numq,num1d3,num2d3,q,qrhs,deltat,    &
                         rho,vel,tmp,prs,qsave,nstep,ctime,performval1
    use dataoper,  only: q2fvar
    !
    real,intent(inout),optional :: comptime
    
    ! local data
    logical,save :: firstcall = .true.
    real(8),save :: rkcoe(3,3)
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
      firstcall=.false.
      !
    endif
    !
    do rkstep=1,3
      !
      !$acc parallel loop collapse(4) present(qrhs)
      do m=1,numq
      do k=0,km
      do j=0,jm
      do i=0,im
        qrhs(i,j,k,m)=0.d0
      enddo
      enddo
      enddo
      enddo
      !
      call bchomo
      
      call rhscal(ctime(7))
      !
      if(rkstep==1) then
        !
        if(mod(nstep,10)==0) call stacal
        !
        !$acc parallel loop collapse(4) present(q,qsave)
        do m=1,numq
        do k=0,km
        do j=0,jm
        do i=0,im
          qsave(i,j,k,m)=q(i,j,k,m)
        enddo
        enddo
        enddo
        enddo
        !
      endif
      !
      !$acc parallel loop collapse(4) present(q,qsave,qrhs)
      do m=1,numq
      do k=0,km
      do j=0,jm
      do i=0,im
        !
        q(i,j,k,m)=rkcoe(1,rkstep)*qsave(i,j,k,m)+   &
                   rkcoe(2,rkstep)*q(i,j,k,m)    +   &
                   rkcoe(3,rkstep)*qrhs(i,j,k,m)*deltat
        !
      enddo
      enddo
      enddo
      enddo
      !
      call bchomovec(q)
      !
      call filterq(ctime(6))
      !
      !$acc parallel loop collapse(3) present(q,rho,vel,prs,tmp)
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
    if (performval1) then
      if(mod(nstep,10)==0) call validation()
    end if
    !
  end subroutine rk3
  !
  subroutine mainloop
    !
    use comvardef, only: time,nstep,deltat,ctime,performval1
    use comvardef, only: qrhs,q,qsave,rho,vel,tmp,prs,dvel,dtmp,sigma,qflux
    !
    !$acc data copy(qrhs,q,qsave,rho,vel,tmp,prs,dvel,dtmp,sigma,qflux )
    !
    if(performval1) then
      open(99, file = "gradcal_oa.txt")
      open(98, file = "convection_oa.txt")
      open(97, file = "diffusion_oa.txt")
      open(96, file = "filterq_oa.txt")
      open(95, file = "q2fvar_oa.txt")
      open(94, file = "bchomo_oa.txt")
    end if
    !
    do while(nstep<101)
      !
      call rk3(ctime(2))
      !
      nstep=nstep+1
      time =time + deltat
      !
      if(performval1 .and. mod(nstep,10)==0) then
        print*,nstep,time
      end if
      !
    enddo
    !
    !$acc end data
    !
    if(performval1) then
      close(99)
      close(98)
      close(97)
      close(96)
      close(95)
      close(94)
    end if
    !
    !
  end subroutine mainloop
  !
end module solver
!
module readwrite
  !
  implicit none
  !
  contains
  !
  subroutine timereport
    !
    use comvardef, only: alloarray,ctime
    !
    real :: dtime(18)
    logical :: lfex
    character(len=64) :: char1,char2,char3
    !
    inquire(file='time_report_oacc.txt',exist=lfex)
    !
    if(lfex) then
      open(15,file='time_report_oacc.txt')
      read(15,'(/)')
      read(15,*)char1,char1,dtime(1)
      print*,dtime(1)
      read(15,*)char1,char1,dtime(2)
      read(15,*)char1,char1,dtime(7)
      read(15,*)char1,char1,dtime(3)
      read(15,*)char1,char1,dtime(4)
      read(15,*)char1,char1,dtime(5)
      read(15,*)char1,char1,dtime(6)
      read(15,*)char1,char1,dtime(8)
      close(15)
      print*,' >> time_report_oacc.txt from a previous calculation'
    endif
    !
    open(14,file='time_report_oacc.txt')
    write(14,'(A)')'------------------cpu time cost------------------'
    write(14,'(A)')'-------------------------------------------------'
    write(14,'(A,2(F13.3,A))')'   total_time |',ctime(1),' s   |',dtime(1),' s'
    write(14,'(A,2(F13.3,A))')'          rk3 |',ctime(2),' s   |',dtime(2),' s'
    write(14,'(A,2(F13.3,A))')'       rhscal |',ctime(7),' s   |',dtime(7),' s'
    write(14,'(A,2(F13.3,A))')'      gradcal |',ctime(3),' s   |',dtime(3),' s'
    write(14,'(A,2(F13.3,A))')'   convection |',ctime(4),' s   |',dtime(4),' s'
    write(14,'(A,2(F13.3,A))')'    diffusion |',ctime(5),' s   |',dtime(5),' s'
    write(14,'(A,2(F13.3,A))')'      filterq |',ctime(6),' s   |',dtime(6),' s'
    write(14,'(A,2(F13.3,A))')'      diff6ec |',ctime(8),' s   |',dtime(8),' s'
    write(14,'(A)')'-------------------------------------------------'
    close(14)
    print*,' << time_report_oacc.txt'
    !
  end subroutine timereport
  !
end module readwrite
!
program boxsolver
    !
    use comvardef, only: alloarray,ctime
    use dataoper, only: tgvinit
    use solver,   only: mainloop
    use readwrite,only: timereport
    !
    implicit none
    !
    real :: tstart,tfinish
    !
    print*,' ** job has started **'
    !
    call cpu_time(tstart)
    !
    call alloarray
    !
    call tgvinit
    !
    call mainloop
    !
    call cpu_time(tfinish)
    !
    ctime(1)=tfinish-tstart
    !
    call timereport
    !
    print*,' ** job has been completed **'
    !
end program boxsolver