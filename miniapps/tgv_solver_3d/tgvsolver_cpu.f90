module comvardef
  !
  implicit none
  !
  integer :: im=128,jm=128,km=128,hm=5,numq=5

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
  real(8),allocatable,dimension(:,:,:) :: rho,prs,tmp
  real(8),allocatable,dimension(:,:,:,:) :: x,q,qrhs,vel,dtmp
  real(8),allocatable,dimension(:,:,:,:,:) :: dvel
  !
  contains
  !
  subroutine alloarray
    !
    allocate( x(0:im,0:jm,0:km,1:3))
    allocate(   q(-hm:im+hm,-hm:jm+hm,-hm:km+hm,1:numq))
    allocate(qrhs(0:im,0:jm,0:km,1:numq))
    !
    allocate( rho(-hm:im+hm,-hm:jm+hm,-hm:km+hm))
    allocate( vel(-hm:im+hm,-hm:jm+hm,-hm:km+hm,1:3))
    allocate( prs(-hm:im+hm,-hm:jm+hm,-hm:km+hm))
    allocate( tmp(-hm:im+hm,-hm:jm+hm,-hm:km+hm))
    allocate(dvel(0:im,0:jm,0:km,1:3,1:3))
    allocate(dtmp(0:im,0:jm,0:km,1:3))
    !
    print*,' ** common array allocated'
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
    real(8),intent(out) :: density
    real(8),intent(out),optional :: velocity(:),pressure,temperature
    !
    density   =q(1)
    !
    if(present(velocity) .or. present(pressure) .or. present(temperature)) then
      velocity(1)=q(2)/density
      velocity(2)=q(3)/density
      velocity(3)=q(4)/density
    endif
    !
    if(present(pressure) .or. present(temperature)) then
      pressure  =( q(5)-0.5d0*density*(velocity(1)**2+velocity(2)**2+  &
                                       velocity(3)**2) )/const6
    endif
    !
    if(present(temperature)) then
      temperature=thermal_scar(pressure=pressure,density=density)
    endif

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
    print*,' ** data initilised'
    !
  end subroutine tgvinit
  !
end module dataoper
!
module numerics
  !
  use comvardef, only : num1d60
  implicit none
  !
  contains
  !
  subroutine diff6ec(vin,dim,n,vout,comptime)
    !
    integer,intent(in) :: dim,n
    real(8),intent(in) :: vin(-n:dim+n)
    real(8) :: vout(0:dim)
    !
    real,intent(inout),optional :: comptime
    
    ! local data
    integer :: i
    real :: tstart,tfinish

    if(present(comptime)) then
        call cpu_time(tstart)
    endif
    !
    do i=0,dim
      vout(i)  =0.75d0 *(vin(i+1)-vin(i-1))- &
                0.15d0 *(vin(i+2)-vin(i-2))+ &
               num1d60 *(vin(i+3)-vin(i-3))
    enddo
    !
    if(present(comptime)) then
      call cpu_time(tfinish)
      !
      comptime=comptime+tfinish-tstart
    endif
    !
  end subroutine diff6ec
    !
  subroutine filter10ec(vin,dim,n,vout)
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
    use comvardef, only: qrhs,ctime
    !
    real,intent(inout),optional :: comptime
    
    real :: tstart,tfinish

    if(present(comptime)) then
        call cpu_time(tstart)
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
      call cpu_time(tfinish)
      !
      comptime=comptime+tfinish-tstart
    endif
  end subroutine rhscal
  !
  subroutine gradcal(comptime)
    !
    use comvardef,only: im,jm,km,hm,vel,dvel,tmp,dtmp,dx,dy,dz,ctime
    use numerics, only: diff6ec
    !
    integer :: i,j,k,n
    !
    real(8),allocatable,dimension(:,:) :: f,df
    real,intent(inout),optional :: comptime
    !
    real :: tstart,tfinish
    !
    !call progress_bar(0,3,'  ** temperature and velocity gradient ',10)
    !
    if(present(comptime)) then
        call cpu_time(tstart)
    endif
    !
    allocate(f(-hm:im+hm,4),df(0:im,4))
    do k=0,km
    do j=0,jm
      !
      do i=-hm,im+hm
        f(i,1)=vel(i,j,k,1)
        f(i,2)=vel(i,j,k,2)
        f(i,3)=vel(i,j,k,3)
        f(i,4)=tmp(i,j,k)
      enddo
      !
      call diff6ec(f(:,1),im,hm,df(:,1))
      call diff6ec(f(:,2),im,hm,df(:,2))
      call diff6ec(f(:,3),im,hm,df(:,3))
      call diff6ec(f(:,4),im,hm,df(:,4))
      !
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
    deallocate(f,df)
    !
    !call progress_bar(1,3,'  ** temperature and velocity gradient ',10)
    !
    allocate(f(-hm:jm+hm,4),df(0:jm,4))
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
      call diff6ec(f(:,1),jm,hm,df(:,1))
      call diff6ec(f(:,2),jm,hm,df(:,2))
      call diff6ec(f(:,3),jm,hm,df(:,3))
      call diff6ec(f(:,4),jm,hm,df(:,4))
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
    deallocate(f,df)
    !
    !call progress_bar(2,3,'  ** temperature and velocity gradient ',10)
    !
    allocate(f(-hm:km+hm,4),df(0:km,4))
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
      call diff6ec(f(:,1),km,hm,df(:,1))
      call diff6ec(f(:,2),km,hm,df(:,2))
      call diff6ec(f(:,3),km,hm,df(:,3))
      call diff6ec(f(:,4),km,hm,df(:,4))
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
    deallocate(f,df)
    !
    if(present(comptime)) then
      call cpu_time(tfinish)
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
    use comvardef, only: im,jm,km,hm,numq,rho,prs,tmp,vel,q,qrhs,dx,dy,dz,ctime
    use dataoper,  only: var2q
    use numerics,  only: diff6ec
    !
    real,intent(inout),optional :: comptime
    !
    real(8),allocatable,dimension(:,:) :: fcs,dfcs
    integer :: i,j,k,n
    !
    real :: tstart,tfinish
    
    if(present(comptime)) then
        call cpu_time(tstart)
    endif
    !
    !call progress_bar(0,3,'  ** convection terms ',10)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! calculating along i direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    allocate(fcs(-hm:im+hm,1:numq),dfcs(0:im,1:numq))
    !
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
      do n=1,numq
        !
        call diff6ec(fcs(:,n),im,hm,dfcs(:,n))
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
      !
      do n=1,numq
        !
        call diff6ec(fcs(:,n),jm,hm,dfcs(:,n))
        !
        do j=0,jm
          !
          qrhs(i,j,k,n)=qrhs(i,j,k,n)+dfcs(j,n)/dy
          !
        enddo
        !
      enddo
      !
    enddo
    enddo
    !
    deallocate(fcs,dfcs)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! end calculating along j direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !call progress_bar(2,3,'  ** convection terms ',10)
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! calculating along k direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    allocate(fcs(-hm:jm+hm,1:numq),dfcs(0:jm,1:numq))
    !
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
      !
      do n=1,numq
        !
        call diff6ec(fcs(:,n),km,hm,dfcs(:,n))
        !
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
    !
    deallocate(fcs,dfcs)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! end calculating along k direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
    use comvardef, only: im,jm,km,hm,vel,dvel,dtmp,qrhs,tmp,ctime
    use dataoper,  only: miucal
    use numerics,  only: diff6ec
    !
    real,intent(inout),optional :: comptime
    !
    real(8),allocatable,dimension(:,:,:,:),save :: sigma,qflux
    real(8),allocatable :: f(:,:),df(:,:)
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

    real(8),allocatable :: phi(:),fph(:)
    !
    integer :: i,j,k,n
    !
    if(present(comptime)) then
        call cpu_time(tstart)
    endif
    !
    !call progress_bar(0,3,'  ** spatial filter ',10)
    !
    allocate(phi(-hm:im+hm),fph(0:im))
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
    deallocate(phi,fph)
    !
    !call progress_bar(1,3,'  ** spatial filter ',10)
    !
    allocate(phi(-hm:jm+hm),fph(0:jm))
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
    deallocate(phi,fph)
    !
    !call progress_bar(2,3,'  ** spatial filter ',10)
    !
    allocate(phi(-hm:km+hm),fph(0:km))
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
    deallocate(phi,fph)
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
  !
  subroutine rk3(comptime)
    !
    use comvardef, only: im,jm,km,numq,num1d3,num2d3,q,qrhs,deltat,    &
                         rho,vel,tmp,prs,nstep,ctime
    use dataoper,  only: q2fvar
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
      !
      call bchomo
      
      call rhscal(ctime(7))
      !
      if(rkstep==1) then
        !
        if(mod(nstep,10)==0) call stacal
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
      call filterq(ctime(6))
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
  !
  subroutine mainloop
    !
    use comvardef, only: time,nstep,deltat,ctime
    !
    do while(nstep<100)
      !
      call rk3(ctime(2))
      !
      nstep=nstep+1
      time =time + deltat
      !
      print*,nstep,time
      !
    enddo
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
    open(14,file='time_report.txt')
    write(14,'(A)')'---------cpu time cost---------'
    write(14,'(A)')'-------------------------------'
    write(14,'(A,F13.3,A)')'   total time |',ctime(1),' s'
    write(14,'(A,F13.3,A)')'          rk3 |',ctime(2),' s'
    write(14,'(A,F13.3,A)')'       rhscal |',ctime(7),' s'
    write(14,'(A,F13.3,A)')'      gradcal |',ctime(3),' s'
    write(14,'(A,F13.3,A)')'   convection |',ctime(4),' s'
    write(14,'(A,F13.3,A)')'    diffusion |',ctime(5),' s'
    write(14,'(A,F13.3,A)')'      filterq |',ctime(6),' s'
    write(14,'(A,F13.3,A)')'      diff6ec |',ctime(8),' s'
    write(14,'(A)')'-------------------------------'
    close(14)
    print*,' << time_report.txt'
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
    print*,' ** the job is done **'
    !
end program boxsolver