module field_init
   
  implicit none

  contains
  
  subroutine flowinit(flowtype)
    
    character(len=*),intent(in) :: flowtype

    if(flowtype=='tgv') then
      call tgvinit
    else
      stop ' !! flowtype not defined !! @flowinit'
    endif

  end subroutine flowinit

  subroutine tgvinit

    !$ use omp_lib
    use constdef
    use comvardef, only: im,jm,km,gamma,mach,reynolds,prandtl,dx,dy,dz, &
                         const1,const2,const3,const4,const5,const6,    &
                         const7,time,nstep,deltat,ref_t,               &
                         x,qrhs,q,vel,rho,prs,tmp,ctime,mthread
    use  fluidynamcs, only: var2q,thermal_scar
    use xdmwrite, only: xdmfwriter
    use tecio, only: tecbin

    integer :: i,j,k
    real(8) :: var1,pinf
    
    const1=1.d0/(gamma*(gamma-1.d0)*mach**2)
    const2=gamma*mach**2
    const3=(gamma-1.d0)/3.d0*prandtl*(mach**2)
    const4=(gamma-1.d0)*mach**2*reynolds*prandtl
    const5=(gamma-1.d0)*mach**2
    const6=1.d0/(gamma-1.d0)
    const7=(gamma-1.d0)*mach**2*Reynolds*prandtl
    
    pinf=1.d0/const2
    
    dx=2.d0*pi/dble(im)
    dy=2.d0*pi/dble(jm)
    dz=2.d0*pi/dble(km)
    
    var1=0.d0

    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k)
    !$OMP DO
    do k=0,km
    do j=0,jm
    do i=0,im
      x(i,j,k,1)  =2.d0*pi/dble(im)*dble(i)
      x(i,j,k,2)  =2.d0*pi/dble(jm)*dble(j)
      x(i,j,k,3)  =2.d0*pi/dble(km)*dble(k)
      
      rho(i,j,k)  =1.d0
      vel(i,j,k,1)= sin(x(i,j,k,1))*cos(x(i,j,k,2))*cos(x(i,j,k,3))
      vel(i,j,k,2)=-cos(x(i,j,k,1))*sin(x(i,j,k,2))*cos(x(i,j,k,3))
      vel(i,j,k,3)=0.d0
      prs(i,j,k)  =pinf+1.d0/16.d0*(cos(2.d0*x(i,j,k,1))+cos(2.d0*x(i,j,k,2)))*(cos(2.d0*x(i,j,k,3))+2.d0)
      
      tmp(i,j,k)  =thermal_scar(density=rho(i,j,k),pressure=prs(i,j,k))
      
      q(i,j,k,:)=var2q(density=rho(i,j,k),velocity=vel(i,j,k,:), pressure=prs(i,j,k))
      
    enddo
    enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    
    qrhs=0.d0
    
    nstep=0
    time =0.d0
    
    ctime=0.0
    
    print*,' ** flow field initilised'
    
     call tecbin(filename='flow000.plt',var1=x(0:im,0:jm,0:km,1),  var1name='x', &
                                        var2=x(0:im,0:jm,0:km,2),  var2name='y', &
                                        var3=x(0:im,0:jm,0:km,3),  var3name='z', &
                                        var4=vel(0:im,0:jm,0:km,1),var4name='u', &
                                        var5=vel(0:im,0:jm,0:km,2),var5name='v', &
                                        var6=prs(0:im,0:jm,0:km),  var6name='p')
    call xdmfwriter(dir='./',filename='flow000',deltax=x(1,0,0,1),  &
                                       var1=vel(0:im,0:jm,0:km,1),   var1name='u',  &
                                       var2=vel(0:im,0:jm,0:km,2),   var2name='v',  &
                                       var3=prs(0:im,0:jm,0:km),     var3name='p',  &
                                         im=im,jm=jm,km=km)

  end subroutine tgvinit

end module field_init