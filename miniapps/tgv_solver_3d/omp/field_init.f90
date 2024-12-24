module field_init
   
  implicit none

  contains
  
  subroutine flowinit
    
    use comvardef, only: flowtype

    if(trim(flowtype)=='tgv') then
      call tgvinit
    elseif(trim(flowtype)=='2dvortex') then
      call vortini
    else
      stop ' !! flowtype not defined !! @flowinit'
    endif

  end subroutine flowinit

  subroutine tgvinit

    !$ use omp_lib
    use constdef
    use comvardef, only: im,jm,km,dx,dy,dz,time,nstep,pinf, &
                         x,qrhs,q,vel,rho,prs,tmp,ctime,file_number
    use fluids, only: var2q,thermal_scar
    use xdmwrite, only: xdmfwriter
    use tecio, only: tecbin

    integer :: i,j,k
    real(rtype) :: var1
    character(len=4) :: file_name
    
    
    dx=2._rtype*pi/real(im,rtype)
    dy=2._rtype*pi/real(jm,rtype)
    dz=2._rtype*pi/real(km,rtype)
    
    var1=0._rtype

    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k)
    !$OMP DO
    do k=0,km
    do j=0,jm
    do i=0,im
      x(i,j,k,1)  =2._rtype*pi/real(im,rtype)*real(i,rtype)
      x(i,j,k,2)  =2._rtype*pi/real(jm,rtype)*real(j,rtype)
      x(i,j,k,3)  =2._rtype*pi/real(km,rtype)*real(k,rtype)
      
      rho(i,j,k)  =1._rtype
      vel(i,j,k,1)= sin(x(i,j,k,1))*cos(x(i,j,k,2))*cos(x(i,j,k,3))
      vel(i,j,k,2)=-cos(x(i,j,k,1))*sin(x(i,j,k,2))*cos(x(i,j,k,3))
      vel(i,j,k,3)=0._rtype
      prs(i,j,k)  =pinf+1._rtype/16._rtype*(cos(2._rtype*x(i,j,k,1))+cos(2._rtype*x(i,j,k,2)))*(cos(2._rtype*x(i,j,k,3))+2._rtype)
      
      tmp(i,j,k)  =thermal_scar(density=rho(i,j,k),pressure=prs(i,j,k))
      
      q(i,j,k,:)=var2q(density=rho(i,j,k),velocity=vel(i,j,k,:), pressure=prs(i,j,k))
      
    enddo
    enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    
    qrhs=0._rtype
    
    nstep=0
    time =0._rtype
    
    ctime=0.0
    
    file_number=0

    write(file_name,'(i4.4)')file_number

    print*,' ** flow field initilised'
    
    call tecbin(filename='flow'//file_name//'.plt',var1=x(0:im,0:jm,0:km,1),  var1name='x', &
                                        var2=x(0:im,0:jm,0:km,2),  var2name='y', &
                                        var3=x(0:im,0:jm,0:km,3),  var3name='z', &
                                        var4=vel(0:im,0:jm,0:km,1),var4name='u', &
                                        var5=vel(0:im,0:jm,0:km,2),var5name='v', &
                                        var6=prs(0:im,0:jm,0:km),  var6name='p')
    ! call xdmfwriter(dir='./',filename='flow000',deltax=x(1,0,0,1),  &
    !                                    var1=vel(0:im,0:jm,0:km,1),   var1name='u',  &
    !                                    var2=vel(0:im,0:jm,0:km,2),   var2name='v',  &
    !                                    var3=prs(0:im,0:jm,0:km),     var3name='p',  &
    !                                      im=im,jm=jm,km=km)

  end subroutine tgvinit

  ! Pirozzoli, S., Stabilized non-dissipative approximations of Euler equations in generalized curvilinear coordinates. Journal of Computational Physics, 2011. 230(8): p. 2997-3014.
  subroutine vortini
    !
    use constdef
    use comvardef, only: im,jm,km,dx,dy,dz,time,nstep,pinf,x,qrhs, &
                         q,vel,rho,prs,tmp,ctime,file_number,mach,const8,gamma
    use fluids, only: var2q,thermal_scar
    use solver, only: gradcal
    use xdmwrite, only: xdmfwriter
    use tecio, only: tecbin
    use io, only: write_field
    !
    ! local data
    integer :: i,j,k
    real(rtype) :: xc,yc,radi2,rvor,mv
    character(len=4) :: file_name
    !
    xc=10._rtype
    yc=5._rtype
    
    rvor=20._rtype/24._rtype
    
    mv=mach

    k=0
    do j=0,jm
    do i=0,im
      x(i,j,k,1)  =20._rtype/real(im,rtype)*real(i,rtype)
      x(i,j,k,2)  =10._rtype/real(jm,rtype)*real(j,rtype)
      x(i,j,k,3)  =0._rtype

      radi2=((x(i,j,k,1)-xc)**2+(x(i,j,k,2)-yc)**2)/rvor/rvor

      tmp(i,j,k)  =1._rtype

      vel(i,j,k,1)=1._rtype-mv/mach*(x(i,j,k,2)-yc)/rvor*exp(0.5_rtype*(1._rtype-radi2))
      vel(i,j,k,2)=         mv/mach*(x(i,j,k,1)-xc)/rvor*exp(0.5_rtype*(1._rtype-radi2))
      vel(i,j,k,3)=0._rtype

      prs(i,j,k)  =pinf*(1._rtype-0.5_rtype*(gamma-1._rtype)*mv*mv*exp(1._rtype-radi2))**const8

      rho(i,j,k)  =thermal_scar(temperature=tmp(i,j,k),pressure=prs(i,j,k))
      
      q(i,j,k,:)=var2q(density=rho(i,j,k),velocity=vel(i,j,k,:), pressure=prs(i,j,k))
      !
    enddo
    enddo

    dx=x(1,0,0,1)-x(0,0,0,1)
    dy=x(0,1,0,2)-x(0,0,0,2)
    dz=0._rtype
    !
    qrhs=0._rtype
    
    nstep=0
    time =0._rtype
    
    ctime=0.0

    file_number=-1

    call gradcal

    call write_field

    ! write(file_name,'(i4.4)')file_number

    ! call tecbin(filename='flow'//file_name//'.plt',var1=x(0:im,0:jm,0:km,1),  var1name='x', &
    !                                    var2=x(0:im,0:jm,0:km,2),  var2name='y', &
    !                                    var3=x(0:im,0:jm,0:km,3),  var3name='z', &
    !                                    var4=vel(0:im,0:jm,0:km,1),var4name='u', &
    !                                    var5=vel(0:im,0:jm,0:km,2),var5name='v', &
    !                                    var6=prs(0:im,0:jm,0:km),  var6name='p')

    write(*,'(A,I1,A)')'  ** 2-D vortical field initialised.'
    !
  end subroutine vortini

end module field_init