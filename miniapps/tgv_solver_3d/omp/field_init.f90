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
    use fluidynamcs, only: var2q,thermal_scar
    use xdmwrite, only: xdmfwriter
    use tecio, only: tecbin

    integer :: i,j,k
    real(8) :: var1
    character(len=4) :: file_name
    
    
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

  subroutine vortini
    !
    use constdef
    use comvardef, only: im,jm,km,dx,dy,dz,time,nstep,pinf,x,qrhs, &
                         q,vel,rho,prs,tmp,ctime,file_number
    use fluidynamcs, only: var2q,thermal_scar
    use xdmwrite, only: xdmfwriter
    use tecio, only: tecbin
    !
    ! local data
    integer :: i,j,k
    real(8) :: xc,yc,radi2,rvor,cvor,var1
    character(len=4) :: file_name
    !
    xc=5.d0
    yc=5.d0
    rvor=0.7d0 
    cvor=0.1d0*rvor
    !
    k=0
    do j=0,jm
    do i=0,im
      x(i,j,k,1)  =10.d0/dble(im)*dble(i)
      x(i,j,k,2)  =10.d0/dble(jm)*dble(j)
      x(i,j,k,3)  =0.d0

      radi2=((x(i,j,k,1)-xc)**2+(x(i,j,k,2)-yc)**2)/rvor/rvor
      var1=cvor/rvor/rvor*exp(-0.5d0*radi2)
      !
      rho(i,j,k)  =1.d0
      vel(i,j,k,1)=1.d0-var1*(x(i,j,k,2)-yc)
      vel(i,j,k,2)=     var1*(x(i,j,k,1)-xc)
      vel(i,j,k,3)=0.d0

      prs(i,j,k)  =pinf-0.5d0*1.d0*cvor**2/rvor**2*exp(-radi2)
      !
      tmp(i,j,k)  =thermal_scar(density=rho(i,j,k),pressure=prs(i,j,k))
      
      q(i,j,k,:)=var2q(density=rho(i,j,k),velocity=vel(i,j,k,:), pressure=prs(i,j,k))
      !
    enddo
    enddo

    dx=x(1,0,0,1)-x(0,0,0,1)
    dy=x(0,1,0,2)-x(0,0,0,2)
    dz=0.d0
    !
    qrhs=0.d0
    
    nstep=0
    time =0.d0
    
    ctime=0.0

    file_number=0

    write(file_name,'(i4.4)')file_number

    call tecbin(filename='flow'//file_name//'.plt',var1=x(0:im,0:jm,0:km,1),  var1name='x', &
                                       var2=x(0:im,0:jm,0:km,2),  var2name='y', &
                                       var3=x(0:im,0:jm,0:km,3),  var3name='z', &
                                       var4=vel(0:im,0:jm,0:km,1),var4name='u', &
                                       var5=vel(0:im,0:jm,0:km,2),var5name='v', &
                                       var6=prs(0:im,0:jm,0:km),  var6name='p')

    write(*,'(A,I1,A)')'  ** 2-D vortical field initialised.'
    !
  end subroutine vortini

end module field_init