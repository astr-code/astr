!+---------------------------------------------------------------------+
!| This module contains subroutines of initialising flow field, should |
!| be highly user defined.                                             |
!| ==============                                                      |
!| CHANGE RECORD                                                       |
!| -------------                                                       |
!| 01-02-2021: Created by J. Fang @ STFC Daresbury Laboratory          |
!+---------------------------------------------------------------------+
module initialisation
  !
  use constdef
  use parallel,only: lio,mpistop,mpirank,mpirankname,bcast
  use commvar, only: im,jm,km,uinf,vinf,pinf,roinf,ndims,num_species
  use tecio
  !
  implicit none
  !
  contains
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is the entrance of flow initialisation.           |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 07-02-2021  | Created by J. Fang @ Warrington                     |
  !+-------------------------------------------------------------------+
  subroutine flowinit
    !
    use commvar,  only: flowtype,nstep,time,filenumb
    use readwrite,only: readcont
    !
    select case(trim(flowtype))
    case('2dvort')
      call vortini
    case('channel')
      call chanini
    case('tgv')
      call tgvini
    case default
      print*,trim(flowtype)
      stop ' !! flowtype not defined @ flowinit'
    end select
    !
    nstep=0
    time=0.d0
    !
    filenumb=0
    !
    call readcont
    !
    if(lio) print*,' ** flowfield initialised.'
    !
  end subroutine flowinit
  !+-------------------------------------------------------------------+
  !| The end of the subroutine flowinit.                               |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to generate an initial field for the      |
  !| 2D vortex transport problem.                                      |
  !| ref: Visbal & Gaitonde, On the Use of Higher-Order Finite-        |
  !|      Difference Schemes on Curvilinear and Deforming Meshes.      |
  !|     Journal of Computational Physics, 2002, 181: 155–185.         |                                      |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 13-Jul-2020: Created by J. Fang @ STFC Daresbury Laboratory       |
  !+-------------------------------------------------------------------+
  subroutine vortini
    !
    use commarray,only: x,vel,rho,prs,spc,tmp,q
    use fludyna,  only: thermal,fvar2q
    !
    ! local data
    integer :: i,j,k
    real(8) :: xc,yc,radi2,rvor,cvor,var1
    !
    xc=10.d0
    yc=5.d0
    rvor=0.7d0 
    cvor=0.1d0*rvor
    !
    do k=0,km
    do j=0,jm
    do i=0,im
      radi2=((x(i,j,k,1)-xc)**2+(x(i,j,k,2)-yc)**2)/rvor/rvor
      var1=cvor/rvor/rvor*exp(-0.5d0*radi2)
      !
      rho(i,j,k)  =roinf
      vel(i,j,k,1)=uinf-var1*(x(i,j,k,2)-yc)
      vel(i,j,k,2)=vinf+var1*(x(i,j,k,1)-xc)
      vel(i,j,k,3)=0.d0
      prs(i,j,k)  =pinf-0.5d0*roinf*cvor**2/rvor**2*exp(-radi2)
      !
      tmp(i,j,k)=thermal(density=rho(i,j,k),pressure=prs(i,j,k))
      !
      if(num_species>1) spc(i,j,k,1)=exp(-0.5d0*radi2)
    enddo
    enddo
    enddo
    !
    call fvar2q(          q=  q(0:im,0:jm,0:km,:),                     &
                    density=rho(0:im,0:jm,0:km),                       &
                   velocity=vel(0:im,0:jm,0:km,:),                     &
                   pressure=prs(0:im,0:jm,0:km),                       &
                    species=spc(0:im,0:jm,0:km,:)                      )
    !
    ! call tecbin('testout/tecinit'//mpirankname//'.plt',                &
    !                                   x(0:im,0:jm,0:km,1),'x',         &
    !                                   x(0:im,0:jm,0:km,2),'y',         &
    !                                   q(0:im,0:jm,0:km,1),'q1',        &
    !                                   q(0:im,0:jm,0:km,2),'q2',        &
    !                                   q(0:im,0:jm,0:km,3),'q3',        &
    !                                   q(0:im,0:jm,0:km,5),'q5',        &
    !                                   q(0:im,0:jm,0:km,6),'q6' )
    !
    if(lio)  write(*,'(A,I1,A)')'  ** 2-D vortical field initialised.'
    !
  end subroutine vortini
  !+-------------------------------------------------------------------+
  !| The end of the subroutine vortini.                                |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to generate 3d TGV initial flow           |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 12-02-2021: Created by J. Fang @ STFC Daresbury Laboratory       |
  !+-------------------------------------------------------------------+
  subroutine tgvini
    !
    use commarray,only: x,vel,rho,prs,spc,tmp,q
    use fludyna,  only: thermal,fvar2q
    !
    ! local data
    integer :: i,j,k
    !
    do k=0,km
    do j=0,jm
    do i=0,im
      rho(i,j,k)  =roinf
      vel(i,j,k,1)= sin(x(i,j,k,1))*cos(x(i,j,k,2))*cos(x(i,j,k,3))
      vel(i,j,k,2)=-cos(x(i,j,k,1))*sin(x(i,j,k,2))*cos(x(i,j,k,3))
      vel(i,j,k,3)=0.d0
      prs(i,j,k)  =pinf+roinf/16.d0*( cos(2.d0*x(i,j,k,1)) +           &
                                      cos(2.d0*x(i,j,k,2)) )*          &
                                             (cos(2.d0*x(i,j,k,3))+2.d0)
      !
      tmp(i,j,k)=thermal(density=rho(i,j,k),pressure=prs(i,j,k))
      !
      if(num_species>1) spc(i,j,k,1)=1.d0
    enddo
    enddo
    enddo
    !
    call fvar2q(          q=  q(0:im,0:jm,0:km,:),                     &
                    density=rho(0:im,0:jm,0:km),                       &
                   velocity=vel(0:im,0:jm,0:km,:),                     &
                   pressure=prs(0:im,0:jm,0:km),                       &
                    species=spc(0:im,0:jm,0:km,:)                      )
    !
    ! call tecbin('testout/tecinit'//mpirankname//'.plt',                &
    !                                   x(0:im,0:jm,0:km,1),'x',         &
    !                                   x(0:im,0:jm,0:km,2),'y',         &
    !                                   x(0:im,0:jm,0:km,3),'z',         &
    !                                rho(0:im,0:jm,0:km)  ,'ro',         &
    !                                vel(0:im,0:jm,0:km,1),'u',          &
    !                                vel(0:im,0:jm,0:km,2),'v',          &
    !                                prs(0:im,0:jm,0:km)  ,'p',          &
    !                                tmp(0:im,0:jm,0:km)  ,'t' )
    !
    if(lio)  write(*,'(A,I1,A)')'  ** 2-D vortical field initialised.'
    !
  end subroutine tgvini
  !+-------------------------------------------------------------------+
  !| The end of the subroutine tgvini.                                 |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to generate an initial field for the      |
  !| simulation of channel flow.                                       |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 27-May-2020: Created by J. Fang @ STFC Daresbury Laboratory       |
  !+-------------------------------------------------------------------+
  subroutine chanini
    !
    use commvar,  only: prandtl,mach,gamma,xmin,xmax,ymin,ymax,zmin,zmax
    use commarray,only: x,vel,rho,prs,spc,tmp,q
    use fludyna,  only: thermal,fvar2q
    !
    ! local data
    integer :: i,j,k,l
    real(8) :: theta,theter,fx,gz,zl,randomv(15)
    !
    theta=0.1d0
    !
    if(ndims==2) then
      !
      do k=0,km
      do j=0,jm
      do i=0,im
        rho(i,j,k)  =roinf
        vel(i,j,k,1)=1.5d0*(1.d0-(x(i,j,k,2)-1.d0)**2)
        vel(i,j,k,2)=0.d0
        vel(i,j,k,3)=0.d0
        tmp(i,j,k)=1.d0+(gamma-1.d0)*prandtl*mach**2/3.d0*             &
                                     1.5d0*(1.d0-((x(i,j,k,2)-1.d0)**4))
        !
        prs(i,j,k)=thermal(density=rho(i,j,k),temperature=tmp(i,j,k))
        !
        if(num_species>0) then
          spc(i,j,k,1)=(tanh((x(i,j,k,2)-1.d0)/theta)+1.d0)*0.5d0
        endif
        !
      enddo
      enddo
      enddo
      !
    elseif(ndims==3) then
      !
      if(mpirank==0) then
        do l=1,15
          call random_number(randomv(l))
        end do
      endif
      !
      call bcast(randomv)
      !
      do k=0,km
      do j=0,jm
      do i=0,im
        !
        theter=x(i,0,k,1)/(xmax-xmin)*2.d0*pi
        fx=4.d0*dsin(theter)*(1.d0-dcos(theter))*0.192450089729875d0
        !
        gz=0.d0
        do l=1,10
          !
          if(l==1) then
            zl=0.2d0/(1.d0-0.8d0**10)
          else
            zl=zl*0.8d0
          end if
          !
          gz=gz+zl*dsin(2.d0*pi*l*(x(i,0,k,3)/(zmax-zmin)+randomv(l)))
          !
        end do
        !
        rho(i,j,k)  =roinf
        vel(i,j,k,1)=1.5d0*(1.d0-(x(i,j,k,2)-1.d0)**2)*(1.d0+0.3d0*fx*gz)
        vel(i,j,k,2)=0.d0
        vel(i,j,k,3)=0.d0
        tmp(i,j,k)=1.d0+(gamma-1.d0)*prandtl*mach**2/3.d0*             &
                                     1.5d0*(1.d0-((x(i,j,k,2)-1.d0)**4))
        !
        prs(i,j,k)=thermal(density=rho(i,j,k),temperature=tmp(i,j,k))
        !
        if(num_species>0) then
          spc(i,j,k,1)=(tanh((x(i,j,k,2)-1.d0)/theta)+1.d0)*0.5d0
        endif
        !
      enddo
      enddo
      enddo
      !
      !
    endif
    !
    call fvar2q(          q=  q(0:im,0:jm,0:km,:),                     &
                    density=rho(0:im,0:jm,0:km),                       &
                   velocity=vel(0:im,0:jm,0:km,:),                     &
                   pressure=prs(0:im,0:jm,0:km),                       &
                    species=spc(0:im,0:jm,0:km,:)                      )
    !
    ! call tecbin('testout/tecinit'//mpirankname//'.plt',                &
    !                                   x(0:im,0:jm,0:km,1),'x',         &
    !                                   x(0:im,0:jm,0:km,2),'y',         &
    !                                rho(0:im,0:jm,0:km)  ,'ro',         &
    !                                vel(0:im,0:jm,0:km,1),'u',          &
    !                                vel(0:im,0:jm,0:km,2),'v',          &
    !                                prs(0:im,0:jm,0:km)  ,'p',          &
    !                                spc(0:im,0:jm,0:km,1),'Y1' )
    
    if(lio)  write(*,'(A,I1,A)')'  ** ',ndims,'-D channel flow initialised.'
    !
  end subroutine chanini
  !+-------------------------------------------------------------------+
  !| The end of the subroutine chanini.                                |
  !+-------------------------------------------------------------------+
  !!
end module initialisation
!+---------------------------------------------------------------------+
!| The end of the module initialisation.                               |
!+---------------------------------------------------------------------+