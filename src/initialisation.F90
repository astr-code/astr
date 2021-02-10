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
  use parallel,only: lio,mpistop,mpirank,mpirankname
  use commvar, only: im,jm,km,uinf,vinf,pinf,roinf
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
    case default
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
      spc(i,j,k,1)=exp(-0.5d0*radi2)
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
end module initialisation
!+---------------------------------------------------------------------+
!| The end of the module initialisation.                               |
!+---------------------------------------------------------------------+