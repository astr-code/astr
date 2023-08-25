!+---------------------------------------------------------------------+
!| This module contains user defined subroutines to interfere  program |
!+---------------------------------------------------------------------+
!| CHANGE RECORD                                                       |
!| -------------                                                       |
!| 18-08-2023  | Created by J. Fang                                    |
!+---------------------------------------------------------------------+
module userdefine
  !
  implicit none
  !
  real(8) :: hsource
  !
  contains
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to set flow environment, such as, incoming     |
  !| free stream variables.                                            |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 18-Aug-2023: created by Jian Fang @ Daresbury                     |
  !+-------------------------------------------------------------------+
  subroutine udf_setflowenv
    !
!     use commvar,  only: roinf,uinf,vinf,winf,pinf,tinf,spcinf,num_species
!     use fludyna,  only: thermal
!     !
! #ifdef COMB
!     use thermchem,only : tranco,spcindex,mixture,convertxiyi
!     use cantera 
!     !
!     real(8) :: specr(num_species)
!     ! 
!     specr(:)=0.d0
!     specr(spcindex('H2'))=0.0173
!     specr(spcindex('O2'))=0.2289
!     specr(spcindex('N2'))=1.d0-sum(specr)
!     !
!     ! pinf=5.d0*pinf
!     uinf=0.97d0
!     vinf=0.d0
!     winf=0.d0
!     tinf=300.d0
!     spcinf(:)=specr(:)
!     roinf=thermal(pressure=pinf,temperature=tinf,species=spcinf(:))
!     !
! #endif
    !
  end subroutine udf_setflowenv
  !+-------------------------------------------------------------------+
  !| The end of the subroutine udf_setflowenv.                         |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to initialise flowfield by a user              |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 25-May-2023: Created by Yifan Xu @ Peking University              |
  !| 18-Aug-2023: Rename and relocated by Jian Fang @ Daresbury        |
  !+-------------------------------------------------------------------+
  subroutine udf_flowinit
    !
!     use commvar,  only: im,jm,km,ndims,roinf,uinf,nondimen,xmax,pinf,  &
!                         ia,num_species
!     use commarray,only: x,vel,rho,prs,spc,tmp,q
!     use parallel, only: lio
!     use fludyna,  only: thermal
!     !
! #ifdef COMB
!     !
!     use thermchem,only : tranco,spcindex,mixture,convertxiyi
!     use cantera 
!     !
!     ! local data
!     integer :: i,j,k
!     real(8) ::  xc,yc,zc,tmpr,tmpp,xloc,xwid,specr(num_species),  &
!       specp(num_species),arg,prgvar,masflx,specx(num_species)
!     real(8) :: pthick
!     !
!     tmpr=300.d0
!     xloc=3.d0*xmax/4.d0
!     xwid=xmax/(12.d0*5.3d0*2.d0)
!     !
!     !reactants
!     specr(:)=0.d0
!     specr(spcindex('H2'))=0.0173
!     specr(spcindex('O2'))=0.2289
!     specr(spcindex('N2'))=1.d0-sum(specr)
!     !
!     !products
!     tmpp=1814.32d0
!     !
!     ! pthick=1.d-4
!     !
!     do k=0,km
!     do j=0,jm
!     do i=0,im
!       !
!       xc=x(i,j,k,1)
!       !
!       !prgvar=0.5d0*(1.d0+tanh(10.d0*(xc-xloc)/xloc))
!       ! if(xc-xloc<xwid*0.5d0*1.2d0) then 
!       !   prgvar=0.d0
!       !   if(xc-xloc>xwid*0.5d0) &
!       !   prgvar=1.d0-(xc-xloc-(xwid*0.5d0))/(xwid*0.5d0*0.2d0)
!       ! else
!       !   prgvar=1.d0
!       ! endif
!       !
!       prgvar=1.d0*exp(-0.5d0*((xc-xloc)/xwid)**2)
!       !
!       spc(i,j,k,:)=specr(:)
!       !
!       vel(i,j,k,1)=uinf
!       !
!       vel(i,j,k,2)=0.d0
!       vel(i,j,k,3)=0.d0
!       !
!       tmp(i,j,k)=tmpr+prgvar*(tmpp-tmpr)
!       !
!       prs(i,j,k)=pinf
!       !
!       rho(i,j,k)=thermal(pressure=prs(i,j,k),temperature=tmp(i,j,k), &
!                           species=spc(i,j,k,:))
!     enddo
!     enddo
!     enddo
!     !
!     !
!     if(lio)  write(*,'(A,I1,A)')'  ** HIT flame initialised.'
!     !
! #endif
    !
  end subroutine udf_flowinit
  !+-------------------------------------------------------------------+
  !| The end of the subroutine udf_flowinit.                           |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to generate grid.                              | 
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 23-Aug-2023: created by Jian Fang @ Daresbury                     |
  !+-------------------------------------------------------------------+
  subroutine udf_grid
  end subroutine udf_grid
  !+-------------------------------------------------------------------+
  !| The end of the subroutine udf_grid.                               |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to list something during a computation.        | 
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 18-Aug-2023: created by Jian Fang @ Daresbury                     |
  !+-------------------------------------------------------------------+
  subroutine udf_stalist
    !
    use constdef
    use commvar,  only : reynolds,lrestart,mach,ia,ja,ka,im,jm,km
    use commarray,only : vel,rho,tmp,dvel,q
    use fludyna,  only : miucal,sos
    use utility,  only : listinit,listwrite
    use parallel, only : psum,lio
    !
    integer :: i,j,k,ns
    real(8) :: s11,s12,s13,s22,s23,s33,div,miu,dissa
    real(8) :: rsamples,miudrho,dudx2,csavg,v2,cs,ufluc,ens,omegax,omegay,omegaz
    real(8) :: urms,taylorlength,kolmoglength,Retaylor,machrms,macht,rhoe
    !
    logical :: fex
    logical,save :: linit=.true.
    integer,save :: hand_fs
    !
    if(linit) then
      !
      if(lio) then
        call listinit(filename='fturbstats.dat',handle=hand_fs, &
                      firstline='nstep time urms enstrophy taylorlength kolmoglength Retaylor machrms macht Tavg hsource')
      endif
      !
      linit=.false.
      !
    endif
    !
    rsamples=dble(ia*ja*ka)
    !
    urms=0.d0
    machrms=0.d0
    !
    miudrho=0.d0
    dudx2=0.d0
    csavg=0.d0
    dissa=0.d0
    ens=0.d0
    do k=1,km
    do j=1,jm
    do i=1,im
      !
      miu=miucal(tmp(i,j,k))/Reynolds
      !
      s11=dvel(i,j,k,1,1)
      s12=0.5d0*(dvel(i,j,k,1,2)+dvel(i,j,k,2,1))
      s13=0.5d0*(dvel(i,j,k,1,3)+dvel(i,j,k,3,1))
      s22=dvel(i,j,k,2,2)
      s23=0.5d0*(dvel(i,j,k,2,3)+dvel(i,j,k,3,2))
      s33=dvel(i,j,k,3,3)
      !
      div=s11+s22+s33
      !
      omegax=dvel(i,j,k,3,2)-dvel(i,j,k,2,3)
      omegay=dvel(i,j,k,1,3)-dvel(i,j,k,3,1)
      omegaz=dvel(i,j,k,2,1)-dvel(i,j,k,1,2)
      !
      v2=vel(i,j,k,1)**2+vel(i,j,k,2)**2+vel(i,j,k,3)**2
      !
      cs=sos(tmp(i,j,k))
      !
      urms=urms+v2
      !
      dudx2=dudx2+s11**2+s22**2+s33**2
      !
      miudrho=miudrho+miu/rho(i,j,k)
      !
      csavg=csavg+cs
      !
      machrms=machrms+v2/(cs*cs)
      !
      rhoe=rhoe+tmp(i,j,k)
      !
      dissa=dissa+2.d0*miu*(s11**2+s22**2+s33**2+2.d0*(s12**2+s13**2+s23**2)-num1d3*div**2)
      !
      ens=ens+(omegax*omegax+omegay*omegay+omegaz*omegaz)
      !
    enddo
    enddo
    enddo
    urms  = sqrt(psum(urms)/rsamples)
    dudx2      = num1d3*psum(dudx2)/rsamples
    miudrho    = psum(miudrho)/rsamples
    csavg      = psum(csavg)/rsamples
    dissa      = psum(dissa)/rsamples
    !
    machrms=sqrt(psum(machrms)/rsamples)
    !
    rhoe=psum(rhoe)/rsamples
    !
    ens=ens/rsamples
    !
    ufluc=urms/sqrt(3.d0)
    !
    macht         = urms/csavg
    taylorlength  = ufluc/sqrt(dudx2)
    retaylor      = ufluc*taylorlength/miudrho
    kolmoglength  = sqrt(sqrt(miudrho**3/dissa))
    ! kolmogvelocity= sqrt(sqrt(dissipation*miudrho))
    ! kolmogtime    = sqrt(miudrho/dissipation)
    !
    if(lio) call listwrite(hand_fs,urms,ens,taylorlength,kolmoglength,Retaylor,machrms,macht,rhoe,hsource)
    !
  end subroutine udf_stalist
  !+-------------------------------------------------------------------+
  !| The end of the subroutine udf_stalist.                            |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to add vortical fluctuations to initial field  | 
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 18-Aug-2023: created by Jian Fang @ Daresbury                     |
  !+-------------------------------------------------------------------+
  subroutine addvortex(xc,yc,radius,amp)
    !
    use commvar,  only: im,jm,km,ndims,roinf,uinf
    use parallel, only: lio
    use commarray,only: x,vel,rho,prs,spc,tmp,q
    use fludyna,  only: thermal
    !
    ! local data
    real(8),intent(in) :: xc,yc,radius,amp
    !
    integer :: i,j,k
    real(8) :: var1,radi2,cvor
    !
    cvor=amp*uinf*radius
    !
    do k=0,km
    do j=0,jm
    do i=0,im
      radi2=((x(i,j,k,1)-xc)**2+(x(i,j,k,2)-yc)**2)/radius/radius
      var1=cvor/radius/radius*exp(-0.5d0*radi2)
      !
      vel(i,j,k,1)=vel(i,j,k,1)-var1*(x(i,j,k,2)-yc)
      if(ndims>=2) vel(i,j,k,2)=vel(i,j,k,2)+var1*(x(i,j,k,1)-xc)
      if(ndims==3) vel(i,j,k,3)=0.d0
      prs(i,j,k)  =prs(i,j,k)-0.5d0*roinf*cvor*cvor/radi2/radi2*exp(-radi2)
      !
      tmp(i,j,k)=thermal(density=rho(i,j,k),pressure=prs(i,j,k),species=spc(i,j,k,:))
      !
    enddo
    enddo
    enddo
    !
  end subroutine addvortex
  !+-------------------------------------------------------------------+
  !| The end of the subroutine addvortex.                              |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine add a source term to the rsd of the equation to   |
  !| hit flame.                                                        |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 13-06-2023: Created by Yifan Xu @ Peking University               |
  !+-------------------------------------------------------------------+
  subroutine udf_src
    !
    use constdef
    use commvar,  only : im,jm,km,ndims,deltat,ia,ja,ka
    use parallel, only : psum 
    use commarray,only : vel,qrhs,x,jacob
    ! !
    ! ! local data
    integer :: i,j,k,k1,k2
    real(8) :: A(3,3),B(3,3),force(3)
    real(8) :: FC,FR,var1,var2
    !
    FR=1.d0/16.d0
    FC=0.d0
    !
    var1=sqrt(num2d3*FC/deltat)
    var2=sqrt(num1d3*FR/deltat)
    !
    do j=1,3
    do i=1,3
      call random_number(A(i,j))
      call random_number(B(i,j))
    end do
    end do
    !
    A=sqrt(3.d0)*(2.d0*A-1.d0)
    B=sqrt(3.d0)*(2.d0*B-1.d0)
    !
    do j=1,3
    do i=1,3
      if(i==j) then
        A(i,j)=var1*A(i,j)
        B(i,j)=var1*B(i,j)
      else
        A(i,j)=var2*A(i,j)
        B(i,j)=var2*B(i,j)
      endif
    end do
    end do
    !
    hsource=0.d0
    !
    do k=0,km
    do j=0,jm
    do i=0,im
      force(1)=A(1,1)*sin(x(i,j,k,1))+A(1,2)*sin(x(i,j,k,2))+      &
               A(1,3)*sin(x(i,j,k,3))+B(1,1)*cos(x(i,j,k,1))+      &
               B(1,2)*cos(x(i,j,k,2))+B(1,3)*cos(x(i,j,k,3))
      !
      force(2)=A(2,1)*sin(x(i,j,k,1))+A(2,2)*sin(x(i,j,k,2))+      &
               A(2,3)*sin(x(i,j,k,3))+B(2,1)*cos(x(i,j,k,1))+      &
               B(2,2)*cos(x(i,j,k,2))+B(2,3)*cos(x(i,j,k,3))
      !                                                          
      force(3)=A(3,1)*sin(x(i,j,k,1))+A(3,2)*sin(x(i,j,k,2))+      &
               A(3,3)*sin(x(i,j,k,3))+B(3,1)*cos(x(i,j,k,1))+      &
               B(3,2)*cos(x(i,j,k,2))+B(3,3)*cos(x(i,j,k,3))
      !
      qrhs(i,j,k,2)=qrhs(i,j,k,2)+force(1)*jacob(i,j,k)
      qrhs(i,j,k,3)=qrhs(i,j,k,3)+force(2)*jacob(i,j,k)
      qrhs(i,j,k,4)=qrhs(i,j,k,4)+force(3)*jacob(i,j,k)
      !
      qrhs(i,j,k,5)=qrhs(i,j,k,5)+( force(1)*vel(i,j,k,1)+force(2)*vel(i,j,k,2)+   &
                                    force(3)*vel(i,j,k,3) )*jacob(i,j,k)
      !
      if(i.ne.0 .and. j.ne.0 .and. k.ne.0) then
        hsource=hsource+force(1)*vel(i,j,k,1)+force(2)*vel(i,j,k,2)+force(3)*vel(i,j,k,3)
      endif
      !
    end do
    end do
    end do
    !
    hsource=psum(hsource)/dble(ia*ja*ka)
    !
    do k=0,km
    do j=0,jm
    do i=0,im
      qrhs(i,j,k,5)=qrhs(i,j,k,5)-hsource*jacob(i,j,k)
    end do
    end do
    end do
    !
    ! real(8) :: dy,u1,u2,u3
    ! !
    ! if(ndims==2) then
    !   k1=0
    !   k2=0
    ! elseif(ndims==3) then
    !   k1=1
    !   k2=km
    ! else
    !   print*,' !! ndims=',ndims
    !   stop ' !! error @ massfluxchan !!'
    ! endif
    ! !
    ! do k=0,km
    ! do j=0,jm
    ! do i=0,im
    !   qrhs(i,j,k,2)=qrhs(i,j,k,2)+force(1)*jacob(i,j,k)
    !   qrhs(i,j,k,3)=qrhs(i,j,k,3)+force(2)*jacob(i,j,k)
    !   qrhs(i,j,k,4)=qrhs(i,j,k,4)+force(3)*jacob(i,j,k)
    !   qrhs(i,j,k,5)=qrhs(i,j,k,5)+( force(1)*u1+force(2)*u2+   &
    !                                 force(3)*u3 )*jacob(i,j,k)
    ! end do
    ! end do
    ! end do
    !
  end subroutine udf_src
  !+-------------------------------------------------------------------+
  !| The end of the subroutine udf_src.                                |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to defined an output by a user.                | 
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 18-Aug-2023: created by Jian Fang @ Daresbury                     |
  !+-------------------------------------------------------------------+
  subroutine udf_write
    !
    
    !
  end subroutine udf_write
  !+-------------------------------------------------------------------+
  !| The end of the subroutine udf_write.                              |
  !+-------------------------------------------------------------------+
  !
end module userdefine
!+---------------------------------------------------------------------+
!| The end of the module userdefine.                                   |
!+---------------------------------------------------------------------+
