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
    use commvar,  only: im,jm,km,ref_den,ref_tem
    use commarray,only: x,vel,rho,prs,tmp
    use parallel, only: lio
    use fludyna,  only: thermal
    !
    ! local data
    integer :: i,j,k
    !
    do k=0,km
    do j=0,jm
    do i=0,im
      !
      rho(i,j,k)=ref_den
      !
      vel(i,j,k,:)=0.d0
      !
      tmp(i,j,k)=ref_tem
      !
      prs(i,j,k)=thermal(density=rho(i,j,k),temperature=tmp(i,j,k))
    enddo
    enddo
    enddo
    !
    !
    if(lio)  write(*,'(A,I1,A)')'  ** HIT flame initialised.'
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
    !
    use commvar,  only : im,jm,km,gridfile,ia,ja,ka,ref_len
    use parallel, only : ig0,jg0,kg0,lio
    use commarray,only : x
    !
    ! local data
    integer :: i,j,k
    !
    do k=0,km
    do j=0,jm
    do i=0,im
      x(i,j,k,1)=ref_len/real(ia,8)*real(i+ig0,8)
      if(ja==0) then
        x(i,j,k,2)=ref_len
      else
        x(i,j,k,2)=ref_len/real(ja,8)*real(j+jg0,8)
      endif
      if(ka==0) then
        x(i,j,k,3)=0.d0
      else
        x(i,j,k,3)=ref_len/real(ka,8)*real(k+kg0,8)
      endif
      !
    enddo
    enddo
    enddo
    !
    if(lio) print*,' ** cubic grid generated'
    !
    !
  end subroutine udf_grid
  !+-------------------------------------------------------------------+
  !| The end of the subroutine udf_grid.                               |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to generate fluctuations for inflow            |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 05-Oct-2023: Created by by Jian Fang @ Daresbury                  |
  !+-------------------------------------------------------------------+
  subroutine udf_inflow_fluc(umean,uinst)
    !
    use commvar, only : jm,km
    !
    real(8),intent(in) ::  umean(0:jm,1:3)  ! inflow mean velocity
    real(8),intent(out) :: uinst(0:jm,0:km,1:3)  ! velocity with fluctuations
    !
  end subroutine udf_inflow_fluc
  !+-------------------------------------------------------------------+
  !| The end of the subroutine udf_inflow_fluc.                        |
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
    use commvar,  only : reynolds,lrestart,mach,ia,ja,ka,im,jm,km,nondimen
    use commarray,only : vel,rho,tmp,dvel,q
    use fludyna,  only : miucal,sos
    use utility,  only : listinit,listwrite
    use parallel, only : psum,lio
    !
    integer :: i,j,k,ns
    real(8) :: s11,s12,s13,s22,s23,s33,div,miu,dissa
    real(8) :: rsamples,miudrho,dudx2,csavg,v2,cs,ufluc,ens,omegax,omegay,omegaz
    real(8) :: urms,taylorlength,kolmoglength,Retaylor,machrms,macht,rhoe,skewness,du2,du3
    !
    logical :: fex
    logical,save :: linit=.true.
    integer,save :: hand_fs
    !
    if(linit) then
      !
      if(lio) then
        call listinit(filename='fturbstats.dat',handle=hand_fs, &
                      firstline='nstep time urms enstrophy taylorlength kolmoglength Retaylor machrms macht Tavg hsource skewness')
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
    du3=0.d0
    du2=0.d0
    do k=1,km
    do j=1,jm
    do i=1,im
      !
      if(nondimen) then
        miu=miucal(tmp(i,j,k))/Reynolds
      else
        miu=miucal(tmp(i,j,k))
      endif
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
      du3=du3+(s11*s11*s11+s22*s22*s22+s33*s33*s33)
      du2=du2+(s11*s11+s22*s22+s33*s33)
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

    skewness      = psum(du3)/(3.d0*rsamples)/sqrt((psum(du2)/(3.d0*rsamples))**3)
    ! kolmogvelocity= sqrt(sqrt(dissipation*miudrho))
    ! kolmogtime    = sqrt(miudrho/dissipation)
    !
    if(lio) call listwrite(hand_fs,urms,ens,taylorlength,kolmoglength, &
                           Retaylor,machrms,macht,rhoe,hsource,skewness)
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
  !| a random force acting like fans to input energy at largest scale  |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 13-06-2023: Created by Yifan Xu @ Peking University               |
  !+-------------------------------------------------------------------+
  subroutine udf_src
    !
    use constdef
    use commvar,  only : im,jm,km,ndims,deltat,ia,ja,ka,rkstep,xmax,ymax,zmax
    use parallel, only : mpirank,psum,bcast,lio
    use commarray,only : rho,tmp,vel,qrhs,x,jacob
    ! !
    ! ! local data
    integer,parameter :: nfan=1 !the number of fans
    !
    integer :: i,j,k,k1,k2,n
    real(8) :: force(3),ampA,ampB
    logical,save :: linit=.true.
    real(8),save :: A(3,3,nfan),B(3,3,nfan),phaseA(3,3,nfan),phaseB(3,3,nfan)
    real(8) :: FC,FR,var1,var2,tavg,at,xs,xe,xx,yy,zz,lwave,ran1,ran2,ran3,ran4
    real(8) :: xs1,xs2,xs3,xs4,xs5
    integer,allocatable :: seed(:)
    !
    if(linit) then
      !
      if(mpirank==0) then
        !
        open(12,file='datin/userinput.txt')
        read(12,*)FC,FR
        close(12)
        print*, ' ** FC =',FC,'FR=',FR
        !
        ! FC=0.d0
        ! FR=1.d0/16.d0
        ! !
        ! var1=sqrt(num2d3*FC/2.d-5)
        ! var2=sqrt(num1d3*FR/2.d-5)
        !
        call random_seed(size=n)
        allocate(seed(n))
        seed = 1    ! putting arbitrary seed to all elements
        call random_seed(put=seed)
        deallocate(seed)
        !
        do n=1,nfan
        do j=1,3
        do i=1,3
          call random_number(A(i,j,n))
          call random_number(B(i,j,n))
        end do
        end do
        end do
        !
        A=(2.d0*A-1.d0)
        B=(2.d0*B-1.d0)
        !
        do j=1,3
        do i=1,3
          if(i==j) then
            A(i,j,:)=0.d0
            B(i,j,:)=0.d0
          else
            A(i,j,:)=FR*A(i,j,:)
            B(i,j,:)=FR*B(i,j,:)
          endif
        end do
        end do
        !
        print*,' ** Amplitude of force, comp:',FC,' sole:',FR
        !
      endif
      !
      call bcast(A)
      call bcast(B)
      call bcast(phaseA)
      call bcast(phaseB)
      !
      linit=.false.
      !
    endif
    !
    lwave=ymax
    !
    xs=0.d0
    xe=xmax
    !
    hsource=0.d0
    tavg=0.d0
    !
    do k=0,km
    do j=0,jm
    do i=0,im
      !
      if(x(i,j,k,1)>=xs .and. x(i,j,k,1)<=xe) then
        !
        ! n=int((x(i,j,k,1)-xs)/lwave)+1
        !
        ! xx=(x(i,j,k,1)-xs)/lwave*2.d0*pi
        n=1
        !
        xx=x(i,j,k,1)/lwave*2.d0*pi
        yy=x(i,j,k,2)/lwave*2.d0*pi
        zz=x(i,j,k,3)/lwave*2.d0*pi
        !
        force(1)=A(1,1,n)*sin(xx)+B(1,1,n)*cos(xx) + &
                 A(1,2,n)*sin(yy)+B(1,2,n)*cos(yy) + &
                 A(1,3,n)*sin(zz)+B(1,3,n)*cos(zz)
        !
        force(2)=A(2,1,n)*sin(xx)+B(2,1,n)*cos(xx) + &
                 A(2,2,n)*sin(yy)+B(2,2,n)*cos(yy) + &
                 A(2,3,n)*sin(zz)+B(2,3,n)*cos(zz)
        !                                                          
        force(3)=A(3,1,n)*sin(xx)+B(3,1,n)*cos(xx) + &
                 A(3,2,n)*sin(yy)+B(3,2,n)*cos(yy) + &
                 A(3,3,n)*sin(zz)+B(3,3,n)*cos(zz)
        !
        qrhs(i,j,k,2)=qrhs(i,j,k,2)+rho(i,j,k)*force(1)*jacob(i,j,k)
        qrhs(i,j,k,3)=qrhs(i,j,k,3)+rho(i,j,k)*force(2)*jacob(i,j,k)
        qrhs(i,j,k,4)=qrhs(i,j,k,4)+rho(i,j,k)*force(3)*jacob(i,j,k)
        !
        qrhs(i,j,k,5)=qrhs(i,j,k,5)+rho(i,j,k)*( force(1)*vel(i,j,k,1) + &
                                                 force(2)*vel(i,j,k,2) + &
                                                 force(3)*vel(i,j,k,3) )*jacob(i,j,k)
        !
        if(i.ne.0 .and. j.ne.0 .and. k.ne.0) then
          hsource=hsource+rho(i,j,k)*(force(1)*vel(i,j,k,1) + &
                                      force(2)*vel(i,j,k,2) + &
                                      force(3)*vel(i,j,k,3))
          tavg=tavg+tmp(i,j,k)
        endif
        !
      endif
      !
    end do
    end do
    end do
    !
    at=psum(hsource)/psum(tavg)
    !
    do k=0,km
    do j=0,jm
    do i=0,im
      !
      if(x(i,j,k,1)>=xs .and. x(i,j,k,1)<=xe) then
        qrhs(i,j,k,5)=qrhs(i,j,k,5)-at*tmp(i,j,k)*jacob(i,j,k)
      endif
      !
    end do
    end do
    end do
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
  !+-------------------------------------------------------------------+
  !| This subroutine is to manipulate data solver as one likes at the  |
  !| end of each loop.                                                 | 
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 30-Oct-2023: created by Jian Fang @ Daresbury                     |
  !+-------------------------------------------------------------------+
  subroutine udf_eom_set
    !
  end subroutine udf_eom_set
  !+-------------------------------------------------------------------+
  !| The end of the subroutine udf_eom_set.                            |
  !+-------------------------------------------------------------------+
  !
end module userdefine
!+---------------------------------------------------------------------+
!| The end of the module userdefine.                                   |
!+---------------------------------------------------------------------+
