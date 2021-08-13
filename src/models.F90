!+---------------------------------------------------------------------+
!| This module contains subroutines and variables related to turbulence|
!| model.                                                              |
!| ==============                                                      |
!| CHANGE RECORD                                                       |
!| -------------                                                       |
!| 09-Aug-2021  | Created by J. Fang @ STFC Daresbury Laboratory       |
!+---------------------------------------------------------------------+
module models
  !
  use parallel, only : mpirank,mpistop,irk,jrk,krk
  use commarray,only : tke,omg,miut,dtke,domg
  use constdef
  !
  implicit none
  !
  type :: komega_coef
    real(8) :: sigma_k1,sigma_k2,sigma_omega1,sigma_omega2,beta1,beta2,&
               beta_star,gamma1,gamma2,a1,b1,c1
    real(8) :: gamma,beta,prt
    real(8),allocatable,dimension(:,:,:) :: sigma_k,sigma_omega
  end type komega_coef
  real(8) :: kamma=0.41d0
  !
  type(komega_coef) :: komega
  !
  contains
  ! 
  !+-------------------------------------------------------------------+
  !| This subroutine is to allocate array for k-omega sst model.       |
  !+-------------------------------------------------------------------+
  !| ref: https://turbmodels.larc.nasa.gov/sst.html
  !|      https://www.openfoam.com/documentation/guides/latest/doc/guide-bcs-wall-turbulence-omegaWallFunction.html
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 09-Aug-2020: Created by J. Fang @ STFC Daresbury Laboratory       |
  !+-------------------------------------------------------------------+
  subroutine init_komegasst
    !
    use commvar,only : im,jm,km,hm
    !
    integer :: lallo
    !
    allocate(   tke(-hm:im+hm,-hm:jm+hm,-hm:km+hm),   &
                omg(-hm:im+hm,-hm:jm+hm,-hm:km+hm),   &
               miut(0:im,0:jm,0:km),stat=lallo)
    if(lallo.ne.0) stop ' !! error at allocating tke,omg,miut'
    !
    allocate( dtke(0:im,0:jm,0:km,1:3),               &
              domg(0:im,0:jm,0:km,1:3),stat=lallo )
    if(lallo.ne.0) stop ' !! error at allocating dtke,domg'
    !
    allocate( komega%sigma_k(0:im,0:jm,0:km),               &
              komega%sigma_omega(0:im,0:jm,0:km),stat=lallo )
    if(lallo.ne.0) stop ' !! error at allocating sigma_k,sigma_omega'
    !
    !
    komega%sigma_k1     = 0.85d0
    komega%sigma_k2     = 1.d0
    komega%sigma_omega1 = 0.5d0
    komega%sigma_omega2 = 0.856d0
    komega%beta1        = 0.075d0
    komega%beta2        = 0.0828d0
    komega%beta_star    = 0.09d0
    komega%a1           = 0.31d0
    komega%b1           = 1.d0
    komega%c1           = 10.d0
    !
    komega%prt          = 0.9d0
    !
    komega%gamma1=5.d0/9.d0
    komega%gamma2=0.44d0
    !
    ! komega%gamma1=komega%beta1/komega%beta_star - &
    !               komega%sigma_omega1*kamma*kamma/sqrt(komega%beta_star)
    ! komega%gamma2=komega%beta2/komega%beta_star - &
    !               komega%sigma_omega2*kamma*kamma/sqrt(komega%beta_star)
    !
  end subroutine init_komegasst
  !+-------------------------------------------------------------------+
  !| The end of the subroutine init_komegasst.                         |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to obtain eddy viscousity of k-omega sst model |
  !+-------------------------------------------------------------------+
  !| ref: https://turbmodels.larc.nasa.gov/sst.html                    |
  !| Menter, F. R. 1994 Two-Equation Eddy-Viscosity Turbulence Models  |
  !|   for Engineering Applications. AIAA J. 32,1598-605.              |
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 09-08-2021: Created by J. Fang @ Warrington.                      |
  !+-------------------------------------------------------------------+
  subroutine src_komega
    !
    use commvar,  only : im,jm,km,Reynolds,num_species
    use commarray,only : dis2wall,rho,tmp,vor,dvel,jacob,qrhs
    use fludyna,  only : miucal
    !
    ! local data
    integer :: i,j,k
    real(8) :: arg1,arg2,f1,f2,cdkomega,miu,sqrtk,dwall,vorti,ro,dkdomg
    real(8) :: var1,var2,var3,var4
    real(8) :: s11,s12,s13,s22,s23,s33,d11,d12,d13,d22,d23,d33,       &
               div,det,tau11,tau12,tau13,tau22,tau23,tau33,s
    real(8) :: produ_k,produ_omega,miueddy
    !
    do k=0,km
    do j=0,jm
    do i=0,im
      !
      miu=miucal(tmp(i,j,k))/Reynolds
      dwall=max(dis2wall(i,j,k),1.d-8)
      sqrtk=sqrt(tke(i,j,k))
      ro=rho(i,j,k)
      !
      dkdomg=dtke(i,j,k,1)*domg(i,j,k,1)+dtke(i,j,k,2)*domg(i,j,k,2) + &
             dtke(i,j,k,3)*domg(i,j,k,3)
      var2=2.d0*ro*komega%sigma_omega2/omg(i,j,k)*dkdomg
      cdkomega=max(var2,10.d-10)
      !
      var1=sqrtk/(komega%beta_star*omg(i,j,k)*dwall)
      var2=500.d0*miu/(ro*dwall*dwall*omg(i,j,k))
      var3=4.d0*ro*komega%sigma_omega2*tke(i,j,k)/(cdkomega*dwall*dwall)
      arg1=min(max(var1,var2),var3)
      arg2=max(2.d0*var1,var2)
      !
      f1=tanh(arg1**4)
      f2=tanh(arg2**2)
      !
      komega%sigma_k(i,j,k)    =komega%sigma_k1*f1+komega%sigma_k2*(1.d0-f1)
      komega%sigma_omega(i,j,k)=  komega%gamma1*f1+  komega%gamma2*(1.d0-f1)
      !
      komega%beta              =   komega%beta1*f1+   komega%beta2*(1.d0-f1)
      komega%gamma             =  komega%gamma1*f1+  komega%gamma2*(1.d0-f1)
      !
      vorti=sqrt(vor(i,j,k,1)**2+vor(i,j,k,2)**2+vor(i,j,k,3)**2)
      !
      s11=dvel(i,j,k,1,1)
      s12=0.5d0*(dvel(i,j,k,1,2)+dvel(i,j,k,2,1))
      s13=0.5d0*(dvel(i,j,k,1,3)+dvel(i,j,k,3,1))
      s22=dvel(i,j,k,2,2)
      s23=0.5d0*(dvel(i,j,k,2,3)+dvel(i,j,k,3,2))
      s33=dvel(i,j,k,3,3)
      !
      s=2.d0*(s11*s11+s22*s22+s33*s33+2.d0*(s12*s12+s13*s13*s23*s23))
      s=sqrt(s)
      !
      var4=max(komega%a1*omg(i,j,k),s*f2)
      !
      miut(i,j,k)=ro*komega%a1*tke(i,j,k)/var4
      miut(i,j,k)=min(miut(i,j,k),100.d0*miu)
      miut(i,j,k)=max(1.d-10,miut(i,j,k))
      !
      div=s11+s22+s33
      !
      det=num2d3*(miut(i,j,k)*div+ro*tke(i,j,k))
      tau11=2.d0*miut(i,j,k)*s11-det
      tau12=2.d0*miut(i,j,k)*s12
      tau13=2.d0*miut(i,j,k)*s13
      tau22=2.d0*miut(i,j,k)*s22-det
      tau23=2.d0*miut(i,j,k)*s23
      tau33=2.d0*miut(i,j,k)*s33-det
      !
      ! production term
      produ_k=      tau11*s11+tau22*s22+tau33*s33  +  &
              2.d0*(tau12*s12+tau13*s13+tau23*s23)
      produ_k=min(produ_k,20.d0*komega%beta_star*ro*tke(i,j,k)*omg(i,j,k))
      !
      var1=-komega%beta_star*ro*omg(i,j,k)*tke(i,j,k)
      !
      qrhs(i,j,k,6+num_species)=qrhs(i,j,k,6+num_species) +            &
                                             (produ_k+var1)*jacob(i,j,k)
      !
      miueddy=max(miut(i,j,k),1.d-10)
      produ_omega=komega%gamma*ro/miueddy*produ_k
      var1=-komega%beta*ro*omg(i,j,k)**2
      var2=2.d0*(1.d0-f1)*ro*komega%sigma_omega2/omg(i,j,k)*dkdomg
      !
      qrhs(i,j,k,7+num_species)=qrhs(i,j,k,7+num_species) +            &
                                    (produ_omega+var1+var2)*jacob(i,j,k)
      !
    enddo
    enddo
    enddo
    !
  end subroutine src_komega
  !+-------------------------------------------------------------------+
  !| The end of the subroutine src_komega.                             |
  !+-------------------------------------------------------------------+
  !
  ! subroutine komega_src
  !   !
  !   ! add source term to the k-omega transport equations
  !   !
  !   ! add production term
  !   do k=ks,ke
  !   do j=js,je
  !   do i=is,ie
  !     !
  !     s11=sigma(i,j,k,1)
  !     s12=sigma(i,j,k,2)
  !     s13=sigma(i,j,k,3)
  !     s22=sigma(i,j,k,4)
  !     s23=sigma(i,j,k,5)
  !     s33=sigma(i,j,k,6)
  !     !
  !     d11=dvel(i,j,k,1,1); d12=dvel(i,j,k,1,2); d13=dvel(i,j,k,1,3)
  !     d21=dvel(i,j,k,2,1); d22=dvel(i,j,k,2,2); d23=dvel(i,j,k,2,3)
  !     d31=dvel(i,j,k,3,1); d32=dvel(i,j,k,3,2); d33=dvel(i,j,k,3,3)
  !     !
  !     var1= s11*d11   +    s12*(d12+d21)  +   s13*(d13+d31) + &
  !                          s22*d22        +   s23*(d23+d32) + &
  !                                             s33*d33
  !     miueddy=max(miut(i,j,k),1.d-9)
  !     var2=var1*komega%gamma(i,j,k)*rho(i,j,k)/miueddy
  !     !
  !     n=5+num_species
  !     !
  !     qrhs(i,j,ks:ke,1+n)=qrhs(i,j,ks:ke,1+n) + var1*jacob(i,j,k)
  !     qrhs(i,j,ks:ke,2+n)=qrhs(i,j,ks:ke,2+n) + var2*jacob(i,j,k)
  !     !
  !   enddo
  !   enddo
  !   enddo
  !   !
  ! end subroutine komega_src
  !
end module models
!+---------------------------------------------------------------------+
!| The end of the module models.                                       |
!+---------------------------------------------------------------------+