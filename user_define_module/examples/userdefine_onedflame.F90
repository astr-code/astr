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
  real(8) :: flamethickness,hsource,equivalence_ratio=0.8d0
  real(8),allocatable :: specr(:)
  !
  logical :: reset_burner=.false.
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
    use commvar,  only: ymax,roinf,uinf,vinf,winf,pinf,tinf,spcinf,num_species
    use fludyna,  only: thermal,sos
    use parallel, only: lio,bcast
    !
#ifdef COMB
    use thermchem,only : tranco,spcindex,mixture,convertxiyi,wmolar,rgcmix,gammarmix
    use cantera 
    !
    real(8) :: cpe,miu,kama,cs,lref
    real(8) :: dispec(num_species)
    real(8) :: airmass,o2frac
    real(8) :: molar_fraction_O2,molar_fraction_N2
    real(8) :: stoichiometric_H2_mol,stoichiometric_O2_mol,stoichiometric_N2_mol, &
               stoichiometric_H2_mas,stoichiometric_O2_mas,stoichiometric_N2_mas, &
               stoichiometric_H2_air_ratio,stoichiometric_O2_air_ratio,           &
               stoichiometric_N2_air_ratio,total_mass,massratio_H2_air
    ! 
    if(lio) then
      open(12,file='datin/userinput.txt')
      read(12,*)flamethickness
      close(12)
      print*, ' ** flamethickness =',flamethickness
    endif
    !
    call bcast(flamethickness)
    !
    stoichiometric_O2_mol=1.d0/(1+0.79d0/0.21d0+2.d0)
    stoichiometric_N2_mol=stoichiometric_O2_mol/0.21d0*0.79d0
    stoichiometric_H2_mol=stoichiometric_O2_mol*2.d0
    !
    stoichiometric_O2_mas=stoichiometric_O2_mol*wmolar(spcindex('O2'))
    stoichiometric_H2_mas=stoichiometric_H2_mol*wmolar(spcindex('H2'))
    stoichiometric_N2_mas=stoichiometric_N2_mol*wmolar(spcindex('N2'))
    !
    total_mass=stoichiometric_O2_mas+stoichiometric_H2_mas+stoichiometric_N2_mas
    !
    stoichiometric_H2_air_ratio=stoichiometric_H2_mas/(stoichiometric_O2_mas+stoichiometric_N2_mas)
    stoichiometric_O2_air_ratio=stoichiometric_O2_mas/(stoichiometric_O2_mas+stoichiometric_N2_mas)
    stoichiometric_N2_air_ratio=stoichiometric_N2_mas/(stoichiometric_O2_mas+stoichiometric_N2_mas)
    !
    if(lio) then
      print*,' ** stoichiometric fuel/air condition'
      write(*,'(1X,A)')'     molar fraction'
      write(*,'(1X,A,1x,F10.5)')'                H2:',stoichiometric_H2_mol
      write(*,'(1X,A,1x,F10.5)')'                O2:',stoichiometric_O2_mol
      write(*,'(1X,A,1x,F10.5)')'                N2:',stoichiometric_N2_mol
      write(*,'(1X,A)')'      mass fraction'
      write(*,'(1X,A,1x,F10.5)')'                H2:',stoichiometric_H2_mas/total_mass
      write(*,'(1X,A,1x,F10.5)')'                O2:',stoichiometric_O2_mas/total_mass
      write(*,'(1X,A,1x,F10.5)')'                N2:',stoichiometric_N2_mas/total_mass
      write(*,'(1X,A)')'         mass ratio'
      write(*,'(1X,A,1x,F10.5)')'            H2/air:',stoichiometric_H2_air_ratio
      write(*,'(1X,A,1x,F10.5)')'            O2/air:',stoichiometric_O2_air_ratio
      write(*,'(1X,A,1x,F10.5)')'            N2/air:',stoichiometric_N2_air_ratio
    endif
    !
    !
    if(.not. allocated(spcinf)) allocate(spcinf(num_species))
    if(.not. allocated(specr)) allocate(specr(num_species))
    !
    specr(:)=0.d0
    !
    massratio_H2_air=equivalence_ratio*stoichiometric_H2_air_ratio
    !
    specr(spcindex('H2'))=massratio_H2_air/(massratio_H2_air+1.d0)
    specr(spcindex('O2'))=1.d0/(massratio_H2_air+1.d0)*stoichiometric_O2_air_ratio
    specr(spcindex('N2'))=1.d0/(massratio_H2_air+1.d0)*stoichiometric_N2_air_ratio
    !
    uinf=0.01d0
    vinf=0.d0
    winf=0.d0
    tinf=300.d0
    pinf=101325.d0
    !
    spcinf(:)=specr(:)
    roinf=thermal(pressure=pinf,temperature=tinf,species=spcinf(:))
    !
    cs=sos(tinf,spcinf)
    !
    lref=flamethickness
    !
    call tranco(den=roinf,tmp=tinf,cp=cpe,mu=miu,lam=kama, &
                spc=specr,rhodi=dispec)

    if(lio) then
      !
      print*,' ---------------------------------------------------------------'
      print*,'                      free stream quatities                     '
      print*,' --------------------------+------------------------------------'
      print*,'                        u∞ | ',uinf,'m/s'
      print*,'                        T∞ | ',tinf,'K'
      print*,'                      rho∞ | ',roinf,'kg/m**3'
      print*,'                        p∞ | ',pinf,'Pa'
      print*,'          reference length | ',lref,'m'
      print*,'                 viscosity | ',miu,'kg/(ms)'
      print*,'                       Re∞ | ',roinf*uinf*lref/miu
      print*,'            speed of sound | ',cs,'m/s'
      print*,'                       Ma∞ | ',uinf/cs
      print*,'         equivalence ratio | ',equivalence_ratio
      print*,'          H2 mass fraction | ',specr(spcindex('H2'))
      print*,'          O2 mass fraction | ',specr(spcindex('O2'))
      print*,'          N2 mass fraction | ',specr(spcindex('N2'))
      print*,'           gas constant, R | ',rgcmix(specr),' J/(kg K)'
      print*,'    heat capacity ratio, γ | ',gammarmix(tinf,specr)
      print*,' --------------------------+------------------------------------'

    endif
    !
#endif
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
    use commvar,  only: im,jm,km,ndims,roinf,uinf,nondimen,xmax,pinf,  &
                        ia,num_species,nstep,time,filenumb
    use commarray,only: x,vel,rho,prs,spc,tmp,q
    use parallel, only: lio
    use fludyna,  only: thermal,updateq
    !
#ifdef COMB
    !
    use thermchem,only : tranco,spcindex,mixture,convertxiyi
    use cantera 
    !
    ! local data
    integer :: i,j,k
    real(8) ::  xc,yc,zc,tmpr,tmpp,xloc,xwid,specp(num_species),arg,prgvar,masflx,specx(num_species)
    real(8) :: pthick
    !
    tmpr=300.d0
    xloc=4.d0*xmax/5.d0
    xwid=xmax/(12.d0*5.3d0)
    !
    !reactants
    ! specr(:)=0.d0
    ! specr(spcindex('H2'))=0.0288
    ! specr(spcindex('O2'))=0.2289
    ! specr(spcindex('N2'))=1.d0-sum(specr)
    !
    !products
    !tmpp=1814.32d0
    tmpp=2814.32d0
    !
    pthick=5.d-4
    !
    do k=0,km
    do j=0,jm
    do i=0,im
      !
      xc=x(i,j,k,1)
      !
      !prgvar=0.5d0*(1.d0+tanh(10.d0*(xc-xloc)/xloc))
      ! if(xc-xloc<xwid*0.5d0*1.2d0) then 
      !   prgvar=0.d0
      !   if(xc-xloc>xwid*0.5d0) &
      !   prgvar=1.d0-(xc-xloc-(xwid*0.5d0))/(xwid*0.5d0*0.2d0)
      ! else
      !   prgvar=1.d0
      ! endif
      !
      ! prgvar=1.d0*exp(-0.5d0*((xc-xloc)/xwid)**2)
      prgvar=1.d0*exp(-0.5d0*((xc-xloc)/pthick)**2)
      !
      spc(i,j,k,:)=specr(:)
      !
      vel(i,j,k,1)=uinf
      !
      vel(i,j,k,2)=0.d0
      vel(i,j,k,3)=0.d0
      !
      tmp(i,j,k)=tmpr+prgvar*(tmpp-tmpr)
      !
      prs(i,j,k)=pinf
      !
      rho(i,j,k)=thermal(pressure=prs(i,j,k),temperature=tmp(i,j,k), &
                          species=spc(i,j,k,:))
    enddo
    enddo
    enddo
    !
    call updateq
    !
    nstep=0
    time=0.d0
    filenumb=0
    !
    if(lio)  write(*,'(A,I1,A)')'  ** HIT flame initialised.'
    !
#endif
    !
  end subroutine udf_flowinit
  !+-------------------------------------------------------------------+
  !| The end of the subroutine udf_flowinit.                           |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to collect statistics.                         | 
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 09-Nov-2023: created by Jian Fang @ Appleton                      |
  !+-------------------------------------------------------------------+
  subroutine udf_meanflow
    !
  end subroutine udf_meanflow
  !+-------------------------------------------------------------------+
  !| The end of the subroutine udf_meanflow.                           |
  !+-------------------------------------------------------------------+
  !!
  !+-------------------------------------------------------------------+
  !| This subroutine is to generate grid.                              | 
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 23-Aug-2023: created by Jian Fang @ Daresbury                     |
  !+-------------------------------------------------------------------+
  subroutine udf_grid(mode)
    !
    use commvar,  only : im,jm,km,gridfile,ia,ja,ka
    use parallel, only : ig0,jg0,kg0,lio,bcast
    use commarray,only : x
    use hdf5io
    !
    ! arguments
    character(len=*),intent(in),optional :: mode
    !
    ! local data
    integer :: i,j,k
    real(8) :: lx,ly,lz
    !
    ! if(mode=='cuboid') then
      ! lx=12.d0*5.3d0*flamethickness
      lx=24.d0*5.3d0*flamethickness
      ly=5.3d0*flamethickness
      lz=5.3d0*flamethickness
    ! elseif(mode=='cubic') then
    !   lx=5.3d0*flamethickness
    !   ly=5.3d0*flamethickness
    !   lz=5.3d0*flamethickness
    ! else
    !   stop ' !! error1 @ gridhitflame'
    ! endif
    !
    do k=0,km
    do j=0,jm
    do i=0,im
      x(i,j,k,1)=lx/real(ia,8)*real(i+ig0,8)
      if(ja==0) then
        x(i,j,k,2)=0.d0
      else
        x(i,j,k,2)=ly/real(ja,8)*real(j+jg0,8)
      endif
      if(ka==0) then
        x(i,j,k,3)=0.d0
      else
        x(i,j,k,3)=lz/real(ka,8)*real(k+kg0,8)
      endif
      !
    enddo
    enddo
    enddo
    !
    if(lio) print*,' ** cubic grid generated, lx=',lx
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
    use commvar,  only : im,jm,km,ia,ja,ka,deltat,time,uinf
    use commarray,only : vel,prs,rho,tmp,spc,cell,x
    use interp,   only : interlinear
    use parallel, only : pmax,psum,irk,irkm,lio
    use utility,  only : listinit,listwrite,get_unit
    !
    use thermchem,only : heatrate
    !
    integer :: i,j,k
    real(8) :: tmpmax,rhomax,umax,qdotmax,poutrt
    real(8) :: qdot,var1,var2,bvelo
    !
    integer,save :: hand_fs,hand_fs2
    real(8),save :: xflame=0.d0,vflame=0.d0,xflame1,xflame2,tflame1,tflame2
    logical,save :: linit=.true.,monitor1=.true.
    !
    if(lio) then
      !
      if(linit) then
        call listinit(filename='flamesta.dat',handle=hand_fs, &
                      firstline='nstep time maxT maxU maxHRR xflame vflame pout')
        !
        hand_fs2=get_unit()
        open(hand_fs2,file='laminar_burnning_speed.txt')
        write(hand_fs2,'(4(1X,A20))')'equivalence_ratio','burning_velocity','maxT','maxHRR'
        !
        linit=.false.
        !
      endif
      !
    endif
    !
    tmpmax=maxval(tmp(0:im,0:jm,0:km))
    tmpmax=pmax(tmpmax)
    !
    rhomax=maxval(rho(0:im,0:jm,0:km))
    rhomax=pmax(rhomax)
    !
    umax=maxval(vel(0:im,0:jm,0:km,1))
    umax=pmax(umax)
    !
    var1=0.d0
    var2=0.d0
    !
    qdotmax=-1.d30
    !
    do i=0,im
      do j=0,jm
        do k=0,km
          !
          qdot=heatrate(rho(i,j,k),tmp(i,j,k),spc(i,j,k,:))
          if(qdot>qdotmax) then
            qdotmax=qdot
          endif
          !
        enddo 
      enddo 
    enddo  
    qdotmax=pmax(qdotmax)
    !
    ! calculate the averaged flame location, set as the T=400K
    var1=0.d0
    !
    j=0
    do i=1,im
      !
      if( (tmp(i-1,j,k)<=400.d0 .and. tmp(i,j,k)>=400.d0) ) then
        var1=var1+interlinear(tmp(i-1,j,k),tmp(i,j,k),        &
                              x(i-1,j,k,1),x(i,j,k,1),400.d0)
        exit
      endif
      !
    enddo
    !
    ! var1=psum(var1)/ja
    var1=psum(var1) !/ja
    !
    ! use xflame to calculate vflame
    if(abs(xflame)<1.d-16) then
      vflame=0.d0
    else
      vflame=(var1-xflame)/deltat
    endif
    !
    xflame=var1
    !
    if(xflame<=0.002d0 .and. xflame>1.d-8 .and. monitor1) then
      !
      xflame1=xflame
      tflame1=time
      !
      monitor1=.false.
      !
    endif
    !
    if(xflame<=0.001d0 .and. xflame>1.d-8) then
      !
      xflame2=xflame
      tflame2=time
      !
      if(lio) then
        !
        bvelo=(xflame1-xflame2)/(tflame2-tflame1)-uinf
        !
        write(hand_fs2,'(4(1X,E20.13E2))')equivalence_ratio,bvelo,tmpmax,qdotmax
        !
        print*,'burning velocity:',equivalence_ratio,bvelo
        !
      endif
      !
      monitor1=.true.
      !
      reset_burner=.true.
      !
      call udf_write
      !
    endif
    !
    ! calculate mean pressure at outflow
    i=im
    k=0
    poutrt=0.d0
    !
    if(irk==irkm) then
      ! do j=1,jm
      !   poutrt=poutrt+prs(i,j,k)
      ! enddo
      poutrt=poutrt+prs(im,0,0)
    endif
    !
    poutrt=psum(poutrt)
    ! poutrt=psum(poutrt)/dble(ja)
    !
    if(lio) call listwrite(hand_fs,tmpmax,umax,qdotmax,xflame,vflame,poutrt)
    !
  end subroutine udf_stalist
  !+-------------------------------------------------------------------+
  !| The end of the subroutine udf_stalist.                            |
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
    use parallel, only: mpistop
    !
    if(reset_burner) then
      !
      ! equivalence_ratio=equivalence_ratio+0.2d0
      ! !
      ! call udf_setflowenv
      ! !
      ! call udf_flowinit
      ! !
      ! reset_burner=.false.
      !
      call mpistop
      !
    endif
    !
  end subroutine udf_eom_set
  !+-------------------------------------------------------------------+
  !| The end of the subroutine udf_eom_set.                            |
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
    ! use constdef
    ! use commvar,  only : im,jm,km,ndims,deltat,ia,ja,ka,rkstep,xmax,ymax,zmax
    ! use parallel, only : mpirank,psum,bcast
    ! use commarray,only : rho,tmp,vel,qrhs,x,jacob
    ! ! !
    ! ! ! local data
    ! integer,parameter :: nfan=5 !the number of fans
    ! !
    ! integer :: i,j,k,k1,k2,n
    ! real(8) :: force(3)
    ! logical,save :: linit=.true.
    ! real(8),save :: A(3,3,nfan),B(3,3,nfan)
    ! real(8) :: FC,FR,var1,var2,tavg,at,xs,xe,xx,yy,zz,lwave
    ! real(8) :: xs1,xs2,xs3,xs4,xs5
    ! integer,allocatable :: seed(:)
    ! !
    ! if(linit) then
    !   !
    !   if(mpirank==0) then
    !     !
    !     FR=1.d0/16.d0
    !     FC=0.d0
    !     !
    !     var1=sqrt(num2d3*FC/deltat)
    !     var2=sqrt(num1d3*FR/deltat)
    !     !
    !     call random_seed(size=n)
    !     allocate(seed(n))
    !     seed = 1    ! putting arbitrary seed to all elements
    !     call random_seed(put=seed)
    !     deallocate(seed)
    !     !
    !     do n=1,nfan
    !     do j=1,3
    !     do i=1,3
    !       call random_number(A(i,j,n))
    !       call random_number(B(i,j,n))
    !     end do
    !     end do
    !     end do
    !     !
    !     A=sqrt(3.d0)*(2.d0*A-1.d0)
    !     B=sqrt(3.d0)*(2.d0*B-1.d0)
    !     !
    !     do j=1,3
    !     do i=1,3
    !       if(i==j) then
    !         A(i,j,:)=var1*A(i,j,:)
    !         B(i,j,:)=var1*B(i,j,:)
    !       else
    !         A(i,j,:)=var2*A(i,j,:)
    !         B(i,j,:)=var2*B(i,j,:)
    !       endif
    !     end do
    !     end do
    !     !
    !   endif
    !   !
    !   call bcast(A)
    !   call bcast(B)
    !   !
    !   A=0.9594d0*A
    !   B=0.9594d0*B
    !   !
    !   linit=.false.
    !   !  
    ! endif
    ! !
    ! lwave=ymax
    ! !
    ! xs=1.1065d-2
    ! xe=xs+4.d0*lwave
    ! !
    ! hsource=0.d0
    ! tavg=0.d0
    ! !
    ! do k=0,km
    ! do j=0,jm
    ! do i=0,im
    !   !
    !   if(x(i,j,k,1)>=xs .and. x(i,j,k,1)<=xe) then
    !     !
    !     n=int((x(i,j,k,1)-xs)/lwave)+1
    !     !
    !     xx=(x(i,j,k,1)-xs)/lwave*2.d0*pi
    !     yy=x(i,j,k,2)/lwave*2.d0*pi
    !     zz=x(i,j,k,3)/lwave*2.d0*pi
    !     !
    !     force(1)=A(1,1,n)*sin(xx)+B(1,1,n)*cos(xx) + &
    !              A(1,2,n)*sin(yy)+B(1,2,n)*cos(yy) + &
    !              A(1,3,n)*sin(zz)+B(1,3,n)*cos(zz)
    !     !
    !     force(2)=A(2,1,n)*sin(xx)+B(2,1,n)*cos(xx) + &
    !              A(2,2,n)*sin(yy)+B(2,2,n)*cos(yy) + &
    !              A(2,3,n)*sin(zz)+B(2,3,n)*cos(zz)
    !     !                                                          
    !     force(3)=A(3,1,n)*sin(xx)+B(3,1,n)*cos(xx) + &
    !              A(3,2,n)*sin(yy)+B(3,2,n)*cos(yy) + &
    !              A(3,3,n)*sin(zz)+B(3,3,n)*cos(zz)
    !     !
    !     qrhs(i,j,k,2)=qrhs(i,j,k,2)+rho(i,j,k)*force(1)*jacob(i,j,k)
    !     qrhs(i,j,k,3)=qrhs(i,j,k,3)+rho(i,j,k)*force(2)*jacob(i,j,k)
    !     qrhs(i,j,k,4)=qrhs(i,j,k,4)+rho(i,j,k)*force(3)*jacob(i,j,k)
    !     !
    !     qrhs(i,j,k,5)=qrhs(i,j,k,5)+rho(i,j,k)*( force(1)*vel(i,j,k,1) + &
    !                                              force(2)*vel(i,j,k,2) + &
    !                                              force(3)*vel(i,j,k,3) )*jacob(i,j,k)
    !     !
    !     if(i.ne.0 .and. j.ne.0 .and. k.ne.0) then
    !       hsource=hsource+rho(i,j,k)*(force(1)*vel(i,j,k,1) + &
    !                                   force(2)*vel(i,j,k,2) + &
    !                                   force(3)*vel(i,j,k,3))
    !       tavg=tavg+tmp(i,j,k)
    !     endif
    !     !
    !   endif
    !   !
    ! end do
    ! end do
    ! end do
    ! !
    ! at=psum(hsource)/psum(tavg)
    ! !
    ! do k=0,km
    ! do j=0,jm
    ! do i=0,im
    !   !
    !   if(x(i,j,k,1)>=xs .and. x(i,j,k,1)<=xe) then
    !     qrhs(i,j,k,5)=qrhs(i,j,k,5)-at*tmp(i,j,k)*jacob(i,j,k)
    !   endif
    !   !
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
    use commvar,  only: ymin,ymax,im,jm,filenumb
    use commarray,only: x,rho,vel,tmp,prs,spc
    use parallel, only: mpistop
    !
#ifdef COMB
    use thermchem,only : heatrate
#endif
    !
    integer :: i,j
    logical :: lwprofile
    real(8) :: ypos
    real(8),allocatable :: hrr(:)
    character(len=4) :: stepname,eqrname
    character(len=128) :: filename
    !
    lwprofile=.true.
    !
    ! ypos=0.5d0*(ymax-ymin)+ymin
    ! !
    ! do j=1,jm
    !   if(x(0,j-1,0,2)<ypos .and. x(0,j,0,2)>=ypos) then
    !     !
    !     lwprofile=.true.
    !     !
    !     exit
    !     !
    !   endif
    ! enddo
    !
#ifdef COMB
    !
    allocate(hrr(0:im))
    do i=0,im
      hrr(i)=heatrate(rho(i,0,0),tmp(i,0,0),spc(i,0,0,:))
    enddo
    !
    write(stepname,'(i4.4)')filenumb
    write(eqrname,'(F4.2)')equivalence_ratio
    !
    !
    if(reset_burner) then
      filename='outdat/profile_equivalence_ratio_'//eqrname//'_end.dat'
    else
      filename='outdat/profile_equivalence_ratio_'//eqrname//'_'//trim(stepname)//'.dat'
    endif
    !
    call writexprofile(profilename=trim(filename),  &
                               var1=rho(0:im,j,0),  var1name='rho', &
                               var2=vel(0:im,j,0,1),var2name='u',   &
                               var3=tmp(0:im,j,0),  var3name='T',   &
                               var4=prs(0:im,j,0),  var4name='P',   &
                               var5=hrr(0:im),      var5name='HRR',truewrite=lwprofile)
    !

    !
#endif
    !
  end subroutine udf_write
  !+-------------------------------------------------------------------+
  !| The end of the subroutine udf_write.                              |
  !+-------------------------------------------------------------------+
  !
  subroutine writexprofile(profilename,var1,var1name, &
                                       var2,var2name, &
                                       var3,var3name, &
                                       var4,var4name, &
                                       var5,var5name,truewrite)
    !
    use commvar,   only : im 
    use parallel,  only : pgather,mpirank
    use commarray, only : x
    !
    character(len=*),intent(in) :: profilename
    !
    real(8),intent(in),optional :: var1(:),var2(:),var3(:),var4(:),var5(:)
    character(len=*),intent(in),optional :: var1name,var2name,var3name,var4name,var5name
    logical,intent(in) :: truewrite
    !
    integer :: nvar,i
    real(8),allocatable :: vdum(:)
    real(8),allocatable :: vout1(:),vout2(:),vout3(:),vout4(:),vout5(:)
    real(8),allocatable,save :: xx(:)
    logical,save :: firstcall=.true.
    !
    if(firstcall) then
      if(truewrite) then
        allocate(vdum(0:im))
        vdum=x(0:im,0,0,1)
      endif
      call pgather(vdum,xx)
      if(truewrite) deallocate(vdum)
      !
      firstcall=.false.
    endif
    !
    nvar=0
    !
    if(truewrite) allocate(vdum(1:size(var1)))
    !
    if(present(var1)) then
      nvar=1
      !
      if(truewrite) then
        vdum=var1
      endif
      !
      call pgather(vdum,vout1)
      !
    endif
    !
    if(present(var2)) then
      nvar=2
      if(truewrite) then
        vdum=var2
      endif
      call pgather(vdum,vout2)
      !
    endif
    !
    if(present(var3)) then
      nvar=3
      if(truewrite) then
        vdum=var3
      endif
      call pgather(vdum,vout3)
      !
    endif
    !
    if(present(var4)) then
      nvar=4
      if(truewrite) then
        vdum=var4
      endif
      call pgather(vdum,vout4)
      !
    endif
    !
    if(present(var5)) then
      nvar=5
      if(truewrite) then
        vdum=var5
      endif
      call pgather(vdum,vout5)
      !
    endif
    !
    if(nvar==0) return
    !
    if(mpirank==0) then
      open(18,file=profilename)
      if(nvar==1) then
        write(18,"(2(1X,A15))")'x',var1name
        write(18,"(2(1X,E15.7E3))")(xx(i),vout1(i),i=1,size(xx))
      elseif(nvar==2) then
        write(18,"(3(1X,A15))")'x',var1name,var2name
        write(18,"(3(1X,E15.7E3))")(xx(i),vout1(i),vout2(i),i=1,size(xx))
      elseif(nvar==3) then
        write(18,"(4(1X,A15))")'x',var1name,var2name,var3name
        write(18,"(4(1X,E15.7E3))")(xx(i),vout1(i),vout2(i),vout3(i),i=1,size(xx))
      elseif(nvar==4) then
        write(18,"(5(1X,A15))")'x',var1name,var2name,var3name,var4name
        write(18,"(5(1X,E15.7E3))")(xx(i),vout1(i),vout2(i),vout3(i),vout4(i),i=1,size(xx))
      elseif(nvar==5) then
        write(18,"(6(1X,A15))")'x',var1name,var2name,var3name,var4name,var5name
        write(18,"(6(1X,E15.7E3))")(xx(i),vout1(i),vout2(i),vout3(i),vout4(i),vout5(i),i=1,size(xx))
      else
        stop ' !! error1 @ writexprofile'
      endif
      close(18)
      print*,' << ',profilename
    endif
    !
  end subroutine writexprofile
  !
end module userdefine
!+---------------------------------------------------------------------+
!| The end of the module userdefine.                                   |
!+---------------------------------------------------------------------+
