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
!     use commvar,  only : im,jm,km,ia,ja,ka,deltat
!     use commarray,only : vel,prs,rho,tmp,spc,cell,x
!     use interp,   only : interlinear
!     use parallel, only : pmax,psum,irk,irkm,lio
!     use utility,  only : listinit,listwrite
!     !
! #ifdef COMB
!     use thermchem,only : heatrate,spcindex
! #endif
!     !
!     integer :: i,j,k
!     real(8) :: tmpmax,rhomax,umax,qdotmax,poutrt
!     real(8) :: qdot,var1,var2
!     !
!     integer,save :: hand_fs
!     real(8),save :: xflame=0.d0,vflame=0.d0
!     logical,save :: linit=.true.
!     !
!     if(lio) then
!       !
!       if(linit) then
!         call listinit(filename='flamesta.dat',handle=hand_fs, &
!                       firstline='nstep time maxT maxU maxHRR xflame vflame pout')
!         linit=.false.
!       endif
!       !
!     endif
!     !
!     tmpmax=maxval(tmp(0:im,0:jm,0:km))
!     tmpmax=pmax(tmpmax)
!     !
!     rhomax=maxval(rho(0:im,0:jm,0:km))
!     rhomax=pmax(rhomax)
!     !
!     umax=maxval(vel(0:im,0:jm,0:km,1))
!     umax=pmax(umax)
!     !
!     var1=0.d0
!     var2=0.d0
!     !
!     qdotmax=-1.d30
!     !
!     do i=0,im
!       do j=0,jm
!         do k=0,km
!           !
!           qdot=heatrate(rho(i,j,k),tmp(i,j,k),spc(i,j,k,:))
!           if(qdot>qdotmax) then
!             qdotmax=qdot
!           endif
!           !
!         enddo 
!       enddo 
!     enddo  
!     qdotmax=pmax(qdotmax)
!     !
!     ! calculate the averaged flame location, set as the T=400K
!     var1=0.d0
!     !
!     do j=1,jm
!       !
!       do i=1,im
!         !
!         if( (tmp(i-1,j,k)<=400.d0 .and. tmp(i,j,k)>=400.d0) ) then
!           var1=var1+interlinear(tmp(i-1,j,k),tmp(i,j,k),        &
!                                 x(i-1,j,k,1),x(i,j,k,1),400.d0)
!           exit
!         endif
!         !
!       enddo
!       !
!     enddo
!     !
!     var1=psum(var1)/ja
!     !
!     ! use xflame to calculate vflame
!     if(abs(xflame)<1.d-16) then
!       vflame=0.d0
!     else
!       vflame=(var1-xflame)/deltat
!     endif
!     !
!     xflame=var1
!     !
!     ! calculate mean pressure at outflow
!     i=im
!     k=0
!     poutrt=0.d0
!     !
!     if(irk==irkm) then
!       do j=1,jm
!         poutrt=poutrt+prs(i,j,k)
!       enddo
!     endif
!     !
!     poutrt=psum(poutrt)/dble(ja)
!     !
!     if(lio) call listwrite(hand_fs,tmpmax,umax,qdotmax,xflame,vflame,poutrt)
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
    ! use commvar,  only : force,im,jm,km,ndims
    ! use parallel, only : psum 
    ! use commarray,only : q,qrhs,x,jacob
    ! !
    ! ! local data
    ! integer :: i,j,k,k1,k2
    ! !
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


  subroutine udf_bc(ndir)
    
    use constdef
    use commarray, only : x,prs,vel,tmp,rho,spc,q
    use commvar,   only : nondimen,im,km,jm,num_species,time
    use parallel,  only : irk,jrk,jrkm,irkm
    use fludyna,   only : thermal,fvar2q

    ! arguments
    integer,intent(in) :: ndir

    ! local data
    integer :: i,j,k
    real(8) :: xs,ys,vss,p2,ro2,u2,t2,a2,p1,ro1,u1,t1,a1

    real(8) :: r00,p00,u00,v00,r01,p01,u01,v01,r10,p10,u10,v10,r11,p11,u11,v11


     r01 = 0.5323d0;  r11 = 1.5d0
     u01 = 1.206d0;   u11 = 0.d0
     v01 = 0.d0;      v11 = 0.d0
     p01 = 0.30;      p11 = 1.5d0

     r00 = 0.138d0;   r10 = 0.5323d0
     u00 = 1.206d0;   u10 = 0.d0
     v00 = 1.206d0;   v10 = 1.206d0
     p00 = 0.029d0;   p10 = 0.3d0
  

    if(ndir==1 .and. irk==0) then
      
        ro1=r01
        p1 =p01
        u1 =v01
        t1 =thermal(pressure=p1,density=ro1)

        ro2=r00
        p2 =p00
        u2 =v00

        vss=shockspeed(p1,t1,p2)

        i=0


        do k=0,km
        do j=0,jm

          ys=0.8d0-vss*time

          if(x(i,j,k,2)>ys) then
            rho(i,j,k)  = ro1
            vel(i,j,k,2)= u1
            prs(i,j,k)  = p1
          else
            rho(i,j,k)  = ro2
            vel(i,j,k,2)= u2
            prs(i,j,k)  = p2
          endif

          vel(i,j,k,1)= u01
          vel(i,j,k,3)= 0.d0

          tmp(i,j,k)  =thermal(pressure=prs(i,j,k),density=rho(i,j,k))

          call fvar2q(        q=  q(i,j,k,:),                      &
                        density=rho(i,j,k),                        &
                       velocity=vel(i,j,k,:),                      &
                       pressure=prs(i,j,k),                        &
                        species=spc(i,j,k,:)                       )
        enddo
        enddo

    elseif(ndir==2 .and. irk==irkm) then
      
        ro1=r11
        p1 =p11
        u1 =v11
        t1 =thermal(pressure=p1,density=ro1)

        ro2=r10
        p2 =p10
        u2 =v10
        vss=shockspeed(p1,t1,p2)

        i=im
        do k=0,km
        do j=0,jm

          ys=0.8d0-vss*time

          ! if(x(i,j,k,2)>ys) then
          !   rho(i,j,k)  = ro1
          !   vel(i,j,k,2)= u1
          !   prs(i,j,k)  = p1 
          ! else
          !   rho(i,j,k)  = ro2
          !   vel(i,j,k,2)= u2
          !   prs(i,j,k)  = p2
          ! endif
          rho(i,j,k)=num1d3*(4.d0*rho(i-1,j,k)-rho(i-2,j,k))
          prs(i,j,k)=num1d3*(4.d0*prs(i-1,j,k)-prs(i-2,j,k))
          vel(i,j,k,2)=num1d3*(4.d0*vel(i-1,j,k,2)-vel(i-2,j,k,2))
          vel(i,j,k,1)= 0.d0
          vel(i,j,k,3)= 0.d0

          tmp(i,j,k)  =thermal(pressure=prs(i,j,k),density=rho(i,j,k))

          call fvar2q(        q=  q(i,j,k,:),                      &
                        density=rho(i,j,k),                        &
                       velocity=vel(i,j,k,:),                      &
                       pressure=prs(i,j,k),                        &
                        species=spc(i,j,k,:)                       )
        enddo
        enddo

        if(jrk==jrkm) then
          i=im
          j=jm
          do k=0,km
  
            rho(i,j,k)=r11
            prs(i,j,k)=p11
            vel(i,j,k,2)=v11
            vel(i,j,k,1)= 0.d0
            vel(i,j,k,3)= 0.d0
  
            tmp(i,j,k)  =thermal(pressure=prs(i,j,k),density=rho(i,j,k))
  
            call fvar2q(        q=  q(i,j,k,:),                      &
                          density=rho(i,j,k),                        &
                         velocity=vel(i,j,k,:),                      &
                         pressure=prs(i,j,k),                        &
                          species=spc(i,j,k,:)                       )
          enddo
        endif

    elseif(ndir==3 .and. jrk==0) then
      
        ro1=r10
        p1 =p10
        u1 =u10
        t1 =thermal(pressure=p1,density=ro1)

        ro2=r00
        p2 =p00
        u2 =u00
        vss=shockspeed(p1,t1,p2)


        j=0
        do k=0,km
        do i=0,im

          xs=0.8d0-vss*time

          if(x(i,j,k,1)>xs) then
            rho(i,j,k)  = ro1
            vel(i,j,k,1)= u1
            prs(i,j,k)  = p1
          else
            rho(i,j,k)  = ro2
            vel(i,j,k,1)= u2
            prs(i,j,k)  = p2
          endif

          vel(i,j,k,2)= v10
          vel(i,j,k,3)= 0.d0

          tmp(i,j,k)  =thermal(pressure=prs(i,j,k),density=rho(i,j,k))

          call fvar2q(        q=  q(i,j,k,:),                      &
                        density=rho(i,j,k),                        &
                       velocity=vel(i,j,k,:),                      &
                       pressure=prs(i,j,k),                        &
                        species=spc(i,j,k,:)                       )
        enddo
        enddo

    elseif(ndir==4 .and. jrk==jrkm) then
      
        ro1=r11
        p1 =p11
        u1 =u11
        t1 =thermal(pressure=p1,density=ro1)

        ro2=r01
        p2 =p01
        u2 =u01
        vss=shockspeed(p1,t1,p2)

        j=jm
        do k=0,km
        do i=0,im

          xs=0.8d0-vss*time

          rho(i,j,k)=num1d3*(4.d0*rho(i,j-1,k)-rho(i,j-2,k))
          prs(i,j,k)=num1d3*(4.d0*prs(i,j-1,k)-prs(i,j-2,k))
          vel(i,j,k,1)= num1d3*(4.d0*vel(i,j-1,k,1)-vel(i,j-2,k,1))

          ! if(x(i,j,k,1)>xs) then
          !   rho(i,j,k)  = ro1
          !   vel(i,j,k,1)= u1
          !   prs(i,j,k)  = p1
          ! else
          !   rho(i,j,k)  = ro2
          !   vel(i,j,k,1)= u2
          !   prs(i,j,k)  = p2
          ! endif

          vel(i,j,k,2)= 0.d0
          vel(i,j,k,3)= 0.d0

          tmp(i,j,k)  =thermal(pressure=prs(i,j,k),density=rho(i,j,k))

          call fvar2q(        q=  q(i,j,k,:),                      &
                        density=rho(i,j,k),                        &
                       velocity=vel(i,j,k,:),                      &
                       pressure=prs(i,j,k),                        &
                        species=spc(i,j,k,:)                       )
        enddo
        enddo

        if(irk==irkm) then
          i=im
          j=jm
          do k=0,km
          
            rho(i,j,k)=r11
            prs(i,j,k)=p11
            vel(i,j,k,2)= 0.d0
            vel(i,j,k,1)= 0.d0
            vel(i,j,k,3)= 0.d0
  
            tmp(i,j,k)  =thermal(pressure=prs(i,j,k),density=rho(i,j,k))
  
            call fvar2q(        q=  q(i,j,k,:),                      &
                          density=rho(i,j,k),                        &
                         velocity=vel(i,j,k,:),                      &
                         pressure=prs(i,j,k),                        &
                          species=spc(i,j,k,:)                       )
          enddo
        endif

    endif

  end subroutine udf_bc

  function shockspeed(p1,t1,p2) result(speed)
    
    use commvar, only : gamma
    use fludyna, only : sos

    real(8) :: speed
    real(8),intent(in) :: p1,t1,p2

    real(8) :: mx2,ss

    mx2=(p2/p1-1.d0)*(gamma+1.d0)/(2.d0*gamma)+1.d0

    ss=sos(t1)

    speed=sqrt(mx2)*ss

  end function shockspeed

end module userdefine
!+---------------------------------------------------------------------+
!| The end of the module userdefine.                                   |
!+---------------------------------------------------------------------+
