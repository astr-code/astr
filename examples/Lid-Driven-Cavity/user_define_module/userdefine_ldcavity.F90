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
    use commvar,  only : im,jm,km,ia,ja,ka,deltat
    use commarray,only : vel,prs,rho,tmp
    use parallel, only : pmax,lio
    use utility,  only : listinit,listwrite

    integer :: i,j,k

    integer,save :: hand_fs
    real(8) :: R_rho,R_u,R_v,R_T,R_p
    real(8) :: A_rho,A_u,A_v,A_T,A_p
    real(8),allocatable,dimension(:,:,:),save :: data_save
    logical,save :: linit=.true.

    if(linit) then

      if(lio) then
        call listinit(filename='residual.dat',handle=hand_fs, &
                      firstline='nstep time R_rho R_u R_v R_T R_p')
      endif

      allocate(data_save(0:im,0:jm,1:5))

      data_save=0.d0

      linit=.false.

    endif

    R_rho=0.d0
    R_u=0.d0
    R_v=0.d0
    R_T=0.d0
    R_p=0.d0

    A_rho=0.d0
    A_u=0.d0
    A_v=0.d0
    A_T=0.d0
    A_p=0.d0

    k=0
    do j=0,jm
    do i=0,im
      A_rho=max(A_rho,abs(rho(i,j,k)  ))
      A_u  =max(A_u,  abs(vel(i,j,k,1)))
      A_v  =max(A_v,  abs(vel(i,j,k,2)))
      A_T  =max(A_T,  abs(tmp(i,j,k)  ))
      A_p  =max(A_p,  abs(prs(i,j,k)  ))

      R_rho=max(R_rho,abs(rho(i,j,k)  -data_save(i,j,1)))
      R_u  =max(R_u,  abs(vel(i,j,k,1)-data_save(i,j,2)))
      R_v  =max(R_v,  abs(vel(i,j,k,2)-data_save(i,j,3)))
      R_T  =max(R_T,  abs(tmp(i,j,k)  -data_save(i,j,4)))
      R_p  =max(R_p,  abs(prs(i,j,k)  -data_save(i,j,5)))
    enddo 
    enddo
    R_rho=pmax(R_rho)/pmax(A_rho)
    R_u  =pmax(R_u  )/pmax(A_u  )
    R_v  =pmax(R_v  )/pmax(A_v  )
    R_T  =pmax(R_T  )/pmax(A_T  )
    R_p  =pmax(R_p  )/pmax(A_p  )

    data_save(0:im,0:jm,1)=rho(0:im,0:jm,0)
    data_save(0:im,0:jm,2)=vel(0:im,0:jm,0,1)
    data_save(0:im,0:jm,3)=vel(0:im,0:jm,0,2)
    data_save(0:im,0:jm,4)=tmp(0:im,0:jm,0)
    data_save(0:im,0:jm,5)=prs(0:im,0:jm,0)

    if(lio) call listwrite(hand_fs,R_rho,R_u,R_v,R_T,R_p)

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
    use commarray, only : prs,vel,tmp,rho,spc,q
    use commvar,   only : nondimen,im,km,jm,num_species
    use parallel,  only : jrk,jrkm
    use fludyna,   only : thermal,fvar2q

    ! arguments
    integer,intent(in) :: ndir

    ! local data
    integer :: i,j,k
    real(8) :: pe
    !

    if(ndir==4 .and. jrk==jrkm) then
      
        j=jm
        do k=0,km
        do i=0,im
          pe=num1d3*(4.d0*prs(i,j-1,k)-prs(i,j-2,k))
          !
          vel(i,j,k,1)=1.d0
          vel(i,j,k,2)=0.d0
          vel(i,j,k,3)=0.d0
          prs(i,j,k)  =pe
          tmp(i,j,k)  =1.d0
          !
          if(num_species>0) then
            spc(i,j,k,1)=num1d3*(4.d0*spc(i,j-1,k,1)-spc(i,j-2,k,1))
          endif
          !
          if(nondimen) then
            !
            rho(i,j,k)  =thermal(pressure=prs(i,j,k),temperature=tmp(i,j,k))
            !
            call fvar2q(      q=  q(i,j,k,:),                            &
                        density=rho(i,j,k),                              &
                      velocity=vel(i,j,k,:),                            &
                      pressure=prs(i,j,k),                              &
                        species=spc(i,j,k,:)                             )
            !
          else
            !
            rho(i,j,k)  =thermal(pressure=prs(i,j,k),temperature=tmp(i,j,k),species=spc(i,j,k,:))
            !
            call fvar2q(      q=  q(i,j,k,:),                            &
                        density=rho(i,j,k),                              &
                        velocity=vel(i,j,k,:),                           &
                        temperature=tmp(i,j,k),                          &
                        species=spc(i,j,k,:)                             )
            !
          endif
          !
        enddo
        enddo

    endif

  end subroutine udf_bc

end module userdefine
!+---------------------------------------------------------------------+
!| The end of the module userdefine.                                   |
!+---------------------------------------------------------------------+
