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
  real(8),allocatable,dimension(:,:) :: denu1f,denu2f,denu3f
  real(8),allocatable,dimension(:,:) :: res
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
    use commvar,  only: jm
    use parallel, only: preadprofile
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
     allocate(res(0:jm,1:4))
     !
     call preadprofile('datin/inlet.res',dir='j',                     &
                               var1=res(:,1),     var2=res(:,2), &
                               var3=res(:,3),     var4=res(:,4),skipline=1)
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
  !| 05-Oct-2023: isotropic digitial filter for mixing layer           |
  !+-------------------------------------------------------------------+
  subroutine udf_inflow_fluc(umean,uinst)
    !
    use hdf5io
    use constdef
    use parallel,only : lio,mpistop,mpirank,jg0,kg0,mpi_imin,mpirankname
    use commvar, only : nstep,deltat,jm,km,ja,ka
    use commarray,only : x
    use tecio
    use utility, only : rnorm_box_muller
    !
    real(8),intent(in) ::  umean(0:jm,1:3)  ! inflow mean velocity
    real(8),intent(out) :: uinst(0:jm,0:km,1:3)  ! velocity with fluctuations
    !
    integer :: ntfjmax,ntfkmax
    !
    real(8),allocatable,save :: btfj1(:),btfk1(:)
    real(8),allocatable,dimension(:,:,:) :: random_array
    real(8),allocatable,dimension(:,:) :: ranf1,ranf2,ranf3
    !
    integer :: i,j,j1,j2,k,k1,k2,n
    real(8) :: var1,var2,var3,var4,var5,var6,var7,varm,taw,ss1,ss2,vran(2)
    real(8) :: u1f,u2f,u3f
    logical,save :: lfirstcal=.true.
    !
    ntfjmax=32
    ntfkmax=26
    ! 
    if(lfirstcal) then
      ! 
      allocate( btfj1(-ntfjmax:ntfjmax),btfk1(-ntfjmax:ntfjmax))
      !
      btfj1=0.d0
      var1=0.d0
      do n=-ntfjmax,ntfjmax
        btfj1(n)=exp(-2.d0*pi*abs(n)/ntfjmax)
        var1=var1+btfj1(n)**2
      end do
      var1=sqrt(var1)
      btfj1=btfj1/var1
      !
      btfk1=0.d0
      var1=0.d0
      do n=-ntfkmax,ntfkmax
        btfk1(n)=exp(-2.d0*pi*abs(n)/ntfkmax)
        var1=var1+btfk1(n)**2
      end do
      var1=sqrt(var1)
      btfk1=btfk1/var1
      !
      allocate( denu1f(0:jm,0:km),denu2f(0:jm,0:km),denu3f(0:jm,0:km)    )
      !
      if(nstep>1) then
        !
        call h5io_init('outdat/denuf.h5',mode='read',comm=mpi_imin)
        !
        call h5read(varname='denu1f',var=denu1f,dir='i')
        call h5read(varname='denu2f',var=denu2f,dir='i')
        call h5read(varname='denu3f',var=denu3f,dir='i')
        !
        call h5io_end
        !
      endif
      !
      lfirstcal=.false.
      !
      if(lio) print*,' ** udf_inflow_fluc initialised'
      !
    end if
    !
    allocate(random_array(-ntfjmax:ja+ntfjmax,-ntfkmax:ka+ntfkmax,1:3))
    !
    ! get random numbers with normal distributions 
    do k=-ntfkmax,ka+ntfkmax
    do j=-ntfjmax,ja+ntfjmax
      !
      vran=rnorm_box_muller(mode='sync')
      random_array(j,k,1:2)=vran(1:2)
      vran=rnorm_box_muller(mode='sync')
      random_array(j,k,3)=vran(1)
      !
    enddo
    enddo
    !
    !
    ! Filter the random number to get random variables with structures
    allocate(ranf1(0:jm,0:km),ranf2(0:jm,0:km),ranf3(0:jm,0:km))
    do k=0,km
    do j=0,jm
      !
      var1=0.0
      var2=0.0
      var3=0.0
      do k1=-ntfkmax,ntfkmax
      do j1=-ntfjmax,ntfjmax
        !
        j2=j+j1+jg0
        k2=k+k1+kg0
        !
        var1=var1+btfj1(j1)*btfk1(k1)*random_array(j2,k2,1)
        var2=var2+btfj1(j1)*btfk1(k1)*random_array(j2,k2,2)
        var3=var3+btfj1(j1)*btfk1(k1)*random_array(j2,k2,3)
        !
      end do
      end do
      !
      ranf1(j,k)=var1
      ranf2(j,k)=var2
      ranf3(j,k)=var3
      !
    enddo
    enddo
    !
    if(nstep==0) then
      !
      do k=0,km
      do j=0,jm
        denu1f(j,k)=ranf1(j,k)
        denu2f(j,k)=ranf2(j,k)
        denu3f(j,k)=ranf3(j,k)
      end do
      end do
      !
    else
      !
      do k=0,km
      do j=0,jm
        !
        taw=0.7d0*0.25d0/umean(j,1)
        !
        denu1f(j,k)=denu1f(j,k)*exp(-0.5d0*pi*deltat/taw)+           &
                    ranf1(j,k)*sqrt(1.d0-exp(-pi*deltat/taw))
        denu2f(j,k)=denu2f(j,k)*exp(-0.5d0*pi*deltat/taw)+           &
                    ranf2(j,k)*sqrt(1.d0-exp(-pi*deltat/taw))
        denu3f(j,k)=denu3f(j,k)*exp(-0.5d0*pi*deltat/taw)+           &
                    ranf3(j,k)*sqrt(1.d0-exp(-pi*deltat/taw))
      end do
      end do
      !
    end if
    !
    do j=0,jm
      !
      if(res(j,1)<1.d-16) then
        var1=0.d0
        var2=0.d0
        var3=0.d0
      else
        var1=sqrt(res(j,1))
        var2=res(j,4)/var1
        var3=sqrt(res(j,3))
      endif
      !
      do k=0,km
       !
       u1f=var1*denu1f(j,k)
       u2f=var2*denu1f(j,k)+sqrt(abs(res(j,2)-var2**2))*denu2f(j,k)
       u3f=var3*denu3f(j,k)
       !
       uinst(j,k,1)=umean(j,1)+u1f
       uinst(j,k,2)=umean(j,2)+u2f
       uinst(j,k,3)=u3f
       !
      end do
      !
    end do
    !
    ! call tecbin('testout/tecinfluc'//mpirankname//'.plt',          &
    !                                   x(0,0:jm,0:km,1),'x',        &
    !                                   x(0,0:jm,0:km,2),'y',        &
    !                                   x(0,0:jm,0:km,3),'z',        &
    !                           denu1f,'denu1f',denu2f,'denu2f',     &
    !                           denu3f,'denu3f',uinst(:,:,1),'u',    &
    !                           uinst(:,:,2),'v',uinst(:,:,3),'w')
    !
    deallocate(random_array,ranf1,ranf2,ranf3)
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
    use commarray,only : vel,prs,rho,tmp,spc,cell,x
    use interp,   only : interlinear
    use parallel, only : pmax,psum,irk,jrk,irkm,jrkm,lio
    use utility,  only : listinit,listwrite
!     !
! #ifdef COMB
!     use thermchem,only : heatrate,spcindex
! #endif
!     !
    integer :: i,j,k
    real(8) :: p_i0,p_j0,p_jm
!     real(8) :: tmpmax,rhomax,umax,qdotmax,poutrt
!     real(8) :: qdot,var1,var2
!     !
    integer,save :: hand_fs
!     real(8),save :: xflame=0.d0,vflame=0.d0
    logical,save :: linit=.true.
!     !
    if(lio) then
      !
      if(linit) then
        call listinit(filename='pressure_nscbc.dat',handle=hand_fs, &
                      firstline='nstep time p_i0 p_j0 p_jm')
        linit=.false.
      endif
      !
    endif
    !
    ! calculate mean pressure at outflow
    p_i0=0.d0
    p_j0=0.d0
    p_jm=0.d0
    !
    if(irk==irkm) then
      !
      i=im
      k=0
      !
      do j=1,jm
        p_i0=p_i0+prs(i,j,k)
      enddo
      !
    endif
    !
    if(jrk==0) then
      j=0
      k=0
      do i=1,im
        p_j0=p_j0+prs(i,j,k)
      enddo
    endif
    !
    if(jrk==jrkm) then
      j=jm
      k=0
      do i=1,im
        p_jm=p_jm+prs(i,j,k)
      enddo
    endif
    !
    p_i0=psum(p_i0)/dble(ja)
    p_j0=psum(p_j0)/dble(ia)
    p_jm=psum(p_jm)/dble(ia)
    !
    if(lio) call listwrite(hand_fs,p_i0,p_j0,p_jm)
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
  !| This subroutine is to defined an output by a user.                | 
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 18-Aug-2023: created by Jian Fang @ Daresbury                     |
  !+-------------------------------------------------------------------+
  subroutine udf_write
    !
    use parallel, only: mpi_imin
    use hdf5io
    !
    if(irk==0) then
      !
      call h5io_init('outdat/denuf.h5',mode='write',comm=mpi_imin)
      !
      call h5write(varname='denu1f',var=denu1f,dir='i')
      call h5write(varname='denu2f',var=denu2f,dir='i')
      call h5write(varname='denu3f',var=denu3f,dir='i')
      !
      call h5io_end
      !
    endif
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
