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
  real(8) :: flamethickness,hsource
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
    use thermchem,only : tranco,spcindex,mixture,convertxiyi
    use cantera 
    !
    real(8) :: cpe,miu,kama,cs,lref
    real(8) :: specr(num_species),dispec(num_species)
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
    specr(:)=0.d0
    specr(spcindex('H2'))=0.0173
    specr(spcindex('O2'))=0.2289
    specr(spcindex('N2'))=1.d0-sum(specr)
    !
    ! pinf=5.d0*pinf
    uinf=0.97d0
    vinf=0.d0
    winf=0.d0
    tinf=300.d0
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
                        ia,num_species
    use commarray,only: x,vel,rho,prs,spc,tmp,q
    use parallel, only: lio
    use fludyna,  only: thermal
    !
#ifdef COMB
    !
    use thermchem,only : tranco,spcindex,mixture,convertxiyi
    use cantera 
    !
    ! local data
    integer :: i,j,k
    real(8) ::  xc,yc,zc,tmpr,tmpp,xloc,xwid,specr(num_species),  &
      specp(num_species),arg,prgvar,masflx,specx(num_species)
    real(8) :: pthick
    !
    tmpr=300.d0
    xloc=3.d0*xmax/4.d0
    xwid=xmax/(12.d0*5.3d0*2.d0)
    !
    !reactants
    specr(:)=0.d0
    specr(spcindex('H2'))=0.0173
    specr(spcindex('O2'))=0.2289
    specr(spcindex('N2'))=1.d0-sum(specr)
    !
    !products
    tmpp=1814.32d0
    !
    ! pthick=1.d-4
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
      prgvar=1.d0*exp(-0.5d0*((xc-xloc)/xwid)**2)
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
    if(mode=='cuboid') then
      lx=12.d0*5.3d0*flamethickness
      ly=5.3d0*flamethickness
      lz=5.3d0*flamethickness
    elseif(mode=='cubic') then
      lx=5.3d0*flamethickness
      ly=5.3d0*flamethickness
      lz=5.3d0*flamethickness
    else
      stop ' !! error1 @ gridhitflame'
    endif
    !
    do k=0,km
    do j=0,jm
    do i=0,im
      x(i,j,k,1)=lx/real(ia,8)*real(i+ig0,8)
      if(ja==0) then
        x(i,j,k,2)=ly
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
    if(lio) print*,' ** cubic grid generated'
    !
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
    use parallel, only : pmax,psum,irk,irkm,lio
    use utility,  only : listinit,listwrite
    !
    use thermchem,only : heatrate
    !
    integer :: i,j,k
    real(8) :: tmpmax,rhomax,umax,qdotmax,poutrt
    real(8) :: qdot,var1,var2
    !
    integer,save :: hand_fs
    real(8),save :: xflame=0.d0,vflame=0.d0
    logical,save :: linit=.true.
    !
    if(lio) then
      !
      if(linit) then
        call listinit(filename='flamesta.dat',handle=hand_fs, &
                      firstline='nstep time maxT maxU maxHRR xflame vflame pout')
        linit=.false.
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
    do j=1,jm
      !
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
    enddo
    !
    var1=psum(var1)/ja
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
    ! calculate mean pressure at outflow
    i=im
    k=0
    poutrt=0.d0
    !
    if(irk==irkm) then
      do j=1,jm
        poutrt=poutrt+prs(i,j,k)
      enddo
    endif
    !
    poutrt=psum(poutrt)/dble(ja)
    !
    if(lio) call listwrite(hand_fs,tmpmax,umax,qdotmax,xflame,vflame,poutrt)
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
    use commvar,  only : im,jm,km,ndims,deltat,ia,ja,ka,rkstep,xmax,ymax,zmax
    use parallel, only : mpirank,psum,bcast
    use commarray,only : rho,tmp,vel,qrhs,x,jacob
    ! !
    ! ! local data
    integer,parameter :: nfan=5 !the number of fans
    !
    integer :: i,j,k,k1,k2,n
    real(8) :: force(3)
    logical,save :: linit=.true.
    real(8),save :: A(3,3,nfan),B(3,3,nfan)
    real(8) :: FC,FR,var1,var2,tavg,at,xs,xe,xx,yy,zz,lwave
    real(8) :: xs1,xs2,xs3,xs4,xs5
    integer,allocatable :: seed(:)
    !
    if(linit) then
      !
      if(mpirank==0) then
        !
        FR=1.d0/16.d0
        FC=0.d0
        !
        var1=sqrt(num2d3*FC/deltat)
        var2=sqrt(num1d3*FR/deltat)
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
        A=sqrt(3.d0)*(2.d0*A-1.d0)
        B=sqrt(3.d0)*(2.d0*B-1.d0)
        !
        do j=1,3
        do i=1,3
          if(i==j) then
            A(i,j,:)=var1*A(i,j,:)
            B(i,j,:)=var1*B(i,j,:)
          else
            A(i,j,:)=var2*A(i,j,:)
            B(i,j,:)=var2*B(i,j,:)
          endif
        end do
        end do
        !
      endif
      !
      call bcast(A)
      call bcast(B)
      !
      A=0.9594d0*A
      B=0.9594d0*B
      !
      linit=.false.
      !  
    endif
    !
    lwave=ymax
    !
    xs=1.1065d-2
    xe=xs+4.d0*lwave
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
        n=int((x(i,j,k,1)-xs)/lwave)+1
        !
        xx=(x(i,j,k,1)-xs)/lwave*2.d0*pi
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
    use commvar,  only: ymin,ymax,im,jm,filenumb
    use commarray,only: x,rho,vel,tmp,prs,spc
    use readwrite,only: writexprofile
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
    character(len=4) :: stepname
    !
    lwprofile=.false.
    ypos=0.5d0*(ymax-ymin)+ymin
    do j=1,jm
      if(x(0,j-1,0,2)<ypos .and. x(0,j,0,2)>=ypos) then
        !
        lwprofile=.true.
        !
        exit
        !
      endif
    enddo
    !
#ifdef COMB
    allocate(hrr(0:im))
    do i=0,im
      hrr(i)=heatrate(rho(i,0,0),tmp(i,0,0),spc(i,0,0,:))
    enddo
    !
    write(stepname,'(i4.4)')filenumb
    call writexprofile(profilename='outdat/profile'//trim(stepname)//'.dat',  &
                               var1=rho(0:im,j,0),  var1name='rho', &
                               var2=vel(0:im,j,0,1),var2name='u',   &
                               var3=tmp(0:im,j,0),  var3name='T',   &
                               var4=prs(0:im,j,0),  var4name='P',   &
                               var5=hrr(0:im),      var5name='HRR',truewrite=lwprofile)
#endif
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
