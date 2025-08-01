!+---------------------------------------------------------------------+
!| This module contains user defined subroutines to interfere  program |
!+---------------------------------------------------------------------+
!| CHANGE RECORD                                                       |
!| -------------                                                       |
!| 18-08-2023  | Created by J. Fang                                    |
!+---------------------------------------------------------------------+
module userdefine

  use constdef

  implicit none

  real(8) :: cavity_length,cavity_depth

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
    
    cavity_depth=2.5d0
    cavity_length=2.d0*cavity_depth

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

    use commvar,  only : im,jm,km,gridfile,ia,ja,ka
    use parallel, only : ig0,jg0,kg0,lio
    use commarray,only : x
    use parallel, only : mpirank
    !
    ! local data
    integer :: jn,js,in
    integer :: i,j,k,ii,jj,kk
    real(8) :: lx,xmax,lz
    real(8),allocatable :: yy(:),xx(:)
    
    allocate(yy(0:ja),xx(0:ia))

    jn=ja/4
    js=ja-40

    yy(0:jn)=cosmesh(jn=jn,dymin=4.d-3,ymax=cavity_depth)-cavity_depth
    yy(jn:js)=quatcosmesh(jn=js-jn,dymin=4.d-3,ymax=5.d0)

    call spongstretch_2(ja-js,10.d0,yy(js-2:ja))

    lx=20.d0+cavity_length+10.d0
    xmax=lx+10.d0

    in=ia-50

    do i=0,in
      xx(i)=lx/real(in,8)*real(i,8)
    enddo
    call spongstretch_2(ia-in,xmax,xx(in-2:ia))

    xx=xx-20.d0

    if(mpirank==0) then
      
      open(18,file='gridy.dat')
      do j=0,ja
        if(j==0) then
          write(18,*)j,yy(j),0.d0
        else
          write(18,*)j,yy(j),yy(j)-yy(j-1)
        endif
      end do
      close(18)
      print*,' << gridy.dat ... done.'
      
      open(18,file='gridx.dat')
      do i=1,ia
      write(18,*)i,xx(i),xx(i)-xx(i-1)
      end do
      close(18)
      print*,' <<< gridx.dat ... done.'

    endif

    lz=0.d0

    do k=0,km
    do j=0,jm
    do i=0,im
      
      ii=i+ig0
      jj=j+jg0
      kk=k+kg0

      x(i,j,k,1)=xx(ii)
      x(i,j,k,2)=yy(jj)

      if(ka==0) then
        x(i,j,k,3)=0.d0
      else
        x(i,j,k,3)=lz/real(ka,8)*real(k+kg0,8)
      endif

    enddo
    enddo
    enddo
    
    if(lio) print*,' ** cubic grid generated'

    ! if(lio) then
    !   open(18,file='tecsolid.plt')
    !   write(18,'(A)')'VARIABLES = "x" "y" '
    !   write(18,'(A)')'ZONE T="solid body 1"'
    !   write(18,'(A)')' I=2, J=2, K=1, ZONETYPE=Ordered DATAPACKING=POINT'
    !   write(18,*)xx(0), -1.d0*cavity_depth
    !   write(18,*)0.d0, -1.d0*cavity_depth
    !   write(18,*)xx(0), 0.d0
    !   write(18,*)0.d0, 0.d0
    !   write(18,'(A)')'ZONE T="solid body 2"'
    !   write(18,'(A)')' I=2, J=2, K=1, ZONETYPE=Ordered DATAPACKING=POINT'
    !   write(18,*)cavity_length, -1.d0*cavity_depth
    !   write(18,*)xx(ia), -1.d0*cavity_depth
    !   write(18,*)cavity_length, 0.d0
    !   write(18,*)xx(ia), 0.d0
    !   close(18)
    !   print*,' << tecsolid.plt'
    ! endif

  end subroutine udf_grid
  !+-------------------------------------------------------------------+
  !| The end of the subroutine udf_grid.                               |
  !+-------------------------------------------------------------------+
  
  function cosmesh(jn,dymin,dymax,ymax) result(y)
    !
    integer,intent(in) :: jn
    real(8),intent(in),optional :: dymin,dymax,ymax
    real(8) :: y(0:jn)
    !
    integer :: j
    real(8) :: var1,dymaxtemp,ymaxtemp,error,ddy,er
    !
    if(present(dymax) .and. present(dymin)) then
      y(0)=0.d0
      do j=1,jn
        var1=2.d0*pi/dble(jn-1)*dble(j-1)
        var1=0.5d0*(dymax-dymin)*(-cos(var1)+1.d0)+dymin
        y(j)=y(j-1)+var1
      enddo
      !
      print*,' ** dymin=',dymin,'dymax=',dymax,'ymax=',y(jn)
      !
    elseif(present(ymax) .and. present(dymin)) then
      ! first guess dymaxtemp
      ddy=0.5d0*dymin
      dymaxtemp=dymin+ddy
      !
      error=1.d10
      do while(abs(error)>1.d-12)
        !
        ymaxtemp=cosmesh_ymax(dymaxtemp)
        !
        if(ymaxtemp<ymax) then
          ! dymax is too small
          dymaxtemp=dymaxtemp+ddy
        else
          dymaxtemp=dymaxtemp-ddy
        endif
        ! !
        er=(ymaxtemp-ymax)
        !
        if(er*error<0) ddy=0.5d0*ddy
        !
        error=er
        !
        ! print*,error,ymaxtemp,dymaxtemp
        !
      enddo
      !
      y(0)=0.d0
      do j=1,jn
        var1=2.d0*pi/dble(jn-1)*dble(j-1)
        var1=0.5d0*(dymaxtemp-dymin)*(-cos(var1)+1.d0)+dymin
        y(j)=y(j-1)+var1
      enddo
      !
      print*,' ** dymin=',dymin,'dymax=',dymaxtemp,'ymax=',y(jn)
      !
    endif
    !
    contains
    !
    function cosmesh_ymax(dymax1) result(ymax1)
      !
      real(8),intent(in) :: dymax1
      real(8) :: ymax1
      !
      real(8) :: varc1
      integer :: j
      !
      ymax1=0.d0
      do j=1,jn
        varc1=2.d0*pi/dble(jn-1)*dble(j-1)
        varc1=0.5d0*(dymax1-dymin)*(-cos(varc1)+1.d0)+dymin
        ymax1=ymax1+varc1
      enddo
      !
      return
      !
    end function cosmesh_ymax
    !
  end function cosmesh

  function quatcosmesh(jn,dymin,dymax,ymax) result(y)
    !
    integer,intent(in) :: jn
    real(8),intent(in),optional :: dymin,dymax,ymax
    real(8) :: y(0:jn)
    !
    integer :: j
    real(8) :: var1,var2,dymaxtemp,ymaxtemp,error,ddy,er
    real(8) :: dymintemp,dymin2
    !
    if(present(dymax) .and. present(dymin)) then
      !
      y(0)=0.d0
      do j=1,jn
        var1=0.5d0*pi/dble(jn-1)*dble(j-1)
        var1=0.5d0*(dymax-dymin)*(-cos(var1)+1.d0)+dymin
        y(j)=y(j-1)+var1
      enddo
      !
      print*,' ** dymin=',dymin,'dymax=',dymax,'ymax=',y(jn)
      !
    elseif(present(ymax) .and. present(dymin)) then
      ! first guess dymaxtemp
      ddy=0.5d0*dymin
      dymaxtemp=dymin+ddy
      !
      error=1.d10
      do while(abs(error)>1.d-12)
        !
        ymaxtemp=quatcosmesh_ymax(dymaxtemp)
        !
        if(ymaxtemp<ymax) then
          ! dymax is too small
          dymaxtemp=dymaxtemp+ddy
        else
          dymaxtemp=dymaxtemp-ddy
        endif
        ! !
        er=(ymaxtemp-ymax)
        !
        if(er*error<0) ddy=0.5d0*ddy
        !
        error=er
        !
        ! print*,error,ymaxtemp,dymaxtemp
        !
      enddo
      !
      y(0)=0.d0
      do j=1,jn
        var1=0.5d0*pi/dble(jn-1)*dble(j-1)
        var1=0.5d0*(dymaxtemp-dymin)*(-cos(var1)+1.d0)+dymin
        y(j)=y(j-1)+var1
      enddo
      !
      print*,' ** dymin=',dymin,'dymax=',dymaxtemp,'ymax=',y(jn)
      !
    elseif(present(ymax) .and. present(dymax)) then
        !
        if( dymax*dble(jn) < ymax) then
          !
          print*,'dymax=',dymax,'jn=',jn,'dymax*jn=',dymax*dble(jn),'ymax=',ymax
          stop ' !! error dymax too small' 
          !
        endif
        !
        var2=0.d0
        do j=1,jn
          var1=0.5d0*pi/dble(jn-1)*dble(j-1)
          var1=0.5d0*(dymax-0.d0)*(cos(var1)+1.d0)+0.d0
          var2=var2+var1
          print*,j,var1
        enddo
        !
        if(var2<=ymax) then
          print*,' minimal ymax can reached:',var2,' but larger than ymax required:',ymax
          stop ' !! increase the number ofnodes' 
        endif
        !
        print*, ' ymax range:',dymax*dble(jn),'~',var2
        print*, ' ymax   set:',ymax,'jn:',jn,'dymax:',dymax
        stop
        !
        dymintemp=0.5d0*dymax
        error=10.d0
        !
        do while(abs(error)>1.d-8)
          var2=0.d0
          ddy=0.5d0*dymintemp
          !
          do j=1,jn
            var1=0.5d0*pi/dble(jn-1)*dble(j-1)
            var1=0.5d0*(dymax-dymintemp)*(cos(var1)+1.d0)+dymintemp
            var2=var2+var1
          enddo
          !
          stop
          !
          if(error*(var2-ymax)<0) ddy=0.5d0*ddy
          !
          if(error>0) then
              dymintemp=dymintemp-ddy
          else
              dymintemp=dymintemp+ddy
          endif
          !
          error=var2-ymax
          !
          print*,error,dymintemp,var2/ymax
          !
        enddo
        !
        if( dymintemp <= 0.d0) then
          !
          print*,dymintemp
          stop ' !! error negative dymin' 
          !
        endif
        !
        y(0)=0.d0
        do j=1,jn
          var1=0.5d0*pi/dble(jn-1)*dble(j-1)
          var1=0.5d0*(dymax-dymintemp)*(cos(var1)+1.d0)+dymintemp
          y(j)=y(j-1)+var1
        enddo
        !
        print*,' ** dymin=',dymintemp,'dymax=',dymax,'ymax=',y(jn)
        !
    endif
    !
    contains
    !
    function quatcosmesh_ymax(dymax1) result(ymax1)
      !
      real(8),intent(in) :: dymax1
      real(8) :: ymax1
      !
      real(8) :: varc1
      integer :: j
      !
      ymax1=0.d0
      do j=1,jn
        varc1=0.5d0*pi/dble(jn-1)*dble(j-1)
        varc1=0.5d0*(dymax1-dymin)*(-cos(varc1)+1.d0)+dymin
        ymax1=ymax1+varc1
      enddo
      !
      return
      !
    end function quatcosmesh_ymax
    !
  end function quatcosmesh

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! This subroutine is used to stretch grid for spong layers
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine spongstretch_2(dim,varmax,var)
    !
    integer :: dim,hh
    real(8) :: var(-2:dim)
    real(8) :: varmax
    real(8) :: err,err2,dra,ddra,var1,var2,dde,ratio
    !
    dra=1.d0
    !
    ddra=0.001d0*dra
    err=1.d10
    err2=0.d0
    !
    do while( dabs(err)>1.d-9 )
      !
      ! ratio=1.d0
      ratio=(var(0)-var(-1))/(var(-1)-var(-2))
      do hh=1,dim
        !
        dde=var(hh-1)-var(hh-2)
        var(hh)=var(hh-1)+dde*ratio
        ratio=ratio*dra
        !
      end do
      !
      err=varmax-var(dim)
      !
      if(err>0) then
        dra=dra+ddra
      else
        dra=dra-ddra
      end if
      !
      if(err*err2<0.d0) ddra=0.5d0*ddra
      !
      err2=err
      !
      !print*,err,ddra,dra,ratio
      !
    end do
    !
    ! print*,' the mesh is strenched to x=',varmax,' from x=',var(0),' the max ratio is:',ratio,' using ',dim,' nodes'
    !
  end subroutine spongstretch_2
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! End of the subroutine spongstretch_2
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
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

  !+-------------------------------------------------------------------+
  !| This subroutine is a user defined boundary condition              |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 02-Juil-2025: created by Jian Fang @ IMech, CAS, Beijing          |
  !+-------------------------------------------------------------------+
  subroutine udf_bc(ndir)
    
    ! arguments
    integer,intent(in) :: ndir

  end subroutine udf_bc
  !+-------------------------------------------------------------------+
  !| The end of the subroutine udf_bc.                                 |
  !+-------------------------------------------------------------------+

  !+-------------------------------------------------------------------+
  !| This subroutine is a user defined immersed solid body.            | 
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 01-Aug-2025: created by Jian Fang @ IMech, CAS, Beijing           |
  !+-------------------------------------------------------------------+
  subroutine udf_immersed_solid(xp,inside,bnode)
    !
    use commtype,  only : sboun
    !
    real(8),intent(in) :: xp(3)
    logical,intent(out),optional :: inside
    type(sboun),intent(out),optional :: bnode
    !
    real(8) :: dis1,dis2
    real(8) :: step_upper_y,step_right_x,step_left_x

    step_left_x=0.d0
    step_right_x=step_left_x+cavity_length
    step_upper_y=0.d0

    if(xp(2)<step_upper_y .and. xp(1)<=step_left_x) then
      !
      if(present(inside)) inside=.true.
      !
      if(present(bnode)) then
        !
        dis1=abs(step_upper_y-xp(2))
        dis2=abs(step_left_x-xp(1))
        !
        if(dis2<dis1) then
          bnode%x(1)=step_left_x
          bnode%x(2)=xp(2)
          bnode%x(3)=xp(3)
          !
          bnode%normdir(1)=1.d0
          bnode%normdir(2)=0.d0
          bnode%normdir(3)=0.d0
        else
          bnode%x(1)=xp(1)
          bnode%x(2)=step_upper_y
          bnode%x(3)=xp(3)
          !
          bnode%normdir(1)=0.d0
          bnode%normdir(2)=1.d0
          bnode%normdir(3)=0.d0
        endif
        !
      endif
      !
    elseif(xp(2)<step_upper_y .and. xp(1)>=step_right_x) then
      !
      if(present(inside)) inside=.true.
      !
      if(present(bnode)) then
        !
        dis1=abs(step_upper_y-xp(2))
        dis2=abs(step_right_x-xp(1))
        !
        if(dis2<dis1) then
          bnode%x(1)=step_right_x
          bnode%x(2)=xp(2)
          bnode%x(3)=xp(3)
          !
          bnode%normdir(1)=-1.d0
          bnode%normdir(2)=0.d0
          bnode%normdir(3)=0.d0
        else
          bnode%x(1)=xp(1)
          bnode%x(2)=step_upper_y
          bnode%x(3)=xp(3)
          !
          bnode%normdir(1)=0.d0
          bnode%normdir(2)=1.d0
          bnode%normdir(3)=0.d0
        endif
        !
      endif
      !
    else
      if(present(inside)) inside=.false.
    endif

    return

  end subroutine udf_immersed_solid
  !+-------------------------------------------------------------------+
  !| The end of the subroutine udf_immersed_solid.                     |
  !+-------------------------------------------------------------------+

end module userdefine
!+---------------------------------------------------------------------+
!| The end of the module userdefine.                                   |
!+---------------------------------------------------------------------+
