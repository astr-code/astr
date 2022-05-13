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
  use constdef
  use parallel,only: lio,mpistop,mpirank,mpirankname,bcast
  use commvar, only: im,jm,km,uinf,vinf,winf,pinf,roinf,tinf,ndims,    &
                     num_species,xmin,xmax,ymin,ymax,zmin,zmax,spcinf
  use tecio
  !
  implicit none
  !
  real(8) :: nomi_thick,disp_thick,mome_thick,fric_velocity
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
    use commvar,  only: flowtype,nstep,time,filenumb,fnumslic,ninit,   &
                        lrestart,lavg,turbmode
    use commarray,only: vel,rho,prs,spc,q,tke,omg
    use readwrite,only: readcont,readflowini3d,readflowini2d,          &
                        readcheckpoint,readmeanflow,readmonc,writeflfed
    use fludyna,  only: updateq
    use statistic,only: nsamples
    use bc,       only: ninflowslice
    !
    if(trim(flowtype)=='bl' .or. trim(flowtype)=='swbli') then
      call blprofile
    endif
    !
    call readcont
    !
    if(lrestart) then
      !
      call readcheckpoint(folder='checkpoint',mode='h')
      !
      call updateq
      !
    else
      !
      if(ninit==3) then
        !
        call readflowini3d
        !
        ! call blcorrect
        ! vel(:,:,:,2)=0.d0
        !
        if(trim(turbmode)=='k-omega') then
          tke=0.d0
          omg=0.d0
        endif
        !
      elseif(ninit==2) then
        !
        call readflowini2d
        !
      else
        !
        select case(trim(flowtype))
        case('2dvort')
          call vortini
        case('channel')
          call chanini
        case('tgv')
          call tgvini
        case('jet')
          call jetini
        case('accutest')
          call accini
        case('cylinder')
          call cylinderini
        case('mixlayer')
          call mixlayerini
        case('shuosher')
          call shuosherini
        case('bl')
          call blini
        case('swbli')
          call blini
        case('windtunn')
          call wtini
        case('0dreactor')
          call reactorini
        case('1dflame')
          call onedflameini
        case('h2supersonic')
          call h2supersonicini
        case('tgvflame')
          call tgvflameini
        case('rti')
          call rtini
        case default
          print*,trim(flowtype)
          stop ' !! flowtype not defined @ flowinit'
        end select
        !
      endif
      !
      call updateq
      !
      nstep=0
      ninflowslice=3
      time=0.d0
      !
      filenumb=0
      fnumslic=0
      !
    endif
    !
    call readmonc
    !
    if(lavg) then
      !
      if(nsamples>0) then
        call readmeanflow(mode='h')
      endif
      !
    endif
    !
    if(lio) print*,' ** flowfield initialised.'
    !
    ! call writeflfed
    ! stop
    !
  end subroutine flowinit
  !+-------------------------------------------------------------------+
  !| The end of the subroutine flowinit.                               |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to initilise sponge layer.                     |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 11-03-2021  | Created by J. Fang @ Warrington                     |
  !|             | (have not consider subdomain situation)             |
  !+-------------------------------------------------------------------+
  subroutine spongelayerini
    !
    use commvar, only : spg_imin,spg_imax,spg_jmin,spg_jmax,           &
                        spg_kmin,spg_kmax,im,jm,km,ia,ja,ka,           &
                        lisponge,ljsponge,lksponge
    use commarray,only: lspg_imin,lspg_imax,lspg_jmin,lspg_jmax,       &
                        lspg_kmin,lspg_kmax,x
    use parallel,only : ig0,jg0,kg0
    !
    ! local data
    integer :: i,j,k
    !
    if(spg_imin>0 .or. spg_imax>0) then
      lisponge=.true.
    else
      lisponge=.false.
    endif
    !
    if(spg_jmin>0 .or. spg_jmax>0) then
      ljsponge=.true.
    else
      ljsponge=.false.
    endif
    !
    if(spg_kmin>0 .or. spg_kmax>0) then
      lksponge=.true.
    else
      lksponge=.false.
    endif
    !
    if(mpirank==0) then
      write(*,'(2X,62A)')('-',i=1,62)
      write(*,'(2X,A)')'                        *** sponge layer ***'
      if(lisponge) then
        write(*,'(22X,3(A,I4))')'i direction:   0 ~',spg_imin,         &
                                        ' ....... ',ia-spg_imax,' ~ ',ia
      endif
      if(ljsponge) then
        write(*,'(22X,3(A,I4))')'j direction:   0 ~',spg_jmin,         &
                                        ' ....... ',ja-spg_jmax,' ~ ',ja
      endif
      if(lksponge) then
        write(*,'(22X,3(A,I4))')'k direction:   0 ~',spg_kmin,         &
                                        ' ....... ',ka-spg_kmax,' ~ ',ka
      endif
      write(*,'(2X,62A)')('-',i=1,62)
    endif
    !
    if(spg_imin>0) then
      !
      if(spg_imin<ig0) then
        spg_imin=0
      elseif(spg_imin>ig0+im) then
        spg_imin=im
      else
        spg_imin=spg_imin-ig0
      endif
      !
    endif
    !
    allocate(lspg_imin(0:jm,0:km))
    lspg_imin=0.d0
    if(spg_imin>0) then
      do k=0,km
      do j=0,jm
        !
        do i=1,spg_imin
          lspg_imin(j,k)=lspg_imin(j,k)+                               &
                                  sqrt((x(i,j,k,1)-x(i-1,j,k,1))**2+   &
                                       (x(i,j,k,2)-x(i-1,j,k,2))**2+   &
                                       (x(i,j,k,3)-x(i-1,j,k,3))**2)
        enddo
        !
      enddo
      enddo
    endif
    !
    if(spg_imax>0) then
      !
      if(ia-spg_imax>ig0+im) then
        spg_imax=0
      elseif(ia-spg_imax<ig0) then
        spg_imax=im
      else
        spg_imax=ig0+im-(ia-spg_imax)
      endif
      !
    endif
    !
    allocate(lspg_imax(0:jm,0:km))
    lspg_imax=0.d0
    if(spg_imax>0) then
      do k=0,km
      do j=0,jm
        !
        do i=im-spg_imax,im-1
          lspg_imax(j,k)=lspg_imax(j,k)+                               &
                                  sqrt((x(i+1,j,k,1)-x(i,j,k,1))**2+   &
                                       (x(i+1,j,k,2)-x(i,j,k,2))**2+   &
                                       (x(i+1,j,k,3)-x(i,j,k,3))**2)
        enddo
        !
      enddo
      enddo
    endif
    !
    if(spg_jmin>0) then
      !
      if(spg_jmin<jg0) then
        spg_jmin=0
      elseif(spg_jmin>jg0+jm) then
        spg_jmin=jm
      else
        spg_jmin=spg_jmin-jg0
      endif
      !
    endif
    !
    allocate(lspg_jmin(0:im,0:km))
    lspg_jmin=0.d0
    if(spg_jmin>0) then
      do k=0,km
      do i=0,im
        !
        do j=1,spg_jmin
          lspg_jmin(i,k)=lspg_jmin(i,k)+                               &
                                  sqrt((x(i,j,k,1)-x(i,j-1,k,1))**2+   &
                                       (x(i,j,k,2)-x(i,j-1,k,2))**2+   &
                                       (x(i,j,k,3)-x(i,j-1,k,3))**2)
        enddo
        !
      enddo
      enddo
      !
    endif
    !
    if(spg_jmax>0) then
      !
      if(ja-spg_jmax>jg0+jm) then
        spg_jmax=0
      elseif(ja-spg_jmax<jg0) then
        spg_jmax=jm
      else
        spg_jmax=jg0+jm-(ja-spg_jmax)
      endif
      !
    endif
    !
    allocate(lspg_jmax(0:im,0:km))
    lspg_jmax=0.d0
    if(spg_jmax>0) then
      do k=0,km
      do i=0,im
        !
        do j=jm-spg_jmax,jm-1
          lspg_jmax(i,k)=lspg_jmax(i,k)+                               &
                                  sqrt((x(i,j+1,k,1)-x(i,j,k,1))**2+   &
                                       (x(i,j+1,k,2)-x(i,j,k,2))**2+   &
                                       (x(i,j+1,k,3)-x(i,j,k,3))**2)
        enddo
        !
      enddo
      enddo
    endif
    !
    if(spg_kmin>0) then
      !
      if(spg_kmin<kg0) then
        spg_kmin=0
      elseif(spg_kmin>kg0+km) then
        spg_kmin=km
      else
        spg_kmin=spg_kmin-kg0
      endif
      !
    endif
    !
    allocate(lspg_kmin(0:im,0:jm))
    lspg_kmin=0.d0
    if(spg_kmin>0) then
      do j=0,jm
      do i=0,im
        !
        do k=1,spg_kmin
          lspg_kmin(i,j)=lspg_kmin(i,j)+                               &
                                  sqrt((x(i,j,k,1)-x(i,j,k-1,1))**2+   &
                                       (x(i,j,k,2)-x(i,j,k-1,2))**2+   &
                                       (x(i,j,k,3)-x(i,j,k-1,3))**2)
        enddo
        !
      enddo
      enddo
    endif
    !
    if(spg_kmax>0) then
      !
      if(ka-spg_kmax>kg0+km) then
        spg_kmax=0
      elseif(ka-spg_kmax<kg0) then
        spg_kmax=km
      else
        spg_kmax=kg0+km-(ka-spg_kmax)
      endif
      !
    endif
    !
    allocate(lspg_kmax(0:im,0:jm))
    lspg_kmax=0.d0
    if(spg_kmax>0) then
      do j=0,jm
      do i=0,im
        !
        do k=1,spg_kmax
          lspg_kmax(i,j)=lspg_kmax(i,j)+                               &
                                  sqrt((x(i,j,k,1)-x(i,j,k-1,1))**2+   &
                                       (x(i,j,k,2)-x(i,j,k-1,2))**2+   &
                                       (x(i,j,k,3)-x(i,j,k-1,3))**2)
        enddo
        !
      enddo
      enddo
    endif
    !
    ! print*,mpirank,'|',lspg_imin(0,0),lspg_imax(0,0),'|',            &
    !                    lspg_jmin(0,0),lspg_jmax(0,0)
    ! !
    ! call mpistop
    !
  end subroutine spongelayerini
  !+-------------------------------------------------------------------+
  !| The end of the subroutine spongelayerini.                         |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to generate an initial field test of      |
  !| accuracy.                                                         |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 01-Mar-2021: Created by J. Fang @ STFC Daresbury Laboratory       |
  !+-------------------------------------------------------------------+
  subroutine accini
    !
    use commarray,only: x,vel,rho,prs,spc,tmp,q,acctest_ref
    use fludyna,  only: thermal
    !
    ! local data
    integer :: i,j,k,l
    real(8) :: xc,yc,zc,radi2,rvor,cvor,var1
    real(8) :: A,slow_wl,fast_wl,damper,randomv(32),gz,zl
    !
    ! if(ndims==1) then
      !
      xc=0.5d0*(xmax-xmin)
      yc=0.5d0*(ymax-ymin)
      zc=0.5d0*(zmax-zmin)
      !
      rvor=0.5d0 
      cvor=0.1d0*rvor
      !
      A = 1.d0
      slow_wl = 4.d0
      fast_wl = 48.d0
      !
      allocate(acctest_ref(0:im))
      !
      if(mpirank==0) then
        do l=1,32
          call random_number(randomv(l))
        end do
      endif
      !
      call bcast(randomv)
      !
      do k=0,km
      do j=0,jm
      do i=0,im
        !
        ! radi2=((x(i,j,k,1)-xc)**2+(x(i,j,k,2)-yc)**2)/rvor/rvor
        ! var1=cvor/rvor/rvor*exp(-0.5d0*radi2)
        radi2=(x(i,j,k,1)-xc)**2
        var1=exp(-5.d0*radi2)
        !
        rho(i,j,k)  =roinf
        vel(i,j,k,1)=1.d0
        vel(i,j,k,2)=0.d0
        vel(i,j,k,3)=0.d0
        prs(i,j,k)  =pinf
        !
        tmp(i,j,k)=thermal(density=rho(i,j,k),pressure=prs(i,j,k))
        !
        ! spc(i,j,k,1)=exp(-0.5d0*radi2)
        ! spc(i,j,k,1)=cos(1.6d0*pi*(x(i,j,k,1)-xc))*var1
        ! spc(i,j,k,1)=cos(1.6d0*pi*(x(i,j,k,1)-xc))*var1
        damper=exp(-0.5d0*(x(i,j,k,1)-xc)**2)
        !
        gz=0.d0
        do l=1,32
          !
          if(l==1) then
            zl=0.2d0/(1.d0-0.95d0**32)
          else
            zl=zl*0.95d0
          end if
          !
          gz=gz+zl*dsin(2.d0*pi*l*(x(i,j,k,1)/(xmax-xmin)+randomv(l)))
          !
        end do
        !
        ! var1= A*cos(abs(sqrt(x(i,j,k,1)-xc)**3)*slow_wl)
        acctest_ref(i)=gz !*damper
        spc(i,j,k,1)=acctest_ref(i)
        !
      enddo
      enddo
      enddo
    !
    ! else
    !   stop ' !! error @ accini'
    ! endif
    !
    !
    call tecbin('testout/tecinit'//mpirankname//'.plt',                &
                                      x(0:im,0:jm,0:km,1),'x',         &
                                      x(0:im,0:jm,0:km,2),'y',         &
                                      x(0:im,0:jm,0:km,3),'z',         &
                                      q(0:im,0:jm,0:km,1),'q1',        &
                                      q(0:im,0:jm,0:km,2),'q2',        &
                                      q(0:im,0:jm,0:km,3),'q3',        &
                                      q(0:im,0:jm,0:km,5),'q5',        &
                                      q(0:im,0:jm,0:km,6),'q6' )
    !
    if(lio)  write(*,'(A,I1,A)')'  ** Gaussian pulse initialised.'
    !
    ! call mpistop
    !
  end subroutine accini
  !+-------------------------------------------------------------------+
  !| The end of the subroutine accini.                                 |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to generate an initial field test of      |
  !| 1D Shu-Osher problem.                                             |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 23-Mar-2021: Created by J. Fang @ STFC Daresbury Laboratory       |
  !+-------------------------------------------------------------------+
  subroutine shuosherini
    !
    use commvar,  only: mach
    use commarray,only: x,vel,rho,prs,spc,tmp,q,acctest_ref
    use fludyna,  only: thermal
    !
    ! local data
    integer :: i,j,k
    real(8) :: xc
    !
    do k=0,km
    do j=0,jm
    do i=0,im
      !
      xc=x(i,j,k,1)
      !
      if(xc<-4.d0) then
        rho(i,j,k)  =3.857143d0
        vel(i,j,k,1)=2.629369d0
        prs(i,j,k)  =10.33333d0
        !
        !
      else
        rho(i,j,k)  =1.d0+0.2d0*sin(5.d0*xc) 
        vel(i,j,k,1)=0.d0
        prs(i,j,k)  =1.d0
      endif
      vel(i,j,k,2)=0.d0
      vel(i,j,k,3)=0.d0
      !
      tmp(i,j,k)=thermal(density=rho(i,j,k),pressure=prs(i,j,k))
      !
    enddo
    enddo
    enddo
    !
    ! call tecbin('testout/tecinit'//mpirankname//'.plt',                &
    !                                   x(0:im,0:jm,0:km,1),'x',         &
    !                                   q(0:im,0:jm,0:km,1),'q1',        &
    !                                   q(0:im,0:jm,0:km,2),'q2',        &
    !                                   q(0:im,0:jm,0:km,3),'q3',        &
    !                                   q(0:im,0:jm,0:km,5),'q5',        &
    !                                   q(0:im,0:jm,0:km,6),'q6' )
    !
    if(lio)  write(*,'(A,I1,A)')'  ** shu-osher profile initialised.'
    !
    ! call mpistop
    !
  end subroutine shuosherini
  !+-------------------------------------------------------------------+
  !| The end of the subroutine shuosherini.                            |
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
    use fludyna,  only: thermal
    !
    ! local data
    integer :: i,j,k,jspc
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
      if(ndims>=2) vel(i,j,k,2)=vinf+var1*(x(i,j,k,1)-xc)
      if(ndims==3) vel(i,j,k,3)=0.d0
      prs(i,j,k)  =pinf-0.5d0*roinf*cvor**2/rvor**2*exp(-radi2)
      !
      tmp(i,j,k)=thermal(density=rho(i,j,k),pressure=prs(i,j,k))
      !
      if(num_species>=1) then
        !
        spc(i,j,k,1)=exp(-0.5d0*radi2)
        do jspc=2,num_species
          spc(i,j,k,jspc)=1.d0-spc(i,j,k,1)
        enddo
        !
      endif
      !
    enddo
    enddo
    enddo
    !
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
  !+-------------------------------------------------------------------+
  !| This subroutine is used to generate 3d TGV initial flow           |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 12-02-2021: Created by J. Fang @ STFC Daresbury Laboratory        |
  !+-------------------------------------------------------------------+
  subroutine tgvini
    !
    use commvar,  only: nondimen
    use commarray,only: x,vel,rho,prs,spc,tmp,q
    use fludyna,  only: thermal
#ifdef COMB
    use thermchem,only: tranco,convertxiyi
#endif
    !
    ! local data
    integer :: i,j,k,jspc
    real(8) :: l_0,miu
    !
    if(nondimen) then
      roinf=1.d0
      tinf=1.d0
      pinf=thermal(temperature=tinf,density=roinf)
      l_0=1.d0
      uinf=1.d0
    else
#ifdef COMB
      tinf=347.d0
      roinf=thermal(temperature=tinf,pressure=pinf,species=spcinf)
      l_0=xmax/(2.d0*pi)
      uinf=40.d0
      !
      ! call tranco(tmp=tinf,spc=spcinf,mu=miu,den=roinf)
      ! if(lio) print*,' ** miu=',miu,'Re=',roinf*uinf*l_0/miu,'pinf=',pinf
#endif
    endif
    !
    do k=0,km
    do j=0,jm
    do i=0,im
      rho(i,j,k)  =roinf
      vel(i,j,k,1)= uinf*sin(x(i,j,k,1)/l_0)*cos(x(i,j,k,2)/l_0)*cos(x(i,j,k,3)/l_0)
      vel(i,j,k,2)=-uinf*cos(x(i,j,k,1)/l_0)*sin(x(i,j,k,2)/l_0)*cos(x(i,j,k,3)/l_0)
      vel(i,j,k,3)=0.d0
      prs(i,j,k)  =pinf+roinf/16.d0*(uinf**2) &
                        *(cos(2.d0*x(i,j,k,1)/l_0)+cos(2.d0*x(i,j,k,2)/l_0)) &
                        *(cos(2.d0*x(i,j,k,3)/l_0)+2.d0)
      !
      if(nondimen) then
        tmp(i,j,k)=thermal(density=rho(i,j,k),pressure=prs(i,j,k))
        if(num_species>1) spc(i,j,k,1)=1.d0
      else 
        spc(i,j,k,:)=spcinf(:)
        tmp(i,j,k)=thermal(density=rho(i,j,k),pressure=prs(i,j,k),species=spc(i,j,k,:))
      endif 
      !
      if(num_species>=1 .and. nondimen) then
        !
        spc(i,j,k,1)=sin(x(i,j,k,1))**2
        do jspc=2,num_species
          spc(i,j,k,jspc)=1.d0-spc(i,j,k,1)
        enddo
        !
      endif
      !
    enddo
    enddo
    enddo
    !
    !
    ! call tecbin('testout/tecinit'//mpirankname//'.plt',                &
    !                                   x(0:im,0:jm,0:km,1),'x',         &
    !                                   x(0:im,0:jm,0:km,2),'y',         &
    !                                   x(0:im,0:jm,0:km,3),'z',         &
    !                                rho(0:im,0:jm,0:km)  ,'ro',         &
    !                                vel(0:im,0:jm,0:km,1),'u',          &
    !                                vel(0:im,0:jm,0:km,2),'v',          &
    !                                prs(0:im,0:jm,0:km)  ,'p',          &
    !                                tmp(0:im,0:jm,0:km)  ,'t' )
    !
    if(lio)  write(*,'(A,I1,A)')'  ** ',ndims,'-D TGV initialised.'
    !
  end subroutine tgvini
  !+-------------------------------------------------------------------+
  !| The end of the subroutine tgvini.                                 |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to generate an initial field for the      |
  !| simulation of Rayleigh–Taylor instability.                        |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 27-May-2020: Created by J. Fang @ STFC Daresbury Laboratory       |
  !+-------------------------------------------------------------------+
  subroutine rtini
    !
    use commarray,only: x,vel,rho,prs,tmp
    use fludyna,  only: thermal,sos
    !
    ! local data
    integer :: i,j,k
    !
    do k=0,km
    do j=0,jm
    do i=0,im
      !
      if(x(i,j,k,2)<0.5d0) then
        rho(i,j,k)=2.d0
        prs(i,j,k)=2.d0*x(i,j,k,2)+1.d0
      else
        rho(i,j,k)=1.d0
        prs(i,j,k)=x(i,j,k,2)+1.5d0
      endif
      !
      tmp(i,j,k)=thermal(density=rho(i,j,k),pressure=prs(i,j,k))
      !
      vel(i,j,k,1)=  0.d0
      vel(i,j,k,2)= -0.025*sos(tmp(i,j,k))*cos(8.d0*pi*x(i,j,k,1))
      vel(i,j,k,3)=  0.d0
      !
    enddo
    enddo
    enddo
    !
    if(lio)  write(*,'(A,I1,A)')'  ** ',ndims,'-D R–T instability initialised.'
    !
  end subroutine rtini
  !+-------------------------------------------------------------------+
  !| The end of the subroutine rtini.                                  |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to generate an initial field for the      |
  !| simulation of channel flow.                                       |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 27-May-2020: Created by J. Fang @ STFC Daresbury Laboratory       |
  !+-------------------------------------------------------------------+
  subroutine chanini
    !
    use commvar,  only: prandtl,mach,gamma,turbmode,                   &
                        xmin,xmax,ymin,ymax,zmin,zmax,Reynolds
    use commarray,only: x,vel,rho,prs,spc,tmp,q,dgrid,tke,omg,miut,res12
    use fludyna,  only: thermal,miucal
    use commfunc, only: dis2point2
    use readwrite,only: readprofile
    !
    ! local data
    integer :: i,j,k,l,ii,jj
    real(8) :: theta,theter,fx,gz,zl,randomv(15),ran
    real(8) :: delta,beta1,miu
    real(8) :: yh(0:jm),nth(0:jm),r12(0:jm)
    integer :: seed1,seed2,seed3,seed11,seed22,seed33
    integer, parameter :: nsemini = 1000
    real(8), dimension(3,nsemini) :: eddy, posvor
    real(8), dimension(3) :: dim_min, dim_max
    real(8) :: volsemini,rrand,ddx,ddy,ddz,lsem,upr,vpr,wpr,rrand1,   &
               init_noise,um,ftent
    !
    theta=0.1d0
    !
    beta1=0.075d0
    !
    if(trim(turbmode)=='udf1') then
      !
      call readprofile('Results/miut.dat',dir='j',                    &
                                              var1=yh,var2=nth,var3=r12)
      ! do j=0,jm
      !   print*,mpirank,'|',yh(j),nth(j),r12(j)
      ! enddo
      !
      allocate(miut(0:im,0:jm,0:km),res12(0:im,0:jm,0:km))
      !
    endif
    !
    if(ndims==2) then
      !
      do k=0,km
      do j=0,jm
      do i=0,im
        rho(i,j,k)  =roinf
        vel(i,j,k,1)=1.5d0*(1.d0-(x(i,j,k,2)-1.d0)**2)
        vel(i,j,k,2)=0.d0
        vel(i,j,k,3)=0.d0
        tmp(i,j,k)=1.d0+(gamma-1.d0)*prandtl*mach**2/3.d0*             &
                                     1.5d0*(1.d0-((x(i,j,k,2)-1.d0)**4))
        !
        prs(i,j,k)=thermal(density=rho(i,j,k),temperature=tmp(i,j,k))
        !
        if(num_species>0) then
          spc(i,j,k,1)=(tanh((x(i,j,k,2)-1.d0)/theta)+1.d0)*0.5d0
        endif
        !
        if(trim(turbmode)=='k-omega') then
          tke(i,j,k)=1.5d0*0.0001d0*10.d0
          !
          ! omg(i,j,k)=sqrt(tke(i,j,k))/(0.09d0)**0.25d0
          delta=dis2point2(x(i,j,k,:),x(i,j+1,k,:))
          miu=miucal(tmp(i,j,k))/Reynolds
          omg(i,j,k)=60.d0*miu/rho(i,j,k)/beta1/delta
          !
        elseif(trim(turbmode)=='udf1') then
          miut(i,j,k)=nth(j)
          res12(i,j,k)=r12(j)
        endif
        !
      enddo
      enddo
      enddo
      !
    elseif(ndims==3) then
      !
      ! if(mpirank==0) then
      !   do l=1,15
      !     call random_number(randomv(l))
      !   end do
      ! endif
      ! !
      ! call bcast(randomv)
      ! !
      ! do k=0,km
      ! do j=0,jm
      ! do i=0,im
      !   !
      !   call random_number(ran)
      !   ran=ran*2.d0-1.d0
      !   !
      !   theter=x(i,0,k,1)/(xmax-xmin)*2.d0*pi
      !   fx=4.d0*dsin(theter)*(1.d0-dcos(theter))*0.192450089729875d0
      !   !
      !   gz=0.d0
      !   do l=1,10
      !     !
      !     if(l==1) then
      !       zl=0.2d0/(1.d0-0.8d0**10)
      !     else
      !       zl=zl*0.8d0
      !     end if
      !     !
      !     gz=gz+zl*dsin(2.d0*pi*l*(x(i,0,k,3)/(zmax-zmin)+randomv(l)))
      !     !
      !   end do
      !   !
      !   rho(i,j,k)  =roinf
      !   vel(i,j,k,1)=1.5d0*(1.d0-(x(i,j,k,2)-1.d0)**2)*(1.d0+0.3d0*fx*gz+0.1d0*ran)
      !   vel(i,j,k,2)=0.d0
      !   vel(i,j,k,3)=0.d0
      !   tmp(i,j,k)=1.d0+(gamma-1.d0)*prandtl*mach**2/3.d0*             &
      !                                1.5d0*(1.d0-((x(i,j,k,2)-1.d0)**4))
      !   !
      !   prs(i,j,k)=thermal(density=rho(i,j,k),temperature=tmp(i,j,k))
      !   !
      !   if(num_species>0) then
      !     spc(i,j,k,1)=(tanh((x(i,j,k,2)-1.d0)/theta)+1.d0)*0.5d0
      !   endif
      !   !
      ! enddo
      ! enddo
      ! enddo
      !
      ! copied from xcompact 
      !! Simplified version of SEM 
      init_noise=0.03d0
      !
      dim_min(1) = 0.d0
      dim_min(2) = 0.d0
      dim_min(3) = 0.d0
      dim_max(1) = xmax
      dim_max(2) = ymax
      dim_max(3) = zmax
      volsemini = xmax * ymax * zmax
      !
      ! 3 int to get different random numbers
      seed1 =  2345
      seed2 = 13456
      seed3 = 24567
      do jj=1,nsemini
        !
        ! Vortex Position
        do ii=1,3
          seed11 = return_30k(seed1+jj*2+ii*379)
          seed22 = return_30k(seed2+jj*5+ii*5250)
          seed33 = return_30k(seed3+jj*3+ii*8170)
          rrand1  = real(r8_random(seed11, seed22, seed33),8)
          call random_number(rrand)
          !write(*,*) ' rr r1 ', rrand, rrand1
          posvor(ii,jj) = dim_min(ii)+(dim_max(ii)-dim_min(ii))*rrand
        enddo
        !
        ! Eddy intensity
        do ii=1,3
           seed11 = return_30k(seed1+jj*7+ii*7924)
           seed22 = return_30k(seed2+jj*11+ii*999)
           seed33 = return_30k(seed3+jj*5+ii*5054)
           rrand1  = real(r8_random(seed11, seed22, seed33),8)
           call random_number(rrand)
           !write(*,*) ' rr r1 ', rrand, rrand1
           if (rrand <= 0.5d0) then
              eddy(ii,jj) = -1.d0
           else
              eddy(ii,jj) =  1.d0
           endif 
        enddo
        !
      enddo
      !
      ! Loops to apply the fluctuations 
      do k=0,km
      do j=0,jm
      do i=0,im
        !
        rho(i,j,k)  =roinf
        vel(i,j,k,1)=1.5d0*(1.d0-(x(i,j,k,2)-1.d0)**2)
        vel(i,j,k,2)=0.d0
        vel(i,j,k,3)=0.d0
        tmp(i,j,k)=1.d0+(gamma-1.d0)*prandtl*mach**2/3.d0*             &
                                     1.5d0*(1.d0-((x(i,j,k,2)-1.d0)**4))
        !
        lsem = 0.15d0 ! For the moment we keep it constant
        upr = 0.d0
        vpr = 0.d0
        wpr = 0.d0
        do jj=1,nsemini
          !
          ddx = abs(x(i,j,k,1)-posvor(1,jj))
          ddy = abs(x(i,j,k,2)-posvor(2,jj))
          ddz = abs(x(i,j,k,3)-posvor(3,jj))
          if (ddx < lsem .and. ddy < lsem .and. ddz < lsem) then
            ! coefficients for the intensity of the fluctuation
            ftent = (1.d0-ddx/lsem)*(1.d0-ddy/lsem)*(1.d0-ddz/lsem)
            ftent = ftent / (sqrt(num2d3*lsem))**3
            upr = upr + eddy(1,jj) * ftent
            vpr = vpr + eddy(2,jj) * ftent
            wpr = wpr + eddy(3,jj) * ftent
          endif
          !
        enddo
        !
        upr = upr * sqrt(volsemini/nsemini)
        vpr = vpr * sqrt(volsemini/nsemini)
        wpr = wpr * sqrt(volsemini/nsemini)
        ! 
        um=vel(i,j,k,1)
        !
        vel(i,j,k,1)=upr*sqrt(num2d3*init_noise*um) + um
        vel(i,j,k,2)=vpr*sqrt(num2d3*init_noise*um)
        vel(i,j,k,3)=wpr*sqrt(num2d3*init_noise*um)
        !
        prs(i,j,k)=thermal(density=rho(i,j,k),temperature=tmp(i,j,k))
        !
      enddo
      enddo
      enddo
      !
    endif
    !
    !
    ! call tecbin('testout/tecinit'//mpirankname//'.plt',                &
    !                                   x(0:im,0:jm,0:km,1),'x',         &
    !                                   x(0:im,0:jm,0:km,2),'y',         &
    !                                   x(0:im,0:jm,0:km,3),'z',         &
    !                                rho(0:im,0:jm,0:km)  ,'ro',         &
    !                                vel(0:im,0:jm,0:km,1),'u',          &
    !                                vel(0:im,0:jm,0:km,2),'v',          &
    !                                prs(0:im,0:jm,0:km)  ,'p' )
    
    if(lio)  write(*,'(A,I1,A)')'  ** ',ndims,'-D channel flow initialised.'
    !
  end subroutine chanini
  !+-------------------------------------------------------------------+
  !| The end of the subroutine chanini.                                |
  !+-------------------------------------------------------------------+
  !!
  !+-------------------------------------------------------------------+
  !| This subroutine is used to generate a initial jet flow.           |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 24-02-2021: Created by J. Fang @ STFC Daresbury Laboratory        |
  !+-------------------------------------------------------------------+
  subroutine jetini
    !
    use commarray,only: x,vel,rho,prs,spc,tmp,q
    use fludyna,  only: thermal,jetvel
    !
    ! local data
    integer :: i,j,k
    real(8) :: radi
    !
    do k=0,km
    do j=0,jm
    do i=0,im
      rho(i,j,k)  = roinf
      !
      radi=sqrt(x(i,j,k,2)**2+x(i,j,k,3)**2)
      ! radi=abs(x(i,j,k,2))
      vel(i,j,k,:)= jetvel(radi)
      ! vel(i,j,k,1)= 0.d0
      vel(i,j,k,2)= 0.d0
      vel(i,j,k,3)= 0.d0
      !
      tmp(i,j,k)  = tinf
      !
      prs(i,j,k)=thermal(density=rho(i,j,k),temperature=tmp(i,j,k))
      !
      if(num_species>1) then
        spc(i,j,k,1)=0.d0
        !
        spc(i,j,k,num_species)=1.d0-sum(spc(i,j,k,1:num_species-1))
        !
      endif
      !
    enddo
    enddo
    enddo
    !
    !
    ! call tecbin('testout/tecinit'//mpirankname//'.plt',                &
    !                                   x(0:im,0:jm,0:km,1),'x',         &
    !                                   x(0:im,0:jm,0:km,2),'y',         &
    !                                   x(0:im,0:jm,0:km,3),'z',         &
    !                                rho(0:im,0:jm,0:km)  ,'ro',         &
    !                                vel(0:im,0:jm,0:km,1),'u',          &
    !                                vel(0:im,0:jm,0:km,2),'v',          &
                                   ! prs(0:im,0:jm,0:km)  ,'p',          &
    !                                tmp(0:im,0:jm,0:km)  ,'t' )
    !
    if(lio)  write(*,'(A,I1,A)')'  ** ',ndims,'-D jet flow initialised.'
    !
  end subroutine jetini
  !+-------------------------------------------------------------------+
  !| The end of the subroutine jetini.                                 |
  !+-------------------------------------------------------------------+
  !!
  !+-------------------------------------------------------------------+
  !| This subroutine is used to generate a initial mixing layer flow.  |
  !+-------------------------------------------------------------------+
  !| ref: Li, Z., Jaberi, F. 2010. Numerical Investigations of         |
  !|      Shock-Turbulence Interaction in a Planar Mixing Layer.       | 
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 24-02-2021: Created by J. Fang @ STFC Daresbury Laboratory        |
  !+-------------------------------------------------------------------+
  subroutine mixlayerini
    !
    use commarray,only: x,vel,rho,prs,spc,tmp,q
    use fludyna,  only: thermal,mixinglayervel
    !
    ! local data
    integer :: i,j,k
    real(8) :: radi
    !
    do k=0,km
    do j=0,jm
    do i=0,im
      rho(i,j,k)  = roinf
      !
      vel(i,j,k,:)= mixinglayervel(x(i,j,k,2))
      !
      tmp(i,j,k)  = tinf
      !
      prs(i,j,k)=thermal(density=rho(i,j,k),temperature=tmp(i,j,k))
      !
      if(num_species>1) then
        spc(i,j,k,1)=0.d0
        !
        spc(i,j,k,num_species)=1.d0-sum(spc(i,j,k,1:num_species-1))
        !
      endif
      !
    enddo
    enddo
    enddo
    !
    !
    ! call tecbin('testout/tecinit'//mpirankname//'.plt',                &
    !                                   x(0:im,0:jm,0:km,1),'x',         &
    !                                   x(0:im,0:jm,0:km,2),'y',         &
    !                                   x(0:im,0:jm,0:km,3),'z',         &
    !                                rho(0:im,0:jm,0:km)  ,'ro',         &
    !                                vel(0:im,0:jm,0:km,1),'u',          &
    !                                vel(0:im,0:jm,0:km,2),'v',          &
                                   ! prs(0:im,0:jm,0:km)  ,'p',          &
    !                                tmp(0:im,0:jm,0:km)  ,'t' )
    !
    if(lio)  write(*,'(A,I1,A)')'  ** ',ndims,'-D jet flow initialised.'
    !
  end subroutine mixlayerini
  !+-------------------------------------------------------------------+
  !| The end of the subroutine mixlayerini.                            |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to correct the boundary layer initial field.   |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 24-10-2021: Created by J. Fang @ STFC Daresbury Laboratory        |
  !+-------------------------------------------------------------------+
  subroutine blcorrect
    !
    use commarray,only: x,vel,rho,prs,spc,tmp,q
    use fludyna,  only: thermal
    use bc,       only: rho_prof,vel_prof,tmp_prof,prs_prof,spc_prof
    !
    ! local data
    integer :: i,j,k
    !
    do k=0,km
    do j=0,jm
    do i=0,im
      !
      if(x(i,j,k,2)>=5.d0*nomi_thick) then
        !
        rho(i,j,k)  = rho_prof(j)
        !
        vel(i,j,k,1)= vel_prof(j,1)
        vel(i,j,k,2)= vel_prof(j,2)
        !
        tmp(i,j,k)  = tmp_prof(j)
        !
        prs(i,j,k)=thermal(density=rho(i,j,k),temperature=tmp(i,j,k))
        ! 
      endif
      !
    enddo
    enddo
    enddo
    !
  end subroutine blcorrect
  !+-------------------------------------------------------------------+
  !| The end of the subroutine blcorrect.                              |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used initilise a boundary layer flow.          |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 27-09-2021: Created by J. Fang @ STFC Daresbury Laboratory        |
  !+-------------------------------------------------------------------+
  subroutine blini
    !
    use commvar,  only: turbmode,Reynolds
    use commarray,only: x,vel,rho,prs,spc,tmp,q,tke,omg
    use fludyna,  only: thermal,miucal
    use bc,       only: rho_prof,vel_prof,tmp_prof,prs_prof,spc_prof
    use commfunc, only: dis2point2
    !
    ! local data
    integer :: i,j,k
    real(8) :: radi,miu,delta,beta1
    !
    beta1=0.075d0
    !
    do k=0,km
    do j=0,jm
    do i=0,im
      !
      rho(i,j,k)  = rho_prof(j)
      !
      vel(i,j,k,1)= vel_prof(j,1)
      vel(i,j,k,2)= vel_prof(j,2)
      !
      tmp(i,j,k)  = tmp_prof(j)
      !
      prs(i,j,k)=thermal(density=rho(i,j,k),temperature=tmp(i,j,k))
      !
      if(num_species>1) then
        spc(i,j,k,1)=0.d0
        !
        spc(i,j,k,num_species)=1.d0-sum(spc(i,j,k,1:num_species-1))
        !
      endif
      !
      if(trim(turbmode)=='k-omega') then
        tke(i,j,k)=1.5d0*0.0001d0
        !
        ! omg(i,j,k)=sqrt(tke(i,j,k))/(0.09d0)**0.25d0
        delta=dis2point2(x(i,j,k,:),x(i,j+1,k,:))
        miu=miucal(tmp(i,j,k))/Reynolds
        omg(i,j,k)=60.d0*miu/rho(i,j,k)/beta1/delta
        !
      endif
      !
    enddo
    enddo
    enddo
    !
    !
    ! call tecbin('testout/tecinit'//mpirankname//'.plt',                &
    !                                   x(0:im,0:jm,0:km,1),'x',         &
    !                                   x(0:im,0:jm,0:km,2),'y',         &
    !                                   x(0:im,0:jm,0:km,3),'z',         &
    !                                rho(0:im,0:jm,0:km)  ,'ro',         &
    !                                vel(0:im,0:jm,0:km,1),'u',          &
    !                                vel(0:im,0:jm,0:km,2),'v',          &
    !                                prs(0:im,0:jm,0:km)  ,'p',          &
    !                                tmp(0:im,0:jm,0:km)  ,'t' )
    !
    if(lio)  write(*,'(A,I1,A)')'  ** ',ndims,'-D BL flow initialised.'
    !
  end subroutine blini
  !+-------------------------------------------------------------------+
  !| The end of the subroutine blini.                                  |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is read the boundary layer profile.               |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 27-09-2021: Created by J. Fang @ STFC Daresbury Laboratory        |
  !+-------------------------------------------------------------------+
  subroutine blprofile
    !
    use readwrite,only: readprofile
    use bc,       only: rho_prof,vel_prof,tmp_prof,prs_prof,spc_prof
    use fludyna,  only: thermal
    use stlaio,  only: get_unit
    !
    ! local data
    integer :: fh
    !
    allocate( rho_prof(0:jm),tmp_prof(0:jm),prs_prof(0:jm),          &
              vel_prof(0:jm,1:3),spc_prof(0:jm,1:num_species) )
    !
    if(lio) then
      fh=get_unit()
      !
      open(fh,file='datin/inlet.prof',action='read',form='formatted')
      read(fh,*)
      read(fh,*)
      read(fh,*)nomi_thick,disp_thick,mome_thick,fric_velocity
      close(fh)
      write(*,"(A)")'  -----------------------------------------------'
      write(*,"(A)")'            δ          δ*           θ        utau'
      write(*,"(1X,4(F12.7))")nomi_thick,disp_thick,mome_thick,fric_velocity
      write(*,"(A)")'  -----------------------------------------------'
      !
    endif
    !
    call bcast(nomi_thick)
    call bcast(disp_thick)
    call bcast(mome_thick)
    call bcast(fric_velocity)
    !
    call readprofile('datin/inlet.prof',dir='j',                     &
                              var1=rho_prof,     var2=vel_prof(:,1), &
                              var3=vel_prof(:,2),var4=tmp_prof,skipline=4)
    !
    prs_prof=thermal(density=rho_prof,temperature=tmp_prof,dim=jm+1)
    !
  end subroutine blprofile
  !+-------------------------------------------------------------------+
  !| The end of the subroutine blprofile.                              |
  !+-------------------------------------------------------------------+
  !  
  !+-------------------------------------------------------------------+
  !| This subroutine is used to generate a initial mixing layer flow.  |
  !+-------------------------------------------------------------------+
  !| ref: Li, Z., Jaberi, F. 2010. Numerical Investigations of         |
  !|      Shock-Turbulence Interaction in a Planar Mixing Layer.       | 
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 24-02-2021: Created by J. Fang @ STFC Daresbury Laboratory        |
  !+-------------------------------------------------------------------+
  subroutine wtini
    !
    use commarray,only: x,vel,rho,prs,spc,tmp,q
    use fludyna,  only: thermal,mixinglayervel
    !
    ! local data
    integer :: i,j,k
    real(8) :: radi
    !
    do k=0,km
    do j=0,jm
    do i=0,im
      rho(i,j,k)  = roinf
      !
      vel(i,j,k,1)= uinf
      vel(i,j,k,2)= vinf
      vel(i,j,k,3)= winf
      !
      tmp(i,j,k)  = tinf
      !
      prs(i,j,k)=thermal(density=rho(i,j,k),temperature=tmp(i,j,k))
      !
      if(num_species>1) then
        spc(i,j,k,1)=0.d0
        !
        spc(i,j,k,num_species)=1.d0-sum(spc(i,j,k,1:num_species-1))
        !
      endif
      !
    enddo
    enddo
    enddo
    !
    !
    ! call tecbin('testout/tecinit'//mpirankname//'.plt',                &
    !                                   x(0:im,0:jm,0:km,1),'x',         &
    !                                   x(0:im,0:jm,0:km,2),'y',         &
    !                                   x(0:im,0:jm,0:km,3),'z',         &
    !                                rho(0:im,0:jm,0:km)  ,'ro',         &
    !                                vel(0:im,0:jm,0:km,1),'u',          &
    !                                vel(0:im,0:jm,0:km,2),'v',          &
                                   ! prs(0:im,0:jm,0:km)  ,'p',          &
    !                                tmp(0:im,0:jm,0:km)  ,'t' )
    !
    if(lio)  write(*,'(A,I1,A)')'  ** ',ndims,'-D wind tunnel initialised.'
    !
  end subroutine wtini
  !+-------------------------------------------------------------------+
  !| The end of the subroutine wtini.                                 |
  !+-------------------------------------------------------------------+
  !  
  !+-------------------------------------------------------------------+
  !| This subroutine is used to generate initial field for flow past a |
  !| cylinder flow.                                                    |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 11-03-2021: Created by J. Fang @ STFC Daresbury Laboratory        |
  !+-------------------------------------------------------------------+
  subroutine cylinderini
    !
    use commarray,only: x,vel,rho,prs,spc,tmp,q
    use fludyna,  only: thermal,jetvel
    !
    ! local data
    integer :: i,j,k
    !
    do k=0,km
    do j=0,jm
    do i=0,im
      rho(i,j,k)  = roinf
      vel(i,j,k,1)= uinf
      vel(i,j,k,2)= 0.d0
      vel(i,j,k,3)= 0.d0
      !
      tmp(i,j,k)  = tinf
      !
      prs(i,j,k)=thermal(density=rho(i,j,k),temperature=tmp(i,j,k))
      !
      if(num_species>=1) then
        spc(i,j,k,1)=0.d0
      endif
      !
    enddo
    enddo
    enddo
    !
    !
    call tecbin('testout/tecinit'//mpirankname//'.plt',                &
                                      x(0:im,0:jm,0:km,1),'x',         &
                                      x(0:im,0:jm,0:km,2),'y',         &
                                      x(0:im,0:jm,0:km,3),'z',         &
                                   rho(0:im,0:jm,0:km)  ,'ro',         &
                                   vel(0:im,0:jm,0:km,1),'u',          &
                                   vel(0:im,0:jm,0:km,2),'v',          &
                                   prs(0:im,0:jm,0:km)  ,'p',          &
                                   tmp(0:im,0:jm,0:km)  ,'t' )
    !
    if(lio)  write(*,'(A,I1,A)')'  ** ',ndims,'-D jet flow initialised.'
    !
  end subroutine cylinderini
  !+-------------------------------------------------------------------+
  !| The end of the subroutine cylinderini.                            |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to generate an initial field for the      |
  !| simulation of 0D reactor.                                                    |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 24-01-2022: Created by Created by Z.X. Chen @ Peking University   |
  !+-------------------------------------------------------------------+
  subroutine reactorini
    !
    use commvar,  only: pinf
    use commarray,only: x,vel,rho,prs,spc,tmp,q
    use fludyna,  only: thermal
    use thermchem,only: convertxiyi,spcindex
    !
#ifdef COMB
    ! local data
    integer :: i,j,k
    real(8) :: tmpr,specr(num_species),specx(num_species)
    !
    ! tmpr=1000.d0/0.8d0
    tmpr=1000.d0
    ! prin=13.5d5
    !reactants
    specx(:)=0.d0
    specx(spcindex('H2'))=0.2d0
    specx(spcindex('O2'))=0.1d0
    !
    ! specx(spcindex('H2'))=0.2867d0
    ! specx(spcindex('O2'))=0.1434d0
    ! specx(spcindex('H2O'))=0.1819d0
    !
    specx(spcindex('N2'))=1.d0-sum(specx)
    call convertxiyi(specx(:),specr(:),'X2Y')
    !
    !!
    ! specr(:)=0.d0
    ! specr(spcindex('nc7h16'))=0.07247482382311918d0
    ! specr(spcindex('o2'))=0.28285951066551174d0
    ! specr(spcindex('he'))=0.08059381448746926d0
    ! specr(spcindex('n2'))=1.d0-sum(specr)
    !
    ! specr(spcindex('CH4'))=0.055d0
    ! specr(spcindex('O2'))=0.220185d0
    ! specr(spcindex('N2'))=1.d0-sum(specr)
    !
    ! print*,specr
    ! stop
    ! specr(1)=0.055d0
    ! specr(2)=0.220185d0
    ! specr(num_species)=1.d0-sum(specr)
    ! specr(1)=0.06218387d0
    ! specr(2)=0.21843332
    ! specr(num_species)=1.d0-sum(specr)
    !
    do k=0,km
    do j=0,jm
    do i=0,im
      !
      vel(i,j,k,:)= 0.d0
      !
      tmp(i,j,k)  = tmpr
      !
      prs(i,j,k)=pinf
      !
      spc(i,j,k,:)=specr(:)
      !
      rho(i,j,k)=thermal(pressure=prs(i,j,k),temperature=tmp(i,j,k), &
                          species=spc(i,j,k,:))
      !
      ! print*,tmp(i,j,k),prs(i,j,k),rho(i,j,k),spc(i,j,k,:)
    enddo
    enddo
    enddo
    !
    if(lio)  write(*,'(A,I1,A)')'  ** reactor initialised.'
    !
#endif
  !
  end subroutine reactorini
  !+-------------------------------------------------------------------+
  !| The end of the subroutine reactorini.                             |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to generate an initial field for the      |
  !| simulation of 1D flame.                                           |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 30-Jan-2022: Created by Z.X. Chen @ Peking University             |
  !+-------------------------------------------------------------------+
  subroutine onedflameini
    !
    use commvar,  only: pinf,ia,num_species
    use commarray,only: x,vel,rho,prs,spc,tmp,q
    use fludyna,  only: thermal
    !
#ifdef COMB
    use thermchem,only: convertxiyi,spcindex
    !
    ! local data
    integer :: i,j,k
    real(8) ::  xc,yc,zc,tmpr,tmpp,xloc,xwid,specr(num_species),  &
      specp(num_species),arg,prgvar,masflx,specx(num_species),yloc
    !
    tmpr=650.d0
    xloc=0.25d-2
    yloc=xloc
    xwid=0.5d-3
    !
    !reactants
    specr(:)=0.d0
    specr(spcindex('H2'))=0.031274  
    specr(spcindex('O2'))=0.22563  
    specr(spcindex('N2'))=1.d0-sum(specr)
    !
    ! specx(:)=0.d0
    ! specx(spcindex('H2'))=0.2957746478873239d0
    ! specx(spcindex('O2'))=0.14788732394366194d0
    ! specx(spcindex('N2'))=1.d0-sum(specx)
    !
    ! specr(spcindex('CH4'))=5.5045872d-2
    ! specr(spcindex('O2'))=2.20183486d-1
    ! call convertxiyi(specx(:),specr(:),'X2Y')
    !
    ! specr(spcindex('nc7h16'))=0.06218387d0
    ! specr(spcindex('o2'))=0.21843332
    ! specr(spcindex('n2'))=1.d0-sum(specr)
    !
    !products
    tmpp=1700.d0
    specp(:)=0.d0
    ! specp(spcindex('CO2'))=1.51376d-1
    ! specp(spcindex('H2O'))=1.23853d-1
    ! specp(spcindex('N2'))=1.d0-sum(specp)
    !
    do k=0,km
    do j=0,jm
    do i=0,im
      !
      xc=x(i,j,k,1)
      yc=x(i,j,k,2)
      !
      if(ndims==3) stop '!!Error - 3D case not configured!!'
      if(abs(xc-xloc)<xwid*0.5d0*1.2d0) then 
        prgvar=1.d0
        if(abs(xc-xloc)>xwid*0.5d0) &
        prgvar=1.d0-(abs(xc-xloc)-(xwid*0.5d0))/(xwid*0.5d0*0.2d0)
      else
        prgvar=0.d0
      endif 
      !
      spc(i,j,k,:)=specr(:)!+prgvar*(specp(js)-specr(js))
      !
      vel(i,j,k,:)=0.d0
      !
      tmp(i,j,k)=tmpr+prgvar*(tmpp-tmpr)
      !
      prs(i,j,k)=pinf
      !
      rho(i,j,k)=thermal(pressure=prs(i,j,k),temperature=tmp(i,j,k), &
                          species=spc(i,j,k,:))
      !
      ! print*,tmp(i,j,k),prs(i,j,k),rho(i,j,k),spc(i,j,k,:)
    enddo
    enddo
    enddo
    !
    if(lio)  write(*,'(A,I1,A)')'  ** 1D flame initialised.'
    !
#endif
    !
  end subroutine onedflameini
  !+-------------------------------------------------------------------+
  !| The end of the subroutine onedflameini.                           |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to generate an initial field for the      |
  !| simulation of H2 supersonic jet flame                             |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 04-Feb-2022: Created by Z.X. Chen @ Peking University             |
  !+-------------------------------------------------------------------+
  subroutine h2supersonicini
#ifdef COMB
    !
    use commvar,  only: &
      pinf,uinf,tinf,num_species,dj_i,dj_o,dco_i,flowtype,ymax
    use commarray,only: x,vel,rho,prs,spc,tmp,q
    use fludyna,  only: thermal,multistream_inflow
    !
    use thermchem,only: convertxiyi,spcindex
    !
    ! local data
    integer :: i,j,k
    real(8) :: rb,ctr,xc,yb,zb,arg,val,xloc,xwid
    real(8) :: vel0(ndims),prs0,rho0,tmp0,spc0(num_species)
    !
    ctr=0.5d0*ymax
    !
    do k=0,km
    do j=0,jm
    do i=0,im
      !
      xc=x(i,j,k,1)
      yb=x(i,j,k,2)
      zb=x(i,j,k,3)
      !
      if(ndims<3) then
        rb=abs(yb-ctr)
      else
        rb=sqrt((yb-ctr)**2+(zb-ctr)**2)
      endif
      !
      if(.true. .and. rb<=0.5d0*dj_i) then
        call multistream_inflow(stream='fuel',rho=rho0,vel=vel0, &
          prs=prs0,tmp=tmp0,spc=spc0,rovd=rb/dj_i)
          !
      elseif(.true. .and. rb>0.5d0*dj_o .and. rb<=0.5d0*dco_i) then
        call multistream_inflow(stream='hotcoflow',rho=rho0,vel=vel0, &
          prs=prs0,tmp=tmp0,spc=spc0, &
          rovd=abs(rb-0.25d0*(dco_i+dj_o))/(0.5d0*(dco_i-dj_o)))
        !
      elseif(.true. .and. rb>0.5d0*dco_i) then
        call multistream_inflow(stream='air',rho=rho0,vel=vel0, &
          prs=prs0,tmp=tmp0,spc=spc0)
        !
      else
        !
        prs0=pinf
        tmp0=tinf
        spc0(:)=spcinf(:)
        rho0=roinf
        vel0(1)=uinf*100.d0
        vel0(2:ndims)=0.d0
        !
      endif 
      !
      ! erf profile in the axial diretion
      xloc=2.d0*dj_i
      xwid=2.d0*dj_i
      arg=-1.d0*(xc-xloc)/xwid
      val=0.5d0*(1.0d0+erf(arg))
      !
      !0.92135039647485750
      ! vel(i,j,k,1)=uinf+val*(vel0(1)-uinf)
      ! vel(i,j,k,2:ndims)=0.d0
      ! prs(i,j,k)=pinf
      ! tmp(i,j,k)=tinf+val*(tmp0-tinf)
      ! spc(i,j,k,:)=spcinf(:)+val*(spc0(:)-spcinf(:))
      vel(i,j,k,1)=uinf+0.92135039647485750*(vel0(1)-uinf)
      vel(i,j,k,2:ndims)=0.d0
      prs(i,j,k)=pinf
      tmp(i,j,k)=tinf+0.92135039647485750*(tmp0-tinf)
      spc(i,j,k,:)=spcinf(:)+0.92135039647485750*(spc0(:)-spcinf(:))
      ! vel(i,j,k,1)=uinf
      ! vel(i,j,k,2:ndims)=0.d0
      ! prs(i,j,k)=pinf
      ! tmp(i,j,k)=tinf
      ! spc(i,j,k,:)=spcinf(:)
      !
      rho(i,j,k)=thermal(pressure=prs(i,j,k),temperature=tmp(i,j,k), &
                          species=spc(i,j,k,:))
      !
      ! print*,tmp(i,j,k),prs(i,j,k),rho(i,j,k),spc(i,j,k,:)
    enddo
    enddo
    enddo
    !
    if(lio)  write(*,'(A,I1,3(A))')  &
      '  ** ',ndims,'-D ',trim(flowtype),' initialised.'
    !
    call tecbin('testout/tecini'//mpirankname//'.plt',                &
                                      x(0:im,0:jm,0:km,1),'x',        &
                                      x(0:im,0:jm,0:km,2),'y',        &
                                      x(0:im,0:jm,0:km,3),'z',        &
                                      rho(0:im,0:jm,0:km),'ro',       &
                                    vel(0:im,0:jm,0:km,1),'u',        &
                                    vel(0:im,0:jm,0:km,2),'v',        &
                                      tmp(0:im,0:jm,0:km),'T',        &
                                    spc(0:im,0:jm,0:km,1),'Y1' )
    !
    !
#endif
    !
  end subroutine h2supersonicini
  !+-------------------------------------------------------------------+
  !| The end of the subroutine h2supersonicini.                        |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to generate an initial field for the      |
  !| simulation of TGV flame.                                          |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 30-Mar-2022: Created by Yifan Xu @ Peking University              |
  !+-------------------------------------------------------------------+
  subroutine tgvflameini
    !
    use commvar,  only: nondimen,xmax
    use commarray,only: x,vel,rho,prs,spc,tmp,q
    use fludyna,  only: thermal
    !
#ifdef COMB
    !
    use thermchem,only : tranco,spcindex,mixture
    use cantera 
    !
    ! local data
    integer :: i,j,k,jspc
    real(8) :: miu,xc,yc,zc,xloc,l_0,xwid,specr(num_species), &
      specp(num_species),prgvar
    !
    tinf=300.d0
    xloc=xmax/2.d0
    xwid=xmax/8.d0
    !
    l_0=xmax/(2.d0*pi)
    uinf=4.d0
    roinf=thermal(temperature=tinf,pressure=pinf,species=spcinf)
    !
    ! nonpremixed reactants include fuel and oxidizer
    specr(:)=0.d0
    specr(spcindex('H2'))=0.0556 
    specr(spcindex('O2'))=0.233  
    specr(spcindex('N2'))=1.d0-sum(specr)
    !
    call tranco(den=roinf,tmp=tinf,mu=miu,spc=spcinf)
    if(lio) print*, &
      ' ** nu=',miu/roinf,'Re=',roinf*uinf*xmax/miu,'pinf=',pinf
    !
    do k=0,km
    do j=0,jm
    do i=0,im
      !
      xc=x(i,j,k,1)
      yc=x(i,j,k,2)
      zc=x(i,j,k,3)
      !
      prgvar=0.5d0*(1.d0+tanh(3.d0*(abs(xc-xloc)-xwid)/xwid))
      specp(:)=0.d0
      specp(spcindex('H2'))=specr(spcindex('H2'))*(1.d0-prgvar)
      specp(spcindex('O2'))=specr(spcindex('O2'))*prgvar
      specp(spcindex('N2'))=1.d0-sum(specp)
      !
      spc(i,j,k,:)=specp(:)
      !
      spc(i,j,k,spcindex('N2'))=1.d0-(sum(spc(i,j,k,:))-spc(i,j,k,spcindex('N2')))
      !
      prs(i,j,k)=pinf
      tmp(i,j,k)=tinf
      ! get initial rho
      roinf=thermal(temperature=tmp(i,j,k),pressure=prs(i,j,k),species=spc(i,j,k,:))
      !
      ! |--CANTERA--|
      call setState_TPY(mixture,tmp(i,j,k),prs(i,j,k),spc(i,j,k,:))
      call equilibrate(mixture,'HP')
      tmp(i,j,k)=temperature(mixture)
      call getMassFractions(mixture,specp(:))
      spc(i,j,k,:)=specp(:)
      !
      rho(i,j,k)=thermal(temperature=tmp(i,j,k),pressure=prs(i,j,k), &
                          species=spc(i,j,k,:))
      ! set velocity and scale vx 
      vel(i,j,k,1)= uinf*sin(xc/l_0)*cos(yc/l_0)*cos(zc/l_0)*roinf/rho(i,j,k)
      vel(i,j,k,2)=-uinf*cos(xc/l_0)*sin(yc/l_0)*cos(zc/l_0)
      vel(i,j,k,3)= 0.d0
      !
      prs(i,j,k)  =pinf+roinf/16.d0*(uinf**2) &
                        *(cos(2.d0*x(i,j,k,1)/l_0)+cos(2.d0*x(i,j,k,2)/l_0)) &
                        *(cos(2.d0*x(i,j,k,3)/l_0)+2.d0)
      !
    enddo
    enddo
    enddo
    !
    if(lio)  write(*,'(A,I1,A)')'  ** ',ndims,'-D tgv flame initialised.'
    !
#endif
    !
  end subroutine tgvflameini
  !+-------------------------------------------------------------------+
  !| The end of the subroutine tgvflameini.                           |
  !+-------------------------------------------------------------------+
  !
  function return_30k(x) result(y)
  
    integer ( kind = 4 ), intent(in) :: x
    integer ( kind = 4 )             :: y
    integer ( kind = 4 ), parameter  :: xmax = 30000
  
    y = iabs(x) - int(iabs(x)/xmax)*xmax
  end function return_30k
  !+-------------------------------------------------------------------+
  !| The end of the function return_30k.                               |
  !+-------------------------------------------------------------------+
  !
  function r8_random ( s1, s2, s3 )
    !*****************************************************************************80
    !
    !! R8_RANDOM returns a pseudorandom number between 0 and 1.
    !
    !  Discussion:
    !
    !    This function returns a pseudo-random number rectangularly distributed
    !    between 0 and 1.   The cycle length is 6.95E+12.  (See page 123
    !    of Applied Statistics (1984) volume 33), not as claimed in the
    !    original article.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    08 July 2008
    !
    !  Author:
    !
    !    FORTRAN77 original version by Brian Wichman, David Hill.
    !    FORTRAN90 version by John Burkardt.
    !
    !  Reference:
    !
    !    Brian Wichman, David Hill,
    !    Algorithm AS 183: An Efficient and Portable Pseudo-Random
    !    Number Generator,
    !    Applied Statistics,
    !    Volume 31, Number 2, 1982, pages 188-190.
    !
    !  Parameters:
    !
    !    Input/output, integer ( kind = 4 ) S1, S2, S3, three values used as the
    !    seed for the sequence.  These values should be positive
    !    integers between 1 and 30,000.
    !
    !    Output, real ( kind = 8 ) R8_RANDOM, the next value in the sequence.
    !
    implicit none

    integer ( kind = 4 ) s1
    integer ( kind = 4 ) s2
    integer ( kind = 4 ) s3
    real ( kind = 8 ) r8_random

    s1 = mod ( 171 * s1, 30269 )
    s2 = mod ( 172 * s2, 30307 )
    s3 = mod ( 170 * s3, 30323 )

    r8_random = mod ( real ( s1, kind = 8 ) / 30269.0D+00 &
                    + real ( s2, kind = 8 ) / 30307.0D+00 &
                    + real ( s3, kind = 8 ) / 30323.0D+00, 1.0D+00 )

    return
  end function r8_random
  !+-------------------------------------------------------------------+
  !| The end of the function r8_random.                                |
  !+-------------------------------------------------------------------+
  !!
end module initialisation
!+---------------------------------------------------------------------+
!| The end of the module initialisation.                               |
!+---------------------------------------------------------------------+