!+---------------------------------------------------------------------+
!| This module contains subroutines of generating grid.                |
!+---------------------------------------------------------------------+
!| CHANGE RECORD                                                       |
!| -------------                                                       |
!| 12-02-2021  | Created by J. Fang                                    |
!+---------------------------------------------------------------------+
module gridgeneration
  !
  use constdef
  !
  implicit none
  !
  contains
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is the main entrance of the grid generation.      |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 12-02-2021  | Created by J. Fang @ Warrington                     |
  !+-------------------------------------------------------------------+
  subroutine gridgen
    !
    use parallel, only : mpirank
    use commvar,  only : flowtype,lreadgrid,nondimen,gridfile
    use readwrite,only : readgrid,writegrid,xdmfwriter
    !
    if(lreadgrid) then
      call readgrid(trim(gridfile))
    else
      if(flowtype(1:3)=='tgv') then
        if(nondimen) then
          call gridcube(2.d0*pi,2.d0*pi,2.d0*pi)
        elseif(trim(flowtype)=='tgv') then
          call gridcube(5.1530915662d-3,5.1530915662d-3,5.1530915662d-3)
        else
          call gridcube(6.283185307d-3,6.283185307d-3,6.283185307d-3)
        endif 
      elseif(trim(flowtype)=='jet') then
        call gridjet
      elseif(trim(flowtype)=='hit') then
        call gridcube(2.d0*pi,2.d0*pi,2.d0*pi)
      elseif(trim(flowtype)=='2dvort') then
        call gridcube(20.d0,10.d0,1.d0)
      elseif(trim(flowtype)=='accutest') then
        call gridcube(10.d0,1.d0,1.d0)
      elseif(trim(flowtype)=='shuosher') then
        call grid1d(-5.d0,5.d0)
      elseif(trim(flowtype)=='windtunn') then
        call gridsandbox
      elseif(trim(flowtype)=='channel') then
        call grichan(2.d0*pi,pi)
      elseif(trim(flowtype)=='0dreactor') then
        call gridcube(1.d0,1.d0,1.d0)
      elseif(trim(flowtype)=='1dflame') then
        call gridcube(0.5d-2,0.05d-3,0.d0)
      elseif(trim(flowtype)=='h2supersonic') then
        call gridsupersonicjet
      elseif(trim(flowtype)=='rti') then
        call gridcube(0.25d0,1.d0,0.d0)
      else
        print*,trim(flowtype),' is not defined @ gridgen'
        stop ' !! error at gridgen' 
      endif
      !
      call writegrid(trim(gridfile))
      !
    endif
    !
    call xdmfwriter(gridh5file=trim(gridfile),mode='grid')
    !
  end subroutine gridgen
  !+-------------------------------------------------------------------+
  !| The end of the subroutine gridgen.                                |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to generate a mesh for jet flow.               |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 24-02-2021  | Created by J. Fang @ Warrington                     |
  !+-------------------------------------------------------------------+
  subroutine gridjet
  !
    use commvar,  only : im,jm,km,gridfile,ia,ja,ka
    use parallel, only : ig0,jg0,kg0,lio,irk,irkm
    use commarray,only : x
    use hdf5io
    !
    ! local data
    real(8) :: lx,ly,lz
    integer :: i,j,k,in
    real(8) :: x1d(0:im)
    !
    lx=50.d0
    ly=15.d0
    lz=15.d0
    !
    do i=0,im
      x1d(i)=lx/real(ia,8)*real(i+ig0,8)
    enddo
    !
    if(irk==irkm) then
      in=im-25
      call spongstretch(im-in,lx+10.d0,x1d(in-2:im))
    endif
    !
    do k=0,km
    do j=0,jm
    do i=0,im
      !
      x(i,j,k,1)=x1d(i)
      x(i,j,k,2)=ly/real(ja,8)*real(j+jg0,8)-0.5d0*ly
      if(ka==0) then
        x(i,j,k,3)=0.d0 !lz/real(ka,8)*real(k+kg0,8)
      else
        x(i,j,k,3)=lz/real(ka,8)*real(k+kg0,8)-0.5d0*lz
      endif
      !
    enddo
    enddo
    enddo
    !
    if(lio) print*,' ** jet grid generated'
    !
  end subroutine gridjet
  !+-------------------------------------------------------------------+
  !| The end of the subroutine gridjet.                                |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to generate a cubic mesh.                      |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 12-02-2021  | Created by J. Fang @ Warrington                     |
  !+-------------------------------------------------------------------+
  subroutine grid1d(xmin,xmax)
    !
    use commvar,  only : im,jm,km,ia,ja,ka
    use parallel, only : ig0,jg0,kg0,lio
    use commarray,only : x
    use hdf5io
    !
    ! arguments
    real(8),intent(in) :: xmin,xmax
    !
    ! local data
    integer :: i,j,k
    !
    do k=0,km
    do j=0,jm
    do i=0,im
      x(i,j,k,1)=(xmax-xmin)/dble(ia)*dble(i+ig0)+xmin
      if(ja==0) then
        x(i,j,k,2)=0.d0
      else
        x(i,j,k,2)=1.d0/real(ja,8)*real(j+jg0,8)
      endif
      x(i,j,k,3)=0.d0
      !
    enddo
    enddo
    enddo
    !
    if(lio) print*,' ** 1-D grid generated'
    !
  end subroutine grid1d
  !
  subroutine gridcube(lx,ly,lz)
    !
    use commvar,  only : im,jm,km,gridfile,ia,ja,ka
    use parallel, only : ig0,jg0,kg0,lio
    use commarray,only : x
    use hdf5io
    !
    ! arguments
    real(8),intent(in) :: lx,ly,lz
    !
    ! local data
    integer :: i,j,k
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
  end subroutine gridcube
  !+-------------------------------------------------------------------+
  !| The end of the subroutine gridcube.                               |
  !+-------------------------------------------------------------------+
  !
  subroutine grichan(lx,lz)
    !
    use commvar,  only : im,jm,km,gridfile,ia,ja,ka
    use parallel, only : ig0,jg0,kg0,lio,pmax,mpistop
    use commarray,only : x
    use commfunc, only : argtanh
    !
    real(8),intent(in) :: lx,lz
    !
    integer :: n,i,j,k
    real(8) :: varc,var1,var2,Retau,dx,yy(0:jm),xmax,ymax,dymax,zmax
    !
    Retau=185.d0
    !
    ! varc=1.02d0
    varc=1.07d0
    do j=0,jm
      !
      var1=argtanh(1.d0/varc)
      var2=2.d0*(j+jg0)/(ja)*1.d0-1.d0
      !
      yy(j)=1.d0*(1.d0+varc*dtanh(var1*var2))
      !
    end do
    !
    do k=0,km
    do j=0,jm
    do i=0,im
      !
      x(i,j,k,1)=lx/real(ia,8)*real(i+ig0,8)
      !
      x(i,j,k,2)=yy(j)
      !
      if(ka==0) then
        x(i,j,k,3)=0.d0
      else
        x(i,j,k,3)=lz/real(ka,8)*real(k+kg0,8)
      endif
      !
    end do
    end do
    end do
    !
    do j=1,jm
      dymax=max(dymax,yy(j)-yy(j-1))
    enddo
    !
    xmax=pmax(x(im,0,0,1))
    ymax=pmax(yy(jm))
    dymax=pmax(dymax)
    zmax=pmax(x(0,0,km,3))
    !
    if(lio) then
      !
      write(*,'(2X,62A)')('-',i=1,62)
      write(*,'(18X,A)')' ** channel grid resolution **'
      write(*,'(2(A,F13.7))')'  ** y1+=',yy(1)*Retau,                  &
                             '   ymax+=',dymax*Retau
      write(*,'(2(A,F13.7))')'  ** dx+=',(x(1,0,0,1)-x(0,0,0,1))*Retau,&
                               '     lx+=',xmax*Retau
      write(*,'(2(A,F13.7))')'  ** dz+=',(x(0,0,1,3)-x(0,0,0,3))*Retau,&
                               '     lz+=',zmax*Retau
      write(*,'(2X,62A)')('-',i=1,62)
      !
      write(*,'(A)')'  ** channel grid generated'
      !
    endif
    !
  end subroutine grichan
  !+-------------------------------------------------------------------+
  !| The end of the subroutine grichan.                                |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to generate grid for jet simulation.           |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 04-Jan-2020: Created by J. Fang @ Warrington.                     |
  !+-------------------------------------------------------------------+
  subroutine gridsupersonicjet
    !
#ifdef COMB
    use commvar, only: dj_i,dco_i,im,jm,km,ia,ja,ka
    use parallel, only : ig0,jg0,kg0
    use commarray,only : x
    !
    ! local data
    real(8),allocatable :: yn(:),xn(:),rn(:)
    real(8) :: varc,var1,var2,varmin,varmax,dy,xmax,ymax,zmax,scale
    integer :: i,j,k,in,nd
    !
    !+--------+
    !| begin  |
    !+--------+
    xmax=118.d-3
    ymax=70.8d-3
    zmax=70.8d-3
    ! NOTE: im, jm ,zm must be even numbers
    ! xmax=44.0d-3
    ! ymax=40.0d-3
    ! zmax=40.0d-3
    !
    allocate(yn(0:ja),xn(0:ia),rn(0:ja))
    !
    ! number of cells in a diameter
    nd=8*24
    ! dy=16.d-3/real(nd,8)
    dy=8.d0*dj_i/real(nd,8)
    in=nd/2
    do j=1,in
      rn(j)=dy*j
    enddo
    !
    call spongstretch(ja/2-in,ymax/2,rn(in-2:ja/2))
    !
    yn(:)=0.d0
    do j=1,ja/2
      yn(j)=yn(j-1)+(rn(ja/2-j+1)-rn(ja/2-j))
    enddo
    !
    do j=ja/2+1,ja
      yn(j)=yn(j-1)+(rn(j-ja/2)-rn(j-ja/2-1))
    enddo
    !
    ! in=nd+22*24
    ! do i=0,in
    !   xn(i)=dy*i
    ! enddo
    ! call spongstretch(ia-in,xmax,xn(in-2:ia))
    do i=0,ia
      xn(i)=xmax/ia*i
    enddo
    !
    !
    do k=0,km
    do j=0,jm
    do i=0,im
      x(i,j,k,1)=xn(i+ig0) 
      x(i,j,k,2)=yn(j+jg0)
      x(i,j,k,3)=yn(k+kg0)
      ! z(i,j,k)=(zmax/km)*k
    end do
    end do
    end do
    !+--------+
    !| end    |
    !+--------+
    !
    return
    !
#endif
  end subroutine gridsupersonicjet
  !+-------------------------------------------------------------------+
  !| The end of the subroutine gridjet.                                |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to generate a sandbox grid for wind tunnel.    |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 23-Jun-2022: Created by J. Fang @ Warrington.                     |
  !+-------------------------------------------------------------------+
  subroutine gridsandbox
    !
    use parallel,  only : mpistop,ig0,jg0,kg0,lio,mpirank
    use commvar,   only : immbody,nsolid,ia,ja,ka,im,jm,km,ndims
    use commarray, only : x
    use stlaio
    use tecio
    !
    integer :: imm,jmm,kmm,im4,jm4,km4,isp
    integer :: i,j,k,jbody
    real(8) :: lx_solid,ly_solid,lz_solid,len_chara,x0,y0,z0
    real(8) :: xmin,ymin,zmin,xmax,ymax,zmax
    real(8),allocatable,dimension(:) :: xx,yy,zz
    real(8),allocatable,dimension(:,:,:) :: xa,ya,za
    !
    xmin=1.d10; xmax=-1.d10
    ymin=1.d10; ymax=-1.d10
    zmin=1.d10; zmax=-1.d10
    do jbody=1,size(immbody)
      xmin=min(xmin,immbody(jbody)%xmin(1))
      ymin=min(ymin,immbody(jbody)%xmin(2))
      zmin=min(zmin,immbody(jbody)%xmin(3))
      !
      xmax=max(xmax,immbody(jbody)%xmax(1))
      ymax=max(ymax,immbody(jbody)%xmax(2))
      zmax=max(zmax,immbody(jbody)%xmax(3))
    enddo
    !
    lx_solid=xmax-xmin
    ly_solid=ymax-ymin
    lz_solid=zmax-zmin
    !
    len_chara=1.d0 !sqrt(lx_solid**2+ly_solid**2)
    ! lx_solid=2.091906d0
    ! ly_solid=2.d0
    ! lz_solid=8.291200d-2
    !
    if(mpirank==0) print*,' **           solid size:',lx_solid,ly_solid,lz_solid
    if(mpirank==0) print*,' **  characteristic size:',len_chara
    if(mpirank==0) print*,' **            mesh size:',ia,ja,ka
    !
    lx_solid=lx_solid
    ly_solid=ly_solid
    lz_solid=lz_solid
    !
    allocate(xx(0:ia),yy(0:ja),zz(0:ka))
    !
    isp=20
    !
    imm=(ia-isp)/2
    jmm=ja/2
    kmm=ka/2
    !
    im4=imm/2
    jm4=jmm/2
    km4=kmm/2
    !
    ! core region 
    do i=imm,imm+im4
      xx(i)=lx_solid/dble(im4)*dble(i-imm)
    enddo
    !
    ! effecive region
    call spongstretch(imm-im4,5.d0*lx_solid,xx(imm+im4-2:imm+imm))
    !
    do i=0,imm-1
      xx(i)=-xx(imm+imm-i)
    enddo
    !
    ! sponge region
    call spongstretch(ia-2*imm,10.d0*lx_solid,xx(imm+imm-2:ia))
    !
    x0=xx(0)
    xx=xx-x0
    !
    if(mpirank==0) then
      print*,' ** -------------------- mesh report --------------------'
      print*,'    x direction'
      print*,'       core zone '
      write(*,'(2(A,1X,F12.8))')'          x:',xx(imm-im4),' ~ ',xx(imm+im4)
      write(*,'(2(A,1X,I0))')   '          i:',imm-im4,' ~ ', imm+im4
      write(*,'(1(A,1X,F12.8))')'         dx:',xx(imm+im4)-xx(imm+im4-1)
      print*,'       effecive zone '
      write(*,'(2(A,1X,F12.8))')'          x:',xx(0),' ~ ',xx(imm+imm)
      write(*,'(2(A,1X,I0))')   '          i:',0,' ~ ', imm+imm
      write(*,'(1(A,1X,F12.8))')'         dx:',xx(imm+imm)-xx(imm+imm-1)
      print*,'       sponge zone '
      write(*,'(2(A,1X,F12.8))')'          x:',xx(imm+imm),' ~ ',xx(ia)
      write(*,'(2(A,1X,I0))')   '          i:',imm+imm,' ~ ', ia
      write(*,'(1(A,1X,F12.8))')'         dx:',xx(ia)-xx(ia-1)
      print*,' ** -----------------------------------------------------'
    endif
    !
    ! core regjon 
    do j=jmm,jmm+jm4
      yy(j)=ly_solid/dble(jm4)*dble(j-jmm)
    enddo
    !
    ! effecive regjon
    call spongstretch(jmm-jm4,5.d0*ly_solid,yy(jmm+jm4-2:jmm+jmm))
    !
    do j=0,jmm-1
      yy(j)=-yy(jmm+jmm-j)
    enddo
    y0=yy(0)
    yy=yy-y0
    !
    if(mpirank==0) then
      print*,'    y djrectjon'
      print*,'       core zone '
      write(*,'(2(A,1x,F12.8))')'          y:',yy(jmm-jm4),' ~ ',yy(jmm+jm4)
      write(*,'(2(A,1x,i0))')   '          j:',jmm-jm4,' ~ ', jmm+jm4
      write(*,'(1(A,1x,F12.8))')'         dy:',yy(jmm+jm4)-yy(jmm+jm4-1)
      print*,'       effecive zone '
      write(*,'(2(A,1x,F12.8))')'          y:',yy(0),' ~ ',yy(jmm+jmm)
      write(*,'(2(A,1x,i0))')   '          j:',0,' ~ ', jmm+jmm
      write(*,'(1(A,1x,F12.8))')'         dy:',yy(jmm+jmm)-yy(jmm+jmm-1)
      print*,' ** -----------------------------------------------------'
    endif
    !
    if(ndims==3) then
      !
      ! core regkon 
      do k=kmm,kmm+km4
        zz(k)=lz_solid/dble(km4)*dble(k-kmm)
      enddo
      !
      ! effeckve regkon
      ! call spongstretch(kmm-km4,5.d0*lz_solid,zz(kmm+km4-2:kmm+kmm))
      call spongstretch(kmm-km4,5.d0,zz(kmm+km4-2:kmm+kmm))
      !
      do k=0,kmm-1
        zz(k)=-zz(kmm+kmm-k)
      enddo
      if(mpirank==0) then
        print*,'    z dkrectkon'
        print*,'       core zone '
        write(*,'(2(A,1x,F12.8))')'          z:',zz(kmm-km4),' ~ ',zz(kmm+km4)
        write(*,'(2(A,1x,i0))')   '          k:',kmm-km4,' ~ ', kmm+km4
        write(*,'(1(A,1x,F12.8))')'         dz:',zz(kmm+km4)-zz(kmm+km4-1)
        print*,'       effecive zone '
        write(*,'(2(A,1x,F12.8))')'          z:',zz(0),' ~ ',zz(kmm+kmm)
        write(*,'(2(A,1x,i0))')   '          k:',0,' ~ ', kmm+kmm
        write(*,'(1(A,1x,F12.8))')'         dz:',zz(kmm+kmm)-zz(kmm+kmm-1)
        print*,' ** -----------------------------------------------------'
      endif
      !
      z0=zz(0)
      zz=zz-z0
      !
    else
      zz=0.d0
    endif
    !
    if(mpirank==0) print*,' ** center of the domain:',xx(imm),yy(jmm),zz(kmm)
    !
    do k=0,km
    do j=0,jm
    do i=0,im
      x(i,j,k,1)=xx(i+ig0)
      x(i,j,k,2)=yy(j+jg0)
      if(ndims==3) then
        x(i,j,k,3)=zz(k+kg0)
      else
        x(i,j,k,3)=0.d0
      endif
    enddo
    enddo
    enddo
    !
    if(mpirank==0) then
      print*,' ** sandbox mesh generated, ranging:'
      write(*,'(2(A,1x,F12.8))')'          x:',xx(0),' ~ ',xx(ia)
      write(*,'(2(A,1x,F12.8))')'          y:',yy(0),' ~ ',yy(ja)
      write(*,'(2(A,1x,F12.8))')'          z:',zz(0),' ~ ',zz(ka)
    endif
    ! allocate(xa(0:ia,0:ja,0:ka),ya(0:ia,0:ja,0:ka),za(0:ia,0:ja,0:ka))
    ! do k=0,ka
    ! do j=0,ja
    ! do i=0,ia
    !   xa(i,j,k)=xx(i)
    !   ya(i,j,k)=yy(j)
    !   za(i,j,k)=zz(k)
    ! enddo
    ! enddo
    ! enddo
    ! !
    ! call tecbin('tecgrid.plt',xa,'x',ya,'y',za,'z')
    !
  end subroutine gridsandbox
  !+-------------------------------------------------------------------+
  !| The end of the subroutine gridjet.                                |
  !+-------------------------------------------------------------------+
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! This subroutine is used to stretch grid for spong layers
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine spongstretch(dim,varmax,var)
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
  end subroutine spongstretch
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! End of the subroutine spongstretch
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
end module gridgeneration
!+---------------------------------------------------------------------+
!| The end of the module gridgeneration.                              |
!+---------------------------------------------------------------------+
