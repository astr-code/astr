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
    use commvar,  only : flowtype,lreadgrid,limmbou,solidfile,nondimen
    use readwrite,only : readgrid,writegrid,readsolid
    !
    if(lreadgrid) then
      call readgrid
    else
      if(trim(flowtype)=='tgv') then
        if(nondimen) then
          call gridcube(2.d0*pi,2.d0*pi,2.d0*pi)
        else 
          call gridcube(5.1530915662d-3,5.1530915662d-3,5.1530915662d-3)
        endif 
      elseif(trim(flowtype)=='jet') then
        call gridjet
      elseif(trim(flowtype)=='2dvort') then
        call gridcube(20.d0,10.d0,1.d0)
      elseif(trim(flowtype)=='accutest') then
        call gridcube(10.d0,1.d0,1.d0)
      elseif(trim(flowtype)=='shuosher') then
        call grid1d(-5.d0,5.d0)
      elseif(trim(flowtype)=='windtunn') then
        call gridcube(20.d0,10.d0,10.d0)
      elseif(trim(flowtype)=='channel') then
        call grichan(2.d0*pi,pi)
      elseif(trim(flowtype)=='0dreactor') then
        call gridcube(1.d0,1.d0,1.d0)
      elseif(trim(flowtype)=='1dflame') then
        call gridcube(0.5d-2,0.05d-3,0.d0)
      elseif(trim(flowtype)=='h2supersonic') then
        call gridsupersonicjet
      else
        print*,trim(flowtype),' is not defined @ gridgen'
        stop ' !! error at gridgen' 
      endif
      !
      ! call writegrid
      !
    endif
    !
    if(limmbou) then
      if(mpirank==0) call readsolid(solidfile)
    endif
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
      ratio=1.d0
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
  end subroutine spongstretch
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! End of the subroutine spongstretch
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
    zmax=pmax(x(im,0,0,3))
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
end module gridgeneration
!+---------------------------------------------------------------------+
!| The end of the module gridgeneration.                              |
!+---------------------------------------------------------------------+