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
    use commvar,  only : flowtype,lreadgrid,limmbou,solidfile
    use readwrite,only : readgrid,writegrid,readsolid
    !
    if(lreadgrid) then
      call readgrid
    else
      if(trim(flowtype)=='tgv') then
        call gridcube(2.d0*pi,2.d0*pi,2.d0*pi)
      elseif(trim(flowtype)=='jet') then
        call gridjet
      elseif(trim(flowtype)=='2dvort') then
        call gridcube(20.d0,10.d0,1.d0)
      elseif(trim(flowtype)=='accutest') then
        call gridcube(10.d0,1.d0,1.d0)
      elseif(trim(flowtype)=='shuosher') then
        call grid1d(-5.d0,5.d0)
      else
        print*,trim(flowtype),'not defined @ gridgen'
        stop ' !! error at gridgen' 
      endif
      !
      call writegrid
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
    if(lio) print*,' ** cubic grid generated'
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
        x(i,j,k,3)=lz
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
  !
end module gridgeneration
!+---------------------------------------------------------------------+
!| The end of the module gridgeneration.                              |
!+---------------------------------------------------------------------+