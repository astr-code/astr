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
    use commvar,  only : flowtype,lreadgrid
    use readwrite,only : readgrid,writegrid
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
      else
        print*,trim(flowtype),'not defined @ gridgen'
        stop ' !! error at gridgen' 
      endif
      !
      call writegrid
      !
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
    use parallel, only : ig0,jg0,kg0,lio
    use commarray,only : x
    use hdf5io
    !
    ! local data
    real(8) :: lx,ly,lz
    integer :: i,j,k
    !
    lx=50.d0
    ly=15.d0
    lz=15.d0
    !
    do k=0,km
    do j=0,jm
    do i=0,im
      !
      x(i,j,k,1)=lx/real(ia,8)*real(i+ig0,8)
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
  !+-------------------------------------------------------------------+
  !| This subroutine is to generate a cubic mesh.                      |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 12-02-2021  | Created by J. Fang @ Warrington                     |
  !+-------------------------------------------------------------------+
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