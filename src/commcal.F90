!+---------------------------------------------------------------------+
!| This module contains some common calculater.                        |
!+---------------------------------------------------------------------+
!| CHANGE RECORD                                                       |
!| -------------                                                       |
!| 19-03-2021  | Created by J. Fang                                    |
!+---------------------------------------------------------------------+
module commcal
  !
  use parallel, only: mpirank,mpistop
  !
  implicit none
  !
  contains
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to calculate CFL number and the           |
  !| corresponding time step.                                          |
  !+-------------------------------------------------------------------+
  !| ref: Adams, N. A., Shariff, K. 1996, J COMPUT PHYS. 127,27-51.    |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 21-03-2021: Created by J. Fang @ STFC Daresbury Laboratory        |
  !+-------------------------------------------------------------------+
  subroutine cflcal(deltat)
    !
    use commvar,  only: im,jm,km
    use commarray,only: vel,rho,prs,tmp,dxi,jacob
    use fludyna,  only: sos
    use parallel, only: pmax
    !
    ! arguments
    real(8),intent(in) :: deltat
    !
    ! local data
    real(8) :: cfl,ubar,vbar,wbar,css,csi,csj,csk
    integer :: i,j,k
    real(8) :: deltai,deltaj,deltak
    !
    deltai=0.d0
    deltaj=0.d0
    deltak=0.d0
    do k=0,km
    do j=0,jm
    do i=0,im
      ubar=dxi(i,j,k,1,1)*vel(i,j,k,1)+dxi(i,j,k,1,2)*vel(i,j,k,2)+    &
           dxi(i,j,k,1,3)*vel(i,j,k,3)
      vbar=dxi(i,j,k,2,1)*vel(i,j,k,1)+dxi(i,j,k,2,2)*vel(i,j,k,2)+    &
           dxi(i,j,k,2,3)*vel(i,j,k,3)
      wbar=dxi(i,j,k,3,1)*vel(i,j,k,1)+dxi(i,j,k,3,2)*vel(i,j,k,2)+    &
           dxi(i,j,k,3,3)*vel(i,j,k,3)
      !
      css=sos(tmp(i,j,k))
      csi=css*sqrt(dxi(i,j,k,1,1)**2+dxi(i,j,k,1,2)**2+dxi(i,j,k,1,3)**2)
      csj=css*sqrt(dxi(i,j,k,2,1)**2+dxi(i,j,k,2,2)**2+dxi(i,j,k,2,3)**2)
      csk=css*sqrt(dxi(i,j,k,3,1)**2+dxi(i,j,k,3,2)**2+dxi(i,j,k,3,3)**2)
      !
      deltai=max(ubar,ubar-csi,ubar+csi)
      deltaj=max(vbar,vbar-csj,vbar+csj)
      deltak=max(wbar,wbar-csk,wbar+csk)
      !
    enddo
    enddo
    enddo
    !
    deltai=pmax(deltai)
    deltaj=pmax(deltaj)
    deltak=pmax(deltak)
    !
    cfl=deltat*(deltai+deltaj+deltak)
    !
    if(mpirank==0) then
      write(*,"(A38)")'  =========== CFL Condition==========='
      write(*,"(A24,1x,F13.7)")'     current time step: ',deltat
      write(*,"(A24,1x,F13.7)")'           current CFL: ',cfl
      write(*,"(A24,1x,F13.7)")'   time step for CFL=1: ',deltat/cfl
      write(*,"(A38)")'  ===================================='
    end if
    !
  end subroutine cflcal
  !+-------------------------------------------------------------------+
  !| The end of the subroutine cflcal.                                 |
  !+-------------------------------------------------------------------+
  !
end module commcal
!+---------------------------------------------------------------------+
!| The end of the module commcal.                                      |
!+---------------------------------------------------------------------+