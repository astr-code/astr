!+---------------------------------------------------------------------+
!| This module contains subroutines of calculating statistics.         |
!+---------------------------------------------------------------------+
!| CHANGE RECORD                                                       |
!| -------------                                                       |
!| 12-02-2021  | Created by J. Fang                                    |
!+---------------------------------------------------------------------+
module statistic
  !
  use parallel, only : mpirank,lio,psum
  !
  implicit none
  !
  real(8) :: enstophy,kenergy
  !
  contains
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to calculate and output instantous status.|
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 12-02-2021  | Created by J. Fang STFC Daresbury Laboratory        |
  !+-------------------------------------------------------------------+
  subroutine statcal
    !
    use commvar, only: flowtype
    !
    ! local data
    !
    if(trim(flowtype)=='tgv') then
      enstophy=enstophycal()
      kenergy =kenergycal()
    endif
    !
  end subroutine statcal
  !+-------------------------------------------------------------------+
  !| The end of the subroutine statcal.                                |
  !+-------------------------------------------------------------------+
  !!
  !+-------------------------------------------------------------------+
  !| This subroutine is used to print state on screen and write to a   |
  !| file.                                                             |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 12-02-2021  | Created by J. Fang STFC Daresbury Laboratory        |
  !+-------------------------------------------------------------------+
  subroutine statout
    !
    use commvar, only : nstep,time,nlstep,maxstep
    !
    ! local data
    logical,save :: linit=.true.
    integer :: nprthead=10
    integer :: i
    ! 
    if(lio) then
      !
      if(linit) then
        open(13,file='flowstate.dat')
        write(13,"(A7,1X,A13,2(1X,A20))")'nstep','time','enstophy','kenergy'
        linit=.false.
      endif
      !
      write(13,"(I7,1X,E13.6E2,2(1X,E20.13E2))")nstep,time,enstophy,kenergy
      !
      if(nstep==maxstep) then
        close(13)
        print*,' << flowstate.dat'
      endif
      !
      if(mod(nstep,nlstep)==0) then
        !
        if(mod(nstep,nprthead*nlstep)==0) then
          write(*,"(2X,A7,3(1X,A13))")'nstep','time','enstophy','kenergy'
           write(*,'(2X,62A)')('-',i=1,62)
        endif
        !
        write(*,"(2X,I7,3(1X,E13.6E2))")nstep,time,enstophy,kenergy
        !
        flush(13)
      endif
    endif
    !
  end subroutine statout
  !+-------------------------------------------------------------------+
  !| The end of the subroutine statprint.                              |
  !+-------------------------------------------------------------------+
  !!
  !+-------------------------------------------------------------------+
  !| This function is to return spatial averaged enstrophy.            |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 12-02-2021  | Created by J. Fang STFC Daresbury Laboratory        |
  !+-------------------------------------------------------------------+
  function enstophycal() result(vout)
    !
    use commvar,   only : im,jm,km,ia,ja,ka
    use commarray, only : vel,celvol,rho,dvel
    use commfunc,  only : volhex
    !
    real(8) :: vout
    !
    ! local data
    integer :: i,j,k
    real(8) :: omega(3),omegam
    !
    vout=0.d0
    do k=1,km
    do j=1,jm
    do i=1,im
      ! dx=
      omega(1)=dvel(i,j,k,3,2)-dvel(i,j,k,2,3)
      omega(2)=dvel(i,j,k,1,3)-dvel(i,j,k,3,1)
      omega(3)=dvel(i,j,k,2,1)-dvel(i,j,k,1,2)
      omegam=omega(1)*omega(1)+omega(2)*omega(2)+omega(3)*omega(3)
      !
      vout=vout+rho(i,j,k)*omegam
    enddo
    enddo
    enddo
    !
    vout=0.5d0*psum(vout)/real(ia*ja*ka,8)
    !
    return
    !
  end function enstophycal
  !+-------------------------------------------------------------------+
  !| The end of the subroutine enstophycal.                            |
  !+-------------------------------------------------------------------+
  !!
  !+-------------------------------------------------------------------+
  !| This function is to return spatial averaged kinetic energy.       |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 12-02-2021  | Created by J. Fang STFC Daresbury Laboratory        |
  !+-------------------------------------------------------------------+
  function kenergycal() result(vout)
    !
    use commvar,   only : im,jm,km,ia,ja,ka
    use commarray, only : vel,celvol,rho,dvel
    use commfunc,  only : volhex
    !
    real(8) :: vout
    !
    ! local data
    integer :: i,j,k
    real(8) :: var1
    !
    vout=0.d0
    do k=1,km
    do j=1,jm
    do i=1,im
      ! 
      var1=vel(i,j,k,1)**2+vel(i,j,k,2)**2+vel(i,j,k,3)**2
      !
      vout=vout+rho(i,j,k)*var1
    enddo
    enddo
    enddo
    !
    vout=0.5d0*psum(vout)/real(ia*ja*ka,8)
    !
    return
    !
  end function kenergycal
  !+-------------------------------------------------------------------+
  !| The end of the subroutine kenergycal.                            |
  !+-------------------------------------------------------------------+
  !!
end module statistic
!+---------------------------------------------------------------------+
!| The end of the module statistic.                                    |
!+---------------------------------------------------------------------+