!+---------------------------------------------------------------------+
!| This module contains subroutines of calculating statistics.         |
!+---------------------------------------------------------------------+
!| CHANGE RECORD                                                       |
!| -------------                                                       |
!| 12-02-2021  | Created by J. Fang                                    |
!+---------------------------------------------------------------------+
module statistic
  !
  use commvar,  only : im,jm,km,ia,ja,ka,ndims,xmin,xmax,ymin,ymax,    &
                       zmin,zmax,nstep,deltat,force,numq
  use parallel, only : mpirank,lio,psum,pmax,pmin,bcast,mpistop
  !
  implicit none
  !
  real(8) :: enstophy,kenergy,fbcx,massflux,massflux_target,wrms
  real(8),allocatable :: max_q(:),min_q(:)
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
    elseif(trim(flowtype)=='channel') then
      fbcx=fbcxchan()
      massflux=massfluxchan()
      wrms=wrmscal()
      !
      if(nstep==0) then
        massflux_target=massflux
        force=0.d0
        !
        if(lio) print*,' ** target mass flux= ',massflux_target
        !
      endif
      !
      force(1)=chanfoce(force(1),massflux,fbcx,massflux_target)
      !
      call bcast(force)
    endif
    !
    call maxmincal
    !
  end subroutine statcal
  !+-------------------------------------------------------------------+
  !| The end of the subroutine statcal.                                |
  !+-------------------------------------------------------------------+
  !!
  !+-------------------------------------------------------------------+
  !| This subroutine is to return min and max variables.               |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 13-02-2021  | Created by J. Fang STFC Daresbury Laboratory        |
  !+-------------------------------------------------------------------+
  subroutine maxmincal
    !
    use commvar,   only : im,jm,km
    use commarray, only : q
    !
    ! local data
    integer :: n
    logical,save :: linit=.true.
    !
    if(linit) then
      allocate(max_q(1:numq),min_q(1:numq))
      linit=.false.
    endif
    !
    do n=1,numq
      max_q(n)=maxval(q(0:im,0:jm,0:km,n))
      min_q(n)=minval(q(0:im,0:jm,0:km,n))
      !
      max_q(n)=pmax(max_q(n))
      min_q(n)=pmin(min_q(n))
      !
    enddo
    !
    return
    !
  end subroutine maxmincal
  !+-------------------------------------------------------------------+
  !| The end of the subroutine maxmincal.                              |
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
    use commvar, only : flowtype,nstep,time,nlstep,maxstep,hand_fs
    !
    ! local data
    logical,save :: linit=.true.
    integer :: nprthead=10
    integer :: i
    ! 
    if(lio) then
      !
      if(linit) then
        open(hand_fs,file='flowstate.dat')
        !
        if(trim(flowtype)=='tgv') then
          write(hand_fs,"(A7,1X,A13,2(1X,A20))")'nstep','time',        &
                                                    'enstophy','kenergy'
        elseif(trim(flowtype)=='channel') then
          write(hand_fs,"(A7,1X,A13,4(1X,A20))")'nstep','time',        &
                                       'massflux','fbcx','forcex','wrms'
        else
          write(hand_fs,"(A7,1X,A13,5(1X,A20))")'nstep','time',        &
                                 'q1max','q2max','q3max','q4max','q5max'
        endif
        !
        linit=.false.
      endif
      !
      if(trim(flowtype)=='tgv') then
        write(hand_fs,"(I7,1X,E13.6E2,2(1X,E20.13E2))")nstep,time,     &
                                                        enstophy,kenergy
      elseif(trim(flowtype)=='channel') then
        write(hand_fs,"(I7,1X,E13.6E2,4(1X,E20.13E2))")nstep,time,     &
                                             massflux,fbcx,force(1),wrms
      else
        ! general flowstate
        write(hand_fs,"(I7,1X,E13.6E2,5(1X,E20.13E2))")nstep,time,     &
                                                        (max_q(i),i=1,5)
      endif
      !
      if(mod(nstep,nlstep)==0) then
        !
        if(mod(nstep,nprthead*nlstep)==0) then
          if(trim(flowtype)=='tgv') then
            write(*,"(2X,A7,3(1X,A13))")'nstep','time','enstophy','kenergy'
          elseif(trim(flowtype)=='channel') then
            write(*,"(2X,A7,5(1X,A13))")'nstep','time','massflux',     &
                                                  'fbcx','forcex','wrms'
          else
            write(*,"(2X,A7,6(1X,A13))")'nstep','time',                &
                                 'q1max','q2max','q3max','q4max','q5max'
          endif
          write(*,'(2X,76A)')('-',i=1,76)
        endif
        !
        if(trim(flowtype)=='tgv') then
          write(*,"(2X,I7,3(1X,E13.6E2))")nstep,time,enstophy,kenergy
        elseif(trim(flowtype)=='channel') then
          write(*,"(2X,I7,5(1X,E13.6E2))")nstep,time,massflux,fbcx,    &
                                          force(1),wrms
        else
          write(*,"(2X,I7,6(1X,E13.6E2))")nstep,time,(max_q(i),i=1,5)
        endif
        !
        flush(hand_fs)
      endif
      !
      if(nstep==maxstep) then
        close(hand_fs)
        print*,' << flowstate.dat'
      endif
      !
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
      !
      !
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
  !| The end of the subroutine kenergycal.                             |
  !+-------------------------------------------------------------------+
  !!
  !+-------------------------------------------------------------------+
  !| This function is to return rms spanwise velocity fluctuation.     |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 12-02-2021  | Created by J. Fang STFC Daresbury Laboratory        |
  !+-------------------------------------------------------------------+
  function wrmscal() result(vout)
    !
    use commvar,   only : im,jm,km,ia,ja,ka
    use commarray, only : vel
    !
    real(8) :: vout
    !
    ! local data
    integer :: i,j,k,k1,k2
    real(8) :: norm
    !
    if(ndims==2) then
      k1=0
      k2=0
      norm=real(ia*ja,8)
    elseif(ndims==3) then
      k1=1
      k2=km
      norm=real(ia*ja*ka,8)
    else
      print*,' !! ndims=',ndims
      stop ' !! error @ fbcxchan !!'
    endif
    !
    vout=0.d0
    do k=1,km
    do j=1,jm
    do i=1,im
      vout=vout+vel(i,j,k,3)*vel(i,j,k,3)
    enddo
    enddo
    enddo
    !
    vout=sqrt(psum(vout)/norm)
    !
    return
    !
  end function wrmscal
  !+-------------------------------------------------------------------+
  !| The end of the subroutine wrmscal.                                |
  !+-------------------------------------------------------------------+
  !!
  !+-------------------------------------------------------------------+
  !| This function is to return wall skin-friction of channel.         |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 13-02-2021  | Created by J. Fang STFC Daresbury Laboratory        |
  !+-------------------------------------------------------------------+
  real(8) function fbcxchan()
    !
    use commvar,  only : reynolds
    use commarray,only : tmp,dvel
    use parallel, only : jrk,jrkm
    use fludyna,  only : miucal
    !
    ! local data
    integer :: i,j,k,k1,k2
    real(8) :: miu,norm
    !
    fbcxchan=0.d0
    !
    if(ndims==2) then
      k1=0
      k2=0
      norm=real(ia,8)
    elseif(ndims==3) then
      k1=1
      k2=km
      norm=real(ia*ka,8)
    else
      print*,' !! ndims=',ndims
      stop ' !! error @ fbcxchan !!'
    endif
    !
    if(jrk==0) then
      j=0
      do k=k1,k2
      do i=1,im
        miu=miucal(tmp(i,j,k))/reynolds
        !
        fbcxchan=fbcxchan+miu*dvel(i,j,k,1,2)
      enddo
      enddo
      !
    endif
    !
    if(jrk==jrkm) then
      !
      j=jm
      do k=k1,k2
      do i=1,im
        miu=miucal(tmp(i,j,k))/reynolds
        !
        fbcxchan=fbcxchan-miu*dvel(i,j,k,1,2)
      enddo
      enddo
    endif
    !
    fbcxchan=psum(fbcxchan)/norm
    !
    return
    !
  end function fbcxchan
  !+-------------------------------------------------------------------+
  !| The end of the subroutine fbcxchan.                               |
  !+-------------------------------------------------------------------+
  !!
  !+-------------------------------------------------------------------+
  !| This function is to return mass flux of channel flow.             |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 13-02-2021  | Created by J. Fang STFC Daresbury Laboratory        |
  !+-------------------------------------------------------------------+
  real(8) function massfluxchan()
    !
    use commarray,only : x,q
    !
    ! local data
    integer :: i,j,k,k1,k2
    real(8) :: dy,var1,norm
    !
    massfluxchan=0.d0
    !
    if(ndims==2) then
      k1=0
      k2=0
      norm=real(ia,8)
    elseif(ndims==3) then
      k1=1
      k2=km
      norm=real(ia*ka,8)
    else
      print*,' !! ndims=',ndims
      stop ' !! error @ massfluxchan !!'
    endif
    !
    do k=k1,k2
    do j=1,jm
    do i=1,im
      dy=x(i,j,k,2)-x(i,j-1,k,2)
      var1=0.5d0*(q(i,j,k,2)+q(i,j-1,k,2))
      !
      massfluxchan=massfluxchan+var1*dy
    enddo
    enddo
    enddo
    !
    massfluxchan=psum(massfluxchan)/norm
    !
    return
    !
  end function massfluxchan
  !+-------------------------------------------------------------------+
  !| The end of the subroutine massfluxchan.                           |
  !+-------------------------------------------------------------------+
  !!
  !+-------------------------------------------------------------------+
  !| This subroutine is used to calcualte the body force in the channel|
  !| to main the target mass flux.                                     |
  !+-------------------------------------------------------------------+
  !| ref: Lenormand, E., Sagaut, P.  and Phuoc, L.T., Large eddy       |
  !|      simulation of subsonic and supersonic channel flow at        |
  !|      moderate Reynolds number. international journal for numerical|
  !|       methods in fluids, 2000, 32: 369-406.                       |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 22-Mar-2019  | Created by J. Fang STFC Daresbury Laboratory       |
  !+-------------------------------------------------------------------+
  real(8) function chanfoce(force,massflux,friction,massfluxtarget)
    ! arguments
    real(8),intent(in) :: force,massflux,friction,massfluxtarget
    !
    ! local data
    real(8) :: gn,qn1,ly
    !
    ly=(ymax-ymin)
    !
    if(nstep==0) then
      chanfoce=friction/ly
    else
      gn=(ly*force-friction)
      qn1=massflux+deltat*gn
      chanfoce=force-(2.d0*(qn1-massfluxtarget)-                       &
               0.2d0*(massflux-massfluxtarget))/ly
    endif
    !
    ! if(lio) print*,massflux,massfluxtarget,'|',qn1,force,chanfoce
    !
  end function chanfoce
  !+-------------------------------------------------------------------+
  !| The end of the subroutine chanfoce.                               |
  !+-------------------------------------------------------------------+
  !!
end module statistic
!+---------------------------------------------------------------------+
!| The end of the module statistic.                                    |
!+---------------------------------------------------------------------+