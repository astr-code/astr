module statistics
   
  use constdef

  implicit none

  contains

  subroutine stacal
    !
    use comvardef, only: im,jm,km,rho,vel,dvel,nstep,time
    !
    integer :: i,j,k
    real(8) :: var1,omega(3),tke,rhom,enst
    !
    logical,save :: firstcall = .true.
    integer,save :: filehand

    !$ save omega
    !$OMP THREADPRIVATE(omega)

    if(firstcall) then
      filehand=13
      open(filehand,file='state.dat')
      firstcall=.false.
    endif
    !
    rhom=0.d0
    tke =0.d0
    enst=0.d0

    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k,var1) REDUCTION(+:tke,rhom,enst)
    !$OMP DO
    do k=1,km
    do j=1,jm
    do i=1,im
      ! 
      rhom=rhom+rho(i,j,k)
      !
      var1=vel(i,j,k,1)**2+vel(i,j,k,2)**2+vel(i,j,k,3)**2
      !
      tke=tke+rho(i,j,k)*var1
      !
      omega(1)=dvel(i,j,k,3,2)-dvel(i,j,k,2,3)
      omega(2)=dvel(i,j,k,1,3)-dvel(i,j,k,3,1)
      omega(3)=dvel(i,j,k,2,1)-dvel(i,j,k,1,2)
      !
      enst=enst+rho(i,j,k)*(omega(1)**2+omega(2)**2+omega(3)**2)
      !
    enddo
    enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL

    rhom=rhom/dble(im*jm*km)
    tke =0.5d0*tke/dble(im*jm*km)
    enst=0.5d0*enst/dble(im*jm*km)
    !
    ! print*,' ** rho=',rhom,', tke=',tke,', enstrophy=',enst
    !
    write(filehand,*)nstep,time,tke,enst
    !
    if(mod(nstep,10)==0) flush(filehand)
    !
  end subroutine stacal

end module statistics