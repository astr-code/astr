module statistics
   
  use constdef

  implicit none

  contains

  subroutine stacal
    
    use comvardef, only: im,jm,km,rho,tmp,vel,dvel,nstep,time,ndims,reynolds
    use fluids,  only: miucal
    
    integer :: i,j,k
    real(rtype) :: var1,omega(3),tke,rhom,enst,dissipation
    real(rtype) :: s11,s12,s13,s22,s23,s33,miu
    
    logical,save :: firstcall = .true.
    integer,save :: filehand

    !$ save omega
    !$OMP THREADPRIVATE(omega)

    if(firstcall) then
      filehand=13
      open(filehand,file='state.dat')
      write(filehand,"(5(1X,A20))")'nstep','time','tke','enst','dissipation'
      firstcall=.false.
    endif
    
    rhom=0._rtype
    tke =0._rtype
    enst=0._rtype
    dissipation = 0._rtype

    if(ndims==3) then

      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k,var1,s11,s12,s13,s22,s23,s33,miu) &
      !$OMP REDUCTION(+:tke,rhom,enst,dissipation)
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
        
        s11=dvel(i,j,k,1,1)
        s12=0.5_rtype*(dvel(i,j,k,1,2)+dvel(i,j,k,2,1))
        s13=0.5_rtype*(dvel(i,j,k,1,3)+dvel(i,j,k,3,1))
        s22=dvel(i,j,k,2,2)
        s23=0.5_rtype*(dvel(i,j,k,2,3)+dvel(i,j,k,3,2))
        s33=dvel(i,j,k,3,3)
 
        miu=miucal(tmp(i,j,k))/reynolds

        var1=2._rtype*miu*(s11**2+s22**2+s33**2+2._rtype*(s12**2+s13**2+s23**2)-num1d3*(s11+s22+s33)**2)

        dissipation = dissipation + var1

      enddo
      enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
  
      rhom        = rhom/real(im*jm*km,rtype)
      tke         = 0.5_rtype*tke/real(im*jm*km,rtype)
      enst        = 0.5_rtype*enst/real(im*jm*km,rtype)
      dissipation = dissipation/real(im*jm*km,rtype)

    elseif(ndims==2) then
      k=0

      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,var1) REDUCTION(+:tke,rhom,enst)
      !$OMP DO
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
      !$OMP END DO
      !$OMP END PARALLEL
  
      rhom=rhom/real(im*jm,rtype)
      tke =0.5_rtype*tke/real(im*jm,rtype)
      enst=0.5_rtype*enst/real(im*jm,rtype)

    endif
    
    print*,' ** rho=',rhom,', tke=',tke,', enstrophy=',enst
    
    write(filehand,"(1X,I20,4(1X,E20.13E2))")nstep,time,tke,enst,dissipation
    
    if(mod(nstep,10)==0) flush(filehand)
    
  end subroutine stacal

end module statistics