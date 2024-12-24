module numericalflux
  
  use constdef

  implicit none

  integer :: lmax=4
  real(rtype) :: dcoe(4,4)
  
  contains

  subroutine flux_solver_init

    dcoe(1,1) =  1._rtype/2._rtype

    dcoe(1,2) =  2._rtype/3._rtype
    dcoe(2,2) = -1._rtype/12._rtype

    dcoe(1,3) =  3._rtype/4._rtype
    dcoe(2,3) = -3._rtype/20._rtype
    dcoe(3,3) =  1._rtype/60._rtype

    dcoe(1,4) =  4._rtype/5._rtype
    dcoe(2,4) = -1._rtype/5._rtype
    dcoe(3,4) =  4._rtype/105._rtype
    dcoe(4,4) = -1._rtype/280._rtype 

  end subroutine flux_solver_init

  
  function flux_sp(dens,velo,pres,enth,dim,dir) result(flux)

    use comvardef, only : hm

    integer,intent(in) :: dim
    character(len=1),intent(in) :: dir
    real(rtype),intent(in) :: dens(-hm:dim+hm),velo(-hm:dim+hm,1:3), &
                              pres(-hm:dim+hm),enth(-hm:dim+hm)
    real(rtype) :: flux(-1:dim,1:5)

    ! local data
    integer :: i,l,m,n
    real(rtype) :: uv(1:6,1:lmax,-hm:dim)
    real(rtype) :: fh(5),ft(6),uvs(6)
    real(rtype) :: u(1:6,-hm:dim+hm),v(1:6,-hm:dim+hm)
    real(rtype) :: rhom
    integer :: ip,iv

    if(dir=='i') then
      iv=1
      ip=2
    elseif(dir=='j') then
      iv=2
      ip=3
    elseif(dir=='k') then
      iv=3
      ip=4
    endif

    do i=-lmax,dim+lmax

      do m=1,5
       u(m,i) = velo(i,iv)
      enddo
      u(6,i) = 1._rtype

      v(1,i) = 1._rtype
      v(2,i) = velo(i,1)
      v(3,i) = velo(i,2)
      v(4,i) = velo(i,3)
      v(5,i) = enth(i)
      v(6,i) = pres(i)

    enddo

    do i=-lmax,dim

      do l=1,lmax
         rhom = dens(i)+dens(i+l)

         uv(1:5,l,i) = (u(1:5,i)+u(1:5,i+l))*(v(1:5,i)+v(1:5,i+l))*rhom

         uv(6,l,i) = (u(6,i)+u(6,i+l))*(v(6,i)+v(6,i+l))
      enddo

    enddo

    do i=-1,dim

      ft = 0._rtype
      do l=1,lmax

        uvs = 0._rtype
        do n=0,l-1
          uvs = uvs + uv(:,l,i-n)
        enddo

        ft = ft + dcoe(l,lmax)*uvs
      enddo

      fh(1:5) = 0.25_rtype*ft(1:5)

      fh(ip) = fh(ip) + 0.5_rtype*ft(6)

      flux(i,1:5) = fh(:)
    enddo

    return

  end function flux_sp

  function flux_div(q,velo,pres,dim,dir) result(flux)

    use comvardef, only : hm

    integer,intent(in) :: dim
    character(len=1),intent(in) :: dir
    real(rtype),intent(in) :: q(-hm:dim+hm,1:5)
    real(rtype),intent(in) :: velo(-hm:dim+hm,1:3),pres(-hm:dim+hm)
    real(rtype) :: flux(-1:dim,1:5)

    ! local data
    integer :: i
    real(rtype) :: pm
    real(rtype) :: qh(-hm:dim+hm,1:5)
    integer :: ip,iv

    if(dir=='i') then
      iv=1
      ip=2
    elseif(dir=='j') then
      iv=2
      ip=3
    elseif(dir=='k') then
      iv=3
      ip=4
    endif

    do i=-hm,dim+hm
      qh(i,1) = q(i,1)*velo(i,iv)
      qh(i,2) = q(i,2)*velo(i,iv)
      qh(i,3) = q(i,3)*velo(i,iv)
      qh(i,4) = q(i,4)*velo(i,iv)
      qh(i,5) = (q(i,5)+pres(i))*velo(i,iv)

      qh(i,ip) = qh(i,ip) + pres(i)
    enddo


    ! forth order central flux reconstruction
    do i=-1,dim
      
      flux(i,1) = num7d12*(qh(i,1)+qh(i+1,1))-num1d12*(qh(i-1,1)+qh(i+2,1))
      flux(i,2) = num7d12*(qh(i,2)+qh(i+1,2))-num1d12*(qh(i-1,2)+qh(i+2,2))
      flux(i,3) = num7d12*(qh(i,3)+qh(i+1,3))-num1d12*(qh(i-1,3)+qh(i+2,3))
      flux(i,4) = num7d12*(qh(i,4)+qh(i+1,4))-num1d12*(qh(i-1,4)+qh(i+2,4))
      flux(i,5) = num7d12*(qh(i,5)+qh(i+1,5))-num1d12*(qh(i-1,5)+qh(i+2,5))

    enddo


    ! compact flux reconstruction
    ! call fluxrecon4cc(qh(:,1),dim,flux(:,1),dir)
    ! call fluxrecon4cc(qh(:,2),dim,flux(:,2),dir)
    ! call fluxrecon4cc(qh(:,3),dim,flux(:,3),dir)
    ! call fluxrecon4cc(qh(:,4),dim,flux(:,4),dir)
    ! call fluxrecon4cc(qh(:,5),dim,flux(:,5),dir)

    return

  end function flux_div

  subroutine fluxrecon4cc(vin,dim,vout,dir)

    use comvardef, only : hm
    use numerics, only: cci,ccj,cck,qtds_solver
    use constdef

    integer,intent(in) :: dim
    real(rtype),intent(in) :: vin(-hm:dim+hm)
    character(len=1),intent(in) :: dir
    real(rtype) :: vout(-1:dim),b(0:dim-1)

    ! local data
    integer :: i

    do i=0,dim-1
      b(i)=0.75_rtype* (vin(i)+vin(i+1))
    end do

    if(dir=='i') then
      call qtds_solver(b,vout(0:dim-1),cci,dim)
    elseif(dir=='j') then
      call qtds_solver(b,vout(0:dim-1),ccj,dim)
    elseif(dir=='k') then
      call qtds_solver(b,vout(0:dim-1),cck,dim)
    endif

    vout(-1)=vout(dim-1)
    vout(dim)=vout(0)

  end subroutine fluxrecon4cc

end module numericalflux