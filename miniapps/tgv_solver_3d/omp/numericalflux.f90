module numericalflux
  
  use constdef
  use comvardef,only: len_sten_max

  implicit none

  real(rtype),allocatable :: cci(:,:),ccj(:,:),cck(:,:)

  real(rtype) :: dcoe(5,5)
  
  contains

  subroutine flux_solver_init

    use comvardef,only: im,jm,km
    use numerics, only: qtds_solver_init

    character(len=6) :: flux_scheme

    select case (len_sten_max)
    case (1)
       flux_scheme='cflux4'
    case (2)
       flux_scheme='cflux6'
    case (3)
       flux_scheme='cflux8'
    case default
       stop 'flux_scheme not defined'
    end select

    call qtds_solver_init(cci,flux_scheme,im)
    call qtds_solver_init(ccj,flux_scheme,jm)
    call qtds_solver_init(cck,flux_scheme,km)

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

    dcoe(1,5) =  5._rtype/6._rtype
    dcoe(2,5) = -5._rtype/21._rtype
    dcoe(3,5) =  5._rtype/84._rtype
    dcoe(4,5) = -5._rtype/504._rtype 
    dcoe(5,5) =  1._rtype/1260._rtype 

  end subroutine flux_solver_init

  ! function flux_sp(dens,velo,pres,enth,dim,dir) result(flux)

  !   use comvardef, only : hm

  !   integer,intent(in) :: dim
  !   character(len=1),intent(in) :: dir
  !   real(rtype),intent(in) :: dens(-hm:dim+hm),velo(-hm:dim+hm,1:3), &
  !                             pres(-hm:dim+hm),enth(-hm:dim+hm)
  !   real(rtype) :: flux(-1:dim,1:5)

  !   ! local data
  !   integer :: i,l,m,n
  !   real(rtype) :: uv(1:6,1:len_sten_max,-hm:dim)
  !   real(rtype) :: fh(5),ft(6),uvs(6)
  !   real(rtype) :: u(1:6,-hm:dim+hm),v(1:6,-hm:dim+hm)
  !   real(rtype) :: rhom
  !   integer :: ip,iv

  !   if(dir=='i') then
  !     iv=1
  !     ip=2
  !   elseif(dir=='j') then
  !     iv=2
  !     ip=3
  !   elseif(dir=='k') then
  !     iv=3
  !     ip=4
  !   endif

  !   do i=-len_sten_max,dim+len_sten_max

  !     do m=1,5
  !      u(m,i) = velo(i,iv)
  !     enddo
  !     u(6,i) = 1._rtype

  !     v(1,i) = 1._rtype
  !     v(2,i) = velo(i,1)
  !     v(3,i) = velo(i,2)
  !     v(4,i) = velo(i,3)
  !     v(5,i) = enth(i)
  !     v(6,i) = pres(i)

  !   enddo

  !   do i=-len_sten_max,dim

  !     do l=1,len_sten_max
  !        rhom = dens(i)+dens(i+l)

  !        uv(1:5,l,i) = (u(1:5,i)+u(1:5,i+l))*(v(1:5,i)+v(1:5,i+l))*rhom

  !        uv(6,l,i) = (u(6,i)+u(6,i+l))*(v(6,i)+v(6,i+l))
  !     enddo

  !   enddo

  !   do i=-1,dim

  !     ft = 0._rtype
  !     do l=1,len_sten_max

  !       uvs = 0._rtype
  !       do n=0,l-1
  !         uvs = uvs + uv(:,l,i-n)
  !       enddo

  !       ft = ft + dcoe(l,len_sten_max)*uvs
  !     enddo

  !     fh(1:5) = 0.25_rtype*ft(1:5)

  !     fh(ip) = fh(ip) + 0.5_rtype*ft(6)

  !     flux(i,1:5) = fh(:)
  !   enddo

  !   return

  ! end function flux_sp

  function flux_stable(dens,velo,pres,enth,dim,dir) result(flux)

    use comvardef, only : hm

    integer,intent(in) :: dim
    character(len=1),intent(in) :: dir
    real(rtype),intent(in) :: dens(-hm:dim+hm),velo(-hm:dim+hm,1:3), &
                              pres(-hm:dim+hm),enth(-hm:dim+hm)
    real(rtype) :: flux(-1:dim,1:5)

    ! local data
    integer :: i,l,m,n
    real(rtype) :: uv(1:6,1:len_sten_max,-hm:dim+hm),uvf(-1:dim,1:6)
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

    do i=-hm,dim+hm

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

    do i=-hm,dim+hm

      do l=1,len_sten_max
         rhom = dens(i)+dens(i+l)

         uv(1:5,l,i) = (u(1:5,i)+u(1:5,i+l))*(v(1:5,i)+v(1:5,i+l))*rhom

         uv(6,l,i) = (u(6,i)+u(6,i+l))*(v(6,i)+v(6,i+l))
      enddo

    enddo

    ! compact flux reconstruction
    ! call flux_recon_subs(uv(1,:,:),dim,uvf(:,1),dir)
    ! call flux_recon_subs(uv(2,:,:),dim,uvf(:,2),dir)
    ! call flux_recon_subs(uv(3,:,:),dim,uvf(:,3),dir)
    ! call flux_recon_subs(uv(4,:,:),dim,uvf(:,4),dir)
    ! call flux_recon_subs(uv(5,:,:),dim,uvf(:,5),dir)
    ! call flux_recon_subs(uv(6,:,:),dim,uvf(:,6),dir)
    call flux_recon_subs_compact(uv(1,:,:),dim,uvf(:,1),dir)
    call flux_recon_subs_compact(uv(2,:,:),dim,uvf(:,2),dir)
    call flux_recon_subs_compact(uv(3,:,:),dim,uvf(:,3),dir)
    call flux_recon_subs_compact(uv(4,:,:),dim,uvf(:,4),dir)
    call flux_recon_subs_compact(uv(5,:,:),dim,uvf(:,5),dir)
    call flux_recon_subs_compact(uv(6,:,:),dim,uvf(:,6),dir)

    do i=-1,dim

      uvf(i,ip) = uvf(i,ip) + 2._rtype*uvf(i,6)

      flux(i,1:5) = 0.25_rtype*uvf(i,1:5)

    enddo

    ! if(dir=='j') stop

      
    ! do i=-1,dim

    !   ft = 0._rtype
    !   do l=1,len_sten_max

    !     uvs = 0._rtype
    !     do n=0,l-1
    !       uvs = uvs + uv(:,l,i-n)
    !     enddo

    !     ft = ft + dcoe(l,len_sten_max)*uvs
    !   enddo

    !   fh(1:5) = 0.25_rtype*ft(1:5)

    !   fh(ip) = fh(ip) + 0.5_rtype*ft(6)

    !   flux(i,1:5) = fh(:)
    ! enddo

    return

  end function flux_stable

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

  subroutine fluxrecon4ec(vin,dim,vout,dir)

    use comvardef, only : hm
    use numerics, only: qtds_solver
    use constdef

    integer,intent(in) :: dim
    real(rtype),intent(in) :: vin(-hm:dim+hm)
    character(len=1),intent(in) :: dir
    real(rtype) :: vout(-1:dim),b(0:dim-1)

    ! local data
    integer :: i

    do i=-1,dim
      vout(i)=num7d12*(vin(i)+vin(i+1))-num1d12*(vin(i-1)+vin(i+2))
    end do

    return

  end subroutine fluxrecon4ec

  subroutine flux_recon_subs(vin,dim,vout,dir)

    use comvardef, only : hm
    use numerics, only: qtds_solver
    use constdef

    integer,intent(in) :: dim
    real(rtype),intent(in) :: vin(1:len_sten_max,-hm:dim+hm)
    character(len=1),intent(in) :: dir
    real(rtype) :: vout(-1:dim),b(0:dim-1),var1,var2

    ! local data
    integer :: i,l,n

    do i=-1,dim

      var1=0._rtype

      do l=1,len_sten_max

        var2 = 0._rtype
        do n=0,l-1
          var2 = var2 + vin(l,i-n)
        enddo

        var1=var1+dcoe(l,len_sten_max)*var2

      enddo

      vout(i)=var1

    end do

    return

  end subroutine flux_recon_subs

  subroutine flux_recon_subs_compact(vin,dim,vout,dir)

    use comvardef, only : hm
    use numerics, only: qtds_solver
    use constdef

    integer,intent(in) :: dim
    real(rtype),intent(in) :: vin(1:len_sten_max,-hm:dim+hm)
    character(len=1),intent(in) :: dir
    real(rtype) :: vout(-1:dim),b(0:dim-1)

    ! local data
    integer :: i

    if(len_sten_max==1) then

      do i=0,dim-1
        b(i)=0.75_rtype* vin(1,i)
      end do

    elseif(len_sten_max==2) then

      do i=0,dim-1
        b(i)=7._rtype/9._rtype*vin(1,i)+1._rtype/36._rtype*(vin(2,i)+vin(2,i-1))
      end do

    elseif(len_sten_max==3) then

      do i=0,dim-1
        b(i)=0.78125_rtype*vin(1,i)+ &
             0.05_rtype*(vin(2,i)+vin(2,i-1))- &
             1._rtype/480._rtype*(vin(3,i)+vin(3,i-1)+vin(3,i-2))
      end do

    endif

    if(dir=='i') then
      call qtds_solver(b,vout(0:dim-1),cci,dim)
    elseif(dir=='j') then
      call qtds_solver(b,vout(0:dim-1),ccj,dim)
    elseif(dir=='k') then
      call qtds_solver(b,vout(0:dim-1),cck,dim)
    endif

    vout(-1)=vout(dim-1)
    vout(dim)=vout(0)

  end subroutine flux_recon_subs_compact

end module numericalflux