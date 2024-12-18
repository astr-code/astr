module bc
   
  implicit none

  contains
  
  subroutine boundarycondition
    
    call bchomo

  end subroutine boundarycondition

  subroutine bchomo
    !
    use comvardef, only: im,jm,km,hm,rho,prs,tmp,vel,q
    use fluidynamcs,  only: var2q
    !
    integer :: i,j,k
    !
    ! b.c. at the i direction
    do k=0,km
    do j=0,jm
      !
      do i=-hm,-1
        rho(i,j,k)  =rho(im+i,j,k)
        vel(i,j,k,:)=vel(im+i,j,k,:)
        tmp(i,j,k)  =tmp(im+i,j,k)
        prs(i,j,k)  =prs(im+i,j,k)
        !
        q(i,j,k,:)  =var2q(density=rho(i,j,k),   &
                          velocity=vel(i,j,k,:), &
                          pressure=prs(i,j,k))
      enddo
      !
      do i=im+1,im+hm
        rho(i,j,k)  =rho(i-im,j,k)
        vel(i,j,k,:)=vel(i-im,j,k,:)
        tmp(i,j,k)  =tmp(i-im,j,k)
        prs(i,j,k)  =prs(i-im,j,k)
        !
        q(i,j,k,:)  =var2q(density=rho(i,j,k),   &
                          velocity=vel(i,j,k,:), &
                          pressure=prs(i,j,k))
      enddo
      !
    enddo
    enddo
    ! end of applying b.c. along i direction
    !
    ! b.c. at the j direction
    do k=0,km
    do i=0,im
      !
      do j=-hm,-1
        rho(i,j,k)  =rho(i,jm+j,k)
        vel(i,j,k,:)=vel(i,jm+j,k,:)
        tmp(i,j,k)  =tmp(i,jm+j,k)
        prs(i,j,k)  =prs(i,jm+j,k)
        !
        q(i,j,k,:)  =var2q(density=rho(i,j,k),   &
                          velocity=vel(i,j,k,:), &
                          pressure=prs(i,j,k))
      enddo
      !
      do j=jm+1,jm+hm
        rho(i,j,k)  =rho(i,j-jm,k)
        vel(i,j,k,:)=vel(i,j-jm,k,:)
        tmp(i,j,k)  =tmp(i,j-jm,k)
        prs(i,j,k)  =prs(i,j-jm,k)
        !
        q(i,j,k,:)  =var2q(density=rho(i,j,k),   &
                          velocity=vel(i,j,k,:), &
                          pressure=prs(i,j,k))
      enddo
      !
    enddo
    enddo
    ! end of applying b.c. along j direction
    !
    ! b.c. at the k direction
    do j=0,jm
    do i=0,im
      !
      do k=-hm,-1
        rho(i,j,k)  =rho(i,j,km+k)
        vel(i,j,k,:)=vel(i,j,km+k,:)
        tmp(i,j,k)  =tmp(i,j,km+k)
        prs(i,j,k)  =prs(i,j,km+k)
        !
        q(i,j,k,:)  =var2q(density=rho(i,j,k),   &
                          velocity=vel(i,j,k,:), &
                          pressure=prs(i,j,k))
      enddo
      !
      do k=km+1,km+hm
        rho(i,j,k)  =rho(i,j,k-km)
        vel(i,j,k,:)=vel(i,j,k-km,:)
        tmp(i,j,k)  =tmp(i,j,k-km)
        prs(i,j,k)  =prs(i,j,k-km)
        !
        q(i,j,k,:)  =var2q(density=rho(i,j,k),   &
                          velocity=vel(i,j,k,:), &
                          pressure=prs(i,j,k))
      enddo
      !
    enddo
    enddo
    ! end of applying b.c. along k direction
    !
  end subroutine bchomo
  !
  subroutine bchomovec(var)
    !
    use comvardef, only: im,jm,km,hm
    !
    real(8),intent(inout) :: var(-hm:,-hm:,-hm:,1:)
    !
    integer :: i,j,k,n,nd4
    !
    nd4=size(var,4)
    !
    ! b.c. at the i direction
    do k=0,km
    do j=0,jm
      !
      do i=-hm,-1
        var(i,j,k,:)  =var(im+i,j,k,:)
      enddo
      !
      do i=im+1,im+hm
        var(i,j,k,:)=var(i-im,j,k,:)
      enddo
      !
    enddo
    enddo
    ! end of applying b.c. along i direction
    !
    ! b.c. at the j direction
    do k=0,km
    do i=0,im
      !
      do j=-hm,-1
        var(i,j,k,:)=var(i,jm+j,k,:)
      enddo
      !
      do j=jm+1,jm+hm
        var(i,j,k,:)=var(i,j-jm,k,:)
      enddo
      !
    enddo
    enddo
    ! end of applying b.c. along j direction
    !
    ! b.c. at the k direction
    do j=0,jm
    do i=0,im
      !
      do k=-hm,-1
        var(i,j,k,:)=var(i,j,km+k,:)
      enddo
      !
      do k=km+1,km+hm
        var(i,j,k,:)=var(i,j,k-km,:)
      enddo
      !
    enddo
    enddo
    ! end of applying b.c. along k direction
    !
  end subroutine bchomovec

end module bc