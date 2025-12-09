module pastr_gradients

    use iso_fortran_env, only: wp => real64
    use pastr_commvar, only: im,jm,km,lihomo,ljhomo,lkhomo,nhalo
    use pastr_constdef

    implicit none

contains
    
    function grad_3d(f,x,y,z) result(df)
      ! 
      ! arguments
      real(wp),intent(in) :: f(0:,0:,0:)
      real(wp),intent(in),optional :: x(0:,0:,0:),y(0:,0:,0:),&
                                     z(0:,0:,0:)
      real(wp) :: df(1:3,0:size(f,1)-1,0:size(f,2)-1,0:size(f,3)-1)
      !
  
      ! local data
      integer :: i,j,k
      real(wp),allocatable :: vtemp(:)
      real(wp),allocatable :: ddi(:,:,:,:),ddj(:,:,:,:),ddk(:,:,:,:)
      save ddi,ddj,ddk
      !$ SAVE vtemp
      !$OMP THREADPRIVATE(vtemp)
      !
      im=size(f,1)-1
      jm=size(f,2)-1
      km=size(f,3)-1

      if(.not. allocated(ddi)) then
        allocate(ddi(1:3,0:im,0:jm,0:km),ddj(1:3,0:im,0:jm,0:km),        &
                 ddk(1:3,0:im,0:jm,0:km) )
        !
        if(present(x) .and. present(y) .and. present(z)) then
          call gridjacobian_3d(x,y,z,ddi,ddj,ddk)
        else
          stop ' !! x,y,z needed in function grad_3d !!'
        endif
        !
      else
        ! print*,' ** Use the saved grid jacobian matrix'
      endif
      !
      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k)
      !
      allocate(vtemp(0:im))
      !$OMP DO
      do k=0,km
      do j=0,jm
        vtemp=dfdi(f(:,j,k))
        df(1,:,j,k)=vtemp(:)*ddi(1,:,j,k)
        df(2,:,j,k)=vtemp(:)*ddi(2,:,j,k)
        df(3,:,j,k)=vtemp(:)*ddi(3,:,j,k)
      end do
      end do
      !$OMP END DO
      deallocate(vtemp)
      !
      allocate(vtemp(0:jm))
      !$OMP DO
      do k=0,km
      do i=0,im
        vtemp=dfdj(f(i,:,k))
        df(1,i,:,k)=df(1,i,:,k)+vtemp(:)*ddj(1,i,:,k)
        df(2,i,:,k)=df(2,i,:,k)+vtemp(:)*ddj(2,i,:,k)
        df(3,i,:,k)=df(3,i,:,k)+vtemp(:)*ddj(3,i,:,k)
      end do
      end do
      !$OMP END DO
      deallocate(vtemp)
      !
      allocate(vtemp(0:km))
      !$OMP DO
      do j=0,jm
      do i=0,im
        vtemp=dfdk(f(i,j,:))
        df(1,i,j,:)=df(1,i,j,:)+vtemp(:)*ddk(1,i,j,:)
        df(2,i,j,:)=df(2,i,j,:)+vtemp(:)*ddk(2,i,j,:)
        df(3,i,j,:)=df(3,i,j,:)+vtemp(:)*ddk(3,i,j,:)
      end do
      end do
      !$OMP END DO
      deallocate(vtemp)
      !
      !$OMP END PARALLEL
      !
      print*,' ** 3D gradient field calculated'
      !
    end function grad_3d

    subroutine gridjacobian_3d(x,y,z,ddi,ddj,ddk)
      !
      ! arguments
      real(wp),intent(in) :: x(0:im,0:jm,0:km),y(0:im,0:jm,0:km),         &
                             z(0:im,0:jm,0:km)
      real(wp),intent(out) :: ddi(1:3,0:im,0:jm,0:km),                    &
                              ddj(1:3,0:im,0:jm,0:km),                    &
                              ddk(1:3,0:im,0:jm,0:km)
      !
      ! local data
      integer :: i,j,k
      real(wp) :: var1
      real(wp),allocatable :: dx(:,:,:,:),dy(:,:,:,:),dz(:,:,:,:)
      !
      allocate(dx(1:3,0:im,0:jm,0:km),dy(1:3,0:im,0:jm,0:km),            &
               dz(1:3,0:im,0:jm,0:km))
      !
      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k,var1)
      !
      !$OMP DO
      do k=0,km
      do j=0,jm
        dx(1,:,j,k)=dfdi(x(:,j,k))
        dy(1,:,j,k)=dfdi(y(:,j,k))
        dz(1,:,j,k)=dfdi(z(:,j,k))
      end do
      end do
      !$OMP END DO
      !
      !$OMP DO
      do k=0,km
      do i=0,im
        dx(2,i,:,k)=dfdj(x(i,:,k))
        dy(2,i,:,k)=dfdj(y(i,:,k))
        dz(2,i,:,k)=dfdj(z(i,:,k))
      end do
      end do
      !$OMP END DO
      !
      !$OMP DO
      do j=0,jm
      do i=0,im
        dx(3,i,j,:)=dfdk(x(i,j,:),geom=.true.)
        dy(3,i,j,:)=dfdk(y(i,j,:),geom=.true.)
        dz(3,i,j,:)=dfdk(z(i,j,:),geom=.true.)
      end do
      end do
      !$OMP END DO

      !$OMP DO
      do k=0,km
      do j=0,jm
      do i=0,im
        var1= dx(1,i,j,k)*dy(2,i,j,k)*dz(3,i,j,k)                        &
             +dx(2,i,j,k)*dy(3,i,j,k)*dz(1,i,j,k)                        &
             +dx(3,i,j,k)*dy(1,i,j,k)*dz(2,i,j,k)                        &
             -dx(3,i,j,k)*dy(2,i,j,k)*dz(1,i,j,k)                        &
             -dx(2,i,j,k)*dy(1,i,j,k)*dz(3,i,j,k)                        &
             -dx(1,i,j,k)*dy(3,i,j,k)*dz(2,i,j,k)
        !
        ddi(1,i,j,k)=dy(2,i,j,k)*dz(3,i,j,k)-dy(3,i,j,k)*dz(2,i,j,k)
        ddi(2,i,j,k)=dx(3,i,j,k)*dz(2,i,j,k)-dx(2,i,j,k)*dz(3,i,j,k)
        ddi(3,i,j,k)=dx(2,i,j,k)*dy(3,i,j,k)-dx(3,i,j,k)*dy(2,i,j,k)
        !
        ddj(1,i,j,k)=dy(3,i,j,k)*dz(1,i,j,k)-dy(1,i,j,k)*dz(3,i,j,k)
        ddj(2,i,j,k)=dx(1,i,j,k)*dz(3,i,j,k)-dx(3,i,j,k)*dz(1,i,j,k)
        ddj(3,i,j,k)=dx(3,i,j,k)*dy(1,i,j,k)-dx(1,i,j,k)*dy(3,i,j,k)
        !
        ddk(1,i,j,k)=dy(1,i,j,k)*dz(2,i,j,k)-dy(2,i,j,k)*dz(1,i,j,k)
        ddk(2,i,j,k)=dx(2,i,j,k)*dz(1,i,j,k)-dx(1,i,j,k)*dz(2,i,j,k)
        ddk(3,i,j,k)=dx(1,i,j,k)*dy(2,i,j,k)-dx(2,i,j,k)*dy(1,i,j,k)
        !
        ddi(1,i,j,k)=ddi(1,i,j,k)/var1
        ddi(2,i,j,k)=ddi(2,i,j,k)/var1
        ddi(3,i,j,k)=ddi(3,i,j,k)/var1
        ddj(1,i,j,k)=ddj(1,i,j,k)/var1
        ddj(2,i,j,k)=ddj(2,i,j,k)/var1
        ddj(3,i,j,k)=ddj(3,i,j,k)/var1
        ddk(1,i,j,k)=ddk(1,i,j,k)/var1
        ddk(2,i,j,k)=ddk(2,i,j,k)/var1
        ddk(3,i,j,k)=ddk(3,i,j,k)/var1
      end do
      end do
      end do
      !$OMP END DO
      !
      !$OMP END PARALLEL

      deallocate(dx,dy,dz)
      !
      print*, ' ** Grid Jacobian matrix is calculated'
      !
    end subroutine gridjacobian_3d

    function grad_xy(f,x,y) result(df)

      ! arguments
      real(wp) :: df(1:2,0:im,0:jm)
      real(wp),intent(in) :: f(0:im,0:jm)
      real(wp),intent(in),optional :: x(0:im,0:jm),y(0:im,0:jm)
      !
      ! local data
      integer :: i,j
      real(wp),allocatable :: vtemp(:)
      real(wp),allocatable :: ddi(:,:,:),ddj(:,:,:)
      save ddi,ddj

      if(.not. allocated(ddi)) then

        allocate(ddi(1:2,0:im,0:jm),ddj(1:2,0:im,0:jm) )

        if(present(x) .and. present(y)) then
          call gridjacobian_xy(x,y,ddi,ddj)
        else
          stop ' !! x,y needed in function grad_xy !!'
        endif
        !
      else
        ! print*,' ** Use the saved grid jacobian matrix'
      endif
      !
      allocate(vtemp(0:im))
      do j=0,jm
        vtemp=dfdi(f(:,j))
        df(1,:,j)=vtemp(:)*ddi(1,:,j)
        df(2,:,j)=vtemp(:)*ddi(2,:,j)
      end do
      deallocate(vtemp)
      !
      allocate(vtemp(0:jm))
      do i=0,im
        vtemp=dfdj(f(i,:))
        df(1,i,:)=df(1,i,:)+vtemp(:)*ddj(1,i,:)
        df(2,i,:)=df(2,i,:)+vtemp(:)*ddj(2,i,:)
      end do
      deallocate(vtemp)
      !
      ! print*,' ** Spatial gradient at x-y plane calculated'
      !
    end function grad_xy

    subroutine gridjacobian_xy(x,y,ddi,ddj,wallnormal)
      !
      !
      ! arguments
      real(wp),intent(in) :: x(0:im,0:jm),y(0:im,0:jm)
      real(wp),intent(out),optional :: ddi(2,0:im,0:jm),ddj(2,0:im,0:jm)
      real(wp),intent(out),optional :: wallnormal(2,0:im)
      !
      ! local data
      integer :: i,j
      real(wp) :: var1
      real(wp),allocatable :: dx(:,:,:),dy(:,:,:)
      real(wp),allocatable :: di(:,:,:),dj(:,:,:)
      !
      allocate(dx(2,0:im,0:jm),dy(2,0:im,0:jm),                          &
               di(2,0:im,0:jm),dj(2,0:im,0:jm))
      !
      do j=0,jm
        dx(1,:,j)=dfdi(x(:,j))
        dy(1,:,j)=dfdi(y(:,j))
      end do
      !
      do i=0,im
        dx(2,i,:)=dfdj(x(i,:))
        dy(2,i,:)=dfdj(y(i,:))
      end do
      !
      do j=0,jm
      do i=0,im
        !
        var1=dx(1,i,j)*dy(2,i,j)-dx(2,i,j)*dy(1,i,j)
        !
        var1=1._wp/var1
        !
        di(1,i,j)= var1*dy(2,i,j)
        di(2,i,j)=-var1*dx(2,i,j)
        !
        dj(1,i,j)=-var1*dy(1,i,j)
        dj(2,i,j)= var1*dx(1,i,j)
        !
      end do
      end do
      !
      if(present(ddi) .and. present(ddj)) then
        ddi=di
        ddj=dj
      endif
      !
      if(present(wallnormal)) then
        do i=0,im
          !
          var1=sqrt(dj(1,i,0)**2+dj(2,i,0)**2)
          !
          wallnormal(1,i)=dj(1,i,0)/var1
          wallnormal(2,i)=dj(2,i,0)/var1
          !
        end do
      endif
      !
      deallocate(dx,dy)
      !
      print*, ' ** Grid Jacobian matrix in a x-y plane is calculated'
      !
    end subroutine gridjacobian_xy

    function dfdi(var)

      real(wp),allocatable :: dfdi(:)
      real(wp),intent(in) :: var(0:im)
      !
      allocate(dfdi(0:im))
      !
      dfdi=diff6ec(var,im,lihomo)
      !
    end function dfdi
    !
    function dfdj(var)

      real(wp),allocatable :: dfdj(:)
      real(wp),intent(in) :: var(0:jm)
      !
      dfdj =diff6ec(var,jm,ljhomo)

    end function dfdj
    !
    function dfdk(var,geom)

      real(wp) ,allocatable :: dfdk(:)
      real(wp),intent(in) :: var(0:km)
      logical,intent(in),optional :: geom

      logical :: lgeom

      if(present(geom)) then
        lgeom=geom
      else
        lgeom=.false.
      endif
      
      if(lgeom) then
        dfdk =diff6ec_geom(var,km,lkhomo)
      else
        dfdk =diff6ec(var,km,lkhomo)
      endif

    end function dfdk

    !+-------------------------------------------------------------------+
    !| This function returns the approximation of derivative of a given  |
    !| array, using 6th-order finite-difference scheme, including the    |
    !| boundary close.                                                   |
    !+-------------------------------------------------------------------+
    !| CHANGE RECORD                                                     |
    !| -------------                                                     |
    !| 31-05-2025  | Moved to here by J. Fang @ Liverpool.               |
    !+-------------------------------------------------------------------+
    pure function diff6ec(vin,dim,homo) result(vout)
      
      integer,intent(in) :: dim
      logical,intent(in) :: homo
      real(8),intent(in) :: vin(0:dim)
      real(8) :: vout(0:dim)
      
      ! local data
      integer :: i

      if(homo) then

        vout(0)=0.75d0 *(vin(1)-vin(dim-1))-           &
                0.15d0 *(vin(2)-vin(dim-2))+           &
                num1d60*(vin(3)-vin(dim-3))       
        vout(1)=0.75d0 *(vin(2)-vin(0))    -           &
                0.15d0 *(vin(3)-vin(dim-1))+           &
                num1d60*(vin(4)-vin(dim-2))       
        vout(2)=0.75d0 *(vin(3)-vin(1))    -           &
                0.15d0 *(vin(4)-vin(0))    +           &
                num1d60*(vin(5)-vin(dim-1))

        vout(dim-2) =0.75d0 *(vin(dim-1)-vin(dim-3))-  &
                     0.15d0 *(vin(dim)  -vin(dim-4))+  &
                     num1d60*(vin(1)    -vin(dim-5))
        vout(dim-1)=0.75d0 *(vin(dim)-vin(dim-2))-     &
                    0.15d0 *(vin(1)  -vin(dim-3))+     &
                    num1d60*(vin(2)  -vin(dim-4))
        vout(dim)  =0.75d0 *(vin(1)-vin(dim-1))-       &
                    0.15d0 *(vin(2)-vin(dim-2))+       &
                    num1d60*(vin(3)-vin(dim-3))
      else

        vout(0)=-0.5d0*vin(2)+2.d0*vin(1)-1.5d0*vin(0)
        vout(1)=0.5d0*(vin(2)-vin(0))
        vout(2)=num2d3*(vin(3)-vin(1))-num1d12*(vin(4)-vin(0))

        vout(dim-2) =num2d3*(vin(dim-1)-vin(dim-3))-      &
                     num1d12*(vin(dim)  -vin(dim-4))
        vout(dim-1)=0.5d0*(vin(dim)-vin(dim-2))
        vout(dim)  =0.5d0*vin(dim-2)-2.d0*vin(dim-1)+1.5d0*vin(dim)

      endif

      do i=3,dim-3
        vout(i)  =0.75d0 *(vin(i+1)-vin(i-1))-            &
                  0.15d0 *(vin(i+2)-vin(i-2))+            &
                  num1d60*(vin(i+3)-vin(i-3))
      enddo
    end function diff6ec
    !+-------------------------------------------------------------------+
    !| The end of the function diff6ec.                                  |
    !+-------------------------------------------------------------------+
    function diff6ec_geom(vin,dim,homo) result(vout)
      
      integer,intent(in) :: dim
      logical,intent(in) :: homo
      real(8),intent(in) :: vin(0:dim)
      real(8) :: vout(0:dim)
      
      ! local data
      integer :: i

      if(homo) then

        vout(0)=0.75d0 *((vin(1)-vin(0))-(vin(dim-1)-vin(dim)))-           &
                0.15d0 *((vin(2)-vin(0))-(vin(dim-2)-vin(dim)))+           &
                num1d60*((vin(3)-vin(0))-(vin(dim-3)-vin(dim)))     
        vout(1)=0.75d0 *((vin(2)-vin(0))-(vin(0)-vin(0)))    -           &
                0.15d0 *((vin(3)-vin(0))-(vin(dim-1)-vin(dim)))+           &
                num1d60*((vin(4)-vin(0))-(vin(dim-2)-vin(dim)))       
        vout(2)=0.75d0 *((vin(3)-vin(0))-(vin(1)-vin(0)))    -           &
                0.15d0 *((vin(4)-vin(0))-(vin(0)-vin(0)))    +           &
                num1d60*((vin(5)-vin(0))-(vin(dim-1)-vin(dim)))

        vout(dim-2) =0.75d0 *(vin(dim-1)-vin(dim-3))-  &
                     0.15d0 *(vin(dim)  -vin(dim-4))+  &
                     num1d60*(vin(1)-vin(0)    -(vin(dim-5)-vin(dim)))
        vout(dim-1)=0.75d0 *(vin(dim)-vin(dim-2))-     &
                    0.15d0 *(vin(1)-vin(0)  +vin(dim)-vin(dim-3))+     &
                    num1d60*(vin(2)-vin(0)  +vin(dim)-vin(dim-4))
        vout(dim)  =0.75d0 *(vin(1)-vin(0)  +vin(dim)-vin(dim-1))-       &
                    0.15d0 *(vin(2)-vin(0)  +vin(dim)-vin(dim-2))+       &
                    num1d60*(vin(3)-vin(0)  +vin(dim)-vin(dim-3))
      else

        vout(0)=-0.5d0*vin(2)+2.d0*vin(1)-1.5d0*vin(0)
        vout(1)=0.5d0*(vin(2)-vin(0))
        vout(2)=num2d3*(vin(3)-vin(1))-num1d12*(vin(4)-vin(0))

        vout(dim-2) =num2d3*(vin(dim-1)-vin(dim-3))-      &
                     num1d12*(vin(dim)  -vin(dim-4))
        vout(dim-1)=0.5d0*(vin(dim)-vin(dim-2))
        vout(dim)  =0.5d0*vin(dim-2)-2.d0*vin(dim-1)+1.5d0*vin(dim)

      endif

      do i=3,dim-3
        vout(i)  =0.75d0 *(vin(i+1)-vin(i-1))-            &
                  0.15d0 *(vin(i+2)-vin(i-2))+            &
                  num1d60*(vin(i+3)-vin(i-3))
      enddo

    end function diff6ec_geom
    !+-------------------------------------------------------------------+
    !| The end of the function diff6ec.                                  |
    !+-------------------------------------------------------------------+

end module pastr_gradients