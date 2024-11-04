program main
  !
  integer :: im=256,jm=256,km=256,hm=3
  real(8),parameter :: pi=4.d0*atan(1.0_8)
  !
  real(8),allocatable,dimension(:,:,:,:) :: x,df,df_ref
  real(8),allocatable,dimension(:,:,:) :: f
  !
  integer :: i,j,k
  real(8) :: dx,dy,dz,error(3)
  real(8) :: time_beg,time_end
  !
  allocate(x(0:im,0:jm,0:km,1:3),f(-hm:im+hm,-hm:jm+hm,-hm:km+hm),df(0:im,0:jm,0:km,1:3),df_ref(0:im,0:jm,0:km,1:3))
  !
  call cpu_time(time_beg)
  ! initial value
  do k=0,km
  do j=0,jm
  do i=0,im
    !
    x(i,j,k,1)  =2.d0*pi/dble(im)*dble(i)
    x(i,j,k,2)  =2.d0*pi/dble(jm)*dble(j)
    x(i,j,k,3)  =2.d0*pi/dble(km)*dble(k)
    !
    f(i,j,k)    = sin(x(i,j,k,1))*cos(x(i,j,k,2))*cos(x(i,j,k,3))
    !
    df_ref(i,j,k,1) =  cos(x(i,j,k,1))*cos(x(i,j,k,2))*cos(x(i,j,k,3))
    df_ref(i,j,k,2) = -sin(x(i,j,k,1))*sin(x(i,j,k,2))*cos(x(i,j,k,3))
    df_ref(i,j,k,3) = -sin(x(i,j,k,1))*cos(x(i,j,k,2))*sin(x(i,j,k,3))
  enddo
  enddo
  enddo
  call cpu_time(time_end)
  print*,' ** CPU time for initialisation :',time_end-time_beg
  !
  dx=x(1,1,1,1)-x(0,0,0,1)
  dy=x(1,1,1,2)-x(0,0,0,2)
  dz=x(1,1,1,3)-x(0,0,0,3)
  !
  call cpu_time(time_beg)
  ! boundary condition
  do k=0,km
  do j=0,jm
     f(-hm:-1,j,k) = f(im-hm:im-1,j,k)
     f(im+1:im+hm,j,k) = f(1:hm,j,k)
  enddo
  enddo
  do k=0,km
  do i=0,im
     f(i,-hm:-1,k) = f(i,jm-hm:jm-1,k)
     f(i,jm+1:jm+hm,k) = f(i,1:hm,k)
  enddo
  enddo
  do j=0,jm
  do i=0,im
     f(i,j,-hm:-1) = f(i,j,km-hm:km-1)
     f(i,j,km+1:km+hm) = f(i,j,1:hm)
  enddo
  enddo
  call cpu_time(time_end)
  print*,' ** CPU time for boundary condition :',time_end-time_beg
  !
  call cpu_time(time_beg)
  ! gradient calculation with 6th-order scheme
  !
  do k=0,km
  do j=0,jm
  do i=0,im
    df(i,j,k,1)  =0.75d0 *(f(i+1,j,k)-f(i-1,j,k))- &
                  0.15d0 *(f(i+2,j,k)-f(i-2,j,k))+ &
                  1.66666666666667d-2*(f(i+3,j,k)-f(i-3,j,k))
    df(i,j,k,2)  =0.75d0 *(f(i,j+1,k)-f(i,j-1,k))- &
                  0.15d0 *(f(i,j+2,k)-f(i,j-2,k))+ &
                  1.66666666666667d-2*(f(i,j+3,k)-f(i,j-3,k))
    df(i,j,k,3)  =0.75d0 *(f(i,j,k+1)-f(i,j,k-1))- &
                  0.15d0 *(f(i,j,k+2)-f(i,j,k-2))+ &
                  1.66666666666667d-2*(f(i,j,k+3)-f(i,j,k-3))
  enddo
  enddo
  enddo
  df(:,:,:,1)=df(:,:,:,1)/dx
  df(:,:,:,2)=df(:,:,:,2)/dy
  df(:,:,:,3)=df(:,:,:,3)/dz
  !
  call cpu_time(time_end)
  print*,' ** CPU time for gradient calculation :',time_end-time_beg
  !
  call cpu_time(time_beg)
  ! validation
  error=0.d0
  do k=1,km
  do j=1,jm
  do i=1,im
    !
    error(1) = error(1) + (df(i,j,k,1)-df_ref(i,j,k,1))**2
    error(2) = error(2) + (df(i,j,k,2)-df_ref(i,j,k,2))**2
    error(3) = error(3) + (df(i,j,k,3)-df_ref(i,j,k,3))**2
    !
    ! if(i==1 .and. j==1) print*,'x-y-z:df3',i,j,k,x(i,j,k,:),df(i,j,k,3)
    ! if(i==1 .and. k==1) print*,'x-y-z:df2',x(i,j,k,:),df(i,j,k,2)
    !
  enddo
  enddo
  enddo
  error=error/dble(im*jm*km)
  !
  print*,' ** errors:',error
  !
  call cpu_time(time_end)
  print*,' ** CPU time for validation :',time_end-time_beg
  !
  print*,' ** job done.'
  !
end program
!