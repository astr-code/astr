!+---------------------------------------------------------------------+
!| This module contains some solver related subroutines.               |
!+---------------------------------------------------------------------+
!| CHANGE RECORD                                                       |
!| -------------                                                       |
!| 08-02-2021  | Created by J. Fang                                    |
!+---------------------------------------------------------------------+
module solver
  !
  use constdef
  use parallel,  only : mpirankname,mpistop,mpirank,lio,dataswap
  !
  implicit none
  !
  real(8) :: alfa_con(3),alfa_dif(3)
  real(8), allocatable, dimension(:,:) :: cci,ccj,cck,dci,dcj,dck
  !
  contains
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to calculate some constant parameters          |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 08-02-2021  | Created by J. Fang @ Warrington                     |
  !+-------------------------------------------------------------------+
  subroutine refcal
    !
    use commvar,   only : numq
    !
    numq=1
    !
  end subroutine refcal
  !+-------------------------------------------------------------------+
  !| The end of the subroutine refcal.                                 |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to initialise solver.                          |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 08-02-2021  | Created by J. Fang @ Warrington                     |
  !+-------------------------------------------------------------------+
  subroutine solvrinit
    !
    use commvar, only : im,jm,km,numq,npdci,npdcj,npdck,conschm,difschm
    use commfunc,  only : coeffcompac,ptds_ini
    !
    ! local data
    integer :: nscheme
    !
    ! convectional term
    if(conschm(4:4)=='c') then
      ! a compact scheme is used
      !
      read(conschm(1:3),*) nscheme
      !
      alfa_con=coeffcompac(nscheme)
      !
      call ptds_ini(cci,alfa_con,nscheme,im,npdci)
      call ptds_ini(ccj,alfa_con,nscheme,jm,npdcj)
      call ptds_ini(cck,alfa_con,nscheme,km,npdck)
      !
    endif
    !
    ! diffusional term
    if(difschm(4:4)=='c') then
      ! a compact scheme is used
      !
      read(difschm(1:3),*) nscheme
      !
      alfa_dif=coeffcompac(nscheme)
      !
      call ptds_ini(dci,alfa_dif,nscheme,im,npdci)
      call ptds_ini(dcj,alfa_dif,nscheme,jm,npdcj)
      call ptds_ini(dck,alfa_dif,nscheme,km,npdck)
      !
    endif
    !
  end subroutine solvrinit
  !+-------------------------------------------------------------------+
  !| The end of the subroutine refcal.                                 |
  !+-------------------------------------------------------------------+
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! This subtoutine is used to calculate geometrical transform matrix
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Writen by Fang Jian, 2008-09-22.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine geomcal
    !
    use commvar,   only : ia,ja,ka,im,jm,km,hm,npdci,npdcj,npdck,hm,   &
                          xmax,xmin,ymax,ymin,zmax,zmin
    use commarray, only : x,jacob,dxi
    use parallel,  only : gridsendrecv,ksize,psum,pmax,pmin
    use commfunc,  only : coeffcompac,ptds_ini,ddf
    use tecio
    !
    ! local data
    character(len=4) :: cscheme
    integer :: nscheme
    integer :: i,j,k,m
    real(8) :: alfa(3)
    real(8), allocatable, dimension(:,:) :: gci,gcj,gck
    real(8), allocatable :: dx(:,:,:,:,:)
    real(8),allocatable :: phi(:),can(:,:,:,:)
    real(8) :: can1av,can2av,can3av,can1var,can2var,can3var
    !
    allocate( dx(-hm:im+hm,-hm:jm+hm,-hm:km+hm,1:3,1:3) )
    !
    call gridsendrecv
    !
    cscheme='642c'
    !
    read(cscheme(1:3),*) nscheme
    !
    if(cscheme(4:4)=='c') then
      ! a compact scheme is used
      !
      alfa=coeffcompac(nscheme)
      !
      call ptds_ini(gci,alfa,nscheme,im,npdci)
      call ptds_ini(gcj,alfa,nscheme,jm,npdcj)
      if(ksize>1) call ptds_ini(gck,alfa,nscheme,km,npdck)
      !
    endif
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Calculating d<x,y,z>/d<i,j,k>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do k=0,km
    do j=0,jm
      do m=1,3
        dx(0:im,j,k,m,1)=ddf(x(:,j,k,m),cscheme,npdci,im,alfa,gci)
      enddo
    enddo
    enddo
    !
    do k=0,km
    do i=0,im
      !
      do m=1,3
        dx(i,0:jm,k,m,2)=ddf(x(i,:,k,m),cscheme,npdcj,jm,alfa,gcj)
      enddo
      !
    enddo
    enddo
    !
    do j=0,jm
    do i=0,im
      !
      do m=1,3
        dx(i,j,0:km,m,3)=ddf(x(i,j,:,m),cscheme,npdck,km,alfa,gck)
      enddo
      !
    enddo
    enddo
    !
    if(lio) print*,' ** dxyz/dijk calculated.'
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! end of cal d<x,y,z>/d<i,j,k>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    call dataswap(dx)
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Calculating geometrical 
    ! Jacobian
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    jacob(0:im,0:jm,0:km)=dx(0:im,0:jm,0:km,1,1)*                      &
                          dx(0:im,0:jm,0:km,2,2)*                      &
                          dx(0:im,0:jm,0:km,3,3)                       &
                        + dx(0:im,0:jm,0:km,1,2)*                      &
                          dx(0:im,0:jm,0:km,2,3)*                      &
                          dx(0:im,0:jm,0:km,3,1)                       &
                        + dx(0:im,0:jm,0:km,1,3)*                      &
                          dx(0:im,0:jm,0:km,2,1)*                      &
                          dx(0:im,0:jm,0:km,3,2)                       &
                        - dx(0:im,0:jm,0:km,1,3)*                      &
                          dx(0:im,0:jm,0:km,2,2)*                      &
                          dx(0:im,0:jm,0:km,3,1)                       &
                        - dx(0:im,0:jm,0:km,1,2)*                      &
                          dx(0:im,0:jm,0:km,2,1)*                      &
                          dx(0:im,0:jm,0:km,3,3)                       &
                        - dx(0:im,0:jm,0:km,1,1)*                      &
                          dx(0:im,0:jm,0:km,2,3)*                      &
                          dx(0:im,0:jm,0:km,3,2)
   !
   call dataswap(jacob)
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! End of Calculating 
   ! geometrical Jacobian
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Calculating d<i,j,k>/d<x,y,z>
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   dxi=0.d0
   !
   allocate( phi(-hm:im+hm)  )
   do k=0,km
   do j=0,jm
     !
     phi(:)=0.5d0*(dx(-hm:im+hm,j,k,2,3)*x(-hm:im+hm,j,k,3)-           &
                   dx(-hm:im+hm,j,k,3,3)*x(-hm:im+hm,j,k,2))
     dxi(0:im,j,k,2,1)=dxi(0:im,j,k,2,1)+ddf(phi,cscheme,npdci,im,alfa,gci)
     !
     phi(:)=0.5d0*(dx(-hm:im+hm,j,k,3,3)*x(-hm:im+hm,j,k,1)-           &
                   dx(-hm:im+hm,j,k,1,3)*x(-hm:im+hm,j,k,3))
     dxi(0:im,j,k,2,2)=dxi(0:im,j,k,2,2)+ddf(phi,cscheme,npdci,im,alfa,gci)
     !
     phi(:)=0.5d0*(dx(-hm:im+hm,j,k,1,3)*x(-hm:im+hm,j,k,2)-           &
                   dx(-hm:im+hm,j,k,2,3)*x(-hm:im+hm,j,k,1))
     dxi(0:im,j,k,2,3)=dxi(0:im,j,k,2,3)+ddf(phi,cscheme,npdci,im,alfa,gci)
     !
     phi(:)=0.5d0*(dx(-hm:im+hm,j,k,3,2)*x(-hm:im+hm,j,k,2)-           &
                   dx(-hm:im+hm,j,k,2,2)*x(-hm:im+hm,j,k,3))
     dxi(0:im,j,k,3,1)=dxi(0:im,j,k,3,1)+ddf(phi,cscheme,npdci,im,alfa,gci)
     !
     phi(:)=0.5d0*(dx(-hm:im+hm,j,k,1,2)*x(-hm:im+hm,j,k,3)-           &
                   dx(-hm:im+hm,j,k,3,2)*x(-hm:im+hm,j,k,1))
     dxi(0:im,j,k,3,2)=dxi(0:im,j,k,3,2)+ddf(phi,cscheme,npdci,im,alfa,gci)
     !
     phi(:)=0.5d0*(dx(-hm:im+hm,j,k,2,2)*x(-hm:im+hm,j,k,1)-           &
                   dx(-hm:im+hm,j,k,1,2)*x(-hm:im+hm,j,k,2))
     dxi(0:im,j,k,3,3)=dxi(0:im,j,k,3,3)+ddf(phi,cscheme,npdci,im,alfa,gci)
     !
   enddo
   enddo
   deallocate( phi )
   !
   allocate( phi(-hm:jm+hm)  )
   do k=0,km
   do i=0,im
     !
     phi(:)=0.5d0*(dx(i,-hm:jm+hm,k,3,3)*x(i,-hm:jm+hm,k,2)-           &
                   dx(i,-hm:jm+hm,k,2,3)*x(i,-hm:jm+hm,k,3))
     dxi(i,0:jm,k,1,1)=dxi(i,0:jm,k,1,1)+ddf(phi,cscheme,npdcj,jm,alfa,gcj)
     !
     phi(:)=0.5d0*(dx(i,-hm:jm+hm,k,1,3)*x(i,-hm:jm+hm,k,3)-           &
                   dx(i,-hm:jm+hm,k,3,3)*x(i,-hm:jm+hm,k,1))
     dxi(i,0:jm,k,1,2)=dxi(i,0:jm,k,1,2)+ddf(phi,cscheme,npdcj,jm,alfa,gcj)
     !
     phi(:)=0.5d0*(dx(i,-hm:jm+hm,k,2,3)*x(i,-hm:jm+hm,k,1)-           &
                   dx(i,-hm:jm+hm,k,1,3)*x(i,-hm:jm+hm,k,2))
     dxi(i,0:jm,k,1,3)=dxi(i,0:jm,k,1,3)+ddf(phi,cscheme,npdcj,jm,alfa,gcj)
     !
     phi(:)=0.5d0*(dx(i,-hm:jm+hm,k,2,1)*x(i,-hm:jm+hm,k,3)-           &
                   dx(i,-hm:jm+hm,k,3,1)*x(i,-hm:jm+hm,k,2))
     dxi(i,0:jm,k,3,1)=dxi(i,0:jm,k,3,1)+ddf(phi,cscheme,npdcj,jm,alfa,gcj)
     !
     phi(:)=0.5d0*(dx(i,-hm:jm+hm,k,3,1)*x(i,-hm:jm+hm,k,1)-           &
                   dx(i,-hm:jm+hm,k,1,1)*x(i,-hm:jm+hm,k,3))
     dxi(i,0:jm,k,3,2)=dxi(i,0:jm,k,3,2)+ddf(phi,cscheme,npdcj,jm,alfa,gcj)
     !
     phi(:)=0.5d0*(dx(i,-hm:jm+hm,k,1,1)*x(i,-hm:jm+hm,k,2)-           &
                   dx(i,-hm:jm+hm,k,2,1)*x(i,-hm:jm+hm,k,1))
     dxi(i,0:jm,k,3,3)=dxi(i,0:jm,k,3,3)+ddf(phi,cscheme,npdcj,jm,alfa,gcj)
   enddo
   enddo
   deallocate( phi )
   !
   allocate( phi(-hm:km+hm)  )
   do j=0,jm
   do i=0,im
     !
     phi(:)=0.5d0*(dx(i,j,-hm:km+hm,2,2)*x(i,j,-hm:km+hm,3)-           &
                   dx(i,j,-hm:km+hm,3,2)*x(i,j,-hm:km+hm,2))
     dxi(i,j,0:km,1,1)=dxi(i,j,0:km,1,1)+ddf(phi,cscheme,npdck,km,alfa,gck)
     !
     phi(:)=0.5d0*(dx(i,j,-hm:km+hm,3,2)*x(i,j,-hm:km+hm,1)-           &
                   dx(i,j,-hm:km+hm,1,2)*x(i,j,-hm:km+hm,3))
     dxi(i,j,0:km,1,2)=dxi(i,j,0:km,1,2)+ddf(phi,cscheme,npdck,km,alfa,gck)
     !
     phi(:)=0.5d0*(dx(i,j,-hm:km+hm,1,2)*x(i,j,-hm:km+hm,2)-           &
                   dx(i,j,-hm:km+hm,2,2)*x(i,j,-hm:km+hm,1))
     dxi(i,j,0:km,1,3)=dxi(i,j,0:km,1,3)+ddf(phi,cscheme,npdck,km,alfa,gck)
     !
     phi(:)=0.5d0*(dx(i,j,-hm:km+hm,3,1)*x(i,j,-hm:km+hm,2)-           &
                   dx(i,j,-hm:km+hm,2,1)*x(i,j,-hm:km+hm,3))
     dxi(i,j,0:km,2,1)=dxi(i,j,0:km,2,1)+ddf(phi,cscheme,npdck,km,alfa,gck)
     !
     phi(:)=0.5d0*(dx(i,j,-hm:km+hm,1,1)*x(i,j,-hm:km+hm,3)-           &
                   dx(i,j,-hm:km+hm,3,1)*x(i,j,-hm:km+hm,1))
     dxi(i,j,0:km,2,2)=dxi(i,j,0:km,2,2)+ddf(phi,cscheme,npdck,km,alfa,gck)
     !
     phi(:)=0.5d0*(dx(i,j,-hm:km+hm,2,1)*x(i,j,-hm:km+hm,1)-           &
                   dx(i,j,-hm:km+hm,1,1)*x(i,j,-hm:km+hm,2))
     dxi(i,j,0:km,2,3)=dxi(i,j,0:km,2,3)+ddf(phi,cscheme,npdck,km,alfa,gck)
     !
   enddo
   enddo
   deallocate( phi )
   !
   call dataswap(dxi)
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Calculating geometrical
   ! metric identities
   !!!!!!!!!!!!!!!!!!!!!!!!!!!
   allocate(can(0:im,0:jm,0:km,1:3))
   !
   can=0.d0
   !
   do k=0,km
   do j=0,jm
    can(:,j,k,1)=can(:,j,k,1)+ddf(dxi(:,j,k,1,1),cscheme,npdci,im,alfa,gci)
    can(:,j,k,2)=can(:,j,k,2)+ddf(dxi(:,j,k,1,2),cscheme,npdci,im,alfa,gci)
    can(:,j,k,3)=can(:,j,k,3)+ddf(dxi(:,j,k,1,3),cscheme,npdci,im,alfa,gci)
   enddo
   enddo
   !
   do k=0,km
   do i=0,im
    can(i,:,k,1)=can(i,:,k,1)+ddf(dxi(i,:,k,2,1),cscheme,npdcj,jm,alfa,gcj)
    can(i,:,k,2)=can(i,:,k,2)+ddf(dxi(i,:,k,2,2),cscheme,npdcj,jm,alfa,gcj)
    can(i,:,k,3)=can(i,:,k,3)+ddf(dxi(i,:,k,2,3),cscheme,npdcj,jm,alfa,gcj)
   enddo
   enddo
   !
   do j=0,jm
   do i=0,im
    can(i,j,:,1)=can(i,j,:,1)+ddf(dxi(i,j,:,3,1),cscheme,npdck,km,alfa,gck)
    can(i,j,:,2)=can(i,j,:,2)+ddf(dxi(i,j,:,3,2),cscheme,npdck,km,alfa,gck)
    can(i,j,:,3)=can(i,j,:,3)+ddf(dxi(i,j,:,3,3),cscheme,npdck,km,alfa,gck)
   enddo
   enddo
   !
   can1av=0.d0
   can2av=0.d0
   can3av=0.d0
   do k=1,km
   do j=1,jm
   do i=1,im
     can1av=can1av+can(i,j,k,1)
     can2av=can2av+can(i,j,k,2)
     can3av=can3av+can(i,j,k,3)
   end do
   end do
   end do
   can1av=psum(can1av)/real(ia*ja*ka,8)
   can2av=psum(can2av)/real(ia*ja*ka,8)
   can3av=psum(can3av)/real(ia*ja*ka,8)
   !
   can1var=0.d0
   can2var=0.d0
   can3var=0.d0
   !
   do k=1,km
   do j=1,jm
   do i=1,im
     can1var=can1var+(can(i,j,k,1)-can1av)**2
     can2var=can2var+(can(i,j,k,2)-can2av)**2
     can3var=can3var+(can(i,j,k,3)-can3av)**2
   end do
   end do
   end do
   can1var=psum(can1var)/real(ia*ja*ka,8)
   can2var=psum(can2var)/real(ia*ja*ka,8)
   can3var=psum(can3var)/real(ia*ja*ka,8)
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! End ofCalculating 
   ! geometrical metric 
   ! identities
   !!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! normlizing 
   ! d<i,j,k>/d<x,y,z> and
   ! geometrical vector.
   !!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   do j=1,3
   do i=1,3
     dxi(0:im,0:jm,0:km,i,j)=dxi(0:im,0:jm,0:km,i,j)/jacob(0:im,0:jm,0:km)
   enddo
   enddo
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! End of normlizing
   ! d<i,j,k>/d<x,y,z> and
   ! geometrical vector.
   !!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   call dataswap(dxi)
   !
   xmax=maxval(x(0:im,0:jm,0:km,1))
   ymax=maxval(x(0:im,0:jm,0:km,2))
   zmax=maxval(x(0:im,0:jm,0:km,3))
   !
   xmin=minval(x(0:im,0:jm,0:km,1))
   ymin=minval(x(0:im,0:jm,0:km,2))
   zmin=minval(x(0:im,0:jm,0:km,3))
   !
   xmax=pmax(xmax)
   ymax=pmax(ymax)
   zmax=pmax(zmax)
   !
   xmin=pmin(xmin)
   ymin=pmin(ymin)
   zmin=pmin(zmin)
   !
   !
   if(lio) then
     !
     write(*,'(2X,62A)')('-',i=1,62)
     write(*,'(2X,A)')'                    *** Grids Information *** '
     write(*,'(2X,62A)')('-',i=1,62)
     write(*,'(3X,62A)')'       xmin      xmax      ymin      ymax      zmin      zmax'
     write(*,"(4X,6(F10.6))")xmin,xmax,ymin,ymax,zmin,zmax
     write(*,'(2X,62A)')('-',i=1,62)
     write(*,'(2X,A)')'                   *** Averaged of Identity ***'
     write(*,"(1X,3(1X,E20.7E3))")can1av,can2av,can3av
     write(*,'(2X,62A)')('-',i=1,62)
     write(*,'(2X,A)')'                   *** Variance of Identity ***'
     write(*,"(1X,3(1X,E20.7E3))")can1var,can2var,can3var
     write(*,'(2X,62A)')('-',i=1,62)
     !
     if(can1av>1d-15 .or. can2av>1d-15 .or. can3av>1d-15) then
       write(*,*)' !! Warning: Averaged Grids Identity is too large'
       write(*,'(2X,62A)')('-',i=1,62)
     end if
     if(can1var>1d-15 .or. can2var>1d-15 .or. can3var>1d-15) then
       write(*,*)' !! Warning: Variance of Grids Identity is too large'
       write(*,'(2X,62A)')('-',i=1,62)
     end if
     !
     print*,' ** Geometrical parameters calculated'
   end if
   !
   deallocate(dx,gci,gcj,gck,can)
   !
    ! if(mpirank==0) then
    !   do k=0,km
    !     print*,x(1,1,k,3),dx(1,1,k,3,3),x(1,1,k,3)-x(1,1,k-1,3)
    !   enddo
    ! endif
    ! call tecbin('testout/tecgrid'//mpirankname//'.plt',                &
    !                                   x(0:im,-hm:jm+hm,0:km,1),'x',    &
    !                                   x(0:im,-hm:jm+hm,0:km,2),'y',    &
    !                                   x(0:im,-hm:jm+hm,0:km,3),'z',    &
    !                                  dx(0:im,-hm:jm+hm,0:km,1,1),'dxdi')
    ! !
  end subroutine geomcal
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! End of the subroutine GeomCal.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  subroutine gradcal
    !
    use commvar,   only : im,jm,km,npdci,npdcj,npdck,conschm
    use commarray, only : x,q,dxi
    use commfunc,  only : ddf
    !
    ! local data
    integer :: i,j,k
    real(8),allocatable :: dq(:)
    !
    do k=0,km
    do j=0,jm
    do i=0,im
      q(i,j,k,1)=sin(x(i,j,k,3)*2.d0)
    enddo
    enddo
    enddo
    !
    call dataswap(q)
    !
    ! allocate(dq(0:im))
    ! dq=ddf(q(:,1,1,1),conschm,npdci,im,alfa_con,cci)*dxi(:,1,1,1,1)
    !
    ! allocate(dq(0:jm))
    ! dq=ddf(q(1,:,1,1),conschm,npdcj,jm,alfa_con,ccj)*dxi(1,0:jm,1,2,2)
    !
    allocate(dq(0:km))
    dq=ddf(q(1,1,:,1),conschm,npdck,km,alfa_con,cck)*dxi(1,1,0:km,3,3)
    !
    open(18,file='testout/profile'//mpirankname//'.dat')
    do k=0,km
      write(18,*)x(1,1,k,3),cos(x(1,1,k,3)*2.d0)*2.d0,dq(k)
    enddo
    close(18)
    print*,' << testout/profile',mpirankname,'.dat'
    !
  end subroutine gradcal
  !
end module solver
!+---------------------------------------------------------------------+
!| The end of the module commarray.                                    |
!+---------------------------------------------------------------------+