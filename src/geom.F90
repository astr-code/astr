!+---------------------------------------------------------------------+
!| This module contains subroutines of calculating geometrically       |
!| parameters                                                          |
!+---------------------------------------------------------------------+
!| CHANGE RECORD                                                       |
!| -------------                                                       |
!| 03-07-2021  | Created by J. Fang                                    |
!+---------------------------------------------------------------------+
module geom
  !
  use constdef
  use parallel, only : mpirankname,mpistop,mpirank,lio,dataswap,       &
                       datasync,ptime,irk,jrk,krk
  use commvar,  only : ndims,ks,ke,hm,hm,lfftk,ctime
  !
  implicit none
  !
  contains
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to calculate geometrical parameters            |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 03-07-2021  | Created by J. Fang @ Warrington                     |
  !+-------------------------------------------------------------------+
  subroutine geomcal
    !
    use commvar,   only : limmbou,immbody
    use parallel,  only : bcast
    use stlaio,    only : stla_write
    !
    call meshgeom
    !
    if(limmbou) then
      call solidgeom
      !
      call bcast(immbody)
      !
      ! call stla_write('solid_real'//mpirankname//'.stl',immbody)
      !
    endif
    !
  end subroutine geomcal
  !+-------------------------------------------------------------------+
  !| The end of the subroutine rhscal.                                 |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to calculate solid's geometrical          |
  !| parameters                                                        |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 02-Jul-2021: Created by J. Fang @ Appleton                        |
  !+-------------------------------------------------------------------+
  subroutine solidgeom
    !
    use commtype,  only : solid,triangle
    use commvar,   only : immbody,nsolid
    !
    ! local data
    integer :: js
    type(solid),pointer :: psolid
    !
    if(mpirank==0) then
      !
      do js=1,nsolid
        !
        call solidrange(immbody(js))
        !
      enddo
      !
      call solidresc(immbody(1),0.025d0)
      !
      call solidshif(immbody(1),x=1.7415d0,y=-0.3838975d0,z=0.9457d0)
      !
      call solidrange(immbody(1))
      !
    endif
    !
    return
    !
  end subroutine solidgeom
  !+-------------------------------------------------------------------+
  !| The end of the subroutine solidgeom.                              |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to obtain the range of the solid          |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 03-Jul-2021: Created by J. Fang @ Appleton                        |
  !+-------------------------------------------------------------------+
  subroutine solidrange(asolid)
    !
    use commtype,  only : solid,triangle
    use commvar,   only : xmax,xmin,ymax,ymin,zmax,zmin
    !
    ! arguments
    type(solid),intent(inout) :: asolid
    !
    ! local data
    integer :: i,jf
    !
    asolid%xmin=1.d10
    asolid%xmax=-1.d10
    do jf=1,asolid%num_face
      do i=1,3
        asolid%xmin(i)=min(asolid%xmin(i),                         &
                           asolid%face(jf)%a(i),                   &
                           asolid%face(jf)%b(i),                   &
                           asolid%face(jf)%c(i) )
        asolid%xmax(i)=max(asolid%xmax(i),                         &
                           asolid%face(jf)%a(i),                   &
                           asolid%face(jf)%b(i),                   &
                           asolid%face(jf)%c(i) )
      enddo
    enddo
    !
    write(*,'(A,A)')      '  solid: ',asolid%name
    write(*,'(2(A,E15.7E3))')'    x: ',asolid%xmin(1),'~',asolid%xmax(1)
    write(*,'(2(A,E15.7E3))')'    y: ',asolid%xmin(2),'~',asolid%xmax(2)
    write(*,'(2(A,E15.7E3))')'    z: ',asolid%xmin(3),'~',asolid%xmax(3)
    !
    if(asolid%xmin(1)<xmin) then
      write(*,'(A,2(E15.7E3))')'  !! WARNING: solid outof computational domain'
      write(*,'(2(A,E15.7E3))')'  !! solid xmin: ',asolid%xmin(1),      &
                                  ' domain xmin: ',xmin
    endif
    if(asolid%xmin(2)<ymin) then
      write(*,'(A,2(E15.7E3))')'  !! WARNING: solid outof computational domain'
      write(*,'(2(A,E15.7E3))')'  !! solid ymin: ',asolid%xmin(2),      &
                                  ' domain ymin: ',ymin
    endif
    if(asolid%xmin(3)<zmin) then
      write(*,'(A,2(E15.7E3))')'  !! WARNING: solid outof computational domain'
      write(*,'(2(A,E15.7E3))')'  !! solid zmin: ',asolid%xmin(3),      &
                                  ' domain zmin: ',zmin
    endif
    !
    if(asolid%xmax(1)>xmax) then
      write(*,'(A,2(E15.7E3))')'  !! WARNING: solid outof computational domain'
      write(*,'(2(A,E15.7E3))')'  !! solid xmax: ',asolid%xmax(1),      &
                                  ' domain xmax: ',xmax
    endif
    if(asolid%xmax(2)>ymax) then
      write(*,'(A,2(E15.7E3))')'  !! WARNING: solid outof computational domain'
      write(*,'(2(A,E15.7E3))')'  !! solid ymax: ',asolid%xmax(2),      &
                                  ' domain ymax: ',ymax
    endif
    if(asolid%xmax(3)>zmax) then
      write(*,'(A,2(E15.7E3))')'  !! WARNING: solid outof computational domain'
      write(*,'(2(A,E15.7E3))')'  !! solid zmax: ',asolid%xmax(3),      &
                                  ' domain zmax: ',zmax
    endif
    !
  end subroutine solidrange
  !+-------------------------------------------------------------------+
  !| The end of the subroutine solidrange.                             |
  !+-------------------------------------------------------------------+
  !!
  !+-------------------------------------------------------------------+
  !| This subroutine is used to rescale size of the solid.             |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 03-Jul-2021: Created by J. Fang @ Appleton                        |
  !+-------------------------------------------------------------------+
  subroutine solidresc(asolid,scale)
    !
    use commtype,  only : solid,triangle
    !
    ! arguments
    type(solid),intent(inout) :: asolid
    real(8),intent(in) :: scale
    !
    ! local data
    integer :: i,jf
    !
    do jf=1,asolid%num_face
      !
      asolid%face(jf)%a(:)=asolid%face(jf)%a(:)*scale
      asolid%face(jf)%b(:)=asolid%face(jf)%b(:)*scale
      asolid%face(jf)%c(:)=asolid%face(jf)%c(:)*scale
      !
    enddo
    !
    write(*,'(3A,E15.7E3)')'  ** solid ',trim(asolid%name),            &
                                                 ' is rescaled by',scale
    !
  end subroutine solidresc
  !+-------------------------------------------------------------------+
  !| The end of the subroutine solidresc.                              |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to rescale size of the solid.             |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 03-Jul-2021: Created by J. Fang @ Appleton                        |
  !+-------------------------------------------------------------------+
  subroutine solidshif(asolid,x,y,z)
    !
    use commtype,  only : solid,triangle
    !
    ! arguments
    type(solid),intent(inout) :: asolid
    real(8),intent(in) :: x,y,z
    !
    ! local data
    integer :: i,jf
    !
    do jf=1,asolid%num_face
      !
      asolid%face(jf)%a(1)=asolid%face(jf)%a(1)+x
      asolid%face(jf)%a(2)=asolid%face(jf)%a(2)+y
      asolid%face(jf)%a(3)=asolid%face(jf)%a(3)+z
      asolid%face(jf)%b(1)=asolid%face(jf)%b(1)+x
      asolid%face(jf)%b(2)=asolid%face(jf)%b(2)+y
      asolid%face(jf)%b(3)=asolid%face(jf)%b(3)+z
      asolid%face(jf)%c(1)=asolid%face(jf)%c(1)+x
      asolid%face(jf)%c(2)=asolid%face(jf)%c(2)+y
      asolid%face(jf)%c(3)=asolid%face(jf)%c(3)+z
      !
    enddo
    !
    write(*,'(3(A),3(E15.7E3))')'  ** solid ',trim(asolid%name),       &
                                                ' is shifted by ',x,y,z
    !
  end subroutine solidshif
  !+-------------------------------------------------------------------+
  !| The end of the subroutine solidshif.                              |
  !+-------------------------------------------------------------------+
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! This subtoutine is used to calculate geometrical transform matrix
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Writen by Fang Jian, 2008-09-22.
  ! Change name by Fang Jian, 2008-09-22.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine meshgeom
    !
    use commvar,   only : ia,ja,ka,im,jm,km,hm,npdci,npdcj,npdck,hm,   &
                          xmax,xmin,ymax,ymin,zmax,zmin,voldom,difschm
    use commarray, only : x,jacob,dxi,celvol
    use parallel,  only : gridsendrecv,jsize,ksize,psum,pmax,pmin
    use commfunc,  only : coeffcompac,ptds_ini,ddfc,volhex,arquad
    use tecio
    use bc,       only : geombc,xyzbc
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
    call xyzbc
    !
    ! cscheme='442e'
    cscheme=difschm
    !
    read(cscheme(1:3),*) nscheme
    !
    if(cscheme(4:4)=='c') then
      ! a compact scheme is used
      !
      alfa=coeffcompac(nscheme)
      !
      call ptds_ini(gci,alfa,im,npdci)
      call ptds_ini(gcj,alfa,jm,npdcj)
      call ptds_ini(gck,alfa,km,npdck)
      !
    endif
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Calculating d<x,y,z>/d<i,j,k>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do k=0,km
    do j=0,jm
      do m=1,3
        dx(0:im,j,k,m,1)=ddfc(x(:,j,k,m),cscheme,npdci,im,alfa,gci)
      enddo
    enddo
    enddo
    !
    do k=0,km
    do i=0,im
      !
      do m=1,3
        dx(i,0:jm,k,m,2)=ddfc(x(i,:,k,m),cscheme,npdcj,jm,alfa,gcj)
      enddo
      !
    enddo
    enddo
    !
    do j=0,jm
    do i=0,im
      !
      if(lfftk) then
        !
        do m=1,2
          dx(i,j,0:km,m,3)=ddfc(x(i,j,:,m),cscheme,npdck,km,alfa,gck,lfft=lfftk)
        enddo
        dx(i,j,0:km,3,3)=x(i,j,1,3)-x(i,j,0,3)
      else
        do m=1,3
          dx(i,j,0:km,m,3)=ddfc(x(i,j,:,m),cscheme,npdck,km,alfa,gck)
        enddo
      endif
      !
    enddo
    enddo
    !
    if(ndims==1) then
      !
      j=0
      k=0
      !
      voldom=0.d0
      do i=1,im
        celvol(i,j,k)=abs(x(i,j,k,1)-x(i-1,j,k,1))
        voldom=voldom+celvol(i,j,k)
      enddo
      voldom=psum(voldom)
      !
    elseif(ndims==2) then
      k=0
      !
      voldom=0.d0
      do j=1,jm
      do i=1,im
        celvol(i,j,k)=arquad( x(i-1,j-1,k,:),   x(i,j-1,k,:),          &
                              x(i,j,k,:),       x(i-1,j,k,:) )
        voldom=voldom+celvol(i,j,k)
      enddo
      enddo
      voldom=psum(voldom)
    elseif(ndims==3) then
      voldom=0.d0
      do k=1,km
      do j=1,jm
      do i=1,im
        celvol(i,j,k)=volhex( x(i-1,j-1,k-1,:), x(i,j-1,k-1,:),        &
                              x(i,j-1,k,:),     x(i-1,j-1,k,:),        &
                              x(i-1,j,k-1,:),   x(i,j,k-1,:)  ,        &
                              x(i,j,k,:),       x(i-1,j,k,:) )
        voldom=voldom+celvol(i,j,k)
      enddo
      enddo
      enddo
      voldom=psum(voldom)
    endif
    !
    if(lio) print*,' ** total volume of the domain is: ',voldom
    !
    if(lio) print*,' ** dxyz/dijk calculated.'
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! end of cal d<x,y,z>/d<i,j,k>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    call dataswap(dx,subtime=ctime(7))
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Calculating geometrical 
    ! Jacobian
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if(ndims==3) then
      jacob(0:im,0:jm,0:km)=dx(0:im,0:jm,0:km,1,1)*                    &
                            dx(0:im,0:jm,0:km,2,2)*                    &
                            dx(0:im,0:jm,0:km,3,3)                     &
                          + dx(0:im,0:jm,0:km,1,2)*                    &
                            dx(0:im,0:jm,0:km,2,3)*                    &
                            dx(0:im,0:jm,0:km,3,1)                     &
                          + dx(0:im,0:jm,0:km,1,3)*                    &
                            dx(0:im,0:jm,0:km,2,1)*                    &
                            dx(0:im,0:jm,0:km,3,2)                     &
                          - dx(0:im,0:jm,0:km,1,3)*                    &
                            dx(0:im,0:jm,0:km,2,2)*                    &
                            dx(0:im,0:jm,0:km,3,1)                     &
                          - dx(0:im,0:jm,0:km,1,2)*                    &
                            dx(0:im,0:jm,0:km,2,1)*                    &
                            dx(0:im,0:jm,0:km,3,3)                     &
                          - dx(0:im,0:jm,0:km,1,1)*                    &
                            dx(0:im,0:jm,0:km,2,3)*                    &
                            dx(0:im,0:jm,0:km,3,2)
    elseif(ndims==2) then  
      jacob(0:im,0:jm,0:km)=dx(0:im,0:jm,0:km,1,1)*                    &
                            dx(0:im,0:jm,0:km,2,2)                     &
                          - dx(0:im,0:jm,0:km,1,2)*                    &
                            dx(0:im,0:jm,0:km,2,1)
    elseif(ndims==1) then  
      jacob(0:im,0:jm,0:km)=dx(0:im,0:jm,0:km,1,1)*                    &
                            dx(0:im,0:jm,0:km,2,2)
    else
      stop ' !! ndimes not defined at jacob calculation.'
    endif
    !
    call dataswap(jacob,subtime=ctime(7))
    !
    call datasync(jacob(0:im,0:jm,0:km),subtime=ctime(7))
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
    if(ndims==3) then
      !
      allocate( phi(-hm:im+hm)  )
      do k=0,km
      do j=0,jm
        !
        phi(:)=0.5d0*(dx(-hm:im+hm,j,k,2,3)*x(-hm:im+hm,j,k,3)-           &
                      dx(-hm:im+hm,j,k,3,3)*x(-hm:im+hm,j,k,2))
        dxi(0:im,j,k,2,1)=dxi(0:im,j,k,2,1)+ddfc(phi,cscheme,npdci,im,alfa,gci)
        !
        phi(:)=0.5d0*(dx(-hm:im+hm,j,k,3,3)*x(-hm:im+hm,j,k,1)-           &
                      dx(-hm:im+hm,j,k,1,3)*x(-hm:im+hm,j,k,3))
        dxi(0:im,j,k,2,2)=dxi(0:im,j,k,2,2)+ddfc(phi,cscheme,npdci,im,alfa,gci)
        !
        phi(:)=0.5d0*(dx(-hm:im+hm,j,k,1,3)*x(-hm:im+hm,j,k,2)-           &
                      dx(-hm:im+hm,j,k,2,3)*x(-hm:im+hm,j,k,1))
        dxi(0:im,j,k,2,3)=dxi(0:im,j,k,2,3)+ddfc(phi,cscheme,npdci,im,alfa,gci)
        !
        phi(:)=0.5d0*(dx(-hm:im+hm,j,k,3,2)*x(-hm:im+hm,j,k,2)-           &
                      dx(-hm:im+hm,j,k,2,2)*x(-hm:im+hm,j,k,3))
        dxi(0:im,j,k,3,1)=dxi(0:im,j,k,3,1)+ddfc(phi,cscheme,npdci,im,alfa,gci)
        !
        phi(:)=0.5d0*(dx(-hm:im+hm,j,k,1,2)*x(-hm:im+hm,j,k,3)-           &
                      dx(-hm:im+hm,j,k,3,2)*x(-hm:im+hm,j,k,1))
        dxi(0:im,j,k,3,2)=dxi(0:im,j,k,3,2)+ddfc(phi,cscheme,npdci,im,alfa,gci)
        !
        phi(:)=0.5d0*(dx(-hm:im+hm,j,k,2,2)*x(-hm:im+hm,j,k,1)-           &
                      dx(-hm:im+hm,j,k,1,2)*x(-hm:im+hm,j,k,2))
        dxi(0:im,j,k,3,3)=dxi(0:im,j,k,3,3)+ddfc(phi,cscheme,npdci,im,alfa,gci)
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
        dxi(i,0:jm,k,1,1)=dxi(i,0:jm,k,1,1)+ddfc(phi,cscheme,npdcj,jm,alfa,gcj)
        !
        phi(:)=0.5d0*(dx(i,-hm:jm+hm,k,1,3)*x(i,-hm:jm+hm,k,3)-           &
                      dx(i,-hm:jm+hm,k,3,3)*x(i,-hm:jm+hm,k,1))
        dxi(i,0:jm,k,1,2)=dxi(i,0:jm,k,1,2)+ddfc(phi,cscheme,npdcj,jm,alfa,gcj)
        !
        phi(:)=0.5d0*(dx(i,-hm:jm+hm,k,2,3)*x(i,-hm:jm+hm,k,1)-           &
                      dx(i,-hm:jm+hm,k,1,3)*x(i,-hm:jm+hm,k,2))
        dxi(i,0:jm,k,1,3)=dxi(i,0:jm,k,1,3)+ddfc(phi,cscheme,npdcj,jm,alfa,gcj)
        !
        phi(:)=0.5d0*(dx(i,-hm:jm+hm,k,2,1)*x(i,-hm:jm+hm,k,3)-           &
                      dx(i,-hm:jm+hm,k,3,1)*x(i,-hm:jm+hm,k,2))
        dxi(i,0:jm,k,3,1)=dxi(i,0:jm,k,3,1)+ddfc(phi,cscheme,npdcj,jm,alfa,gcj)
        !
        phi(:)=0.5d0*(dx(i,-hm:jm+hm,k,3,1)*x(i,-hm:jm+hm,k,1)-           &
                      dx(i,-hm:jm+hm,k,1,1)*x(i,-hm:jm+hm,k,3))
        dxi(i,0:jm,k,3,2)=dxi(i,0:jm,k,3,2)+ddfc(phi,cscheme,npdcj,jm,alfa,gcj)
        !
        phi(:)=0.5d0*(dx(i,-hm:jm+hm,k,1,1)*x(i,-hm:jm+hm,k,2)-           &
                      dx(i,-hm:jm+hm,k,2,1)*x(i,-hm:jm+hm,k,1))
        dxi(i,0:jm,k,3,3)=dxi(i,0:jm,k,3,3)+ddfc(phi,cscheme,npdcj,jm,alfa,gcj)
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
        if(lfftk) then
          dxi(i,j,0:km,1,1)=dxi(i,j,0:km,1,1)+phi(1)-phi(0)
        else
          dxi(i,j,0:km,1,1)=dxi(i,j,0:km,1,1)+ddfc(phi,cscheme,npdck,km,alfa,gck,lfft=lfftk)
        endif
        !
        phi(:)=0.5d0*(dx(i,j,-hm:km+hm,3,2)*x(i,j,-hm:km+hm,1)-           &
                      dx(i,j,-hm:km+hm,1,2)*x(i,j,-hm:km+hm,3))
        if(lfftk) then
          dxi(i,j,0:km,1,2)=dxi(i,j,0:km,1,2)+phi(1)-phi(0)
        else
          dxi(i,j,0:km,1,2)=dxi(i,j,0:km,1,2)+ddfc(phi,cscheme,npdck,km,alfa,gck,lfft=lfftk)
        endif
        !
        phi(:)=0.5d0*(dx(i,j,-hm:km+hm,1,2)*x(i,j,-hm:km+hm,2)-           &
                      dx(i,j,-hm:km+hm,2,2)*x(i,j,-hm:km+hm,1))
        if(lfftk) then
          dxi(i,j,0:km,1,3)=dxi(i,j,0:km,1,3)+phi(1)-phi(0)
        else
          dxi(i,j,0:km,1,3)=dxi(i,j,0:km,1,3)+ddfc(phi,cscheme,npdck,km,alfa,gck,lfft=lfftk)
        endif
        !
        phi(:)=0.5d0*(dx(i,j,-hm:km+hm,3,1)*x(i,j,-hm:km+hm,2)-           &
                      dx(i,j,-hm:km+hm,2,1)*x(i,j,-hm:km+hm,3))
        if(lfftk) then
          dxi(i,j,0:km,2,1)=dxi(i,j,0:km,2,1)+phi(1)-phi(0)
        else
          dxi(i,j,0:km,2,1)=dxi(i,j,0:km,2,1)+ddfc(phi,cscheme,npdck,km,alfa,gck,lfft=lfftk)
        endif
        !
        phi(:)=0.5d0*(dx(i,j,-hm:km+hm,1,1)*x(i,j,-hm:km+hm,3)-           &
                      dx(i,j,-hm:km+hm,3,1)*x(i,j,-hm:km+hm,1))
        if(lfftk) then
          dxi(i,j,0:km,2,2)=dxi(i,j,0:km,2,2)+phi(1)-phi(0)
        else
          dxi(i,j,0:km,2,2)=dxi(i,j,0:km,2,2)+ddfc(phi,cscheme,npdck,km,alfa,gck,lfft=lfftk)
        endif
        !
        phi(:)=0.5d0*(dx(i,j,-hm:km+hm,2,1)*x(i,j,-hm:km+hm,1)-           &
                      dx(i,j,-hm:km+hm,1,1)*x(i,j,-hm:km+hm,2))
        if(lfftk) then
          dxi(i,j,0:km,2,3)=dxi(i,j,0:km,2,3)+phi(1)-phi(0)
        else
          dxi(i,j,0:km,2,3)=dxi(i,j,0:km,2,3)+ddfc(phi,cscheme,npdck,km,alfa,gck,lfft=lfftk)
        endif
        !
      enddo
      enddo
      deallocate( phi )
      !
    elseif(ndims==2) then
      !
      dxi(0:im,0:jm,0:km,1,1)= dx(0:im,0:jm,0:km,2,2)  
      dxi(0:im,0:jm,0:km,1,2)=-dx(0:im,0:jm,0:km,1,2)  
      dxi(0:im,0:jm,0:km,2,1)=-dx(0:im,0:jm,0:km,2,1)  
      dxi(0:im,0:jm,0:km,2,2)= dx(0:im,0:jm,0:km,1,1) 
      !
      dxi(0:im,0:jm,0:km,1,3)=0.d0
      dxi(0:im,0:jm,0:km,2,3)=0.d0
      dxi(0:im,0:jm,0:km,3,:)=0.d0
      !
    elseif(ndims==1) then
      !
      dxi(0:im,0:jm,0:km,1,1)= dx(0:im,0:jm,0:km,2,2)    
      dxi(0:im,0:jm,0:km,2,2)= dx(0:im,0:jm,0:km,1,1) 
      !
      dxi(0:im,0:jm,0:km,1,2)=0.d0
      dxi(0:im,0:jm,0:km,2,1)=0.d0
      dxi(0:im,0:jm,0:km,1,3)=0.d0
      dxi(0:im,0:jm,0:km,2,3)=0.d0
      dxi(0:im,0:jm,0:km,3,:)=0.d0
      !
    endif
    !
    call dataswap(dxi,subtime=ctime(7))
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
      can(:,j,k,1)=can(:,j,k,1)+ddfc(dxi(:,j,k,1,1),cscheme,npdci,im,alfa,gci)
      can(:,j,k,2)=can(:,j,k,2)+ddfc(dxi(:,j,k,1,2),cscheme,npdci,im,alfa,gci)
      can(:,j,k,3)=can(:,j,k,3)+ddfc(dxi(:,j,k,1,3),cscheme,npdci,im,alfa,gci)
    enddo
    enddo
    !
    do k=0,km
    do i=0,im
      can(i,:,k,1)=can(i,:,k,1)+ddfc(dxi(i,:,k,2,1),cscheme,npdcj,jm,alfa,gcj)
      can(i,:,k,2)=can(i,:,k,2)+ddfc(dxi(i,:,k,2,2),cscheme,npdcj,jm,alfa,gcj)
      can(i,:,k,3)=can(i,:,k,3)+ddfc(dxi(i,:,k,2,3),cscheme,npdcj,jm,alfa,gcj)
    enddo
    enddo
    !
    do j=0,jm
    do i=0,im
      can(i,j,:,1)=can(i,j,:,1)+ddfc(dxi(i,j,:,3,1),cscheme,npdck,km,alfa,gck,lfft=lfftk)
      can(i,j,:,2)=can(i,j,:,2)+ddfc(dxi(i,j,:,3,2),cscheme,npdck,km,alfa,gck,lfft=lfftk)
      can(i,j,:,3)=can(i,j,:,3)+ddfc(dxi(i,j,:,3,3),cscheme,npdck,km,alfa,gck,lfft=lfftk)
      !
    enddo
    enddo
    !
    can1av=0.d0
    can2av=0.d0
    can3av=0.d0
    if(ndims==3) then
      do k=1,km
      do j=1,jm
      do i=1,im
        can1av=can1av+can(i,j,k,1)
        can2av=can2av+can(i,j,k,2)
        can3av=can3av+can(i,j,k,3)
        !
        if(isnan(can(i,j,k,1))) print*,mpirank,'|-1',i,j,k,can(i,j,k,:)
        if(isnan(can(i,j,k,2))) print*,mpirank,'|-2',i,j,k,can(i,j,k,:)
        if(isnan(can(i,j,k,3))) print*,mpirank,'|-3',i,j,k,can(i,j,k,:)
        !
      end do
      end do
      end do
      can1av=psum(can1av)/real(ia*ja*ka,8)
      can2av=psum(can2av)/real(ia*ja*ka,8)
      can3av=psum(can3av)/real(ia*ja*ka,8)
    elseif(ndims==2) then
      k=0
      do j=1,jm
      do i=1,im
        can1av=can1av+can(i,j,k,1)
        can2av=can2av+can(i,j,k,2)
        can3av=can3av+can(i,j,k,3)
        !
        if(isnan(can(i,j,k,1))) print*,mpirank,'|-1',i,j,k,can(i,j,k,:)
        if(isnan(can(i,j,k,2))) print*,mpirank,'|-2',i,j,k,can(i,j,k,:)
        if(isnan(can(i,j,k,3))) print*,mpirank,'|-3',i,j,k,can(i,j,k,:)
        !
      end do
      end do
      can1av=psum(can1av)/real(ia*ja,8)
      can2av=psum(can2av)/real(ia*ja,8)
      can3av=psum(can3av)/real(ia*ja,8)
    endif
    !
    can1var=0.d0
    can2var=0.d0
    can3var=0.d0
    !
    if(ndims==3) then
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
    elseif(ndims==2) then
      do j=1,jm
      do i=1,im
        can1var=can1var+(can(i,j,k,1)-can1av)**2
        can2var=can2var+(can(i,j,k,2)-can2av)**2
        can3var=can3var+(can(i,j,k,3)-can3av)**2
      end do
      end do
      can1var=psum(can1var)/real(ia*ja,8)
      can2var=psum(can2var)/real(ia*ja,8)
      can3var=psum(can3var)/real(ia*ja,8)
    endif
    !
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
    call dataswap(dxi,subtime=ctime(7))
    !
    do j=1,3
    do i=1,3
      call datasync(dxi(0:im,0:jm,0:km,i,j),subtime=ctime(7))
    enddo
    enddo
    !
    call geombc
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
    if(lio) then
      !
      write(*,'(2X,62A)')('-',i=1,62)
      write(*,'(2X,A)')'                    *** Grids Information *** '
      write(*,'(2X,62A)')('-',i=1,62)
      write(*,'(3X,62A)')'       xmin      xmax      ymin      ymax      zmin      zmax'
      write(*,"(4X,6(F10.3))")xmin,xmax,ymin,ymax,zmin,zmax
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
      print*,' ** geometrical parameters calculated'
    end if
    !
    ! if(mpirank==0) then
    !   do k=0,km
    !     print*,x(1,1,k,3),dx(1,1,k,3,3),x(1,1,k,3)-x(1,1,k-1,3)
    !   enddo
    ! ! endif
    ! call tecbin('testout/tecgrid'//mpirankname//'.plt',                &
    !                                   x(0:im,0:jm,0:km,1),'x',    &
    !                                   x(0:im,0:jm,0:km,2),'y',    &
    !                                   x(0:im,0:jm,0:km,3),'z',    &
    !                                 dxi(0:im,0:jm,0:km,1,1),'dxdi',    &
    !                                 dxi(0:im,0:jm,0:km,2,2),'dydj',    &
    !                                 dxi(0:im,0:jm,0:km,3,3),'dzdk')
    ! ! ! 
    deallocate(dx,can)
    !
    !
    if(allocated(gci)) deallocate(gci)
    if(allocated(gcj)) deallocate(gcj)
    if(allocated(gck)) deallocate(gck)
    !
    
    return
    !
    ! call mpistop
    !
  end subroutine meshgeom
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! End of the subroutine meshgeom.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
end module geom