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
                       datasync,ptime,irk,jrk,krk,ig0,jg0,kg0,         &
                       mpirankmax,psum,pmax
  use commvar,  only : ndims,ks,ke,hm,hm,lfftk,ctime,im,jm,km,         &
                       npdci,npdcj,npdck
  use tecio
  use stlaio,  only: get_unit
  use commcal, only: ijkcellin,ijkin
  !
  implicit none
  !
  interface pointintriangle
    module procedure pointintriangle_tri
    module procedure pointintriangle_nodes
  end interface
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
    use commvar,   only : limmbou,immbody,nsolid
    use parallel,  only : bcast
    use stlaio,    only : stla_write
    !
    call gridgeom
    !
    if(limmbou) then
      !
      call solidgeom
      !
      call bcast(nsolid)
      !
      call bcast(immbody)
      !
      call immsgrid
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
    use commvar,   only : xmax,xmin,ymax,ymin,zmax,zmin,               &
                          immbody,nsolid,ndims
    use tecio,     only : tecsolid
    !
    ! local data
    integer :: js,fh,ios,i
    type(solid),pointer :: psolid
    character(len=64) :: infile,head
    logical :: lexist,lread
    logical :: lshift=.false.,lscale=.false.,lrotate=.false.
    real(8) :: xcen(3),scale,theta,rot_vec(3)
    !
    if(mpirank==0) then
      !
      infile='datin/stltransform.dat'
      inquire(file=trim(infile), exist=lexist)
      !
      if(lexist) then
        !
        fh=get_unit()
        open(fh,file=trim(infile),action='read')
        lread=.true.
        ios=0
        do while(lread .and. ios==0)
          !
          read(fh,*,iostat=ios)head
          !
          if(ios.ne.0) exit
          !
          select case(trim(head))
            !
            case('center')
              !
              backspace(fh)
              read(fh,*,iostat=ios)head,(xcen(i),i=1,3)
              lshift=.true.
              !
            case('rescale')
              !
              backspace(fh)
              read(fh,*,iostat=ios)head,scale
              lscale=.true.
              !
            case('rotate')
              !
              backspace(fh)
              read(fh,*,iostat=ios)head,theta,(rot_vec(i),i=1,3)
              lrotate=.true.
              !
          case default
            print*,' ERROR: head not recognised: ',head
            stop
          end select
          !
        enddo
        close(fh)
        print*,' >> ',trim(infile)
        !
      endif
      !
      do js=1,nsolid
        !
        call solidrange(immbody(js),inputcmd='checkdomain')
        !
        if(lscale) then
          call solidresc(immbody(js),scale)
        endif
        !
        if(lrotate) then
          call solidrota(immbody(js),theta,rot_vec)
        endif
        !
        if(lshift) then
          call solidshif(immbody(js),x=xcen(1)-immbody(js)%xcen(1),  &
                                     y=xcen(2)-immbody(js)%xcen(2),  &
                                     z=xcen(3)-immbody(js)%xcen(3))
        endif
        !
        if(ndims==2) then
          ! call solidreduc(immbody(js))
          call solidsilce(immbody(js),zsec=zmin)
        else
          immbody(js)%num_edge=0
        endif
        !
      enddo
      !
      call tecsolid('tecsolid.plt',immbody,dim=3)
      !
    endif
    !
    return
    !
  end subroutine solidgeom
  !+-------------------------------------------------------------------+
  !| The end of the subroutine solidgeom.                              |
  !+-------------------------------------------------------------------+
  !!
  !+-------------------------------------------------------------------+
  !| This subroutine is used to build a immersed body in the grid.     |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 23-Aug-2021: Created by J. Fang @ Appleton                        |
  !+-------------------------------------------------------------------+
  subroutine immsgrid
    !
    ! local data
    real(8) :: time_beg,subtime
    !
    time_beg=ptime() 
    !
    call nodestatecal
    !
    call boundnodecal
    !
    call icellextent
    !
    call icellijkcal
    !
    call cellbilincoef
    !
    call immblocal
    !
    subtime=ptime()-time_beg 
    !
    if(lio) print*,' ** grid in solid calculated'
    !
    if(lio) write(*,'(A,F12.8)')'  ** time cost in processing immersed grid :',subtime
    !
  end subroutine immsgrid
  !+-------------------------------------------------------------------+
  !| The end of the subroutine immsgrid.                               |
  !+-------------------------------------------------------------------+
  !!
  !+-------------------------------------------------------------------+
  !| This subroutine is used to calculate which nodes in/out of solid. |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 04-Jul-2021: Created by J. Fang @ Appleton                        |
  !+-------------------------------------------------------------------+
  subroutine gridinsolid
    !
    use commtype,  only : solid,triangle,sboun,nodcel
    use commvar,   only : immbody,nsolid,immbond,dxyzmax,dxyzmin,      &
                          npdci,npdcj,npdck,ndims,imb_node_have,       &
                          imb_node_need,num_icell_rank,num_ighost_rank
    use commarray, only : x,nodestat,cell
    use commfunc,  only : dis2point,dis2point2,matinv4,matinv,isidenmar
    use commcal,   only : ijkin
    use parallel,  only : ig0,jg0,kg0,pmerg,syncinmg,syncisup,psum,    &
                          syncweig,pgeticell,pcollecicell,pgather
    use readwrite, only : write_sboun
    !
    ! local data
    integer :: i,j,k,jsd,jfc,counter,ninters,n,n1,m,ii,jj,kk,bignum,   &
               jb,kb,jdir,jx,iss,jss,kss,ks1,nc_f,nc_g,nc_b,k1,i0,j0,k0
    integer :: fh,ncou,ke1,icell_counter,ighost_counter
    integer,allocatable :: icell_order(:),ighost_order(:)
    type(solid),pointer :: pso
    logical :: crossface
    integer :: snodes(27,3)
    integer,allocatable :: i_cell(:,:)
    real(8) :: epsilon
    real(8) :: dist,distmin,var1,var2
    real(8) :: xmin(3),xmax(3),xinc(3)
    real(8),allocatable :: Tm1(:,:),Tm2(:,:),Ti1(:,:),Ti2(:,:),        &
                           xvec(:),xcell(:,:),xnorm(:,:)
    !
    type(sboun),allocatable :: bnodes(:)
    logical :: liout,ljout,lkout,lin
    logical,allocatable :: marker(:,:,:)
    real(8) :: time_beg,time_loc,subtime,subtime1,subtime2,subtime3
    !
    time_beg=ptime() 
    !
    if(lio) print*,' ** grid in solid calculating ...'
    !
    epsilon=1.d-10
    !
    ninters=10
    !
    if(npdci==1) then
      iss=0
    else
      iss=1
    endif
    if(npdcj==1) then
      jss=0
    else
      jss=1
    endif
    if(npdck==1 .or. ndims==2) then
      kss=0
    else
      kss=1
    endif
    !
    if(ndims==2) then
      k1=0
      ks1=0
      !
      allocate( Tm1(4,4),Tm2(4,4),Ti1(4,4),Ti2(4,4),xvec(4),          &
                xcell(4,3),xnorm(4,3))
      !
    elseif(ndims==3) then
      k1=-1
      ks1=1
      !
      allocate( Tm1(8,8),Tm2(8,8),Ti1(8,8),Ti2(8,8),xvec(8),          &
                xcell(8,3),xnorm(8,3))
    endif
    !
    ! allocate( icell_order(1:size(immbond)),                            &
    !           ighost_order(1:size(immbond)*8),                         &
    !           num_icell_rank(0:mpirankmax),                            &
    !           num_ighost_rank(0:mpirankmax),                           &
    !           i_cell(size(immbond),3) )
    ! !
    ! num_ighost_rank=0
    ! num_icell_rank=0
    ! icell_order=0
    ! ighost_order=0
    ! !
    ! icell_counter=0
    ! ighost_counter=0
    ! !
    ! if(lio) print*,'    ** searching cells containing image nodes ... '
    ! !
    ! i_cell=0
    ! !
    ! time_loc=ptime()
    ! !
    ! do jb=1,size(immbond)
    !   !
    !   ! search the cell that contains the image node
    !   !
    !   if(ndims==2) then
    !     !
    !     k=0
    !     loopj: do j=1,jm
    !     loopi: do i=1,im
    !       !
    !       !
    !       lin=nodeincell(cell(i,j,k),immbond(jb)%ximag)
    !       !
    !       if(lin) then
    !         !
    !         i_cell(jb,1)=i+ig0
    !         i_cell(jb,2)=j+jg0
    !         i_cell(jb,3)=k+kg0
    !         !
    !         ! print*,' ** icell found for jb=',jb,' rank=',mpirank
    !         !
    !         exit loopj
    !         !
    !       endif
    !       !
    !     enddo loopi
    !     enddo loopj
    !     !
    !     if(cell(i,j,k)%celltype=='i') then
    !       ! it is not a pure fluid or pure solid cell.
    !       ! extent the ximage
    !     endif
    !     !
    !   elseif(ndims==3) then
    !     !
    !     i_cell(jb,:)=pngrid(immbond(jb)%ximag)
    !     !
    !   endif
    !   !
    ! enddo
    ! !
    ! subtime1=subtime1+ptime()-time_loc 
    ! !
    ! time_loc=ptime()
    ! !
    ! ! put icell to immbond
    ! call pcollecicell(i_cell,immbond)
    ! !
    ! subtime2=ptime()-time_loc 
    ! !
    ! ! to check if icell contain a ghost node
    ! do jb=1,size(immbond)
    !   !
    !   !
    !   ! time_loc=ptime()
    !   ! !
    !   ! immbond(jb)%icell=pgeticell(i_cell(jb,:))
    !   ! !
    !   ! subtime2=subtime2+ptime()-time_loc 
    !   !
    !   immbond(jb)%icell_bnode=0
    !   immbond(jb)%icell_ijk=-1
    !   !
    !   time_loc=ptime()
    !   !
    !   ncou=0
    !   !
    !   i=immbond(jb)%icell(1)-ig0
    !   j=immbond(jb)%icell(2)-jg0
    !   k=immbond(jb)%icell(3)-kg0
    !   !
    !   if(i>=1 .and. i<=im .and. j>=1 .and. j<=jm .and. k>=ks1 .and. k<=km) then
    !     !
    !     ! it is not a pure fluid or pure solid cell.
    !     if(cell(i,j,k)%celltype=='i') then
    !       ! if(immbond(jb)%dis2image==0.d0) then
    !       !   print*,mpirank,'|',immbond(jb)%dis2image
    !       !   print*,mpirank,'|',immbond(jb)%ximag
    !       !   print*,mpirank,'|',immbond(jb)%x,immbond(jb)%nodetype
    !       ! endif
    !       ! now to get the intercepting point of the ray and cell
    !       !
    !       xinc=cellintercep(cell(i,j,k),immbond(jb))
    !       !
    !       ! immbond(jb)%ximag=xinc+immbond(jb)%normdir*dxyzmin*0.1d0
    !       ! immbond(jb)%dis2image=immbond(jb)%dis2image+dxyzmin*0.1d0
    !       ! xinc=xinc-immbond(jb)%x
    !       ! xinc=xinc/sqrt(xinc(1)**2+xinc(2)**2)
    !       print*,mpirank,'|',xinc(1:2),':',dis2point(xinc,immbond(jb)%ximag)
    !       !
    !     endif
    !     !
    !   endif
    !   !
    !   ! do kk=k1,0
    !   ! do jj=-1,0
    !   ! do ii=-1,0
    !   !   !
    !   !   ncou=ncou+1
    !   !   !
    !   !   immbond(jb)%icell_ijk(ncou,1)=i
    !   !   immbond(jb)%icell_ijk(ncou,2)=j
    !   !   immbond(jb)%icell_ijk(ncou,3)=k
    !   !   !
    !   !   if(i>=0 .and. i<=im .and. j>=0 .and. j<=jm .and. k>=0 .and. k<=km) then
    !   !     !
    !   !     if(nodestat(i,j,k)>0) then
    !   !       ! icell contain a solid node or a boundary node
    !   !       !
    !   !       ! ! use the boundary node to instead the ghost node
    !   !       ! do kb=1,size(immbond)
    !   !       !   !
    !   !       !   if( immbond(kb)%igh(1)==i+ig0 .and. &
    !   !       !       immbond(kb)%igh(2)==j+jg0 .and. &
    !   !       !       immbond(kb)%igh(3)==k+kg0 ) then
    !   !       !     !
    !   !       !     immbond(jb)%icell_bnode(ncou)=kb
    !   !       !     !
    !   !       !     ! reset ijk that is effective
    !   !       !     immbond(jb)%icell_ijk(ncou,1)=-1
    !   !       !     immbond(jb)%icell_ijk(ncou,2)=-1
    !   !       !     immbond(jb)%icell_ijk(ncou,3)=-1
    !   !       !     !
    !   !       !     ! if(mpirank==0) then
    !   !       !     !   print*,mpirank,'|',immbond(jb)%icell
    !   !       !     !   print*,mpirank,'|',kb
    !   !       !     ! endif
    !   !       !     !
    !   !       !     exit
    !   !       !     !
    !   !       !   endif
    !   !       !   !
    !   !       ! enddo
    !   !       !
    !   !       ! extent ximag
    !   !       !
    !   !       immbond(jb)%ximag=immbond(jb)%ximag
    !   !       !
    !   !     endif
    !   !     !
    !   !   endif
    !   !   !
    !   ! enddo
    !   ! enddo
    !   ! enddo
    !   !
    !   subtime3=subtime3+ptime()-time_loc 
    !   !
    ! enddo
    ! !
    ! if(lio) print*,'    ** supporting cell established'
    ! !
    ! if(lio) write(*,'(A,F12.8)')'    ** time cost in search        :',subtime1
    ! if(lio) write(*,'(A,F12.8)')'    ** time cost in pgeticell     :',subtime2
    ! if(lio) write(*,'(A,F12.8)')'    ** time cost in checking icell:',subtime3
    ! if(lio) write(*,'(A,F12.8)')'    ** total time cost :',subtime1+subtime2+subtime3
    ! !
    ! call mpistop
    ! !
    ! do jb=1,size(immbond)
    !   !
    !   ! determine interpolation coefficient
    !   !
    !   i=immbond(jb)%icell(1)-ig0
    !   j=immbond(jb)%icell(2)-jg0
    !   k=immbond(jb)%icell(3)-kg0
    !   !
    !   if(ndims==2) then
    !     !
    !     !
    !     if(i>=1 .and. i<=im .and. j>=1 .and. j<=jm ) then
    !       ! icell is in the domain
    !       !
    !       do m=1,4
    !         !
    !         if(immbond(jb)%icell_ijk(m,1)>=0) then
    !           i=immbond(jb)%icell_ijk(m,1)
    !           j=immbond(jb)%icell_ijk(m,2)
    !           k=immbond(jb)%icell_ijk(m,3)
    !           !
    !           xcell(m,:)=x(i,j,k,:)
    !           !
    !           Tm1(m,1)=xcell(m,1)*xcell(m,2)
    !           Tm1(m,2)=xcell(m,1)
    !           Tm1(m,3)=xcell(m,2)
    !           Tm1(m,4)=1.d0
    !           !
    !           Tm2(m,1)=xcell(m,1)*xcell(m,2)
    !           Tm2(m,2)=xcell(m,1)
    !           Tm2(m,3)=xcell(m,2)
    !           Tm2(m,4)=1.d0
    !           !
    !         elseif(immbond(jb)%icell_bnode(m)>0) then
    !           kb=immbond(jb)%icell_bnode(m)
    !           !
    !           xcell(m,:)=immbond(kb)%x(:)
    !           xnorm(m,:)=immbond(kb)%normdir(:)
    !           !
    !           Tm1(m,1)=xcell(m,1)*xcell(m,2)
    !           Tm1(m,2)=xcell(m,1)
    !           Tm1(m,3)=xcell(m,2)
    !           Tm1(m,4)=1.d0
    !           !
    !           Tm2(m,1)=xcell(m,1)*xnorm(m,2)+xcell(m,2)*xnorm(m,1)
    !           Tm2(m,2)=xnorm(m,1)
    !           Tm2(m,3)=xnorm(m,2)
    !           Tm2(m,4)=0.d0
    !           !
    !         else
    !           stop ' !! ERROR in determining interpolation coefficient'
    !         endif
    !         !
    !       enddo
    !       !
    !       allocate( immbond(jb)%coef_dirichlet(4),                   &
    !                 immbond(jb)%coef_neumann(4)  )
    !       !
    !       Ti1=matinv(Tm1,4)
    !       Ti2=matinv(Tm2,4)
    !       !
    !       xvec(1)=immbond(jb)%ximag(1)*immbond(jb)%ximag(2)
    !       xvec(2)=immbond(jb)%ximag(1)
    !       xvec(3)=immbond(jb)%ximag(2)
    !       xvec(4)=1.d0
    !       !
    !       do m=1,4
    !         immbond(jb)%coef_dirichlet(m)=dot_product(xvec,Ti1(:,m))
    !         immbond(jb)%coef_neumann(m)  =dot_product(xvec,Ti2(:,m))
    !       enddo
    !       !
    !       ! immbond(jb)%coef_dirichlet=matinv4(Tm1)
    !       ! immbond(jb)%coef_neumann  =matinv4(Tm2)
    !       !
    !       ! Ti1=matmul(Tm1,Ti1)
    !       ! Ti2=matmul(Tm2,Ti2)
    !       ! !
    !       ! write(*,"(I0,A,4(1X,F16.12))")mpirank,'|',Ti2(1,:)
    !       ! write(*,"(I0,A,4(1X,F16.12))")mpirank,'|',Ti2(2,:)
    !       ! write(*,"(I0,A,4(1X,F16.12))")mpirank,'|',Ti2(3,:)
    !       ! write(*,"(I0,A,4(1X,F16.12))")mpirank,'|',Ti2(4,:)
    !       ! print*,'---------------------------------------------------------'
    !       !
    !     endif
    !     !
    !   elseif(ndims==3) then
    !     !
    !     if(i>=1 .and. i<=im .and. j>=1 .and. j<=jm .and. k>=1 .and. k<=km ) then
    !       ! icell is in the domain
    !       !
    !       !
    !       do m=1,8
    !         !
    !         if(immbond(jb)%icell_ijk(m,1)>=0) then
    !           !
    !           i=immbond(jb)%icell_ijk(m,1)
    !           j=immbond(jb)%icell_ijk(m,2)
    !           k=immbond(jb)%icell_ijk(m,3)
    !           !
    !           xcell(m,:)=x(i,j,k,:)
    !           !
    !           Tm1(m,1)=xcell(m,1)*xcell(m,2)*xcell(m,3)
    !           Tm1(m,2)=xcell(m,1)*xcell(m,2)
    !           Tm1(m,3)=xcell(m,1)*xcell(m,3)
    !           Tm1(m,4)=xcell(m,2)*xcell(m,3)
    !           Tm1(m,5)=xcell(m,1)
    !           Tm1(m,6)=xcell(m,2)
    !           Tm1(m,7)=xcell(m,3)
    !           Tm1(m,8)=1.d0
    !           !
    !           Tm2(m,1)=xcell(m,1)*xcell(m,2)*xcell(m,3)
    !           Tm2(m,2)=xcell(m,1)*xcell(m,2)
    !           Tm2(m,3)=xcell(m,1)*xcell(m,3)
    !           Tm2(m,4)=xcell(m,2)*xcell(m,3)
    !           Tm2(m,5)=xcell(m,1)
    !           Tm2(m,6)=xcell(m,2)
    !           Tm2(m,7)=xcell(m,3)
    !           Tm2(m,8)=1.d0
    !           !
    !         elseif(immbond(jb)%icell_bnode(m)>0) then
    !           !
    !           kb=immbond(jb)%icell_bnode(m)
    !           !
    !           xcell(m,:)=immbond(kb)%x(:)
    !           xnorm(m,:)=immbond(kb)%normdir(:)
    !           !
    !           Tm1(m,1)=xcell(m,1)*xcell(m,2)*xcell(m,3)
    !           Tm1(m,2)=xcell(m,1)*xcell(m,2)
    !           Tm1(m,3)=xcell(m,1)*xcell(m,3)
    !           Tm1(m,4)=xcell(m,2)*xcell(m,3)
    !           Tm1(m,5)=xcell(m,1)
    !           Tm1(m,6)=xcell(m,2)
    !           Tm1(m,7)=xcell(m,3)
    !           Tm1(m,8)=1.d0
    !           !
    !           Tm2(m,1)=xnorm(m,1)*xcell(m,2)*xcell(m,3) +  &
    !                    xcell(m,1)*xnorm(m,2)*xcell(m,3) +  &
    !                    xcell(m,1)*xcell(m,2)*xnorm(m,3)
    !           Tm2(m,2)=xnorm(m,1)*xcell(m,2)+xcell(m,1)*xnorm(m,2)
    !           Tm2(m,3)=xnorm(m,1)*xcell(m,3)+xcell(m,1)*xnorm(m,3)
    !           Tm2(m,4)=xnorm(m,2)*xcell(m,3)+xcell(m,2)*xnorm(m,3)
    !           Tm2(m,5)=xnorm(m,1)
    !           Tm2(m,6)=xnorm(m,2)
    !           Tm2(m,7)=xnorm(m,3)
    !           Tm2(m,8)=0.d0
    !           !
    !         else
    !           print*,' ** mpirank=',mpirank
    !           print*,' ** immbond(jb)%icell      :',immbond(jb)%icell
    !           print*,' ** immbond(jb)%icell_ijk  :',immbond(jb)%icell_ijk(m,:)
    !           print*,' ** immbond(jb)%icell_bnode:',immbond(jb)%icell_bnode(m)
    !           stop ' !! ERROR in determining interpolation coefficient'
    !         endif
    !         !
    !       enddo
    !       !
    !       allocate( immbond(jb)%coef_dirichlet(8),                   &
    !                 immbond(jb)%coef_neumann(8),immbond(jb)%invmatrx(8,8)  )
    !       !
    !       Ti1=matinv(Tm1,8)
    !       Ti2=matinv(Tm2,8)
    !       !
    !       xvec(1)=immbond(jb)%ximag(1)*immbond(jb)%ximag(2)*immbond(jb)%ximag(3)
    !       xvec(2)=immbond(jb)%ximag(1)*immbond(jb)%ximag(2)
    !       xvec(3)=immbond(jb)%ximag(1)*immbond(jb)%ximag(3)
    !       xvec(4)=immbond(jb)%ximag(2)*immbond(jb)%ximag(3)
    !       xvec(5)=immbond(jb)%ximag(1)
    !       xvec(6)=immbond(jb)%ximag(2)
    !       xvec(7)=immbond(jb)%ximag(3)
    !       xvec(8)=1.d0
    !       !
    !       do m=1,8
    !         immbond(jb)%coef_dirichlet(m)=dot_product(xvec,Ti1(:,m))
    !         immbond(jb)%coef_neumann(m)  =dot_product(xvec,Ti2(:,m))
    !       enddo
    !       ! immbond(jb)%invmatrx=Ti2
    !       !
    !       ! immbond(jb)%coef_dirichlet=matinv4(Tm1)
    !       ! immbond(jb)%coef_neumann  =matinv4(Tm2)
    !       !
    !       ! Ti1=matmul(Tm1,Ti1)
    !       ! Ti2=matmul(Tm2,Ti2)
    !       ! !
    !       ! if(.not. isidenmar(Ti1,1.d-5)) then
    !       !   write(*,"(8(1X,F10.6))")Ti1(1,:)
    !       !   write(*,"(8(1X,F10.6))")Ti1(2,:)
    !       !   write(*,"(8(1X,F10.6))")Ti1(3,:)
    !       !   write(*,"(8(1X,F10.6))")Ti1(4,:)
    !       !   write(*,"(8(1X,F10.6))")Ti1(5,:)
    !       !   write(*,"(8(1X,F10.6))")Ti1(6,:)
    !       !   write(*,"(8(1X,F10.6))")Ti1(7,:)
    !       !   write(*,"(8(1X,F10.6))")Ti1(8,:)
    !       !   print*,'---------------------------------------------------------'
    !       ! endif
    !       ! if(.not. isidenmar(Ti2,1.d-5)) then
    !       !   write(*,"(8(1X,F10.6))")Ti2(1,:)
    !       !   write(*,"(8(1X,F10.6))")Ti2(2,:)
    !       !   write(*,"(8(1X,F10.6))")Ti2(3,:)
    !       !   write(*,"(8(1X,F10.6))")Ti2(4,:)
    !       !   write(*,"(8(1X,F10.6))")Ti2(5,:)
    !       !   write(*,"(8(1X,F10.6))")Ti2(6,:)
    !       !   write(*,"(8(1X,F10.6))")Ti2(7,:)
    !       !   write(*,"(8(1X,F10.6))")Ti2(8,:)
    !       !   print*,'---------------------------------------------------------'
    !       ! endif
    !       !
    !     endif
    !     !
    !   endif
    !   !
    ! enddo
    ! if(lio) print*,'    ** interpolating coefficients established.'
    ! !
    ! do jb=1,size(immbond)
    !   !
    !   ! setting locality
    !   immbond(jb)%localin=.false.
    !   !
    !   i=immbond(jb)%igh(1)-ig0
    !   j=immbond(jb)%igh(2)-jg0
    !   k=immbond(jb)%igh(3)-kg0
    !   !
    !   if(i>=0 .and. i<=im .and. &
    !      j>=0 .and. j<=jm .and. &
    !      k>=0 .and. k<=km ) then
    !     !
    !     immbond(jb)%localin=.true.
    !     !
    !     ighost_counter=ighost_counter+1
    !     !
    !     ighost_order(ighost_counter)=jb
    !     !
    !     num_ighost_rank(mpirank)=num_ighost_rank(mpirank)+1
    !     !
    !   endif
    !   !
    !   ! checking icell 
    !   i=immbond(jb)%icell(1)-ig0
    !   j=immbond(jb)%icell(2)-jg0
    !   k=immbond(jb)%icell(3)-kg0
    !   !
    !   if(i>=1 .and. i<=im .and. &
    !      j>=1 .and. j<=jm .and. &
    !      k>=ks1 .and. k<=km ) then
    !     !
    !     immbond(jb)%cellin=.true.
    !     !
    !     icell_counter=icell_counter+1
    !     !
    !     icell_order(icell_counter)=jb
    !     !
    !     num_icell_rank(mpirank)=num_icell_rank(mpirank)+1
    !     !
    !   else
    !     immbond(jb)%cellin=.false.
    !   endif
    !   !
    ! enddo
    ! !
    ! call pgather(icell_order(1:icell_counter),imb_node_have)
    ! !
    ! call pgather(ighost_order(1:ighost_counter),imb_node_need)
    ! !
    ! num_icell_rank=psum(num_icell_rank)
    ! !
    ! num_ighost_rank=psum(num_ighost_rank)
    ! !
    ! if(lio) print*,'    ** boundary nodes re-ordered'
    ! !
    ! do k=0,km
    ! do j=0,jm
    ! do i=0,im
    !   rnodestat(i,j,k)=dble(nodestat(i,j,k))
    ! enddo
    ! enddo
    ! enddo
    ! !
    ! call tecbin('testout/tecgrid'//mpirankname//'.plt',           &
    !                                   x(0:im,0:jm,0:km,1),'x',    &
    !                                   x(0:im,0:jm,0:km,2),'y',    &
    !                                   x(0:im,0:jm,0:km,3),'z',    &
    !                           rnodestat(0:im,0:jm,0:km),'ns' )
    ! ! !
    ! ! ! 
    ! ! !
    ! call write_sboun(immbond,'image')
    ! call write_sboun(immbond,'icell')
    ! call write_sboun(immbond,'ghost')
    ! ! !
    ! deallocate(rnodestat,bnodes,marker)
    ! deallocate(Tm1,Tm2,Ti1,Ti2,xvec,xcell,xnorm)
    ! !
    ! subtime=ptime()-time_beg 
    ! !
    ! if(lio) print*,' ** grid in solid calculated'
    ! if(lio) write(*,'(A,F12.8)')'    ** time cost:',subtime
    !
    call mpistop
    !
    return
    !
  end subroutine gridinsolid
  !+-------------------------------------------------------------------+
  !| The end of the subroutine gridinsolid.                            |
  !+-------------------------------------------------------------------+
  !!
  !+-------------------------------------------------------------------+
  !| This subroutine is to calculate nodestat.                         |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 23-Aug-2021: Created by J. Fang @ Appleton                        |
  !+-------------------------------------------------------------------+
  subroutine nodestatecal
    !
    use commtype,  only : solid
    use commvar,   only : immbody,npdci,npdcj,npdck,nsolid,dxyzmin
    use commarray, only : x,nodestat,cell
    use commfunc,  only : dis2point
    !
    ! local data
    integer :: n,ncou,jsd,i,j,k,ii,jj,kk
    type(solid),pointer :: pso
    logical,allocatable :: marker(:,:,:)
    logical :: linsold
    real(8) :: xp(3),pcdir(3),dist
    !
    if(lio) print*,' ** identifying solid nodes ...'
    !
    nodestat=0
    ncou=0
    do jsd=1,nsolid
      !
      pso=>immbody(jsd)
      !
      do k=0,km
      do j=0,jm
      do i=0,im
        !
        xp=x(i,j,k,:)
        !
        pcdir=pso%xcen(:)-xp
        dist=pcdir(1)**2+pcdir(2)**2+pcdir(3)**2
        if(dist>1.d-10) then
          xp=xp+1.d-4*dxyzmin*pcdir/sqrt(dist)
        endif
        !
        if(xp(1)<pso%xmin(1) .or.  xp(1)>pso%xmax(1) .or.    &
           xp(2)<pso%xmin(2) .or.  xp(2)>pso%xmax(2) .or.    &
           xp(3)<pso%xmin(3) .or.  xp(3)>pso%xmax(3) ) then
          !
          nodestat(i,j,k)=0
          ! fluids 
          !
        else
          !
          ! to calculation the intersection between nodes and face
          ! !
          if(ndims==2) then
            linsold=polyhedron_contains_point_2d(pso,xp)
          elseif(ndims==3) then
            linsold=polyhedron_contains_point_3d(pso,xp)
          else
            stop ' !! ERROR1 @ nodestatecal' 
          endif
          !
          if(linsold) then
            ! ths point is in the solid
            nodestat(i,j,k)=5
            !
            ncou=ncou+1
            !
          else
            ! fluids
            nodestat(i,j,k)=0
            !
          endif
          !
        endif
        !
      enddo
      enddo
      enddo
      !
    enddo
    !
    ncou=psum(ncou)
    if(lio) write(*,'(A,I0)')'  ** total number of solide nodes: ',ncou
    !
    call dataswap(nodestat)
    !
    if(lio) print*,' ** state of grid nodes calculated'
    !
    ! search for near-boundary ghost nodes (1-5)
    allocate(marker(0:im,0:jm,0:km))
    n=0
    do while(n<5)
      !
      marker=.false.
      !
      do k=0,km
      do j=0,jm
      do i=0,im
        !
        if(nodestat(i,j,k)==5) then
          ! solid nodes
          !
          do ii=-1,1,2
            if(nodestat(i+ii,j,k)==n) then
              marker(i,j,k)=.true.
              exit
            endif
          enddo
          !
          do jj=-1,1,2
            if(nodestat(i,j+jj,k)==n) then
              marker(i,j,k)=.true.
              exit
            endif
          enddo
          !
          if(ndims==3) then
            do kk=-1,1,2
              if(nodestat(i,j,k+kk)==n) then
                marker(i,j,k)=.true.
                exit
              endif
            enddo
          endif
          !
        endif
        !
      enddo
      enddo
      enddo
      !
      do k=0,km
      do j=0,jm
      do i=0,im
        !
        if(marker(i,j,k)) nodestat(i,j,k)=n+1
        !
      enddo
      enddo
      enddo
      !
      call dataswap(nodestat)
      !
      n=n+1
      !
    enddo
    !
    ! call tecbin('testout/tecgrid'//mpirankname//'.plt',           &
    !                              nodestat(0:im,0:jm,0:km),'ns' )
    !
    ! search for near-boundary force nodes (1-5)
    !
    marker=.false.
    !
    do k=0,km
    do j=0,jm
    do i=0,im
      !
      if(nodestat(i,j,k)==0) then
        ! for fluids nodes
        !
        do ii=-1,1,2
          if(nodestat(i+ii,j,k)>0) then
            marker(i,j,k)=.true.
            exit
          endif
        enddo
        !
        do jj=-1,1,2
          if(nodestat(i,j+jj,k)>0) then
            marker(i,j,k)=.true.
            exit
          endif
        enddo
        !
        if(ndims==3) then
          do kk=-1,1,2
            if(nodestat(i,j,k+kk)>0) then
              marker(i,j,k)=.true.
              exit
            endif
          enddo
        endif
        !
      endif
      !
    enddo
    enddo
    enddo
    !
    do k=0,km
    do j=0,jm
    do i=0,im
      !
      if(marker(i,j,k)) nodestat(i,j,k)=-1
      !
    enddo
    enddo
    enddo
    !
    call dataswap(nodestat)
    !
    ! set bc for nodes, no need to go to the dummy edge or corner
    nodestat(-hm:-1,     -hm:-1,   :)=10
    nodestat(im+1:im+hm, -hm:-1,   :)=10
    nodestat(-hm:-1,    jm+1:jm+hm,:)=10
    nodestat(im+1:im+hm,jm+1:jm+hm,:)=10
    !
    nodestat(-hm:-1,    :,-hm:-1)=10
    nodestat(im+1:im+hm,:,-hm:-1)=10
    nodestat(-hm:-1,    :,km+1:km+hm)=10
    nodestat(im+1:im+hm,:,km+1:km+hm)=10
    !
    nodestat(:,-hm:-1,-hm:-1)=10
    nodestat(:,jm+1:jm+hm,-hm:-1)=10
    nodestat(:,-hm:-1,km+1:km+hm)=10
    nodestat(:,jm+1:jm+hm,km+1:km+hm)=10
    !
    if(npdci==1) then
      nodestat(-hm:-1,:,:)=10
    elseif(npdci==2) then
      nodestat(im+1:im+hm,:,:)=10
    endif
    !
    if(npdcj==1) then
      nodestat(:,-hm:-1,:)=10
    elseif(npdcj==2) then
      nodestat(:,jm+1:jm+hm,:)=10
    endif
    !
    if(npdck==1) then
      nodestat(:,:,-hm:-1)=10
    elseif(npdck==2) then
      nodestat(:,:,km+1:km+hm)=10
    endif
    !
    !
    if(lio) print*,' ** near-boundary ghost nodes identified'
    !
    ! set cell state.
    if(ndims==2) then
      !
      k=0
      do j=1,jm
      do i=1,im
        !
        if( nodestat(i-1,j-1,k)<=0 .and. &
            nodestat(i,  j-1,k)<=0 .and. &
            nodestat(i,  j,  k)<=0 .and. &
            nodestat(i-1,j,  k)<=0 ) then
          !
          cell(i,j,k)%celltype='f'
          ! fluids
          !
        elseif( nodestat(i-1,j-1,k)>0 .and. &
                nodestat(i,  j-1,k)>0 .and. &
                nodestat(i,  j,  k)>0 .and. &
                nodestat(i-1,j,  k)>0 ) then
          !
          cell(i,j,k)%celltype='s'
          ! solid
          !
        else
          !
          cell(i,j,k)%celltype='i'
          ! interface
          !
        endif
        !
      enddo
      enddo
      !
    elseif(ndims==3) then
      !
      do k=1,km
      do j=1,jm
      do i=1,im
        !
        if( nodestat(i-1,j-1,k)  <=0 .and. &
            nodestat(i,  j-1,k)  <=0 .and. &
            nodestat(i,  j,  k)  <=0 .and. &
            nodestat(i-1,j,  k)  <=0 .and. &
            nodestat(i-1,j-1,k-1)<=0 .and. &
            nodestat(i,  j-1,k-1)<=0 .and. &
            nodestat(i,  j,  k-1)<=0 .and. &
            nodestat(i-1,j,  k-1)<=0 ) then
          !
          cell(i,j,k)%celltype='f'
          ! fluids
          !
        elseif( nodestat(i-1,j-1,k)  >0 .and. &
                nodestat(i,  j-1,k)  >0 .and. &
                nodestat(i,  j,  k)  >0 .and. &
                nodestat(i-1,j,  k)  >0 .and. &
                nodestat(i-1,j-1,k-1)>0 .and. &
                nodestat(i,  j-1,k-1)>0 .and. &
                nodestat(i,  j,  k-1)>0 .and. &
                nodestat(i-1,j,  k-1)>0) then
          !
          cell(i,j,k)%celltype='s'
          ! solid
          !
        else
          !
          cell(i,j,k)%celltype='i'
          ! interface
          !
        endif
        !
      enddo
      enddo
      enddo
      !
    else 
      stop ' !! ERROR in ndims @ solidgeom'
    endif
    !
    !
    if(lio) print*,' ** cell state set'
    !
  end subroutine nodestatecal
  !+-------------------------------------------------------------------+
  !| The end of the subroutine nodestatecal.                           |
  !+-------------------------------------------------------------------+
  !!
  !+-------------------------------------------------------------------+
  !| This subroutine is to establish the boundary node.                |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 05-Jul-2021: Created by J. Fang @ Appleton                        |
  !+-------------------------------------------------------------------+
  subroutine boundnodecal
    !
    use commvar,   only : nsolid,immbody,immbond,dxyzmin
    use commtype,  only : solid,sboun
    use commarray, only : nodestat,x
    use parallel,  only : ig0,jg0,kg0,pmerg
    use commfunc,  only : dis2point
    use readwrite, only : write_sboun
    !
    ! local data
    integer :: counter,nc_f,nc_b,nc_g,bignum,jsd,iss,jss,kss,i,j,k
    type(solid),pointer :: pso
    type(sboun),allocatable :: bnodes(:)
    !
    bignum=im*jm*(km+1)
    allocate(bnodes(bignum))
    !
    if(npdci==1) then
      iss=0
    else
      iss=1
    endif
    if(npdcj==1) then
      jss=0
    else
      jss=1
    endif
    if(npdck==1 .or. ndims==2) then
      kss=0
    else
      kss=1
    endif
    !
    ! to get the nodes and distance of inner solid nodes to boundary 
    counter=0
    nc_f=0
    nc_b=0
    nc_g=0
    !
    do jsd=1,nsolid
      !
      pso=>immbody(jsd)
      !
      do k=kss,km
      do j=jss,jm
      do i=iss,im
        !
        if(nodestat(i,j,k)>0 .and. nodestat(i,j,k)<5) then
          ! ghost point
          !
          counter=counter+1
          !
          call polyhedron_bound_search(pso,x(i,j,k,:),bnodes(counter),dir='+')
          !
          bnodes(counter)%igh(1)=i+ig0
          bnodes(counter)%igh(2)=j+jg0
          bnodes(counter)%igh(3)=k+kg0
          !
          bnodes(counter)%ximag(:)=2.d0*bnodes(counter)%x(:)-x(i,j,k,:)
          !
          bnodes(counter)%dis2image=dis2point(bnodes(counter)%x,bnodes(counter)%ximag)
          bnodes(counter)%dis2ghost=dis2point(bnodes(counter)%x,x(i,j,k,:))
          !
          if(dis2point(bnodes(counter)%x,x(i,j,k,:))<dxyzmin*0.01d0) then
            ! the distance of node to wall is less thant 1/100 of 
            ! the min grid spacing
            bnodes(counter)%nodetype='b'
            nc_b=nc_b+1
            !
            ! boundary marker
          else
            bnodes(counter)%nodetype='g'
            nc_g=nc_g+1
          endif
          !
          !
        ! elseif(nodestat(i,j,k)==-1) then
        !   ! force point
        !   !
        !   counter=counter+1
        !   !
        !   call polyhedron_bound_search(pso,x(i,j,k,:),bnodes(counter),dir='-')
        !   !
        !   bnodes(counter)%igh(1)=i+ig0
        !   bnodes(counter)%igh(2)=j+jg0
        !   bnodes(counter)%igh(3)=k+kg0
        !   !
        !   bnodes(counter)%ximag(:)=2.d0*x(i,j,k,:)-bnodes(counter)%x(:)
        !   !
        !   bnodes(counter)%nodetype='f'
        !   !
        !   nc_f=nc_f+1
        endif
        !
      enddo
      enddo
      enddo
      !
    enddo
    !
    call pmerg(var=bnodes,nvar=counter,vmerg=immbond)
    !
    if(lio) print*,' ** ghost, boundary and image nodes collected'
    !
    nc_b=psum(nc_b)
    nc_f=psum(nc_f)
    nc_g=psum(nc_g)
    if(lio) then
      write(*,'(A,I0)')'     ** number of boundary nodes: ',size(immbond)
      write(*,'(A,I0)')'        **    ghost nodes: ',nc_g
      write(*,'(A,I0)')'        **  forcing nodes: ',nc_f
      write(*,'(A,I0)')'        ** boundary nodes: ',nc_b
    endif
    !
    return
    !
  end subroutine boundnodecal
  !+-------------------------------------------------------------------+
  !| The end of the subroutine boundnodecal.                           |
  !+-------------------------------------------------------------------+
  !!
  !+-------------------------------------------------------------------+
  !| This subroutine is to build the cell that contains image nodes.   |
  !| for those cell in the sold, extent the image node.                |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 23-Aug-2021: Created by J. Fang @ Appleton                        |
  !+-------------------------------------------------------------------+
  subroutine icellextent
    !
    use commtype,  only : sboun
    use commvar,   only : immbond,dxyzmin
    use commarray, only : x,nodestat,cell
    use commfunc,  only : dis2point
    use parallel,  only : pcollecicell,por
    use readwrite, only : write_sboun
    !
    ! local data
    integer,allocatable :: i_cell(:,:)
    integer :: jb,i,j,k,ncou
    real(8) :: time_beg,subtime
    real(8) :: deltamv,xintcell(3)
    logical,allocatable :: icell_marker(:)
    type(sboun),pointer :: pbon
    !
    ! start ... 
    time_beg=ptime() 
    !
    allocate( i_cell(size(immbond),3),icell_marker(size(immbond)) )
    !
    if(lio) print*,' ** searching cells containing image nodes ... '
    !
    icell_marker=.false.
    !
    i_cell=0
    !
    ncou=0
    do jb=1,size(immbond)
      ! go through all boundary nodes
      !
      pbon=>immbond(jb)
      ! search the cell that contains the image node
      !
      if(ndims==2) then
        !
        i_cell(jb,:)=pngrid2d(pbon%ximag)
        !
      elseif(ndims==3) then
        !
        i_cell(jb,:)=pngrid3d(pbon%ximag)
        !
      endif
      !
      i=i_cell(jb,1)-ig0
      j=i_cell(jb,2)-jg0
      k=i_cell(jb,3)-kg0
      !
      if(ijkcellin(i,j,k)) then
        !
        if(cell(i,j,k)%celltype=='i' .or. cell(i,j,k)%celltype=='s') then
          ! it is not a pure fluid.
          ! extent the ximage
          ncou=ncou+1
          !
          icell_marker(jb)=.true.
          !
        else
        endif
        !
      endif
      !
    enddo
    !
    ncou=psum(ncou)
    icell_marker=por(icell_marker)
    !
    if(lio) print*,' ** number of cut cell:',ncou
    !
    ! go through the marked boundary nodes
    !
    do while(ncou>0)
      !
      !
      do jb=1,size(immbond)
        !
        deltamv=0.d0
        !
        pbon=>immbond(jb)
        !
        if(icell_marker(jb)) then
          !
          i=i_cell(jb,1)-ig0
          j=i_cell(jb,2)-jg0
          k=i_cell(jb,3)-kg0
          !
          if(ijkcellin(i,j,k)) then
            ! calculate the intersection the boundary extension 
            ! ray and the cell
            xintcell=cellintersec(cell(i,j,k),pbon)
            !
            deltamv=dis2point(xintcell,pbon%ximag)
            !
          endif
          !
          ! extend the ximag
          deltamv=pmax(deltamv)
          !
          xintcell=pbon%ximag+pbon%normdir*(deltamv+0.01d0*dxyzmin)
          pbon%ximag=xintcell
          !
          pbon%dis2image=dis2point(pbon%x,pbon%ximag)
          !
        endif
        !
      enddo
      !
      ! rebuilt icell
      ncou=0
      do jb=1,size(immbond)
        !
        pbon=>immbond(jb)
        !
        if(icell_marker(jb)) then
          !
          icell_marker(jb)=.false.
          !
          if(ndims==2) then
            !
            i_cell(jb,:)=pngrid2d(pbon%ximag)
            !
          elseif(ndims==3) then
            !
            i_cell(jb,:)=pngrid3d(pbon%ximag)
            !
          endif
          !
          i=i_cell(jb,1)-ig0
          j=i_cell(jb,2)-jg0
          k=i_cell(jb,3)-kg0
          !
          if(ijkcellin(i,j,k)) then
            !
            !
            if(cell(i,j,k)%celltype=='i' .or. cell(i,j,k)%celltype=='s') then
              ! it is not a pure fluid or pure solid cell.
              ! extent the ximage
              ncou=ncou+1
              !
              icell_marker(jb)=.true.
              !
            else
            endif
            !
          endif
          !
        endif
        !
      enddo
      !
      ncou=psum(ncou)
      icell_marker=por(icell_marker)
      !
      if(lio) print*,' ** number of cut cell:',ncou
      !
    enddo
    !
    !
    ! put icell to immbond
    call pcollecicell(i_cell,immbond)
    !
    call write_sboun(immbond,'image')
    !
    call write_sboun(immbond,'ghost')
    !
    call write_sboun(immbond,'icell')
    !
    subtime=ptime()-time_beg 
    !
    if(lio) write(*,'(A,F12.8)')'  ** time cost in searching icell :',subtime
    !
  end subroutine icellextent
  !+-------------------------------------------------------------------+
  !| The end of the subroutine icellextent.                            |
  !+-------------------------------------------------------------------+
  !!
  !+-------------------------------------------------------------------+
  !| This subroutine is to calculate the ijk of icell.                 |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 24-Aug-2021: Created by J. Fang @ Appleton                        |
  !+-------------------------------------------------------------------+
  subroutine icellijkcal
    !
    use commtype,  only : sboun
    use commvar,   only : immbond,ndims
    use commarray, only : nodestat,cell
    !
    integer :: jb,ncou,k1,i,j,k,ic,jc,kc,ii,jj,kk
    type(sboun),pointer :: pbon
    real(8) :: epsilon=1.d-12
    !
    if(ndims==2) then
      k1=0
    elseif(ndims==3) then
      k1=-1
    endif
    !
    do jb=1,size(immbond)
      !
      pbon=>immbond(jb)
      !
      pbon%icell_bnode=0
      pbon%icell_ijk=-1
      !
      ncou=0
      do kk=k1,0
      do jj=-1,0
      do ii=-1,0
        !
        ncou=ncou+1
        !
        ic=pbon%icell(1)-ig0
        jc=pbon%icell(2)-jg0
        kc=pbon%icell(3)-kg0
        !
        i=ic+ii
        j=jc+jj
        k=kc+kk
        !
        pbon%icell_ijk(ncou,1)=i
        pbon%icell_ijk(ncou,2)=j
        pbon%icell_ijk(ncou,3)=k
        !
        if(ijkin(i,j,k)) then
          !
          if(nodestat(i,j,k)>0) then
            ! icell contain a solid node or a boundary node
            ! should not happen now
            print*,mpirank,'| cell:',ic,jc,kc,cell(ic,jc,kc)%celltype
            print*,mpirank,'| node:',i,j,k
            print*,' !! ERROR: icell contain a solid node or a boundary node'
            stop 
          endif
          !
          if(pbon%dis2image<epsilon) then
            print*,mpirank,'|',i,j,k
            print*,' !! WARNING the distance between image node to boundary is too small',pbon%dis2image
            stop 
          endif
          !
        endif
        !
      enddo
      enddo
      enddo
      !
    enddo
    !
    if(lio) write(*,'(A)')'  ** icell ijk built.'
    !
  end subroutine icellijkcal
  !+-------------------------------------------------------------------+
  !| The end of the subroutine icellijkcal.                            |
  !+-------------------------------------------------------------------+
  !!
  !+-------------------------------------------------------------------+
  !| This subroutine is to calcualte coefficient for bilinear          |
  !| interpolation                                                     |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 05-Jul-2021: Created by J. Fang @ Appleton                        |
  !+-------------------------------------------------------------------+
  subroutine cellbilincoef
    !
    use commtype,  only : sboun
    use commvar,   only : immbond,ndims
    use commarray, only : x
    use commfunc,  only : matinv
    !
    integer :: jb,i,j,k,ic,jc,kc,m,ncord
    type(sboun),pointer :: pbon
    real(8),allocatable :: xcell(:,:),Tm1(:,:),Ti1(:,:),xvec(:)
    !
    ncord=2**ndims
    allocate(Tm1(1:ncord,1:ncord),xcell(1:ncord,1:3),xvec(1:ncord))
    !
    do jb=1,size(immbond)
      !
      pbon=>immbond(jb)
      !
      ic=pbon%icell(1)-ig0
      jc=pbon%icell(2)-jg0
      kc=pbon%icell(3)-kg0
      !
      if(ijkcellin(ic,jc,kc)) then
        !
        if(ndims==2) then
          !
          do m=1,4
            !
            i=pbon%icell_ijk(m,1)
            j=pbon%icell_ijk(m,2)
            k=pbon%icell_ijk(m,3)
            !
            xcell(m,:)=x(i,j,k,:)
            !
            Tm1(m,1)=xcell(m,1)*xcell(m,2)
            Tm1(m,2)=xcell(m,1)
            Tm1(m,3)=xcell(m,2)
            Tm1(m,4)=1.d0
            !
          enddo
          !
          allocate( pbon%coef_dirichlet(4) )
          !
          Ti1=matinv(Tm1,4)
          !
          xvec(1)=pbon%ximag(1)*pbon%ximag(2)
          xvec(2)=pbon%ximag(1)
          xvec(3)=pbon%ximag(2)
          xvec(4)=1.d0
          !
          do m=1,4
            pbon%coef_dirichlet(m)=dot_product(xvec,Ti1(:,m))
          enddo
          !
        elseif(ndims==3) then
          !
          do m=1,8
            !
            i=pbon%icell_ijk(m,1)
            j=pbon%icell_ijk(m,2)
            k=pbon%icell_ijk(m,3)
            !
            xcell(m,:)=x(i,j,k,:)
            !
            Tm1(m,1)=xcell(m,1)*xcell(m,2)*xcell(m,3)
            Tm1(m,2)=xcell(m,1)*xcell(m,2)
            Tm1(m,3)=xcell(m,1)*xcell(m,3)
            Tm1(m,4)=xcell(m,2)*xcell(m,3)
            Tm1(m,5)=xcell(m,1)
            Tm1(m,6)=xcell(m,2)
            Tm1(m,7)=xcell(m,3)
            Tm1(m,8)=1.d0
            !
          enddo
          !
          allocate( pbon%coef_dirichlet(8) )
          !
          Ti1=matinv(Tm1,8)
          !
          xvec(1)=immbond(jb)%ximag(1)*immbond(jb)%ximag(2)*immbond(jb)%ximag(3)
          xvec(2)=immbond(jb)%ximag(1)*immbond(jb)%ximag(2)
          xvec(3)=immbond(jb)%ximag(1)*immbond(jb)%ximag(3)
          xvec(4)=immbond(jb)%ximag(2)*immbond(jb)%ximag(3)
          xvec(5)=immbond(jb)%ximag(1)
          xvec(6)=immbond(jb)%ximag(2)
          xvec(7)=immbond(jb)%ximag(3)
          xvec(8)=1.d0
          !
          do m=1,8
            pbon%coef_dirichlet(m)=dot_product(xvec,Ti1(:,m))
          enddo
          !
        endif
        !
        ! Ti1=matmul(Tm1,Ti1)
        ! !
        ! write(*,"(I0,A,8(1X,F8.4))")mpirank,'|',Ti1(1,:)
        ! write(*,"(I0,A,8(1X,F8.4))")mpirank,'|',Ti1(2,:)
        ! write(*,"(I0,A,8(1X,F8.4))")mpirank,'|',Ti1(3,:)
        ! write(*,"(I0,A,8(1X,F8.4))")mpirank,'|',Ti1(4,:)
        ! write(*,"(I0,A,8(1X,F8.4))")mpirank,'|',Ti1(5,:)
        ! write(*,"(I0,A,8(1X,F8.4))")mpirank,'|',Ti1(6,:)
        ! write(*,"(I0,A,8(1X,F8.4))")mpirank,'|',Ti1(7,:)
        ! write(*,"(I0,A,8(1X,F8.4))")mpirank,'|',Ti1(8,:)
        ! print*,'---------------------------------------------------------'
        !
      endif
      !
    enddo
    !
    deallocate(Tm1,xcell)
    !
    if(lio) write(*,'(A)')'  ** cell interpolation coefficients calculated'
    !
  end subroutine cellbilincoef
  !+-------------------------------------------------------------------+
  !| The end of the subroutine cellbilincoef.                          |
  !+-------------------------------------------------------------------+
  !!
  !+-------------------------------------------------------------------+
  !| This subroutine is to collect ib information.                     |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 25-Agu-2021: Created by J. Fang @ Appleton                        |
  !+-------------------------------------------------------------------+
  subroutine immblocal
    !!
    use commtype,  only : sboun
    use commvar,   only : immbond,ndims,imb_node_have,imb_node_need,   &
                          num_icell_rank,num_ighost_rank
    use parallel,  only : pgather
    !
    integer :: jb,i,j,k,icell_counter,ighost_counter
    integer,allocatable :: ighost_order(:),icell_order(:)
    type(sboun),pointer :: pbon
    !
    allocate( icell_order(1:size(immbond)),                            &
              ighost_order(1:size(immbond)*8),                         &
              num_icell_rank(0:mpirankmax),                            &
              num_ighost_rank(0:mpirankmax) )
    icell_counter=0
    ighost_counter=0
    !
    num_ighost_rank=0
    num_icell_rank=0
    icell_order=0
    ighost_order=0
    
    do jb=1,size(immbond)
      !
      pbon=>immbond(jb)
      !
      ! setting locality
      pbon%localin=.false.
      !
      i=pbon%igh(1)-ig0
      j=pbon%igh(2)-jg0
      k=pbon%igh(3)-kg0
      !
      if( ijkin(i,j,k) ) then
        !
        pbon%localin=.true.
        !
        ighost_counter=ighost_counter+1
        !
        ighost_order(ighost_counter)=jb
        !
        num_ighost_rank(mpirank)=num_ighost_rank(mpirank)+1
        !
      endif
      !
      ! checking icell 
      pbon%cellin=.false.
      !
      i=pbon%icell(1)-ig0
      j=pbon%icell(2)-jg0
      k=pbon%icell(3)-kg0
      !
      if( ijkcellin(i,j,k) ) then
        !
        pbon%cellin=.true.
        !
        icell_counter=icell_counter+1
        !
        icell_order(icell_counter)=jb
        !
        num_icell_rank(mpirank)=num_icell_rank(mpirank)+1
        !
      endif
      !
    enddo
    !
    call pgather(icell_order(1:icell_counter),imb_node_have)
    !
    call pgather(ighost_order(1:ighost_counter),imb_node_need)
    !
    num_icell_rank=psum(num_icell_rank)
    !
    num_ighost_rank=psum(num_ighost_rank)
    !
    if(lio) print*,' ** boundary nodes re-ordered'
    !
  end subroutine immblocal
  !+-------------------------------------------------------------------+
  !| The end of the subroutine immblocal.                              |
  !+-------------------------------------------------------------------+
  !!
  !!
  !+-------------------------------------------------------------------+
  !| This function is to judge if a point in a cell.                   |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 05-Jul-2021: Created by J. Fang @ Appleton                        |
  !+-------------------------------------------------------------------+
  function nodeincell(acell,p) result(lin)
    !
    use commtype, only : nodcel,solid
    use commvar,  only : ndims
    use commfunc,  only : cross_product
    !
    ! arguments
    type(nodcel),intent(in) :: acell
    type(solid) :: scell
    real(8),intent(in) :: p(3)
    logical :: lin
    !
    real(8) :: norm1(3),norm2(3),var1(3)
    integer :: n
    !
    if(p(1)<acell%xmin(1) .or. p(1)>acell%xmax(1) .or. &
       p(2)<acell%xmin(2) .or. p(2)>acell%xmax(2) .or. &
       p(3)<acell%xmin(3) .or. p(3)>acell%xmax(3) ) then
      !
      lin=.false.
      !
    else
      !
      if(ndims==2) then
        !
        ! checking in the point p is in any of the two triangles.
        lin=pointintriangle(acell%x(1,:),acell%x(2,:),acell%x(3,:),p)
        !
        if(lin) then
          return
        else
          lin=pointintriangle(acell%x(1,:),acell%x(4,:),acell%x(3,:),p)
          return
        endif
        !
      elseif(ndims==3) then
        !
        scell%num_face=12
        !
        call scell%alloface()
        !
        scell%xcen(:)=0.125d0*( acell%x(1,:)+acell%x(2,:) + & 
                                acell%x(3,:)+acell%x(4,:) + &
                                acell%x(5,:)+acell%x(6,:) + &
                                acell%x(7,:)+acell%x(8,:)   )
        !
        ! k-1 face
        scell%face(1)%a=acell%x(1,:)
        scell%face(1)%b=acell%x(4,:)
        scell%face(1)%c=acell%x(3,:)
        !
        scell%face(2)%a=acell%x(3,:)
        scell%face(2)%b=acell%x(2,:)
        scell%face(2)%c=acell%x(1,:)
        !
        ! k face
        scell%face(3)%a=acell%x(6,:)
        scell%face(3)%b=acell%x(7,:)
        scell%face(3)%c=acell%x(8,:)
        !
        scell%face(4)%a=acell%x(8,:)
        scell%face(4)%b=acell%x(5,:)
        scell%face(4)%c=acell%x(6,:)
        !
        ! j-1 face
        scell%face(5)%a=acell%x(5,:)
        scell%face(5)%b=acell%x(1,:)
        scell%face(5)%c=acell%x(2,:)
        !
        scell%face(6)%a=acell%x(2,:)
        scell%face(6)%b=acell%x(6,:)
        scell%face(6)%c=acell%x(5,:)
        !
        ! j face
        scell%face(7)%a=acell%x(7,:)
        scell%face(7)%b=acell%x(3,:)
        scell%face(7)%c=acell%x(4,:)
        !
        scell%face(8)%a=acell%x(4,:)
        scell%face(8)%b=acell%x(8,:)
        scell%face(8)%c=acell%x(7,:)
        !
        ! i-1 face
        scell%face(9)%a=acell%x(4,:)
        scell%face(9)%b=acell%x(1,:)
        scell%face(9)%c=acell%x(5,:)
        !
        scell%face(10)%a=acell%x(5,:)
        scell%face(10)%b=acell%x(8,:)
        scell%face(10)%c=acell%x(4,:)
        !
        ! i face
        scell%face(11)%a=acell%x(2,:)
        scell%face(11)%b=acell%x(3,:)
        scell%face(11)%c=acell%x(7,:)
        !
        scell%face(12)%a=acell%x(7,:)
        scell%face(12)%b=acell%x(6,:)
        scell%face(12)%c=acell%x(2,:)
        !
        do n=1,12
          norm1=cross_product(scell%face(n)%b-scell%face(n)%a,  &
                              scell%face(n)%c-scell%face(n)%a  )
          norm2=scell%face(n)%b-scell%xcen
          !
          if(dot_product(norm1,norm2)<0.d0) then
            var1=scell%face(n)%a
            scell%face(n)%a=scell%face(n)%c
            scell%face(n)%c=var1
            print*,norm1,'-',norm2
          endif
          !
        enddo
        !
        lin=polyhedron_contains_point_3d(scell,p)
        !
        ! if( p(1)>=acell%xmin(1) .and. p(1)<=acell%xmax(1) .and. &
        !     p(2)>=acell%xmin(2) .and. p(2)<=acell%xmax(2) .and. &
        !     p(3)>=acell%xmin(3) .and. p(3)<=acell%xmax(3) ) then
        !   !
        !   lin=.true.
        ! else
        !   lin=.false.
        ! endif
        !
      else
        stop ' !! dimension not set yet @ nodeincell!!'
      endif
      !
    endif
    !
    return
    !
  end function nodeincell
  !+-------------------------------------------------------------------+
  !| The end of the function nodeincell.                               |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to improve quality of a solid,                 |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 05-Jul-2021: Created by J. Fang @ Appleton                        |
  !+-------------------------------------------------------------------+
  subroutine solidimpro(asolid)
    !
    use commtype,  only : solid,triangle
    use commvar,   only : xmax,xmin,ymax,ymin,zmax,zmin
    use commfunc,  only : dis2point,cross_product
    !
    ! arguments
    type(solid),intent(inout) :: asolid
    !
    ! local data
    integer :: jf,n,counter,largn
    real(8) :: var1
    real(8),dimension(3) :: norm1,a,b,c
    real(8),allocatable ::  psav(:,:)
    real(8) :: epsilon
    !
    epsilon=2.5d-6
    largn=10000
    !
    allocate(psav(3,largn))
    !
    psav(:,1)=asolid%face(1)%a
    psav(:,2)=asolid%face(1)%b
    psav(:,3)=asolid%face(1)%c
    !
    counter=3
    !
    do jf=2,asolid%num_face
      !
      a=asolid%face(jf)%a
      b=asolid%face(jf)%b
      c=asolid%face(jf)%c
      !
      do n=1,counter
        var1=dis2point(a,psav(:,n))
        if(var1>0 .and. var1<=epsilon) then
          print*,a,'-',psav(:,n),':',var1
          asolid%face(jf)%a=psav(:,n)
          exit
        endif
      enddo
      !
      if(n>counter) then
        ! no record in the psav
        counter=counter+1
        psav(:,counter)=a
      endif
      !
      do n=1,counter
        var1=dis2point(b,psav(:,n))
        if(var1>0 .and. var1<=epsilon) then
          print*,b,'-',psav(:,n),':',var1
          asolid%face(jf)%b=psav(:,n)
          exit
        endif
      enddo
      !
      if(n>counter) then
        ! no record in the psav
        counter=counter+1
        psav(:,counter)=b
      endif
      !
      do n=1,counter
        var1=dis2point(c,psav(:,n))
        if(var1>0 .and. var1<=epsilon) then
          print*,c,'-',psav(:,n),':',var1
          asolid%face(jf)%c=psav(:,n)
          exit
        endif
      enddo
      !
      if(n>counter) then
        ! no record in the psav
        counter=counter+1
        psav(:,counter)=c
      endif
      !
      !
    enddo
    !
    do jf=1,asolid%num_face
      norm1=cross_product(asolid%face(jf)%b-asolid%face(jf)%a,         &
                          asolid%face(jf)%c-asolid%face(jf)%a )
      !
      var1=sqrt(norm1(1)**2+norm1(2)**2+norm1(3)**2)
      norm1=norm1/var1
      !
      if(dot_product(norm1,asolid%face(jf)%normdir)>0.d0) then
        asolid%face(jf)%normdir=norm1
      elseif(dot_product(norm1,asolid%face(jf)%normdir)<0.d0) then
        asolid%face(jf)%normdir=-norm1
      else
        print*,'face number=',jf
        print*,asolid%face(jf)%a
        print*,asolid%face(jf)%b
        print*,asolid%face(jf)%c
        print*,norm1,'-',asolid%face(jf)%normdir
        stop ' !! ERROR 1 @ solidimpro'
      endif
      !
    enddo
    !
  end subroutine solidimpro
  !+-------------------------------------------------------------------+
  !| The end of the subroutine solidimpro.                             |
  !+-------------------------------------------------------------------+
  !

  !+-------------------------------------------------------------------+
  !| This subroutine is used to cut a 2d slice from a 3D geometry.     |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 07-Jul-2021: Created by J. Fang @ Appleton                        |
  !+-------------------------------------------------------------------+
  subroutine solidsilce(asolid,zsec)
    !
    use commtype,  only : solid,triangle,lsegment
    use commvar,   only : xmax,xmin,ymax,ymin,zmax,zmin
    use commfunc,  only : areatriangle,cross_product,dis2point
    !
    ! arguments
    type(solid),intent(inout),target :: asolid
    real(8),intent(in) :: zsec
    !
    ! local data
    integer :: i,jf,je,nedge,nemax
    type(lsegment),allocatable :: edge_temp(:)
    type(triangle),pointer :: aface
    real(8) :: dz2,dz1,epsilon,tange(2),norma(2),n2cen(2)
    logical :: lbpoint
    !
    epsilon=1.d-12
    !
    nemax=asolid%num_face
    !
    allocate(edge_temp(nemax))
    !
    nedge=0
    do jf=1,asolid%num_face
      !
      aface=>asolid%face(jf)
      !
      lbpoint=.false.
      !
      if( (aface%a(3)>=zsec .and. aface%b(3)<=zsec) .or. &
          (aface%a(3)<=zsec .and. aface%b(3)>=zsec)  ) then
        !
        dz2=aface%b(3)-aface%a(3)
        dz1=zsec-aface%a(3)
        !
        if(abs(dz2)<epsilon) then
          !
          nedge=nedge+1
          !
          edge_temp(nedge)%a=aface%a(1:2)
          edge_temp(nedge)%b=aface%b(1:2)
          !
          if(dis2point(edge_temp(nedge)%a,edge_temp(nedge)%b)>=epsilon) then
            cycle
          else
            nedge=nedge-1
          endif
          !
        else
          nedge=nedge+1
          !
          edge_temp(nedge)%a=dz1/dz2*(aface%b(1:2)-aface%a(1:2))+aface%a(1:2)
          lbpoint=.true.
          !
        endif
        !
      endif
      !
      if( (aface%b(3)>=zsec .and. aface%c(3)<=zsec) .or. &
          (aface%b(3)<=zsec .and. aface%c(3)>=zsec)  ) then
        !
        dz2=aface%c(3)-aface%b(3)
        dz1=zsec-aface%b(3)
        !
        if(abs(dz2)<epsilon) then
          !
          if(.not. lbpoint) nedge=nedge+1
          !
          edge_temp(nedge)%a=aface%b(1:2)
          edge_temp(nedge)%b=aface%c(1:2)
          !
          if(dis2point(edge_temp(nedge)%a,edge_temp(nedge)%b)>=epsilon) then
            cycle
          else
            nedge=nedge-1
          endif
          !
        else
          !
          if(lbpoint) then
            !
            edge_temp(nedge)%b=dz1/dz2*(aface%c(1:2)-aface%b(1:2))+aface%b(1:2)
            !
            if(dis2point(edge_temp(nedge)%a,edge_temp(nedge)%b)>=epsilon) then
              cycle
            else
              nedge=nedge-1
            endif
            !
          else
            nedge=nedge+1
            edge_temp(nedge)%a=dz1/dz2*(aface%c(1:2)-aface%b(1:2))+aface%b(1:2)
            lbpoint=.true.
          endif
          !
        endif
        !
      endif
      !
      if( (aface%c(3)>=zsec .and. aface%a(3)<=zsec) .or. &
          (aface%c(3)<=zsec .and. aface%a(3)>=zsec)  ) then
        !
        dz2=aface%a(3)-aface%c(3)
        dz1=zsec-aface%c(3)
        !
        if(abs(dz2)<epsilon) then
          !
          if(.not. lbpoint) nedge=nedge+1
          !
          edge_temp(nedge)%a=aface%c(1:2)
          edge_temp(nedge)%b=aface%a(1:2)
          !
          !
          if(dis2point(edge_temp(nedge)%a,edge_temp(nedge)%b)>=epsilon) then
            cycle
          else
            nedge=nedge-1
          endif
          !
        else
          !
          if(lbpoint) then
            !
            edge_temp(nedge)%b=dz1/dz2*(aface%a(1:2)-aface%c(1:2))+aface%c(1:2)
            !
            if(dis2point(edge_temp(nedge)%a,edge_temp(nedge)%b)>=epsilon) then
              cycle
            else
              nedge=nedge-1
            endif
            !
          else
            stop ' !! ERROR @ solidsilce'
          endif
          !
        endif
        !
      endif
      !
    enddo
    !
    if(nedge==0) then
      print*,' !! ERROR, failed to cut a slice from the solid'
      print*,' ** zsec=',zsec
      print*,' ** solid extent:',asolid%xmin,'~',asolid%xmax
    endif
    !
    asolid%num_edge=nedge
    call asolid%alloedge()
    !
    asolid%edge(1:nedge)=edge_temp(1:nedge)
    !
    asolid%xcen=0.d0
    do je=1,asolid%num_edge
      asolid%xcen(1:2)=asolid%xcen(1:2)+0.5d0*(asolid%edge(je)%a+asolid%edge(je)%b)
    enddo
    asolid%xcen=asolid%xcen/dble(asolid%num_edge)
    asolid%xcen(3)=zsec
    !
    print*,' ** center of the solid: ',asolid%xcen
    !
    do je=1,asolid%num_edge
      !
      tange=asolid%edge(je)%b-asolid%edge(je)%a
      norma(1)=-tange(2)/sqrt(tange(1)**2+tange(2)**2)
      norma(2)= tange(1)/sqrt(tange(1)**2+tange(2)**2)
      !
      if(sqrt(tange(1)**2+tange(2)**2)<1.d-16) then
        print*,' !! ERROR: the two ends of the edge is two close'
        print*,je,asolid%edge(je)%a,':',asolid%edge(je)%b
      endif
      !
      n2cen=0.5d0*(asolid%edge(je)%a+asolid%edge(je)%b)-asolid%xcen(1:2)
      !
      if(dot_product(norma,n2cen)>=0.d0) then
        !
        asolid%edge(je)%normdir=norma
        !
      else
        !
        asolid%edge(je)%normdir=-norma
        !
      endif
      !
    enddo
    !
    if(lio) print*,' ** the body is slided at z=',zsec
    !
  end subroutine solidsilce
  !+-------------------------------------------------------------------+
  !| The end of the subroutine solidsilce.                             |
  !+-------------------------------------------------------------------+
  !!
  !+-------------------------------------------------------------------+
  !| This subroutine is used to reduce a 3D polyhedron to a 2D polygon |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 07-Jul-2021: Created by J. Fang @ Appleton                        |
  !+-------------------------------------------------------------------+
  subroutine solidreduc(asolid)
    !
    use commtype,  only : solid,triangle,lsegment
    use commvar,   only : xmax,xmin,ymax,ymin,zmax,zmin
    use commfunc,  only : areatriangle,cross_product
    !
    ! arguments
    type(solid),intent(inout) :: asolid
    !
    ! local data
    integer :: i,jf,je,nedge,nemax
    type(lsegment),allocatable :: edge_temp(:)
    real(8) :: dz,epsilon,nz
    logical :: lbpoint
    !
    epsilon=1.d-12
    !
    nemax=asolid%num_face*3
    !
    allocate(edge_temp(nemax))
    !
    nedge=0
    do jf=1,asolid%num_face
      !
      dz=abs(asolid%face(jf)%a(3)-asolid%xmin(3))
      nz=abs(abs(asolid%face(jf)%normdir(3))-1.d0)
      !
      if(dz<=epsilon .and. nz>epsilon) then
        nedge=nedge+1
        edge_temp(nedge)%a(1:2)=asolid%face(jf)%a(1:2)
        lbpoint=.true.
      else
        lbpoint=.false.
      endif
      !
      dz=abs(asolid%face(jf)%b(3)-asolid%xmin(3))
      nz=abs(abs(asolid%face(jf)%normdir(3))-1.d0)
      !
      if(dz<=epsilon .and. nz>epsilon) then
        !
        if(lbpoint) then
          edge_temp(nedge)%b(1:2)=asolid%face(jf)%b(1:2)
          edge_temp(nedge)%normdir(1:2)=asolid%face(jf)%normdir(1:2)
          !
          lbpoint=.false.
        else
          nedge=nedge+1
          edge_temp(nedge)%a(1:2)=asolid%face(jf)%b(1:2)
          lbpoint=.true.
        endif
        !
      endif
      !
      if(.not. lbpoint) cycle
      !
      dz=abs(asolid%face(jf)%c(3)-asolid%xmin(3))
      nz=abs(abs(asolid%face(jf)%normdir(3))-1.d0)
      !
      if(dz<=epsilon .and. nz>epsilon) then
        !
        if(lbpoint) then
          edge_temp(nedge)%b(1:2)=asolid%face(jf)%c(1:2)
          edge_temp(nedge)%normdir(1:2)=asolid%face(jf)%normdir(1:2)
          lbpoint=.false.
          !
        else
          stop ' !! ERROR @ solidreduc'
          !
        endif
        !
      elseif(lbpoint) then
        nedge=nedge-1
      endif
      
    enddo
    !
    asolid%num_edge=nedge
    call asolid%alloedge()
    !
    asolid%edge(1:nedge)=edge_temp(1:nedge)
    !
    ! do jf=1,nedge
    !   print*,jf,'|',asolid%edge(jf)%normdir
    ! enddo
    !
    ! open(18,file='test.dat')
    ! do je=1,nedge
    ! enddo
    ! close(18)
    !
  end subroutine solidreduc
  !+-------------------------------------------------------------------+
  !| The end of the subroutine solidimpro.                             |
  !+-------------------------------------------------------------------+
  !+-------------------------------------------------------------------+
  !| This subroutine is used to obtain the range of the solid          |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 03-Jul-2021: Created by J. Fang @ Appleton                        |
  !+-------------------------------------------------------------------+
  subroutine solidrange(asolid,inputcmd)
    !
    use commtype,  only : solid,triangle
    use commvar,   only : xmax,xmin,ymax,ymin,zmax,zmin
    use commfunc,  only : areatriangle,cross_product
    !
    ! arguments
    type(solid),intent(inout) :: asolid
    character(len=*),intent(in),optional :: inputcmd
    !
    ! local data
    integer :: i,jf
    real(8) :: var1
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
      !
      asolid%face(jf)%area=areatriangle(asolid%face(jf)%a,           &
                                        asolid%face(jf)%b,           &
                                        asolid%face(jf)%c )
      !
      var1= asolid%face(jf)%normdir(1)**2 +  &
            asolid%face(jf)%normdir(2)**2 +  &
            asolid%face(jf)%normdir(3)**2
      !
      asolid%face(jf)%normdir=asolid%face(jf)%normdir/sqrt(var1)
      !
      if(asolid%face(jf)%area<1.d-16) then
        print*,' !! WARNING @ solidrange'
        print*,jf
        print*,asolid%face(jf)%a,'-',asolid%face(jf)%b,'-',asolid%face(jf)%c
        print*,asolid%face(jf)%area
        stop
      endif
      !
      asolid%face(jf)%cen(:)=num1d3*(asolid%face(jf)%a(:)+           &
                                     asolid%face(jf)%b(:)+           &
                                     asolid%face(jf)%c(:) )
    enddo
    !
    asolid%xcen(:)=0.5d0*(asolid%xmin(:)+asolid%xmax(:))
    !
    write(*,'(A,A)')      '  solid: ',asolid%name
    write(*,'(2(A,E15.7E3))')'     x: ',asolid%xmin(1),'~',asolid%xmax(1)
    write(*,'(2(A,E15.7E3))')'     y: ',asolid%xmin(2),'~',asolid%xmax(2)
    write(*,'(2(A,E15.7E3))')'     z: ',asolid%xmin(3),'~',asolid%xmax(3)
    write(*,'(3(A,E15.7E3))')'  size: ',asolid%xmax(1)-asolid%xmin(1),':', &
                                        asolid%xmax(2)-asolid%xmin(2),':', &
                                        asolid%xmax(3)-asolid%xmin(3)
    write(*,'(A,3(E15.7E3))')'  cent: ',asolid%xcen(:)
    !
    if(trim(inputcmd)=='checkdomain') then
      !
      if(asolid%xmin(1)<xmin) then
        write(*,'(A,2(E15.7E3))')'  !! WARNING: solid outof domain'
        write(*,'(2(A,E15.7E3))')'  !! solid xmin: ',asolid%xmin(1),   &
                                    ' domain xmin: ',xmin
      endif
      if(asolid%xmax(1)>xmax) then
        write(*,'(A,2(E15.7E3))')'  !! WARNING: solid outof domain'
        write(*,'(2(A,E15.7E3))')'  !! solid xmax: ',asolid%xmax(1),   &
                                    ' domain xmax: ',xmax
      endif
      if(asolid%xmin(2)<ymin) then
        write(*,'(A,2(E15.7E3))')'  !! WARNING: solid outof domain'
        write(*,'(2(A,E15.7E3))')'  !! solid ymin: ',asolid%xmin(2),   &
                                    ' domain ymin: ',ymin
      endif
      if(asolid%xmax(2)>ymax) then
        write(*,'(A,2(E15.7E3))')'  !! WARNING: solid outof domain'
        write(*,'(2(A,E15.7E3))')'  !! solid ymax: ',asolid%xmax(2),  &
                                    ' domain ymax: ',ymax
      endif
      !
      if(ndims==3) then
        if(asolid%xmin(3)<zmin) then
          write(*,'(A,2(E15.7E3))')'  !! WARNING: solid outof domain'
          write(*,'(2(A,E15.7E3))')'  !! solid zmin: ',asolid%xmin(3), &
                                      ' domain zmin: ',zmin
        endif
        !
        if(asolid%xmax(3)>zmax) then
          write(*,'(A,2(E15.7E3))')'  !! WARNING: solid outof domain'
          write(*,'(2(A,E15.7E3))')'  !! solid zmax: ',asolid%xmax(3), &
                                      ' domain zmax: ',zmax
        endif
      endif
      !
    endif
    !
    write(*,'(2X,62A)')('-',i=1,62)
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
    asolid%xmin=asolid%xmin*scale
    asolid%xmax=asolid%xmax*scale
    asolid%xcen=asolid%xcen*scale
    !
    write(*,'(3A,E15.7E3)')'  ** solid ',trim(asolid%name),            &
                                                 ' is rescaled by',scale
    !
    call solidrange(asolid)
    !
  end subroutine solidresc
  !+-------------------------------------------------------------------+
  !| The end of the subroutine solidresc.                              |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to shift size of the solid.               |
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
    call solidrange(asolid)
    !
  end subroutine solidshif
  !+-------------------------------------------------------------------+
  !| The end of the subroutine solidshif.                              |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to rotate size of the solid.              |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 22-Aug-2021: Created by J. Fang @ Appleton                        |
  !+-------------------------------------------------------------------+
  !| ref: https://zh.wikipedia.org/wiki/%E6%97%8B%E8%BD%AC%E7%9F%A9%E9%98%B5
  !+-------------------------------------------------------------------+
  subroutine solidrota(asolid,theta,vec)
    !
    use commtype,  only : solid,triangle
    !
    ! arguments
    type(solid),intent(inout) :: asolid
    real(8),intent(in) :: theta,vec(3)
    !
    ! local data
    integer :: i,jf
    real(8) :: rotmat(3,3),a,vec_xyz(3,3),vec_abc(3,3)
    !
    a=theta/180.d0*pi
    !
    rotmat(1,1)=(1.d0-cos(a))*vec(1)*vec(1)+cos(a)
    rotmat(1,2)=(1.d0-cos(a))*vec(1)*vec(2)-sin(a)*vec(3)
    rotmat(1,3)=(1.d0-cos(a))*vec(1)*vec(3)+sin(a)*vec(2)
    !
    rotmat(2,1)=(1.d0-cos(a))*vec(1)*vec(2)+sin(a)*vec(3)
    rotmat(2,2)=(1.d0-cos(a))*vec(2)*vec(2)+cos(a)
    rotmat(2,3)=(1.d0-cos(a))*vec(2)*vec(3)-sin(a)*vec(1)
    !
    rotmat(3,1)=(1.d0-cos(a))*vec(1)*vec(3)-sin(a)*vec(2)
    rotmat(3,2)=(1.d0-cos(a))*vec(2)*vec(3)+sin(a)*vec(1)
    rotmat(3,3)=(1.d0-cos(a))*vec(3)*vec(3)+cos(a)
    !
    ! print*,rotmat(1,:)
    ! print*,rotmat(2,:)
    ! print*,rotmat(3,:)
    !
    do jf=1,asolid%num_face
      !
      vec_xyz(1,:)=asolid%face(jf)%a
      vec_xyz(2,:)=asolid%face(jf)%b
      vec_xyz(3,:)=asolid%face(jf)%c
      !
      do i=1,3
        vec_abc(1,i)=dot_product(rotmat(i,:),vec_xyz(1,:))
        vec_abc(2,i)=dot_product(rotmat(i,:),vec_xyz(2,:))
        vec_abc(3,i)=dot_product(rotmat(i,:),vec_xyz(3,:))
      enddo
      !
      asolid%face(jf)%a=vec_abc(1,:)
      asolid%face(jf)%b=vec_abc(2,:)
      asolid%face(jf)%c=vec_abc(3,:)
      !
      vec_xyz(1,:)=asolid%face(jf)%normdir
      !
      asolid%face(jf)%normdir(1)=dot_product(rotmat(1,:),vec_xyz(1,:))
      asolid%face(jf)%normdir(2)=dot_product(rotmat(2,:),vec_xyz(1,:))
      asolid%face(jf)%normdir(3)=dot_product(rotmat(3,:),vec_xyz(1,:))
      !
    enddo
    !
    write(*,'(3(A),F10.6,A,3(F10.6))')'  ** solid ',trim(asolid%name), &
                                 ' roates ',theta,'deg along a vec:',vec
    !
    call solidrange(asolid)
    !
  end subroutine solidrota
  !+-------------------------------------------------------------------+
  !| The end of the subroutine solidrota.                              |
  !+-------------------------------------------------------------------+
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! This subtoutine is used to calculate geometrical transform matrix
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Writen by Fang Jian, 2008-09-22.
  ! Change name by Fang Jian, 2008-09-22.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine gridgeom
    !
    use commvar,   only : ia,ja,ka,hm,npdci,npdcj,npdck,               &
                          xmax,xmin,ymax,ymin,zmax,zmin,voldom,difschm,&
                          dxyzmax,dxyzmin
    use commarray, only : x,jacob,dxi,cell,dgrid,dis2wall
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
    real(8) :: var1,var2,var3
    !
    allocate( dx(-hm:im+hm,-hm:jm+hm,-hm:km+hm,1:3,1:3) )
    !
    call gridsendrecv
    !
    ! call xyzbc
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
    dxyzmax=0.d0
    dxyzmin=1.d10
    !
    if(ndims==1) then
      !
      j=0
      k=0
      !
      voldom=0.d0
      do i=1,im
        allocate(cell(i,j,k)%x(2,3))
        cell(i,j,k)%x(1,:)=x(i-1,j,k,:)
        cell(i,j,k)%x(2,:)=x(i,j,k,:)
        !
        cell(i,j,k)%xmin(1)=min(x(i-1,j,k,1),x(i,j,k,1))
        cell(i,j,k)%xmin(2)=min(x(i-1,j,k,1),x(i,j,k,2))
        cell(i,j,k)%xmin(3)=min(x(i-1,j,k,1),x(i,j,k,3))
        !
        cell(i,j,k)%xmax(1)=max(x(i-1,j,k,1),x(i,j,k,1))
        cell(i,j,k)%xmax(2)=max(x(i-1,j,k,1),x(i,j,k,2))
        cell(i,j,k)%xmax(3)=max(x(i-1,j,k,1),x(i,j,k,3))
        !
        cell(i,j,k)%vol=abs(x(i,j,k,1)-x(i-1,j,k,1))
        voldom=voldom+cell(i,j,k)%vol
      enddo
      voldom=psum(voldom)
      !
      do i=0,im
        dxyzmax=max(dxyzmax,abs(dx(i,j,k,1,1)))
        dxyzmin=min(dxyzmin,abs(dx(i,j,k,1,1)))
        !
        dgrid(i,j,k,1)=abs(dx(i,j,k,1,1))
        dgrid(i,j,k,2)=0.d0
        dgrid(i,j,k,3)=0.d0
        !
      enddo
      dxyzmax=pmax(dxyzmax)
      dxyzmin=pmin(dxyzmin)
      !
    elseif(ndims==2) then
      k=0
      !
      voldom=0.d0
      do j=1,jm
      do i=1,im
        allocate(cell(i,j,k)%x(4,3))
        cell(i,j,k)%x(1,:)=x(i-1,j-1,k,:)
        cell(i,j,k)%x(2,:)=x(i,  j-1,k,:)
        cell(i,j,k)%x(3,:)=x(i,  j,  k,:)
        cell(i,j,k)%x(4,:)=x(i-1,j,  k,:)
        !
        cell(i,j,k)%xmin(1)=min(x(i-1,j-1,k,1),x(i,  j-1,k,1),         &
                                x(i,  j,  k,1),x(i-1,j,  k,1))
        cell(i,j,k)%xmin(2)=min(x(i-1,j-1,k,2),x(i,  j-1,k,2),         &
                                x(i,  j,  k,2),x(i-1,j,  k,2))
        cell(i,j,k)%xmin(3)=min(x(i-1,j-1,k,3),x(i,  j-1,k,3),         &
                                x(i,  j,  k,3),x(i-1,j,  k,3))
        !
        cell(i,j,k)%xmax(1)=max(x(i-1,j-1,k,1),x(i,  j-1,k,1),         &
                                x(i,  j,  k,1),x(i-1,j,  k,1))
        cell(i,j,k)%xmax(2)=max(x(i-1,j-1,k,2),x(i,  j-1,k,2),         &
                                x(i,  j,  k,2),x(i-1,j,  k,2))
        cell(i,j,k)%xmax(3)=max(x(i-1,j-1,k,3),x(i,  j-1,k,3),         &
                                x(i,  j,  k,3),x(i-1,j,  k,3))
        !
        cell(i,j,k)%vol=arquad( x(i-1,j-1,k,:),   x(i,j-1,k,:),        &
                                x(i,j,k,:),       x(i-1,j,k,:) )
        voldom=voldom+cell(i,j,k)%vol
      enddo
      enddo
      voldom=psum(voldom)
      !
      do j=0,jm
      do i=0,im
        var1=sqrt(dx(i,j,k,1,1)**2+dx(i,j,k,2,1)**2)
        var2=sqrt(dx(i,j,k,1,2)**2+dx(i,j,k,2,2)**2)
        !
        dgrid(i,j,k,1)=var1
        dgrid(i,j,k,2)=var2
        dgrid(i,j,k,3)=0.d0
        !
        dxyzmax=max(dxyzmax,var1,var2)
        dxyzmin=min(dxyzmin,var1,var2)
      enddo
      enddo
      dxyzmax=pmax(dxyzmax)
      dxyzmin=pmin(dxyzmin)
      !
    elseif(ndims==3) then
      voldom=0.d0
      do k=1,km
      do j=1,jm
      do i=1,im
        allocate(cell(i,j,k)%x(8,3))
        cell(i,j,k)%x(1,:)=x(i-1,j-1,k-1,:)
        cell(i,j,k)%x(2,:)=x(i,  j-1,k-1,:)
        cell(i,j,k)%x(3,:)=x(i,  j,  k-1,:)
        cell(i,j,k)%x(4,:)=x(i-1,j,  k-1,:)
        cell(i,j,k)%x(5,:)=x(i-1,j-1,k,:)
        cell(i,j,k)%x(6,:)=x(i,  j-1,k,:)
        cell(i,j,k)%x(7,:)=x(i,  j,  k,:)
        cell(i,j,k)%x(8,:)=x(i-1,j,  k,:)
        !
        cell(i,j,k)%xmin(1)=min(x(i-1,j-1,k,1),  x(i,  j-1,k,1),       &
                                x(i,  j,  k,1),  x(i-1,j,  k,1),       &
                                x(i-1,j-1,k-1,1),x(i,  j-1,k-1,1),     &
                                x(i,  j,  k-1,1),x(i-1,j,  k-1,1))
        cell(i,j,k)%xmin(2)=min(x(i-1,j-1,k,2),  x(i,  j-1,k,2),       &
                                x(i,  j,  k,2),  x(i-1,j,  k,2),       &
                                x(i-1,j-1,k-1,2),x(i,  j-1,k-1,2),     &
                                x(i,  j,  k-1,2),x(i-1,j,  k-1,2))
        cell(i,j,k)%xmin(3)=min(x(i-1,j-1,k,3),  x(i,  j-1,k,3),       &
                                x(i,  j,  k,3),  x(i-1,j,  k,3),      &
                                x(i-1,j-1,k-1,3),x(i,  j-1,k-1,3),     &
                                x(i,  j,  k-1,3),x(i-1,j,  k-1,3))
        !
        cell(i,j,k)%xmax(1)=max(x(i-1,j-1,k,1),  x(i,  j-1,k,1),       &
                                x(i,  j,  k,1),  x(i-1,j,  k,1),       &
                                x(i-1,j-1,k-1,1),x(i,  j-1,k-1,1),     &
                                x(i,  j,  k-1,1),x(i-1,j,  k-1,1))
        cell(i,j,k)%xmax(2)=max(x(i-1,j-1,k,2),  x(i,  j-1,k,2),       &
                                x(i,  j,  k,2),  x(i-1,j,  k,2),       &
                                x(i-1,j-1,k-1,2),x(i,  j-1,k-1,2),     &
                                x(i,  j,  k-1,2),x(i-1,j,  k-1,2))
        cell(i,j,k)%xmax(3)=max(x(i-1,j-1,k,3),  x(i,  j-1,k,3),       &
                                x(i,  j,  k,3),  x(i-1,j,  k,3),      &
                                x(i-1,j-1,k-1,3),x(i,  j-1,k-1,3),     &
                                x(i,  j,  k-1,3),x(i-1,j,  k-1,3))
        !
        cell(i,j,k)%vol=volhex( x(i-1,j-1,k-1,:), x(i,j-1,k-1,:),      &
                                x(i,j-1,k,:),     x(i-1,j-1,k,:),      &
                                x(i-1,j,k-1,:),   x(i,j,k-1,:)  ,      &
                                x(i,j,k,:),       x(i-1,j,k,:) )
        voldom=voldom+cell(i,j,k)%vol
      enddo
      enddo
      enddo
      voldom=psum(voldom)
      !
      do k=0,km
      do j=0,jm
      do i=0,im
        var1=sqrt(dx(i,j,k,1,1)**2+dx(i,j,k,2,1)**2+dx(i,j,k,3,1)**2)
        var2=sqrt(dx(i,j,k,1,2)**2+dx(i,j,k,2,2)**2+dx(i,j,k,3,2)**2)
        var3=sqrt(dx(i,j,k,1,3)**2+dx(i,j,k,2,3)**2+dx(i,j,k,3,3)**2)
        !
        dgrid(i,j,k,1)=var1
        dgrid(i,j,k,2)=var2
        dgrid(i,j,k,3)=var3
        !
        dxyzmax=max(dxyzmax,var1,var2,var3)
        dxyzmin=min(dxyzmin,var1,var2,var3)
      enddo
      enddo
      enddo
      dxyzmax=pmax(dxyzmax)
      dxyzmin=pmin(dxyzmin)
      !
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
    do k=0,km
    do j=0,jm
    do i=0,im
      ! only for channel flow.
      dis2wall(i,j,k)=min(x(i,j,k,2),2.d0-x(i,j,k,2))
    enddo
    enddo
    enddo
    !
    if(lio) then
      !
      write(*,'(2X,62A)')('-',i=1,62)
      write(*,'(2X,A)')'                    *** Grids Information *** '
      write(*,'(2X,62A)')('-',i=1,62)
      write(*,'(3X,62A)')'       xmin      xmax      ymin      ymax      zmin      zmax'
      write(*,"(4X,6(F10.3))")xmin,xmax,ymin,ymax,zmin,zmax
      write(*,'(2X,62A)')('-',i=1,62)
      write(*,'(7X,2(A,E20.7E3))')'  grid spacing:',dxyzmin,' ~',dxyzmax
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
  end subroutine gridgeom
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! End of the subroutine gridgeom.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !+-------------------------------------------------------------------+
  !| This function is to find the intercept between a ray and a cell.  |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 23-08-2021  | Created by J. Fang @ Warrington                     |
  !+-------------------------------------------------------------------+
  function cellintersec(acell,abond,debug) result(xincep)
    !
    use commtype, only : nodcel,sboun,triangle
    use commfunc, only : dis2point,cross_product
    !
    ! arguments
    type(nodcel),intent(in) :: acell
    type(sboun),intent(in) :: abond
    logical,intent(in),optional :: debug
    real(8) :: xincep(3)
    !
    ! local data
    integer :: n
    real(8) :: x1(4,3),x2(4,3),p(3),dir(3),xint(3),dist,distmax
    type(triangle) :: atri(12)
    logical :: lintercept
    !
    if(ndims==2) then
      !
      p=abond%ximag
      dir=abond%normdir
      !
      ! edge 1
      x1(1,:)=acell%x(1,:); x2(1,:)=acell%x(2,:)
      x1(2,:)=acell%x(2,:); x2(2,:)=acell%x(3,:)
      x1(3,:)=acell%x(3,:); x2(3,:)=acell%x(4,:)
      x1(4,:)=acell%x(4,:); x2(4,:)=acell%x(1,:)
      !
      distmax=-1.d10
      do n=1,4
        !
        call ray2segment(x1(n,:),x2(n,:),p,dir,xint,lintercept)
        !
        if(lintercept) then
          !
          dist=dis2point(xint,p)
          !
          if(dist>distmax) then
            distmax=dist
            xincep=xint
          endif
          !
        endif
        !
      enddo
      !
    elseif(ndims==3) then
      !
      p=abond%ximag
      dir=abond%normdir
      !
      ! establish a hexahedron
      ! k-1 face
      atri(1)%a=acell%x(1,:)
      atri(1)%b=acell%x(4,:)
      atri(1)%c=acell%x(3,:)
      !
      atri(2)%a=acell%x(3,:)
      atri(2)%b=acell%x(2,:)
      atri(2)%c=acell%x(1,:)
      !
      ! k face
      atri(3)%a=acell%x(6,:)
      atri(3)%b=acell%x(7,:)
      atri(3)%c=acell%x(8,:)
      !
      atri(4)%a=acell%x(8,:)
      atri(4)%b=acell%x(5,:)
      atri(4)%c=acell%x(6,:)
      !
      ! j-1 face
      atri(5)%a=acell%x(5,:)
      atri(5)%b=acell%x(1,:)
      atri(5)%c=acell%x(2,:)
      !
      atri(6)%a=acell%x(2,:)
      atri(6)%b=acell%x(6,:)
      atri(6)%c=acell%x(5,:)
      !
      ! j face
      atri(7)%a=acell%x(7,:)
      atri(7)%b=acell%x(3,:)
      atri(7)%c=acell%x(4,:)
      !
      atri(8)%a=acell%x(4,:)
      atri(8)%b=acell%x(8,:)
      atri(8)%c=acell%x(7,:)
      !
      ! i-1 face
      atri(9)%a=acell%x(4,:)
      atri(9)%b=acell%x(1,:)
      atri(9)%c=acell%x(5,:)
      !
      atri(10)%a=acell%x(5,:)
      atri(10)%b=acell%x(8,:)
      atri(10)%c=acell%x(4,:)
      !
      ! i face
      atri(11)%a=acell%x(2,:)
      atri(11)%b=acell%x(3,:)
      atri(11)%c=acell%x(7,:)
      !
      atri(12)%a=acell%x(7,:)
      atri(12)%b=acell%x(6,:)
      atri(12)%c=acell%x(2,:)
      !
      do n=1,12
        atri(n)%normdir=cross_product(atri(n)%b-atri(n)%a,  &
                                      atri(n)%c-atri(n)%a  )
      enddo
      !
      distmax=-1.d10
      do n=1,12
        call ray2triangle(atri(n),p,dir,xint,lintercept)
        !
        if(lintercept) then
          !
          dist=dis2point(xint,p)
          !
          if(dist>distmax) then
            distmax=dist
            xincep=xint
          endif
          !
        endif
        !
      enddo
      !
    else
      print*,' !! ERROR 1 dimension @ cellintersec'
      stop
    endif
    !
    return
    !
  end function cellintersec
  !+-------------------------------------------------------------------+
  !| The end of the function cellintersec.                             |
  !+-------------------------------------------------------------------+
  !!
  !+-------------------------------------------------------------------+
  !| This subroutine is to find the intersection between a ray and a   |
  !| segment                                                           |
  !+-------------------------------------------------------------------+
  ! https://en.wikipedia.org/wiki/Line%E2%80%93line_intersection#Given_two_points_on_each_line_segment
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 24-08-2021  | Created by J. Fang @ Warrington                     |
  !+-------------------------------------------------------------------+
  subroutine ray2segment(x1,x2,p,dir,xint,lint)
    !
    use commfunc, only : dis2point2
    ! arguments
    real(8),intent(in)  :: x1(3),x2(3),p(3),dir(3)
    real(8),intent(out) :: xint(3)
    logical :: lint
    !
    ! local data
    real(8) :: t,u,var1,dis12,dis1p,disp2
    if(ndims==2) then
      !
      var1=-(x1(1)-x2(1))*dir(2)+(x1(2)-x2(2))*dir(1)
      !
      if(abs(var1)<1.d-16) then
        ! The segment is parallel to the ray, check if p is on x1-x2
        !
        lint=.false.
        xint(1)=0.d0
        xint(2)=0.d0
        xint(3)=0.d0
        !
      else
        t=(-(x1(1)-p(1)) *dir(2)+(x1(2)-p(2))*dir(1))/var1
        u=( (x2(1)-x1(1))*(x1(2)-p(2)) - (x2(2)-x1(2))*(x1(1)-p(1)))/var1
        !
        if(u>=0.d0 .and. t>=0.d0 .and. t<=1.d0 ) then
          lint=.true.
          !
          xint(1)=p(1)+dir(1)*u
          xint(2)=p(2)+dir(2)*u
          xint(3)=0.d0
        else
          lint=.false.
          xint(1)=0.d0
          xint(2)=0.d0
          xint(3)=0.d0
        endif
        !
      endif
      !
    else
      stop ' !! ERROR in dimension @ ray2segment'
    endif
    !
    return
    !
  end subroutine ray2segment
  !+-------------------------------------------------------------------+
  !| The end of the subroutine ray2segment.                            |
  !+-------------------------------------------------------------------+
  !!
  !+-------------------------------------------------------------------+
  !| This function is to find the point in 3d grid.                    |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 19-08-2021  | Created by J. Fang @ Warrington                     |
  !+-------------------------------------------------------------------+
  function pngrid2d(p) result(ijk)
    !
    use commvar,   only : im,jm,km
    use commarray, only : x,cell
    !
    ! arguments
    real(8),intent(in) :: p(3)
    integer :: ijk(3)
    !
    integer :: i,j,k
    logical :: lin
    !
    ijk=0
    !
    k=0
    loopj: do j=1,jm
    loopi: do i=1,im
      !
      if(nodeincell(cell(i,j,k),p)) then
        !
        ijk(1)=i+ig0
        ijk(2)=j+jg0
        ijk(3)=k+kg0
        !
        exit loopj
        !
      endif
      !
    enddo loopi
    enddo loopj
    !
  end function pngrid2d
  !+-------------------------------------------------------------------+
  !| The end of the function pngrid2d.                                 |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  ! !| This function is to find the point in 3d grid.                    |
  ! !+-------------------------------------------------------------------+
  ! !| CHANGE RECORD                                                     |
  ! !| -------------                                                     |
  ! !| 19-08-2021  | Created by J. Fang @ Warrington                     |
  ! !+-------------------------------------------------------------------+
  ! function pngrid3d(p) result(ijk)
  !   !
  !   use commvar,   only : im,jm,km
  !   use commarray, only : x,cell
  !   !
  !   ! arguments
  !   real(8),intent(in) :: p(3)
  !   integer :: ijk(3)
  !   !
  !   ! local data
  !   integer :: nlevel,ilevel,idomai,inxt,ncou,nco2
  !   integer :: i,j,k,n,m,m20,is,ie,js,je,ks,ke
  !   integer :: idomai_save(1024),idomai_save2(1024)
  !   real(8) :: xmin(3),xmax(3)
  !   !
  !   logical,save :: lfirstcal=.true.
  !   type :: subdomaim
  !     integer :: i0,j0,k0
  !     integer :: im,jm,km
  !     integer :: next(8)
  !     real(8) :: xmin(3),xmax(3)
  !     type(subdomaim),allocatable :: subdom(:,:,:)
  !   end type subdomaim
  !   type :: level
  !     integer :: numdom
  !     type(subdomaim),allocatable :: domain(:)
  !   end type level
  !   !
  !   type(level),allocatable,save :: dolevl(:)
  !   !
  !   if(lfirstcal) then
  !     !
  !     nlevel=1
  !     do while(2**nlevel<max(im,jm,km))
  !       nlevel=nlevel+1
  !     enddo
  !     ! write(*,'(3(A,I0))')'     ** nlevel= ',nlevel,         &
  !     !                          ',  max nodes= ',2**nlevel, &
  !     !                         '/',max(im,jm,km)
  !     !
  !     allocate(dolevl(1:nlevel))
  !     !
  !     ! the master level
  !     dolevl(1)%numdom=1
  !     allocate(dolevl(1)%domain(dolevl(1)%numdom))
  !     !
  !     do ilevel=2,nlevel
  !       dolevl(ilevel)%numdom=dolevl(ilevel-1)%numdom*8
  !       allocate(dolevl(ilevel)%domain(dolevl(ilevel)%numdom))
  !     enddo
  !     !
  !     dolevl(1)%domain(1)%i0=0
  !     dolevl(1)%domain(1)%j0=0
  !     dolevl(1)%domain(1)%k0=0
  !     dolevl(1)%domain(1)%im=im
  !     dolevl(1)%domain(1)%jm=jm
  !     dolevl(1)%domain(1)%km=km
  !     !
  !     do ilevel=1,nlevel-1
  !       !
  !       ncou=0
  !       do n=1,dolevl(ilevel)%numdom
  !         !
  !         allocate(dolevl(ilevel)%domain(n)%subdom(2,2,2))
  !         !
  !         m20=(dolevl(ilevel)%domain(n)%im-dolevl(ilevel)%domain(n)%i0)/2 + &
  !              dolevl(ilevel)%domain(n)%i0 
  !         dolevl(ilevel)%domain(n)%subdom(1,:,:)%i0=dolevl(ilevel)%domain(n)%i0
  !         dolevl(ilevel)%domain(n)%subdom(1,:,:)%im=m20
  !         dolevl(ilevel)%domain(n)%subdom(2,:,:)%i0=m20
  !         dolevl(ilevel)%domain(n)%subdom(2,:,:)%im=dolevl(ilevel)%domain(n)%im
  !         !
  !         m20=(dolevl(ilevel)%domain(n)%jm-dolevl(ilevel)%domain(n)%j0)/2 + &
  !              dolevl(ilevel)%domain(n)%j0 
  !         dolevl(ilevel)%domain(n)%subdom(:,1,:)%j0=dolevl(ilevel)%domain(n)%j0
  !         dolevl(ilevel)%domain(n)%subdom(:,1,:)%jm=m20
  !         dolevl(ilevel)%domain(n)%subdom(:,2,:)%j0=m20
  !         dolevl(ilevel)%domain(n)%subdom(:,2,:)%jm=dolevl(ilevel)%domain(n)%jm
  !         !
  !         m20=(dolevl(ilevel)%domain(n)%km-dolevl(ilevel)%domain(n)%k0)/2 + &
  !              dolevl(ilevel)%domain(n)%k0 
  !         dolevl(ilevel)%domain(n)%subdom(:,:,1)%k0=dolevl(ilevel)%domain(n)%k0
  !         dolevl(ilevel)%domain(n)%subdom(:,:,1)%km=m20
  !         dolevl(ilevel)%domain(n)%subdom(:,:,2)%k0=m20
  !         dolevl(ilevel)%domain(n)%subdom(:,:,2)%km=dolevl(ilevel)%domain(n)%km
  !         !
  !         nco2=0
  !         do k=1,2
  !         do j=1,2
  !         do i=1,2
  !           ncou=ncou+1
  !           nco2=nco2+1
  !           dolevl(ilevel+1)%domain(ncou)=dolevl(ilevel)%domain(n)%subdom(i,j,k)
  !           dolevl(ilevel)%domain(n)%next(nco2)=ncou
  !         enddo
  !         enddo
  !         enddo
  !         !
  !         ! if(mpirank==0) then
  !         !   write(*,'(2(A,I0))')'     ** domain at level: ',ilevel,' n= ',n
  !         !   write(*,'(12X,5(I3,A),I0)') dolevl(ilevel)%domain(n)%i0,' ~ ', &
  !         !                              dolevl(ilevel)%domain(n)%im,' | ', &
  !         !                              dolevl(ilevel)%domain(n)%j0,' ~ ', &
  !         !                              dolevl(ilevel)%domain(n)%jm,' | ', &
  !         !                              dolevl(ilevel)%domain(n)%k0,' ~ ', &
  !         !                              dolevl(ilevel)%domain(n)%km
  !         !   do k=1,2
  !         !   do j=1,2
  !         !   do i=1,2
  !         !     write(*,'((A,3I0))')'        ** sub-domain: ',i,j,k
  !         !     write(*,'(12X,5(I3,A),I0)')dolevl(ilevel)%domain(n)%subdom(i,j,k)%i0,' ~ ', &
  !         !                               dolevl(ilevel)%domain(n)%subdom(i,j,k)%im,' | ', &
  !         !                               dolevl(ilevel)%domain(n)%subdom(i,j,k)%j0,' ~ ', &
  !         !                               dolevl(ilevel)%domain(n)%subdom(i,j,k)%jm,' | ', &
  !         !                               dolevl(ilevel)%domain(n)%subdom(i,j,k)%k0,' ~ ', &
  !         !                               dolevl(ilevel)%domain(n)%subdom(i,j,k)%km
  !         !   enddo
  !         !   enddo
  !         !   enddo
  !         ! endif
  !         !
  !       enddo
  !       !
  !     enddo
  !     !
  !     ! search the extent
  !     do ilevel=1,nlevel
  !       ! 
  !       do n=1,dolevl(ilevel)%numdom
  !         !
  !         is=dolevl(ilevel)%domain(n)%i0
  !         ie=dolevl(ilevel)%domain(n)%im
  !         js=dolevl(ilevel)%domain(n)%j0
  !         je=dolevl(ilevel)%domain(n)%jm
  !         ks=dolevl(ilevel)%domain(n)%k0
  !         ke=dolevl(ilevel)%domain(n)%km
  !         !
  !         dolevl(ilevel)%domain(n)%xmin=1.d10
  !         dolevl(ilevel)%domain(n)%xmax=-1.d10
  !         do k=ks,ke
  !         do j=js,je
  !         do i=is,ie
  !           !
  !           if(i==is .or. j==js .or. k==ks .or. i==ie .or. j==je .or. k==ke) then
  !             ! only search the outerlayer
  !             dolevl(ilevel)%domain(n)%xmin(1)=min(dolevl(ilevel)%domain(n)%xmin(1),x(i,j,k,1))
  !             dolevl(ilevel)%domain(n)%xmin(2)=min(dolevl(ilevel)%domain(n)%xmin(2),x(i,j,k,2))
  !             dolevl(ilevel)%domain(n)%xmin(3)=min(dolevl(ilevel)%domain(n)%xmin(3),x(i,j,k,3))
  !             !
  !             dolevl(ilevel)%domain(n)%xmax(1)=max(dolevl(ilevel)%domain(n)%xmax(1),x(i,j,k,1))
  !             dolevl(ilevel)%domain(n)%xmax(2)=max(dolevl(ilevel)%domain(n)%xmax(2),x(i,j,k,2))
  !             dolevl(ilevel)%domain(n)%xmax(3)=max(dolevl(ilevel)%domain(n)%xmax(3),x(i,j,k,3))
  !           endif
  !           !
  !         enddo
  !         enddo
  !         enddo
  !         !
  !         ! if(mpirank==0) then
  !         !   write(*,'(2(A,I0))')'     ** domain at level: ',ilevel,' n= ',n
  !         !   write(*,'(12X,6(F10.6,A))') dolevl(ilevel)%domain(n)%xmin(1),' ~ ', &
  !         !                               dolevl(ilevel)%domain(n)%xmax(1),' | ', &
  !         !                               dolevl(ilevel)%domain(n)%xmin(2),' ~ ', &
  !         !                               dolevl(ilevel)%domain(n)%xmax(2),' | ', &
  !         !                               dolevl(ilevel)%domain(n)%xmin(3),' ~ ', &
  !         !                               dolevl(ilevel)%domain(n)%xmax(3),' '
  !         ! endif
  !         !
  !       enddo
  !       !
  !     enddo
  !     !
  !     lfirstcal=.false.
  !   endif
  !   !
  !   ! if(mpirank==0) then
  !   !   write(*,'(2(A,I0))')'     ** domain at level: ',ilevel,' n= ',n
  !   !   write(*,'(12X,6(F10.6,A))') dolevl(1)%domain(1)%xmin(1),' ~ ', &
  !   !                               dolevl(1)%domain(1)%xmax(1),' | ', &
  !   !                               dolevl(1)%domain(1)%xmin(2),' ~ ', &
  !   !                               dolevl(1)%domain(1)%xmax(2),' | ', &
  !   !                               dolevl(1)%domain(1)%xmin(3),' ~ ', &
  !   !                               dolevl(1)%domain(1)%xmax(3),' '
  !   ! endif
  !   !
  !   ! Now do the search work
  !   ijk=0
  !   !
  !   idomai_save=0
  !   idomai_save2=0
  !   !
  !   ilevel=1
  !   n=1
  !   !
  !   if( p(1)>=dolevl(ilevel)%domain(n)%xmin(1) .and. &
  !       p(2)>=dolevl(ilevel)%domain(n)%xmin(2) .and. &
  !       p(3)>=dolevl(ilevel)%domain(n)%xmin(3) .and. &
  !       p(1)<=dolevl(ilevel)%domain(n)%xmax(1) .and. &
  !       p(2)<=dolevl(ilevel)%domain(n)%xmax(2) .and. &
  !       p(3)<=dolevl(ilevel)%domain(n)%xmax(3) ) then
  !     ! first search the whole domain 
  !     !
  !     ncou=1
  !     idomai_save(ncou)=n
  !     !
  !     do while(ilevel<nlevel)
  !       !
  !       nco2=0
  !       do n=1,ncou
  !         !
  !         idomai=idomai_save(n)
  !         !
  !         ! go through branches
  !         do m=1,8
  !           !
  !           inxt=dolevl(ilevel)%domain(idomai)%next(m)
  !           !
  !           xmin=dolevl(ilevel+1)%domain(inxt)%xmin
  !           xmax=dolevl(ilevel+1)%domain(inxt)%xmax
  !           !
  !           if( p(1)>=xmin(1) .and. p(1)<=xmax(1) .and. &
  !               p(2)>=xmin(2) .and. p(2)<=xmax(2) .and. &
  !               p(3)>=xmin(3) .and. p(3)<=xmax(3) ) then
  !             !
  !             nco2=nco2+1
  !             idomai_save2(nco2)=inxt
  !             !
  !             ! write(*,'(4(A,I0))')'  ** rank: ',mpirank,             &
  !             !                     ' | locate p at level: ',ilevel+1, &
  !             !                     ' domain: ',inxt,' numb: ',nco2
  !             !
  !           endif
  !           !
  !         enddo
  !         !
  !       enddo
  !       !
  !       ilevel=ilevel+1
  !       idomai_save=idomai_save2
  !       ncou=nco2
  !     enddo
  !     !
  !     do n=1,ncou
  !       !
  !       idomai=idomai_save(n)
  !       !
  !       do k=dolevl(ilevel)%domain(idomai)%k0+1,dolevl(ilevel)%domain(idomai)%km
  !       do j=dolevl(ilevel)%domain(idomai)%j0+1,dolevl(ilevel)%domain(idomai)%jm
  !       do i=dolevl(ilevel)%domain(idomai)%i0+1,dolevl(ilevel)%domain(idomai)%im
  !         !
  !         if(nodeincell(cell(i,j,k),p)) then
  !           !
  !           ijk(1)=i+ig0
  !           ijk(2)=j+jg0
  !           ijk(3)=k+kg0
  !           !
  !         endif
  !         !
  !       enddo
  !       enddo
  !       enddo
  !       !
  !     enddo
  !     !
  !   endif
  !   !
  !   return
  !   !
  ! end function pngrid3d
  !+-------------------------------------------------------------------+
  !| The end of the function pngrid.                                   |
  !+-------------------------------------------------------------------+
  !
  !!+-------------------------------------------------------------------+
  !| This subroutine is to search a point in a structure mesh.         |
  !| return the number of i,j,k where the point is belong.             |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 27-Aug-2021: Created by J. Fang @ Appleton                        |
  !+-------------------------------------------------------------------+
  function pngrid3d(p) result(ijk)
    !
    use commvar,   only : im,jm,km
    use commarray, only : x,cell
    !
    ! arguments
    real(8),intent(in) :: p(3)
    integer :: ijk(3)
    !
    !
    ! local data
    type :: subdomaim
      integer :: level,i0,j0,k0,im,jm,km
      integer :: next(8)
      real(8) :: xmin(3),xmax(3)
    end type subdomaim
    !
    type(subdomaim),allocatable,save :: domain(:)
    integer,allocatable,save :: ndomai(:)
    integer,save :: nlevel,nalldom
    logical,save :: initial=.true.
    !
    integer :: ilevel,ic,id,i,j,k,in,jn,kn,                            &
               ncou,nco2,imm,jmm,kmm,is,js,ks,ie,je,ke,isub,n
    integer,allocatable :: idsaver(:),idstemp(:)
    logical :: lfound
    !
    if(initial) then
      ! divide the domain into many subdomains using a tree structure
      !
      nlevel=1
      do while(2**nlevel<max(im,jm,km))
        nlevel=nlevel+1
      enddo
      !
      nlevel=nlevel-2
      !
      allocate(ndomai(1:nlevel))
      ndomai(1)=1
      nalldom=1
      do ilevel=2,nlevel
        ndomai(ilevel)=ndomai(ilevel-1)*8
        nalldom=nalldom+ndomai(ilevel)
      enddo
      !
      ! print*,' ** number of levels  :',nlevel
      ! print*,' ** number of domains :',nalldom
      !
      allocate(domain(nalldom))
      !
      domain(1)%level=1
      domain(1)%i0=0; domain(1)%im=im
      domain(1)%j0=0; domain(1)%jm=jm
      domain(1)%k0=0; domain(1)%km=km
      !
      ic=1
      id=1
      ncou=1
      !
      do while(ncou<nalldom)
        !
        ! go through allocated domains
        do i=ic,id
          !
          nco2=0
          do kn=1,2
          do jn=1,2
          do in=1,2
            !
            ncou=ncou+1
            !
            nco2=nco2+1
            !
            ! assgin the children domains
            domain(i)%next(nco2)=ncou
            !
            imm=(domain(i)%im-domain(i)%i0)/2
            jmm=(domain(i)%jm-domain(i)%j0)/2
            kmm=(domain(i)%km-domain(i)%k0)/2
            !
            ! set the children domains
            domain(ncou)%level=domain(i)%level+1
            domain(ncou)%next=-1
            !
            if(in==1) then
              domain(ncou)%i0=domain(i)%i0
              domain(ncou)%im=domain(i)%i0 + imm
            elseif(in==2) then
              domain(ncou)%i0=domain(i)%i0 + imm
              domain(ncou)%im=domain(i)%im
            endif
            !
            if(jn==1) then
              domain(ncou)%j0=domain(i)%j0
              domain(ncou)%jm=domain(i)%j0 + jmm
            elseif(jn==2) then
              domain(ncou)%j0=domain(i)%j0 + jmm
              domain(ncou)%jm=domain(i)%jm
            endif
            !
            if(kn==1) then
              domain(ncou)%k0=domain(i)%k0
              domain(ncou)%km=domain(i)%k0 + kmm
            elseif(kn==2) then
              domain(ncou)%k0=domain(i)%k0 + kmm
              domain(ncou)%km=domain(i)%km
            endif
            !
            !
          enddo
          enddo
          enddo
          !
        enddo
        !
        ic=id+1
        id=ncou
        !
        ! print*,' ** search domain: ',ic,id
        !
      enddo
      !
      ! set search x,y,z extent
      do id=1,nalldom
        !
        ! only check the outer layer of the domain
        !
        is=domain(id)%i0
        ie=domain(id)%im
        js=domain(id)%j0
        je=domain(id)%jm
        ks=domain(id)%k0
        ke=domain(id)%km
        !
        domain(id)%xmin= 1.d10
        domain(id)%xmax=-1.d10
        do k=ks,ke
        do j=js,je
        do i=is,ie
          !
          if(i==is .or. j==js .or. k==ks .or. i==ie .or. j==je .or. k==ke) then
            ! only search the outerlayer
            domain(id)%xmin(1)=min(domain(id)%xmin(1),x(i,j,k,1))
            domain(id)%xmin(2)=min(domain(id)%xmin(2),x(i,j,k,2))
            domain(id)%xmin(3)=min(domain(id)%xmin(3),x(i,j,k,3))
            !
            domain(id)%xmax(1)=max(domain(id)%xmax(1),x(i,j,k,1))
            domain(id)%xmax(2)=max(domain(id)%xmax(2),x(i,j,k,2))
            domain(id)%xmax(3)=max(domain(id)%xmax(3),x(i,j,k,3))
          endif
          !
        enddo
        enddo
        enddo
        !
        ! if(id==1) then
        !   print*,' ** root domain id:',id,' level: ',domain(id)%level
        !   print*,' **  i extent:',domain(id)%i0,' ~ ',domain(id)%im
        !   print*,' **  j extent:',domain(id)%j0,' ~ ',domain(id)%jm
        !   print*,' **  k extent:',domain(id)%k0,' ~ ',domain(id)%km
        !   print*,' **  x extent:',domain(id)%xmin(1),' ~ ',domain(id)%xmax(1)
        !   print*,' **  y extent:',domain(id)%xmin(2),' ~ ',domain(id)%xmax(2)
        !   print*,' **  z extent:',domain(id)%xmin(3),' ~ ',domain(id)%xmax(3)
        ! endif
        ! if(id==nalldom) then
        !   print*,' ** highest domain id:',id,' level: ',domain(id)%level
        !   print*,' **  i extent:',domain(id)%i0,' ~ ',domain(id)%im
        !   print*,' **  j extent:',domain(id)%j0,' ~ ',domain(id)%jm
        !   print*,' **  k extent:',domain(id)%k0,' ~ ',domain(id)%km
        !   print*,' **  x extent:',domain(id)%xmin(1),' ~ ',domain(id)%xmax(1)
        !   print*,' **  y extent:',domain(id)%xmin(2),' ~ ',domain(id)%xmax(2)
        !   print*,' **  z extent:',domain(id)%xmin(3),' ~ ',domain(id)%xmax(3)
        ! endif
        !
      enddo
      !
      initial=.false.
      !
    endif
    !
    ! now the search the subdomain that possibly contain p
    !
    allocate(idsaver(ndomai(nlevel)),idstemp(ndomai(nlevel)))
    !
    id=1
    ilevel=1
    ncou=1
    idsaver(1)=1
    !
    if( p(1)>=domain(id)%xmin(1) .and. p(1)<=domain(id)%xmax(1) .and. &
        p(2)>=domain(id)%xmin(2) .and. p(2)<=domain(id)%xmax(2) .and. &
        p(3)>=domain(id)%xmin(3) .and. p(3)<=domain(id)%xmax(3) ) then
      !
      ! write(*,'(2(A,I0))')'  **  found p at domain: ',id,', level: ',domain(id)%level
      ! !
      ! print*,' **  x extent:',domain(id)%xmin(1),' ~ ',domain(id)%xmax(1)
      ! print*,' **  y extent:',domain(id)%xmin(2),' ~ ',domain(id)%xmax(2)
      ! print*,' **  z extent:',domain(id)%xmin(3),' ~ ',domain(id)%xmax(3)
      !
      lfound=.true.
      !
    else
      !
      ! print*,' ** the point is not in the main domain.'
      ijk=0
      !
      return
      !
    endif
    !
    do while(ilevel<nlevel .and. lfound)
      !
      ! check the extent of the subdomain
      lfound=.false.
      !
      nco2=0
      idstemp=0
      !
      do n=1,ncou
        !
        id=idsaver(n)
        !
        do isub=1,8
          !
          ic=domain(id)%next(isub)
          !
          if( p(1)>=domain(ic)%xmin(1) .and. p(1)<=domain(ic)%xmax(1) .and. &
              p(2)>=domain(ic)%xmin(2) .and. p(2)<=domain(ic)%xmax(2) .and. &
              p(3)>=domain(ic)%xmin(3) .and. p(3)<=domain(ic)%xmax(3) ) then
            !
            ! write(*,'(2(A,I0))')'  **  found p at domain: ',ic,', level: ',domain(ic)%level
            ! !
            ! print*,' **  x extent:',domain(ic)%xmin(1),' ~ ',domain(ic)%xmax(1)
            ! print*,' **  y extent:',domain(ic)%xmin(2),' ~ ',domain(ic)%xmax(2)
            ! print*,' **  z extent:',domain(ic)%xmin(3),' ~ ',domain(ic)%xmax(3)
            !
            nco2=nco2+1
            !
            idstemp(nco2)=ic
            !
            lfound=.true.
            ilevel=domain(ic)%level
            !
          endif
          !
        enddo
        !
      enddo
      !
      ncou=nco2
      idsaver=idstemp
      !
    enddo
    !
    ! now find p
    do n=1,ncou
      !
      id=idsaver(n)
      !
      is=domain(id)%i0
      ie=domain(id)%im
      js=domain(id)%j0
      je=domain(id)%jm
      ks=domain(id)%k0
      ke=domain(id)%km
      !
      do k=domain(id)%k0+1,domain(id)%km
      do j=domain(id)%j0+1,domain(id)%jm
      do i=domain(id)%i0+1,domain(id)%im
        !
        if(nodeincell(cell(i,j,k),p)) then
          !
          ijk(1)=i+ig0
          ijk(2)=j+jg0
          ijk(3)=k+kg0
          !
          ! print*,' ** found p at domain:',id,'i,j,k',i,j,k
          !
          return
          !
        endif
        !
      enddo
      enddo
      enddo
      !
    enddo
    !
  end function pngrid3d
  !+-------------------------------------------------------------------+
  !| The end of the function pngrid3d.                                 |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to decide is a line from a point intersects    |
  !| with a triangle.                                                  |
  !+-------------------------------------------------------------------+
  !| ref: https://zh.wikipedia.org/wiki/%E7%BA%BF%E9%9D%A2%E4%BA%A4%E7%82%B9
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 04-07-2021  | Created by J. Fang @ Warrington                     |
  !+-------------------------------------------------------------------+
  subroutine ray2triangle(tria,p,dir,pintersection,intertri)
    !
    use commtype, only :  triangle
    use commfunc,  only : dis2point,areatriangle
    !
    ! arguments
    type(triangle),intent(in) :: tria
    real(8),intent(in) :: p(3),dir(3)
    logical,intent(out) :: intertri
    real(8),intent(out) :: pintersection(3)
    !
    ! local data
    integer :: i,j
    real(8) :: vec1(3),vec2(3),intpoint(3),d,d2,ldn
    real(8) :: epsilon=1.d-12
    real(8) :: trxmin(3),trxmax(3)
    !
    !
    ldn=dot_product(dir,tria%normdir)
    !
    ! vec1=num1d3*(tria%a+tria%b+tria%c)-p
    vec1=tria%a-p
    d=dot_product(vec1,tria%normdir)
    !
    !
    if(abs(ldn)<=epsilon) then
      ! the ray line is parallel to the plane
      !
      intertri=.false.
      pintersection=(/1.d10,1.d10,1.d10/)
      !
    else
      !
      d=d/ldn
      !
      intpoint=abs(d)*dir+p
      !
      intertri=pointintriangle(tria,intpoint)
      !
      if(intertri) then
        pintersection=intpoint
      else
        pintersection=(/1.d10,1.d10,1.d10/)
      endif
      !
    endif
    !
    !
    return
    !
  end subroutine ray2triangle
  !+-------------------------------------------------------------------+
  !| The end of the subroutine ray2triangle.                           |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !|This function returns the distance of a point to a triangle.       |
  !+-------------------------------------------------------------------+
  !| ref: https://www.cnblogs.com/tenosdoit/p/4024413.html
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 15-07-2021  | Created by J. Fang @ Warrington                     |
  !+-------------------------------------------------------------------+
  subroutine dis2tri(tria,p,distance,intpoint,nodir,dir,debug)
    !
    use commtype, only : triangle
    use commfunc, only : dis2point
    !
    type(triangle),intent(in) :: tria
    real(8),intent(in) :: p(3)
    character(len=1),intent(in) :: dir
    logical,intent(in),optional :: debug
    real(8),intent(out) :: distance
    real(8),intent(out) :: intpoint(3),nodir(3)
    !
    ! local data
    real(8) :: unitvect(3),vectoa(3)
    real(8) :: var1
    logical :: ldeg,linout
    !
    if(present(debug)) then
      ldeg=debug
    else
      ldeg=.false.
    endif
    !
    if(dir=='+') then
      unitvect = tria%normdir 
    elseif(dir=='-') then
      unitvect =-tria%normdir 
    else
      stop ' ERROR @ dis2tri'
    endif
    !
    ! The distance from the plane to the point is the projection of the
    ! vector created between the point of interest off of the plane
    ! and any point on the plane onto the unit normal vector. 
    ! The first point in the plane [x(1), y(1), z(1)] is used.
    !
    vectoa=tria%a-p
    !
    var1 = dot_product(unitvect,vectoa)
    !
    intpoint=p+var1*unitvect
    !
    linout=pointintriangle(tria,intpoint)
    !
    if(linout) then
      ! if the intersection point is in the triangle
      distance=abs(var1)
      nodir   =tria%normdir 
    else
      ! calculate the distance between corners and p
      !
      distance=dis2point(p,tria%a)
      intpoint=tria%a
      !
      var1=dis2point(p,tria%b)
      if(var1<distance) then
        distance=var1
        intpoint=tria%b
      endif
      !
      var1=dis2point(p,tria%c)
      if(var1<distance) then
        distance=var1
        intpoint=tria%c
      endif
      !
      nodir=p-intpoint
      nodir=nodir/(sqrt(nodir(1)**2+nodir(2)**2+nodir(3)**2))
      !
    endif
    !
    return
    !
  end subroutine dis2tri
  !+-------------------------------------------------------------------+
  !| The end of the subroutine dis2tri.                                |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !|This function returns the distance of a point to a segment.        |
  !+-------------------------------------------------------------------+
  !| https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line     |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 18-07-2021  | Created by J. Fang @ Warrington                     |
  !+-------------------------------------------------------------------+
  subroutine dis2edge(edge,p,distance,intpoint,nodir,dir,debug)
    !
    use commtype, only : lsegment
    use commfunc, only : dis2point
    !
    type(lsegment),intent(in) :: edge
    real(8),intent(in) :: p(3)
    character(len=1),intent(in) :: dir
    logical,intent(in),optional :: debug
    real(8),intent(out) :: distance
    real(8),intent(out) :: intpoint(3),nodir(3)
    !
    ! local data
    real(8) :: unitvect(2),ip(2)
    real(8) :: var1,var2,var3,d1
    logical :: ldeg,linout
    real(8) :: epsilon=1.d-12
    !
    if(present(debug)) then
      ldeg=debug
    else
      ldeg=.false.
    endif
    !
    if(dir=='+') then
      unitvect = edge%normdir 
    elseif(dir=='-') then
      unitvect =-edge%normdir 
    else
      stop ' ERROR @ dis2tri'
    endif
    !
    !
    var1=(edge%b(1)-edge%a(1))*(edge%a(2)-p(2))
    var2=(edge%a(1)-p(1))     *(edge%b(2)-edge%a(2))
    var3=sqrt((edge%b(1)-edge%a(1))**2+(edge%b(2)-edge%a(2))**2)
    !
    d1=abs(var1-var2)/var3 
    !
    ip=p(1:2)+d1*unitvect
    !
    nodir(1:2)=unitvect
    nodir(3)=0.d0
    !
    var1=sqrt((p(1)-edge%a(1))**2+(p(2)-edge%a(2))**2)
    var2=sqrt((p(1)-edge%b(1))**2+(p(2)-edge%b(2))**2)
    !
    if( (ip(1)-edge%a(1))*(ip(1)-edge%b(1))<=epsilon .and. &
        (ip(2)-edge%a(2))*(ip(2)-edge%b(2))<=epsilon   ) then
      linout=.true.
      distance=d1
    else
      linout=.false.
      !
      if(var1<=var2) then
        ip=edge%a
        distance=var1
        !
      else
        ip=edge%b
        distance=var2
      endif
      !
      nodir(1:2)=p(1:2)-ip
      nodir(3)=0.d0
      !
      nodir=nodir/sqrt(nodir(1)**2+nodir(2)**2)
      !
    endif
    !
    intpoint(1:2)=ip(1:2)
    intpoint(3)  =0.d0
    !
    return
    !
  end subroutine dis2edge
  !+-------------------------------------------------------------------+
  !| The end of the subroutine dis2edge.                               |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to decide if a point is inside of a triangle.  |
  !+-------------------------------------------------------------------+
  !| ref: https://www.cnblogs.com/tenosdoit/p/4024413.html
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 04-07-2021  | Created by J. Fang @ Warrington                     |
  !+-------------------------------------------------------------------+
  function pointintriangle_tri(tria,p) result(lin)
    !
    use commtype, only :  triangle
    use commfunc,  only : areatriangle,cross_product
    !
    ! arguments
    type(triangle),intent(in) :: tria
    real(8),intent(in) :: p(3)
    logical :: lin
    !
    ! local data
    real(8) :: pa(3),pb(3),pc(3),t1(3),t2(3),t3(3),v0(3),v1(3),v2(3)
    real(8) :: a1,a2,a3,error,var1,var2,dot00,dot01,dot02,dot11,dot12, &
               inverDeno,u,v
    real(8) :: epsilon
    !
    ! epsilon=1.d-12
    ! !
    ! a1=areatriangle(tria%a,tria%b,p)
    ! a2=areatriangle(tria%a,tria%c,p)
    ! a3=areatriangle(tria%b,tria%c,p)
    ! !
    ! error=abs(a1+a2+a3-tria%area)
    ! if(error<epsilon) then
    !   ! the point is in the triangle
    !   pointintriangle=.true.
    ! elseif(a1+a2+a3>tria%area) then
    !   ! the point is outof the triangle
    !   pointintriangle=.false.
    ! else
    !   print*,' !! error in pointintriangle'
    !   print*,a1+a2+a3,tria%area,error
    ! endif
    !
    ! pointintriangle=sameside(tria%a,tria%b,tria%c,p) .and.             &
    !                 sameside(tria%b,tria%c,tria%a,p) .and.             &
    !                 sameside(tria%c,tria%a,tria%b,p)
    !
    ! pa=tria%a-p
    ! pb=tria%b-p
    ! pc=tria%c-p
    ! !
    ! t1=cross_product(pa,pb)
    ! t2=cross_product(pb,pc)
    ! t3=cross_product(pc,pa)
    ! !
    ! var1=dot_product(t1,t2)
    ! var2=dot_product(t1,t3)
    ! !
    ! if(var1>=0.d0 .and. var2>=0.d0) then
    !   pointintriangle=.true.
    ! else
    !   pointintriangle=.false.
    ! endif
    !PointinTriangle(Vector3 A, Vector3 B, Vector3 C, Vector3 P)
    
    v0 = tria%c - tria%a 
    v1 = tria%b - tria%a 
    v2 =      p - tria%a 

    dot00 = dot_product(v0,v0)
    dot01 = dot_product(v0,v1)
    dot02 = dot_product(v0,v2)
    dot11 = dot_product(v1,v1)
    dot12 = dot_product(v1,v2)

    inverDeno = 1.d0 / (dot00 * dot11 - dot01 * dot01) 
    
    u = (dot11 * dot02 - dot01 * dot12) * inverDeno

    if(isnan(u)) then
      print*,tria%a 
      print*,tria%b
      print*,tria%c 
      print*,tria%area
      stop ' !! ERROR @ pointintriangle_tri'
    endif
    ! print*,' ** u=',u
    !
    if (u < 0.d0 .or. u > 1.d0) then
      ! if u out of range, return directly
      lin=.false.
      return
    endif
    !
    v = (dot00 * dot12 - dot01 * dot02) * inverDeno
    ! print*,' ** v=',v
    if (v < 0.d0 .or. v > 1.d0) then
      ! if v out of range, return directly
      lin=.false.
      return 
    endif
    !
    if(u + v <= 1.d0) then
      lin=.true.
    else
      lin=.false.
    endif
    !
    return
    !
  end function pointintriangle_tri
  !
  function pointintriangle_nodes(a,b,c,p) result(lin)
    !
    use commfunc,  only : areatriangle,cross_product
    !
    ! arguments
    real(8),intent(in) :: a(3),b(3),c(3),p(3)
    logical :: lin
    !
    ! local data
    real(8) :: pa(3),pb(3),pc(3),t1(3),t2(3),t3(3),v0(3),v1(3),v2(3)
    real(8) :: a1,a2,a3,error,var1,var2,dot00,dot01,dot02,dot11,dot12, &
               inverDeno,u,v
    real(8) :: epsilon
    !
    v0 = c - a 
    v1 = b - a 
    v2 = p - a 

    dot00 = dot_product(v0,v0)
    dot01 = dot_product(v0,v1)
    dot02 = dot_product(v0,v2)
    dot11 = dot_product(v1,v1)
    dot12 = dot_product(v1,v2)

    inverDeno = 1.d0 / (dot00 * dot11 - dot01 * dot01) 
    
    u = (dot11 * dot02 - dot01 * dot12) * inverDeno

    if(isnan(u)) then
      print*,a 
      print*,b
      print*,c
      stop ' !! ERROR @ pointintriangle_nodes'
    endif
    ! print*,' ** u=',u
    !
    if (u < 0.d0 .or. u > 1.d0) then
      ! if u out of range, return directly
      lin=.false.
      return
    endif
    !
    v = (dot00 * dot12 - dot01 * dot02) * inverDeno
    ! print*,' ** v=',v
    if (v < 0.d0 .or. v > 1.d0) then
      ! if v out of range, return directly
      lin=.false.
      return 
    endif
    !
    if(u + v <= 1.d0) then
      lin=.true.
    else
      lin=.false.
    endif
    !
    return
    !
  end function pointintriangle_nodes
  !+-------------------------------------------------------------------+
  !| The end of the subroutine pointintriangle.                        |
  !+-------------------------------------------------------------------+
  !!
  !+-------------------------------------------------------------------+
  !| This subroutine is to whether two vectors v1 and v2 point to the  |
  !| same direction.                                                   |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 04-07-2021  | Added by J. Fang @ Warrington                       |
  !+-------------------------------------------------------------------+
  pure logical function sameside(a,b,c,p)
    !
    use commfunc, only: cross_product
    !
    ! arguments
    real(8),intent(in) :: a(3),b(3),c(3),p(3)
    !
    ! local data
    real(8) :: ab(3),ac(3),ap(3),v1(3),v2(3),var1
    !
    ab = b - a
    ac = c - a
    ap = p - a
    !
    v1 = cross_product(ab,ac)
    v2 = cross_product(ab,ap)
    !
    ! v1 and v2 should point to the same direction
    var1=dot_product(v1,v2)
    !
    if(var1>=0.d0) then
      sameside=.true.
    else
      sameside=.false.
    endif  
    !
  end function sameside
  !+-------------------------------------------------------------------+
  !| The end of the function sameside.                                 |
  !+-------------------------------------------------------------------+
  !!
  !+-------------------------------------------------------------------+
  !| This function is to return the boundary node. which is the.       |
  !| shortest distance to a solid point.                               |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 04-07-2021  | Added by J. Fang @ Warrington                       |
  !+-------------------------------------------------------------------+
  subroutine polyhedron_bound_search(asolid,p,bnode,dir,debug)
    !
    use commtype, only : solid,sboun
    !
    ! arguments
    type(solid),intent(in) :: asolid
    real(kind=8),intent(in) :: p(3)
    character(len=1),intent(in) :: dir
    logical,intent(in),optional :: debug
    type(sboun),intent(out) :: bnode
    !
    ! local data
    integer :: jface,jedge,jfsave
    real(8) :: dismin,dis,var1
    real(8) :: nodeb(3),nodesave(3),vec(3),bnx(3)
    logical :: ldeg
    !
    if(present(debug)) then
      ldeg=debug
    else
      ldeg=.false.
    endif
    !
    dismin=1.d10
    !
    if(ndims==2) then
      !
      do jedge=1,asolid%num_edge
        !
        call dis2edge(asolid%edge(jedge),p,dis,nodeb,bnx,dir,debug=ldeg)
        !
        if(dis<dismin) then
          !
          bnode%x=nodeb
          bnode%normdir=bnx
          !
          dismin=dis
          !
        endif
        !
      enddo
      !
    elseif(ndims==3) then
      !
      do jface = 1, asolid%num_face
        !
        call dis2tri(asolid%face(jface),p,dis,nodeb,bnx,dir,debug=ldeg)
        !
        if(dis<dismin) then
          !
          bnode%x=nodeb
          bnode%normdir=bnx
          !
          dismin=dis
          !
        endif
        !
      end do
      !
    endif
    !
    return
    !
  end subroutine polyhedron_bound_search
  !+-------------------------------------------------------------------+
  !| The end of the function polyhedron_bound_search.                  |
  !+-------------------------------------------------------------------+
  !!
  !+-------------------------------------------------------------------+
  !
  ! POINT_IN_POLYGON determines if a point is inside a polygon
  !
  ! Discussion:
  !
  ! If the points ( x(i), y(i) ) ( i = 1, 2, ..., n ) are,
  ! in this cyclic order, the vertices of a simple closed polygon and
  ! (x0,y0) is a point not on any side of the polygon, then the
  ! procedure determines, by setting "point_in_polygon" to TRUE or FALSE,
  ! whether (x0,y0) lies in the interior of the polygon.
  !
  ! Licensing:
  !
  !   This code is distributed under the GNU LGPL license. 
  !
  ! Modified:
  !
  !   07 November 2016
  !
  ! Author:
  !
  !   John Burkardt
  !
  ! Reference:
  !
  !   Moshe Shimrat,
  !   ACM Algorithm 112,
  !   Position of Point Relative to Polygon,
  !   Communications of the ACM,
  !   Volume 5, Number 8, page 434, August 1962.
  !
  !   Richard Hacker,
  !   Certification of Algorithm 112,
  !   Communications of the ACM,
  !   Volume 5, Number 12, page  606, December 1962.
  !
  ! Parameters:
  !
  !   Input, integer ( kind = 4 ) N, the number of nodes or vertices in 
  !   the polygon.  N must be at least 3.
  !
  !   Input, real ( kind = 8 ) V(2,N), the vertices of the polygon.
  !
  !   Input, real ( kind = 8 ) P(2), the coordinates of the point to be tested.
  !
  !   Output, logical ( kind = 4 ) INSIDE, is TRUE if the point is
  !   inside the polygon.
  !
  !+-------------------------------------------------------------------+
  function polyhedron_contains_point_2d ( asolid, p ) result (inside)
    !
    use commtype, only : solid
    !
    ! arguments
    type(solid),intent(in) :: asolid
    real(kind=8),intent(in) :: p(3)
    logical :: inside
    !
    ! local data
    integer :: i
    real(8) :: t,var1,var2,var3,epsilon
    real(8) :: a(2),b(2)
    
    epsilon=1.d-16
    !
    inside = .false.

    do i = 1, asolid%num_edge
      !
      a(:)=asolid%edge(i)%a(1:2)
      b(:)=asolid%edge(i)%b(1:2)
      !
      var1=sqrt((p(1)-a(1))**2+(p(2)-a(2))**2)
      var2=sqrt((p(1)-b(1))**2+(p(2)-b(2))**2)
      var3=sqrt((a(1)-b(1))**2+(a(2)-b(2))**2)
      !
      if(var1+var2-var3<=epsilon) then
        inside=.true.
      else
        !
        if(p(2)>b(2) .eqv. p(2)<=a(2)) then
          !
          if(abs(b(2)-p(2))>epsilon) then
            t=(p(1)-a(1))-(p(2)-a(2))*(b(1)-p(1))/(b(2)-p(2))
          else
            t=(p(1)-b(1))-(p(2)-b(2))*(a(1)-p(1))/(a(2)-p(2))
          endif
          !
          if ( t < 0.d0 ) then
            inside = .not. inside
          end if
        end if
        !
      endif
      !
    end do
    !
    return

  end function polyhedron_contains_point_2d
  !+-------------------------------------------------------------------+
  !| The end of the function polyhedron_contains_point_2d.             |
  !+-------------------------------------------------------------------+
  !!
  !+-------------------------------------------------------------------+
  ! ref: https://searchcode.com/codesearch/view/13427361/
  !
  !! POLYHEDRON_CONTAINS_POINT_3D determines if a point is inside a polyhedron.
  !  R. Schuh - October 2006
  !
  !  PIP (point inside polyhedron)
  !
  !  Discussion:
  !
  !    The reference states that the polyhedron should be simple (that
  !    is, the faces should form a single connected surface), and that 
  !    the individual faces should be consistently oriented.
  !
  !    However, the polyhedron does not, apparently, need to be convex.
  !
  !  Modified:
  !
  !    30 August 2005
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Reference:
  !
  !    Paulo Cezar Pinto Carvalho, Paulo Roma Cavalcanti,
  !    Point in Polyhedron Testing Using Spherical Polygons,
  !    in Graphics Gems V,
  !    edited by Alan Paeth,
  !    Academic Press, 1995, T385.G6975.
  !
  !  Parameters:
  !
  !    input, integer node_num, the number of vertices.
  !
  !    input, integer face_num, the number of faces.
  !
  !    input, integer face_order_max, the maximum order of any face.
  !
  !    input, real(kind=8) v(3,node_num), the coordinates of the vertices.
  !
  !    input, integer face_order(face_num), the order of each face.
  !
  !    input, integer face_point(face_order_max,face_num), the indices of the
  !    nodes that make up each face.
  !
  !    input, real(kind=8) p(3), the point to be tested.
  !
  !    output, logical inside, is true if the point 
  !    is inside the polyhedron.
  !+-------------------------------------------------------------------+
  function polyhedron_contains_point_3d ( asolid, p,debug ) result (inside)
    !
    use commtype, only : solid
    !
    ! arguments
    type(solid),intent(in) :: asolid
    real(kind=8),intent(in) :: p(3)
    logical,intent(in),optional :: debug
    !
    ! local data
    integer, parameter :: dim_num = 3, face_order_max = 3
    !
    real(kind=8) area
    integer jface
    integer face_order
    logical inside
    real(kind=8), parameter :: pi = 3.141592653589793d+00
    real(kind=8), parameter :: epsilon = 1.d-12
    real(kind=8) solid_angle
    real(kind=8) v_face(dim_num,face_order_max)
    real(8) :: p2(3),norm(3),var1
    
    norm=p-asolid%xcen
    var1=sqrt(norm(1)**2+norm(2)**2+norm(3)**2)
    if(var1>epsilon*10.d0) then
      norm=norm/var1
      p2=p-norm*epsilon
    else
      p2=p
    endif

    area = 0.0d+00
    face_order=3
    
    do jface = 1, asolid%num_face
  
      v_face(:,1) = asolid%face(jface)%a
      v_face(:,2) = asolid%face(jface)%b
      v_face(:,3) = asolid%face(jface)%c
  
      call geo_polygon_solid_angle_3d(face_order,v_face,p2,solid_angle)
  
      area = area + solid_angle
  
    end do
    !
    !  area should be -4*pi, 0, or 4*pi.
    !  so this test should be quite safe!
    !
    if(present(debug)) then
      print*,'area=',area
    endif
    !
    if ( area < -2.0d+00 * pi .or. 2.0d+00 * pi < area ) then
      inside = .true.
    else
      inside = .false.
    end if
  
    return

  end function polyhedron_contains_point_3d
  !+-------------------------------------------------------------------+
  !| The end of the subroutine polyhedron_contains_point_3d            |
  !+-------------------------------------------------------------------+
  !!
  !+-------------------------------------------------------------------+
  !
  !! POLYGON_SOLID_ANGLE_3D: projected solid angle of a 3D plane polygon.
  !
  !  Discussion:
  !
  !    A point P is at the center of the unit sphere.  A planar polygon
  !    is to be projected onto the surface of the unit sphere, by drawing
  !    the ray from P to each polygonal vertex, and noting where this ray
  !    intersects the unit sphere.  The area of the projected polygon is
  !    equal to the solid angle, since we are considering the unit sphere.
  !
  !  Modified:
  !
  !    30 August 2005
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Reference:
  !
  !    Paulo Cezar Pinto Carvalho, Paulo Roma Cavalcanti,
  !    Point in Polyhedron Testing Using Spherical Polygons,
  !    in Graphics Gems V,
  !    edited by Alan Paeth,
  !    Academic Press, 1995, T385.G6975.
  !
  !  Parameters:
  !
  !    Input, integer N, the number of vertices.
  !
  !    Input, real(kind=8) V(3,N), the coordinates of the vertices.
  !
  !    Input, real(kind=8) P(3), the point at the center of the unit sphere.
  !
  !    Output, double SOLID_ANGLE, the solid angle subtended
  !    by the polygon, as projected onto the unit sphere around the point P.
  !
  !+-------------------------------------------------------------------+
  subroutine geo_polygon_solid_angle_3d ( n, v, p, solid_angle )
  
    implicit none
  
    integer, parameter :: dim_num = 3
    integer n
  
    real(kind=8) a(dim_num)
    real(kind=8) angle
    real(kind=8) area
    real(kind=8) b(dim_num)
    integer i
    integer j
    integer jp1
    real(kind=8) normal1(dim_num)
    real(kind=8) normal1_norm
    real(kind=8) normal2(dim_num)
    real(kind=8) normal2_norm
    real(kind=8) p(dim_num)
    real(kind=8), parameter :: pi = 3.141592653589793D+00
    real(kind=8) plane(dim_num)
    real(kind=8) r1(dim_num)
    real(kind=8) s
    real(kind=8) solid_angle
    real(kind=8) v(dim_num,n)
  
    if ( n < 3 ) then
      solid_angle = 0.0D+00
      return
    end if
  
    call geo_polygon_normal_3d  ( n, v, plane )
   
    a(1:dim_num) = v(1:dim_num,n) - v(1:dim_num,1)
  
    area = 0.0D+00
  
    do j = 1, n
  
      r1(1:dim_num) = v(1:dim_num,j) - p(1:dim_num)
  
      jp1 = i_wrap ( j + 1, 1, n )
  
      b(1:dim_num) = v(1:dim_num,jp1) - v(1:dim_num,j)
  
      call geo_dvec_cross_3d ( a, r1, normal1 )
  
      normal1_norm = dvec_length ( dim_num, normal1 )
  
      call geo_dvec_cross_3d ( r1, b, normal2 )
  
      normal2_norm = dvec_length ( dim_num, normal2 )
      
      s = dot_product ( normal1(1:dim_num), normal2(1:dim_num) ) &
        / ( normal1_norm * normal2_norm )
  
      angle = arc_cosine ( s )
  
      s = dvec_triple_product ( b, a, plane )
  
      if ( 0.0D+00 < s ) then
        area = area + pi - angle
      else
        area = area + pi + angle
      end if
  
      a(1:dim_num) = -b(1:dim_num)
  
    end do
  
    area = area - pi * real ( n - 2, kind = 8 )
  
    if ( 0.0D+00 < dot_product ( plane(1:dim_num), r1(1:dim_num) ) ) then
      solid_angle = -area
    else
      solid_angle = area
    end if
  
    return

  end subroutine geo_polygon_solid_angle_3d
  !+-------------------------------------------------------------------+
  !| The end of subroutine geo_polygon_solid_angle_3d                  |
  !+-------------------------------------------------------------------+
  !!
  !+-------------------------------------------------------------------+
  !
  !! ARC_COSINE computes the arc cosine function, with argument truncation.
  !
  !  Discussion:
  !
  !    If you call geo_your system ACOS routine with an input argument that is
  !    even slightly outside the range [-1.0, 1.0 ], you may get an unpleasant 
  !    surprise (I did).
  !
  !    This routine simply truncates arguments outside the range.
  !
  !  Modified:
  !
  !    02 December 2000
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, real(kind=8) C, the argument.
  !
  !    Output, real(kind=8) ARC_COSINE, an angle whose cosine is C.
  !+-------------------------------------------------------------------+
  real(kind=8) function arc_cosine ( c )
  
    !
    implicit none
  
    real(kind=8) c
    real(kind=8) c2
  
    c2 = c
    c2 = max ( c2, -1.0D+00 )
    c2 = min ( c2, +1.0D+00 )
  
    arc_cosine = acos ( c2 )
  
    return

  end function arc_cosine
  !+-------------------------------------------------------------------+
  !| The end of the function arc_cosine                                |
  !+-------------------------------------------------------------------+
  !!
  !+-------------------------------------------------------------------+
  !
  !! DVEC_CROSS_3D computes the cross product of two vectors in 3D.
  !
  !  Discussion:
  !
  !    The cross product in 3D can be regarded as the determinant of the
  !    symbolic matrix:
  !
  !          |  i  j  k |
  !      det | x1 y1 z1 |
  !          | x2 y2 z2 |
  !
  !      = ( y1 * z2 - z1 * y2 ) * i
  !      + ( z1 * x2 - x1 * z2 ) * j
  !      + ( x1 * y2 - y1 * x2 ) * k
  !
  !  Modified:
  !
  !    07 August 2005
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, real(kind=8) V1(3), V2(3), the two vectors.
  !
  !    Output, real(kind=8) V3(3), the cross product vector.
  !
  !+-------------------------------------------------------------------+
  subroutine geo_dvec_cross_3d ( v1, v2, v3 )
  
    implicit none
  
    integer, parameter :: dim_num = 3
  
    real(kind=8) v1(dim_num)
    real(kind=8) v2(dim_num)
    real(kind=8) v3(dim_num)
  
    v3(1) = v1(2) * v2(3) - v1(3) * v2(2)
    v3(2) = v1(3) * v2(1) - v1(1) * v2(3)
    v3(3) = v1(1) * v2(2) - v1(2) * v2(1)
  
    return

  end subroutine geo_dvec_cross_3d
  !+-------------------------------------------------------------------+
  !| The end of the subroutine geo_dvec_cross_3d                       |
  !+-------------------------------------------------------------------+
  !!
  !+-------------------------------------------------------------------+
  !
  !! DVEC_LENGTH returns the Euclidean length of a vector.
  !
  !  Modified:
  !
  !    08 August 2005
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer DIM_NUM, the spatial dimension.
  !
  !    Input, real(kind=8) X(DIM_NUM), the vector.
  !
  !    Output, real(kind=8) DVEC_LENGTH, the Euclidean length of the vector.
  !+-------------------------------------------------------------------+
  real(kind=8) function dvec_length ( dim_num, x )
  
    implicit none
  
    integer dim_num
     
    real(kind=8) x(dim_num)
  
    dvec_length = sqrt ( sum ( ( x(1:dim_num) )**2 ) )
  
    return

  end function dvec_length
  !+-------------------------------------------------------------------+
  !| The end of the function dvec_length                               |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !
  !! DVEC_TRIPLE_PRODUCT finds the triple product in 3D.
  !
  !  Discussion:
  !
  !    [A,B,C] = A dot ( B cross C ) 
  !            = B dot ( C cross A )
  !            = C dot ( A cross B )
  !
  !    The volume of a parallelepiped, whose sides are given by
  !    vectors A, B, and C, is abs ( A dot ( B cross C ) ).
  !
  !    Three vectors are coplanar if and only if their triple product vanishes.
  !
  !  Modified:
  !
  !    07 August 2005
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Reference:
  !
  !    Eric Weisstein,
  !    "Scalar Triple Product",
  !    CRC Concise Encyclopedia of Mathematics, 1999
  !
  !  Parameters:
  !
  !    Input, real(kind=8) V1(3), V2(3), V3(3), the vectors.
  !
  !    Output, real(kind=8) DVEC_TRIPLE_PRODUCT, the triple product.
  !+-------------------------------------------------------------------+
  real(kind=8) function dvec_triple_product ( v1, v2, v3 )
  
    implicit none
  
    integer, parameter :: dim_num = 3
  
    real(kind=8) v1(dim_num)
    real(kind=8) v2(dim_num)
    real(kind=8) v3(dim_num)
    real(kind=8) v4(dim_num)
  
    call geo_dvec_cross_3d ( v2, v3, v4 )
  
    dvec_triple_product = dot_product ( v1(1:dim_num), v4(1:dim_num) )
  
    return

  end function dvec_triple_product
  !+-------------------------------------------------------------------+
  !| The end of the function dvec_triple_product                       |
  !+-------------------------------------------------------------------+
  !!
  !+-------------------------------------------------------------------+
  !
  !! geo_i_modp returns the nonnegative remainder of integer division.
  !
  !  Discussion:
  !
  !    If
  !      NREM = geo_i_modp ( I, J )
  !      NMULT = ( I - NREM ) / J
  !    then
  !      I = J * NMULT + NREM
  !    where NREM is always nonnegative.
  !
  !    The MOD function computes a result with the same sign as the
  !    quantity being divided.  Thus, suppose you had an angle A,
  !    and you wanted to ensure that it was between 0 and 360.
  !    Then mod(A,360) would do, if A was positive, but if A
  !    was negative, your result would be between -360 and 0.
  !
  !    On the other hand, geo_i_modp(A,360) is between 0 and 360, always.
  !
  !  Examples:
  !
  !        I     J     MOD  geo_i_modp    Factorization
  !
  !      107    50       7       7    107 =  2 *  50 + 7
  !      107   -50       7       7    107 = -2 * -50 + 7
  !     -107    50      -7      43   -107 = -3 *  50 + 43
  !     -107   -50      -7      43   -107 =  3 * -50 + 43
  !
  !  Modified:
  !
  !    02 March 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer I, the number to be divided.
  !
  !    Input, integer J, the number that divides I.
  !
  !    Output, integer geo_i_modp, the nonnegative remainder when I is
  !    divided by J.
  !
  !+-------------------------------------------------------------------+
  integer function geo_i_modp ( i, j )
  
    implicit none
  
    integer i
    integer j
  
    if ( j == 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'geo_i_modp - Fatal error!'
      write ( *, '(a,i8)' ) '  geo_i_modp ( I, J ) called with J = ', j
      stop
    end if
  
    geo_i_modp = mod ( i, j )
  
    if ( geo_i_modp < 0 ) then
      geo_i_modp = geo_i_modp + abs ( j )
    end if
  
    return

  end function geo_i_modp
  !+-------------------------------------------------------------------+
  !| The end of the function geo_i_modp                                |
  !+-------------------------------------------------------------------+
  !!
  !+-------------------------------------------------------------------+
  !
  !! I_WRAP forces an integer to lie between given limits by wrapping.
  !
  !  Example:
  !
  !    ILO = 4, IHI = 8
  !
  !    I  I_WRAP
  !
  !    -2     8
  !    -1     4
  !     0     5
  !     1     6
  !     2     7
  !     3     8
  !     4     4
  !     5     5
  !     6     6
  !     7     7
  !     8     8
  !     9     4
  !    10     5
  !    11     6
  !    12     7
  !    13     8
  !    14     4
  !
  !  Modified:
  !
  !    19 August 2003
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer IVAL, an integer value.
  !
  !    Input, integer ILO, IHI, the desired bounds for the integer value.
  !
  !    Output, integer I_WRAP, a "wrapped" version of IVAL.
  !
  !+-------------------------------------------------------------------+
  integer function i_wrap ( ival, ilo, ihi )
  
    implicit none
  
    integer ihi
    integer ilo
    integer ival
    integer jhi
    integer jlo
    integer wide
  
    jlo = min ( ilo, ihi )
    jhi = max ( ilo, ihi )
  
    wide = jhi - jlo + 1
  
    if ( wide == 1 ) then
      i_wrap = jlo
    else
      i_wrap = jlo + geo_i_modp ( ival - jlo, wide )
    end if
  
    return

  end function i_wrap
  !+-------------------------------------------------------------------+
  !| The end of the function i_wrap                                    |
  !+-------------------------------------------------------------------+
  !!
  !+-------------------------------------------------------------------+
  !
  !! POLYGON_NORMAL_3D computes the normal vector to a polygon in 3D.
  !
  !  Discussion:
  !
  !    If the polygon is planar, then this calculation is correct.
  !
  !    Otherwise, the normal vector calculated is the simple average
  !    of the normals defined by the planes of successive triples
  !    of vertices.
  !
  !    If the polygon is "almost" planar, this is still acceptable.
  !    But as the polygon is less and less planar, so this averaged normal
  !    vector becomes more and more meaningless.
  !
  !  Modified:
  !
  !    12 August 2005
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Reference:
  !
  !    Paulo Cezar Pinto Carvalho, Paulo Roma Cavalcanti,
  !    Point in Polyhedron Testing Using Spherical Polygons,
  !    in Graphics Gems V,
  !    edited by Alan Paeth,
  !    Academic Press, 1995, T385.G6975.
  !
  !  Parameters:
  !
  !    Input, integer N, the number of vertices.
  !
  !    Input, real(kind=8) V(3,N), the coordinates of the vertices.
  !
  !    Output, real(kind=8) NORMAL(3), the averaged normal vector
  !    to the polygon. 
  !
  !+-------------------------------------------------------------------+
  subroutine geo_polygon_normal_3d ( n, v, normal ) 
  
    implicit none
  
    integer, parameter :: dim_num = 3
    integer n
  
    integer i
    integer j
    real(kind=8) normal(dim_num)
    real(kind=8) normal_norm
    real(kind=8) p(dim_num)
    real(kind=8) v(dim_num,n)
    real(kind=8) v1(dim_num)
    real(kind=8) v2(dim_num)
  
    normal(1:dim_num) = 0.0D+00
  
    v1(1:dim_num) = v(1:dim_num,2) - v(1:dim_num,1)
  
    do j = 3, n
  
      v2(1:dim_num) = v(1:dim_num,j) - v(1:dim_num,1)
  
      call geo_dvec_cross_3d ( v1, v2, p )
  
      normal(1:dim_num) = normal(1:dim_num) + p(1:dim_num)
  
      v1(1:dim_num) = v2(1:dim_num)
  
    end do
  !
  !  Normalize.
  !
    normal_norm = dvec_length ( dim_num, normal )
  
    if ( normal_norm == 0.0D+00 ) then
      return
    end if
  
    normal(1:dim_num) = normal(1:dim_num) / normal_norm
  
    return

  end subroutine geo_polygon_normal_3d
  !+-------------------------------------------------------------------+
  !| The end of subroutine geo_polygon_normal_3d                       |
  !+-------------------------------------------------------------------+
  !!
end module geom