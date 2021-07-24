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
                       datasync,ptime,irk,jrk,krk,ig0,jg0,kg0
  use commvar,  only : ndims,ks,ke,hm,hm,lfftk,ctime,im,jm,km
  use tecio
  use stlaio,  only: get_unit
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
      call gridinsolid
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
    integer :: js
    type(solid),pointer :: psolid
    !
    if(mpirank==0) then
      !
      do js=1,nsolid
        !
        call solidrange(immbody(js))
        !
        if(ndims==2) then
          call solidreduc(immbody(js))
        endif
        !
      enddo
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
  !| This subroutine is used to calculate which nodes in/out of solid. |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 04-Jul-2021: Created by J. Fang @ Appleton                        |
  !+-------------------------------------------------------------------+
  subroutine gridinsolid
    !
    use commtype,  only : solid,triangle,sboun,nodcel
    use commvar,   only : immbody,nsolid,immbnod,dxyzmax,dxyzmin,      &
                          npdci,npdcj,npdck,ndims
    use commarray, only : x,nodestat,cell
    use commfunc,  only : dis2point,dis2point2,matinv4
    use commcal,   only : ijkin
    use parallel,  only : ig0,jg0,kg0,pmerg,syncinmg,syncisup,psum,    &
                          syncweig
    !
    ! local data
    integer :: i,j,k,jsd,jfc,counter,ninters,n,n1,m,ii,jj,kk,bignum,   &
               jb,kb,jdir,jx,iss,jss,kss,nc_f,nc_g,nc_b
    integer :: fh,ncou,ke1
    type(solid),pointer :: pso
    logical :: crossface
    real(8),allocatable :: rnodestat(:,:,:)
    integer :: snodes(27,3),i_cell(3)
    real(8) :: epsilon
    real(8) :: dist,distmin,var1,var2
    real(8) :: xmin(3),xmax(3),xcell(4,3),xnorm(4,3)
    real(8) :: Tm1(4,4),Tm2(4,4),Ti1(4,4),Ti2(4,4)
    !
    type(sboun),allocatable :: bnodes(:)
    logical :: liout,ljout,lkout,lin
    logical,allocatable :: marker(:,:,:)
    !
    if(lio) print*,' ** calculating grid in solid'
    !
    epsilon=1.d-10
    !
    ninters=10
    !
    allocate(rnodestat(0:im,0:jm,0:km))
    !
    nodestat=0
    !
    rnodestat=0.d0
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
    do jsd=1,nsolid
      !
      pso=>immbody(jsd)
      !
      do k=0,km
      do j=0,jm
      do i=0,im
        !
        if(x(i,j,k,1)<pso%xmin(1) .or.  x(i,j,k,1)>pso%xmax(1) .or.    &
           x(i,j,k,2)<pso%xmin(2) .or.  x(i,j,k,2)>pso%xmax(2) .or.    &
           x(i,j,k,3)<pso%xmin(3) .or.  x(i,j,k,3)>pso%xmax(3) ) then
          !
          nodestat(i,j,k)=0
          ! fluids 
          !
        else
          !
          ! to calculation the intersection between nodes and face
          ! !
          if(ndims==2) then
            crossface=polyhedron_contains_point_2d(pso,x(i,j,k,:))
          elseif(ndims==3) then
            crossface=polyhedron_contains_point_3d(pso,x(i,j,k,:))
          else
            stop ' !! ERROR1 @ gridinsolid' 
          endif
          !
          if(crossface) then
            ! ths point is in the solid
            nodestat(i,j,k)=5
          else
            ! fluids
            nodestat(i,j,k)=0
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
    call dataswap(nodestat)
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
          ! for solid nodes
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
          do kk=-1,1,2
            if(nodestat(i,j,k+kk)==n) then
              marker(i,j,k)=.true.
              exit
            endif
          enddo
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
        do kk=-1,1,2
          if(nodestat(i,j,k+kk)>0) then
            marker(i,j,k)=.true.
            exit
          endif
        enddo
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
    ! set cell state. currently, 2D only
    do k=kss,km
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
    enddo
    !
    bignum=im*jm*(km+1)
    allocate(bnodes(bignum))
    !
    ! to get the nodes and distance of inner solid nodes to boundary 
    counter=0
    nc_f=0
    nc_b=0
    nc_g=0
    !
    marker=.false.
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
          if(dis2point(bnodes(counter)%x,x(i,j,k,:))<1.d-6) then
            bnodes(counter)%nodetype='b'
            nc_b=nc_b+1
            !
            marker(i,j,k)=.true.
            ! boundary marker
          else
            bnodes(counter)%nodetype='g'
            nc_g=nc_g+1
          endif
          !
          !
        elseif(nodestat(i,j,k)==-1) then
          ! force point
          !
          ! counter=counter+1
          ! !
          ! call polyhedron_bound_search(pso,x(i,j,k,:),bnodes(counter),dir='-')
          ! !
          ! bnodes(counter)%igh(1)=i+ig0
          ! bnodes(counter)%igh(2)=j+jg0
          ! bnodes(counter)%igh(3)=k+kg0
          ! !
          ! bnodes(counter)%ximag(:)=2.d0*x(i,j,k,:)-bnodes(counter)%x(:)
          ! !
          ! bnodes(counter)%nodetype='f'
          ! !
          ! nc_f=nc_f+1
        endif
        !
      enddo
      enddo
      enddo
      !
    enddo
    !
    call pmerg(var=bnodes,nvar=counter,vmerg=immbnod)
    !
    nc_b=psum(nc_b)
    nc_f=psum(nc_f)
    nc_g=psum(nc_g)
    if(lio) then
      write(*,'(A,I0)')'  ** number of boundary nodes: ',size(immbnod)
      write(*,'(A,I0)')'     **   ghost nodes: ',nc_g
      write(*,'(A,I0)')'     ** forcing nodes: ',nc_f
      write(*,'(A,I0)')'    ** boundary nodes: ',nc_b
    endif
    !
    if(ndims==2) then
      !
      do jb=1,size(immbnod)
        !
        ! search the cell that contains the image node
        i_cell=0
        !
        loopk: do k=kss,km
        loopj: do j=1,jm
        loopi: do i=1,im
          !
          lin=nodeincell(cell(i,j,k),immbnod(jb)%ximag)
          !
          if(lin) then
            !
            i_cell(1)=i+ig0
            i_cell(2)=j+jg0
            i_cell(3)=k+kg0
            !
            exit loopk
            !
          endif
          !
        enddo loopi
        enddo loopj
        enddo loopk
        !
        i_cell=psum(i_cell)
        immbnod(jb)%icell=i_cell
        !
        ! to check if icell contain a ghost node
        immbnod(jb)%icell_bnode=0
        immbnod(jb)%icell_ijk=-1
        !
        ncou=0
        do jj=-1,0
        do ii=-1,0
          !
          ncou=ncou+1
          !
          i=immbnod(jb)%icell(1)-ig0+ii
          j=immbnod(jb)%icell(2)-jg0+jj
          k=immbnod(jb)%icell(3)-kg0
          !
          immbnod(jb)%icell_ijk(ncou,1)=i
          immbnod(jb)%icell_ijk(ncou,2)=j
          immbnod(jb)%icell_ijk(ncou,3)=k
          !
          if(i>=0 .and. i<=im .and. j>=0 .and. j<=jm) then
            !
            if(nodestat(i,j,k)>0) then
              ! icell contain a solid node or a boundary node
              !
              do kb=1,size(immbnod)
                !
                if( immbnod(kb)%igh(1)==i+ig0 .and. &
                    immbnod(kb)%igh(2)==j+jg0 .and. &
                    immbnod(kb)%igh(3)==k+kg0 ) then
                  !
                  immbnod(jb)%icell_bnode(ncou)=kb
                  !
                  ! reset ijk that is effective
                  immbnod(jb)%icell_ijk(ncou,1)=-1
                  immbnod(jb)%icell_ijk(ncou,2)=-1
                  immbnod(jb)%icell_ijk(ncou,3)=-1
                  !
                  exit
                  !
                endif
                !
              enddo
              !
            endif
            !
          endif
          !
        enddo
        enddo
        !
        ! determine interpolation coefficient
        !
        i=immbnod(jb)%icell(1)-ig0
        j=immbnod(jb)%icell(2)-jg0
        !
        if(i>=1 .and. i<=im .and. j>=1 .and. j<=jm ) then
          ! icell is in the domain
          !
          do m=1,4
            !
            if(immbnod(jb)%icell_ijk(m,1)>=0) then
              i=immbnod(jb)%icell_ijk(m,1)
              j=immbnod(jb)%icell_ijk(m,2)
              k=immbnod(jb)%icell_ijk(m,3)
              !
              xcell(m,:)=x(i,j,k,:)
              !
              Tm1(m,1)=xcell(m,1)*xcell(m,2)
              Tm1(m,2)=xcell(m,1)
              Tm1(m,3)=xcell(m,2)
              Tm1(m,4)=1.d0
              !
              Tm2(m,1)=xcell(m,1)*xcell(m,2)
              Tm2(m,2)=xcell(m,1)
              Tm2(m,3)=xcell(m,2)
              Tm2(m,4)=1.d0
              !
            elseif(immbnod(jb)%icell_bnode(m)>0) then
              kb=immbnod(jb)%icell_bnode(m)
              !
              xcell(m,:)=immbnod(kb)%x(:)
              xnorm(m,:)=immbnod(kb)%normdir(:)
              !
              Tm1(m,1)=xcell(m,1)*xcell(m,2)
              Tm1(m,2)=xcell(m,1)
              Tm1(m,3)=xcell(m,2)
              Tm1(m,4)=1.d0
              !
              Tm2(m,1)=xcell(m,1)*xnorm(m,2)+xcell(m,2)*xnorm(m,1)
              Tm2(m,2)=xnorm(m,1)
              Tm2(m,3)=xnorm(m,2)
              Tm2(m,4)=0.d0
              !
            else
              stop ' !! ERROR in determining interpolation coefficient'
            endif
            !
          enddo
          !
          allocate( immbnod(jb)%coef_dirichlet(4,4),                   &
                    immbnod(jb)%coef_neumann(4,4)  )
          !
          immbnod(jb)%coef_dirichlet=matinv4(Tm1)
          immbnod(jb)%coef_neumann  =matinv4(Tm2)
          !
          ! Ti1=matmul(Tm1,immbnod(jb)%coef_dirichlet)
          ! Ti2=matmul(Tm2,immbnod(jb)%coef_neumann)
          ! !
          ! write(*,"(I0,A,4(1X,F16.12))")mpirank,'|',Ti2(1,:)
          ! write(*,"(I0,A,4(1X,F16.12))")mpirank,'|',Ti2(2,:)
          ! write(*,"(I0,A,4(1X,F16.12))")mpirank,'|',Ti2(3,:)
          ! write(*,"(I0,A,4(1X,F16.12))")mpirank,'|',Ti2(4,:)
          ! print*,'---------------------------------------------------------'
          !
        endif
        !
        ! setting locality
        immbnod(jb)%localin=.false.
        !
        i=immbnod(jb)%igh(1)-ig0
        j=immbnod(jb)%igh(2)-jg0
        k=immbnod(jb)%igh(3)-kg0
        !
        if(i>=0 .and. i<=im .and. &
           j>=0 .and. j<=jm ) then
          immbnod(jb)%localin=.true.
          !
        endif
        !
        ! checking icell 
        i=immbnod(jb)%icell(1)-ig0
        j=immbnod(jb)%icell(2)-jg0
        k=immbnod(jb)%icell(3)-kg0
        !
        if(i>=1 .and. i<=im .and. &
           j>=1 .and. j<=jm ) then
          immbnod(jb)%localin=.true.
        endif
        !
      enddo
      !
    else
      stop ' !! ERROR DIMENSION NOT SET @ gridinsolid'
    endif
    !
    ! jb=1098
    ! print*,mpirank,'|',immbnod(jb)%igh(1)-ig0,immbnod(jb)%igh(2)-jg0
    !
    ! search the supporting points
    ! if(ndims==2) then
    !   !
    !   do jb=1,size(immbnod)
    !     !
    !     immbnod(jb)%dis_imga_inmg=1.d10
    !     !
    !     do k=0,ke1
    !     do j=0,jm-1
    !     do i=0,im-1
    !       !
    !       if(any(nodestat(i:i+1,j:j+1,k)==10)) cycle
    !       !
    !       do jdir=1,2
    !         xmin(jdir)=min(x(i,  j,k,jdir),x(i,  j+1,k,jdir),          &
    !                        x(i+1,j,k,jdir),x(i+1,j+1,k,jdir)  )
    !         xmax(jdir)=max(x(i,  j,k,jdir),x(i,  j+1,k,jdir),          &
    !                        x(i+1,j,k,jdir),x(i+1,j+1,k,jdir) )
    !       enddo
    !       !
    !       if( immbnod(jb)%ximag(1)>=xmin(1) .and.                      &
    !           immbnod(jb)%ximag(1)<=xmax(1) .and.                      &
    !           immbnod(jb)%ximag(2)>=xmin(2) .and.                      &
    !           immbnod(jb)%ximag(2)<=xmax(2) ) then
    !         !
    !         do jj=0,1
    !         do ii=0,1
    !           !
    !           if(nodestat(i+ii,j+jj,k)>1) cycle ! keep the boundary node
    !           !
    !           dist=dis2point(immbnod(jb)%ximag,x(i+ii,j+jj,k,:))
    !           !
    !           if(dist<=immbnod(jb)%dis_imga_inmg) then
    !             !
    !             immbnod(jb)%dis_imga_inmg=dist
    !             !
    !             immbnod(jb)%inmg(1)=ig0+i+ii
    !             immbnod(jb)%inmg(2)=jg0+j+jj
    !             immbnod(jb)%inmg(3)=kg0+k
    !             !
    !           endif
    !           !
    !         enddo
    !         enddo
    !         !
    !       endif
    !       !
    !     enddo
    !     enddo
    !     enddo
    !     !
    !     ! to get real inmg
    !     call syncinmg(immbnod(jb))
    !     !
    !     ! checking the distance between the ximag and inmg
    !     i=immbnod(jb)%inmg(1)-ig0
    !     j=immbnod(jb)%inmg(2)-jg0
    !     k=immbnod(jb)%inmg(3)-kg0
    !     !
    !     if(i>=0 .and. i<=im .and. j>=0 .and. j<=jm .and. k>=0 .and. k<=km) then
    !       var1=dis2point(x(i,j,k,:),immbnod(jb)%ximag)
    !       !
    !       ! print*,var1,'vs',immbnod(jb)%dis_imga_inmg
    !       if(var1>dxyzmax) then
    !         print*,mpirank,'|',jb
    !         print*,var1,immbnod(jb)%dis_imga_inmg<1.d5
    !         print*,' !! distance too large between the ximag and x(inmg) !!'
    !         stop
    !       endif
    !       !
    !     endif
    !     !
    !     ! build supporting nodes
    !     ncou=0
    !     do ii=-2,2
    !     do jj=-2,2
    !       !
    !       i=immbnod(jb)%inmg(1)-ig0+ii
    !       j=immbnod(jb)%inmg(2)-jg0+jj
    !       k=immbnod(jb)%inmg(3)-kg0
    !       !
    !       if(i<0 .or. i>im .or. j<0 .or. j>jm) cycle
    !       ! only search for the nodes in the domain, not dummy nodes
    !       !
    !       if(nodestat(i,j,k)==0) then
    !         ! supporting nodes only in the fluids domain
    !         !
    !         if(dis2point2(x(i,j,k,:),immbnod(jb)%x)<epsilon) then
    !           ! for the supporting nodes is too close to the boundary
    !           cycle
    !           !
    !         endif
    !         !
    !         ncou=ncou+1
    !         !
    !         snodes(ncou,1)=i+ig0
    !         snodes(ncou,2)=j+jg0
    !         snodes(ncou,3)=k+kg0
    !         !
    !       endif
    !       !
    !     enddo
    !     enddo
    !     !
    !     ! to get assumble supporting nodes
    !     call syncisup(immbnod(jb),snodes,ncou)
    !     !
    !     ! if(mpirank==0) then
    !     !   ! print supporting nodes.
    !     !   m=size(immbnod(jb)%isup,1)
    !     !   write(*,'(I0,A,2(E15.7E3),A,9(I0,A,I0,1X))')jb,'|',          &
    !     !                                    immbnod(jb)%ximag(1:2),'|', &
    !     !       (immbnod(jb)%isup(jx,1),'-',immbnod(jb)%isup(jx,2),jx=1,m)
    !     ! endif
    !     !
    !     ! to calculate weight function from supporting nodes
    !     m=size(immbnod(jb)%isup,1)
    !     !
    !     allocate(immbnod(jb)%weig(m))
    !     !
    !     immbnod(jb)%weig=0.d0
    !     !
    !     do jx=1,m
    !       !
    !       i=immbnod(jb)%isup(jx,1)-ig0
    !       j=immbnod(jb)%isup(jx,2)-jg0
    !       k=immbnod(jb)%isup(jx,3)-kg0
    !       !
    !       if(i>=iss .and. i<=im .and. &
    !          j>=jss .and. j<=jm .and. &
    !          k>=kss .and. k<=km ) then
    !         !
    !         immbnod(jb)%weig(jx)=dis2point2(immbnod(jb)%ximag,x(i,j,k,:))
    !         !
    !       endif
    !       !
    !     enddo
    !     !
    !     ! to get the weight using distance
    !     call syncweig(immbnod(jb))
    !     !
      ! do jb=1,size(immbnod)
      !   !
      !   ! to get the localisation preoperty
      !   immbnod(jb)%localin=.false.
      !   !
      !   ! checking ghost nodes
      !   i=immbnod(jb)%igh(1)-ig0
      !   j=immbnod(jb)%igh(2)-jg0
      !   k=immbnod(jb)%igh(3)-kg0
      !   !
      !   if(i>=0 .and. i<=im .and. &
      !      j>=0 .and. j<=jm .and. &
      !      k>=0 .and. k<=km  ) then
      !     immbnod(jb)%localin=.true.
      !     !
      !   endif
      !   !
      !   ! checking icell 
      !   i=immbnod(jb)%icell(1)-ig0
      !   j=immbnod(jb)%icell(2)-jg0
      !   k=immbnod(jb)%icell(3)-kg0
      !   !
      !   if(i>=0 .and. i<=im .and. &
      !      j>=0 .and. j<=jm .and. &
      !      k>=0 .and. k<=km  ) then
      !     immbnod(jb)%localin=.true.
      !   endif
        !
        ! ! checking supporting nodes
        ! do jx=1,m
        !   !
        !   i=immbnod(jb)%isup(jx,1)-ig0
        !   j=immbnod(jb)%isup(jx,2)-jg0
        !   k=immbnod(jb)%isup(jx,3)-kg0
        !   !
        !   if(i>=iss .and. i<=im .and. &
        !      j>=jss .and. j<=jm .and. &
        !      k>=kss .and. k<=km ) then
        !     !
        !     immbnod(jb)%localin=.true.
        !     !
        !     exit
        !     !
        !   endif
        !   !
        ! enddo
        !
      ! enddo
    !   !
    ! elseif(ndims==3) then
    !   stop ' !! 3D not setup yet @ sub. gridinsolid' 
    ! endif
    !
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
    ! fh=get_unit()
    ! open(fh,file='tecbound_f'//mpirankname//'.dat')
    ! do n=1,size(immbnod)
    !   if(immbnod(n)%localin .and. immbnod(n)%nodetype=='f') then
    !     write(fh,'(3(1X,E20.13E2))')immbnod(n)%x(:)
    !   endif
    ! enddo
    ! close(fh)
    ! print*,' << tecbound_f.dat'
    ! !
    ! fh=get_unit()
    ! open(fh,file='tecbound_g'//mpirankname//'.dat')
    ! do n=1,size(immbnod)
    !   if(immbnod(n)%localin .and. immbnod(n)%nodetype=='g') then
    !     write(fh,'(3(1X,E20.13E2))')immbnod(n)%x(:)
    !   endif
    ! enddo
    ! close(fh)
    ! print*,' << tecbound_g.dat'
    ! !
    ! fh=get_unit()
    ! open(fh,file='tecimag_f'//mpirankname//'.dat')
    ! do n=1,size(immbnod)
    !   if(immbnod(n)%localin .and. immbnod(n)%nodetype=='f') then
    !     write(fh,'(3(1X,E20.13E2))')immbnod(n)%ximag(:)
    !   endif
    ! enddo
    ! close(fh)
    ! print*,' << tecimag_f.dat'
    ! fh=get_unit()
    ! open(fh,file='tecimag_g'//mpirankname//'.dat')
    ! do n=1,size(immbnod)
    !   if(immbnod(n)%localin .and. immbnod(n)%nodetype=='g') then
    !     write(fh,'(3(1X,E20.13E2))')immbnod(n)%ximag(:)
    !   endif
    ! enddo
    ! close(fh)
    ! print*,' << tecimag_g.dat'
    ! !
    ! fh=get_unit()
    ! open(fh,file='tecimagijk_f'//mpirankname//'.dat')
    ! do n=1,size(immbnod)
    !   if(immbnod(n)%localin .and. immbnod(n)%nodetype=='f') then
    !     i=immbnod(n)%inmg(1)-ig0
    !     j=immbnod(n)%inmg(2)-jg0
    !     k=immbnod(n)%inmg(3)-kg0
    !     write(fh,'(3(1X,E20.13E2))')x(i,j,k,:)
    !   endif
    ! enddo
    ! close(fh)
    ! print*,' << tecimagijk_f.dat'
    ! !
    ! fh=get_unit()
    ! open(fh,file='tecimagijk_g'//mpirankname//'.dat')
    ! do n=1,size(immbnod)
    !   if(immbnod(n)%localin .and. immbnod(n)%nodetype=='g') then
    !     i=immbnod(n)%inmg(1)-ig0
    !     j=immbnod(n)%inmg(2)-jg0
    !     k=immbnod(n)%inmg(3)-kg0
    !     write(fh,'(3(1X,E20.13E2))')x(i,j,k,:)
    !   endif
    ! enddo
    ! close(fh)
    ! print*,' << tecimagijk_g.dat'
    ! !
    ! fh=get_unit()
    ! open(fh,file='tecgho'//mpirankname//'.dat')
    ! do n=1,size(immbnod)
    !   if(immbnod(n)%localin .and. immbnod(n)%nodetype=='g' ) then
    !     i=immbnod(n)%igh(1)-ig0
    !     j=immbnod(n)%igh(2)-jg0
    !     k=immbnod(n)%igh(3)-kg0
    !     write(fh,'(3(1X,E20.13E2))')x(i,j,k,:)
    !     if(x(i,j,k,1)>5.6d0) print*,mpirank,'|',n,i,j,k
    !   endif
    ! enddo
    ! close(fh)
    ! print*,' << tecgho.dat'
    ! !
    ! fh=get_unit()
    ! open(fh,file='tecforce'//mpirankname//'.dat')
    ! do n=1,size(immbnod)
    !   if(immbnod(n)%localin .and. immbnod(n)%nodetype=='f' ) then
    !     i=immbnod(n)%igh(1)-ig0
    !     j=immbnod(n)%igh(2)-jg0
    !     k=immbnod(n)%igh(3)-kg0
    !     write(fh,'(3(1X,E20.13E2))')x(i,j,k,:)
    !     if(x(i,j,k,1)>5.6d0) print*,mpirank,'|',n,i,j,k
    !     if(x(i,j,k,2)>5.6d0) print*,mpirank,'|',n,i,j,k
    !   endif
    ! enddo
    ! close(fh)
    ! print*,' << tecforce.dat'
    ! !
    ! fh=get_unit()
    ! open(fh,file='tecsupp_f'//mpirankname//'.dat')
    ! do n=1,size(immbnod)
    !   if(immbnod(n)%localin .and. immbnod(n)%nodetype=='f' ) then
    !     do jb=1,size(immbnod(n)%isup,1)
    !       i=immbnod(n)%isup(jb,1)-ig0
    !       j=immbnod(n)%isup(jb,2)-jg0
    !       k=immbnod(n)%isup(jb,3)-kg0
    !       write(fh,'(3(1X,E20.13E2))')x(i,j,k,:)
    !     enddo
    !   endif
    ! enddo
    ! close(fh)
    ! print*,' << tecsupp_f.dat'
    ! !
    ! fh=get_unit()
    ! open(fh,file='tecsupp_g'//mpirankname//'.dat')
    ! do n=1,size(immbnod)
    !   if(immbnod(n)%localin .and. immbnod(n)%nodetype=='g' ) then
    !     do jb=1,size(immbnod(n)%isup,1)
    !       i=immbnod(n)%isup(jb,1)-ig0
    !       j=immbnod(n)%isup(jb,2)-jg0
    !       k=immbnod(n)%isup(jb,3)-kg0
    !       write(fh,'(3(1X,E20.13E2))')x(i,j,k,:)
    !     enddo
    !   endif
    ! enddo
    ! close(fh)
    ! print*,' << tecsupp_g.dat'
    !
    ! fh=get_unit()
    ! open(fh,file='tec_icellin_g'//mpirankname//'.dat')
    ! do n=1,size(immbnod)
    !   i=immbnod(n)%icell(1)-ig0
    !   j=immbnod(n)%icell(2)-jg0
    !   k=immbnod(n)%icell(3)-kg0
    !   if(ijkin(i,j,k) .and. ijkin(i-1,j-1,k))  then
    !     !
    !     do m=1,4
    !       !
    !       if(immbnod(n)%icell_ijk(m,1)>=0) then
    !         i=immbnod(n)%icell_ijk(m,1)
    !         j=immbnod(n)%icell_ijk(m,2)
    !         k=immbnod(n)%icell_ijk(m,3)
    !         xcell(m,:)=x(i,j,k,:)
    !       elseif(immbnod(n)%icell_bnode(m)>0) then
    !         n1=immbnod(n)%icell_bnode(m)
    !         xcell(m,:)=immbnod(n1)%x(:)
    !       endif
    !       !
    !     enddo
    !     !
    !     write(fh,'(A)')'TITLE = "FE-Volume QUADRILATERAL Data"'
    !     write(fh,'(a)')'variables = "x", "y", "z"'
    !     write(fh,'(A,I0,A)')'ZONE T="P_',n,'", DATAPACKING=BLOCK, NODES=4, ELEMENTS= 1, ZONETYPE=FEQUADRILATERAL'

    !     write(fh,'(4(1X,E20.13E2))')xcell(1,1),xcell(2,1),xcell(4,1),xcell(3,1)
    !     write(fh,'(4(1X,E20.13E2))')xcell(1,2),xcell(2,2),xcell(4,2),xcell(3,2)
    !     write(fh,'(4(1X,E20.13E2))')xcell(1,3),xcell(2,3),xcell(4,3),xcell(3,3)
    !     write(fh,'(A)')'1,2,3,4'
    !   endif
    ! enddo
    ! close(fh)
    ! print*,' << tec_icellin_g.dat'
    !
    deallocate(rnodestat,bnodes,marker)
    !
    ! call mpistop
    !
    if(lio) print*,' ** grid in solid calculated.'
    !
    return
    !
  end subroutine gridinsolid
  !+-------------------------------------------------------------------+
  !| The end of the subroutine solidingrid.                            |
  !+-------------------------------------------------------------------+
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
    use commtype, only : nodcel
    use commvar,  only : ndims
    !
    ! arguments
    type(nodcel),intent(in) :: acell
    real(8),intent(in) :: p(3)
    logical :: lin
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
        stop ' !! ERROR 1 @ solidrange'
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
  subroutine solidrange(asolid)
    !
    use commtype,  only : solid,triangle
    use commvar,   only : xmax,xmin,ymax,ymin,zmax,zmin
    use commfunc,  only : areatriangle,cross_product
    !
    ! arguments
    type(solid),intent(inout) :: asolid
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
    write(*,'(2(A,E15.7E3))')'    x: ',asolid%xmin(1),'~',asolid%xmax(1)
    write(*,'(2(A,E15.7E3))')'    y: ',asolid%xmin(2),'~',asolid%xmax(2)
    write(*,'(2(A,E15.7E3))')'    z: ',asolid%xmin(3),'~',asolid%xmax(3)
    write(*,'(A,3(E15.7E3))')' cent: ',asolid%xcen(:)
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
    asolid%xmin=asolid%xmin*scale
    asolid%xmax=asolid%xmax*scale
    asolid%xcen=asolid%xcen*scale
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
  subroutine gridgeom
    !
    use commvar,   only : ia,ja,ka,hm,npdci,npdcj,npdck,               &
                          xmax,xmin,ymax,ymin,zmax,zmin,voldom,difschm,&
                          dxyzmax,dxyzmin
    use commarray, only : x,jacob,dxi,cell
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
        cell(i,j,k)%x(3,:)=x(i,  j-1,k,:)
        cell(i,j,k)%x(4,:)=x(i-1,j-1,k,:)
        cell(i,j,k)%x(5,:)=x(i-1,j,  k-1,:)
        cell(i,j,k)%x(6,:)=x(i,  j,  k-1,:)
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
  !| This subroutine is to decide is a line from a point intersects    |
  !| with a triangle.                                                  |
  !+-------------------------------------------------------------------+
  !| ref: https://zh.wikipedia.org/wiki/%E7%BA%BF%E9%9D%A2%E4%BA%A4%E7%82%B9
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 04-07-2021  | Created by J. Fang @ Warrington                     |
  !+-------------------------------------------------------------------+
  subroutine ray2triangle(tria,p,intertri,pintersection)
    !
    use commtype, only :  triangle
    use commfunc,  only : dis2point,areatriangle
    !
    ! arguments
    type(triangle),intent(in) :: tria
    real(8),intent(in) :: p(3)
    logical,intent(out) :: intertri
    real(8),intent(out) :: pintersection(3)
    !
    ! local data
    integer :: i,j
    real(8) :: slop(3),vec1(3),vec2(3),intpoint(3),d,d2,ldn
    real(8) :: epsilon=1.d-10
    real(8) :: trxmin(3),trxmax(3)
    !
    ! establish a random slop
    slop=(/1.d0,0.d0,0.d0/)
    !
    ldn=dot_product(slop,tria%normdir)
    ! vec1=num1d3*(tria%a+tria%b+tria%c)-p
    vec1=tria%a-p
    d=dot_product(vec1,tria%normdir)
    !
    !
    if(abs(ldn)<epsilon) then
      ! the ray line is parallel to the plane
      !
      if(abs(d)<epsilon) then
        ! the line is in the plane
        intertri=pointintriangle(tria,p)
        if(intertri) pintersection=intpoint
        !
      else
        intertri=.false.
        pintersection=(/1.d10,1.d10,1.d10/)
      endif
      !
    else
      !
      d=d/ldn
      intpoint=d*slop+p
      !
      do i=1,3
        trxmin(i)=min(tria%a(i),tria%b(i),tria%c(i))
        trxmax(i)=max(tria%a(i),tria%b(i),tria%c(i))
      enddo
      !
      ! use the ray line condition
      vec2=intpoint-p
      d2=dot_product(vec2,slop)
      !
      if(d2>=0 .or. dis2point(intpoint,p)<epsilon) then
        !
        intertri=pointintriangle(tria,intpoint)
        if(intertri)  pintersection=intpoint
        !
        ! if(intertri) then
        !   print*,tria%normdir,':',intpoint
        !   ! print*,tria%a,'-',tria%b,'-',tria%c
        ! endif
        ! if(intpoint(1)>=trxmin(1) .and. intpoint(1)<=trxmax(1) .and.     &
        !    intpoint(2)>=trxmin(2) .and. intpoint(2)<=trxmax(2) .and.     & 
        !    intpoint(3)>=trxmin(3) .and. intpoint(3)<=trxmax(3) ) then
        !   !
        !   intertri=pointintriangle(tria,intpoint)
          write(*,'(A)')'----------------------------------------------'
          write(*,'(3(1X,E20.12E3))')tria%normdir
          write(*,'(3(1X,E20.12E3))')tria%a
          write(*,'(3(1X,E20.12E3))')tria%b
          write(*,'(3(1X,E20.12E3))')tria%c
          write(*,'(A)')'----------------------------------------------'
          write(*,'(3(1X,E20.12E3))')intpoint
        ! !   !
        ! !   print*,areatriangle(tria%a,tria%b,intpoint) + &
        ! !          areatriangle(tria%a,tria%c,intpoint) + &
        ! !          areatriangle(tria%b,tria%c,intpoint),tria%area,intertri
        !   !
        !   ! open(18,file='test.plt',position='append')
        !   ! write(18,'(a)')'variables = "x", "y", "z"'
        !   ! write(18,'(2(A,I0),A)')'ZONE T="P_1", DATAPACKING=POINT, NODES=', &
        !   !                      4,', ELEMENTS=',4,', ZONETYPE=FETRIANGLE'
        !   ! write(18,'(3(1X,E15.7E3))')tria%a
        !   ! write(18,'(3(1X,E15.7E3))')tria%b
        !   ! write(18,'(3(1X,E15.7E3))')tria%c
        !   ! write(18,'(3(1X,E15.7E3))')intpoint
        !   ! write(18,'(3(1X,I8))')1,2,3
        !   ! write(18,'(3(1X,I8))')1,2,4
        !   ! write(18,'(3(1X,I8))')1,3,4
        !   ! write(18,'(3(1X,I8))')2,3,4
        !   ! !
        !   ! close(18)
        !   ! print*,' << test.plt'
        !   !
        ! endif
        !
      else
        ! the point p is opposite of the ray's direction
        intertri=.false.
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
    logical b
    integer ( kind = 4 ) i
    real ( kind = 8 ) t
  
    b = .false.

    do i = 1, asolid%num_edge

      if(asolid%edge(i)%b(2)<p(2) .eqv. p(2)<=asolid%edge(i)%a(2)) then
        t=p(1)-asolid%edge(i)%a(1)-(p(2)-asolid%edge(i)%a(2))*         &
                   (asolid%edge(i)%b(1)-p(1))/(asolid%edge(i)%b(2)-p(2))
        if ( t < 0.0D+00 ) then
          b = .not. b
        end if
      end if
    end do

    inside = b

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
  function polyhedron_contains_point_3d ( asolid, p ) result (inside)
    !
    use commtype, only : solid
    !
    ! arguments
    type(solid),intent(in) :: asolid
    real(kind=8),intent(in) :: p(3)
    !
    ! local data
    integer, parameter :: dim_num = 3, face_order_max = 3
    !
    real(kind=8) area
    integer jface
    integer face_order
    logical inside
    real(kind=8), parameter :: pi = 3.141592653589793D+00
    real(kind=8) solid_angle
    real(kind=8) v_face(dim_num,face_order_max)
  
    area = 0.0D+00
    face_order=3
  
    do jface = 1, asolid%num_face
  
      v_face(:,1) = asolid%face(jface)%a
      v_face(:,2) = asolid%face(jface)%b
      v_face(:,3) = asolid%face(jface)%c
  
      call geo_polygon_solid_angle_3d(face_order,v_face,p,solid_angle)
  
      area = area + solid_angle
  
    end do
    !
    !  AREA should be -4*PI, 0, or 4*PI.
    !  So this test should be quite safe!
    !
    if ( area < -2.0D+00 * pi .or. 2.0D+00 * pi < area ) then
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