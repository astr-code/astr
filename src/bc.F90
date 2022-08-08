!+---------------------------------------------------------------------+
!| This module contains subroutines of applying boundary conditions.   |
!| ==============                                                      |
!| CHANGE RECORD                                                       |
!| -------------                                                       |
!| 13-02-2021: Created by J. Fang @ Warrington                         |
!+---------------------------------------------------------------------+
module bc
  !
  use constdef
  use parallel,only: lio,mpistop,mpirank,mpirankname,irk,jrk,krk,      &
                     irkm,jrkm,krkm,pmax,ptime
  use commvar, only: hm,im,jm,km,uinf,vinf,winf,pinf,roinf,tinf,ndims, &
                     num_species,flowtype,gamma,numq,npdci,npdcj,      &
                     npdck,is,ie,js,je,ks,ke,xmin,xmax,ymin,ymax,      &
                     zmin,zmax,time,num_modequ
  use commarray, only : x,dxi,jacob,prs,vel,tmp,rho,spc,q,qrhs
  use tecio
  use stlaio,  only: get_unit
  !
  implicit none
  !
  integer :: ninflowslice
  real(8) :: pout
  real(8) :: uinf_j0,uinf_jm
  integer :: bctype(6)
  real(8) :: twall(6),xrhjump,angshk,xslip
  character(len=4) :: turbinf
  !+---------------------+---------------------------------------------+
  !|              bctype | define type of boundary condition.          |
  !|               twall | wall temperature.                           |
  !|             turbinf | method of generating inflow turbulence.     |
  !+---------------------+---------------------------------------------+
  real(8),allocatable :: rho_in(:,:),vel_in(:,:,:),tmp_in(:,:),        &
                         prs_in(:,:),spc_in(:,:,:)
  real(8),allocatable :: rho_prof(:),vel_prof(:,:),tmp_prof(:),        &
                         prs_prof(:),spc_prof(:,:)
  real(8),allocatable :: rho_far(:),vel_far(:,:),tmp_far(:),           & 
                         prs_far(:),spc_far(:,:)
  real(8),allocatable,dimension(:,:,:) :: bvec_i0,bvec_im,             &
                                          bvec_j0,bvec_jm,             &
                                          bvec_k0,bvec_km
  real(8),allocatable,target :: flowvarins(:,:,:,:),timeins(:)
  !
  contains
  !
  ! Table of values for parameter (ibc(ilat))
  ! ilat = 1,6
  !
  !     ____________           
  !    /|         /|           
  !   / |  4     / |           
  !  /  |    5  /  |  j        
  ! /__________/   |  |        
  ! | 1 |______|_2 |  |____ i  
  ! |  / 6     |  /  /         
  ! | /    3   | /  /          
  ! |/         |/  k           
  ! /----------/
  !+-------------------------------------------------------------------+
  !| This subroutine is to allocate inflow array.                      |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 24-02-2021: Created by J. Fang @ Warrington                       |
  !+-------------------------------------------------------------------+
  subroutine alloinflow(ndir)
    !
    ! arguments
    integer,intent(in) :: ndir
    !
    if(ndir==1) then
      allocate( rho_in(0:jm,0:km),tmp_in(0:jm,0:km),prs_in(0:jm,0:km), &
                vel_in(0:jm,0:km,1:3),spc_in(0:jm,0:km,1:num_species)  )
    elseif(ndir==3) then
      allocate( rho_in(0:im,0:km),vel_in(0:im,0:km,1:3),               &
                tmp_in(0:im,0:km),prs_in(0:im,0:km),                   &
                spc_in(0:im,0:km,1:num_species) )
    endif
    !
  end subroutine alloinflow
  !+-------------------------------------------------------------------+
  !| The end of the subroutine alloinflow.                             |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to apply bounday conditions to geometrical     |
  !| variables.                                                        |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 14-02-2021: Created by J. Fang @ Warrington                       |
  !+-------------------------------------------------------------------+
  subroutine geombc
    !
    use commfunc,  only : gradextrp
    !
    ! local data
    integer :: i,j,k,l
    real(8) :: fex(4,10)
    !
    if(jrk==0) then
      !
      j=0
      !
      if(bctype(3)==41) then
        do k=0,km
        do i=0,im
          !
          fex(:,1) =gradextrp(qbou=dxi(i,0,k,1,1),q1st=dxi(i,1,k,1,1))
          fex(:,2) =gradextrp(qbou=dxi(i,0,k,1,2),q1st=dxi(i,1,k,1,2))
          fex(:,3) =gradextrp(qbou=dxi(i,0,k,1,3),q1st=dxi(i,1,k,1,3))
          fex(:,4) =gradextrp(qbou=dxi(i,0,k,2,1),q1st=dxi(i,1,k,2,1))
          fex(:,5) =gradextrp(qbou=dxi(i,0,k,2,2),q1st=dxi(i,1,k,2,2))
          fex(:,6) =gradextrp(qbou=dxi(i,0,k,2,3),q1st=dxi(i,1,k,2,3))
          fex(:,7) =gradextrp(qbou=dxi(i,0,k,3,1),q1st=dxi(i,1,k,3,1))
          fex(:,8) =gradextrp(qbou=dxi(i,0,k,3,2),q1st=dxi(i,1,k,3,2))
          fex(:,9) =gradextrp(qbou=dxi(i,0,k,3,3),q1st=dxi(i,1,k,3,3))
          fex(:,10)=gradextrp(qbou=jacob(i,0,k),  q1st=jacob(i,1,k))
          !
          do l=1,4
            dxi(i,-l,k,1,1)=fex(l,1) 
            dxi(i,-l,k,1,2)=fex(l,2) 
            dxi(i,-l,k,1,3)=fex(l,3) 
            dxi(i,-l,k,2,1)=fex(l,4) 
            dxi(i,-l,k,2,2)=fex(l,5) 
            dxi(i,-l,k,2,3)=fex(l,6) 
            dxi(i,-l,k,3,1)=fex(l,7) 
            dxi(i,-l,k,3,2)=fex(l,8) 
            dxi(i,-l,k,3,3)=fex(l,9) 
            jacob(i,-l,k)  =fex(l,10)
          enddo
          !
        enddo
        enddo
      endif
      !
    endif
    !
    if(jrk==jrkm) then
      !
      j=0
      !
      if(bctype(4)==41) then
        do k=0,km
        do i=0,im
          !
          dxi(i,jm+1:jm+4,k,1,1)=gradextrp(qbou=dxi(i,jm,k,1,1),q1st=dxi(i,jm-1,k,1,1))
          dxi(i,jm+1:jm+4,k,1,2)=gradextrp(qbou=dxi(i,jm,k,1,2),q1st=dxi(i,jm-1,k,1,2))
          dxi(i,jm+1:jm+4,k,1,3)=gradextrp(qbou=dxi(i,jm,k,1,3),q1st=dxi(i,jm-1,k,1,3))
          dxi(i,jm+1:jm+4,k,2,1)=gradextrp(qbou=dxi(i,jm,k,2,1),q1st=dxi(i,jm-1,k,2,1))
          dxi(i,jm+1:jm+4,k,2,2)=gradextrp(qbou=dxi(i,jm,k,2,2),q1st=dxi(i,jm-1,k,2,2))
          dxi(i,jm+1:jm+4,k,2,3)=gradextrp(qbou=dxi(i,jm,k,2,3),q1st=dxi(i,jm-1,k,2,3))
          dxi(i,jm+1:jm+4,k,3,1)=gradextrp(qbou=dxi(i,jm,k,3,1),q1st=dxi(i,jm-1,k,3,1))
          dxi(i,jm+1:jm+4,k,3,2)=gradextrp(qbou=dxi(i,jm,k,3,2),q1st=dxi(i,jm-1,k,3,2))
          dxi(i,jm+1:jm+4,k,3,3)=gradextrp(qbou=dxi(i,jm,k,3,3),q1st=dxi(i,jm-1,k,3,3))
          jacob(1,jm+1:jm+4,k)  =gradextrp(qbou=jacob(i,jm,k),  q1st=jacob(i,jm-1,k))
          !
        enddo
        enddo
      endif
      !
    endif
    !
  end subroutine geombc
  !
  subroutine xyzbc
    !
    use commvar,  only : npdci,npdcj,npdck
    use commfunc, only : gradextrp
    !
    ! local data
    integer :: i,j,k,l
    real(8) :: fex(4,3)
    !
    if(npdci==1) then
      do k=0,km
      do j=0,jm
        !
        fex(:,1) =gradextrp(qbou=x(0,j,k,1),q1st=x(1,j,k,1))
        fex(:,2) =gradextrp(qbou=x(0,j,k,2),q1st=x(1,j,k,2))
        fex(:,3) =gradextrp(qbou=x(0,j,k,3),q1st=x(1,j,k,3))
        !
        do l=1,4
          x(-l,j,k,1)=fex(l,1) 
          x(-l,j,k,2)=fex(l,2) 
          x(-l,j,k,3)=fex(l,3) 
        enddo
        !
      enddo
      enddo
    elseif(npdci==2) then
      do k=0,km
      do j=0,jm
        !
        x(im+1:im+4,j,k,1)=gradextrp(qbou=x(im,j,k,1),q1st=x(im-1,j,k,1))
        x(im+1:im+4,j,k,2)=gradextrp(qbou=x(im,j,k,2),q1st=x(im-1,j,k,2))
        x(im+1:im+4,j,k,3)=gradextrp(qbou=x(im,j,k,3),q1st=x(im-1,j,k,3))
        !
      enddo
      enddo
    endif
    !
    if(npdcj==1) then
      do k=0,km
      do i=0,im
        !
        fex(:,1) =gradextrp(qbou=x(i,0,k,1),q1st=x(i,1,k,1))
        fex(:,2) =gradextrp(qbou=x(i,0,k,2),q1st=x(i,1,k,2))
        fex(:,3) =gradextrp(qbou=x(i,0,k,3),q1st=x(i,1,k,3))
        !
        do l=1,4
          x(i,-l,k,1)=fex(l,1) 
          x(i,-l,k,2)=fex(l,2) 
          x(i,-l,k,3)=fex(l,3) 
        enddo
        !
      enddo
      enddo
    elseif(npdcj==2) then
      do k=0,km
      do i=0,im
        !
        x(i,jm+1:jm+4,k,1)=gradextrp(qbou=x(i,jm,k,1),q1st=x(i,jm-1,k,1))
        x(i,jm+1:jm+4,k,2)=gradextrp(qbou=x(i,jm,k,2),q1st=x(i,jm-1,k,2))
        x(i,jm+1:jm+4,k,3)=gradextrp(qbou=x(i,jm,k,3),q1st=x(i,jm-1,k,3))
        !
      enddo
      enddo
    endif
    !
  end subroutine xyzbc
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to apply bounday conditions to node state.     |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 21-04-2022: Created by J. Fang @ Warrington                       |
  !+-------------------------------------------------------------------+
  subroutine nodebc(nodestat)
    !
    integer,intent(inout) :: nodestat(-hm:im+hm,-hm:jm+hm,-hm:km+hm)
    !
    ! local data
    integer :: n
    !
    do n=1,6
      !
      if(bctype(n)==41 .or. bctype(n)==411 .or. bctype(n)==421) then
        call solidnode(n,nodestat)
      endif
      !
    enddo
    !
  end subroutine nodebc
  !+-------------------------------------------------------------------+
  !| The end of the subroutine nodebc.                                 |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to apply solid bc to nodes.                    |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 21-04-2021: Created by J. Fang @ Warrington                       |
  !+-------------------------------------------------------------------+
  subroutine solidnode(ndir,nodestat)
    !
    ! arguments
    integer,intent(in) :: ndir
    integer,intent(inout) :: nodestat(-hm:im+hm,-hm:jm+hm,-hm:km+hm)
    !
    ! local data
    integer :: i,j,k
    !
    if(ndir==3) then
      !
      if(jrk==0) then
        !
        j=0
        !
        do k=0,km
        do i=0,im
          !
          nodestat(i,-hm:0,k)=5
          !
        enddo
        enddo
        !
      endif
      !
    elseif(ndir==4) then
      !
      if(jrk==jrkm) then
        !
        j=jm
        !
        do k=0,km
        do i=0,im
          !
          nodestat(i,jm:jm+hm,k)=5
          !
        enddo
        enddo
        !
      endif
      !
    else
      stop ' !! ndir not defined @ solidnode'
    endif
    !
  end subroutine solidnode
  !+-------------------------------------------------------------------+
  !| The end of the subroutine solidnode.                              |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to apply bounday conditions.                   |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 13-02-2021: Created by J. Fang @ Warrington                       |
  !+-------------------------------------------------------------------+
  subroutine boucon(subtime)
    !
    use commvar, only : limmbou
    !
    ! arguments
    real(8),intent(inout),optional :: subtime
    !
    ! local data
    integer :: n
    real(8) :: time_beg
    !
    do n=1,6
      !
      if(bctype(n)==12) then
        call inflow_nscbc(n)
      endif
      !
      if(bctype(n)==22) then
        call outflow_nscbc(n)
      endif
      !
      if(bctype(n)==52) then
        call farfield_nscbc(n)
      endif
      !
    enddo
    !
    do n=1,6
      !
      if(bctype(n)==41) then
        call noslip(n,twall(n))
      endif
      !
      if(bctype(n)==42) then
        call noslip_adibatic(n)
      endif
      !
      if(bctype(n)==411) then
        call slipisotwall(n,twall(n))
      endif
      !
      if(bctype(n)==421) then
        call slipadibwall(n)
      endif
      !
      if(bctype(n)==51) then
        call farfield(n)
      endif
      !
      if(bctype(n)==23) then
        call gcnscbc(n,prs_t=pinf)
      endif
      !
      if(bctype(n)==11) then
        call inflow(n)
      endif
      !
      if(bctype(n)==21) then
        call outflow(n)
      endif
      !
      if(bctype(n)==31) then
        call fixbc(n)
      endif
      !
    enddo
    !
    ! if(present(subtime)) time_beg=ptime()
    !
    ! if(present(subtime)) subtime=subtime+ptime()-time_beg
    !
    ! call mpistop
    !
    return
    !
  end subroutine boucon
  !+-------------------------------------------------------------------+
  !| The end of the subroutine boucon.                                 |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to apply a immersed boundary condition to the  |
  !| sloid body nodes.                                                 |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 30-06-2021: Created by J. Fang @ Warrington                       |
  !+-------------------------------------------------------------------+
  subroutine immbody(subtime1,subtime2,subtime3,subtime4)
    !
    use commtype,  only : sboun
    !
    use commvar,   only : pinf,immbond,num_species,              &
                          num_icell_rank,num_ighost_rank
    use commarray, only : nodestat,vel,x
    use commcal,   only : ijkin,ijkcellin
    use fludyna,   only : fvar2q,q2fvar,thermal
    use parallel,  only : ig0,jg0,kg0,npdci,npdcj,npdck,qswap,msize
    use commfunc,  only : median
    !
    ! arguments
    real(8),intent(inout),optional :: subtime1,subtime2,subtime3,subtime4
    !
    ! local data
    integer :: i,j,k,m,kb,iss,jss,kss,n,jspec,jq,jb,jc,jsup,ncube,counter
    real(8) :: var_ro,var_u(3),var_t,var_p,var_sp(num_species)
    real(8) :: vel_bou(3),prs_bou,rho_bou,tmp_bou,spc_bou(num_species)
    real(8),allocatable :: vel_icell(:,:),tmp_icell(:),prs_icell(:),c_dir(:),c_neu(:)
    real(8) :: vel_image(3),tmp_image,prs_image,rho_image,spc_image(num_species)
    real(8),allocatable :: qimag(:,:)
    real(8) :: prslmin,prslmax
    type(sboun),pointer :: pb
    !
    logical,save :: lfirstcal=.true.
    integer,allocatable,save :: nbnode(:),ngnode(:),ijkcub(:,:,:),ijkgho(:,:)
    !
    real(8) :: time_beg,time_tmp
    !
    if(present(subtime1) .or. present(subtime2) ) then
      time_beg=ptime()
    endif
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
    if(ndims==2) then
      kss=0
    else
      kss=1
    endif
    !
    ncube=2**ndims
    !
    ! go over the boundary nodes
    !
    if(lfirstcal) then
      !
      counter=0
      do jb=1,size(immbond)
        !
        pb=>immbond(jb)
        !
        if(.not. allocated(pb%qimag)) allocate(pb%qimag(numq))
        !
        pb%qimag=0.d0
        !
        i=pb%icell(1)-ig0
        j=pb%icell(2)-jg0
        k=pb%icell(3)-kg0
        !
        if( ijkcellin(i,j,k)) then
          !
          counter=counter+1
          !
        endif
        !
      enddo
      !
      allocate(ijkcub(3,ncube,counter))
      allocate(nbnode(counter))
      !
      counter=0
      !
      do jb=1,size(immbond)
        !
        pb=>immbond(jb)
        !
        i=pb%icell(1)-ig0
        j=pb%icell(2)-jg0
        k=pb%icell(3)-kg0
        !
        if( ijkcellin(i,j,k)) then
          !
          counter=counter+1
          !
          nbnode(counter)=jb
          !
          do m=1,ncube
            if(pb%icell_ijk(m,1)>=0) then
              ijkcub(1,m,counter)=pb%icell_ijk(m,1)
              ijkcub(2,m,counter)=pb%icell_ijk(m,2)
              ijkcub(3,m,counter)=pb%icell_ijk(m,3)
            else
              stop ' !! ERROR in determining pb%icell_ijk @ immbody'
            endif
          enddo
        endif
        !
      enddo
      !
      counter=0
      do jb=1,size(immbond)
        !
        pb=>immbond(jb)
        !
        i=pb%igh(1)-ig0
        j=pb%igh(2)-jg0
        k=pb%igh(3)-kg0
        !
        if(i>=0 .and. i<=im .and. &
           j>=0 .and. j<=jm .and. &
           k>=0 .and. k<=km ) then
          !
          counter=counter+1
          !
        endif
        !
      enddo
      !
      allocate(ngnode(counter),ijkgho(3,counter))
      !
      counter=0
      do jb=1,size(immbond)
        !
        pb=>immbond(jb)
        !
        i=pb%igh(1)-ig0
        j=pb%igh(2)-jg0
        k=pb%igh(3)-kg0
        !
        if(i>=0 .and. i<=im .and. &
           j>=0 .and. j<=jm .and. &
           k>=0 .and. k<=km ) then
          !
          counter=counter+1
          !
          ngnode(counter)=jb
          !
          ijkgho(1,counter)=i
          ijkgho(2,counter)=j
          ijkgho(3,counter)=k
          !
        endif
        !
      enddo
      !
      lfirstcal=.false.
      !
    endif
    !
    allocate( vel_icell(ncube,3),tmp_icell(ncube),prs_icell(ncube), &
              c_dir(ncube),c_neu(ncube))
    allocate(qimag(numq,1:size(immbond)))
    !
    do jc=1,msize(nbnode)
      !
      jb=nbnode(jc)
      !
      pb=>immbond(jb)
      !
      prslmin= 1.d10
      prslmax=-1.d10
      do m=1,ncube
        !
        i=ijkcub(1,m,jc)
        j=ijkcub(2,m,jc)
        k=ijkcub(3,m,jc)
        !
        vel_icell(m,:)=vel(i,j,k,:)
        tmp_icell(m)  =tmp(i,j,k)
        prs_icell(m)  =prs(i,j,k)
        !
        prslmin=min(prslmin,prs(i,j,k))
        prslmax=max(prslmax,prs(i,j,k))
        !
      enddo
      !
      c_dir=pb%coef_dirichlet
      ! c_neu=pb%coef_neumann
      !
      vel_image(1)=dot_product(c_dir,vel_icell(:,1))
      vel_image(2)=dot_product(c_dir,vel_icell(:,2))
      vel_image(3)=dot_product(c_dir,vel_icell(:,3))
      !
      tmp_image   =dot_product(c_dir,tmp_icell)
      !
      ! prs_image   =dot_product(c_neu,prs_icell)
      prs_image   =dot_product(c_dir,prs_icell)
      prs_image   =median(prs_image,prslmin,prslmin)
      !
      rho_image=thermal(pressure=prs_image,temperature=tmp_image)
      !
      spc_image=1.d0
      !
      call fvar2q(       q=pb%qimag(:),                         &
                   density=rho_image,  velocity=vel_image,      &
                  pressure=prs_image,  species=spc_image         )
      !
    enddo
    !
    subtime2=subtime2+ptime()-time_beg
    time_tmp=ptime()
    !
    call syncqimag_nonlocal
    ! call syncqimag_nonlocal(subtime1=subtime2,subtime2=subtime3, &
    !                         subtime3=subtime4)
    !
    subtime3=subtime3+ptime()-time_tmp
    time_tmp=ptime()
    !
    ! call syncqimag(qimag(:,1:counter))
    !
    do jc=1,msize(ngnode)
      !
      jb=ngnode(jc)
      !
      pb=>immbond(jb)
      !
      i=ijkgho(1,jc)
      j=ijkgho(2,jc)
      k=ijkgho(3,jc)
      !
      call q2fvar(        q=pb%qimag,                   &
                    density=var_ro,                     &
                   velocity=var_u(:),                   &
                   pressure=var_p,                      &
                temperature=var_t                        )
      !
      if(pb%nodetype=='g') then
        ! for ghost nodes
        !
        vel(i,j,k,1)= -1.d0*var_u(1)*pb%dis2ghost/pb%dis2image
        vel(i,j,k,2)= -1.d0*var_u(2)*pb%dis2ghost/pb%dis2image
        vel(i,j,k,3)= -1.d0*var_u(3)*pb%dis2ghost/pb%dis2image
          ! tmp(i,j,k)=tinf-(var_t-tinf)*pb%dis2ghost/pb%dis2image
          tmp(i,j,k)=var_t
          prs(i,j,k)=var_p
          rho(i,j,k)=thermal(pressure=prs(i,j,k),temperature=tmp(i,j,k))
        !
        spc(i,j,k,:)=1.d0
        !
        ! print*,i,j,k,vel(i,j,k,:)
        !
      elseif(pb%nodetype=='b') then
        ! for boundary nodes
        !
        vel(i,j,k,1)=0.d0
        vel(i,j,k,2)=0.d0
        vel(i,j,k,3)=0.d0
          tmp(i,j,k)=var_t
          ! tmp(i,j,k)=tinf
          prs(i,j,k)=var_p
          rho(i,j,k)=thermal(pressure=prs(i,j,k),temperature=tmp(i,j,k))
        !
        spc(i,j,k,:)=1.d0
        !
      ! elseif(pb%nodetype=='f') then
      !   ! for force nodes
      !   vel(i,j,k,:)=0.5d0*var_u(:)
      !     tmp(i,j,k)=0.5d0*(tinf+var_t)
      !     prs(i,j,k)=var_p
      !   rho(i,j,k)=thermal(pressure=prs(i,j,k),temperature=tmp(i,j,k))
      !   !
      !   spc(i,j,k,:)=1.d0
      !   !
      ! else
      !   stop ' ERROR @ immbody'
      endif
      
      call fvar2q(      q=  q(i,j,k,:),   density=rho(i,j,k),        &
                 velocity=vel(i,j,k,:),  pressure=prs(i,j,k),        &
                  species=spc(i,j,k,:)                               )
      !
      ! if(jb==1)  then
      !   print*,mpirank,'|',jb,'-',pb%icell_rank,pb%ighst_rank
      ! endif
      !
    enddo
    !
    do k=0,km
    do j=0,jm
    do i=0,im
      !
      if(nodestat(i,j,k)==5) then
        ! inner solid
      ! if(nodestat(i,j,k)>0) then
        ! all solid nodes
        !
        vel(i,j,k,1)=0.d0
        vel(i,j,k,2)=0.d0
        vel(i,j,k,3)=0.d0
        ! tmp(i,j,k)  =twall(3)
        tmp(i,j,k)  =tinf
        prs(i,j,k)  =pinf
        rho(i,j,k)  =thermal(pressure=prs(i,j,k),temperature=tmp(i,j,k))
        !
        call fvar2q(      q=  q(i,j,k,:),   density=rho(i,j,k),        &
                   velocity=vel(i,j,k,:),  pressure=prs(i,j,k),        &
                    species=spc(i,j,k,:)                               )
      endif
      !
    enddo
    enddo
    enddo
    !
    subtime4=subtime4+ptime()-time_tmp
    if(present(subtime1)) subtime1=subtime1+ptime()-time_beg
    !
    ! call mpistop
    !
    return
    !
  end subroutine immbody
  !+-------------------------------------------------------------------+
  !| The end of the subroutine immbody.                                |
  !+-------------------------------------------------------------------+
  !!
  !+-------------------------------------------------------------------+
  !| This subroutine is used to synconize immbond(jb)%qimag for only   |
  !| nonelocal element.                                                |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 21-Jun-2022: Created by J. Fang @ STFC Daresbury Laboratory       |
  !+-------------------------------------------------------------------+
  subroutine syncqimag_nonlocal(subtime1,subtime2,subtime3)
    !
    use commtype, only : sboun
    use parallel, only : mpirankmax,pswapall,ptabupd,msize
    use commvar,  only : immbond
    !
    ! arguments
    real(8),intent(inout),optional :: subtime1,subtime2,subtime3
    !
    ! local data
    integer :: nsize,jb,js,kb,qsize,jrank,qzmax,offset,bignumber,jsend,n,nnl
    integer,allocatable :: recorder(:,:)
    integer,allocatable,save :: isendtable(:),irecvtable(:), &
                                num_bound_elemt_send(:),num_bound_elemt_recv(:)
    real(8),allocatable :: q2send(:,:),q2recv(:,:)
    type(sboun),pointer :: pbon
    !
    real(8) :: time_beg,time_tmp
    !
    integer,save :: nsend,nrecv
    !
    logical,save :: lfirstcal=.true.
    !
    if(present(subtime1) .or. present(subtime2) .or. present(subtime3) ) then
      time_beg=ptime()
    endif
    !
    if(lfirstcal) then
      !
      nnl=0
      do jb=1,size(immbond)
        !
        pbon=>immbond(jb)
        !
        if(pbon%cellin) then
          if(pbon%locality=='n') then
            do js=1,size(pbon%ighst_rank_all)
              jrank=pbon%ighst_rank_all(js)
              if(jrank==mpirank) cycle
              nnl=nnl+1
            enddo
          endif
        endif
        !
      enddo
      !
      ! print*,mpirank,'|',nnl
      !
      allocate(isendtable(0:mpirankmax),irecvtable(0:mpirankmax), &
               recorder(0:mpirankmax,nnl))
      !
      isendtable=0
      recorder=0
      nsend=0
      !
      ! call mpistop
      !
      ! first establish the sending table for nonlocal elements.
      do jb=1,size(immbond)
        !
        pbon=>immbond(jb)
        !
        if(pbon%cellin) then
          !
          if(pbon%locality=='n') then
            ! non-local element identified.
            !
            do js=1,size(pbon%ighst_rank_all)
              !
              jrank=pbon%ighst_rank_all(js)
              !
              if(jrank==mpirank) cycle
              !
              isendtable(jrank)=isendtable(jrank)+1
              !
              recorder(jrank,isendtable(jrank))=jb
              !
              nsend=nsend+1
              !
            enddo
            !
            ! if(size(pbon%ighst_rank_all)>1) then
            !   print*,mpirank,'|',pbon%ighst_rank_all
            ! endif
            !
          elseif(pbon%locality=='l') then
          else
            print*,'jb=',jb,'rank=',mpirank
            stop ' !! ERROR !! in pbon%locality @ syncqimag_nonlocal '
          endif
          !
        endif
        !
      enddo
      !
      ! print*,mpirank,'|',nsend
      !
      irecvtable=pswapall(isendtable)
      !
      if(nsend>0) then
        !
        allocate(num_bound_elemt_send(1:nsend))
        !
        n=0
        do jrank=0,mpirankmax
          !
          do jsend=1,isendtable(jrank)
            !
            n=n+1
            !
            num_bound_elemt_send(n)=recorder(jrank,jsend)
            !
              ! if(num_bound_elemt_send(n)==352) then
              !   print*,mpirank,'~',jrank,jsend
              ! endif
              !
            ! if(mpirank==1) then
            !   print*,n,num_bound_elemt_send(n),'|',mpirank,'->',jrank
            ! endif
            !
          enddo
          !
        enddo
        !
      endif
      !
      call ptabupd(num_bound_elemt_send,num_bound_elemt_recv,isendtable,irecvtable)
      !
      nrecv=msize(num_bound_elemt_recv)
      !
      if(sum(irecvtable) .ne. size(num_bound_elemt_recv)) then
        print*,mpirank,'-',sum(irecvtable),size(num_bound_elemt_recv)
        print*,' !! WARNING !! sum of  irecvtable not compatable with size of num_bound_elemt_recv'
      endif
      !
      ! if(allocated(irecvtable) .and. size(irecvtable)>0) then
      !   !
      !   n=0
      !   do jrank=0,mpirankmax
      !     do jsend=1,irecvtable(jrank)
      !       !
      !       n=n+1
      !       !
      !       if(mpirank==1) then
      !         print*,' ** rank ',mpirank,'recv element',num_bound_elemt_recv(n),' from ',jrank
      !         !
      !         jb=num_bound_elemt_recv(n)
      !         !
      !         print*,' ** image node is in rank: ',immbond(jb)%icell_rank,  &
      !                    'ghost node is in rank: ',immbond(jb)%ighst_rank
      !         !
      !       endif
      !       !
      !     enddo
      !   enddo
      !   !
      ! endif
      !
      deallocate(recorder)
      !
      lfirstcal=.false.
      !
    endif
    !
    allocate(q2send(numq,nsend),q2recv(numq,nrecv))
    !
    ! pack the data to send
    do n=1,nsend
      !
      jb=num_bound_elemt_send(n)
      !
      q2send(:,n)=immbond(jb)%qimag(:)
      !
    enddo
    !
    if(present(subtime1)) then
      subtime1=subtime1+ptime()-time_beg
      time_tmp=ptime()
    endif
    !
    ! data passing via alltoall
    call ptabupd(q2send,q2recv,isendtable,irecvtable)
    !
    if(present(subtime2)) then
      subtime2=subtime2+ptime()-time_tmp
      time_tmp=ptime()
    endif
    !
    ! unpack received data
    do n=1,nrecv
      !
      jb=num_bound_elemt_recv(n)
      immbond(jb)%qimag(:)=q2recv(:,n)
      !
    enddo
    !
    deallocate(q2send,q2recv)
    !
    if(present(subtime3)) then
      subtime3=subtime3+ptime()-time_tmp
    endif
    ! print*,mpirank,'|',isendtable
    ! print*,mpirank,'-',sum(irecvtable),size(num_bound_elemt_recv)
    !
  end subroutine syncqimag_nonlocal
  !+-------------------------------------------------------------------+
  !| The end of the subroutine syncqimag_nonlocal.                     |
  !+-------------------------------------------------------------------+
  !!
  !+-------------------------------------------------------------------+
  !| This subroutine is used to synconize immbond(jb)%qimag            |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 16-Jul-2021: Created by J. Fang @ STFC Daresbury Laboratory       |
  !+-------------------------------------------------------------------+
  subroutine syncqimag(qin)
    !
    use commtype, only : sboun
    use parallel, only : psum,ptabupd,mpirankmax,pgather,bcast,pscatter
    use commvar,  only : immbond,imb_node_have,imb_node_need,          &
                         num_ighost_rank,numq,imbroot
    !
    ! arguments
    real(8),intent(in) :: qin(:,:)
    !
    ! local data
    integer :: nsize,nsend,jb,kb,qsize,jrank,qzmax,offset
    integer :: qztable(0:mpirankmax)
    real(8),allocatable :: buf(:,:),qimmb(:,:),q2send(:,:),qboud(:,:)
    integer,allocatable :: icell_order(:)
    real(8) :: epslion=1.d-10
    !
    logical,save :: lfirstcal=.true.
    !
    nsize=size(immbond)
    nsend=size(imb_node_need)
    !
    if(lfirstcal) then
      !
      qsize=size(qin,2)
      !
      call ptabupd(qsize,qztable)
      !
      qzmax=0
      do jrank=0,mpirankmax
        if(qztable(jrank)>qzmax) then
          imbroot=jrank
          qzmax=qztable(jrank)
        endif
      enddo
      !
      if(mpirank==imbroot) then
        write(*,'(3(A,I0))')'  ** imb root rank= ',imbroot, &
                            ', local size: ',qsize,', total imb size: ',nsize
      endif
      !
      lfirstcal=.false.
      !
    endif
    !
    call pgather(qin,buf,imbroot)
    !
    allocate(qimmb(numq,nsize),q2send(numq,nsend))
    !
    if(mpirank==imbroot) then
      !
      do jb=1,nsize
        !
        kb=imb_node_have(jb)
        !
        qimmb(:,kb)=buf(:,jb)
        !
        if(.not. allocated(immbond(kb)%qimag)) allocate(immbond(kb)%qimag(numq))
        !
        immbond(kb)%qimag(:)=qimmb(:,kb)
        !
      enddo
      !
      ! pack data send to each ranks
      do jb=1,nsend
        !
        kb=imb_node_need(jb)
        !
        q2send(:,jb)=qimmb(:,kb)
        !
      enddo
      !
    endif
    !
    call pscatter(q2send,qboud,num_ighost_rank,rank=imbroot)
    !
    if(mpirank==0) then
      offset=0
    else
      offset=0
      do jrank=1,mpirank
        offset=offset+num_ighost_rank(jrank-1)
      enddo
    endif
    !
    ! print*,mpirank,'|',offset
    !
    if(mpirank.ne.imbroot) then
      !
      do jb=1,size(qboud,2)
        !
        kb=imb_node_need(jb+offset)
        !
        if(.not. allocated(immbond(kb)%qimag)) allocate(immbond(kb)%qimag(numq))
        !
        immbond(kb)%qimag(:)=qboud(:,jb)
        !
        ! if(jb==size(qboud,2)) then
        !   print*,mpirank,'|',jb,kb
        ! endif
        !
      enddo
      !
    endif
    !
    ! nsize=size(asboun)
    ! !
    ! allocate(buf(1:numq,1:nsize))
    ! !
    ! do jb=1,nsize
    !   buf(1:numq,jb)=asboun(jb)%qimag(1:numq)
    ! enddo
    ! !
    ! buf=psum(buf)
    ! !
    ! do jb=1,nsize
    !   asboun(jb)%qimag(1:numq)=buf(1:numq,jb)
    ! enddo
    ! !
    ! deallocate(buf)
    !
    return
    !
  end subroutine syncqimag
  !+-------------------------------------------------------------------+
  !| The end of the subroutine syncisup.                               |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to apply corner bc.                            |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 06-03-2021: Created by J. Fang @ Warrington                       |
  !+-------------------------------------------------------------------+
  subroutine corner4nscbc
    !
    use fludyna,   only : fvar2q,q2fvar,thermal,sos
    use commfunc,  only : extrapolate
    !
    integer :: i,j,k,n,jspec
    real(8) :: css,ue,pe,roe
    !
    ! 
    !
    !+---------+
    !| edge 1  |
    !+---------+
    if(irk==0 .and. jrk==0) then
      i=0
      j=0
      !
      do k=ks,ke
        !
        qrhs(i,j,k,:)=0.d0
        !
        css=sos(tmp(i,j,k))
        ue  =extrapolate(vel(i+1,j,k,1),vel(i+2,j,k,1),dv=0.d0)
        pe  =extrapolate(prs(i+1,j,k),  prs(i+2,j,k),  dv=0.d0)
        roe =extrapolate(rho(i+1,j,k),  rho(i+2,j,k),  dv=0.d0)
        !
        vel(i,j,k,1)=0.5d0*(prs_in(j,k)-pe)/(rho(i,j,k)*css)+0.5d0*(vel_in(j,k,1)+ue)
        vel(i,j,k,2)=vel_in(j,k,2)
        vel(i,j,k,3)=vel_in(j,k,3)
        prs(i,j,k)  =0.5d0*(prs_in(j,k)+pe)+0.5d0*rho(i,j,k)*css*(vel_in(j,k,1)-ue)
        rho(i,j,k)  =rho_in(j,k)*(prs(i,j,k)/prs_in(j,k))**(1.d0/gamma)
        !
        tmp(i,j,k)  =thermal(pressure=prs(i,j,k),density=rho(i,j,k))
        !
        spc(i,j,k,:)=spc_in(j,k,:)
        !
        call fvar2q(      q=  q(i,j,k,:),   density=rho(i,j,k),        &
                   velocity=vel(i,j,k,:),  pressure=prs(i,j,k),        &
                    species=spc(i,j,k,:)                               )
      enddo
      !
    endif
    !
    !+---------+
    !| edge 2  |
    !+---------+
    if(irk==0 .and. jrk==jrkm) then
      i=0
      j=jm
      !
      do k=ks,ke
        !
        qrhs(i,j,k,:)=0.d0
        !
        css=sos(tmp(i,j,k))
        ue  =extrapolate(vel(i+1,j,k,1),vel(i+2,j,k,1),dv=0.d0)
        pe  =extrapolate(prs(i+1,j,k),  prs(i+2,j,k),  dv=0.d0)
        roe =extrapolate(rho(i+1,j,k),  rho(i+2,j,k),  dv=0.d0)
        !
        vel(i,j,k,1)=0.5d0*(prs_in(j,k)-pe)/(rho(i,j,k)*css)+0.5d0*(vel_in(j,k,1)+ue)
        vel(i,j,k,2)=vel_in(j,k,2)
        vel(i,j,k,3)=vel_in(j,k,3)
        prs(i,j,k)  =0.5d0*(prs_in(j,k)+pe)+0.5d0*rho(i,j,k)*css*(vel_in(j,k,1)-ue)
        rho(i,j,k)  =rho_in(j,k)*(prs(i,j,k)/prs_in(j,k))**(1.d0/gamma)
        !
        tmp(i,j,k)  =thermal(pressure=prs(i,j,k),density=rho(i,j,k))
        !
        spc(i,j,k,:)=spc_in(j,k,:)
        !
        call fvar2q(      q=  q(i,j,k,:),   density=rho(i,j,k),        &
                   velocity=vel(i,j,k,:),  pressure=prs(i,j,k),        &
                    species=spc(i,j,k,:)                               )
      enddo
      !
    endif
    !
    !+---------+
    !| edge 3  |
    !+---------+
    if(irk==irkm .and. jrk==0) then
      i=im
      j=0
      !
      do k=ks,ke
        !
        qrhs(i,j,k,:)=0.d0
        !
        do n=1,numq
          q(i,j,k,n)=0.5d0*(extrapolate(q(i,j-1,k,n),q(i,j-2,k,n),0.d0)+ &
                            extrapolate(q(i,j+1,k,n),q(i,j+2,k,n),0.d0))
        enddo
        call fvar2q(          q=  q(i,j,k,:),                          &
                        density=rho(i,j,k),                            &
                       velocity=vel(i,j,k,:),                          &
                       pressure=prs(i,j,k),                            &
                        species=spc(i,j,k,:)                           )
      enddo
      !
    endif
    !
    !+---------+
    !| edge 4  |
    !+---------+
    if(irk==irkm .and. jrk==jrkm) then
      i=im
      j=jm
      !
      do k=ks,ke
        !
        qrhs(i,j,k,:)=0.d0
        !
        do n=1,numq
          q(i,j,k,n)=0.5d0*(extrapolate(q(i,j-1,k,n),q(i,j-2,k,n),0.d0)+  &
                            extrapolate(q(i,j-1,k,n),q(i,j-2,k,n),0.d0))
        enddo
        call fvar2q(          q=  q(i,j,k,:),                          &
                        density=rho(i,j,k),                            &
                       velocity=vel(i,j,k,:),                          &
                       pressure=prs(i,j,k),                            &
                        species=spc(i,j,k,:)                           )
      enddo
      !
    endif
    !
  end subroutine corner4nscbc
  !+-------------------------------------------------------------------+
  !| The end of the subroutine corner4nscbc.                           |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to apply inflow bc.                            |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 24-02-2021: Created by J. Fang @ Warrington                       |
  !+-------------------------------------------------------------------+
  subroutine inflow(ndir)
    !
    use commvar,   only : turbmode
    use fludyna,   only : thermal,fvar2q,q2fvar,sos
    use commfunc,  only : extrapolate
    use commarray, only : tke,omg
    !
    ! arguments
    integer,intent(in) :: ndir
    !
    ! local data
    integer :: i,j,k,l,jspc
    real(8) :: css,csse,ub,pe,roe,ue,pwave_out,pwave_in,malo
    !
    logical,save :: lfirstcal=.true.
    !
    if(ndir==1 .and. irk==0) then
      !
      if(lfirstcal) then
        !
        call alloinflow(ndir)
        !
        if(trim(flowtype)=='jet') then
          call jetinflow
        elseif(trim(flowtype)=='mixlayer') then
          call mixlayerinflow
        elseif(trim(flowtype)=='bl' .or. trim(flowtype)=='swbli') then
          !
          call profileinflow
          !
        else
          call freestreaminflow
        endif
        !
        lfirstcal=.false.
        !
      endif
      !
      i=0
      !
      if(trim(flowtype)=='bl' .or. trim(flowtype)=='swbli') then
        !
        if(turbinf=='intp')  call inflowintp
        !
      endif
      !
      do k=0,km
      do j=0,jm
        !
        ! css=sos(tmp(i,j,k))
        ! ub =vel_in(j,k,1)
        css=sos(tmp_prof(j))
        ub=vel_prof(j,1)
        ! ub =vel(i,j,k,1)
        !
        ! if(.true.) then
        if(ub>=css) then
          ! supersonic inlet
          !
          vel(i,j,k,:)=vel_in(j,k,:)
          tmp(i,j,k)  =tmp_in(j,k)
          prs(i,j,k)  =prs_in(j,k)
          rho(i,j,k)  =rho_in(j,k)
          !
          spc(i,j,k,:)=spc_in(j,k,:)
          !
        else
        ! elseif(ub<css .and. ub>=0.d0) then
          ! subsonic inlet
          ue  =extrapolate(vel(i+1,j,k,1),vel(i+2,j,k,1),dv=0.d0)
          pe  =extrapolate(prs(i+1,j,k),  prs(i+2,j,k),  dv=0.d0)
          roe =extrapolate(rho(i+1,j,k),  rho(i+2,j,k),  dv=0.d0)
          ! csse=extrapolate(sos(tmp(i+1,j,k)),sos(tmp(i+2,j,k)),dv=0.d0)
          !
          vel(i,j,k,1)=vel_in(j,k,1)
          vel(i,j,k,2)=vel_in(j,k,2)
          vel(i,j,k,3)=vel_in(j,k,3)
          tmp(i,j,k)  =tmp_in(j,k)
          !
          pwave_out  = prs(i+1,j,k)-rho(i,j,k)*css*(vel(i+1,j,k,1)-vel(i,j,k,1))
          ! ! Rudy & Strikwerda, 1980
          ! pwave_out  = pe
          ! !
          ! pwave_in   = prs_in(j,k)
          ! ! !
          ! malo=(ub/css)
          ! !
          ! prs(i,j,k)=malo*pwave_in+(1.d0-malo)*pwave_out
          prs(i,j,k) = pwave_out
          ! prs(i,j,k) = pwave_out
          !
          rho(i,j,k) =thermal(pressure=prs(i,j,k),temperature=tmp(i,j,k))
          !
          ! vel(i,j,k,1)=0.5d0*(prs_in(j,k)-pe)/(rho(i,j,k)*css)+0.5d0*(vel_in(j,k,1)+ue)
          ! prs(i,j,k)  =0.5d0*(prs_in(j,k)+pe)+0.5d0*rho(i,j,k)*css*(vel_in(j,k,1)-ue)
          ! rho(i,j,k)  =rho_in(j,k)*(prs(i,j,k)/prs_in(j,k))**(1.d0/gamma)
          ! !
          ! vel(i,j,k,2)=vel_in(j,k,2)
          ! vel(i,j,k,3)=vel_in(j,k,3)
          ! !
          ! tmp(i,j,k)  =thermal(pressure=prs(i,j,k),density=rho(i,j,k))
          !
          do jspc=1,num_species
            spc(i,j,k,jspc)=spc_in(j,k,jspc) 
          enddo
          !
          ! print*,vel_in(j,k,1),rho(i,j,k)
          !
        ! else
        !   print*,mpirank,'|',i,j,k
        !   print*,'ub=',ub,'vel_in=',vel_in(j,k,:),'css=',css
        !   stop ' !! velocity at inflow error 1 !! @ inflow'
        endif
        !
        call fvar2q(      q=  q(i,j,k,:),   density=rho(i,j,k),        &
                   velocity=vel(i,j,k,:),  pressure=prs(i,j,k),        &
                    species=spc(i,j,k,:)                               )
        !
        if(trim(turbmode)=='k-omega') then
          tke(i,j,k)=0.d0
          omg(i,j,k)=0.d0
        endif
        !
        qrhs(i,j,k,:)=0.d0
        !
      enddo
      enddo
      !
    endif
    !
    ! if(ndir==3 .and. jrk==0) then
    !   !
    !   if(lfirstcal) then
    !     !
    !     call alloinflow(ndir)
    !     !
    !     if(trim(flowtype)=='jet') then
    !       call jetinflow
    !     else
    !       call freestreaminflow
    !     endif
    !     !
    !     lfirstcal=.false.
    !     !
    !   endif
    !   !
    !   i=0
    !   do k=0,km
    !   do j=0,jm
    !     !
    !     css=sos(tmp(i,j,k))
    !     ub =vel(i,j,k,1)
    !     !
    !     if(ub>=css) then
    !       ! supersonic inlet
    !       !
    !       vel(i,j,k,:)=vel_in(j,k,:)
    !       tmp(i,j,k)  =tmp_in(j,k)
    !       prs(i,j,k)  =prs_in(j,k)
    !       rho(i,j,k)  =rho_in(j,k)
    !       !
    !       spc(i,j,k,:)=spc_in(j,k,:)
    !       !
    !     elseif(ub<css .and. ub>=0.d0) then
    !       ! subsonic inlet
    !       ue  =extrapolate(vel(i+1,j,k,1),vel(i+2,j,k,1),dv=0.d0)
    !       pe  =extrapolate(prs(i+1,j,k),  prs(i+2,j,k),  dv=0.d0)
    !       roe =extrapolate(rho(i+1,j,k),  rho(i+2,j,k),  dv=0.d0)
    !       ! csse=extrapolate(sos(tmp(i+1,j,k)),sos(tmp(i+2,j,k)),dv=0.d0)
    !       !
    !       vel(i,j,k,1)=0.5d0*(prs_in(j,k)-pe)/(rho(i,j,k)*css)+0.5d0*(vel_in(j,k,1)+ue)
    !       vel(i,j,k,2)=vel_in(j,k,2)
    !       vel(i,j,k,3)=vel_in(j,k,3)
    !       prs(i,j,k)  =0.5d0*(prs_in(j,k)+pe)+0.5d0*rho(i,j,k)*css*(vel_in(j,k,1)-ue)
    !       rho(i,j,k)  =rho_in(j,k)*(prs(i,j,k)/prs_in(j,k))**(1.d0/gamma)
    !       !
    !       tmp(i,j,k)  =thermal(pressure=prs(i,j,k),density=rho(i,j,k))
    !       !
    !       spc(i,j,k,:)=spc_in(j,k,:)
    !       !
    !       ! print*,vel_in(j,k,1),rho(i,j,k)
    !       !
    !     else
    !       stop ' !! velocity at inflow error 2 !! @ inflow'
    !     endif
    !     !
    !     call fvar2q(      q=  q(i,j,k,:),   density=rho(i,j,k),        &
    !                velocity=vel(i,j,k,:),  pressure=prs(i,j,k),        &
    !                 species=spc(i,j,k,:)                               )
    !     !
    !     qrhs(i,j,k,:)=0.d0
    !     !
    !   enddo
    !   enddo
    !   !
    ! endif
    !
    ! call mpistop
    !
  end subroutine inflow
  !+-------------------------------------------------------------------+
  !| The end of the subroutine inflow.                                 |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to interpolate the inflow turbulence.          |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 20-10-2021: Created by J. Fang @ Warrington                       |
  !+-------------------------------------------------------------------+
  subroutine inflowintp
    !
    use commvar, only : time,nstep
    use commfunc,only : cubic
    use parallel,only : mpi_imin
    use fludyna, only : thermal
    use hdf5io
    use tecio
    !
    ! local data
    real(8),pointer,dimension(:,:,:) :: vp1,vp2,vp3,vp4
    integer,save :: nstep_save=-1
    logical,save :: loadiniflow=.true.
    integer :: n,n2,n0,i,j,k
    character(len=5) :: isname
    integer :: nvp(0:3)
    real(8) :: vartime(0:3)
    real(8) :: time0,deltime,var0
    !
    if(nstep<nstep_save) then
      ! this can happen when the solver is recovered from a bakup
      loadiniflow=.true.
      nstep_save=-1
      !
      deallocate(flowvarins,timeins )
    endif
    !
    if(nstep>nstep_save) then
      !
      nstep_save=nstep
      !
      if(loadiniflow) then
        !
        allocate(flowvarins(0:jm,0:km,1:5,0:3),timeins(0:3) )
        !
        do n=0,3
          !
          if(ninflowslice<=3) then
            write(isname,'(i5.5)')n
          else 
            write(isname,'(i5.5)')n+ninflowslice-3
          endif
          !
          call h5io_init(filename='inflow/islice'//isname//'.h5',mode='read',comm=mpi_imin)
          !
          call h5read(varname='time',var=timeins(n))
          call h5read(varname='ro', var=flowvarins(0:jm,0:km,1,n),dir='i',display=.false.)
          call h5read(varname='u1', var=flowvarins(0:jm,0:km,2,n),dir='i',display=.false.)
          call h5read(varname='u2', var=flowvarins(0:jm,0:km,3,n),dir='i',display=.false.)
          call h5read(varname='u3', var=flowvarins(0:jm,0:km,4,n),dir='i',display=.false.)
          call h5read(varname='t',  var=flowvarins(0:jm,0:km,5,n),dir='i',display=.false.)
          !
          call h5io_end
          ! 
        enddo
        !
        time0=timeins(0)
        deltime=timeins(1)-timeins(0)
        !
        vp1=>flowvarins(0:jm,0:km,1:5,0)
        vp2=>flowvarins(0:jm,0:km,1:5,1)
        vp3=>flowvarins(0:jm,0:km,1:5,2)
        vp4=>flowvarins(0:jm,0:km,1:5,3)
        !
        do k=0,km
        do j=0,jm
          !
          rho_in(j,k)  =cubic(vp1(j+1,k+1,1),vp2(j+1,k+1,1),           &
                              vp3(j+1,k+1,1),vp4(j+1,k+1,1),           &
                                         time0,deltime,time) + rho_prof(j)
          vel_in(j,k,1)=cubic(vp1(j+1,k+1,2),vp2(j+1,k+1,2),           &
                              vp3(j+1,k+1,2),vp4(j+1,k+1,2),           &
                                         time0,deltime,time)
          !
          if(vel_in(j,k,1)+vel_prof(j,1)<0.d0) then
            vel_in(j,k,1)=vel_prof(j,1)
          else
            vel_in(j,k,1)=vel_prof(j,1)+vel_in(j,k,1)
          endif
          !
          vel_in(j,k,2)=cubic(vp1(j+1,k+1,3),vp2(j+1,k+1,3),           &
                              vp3(j+1,k+1,3),vp4(j+1,k+1,3),           &
                                         time0,deltime,time) + vel_prof(j,2)
          vel_in(j,k,3)=cubic(vp1(j+1,k+1,4),vp2(j+1,k+1,4),           &
                              vp3(j+1,k+1,4),vp4(j+1,k+1,4),           &
                                         time0,deltime,time)
          tmp_in(j,k)  =cubic(vp1(j+1,k+1,5),vp2(j+1,k+1,5),           &
                              vp3(j+1,k+1,5),vp4(j+1,k+1,5),           &
                                         time0,deltime,time) + tmp_prof(j)
          prs_in(j,k)  =thermal(density=rho_in(j,k),temperature=tmp_in(j,k))
          !
        end do
        end do
        !
        ! call tecbin('Testout/tecyz'//mpirankname//'.plt',                &
        !                                   x(0,0:jm,0:km,1),'x',          &
        !                                   x(0,0:jm,0:km,2),'y',          &
        !                                   x(0,0:jm,0:km,3),'z',          &
        !                                   rho_in(0:jm,0:km),'ro',        &
        !                                   vel_in(0:jm,0:km,1),'u',       &
        !                                   vel_in(0:jm,0:km,2),'v',       &
        !                                   vel_in(0:jm,0:km,3),'w',       &
        !                                   tmp_in(0:jm,0:km),'t')
        ! !
        ! if(mpirank==0) then
        !   !
        !   print*,'time0',time0
        !   print*,'deltime',deltime
        !   print*,'timeins',timeins
        !   !
        ! end if
        !
        loadiniflow=.false.
        !
      else
        !
        do n=0,3
          !
          nvp(n)=n
          vartime(n)=timeins(n)
          !
        end do
        !
        do n=0,3
          do n2=n,3
            !
            if(vartime(n)>vartime(n2)) then
              var0=vartime(n)
              vartime(n)=vartime(n2)
              vartime(n2)=var0
              !
              n0=nvp(n)
              nvp(n)=nvp(n2)
              nvp(n2)=n0
            end if
            !
          end do
        end do
        !
        vp1=>flowvarins(:,:,:,nvp(0))
        vp2=>flowvarins(:,:,:,nvp(1))
        vp3=>flowvarins(:,:,:,nvp(2))
        vp4=>flowvarins(:,:,:,nvp(3))
        !
        time0=timeins(nvp(0))
        deltime=timeins(nvp(1))-timeins(nvp(0))
        !
        ! to load next inflow slice
        if(time>timeins(nvp(2))) then
          !
          ninflowslice=ninflowslice+1
          !
          write(isname,'(i5.5)')ninflowslice
          !
          call h5io_init(filename='inflow/islice'//isname//'.h5',      &
                                              mode='read',comm=mpi_imin)
          !
          n=nvp(0)
          call h5read(varname='time',var=timeins(n))
          call h5read(varname='ro', var=flowvarins(0:jm,0:km,1,n),dir='i',display=.false.)
          call h5read(varname='u1', var=flowvarins(0:jm,0:km,2,n),dir='i',display=.false.)
          call h5read(varname='u2', var=flowvarins(0:jm,0:km,3,n),dir='i',display=.false.)
          call h5read(varname='u3', var=flowvarins(0:jm,0:km,4,n),dir='i',display=.false.)
          call h5read(varname='t',  var=flowvarins(0:jm,0:km,5,n),dir='i',display=.false.)
          !
          call h5io_end
          ! 
          vp1=>flowvarins(:,:,:,nvp(1))
          vp2=>flowvarins(:,:,:,nvp(2))
          vp3=>flowvarins(:,:,:,nvp(3))
          vp4=>flowvarins(:,:,:,nvp(0))
          !
          time0=timeins(nvp(1))
          !
          ! if(lio) then
          !   write(*,'(5(F12.6,A))')timeins(nvp(1)),' <',               &
          !                          timeins(nvp(2)),' < |',time,'| <',  &
          !                          timeins(nvp(3)),' <',timeins(nvp(0)),''
          ! endif
          !
        endif
        !
        do k=0,km
        do j=0,jm
          !
          rho_in(j,k)  =cubic(vp1(j+1,k+1,1),vp2(j+1,k+1,1),           &
                              vp3(j+1,k+1,1),vp4(j+1,k+1,1),           &
                                         time0,deltime,time) + rho_prof(j)
          vel_in(j,k,1)=cubic(vp1(j+1,k+1,2),vp2(j+1,k+1,2),           &
                              vp3(j+1,k+1,2),vp4(j+1,k+1,2),           &
                                         time0,deltime,time)
          !
          if(vel_in(j,k,1)+vel_prof(j,1)<0.d0) then
            vel_in(j,k,1)=vel_prof(j,1)
          else
            vel_in(j,k,1)=vel_prof(j,1)+vel_in(j,k,1)
          endif
          !
          vel_in(j,k,2)=cubic(vp1(j+1,k+1,3),vp2(j+1,k+1,3),           &
                              vp3(j+1,k+1,3),vp4(j+1,k+1,3),           &
                                         time0,deltime,time) + vel_prof(j,2)
          vel_in(j,k,3)=cubic(vp1(j+1,k+1,4),vp2(j+1,k+1,4),           &
                              vp3(j+1,k+1,4),vp4(j+1,k+1,4),           &
                                         time0,deltime,time)
          tmp_in(j,k)  =cubic(vp1(j+1,k+1,5),vp2(j+1,k+1,5),           &
                              vp3(j+1,k+1,5),vp4(j+1,k+1,5),           &
                                         time0,deltime,time) + tmp_prof(j)
          !
          prs_in(j,k)  =thermal(density=rho_in(j,k),temperature=tmp_in(j,k))
          !
        end do
        end do
        !
      endif
      !
    else
      !
      return
      !
    endif
    !
  end subroutine inflowintp
  !+-------------------------------------------------------------------+
  !| The end of the subroutine inflowintp.                             |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to apply inflow bc using nscbc.                |
  !+-------------------------------------------------------------------+
  !| ref: Jae Wook Kim, AIAA JOURNAL Vol. 38, No. 11, November 2000    |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 25-02-2021: Created by J. Fang @ Warrington                       |
  !+-------------------------------------------------------------------+
  subroutine inflow_nscbc(ndir)
    !
    use fludyna,   only : thermal,fvar2q,q2fvar,sos
    use commfunc,  only : deriv,ddfc
    !
    ! arguments
    integer,intent(in) :: ndir
    !
    ! local data
    integer :: i,j,k,l,jspc,ii,n,m
    real(8) :: pinv(5,5),pnor(5,5),Pmult(5,5),E(5),F(5),G(5),Rest(5),  &
               jcbi(3),LODi1(5),LODi(5)
    real(8),allocatable :: Ecs(:,:),dEcs(:),fcs(:,:),dfcs(:,:)
    real(8) :: uu,css,gmachmax2,kin,var1
    !
    logical,save :: lfirstcal=.true.
    !
    gmachmax2=0.d0
    do k=0,km
    do j=0,jm
    do i=0,im
      var1=vel(i,j,k,1)**2+vel(i,j,k,2)**2+vel(i,j,k,3)**2
      css=sos(tmp(i,j,k))
      gmachmax2=max(gmachmax2,var1/css/css)
    enddo
    enddo
    enddo
    gmachmax2=pmax(gmachmax2)
    !
    if(ndir==1 .and. irk==0) then
      !
      if(lfirstcal) then
        !
        call alloinflow(ndir)
        !
        if(trim(flowtype)=='jet') then
          call jetinflow
        else
          call freestreaminflow
        endif
        !
        lfirstcal=.false.
        !
      endif
      !
      i=0
      !
      allocate(Ecs(0:2,1:numq),dEcs(1:numq))
      ! do k=ks,ke
      ! do j=js,je
      do k=0,km
      do j=0,jm
        !
        pnor=pmatrix(rho(i,j,k),vel(i,j,k,1),vel(i,j,k,2),            &
                     vel(i,j,k,3),tmp(i,j,k),dxi(i,j,k,1,:),inv=.false.)
        pinv=pmatrix(rho(i,j,k),vel(i,j,k,1),vel(i,j,k,2),            &
                     vel(i,j,k,3),tmp(i,j,k),dxi(i,j,k,1,:),inv=.true.)
        !
        ! Pmult=MatMul(Pinv,pnor)
        ! if(jrk==0) then
        !   print*,'---------------------------------------------------------'
        !   write(*,"(5(F7.4))")Pmult(1,:)
        !   write(*,"(5(F7.4))")Pmult(2,:)
        !   write(*,"(5(F7.4))")Pmult(3,:)
        !   write(*,"(5(F7.4))")Pmult(4,:)
        !   write(*,"(5(F7.4))")Pmult(5,:)
        ! end if
        !
        do ii=0,2
          uu=dxi(i+ii,j,k,1,1)*vel(i+ii,j,k,1) +                       &
             dxi(i+ii,j,k,1,2)*vel(i+ii,j,k,2) +                       &
             dxi(i+ii,j,k,1,3)*vel(i+ii,j,k,3)
          !
          Ecs(ii,1)=jacob(i+ii,j,k)*  q(i+ii,j,k,1)*uu
          Ecs(ii,2)=jacob(i+ii,j,k)*( q(i+ii,j,k,2)*uu+dxi(i+ii,j,k,1,1)*prs(i+ii,j,k) )
          Ecs(ii,3)=jacob(i+ii,j,k)*( q(i+ii,j,k,3)*uu+dxi(i+ii,j,k,1,2)*prs(i+ii,j,k) )
          Ecs(ii,4)=jacob(i+ii,j,k)*( q(i+ii,j,k,4)*uu+dxi(i+ii,j,k,1,3)*prs(i+ii,j,k) )
          Ecs(ii,5)=jacob(i+ii,j,k)*( q(i+ii,j,k,5)+prs(i+ii,j,k) )*uu
          do jspc=1,num_species
            Ecs(ii,5+jspc)=jacob(i+ii,j,k)*q(i+ii,j,k,5+jspc)*uu
          enddo
        enddo
        !
        do n=1,numq
          dEcs(n)= deriv( Ecs(0,n),Ecs(1,n) ) !,Ecs(2,n) )
        enddo
        !
        E(1)= q(i,j,k,2)
        E(2)= q(i,j,k,2)*vel(i,j,k,1)+prs(i,j,k)
        E(3)= q(i,j,k,3)*vel(i,j,k,1)
        E(4)= q(i,j,k,4)*vel(i,j,k,1)
        E(5)=(q(i,j,k,5)+prs(i,j,k))*vel(i,j,k,1)
        !
        F(1)= q(i,j,k,3)
        F(2)= q(i,j,k,2)*vel(i,j,k,2)
        F(3)= q(i,j,k,3)*vel(i,j,k,2)+prs(i,j,k)
        F(4)= q(i,j,k,4)*vel(i,j,k,2)
        F(5)=(q(i,j,k,5)+prs(i,j,k))*vel(i,j,k,2)
        !
        G(1)= q(i,j,k,4)
        G(2)= q(i,j,k,2)*vel(i,j,k,3)
        G(3)= q(i,j,k,3)*vel(i,j,k,3)
        G(4)= q(i,j,k,4)*vel(i,j,k,3)+prs(i,j,k)
        G(5)=(q(i,j,k,5)+prs(i,j,k))*vel(i,j,k,3)
        !
        jcbi(1)= deriv( dxi(i,j,k,1,1)  *jacob(i,j,k),                 &
                        dxi(i+1,j,k,1,1)*jacob(i+1,j,k) )
        jcbi(2)= deriv( dxi(i,j,k,1,2)  *jacob(i,j,k),                 &
                        dxi(i+1,j,k,1,2)*jacob(i+1,j,k) )
        jcbi(3)= deriv( dxi(i,j,k,1,3)  *jacob(i,j,k),                 &
                        dxi(i+1,j,k,1,3)*jacob(i+1,j,k) )
        !
        Rest(1)=E(1)*Jcbi(1)+F(1)*Jcbi(2)+G(1)*Jcbi(3)
        Rest(2)=E(2)*Jcbi(1)+F(2)*Jcbi(2)+G(2)*Jcbi(3)
        Rest(3)=E(3)*Jcbi(1)+F(3)*Jcbi(2)+G(3)*Jcbi(3)
        Rest(4)=E(4)*Jcbi(1)+F(4)*Jcbi(2)+G(4)*Jcbi(3)
        Rest(5)=E(5)*Jcbi(1)+F(5)*Jcbi(2)+G(5)*Jcbi(3)
        !
        LODi1(1)=dEcs(1)-Rest(1)
        LODi1(2)=dEcs(2)-Rest(2)
        LODi1(3)=dEcs(3)-Rest(3)
        LODi1(4)=dEcs(4)-Rest(4)
        LODi1(5)=dEcs(5)-Rest(5)
        !
        LODi=MatMul(pinv,LODi1)/jacob(i,j,k)
        !
        css=sos(tmp(i,j,k))
        !
        kin=0.25d0*(1.d0-gmachmax2)*css/(xmax-xmin)
        !
        var1=1.d0/sqrt( dxi(i,j,k,1,1)**2+dxi(i,j,k,1,2)**2+           &
                        dxi(i,j,k,1,3)**2 )
        !
        LODi(1)=0.d0
        LODi(2)=kin*0.5d0*                                             &
                   ( dxi(i,j,k,1,1)*var1*(vel(i,j,k,3)-vel_in(j,k,3))- &
                     dxi(i,j,k,1,3)*var1*(vel(i,j,k,1)-vel_in(j,k,1)) )
        LODi(3)=kin*0.5d0*                                             &
                   (-dxi(i,j,k,1,1)*var1*(vel(i,j,k,2)-vel_in(j,k,2))+ &
                     dxi(i,j,k,1,2)*var1*(vel(i,j,k,1)-vel_in(j,k,1)) )
        LODi(4)=kin*(dxi(i,j,k,1,1)*var1*(vel(i,j,k,1)-vel_in(j,k,1))+ &
                     dxi(i,j,k,1,2)*var1*(vel(i,j,k,2)-vel_in(j,k,2))+ &
                     dxi(i,j,k,1,3)*var1*(vel(i,j,k,3)-vel_in(j,k,3))+ &
                     (prs(i,j,k)-prs_in(j,k))/rho(i,j,k)/css )
        !
        LODi1=MatMul(pnor,LODi)*jacob(i,j,k)
        !
        dEcs(1)=LODi1(1)+Rest(1)
        dEcs(2)=LODi1(2)+Rest(2)
        dEcs(3)=LODi1(3)+Rest(3)
        dEcs(4)=LODi1(4)+Rest(4)
        dEcs(5)=LODi1(5)+Rest(5)
        !
        uu=dxi(i,j,k,1,1)*vel(i,j,k,1) +                               &
           dxi(i,j,k,1,2)*vel(i,j,k,2) +                               &
           dxi(i,j,k,1,3)*vel(i,j,k,3)
        do jspc=1,num_species
          var1=kin*0.5d0*(spc(i,j,k,jspc)-spc_in(j,k,jspc))
          !
          dEcs(5+jspc)=jacob(i,j,k)*uu*var1
        enddo
        !
        qrhs(i,j,k,:)=qrhs(i,j,k,:)+dEcs(:)
        !
      enddo
      enddo
      !
      deallocate(Ecs,dEcs)
      !
      allocate(fcs(-hm:jm+hm,1:numq),dfcs(0:jm,1:numq))
      do k=ks,ke
        !
        do j=-hm,jm+hm
          !
          uu=dxi(i,j,k,2,1)*vel(i,j,k,1)+dxi(i,j,k,2,2)*vel(i,j,k,2) + &
             dxi(i,j,k,2,3)*vel(i,j,k,3)
          fcs(j,1)=jacob(i,j,k)*  q(i,j,k,1)*uu
          fcs(j,2)=jacob(i,j,k)*( q(i,j,k,2)*uu+dxi(i,j,k,2,1)*prs(i,j,k) )
          fcs(j,3)=jacob(i,j,k)*( q(i,j,k,3)*uu+dxi(i,j,k,2,2)*prs(i,j,k) )
          fcs(j,4)=jacob(i,j,k)*( q(i,j,k,4)*uu+dxi(i,j,k,2,3)*prs(i,j,k) )
          fcs(j,5)=jacob(i,j,k)*( q(i,j,k,5)+prs(i,j,k) )*uu
          do jspc=1,num_species
            fcs(j,5+jspc)=jacob(i,j,k)*q(i,j,k,5+jspc)*uu
          enddo
          !
        enddo
        !
        do n=1,numq
          dfcs(:,n)=ddfc(fcs(:,n),'222e',npdcj,jm)
        enddo
        !
        qrhs(i,js:je,k,:)=qrhs(i,js:je,k,:)+dfcs(js:je,:)
        !
      enddo
      !
      deallocate(fcs,dfcs)
      ! !
      if(ndims==3) then
        allocate(fcs(-hm:km+hm,1:numq),dfcs(0:km,1:numq))
        do j=js,je
          !
          do k=-hm,km+hm
            !
            uu=dxi(i,j,k,3,1)*vel(i,j,k,1)+dxi(i,j,k,3,2)*vel(i,j,k,2) + &
               dxi(i,j,k,3,3)*vel(i,j,k,3)
            fcs(k,1)=jacob(i,j,k)*  q(i,j,k,1)*uu
            fcs(k,2)=jacob(i,j,k)*( q(i,j,k,2)*uu+dxi(i,j,k,3,1)*prs(i,j,k) )
            fcs(k,3)=jacob(i,j,k)*( q(i,j,k,3)*uu+dxi(i,j,k,3,2)*prs(i,j,k) )
            fcs(k,4)=jacob(i,j,k)*( q(i,j,k,4)*uu+dxi(i,j,k,3,3)*prs(i,j,k) )
            fcs(k,5)=jacob(i,j,k)*( q(i,j,k,5)+prs(i,j,k) )*uu
            do jspc=1,num_species
              fcs(k,5+jspc)=jacob(i,j,k)*q(i,j,k,5+jspc)*uu
            enddo
            !
          enddo
          !
          do n=1,numq
            dfcs(:,n)=ddfc(fcs(:,n),'222e',npdck,km)
          enddo
          !
          qrhs(i,j,ks:ke,:)=qrhs(i,j,ks:ke,:)+dfcs(ks:ke,:)
          !
        enddo
        !
        deallocate(fcs,dfcs)
      endif
      !
    endif
    !
  end subroutine inflow_nscbc
  !+-------------------------------------------------------------------+
  !| The end of the subroutine inflow_nscbc.                           |
  !+-------------------------------------------------------------------+
  !!
  !+-------------------------------------------------------------------+
  !| This subroutine is to apply far-field bc.                         |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 04-04-2021: Created by J. Fang @ Warrington                       |
  !+-------------------------------------------------------------------+
  subroutine farfield(ndir)
    !
    use fludyna,   only : thermal,fvar2q,q2fvar,sos,postshock
    use commfunc,  only : extrapolate
    !
    ! arguments
    integer,intent(in) :: ndir
    !
    ! local data
    integer :: i,j,k,l,jspec
    real(8) :: css,csse,ub,pe,roe,ue,ve,we,spce(1:num_species),        &
               vnb,vtb,vne,vte
    real(8) :: var1
    !
    logical,save :: lfirstcal=.true.
    !
    if(ndir==3 .and. jrk==0) then
      !
      j=0
      !
      ! if(lfirstcal) then
      !   !
      !   allocate(bvec_jm(0:im,0:km,1:3))
      !   !
      !   do k=0,km
      !   do i=0,im
      !     var1=sqrt( dxi(i,j,k,2,1)**2+dxi(i,j,k,2,2)**2+              &
      !                dxi(i,j,k,2,3)**2 )
      !     !
      !     bvec_jm(i,k,1)=dxi(i,j,k,2,1)/var1
      !     bvec_jm(i,k,2)=dxi(i,j,k,2,2)/var1
      !     bvec_jm(i,k,3)=dxi(i,j,k,2,3)/var1
      !     !
      !     ! print*,bvec_im(:,j,k)
      !   enddo
      !   enddo
      !   !
      !   lfirstcal=.false.
      !   !
      ! endif
      !
      do k=0,km
      do i=0,im
        !
        css=sos(tmp(i,j,k))
        ! ub =vel(i,j,k,1)*bvec_jm(i,k,1)+vel(i,j,k,2)*bvec_jm(i,k,2)+   &
        !     vel(i,j,k,3)*bvec_jm(i,k,3)
        !
        ue  =extrapolate(vel(i,j+1,k,1),vel(i,j+2,k,1),dv=0.d0)
        ve  =extrapolate(vel(i,j+1,k,2),vel(i,j+2,k,2),dv=0.d0)
        we  =extrapolate(vel(i,j+1,k,3),vel(i,j+2,k,3),dv=0.d0)
        pe  =extrapolate(prs(i,j+1,k),  prs(i,j+2,k),dv=0.d0)
        roe =extrapolate(rho(i,j+1,k),  rho(i,j+2,k),dv=0.d0)
        csse=extrapolate(sos(tmp(i,j+1,k)),sos(tmp(i,j+2,k)),dv=0.d0)
        !
        do jspec=1,num_species
          spce(jspec)=extrapolate(spc(i,j+1,k,jspec),                  &
                                  spc(i,j+2,k,jspec),dv=0.d0)
        enddo
        !
        if(vel(i,j,k,2)>=0.d0) then
          ! subsonic inflow
          vel(i,j,k,1)=uinf
          vel(i,j,k,2)=0.5d0*(pinf-pe)/(rho(i,j,k)*css)+0.5d0*(vinf+ve)
          vel(i,j,k,3)=winf
          prs(i,j,k)  =0.5d0*(pinf+pe)+0.5d0*rho(i,j,k)*css*(vinf-ve)
          rho(i,j,k)  =roinf*(prs(i,j,k)/pinf)**(1.d0/gamma)
          !
          tmp(i,j,k)  =thermal(pressure=prs(i,j,k),density=rho(i,j,k))
          !
          spc(i,j,k,:)=0.d0
          !
        else
          ! subsonic outflow
          !
          prs(i,j,k)=pinf
          rho(i,j,k)=roe+(prs(i,j,k)-pe)/csse/csse
          !
          vel(i,j,k,1)=ue
          vel(i,j,k,2)=ve-(pe-prs(i,j,k))/roe/Csse
          vel(i,j,k,3)=we
          tmp(i,j,k)  =thermal(pressure=prs(i,j,k),density=rho(i,j,k))
          spc(i,j,k,:)=spce(:)
          !
        endif
        !
        call fvar2q(      q=  q(i,j,k,:),   density=rho(i,j,k),        &
                   velocity=vel(i,j,k,:),  pressure=prs(i,j,k),        &
                    species=spc(i,j,k,:)                               )
        !
        qrhs(i,j,k,:)=0.d0
        !
      enddo
      enddo
      !
    endif
    !
    if(ndir==4 .and. jrk==jrkm) then
      !
      j=jm
      !
      if(lfirstcal) then
        !
        allocate(rho_far(0:im),vel_far(0:im,1:3),tmp_far(0:im),        & 
                 prs_far(0:im))
        !
        do i=0,im
          !
          if(trim(flowtype)=='swbli') then
            !
            if(x(i,j,0,1)<=xrhjump) then
              rho_far(i)=rho_prof(jm)
              vel_far(i,:)=vel_prof(jm,:)
              tmp_far(i)=tmp_prof(jm)
              prs_far(i)=prs_prof(jm)
            else
              call postshock(rho_prof(jm),vel_prof(jm,1),vel_prof(jm,2), &
                             prs_prof(jm),tmp_prof(jm),                  &
                             rho_far(i),vel_far(i,1),vel_far(i,2),       &
                             prs_far(i),tmp_far(i),angshk)
            endif
            !
          elseif(trim(flowtype)=='bl') then
            !
            rho_far(i)=rho_prof(jm)
            vel_far(i,:)=vel_prof(jm,:)
            tmp_far(i)=tmp_prof(jm)
            prs_far(i)=prs_prof(jm)
            !
          else
            !
            rho_far(i)  =roinf
            vel_far(i,1)=uinf
            vel_far(i,2)=vinf
            vel_far(i,3)=winf
            tmp_far(i)  =tinf
            prs_far(i)  =pinf
            !
          endif
          !
        enddo
        !
        lfirstcal=.false.
        !
      endif
      !
      do k=0,km
      do i=0,im
        !
        css=sos(tmp(i,j,k))
        ! ub =vel(i,j,k,1)*bvec_jm(i,k,1)+vel(i,j,k,2)*bvec_jm(i,k,2)+   &
        !     vel(i,j,k,3)*bvec_jm(i,k,3)
        !
        ue  =extrapolate(vel(i,j-1,k,1),vel(i,j-2,k,1),dv=0.d0)
        ve  =extrapolate(vel(i,j-1,k,2),vel(i,j-2,k,2),dv=0.d0)
        we  =extrapolate(vel(i,j-1,k,3),vel(i,j-2,k,3),dv=0.d0)
        pe  =extrapolate(prs(i,j-1,k),  prs(i,j-2,k),dv=0.d0)
        roe =extrapolate(rho(i,j-1,k),  rho(i,j-2,k),dv=0.d0)
        csse=extrapolate(sos(tmp(i,j-1,k)),sos(tmp(i,j-2,k)),dv=0.d0)
        !
        do jspec=1,num_species
          spce(jspec)=extrapolate(spc(i,j-1,k,jspec),                  &
                                  spc(i,j-2,k,jspec),dv=0.d0)
        enddo
        !
        ! if(vel(i,j,k,2)<=0.d0) then
        !   ! subsonic inflow
        !   ! vel(i,j,k,1)=uinf
        !   ! vel(i,j,k,2)=-0.5d0*(pinf-pe)/(rho(i,j,k)*css)+0.5d0*(vinf+ve)
        !   ! vel(i,j,k,3)=winf
        !   ! prs(i,j,k)  =0.5d0*(pinf+pe)-0.5d0*rho(i,j,k)*css*(vinf-ve)
        !   ! rho(i,j,k)  =roinf*(prs(i,j,k)/pinf)**(1.d0/gamma)
        !   ! !
        !   ! tmp(i,j,k)  =thermal(pressure=prs(i,j,k),density=rho(i,j,k))
        !   ! !
        !   ! spc(i,j,k,:)=0.d0
        !   !
        !   vel(i,j,k,1)=vel_far(i,1)
        !   vel(i,j,k,2)=-0.5d0*(prs_far(i)-pe)/(rho(i,j,k)*css) +       &
        !                                          0.5d0*(vel_far(i,2)+ve)
        !   vel(i,j,k,3)=0.d0
        !   prs(i,j,k)  =0.5d0*(prs_far(i)+pe) -                         &
        !                           0.5d0*rho(i,j,k)*css*(vel_far(i,2)-ve)
        !   rho(i,j,k)  =rho_far(i)*(prs(i,j,k)/prs_far(i))**(1.d0/gamma)
        !   !
        !   tmp(i,j,k)  =thermal(pressure=prs(i,j,k),density=rho(i,j,k))
        !   !
        !   spc(i,j,k,:)=0.d0
        !   !
        ! else
        !   ! subsonic outflow
        !   !
        !   ! prs(i,j,k)=pinf
        !   prs(i,j,k)=prs_far(i)
        !   rho(i,j,k)=roe+(prs(i,j,k)-pe)/csse/csse
        !   !
        !   vel(i,j,k,1)=ue
        !   vel(i,j,k,2)=ve+(pe-prs(i,j,k))/roe/Csse
        !   vel(i,j,k,3)=we
        !   tmp(i,j,k)  =thermal(pressure=prs(i,j,k),density=rho(i,j,k))
        !   spc(i,j,k,:)=spce(:)
        !   !
        ! endif
          prs(i,j,k)=pe
          rho(i,j,k)=roe
          !
          vel(i,j,k,1)=ue
          vel(i,j,k,2)=ve
          vel(i,j,k,3)=we
          tmp(i,j,k)  =thermal(pressure=prs(i,j,k),density=rho(i,j,k))
          spc(i,j,k,:)=spce(:)
        !
        call fvar2q(      q=  q(i,j,k,:),   density=rho(i,j,k),        &
                   velocity=vel(i,j,k,:),  pressure=prs(i,j,k),        &
                    species=spc(i,j,k,:)                               )
        !
        qrhs(i,j,k,:)=0.d0
        !
      enddo
      enddo
      !
    endif
    !
    if(ndir==5 .and. krk==0) then
      !
      k=0
      !
      ! if(lfirstcal) then
      !   !
      !   allocate(bvec_jm(0:im,0:km,1:3))
      !   !
      !   do k=0,km
      !   do i=0,im
      !     var1=sqrt( dxi(i,j,k,2,1)**2+dxi(i,j,k,2,2)**2+              &
      !                dxi(i,j,k,2,3)**2 )
      !     !
      !     bvec_jm(i,k,1)=dxi(i,j,k,2,1)/var1
      !     bvec_jm(i,k,2)=dxi(i,j,k,2,2)/var1
      !     bvec_jm(i,k,3)=dxi(i,j,k,2,3)/var1
      !     !
      !     ! print*,bvec_im(:,j,k)
      !   enddo
      !   enddo
      !   !
      !   lfirstcal=.false.
      !   !
      ! endif
      !
      do j=0,jm
      do i=0,im
        !
        css=sos(tmp(i,j,k))
        ! ub =vel(i,j,k,1)*bvec_jm(i,k,1)+vel(i,j,k,2)*bvec_jm(i,k,2)+   &
        !     vel(i,j,k,3)*bvec_jm(i,k,3)
        !
        ue  =extrapolate(vel(i,j,k+1,1),vel(i,j,k+2,1),dv=0.d0)
        ve  =extrapolate(vel(i,j,k+1,2),vel(i,j,k+2,2),dv=0.d0)
        we  =extrapolate(vel(i,j,k+1,3),vel(i,j,k+2,3),dv=0.d0)
        pe  =extrapolate(prs(i,j,k+1),  prs(i,j,k+2),dv=0.d0)
        roe =extrapolate(rho(i,j,k+1),  rho(i,j,k+2),dv=0.d0)
        csse=extrapolate(sos(tmp(i,j,k+1)),sos(tmp(i,j,k+2)),dv=0.d0)
        !
        do jspec=1,num_species
          spce(jspec)=extrapolate(spc(i,j,k+1,jspec),                  &
                                  spc(i,j,k+2,jspec),dv=0.d0)
        enddo
        !
        if(vel(i,j,k,3)>=0.d0) then
          ! subsonic inflow
          vel(i,j,k,1)=uinf
          vel(i,j,k,2)=vinf
          vel(i,j,k,3)=0.5d0*(pinf-pe)/(rho(i,j,k)*css)+0.5d0*(winf+we)
          prs(i,j,k)  =0.5d0*(pinf+pe)+0.5d0*rho(i,j,k)*css*(winf-we)
          rho(i,j,k)  =roinf*(prs(i,j,k)/pinf)**(1.d0/gamma)
          !
          tmp(i,j,k)  =thermal(pressure=prs(i,j,k),density=rho(i,j,k))
          !
          spc(i,j,k,:)=0.d0
          !
        else
          ! subsonic outflow
          !
          prs(i,j,k)=pinf
          rho(i,j,k)=roe+(prs(i,j,k)-pe)/csse/csse
          !
          vel(i,j,k,1)=ue
          vel(i,j,k,2)=ve
          vel(i,j,k,3)=we-(pe-prs(i,j,k))/roe/Csse
          tmp(i,j,k)  =thermal(pressure=prs(i,j,k),density=rho(i,j,k))
          spc(i,j,k,:)=spce(:)
          !
        endif
        !
        call fvar2q(      q=  q(i,j,k,:),   density=rho(i,j,k),        &
                   velocity=vel(i,j,k,:),  pressure=prs(i,j,k),        &
                    species=spc(i,j,k,:)                               )
        !
      enddo
      enddo
      !
    endif
    !
    if(ndir==6 .and. krk==krkm) then
      !
      k=km
      !
      ! if(lfirstcal) then
      !   !
      !   allocate(bvec_jm(0:im,0:km,1:3))
      !   !
      !   do k=0,km
      !   do i=0,im
      !     var1=sqrt( dxi(i,j,k,2,1)**2+dxi(i,j,k,2,2)**2+              &
      !                dxi(i,j,k,2,3)**2 )
      !     !
      !     bvec_jm(i,k,1)=dxi(i,j,k,2,1)/var1
      !     bvec_jm(i,k,2)=dxi(i,j,k,2,2)/var1
      !     bvec_jm(i,k,3)=dxi(i,j,k,2,3)/var1
      !     !
      !     ! print*,bvec_im(:,j,k)
      !   enddo
      !   enddo
      !   !
      !   lfirstcal=.false.
      !   !
      ! endif
      !
      do j=0,jm
      do i=0,im
        !
        css=sos(tmp(i,j,k))
        ! ub =vel(i,j,k,1)*bvec_jm(i,k,1)+vel(i,j,k,2)*bvec_jm(i,k,2)+   &
        !     vel(i,j,k,3)*bvec_jm(i,k,3)
        !
        ue  =extrapolate(vel(i,j,k-1,1),vel(i,j,k-2,1),dv=0.d0)
        ve  =extrapolate(vel(i,j,k-1,2),vel(i,j,k-2,2),dv=0.d0)
        we  =extrapolate(vel(i,j,k-1,3),vel(i,j,k-2,3),dv=0.d0)
        pe  =extrapolate(prs(i,j,k-1),  prs(i,j,k-2),dv=0.d0)
        roe =extrapolate(rho(i,j,k-1),  rho(i,j,k-2),dv=0.d0)
        csse=extrapolate(sos(tmp(i,j,k-1)),sos(tmp(i,j,k-2)),dv=0.d0)
        !
        do jspec=1,num_species
          spce(jspec)=extrapolate(spc(i,j,k-1,jspec),                  &
                                  spc(i,j,k-2,jspec),dv=0.d0)
        enddo
        !
        if(vel(i,j,k,3)<=0.d0) then
          ! subsonic inflow
          vel(i,j,k,1)=uinf
          vel(i,j,k,2)=vinf
          vel(i,j,k,3)=-0.5d0*(pinf-pe)/(rho(i,j,k)*css)+0.5d0*(winf+we)
          prs(i,j,k)  = 0.5d0*(pinf+pe)-0.5d0*rho(i,j,k)*css*(winf-we)
          rho(i,j,k)  =roinf*(prs(i,j,k)/pinf)**(1.d0/gamma)
          !
          tmp(i,j,k)  =thermal(pressure=prs(i,j,k),density=rho(i,j,k))
          !
          spc(i,j,k,:)=0.d0
          !
        else
          ! subsonic outflow
          prs(i,j,k)=pinf
          rho(i,j,k)=roe+(prs(i,j,k)-pe)/csse/csse
          !
          vel(i,j,k,1)=ue
          vel(i,j,k,2)=ve
          vel(i,j,k,3)=we+(pe-prs(i,j,k))/roe/Csse
          tmp(i,j,k)  =thermal(pressure=prs(i,j,k),density=rho(i,j,k))
          spc(i,j,k,:)=spce(:)
          !
        endif
        !
        call fvar2q(      q=  q(i,j,k,:),   density=rho(i,j,k),        &
                   velocity=vel(i,j,k,:),  pressure=prs(i,j,k),        &
                    species=spc(i,j,k,:)                               )
        !
      enddo
      enddo
      !
    endif
    !
  end subroutine farfield
  !+-------------------------------------------------------------------+
  !| The end of the subroutine farfield.                               |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to apply outflow bc.                           |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 25-02-2021: Created by J. Fang @ Warrington                       |
  !+-------------------------------------------------------------------+
  subroutine outflow(ndir)
    !
    use commvar,   only : turbmode,deltat
    use fludyna,   only : thermal,fvar2q,q2fvar,sos
    use commfunc,  only : extrapolate
    use commarray, only : tke,omg
    !
    ! arguments
    integer,intent(in) :: ndir
    !
    ! local data
    integer :: i,j,k,l,jspec
    real(8) :: css,csse,ub,pe,roe,ue,ve,we,spce(1:num_species),        &
               vnb,vtb,vne,vte,alpha,pwave_in,pwave_out
    real(8) :: var1
    !
    logical,save :: lfirstcal=.true.
    !
    alpha=0.3d0
    !
    if(ndir==2 .and. irk==irkm) then
      !
      i=im
      !
      if(lfirstcal) then
        !
        allocate(bvec_im(0:jm,0:km,1:3))
        !
        do k=0,km
        do j=0,jm
          var1=sqrt( dxi(i,j,k,1,1)**2+dxi(i,j,k,1,2)**2+              &
                     dxi(i,j,k,1,3)**2 )
          !
          bvec_im(j,k,1)=dxi(i,j,k,1,1)/var1
          bvec_im(j,k,2)=dxi(i,j,k,1,2)/var1
          bvec_im(j,k,3)=dxi(i,j,k,1,3)/var1
          !
          ! print*,bvec_im(:,j,k)
        enddo
        enddo
        !
        lfirstcal=.false.
        !
      endif
      !
      do k=0,km
      do j=0,jm
        !
        css=sos(tmp(i,j,k))
        ub =vel(i,j,k,1)*bvec_im(j,k,1)+vel(i,j,k,2)*bvec_im(j,k,2)+   &
            vel(i,j,k,3)*bvec_im(j,k,3)
        !
        ue  =extrapolate(vel(i-1,j,k,1),vel(i-2,j,k,1),dv=0.d0)
        ve  =extrapolate(vel(i-1,j,k,2),vel(i-2,j,k,2),dv=0.d0)
        we  =extrapolate(vel(i-1,j,k,3),vel(i-2,j,k,3),dv=0.d0)
        pe  =extrapolate(prs(i-1,j,k),  prs(i-2,j,k),dv=0.d0)
        roe =extrapolate(rho(i-1,j,k),  rho(i-2,j,k),dv=0.d0)
        csse=extrapolate(sos(tmp(i-1,j,k)),sos(tmp(i-2,j,k)),dv=0.d0)
        !
        do jspec=1,num_species
          spce(jspec)=extrapolate(spc(i-1,j,k,jspec),                  &
                                  spc(i-2,j,k,jspec),dv=0.d0)
        enddo
        !
        vne=ue*bvec_im(j,k,1)+ve*bvec_im(j,k,2)
        vte=ue*bvec_im(j,k,2)-ve*bvec_im(j,k,1)
        if(.true.) then
        ! if(ub>=css) then
          ! supersonic outlet
          !
          vel(i,j,k,1)=ue 
          vel(i,j,k,2)=ve 
          vel(i,j,k,3)=we
          prs(i,j,k)  =pe
          rho(i,j,k)  =roe
          !
        else !if(ub<css .and. ub>=0.d0) then
          ! subsonic outlet
          !
          ! prs(i,j,k)= pinf
          ! rho(i,j,k)= roe+(prs(i,j,k)-pe)/csse/csse
          ! vel(i,j,k,1)= ue+ (pe-prs(i,j,k))/roe/csse
          ! vel(i,j,k,2)= ve
          ! vel(i,j,k,3)= we
          !
          ! pwave_in = (prs(i,j,k)+alpha*deltat*pinf+rho(i,j,k)*css*(ue-vel(i,j,k,1)))/(1.d0+alpha*deltat)
          ! RUDY & STRIKWERDA, 1981
          pwave_in = pinf
          pwave_out=pe
          !
          var1=ub/css
          prs(i,j,k)=var1*pwave_out+(1.d0-var1)*pwave_in
          !
          vel(i,j,k,1)=ue 
          vel(i,j,k,2)=ve 
          vel(i,j,k,3)=we
          !
          rho(i,j,k)  =roe
          !
        ! else
        !   stop ' !! velocity at outflow error !! @ outflow'
        endif
        !
        tmp(i,j,k)  =thermal(pressure=prs(i,j,k),density=rho(i,j,k))
        spc(i,j,k,:)=spce(:)
        !
        call fvar2q(      q=  q(i,j,k,:),   density=rho(i,j,k),        &
                   velocity=vel(i,j,k,:),  pressure=prs(i,j,k),        &
                    species=spc(i,j,k,:)                               )
        !
        qrhs(i,j,k,:)=0.d0
        !
        if(trim(turbmode)=='k-omega') then
          tke(i,j,k)=tke(i-1,j,k)
          omg(i,j,k)=omg(i-1,j,k)
        endif
        !
      enddo
      enddo
      !
    endif
    !
  end subroutine outflow
  !+-------------------------------------------------------------------+
  !| The end of the subroutine outflow.                                |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to apply farfield bc using nscbc.              |
  !+-------------------------------------------------------------------+
  !| ref: Jae Wook Kim, AIAA JOURNAL Vol. 38, No. 11, November 2000    |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 05-03-2021: Created by J. Fang @ Warrington                       |
  !+-------------------------------------------------------------------+
  subroutine farfield_nscbc(ndir)
    !
    use commvar,   only : xmin,xmax,ymin,ymax,mach
    use fludyna,   only : thermal,fvar2q,q2fvar,sos
    use commfunc,  only : deriv,ddfc,spafilter6exp
    !
    ! arguments
    integer,intent(in) :: ndir
    !
    ! local data
    integer :: i,j,k,l,jspc,jmod,ii,n,m,jq
    real(8) :: pinv(5,5),pnor(5,5),Pmult(5,5),E(5),F(5),G(5),Rest(5),  &
               jcbi(3),LODi1(5),LODi(5)
    real(8),allocatable :: Ecs(:,:),dEcs(:),fcs(:,:),dfcs(:,:),qfilt(:,:)
    real(8) :: uu,css,gmachmax2,kinout,kin,var1,var2
    !
    logical,save :: lfirstcal=.true.
    !
    if(lfirstcal) then
      !
      uinf_j0=uinf
      uinf_jm=uinf
      !
      lfirstcal=.false.
      !
    endif
    !
    gmachmax2=0.d0
    if(ndir==3 .and. jrk==0) then
      j=0
      do k=0,km
      do i=0,im
        var1=1.d0/( dxi(i,j,k,2,1)**2+dxi(i,j,k,2,2)**2+               &
                    dxi(i,j,k,2,3)**2 )
        var2=vel(i,j,k,1)*dxi(i,j,k,2,1)+vel(i,j,k,2)*dxi(i,j,k,2,2)+  &
             vel(i,j,k,3)*dxi(i,j,k,2,3)
        css=sos(tmp(i,j,k))
        gmachmax2=max(gmachmax2,var2*var2*var1/css/css)
      enddo
      enddo
    endif
    if(ndir==4 .and. jrk==jrkm) then
      j=jm
      do k=0,km
      do i=0,im
        var1=1.d0/( dxi(i,j,k,2,1)**2+dxi(i,j,k,2,2)**2+               &
                    dxi(i,j,k,2,3)**2 )
        var2=vel(i,j,k,1)*dxi(i,j,k,2,1)+vel(i,j,k,2)*dxi(i,j,k,2,2)+  &
             vel(i,j,k,3)*dxi(i,j,k,2,3)
        css=sos(tmp(i,j,k))
        gmachmax2=max(gmachmax2,var2*var2*var1/css/css)
      enddo
      enddo
    endif
    if(ndir==5 .and. krk==0) then
      k=0
      do j=0,jm
      do i=0,im
        var1=1.d0/( dxi(i,j,k,3,1)**2+dxi(i,j,k,3,2)**2+               &
                    dxi(i,j,k,3,3)**2 )
        var2=vel(i,j,k,1)*dxi(i,j,k,3,1)+vel(i,j,k,2)*dxi(i,j,k,3,2)+  &
             vel(i,j,k,3)*dxi(i,j,k,3,3)
        css=sos(tmp(i,j,k))
        gmachmax2=max(gmachmax2,var2*var2*var1/css/css)
      enddo
      enddo
    endif
    if(ndir==6 .and. krk==krkm) then
      k=km
      do j=0,jm
      do i=0,im
        var1=1.d0/( dxi(i,j,k,3,1)**2+dxi(i,j,k,3,2)**2+               &
                    dxi(i,j,k,3,3)**2 )
        var2=vel(i,j,k,1)*dxi(i,j,k,3,1)+vel(i,j,k,2)*dxi(i,j,k,3,2)+  &
             vel(i,j,k,3)*dxi(i,j,k,3,3)
        css=sos(tmp(i,j,k))
        gmachmax2=max(gmachmax2,var2*var2*var1/css/css)
      enddo
      enddo
    endif
    gmachmax2=pmax(gmachmax2)
    !
    if(ndir==3 .and. jrk==0) then
      !
      j=0
      !
      allocate(qfilt(0:im,1:numq))
      do k=0,km
        !
        do jq=1,numq
          qfilt(:,jq)=spafilter6exp(q(:,j,k,jq),npdci,im)
          q(0:im,j,k,jq)=qfilt(:,jq)
        enddo
        call q2fvar(      q=  q(0:im,j,k,:),                   &
                    density=rho(0:im,j,k),                     &
                   velocity=vel(0:im,j,k,:),                   &
                   pressure=prs(0:im,j,k),                     &
                temperature=tmp(0:im,j,k),                     &
                    species=spc(0:im,j,k,:)                    )
        !
      enddo
      deallocate(qfilt)
      !
      if(ndims==3) then 
        !
        allocate(qfilt(0:km,1:numq))
        do i=0,im
          !
          do jq=1,numq
            qfilt(:,jq)=spafilter6exp(q(i,j,:,jq),npdck,km)
            q(i,j,0:km,jq)=qfilt(:,jq)
          enddo
          call q2fvar(      q=  q(i,j,0:km,:),                 &
                      density=rho(i,j,0:km),                   &
                     velocity=vel(i,j,0:km,:),                 &
                     pressure=prs(i,j,0:km),                   &
                  temperature=tmp(i,j,0:km),                   &
                      species=spc(i,j,0:km,:)                  )
          !
        enddo
        deallocate(qfilt)
        !
      endif
      !
      allocate(Ecs(0:2,1:numq),dEcs(1:numq))
      ! do k=ks,ke
      ! do i=is,ie
      do k=0,km
      do i=0,im
        !
        pnor=pmatrix(rho(i,j,k),vel(i,j,k,1),vel(i,j,k,2),            &
                     vel(i,j,k,3),tmp(i,j,k),dxi(i,j,k,2,:),inv=.false.)
        pinv=pmatrix(rho(i,j,k),vel(i,j,k,1),vel(i,j,k,2),            &
                     vel(i,j,k,3),tmp(i,j,k),dxi(i,j,k,2,:),inv=.true.)
        !
        ! Pmult=MatMul(Pinv,pnor)
        ! if(irk==0) then
        !   print*,'---------------------------------------------------------'
        !   write(*,"(5(F7.4))")Pmult(1,:)
        !   write(*,"(5(F7.4))")Pmult(2,:)
        !   write(*,"(5(F7.4))")Pmult(3,:)
        !   write(*,"(5(F7.4))")Pmult(4,:)
        !   write(*,"(5(F7.4))")Pmult(5,:)
        ! end if
        !
        do ii=0,2
          uu=dxi(i,j+ii,k,2,1)*vel(i,j+ii,k,1) +                       &
             dxi(i,j+ii,k,2,2)*vel(i,j+ii,k,2) +                       &
             dxi(i,j+ii,k,2,3)*vel(i,j+ii,k,3)
          !
          Ecs(ii,1)=jacob(i,j+ii,k)*  q(i,j+ii,k,1)*uu
          Ecs(ii,2)=jacob(i,j+ii,k)*( q(i,j+ii,k,2)*uu+dxi(i,j+ii,k,2,1)*prs(i,j+ii,k) )
          Ecs(ii,3)=jacob(i,j+ii,k)*( q(i,j+ii,k,3)*uu+dxi(i,j+ii,k,2,2)*prs(i,j+ii,k) )
          Ecs(ii,4)=jacob(i,j+ii,k)*( q(i,j+ii,k,4)*uu+dxi(i,j+ii,k,2,3)*prs(i,j+ii,k) )
          Ecs(ii,5)=jacob(i,j+ii,k)*( q(i,j+ii,k,5)+prs(i,j+ii,k) )*uu
          do jspc=1,num_species
            Ecs(ii,5+jspc)=jacob(i,j+ii,k)*q(i,j+ii,k,5+jspc)*uu
          enddo
          !
          if(num_modequ>0) then
            n=5+num_species
            do jmod=1,num_modequ
              Ecs(ii,n+jmod)=jacob(i,j+ii,k)*q(i,j+ii,k,n+jmod)*uu
            enddo
          endif
          !
        enddo
        !
        do n=1,numq
          dEcs(n)= deriv( Ecs(0,n),Ecs(1,n),Ecs(2,n) )
        enddo
        !
        E(1)= q(i,j,k,2)
        E(2)= q(i,j,k,2)*vel(i,j,k,1)+prs(i,j,k)
        E(3)= q(i,j,k,3)*vel(i,j,k,1)
        E(4)= q(i,j,k,4)*vel(i,j,k,1)
        E(5)=(q(i,j,k,5)+prs(i,j,k))*vel(i,j,k,1)
        !
        F(1)= q(i,j,k,3)
        F(2)= q(i,j,k,2)*vel(i,j,k,2)
        F(3)= q(i,j,k,3)*vel(i,j,k,2)+prs(i,j,k)
        F(4)= q(i,j,k,4)*vel(i,j,k,2)
        F(5)=(q(i,j,k,5)+prs(i,j,k))*vel(i,j,k,2)
        !
        G(1)= q(i,j,k,4)
        G(2)= q(i,j,k,2)*vel(i,j,k,3)
        G(3)= q(i,j,k,3)*vel(i,j,k,3)
        G(4)= q(i,j,k,4)*vel(i,j,k,3)+prs(i,j,k)
        G(5)=(q(i,j,k,5)+prs(i,j,k))*vel(i,j,k,3)
        !
        jcbi(1)= deriv( dxi(i,j,k,2,1)  *jacob(i,j,k),                 &
                        dxi(i,j+1,k,2,1)*jacob(i,j+1,k),                 &
                        dxi(i,j+2,k,2,1)*jacob(i,j+2,k) )
        jcbi(2)= deriv( dxi(i,j,k,2,2)  *jacob(i,j,k),                 &
                        dxi(i,j+1,k,2,2)*jacob(i,j+1,k),                 &
                        dxi(i,j+2,k,2,2)*jacob(i,j+2,k) )
        jcbi(3)= deriv( dxi(i,j,k,2,3)  *jacob(i,j,k),                 &
                        dxi(i,j+1,k,2,3)*jacob(i,j+1,k),                 &
                        dxi(i,j+3,k,2,3)*jacob(i,j+3,k) )
        !
        Rest(1)=E(1)*Jcbi(1)+F(1)*Jcbi(2)+G(1)*Jcbi(3)
        Rest(2)=E(2)*Jcbi(1)+F(2)*Jcbi(2)+G(2)*Jcbi(3)
        Rest(3)=E(3)*Jcbi(1)+F(3)*Jcbi(2)+G(3)*Jcbi(3)
        Rest(4)=E(4)*Jcbi(1)+F(4)*Jcbi(2)+G(4)*Jcbi(3)
        Rest(5)=E(5)*Jcbi(1)+F(5)*Jcbi(2)+G(5)*Jcbi(3)
        !
        LODi1(1)=dEcs(1)-Rest(1)
        LODi1(2)=dEcs(2)-Rest(2)
        LODi1(3)=dEcs(3)-Rest(3)
        LODi1(4)=dEcs(4)-Rest(4)
        LODi1(5)=dEcs(5)-Rest(5)
        !
        LODi=MatMul(pinv,LODi1)/jacob(i,j,k)
        !
        css=sos(tmp(i,j,k))
        !
        uu=-(dxi(i,j,k,2,1)*vel(i,j,k,1) +                            &
             dxi(i,j,k,2,2)*vel(i,j,k,2) +                            &
             dxi(i,j,k,2,3)*vel(i,j,k,3))
        ! if(uu>=0.d0) then
          ! kinout=0.25d0*(1.d0-gmachmax2)*css/(ymax-ymin)
          ! LODi(4)=kinout*(pinf-prs(i,j,k))/rho(i,j,k)/css
          ! LODi(4)=kinout*(prs(i,j,k)-pinf)/rho(i,j,k)/css
        ! else
        !   var1=1.d0/sqrt( dxi(i,j,k,2,1)**2+dxi(i,j,k,2,2)**2+         &
        !                   dxi(i,j,k,2,3)**2 )
        !   kin=0.25d0*(1.d0-gmachmax2)*css/(ymax-ymin)
        !   LODi(1)=0.d0
        !   LODi(2)=kin*0.5d0*                                           &
        !                     ( dxi(i,j,k,2,3)*var1*(vel(i,j,k,2)-vinf)- &
        !                       dxi(i,j,k,2,2)*var1*(vel(i,j,k,3)-winf) )
        !   LODi(3)=kin*0.5d0*                                           &
        !                    (-dxi(i,j,k,2,1)*var1*(vel(i,j,k,2)-vinf)+  &
        !                      dxi(i,j,k,2,2)*var1*(vel(i,j,k,1)-uinf_j0) )
        !   LODi(4)=kin*(dxi(i,j,k,2,1)*var1*(vel(i,j,k,1)-uinf_j0)+        &
        !                dxi(i,j,k,2,2)*var1*(vel(i,j,k,2)-vinf)+        &
        !                dxi(i,j,k,2,3)*var1*(vel(i,j,k,3)-winf)+        &
        !                (prs(i,j,k)-pinf)/rho(i,j,k)/css )
        ! endif
        !
        LODi1=MatMul(pnor,LODi)*jacob(i,j,k)
        !
        dEcs(1)=LODi1(1)+Rest(1)
        dEcs(2)=LODi1(2)+Rest(2)
        dEcs(3)=LODi1(3)+Rest(3)
        dEcs(4)=LODi1(4)+Rest(4)
        dEcs(5)=LODi1(5)+Rest(5)
        !
        qrhs(i,j,k,:)=qrhs(i,j,k,:)+dEcs(:)
        !
      enddo
      enddo
      !
      deallocate(Ecs,dEcs)
      !
      allocate(fcs(-hm:im+hm,1:numq),dfcs(0:im,1:numq))
      do k=ks,ke
        !
        do i=-hm,im+hm
          !
          uu=dxi(i,j,k,1,1)*vel(i,j,k,1)+dxi(i,j,k,1,2)*vel(i,j,k,2) + &
             dxi(i,j,k,1,3)*vel(i,j,k,3)
          fcs(i,1)=jacob(i,j,k)*  q(i,j,k,1)*uu
          fcs(i,2)=jacob(i,j,k)*( q(i,j,k,2)*uu+dxi(i,j,k,1,1)*prs(i,j,k) )
          fcs(i,3)=jacob(i,j,k)*( q(i,j,k,3)*uu+dxi(i,j,k,1,2)*prs(i,j,k) )
          fcs(i,4)=jacob(i,j,k)*( q(i,j,k,4)*uu+dxi(i,j,k,1,3)*prs(i,j,k) )
          fcs(i,5)=jacob(i,j,k)*( q(i,j,k,5)+prs(i,j,k) )*uu
          do jspc=1,num_species
            fcs(i,5+jspc)=jacob(i,j,k)*q(i,j,k,5+jspc)*uu
          enddo
          !
        enddo
        !
        do n=1,numq
          dfcs(:,n)=ddfc(fcs(:,n),'222e',npdci,im)
        enddo
        !
        qrhs(is:ie,j,k,:)=qrhs(is:ie,j,k,:)+dfcs(is:ie,:)
        !
      enddo
      !
      deallocate(fcs,dfcs)
      ! !
      if(ndims==3) then
        !
        allocate(fcs(-hm:km+hm,1:numq),dfcs(0:km,1:numq))
        do i=is,ie
          !
          do k=-hm,km+hm
            !
            uu=dxi(i,j,k,3,1)*vel(i,j,k,1)+dxi(i,j,k,3,2)*vel(i,j,k,2) + &
               dxi(i,j,k,3,3)*vel(i,j,k,3)
            fcs(k,1)=jacob(i,j,k)*  q(i,j,k,1)*uu
            fcs(k,2)=jacob(i,j,k)*( q(i,j,k,2)*uu+dxi(i,j,k,3,1)*prs(i,j,k) )
            fcs(k,3)=jacob(i,j,k)*( q(i,j,k,3)*uu+dxi(i,j,k,3,2)*prs(i,j,k) )
            fcs(k,4)=jacob(i,j,k)*( q(i,j,k,4)*uu+dxi(i,j,k,3,3)*prs(i,j,k) )
            fcs(k,5)=jacob(i,j,k)*( q(i,j,k,5)+prs(i,j,k) )*uu
            do jspc=1,num_species
              fcs(k,5+jspc)=jacob(i,j,k)*q(i,j,k,5+jspc)*uu
            enddo
            !
          enddo
          !
          do n=1,numq
            dfcs(:,n)=ddfc(fcs(:,n),'222e',npdck,km)
          enddo
          !
          qrhs(i,j,ks:ke,:)=qrhs(i,j,ks:ke,:)+dfcs(ks:ke,:)
          !
        enddo
        !
        deallocate(fcs,dfcs)
      endif
      !
    endif
    !
    if(ndir==4 .and. jrk==jrkm) then
      !
      j=jm
      !
      allocate(qfilt(0:im,1:numq))
      do k=0,km
        !
        do jq=1,numq
          qfilt(:,jq)=spafilter6exp(q(:,j,k,jq),npdci,im)
          q(0:im,j,k,jq)=qfilt(:,jq)
        enddo
        call q2fvar(      q=  q(0:im,j,k,:),                   &
                    density=rho(0:im,j,k),                     &
                   velocity=vel(0:im,j,k,:),                   &
                   pressure=prs(0:im,j,k),                     &
                temperature=tmp(0:im,j,k),                     &
                    species=spc(0:im,j,k,:)                    )
        !
      enddo
      deallocate(qfilt)
      !
      if(ndims==3) then 
        !
        allocate(qfilt(0:km,1:numq))
        do i=0,im
          !
          do jq=1,numq
            qfilt(:,jq)=spafilter6exp(q(i,j,:,jq),npdck,km)
            q(i,j,0:km,jq)=qfilt(:,jq)
          enddo
          call q2fvar(      q=  q(i,j,0:km,:),                 &
                      density=rho(i,j,0:km),                   &
                     velocity=vel(i,j,0:km,:),                 &
                     pressure=prs(i,j,0:km),                   &
                  temperature=tmp(i,j,0:km),                   &
                      species=spc(i,j,0:km,:)                  )
          !
        enddo
        deallocate(qfilt)
        !
      endif
      !
      allocate(Ecs(0:2,1:numq),dEcs(1:numq))
      ! do k=ks,ke
      ! do i=is,ie
      do k=0,km
      do i=0,im
        !
        pnor=pmatrix(rho(i,j,k),vel(i,j,k,1),vel(i,j,k,2),            &
                     vel(i,j,k,3),tmp(i,j,k),dxi(i,j,k,2,:),inv=.false.)
        pinv=pmatrix(rho(i,j,k),vel(i,j,k,1),vel(i,j,k,2),            &
                     vel(i,j,k,3),tmp(i,j,k),dxi(i,j,k,2,:),inv=.true.)
        !
        ! Pmult=MatMul(Pinv,pnor)
        ! if(irk==0) then
        !   print*,'---------------------------------------------------------'
        !   write(*,"(5(F7.4))")Pmult(1,:)
        !   write(*,"(5(F7.4))")Pmult(2,:)
        !   write(*,"(5(F7.4))")Pmult(3,:)
        !   write(*,"(5(F7.4))")Pmult(4,:)
        !   write(*,"(5(F7.4))")Pmult(5,:)
        ! end if
        !
        do ii=0,2
          !
          uu=dxi(i,j-ii,k,2,1)*vel(i,j-ii,k,1) +                       &
             dxi(i,j-ii,k,2,2)*vel(i,j-ii,k,2) +                       &
             dxi(i,j-ii,k,2,3)*vel(i,j-ii,k,3)
          !
          Ecs(ii,1)=jacob(i,j-ii,k)*  q(i,j-ii,k,1)*uu
          Ecs(ii,2)=jacob(i,j-ii,k)*( q(i,j-ii,k,2)*uu+dxi(i,j-ii,k,2,1)*prs(i,j-ii,k) )
          Ecs(ii,3)=jacob(i,j-ii,k)*( q(i,j-ii,k,3)*uu+dxi(i,j-ii,k,2,2)*prs(i,j-ii,k) )
          Ecs(ii,4)=jacob(i,j-ii,k)*( q(i,j-ii,k,4)*uu+dxi(i,j-ii,k,2,3)*prs(i,j-ii,k) )
          Ecs(ii,5)=jacob(i,j-ii,k)*( q(i,j-ii,k,5)+prs(i,j-ii,k) )*uu
          do jspc=1,num_species
            Ecs(ii,5+jspc)=jacob(i,j-ii,k)*q(i,j-ii,k,5+jspc)*uu
          enddo
          !
          if(num_modequ>0) then
            n=5+num_species
            do jmod=1,num_modequ
              Ecs(ii,n+jmod)=jacob(i,j-ii,k)*q(i,j-ii,k,n+jmod)*uu
            enddo
          endif
          !
        enddo
        !
        do n=1,numq
          dEcs(n)=-deriv( Ecs(0,n),Ecs(1,n),Ecs(2,n) )
        enddo
        !
        E(1)= q(i,j,k,2)
        E(2)= q(i,j,k,2)*vel(i,j,k,1)+prs(i,j,k)
        E(3)= q(i,j,k,3)*vel(i,j,k,1)
        E(4)= q(i,j,k,4)*vel(i,j,k,1)
        E(5)=(q(i,j,k,5)+prs(i,j,k))*vel(i,j,k,1)
        !
        F(1)= q(i,j,k,3)
        F(2)= q(i,j,k,2)*vel(i,j,k,2)
        F(3)= q(i,j,k,3)*vel(i,j,k,2)+prs(i,j,k)
        F(4)= q(i,j,k,4)*vel(i,j,k,2)
        F(5)=(q(i,j,k,5)+prs(i,j,k))*vel(i,j,k,2)
        !
        G(1)= q(i,j,k,4)
        G(2)= q(i,j,k,2)*vel(i,j,k,3)
        G(3)= q(i,j,k,3)*vel(i,j,k,3)
        G(4)= q(i,j,k,4)*vel(i,j,k,3)+prs(i,j,k)
        G(5)=(q(i,j,k,5)+prs(i,j,k))*vel(i,j,k,3)
        !
        jcbi(1)=-deriv( dxi(i,j,k,2,1)  *jacob(i,j,k),                 &
                        dxi(i,j-1,k,2,1)*jacob(i,j-1,k),               &
                        dxi(i,j-2,k,2,1)*jacob(i,j-2,k) )
        jcbi(2)=-deriv( dxi(i,j,k,2,2)  *jacob(i,j,k),                 &
                        dxi(i,j-1,k,2,2)*jacob(i,j-1,k) ,              &
                        dxi(i,j-2,k,2,2)*jacob(i,j-2,k))
        jcbi(3)=-deriv( dxi(i,j,k,2,3)  *jacob(i,j,k),                 &
                        dxi(i,j-1,k,2,3)*jacob(i,j-1,k),               &
                        dxi(i,j-2,k,2,3)*jacob(i,j-2,k) )
        !
        Rest(1)=E(1)*Jcbi(1)+F(1)*Jcbi(2)+G(1)*Jcbi(3)
        Rest(2)=E(2)*Jcbi(1)+F(2)*Jcbi(2)+G(2)*Jcbi(3)
        Rest(3)=E(3)*Jcbi(1)+F(3)*Jcbi(2)+G(3)*Jcbi(3)
        Rest(4)=E(4)*Jcbi(1)+F(4)*Jcbi(2)+G(4)*Jcbi(3)
        Rest(5)=E(5)*Jcbi(1)+F(5)*Jcbi(2)+G(5)*Jcbi(3)
        !
        LODi1(1)=dEcs(1)-Rest(1)
        LODi1(2)=dEcs(2)-Rest(2)
        LODi1(3)=dEcs(3)-Rest(3)
        LODi1(4)=dEcs(4)-Rest(4)
        LODi1(5)=dEcs(5)-Rest(5)
        !
        LODi=MatMul(pinv,LODi1)/jacob(i,j,k)
        !
        ! if(irk==0 .and. i==0) then
        !   print*,LODi(2),LODi(3),LODi(5)
        ! endif
        !
        css=sos(tmp(i,j,k))
        !
        ! uu=dxi(i,j,k,2,1)*vel(i,j,k,1) +                       &
        !    dxi(i,j,k,2,2)*vel(i,j,k,2) +                       &
        !    dxi(i,j,k,2,3)*vel(i,j,k,3)
        ! if(uu>=0.d0) then
          kinout=0.25d0*(1.d0-gmachmax2)*css/(ymax-ymin)
          LODi(5)=kinout*(prs(i,j,k)-pinf)/rho(i,j,k)/css
          ! LODi(5)=kinout*(prs(i,j,k)-prs_prof(jm))/rho(i,j,k)/css
        ! else
        !   var1=1.d0/sqrt( dxi(i,j,k,2,1)**2+dxi(i,j,k,2,2)**2+         &
        !                   dxi(i,j,k,2,3)**2 )
        !   kin=-0.250*(1.d0-gmachmax2)*css/(ymax-ymin)
        !   !
        !   ! kin=1.d0
        !   !
        !   LODi(1)=0.d0 !0.25d0*rho(i,j,k)*css*css/gamma/(ymax-ymin)*(tmp(i,j,k)-tinf)
        !   !
        !   LODi(2)=kin*0.5d0*                                           &
        !                     ( dxi(i,j,k,2,3)*var1*(vel(i,j,k,2)-vinf)- &
        !                       dxi(i,j,k,2,2)*var1*(vel(i,j,k,3)-winf) )
        !   LODi(3)=kin*0.5d0*                                           &
        !                    (-dxi(i,j,k,2,1)*var1*(vel(i,j,k,2)-vinf)+  &
        !                      dxi(i,j,k,2,2)*var1*(vel(i,j,k,1)-uinf_jm) )
        !   LODi(5)=kin*(dxi(i,j,k,2,1)*var1*(vel(i,j,k,1)-uinf_jm)+     &
        !                dxi(i,j,k,2,2)*var1*(vel(i,j,k,2)-vinf)+        &
        !                dxi(i,j,k,2,3)*var1*(vel(i,j,k,3)-winf)+        &
        !                (prs(i,j,k)-pinf)/rho(i,j,k)/css )
        !   ! if(irk==0 .and. i==0) then
        !   !   print*,dxi(i,j,k,2,1)*var1,gmachmax2
        !   ! endif
        !   !
        ! endif
        ! LODi(5)=0.d0
        !
        LODi1=MatMul(pnor,LODi)*jacob(i,j,k)
        !
        dEcs(1)=LODi1(1)+Rest(1)
        dEcs(2)=LODi1(2)+Rest(2)
        dEcs(3)=LODi1(3)+Rest(3)
        dEcs(4)=LODi1(4)+Rest(4)
        dEcs(5)=LODi1(5)+Rest(5)
        !
        qrhs(i,j,k,:)=qrhs(i,j,k,:)+dEcs(:)
        ! !
        ! if(irk==0 .and. i==0) then
        !   print*,LODi(2),LODi(3),LODi(5)
        !   print*,'-------------------------'
        ! endif
        !
      enddo
      enddo
      !
      !
      deallocate(Ecs,dEcs)
      !
      allocate(fcs(-hm:im+hm,1:numq),dfcs(0:im,1:numq))
      do k=ks,ke
        !
        do i=-hm,im+hm
          !
          uu=dxi(i,j,k,1,1)*vel(i,j,k,1)+dxi(i,j,k,1,2)*vel(i,j,k,2) + &
             dxi(i,j,k,1,3)*vel(i,j,k,3)
          fcs(i,1)=jacob(i,j,k)*  q(i,j,k,1)*uu
          fcs(i,2)=jacob(i,j,k)*( q(i,j,k,2)*uu+dxi(i,j,k,1,1)*prs(i,j,k) )
          fcs(i,3)=jacob(i,j,k)*( q(i,j,k,3)*uu+dxi(i,j,k,1,2)*prs(i,j,k) )
          fcs(i,4)=jacob(i,j,k)*( q(i,j,k,4)*uu+dxi(i,j,k,1,3)*prs(i,j,k) )
          fcs(i,5)=jacob(i,j,k)*( q(i,j,k,5)+prs(i,j,k) )*uu
          do jspc=1,num_species
            fcs(i,5+jspc)=jacob(i,j,k)*q(i,j,k,5+jspc)*uu
          enddo
          !
        enddo
        !
        do n=1,numq
          dfcs(:,n)=ddfc(fcs(:,n),'222e',npdci,im)
        enddo
        !
        qrhs(is:ie,j,k,:)=qrhs(is:ie,j,k,:)+dfcs(is:ie,:)
        !
      enddo
      !
      deallocate(fcs,dfcs)
      !
      if(ndims==3) then
        !
        allocate(fcs(-hm:km+hm,1:numq),dfcs(0:km,1:numq))
        do i=is,ie
          !
          do k=-hm,km+hm
            !
            uu=dxi(i,j,k,3,1)*vel(i,j,k,1)+dxi(i,j,k,3,2)*vel(i,j,k,2) + &
               dxi(i,j,k,3,3)*vel(i,j,k,3)
            fcs(k,1)=jacob(i,j,k)*  q(i,j,k,1)*uu
            fcs(k,2)=jacob(i,j,k)*( q(i,j,k,2)*uu+dxi(i,j,k,3,1)*prs(i,j,k) )
            fcs(k,3)=jacob(i,j,k)*( q(i,j,k,3)*uu+dxi(i,j,k,3,2)*prs(i,j,k) )
            fcs(k,4)=jacob(i,j,k)*( q(i,j,k,4)*uu+dxi(i,j,k,3,3)*prs(i,j,k) )
            fcs(k,5)=jacob(i,j,k)*( q(i,j,k,5)+prs(i,j,k) )*uu
            do jspc=1,num_species
              fcs(k,5+jspc)=jacob(i,j,k)*q(i,j,k,5+jspc)*uu
            enddo
            !
          enddo
          !
          do n=1,numq
            dfcs(:,n)=ddfc(fcs(:,n),'222e',npdck,km)
          enddo
          !
          qrhs(i,j,ks:ke,:)=qrhs(i,j,ks:ke,:)+dfcs(ks:ke,:)
          !
        enddo
        !
        deallocate(fcs,dfcs)
        !
      endif
      !
    endif
    !
    if(ndir==5 .and. krk==0) then
      !
      k=0
      !
      allocate(qfilt(0:im,1:numq))
      do j=0,jm
        !
        do jq=1,numq
          qfilt(:,jq)=spafilter6exp(q(:,j,k,jq),npdci,im)
          q(0:im,j,k,jq)=qfilt(:,jq)
        enddo
        call q2fvar(      q=  q(0:im,j,k,:),                   &
                    density=rho(0:im,j,k),                     &
                   velocity=vel(0:im,j,k,:),                   &
                   pressure=prs(0:im,j,k),                     &
                temperature=tmp(0:im,j,k),                     &
                    species=spc(0:im,j,k,:)                    )
        !
      enddo
      deallocate(qfilt)
      !
      allocate(qfilt(0:jm,1:numq))
      do i=0,im
        !
        do jq=1,numq
          qfilt(:,jq)=spafilter6exp(q(i,:,k,jq),npdcj,jm)
          q(i,0:jm,k,jq)=qfilt(:,jq)
        enddo
        call q2fvar(      q=  q(i,0:jm,k,:),                   &
                    density=rho(i,0:jm,k),                     &
                   velocity=vel(i,0:jm,k,:),                   &
                   pressure=prs(i,0:jm,k),                     &
                temperature=tmp(i,0:jm,k),                     &
                    species=spc(i,0:jm,k,:)                    )
        !
      enddo
      deallocate(qfilt)
      !
      allocate(Ecs(0:2,1:numq),dEcs(1:numq))
      do j=0,jm
      do i=0,im
        ! do j=js,je
        ! do i=is,ie
        !
        pnor=pmatrix(rho(i,j,k),vel(i,j,k,1),vel(i,j,k,2),            &
                     vel(i,j,k,3),tmp(i,j,k),dxi(i,j,k,3,:),inv=.false.)
        pinv=pmatrix(rho(i,j,k),vel(i,j,k,1),vel(i,j,k,2),            &
                     vel(i,j,k,3),tmp(i,j,k),dxi(i,j,k,3,:),inv=.true.)
        !
        ! Pmult=MatMul(Pinv,pnor)
        ! if(irk==0 .and. jrk==0) then
        !   print*,'---------------------------------------------------------'
        !   write(*,"(5(F7.4))")Pmult(1,:)
        !   write(*,"(5(F7.4))")Pmult(2,:)
        !   write(*,"(5(F7.4))")Pmult(3,:)
        !   write(*,"(5(F7.4))")Pmult(4,:)
        !   write(*,"(5(F7.4))")Pmult(5,:)
        ! end if
        !
        do ii=0,2
          uu=dxi(i,j,k+ii,3,1)*vel(i,j,k+ii,1) +                       &
             dxi(i,j,k+ii,3,2)*vel(i,j,k+ii,2) +                       &
             dxi(i,j,k+ii,3,3)*vel(i,j,k+ii,3)
          !
          Ecs(ii,1)=jacob(i,j,k+ii)*  q(i,j,k+ii,1)*uu
          Ecs(ii,2)=jacob(i,j,k+ii)*( q(i,j,k+ii,2)*uu+dxi(i,j,k+ii,3,1)*prs(i,j,k+ii) )
          Ecs(ii,3)=jacob(i,j,k+ii)*( q(i,j,k+ii,3)*uu+dxi(i,j,k+ii,3,2)*prs(i,j,k+ii) )
          Ecs(ii,4)=jacob(i,j,k+ii)*( q(i,j,k+ii,4)*uu+dxi(i,j,k+ii,3,3)*prs(i,j,k+ii) )
          Ecs(ii,5)=jacob(i,j,k+ii)*( q(i,j,k+ii,5)+prs(i,j,k+ii) )*uu
          do jspc=1,num_species
            Ecs(ii,5+jspc)=jacob(i,j,k+ii)*q(i,j,k+ii,5+jspc)*uu
          enddo
        enddo
        !
        do n=1,numq
          dEcs(n)= deriv( Ecs(0,n),Ecs(1,n),Ecs(2,n) )
        enddo
        !
        E(1)= q(i,j,k,2)
        E(2)= q(i,j,k,2)*vel(i,j,k,1)+prs(i,j,k)
        E(3)= q(i,j,k,3)*vel(i,j,k,1)
        E(4)= q(i,j,k,4)*vel(i,j,k,1)
        E(5)=(q(i,j,k,5)+prs(i,j,k))*vel(i,j,k,1)
        !
        F(1)= q(i,j,k,3)
        F(2)= q(i,j,k,2)*vel(i,j,k,2)
        F(3)= q(i,j,k,3)*vel(i,j,k,2)+prs(i,j,k)
        F(4)= q(i,j,k,4)*vel(i,j,k,2)
        F(5)=(q(i,j,k,5)+prs(i,j,k))*vel(i,j,k,2)
        !
        G(1)= q(i,j,k,4)
        G(2)= q(i,j,k,2)*vel(i,j,k,3)
        G(3)= q(i,j,k,3)*vel(i,j,k,3)
        G(4)= q(i,j,k,4)*vel(i,j,k,3)+prs(i,j,k)
        G(5)=(q(i,j,k,5)+prs(i,j,k))*vel(i,j,k,3)
        !
        jcbi(1)= deriv( dxi(i,j,k,  3,1)*jacob(i,j,k),                 &
                        dxi(i,j,k+1,3,1)*jacob(i,j,k+1),               &
                        dxi(i,j,k+2,3,1)*jacob(i,j,k+2) )
        jcbi(2)= deriv( dxi(i,j,k,  3,2)*jacob(i,j,k),                 &
                        dxi(i,j,k+1,3,2)*jacob(i,j,k+1),               &
                        dxi(i,j,k+2,3,2)*jacob(i,j,k+2) )
        jcbi(3)= deriv( dxi(i,j,k,  3,3)*jacob(i,j,k),                 &
                        dxi(i,j,k+1,3,3)*jacob(i,j,k+1),               &
                        dxi(i,j,k+2,3,3)*jacob(i,j,k+2) )
        !
        Rest(1)=E(1)*Jcbi(1)+F(1)*Jcbi(2)+G(1)*Jcbi(3)
        Rest(2)=E(2)*Jcbi(1)+F(2)*Jcbi(2)+G(2)*Jcbi(3)
        Rest(3)=E(3)*Jcbi(1)+F(3)*Jcbi(2)+G(3)*Jcbi(3)
        Rest(4)=E(4)*Jcbi(1)+F(4)*Jcbi(2)+G(4)*Jcbi(3)
        Rest(5)=E(5)*Jcbi(1)+F(5)*Jcbi(2)+G(5)*Jcbi(3)
        !
        LODi1(1)=dEcs(1)-Rest(1)
        LODi1(2)=dEcs(2)-Rest(2)
        LODi1(3)=dEcs(3)-Rest(3)
        LODi1(4)=dEcs(4)-Rest(4)
        LODi1(5)=dEcs(5)-Rest(5)
        !
        LODi=MatMul(pinv,LODi1)/jacob(i,j,k)
        !
        css=sos(tmp(i,j,k))
        !
        ! uu=-(dxi(i,j,k,3,1)*vel(i,j,k,1) +                            &
        !      dxi(i,j,k,3,2)*vel(i,j,k,2) +                            &
        !      dxi(i,j,k,3,3)*vel(i,j,k,3))
        ! if(uu>=0.d0) then
          kinout=0.25d0*(1.d0-gmachmax2)*css/(zmax-zmin)
          LODi(4)=kinout*(prs(i,j,k)-pinf)/rho(i,j,k)/css
        ! else
        !   var1=1.d0/sqrt( dxi(i,j,k,3,1)**2+dxi(i,j,k,3,2)**2+         &
        !                   dxi(i,j,k,3,3)**2 )
        !   kin=0.25d0*(1.d0-gmachmax2)*css/(zmax-zmin)
        !   LODi(1)=0.d0
        !   LODi(2)=kin*0.5d0*                                           &
        !                     (-dxi(i,j,k,3,3)*var1*(vel(i,j,k,2)-vinf)+ &
        !                       dxi(i,j,k,3,2)*var1*(vel(i,j,k,3)-winf) )
        !   LODi(3)=kin*0.5d0*                                           &
        !                    ( dxi(i,j,k,3,3)*var1*(vel(i,j,k,1)-uinf)+  &
        !                      dxi(i,j,k,3,1)*var1*(vel(i,j,k,3)-winf) )
        !   LODi(4)=kin*(dxi(i,j,k,3,1)*var1*(vel(i,j,k,1)-uinf)+        &
        !                dxi(i,j,k,3,2)*var1*(vel(i,j,k,2)-vinf)+        &
        !                dxi(i,j,k,3,3)*var1*(vel(i,j,k,3)-winf)+        &
        !                (prs(i,j,k)-pinf)/rho(i,j,k)/css )
        ! endif
        !
        LODi1=MatMul(pnor,LODi)*jacob(i,j,k)
        !
        dEcs(1)=LODi1(1)+Rest(1)
        dEcs(2)=LODi1(2)+Rest(2)
        dEcs(3)=LODi1(3)+Rest(3)
        dEcs(4)=LODi1(4)+Rest(4)
        dEcs(5)=LODi1(5)+Rest(5)
        !
        qrhs(i,j,k,:)=qrhs(i,j,k,:)+dEcs(:)
        !
      enddo
      enddo
      !
      deallocate(Ecs,dEcs)
      !
      allocate(fcs(-hm:im+hm,1:numq),dfcs(0:im,1:numq))
      do j=js,je
        !
        do i=-hm,im+hm
          !
          uu=dxi(i,j,k,1,1)*vel(i,j,k,1)+dxi(i,j,k,1,2)*vel(i,j,k,2) + &
             dxi(i,j,k,1,3)*vel(i,j,k,3)
          fcs(i,1)=jacob(i,j,k)*  q(i,j,k,1)*uu
          fcs(i,2)=jacob(i,j,k)*( q(i,j,k,2)*uu+dxi(i,j,k,1,1)*prs(i,j,k) )
          fcs(i,3)=jacob(i,j,k)*( q(i,j,k,3)*uu+dxi(i,j,k,1,2)*prs(i,j,k) )
          fcs(i,4)=jacob(i,j,k)*( q(i,j,k,4)*uu+dxi(i,j,k,1,3)*prs(i,j,k) )
          fcs(i,5)=jacob(i,j,k)*( q(i,j,k,5)+prs(i,j,k) )*uu
          do jspc=1,num_species
            fcs(i,5+jspc)=jacob(i,j,k)*q(i,j,k,5+jspc)*uu
          enddo
          !
        enddo
        !
        do n=1,numq
          dfcs(:,n)=ddfc(fcs(:,n),'222e',npdci,im)
        enddo
        !
        qrhs(is:ie,j,k,:)=qrhs(is:ie,j,k,:)+dfcs(is:ie,:)
        !
      enddo
      !
      deallocate(fcs,dfcs)
      !
      allocate(fcs(-hm:jm+hm,1:numq),dfcs(0:jm,1:numq))
      do i=is,ie
        !
        do j=-hm,jm+hm
          !
          uu=dxi(i,j,k,2,1)*vel(i,j,k,1)+dxi(i,j,k,2,2)*vel(i,j,k,2) + &
             dxi(i,j,k,2,3)*vel(i,j,k,3)
          fcs(j,1)=jacob(i,j,k)*  q(i,j,k,1)*uu
          fcs(j,2)=jacob(i,j,k)*( q(i,j,k,2)*uu+dxi(i,j,k,2,1)*prs(i,j,k) )
          fcs(j,3)=jacob(i,j,k)*( q(i,j,k,3)*uu+dxi(i,j,k,2,2)*prs(i,j,k) )
          fcs(j,4)=jacob(i,j,k)*( q(i,j,k,4)*uu+dxi(i,j,k,2,3)*prs(i,j,k) )
          fcs(j,5)=jacob(i,j,k)*( q(i,j,k,5)+prs(i,j,k) )*uu
          do jspc=1,num_species
            fcs(j,5+jspc)=jacob(i,j,k)*q(i,j,k,5+jspc)*uu
          enddo
          !
        enddo
        !
        do n=1,numq
          dfcs(:,n)=ddfc(fcs(:,n),'222e',npdcj,jm)
        enddo
        !
        qrhs(i,js:je,k,:)=qrhs(i,js:je,k,:)+dfcs(js:je,:)
        !
      enddo
      !
      deallocate(fcs,dfcs)
      !
    endif
    !
    if(ndir==6 .and. krk==krkm) then
      !
      k=km
      !
      allocate(qfilt(0:im,1:numq))
      do j=0,jm
        !
        do jq=1,numq
          qfilt(:,jq)=spafilter6exp(q(:,j,k,jq),npdci,im)
          q(0:im,j,k,jq)=qfilt(:,jq)
        enddo
        call q2fvar(      q=  q(0:im,j,k,:),                   &
                    density=rho(0:im,j,k),                     &
                   velocity=vel(0:im,j,k,:),                   &
                   pressure=prs(0:im,j,k),                     &
                temperature=tmp(0:im,j,k),                     &
                    species=spc(0:im,j,k,:)                    )
        !
      enddo
      deallocate(qfilt)
      !
      allocate(qfilt(0:jm,1:numq))
      do i=0,im
        !
        do jq=1,numq
          qfilt(:,jq)=spafilter6exp(q(i,:,k,jq),npdcj,jm)
          q(i,0:jm,k,jq)=qfilt(:,jq)
        enddo
        call q2fvar(      q=  q(i,0:jm,k,:),                   &
                    density=rho(i,0:jm,k),                     &
                   velocity=vel(i,0:jm,k,:),                   &
                   pressure=prs(i,0:jm,k),                     &
                temperature=tmp(i,0:jm,k),                     &
                    species=spc(i,0:jm,k,:)                    )
        !
      enddo
      deallocate(qfilt)
      !
      allocate(Ecs(0:2,1:numq),dEcs(1:numq))
      do j=0,jm
      do i=0,im
      ! do j=js,je
      ! do i=is,ie
        !
        pnor=pmatrix(rho(i,j,k),vel(i,j,k,1),vel(i,j,k,2),            &
                     vel(i,j,k,3),tmp(i,j,k),dxi(i,j,k,3,:),inv=.false.)
        pinv=pmatrix(rho(i,j,k),vel(i,j,k,1),vel(i,j,k,2),            &
                     vel(i,j,k,3),tmp(i,j,k),dxi(i,j,k,3,:),inv=.true.)
        !
        ! Pmult=MatMul(Pinv,pnor)
        ! if(irk==0) then
        !   print*,'-----------------------------------------------------'
        !   write(*,"(5(F7.4))")Pmult(1,:)
        !   write(*,"(5(F7.4))")Pmult(2,:)
        !   write(*,"(5(F7.4))")Pmult(3,:)
        !   write(*,"(5(F7.4))")Pmult(4,:)
        !   write(*,"(5(F7.4))")Pmult(5,:)
        ! end if
        !
        do ii=0,2
          uu=dxi(i,j,k-ii,3,1)*vel(i,j,k-ii,1) +                       &
             dxi(i,j,k-ii,3,2)*vel(i,j,k-ii,2) +                       &
             dxi(i,j,k-ii,3,3)*vel(i,j,k-ii,3)
          !
          Ecs(ii,1)=jacob(i,j,k-ii)*  q(i,j,k-ii,1)*uu
          Ecs(ii,2)=jacob(i,j,k-ii)*( q(i,j,k-ii,2)*uu+dxi(i,j,k-ii,3,1)*prs(i,j,k-ii) )
          Ecs(ii,3)=jacob(i,j,k-ii)*( q(i,j,k-ii,3)*uu+dxi(i,j,k-ii,3,2)*prs(i,j,k-ii) )
          Ecs(ii,4)=jacob(i,j,k-ii)*( q(i,j,k-ii,4)*uu+dxi(i,j,k-ii,3,3)*prs(i,j,k-ii) )
          Ecs(ii,5)=jacob(i,j,k-ii)*( q(i,j,k-ii,5)+prs(i,j,k-ii) )*uu
          do jspc=1,num_species
            Ecs(ii,5+jspc)=jacob(i,j,k-ii)*q(i,j,k-ii,5+jspc)*uu
          enddo
        enddo
        !
        do n=1,numq
          dEcs(n)=-deriv( Ecs(0,n),Ecs(1,n),Ecs(2,n) )
        enddo
        !
        E(1)= q(i,j,k,2)
        E(2)= q(i,j,k,2)*vel(i,j,k,1)+prs(i,j,k)
        E(3)= q(i,j,k,3)*vel(i,j,k,1)
        E(4)= q(i,j,k,4)*vel(i,j,k,1)
        E(5)=(q(i,j,k,5)+prs(i,j,k))*vel(i,j,k,1)
        !
        F(1)= q(i,j,k,3)
        F(2)= q(i,j,k,2)*vel(i,j,k,2)
        F(3)= q(i,j,k,3)*vel(i,j,k,2)+prs(i,j,k)
        F(4)= q(i,j,k,4)*vel(i,j,k,2)
        F(5)=(q(i,j,k,5)+prs(i,j,k))*vel(i,j,k,2)
        !
        G(1)= q(i,j,k,4)
        G(2)= q(i,j,k,2)*vel(i,j,k,3)
        G(3)= q(i,j,k,3)*vel(i,j,k,3)
        G(4)= q(i,j,k,4)*vel(i,j,k,3)+prs(i,j,k)
        G(5)=(q(i,j,k,5)+prs(i,j,k))*vel(i,j,k,3)
        !
        jcbi(1)=-deriv( dxi(i,j,k,  3,1)*jacob(i,j,k),                 &
                        dxi(i,j,k-1,3,1)*jacob(i,j,k-1),               &
                        dxi(i,j,k-2,3,1)*jacob(i,j,k-2) )
        jcbi(2)=-deriv( dxi(i,j,k,  3,2)*jacob(i,j,k),                 &
                        dxi(i,j,k-1,3,2)*jacob(i,j,k-1),               &
                        dxi(i,j,k-2,3,2)*jacob(i,j,k-2) )
        jcbi(3)=-deriv( dxi(i,j,k,  3,3)*jacob(i,j,k),                 &
                        dxi(i,j,k-1,3,3)*jacob(i,j,k-1),               &
                        dxi(i,j,k-2,3,3)*jacob(i,j,k-2) )
        !
        Rest(1)=E(1)*Jcbi(1)+F(1)*Jcbi(2)+G(1)*Jcbi(3)
        Rest(2)=E(2)*Jcbi(1)+F(2)*Jcbi(2)+G(2)*Jcbi(3)
        Rest(3)=E(3)*Jcbi(1)+F(3)*Jcbi(2)+G(3)*Jcbi(3)
        Rest(4)=E(4)*Jcbi(1)+F(4)*Jcbi(2)+G(4)*Jcbi(3)
        Rest(5)=E(5)*Jcbi(1)+F(5)*Jcbi(2)+G(5)*Jcbi(3)
        !
        LODi1(1)=dEcs(1)-Rest(1)
        LODi1(2)=dEcs(2)-Rest(2)
        LODi1(3)=dEcs(3)-Rest(3)
        LODi1(4)=dEcs(4)-Rest(4)
        LODi1(5)=dEcs(5)-Rest(5)
        !
        LODi=MatMul(pinv,LODi1)/jacob(i,j,k)
        !
        css=sos(tmp(i,j,k))
        ! !
        ! uu= (dxi(i,j,k,3,1)*vel(i,j,k,1) +                            &
        !      dxi(i,j,k,3,2)*vel(i,j,k,2) +                            &
        !      dxi(i,j,k,3,3)*vel(i,j,k,3))
        ! if(uu>=0.d0) then
          kinout=0.25d0*(1.d0-gmachmax2)*css/(zmax-zmin)
          LODi(5)=kinout*(prs(i,j,k)-pinf)/rho(i,j,k)/css
        ! else
        !   var1=1.d0/sqrt( dxi(i,j,k,3,1)**2+dxi(i,j,k,3,2)**2+         &
        !                   dxi(i,j,k,3,3)**2 )
        !   kin=-0.25d0*(1.d0-gmachmax2)*css/(zmax-zmin)
        !   LODi(1)=0.d0
        !   LODi(2)=kin*0.5d0*                                           &
        !                     (-dxi(i,j,k,3,3)*var1*(vel(i,j,k,2)-vinf)+ &
        !                       dxi(i,j,k,3,2)*var1*(vel(i,j,k,3)-winf) )
        !   LODi(3)=kin*0.5d0*                                           &
        !                    ( dxi(i,j,k,3,3)*var1*(vel(i,j,k,1)-uinf)+  &
        !                      dxi(i,j,k,3,1)*var1*(vel(i,j,k,3)-winf) )
        !   LODi(5)=kin*(dxi(i,j,k,3,1)*var1*(vel(i,j,k,1)-uinf)+        &
        !                dxi(i,j,k,3,2)*var1*(vel(i,j,k,2)-vinf)+        &
        !                dxi(i,j,k,3,3)*var1*(vel(i,j,k,3)-winf)+        &
        !                (prs(i,j,k)-pinf)/rho(i,j,k)/css )
        ! endif
        !
        LODi1=MatMul(pnor,LODi)*jacob(i,j,k)
        !
        dEcs(1)=LODi1(1)+Rest(1)
        dEcs(2)=LODi1(2)+Rest(2)
        dEcs(3)=LODi1(3)+Rest(3)
        dEcs(4)=LODi1(4)+Rest(4)
        dEcs(5)=LODi1(5)+Rest(5)
        !
        qrhs(i,j,k,:)=qrhs(i,j,k,:)+dEcs(:)
        !
      enddo
      enddo
      !
      deallocate(Ecs,dEcs)
      !
      allocate(fcs(-hm:im+hm,1:numq),dfcs(0:im,1:numq))
      do j=js,je
        !
        do i=-hm,im+hm
          !
          uu=dxi(i,j,k,1,1)*vel(i,j,k,1)+dxi(i,j,k,1,2)*vel(i,j,k,2) + &
             dxi(i,j,k,1,3)*vel(i,j,k,3)
          fcs(i,1)=jacob(i,j,k)*  q(i,j,k,1)*uu
          fcs(i,2)=jacob(i,j,k)*( q(i,j,k,2)*uu+dxi(i,j,k,1,1)*prs(i,j,k) )
          fcs(i,3)=jacob(i,j,k)*( q(i,j,k,3)*uu+dxi(i,j,k,1,2)*prs(i,j,k) )
          fcs(i,4)=jacob(i,j,k)*( q(i,j,k,4)*uu+dxi(i,j,k,1,3)*prs(i,j,k) )
          fcs(i,5)=jacob(i,j,k)*( q(i,j,k,5)+prs(i,j,k) )*uu
          do jspc=1,num_species
            fcs(i,5+jspc)=jacob(i,j,k)*q(i,j,k,5+jspc)*uu
          enddo
          !
        enddo
        !
        do n=1,numq
          dfcs(:,n)=ddfc(fcs(:,n),'222e',npdci,im)
        enddo
        !
        qrhs(is:ie,j,k,:)=qrhs(is:ie,j,k,:)+dfcs(is:ie,:)
        !
      enddo
      !
      deallocate(fcs,dfcs)
      !
      allocate(fcs(-hm:jm+hm,1:numq),dfcs(0:jm,1:numq))
      do i=is,ie
        !
        do j=-hm,jm+hm
          !
          uu=dxi(i,j,k,2,1)*vel(i,j,k,1)+dxi(i,j,k,2,2)*vel(i,j,k,2) + &
             dxi(i,j,k,2,3)*vel(i,j,k,3)
          fcs(j,1)=jacob(i,j,k)*  q(i,j,k,1)*uu
          fcs(j,2)=jacob(i,j,k)*( q(i,j,k,2)*uu+dxi(i,j,k,2,1)*prs(i,j,k) )
          fcs(j,3)=jacob(i,j,k)*( q(i,j,k,3)*uu+dxi(i,j,k,2,2)*prs(i,j,k) )
          fcs(j,4)=jacob(i,j,k)*( q(i,j,k,4)*uu+dxi(i,j,k,2,3)*prs(i,j,k) )
          fcs(j,5)=jacob(i,j,k)*( q(i,j,k,5)+prs(i,j,k) )*uu
          do jspc=1,num_species
            fcs(j,5+jspc)=jacob(i,j,k)*q(i,j,k,5+jspc)*uu
          enddo
          !
        enddo
        !
        do n=1,numq
          dfcs(:,n)=ddfc(fcs(:,n),'222e',npdcj,jm)
        enddo
        !
        qrhs(i,js:je,k,:)=qrhs(i,js:je,k,:)+dfcs(js:je,:)
        !
      enddo
      !
      deallocate(fcs,dfcs)
      !
    endif
    !
  end subroutine farfield_nscbc
  !+-------------------------------------------------------------------+
  !| The end of the subroutine farfield_nscbc.                         |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to apply outflow bc using nscbc.               |
  !+-------------------------------------------------------------------+
  !| ref: Jae Wook Kim, AIAA JOURNAL Vol. 38, No. 11, November 2000    |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 25-02-2021: Created by J. Fang @ Warrington                       |
  !+-------------------------------------------------------------------+
  subroutine outflow_nscbc(ndir)
    !
    use commvar,   only : xmin,xmax,mach
    use fludyna,   only : thermal,fvar2q,q2fvar,sos
    use commfunc,  only : deriv,ddfc,spafilter6exp
    use parallel,  only : psum,mpi_imax
    !
    ! arguments
    integer,intent(in) :: ndir
    !
    ! local data
    integer :: i,j,k,l,jspc,jq,ii,n,m
    real(8) :: pinv(5,5),pnor(5,5),Pmult(5,5),E(5),F(5),G(5),Rest(5),  &
               jcbi(3),LODi1(5),LODi(5)
    real(8),allocatable :: Ecs(:,:),dEcs(:),fcs(:,:),dfcs(:,:),qfilt(:,:)
    real(8) :: uu,css,gmachmax2,kinout,kin,var1,var2,psonic,rncout
    !
    gmachmax2=0.d0
    if(ndir==2 .and. irk==irkm) then
      i=im
      do k=0,km
      do j=0,jm
        css=sos(tmp(i,j,k)) 
        var1=1.d0/( dxi(i,j,k,1,1)**2+dxi(i,j,k,1,2)**2+               &
                    dxi(i,j,k,1,3)**2 )
        var2=vel(i,j,k,1)*dxi(i,j,k,1,1)+vel(i,j,k,2)*dxi(i,j,k,1,2)+  &
             vel(i,j,k,3)*dxi(i,j,k,1,3)
        gmachmax2=max(gmachmax2,var2*var2/(var1*css*css))
      enddo
      enddo
    endif
    gmachmax2=pmax(gmachmax2)
    !
    if(ndir==2 .and. irk==irkm) then
      !
      i=im
      !
      allocate(qfilt(0:jm,1:numq))
      do k=0,km
        !
        do jq=1,numq
          qfilt(:,jq)=spafilter6exp(q(i,:,k,jq),npdcj,jm)
          q(i,0:jm,k,jq)=qfilt(:,jq)
        enddo
        !
        ! open(18,file='profileq'//mpirankname//'.dat')
        ! write(18,"(3(1X,A15))")'y','q1','q1f'
        ! write(18,"(3(1X,E15.7E3))")(x(i,j,k,2),q(i,j,k,1),qfilt(j,1),j=0,jm)
        ! close(18)
        ! print*,' << profileq',mpirankname,'.dat'
        !
        call q2fvar(      q=  q(i,0:jm,k,:),                   &
                    density=rho(i,0:jm,k),                     &
                   velocity=vel(i,0:jm,k,:),                   &
                   pressure=prs(i,0:jm,k),                     &
                temperature=tmp(i,0:jm,k),                     &
                    species=spc(i,0:jm,k,:)                    )
      enddo
      deallocate(qfilt)
      !
      allocate(Ecs(0:2,1:numq),dEcs(1:numq))
      ! do k=ks,ke
      ! do j=js,je
      !
      ! psonic=0.d0
      ! rncout=0.d0
      ! do j=1,jm
      !   do k=0,0
      !     css=sos(tmp(i,j,k))
      !     if(vel(i,j,k,1)>css .and. vel(i,j-1,k,1)<=css) then
      !       psonic=psonic+prs(i,j,k)
      !       rncout=rncout+1.d0
      !     endif
      !   enddo
      ! enddo
      ! psonic=psum(psonic,comm=mpi_imax)/psum(rncout,comm=mpi_imax)
      ! print*,' ** psonic',psonic
      !
      do k=0,km
      do j=0,jm
        !
        pnor=pmatrix(rho(i,j,k),vel(i,j,k,1),vel(i,j,k,2),            &
                     vel(i,j,k,3),tmp(i,j,k),dxi(i,j,k,1,:),inv=.false.)
        pinv=pmatrix(rho(i,j,k),vel(i,j,k,1),vel(i,j,k,2),            &
                     vel(i,j,k,3),tmp(i,j,k),dxi(i,j,k,1,:),inv=.true.)
        !
        ! Pmult=MatMul(pnor,Pinv)
        ! if(jrk==0) then
        !   print*,'---------------------------------------------------------'
        !   write(*,"(5(1X,E15.7E3))")Pmult(1,:)
        !   write(*,"(5(1X,E15.7E3))")Pmult(2,:)
        !   write(*,"(5(1X,E15.7E3))")Pmult(3,:)
        !   write(*,"(5(1X,E15.7E3))")Pmult(4,:)
        !   write(*,"(5(1X,E15.7E3))")Pmult(5,:)
        ! end if
        !
        do ii=0,2
          uu=dxi(i-ii,j,k,1,1)*vel(i-ii,j,k,1) +                       &
             dxi(i-ii,j,k,1,2)*vel(i-ii,j,k,2) +                       &
             dxi(i-ii,j,k,1,3)*vel(i-ii,j,k,3)
          !
          Ecs(ii,1)=jacob(i-ii,j,k)*  q(i-ii,j,k,1)*uu
          Ecs(ii,2)=jacob(i-ii,j,k)*( q(i-ii,j,k,2)*uu+dxi(i-ii,j,k,1,1)*prs(i-ii,j,k) )
          Ecs(ii,3)=jacob(i-ii,j,k)*( q(i-ii,j,k,3)*uu+dxi(i-ii,j,k,1,2)*prs(i-ii,j,k) )
          Ecs(ii,4)=jacob(i-ii,j,k)*( q(i-ii,j,k,4)*uu+dxi(i-ii,j,k,1,3)*prs(i-ii,j,k) )
          Ecs(ii,5)=jacob(i-ii,j,k)*( q(i-ii,j,k,5)+prs(i-ii,j,k) )*uu
          do jspc=1,num_species
            Ecs(ii,5+jspc)=jacob(i-ii,j,k)*q(i-ii,j,k,5+jspc)*uu
          enddo
        enddo
        !
        do n=1,numq
          dEcs(n)=-deriv( Ecs(0,n),Ecs(1,n),Ecs(2,n) )
        enddo
        !
        E(1)= q(i,j,k,2)
        E(2)= q(i,j,k,2)*vel(i,j,k,1)+prs(i,j,k)
        E(3)= q(i,j,k,3)*vel(i,j,k,1)
        E(4)= q(i,j,k,4)*vel(i,j,k,1)
        E(5)=(q(i,j,k,5)+prs(i,j,k))*vel(i,j,k,1)
        !
        F(1)= q(i,j,k,3)
        F(2)= q(i,j,k,2)*vel(i,j,k,2)
        F(3)= q(i,j,k,3)*vel(i,j,k,2)+prs(i,j,k)
        F(4)= q(i,j,k,4)*vel(i,j,k,2)
        F(5)=(q(i,j,k,5)+prs(i,j,k))*vel(i,j,k,2)
        !
        G(1)= q(i,j,k,4)
        G(2)= q(i,j,k,2)*vel(i,j,k,3)
        G(3)= q(i,j,k,3)*vel(i,j,k,3)
        G(4)= q(i,j,k,4)*vel(i,j,k,3)+prs(i,j,k)
        G(5)=(q(i,j,k,5)+prs(i,j,k))*vel(i,j,k,3)
        !
        jcbi(1)=-deriv( dxi(i,j,k,1,1)  *jacob(i,j,k),                 &
                        dxi(i-1,j,k,1,1)*jacob(i-1,j,k),               &
                        dxi(i-2,j,k,1,1)*jacob(i-2,j,k) )
        jcbi(2)=-deriv( dxi(i,j,k,1,2)  *jacob(i,j,k),                 &
                        dxi(i-1,j,k,1,2)*jacob(i-1,j,k),               &
                        dxi(i-2,j,k,1,2)*jacob(i-2,j,k) )
        jcbi(3)=-deriv( dxi(i,j,k,1,3)  *jacob(i,j,k),                 &
                        dxi(i-1,j,k,1,3)*jacob(i-1,j,k),               &
                        dxi(i-2,j,k,1,3)*jacob(i-2,j,k) )
        !
        Rest(1)=E(1)*Jcbi(1)+F(1)*Jcbi(2)+G(1)*Jcbi(3)
        Rest(2)=E(2)*Jcbi(1)+F(2)*Jcbi(2)+G(2)*Jcbi(3)
        Rest(3)=E(3)*Jcbi(1)+F(3)*Jcbi(2)+G(3)*Jcbi(3)
        Rest(4)=E(4)*Jcbi(1)+F(4)*Jcbi(2)+G(4)*Jcbi(3)
        Rest(5)=E(5)*Jcbi(1)+F(5)*Jcbi(2)+G(5)*Jcbi(3)
        !
        LODi1(1)=dEcs(1)-Rest(1)
        LODi1(2)=dEcs(2)-Rest(2)
        LODi1(3)=dEcs(3)-Rest(3)
        LODi1(4)=dEcs(4)-Rest(4)
        LODi1(5)=dEcs(5)-Rest(5)
        !
        LODi=MatMul(pinv,LODi1)/jacob(i,j,k)
        ! do m=1,5
        !   LODi(m)=0.d0
        !   do n=1,5
        !     LODi(m)=pinv(m,n)*LODi1(n)
        !   enddo
        ! enddo
        ! LODi=LODi/jacob(i,j,k)
        !
        uu=dxi(i,j,k,1,1)*vel(i,j,k,1)+dxi(i,j,k,1,2)*vel(i,j,k,2) + &
           dxi(i,j,k,1,3)*vel(i,j,k,3)
        !
        css=sos(tmp(i,j,k))
        !
        ! if(uu>=0.d0) then
        !   !
        !   ! if(uu<css) then
        !   !   kinout=0.25d0*(1.d0-gmachmax2)*css/(xmax-xmin)
        !   !   LODi(5)=kinout*(prs(i,j,k)-pinf)/rho(i,j,k)/css
        !   ! endif
        !   ! else
        !   !   LODi(5)=0.d0
        !   ! endif
        !   !
        !   ! kinout=0.25d0*(1.d0-gmachmax2)*css/(ymax-ymin)
        !   ! LODi(5)=kinout*(pinf-prs(i,j,k))/rho(i,j,k)/css
        !   ! LODi(5)=kinout*(prs(i,j,k)-pinf)/rho(i,j,k)/css
        ! else
        !   ! back flow
        !   var1=1.d0/sqrt( dxi(i,j,k,1,1)**2+dxi(i,j,k,1,2)**2+         &
        !                   dxi(i,j,k,1,3)**2 )
        !   kin=-0.250*(1.d0-gmachmax2)*css/(xmax-xmin)
        !   !
        !   LODi(1)=0.d0
        !   !
        !   LODi(2)=kin*0.5d0*( dxi(i,j,k,1,1)*var1*(vel(i,j,k,3)-winf)- &
        !                       dxi(i,j,k,1,3)*var1*(vel(i,j,k,1)-uinf) )
        !   LODi(3)=kin*0.5d0*(-dxi(i,j,k,2,1)*var1*(vel(i,j,k,2)-vinf)+ &
        !                       dxi(i,j,k,2,2)*var1*(vel(i,j,k,1)-uinf) )
        !   LODi(5)=kin*(dxi(i,j,k,1,1)*var1*(vel(i,j,k,1)-uinf)+        &
        !                dxi(i,j,k,1,2)*var1*(vel(i,j,k,2)-vinf)+        &
        !                dxi(i,j,k,1,3)*var1*(vel(i,j,k,3)-winf)+        &
        !                (prs(i,j,k)-pinf)/rho(i,j,k)/css )
        ! endif  
        !
        LODi1=MatMul(pnor,LODi)*jacob(i,j,k)
        !
        dEcs(1)=LODi1(1)+Rest(1)
        dEcs(2)=LODi1(2)+Rest(2)
        dEcs(3)=LODi1(3)+Rest(3)
        dEcs(4)=LODi1(4)+Rest(4)
        dEcs(5)=LODi1(5)+Rest(5)
        !
        qrhs(i,j,k,:)=qrhs(i,j,k,:)+dEcs(:)
        !
      enddo
      enddo
      !
      deallocate(Ecs,dEcs)
      !
      allocate(fcs(-hm:jm+hm,1:numq),dfcs(0:jm,1:numq))
      do k=ks,ke
        !
        do j=-hm,jm+hm
          !
          uu=dxi(i,j,k,2,1)*vel(i,j,k,1)+dxi(i,j,k,2,2)*vel(i,j,k,2) + &
             dxi(i,j,k,2,3)*vel(i,j,k,3)
          fcs(j,1)=jacob(i,j,k)*  q(i,j,k,1)*uu
          fcs(j,2)=jacob(i,j,k)*( q(i,j,k,2)*uu+dxi(i,j,k,2,1)*prs(i,j,k) )
          fcs(j,3)=jacob(i,j,k)*( q(i,j,k,3)*uu+dxi(i,j,k,2,2)*prs(i,j,k) )
          fcs(j,4)=jacob(i,j,k)*( q(i,j,k,4)*uu+dxi(i,j,k,2,3)*prs(i,j,k) )
          fcs(j,5)=jacob(i,j,k)*( q(i,j,k,5)+prs(i,j,k) )*uu
          do jspc=1,num_species
            fcs(j,5+jspc)=jacob(i,j,k)*q(i,j,k,5+jspc)*uu
          enddo
          !
        enddo
        !
        do n=1,numq
          dfcs(:,n)=ddfc(fcs(:,n),'222e',npdcj,jm)
        enddo
        !
        qrhs(i,js:je,k,:)=qrhs(i,js:je,k,:)+dfcs(js:je,:)
        !
      enddo
      !
      deallocate(fcs,dfcs)
      !
      if(ndims==3) then
        !
        allocate(fcs(-hm:km+hm,1:numq),dfcs(0:km,1:numq))
        do j=js,je
          !
          do k=-hm,km+hm
            !
            uu=dxi(i,j,k,3,1)*vel(i,j,k,1)+dxi(i,j,k,3,2)*vel(i,j,k,2) + &
               dxi(i,j,k,3,3)*vel(i,j,k,3)
            fcs(k,1)=jacob(i,j,k)*  q(i,j,k,1)*uu
            fcs(k,2)=jacob(i,j,k)*( q(i,j,k,2)*uu+dxi(i,j,k,3,1)*prs(i,j,k) )
            fcs(k,3)=jacob(i,j,k)*( q(i,j,k,3)*uu+dxi(i,j,k,3,2)*prs(i,j,k) )
            fcs(k,4)=jacob(i,j,k)*( q(i,j,k,4)*uu+dxi(i,j,k,3,3)*prs(i,j,k) )
            fcs(k,5)=jacob(i,j,k)*( q(i,j,k,5)+prs(i,j,k) )*uu
            do jspc=1,num_species
              fcs(k,5+jspc)=jacob(i,j,k)*q(i,j,k,5+jspc)*uu
            enddo
            !
          enddo
          !
          do n=1,numq
            dfcs(:,n)=ddfc(fcs(:,n),'222e',npdck,km)
          enddo
          !
          qrhs(i,j,ks:ke,:)=qrhs(i,j,ks:ke,:)+dfcs(ks:ke,:)
          !
        enddo
        !
        deallocate(fcs,dfcs)
        !
      endif
      !
    endif
    !
  end subroutine outflow_nscbc
  !
  function pmatrix(rho,u,v,w,t,ddi,inv)
    !
    use commvar, only : gamma
    use fludyna, only : sos
    !
    ! arguments
    real(8) :: pmatrix(5,5)
    logical,intent(in) :: inv
    real(8),intent(in) :: rho,u,v,w,t,ddi(3)
    !
    ! local data
    real(8) :: var1,nx,ny,nz,a1,a2,a3,a4,a5,a6,c,VV,phi
    !
    c=sos(t)
    var1=1.d0/sqrt(ddi(1)*ddi(1)+ddi(2)*ddi(2)+ddi(3)*ddi(3))
    nx=ddi(1)*var1
    ny=ddi(2)*var1
    nz=ddi(3)*var1
    phi=0.5d0*(gamma-1.d0)*(u*u+v*v+w*w)
    VV=nx*u+ny*v+nz*w
    a1=gamma-1.d0
    a2=1.d0/(sqrt(2.d0)*rho*c)
    a3=rho/(c*sqrt(2.d0))
    a4=(phi+c*c)/(gamma-1.d0)
    a5=1.d0-phi/(c*c)
    a6=phi/(gamma-1.d0)
    !
    ! pmatrix=0.d0
    !
    if(inv) then
      !
      pmatrix(1,1)= nx*a5-(nz*v-ny*w)/rho
      pmatrix(1,2)= nx*a1*u/(c*c)
      pmatrix(1,3)= nx*a1*v/(c*c)+nz/rho
      pmatrix(1,4)= nx*a1*w/(c*c)-ny/rho
      pmatrix(1,5)=-nx*a1/(c*c)
      ! 
      pmatrix(2,1)= ny*a5-(nx*w-nz*u)/rho
      pmatrix(2,2)= ny*a1*u/(c*c)-nz/rho
      pmatrix(2,3)= ny*a1*v/(c*c)
      pmatrix(2,4)= ny*a1*w/(c*c)+nx/rho
      pmatrix(2,5)=-ny*a1/(c*c)
      ! 
      pmatrix(3,1)= nz*a5-(ny*u-nx*v)/rho
      pmatrix(3,2)= nz*a1*u/(c*c)+ny/rho
      pmatrix(3,3)= nz*a1*v/(c*c)-nx/rho
      pmatrix(3,4)= nz*a1*w/(c*c)
      pmatrix(3,5)=-nz*a1/(c*c)
      ! 
      pmatrix(4,1)= a2*(phi-c*VV)
      pmatrix(4,2)=-a2*(a1*u-nx*c)
      pmatrix(4,3)=-a2*(a1*v-ny*c)
      pmatrix(4,4)=-a2*(a1*w-nz*c)
      pmatrix(4,5)= a1*a2
      ! 
      pmatrix(5,1)= a2*(phi+c*VV)
      pmatrix(5,2)=-a2*(a1*u+nx*c)
      pmatrix(5,3)=-a2*(a1*v+ny*c)
      pmatrix(5,4)=-a2*(a1*w+nz*c)
      pmatrix(5,5)= a1*a2
      ! 
    else
      !
      pmatrix(1,1)=nx
      pmatrix(1,2)=ny
      pmatrix(1,3)=nz
      pmatrix(1,4)=a3
      pmatrix(1,5)=a3
      !
      pmatrix(2,1)=nx*u
      pmatrix(2,2)=ny*u-nz*rho
      pmatrix(2,3)=nz*u+ny*rho
      pmatrix(2,4)=a3*(u+nx*c)
      pmatrix(2,5)=a3*(u-nx*c)
      ! 
      pmatrix(3,1)=nx*v+nz*rho
      pmatrix(3,2)=ny*v
      pmatrix(3,3)=nz*v-nx*rho
      pmatrix(3,4)=a3*(v+ny*c)
      pmatrix(3,5)=a3*(v-ny*c)
      !
      pmatrix(4,1)=nx*w-ny*rho
      pmatrix(4,2)=ny*w+nx*rho
      pmatrix(4,3)=nz*w
      pmatrix(4,4)=a3*(w+nz*c)
      pmatrix(4,5)=a3*(w-nz*c)
      !
      pmatrix(5,1)=nx*a6+rho*(nz*v-ny*w)
      pmatrix(5,2)=ny*a6+rho*(nx*w-nz*u)
      pmatrix(5,3)=nz*a6+rho*(ny*u-nx*v)
      pmatrix(5,4)=a3*(a4+c*VV)
      pmatrix(5,5)=a3*(a4-c*VV)
      !
    endif
    !
  end function pmatrix
  !
  ! function pmatrix2(rho,u,v,w,t,ddi,inv) result(pmatrix)
  !   !
  !   use commvar, only : gamma
  !   use fludyna, only : sos
  !   !
  !   ! arguments
  !   real(8) :: pmatrix(5,5)
  !   logical,intent(in) :: inv
  !   real(8),intent(in) :: rho,u,v,w,t,ddi(3)
  !   !
  !   ! local data
  !   real(8) :: ke,ma2,css,cs2,b0(3),lxi(3),vlxi(3),cplus(3),cminu(3),  &
  !              b(3)
  !   real(8) :: var1,gamm1,rhi,h,cvlxi
  !   !
  !   gamm1=gamma-1.d0
  !   rhi  =1.d0/rho
  !   !
  !   ke=(u**2+v**2+w**2)
  !   css=sos(t)
  !   cs2=css*css
  !   ma2=ke/cs2
  !   !
  !   var1=1.d0/sqrt(ddi(1)*ddi(1)+ddi(2)*ddi(2)+ddi(3)*ddi(3))
  !   lxi(1)=ddi(1)*var1
  !   lxi(2)=ddi(2)*var1
  !   lxi(3)=ddi(3)*var1
  !   !
  !   vlxi(1)=lxi(3)*v-lxi(2)*w
  !   vlxi(2)=lxi(1)*w-lxi(3)*u
  !   vlxi(3)=lxi(2)*u-lxi(1)*v
  !   !
  !   pmatrix=0.d0
  !   !
  !   if(inv) then
  !     !
  !     b0(1)=(1.d0-0.5d0*gamm1*ma2)*lxi(1)-vlxi(1)*rhi
  !     b0(2)=(1.d0-0.5d0*gamm1*ma2)*lxi(2)-vlxi(2)*rhi
  !     b0(3)=(1.d0-0.5d0*gamm1*ma2)*lxi(3)-vlxi(3)*rhi
  !     !
  !     cplus(1)=(lxi(1)-gamm1/css*u)*rhi
  !     cplus(2)=(lxi(2)-gamm1/css*v)*rhi
  !     cplus(3)=(lxi(3)-gamm1/css*w)*rhi
  !     !
  !     cminu(1)=(-lxi(1)-gamm1/css*u)*rhi
  !     cminu(2)=(-lxi(2)-gamm1/css*v)*rhi
  !     cminu(3)=(-lxi(3)-gamm1/css*w)*rhi
  !     !
  !     pmatrix(1,1)= b0(1)
  !     pmatrix(1,2)= gamm1*u/cs2*lxi(1)
  !     pmatrix(1,3)= gamm1*v/cs2*lxi(1)+lxi(3)*rhi
  !     pmatrix(1,4)= gamm1*w/cs2*lxi(1)-lxi(2)*rhi
  !     pmatrix(1,5)=-gamm1/cs2*lxi(1)
  !     !
  !     pmatrix(2,1)= b0(2)
  !     pmatrix(2,2)= gamm1*u/cs2*lxi(2)-lxi(3)*rhi
  !     pmatrix(2,3)= gamm1*v/cs2*lxi(2)
  !     pmatrix(2,4)= gamm1*w/cs2*lxi(2)+lxi(1)*rhi
  !     pmatrix(2,5)=-gamm1/cs2*lxi(2)
  !     !
  !     pmatrix(3,1)= b0(3)
  !     pmatrix(3,2)= gamm1*u/cs2*lxi(3)+lxi(2)*rhi
  !     pmatrix(3,3)= gamm1*v/cs2*lxi(3)-lxi(1)*rhi
  !     pmatrix(3,4)= gamm1*w/cs2*lxi(3)
  !     pmatrix(3,5)=-gamm1/cs2*lxi(3)
  !     !
  !     pmatrix(4,1)= css*rhi*(0.5d0*gamm1*ma2-(u*lxi(1)+v*lxi(2)+w*lxi(3))/css)
  !     pmatrix(4,2)= cplus(1)
  !     pmatrix(4,3)= cplus(2)
  !     pmatrix(4,4)= cplus(3)
  !     pmatrix(4,5)= gamm1/css*rhi
  !     !
  !     pmatrix(5,1)= css*rhi*(0.5d0*gamm1*ma2+(u*lxi(1)+v*lxi(2)+w*lxi(3))/css)
  !     pmatrix(5,2)= cminu(1)
  !     pmatrix(5,3)= cminu(2)
  !     pmatrix(5,4)= cminu(3)
  !     pmatrix(5,5)= gamm1/css*rhi
  !     !
  !   else
  !     !
  !     pmatrix(1,1)=lxi(1)
  !     pmatrix(1,2)=lxi(2)
  !     pmatrix(1,3)=lxi(3)
  !     pmatrix(1,4)=0.5d0*rho/css
  !     pmatrix(1,5)=0.5d0*rho/css
  !     !
  !     pmatrix(2,1)=lxi(1)*u
  !     pmatrix(2,2)=lxi(2)*u-lxi(3)*rho
  !     pmatrix(2,3)=lxi(3)*u+lxi(2)*rho
  !     pmatrix(2,4)=0.5d0*rho/css*(u+lxi(1)*css)
  !     pmatrix(2,5)=0.5d0*rho/css*(u-lxi(1)*css)
  !     !
  !     pmatrix(3,1)=lxi(1)*v+lxi(3)*rho
  !     pmatrix(3,2)=lxi(2)*v
  !     pmatrix(3,3)=lxi(3)*v-lxi(1)*rho
  !     pmatrix(3,4)=0.5d0*rho/css*(v+lxi(2)*css)
  !     pmatrix(3,5)=0.5d0*rho/css*(v-lxi(2)*css)
  !     !
  !     pmatrix(4,1)=lxi(1)*w-lxi(2)*rho
  !     pmatrix(4,2)=lxi(2)*w+lxi(1)*rho
  !     pmatrix(4,3)=lxi(3)*w
  !     pmatrix(4,4)=0.5d0*rho/css*(w+lxi(3)*css)
  !     pmatrix(4,5)=0.5d0*rho/css*(w-lxi(3)*css)
  !     !
  !     b(1)=0.5d0*ke*lxi(1)+rho*vlxi(1)
  !     b(2)=0.5d0*ke*lxi(2)+rho*vlxi(2)
  !     b(3)=0.5d0*ke*lxi(3)+rho*vlxi(3)
  !     !
  !     h=0.5d0*ke+css*css/gamm1
  !     cvlxi=css*(lxi(1)*u+lxi(2)*v+lxi(3)*w)
  !     !
  !     pmatrix(5,1)=b(1)
  !     pmatrix(5,2)=b(2)
  !     pmatrix(5,3)=b(3)
  !     pmatrix(5,4)=0.5d0*rho/css*(h+cvlxi)
  !     pmatrix(5,5)=0.5d0*rho/css*(h-cvlxi)
  !     !
  !   endif
  !   !
  ! end function pmatrix2
  !
  ! subroutine outflow_nscbc(ndir)
  !   !
  !   use commarray, only : prs,vel,tmp,rho,spc,q,qrhs,dxi,jacob
  !   use fludyna,   only : thermal,fvar2q,q2fvar,sos
  !   use commfunc,  only : deriv
  !   !
  !   ! arguments
  !   integer,intent(in) :: ndir
  !   !
  !   ! local data
  !   integer :: i,j,k,l,jspec
  !   real(8) :: css,ri,uu,vv,ww,qq,gm
  !   real(8),allocatable :: qinf(:),vel_inf(:),spc_inf(:)
  !   real(8) :: jac(5,5),jacinv(5,5),el(5,5),er(5,5),ev(5),             &
  !              dwdxi(5),dwdxo(5),dwcdxi(5),dwcdxo(5),dwcdx(5),df(5)
  !   !
  !   if(ndir==2 .and. irk==irkm) then
  !     !
  !     i=im
  !     !
  !     allocate(qinf(numq),vel_inf(3),spc_inf(num_species))
  !     !
  !     vel_inf(1)=uinf
  !     vel_inf(2)=vinf
  !     vel_inf(3)=winf
  !     spc_inf(1)=0.d0
  !     !
  !     call fvar2q(      q=  qinf, density=roinf, velocity=vel_inf,     &
  !                                 pressure=pinf,  species=spc_inf      )
  !     !
  !     do k=0,km
  !     do j=0,jm
  !       !
  !       css=sos(tmp(i,j,k))
  !       !
  !       ri  =  1.d0/rho(i,j,k)
  !       uu  =  vel(i,j,k,1)
  !       vv  =  vel(i,j,k,2)
  !       ww  =  vel(i,j,k,3)
  !       qq  =  0.5d0 * (uu*uu  + vv*vv + ww*ww)
  !       !
  !       !   Jacobian of conservative/primitive transformation          
  !       !   (Eqn. (A.5) of Lodato et al, JCP 2008)
  !       jac(1,1)  =  1.d0
  !       jac(1,2)  =  0.d0
  !       jac(1,3)  =  0.d0
  !       jac(1,4)  =  0.d0
  !       jac(1,5)  =  0.d0
  !       jac(2,1)  =  uu
  !       jac(2,2)  =  rho(i,j,k)
  !       jac(2,3)  =  0.d0
  !       jac(2,4)  =  0.d0
  !       jac(2,5)  =  0.d0
  !       jac(3,1)  =  vv
  !       jac(3,2)  =  0.d0
  !       jac(3,3)  =  rho(i,j,k)
  !       jac(3,4)  =  0.d0
  !       jac(3,5)  =  0.d0
  !       jac(4,1)  =  ww
  !       jac(4,2)  =  0.d0
  !       jac(4,3)  =  0.d0
  !       jac(4,4)  =  rho(i,j,k)
  !       jac(4,5)  =  0.d0
  !       jac(5,1)  =  qq
  !       jac(5,2)  =  q(i,j,k,2)
  !       jac(5,3)  =  q(i,j,k,3)
  !       jac(5,4)  =  q(i,j,k,4)
  !       jac(5,5)  =  1.d0/(gamma-1.d0)
  !       !
  !       ! Jacobian of inverse conservative/primitive transformation 
  !       ! (Eqn. (A.5) of Lodato et al, JCP 2008)
  !       jacinv(1,1) =  1.d0
  !       jacinv(1,2) =  0.d0
  !       jacinv(1,3) =  0.d0
  !       jacinv(1,4) =  0.d0
  !       jacinv(1,5) =  0.d0
  !       jacinv(2,1) = -uu*ri
  !       jacinv(2,2) =  ri
  !       jacinv(2,3) =  0.d0
  !       jacinv(2,4) =  0.d0
  !       jacinv(2,5) =  0.d0
  !       jacinv(3,1) = -vv*ri
  !       jacinv(3,2) =  0.d0
  !       jacinv(3,3) =  ri
  !       jacinv(3,4) =  0.d0
  !       jacinv(3,5) =  0.d0
  !       jacinv(4,1) = -ww*ri
  !       jacinv(4,2) =  0.d0
  !       jacinv(4,3) =  0.d0
  !       jacinv(4,4) =  ri
  !       jacinv(4,5) =  0.d0
  !       jacinv(5,1) =  (gamma-1.d0)*qq
  !       jacinv(5,2) = -(gamma-1.d0)*uu
  !       jacinv(5,3) = -(gamma-1.d0)*vv
  !       jacinv(5,4) = -(gamma-1.d0)*ww
  !       jacinv(5,5) =  (gamma-1.d0)
  !       !
  !       ! left eigenvectors matrix (Eqn.d0 (A.12d0) of 
  !       ! Lodato et al, JCP 2008)
  !       el(1,1) =  0.d0
  !       el(1,2) = -rho(i,j,k)*css
  !       el(1,3) =  0.d0
  !       el(1,4) =  0.d0
  !       el(1,5) =  1.d0
  !       el(2,1) =  css*css
  !       el(2,2) =  0.d0
  !       el(2,3) =  0.d0
  !       el(2,4) =  0.d0
  !       el(2,5) = -1.d0
  !       el(3,1) =  0.d0
  !       el(3,2) =  0.d0
  !       el(3,3) =  1.d0
  !       el(3,4) =  0.d0
  !       el(3,5) =  0.d0
  !       el(4,1) =  0.d0
  !       el(4,2) =  0.d0
  !       el(4,3) =  0.d0
  !       el(4,4) =  1.d0
  !       el(4,5) =  0.d0
  !       el(5,1) =  0.d0
  !       el(5,2) =  rho(i,j,k)*css
  !       el(5,3) =  0.d0
  !       el(5,4) =  0.d0
  !       el(5,5) =  1.d0
  !       !
  !       ! left eigenvectors matrix (Eqn.d0 (A.11d0) of Lodato et al, JCP 2008)
  !       er(1,1) =  0.5d0/css/css
  !       er(2,1) = -0.5d0*ri/css
  !       er(3,1) =  0.d0
  !       er(4,1) =  0.d0
  !       er(5,1) =  0.5d0
  !       er(1,2) =  1.d0/css/css
  !       er(2,2) =  0.d0
  !       er(3,2) =  0.d0
  !       er(4,2) =  0.d0
  !       er(5,2) =  0.d0
  !       er(1,3) =  0.d0
  !       er(2,3) =  0.d0
  !       er(3,3) =  1.d0
  !       er(4,3) =  0.d0
  !       er(5,3) =  0.d0
  !       er(1,4) =  0.d0
  !       er(2,4) =  0.d0
  !       er(3,4) =  0.d0
  !       er(4,4) =  1.d0
  !       er(5,4) =  0.d0
  !       er(1,5) =  0.5d0/css/css
  !       er(2,5) =  0.5d0*ri/css
  !       er(3,5) =  0.d0
  !       er(4,5) =  0.d0
  !       er(5,5) =  0.5d0
  !       !
  !       !  Eigenvalues 
  !       ev(1)    =  uu-css
  !       ev(2)    =  uu
  !       ev(3)    =  uu
  !       ev(4)    =  uu
  !       ev(5)    =  uu+css
  !       !
  !       ! Derivatives of conservative variables
  !       ! Inner derivatives
  !       dwdxi(1) = -deriv(q(i,j,k,1),q(i-1,j,k,1),q(i-2,j,k,1))
  !       dwdxi(2) = -deriv(q(i,j,k,2),q(i-1,j,k,2),q(i-2,j,k,2))
  !       dwdxi(3) = -deriv(q(i,j,k,3),q(i-1,j,k,3),q(i-2,j,k,3))
  !       dwdxi(4) = -deriv(q(i,j,k,4),q(i-1,j,k,4),q(i-2,j,k,4))
  !       dwdxi(5) = -deriv(q(i,j,k,5),q(i-1,j,k,5),q(i-2,j,k,5))
  !       !
  !       ! Outer derivatives
  !       dwdxo(1) = -deriv(qinf(1),q(i,j,k,1))
  !       dwdxo(2) = -deriv(qinf(2),q(i,j,k,2))
  !       dwdxo(3) = -deriv(qinf(3),q(i,j,k,3))
  !       dwdxo(4) = -deriv(qinf(4),q(i,j,k,4))
  !       dwdxo(5) = -deriv(qinf(5),q(i,j,k,5))
  !       !
  !       !   Derivatives of characteristic variables
  !       dwcdxi = 0.d0
  !       dwcdxo = 0.d0
  !       !m=1
  !       !mm=1
  !       dwcdxi(1) = dwcdxi(1) +el(1,1)*jacinv(1,1)*dwdxi(1) !  inner
  !       dwcdxi(1) = dwcdxi(1) +el(1,1)*jacinv(1,2)*dwdxi(2) !  inner
  !       dwcdxi(1) = dwcdxi(1) +el(1,1)*jacinv(1,3)*dwdxi(3) !  inner
  !       dwcdxi(1) = dwcdxi(1) +el(1,1)*jacinv(1,4)*dwdxi(4) !  inner
  !       dwcdxi(1) = dwcdxi(1) +el(1,1)*jacinv(1,5)*dwdxi(5) !  inner
  !       !mm=2
  !       dwcdxi(1) =dwcdxi(1) +el(1,2)*jacinv(2,1)*dwdxi(1) !  inner
  !       dwcdxi(1) =dwcdxi(1) +el(1,2)*jacinv(2,2)*dwdxi(2) !  inner
  !       dwcdxi(1) =dwcdxi(1) +el(1,2)*jacinv(2,3)*dwdxi(3) !  inner
  !       dwcdxi(1) =dwcdxi(1) +el(1,2)*jacinv(2,4)*dwdxi(4) !  inner
  !       dwcdxi(1) =dwcdxi(1) +el(1,2)*jacinv(2,5)*dwdxi(5) !  inner
  !       !mm=3
  !       dwcdxi(1) =dwcdxi(1) +el(1,3)*jacinv(3,1)*dwdxi(1) !  inner
  !       dwcdxi(1) =dwcdxi(1) +el(1,3)*jacinv(3,2)*dwdxi(2) !  inner
  !       dwcdxi(1) =dwcdxi(1) +el(1,3)*jacinv(3,3)*dwdxi(3) !  inner
  !       dwcdxi(1) =dwcdxi(1) +el(1,3)*jacinv(3,4)*dwdxi(4) !  inner
  !       dwcdxi(1) =dwcdxi(1) +el(1,3)*jacinv(3,5)*dwdxi(5) !  inner
  !       !mm=4
  !       dwcdxi(1) =dwcdxi(1) +el(1,4)*jacinv(4,1)*dwdxi(1) !  inner
  !       dwcdxi(1) =dwcdxi(1) +el(1,4)*jacinv(4,2)*dwdxi(2) !  inner
  !       dwcdxi(1) =dwcdxi(1) +el(1,4)*jacinv(4,3)*dwdxi(3) !  inner
  !       dwcdxi(1) =dwcdxi(1) +el(1,4)*jacinv(4,4)*dwdxi(4) !  inner
  !       dwcdxi(1) =dwcdxi(1) +el(1,4)*jacinv(4,5)*dwdxi(5) !  inner
  !       !mm=5
  !       dwcdxi(1) =dwcdxi(1) +el(1,5)*jacinv(5,1)*dwdxi(1) !  inner
  !       dwcdxi(1) =dwcdxi(1) +el(1,5)*jacinv(5,2)*dwdxi(2) !  inner
  !       dwcdxi(1) =dwcdxi(1) +el(1,5)*jacinv(5,3)*dwdxi(3) !  inner
  !       dwcdxi(1) =dwcdxi(1) +el(1,5)*jacinv(5,4)*dwdxi(4) !  inner
  !       dwcdxi(1) =dwcdxi(1) +el(1,5)*jacinv(5,5)*dwdxi(5) !  inner
  !       !
  !       dwcdxo(1) =dwcdxo(1) +el(1,5)*jacinv(5,1)*dwdxo(1) !  outer
  !       dwcdxo(1) =dwcdxo(1) +el(1,5)*jacinv(5,2)*dwdxo(2) !  outer
  !       dwcdxo(1) =dwcdxo(1) +el(1,5)*jacinv(5,3)*dwdxo(3) !  outer
  !       dwcdxo(1) =dwcdxo(1) +el(1,5)*jacinv(5,4)*dwdxo(4) !  outer
  !       dwcdxo(1) =dwcdxo(1) +el(1,5)*jacinv(5,5)*dwdxo(5) !  outer
  !       !
  !       !m=2
  !       !
  !       !mm=1
  !       dwcdxi(2) =dwcdxi(2) +el(2,1)*jacinv(1,1)*dwdxi(1) !  inner
  !       dwcdxi(2) =dwcdxi(2) +el(2,1)*jacinv(1,2)*dwdxi(2) !  inner
  !       dwcdxi(2) =dwcdxi(2) +el(2,1)*jacinv(1,3)*dwdxi(3) !  inner
  !       dwcdxi(2) =dwcdxi(2) +el(2,1)*jacinv(1,4)*dwdxi(4) !  inner
  !       dwcdxi(2) =dwcdxi(2) +el(2,1)*jacinv(1,5)*dwdxi(5) !  inner
  !       !
  !       dwcdxi(2) =dwcdxi(2) +el(2,2)*jacinv(2,1)*dwdxi(1) !  inner
  !       dwcdxi(2) =dwcdxi(2) +el(2,2)*jacinv(2,2)*dwdxi(2) !  inner
  !       dwcdxi(2) =dwcdxi(2) +el(2,2)*jacinv(2,3)*dwdxi(3) !  inner
  !       dwcdxi(2) =dwcdxi(2) +el(2,2)*jacinv(2,4)*dwdxi(4) !  inner
  !       dwcdxi(2) =dwcdxi(2) +el(2,2)*jacinv(2,5)*dwdxi(5) !  inner
  !       !
  !       dwcdxi(2) =dwcdxi(2) +el(2,3)*jacinv(3,1)*dwdxi(1) !  inner
  !       dwcdxi(2) =dwcdxi(2) +el(2,3)*jacinv(3,2)*dwdxi(2) !  inner
  !       dwcdxi(2) =dwcdxi(2) +el(2,3)*jacinv(3,3)*dwdxi(3) !  inner
  !       dwcdxi(2) =dwcdxi(2) +el(2,3)*jacinv(3,4)*dwdxi(4) !  inner
  !       dwcdxi(2) =dwcdxi(2) +el(2,3)*jacinv(3,5)*dwdxi(5) !  inner
  !       !
  !       dwcdxi(2) =dwcdxi(2) +el(2,4)*jacinv(4,1)*dwdxi(1) !  inner
  !       dwcdxi(2) =dwcdxi(2) +el(2,4)*jacinv(4,2)*dwdxi(2) !  inner
  !       dwcdxi(2) =dwcdxi(2) +el(2,4)*jacinv(4,3)*dwdxi(3) !  inner
  !       dwcdxi(2) =dwcdxi(2) +el(2,4)*jacinv(4,4)*dwdxi(4) !  inner
  !       dwcdxi(2) =dwcdxi(2) +el(2,4)*jacinv(4,5)*dwdxi(5) !  inner
  !       !
  !       dwcdxi(2) =dwcdxi(2) +el(2,5)*jacinv(5,1)*dwdxi(1) !  inner
  !       dwcdxi(2) =dwcdxi(2) +el(2,5)*jacinv(5,2)*dwdxi(2) !  inner
  !       dwcdxi(2) =dwcdxi(2) +el(2,5)*jacinv(5,3)*dwdxi(3) !  inner
  !       dwcdxi(2) =dwcdxi(2) +el(2,5)*jacinv(5,4)*dwdxi(4) !  inner
  !       dwcdxi(2) =dwcdxi(2) +el(2,5)*jacinv(5,5)*dwdxi(5) !  inner
  !       !
  !       dwcdxo(2) =dwcdxo(2) +el(2,5)*jacinv(5,1)*dwdxo(1) !  outer
  !       dwcdxo(2) =dwcdxo(2) +el(2,5)*jacinv(5,2)*dwdxo(2) !  outer
  !       dwcdxo(2) =dwcdxo(2) +el(2,5)*jacinv(5,3)*dwdxo(3) !  outer
  !       dwcdxo(2) =dwcdxo(2) +el(2,5)*jacinv(5,4)*dwdxo(4) !  outer
  !       dwcdxo(2) =dwcdxo(2) +el(2,5)*jacinv(5,5)*dwdxo(5) !  outer
  !       !
  !       !m=3
  !       !
  !       !mm=1
  !       dwcdxi(3) =dwcdxi(3) +el(3,1)*jacinv(1,1)*dwdxi(1) !  inner
  !       dwcdxi(3) =dwcdxi(3) +el(3,1)*jacinv(1,2)*dwdxi(2) !  inner
  !       dwcdxi(3) =dwcdxi(3) +el(3,1)*jacinv(1,3)*dwdxi(3) !  inner
  !       dwcdxi(3) =dwcdxi(3) +el(3,1)*jacinv(1,4)*dwdxi(4) !  inner
  !       dwcdxi(3) =dwcdxi(3) +el(3,1)*jacinv(1,5)*dwdxi(5) !  inner
  !       !
  !       dwcdxi(3) =dwcdxi(3) +el(3,2)*jacinv(2,1)*dwdxi(1) !  inner
  !       dwcdxi(3) =dwcdxi(3) +el(3,2)*jacinv(2,2)*dwdxi(2) !  inner
  !       dwcdxi(3) =dwcdxi(3) +el(3,2)*jacinv(2,3)*dwdxi(3) !  inner
  !       dwcdxi(3) =dwcdxi(3) +el(3,2)*jacinv(2,4)*dwdxi(4) !  inner
  !       dwcdxi(3) =dwcdxi(3) +el(3,2)*jacinv(2,5)*dwdxi(5) !  inner
  !       !
  !       dwcdxi(3) =dwcdxi(3) +el(3,3)*jacinv(3,1)*dwdxi(1) !  inner
  !       dwcdxi(3) =dwcdxi(3) +el(3,3)*jacinv(3,2)*dwdxi(2) !  inner
  !       dwcdxi(3) =dwcdxi(3) +el(3,3)*jacinv(3,3)*dwdxi(3) !  inner
  !       dwcdxi(3) =dwcdxi(3) +el(3,3)*jacinv(3,4)*dwdxi(4) !  inner
  !       dwcdxi(3) =dwcdxi(3) +el(3,3)*jacinv(3,5)*dwdxi(5) !  inner
  !       !
  !       dwcdxi(3) =dwcdxi(3) +el(3,4)*jacinv(4,1)*dwdxi(1) !  inner
  !       dwcdxi(3) =dwcdxi(3) +el(3,4)*jacinv(4,2)*dwdxi(2) !  inner
  !       dwcdxi(3) =dwcdxi(3) +el(3,4)*jacinv(4,3)*dwdxi(3) !  inner
  !       dwcdxi(3) =dwcdxi(3) +el(3,4)*jacinv(4,4)*dwdxi(4) !  inner
  !       dwcdxi(3) =dwcdxi(3) +el(3,4)*jacinv(4,5)*dwdxi(5) !  inner
  !       !
  !       dwcdxi(3) =dwcdxi(3) +el(3,5)*jacinv(5,1)*dwdxi(1) !  inner
  !       dwcdxi(3) =dwcdxi(3) +el(3,5)*jacinv(5,2)*dwdxi(2) !  inner
  !       dwcdxi(3) =dwcdxi(3) +el(3,5)*jacinv(5,3)*dwdxi(3) !  inner
  !       dwcdxi(3) =dwcdxi(3) +el(3,5)*jacinv(5,4)*dwdxi(4) !  inner
  !       dwcdxi(3) =dwcdxi(3) +el(3,5)*jacinv(5,5)*dwdxi(5) !  inner
  !       !
  !       dwcdxo(3) =dwcdxo(3) +el(3,5)*jacinv(5,1)*dwdxo(1) !  outer
  !       dwcdxo(3) =dwcdxo(3) +el(3,5)*jacinv(5,2)*dwdxo(2) !  outer
  !       dwcdxo(3) =dwcdxo(3) +el(3,5)*jacinv(5,3)*dwdxo(3) !  outer
  !       dwcdxo(3) =dwcdxo(3) +el(3,5)*jacinv(5,4)*dwdxo(4) !  outer
  !       dwcdxo(3) =dwcdxo(3) +el(3,5)*jacinv(5,5)*dwdxo(5) !  outer
  !       !
  !       !m=4
  !       !
  !       !mm=1
  !       dwcdxi(4) =dwcdxi(4) +el(4,1)*jacinv(1,1)*dwdxi(1) !  inner
  !       dwcdxi(4) =dwcdxi(4) +el(4,1)*jacinv(1,2)*dwdxi(2) !  inner
  !       dwcdxi(4) =dwcdxi(4) +el(4,1)*jacinv(1,3)*dwdxi(3) !  inner
  !       dwcdxi(4) =dwcdxi(4) +el(4,1)*jacinv(1,4)*dwdxi(4) !  inner
  !       dwcdxi(4) =dwcdxi(4) +el(4,1)*jacinv(1,5)*dwdxi(5) !  inner
  !       !
  !       dwcdxi(4) =dwcdxi(4) +el(4,2)*jacinv(2,1)*dwdxi(1) !  inner
  !       dwcdxi(4) =dwcdxi(4) +el(4,2)*jacinv(2,2)*dwdxi(2) !  inner
  !       dwcdxi(4) =dwcdxi(4) +el(4,2)*jacinv(2,3)*dwdxi(3) !  inner
  !       dwcdxi(4) =dwcdxi(4) +el(4,2)*jacinv(2,4)*dwdxi(4) !  inner
  !       dwcdxi(4) =dwcdxi(4) +el(4,2)*jacinv(2,5)*dwdxi(5) !  inner
  !       !
  !       dwcdxi(4) =dwcdxi(4) +el(4,3)*jacinv(3,1)*dwdxi(1) !  inner
  !       dwcdxi(4) =dwcdxi(4) +el(4,3)*jacinv(3,2)*dwdxi(2) !  inner
  !       dwcdxi(4) =dwcdxi(4) +el(4,3)*jacinv(3,3)*dwdxi(3) !  inner
  !       dwcdxi(4) =dwcdxi(4) +el(4,3)*jacinv(3,4)*dwdxi(4) !  inner
  !       dwcdxi(4) =dwcdxi(4) +el(4,3)*jacinv(3,5)*dwdxi(5) !  inner
  !       !
  !       dwcdxi(4) =dwcdxi(4) +el(4,4)*jacinv(4,1)*dwdxi(1) !  inner
  !       dwcdxi(4) =dwcdxi(4) +el(4,4)*jacinv(4,2)*dwdxi(2) !  inner
  !       dwcdxi(4) =dwcdxi(4) +el(4,4)*jacinv(4,3)*dwdxi(3) !  inner
  !       dwcdxi(4) =dwcdxi(4) +el(4,4)*jacinv(4,4)*dwdxi(4) !  inner
  !       dwcdxi(4) =dwcdxi(4) +el(4,4)*jacinv(4,5)*dwdxi(5) !  inner
  !       !
  !       dwcdxi(4) =dwcdxi(4) +el(4,5)*jacinv(5,1)*dwdxi(1) !  inner
  !       dwcdxi(4) =dwcdxi(4) +el(4,5)*jacinv(5,2)*dwdxi(2) !  inner
  !       dwcdxi(4) =dwcdxi(4) +el(4,5)*jacinv(5,3)*dwdxi(3) !  inner
  !       dwcdxi(4) =dwcdxi(4) +el(4,5)*jacinv(5,4)*dwdxi(4) !  inner
  !       dwcdxi(4) =dwcdxi(4) +el(4,5)*jacinv(5,5)*dwdxi(5) !  inner
  !       !
  !       dwcdxo(4) =dwcdxo(4) +el(4,5)*jacinv(5,1)*dwdxo(1) !  outer
  !       dwcdxo(4) =dwcdxo(4) +el(4,5)*jacinv(5,2)*dwdxo(2) !  outer
  !       dwcdxo(4) =dwcdxo(4) +el(4,5)*jacinv(5,3)*dwdxo(3) !  outer
  !       dwcdxo(4) =dwcdxo(4) +el(4,5)*jacinv(5,4)*dwdxo(4) !  outer
  !       dwcdxo(4) =dwcdxo(4) +el(4,5)*jacinv(5,5)*dwdxo(5) !  outer
  !       !
  !       !m=5
  !       !
  !       !mm=1
  !       dwcdxi(5) =dwcdxi(5) +el(5,1)*jacinv(1,1)*dwdxi(1) !  inner
  !       dwcdxi(5) =dwcdxi(5) +el(5,1)*jacinv(1,2)*dwdxi(2) !  inner
  !       dwcdxi(5) =dwcdxi(5) +el(5,1)*jacinv(1,3)*dwdxi(3) !  inner
  !       dwcdxi(5) =dwcdxi(5) +el(5,1)*jacinv(1,4)*dwdxi(4) !  inner
  !       dwcdxi(5) =dwcdxi(5) +el(5,1)*jacinv(1,5)*dwdxi(5) !  inner
  !       !
  !       dwcdxi(5) =dwcdxi(5) +el(5,2)*jacinv(2,1)*dwdxi(1) !  inner
  !       dwcdxi(5) =dwcdxi(5) +el(5,2)*jacinv(2,2)*dwdxi(2) !  inner
  !       dwcdxi(5) =dwcdxi(5) +el(5,2)*jacinv(2,3)*dwdxi(3) !  inner
  !       dwcdxi(5) =dwcdxi(5) +el(5,2)*jacinv(2,4)*dwdxi(4) !  inner
  !       dwcdxi(5) =dwcdxi(5) +el(5,2)*jacinv(2,5)*dwdxi(5) !  inner
  !       !
  !       dwcdxi(5) =dwcdxi(5) +el(5,3)*jacinv(3,1)*dwdxi(1) !  inner
  !       dwcdxi(5) =dwcdxi(5) +el(5,3)*jacinv(3,2)*dwdxi(2) !  inner
  !       dwcdxi(5) =dwcdxi(5) +el(5,3)*jacinv(3,3)*dwdxi(3) !  inner
  !       dwcdxi(5) =dwcdxi(5) +el(5,3)*jacinv(3,4)*dwdxi(4) !  inner
  !       dwcdxi(5) =dwcdxi(5) +el(5,3)*jacinv(3,5)*dwdxi(5) !  inner
  !       !
  !       dwcdxi(5) =dwcdxi(5) +el(5,4)*jacinv(4,1)*dwdxi(1) !  inner
  !       dwcdxi(5) =dwcdxi(5) +el(5,4)*jacinv(4,2)*dwdxi(2) !  inner
  !       dwcdxi(5) =dwcdxi(5) +el(5,4)*jacinv(4,3)*dwdxi(3) !  inner
  !       dwcdxi(5) =dwcdxi(5) +el(5,4)*jacinv(4,4)*dwdxi(4) !  inner
  !       dwcdxi(5) =dwcdxi(5) +el(5,4)*jacinv(4,5)*dwdxi(5) !  inner
  !       !
  !       dwcdxi(5) =dwcdxi(5) +el(5,5)*jacinv(5,1)*dwdxi(1) !  inner
  !       dwcdxi(5) =dwcdxi(5) +el(5,5)*jacinv(5,2)*dwdxi(2) !  inner
  !       dwcdxi(5) =dwcdxi(5) +el(5,5)*jacinv(5,3)*dwdxi(3) !  inner
  !       dwcdxi(5) =dwcdxi(5) +el(5,5)*jacinv(5,4)*dwdxi(4) !  inner
  !       dwcdxi(5) =dwcdxi(5) +el(5,5)*jacinv(5,5)*dwdxi(5) !  inner
  !       !
  !       dwcdxo(5) =dwcdxo(5) +el(5,5)*jacinv(5,1)*dwdxo(1) !  outer
  !       dwcdxo(5) =dwcdxo(5) +el(5,5)*jacinv(5,2)*dwdxo(2) !  outer
  !       dwcdxo(5) =dwcdxo(5) +el(5,5)*jacinv(5,3)*dwdxo(3) !  outer
  !       dwcdxo(5) =dwcdxo(5) +el(5,5)*jacinv(5,4)*dwdxo(4) !  outer
  !       dwcdxo(5) =dwcdxo(5) +el(5,5)*jacinv(5,5)*dwdxo(5) !  outer
  !       !
  !       !   Enforce LODI relations
  !       !m=1
  !       if (ev(1)>0.d0) then
  !         ! Waves pointing out of the domain
  !         dwcdx(1) = dwcdxi(1)
  !       else
  !         ! Waves entering the domain
  !         dwcdx(1) = dwcdxo(1) ! n.r with or without relaxation
  !       endif
  !       !m=2
  !       if (ev(2)>0.d0) then
  !         ! Waves pointing out of the domain
  !         dwcdx(2) = dwcdxi(2)
  !       else
  !         ! Waves entering the domain
  !         dwcdx(2) = dwcdxo(2) ! n.r. with or without relaxation
  !       endif
  !       !m=3
  !       if (ev(3)>0.d0) then
  !         ! Waves pointing out of the domain
  !         dwcdx(3) = dwcdxi(3)
  !       else
  !         ! Waves entering the domain
  !         dwcdx(3) = dwcdxo(3) ! n.r with or without relaxation
  !       endif
  !       !m=4
  !       if (ev(4)>0.d0) then
  !         ! Waves pointing out of the domain
  !         dwcdx(4) = dwcdxi(4)
  !       else
  !         ! Waves entering the domain
  !         dwcdx(4) = dwcdxo(4) ! n.r with or without relaxation
  !       endif
  !       !m=5
  !       if (ev(5)>0.d0) then
  !         ! Waves pointing out of the domain
  !         dwcdx(5) = dwcdxi(5)
  !       else
  !         ! Waves entering the domain
  !         dwcdx(5) = dwcdxo(5) ! n.r with or without relaxation
  !       endif
  !       !
  !       !   Amplitude of characteristic waves
  !       dwcdx = dwcdx * ev
  !       !
  !       !  Return to conservative variables 
  !       df=0.d0
  !       !
  !       !mm=1
  !       df(1) = df(1) + jac(1,1)*er(1,1)*dwcdx(1)
  !       df(1) = df(1) + jac(1,1)*er(1,2)*dwcdx(2)
  !       df(1) = df(1) + jac(1,1)*er(1,3)*dwcdx(3)
  !       df(1) = df(1) + jac(1,1)*er(1,4)*dwcdx(4)
  !       df(1) = df(1) + jac(1,1)*er(1,5)*dwcdx(5)
  !       df(1) = df(1) + jac(1,2)*er(2,1)*dwcdx(1)
  !       df(1) = df(1) + jac(1,2)*er(2,2)*dwcdx(2)
  !       df(1) = df(1) + jac(1,2)*er(2,3)*dwcdx(3)
  !       df(1) = df(1) + jac(1,2)*er(2,4)*dwcdx(4)
  !       df(1) = df(1) + jac(1,2)*er(2,5)*dwcdx(5)
  !       df(1) = df(1) + jac(1,3)*er(3,1)*dwcdx(1)
  !       df(1) = df(1) + jac(1,3)*er(3,2)*dwcdx(2)
  !       df(1) = df(1) + jac(1,3)*er(3,3)*dwcdx(3)
  !       df(1) = df(1) + jac(1,3)*er(3,4)*dwcdx(4)
  !       df(1) = df(1) + jac(1,3)*er(3,5)*dwcdx(5)
  !       df(1) = df(1) + jac(1,4)*er(4,1)*dwcdx(1)
  !       df(1) = df(1) + jac(1,4)*er(4,2)*dwcdx(2)
  !       df(1) = df(1) + jac(1,4)*er(4,3)*dwcdx(3)
  !       df(1) = df(1) + jac(1,4)*er(4,4)*dwcdx(4)
  !       df(1) = df(1) + jac(1,4)*er(4,5)*dwcdx(5)
  !       df(1) = df(1) + jac(1,5)*er(5,1)*dwcdx(1)
  !       df(1) = df(1) + jac(1,5)*er(5,2)*dwcdx(2)
  !       df(1) = df(1) + jac(1,5)*er(5,3)*dwcdx(3)
  !       df(1) = df(1) + jac(1,5)*er(5,4)*dwcdx(4)
  !       df(1) = df(1) + jac(1,5)*er(5,5)*dwcdx(5)
  !       !
  !       df(2) = df(2) + jac(2,1)*er(1,1)*dwcdx(1)
  !       df(2) = df(2) + jac(2,1)*er(1,2)*dwcdx(2)
  !       df(2) = df(2) + jac(2,1)*er(1,3)*dwcdx(3)
  !       df(2) = df(2) + jac(2,1)*er(1,4)*dwcdx(4)
  !       df(2) = df(2) + jac(2,1)*er(1,5)*dwcdx(5)
  !       df(2) = df(2) + jac(2,2)*er(2,1)*dwcdx(1)
  !       df(2) = df(2) + jac(2,2)*er(2,2)*dwcdx(2)
  !       df(2) = df(2) + jac(2,2)*er(2,3)*dwcdx(3)
  !       df(2) = df(2) + jac(2,2)*er(2,4)*dwcdx(4)
  !       df(2) = df(2) + jac(2,2)*er(2,5)*dwcdx(5)
  !       df(2) = df(2) + jac(2,3)*er(3,1)*dwcdx(1)
  !       df(2) = df(2) + jac(2,3)*er(3,2)*dwcdx(2)
  !       df(2) = df(2) + jac(2,3)*er(3,3)*dwcdx(3)
  !       df(2) = df(2) + jac(2,3)*er(3,4)*dwcdx(4)
  !       df(2) = df(2) + jac(2,3)*er(3,5)*dwcdx(5)
  !       df(2) = df(2) + jac(2,4)*er(4,1)*dwcdx(1)
  !       df(2) = df(2) + jac(2,4)*er(4,2)*dwcdx(2)
  !       df(2) = df(2) + jac(2,4)*er(4,3)*dwcdx(3)
  !       df(2) = df(2) + jac(2,4)*er(4,4)*dwcdx(4)
  !       df(2) = df(2) + jac(2,4)*er(4,5)*dwcdx(5)
  !       df(2) = df(2) + jac(2,5)*er(5,1)*dwcdx(1)
  !       df(2) = df(2) + jac(2,5)*er(5,2)*dwcdx(2)
  !       df(2) = df(2) + jac(2,5)*er(5,3)*dwcdx(3)
  !       df(2) = df(2) + jac(2,5)*er(5,4)*dwcdx(4)
  !       df(2) = df(2) + jac(2,5)*er(5,5)*dwcdx(5)

  !       df(3) = df(3) + jac(3,1)*er(1,1)*dwcdx(1)
  !       df(3) = df(3) + jac(3,1)*er(1,2)*dwcdx(2)
  !       df(3) = df(3) + jac(3,1)*er(1,3)*dwcdx(3)
  !       df(3) = df(3) + jac(3,1)*er(1,4)*dwcdx(4)
  !       df(3) = df(3) + jac(3,1)*er(1,5)*dwcdx(5)
  !       df(3) = df(3) + jac(3,2)*er(2,1)*dwcdx(1)
  !       df(3) = df(3) + jac(3,2)*er(2,2)*dwcdx(2)
  !       df(3) = df(3) + jac(3,2)*er(2,3)*dwcdx(3)
  !       df(3) = df(3) + jac(3,2)*er(2,4)*dwcdx(4)
  !       df(3) = df(3) + jac(3,2)*er(2,5)*dwcdx(5)
  !       df(3) = df(3) + jac(3,3)*er(3,1)*dwcdx(1)
  !       df(3) = df(3) + jac(3,3)*er(3,2)*dwcdx(2)
  !       df(3) = df(3) + jac(3,3)*er(3,3)*dwcdx(3)
  !       df(3) = df(3) + jac(3,3)*er(3,4)*dwcdx(4)
  !       df(3) = df(3) + jac(3,3)*er(3,5)*dwcdx(5)
  !       df(3) = df(3) + jac(3,4)*er(4,1)*dwcdx(1)
  !       df(3) = df(3) + jac(3,4)*er(4,2)*dwcdx(2)
  !       df(3) = df(3) + jac(3,4)*er(4,3)*dwcdx(3)
  !       df(3) = df(3) + jac(3,4)*er(4,4)*dwcdx(4)
  !       df(3) = df(3) + jac(3,4)*er(4,5)*dwcdx(5)
  !       df(3) = df(3) + jac(3,5)*er(5,1)*dwcdx(1)
  !       df(3) = df(3) + jac(3,5)*er(5,2)*dwcdx(2)
  !       df(3) = df(3) + jac(3,5)*er(5,3)*dwcdx(3)
  !       df(3) = df(3) + jac(3,5)*er(5,4)*dwcdx(4)
  !       df(3) = df(3) + jac(3,5)*er(5,5)*dwcdx(5)
  !       !
  !       df(4) = df(4) + jac(4,1)*er(1,1)*dwcdx(1)
  !       df(4) = df(4) + jac(4,1)*er(1,2)*dwcdx(2)
  !       df(4) = df(4) + jac(4,1)*er(1,3)*dwcdx(3)
  !       df(4) = df(4) + jac(4,1)*er(1,4)*dwcdx(4)
  !       df(4) = df(4) + jac(4,1)*er(1,5)*dwcdx(5)
  !       df(4) = df(4) + jac(4,2)*er(2,1)*dwcdx(1)
  !       df(4) = df(4) + jac(4,2)*er(2,2)*dwcdx(2)
  !       df(4) = df(4) + jac(4,2)*er(2,3)*dwcdx(3)
  !       df(4) = df(4) + jac(4,2)*er(2,4)*dwcdx(4)
  !       df(4) = df(4) + jac(4,2)*er(2,5)*dwcdx(5)
  !       df(4) = df(4) + jac(4,3)*er(3,1)*dwcdx(1)
  !       df(4) = df(4) + jac(4,3)*er(3,2)*dwcdx(2)
  !       df(4) = df(4) + jac(4,3)*er(3,3)*dwcdx(3)
  !       df(4) = df(4) + jac(4,3)*er(3,4)*dwcdx(4)
  !       df(4) = df(4) + jac(4,3)*er(3,5)*dwcdx(5)
  !       df(4) = df(4) + jac(4,4)*er(4,1)*dwcdx(1)
  !       df(4) = df(4) + jac(4,4)*er(4,2)*dwcdx(2)
  !       df(4) = df(4) + jac(4,4)*er(4,3)*dwcdx(3)
  !       df(4) = df(4) + jac(4,4)*er(4,4)*dwcdx(4)
  !       df(4) = df(4) + jac(4,4)*er(4,5)*dwcdx(5)
  !       df(4) = df(4) + jac(4,5)*er(5,1)*dwcdx(1)
  !       df(4) = df(4) + jac(4,5)*er(5,2)*dwcdx(2)
  !       df(4) = df(4) + jac(4,5)*er(5,3)*dwcdx(3)
  !       df(4) = df(4) + jac(4,5)*er(5,4)*dwcdx(4)
  !       df(4) = df(4) + jac(4,5)*er(5,5)*dwcdx(5)
  !       !
  !       df(5) = df(5) + jac(5,1)*er(1,1)*dwcdx(1)
  !       df(5) = df(5) + jac(5,1)*er(1,2)*dwcdx(2)
  !       df(5) = df(5) + jac(5,1)*er(1,3)*dwcdx(3)
  !       df(5) = df(5) + jac(5,1)*er(1,4)*dwcdx(4)
  !       df(5) = df(5) + jac(5,1)*er(1,5)*dwcdx(5)
  !       df(5) = df(5) + jac(5,2)*er(2,1)*dwcdx(1)
  !       df(5) = df(5) + jac(5,2)*er(2,2)*dwcdx(2)
  !       df(5) = df(5) + jac(5,2)*er(2,3)*dwcdx(3)
  !       df(5) = df(5) + jac(5,2)*er(2,4)*dwcdx(4)
  !       df(5) = df(5) + jac(5,2)*er(2,5)*dwcdx(5)
  !       df(5) = df(5) + jac(5,3)*er(3,1)*dwcdx(1)
  !       df(5) = df(5) + jac(5,3)*er(3,2)*dwcdx(2)
  !       df(5) = df(5) + jac(5,3)*er(3,3)*dwcdx(3)
  !       df(5) = df(5) + jac(5,3)*er(3,4)*dwcdx(4)
  !       df(5) = df(5) + jac(5,3)*er(3,5)*dwcdx(5)
  !       df(5) = df(5) + jac(5,4)*er(4,1)*dwcdx(1)
  !       df(5) = df(5) + jac(5,4)*er(4,2)*dwcdx(2)
  !       df(5) = df(5) + jac(5,4)*er(4,3)*dwcdx(3)
  !       df(5) = df(5) + jac(5,4)*er(4,4)*dwcdx(4)
  !       df(5) = df(5) + jac(5,5)*er(5,1)*dwcdx(1)
  !       df(5) = df(5) + jac(5,5)*er(5,2)*dwcdx(2)
  !       df(5) = df(5) + jac(5,5)*er(5,3)*dwcdx(3)
  !       df(5) = df(5) + jac(5,5)*er(5,4)*dwcdx(4)
  !       df(5) = df(5) + jac(5,5)*er(5,5)*dwcdx(5)
  !       !
  !       qrhs(i,j,k,1) = df(1)*dxi(i,j,k,1,1)/jacob(i,j,k)
  !       qrhs(i,j,k,2) = df(2)*dxi(i,j,k,1,1)/jacob(i,j,k)
  !       qrhs(i,j,k,3) = df(3)*dxi(i,j,k,1,1)/jacob(i,j,k)
  !       qrhs(i,j,k,4) = df(4)*dxi(i,j,k,1,1)/jacob(i,j,k)
  !       qrhs(i,j,k,5) = df(5)*dxi(i,j,k,1,1)/jacob(i,j,k)
  !       !
  !     enddo
  !     enddo
  !     !
  !   endif
  !   !
  ! end subroutine outflow_nscbc
  !+-------------------------------------------------------------------+
  !| The end of the subroutine outflow_nscbc.                          |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to apply nonslip bc.                           |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 13-02-2021: Created by J. Fang @ Warrington                       |
  !+-------------------------------------------------------------------+
  subroutine noslip(ndir,tw)
    !
    use commvar,   only : Reynolds,turbmode,lreport
    use commarray, only : tke,omg
    use fludyna,   only : thermal,fvar2q,q2fvar,miucal
    use commfunc,  only : dis2point2
    use parallel,  only : mpi_jmin,bcast
    !
    ! arguments
    integer,intent(in) :: ndir
    real(8),intent(in) :: tw
    !
    ! local data
    integer :: i,j,k,l,jspec,fh
    real(8) :: pe,delta_d1,miu,beta1
    real(8) :: vwall(0:im,0:km)
    character(len=64) :: filewbs
    logical :: lexist
    !
    real(8),save :: beter,wallamplit,xa,xb,xc
    integer,save :: nmod_t,nmod_z
    logical,save :: lfirstcal=.true.
    !
    beta1=0.075d0
    !
    if(ndir==3) then
      !
      if(jrk==0) then
        !
        if(lreport) lfirstcal=.true.
        !
        if(lfirstcal) then
          !
          if(mpirank==0) then
            !
            filewbs='datin/wallbs.dat'
            !
            inquire(file=trim(filewbs), exist=lexist)
            !
            if(lexist) then
              fh=get_unit()
              open(fh,file=trim(filewbs),action='read')
              read(fh,*)
              read(fh,*)
              read(fh,*)wallamplit,beter,xa,xb,xc
              read(fh,*)
              read(fh,*)nmod_t,nmod_z
              close(fh)
              print*,' >> ',trim(filewbs)
              !
              print*,' ------------- wall blowing & suction parameters -------------'
              write(*,"(38x,A,1X,F12.5)")   '  amplitude:',wallamplit
              write(*,"(38x,A,1X,F12.5)")   '  frequency:',beter
              write(*,"(23x,3(A,1X,F12.5))")'   x extent:',xa,' ~',xb,' ~',xc
              write(*,"(42x,(A,1X,I0))")'    temporal modes:',nmod_t
              write(*,"(42x,(A,1X,I0))")'    spanwise modes:',nmod_z
              print*,' --------------------------------------------------------------'
            else
              wallamplit=0.d0
              beter=1.d0
              xa=0.d0
              xb=0.d0
              xc=0.d0
              nmod_t=0
              nmod_z=0
            endif
            !
          endif
          !
          call bcast(wallamplit,comm=mpi_jmin)
          call bcast(beter,comm=mpi_jmin)
          call bcast(xa,comm=mpi_jmin)
          call bcast(xb,comm=mpi_jmin)
          call bcast(xc,comm=mpi_jmin)
          call bcast(nmod_t,comm=mpi_jmin)
          call bcast(nmod_z,comm=mpi_jmin)
          !
          !
          lfirstcal=.false.
          !
        endif
        !
        if(ndims==3 .and. wallamplit>1.d-10) then
          vwall=wallbs_rand(beter,wallamplit,xa,xb,xc,nmod_t,nmod_z)
        else
          vwall=0.d0
        endif
        !
        j=0
        !
        do k=0,km
        do i=0,im
          !
          pe=num1d3*(4.d0*prs(i,1,k)-prs(i,2,k))
          !
          vel(i,j,k,1)=0.d0
          vel(i,j,k,2)=vwall(i,k)
          vel(i,j,k,3)=0.d0
          prs(i,j,k)  =pe
          tmp(i,j,k)  =tw
          !
          do jspec=1,num_species
            spc(i,j,k,jspec)=num1d3*(4.d0*spc(i,1,k,jspec)-spc(i,2,k,jspec))
          enddo
          !
          rho(i,j,k)  =thermal(pressure=prs(i,j,k),temperature=tmp(i,j,k))
          !
          call fvar2q(      q=  q(i,j,k,:),                            &
                      density=rho(i,j,k),                              &
                     velocity=vel(i,j,k,:),                            &
                     pressure=prs(i,j,k),                              &
                      species=spc(i,j,k,:)                             )
          !
          if(trim(turbmode)=='k-omega') then
            tke(i,j,k)=0.d0
            !
            delta_d1=dis2point2(x(i,j,k,:),x(i,j+1,k,:))
            miu=miucal(tmp(i,j,k))/Reynolds
            omg(i,j,k)=50.d0*miu/rho(i,j,k)/beta1/delta_d1
            !
            q(i,j,k,5+num_species+1)=tke(i,j,k)*rho(i,j,k)
            q(i,j,k,5+num_species+2)=omg(i,j,k)*rho(i,j,k)
            !
          endif
          !
        enddo
        enddo
        !
      endif
      !
    elseif(ndir==4) then
      !
      if(jrk==jrkm) then
        !
        j=jm
        do k=0,km
        do i=0,im
          pe=num1d3*(4.d0*prs(i,j-1,k)-prs(i,j-2,k))
          !
          vel(i,j,k,1)=0.d0
          vel(i,j,k,2)=0.d0
          vel(i,j,k,3)=0.d0
          prs(i,j,k)  =pe
          tmp(i,j,k)  =tw
          !
          if(num_species>0) then
            spc(i,j,k,1)=num1d3*(4.d0*spc(i,j-1,k,1)-spc(i,j-2,k,1))
          endif
          !
          rho(i,j,k)  =thermal(pressure=prs(i,j,k),temperature=tmp(i,j,k))
          !
          call fvar2q(      q=  q(i,j,k,:),                            &
                      density=rho(i,j,k),                              &
                     velocity=vel(i,j,k,:),                            &
                     pressure=prs(i,j,k),                              &
                      species=spc(i,j,k,:)                             )
          !
          if(trim(turbmode)=='k-omega') then
            tke(i,j,k)=0.d0
            !
            delta_d1=dis2point2(x(i,j,k,:),x(i,j-1,k,:))
            miu=miucal(tmp(i,j,k))/Reynolds
            omg(i,j,k)=50.d0*miu/rho(i,j,k)/beta1/delta_d1
            !
            q(i,j,k,5+num_species+1)=tke(i,j,k)*rho(i,j,k)
            q(i,j,k,5+num_species+2)=omg(i,j,k)*rho(i,j,k)
            !
          endif
          !
        enddo
        enddo
        !
      endif
      !
    else
      stop ' !! ndir not defined @ noslip'
    endif
    !
  end subroutine noslip
  !!
  subroutine noslip_adibatic(ndir)
    !
    use commvar,   only : Reynolds,turbmode,lreport
    use commarray, only : tke,omg
    use fludyna,   only : thermal,fvar2q,q2fvar,miucal
    use commfunc,  only : dis2point2
    use parallel,  only : mpi_jmin,bcast
    !
    ! arguments
    integer,intent(in) :: ndir
    !
    ! local data
    integer :: i,j,k,l,jspec,fh
    real(8) :: pe,te,delta_d1,miu,beta1
    real(8) :: vwall(0:im,0:km)
    character(len=64) :: filewbs
    logical :: lexist
    !
    real(8),save :: beter,wallamplit,xa,xb,xc
    integer,save :: nmod_t,nmod_z
    logical,save :: lfirstcal=.true.
    !
    beta1=0.075d0
    !
    if(ndir==3) then
      !
      if(jrk==0) then
        !
        if(lreport) lfirstcal=.true.
        !
        if(lfirstcal) then
          !
          if(mpirank==0) then
            !
            filewbs='datin/wallbs.dat'
            !
            inquire(file=trim(filewbs), exist=lexist)
            !
            if(lexist) then
              fh=get_unit()
              open(fh,file=trim(filewbs),action='read')
              read(fh,*)
              read(fh,*)
              read(fh,*)wallamplit,beter,xa,xb,xc
              read(fh,*)
              read(fh,*)nmod_t,nmod_z
              close(fh)
              print*,' >> ',trim(filewbs)
              !
              print*,' ------------- wall blowing & suction parameters -------------'
              write(*,"(38x,A,1X,F12.5)")   '  amplitude:',wallamplit
              write(*,"(38x,A,1X,F12.5)")   '  frequency:',beter
              write(*,"(23x,3(A,1X,F12.5))")'   x extent:',xa,' ~',xb,' ~',xc
              write(*,"(42x,(A,1X,I0))")'    temporal modes:',nmod_t
              write(*,"(42x,(A,1X,I0))")'    spanwise modes:',nmod_z
              print*,' --------------------------------------------------------------'
            else
              wallamplit=0.d0
              beter=1.d0
              xa=0.d0
              xb=0.d0
              xc=0.d0
              nmod_t=0
              nmod_z=0
            endif
            !
          endif
          !
          call bcast(wallamplit,comm=mpi_jmin)
          call bcast(beter,comm=mpi_jmin)
          call bcast(xa,comm=mpi_jmin)
          call bcast(xb,comm=mpi_jmin)
          call bcast(xc,comm=mpi_jmin)
          call bcast(nmod_t,comm=mpi_jmin)
          call bcast(nmod_z,comm=mpi_jmin)
          !
          !
          lfirstcal=.false.
          !
        endif
        !
        if(ndims==3 .and. wallamplit>1.d-10) then
          vwall=wallbs_rand(beter,wallamplit,xa,xb,xc,nmod_t,nmod_z)
        else
          vwall=0.d0
        endif
        !
        j=0
        !
        do k=0,km
        do i=0,im
          !
          pe=num1d3*(4.d0*prs(i,1,k)-prs(i,2,k))
          te=num1d3*(4.d0*tmp(i,1,k)-tmp(i,2,k))
          !
          vel(i,j,k,1)=0.d0
          vel(i,j,k,2)=vwall(i,k)
          vel(i,j,k,3)=0.d0
          prs(i,j,k)  =pe
          tmp(i,j,k)  =te
          !
          do jspec=1,num_species
            spc(i,j,k,jspec)=num1d3*(4.d0*spc(i,1,k,jspec)-spc(i,2,k,jspec))
          enddo
          !
          rho(i,j,k)  =thermal(pressure=prs(i,j,k),temperature=tmp(i,j,k))
          !
          call fvar2q(      q=  q(i,j,k,:),                            &
                      density=rho(i,j,k),                              &
                     velocity=vel(i,j,k,:),                            &
                     pressure=prs(i,j,k),                              &
                      species=spc(i,j,k,:)                             )
          !
          if(trim(turbmode)=='k-omega') then
            tke(i,j,k)=0.d0
            !
            delta_d1=dis2point2(x(i,j,k,:),x(i,j+1,k,:))
            miu=miucal(tmp(i,j,k))/Reynolds
            omg(i,j,k)=50.d0*miu/rho(i,j,k)/beta1/delta_d1
            !
            q(i,j,k,5+num_species+1)=tke(i,j,k)*rho(i,j,k)
            q(i,j,k,5+num_species+2)=omg(i,j,k)*rho(i,j,k)
            !
          endif
          !
        enddo
        enddo
        !
      endif
      !
    elseif(ndir==4) then
      !
      if(jrk==jrkm) then
        !
        j=jm
        do k=0,km
        do i=0,im
          pe=num1d3*(4.d0*prs(i,j-1,k)-prs(i,j-2,k))
          te=num1d3*(4.d0*tmp(i,j-1,k)-tmp(i,j-2,k))
          !
          vel(i,j,k,1)=0.d0
          vel(i,j,k,2)=0.d0
          vel(i,j,k,3)=0.d0
          prs(i,j,k)  =pe
          tmp(i,j,k)  =te
          !
          if(num_species>0) then
            spc(i,j,k,1)=num1d3*(4.d0*spc(i,j-1,k,1)-spc(i,j-2,k,1))
          endif
          !
          rho(i,j,k)  =thermal(pressure=prs(i,j,k),temperature=tmp(i,j,k))
          !
          call fvar2q(      q=  q(i,j,k,:),                            &
                      density=rho(i,j,k),                              &
                     velocity=vel(i,j,k,:),                            &
                     pressure=prs(i,j,k),                              &
                      species=spc(i,j,k,:)                             )
          !
          if(trim(turbmode)=='k-omega') then
            tke(i,j,k)=0.d0
            !
            delta_d1=dis2point2(x(i,j,k,:),x(i,j-1,k,:))
            miu=miucal(tmp(i,j,k))/Reynolds
            omg(i,j,k)=50.d0*miu/rho(i,j,k)/beta1/delta_d1
            !
            q(i,j,k,5+num_species+1)=tke(i,j,k)*rho(i,j,k)
            q(i,j,k,5+num_species+2)=omg(i,j,k)*rho(i,j,k)
            !
          endif
          !
        enddo
        enddo
        !
      endif
      !
    else
      stop ' !! ndir not defined @ noslip'
    endif
    !
  end subroutine noslip_adibatic
  !
  subroutine slipisotwall(ndir,tw)
    !
    use commvar,   only : Reynolds,turbmode,lreport
    use commarray, only : tke,omg
    use fludyna,   only : thermal,fvar2q,q2fvar,miucal
    use commfunc,  only : dis2point2
    use parallel,  only : mpi_jmin,bcast
    !
    ! arguments
    integer,intent(in) :: ndir
    real(8),intent(in) :: tw
    !
    ! local data
    integer :: i,j,k,l,jspec,fh
    real(8) :: pe,ue,te,delta_d1,miu,beta1
    real(8) :: vwall(0:im,0:km)
    character(len=64) :: filewbs
    logical :: lexist
    !
    real(8),save :: beter,wallamplit,xa,xb,xc
    integer,save :: nmod_t,nmod_z
    logical,save :: lfirstcal=.true.
    !
    beta1=0.075d0
    !
    if(ndir==3) then
      !
      if(jrk==0) then
        !
        if(lreport) lfirstcal=.true.
        !
        if(lfirstcal) then
          !
          if(mpirank==0) then
            !
            filewbs='datin/wallbs.dat'
            !
            inquire(file=trim(filewbs), exist=lexist)
            !
            if(lexist) then
              fh=get_unit()
              open(fh,file=trim(filewbs),action='read')
              read(fh,*)
              read(fh,*)
              read(fh,*)wallamplit,beter,xa,xb,xc
              read(fh,*)
              read(fh,*)nmod_t,nmod_z
              close(fh)
              print*,' >> ',trim(filewbs)
              !
              print*,' ------------- wall blowing & suction parameters -------------'
              write(*,"(38x,A,1X,F12.5)")   '  amplitude:',wallamplit
              write(*,"(38x,A,1X,F12.5)")   '  frequency:',beter
              write(*,"(23x,3(A,1X,F12.5))")'   x extent:',xa,' ~',xb,' ~',xc
              write(*,"(42x,(A,1X,I0))")'    temporal modes:',nmod_t
              write(*,"(42x,(A,1X,I0))")'    spanwise modes:',nmod_z
              print*,' --------------------------------------------------------------'
            else
              wallamplit=0.d0
              beter=1.d0
              xa=0.d0
              xb=0.d0
              xc=0.d0
              nmod_t=0
              nmod_z=0
            endif
          endif
          !
          call bcast(wallamplit,comm=mpi_jmin)
          call bcast(beter,comm=mpi_jmin)
          call bcast(xa,comm=mpi_jmin)
          call bcast(xb,comm=mpi_jmin)
          call bcast(xc,comm=mpi_jmin)
          call bcast(nmod_t,comm=mpi_jmin)
          call bcast(nmod_z,comm=mpi_jmin)
          !
          !
          lfirstcal=.false.
          !
        endif
        !
        if(ndims==3 .and. wallamplit>1.d-10) then
          vwall=wallbs_rand(beter,wallamplit,xa,xb,xc,nmod_t,nmod_z)
        else
          vwall=0.d0
        endif
        !
        j=0
        !
        do k=0,km
        do i=0,im
          !
          pe=num1d3*(4.d0*prs(i,1,k)-prs(i,2,k))
          !
          if(x(i,j,k,1)<=xslip) then
            ue=num1d3*(4.d0*vel(i,1,k,1)-vel(i,2,k,1))
            te=num1d3*(4.d0*tmp(i,1,k)-tmp(i,2,k))
            !
            vel(i,j,k,1)=ue
            vel(i,j,k,2)=0.d0
            tmp(i,j,k)  =te
            !
          else
            !
            vel(i,j,k,1)=0.d0
            vel(i,j,k,2)=vwall(i,k)
            tmp(i,j,k)  =tw
            !
          endif
          !
          !
          vel(i,j,k,3)=0.d0
          prs(i,j,k)  =pe
          !
          do jspec=1,num_species
            spc(i,j,k,jspec)=num1d3*(4.d0*spc(i,1,k,jspec)-spc(i,2,k,jspec))
          enddo
          !
          rho(i,j,k)  =thermal(pressure=prs(i,j,k),temperature=tmp(i,j,k))
          !
          call fvar2q(      q=  q(i,j,k,:),                            &
                      density=rho(i,j,k),                              &
                     velocity=vel(i,j,k,:),                            &
                     pressure=prs(i,j,k),                              &
                      species=spc(i,j,k,:)                             )
          !
          if(trim(turbmode)=='k-omega') then
            tke(i,j,k)=0.d0
            !
            delta_d1=dis2point2(x(i,j,k,:),x(i,j+1,k,:))
            miu=miucal(tmp(i,j,k))/Reynolds
            omg(i,j,k)=50.d0*miu/rho(i,j,k)/beta1/delta_d1
            !
            q(i,j,k,5+num_species+1)=tke(i,j,k)*rho(i,j,k)
            q(i,j,k,5+num_species+2)=omg(i,j,k)*rho(i,j,k)
            !
          endif
            !
        enddo
        enddo
        !
      endif
      !
    elseif(ndir==4) then
      !
      if(jrk==jrkm) then
        !
        j=jm
        do k=0,km
        do i=0,im
          pe=num1d3*(4.d0*prs(i,j-1,k)-prs(i,j-2,k))
          !
          vel(i,j,k,1)=0.d0
          vel(i,j,k,2)=0.d0
          vel(i,j,k,3)=0.d0
          prs(i,j,k)  =pe
          tmp(i,j,k)  =tw
          !
          if(num_species>0) then
            spc(i,j,k,1)=num1d3*(4.d0*spc(i,j-1,k,1)-spc(i,j-2,k,1))
          endif
          !
          rho(i,j,k)  =thermal(pressure=prs(i,j,k),temperature=tmp(i,j,k))
          !
          call fvar2q(      q=  q(i,j,k,:),                            &
                      density=rho(i,j,k),                              &
                     velocity=vel(i,j,k,:),                            &
                     pressure=prs(i,j,k),                              &
                      species=spc(i,j,k,:)                             )
          !
          if(trim(turbmode)=='k-omega') then
            tke(i,j,k)=0.d0
            !
            delta_d1=dis2point2(x(i,j,k,:),x(i,j-1,k,:))
            miu=miucal(tmp(i,j,k))/Reynolds
            omg(i,j,k)=50.d0*miu/rho(i,j,k)/beta1/delta_d1
            !
            q(i,j,k,5+num_species+1)=tke(i,j,k)*rho(i,j,k)
            q(i,j,k,5+num_species+2)=omg(i,j,k)*rho(i,j,k)
            !
          endif
          !
        enddo
        enddo
        !
      endif
      !
    else
      stop ' !! ndir not defined @ noslip'
    endif
    !
  end subroutine slipisotwall
  !
  subroutine slipadibwall(ndir)
    !
    use commvar,   only : Reynolds,turbmode,lreport
    use commarray, only : tke,omg
    use fludyna,   only : thermal,fvar2q,q2fvar,miucal
    use commfunc,  only : dis2point2
    use parallel,  only : mpi_jmin,bcast
    !
    ! arguments
    integer,intent(in) :: ndir
    !
    ! local data
    integer :: i,j,k,l,jspec,fh
    real(8) :: pe,ue,ve,te,delta_d1,miu,beta1
    real(8) :: vwall(0:im,0:km)
    character(len=64) :: filewbs
    logical :: lexist
    !
    real(8),save :: beter,wallamplit,xa,xb,xc
    integer,save :: nmod_t,nmod_z
    logical,save :: lfirstcal=.true.
    !
    beta1=0.075d0
    !
    if(ndir==3) then
      !
      if(jrk==0) then
        !
        if(lreport) lfirstcal=.true.
        !
        if(lfirstcal) then
          !
          if(mpirank==0) then
            !
            filewbs='datin/wallbs.dat'
            !
            inquire(file=trim(filewbs), exist=lexist)
            !
            if(lexist) then
              fh=get_unit()
              open(fh,file=trim(filewbs),action='read')
              read(fh,*)
              read(fh,*)
              read(fh,*)wallamplit,beter,xa,xb,xc
              read(fh,*)
              read(fh,*)nmod_t,nmod_z
              close(fh)
              print*,' >> ',trim(filewbs)
              !
              print*,' ------------- wall blowing & suction parameters -------------'
              write(*,"(38x,A,1X,F12.5)")   '  amplitude:',wallamplit
              write(*,"(38x,A,1X,F12.5)")   '  frequency:',beter
              write(*,"(23x,3(A,1X,F12.5))")'   x extent:',xa,' ~',xb,' ~',xc
              write(*,"(42x,(A,1X,I0))")'    temporal modes:',nmod_t
              write(*,"(42x,(A,1X,I0))")'    spanwise modes:',nmod_z
              print*,' --------------------------------------------------------------'
            else
              wallamplit=0.d0
              beter=1.d0
              xa=0.d0
              xb=0.d0
              xc=0.d0
              nmod_t=0
              nmod_z=0
            endif
          endif
          !
          call bcast(wallamplit,comm=mpi_jmin)
          call bcast(beter,comm=mpi_jmin)
          call bcast(xa,comm=mpi_jmin)
          call bcast(xb,comm=mpi_jmin)
          call bcast(xc,comm=mpi_jmin)
          call bcast(nmod_t,comm=mpi_jmin)
          call bcast(nmod_z,comm=mpi_jmin)
          !
          !
          lfirstcal=.false.
          !
        endif
        !
        if(ndims==3 .and. wallamplit>1.d-10) then
          vwall=wallbs_rand(beter,wallamplit,xa,xb,xc,nmod_t,nmod_z)
        else
          vwall=0.d0
        endif
        !
        j=0
        !
        do k=0,km
        do i=0,im
          !
          pe=num1d3*(4.d0*prs(i,1,k)-prs(i,2,k))
          te=num1d3*(4.d0*tmp(i,1,k)-tmp(i,2,k))
          !
          ! if(x(i,j,k,1)<=xslip) then
          !   ue=num1d3*(4.d0*vel(i,1,k,1)-vel(i,2,k,1))
          !   !
          !   vel(i,j,k,1)=ue
          !   vel(i,j,k,2)=0.d0
          !   !
          ! else
          !   !
          !   vel(i,j,k,1)=0.d0
          !   vel(i,j,k,2)=0.d0
          !   !
          ! endif
          !
          ue=num1d3*(4.d0*vel(i,1,k,1)-vel(i,2,k,1))
          vel(i,j,k,1)=ue
          vel(i,j,k,2)=0.d0
          !
          vel(i,j,k,3)=0.d0
          tmp(i,j,k)  =te
          prs(i,j,k)  =pe
          !
          do jspec=1,num_species
            spc(i,j,k,jspec)=num1d3*(4.d0*spc(i,1,k,jspec)-spc(i,2,k,jspec))
          enddo
          !
          rho(i,j,k)  =thermal(pressure=prs(i,j,k),temperature=tmp(i,j,k))
          !
          call fvar2q(      q=  q(i,j,k,:),                            &
                      density=rho(i,j,k),                              &
                     velocity=vel(i,j,k,:),                            &
                     pressure=prs(i,j,k),                              &
                      species=spc(i,j,k,:)                             )
          !
          if(trim(turbmode)=='k-omega') then
            tke(i,j,k)=0.d0
            !
            delta_d1=dis2point2(x(i,j,k,:),x(i,j+1,k,:))
            miu=miucal(tmp(i,j,k))/Reynolds
            omg(i,j,k)=50.d0*miu/rho(i,j,k)/beta1/delta_d1
            !
            q(i,j,k,5+num_species+1)=tke(i,j,k)*rho(i,j,k)
            q(i,j,k,5+num_species+2)=omg(i,j,k)*rho(i,j,k)
            !
          endif
            !
        enddo
        enddo
        !
      endif
      !
    elseif(ndir==4) then
      !
      if(jrk==jrkm) then
        !
        j=jm
        !
        do k=0,km
        do i=0,im
          pe=num1d3*(4.d0*prs(i,j-1,k)-prs(i,j-2,k))
          te=num1d3*(4.d0*tmp(i,j-1,k)-tmp(i,j-2,k))
          !
          ue=num1d3*(4.d0*vel(i,j-1,k,1)-vel(i,j-2,k,1))
          ve=num1d3*(4.d0*vel(i,j-1,k,2)-vel(i,j-2,k,2))
          !
          ! vel(i,j,k,1)=0.d0
          vel(i,j,k,1)=ue
          vel(i,j,k,2)=ve
          !
          vel(i,j,k,3)=0.d0
          prs(i,j,k)  =pe
          tmp(i,j,k)  =te
          !
          if(num_species>0) then
            spc(i,j,k,1)=num1d3*(4.d0*spc(i,j-1,k,1)-spc(i,j-2,k,1))
          endif
          !
          rho(i,j,k)  =thermal(pressure=prs(i,j,k),temperature=tmp(i,j,k))
          !
          call fvar2q(      q=  q(i,j,k,:),                            &
                      density=rho(i,j,k),                              &
                     velocity=vel(i,j,k,:),                            &
                     pressure=prs(i,j,k),                              &
                      species=spc(i,j,k,:)                             )
          !
          if(trim(turbmode)=='k-omega') then
            tke(i,j,k)=0.d0
            !
            delta_d1=dis2point2(x(i,j,k,:),x(i,j-1,k,:))
            miu=miucal(tmp(i,j,k))/Reynolds
            omg(i,j,k)=50.d0*miu/rho(i,j,k)/beta1/delta_d1
            !
            q(i,j,k,5+num_species+1)=tke(i,j,k)*rho(i,j,k)
            q(i,j,k,5+num_species+2)=omg(i,j,k)*rho(i,j,k)
            !
          endif
          !
        enddo
        enddo
        !
      endif
      !
    else
      stop ' !! ndir not defined @ noslip'
    endif
    !
  end subroutine slipadibwall
  !+-------------------------------------------------------------------+
  !| The end of the subroutine noslip.                                 |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This function is to get wall normall blowing and suction velocity.|
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 27-09-2021: Created by J. Fang @ Warrington                       |
  !+-------------------------------------------------------------------+
  function wallbs(beter,wallamplit,xa,xb,nmod_t,nmod_z) result(vwall)
    !
    ! arguments
    real(8),intent(in) :: beter,wallamplit,xa,xb
    integer,intent(in) :: nmod_t,nmod_z
    real(8) :: vwall(0:im,0:km)
    !
    ! local data
    integer :: i,k,l,m,seed_size
    real(8) :: theter,fx,gz,ht,zl,tm
    integer,allocatable :: seed(:)
    !
    real(8),save :: sqrt27,z0,t0,lz,rfluc,rampl
    real(8),save :: randomv(15)
    logical,save :: lfirstcal=.true.
    !
    if(lfirstcal) then
      !
      ! beter=0.02d0
      ! !
      ! xa=5.d0
      ! xb=40.d0
      ! !
      ! nmod_z=3
      ! nmod_t=2
      !
      lz=zmax-zmin
      !
      sqrt27=1.d0/sqrt(27.d0)
      z0=0.2d0/(1.d0-0.8d0**nmod_z)
      !
      if(nmod_t==0) then
        t0=0
      else
        t0=0.2d0/(1.d0-0.8d0**nmod_t)
      endif
      !
      ! call random_seed() ! initialize with system generated seed
      call random_seed(size=seed_size) ! find out size of seed
      allocate(seed(seed_size))
      ! call random_seed(get=seed) ! get system generated seed
      ! write(*,*) seed            ! writes system generated seed
      seed=0
      call random_seed(put=seed) ! set current seed
      ! call random_seed(get=seed) ! get current seed
      ! write(*,*) seed            ! writes 0
      deallocate(seed)           ! safe
      !
      do m=1,15
        call random_number(randomv(m))
        ! write(*,*)mpirank,'|',m,randomv(m)
      end do
      !
      rampl=0.05d0
      !
      call random_seed(size=seed_size) ! find out size of seed
      allocate(seed(seed_size))
      seed=mpirank
      call random_seed(put=seed) 
      deallocate(seed)
      !
      lfirstcal=.false.
      !
    endif
    !
    do k=0,km
    do i=0,im
      !
      if(x(i,0,k,1)<=xb .and. x(i,0,k,1)>=xa) then
        !
        theter=2.d0*pi*(x(i,0,k,1)-xa)/(xb-xa)
        fx=4.d0*dsin(theter)*(1.d0-dcos(theter))*sqrt27
        !
        if(km==0) then
          gz=1.d0
        else
          gz=0.d0
          do l=1,nmod_z
            !
            if(l==1) then
              zl=z0
            else
              zl=zl*0.8d0
            end if
            !
            gz=gz+zl*dsin(2.d0*pi*l*(x(i,0,k,3)/lz+randomv(l)))
            !
          end do
        endif
        !
        if(nmod_t==0) then
          ht=1.d0
        else
          ht=0.d0
          do m=1,nmod_t
            !
            if(m==1) then
              tm=t0
            else
              tm=tm*0.8d0
            end if
            !
            ht=ht+tm*dsin(2.d0*pi*m*(beter*time+randomv(m+10)))
            !
          end do
        endif
        !
        vwall(i,k)=wallamplit*uinf*fx*gz*ht
        !
      else
        vwall(i,k)=0.d0
      end if
      !
    end do
    end do
    !
    return
    !
  end function wallbs
  !
  function wallbs_rand(beter,wallamplit,xa,xb,xc,nmod_t,nmod_z) result(vwall)
    !
    ! arguments
    real(8),intent(in) :: beter,wallamplit,xa,xb,xc
    integer,intent(in) :: nmod_t,nmod_z
    real(8) :: vwall(0:im,0:km)
    !
    ! local data
    integer :: i,k,l,m,seed_size
    real(8) :: theter,fx,gz,ht,zl,tm
    integer,allocatable :: seed(:)
    !
    real(8),save :: sqrt27,z0,t0,lz,rfluc,rampl
    real(8),save :: randomv(15)
    logical,save :: lfirstcal=.true.
    !
    if(lfirstcal) then
      !
      ! beter=0.02d0
      ! !
      ! xa=5.d0
      ! xb=40.d0
      ! !
      ! nmod_z=3
      ! nmod_t=2
      !
      lz=zmax-zmin
      !
      sqrt27=1.d0/sqrt(27.d0)
      z0=0.2d0/(1.d0-0.8d0**nmod_z)
      !
      if(nmod_t==0) then
        t0=0
      else
        t0=0.2d0/(1.d0-0.8d0**nmod_t)
      endif
      !
      ! call random_seed() ! initialize with system generated seed
      call random_seed(size=seed_size) ! find out size of seed
      allocate(seed(seed_size))
      ! call random_seed(get=seed) ! get system generated seed
      ! write(*,*) seed            ! writes system generated seed
      seed=0
      call random_seed(put=seed) ! set current seed
      ! call random_seed(get=seed) ! get current seed
      ! write(*,*) seed            ! writes 0
      deallocate(seed)           ! safe
      !
      do m=1,15
        call random_number(randomv(m))
        ! write(*,*)mpirank,'|',m,randomv(m)
      end do
      !
      rampl=0.1d0
      !
      call random_seed(size=seed_size) ! find out size of seed
      allocate(seed(seed_size))
      seed=mpirank
      call random_seed(put=seed) 
      deallocate(seed)
      !
      lfirstcal=.false.
      !
    endif
    !
    do k=0,km
    do i=0,im
      !
      if(x(i,0,k,1)<=xb .and. x(i,0,k,1)>=xa) then
        !
        theter=2.d0*pi*(x(i,0,k,1)-xa)/(xb-xa)
        fx=4.d0*dsin(theter)*(1.d0-dcos(theter))*sqrt27
        !
        gz=dsin(2.d0*nmod_z*pi*(x(i,0,k,3)/lz))
        !
        ht=1.d0
        !
        call random_number(rfluc)
        rfluc=(rfluc*2.d0-1.d0)*rampl
        !
        vwall(i,k)=wallamplit*uinf*fx*gz*ht*(1.d0+rfluc)
        !
      elseif(x(i,0,k,1)<=xc .and. x(i,0,k,1)>=xb) then
        !
        theter=2.d0*pi*(x(i,0,k,1)-xb)/(xc-xb)
        fx=4.d0*dsin(theter)*(1.d0-dcos(theter))*sqrt27
        !
        gz=dsin(2.d0*nmod_z*pi*(x(i,0,k,3)/lz)+0.5d0*pi)
        !
        ht=1.d0
        !
        call random_number(rfluc)
        rfluc=(rfluc*2.d0-1.d0)*rampl
        !
        vwall(i,k)=wallamplit*uinf*fx*gz*ht*(1.d0+rfluc)
        !
      else
        vwall(i,k)=0.d0
      end if
      !
    end do
    end do
    !
    return
    !
  end function wallbs_rand
  !+-------------------------------------------------------------------+
  !| The end of the subroutine wallbs.                                 |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to apply a fixed b.c.                          |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 09-05-2022: Created by J. Fang @ Warrington                       |
  !+-------------------------------------------------------------------+
  subroutine fixbc(ndir)
    !
    use fludyna, only : fvar2q,thermal
    !
    ! arguments
    integer,intent(in) :: ndir
    !
    ! local data
    integer :: i,j,k
    !
    if(ndir==1 .and. irk==0) then
      !
      i=0
      !
      do k=0,km
      do j=0,jm
      enddo
      enddo
      !
    elseif(ndir==2 .and. irk==irkm) then
      !
      i=im
      !
      do k=0,km
      do j=0,jm
      enddo
      enddo
      !
    elseif(ndir==3 .and. jrk==0) then
      !
      j=0
      !
      do k=0,km
      do i=0,im
        !
        rho(i,j,k)  =2.d0
        vel(i,j,k,:)=0.d0
        prs(i,j,k)  =1.d0
        !
        tmp(i,j,k)  =thermal(density=rho(i,j,k),pressure=prs(i,j,k))
        !
        call fvar2q(      q=  q(i,j,k,:),   density=rho(i,j,k),        &
                   velocity=vel(i,j,k,:),  pressure=prs(i,j,k)         )
        !
      enddo
      enddo
      !
    elseif(ndir==4 .and. jrk==jrkm) then
      !
      j=jm
      !
      do k=0,km
      do i=0,im
        !
        rho(i,j,k)  =1.d0
        vel(i,j,k,:)=0.d0
        prs(i,j,k)  =2.5d0
        !
        tmp(i,j,k)  =thermal(density=rho(i,j,k),pressure=prs(i,j,k))
        !
        call fvar2q(      q=  q(i,j,k,:),   density=rho(i,j,k),        &
                   velocity=vel(i,j,k,:),  pressure=prs(i,j,k)         )
        !
      enddo
      enddo
      !
    elseif(ndir==5 .and. krk==0) then
      !
      k=0
      !
      do j=0,jm
      do i=0,im
      enddo
      enddo
      !
    elseif(ndir==6 .and. krk==krkm) then
      !
      k=km
      !
      do j=0,jm
      do i=0,im
      enddo
      enddo
      !
    endif
    !
  end subroutine fixbc
  !+-------------------------------------------------------------------+
  !| The end of the subroutine fixbc.                                  |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to obtain the inlet flow for jet.              |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 24-02-2021: Created by J. Fang @ Warrington                       |
  !+-------------------------------------------------------------------+
  subroutine jetinflow
    !
    use fludyna, only : jetvel,thermal
    !
    ! local data
    integer :: i,j,k
    real(8) :: radi
    !
    i=0
    do k=0,km
    do j=0,jm
      !
      rho_in(j,k)  =roinf
      !
      radi=sqrt(x(i,j,k,2)**2+x(i,j,k,3)**2)
      ! radi=abs(x(i,j,k,2))
      !
      vel_in(j,k,:)=jetvel(radi)
      tmp_in(j,k)  =tinf
      prs_in(j,k)  =thermal(density=rho_in(j,k),temperature=tmp_in(j,k))
      !
      if(num_species>0) then
        if(vel_in(j,k,1)>(1.d0+1.d-10)*uinf) then
          spc_in(j,k,1)=1.d0
        else
          spc_in(j,k,1)=0.d0
        endif
        !
        spc_in(j,k,2)=1.d0-spc_in(j,k,1)
      endif
      !
    enddo
    enddo
    !
  end subroutine jetinflow
  !+-------------------------------------------------------------------+
  !| The end of the subroutine jetinflow.                              |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to obtain the inlet flow for mixing layer.     |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 19-03-2021: Created by J. Fang @ Warrington                       |
  !+-------------------------------------------------------------------+
  subroutine mixlayerinflow
    !
    use fludyna, only : mixinglayervel,thermal
    !
    ! local data
    integer :: i,j,k
    real(8) :: radi
    !
    i=0
    do k=0,km
    do j=0,jm
      !
      rho_in(j,k)  =roinf
      vel_in(j,k,:)=mixinglayervel(x(i,j,k,2))
      tmp_in(j,k)  =tinf
      prs_in(j,k)  =thermal(density=rho_in(j,k),temperature=tmp_in(j,k))
      !
      if(num_species>0) then
        if(vel_in(j,k,1)>(1.d0+1.d-10)*uinf) then
          spc_in(j,k,1)=1.d0
        else
          spc_in(j,k,1)=0.d0
        endif
        !
        spc_in(j,k,2)=1.d0-spc_in(j,k,1)
      endif
      !
    enddo
    enddo
    !
  end subroutine mixlayerinflow
  !+-------------------------------------------------------------------+
  !| The end of the subroutine mixlayerinflow.                         |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to obtain the inlet flow by reading a inlet    |
  !| profile.                                                          |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 27-09-2021: Created by J. Fang @ Warrington                       |
  !+-------------------------------------------------------------------+
  subroutine profileinflow
    !
    use fludyna, only : mixinglayervel,thermal
    !
    ! local data
    integer :: i,j,k
    real(8) :: radi
    ! real(8),allocatable
    !
    i=0
    do k=0,km
    do j=0,jm
      !
      rho_in(j,k)  =rho_prof(j)
      vel_in(j,k,1)=vel_prof(j,1)
      vel_in(j,k,2)=vel_prof(j,2)
      vel_in(j,k,3)=0.d0
      tmp_in(j,k)  =tmp_prof(j)
      prs_in(j,k)  =thermal(density=rho_in(j,k),temperature=tmp_in(j,k))
      !
      if(num_species>0) then
        spc_in(j,k,:)=1.d0
      endif
      !
    enddo
    enddo
    !
  end subroutine profileinflow
  !+-------------------------------------------------------------------+
  !| The end of the subroutine profileinflow.                          |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to obtain the inlet of freestream flow.        |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 25-02-2021: Created by J. Fang @ Warrington                       |
  !+-------------------------------------------------------------------+
  subroutine freestreaminflow
    !
    use fludyna, only : jetvel,thermal
    !
    ! local data
    integer :: i,j,k,jspec
    !
    i=0
    do k=0,km
    do j=0,jm
      !
      rho_in(j,k)  =roinf
      vel_in(j,k,1)=uinf
      vel_in(j,k,2)=vinf
      vel_in(j,k,3)=winf
      tmp_in(j,k)  =tinf
      prs_in(j,k)  =thermal(density=rho_in(j,k),temperature=tmp_in(j,k))
      !
      do jspec=1,num_species
        spc_in(j,k,jspec)=0.d0
      enddo
      !
    enddo
    enddo
    !
  end subroutine freestreaminflow
  !+-------------------------------------------------------------------+
  !| The end of the subroutine freestreaminflow.                       |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to implement NSCBC for ghost cells.            |
  !| pq1-----pq2---|---gq1-----gq2                                     |
  !+-------------------------------------------------------------------+
  !| ref: Motheau, Almgreny & Bell, Navier-StokesCharacteristic        |
  !|      Boundary Conditions Using Ghost Cells, AIAA J., 2017, 55(10) |
  !| ref: Gross & Fasel, Characteristic Ghost-Cell Boundary Condition, |
  !|      AIAA J., 2007, 45(1).                                        |
  !| ref: Yoo,et al. Characteristic boundary conditions for direct     |
  !|      simulations of turbulent counterflow flames. Combustion      |
  !|      Theory and Modelling, 2005, 9(4): 617-646.                   |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 14-Jul-2020  | Created by J. Fang @ STFC Daresbury Laboratory     |
  !+-------------------------------------------------------------------+
  subroutine gcnscbc(ndir,prs_t)
    !
    use fludyna,   only : fvar2q,q2fvar,thermal,sos
    use commfunc,  only : deriv
    !
    ! arguments
    integer,intent(in) :: ndir
    real(8),intent(in) :: prs_t
    !
    ! local data
    real(8),parameter :: eps=1.d-10
    integer :: vnorm,css,css2,i,j,k,jspc
    real(8) :: mavg,lambda(numq)
    real(8) :: dvel(3),dprs,drho,dtmp,dspc(num_species),              &
               lodi(5+num_species),tran(5)
    real(8) :: dx,sigma,beter,kamma,dpdy,dudy,dvdy,dy,eta
    !
    if(ndir==2 .and. irk==irkm) then
      !
      i=im
      do k=0,km
      do j=0,jm
        !
        css=sos(tmp(i-1,j,k))
        css2=css*css
        vnorm=vel(i-1,j,k,1)
        !
        mavg=vnorm/css
        !
        ! characteristic velocities
        lambda(1)        =vnorm-css
        lambda(2:4)      =vnorm
        lambda(5)        =vnorm+css
        lambda(5:numq)   =vnorm ! for scalars
        !
        lambda = lambda + eps
        !
        dx=(x(i-1,j,k,1)-x(i-2,j,k,1))
        dvel(1)=-deriv(vel(i-1,j,k,1),vel(i-2,j,k,1) )/dx
        dvel(2)=-deriv(vel(i-1,j,k,2),vel(i-2,j,k,2) )/dx
        dvel(3)=-deriv(vel(i-1,j,k,3),vel(i-2,j,k,3) )/dx
        !
        dprs   =-deriv(prs(i-1,j,k),prs(i-2,j,k) )/dx
        drho   =-deriv(rho(i-1,j,k),rho(i-2,j,k) )/dx
        dtmp   =-deriv(tmp(i-1,j,k),tmp(i-2,j,k) )/dx
        !
        if(num_species>0) then
          do jspc=1,num_species
            dspc(jspc)=-deriv(spc(i-1,j,k,jspc),spc(i-2,j,k,jspc) )/dx
          enddo
        endif
        !
        ! if(jrk==0 .and. j==0) then
        !   dy=x(i-1,j+1,k,2)-x(i-1,j,k,2)
        !   dpdy=(prs(i-1,j+1,k)  -prs(i-1,j,k))/dy
        !   dudy=(vel(i-1,j+1,k,1)-vel(i-1,j,k,1))/dy
        !   dvdy=(vel(i-1,j+1,k,2)-vel(i-1,j,k,2))/dy
        ! elseif(jrk==jrkm .and. j==jm) then
        !   dy=x(i-1,j,k,2)-x(i-1,j-1,k,2)
        !   dpdy=(prs(i-1,j,k)-prs(i-1,j-1,k))/dy
        !   dudy=(vel(i-1,j,k,1)-vel(i-1,j-1,k,1))/dy
        !   dvdy=(vel(i-1,j,k,2)-vel(i-1,j-1,k,2))/dy
        ! else
        !   dy=0.5d0*(x(i-1,j+1,k,2)-x(i-1,j-1,k,2))
        !   dpdy=0.5d0*(prs(i-1,j+1,k)-prs(i-1,j-1,k))/dy
        !   dudy=0.5d0*(vel(i-1,j+1,k,1)-vel(i-1,j-1,k,1))/dy
        !   dvdy=0.5d0*(vel(i-1,j+1,k,2)-vel(i-1,j-1,k,2))/dy
        ! endif
        !
        ! LODI with 1-side difference
        lodi(1) = lambda(1) *(dprs-rho(i-1,j,k)*css*dvel(1))
        lodi(2) = lambda(2) *(css2*drho-dprs)
        lodi(3) = lambda(3) *dvel(2)
        lodi(4) = lambda(4) *dvel(3)
        lodi(5) = lambda(5) *(dprs+rho(i-1,j,k)*css*dvel(1))
        !
        if(num_species>0) then
          do jspc=1,num_species
            lodi(5+jspc)=lambda(5+jspc)*dspc(jspc)
          enddo
        endif
        !
        ! soft outflow condition
        sigma=0.25d0
        beter=0.5d0
        kamma=sigma*css*(1.d0-mavg**2)/(xmax-xmin)
        !
        ! tran(1)=vel(i-1,j,k,2)*(dpdy-rho(i-1,j,k)*css*dudy)+           &
        !                                          gamma*prs(i-1,j,k)*dvdy
        !
        lodi(1)=kamma*(prs(i-1,j,k)-prs_t) !-(1.d0-beter)*tran(1)
        !
        ! project back to physical gradient
        ! ref: Eqs(33-36), Poinsot, T. J., Lele, S. K. 1992 Boundary
        !|Conditions for Direct Simulations of Compressible Viscous Flows.
        ! J COMPUT PHYS. 101.
        !
        lodi(:) = lodi(:)/lambda(:)
        !
        drho = 1.d0/css**2*(lodi(2)+0.5d0*(lodi(1)+lodi(5)))
        dtmp = tmp(i-1,j,k)/rho(i-1,j,k)/css*(-lodi(2)+0.5d0*(gamma-1.d0)*(lodi(1)+lodi(5)))
        dprs = 0.5d0*(lodi(1)+lodi(5))
        dvel(1) = 1.d0/(2.d0*rho(i-1,j,k)*css)*(lodi(5)-lodi(1))
        dvel(2) = lodi(3)
        dvel(3) = lodi(4)
        dspc(:) = lodi(5:5+num_species)
        !
        tmp(i,j,k)  =tmp(i-1,j,k)  +2.d0*dx*dtmp
        rho(i,j,k)  =rho(i-1,j,k)  +2.d0*dx*drho
        prs(i,j,k)  =prs(i-1,j,k)  +2.d0*dx*dprs
        vel(i,j,k,:)=vel(i-1,j,k,:)+2.d0*dx*dvel(:)
        spc(i,j,k,:)=spc(i-1,j,k,:)+2.d0*dx*dspc(:)
        !
        rho(i,j,k)  =thermal(pressure=prs(i,j,k),temperature=tmp(i,j,k))
        !
        call fvar2q(      q=  q(i,j,k,:),   density=rho(i,j,k),        &
                   velocity=vel(i,j,k,:),  pressure=prs(i,j,k),        &
                    species=spc(i,j,k,:)                               )
        !
        qrhs(i,j,k,:)=0.d0
        !
      enddo
      enddo
      !
    endif
    !
    if(ndir==4 .and. jrk==jrkm) then
      !
      j=jm
      do k=0,km
      do i=0,im
        !
        css=sos(tmp(i,j-1,k))
        css2=css*css
        vnorm=vel(i,j-1,k,2)
        !
        mavg=vnorm/css
        !
        ! characteristic velocities
        lambda(1)        =vnorm-css
        lambda(2:4)      =vnorm
        lambda(5)        =vnorm+css
        lambda(5:numq)   =vnorm ! for scalars
        !
        lambda = lambda + eps
        !
        dy=(x(i,j-1,k,2)-x(i,j-2,k,2))
        dvel(1)=-deriv(vel(i,j-1,k,1),vel(i,j-2,k,1) )/dy
        dvel(2)=-deriv(vel(i,j-1,k,2),vel(i,j-2,k,2) )/dy
        dvel(3)=-deriv(vel(i,j-1,k,3),vel(i,j-2,k,3) )/dy
        !
        dprs   =-deriv(prs(i,j-1,k),prs(i,j-2,k) )/dy
        drho   =-deriv(rho(i,j-1,k),rho(i,j-2,k) )/dy
        dtmp   =-deriv(tmp(i,j-1,k),tmp(i,j-2,k) )/dy
        !
        if(num_species>0) then
          do jspc=1,num_species
            dspc(jspc)=-deriv(spc(i,j-1,k,jspc),spc(i,j-2,k,jspc) )/dy
          enddo
        endif
        !
        ! LODI with 1-side difference
        lodi(1) = lambda(1) *(dprs-rho(i,j-1,k)*css*dvel(2))
        lodi(2) = lambda(2) *(css2*drho-dprs)
        lodi(3) = lambda(3) *dvel(1)
        lodi(4) = lambda(4) *dvel(3)
        lodi(5) = lambda(5) *(dprs+rho(i,j-1,k)*css*dvel(2))
        !
        if(num_species>0) then
          do jspc=1,num_species
            lodi(5+jspc)=lambda(5+jspc)*dspc(jspc)
          enddo
        endif
        !
        ! soft outflow condition
        ! tran(1)=vel(i-1,j,k,2)*(dpdy-rho(i-1,j,k)*css*dudy)+           &
        !                                          gamma*prs(i-1,j,k)*dvdy
        !
        if(vnorm>=0.d0) then
          sigma=0.25d0
          beter=0.5d0
          kamma=sigma*css*(1.d0-mavg**2)/(ymax-ymin)
          !
          lodi(1)=kamma*(prs(i,j-1,k)-prs_t) !-(1.d0-beter)*tran(1)
        else
          eta=0.25d0
          !
          lodi(1)=eta*rho(i,j-1,k)*css2*(1.d0-mavg**2)/(ymax-ymin)*(vel(i,j-1,k,2)-vinf)
          lodi(2)=eta*rho(i,j-1,k)*css2/gamma/(ymax-ymin)*(tmp(i,j-1,k)-tinf)
          lodi(3)=eta*css/(ymax-ymin)*(vel(i,j-1,k,1)-uinf)
          lodi(4)=eta*css/(ymax-ymin)*(vel(i,j-1,k,3)-winf)
          do jspc=1,num_species
            lodi(5+jspc)=eta*css/(ymax-ymin)*(spc(i,j-1,k,jspc)-0.d0)
          enddo
        endif
        !
        ! project back to physical gradient
        ! ref: Eqs(33-36), Poinsot, T. J., Lele, S. K. 1992 Boundary
        !|Conditions for Direct Simulations of Compressible Viscous Flows.
        ! J COMPUT PHYS. 101.
        !
        lodi(:) = lodi(:)/lambda(:)
        !
        drho = 1.d0/css**2*(lodi(2)+0.5d0*(lodi(1)+lodi(5)))
        dtmp = tmp(i,j-1,k)/rho(i,j-1,k)/css*(-lodi(2)+0.5d0*(gamma-1.d0)*(lodi(1)+lodi(5)))
        dprs = 0.5d0*(lodi(1)+lodi(5))
        dvel(2) = 1.d0/(2.d0*rho(i,j-1,k)*css)*(lodi(5)-lodi(1))
        dvel(1) = lodi(3)
        dvel(3) = lodi(4)
        dspc(:) = lodi(5:5+num_species)
        !
        tmp(i,j,k)  =tmp(i,j-1,k)  +2.d0*dy*dtmp
        rho(i,j,k)  =rho(i,j-1,k)  +2.d0*dy*drho
        prs(i,j,k)  =prs(i,j-1,k)  +2.d0*dy*dprs
        vel(i,j,k,:)=vel(i,j-1,k,:)+2.d0*dy*dvel(:)
        spc(i,j,k,:)=spc(i,j-1,k,:)+2.d0*dy*dspc(:)
        !
        rho(i,j,k)  =thermal(pressure=prs(i,j,k),temperature=tmp(i,j,k))
        !
        call fvar2q(      q=  q(i,j,k,:),   density=rho(i,j,k),        &
                   velocity=vel(i,j,k,:),  pressure=prs(i,j,k),        &
                    species=spc(i,j,k,:)                               )
        !
        qrhs(i,j,k,:)=0.d0
        !
      enddo
      enddo
      !
    endif
    !
  end subroutine gcnscbc
  !+-------------------------------------------------------------------+
  !| The end of the subroutine gcnscbc.                                |
  !+-------------------------------------------------------------------+
  !!
  !+-------------------------------------------------------------------+
  !| This subroutine is a templet of designing b.c.                    |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 09-05-2022: Created by J. Fang @ Warrington                       |
  !+-------------------------------------------------------------------+
  subroutine bc_template(ndir)
    !
    ! arguments
    integer,intent(in) :: ndir
    !
    ! local data
    integer :: i,j,k
    !
    if(ndir==1 .and. irk==0) then
      !
      i=0
      !
      do k=0,km
      do j=0,jm
      enddo
      enddo
      !
    elseif(ndir==2 .and. irk==irkm) then
      !
      i=im
      !
      do k=0,km
      do j=0,jm
      enddo
      enddo
      !
    elseif(ndir==3 .and. jrk==0) then
      !
      j=0
      !
      do k=0,km
      do i=0,im
      enddo
      enddo
      !
    elseif(ndir==4 .and. jrk==jrkm) then
      !
      j=jm
      !
      do k=0,km
      do i=0,im
      enddo
      enddo
      !
    elseif(ndir==5 .and. krk==0) then
      !
      k=0
      !
      do j=0,jm
      do i=0,im
      enddo
      enddo
      !
    elseif(ndir==6 .and. krk==krkm) then
      !
      k=km
      !
      do j=0,jm
      do i=0,im
      enddo
      enddo
      !
    endif
    !
  end subroutine bc_template
  !+-------------------------------------------------------------------+
  !| The end of the subroutine bc_template.                            |
  !+-------------------------------------------------------------------+
  !!
end module bc
!+---------------------------------------------------------------------+
!| The end of the module bc.                                           |
!+---------------------------------------------------------------------+
