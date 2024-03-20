!+---------------------------------------------------------------------+
!| This module contains subroutines of dealing with immersed boundary  |
!| method geometrically                                                |
!+---------------------------------------------------------------------+
!| CHANGE RECORD                                                       |
!| -------------                                                       |
!| 23-06-2022  | Created by J. Fang                                    |
!+---------------------------------------------------------------------+
module ibmethod
  !
  use parallel, only : mpirankname,mpistop,mpirank,lio,ptime,ig0,jg0,kg0
  use stlaio,  only: get_unit
  use commcal, only: ijkcellin,ijkin
  use stlaio,  only: get_unit
  !
  implicit none
  !
  contains
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to pre-process immersed solid for ib metho     |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 23-06-2022  | Created by J. Fang @ Warrington                     |
  !+-------------------------------------------------------------------+
  subroutine ibprocess
    !
    use commvar,   only : limmbou,solidfile,ibmode
    use readwrite, only : readsolid
    !
    if(limmbou .and. ibmode=='stl') then
      !
      if(mpirank==0) call readsolid(solidfile)
      !
      call solidgeom
      !
    endif
    !
  end subroutine ibprocess
  !+-------------------------------------------------------------------+
  !| The end of the subroutine ibprocess.                              |
  !+-------------------------------------------------------------------+
  !!
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
    use parallel,  only : bcast
    use tecio,     only : tecsolid
    use geom,      only : solidsilce,solidrange,solidresc,solidrota,solidshif
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
        call solidrange(immbody(js))
        ! call solidrange(immbody(js),inputcmd='checkdomain')
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
          call solidsilce(immbody(js),zsec=0.d0)
          !
          call solidrange(immbody(js),inputcmd='edge')
          !
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
    call bcast(nsolid)
    !
    call bcast(immbody)
    !
    return
    !
  end subroutine solidgeom
  !+-------------------------------------------------------------------+
  !| The end of the subroutine solidgeom.                              |
  !+-------------------------------------------------------------------+
  !!
  !+-------------------------------------------------------------------+
  !| This subroutine is to interpolate flowfield to the solid surface. |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 06-Jul-2022: Created by J. Fang @ Warrington.                     |
  !+-------------------------------------------------------------------+
  subroutine ibforce
    !
    use constdef
    use commtype,  only : solid,triangle,lsegment
    use commvar,   only : xmax,xmin,ymax,ymin,zmax,zmin,               &
                          immbody,nsolid,ndims,numq,im,jm,km,reynolds, &
                          uinf,pinf,roinf,nstep,time
    use geom,      only : pngrid2d
    use parallel,  only : psum,pmax
    use fludyna,   only : q2fvar,miucal
    !
    ! local data
    type(solid),pointer :: pso
    type(lsegment),pointer :: pedge
    real(8) :: xnode(3) 
    integer :: jsd,jedge,ijk(3)
    real(8),allocatable :: qib(:,:),dv(:,:,:),vel(:,:),prs(:),rho(:),tmp(:)
    real(8) :: drag,lift,pres_drag,fric_drag,pres_lift,fric_lift,aera,     &
               div,s11,s22,s12,miu,fx,fy,px,py
    integer,save :: counter=0,fh
    !
    do jsd=1,nsolid
      !
      pso=>immbody(jsd)
      !
      allocate( qib(size(pso%edge),numq),dv(size(pso%edge),1:ndims,1:ndims) )
      !
      if(ndims==2) then
        !
        do jedge=1,size(pso%edge)
          !
          xnode(1:2)=pso%edge(jedge)%cen
          xnode(3)  =0.d0
          !
          call pngrid2d(xnode,ijk)
          !
          ijk(1)=ijk(1)-ig0
          ijk(2)=ijk(2)-jg0
          ijk(3)=ijk(3)-kg0
          !
          ! if(ijk(1)>im .or. ijk(2)>jm) then
          !   print*,mpirank,'|',xnode,'-',ijk(1:2)
          ! endif
          !
          if(ijk(1)>0 .and. ijk(2)>0) then
            ! only do calculation for a phyiscal cel
            call interp_cell(xnode,ijk,qint=qib(jedge,:),dvint=dv(jedge,:,:))
          else
            qib(jedge,:) =-1.d10
            dv(jedge,:,:)=-1.d10
          endif
          !
          ! if(ijk(1)>0 .and. ijk(2)>0) then
          !   print*,mpirank,'|',xnode,'-',ijk(1:2)
          ! endif
          !
          ! if(jedge==217) then
          !   print*,mpirank,'|',pso%edge(jedge)%cen,ijk
          ! endif
          !
        enddo
        !
      elseif(ndims==3) then
        stop ' !! nothing yet @ ibforce'
      endif
      !
      qib=pmax(qib)
      dv =pmax(dv)
      !
      if(mpirank==0) then
        !
        allocate( vel(size(pso%edge),3),prs(size(pso%edge)),     &
                  rho(size(pso%edge)),  tmp(size(pso%edge)))
        !
        do jedge=1,size(pso%edge)
          call q2fvar(q=qib(jedge,:),    density =rho(jedge),    &
                                         velocity=vel(jedge,:),  &
                                         pressure=prs(jedge),    &
                                      temperature=tmp(jedge)      )
        enddo
        !
        ! not integrate the force
        drag=0.d0
        lift=0.d0
        pres_drag=0.d0
        pres_lift=0.d0
        fric_drag=0.d0
        fric_lift=0.d0
        aera=0.d0
        !
        do jedge=1,size(pso%edge)
          !
          pedge=>pso%edge(jedge)
          !
          aera=aera+pedge%length
          !
          div=num1d3*(dv(jedge,1,1)+dv(jedge,2,2))
          s11=2.d0*dv(jedge,1,1)-div
          s22=2.d0*dv(jedge,2,2)-div
          s12=(dv(jedge,1,2)+dv(jedge,2,1))
          !
          miu=miucal(tmp(jedge))/reynolds
          !
          ! print*,jedge,s11,s22,s12,prs(jedge)
          !
          fx=miu*pedge%length*(s11*pedge%normdir(1)+s12*pedge%normdir(2))
          fy=miu*pedge%length*(s12*pedge%normdir(1)+s22*pedge%normdir(2))
          px=-prs(jedge)*pedge%length*pedge%normdir(1)
          py=-prs(jedge)*pedge%length*pedge%normdir(2)
          !
          pres_drag=pres_drag+px
          pres_lift=pres_lift+py
          fric_drag=fric_drag+fx
          fric_lift=fric_lift+fy
          !
        enddo
        !
        pres_drag=pres_drag/(0.5d0*roinf*uinf**2)
        pres_lift=pres_lift/(0.5d0*roinf*uinf**2)
        fric_drag=fric_drag/(0.5d0*roinf*uinf**2)
        fric_lift=fric_lift/(0.5d0*roinf*uinf**2)
        !
        drag=pres_drag+fric_drag
        lift=pres_lift+fric_lift
        !
        if(counter==0 .and. nstep==0) then
          fh=get_unit()
          open(fh,file='surface_variables.dat')
          write(fh,"(A7,1X,A13,2(1X,A20))")'nstep','time','lift','drag'
        else
          fh=get_unit()
          open(fh,file='surface_variables.dat',access='append')
        endif
        !
        write(fh,"(I7,1X,E13.6E2,2(1X,E20.13E2))")nstep,time,lift,drag
        !
        close(fh)
        !
        ! print*,' ** aera of airfoil:',aera
        ! print*,' ** drag on airfoil:',drag
        ! print*,' **   friction drag:',fric_drag
        ! print*,' **   pressure drag:',pres_drag
        ! print*,' ** lift on airfoil:',lift
        ! print*,' **   friction lift:',fric_lift
        ! print*,' **   pressure lift:',pres_lift
        ! !
        ! open(18,file='airfoil.dat')
        ! do jedge=1,size(pso%edge)
        !   write(18,"(7(1X,E15.7E3))")pso%edge(jedge)%cen(1),     &
        !                              pso%edge(jedge)%cen(2),     &
        !                              pso%edge(jedge)%normdir(1), &
        !                              pso%edge(jedge)%normdir(2), &
        !                              prs(jedge),vel(jedge,1),  &
        !                              vel(jedge,2)
        ! enddo
        ! close(18)
        ! print*,' << airfoil.dat'
        !
        deallocate(vel,prs,rho,tmp)
        !
        counter=counter+1
        !
      endif
      !
      deallocate( qib,dv )
      !
    enddo
    !
  end subroutine ibforce
  !+-------------------------------------------------------------------+
  !| The end of the subroutine ibforce.                                |
  !+-------------------------------------------------------------------+
  !!
  !+-------------------------------------------------------------------+
  !| This function is to intepolate flow field from a cell             |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 06-Jul-2022: Created by J. Fang @ Warrington.                     |
  !+-------------------------------------------------------------------+
  subroutine interp_cell(xn,ic,qint,dvint)
    !
    use commvar,   only : ndims,numq
    use commarray, only : x,q,dvel
    use commfunc,  only : matinv
    !
    real(8),intent(in) :: xn(3)
    integer,intent(in) :: ic(3)
    real(8),intent(out),optional :: qint(numq),dvint(1:ndims,1:ndims)
    !
    real(8) :: qcell(4,numq),dvcell(4,1:ndims,1:ndims)
    real(8) :: Tm1(4,4),xcell(4,3),coef_dirichlet(4),xvec(4),Ti1(4,4)
    integer :: i,j,k,m,jq
    !
    if(ndims==2) then
      !
      do m=1,4
        !
        i=icell(ic(1),m)
        j=jcell(ic(2),m)
        k=ic(3)
        !
        xcell(m,:)=x(i,j,k,:)
        !
        Tm1(m,1)=xcell(m,1)*xcell(m,2)
        Tm1(m,2)=xcell(m,1)
        Tm1(m,3)=xcell(m,2)
        Tm1(m,4)=1.d0
        !
        if(present(qint)) qcell(m,:)=q(i,j,k,:)
        !
        if(present(dvint)) dvcell(m,1:ndims,1:ndims)=dvel(i,j,k,1:ndims,1:ndims)
      enddo
      !
      Ti1=matinv(Tm1,4)
      !
      xvec(1)=xn(1)*xn(2)
      xvec(2)=xn(1)
      xvec(3)=xn(2)
      xvec(4)=1.d0
      !
      do m=1,4
        coef_dirichlet(m)=dot_product(xvec,Ti1(:,m))
      enddo
      !
      if(present(qint)) then
        do jq=1,numq
          qint(jq)=dot_product(coef_dirichlet,qcell(:,jq))
        enddo
      endif
      !
      if(present(dvint)) then
        do j=1,ndims
        do i=1,ndims
          dvint(i,j)=dot_product(coef_dirichlet,dvcell(:,i,j))
        enddo
        enddo
      endif
      !
    elseif(ndims==3) then
      stop ' !! not yet @ interp_cell'
    endif
    !
    contains
    !
    integer function icell(i,m)
      !
      integer,intent(in) :: i,m
      !
      if(m==1) then
        icell=i-1
      elseif(m==2) then
        icell=i
      elseif(m==3) then
        icell=i
      elseif(m==4) then
        icell=i-1
      endif
      !
    end function icell
    !
    integer function jcell(j,m)
      !
      integer,intent(in) :: j,m
      !
      if(m==1) then
        jcell=j-1
      elseif(m==2) then
        jcell=j-1
      elseif(m==3) then
        jcell=j
      elseif(m==4) then
        jcell=j
      endif
      !
    end function jcell
    !
  end subroutine interp_cell
  !+-------------------------------------------------------------------+
  !| The end of the subroutine interp_cell.                            |
  !+-------------------------------------------------------------------+
  !
end module ibmethod
!+---------------------------------------------------------------------+
!| The end of the module ibmethod                                      |
!+---------------------------------------------------------------------+