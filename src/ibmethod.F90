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
end module ibmethod
!+---------------------------------------------------------------------+
!| The end of the module ibmethod                                      |
!+---------------------------------------------------------------------+