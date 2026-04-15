module pastr_grid

    use iso_fortran_env, only: wp => real64
    use pastr_constdef

    implicit none

    real(wp) :: Pr,gamma,Me,Re,f1inf,g0inf,f00,f10,g00,h,ref_t
    integer :: n

contains
    
    subroutine gridgen

      use pastr_io,     only: parse_command_line

      character(len=32) :: gtype

      call parse_command_line(  string=gtype )

      if(trim(gtype)=='bl') then
        print*,' ** to generate a boundary layer grid'
        call grid_bl()
      else
        stop '  !! gtype not defined !! @ gridgen'
      endif

    end subroutine gridgen

    subroutine grid_bl

      use pastr_io,     only: parse_command_line
      use pastr_gradients, only: grad_y
      use pastr_h5io
      use pastr_tecio

      real(wp) :: bl_thickness,xmax,ymax,zmax,xlen,xlen1,lvis,dxmin
      real(wp),allocatable :: x1(:),y1(:),z1(:),dx(:),dy(:)
      real(wp),allocatable :: x(:,:,:),y(:,:,:),z(:,:,:)
      integer :: im,jm,km,in,in1,in2,jn,i,j,k
      real(wp) :: vtemp1

      call parse_command_line( inumber=im )
      call parse_command_line( inumber=jm )
      call parse_command_line( inumber=km )

      print*,' ** mesh dimenson:  ',im,'x',jm,'x',km

      open(11,file='tbl_stat.dat')
      read(11,*)lvis
      read(11,*)bl_thickness
      close(11)
      print*,' >> tbl_stat.dat'

      print*,' ** TBL thickness:   ',bl_thickness
      print*,' ** viscous length:  ',1.d0/lvis

      xlen=bl_thickness*40.d0
      xmax=50.d0
      ymax=real(int(bl_thickness*10.d0))
      zmax=real(int(bl_thickness*3.d0))

      print*,' ** xmax:',xmax
      print*,' ** ymax:',ymax
      print*,' ** zmax:',zmax

      allocate(x1(0:im),y1(0:jm),z1(0:km),dx(0:im),dy(0:jm))

      in1=im/10
      xlen1=bl_thickness*15.d0
      do i=0,in1
        x1(i)=xlen1/dble(in1)*dble(i)
      enddo
      print*,' ** Delta_x0+=',(x1(in1)-x1(in1-1))*lvis,'x',x1(in1)
      vtemp1=x1(in1)

      dxmin=7.5d0/lvis
      print*,dxmin

      in2=in1+50

      x1(in1:in2)=halfcosmesh(jn=in2-in1,dymin=dxmin,dymax=x1(in1)-x1(in1-1),mode='>')+vtemp1
      print*,' ** Delta_x+=',(x1(in2)-x1(in2-1))*lvis,'x:',x1(in2)

      in=im-50
      vtemp1=x1(in2)-x1(in2-1)
      do i=in2,in
        x1(i)=x1(i-1)+vtemp1
      enddo
      print*,' ** Delta_x+=',(x1(in)-x1(in-1))*lvis,'x:',x1(in)

      call spongstretch(im-in,xmax,x1(in-2:im))

      dx=grad_y(x1)
      open(18,file='gridx.dat')
      do i=0,im
        write(18,*)i,x1(i),dx(i)
      end do
      close(18)
      print*,' << gridx.dat ... done.'

      y1=quatcosmesh(jn=jm,dymin=0.7d0/lvis,ymax=ymax)
      print*,' ** Delta_y1+=',y1(1)*lvis
      do j=1,jm
        if(y1(j)>=bl_thickness) then
          print*,' ** Delta_yd+=',(y1(j)-y1(j-1))*lvis
          exit
        endif
      enddo
      print*,' ** Delta_ymax+=',(y1(jm)-y1(jm-1))*lvis

      dy=grad_y(y1)
      open(18,file='gridy.dat')
      do j=0,jm
        write(18,*)j,y1(j),dy(j)
      end do
      close(18)
      print*,' << gridy.dat ... done.'

      if(km==0) then
        z1=0.d0
      else
        do k=0,km
          z1(k)=zmax/dble(km)*dble(k)
        enddo
        print*,' ** Delta_z+=',z1(1)*lvis
      endif

      allocate(x(0:im,0:jm,0:km),y(0:im,0:jm,0:km),z(0:im,0:jm,0:km))
      do k=0,km
      do j=0,jm
      do i=0,im
        x(i,j,k)=x1(i)
        y(i,j,k)=y1(j)
        z(i,j,k)=z1(k)
      enddo
      enddo
      enddo

      call tecbin('tecgrid_2d.plt',x(:,:,0),'x',y(:,:,0),'y')

      call h5_writearray3d(x,im,jm,km,'x','grid.h5')
      call h5_writearray3d(y,im,jm,km,'y','grid.h5')
      call h5_writearray3d(z,im,jm,km,'z','grid.h5')

    end subroutine grid_bl

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! This subroutine is used to stretch grid for spong layers
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine spongstretch(dim,varmax,var)
      !
      integer :: dim,hh
      real(wp) :: var(-2:dim)
      real(wp) :: varmax
      real(wp) :: err,err2,dra,ddra,var1,var2,dde,ratio
      !
      dra=1._wp
      !
      ddra=0.001_wp*dra
      err=1.d10
      err2=0._wp
      !
      do while( abs(err)>1.d-9 )
        !
        ! ratio=1._wp
        ratio=(var(0)-var(-1))/(var(-1)-var(-2))
        do hh=1,dim
          !
          dde=var(hh-1)-var(hh-2)
          var(hh)=var(hh-1)+dde*ratio
          ratio=ratio*dra
          !
        end do
        !
        err=varmax-var(dim)
        !
        if(err>0) then
          dra=dra+ddra
        else
          dra=dra-ddra
        end if
        !
        if(err*err2<0._wp) ddra=0.5_wp*ddra
        !
        err2=err
        !
        !print*,err,ddra,dra,ratio
        !
      end do
      !
      print*,' the mesh is strenched to x=',varmax,' from x=',var(0),' the max ratio is:',ratio,' using ',dim,' nodes'
      !
    end subroutine spongstretch
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! End of the subroutine spongstretch
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    function quatcosmesh(jn,dymin,dymax,ymax) result(y)
      !
      integer,intent(in) :: jn
      real(wp),intent(in),optional :: dymin,dymax,ymax
      real(wp) :: y(0:jn)
      !
      integer :: j
      real(wp) :: var1,var2,dymaxtemp,ymaxtemp,error,ddy,er
      real(wp) :: dymintemp,dymin2
      !
      if(present(dymax) .and. present(dymin)) then
        !
        y(0)=0._wp
        do j=1,jn
          var1=0.5_wp*pi/dble(jn-1)*dble(j-1)
          var1=0.5_wp*(dymax-dymin)*(-cos(var1)+1._wp)+dymin
          y(j)=y(j-1)+var1
        enddo
        !
        print*,' ** dymin=',dymin,'dymax=',dymax,'ymax=',y(jn)
        !
      elseif(present(ymax) .and. present(dymin)) then
        ! first guess dymaxtemp
        ddy=0.5_wp*dymin
        dymaxtemp=dymin+ddy
        !
        error=1.d10
        do while(abs(error)>1.d-12)
          !
          ymaxtemp=quatcosmesh_ymax(dymaxtemp)
          !
          if(ymaxtemp<ymax) then
            ! dymax is too small
            dymaxtemp=dymaxtemp+ddy
          else
            dymaxtemp=dymaxtemp-ddy
          endif
          ! !
          er=(ymaxtemp-ymax)
          !
          if(er*error<0) ddy=0.5_wp*ddy
          !
          error=er
          !
          ! print*,error,ymaxtemp,dymaxtemp
          !
        enddo
        !
        y(0)=0._wp
        do j=1,jn
          var1=0.5_wp*pi/dble(jn-1)*dble(j-1)
          var1=0.5_wp*(dymaxtemp-dymin)*(-cos(var1)+1._wp)+dymin
          y(j)=y(j-1)+var1
        enddo
        !
        print*,' ** dymin=',dymin,'dymax=',dymaxtemp,'ymax=',y(jn)
        !
      elseif(present(ymax) .and. present(dymax)) then
          !
          if( dymax*dble(jn) < ymax) then
            !
            print*,'dymax=',dymax,'jn=',jn,'dymax*jn=',dymax*dble(jn),'ymax=',ymax
            stop ' !! error dymax too small' 
            !
          endif
          !
          var2=0._wp
          do j=1,jn
            var1=0.5_wp*pi/dble(jn-1)*dble(j-1)
            var1=0.5_wp*(dymax-0._wp)*(cos(var1)+1._wp)+0._wp
            var2=var2+var1
            print*,j,var1
          enddo
          !
          if(var2<=ymax) then
            print*,' minimal ymax can reached:',var2,' but larger than ymax required:',ymax
            stop ' !! increase the number ofnodes' 
          endif
          !
          print*, ' ymax range:',dymax*dble(jn),'~',var2
          print*, ' ymax   set:',ymax,'jn:',jn,'dymax:',dymax
          stop
          !
          dymintemp=0.5_wp*dymax
          error=10._wp
          !
          do while(abs(error)>1.d-8)
            var2=0._wp
            ddy=0.5_wp*dymintemp
            !
            do j=1,jn
              var1=0.5_wp*pi/dble(jn-1)*dble(j-1)
              var1=0.5_wp*(dymax-dymintemp)*(cos(var1)+1._wp)+dymintemp
              var2=var2+var1
            enddo
            !
            stop
            !
            if(error*(var2-ymax)<0) ddy=0.5_wp*ddy
            !
            if(error>0) then
                dymintemp=dymintemp-ddy
            else
                dymintemp=dymintemp+ddy
            endif
            !
            error=var2-ymax
            !
            print*,error,dymintemp,var2/ymax
            !
          enddo
          !
          if( dymintemp <= 0._wp) then
            !
            print*,dymintemp
            stop ' !! error negative dymin' 
            !
          endif
          !
          y(0)=0._wp
          do j=1,jn
            var1=0.5_wp*pi/dble(jn-1)*dble(j-1)
            var1=0.5_wp*(dymax-dymintemp)*(cos(var1)+1._wp)+dymintemp
            y(j)=y(j-1)+var1
          enddo
          !
          print*,' ** dymin=',dymintemp,'dymax=',dymax,'ymax=',y(jn)
          !
      endif
      !
      contains
      !
      function quatcosmesh_ymax(dymax1) result(ymax1)
        !
        real(wp),intent(in) :: dymax1
        real(wp) :: ymax1
        !
        real(wp) :: varc1
        integer :: j
        !
        ymax1=0._wp
        do j=1,jn
          varc1=0.5_wp*pi/dble(jn-1)*dble(j-1)
          varc1=0.5_wp*(dymax1-dymin)*(-cos(varc1)+1._wp)+dymin
          ymax1=ymax1+varc1
        enddo
        !
        return
        !
      end function quatcosmesh_ymax
      !
    end function quatcosmesh

    function halfcosmesh(jn,dymin,dymax,ymax,mode) result(y)
      !
      use pastr_fsolver,only: fxsolver
      !
      integer,intent(in) :: jn
      real(wp),intent(in),optional :: dymax,dymin,ymax
      character(len=*),intent(in) :: mode
      real(wp) :: y(0:jn),power
      !
      integer :: j
      real(wp) :: var1,var2,dymaxtemp,ymaxtemp,dymintemp,error,ddy,er
      !
      if(present(dymax) .and. present(dymin) .and. present(ymax)) then
        !
        var1=dymax*dble(jn)
        if(var1<ymax) then
          print*,' !! dymax*jn=',var1,'< ymax:',ymax
          print*,' Even with max mesh resolution required, still cannot fill the length'
          print*,' The least jn suggest here is:',int(ymax/dymax)+1
          print*,' !! halfcosmesh quite !!'
          stop
        endif
        !
        y(0)=0._wp
        !
        if(mode=='<') then
          power=fxsolver(xmin=0.1_wp,xmax=2.5_wp,ftarget=ymax,fx=getymax_dec)
          do j=1,jn
            var1=pi*(dble(j-1)/dble(jn-1))**power
            var1=0.5_wp*(dymax-dymin)*(-cos(var1)+1._wp)+dymin
            y(j)=y(j-1)+var1
          enddo
        elseif(mode=='>') then
          power=fxsolver(xmin=0.1_wp,xmax=2.5_wp,ftarget=ymax,fx=getymax_inc)
          do j=1,jn
            var1=pi*(dble(j-1)/dble(jn-1))**power
            var1=0.5_wp*(dymax-dymin)*(cos(var1)+1._wp)+dymin
            y(j)=y(j-1)+var1
          enddo
        else
          stop ' mode error @ function halfcosmesh'
        endif
        !
      elseif(present(dymax) .and. present(dymin)) then
        !
        if(mode=='<') then
          y(0)=0._wp
          do j=1,jn
            var1=pi/dble(jn-1)*dble(j-1)
            var1=0.5_wp*(dymax-dymin)*(-cos(var1)+1._wp)+dymin
            y(j)=y(j-1)+var1
          enddo
        elseif(mode=='>') then
          y(0)=0._wp
          do j=1,jn
            var1=pi/dble(jn-1)*dble(j-1)
            var1=0.5_wp*(dymax-dymin)*(cos(var1)+1._wp)+dymin
            y(j)=y(j-1)+var1
          enddo
        else
          stop ' mode error @ function halfcosmesh'
        endif
        ! 
      elseif(present(ymax) .and. present(dymin)) then
        ! first guess dymaxtemp
        ddy=0.5_wp*dymin
        dymaxtemp=dymin+ddy
        !
        error=1.d10
        do while(abs(error)>1.d-12)
          !
          ymaxtemp=halfcosmesh_ymax(dymaxtemp)
          !
          if(ymaxtemp<ymax) then
            ! dymax is too small
            dymaxtemp=dymaxtemp+ddy
          else
            dymaxtemp=dymaxtemp-ddy
          endif
          ! !
          er=(ymaxtemp-ymax)
          !
          if(er*error<0) ddy=0.5_wp*ddy
          !
          error=er
          !
          ! print*,error,ymaxtemp,dymaxtemp
          !
        enddo
        !
        y(0)=0._wp
        do j=1,jn
          var1=1._wp*pi/dble(jn-1)*dble(j-1)
          var1=0.5_wp*(dymaxtemp-dymin)*(-cos(var1)+1._wp)+dymin
          y(j)=y(j-1)+var1
        enddo
        !
      elseif(present(ymax) .and. present(dymax)) then
          !
          if( dymax*dble(jn) < ymax) then
            !
            print*,'dymax=',dymax,'jn=',jn,'dymax*jn=',dymax*dble(jn),'ymax=',ymax
            stop ' !! error dymax too small' 
            !
          endif
          !
          var2=0._wp
          do j=1,jn
            var1=pi/dble(jn-1)*dble(j-1)
            var1=0.5_wp*(dymax-0._wp)*(cos(var1)+1._wp)+0._wp
            var2=var2+var1
            ! print*,j,var1
          enddo
          !
          if(var2>=ymax) then
            print*,' minimal ymax can reached:',var2,' but larger than ymax required:',ymax
            stop ' !! reduce the number ofnodes' 
          endif
          !
          print*, ' ymax range:',dymax*dble(jn),'~',var2
          print*, ' ymax   set:',ymax,'jn:',jn,'dymax:',dymax
          !
          dymintemp=0.5_wp*dymax
          error=10._wp
          !
          ddy=0.5_wp*dymintemp
          !
          do while(.true.)
            !
            var2=0._wp
            !
            do j=1,jn
              var1=pi/dble(jn-1)*dble(j-1)
              var1=0.5_wp*(dymax-dymintemp)*(cos(var1)+1._wp)+dymintemp
              var2=var2+var1
            enddo
            !
            ! print*,error,ddy,dymintemp,var2,ymax
            !
            if(abs(error)<1.d-12) exit
            !
            if(error*(var2-ymax)<0) ddy=0.5_wp*ddy
            !
            if(error>0) then
                dymintemp=dymintemp-ddy
            else
                dymintemp=dymintemp+ddy
            endif
            !
            error=(var2-ymax)
            !
            ! stop
            !
          enddo
          !
          if( dymintemp <= 0._wp) then
            !
            print*,dymintemp
            stop ' !! error negative dymin' 
            !
          endif
          !
          y(0)=0._wp
          do j=1,jn
            var1=pi/dble(jn-1)*dble(j-1)
            var1=0.5_wp*(dymax-dymintemp)*(cos(var1)+1._wp)+dymintemp
            y(j)=y(j-1)+var1
          enddo
          !
          print*,' ** dymin=',dymintemp,'dymax=',dymax,'ymax=',y(jn)
          !
      endif
      !
      print*
      print*,' ** mesh generated by halfcosmesh  **'
      print*,'               size of mesh:',size(y)
      print*,'    mesh reso in the middle:',y(size(y)/2)-y(size(y)/2-1)
      print*,'    mesh reso at  left  end:',y(1)-y(0)
      print*,'    mesh reso at  right end:',y(size(y)-1)-y(size(y)-2)
      print*,' ** ****************************** **'
      print*
      !
      contains
      !
      function halfcosmesh_ymax(dymax1) result(ymax1)
        !
        real(wp),intent(in) :: dymax1
        real(wp) :: ymax1
        !
        real(wp) :: varc1
        integer :: j
        !
        ymax1=0._wp
        do j=1,jn
          varc1=1._wp*pi/dble(jn-1)*dble(j-1)
          varc1=0.5_wp*(dymax1-dymin)*(-cos(varc1)+1._wp)+dymin
          ymax1=ymax1+varc1
        enddo
        !
        return
        !
      end function halfcosmesh_ymax
      !
      ! return x from 0,pi
      function x0_pi(n,nmax,np) result(xcor)
        !
        integer,intent(in) :: n,nmax
        real(wp),intent(in) :: np
        real(wp) :: xcor
        !
        xcor=2._wp*dble(n)/dble(nmax)-1._wp ! -1,1
        !
        xcor=signpower(xcor,np) ! -1,1
        !
        xcor=0.5_wp*pi*xcor+0.5_wp*pi ! 0-pi
        !
      end function x0_pi
      !
      real(wp) function getymax_dec(npower)
        !
        real(wp) :: npower
        !
        integer :: j
        real(wp) :: xco
        !
        getymax_dec=0._wp
        !
        do j=1,jn
          xco=pi*(dble(j-1)/dble(jn-1))**npower
          var1=0.5_wp*(dymax-dymin)*(-cos(xco)+1._wp)+dymin
          getymax_dec=getymax_dec+var1
          ! print*,xco,var1,getymax_dec
        enddo
        !
        !
        return
        !
      end function getymax_dec
      !
      real(wp) function getymax_inc(npower)
        !
        real(wp) :: npower
        !
        integer :: j
        real(wp) :: xco
        !
        getymax_inc=0._wp
        !
        do j=1,jn
          xco=pi*(dble(j-1)/dble(jn-1))**npower
          var1=0.5_wp*(dymax-dymin)*(cos(xco)+1._wp)+dymin
          getymax_inc=getymax_inc+var1
        enddo
        !
        return
        !
      end function getymax_inc
      !
    end function halfcosmesh

    function signpower(x,n) result(xp)
      !
      real(wp),intent(in) :: x,n
      real(wp) :: xp
      !
      xp = sign(1._wp,x)*(abs(x)**n)
      !
      return
      !
    end function signpower

end module pastr_grid