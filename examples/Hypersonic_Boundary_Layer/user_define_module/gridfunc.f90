
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! This function is the the inverse function of the function tanh(x).
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function argtanh(var)
    !
    real(8) :: var,argtanh
    !
    argtanh=0.5d0*log((1.d0+var)/(1.d0-var))
    !
    return
    !
  end function
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! End of the function argtanh.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine hyperb(ymin,ymax,y,jn)
    !
    integer,intent(in) :: jn
    real(8),intent(in) :: ymin,ymax
    !
    real(8) :: y(0:jn)
    !
    real(8) :: var1,var2,varc,vmin,dvc,delta,error,err
    integer :: j
    !
    dvc=1.e-4
    delta=0.01d0*dvc
    !
    error=1.e10
    do while(abs(error)>1.d-5)
      varc=1.d0+dvc
      !
      var1=argtanh(1.d0/varc)
      var2=1.d0/real(jn,8)*1.d0-1.d0
      !
      vmin=ymax*(1.d0+varc*tanh(var1*var2))
      !
      err=vmin-ymin
      !
      if(err>0.d0) then
        dvc=dvc-delta
      else
        dvc=dvc+delta
      endif
      !
      if(error*err<0.d0) delta=delta*0.5d0
      error=err/ymin
      !
      ! print*,vmin,ymax,error,err,varc
    enddo
    ! stop
    !
    do j=0,jn
      !
      var1=argtanh(1.d0/varc)
      var2=1.d0*(j)/(jn)*1.d0-1.d0
      !
      y(j)=ymax*(1.d0+varc*tanh(var1*var2))
      !
    end do
    !
    ! write(*,'(A)')'  |              y0              y1              ym|'
    ! write(*,'(A)')'  +------------------------------------------------+'
    ! write(*,'(A,3(1X,E15.7E3),A)')'  |',y(0),y(1),y(jm),'|'
    !
    print*,' ** the mesh is hyperbolic strenched from ',ymin,' to ',ymax,' with ',jn,' nodes'
    !
  end subroutine hyperb

function cosmesh(jn,dymin,dymax,ymax) result(y)
    !
    integer,intent(in) :: jn
    real(8),intent(in),optional :: dymin,dymax,ymax
    real(8) :: y(0:jn)
    !
    integer :: j
    real(8) :: var1,dymaxtemp,ymaxtemp,error,ddy,er
    !
    if(present(dymax) .and. present(dymin)) then
      y(0)=0.d0
      do j=1,jn
        var1=2.d0*pi/dble(jn-1)*dble(j-1)
        var1=0.5d0*(dymax-dymin)*(-cos(var1)+1.d0)+dymin
        y(j)=y(j-1)+var1
      enddo
      !
      print*,' ** dymin=',dymin,'dymax=',dymax,'ymax=',y(jn)
      !
    elseif(present(ymax) .and. present(dymin)) then
      ! first guess dymaxtemp
      ddy=0.5d0*dymin
      dymaxtemp=dymin+ddy
      !
      error=1.d10
      do while(abs(error)>1.d-12)
        !
        ymaxtemp=cosmesh_ymax(dymaxtemp)
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
        if(er*error<0) ddy=0.5d0*ddy
        !
        error=er
        !
        ! print*,error,ymaxtemp,dymaxtemp
        !
      enddo
      !
      y(0)=0.d0
      do j=1,jn
        var1=2.d0*pi/dble(jn-1)*dble(j-1)
        var1=0.5d0*(dymaxtemp-dymin)*(-cos(var1)+1.d0)+dymin
        y(j)=y(j-1)+var1
      enddo
      !
      print*,' ** dymin=',dymin,'dymax=',dymaxtemp,'ymax=',y(jn)
      !
    endif
    !
    contains
    !
    function cosmesh_ymax(dymax1) result(ymax1)
      !
      real(8),intent(in) :: dymax1
      real(8) :: ymax1
      !
      real(8) :: varc1
      integer :: j
      !
      ymax1=0.d0
      do j=1,jn
        varc1=2.d0*pi/dble(jn-1)*dble(j-1)
        varc1=0.5d0*(dymax1-dymin)*(-cos(varc1)+1.d0)+dymin
        ymax1=ymax1+varc1
      enddo
      !
      return
      !
    end function cosmesh_ymax
    !
  end function cosmesh

  function quatcosmesh(jn,dymin,dymax,ymax) result(y)
    !
    integer,intent(in) :: jn
    real(8),intent(in),optional :: dymin,dymax,ymax
    real(8) :: y(0:jn)
    !
    integer :: j
    real(8) :: var1,var2,dymaxtemp,ymaxtemp,error,ddy,er
    real(8) :: dymintemp,dymin2
    !
    if(present(dymax) .and. present(dymin)) then
      !
      y(0)=0.d0
      do j=1,jn
        var1=0.5d0*pi/dble(jn-1)*dble(j-1)
        var1=0.5d0*(dymax-dymin)*(-cos(var1)+1.d0)+dymin
        y(j)=y(j-1)+var1
      enddo
      !
      print*,' ** dymin=',dymin,'dymax=',dymax,'ymax=',y(jn)
      !
    elseif(present(ymax) .and. present(dymin)) then
      ! first guess dymaxtemp
      ddy=0.5d0*dymin
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
        if(er*error<0) ddy=0.5d0*ddy
        !
        error=er
        !
        ! print*,error,ymaxtemp,dymaxtemp
        !
      enddo
      !
      y(0)=0.d0
      do j=1,jn
        var1=0.5d0*pi/dble(jn-1)*dble(j-1)
        var1=0.5d0*(dymaxtemp-dymin)*(-cos(var1)+1.d0)+dymin
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
        var2=0.d0
        do j=1,jn
          var1=0.5d0*pi/dble(jn-1)*dble(j-1)
          var1=0.5d0*(dymax-0.d0)*(cos(var1)+1.d0)+0.d0
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
        dymintemp=0.5d0*dymax
        error=10.d0
        !
        do while(abs(error)>1.d-8)
          var2=0.d0
          ddy=0.5d0*dymintemp
          !
          do j=1,jn
            var1=0.5d0*pi/dble(jn-1)*dble(j-1)
            var1=0.5d0*(dymax-dymintemp)*(cos(var1)+1.d0)+dymintemp
            var2=var2+var1
          enddo
          !
          stop
          !
          if(error*(var2-ymax)<0) ddy=0.5d0*ddy
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
        if( dymintemp <= 0.d0) then
          !
          print*,dymintemp
          stop ' !! error negative dymin' 
          !
        endif
        !
        y(0)=0.d0
        do j=1,jn
          var1=0.5d0*pi/dble(jn-1)*dble(j-1)
          var1=0.5d0*(dymax-dymintemp)*(cos(var1)+1.d0)+dymintemp
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
      real(8),intent(in) :: dymax1
      real(8) :: ymax1
      !
      real(8) :: varc1
      integer :: j
      !
      ymax1=0.d0
      do j=1,jn
        varc1=0.5d0*pi/dble(jn-1)*dble(j-1)
        varc1=0.5d0*(dymax1-dymin)*(-cos(varc1)+1.d0)+dymin
        ymax1=ymax1+varc1
      enddo
      !
      return
      !
    end function quatcosmesh_ymax
    !
  end function quatcosmesh

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! This subroutine is used to stretch grid for spong layers
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine spongstretch_2(dim,varmax,var)
    !
    integer :: dim,hh
    real(8) :: var(-2:dim)
    real(8) :: varmax
    real(8) :: err,err2,dra,ddra,var1,var2,dde,ratio
    !
    dra=1.d0
    !
    ddra=0.001d0*dra
    err=1.d10
    err2=0.d0
    !
    do while( dabs(err)>1.d-9 )
      !
      ! ratio=1.d0
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
      if(err*err2<0.d0) ddra=0.5d0*ddra
      !
      err2=err
      !
      !print*,err,ddra,dra,ratio
      !
    end do
    !
    ! print*,' the mesh is strenched to x=',varmax,' from x=',var(0),' the max ratio is:',ratio,' using ',dim,' nodes'
    !
  end subroutine spongstretch_2
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! End of the subroutine spongstretch_2
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!