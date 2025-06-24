!+---------------------------------------------------------------------+
!| This module contains routines related to calculate numerical flux   |
!| calculations.                                                       |
!+---------------------------------------------------------------------+
!| CHANGE RECORD                                                       |
!| -------------                                                       |
!| 03-06-2025  | Created by J. Fang @ Liverpool                        |
!+---------------------------------------------------------------------+
module flux

    use parallel, only : mpirank,mpistop

    implicit none

    real(8),allocatable,dimension(:,:),target :: cfluxi_uwd,cfluxj_uwd,cfluxk_uwd, & !cflux: compact flux calculation coefficient
                                                 cfluxi_dwd,cfluxj_dwd,cfluxk_dwd
    integer :: ifxs,ifxe,jfxs,jfxe,kfxs,kfxe ! starting and ending nodes of finite difference schemes

    contains

    !+-------------------------------------------------------------------+
    !| This subroutine is to initiate the comapct flux scheme.           |
    !+-------------------------------------------------------------------+
    !| CHANGE RECORD                                                     |
    !| -------------                                                     |
    !| 03-06-2025  | Rewrite by J. Fang @ Liverpool.                     |
    !+-------------------------------------------------------------------+
    subroutine compact_flux_initiate(scheme,ntype,dim,dir)

        use constdef
        use commfunc,only : tridiagonal_thomas_proprocess

        integer,intent(in) :: scheme,ntype,dim,dir

        integer :: i_0,i_m,j
        real(8),allocatable :: a_uwd(:),c_uwd(:),a_dwd(:),c_dwd(:)

        select case (ntype)
        case (1)
          ! the block with boundary at i==0
          i_0=-1
          i_m=dim+1
        case (2)
          ! the block with boundary at i==im
          i_0=-2
          i_m=dim
        case (3)
          ! inner block
          i_0=-2
          i_m=dim+1
        case (4)
          ! the block with boundary at i=0 and i=im
          i_0=-1
          i_m=dim
        case default
          print*, ' !! error 1 in subroutine compact_filter_initiate !'
        end select

        allocate(a_uwd(i_0:i_m),c_uwd(i_0:i_m),a_dwd(i_0:i_m),c_dwd(i_0:i_m))

        if(scheme/100==5) then

          a_uwd(:)=0.5d0;  c_uwd(:)=num1d6
          a_dwd(:)=num1d6; c_dwd(:)=0.5d0

          ! default setup, interface
          a_uwd(i_0)=0.d0;     c_uwd(i_0)=0.d0 ! explicit central scheme at interface
          a_uwd(i_m)=0.d0;     c_uwd(i_m)=0.d0 ! explicit central scheme at interface

          a_dwd(i_0)=0.d0;     c_dwd(i_0)=0.d0 ! explicit central scheme at interface
          a_dwd(i_m)=0.d0;     c_dwd(i_m)=0.d0 ! explicit central scheme at interface

          if(scheme==543) then

            ! set near-boundary/interface schemes
            if(ntype==1 .or. ntype==4) then
              a_uwd(i_0)   =2.d0;   c_uwd(i_0)  =2.d0     ! 3nd-order downwind-biased scheme
              a_uwd(i_0+1) =0.25d0; c_uwd(i_0+1)=0.25d0   ! 4th-order central scheme
              a_dwd(i_0)   =2.d0;   c_dwd(i_0)  =2.d0     ! 3nd-order upwind-biased scheme
              a_dwd(i_0+1) =0.25d0; c_dwd(i_0+1)=0.25d0   ! 4th-order central scheme
            endif
    
            if(ntype==2 .or. ntype==4) then
              a_uwd(i_m)   =2.d0;   c_uwd(i_m)  =2.d0    ! 3nd-order upwind-biased scheme
              a_uwd(i_m-1) =0.25d0; c_uwd(i_m-1)=0.25d0  ! 4th-order central scheme
              a_dwd(i_m)   =2.d0;   c_dwd(i_m)  =2.d0    ! 3nd-order down-biased scheme
              a_dwd(i_m-1) =0.25d0; c_dwd(i_m-1)=0.25d0  ! 4th-order central scheme
            endif
          else
            print*,' scheme:',scheme
            stop ' !! scheme not defined @ coef_diffcompac'
          endif

        else
          stop ' !! scheme not defined @ coef_diffcompac'
        endif

        if(dir==1) then
            ifxs=i_0
            ifxe=i_m
            allocate(cfluxi_uwd(1:3,i_0:i_m),cfluxi_dwd(1:3,i_0:i_m))

            call tridiagonal_thomas_proprocess(a_uwd,c_uwd,cfluxi_uwd)
            call tridiagonal_thomas_proprocess(a_dwd,c_dwd,cfluxi_dwd)
        elseif(dir==2) then
            jfxs=i_0
            jfxe=i_m
            allocate(cfluxj_uwd(1:3,i_0:i_m),cfluxj_dwd(1:3,i_0:i_m))

            call tridiagonal_thomas_proprocess(a_uwd,c_uwd,cfluxj_uwd)
            call tridiagonal_thomas_proprocess(a_dwd,c_dwd,cfluxj_dwd)
        elseif(dir==3) then
            kfxs=i_0
            kfxe=i_m
            allocate(cfluxk_uwd(1:3,i_0:i_m),cfluxk_dwd(1:3,i_0:i_m))

            call tridiagonal_thomas_proprocess(a_uwd,c_uwd,cfluxk_uwd)
            call tridiagonal_thomas_proprocess(a_dwd,c_dwd,cfluxk_dwd)
        else
            stop ' !! error 2 in subroutine compact_filter_initiate !'
        endif

    end subroutine compact_flux_initiate
    !+-------------------------------------------------------------------+
    !| The end of the subroutine compact_flux_initiate.                  |
    !+-------------------------------------------------------------------+

    !+-------------------------------------------------------------------+
    !| This function is to calculate the interface flux of a input array | 
    !| using compact scheme                                              |
    !+-------------------------------------------------------------------+
    !| CHANGE RECORD                                                     |
    !| -------------                                                     |
    !| 03-06-2025  | Rewrite by J. Fang @ Liverpool.                     |
    !+-------------------------------------------------------------------+
    function flux_compact(f,ntype,dim,dir,wind) result(fh)

        use commvar, only : hm
        use commfunc,only : tridiagonal_thomas_solver

        ! arguments
        integer,intent(in) :: ntype,dim,dir
        real(8),intent(in) :: f(-hm:dim+hm)
        character(len=1),intent(in) :: wind
        real(8) :: fh(-1:dim)

        ! local data
        integer :: l,i_0,i_m
        real(8) :: alpha,var0,var1,var2,var3,var4,var5
        real(8),pointer :: ac(:,:)
        real(8),allocatable :: d(:),xx(:)

        if(dir==1) then

            i_0=ifxs
            i_m=ifxe
            
            if(wind=='+') then
                ac=>cfluxi_uwd
            elseif(wind=='-') then
                ac=>cfluxi_dwd
            endif

        elseif(dir==2) then

            i_0=jfxs
            i_m=jfxe

            if(wind=='+') then
                ac=>cfluxj_uwd
            elseif(wind=='-') then
                ac=>cfluxj_dwd
            endif

        elseif(dir==3) then

            i_0=kfxs
            i_m=kfxe

            if(wind=='+') then
                ac=>cfluxk_uwd
            elseif(wind=='-') then
                ac=>cfluxk_dwd
            endif

        endif

        d=compact_flux_rhs(f,ntype,dim,dir,wind)

        allocate(xx(i_0:i_m))

        xx=tridiagonal_thomas_solver(ac,d)

        fh(-1:dim)=xx(-1:dim)

        return

    end function flux_compact
    !+-------------------------------------------------------------------+
    !| The end of the function flux_compact.                             |
    !+-------------------------------------------------------------------+

    !+-------------------------------------------------------------------+
    !| This function is to calculate R.H.S of a compact flux             |
    !| reconstruction                                                    | 
    !+-------------------------------------------------------------------+
    !| CHANGE RECORD                                                     |
    !| -------------                                                     |
    !| 03-06-2025  | Rewrite by J. Fang @ Liverpool.                     |
    !+-------------------------------------------------------------------+
    function compact_flux_rhs(f,ntype,dim,dir,wind) result(d)

        use constdef
        use commvar, only : hm

        integer,intent(in) :: ntype,dim,dir
        real(8),intent(in) :: f(-hm:dim+hm)
        character(len=1),intent(in) :: wind
        real(8),allocatable :: d(:)

        integer :: i_0,i_m,j,k,i_s,i_e
        real(8) :: var1,var2,var3

        if(dir==1) then
            i_0=ifxs
            i_m=ifxe
        elseif(dir==2) then
            i_0=jfxs
            i_m=jfxe
        elseif(dir==3) then
            i_0=kfxs
            i_m=kfxe
        endif

        allocate(d(i_0:i_m))

        if(ntype==1 .or. ntype==4) then

          ! physical boundary

          i_s=i_0+2

          ! 4th-order compact
          j=i_0+1 ! i=0
          var1=f(j)+f(j+1)
          d(j)=0.75d0*var1

          ! 3rd-order flux
          j=i_0 ! i=-1
          d(j)=2.5d0*var1*f(j+1)+0.5d0*f(j+2)

        else

          ! interfaces
          i_s=i_0+1

          ! 6th-order explicit
          j=i_0
          var1=f(j)  +f(j+1)
          var2=f(j-1)+f(j+2)
          var3=f(j-2)+f(j+3)
          d(j)=num37d60*var1 - num2d15*var2 + num1d60*var3

        endif

        if(ntype==2 .or. ntype==4) then

          ! physical boundary
          i_e=i_m-2

          ! 4th-order compact i=im-1
          j=i_m-1
          var1=f(j)+f(j+1)
          d(j)=0.75d0*var1
    
          ! 3rd-order flux i=im
          j=i_m
          d(j)=0.5d0*f(j-1) + 2.5d0*f(j)

        else
          ! interfaces
          i_e=i_m-1

          ! 6th-order explicit
          j=i_m
          var1=f(j)  +f(j+1)
          var2=f(j-1)+f(j+2)
          var3=f(j-2)+f(j+3)
          d(j)=num37d60*var1 - num2d15*var2 + num1d60*var3

        endif

        if(wind=='+') then
          ! 5th-order compact upwind
          do j=i_s,i_e
            d(j)=num1d18*f(j-1) + num19d18*f(j) + num5d9*f(j+1)
          end do
        elseif(wind=='-') then
          do j=i_s,i_e
            d(j)=num5d9*f(j) + num19d18*f(j+1) + num1d18*f(j+2)
          end do
        endif
    
    end function compact_flux_rhs
    !+-------------------------------------------------------------------+
    !| The end of the function compact_flux_rhs.                         |
    !+-------------------------------------------------------------------+

    !+-------------------------------------------------------------------+
    !| This function is to return the explicit reconstructed value of a  !
    !| input  function.                                                  |
    !+-------------------------------------------------------------------+
    !| CHANGE RECORD                                                     |
    !| -------------                                                     |
    !| 08-02-2021  | Created by J. Fang @ Warrington.                    |
    !+-------------------------------------------------------------------+
    function recons_exp(f,inode,dim,ntype,reschem,shock,solid) result(fc)
      !
      use commvar,  only: bfacmpld
      !
      ! arguments
      real(8),intent(in) :: f(1:8)
      integer,intent(in) :: inode,dim,ntype,reschem
      logical,intent(in) :: shock,solid
      real(8) :: fc
      !
      ! local data
      !
      if((ntype==1 .and. inode==0)    .or.(ntype==2 .and. inode==dim-1)) then
        !
        select case(reschem)
        case(-1)
          fc=f(4)
        case default
          fc=0.5d0*(f(4)+f(5))
        end select
        !
      elseif((ntype==1 .and. inode==1).or.(ntype==2 .and. inode==dim-2)) then
        !
        select case(reschem)
        case(-1)
          fc=f(4)
        case default
          fc=SUW3(f(3:5))
        end select
        !
      ! else
      elseif((ntype==1 .and. inode==2).or.(ntype==2 .and. inode==dim-3)) then
        !
        select case(reschem)
        case(-1)
          fc=f(4)
        case(0)
          fc=suw5(f(2:6))
        case(1)
          fc=weno5(f(2:6))
        case(2)
          fc=weno5z(f(2:6))
        case(3)
          fc=mp5(f(2:6))
        ! case(4)
        !   call weno7sym(f(1:8),var1)
        case(5)
          fc=mp5ld(f(2:7),bfacmpld,shock,solid)
        case(6) !xi
          fc=round(f(3:5))  
        case default
          print*,' !! 1 Reconstruction scheme not defined @ recons_exp !!'
          stop
        end select
        !
      else
        !
        select case(reschem)
        case(-1)
          fc=f(4)
        case(0)
          fc=suw7(f(1:7))
        case(1)
          fc=weno7(f(1:7))
        case(2)
          fc=weno7z(f(1:7))
        case(3)
          fc=mp7(f(1:7))
        ! case(4)
        !   call weno7sym(f(1:8),var1)
        case(5)
          fc=mp7ld(f(1:8),bfacmpld,shock,solid)
        case(6) !xi
          fc=round(f(3:5))
        case default
          print*,' !! 2 Reconstruction scheme not defined @ recons_exp !!'
          stop
        end select
        !
      endif
      !
      return
      !
    end function recons_exp
    !+-------------------------------------------------------------------+
    !| The end of the function recons.                                   |
    !+-------------------------------------------------------------------+

    !+-------------------------------------------------------------------+
    !| This function is to apply the MP limiter to reconstructured       |
    !| interface value                                                   |
    !+-------------------------------------------------------------------+
    !| CHANGE RECORD                                                     |
    !| -------------                                                     |
    !| 17-05-2022  | Created by J. Fang @ Warrington.                    |
    !+-------------------------------------------------------------------+
    pure function mplimiter(f,fl,shock,inode,dim,ntype)  result(fc)
      !
      real(8),intent(in) :: f(1:5),fl
      integer,intent(in) :: inode,dim,ntype
      logical,intent(in) :: shock
      real(8) :: fc
      !
      if((ntype==1 .and. inode==0) .or. (ntype==2 .and. inode==dim-1) .or. &
         (ntype==1 .and. inode==1) .or. (ntype==2 .and. inode==dim-2) ) then
        fc=fl
      else
        fc=mp5(f,fl,discont=shock)
      endif
      !
      return
      !
    end function mplimiter
    !+-------------------------------------------------------------------+
    !| The end of the function mplimiter.                                |
    !+-------------------------------------------------------------------+

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! This subroutine offers the interface's values with standard  upwind 
    ! scheme.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    pure function suw3(u) result(uh)
      
      use constdef

      real(8),intent(in) :: u(1:3)
      real(8) :: uh
      !
      uh=-num1d6*u(1)+num5d6*u(2)+num1d3*u(3) 
      !
      return
      !
    end function suw3
    !
    pure function suw5(u) result(uh)
      !
      real(8),intent(in) :: u(1:5)
      real(8) :: uh
      !
      uh=(2.d0*u(1)-13.d0*u(2)+47.d0*u(3)+27.d0*u(4)-3.d0*u(5))/60.d0
      !
      return
      !
    end function suw5
    !
    pure function suw7(u) result(uh)
      !
      real(8),intent(in) :: u(1:7)
      real(8) :: uh
      !
      uh=(-3.d0*u(1)+25.d0*u(2)-101.d0*u(3)+319.d0*u(4)+214.d0*u(5)-     &
          38.d0*u(6)+ 4.d0*u(7))/420.d0
      !
      return
      !
    end function suw7
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! end of the subroutine suw7
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! This subroutine is used to calculate the values of cell's interface
    ! using 5th-order Monotonicity-Preserving method.
    ! Ref. 1 : Suresh A. and Huynh H. T., Accurate Monotonicity-Preserving 
    ! Schemes with Runge¨CKutta Time Stepping, J. C. P., 1997, 136: 83¨C99
    ! .
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    pure function MP5(u,ul,discont)  result(uh)
      
      use commfunc, only : minmod2,minmod4

      real(8),intent(in) :: u(1:5)
      real(8),intent(in),optional :: ul
      logical,intent(in),optional :: discont
      real(8) :: uh
      !
      real(8) :: ulinear,uMP,uUL,uAV,uMD,uLC,uMIN,uMAX
      real(8) :: var1,var2,dm1,d0,d1,dhm1,dh0
      logical :: lskt
      !
      if(present(ul)) then
        ulinear=ul
      else
        ulinear=(2.d0*u(1)-13.d0*u(2)+47.d0*u(3)+27.d0*u(4)-3.d0*u(5))/60.d0
      endif
      !
      if(present(discont)) then
        lskt=discont
      else
        lskt=.true.
      endif
      !
      var1=u(4)-u(3)
      var2=4.d0*(u(3)-u(2))
      uMP=u(3)+minmod2(var1,var2)
      !
      var1=(ulinear-u(3))*(ulinear-uMP)
      if(lskt .and. var1>=1.d-10) then
        !
        dm1=u(1)-2.d0*u(2)+u(3)
        d0= u(2)-2.d0*u(3)+u(4)
        d1= u(3)-2.d0*u(4)+u(5)
        !
        dhm1=minmod4( 4.d0*dm1-d0,4.d0*d0-dm1,dm1,d0 )
        dh0 =minmod4(  4.d0*d0-d1, 4.d0*d1-d0, d0,d1 )
        !
        uUL=u(3)+4.d0*(u(3)-u(2))
        uAV=0.5d0*(u(3)+u(4))
        uMD=uAV-0.5d0*dh0
        uLC=u(3)+0.5d0*(u(3)-u(2))+1.333333333333333d0*dhm1
        !
        var1=min(u(3),u(4),uMD)
        var2=min(u(3),uUL,uLC)
        uMIN=max(var1,var2)
        !
        var1=max(u(3),u(4),uMD)
        var2=max(u(3),uUL,uLC)
        uMAX=min(var1,var2)
        !
        var1=uMIN-ulinear
        var2=uMAX-ulinear
        uh=ulinear+minmod2(var1,var2)
      else
        uh=ulinear
        ! No limiter is needed
      end if
      !
      return
      !
    end function MP5
    !
    pure function MP7(u,ul) result(uh)

      use commfunc, only : minmod2,minmod4

      real(8),intent(in) :: u(0:6)
      real(8),intent(in),optional :: ul
      real(8) :: uh
      !
      real(8) :: ulinear,uMP,uUL,uAV,uMD,uLC,uMIN,uMAX
      real(8) :: var1,var2,dm1,d0,d1,dhm1,dh0
      !
      if(present(ul)) then
        ulinear=ul
      else
        ulinear=(-3.d0*u(0)+25.d0*u(1)-101.d0*u(2)+319.d0*u(3)+          &
                 214.d0*u(4)-38.d0*u(5)+ 4.d0*u(6))/420.d0
      endif
      !
      var1=u(4)-u(3)
      var2=4.d0*(u(3)-u(2))
      uMP=u(3)+minmod2(var1,var2)
      !
      var1=(ulinear-u(3))*(ulinear-uMP)
      if(var1>=1.d-10) then
        dm1=u(1)-2.d0*u(2)+u(3)
        d0= u(2)-2.d0*u(3)+u(4)
        d1= u(3)-2.d0*u(4)+u(5)
        !
        dhm1=minmod4( 4.d0*dm1-d0,4.d0*d0-dm1,dm1,d0 )
        dh0 =minmod4(  4.d0*d0-d1, 4.d0*d1-d0, d0,d1 )
        !
        uUL=u(3)+4.d0*(u(3)-u(2))
        uAV=0.5d0*(u(3)+u(4))
        uMD=uAV-0.5d0*dh0
        uLC=u(3)+0.5d0*(u(3)-u(2))+1.333333333333333d0*dhm1
        !
        var1=min(u(3),u(4),uMD)
        var2=min(u(3),uUL,uLC)
        uMIN=max(var1,var2)
        !
        var1=max(u(3),u(4),uMD)
        var2=max(u(3),uUL,uLC)
        uMAX=min(var1,var2)
        !
        var1=uMIN-ulinear
        var2=uMAX-ulinear
        uh=ulinear+minmod2(var1,var2)
      else
        uh=ulinear
      end if
      !
      return
      !
    end function MP7
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! end of the function MP
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    function MP5LD(u,weightBW,lskt,lsod) result(uh)
      
      use commfunc, only : minmod2,minmod4

      real(8),intent(in) :: u(1:6),weightBW
      logical,intent(in) :: lskt,lsod
      real(8) :: uh
      !
      ! local data
      real(8) :: ulinear,uMP,uUL,uAV,uMD,uLC,uMIN,uMAX
      real(8) :: var1,var2,dm1,d0,d1,dhm1,dh0
      !
      real(8) :: b1,b2,b3,b4,b5,b6,vadp
      !
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! vadp: optimize parameter, to control the dissipation
      ! of the linear scheme. 
      ! vadp from 0 to 1/60
      ! vadp=0: standard 5 order upwind. Max dissipation
      ! vadp=1/60: standard 6 order center. 0 dissipation
      ! lskt: shock or not
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      ! if(lsod) then
      !   ulinear=-num1d6*u(2)+num5d6*u(3)+num1d3*u(4)
      !   ! ulinear=u(3)
      ! else
        vadp=1.666666666666667d-2*weightBW
        !vadp=0.0d0*1.666666666666667d-2
        !
        b1=      -vadp+3.333333333333333d-2
        b2=  5.d0*vadp-2.166666666666667d-1
        b3=-10.d0*vadp+7.833333333333333d-1
        b4= 10.d0*vadp+0.45d0
        b5= -5.d0*vadp-5.d-2
        b6=       vadp
        !
        ulinear=b1*u(1)+b2*u(2)+b3*u(3)+b4*u(4)+b5*u(5)+b6*u(6)
      ! endif
      !!
      var1=u(4)-u(3)
      var2=4.d0*(u(3)-u(2))
      uMP=u(3)+minmod2(var1,var2)
      !
      var1=(ulinear-u(3))*(ulinear-uMP)
      ! if(lskt) then
      if(var1>=1.d-10) then
        dm1=u(1)-2.d0*u(2)+u(3)
        d0= u(2)-2.d0*u(3)+u(4)
        d1= u(3)-2.d0*u(4)+u(5)
        !
        dhm1=minmod4( 4.d0*dm1-d0,4.d0*d0-dm1,dm1,d0 )
        dh0 =minmod4(  4.d0*d0-d1, 4.d0*d1-d0, d0,d1 )
        !
        uUL=u(3)+4.d0*(u(3)-u(2))
        uAV=0.5d0*(u(3)+u(4))
        uMD=uAV-0.5d0*dh0
        uLC=u(3)+0.5d0*(u(3)-u(2))+1.333333333333333d0*dhm1
        !
        var1=min(u(3),u(4),uMD)
        var2=min(u(3),uUL,uLC)
        uMIN=max(var1,var2)
        !
        var1=max(u(3),u(4),uMD)
        var2=max(u(3),uUL,uLC)
        uMAX=min(var1,var2)
        !
        var1=uMIN-ulinear
        var2=uMAX-ulinear
        uh=ulinear+minmod2(var1,var2)
      else
        uh=ulinear
        !
      end if
      !
      return
      !
    end function MP5LD
    !
    function MP7LD(u,weightBW,lskt,lsod) result(uh)
        
      use commfunc, only : minmod2,minmod4

      real(8),intent(in) :: u(0:7),weightBW
      logical,intent(in) :: lskt,lsod
      real(8) :: uh
      !
      ! local data
      real(8) :: ulinear,uMP,uUL,uAV,uMD,uLC,uMIN,uMAX
      real(8) :: var1,var2,dm1,d0,d1,dhm1,dh0
      real(8) :: b(0:7)
      real(8) :: vadp
      !
      ! if(lsod) then
      !   ulinear=-num1d6*u(2)+num5d6*u(3)+num1d3*u(4) 
      !   ! ulinear=u(3)
      ! else
        !
        vadp=3.571428571428571d-3*weightBW
        !
        b(0)=  1.d0*vadp-7.142857142857143d-3
        b(1)= -7.d0*vadp+5.952380952380952d-2
        b(2)= 21.d0*vadp-0.240476190476190d0 
        b(3)=-35.d0*vadp+0.759523809523809d0 
        b(4)= 35.d0*vadp+0.509523809523809d0 
        b(5)=-21.d0*vadp-9.047619047619047d-2
        b(6)=  7.d0*vadp+9.523809523809525d-3
        b(7)= -1.d0*vadp
        !
        ulinear= b(0)*u(0)+b(1)*u(1)+b(2)*u(2)+b(3)*u(3)+b(4)*u(4)         &
                +b(5)*u(5)+b(6)*u(6)+b(7)*u(7)
      ! endif
      !
      var1=u(4)-u(3)
      var2=4.d0*(u(3)-u(2))
      uMP=u(3)+minmod2(var1,var2)
      !
      var1=(ulinear-u(3))*(ulinear-uMP)
      if(lskt) then
        !
        dm1=u(1)-2.d0*u(2)+u(3)
        d0= u(2)-2.d0*u(3)+u(4)
        d1= u(3)-2.d0*u(4)+u(5)
        !
        dhm1=minmod4( 4.d0*dm1-d0,4.d0*d0-dm1,dm1,d0 )
        dh0 =minmod4(  4.d0*d0-d1, 4.d0*d1-d0, d0,d1 )
        !
        uUL=u(3)+4.d0*(u(3)-u(2))
        uAV=0.5d0*(u(3)+u(4))
        uMD=uAV-0.5d0*dh0
        uLC=u(3)+0.5d0*(u(3)-u(2))+1.333333333333333d0*dhm1
        !
        var1=min(u(3),u(4),uMD)
        var2=min(u(3),uUL,uLC)
        uMIN=max(var1,var2)
        !
        var1=max(u(3),u(4),uMD)
        var2=max(u(3),uUL,uLC)
        uMAX=min(var1,var2)
        !
        var1=uMIN-ulinear
        var2=uMAX-ulinear
        uh=ulinear+minmod2(var1,var2)
      else
        uh=ulinear
        ! No limiter is needed
        !
      end if
      !
      return
      !
    end function MP7LD
    !
  
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Three point ROUND schemes: Reconstruction Operator on Unified    
    ! Normalized-variable Diagram
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! added by Xi Deng, 06/02/2022
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    pure function round(u) result(uh)
      !
      real(8),intent(in) :: u(1:3)
      real(8) :: uh
      !!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! u: 3 points grid value
      ! uh: final interface value
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real(8) :: eps
      
      real(8) :: z0c,a1c,a2c,plc,prc,p1c,p2c,p3c
      
      real(8) :: wc1c,wc2c,gc
      
      eps=1.d-16  
      
      z0c=(u(2)-u(1)+eps)/(u(3)-u(1)+eps);
      
      
      a1c=1.0d0+12.0d0*z0c*z0c;
      a2c=1.0d0+5.0d0*(z0c-1.0d0)*(z0c-1.0d0);
      plc=1100.0d0*(z0c-0.05d0)*(z0c-0.05d0)*(z0c-0.05d0)*(0.47d0-z0c)*(0.47d0-z0c)*(0.47d0-z0c)
      prc=18000.0d0*(z0c-0.55d0)*(z0c-0.55d0)*(z0c-0.55d0)*(0.97d0-z0c)*(0.97d0-z0c)*(0.97d0-z0c)*(0.97d0-z0c)*(0.97d0-z0c)
      p1c=0.833333333333333d0*z0c+0.333333333333333d0+max(plc,0.0d0)+max(prc,0.0d0)
      p2c=1.5d0*z0c
      p3c=0.5d0*z0c+0.5d0
      wc1c=1.0d0/a1c/a1c/a1c/a1c
      wc2c=1.0d0/a2c/a2c/a2c/a2c/a2c/a2c/a2c/a2c
      
      gc=(p1c*(1.0d0-wc1c)+p2c*wc1c)*(1.0d0-wc2c)+p3c*wc2c
      
      uh=gc*(u(3)-u(1))+u(1)
      
      return
      
    end function round
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! End of the subroutine ROUND.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! this subroutine offers the interface's values with 5 point WENO
    ! Scheme according the points' value you give.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    pure function WENO5(u) result(uh)
      !
      real(8),intent(in) :: u(1:5)
      real(8) :: uh
      !!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! u: 5 points grid value
      ! uh: final interface value
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real(8) :: uh1,uh2,uh3,beter1,beter2,beter3,alfa1,alfa2,alfa3,     &
                 WT1,WT2,WT3,C1,C2,C3,alfaS
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! uh1: i+1/2 value for stencil 1 
      !      with 3 pints upwind shmeme
      ! uh2: i+1/2 value for stencil 2 
      !      with 3 pints upwind shmeme 
      ! uh3: i+1/2 value for stencil 3 
      !      with 3 pints upwind shmeme
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real(8) :: eps
      !
      eps=1.d-6
      !
      uh1=0.333333333333333d0*u(1)-1.16666666666667d0*u(2)+            &
          1.83333333333333d0*u(3)                                      
      uh2=-0.166666666666667d0*u(2)+0.833333333333333d0*u(3)+          &
           0.333333333333333d0*u(4)                                    
      uh3=0.333333333333333d0*u(3)+0.833333333333333d0*u(4)           &
          -0.166666666666667d0*u(5)                                    
      !                                                                
      beter1=1.08333333333333d0*(u(1)-2.d0*u(2)+u(3))**2+              &
             0.25d0*(u(1)-4.d0*u(2)+3.d0*u(3))**2                      
      beter2=1.08333333333333d0*(u(2)-2.d0*u(3)+u(4))**2+              &
             0.25d0*(u(2)-u(4))**2                                     
      beter3=1.08333333333333d0*(u(3)-2.d0*u(4)+u(5))**2+              &
             0.25d0*(3.d0*u(3)-4.d0*u(4)+u(5))**2
      !
      C1=0.1d0
      C2=0.6d0
      C3=0.3d0
      !
      alfa1=C1/(beter1+eps)**2
      alfa2=C2/(beter2+eps)**2
      alfa3=C3/(beter3+eps)**2
      !
      alfaS=alfa1+alfa2+alfa3
      !
      WT1=alfa1/alfaS
      WT2=alfa2/alfaS
      WT3=alfa3/alfaS
      !
      uh=WT1*uh1+WT2*uh2+WT3*uh3
      !
      return
      !
    end function WENO5
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! End of the subroutine WENO5.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! this subroutine offers the interface's values with 7 point WENO
    ! Scheme according the points' value you give.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    pure function WENO7(u) result(uh)
      !
      real(8),intent(in) :: u(1:7)
      real(8) :: uh
      !!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! dir=1: upwind
      ! dir=-1: backwind
      ! u: 5 points grid value
      ! uh: final interface value
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real(8) :: uh1,uh2,uh3,uh4
      real(8) :: beter1,beter2,beter3,beter4,alfa1,alfa2,alfa3,alfa4,  &
                 WT1,WT2,WT3,WT4,C1,C2,C3,C4,alfaS
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! uh1: i+1/2 value for stencil 1 
      !      with 3 pints upwind shmeme
      ! uh2: i+1/2 value for stencil 2 
      !      with 3 pints upwind shmeme 
      ! uh3: i+1/2 value for stencil 3 
      !      with 3 pints upwind shmeme
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real(8) :: df1,df2,df3
      real(8) :: eps
      !
      eps=1.d-6
      !
      uh1=-0.25d0              *u(1)+1.083333333333333d0 *u(2)       &
          -1.916666666666667d0 *u(3)+2.083333333333333d0 *u(4)
      uh2= 8.333333333333333d-2*u(2)-4.166666666666667d-1*u(3)       &
          +1.083333333333333d0 *u(4)+0.25d0              *u(5)
      uh3=-8.333333333333333d-2*u(3)+5.833333333333333d-1*u(4)       &
          +5.833333333333333d-1*u(5)-8.333333333333333d-2*u(6)
      uh4= 0.25d0              *u(4)+1.083333333333333d0 *u(5)       &
          -4.166666666666667d-1*u(6)+8.333333333333333d-2*u(7)
      !
      df1=2.777777777777778d-2*(-2.d0*u(1)+9.d0*u(2)-18.d0*u(3)+11d0*u(4))**2
      df2= 1.083333333333333d0*(-1.d0*u(1)+4.d0*u(2) -5.d0*u(3) +2d0*u(4))**2
      df3= 1.084722222222222d0*(-1.d0*u(1)+3.d0*u(2) -3.d0*u(3) +1d0*u(4))**2
      beter1=df1+df2+df3
      !
      df1=2.777777777777778d-2*(      u(2)-6.d0*u(3)+3.d0*u(4)+2.d0*u(5))**2
      df2= 1.083333333333333d0*(                u(3)-2.d0*u(4)+     u(5))**2
      df3= 1.084722222222222d0*(-1.d0*u(2)+3.d0*u(3)-3.d0*u(4)+1.d0*u(5))**2
      beter2=df1+df2+df3
      !
      df1=2.777777777777778d-2*(-2.d0*u(3)-3.d0*u(4)+6.d0*u(5)-1.d0*u(6))**2
      df2= 1.083333333333333d0*(      u(3)-2.d0*u(4)     +u(5)          )**2
      df3= 1.084722222222222d0*(-1.d0*u(3)+3.d0*u(4)-3.d0*u(5)+1.d0*u(6))**2
      beter3=df1+df2+df3
      !
      df1=2.777777777777778d-2*(-11.d0*u(4)+18.d0*u(5)-9.d0*u(6)+2.d0*u(7))**2
      df2= 1.083333333333333d0*(  2.d0*u(4) -5.d0*u(5)+4.d0*u(6)     -u(7))**2
      df3= 1.084722222222222d0*( -1.d0*u(4) +3.d0*u(5)-3.d0*u(6)+1.d0*u(7))**2
      beter4=df1+df2+df3
      !
      C1=2.857142857142857d-2
      C2=3.428571428571429d-1
      C3=5.142857142857143d-1
      C4=1.142857142857143d-1
      !
      alfa1=C1/(beter1+eps)**2
      alfa2=C2/(beter2+eps)**2
      alfa3=C3/(beter3+eps)**2
      alfa4=C4/(beter4+eps)**2
      !
      alfaS=alfa1+alfa2+alfa3+alfa4
      !
      WT1=alfa1/alfaS
      WT2=alfa2/alfaS
      WT3=alfa3/alfaS
      WT4=alfa4/alfaS
      !
      uh=WT1*uh1+WT2*uh2+WT3*uh3+WT4*uh4
      !
      return
      !
    end function WENO7
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !End of the subroutine WENO7.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! this subroutine offers the interface's values with 5 order WENOZ 
    ! scheme of Borges, 2008. that modified wight function.
    ! Ref:  R. Borges, M. Carmona, B. Costa B and et al. An improved 
    ! weighted essentially non-oscillatory scheme for hyperbolic 
    ! conservation laws, J. Comput. Phys. 227 (2008) 3191¨C3211.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    pure function WENO5Z(u) result(uh)
      !
      real(8),intent(in) :: u(1:5)
      real(8) :: uh
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! u: 5 points grid value
      ! uh: final interface value
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real(8) :: uh1,uh2,uh3,beter1,beter2,beter3,alfa1,alfa2,alfa3,     &
                 WT1,WT2,WT3,C1,C2,C3,alfaS,tau5
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! uh1: i+1/2 value for stencil 1 
      !      with 3 pints upwind shmeme
      ! uh2: i+1/2 value for stencil 2 
      !      with 3 pints upwind shmeme 
      ! uh3: i+1/2 value for stencil 3 
      !      with 3 pints upwind shmeme
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real(8) :: eps
      !
      eps=1.d-6
      !
      uh1=0.333333333333333d0*u(1)-1.16666666666667d0*u(2)+            &
          1.83333333333333d0*u(3)                                      
      uh2=-0.166666666666667d0*u(2)+0.833333333333333d0*u(3)+          &
           0.333333333333333d0*u(4)                                    
      uh3=0.333333333333333d0*u(3)+0.833333333333333d0*u(4)           &
          -0.166666666666667d0*u(5)                                    
      !                                                                
      beter1=1.08333333333333d0*(u(1)-2.d0*u(2)+u(3))**2+              &
             0.25d0*(u(1)-4.d0*u(2)+3.d0*u(3))**2                      
      beter2=1.08333333333333d0*(u(2)-2.d0*u(3)+u(4))**2+              &
             0.25d0*(u(2)-u(4))**2                                     
      beter3=1.08333333333333d0*(u(3)-2.d0*u(4)+u(5))**2+              &
             0.25d0*(3.d0*u(3)-4.d0*u(4)+u(5))**2
      !
      C1=0.1d0
      C2=0.6d0
      C3=0.3d0
      !
      tau5=dabs(beter3-beter1)
      !
      alfa1=C1+C1*(tau5/(beter1+eps)**2)
      alfa2=C2+C2*(tau5/(beter2+eps)**2)
      alfa3=C3+C3*(tau5/(beter3+eps)**2)
      !
      alfaS=alfa1+alfa2+alfa3
      !
      WT1=alfa1/alfaS
      WT2=alfa2/alfaS
      WT3=alfa3/alfaS
      !
      uh=WT1*uh1+WT2*uh2+WT3*uh3
      !
      return
      !
    end function WENO5Z
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! End of the subroutine WENO5Z.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! this subroutine offers the interface's values with 7 point WENO
    ! Scheme according the points' value you give.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    pure function WENO7Z(u) result(uh)
      !
      real(8),intent(in) :: u(1:7)
      real(8) :: uh
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! dir=1: upwind
      ! dir=-1: backwind
      ! u: 5 points grid value
      ! uh: final interface value
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real(8) :: beter1,beter2,beter3,beter4,alfa1,alfa2,alfa3,alfa4,  &
                 WT1,WT2,WT3,WT4,C1,C2,C3,C4,alfaS,tau7,uh1,uh2,uh3,uh4
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! uh1: i+1/2 value for stencil 1 
      !      with 3 pints upwind shmeme
      ! uh2: i+1/2 value for stencil 2 
      !      with 3 pints upwind shmeme 
      ! uh3: i+1/2 value for stencil 3 
      !      with 3 pints upwind shmeme
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real(8) :: df1,df2,df3
      real(8) :: eps
      !
      eps=1.d-6
      !
      uh1=-0.25d0              *u(1)+1.083333333333333d0 *u(2)       &
          -1.916666666666667d0 *u(3)+2.083333333333333d0 *u(4)
      uh2= 8.333333333333333d-2*u(2)-4.166666666666667d-1*u(3)       &
          +1.083333333333333d0 *u(4)+0.25d0              *u(5)
      uh3=-8.333333333333333d-2*u(3)+5.833333333333333d-1*u(4)       &
          +5.833333333333333d-1*u(5)-8.333333333333333d-2*u(6)
      uh4= 0.25d0              *u(4)+1.083333333333333d0 *u(5)       &
          -4.166666666666667d-1*u(6)+8.333333333333333d-2*u(7)
      !
      df1=2.777777777777778d-2*(-2.d0*u(1)+9.d0*u(2)-18.d0*u(3)+11d0*u(4))**2
      df2= 1.083333333333333d0*(-1.d0*u(1)+4.d0*u(2) -5.d0*u(3) +2d0*u(4))**2
      df3= 1.084722222222222d0*(-1.d0*u(1)+3.d0*u(2) -3.d0*u(3) +1d0*u(4))**2
      beter1=df1+df2+df3
      !
      df1=2.777777777777778d-2*(      u(2)-6.d0*u(3)+3.d0*u(4)+2.d0*u(5))**2
      df2= 1.083333333333333d0*(                u(3)-2.d0*u(4)+     u(5))**2
      df3= 1.084722222222222d0*(-1.d0*u(2)+3.d0*u(3)-3.d0*u(4)+1.d0*u(5))**2
      beter2=df1+df2+df3
      !
      df1=2.777777777777778d-2*(-2.d0*u(3)-3.d0*u(4)+6.d0*u(5)-1.d0*u(6))**2
      df2= 1.083333333333333d0*(      u(3)-2.d0*u(4)     +u(5)          )**2
      df3= 1.084722222222222d0*(-1.d0*u(3)+3.d0*u(4)-3.d0*u(5)+1.d0*u(6))**2
      beter3=df1+df2+df3
      !
      df1=2.777777777777778d-2*(-11.d0*u(4)+18.d0*u(5)-9.d0*u(6)+2.d0*u(7))**2
      df2= 1.083333333333333d0*(  2.d0*u(4) -5.d0*u(5)+4.d0*u(6)     -u(7))**2
      df3= 1.084722222222222d0*( -1.d0*u(4) +3.d0*u(5)-3.d0*u(6)+1.d0*u(7))**2
      beter4=df1+df2+df3
      !
      C1=2.857142857142857d-2
      C2=3.428571428571429d-1
      C3=5.142857142857143d-1
      C4=1.142857142857143d-1
      !
      tau7=dabs(beter4-beter1)
      !
      alfa1=C1+C1*(tau7/(beter1+eps)**2)
      alfa2=C2+C2*(tau7/(beter2+eps)**2)
      alfa3=C3+C3*(tau7/(beter3+eps)**2)
      alfa4=C4+C4*(tau7/(beter4+eps)**2)
      !
      alfaS=alfa1+alfa2+alfa3+alfa4
      !
      WT1=alfa1/alfaS
      WT2=alfa2/alfaS
      WT3=alfa3/alfaS
      WT4=alfa4/alfaS
      !
      uh=WT1*uh1+WT2*uh2+WT3*uh3+WT4*uh4
      !
      return
      !
    end function WENO7Z
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !End of the subroutine WENO7Z.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module flux
