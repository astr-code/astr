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
    use commtype, only : compact_scheme

    implicit none

    real(8),allocatable,dimension(:,:),target :: cfluxi_uwd,cfluxj_uwd,cfluxk_uwd, & !cflux: compact flux calculation coefficient
                                                 cfluxi_dwd,cfluxj_dwd,cfluxk_dwd
    integer :: ifxs,ifxe,jfxs,jfxe,kfxs,kfxe ! starting and ending nodes of finite difference schemes
    type(compact_scheme) :: flux_uw_i,flux_dw_i,flux_uw_j,flux_dw_j, &
                            flux_uw_k,flux_dw_k

    contains

    !+-------------------------------------------------------------------+
    !| This subroutine is to initiate the comapct flux scheme.           |
    !+-------------------------------------------------------------------+
    !| CHANGE RECORD                                                     |
    !| -------------                                                     |
    !| 03-06-2025  | Rewrite by J. Fang @ Liverpool.                     |
    !+-------------------------------------------------------------------+
    subroutine compact_flux_initiate(asolver,scheme,ntype,dim,wind)

        use constdef
        use commfunc,only : tridiagonal_thomas_proprocess

        type(compact_scheme) :: asolver
        integer,intent(in) :: scheme,ntype,dim
        character(len=1),intent(in) :: wind

        integer :: i_0,i_m,j

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

        asolver%first_node=i_0
        asolver%last_node =i_m
        asolver%dimension =dim
        asolver%nbctype   =ntype
        asolver%wind      =wind

        call asolver%init()

        if(scheme/100==5) then

          if(asolver%wind=='+') then
            asolver%a(:)=0.5d0;  asolver%c(:)=num1d6
          elseif(asolver%wind=='-') then
            asolver%a(:)=num1d6; asolver%c(:)=0.5d0
          else
            stop 'wind error @ compact_flux_initiate'
          endif

          ! default setup, interface
          asolver%a(i_0)=0.d0;     asolver%c(i_0)=0.d0 ! explicit central scheme at interface
          asolver%a(i_m)=0.d0;     asolver%c(i_m)=0.d0 ! explicit central scheme at interface

          if(scheme==543) then  
            ! set near-boundary/interface schemes
            if(ntype==1 .or. ntype==4) then
              asolver%a(i_0)   =2.d0;    asolver%c(i_0)  =2.d0  ! 3nd-order downwind-biased scheme
              asolver%a(i_0+1) =0.25d0;  asolver%c(i_0+1)=0.25d0  ! 3nd-order downwind-biased scheme
            endif
    
            if(ntype==2 .or. ntype==4) then
              asolver%a(i_m)   =2.d0;    asolver%c(i_m)  =2.d0    ! 3nd-order upwind-biased scheme        
              asolver%a(i_m-1) =0.25d0;  asolver%c(i_m-1)=0.25d0    ! 3nd-order upwind-biased scheme
            endif
          else
            print*,' scheme:',scheme
            stop ' !! scheme not defined @ compact_flux_initiate'
          endif

        else

          stop ' !! scheme not defined @ compact_flux_initiate'

        endif

        asolver%ac=tridiagonal_thomas_proprocess(asolver%a,asolver%c)

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
    function flux_compact(asolver,f,dim) result(fh)

        use commvar, only : hm
        use commfunc,only : tridiagonal_thomas_solver

        ! arguments
        type(compact_scheme),intent(in) :: asolver
        integer,intent(in) :: dim
        real(8),intent(in) :: f(-hm:dim+hm)
        real(8) :: fh(-1:dim)

        ! local data
        integer :: l
        real(8),allocatable :: d(:),xx(:)

        d=compact_flux_rhs(asolver,f,dim)

        allocate(xx(-2:dim+1))

        xx(asolver%first_node:asolver%last_node)=tridiagonal_thomas_solver(asolver%ac,d)

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
    function compact_flux_rhs(asolver,f,dim) result(d)

        use constdef
        use commvar, only : hm

        type(compact_scheme),intent(in) :: asolver
        integer,intent(in) :: dim
        real(8),intent(in) :: f(-hm:dim+hm)
        real(8),allocatable :: d(:)

        integer :: i_0,i_m,j,k,i_s,i_e,ntype
        real(8) :: var1,var2,var3
        character(len=1) :: wind

        i_0=asolver%first_node
        i_m=asolver%last_node

        ntype=asolver%nbctype

        wind=asolver%wind

        allocate(d(i_0:i_m))

        i_s=i_0+1
        i_e=i_m-1

        if(ntype==1 .or. ntype==4) then

          ! physical boundary

          ! 3rd-order flux
          j=i_0 ! i=-1
          d(j)=2.5d0*f(j+1)+0.5d0*f(j+2)

          ! 4th-order flux
          j=i_0+1 ! i=0
          d(j)=0.75d0*f(j)+0.75d0*f(j+1)

          i_s=i_0+2

        else

          ! interfaces

          ! 6th-order explicit
          j=i_0
          var1=f(j)  +f(j+1)
          var2=f(j-1)+f(j+2)
          var3=f(j-2)+f(j+3)
          d(j)=num37d60*var1 - num2d15*var2 + num1d60*var3

        endif

        if(ntype==2 .or. ntype==4) then

          ! physical boundary

          ! 3rd-order flux i=im-1
          j=i_m-1
          d(j)=0.75d0*f(j+1)+0.75d0*f(j)

          ! 4th-order flux
          j=i_m
          d(j)=2.5d0*f(j)+0.5d0*f(j-1)

          i_e=i_m-2
        else
          ! interfaces

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
        fc=mp5(f,fl,discont=.true.)
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

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! This subroutine is used to initial for solving the tridiagonal 
  ! martix with two layer of boundary scheme:, considering asymmetry 
  ! upwind scheme
  ! A*x=b
  !   |1,af1,...............|
  !   |af2,1,af3,.... ......|
  !   |..af4,1,af5,.........|
  !   |.....................|
  ! A=|...,af4,1,af5,.......|
  !   |.....................|
  !   |.........,af4,1,af5..|
  !   |...........,af2,1,af3|
  !   |...............,af1,1|
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Writen by Fang Jian, 2008-11-04.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function recons(f,stype,ntype,dim,af,cc,windir) result(fc)
    !
    use commvar,  only: hm
    ! arguments
    character(len=4),intent(in) :: stype
    integer,intent(in) :: ntype,dim
    real(8),intent(in) :: f(-hm:dim+hm)
    real(8),intent(in) :: af(1:5),cc(1:2,-2:dim+1)
    character(len=1) :: windir
    real(8) :: fc(-1:dim)
    !
    ! local data
    integer :: nscheme,k
    real(8) :: b(-2:dim+1)
    ! integer(8),save :: plan_f,plan_b
    !
    complex(8),allocatable :: cf(:)
    real(8),allocatable :: kama(:)
    !
    ! print*,mpirank,'|',dim,im
    !
    read(stype(1:3),*) nscheme
    !
    if(dim==0) then
      fc=f(0)
    else
      b =ptds_recon_rhs(f,dim,nscheme,ntype,windir)
      fc=ptds_recon_cal(b,af,cc,dim,ntype,windir)
    endif
    !
    return
    !
  end function recons

  function coeffcompac(scheme) result(alfa)
    !
    use constdef
    use commvar,  only: bfacmpld
    !
    integer,intent(in) :: scheme
    real(8),allocatable :: alfa(:)
    !
    if(scheme==642) then
      allocate(alfa(3))
      alfa(3)=num1d3
      alfa(2)=0.25d0
      alfa(1)=1.d0
    elseif(scheme==644) then
      allocate(alfa(3))
      alfa(3)=num1d3
      alfa(2)=0.25d0
      alfa(1)=0.d0
    elseif(scheme==553) then
      allocate(alfa(5))
      alfa(1)=0.5d0
      alfa(2)=0.5d0 -1.d0/6.d0*bfacmpld
      alfa(3)=num1d6+1.d0/6.d0*bfacmpld
      alfa(4)=0.5d0 -1.d0/6.d0*bfacmpld
      alfa(5)=num1d6+1.d0/6.d0*bfacmpld
    elseif(scheme==543) then
      allocate(alfa(5))
      alfa(1)=0.5d0
      alfa(2)=0.25d0
      alfa(3)=0.25d0
      alfa(4)=0.5d0
      alfa(5)=num1d6
    elseif(scheme==753) then
      allocate(alfa(5))
      alfa(1)=0.5d0
      alfa(2)=0.5d0 -1.d0/6.d0*bfacmpld
      alfa(3)=num1d6+1.d0/6.d0*bfacmpld
      alfa(4)=0.5d0 -1.d0/8.d0*bfacmpld
      alfa(5)=0.25d0+1.d0/8.d0*bfacmpld
    else
      stop ' !! scheme not defined @ coef_diffcompac'
    endif
    !
    return
    !
  end function coeffcompac

  subroutine ptds_aym_ini(cc,af,dim,ntype,windir)
    !
    integer,intent(in) :: dim,ntype
    real(8),intent(in) :: af(1:5)
    real(8),allocatable,intent(out) :: cc(:,:)
    character(len=1),intent(in) :: windir
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! af(1),af2,af: input dat
    ! cc: output array
    ! dim: input dat
    ! l: temporary variable
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! local data
    integer :: l
    !
    allocate(cc(1:2,-2:dim+1))
    !
    if(dim==0) then
      cc=0.d0
      return
    endif
    !
    if(ntype==1) then
      ! the block with boundary at i==0
      !
      cc(1,0)=1.d0
      cc(2,0)=af(1)/cc(1,0)
      !
      if(windir=='+') then
        cc(1,1)=1.d0-af(2)*cc(2,0)
        cc(2,1)=af(3)/cc(1,1)
        do l=2,dim+1
          cc(1,l)=1.d0-af(4)*cc(2,l-1)
          cc(2,l)=af(5)/cc(1,l)
        end do
      elseif(windir=='-') then
        cc(1,1)=1.d0-af(3)*cc(2,0)
        cc(2,1)=af(2)/cc(1,1)
        do l=2,dim+1
          cc(1,l)=1.d0-af(5)*cc(2,l-1)
          cc(2,l)=af(4)/cc(1,l)
        end do
      endif
      !
      cc(1,dim+1)=1.d0
      !
    elseif(ntype==2) then
      ! the block with boundary at i==im
      !
      cc(1,-2)=1.d0
      cc(2,-2)=0.d0
      !
      if(windir=='+') then
        do l=-1,dim-3
          cc(1,l)=1.d0-af(4)*cc(2,l-1)
          cc(2,l)=af(5)/cc(1,l)
        end do
        l=dim-2
        cc(1,l)=1.d0-af(2)*cc(2,l-1)
        cc(2,l)=af(3)/cc(1,l)
        l=dim-1
        cc(1,l)=1.d0-af(1)*cc(2,l-1)
      elseif(windir=='-') then
        do l=-1,dim+1
          cc(1,l)=1.d0-af(5)*cc(2,l-1)
          cc(2,l)=af(4)/cc(1,l)
        end do
        l=dim-2
        cc(1,l)=1.d0-af(3)*cc(2,l-1)
        cc(2,l)=af(2)/cc(1,l)
        l=dim-1
        cc(1,l)=1.d0-af(1)*cc(2,l-1)
      endif
      !
    elseif(ntype==3) then
      ! inner block
      !
      cc(1,-2)=1.d0
      cc(2,-2)=0.d0
      !
      if(windir=='+') then
        do l=-1,dim+1
          cc(1,l)=1.d0-af(4)*cc(2,l-1)
          cc(2,l)=af(5)/cc(1,l)
        end do
      elseif(windir=='-') then
        do l=-1,dim+1
          cc(1,l)=1.d0-af(5)*cc(2,l-1)
          cc(2,l)=af(4)/cc(1,l)
        end do
      endif
      !
      cc(1,dim+1)=1.d0
      !
    else
      print*, ' !! error in subroutine ptds_aym_ini !'
      stop
    end if
    !
    cc(1,:)=1.d0/cc(1,:)
    !
    return
    !
  end subroutine ptds_aym_ini
  !
  function ptds_recon_rhs(vin,dim,ns,ntype,windir) result(vout)
    !
    use constdef
    use commvar,  only: bfacmpld,hm
    !
    integer,intent(in) :: dim,ns,ntype
    real(8),intent(in) :: vin(-hm:dim+hm)
    character(len=1),intent(in) :: windir
    real(8) :: vout(-2:dim+1)
    !
    ! local data
    integer :: l
    !
    if(ntype==1) then
      ! the block with boundary at i==0
      !
      if(ns==553) then
        ! ns==553 3-5-5-...-5-5-3
        l=0
        vout(0)=0.25d0*vin(l)+1.25d0*vin(l+1)
        if(windir=='+') then
          do l=1,dim
            vout(l)=(num1d18 -num1d36*bfacmpld)*vin(l-1) + &
                    (num19d18-num9d36*bfacmpld)*vin(l)   + &
                    (num5d9  +num9d36*bfacmpld)*vin(l+1) + &
                              num1d36*bfacmpld *vin(l+2)
          end do
          l=dim+1
          vout(l)=(2.d0*vin(l-2)-13.d0*vin(l-1)+47.d0*vin(l)+          &
                   27.d0*vin(l+1)-3.d0*vin(l+2))/60.d0
        elseif(windir=='-') then
          do l=1,dim
            vout(l)=(num1d18 -num1d36*bfacmpld)*vin(l+2) + &
                    (num19d18-num9d36*bfacmpld)*vin(l+1) + &
                    (num5d9  +num9d36*bfacmpld)*vin(l)   + &
                              num1d36*bfacmpld *vin(l-1)
          end do
          l=dim+1
          vout(l)=(2.d0*vin(l+3)-13.d0*vin(l+2)+47.d0*vin(l+1)+        &
                   27.d0*vin(l)  -3.d0*vin(l-1))/60.d0
        endif
      elseif(ns==543) then

        ! ns==543 3-5-5-...-5-5-3
        l=0
        vout(l)=0.25d0*vin(l)+1.25d0*vin(l+1)
        l=1
        vout(l)=0.75d0*vin(l)+0.75d0*vin(l+1)
        if(windir=='+') then
          do l=2,dim
            vout(l)=(num1d18 )*vin(l-1) + &
                    (num19d18)*vin(l)   + &
                    (num5d9  )*vin(l+1)
          end do
        elseif(windir=='-') then
          do l=2,dim
            vout(l)=(num1d18)*vin(l+2) + &
                    (num19d18)*vin(l+1) + &
                    (num5d9)*vin(l)
          end do
        endif
        l=dim+1
        vout(l)=(vin(l-2)-8.d0*vin(l-1)+37.d0*vin(l)+          &
                 37.d0*vin(l+1)-8.d0*vin(l+2)+vin(l+3))/60.d0
        !
      elseif(ns==753) then
        ! ns==553 3-5-7-...-7-5-3
        !
        ! 3rd-order compact upwind
        l=0
        vout(0)=0.25d0*vin(l)+1.25d0*vin(l+1)
        !
        if(windir=='+') then
          !
          ! 5th-order compact upwind near boundary
          l=1
          vout(l)=(num1d18 -num1d36*bfacmpld)*vin(l-1) + &
                  (num19d18-num9d36*bfacmpld)*vin(l)   + &
                  (num5d9  +num9d36*bfacmpld)*vin(l+1) + &
                            num1d36*bfacmpld *vin(l+2)
          !
          ! 7th-order compact upwind
          do l=2,dim
            vout(l)=(( -1.d0 + 0.5d0*bfacmpld)*vin(l-2) +   &
                     ( 19.d0 - 7.5d0*bfacmpld)*vin(l-1) +   &
                     (239.d0 - 40.d0*bfacmpld)*vin(l)   +   &
                     (159.d0 + 40.d0*bfacmpld)*vin(l+1) +   &
                     (  4.d0 + 7.5d0*bfacmpld)*vin(l+2)     &
                             - 0.5d0*bfacmpld *vin(l+3))/240.d0
          end do
          !
          ! 7th-order explicit upwind
          l=dim+1
          vout(l)=(   -3.d0*vin(l-3)+  &
                      25.d0*vin(l-2)-  &
                     101.d0*vin(l-1)+  &
                     319.d0*vin(l)+    &
                     214.d0*vin(l+1)-  &
                      38.d0*vin(l+2)+  &
                       4.d0*vin(l+3) )/420.d0
          !
        elseif(windir=='-') then
          !
          ! 5th-order compact upwind near boundary
          l=1
          vout(l)=(num1d18 -num1d36*bfacmpld)*vin(l+2) + &
                  (num19d18-num9d36*bfacmpld)*vin(l+1) + &
                  (num5d9  +num9d36*bfacmpld)*vin(l)   + &
                            num1d36*bfacmpld *vin(l-1)
          !
          ! 7th-order compact upwind
          do l=2,dim
            vout(l)=(( -1.d0 + 0.5d0*bfacmpld)*vin(l+3) +   &
                     ( 19.d0 - 7.5d0*bfacmpld)*vin(l+2) +   &
                     (239.d0 - 40.d0*bfacmpld)*vin(l+1) +   &
                     (159.d0 + 40.d0*bfacmpld)*vin(l)   +   &
                     (  4.d0 + 7.5d0*bfacmpld)*vin(l-1)     &
                             - 0.5d0*bfacmpld *vin(l-2))/240.d0
          end do
          !
          ! 7th-order explicit upwind for interface
          l=dim+1
          vout(l)=(   -3.d0*vin(l+4)+  &
                      25.d0*vin(l+3)-  &
                     101.d0*vin(l+2)+  &
                     319.d0*vin(l+1)+  &
                     214.d0*vin(l)-    &
                      38.d0*vin(l-1)+  &
                       4.d0*vin(l-2) )/420.d0
        endif
        !
      end if
      !
    elseif(ntype==2) then
      ! the block with boundary at i==im
      !
      if(ns==553) then
        ! ns==553 3-5-5-...-5-5-3
        !
        if(windir=='+') then
          l=-2
          vout(l)=(2.d0*vin(l-2)-13.d0*vin(l-1)+47.d0*vin(l)+          &
                   27.d0*vin(l+1)-3.d0*vin(l+2))/60.d0
          do l=-1,dim-2
            vout(l)=(num1d18 -num1d36*bfacmpld)*vin(l-1) + &
                    (num19d18-num9d36*bfacmpld)*vin(l)   + &
                    (num5d9  +num9d36*bfacmpld)*vin(l+1) + &
                              num1d36*bfacmpld *vin(l+2)
          end do
        elseif(windir=='-') then
          l=-2
          vout(l)=(2.d0*vin(l+3)-13.d0*vin(l+2)+47.d0*vin(l+1)+        &
                   27.d0*vin(l)  -3.d0*vin(l-1))/60.d0
          do l=-1,dim-2
            vout(l)=(num1d18 -num1d36*bfacmpld)*vin(l+2) + &
                    (num19d18-num9d36*bfacmpld)*vin(l+1) + &
                    (num5d9  +num9d36*bfacmpld)*vin(l)   + &
                              num1d36*bfacmpld *vin(l-1)
          end do
        endif
        l=dim-1
        vout(l)=0.25d0*vin(l+1)+1.25d0*vin(l)
        !
      elseif(ns==543) then
        ! ns==543 3-5-5-...-5-4-3
        !
        l=-2
        vout(l)=(vin(l-2)-8.d0*vin(l-1)+37.d0*vin(l)+          &
                 37.d0*vin(l+1)-8.d0*vin(l+2)+vin(l+3))/60.d0
        if(windir=='+') then
          do l=-1,dim-3
            vout(l)=(num1d18 )*vin(l-1) + &
                    (num19d18)*vin(l)   + &
                    (num5d9  )*vin(l+1)
          end do

        elseif(windir=='-') then
          do l=-1,dim-3
            vout(l)=(num1d18 )*vin(l+2) + &
                    (num19d18)*vin(l+1) + &
                    (num5d9  )*vin(l) 
          end do
        endif
        l=dim-2
        vout(l)=0.75d0*vin(l)+0.75d0*vin(l+1)
        l=dim-1
        vout(l)=0.25d0*vin(l+1)+1.25d0*vin(l)
        !
      elseif(ns==753) then
        ! ns==753 3-5-7-...-7-5-3
        !
        if(windir=='+') then
          !
          l=-2
          vout(l)=(   -3.d0*vin(l-3)+  &
                      25.d0*vin(l-2)-  &
                     101.d0*vin(l-1)+  &
                     319.d0*vin(l)+    &
                     214.d0*vin(l+1)-  &
                      38.d0*vin(l+2)+  &
                       4.d0*vin(l+3) )/420.d0
          ! 7th-order compact upwind
          do l=-1,dim-3
            vout(l)=(( -1.d0 + 0.5d0*bfacmpld)*vin(l-2) +   &
                     ( 19.d0 - 7.5d0*bfacmpld)*vin(l-1) +   &
                     (239.d0 - 40.d0*bfacmpld)*vin(l)   +   &
                     (159.d0 + 40.d0*bfacmpld)*vin(l+1) +   &
                     (  4.d0 + 7.5d0*bfacmpld)*vin(l+2)     &
                             - 0.5d0*bfacmpld *vin(l+3))/240.d0
          end do
          !
          ! 5th-order compact upwind near the boundary
          l=dim-2
          vout(l)=(num1d18 -num1d36*bfacmpld)*vin(l-1) + &
                  (num19d18-num9d36*bfacmpld)*vin(l)   + &
                  (num5d9  +num9d36*bfacmpld)*vin(l+1) + &
                            num1d36*bfacmpld *vin(l+2)
          !
        elseif(windir=='-') then
          !
          ! 7th-order explicit upwind for interface
          l=-2
          vout(l)=(   -3.d0*vin(l+4)+  &
                      25.d0*vin(l+3)-  &
                     101.d0*vin(l+2)+  &
                     319.d0*vin(l+1)+  &
                     214.d0*vin(l)-    &
                      38.d0*vin(l-1)+  &
                       4.d0*vin(l-2) )/420.d0
          ! 7th-order compact upwind
          do l=-1,dim-3
            vout(l)=(( -1.d0 + 0.5d0*bfacmpld)*vin(l+3) +   &
                     ( 19.d0 - 7.5d0*bfacmpld)*vin(l+2) +   &
                     (239.d0 - 40.d0*bfacmpld)*vin(l+1) +   &
                     (159.d0 + 40.d0*bfacmpld)*vin(l)   +   &
                     (  4.d0 + 7.5d0*bfacmpld)*vin(l-1)     &
                             - 0.5d0*bfacmpld *vin(l-2))/240.d0
          end do
          !
          ! 5th-order compact upwind near the boundary
          l=dim-2
          vout(l)=(num1d18 -num1d36*bfacmpld)*vin(l+2) + &
                  (num19d18-num9d36*bfacmpld)*vin(l+1) + &
                  (num5d9  +num9d36*bfacmpld)*vin(l)   + &
                            num1d36*bfacmpld *vin(l-1)
          !
        endif
        !
        ! 3rd-order compact upwind at the boundary
        l=dim-1
        vout(l)=0.25d0*vin(l+1)+1.25d0*vin(l)
        !
      else
        print*,' ** ntype: ',ntype
        print*, ' !! error 2 in subroutine ptds_recon_rhs !'
        stop
      endif
      !
    elseif(ntype==3) then
      !
      ! inner block
      if(ns/100==5) then
        !
        if(windir=='+') then
          !
          l=-2
          vout(l)=(vin(l-2)-8.d0*vin(l-1)+37.d0*vin(l)+          &
                   37.d0*vin(l+1)-8.d0*vin(l+2)+vin(l+3))/60.d0
          do l=-1,dim
            vout(l)=(num1d18 )*vin(l-1) + &
                    (num19d18)*vin(l)   + &
                    (num5d9  )*vin(l+1)
          end do
          l=dim+1
          vout(l)=(vin(l-2)-8.d0*vin(l-1)+37.d0*vin(l)+          &
                   37.d0*vin(l+1)-8.d0*vin(l+2)+vin(l+3))/60.d0
        elseif(windir=='-') then
          l=-2
          vout(l)=(vin(l-2)-8.d0*vin(l-1)+37.d0*vin(l)+          &
                   37.d0*vin(l+1)-8.d0*vin(l+2)+vin(l+3))/60.d0
          do l=-1,dim
            vout(l)=(num1d18 )*vin(l+2) + &
                    (num19d18)*vin(l+1) + &
                    (num5d9  )*vin(l)  
          end do
          l=dim+1
          vout(l)=(vin(l-2)-8.d0*vin(l-1)+37.d0*vin(l)+          &
                   37.d0*vin(l+1)-8.d0*vin(l+2)+vin(l+3))/60.d0
        endif
        !
      elseif(ns/100==7) then
        !
        if(windir=='+') then
          !
          ! 7th-order explicit upwind for interface
          l=-2
          vout(l)=(   -3.d0*vin(l-3)+  &
                      25.d0*vin(l-2)-  &
                     101.d0*vin(l-1)+  &
                     319.d0*vin(l)+    &
                     214.d0*vin(l+1)-  &
                      38.d0*vin(l+2)+  &
                       4.d0*vin(l+3) )/420.d0
          ! 
          ! 7th-order compact scheme
          do l=-1,dim
            vout(l)=(( -1.d0 + 0.5d0*bfacmpld)*vin(l-2) +   &
                     ( 19.d0 - 7.5d0*bfacmpld)*vin(l-1) +   &
                     (239.d0 - 40.d0*bfacmpld)*vin(l)   +   &
                     (159.d0 + 40.d0*bfacmpld)*vin(l+1) +   &
                     (  4.d0 + 7.5d0*bfacmpld)*vin(l+2)     &
                             - 0.5d0*bfacmpld *vin(l+3))/240.d0
          end do
          !
          ! 7th-order explicit upwind for interface
          l=dim+1
          vout(l)=(   -3.d0*vin(l-3)+  &
                      25.d0*vin(l-2)-  &
                     101.d0*vin(l-1)+  &
                     319.d0*vin(l)+    &
                     214.d0*vin(l+1)-  &
                      38.d0*vin(l+2)+  &
                       4.d0*vin(l+3) )/420.d0
          !
        elseif(windir=='-') then
          !
          ! 7th-order explicit upwind for interface
          l=-2
          vout(l)=(   -3.d0*vin(l+4)+  &
                      25.d0*vin(l+3)-  &
                     101.d0*vin(l+2)+  &
                     319.d0*vin(l+1)+  &
                     214.d0*vin(l)-    &
                      38.d0*vin(l-1)+  &
                       4.d0*vin(l-2) )/420.d0
          ! 7th-order compact scheme
          do l=-1,dim
            vout(l)=(( -1.d0 + 0.5d0*bfacmpld)*vin(l+3) +   &
                     ( 19.d0 - 7.5d0*bfacmpld)*vin(l+2) +   &
                     (239.d0 - 40.d0*bfacmpld)*vin(l+1) +   &
                     (159.d0 + 40.d0*bfacmpld)*vin(l)   +   &
                     (  4.d0 + 7.5d0*bfacmpld)*vin(l-1)     &
                             - 0.5d0*bfacmpld *vin(l-2))/240.d0
          end do
          ! 7th-order explicit upwind for interface
          l=dim+1
          vout(l)=(   -3.d0*vin(l+4)+  &
                      25.d0*vin(l+3)-  &
                     101.d0*vin(l+2)+  &
                     319.d0*vin(l+1)+  &
                     214.d0*vin(l)-    &
                      38.d0*vin(l-1)+  &
                       4.d0*vin(l-2) )/420.d0
        endif
        !
      endif
      !
    else
      print*,' ** ntype: ',ntype
      print*, ' !! error in subroutine ptds_recon_rhs !'
      stop
    end if
  end function ptds_recon_rhs
  !
  function ptds_recon_cal(bd,af,cc,dim,ntype,windir) result(xd)
    !
    integer,intent(in) :: dim,ntype
    real(8),intent(in) :: af(5),bd(-2:dim+1),cc(1:2,-2:dim+1)
    character(len=1),intent(in) :: windir
    real(8) :: xd(-1:dim)
    !
    ! local data
    integer :: l
    real(8),allocatable,dimension(:) :: yd,md
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! af(3): input dat
    ! bd: input array
    ! xd: output array
    ! cc: input array
    ! dim: input dat
    ! l, yd: temporary variable
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    allocate ( yd(-2:dim+1),md(-2:dim+1) )
    !
    if(ntype==1) then
      ! the block with boundary at i==0
      !
      yd(0)=bd(0)*cc(1,0)
      if(windir=='+') then
        yd(1)=(bd(1)-af(2)*yd(0))*cc(1,1)
        do l=2,dim
          yd(l)=(bd(l)-af(4)*yd(l-1))*cc(1,l)
        end do
      elseif(windir=='-') then
        yd(1)=(bd(1)-af(3)*yd(0))*cc(1,1)
        do l=2,dim
          yd(l)=(bd(l)-af(5)*yd(l-1))*cc(1,l)
        end do
      endif
      yd(dim+1)=bd(dim+1)*cc(1,dim+1)
      !
      md(dim+1)=yd(dim+1)
      do l=dim,-2,-1
        md(l)=yd(l)-cc(2,l)*md(l+1)
      end do
      !
      xd(0:dim)=md(0:dim)
      !
    elseif(ntype==2) then
      ! the block with boundary at i==im
      !
      yd(-2)=bd(-2)*cc(1,-2)
      if(windir=='+') then
        do l=-1,dim-3
          yd(l)=(bd(l)-af(4)*yd(l-1))*cc(1,l)
        end do
        l=dim-2
        yd(l)=(bd(l)-af(2)*yd(l-1))*cc(1,l)
      elseif(windir=='-') then
        do l=-1,dim-3
          yd(l)=(bd(l)-af(5)*yd(l-1))*cc(1,l)
        end do
        l=dim-2
        yd(l)=(bd(l)-af(3)*yd(l-1))*cc(1,l)
      endif
      l=dim-1
      yd(l)=(bd(l)-af(1)*yd(l-1))*cc(1,l)
      !
      md(dim-1)=yd(dim-1)
      do l=dim-2,-2,-1
        md(l)=yd(l)-cc(2,l)*md(l+1)
      end do
      !
      xd(-1:dim-1)=md(-1:dim-1)
      !
    elseif(ntype==3) then
      !
      ! inner block
      !
      yd(-2)=bd(-2)*cc(1,-2)
      if(windir=='+') then
        do l=-1,dim
          yd(l)=(bd(l)-af(4)*yd(l-1))*cc(1,l)
        end do
      elseif(windir=='-') then
        do l=-1,dim
          yd(l)=(bd(l)-af(5)*yd(l-1))*cc(1,l)
        end do
      endif
      yd(dim+1)=bd(dim+1)*cc(1,dim+1)
      !
      md(dim+1)=yd(dim+1)
      do l=dim,-2,-1
        md(l)=yd(l)-cc(2,l)*md(l+1)
      end do
      !
      xd(-1:dim)=md(-1:dim)
      !
    else
      print*,' ** ntype: ',ntype
      print*, ' !! error in subroutine ptds_recon_cal !'
      stop
    end if
    !
    deallocate( yd )
    !
    return
    !
  end function ptds_recon_cal
  !+-------------------------------------------------------------------+
  !| The end of the function ptds_recon_.                              |
  !+-------------------------------------------------------------------+
  !
end module flux
