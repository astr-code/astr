module pastr_blasius

    use iso_fortran_env, only: wp => real64

    implicit none

    real(wp) :: Pr,gamma,Me,Re,f1inf,g0inf,f00,f10,g00,h,ref_t
    integer :: n

contains
    
    ! reference of defining compressible blasius equation: 
    ! https://www.mdpi.com/2311-5521/6/11/400
    subroutine blasius_solution

      use pastr_interpolation,only: interpolat
      use pastr_fsolver, only: fsolve
      use pastr_io,     only: parse_command_line

      real(wp) :: omega,theta,Twall
      real(wp) :: a,b,r,xin,rex,f20_guess,g0_guess,var1,var2
      real(wp),allocatable,dimension(:) :: eta,Me_study,f20_guess_store, &
                                           g0_guess_store,uin,uout
      real(wp),allocatable :: f(:,:),g(:,:),t(:),ro(:),u(:),v(:),y(:)
      real(wp) :: thick1,thick2,thick3,delta_target

      integer :: m,nguess
      integer :: i,j,io
      character(len=16) :: ecomm

      a = 0._wp
      b = 12._wp
      h = 0.01_wp
      n = (b-a)/h
      n = int(n)
      m = 3
      allocate(eta(0:n))
      do j=0,n
        eta(j) = (b-a)/dble(n)*dble(j)
      enddo
      print*,' ** number of nodes:',n
      print*,' ** range of eta:',eta(0),eta(n)
      ! k = ((2.0*(gamma-1.0)*ma**2.0)/(2.0+(gamma-1.0)*ma**2.0)) 
      ! coefficient ue^2/(h0)e^2 in the energy equation
      ! Flow parameters
      Pr = 0.72_wp
      gamma = 1.4_wp

      ! Boundary conditions
      f00 = 0._wp
      f10 = 0._wp
      f1inf = 1.0_wp
      g0inf = 1.0_wp

      Me=0.1_wp
      Re=100._wp
      ref_t=273.15_wp
      xin=1.d0
      delta_target=0._wp
      
      open(12,file='bldef.txt')
      read(12,*,iostat=io)Me,Re,ref_t,xin,delta_target
      read(12,*,iostat=io)Twall
      close(12)
      print*,' >> bldef.txt'
      if(io==0 .and. Twall > 0._wp) then
      else
        r=sqrt(Pr) !https://www.thermopedia.com/content/291/
        Twall=1.0_wp+0.5_wp*r*(gamma-1._wp)*Me*Me
      endif
      
      if(delta_target>0._wp) then
        xin=(delta_target/5.29_wp)**2*Re
      else
        delta_target=5.29_wp*sqrt(xin/Re)
      endif

      rex=Re*xin


      print*,' ** Mach number           :',Me
      print*,' ** Reynolds number       :',Re
      print*,' ** Distance to LE        :',xin
      print*,' ** Target BL thickness   :',delta_target
      print*,' ** Rex                   :',rex
      print*,' ** Reference temperature :',ref_t
      print*,' ** Wall temperature      :',Twall

      nguess=13
      allocate(Me_study(0:nguess),f20_guess_store(0:nguess),g0_guess_store(0:nguess))
      Me_study        = (/   0.0,   2.0,   3.1,   4.5,   5.0,   6.0,   7.0,   8.0,   9.0,  10.0,  10.3,  10.6,  11.0, 11.25/)
      f20_guess_store = (/0.4696,0.4746,0.4803,0.4924,0.4976,0.5094,0.5225,0.5367,0.5517,0.5671,0.5718,0.5766,0.5829,0.5869/)
      g0_guess_store  = (/   0.0,0.1765,0.2471, 0.309,0.3246,0.3503,0.3716,0.3904,0.4077,0.424, 0.4288,0.4335,0.4397,0.4436/)

      g00 = (Twall+0.5_wp*(Me**2)*(gamma-1.0_wp)*(f10**2))/(1.0+0.5*(Me**2)*(gamma-1._wp))

      f20_guess = interpolat(Me_study,f20_guess_store,Me)
      g0_guess  = interpolat(Me_study,g0_guess_store,Me)

      print*,' ** f20_guess',f20_guess,'g0_guess',g0_guess
      allocate(uin(2),uout(2))
      uin(1)=f20_guess
      uin(2)=g0_guess

      call fsolve(fcn=nonlinear_FG,n=2,x=uin,fvec=uout,tol=1.d-12,info=io )

      allocate(f(0:n,0:2),g(0:n,0:1))

      call update_FG(uin(1),uin(2),f,g)
      write(*,*)' ** f20_error = ',f(n,1),f1inf
      write(*,*)' ** g00_error = ',g(n,0),g0inf

      print*,' ** Blasius equation solved.'

      allocate(t(0:n),ro(0:n),u(0:n),v(0:n),y(0:n))

      t(:) = g(:,0) + 0.5_wp*(Me**2)*(gamma-1._wp)*(g(:,0)-(f(:,1))**2)
      ro   =1._wp/t
      u    =f(:,1)


      ! Howarth–Dorodnitsyn transformation
      ! https://en.wikipedia.org/wiki/Blasius_boundary_layer#Compressible_Blasius_boundary_layer
      y(0)=0._wp
      do j=1,n
        y(j)=0._wp
        do i=1,j
          y(j)=y(j)+0.5_wp*(t(i)+t(i-1))*h
        end do
      end do
      y=y*sqrt(2._wp*xin/Re)

      ! var1=0._wp
      ! do j=1,n
      !   var1=var1+0.5_wp*(ro(j)+ro(j-1))*(y(j)-y(j-1))*sqrt(Re/2._wp/xin)
      !   print*,j,var1,eta(j),y(j)
      ! enddo

      do j=0,n
        v(j)=-1._wp*t(j)*(f(j,0)-f(j,1)*eta(j))
      end do
      v=v/sqrt(2._wp*rex)

      thick1=0._wp
      thick2=0._wp
      thick3=0._wp

      thick1=interpolat(u,y,0.99_wp)

      do j=1,n
        var1=0.5_wp*(ro(j)+ro(j-1))
        var2=0.5_wp*(u(j)+u(j-1))
        thick2=thick2+(1._wp-var1*var2)*(y(j)-y(j-1))
      end do

      do j=1,n
        var1=0.5_wp*(ro(j)+ro(j-1))
        var2=0.5_wp*(u(j)+u(j-1))
        thick3=thick3+(var1*var2*(u(n)-var2))*(y(j)-y(j-1))
      end do

      print*,' ** norminal thickness     : ',thick1
      print*,' ** displacement thickness : ',thick2
      print*,' ** momentum thickness     : ',thick3

      print*,' ** Reynolds Based on x    : ',rex
      print*,' ** Reynolds Based on δ    : ',re*thick1
      print*,' ** Reynolds Based on δ*   : ',re*thick2
      print*,' ** Reynolds Based on θ    : ',re*thick3

      open(18,file='profile_blasius.dat')
      write(18,'(5(1X,A20))')'y','ro','u','v','T'
      do i=0,n
        write(18,'(5(1X,E20.13E2))')y(i),ro(i),u(i),v(i),t(i)
      enddo
      close(18)
      print*,' << profile_blasius.dat'

      call parse_command_line(  string=ecomm )
      if(trim(ecomm)=='inletgen') then
        print*,' ** generating inlet.prof'
      endif

    end subroutine blasius_solution

    subroutine rk4_fg(y)
      
      use pastr_constdef

      real(wp) :: y(0:n,0:4)

      real(wp) :: k1(0:4),k2(0:4),k3(0:4),k4(0:4)

      integer :: i

      do i=0,n-1
        k1=h*sysFG(y(i,:))
        k2=h*sysFG(y(i,:)+0.5_wp*k1)
        k3=h*sysFG(y(i,:)+0.5_wp*k2)
        k4=h*sysFG(y(i,:)+k3)
        y(i+1,:)=y(i,:) + num1d6*(k1+2.0_wp*k2+2.0_wp*k3+k4)
      enddo

    end subroutine rk4_fg

    function sysFG(u) result(dfsysFG)

      use pastr_thermo_phys, only: viscosity_suth

      real(wp) :: u(0:)
      real(wp) :: dfsysFG(0:4)

      real(wp) :: f(0:2),g(0:1)
      real(wp) :: k,invC,T_ratio

      f(0:2)=u(0:2)
      g(0:1)=u(3:4)

      k = ((2._wp*(gamma-1._wp)*Me**2)/(2._wp+(gamma-1._wp)*Me**2)) 
      T_ratio = g(0) + 0.5_wp*(Me**2)*(gamma-1._wp)*(g(0)-f(1)**2)

      invC=1._wp/viscosity_suth(T=1._wp,nondim=.true.,ref_temperature=ref_t)

      dfsysFG(0) = f(1)
      dfsysFG(1) = f(2)
      dfsysFG(2) = -invC*f(0)*f(2)
      dfsysFG(3) = g(1)
      dfsysFG(4) = -Pr*invC*f(0)*g(1) - K*(Pr-1._wp)*f(2)*(f(2)-invC*f(0)*f(1))

    end function sysFG

    subroutine nonlinear_FG( nx, x, fvec )
      
      integer :: nx
      real(wp) :: x(1:nx),fvec(1:nx)

      real(wp) :: u(0:n,0:4)

      u(0,0)=f00
      u(0,1)=f10
      u(0,2)=x(1)
      u(0,3)=g00
      u(0,4)=x(2)

      call rk4_fg(u)

      fvec(1)=u(n,1)-f1inf
      fvec(2)=u(n,3)-g0inf

      return

    end subroutine nonlinear_FG

    subroutine update_FG(f20,g0,f,g)
      
      real(wp),intent(in) :: f20,g0
      real(wp) :: f(0:n,0:2),g(0:n,0:1)

      real(wp) :: u(0:n,0:4)

      integer :: i
      real(wp) :: var

      u(0,0)=f00
      u(0,1)=f10
      u(0,2)=f20
      u(0,3)=g00
      u(0,4)=g0

      call rk4_fg(u)

      do i=0,n
        f(i,0)=u(i,0)
        f(i,1)=u(i,1)
        f(i,2)=u(i,2)
        g(i,0)=u(i,3)
        g(i,1)=u(i,4)
      enddo

      return

    end subroutine update_FG

end module pastr_blasius