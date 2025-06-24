!+---------------------------------------------------------------------+
!| This module contains routines related to numerical derivative       |
!| calculations.                                                       |
!+---------------------------------------------------------------------+
!| CHANGE RECORD                                                       |
!| -------------                                                       |
!| 30-05-2025  | Created by J. Fang @ Liverpool                        |
!+---------------------------------------------------------------------+
module derivative

    use commtype,only : compact_scheme

    implicit none

    ! Abstract base type
    type, abstract :: finite_difference
    contains
        procedure(derivative_interface), deferred :: central
    end type finite_difference
  
    ! Abstract interface for method
    abstract interface
        function derivative_interface(this,asolver,f,ntype,dim,dir) result(df)
            use commtype,only : compact_scheme
            use commvar, only : hm
            import :: finite_difference
            implicit none
            class(finite_difference), intent(in) :: this
            type(compact_scheme),intent(in) :: asolver
            real(8), intent(in) :: f(-hm:dim+hm)
            integer, intent(in) :: ntype,dim,dir
            real(8) :: df(0:dim)
        end function derivative_interface
    end interface
  
    class(finite_difference), allocatable :: fds !finite-difference-solver

    type, extends(finite_difference) :: compact_central
    contains
        procedure :: central => df_compact
    end type compact_central

    type, extends(finite_difference) :: explicit_central
    contains
        procedure :: central => df_explicit
    end type explicit_central

    type(compact_scheme) :: fds_compact_i,fds_compact_j,fds_compact_k

    contains

    !+-------------------------------------------------------------------+
    !| This subroutine is to initiate the comapct scheme.                |
    !+-------------------------------------------------------------------+
    !| CHANGE RECORD                                                     |
    !| -------------                                                     |
    !| 30-05-2025  | Rewrite by J. Fang @ Liverpool.                     |
    !+-------------------------------------------------------------------+
    subroutine fd_scheme_initiate(asolver,scheme,ntype,dim,dir)

        use constdef
        use commfunc,only : tridiagonal_thomas_proprocess

        type(compact_scheme) :: asolver
        character(len=4),intent(in) :: scheme
        integer,intent(in) :: ntype,dim,dir

        integer :: i_0,i_m,nscheme

        select case (ntype)
        case (1)
          ! the block with boundary at i==0
          i_0=0
          i_m=dim+1
        case (2)
          ! the block with boundary at i==im
          i_0=-1
          i_m=dim
        case (3)
          ! inner block
          i_0=-1
          i_m=dim+1
        case (4)
          ! the block with boundary at i=0 and i=im
          i_0=0
          i_m=dim
        case default
          print*, ' !! error 1 in subroutine compact_filter_initiate !'
        end select

        asolver%first_node=i_0
        asolver%last_node=i_m
        call asolver%init()

        asolver%nbctype=ntype

        read(scheme(1:3),*) nscheme

        if(scheme(4:4)=='c') then

          if(nscheme/100==6) then
  
            asolver%a(:)=num1d3
            asolver%c(:)=num1d3
  
            ! default setup, interface
            asolver%a(i_0)=0.d0;     asolver%c(i_0)=0.d0 ! explicit central scheme at interface
            asolver%a(i_m)=0.d0;     asolver%c(i_m)=0.d0 ! explicit central scheme at interface

            if(nscheme==642) then
  
              ! set near-boundary/interface schemes
              if(ntype==1 .or. ntype==4) then
                asolver%a(i_0)  =1.d0;   asolver%c(i_0)  =1.d0   ! 2nd-order compact biased scheme
                asolver%a(i_0+1)=0.25d0; asolver%c(i_0+1)=0.25d0 ! 4th-order compact central scheme
              endif
      
              if(ntype==2 .or. ntype==4) then
                asolver%a(i_m)  =1.d0;   asolver%c(i_m)  =1.d0   ! 2nd-order compact biased scheme
                asolver%a(i_m-1)=0.25d0; asolver%c(i_m-1)=0.25d0 ! 4th-order compact central scheme
              endif
  
            endif
  
          else
            stop ' !! scheme not defined @ coef_diffcompac'
          endif
  
          call tridiagonal_thomas_proprocess(asolver%a,asolver%c,asolver%ac)

        endif

    end subroutine fd_scheme_initiate
    !+-------------------------------------------------------------------+
    !| The end of the subroutine fd_scheme_initiate.                     |
    !+-------------------------------------------------------------------+

    !+-------------------------------------------------------------------+
    !| This function is to calculate the difference of a input array     | 
    !| using compact scheme                                              |
    !+-------------------------------------------------------------------+
    !| CHANGE RECORD                                                     |
    !| -------------                                                     |
    !| 30-05-2025  | Rewrite by J. Fang @ Liverpool.                     |
    !+-------------------------------------------------------------------+
    function df_compact(this,asolver,f,ntype,dim,dir) result(df)

        use commvar, only : hm
        use commfunc,only : tridiagonal_thomas_solver

        ! arguments
        class(compact_central), intent(in) :: this
        type(compact_scheme),intent(in) :: asolver
        integer,intent(in) :: ntype,dim,dir
        real(8),intent(in) :: f(-hm:dim+hm)
        real(8) :: df(0:dim)

        ! local data
        integer :: l
        real(8) :: alpha,var0,var1,var2,var3,var4,var5
        real(8),allocatable :: d(:),xx(:)

        d=compact_fd_rhs(asolver,f,ntype,dim,dir)

        allocate(xx(asolver%first_node:asolver%last_node))

        xx=tridiagonal_thomas_solver(asolver%ac,d)

        df(0:dim)=xx(0:dim)

        return

    end function df_compact
    !+-------------------------------------------------------------------+
    !| The end of the function df_compact.                               |
    !+-------------------------------------------------------------------+

    !+-------------------------------------------------------------------+
    !| This function is to calculate R.H.S of a compact finite-difference|
    !+-------------------------------------------------------------------+
    !| CHANGE RECORD                                                     |
    !| -------------                                                     |
    !| 30-05-2025  | Rewrite by J. Fang @ Liverpool.                     |
    !+-------------------------------------------------------------------+
    function compact_fd_rhs(asolver,f,ntype,dim,dir) result(d)

        use constdef
        use commvar, only : hm

        type(compact_scheme),intent(in) :: asolver
        integer,intent(in) :: ntype,dim,dir
        real(8),intent(in) :: f(-hm:dim+hm)
        real(8),allocatable :: d(:)

        integer :: i_0,i_m,j,k,i_s,i_e
        real(8) :: var1,var2,var3

        i_0=asolver%first_node
        i_m=asolver%last_node

        allocate(d(i_0:i_m))

        if(ntype==1 .or. ntype==4) then

          ! physical boundary

          i_s=i_0+2

          ! 4th-order compact
          j=i_0+1
          var1=f(j+1)-f(j-1)
          d(j)=0.75d0*var1
    
          ! 2nd-order biased compact
          j=i_0
          var1=f(j+1)-f(j)
          d(j)=2.d0*var1


        else

          ! interfaces
          i_s=i_0+1

          ! 6th-order explicit
          j=i_0
          var1=f(j+1)-f(j-1)
          var2=f(j+2)-f(j-2)
          var3=f(j+3)-f(j-3)
          d(j)=0.75d0*var1 - 0.15d0*var2 + num1d60*var3

        endif

        if(ntype==2 .or. ntype==4) then

          ! physical boundary
          i_e=i_m-2

          ! 4th-order compact
          j=i_m-1
          var1=f(j+1)-f(j-1)
          d(j)=0.75d0*var1
    
          ! 2nd-order biased compact
          j=i_m
          var1=f(j)-f(j-1)
          d(j)=2.d0*var1

        else
          ! interfaces
          i_e=i_m-1

          ! 6th-order explicit
          j=i_m
          var1=f(j+1)-f(j-1)
          var2=f(j+2)-f(j-2)
          var3=f(j+3)-f(j-3)
          d(j)=0.75d0*var1 - 0.15d0*var2 + num1d60*var3

        endif

        do j=i_s,i_e

          ! 6th-order compact
          var1=f(j+1)-f(j-1)
          var2=f(j+2)-f(j-2)

          d(j)=num7d9*var1+   num1d36*var2
    
        end do
    
    end function compact_fd_rhs
    !+-------------------------------------------------------------------+
    !| The end of the function compact_fd_rhs.                           |
    !+-------------------------------------------------------------------+

    !+-------------------------------------------------------------------+
    !| This function is to calculate the difference of a input array     | 
    !| using explicit scheme                                             |
    !+-------------------------------------------------------------------+
    !| CHANGE RECORD                                                     |
    !| -------------                                                     |
    !| 30-05-2025  | Rewrite by J. Fang @ Liverpool.                     |
    !+-------------------------------------------------------------------+
    function df_explicit(this,asolver,f,ntype,dim,dir) result(df)
  
        use commvar, only : hm
  
        ! arguments
        class(explicit_central), intent(in) :: this
        type(compact_scheme),intent(in) :: asolver
        integer,intent(in) :: ntype,dim,dir
        real(8),intent(in) :: f(-hm:dim+hm)
        real(8) :: df(0:dim)
  
        integer :: i
  
        df=diff6ec(f,dim,ntype)
        
        return
  
    end function df_explicit
    !+-------------------------------------------------------------------+
    !| The end of the function df_explicit.                              |
    !+-------------------------------------------------------------------+
  
    !+-------------------------------------------------------------------+
    !| This function returns the approximation of derivative of a given  |
    !| array, using 6th-order finite-difference scheme, including the    |
    !| boundary close.                                                   |
    !+-------------------------------------------------------------------+
    !| CHANGE RECORD                                                     |
    !| -------------                                                     |
    !| 31-05-2025  | Moved to here by J. Fang @ Liverpool.              |
    !+-------------------------------------------------------------------+
    pure function diff6ec(vin,dim,ntype) result(vout)
      
      use constdef
      use commvar, only : hm

      integer,intent(in) :: dim,ntype
      real(8),intent(in) :: vin(-hm:dim+hm)
      real(8) :: vout(0:dim)
      
      ! local data
      integer :: i

      select case (ntype)
      case (1)

        vout(0)=-0.5d0*vin(2)+2.d0*vin(1)-1.5d0*vin(0)
        vout(1)=0.5d0*(vin(2)-vin(0))
        vout(2)=num2d3*(vin(3)-vin(1))-num1d12*(vin(4)-vin(0))

        do i=3,dim
          vout(i)  =0.75d0 *(vin(i+1)-vin(i-1))-                         &
                    0.15d0 *(vin(i+2)-vin(i-2))+                         &
                    num1d60*(vin(i+3)-vin(i-3))
        enddo

      case (2)

        do i=0,dim-3
          vout(i)  =0.75d0 *(vin(i+1)-vin(i-1))-                         &
                    0.15d0 *(vin(i+2)-vin(i-2))+                         &
                    num1d60*(vin(i+3)-vin(i-3))
        enddo

        vout(dim-2) =num2d3*(vin(dim-1)-vin(dim-3))-                   &
                    num1d12*(vin(dim)  -vin(dim-4))
        vout(dim-1)=0.5d0*(vin(dim)-vin(dim-2))
        vout(dim)  =0.5d0*vin(dim-2)-2.d0*vin(dim-1)+1.5d0*vin(dim)
        
      case (3)

        do i=0,dim
          vout(i)  =0.75d0 *(vin(i+1)-vin(i-1))-                         &
                    0.15d0 *(vin(i+2)-vin(i-2))+                         &
                    num1d60*(vin(i+3)-vin(i-3))
        enddo

      case (4)

        vout(0)=-0.5d0*vin(2)+2.d0*vin(1)-1.5d0*vin(0)
        vout(1)=0.5d0*(vin(2)-vin(0))
        vout(2)=num2d3*(vin(3)-vin(1))-num1d12*(vin(4)-vin(0))
        do i=3,dim-3
          vout(i)  =0.75d0 *(vin(i+1)-vin(i-1))-                         &
                    0.15d0 *(vin(i+2)-vin(i-2))+                         &
                    num1d60*(vin(i+3)-vin(i-3))
        enddo
        vout(dim-2) =num2d3*(vin(dim-1)-vin(dim-3))-                   &
                    num1d12*(vin(dim)  -vin(dim-4))
        vout(dim-1)=0.5d0*(vin(dim)-vin(dim-2))
        vout(dim)  =0.5d0*vin(dim-2)-2.d0*vin(dim-1)+1.5d0*vin(dim)

      end select
      !
    end function diff6ec
    !+-------------------------------------------------------------------+
    !| The end of the function diff6ec.                                  |
    !+-------------------------------------------------------------------+

end module derivative
!+---------------------------------------------------------------------+
!| The end of the module derivative.                                   |
!+---------------------------------------------------------------------+