module pastr_fdm

    use iso_fortran_env, only: wp => real64
    use pastr_constdef

    implicit none

    type :: stencil_type
      real(wp),pointer :: a(:)
    end type stencil_type
    
    type :: fdm_type

      integer :: first_node,last_node,dim
      type(stencil_type),allocatable :: s(:)

      contains

      procedure :: init =>fdm_solver_init
      procedure :: init2 =>fdm_solver_init2
      procedure :: cal=>fdm_solver_operator
      procedure :: info=>fdm_solver_print

    end type fdm_type

contains
    
    function fdm_solver_operator(asolver,f,nhalo) result(df)

      class(fdm_type),target :: asolver
      integer,intent(in) :: nhalo
      real(wp),intent(in) :: f(-nhalo:asolver%dim+nhalo)
      real(wp) :: df(0:asolver%dim)

      integer :: n,j

      df=0.d0

      do n=0,asolver%dim
        do j=-3,3
          df(n)=df(n)+asolver%s(n)%a(j)*f(n+j)
        enddo
      enddo

    end function fdm_solver_operator

    subroutine fdm_solver_init2(asolver,dim,nhalo,node_state)

      class(fdm_type),target :: asolver
      integer,intent(in) :: dim,nhalo

      character(len=*),intent(in) :: node_state

      integer :: n
      character(len=9) :: cstencil
      
      asolver%dim=dim
      allocate(asolver%s(0:dim))

      do n=0,dim
        cstencil=node_state(n+1:n+3)//'-'// &
                 node_state(n+4:n+4)//'-'// &
                 node_state(n+5:n+7)

        allocate(asolver%s(n)%a(-3:3))
        asolver%s(n)%a=0.d0
        
        if(cstencil(1:9)=='fff-f-fff' .or. &
           cstencil(1:9)=='bff-f-fff' .or. &
           cstencil(1:9)=='fff-f-ffb' .or. &
           cstencil(1:9)=='bff-f-ffb' ) then
          ! sixth-order central scheme
          asolver%s(n)%a(-3) =-num1d60
          asolver%s(n)%a(-2) = 0.15d0
          asolver%s(n)%a(-1) =-0.75d0
          
          asolver%s(n)%a( 1) = 0.75d0
          asolver%s(n)%a( 2) =-0.15d0
          asolver%s(n)%a( 3) = num1d60
        elseif(cstencil(2:8)=='ff-f-ff' .or. &
               cstencil(2:8)=='bf-f-ff' .or. &
               cstencil(2:8)=='ff-f-fb' .or. &
               cstencil(2:8)=='bf-f-fb' ) then
          ! fourth-order central scheme
          asolver%s(n)%a(-2) = num1d12
          asolver%s(n)%a(-1) =-num2d3

          asolver%s(n)%a( 1) = num2d3
          asolver%s(n)%a( 2) =-num1d12
        elseif(cstencil(3:7)=='f-f-f' .or. &
               cstencil(3:7)=='b-f-f' .or. &
               cstencil(3:7)=='f-f-b' .or. &
               cstencil(3:7)=='b-f-b'  ) then
         ! second-order central scheme
          asolver%s(n)%a(-1) =-0.5d0

          asolver%s(n)%a( 1) = 0.5d0
        elseif(cstencil(4:7)=='-f-f' .or. &
               cstencil(4:7)=='-b-f'  ) then
          ! second-order biased scheme, i- direction
          asolver%s(n)%a( 0) =-1.5d0
          asolver%s(n)%a( 1) = 2.d0
          asolver%s(n)%a( 2) =-0.5d0
        elseif(cstencil(3:6)=='f-f-' .or. &
               cstencil(3:6)=='f-b-' ) then
          ! second-order biased scheme, i+ direction
          asolver%s(n)%a(-2) =+0.5d0
          asolver%s(n)%a(-1) =-2.d0
          asolver%s(n)%a( 0) = 1.5d0
        elseif(cstencil(1:9)=='bbb-b-bbb') then
          asolver%s(n)%a(-3) =-num1d60
          asolver%s(n)%a(-2) = 0.15d0
          asolver%s(n)%a(-1) =-0.75d0
          
          asolver%s(n)%a( 1) = 0.75d0
          asolver%s(n)%a( 2) =-0.15d0
          asolver%s(n)%a( 3) = num1d60
        elseif(cstencil(2:8)=='bb-b-bb') then
          asolver%s(n)%a(-2) = num1d12
          asolver%s(n)%a(-1) =-num2d3

          asolver%s(n)%a( 1) = num2d3
          asolver%s(n)%a( 2) =-num1d12
        elseif(cstencil(3:7)=='b-b-b') then
          asolver%s(n)%a(-1) =-0.5d0

          asolver%s(n)%a( 1) = 0.5d0
        elseif(cstencil(4:7)=='-b-b') then
          asolver%s(n)%a( 0) =-1.5d0
          asolver%s(n)%a( 1) = 2.d0
          asolver%s(n)%a( 2) =-0.5d0
        elseif(cstencil(3:6)=='b-b-') then
          asolver%s(n)%a(-2) =+0.5d0
          asolver%s(n)%a(-1) =-2.d0
          asolver%s(n)%a( 0) = 1.5d0
        endif

      enddo

    end subroutine fdm_solver_init2

    subroutine fdm_solver_init(asolver,dim,nhalo,node_state)

      class(fdm_type),target :: asolver
      integer,intent(in) :: dim,nhalo

      character(len=1),intent(in) :: node_state(-nhalo:dim+nhalo)

      integer :: n
      character(len=9) :: cstencil
      
      asolver%dim=dim
      allocate(asolver%s(0:dim))

      do n=0,dim
        cstencil=node_state(n-3)//node_state(n-2)//node_state(n-1)// &
                             '-'//node_state(n)//'-'//               &
                 node_state(n+1)//node_state(n+2)//node_state(n+3)

        allocate(asolver%s(n)%a(-3:3))
        asolver%s(n)%a=0.d0
        
        if(cstencil(1:9)=='fff-f-fff' .or. &
           cstencil(1:9)=='bff-f-fff' .or. &
           cstencil(1:9)=='fff-f-ffb' .or. &
           cstencil(1:9)=='bff-f-ffb' ) then
          ! sixth-order central scheme
          asolver%s(n)%a(-3) =-num1d60
          asolver%s(n)%a(-2) = 0.15d0
          asolver%s(n)%a(-1) =-0.75d0
          
          asolver%s(n)%a( 1) = 0.75d0
          asolver%s(n)%a( 2) =-0.15d0
          asolver%s(n)%a( 3) = num1d60
        elseif(cstencil(2:8)=='ff-f-ff' .or. &
               cstencil(2:8)=='bf-f-ff' .or. &
               cstencil(2:8)=='ff-f-fb' .or. &
               cstencil(2:8)=='bf-f-fb' ) then
          ! fourth-order central scheme
          asolver%s(n)%a(-2) = num1d12
          asolver%s(n)%a(-1) =-num2d3

          asolver%s(n)%a( 1) = num2d3
          asolver%s(n)%a( 2) =-num1d12
        elseif(cstencil(3:7)=='f-f-f' .or. &
               cstencil(3:7)=='b-f-f' .or. &
               cstencil(3:7)=='f-f-b' .or. &
               cstencil(3:7)=='b-f-b'  ) then
         ! second-order central scheme
          asolver%s(n)%a(-1) =-0.5d0

          asolver%s(n)%a( 1) = 0.5d0
        elseif(cstencil(4:7)=='-f-f' .or. &
               cstencil(4:7)=='-b-f'  ) then
          ! second-order biased scheme, i- direction
          asolver%s(n)%a( 0) =-1.5d0
          asolver%s(n)%a( 1) = 2.d0
          asolver%s(n)%a( 2) =-0.5d0
        elseif(cstencil(3:6)=='f-f-' .or. &
               cstencil(3:6)=='f-b-' ) then
          ! second-order biased scheme, i+ direction
          asolver%s(n)%a(-2) =+0.5d0
          asolver%s(n)%a(-1) =-2.d0
          asolver%s(n)%a( 0) = 1.5d0
        elseif(cstencil(1:9)=='bbb-b-bbb') then
          asolver%s(n)%a(-3) =-num1d60
          asolver%s(n)%a(-2) = 0.15d0
          asolver%s(n)%a(-1) =-0.75d0
          
          asolver%s(n)%a( 1) = 0.75d0
          asolver%s(n)%a( 2) =-0.15d0
          asolver%s(n)%a( 3) = num1d60
        elseif(cstencil(2:8)=='bb-b-bb') then
          asolver%s(n)%a(-2) = num1d12
          asolver%s(n)%a(-1) =-num2d3

          asolver%s(n)%a( 1) = num2d3
          asolver%s(n)%a( 2) =-num1d12
        elseif(cstencil(3:7)=='b-b-b') then
          asolver%s(n)%a(-1) =-0.5d0

          asolver%s(n)%a( 1) = 0.5d0
        elseif(cstencil(4:7)=='-b-b') then
          asolver%s(n)%a( 0) =-1.5d0
          asolver%s(n)%a( 1) = 2.d0
          asolver%s(n)%a( 2) =-0.5d0
        elseif(cstencil(3:6)=='b-b-') then
          asolver%s(n)%a(-2) =+0.5d0
          asolver%s(n)%a(-1) =-2.d0
          asolver%s(n)%a( 0) = 1.5d0
        endif

      enddo

    end subroutine fdm_solver_init

    subroutine fdm_solver_print(asolver)
      class(fdm_type),target :: asolver
      integer :: n

      do n=0,size(asolver%s)-1
        print*,asolver%s(n)%a(:)
      enddo

    end subroutine fdm_solver_print

end module pastr_fdm