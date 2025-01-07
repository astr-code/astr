module schemes
  
  use constdef
  use commvar, only : hm

  implicit none

  integer :: len_sten_max=3
  real(8) :: dcoe(5,5)
  
    !   dcoe(1,1) =  1.d0/2.d0
  
    !   dcoe(1,2) =  2.d0/3.d0
    !   dcoe(2,2) = -1.d0/12.d0
  
    !   dcoe(1,3) =  3.d0/4.d0
    !   dcoe(2,3) = -3.d0/20.d0
    !   dcoe(3,3) =  1.d0/60.d0
  
    !   dcoe(1,4) =  4.d0/5.d0
    !   dcoe(2,4) = -1.d0/5.d0
    !   dcoe(3,4) =  4.d0/105.d0
    !   dcoe(4,4) = -1.d0/280.d0 
  
    !   dcoe(1,5) =  5.d0/6.d0
    !   dcoe(2,5) = -5.d0/21.d0
    !   dcoe(3,5) =  5.d0/84.d0
    !   dcoe(4,5) = -5.d0/504.d0 
    !   dcoe(5,5) =  1.d0/1260.d0 

  contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! This subroutine is used to initial for solving the tridiagonal 
  ! martix with two layer of boundary scheme:
  ! A*x=b
  !   |1,c10,.........................|
  !   |c21,1,c11,.....................|
  !   |..c22,1,c12,...................|
  !   |...............................|
  ! A=|...,c2i,1,c1i..................|
  !   |...............................|
  !   |.........,c(2,m-2),1,c(1,m-2)..|
  !   |...........,c(2,m-1),1,c(1,m-1)|
  !   |.........................,c2m,1|
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function compact_scheme_lhs_init(m,ntype,scheme,addata) result(c)

    integer,intent(in) :: m,ntype
    character(len=4),intent(in) :: scheme
    real(8),intent(in),optional :: addata

    real(8),allocatable :: c(:,:)

    allocate(c(1:4,0:m))
    
    c=0.d0

    if(m==0) return
    
    if(ntype==1) then
      ! the block with boundary at i==0

      if(scheme=='642c') then
        c(1,0)=1.d0
        c(1,1)=0.25d0;       c(2,1)=0.25d0
        
        c(1,2:m-1)=num1d3; c(2,2:m-1)=num1d3
      elseif(scheme=='442c') then
        c(1,0)=1.d0
        c(1,1:m-1)=0.25d0; c(2,1:m-1)=0.25d0
      elseif(scheme=='filt') then
        c(1,0:m-1)=addata
        c(2,1:m)=addata
      endif

    elseif(ntype==2) then
      ! the block with boundary at i==im

      if(scheme=='642c') then
        c(1,1:m-2)=num1d3;  c(2,1:m-2)=num1d3
        c(1,m-1)  =0.25d0;  c(2,m-1)  =0.25d0
        c(2,m)  =1.d0
      elseif(scheme=='442c') then
        c(1,1:m-1)  =0.25d0;c(2,1:m-1)  =0.25d0
        c(2,m)  =1.d0
      elseif(scheme=='filt') then
        c(1,1:m-1)=addata
        c(2,2:m)=addata
      endif

    elseif(ntype==3) then
      ! inner block

      if(scheme=='642c') then
        c(1,1:m-1)=num1d3
        c(2,1:m-1)=num1d3
      elseif(scheme=='442c') then
        c(1,1:m-1)=0.25d0
        c(2,1:m-1)=0.25d0
      elseif(scheme=='filt') then
        c(1,1:m-1)=addata
        c(2,1:m-1)=addata

      endif

    elseif(ntype==4) then
      ! the block with boundary at i==0 .and. i==im

      if(scheme=='642c') then
        c(1,0)=1.d0
        c(1,1)=0.25d0;       c(2,1)=0.25d0
        
        c(1,2:m-2)=num1d3; c(2,2:m-2)=num1d3

        c(1,m-1)  =0.25d0;  c(2,m-1)  =0.25d0
        c(2,m)  =1.d0
      elseif(scheme=='642c') then
        c(1,0)=1.d0
        
        c(1,1:m-1)=0.25d0; c(2,1:m-1)=0.25d0

        c(2,m)  =1.d0
      elseif(scheme=='filt') then
        c(1,0:m-1)=addata
        c(2,1:m)  =addata
      endif
      !
    else
      print*,'ntype=',ntype
      print*, ' !! error in subroutine ptds_ini_com !'
      stop
    end if

  end function compact_scheme_lhs_init

  function compact_scheme_rhs(vin,m,hm,ntype,scheme) result(d)

    integer,intent(in) :: m,hm,ntype
    real(8),intent(in) :: vin(-hm:m+hm)
    character(len=4),intent(in) :: scheme
    real(8) :: d(0:m)

    ! local data
    integer :: i
    
    if(scheme=='642c') then
      do i=1,m-1
        d(i)=num7d9*(vin(i+1)-vin(i-1))+ &
            num1d36*(vin(i+2)-vin(i-2))
      end do
    endif

    if(scheme=='442c') then
      do i=1,m-1
        d(i)=0.75d0*(vin(i+1)-vin(i-1))
      end do
    endif

    if(ntype==1) then
      ! the block with boundary at i==0

      if(scheme=='642c') then
        
        d(0)  =2.d0*  (-vin(0)+vin(1))
        d(1)  =0.75d0*( vin(2)-vin(0))

        d(m)=0.75d0* (vin(m+1)-vin(m-1))-  &
             0.15d0* (vin(m+2)-vin(m-2))+  &
             num1d60*(vin(m+3)-vin(m-3))

      
      elseif(scheme=='442c') then
        
        d(0)  =2.d0*  (-vin(0)+vin(1))

        d(m)= num2d3* (vin(m+1)-vin(m-1))-  &
             num1d12* (vin(m+2)-vin(m-2))

      endif

    elseif(ntype==2) then
      ! the block with boundary at i==im

      if(scheme=='642c') then

        d(0)=0.75d0* (vin(1)-vin(-1)) -      &
             0.15d0* (vin(2)-vin(-2))+      &
             num1d60*(vin(3)-vin(-3))

        d(m-1)=0.75d0*( vin(m)-vin(m-2))
        d(m)  =  2.d0*(-vin(m-1)+vin(m))

      elseif(scheme=='442c') then

        d(0)= num2d3* (vin(1)-vin(-1))-  &
             num1d12* (vin(2)-vin(-2))

        d(m)  =  2.d0*(-vin(m-1)+vin(m))


      endif

    elseif(ntype==3) then
      ! inner block

      if(scheme=='642c') then

        d(0)=0.75d0* (vin(1)-vin(-1)) -      &
             0.15d0* (vin(2)-vin(-2))+      &
             num1d60*(vin(3)-vin(-3))

        d(m)=0.75d0* (vin(m+1)-vin(m-1))-  &
             0.15d0* (vin(m+2)-vin(m-2))+  &
             num1d60*(vin(m+3)-vin(m-3))

      elseif(scheme=='442c') then

        d(0)= num2d3* (vin(1)-vin(-1))-  &
             num1d12* (vin(2)-vin(-2))

        d(m)=num2d3* (vin(m+1)-vin(m-1))-  &
             num1d12* (vin(m+2)-vin(m-2))

      endif

    elseif(ntype==4) then
      ! the block with boundary at i==0 .and. i==im

      if(scheme=='642c') then

        d(0)  =2.d0*  (-vin(0)+vin(1))
        d(1)  =0.75d0*( vin(2)-vin(0))
        d(m-1)=0.75d0*( vin(m)-vin(m-2))
        d(m)  =  2.d0*(-vin(m-1)+vin(m))
      elseif(scheme=='442c') then

        d(0)  =2.d0*  (-vin(0)+vin(1))
        d(m)  =  2.d0*(-vin(m-1)+vin(m))
      endif
      !
    else
      print*,'ntype=',ntype
      print*, ' !! error in subroutine ptds_ini_com !'
      stop
    end if

  end function compact_scheme_rhs

  !----------------------------------------
  ! A*x=b
  !   |1,c10,.........................|
  !   |c21,1,c11,.....................|
  !   |..c22,1,c12,...................|
  !   |...............................|
  ! A=|...,c2i,1,c1i..................|
  !   |...............................|
  !   |.........,c(2,m-2),1,c(1,m-2)..|
  !   |...........,c(2,m-1),1,c(1,m-1)|
  !   |.........................,c2m,1|
  !----------------------------------------

  function compact_flux_lhs_init(m,ntype,scheme,addata) result(c)

    integer,intent(in) :: m,ntype
    character(len=4),intent(in) :: scheme
    real(8),intent(in),optional :: addata

    real(8),allocatable :: c(:,:)

    allocate(c(1:4,-1:m))
    
    c=0.d0

    if(m==0) return
    
    if(ntype==1) then
      ! the block with boundary at i==0

      if(scheme=='642c') then
        ! c(1,0) = 0.25d0
        ! c(2,0) = 0.25d0

        c(1,1:m-1)=num1d3
        c(2,1:m-1)=num1d3

      endif

      if(scheme=='442c') then
        c(1,1:m-1)=0.25d0
        c(2,1:m-1)=0.25d0
      endif

    elseif(ntype==2) then
      ! the block with boundary at i==im

      if(scheme=='642c') then
        c(1,0:m-2)=num1d3
        c(2,0:m-2)=num1d3
      endif

      if(scheme=='442c') then
        c(1,0:m-2)=0.25d0
        c(2,0:m-2)=0.25d0
      endif

    elseif(ntype==3) then
      ! inner block

      if(scheme=='642c') then
        c(1,0:m-1)=num1d3
        c(2,0:m-1)=num1d3
      endif

      if(scheme=='442c') then
        c(1,0:m-1)=0.25d0
        c(2,0:m-1)=0.25d0
      endif

    elseif(ntype==4) then
      ! the block with boundary at i==0 .and. i==im
      stop '44444'
      ! if(scheme=='642c') then

      !   c(1,1:m-2)=num1d3
      !   c(2,1:m-2)=num1d3
  
      ! endif
      ! if(scheme=='442c') then

      !   c(1,1:m-2)==0.25d0
      !   c(2,1:m-2)==0.25d0
  
      ! endif

    else
      print*,'ntype=',ntype
      print*, ' !! error in subroutine compact_flux_lhs_init !'
      stop
    end if

  end function compact_flux_lhs_init

  function flux_rhs_642e(vin,m,hm,ntype) result(d)

    integer,intent(in) :: m,hm,ntype
    real(8),intent(in) :: vin(0:len_sten_max,-hm:m+hm)
    real(8) :: d(-1:m)

    ! local data
    integer :: i

    select case (ntype)
       case (1)
       ! boundary at i==0, interface at i==im
          i=-1
            ! v-1=1.5*v0-0.5v1
            d(i)= vin(1,-1)+0.5d0*vin(1,0)-vin(2,-1)
          i=0
            d(i)= 0.5d0*vin(1,i) 
          
          i=1
            d(i)=  num2d3*vin(1,i)-num1d12*(vin(2,i)+vin(2,i-1))
          
          do i=2,m
            d(i)=  num3d4*vin(1,i) - &
                  num3d20*(vin(2,i)+vin(2,i-1)) + &
                  num1d60*(vin(3,i)+vin(3,i-1)+vin(3,i-2))
          enddo
       case (2)
       ! boundary at i==im, interface at i==0
          do i=-1,m-3
            d(i)=  num3d4*vin(1,i) - &
                  num3d20*(vin(2,i)+vin(2,i-1)) + &
                  num1d60*(vin(3,i)+vin(3,i-1)+vin(3,i-2))
          enddo
          i=m-2
            d(i)=  num2d3*vin(1,i)-num1d12*(vin(2,i)+vin(2,i-1))
          
          i=m-1
            d(i)= 0.5d0*vin(1,i) 

          i=m
            d(i)= vin(1,m)+0.5d0*vin(1,m-1)-vin(2,m-1)

       case (3)
       ! interface at i==im and i==0
          do i=-1,m
            d(i)=  num3d4*vin(1,i) - &
                  num3d20*(vin(2,i)+vin(2,i-1)) + &
                  num1d60*(vin(3,i)+vin(3,i-1)+vin(3,i-2))
          enddo
       case (4)
       ! boundary at i==im and i==0
          i=-1
            ! v-1=2*v0-v1
            d(i)= vin(1,-1)+0.5d0*vin(1,0)-vin(2,-1)
          i=0
            d(i)= 0.5d0*vin(1,i) 
          
          i=1
            d(i)=  num2d3*vin(1,i)-num1d12*(vin(2,i)+vin(2,i-1))
          do i=2,m-3
            d(i)=  num3d4*vin(1,i) - &
                  num3d20*(vin(2,i)+vin(2,i-1)) + &
                  num1d60*(vin(3,i)+vin(3,i-1)+vin(3,i-2))
          enddo
          i=m-2
            d(i)=  num2d3*vin(1,i)-num1d12*(vin(2,i)+vin(2,i-1))
          
          i=m-1
            d(i)= 0.5d0*vin(1,i) 

          i=m
            d(i)= vin(1,m)+0.5d0*vin(1,m-1)-vin(2,m-1)
       case default
         print*, ' !! ntype error',ntype,' @ flux_rhs_642e !'
         stop
    end select

    return

  end function flux_rhs_642e

  function flux_rhs_442e(vin,m,hm,ntype) result(d)

    integer,intent(in) :: m,hm,ntype
    real(8),intent(in) :: vin(0:len_sten_max,-hm:m+hm)
    real(8) :: d(-1:m)

    ! local data
    integer :: i

    select case (ntype)
       case (1)
       ! boundary at i==0, interface at i==im
          i=-1
            ! v-1=1.5*v0-0.5v1
            d(i)= vin(1,-1)+0.5d0*vin(1,0)-vin(2,-1)
          i=0
            d(i)= 0.5d0*vin(1,i) 
          
          do i=1,m
            d(i)=  num2d3*vin(1,i)-num1d12*(vin(2,i)+vin(2,i-1))
          enddo
       case (2)
       ! boundary at i==im, interface at i==0
          do i=-1,m-2
            d(i)=  num2d3*vin(1,i)-num1d12*(vin(2,i)+vin(2,i-1))
          enddo
          i=m-1
            d(i)= 0.5d0*vin(1,i) 

          i=m
            d(i)= vin(1,m)+0.5d0*vin(1,m-1)-vin(2,m-1)

       case (3)
       ! interface at i==im and i==0
          do i=-1,m
            d(i)=  num2d3*vin(1,i)-num1d12*(vin(2,i)+vin(2,i-1))
          enddo
       case (4)
       ! boundary at i==im and i==0
          i=-1
            ! v-1=2*v0-v1
            d(i)= vin(1,-1)+0.5d0*vin(1,0)-vin(2,-1)
          i=0
            d(i)= 0.5d0*vin(1,i) 

          do i=1,m-2
            d(i)=  num2d3*vin(1,i)-num1d12*(vin(2,i)+vin(2,i-1))
          enddo
        
          i=m-1
            d(i)= 0.5d0*vin(1,i) 

          i=m
            d(i)= vin(1,m)+0.5d0*vin(1,m-1)-vin(2,m-1)
       case default
         print*, ' !! ntype error',ntype,' @ flux_rhs_442e !'
         stop
    end select

    return

  end function flux_rhs_442e

  function flux_rhs_642c(vin,m,hm,ntype) result(d)

    integer,intent(in) :: m,hm,ntype
    real(8),intent(in) :: vin(0:len_sten_max,-hm:m+hm)
    real(8) :: d(-1:m)

    ! local data
    integer :: i

    if(ntype==1) then
      ! the block with boundary at i==0

      i=-1
        ! v-1=1.5*v0-0.5v1
        d(i)= vin(1,-1)+0.5d0*vin(1,0)-vin(2,-1)
      i=0
        ! d(i)= 0.5d0*vin(1,i) 
        d(i)= 1.25d0*vin(1,i)-0.5d0*vin(0,i)

      do i=1,m-1
        d(i)=num7d9*vin(1,i)+num1d36*(vin(2,i)+vin(2,i-1))
      end do
      i=m
      d(i)=  num3d4*vin(1,i) - &
            num3d20*(vin(2,i)+vin(2,i-1)) + &
            num1d60*(vin(3,i)+vin(3,i-1)+vin(3,i-2))

    elseif(ntype==2) then
      ! the block with boundary at i==im
      
      i=-1
      d(i)=  num3d4*vin(1,i) - &
            num3d20*(vin(2,i)+vin(2,i-1)) + &
            num1d60*(vin(3,i)+vin(3,i-1)+vin(3,i-2))
      do i=0,m-2
        d(i)=num7d9*vin(1,i)+num1d36*(vin(2,i)+vin(2,i-1))
      end do
      i=m-1
        ! d(i)= 0.5d0*vin(1,i)
        d(i)=0.25d0*vin(1,i)+0.5d0*vin(0,i)

      i=m
        d(i)= vin(1,m)+0.5d0*vin(1,m-1)-vin(2,m-1)

    elseif(ntype==3) then
      ! inner block

      i=-1
      d(i)=  num3d4*vin(1,i) - &
            num3d20*(vin(2,i)+vin(2,i-1)) + &
            num1d60*(vin(3,i)+vin(3,i-1)+vin(3,i-2))
      do i=0,m-1
        d(i)=num7d9*vin(1,i)+num1d36*(vin(2,i)+vin(2,i-1))
      end do
      i=m
      d(i)=  num3d4*vin(1,i) - &
            num3d20*(vin(2,i)+vin(2,i-1)) + &
            num1d60*(vin(3,i)+vin(3,i-1)+vin(3,i-2))
   
    elseif(ntype==4) then
      ! the block with boundary at i==0 .and. i==im

      i=-1
        ! v-1=1.5*v0-0.5v1
        d(i)= vin(1,-1)+0.5d0*vin(1,0)-vin(2,-1)
      i=0
        ! d(i)= 0.5d0*vin(1,i)
        d(i)= 1.25d0*vin(1,i)-0.5d0*vin(0,i)
      do i=1,m-2
        d(i)=num7d9*vin(1,i)+num1d36*(vin(2,i)+vin(2,i-1))
      end do
      i=m-1
        ! d(i)= 0.5d0*vin(1,i) 
        d(i)=0.25d0*vin(1,i)+0.5d0*vin(0,i)

      i=m
        d(i)= vin(1,m)+0.5d0*vin(1,m-1)-vin(2,m-1)

    else
      print*,'ntype=',ntype
      print*, ' !! error in subroutine flux_rhs_642c !'
      stop
    end if

  end function flux_rhs_642c

  function flux_rhs_442c(vin,m,hm,ntype) result(d)

    integer,intent(in) :: m,hm,ntype
    real(8),intent(in) :: vin(0:len_sten_max,-hm:m+hm)
    real(8) :: d(-1:m)

    ! local data
    integer :: i

    if(ntype==1) then
      ! the block with boundary at i==0

      i=-1
        ! v-1=1.5*v0-0.5v1
        d(i)= vin(1,-1)+0.5d0*vin(1,0)-vin(2,-1)
      i=0
        d(i)= 0.5d0*vin(1,i) 
        ! d(i)= 1.25d0*vin(1,i)-0.5d0*vin(0,i)

      do i=1,m-1
        d(i)=0.75d0*vin(1,i)
      end do
      i=m
      d(i)=  num2d3*vin(1,i) - &
            num1d12*(vin(2,i)+vin(2,i-1))

    elseif(ntype==2) then
      ! the block with boundary at i==im
      
      i=-1
      d(i)=  num2d3*vin(1,i) - &
            num1d12*(vin(2,i)+vin(2,i-1))

      do i=0,m-2
        d(i)=0.75d0*vin(1,i)
      end do
      i=m-1
        d(i)= 0.5d0*vin(1,i)
        ! d(i)=0.25d0*vin(1,i)+0.5d0*vin(0,i)

      i=m
        d(i)= vin(1,m)+0.5d0*vin(1,m-1)-vin(2,m-1)

    elseif(ntype==3) then
      ! inner block

      i=-1
      d(i)=  num2d3*vin(1,i) - &
            num1d12*(vin(2,i)+vin(2,i-1))
      do i=0,m-1
        d(i)=0.75d0*vin(1,i)
      end do
      i=m
      d(i)=  num2d3*vin(1,i) - &
            num1d12*(vin(2,i)+vin(2,i-1))
   
    elseif(ntype==4) then
      ! the block with boundary at i==0 .and. i==im

      stop '4444'
      ! i=-1
      !   ! v-1=1.5*v0-0.5v1
      !   d(i)= vin(1,-1)+0.5d0*vin(1,0)-vin(2,-1)
      ! i=0
      !   ! d(i)= 0.5d0*vin(1,i)
      !   d(i)= 1.25d0*vin(1,i)-0.5d0*vin(0,i)
      ! do i=1,m-2
      !   d(i)=num7d9*vin(1,i)+num1d36*(vin(2,i)+vin(2,i-1))
      ! end do
      ! i=m-1
      !   ! d(i)= 0.5d0*vin(1,i) 
      !   d(i)=0.25d0*vin(1,i)+0.5d0*vin(0,i)

      ! i=m
      !   d(i)= vin(1,m)+0.5d0*vin(1,m-1)-vin(2,m-1)

    else
      print*,'ntype=',ntype
      print*, ' !! error in subroutine flux_rhs_442c !'
      stop
    end if

  end function flux_rhs_442c
  !+-------------------------------------------------------------------+
  !| This function is to return the finit difference value of 6th-order|
  !| central scheme, inculding boundary closure.                       |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 19-02-2021  | Created by J. Fang @ Warrington.                    |
  !+-------------------------------------------------------------------+
  function diff8ec(vin,dim,ns,ntype) result(vout)
    !
    integer,intent(in) :: dim,ns,ntype
    real(8),intent(in) :: vin(-hm:dim+hm)
    real(8) :: vout(0:dim)
    !
    ! local data
    integer :: i
    !
    if(ntype==1) then
      !
      if(ns==642) then
        ! ns==642: 2-4-6-6-6-...-6-6-6-4-2
        vout(0)=-0.5d0*vin(2)+2.d0*vin(1)-1.5d0*vin(0)
        vout(1)=0.5d0*(vin(2)-vin(0))
        vout(2)=num2d3*(vin(3)-vin(1))-num1d12*(vin(4)-vin(0))
      else
        print*,' !! ns=',ns
        stop ' error 1 @ diff6c'
      endif
      !
      do i=3,dim
        vout(i)  =0.75d0 *(vin(i+1)-vin(i-1))-                         &
                  0.15d0 *(vin(i+2)-vin(i-2))+                         &
                  num1d60*(vin(i+3)-vin(i-3))
      enddo
      !
    elseif(ntype==2) then
      !
      do i=0,dim-3
        vout(i)  =0.75d0 *(vin(i+1)-vin(i-1))-                         &
                  0.15d0 *(vin(i+2)-vin(i-2))+                         &
                  num1d60*(vin(i+3)-vin(i-3))
      enddo
      !
      if(ns==642) then
        ! ns==642: 2-4-6-6-6-...-6-6-6-4-2
        vout(dim-2) =num2d3*(vin(dim-1)-vin(dim-3))-                   &
                    num1d12*(vin(dim)  -vin(dim-4))
        vout(dim-1)=0.5d0*(vin(dim)-vin(dim-2))
        vout(dim)  =0.5d0*vin(dim-2)-2.d0*vin(dim-1)+1.5d0*vin(dim)
      else
        print*,' !! ns=',ns
        stop ' error 2 @ diff6c'
      endif
      !
    elseif(ntype==3) then
      do i=0,dim
        vout(i)  =(672.d0*(vin(i+1)-vin(i-1))-                         &
                   168.d0*(vin(i+2)-vin(i-2))+                         &
                    32.d0*(vin(i+3)-vin(i-3))-                         &
                     3.d0*(vin(i+4)-vin(i-4)))*num1d840
      enddo
    else
      print*,' !! ntype=',ntype
      stop ' !! errpr 3 @ diff6c'
    endif
    !
  end function diff8ec
  !
  function diff6ec(vin,dim,ns,ntype) result(vout)
    !
    integer,intent(in) :: dim,ns,ntype
    real(8),intent(in) :: vin(-hm:dim+hm)
    real(8) :: vout(0:dim)
    !
    ! local data
    integer :: i
    !
    if(ntype==1) then
      !
      if(ns==642) then
        ! ns==642: 2-4-6-6-6-...-6-6-6-4-2
        vout(0)=-0.5d0*vin(2)+2.d0*vin(1)-1.5d0*vin(0)
        vout(1)=0.5d0*(vin(2)-vin(0))
        vout(2)=num2d3*(vin(3)-vin(1))-num1d12*(vin(4)-vin(0))
      else
        print*,' !! ns=',ns
        stop ' error 1 @ diff6c'
      endif
      !
      do i=3,dim
        vout(i)  =0.75d0 *(vin(i+1)-vin(i-1))-                         &
                  0.15d0 *(vin(i+2)-vin(i-2))+                         &
                  num1d60*(vin(i+3)-vin(i-3))
      enddo
      !
    elseif(ntype==2) then
      !
      do i=0,dim-3
        vout(i)  =0.75d0 *(vin(i+1)-vin(i-1))-                         &
                  0.15d0 *(vin(i+2)-vin(i-2))+                         &
                  num1d60*(vin(i+3)-vin(i-3))
      enddo
      !
      if(ns==642) then
        ! ns==642: 2-4-6-6-6-...-6-6-6-4-2
        vout(dim-2) =num2d3*(vin(dim-1)-vin(dim-3))-                   &
                    num1d12*(vin(dim)  -vin(dim-4))
        vout(dim-1)=0.5d0*(vin(dim)-vin(dim-2))
        vout(dim)  =0.5d0*vin(dim-2)-2.d0*vin(dim-1)+1.5d0*vin(dim)
      else
        print*,' !! ns=',ns
        stop ' error 2 @ diff6c'
      endif
      !
    elseif(ntype==3) then
      do i=0,dim
        vout(i)  =0.75d0 *(vin(i+1)-vin(i-1))-                         &
                  0.15d0 *(vin(i+2)-vin(i-2))+                         &
                  num1d60*(vin(i+3)-vin(i-3))
      enddo
    else
      print*,' !! ntype=',ntype
      stop ' !! errpr 3 @ diff6c'
    endif
    !
  end function diff6ec
  !
  function diff4ec(vin,dim,ns,ntype) result(vout)
    !
    integer,intent(in) :: dim,ns,ntype
    real(8),intent(in) :: vin(-hm:dim+hm)
    real(8) :: vout(0:dim)
    !
    ! local data
    integer :: i
    !
    if(ntype==1) then
      !
      if(ns==442) then
        ! ns==422: 2-4-4-...-4-4-2
        vout(0)=-0.5d0*vin(2)+2.d0*vin(1)-1.5d0*vin(0)
        vout(1)=0.5d0*(vin(2)-vin(0))
      else
        print*,' !! ns=',ns
        stop ' error 1 @ diff6c'
      endif
      !
      do i=2,dim
        vout(i)  =num2d3*(vin(i+1)-vin(i-1))-                         &
                  num1d12*(vin(i+2)-vin(i-2))
      enddo
      !
    elseif(ntype==2) then
      !
      do i=0,dim-2
        vout(i)  =num2d3*(vin(i+1)-vin(i-1))-                         &
                  num1d12*(vin(i+2)-vin(i-2))
      enddo
      !
      if(ns==442) then
        ! ns==422: 2-2-4-...-4-2-2
        vout(dim-1)=0.5d0*(vin(dim)-vin(dim-2))
        vout(dim)  =0.5d0*vin(dim-2)-2.d0*vin(dim-1)+1.5d0*vin(dim)
      else
        print*,' !! ns=',ns
        stop ' error 2 @ diff6c'
      endif
      !
    elseif(ntype==3) then
      do i=0,dim
        vout(i)  =num2d3*(vin(i+1)-vin(i-1))-                         &
                  num1d12*(vin(i+2)-vin(i-2))
      enddo
    else
      print*,' !! ntype=',ntype
      stop ' !! errpr 3 @ diff6c'
    endif
    !
  end function diff4ec
  !
  function diff2ec(vin,dim,ns,ntype) result(vout)
    !
    integer,intent(in) :: dim,ns,ntype
    real(8),intent(in) :: vin(-hm:dim+hm)
    real(8) :: vout(0:dim)
    !
    ! local data
    integer :: i
    !
    if(ntype==1) then
      !
      ! vout(0)=-0.5d0*vin(2)+2.d0*vin(1)-1.5d0*vin(0)
      vout(0)=vin(1)-vin(0)
      !
      do i=1,dim
        vout(i)  =0.5d0*(vin(i+1)-vin(i-1))
      enddo
      !
    elseif(ntype==2) then
      !
      do i=0,dim-1
        vout(i)  =0.5d0*(vin(i+1)-vin(i-1))
      enddo
      !
      ! vout(dim)  =0.5d0*vin(dim-2)-2.d0*vin(dim-1)+1.5d0*vin(dim)
      vout(dim)  =-vin(dim-1)+vin(dim)
      !
    elseif(ntype==3) then
      do i=0,dim
        vout(i)  =0.5d0*(vin(i+1)-vin(i-1))
      enddo
    else
      print*,' !! ntype=',ntype
      stop ' !! errpr 3 @ diff6c'
    endif
    !
  end function diff2ec
  !+-------------------------------------------------------------------+
  !| The end of the diffcen ddf.                                       |
  !+-------------------------------------------------------------------+


end module schemes