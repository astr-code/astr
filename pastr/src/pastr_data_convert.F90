module pastr_data_convert

    use iso_fortran_env, only: wp => real64

    implicit none


    private


    public :: dat_out_cal

    
    Interface dat_out_cal
      module procedure dat_out_cal_2d
      module procedure dat_out_cal_3d
    end Interface dat_out_cal

contains

    function dat_out_cal_2d(data_in,name_in,name_out,x,y) result(data_out)
      
      use pastr_commvar, only : im,jm,km,const2
      use pastr_flowvis, only : schlieren_xy
      use pastr_gradients, only : grad_xy
      use pastr_thermo_phys, only : sos

      real(wp),intent(in) :: data_in(:,:,:)
      character(len=6),intent(in) :: name_in(:),name_out(:)
      real(wp),intent(in),optional :: x(:,:),y(:,:)
      real(wp) :: data_out(0:im,0:jm,size(name_out))

      real(wp),allocatable :: var1(:,:),var2(:,:),dvar1(:,:,:),  &
                              dvar2(:,:,:),du1(:,:,:),du2(:,:,:)

      integer :: i,j,i1,j1,k,nvar_in,nvar_out
      integer,allocatable :: is(:)

      nvar_in =size(name_in)
      nvar_out=size(name_out)

      do i=1,nvar_out

        do j=1,nvar_in
          
          if(name_out(i)==name_in(j)) then
            data_out(:,:,i)=data_in(:,:,j)
          endif

        enddo

        if(name_out(i)=='Ma') then
          allocate(is(3))
          is=0
          do j=1,nvar_in
            if(name_in(j)=='t')  is(1)=j
            if(name_in(j)=='u1') is(2)=j
            if(name_in(j)=='u2') is(3)=j
          enddo

          allocate(var1(0:im,0:jm))
          do j1=0,jm
          do i1=0,im
            var1(i1,j1)=sos(data_in(i1,j1,is(1)))
          enddo
          enddo
          data_out(:,:,i)=sqrt(data_in(:,:,is(2))**2+data_in(:,:,is(3))**2)/var1

          deallocate(var1,is)
        endif

        if(name_out(i)=='div') then

          allocate(is(3))
          is=0
          do j=1,nvar_in
            if(name_in(j)=='dudx') is(1)=j
            if(name_in(j)=='dvdy') is(2)=j
            if(name_in(j)=='dwdz') is(3)=j
          enddo

          if(is(1)>0) then
            data_out(:,:,i)=data_in(:,:,is(1))+data_in(:,:,is(2))+data_in(:,:,is(3))
          else
            do j=1,nvar_in
              if(name_in(j)=='u1') is(1)=j
              if(name_in(j)=='u2') is(2)=j
            enddo

            if(.not. allocated(du1)) then
              allocate(du1(1:2,0:im,0:jm),du2(1:2,0:im,0:jm))
              du1=grad_xy(data_in(:,:,is(1)),x,y)
              du2=grad_xy(data_in(:,:,is(2)))
            endif

            data_out(:,:,i)=du1(1,:,:)+du2(2,:,:)

          endif

          
          deallocate(is)

        endif

        if(name_out(i)=='omegaz') then

          allocate(is(2))
          is=0
          do j=1,nvar_in
            if(name_in(j)=='dudy') is(1)=j
            if(name_in(j)=='dvdx') is(2)=j
          enddo
          if(is(1)>0) then
            data_out(:,:,i)=-data_in(:,:,is(1))+data_in(:,:,is(2))
          else
           do j=1,nvar_in
              if(name_in(j)=='u1') is(1)=j
              if(name_in(j)=='u2') is(2)=j
            enddo

            if(.not. allocated(du1)) then
              allocate(du1(1:2,0:im,0:jm),du2(1:2,0:im,0:jm))
              du1=grad_xy(data_in(:,:,is(1)),x,y)
              du2=grad_xy(data_in(:,:,is(2)))
            endif

            data_out(:,:,i)=-du1(2,:,:)+du2(1,:,:)
          endif

          deallocate(is)

        endif

        if(name_out(i)=='rns') then

          allocate(is(3))
          is=0
          do j=1,nvar_in
            if(name_in(j)=='ro') is(1)=j
          enddo

          if(is(1)>0) then
            data_out(:,:,i)=schlieren_xy(data_in(:,:,is(1)),x,y)
          else
            do j=1,nvar_in
              if(name_in(j)=='p') is(2)=j
              if(name_in(j)=='t') is(3)=j
            enddo
            allocate(var1(0:im,0:jm))
            var1=data_in(:,:,is(2))/data_in(:,:,is(3))*const2

            data_out(:,:,i)=schlieren_xy(var1,x,y)
          endif
          deallocate(var1)
          deallocate(is)

        endif

      enddo

    end function dat_out_cal_2d

    function dat_out_cal_3d(data_in,name_in,name_out,x,y,z) result(data_out)
      
      use pastr_commvar, only : im,jm,km,const2
      use pastr_flowvis, only : schlieren_xy
      use pastr_gradients, only : grad_xy
      use pastr_thermo_phys, only : sos

      real(wp),intent(in) :: data_in(:,:,:,:)
      character(len=6),intent(in) :: name_in(:),name_out(:)
      real(wp),intent(in),optional :: x(:,:,:),y(:,:,:),z(:,:,:)
      real(wp) :: data_out(0:im,0:jm,0:km,size(name_out))

      real(wp),allocatable :: var1(:,:,:),var2(:,:,:),dvar1(:,:,:,:),  &
                              dvar2(:,:,:,:),du1(:,:,:,:),du2(:,:,:,:),du3(:,:,:,:)

      integer :: i,j,k,i1,j1,k1,nvar_in,nvar_out
      integer,allocatable :: is(:)

      nvar_in =size(name_in)
      nvar_out=size(name_out)

      do i=1,nvar_out

        do j=1,nvar_in

          if(name_out(i)==name_in(j)) then
            data_out(:,:,:,i)=data_in(:,:,:,j)
          endif

        enddo

      enddo

    end function dat_out_cal_3d

end module pastr_data_convert