module pastr_data_process

    use iso_fortran_env, only: wp => real64

    implicit none

    real(wp),allocatable,dimension(:,:,:) :: ro,u1,u2,u3,p,t
    real(wp),allocatable,dimension(:,:) :: ro_xy,u1_xy,u2_xy,u3_xy,p_xy,t_xy

contains

    subroutine stats_read_process(mode)

      use pastr_commvar, only: im,jm,km
      use pastr_io, only : read_stats
      
      character(len=*),intent(in) :: mode

      if(mode=='meanflow' .or. mode=='all') then

        allocate( ro(0:im,0:jm,0:km),u1(0:im,0:jm,0:km),  &
                  u2(0:im,0:jm,0:km),u3(0:im,0:jm,0:km),  &
                  p(0:im,0:jm,0:km),t(0:im,0:jm,0:km) )

        call read_stats(var=ro,varname='rom')
        call read_stats(var=u1,varname='u1m')
        call read_stats(var=u2,varname='u2m')
        call read_stats(var=u3,varname='u3m')
        call read_stats(var= p,varname= 'pm')
        call read_stats(var= t,varname= 'tm')

        allocate( ro_xy(0:im,0:jm),u1_xy(0:im,0:jm),  &
                  u2_xy(0:im,0:jm),u3_xy(0:im,0:jm),  &
                   p_xy(0:im,0:jm), t_xy(0:im,0:jm) )

        ro_xy=average(ro)
        u1_xy=average(u1)/ro_xy
        u2_xy=average(u2)/ro_xy
        u3_xy=average(u3)/ro_xy
         p_xy=average(p)
         t_xy=average(t)/ro_xy

        call write_stats()

        call plot_meanflow_xy()

      elseif(mode=='2ndsta' .or. mode=='all') then
      elseif(mode=='3rdsta' .or. mode=='all') then
      elseif(mode=='budget' .or. mode=='all') then
      else
       stop 1
      endif

    end subroutine stats_read_process

    subroutine write_stats

      use pastr_commtype, only : tblock
      use pastr_commvar, only: im,jm,km
      use pastr_multiblock_type, only: block_define
      use pastr_h5io

      type(tblock),allocatable,target :: pblocks(:)
      integer :: nblocks
      logical :: multi_block
      integer :: i
      type(tblock),pointer :: b

      call block_define(multi_block,nblocks,pblocks)

      do i=1,nblocks
          b=>pblocks(i)
          b%nvar=6
          write(b%name,'(A,I5.5)')'b',i
          call b%init_data()
          b%var(0:b%im,0:b%jm,0,1)=ro_xy(b%ilo:b%ihi,b%jlo:b%jhi)
          b%var(0:b%im,0:b%jm,0,2)=u1_xy(b%ilo:b%ihi,b%jlo:b%jhi)
          b%var(0:b%im,0:b%jm,0,3)=u2_xy(b%ilo:b%ihi,b%jlo:b%jhi)
          b%var(0:b%im,0:b%jm,0,4)=u3_xy(b%ilo:b%ihi,b%jlo:b%jhi)
          b%var(0:b%im,0:b%jm,0,5)= p_xy(b%ilo:b%ihi,b%jlo:b%jhi)
          b%var(0:b%im,0:b%jm,0,6)= t_xy(b%ilo:b%ihi,b%jlo:b%jhi)
          b%varname(1) = 'ro'
          b%varname(2) = 'u'
          b%varname(3) = 'v'
          b%varname(4) = 'w'
          b%varname(5) = 'p'
          b%varname(6) = 't'
      enddo

      call h5_write_blocks(blocks=pblocks,fname='Results/mean.fav.zm.h5')

      ! call H5WriteArray(ro_xy,im,jm,'ro','Results/mean.fav.zm.h5')
      ! call H5WriteArray(u1_xy,im,jm,'u1','Results/mean.fav.zm.h5')
      ! call H5WriteArray(u2_xy,im,jm,'u2','Results/mean.fav.zm.h5')
      ! call H5WriteArray(u3_xy,im,jm,'u3','Results/mean.fav.zm.h5')
      ! call H5WriteArray( p_xy,im,jm, 'p','Results/mean.fav.zm.h5')
      ! call H5WriteArray( t_xy,im,jm, 't','Results/mean.fav.zm.h5')

    end subroutine write_stats

    subroutine plot_meanflow_xy

      use pastr_commvar, only: im,jm,km
      use pastr_io, only : read_grid
      use pastr_multiblock_type, only: block_define
      use pastr_tecio
      use pastr_commtype, only : tblock

      real(wp),allocatable :: x(:,:),y(:,:),z(:,:)
      type(tblock),allocatable,target :: pblocks(:)
      integer :: nblocks
      logical :: multi_block
      integer :: i
      type(tblock),pointer :: b

      allocate(x(0:im,0:jm),y(0:im,0:jm),z(0:im,0:jm))

      call read_grid(x=x,y=y,z=z,kslice=0)

      call block_define(multi_block,nblocks,pblocks)

      do i=1,nblocks
          b=>pblocks(i)
          b%nvar=6
          call b%init_data()
          b%x(0:b%im,0:b%jm,0)=x(b%ilo:b%ihi,b%jlo:b%jhi)
          b%y(0:b%im,0:b%jm,0)=y(b%ilo:b%ihi,b%jlo:b%jhi)
          b%var(0:b%im,0:b%jm,0,1)=ro_xy(b%ilo:b%ihi,b%jlo:b%jhi)
          b%var(0:b%im,0:b%jm,0,2)=u1_xy(b%ilo:b%ihi,b%jlo:b%jhi)
          b%var(0:b%im,0:b%jm,0,3)=u2_xy(b%ilo:b%ihi,b%jlo:b%jhi)
          b%var(0:b%im,0:b%jm,0,4)=u3_xy(b%ilo:b%ihi,b%jlo:b%jhi)
          b%var(0:b%im,0:b%jm,0,5)= p_xy(b%ilo:b%ihi,b%jlo:b%jhi)
          b%var(0:b%im,0:b%jm,0,6)= t_xy(b%ilo:b%ihi,b%jlo:b%jhi)
          b%varname(1) = 'ro'
          b%varname(2) = 'u'
          b%varname(3) = 'v'
          b%varname(4) = 'w'
          b%varname(5) = 'p'
          b%varname(6) = 't'

      enddo

      call tecbin(filename='Results/tecmean',block=pblocks)

    end subroutine plot_meanflow_xy

    function average(var) result(varm)

      use pastr_commvar, only: im,jm,km

      real(wp),intent(in) :: var(0:im,0:jm,0:km)
      real(wp) :: varm(0:im,0:jm)

      integer :: k

      varm=0._wp
      do k=1,km
        varm(:,:)=varm(:,:)+var(:,:,k)
      enddo
      varm=varm/dble(km)

    end function average

end module pastr_data_process