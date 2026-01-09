module pastr_data_process

    use iso_fortran_env, only: wp => real64
    use pastr_commtype, only : tblock

    implicit none


contains

    subroutine stats_read_process(mode)

      use pastr_commvar, only: im,jm,km
      use pastr_io, only : read_stats
      use pastr_h5io
      use pastr_multiblock_type, only: block_define
      
      character(len=*),intent(in) :: mode

      integer :: nblocks
      logical :: multi_block
      type(tblock),allocatable,target :: pblocks(:)
      real(wp),allocatable :: var(:,:,:),varxy(:,:,:)
      type(tblock),pointer :: b,c
      integer :: i
      character(len=32), allocatable :: names(:)

      type(tblock),allocatable,target :: block_mean(:),block_2nd(:), &
                                         block_3rd(:)

      call block_define(multi_block,nblocks,pblocks)

      do i=1,size(pblocks)
        call pblocks(i)%patch_info
      enddo

      stop

      if(mode=='meanflow' .or. mode=='all') then

        allocate(var(0:im,0:jm,0:km))
        allocate( varxy(0:im,0:jm,1:6) )

        var=read_stats(varname='rom',filename='outdat/meanflow.h5')
        varxy(:,:,1)=average_k(var)
        var=read_stats(varname='u1m',filename='outdat/meanflow.h5')
        varxy(:,:,2)=average_k(var)/varxy(:,:,1)
        var=read_stats(varname='u2m',filename='outdat/meanflow.h5')
        varxy(:,:,3)=average_k(var)/varxy(:,:,1)
        var=read_stats(varname='u3m',filename='outdat/meanflow.h5')
        varxy(:,:,4)=average_k(var)/varxy(:,:,1)
        var=read_stats(varname= 'pm',filename='outdat/meanflow.h5')
        varxy(:,:,5)=average_k(var)
        var=read_stats(varname= 'tm',filename='outdat/meanflow.h5')
        varxy(:,:,6)=average_k(var)/varxy(:,:,1)

        allocate(block_mean, source=pblocks)

        do i=1,nblocks
          b=>block_mean(i)
          b%nvar=6
          call b%init_data()
          b%var(0:b%im,0:b%jm,0,1:6)=varxy(b%ilo:b%ihi,b%jlo:b%jhi,1:6)
          b%varname(1) = 'ro'
          b%varname(2) = 'u'
          b%varname(3) = 'v'
          b%varname(4) = 'w'
          b%varname(5) = 'p'
          b%varname(6) = 't'
        enddo

        call h5_write_blocks(blocks=block_mean,fname='Results/mean.fav.zm.h5')

        call plot_block_xy(blocks=block_mean,filename='Results/tecmean')

        deallocate(varxy,var)
      endif
      
      if(mode=='2ndsta' .or. mode=='all') then

        if(allocated(block_mean)) deallocate(block_mean)

        allocate(block_mean, source=pblocks)
        allocate(names(6))
        names(1)='ro'
        names(2)='u'
        names(3)='v'
        names(4)='w'
        names(5)='t'
        names(6)='p'

        call h5_read_blocks(blocks=block_mean,fname='Results/mean.fav.zm.h5',namelist=names)

        allocate( var(0:im,0:jm,0:km))
        allocate( varxy(0:im,0:jm,1:11) )        

        var=read_stats(varname='tt',filename='outdat/2ndsta.h5')
        varxy(:,:,1)=average_k(var)
        var=read_stats(varname='tu1',filename='outdat/2ndsta.h5')
        varxy(:,:,2)=average_k(var)
        var=read_stats(varname='tu2',filename='outdat/2ndsta.h5')
        varxy(:,:,3)=average_k(var)
        var=read_stats(varname='tu3',filename='outdat/2ndsta.h5')
        varxy(:,:,4)=average_k(var)
        var=read_stats(varname='u11',filename='outdat/2ndsta.h5')
        varxy(:,:,5)=average_k(var)
        var=read_stats(varname='u12',filename='outdat/2ndsta.h5')
        varxy(:,:,6)=average_k(var)
        var=read_stats(varname='u13',filename='outdat/2ndsta.h5')
        varxy(:,:,7)=average_k(var)
        var=read_stats(varname='u22',filename='outdat/2ndsta.h5')
        varxy(:,:,8)=average_k(var)
        var=read_stats(varname='u23',filename='outdat/2ndsta.h5')
        varxy(:,:,9)=average_k(var)
        var=read_stats(varname='u33',filename='outdat/2ndsta.h5')
        varxy(:,:,10)=average_k(var)
        var=read_stats(varname='pp',filename='outdat/2ndsta.h5')
        varxy(:,:,11)=average_k(var)

        allocate(block_2nd, source=pblocks)

        do i=1,nblocks
          b=>block_2nd(i)
          c=>block_mean(i)
          b%nvar=11
          call b%init_data()
          b%var(0:b%im,0:b%jm,0,1:11)=varxy(b%ilo:b%ihi,b%jlo:b%jhi,1:11)

          b%var(:,:,:,1) =b%var(:,:,:,1) /c%var(:,:,:,1)-c%var(:,:,:,5)*c%var(:,:,:,5)
          b%var(:,:,:,2) =b%var(:,:,:,2) /c%var(:,:,:,1)-c%var(:,:,:,2)*c%var(:,:,:,5)
          b%var(:,:,:,3) =b%var(:,:,:,3) /c%var(:,:,:,1)-c%var(:,:,:,3)*c%var(:,:,:,5)
          b%var(:,:,:,4) =b%var(:,:,:,4) /c%var(:,:,:,1)-c%var(:,:,:,4)*c%var(:,:,:,5)
          b%var(:,:,:,5) =b%var(:,:,:,5) /c%var(:,:,:,1)-c%var(:,:,:,2)*c%var(:,:,:,2)
          b%var(:,:,:,6) =b%var(:,:,:,6) /c%var(:,:,:,1)-c%var(:,:,:,2)*c%var(:,:,:,3)
          b%var(:,:,:,7) =b%var(:,:,:,7) /c%var(:,:,:,1)-c%var(:,:,:,2)*c%var(:,:,:,4)
          b%var(:,:,:,8) =b%var(:,:,:,8) /c%var(:,:,:,1)-c%var(:,:,:,3)*c%var(:,:,:,3)
          b%var(:,:,:,9) =b%var(:,:,:,9) /c%var(:,:,:,1)-c%var(:,:,:,3)*c%var(:,:,:,4)
          b%var(:,:,:,10)=b%var(:,:,:,10)/c%var(:,:,:,1)-c%var(:,:,:,4)*c%var(:,:,:,4)
          b%var(:,:,:,11)=b%var(:,:,:,11)-c%var(:,:,:,6)*c%var(:,:,:,6)

          b%varname(1) = 'tt'
          b%varname(2) = 'tu'
          b%varname(3) = 'tv'
          b%varname(4) = 'tw'
          b%varname(5) = 'uu'
          b%varname(6) = 'uv'
          b%varname(7) = 'uw'
          b%varname(8) = 'vv'
          b%varname(9) = 'vw'
          b%varname(10)= 'ww'
          b%varname(11)= 'pp'

        enddo

        call h5_write_blocks(blocks=block_2nd,fname='Results/2nd.fav.zm.h5')

        call plot_block_xy(blocks=block_2nd,filename='Results/tec2nd')

        deallocate(varxy,var)

      endif
      
      if(mode=='3rdsta' .or. mode=='all') then

        if(allocated(block_mean)) deallocate(block_mean)
        allocate(block_mean, source=pblocks)
        if(allocated(names)) deallocate(names)
        allocate(names(4))
        names(1)='ro'
        names(2)='u'
        names(3)='v'
        names(4)='w'

        call h5_read_blocks(blocks=block_mean,fname='Results/mean.fav.zm.h5',namelist=names)

        allocate( var(0:im,0:jm,0:km))
        allocate( varxy(0:im,0:jm,1:10) )        

        var=read_stats(varname='u111',filename='outdat/3rdsta.h5')
        varxy(:,:,1)=average_k(var)
        var=read_stats(varname='u112',filename='outdat/3rdsta.h5')
        varxy(:,:,2)=average_k(var)
        var=read_stats(varname='u113',filename='outdat/3rdsta.h5')
        varxy(:,:,3)=average_k(var)
        var=read_stats(varname='u122',filename='outdat/3rdsta.h5')
        varxy(:,:,4)=average_k(var)
        var=read_stats(varname='u123',filename='outdat/3rdsta.h5')
        varxy(:,:,5)=average_k(var)
        var=read_stats(varname='u133',filename='outdat/3rdsta.h5')
        varxy(:,:,6)=average_k(var)
        var=read_stats(varname='u222',filename='outdat/3rdsta.h5')
        varxy(:,:,7)=average_k(var)
        var=read_stats(varname='u223',filename='outdat/3rdsta.h5')
        varxy(:,:,8)=average_k(var)
        var=read_stats(varname='u233',filename='outdat/3rdsta.h5')
        varxy(:,:,9)=average_k(var)
        var=read_stats(varname='u333',filename='outdat/3rdsta.h5')
        varxy(:,:,10)=average_k(var)

        allocate(block_3rd, source=pblocks)

        do i=1,nblocks
          b=>block_3rd(i)
          c=>block_mean(i)
          b%nvar=10
          call b%init_data()
          b%var(0:b%im,0:b%jm,0,1:10)=varxy(b%ilo:b%ihi,b%jlo:b%jhi,1:10)

          b%var(:,:,:,1) =b%var(:,:,:,1) /c%var(:,:,:,1)-c%var(:,:,:,2)*c%var(:,:,:,2)*c%var(:,:,:,2)
          b%var(:,:,:,2) =b%var(:,:,:,2) /c%var(:,:,:,1)-c%var(:,:,:,2)*c%var(:,:,:,2)*c%var(:,:,:,3)
          b%var(:,:,:,3) =b%var(:,:,:,3) /c%var(:,:,:,1)-c%var(:,:,:,2)*c%var(:,:,:,2)*c%var(:,:,:,4)
          b%var(:,:,:,4) =b%var(:,:,:,4) /c%var(:,:,:,1)-c%var(:,:,:,2)*c%var(:,:,:,3)*c%var(:,:,:,3)
          b%var(:,:,:,5) =b%var(:,:,:,5) /c%var(:,:,:,1)-c%var(:,:,:,2)*c%var(:,:,:,3)*c%var(:,:,:,4)
          b%var(:,:,:,6) =b%var(:,:,:,6) /c%var(:,:,:,1)-c%var(:,:,:,2)*c%var(:,:,:,4)*c%var(:,:,:,4)
          b%var(:,:,:,7) =b%var(:,:,:,7) /c%var(:,:,:,1)-c%var(:,:,:,3)*c%var(:,:,:,3)*c%var(:,:,:,3)
          b%var(:,:,:,8) =b%var(:,:,:,8) /c%var(:,:,:,1)-c%var(:,:,:,3)*c%var(:,:,:,3)*c%var(:,:,:,4)
          b%var(:,:,:,9) =b%var(:,:,:,9) /c%var(:,:,:,1)-c%var(:,:,:,3)*c%var(:,:,:,4)*c%var(:,:,:,4)
          b%var(:,:,:,10)=b%var(:,:,:,10)/c%var(:,:,:,1)-c%var(:,:,:,4)*c%var(:,:,:,4)*c%var(:,:,:,4)

          b%varname(1) = 'uuu'
          b%varname(2) = 'uuv'
          b%varname(3) = 'uuw'
          b%varname(4) = 'uvv'
          b%varname(5) = 'uvw'
          b%varname(6) = 'uww'
          b%varname(7) = 'vvv'
          b%varname(8) = 'vvw'
          b%varname(9) = 'vww'
          b%varname(10)= 'www'
        enddo

        call h5_write_blocks(blocks=block_3rd,fname='Results/3rd.fav.zm.h5')

        call plot_block_xy(blocks=block_3rd,filename='Results/tec3rd')

        deallocate(varxy,var)

      endif
      
      if(mode=='budget' .or. mode=='all') then
      endif

    end subroutine stats_read_process

    subroutine plot_block_xy(blocks,filename)

      use pastr_commvar, only: im,jm,km
      use pastr_io, only : read_grid
      use pastr_multiblock_type, only: block_define
      use pastr_tecio
      use pastr_commtype, only : tblock

      type(tblock),intent(in), target :: blocks(:)
      character(len=*),intent(in) :: filename

      real(wp),allocatable :: x(:,:),y(:,:),z(:,:)
      integer :: i
      type(tblock),pointer :: b

      allocate(x(0:im,0:jm),y(0:im,0:jm),z(0:im,0:jm))

      call read_grid(x=x,y=y,z=z,kslice=0)

      do i=1,size(blocks)
        b=>blocks(i)

        b%x(0:b%im,0:b%jm,0)=x(b%ilo:b%ihi,b%jlo:b%jhi)
        b%y(0:b%im,0:b%jm,0)=y(b%ilo:b%ihi,b%jlo:b%jhi)

      enddo

      call tecbin(filename=filename,block=blocks)

    end subroutine plot_block_xy

    function average_k(var) result(varm)

      use pastr_commvar, only: im,jm,km

      real(wp),intent(in) :: var(0:im,0:jm,0:km)
      real(wp) :: varm(0:im,0:jm)

      integer :: k

      varm=0._wp
      do k=1,km
        varm(:,:)=varm(:,:)+var(:,:,k)
      enddo
      varm=varm/dble(km)

    end function average_k

end module pastr_data_process