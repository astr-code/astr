module pastr_commtype

    use iso_fortran_env, only: wp => real64
    use pastr_constdef
    use pastr_fdm
    use module_tdata
    use module_tpatch

    implicit none

    type :: montype
      integer :: npoints,nvariables
      character(len=16),allocatable :: varname(:)
      integer,allocatable :: nstep(:)
      real(wp),allocatable :: time(:),data(:,:)
      contains
      procedure :: init => alloc_monitor
    end type montype


    type :: tblock
      character(len=6) :: name
      integer :: id
      integer :: im,jm,km,nvar,nres,nswp
      integer :: ilo,ihi,jlo,jhi,klo,khi
      character(len=16),allocatable :: varname(:),resname(:)
      character(len=1),allocatable,dimension(:,:,:) :: nds
      real(wp), allocatable, dimension(:,:,:) :: x,y,z
      type(tdata) :: xt,yt,zt

      real(wp), allocatable, dimension(:,:,:,:) :: ddi,ddj,ddk
      real(wp), allocatable, dimension(:,:,:,:) :: var
      real(wp), allocatable, dimension(:,:,:,:) :: res
      real(wp), allocatable, dimension(:,:,:,:) :: rswap
      integer :: npatches=0
      type(tpatch),allocatable :: patches(:)

      type(fdm_type),allocatable :: fdi(:),fdj(:),fdk(:)
      integer,allocatable :: fdi_map(:,:),fdj_map(:,:),fdk_map(:,:)

      contains

      procedure :: init_grid  => alloc_block_grid
      procedure :: init_var_data  => alloc_block_var
      procedure :: init_res_data  => alloc_block_res
      procedure :: init_rswap  => alloc_block_rswap
      procedure :: data_map   => block_data_map
      procedure :: add_patch  => add_block_patch
      procedure :: patch_info => print_patch_info
    end type tblock


contains


    subroutine alloc_monitor(amonitor)

      class(montype),target :: amonitor

      allocate(amonitor%varname(amonitor%nvariables))
      allocate(amonitor%time(amonitor%npoints))
      allocate(amonitor%nstep(amonitor%npoints))
      allocate(amonitor%data(amonitor%nvariables,amonitor%npoints))

    end subroutine alloc_monitor

    subroutine print_patch_info(ablock)

      class(tblock),target :: ablock

      integer :: i

      print*,' ** number of patches:',ablock%npatches,' in block ',ablock%id
      do i=1,ablock%npatches
        print*,' ** patch',i,' from block',ablock%id,' -dir: ',ablock%patches(i)%dir
        print*,' ** dimension: ',ablock%patches(i)%is,ablock%patches(i)%ie, &
                                 ablock%patches(i)%js,ablock%patches(i)%je, &
                                 ablock%patches(i)%ks,ablock%patches(i)%ke
        print*,' ** connected to patch ',ablock%patches(i)%nbr_patch,' from block ',ablock%patches(i)%nbr_block
      enddo
      print*,' ** --------------------------------------------------------------------------- **'

    end subroutine print_patch_info

    subroutine alloc_block_grid(ablock,nhalo)

      class(tblock),target :: ablock
      integer,intent(in) :: nhalo

      if(.not.allocated(ablock%x)) allocate(ablock%x(-nhalo:ablock%im+nhalo, &
                                                     -nhalo:ablock%jm+nhalo, &
                                                     -nhalo:ablock%km+nhalo))
      if(.not.allocated(ablock%y)) allocate(ablock%y(-nhalo:ablock%im+nhalo, &
                                                     -nhalo:ablock%jm+nhalo, &
                                                     -nhalo:ablock%km+nhalo))
      if(.not.allocated(ablock%z)) allocate(ablock%z(-nhalo:ablock%im+nhalo, &
                                                     -nhalo:ablock%jm+nhalo, &
                                                     -nhalo:ablock%km+nhalo))
      
      call ablock%xt%alloc_tdata(ablock%im,ablock%jm,ablock%km,nhalo)
      call ablock%yt%alloc_tdata(ablock%im,ablock%jm,ablock%km,nhalo)
      call ablock%zt%alloc_tdata(ablock%im,ablock%jm,ablock%km,nhalo)

    end subroutine alloc_block_grid

    subroutine alloc_block_var(ablock)

      class(tblock),target :: ablock

      allocate(ablock%varname(ablock%nvar))
      allocate(ablock%var(0:ablock%im,0:ablock%jm,0:ablock%km,1:ablock%nvar))

    end subroutine alloc_block_var

    subroutine alloc_block_res(ablock)

      class(tblock),target :: ablock

      allocate(ablock%resname(ablock%nres))
      allocate(ablock%res(0:ablock%im,0:ablock%jm,0:ablock%km,1:ablock%nres))

    end subroutine alloc_block_res

    subroutine alloc_block_rswap(ablock,n,nhalo)

      class(tblock),target :: ablock
      integer,intent(in) :: n,nhalo

      ablock%nswp=n
      allocate(ablock%rswap(-nhalo:ablock%im+nhalo, &
                            -nhalo:ablock%jm+nhalo, &
                            -nhalo:ablock%km+nhalo,1:ablock%nswp))

    end subroutine alloc_block_rswap

    subroutine block_data_map(ablock,datain)

      class(tblock),target :: ablock

      real(wp),intent(in) :: datain(0:,0:,0:,1:)

      ablock%var(0:ablock%im, &
                 0:ablock%jm, &
                 0:ablock%km, &
                 1:ablock%nvar)=datain(ablock%ilo:ablock%ihi, &
                                       ablock%jlo:ablock%jhi, &
                                       ablock%klo:ablock%khi, &
                                       1:ablock%nvar   )

    end subroutine block_data_map

    subroutine add_block_patch(ablock,patch_in)

      class(tblock),target :: ablock
      type(tpatch),intent(in) :: patch_in

      type(tpatch), allocatable :: tmp(:)

      if(.not. allocated(ablock%patches)) then
        ablock%npatches=1
        allocate(ablock%patches(ablock%npatches))
        ablock%patches(1)=patch_in
        ablock%patches(1)%id=ablock%npatches
        ablock%patches(1)%blk_id=ablock%id
      else
        allocate(tmp(ablock%npatches))
        tmp(:) = ablock%patches(:)
        deallocate(ablock%patches)

        ablock%npatches=ablock%npatches+1
        allocate(ablock%patches(ablock%npatches))

        ablock%patches(1:ablock%npatches-1)=tmp(:)
        ablock%patches(ablock%npatches)=patch_in
        ablock%patches(ablock%npatches)%id=ablock%npatches
        ablock%patches(ablock%npatches)%blk_id=ablock%id

        deallocate(tmp)
      endif

    end subroutine add_block_patch

    subroutine pack_data_slap(adata,apatch,nrswap,nhalo,dtype,etype)

      type(tdata),intent(inout) :: adata
      type(tpatch),intent(in) :: apatch
      character(len=*),intent(in) :: dtype,etype
      integer,intent(in) :: nrswap,nhalo

      integer :: i,j,k

      if(dtype=='xyz') then
        
        select case(apatch%dir)
        case('i-')
          if(etype=='pack') then
            do i=1,nhalo
              adata%spi0(-i,apatch%js:apatch%je,apatch%ks:apatch%ke)= &
              adata%core(i,apatch%js:apatch%je,apatch%ks:apatch%ke)- &
              adata%core(0,apatch%js:apatch%je,apatch%ks:apatch%ke)
            end do
          elseif(etype=='unpack') then
            do i=1,nhalo
                adata%spi0(-i,apatch%js:apatch%je,apatch%ks:apatch%ke)= &
              apatch%rswap(nhalo-i,apatch%js:apatch%je,apatch%ks:apatch%ke,nrswap)+ &
                adata%core( 0,apatch%js:apatch%je,apatch%ks:apatch%ke)
            end do
          endif
        case('i+')
          if(etype=='pack') then
            do i=1,nhalo
              adata%spim(          i,apatch%js:apatch%je,apatch%ks:apatch%ke)= &
              adata%core(apatch%is-i,apatch%js:apatch%je,apatch%ks:apatch%ke)- &
              adata%core(  apatch%is,apatch%js:apatch%je,apatch%ks:apatch%ke)
            end do
          elseif(etype=='unpack') then
            do i=1,nhalo
                adata%spim(        i,apatch%js:apatch%je,apatch%ks:apatch%ke)= &
              apatch%rswap(        i,apatch%js:apatch%je,apatch%ks:apatch%ke,nrswap)+ &
                adata%core(apatch%is,apatch%js:apatch%je,apatch%ks:apatch%ke)
            end do
          endif
        case('j-')
          if(etype=='pack') then
            do j=1,nhalo
              adata%spj0(apatch%is:apatch%ie,-j,apatch%ks:apatch%ke)= &
              adata%core(apatch%is:apatch%ie, j,apatch%ks:apatch%ke)- &
              adata%core(apatch%is:apatch%ie, 0,apatch%ks:apatch%ke)
            end do
          elseif(etype=='unpack') then
            do j=1,nhalo
                adata%spj0(apatch%is:apatch%ie,     -j,apatch%ks:apatch%ke)= &
              apatch%rswap(apatch%is:apatch%ie,nhalo-j,apatch%ks:apatch%ke,nrswap) + &
                adata%core(apatch%is:apatch%ie,      0,apatch%ks:apatch%ke)
            end do
          endif
        case('j+')
          if(etype=='pack') then
            do j=1,nhalo
              adata%spjm(apatch%is:apatch%ie,          j,apatch%ks:apatch%ke)= &
              adata%core(apatch%is:apatch%ie,apatch%js-j,apatch%ks:apatch%ke)- &
              adata%core(apatch%is:apatch%ie,  apatch%js,apatch%ks:apatch%ke)
            end do
          elseif(etype=='unpack') then
            do j=1,nhalo
                adata%spjm(apatch%is:apatch%ie,        j,apatch%ks:apatch%ke)= &
              apatch%rswap(apatch%is:apatch%ie,        j,apatch%ks:apatch%ke,nrswap)+ &
               adata%core(apatch%is:apatch%ie, apatch%js,apatch%ks:apatch%ke)
            end do
          endif
        case('k-')
          if(etype=='pack') then
            do k=1,nhalo
              adata%spk0(apatch%is:apatch%ie,apatch%js:apatch%je,-k)= &
              adata%core(apatch%is:apatch%ie,apatch%js:apatch%je, k)- &
              adata%core(apatch%is:apatch%ie,apatch%js:apatch%je, 0)
            end do
          elseif(etype=='unpack') then
            do k=1,nhalo
                adata%spk0(apatch%is:apatch%ie,apatch%js:apatch%je,    -k)= &
              apatch%rswap(apatch%is:apatch%ie,apatch%js:apatch%je,nhalo-k,nrswap)+ &
                adata%core(apatch%is:apatch%ie,apatch%js:apatch%je,0)
            end do
          endif
        case('k+')
          if(etype=='pack') then
            do k=1,nhalo
              adata%spkm(apatch%is:apatch%ie,apatch%js:apatch%je,k)= &
              adata%core(apatch%is:apatch%ie,apatch%js:apatch%je,apatch%ks-k)- &
              adata%core(apatch%is:apatch%ie,apatch%js:apatch%je,apatch%ks)
            end do
          elseif(etype=='unpack') then
            do k=1,nhalo
                adata%spkm(apatch%is:apatch%ie,apatch%js:apatch%je,k)= &
              apatch%rswap(apatch%is:apatch%ie,apatch%js:apatch%je,k,nrswap)+ &
                adata%core(apatch%is:apatch%ie,apatch%js:apatch%je,apatch%ks)
            end do
          endif
        case default
          stop '1 @ pack_data_slap'
        end select

      endif

    end subroutine pack_data_slap

end module pastr_commtype