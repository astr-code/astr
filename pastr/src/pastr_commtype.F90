module pastr_commtype

    use iso_fortran_env, only: wp => real64
    use pastr_constdef

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
      integer :: im,jm,km,nvar
      integer :: ilo,ihi,jlo,jhi,klo,khi
      character(len=16),allocatable :: varname(:)
      real(wp), allocatable, dimension(:,:,:) :: x,y,z
      real(wp), allocatable, dimension(:,:,:,:) :: var
      integer :: npatches=0
      type(tpatch),allocatable :: patches(:)

      contains

      procedure :: init_data => alloc_block_data
      procedure :: add_patch => add_block_patch
      procedure :: patch_info => print_patch_info
    end type tblock

    type :: tpatch

      character(len=2) :: dir
      integer :: id
      integer :: nbr_block,nbr_patch
      integer :: ilo,ihi,jlo,jhi,klo,khi

    end type tpatch

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
        print*,' ** patch',i,'  --  ',ablock%patches(i)%dir
        print*,ablock%patches(i)%ilo,ablock%patches(i)%ihi, &
               ablock%patches(i)%jlo,ablock%patches(i)%jhi, &
               ablock%patches(i)%klo,ablock%patches(i)%khi
        print*,' ** connected to patch ',ablock%patches(i)%nbr_patch,' from block ',ablock%patches(i)%nbr_block
      enddo
      print*,' ** --------------------------------------------------------------------------- **'

    end subroutine print_patch_info

    subroutine alloc_block_data(ablock)

      class(tblock),target :: ablock


      if(.not.allocated(ablock%x)) allocate(ablock%x(0:ablock%im,0:ablock%jm,0:ablock%km))
      if(.not.allocated(ablock%y)) allocate(ablock%y(0:ablock%im,0:ablock%jm,0:ablock%km))
      if(.not.allocated(ablock%z)) allocate(ablock%z(0:ablock%im,0:ablock%jm,0:ablock%km))
      allocate(ablock%varname(ablock%nvar))
      allocate(ablock%var(0:ablock%im,0:ablock%jm,0:ablock%km,1:ablock%nvar))

    end subroutine alloc_block_data

    subroutine add_block_patch(ablock,patch_in)

      class(tblock),target :: ablock
      type(tpatch),intent(in) :: patch_in

      type(tpatch), allocatable :: tmp(:)

      if(.not. allocated(ablock%patches)) then
        ablock%npatches=1
        allocate(ablock%patches(ablock%npatches))
        ablock%patches(1)=patch_in
        ablock%patches(1)%id=ablock%npatches
      else
        allocate(tmp(ablock%npatches))
        tmp(:) = ablock%patches(:)
        deallocate(ablock%patches)

        ablock%npatches=ablock%npatches+1
        allocate(ablock%patches(ablock%npatches))

        ablock%patches(1:ablock%npatches-1)=tmp(:)
        ablock%patches(ablock%npatches)=patch_in
        ablock%patches(ablock%npatches)%id=ablock%npatches

        deallocate(tmp)
      endif

    end subroutine add_block_patch

end module pastr_commtype