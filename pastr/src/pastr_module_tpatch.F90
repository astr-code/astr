module module_tpatch

    use iso_fortran_env, only: wp => real64

    implicit none

    type :: tpatch

      character(len=2) :: dir
      character(len=3) :: bct
      logical :: tosplit=.false.
      integer :: split_ijk(3)
      integer :: id,blk_id
      integer :: nbr_block,nbr_patch
      integer :: ilo,ihi,jlo,jhi,klo,khi
      integer :: is,ie,js,je,ks,ke

      logical :: data_not_swaped=.true.

      character(len=1),allocatable,dimension(:,:,:) :: cswap
      real(wp),allocatable,dimension(:,:,:,:) :: rswap

      contains

      procedure :: alloc_swap  => alloc_patch_swap

    end type tpatch

contains

    subroutine alloc_patch_swap(apatch,nhalo,dtype,nswp)
      
      class(tpatch),target :: apatch
      integer,intent(in) :: nhalo
      character(len=*),intent(in) :: dtype
      integer,intent(in) :: nswp

      if(apatch%dir(1:1)=='i') then
        if(dtype=='nds') then
          if(allocated(apatch%cswap)) deallocate(apatch%cswap)
          allocate(apatch%cswap(0:nhalo,  &
                                apatch%js:apatch%je,  &
                                apatch%ks:apatch%ke))
        elseif(dtype=='xyz' .or. dtype=='rswap') then
          if(allocated(apatch%rswap)) deallocate(apatch%rswap)
          allocate(apatch%rswap(0:nhalo,  &
                                apatch%js:apatch%je,  &
                                apatch%ks:apatch%ke,nswp))
        endif
      elseif(apatch%dir(1:1)=='j') then
        if(dtype=='nds') then
          if(allocated(apatch%cswap)) deallocate(apatch%cswap)
          allocate(apatch%cswap(apatch%is:apatch%ie,  &
                                0:nhalo,  &
                                apatch%ks:apatch%ke))
        elseif(dtype=='xyz' .or. dtype=='rswap') then
          if(allocated(apatch%rswap)) deallocate(apatch%rswap)
          allocate(apatch%rswap(apatch%is:apatch%ie,  &
                                0:nhalo,  &
                                apatch%ks:apatch%ke,nswp))
        endif
      elseif(apatch%dir(1:1)=='k') then
        if(dtype=='nds') then
          if(allocated(apatch%cswap)) deallocate(apatch%cswap)
          allocate(apatch%cswap(apatch%is:apatch%ie,  &
                                apatch%js:apatch%je,  &
                                0:nhalo))
        elseif(dtype=='xyz' .or. dtype=='rswap') then
          if(allocated(apatch%rswap)) deallocate(apatch%rswap)
          allocate(apatch%rswap(apatch%is:apatch%ie,  &
                                apatch%js:apatch%je,  &
                                0:nhalo,nswp))
        endif
      endif

    end subroutine alloc_patch_swap

end module module_tpatch