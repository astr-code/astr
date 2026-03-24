module module_tdata

    use iso_fortran_env, only: wp => real64

    implicit none

    type :: tdata
      character(len=6) :: name
      real(wp), allocatable, dimension(:,:,:) :: core
      real(wp), allocatable, dimension(:,:,:) :: spi0,spim,  &
                                                 spj0,spjm,  &
                                                 spk0,spkm
      contains

      procedure :: alloc_tdata  => alloc_tdata_data
    end type tdata

contains

    subroutine alloc_tdata_data(adata,im,jm,km,nhalo)

      class(tdata),target :: adata
      integer,intent(in) :: im,jm,km,nhalo

      allocate(adata%core(0:im,0:jm,0:km))
      allocate(adata%spi0(-nhalo:-1,0:jm,0:km))
      allocate(adata%spim(1:nhalo,0:jm,0:km))
      allocate(adata%spj0(0:im,-nhalo:-1,0:km))
      allocate(adata%spjm(0:im,1:nhalo,0:km))
      allocate(adata%spk0(0:im,0:jm,-nhalo:-1))
      allocate(adata%spkm(0:im,0:jm,1:nhalo))

    end subroutine alloc_tdata_data

end module module_tdata