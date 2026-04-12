module vtkio
  use commvar, only: im,jm,km,hm
  use commarray, only:dx,dy,dz
  implicit none

contains

  subroutine write_vtk_binary(filename, rho, vel, prs, tmp)
    implicit none
    character(len=*), intent(in) :: filename
    real(8), intent(in) :: rho(-hm:im+hm,-hm:jm+hm,-hm:km+hm)
    real(8), intent(in) :: vel(-hm:im+hm,-hm:jm+hm,-hm:km+hm,1:3)
    real(8), intent(in) :: prs(-hm:im+hm,-hm:jm+hm,-hm:km+hm)
    real(8), intent(in) :: tmp(-hm:im+hm,-hm:jm+hm,-hm:km+hm)

    integer :: i,j,k
    integer :: nx,ny,nz
    integer :: unit
    character(len=200) :: line

    nx = im + 1
    ny = jm + 1
    nz = km + 1

    open(newunit=unit, file=filename, access='stream', form='unformatted', &
         status='replace', convert='big_endian')

    write(unit) '# vtk DataFile Version 3.0'//char(10)
    write(unit) 'CFD binary output'//char(10)
    write(unit) 'BINARY'//char(10)
    write(unit) 'DATASET STRUCTURED_POINTS'//char(10)

    ! ===== DIMENSIONS =====
    write(line,'(A,3(I0,1X))') 'DIMENSIONS ', nx, ny, nz
    write(unit) trim(line)//char(10)

    ! ===== ORIGIN =====
    write(unit) 'ORIGIN 0.0 0.0 0.0'//char(10)

    ! ===== SPACING =====
    write(line,'(A,3(F15.10,1X))') 'SPACING ', dx, dy, dz
    write(unit) trim(line)//char(10)

    ! ===== POINT_DATA =====
    write(line,'(A,I0)') 'POINT_DATA ', nx*ny*nz
    write(unit) trim(line)//char(10)

    ! =========================
    ! density
    ! =========================
    write(unit) 'SCALARS density double'//char(10)
    write(unit) 'LOOKUP_TABLE default'//char(10)

    do k=0,km
      do j=0,jm
        write(unit) (rho(i,j,k), i=0,im)
      end do
    end do
    write(unit) char(10)

    ! =========================
    ! pressure
    ! =========================
    write(unit) 'SCALARS pressure double'//char(10)
    write(unit) 'LOOKUP_TABLE default'//char(10)

    do k=0,km
      do j=0,jm
        write(unit) (prs(i,j,k), i=0,im)
      end do
    end do
    write(unit) char(10)

    ! =========================
    ! temperature
    ! =========================
    write(unit) 'SCALARS temperature double'//char(10)
    write(unit) 'LOOKUP_TABLE default'//char(10)

    do k=0,km
      do j=0,jm
        write(unit) (tmp(i,j,k), i=0,im)
      end do
    end do
    write(unit) char(10)

    ! =========================
    ! velocity
    ! =========================
    write(unit) 'VECTORS velocity double'//char(10)

    do k=0,km
      do j=0,jm
        do i=0,im
          write(unit) vel(i,j,k,1), vel(i,j,k,2), vel(i,j,k,3)
        end do
      end do
    end do

    close(unit)

  end subroutine write_vtk_binary

end module vtkio