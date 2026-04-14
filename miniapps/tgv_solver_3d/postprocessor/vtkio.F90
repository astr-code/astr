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

  subroutine write_vtk_ascii(filename, rho, vel, prs, tmp)
      implicit none
      character(len=*), intent(in) :: filename
      real(8), intent(in) :: rho(-hm:im + hm, -hm:jm + hm, -hm:km + hm)
      real(8), intent(in) :: vel(-hm:im + hm, -hm:jm + hm, -hm:km + hm, 1:3)
      real(8), intent(in) :: prs(-hm:im + hm, -hm:jm + hm, -hm:km + hm)
      real(8), intent(in) :: tmp(-hm:im + hm, -hm:jm + hm, -hm:km + hm)

      integer :: i, j, k, m
      integer :: nx, ny, nz
      integer :: unit
      character(len=200) :: line

      nx = im + 1
      ny = jm + 1
      nz = km + 1

      open (newunit=unit, file=filename, form='formatted', status='replace')

      write (unit, '(A)') '# vtk DataFile Version 3.0'
      write (unit, '(A)') 'CFD ASCII output'
      write (unit, '(A)') 'ASCII'
      write (unit, '(A)') 'DATASET STRUCTURED_POINTS'

      ! ===== DIMENSIONS =====
      write (line, '(A,3(I0,1X))') 'DIMENSIONS ', nx, ny, nz
      write (unit, '(A)') trim(line)

      ! ===== ORIGIN =====
      write (unit, '(A)') 'ORIGIN 0.0 0.0 0.0'
      ! ===== SPACING =====
      write (line, '(A,3(F16.9,1X))') 'SPACING ', dx, dy, dz
      write (unit, '(A)') trim(line)

      ! ===== POINT_DATA =====
      write (line, '(A,I0)') 'POINT_DATA ', nx*ny*nz
      write (unit, '(A)') trim(line)

      ! =========================
      ! density
      ! =========================
      write (unit, '(A)') 'SCALARS density double'
      write (unit, '(A)') 'LOOKUP_TABLE default'

      do k = 0, km
         do j = 0, jm
            write (unit, '(*(ES25.15E3,1X))') (rho(i, j, k), i=0, im)
         end do
      end do

      ! =========================
      ! pressure
      ! =========================
      write (unit, '(A)') 'SCALARS pressure double'
      write (unit, '(A)') 'LOOKUP_TABLE default'

      do k = 0, km
         do j = 0, jm
            write (unit, '(*(ES25.15E3,1X))') (prs(i, j, k), i=0, im)
         end do
      end do

      ! =========================
      ! temperature
      ! =========================
      write (unit, '(A)') 'SCALARS temperature double'
      write (unit, '(A)') 'LOOKUP_TABLE default'

      do k = 0, km
         do j = 0, jm
            write (unit, '(*(ES25.15E3,1X))') (tmp(i, j, k), i=0, im)
         end do
      end do

      ! =========================
      ! velocity
      ! =========================
      write (unit, '(A)') 'VECTORS velocity double'

      do k = 0, km
         do j = 0, jm
            write (unit, '(*(ES25.15E3,1X))') ((vel(i, j, k, m), m=1, 3), i=0, im)
         end do
      end do

      close (unit)

  end subroutine write_vtk_ascii

  subroutine read_vtk_binary(filename, rho, vel, prs, tmp)
      implicit none
      character(len=*), intent(in) :: filename
      real(8), intent(out) :: rho(0:im, 0:jm, 0:km)
      real(8), intent(out) :: vel(0:im, 0:jm, 0:km, 1:3)
      real(8), intent(out) :: prs(0:im, 0:jm, 0:km)
      real(8), intent(out) :: tmp(0:im, 0:jm, 0:km)

      integer :: i, j, k, m
      integer :: nx, ny, nz, nx_read, ny_read, nz_read
      integer :: unit, ios, point_data_size
      character(len=200) :: line
      character(len=50) :: keyword
      character :: c
      real(8) :: dx_read, dy_read, dz_read, origin_x, origin_y, origin_z

      nx = im + 1
      ny = jm + 1
      nz = km + 1

      open (newunit=unit, file=filename, access='stream', form='unformatted', &
            status='old', convert='big_endian', iostat=ios)

      if (ios /= 0) then
         print *, 'Error opening file: ', filename
         return
      end if

      ! Read header lines
      call read_line(unit, line)  ! '# vtk DataFile Version 3.0'
      call read_line(unit, line)  ! 'CFD binary output'
      call read_line(unit, line)  ! 'BINARY'
      call read_line(unit, line)  ! 'DATASET STRUCTURED_POINTS'

      ! Read DIMENSIONS
      call read_line(unit, line)
      read (line, *) keyword, nx_read, ny_read, nz_read
      if (nx_read /= nx .or. ny_read /= ny .or. nz_read /= nz) then
         print *, 'Warning: Dimensions mismatch!'
         print *, 'Expected:', nx, ny, nz
         print *, 'Read:', nx_read, ny_read, nz_read
      end if

      ! Read ORIGIN
      call read_line(unit, line)
      read (line, *) keyword, origin_x, origin_y, origin_z

      ! Read SPACING
      call read_line(unit, line)
      read (line, *) keyword, dx_read, dy_read, dz_read

      ! Read POINT_DATA
      call read_line(unit, line)
      read (line, *) keyword, point_data_size

      ! =========================
      ! Read density
      ! =========================
      call read_line(unit, line)  ! 'SCALARS density double'
      call read_line(unit, line)  ! 'LOOKUP_TABLE default'

      do k = 0, km
         do j = 0, jm
            read (unit) (rho(i, j, k), i=0, im)
         end do
      end do
      read (unit) c  ! consume the trailing newline

      ! =========================
      ! Read pressure
      ! =========================
      call read_line(unit, line)  ! 'SCALARS pressure double'
      call read_line(unit, line)  ! 'LOOKUP_TABLE default'

      do k = 0, km
         do j = 0, jm
            read (unit) (prs(i, j, k), i=0, im)
         end do
      end do
      read (unit) c  ! consume the trailing newline

      ! =========================
      ! Read temperature
      ! =========================
      call read_line(unit, line)  ! 'SCALARS temperature double'
      call read_line(unit, line)  ! 'LOOKUP_TABLE default'

      do k = 0, km
         do j = 0, jm
            read (unit) (tmp(i, j, k), i=0, im)
         end do
      end do
      read (unit) c  ! consume the trailing newline

      ! =========================
      ! Read velocity
      ! =========================
      call read_line(unit, line)  ! 'VECTORS velocity double'

      do k = 0, km
         do j = 0, jm
            read (unit) ((vel(i, j, k, m), m=1, 3), i=0, im)
         end do
      end do

      close (unit)

  end subroutine read_vtk_binary

! Helper subroutine to read a line from binary stream file
  subroutine read_line(unit, line)
      implicit none
      integer, intent(in) :: unit
      character(len=*), intent(out) :: line
      character :: c
      integer :: i

      line = ' '
      i = 1
      do
         read (unit) c
         if (c == char(10)) exit  ! newline character
         if (i <= len(line)) then
            line(i:i) = c
            i = i + 1
         end if
      end do
  end subroutine read_line

end module vtkio