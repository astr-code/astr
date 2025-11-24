module pastr_input

    use iso_fortran_env, only: wp => real64

    implicit none

    private

    public ::  read_astr_input

contains

    !----------------------------------------------------------------------
    !> Read ASTR input file and fill astr_input_type
    !!  Format assumed:
    !!     key = value
    !!  Example:
    !!     nx = 200
    !!     ny = 150
    !!     dt = 0.0005
    !!  Missing keys use defaults
    !----------------------------------------------------------------------
    subroutine read_astr_input()

        use pastr_commvar
        use pastr_io, only: parse_command_line

        character(len=512) :: input_file,line, key, val
        integer :: ios, iunit
        logical :: eof_flag
        integer :: eq_pos
        character(len=6) :: iden
        logical :: file_exist

        ! Initialize with defaults

        call parse_command_line(string=iden)
        if(iden=='-input') then
          call parse_command_line(string=input_file)
        endif

        inquire(file=trim(input_file), exist=file_exist)
        if (.not. file_exist) then
            write(*,*) "Error: ASTR input file '", trim(input_file), "' does not exist!"
            stop 1
        end if

        ! Open file
        open(newunit=iunit, file=trim(input_file), status='old', action='read', iostat=ios)
        read(iunit,'(///////)')
        read(iunit,*)im,jm,km
        read(iunit,"(/)")
        read(iunit,*)lihomo,ljhomo,lkhomo
        read(iunit,"(/)")
        read(iunit,*)nondimen
        read(iunit,'(///////)')
        if(nondimen) then
          read(iunit,*)ref_t,reynolds,mach
        else
          read(iunit,*)ref_t,ref_vel,ref_len,ref_den
        endif
        read(iunit,'(///////////////////////////)')
        read(iunit,'(A)')gridfile
        close(iunit)
        write(*,*)' >> ',trim(input_file)

        write(*,*)' ==========================informout=========================='
        write(*,*)'                    im                  jm                  km'
        write(*,"(3X,3(I20))")im,jm,km
        write(*,*)' -------------------------------------------------------------'
        write(*,*)'              Reynolds                Mach               ref_t'
        write(*,"(3X,3(F20.5))")Reynolds,mach,ref_t
        write(*,*)' -------------------------------------------------------------'
        write(*,"(3X,A40,A20)")'grid: ',trim(gridfile)
        write(*,*)' ==========================informout=========================='

    end subroutine read_astr_input

end module pastr_input
