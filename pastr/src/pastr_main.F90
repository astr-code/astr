module pastr_main_mod

    use iso_fortran_env,  only: wp => real64

    ! use pastr_config_mod, only: pastr_settings, init_default_settings
    ! use pastr_filter_mod, only: run_post_processing
    ! use astr_input_mod,   only: astr_input_type, read_astr_input

    implicit none

    private

    public :: pastr_main

contains

    !====================================================================
    !> Main driver for PASTR post-processing tool
    !!  - Reads ASTR input file (required)
    !!  - Reads PASTR config file (optional)
    !!  - Runs post-processing routines
    !====================================================================
    subroutine pastr_main()

        use pastr_io,     only: parse_command_line
        use pastr_process,only: run_process_entry

        character(len=16) :: pastr_command

        character(len=:), allocatable :: pastr_cfg_file

        logical :: cfg_found

        !----------------------------------------------------------
        ! 1. Command line arguments
        !----------------------------------------------------------
        call parse_command_line(pastr_command)

        print*,' **             command: ',pastr_command

        !----------------------------------------------------------
        ! 2. Run the post-processing
        !----------------------------------------------------------
        call run_process_entry(pastr_command)

        !----------------------------------------------------------
        ! 3. Finish
        !----------------------------------------------------------
        write(*,*)
        write(*,*) " ******************************************* "
        write(*,*) " **      pastr completed successfully     ** "
        write(*,*) " ******************************************* "
        write(*,*)

    end subroutine pastr_main

end module pastr_main_mod

program pastr

    use pastr_main_mod

    implicit none

    ! Call the main driver subroutine
    call pastr_main()

end program pastr
