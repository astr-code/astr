module pastr_io

    use iso_fortran_env, only: wp => real64

    implicit none

    private

    public :: parse_command_line,read_grid
    public :: global_n_arguments

    integer :: global_n_arguments=0

    Interface read_grid
      module procedure read_grid_2d
    end Interface read_grid
    !
contains

    !------------------------------------------------------------------
    !> Automatically parse command-line arguments into strings and integers
    !!  - Arguments that can be read as integers go to `numbers`
    !!  - All other arguments go to `strings`
    !------------------------------------------------------------------
    subroutine parse_command_line(string, inumber,rnumber)
        
        use pastr_utility

        character(len=*), intent(inout), optional :: string
        integer, intent(inout), optional :: inumber
        real(wp), intent(inout), optional :: rnumber

        integer :: n_args, i, n_count, s_count
        character(len=255) :: arg
        integer :: num
        logical :: proper_assign

        proper_assign=.false.

        global_n_arguments=global_n_arguments+1

        call get_cmd_argument(global_n_arguments,arg)
        
        if(present(string)) string=''
        if(present(inumber)) inumber=-1
        if(present(rnumber)) rnumber=0.d0

        if(isnum(arg)==0) then
          ! non-numeric
          if(present(string)) then
            string=trim(arg)
            proper_assign=.true.
          endif

        elseif(isnum(arg)==1) then
          ! integer
          if(present(inumber)) then
            read(arg,*)inumber
            proper_assign=.true.
          endif

        elseif(isnum(arg)>1) then
          ! real number
          if(present(rnumber)) then
            read(arg,*)rnumber
            proper_assign=.true.
          endif
        else
          stop ' !!  error1 @ parse_command_line !!'
        endif

        if(.not.  proper_assign) global_n_arguments=global_n_arguments-1
        

    end subroutine parse_command_line

    subroutine read_grid_2d(x,y,z,islice,jslice,kslice)
      
      use pastr_commvar, only: gridfile,im,jm,km
      use pastr_h5io

      real(wp),intent(inout),allocatable,optional :: x(:,:),y(:,:),z(:,:)
      integer,intent(in),optional :: islice,jslice,kslice

      if(present(islice)) then
        if(present(x)) call h5_read2dfrom3d(x,im,jm,km,'x',trim(gridfile),islice=islice)
        if(present(y)) call h5_read2dfrom3d(y,im,jm,km,'y',trim(gridfile),islice=islice)
        if(present(z)) call h5_read2dfrom3d(z,im,jm,km,'z',trim(gridfile),islice=islice)
      elseif(present(jslice)) then
        if(present(x)) call h5_read2dfrom3d(x,im,jm,km,'x',trim(gridfile),jslice=jslice)
        if(present(y)) call h5_read2dfrom3d(y,im,jm,km,'y',trim(gridfile),jslice=jslice)
        if(present(z)) call h5_read2dfrom3d(z,im,jm,km,'z',trim(gridfile),jslice=jslice)
      elseif(present(kslice)) then
        if(present(x)) call h5_read2dfrom3d(x,im,jm,km,'x',trim(gridfile),kslice=kslice)
        if(present(y)) call h5_read2dfrom3d(y,im,jm,km,'y',trim(gridfile),kslice=kslice)
        if(present(z)) call h5_read2dfrom3d(z,im,jm,km,'z',trim(gridfile),kslice=kslice)
      endif

    end subroutine read_grid_2d

end module pastr_io
