module pastr_udf

    use iso_fortran_env, only: wp => real64
    implicit none

contains

    subroutine udf_func1

      print*,' ** udf_func1 is executed ** '

    end subroutine udf_func1

end module pastr_udf