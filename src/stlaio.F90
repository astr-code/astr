!+---------------------------------------------------------------------+
!| This module contains subroutines and functions to read and write    |
!| sla file.                                                           |
!+---------------------------------------------------------------------+
!| src: https://people.sc.fsu.edu/~jburkardt/f_src/stla_io/stla_io.html
!+---------------------------------------------------------------------+
!| CHANGE RECORD                                                       |
!| -------------                                                       |
!| 01-07-2021  | Added by J. Fang                                      |
!+---------------------------------------------------------------------+
module stlaio
!
implicit none
!
contains
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to get the number of stl bodies           !
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 02-07-2021  | Created by J. Fang @ Warrington                     |
  !+-------------------------------------------------------------------+
  integer function num_solid(input_file_name)
    !
    ! argument
    character(len=*),intent(in) :: input_file_name
    !
    ! local data
    integer :: iunit,ios,ierror,state,count
    character(len=255) :: word1,text
    logical :: done
    !
    iunit=get_unit()
    !
    state=0
    num_solid=0
    count=0
    !
    open(unit=iunit,file=input_file_name,status='old',iostat=ios )
    !
    if ( ios /= 0 ) then
      write ( *, '(2A)' ) '  Could not open the file "',               &
                                                  trim (input_file_name)
      ierror = 1
      stop
    end if
    ! !
    do while(ios==0)
      !
      read(iunit,'(a)',iostat=ios) text
      !
      if( ios /= 0 ) then
        exit
      end if
      !
      done = .true.
      !  Read the first word in the line.
      call word_next_read ( text, word1, done )
      if(s_eqi(word1,'SOLID')) then
        if ( state /= 0 ) exit
        state = 1
      elseif(s_eqi(word1,'ENDSOLID')) then
        if ( state /= 1 ) exit
        !
        state = 0
        !
        num_solid = num_solid + 1
      endif
      !
    enddo
    !
    close(iunit)
    !
    return
    !
  end function num_solid
  !+-------------------------------------------------------------------+
  !| The end of the function num_solid.                                |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to get the number of face within 1 solid. !
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 02-07-2021  | Created by J. Fang @ Warrington                     |
  !+-------------------------------------------------------------------+
  integer function num_face_in_solid(input_file_name,nsolid)
    !
    ! argument
    character(len=*),intent(in) :: input_file_name
    integer,intent(in) :: nsolid
    !
    ! local data
    integer :: iunit,ios,ierror,solidstate,num_solid
    character(len=255) word1,word2,text
    logical :: done
    !
    !
    solidstate=0
    num_solid=0
    num_face_in_solid=0
    !
    iunit=get_unit()
    !
    open(unit=iunit,file=input_file_name,status='old',iostat=ios )
    !
    !
    if ( ios /= 0 ) then
      write ( *, '(2A)' ) '  Could not open the file "',               &
                                                  trim (input_file_name)
      ierror = 1
      stop
    end if
    !
    do while(ios==0)
      !
      if( ios /= 0 ) exit
      !
      done = .true.
      !
      !  Read the first word in the line.
      read(iunit,'(a)',iostat=ios) text
      !
      if( ios /= 0 ) then
        exit
      end if
      !
      done = .true.
      !
      call word_next_read ( text, word1, done )
      !
      if(s_eqi(word1,'SOLID')) then
        if ( solidstate /= 0 ) exit
        !
        solidstate = 1
        !
        num_solid = num_solid + 1
        !
      elseif(s_eqi(word1,'ENDSOLID')) then
        if ( solidstate /= 1 ) exit
        !
        solidstate = 0
        !
      elseif(s_eqi(word1,'FACET')) then
        !
        call word_next_read ( text, word2, done )
        !
        if(s_eqi ( word2, 'NORMAL' )) then
          !
          if(solidstate==1 .and. num_solid==nsolid) then
            num_face_in_solid=num_face_in_solid+1
          endif
          !
        endif
        !
      endif
      !
    enddo
    !
    if(num_face_in_solid==0) then
      print*,' !! no face can be found for solid: ',nsolid
    endif
    !
    close(iunit)
    !
    return
    !
  end function num_face_in_solid
  !+-------------------------------------------------------------------+
  !| The end of the function num_face_in_solid.                        |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to read a solid data in STL file.         !
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 02-07-2021  | Created by J. Fang @ Warrington                     |
  !+-------------------------------------------------------------------+
  subroutine stla_read(input_file_name,asolid)
    !
    use commtype,  only : solid,triangle
    !
    ! argument
    character(len=*),intent(in) :: input_file_name
    type(solid),intent(inout) :: asolid(:)
    !
    ! local data
    integer :: iunit,ios,ierror,solidstate,lchar,i
    integer :: num_solid,num_face,num_vertex
    character(len=255) word1,word2,word3,text
    logical :: done
    real(8) :: dval
    !
    !
    solidstate=0
    num_solid=0
    !
    iunit=get_unit()
    !
    open(unit=iunit,file=input_file_name,status='old',iostat=ios )
    !
    !
    if ( ios /= 0 ) then
      write ( *, '(2A)' ) '  Could not open the file "',               &
                                                  trim (input_file_name)
      ierror = 1
      stop
    end if
    !
    do while(ios==0)
      !
      if( ios /= 0 ) exit
      !
      done = .true.
      !
      !  Read the first word in the line.
      read(iunit,'(a)',iostat=ios) text
      !
      if( ios /= 0 ) then
        exit
      end if
      !
      done = .true.
      !
      call word_next_read(text,word1,done)
      !
      if(s_eqi(word1,'SOLID')) then
        if ( solidstate /= 0 ) exit
        !
        solidstate = 1
        !
        num_solid = num_solid + 1
        num_face=0
        !
        call word_next_read(text,word2,done)
        !
        asolid(num_solid)%name=trim(word2)
        !
      elseif(s_eqi(word1,'ENDSOLID')) then
        if ( solidstate /= 1 ) exit
        !
        solidstate = 0
        !
      elseif(s_eqi(word1,'FACET')) then
        !
        call word_next_read(text,word2,done)
        !
        if(s_eqi ( word2, 'NORMAL' )) then
          !
          if(solidstate==1) then
            !
            num_face=num_face+1
            num_vertex=0
            !
            do i=1,3
              !
              call word_next_read(text,word3,done)
              call s_to_r8(word3,dval,ierror,lchar)
              !
              asolid(num_solid)%face(num_face)%normdir(i)=dval
              !
            enddo
            ! print*,num_solid,num_face,'-',asolid(num_solid)%face(num_face)%normdir(:)
            !
          else
            stop ' !! ERROR 1 @ stla_read'
          endif
          !
        endif
        !
      elseif(s_eqi(word1,'VERTEX')) then
        !
        num_vertex=num_vertex+1
        !
        do i=1,3
          !
          call word_next_read(text,word2,done)
          call s_to_r8(word2,dval,ierror,lchar)
          !
          select case(num_vertex)
          case(1)
            asolid(num_solid)%face(num_face)%a(i)=dval
          case(2)
            asolid(num_solid)%face(num_face)%b(i)=dval
          case(3)
            asolid(num_solid)%face(num_face)%c(i)=dval
          case default
            stop ' !! ERROR 2 @ stla_read'
          end select
        enddo
        !
      endif
      !
    enddo
    !
    close(iunit)
    print*,' >> ',input_file_name
    !
    return
    !
  end subroutine stla_read
  !+-------------------------------------------------------------------+
  !| The end of the function stla_read.                                |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to write a STL file.                      !
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 02-07-2021  | Created by J. Fang @ Warrington                     |
  !+-------------------------------------------------------------------+
  subroutine stla_write(out_file_name,asolid)
    !
    use commtype,  only : solid,triangle
    !
    ! argument
    character(len=*),intent(in) :: out_file_name
    type(solid),intent(in) :: asolid(:)
    !
    ! local data
    integer :: iunit,ios,ierror,solidstate,lchar,i
    integer :: num_solid,nsolid,nface
    character(len=255) word1,word2,word3,text
    logical :: done
    real(8) :: dval
    !
    num_solid=size(asolid)
    !
    iunit=get_unit()
    !
    open(unit=iunit,file=out_file_name,iostat=ios)
    do nsolid=1,num_solid
      write(iunit,'(A,I0)')'solid ',nsolid
      do nface=1,asolid(nsolid)%num_face
        write(iunit,'(A,3(1X,E15.7E3))')'  facet normal ',             &
                                   asolid(nsolid)%face(nface)%normdir(:)
        write(iunit,'(A)')'    outer loop'
        write(iunit,'(A,3(1X,E15.7E3))')'      vertex ',              &
                                         asolid(nsolid)%face(nface)%a(:)
        write(iunit,'(A,3(1X,E15.7E3))')'      vertex ',              &
                                         asolid(nsolid)%face(nface)%b(:)
        write(iunit,'(A,3(1X,E15.7E3))')'      vertex ',              &
                                         asolid(nsolid)%face(nface)%c(:)
        write(iunit,'(A)')'    endloop'
        write(iunit,'(A)')'  endfacet'
      enddo
      write(iunit,'(A,I0)')'endsolid ',nsolid
    enddo
    close(iunit)
    print*,' << ',out_file_name
    !
    return
    !
  end subroutine stla_write
  !+-------------------------------------------------------------------+
  !| The end of the function stla_write.                               |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !
  !! GET_UNIT returns a free FORTRAN unit number.
  !
  !  Discussion:
  !
  !    A "free" FORTRAN unit number is an integer between 1 and 99 which
  !    is not currently associated with an I/O device.  A free FORTRAN unit
  !    number is needed in order to open a file with the OPEN command.
  !
  !    If IUNIT = 0, then no free FORTRAN unit could be found, although
  !    all 99 units were checked (except for units 5, 6 and 9, which
  !    are commonly reserved for console I/O).
  !
  !    Otherwise, IUNIT is an integer between 1 and 99, representing a
  !    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
  !    are special, and will never return those values.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    18 September 2005
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Output, integer ( kind = 4 ) IUNIT, the free unit number.
  !+-------------------------------------------------------------------+
  integer function get_unit()
    !
    implicit none

    integer ::  i
    integer ::  ios
    integer ::  iunit
    logical :: lopen

    get_unit = 0
    !
    i=20
    lopen=.true.
    do while(lopen)
      !
      inquire ( unit = i, opened = lopen, iostat = ios )

      if ( ios == 0  .and. (.not. lopen) ) then
        get_unit = i
      elseif(ios .ne. 0) then
        print*,' !! error with opening file unit:',i
      else
        i=i+1
      end if

    end do
    !
    ! write(*,'(A,I0,A)')'  ** file unit ',get_unit,' is activated'
    !
    return
    !
  end function get_unit
  !+-------------------------------------------------------------------+
  !| The end of the function get_unit.                                 |
  !+-------------------------------------------------------------------+
  !!
  !+-------------------------------------------------------------------+
  !
  !! CH_CAP capitalizes a single character.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    19 July 1998
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input/output, character C, the character to capitalize.
  !+-------------------------------------------------------------------+
  subroutine ch_cap ( c )
    !
    implicit none

    character c
    integer ( kind = 4 ) itemp

    itemp = ichar ( c )

    if ( 97 <= itemp .and. itemp <= 122 ) then
      c = char ( itemp - 32 )
    end if

    return
  end subroutine ch_cap
  !+-------------------------------------------------------------------+
  !| The end of the subroutine ch_cap.                                 |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !
  !! CH_EQI is a case insensitive comparison of two characters for equality.
  !
  !  Example:
  !
  !    CH_EQI ( 'A', 'a' ) is TRUE.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    28 July 2000
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character C1, C2, the characters to compare.
  !
  !    Output, logical CH_EQI, the result of the comparison.
  !+-------------------------------------------------------------------+
  logical function ch_eqi ( c1, c2 ) 
    !
    implicit none

    character c1
    character c1_cap
    character c2
    character c2_cap

    c1_cap = c1
    c2_cap = c2

    call ch_cap ( c1_cap )
    call ch_cap ( c2_cap )

    if ( c1_cap == c2_cap ) then
      ch_eqi = .true.
    else
      ch_eqi = .false.
    end if

    return
    !
  end function ch_eqi
  !+-------------------------------------------------------------------+
  !| The end of the subroutine ch_eqi.                                 |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !
  !! S_EQI is a case insensitive comparison of two strings for equality.
  !
  !  Example:
  !
  !    S_EQI ( 'Anjana', 'ANJANA' ) is TRUE.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    14 April 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S1, S2, the strings to compare.
  !
  !    Output, logical S_EQI, the result of the comparison.
  !+-------------------------------------------------------------------+
  logical function s_eqi ( s1, s2 )
    !
    implicit none

    character c1
    character c2
    integer ( kind = 4 ) i
    integer ( kind = 4 ) len1
    integer ( kind = 4 ) len2
    integer ( kind = 4 ) lenc
    character ( len = * ) s1
    character ( len = * ) s2

    len1 = len ( s1 )
    len2 = len ( s2 )
    lenc = min ( len1, len2 )

    s_eqi = .false.

    do i = 1, lenc

      c1 = s1(i:i)
      c2 = s2(i:i)
      call ch_cap ( c1 )
      call ch_cap ( c2 )

      if ( c1 /= c2 ) then
        return
      end if

    end do

    do i = lenc + 1, len1
      if ( s1(i:i) /= ' ' ) then
        return
      end if
    end do

    do i = lenc + 1, len2
      if ( s2(i:i) /= ' ' ) then
        return
      end if
    end do

    s_eqi = .true.

    return
  end function s_eqi
  !+-------------------------------------------------------------------+
  !| The end of the function s_eqi.                                    |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !
  !! WORD_NEXT_READ "reads" words from a string, one at a time.
  !
  !  Special cases:
  !
  !    The following characters are considered to be a single word,
  !    whether surrounded by spaces or not:
  !
  !      " ( ) { } [ ]
  !
  !    Also, if there is a trailing comma on the word, it is stripped off.
  !    This is to facilitate the reading of lists.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    23 May 2001
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S, a string, presumably containing words
  !    separated by spaces.
  !
  !    Output, character ( len = * ) WORD.
  !    If DONE is FALSE, then WORD contains the "next" word read.
  !    If DONE is TRUE, then WORD is blank, because there was no more to read.
  !
  !    Input/output, logical DONE.
  !    On input with a fresh string, set DONE to TRUE.
  !    On output, the routine sets DONE:
  !      FALSE if another word was read,
  !      TRUE if no more words could be read.
  !+-------------------------------------------------------------------+
  subroutine word_next_read ( s, word, done )

    implicit none

    logical done
    integer ( kind = 4 ) ilo
    integer ( kind = 4 ), save :: lenc = 0
    integer ( kind = 4 ), save :: next = 1
    character ( len = * ) s
    character, parameter :: TAB = char ( 9 )
    character ( len = * ) word
  !
  !  We "remember" LENC and NEXT from the previous call.
  !
  !  An input value of DONE = TRUE signals a new line of text to examine.
  !
    if ( done ) then

      next = 1
      done = .false.
      lenc = len_trim ( s )

      if ( lenc <= 0 ) then
        done = .true.
        word = ' '
        return
      end if

    end if
  !
  !  Beginning at index NEXT, search the string for the next nonblank,
  !  which signals the beginning of a word.
  !
    ilo = next
  !
  !  ...S(NEXT:) is blank.  Return with WORD = ' ' and DONE = TRUE.
  !
    do

      if ( lenc < ilo ) then
        word = ' '
        done = .true.
        next = lenc + 1
        return
      end if
  !
  !  If the current character is blank, skip to the next one.
  !
      if ( s(ilo:ilo) /= ' ' .and. s(ilo:ilo) /= TAB ) then
        exit
      end if

      ilo = ilo + 1

    end do
  !
  !  ILO is the index of the next nonblank character in the string.
  !
  !  If this initial nonblank is a special character,
  !  then that's the whole word as far as we're concerned,
  !  so return immediately.
  !
    if ( s(ilo:ilo) == '"' .or. &
         s(ilo:ilo) == '(' .or. &
         s(ilo:ilo) == ')' .or. &
         s(ilo:ilo) == '{' .or. &
         s(ilo:ilo) == '}' .or. &
         s(ilo:ilo) == '[' .or. &
         s(ilo:ilo) == ']' ) then

      word = s(ilo:ilo)
      next = ilo + 1
      return

    end if
  !
  !  Now search for the last contiguous character that is not a
  !  blank, TAB, or special character.
  !
    next = ilo + 1

    do while ( next <= lenc )

      if ( s(next:next) == ' ' ) then
        exit
      else if ( s(next:next) == TAB ) then
        exit
      else if ( s(next:next) == '"' ) then
        exit
      else if ( s(next:next) == '(' ) then
        exit
      else if ( s(next:next) == ')' ) then
        exit
      else if ( s(next:next) == '{' ) then
        exit
      else if ( s(next:next) == '}' ) then
        exit
      else if ( s(next:next) == '[' ) then
        exit
      else if ( s(next:next) == ']' ) then
        exit
      end if

      next = next + 1

    end do

    if ( s(next-1:next-1) == ',' ) then
      word = s(ilo:next-2)
    else
      word = s(ilo:next-1)
    end if

    return
    !
  end subroutine word_next_read
  !+-------------------------------------------------------------------+
  !| The end of the function word_next_read.                           |
  !+-------------------------------------------------------------------+
  !
  subroutine s_cat ( s1, s2, s3 )

  !*****************************************************************************80
  !
  !! S_CAT concatenates two strings to make a third string.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    18 September 2000
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S1, the "prefix" string.
  !
  !    Input, character ( len = * ) S2, the "postfix" string.
  !
  !    Output, character ( len = * ) S3, the string made by
  !    concatenating S1 and S2, ignoring any trailing blanks.
  !
    implicit none

    character ( len = * ) s1
    character ( len = * ) s2
    character ( len = * ) s3

    if ( s1 == ' ' .and. s2 == ' ' ) then
      s3 = ' '
    else if ( s1 == ' ' ) then
      s3 = s2
    else if ( s2 == ' ' ) then
      s3 = s1
    else
      s3 = trim ( s1 ) // trim ( s2 )
    end if

    return
  end
  !

  subroutine ch_to_digit ( c, digit )

  !*****************************************************************************80
  !
  !! CH_TO_DIGIT returns the integer value of a base 10 digit.
  !
  !  Example:
  !
  !     C   DIGIT
  !    ---  -----
  !    '0'    0
  !    '1'    1
  !    ...  ...
  !    '9'    9
  !    ' '    0
  !    'X'   -1
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    04 August 1999
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character C, the decimal digit, '0' through '9' or blank
  !    are legal.
  !
  !    Output, integer ( kind = 4 ) DIGIT, the corresponding integer value.  
  !    If C was 'illegal', then DIGIT is -1.
  !
    implicit none

    character c
    integer ( kind = 4 ) digit

    if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then

      digit = ichar ( c ) - 48

    else if ( c == ' ' ) then

      digit = 0

    else

      digit = -1

    end if

    return
  end

  subroutine s_to_r8 ( s, dval, ierror, length )

  !*****************************************************************************80
  !
  !! S_TO_R8 reads an R8 from a string.
  !
  !  Discussion:
  !
  !    The routine will read as many characters as possible until it reaches
  !    the end of the string, or encounters a character which cannot be
  !    part of the number.
  !
  !    Legal input is:
  !
  !       1 blanks,
  !       2 '+' or '-' sign,
  !       2.5 blanks
  !       3 integer part,
  !       4 decimal point,
  !       5 fraction part,
  !       6 'E' or 'e' or 'D' or 'd', exponent marker,
  !       7 exponent sign,
  !       8 exponent integer part,
  !       9 exponent decimal point,
  !      10 exponent fraction part,
  !      11 blanks,
  !      12 final comma or semicolon,
  !
  !    with most quantities optional.
  !
  !  Example:
  !
  !    S                 DVAL
  !
  !    '1'               1.0
  !    '     1   '       1.0
  !    '1A'              1.0
  !    '12,34,56'        12.0
  !    '  34 7'          34.0
  !    '-1E2ABCD'        -100.0
  !    '-1X2ABCD'        -1.0
  !    ' 2E-1'           0.2
  !    '23.45'           23.45
  !    '-4.2E+2'         -420.0
  !    '17d2'            1700.0
  !    '-14e-2'         -0.14
  !    'e2'              100.0
  !    '-12.73e-9.23'   -12.73 * 10.0**(-9.23)
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    07 September 2004
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) S, the string containing the
  !    data to be read.  Reading will begin at position 1 and
  !    terminate at the end of the string, or when no more
  !    characters can be read to form a legal real.  Blanks,
  !    commas, or other nonnumeric data will, in particular,
  !    cause the conversion to halt.
  !
  !    Output, real ( kind = 8 ) DVAL, the value read from the string.
  !
  !    Output, integer ( kind = 4 ) IERROR, error flag.
  !    0, no errors occurred.
  !    1, 2, 6 or 7, the input number was garbled.  The
  !    value of IERROR is the last type of input successfully
  !    read.  For instance, 1 means initial blanks, 2 means
  !    a plus or minus sign, and so on.
  !
  !    Output, integer ( kind = 4 ) LENGTH, the number of characters read
  !    to form the number, including any terminating
  !    characters such as a trailing comma or blanks.
  !
    implicit none

    character c
    real ( kind = 8 ) dval
    integer ( kind = 4 ) ierror
    integer ( kind = 4 ) ihave
    integer ( kind = 4 ) isgn
    integer ( kind = 4 ) iterm
    integer ( kind = 4 ) jbot
    integer ( kind = 4 ) jsgn
    integer ( kind = 4 ) jtop
    integer ( kind = 4 ) length
    integer ( kind = 4 ) nchar
    integer ( kind = 4 ) ndig
    real ( kind = 8 ) rbot
    real ( kind = 8 ) rexp
    real ( kind = 8 ) rtop
    character ( len = * ) s

    nchar = len_trim ( s )

    ierror = 0
    dval = 0.0D+00
    length = -1
    isgn = 1
    rtop = 0
    rbot = 1
    jsgn = 1
    jtop = 0
    jbot = 1
    ihave = 1
    iterm = 0

    do

      length = length + 1

      if ( nchar < length+1 ) then
        exit
      end if

      c = s(length+1:length+1)
  !
  !  Blank character.
  !
      if ( c == ' ' ) then

        if ( ihave == 2 ) then

        else if ( ihave == 6 .or. ihave == 7 ) then
          iterm = 1
        else if ( 1 < ihave ) then
          ihave = 11
        end if
  !
  !  Comma.
  !
      else if ( c == ',' .or. c == ';' ) then

        if ( ihave /= 1 ) then
          iterm = 1
          ihave = 12
          length = length + 1
        end if
  !
  !  Minus sign.
  !
      else if ( c == '-' ) then

        if ( ihave == 1 ) then
          ihave = 2
          isgn = -1
        else if ( ihave == 6 ) then
          ihave = 7
          jsgn = -1
        else
          iterm = 1
        end if
  !
  !  Plus sign.
  !
      else if ( c == '+' ) then

        if ( ihave == 1 ) then
          ihave = 2
        else if ( ihave == 6 ) then
          ihave = 7
        else
          iterm = 1
        end if
  !
  !  Decimal point.
  !
      else if ( c == '.' ) then

        if ( ihave < 4 ) then
          ihave = 4
        else if ( 6 <= ihave .and. ihave <= 8 ) then
          ihave = 9
        else
          iterm = 1
        end if
  !
  !  Scientific notation exponent marker.
  !
      else if ( ch_eqi ( c, 'E' ) .or. ch_eqi ( c, 'D' ) ) then

        if ( ihave < 6 ) then
          ihave = 6
        else
          iterm = 1
        end if
  !
  !  Digit.
  !
      else if (  ihave < 11 .and. lle ( '0', c ) .and. lle ( c, '9' ) ) then

        if ( ihave <= 2 ) then
          ihave = 3
        else if ( ihave == 4 ) then
          ihave = 5
        else if ( ihave == 6 .or. ihave == 7 ) then
          ihave = 8
        else if ( ihave == 9 ) then
          ihave = 10
        end if

        call ch_to_digit ( c, ndig )

        if ( ihave == 3 ) then
          rtop = 10.0D+00 * rtop + real ( ndig, kind = 8 )
        else if ( ihave == 5 ) then
          rtop = 10.0D+00 * rtop + real ( ndig, kind = 8 )
          rbot = 10.0D+00 * rbot
        else if ( ihave == 8 ) then
          jtop = 10 * jtop + ndig
        else if ( ihave == 10 ) then
          jtop = 10 * jtop + ndig
          jbot = 10 * jbot
        end if
  !
  !  Anything else is regarded as a terminator.
  !
      else
        iterm = 1
      end if
  !
  !  If we haven't seen a terminator, and we haven't examined the
  !  entire string, go get the next character.
  !
      if ( iterm == 1 ) then
        exit
      end if

    end do
  !
  !  If we haven't seen a terminator, and we have examined the
  !  entire string, then we're done, and LENGTH is equal to NCHAR.
  !
    if ( iterm /= 1 .and. length+1 == nchar ) then
      length = nchar
    end if
  !
  !  Number seems to have terminated.  Have we got a legal number?
  !  Not if we terminated in states 1, 2, 6 or 7!
  !
    if ( ihave == 1 .or. ihave == 2 .or. ihave == 6 .or. ihave == 7 ) then
      ierror = ihave
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'S_TO_R8 - Serious error!'
      write ( *, '(a)' ) '  Illegal or nonnumeric input:'
      write ( *, '(a)' ) '    ' // trim ( s )
      return
    end if
  !
  !  Number seems OK.  Form it.
  !
    if ( jtop == 0 ) then
      rexp = 1.0D+00
    else
      if ( jbot == 1 ) then
        rexp = 10.0D+00 ** ( jsgn * jtop )
      else
        rexp = 10.0D+00 ** ( real ( jsgn * jtop, kind = 8 ) &
          / real ( jbot, kind = 8 ) )
      end if
    end if

    dval = real ( isgn, kind = 8 ) * rexp * rtop / rbot

    return
  end
  subroutine stla_size ( input_file_name, solid_num, node_num, face_num, &
    text_num )

  !*****************************************************************************80
  !
  !! STLA_SIZE determines sizes associated with an STLA file.
  !
  !  Discussion:
  !
  !    This routine assumes that the file is a legal STLA file.
  !
  !    To perform checks on the file, call STLA_CHECK first.
  !
  !    Note that the counts for the number of nodes and edges are
  !    overestimates, since presumably, most nodes will be defined several
  !    times, once for each face they are part of, and most edges will
  !    be defined twice.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    15 February 2007
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Reference:
  !
  !    3D Systems, Inc,
  !    Stereolithography Interface Specification,
  !    October 1989.
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) INPUT_FILE_NAME, the name of the input file.
  !
  !    Output, integer ( kind = 4 ) SOLID_NUM, the number of solids defined.
  !    Presumably, this is 1.
  !
  !    Output, integer ( kind = 4 ) NODE_NUM, the number of vertices defined.
  !
  !    Output, integer ( kind = 4 ) FACE_NUM, the number of faces defined.
  !
  !    Output, integer ( kind = 4 ) TEXT_NUM, the number of lines of text 
  !    in the file.
  !
    implicit none

    logical done
    real ( kind = 8 ) dval
    integer ( kind = 4 ) face_num
    integer ( kind = 4 ) i
    integer ( kind = 4 ) ierror
    character ( len = * ) input_file_name
    integer ( kind = 4 ) ios
    integer ( kind = 4 ) iunit
    integer ( kind = 4 ) lchar
    integer ( kind = 4 ) node_num
    integer ( kind = 4 ) solid_num
    integer ( kind = 4 ) state
    character ( len = 255 ) text
    integer ( kind = 4 ) text_num
    integer ( kind = 4 ) vertex
    character ( len = 255 ) word1
    character ( len = 255 ) word2

    ierror = 0

    state = 0
    text_num = 0

    solid_num = 0
    node_num = 0
    face_num = 0
  !
  !  Open the file.
  !
    iunit=get_unit()

    open ( unit = iunit, file = input_file_name, status = 'old', iostat = ios )

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'STLA_SIZE - Fatal error!'
      write ( *, '(a)' ) '  Could not open the file "' &
        // trim ( input_file_name ) // '".'
      ierror = 1
      return
    end if
  !
  !  Read the next line of text.
  !
    do

      read ( iunit, '(a)', iostat = ios ) text

      if ( ios /= 0 ) then
        if ( state /= 0 .and. state /= 1 ) then
          return
        end if
        exit
      end if

      text_num = text_num + 1

      done = .true.
  !
  !  Read the first word in the line.
  !
      call word_next_read ( text, word1, done )

      if ( done ) then
        return
      end if
  !
  !  "Doctor" the text, changing a beginning occurrence of:
  !
  !      END FACET to ENDFACET
  !      END LOOP to ENDLOOP
  !      END SOLID to ENDSOLID
  !      FACET NORMAL to FACETNORMAL
  !      OUTER LOOP to OUTERLOOP
  !
      if ( s_eqi ( word1, 'END' ) ) then

        call word_next_read ( text, word2, done )

        if ( .not. s_eqi ( word2, 'FACET' ) .and. &
             .not. s_eqi ( word2, 'LOOP' ) .and. &
             .not. s_eqi ( word2, 'SOLID' ) ) then
          return
        end if

        call s_cat ( word1, word2, word1 )

      else if ( s_eqi ( word1, 'FACET' ) ) then

        call word_next_read ( text, word2, done )

        if ( .not. s_eqi ( word2, 'NORMAL' ) ) then
          return
        end if

        call s_cat ( word1, word2, word1 )

      else if ( s_eqi ( word1, 'OUTER' ) ) then

        call word_next_read ( text, word2, done )

        if ( .not. s_eqi ( word2, 'LOOP' ) ) then
          return
        end if

        call s_cat ( word1, word2, word1 )

      end if
  !
  !  This first word tells us what to do.
  !
  !  SOLID - begin a new solid.
  !    Valid in state 0, moves to state 1.
  !  ENDSOLID - end current solid.
  !    Valid in state 1, moves to state 0.
  !
  !  FACETNORMAL - begin a new facet.
  !    Valid in state 0 or 1, moves to state 2.
  !  ENDFACET - end current facet.
  !    Valid in state 2, moves to state 1.
  !
  !  OUTERLOOP - begin a list of vertices.
  !    Valid in state 2, moves to state 3.
  !  ENDLOOP - end vertex list.
  !    Valid in state 3, moves to state 2.
  !
  !  VERTEX - give coordinates of next vertex.
  !    Valid in state 3.
  !
  !  End of file -
  !    Valid in state 0 or 1.
  !
      if ( s_eqi ( word1, 'SOLID' ) ) then

        if ( state /= 0 ) then
          return
        end if

        state = 1

      else if ( s_eqi ( word1, 'ENDSOLID' ) ) then

        if ( state /= 1 ) then
          return
        end if

        state = 0

        solid_num = solid_num + 1

      else if ( s_eqi ( word1, 'FACETNORMAL' ) ) then

        if ( state /= 0 .and. state /= 1 ) then
          return
        end if

        state = 2

        do i = 1, 3

          call word_next_read ( text, word2, done )

          if ( done ) then
            return
          end if

          call s_to_r8 ( word2, dval, ierror, lchar )

          if ( ierror /= 0 ) then
            return
          end if

        end do

      else if ( s_eqi ( word1, 'ENDFACET' ) ) then

        if ( state /= 2 ) then
          return
        end if

        state = 1

        face_num = face_num + 1

      else if ( s_eqi ( word1, 'OUTERLOOP' ) ) then

        if ( state /= 2 ) then
          return
        end if

        state = 3
        vertex = 0

      else if ( s_eqi ( word1, 'ENDLOOP' ) ) then

        if ( state /= 3 ) then
          return
        end if

        state = 2

      else if ( s_eqi ( word1, 'VERTEX' ) ) then

        if ( state /= 3 ) then
          return
        end if

        if ( 3 <= vertex ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'STLA_SIZE - Fatal error!'
          write ( *, '(a)' ) '  Too many vertices for a face.'
          ierror = 1
          return
        end if

        do i = 1, 3

          call word_next_read ( text, word2, done )

          if ( done ) then
            return
          end if

          call s_to_r8 ( word2, dval, ierror, lchar )

          if ( ierror /= 0 ) then
            return
          end if

        end do

        vertex = vertex + 1
        node_num = node_num + 1

      else

        return

      end if

    end do
  !
  !  Close the file.
  !
    close ( unit = iunit )

    return
  end
  subroutine stla_size_print ( input_file_name, solid_num, node_num, face_num, &
    text_num )

  !*****************************************************************************80
  !
  !! STLA_SIZE_PRINT prints sizes associated with an STLA file.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  Modified:
  !
  !    15 February 2007
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Reference:
  !
  !    3D Systems, Inc,
  !    Stereolithography Interface Specification,
  !    October 1989.
  !
  !  Parameters:
  !
  !    Input, character ( len = * ) INPUT_FILE_NAME, the name of the input file.
  !
  !    Input, integer ( kind = 4 ) SOLID_NUM, the number of solids defined.
  !
  !    Input, integer ( kind = 4 ) NODE_NUM, the number of vertices defined.
  !
  !    Input, integer ( kind = 4 ) FACE_NUM, the number of faces defined.
  !
  !    Input, integer ( kind = 4 ) TEXT_NUM, the number of lines of text 
  !    in the file.
  !
    implicit none

    integer ( kind = 4 ) face_num
    character ( len = * ) input_file_name
    integer ( kind = 4 ) node_num
    integer ( kind = 4 ) solid_num
    integer ( kind = 4 ) text_num

    write ( *, '(a)'    ) ' '
    write ( *, '(a)'    ) '  Object sizes for STLA file "' // &
      trim ( input_file_name ) // '":'
    write ( *, '(a)'    ) ' '
    write ( *, '(a,i8)' ) '  Solids =                   ', solid_num
    write ( *, '(a,i8)' ) '  Nodes (may be repeated) =  ', node_num
    write ( *, '(a,i8)' ) '  Faces (triangular only) =  ', face_num
    write ( *, '(a)'    ) ' '
    write ( *, '(a,i8)' ) '  Number of lines of text =  ', text_num

    return
  end
  !
end module stlaio
!+---------------------------------------------------------------------+
!| The end of the module stlaio.                                       |
!+---------------------------------------------------------------------+