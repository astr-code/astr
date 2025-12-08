module pastr_utility

    implicit none

    private

    public :: get_cmd_argument,isnum,make_dir,sort_unique,add_new_element,remove_element

    interface sort_unique
      module procedure  :: sort_unique_int
    end interface sort_unique
    interface add_new_element
      module procedure  :: add_new_element_int
      module procedure  :: add_new_element_patch
      module procedure  :: add_new_element_linker
    end interface add_new_element
    interface remove_element
      module procedure  :: remove_element_linker
    end interface remove_element

contains

  subroutine make_dir(dirname)

      implicit none
      character(*), intent(in) :: dirname
      integer :: status
      logical :: exists
  
      inquire(file=trim(dirname), exist=exists)
      if (.not. exists) then
      
        call execute_command_line("mkdir -p " // trim(dirname), exitstat=status)
  
        if (status /= 0) then
            print *, "Error: could not create directory: ", trim(dirname)
        end if

      endif

  end subroutine make_dir

  !+-------------------------------------------------------------------+
  !|Description. Returns a command argument.                           |
  !|                                                                   |
  !|Class. Subroutine.                                                 |
  !|                                                                   |
  !|Arguments.                                                         |
  !|NUMBER shall be scalar and of type default integer. It is an       |
  !|  INTENT(IN) argument. It specifies the number of the command      |
  !|  argument that the other arguments give information about. Useful |
  !|  values of NUMBER are those between 0 and the argument count      |
  !|  returned by the COMMAND_ARGUMENT_COUNT intrinsic.                |
  !|  Other values are allowed, but will result in error status return |
  !|  (see below).  Command argument 0 is defined to be the command    |
  !|  name by which the program was invoked if the processor has such  |
  !|  a concept. It is allowed to call the get_cmd_argument        |
  !|  procedure for command argument number 0, even if the processor   |
  !|  does not define command names or other command arguments.        |
  !|  The remaining command arguments are numbered consecutively from  |
  !|  1 to the argument count in an order determined by the processor. |
  !|VALUE (optional) shall be scalar and of type default character.    |
  !|  It is an INTENT(OUT) argument. It is assigned the value of the   |
  !|  command argument specified by NUMBER. If the command argument    |
  !|  value cannot be determined, VALUE is assigned all blanks.        |
  !|LENGTH (optional) shall be scalar and of type default integer.     |
  !|  It is an INTENT(OUT) argument. It is assigned the significant    |
  !|  length of the command argument specified by NUMBER. The          |
  !|  significant length may include trailing blanks if the processor  |
  !|   allows command arguments with significant trailing blanks. This |
  !|  length does not consider any possible truncation or padding in   |
  !|   assigning the command argument value to the VALUE argument; in  |
  !|  fact the VALUE argument need not even be present. If the command |
  !|  argument length cannot be determined, a length of 0 is assigned. |
  !|STATUS (optional) shall be scalar and of type default integer.     |
  !|  It is an INTENT(OUT) argument. It is assigned the value 0 if     |
  !|  the argument retrieval is sucessful. It is assigned a            |
  !|  processor-dependent non-zero value if the argument retrieval.    |
  !|  fails                                                            |
  !|NOTE                                                               |
  !|  One possible reason for failure is that NUMBER is negative or    |
  !|  greater than COMMAND_ARGUMENT_COUNT().                           |
  !+-------------------------------------------------------------------+
  subroutine get_cmd_argument(number,value,length,status)
    integer         , intent(in)            :: number
    character(len=*), intent(out), optional :: value
    integer         , intent(out), optional :: length
    integer         , intent(out), optional :: status
    !
    !  A temporary variable for the rare case case where LENGTH is
    !  specified but VALUE is not. An arbitrary maximum argument length
    !  of 1000 characters should cover virtually all situations.
    character(len=1000) :: tmpval
    !
    integer :: iargc
    ! Possible error codes:
    ! 1 = Argument number is less than minimum
    ! 2 = Argument number exceeds maximum
    if (number < 0) then
      if (present(value )) value  = ' '
      if (present(length)) length = 0
      if (present(status)) status = 1
      return
    else if (number > iargc()) then
      if (present(value )) value  = ' '
      if (present(length)) length = 0
      if (present(status)) status = 2
      return
    end if
    !
    ! Get the argument if VALUE is present
    !
    if (present(value)) call getarg(number,value)
    !
    ! The LENGTH option is fairly pointless under Unix.
    ! Trailing spaces can only be specified using quotes.
    ! Since the command line has already been processed by the
    ! shell before the application sees it, we have no way of
    ! knowing the true length of any quoted arguments. LEN_TRIM
    ! is used to ensure at least some sort of meaningful result.
    !
    if (present(length)) then
      if (present(value)) then
        length = len_trim(value)
      else
        call getarg(number,tmpval)
        length = len_trim(tmpval)
      end if
    end if
    !
    ! Since GETARG does not return a result code, assume success
    !
    if (present(status)) status = 0
    return
    !
  end subroutine get_cmd_argument
  !+-------------------------------------------------------------------+
  !| The end of the subroutine get_cmd_argument                        |
  !+-------------------------------------------------------------------+
  !!
  !+-------------------------------------------------------------------+
  !| This function is to to determine if a string is numeric.          |
  !+-------------------------------------------------------------------+
  !| ref: http://fcode.cn/code_gen-115-1.html                          |
  !+-------------------------------------------------------------------+
  integer function isnum(zval)
    !
    ! Verify that a character string represents a numerical value
    ! 确定字符是否是数值类型：
    !   0-非数值的字符串
    !   1-整数(integer)
    !   2-小数(fixed point real)
    !   3-指数类型实数(exponent type real)
    !   4-双精度实数指数形式(exponent type double)
    !
    character(len=*),intent (in) :: zval
    integer :: num, nmts, nexp, kmts, ifexp, ichr
    integer, parameter :: kint = 1 ! integer
    integer, parameter :: kfix = 2 ! fixed point real
    integer, parameter :: kexp = 3 ! exponent type real
    integer, parameter :: kdbl = 4 ! exponent type double
    ! 
    ! initialise
    !
    num = 0  ! 数字的格式，最后传递给isnum返回
    nmts = 0 ! 整数或浮点数的数字个数
    nexp = 0 ! 指数形式的数字个数
    kmts = 0 ! 有+-号为1，否则为0
    ifexp = 0! 似乎没用
    !
    ! loop over characters
    !
    ichr = 0
    !
    do
      !
      if(ichr>=len(zval)) then
        !
        ! last check
        !
        if (nmts==0) exit
        !
        if (num>=kexp .and. nexp==0) exit
        !
        isnum = num
        !
        return
        !
      end if
      !
      ichr = ichr + 1
      !
      select case (zval(ichr:ichr))
        ! process blanks
      case (' ')
        continue
      ! process digits
      case ('0', '1', '2', '3', '4', '5', '6', '7', '8', '9')
        if (num==0) num = kint
        if (num<kexp) then
          nmts = nmts + 1
          ! 整数或浮点数+1
        else
          nexp = nexp + 1
          ! 指数形式+1
        end if
      ! process signs
      case ('+', '-')
        !
        if (num==0) then
          if (kmts>0) exit
          ! 出现2个符号，非数字
          kmts = 1
          num = kint
        else
          if (num<kexp) exit
          if (ifexp>0) exit
          ifexp = 1
        end if
        ! process decimal point
      case ('.')
        if (num/=kint .and. ichr/=1) exit
        ! 前面不是整数，小数点也不是第一个字符，则非数字
        num = kfix
        ! process exponent
      case ('e', 'E')
        if (num>=kexp) exit
        if (nmts==0) exit
        num = kexp
      case ('d', 'D')
        if (num>=kexp) exit
        if (nmts==0) exit
        num = kdbl
        ! any other character means the string is non-numeric
      case default
        exit
      end select
      !
    end do
    !
    ! if this point is reached, the string is non-numeric
    !
    isnum = 0
    !  
    return
    !
  end function isnum
  !+-------------------------------------------------------------------+
  !| The end of the function isnum                                     |
  !+-------------------------------------------------------------------+
  
  !+-------------------------------------------------------------------+
  !| Progress indicators library.                                      |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 31-01-2024  | copied by J. Fang via:                              |
  !| https://github.com/macie/fortran-libs                             |
  !|  Maciej Żok, 2010 MIT License                                     |
  !+-------------------------------------------------------------------+
  subroutine progress_bar(iteration,maximum,info2show,barlength)
     !
     ! Prints progress bar.
     !
     ! Args: 
     !     iteration - iteration number
     !     maximum - total iterations
     !     barlength - length of the bar
     !   
     ! use iso_fortran_env
     integer,intent(in) :: iteration,maximum
     character(len=*),intent(in),optional :: info2show
     integer,intent(in),optional :: barlength
     integer :: counter,nlength
     integer :: done
     real(4) :: perc
     !
     if(present(barlength)) then
         nlength=barlength
     else
         nlength=10
     endif
     !
     perc = 100.0*real(iteration)/real(maximum)
     done = floor(perc/(100.0/real(nlength)))  ! mark length
     !
     write(6,'(1A1,A,A)',advance='no')char(13),info2show,'['
     if (done .LE. 0) then
         do counter = 1, nlength
             write(6,'(1A1,A)',advance='no')'='
         end do
     else if ((done .GT. 0) .and. (done .LT. nlength)) then
         do counter = 1, done
             write(6,'(1A1,A)',advance='no')'>'
         end do
         do counter = done+1, nlength
             write(6,'(1A1,A)',advance='no')'='
         end do 
     else
         do counter = 1, nlength
             write(6,'(1A1,A)',advance='no')'>'
         end do
     end if
     write(6,'(A,F5.1,A)',advance='no')'] ',perc,'%'
     !
     if(iteration==maximum) write(6,*)
     !
  end subroutine progress_bar
  !+-------------------------------------------------------------------+
  !| The end of the subroutine progress_bar.                           |
  !+-------------------------------------------------------------------+


  function line_2_strings(line) result(strings)

    ! arguments
    character(len=*), intent(in) :: line

    ! local data
    character(len=20), allocatable :: strings(:)
    character(len=20) :: token

    integer :: n,i,istart

    ! Initialize counters
  
    ! Count number of tokens for allocation
    n = count_tokens(line)
  
    ! Allocate arrays based on the number of tokens
    allocate(strings(n))

    ! Split line into tokens and categorize them
    istart=1
    do i = 1, n
      call get_token(line, istart, token)
      strings(i) = token
    end do

    return

  end function line_2_strings

  ! Function to count tokens in a line
  integer function count_tokens(line)
    character(len=*), intent(in) :: line
    integer :: idx, len_line, count

    count = 0
    idx = 1
    len_line = len_trim(line)

    do while (idx <= len_line)
      ! Skip delimiters
      if (line(idx:idx) == ' ' .or. line(idx:idx) == ',') then
        idx = idx + 1
      else
        count = count + 1
        ! Skip the current token
        do while (idx <= len_line .and. &
                  line(idx:idx) /= ' ' .and. line(idx:idx) /= ',')
          idx = idx + 1
        end do
      end if
    end do

    count_tokens = count

  end function count_tokens


  ! Function to get each token from the line
  subroutine get_token(line, istart, token)
    character(len=*), intent(in) :: line
    integer, intent(inout) :: istart
    character(len=*), intent(out) :: token
    integer :: i, len_line

    len_line = len_trim(line)
    token = ''
    i = istart
    do while (i <= len_line .and. (line(i:i) == ' ' .or. line(i:i) == ','))
      i = i + 1
    end do
    istart = i
    do while (i <= len_line .and. line(i:i) /= ' ' .and. line(i:i) /= ',')
      token = trim(adjustl(token)) // line(i:i)
      i = i + 1
    end do
    istart = i
  end subroutine get_token

  ! Function to check if a string is numeric
  logical function is_numeric(str)
    character(len=*), intent(in) :: str
    integer :: i, dot_count
    dot_count = 0
    is_numeric = .true.
    do i = 1, len_trim(str)
      if (str(i:i) == '.') then
        dot_count = dot_count + 1
        if (dot_count > 1) then
          is_numeric = .false.
          return
        end if
      else if (.not. (str(i:i) >= '0' .and. str(i:i) <= '9')) then
        is_numeric = .false.
        return
      end if
    end do
  end function is_numeric

  subroutine sort_unique_int(arr)
    !----------------------------------------------------------------------
    ! Sort an integer array and remove duplicated elements.
    !
    ! Input/Output:
    !   arr  : integer array, allocated; will be overwritten with unique sorted list
    !   n    : integer; input = original size, output = size after removing duplicates
    !
    ! Usage:
    !   integer, allocatable :: a(:)
    !   integer :: n
    !
    !   a = [5,2,3,3,2,10]
    !   n = size(a)
    !   call sort_unique_int(a, n)
    !   ! Now: a = [2,3,5,10], n = 4
    !----------------------------------------------------------------------
    implicit none
    integer, allocatable, intent(inout) :: arr(:)

    integer :: i, m, n
    integer, allocatable :: tmp(:)

    if(.not. allocated(arr)) return
    n=size(arr)
    if(n<2) return

    ! ---- Sort array (Insertion Sort: simple & stable) ----
    do i = 2, n
        m = arr(i)
        call insert_position(arr, i, m)
    end do

    ! ---- Remove duplicates ----
    m = 1
    do i = 2, n
        if (arr(i) /= arr(m)) then
            m = m + 1
            arr(m) = arr(i)
        end if
    end do

    ! ---- Resize array to new length ----
    if (m < n) then
        allocate(tmp(m))
        tmp = arr(1:m)
        call move_alloc(tmp, arr)
    end if

contains

    ! insertion into sorted part arr(1:i-1)
    subroutine insert_position(a, iend, value)
        integer, intent(inout) :: a(:)
        integer, intent(in) :: iend, value
        integer :: j
        j = iend - 1
        do while (j >= 1 .and. a(j) > value)
            a(j+1) = a(j)
            j = j - 1
        end do
        a(j+1) = value
    end subroutine insert_position

  end subroutine sort_unique_int

  subroutine add_new_element_int(array,element)
      
      integer,intent(inout),allocatable :: array(:)
      integer,intent(in) :: element
       
      if(.not.allocated(array)) then
        allocate(array(1))
        array(1)=element
      else
        array=[array,element]
      endif

  end subroutine add_new_element_int
  subroutine add_new_element_patch(array,element)

      use pastr_multiblock_type, only: patch_type
      
      type(patch_type),intent(inout),allocatable :: array(:)
      type(patch_type),intent(in) :: element
       
      if(.not.allocated(array)) then
        allocate(array(1))
        array(1)=element
      else
        array=[array,element]
      endif

  end subroutine add_new_element_patch
  subroutine add_new_element_linker(array,element)

      use pastr_multiblock_type, only: link_type
      
      type(link_type),intent(inout),allocatable :: array(:)
      type(link_type),intent(in) :: element
       
      if(.not.allocated(array)) then
        allocate(array(1))
        array(1)=element
      else
        array=[array,element]
      endif

  end subroutine add_new_element_linker
  
  subroutine remove_element_linker(array,idx)

    use pastr_multiblock_type, only: link_type

    type(link_type),intent(inout), allocatable :: array(:)
    integer,intent(in) :: idx
    type(link_type),allocatable :: tmp(:)

    integer :: n
    
    n = size(array)
    allocate(tmp(n-1))
    
    if (idx > 1) tmp(1:idx-1) = array(1:idx-1)
    if (idx < n) tmp(idx:n-1) = array(idx+1:n)
    
    call move_alloc(tmp, array)

  end subroutine remove_element_linker

end module pastr_utility