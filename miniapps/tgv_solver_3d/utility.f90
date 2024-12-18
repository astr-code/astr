!+---------------------------------------------------------------------+
!| This module contains utility subroutines                            |
!+---------------------------------------------------------------------+
!| CHANGE RECORD                                                       |
!| -------------                                                       |
!| 10-08-2022  | Created by J. Fang                                    |
!+---------------------------------------------------------------------+
module utility

   implicit none

   contains

   !+-------------------------------------------------------------------+
   !| This subroutine is used to print welcome infomation.              |
   !+-------------------------------------------------------------------+
   !| CHANGE RECORD                                                     |
   !| -------------                                                     |
   !| 08-Oct-2018  | Created by J. Fang @ Warrington                    |
   !+-------------------------------------------------------------------+
   subroutine statement
      !
      write (*, *)
      write (*, *)
      print *, ' +------------------------ Statement -------------------------+'
      print *, ' |                                                            |'
      print *, ' |      A mini application of astr by Jian Fang since 2008    |'
      print *, ' |               <ASTR> Copyright Resvered <2008>             |'
      print *, ' |         Advanced Simulatior for Turbulence Research        |'
      print *, ' |                                                            |'
      print *, ' +------------------------------------------------------------+'
      print *, ' |                                                            |'
      print *, ' | Copyright 2008 Jian Fang                                   |'
      print *, ' |                                                            |'
      print *, ' | Licensed under the Apache License, Version 2.0             |'
      print *, ' | you may not use this file except in compliance with the    |'
      print *, ' | license. You may obtain a copy of the License at           |'
      print *, ' |                                                            |'
      print *, ' |        http://www.apache.org/licenses/LICENSE-2.0          |'
      print *, ' |                                                            |'
      print *, ' | Unless required by applicable law or agreed to in writing, |'
      print *, ' | software distributed under the License is distributed on an|'
      print *, ' | "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY     |'
      print *, ' | KIND, either express or implied. See the License for the   |'
      print *, ' | specific language governing permissions and limitations    |'
      print *, ' | under the License.                                         |'
      print *, ' |                                                            |'
      print *, ' +--------------------- End of Statement ---------------------+'
      write (*, *)
      !
      write (*, *) '                        ___   _____________  '
      write (*, *) '                       / _ | / __/_  __/ _ \ '
      write (*, *) '                      / __ |_\ \  / / / , _/ '
      write (*, *) '                     /_/ |_/___/ /_/ /_/|_|  '
      write (*, *)
      write (*, *)
      !
      return
      !
   end subroutine statement
   !+-------------------------------------------------------------------+
   !| The end of the subroutine statement.                              |
   !+-------------------------------------------------------------------+

  
  !+-------------------------------------------------------------------+
  !| This subroutine is used to input command from keyboard.           |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 11-02-2021  | Created by J. Fang @ Warrington                     |
  !| 28-05-2021  | Moved to the module command by J. Fang @ Warrington |
  !+-------------------------------------------------------------------+
  subroutine readkeyboad(keyin)
    !
    character(len=*),intent(out),optional :: keyin
    !
    ! local data
    integer :: ierr,cli_len,nlen,arg_count
    integer,save :: nkey=0
    !
    nkey=nkey+1
    call get_command_argument(nkey,keyin,cli_len,ierr)
    !
  end subroutine readkeyboad
  !+-------------------------------------------------------------------+
  !| The end of the subroutine readkeyboad.                            |
  !+-------------------------------------------------------------------+
  !
  subroutine timereport
    !
    use comvardef, only: ctime
    !
    open(14,file='time_report.txt')
    write(14,'(A)')'---------cpu time cost---------'
    write(14,'(A)')'-------------------------------'
    write(14,'(A,F13.3,A)')'   total time |',ctime(1),' s'
    write(14,'(A,F13.3,A)')'          rk3 |',ctime(2),' s'
    write(14,'(A,F13.3,A)')'       rhscal |',ctime(7),' s'
    write(14,'(A,F13.3,A)')'      gradcal |',ctime(3),' s'
    write(14,'(A,F13.3,A)')'   convection |',ctime(4),' s'
    write(14,'(A,F13.3,A)')'    diffusion |',ctime(5),' s'
    write(14,'(A,F13.3,A)')'      filterq |',ctime(6),' s'
    write(14,'(A,F13.3,A)')'      diff6ec |',ctime(8),' s'
    write(14,'(A)')'-------------------------------'
    close(14)
    print*,' << time_report.txt'
    !
  end subroutine timereport

   !+-------------------------------------------------------------------+
   !| This function Verifies that a character string represents a       |
   !|  numerical value                                                  |
   !+-------------------------------------------------------------------+
   ! ref: http://fcode.cn/code_gen-115-1.html
   !+-------------------------------------------------------------------+
   Integer Function IsNum(zval)
      ! 确定字符是否是数值类型：
      ! 0-非数值的字符串
      ! 1-整数(integer)
      ! 2-小数(fixed point real)
      ! 3-指数类型实数(exponent type real)
      ! 4-双精度实数指数形式(exponent type double)
      Character(Len=*), Intent(In) :: zval
      !
      Integer :: num, nmts, nexp, kmts, ifexp, ichr
      !
      Integer, Parameter :: kint = 1 ! integer
      Integer, Parameter :: kfix = 2 ! fixed point real
      Integer, Parameter :: kexp = 3 ! exponent type real
      Integer, Parameter :: kdbl = 4 ! exponent type double
      !
      ! initialise
      num = 0  ! 数字的格式，最后传递给ISNUM返回
      nmts = 0 ! 整数或浮点数的数字个数
      nexp = 0 ! 指数形式的数字个数
      kmts = 0 ! 有+-号为1，否则为0
      ifexp = 0! 似乎没用
      ! loop over characters
      ichr = 0
      !
      Do

         If (ichr >= len(zval)) Then

            ! last check

            If (nmts == 0) Exit

            If (num >= kexp .And. nexp == 0) Exit

            isnum = num

            Return

         End If

         ichr = ichr + 1

         Select Case (zval(ichr:ichr))

            ! process blanks

         Case (' ')

            Continue

            ! process digits

         Case ('0', '1', '2', '3', '4', '5', '6', '7', '8', '9')

            If (num == 0) num = kint

            If (num < kexp) Then

               nmts = nmts + 1

               ! 整数或浮点数+1

            Else

               nexp = nexp + 1

               ! 指数形式+1

            End If

            ! process signs

         Case ('+', '-')

            If (num == 0) Then

               If (kmts > 0) Exit

               ! 出现2个符号，非数字

               kmts = 1

               num = kint

            Else

               If (num < kexp) Exit

               If (ifexp > 0) Exit

               ifexp = 1

            End If

            ! process decimal point

         Case ('.')

            If (num /= kint .And. ichr /= 1) Exit

            ! 前面不是整数，小数点也不是第一个字符，则非数字

            num = kfix

            ! process exponent

         Case ('e', 'E')

            If (num >= kexp) Exit

            If (nmts == 0) Exit

            num = kexp
         Case ('d', 'D')

            If (num >= kexp) Exit

            If (nmts == 0) Exit

            num = kdbl

            ! any other character means the string is non-numeric

         Case Default

            Exit

         End Select

      End Do

      ! if this point is reached, the string is non-numeric

      isnum = 0

      Return

   End Function IsNum
   !+-------------------------------------------------------------------+
   !| The end of the Function IsNum.                                    |
   !+-------------------------------------------------------------------+
   !
   function rnorm_box_muller(mode) result(variates)
      !
      use constdef
      ! coded formulas from https://en.wikipedia.org/wiki/Box%E2%80%93Muller_transform
      !
      character(len=*), intent(in) :: mode
      !
      ! return two uncorrelated standard normal variates
      integer, allocatable :: seed(:)
      integer :: rantime(8)
      !
      integer :: n
      real(8) :: variates(2)
      real(8) :: u(2), factor, arg
      !
      logical, save :: firstcall = .true.
      !
      if (mode == 'sync') then
         ! all processor generate same random numbers
         if (firstcall) then
            !
            call random_seed(size=n)
            allocate (seed(n))
            call date_and_time(values=rantime)
            !  use date and minutes for synthetisation
            !    1      2    3    4     5      6      7    8
            !-----------------------------------------------
            ! 2023     10    5   60    23      6     28  962
            !-----------------------------------------------
            ! year  month date   ??  hour minute second msec
            !
            seed = 0
            seed(1:6) = rantime(1:6)
            !
            call random_seed(put=seed)
            !
            deallocate (seed)
            !
            firstcall = .false.
         end if
         !
      end if
      !
      do
         call random_number(u)
         if (u(1) > 0.d0) exit
      end do
      factor = sqrt(-2*log(u(1)))
      arg = 2.d0*pi*u(2)
      variates = factor*[cos(arg), sin(arg)]
      !
   end function rnorm_box_muller
   !
   !+-------------------------------------------------------------------+
   !| Progress indicators library.                                      |
   !+-------------------------------------------------------------------+
   !| CHANGE RECORD                                                     |
   !| -------------                                                     |
   !| 31-01-2024  | copied by J. Fang via:                              |
   !| https://github.com/macie/fortran-libs                             |
   !|  Maciej Żok, 2010 MIT License                                     |
   !+-------------------------------------------------------------------+
   subroutine progress_bar(iteration, maximum, info2show, barlength)
      !
      ! Prints progress bar.
      !
      ! Args:
      !     iteration - iteration number
      !     maximum - total iterations
      !     barlength - length of the bar
      !
      ! use iso_fortran_env
      integer, intent(in) :: iteration, maximum
      character(len=*), intent(in), optional :: info2show
      integer, intent(in), optional :: barlength
      integer :: counter, nlength
      integer :: done
      real(4) :: perc
      !
      if (present(barlength)) then
         nlength = barlength
      else
         nlength = 10
      end if
      !
      perc = 100.0*real(iteration)/real(maximum)
      done = floor(perc/(100.0/real(nlength)))  ! mark length
      !
      write (6, '(1A1,A,A)', advance='no') char(13), info2show, '['
      if (done .LE. 0) then
         do counter = 1, nlength
            write (6, '(1A1,A)', advance='no') '='
         end do
      else if ((done .GT. 0) .and. (done .LT. nlength)) then
         do counter = 1, done
            write (6, '(1A1,A)', advance='no') '>'
         end do
         do counter = done + 1, nlength
            write (6, '(1A1,A)', advance='no') '='
         end do
      else
         do counter = 1, nlength
            write (6, '(1A1,A)', advance='no') '>'
         end do
      end if
      write (6, '(A,F5.1,A)', advance='no') '] ', perc, '%'
      !
      if (iteration == maximum) write (6, *)
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

      integer :: n, i, istart

      ! Initialize counters

      ! Count number of tokens for allocation
      n = count_tokens(line)

      ! Allocate arrays based on the number of tokens
      allocate (strings(n))

      ! Split line into tokens and categorize them
      istart = 1
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
         token = trim(adjustl(token))//line(i:i)
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
      i = 20
      lopen = .true.
      do while (lopen)
         !
         inquire (unit=i, opened=lopen, iostat=ios)

         if (ios == 0 .and. (.not. lopen)) then
            get_unit = i
         elseif (ios .ne. 0) then
            print *, ' !! error with opening file unit:', i
         else
            i = i + 1
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
end module utility
!+---------------------------------------------------------------------+
!| The end of the module utility                                       |
!+---------------------------------------------------------------------+
