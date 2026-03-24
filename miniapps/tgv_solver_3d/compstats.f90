program merge_files
    use ieee_arithmetic
    ! 
    implicit none
    ! 
    character(len=40) :: line1, line2
    integer :: ios1, ios2, i1, i2
    real :: num1, num2, num3
    logical :: flag = .true.
    character(len=10) :: second_last_word
    integer :: pos_last_space, pos_second_last_space, ln

    ln = 0

    open(unit=10, file='time_report_cpu.txt', status='old', action='read')
    open(unit=20, file='time_report_gpu.txt', status='old', action='read')
    open(unit=30, file='time_report.txt', status='replace', action='write')

    do
        read(10, '(A)', iostat=ios1) line1
        read(20, '(A)', iostat=ios2) line2

        if (ios1 /= 0 .or. ios2 /= 0) then
            write(30, '(A, A, A)') line1, line2, '---------------------'
            exit
        end if

        ln = ln + 1
        flag = .true.
        
        ! read(line1, '(F13.3)', IOSTAT=i) num1
        ! read(line1, '(F13.3)', IOSTAT=i) num2
        ! num3 = num1/num2

        pos_last_space = len_trim(line1)
        do while (pos_last_space > 1 .and. line1(pos_last_space:pos_last_space) /= ' ')
          pos_last_space = pos_last_space - 1
        end do
      
        pos_second_last_space = pos_last_space - 1
        do while (pos_second_last_space > 1 .and. line1(pos_second_last_space:pos_second_last_space) /= ' ')
          pos_second_last_space = pos_second_last_space - 1
        end do
      
        second_last_word = adjustl(line1(pos_second_last_space+1:pos_last_space-1))
        read(second_last_word, '(F13.3)', iostat=i1) num1

        ! print *, line1, second_last_word

        pos_last_space = len_trim(line2)
        do while (pos_last_space > 1 .and. line2(pos_last_space:pos_last_space) /= ' ')
          pos_last_space = pos_last_space - 1
        end do
      
        pos_second_last_space = pos_last_space - 1
        do while (pos_second_last_space > 1 .and. line2(pos_second_last_space:pos_second_last_space) /= ' ')
          pos_second_last_space = pos_second_last_space - 1
        end do
      
        second_last_word = adjustl(line2(pos_second_last_space+1:pos_last_space-1))
        read(second_last_word, '(F13.3)', iostat=i2) num2

        if (ln == 1) then
            write(30, '(A, A, A)') line1, line2, '-----acceleration----'
            cycle
        end if

        if (ln == 2) then
            write(30, '(A, A, A)') line1, line2, '---------------------'
            cycle
        end if

        if(num1 == 0 .or. num2 == 0) flag = .false.
        ! print *, num1, num2, flag

        if (flag) then
            num3 = num1/num2
            write(30, '(A, A, F13.1)') line1, line2, num3
        else 
            write(30, '(A, A, A)') line1, line2
        end if

        if (line1 == '.end' .or. line2 == '.end') print *, "end", line1, line2
    end do
    
    close(10)
    close(20)
    close(30)

    print*,'-- acceleration statistiscs report generated : time_report.txt --'
end program merge_files
