! kernels.f90
module kernels

  use iso_fortran_env, only: wp => real64
  
  implicit none

contains

  subroutine inout(arg1)
  
    use pastr_constdef

    character(len=*), intent(in) :: arg1
  
    character(len=30) :: pastr_command
  
    print *, 'Fortran received:', arg1
  
    print*,' pi=',pi
  
  end subroutine inout

  subroutine monitor_data_count(filename,n,recl_size)
    
    character(len=*),intent(in) :: filename
    integer, intent(out)        :: n,recl_size

    integer :: n1,n2,n3,nr,n_col,n_data,stat
    real(wp),allocatable :: data_test(:)

    recl_size=0
    do while(.true.)
      recl_size=recl_size+8

      open(12,file=filename,access='direct',recl=recl_size,action='read')
      read(12,rec=1,iostat=stat)n1
      read(12,rec=2,iostat=stat)n2
      read(12,rec=3,iostat=stat)n3
      close(12)
      print*,n1,n2,n3
      if(n2-n1 == n3-n2 .and. n2>n1) exit

    enddo
    n_data=(recl_size-1)/8
    print*,' **     recl_size:',recl_size
    print*,' **  data numbers:',n_data

    allocate(data_test(n_data))
    nr=1
    n_col=1
    open(12,file=filename,access='direct',recl=recl_size,action='read')
    do while(.true.)
      
      read(12,rec=nr,iostat=stat)n,data_test
      ! print*,n,data_test(1)
      if(stat==0) then
        nr = nr + 1
        n_col=n_col+1
      else
        exit
      endif
      !
    enddo
    close(12)
    print*,' ** last time step:',n,data_test(1)

  end subroutine monitor_data_count

  subroutine read_monitor_data(filename,recl_size,n_col,n_row,data_x,data_y)

    character(len=*),intent(in) :: filename
    integer, intent(in)         :: n_col,recl_size,n_row
    real(8), intent(out) :: data_x(n_col),data_y(n_col)

    integer :: num_first_file,num_last_file
    integer :: nr,n,n1,n2,n3,stat,n_data,i

    real(8),allocatable :: time(:),data(:,:)
    integer,allocatable :: nstep(:)

    n_data=(recl_size-1)/8-1
    allocate(nstep(n_col))
    allocate(time(n_col),data(n_data,n_col))

    print*,n_col,recl_size,n_data

    open(12,file=filename,access='direct',recl=recl_size,action='read')
    do i=1,n_col
      read(12,rec=i,iostat=stat)nstep(i),time(i),data(:,i)
      ! if(i<=2) print*,i,nstep(i),time(i)
    enddo
    close(12)
    print*,' >> ',filename

    data_x(:)=time(:)

    data_y(:)=data(n_row,:)

  end subroutine read_monitor_data

end module kernels