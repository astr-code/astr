program test_compact_filter_cpu
    use commvar
    use commtype
    use filter
    implicit none

    integer, parameter :: n = 64
    type(compact_scheme) :: afilter
    real(8) :: f(0:n)
    real(8), allocatable :: f_filtered(:)
    integer :: i
    real(8) :: alfa, beter_halo, beter_bouond

    !-----------------------------
    ! 初始化滤波器系数
    !-----------------------------
    alfa = 0.49d0
    beter_halo = 1.11d0
    beter_bouond = 1.09d0
    call filter_coefficient_cal(alfa, beter_halo, beter_bouond)

    !-----------------------------
    ! 初始化紧致滤波器对象
    !-----------------------------
    afilter%first_node = 0
    afilter%last_node  = n
    afilter%nbctype    = 1  ! 物理边界类型

    !-----------------------------
    ! 初始化测试函数
    !-----------------------------
    do i = 0, n
        f(i) = sin(2.d0*3.141592653589793d0*i/n)
    end do

    !-----------------------------
    ! 调用 CPU 紧致滤波
    !-----------------------------
    f_filtered = compact_filter_rhs(afilter, f, n)

    !-----------------------------
    ! 输出部分结果
    !-----------------------------
    print *, ' i    f_filtered_cpu'
    do i = 0, 10
        print '(I4,E20.10)', i, f_filtered(i)
    end do

end program test_compact_filter_cpu