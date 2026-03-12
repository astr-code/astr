program compact_filter_cpu
    use commvar
    use commtype
    use filter
    implicit none

    integer, parameter :: dim = 64
    type(compact_scheme) :: afilter
    real(8) :: f(-hm:dim+hm)
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
    afilter%last_node  = dim
    afilter%nbctype    = 1  ! 物理边界类型
    allocate(afilter%ac(3, dim+1))
    do i=0,dim
        afilter%ac(1,i) = 0.25d0
        afilter%ac(2,i) =  1.0d0
        afilter%ac(3,i) = 0.25d0
end do
    !-----------------------------
    ! 初始化测试函数
    !-----------------------------
    f = 0.d0
    do i = -hm, dim+hm
        f(i) = sin(2.d0*3.141592653589793d0*i/dim)
    end do

    !-----------------------------
    ! 调用 CPU 紧致滤波
    !-----------------------------
    f_filtered = compact_filter(afilter,f,dim)

    !-----------------------------
    ! 输出部分结果
    !-----------------------------
    print *, ' i    f_filtered_cpu'
    do i = 0, 10
        print '(I4, 1X, E25.16)', i, f_filtered(i)
    end do

end program compact_filter_cpu