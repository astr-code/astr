program comparison
    implicit none
    character(len=200) :: line1, line2
    integer :: i, ios1, ios2, choice
    logical :: flag1=.true.

    print *, "Enter your choice: "
    read(*,*) choice

    if(choice == 1) then
        open(unit=10, file='gradcal_cpu.txt', status='old', action='read')
        open(unit=20, file='gradcal_gpu.txt', status='old', action='read')
        open(unit=30, file='gradcal.txt', status='unknown', action='write')
        open(unit=40, file='gradcal_ne.txt', status='unknown', action='write')
    else if(choice == 2) then
        open(unit=10, file='convection_cpu.txt', status='old', action='read')
        open(unit=20, file='convection_gpu.txt', status='old', action='read')
        open(unit=30, file='convection.txt', status='unknown', action='write')
        open(unit=40, file='convection_ne.txt', status='unknown', action='write')
    else if(choice == 3) then
        open(unit=10, file='diffusion_cpu.txt', status='old', action='read')
        open(unit=20, file='diffusion_gpu.txt', status='old', action='read')
        open(unit=30, file='diffusion.txt', status='unknown', action='write')
        open(unit=40, file='diffusion_ne.txt', status='unknown', action='write')
    else if(choice == 4) then
        open(unit=10, file='filterq_cpu.txt', status='old', action='read')
        open(unit=20, file='filterq_gpu.txt', status='old', action='read')
        open(unit=30, file='filterq.txt', status='unknown', action='write')
        open(unit=40, file='filterq_ne.txt', status='unknown', action='write')
    else if(choice == 5) then
        open(unit=10, file='q2fvar_cpu.txt', status='old', action='read')
        open(unit=20, file='q2fvar_gpu.txt', status='old', action='read')
        open(unit=30, file='q2fvar.txt', status='unknown', action='write')
        open(unit=40, file='q2fvar_ne.txt', status='unknown', action='write')
    else if(choice == 6) then
        open(unit=10, file='bchomo_cpu.txt', status='old', action='read')
        open(unit=20, file='bchomo_gpu.txt', status='old', action='read')
        open(unit=30, file='bchomo.txt', status='unknown', action='write')
        open(unit=40, file='bchomo_ne.txt', status='unknown', action='write')
    end if
    
    
    do
        read(10, '(A)', iostat=ios1) line1
        read(20, '(A)', iostat=ios2) line2

        if (ios1 /= 0 .and. ios2 /= 0) exit

        write(30, '(A)') trim(line1)
        write(30, '(A)') trim(line2)
        write(30, '(A)') '----------------------------------------------------------------------'

        if(line1/=line2) then
            flag1 = .false.
            write(40, '(A)') trim(line1)
            write(40, '(A)') trim(line2)
            write(40, '(A)') '----------------------------------------------------------------------'
        else
        end if

    end do

    if(flag1) then
        write(30, '(A)') 'Data is same'
        write(40, '(A)') '----- NO DATA -----'
    else
        write(30, '(A)') 'Data is not same'
    end if
    write(30, '(A)') '----------------------------------------------------------------------'

    close(10)
    close(20)
    close(30)
    close(40)
end program comparison
