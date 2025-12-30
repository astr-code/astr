module pastr_process

    use iso_fortran_env, only: wp => real64

    implicit none

    private

    public :: run_process_entry

contains

    !----------------------------------------------------------------------
    !> Central driver for post-processing
    !!  Calls modular processing routines: statistics, visualization, sampling
    !----------------------------------------------------------------------
    subroutine run_process_entry(command)

        use pastr_udf
        use pastr_blasius, only: blasius_solution

        character(len=*), intent(in) :: command

        select case (command)
           case ('viewxy')
              call field_view('xy')
           case ('viewxz')
              call field_view('xz')
           case ('view3d')
              call field_view('3d')
           case ('monitor')
              call monitor_data_process()
           case ('blasius')
              call blasius_solution()
           case ('stats')
              call stats_cal()
           case ('udf')
              call udf_func1()
           case default
              print*,' !! command not defined @ run_process_entry: ',command
              stop 1
        end select 

        ! !----------------------------------------------------------
        ! ! 1. Compute statistics
        ! !----------------------------------------------------------
        ! call compute_statistics(astrin, cfg)

        ! !----------------------------------------------------------
        ! ! 2. Sample flow data (e.g., probe points, slices)
        ! !----------------------------------------------------------
        ! call sample_flow_data(astrin, cfg)

        ! !----------------------------------------------------------
        ! ! 3. Generate visualization files
        ! !----------------------------------------------------------
        ! call generate_visualization(astrin, cfg)

        ! !----------------------------------------------------------
        ! ! Finished
        ! !----------------------------------------------------------
        ! write(*,*) "Post-processing completed successfully."

    end subroutine run_process_entry

    subroutine stats_cal

      use pastr_io,     only: parse_command_line
      use pastr_input,  only: read_astr_input
      use pastr_data_process, only: stats_read_process

      integer :: num_first_file,num_last_file

      call parse_command_line( inumber=num_first_file )
      call parse_command_line( inumber=num_last_file )

      call read_astr_input()

      if(num_first_file>0 .and. num_last_file>=num_first_file) then
      else
        call stats_read_process(mode='meanflow')
      endif

    end subroutine stats_cal

    subroutine field_view(mode)

        use pastr_io,     only: parse_command_line
        use pastr_input,  only: read_astr_input
        use pastr_field_view, only: write_xy_slice,write_xz_slice,write_3d_field

        character(len=*),intent(in) :: mode

        character(len=128) :: file_in,file_out
        character(len=4) :: visu_format,vform
        integer :: num_first_file,num_last_file,slice

        visu_format='plt'

        if(mode=='xy') then 
          write(*,*)' ** visulise flowfield from a xy slice'
        elseif(mode=='xz') then 
          write(*,*)' ** visulise flowfield from a xz slice'
        elseif(mode=='3d') then 
          write(*,*)' ** visulise 3d flowfield '
        endif

        write(*,*)' ** input either file numbers or names of files.'
        write(*,*)' ** examples: ./pastr viewxy                   0         2 110 -input datin/input.3d xdmf kslice'
        write(*,*)'              ./pastr viewxy outdat/flowfield.h5 tecxy.plt 110 -input datin/input.3d'

        call parse_command_line(  string=file_in )
        call parse_command_line(  string=file_out )
        call parse_command_line( inumber=num_first_file )
        call parse_command_line( inumber=num_last_file )
        call parse_command_line( inumber=slice )

        call read_astr_input()

        call refconst()

        call parse_command_line(  string=vform )

        if(len(trim(vform))>1) visu_format=vform
        if(trim(visu_format)=='plt') then
          write(*,*)' **          visu fomat: tecplot'
        elseif(trim(visu_format)=='xdmf') then
          write(*,*)' **          visu fomat: xdmf'
        else
          stop ' error 1 @ field_view'
        endif
        
        if(num_first_file>=0 .and. num_last_file>=0) then

          write(*,'(A,I0)')'  **      first file num: ',num_first_file
          write(*,'(A,I0)')'  **       last file num: ',num_last_file
          write(*,'(A,I0)')'  **      slice position: ',slice

          if(mode=='xy') then
            call write_xy_slice( nfirst=num_first_file, &
                                  nlast=num_last_file,  &
                                  slice=slice,format=trim(visu_format))
          elseif(mode=='xz') then
            call write_xz_slice( nfirst=num_first_file, &
                                  nlast=num_last_file,  &
                                  slice=slice,format=trim(visu_format))
          elseif(mode=='3d') then
          else
            print*,' !! mode not defined @ field_view ',mode
            stop 2
          endif

        elseif(len(trim(file_in))>3 .and. len(trim(file_in))>3) then

          write(*,*)' **    input field file: ',trim(file_in)
          write(*,*)' **    output visu file: ',trim(file_out)
          slice=num_first_file
          write(*,'(A,I0)')'  **         slice pos k: ',slice

          if(mode=='xy') then
            call write_xy_slice( filein=trim(file_in),  &
                                fileout=trim(file_out), &
                                  slice=slice,format=trim(visu_format))
          elseif(mode=='xz') then
            call write_xz_slice( filein=trim(file_in),  &
                                fileout=trim(file_out), &
                                  slice=slice,format=trim(visu_format))
          elseif(mode=='3d') then
            call write_3d_field( filein=trim(file_in),  &
                                fileout=trim(file_out), &
                                 format=trim(visu_format))
          else
            print*,' !! mode not defined @ field_view '
            stop 3
          endif

        else
          stop ' error 4 @ field_view'
        endif

    end subroutine field_view

    subroutine monitor_data_process

        use pastr_io,     only: parse_command_line,read_monitor_data,write_monitor_data
        use pastr_input,  only: read_astr_input,read_moniter_input
        use pastr_commvar,only: montype

        integer :: num_first_file,num_last_file,recl_size
        integer :: nr,n,n1,n2,n3,stat,n_data,n_col

        integer,allocatable :: imon(:),jmon(:),kmon(:)
        real(wp),allocatable :: data_test(:)
        logical :: ltemp
        character(len=16),allocatable :: data_name(:)
        type(montype),allocatable :: data_monitor(:)
        
        call read_astr_input()

        call read_moniter_input(imon,jmon,kmon)

        call parse_command_line( inumber=num_first_file )
        call parse_command_line( inumber=num_last_file )
        
        recl_size=0
        do while(.true.)
          recl_size=recl_size+8

          open(12,file='monitor/monitor0001.dat',access='direct',recl=recl_size,action='read')
          read(12,rec=1,iostat=stat)n1
          read(12,rec=2,iostat=stat)n2
          read(12,rec=3,iostat=stat)n3
          close(12)
          print*,n1,n2,n3
          if(n2-n1 == n3-n2 .and. n2>n1) exit

        enddo
        n_data=(recl_size-1)/8
        print*,' ** recl_size= ',recl_size,'data numbers:',n_data

        allocate(data_test(n_data))
        nr=1
        n_col=1
        open(12,file='monitor/monitor0001.dat',access='direct',recl=recl_size,action='read')
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

        n_col=n_col-1

        call read_monitor_data(num_first_file,num_last_file,n_data,n_col,data_monitor)

        call write_monitor_data(data_monitor)

    end subroutine monitor_data_process

    subroutine refconst
      !
      use pastr_commvar
      use pastr_thermo_phys, only : sos
      !
      prandtl=0.72d0
      gamma=1.4d0
      rgas=287.1d0
      !
      ref_miu=1.458d-6*ref_t**1.5d0/(ref_t+110.4d0)
      !
      tempconst=110.4d0/ref_t
      tempconst1=1.d0+tempconst
      !
      !
      const1=gamma*(gamma-1.d0)*mach**2
      const2=gamma*mach**2
      const3=(gamma-1.d0)/3.d0*prandtl*(mach**2)
      const4=(gamma-1.d0)*mach**2*reynolds*prandtl
      const5=(gamma-1.d0)*mach**2
      const6=gamma-1.d0
      const7=(gamma-1.d0)*mach**2*reynolds*prandtl
      const8=1.d0/reynolds
      
      if(nondimen) then
        uinf=1.d0
        vinf=0.d0
        tinf=1.d0
        roinf=1.d0
        
        pinf=roinf*tinf/const2
      else
        uinf =ref_vel
        vinf =0.d0
        tinf =ref_t
        roinf=ref_den
  
        pinf=roinf*tinf*rgas
  
        reynolds=ref_den*ref_vel*ref_len/ref_miu
        mach    =ref_vel/sos(ref_t)
      endif
      !
      print*,' ** Reference parameters calculation done.'
      !
    end subroutine refconst

    !----------------------------------------------------------------------
    ! Example subroutine: compute statistics
    !----------------------------------------------------------------------
    subroutine compute_statistics()

        ! Placeholder: compute mean, RMS, max/min of flow variables
        write(*,*) "Computing statistics..."
        ! TODO: implement actual statistics computations
    end subroutine compute_statistics

    !----------------------------------------------------------------------
    ! Example subroutine: sample flow data
    !----------------------------------------------------------------------
    subroutine sample_flow_data()

        write(*,*) "Sampling flow data..."
        ! TODO: implement point probes, slices, or line sampling
    end subroutine sample_flow_data

    !----------------------------------------------------------------------
    ! Example subroutine: generate visualization
    !----------------------------------------------------------------------
    subroutine generate_visualization()

        write(*,*) "Generating visualization files..."
        ! TODO: write VTK, HDF5, or other format files
    end subroutine generate_visualization

end module pastr_process