module udf_postprocess
  !
  implicit none
  !
  contains
  !
  subroutine flame2d_pp(inputfile)
    !
    use cmdefne
    use hdf5io
    use tecio
    use WriteVTK
    use commvar, only: num_species
    !
#ifdef COMB
    use cantera
    use thermchem,only: chemrep,chemread,thermdyn,heatrate
    !
    type(phase_t) :: mixture
    !
#endif
    !
    character(len=*),intent(in) :: inputfile
    !
    ! local data
    integer :: nfirst,nlast,n,im,jm,km,jsp,i,j
    logical :: lihomo,ljhomo,lkhomo
    real(8) :: ref_tem,reynolds,mach,qdotmax
    character(len=32) :: gridfile,chemfile
    !
    real(8),allocatable,dimension(:,:) :: x_2d,y_2d,z_2d,ro_2d,u_2d,   &
                                          v_2d,w_2d,p_2d,t_2d
    real(8),allocatable,dimension(:,:,:) :: spec
    real(8),allocatable,dimension(:,:) :: vtmp,qdot
    !
    character(len=4) :: stepname
    character(len=5) :: cnmuber1,cnmuber2
    character(len=3) :: spname
    character(len=10),allocatable :: spcsym(:)
    character(len=10) :: phase_id
    character(len=128) :: filename,tecname
    !
    call readkeyboad(cnmuber1)
    call readkeyboad(cnmuber2)
    !
    read(cnmuber1,*)nfirst
    read(cnmuber2,*)nlast
    !
    print*,' ==========================readinput=========================='
    !
    open(11,file=inputfile,form='formatted',status='old')
    read(11,'(///////)')
    read(11,*)im,jm,km
    read(11,"(/)")
    read(11,*)lihomo,ljhomo,lkhomo
    read(11,'(//////////)')
    read(11,*)ref_tem,reynolds,mach
    read(11,'(///////)')
    read(11,*)num_species
    read(11,'(//////////////////)')
    read(11,'(A)')gridfile
#ifdef COMB
    read(11,"(/)")
    read(11,'(A)')chemfile
#endif
    close(11)
    print*,' >> ',inputfile
    !
    print*,' ** grid file: ',trim(gridfile)
    !
    allocate(x_2d(0:im,0:jm),y_2d(0:im,0:jm))
    call H5ReadSubset(x_2d,im,jm,km,'x',gridfile,kslice=0)
    call H5ReadSubset(y_2d,im,jm,km,'y',gridfile,kslice=0)
    !
#ifdef COMB
    call chemread(trim(chemfile))
    call thermdyn
    print*,' ** num_species: ',num_species
    !
    phase_id='gas'
    !---CANTERA---
    mixture=importPhase(trim(chemfile),trim(phase_id))
    !
    allocate(spcsym(num_species))
    !
    write(*,'(2X,A12)')'------------'
    write(*,'(2X,A12)')'   chemistry'
    write(*,'(2X,A12)')'----+-------'
    do jsp=1,num_species
      
      call getSpeciesName(mixture,jsp,spcsym(jsp))
      
      write(*,'(2X,I3,A3,A5)')jsp,' | ',spcsym(jsp)
    enddo
    write(*,'(2X,A12)')'----+-------'
    !
    !
#endif
    !
    allocate(ro_2d(0:im,0:jm),u_2d(0:im,0:jm), v_2d(0:im,0:jm),      &
              w_2d(0:im,0:jm),p_2d(0:im,0:jm), t_2d(0:im,0:jm)       )
    !
    if(num_species>0) then

      allocate(spec(0:im,0:jm,1:num_species) )
      allocate(vtmp(0:im,0:jm),qdot(0:im,0:jm) )
    endif
    !
    do n=nfirst,nlast
      !
      write(stepname,'(i4.4)')n
      filename='outdat/flowfield'//stepname//'.h5'
      !
      !
      call h5_read2dfrom3d(ro_2d,im,jm,km,'ro',trim(filename),kslice=0)
      call h5_read2dfrom3d( u_2d,im,jm,km,'u1',trim(filename),kslice=0)
      call h5_read2dfrom3d( v_2d,im,jm,km,'u2',trim(filename),kslice=0)
      call h5_read2dfrom3d( p_2d,im,jm,km, 'p',trim(filename),kslice=0)
      call h5_read2dfrom3d( t_2d,im,jm,km, 't',trim(filename),kslice=0)
      !
      tecname='Results/tecflow'//stepname//'.plt'
      !
      do jsp=1,num_species
        write(spname,'(i3.3)') jsp
        call h5_read2dfrom3d(vtmp,im,jm,km, 'sp'//spname,trim(filename),kslice=0)
        spec(:,:,jsp)=vtmp
      enddo
      !
#ifdef COMB
      qdotmax=0.d0
      do j=0,jm
      do i=0,im
        qdot(i,j)=heatrate(ro_2d(i,j),t_2d(i,j),spec(i,j,:))
        !
        if(qdot(i,j)>qdotmax)qdotmax=qdot(i,j)
      enddo
      enddo
      !
      print*,' ** qdotmax=',qdotmax
#endif
      !
      call tecbin(trim(tecname),x_2d,'x',y_2d,'y',ro_2d,'ro',  &
      	                                           u_2d,'u',   &
      	                                           v_2d,'v',   &
      	                                           p_2d,'p',   &
      	                                           t_2d,'t',   &
      	                                        spec(:,:,1),spcsym(1), &
      	                                        spec(:,:,2),spcsym(2), &
      	                                        spec(:,:,5),spcsym(5), &
      	                                        qdot,'HRR')
     enddo

    !   ! call writeprvbin(outputfile,x_2d,'x',y_2d,'y',ro_2d,'ro',u_2d,'u',   &
    !   !                                v_2d,'v',p_2d,'p',t_2d,'t',im,jm)
    !    call tecbin(outputfile,x_2d,'x',y_2d,'y',ro_2d,'ro',u_2d,'u',   &
    !                                      v_2d,'v',p_2d,'p',t_2d,'t')
    ! elseif(viewmode=='3d') then
    !   allocate(x(0:im,0:jm,0:km),y(0:im,0:jm,0:km),z(0:im,0:jm,0:km))
    ! elseif(viewmode=='1d') then
    !   !
    !   allocate(x_1d(0:im),ro_1d(0:im),u_1d(0:im),p_1d(0:im),t_1d(0:im))
    !   !
    !   call H5ReadSubset( x_1d,im,jm,km, 'x',gridfile,jslice=jm/2-1,kslice=-1)
    !   !
    !   call H5ReadSubset(ro_1d,im,jm,km,'ro',flowfile,jslice=jm/2-1,kslice=-1)
    !   call H5ReadSubset( u_1d,im,jm,km,'u1',flowfile,jslice=jm/2-1,kslice=-1)
    !   call H5ReadSubset( p_1d,im,jm,km, 'p',flowfile,jslice=jm/2-1,kslice=-1)
    !   call H5ReadSubset( t_1d,im,jm,km, 't',flowfile,jslice=jm/2-1,kslice=-1)
    !   !
    !   allocate(spc_1d(0:im,1:num_species))
    !   allocate(var1d(0:im))
    !   do jsp=1,num_species
    !     write(spname,'(i3.3)')jsp
    !     call H5ReadSubset(var1d,im,jm,km,'sp'//spname,flowfile,jslice=jm/2-1,kslice=-1)
    !     spc_1d(:,jsp)=var1d
    !   enddo
    !   !
    !   open(18,file=outputfile)
    !   write(18,"(8(1X,A15))")'x','ro','u','p','t','Y1','Y2','Y3'
    !   write(18,"(8(1X,E15.7E3))")(x_1d(i),ro_1d(i),u_1d(i),p_1d(i), &
    !               t_1d(i),spc_1d(i,1),spc_1d(i,2),spc_1d(i,3),i=0,im)
    !   close(18)
    !   print*,' << ',outputfile
    !   !
    !   call h5srite(var=ro_1d,varname='ro',filename='flowini1d.h5',explicit=.true.,newfile=.true.)
    !   call h5srite(var= u_1d,varname='u1',filename='flowini1d.h5',explicit=.true.,newfile=.false.)
    !   call h5srite(var=t_1d, varname= 't',filename='flowini1d.h5',explicit=.true.,newfile=.false.)
    !   call h5srite(var=p_1d, varname= 'p',filename='flowini1d.h5',explicit=.true.,newfile=.false.)
    !   !
    !   do jsp=1,num_species
    !     write(spname,'(i3.3)')jsp
    !     call h5srite(var=spc_1d(:,jsp),varname='sp'//spname,filename='flowini1d.h5',explicit=.true.,newfile=.false.)
    !   enddo
    ! else
    !   print*,viewmode
    !   stop ' !! mode is not defined @ fieldview'
    ! endif
    !
  end subroutine flame2d_pp
  !
end module udf_postprocess