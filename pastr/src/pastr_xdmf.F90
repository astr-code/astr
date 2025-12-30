module pastr_xdmf

  use iso_fortran_env, only: wp => real64

  implicit none
  !
  Interface xdmfwriter
    module procedure xdmfwriter_2d
    module procedure xdmfwriter_2drec
    module procedure xdmfwriter_2drec_block
    module procedure xdmfwriter_2drec_list_xy
    module procedure xdmfwriter_3d
    module procedure xdmfwriter_3drec
    module procedure xdmfwriter_3drec_block
    module procedure xdmfwriter_3drec_list
    module procedure xdmfwriter_3dbox
  end Interface xdmfwriter

contains

  !+-------------------------------------------------------------------+
  !| This subroutine is used to write a xdmf head file for             |
  !| visulisation of flow field.                                       |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 01-07-2022  | Created by J. Fang @ Warrington                     |
  !+-------------------------------------------------------------------+
  subroutine xdmfwriter_3d(dir,filename,gridfile,var1,var1name,var2,var2name, &
                                              var3,var3name,var4,var4name, &
                                              var5,var5name,var6,var6name, im,jm,km)
    
    use pastr_h5io
    use pastr_utility, only : make_dir

    ! arguments
    character(len=*),intent(in) :: dir,filename,gridfile
    real(wp),intent(in),dimension(:,:,:) :: var1
    character(len=*),intent(in) :: var1name
    integer :: im,jm,km
    real(wp),intent(in),dimension(:,:,:),optional :: var2,var3,var4,var5,var6
    character(len=*),intent(in),optional :: var2name,var3name,var4name,var5name,var6name
    !
    real(4),allocatable :: bufr4(:,:,:)
    !
    ! local data
    integer :: fh,i
    !
    allocate(bufr4(0:im,0:jm,0:km))
    
    call make_dir(dir)

    fh=18
    !
    ! write the head of xdmf and grid
    open(fh,file=dir//filename//'.xdmf',form='formatted')
    write(fh,'(A)')'<?xml version="1.0" ?>'
    write(fh,'(A)')'<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
    write(fh,'(A)')'<Xdmf xmlns:xi="http://www.w3.org/2001/XInclude" Version="2.0">'
    write(fh,'(A)')'  <Domain>'
    !
    write(fh,'(A,3(1X,I0),A)')'    <Topology name="topo" TopologyType="3DSMESH" Dimensions="',  &
                              km+1,jm+1,im+1,'"> </Topology>'
    write(fh,'(A)')'    <Geometry name="geo" Type="X_Y_Z">'
    write(fh,'(A)')'      <DataItem Format="HDF" DataType="Float" Precision="4" Endian="little" Seek="0"'
    write(fh,'(A,3(1X,I0),3(A))')'                Dimensions="',km+1,jm+1,im+1,'"> ',gridfile,':x </DataItem>'
    write(fh,'(A)')'      <DataItem Format="HDF" DataType="Float" Precision="4" Endian="little" Seek="0"'
    write(fh,'(A,3(1X,I0),3(A))')'                Dimensions="',km+1,jm+1,im+1,'"> ',gridfile,':y </DataItem>'
    write(fh,'(A)')'      <DataItem Format="HDF" DataType="Float" Precision="4" Endian="little" Seek="0"'
    write(fh,'(A,3(1X,I0),3(A))')'                Dimensions="',km+1,jm+1,im+1,'"> ',gridfile,':z </DataItem>'
    write(fh,'(A)')'    </Geometry>'
    !
    write(fh,'(A)')'    <Grid Name="001" GridType="Uniform">'
    write(fh,'(A)')'      <Topology Reference="/Xdmf/Domain/Topology[1]"/>'
    write(fh,'(A)')'      <Geometry Reference="/Xdmf/Domain/Geometry[1]"/>'
    !
    write(fh,'(3(A))')'      <Attribute Name="',var1name,'" Center="Node">'
    write(fh,'(A)')'        <DataItem Format="HDF" DataType="Float" Precision="4" Endian="little" Seek="0"'
    write(fh,'(A,3(1X,I0),5(A))')'                   Dimensions="',km+1,jm+1,im+1,'"> ',filename//'.h5',':',var1name,'</DataItem>'
    write(fh,'(A)')'      </Attribute>'
    bufr4=var1
    call H5WriteArray(bufr4,im,jm,km,var1name,dir//filename//'.h5')

    if(present(var2)) then
      write(fh,'(3(A))')'      <Attribute Name="',var2name,'" Center="Node">'
      write(fh,'(A)')'        <DataItem Format="HDF" DataType="Float" Precision="4" Endian="little" Seek="0"'
      write(fh,'(A,3(1X,I0),5(A))')'                   Dimensions="',km+1,jm+1,im+1,'"> ',filename//'.h5',':',var2name,'</DataItem>'
      write(fh,'(A)')'      </Attribute>'
      bufr4=var2
      call H5WriteArray(bufr4,im,jm,km,var2name,dir//filename//'.h5')
    endif
    if(present(var3)) then
      write(fh,'(3(A))')'      <Attribute Name="',var3name,'" Center="Node">'
      write(fh,'(A)')'        <DataItem Format="HDF" DataType="Float" Precision="4" Endian="little" Seek="0"'
      write(fh,'(A,3(1X,I0),5(A))')'                   Dimensions="',km+1,jm+1,im+1,'"> ',filename//'.h5',':',var3name,'</DataItem>'
      write(fh,'(A)')'      </Attribute>'
      bufr4=var3
      call H5WriteArray(bufr4,im,jm,km,var3name,dir//filename//'.h5')
    endif
    if(present(var4)) then
      write(fh,'(3(A))')'      <Attribute Name="',var4name,'" Center="Node">'
      write(fh,'(A)')'        <DataItem Format="HDF" DataType="Float" Precision="4" Endian="little" Seek="0"'
      write(fh,'(A,3(1X,I0),5(A))')'                   Dimensions="',km+1,jm+1,im+1,'"> ',filename//'.h5',':',var4name,'</DataItem>'
      write(fh,'(A)')'      </Attribute>'
      bufr4=var4
      call H5WriteArray(bufr4,im,jm,km,var4name,dir//filename//'.h5')
    endif
    if(present(var5)) then
      write(fh,'(3(A))')'      <Attribute Name="',var5name,'" Center="Node">'
      write(fh,'(A)')'        <DataItem Format="HDF" DataType="Float" Precision="4" Endian="little" Seek="0"'
      write(fh,'(A,3(1X,I0),5(A))')'                   Dimensions="',km+1,jm+1,im+1,'"> ',filename//'.h5',':',var5name,'</DataItem>'
      write(fh,'(A)')'      </Attribute>'
      bufr4=var5
      call H5WriteArray(bufr4,im,jm,km,var5name,dir//filename//'.h5')
    endif
    if(present(var6)) then
      write(fh,'(3(A))')'      <Attribute Name="',var6name,'" Center="Node">'
      write(fh,'(A)')'        <DataItem Format="HDF" DataType="Float" Precision="4" Endian="little" Seek="0"'
      write(fh,'(A,3(1X,I0),5(A))')'                   Dimensions="',km+1,jm+1,im+1,'"> ',filename//'.h5',':',var6name,'</DataItem>'
      write(fh,'(A)')'      </Attribute>'
      bufr4=var6
      call H5WriteArray(bufr4,im,jm,km,var6name,dir//filename//'.h5')
    endif

    write(fh,'(A)')'     </Grid>'
    !
    write(fh,'(A)')'  </Domain>'
    write(fh,'(A)')'</Xdmf>'
    !
    close(fh)
    !
    !
  end subroutine xdmfwriter_3d
  !
  subroutine xdmfwriter_3dbox(dir,filename,deltax,var1,var1name,var2,var2name, &
                                                  var3,var3name,var4,var4name, &
                                                  var5,var5name,var6,var6name, im,jm,km)
    !
    ! arguments
    character(len=*),intent(in) :: dir,filename
    real(wp),intent(in) :: deltax
    real(wp),intent(in),dimension(:,:,:) :: var1
    character(len=*),intent(in) :: var1name
    integer :: im,jm,km
    real(wp),intent(in),dimension(:,:,:),optional :: var2,var3,var4,var5,var6
    character(len=*),intent(in),optional :: var2name,var3name,var4name,var5name,var6name
    !
    real(4),allocatable :: bufr4(:,:,:)
    !
    ! local data
    integer :: fh,i
    !
    allocate(bufr4(0:im,0:jm,0:km))
    !
    fh=18
    !
    ! write the head of xdmf and grid
    open(fh,file=dir//filename//'.xdmf',form='formatted')
    write(fh,'(A)')'<?xml version="1.0" ?>'
    write(fh,'(A)')'<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
    write(fh,'(A)')'<Xdmf xmlns:xi="http://www.w3.org/2001/XInclude" Version="2.0">'
    write(fh,'(A)')'  <Domain>'
    !
    write(fh,'(A)')'    <Topology name="topo" TopologyType="3DCoRectMesh"'
    write(fh,'(A,3(1X,I0),A)')'        Dimensions="',km+1,jm+1,im+1,'">'
    write(fh,'(A)')'    </Topology>'
    write(fh,'(A)')'     <Geometry name="geo" Type="ORIGIN_DXDYDZ">'
    write(fh,'(A)')'         <!-- Origin -->'
    write(fh,'(A)')'         <DataItem Format="XML" Dimensions="3">'
    write(fh,'(A)')'         0.0 0.0 0.0'
    write(fh,'(A)')'         </DataItem>'
    write(fh,'(A)')'         <!-- DxDyDz -->'
    write(fh,'(A)')'         <DataItem Format="XML" Dimensions="3">'
    write(fh,'(A,3(1X,E15.7E3))')'         ',deltax,deltax,deltax
    write(fh,'(A)')'         </DataItem>'
    write(fh,'(A)')'     </Geometry>'
    !
    write(fh,'(A)')'    <Grid Name="001" GridType="Uniform">'
    write(fh,'(A)')'      <Topology Reference="/Xdmf/Domain/Topology[1]"/>'
    write(fh,'(A)')'      <Geometry Reference="/Xdmf/Domain/Geometry[1]"/>'
    !
    write(fh,'(3(A))')'      <Attribute Name="',var1name,'" Center="Node">'
    write(fh,'(A)')'        <DataItem Format="Binary" DataType="Float" Precision="4" Endian="little" Seek="0"'
    write(fh,'(A,3(1X,I0),3(A))')'                   Dimensions="',km+1,jm+1,im+1,'"> ',filename//'-'//var1name,'</DataItem>'
    write(fh,'(A)')'      </Attribute>'

    bufr4=var1
    open(fh+1,file=dir//filename//'-'//var1name,access="stream")
    write(fh+1)bufr4
    close(fh+1)
    print*,' << ',dir//filename//'-'//var1name

    if(present(var2)) then
      write(fh,'(3(A))')'      <Attribute Name="',var2name,'" Center="Node">'
      write(fh,'(A)')'        <DataItem Format="Binary" DataType="Float" Precision="4" Endian="little" Seek="0"'
      write(fh,'(A,3(1X,I0),3(A))')'                   Dimensions="',km+1,jm+1,im+1,'"> ',filename//'-'//var2name,'</DataItem>'
      write(fh,'(A)')'      </Attribute>'
  
      bufr4=var2
      open(fh+1,file=dir//filename//'-'//var2name,access="stream")
      write(fh+1)bufr4
      close(fh+1)
      print*,' << ',dir//filename//'-'//var2name
    endif
    if(present(var3)) then
      write(fh,'(3(A))')'      <Attribute Name="',var3name,'" Center="Node">'
      write(fh,'(A)')'        <DataItem Format="Binary" DataType="Float" Precision="4" Endian="little" Seek="0"'
      write(fh,'(A,3(1X,I0),3(A))')'                   Dimensions="',km+1,jm+1,im+1,'"> ',filename//'-'//var3name,'</DataItem>'
      write(fh,'(A)')'      </Attribute>'
  
      bufr4=var3
      open(fh+1,file=dir//filename//'-'//var3name,access="stream")
      write(fh+1)bufr4
      close(fh+1)
      print*,' << ',dir//filename//'-'//var3name
    endif
    if(present(var4)) then
      write(fh,'(3(A))')'      <Attribute Name="',var4name,'" Center="Node">'
      write(fh,'(A)')'        <DataItem Format="Binary" DataType="Float" Precision="4" Endian="little" Seek="0"'
      write(fh,'(A,3(1X,I0),3(A))')'                   Dimensions="',km+1,jm+1,im+1,'"> ',filename//'-'//var4name,'</DataItem>'
      write(fh,'(A)')'      </Attribute>'
  
      bufr4=var4
      open(fh+1,file=dir//filename//'-'//var4name,access="stream")
      write(fh+1)bufr4
      close(fh+1)
      print*,' << ',dir//filename//'-'//var4name
    endif
    if(present(var5)) then
      write(fh,'(3(A))')'      <Attribute Name="',var5name,'" Center="Node">'
      write(fh,'(A)')'        <DataItem Format="Binary" DataType="Float" Precision="4" Endian="little" Seek="0"'
      write(fh,'(A,3(1X,I0),3(A))')'                   Dimensions="',km+1,jm+1,im+1,'"> ',filename//'-'//var5name,'</DataItem>'
      write(fh,'(A)')'      </Attribute>'
  
      bufr4=var5
      open(fh+1,file=dir//filename//'-'//var5name,access="stream")
      write(fh+1)bufr4
      close(fh+1)
      print*,' << ',dir//filename//'-'//var5name
    endif
    if(present(var6)) then
      write(fh,'(3(A))')'      <Attribute Name="',var6name,'" Center="Node">'
      write(fh,'(A)')'        <DataItem Format="Binary" DataType="Float" Precision="4" Endian="little" Seek="0"'
      write(fh,'(A,3(1X,I0),3(A))')'                   Dimensions="',km+1,jm+1,im+1,'"> ',filename//'-'//var6name,'</DataItem>'
      write(fh,'(A)')'      </Attribute>'
  
      bufr4=var6
      open(fh+1,file=dir//filename//'-'//var6name,access="stream")
      write(fh+1)bufr4
      close(fh+1)
      print*,' << ',dir//filename//'-'//var6name
    endif

    write(fh,'(A)')'     </Grid>'
    !
    write(fh,'(A)')'  </Domain>'
    write(fh,'(A)')'</Xdmf>'
    !
    close(fh)
    !
    !
  end subroutine xdmfwriter_3dbox
  
  subroutine xdmfwriter_3drec(dir,filename,x,y,z,var1,var1name,var2,var2name, &
                                                 var3,var3name,var4,var4name, &
                                                 var5,var5name,var6,var6name)
    !
    ! arguments
    character(len=*),intent(in) :: dir,filename
    real(wp),intent(in),dimension(:) :: x,y,z
    real(wp),intent(in),dimension(:,:,:) :: var1
    character(len=*),intent(in) :: var1name

    real(wp),intent(in),dimension(:,:,:),optional :: var2,var3,var4,var5,var6
    character(len=*),intent(in),optional :: var2name,var3name,var4name,var5name,var6name
    !
    real(4),allocatable :: bufr4(:,:,:)
    !
    ! local data
    integer :: fh,i
    integer :: im,jm,km
    logical :: lfex
    !
    im=size(x)
    jm=size(y)
    km=size(z)

    inquire(file=dir//filename//'-grid-x.bin',exist=lfex)
    if(.not. lfex) then
      open(17,file=dir//filename//'-grid-x.bin',access="stream")
      write(17)sngl(x)
      close(17)
      print*,' << ',dir//filename//'-grid-x.bin'
    endif
    inquire(file=dir//filename//'-grid-y.bin',exist=lfex)
    if(.not. lfex) then
      open(17,file=dir//filename//'-grid-y.bin',access="stream")
      write(17)sngl(y)
      close(17)
      print*,' << ',dir//filename//'-grid-y.bin'
    endif
    inquire(file=dir//filename//'-grid-z.bin',exist=lfex)
    if(.not. lfex) then
      open(17,file=dir//filename//'-grid-z.bin',access="stream")
      write(17)sngl(z)
      close(17)
      print*,' << ',dir//filename//'-grid-z.bin'
    endif

    !
    allocate(bufr4(im,jm,km))
    !
    fh=18
    !
    ! write the head of xdmf and grid
    open(fh,file=dir//filename//'.xdmf',form='formatted')
    write(fh,'(A)')'<?xml version="1.0" ?>'
    write(fh,'(A)')'<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
    write(fh,'(A)')'<Xdmf xmlns:xi="http://www.w3.org/2001/XInclude" Version="2.0">'
    write(fh,'(A)')'  <Domain>'
    !
    write(fh,'(A,3(1X,I0),A)')'    <Topology name="topo" TopologyType="3DRectMesh" Dimensions="',  &
                              km,jm,im,'"> </Topology>'
    write(fh,'(A)')'    <Geometry name="geo" Type="VXVYVZ">'
    write(fh,'(A)')'      <DataItem Format="Binary" DataType="Float" Precision="4" Endian="little" Seek="0" '
    write(fh,'(A,1X,I0,3A)') '                Dimensions=" ',im,'"> '//filename//'-grid-x.bin  </DataItem>'
    write(fh,'(A)')'      <DataItem Format="Binary" DataType="Float" Precision="4" Endian="little" Seek="0" '
    write(fh,'(A,1X,I0,3A)') '                Dimensions=" ',jm,'"> '//filename//'-grid-y.bin  </DataItem>'
    write(fh,'(A)')'      <DataItem Format="Binary" DataType="Float" Precision="4" Endian="little" Seek="0" '
    write(fh,'(A,1X,I0,3A)') '                Dimensions=" ',km,'"> '//filename//'-grid-z.bin  </DataItem>'
    write(fh,'(A)')'    </Geometry>'
    write(fh,*)
    !
    write(fh,'(A)')'    <Grid Name="001" GridType="Uniform">'
    write(fh,'(A)')'      <Topology Reference="/Xdmf/Domain/Topology[1]"/>'
    write(fh,'(A)')'      <Geometry Reference="/Xdmf/Domain/Geometry[1]"/>'
    !
    write(fh,'(3(A))')'      <Attribute Name="',var1name,'" Center="Node">'
    write(fh,'(A)')'        <DataItem Format="Binary" DataType="Float" Precision="4" Endian="little" Seek="0"'
    write(fh,'(A,3(1X,I0),3(A))')'                   Dimensions="',km,jm,im,'"> ',filename//'-'//var1name,'</DataItem>'
    write(fh,'(A)')'      </Attribute>'

    bufr4=var1
    open(fh+1,file=dir//filename//'-'//var1name,access="stream")
    write(fh+1)bufr4
    close(fh+1)
    print*,' << ',dir//filename//'-'//var1name

    if(present(var2)) then
      write(fh,'(3(A))')'      <Attribute Name="',var2name,'" Center="Node">'
      write(fh,'(A)')'        <DataItem Format="Binary" DataType="Float" Precision="4" Endian="little" Seek="0"'
      write(fh,'(A,3(1X,I0),3(A))')'                   Dimensions="',km,jm,im,'"> ',filename//'-'//var2name,'</DataItem>'
      write(fh,'(A)')'      </Attribute>'
  
      bufr4=var2
      open(fh+1,file=dir//filename//'-'//var2name,access="stream")
      write(fh+1)bufr4
      close(fh+1)
      print*,' << ',dir//filename//'-'//var2name
    endif
    if(present(var3)) then
      write(fh,'(3(A))')'      <Attribute Name="',var3name,'" Center="Node">'
      write(fh,'(A)')'        <DataItem Format="Binary" DataType="Float" Precision="4" Endian="little" Seek="0"'
      write(fh,'(A,3(1X,I0),3(A))')'                   Dimensions="',km,jm,im,'"> ',filename//'-'//var3name,'</DataItem>'
      write(fh,'(A)')'      </Attribute>'
  
      bufr4=var3
      open(fh+1,file=dir//filename//'-'//var3name,access="stream")
      write(fh+1)bufr4
      close(fh+1)
      print*,' << ',dir//filename//'-'//var3name
    endif
    if(present(var4)) then
      write(fh,'(3(A))')'      <Attribute Name="',var4name,'" Center="Node">'
      write(fh,'(A)')'        <DataItem Format="Binary" DataType="Float" Precision="4" Endian="little" Seek="0"'
      write(fh,'(A,3(1X,I0),3(A))')'                   Dimensions="',km,jm,im,'"> ',filename//'-'//var4name,'</DataItem>'
      write(fh,'(A)')'      </Attribute>'
  
      bufr4=var4
      open(fh+1,file=dir//filename//'-'//var4name,access="stream")
      write(fh+1)bufr4
      close(fh+1)
      print*,' << ',dir//filename//'-'//var4name
    endif
    if(present(var5)) then
      write(fh,'(3(A))')'      <Attribute Name="',var5name,'" Center="Node">'
      write(fh,'(A)')'        <DataItem Format="Binary" DataType="Float" Precision="4" Endian="little" Seek="0"'
      write(fh,'(A,3(1X,I0),3(A))')'                   Dimensions="',km,jm,im,'"> ',filename//'-'//var5name,'</DataItem>'
      write(fh,'(A)')'      </Attribute>'
  
      bufr4=var5
      open(fh+1,file=dir//filename//'-'//var5name,access="stream")
      write(fh+1)bufr4
      close(fh+1)
      print*,' << ',dir//filename//'-'//var5name
    endif
    if(present(var6)) then
      write(fh,'(3(A))')'      <Attribute Name="',var6name,'" Center="Node">'
      write(fh,'(A)')'        <DataItem Format="Binary" DataType="Float" Precision="4" Endian="little" Seek="0"'
      write(fh,'(A,3(1X,I0),3(A))')'                   Dimensions="',km,jm,im,'"> ',filename//'-'//var6name,'</DataItem>'
      write(fh,'(A)')'      </Attribute>'
  
      bufr4=var6
      open(fh+1,file=dir//filename//'-'//var6name,access="stream")
      write(fh+1)bufr4
      close(fh+1)
      print*,' << ',dir//filename//'-'//var6name
    endif

    write(fh,'(A)')'     </Grid>'
    !
    write(fh,'(A)')'  </Domain>'
    write(fh,'(A)')'</Xdmf>'
    !
    close(fh)
    !
    !
  end subroutine xdmfwriter_3drec


  subroutine xdmfwriter_3drec_list(dir,filename,x,y,z,var,varname)
    
    use pastr_utility, only : make_dir
    
    ! arguments
    character(len=*),intent(in) :: dir,filename
    real(wp),intent(in),dimension(:) :: x,y,z
    real(wp),intent(in),dimension(:,:,:,:) :: var
    character(len=*),intent(in) :: varname(:)
    !
    real(4),allocatable :: bufr4(:,:,:)
    !
    ! local data
    integer :: fh,i,n,nvar
    integer :: im,jm,km
    logical :: lfex
    
    call make_dir(dir)

    im=size(var,1)
    jm=size(var,2)
    km=size(var,3)
    nvar=size(var,4)

    inquire(file=dir//'/grid-x.bin',exist=lfex)
    if(.not. lfex) then
      open(17,file=dir//'/grid-x.bin',access="stream")
      write(17)sngl(x)
      close(17)
      print*,' << ',dir//'/grid-x.bin'
    endif
    inquire(file=dir//'/grid-y.bin',exist=lfex)
    if(.not. lfex) then
      open(17,file=dir//'/grid-y.bin',access="stream")
      write(17)sngl(y)
      close(17)
      print*,' << ',dir//'/grid-y.bin'
    endif
    inquire(file=dir//'/grid-z.bin',exist=lfex)
    if(.not. lfex) then
      open(17,file=dir//'/grid-z.bin',access="stream")
      write(17)sngl(z)
      close(17)
      print*,' << ',dir//'/grid-z.bin'
    endif

    !
    allocate(bufr4(im,jm,km))
    !
    fh=18
    !
    ! write the head of xdmf and grid
    open(fh,file=dir//filename//'.xdmf',form='formatted')
    write(fh,'(A)')'<?xml version="1.0" ?>'
    write(fh,'(A)')'<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
    write(fh,'(A)')'<Xdmf xmlns:xi="http://www.w3.org/2001/XInclude" Version="2.0">'
    write(fh,'(A)')'  <Domain>'
    !
    write(fh,'(A,3(1X,I0),A)')'    <Topology name="topo" TopologyType="3DRectMesh" Dimensions="',  &
                              km,jm,im,'"> </Topology>'
    write(fh,'(A)')'    <Geometry name="geo" Type="VXVYVZ">'
    write(fh,'(A)')'      <DataItem Format="Binary" DataType="Float" Precision="4" Endian="little" Seek="0" '
    write(fh,'(A,1X,I0,A)') '                Dimensions=" ',im,'">grid-x.bin  </DataItem>'
    write(fh,'(A)')'      <DataItem Format="Binary" DataType="Float" Precision="4" Endian="little" Seek="0" '
    write(fh,'(A,1X,I0,A)') '                Dimensions=" ',jm,'">grid-y.bin  </DataItem>'
    write(fh,'(A)')'      <DataItem Format="Binary" DataType="Float" Precision="4" Endian="little" Seek="0" '
    write(fh,'(A,1X,I0,A)') '                Dimensions=" ',km,'">grid-z.bin  </DataItem>'
    write(fh,'(A)')'    </Geometry>'
    write(fh,*)
    !
    write(fh,'(A)')'    <Grid Name="001" GridType="Uniform">'
    write(fh,'(A)')'      <Topology Reference="/Xdmf/Domain/Topology[1]"/>'
    write(fh,'(A)')'      <Geometry Reference="/Xdmf/Domain/Geometry[1]"/>'
    

    do n=1,nvar
      write(fh,'(3(A))')'      <Attribute Name="',varname(n),'" Center="Node">'
      write(fh,'(A)')'        <DataItem Format="Binary" DataType="Float" Precision="4" Endian="little" Seek="0"'
      write(fh,'(A,3(1X,I0),3(A))')'                   Dimensions="',km,jm,im,'"> ',filename//'-'//varname(n),'</DataItem>'
      write(fh,'(A)')'      </Attribute>'

      bufr4=var(:,:,:,n)
      open(fh+1,file=dir//filename//'-'//varname(n),access="stream")
      write(fh+1)bufr4
      close(fh+1)
      print*,' << ',dir//filename//'-'//varname(n)
    enddo

    write(fh,'(A)')'     </Grid>'
    !
    write(fh,'(A)')'  </Domain>'
    write(fh,'(A)')'</Xdmf>'
    !
    close(fh)
    !
    !
  end subroutine xdmfwriter_3drec_list

  subroutine xdmfwriter_3drec_block(dir,filename,block)
    
    use pastr_utility, only : make_dir
    use pastr_multiblock_type, only : block_type
    
    ! arguments
    character(len=*),intent(in) :: dir,filename
    type(block_type),intent(in),dimension(:) :: block
    !
    real(4),allocatable :: bufr4(:,:,:)
    !
    ! local data
    integer :: fh,i,n,m,nblk,nb
    integer :: im,jm,km
    character(len=4) :: blkname
    character(len=64) :: gridnamex,gridnamey,gridnamez,filavarname
    logical :: lfex
    
    call make_dir(dir)

    nblk=size(block)

    fh=18

    open(fh,file=dir//filename//'.xdmf',form='formatted')
    write(fh,'(A)')'<?xml version="1.0" ?>'
    write(fh,'(A)')'<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
    write(fh,'(A)')'<Xdmf xmlns:xi="http://www.w3.org/2001/XInclude" Version="2.0">'
    write(fh,'(A)')'  <Domain>'

    do n=1,nblk

      im=block(n)%im
      jm=block(n)%jm
      km=block(n)%km

      write(blkname,'(i4.4)')n

      gridnamex=filename//'-block'//blkname//'-gridx.bin'
      inquire(file=dir//trim(gridnamex),exist=lfex)
      if(.not. lfex) then
        open(17,file=dir//trim(gridnamex),access="stream")
        write(17)sngl(block(n)%x(0:im,0,0,1))
        close(17)
        print*,' << ',dir//trim(gridnamex)
      endif

      gridnamey=filename//'-block'//blkname//'-gridy.bin'
      inquire(file=trim(gridnamey),exist=lfex)
      if(.not. lfex) then
        open(17,file=dir//trim(gridnamey),access="stream")
        write(17)sngl(block(n)%x(0,0:jm,0,2))
        close(17)
        print*,' << ',dir//trim(gridnamey)
      endif

      gridnamez=filename//'-block'//blkname//'-gridz.bin'
      inquire(file=trim(gridnamez),exist=lfex)
      if(.not. lfex) then
        open(17,file=dir//trim(gridnamez),access="stream")
        write(17)sngl(block(n)%x(0,0,0:km,3))
        close(17)
        print*,' << ',dir//trim(gridnamez)
      endif

      ! write the head of xdmf and grid
      
      write(fh,'(3(A))')'    <Grid Name="',trim(block(n)%name),'" GridType="Uniform">'

      write(fh,'(A,3(1X,I0),A)')'    <Topology name="topo" TopologyType="3DRectMesh" Dimensions="',  &
                                km,jm,im,'"> </Topology>'
      write(fh,'(A)')'    <Geometry name="geo" Type="VXVYVZ">'
      write(fh,'(A)')'      <DataItem Format="Binary" DataType="Float" Precision="4" Endian="little" Seek="0" '
      write(fh,'(A,1X,I0,3(A))') '                Dimensions=" ',im+1,'"> ',gridnamex,'  </DataItem>'
      write(fh,'(A)')'      <DataItem Format="Binary" DataType="Float" Precision="4" Endian="little" Seek="0" '
      write(fh,'(A,1X,I0,3(A))') '                Dimensions=" ',jm+1,'"> ',gridnamey,' </DataItem>'
      write(fh,'(A)')'      <DataItem Format="Binary" DataType="Float" Precision="4" Endian="little" Seek="0" '
      write(fh,'(A,1X,I0,3(A))') '                Dimensions=" ',km+1,'"> ',gridnamez,'  </DataItem>'
      write(fh,'(A)')'    </Geometry>'
      write(fh,*)    

      allocate(bufr4(0:im,0:jm,0:km))

      do m=1,block(n)%nvar

        filavarname=filename//'-'//'block'//blkname//'-'//trim(block(n)%varname(m))

        write(fh,'(3(A))')'      <Attribute Name="',trim(block(n)%varname(m)),'" Center="Node">'
        write(fh,'(A)')'        <DataItem Format="Binary" DataType="Float" Precision="4" Endian="little" Seek="0"'
        write(fh,'(A,3(1X,I0),3(A))')'                   Dimensions="',km+1,jm+1,im+1,'"> ',trim(filavarname),'</DataItem>'
        write(fh,'(A)')'      </Attribute>'
  
        bufr4=sngl(block(n)%var(0:im,0:jm,0:km,m))
        open(fh+1,file=dir//trim(filavarname),access="stream")
        write(fh+1)bufr4
        close(fh+1)
        print*,' << ',dir//filename//'-'//trim(block(n)%varname(m))

        deallocate(bufr4)
      enddo

      write(fh,'(A)')'     </Grid>'
    enddo
    !
    write(fh,'(A)')'  </Domain>'
    write(fh,'(A)')'</Xdmf>'
    !
    close(fh)
    print*,' <<',dir//filename//'.xdmf done.'

  end subroutine xdmfwriter_3drec_block

  subroutine xdmfwriter_2d(dir,filename,gridfile,var1,var1name,var2,var2name, &
                                              var3,var3name,var4,var4name, &
                                              var5,var5name,var6,var6name, im,jm)
    !
    use pastr_h5io
    !
    ! arguments
    character(len=*),intent(in) :: dir,filename,gridfile
    real(wp),intent(in),dimension(:,:) :: var1
    character(len=*),intent(in) :: var1name
    integer :: im,jm
    real(wp),intent(in),dimension(:,:),optional :: var2,var3,var4,var5,var6
    character(len=*),intent(in),optional :: var2name,var3name,var4name,var5name,var6name
    !
    real(4),allocatable :: bufr4(:,:)
    !
    ! local data
    integer :: fh,i
    !
    allocate(bufr4(0:im,0:jm))
    !
    fh=18
    !
    ! write the head of xdmf and grid
    open(fh,file=dir//filename//'.xdmf',form='formatted')
    write(fh,'(A)')'<?xml version="1.0" ?>'
    write(fh,'(A)')'<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
    write(fh,'(A)')'<Xdmf xmlns:xi="http://www.w3.org/2001/XInclude" Version="2.0">'
    write(fh,'(A)')'  <Domain>'
    !
    write(fh,'(A,2(1X,I0),A)')'    <Topology name="topo" TopologyType="2DSMESH" Dimensions="',  &
                              jm+1,im+1,'"> </Topology>'
    write(fh,'(A)')'    <Geometry name="geo" Type="X_Y">'
    write(fh,'(A)')'      <DataItem Format="HDF" DataType="Float" Precision="4" Endian="little" Seek="0"'
    write(fh,'(A,2(1X,I0),3(A))')'                Dimensions="',jm+1,im+1,'"> ',gridfile,':x </DataItem>'
    write(fh,'(A)')'      <DataItem Format="HDF" DataType="Float" Precision="4" Endian="little" Seek="0"'
    write(fh,'(A,2(1X,I0),3(A))')'                Dimensions="',jm+1,im+1,'"> ',gridfile,':y </DataItem>'
    write(fh,'(A)')'    </Geometry>'
    !
    write(fh,'(A)')'    <Grid Name="001" GridType="Uniform">'
    write(fh,'(A)')'      <Topology Reference="/Xdmf/Domain/Topology[1]"/>'
    write(fh,'(A)')'      <Geometry Reference="/Xdmf/Domain/Geometry[1]"/>'
    !
    write(fh,'(3(A))')'      <Attribute Name="',var1name,'" Center="Node">'
    write(fh,'(A)')'        <DataItem Format="HDF" DataType="Float" Precision="4" Endian="little" Seek="0"'
    write(fh,'(A,2(1X,I0),5(A))')'                   Dimensions="',jm+1,im+1,'"> ',filename//'.h5',':',var1name,'</DataItem>'
    write(fh,'(A)')'      </Attribute>'
    bufr4=var1
    call h5_writearray2d_r4(bufr4,im,jm,var1name,dir//filename//'.h5')

    if(present(var2)) then
      write(fh,'(3(A))')'      <Attribute Name="',var2name,'" Center="Node">'
      write(fh,'(A)')'        <DataItem Format="HDF" DataType="Float" Precision="4" Endian="little" Seek="0"'
      write(fh,'(A,2(1X,I0),5(A))')'                   Dimensions="',jm+1,im+1,'"> ',filename//'.h5',':',var2name,'</DataItem>'
      write(fh,'(A)')'      </Attribute>'
      bufr4=var2
      call h5_writearray2d_r4(bufr4,im,jm,var2name,dir//filename//'.h5')
    endif
    if(present(var3)) then
      write(fh,'(3(A))')'      <Attribute Name="',var3name,'" Center="Node">'
      write(fh,'(A)')'        <DataItem Format="HDF" DataType="Float" Precision="4" Endian="little" Seek="0"'
      write(fh,'(A,2(1X,I0),5(A))')'                   Dimensions="',jm+1,im+1,'"> ',filename//'.h5',':',var3name,'</DataItem>'
      write(fh,'(A)')'      </Attribute>'
      bufr4=var3
      call h5_writearray2d_r4(bufr4,im,jm,var3name,dir//filename//'.h5')
    endif
    if(present(var4)) then
      write(fh,'(3(A))')'      <Attribute Name="',var4name,'" Center="Node">'
      write(fh,'(A)')'        <DataItem Format="HDF" DataType="Float" Precision="4" Endian="little" Seek="0"'
      write(fh,'(A,2(1X,I0),5(A))')'                   Dimensions="',jm+1,im+1,'"> ',filename//'.h5',':',var4name,'</DataItem>'
      write(fh,'(A)')'      </Attribute>'
      bufr4=var4
      call h5_writearray2d_r4(bufr4,im,jm,var4name,dir//filename//'.h5')
    endif
    if(present(var5)) then
      write(fh,'(3(A))')'      <Attribute Name="',var5name,'" Center="Node">'
      write(fh,'(A)')'        <DataItem Format="HDF" DataType="Float" Precision="4" Endian="little" Seek="0"'
      write(fh,'(A,2(1X,I0),5(A))')'                   Dimensions="',jm+1,im+1,'"> ',filename//'.h5',':',var5name,'</DataItem>'
      write(fh,'(A)')'      </Attribute>'
      bufr4=var5
      call h5_writearray2d_r4(bufr4,im,jm,var5name,dir//filename//'.h5')
    endif
    if(present(var6)) then
      write(fh,'(3(A))')'      <Attribute Name="',var6name,'" Center="Node">'
      write(fh,'(A)')'        <DataItem Format="HDF" DataType="Float" Precision="4" Endian="little" Seek="0"'
      write(fh,'(A,2(1X,I0),5(A))')'                   Dimensions="',jm+1,im+1,'"> ',filename//'.h5',':',var6name,'</DataItem>'
      write(fh,'(A)')'      </Attribute>'
      bufr4=var6
      call h5_writearray2d_r4(bufr4,im,jm,var6name,dir//filename//'.h5')
    endif

    write(fh,'(A)')'     </Grid>'
    !
    write(fh,'(A)')'  </Domain>'
    write(fh,'(A)')'</Xdmf>'
    !
    close(fh)
    !
    !
  end subroutine xdmfwriter_2d

  subroutine xdmfwriter_2drec(dir,filename,x,y,var1,var1name,var2,var2name, &
                                               var3,var3name,var4,var4name, &
                                               var5,var5name,var6,var6name, &
                                               var7,var7name,var8,var8name, &
                                               var9,var9name,var10,var10name)
    !
    ! arguments
    character(len=*),intent(in) :: dir,filename
    real(wp),intent(in),dimension(:) :: x,y
    real(wp),intent(in),dimension(:,:) :: var1
    character(len=*),intent(in) :: var1name

    real(wp),intent(in),dimension(:,:),optional :: var2,var3,var4,var5,var6,var7,var8,var9,var10
    character(len=*),intent(in),optional :: var2name,var3name,var4name,var5name,var6name, &
                                            var7name,var8name,var9name,var10name
    !
    real(4),allocatable :: bufr4(:,:,:)
    !
    ! local data
    integer :: fh,i
    integer :: im,jm
    logical :: lfex
    !
    im=size(var1,1)
    jm=size(var1,2)

    inquire(file=dir//'/grid-x.bin',exist=lfex)
    if(.not. lfex) then
      open(17,file=dir//'/grid-x.bin',access="stream")
      write(17)sngl(x)
      close(17)
      print*,' << ',dir//'/grid-x.bin'
    endif
    inquire(file=dir//'/grid-y.bin',exist=lfex)
    if(.not. lfex) then
      open(17,file=dir//'/grid-y.bin',access="stream")
      write(17)sngl(y)
      close(17)
      print*,' << ',dir//'/grid-y.bin'
    endif

    !
    allocate(bufr4(im,jm,1))
    !
    fh=18
    !
    ! write the head of xdmf and grid
    open(fh,file=dir//filename//'.xdmf',form='formatted')
    write(fh,'(A)')'<?xml version="1.0" ?>'
    write(fh,'(A)')'<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
    write(fh,'(A)')'<Xdmf xmlns:xi="http://www.w3.org/2001/XInclude" Version="2.0">'
    write(fh,'(A)')'  <Domain>'
    !
    write(fh,'(A,2(1X,I0),A)')'    <Topology name="topo" TopologyType="3DRectMesh" Dimensions=" 1 ',  &
                              jm,im,'"> </Topology>'
    write(fh,'(A)')'    <Geometry name="geo" Type="VXVYVZ">'
    write(fh,'(A)')'      <DataItem Format="Binary" DataType="Float" Precision="4" Endian="little" Seek="0" '
    write(fh,'(A,1X,I0,A)') '                Dimensions=" ',im,'">grid-x.bin  </DataItem>'
    write(fh,'(A)')'      <DataItem Format="Binary" DataType="Float" Precision="4" Endian="little" Seek="0" '
    write(fh,'(A,1X,I0,A)') '                Dimensions=" ',jm,'">grid-y.bin  </DataItem>'
    write(fh,'(A)')'      <DataItem Dimensions=" 1 " Numbewp="Float" Precision="4" Format="XML" > 0.0 </DataItem>'
    write(fh,'(A)')'    </Geometry>'
    write(fh,*)
    !
    write(fh,'(A)')'    <Grid Name="001" GridType="Uniform">'
    write(fh,'(A)')'      <Topology Reference="/Xdmf/Domain/Topology[1]"/>'
    write(fh,'(A)')'      <Geometry Reference="/Xdmf/Domain/Geometry[1]"/>'
    !
    write(fh,'(3(A))')'      <Attribute Name="',var1name,'" Center="Node">'
    write(fh,'(A)')'        <DataItem Format="Binary" DataType="Float" Precision="4" Endian="little" Seek="0"'
    write(fh,'(A,2(1X,I0),3(A))')'                   Dimensions=" 1 ',jm,im,'"> ',filename//'-'//var1name,'</DataItem>'
    write(fh,'(A)')'      </Attribute>'

    bufr4(:,:,1)=var1(:,:)
    open(fh+1,file=dir//filename//'-'//var1name,access="stream")
    write(fh+1)bufr4
    close(fh+1)
    print*,' << ',dir//filename//'-'//var1name

    if(present(var2)) then
      write(fh,'(3(A))')'      <Attribute Name="',var2name,'" Center="Node">'
      write(fh,'(A)')'        <DataItem Format="Binary" DataType="Float" Precision="4" Endian="little" Seek="0"'
      write(fh,'(A,2(1X,I0),3(A))')'                   Dimensions=" 1 ',jm,im,'"> ',filename//'-'//var2name,'</DataItem>'
      write(fh,'(A)')'      </Attribute>'
  
      bufr4(:,:,1)=var2(:,:)
      open(fh+1,file=dir//filename//'-'//var2name,access="stream")
      write(fh+1)bufr4
      close(fh+1)
      print*,' << ',dir//filename//'-'//var2name
    endif
    if(present(var3)) then
      write(fh,'(3(A))')'      <Attribute Name="',var3name,'" Center="Node">'
      write(fh,'(A)')'        <DataItem Format="Binary" DataType="Float" Precision="4" Endian="little" Seek="0"'
      write(fh,'(A,2(1X,I0),3(A))')'                   Dimensions=" 1 ',jm,im,'"> ',filename//'-'//var3name,'</DataItem>'
      write(fh,'(A)')'      </Attribute>'
  
      bufr4(:,:,1)=var3(:,:)
      open(fh+1,file=dir//filename//'-'//var3name,access="stream")
      write(fh+1)bufr4
      close(fh+1)
      print*,' << ',dir//filename//'-'//var3name
    endif
    if(present(var4)) then
      write(fh,'(3(A))')'      <Attribute Name="',var4name,'" Center="Node">'
      write(fh,'(A)')'        <DataItem Format="Binary" DataType="Float" Precision="4" Endian="little" Seek="0"'
      write(fh,'(A,2(1X,I0),3(A))')'                   Dimensions=" 1 ',jm,im,'"> ',filename//'-'//var4name,'</DataItem>'
      write(fh,'(A)')'      </Attribute>'
  
      bufr4(:,:,1)=var4(:,:)
      open(fh+1,file=dir//filename//'-'//var4name,access="stream")
      write(fh+1)bufr4
      close(fh+1)
      print*,' << ',dir//filename//'-'//var4name
    endif
    if(present(var5)) then
      write(fh,'(3(A))')'      <Attribute Name="',var5name,'" Center="Node">'
      write(fh,'(A)')'        <DataItem Format="Binary" DataType="Float" Precision="4" Endian="little" Seek="0"'
      write(fh,'(A,2(1X,I0),3(A))')'                   Dimensions=" 1 ',jm,im,'"> ',filename//'-'//var5name,'</DataItem>'
      write(fh,'(A)')'      </Attribute>'
  
      bufr4(:,:,1)=var5(:,:)
      open(fh+1,file=dir//filename//'-'//var5name,access="stream")
      write(fh+1)bufr4
      close(fh+1)
      print*,' << ',dir//filename//'-'//var5name
    endif
    if(present(var6)) then
      write(fh,'(3(A))')'      <Attribute Name="',var6name,'" Center="Node">'
      write(fh,'(A)')'        <DataItem Format="Binary" DataType="Float" Precision="4" Endian="little" Seek="0"'
      write(fh,'(A,2(1X,I0),3(A))')'                   Dimensions=" 1 ',jm,im,'"> ',filename//'-'//var6name,'</DataItem>'
      write(fh,'(A)')'      </Attribute>'
  
      bufr4(:,:,1)=var6(:,:)
      open(fh+1,file=dir//filename//'-'//var6name,access="stream")
      write(fh+1)bufr4
      close(fh+1)
      print*,' << ',dir//filename//'-'//var6name
    endif
    if(present(var7)) then
      write(fh,'(3(A))')'      <Attribute Name="',var7name,'" Center="Node">'
      write(fh,'(A)')'        <DataItem Format="Binary" DataType="Float" Precision="4" Endian="little" Seek="0"'
      write(fh,'(A,2(1X,I0),3(A))')'                   Dimensions=" 1 ',jm,im,'"> ',filename//'-'//var7name,'</DataItem>'
      write(fh,'(A)')'      </Attribute>'
  
      bufr4(:,:,1)=var7(:,:)
      open(fh+1,file=dir//filename//'-'//var7name,access="stream")
      write(fh+1)bufr4
      close(fh+1)
      print*,' << ',dir//filename//'-'//var7name
    endif
    if(present(var8)) then
      write(fh,'(3(A))')'      <Attribute Name="',var8name,'" Center="Node">'
      write(fh,'(A)')'        <DataItem Format="Binary" DataType="Float" Precision="4" Endian="little" Seek="0"'
      write(fh,'(A,2(1X,I0),3(A))')'                   Dimensions=" 1 ',jm,im,'"> ',filename//'-'//var8name,'</DataItem>'
      write(fh,'(A)')'      </Attribute>'
  
      bufr4(:,:,1)=var8(:,:)
      open(fh+1,file=dir//filename//'-'//var8name,access="stream")
      write(fh+1)bufr4
      close(fh+1)
      print*,' << ',dir//filename//'-'//var8name
    endif
    if(present(var9)) then
      write(fh,'(3(A))')'      <Attribute Name="',var9name,'" Center="Node">'
      write(fh,'(A)')'        <DataItem Format="Binary" DataType="Float" Precision="4" Endian="little" Seek="0"'
      write(fh,'(A,2(1X,I0),3(A))')'                   Dimensions=" 1 ',jm,im,'"> ',filename//'-'//var9name,'</DataItem>'
      write(fh,'(A)')'      </Attribute>'
  
      bufr4(:,:,1)=var9(:,:)
      open(fh+1,file=dir//filename//'-'//var9name,access="stream")
      write(fh+1)bufr4
      close(fh+1)
      print*,' << ',dir//filename//'-'//var9name
    endif
    if(present(var10)) then
      write(fh,'(3(A))')'      <Attribute Name="',var10name,'" Center="Node">'
      write(fh,'(A)')'        <DataItem Format="Binary" DataType="Float" Precision="4" Endian="little" Seek="0"'
      write(fh,'(A,2(1X,I0),3(A))')'                   Dimensions=" 1 ',jm,im,'"> ',filename//'-'//var10name,'</DataItem>'
      write(fh,'(A)')'      </Attribute>'
  
      bufr4(:,:,1)=var10(:,:)
      open(fh+1,file=dir//filename//'-'//var10name,access="stream")
      write(fh+1)bufr4
      close(fh+1)
      print*,' << ',dir//filename//'-'//var10name
    endif

    write(fh,'(A)')'     </Grid>'
    !
    write(fh,'(A)')'  </Domain>'
    write(fh,'(A)')'</Xdmf>'
    !
    close(fh)
    !
    !
  end subroutine xdmfwriter_2drec

  subroutine xdmfwriter_2drec_list_xy(dir,filename,x,y,var,varname)
    
    use pastr_utility, only : make_dir
    ! arguments
    character(len=*),intent(in) :: dir,filename
    real(wp),intent(in),dimension(:) :: x,y
    real(wp),intent(in) :: var(:,:,:)
    character(len=*),intent(in) :: varname(:)

    real(4),allocatable :: bufr4(:,:,:)
    !
    ! local data
    integer :: fh,i,n,nvar
    integer :: im,jm
    logical :: lfex
    
    call make_dir(dir)

    im=size(var,1)
    jm=size(var,2)

    inquire(file=dir//'/grid-x.bin',exist=lfex)
    if(.not. lfex) then
      open(17,file=dir//'/grid-x.bin',access="stream")
      write(17)sngl(x)
      close(17)
      print*,' << ',dir//'/grid-x.bin'
    endif
    inquire(file=dir//'/grid-y.bin',exist=lfex)
    if(.not. lfex) then
      open(17,file=dir//'/grid-y.bin',access="stream")
      write(17)sngl(y)
      close(17)
      print*,' << ',dir//'/grid-y.bin'
    endif

    allocate(bufr4(im,jm,1))
    
    fh=18
    !
    ! write the head of xdmf and grid
    open(fh,file=dir//filename//'.xdmf',form='formatted')
    write(fh,'(A)')'<?xml version="1.0" ?>'
    write(fh,'(A)')'<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
    write(fh,'(A)')'<Xdmf xmlns:xi="http://www.w3.org/2001/XInclude" Version="2.0">'
    write(fh,'(A)')'  <Domain>'
    !
    write(fh,'(A,2(1X,I0),A)')'    <Topology name="topo" TopologyType="3DRectMesh" Dimensions=" 1 ',  &
                              jm,im,'"> </Topology>'
    write(fh,'(A)')'    <Geometry name="geo" Type="VXVYVZ">'
    write(fh,'(A)')'      <DataItem Format="Binary" DataType="Float" Precision="4" Endian="little" Seek="0" '
    write(fh,'(A,1X,I0,A)') '                Dimensions=" ',im,'">grid-x.bin  </DataItem>'
    write(fh,'(A)')'      <DataItem Format="Binary" DataType="Float" Precision="4" Endian="little" Seek="0" '
    write(fh,'(A,1X,I0,A)') '                Dimensions=" ',jm,'">grid-y.bin  </DataItem>'
    write(fh,'(A)')'      <DataItem Dimensions=" 1 " Numbewp="Float" Precision="4" Format="XML" > 0.0 </DataItem>'
    write(fh,'(A)')'    </Geometry>'
    write(fh,*)
    !
    write(fh,'(A)')'    <Grid Name="001" GridType="Uniform">'
    write(fh,'(A)')'      <Topology Reference="/Xdmf/Domain/Topology[1]"/>'
    write(fh,'(A)')'      <Geometry Reference="/Xdmf/Domain/Geometry[1]"/>'
    
    nvar=size(var,3)
    do n=1,nvar
      write(fh,'(3(A))')'      <Attribute Name="',varname(n),'" Center="Node">'
      write(fh,'(A)')'        <DataItem Format="Binary" DataType="Float" Precision="4" Endian="little" Seek="0"'
      write(fh,'(A,2(1X,I0),3(A))')'                   Dimensions=" 1 ',jm,im,'"> ',filename//'-'//varname(n),'</DataItem>'
      write(fh,'(A)')'      </Attribute>'

      bufr4(:,:,1)=var(:,:,n)
      open(fh+1,file=dir//filename//'-'//varname(n),access="stream")
      write(fh+1)bufr4
      close(fh+1)
      print*,' << ',dir//filename//'-'//varname(n)
    enddo

    write(fh,'(A)')'     </Grid>'
    !
    write(fh,'(A)')'  </Domain>'
    write(fh,'(A)')'</Xdmf>'
    !
    close(fh)
    !
    !
  end subroutine xdmfwriter_2drec_list_xy

  subroutine xdmfwriter_2drec_block(dir,filename,blocks)
    
    use pastr_utility, only : make_dir
    use pastr_commtype,only : tblock

    ! arguments
    character(len=*),intent(in) :: dir,filename
    type(tblock),intent(in),dimension(:) :: blocks

    real(4),allocatable :: bufr4(:,:)
    character(len=3) :: bname
    character(len=64) :: gridnamex,gridnamey,filavarname
    !
    ! local data
    integer :: fh,i,m,nvar
    integer :: im,jm
    logical :: lfex
    
    call make_dir(dir)

    fh=18
    open(fh,file=dir//filename//'.xdmf',form='formatted')
    write(fh,'(A)')'<?xml version="1.0" ?>'
    write(fh,'(A)')'<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
    write(fh,'(A)')'<Xdmf xmlns:xi="http://www.w3.org/2001/XInclude" Version="2.0">'
    write(fh,'(A)')'  <Domain>'

    do i=1,size(blocks)

      im=blocks(i)%im
      jm=blocks(i)%jm
      
      write(bname,'(i3.3)') i
      gridnamex='block'//bname//'-grid-x.bin'
      gridnamey='block'//bname//'-grid-y.bin'

      inquire(file=dir//trim(gridnamex),exist=lfex)
      if(.not. lfex) then
        open(17,file=dir//trim(gridnamex),access="stream")
        write(17)sngl(blocks(i)%x(0:im,0,0))
        close(17)
        print*,' << ',dir//trim(gridnamex)
      endif
      inquire(file=dir//trim(gridnamey),exist=lfex)
      if(.not. lfex) then
        open(17,file=dir//trim(gridnamey),access="stream")
        write(17)sngl(blocks(i)%y(0,0:jm,0))
        close(17)
        print*,' << ',dir//trim(gridnamey)
      endif

      write(fh,'(3(A))')'    <Grid Name=" block-',bname,'" GridType="Uniform">'

      write(fh,'(A,2(1X,I0),A)')'    <Topology name="topo" TopologyType="3DRectMesh" Dimensions=" 1',  &
                                jm+1,im+1,'"> </Topology>'
      write(fh,'(A)')'    <Geometry name="geo" Type="VXVYVZ">'
      write(fh,'(A)')'      <DataItem Format="Binary" DataType="Float" Precision="4" Endian="little" Seek="0" '
      write(fh,'(A,1X,I0,3(A))') '                Dimensions=" ',im+1,'"> ',gridnamex,'  </DataItem>'
      write(fh,'(A)')'      <DataItem Format="Binary" DataType="Float" Precision="4" Endian="little" Seek="0" '
      write(fh,'(A,1X,I0,3(A))') '                Dimensions=" ',jm+1,'"> ',gridnamey,' </DataItem>'
      write(fh,'(A)')'      <DataItem Dimensions=" 1 " Numbewp="Float" Precision="4" Format="XML" > 0.0 </DataItem>'
      write(fh,'(A)')'    </Geometry>'
      write(fh,*)    

      allocate(bufr4(0:im,0:jm))

      do m=1,blocks(i)%nvar

        filavarname=filename//'-'//'block'//bname//'-'//trim(blocks(i)%varname(m))

        write(fh,'(3(A))')'      <Attribute Name="',trim(blocks(i)%varname(m)),'" Center="Node">'
        write(fh,'(A)')'        <DataItem Format="Binary" DataType="Float" Precision="4" Endian="little" Seek="0"'
        write(fh,'(A,2(1X,I0),3(A))')'                   Dimensions=" 1 ',jm+1,im+1,'"> ',trim(filavarname),'</DataItem>'

        write(fh,'(A)')'      </Attribute>'
        
        bufr4=sngl(blocks(i)%var(0:im,0:jm,0,m))
        open(fh+1,file=dir//trim(filavarname),access="stream")
        write(fh+1)bufr4
        close(fh+1)
        print*,' << ',dir//filename//'-'//trim(blocks(i)%varname(m))

        deallocate(bufr4)
      enddo

      write(fh,'(A)')'     </Grid>'
    enddo
    !
    write(fh,'(A)')'  </Domain>'
    write(fh,'(A)')'</Xdmf>'
    !
    close(fh)
    print*,' <<',dir//filename//'.xdmf done.'

  end subroutine xdmfwriter_2drec_block

end module pastr_xdmf