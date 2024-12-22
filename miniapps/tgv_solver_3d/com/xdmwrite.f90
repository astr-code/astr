module xdmwrite
   
  implicit none

  Interface xdmfwriter
    module procedure xdmfwriter_2d
    module procedure xdmfwriter_2drec
    module procedure xdmfwriter_3d
    module procedure xdmfwriter_3drec
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
    !
    ! arguments
    character(len=*),intent(in) :: dir,filename,gridfile
    real(8),intent(in),dimension(:,:,:) :: var1
    character(len=*),intent(in) :: var1name
    integer :: im,jm,km
    real(8),intent(in),dimension(:,:,:),optional :: var2,var3,var4,var5,var6
    character(len=*),intent(in),optional :: var2name,var3name,var4name,var5name,var6name
    !
  end subroutine xdmfwriter_3d
  !
  subroutine xdmfwriter_3dbox(dir,filename,deltax,var1,var1name,var2,var2name, &
                                                  var3,var3name,var4,var4name, &
                                                  var5,var5name,var6,var6name, im,jm,km)
    !
    ! arguments
    character(len=*),intent(in) :: dir,filename
    real(8),intent(in) :: deltax
    real(8),intent(in),dimension(:,:,:) :: var1
    character(len=*),intent(in) :: var1name
    integer :: im,jm,km
    real(8),intent(in),dimension(:,:,:),optional :: var2,var3,var4,var5,var6
    character(len=*),intent(in),optional :: var2name,var3name,var4name,var5name,var6name
    !
    real(4),allocatable :: bufr4(:,:,:)
    !
    ! local data
    integer :: fh
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

    bufr4=sngl(var1)
    open(fh+1,file=dir//filename//'-'//var1name,access="stream")
    write(fh+1)bufr4
    close(fh+1)
    print*,' << ',dir//filename//'-'//var1name

    if(present(var2)) then
      write(fh,'(3(A))')'      <Attribute Name="',var2name,'" Center="Node">'
      write(fh,'(A)')'        <DataItem Format="Binary" DataType="Float" Precision="4" Endian="little" Seek="0"'
      write(fh,'(A,3(1X,I0),3(A))')'                   Dimensions="',km+1,jm+1,im+1,'"> ',filename//'-'//var2name,'</DataItem>'
      write(fh,'(A)')'      </Attribute>'
  
      bufr4=sngl(var2)
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
  
      bufr4=sngl(var3)
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
  
      bufr4=sngl(var4)
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
  
      bufr4=sngl(var5)
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
  
      bufr4=sngl(var6)
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
    real(8),intent(in),dimension(:) :: x,y,z
    real(8),intent(in),dimension(:,:,:) :: var1
    character(len=*),intent(in) :: var1name

    real(8),intent(in),dimension(:,:,:),optional :: var2,var3,var4,var5,var6
    character(len=*),intent(in),optional :: var2name,var3name,var4name,var5name,var6name
    !
    real(4),allocatable :: bufr4(:,:,:)
    !
    ! local data
    integer :: fh
    integer :: im,jm,km
    logical :: lfex
    !
    im=size(var1,1)
    jm=size(var1,2)
    km=size(var1,3)

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
    !
    write(fh,'(3(A))')'      <Attribute Name="',var1name,'" Center="Node">'
    write(fh,'(A)')'        <DataItem Format="Binary" DataType="Float" Precision="4" Endian="little" Seek="0"'
    write(fh,'(A,3(1X,I0),3(A))')'                   Dimensions="',km,jm,im,'"> ',filename//'-'//var1name,'</DataItem>'
    write(fh,'(A)')'      </Attribute>'

    bufr4=sngl(var1)
    open(fh+1,file=dir//filename//'-'//var1name,access="stream")
    write(fh+1)bufr4
    close(fh+1)
    print*,' << ',dir//filename//'-'//var1name

    if(present(var2)) then
      write(fh,'(3(A))')'      <Attribute Name="',var2name,'" Center="Node">'
      write(fh,'(A)')'        <DataItem Format="Binary" DataType="Float" Precision="4" Endian="little" Seek="0"'
      write(fh,'(A,3(1X,I0),3(A))')'                   Dimensions="',km,jm,im,'"> ',filename//'-'//var2name,'</DataItem>'
      write(fh,'(A)')'      </Attribute>'
  
      bufr4=sngl(var2)
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
  
      bufr4=sngl(var3)
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
  
      bufr4=sngl(var4)
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
  
      bufr4=sngl(var5)
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
  
      bufr4=sngl(var6)
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

  subroutine xdmfwriter_2d(dir,filename,gridfile,var1,var1name,var2,var2name, &
                                              var3,var3name,var4,var4name, &
                                              var5,var5name,var6,var6name, im,jm)
    !
    ! arguments
    character(len=*),intent(in) :: dir,filename,gridfile
    real(8),intent(in),dimension(:,:) :: var1
    character(len=*),intent(in) :: var1name
    integer :: im,jm
    real(8),intent(in),dimension(:,:),optional :: var2,var3,var4,var5,var6
    character(len=*),intent(in),optional :: var2name,var3name,var4name,var5name,var6name
    !
    !
  end subroutine xdmfwriter_2d

  subroutine xdmfwriter_2drec(dir,filename,x,y,var1,var1name,var2,var2name, &
                                               var3,var3name,var4,var4name, &
                                               var5,var5name,var6,var6name)
    !
    ! arguments
    character(len=*),intent(in) :: dir,filename
    real(8),intent(in),dimension(:) :: x,y
    real(8),intent(in),dimension(:,:) :: var1
    character(len=*),intent(in) :: var1name

    real(8),intent(in),dimension(:,:),optional :: var2,var3,var4,var5,var6
    character(len=*),intent(in),optional :: var2name,var3name,var4name,var5name,var6name
    !
    real(4),allocatable :: bufr4(:,:,:)
    !
    ! local data
    integer :: fh
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
    write(fh,'(A)')'      <DataItem Dimensions=" 1 " NumberType="Float" Precision="4" Format="XML" > 0.0 </DataItem>'
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

    bufr4(:,:,1)=sngl(var1(:,:))
    open(fh+1,file=dir//filename//'-'//var1name,access="stream")
    write(fh+1)bufr4
    close(fh+1)
    print*,' << ',dir//filename//'-'//var1name

    if(present(var2)) then
      write(fh,'(3(A))')'      <Attribute Name="',var2name,'" Center="Node">'
      write(fh,'(A)')'        <DataItem Format="Binary" DataType="Float" Precision="4" Endian="little" Seek="0"'
      write(fh,'(A,2(1X,I0),3(A))')'                   Dimensions=" 1 ',jm,im,'"> ',filename//'-'//var2name,'</DataItem>'
      write(fh,'(A)')'      </Attribute>'
  
      bufr4(:,:,1)=sngl(var2(:,:))
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
  
      bufr4(:,:,1)=sngl(var3(:,:))
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
  
      bufr4(:,:,1)=sngl(var4(:,:))
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
  
      bufr4(:,:,1)=sngl(var5(:,:))
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
  
      bufr4(:,:,1)=sngl(var6(:,:))
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
  end subroutine xdmfwriter_2drec
  !+-------------------------------------------------------------------+
  !| The end of the subroutine xdmfwriter.                             |
  !+-------------------------------------------------------------------+

end module xdmwrite