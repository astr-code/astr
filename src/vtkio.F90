!+---------------------------------------------------------------------+
!|This module is used to write field in tecplot formate.                |
!+---------------------------------------------------------------------+
module WriteVTK
  !
  implicit none
  !
  Interface writeprvbin
    !
    module procedure writeprvbin2d7var
        !
    module procedure writeprvbin3d9var
    module procedure writeprvbin3d4var
    module procedure writeprvbin3d6var
    !
  end Interface writeprvbin
  !
  contains
  !
  !+-------------------------------------------------------------------+
  !|This subroutine is used to write bin file for 3d paraview field.   |
  !|ifort compiler only                                                |
  !+-------------------------------------------------------------------+
  subroutine writeprvbin3d4var(filename,var1,var1name,var2,var2name,   &
    var3,var3name,var4,var4name,imax,jmax,kmax)
    !
    character(len=*),intent(in) :: filename
    integer,intent(in) :: imax,jmax,kmax
    real(8),intent(in) :: var1(0:imax,0:jmax,0:kmax),                 &
    var2(0:imax,0:jmax,0:kmax),                 &
    var3(0:imax,0:jmax,0:kmax),                 &
    var4(0:imax,0:jmax,0:kmax)
    character(len=*),intent(in) :: var1name,var2name,var3name,var4name
    !
    integer :: int32,nbrvar,jvar
    integer :: int_bytesize,real_bytesize,nbt_scal,nbt_vec,nbt_xyz,nbt_ien,nbt_offset,nbt_etype
    integer(16) :: ioff
    ! ip : le point actuel
    !
    real(4) :: float32
    !
    real(4),allocatable,dimension(:,:,:,:) :: var,points
    character(256),allocatable,dimension(:) :: vname
    !
    character(15) :: charxmax,charymax,charzmax
    character :: lf*1,offset*16
    !        
    nbrvar=4
    !
    allocate(var(0:imax,0:jmax,0:kmax,nbrvar),points(3,0:imax,0:jmax,0:kmax))
    allocate(vname(nbrvar))
    !
    var(0:imax,0:jmax,0:kmax,1)=real(var1(0:imax,0:jmax,0:kmax))
    var(0:imax,0:jmax,0:kmax,2)=real(var2(0:imax,0:jmax,0:kmax))
    var(0:imax,0:jmax,0:kmax,3)=real(var3(0:imax,0:jmax,0:kmax))
    var(0:imax,0:jmax,0:kmax,4)=real(var4(0:imax,0:jmax,0:kmax))
    !
    points(1,0:imax,0:jmax,0:kmax)=var(:,:,:,1)
    points(2,0:imax,0:jmax,0:kmax)=var(:,:,:,2)
    points(3,0:imax,0:jmax,0:kmax)=var(:,:,:,3)
    !
    vname(1)=var1name
    vname(2)=var2name
    vname(3)=var3name
    vname(4)=var4name
    !
    write(charxmax,'(I4)') imax
    write(charymax,'(I4)') jmax
    write(charzmax,'(I4)') kmax
    ! print*,charxmin
    ! stop

    ! open(18,file=filename)
    ! write(18,*)'<VTKFile type="StructuredGrid">'
    ! write(18,*)'  <StructuredGrid WholeExtent="' &
    !               //'0'//trim(charxmax)//' 0' &
    !               //trim(charymax)//' 0'//trim(charzmax)//'">'
    ! write(18,*)'    <Piece Extent="' &
    !               //'0'//trim(charxmax)//' 0' &
    !               //trim(charymax)//' 0'//trim(charzmax)//'">'
    ! write(18,*)'      <PointData>'
    ! write(18,*)'        <DataArray'
    ! write(18,*)'          type="Float32"'
    ! write(18,*)'          Name="'//trim(vname(4))//'"'
    ! write(18,*)'          format="'//'ascii'//'"'
    ! write(18,*)'        >'
    ! write(18,'(9(1PE12.4))')(var(:,:,:,4))
    ! write(18,*)'        </DataArray>'
    ! write(18,*)'      </PointData>'
    ! write(18,*)'      <CellData>'
    ! write(18,*)'      </CellData>'
    ! write(18,*)'      <Points>'
    ! write(18,*)'        <DataArray NumberOfComponents="3"'
    ! write(18,*)'          type="Float32"'
    ! write(18,*)'          format="'//'ascii'//'"'
    ! write(18,*)'        >'
    ! write(18,'(9(1PE12.4))')(points(:,:,:,:))
    ! write(18,*)'        </DataArray>'
    ! write(18,*)'      </Points>'
    ! write(18,*)'    </Piece>'
    ! write(18,*)'  </StructuredGrid>'
    ! write(18,*)'</VTKFile>'
    !
    ! Layout for the Appended Binary Data
    lf = char(10)
    inquire(iolength=real_bytesize) float32
    nbt_scal   = (imax+1) * (jmax+1) * (kmax+1) * real_bytesize
    nbt_vec    = 3 * (imax+1) * (jmax+1) * (kmax+1) * real_bytesize
    inquire(iolength=int_bytesize) int32
    ! nbt_ien    = 8 * 8  * int_bytesize
    ! nbt_offset = 8  * int_bytesize
    ! nbt_etype  = 8  * int_bytesize
    !
    open(18,file=filename,access='stream',convert='LITTLE_ENDIAN')
    write(18)'<VTKFile type="StructuredGrid" byte_order="LittleEndian">'//lf
    write(18)'  <StructuredGrid WholeExtent="' &
    //'0'//trim(charxmax)//' 0' &
    //trim(charymax)//' 0'//trim(charzmax)//'">'//lf
    write(18)'    <Piece Extent="' &
    //'0'//trim(charxmax)//' 0' &
    //trim(charymax)//' 0'//trim(charzmax)//'">'//lf
    write(18)'      <PointData>'//lf
    !
    ioff=0
    do jvar=4,nbrvar
    if(jvar>4)ioff=ioff + int_bytesize + nbt_scal
    write(offset(1:16),'(I16)') ioff
    write(18)'        <DataArray'//lf
    write(18)'          type="Float32"'//lf
    write(18)'          Name="'//trim(vname(jvar))//'"'//lf
    write(18)'          format="appended"'//lf
    write(18)'          offset="'//offset//'"'//lf
    write(18)'        />'//lf
    enddo
    write(18)'      </PointData>'//lf
    write(18)'      <CellData>'//lf
    write(18)'      </CellData>'//lf 
    !
    write(18)'      <Points>'//lf


    ioff=ioff + int_bytesize + nbt_scal
    write(offset(1:16),'(I16)') ioff
    !
    write(18)'        <DataArray NumberOfComponents="3"'//lf
    write(18)'          type="Float32"'//lf
    write(18)'          format="appended"'//lf
    write(18)'          offset="'//offset//'"'//lf
    write(18)'        />'//lf
    write(18)'      </Points>'//lf
    write(18)'    </Piece>'//lf
    write(18)'  </StructuredGrid>'//lf
    write(18)'  <AppendedData encoding="raw">'//lf
    write(18)'_'
    !
    do jvar=4,nbrvar
    write(18)nbt_scal,var(:,:,:,jvar)
    enddo 
    !
    write(18)nbt_vec,points(:,:,:,:)
    !
    write(18) lf//'  </AppendedData>'//lf
    write(18)'</VTKFile>'//lf
    !
    close(18)
    !
    print*,' << ',filename
    !
    deallocate(var)
  end subroutine writeprvbin3d4var
  !
  subroutine writeprvbin3d6var(filename,var1,var1name,var2,var2name,   &
    var3,var3name,var4,var4name, &
    var5,var5name,var6,var6name, &
    imax,jmax,kmax)
    !
    character(len=*),intent(in) :: filename
    integer,intent(in) :: imax,jmax,kmax
    real(8),intent(in) :: var1(0:imax,0:jmax,0:kmax),                 &
    var2(0:imax,0:jmax,0:kmax),                 &
    var3(0:imax,0:jmax,0:kmax),                 &
    var4(0:imax,0:jmax,0:kmax),                 &
    var5(0:imax,0:jmax,0:kmax),                 &
    var6(0:imax,0:jmax,0:kmax)
    character(len=*),intent(in) :: var1name,var2name,var3name,var4name,var5name,var6name
    !
    integer :: int32,nbrvar,jvar
    integer :: int_bytesize,real_bytesize,nbt_scal,nbt_vec!,nbt_xyz,nbt_ien,nbt_offset,nbt_etype
    integer(16) :: ioff
    ! ip : le point actuel
    !
    real(4) :: float32
    !
    real(4),allocatable,dimension(:,:,:,:) :: var,points
    character(256),allocatable,dimension(:) :: vname
    !
    character(15) :: charxmax,charymax,charzmax
    character :: lf*1,offset*16
    !        
    nbrvar=6
    !
    allocate(var(0:imax,0:jmax,0:kmax,nbrvar),points(3,0:imax,0:jmax,0:kmax))
    allocate(vname(nbrvar))
    !
    var(0:imax,0:jmax,0:kmax,1)=real(var1(0:imax,0:jmax,0:kmax))
    var(0:imax,0:jmax,0:kmax,2)=real(var2(0:imax,0:jmax,0:kmax))
    var(0:imax,0:jmax,0:kmax,3)=real(var3(0:imax,0:jmax,0:kmax))
    var(0:imax,0:jmax,0:kmax,4)=real(var4(0:imax,0:jmax,0:kmax))
    var(0:imax,0:jmax,0:kmax,5)=real(var5(0:imax,0:jmax,0:kmax))
    var(0:imax,0:jmax,0:kmax,6)=real(var6(0:imax,0:jmax,0:kmax))
    !
    points(1,0:imax,0:jmax,0:kmax)=var(:,:,:,1)
    points(2,0:imax,0:jmax,0:kmax)=var(:,:,:,2)
    points(3,0:imax,0:jmax,0:kmax)=var(:,:,:,3)
    !
    vname(1)=var1name
    vname(2)=var2name
    vname(3)=var3name
    vname(4)=var4name
    vname(5)=var5name
    vname(6)=var6name
    !
    write(charxmax,'(I4)') imax
    write(charymax,'(I4)') jmax
    write(charzmax,'(I4)') kmax
    ! print*,charxmin
    ! stop

    ! open(18,file=filename)
    ! write(18,*)'<VTKFile type="StructuredGrid">'
    ! write(18,*)'  <StructuredGrid WholeExtent="' &
    !               //'0'//trim(charxmax)//' 0' &
    !               //trim(charymax)//' 0'//trim(charzmax)//'">'
    ! write(18,*)'    <Piece Extent="' &
    !               //'0'//trim(charxmax)//' 0' &
    !               //trim(charymax)//' 0'//trim(charzmax)//'">'
    ! write(18,*)'      <PointData>'
    ! write(18,*)'        <DataArray'
    ! write(18,*)'          type="Float32"'
    ! write(18,*)'          Name="'//trim(vname(4))//'"'
    ! write(18,*)'          format="'//'ascii'//'"'
    ! write(18,*)'        >'
    ! write(18,'(9(1PE12.4))')(var(:,:,:,4))
    ! write(18,*)'        </DataArray>'
    ! write(18,*)'      </PointData>'
    ! write(18,*)'      <CellData>'
    ! write(18,*)'      </CellData>'
    ! write(18,*)'      <Points>'
    ! write(18,*)'        <DataArray NumberOfComponents="3"'
    ! write(18,*)'          type="Float32"'
    ! write(18,*)'          format="'//'ascii'//'"'
    ! write(18,*)'        >'
    ! write(18,'(9(1PE12.4))')(points(:,:,:,:))
    ! write(18,*)'        </DataArray>'
    ! write(18,*)'      </Points>'
    ! write(18,*)'    </Piece>'
    ! write(18,*)'  </StructuredGrid>'
    ! write(18,*)'</VTKFile>'
    !
    ! Layout for the Appended Binary Data
    lf = char(10)
    inquire(iolength=real_bytesize) float32
    nbt_scal   = (imax+1) * (jmax+1) * (kmax+1) * real_bytesize
    nbt_vec    = 3 * (imax+1) * (jmax+1) * (kmax+1) * real_bytesize
    inquire(iolength=int_bytesize) int32
    ! nbt_ien    = 8 * 8  * int_bytesize
    ! nbt_offset = 8  * int_bytesize
    ! nbt_etype  = 8  * int_bytesize
    !
    open(18,file=filename,access='stream',convert='LITTLE_ENDIAN')
    write(18)'<VTKFile type="StructuredGrid" byte_order="LittleEndian">'//lf
    write(18)'  <StructuredGrid WholeExtent="' &
    //'0'//trim(charxmax)//' 0' &
    //trim(charymax)//' 0'//trim(charzmax)//'">'//lf
    write(18)'    <Piece Extent="' &
    //'0'//trim(charxmax)//' 0' &
    //trim(charymax)//' 0'//trim(charzmax)//'">'//lf
    write(18)'      <PointData>'//lf
    !
    ioff=0
    do jvar=4,nbrvar
      write(*,*) '** writing '//trim(vname(jvar))
      if(jvar>4)ioff=ioff + int_bytesize + nbt_scal
      write(offset(1:16),'(I16)') ioff
      write(18)'        <DataArray'//lf
      write(18)'          type="Float32"'//lf
      write(18)'          Name="'//trim(vname(jvar))//'"'//lf
      write(18)'          format="appended"'//lf
      write(18)'          offset="'//offset//'"'//lf
      write(18)'        />'//lf
    enddo
    !
    write(18)'      </PointData>'//lf
    write(18)'      <CellData>'//lf
    write(18)'      </CellData>'//lf 
    !
    write(18)'      <Points>'//lf


    ioff=ioff + int_bytesize + nbt_scal
    write(offset(1:16),'(I16)') ioff
    !
    write(18)'        <DataArray NumberOfComponents="3"'//lf
    write(18)'          type="Float32"'//lf
    write(18)'          format="appended"'//lf
    write(18)'          offset="'//offset//'"'//lf
    write(18)'        />'//lf
    write(18)'      </Points>'//lf
    write(18)'    </Piece>'//lf
    write(18)'  </StructuredGrid>'//lf
    write(18)'  <AppendedData encoding="raw">'//lf
    write(18)'_'
    !
    write(*,*) '** appending data...'
    do jvar=4,nbrvar
    write(18)nbt_scal,var(:,:,:,jvar)
    enddo 
    !
    write(18)nbt_vec,points(:,:,:,:)
    !
    write(18) lf//'  </AppendedData>'//lf
    write(18)'</VTKFile>'//lf
    !
    close(18)
    !
    print*,' << ',filename
    !
    deallocate(var)
  end subroutine writeprvbin3d6var
  !
  subroutine writeprvbin3d9var(filename,var1,var1name,var2,var2name,   &
                                        var3,var3name,var4,var4name,   &
                                        var5,var5name,var6,var6name,   &
                                        var7,var7name,var8,var8name,   &
                                        var9,var9name,imax,jmax,kmax)
    !
    character(len=*),intent(in) :: filename
    integer,intent(in) :: imax,jmax,kmax
    real(8),intent(in) :: var1(0:imax,0:jmax,0:kmax),                 &
                          var2(0:imax,0:jmax,0:kmax),                 &
                          var3(0:imax,0:jmax,0:kmax),                 &
                          var4(0:imax,0:jmax,0:kmax),                 &
                          var5(0:imax,0:jmax,0:kmax),                 &
                          var6(0:imax,0:jmax,0:kmax),                 &
                          var7(0:imax,0:jmax,0:kmax),                 &
                          var8(0:imax,0:jmax,0:kmax),                 &
                          var9(0:imax,0:jmax,0:kmax)
    character(len=*),intent(in) :: var1name,var2name,var3name,var4name,&
                                    var5name,var6name,var7name,var8name,&
                                    var9name
    real(8) :: solutiontime1
    integer :: zonenumber1
    character(256) :: title1
    !
    integer :: int32,unitf,nbrvar,n,jvar
    integer :: int_bytesize,real_bytesize,nbt_scal,nbt_vec,nbt_xyz,nbt_ien,nbt_offset,nbt_etype
    integer(16) :: ioff
    ! ip : le point actuel
    !
    real(4) :: float32
    real(8) :: float64
    !
    real(4),allocatable,dimension(:,:,:,:) :: var,points
    character(256),allocatable,dimension(:) :: vname
    !
    character(40) :: zonename1
    character(256) :: ligne
    character(15) :: charxmin,charxmax,charymin,charymax,charzmin,charzmax
    character :: lf*1,str_npoint*15,str_ncell*15,offset*16
    !        
    nbrvar=9
    !
    allocate(var(0:imax,0:jmax,0:kmax,nbrvar),points(3,0:imax,0:jmax,0:kmax))
    allocate(vname(nbrvar))
    !
    var(0:imax,0:jmax,0:kmax,1)=real(var1(0:imax,0:jmax,0:kmax))
    var(0:imax,0:jmax,0:kmax,2)=real(var2(0:imax,0:jmax,0:kmax))
    var(0:imax,0:jmax,0:kmax,3)=real(var3(0:imax,0:jmax,0:kmax))
    var(0:imax,0:jmax,0:kmax,4)=real(var4(0:imax,0:jmax,0:kmax))
    var(0:imax,0:jmax,0:kmax,5)=real(var5(0:imax,0:jmax,0:kmax))
    var(0:imax,0:jmax,0:kmax,6)=real(var6(0:imax,0:jmax,0:kmax))
    var(0:imax,0:jmax,0:kmax,7)=real(var7(0:imax,0:jmax,0:kmax))
    var(0:imax,0:jmax,0:kmax,8)=real(var8(0:imax,0:jmax,0:kmax))
    var(0:imax,0:jmax,0:kmax,9)=real(var9(0:imax,0:jmax,0:kmax))
    !
    points(1,0:imax,0:jmax,0:kmax)=var(:,:,:,1)
    points(2,0:imax,0:jmax,0:kmax)=var(:,:,:,2)
    points(3,0:imax,0:jmax,0:kmax)=var(:,:,:,3)
    !
      vname(1)=var1name
      vname(2)=var2name
      vname(3)=var3name
      vname(4)=var4name
      vname(5)=var5name
      vname(6)=var6name
      vname(7)=var7name
      vname(8)=var8name
      vname(9)=var9name
      !
      !
      ! if(tecinfout) then
      !   write(*,'(1x,A31,1X,A50)')' ** Title of the tecplot file: ',title1
      !   write(*,'(1x,A19,1X,F8.3)')' ** Solution time: ',solutiontime1
      !   write(*,'(1x,A33,1X,I2)')' ** Number of variable to write: ',nbrvar
      ! endif
      !
    write(charxmax,'(I4)') imax
    write(charymax,'(I4)') jmax
    write(charzmax,'(I4)') kmax
    ! print*,charxmin
    ! stop

      ! open(18,file=filename)
    ! write(18,*)'<VTKFile type="StructuredGrid">'
    ! write(18,*)'  <StructuredGrid WholeExtent="' &
    !               //'0'//trim(charxmax)//' 0' &
    !               //trim(charymax)//' 0'//trim(charzmax)//'">'
    ! write(18,*)'    <Piece Extent="' &
    !               //'0'//trim(charxmax)//' 0' &
    !               //trim(charymax)//' 0'//trim(charzmax)//'">'
    ! write(18,*)'      <PointData>'
    ! write(18,*)'        <DataArray'
    ! write(18,*)'          type="Float32"'
    ! write(18,*)'          Name="Density"'
    ! write(18,*)'          format="'//'ascii'//'"'
    ! write(18,*)'        >'
    ! write(18,'(9(1PE12.4))')(var(:,:,:,5))
    ! write(18,*)'        </DataArray>'
    ! write(18,*)'      </PointData>'
    ! write(18,*)'      <CellData>'
    ! write(18,*)'      </CellData>'
    ! write(18,*)'      <Points>'
    ! write(18,*)'        <DataArray NumberOfComponents="3"'
    ! write(18,*)'          type="Float32"'
    ! write(18,*)'          format="'//'ascii'//'"'
    ! write(18,*)'        >'
    ! write(18,'(9(1PE12.4))')(points(:,:,:,:))
    ! write(18,*)'        </DataArray>'
    ! write(18,*)'      </Points>'
    ! write(18,*)'    </Piece>'
    ! write(18,*)'  </StructuredGrid>'
    ! write(18,*)'</VTKFile>'
    !
      close(18)
      !
      ! Layout for the Appended Binary Data
    lf = char(10)
    inquire(iolength=real_bytesize) float32
    nbt_scal   = (imax+1) * (jmax+1) * (kmax+1) * real_bytesize
    nbt_vec    = 3 * (imax+1) * (jmax+1) * (kmax+1) * real_bytesize
    inquire(iolength=int_bytesize) int32
    ! nbt_ien    = 8 * 8  * int_bytesize
    ! nbt_offset = 8  * int_bytesize
    ! nbt_etype  = 8  * int_bytesize
    !
    open(18,file=filename,access='stream',convert='LITTLE_ENDIAN')
    write(18)'<VTKFile type="StructuredGrid" byte_order="LittleEndian">'//lf
    write(18)'  <StructuredGrid WholeExtent="' &
                  //'0'//trim(charxmax)//' 0' &
                  //trim(charymax)//' 0'//trim(charzmax)//'">'//lf
    write(18)'    <Piece Extent="' &
                  //'0'//trim(charxmax)//' 0' &
                  //trim(charymax)//' 0'//trim(charzmax)//'">'//lf
    write(18)'      <PointData>'//lf
    !
    ioff=0
    do jvar=4,nbrvar
      if(jvar>4)ioff=ioff + int_bytesize + nbt_scal
      write(offset(1:16),'(I16)') ioff
      write(18)'        <DataArray'//lf
      write(18)'          type="Float32"'//lf
      write(18)'          Name="'//trim(vname(jvar))//'"'//lf
      write(18)'          format="appended"'//lf
      write(18)'          offset="'//offset//'"'//lf
      write(18)'        />'//lf
    enddo
    write(18)'      </PointData>'//lf
    write(18)'      <CellData>'//lf
    write(18)'      </CellData>'//lf 
    !
    write(18)'      <Points>'//lf


    ioff=ioff + int_bytesize + nbt_scal
    write(offset(1:16),'(I16)') ioff
    !
    write(18)'        <DataArray NumberOfComponents="3"'//lf
    write(18)'          type="Float32"'//lf
    write(18)'          format="appended"'//lf
    write(18)'          offset="'//offset//'"'//lf
    write(18)'        />'//lf
    write(18)'      </Points>'//lf
    write(18)'    </Piece>'//lf
    write(18)'  </StructuredGrid>'//lf
    write(18)'  <AppendedData encoding="raw">'//lf
    write(18)'_'
    !
    do jvar=4,nbrvar
      write(18)nbt_scal,var(:,:,:,jvar)
    enddo 
    !
    write(18)nbt_vec,points(:,:,:,:)
    !
    write(18) lf//'  </AppendedData>'//lf
    write(18)'</VTKFile>'//lf
    !
      close(18)
    !
      print*,' << ',filename
      !
      deallocate(var)
  end subroutine writeprvbin3d9var
    !
  subroutine writeprvbin2d7var(filename,var1,var1name,var2,var2name,   &
                                        var3,var3name,var4,var4name,   &
                                        var5,var5name,var6,var6name,   &
                                        var7,var7name,imax,jmax)
    !
    character(len=*),intent(in) :: filename
    integer,intent(in) :: imax,jmax
    real(8),intent(in) :: var1(0:imax,0:jmax),                 &
                          var2(0:imax,0:jmax),                 &
                          var3(0:imax,0:jmax),                 &
                          var4(0:imax,0:jmax),                 &
                          var5(0:imax,0:jmax),                 &
                          var6(0:imax,0:jmax),                 &
                          var7(0:imax,0:jmax)
    character(len=*),intent(in) :: var1name,var2name,var3name, &
                                   var4name,var5name,var6name, &
                                   var7name
    !
    integer :: int32,nbrvar,jvar,kmax
    integer :: int_bytesize,real_bytesize,nbt_scal,nbt_vec!,nbt_xyz,nbt_ien,nbt_offset,nbt_etype
    integer(16) :: ioff
    ! ip : le point actuel
    !
    real(4) :: float32
    !
    real(4),allocatable,dimension(:,:,:,:) :: var,points
    character(256),allocatable,dimension(:) :: vname
    !
    character(15) :: charxmax,charymax,charzmax
    character :: lf*1,offset*16
    !        
    nbrvar=7
    kmax=0
    !
    allocate(var(0:imax,0:jmax,0:kmax,nbrvar),points(3,0:imax,0:jmax,0:kmax))
    allocate(vname(nbrvar))
    !
    var(0:imax,0:jmax,0,1)=real(var1(0:imax,0:jmax))
    var(0:imax,0:jmax,0,2)=real(var2(0:imax,0:jmax))
    var(0:imax,0:jmax,0,3)=real(var3(0:imax,0:jmax))
    var(0:imax,0:jmax,0,4)=real(var4(0:imax,0:jmax))
    var(0:imax,0:jmax,0,5)=real(var5(0:imax,0:jmax))
    var(0:imax,0:jmax,0,6)=real(var6(0:imax,0:jmax))
    var(0:imax,0:jmax,0,7)=real(var7(0:imax,0:jmax))
    !
    points(1,0:imax,0:jmax,0:kmax)=var(:,:,:,1)
    points(2,0:imax,0:jmax,0:kmax)=var(:,:,:,2)
    points(3,0:imax,0:jmax,0:kmax)=var(:,:,:,3)
    !
    vname(1)=var1name
    vname(2)=var2name
    vname(3)=var3name
    vname(4)=var4name
    vname(5)=var5name
    vname(6)=var6name
    vname(7)=var7name
    !
    write(charxmax,'(I4)') imax
    write(charymax,'(I4)') jmax
    write(charzmax,'(I4)') kmax
    !
    ! Layout for the Appended Binary Data
    lf = char(10)
    inquire(iolength=real_bytesize) float32
    nbt_scal   = (imax+1) * (jmax+1) * (kmax+1) * real_bytesize
    nbt_vec    = 3 * (imax+1) * (jmax+1) * (kmax+1) * real_bytesize
    inquire(iolength=int_bytesize) int32
    ! nbt_ien    = 8 * 8  * int_bytesize
    ! nbt_offset = 8  * int_bytesize
    ! nbt_etype  = 8  * int_bytesize
    !
    open(18,file=filename,access='stream',convert='LITTLE_ENDIAN')
    write(18)'<VTKFile type="StructuredGrid" byte_order="LittleEndian">'//lf
    write(18)'  <StructuredGrid WholeExtent="' &
    //'0'//trim(charxmax)//' 0' &
    //trim(charymax)//' 0'//trim(charzmax)//'">'//lf
    write(18)'    <Piece Extent="' &
    //'0'//trim(charxmax)//' 0' &
    //trim(charymax)//' 0'//trim(charzmax)//'">'//lf
    write(18)'      <PointData>'//lf
    !
    ioff=0
    do jvar=4,nbrvar
      write(*,*) '** writing '//trim(vname(jvar))
      if(jvar>4)ioff=ioff + int_bytesize + nbt_scal
      write(offset(1:16),'(I16)') ioff
      write(18)'        <DataArray'//lf
      write(18)'          type="Float32"'//lf
      write(18)'          Name="'//trim(vname(jvar))//'"'//lf
      write(18)'          format="appended"'//lf
      write(18)'          offset="'//offset//'"'//lf
      write(18)'        />'//lf
    enddo
    !
    write(18)'      </PointData>'//lf
    write(18)'      <CellData>'//lf
    write(18)'      </CellData>'//lf 
    !
    write(18)'      <Points>'//lf


    ioff=ioff + int_bytesize + nbt_scal
    write(offset(1:16),'(I16)') ioff
    !
    write(18)'        <DataArray NumberOfComponents="3"'//lf
    write(18)'          type="Float32"'//lf
    write(18)'          format="appended"'//lf
    write(18)'          offset="'//offset//'"'//lf
    write(18)'        />'//lf
    write(18)'      </Points>'//lf
    write(18)'    </Piece>'//lf
    write(18)'  </StructuredGrid>'//lf
    write(18)'  <AppendedData encoding="raw">'//lf
    write(18)'_'
    !
    write(*,*) '** appending data...'
    do jvar=4,nbrvar
    write(18)nbt_scal,var(:,:,:,jvar)
    enddo 
    !
    write(18)nbt_vec,points(:,:,:,:)
    !
    write(18) lf//'  </AppendedData>'//lf
    write(18)'</VTKFile>'//lf
    !
    close(18)
    !
    print*,' << ',filename
    !
    deallocate(var)
  end subroutine writeprvbin2d7var
  !
end module WriteVTK
  !+---------------------------------------------------------------------+
  !|The end of the module WriteTec                                       |
  !+---------------------------------------------------------------------+
  