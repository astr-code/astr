!+---------------------------------------------------------------------+
!| This module contains subroutines of writing file to view via tecplot|
!| ==============                                                      |
!| CHANGE RECORD                                                       |
!| -------------                                                       |
!| 07-Feb-2020  | Created by J. Fang @ Warrington                      |
!+---------------------------------------------------------------------+
module tecio
  !
  implicit none
  !
  real(4) :: eohmarker         = 357.0,                                &
             zonemarker        = 299.0,                                &
             geometrymarker    = 399.0,                                &
             textmarker        = 499.0,                                &
             customlabelmarker = 599.0,                                &
             userrecmarker     = 699.0,                                &
             datasetauxmarker  = 799.0,                                &
             varauxmarker      = 899.0
  logical :: tecinfout=.false.
  !
  Interface tecbin
    !
    module procedure writetecbin3d3var
    module procedure writetecbin3d4var
    module procedure writetecbin3d6var
    !
  end Interface tecbin
  !
  contains
  !
  !+-------------------------------------------------------------------+
  !|This subroutine is used to write bin file for 3d tecplot field.    |
  !|ifort compiler only                                                |
  !+-------------------------------------------------------------------+
  subroutine writetecbin3d3var(filename,var1,var1name,var2,var2name,   &
                                                      var3,var3name)
    ! 
    character(len=*),intent(in) :: filename
    real(8),intent(in) :: var1(:,:,:),var2(:,:,:),var3(:,:,:)
    character(len=*),intent(in) :: var1name,var2name,var3name
    !
    ! local data
    !
    integer :: imax,jmax,kmax
    real(8) :: solutiontime1
    integer :: zonenumber1
    character(256) :: title1
    !
    integer :: int32,unitf,nbrvar,n
    ! ip : le point actuel
    !
    real(4) :: float32
    real(8) :: float64
    !
    real(4),allocatable,dimension(:,:,:,:) :: var
    character(256),allocatable,dimension(:) :: vname
    !
    character(40) :: zonename1
    character(256) :: ligne
    !
    solutiontime1=0.d0
    zonenumber1=1
    title1="Bin field for tecplot"
    !
    nbrvar=3
    !
    imax=size(var1,1)
    jmax=size(var1,2)
    kmax=size(var1,3)
    !
    allocate(var(1:imax,1:jmax,1:kmax,nbrvar))
    allocate(vname(nbrvar))
    !
    var(1:imax,1:jmax,1:kmax,1)=real(var1(1:imax,1:jmax,1:kmax),4)
    var(1:imax,1:jmax,1:kmax,2)=real(var2(1:imax,1:jmax,1:kmax),4)
    var(1:imax,1:jmax,1:kmax,3)=real(var3(1:imax,1:jmax,1:kmax),4)
    !
    vname(1)=var1name
    vname(2)=var2name
    vname(3)=var3name
    !
    unitf=18
    !
    open(18,file=filename,form='unformatted',access='stream')
    !
      !i. header section
    ! i. magic number, version number
    ! +------------+
    ! | "#!tdv112" |
    ! +------------+
    write(unitf)"#!TDV112"
    ! ii. integer value of 1
    ! +------------+
    ! | int32      |
    ! +------------+
    int32=1
    write(unitf)int32
    ! iii. title and variable names
    ! +------------+
    ! | int32      | filetype: 0=full, 1=grid, 2=solution
    ! +------------+
    int32=0
    write(unitf)int32
    ! +------------+
    ! | int32*n    | the title
    ! +------------+?
    call ecrirebin(unitf,title1)
    ! +------------+
    ! | int32      | number of variables in the datafile
    ! +------------+
    write(unitf)nbrvar
    ! +------------+
    ! | int32*n    | variable names
    ! +------------+
    do n=1,nbrvar
      call ecrirebin(unitf,vname(n))
        if(tecinfout) write(*,'(1x,A12,I2,A4,A10)')' ** Variable',n,' is ',vname(n)
    enddo
    ! iv. zones
    ! +------------+
    ! | float32    | zone marker. value = 299.0
    ! +------------+
    write(unitf)zonemarker
    ! +------------+
    ! | int32*n    | zone name
    ! +------------+ 
    Ligne=""
    write(Ligne,"(A,I3.3)")"Zone",zonenumber1
    call EcrireBin(UnitF,Ligne)
      if(tecinfout) write(*,'(1x,A15,1X,A10)')' ** Zone name: ',Ligne
    ! +------------+
    ! | int32      | parentzone
    ! +------------+
    int32=-1
    write(unitf)int32
    ! +------------+
    ! | int32      | strandid
    ! +------------+
    int32=-1
    write(unitf)int32
    ! +------------+
    ! | float64    | solution time
    ! +------------+
    write(unitf)solutiontime1
    ! +------------+
    ! | int32      | not used. set to -1
    ! +------------+
    int32=-1
    write(unitf)int32
    ! +------------+
    ! | int32      | zonetype 0=ordered,       1=felineseg,
    ! +------------+          2=fetriangle,    3=fequadrilateral,
    !                         4=fetetrahedron, 5=febrick,
    !                         6=fepolygon,     7=fepolyhedron
    int32=0
    write(unitf)int32
    ! +------------+
    ! | int32      | specify var location
    ! +------------+    0 = don't specify, 1 = specify
    int32=0
    write(unitf)int32
    ! +------------+
    ! | int32      | are raw local 1-to-1 face neighbors supplied?
    ! +------------+
    int32=0
    write(unitf)int32
    ! +------------+
    ! | int32      | number of miscellaneous user-defined face neighbor connections
    ! +------------+
    int32=0
    write(unitf)int32
    ! +------------+
    ! | int32*3    | imax,jmax,kmax
    ! +------------+
    write(unitf)imax
    write(unitf)jmax
    write(unitf)kmax
    ! +------------+
    ! | int32      | 1=auxiliary name/value pair to follow
    ! +------------+ 0=no more auxiliary name/value pairs
    int32=0
    write(unitf)int32
    ! +------------+
    ! | float32    | eohmarker, value = 357.0, end of header section
    ! +------------+
    write(unitf)eohmarker
    !ii. data section
    ! i. for both ordered and fe zones
    ! +------------+
    ! | float32    | zone marker value = 299.0
    ! +------------+
    write(unitf)zonemarker
    ! +------------+
    ! | int32*n    | variable data format, n=total number of vars
    ! +------------+     1=float,    2=double, 3=longint
    !                    4=shortint, 5=byte,   6=bit
    do n=1,nbrvar
      int32=1
      write(18)int32
    enddo
    ! +------------+
    ! | int32      | has passive variables: 0=no, 1=yes
    ! +------------+
    int32=0
    write(18)int32
    ! +------------+
    ! | int32      | has variable sharing: 0=no, 1=yes
    ! +------------+
    int32=0
    write(18)int32
    ! +------------+
    ! | int32      | zero based zone number to share connectivity list with (-1 = no sharing)
    ! +------------+
    int32=-1
    write(18)int32
    !
    do n=1,nbrvar
      !+------------+
      !| float64    | min value
      !+------------+
      float32=minval(var(:,:,:,n))
      float64=real(float32,8)
      write(18)float64
      !+------------+
      !| float64    | max value
      !+------------+
      float32=maxval(var(:,:,:,n))
      float64=real(float32,8)
      write(18)float64
    enddo
    ! +------------+
    ! | xxxxxxxxxx | zone data
    ! +------------+
    write(18)var
    !
    !write(18)((((var(n,i,j,k),n=1,nbrvar),i=0,im),j=0,jm),k=0,km)
    !do n=1,nbrvar
    !write(18)((((var(i,j,k,n),i=0,imax),j=0,jmax),k=0,kmax),n=1,nbrvar)
      
    !end do
      !
    close(18)
    !
    print*,' << ',filename
    !
    deallocate(var)
  end subroutine writetecbin3d3var
  !
  subroutine writetecbin3d4var(filename,var1,var1name,var2,var2name,   &
                                        var3,var3name,var4,var4name )
    ! 
    character(len=*),intent(in) :: filename
    real(8),dimension(:,:,:),intent(in) :: var1,var2,var3,var4
    character(len=*),intent(in) :: var1name,var2name,var3name,         &
                                   var4name
    !
    ! local data
    !
    integer :: imax,jmax,kmax
    real(8) :: solutiontime1
    integer :: zonenumber1
    character(256) :: title1
    !
    integer :: int32,unitf,nbrvar,n
    ! ip : le point actuel
    !
    real(4) :: float32
    real(8) :: float64
    !
    real(4),allocatable,dimension(:,:,:,:) :: var
    character(256),allocatable,dimension(:) :: vname
    !
    character(40) :: zonename1
    character(256) :: ligne
    !
    solutiontime1=0.d0
    zonenumber1=1
    title1="Bin field for tecplot"
    !
    nbrvar=4
    !
    imax=size(var1,1)
    jmax=size(var1,2)
    kmax=size(var1,3)
    !
    allocate(var(1:imax,1:jmax,1:kmax,nbrvar))
    allocate(vname(nbrvar))
    !
    var(1:imax,1:jmax,1:kmax,1)=real(var1(1:imax,1:jmax,1:kmax),4)
    var(1:imax,1:jmax,1:kmax,2)=real(var2(1:imax,1:jmax,1:kmax),4)
    var(1:imax,1:jmax,1:kmax,3)=real(var3(1:imax,1:jmax,1:kmax),4)
    var(1:imax,1:jmax,1:kmax,4)=real(var4(1:imax,1:jmax,1:kmax),4)
    !
    vname(1)=var1name
    vname(2)=var2name
    vname(3)=var3name
    vname(4)=var4name
    !
    unitf=18
    !
    open(18,file=filename,form='unformatted',access='stream')
    !
      !i. header section
    ! i. magic number, version number
    ! +------------+
    ! | "#!tdv112" |
    ! +------------+
    write(unitf)"#!TDV112"
    ! ii. integer value of 1
    ! +------------+
    ! | int32      |
    ! +------------+
    int32=1
    write(unitf)int32
    ! iii. title and variable names
    ! +------------+
    ! | int32      | filetype: 0=full, 1=grid, 2=solution
    ! +------------+
    int32=0
    write(unitf)int32
    ! +------------+
    ! | int32*n    | the title
    ! +------------+?
    call ecrirebin(unitf,title1)
    ! +------------+
    ! | int32      | number of variables in the datafile
    ! +------------+
    write(unitf)nbrvar
    ! +------------+
    ! | int32*n    | variable names
    ! +------------+
    do n=1,nbrvar
      call ecrirebin(unitf,vname(n))
        if(tecinfout) write(*,'(1x,A12,I2,A4,A10)')' ** Variable',n,' is ',vname(n)
    enddo
    ! iv. zones
    ! +------------+
    ! | float32    | zone marker. value = 299.0
    ! +------------+
    write(unitf)zonemarker
    ! +------------+
    ! | int32*n    | zone name
    ! +------------+ 
    Ligne=""
    write(Ligne,"(A,I3.3)")"Zone",zonenumber1
    call EcrireBin(UnitF,Ligne)
      if(tecinfout) write(*,'(1x,A15,1X,A10)')' ** Zone name: ',Ligne
    ! +------------+
    ! | int32      | parentzone
    ! +------------+
    int32=-1
    write(unitf)int32
    ! +------------+
    ! | int32      | strandid
    ! +------------+
    int32=-1
    write(unitf)int32
    ! +------------+
    ! | float64    | solution time
    ! +------------+
    write(unitf)solutiontime1
    ! +------------+
    ! | int32      | not used. set to -1
    ! +------------+
    int32=-1
    write(unitf)int32
    ! +------------+
    ! | int32      | zonetype 0=ordered,       1=felineseg,
    ! +------------+          2=fetriangle,    3=fequadrilateral,
    !                         4=fetetrahedron, 5=febrick,
    !                         6=fepolygon,     7=fepolyhedron
    int32=0
    write(unitf)int32
    ! +------------+
    ! | int32      | specify var location
    ! +------------+    0 = don't specify, 1 = specify
    int32=0
    write(unitf)int32
    ! +------------+
    ! | int32      | are raw local 1-to-1 face neighbors supplied?
    ! +------------+
    int32=0
    write(unitf)int32
    ! +------------+
    ! | int32      | number of miscellaneous user-defined face neighbor connections
    ! +------------+
    int32=0
    write(unitf)int32
    ! +------------+
    ! | int32*3    | imax,jmax,kmax
    ! +------------+
    write(unitf)imax
    write(unitf)jmax
    write(unitf)kmax
    ! +------------+
    ! | int32      | 1=auxiliary name/value pair to follow
    ! +------------+ 0=no more auxiliary name/value pairs
    int32=0
    write(unitf)int32
    ! +------------+
    ! | float32    | eohmarker, value = 357.0, end of header section
    ! +------------+
    write(unitf)eohmarker
    !ii. data section
    ! i. for both ordered and fe zones
    ! +------------+
    ! | float32    | zone marker value = 299.0
    ! +------------+
    write(unitf)zonemarker
    ! +------------+
    ! | int32*n    | variable data format, n=total number of vars
    ! +------------+     1=float,    2=double, 3=longint
    !                    4=shortint, 5=byte,   6=bit
    do n=1,nbrvar
      int32=1
      write(18)int32
    enddo
    ! +------------+
    ! | int32      | has passive variables: 0=no, 1=yes
    ! +------------+
    int32=0
    write(18)int32
    ! +------------+
    ! | int32      | has variable sharing: 0=no, 1=yes
    ! +------------+
    int32=0
    write(18)int32
    ! +------------+
    ! | int32      | zero based zone number to share connectivity list with (-1 = no sharing)
    ! +------------+
    int32=-1
    write(18)int32
    !
    do n=1,nbrvar
      !+------------+
      !| float64    | min value
      !+------------+
      float32=minval(var(:,:,:,n))
      float64=real(float32,8)
      write(18)float64
      !+------------+
      !| float64    | max value
      !+------------+
      float32=maxval(var(:,:,:,n))
      float64=real(float32,8)
      write(18)float64
    enddo
    ! +------------+
    ! | xxxxxxxxxx | zone data
    ! +------------+
    write(18)var
    !
    !write(18)((((var(n,i,j,k),n=1,nbrvar),i=0,im),j=0,jm),k=0,km)
    !do n=1,nbrvar
    !write(18)((((var(i,j,k,n),i=0,imax),j=0,jmax),k=0,kmax),n=1,nbrvar)
      
    !end do
      !
    close(18)
    !
    print*,' << ',filename
    !
    deallocate(var)
  end subroutine writetecbin3d4var
  !
  subroutine writetecbin3d6var(filename,var1,var1name,var2,var2name,   &
                                        var3,var3name,var4,var4name,   &
                                        var5,var5name,var6,var6name )
    ! 
    character(len=*),intent(in) :: filename
    real(8),dimension(:,:,:),intent(in) :: var1,var2,var3,var4,var5,var6
    character(len=*),intent(in) :: var1name,var2name,var3name,         &
                                   var4name,var5name,var6name
    !
    ! local data
    !
    integer :: imax,jmax,kmax
    real(8) :: solutiontime1
    integer :: zonenumber1
    character(256) :: title1
    !
    integer :: int32,unitf,nbrvar,n
    ! ip : le point actuel
    !
    real(4) :: float32
    real(8) :: float64
    !
    real(4),allocatable,dimension(:,:,:,:) :: var
    character(256),allocatable,dimension(:) :: vname
    !
    character(40) :: zonename1
    character(256) :: ligne
    !
    solutiontime1=0.d0
    zonenumber1=1
    title1="Bin field for tecplot"
    !
    nbrvar=6
    !
    imax=size(var1,1)
    jmax=size(var1,2)
    kmax=size(var1,3)
    !
    allocate(var(1:imax,1:jmax,1:kmax,nbrvar))
    allocate(vname(nbrvar))
    !
    var(1:imax,1:jmax,1:kmax,1)=real(var1(1:imax,1:jmax,1:kmax),4)
    var(1:imax,1:jmax,1:kmax,2)=real(var2(1:imax,1:jmax,1:kmax),4)
    var(1:imax,1:jmax,1:kmax,3)=real(var3(1:imax,1:jmax,1:kmax),4)
    var(1:imax,1:jmax,1:kmax,4)=real(var4(1:imax,1:jmax,1:kmax),4)
    var(1:imax,1:jmax,1:kmax,5)=real(var5(1:imax,1:jmax,1:kmax),4)
    var(1:imax,1:jmax,1:kmax,6)=real(var6(1:imax,1:jmax,1:kmax),4)
    !
    vname(1)=var1name
    vname(2)=var2name
    vname(3)=var3name
    vname(4)=var4name
    vname(5)=var5name
    vname(6)=var6name
    !
    unitf=18
    !
    open(18,file=filename,form='unformatted',access='stream')
    !
      !i. header section
    ! i. magic number, version number
    ! +------------+
    ! | "#!tdv112" |
    ! +------------+
    write(unitf)"#!TDV112"
    ! ii. integer value of 1
    ! +------------+
    ! | int32      |
    ! +------------+
    int32=1
    write(unitf)int32
    ! iii. title and variable names
    ! +------------+
    ! | int32      | filetype: 0=full, 1=grid, 2=solution
    ! +------------+
    int32=0
    write(unitf)int32
    ! +------------+
    ! | int32*n    | the title
    ! +------------+?
    call ecrirebin(unitf,title1)
    ! +------------+
    ! | int32      | number of variables in the datafile
    ! +------------+
    write(unitf)nbrvar
    ! +------------+
    ! | int32*n    | variable names
    ! +------------+
    do n=1,nbrvar
      call ecrirebin(unitf,vname(n))
        if(tecinfout) write(*,'(1x,A12,I2,A4,A10)')' ** Variable',n,' is ',vname(n)
    enddo
    ! iv. zones
    ! +------------+
    ! | float32    | zone marker. value = 299.0
    ! +------------+
    write(unitf)zonemarker
    ! +------------+
    ! | int32*n    | zone name
    ! +------------+ 
    Ligne=""
    write(Ligne,"(A,I3.3)")"Zone",zonenumber1
    call EcrireBin(UnitF,Ligne)
      if(tecinfout) write(*,'(1x,A15,1X,A10)')' ** Zone name: ',Ligne
    ! +------------+
    ! | int32      | parentzone
    ! +------------+
    int32=-1
    write(unitf)int32
    ! +------------+
    ! | int32      | strandid
    ! +------------+
    int32=-1
    write(unitf)int32
    ! +------------+
    ! | float64    | solution time
    ! +------------+
    write(unitf)solutiontime1
    ! +------------+
    ! | int32      | not used. set to -1
    ! +------------+
    int32=-1
    write(unitf)int32
    ! +------------+
    ! | int32      | zonetype 0=ordered,       1=felineseg,
    ! +------------+          2=fetriangle,    3=fequadrilateral,
    !                         4=fetetrahedron, 5=febrick,
    !                         6=fepolygon,     7=fepolyhedron
    int32=0
    write(unitf)int32
    ! +------------+
    ! | int32      | specify var location
    ! +------------+    0 = don't specify, 1 = specify
    int32=0
    write(unitf)int32
    ! +------------+
    ! | int32      | are raw local 1-to-1 face neighbors supplied?
    ! +------------+
    int32=0
    write(unitf)int32
    ! +------------+
    ! | int32      | number of miscellaneous user-defined face neighbor connections
    ! +------------+
    int32=0
    write(unitf)int32
    ! +------------+
    ! | int32*3    | imax,jmax,kmax
    ! +------------+
    write(unitf)imax
    write(unitf)jmax
    write(unitf)kmax
    ! +------------+
    ! | int32      | 1=auxiliary name/value pair to follow
    ! +------------+ 0=no more auxiliary name/value pairs
    int32=0
    write(unitf)int32
    ! +------------+
    ! | float32    | eohmarker, value = 357.0, end of header section
    ! +------------+
    write(unitf)eohmarker
    !ii. data section
    ! i. for both ordered and fe zones
    ! +------------+
    ! | float32    | zone marker value = 299.0
    ! +------------+
    write(unitf)zonemarker
    ! +------------+
    ! | int32*n    | variable data format, n=total number of vars
    ! +------------+     1=float,    2=double, 3=longint
    !                    4=shortint, 5=byte,   6=bit
    do n=1,nbrvar
      int32=1
      write(18)int32
    enddo
    ! +------------+
    ! | int32      | has passive variables: 0=no, 1=yes
    ! +------------+
    int32=0
    write(18)int32
    ! +------------+
    ! | int32      | has variable sharing: 0=no, 1=yes
    ! +------------+
    int32=0
    write(18)int32
    ! +------------+
    ! | int32      | zero based zone number to share connectivity list with (-1 = no sharing)
    ! +------------+
    int32=-1
    write(18)int32
    !
    do n=1,nbrvar
      !+------------+
      !| float64    | min value
      !+------------+
      float32=minval(var(:,:,:,n))
      float64=real(float32,8)
      write(18)float64
      !+------------+
      !| float64    | max value
      !+------------+
      float32=maxval(var(:,:,:,n))
      float64=real(float32,8)
      write(18)float64
    enddo
    ! +------------+
    ! | xxxxxxxxxx | zone data
    ! +------------+
    write(18)var
    !
    !write(18)((((var(n,i,j,k),n=1,nbrvar),i=0,im),j=0,jm),k=0,km)
    !do n=1,nbrvar
    !write(18)((((var(i,j,k,n),i=0,imax),j=0,jmax),k=0,kmax),n=1,nbrvar)
      
    !end do
      !
    close(18)
    !
    print*,' << ',filename
    !
    deallocate(var)
  end subroutine writetecbin3d6var
  !+-------------------------------------------------------------------+
  !| The end of the subroutine writetecbin3d.                          |
  !+-------------------------------------------------------------------+
  !
  subroutine ecrirebin(unitf,string)
    !
    integer unitf
    integer i
    character*80 string
    !
    do i=1,len_trim(string)
      write(unitf)ichar(string(i:i))
    enddo
    write(unitf)0
    !
  end subroutine ecrirebin
  !
end module tecio
!+---------------------------------------------------------------------+
!| The end of the module tecio.                                        |
!+---------------------------------------------------------------------+