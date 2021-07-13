!+---------------------------------------------------------------------+
!| The module declares new types.                                      |
!+---------------------------------------------------------------------+
!| ==============                                                      |
!| CHANGE RECORD                                                       |
!| -------------                                                       |
!| 02-Jul-2021  | Created by J. Fang STFC                              |
!+---------------------------------------------------------------------+
module commtype
  !
  implicit none
  !
  type :: triangle
    real(8) :: a(3),b(3),c(3),normdir(3),area,cen(3)
    !+-------------------+---------------------------------------------+
    !|            a,b,c  | coordinates of the 3 vertex of the triangle.|
    !|           normdir | normal direction of the face.               |
    !|              area | area of of the triangle.                    |
    !+-------------------+---------------------------------------------+
  end type triangle
  !
  type :: lsegment
    real(8) :: a(2),b(2),normdir(2),length,cen(2)
    !+-------------------+---------------------------------------------+
    !|            a,b,c  | coordinates of the 3 vertex of the triangle.|
    !|           normdir | normal direction of the face.               |
    !|              area | area of of the triangle.                    |
    !+-------------------+---------------------------------------------+
  end type lsegment
  !
  type :: solid
    !
    character(len=32) :: name
    real(8) :: xmin(3),xmax(3),xref(3),xcen(3)
    integer :: num_face,num_edge
    type(triangle),allocatable :: face(:)
    type(lsegment),allocatable :: edge(:)
    !
    contains
    !
    procedure :: alloface
    procedure :: alloedge
    !
  end type solid
  !
  !+---------------------+---------------------------------------------+
  !|           solidbody | a type of describing immersed solid body    |
  !|            triangle | a type of describing a  triangle            |
  !+---------------------+---------------------------------------------+
  !
  contains 
  !
  !+-------------------------------------------------------------------+
  !| This subroutine initialise faces of a solid.                      |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 02-Jul-2021  | Created by J. Fang STFC                            |
  !+-------------------------------------------------------------------+
  subroutine alloface(asolid)
    !
    class(solid),target :: asolid
    !
    allocate(asolid%face(asolid%num_face))
    !
  end subroutine alloface
  !+-------------------------------------------------------------------+
  !| The end of the subroutine alloface                                |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine initialise edges of a solid.                      |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 07-Jul-2021  | Created by J. Fang STFC                            |
  !+-------------------------------------------------------------------+
  subroutine alloedge(asolid)
    !
    class(solid),target :: asolid
    !
    allocate(asolid%edge(asolid%num_edge))
    !
  end subroutine alloedge
  !+-------------------------------------------------------------------+
  !| The end of the subroutine alloedge                                |
  !+-------------------------------------------------------------------+
  !
end module commtype
!+---------------------------------------------------------------------+
!| The end of the module commtype.                                     |
!+---------------------------------------------------------------------+
