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
    real(8) :: a(3),b(3),c(3),normdir(3),area
    !+-------------------+---------------------------------------------+
    !|            a,b,c  | coordinates of the 3 vertex of the triangle.|
    !|           normdir | normal direction of the face.               |
    !|              area | area of of the triangle.                    |
    !+-------------------+---------------------------------------------+
  end type triangle
  !
  type :: solid
    !
    character(len=32) :: name
    real(8) :: xmin(3),xmax(3),xref(3),xcen(3)
    integer :: num_face
    type(triangle),allocatable :: face(:)
    !
    contains
    !
    procedure :: alloface
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
end module commtype
!+---------------------------------------------------------------------+
!| The end of the module commtype.                                     |
!+---------------------------------------------------------------------+
