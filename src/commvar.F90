!+---------------------------------------------------------------------+
!| This module is to define common variables.                          |
!+---------------------------------------------------------------------+
!| CHANGE RECORD                                                       |
!| -------------                                                       |
!| 06-02-2021  | Created by J. Fang                                    |
!+---------------------------------------------------------------------+
module commvar
  !
  implicit none
  !
  integer :: ia,ja,ka,im,jm,km,is,ie,js,je,ks,ke
  integer :: nh
  !+---------------------+---------------------------------------------+
  !|            ia,ja,ka | the dimension of the entire domain          | 
  !|            im,jm,km | the dimension of the local domain           | 
  !|   is,ie,js,je,ks,ke | start and end nodes number.                 |
  !|                  nh | level of halo nodes.                        |
  !+---------------------+---------------------------------------------+
  logical :: lihomo,ljhomo,lkhomo
  !+---------------------+---------------------------------------------+
  !| lihomo,ljhomo,lkhomo| to define homogeneous direction.            |
  !+---------------------+---------------------------------------------+
  !
  integer :: npdci,npdcj,npdck
  !+---------------------+---------------------------------------------+
  !|               npdci | to control scheme at boundary.              |
  !+---------------------+---------------------------------------------+
  !
  parameter(nh=4)
  !
end module commvar
!+---------------------------------------------------------------------+
!| The end of the module commvar.                                      |
!+---------------------------------------------------------------------+