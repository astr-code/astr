!+---------------------------------------------------------------------+
!| This module is to define constants.                                 |
!+---------------------------------------------------------------------+
!| CHANGE RECORD                                                       |
!| -------------                                                       |
!| 06-02-2021  | Created by J. Fang                                    |
!+---------------------------------------------------------------------+
module constdef
   
   use, intrinsic :: iso_fortran_env, only: real32, real64

   implicit none
   
   ! integer, parameter, public :: rtype = KIND(0._real32)

   integer, parameter, public :: rtype = KIND(0._real64)

   real(rtype), parameter :: pi = 4._rtype*atan(1.0_rtype), &
                       num1d35  = 1._rtype /35._rtype,  num1d3   = 1._rtype /3._rtype,   &
                       num2d3   = 2._rtype /3._rtype,   num1d24  = 1._rtype /24._rtype,  &
                       num4d3   = 4._rtype /3._rtype,   num1d6   = 1._rtype /6._rtype,   &
                       num1d12  = 1._rtype /12._rtype,  num7d12  = 7._rtype /12._rtype,  &
                       num7d9   = 7._rtype /9._rtype,   num1d36  = 1._rtype /36._rtype,  &
                       num1d60  = 1._rtype /60._rtype,  num65d3  = 65._rtype/3._rtype,   &
                       num20d3  = 20._rtype/3._rtype,   num1d11  = 1._rtype /11._rtype,  &
                       num25d12 = 25._rtype/12._rtype,  num11d6  = 11._rtype/6._rtype,   &
                       num1d840 = 1._rtype /840._rtype, num1d280 = 1._rtype /280._rtype, &
                       num4d105 = 4._rtype /105._rtype, num13d60 = 13._rtype/60._rtype,  &
                       num1d30  = 1._rtype /30._rtype,  num47d60 = 47._rtype/60._rtype,  &
                       num5d6   = 5._rtype /6._rtype,   num1d18  = 1._rtype /18._rtype,  &
                       num19d18 = 19._rtype/18._rtype,  num5d9   = 5._rtype /9._rtype,   &
                       num3d44  = 3._rtype/44._rtype,   num12d11 = 12._rtype/11._rtype,  &
                       num2d11  = 2._rtype/11._rtype
   !+-------------------------------------------------------------------+
   !| constants.                                                        |
   !+-------------------------------------------------------------------+
   !
end module constdef
!+---------------------------------------------------------------------+
!| The end of the module constdef.                                     |
!+---------------------------------------------------------------------+