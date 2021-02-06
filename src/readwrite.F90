!+---------------------------------------------------------------------+
!| This module contains subroutines of reading and writing files.      |
!| ==============                                                      |
!| CHANGE RECORD                                                       |
!| -------------                                                       |
!| 06-Oct-2018  | Created by J. Fang STFC Daresbury Laboratory         |
!+---------------------------------------------------------------------+
module readwrite
  !
  contains
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used to print welcome infomation.              |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 08-Oct-2018  | Created by J. Fang STFC Daresbury Laboratory       |
  !+-------------------------------------------------------------------+
  subroutine statement
    !
    write(*,*)
    write(*,*)'                ___   _____________            ___ '
    write(*,*)'               / _ | / __/_  __/ _ \   ____   / _\\'
    write(*,*)'              / __ |_\ \  / / / , _/  /___/  / , _/'
    write(*,*)'             /_/ |_/___/ /_/ /_/|_|         /_/|_| '
    write(*,*)
    write(*,*)

    print*,' +----------------------- Statement --------------------------+'
    print*,' |                                                            |'
    print*,' |                   Developed by Jian Fang                   |'
    print*,' |             <ASTR-R> Copyright Resvered <2020>             |'
    print*,' |                 ASTR code for reacting flow                |'
    print*,' |                                                            |'
    print*,' +------------------------------------------------------------+'
    print*,' |                                                            |'
    print*,' | Copyright 2020 Jian Fang                                   |'
    print*,' |                                                            |'
    print*,' | Licensed under the Apache License, Version 2.0             |'
    print*,' | you may not use this file except in compliance with the    |'
    print*,' | license. You may obtain a copy of the License at           |'
    print*,' |                                                            |'
    print*,' |        http://www.apache.org/licenses/LICENSE-2.0          |'
    print*,' |                                                            |'
    print*,' | Unless required by applicable law or agreed to in writing, |'
    print*,' | software distributed under the License is distributed on an|'
    print*,' | "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY     |'
    print*,' | KIND, either express or implied. See the License for the   |'
    print*,' | specific language governing permissions and limitations    |'
    print*,' | under the License.                                         |'
    print*,' |                                                            |'
    print*,' +--------------------- End of Statement ---------------------+'
    !
  end subroutine statement
  !+-------------------------------------------------------------------+
  !| The end of the subroutine statement.                              |
  !+-------------------------------------------------------------------+
  !
end module readwrite
!+---------------------------------------------------------------------+
!| The end of the module readwrite.                                    |
!+---------------------------------------------------------------------+
