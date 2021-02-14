!+---------------------------------------------------------------------+
!| This module contains subroutines of applying boundary conditions.   |
!| ==============                                                      |
!| CHANGE RECORD                                                       |
!| -------------                                                       |
!| 13-02-2021: Created by J. Fang @ Warrington                         |
!+---------------------------------------------------------------------+
module bc
  !
  use constdef
  use parallel,only: lio,mpistop,mpirank,mpirankname,irk,jrk,krk,      &
                     irkm,jrkm,krkm
  use commvar, only: hm,im,jm,km,uinf,vinf,pinf,roinf,ndims,num_species
  use tecio
  !
  implicit none
  !
  contains
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to apply bounday conditions.                   |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 13-02-2021: Created by J. Fang @ Warrington                       |
  !+-------------------------------------------------------------------+
  subroutine boucon
    !
    use commvar, only : bctype,twall
    !
    ! local data
    integer :: n
    !
    do n=1,6
      !
      if(bctype(n)==41) then
        call noslip(n,twall(n))
      endif
      !
    enddo
    !
  end subroutine boucon
  !+-------------------------------------------------------------------+
  !| The end of the subroutine boucon.                                 |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to apply nonslip bc.                           |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 13-02-2021: Created by J. Fang @ Warrington                       |
  !+-------------------------------------------------------------------+
  subroutine noslip(ndir,tw)
    !
    use commarray, only : prs,vel,tmp,rho,spc,q
    use fludyna,   only : thermal,fvar2q
    !
    ! arguments
    integer,intent(in) :: ndir
    real(8),intent(in) :: tw
    !
    ! local data
    integer :: i,j,k,l,jspec
    real(8) :: pe
    !
    if(ndir==3) then
      !
      if(jrk==0) then
        !
        j=0
        do k=0,km
        do i=0,im
          !
          pe=num1d3*(4.d0*prs(i,1,k)-prs(i,2,k))
          !
          vel(i,j,k,1)=0.d0
          vel(i,j,k,2)=0.d0
          vel(i,j,k,3)=0.d0
          prs(i,j,k)  =pe
          tmp(i,j,k)  =tw
          !
          do jspec=1,num_species
            spc(i,j,k,jspec)=num1d3*(4.d0*spc(i,1,k,jspec)-spc(i,2,k,jspec))
          enddo
          !
          rho(i,j,k)  =thermal(pressure=prs(i,j,k),temperature=tmp(i,j,k))
          !
          call fvar2q(      q=  q(i,j,k,:),                            &
                      density=rho(i,j,k),                              &
                     velocity=vel(i,j,k,:),                            &
                     pressure=prs(i,j,k),                              &
                      species=spc(i,j,k,:)                             )
          !
          do l=1,hm
            q(i,j-l,k,1)= q(i,j+l,k,1) ! rho   is even
            q(i,j-l,k,2)=-q(i,j+l,k,2) ! rho*u is odd 
            q(i,j-l,k,3)=-q(i,j+l,k,3) ! rho*v is odd 
            q(i,j-l,k,4)=-q(i,j+l,k,4) ! rho*w is odd
            q(i,j-l,k,5)= q(i,j+l,k,5) ! rho*E is even
            !
            do jspec=1,num_species
              q(i,j-l,k,5+jspec)= q(i,j+l,k,5+jspec) ! rho*Yj is even
            enddo
            !
          enddo
          !
        enddo
        enddo
        !
      endif
      !
    elseif(ndir==4) then
      !
      if(jrk==jrkm) then
        !
        j=jm
        do k=0,km
        do i=0,im
          pe=num1d3*(4.d0*prs(i,jm-1,k)-prs(i,jm-2,k))
          !
          vel(i,j,k,1)=0.d0
          vel(i,j,k,2)=0.d0
          vel(i,j,k,3)=0.d0
          prs(i,j,k)  =pe
          tmp(i,j,k)  =tw
          !
          if(num_species>0) then
            spc(i,j,k,1)=num1d3*(4.d0*spc(i,jm-1,k,1)-spc(i,jm-2,k,1))
          endif
          !
          rho(i,j,k)  =thermal(pressure=prs(i,j,k),temperature=tmp(i,j,k))
          !
          call fvar2q(      q=  q(i,j,k,:),                            &
                      density=rho(i,j,k),                              &
                     velocity=vel(i,j,k,:),                            &
                     pressure=prs(i,j,k),                              &
                      species=spc(i,j,k,:)                             )
          !
          do l=1,hm
            q(i,j+l,k,1)= q(i,j-l,k,1) ! rho   is even
            q(i,j+l,k,2)=-q(i,j-l,k,2) ! rho*u is odd 
            q(i,j+l,k,3)=-q(i,j-l,k,3) ! rho*v is odd 
            q(i,j+l,k,4)=-q(i,j-l,k,4) ! rho*w is odd
            q(i,j+l,k,5)= q(i,j-l,k,5) ! rho*E is even
            !
            do jspec=1,num_species
              q(i,j+l,k,5+jspec)= q(i,j-l,k,5+jspec) ! rho*Yj is even
            enddo
            !
          enddo
          !
        enddo
        enddo
        !
      endif
      !
    else
      stop ' !! ndir not defined @ noslip'
    endif
    !
  end subroutine noslip
  !+-------------------------------------------------------------------+
  !| The end of the subroutine noslip.                                 |
  !+-------------------------------------------------------------------+
    !
  !
end module bc
!+---------------------------------------------------------------------+
!| The end of the module bc.                                           |
!+---------------------------------------------------------------------+