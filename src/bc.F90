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
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to apply bounday conditions to geometrical     |
  !| variables.                                                        |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 14-02-2021: Created by J. Fang @ Warrington                       |
  !+-------------------------------------------------------------------+
  subroutine geombc
    !
    use commvar, only : bctype
    use commarray, only : jacob,dxi
    use commfunc,  only : gradextrp
    !
    ! local data
    integer :: i,j,k,l
    real(8) :: fex(4,10)
    !
    if(jrk==0) then
      !
      j=0
      !
      if(bctype(3)==41) then
        do k=0,km
        do i=0,im
          !
          fex(:,1) =gradextrp(qbou=dxi(i,0,k,1,1),q1st=dxi(i,1,k,1,1))
          fex(:,2) =gradextrp(qbou=dxi(i,0,k,1,2),q1st=dxi(i,1,k,1,2))
          fex(:,3) =gradextrp(qbou=dxi(i,0,k,1,3),q1st=dxi(i,1,k,1,3))
          fex(:,4) =gradextrp(qbou=dxi(i,0,k,2,1),q1st=dxi(i,1,k,2,1))
          fex(:,5) =gradextrp(qbou=dxi(i,0,k,2,2),q1st=dxi(i,1,k,2,2))
          fex(:,6) =gradextrp(qbou=dxi(i,0,k,2,3),q1st=dxi(i,1,k,2,3))
          fex(:,7) =gradextrp(qbou=dxi(i,0,k,3,1),q1st=dxi(i,1,k,3,1))
          fex(:,8) =gradextrp(qbou=dxi(i,0,k,3,2),q1st=dxi(i,1,k,3,2))
          fex(:,9) =gradextrp(qbou=dxi(i,0,k,3,3),q1st=dxi(i,1,k,3,3))
          fex(:,10)=gradextrp(qbou=jacob(i,0,k),  q1st=jacob(i,1,k))
          !
          do l=1,4
            dxi(i,-l,k,1,1)=fex(l,1) 
            dxi(i,-l,k,1,2)=fex(l,2) 
            dxi(i,-l,k,1,3)=fex(l,3) 
            dxi(i,-l,k,2,1)=fex(l,4) 
            dxi(i,-l,k,2,2)=fex(l,5) 
            dxi(i,-l,k,2,3)=fex(l,6) 
            dxi(i,-l,k,3,1)=fex(l,7) 
            dxi(i,-l,k,3,2)=fex(l,8) 
            dxi(i,-l,k,3,3)=fex(l,9) 
            jacob(i,-l,k)  =fex(l,10)
          enddo
          !
        enddo
        enddo
      endif
      !
    endif
    !
    if(jrk==jrkm) then
      !
      j=0
      !
      if(bctype(4)==41) then
        do k=0,km
        do i=0,im
          !
          dxi(i,jm+1:jm+4,k,1,1)=gradextrp(qbou=dxi(i,jm,k,1,1),q1st=dxi(i,jm-1,k,1,1))
          dxi(i,jm+1:jm+4,k,1,2)=gradextrp(qbou=dxi(i,jm,k,1,2),q1st=dxi(i,jm-1,k,1,2))
          dxi(i,jm+1:jm+4,k,1,3)=gradextrp(qbou=dxi(i,jm,k,1,3),q1st=dxi(i,jm-1,k,1,3))
          dxi(i,jm+1:jm+4,k,2,1)=gradextrp(qbou=dxi(i,jm,k,2,1),q1st=dxi(i,jm-1,k,2,1))
          dxi(i,jm+1:jm+4,k,2,2)=gradextrp(qbou=dxi(i,jm,k,2,2),q1st=dxi(i,jm-1,k,2,2))
          dxi(i,jm+1:jm+4,k,2,3)=gradextrp(qbou=dxi(i,jm,k,2,3),q1st=dxi(i,jm-1,k,2,3))
          dxi(i,jm+1:jm+4,k,3,1)=gradextrp(qbou=dxi(i,jm,k,3,1),q1st=dxi(i,jm-1,k,3,1))
          dxi(i,jm+1:jm+4,k,3,2)=gradextrp(qbou=dxi(i,jm,k,3,2),q1st=dxi(i,jm-1,k,3,2))
          dxi(i,jm+1:jm+4,k,3,3)=gradextrp(qbou=dxi(i,jm,k,3,3),q1st=dxi(i,jm-1,k,3,3))
          jacob(1,jm+1:jm+4,k)  =gradextrp(qbou=jacob(i,jm,k),  q1st=jacob(i,jm-1,k))
          !
        enddo
        enddo
      endif
      !
    endif
    !
  end subroutine geombc
  !
  subroutine xyzbc
    !
    use commvar,  only : bctype
    use commarray,only : x
    use commfunc, only : gradextrp
    !
    ! local data
    integer :: i,j,k,l
    real(8) :: fex(4,3)
    !
    if(jrk==0) then
      !
      j=0
      !
      if(bctype(3)==41) then
        do k=0,km
        do i=0,im
          !
          fex(:,1) =gradextrp(qbou=x(i,0,k,1),q1st=x(i,1,k,1))
          fex(:,2) =gradextrp(qbou=x(i,0,k,2),q1st=x(i,1,k,2))
          fex(:,3) =gradextrp(qbou=x(i,0,k,3),q1st=x(i,1,k,3))
          !
          do l=1,4
            x(i,-l,k,1)=fex(l,1) 
            x(i,-l,k,2)=fex(l,2) 
            x(i,-l,k,3)=fex(l,3) 
          enddo
          !
        enddo
        enddo
      endif
      !
    endif
    !
    if(jrk==jrkm) then
      !
      j=0
      !
      if(bctype(4)==41) then
        do k=0,km
        do i=0,im
          !
          x(i,jm+1:jm+4,k,1)=gradextrp(qbou=x(i,jm,k,1),q1st=x(i,jm-1,k,1))
          x(i,jm+1:jm+4,k,2)=gradextrp(qbou=x(i,jm,k,2),q1st=x(i,jm-1,k,2))
          x(i,jm+1:jm+4,k,3)=gradextrp(qbou=x(i,jm,k,3),q1st=x(i,jm-1,k,3))
          !
        enddo
        enddo
      endif
      !
    endif
    !
  end subroutine xyzbc
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
    use fludyna,   only : thermal,fvar2q,q2fvar
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
            call q2fvar(q=q(i,j-l,k,:),density=rho(i,j-l,k),           &
                                      velocity=vel(i,j-l,k,:),         &
                                      pressure=prs(i,j-l,k),           &
                                   temperature=tmp(i,j-l,k),           &
                                       species=spc(i,j-l,k,:))
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
          pe=num1d3*(4.d0*prs(i,j-1,k)-prs(i,j-2,k))
          !
          vel(i,j,k,1)=0.d0
          vel(i,j,k,2)=0.d0
          vel(i,j,k,3)=0.d0
          prs(i,j,k)  =pe
          tmp(i,j,k)  =tw
          !
          if(num_species>0) then
            spc(i,j,k,1)=num1d3*(4.d0*spc(i,j-1,k,1)-spc(i,j-2,k,1))
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
            call q2fvar(q=q(i,j+l,k,:),density=rho(i,j+l,k),           &
                                      velocity=vel(i,j+l,k,:),         &
                                      pressure=prs(i,j+l,k),           &
                                   temperature=tmp(i,j+l,k),           &
                                       species=spc(i,j+l,k,:))
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