!+---------------------------------------------------------------------+
!| This module contains subroutines related to Riemann solver.         |
!+---------------------------------------------------------------------+
!| CHANGE RECORD                                                       |
!| -------------                                                       |
!| 10-02-2022  | Created by J. Fang                                    |
!+---------------------------------------------------------------------+
module riemann
  !
  use parallel, only : mpirank
  !
  implicit none
  !
  contains
  !
  !+-------------------------------------------------------------------+
  !| this subroutine is to give Steger-Warming flux vector split       |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 10-02-2022: Created by J. Fang @ Warrington.                      |
  !+-------------------------------------------------------------------+
  subroutine flux_steger_warming(fplus,fmius,rho,vel,prs,tmp,spc,q,dxi,jacob)
    !
    use commvar,  only: numq,gamma,nondimen
    use fludyna,  only: sos
#ifdef COMB
    use thermchem,only:aceval,gammarmix
#endif
    !
    real(8),intent(out) :: fplus(:,:),fmius(:,:)
    real(8),intent(in) ::  rho(:),vel(:,:),prs(:),tmp(:),spc(:,:),   &
                           q(:,:),dxi(:,:),jacob(:)
    !
    ! local data
    real(8) :: uu,eps,gm2,css,csa,lmach,fhi,jro
    real(8) :: gpd(3),lmda(5),lmdap(5),lmdam(5)
    real(8) :: var0,var1,var2,var3,var4
    integer :: nsize,i
    !
    eps=0.04d0
    ! gm2=0.5d0/gamma
    nsize=size(rho,1)
    !
    do i=1,nsize
      !
      uu=dxi(i,1)*vel(i,1)+dxi(i,2)*vel(i,2)+dxi(i,3)*vel(i,3)
      var0=1.d0/sqrt(dxi(i,1)**2+dxi(i,2)**2+dxi(i,3)**2)
      !
      gpd(1)=dxi(i,1)*var0
      gpd(2)=dxi(i,2)*var0
      gpd(3)=dxi(i,3)*var0
      !
#ifdef COMB
      gamma = gammarmix(tmp(i),spc(i,:))
      gm2=0.5d0/gamma
      call aceval(tmp(i),spc(i,:),css)
#else
      gm2=0.5d0/gamma
      css=sos(tmp(i))
#endif
      csa=css/var0
      lmach=uu/csa
      !
      lmda(1)=uu
      lmda(2)=uu
      lmda(3)=uu
      lmda(4)=uu+csa
      lmda(5)=uu-csa
      !
      lmdap(1)=0.5d0*(lmda(1)+sqrt(lmda(1)**2+eps**2))
      lmdap(2)=0.5d0*(lmda(2)+sqrt(lmda(2)**2+eps**2))
      lmdap(3)=0.5d0*(lmda(3)+sqrt(lmda(3)**2+eps**2))
      lmdap(4)=0.5d0*(lmda(4)+sqrt(lmda(4)**2+eps**2))
      lmdap(5)=0.5d0*(lmda(5)+sqrt(lmda(5)**2+eps**2))
      !
      lmdam(1)=lmda(1)-lmdap(1)
      lmdam(2)=lmda(2)-lmdap(2)
      lmdam(3)=lmda(3)-lmdap(3)
      lmdam(4)=lmda(4)-lmdap(4)
      lmdam(5)=lmda(5)-lmdap(5)
      !
      if(lmach>=1.d0) then
        fplus(i,1)=jacob(i)*  q(i,1)*uu
        fplus(i,2)=jacob(i)*( q(i,2)*uu+dxi(i,1)*prs(i) )
        fplus(i,3)=jacob(i)*( q(i,3)*uu+dxi(i,2)*prs(i) )
        fplus(i,4)=jacob(i)*( q(i,4)*uu+dxi(i,3)*prs(i) )
        fplus(i,5)=jacob(i)*( q(i,5)+prs(i) )*uu
        !
        fmius(i,1)=0.d0
        fmius(i,2)=0.d0
        fmius(i,3)=0.d0
        fmius(i,4)=0.d0
        fmius(i,5)=0.d0  
        !
        if(numq>5) then
          fplus(i,6:numq)=jacob(i)*q(i,6:numq)*uu
          fmius(i,6:numq)=0.d0
        endif
        !
      elseif(lmach<=-1.d0) then
        fplus(i,1)=0.d0
        fplus(i,2)=0.d0
        fplus(i,3)=0.d0
        fplus(i,4)=0.d0
        fplus(i,5)=0.d0
        !
        fmius(i,1)=jacob(i)*  q(i,1)*uu
        fmius(i,2)=jacob(i)*( q(i,2)*uu+dxi(i,1)*prs(i) )
        fmius(i,3)=jacob(i)*( q(i,3)*uu+dxi(i,2)*prs(i) )
        fmius(i,4)=jacob(i)*( q(i,4)*uu+dxi(i,3)*prs(i) )
        fmius(i,5)=jacob(i)*( q(i,5)+prs(i) )*uu
        !
        if(numq>5) then
          fplus(i,6:numq)=0.d0
          fmius(i,6:numq)=jacob(i)*q(i,6:numq)*uu
        endif
        !
      else
        !
        fhi=0.5d0*(gamma-1.d0)*(vel(i,1)**2+vel(i,2)**2+vel(i,3)**2)
        !
        var1=lmdap(1)
        var2=lmdap(4)-lmdap(5)
        var3=2.d0*lmdap(1)-lmdap(4)-lmdap(5)
        var4=var1-var3*gm2
        !
        jro=jacob(i)*rho(i)
        !
        fplus(i,1)=jro*var4
        fplus(i,2)=jro*(var4*vel(i,1)+var2*css*gpd(1)*gm2)
        fplus(i,3)=jro*(var4*vel(i,2)+var2*css*gpd(2)*gm2)
        fplus(i,4)=jro*(var4*vel(i,3)+var2*css*gpd(3)*gm2)
        fplus(i,5)=jacob(i)*(var1*q(i,5)+rho(i)*(var2*uu*var0*css*gm2-var3*(fhi+css**2)*gm2/(gamma-1.d0)))
        !
        if(numq>5) then
          fplus(i,6:numq)=jacob(i)*q(i,6:numq)*var4
        endif
        !
        var1=lmdam(1)
        var2=lmdam(4)-lmdam(5)
        var3=2.d0*lmdam(1)-lmdam(4)-lmdam(5)
        var4=var1-var3*gm2
        !
        fmius(i,1)=jro*var4                                                         
        fmius(i,2)=jro*(var4*vel(i,1)+var2*css*gpd(1)*gm2)                 
        fmius(i,3)=jro*(var4*vel(i,2)+var2*css*gpd(2)*gm2)                 
        fmius(i,4)=jro*(var4*vel(i,3)+var2*css*gpd(3)*gm2)                 
        fmius(i,5)=jacob(i)*(var1*q(i,5)+rho(i)*(var2*uu*var0*css*gm2-var3*(fhi+css**2)*gm2/(gamma-1.d0)))
        !
        if(numq>5) then
          fmius(i,6:numq)=jacob(i)*q(i,6:numq)*var4
        endif
        !
        !
      end if
    enddo
    !
    return
    !
  end subroutine flux_steger_warming
  !+-------------------------------------------------------------------+
  !| The end of the subroutine flux_steger_warming.                    |
  !+-------------------------------------------------------------------+
  !!
  !
end module riemann