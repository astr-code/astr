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
  use commvar,only: mach,gamma,const2

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
  
  function ausmpwp3d(fhc_l,fhc_r,p1_l,p2_l,p1_r,p2_r) result(fh)
      !
      real(8) :: fh(5)
      real(8),intent(in) :: fhc_l(9),fhc_r(9)
      real(8),intent(in),optional :: p1_l,p2_l,p1_r,p2_r
      !
      real(8) :: ro_l,rou_l,e_l,p_l,u_l,v_l,w_l,un_l,h_l,cs_l,ma_l,ma_l_p,f_l,plp,ut_l2
      real(8) :: ro_r,rou_r,e_r,p_r,u_r,v_r,w_r,un_r,h_r,cs_r,ma_r,ma_r_m,f_r,prm,ut_r2
      real(8) :: cs,cs_m,ma_m,h_m,det_m,jac_m
      real(8) :: var1,var2,det_l,det_r
      real(8) :: w,ps
      real(8) :: mass_flux_l,mass_flux_r,m_l_plus,m_r_minus,p_r_minus,p_l_plus
      real(8) :: m_bar_l_plus,m_bar_r_minus
      real(8) :: xi_l(3),xi_r(3),ja_l,ja_r
      !
      ! left                                         ! right
      ro_l=fhc_l(1);                                 ro_r=fhc_r(1)
      u_l =fhc_l(2)/ro_l;                            u_r =fhc_r(2)/ro_r
      v_l =fhc_l(3)/ro_l;                            v_r =fhc_r(3)/ro_r
      w_l =fhc_l(4)/ro_l;                            w_r =fhc_r(4)/ro_r
      e_l =fhc_l(5);                                 e_r =fhc_r(5)   
      p_l=(e_l-0.5d0*ro_l*(u_l**2+v_l**2+w_l**2))*(gamma-1.d0)
      p_r=(e_r-0.5d0*ro_r*(u_r**2+v_r**2+w_r**2))*(gamma-1.d0)
      h_l =(e_l+p_l)/ro_l;                           h_r =(e_r+p_r)/ro_r
      xi_l(1:3)=fhc_l(6:8);                          xi_r(1:3)=fhc_r(6:8)
      ja_l=fhc_l(9);                                 ja_r=fhc_r(9)
      !
      un_l=u_l*xi_l(1)+v_l*xi_l(2)+w_l*xi_l(3);      un_r=u_r*xi_r(1)+v_r*xi_r(2)+w_r*xi_r(3)
      det_l=sqrt(xi_l(1)**2+xi_l(2)**2+xi_l(3)**2);  det_r=sqrt(xi_r(1)**2+xi_r(2)**2+xi_r(3)**2)
      un_l=un_l/det_l;                               un_r=un_r/det_r
      ! 
      ut_l2=u_l**2+v_l**2+w_l**2-un_l**2;            ut_r2=u_r**2+v_r**2+w_r**2-un_r**2
      ! det_m=0.5d0*(det_l+det_r)
      ! jac_m=0.5d0*(ja_l +ja_r)
      !               
      h_m=0.5d0*(h_l-0.5d0*ut_l2+h_r-0.5d0*ut_r2)
      ! h_m=min(h_l-0.5d0*ut_l2,h_r-0.5d0*ut_r2)
      cs=sqrt(2.d0*h_m*(gamma-1.d0)/(gamma+1.d0))
      ! !
      if(0.5d0*(un_l+un_r)<0.d0) then
        cs_m=cs**2/max(cs,abs(un_r))  
      else
        cs_m=cs**2/max(cs,abs(un_l))  
      endif
      ! cs_l=sqrt(2.d0*h_l*(gamma-1.d0)/(gamma+1.d0)); cs_r=sqrt(2.d0*h_r*(gamma-1.d0)/(gamma+1.d0))
      ! !
      ! cs_l=cs_l**2/max(abs(un_l),cs_l);              cs_r=cs_r**2/max(abs(un_r),cs_r)
      ! !
      ! cs_m=min(cs_l,cs_r)
      !
      ma_l=un_l/cs_m;                                  ma_r=un_r/cs_m
      !
      ! print*,xi_l(2),xi_r(2)
      ! print*,'================================'
      ! mach splitting functions
      m_l_plus  = mach_split_plus(ma_l)
      m_r_minus = mach_split_minus(ma_r)
      ! pressure splitting functions
      p_r_minus = pressure_split_minus(ma_r)
      p_l_plus  = pressure_split_plus(ma_l)
      !
      w = get_w(p_l=p_l, p_r=p_r)
      !
      ps=p_l * p_l_plus + p_r * p_r_minus
      !
      if(present(p1_l) .and. present(p2_l) .and. present(p1_r) .and. present(p2_r)) then
        call get_f(p_l=p_l,p_r=p_r,p_s=ps,p1_l=p1_l,p1_r=p1_r,p2_l=p2_l,p2_r=p2_r,f_l=f_l,f_r=f_r)
      else
        call get_f(p_l=p_l,p_r=p_r,p_s=ps,p1_l=p_l,p1_r=p_r,p2_l=p_l,p2_r=p_r,f_l=f_l,f_r=f_r)
      endif
      !

      if(m_l_plus + m_r_minus < 0.0d0) then
        m_bar_l_plus = m_l_plus * w * (1.0d0 + f_l)
        m_bar_r_minus = m_r_minus + m_l_plus * ((1.0d0 - w) * (1.0d0 + f_l) - f_r)
      else ! from eq. 24 (i) in ref [1]
        m_bar_l_plus = m_l_plus + m_r_minus * ((1.0d0 - w) * (1.0d0 + f_r) - f_l)
        m_bar_r_minus = m_r_minus * w * (1.0d0 + f_r)
      end if
      !
      mass_flux_l = ja_l*det_l*m_bar_l_plus*cs_m*ro_l
      mass_flux_r = ja_r*det_r*m_bar_r_minus*cs_m*ro_r
      !
      p_l_plus=p_l_plus*ja_l*p_l
      p_r_minus=p_r_minus*ja_r*p_r
      !
      !
      fh(1)=mass_flux_l+mass_flux_r
      fh(2)=mass_flux_l*u_l+mass_flux_r*u_r+p_l_plus*xi_l(1)+p_r_minus*xi_r(1)
      fh(3)=mass_flux_l*v_l+mass_flux_r*v_r+p_l_plus*xi_l(2)+p_r_minus*xi_r(2)
      fh(4)=mass_flux_l*w_l+mass_flux_r*w_r+p_l_plus*xi_l(3)+p_r_minus*xi_r(3)
      fh(5)=(mass_flux_l*h_l)+(mass_flux_r*h_r)
      !
    end function ausmpwp3d

    pure subroutine get_f(p_L, p_R, p_s, p1_L, p1_R, p2_L, p2_R, f_L, f_R)
      !< Implementation of the f function (Eq. 26 in Ref [1]). This is used to help determine
      !< where shock discontinuities are
      real(8), intent(in) :: p_L    !< L pressure
      real(8), intent(in) :: p_R    !< R pressure
      real(8), intent(in) :: p1_L   !< L pressure, j+1 for i direction, i-1 for j direction
      real(8), intent(in) :: p1_R   !< R pressure, j-1 for i direction, i+1 for j direction
      real(8), intent(in) :: p2_L   !< L pressure, j+1 for i direction, i-1 for j direction
      real(8), intent(in) :: p2_R   !< R pressure, j-1 for i direction, i+1 for j direction
      real(8), intent(in) :: p_s !< p_s = p_L * P_L(+) + p_R * P_R(-)

      real(8), intent(out) :: f_L !< L "f" function
      real(8), intent(out) :: f_R !< R "f" function

      real(8) :: min_term

      if(abs(p_s) < tiny(1.0d0)) then ! p_s == 0
        f_L = 0.0d0
        f_R = 0.0d0
      else
        min_term = min(1.0d0, min(p1_L, p1_R, p2_L, p2_R) / min(p_L, p_R))**2
        f_L = ((p_L / p_s) - 1.0d0) * min_term
        f_R = ((p_R / p_s) - 1.0d0) * min_term
      end if
    end subroutine get_f
    !
    pure real(8) function get_w(p_L, p_R) result(w)
      !< Discontinuity sensor w_1 (Eq. 25 in Ref [1]). This detects whether a shock
      !< exists in the normal direction to the cell-interface or not
      real(8), intent(in) :: p_L !< left interface pressure
      real(8), intent(in) :: p_R !< right interface pressure

      real(8) :: min_term

      min_term = min((p_L / p_R),(p_R / p_L))
      w = 1.0d0 - min_term * min_term * min_term
    end function get_w
    !
    pure real(8) function mach_split_plus(M) result(M_plus)
      !< The Mach splitting function (Eq. 9 in Ref[1]). This is the M+ version. This is kept
      !< simple for inlining
      real(8), intent(in) :: M !< interface Mach number

      if(abs(M) > 1.0d0) then
        M_plus = 0.5d0 * (M + abs(M))
      else ! |M| <= 1
        M_plus = 0.25d0 * (M + 1.0d0)**2
      endif
    end function mach_split_plus
    !
    pure real(8) function mach_split_minus(M) result(M_minus)
      !< The Mach splitting function (Eq. 9 in Ref [1]). This is the M- version. This is kept
      !< simple for inlining
      real(8), intent(in) :: M !< interface Mach number
      !
      if(abs(M) > 1.0d0) then
        M_minus = 0.5d0 * (M - abs(M))
      else ! |M| <= 1
        M_minus = -0.25d0 * (M - 1.0d0)**2
      endif
    end function mach_split_minus
    !
    pure real(8) function pressure_split_plus(M) result(P_plus)
      !< The pressure splitting function (Eq. 10 in Ref[1]). This is the P+ version
      real(8), intent(in) :: M !< interface Mach number

      ! real(8) :: alfa=0.d0
  
      if(abs(M) > 1.0d0) then
        P_plus = 0.5d0 * (1.0d0 + sign(1.0d0, M))
      else ! |M| <= 1
        P_plus = 0.25d0 * (M + 1.0d0)**2 * (2.0d0 - M) + 0.d0 * M * (M**2 - 1.0d0)**2
      endif
    end function pressure_split_plus
    !
    pure real(8) function pressure_split_minus(M) result(P_minus)
      !< The pressure splitting function (Eq. 10 in Ref[1]). This is the P- version
      real(8), intent(in) :: M !< interface Mach number
  
      if(abs(M) > 1.0d0) then
        P_minus = 0.5d0 * (1.0d0 - sign(1.0d0, M))
      else ! |M| <= 1
        P_minus = 0.25d0 * (M - 1.0d0)**2 * (2.0d0 + M) - 0.d0 * M * (M**2 - 1.0d0)**2
      endif
    end function pressure_split_minus
    
end module riemann