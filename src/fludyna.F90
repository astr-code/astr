!+---------------------------------------------------------------------+
!| This module contains subroutines and functions related to fluid     |
!| dynamicsof testing hanmish                                          |
!| ==============                                                      |
!| CHANGE RECORD                                                       |
!| -------------                                                       |
!| 09-Jnu-2020  | Created by J. Fang @ STFC Daresbury Laboratory       |
!+---------------------------------------------------------------------+
module fludyna
  !
  ! use parallel, only: mpirank,mpistop
  use commvar, only: nondimen
  use thermchem, only: rgcmix,cpeval,temperature_calc
  !
  implicit none
  !
  interface thermal
     module procedure thermal_scar
     module procedure thermal_1d
     module procedure thermal_3d
  end interface
  !
  interface fvar2q
     module procedure fvar2q_sca
     module procedure fvar2q_1da
     module procedure fvar2q_3da
  end interface
  !
  interface q2fvar
     module procedure q2fvar_sca
     module procedure q2fvar_1da
     module procedure q2fvar_3da
  end interface
  !
  contains
  !
  !+-------------------------------------------------------------------+
  !| This function is used to get thermal relation among density,      |
  !| pressure and temperature                                          |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 13-May-2020: Created by J. Fang @ STFC Daresbury Laboratory       |
  !+-------------------------------------------------------------------+
  function thermal_scar(density,pressure,temperature,species) result(vout)
    !
    use commvar,only : const2,rgas
    !
    ! arguments
    real(8) :: vout
    real(8),intent(in) ,optional :: density,pressure,temperature,species(:)
    !
    real(8) :: rloc

    
    if(nondimen) then
      !
      if(present(density) .and. present(temperature)) then
        vout=density*temperature/const2
      elseif(present(density) .and. present(pressure)) then
        vout=pressure/density*const2
      elseif(present(temperature) .and. present(pressure)) then
        vout=pressure/temperature*const2
      else
        stop ' !! unable to get thermal variable  @ thermal !!'
      endif
      !
    else 
      !
#ifdef COMB
      rloc=rgcmix(species)
#else
      rloc=rgas
#endif
      !
      if(present(density).and.present(temperature)) then
        vout = density*temperature*rloc
      elseif(present(density).and.present(pressure)) then
        vout = pressure/density/rloc
      elseif(present(temperature).and.present(pressure)) then
        vout = pressure/temperature/rloc
      else
        stop ' !! unable to use dimensional EoS @ thermal !!'
      endif
      !
    endif 
    
  end function thermal_scar
  !
  function thermal_1d(density,pressure,temperature,species,dim) result(vout)
    !
    use commvar,only : const2,rgas
    !
    ! arguments
    integer,intent(in) :: dim
    real(8) :: vout(dim)
    real(8),intent(in),optional :: density(:),pressure(:),        &
                                   temperature(:),species(:,:)
    !
    real(8) :: rloc(dim)
    !
    if(nondimen) then
      !
      if(present(density) .and. present(temperature)) then
        vout=density*temperature/const2
      elseif(present(density) .and. present(pressure)) then
        vout=pressure/density*const2
      elseif(present(temperature) .and. present(pressure)) then
        vout=pressure/temperature*const2
      else
        stop ' !! unable to get thermal variable  @ thermal !!'
      endif
      !
    else 
      !
#ifdef COMB
      rloc=rgcmix(species,dim)
#else
      rloc=rgas
#endif

      if(present(density).and.present(temperature)) then
        vout = density*temperature*rloc
      elseif(present(density).and.present(pressure)) then
        vout = pressure/density/rloc
      elseif(present(temperature).and.present(pressure)) then
        vout = pressure/temperature/rloc
      else
        stop ' !! unable to use dimensional EoS @ thermal !!'
      endif
      !
    endif 
    !
  end function thermal_1d
  !
  function thermal_3d(density,pressure,temperature,species,dim) result(vout)
    !
    use commvar,only : const2,rgas
    !
    ! arguments
    integer,intent(in) :: dim(3)
    real(8) :: vout(dim(1),dim(2),dim(3))
    real(8),intent(in),optional :: density(:,:,:),pressure(:,:,:),    &
                                   temperature(:,:,:),species(:,:,:,:)
    !
    real(8) :: rloc(dim(1),dim(2),dim(3))
    !
    if(nondimen) then
      !
      if(present(density) .and. present(temperature)) then
        vout=density*temperature/const2
      elseif(present(density) .and. present(pressure)) then
        vout=pressure/density*const2
      elseif(present(temperature) .and. present(pressure)) then
        vout=pressure/temperature*const2
      else
        stop ' !! unable to get thermal variable  @ thermal_3d !!'
      endif
      !
    else
      !
#ifdef COMB
      rloc=rgcmix(species,dim)
#else
      rloc=rgas
#endif
      if(present(density).and.present(temperature)) then
        vout = density*temperature*rloc
      elseif(present(density).and.present(pressure)) then
        vout = pressure/density/rloc
      elseif(present(temperature).and.present(pressure)) then
        vout = pressure/temperature/rloc
      else
        stop ' !! unable to use dimensional EoS @ thermal_3d !!'
      endif
      !
    endif
    !
  end function thermal_3d
  !+-------------------------------------------------------------------+
  !| The end of the subroutine thermal.                                |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to update flow variables from q.               |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 04-Aug-2018: Created by J. Fang @ STFC Daresbury Laboratory       |
  !+-------------------------------------------------------------------+
  subroutine updatefvar
    !
    use commarray,only : q,rho,vel,prs,tmp,spc,tke,omg
    use commvar,  only : im,jm,km,num_species,num_modequ,turbmode
    !
    ! local data
    integer :: i,j,k
    !
    if(trim(turbmode)=='k-omega') then
      !
      call q2fvar(q=q(0:im,0:jm,0:km,:),                               &
                                     density=rho(0:im,0:jm,0:km),      &
                                    velocity=vel(0:im,0:jm,0:km,:),    &
                                    pressure=prs(0:im,0:jm,0:km),      &
                                 temperature=tmp(0:im,0:jm,0:km),      &
                                     species=spc(0:im,0:jm,0:km,:),    &
                                         tke=tke(0:im,0:jm,0:km),      &
                                       omega=omg(0:im,0:jm,0:km) )
      !
      do k=0,km
      do j=0,jm
      do i=0,im
        !
        if(tke(i,j,k)<0.d0) then
          tke(i,j,k)=0.d0
          q(i,j,k,6+num_species)=rho(i,j,k)*tke(i,j,k)
        endif
        !
        if(omg(i,j,k)<0.01d0) then
          omg(i,j,k)=0.01d0
          q(i,j,k,7+num_species)=rho(i,j,k)*omg(i,j,k)
        endif
        !
      enddo
      enddo
      enddo
      !
    elseif(trim(turbmode)=='none' .or. trim(turbmode)=='udf1') then
      !
      call q2fvar(q=q(0:im,0:jm,0:km,:),                               &
                                     density=rho(0:im,0:jm,0:km),      &
                                    velocity=vel(0:im,0:jm,0:km,:),    &
                                    pressure=prs(0:im,0:jm,0:km),      &
                                 temperature=tmp(0:im,0:jm,0:km),      &
                                     species=spc(0:im,0:jm,0:km,:) )
      !
    else
      print*,' !! ERROR @ updatefvar'
      stop
    endif
    !
  end subroutine updatefvar
  !+-------------------------------------------------------------------+
  !| The end of the subroutine updatefvar.                             |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to update q from flow variables.               |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 06-Aug-2018: Created by J. Fang @ STFC Daresbury Laboratory       |
  !+-------------------------------------------------------------------+
  subroutine updateq
    !
    use commarray,only : q,rho,vel,prs,tmp,spc,tke,omg
    use commvar,  only : im,jm,km,num_species,num_modequ,turbmode,numq
    !
    integer :: i,j,k,n
    !
    if(trim(turbmode)=='k-omega') then
      !
      call fvar2q(          q=  q(0:im,0:jm,0:km,:),                   &
                      density=rho(0:im,0:jm,0:km),                     &
                     velocity=vel(0:im,0:jm,0:km,:),                   &
                     pressure=prs(0:im,0:jm,0:km),                     &
                      species=spc(0:im,0:jm,0:km,:),                   &
                          tke=tke(0:im,0:jm,0:km),                     &
                        omega=omg(0:im,0:jm,0:km)                      )
      !
    elseif(trim(turbmode)=='none' .or. trim(turbmode)=='udf1') then
      !
      call fvar2q(          q=  q(0:im,0:jm,0:km,:),                   &
                      density=rho(0:im,0:jm,0:km),                     &
                     velocity=vel(0:im,0:jm,0:km,:),                   &
                     pressure=prs(0:im,0:jm,0:km),                     &
                     temperature=tmp(0:im,0:jm,0:km),                  &
                      species=spc(0:im,0:jm,0:km,:)                    )
      !
      do k=0,km
      do j=0,jm
      do i=0,im
        !
        do n=1,numq
          if(isnan(q(i,j,k,n))) then
            print*,i,j,k,n
            print*,'q:',q(i,j,k,:),'@ updateq'
          endif
        enddo
        !
      enddo
      enddo
      enddo
      !
    else
      print*,' !! ERROR @ updatefvar'
      stop
    endif
    !
  end subroutine updateq
  !+-------------------------------------------------------------------+
  !| The end of the subroutine updateq.                                |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This function is used to convert flow variables to q.             |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 09-02-2021: Created by J. Fang @ Warrington.                      |
  !+-------------------------------------------------------------------+
  subroutine fvar2q_sca(q,density,velocity,pressure,temperature,       &
                       species,tke,omega)
    !
    use commvar, only: numq,ndims,num_species,const1,const6,cv
    !
    real(8),intent(in) :: density,velocity(:)
    real(8),intent(in),optional :: pressure,temperature,species(:),tke,omega
    real(8),intent(out) :: q(:)
    !
    ! local data
    integer :: jspec,j
    real(8) :: var1,var2,cotem
    !
    q(1)=density
    q(2)=density*velocity(1)
    q(3)=density*velocity(2)
    q(4)=density*velocity(3)
    !
#ifdef COMB
    if(.not.present(species))  &
      stop ' !! error @ fvar2q - species not set for energy calculation !!'
      !
    if(present(temperature)) then
      var1=0.5d0*sum(velocity(:)*velocity(:))
      call cpeval(tmp=temperature,spc=species,ke=var1,eng=var2)
      q(5)=var2*density
    else
      print*,' !! temperature required for energy calculation !!'
      stop ' !! error @ fvar2q'
    endif 
    !
#else
    if(nondimen) then 
      cotem=const1
    else
      cotem=cv
    endif
    !
    var1=0.5d0*sum(velocity(:)*velocity(:))
    if(present(temperature)) then
      q(5)=density*(temperature*cotem+var1)
    elseif(present(pressure)) then
      q(5)=pressure*const6+density*var1
    else
      print*,' !! pressure or temperature required for energy calculation !!'
      stop ' !! error @ fvar2q'
    endif
#endif
    !
    if(num_species>0) then
      !
      do jspec=1,num_species
        q(5+jspec)=density*species(jspec)
      enddo
      !
    endif
    !
    if(present(tke) .and. present(omega)) then
      !
      q(5+num_species+1)=density*tke
      q(5+num_species+2)=density*omega
      !
    endif
    !
  end subroutine fvar2q_sca
  !
  subroutine fvar2q_1da(q,density,velocity,pressure,temperature,species)
    !
    use commvar, only: numq,ndims,num_species,const1,const6,cv
    !
    real(8),intent(in) :: density(:),velocity(:,:)
    real(8),intent(in),optional :: pressure(:),temperature(:), &
                                   species(:,:)
    real(8),intent(out) :: q(:,:)
    !
    ! local data
    integer :: jspec,j
    real(8) :: var1,var2,cotem
    !
    q(:,1)=density(:)
    q(:,2)=density(:)*velocity(:,1)
    q(:,3)=density(:)*velocity(:,2)
    q(:,4)=density(:)*velocity(:,3)
    !
#ifdef COMB

    if(.not.present(species))  &
      stop ' !! error @ fvar2q - species not set for energy calculation !!'
    !
    if(present(temperature)) then
      !
      do j=1,size(density)
        var1=0.5d0*sum(velocity(j,:)*velocity(j,:))
        call cpeval(tmp=temperature(j),spc=species(j,:),ke=var1,eng=var2)
        q(j,5)=var2*density(j)
      enddo   
      !
    else
      print*,' !! temperature required for energy calculation !!'
      stop ' !! error @ fvar2q'
    endif 

#else
    
    if(nondimen) then 
      cotem=const1
    else
      cotem=cv
    endif
    !
    if(present(temperature)) then
      q(:,5)=density(:)*( temperature(:)*cotem +         &
                          0.5d0*(velocity(:,1)**2 +      &
                                 velocity(:,2)**2 +      &
                                 velocity(:,3)**2) )
    elseif(present(pressure)) then
      q(:,5)=pressure(:)*const6+0.5d0*density(:)*(       &
                                    velocity(:,1)**2 +   &
                                    velocity(:,2)**2 +   &
                                    velocity(:,3)**2 )
    else
      print*,' !! pressure or temperature required for energy calculation !!'
      stop ' !! error @ fvar2q'
    endif

#endif
    !
    if(num_species>0) then
      !
      do jspec=1,num_species
        q(:,5+jspec)=density(:)*species(:,jspec)
      enddo
      !
    endif
    !
  end subroutine fvar2q_1da
  !
  subroutine fvar2q_3da(q,density,velocity,pressure,temperature,species,tke,omega)
    !
    use commvar, only: numq,ndims,num_species,const1,const6,cv
    !
    real(8),intent(in) :: density(:,:,:),velocity(:,:,:,:)
    real(8),intent(in),optional :: pressure(:,:,:),temperature(:,:,:), &
                                   species(:,:,:,:),tke(:,:,:),        &
                                   omega(:,:,:)
    real(8),intent(out) :: q(:,:,:,:)
    !
    ! local data
    integer :: jspec,i,j,k
    real(8) :: var1,var2,cotem
    !
    q(:,:,:,1)=density(:,:,:)
    q(:,:,:,2)=density(:,:,:)*velocity(:,:,:,1)
    q(:,:,:,3)=density(:,:,:)*velocity(:,:,:,2)
    q(:,:,:,4)=density(:,:,:)*velocity(:,:,:,3)
    !
#ifdef COMB
    
    if(.not.present(species))  &
      stop ' !! error @ fvar2q - species not set for energy calculation !!'
      !
    if(present(temperature)) then
      !
      do i=1,size(density,1)
        do j=1,size(density,2)
          do k=1,size(density,3)
            !
            var1=0.5d0*sum(velocity(i,j,k,:)*velocity(i,j,k,:))
            call cpeval(tmp=temperature(i,j,k),spc=species(i,j,k,:), &
                        ke=var1,eng=var2)
            q(i,j,k,5)=var2*density(i,j,k)
            !
          enddo   
        enddo 
      enddo 
      !
    else
      print*,' !! temperature required for energy calculation !!'
      stop ' !! error @ fvar2q'
      endif

#else

    if(nondimen) then 
      cotem=const1
    else
      cotem=cv
    endif
    !
    if(present(temperature)) then
      q(:,:,:,5)=density(:,:,:)*( temperature(:,:,:)*cotem +     &
                           0.5d0*(velocity(:,:,:,1)**2 +         &
                                  velocity(:,:,:,2)**2 +         &
                                  velocity(:,:,:,3)**2) )
    elseif(present(pressure)) then
        q(:,:,:,5)=pressure(:,:,:)*const6+0.5d0*density(:,:,:)*(       &
                                    velocity(:,:,:,1)**2 +             &
                                    velocity(:,:,:,2)**2 +             &
                                    velocity(:,:,:,3)**2 )
    else
      print*,' !! pressure or temperature required !!'
      stop ' !! error @ fvar2q'
    endif

#endif
    !
    if(num_species>0) then
      !
      do jspec=1,num_species
        q(:,:,:,5+jspec)=density(:,:,:)*species(:,:,:,jspec)
      enddo
      !
    endif
    !
    if(present(tke) .and. present(omega)) then
      !
      q(:,:,:,5+num_species+1)=tke(:,:,:)*density(:,:,:)
      q(:,:,:,5+num_species+2)=omega(:,:,:)*density(:,:,:)
      !
    endif
    !
  end subroutine fvar2q_3da
  !+-------------------------------------------------------------------+
  !| The end of the subroutine fvar2q.                                 |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine calcualtes field variables from q.                |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 09-02-2021: Created by J. Fang @ Warrington.                      |
  !+-------------------------------------------------------------------+
  subroutine q2fvar_3da(q,density,velocity,pressure,temperature,species,tke,omega)
    !
    use commvar, only: numq,ndims,num_species,const1,const6,tinf
    !
    real(8),intent(inout) :: q(:,:,:,:)
    real(8),intent(out) :: density(:,:,:)
    real(8),intent(out),optional :: temperature(:,:,:),pressure(:,:,:),&
                                    species(:,:,:,:),velocity(:,:,:,:),&
                                    tke(:,:,:),omega(:,:,:)
    !
    ! local data
    integer :: jspec,i,j,k
    integer :: dim(3)
    real(8) :: var1
    !
    density(:,:,:)   =q(:,:,:,1)
    !
    if(present(velocity) .or. present(pressure) .or. present(temperature)) then
      velocity(:,:,:,1)=q(:,:,:,2)/density
      velocity(:,:,:,2)=q(:,:,:,3)/density
      velocity(:,:,:,3)=q(:,:,:,4)/density
    endif
    !
    if(num_species>0 .and. present(species)) then
      !
      do jspec=1,num_species
        species(:,:,:,jspec)=q(:,:,:,5+jspec)/density(:,:,:)
      enddo
      !
    endif
    !
#ifdef COMB

    if(.not.present(species))  &
      stop ' !! error @ q2fvar - species not set for T calculation !!'
      !
    if(present(pressure) .or. present(temperature)) then
      !
      do i=1,size(q,1)
        do j=1,size(q,2)
          do k=1,size(q,3)
            !
            var1=q(i,j,k,5)/density(i,j,k) &
                 -0.5d0*sum(velocity(i,j,k,:)*velocity(i,j,k,:))
            ! temperature(i,j,k)=tinf
            call temperature_calc(tmp=temperature(i,j,k),den=density(i,j,k), &
                                  spc=species(i,j,k,:),eint=var1)
          enddo
        enddo 
      enddo          
      !
    endif 
    !
    if(present(pressure)) then
      dim(1)=size(q,1)
      dim(2)=size(q,2)
      dim(3)=size(q,3)
      !
      pressure=thermal(temperature=temperature,density=density,species=species,dim=dim)

    endif
    !  
#else

    if(present(pressure) .or. present(temperature)) then
      pressure(:,:,:)  =( q(:,:,:,5)-0.5d0*density(:,:,:)*(              &
                                          velocity(:,:,:,1)**2+         &
                                          velocity(:,:,:,2)**2+         &
                                          velocity(:,:,:,3)**2) )/const6
    endif
    !
    if(present(temperature)) then
      dim(1)=size(q,1)
      dim(2)=size(q,2)
      dim(3)=size(q,3)
      !
      temperature=thermal(pressure=pressure,density=density,dim=dim)
    endif
    !
#endif
    !
    if(present(tke)) then
      tke(:,:,:)=q(:,:,:,5+num_species+1)/density
    endif
    !
    if(present(omega)) then
      omega(:,:,:)=q(:,:,:,5+num_species+2)/density
    endif
    !
  end subroutine q2fvar_3da
  !
  subroutine q2fvar_1da(q,density,velocity,pressure,temperature,species,tke,omega)
    !
    use commvar, only: numq,ndims,num_species,const1,const6
    !
    real(8),intent(in) :: q(:,:)
    real(8),intent(out) :: density(:)
    real(8),intent(out),optional :: velocity(:,:),pressure(:),         &
                                    temperature(:),species(:,:),       &
                                    tke(:),omega(:)
    !
    ! local data
    integer :: jspec,j
    integer :: dim
    real(8) :: var1
    !
    density(:)   =q(:,1)
    !
    if(present(velocity) .or. present(pressure) .or. present(temperature)) then
      velocity(:,1)=q(:,2)/density
      velocity(:,2)=q(:,3)/density
      velocity(:,3)=q(:,4)/density
    endif
    !
    if(num_species>0 .and. present(species)) then
      !
      do jspec=1,num_species
        species(:,jspec)=q(:,5+jspec)/density(:)
      enddo
      !
    endif
    !
#ifdef COMB

    if(.not.present(species))  &
      stop ' !! error @ q2fvar - species not set for T calculation !!'
      !
    if(present(pressure) .or. present(temperature)) then
      !
      do j=1,size(density)
        !
        var1=q(j,5)/density(j)-0.5d0*sum(velocity(j,:)*velocity(j,:))
        call temperature_calc(tmp=temperature(j),den=density(j), &
                                  spc=species(j,:),eint=var1)
      enddo          
      !
    endif 
    !
    if(present(pressure)) then
      dim=size(q,1)
      !
      pressure=thermal(temperature=temperature,density=density, &
                        species=species,dim=dim)
    endif

#else

    if(present(pressure) .or. present(temperature)) then
      pressure(:)  =( q(:,5)-0.5d0*density(:)*(                   &
                                          velocity(:,1)**2+      &
                                          velocity(:,2)**2+      &
                                          velocity(:,3)**2) )/const6
    endif
    !
    if(present(temperature)) then
      dim=size(q,1)
      !
      temperature=thermal(pressure=pressure,density=density,dim=dim)
    endif

#endif
    !
    if(present(tke)) then
      tke(:)=q(:,5+num_species+1)/density
    endif
    !
    if(present(omega)) then
      omega(:)=q(:,5+num_species+2)/density
    endif
    !
  end subroutine q2fvar_1da
  !
  subroutine q2fvar_sca(q,density,velocity,pressure,temperature,species,tke,omega)
    !
    use commvar, only: numq,ndims,num_species,const1,const6
    !
    real(8),intent(in) :: q(:)
    real(8),intent(out) :: density
    real(8),intent(out),optional :: velocity(:),pressure,temperature,  &
                                    species(:),tke,omega
    !
    ! local data
    integer :: jspec
    real(8) :: var1
    !
    density   =q(1)
    !
    if(present(velocity) .or. present(pressure) .or. present(temperature)) then
      velocity(1)=q(2)/density
      velocity(2)=q(3)/density
      velocity(3)=q(4)/density
    endif
    !
    if(num_species>0 .and. present(species)) then
      !
      do jspec=1,num_species
        species(jspec)=q(5+jspec)/density
      enddo
      !
    endif
    !
#ifdef COMB
    if(.not.present(species))  &
      stop ' !! error @ q2fvar - species not set for T calculation !!'
      !
    if(present(pressure) .or. present(temperature)) then
      !
      var1=q(5)/density-0.5d0*sum(velocity(:)*velocity(:))
      call temperature_calc(tmp=temperature,den=density, &
                            spc=species(:),eint=var1)
    endif 
    !
    if(present(pressure)) then
      pressure=thermal(temperature=temperature,density=density, &
                        species=species)
    endif
#else

    if(present(pressure) .or. present(temperature)) then
      pressure  =( q(5)-0.5d0*density*(velocity(1)**2+velocity(2)**2+    &
                                       velocity(3)**2) )/const6
    endif
    !
    if(present(temperature)) then
      temperature=thermal(pressure=pressure,density=density)
    endif

#endif
    !
    if(present(tke)) then
      tke=q(5+num_species+1)/density
    endif
    !
    if(present(omega)) then
      omega=q(5+num_species+2)/density
    endif
    !
  end subroutine q2fvar_sca
  !+-------------------------------------------------------------------+
  !| The end of the subroutine q2fvar.                                 |
  !+-------------------------------------------------------------------+
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! calculate dynamic viscosity coefficient at different temperature.
  ! ref: https://www.cfd-online.com/Wiki/Sutherland%27s_law
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  real(8) function miucal(temper)
    !
    use commvar, only :  tempconst,tempconst1
    !
    real(8),intent(in) :: temper
    ! temper represent temperature, dimensionless
    ! below calculate miucal using sutherland's law
    ! tempconst=110.4d0/ref_t
    ! tempconst1=1.d0+tempconst
    !
    real(8) :: s,tmpers,tmper0,miu0,tnondim
    !
    if(nondimen) then
      miucal=temper*sqrt(temper)*tempconst1/(temper+tempconst)
    else
      ! for air, kg/(m s)
      tnondim=temper/273.15d0
      !
      miucal=1.716d-5*tnondim*sqrt(tnondim)*(273.15d0+110.4d0)/(temper+110.4d0)
      !
      ! tmpers=110.d0
      ! tmper0=273.125d0
      ! miu0=1.845d-5
      ! miucal=miu0*(tmper0+tmpers)/(temper+tmpers) &
      !        *((temper/tmper0)**1.5d0)
    endif 
    !
    return
    !
  end function miucal
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! End of function MiuCal.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !+-------------------------------------------------------------------+
  !| This function is used to calculate speed of sound.                |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 24-Feb-2021: Created by J. Fang @ STFC Daresbury Laboratory       |
  !+-------------------------------------------------------------------+
  real(8) function sos(tmp,spc)
    !
    use commvar,   only: nondimen,mach,num_species,gamma,rgas
    use thermchem, only: aceval
    !
    real(8),intent(in) :: tmp
    real(8),intent(in),optional :: spc(:)
    !
    ! local data
    real(8) :: cpcmix,gamrgc
    !
#ifdef COMB
    if(.not. present(spc)) then
      stop ' !! error, species needed to calculate speed of sound @ function sos'
    else
      call aceval(tmp,spc,sos)
    endif
#else
    if(nondimen) then
      sos=sqrt(tmp)/mach
    else
      sos=sqrt(gamma*rgas*tmp)
    endif
#endif
    !
    return
    !
  end function sos
  !+-------------------------------------------------------------------+
  !| The end of the function sos.                                      |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This function is used to assign the velocity profile of the jet   |
  !| flow.                                                             |
  !+-------------------------------------------------------------------+
  !| ref: Bogey, C., Bailly C., and Juve,D., Noise Investigation of a  |
  !|      High Subsonic, Moderate Reynolds Number Jet Using a          |
  !|      Compressible Large Eddy Simulation. Theoret. Comput. Fluid   |
  !|      Dynamics, 2003, 16, 273–297.                                 |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 09-Jun-2020: Created by J. Fang @ STFC Daresbury Laboratory       |
  !+-------------------------------------------------------------------+
  function jetvel(r,time) result(u)
    !
    use constdef
    use commvar, only: uinf,ymin,ymax
    !
    ! arguments
    real(8) :: u(3)
    real(8),intent(in) :: r
    real(8),intent(in),optional :: time
    !
    ! local data
    real(8) :: r0,uc,delta,ujet,alfa,theter,Tperi,var1
    !+-----------------------------+
    !| delta: the initial momentum |
    !| thickness of the shear layer|
    !| uc: jet centerline velocity |
    !+-----------------------------+
    !
    !
    r0=1.d0
    uc=1.67d0*uinf    ! ref: Reichert & Biringen, Mechanics Research
                      !      Communications 34 (2007) 249–259
    delta=0.05d0*r0   ! ref: Bogey
    !
    var1=0.5d0+0.5d0*tanh((r0-abs(r))/(2.d0*delta)) ! profile 0~1
    ! ujet=var1*(uc-uinf)+uinf
    ujet=var1
    !
    if(present(time)) then
      !
      ! alfa=2.5d0/180.d0*pi
      ! theter=sin(2.d0*pi/Tperi*time)*alfa
      ! !
      ! u(1)=ujet*cos(theter)+uinf
      ! u(2)=ujet*sin(theter)
      !
    else
      u(1)=ujet !+uinf
      u(2)=0.d0
      u(3)=0.d0
    endif
    !
    return
    !
  end function jetvel
  !+-------------------------------------------------------------------+
  !| The end of the function jetvel.                                   |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This function is used to assign velocity profile of mixing layer. |
  !+-------------------------------------------------------------------+
  !| ref: Li, Z., Jaberi, F. 2010. Numerical Investigations of         |
  !|      Shock-Turbulence Interaction in a Planar Mixing Layer.       |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 19-Mar-2021: Created by J. Fang @ STFC Daresbury Laboratory       |
  !+-------------------------------------------------------------------+
  subroutine mixinglayervel(y,u,t,ro,info) 
    !
    use constdef
    use commvar, only: ymin,ymax,gamma,mach
    !
    ! arguments
    real(8),intent(in) :: y
    real(8),intent(out) :: u(3),t,ro
    logical,intent(in) :: info
    !
    ! local data
    real(8) :: uref,u1,u2,delta,mc,ustar
    logical,save :: firstcall=.true.
    !+-----------------------------------+
    !| delta: inflow vorticity thickness |
    !+-----------------------------------+
    u1=1.5d0
    u2=0.5d0
    uref = 0.5d0*(u1-u2)
    !
    delta=1.d0
    mc=(u1-u2)/2.d0*mach
    !
    u(1)=0.5d0*(u1+u2)+uref*tanh(2.d0*y/delta)
    u(2)=0.d0
    u(3)=0.d0
    !
    ustar=u(1)-0.5d0*(u1+u2)
    t=1.d0*(1.d0+mc**2*0.5d0*(gamma-1.d0)*(1.d0-(ustar/uref)**2))
    !
    ro=1.d0/t
    !
    if(info .and. firstcall) then
      print*,' ** High-speed Mach number=',u1*mach
      print*,' **  Low-speed Mach number=',u2*mach
      print*,' ** Convective Mach number=',mc
      print*,' **  Shear-layer thickness=',delta
      firstcall=.false.
    endif
    !
    return
    !
  end subroutine mixinglayervel
  !+-------------------------------------------------------------------+
  !| The end of the function mixinglayervel.                           |
  !+-------------------------------------------------------------------+
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! This subroutine is used to cauculate the post-shock flow paramters,
  ! by using pre-shock flow paramters and shock angle and invisid R-H
  ! Relation.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Writen by Fang Jian, 2010-10-23.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine postshock(ropre,upre,vpre,ppre,tpre,ropos,upos,vpos,      &
                                                          ppos,tpos,ang)
    !
    use commvar, only: const2,gamma,mach
    use constdef,only: pi
    !
    real(8),intent(in)  :: ropre,upre,vpre,ppre,tpre,ang
    real(8),intent(out) :: ropos,upos,vpos,ppos,tpos
    !
    real(8) :: csspre,unpre,utpre,machnpre,unpos,utpos,sratio,ang2
    !
    ang2=ang/180.d0*pi
    !
    csspre=dsqrt(tpre)/mach
    !
    unpre=upre*dsin(ang2)-vpre*dcos(ang2)
    utpre=upre*dcos(ang2)+vpre*dsin(ang2)
    machnpre=unpre/csspre
    ! print*,'machnpre=',machnpre
    ! print*,'unpre=',unpre
    ! print*,'csspre=',csspre
    sratio=2.d0*gamma/(gamma+1.d0)*machnpre**2-(gamma-1.d0)/(gamma+1.d0)
    !
    ppos=ppre*sratio
    ropos=(sratio*(gamma+1.d0)+gamma-1.d0)/                            &
                                  (gamma+1.d0+sratio*(gamma-1.d0))*ropre
    tpos=ppos/ropos*const2
    !
    unpos=ropre*unpre/ropos
    utpos=utpre
    upos=unpos*dsin(ang2)+utpos*dcos(ang2)
    vpos=-unpos*dcos(ang2)+utpos*dsin(ang2)
    !
  end subroutine postshock
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! end of the subroutine postshockcal.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !+-------------------------------------------------------------------+
  !| This function is used to assign the inflow condition for the three|
  !| stream coaxial jet flow.                                          |
  !+-------------------------------------------------------------------+
  !| ref: Bouheraoua, L., Domingo, P., Ribert, G. Large-eddy simulation|
  !| of a supersonic lifted jet flame: Analysis of the turbulent flame |
  !| base. Combust. Flame, 2017, 179, 199-218.                         |                                               |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 04-Feb-2022: Created by Z.X. Chen @ Hengdong                      |
  !+-------------------------------------------------------------------+
  subroutine multistream_inflow(stream,rovd,rho,vel,prs,tmp,spc,ubulk)
    !
    use commvar
#ifdef COMB
    use thermchem, only: spcindex,convertxiyi
#endif
    !
    ! arguments
    character(len=*),intent(in) :: stream
    real(8),intent(in),optional :: rovd
    real(8),intent(out),optional :: &
      rho,prs,tmp,spc(num_species),vel(ndims),ubulk
    !
#ifdef COMB
    !local
    real(8) :: vel_ref(ndims),vel0(ndims),ub,spcx(num_species),h,delta &
              ,spc_ref(num_species),tmp_ref
    !
    vel_ref(:)=0.d0
    spc_ref(:)=0.d0
    spcx(:)=0.d0
    tmp=tinf
    !
    if(trim(stream)=='fuel') then
      !
      tmp_ref=545.d0
      spcx(spcindex('H2'))=1.d0
      ! spcx(spcindex('N2'))=1.d0-sum(spcx(:))
      call convertxiyi(spcx(:),spc_ref(:),'X2Y')
      spc(:)=spc_ref(:)
      ! tmp=tmp_ref
      !
      ub=1.78d2
      ! ub=680.d0!14.1d0
      if(.true. .and. present(rovd)) then 
        ! vel_ref(1)=ub*1.128d0*((1.d0-2.d0*rovd)**(1.d0/7.d0))
        h=0.5d0
        delta=0.05d0
        vel_ref(1)=100.d0+0.5d0*(ub-100.d0)*(1.d0-tanh(2.d0*(rovd-h)/delta))
        ! spc(spcindex('H2'))=spcinf(spcindex('H2'))+0.5d0*(spc_ref(spcindex('H2')) &
        !                 -spcinf(spcindex('H2')))*(1.d0-tanh(2.d0*(rovd-h)/delta))
        ! spc(spcindex('N2'))=1.d0-sum(spc(:))
        ! vel_ref(1)=vel_ref(1)*1.27767d0*((1.d0-2.d0*rovd)**(1.d0/3.d0))
        ! vel_ref(2:ndims)=0.d0*vel_ref(1)
        !
        tmp=tinf+0.5d0*(tmp_ref-tinf)*(1.d0-tanh(2.d0*(rovd-h)/delta))               
        vel(:)=vel_ref(:)
        ! if(nstep==0) vel(:)=whitenoise(vel_ref(:),vel(:),0.1d0)
      else 
        vel(1)=ub
      endif 
      !
    elseif(trim(stream)=='hotcoflow') then
      !
      tmp_ref=1250.d0
      spc_ref(spcindex('O2'))=0.24488d0
      spc_ref(spcindex('H2O'))=0.17491d0
      spc_ref(spcindex('N2'))=1.d0-sum(spc_ref(:))
      ! call convertxiyi(spcx(:),spc_ref(:),'X2Y')
      spc(:)=spc_ref(:)
      ! tmp=tmp_ref
      ! spc(spcindex('O2'))=0.233d0
      ! spc(spcindex('N2'))=1.d0-sum(spc(:))
      !
      ub=1.42d2
      ! ub=1.d0
      if(.true. .and. present(rovd)) then 
        h=0.5d0
        delta=0.05d0
        vel_ref(1)=100.d0+0.5d0*(ub-100.d0)*(1.d0-tanh(2.d0*(rovd-h)/delta))
        ! spc(spcindex('O2'))=spcinf(spcindex('O2'))+0.5d0*(spc_ref(spcindex('H2')) &
        !                 -spcinf(spcindex('H2')))*(1.d0-tanh(2.d0*(rovd-h)/delta))
        ! vel_ref(1)=ub*1.128d0*((1.d0-2.d0*rovd)**(1.d0/7.d0))
        ! vel_ref(1)=vel_ref(1)*1.27767d0*((1.d0-2.d0*rovd)**(1.d0/3.d0))
        ! vel_ref(2:ndims)=0.d0*vel_ref(1)
        !
        tmp=tinf+0.5d0*(tmp_ref-tinf)*(1.d0-tanh(2.d0*(rovd-h)/delta))
        vel(:)=vel_ref(:)
        ! if(nstep==0) vel(:)=whitenoise(vel_ref(:),vel(:),0.1d0)
      else 
        vel(1)=ub
      endif
      !
    elseif(trim(stream)=='air') then
      !
      tmp=tinf
      spc(:)=spcinf(:)
      !
      ub=100.d0
      vel(:)=0.d0
      vel(1)=ub
      !
    else
      stop ' !! stream type not set for the inlet boundary !!'
    endif
    !
    prs=pinf
    rho=thermal(pressure=prs,temperature=tmp,species=spc(:))
    !
    if(present(ubulk))ubulk=ub
    !
#endif
    !
  end subroutine multistream_inflow
  !+-------------------------------------------------------------------+
  !| The end of the subroutine multistream_inflow.                     |
  !+-------------------------------------------------------------------+
  !
  !!
end module fludyna
!+---------------------------------------------------------------------+
!| The end of the module fludyna.                                      |
!+---------------------------------------------------------------------+
