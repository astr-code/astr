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
  function thermal_scar(density,pressure,temperature) result(vout)
    !
    use commvar,only : const2
    !
    ! arguments
    real(8) :: vout
    real(8),intent(in) ,optional :: density,pressure,temperature
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
  end function thermal_scar
  !
  function thermal_1d(density,pressure,temperature,dim) result(vout)
    !
    use commvar,only : const2
    !
    ! arguments
    integer,intent(in) :: dim
    real(8) :: vout(dim)
    real(8),intent(in),optional :: density(:),pressure(:),        &
                                   temperature(:)
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
  end function thermal_1d
  !
  function thermal_3d(density,pressure,temperature,dim) result(vout)
    !
    use commvar,only : const2
    !
    ! arguments
    integer,intent(in) :: dim(3)
    real(8) :: vout(dim(1),dim(2),dim(3))
    real(8),intent(in),optional :: density(:,:,:),pressure(:,:,:),    &
                                    temperature(:,:,:)
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
    use commvar,  only : im,jm,km,num_species,num_modequ,turbmode
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
                      species=spc(0:im,0:jm,0:km,:)                    )
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
  subroutine fvar2q_sca(q,density,velocity,pressure,temperature,species)
    !
    use commvar, only: numq,ndims,num_species,const1,const6
    !
    real(8),intent(in) :: density,velocity(:)
    real(8),intent(in),optional :: pressure,temperature,species(:)
    real(8),intent(out) :: q(:)
    !
    ! local data
    integer :: jspec
    !
    q(1)=density
    q(2)=density*velocity(1)
    q(3)=density*velocity(2)
    q(4)=density*velocity(3)
    !
    if(present(temperature)) then
        q(5)=density*( temperature*const1 + 0.5d0*(velocity(1)**2 +    &
                                                   velocity(2)**2 +    &
                                                   velocity(3)**2) )
    elseif(present(pressure)) then
        q(5)=pressure*const6+0.5d0*density*( velocity(1)**2 +          &
                                             velocity(2)**2 +          &
                                             velocity(3)**2 )
    else
      print*,' !! pressure or temperature required !!'
      stop ' !! error @ fvar2q'
    endif
    !
    if(num_species>0) then
      !
      do jspec=1,num_species
        q(5+jspec)=density*species(jspec)
      enddo
      !
    endif
    !
  end subroutine fvar2q_sca
  !
  subroutine fvar2q_1da(q,density,velocity,pressure,temperature,species)
    !
    use commvar, only: numq,ndims,num_species,const1,const6
    !
    real(8),intent(in) :: density(:),velocity(:,:)
    real(8),intent(in),optional :: pressure(:),temperature(:), &
                                   species(:,:)
    real(8),intent(out) :: q(:,:)
    !
    ! local data
    integer :: jspec
    !
    q(:,1)=density(:)
    q(:,2)=density(:)*velocity(:,1)
    q(:,3)=density(:)*velocity(:,2)
    q(:,4)=density(:)*velocity(:,3)
    !
    if(present(temperature)) then
        q(:,5)=density(:)*( temperature(:)*const1 +        &
                             0.5d0*(velocity(:,1)**2 +             &
                                    velocity(:,2)**2 +             &
                                    velocity(:,3)**2) )
    elseif(present(pressure)) then
        q(:,5)=pressure(:)*const6+0.5d0*density(:)*(       &
                                    velocity(:,1)**2 +             &
                                    velocity(:,2)**2 +             &
                                    velocity(:,3)**2 )
    else
      print*,' !! pressure or temperature required !!'
      stop ' !! error @ fvar2q'
    endif
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
    use commvar, only: numq,ndims,num_species,const1,const6
    !
    real(8),intent(in) :: density(:,:,:),velocity(:,:,:,:)
    real(8),intent(in),optional :: pressure(:,:,:),temperature(:,:,:), &
                                   species(:,:,:,:),tke(:,:,:),        &
                                   omega(:,:,:)
    real(8),intent(out) :: q(:,:,:,:)
    !
    ! local data
    integer :: jspec
    !
    q(:,:,:,1)=density(:,:,:)
    q(:,:,:,2)=density(:,:,:)*velocity(:,:,:,1)
    q(:,:,:,3)=density(:,:,:)*velocity(:,:,:,2)
    q(:,:,:,4)=density(:,:,:)*velocity(:,:,:,3)
    !
    if(present(temperature)) then
        q(:,:,:,5)=density(:,:,:)*( temperature(:,:,:)*const1 +        &
                             0.5d0*(velocity(:,:,:,1)**2 +             &
                                    velocity(:,:,:,2)**2 +             &
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
    use commvar, only: numq,ndims,num_species,const1,const6
    !
    real(8),intent(in) :: q(:,:,:,:)
    real(8),intent(out) :: density(:,:,:)
    real(8),intent(out),optional :: temperature(:,:,:),pressure(:,:,:),&
                                    species(:,:,:,:),velocity(:,:,:,:),&
                                    tke(:,:,:),omega(:,:,:)
    !
    ! local data
    integer :: jspec
    integer :: dim(3)
    !
    density(:,:,:)   =q(:,:,:,1)
    !
    if(present(velocity) .or. present(pressure) .or. present(temperature)) then
      velocity(:,:,:,1)=q(:,:,:,2)/density
      velocity(:,:,:,2)=q(:,:,:,3)/density
      velocity(:,:,:,3)=q(:,:,:,4)/density
    endif
    !
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
    if(num_species>0 .and. present(species)) then
      !
      do jspec=1,num_species
        species(:,:,:,jspec)=q(:,:,:,5+jspec)/density(:,:,:)
      enddo
      !
    endif
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
    integer :: jspec
    integer :: dim
    !
    density(:)   =q(:,1)
    !
    if(present(velocity) .or. present(pressure) .or. present(temperature)) then
      velocity(:,1)=q(:,2)/density
      velocity(:,2)=q(:,3)/density
      velocity(:,3)=q(:,4)/density
    endif
    !
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
    !
    if(num_species>0 .and. present(species)) then
      !
      do jspec=1,num_species
        species(:,jspec)=q(:,5+jspec)/density(:)
      enddo
      !
    endif
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
    !
    density   =q(1)
    !
    if(present(velocity) .or. present(pressure) .or. present(temperature)) then
      velocity(1)=q(2)/density
      velocity(2)=q(3)/density
      velocity(3)=q(4)/density
    endif
    !
    if(present(pressure) .or. present(temperature)) then
      pressure  =( q(5)-0.5d0*density*(velocity(1)**2+velocity(2)**2+    &
                                       velocity(3)**2) )/const6
    endif
    !
    if(present(temperature)) then
      temperature=thermal(pressure=pressure,density=density)
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
    miucal=temper*sqrt(temper)*tempconst1/(temper+tempconst)
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
  pure real(8) function sos(tmp)
    !
    use commvar, only: mach
    !
    real(8),intent(in) :: tmp
    !
    sos=sqrt(tmp)/mach
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
  function mixinglayervel(y) result(u)
    !
    use constdef
    use commvar, only: uinf,ymin,ymax
    !
    ! arguments
    real(8) :: u(3)
    real(8),intent(in) :: y
    !
    ! local data
    real(8) :: u1,u2,delta
    !+-----------------------------------+
    !| delta: inflow vorticity thickness |
    !+-----------------------------------+
    u1=3.5d0
    u2=1.5d0
    delta=1.d0
    !
    u(1)=0.5d0*((u1+u2)+(u1-u2)*tanh(2.d0*y/delta))
    u(2)=0.d0
    u(3)=0.d0
    !
    return
    !
  end function mixinglayervel
  !+-------------------------------------------------------------------+
  !| The end of the function mixinglayervel.                           |
  !+-------------------------------------------------------------------+
  !
  !!
end module fludyna
!+---------------------------------------------------------------------+
!| The end of the module fludyna.                                      |
!+---------------------------------------------------------------------+
