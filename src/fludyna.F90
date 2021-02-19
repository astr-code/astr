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
     module procedure thermal_3d
  end interface
  !
  interface fvar2q
     module procedure fvar2q_sca
     module procedure fvar2q_3da
  end interface
  !
  interface q2fvar
     module procedure q2fvar_sca
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
  subroutine fvar2q_3da(q,density,velocity,pressure,temperature,species)
    !
    use commvar, only: numq,ndims,num_species,const1,const6
    !
    real(8),intent(in) :: density(:,:,:),velocity(:,:,:,:)
    real(8),intent(in),optional :: pressure(:,:,:),temperature(:,:,:), &
                                   species(:,:,:,:)
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
  subroutine q2fvar_3da(q,density,velocity,pressure,temperature,species)
    !
    use commvar, only: numq,ndims,num_species,const1,const6
    !
    real(8),intent(in) :: q(:,:,:,:)
    real(8),intent(out) :: density(:,:,:),velocity(:,:,:,:),           &
                           pressure(:,:,:),temperature(:,:,:) 
    real(8),intent(out),optional :: species(:,:,:,:)
    !
    ! local data
    integer :: jspec
    integer :: dim(3)
    !
    density(:,:,:)   =q(:,:,:,1)
    !
    velocity(:,:,:,1)=q(:,:,:,2)/density
    velocity(:,:,:,2)=q(:,:,:,3)/density
    velocity(:,:,:,3)=q(:,:,:,4)/density
    !
    pressure(:,:,:)  =( q(:,:,:,5)-0.5d0*density(:,:,:)*(              &
                                         velocity(:,:,:,1)**2+         &
                                         velocity(:,:,:,2)**2+         &
                                         velocity(:,:,:,3)**2) )/const6
    ! !
    !
    dim(1)=size(q,1)
    dim(2)=size(q,2)
    dim(3)=size(q,3)
    !
    temperature=thermal(pressure=pressure,density=density,dim=dim)
    !
    if(num_species>0 .and. present(species)) then
      !
      do jspec=1,num_species
        species(:,:,:,jspec)=q(:,:,:,5+jspec)/density(:,:,:)
      enddo
      !
    endif
    !
  end subroutine q2fvar_3da
  !
  subroutine q2fvar_sca(q,density,velocity,pressure,temperature,species)
    !
    use commvar, only: numq,ndims,num_species,const1,const6
    !
    real(8),intent(in) :: q(:)
    real(8),intent(out) :: density,velocity(:),           &
                           pressure,temperature 
    real(8),intent(out),optional :: species(:)
    !
    ! local data
    integer :: jspec
    !
    density   =q(1)
    !
    velocity(1)=q(2)/density
    velocity(2)=q(3)/density
    velocity(3)=q(4)/density
    !
    pressure  =( q(5)-0.5d0*density*(velocity(1)**2+velocity(2)**2+    &
                                     velocity(3)**2) )/const6
    !
    temperature=thermal(pressure=pressure,density=density)
    !
    if(num_species>0 .and. present(species)) then
      !
      do jspec=1,num_species
        species(jspec)=q(5+jspec)/density
      enddo
      !
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
end module fludyna
!+---------------------------------------------------------------------+
!| The end of the module fludyna.                                      |
!+---------------------------------------------------------------------+
