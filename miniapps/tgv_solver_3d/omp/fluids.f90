module fluids

  use constdef
  
  implicit none
  
  contains
  
  function thermal_scar(density,pressure,temperature) result(vout)
    !
    use comvardef,only : const2
    !
    ! arguments
    real(rtype) :: vout
    real(rtype),intent(in) ,optional :: density,pressure,temperature
    !
    if(present(density) .and. present(temperature)) then
      vout=density*temperature/const2
    elseif(present(density) .and. present(pressure)) then
      vout=pressure/density*const2
    elseif(present(temperature) .and. present(pressure)) then
      vout=pressure/temperature*const2
    else
      stop ' !! unable to get thermal variable  @ thermal_scar !!'
    endif
    !
  end function thermal_scar
  
  function var2q(density,velocity,pressure,temperature) result(q)
    !
    use comvardef, only: const1,const6,numq
    !
    real(rtype) :: q(1:numq)
    real(rtype),intent(in) :: density,velocity(3)
    real(rtype),intent(in),optional :: pressure,temperature
    !
    ! local data
    real(rtype) :: var1
    !
    q(1)=density
    q(2)=density*velocity(1)
    q(3)=density*velocity(2)
    q(4)=density*velocity(3)
    !
    var1=0.5_rtype*sum(velocity(:)*velocity(:))
    if(present(temperature)) then
      q(5)=density*(temperature*const1+var1)
    elseif(present(pressure)) then
      q(5)=pressure*const6+density*var1
    endif
    !
  end function var2q
  
  subroutine q2fvar(q,density,velocity,pressure,temperature)
    !
    use comvardef, only: const6
    !
    real(rtype),intent(in) :: q(:)
    real(rtype),intent(out) :: density
    real(rtype),intent(out),optional :: velocity(:),pressure,temperature
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
      pressure  =( q(5)-0.5_rtype*density*(velocity(1)**2+velocity(2)**2+  &
                                                          velocity(3)**2) )/const6
    endif
    !
    if(present(temperature)) then
      temperature=thermal_scar(pressure=pressure,density=density)
    endif

  end subroutine q2fvar
  
  pure real(rtype) function miucal(temper)
    !
    use comvardef, only :  ref_t
    !
    real(rtype),intent(in) :: temper
    ! temper represent temperature, dimensionless
    ! below calculate miucal using sutherland's law
    !
    real(rtype) :: tempconst,tempconst1
    ! 
    tempconst=110.4_rtype/ref_t
    tempconst1=1._rtype+tempconst
    !
    miucal=temper*sqrt(temper)*tempconst1/(temper+tempconst)
    !
    return
    !
  end function miucal
  
end module fluids