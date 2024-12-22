module fluidynamcs
  
  implicit none
  
  contains
  
  function thermal_scar(density,pressure,temperature) result(vout)
    !
    use comvardef,only : const2
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
      stop ' !! unable to get thermal variable  @ thermal_scar !!'
    endif
    !
  end function thermal_scar
  
  function var2q(density,velocity,pressure,temperature) result(q)
    !
    use comvardef, only: const1,const6,numq
    !
    real(8) :: q(1:numq)
    real(8),intent(in) :: density,velocity(3)
    real(8),intent(in),optional :: pressure,temperature
    !
    ! local data
    real(8) :: var1
    !
    q(1)=density
    q(2)=density*velocity(1)
    q(3)=density*velocity(2)
    q(4)=density*velocity(3)
    !
    var1=0.5d0*sum(velocity(:)*velocity(:))
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
    real(8),intent(in) :: q(:)
    real(8),intent(out) :: density
    real(8),intent(out),optional :: velocity(:),pressure,temperature
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
      pressure  =( q(5)-0.5d0*density*(velocity(1)**2+velocity(2)**2+  &
                                       velocity(3)**2) )/const6
    endif
    !
    if(present(temperature)) then
      temperature=thermal_scar(pressure=pressure,density=density)
    endif

  end subroutine q2fvar
  
  pure real(8) function miucal(temper)
    !
    use comvardef, only :  ref_t
    !
    real(8),intent(in) :: temper
    ! temper represent temperature, dimensionless
    ! below calculate miucal using sutherland's law
    !
    real(8) :: tempconst,tempconst1
    ! 
    tempconst=110.4d0/ref_t
    tempconst1=1.d0+tempconst
    !
    miucal=temper*sqrt(temper)*tempconst1/(temper+tempconst)
    !
    return
    !
  end function miucal
  
end module fluidynamcs