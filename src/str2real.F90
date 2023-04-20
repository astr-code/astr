module str2real_m
  !
  implicit none
  private
  public :: str2real
  
  integer, parameter :: wp    = kind(1.0d0)
  integer, parameter :: ikind = selected_int_kind(2)
  integer(kind=ikind), parameter :: digit_0    = ichar('0')
  integer(kind=ikind), parameter :: period     = ichar('.') - digit_0
  integer(kind=ikind), parameter :: minus_sign = ichar('-') - digit_0
  
  real(wp), parameter :: whole_number_base(*) =                    &
          [1d15,  1d14,  1d13,  1d12,  1d11,  1d10,  1d9,   1d8,   &
           1d7,   1d6,   1d5,   1d4,   1d3,   1d2,   1d1,   1d0]
  real(wp), parameter :: fractional_base(*)   =                    &
          [1d-1,  1d-2,  1d-3,  1d-4,  1d-5,  1d-6,  1d-7,  1d-8,  &
           1d-9,  1d-10, 1d-11, 1d-12, 1d-13, 1d-14, 1d-15, 1d-16, &
           1d-17, 1d-18, 1d-19, 1d-20, 1d-21, 1d-22, 1d-23, 1d-24]
  real(wp), parameter :: period_skip = 0d0
  real(wp), parameter :: base(*) = &
          [whole_number_base, period_skip, fractional_base]
  
  contains
  
  
  
  function str2real(s) result(r)
  character(*), intent(in) :: s
  real(wp) :: r
  
  real(wp) :: r_coefficient
  real(wp) :: r_exponent
  
  integer, parameter :: N = 32
  character(N) :: factors_char
  integer(kind=ikind)    :: factors(N)
  integer(kind=ikind)    :: mask(N)
  integer      :: period_loc
  integer      :: exponent_loc
  integer      :: mask_from
  integer      :: mask_till
  
  equivalence(factors, factors_char)
  factors_char = s
  factors = factors - digit_0
  
  period_loc   = findloc(factors, period, 1)
  exponent_loc = scan(s, 'eE', back=.true.)
  if(exponent_loc == 0) exponent_loc = len(s) + 1
  if(period_loc   == 0) period_loc   = exponent_loc
  
  where (0 <= factors .and. factors <= 9)
      mask = 1
  elsewhere
      mask = 0
  end where
  
  mask_from = 18 - period_loc
  mask_till = mask_from + exponent_loc - 2
  
  r_coefficient = sum( &
          factors(:exponent_loc - 1)  * &
          base(mask_from:mask_till) * &
          mask(:exponent_loc - 1))
  r_exponent = sum( &
          factors(exponent_loc+1:len(s)) * &
          mask(exponent_loc+1:len(s))  * &
          base(17-(len(s)-exponent_loc):16))
  if(factors(exponent_loc+1) == minus_sign) r_exponent    = -r_exponent
  if(factors(1)              == minus_sign) r_coefficient = -r_coefficient
  r = r_coefficient * 10 ** r_exponent
  end function str2real



end module str2real_m