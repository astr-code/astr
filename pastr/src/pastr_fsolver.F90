! fsolve, a Fortran90 code which solves systems of nonlinear equations, 
! inspired by the fsolve() function in minpack(), with special interfaces 
! fsolve_bdf2(), fsolve_bdf3, fsolve_be() and fsolve_tr() for handling 
! systems associated with implicit ODE solvers of type bdf2, bdf3, backward Euler, 
! midpoint, or trapezoidal.
! Author:
! Original Fortran77 version by Jorge More, Danny Sorenson, Burton Garbow, 
! Kenneth Hillstrom. This version by John Burkardt.
! 
! Reference:
! 
! Jorge More, Burton Garbow, Kenneth Hillstrom,
! User Guide for MINPACK-1,
! Technical Report ANL-80-74,
! Argonne National Laboratory, 1980.
! Jorge More, Danny Sorenson, Burton Garbow, Kenneth Hillstrom,
! The MINPACK Project,
! in Sources and Development of Mathematical Software,
! edited by Wayne Cowell,
! Prentice-Hall, 1984,
! ISBN: 0-13-823501-5,
! LC: QA76.95.S68.
! https://people.sc.fsu.edu/~jburkardt/f_src/fsolve/fsolve.html
module pastr_fsolver

    use iso_fortran_env, only: wp => real64

    implicit none

     
contains


    subroutine backward_euler_residual ( dydt, n, to, yo, tm, ym, fm )

    !*****************************************************************************80
    !
    !! backward_euler_residual() evaluates the backward Euler residual.
    !
    !  Discussion:
    !
    !    Let to and tm be two times, with yo and ym the associated ODE
    !    solution values there.  If ym satisfies the backward Euler condition,
    !    then
    !
    !      dydt(tm,ym) = ( ym - yo ) / ( tm - to )
    !
    !    This can be rewritten as
    !
    !      residual = ym - yo - ( tm - to ) * dydt(tm,ym)
    !
    !    Given the other information, a nonlinear equation solver can be used
    !    to estimate the value ym that makes the residual zero.
    !
    !  Licensing:
    !
    !    This code is distributed under the MIT license.
    !
    !  Modified:
    !
    !    08 November 2023
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Input:
    !
    !    external dydt(): the name of the user-supplied code which
    !    evaluates the right hand side of the ODE, of the form:
    !      subroutine dydt ( t, y, dy )
    !      real(wp) dy(*)
    !      real(wp) t
    !      real(wp) y(*)
    !
    !    integer n: the vector size.
    !
    !    real(wp) to, yo(n): the old time and solution.
    !
    !    real(wp) tm, ym(n): the new time and tentative solution.
    !
    !  Output:
    !
    !    real(wp) fm(n): the backward Euler residual.
    !
           
       
      integer n
    
      real(wp) dydtm(n)
      external dydt
      real(wp) fm(n)
      real(wp) tm
      real(wp) to
      real(wp) ym(n)
      real(wp) yo(n)
    
      call dydt ( tm, ym, dydtm )
    
      fm = ym - yo - ( tm - to ) * dydtm
    
      return
    end
    subroutine bdf2_residual ( dydt, n, t1, y1, t2, y2, t3, y3, fm )
    
    !*****************************************************************************80
    !
    !! bdf2_residual() evaluates the backward difference order 2 residual.
    !
    !  Discussion:
    !
    !    Let t1, t2 and t3 be three times, with y1, y2 and y3 the associated ODE
    !    solution values there.  Assume only the y3 value may be varied.
    !
    !    The BDF2 condition is:
    !
    !      w = ( t3 - t2 ) / ( t2 - t1 )
    !      b = ( 1 + w )^2 / ( 1 + 2 w )
    !      c = w^2 / ( 1 + 2 w )
    !      d = ( 1 + w ) / ( 1 + 2 w )
    !
    !      y3 - b y2 + c y1 = ( t3 - t2 ) * dydt( t3, y3 )
    !
    !    but if (t3-t2) = (t2-t1), we have:
    !
    !      w = 1
    !      b = 4/3
    !      c = 1/3
    !      d = 2/3
    !      y3 - 4/3 y2 + 1/3 y1 = 2 dt * dydt( t3, y3 )
    !
    !    This can be rewritten as
    !
    !      residual = y3 - b y2 + c y1 - ( t3 - t2 ) * dydt(t3,y3)
    !
    !    This is the BDF2 residual to be evaluated.
    !
    !  Licensing:
    !
    !    This code is distributed under the MIT license.
    !
    !  Modified:
    !
    !    17 November 2023
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Input:
    !
    !    external dydt(): the name of the user-supplied code which
    !    evaluates the right hand side of the ODE, of the form:
    !      subroutine dydt ( t, y, dy )
    !      real(wp) dy(*)
    !      real(wp) t
    !      real(wp) y(*)
    !
    !    integer n: the vector size.
    !
    !    real(wp) t1, y1(n), t2, y2(n), t3, y3(n): three sets of
    !    data at a sequence of times.
    !
    !  Output:
    !
    !    real(wp) fm(n): the residual.
    !
           
      integer n
    
      real(wp) b
      real(wp) c
      real(wp) d
      real(wp) dydt3(n)
      external dydt
      real(wp) fm(n)
      real(wp) t1
      real(wp) t2
      real(wp) t3
      real(wp) w
      real(wp) y1(n)
      real(wp) y2(n)
      real(wp) y3(n)
    
      w = ( t3 - t2 ) / ( t2 - t1 )
      b = ( 1.0D+00 + w )**2 / ( 1.0D+00 + 2.0D+00 * w )
      c = w**2 / ( 1.0D+00 + 2.0D+00 * w )
      d = ( 1.0D+00 + w ) / ( 1.0D+00 + 2.0D+00 * w )
    
      call dydt ( t3, y3, dydt3 )
    
      fm = y3 - b * y2 + c * y1 - d * ( t3 - t2 ) * dydt3
    
      return
    end
    subroutine bdf3_residual ( dydt, n, t1, y1, t2, y2, t3, y3, t4, y4, fm )
    
    !*****************************************************************************80
    !
    !! bdf3_residual() evaluates the backward difference order 3 residual.
    !
    !  Discussion:
    !
    !    We seek Y4 defined by the implicit equation:
    !
    !      11 Y4 - 18 Y3 + 9 Y2 - 2 Y1 = 6 * DT * F ( T4, Y4 )
    !
    !  Licensing:
    !
    !    This code is distributed under the MIT license.
    !
    !  Modified:
    !
    !    27 May 2025
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Input:
    !
    !    external dydt(): the name of the user-supplied code which
    !    evaluates the right hand side of the ODE, of the form:
    !      subroutine dydt ( t, y, dy )
    !      real(wp) dy(*)
    !      real(wp) t
    !      real(wp) y(*)
    !
    !    integer n: the vector size.
    !
    !    real(wp) t1, y1(n), t2, y2(n), t3, y3(n), t4, y4(n): 
    !    sets of data at a sequence of times.
    !
    !  Output:
    !
    !    real(wp) fm(n): the residual.
    !
           
      integer n
    
      real(wp) dydt4(n)
      external dydt
      real(wp) fm(n)
      real(wp) t1
      real(wp) t2
      real(wp) t3
      real(wp) t4
      real(wp) y1(n)
      real(wp) y2(n)
      real(wp) y3(n)
      real(wp) y4(n)
    
      ! call r8_fake_use ( t1 )
      ! call r8_fake_use ( t2 )
    
      call dydt ( t4, y4, dydt4 )
    
      fm = 11.0 * y4 - 18.0 * y3 + 9.0 * y2 - 2.0 * y1 - 6.0 * ( t4 - t3 ) * dydt4
    
      return
    end
    subroutine dogleg ( n, r, lr, diag, qtb, delta, x )
    
    !*****************************************************************************80
    !
    !! dogleg() finds the minimizing combination of Gauss-Newton and gradient steps.
    !
    !  Discussion:
    !
    !    Given an M by N matrix A, an N by N nonsingular diagonal
    !    matrix D, an M-vector B, and a positive number DELTA, the
    !    problem is to determine the convex combination X of the
    !    Gauss-Newton and scaled gradient directions that minimizes
    !    (A*X - B) in the least squares sense, subject to the
    !    restriction that the euclidean norm of D*X be at most DELTA.
    !
    !    This function completes the solution of the problem
    !    if it is provided with the necessary information from the
    !    QR factorization of A.  That is, if A = Q*R, where Q has
    !    orthogonal columns and R is an upper triangular matrix,
    !    then DOGLEG expects the full upper triangle of R and
    !    the first N components of Q'*B.
    !
    !  Licensing:
    !
    !    This code is distributed under the MIT license.
    !
    !  Modified:
    !
    !    06 April 2010
    !
    !  Author:
    !
    !    Original Fortran77 version by Jorge More, Burton Garbow, Kenneth Hillstrom.
    !    This version by John Burkardt.
    !
    !  Reference:
    !
    !    Jorge More, Burton Garbow, Kenneth Hillstrom,
    !    User Guide for MINPACK-1,
    !    Technical Report ANL-80-74,
    !    Argonne National Laboratory, 1980.
    !
    !  Input:
    !
    !    integer N, the order of the matrix R.
    !
    !    real(wp) R(LR), the upper triangular matrix R stored
    !    by rows.
    !
    !    integer LR, the size of the R array, which must be 
    !    no less than (N*(N+1))/2.
    !
    !    real(wp) DIAG(N), the diagonal elements of the matrix D.
    !
    !    real(wp) QTB(N), the first N elements of the vector Q'* B.
    !
    !    real(wp) DELTA, is a positive upper bound on the
    !    euclidean norm of D*X(1:N).
    !
    !  Output:
    !
    !    real(wp) X(N), the desired convex combination of the
    !    Gauss-Newton direction and the scaled gradient direction.
    !
           
      integer lr
      integer n
    
      real(wp) alpha
      real(wp) bnorm
      real(wp) delta
      real(wp) diag(n) 
      real(wp) epsmch
      real(wp) gnorm
      integer i
      integer j
      integer jj
      integer k
      integer l
      real(wp) qnorm
      real(wp) qtb(n)
      real(wp) r(lr)
      real(wp) sgnorm
      real(wp) sum2
      real(wp) temp
      real(wp) wa1(n)
      real(wp) wa2(n)
      real(wp) x(n)
    
      epsmch = epsilon ( epsmch )
    !
    !  Calculate the Gauss-Newton direction.
    !
      jj = ( n * ( n + 1 ) ) / 2 + 1
    
      do k = 1, n
    
         j = n - k + 1
         jj = jj - k
         l = jj + 1
         sum2 = 0.0D+00
    
         do i = j + 1, n
           sum2 = sum2 + r(l) * x(i)
           l = l + 1
         end do
    
         temp = r(jj)
    
         if ( temp == 0.0D+00 ) then
    
           l = j
           do i = 1, j
             temp = max ( temp, abs ( r(l)) )
             l = l + n - i
           end do
    
           if ( temp == 0.0D+00 ) then
             temp = epsmch
           else
             temp = epsmch * temp
           end if
    
         end if
    
         x(j) = ( qtb(j) - sum2 ) / temp
    
      end do
    !
    !  Test whether the Gauss-Newton direction is acceptable.
    !
      wa1(1:n) = 0.0D+00
      wa2(1:n) = diag(1:n) * x(1:n)
      qnorm = enorm ( n, wa2 )
    
      if ( qnorm <= delta ) then
        return
      end if
    !
    !  The Gauss-Newton direction is not acceptable.
    !  Calculate the scaled gradient direction.
    !
      l = 1
      do j = 1, n
         temp = qtb(j)
         do i = j, n
           wa1(i) = wa1(i) + r(l) * temp
           l = l + 1
         end do
         wa1(j) = wa1(j) / diag(j)
      end do
    !
    !  Calculate the norm of the scaled gradient.
    !  Test for the special case in which the scaled gradient is zero.
    !
      gnorm = enorm ( n, wa1 )
      sgnorm = 0.0D+00
      alpha = delta / qnorm
    
      if ( gnorm /= 0.0D+00 ) then
    !
    !  Calculate the point along the scaled gradient which minimizes the quadratic.
    !
        wa1(1:n) = ( wa1(1:n) / gnorm ) / diag(1:n)
    
        l = 1
        do j = 1, n
          sum2 = 0.0D+00
          do i = j, n
            sum2 = sum2 + r(l) * wa1(i)
            l = l + 1
          end do
          wa2(j) = sum2
        end do
    
        temp = enorm ( n, wa2 )
        sgnorm = ( gnorm / temp ) / temp
    !
    !  Test whether the scaled gradient direction is acceptable.
    !
        alpha = 0.0D+00
    !
    !  The scaled gradient direction is not acceptable.
    !  Calculate the point along the dogleg at which the quadratic is minimized.
    !
        if ( sgnorm < delta ) then
    
          bnorm = enorm ( n, qtb )
          temp = ( bnorm / gnorm ) * ( bnorm / qnorm ) * ( sgnorm / delta )
          temp = temp - ( delta / qnorm ) * ( sgnorm / delta) ** 2 &
            + sqrt ( ( temp - ( delta / qnorm ) ) ** 2 &
            + ( 1.0D+00 - ( delta / qnorm ) ** 2 ) &
            * ( 1.0D+00 - ( sgnorm / delta ) ** 2 ) )
    
          alpha = ( ( delta / qnorm ) * ( 1.0D+00 - ( sgnorm / delta ) ** 2 ) ) &
            / temp
    
        end if
    
      end if
    !
    !  Form appropriate convex combination of the Gauss-Newton
    !  direction and the scaled gradient direction.
    !
      temp = ( 1.0D+00 - alpha ) * min ( sgnorm, delta )
    
      x(1:n) = temp * wa1(1:n) + alpha * x(1:n)
    
      return
    end
    function enorm ( n, x )
    
    !*****************************************************************************80
    !
    !! enorm() computes the Euclidean norm of a vector.
    !
    !  Discussion:
    !
    !    The Euclidean norm is computed by accumulating the sum of
    !    squares in three different sums.  The sums of squares for the
    !    small and large components are scaled so that no overflows
    !    occur.  Non-destructive underflows are permitted.  Underflows
    !    and overflows do not occur in the computation of the unscaled
    !    sum of squares for the intermediate components.
    !
    !    The definitions of small, intermediate and large components
    !    depend on two constants, RDWARF and RGIANT.  The main
    !    restrictions on these constants are that RDWARF^2 not
    !    underflow and RGIANT^2 not overflow.
    !
    !  Licensing:
    !
    !    This code is distributed under the MIT license.
    !
    !  Modified:
    !
    !    06 April 2010
    !
    !  Author:
    !
    !    Original Fortran77 version by Jorge More, Burton Garbow, Kenneth Hillstrom.
    !    This version by John Burkardt.
    !
    !  Reference:
    !
    !    Jorge More, Burton Garbow, Kenneth Hillstrom,
    !    User Guide for MINPACK-1
    !    Argonne National Laboratory,
    !    Argonne, Illinois.
    !
    !  Input:
    !
    !    integer N, is the length of the vector.
    !
    !    real(wp) X(N), the vector whose norm is desired.
    !
    !  Output:
    !
    !    real(wp) ENORM, the Euclidean norm of the vector.
    !
           
       
      integer n
    
      real(wp) agiant 
      real(wp) enorm
      integer i
      real(wp) rdwarf
      real(wp) rgiant
      real(wp) s1
      real(wp) s2
      real(wp) s3
      real(wp) x(n)
      real(wp) xabs
      real(wp) x1max
      real(wp) x3max
    
      rdwarf = sqrt ( tiny ( rdwarf ) )
      rgiant = sqrt ( huge ( rgiant ) )
    
      s1 = 0.0D+00
      s2 = 0.0D+00
      s3 = 0.0D+00
      x1max = 0.0D+00
      x3max = 0.0D+00
      agiant = rgiant / real ( n, kind = wp )
    
      do i = 1, n
    
        xabs = abs ( x(i) )
    
        if ( xabs <= rdwarf ) then
    
          if ( x3max < xabs ) then
            s3 = 1.0D+00 + s3 * ( x3max / xabs ) ** 2
            x3max = xabs
          else if ( xabs /= 0.0D+00 ) then
            s3 = s3 + ( xabs / x3max ) ** 2
          end if
    
        else if ( agiant <= xabs ) then
    
          if ( x1max < xabs ) then
            s1 = 1.0D+00 + s1 * ( x1max / xabs ) ** 2
            x1max = xabs
          else
            s1 = s1 + ( xabs / x1max ) ** 2
          end if
    
        else
    
          s2 = s2 + xabs ** 2
    
        end if
    
      end do
    !
    !  Calculation of norm.
    !
      if ( s1 /= 0.0D+00 ) then
    
        enorm = x1max * sqrt ( s1 + ( s2 / x1max ) / x1max )
    
      else if ( s2 /= 0.0D+00 ) then
    
        if ( x3max <= s2 ) then
          enorm = sqrt ( s2 * ( 1.0D+00 + ( x3max / s2 ) * ( x3max * s3 ) ) )
        else
          enorm = sqrt ( x3max * ( ( s2 / x3max ) + ( x3max * s3 ) ) )
        end if
    
      else
    
        enorm = x3max * sqrt ( s3 )
    
      end if
    
      return
    end
    subroutine fdjac1 ( fcn, n, x, fvec, fjac, ldfjac, ml, mu, epsfcn )
    
    !*****************************************************************************80
    !
    !! fdjac1() estimates a jacobian matrix using forward differences.
    !
    !  Discussion:
    !
    !    This function computes a forward-difference approximation
    !    to the N by N jacobian matrix associated with a specified
    !    problem of N functions in N variables. If the jacobian has
    !    a banded form, then function evaluations are saved by only
    !    approximating the nonzero terms.
    !
    !  Licensing:
    !
    !    This code is distributed under the MIT license.
    !
    !  Modified:
    !
    !    06 April 2010
    !
    !  Author:
    !
    !    Original Fortran77 version by Jorge More, Burton Garbow, Kenneth Hillstrom.
    !    This version by John Burkardt.
    !
    !  Reference:
    !
    !    Jorge More, Burton Garbow, Kenneth Hillstrom,
    !    User Guide for MINPACK-1,
    !    Technical Report ANL-80-74,
    !    Argonne National Laboratory, 1980.
    !
    !  Input:
    !
    !    external FCN, the name of the user-supplied subroutine which
    !    calculates the functions.  The routine should have the form:
    !      subroutine fcn ( n, x, fvec )
    !      integer n
    !      real(wp) fvec(n)
    !      real(wp) x(n)
    !
    !    integer N, the number of functions and variables.
    !
    !    real(wp) X(N), the point where the jacobian is evaluated.
    !
    !    real(wp) FVEC(N), the functions evaluated at X.
    !
    !    integer LDFJAC, the leading dimension of FJAC, which
    !    must not be less than N.
    !
    !    integer ML, MU, specify the number of subdiagonals and
    !    superdiagonals within the band of the jacobian matrix.  If the
    !    jacobian is not banded, set ML and MU to N-1.
    !
    !    real(wp) EPSFCN, is used in determining a suitable step
    !    length for the forward-difference approximation.  This approximation
    !    assumes that the relative errors in the functions are of the order of
    !    EPSFCN.  If EPSFCN is less than the machine precision, it is assumed that
    !    the relative errors in the functions are of the order of the machine
    !    precision.
    !
    !  Output:
    !
    !    real(wp) FJAC(LDFJAC,N), the N by N approximate
    !    jacobian matrix.
    !
           
       
      integer ldfjac
      integer n
    
      real(wp) eps
      real(wp) epsfcn
      real(wp) epsmch
      external fcn
      real(wp) fjac(ldfjac,n)
      real(wp) fvec(n)
      real(wp) h
      integer i
      integer j
      integer k
      integer ml
      integer msum
      integer mu
      real(wp) temp
      real(wp) wa1(n)
      real(wp) wa2(n)
      real(wp) x(n)
    
      epsmch = epsilon ( epsmch )
    
      eps = sqrt ( max ( epsfcn, epsmch ) )
      msum = ml + mu + 1
    !
    !  Computation of dense approximate jacobian.
    !
      if ( n <= msum ) then
    
         do j = 1, n
    
            temp = x(j)
            h = eps * abs ( temp )
            if ( h == 0.0D+00 ) then
              h = eps
            end if
    
            x(j) = temp + h
            call fcn ( n, x, wa1 )
    
            x(j) = temp
            fjac(1:n,j) = ( wa1(1:n) - fvec(1:n) ) / h
    
         end do
    
      else
    !
    !  Computation of banded approximate jacobian.
    !
         do k = 1, msum
    
            do j = k, n, msum
              wa2(j) = x(j)
              h = eps * abs ( wa2(j) )
              if ( h == 0.0D+00 ) then
                h = eps
              end if
              x(j) = wa2(j) + h
            end do
    
            call fcn ( n, x, wa1 )
    
            do j = k, n, msum
    
              x(j) = wa2(j)
    
              h = eps * abs ( wa2(j) )
              if ( h == 0.0D+00 ) then
                h = eps
              end if
    
              fjac(1:n,j) = 0.0D+00
    
              do i = 1, n
                if ( j - mu <= i .and. i <= j + ml ) then
                  fjac(i,j) = ( wa1(i) - fvec(i) ) / h
                end if
              end do
    
            end do
    
         end do
    
      end if
    
      return
    end
    subroutine fdjac_bdf2 ( dydt, n, t1, x1, t2, x2, t3, x3, fvec, fjac, ldfjac, &
      ml, mu, epsfcn )
    
    !*****************************************************************************80
    !
    !! fdjac_bdf2() estimates a jacobian matrix using forward differences.
    !
    !  Discussion:
    !
    !    This function computes a forward-difference approximation
    !    to the N by N jacobian matrix associated with a specified
    !    problem of N functions in N variables.  If the jacobian has
    !    a banded form, then function evaluations are saved by only
    !    approximating the nonzero terms.
    !
    !  Licensing:
    !
    !    This code is distributed under the MIT license.
    !
    !  Modified:
    !
    !    17 November 2023
    !
    !  Author:
    !
    !    Original Fortran77 version by Jorge More, Burton Garbow, Kenneth Hillstrom.
    !    This version by John Burkardt.
    !
    !  Reference:
    !
    !    Jorge More, Burton Garbow, Kenneth Hillstrom,
    !    User Guide for MINPACK-1,
    !    Technical Report ANL-80-74,
    !    Argonne National Laboratory, 1980.
    !
    !  Input:
    !
    !    external dydt(): the name of the user-supplied code which
    !    evaluates the right hand side of the ODE, of the form:
    !      subroutine dydt ( t, y, dy )
    !      real(wp) dy(*)
    !      real(wp) t
    !      real(wp) y(*)
    !
    !    integer n: the number of functions and variables.
    !
    !    real(wp) t1, x1(n), t2, x2(n), t3, x3(n): 
    !    a sequence of three times and solution estimates.
    !   
    !    real(wp) fvec(n): the functions evaluated at x3.
    !
    !    integer ldfjac: the leading dimension of FJAC, which
    !    must not be less than N.
    !
    !    integer ml, mu: specify the number of subdiagonals and
    !    superdiagonals within the band of the jacobian matrix.  If the
    !    jacobian is not banded, set ML and MU to N-1.
    !
    !    real(wp) epsfcn: is used in determining a suitable step
    !    length for the forward-difference approximation.  This approximation
    !    assumes that the relative errors in the functions are of the order of
    !    EPSFCN.  If EPSFCN is less than the machine precision, it is assumed that
    !    the relative errors in the functions are of the order of the machine
    !    precision.
    !
    !  Output:
    !
    !    real(wp) fjac(ldfjac,n): the N by N approximate
    !    jacobian matrix.
    !
           
       
      integer ldfjac
      integer n
    
      external dydt
      real(wp) eps
      real(wp) epsfcn
      real(wp) epsmch
      real(wp) fjac(ldfjac,n)
      real(wp) fvec(n)
      real(wp) h
      integer i
      integer j
      integer k
      integer ml
      integer msum
      integer mu
      real(wp) t1
      real(wp) t2
      real(wp) t3
      real(wp) temp
      real(wp) wa1(n)
      real(wp) wa2(n)
      real(wp) x1(n)
      real(wp) x2(n)
      real(wp) x3(n)
    
      epsmch = epsilon ( epsmch )
    
      eps = sqrt ( max ( epsfcn, epsmch ) )
      msum = ml + mu + 1
    !
    !  Computation of dense approximate jacobian.
    !
      if ( n <= msum ) then
    
        do j = 1, n
    
          temp = x3(j)
          h = eps * abs ( temp )
          if ( h == 0.0D+00 ) then
            h = eps
          end if
    
          x3(j) = temp + h
          call bdf2_residual ( dydt, n, t1, x1, t2, x2, t3, x3, wa1 )
    
          x3(j) = temp
          fjac(1:n,j) = ( wa1(1:n) - fvec(1:n) ) / h
    
        end do
    !
    !  Computation of banded approximate jacobian.
    !
      else
    
        do k = 1, msum
    
          do j = k, n, msum
            wa2(j) = x3(j)
            h = eps * abs ( wa2(j) )
            if ( h == 0.0D+00 ) then
              h = eps
            end if
            x3(j) = wa2(j) + h
          end do
    
          call bdf2_residual ( dydt, n, t1, x1, t2, x2, t3, x3, wa1 )
    
          do j = k, n, msum
    
            x3(j) = wa2(j)
    
            h = eps * abs ( wa2(j) )
            if ( h == 0.0D+00 ) then
              h = eps
            end if
    
            fjac(1:n,j) = 0.0D+00
    
            do i = 1, n
              if ( j - mu <= i .and. i <= j + ml ) then
                fjac(i,j) = ( wa1(i) - fvec(i) ) / h
              end if
            end do
    
          end do
    
        end do
    
      end if
    
      return
    end
    subroutine fdjac_bdf3 ( dydt, n, t1, x1, t2, x2, t3, x3, t4, x4, fvec, fjac, &
      ldfjac, ml, mu, epsfcn )
    
    !*****************************************************************************80
    !
    !! fdjac_bdf3() estimates a jacobian matrix using forward differences.
    !
    !  Discussion:
    !
    !    This function computes a forward-difference approximation
    !    to the N by N jacobian matrix associated with a specified
    !    problem of N functions in N variables.  If the jacobian has
    !    a banded form, then function evaluations are saved by only
    !    approximating the nonzero terms.
    !
    !  Licensing:
    !
    !    This code is distributed under the MIT license.
    !
    !  Modified:
    !
    !    27 May 2025
    !
    !  Author:
    !
    !    Original Fortran77 version by Jorge More, Burton Garbow, Kenneth Hillstrom.
    !    This version by John Burkardt.
    !
    !  Reference:
    !
    !    Jorge More, Burton Garbow, Kenneth Hillstrom,
    !    User Guide for MINPACK-1,
    !    Technical Report ANL-80-74,
    !    Argonne National Laboratory, 1980.
    !
    !  Input:
    !
    !    external dydt(): the name of the user-supplied code which
    !    evaluates the right hand side of the ODE, of the form:
    !      subroutine dydt ( t, y, dy )
    !      real(wp) dy(*)
    !      real(wp) t
    !      real(wp) y(*)
    !
    !    integer n: the number of functions and variables.
    !
    !    real(wp) t1, x1(n), t2, x2(n), t3, x3(n), t4, x4(n): 
    !    a sequence of three times and solution estimates.
    !   
    !    real(wp) fvec(n): the functions evaluated at x3.
    !
    !    integer ldfjac: the leading dimension of FJAC, which
    !    must not be less than N.
    !
    !    integer ml, mu: specify the number of subdiagonals and
    !    superdiagonals within the band of the jacobian matrix.  If the
    !    jacobian is not banded, set ML and MU to N-1.
    !
    !    real(wp) epsfcn: is used in determining a suitable step
    !    length for the forward-difference approximation.  This approximation
    !    assumes that the relative errors in the functions are of the order of
    !    EPSFCN.  If EPSFCN is less than the machine precision, it is assumed that
    !    the relative errors in the functions are of the order of the machine
    !    precision.
    !
    !  Output:
    !
    !    real(wp) fjac(ldfjac,n): the N by N approximate
    !    jacobian matrix.
    !
           
       
      integer ldfjac
      integer n
    
      external dydt
      real(wp) eps
      real(wp) epsfcn
      real(wp) epsmch
      real(wp) fjac(ldfjac,n)
      real(wp) fvec(n)
      real(wp) h
      integer i
      integer j
      integer k
      integer ml
      integer msum
      integer mu
      real(wp) t1
      real(wp) t2
      real(wp) t3
      real(wp) t4
      real(wp) temp
      real(wp) wa1(n)
      real(wp) wa2(n)
      real(wp) x1(n)
      real(wp) x2(n)
      real(wp) x3(n)
      real(wp) x4(n)
    
      epsmch = epsilon ( epsmch )
    
      eps = sqrt ( max ( epsfcn, epsmch ) )
      msum = ml + mu + 1
    !
    !  Computation of dense approximate jacobian.
    !
      if ( n <= msum ) then
    
        do j = 1, n
    
          temp = x4(j)
          h = eps * abs ( temp )
          if ( h == 0.0D+00 ) then
            h = eps
          end if
    
          x4(j) = temp + h
          call bdf3_residual ( dydt, n, t1, x1, t2, x2, t3, x3, t4, x4, wa1 )
    
          x4(j) = temp
          fjac(1:n,j) = ( wa1(1:n) - fvec(1:n) ) / h
    
        end do
    !
    !  Computation of banded approximate jacobian.
    !
      else
    
        do k = 1, msum
    
          do j = k, n, msum
            wa2(j) = x4(j)
            h = eps * abs ( wa2(j) )
            if ( h == 0.0D+00 ) then
              h = eps
            end if
            x4(j) = wa2(j) + h
          end do
    
          call bdf3_residual ( dydt, n, t1, x1, t2, x2, t3, x3, t4, x4, wa1 )
    
          do j = k, n, msum
    
            x4(j) = wa2(j)
    
            h = eps * abs ( wa2(j) )
            if ( h == 0.0D+00 ) then
              h = eps
            end if
    
            fjac(1:n,j) = 0.0D+00
    
            do i = 1, n
              if ( j - mu <= i .and. i <= j + ml ) then
                fjac(i,j) = ( wa1(i) - fvec(i) ) / h
              end if
            end do
    
          end do
    
        end do
    
      end if
    
      return
    end
    subroutine fdjac_be ( dydt, n, to, xo, t, x, fvec, fjac, ldfjac, ml, &
      mu, epsfcn )
    
    !*****************************************************************************80
    !
    !! fdjac_be() estimates a jacobian matrix using forward differences.
    !
    !  Discussion:
    !
    !    This function computes a forward-difference approximation
    !    to the N by N jacobian matrix associated with a specified
    !    problem of N functions in N variables.  If the jacobian has
    !    a banded form, then function evaluations are saved by only
    !    approximating the nonzero terms.
    !
    !    The original code fdjac1() was modified to deal with problems
    !    involving a backward Euler residual.
    !
    !  Licensing:
    !
    !    This code is distributed under the MIT license.
    !
    !  Modified:
    !
    !    08 November 2023
    !
    !  Author:
    !
    !    Original Fortran77 version by Jorge More, Burton Garbow, Kenneth Hillstrom.
    !    This version by John Burkardt.
    !
    !  Reference:
    !
    !    Jorge More, Burton Garbow, Kenneth Hillstrom,
    !    User Guide for MINPACK-1,
    !    Technical Report ANL-80-74,
    !    Argonne National Laboratory, 1980.
    !
    !  Input:
    !
    !    external dydt(): the name of the user-supplied code which
    !    evaluates the right hand side of the ODE, of the form:
    !      subroutine dydt ( t, y, dy )
    !      real(wp) dy(*)
    !      real(wp) t
    !      real(wp) y(*)
    !
    !    integer n: the number of functions and variables.
    !
    !    real(wp) to, xo(n): the old time and solution.
    !   
    !    real(wp) t, x(n): the new time and current solution estimate. 
    !
    !    real(wp) fvec(n): the functions evaluated at X.
    !
    !    integer ldfjac: the leading dimension of FJAC, which
    !    must not be less than N.
    !
    !    integer ml, mu: specify the number of subdiagonals and
    !    superdiagonals within the band of the jacobian matrix.  If the
    !    jacobian is not banded, set ML and MU to N-1.
    !
    !    real(wp) epsfcn: is used in determining a suitable step
    !    length for the forward-difference approximation.  This approximation
    !    assumes that the relative errors in the functions are of the order of
    !    EPSFCN.  If EPSFCN is less than the machine precision, it is assumed that
    !    the relative errors in the functions are of the order of the machine
    !    precision.
    !
    !  Output:
    !
    !    real(wp) fjac(ldfjac,n): the N by N approximate
    !    jacobian matrix.
    !
           
       
      integer ldfjac
      integer n
    
      external dydt
      real(wp) eps
      real(wp) epsfcn
      real(wp) epsmch
      real(wp) fjac(ldfjac,n)
      real(wp) fvec(n)
      real(wp) h
      integer i
      integer j
      integer k
      integer ml
      integer msum
      integer mu
      real(wp) t
      real(wp) temp
      real(wp) to
      real(wp) wa1(n)
      real(wp) wa2(n)
      real(wp) x(n)
      real(wp) xo(n)
    
      epsmch = epsilon ( epsmch )
    
      eps = sqrt ( max ( epsfcn, epsmch ) )
      msum = ml + mu + 1
    !
    !  Computation of dense approximate jacobian.
    !
      if ( n <= msum ) then
    
        do j = 1, n
    
          temp = x(j)
          h = eps * abs ( temp )
          if ( h == 0.0D+00 ) then
            h = eps
          end if
    
          x(j) = temp + h
          call backward_euler_residual ( dydt, n, to, xo, t, x, wa1 )
    
          x(j) = temp
          fjac(1:n,j) = ( wa1(1:n) - fvec(1:n) ) / h
    
        end do
    !
    !  Computation of banded approximate jacobian.
    !
      else
    
        do k = 1, msum
    
          do j = k, n, msum
            wa2(j) = x(j)
            h = eps * abs ( wa2(j) )
            if ( h == 0.0D+00 ) then
              h = eps
            end if
            x(j) = wa2(j) + h
          end do
    
          call backward_euler_residual ( dydt, n, to, xo, t, x, wa1 )
    
          do j = k, n, msum
    
            x(j) = wa2(j)
    
            h = eps * abs ( wa2(j) )
            if ( h == 0.0D+00 ) then
              h = eps
            end if
    
            fjac(1:n,j) = 0.0D+00
    
            do i = 1, n
              if ( j - mu <= i .and. i <= j + ml ) then
                fjac(i,j) = ( wa1(i) - fvec(i) ) / h
              end if
            end do
    
          end do
    
        end do
    
      end if
    
      return
    end
    subroutine fdjac_tr ( dydt, n, to, xo, tn, xn, fvec, fjac, ldfjac, ml, &
      mu, epsfcn )
    
    !*****************************************************************************80
    !
    !! fdjac_tr() estimates a jacobian matrix using forward differences.
    !
    !  Discussion:
    !
    !    This function computes a forward-difference approximation
    !    to the N by N jacobian matrix associated with a specified
    !    problem of N functions in N variables.  If the jacobian has
    !    a banded form, then function evaluations are saved by only
    !    approximating the nonzero terms.
    !
    !    The original code fdjac1() was modified to deal with problems
    !    involving an implicit trapezoidal ODE residual.
    !
    !  Licensing:
    !
    !    This code is distributed under the MIT license.
    !
    !  Modified:
    !
    !    15 November 2023
    !
    !  Author:
    !
    !    Original Fortran77 version by Jorge More, Burton Garbow, Kenneth Hillstrom.
    !    This version by John Burkardt.
    !
    !  Reference:
    !
    !    Jorge More, Burton Garbow, Kenneth Hillstrom,
    !    User Guide for MINPACK-1,
    !    Technical Report ANL-80-74,
    !    Argonne National Laboratory, 1980.
    !
    !  Input:
    !
    !    external dydt(): the name of the user-supplied code which
    !    evaluates the right hand side of the ODE, of the form:
    !      subroutine dydt ( t, y, dy )
    !      real(wp) dy(*)
    !      real(wp) t
    !      real(wp) y(*)
    !
    !    integer n: the number of functions and variables.
    !
    !    real(wp) to, xo(n): the old time and solution.
    !   
    !    real(wp) tn, xn(n): the new time and current solution estimate. 
    !
    !    real(wp) fvec(n): the functions evaluated at X.
    !
    !    integer ldfjac: the leading dimension of FJAC, which
    !    must not be less than N.
    !
    !    integer ml, mu: specify the number of subdiagonals and
    !    superdiagonals within the band of the jacobian matrix.  If the
    !    jacobian is not banded, set ML and MU to N-1.
    !
    !    real(wp) epsfcn: is used in determining a suitable step
    !    length for the forward-difference approximation.  This approximation
    !    assumes that the relative errors in the functions are of the order of
    !    EPSFCN.  If EPSFCN is less than the machine precision, it is assumed that
    !    the relative errors in the functions are of the order of the machine
    !    precision.
    !
    !  Output:
    !
    !    real(wp) fjac(ldfjac,n): the N by N approximate
    !    jacobian matrix.
    !
           
       
      integer ldfjac
      integer n
    
      external dydt
      real(wp) eps
      real(wp) epsfcn
      real(wp) epsmch
      real(wp) fjac(ldfjac,n)
      real(wp) fvec(n)
      real(wp) h
      integer i
      integer j
      integer k
      integer ml
      integer msum
      integer mu
      real(wp) temp
      real(wp) tn
      real(wp) to
      real(wp) wa1(n)
      real(wp) wa2(n)
      real(wp) xn(n)
      real(wp) xo(n)
    
      epsmch = epsilon ( epsmch )
    
      eps = sqrt ( max ( epsfcn, epsmch ) )
      msum = ml + mu + 1
    !
    !  Computation of dense approximate jacobian.
    !
      if ( n <= msum ) then
    
        do j = 1, n
    
          temp = xn(j)
          h = eps * abs ( temp )
          if ( h == 0.0D+00 ) then
            h = eps
          end if
    
          xn(j) = temp + h
          call trapezoidal_residual ( dydt, n, to, xo, tn, xn, wa1 )
    
          xn(j) = temp
          fjac(1:n,j) = ( wa1(1:n) - fvec(1:n) ) / h
    
        end do
    !
    !  Computation of banded approximate jacobian.
    !
      else
    
        do k = 1, msum
    
          do j = k, n, msum
            wa2(j) = xn(j)
            h = eps * abs ( wa2(j) )
            if ( h == 0.0D+00 ) then
              h = eps
            end if
            xn(j) = wa2(j) + h
          end do
    
          call trapezoidal_residual ( dydt, n, to, xo, tn, xn, wa1 )
    
          do j = k, n, msum
    
            xn(j) = wa2(j)
    
            h = eps * abs ( wa2(j) )
            if ( h == 0.0D+00 ) then
              h = eps
            end if
    
            fjac(1:n,j) = 0.0D+00
    
            do i = 1, n
              if ( j - mu <= i .and. i <= j + ml ) then
                fjac(i,j) = ( wa1(i) - fvec(i) ) / h
              end if
            end do
    
          end do
    
        end do
    
      end if
    
      return
    end
    subroutine fsolve ( fcn, n, x, fvec, tol, info )
    
    !*****************************************************************************80
    !
    !! fsolve() seeks a zero of N nonlinear equations in N variables.
    !
    !  Discussion:
    !
    !    The code finds a zero of a system of N nonlinear functions in N variables
    !    by a modification of the Powell hybrid method.  This is done by using the
    !    more general nonlinear equation solver HYBRD.  The user provides a
    !    subroutine which calculates the functions.  
    !
    !    The jacobian is calculated by a forward-difference approximation.
    !
    !  Licensing:
    !
    !    This code is distributed under the MIT license.
    !
    !  Modified:
    !
    !    07 April 2021
    !
    !  Author:
    !
    !    Original Fortran77 version by Jorge More, Burton Garbow, Kenneth Hillstrom.
    !    This version by John Burkardt.
    !
    !  Reference:
    !
    !    Jorge More, Burton Garbow, Kenneth Hillstrom,
    !    User Guide for MINPACK-1,
    !    Technical Report ANL-80-74,
    !    Argonne National Laboratory, 1980.
    !
    !  Input:
    !
    !    external FCN, the user subroutine which calculates the functions.  
    !    The routine should have the form:
    !      subroutine fcn ( n, x, fvec )
    !      integer n
    !      real(wp) fvec(n)
    !      real(wp) x(n)
    !
    !    integer N, the number of functions and variables.
    !
    !    real(wp) X(N), an initial estimate of the solution vector.  
    !
    !    real(wp) TOL.  Satisfactory termination occurs when the algorithm
    !    estimates that the relative error between X and the solution is at
    !    most TOL.  TOL should be nonnegative.
    !
    !  Output:
    !
    !    real(wp) X(N), the estimate of the solution vector.
    !
    !    real(wp) FVEC(N), the functions evaluated at the output X.
    !
    !    integer INFO, error flag.
    !    0, improper input parameter values.
    !    1, algorithm estimates that the relative error between X and the
    !       solution is at most TOL.
    !    2, number of calls to FCN has reached or exceeded 200*(N+1).
    !    3, TOL is too small.  No further improvement in the approximate
    !       solution X is possible.
    !    4, the iteration is not making good progress.
    !
           
       
      integer n
    
      real(wp) diag(n)
      real(wp) epsfcn
      real(wp) factor
      external fcn
      real(wp) fjac(n,n)
      real(wp) fvec(n)
      integer info
      integer ldfjac
      integer lr
      integer maxfev
      integer ml
      integer mode
      integer mu
      integer nfev
      real(wp) qtf(n)
      real(wp) r((n*(n+1))/2)
      real(wp) tol
      real(wp) x(n)
      real(wp) xtol
    
      if ( n <= 0 ) then
        info = 0
        return
      end if
    
      if ( tol < 0.0D+00 ) then
        info = 0
        return
      end if
    
      xtol = tol
      maxfev = 200 * ( n + 1 )
      ml = n - 1
      mu = n - 1
      epsfcn = 0.0D+00
      diag(1:n) = 1.0D+00
      mode = 2
      factor = 100.0D+00
      info = 0
      nfev = 0
      fjac(1:n,1:n) = 0.0D+00
      ldfjac = n
      r(1:(n*(n+1))/2) = 0.0D+00
      lr = ( n * ( n + 1 ) ) / 2
      qtf(1:n) = 0.0D+00
    
      call hybrd ( fcn, n, x, fvec, xtol, maxfev, ml, mu, epsfcn, diag, mode, &
        factor, info, nfev, fjac, ldfjac, r, lr, qtf )
    
      if ( info == 5 ) then
        info = 4
      end if
    
      return
    end
    subroutine fsolve_bdf2 ( dydt, n, t1, x1, t2, x2, t3, x3, fvec, tol, info )
    
    !*****************************************************************************80
    !
    !! fsolve_bdf2() seeks a zero of N nonlinear equations in N variables.
    !
    !  Discussion:
    !
    !    The code finds a zero of a system of N nonlinear functions in N variables
    !    by a modification of the Powell hybrid method.  This is done by using the
    !    more general nonlinear equation solver hybrd_bdf2().  The user provides a
    !    subroutine which calculates the functions.  
    !
    !    The jacobian is calculated by a forward-difference approximation.
    !
    !  Licensing:
    !
    !    This code is distributed under the MIT license.
    !
    !  Modified:
    !
    !    17 November 2023
    !
    !  Author:
    !
    !    Original Fortran77 version by Jorge More, Burton Garbow, Kenneth Hillstrom.
    !    This version by John Burkardt.
    !
    !  Reference:
    !
    !    Jorge More, Burton Garbow, Kenneth Hillstrom,
    !    User Guide for MINPACK-1,
    !    Technical Report ANL-80-74,
    !    Argonne National Laboratory, 1980.
    !
    !  Input:
    !
    !    external dydt(), the name of the user-supplied code which
    !    evaluates the right hand side of the ODE, of the form:
    !      subroutine dydt ( t, y, dy )
    !      real(wp) dy(*)
    !      real(wp) t
    !      real(wp) y(*)
    !
    !    integer n: the number of functions and variables.
    !
    !    real(wp) t1, x1(n), t2, x2(n), t3, x3(n): 
    !    a sequence of three times and solution estimates.
    !
    !    real(wp) tol:  Satisfactory termination occurs when the algorithm
    !    estimates that the relative error between X and the solution is at
    !    most TOL.  TOL should be nonnegative.
    !
    !  Output:
    !
    !    real(wp) x3(n): the new estimate of the solution vector.
    !
    !    real(wp) fvec(n): the residuals evaluated at the output x.
    !
    !    integer info: a status flag.  A value of 1 represents success.
    !    0, improper input parameter values.
    !    1, algorithm estimates that the relative error between X and the
    !       solution is at most TOL.
    !    2, number of calls to derivative has reached or exceeded 200*(N+1).
    !    3, TOL is too small.  No further improvement in the approximate
    !       solution X is possible.
    !    4, the iteration is not making good progress.
    !
           
       
      integer n
    
      real(wp) diag(n)
      external dydt
      real(wp) epsfcn
      real(wp) factor
      real(wp) fjac(n,n)
      real(wp) fvec(n)
      integer info
      integer ldfjac
      integer lr
      integer maxfev
      integer ml
      integer mode
      integer mu
      integer nfev
      real(wp) qtf(n)
      real(wp) r((n*(n+1))/2)
      real(wp) t1
      real(wp) t2
      real(wp) t3
      real(wp) tol
      real(wp) x1(n)
      real(wp) x2(n)
      real(wp) x3(n)
      real(wp) xtol
    
      if ( n <= 0 ) then
        info = 0
        return
      end if
    
      if ( tol < 0.0D+00 ) then
        info = 0
        return
      end if
    
      xtol = tol
      maxfev = 200 * ( n + 1 )
      ml = n - 1
      mu = n - 1
      epsfcn = 0.0D+00
      diag(1:n) = 1.0D+00
      mode = 2
      factor = 100.0D+00
      info = 0
      nfev = 0
      fjac(1:n,1:n) = 0.0D+00
      ldfjac = n
      r(1:(n*(n+1))/2) = 0.0D+00
      lr = ( n * ( n + 1 ) ) / 2
      qtf(1:n) = 0.0D+00
    
      call hybrd_bdf2 ( dydt, n, t1, x1, t2, x2, t3, x3, fvec, xtol, maxfev, &
        ml, mu, epsfcn, diag, mode, factor, info, nfev, fjac, ldfjac, r, lr, qtf )
    
      if ( info == 5 ) then
        info = 4
      end if
    
      return
    end
    subroutine fsolve_bdf3 ( dydt, n, t1, x1, t2, x2, t3, x3, t4, x4, fvec, &
      tol, info )
    
    !*****************************************************************************80
    !
    !! fsolve_bdf3() seeks a zero of N nonlinear equations in N variables.
    !
    !  Discussion:
    !
    !    The code finds a zero of a system of N nonlinear functions in N variables
    !    by a modification of the Powell hybrid method.  This is done by using the
    !    more general nonlinear equation solver hybrd_bdf23().  The user provides a
    !    subroutine which calculates the functions.  
    !
    !    The jacobian is calculated by a forward-difference approximation.
    !
    !  Licensing:
    !
    !    This code is distributed under the MIT license.
    !
    !  Modified:
    !
    !    27 May 2025
    !
    !  Author:
    !
    !    Original Fortran77 version by Jorge More, Burton Garbow, Kenneth Hillstrom.
    !    This version by John Burkardt.
    !
    !  Reference:
    !
    !    Jorge More, Burton Garbow, Kenneth Hillstrom,
    !    User Guide for MINPACK-1,
    !    Technical Report ANL-80-74,
    !    Argonne National Laboratory, 1980.
    !
    !  Input:
    !
    !    external dydt(), the name of the user-supplied code which
    !    evaluates the right hand side of the ODE, of the form:
    !      subroutine dydt ( t, y, dy )
    !      real(wp) dy(*)
    !      real(wp) t
    !      real(wp) y(*)
    !
    !    integer n: the number of functions and variables.
    !
    !    real(wp) t1, x1(n), t2, x2(n), t3, x3(n), t4, x4(n): 
    !    a sequence of three times and solution estimates.
    !
    !    real(wp) tol:  Satisfactory termination occurs when the algorithm
    !    estimates that the relative error between X and the solution is at
    !    most TOL.  TOL should be nonnegative.
    !
    !  Output:
    !
    !    real(wp) x4(n): the new estimate of the solution vector.
    !
    !    real(wp) fvec(n): the residuals evaluated at the output x.
    !
    !    integer info: a status flag.  A value of 1 represents success.
    !    0, improper input parameter values.
    !    1, algorithm estimates that the relative error between X and the
    !       solution is at most TOL.
    !    2, number of calls to derivative has reached or exceeded 200*(N+1).
    !    3, TOL is too small.  No further improvement in the approximate
    !       solution X is possible.
    !    4, the iteration is not making good progress.
    !
           
       
      integer n
    
      real(wp) diag(n)
      external dydt
      real(wp) epsfcn
      real(wp) factor
      real(wp) fjac(n,n)
      real(wp) fvec(n)
      integer info
      integer ldfjac
      integer lr
      integer maxfev
      integer ml
      integer mode
      integer mu
      integer nfev
      real(wp) qtf(n)
      real(wp) r((n*(n+1))/2)
      real(wp) t1
      real(wp) t2
      real(wp) t3
      real(wp) t4
      real(wp) tol
      real(wp) x1(n)
      real(wp) x2(n)
      real(wp) x3(n)
      real(wp) x4(n)
      real(wp) xtol
    
      if ( n <= 0 ) then
        info = 0
        return
      end if
    
      if ( tol < 0.0D+00 ) then
        info = 0
        return
      end if
    
      xtol = tol
      maxfev = 200 * ( n + 1 )
      ml = n - 1
      mu = n - 1
      epsfcn = 0.0D+00
      diag(1:n) = 1.0D+00
      mode = 2
      factor = 100.0D+00
      info = 0
      nfev = 0
      fjac(1:n,1:n) = 0.0D+00
      ldfjac = n
      r(1:(n*(n+1))/2) = 0.0D+00
      lr = ( n * ( n + 1 ) ) / 2
      qtf(1:n) = 0.0D+00
    
      call hybrd_bdf3 ( dydt, n, t1, x1, t2, x2, t3, x3, t4, x4, fvec, xtol, &
        maxfev, ml, mu, epsfcn, diag, mode, factor, info, nfev, fjac, ldfjac, &
        r, lr, qtf )
    
      if ( info == 5 ) then
        info = 4
      end if
    
      return
    end
    subroutine fsolve_be ( dydt, n, to, xo, t, x, fvec, tol, info )
    
    !*****************************************************************************80
    !
    !! fsolve_be() seeks a zero of N nonlinear equations in N variables.
    !
    !  Discussion:
    !
    !    The code finds a zero of a system of N nonlinear functions in N variables
    !    by a modification of the Powell hybrid method.  This is done by using the
    !    more general nonlinear equation solver hybrd_be().  The user provides a
    !    subroutine which calculates the functions.  
    !
    !    The jacobian is calculated by a forward-difference approximation.
    !
    !    The original code fsolve() was modified to deal with problems
    !    involving a backward Euler residual.
    !
    !  Licensing:
    !
    !    This code is distributed under the MIT license.
    !
    !  Modified:
    !
    !    08 November 2023
    !
    !  Author:
    !
    !    Original Fortran77 version by Jorge More, Burton Garbow, Kenneth Hillstrom.
    !    This version by John Burkardt.
    !
    !  Reference:
    !
    !    Jorge More, Burton Garbow, Kenneth Hillstrom,
    !    User Guide for MINPACK-1,
    !    Technical Report ANL-80-74,
    !    Argonne National Laboratory, 1980.
    !
    !  Input:
    !
    !    external dydt(), the name of the user-supplied code which
    !    evaluates the right hand side of the ODE, of the form:
    !      subroutine dydt ( t, y, dy )
    !      real(wp) dy(*)
    !      real(wp) t
    !      real(wp) y(*)
    !
    !    integer n: the number of functions and variables.
    !
    !    real(wp) to, xo(n): the old time and solution.
    !   
    !    real(wp) t, x(n): the new time and current solution estimate. 
    !
    !    real(wp) tol:  Satisfactory termination occurs when the algorithm
    !    estimates that the relative error between X and the solution is at
    !    most TOL.  TOL should be nonnegative.
    !
    !  Output:
    !
    !    real(wp) x(n): the new estimate of the solution vector.
    !
    !    real(wp) fvec(n): the residuals evaluated at the output x.
    !
    !    integer info: a status flag.  A value of 1 represents success.
    !    0, improper input parameter values.
    !    1, algorithm estimates that the relative error between X and the
    !       solution is at most TOL.
    !    2, number of calls to derivative has reached or exceeded 200*(N+1).
    !    3, TOL is too small.  No further improvement in the approximate
    !       solution X is possible.
    !    4, the iteration is not making good progress.
    !
           
       
      integer n
    
      real(wp) diag(n)
      external dydt
      real(wp) epsfcn
      real(wp) factor
      real(wp) fjac(n,n)
      real(wp) fvec(n)
      integer info
      integer ldfjac
      integer lr
      integer maxfev
      integer ml
      integer mode
      integer mu
      integer nfev
      real(wp) qtf(n)
      real(wp) r((n*(n+1))/2)
      real(wp) t
      real(wp) to
      real(wp) tol
      real(wp) x(n)
      real(wp) xo(n)
      real(wp) xtol
    
      info = 1
    
      if ( n <= 0 ) then
        info = 0
        write ( *, '(a)' ) ''
        write ( *, '(a)' ) 'fsolve_be(): Fatal error!'
        write ( *, '(a)' ) '  n <= 0'
        return
      end if
    
      if ( tol < 0.0D+00 ) then
        info = 0
        write ( *, '(a)' ) ''
        write ( *, '(a)' ) 'fsolve_be(): Fatal error!'
        write ( *, '(a)' ) '  tol < 0'
        return
      end if

      xtol = tol
      maxfev = 200 * ( n + 1 )
      ml = n - 1
      mu = n - 1
      epsfcn = 0.0D+00
      diag(1:n) = 1.0D+00
      mode = 2
      factor = 100.0D+00
      info = 0
      nfev = 0
      fjac(1:n,1:n) = 0.0D+00
      ldfjac = n
      r(1:(n*(n+1))/2) = 0.0D+00
      lr = ( n * ( n + 1 ) ) / 2
      qtf(1:n) = 0.0D+00

      call hybrd_be ( dydt, n, to, xo, t, x, fvec, xtol, maxfev, ml, mu, &
        epsfcn, diag, mode, factor, info, nfev, fjac, ldfjac, r, lr, qtf )

      if ( info == 5 ) then
        info = 4
      end if

      return
    end
    subroutine fsolve_tr ( dydt, n, to, xo, tn, xn, fvec, tol, info )

    !*****************************************************************************80
    !
    !! fsolve_tr() seeks a zero of N nonlinear equations in N variables.
    !
    !  Discussion:
    !
    !    The code finds a zero of a system of N nonlinear functions in N variables
    !    by a modification of the Powell hybrid method.  This is done by using the
    !    more general nonlinear equation solver hybrd_be().  The user provides a
    !    subroutine which calculates the functions.  
    !
    !    The jacobian is calculated by a forward-difference approximation.
    !
    !    The original code fsolve() was modified to deal with problems
    !    involving an implicit trapezoidal ODE residual.
    !
    !  Licensing:
    !
    !    This code is distributed under the MIT license.
    !
    !  Modified:
    !
    !    15 November 2023
    !
    !  Author:
    !
    !    Original Fortran77 version by Jorge More, Burton Garbow, Kenneth Hillstrom.
    !    This version by John Burkardt.
    !
    !  Reference:
    !
    !    Jorge More, Burton Garbow, Kenneth Hillstrom,
    !    User Guide for MINPACK-1,
    !    Technical Report ANL-80-74,
    !    Argonne National Laboratory, 1980.
    !
    !  Input:
    !
    !    external dydt(), the name of the user-supplied code which
    !    evaluates the right hand side of the ODE, of the form:
    !      subroutine dydt ( t, y, dy )
    !      real(wp) dy(*)
    !      real(wp) t
    !      real(wp) y(*)
    !
    !    integer n: the number of functions and variables.
    !
    !    real(wp) to, xo(n): the old time and solution.
    !   
    !    real(wp) tn, xn(n): the new time and current solution estimate. 
    !
    !    real(wp) tol:  Satisfactory termination occurs when the algorithm
    !    estimates that the relative error between X and the solution is at
    !    most TOL.  TOL should be nonnegative.
    !
    !  Output:
    !
    !    real(wp) xn(n): the new estimate of the solution vector.
    !
    !    real(wp) fvec(n): the residuals evaluated at the output xn.
    !
    !    integer info: a status flag.  A value of 1 represents success.
    !    0, improper input parameter values.
    !    1, algorithm estimates that the relative error between X and the
    !       solution is at most TOL.
    !    2, number of calls to derivative has reached or exceeded 200*(N+1).
    !    3, TOL is too small.  No further improvement in the approximate
    !       solution X is possible.
    !    4, the iteration is not making good progress.
    !
       
      integer n

      real(wp) diag(n)
      external dydt
      real(wp) epsfcn
      real(wp) factor
      real(wp) fjac(n,n)
      real(wp) fvec(n)
      integer info
      integer ldfjac
      integer lr
      integer maxfev
      integer ml
      integer mode
      integer mu
      integer nfev
      real(wp) qtf(n)
      real(wp) r((n*(n+1))/2)
      real(wp) tn
      real(wp) to
      real(wp) tol
      real(wp) xn(n)
      real(wp) xo(n)
      real(wp) xtol

      if ( n <= 0 ) then
        info = 0
        return
      end if

      if ( tol < 0.0D+00 ) then
        info = 0
        return
      end if

      xtol = tol
      maxfev = 200 * ( n + 1 )
      ml = n - 1
      mu = n - 1
      epsfcn = 0.0D+00
      diag(1:n) = 1.0D+00
      mode = 2
      factor = 100.0D+00
      info = 0
      nfev = 0
      fjac(1:n,1:n) = 0.0D+00
      ldfjac = n
      r(1:(n*(n+1))/2) = 0.0D+00
      lr = ( n * ( n + 1 ) ) / 2
      qtf(1:n) = 0.0D+00

      call hybrd_tr ( dydt, n, to, xo, tn, xn, fvec, xtol, maxfev, ml, mu, &
        epsfcn, diag, mode, factor, info, nfev, fjac, ldfjac, r, lr, qtf )

      if ( info == 5 ) then
        info = 4
      end if

      return
    end
    subroutine hybrd ( fcn, n, x, fvec, xtol, maxfev, ml, mu, epsfcn, diag, mode, &
      factor, info, nfev, fjac, ldfjac, r, lr, qtf )

    !*****************************************************************************80
    !
    !! hybrd() seeks a zero of N nonlinear equations in N variables.
    !
    !  Discussion:
    !
    !    The code finds a zero of a system of N nonlinear functions in N variables
    !    by a modification of the Powell hybrid method.  The user must provide a
    !    subroutine which calculates the functions.  
    !
    !    The jacobian is then calculated by a forward-difference approximation.
    !
    !  Licensing:
    !
    !    This code is distributed under the MIT license.
    !
    !  Modified:
    !
    !    06 April 2010
    !
    !  Author:
    !
    !    Original Fortran77 version by Jorge More, Burton Garbow, Kenneth Hillstrom.
    !    This version by John Burkardt.
    !
    !  Reference:
    !
    !    Jorge More, Burton Garbow, Kenneth Hillstrom,
    !    User Guide for MINPACK-1,
    !    Technical Report ANL-80-74,
    !    Argonne National Laboratory, 1980.
    !
    !  Input:
    !
    !    external FCN, the name of the user-supplied subroutine which
    !    calculates the functions.  The routine should have the form:
    !      subroutine fcn ( n, x, fvec )
    !      integer n
    !      real(wp) fvec(n)
    !      real(wp) x(n)
    !
    !    integer N, the number of functions and variables.
    !
    !    real(wp) X(N): an initial estimate of the solution vector.  
    !
    !    real(wp) XTOL.  Termination occurs when the relative error
    !    between two consecutive iterates is at most XTOL.  XTOL should be
    !    nonnegative.
    !
    !    integer MAXFEV.  Termination occurs when the number of
    !    calls to FCN is at least MAXFEV by the end of an iteration.
    !
    !    integer ML, MU, specify the number of subdiagonals and
    !    superdiagonals within the band of the jacobian matrix.  If the jacobian
    !    is not banded, set ML and MU to at least n - 1.
    !
    !    real(wp) EPSFCN, is used in determining a suitable step
    !    length for the forward-difference approximation.  This approximation
    !    assumes that the relative errors in the functions are of the order of
    !    EPSFCN.  If EPSFCN is less than the machine precision, it is assumed that
    !    the relative errors in the functions are of the order of the machine
    !    precision.
    !
    !    real(wp) DIAG(N).  If MODE = 2, then DIAG must contain 
    !    positive entries that serve as multiplicative scale factors 
    !    for the variables.
    !
    !    integer MODE, scaling option.
    !    1, variables will be scaled internally.
    !    2, scaling is specified by the input DIAG vector.
    !
    !    real(wp) FACTOR, determines the initial step bound.  This
    !    bound is set to the product of FACTOR and the euclidean norm of DIAG*X if
    !    nonzero, or else to FACTOR itself.  In most cases, FACTOR should lie
    !    in the interval (0.1, 100) with 100 the recommended value.
    !
    !    integer LDFJAC, the leading dimension of FJAC.
    !    LDFJAC must be at least N.
    !
    !    integer LR, the size of the R array, which must be no
    !    less than (N*(N+1))/2.
    !
    !  Output:
    !
    !    real(wp) X(N): the final estimate of the solution vector.
    !
    !    real(wp) FVEC(N), the functions evaluated at the output X.
    !
    !    real(wp) DIAG(N): positive entries that
    !    served as multiplicative scale factors for the variables.
    !
    !    integer INFO, error flag. 
    !    0, improper input parameter values.
    !    1, relative error between two consecutive iterates is at most XTOL.
    !    2, number of calls to FCN has reached or exceeded MAXFEV.
    !    3, XTOL is too small.  No further improvement in the approximate
    !       solution X is possible.
    !    4, iteration is not making good progress, as measured by the improvement
    !       from the last five jacobian evaluations.
    !    5, iteration is not making good progress, as measured by the improvement
    !       from the last ten iterations.
    !
    !    integer NFEV, the number of calls to FCN.
    !
    !    real(wp) FJAC(LDFJAC,N), an N by N array which contains
    !    the orthogonal matrix Q produced by the QR factorization of the final
    !    approximate jacobian.
    !
    !    real(wp) R(LR), the upper triangular matrix produced by
    !    the QR factorization of the final approximate jacobian, stored rowwise.
    !
    !    real(wp) QTF(N), contains the vector Q'*FVEC.
    !
       
      integer ldfjac
      integer lr
      integer n

      real(wp) actred
      real(wp) delta
      real(wp) diag(n) 
      real(wp) epsfcn
      real(wp) epsmch
      real(wp) factor
      external fcn
      real(wp) fjac(ldfjac,n)
      real(wp) fnorm
      real(wp) fnorm1
      real(wp) fvec(n)
      integer i
      integer info
      integer iter
      integer iwa(1)
      integer j
      logical jeval
      integer l
      integer maxfev
      integer ml
      integer mode
      integer msum
      integer mu
      integer ncfail
      integer nslow1
      integer nslow2
      integer ncsuc
      integer nfev
      logical pivot
      real(wp) pnorm
      real(wp) prered
      real(wp) qtf(n)
      real(wp) r(lr)
      real(wp) ratio
      logical sing
      real(wp) sum2
      real(wp) temp
      real(wp) wa1(n)
      real(wp) wa2(n)
      real(wp) wa3(n)
      real(wp) wa4(n)
      real(wp) x(n)
      real(wp) xnorm
      real(wp) xtol

      epsmch = epsilon ( epsmch )

      info = 0
      nfev = 0
    !
    !  Check the input for errors.
    !
      if ( n <= 0 ) then
        return
      else if ( xtol < 0.0D+00 ) then
        return
      else if ( maxfev <= 0 ) then
        return
      else if ( ml < 0 ) then
        return
      else if ( mu < 0 ) then
        return
      else if ( factor <= 0.0D+00 ) then
        return
      else if ( ldfjac < n ) then
        return
      else if ( lr < ( n * ( n + 1 ) ) / 2 ) then
        return
      end if

      if ( mode == 2 ) then

        do j = 1, n
          if ( diag(j) <= 0.0D+00 ) then
            return
          end if
        end do

      end if
    !
    !  Evaluate the function at the starting point
    !  and calculate its norm.
    !
      call fcn ( n, x, fvec )
      nfev = 1

      fnorm = enorm ( n, fvec )
    !
    !  Determine the number of calls to FCN needed to compute the jacobian matrix.
    !
      msum = min ( ml + mu + 1, n )
    !
    !  Initialize iteration counter and monitors.
    !
      iter = 1
      ncsuc = 0
      ncfail = 0
      nslow1 = 0
      nslow2 = 0
    !
    !  Beginning of the outer loop.
    !
    30 continue

        jeval = .true.
    !
    !  Calculate the jacobian matrix.
    !
        call fdjac1 ( fcn, n, x, fvec, fjac, ldfjac, ml, mu, epsfcn )

        nfev = nfev + msum
    !
    !  Compute the QR factorization of the jacobian.
    !
        pivot = .false.
        call qrfac ( n, n, fjac, ldfjac, pivot, iwa, 1, wa1, wa2 )
    !
    !  On the first iteration, if MODE is 1, scale according
    !  to the norms of the columns of the initial jacobian.
    !
        if ( iter == 1 ) then

          if ( mode /= 2 ) then

            diag(1:n) = wa2(1:n)
            do j = 1, n
              if ( wa2(j) == 0.0D+00 ) then
                diag(j) = 1.0D+00
              end if
            end do

          end if
    !
    !  On the first iteration, calculate the norm of the scaled X
    !  and initialize the step bound DELTA.
    !
          wa3(1:n) = diag(1:n) * x(1:n)
          xnorm = enorm ( n, wa3 )
          delta = factor * xnorm
          if ( delta == 0.0D+00 ) then
            delta = factor
          end if

        end if
    !
    !  Form Q' * FVEC and store in QTF.
    !
         qtf(1:n) = fvec(1:n)

         do j = 1, n

           if ( fjac(j,j) /= 0.0D+00 ) then
             temp = - dot_product ( qtf(j:n), fjac(j:n,j) ) / fjac(j,j)
             qtf(j:n) = qtf(j:n) + fjac(j:n,j) * temp
           end if

         end do
    !
    !  Copy the triangular factor of the QR factorization into R.
    !
         sing = .false.

         do j = 1, n
            l = j
            do i = 1, j - 1
              r(l) = fjac(i,j)
              l = l + n - i
            end do
            r(l) = wa1(j)
            if ( wa1(j) == 0.0D+00 ) then
              sing = .true.
            end if
         end do
    !
    !  Accumulate the orthogonal factor in FJAC.
    !
         call qform ( n, n, fjac, ldfjac )
    !
    !  Rescale if necessary.
    !
         if ( mode /= 2 ) then
           do j = 1, n
             diag(j) = max ( diag(j), wa2(j) )
           end do
         end if
    !
    !  Beginning of the inner loop.
    !
    180    continue
    !
    !  Determine the direction P.
    !
            call dogleg ( n, r, lr, diag, qtf, delta, wa1 )
    !
    !  Store the direction P and X + P.
    !  Calculate the norm of P.
    !
            wa1(1:n) = - wa1(1:n)
            wa2(1:n) = x(1:n) + wa1(1:n)
            wa3(1:n) = diag(1:n) * wa1(1:n)

            pnorm = enorm ( n, wa3 )
    !
    !  On the first iteration, adjust the initial step bound.
    !
            if ( iter == 1 ) then
              delta = min ( delta, pnorm )
            end if
    !
    !  Evaluate the function at X + P and calculate its norm.
    !
            call fcn ( n, wa2, wa4 )
            nfev = nfev + 1
            fnorm1 = enorm ( n, wa4 )
    !
    !  Compute the scaled actual reduction.
    !
            actred = -1.0D+00
            if ( fnorm1 < fnorm ) then
              actred = 1.0D+00 - ( fnorm1 / fnorm ) ** 2
            endif
    !
    !  Compute the scaled predicted reduction.
    !
            l = 1
            do i = 1, n
              sum2 = 0.0D+00
              do j = i, n
                sum2 = sum2 + r(l) * wa1(j)
                l = l + 1
              end do
              wa3(i) = qtf(i) + sum2
            end do

            temp = enorm ( n, wa3 )
            prered = 0.0D+00
            if ( temp < fnorm ) then
              prered = 1.0D+00 - ( temp / fnorm ) ** 2
            end if
    !
    !  Compute the ratio of the actual to the predicted reduction.
    !
            ratio = 0.0D+00
            if ( 0.0D+00 < prered ) then
              ratio = actred / prered
            end if
    !
    !  Update the step bound.
    !
            if ( ratio < 0.1D+00 ) then

              ncsuc = 0
              ncfail = ncfail + 1
              delta = 0.5D+00 * delta

            else

              ncfail = 0
              ncsuc = ncsuc + 1

              if ( 0.5D+00 <= ratio .or. 1 < ncsuc ) then
                delta = max ( delta, pnorm / 0.5D+00 )
              end if

              if ( abs ( ratio - 1.0D+00 ) <= 0.1D+00 ) then
                delta = pnorm / 0.5D+00
              end if

            end if
    !
    !  Test for successful iteration.
    !
    !  Successful iteration.
    !  Update X, FVEC, and their norms.
    !
            if ( 0.0001D+00 <= ratio ) then
              x(1:n) = wa2(1:n)
              wa2(1:n) = diag(1:n) * x(1:n)
              fvec(1:n) = wa4(1:n)
              xnorm = enorm ( n, wa2 )
              fnorm = fnorm1
              iter = iter + 1
            end if
    !
    !  Determine the progress of the iteration.
    !
            nslow1 = nslow1 + 1
            if ( 0.001D+00 <= actred ) then
              nslow1 = 0
            end if

            if ( jeval ) then
              nslow2 = nslow2 + 1
            end if

            if ( 0.1D+00 <= actred ) then
              nslow2 = 0
            end if
    !
    !  Test for convergence.
    !
            if ( delta <= xtol * xnorm .or. fnorm == 0.0D+00 ) then
              info = 1
            end if

            if ( info /= 0 ) then
              return
            end if
    !
    !  Tests for termination and stringent tolerances.
    !
            if ( maxfev <= nfev ) then
              info = 2
            end if

            if ( 0.1D+00 * max ( 0.1D+00 * delta, pnorm ) <= epsmch * xnorm ) then
              info = 3
            end if

            if ( nslow2 == 5 ) then
              info = 4
            end if

            if ( nslow1 == 10 ) then
              info = 5
            end if

            if ( info /= 0 ) then
              return
            end if
    !
    !  Criterion for recalculating jacobian approximation
    !  by forward differences.
    !
            if ( ncfail == 2 ) then
              go to 290
            end if
    !
    !  Calculate the rank one modification to the jacobian
    !  and update QTF if necessary.
    !
            do j = 1, n
              sum2 = dot_product ( wa4(1:n), fjac(1:n,j) )
              wa2(j) = ( sum2 - wa3(j) ) / pnorm
              wa1(j) = diag(j) * ( ( diag(j) * wa1(j) ) / pnorm )
              if ( 0.0001D+00 <= ratio ) then
                qtf(j) = sum2
              end if
            end do
    !
    !  Compute the QR factorization of the updated jacobian.
    !
            call r1updt ( n, n, r, lr, wa1, wa2, wa3, sing )
            call r1mpyq ( n, n, fjac, ldfjac, wa2, wa3 )
            call r1mpyq ( 1, n, qtf, 1, wa2, wa3 )
    !
    !  End of the inner loop.
    !
            jeval = .false.
            go to 180

      290   continue
    !
    !  End of the outer loop.
    !
         go to 30

      return
    end
    subroutine hybrd_bdf2 ( dydt, n, t1, x1, t2, x2, t3, x3, fvec, xtol, maxfev, &
      ml, mu, epsfcn, diag, mode, factor, info, nfev, fjac, ldfjac, r, lr, qtf )

    !*****************************************************************************80
    !
    !! hybrd_bdf2() seeks a zero of N nonlinear equations in N variables.
    !
    !  Discussion:
    !
    !    The code finds a zero of a system of N nonlinear functions in N variables
    !    by a modification of the Powell hybrid method.  The user must provide a
    !    subroutine which calculates the functions.  
    !
    !    The jacobian is then calculated by a forward-difference approximation.
    !
    !  Licensing:
    !
    !    This code is distributed under the MIT license.
    !
    !  Modified:
    !
    !    17 November 2023
    !
    !  Author:
    !
    !    Original Fortran77 version by Jorge More, Burton Garbow, Kenneth Hillstrom.
    !    This version by John Burkardt.
    !
    !  Reference:
    !
    !    Jorge More, Burton Garbow, Kenneth Hillstrom,
    !    User Guide for MINPACK-1,
    !    Technical Report ANL-80-74,
    !    Argonne National Laboratory, 1980.
    !
    !  Input:
    !
    !    external dydt(), the name of the user-supplied code which
    !    evaluates the right hand side of the ODE, of the form:
    !      subroutine dydt ( t, y, dy )
    !      real(wp) dy(n)
    !      real(wp) t
    !      real(wp) y(n)
    !
    !    integer n: the number of functions and variables.
    !
    !    real(wp) t1, x1(n), t2, x2(n), t3, x3(n): 
    !    a sequence of three times and solution estimates.
    !
    !    real(wp) xtol: Termination occurs when the relative error
    !    between two consecutive iterates is at most XTOL.  XTOL should be
    !    nonnegative.
    !
    !    integer maxfev: Termination occurs when the number of
    !    calls to the derivative code is at least MAXFEV by the end of an iteration.
    !
    !    integer ml, mu: specify the number of subdiagonals and
    !    superdiagonals within the band of the jacobian matrix.  If the jacobian
    !    is not banded, set ML and MU to at least n - 1.
    !
    !    real(wp) epsfcn: is used in determining a suitable step
    !    length for the forward-difference approximation.  This approximation
    !    assumes that the relative errors in the functions are of the order of
    !    EPSFCN.  If EPSFCN is less than the machine precision, it is assumed that
    !    the relative errors in the functions are of the order of the machine
    !    precision.
    !
    !    real(wp) diag(n): multiplicative scale factors for the 
    !    variables.  Only needed as input if MODE = 2.  
    !
    !    integer mode: scaling option.
    !    1, variables will be scaled internally.
    !    2, scaling is specified by the input DIAG vector.
    !
    !    real(wp) factor: determines the initial step bound.  This
    !    bound is set to the product of FACTOR and the euclidean norm of DIAG*X if
    !    nonzero, or else to FACTOR itself.  In most cases, FACTOR should lie
    !    in the interval (0.1, 100) with 100 the recommended value.
    !
    !    integer ldfjac: the leading dimension of FJAC.
    !    LDFJAC must be at least N.
    !
    !    integer lr: the size of the R array, which must be no
    !    less than (N*(N+1))/2.
    !
    !  Output:
    !
    !    real(wp) x3(n): the final estimate of the solution vector.
    !
    !    real(wp) fvec(n): the functions evaluated at the output x.
    !
    !    real(wp) diag(n): multiplicative scale factors for the 
    !    variables.  Set internally if MODE = 1.
    !    
    !    integer info: status flag.  A value of 1 indicates success. 
    !    0, improper input parameter values.
    !    1, relative error between two consecutive iterates is at most XTOL.
    !    2, number of calls to derivative has reached or exceeded MAXFEV.
    !    3, XTOL is too small.  No further improvement in the approximate
    !       solution X is possible.
    !    4, iteration is not making good progress, as measured by the improvement
    !       from the last five jacobian evaluations.
    !    5, iteration is not making good progress, as measured by the improvement
    !       from the last ten iterations.
    !
    !    integer nfev: the number of calls to the derivative function.
    !
    !    real(wp) fjac(ldfjac,n): an N by N array which contains
    !    the orthogonal matrix Q produced by the QR factorization of the final
    !    approximate jacobian.
    !
    !    real(wp) r(lr): the upper triangular matrix produced by
    !    the QR factorization of the final approximate jacobian, stored rowwise.
    !
    !    real(wp) qtf(n): contains the vector Q'*FVEC.
    !
       
      integer ldfjac
      integer lr
      integer n

      real(wp) actred
      real(wp) delta
      real(wp) diag(n)
      external dydt 
      real(wp) epsfcn
      real(wp) epsmch
      real(wp) factor
      real(wp) fjac(ldfjac,n)
      real(wp) fnorm
      real(wp) fnorm1
      real(wp) fvec(n)
      integer i
      integer info
      integer iter
      integer iwa(1)
      integer j
      logical jeval
      integer l
      integer maxfev
      integer ml
      integer mode
      integer msum
      integer mu
      integer ncfail
      integer ncsuc
      integer nfev
      integer nslow1
      integer nslow2
      logical pivot
      real(wp) pnorm
      real(wp) prered
      real(wp) qtf(n)
      real(wp) r(lr)
      real(wp) ratio
      logical sing
      real(wp) sum2
      real(wp) t1
      real(wp) t2
      real(wp) t3
      real(wp) temp
      real(wp) wa1(n)
      real(wp) wa2(n)
      real(wp) wa3(n)
      real(wp) wa4(n)
      real(wp) x1(n)
      real(wp) x2(n)
      real(wp) x3(n)
      real(wp) xnorm
      real(wp) xtol

      epsmch = epsilon ( epsmch )

      info = 0
      nfev = 0
    !
    !  Check the input for errors.
    !
      if ( n <= 0 ) then
        return
      else if ( xtol < 0.0D+00 ) then
        return
      else if ( maxfev <= 0 ) then
        return
      else if ( ml < 0 ) then
        return
      else if ( mu < 0 ) then
        return
      else if ( factor <= 0.0D+00 ) then
        return
      else if ( ldfjac < n ) then
        return
      else if ( lr < ( n * ( n + 1 ) ) / 2 ) then
        return
      end if

      if ( mode == 2 ) then

        do j = 1, n
          if ( diag(j) <= 0.0D+00 ) then
            return
          end if
        end do

      end if
    !
    !  Evaluate the function at the starting point
    !  and calculate its norm.
    !
      call bdf2_residual ( dydt, n, t1, x1, t2, x2, t3, x3, fvec )
      nfev = 1

      fnorm = enorm ( n, fvec )
    !
    !  Determine the number of calls needed to compute the jacobian matrix.
    !
      msum = min ( ml + mu + 1, n )
    !
    !  Initialize iteration counter and monitors.
    !
      iter = 1
      ncsuc = 0
      ncfail = 0
      nslow1 = 0
      nslow2 = 0
    !
    !  Beginning of the outer loop.
    !
      do

        jeval = .true.
    !
    !  Calculate the jacobian matrix.
    !
        call fdjac_bdf2 ( dydt, n, t1, x1, t2, x2, t3, x3, fvec, fjac, ldfjac, &
          ml, mu, epsfcn )

        nfev = nfev + msum
    !
    !  Compute the QR factorization of the jacobian.
    !
        pivot = .false.
        call qrfac ( n, n, fjac, ldfjac, pivot, iwa, 1, wa1, wa2 )
    !
    !  On the first iteration, if MODE is 1, scale according
    !  to the norms of the columns of the initial jacobian.
    !
        if ( iter == 1 ) then

          if ( mode /= 2 ) then

            diag(1:n) = wa2(1:n)
            do j = 1, n
              if ( wa2(j) == 0.0D+00 ) then
                diag(j) = 1.0D+00
              end if
            end do

          end if
    !
    !  On the first iteration, calculate the norm of the scaled X
    !  and initialize the step bound DELTA.
    !
          wa3(1:n) = diag(1:n) * x3(1:n)
          xnorm = enorm ( n, wa3 )
          delta = factor * xnorm
          if ( delta == 0.0D+00 ) then
            delta = factor
          end if

        end if
    !
    !  Form Q' * FVEC and store in QTF.
    !
        qtf(1:n) = fvec(1:n)

        do j = 1, n

          if ( fjac(j,j) /= 0.0D+00 ) then
            temp = - dot_product ( qtf(j:n), fjac(j:n,j) ) / fjac(j,j)
            qtf(j:n) = qtf(j:n) + fjac(j:n,j) * temp
          end if

        end do
    !
    !  Copy the triangular factor of the QR factorization into R.
    !
        sing = .false.

        do j = 1, n
          l = j
          do i = 1, j - 1
            r(l) = fjac(i,j)
            l = l + n - i
          end do
          r(l) = wa1(j)
          if ( wa1(j) == 0.0D+00 ) then
            sing = .true.
          end if
        end do
    !
    !  Accumulate the orthogonal factor in FJAC.
    !
        call qform ( n, n, fjac, ldfjac )
    !
    !  Rescale if necessary.
    !
        if ( mode /= 2 ) then
          do j = 1, n
            diag(j) = max ( diag(j), wa2(j) )
          end do
        end if
    !
    !  Beginning of the inner loop.
    !
        do
    !
    !  Determine the direction P.
    !
          call dogleg ( n, r, lr, diag, qtf, delta, wa1 )
    !
    !  Store the direction P and X + P.
    !  Calculate the norm of P.
    !
          wa1(1:n) = - wa1(1:n)
          wa2(1:n) = x3(1:n) + wa1(1:n)
          wa3(1:n) = diag(1:n) * wa1(1:n)

          pnorm = enorm ( n, wa3 )
    !
    !  On the first iteration, adjust the initial step bound.
    !
          if ( iter == 1 ) then
            delta = min ( delta, pnorm )
          end if
    !
    !  Evaluate the function at X + P and calculate its norm.
    !
          call bdf2_residual ( dydt, n, t1, x1, t2, x2, t3, wa2, wa4 )
          nfev = nfev + 1
          fnorm1 = enorm ( n, wa4 )
    !
    !  Compute the scaled actual reduction.
    !
          actred = -1.0D+00
          if ( fnorm1 < fnorm ) then
            actred = 1.0D+00 - ( fnorm1 / fnorm ) ** 2
          end if
    !
    !  Compute the scaled predicted reduction.
    !
          l = 1
          do i = 1, n
            sum2 = 0.0D+00
            do j = i, n
              sum2 = sum2 + r(l) * wa1(j)
              l = l + 1
            end do
            wa3(i) = qtf(i) + sum2
          end do

          temp = enorm ( n, wa3 )
          prered = 0.0D+00
          if ( temp < fnorm ) then
            prered = 1.0D+00 - ( temp / fnorm ) ** 2
          end if
    !
    !  Compute the ratio of the actual to the predicted reduction.
    !
          ratio = 0.0D+00
          if ( 0.0D+00 < prered ) then
            ratio = actred / prered
          end if
    !
    !  Update the step bound.
    !
          if ( ratio < 0.1D+00 ) then

            ncsuc = 0
            ncfail = ncfail + 1
            delta = 0.5D+00 * delta

          else

            ncfail = 0
            ncsuc = ncsuc + 1

            if ( 0.5D+00 <= ratio .or. 1 < ncsuc ) then
              delta = max ( delta, pnorm / 0.5D+00 )
            end if

            if ( abs ( ratio - 1.0D+00 ) <= 0.1D+00 ) then
              delta = pnorm / 0.5D+00
            end if

          end if
    !
    !  Successful iteration.
    !  Update X, FVEC, and their norms.
    !
          if ( 0.0001D+00 <= ratio ) then
            x3(1:n) = wa2(1:n)
            wa2(1:n) = diag(1:n) * x3(1:n)
            fvec(1:n) = wa4(1:n)
            xnorm = enorm ( n, wa2 )
            fnorm = fnorm1
            iter = iter + 1
          end if
    !
    !  Determine the progress of the iteration.
    !
          nslow1 = nslow1 + 1
          if ( 0.001D+00 <= actred ) then
            nslow1 = 0
          end if

          if ( jeval ) then
            nslow2 = nslow2 + 1
          end if

          if ( 0.1D+00 <= actred ) then
            nslow2 = 0
          end if
    !
    !  Test for convergence.
    !
          if ( delta <= xtol * xnorm .or. fnorm == 0.0D+00 ) then
            info = 1
          end if

          if ( info /= 0 ) then
            return
          end if
    !
    !  Tests for termination and stringent tolerances.
    !
          if ( maxfev <= nfev ) then
            info = 2
          end if

          if ( 0.1D+00 * max ( 0.1D+00 * delta, pnorm ) <= epsmch * xnorm ) then
            info = 3
          end if

          if ( nslow2 == 5 ) then
            info = 4
          end if

          if ( nslow1 == 10 ) then
            info = 5
          end if

          if ( info /= 0 ) then
            return
          end if
    !
    !  Criterion for recalculating jacobian approximation
    !  by forward differences.
    !
          if ( ncfail == 2 ) then
            exit
          end if
    !
    !  Calculate the rank one modification to the jacobian
    !  and update QTF if necessary.
    !
          do j = 1, n
            sum2 = dot_product ( wa4(1:n), fjac(1:n,j) )
            wa2(j) = ( sum2 - wa3(j) ) / pnorm
            wa1(j) = diag(j) * ( ( diag(j) * wa1(j) ) / pnorm )
            if ( 0.0001D+00 <= ratio ) then
              qtf(j) = sum2
            end if
          end do
    !
    !  Compute the QR factorization of the updated jacobian.
    !
          call r1updt ( n, n, r, lr, wa1, wa2, wa3, sing )
          call r1mpyq ( n, n, fjac, ldfjac, wa2, wa3 )
          call r1mpyq ( 1, n, qtf, 1, wa2, wa3 )
    !
    !  End of the inner loop.
    !
          jeval = .false.

        end do
    !
    !  End of the outer loop.
    !
      end do

      return
    end
    subroutine hybrd_bdf3 ( dydt, n, t1, x1, t2, x2, t3, x3, t4, x4, fvec, xtol, &
      maxfev, ml, mu, epsfcn, diag, mode, factor, info, nfev, fjac, ldfjac, r, &
      lr, qtf )

    !*****************************************************************************80
    !
    !! hybrd_bdf3() seeks a zero of N nonlinear equations in N variables.
    !
    !  Discussion:
    !
    !    The code finds a zero of a system of N nonlinear functions in N variables
    !    by a modification of the Powell hybrid method.  The user must provide a
    !    subroutine which calculates the functions.  
    !
    !    The jacobian is then calculated by a forward-difference approximation.
    !
    !  Licensing:
    !
    !    This code is distributed under the MIT license.
    !
    !  Modified:
    !
    !    27 May 2025
    !
    !  Author:
    !
    !    Original Fortran77 version by Jorge More, Burton Garbow, Kenneth Hillstrom.
    !    This version by John Burkardt.
    !
    !  Reference:
    !
    !    Jorge More, Burton Garbow, Kenneth Hillstrom,
    !    User Guide for MINPACK-1,
    !    Technical Report ANL-80-74,
    !    Argonne National Laboratory, 1980.
    !
    !  Input:
    !
    !    external dydt(), the name of the user-supplied code which
    !    evaluates the right hand side of the ODE, of the form:
    !      subroutine dydt ( t, y, dy )
    !      real(wp) dy(n)
    !      real(wp) t
    !      real(wp) y(n)
    !
    !    integer n: the number of functions and variables.
    !
    !    real(wp) t1, x1(n), t2, x2(n), t3, x3(n), t4, x4(n): 
    !    a sequence of four times and solution estimates.
    !
    !    real(wp) xtol: Termination occurs when the relative error
    !    between two consecutive iterates is at most XTOL.  XTOL should be
    !    nonnegative.
    !
    !    integer maxfev: Termination occurs when the number of
    !    calls to the derivative code is at least MAXFEV by the end of an iteration.
    !
    !    integer ml, mu: specify the number of subdiagonals and
    !    superdiagonals within the band of the jacobian matrix.  If the jacobian
    !    is not banded, set ML and MU to at least n - 1.
    !
    !    real(wp) epsfcn: is used in determining a suitable step
    !    length for the forward-difference approximation.  This approximation
    !    assumes that the relative errors in the functions are of the order of
    !    EPSFCN.  If EPSFCN is less than the machine precision, it is assumed that
    !    the relative errors in the functions are of the order of the machine
    !    precision.
    !
    !    real(wp) diag(n): multiplicative scale factors for the 
    !    variables.  Only needed as input if MODE = 2.  
    !
    !    integer mode: scaling option.
    !    1, variables will be scaled internally.
    !    2, scaling is specified by the input DIAG vector.
    !
    !    real(wp) factor: determines the initial step bound.  This
    !    bound is set to the product of FACTOR and the euclidean norm of DIAG*X if
    !    nonzero, or else to FACTOR itself.  In most cases, FACTOR should lie
    !    in the interval (0.1, 100) with 100 the recommended value.
    !
    !    integer ldfjac: the leading dimension of FJAC.
    !    LDFJAC must be at least N.
    !
    !    integer lr: the size of the R array, which must be no
    !    less than (N*(N+1))/2.
    !
    !  Output:
    !
    !    real(wp) x4(n): the final estimate of the solution vector.
    !
    !    real(wp) fvec(n): the functions evaluated at the output x.
    !
    !    real(wp) diag(n): multiplicative scale factors for the 
    !    variables.  Set internally if MODE = 1.
    !    
    !    integer info: status flag.  A value of 1 indicates success. 
    !    0, improper input parameter values.
    !    1, relative error between two consecutive iterates is at most XTOL.
    !    2, number of calls to derivative has reached or exceeded MAXFEV.
    !    3, XTOL is too small.  No further improvement in the approximate
    !       solution X is possible.
    !    4, iteration is not making good progress, as measured by the improvement
    !       from the last five jacobian evaluations.
    !    5, iteration is not making good progress, as measured by the improvement
    !       from the last ten iterations.
    !
    !    integer nfev: the number of calls to the derivative function.
    !
    !    real(wp) fjac(ldfjac,n): an N by N array which contains
    !    the orthogonal matrix Q produced by the QR factorization of the final
    !    approximate jacobian.
    !
    !    real(wp) r(lr): the upper triangular matrix produced by
    !    the QR factorization of the final approximate jacobian, stored rowwise.
    !
    !    real(wp) qtf(n): contains the vector Q'*FVEC.
    !
       
      integer ldfjac
      integer lr
      integer n

      real(wp) actred
      real(wp) delta
      real(wp) diag(n)
      external dydt 
      real(wp) epsfcn
      real(wp) epsmch
      real(wp) factor
      real(wp) fjac(ldfjac,n)
      real(wp) fnorm
      real(wp) fnorm1
      real(wp) fvec(n)
      integer i
      integer info
      integer iter
      integer iwa(1)
      integer j
      logical jeval
      integer l
      integer maxfev
      integer ml
      integer mode
      integer msum
      integer mu
      integer ncfail
      integer ncsuc
      integer nfev
      integer nslow1
      integer nslow2
      logical pivot
      real(wp) pnorm
      real(wp) prered
      real(wp) qtf(n)
      real(wp) r(lr)
      real(wp) ratio
      logical sing
      real(wp) sum2
      real(wp) t1
      real(wp) t2
      real(wp) t3
      real(wp) t4
      real(wp) temp
      real(wp) wa1(n)
      real(wp) wa2(n)
      real(wp) wa3(n)
      real(wp) wa4(n)
      real(wp) x1(n)
      real(wp) x2(n)
      real(wp) x3(n)
      real(wp) x4(n)
      real(wp) xnorm
      real(wp) xtol

      epsmch = epsilon ( epsmch )

      info = 0
      nfev = 0
    !
    !  Check the input for errors.
    !
      if ( n <= 0 ) then
        return
      else if ( xtol < 0.0D+00 ) then
        return
      else if ( maxfev <= 0 ) then
        return
      else if ( ml < 0 ) then
        return
      else if ( mu < 0 ) then
        return
      else if ( factor <= 0.0D+00 ) then
        return
      else if ( ldfjac < n ) then
        return
      else if ( lr < ( n * ( n + 1 ) ) / 2 ) then
        return
      end if

      if ( mode == 2 ) then

        do j = 1, n
          if ( diag(j) <= 0.0D+00 ) then
            return
          end if
        end do

      end if
    !
    !  Evaluate the function at the starting point
    !  and calculate its norm.
    !
      call bdf3_residual ( dydt, n, t1, x1, t2, x2, t3, x3, t4, x4, fvec )
      nfev = 1

      fnorm = enorm ( n, fvec )
    !
    !  Determine the number of calls needed to compute the jacobian matrix.
    !
      msum = min ( ml + mu + 1, n )
    !
    !  Initialize iteration counter and monitors.
    !
      iter = 1
      ncsuc = 0
      ncfail = 0
      nslow1 = 0
      nslow2 = 0
    !
    !  Beginning of the outer loop.
    !
      do

        jeval = .true.
    !
    !  Calculate the jacobian matrix.
    !
        call fdjac_bdf3 ( dydt, n, t1, x1, t2, x2, t3, x3, t4, x4, fvec, fjac, &
          ldfjac, ml, mu, epsfcn )

        nfev = nfev + msum
    !
    !  Compute the QR factorization of the jacobian.
    !
        pivot = .false.
        call qrfac ( n, n, fjac, ldfjac, pivot, iwa, 1, wa1, wa2 )
    !
    !  On the first iteration, if MODE is 1, scale according
    !  to the norms of the columns of the initial jacobian.
    !
        if ( iter == 1 ) then

          if ( mode /= 2 ) then

            diag(1:n) = wa2(1:n)
            do j = 1, n
              if ( wa2(j) == 0.0D+00 ) then
                diag(j) = 1.0D+00
              end if
            end do

          end if
    !
    !  On the first iteration, calculate the norm of the scaled X
    !  and initialize the step bound DELTA.
    !
          wa3(1:n) = diag(1:n) * x4(1:n)
          xnorm = enorm ( n, wa3 )
          delta = factor * xnorm
          if ( delta == 0.0D+00 ) then
            delta = factor
          end if

        end if
    !
    !  Form Q' * FVEC and store in QTF.
    !
        qtf(1:n) = fvec(1:n)

        do j = 1, n

          if ( fjac(j,j) /= 0.0D+00 ) then
            temp = - dot_product ( qtf(j:n), fjac(j:n,j) ) / fjac(j,j)
            qtf(j:n) = qtf(j:n) + fjac(j:n,j) * temp
          end if

        end do
    !
    !  Copy the triangular factor of the QR factorization into R.
    !
        sing = .false.

        do j = 1, n
          l = j
          do i = 1, j - 1
            r(l) = fjac(i,j)
            l = l + n - i
          end do
          r(l) = wa1(j)
          if ( wa1(j) == 0.0D+00 ) then
            sing = .true.
          end if
        end do
    !
    !  Accumulate the orthogonal factor in FJAC.
    !
        call qform ( n, n, fjac, ldfjac )
    !
    !  Rescale if necessary.
    !
        if ( mode /= 2 ) then
          do j = 1, n
            diag(j) = max ( diag(j), wa2(j) )
          end do
        end if
    !
    !  Beginning of the inner loop.
    !
        do
    !
    !  Determine the direction P.
    !
          call dogleg ( n, r, lr, diag, qtf, delta, wa1 )
    !
    !  Store the direction P and X + P.
    !  Calculate the norm of P.
    !
          wa1(1:n) = - wa1(1:n)
          wa2(1:n) = x4(1:n) + wa1(1:n)
          wa3(1:n) = diag(1:n) * wa1(1:n)

          pnorm = enorm ( n, wa3 )
    !
    !  On the first iteration, adjust the initial step bound.
    !
          if ( iter == 1 ) then
            delta = min ( delta, pnorm )
          end if
    !
    !  Evaluate the function at X + P and calculate its norm.
    !
          call bdf3_residual ( dydt, n, t1, x1, t2, x2, t3, x3, t4, wa2, wa4 )
          nfev = nfev + 1
          fnorm1 = enorm ( n, wa4 )
    !
    !  Compute the scaled actual reduction.
    !
          actred = -1.0D+00
          if ( fnorm1 < fnorm ) then
            actred = 1.0D+00 - ( fnorm1 / fnorm ) ** 2
          end if
    !
    !  Compute the scaled predicted reduction.
    !
          l = 1
          do i = 1, n
            sum2 = 0.0D+00
            do j = i, n
              sum2 = sum2 + r(l) * wa1(j)
              l = l + 1
            end do
            wa3(i) = qtf(i) + sum2
          end do

          temp = enorm ( n, wa3 )
          prered = 0.0D+00
          if ( temp < fnorm ) then
            prered = 1.0D+00 - ( temp / fnorm ) ** 2
          end if
    !
    !  Compute the ratio of the actual to the predicted reduction.
    !
          ratio = 0.0D+00
          if ( 0.0D+00 < prered ) then
            ratio = actred / prered
          end if
    !
    !  Update the step bound.
    !
          if ( ratio < 0.1D+00 ) then

            ncsuc = 0
            ncfail = ncfail + 1
            delta = 0.5D+00 * delta

          else

            ncfail = 0
            ncsuc = ncsuc + 1

            if ( 0.5D+00 <= ratio .or. 1 < ncsuc ) then
              delta = max ( delta, pnorm / 0.5D+00 )
            end if

            if ( abs ( ratio - 1.0D+00 ) <= 0.1D+00 ) then
              delta = pnorm / 0.5D+00
            end if

          end if
    !
    !  Successful iteration.
    !  Update X, FVEC, and their norms.
    !
          if ( 0.0001D+00 <= ratio ) then
            x4(1:n) = wa2(1:n)
            wa2(1:n) = diag(1:n) * x4(1:n)
            fvec(1:n) = wa4(1:n)
            xnorm = enorm ( n, wa2 )
            fnorm = fnorm1
            iter = iter + 1
          end if
    !
    !  Determine the progress of the iteration.
    !
          nslow1 = nslow1 + 1
          if ( 0.001D+00 <= actred ) then
            nslow1 = 0
          end if

          if ( jeval ) then
            nslow2 = nslow2 + 1
          end if

          if ( 0.1D+00 <= actred ) then
            nslow2 = 0
          end if
    !
    !  Test for convergence.
    !
          if ( delta <= xtol * xnorm .or. fnorm == 0.0D+00 ) then
            info = 1
          end if

          if ( info /= 0 ) then
            return
          end if
    !
    !  Tests for termination and stringent tolerances.
    !
          if ( maxfev <= nfev ) then
            info = 2
          end if

          if ( 0.1D+00 * max ( 0.1D+00 * delta, pnorm ) <= epsmch * xnorm ) then
            info = 3
          end if

          if ( nslow2 == 5 ) then
            info = 4
          end if

          if ( nslow1 == 10 ) then
            info = 5
          end if

          if ( info /= 0 ) then
            return
          end if
    !
    !  Criterion for recalculating jacobian approximation
    !  by forward differences.
    !
          if ( ncfail == 2 ) then
            exit
          end if
    !
    !  Calculate the rank one modification to the jacobian
    !  and update QTF if necessary.
    !
          do j = 1, n
            sum2 = dot_product ( wa4(1:n), fjac(1:n,j) )
            wa2(j) = ( sum2 - wa3(j) ) / pnorm
            wa1(j) = diag(j) * ( ( diag(j) * wa1(j) ) / pnorm )
            if ( 0.0001D+00 <= ratio ) then
              qtf(j) = sum2
            end if
          end do
    !
    !  Compute the QR factorization of the updated jacobian.
    !
          call r1updt ( n, n, r, lr, wa1, wa2, wa3, sing )
          call r1mpyq ( n, n, fjac, ldfjac, wa2, wa3 )
          call r1mpyq ( 1, n, qtf, 1, wa2, wa3 )
    !
    !  End of the inner loop.
    !
          jeval = .false.

        end do
    !
    !  End of the outer loop.
    !
      end do

      return
    end
    subroutine hybrd_be ( dydt, n, to, xo, t, x, fvec, xtol, maxfev, ml, mu, &
      epsfcn, diag, mode, factor, info, nfev, fjac, ldfjac, r, lr, qtf )

    !*****************************************************************************80
    !
    !! hybrd_be() seeks a zero of N nonlinear equations in N variables.
    !
    !  Discussion:
    !
    !    The code finds a zero of a system of N nonlinear functions in N variables
    !    by a modification of the Powell hybrid method.  The user must provide a
    !    subroutine which calculates the functions.  
    !
    !    The jacobian is then calculated by a forward-difference approximation.
    !
    !    The original code hybrd() was modified to deal with problems
    !    involving a backward Euler residual.
    !
    !  Licensing:
    !
    !    This code is distributed under the MIT license.
    !
    !  Modified:
    !
    !    08 November 2023
    !
    !  Author:
    !
    !    Original Fortran77 version by Jorge More, Burton Garbow, Kenneth Hillstrom.
    !    This version by John Burkardt.
    !
    !  Reference:
    !
    !    Jorge More, Burton Garbow, Kenneth Hillstrom,
    !    User Guide for MINPACK-1,
    !    Technical Report ANL-80-74,
    !    Argonne National Laboratory, 1980.
    !
    !  Input:
    !
    !    external dydt(), the name of the user-supplied code which
    !    evaluates the right hand side of the ODE, of the form:
    !      subroutine dydt ( t, y, dy )
    !      real(wp) dy(n)
    !      real(wp) t
    !      real(wp) y(n)
    !
    !    integer n: the number of functions and variables.
    !
    !    real(wp) to, xo(n): the old time and solution.
    !   
    !    real(wp) t, x(n): the new time and current solution estimate. 
    !
    !    real(wp) xtol: Termination occurs when the relative error
    !    between two consecutive iterates is at most XTOL.  XTOL should be
    !    nonnegative.
    !
    !    integer maxfev: Termination occurs when the number of
    !    calls to the derivative code is at least MAXFEV by the end of an iteration.
    !
    !    integer ml, mu: specify the number of subdiagonals and
    !    superdiagonals within the band of the jacobian matrix.  If the jacobian
    !    is not banded, set ML and MU to at least n - 1.
    !
    !    real(wp) epsfcn: is used in determining a suitable step
    !    length for the forward-difference approximation.  This approximation
    !    assumes that the relative errors in the functions are of the order of
    !    EPSFCN.  If EPSFCN is less than the machine precision, it is assumed that
    !    the relative errors in the functions are of the order of the machine
    !    precision.
    !
    !    real(wp) diag(n): multiplicative scale factors for the 
    !    variables.  Only needed as input if MODE = 2.  
    !
    !    integer mode: scaling option.
    !    1, variables will be scaled internally.
    !    2, scaling is specified by the input DIAG vector.
    !
    !    real(wp) factor: determines the initial step bound.  This
    !    bound is set to the product of FACTOR and the euclidean norm of DIAG*X if
    !    nonzero, or else to FACTOR itself.  In most cases, FACTOR should lie
    !    in the interval (0.1, 100) with 100 the recommended value.
    !
    !    integer ldfjac: the leading dimension of FJAC.
    !    LDFJAC must be at least N.
    !
    !    integer lr: the size of the R array, which must be no
    !    less than (N*(N+1))/2.
    !
    !  Output:
    !
    !    real(wp) x(n): the final estimate of the solution vector.
    !
    !    real(wp) fvec(n): the functions evaluated at the output x.
    !
    !    real(wp) diag(n): multiplicative scale factors for the 
    !    variables.  Set internally if MODE = 1.
    !    
    !    integer info: status flag.  A value of 1 indicates success. 
    !    0, improper input parameter values.
    !    1, relative error between two consecutive iterates is at most XTOL.
    !    2, number of calls to derivative has reached or exceeded MAXFEV.
    !    3, XTOL is too small.  No further improvement in the approximate
    !       solution X is possible.
    !    4, iteration is not making good progress, as measured by the improvement
    !       from the last five jacobian evaluations.
    !    5, iteration is not making good progress, as measured by the improvement
    !       from the last ten iterations.
    !
    !    integer nfev: the number of calls to the derivative function.
    !
    !    real(wp) fjac(ldfjac,n): an N by N array which contains
    !    the orthogonal matrix Q produced by the QR factorization of the final
    !    approximate jacobian.
    !
    !    real(wp) r(lr): the upper triangular matrix produced by
    !    the QR factorization of the final approximate jacobian, stored rowwise.
    !
    !    real(wp) qtf(n): contains the vector Q'*FVEC.
    !
       
      integer ldfjac
      integer lr
      integer n

      real(wp) actred
      real(wp) delta
      real(wp) diag(n)
      external dydt 
      real(wp) epsfcn
      real(wp) epsmch
      real(wp) factor
      real(wp) fjac(ldfjac,n)
      real(wp) fnorm
      real(wp) fnorm1
      real(wp) fvec(n)
      integer i
      integer info
      integer iter
      integer iwa(1)
      integer j
      logical jeval
      integer l
      integer maxfev
      integer ml
      integer mode
      integer msum
      integer mu
      integer ncfail
      integer ncsuc
      integer nfev
      integer nslow1
      integer nslow2
      logical pivot
      real(wp) pnorm
      real(wp) prered
      real(wp) qtf(n)
      real(wp) r(lr)
      real(wp) ratio
      logical sing
      real(wp) sum2
      real(wp) t
      real(wp) temp
      real(wp) to
      real(wp) wa1(n)
      real(wp) wa2(n)
      real(wp) wa3(n)
      real(wp) wa4(n)
      real(wp) x(n)
      real(wp) xnorm
      real(wp) xo(n)
      real(wp) xtol

      epsmch = epsilon ( epsmch )

      info = 0
      nfev = 0
    !
    !  Check the input for errors.
    !
      if ( n <= 0 ) then
        return
      else if ( xtol < 0.0D+00 ) then
        return
      else if ( maxfev <= 0 ) then
        return
      else if ( ml < 0 ) then
        return
      else if ( mu < 0 ) then
        return
      else if ( factor <= 0.0D+00 ) then
        return
      else if ( ldfjac < n ) then
        return
      else if ( lr < ( n * ( n + 1 ) ) / 2 ) then
        return
      end if

      if ( mode == 2 ) then

        do j = 1, n
          if ( diag(j) <= 0.0D+00 ) then
            return
          end if
        end do

      end if
    !
    !  Evaluate the function at the starting point
    !  and calculate its norm.
    !
      call backward_euler_residual ( dydt, n, to, xo, t, x, fvec )
      nfev = 1

      fnorm = enorm ( n, fvec )
    !
    !  Determine the number of calls needed to compute the jacobian matrix.
    !
      msum = min ( ml + mu + 1, n )
    !
    !  Initialize iteration counter and monitors.
    !
      iter = 1
      ncsuc = 0
      ncfail = 0
      nslow1 = 0
      nslow2 = 0
    !
    !  Beginning of the outer loop.
    !
      do

        jeval = .true.
    !
    !  Calculate the jacobian matrix.
    !
        call fdjac_be ( dydt, n, to, xo, t, x, fvec, fjac, ldfjac, ml, mu, epsfcn )

        nfev = nfev + msum
    !
    !  Compute the QR factorization of the jacobian.
    !
        pivot = .false.
        call qrfac ( n, n, fjac, ldfjac, pivot, iwa, 1, wa1, wa2 )
    !
    !  On the first iteration, if MODE is 1, scale according
    !  to the norms of the columns of the initial jacobian.
    !
        if ( iter == 1 ) then

          if ( mode /= 2 ) then

            diag(1:n) = wa2(1:n)
            do j = 1, n
              if ( wa2(j) == 0.0D+00 ) then
                diag(j) = 1.0D+00
              end if
            end do

          end if
    !
    !  On the first iteration, calculate the norm of the scaled X
    !  and initialize the step bound DELTA.
    !
          wa3(1:n) = diag(1:n) * x(1:n)
          xnorm = enorm ( n, wa3 )
          delta = factor * xnorm
          if ( delta == 0.0D+00 ) then
            delta = factor
          end if

        end if
    !
    !  Form Q' * FVEC and store in QTF.
    !
        qtf(1:n) = fvec(1:n)

        do j = 1, n

          if ( fjac(j,j) /= 0.0D+00 ) then
            temp = - dot_product ( qtf(j:n), fjac(j:n,j) ) / fjac(j,j)
            qtf(j:n) = qtf(j:n) + fjac(j:n,j) * temp
          end if

        end do
    !
    !  Copy the triangular factor of the QR factorization into R.
    !
        sing = .false.

        do j = 1, n
          l = j
          do i = 1, j - 1
            r(l) = fjac(i,j)
            l = l + n - i
          end do
          r(l) = wa1(j)
          if ( wa1(j) == 0.0D+00 ) then
            sing = .true.
          end if
        end do
    !
    !  Accumulate the orthogonal factor in FJAC.
    !
        call qform ( n, n, fjac, ldfjac )
    !
    !  Rescale if necessary.
    !
        if ( mode /= 2 ) then
          do j = 1, n
            diag(j) = max ( diag(j), wa2(j) )
          end do
        end if
    !
    !  Beginning of the inner loop.
    !
        do
    !
    !  Determine the direction P.
    !
          call dogleg ( n, r, lr, diag, qtf, delta, wa1 )
    !
    !  Store the direction P and X + P.
    !  Calculate the norm of P.
    !
          wa1(1:n) = - wa1(1:n)
          wa2(1:n) = x(1:n) + wa1(1:n)
          wa3(1:n) = diag(1:n) * wa1(1:n)

          pnorm = enorm ( n, wa3 )
    !
    !  On the first iteration, adjust the initial step bound.
    !
          if ( iter == 1 ) then
            delta = min ( delta, pnorm )
          end if
    !
    !  Evaluate the function at X + P and calculate its norm.
    !
          call backward_euler_residual ( dydt, n, to, xo, t, wa2, wa4 )
          nfev = nfev + 1
          fnorm1 = enorm ( n, wa4 )
    !
    !  Compute the scaled actual reduction.
    !
          actred = -1.0D+00
          if ( fnorm1 < fnorm ) then
            actred = 1.0D+00 - ( fnorm1 / fnorm ) ** 2
          end if
    !
    !  Compute the scaled predicted reduction.
    !
          l = 1
          do i = 1, n
            sum2 = 0.0D+00
            do j = i, n
              sum2 = sum2 + r(l) * wa1(j)
              l = l + 1
            end do
            wa3(i) = qtf(i) + sum2
          end do

          temp = enorm ( n, wa3 )
          prered = 0.0D+00
          if ( temp < fnorm ) then
            prered = 1.0D+00 - ( temp / fnorm ) ** 2
          end if
    !
    !  Compute the ratio of the actual to the predicted reduction.
    !
          ratio = 0.0D+00
          if ( 0.0D+00 < prered ) then
            ratio = actred / prered
          end if
    !
    !  Update the step bound.
    !
          if ( ratio < 0.1D+00 ) then

            ncsuc = 0
            ncfail = ncfail + 1
            delta = 0.5D+00 * delta

          else

            ncfail = 0
            ncsuc = ncsuc + 1

            if ( 0.5D+00 <= ratio .or. 1 < ncsuc ) then
              delta = max ( delta, pnorm / 0.5D+00 )
            end if

            if ( abs ( ratio - 1.0D+00 ) <= 0.1D+00 ) then
              delta = pnorm / 0.5D+00
            end if

          end if
    !
    !  Successful iteration.
    !  Update X, FVEC, and their norms.
    !
          if ( 0.0001D+00 <= ratio ) then
            x(1:n) = wa2(1:n)
            wa2(1:n) = diag(1:n) * x(1:n)
            fvec(1:n) = wa4(1:n)
            xnorm = enorm ( n, wa2 )
            fnorm = fnorm1
            iter = iter + 1
          end if
    !
    !  Determine the progress of the iteration.
    !
          nslow1 = nslow1 + 1
          if ( 0.001D+00 <= actred ) then
            nslow1 = 0
          end if

          if ( jeval ) then
            nslow2 = nslow2 + 1
          end if

          if ( 0.1D+00 <= actred ) then
            nslow2 = 0
          end if
    !
    !  Test for convergence.
    !
          if ( delta <= xtol * xnorm .or. fnorm == 0.0D+00 ) then
            info = 1
          end if

          if ( info /= 0 ) then
            return
          end if
    !
    !  Tests for termination and stringent tolerances.
    !
          if ( maxfev <= nfev ) then
            info = 2
          end if

          if ( 0.1D+00 * max ( 0.1D+00 * delta, pnorm ) <= epsmch * xnorm ) then
            info = 3
          end if

          if ( nslow2 == 5 ) then
            info = 4
          end if

          if ( nslow1 == 10 ) then
            info = 5
          end if

          if ( info /= 0 ) then
            return
          end if
    !
    !  Criterion for recalculating jacobian approximation
    !  by forward differences.
    !
          if ( ncfail == 2 ) then
            exit
          end if
    !
    !  Calculate the rank one modification to the jacobian
    !  and update QTF if necessary.
    !
          do j = 1, n
            sum2 = dot_product ( wa4(1:n), fjac(1:n,j) )
            wa2(j) = ( sum2 - wa3(j) ) / pnorm
            wa1(j) = diag(j) * ( ( diag(j) * wa1(j) ) / pnorm )
            if ( 0.0001D+00 <= ratio ) then
              qtf(j) = sum2
            end if
          end do
    !
    !  Compute the QR factorization of the updated jacobian.
    !
          call r1updt ( n, n, r, lr, wa1, wa2, wa3, sing )
          call r1mpyq ( n, n, fjac, ldfjac, wa2, wa3 )
          call r1mpyq ( 1, n, qtf, 1, wa2, wa3 )
    !
    !  End of the inner loop.
    !
          jeval = .false.

        end do
    !
    !  End of the outer loop.
    !
      end do

      return
    end
    subroutine hybrd_tr ( dydt, n, to, xo, tn, xn, fvec, xtol, maxfev, ml, mu, &
      epsfcn, diag, mode, factor, info, nfev, fjac, ldfjac, r, lr, qtf )

    !*****************************************************************************80
    !
    !! hybrd_tr() seeks a zero of N nonlinear equations in N variables.
    !
    !  Discussion:
    !
    !    The code finds a zero of a system of N nonlinear functions in N variables
    !    by a modification of the Powell hybrid method.  The user must provide a
    !    subroutine which calculates the functions.  
    !
    !    The jacobian is then calculated by a forward-difference approximation.
    !
    !    The original code hybrd() was modified to deal with problems
    !    involving an implicit trapezoidal ODE residual.
    !
    !  Licensing:
    !
    !    This code is distributed under the MIT license.
    !
    !  Modified:
    !
    !    15 November 2023
    !
    !  Author:
    !
    !    Original Fortran77 version by Jorge More, Burton Garbow, Kenneth Hillstrom.
    !    This version by John Burkardt.
    !
    !  Reference:
    !
    !    Jorge More, Burton Garbow, Kenneth Hillstrom,
    !    User Guide for MINPACK-1,
    !    Technical Report ANL-80-74,
    !    Argonne National Laboratory, 1980.
    !
    !  Input:
    !
    !    external dydt(), the name of the user-supplied code which
    !    evaluates the right hand side of the ODE, of the form:
    !      subroutine dydt ( t, y, dy )
    !      real(wp) dy(n)
    !      real(wp) t
    !      real(wp) y(n)
    !
    !    integer n: the number of functions and variables.
    !
    !    real(wp) to, xo(n): the old time and solution.
    !   
    !    real(wp) tn, xn(n): the new time and current solution estimate. 
    !
    !    real(wp) xtol: Termination occurs when the relative error
    !    between two consecutive iterates is at most XTOL.  XTOL should be
    !    nonnegative.
    !
    !    integer maxfev: Termination occurs when the number of
    !    calls to the derivative code is at least MAXFEV by the end of an iteration.
    !
    !    integer ml, mu: specify the number of subdiagonals and
    !    superdiagonals within the band of the jacobian matrix.  If the jacobian
    !    is not banded, set ML and MU to at least n - 1.
    !
    !    real(wp) epsfcn: is used in determining a suitable step
    !    length for the forward-difference approximation.  This approximation
    !    assumes that the relative errors in the functions are of the order of
    !    EPSFCN.  If EPSFCN is less than the machine precision, it is assumed that
    !    the relative errors in the functions are of the order of the machine
    !    precision.
    !
    !    real(wp) diag(n): multiplicative scale factors for the 
    !    variables.  Only needed as input if MODE = 2.  
    !
    !    integer mode: scaling option.
    !    1, variables will be scaled internally.
    !    2, scaling is specified by the input DIAG vector.
    !
    !    real(wp) factor: determines the initial step bound.  This
    !    bound is set to the product of FACTOR and the euclidean norm of DIAG*X if
    !    nonzero, or else to FACTOR itself.  In most cases, FACTOR should lie
    !    in the interval (0.1, 100) with 100 the recommended value.
    !
    !    integer ldfjac: the leading dimension of FJAC.
    !    LDFJAC must be at least N.
    !
    !    integer lr: the size of the R array, which must be no
    !    less than (N*(N+1))/2.
    !
    !  Output:
    !
    !    real(wp) xn(n): the final estimate of the solution vector.
    !
    !    real(wp) fvec(n): the functions evaluated at the output xn.
    !
    !    real(wp) diag(n): multiplicative scale factors for the 
    !    variables.  Set internally if MODE = 1.
    !    
    !    integer info: status flag.  A value of 1 indicates success. 
    !    0, improper input parameter values.
    !    1, relative error between two consecutive iterates is at most XTOL.
    !    2, number of calls to derivative has reached or exceeded MAXFEV.
    !    3, XTOL is too small.  No further improvement in the approximate
    !       solution X is possible.
    !    4, iteration is not making good progress, as measured by the improvement
    !       from the last five jacobian evaluations.
    !    5, iteration is not making good progress, as measured by the improvement
    !       from the last ten iterations.
    !
    !    integer nfev: the number of calls to the derivative function.
    !
    !    real(wp) fjac(ldfjac,n): an N by N array which contains
    !    the orthogonal matrix Q produced by the QR factorization of the final
    !    approximate jacobian.
    !
    !    real(wp) r(lr): the upper triangular matrix produced by
    !    the QR factorization of the final approximate jacobian, stored rowwise.
    !
    !    real(wp) qtf(n): contains the vector Q'*FVEC.
    !
           
       
      integer ldfjac
      integer lr
      integer n
    
      real(wp) actred
      real(wp) delta
      real(wp) diag(n)
      external dydt 
      real(wp) epsfcn
      real(wp) epsmch
      real(wp) factor
      real(wp) fjac(ldfjac,n)
      real(wp) fnorm
      real(wp) fnorm1
      real(wp) fvec(n)
      integer i
      integer info
      integer iter
      integer iwa(1)
      integer j
      logical jeval
      integer l
      integer maxfev
      integer ml
      integer mode
      integer msum
      integer mu
      integer ncfail
      integer ncsuc
      integer nfev
      integer nslow1
      integer nslow2
      logical pivot
      real(wp) pnorm
      real(wp) prered
      real(wp) qtf(n)
      real(wp) r(lr)
      real(wp) ratio
      logical sing
      real(wp) sum2
      real(wp) temp
      real(wp) tn
      real(wp) to
      real(wp) wa1(n)
      real(wp) wa2(n)
      real(wp) wa3(n)
      real(wp) wa4(n)
      real(wp) xn(n)
      real(wp) xnorm
      real(wp) xo(n)
      real(wp) xtol
    
      epsmch = epsilon ( epsmch )
    
      info = 0
      nfev = 0
    !
    !  Check the input for errors.
    !
      if ( n <= 0 ) then
        return
      else if ( xtol < 0.0D+00 ) then
        return
      else if ( maxfev <= 0 ) then
        return
      else if ( ml < 0 ) then
        return
      else if ( mu < 0 ) then
        return
      else if ( factor <= 0.0D+00 ) then
        return
      else if ( ldfjac < n ) then
        return
      else if ( lr < ( n * ( n + 1 ) ) / 2 ) then
        return
      end if
    
      if ( mode == 2 ) then
    
        do j = 1, n
          if ( diag(j) <= 0.0D+00 ) then
            return
          end if
        end do
    
      end if
    !
    !  Evaluate the function at the starting point
    !  and calculate its norm.
    !
      call trapezoidal_residual ( dydt, n, to, xo, tn, xn, fvec )
      nfev = 1
    
      fnorm = enorm ( n, fvec )
    !
    !  Determine the number of calls needed to compute the jacobian matrix.
    !
      msum = min ( ml + mu + 1, n )
    !
    !  Initialize iteration counter and monitors.
    !
      iter = 1
      ncsuc = 0
      ncfail = 0
      nslow1 = 0
      nslow2 = 0
    !
    !  Beginning of the outer loop.
    !
      do
    
        jeval = .true.
    !
    !  Calculate the jacobian matrix.
    !
        call fdjac_tr ( dydt, n, to, xo, tn, xn, fvec, fjac, ldfjac, &
          ml, mu, epsfcn )
    
        nfev = nfev + msum
    !
    !  Compute the QR factorization of the jacobian.
    !
        pivot = .false.
        call qrfac ( n, n, fjac, ldfjac, pivot, iwa, 1, wa1, wa2 )
    !
    !  On the first iteration, if MODE is 1, scale according
    !  to the norms of the columns of the initial jacobian.
    !
        if ( iter == 1 ) then
    
          if ( mode /= 2 ) then
    
            diag(1:n) = wa2(1:n)
            do j = 1, n
              if ( wa2(j) == 0.0D+00 ) then
                diag(j) = 1.0D+00
              end if
            end do
    
          end if
    !
    !  On the first iteration, calculate the norm of the scaled X
    !  and initialize the step bound DELTA.
    !
          wa3(1:n) = diag(1:n) * xn(1:n)
          xnorm = enorm ( n, wa3 )
          delta = factor * xnorm
          if ( delta == 0.0D+00 ) then
            delta = factor
          end if
    
        end if
    !
    !  Form Q' * FVEC and store in QTF.
    !
        qtf(1:n) = fvec(1:n)
    
        do j = 1, n
    
          if ( fjac(j,j) /= 0.0D+00 ) then
            temp = - dot_product ( qtf(j:n), fjac(j:n,j) ) / fjac(j,j)
            qtf(j:n) = qtf(j:n) + fjac(j:n,j) * temp
          end if
    
        end do
    !
    !  Copy the triangular factor of the QR factorization into R.
    !
        sing = .false.
    
        do j = 1, n
          l = j
          do i = 1, j - 1
            r(l) = fjac(i,j)
            l = l + n - i
          end do
          r(l) = wa1(j)
          if ( wa1(j) == 0.0D+00 ) then
            sing = .true.
          end if
        end do
    !
    !  Accumulate the orthogonal factor in FJAC.
    !
        call qform ( n, n, fjac, ldfjac )
    !
    !  Rescale if necessary.
    !
        if ( mode /= 2 ) then
          do j = 1, n
            diag(j) = max ( diag(j), wa2(j) )
          end do
        end if
    !
    !  Beginning of the inner loop.
    !
        do
    !
    !  Determine the direction P.
    !
          call dogleg ( n, r, lr, diag, qtf, delta, wa1 )
    !
    !  Store the direction P and X + P.
    !  Calculate the norm of P.
    !
          wa1(1:n) = - wa1(1:n)
          wa2(1:n) = xn(1:n) + wa1(1:n)
          wa3(1:n) = diag(1:n) * wa1(1:n)
    
          pnorm = enorm ( n, wa3 )
    !
    !  On the first iteration, adjust the initial step bound.
    !
          if ( iter == 1 ) then
            delta = min ( delta, pnorm )
          end if
    !
    !  Evaluate the function at X + P and calculate its norm.
    !
          call trapezoidal_residual ( dydt, n, to, xo, tn, wa2, wa4 )
          nfev = nfev + 1
          fnorm1 = enorm ( n, wa4 )
    !
    !  Compute the scaled actual reduction.
    !
          actred = -1.0D+00
          if ( fnorm1 < fnorm ) then
            actred = 1.0D+00 - ( fnorm1 / fnorm ) ** 2
          end if
    !
    !  Compute the scaled predicted reduction.
    !
          l = 1
          do i = 1, n
            sum2 = 0.0D+00
            do j = i, n
              sum2 = sum2 + r(l) * wa1(j)
              l = l + 1
            end do
            wa3(i) = qtf(i) + sum2
          end do
    
          temp = enorm ( n, wa3 )
          prered = 0.0D+00
          if ( temp < fnorm ) then
            prered = 1.0D+00 - ( temp / fnorm ) ** 2
          end if
    !
    !  Compute the ratio of the actual to the predicted reduction.
    !
          ratio = 0.0D+00
          if ( 0.0D+00 < prered ) then
            ratio = actred / prered
          end if
    !
    !  Update the step bound.
    !
          if ( ratio < 0.1D+00 ) then
    
            ncsuc = 0
            ncfail = ncfail + 1
            delta = 0.5D+00 * delta
    
          else
    
            ncfail = 0
            ncsuc = ncsuc + 1
    
            if ( 0.5D+00 <= ratio .or. 1 < ncsuc ) then
              delta = max ( delta, pnorm / 0.5D+00 )
            end if
    
            if ( abs ( ratio - 1.0D+00 ) <= 0.1D+00 ) then
              delta = pnorm / 0.5D+00
            end if
    
          end if
    !
    !  Successful iteration.
    !  Update X, FVEC, and their norms.
    !
          if ( 0.0001D+00 <= ratio ) then
            xn(1:n) = wa2(1:n)
            wa2(1:n) = diag(1:n) * xn(1:n)
            fvec(1:n) = wa4(1:n)
            xnorm = enorm ( n, wa2 )
            fnorm = fnorm1
            iter = iter + 1
          end if
    !
    !  Determine the progress of the iteration.
    !
          nslow1 = nslow1 + 1
          if ( 0.001D+00 <= actred ) then
            nslow1 = 0
          end if
    
          if ( jeval ) then
            nslow2 = nslow2 + 1
          end if
    
          if ( 0.1D+00 <= actred ) then
            nslow2 = 0
          end if
    !
    !  Test for convergence.
    !
          if ( delta <= xtol * xnorm .or. fnorm == 0.0D+00 ) then
            info = 1
          end if
    
          if ( info /= 0 ) then
            return
          end if
    !
    !  Tests for termination and stringent tolerances.
    !
          if ( maxfev <= nfev ) then
            info = 2
          end if
    
          if ( 0.1D+00 * max ( 0.1D+00 * delta, pnorm ) <= epsmch * xnorm ) then
            info = 3
          end if
    
          if ( nslow2 == 5 ) then
            info = 4
          end if
    
          if ( nslow1 == 10 ) then
            info = 5
          end if
    
          if ( info /= 0 ) then
            return
          end if
    !
    !  Criterion for recalculating jacobian approximation
    !  by forward differences.
    !
          if ( ncfail == 2 ) then
            exit
          end if
    !
    !  Calculate the rank one modification to the jacobian
    !  and update QTF if necessary.
    !
          do j = 1, n
            sum2 = dot_product ( wa4(1:n), fjac(1:n,j) )
            wa2(j) = ( sum2 - wa3(j) ) / pnorm
            wa1(j) = diag(j) * ( ( diag(j) * wa1(j) ) / pnorm )
            if ( 0.0001D+00 <= ratio ) then
              qtf(j) = sum2
            end if
          end do
    !
    !  Compute the QR factorization of the updated jacobian.
    !
          call r1updt ( n, n, r, lr, wa1, wa2, wa3, sing )
          call r1mpyq ( n, n, fjac, ldfjac, wa2, wa3 )
          call r1mpyq ( 1, n, qtf, 1, wa2, wa3 )
    !
    !  End of the inner loop.
    !
          jeval = .false.
    
        end do
    !
    !  End of the outer loop.
    !
      end do
    
      return
    end
    subroutine qform ( m, n, q, ldq )
    
    !*****************************************************************************80
    !
    !! qform() produces the explicit QR factorization of a matrix.
    !
    !  Discussion:
    !
    !    The QR factorization of a matrix is usually accumulated in implicit
    !    form, that is, as a series of orthogonal transformations of the
    !    original matrix.  This routine carries out those transformations,
    !    to explicitly exhibit the factorization constructed by QRFAC.
    !
    !  Licensing:
    !
    !    This code is distributed under the MIT license.
    !
    !  Modified:
    !
    !    06 April 2010
    !
    !  Author:
    !
    !    Original Fortran77 version by Jorge More, Burton Garbow, Kenneth Hillstrom.
    !    This version by John Burkardt.
    !
    !  Reference:
    !
    !    Jorge More, Burton Garbow, Kenneth Hillstrom,
    !    User Guide for MINPACK-1,
    !    Technical Report ANL-80-74,
    !    Argonne National Laboratory, 1980.
    !
    !  Input:
    !
    !    integer M, is a positive integer input variable set
    !    to the number of rows of A and the order of Q.
    !
    !    integer N, is a positive integer input variable set
    !    to the number of columns of A.
    !
    !    real(wp) Q(LDQ,M).  Q is an M by M array.
    !    the full lower trapezoid in the first min(M,N) columns of Q
    !    contains the factored form.
    !
    !    integer LDQ, is a positive integer input variable 
    !    not less than M which specifies the leading dimension of the array Q.
    !
    !  Output:
    !
    !    real(wp) Q(LDQ,M): an M by M array.
    !    Q has been accumulated into a square matrix.
    !
           
       
      integer ldq
      integer m
      integer n
    
      integer j
      integer k
      integer l
      integer minmn
      real(wp) q(ldq,m)
      real(wp) temp
      real(wp) wa(m)
    
      minmn = min ( m, n )
    
      do j = 2, minmn
        q(1:j-1,j) = 0.0D+00
      end do
    !
    !  Initialize remaining columns to those of the identity matrix.
    !
      q(1:m,n+1:m) = 0.0D+00
    
      do j = n + 1, m
        q(j,j) = 1.0D+00
      end do
    !
    !  Accumulate Q from its factored form.
    !
      do l = 1, minmn
    
        k = minmn - l + 1
    
        wa(k:m) = q(k:m,k)
    
        q(k:m,k) = 0.0D+00
        q(k,k) = 1.0D+00
    
        if ( wa(k) /= 0.0D+00 ) then
    
          do j = k, m
            temp = dot_product ( wa(k:m), q(k:m,j) ) / wa(k)
            q(k:m,j) = q(k:m,j) - temp * wa(k:m)
          end do
    
        end if
    
      end do
    
      return
    end
    subroutine qrfac ( m, n, a, lda, pivot, ipvt, lipvt, rdiag, acnorm )
    
    !*****************************************************************************80
    !
    !! qrfac() computes a QR factorization using Householder transformations.
    !
    !  Discussion:
    !
    !    This function uses Householder transformations with optional column
    !    pivoting to compute a QR factorization of the
    !    M by N matrix A.  That is, QRFAC determines an orthogonal
    !    matrix Q, a permutation matrix P, and an upper trapezoidal
    !    matrix R with diagonal elements of nonincreasing magnitude,
    !    such that A*P = Q*R.  
    !
    !    The Householder transformation for column K, K = 1,2,...,min(M,N), 
    !    is of the form
    !
    !      I - ( 1 / U(K) ) * U * U'
    !
    !    where U has zeros in the first K-1 positions.  
    !
    !    The form of this transformation and the method of pivoting first
    !    appeared in the corresponding LINPACK routine.
    !
    !  Licensing:
    !
    !    This code is distributed under the MIT license.
    !
    !  Modified:
    !
    !    06 April 2010
    !
    !  Author:
    !
    !    Original Fortran77 version by Jorge More, Burton Garbow, Kenneth Hillstrom.
    !    This version by John Burkardt.
    !
    !  Reference:
    !
    !    Jorge More, Burton Garbow, Kenneth Hillstrom,
    !    User Guide for MINPACK-1,
    !    Technical Report ANL-80-74,
    !    Argonne National Laboratory, 1980.
    !
    !  Input:
    !
    !    integer M, the number of rows of A.
    !
    !    integer N, the number of columns of A.
    !
    !    real(wp) A(LDA,N): the M by N matrix for which the QR 
    !    factorization is to be computed.  
    !
    !    integer LDA, the leading dimension of A, which must
    !    be no less than M.
    !
    !    logical PIVOT, is TRUE if column pivoting is to be carried out.
    !
    !    integer LIPVT, the dimension of IPVT, which should 
    !    be N if pivoting is used.
    !
    !  Output:
    !
    !    real(wp) A(LDA,N), the M by N array.
    !    The strict upper trapezoidal part of A contains
    !    the strict upper trapezoidal part of R, and the lower trapezoidal
    !    part of A contains a factored form of Q, the non-trivial elements of
    !    the U vectors described above.
    !
    !    integer IPVT(LIPVT), defines the permutation matrix P 
    !    such that A*P = Q*R.  Column J of P is column IPVT(J) of the identity 
    !    matrix.  If PIVOT is false, IPVT is not referenced.
    !
    !    real(wp) RDIAG(N), contains the diagonal elements of R.
    !
    !    real(wp) ACNORM(N), the norms of the corresponding
    !    columns of the input matrix A.  If this information is not needed,
    !    then ACNORM can coincide with RDIAG.
    !
           
       
      integer lda
      integer lipvt
      integer m
      integer n
    
      real(wp) a(lda,n)
      real(wp) acnorm(n)
      real(wp) ajnorm 
      real(wp) epsmch
      integer i4_temp
      integer ipvt(lipvt)
      integer j
      integer k
      integer kmax
      integer minmn
      logical pivot
      real(wp) r8_temp(m)
      real(wp) rdiag(n)
      real(wp) temp
      real(wp) wa(n)
    
      epsmch = epsilon ( epsmch )
    !
    !  Compute the initial column norms and initialize several arrays.
    !
      do j = 1, n
        acnorm(j) = enorm ( m, a(1:m,j) )
      end do
    
      rdiag(1:n) = acnorm(1:n)
      wa(1:n) = acnorm(1:n)
    
      if ( pivot ) then
        do j = 1, n
          ipvt(j) = j
        end do
      end if
    !
    !  Reduce A to R with Householder transformations.
    !
      minmn = min ( m, n )
    
      do j = 1, minmn
    !
    !  Bring the column of largest norm into the pivot position.
    !
        if ( pivot ) then
    
          kmax = j
    
          do k = j, n
            if ( rdiag(kmax) < rdiag(k) ) then
              kmax = k
            end if
          end do
    
          if ( kmax /= j ) then
    
            r8_temp(1:m) = a(1:m,j)
            a(1:m,j)     = a(1:m,kmax)
            a(1:m,kmax)  = r8_temp(1:m)
    
            rdiag(kmax) = rdiag(j)
            wa(kmax) = wa(j)
    
            i4_temp    = ipvt(j)
            ipvt(j)    = ipvt(kmax)
            ipvt(kmax) = i4_temp
    
          end if
    
        end if
    !
    !  Compute the Householder transformation to reduce the
    !  J-th column of A to a multiple of the J-th unit vector.
    !
        ajnorm = enorm ( m-j+1, a(j,j) )
    
        if ( ajnorm /= 0.0D+00 ) then
    
          if ( a(j,j) < 0.0D+00 ) then
            ajnorm = -ajnorm
          end if
    
          a(j:m,j) = a(j:m,j) / ajnorm
          a(j,j) = a(j,j) + 1.0D+00
    !
    !  Apply the transformation to the remaining columns and update the norms.
    !
          do k = j + 1, n
    
            temp = dot_product ( a(j:m,j), a(j:m,k) ) / a(j,j)
    
            a(j:m,k) = a(j:m,k) - temp * a(j:m,j)
    
            if ( pivot .and. rdiag(k) /= 0.0D+00 ) then
    
              temp = a(j,k) / rdiag(k)
              rdiag(k) = rdiag(k) * sqrt ( max ( 0.0D+00, 1.0D+00 - temp ** 2 ) )
    
              if ( 0.05D+00 * ( rdiag(k) / wa(k) ) ** 2 <= epsmch ) then
                rdiag(k) = enorm ( m-j, a(j+1,k) )
                wa(k) = rdiag(k)
              end if
    
            end if
    
          end do
    
        end if
    
        rdiag(j) = - ajnorm
    
      end do
    
      return
    end
    subroutine r1mpyq ( m, n, a, lda, v, w )
    
    !*****************************************************************************80
    !
    !! r1mpyq() computes A*Q, where Q is the product of Householder transformations.
    !
    !  Discussion:
    !
    !    Given an M by N matrix A, this function computes A*Q where
    !    Q is the product of 2*(N - 1) transformations
    !
    !      GV(N-1)*...*GV(1)*GW(1)*...*GW(N-1)
    !
    !    and GV(I), GW(I) are Givens rotations in the (I,N) plane which
    !    eliminate elements in the I-th and N-th planes, respectively.
    !    Q itself is not given, rather the information to recover the
    !    GV, GW rotations is supplied.
    !
    !  Licensing:
    !
    !    This code is distributed under the MIT license.
    !
    !  Modified:
    !
    !    06 April 2010
    !
    !  Author:
    !
    !    Original Fortran77 version by Jorge More, Burton Garbow, Kenneth Hillstrom.
    !    This version by John Burkardt.
    !
    !  Reference:
    !
    !    Jorge More, Burton Garbow, Kenneth Hillstrom,
    !    User Guide for MINPACK-1,
    !    Technical Report ANL-80-74,
    !    Argonne National Laboratory, 1980.
    !
    !  Input:
    !
    !    integer M, the number of rows of A.
    !
    !    integer N, the number of columns of A.
    !
    !    real(wp) A(LDA,N), the M by N matrix A to be 
    !    postmultiplied by the orthogonal matrix Q.
    !
    !    integer LDA, the leading dimension of A, which must not
    !    be less than M.
    !
    !    real(wp) V(N), W(N), contain the information necessary
    !    to recover the Givens rotations GV and GW.
    !
    !  Output:
    !
    !    real(wp) A(LDA,N): the value of A*Q.
    !
           
       
      integer lda
      integer m
      integer n
    
      real(wp) a(lda,n)
      real(wp) c
      integer i
      integer j
      real(wp) s
      real(wp) temp
      real(wp) v(n)
      real(wp) w(n)
    !
    !  Apply the first set of Givens rotations to A.
    !
      do j = n - 1, 1, -1
    
        if ( 1.0D+00 < abs ( v(j) ) ) then
          c = 1.0D+00 / v(j)
          s = sqrt ( 1.0D+00 - c ** 2 )
        else
          s = v(j)
          c = sqrt ( 1.0D+00 - s ** 2 )
        end if
    
        do i = 1, m
          temp =   c * a(i,j) - s * a(i,n)
          a(i,n) = s * a(i,j) + c * a(i,n)
          a(i,j) = temp
        end do
    
      end do
    !
    !  Apply the second set of Givens rotations to A.
    !
      do j = 1, n - 1
    
        if ( 1.0D+00 < abs ( w(j) ) ) then
          c = 1.0D+00 / w(j)
          s = sqrt ( 1.0D+00 - c ** 2 )
        else
          s = w(j)
          c = sqrt ( 1.0D+00 - s ** 2 )
        end if
    
        do i = 1, m
          temp =     c * a(i,j) + s * a(i,n)
          a(i,n) = - s * a(i,j) + c * a(i,n)
          a(i,j) = temp
        end do
    
      end do
    
      return
    end
    subroutine r1updt ( m, n, s, ls, u, v, w, sing )
    
    !*****************************************************************************80
    !
    !! r1updt() re-triangularizes a matrix after a rank one update.
    !
    !  Discussion:
    !
    !    Given an M by N lower trapezoidal matrix S, an M-vector U, and an
    !    N-vector V, the problem is to determine an orthogonal matrix Q such that
    !
    !      (S + U * V' ) * Q
    !
    !    is again lower trapezoidal.
    !
    !    This function determines Q as the product of 2 * (N - 1)
    !    transformations
    !
    !      GV(N-1)*...*GV(1)*GW(1)*...*GW(N-1)
    !
    !    where GV(I), GW(I) are Givens rotations in the (I,N) plane
    !    which eliminate elements in the I-th and N-th planes,
    !    respectively.  Q itself is not accumulated, rather the
    !    information to recover the GV and GW rotations is returned.
    !
    !  Licensing:
    !
    !    This code is distributed under the MIT license.
    !
    !  Modified:
    !
    !    06 April 2010
    !
    !  Author:
    !
    !    Original Fortran77 version by Jorge More, Burton Garbow, Kenneth Hillstrom.
    !    This version by John Burkardt.
    !
    !  Reference:
    !
    !    Jorge More, Burton Garbow, Kenneth Hillstrom,
    !    User Guide for MINPACK-1,
    !    Technical Report ANL-80-74,
    !    Argonne National Laboratory, 1980.
    !
    !  Input:
    !
    !    integer M, the number of rows of S.
    !
    !    integer N, the number of columns of S.  
    !    N must not exceed M.
    !
    !    real(wp) S(LS): the lower trapezoidal matrix S stored 
    !    by columns.  
    !
    !    integer LS, the length of the S array.  LS must be at
    !    least (N*(2*M-N+1))/2.
    !
    !    real(wp) U(M), the U vector.
    !
    !    real(wp) V(N): the vector V.
    !
    !  Output:
    !
    !    real(wp) S(LS): the lower trapezoidal
    !    matrix produced as described above.
    !
    !    real(wp) V(N): the information necessary to recover the
    !    Givens rotations GV described above.
    !
    !    real(wp) W(M), contains information necessary to
    !    recover the Givens rotations GW described above.
    !
    !    logical SING, is set to TRUE if any of the diagonal elements
    !    of the output S are zero.  Otherwise SING is set FALSE.
    !
           
       
      integer ls
      integer m
      integer n
    
      real(wp) cos
      real(wp) cotan
      real(wp) giant
      integer i
      integer j
      integer jj
      integer l
      real(wp) s(ls)
      real(wp) sin
      logical sing
      real(wp) tan
      real(wp) tau
      real(wp) temp
      real(wp) u(m)
      real(wp) v(n)
      real(wp) w(m)
    !
    !  GIANT is the largest magnitude.
    !
      giant = huge ( giant )
    !
    !  Initialize the diagonal element pointer.
    !
      jj = ( n * ( 2 * m - n + 1 ) ) / 2 - ( m - n )
    !
    !  Move the nontrivial part of the last column of S into W.
    !
      l = jj
      do i = n, m
        w(i) = s(l)
        l = l + 1
      end do
    !
    !  Rotate the vector V into a multiple of the N-th unit vector
    !  in such a way that a spike is introduced into W.
    !
      do j = n - 1, 1, -1
    
        jj = jj - ( m - j + 1 )
        w(j) = 0.0D+00
    
        if ( v(j) /= 0.0D+00 ) then
    !
    !  Determine a Givens rotation which eliminates the J-th element of V.
    !
          if ( abs ( v(n) ) < abs ( v(j) ) ) then
            cotan = v(n) / v(j)
            sin = 0.5D+00 / sqrt ( 0.25D+00 + 0.25D+00 * cotan ** 2 )
            cos = sin * cotan
            tau = 1.0D+00
            if ( abs ( cos ) * giant > 1.0D+00 ) then
              tau = 1.0D+00 / cos
            end if
          else
            tan = v(j) / v(n)
            cos = 0.5D+00 / sqrt ( 0.25D+00 + 0.25D+00 * tan ** 2 )
            sin = cos * tan
            tau = sin
          end if
    !
    !  Apply the transformation to V and store the information
    !  necessary to recover the Givens rotation.
    !
          v(n) = sin * v(j) + cos * v(n)
          v(j) = tau
    !
    !  Apply the transformation to S and extend the spike in W.
    !
          l = jj
          do i = j, m
            temp = cos * s(l) - sin * w(i)
            w(i) = sin * s(l) + cos * w(i)
            s(l) = temp
            l = l + 1
          end do
    
        end if
    
      end do
    !
    !  Add the spike from the rank 1 update to W.
    !
       w(1:m) = w(1:m) + v(n) * u(1:m)
    !
    !  Eliminate the spike.
    !
      sing = .false.
    
      do j = 1, n-1
    
        if ( w(j) /= 0.0D+00 ) then
    !
    !  Determine a Givens rotation which eliminates the
    !  J-th element of the spike.
    !
          if ( abs ( s(jj) ) < abs ( w(j) ) ) then
    
            cotan = s(jj) / w(j)
            sin = 0.5D+00 / sqrt ( 0.25D+00 + 0.25D+00 * cotan ** 2 )
            cos = sin * cotan
    
            if ( 1.0D+00 < abs ( cos ) * giant ) then
              tau = 1.0D+00 / cos
            else
              tau = 1.0D+00
            end if
    
          else
    
            tan = w(j) / s(jj)
            cos = 0.5D+00 / sqrt ( 0.25D+00 + 0.25D+00 * tan ** 2 )
            sin = cos * tan
            tau = sin
    
          end if
    !
    !  Apply the transformation to S and reduce the spike in W.
    !
          l = jj
          do i = j, m
            temp = cos * s(l) + sin * w(i)
            w(i) = - sin * s(l) + cos * w(i)
            s(l) = temp
            l = l + 1
          end do
    !
    !  Store the information necessary to recover the Givens rotation.
    !
          w(j) = tau
    
        end if
    !
    !  Test for zero diagonal elements in the output S.
    !
        if ( s(jj) == 0.0D+00 ) then
          sing = .true.
        end if
    
        jj = jj + ( m - j + 1 )
    
      end do
    !
    !  Move W back into the last column of the output S.
    !
      l = jj
      do i = n, m
        s(l) = w(i)
        l = l + 1
      end do
    
      if ( s(jj) == 0.0D+00 ) then
        sing = .true.
      end if
    
      return
    end
    subroutine trapezoidal_residual ( dydt, n, to, yo, tn, yn, ft )
    
    !*****************************************************************************80
    !
    !! trapezoidal_residual() evaluates the trapezoidal ODE solver residual.
    !
    !  Discussion:
    !
    !    Let to and tn be two times, with yo and yn the associated ODE
    !    solution values there.  If yn satisfies the implicit trapezoidal
    !    ODE solver condiion, then:
    !
    !      0.5 * ( dydt(to,yo) + dydt(tn,yn) ) = ( yn - yo ) / ( tn - to )
    !
    !    This can be rewritten as
    !
    !      residual = yn - yo - 0.5 * ( tn - to ) * ( dydt(to,yo) + dydt(tn,yn) )
    !
    !    Given the other information as fixed values, a nonlinear equation 
    !    solver can be used to estimate the value yn that makes the residual zero.
    !
    !  Licensing:
    !
    !    This code is distributed under the MIT license.
    !
    !  Modified:
    !
    !    15 November 2023
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Input:
    !
    !    external dydt(): the name of the user-supplied code which
    !    evaluates the right hand side of the ODE, of the form:
    !      subroutine dydt ( t, y, dy )
    !      real(wp) dy(*)
    !      real(wp) t
    !      real(wp) y(*)
    !
    !    integer n: the vector size.
    !
    !    real(wp) to, yo(n): the old time and solution.
    !
    !    real(wp) tn, yn(n): the new time and tentative solution.
    !
    !  Output:
    !
    !    real(wp) ft(n): the trapezoidal residual.
    !
           
       
      integer n
    
      real(wp) dydtn(n)
      real(wp) dydto(n)
      external dydt
      real(wp) ft(n)
      real(wp) tn
      real(wp) to
      real(wp) yn(n)
      real(wp) yo(n)
    
      call dydt ( to, yo, dydto )
    
      call dydt ( tn, yn, dydtn )
    
      ft = yn - yo - 0.5 * ( tn - to ) * ( dydto + dydtn )
    
      return
    end

end module pastr_fsolver    