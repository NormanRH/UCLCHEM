PROGRAM nan_mask
use, intrinsic :: ieee_arithmetic, only: IEEE_Value, IEEE_QUIET_NAN
use, intrinsic :: iso_fortran_env, only: real32
real(real32) :: nan
double precision ::  values(4)=(/ 2.0,1.0e+01,3.0e+00,1.0e+01/)
double precision :: a,b


    nan = IEEE_VALUE(nan, IEEE_QUIET_NAN)
    values(1)=nan
    a=SUM(values,MASK=.NOT.ISNAN(values))
    b=SUM(values)
    write(*,*) a,b,nan,HUGE(1.0)

END PROGRAM nan_mask