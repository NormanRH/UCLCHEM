MODULE CONSTANTS
   use, intrinsic :: iso_fortran_env, dp=>real64

   REAL(dp), parameter :: C  = 2.99792458D+10 !Speed of light in cgs
   REAL(dp), PARAMETER :: K_BOLTZ = 1.38065040D-16 ! Boltzmann constant cgs
   REAL(dp), PARAMETER :: HP = 6.62606896D-27 !Planck constant in cgs
   REAL(dp), PARAMETER :: REDUCED_PLANCK=1.054571628d-27
   REAL(dp), PARAMETER :: MH = 1.67262164D-24 !H nucleus mass in cgs
   REAL(dp), PARAMETER :: PI = 3.141592654
   REAL(dp), PARAMETER :: K_BOLTZ_SI=1.38d-23 !Boltzmann constant SI
   REAL(dp), PARAMETER :: PC=3.086d18 !parsec in cgs
   REAL(dp), PARAMETER :: SECONDS_PER_YEAR=3.16d7
   REAL(dp), PARAMETER :: T_CMB=2.73
   REAL(dp), PARAMETER :: EV = 1.60217646D-12 ! electron volt in erg
   REAL(dp), PARAMETER :: GRAV_G = 6.674d-8 !gravitational constant in cgs
   REAL(dp), PARAMETER :: SB_CONST=5.6704d-5 !Stefan Boltzmann constant in cgs
END MODULE CONSTANTS