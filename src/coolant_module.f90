MODULE COOLANT_MODULE
   USE constants
   USE network
   use, intrinsic :: iso_fortran_env, dp=>real64
   IMPLICIT NONE

!  Specify the properties that define each coolant species
   TYPE COOLANT_TYPE

      CHARACTER(LEN=10) :: NAME ! Coolant species name
      CHARACTER(LEN=256):: FILENAME ! Name of the coolant data file

      INTEGER :: INDEX  ! Index number of the coolant species
      INTEGER :: NLEVEL ! Number of levels in the system
      INTEGER :: NTEMP  ! Number of temperatures for which collisional rates are available

      REAL(dp) :: MOLECULAR_MASS, previousCooling ! Molecular mass of the coolant species

      REAL(dp), ALLOCATABLE :: ENERGY(:),WEIGHT(:) ! Energy (K) and statistical weight of each level
      REAL(dp), ALLOCATABLE :: A_COEFF(:,:),B_COEFF(:,:) ! Einstein A and B coefficients for each transition between levels
      REAL(dp), ALLOCATABLE :: FREQUENCY(:,:) ! Frequency (Hz) of each transition between levels
      REAL(dp), ALLOCATABLE :: TEMPERATURE(:,:) ! Temperatures (K) at which collisional rates are available for each collision partner
      REAL(dp), ALLOCATABLE :: C_COEFF(:,:,:,:) ! Collisional rate coefficient (cm^3 s^-1) for each transition at each specified temperature for each collision partner

      LOGICAL :: CONVERGED ! Flag indicating whether the level populations of all particles have converged

      !Below were stolen from population type since I think we just want one coolant.
      REAL(dp) :: DENSITY ! Total number density of the species (cm^-3)
      REAL(dp) :: LINEWIDTH ! Doppler line width of the species (cm s^-1)


      REAL(dp), ALLOCATABLE :: POPULATION(:) ! Population density (cm^-3) of each level
      REAL(dp), ALLOCATABLE :: PREVIOUS_POPULATION(:) ! Population density calculated at the previous iteration step
      REAL(dp), ALLOCATABLE :: EMISSIVITY(:,:) ! Local emissivity (erg cm^-3 s^-1) of each transition
      REAL(dp), ALLOCATABLE :: OPACITY(:,:) ! Optical depth of each transition along each HEALPix ray to the PDR surface (or simulation boundary)
      REAL(dp), ALLOCATABLE :: LAMBDA(:,:) ! Lambda operator value for each transition


   END TYPE COOLANT_TYPE



   TYPE(COOLANT_TYPE), allocatable :: coolants(:)
   integer,PARAMETER :: NCOOL=7
   CHARACTER(*), PARAMETER :: coolantFiles(NCOOL)=(/"ly-a.dat       ","12c+_nometa.dat","16o.dat        ","12c.dat        ",&
            &"12co.dat       ","p-h2.dat       ","o-h2.dat       "/)
   CHARACTER(*), PARAMETER :: coolantNames(NCOOL)=(/"H  ","C+ ","O  ","C  ","CO ","H2 ","H2 "/)
   ! CHARACTER(*), PARAMETER :: coolantFiles(NCOOL)=(/"12co.dat       ","p-h2.dat       ","o-h2.dat       ","p-h2o.dat      ","o-h2o.dat      "/)
   ! CHARACTER(*), PARAMETER :: coolantNames(NCOOL)=(/"CO ","H2 ","H2 ","H2O","H2O"/)
   INTEGER :: coolantIndices(NCOOL)
   REAL(dp) :: CLOUD_SIZE


CONTAINS
   !=======================================================================
   !  STOLEN FROM UCL-PDR
   !  Read in the coolant species, their energy level structure, transition
   !  properties and collisional rate coefficients. The specified files are
   !  assumed to contain entries in the LAMDA/RADEX format, allowing files
   !  to be downloaded directly from the online database.
   !
   !-----------------------------------------------------------------------
   SUBROUTINE READ_COOLANTS()
      IMPLICIT NONE
      INTEGER :: I,J,K,L,M,N,INDEX,IER
      INTEGER :: NLEVEL,NLINE,NTEMP,NPARTNER,NCOLL,PARTNER_ID,MAX_NTEMP

      IF (ALLOCATED(coolants)) DEALLOCATE(coolants)
      ALLOCATE(coolants(NCOOL))
      DO N=1,NCOOL ! Loop over coolants
         coolants(N)%FILENAME=coolantFiles(N)
   !     Open the input file
         OPEN(UNIT=1,FILE='Datafiles/Collisional-Rates/'//coolants(N)%FILENAME,IOSTAT=IER,ACTION='READ',STATUS='OLD')
         READ(1,*,IOSTAT=IER) ! Skip the first comment line

   !     Produce an error message if the file does not exist (or cannot be opened for whatever reason)
         IF(IER.NE.0) THEN
            WRITE(6,*) 'ERROR! Cannot open coolant data file ',TRIM(coolants(N)%FILENAME),' for input'
            WRITE(6,*)
            CLOSE(1)
            STOP
         END IF

         READ(1,*,IOSTAT=IER) coolants(N)%NAME ! Read the name of the coolant
         READ(1,*,IOSTAT=IER)
         READ(1,*,IOSTAT=IER) coolants(N)%MOLECULAR_MASS ! Read the molecular mass
         READ(1,*,IOSTAT=IER)
         coolants(N)%INDEX=0 ! Initialize the coolant species index (assigned later)

   !     Read the number of levels and allocate the energy, statistical weight,
   !     Einstein A & B coefficient and transition frequency arrays accordingly
         READ(1,*,IOSTAT=IER) NLEVEL
         READ(1,*,IOSTAT=IER)
         IF(NLEVEL.LT.2) THEN
            WRITE(6,*) 'ERROR! Incorrect number of energy levels in coolant data file ',&
                        &TRIM(coolants(N)%FILENAME),' (NLEVEL=',NLEVEL,')'
            CLOSE(1)
            STOP
         END IF
         coolants(N)%NLEVEL=NLEVEL
         ALLOCATE(coolants(N)%ENERGY(1:NLEVEL))
         ALLOCATE(coolants(N)%WEIGHT(1:NLEVEL))
         ALLOCATE(coolants(N)%A_COEFF(1:NLEVEL,1:NLEVEL))
         ALLOCATE(coolants(N)%B_COEFF(1:NLEVEL,1:NLEVEL))
         ALLOCATE(coolants(N)%FREQUENCY(1:NLEVEL,1:NLEVEL))

   !     Initialize the level energies, statistical weights,
   !     Einstein coefficients and transition frequencies
         coolants(N)%ENERGY=0.0D0
         coolants(N)%WEIGHT=0.0D0
         coolants(N)%A_COEFF=0.0D0
         coolants(N)%B_COEFF=0.0D0
         coolants(N)%FREQUENCY=0.0D0

   !     Read the energy (cm^-1) and statistical weight of each level
         DO L=1,NLEVEL ! Loop over levels
            READ(1,*,IOSTAT=IER) I,coolants(N)%ENERGY(I),coolants(N)%WEIGHT(I)
            coolants(N)%ENERGY(I)=coolants(N)%ENERGY(I)*C*HP ! Convert from cm^-1 to erg
         END DO ! End of loop over levels
         READ(1,*,IOSTAT=IER)

   !     Read the Einstein A coefficient (s^-1) and frequency (GHz) of each radiative transition
         READ(1,*,IOSTAT=IER) NLINE
         READ(1,*,IOSTAT=IER)
         DO L=1,NLINE ! Loop over radiative transitions
            READ(1,*,IOSTAT=IER) INDEX,I,J,coolants(N)%A_COEFF(I,J),coolants(N)%FREQUENCY(I,J)
            coolants(N)%FREQUENCY(I,J)=coolants(N)%FREQUENCY(I,J)*1.0D9 ! Convert from GHz to Hz
   !        Calculate the Einstein B coefficient using B_ij = A_ij/(2.h.nu^3/c^2)
            coolants(N)%B_COEFF(I,J)=coolants(N)%A_COEFF(I,J)/(2*HP*(coolants(N)%FREQUENCY(I,J)**3)/(C**2))
   !        Calculate the Einstein B coefficient for the reverse transition from detailed balance
            coolants(N)%B_COEFF(J,I)=coolants(N)%B_COEFF(I,J)*(coolants(N)%WEIGHT(I)/coolants(N)%WEIGHT(J))
         END DO ! End of loop over radiative transitions
         READ(1,*,IOSTAT=IER)

   !     Calculate the transition frequencies between all levels (even if forbidden)
         DO I=1,NLEVEL
            DO J=1,NLEVEL
   !           Check that the calculated and measured frequencies differ by <0.1%
   !           Produce an error message if the difference between them is greater
               IF(coolants(N)%FREQUENCY(I,J).NE.0.0D0) THEN
                  IF(ABS(coolants(N)%FREQUENCY(I,J)-ABS(coolants(N)%ENERGY(I)-coolants(N)%ENERGY(J))/HP) &
                      & /coolants(N)%FREQUENCY(I,J).GT.1.0D-3) THEN
                     WRITE(6,*) 'ERROR! Calculated frequency differs from measured frequency by >0.1%'
                     WRITE(6,"(1PD12.5,'Hz vs',1PD12.5,'Hz')") ABS(coolants(N)%ENERGY(I)-coolants(N)%ENERGY(J))/HP, &
                                                             & coolants(N)%FREQUENCY(I,J)
                     CLOSE(1)
                     STOP
                  END IF
               ELSE
                  coolants(N)%FREQUENCY(I,J)=ABS(coolants(N)%ENERGY(I)-coolants(N)%ENERGY(J))/HP
               END IF
            END DO
         END DO

   !     Allocate and initialize the collisional rate coefficient arrays
   !     allowing a maximum of 1000 temperature values per collision partner
         MAX_NTEMP=1000
         coolants(N)%NTEMP=MAX_NTEMP
         ALLOCATE(coolants(N)%TEMPERATURE(1:7,1:MAX_NTEMP))
         ALLOCATE(coolants(N)%C_COEFF(1:7,1:NLEVEL,1:NLEVEL,1:MAX_NTEMP))
         coolants(N)%TEMPERATURE=0.0D0
         coolants(N)%C_COEFF=0.0D0

   !     Read the collisional rate coefficients (cm^3 s^-1) for each collision partner
         MAX_NTEMP=0
         READ(1,*,IOSTAT=IER) NPARTNER
         DO L=1,NPARTNER ! Loop over collision partners
            READ(1,*,IOSTAT=IER)
            READ(1,*,IOSTAT=IER) PARTNER_ID
            IF(PARTNER_ID.LT.1 .OR. PARTNER_ID.GT.7) THEN
               WRITE(6,*) 'ERROR! Unrecognized collision partner ID in coolant data file ',&
                           &TRIM(coolants(N)%FILENAME),' (ID=',PARTNER_ID,')'
               CLOSE(1)
               STOP
            END IF
            READ(1,*,IOSTAT=IER)
            READ(1,*,IOSTAT=IER) NCOLL
            READ(1,*,IOSTAT=IER)
            READ(1,*,IOSTAT=IER) NTEMP ; MAX_NTEMP=MAX(MAX_NTEMP,NTEMP)
            IF(NTEMP.GT.1000) THEN
               WRITE(6,*) 'ERROR! Number of temperature values exceeds limit in coolant data file ',&
                           &TRIM(coolants(N)%FILENAME),' (NTEMP=',NTEMP,')'
               CLOSE(1)
               STOP
            END IF
            READ(1,*,IOSTAT=IER)
            READ(1,*,IOSTAT=IER) (coolants(N)%TEMPERATURE(PARTNER_ID,K),K=1,NTEMP)
            READ(1,*,IOSTAT=IER)
            DO M=1,NCOLL ! Loop over collisional transitions
               READ(1,*,IOSTAT=IER) INDEX,I,J,(coolants(N)%C_COEFF(PARTNER_ID,I,J,K),K=1,NTEMP)
   !           Calculate the reverse (excitation) rate coefficients from
   !           detailed balance: C_ji = C_ij*gi/gj*exp(-(Ei-Ej)/kT)
               DO K=1,NTEMP ! Loop over temperatures
                  IF(coolants(N)%C_COEFF(PARTNER_ID,I,J,K).NE.0.0D0 .AND. &
                   & coolants(N)%C_COEFF(PARTNER_ID,J,I,K).EQ.0.0D0) THEN
                     coolants(N)%C_COEFF(PARTNER_ID,J,I,K)=coolants(N)%C_COEFF(PARTNER_ID,I,J,K) &
                                                       & *(coolants(N)%WEIGHT(I)/coolants(N)%WEIGHT(J)) &
                                                       & *EXP(-(coolants(N)%ENERGY(I)-coolants(N)%ENERGY(J)) &
                                                       &      /(K_BOLTZ*coolants(N)%TEMPERATURE(PARTNER_ID,K)))
                  END IF
               END DO ! End of loop over temperatures
            END DO ! End of loop over collisional transitions
         END DO ! End of loop over collision partners

         coolants(N)%NTEMP=MAX_NTEMP
         coolants(N)%previousCooling=0.0d0
         CLOSE(1)

      END DO ! End of loop over coolants
      

       DO N=1,NCOOL
         ALLOCATE(coolants(N)%POPULATION(1:coolants(N)%NLEVEL))
         ALLOCATE(coolants(N)%PREVIOUS_POPULATION(1:coolants(N)%NLEVEL))
         ALLOCATE(coolants(N)%EMISSIVITY(1:coolants(N)%NLEVEL,1:coolants(N)%NLEVEL))
         !ALLOCATE(coolants(N)%OPACITY(0:NRAYS-1,1:coolants(N)%NLEVEL,1:coolants(N)%NLEVEL))
         ALLOCATE(coolants(N)%OPACITY(1:coolants(N)%NLEVEL,1:coolants(N)%NLEVEL))
         ALLOCATE(coolants(N)%LAMBDA(1:coolants(N)%NLEVEL,1:coolants(N)%NLEVEL))
         coolants(N)%POPULATION=0.0
         coolants(N)%PREVIOUS_POPULATION=0.0
         coolants(N)%EMISSIVITY=0.0
         coolants(N)%OPACITY=0.0
         coolants(N)%LAMBDA=0.0
      END DO

      RETURN
   END SUBROUTINE READ_COOLANTS   

   !=======================================================================
   !
   !  Calculate the Doppler line width of each coolant species for all
   !  particles. Contributions from both thermal and turbulent motions
   !  are included.
   !
   !-----------------------------------------------------------------------
   SUBROUTINE UPDATE_COOLANT_LINEWIDTHS(GasTemperature)
      IMPLICIT NONE
      REAL(dp), INTENT(IN) :: GasTemperature
      INTEGER :: N

   !  Calculate the mean thermal velocity of each coolant species at the relevant
   !  gas temperature for each particle. Update the Doppler line widths (cm s^-1)
   !  by adding the thermal and turbulent velocities in quadrature.
     DO N=1,NCOOL ! Loop over coolants
      COOLANTS(N)%LINEWIDTH = SQRT(2*K_BOLTZ*GasTemperature/(coolants(N)%MOLECULAR_MASS*MH)) ! v_thermal = (2kT/m)^1/2
     END DO ! End of loop over coolants

      RETURN
   END SUBROUTINE UPDATE_COOLANT_LINEWIDTHS

   !=======================================================================
   !
   !  Calculate the level populations at LTE for the given species.
   !
   !-----------------------------------------------------------------------
   SUBROUTINE CALCULATE_LTE_POPULATIONS(NLEVEL,ENERGY,WEIGHT,DENSITY, &
                                      & TEMPERATURE,POPULATION)

       IMPLICIT NONE

      INTEGER, INTENT(IN)  :: NLEVEL
      REAL(dp),     INTENT(IN)  :: ENERGY(:),WEIGHT(:)
      REAL(dp),     INTENT(IN)  :: DENSITY,TEMPERATURE
      REAL(dp),     INTENT(OUT) :: POPULATION(:)

      INTEGER :: ILEVEL
      REAL(dp) :: PARTITION_FUNCTION

   !  Initialize the level populations
      POPULATION=0.0D0

      PARTITION_FUNCTION=0.0D0
      DO ILEVEL=1,NLEVEL
         POPULATION(ILEVEL)=WEIGHT(ILEVEL)*EXP(-ENERGY(ILEVEL)/(K_BOLTZ*TEMPERATURE))
         IF (isnan(population(ilevel))) write(*,*) "LTE",temperature
         PARTITION_FUNCTION=PARTITION_FUNCTION+POPULATION(ILEVEL)
         IF (isnan(PARTITION_FUNCTION)) write(*,*) density,PARTITION_FUNCTION
      END DO
      POPULATION=POPULATION*DENSITY/PARTITION_FUNCTION
      !WRITE(*,*) "LTE"
      !WRITE(*,*) population
      DO ILEVEL=1,NLEVEL
         IF (isnan(population(ilevel))) write(*,*) ilevel,PARTITION_FUNCTION,density
      END DO
   !  Check that the sum of the level populations matches the total density to within 0.1%
      IF(ABS(DENSITY-SUM(POPULATION))/DENSITY.GT.1.0D-3) THEN
         WRITE(*,"('ERROR! Sum of LTE level populations differs from the total density by',F4.1,'%')") &
            & 1.0D2*ABS(SUM(POPULATION)-DENSITY)/DENSITY
         STOP
      END IF

      RETURN
   END SUBROUTINE CALCULATE_LTE_POPULATIONS

   !=======================================================================
   !
   !  Calculate the line opacity of each coolant transition along each
   !  HEALPix ray to the PDR surface (or simulation boundary) for all
   !  particles. Calculations can be performed using either the Euler
   !  or Trapezium integration scheme (specified by a compiler flag).
   !
   !  UCLPDR major difference is we assume a single ray and instead of looking
   !  at particles along the ray, we assume homogenous medium at a distance
   !  of STEP_SIZE from edge
   !-----------------------------------------------------------------------
   SUBROUTINE CALCULATE_LINE_OPACITIES()
      IMPLICIT NONE
      !INTEGER    :: NRAYS=1 !hard coding 1 ray
      INTEGER :: N,ILEVEL,JLEVEL
      REAL(dp) :: STEP_SIZE,FACTOR1,FACTOR2,FACTOR3

      DO N=1,NCOOL ! Loop over coolants
         IF(coolants(N)%CONVERGED) CYCLE
         coolants(N)%OPACITY = 0.0D0
         DO ILEVEL=1,coolants(N)%NLEVEL ! Loop over levels (i)
            DO JLEVEL=1,coolants(N)%NLEVEL ! Loop over levels (j)
               IF(coolants(N)%A_COEFF(ILEVEL,JLEVEL).EQ.0) CYCLE
               !Factor 1 combines constants: = A_ij.c^3/8π.nu_ij^3
               FACTOR1 = (coolants(N)%A_COEFF(ILEVEL,JLEVEL)*C**3)/(8*PI*coolants(N)%FREQUENCY(ILEVEL,JLEVEL)**3) 
               
               !geometric distance between two particles
               !UCLPDR does this with particles along ray and get euler distance between them
               STEP_SIZE = CLOUD_SIZE
               ! STEP_SIZE=5.0d16
!                    Line opacity of the coolant transition (i,j) using the Trapezium integration method (default)

               !Factor 2 - inverse average LINEWIDTH between the two points
               !FACTOR2 = 2.0D0/(PARTICLE(L)%coolants(N)%LINEWIDTH+PARTICLE(M)%coolants(N)%LINEWIDTH) ! = 1/δv_D
               FACTOR2 = 1.0/coolants(N)%LINEWIDTH
  
               !Difference between average weight of ith level and average weight of jth level
               !divided by weight of jth level.
               ! = (n_j.g_i - n_i.g_j)/g_j = n_i.(n_j.g_i/n_i.g_j - 1)
               FACTOR3 = (coolants(N)%POPULATION(JLEVEL)*coolants(N)%WEIGHT(ILEVEL)  &
                     & -  coolants(N)%POPULATION(ILEVEL)*coolants(N)%WEIGHT(JLEVEL)) &
                     &  /coolants(N)%WEIGHT(JLEVEL)
                     
               coolants(N)%OPACITY(ILEVEL,JLEVEL) = coolants(N)%OPACITY(ILEVEL,JLEVEL) &
                                                       & + FACTOR1*FACTOR2*FACTOR3*STEP_SIZE 
               ! dtau_ij = A_ij.c^3/8π.nu_ij^3 * 1/δv_D * n_i.(n_j.g_i/n_i.g_j - 1) * dr
            END DO ! End of loop over levels (j)
         END DO ! End of loop over levels (i)
      END DO ! End of loop over coolants
   RETURN
   END SUBROUTINE CALCULATE_LINE_OPACITIES

!=======================================================================
!
!  Calculate the lambda operator for each allowed coolant transition,
!  needed to solve the level populations using the Accelerated Lambda
!  Iteration (ALI) method.
!
!  UCLPDR lambda operator seems to only take into accoutn change in opacity
!  between a point and surface. If we have homogenous medium then you get
!  zero. Do we need a different way to calculate this?
!-----------------------------------------------------------------------
SUBROUTINE CALCULATE_LAMBDA_OPERATOR()
   IMPLICIT NONE
   INTEGER :: N,i,j
   REAL(kind=dp) :: dTau_1=0.0,ALI_ij
   DO N=1,NCOOL ! Loop over coolants
      IF(coolants(N)%CONVERGED) CYCLE

      DO i=1,coolants(N)%nLevel
         DO j=1,coolants(N)%nLevel
            IF (coolants(N)%A_COEFF(i,j) .eq. 0.0) CYCLE

            !We want ifference in opacity between current "particle" and next
            !But I'm dealing with 1 particle so let's just set full opacity as difference between here and surface
            dTau_1=coolants(N)%OPACITY(I,J)

            IF (dTau_1 .ne.0.0) THEN   
               ALI_ij=2.0D0*(1.0D0-exp(-dTau_1))/dTau_1
               ALI_ij=1.0D0/(1.0D0+(ALI_ij+2.0)/(dTau_1)) - 1.0D0
            ELSE
               ALI_ij=0.0
            END IF
            coolants(N)%lambda(i,j)=ALI_ij+1.0D0
         END DO
      END DO
   END DO
   RETURN
END SUBROUTINE CALCULATE_LAMBDA_OPERATOR



!=======================================================================
!
!  Calculate the level populations of a given coolant by constructing
!  the matrix of transition rates and solving to find the populations
!  assuming statistical equilibrium. The resulting set of statistical
!  equilibrium equations take the form:
!
!     n_i.∑_j R_ij = ∑_j n_j.R_ji
!  or:
!     n_i.∑_j R_ij - ∑_j n_j.R_ji = 0
!
!  where n_i is the population density (cm^-3) of level i and R_ij is
!  the transition rate (s^-1) from level i to level j. By rearranging
!  these equations, they can then be put into the matrix form: A.n=0,
!  where A is a coefficient matrix and n is a vector containing the N
!  population densities of all the levels. The right-hand side of the
!  matrix equation is then a vector of length N, composed of zeroes.
!
!  The elements of the coefficient matrix A are specified as follows:
!
!     A_ij = -R_ji    (j≠i)
!     A_ii = ∑_j R_ij (j≠i)
!
!  Since this set of equilibrium equations is not independent, one of
!  the equations has to be replaced by the conservation equation:
!
!     ∑_j n_j = n_tot
!
!  where n_tot is the density of the coolant species in all levels.
!
!  Therefore, the last row of the coefficient matrix is replaced with
!  this summation over all levels, and the last right-hand-side value
!  is set to the total density of the coolant species (cm^-3).
!
!  This system of linear equations is then solved using Gauss-Jordan
!  elimination to determine the population densities of all N levels.
!
!-----------------------------------------------------------------------
SUBROUTINE CALCULATE_LEVEL_POPULATIONS(COOLANT,GasTemperature,gasDensity,abundances)
   IMPLICIT NONE
   REAL(dp), INTENT(IN) :: gasDensity,gasTemperature,abundances(:)
   TYPE(COOLANT_TYPE),  INTENT(INOUT)    :: COOLANT
   INTEGER:: I,J,N
   REAL(dp)     :: SUM
   REAL(dp), ALLOCATABLE :: R(:,:),A(:,:),B(:)

!  Allocate the arrays
   ALLOCATE(R(1:COOLANT%NLEVEL,1:COOLANT%NLEVEL))
   ALLOCATE(A(1:COOLANT%NLEVEL,1:COOLANT%NLEVEL))
   ALLOCATE(B(1:COOLANT%NLEVEL))
   A=0.0D0
   B=0.0D0
   R=0.0D0
   ! write(*,*) "temp",gasTemperature
   ! write(*,*) "pop",coolant%population

!  Construct the matrix of transition rates R_ij (s^-1)
   CALL CONSTRUCT_TRANSITION_MATRIX(COOLANT,R,gasTemperature,gasDensity&
      &,abundances)


!  Fill the coefficient matrix A and the right-hand-side vector b
   DO I=1,COOLANT%NLEVEL
      SUM=0.0D0
      DO J=1,COOLANT%NLEVEL
         IF(J.EQ.I) CYCLE
         IF(ISNAN(R(J,I))) THEN
            write(*,*) coolant%name
             DO N=1,COOLANT%NLEVEL
               write(*,*) R(N,:)
            END DO
            STOP
         END IF
         ! IF (R(J,I) .eq. 0.0) write(*,*) coolant%name,I,J
         A(I,J) = -R(J,I)
         SUM=SUM+R(I,J)
      END DO
      A(I,I)=SUM
      IF (A(I,I) .eq. 0.0) write(*,*) coolant%name,I
   END DO
   B=0.0D0

!  Replace the last equilibrium equation in the transition matrix with
!  the conservation equation (i.e. the sum of the population densities
!  over all levels), and replace the last entry in the right-hand-side
!  vector with the total density of the coolant species.
   A(COOLANT%NLEVEL,:)=1 ! Sum over all levels, ∑_j n_j
   B(COOLANT%NLEVEL)=coolant%density

!  Call the Gauss-Jordan solver (the solution is returned in vector b)
   CALL GAUSS_JORDAN(COOLANT%NLEVEL,A,B)

!  Replace negative or NaN values caused by numerical noise around zero
   DO I=1,COOLANT%NLEVEL
      IF(.NOT.B(I).GE.0) B(I)=0.0D0
      IF (ISNAN(B(I))) THEN
         write(*,*) I, B(I),COOLANT%NLEVEL
         write(*,*) A
      END IF
   END DO

!  Store the previously-calculated population densities for comparison
   coolant%PREVIOUS_POPULATION=coolant%POPULATION

!  Store the new population densities
   coolant%POPULATION=B

   DEALLOCATE(R,A,B)
         
   RETURN
END SUBROUTINE CALCULATE_LEVEL_POPULATIONS
!=======================================================================

!=======================================================================
!
!  Standard linear equation solver using Gauss-Jordon elimination taken
!  directly from Numerical Recipes (Chapter B2).
!
!  A is an NxN input coefficient matrix, B is an input vector of size N
!  containing the right-hand-side values. On output, A is replaced by its
!  matrix inverse and B is replaced by the corresponding set of solution
!  values.
!
!-----------------------------------------------------------------------
SUBROUTINE GAUSS_JORDAN(N,A,B)
   USE SWAP_FUNCTION
   IMPLICIT NONE

   INTEGER, INTENT(IN)    :: N
   REAL(dp),     INTENT(INOUT) :: A(:,:)
   REAL(dp),     INTENT(INOUT) :: B(:)

   INTEGER :: I,J,K,L,IROW,ICOL
   INTEGER, ALLOCATABLE :: IPIV(:),INDEX_ROW(:),INDEX_COL(:)
   REAL(dp) :: MAX,DUMMY,PIVINV

   ALLOCATE(IPIV(1:N),INDEX_ROW(1:N),INDEX_COL(1:N))

   ICOL=0
   IROW=0
   IPIV=0

   DO I=1,N ! Main loop over columns to be reduced
      MAX=0.0D0
      DO J=1,N
         IF(IPIV(J).NE.1) THEN
            DO K=1,N
               IF(IPIV(K).EQ.0) THEN
                  IF(ISNAN(A(J,K))) THEN
                     WRITE(6,*)
                     WRITE(6,*) 'ERROR! NaN found in coefficient matrix A of Gauss-Jordan routine'
                     WRITE(6,"(' A(',I3,',',I3,') =',F4.1)") J,K,A(J,K)
                     WRITE(6,*) 'A ='
                     DO L=1,N
                        WRITE(6,"(100(0PD9.1))") A(L,:)
                     END DO
                     WRITE(6,*)
                     STOP
                  END IF
                  IF(ABS(A(J,K)).GE.MAX) THEN
                     MAX=ABS(A(J,K))
                     IROW=J
                     ICOL=K
                  END IF
               ELSE IF(IPIV(K).GT.1) THEN
                  WRITE(6,*)
                  WRITE(6,*) 'ERROR! Singular matrix found in Gauss-Jordan routine (#1)'
                  WRITE(6,*)
                  STOP
               END IF
            END DO
         END IF
      END DO
      IPIV(ICOL)=IPIV(ICOL)+1
      IF(IROW.NE.ICOL) THEN
         CALL SWAP(A(IROW,:),A(ICOL,:))
         CALL SWAP(B(IROW),B(ICOL))
      END IF
      INDEX_ROW(I)=IROW
      INDEX_COL(I)=ICOL
      IF(A(ICOL,ICOL).EQ.0.0D0) THEN
         WRITE(6,*)
         WRITE(6,*) 'ERROR! Singular matrix found in Gauss-Jordan routine (#2)'
         WRITE(6,*)
         STOP
      END IF
      PIVINV=1.0D0/A(ICOL,ICOL)
      A(ICOL,ICOL)=1.0D0
      A(ICOL,:)=A(ICOL,:)*PIVINV
      B(ICOL)=B(ICOL)*PIVINV
      DO L=1,N
         IF(L.NE.ICOL) THEN
            DUMMY=A(L,ICOL)
            A(L,ICOL)=0.0D0
            A(L,:)=A(L,:)-A(ICOL,:)*DUMMY
            B(L)=B(L)-B(ICOL)*DUMMY
         END IF
      END DO
   END DO ! End of main loop over columns

!  Unscramble the solution by interchanging pairs of columns
!  in the reverse order to which the permutation was built up
   DO L=N,1,-1
      CALL SWAP(A(:,INDEX_ROW(L)),A(:,INDEX_COL(L)))
   END DO

   DEALLOCATE(IPIV,INDEX_ROW,INDEX_COL)

   RETURN
END SUBROUTINE GAUSS_JORDAN

!=======================================================================
!=======================================================================
!
!  Construct the matrix of transition rates for a given coolant species
!  using the escape probability formalism to determine the local field.
!
!-----------------------------------------------------------------------
SUBROUTINE CONSTRUCT_TRANSITION_MATRIX(COOLANT,TRANSITION_MATRIX,GasTemperature,gasDensity,abundances)
   IMPLICIT NONE
   TYPE(COOLANT_TYPE), INTENT(INOUT)    :: COOLANT
   REAL(dp), INTENT(IN)    :: GasTemperature,gasDensity,abundances(:)
   REAL(dp), INTENT(OUT)   :: TRANSITION_MATRIX(:,:)

   INTEGER :: I,J,K
   REAL(dp) :: TPOP,RHO_GRAIN,DUST_EMISSIVITY
   REAL(dp) :: S_ij,B_ij,B_ij_CMB,B_ij_DUST,BETA_ij
   REAL(dp), ALLOCATABLE :: RADIATION_FIELD(:,:)
   REAL(dp), ALLOCATABLE :: COLLISIONAL_RATE(:,:)

   REAL(dp) :: S_ij_PREVIOUS,LAMBDA_ij

!  Allocate and initialize the mean integrated radiation field
   ALLOCATE(RADIATION_FIELD(1:COOLANT%NLEVEL,1:COOLANT%NLEVEL))
   RADIATION_FIELD=0.0D0
   ! TRANSITION_MATRIX=0.0D0
!  Initialize the local emissivities
   COOLANT%EMISSIVITY=0.0D0

   DO I=1,COOLANT%NLEVEL
      DO J=1,COOLANT%NLEVEL
         IF(J.GE.I) EXIT

!        Calculate the background radiation field including contributions
!        from CMB blackbody emission and dust modified blackbody emission
         IF(COOLANT%FREQUENCY(I,J).LT.1.0D15) THEN
            B_ij_CMB=(2*HP*COOLANT%FREQUENCY(I,J)**3)/(C**2)/(EXP(HP*COOLANT%FREQUENCY(I,J)/(K_BOLTZ*T_CMB))-1.0D0)
!#ifdef DUST
!            RHO_GRAIN=2.0D0 ! Grain mass density (g cm^-3)
!            DUST_EMISSIVITY=(RHO_GRAIN*DUST_DENSITY)*(0.01*(1.3*COOLANT%FREQUENCY(I,J)/3.0D11))
!            B_ij_DUST=(2*HP*COOLANT%FREQUENCY(I,J)**3)/(C**2)/(EXP(HP*COOLANT%FREQUENCY(I,J)/(KB*DUST_TEMPERATURE))-1.D0)*DUST_EMISSIVITY
!#else
            B_ij_DUST=0.0D0
!#endif
            B_ij=B_ij_CMB+B_ij_DUST
         ELSE
            B_ij=0.0D0
         END IF

!        If the population n_i is zero then the source function is zero and the
!        mean integrated radiation field is just the background radiation field
         IF(COOLANT%POPULATION(I).lt. 1.0d-20) THEN
            S_ij=0.0D0
            BETA_ij=1.0D0

!        If the difference between n_i.g_j and n_j.g_i is vanishingly small, then
!        calculate the source function directly and set the escape probability to 1
         ELSE IF(ABS(COOLANT%POPULATION(I)*COOLANT%WEIGHT(J)-COOLANT%POPULATION(J)*COOLANT%WEIGHT(I)).EQ.0) THEN
            S_ij=HP*COOLANT%FREQUENCY(I,J)*COOLANT%POPULATION(I)*COOLANT%A_COEFF(I,J)/(4.0*PI)
            BETA_ij=1.0D0
           ! write(*,*) "vanishingly",HP,COOLANT%FREQUENCY(I,J),COOLANT%POPULATION(I),COOLANT%A_COEFF(I,J),(4.0*PI)
         ELSE
!           Calculate the source function
            S_ij=(2*HP*COOLANT%FREQUENCY(I,J)**3)/(C**2) &
              & /((COOLANT%POPULATION(J)*COOLANT%WEIGHT(I)) &
              &  /(COOLANT%POPULATION(I)*COOLANT%WEIGHT(J))-1.0D0)

!           Calculate the escape probability
            BETA_ij=ESCAPE_PROBABILITY(COOLANT%OPACITY(I,J))
            !write(*,*) "regular",HP,COOLANT%FREQUENCY(I,J),COOLANT%POPULATION(J),COOLANT%WEIGHT(I),COOLANT%POPULATION(I),COOLANT%WEIGHT(J)
         END IF

!        Calculate the local emissivity (erg cm^-3 s^-1) for the transition line
         IF(COOLANT%POPULATION(I).GT.0) THEN
            COOLANT%EMISSIVITY(I,J)=COOLANT%POPULATION(I)*COOLANT%A_COEFF(I,J) &
                                             & *HP*COOLANT%FREQUENCY(I,J)*BETA_ij*(S_ij-B_ij)/S_ij
         END IF

! !        Calculate the mean integrated radiation field <J_ij>
!          RADIATION_FIELD(I,J)=(1.0D0-BETA_ij)*S_ij+BETA_ij*B_ij
!          RADIATION_FIELD(J,I)=RADIATION_FIELD(I,J)
!          if (isnan(RADIATION_FIELD(I,J))) THEN
!             WRITE(*,*) "beta_ij", BETA_ij
!             write(*,*) "s_ij",S_ij
!             write(*,*) "b_ij",B_ij
!             write(*,*) "freq",COOLANT%FREQUENCY(I,J)
!             write(*,*) "pop j" ,COOLANT%POPULATION(J)
!             write(*,*) "weight i",COOLANT%WEIGHT(I)
!             write(*,*) "pop i",COOLANT%POPULATION(I)
!             write(*,*) "weight j",COOLANT%WEIGHT(J)
!          END IF
!        Lambda operator keeps breaking so lets not use ALI for now
! !        Calculate the source function in the same manner for
! !        the populations determined on the previous iteration
         IF(COOLANT%PREVIOUS_POPULATION(I).lt.1.d-20) THEN
            S_ij_PREVIOUS=0.0D0
         ELSE IF(ABS(COOLANT%PREVIOUS_POPULATION(I)*COOLANT%WEIGHT(J)-COOLANT%PREVIOUS_POPULATION(J)*COOLANT%WEIGHT(I)).EQ.0) THEN
            S_ij_PREVIOUS=HP*COOLANT%FREQUENCY(I,J)*COOLANT%PREVIOUS_POPULATION(I)*COOLANT%A_COEFF(I,J)/(4.0*PI)
         ELSE
            S_ij_PREVIOUS=(2*HP*COOLANT%FREQUENCY(I,J)**3)/(C**2) &
              & /((COOLANT%PREVIOUS_POPULATION(J)*COOLANT%WEIGHT(I)) &
              &  /(COOLANT%PREVIOUS_POPULATION(I)*COOLANT%WEIGHT(J))-1.0D0)
         END IF

!        Use the Accelerated Lambda Iteration method to speed up convergence
!        by amplifying the incremental difference of the new source function
         LAMBDA_ij=COOLANT%LAMBDA(I,J)
         RADIATION_FIELD(I,J)=(1.0D0-BETA_ij)*(LAMBDA_ij*S_ij+(1.0D0-LAMBDA_ij)*S_ij_PREVIOUS)+BETA_ij*B_ij
         ! if (isnan(RADIATION_FIELD(I,J))) THEN
         !    WRITE(*,*) "beta",BETA_ij
         !    WRITE(*,*) "lambda",LAMBDA_ij
         !    WRITE(*,*) "source and prev",S_ij,S_ij_PREVIOUS
         !    WRITE(*,*) "B",B_ij
         !    WRITE(*,*) "opacity",COOLANT%OPACITY(I,J)
         !    write(*,*) "pop i",COOLANT%POPULATION(I)

         ! END IF

         RADIATION_FIELD(J,I)=RADIATION_FIELD(I,J)
      END DO ! End of loop over levels (j)
   END DO ! End of loop over levels (i)

!  Allocate and calculate the collisional rates
   ALLOCATE(COLLISIONAL_RATE(1:COOLANT%NLEVEL,1:COOLANT%NLEVEL))
   CALL CALCULATE_COLLISIONAL_RATES(COOLANT,gasDensity,GasTemperature,abundances,COLLISIONAL_RATE)

!  Construct the transition matrix: R_ij = A_ij + B_ij.<J> + C_ij
   DO I=1,COOLANT%NLEVEL
      DO J=1,COOLANT%NLEVEL
         if (isnan(coolant%A_COEFF(I,J))) WRITE(*,*) "A_COEFF!"
         if (isnan(coolant%b_COEFF(I,J))) WRITE(*,*) "b_COEFF!"
         if (isnan(RADIATION_FIELD(I,J))) WRITE(*,*) "RADFIELD!"
         if (isnan(COLLISIONAL_RATE(I,J))) WRITE(*,*) "COLLS!"
         TRANSITION_MATRIX(I,J)=COOLANT%A_COEFF(I,J) &
                            & + COOLANT%B_COEFF(I,J)*RADIATION_FIELD(I,J) &
                            & + COLLISIONAL_RATE(I,J)
      END DO
   END DO

   DEALLOCATE(RADIATION_FIELD,COLLISIONAL_RATE)

   RETURN
END SUBROUTINE CONSTRUCT_TRANSITION_MATRIX
!=======================================================================


!=======================================================================
!
!  Calculate the total collisional rates (s^-1) for a given coolant at
!  the specified temperature by summing the individual rates from each
!  of the available collision partners by linear interpolation between
!  the available rate coefficients at specific temperature values.
!
!-----------------------------------------------------------------------
SUBROUTINE CALCULATE_COLLISIONAL_RATES(COOLANT,DENSITY,TEMPERATURE, &
                                     & ABUNDANCE,COLLISIONAL_RATE)
  ! USE FUNCTIONS_MODULE

   IMPLICIT NONE

   TYPE(COOLANT_TYPE), INTENT(IN)  :: COOLANT
   REAL(dp),      INTENT(IN)  :: DENSITY,TEMPERATURE
   REAL(dp),      INTENT(IN)  :: ABUNDANCE(:)
   REAL(dp),      INTENT(OUT) :: COLLISIONAL_RATE(:,:)

   INTEGER:: I,J,K,KLO,KHI,PARTNER_ID
   REAL(dp) :: PARA_FRACTION,ORTHO_FRACTION
   REAL(dp) :: STEP,C_COEFF

!  Initialize the collisional rates
   COLLISIONAL_RATE=0.0D0

!  Calculate the H2 ortho/para ratio at equilibrium for the specified
!  temperature and the resulting fractions of H2 in para & ortho form
   PARA_FRACTION=1.0D0/(1.0D0+ORTHO_PARA_RATIO(TEMPERATURE))
   ORTHO_FRACTION=1.0D0-PARA_FRACTION


   DO PARTNER_ID=1,7 ! Loop over collision partners
!     Skip the collision partner if no rates are available
      IF(COOLANT%TEMPERATURE(PARTNER_ID,1).EQ.0.0D0) CYCLE

!     Determine the two nearest temperature values
!     present within the list of collisional rates
      KLO=0; KHI=0
      DO K=1,COOLANT%NTEMP ! Loop over temperatures
         IF(COOLANT%TEMPERATURE(PARTNER_ID,K).GT.TEMPERATURE) THEN
            KLO=K-1
            KHI=K
            EXIT
         ELSE IF(COOLANT%TEMPERATURE(PARTNER_ID,K).EQ.0.0D0) THEN
            KLO=K-1
            KHI=K-1
            EXIT
         END IF
      END DO

!     If the required temperature is above or below the range of available
!     temperature values then use the highest or lowest value in the range
      IF(KHI.EQ.0) THEN
         KLO=COOLANT%NTEMP
         KHI=COOLANT%NTEMP
      ELSE IF(KHI.EQ.1) THEN
         KLO=1
         KHI=1
      END IF

!     Calculate the "distance" between the two temperature
!     values, to be used in the linear interpolation below
      IF(KLO.EQ.KHI) THEN
         STEP=0.0D0
      ELSE
         STEP=(TEMPERATURE-COOLANT%TEMPERATURE(PARTNER_ID,KLO)) &
           & /(COOLANT%TEMPERATURE(PARTNER_ID,KHI)-COOLANT%TEMPERATURE(PARTNER_ID,KLO))
      END IF

!     Linearly interpolate the collisional rate coefficients
!     for each collision partner at the required temperature
      IF(PARTNER_ID.EQ.1) THEN ! Collisions with H2
         DO I=1,COOLANT%NLEVEL
            DO J=1,COOLANT%NLEVEL
               C_COEFF=COOLANT%C_COEFF(PARTNER_ID,I,J,KLO)+(COOLANT%C_COEFF(PARTNER_ID,I,J,KHI)-COOLANT%C_COEFF(PARTNER_ID,I,J,KLO))*STEP
               COLLISIONAL_RATE(I,J)=COLLISIONAL_RATE(I,J)+C_COEFF*DENSITY*ABUNDANCE(nH2)
            END DO
         END DO
      END IF
      IF(PARTNER_ID.EQ.2) THEN ! Collisions with para-H2
         DO I=1,COOLANT%NLEVEL
            DO J=1,COOLANT%NLEVEL
               C_COEFF=COOLANT%C_COEFF(PARTNER_ID,I,J,KLO)+(COOLANT%C_COEFF(PARTNER_ID,I,J,KHI)-COOLANT%C_COEFF(PARTNER_ID,I,J,KLO))*STEP
               COLLISIONAL_RATE(I,J)=COLLISIONAL_RATE(I,J)+C_COEFF*DENSITY*ABUNDANCE(nH2)*PARA_FRACTION
            END DO
         END DO
      END IF
      IF(PARTNER_ID.EQ.3) THEN ! Collisions with ortho-H2
         DO I=1,COOLANT%NLEVEL
            DO J=1,COOLANT%NLEVEL
               C_COEFF=COOLANT%C_COEFF(PARTNER_ID,I,J,KLO)+(COOLANT%C_COEFF(PARTNER_ID,I,J,KHI)-COOLANT%C_COEFF(PARTNER_ID,I,J,KLO))*STEP
               COLLISIONAL_RATE(I,J)=COLLISIONAL_RATE(I,J)+C_COEFF*DENSITY*ABUNDANCE(nH2)*ORTHO_FRACTION
            END DO
         END DO
      END IF
      IF(PARTNER_ID.EQ.4) THEN ! Collisions with electrons
         DO I=1,COOLANT%NLEVEL
            DO J=1,COOLANT%NLEVEL
               C_COEFF=COOLANT%C_COEFF(PARTNER_ID,I,J,KLO)+(COOLANT%C_COEFF(PARTNER_ID,I,J,KHI)-COOLANT%C_COEFF(PARTNER_ID,I,J,KLO))*STEP
               COLLISIONAL_RATE(I,J)=COLLISIONAL_RATE(I,J)+C_COEFF*DENSITY*ABUNDANCE(nelec)
            END DO
         END DO
      END IF
      IF(PARTNER_ID.EQ.5) THEN ! Collisions with H
         DO I=1,COOLANT%NLEVEL
            DO J=1,COOLANT%NLEVEL
               C_COEFF=COOLANT%C_COEFF(PARTNER_ID,I,J,KLO)+(COOLANT%C_COEFF(PARTNER_ID,I,J,KHI)-COOLANT%C_COEFF(PARTNER_ID,I,J,KLO))*STEP
               COLLISIONAL_RATE(I,J)=COLLISIONAL_RATE(I,J)+C_COEFF*DENSITY*ABUNDANCE(nH)
            END DO
         END DO
      END IF
      IF(PARTNER_ID.EQ.6) THEN ! Collisions with He
         DO I=1,COOLANT%NLEVEL
            DO J=1,COOLANT%NLEVEL
               C_COEFF=COOLANT%C_COEFF(PARTNER_ID,I,J,KLO)+(COOLANT%C_COEFF(PARTNER_ID,I,J,KHI)-COOLANT%C_COEFF(PARTNER_ID,I,J,KLO))*STEP
               COLLISIONAL_RATE(I,J)=COLLISIONAL_RATE(I,J)+C_COEFF*DENSITY*ABUNDANCE(nHe)
            END DO
         END DO
      END IF
      IF(PARTNER_ID.EQ.7) THEN ! Collisions with protons
         DO I=1,COOLANT%NLEVEL
            DO J=1,COOLANT%NLEVEL
               C_COEFF=COOLANT%C_COEFF(PARTNER_ID,I,J,KLO)+(COOLANT%C_COEFF(PARTNER_ID,I,J,KHI)-COOLANT%C_COEFF(PARTNER_ID,I,J,KLO))*STEP
               COLLISIONAL_RATE(I,J)=COLLISIONAL_RATE(I,J)+C_COEFF*DENSITY*ABUNDANCE(nHx)
            END DO
         END DO
      END IF

   END DO ! End of loop over collision partners

   RETURN
END SUBROUTINE CALCULATE_COLLISIONAL_RATES

!=======================================================================
!
!  Check for convergence in the population densities of all levels of
!  each coolant of each particle. Set the relevant convergence flags.
!  Calculate the percentage of particles that have converged for each
!  coolant.
!
!-----------------------------------------------------------------------
LOGICAL FUNCTION CHECK_CONVERGENCE()
   INTEGER :: I,N
   REAL(dp) :: RELATIVE_CHANGE
   REAL(dp), PARAMETER :: POPULATION_LIMIT=1.0D-12,POPULATION_CONVERGENCE_CRITERION=0.0001!pretty sure this is 0.1 in uclpdr
   LOGICAL :: convergence(NCOOL)

      DO N=1,NCOOL ! Loop over coolants
         coolants(N)%CONVERGED=.True.
         DO I=1,coolants(N)%NLEVEL ! Loop over levels

!           Skip this level if its population density is below the cut-off
            IF(coolants(N)%POPULATION(I).LT.POPULATION_LIMIT*coolants(N)%DENSITY) CYCLE

!           Skip this level if its population density has not changed
            IF(coolants(N)%POPULATION(I).EQ.coolants(N)%PREVIOUS_POPULATION(I)) CYCLE

!           Calculate the relative change in population density between this iteration and the previous
            RELATIVE_CHANGE=ABS(coolants(N)%POPULATION(I)-coolants(N)%PREVIOUS_POPULATION(I)) &
                          & *2/(coolants(N)%POPULATION(I)+coolants(N)%PREVIOUS_POPULATION(I))

!           If the relative change is greater than the criterion for convergence, set the flag to false
            IF(RELATIVE_CHANGE.GT.POPULATION_CONVERGENCE_CRITERION) THEN
               coolants(N)%CONVERGED=.FALSE.
               EXIT
            END IF
         END DO ! End of loop over levels
         convergence(N)=coolants(N)%converged
      END DO ! End of loop over coolant
   CHECK_CONVERGENCE=ALL(convergence)
END FUNCTION CHECK_CONVERGENCE


!=======================================================================
!
!  Calculate the ortho-to-para ratio of H2 at thermal equilibrium
!  for the specified temperature, making use of energy level data
!  if available, or an approximation if not.
!
!-----------------------------------------------------------------------
FUNCTION ORTHO_PARA_RATIO(TEMPERATURE)
   IMPLICIT NONE
   REAL(dp) :: ORTHO_PARA_RATIO
   REAL(dp), INTENT(IN) :: TEMPERATURE

   INTEGER :: I,J,N,ORTHO_INDEX,PARA_INDEX
   REAL(dp) :: I_ORTHO,I_PARA,ORTHO_FRACTION,PARA_FRACTION

!  Check if coolant data is available for the ortho and para forms
   ORTHO_INDEX=0; PARA_INDEX=0
   DO N=1,NCOOL
      IF(COOLANTS(N)%NAME.EQ."o-H2") ORTHO_INDEX=N
      IF(COOLANTS(N)%NAME.EQ."p-H2") PARA_INDEX=N
   END DO

!  Calculate the exact ortho/para ratio if molecular data is available
   IF(ORTHO_INDEX.NE.0 .AND. PARA_INDEX.NE.0) THEN

      I_ORTHO=1.0D0; I_PARA=0.0D0 ! Total nuclear spins of the two forms

!     Calculate the ortho/para ratio of H2 using the expression
!     from Poelman & Spaans (2005, A&A, 440, 559, equation 11)
      ORTHO_FRACTION=0.0D0
      DO I=1,COOLANTS(ORTHO_INDEX)%NLEVEL
         ORTHO_FRACTION=ORTHO_FRACTION+COOLANTS(ORTHO_INDEX)%WEIGHT(I) &
                     & *EXP(-COOLANTS(ORTHO_INDEX)%ENERGY(I)/(K_BOLTZ*TEMPERATURE))
      END DO

      PARA_FRACTION=0.0D0
      DO I=1,COOLANTS(PARA_INDEX)%NLEVEL
         PARA_FRACTION=PARA_FRACTION+COOLANTS(PARA_INDEX)%WEIGHT(I) &
                    & *EXP(-COOLANTS(PARA_INDEX)%ENERGY(I)/(K_BOLTZ*TEMPERATURE))
      END DO

      ORTHO_PARA_RATIO=(2*I_ORTHO+1)/(2*I_PARA+1)*(ORTHO_FRACTION/PARA_FRACTION)

   ELSE

!     Approximate the ortho/para ratio of H2 using the expression
!     given by Flower & Watt (1985, MNRAS, 213, 991, equation 2)
      ORTHO_PARA_RATIO=9.0D0*EXP(-170.5D0/TEMPERATURE)

!     Limit the ortho/para ratio to its statistical limit
      IF(ORTHO_PARA_RATIO.GT.3.0D0) ORTHO_PARA_RATIO=3.0D0

   END IF

END FUNCTION ORTHO_PARA_RATIO

FUNCTION ESCAPE_PROBABILITY(TAU) RESULT(BETA)
   IMPLICIT NONE
   REAL(dp),  INTENT(IN) :: TAU
   INTEGER :: K
   REAL(dp)  :: BETA

!  Initialize the escape probability values along each ray
   BETA=0.0D0

!     Limit the escape probability to unity for masing transitions (tau <= 0)
   IF(TAU.LE.0) THEN
      BETA=1.0D0

!     Prevent floating point overflow caused by very low opacity (tau < 1E-8)
   ELSE IF(ABS(TAU).LT.1.0D-8) THEN
      BETA=1.0D0

!     For all other cases use the standard escape probability formalism
   ELSE
      BETA=(1.0D0-EXP(-TAU))/TAU
   END IF

!  The total escape probability must be divided by the number of rays to
!  account for the fraction of the total solid angle covered by each ray
!  (assuming that each ray covers the same fraction of the total 4π sr).
!  In the case of only 1 ray (i.e., semi-infinite slab geometry) the ray
!  subtends a solid angle of 2π sr, since the photons escape through the
!  hemisphere in the outward direction, so its escape probability should
!  be divided by two.
   BETA=0.5*BETA
   RETURN
END FUNCTION ESCAPE_PROBABILITY

END MODULE COOLANT_MODULE