!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Module that provides heating and cooling rates		  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE heating
    USE NETWORK
    USE COOLANT_MODULE
    USE CONSTANTS
IMPLICIT NONE



CONTAINS

    SUBROUTINE initializeHeating(gasTemperature, gasDensity,abundances,size)
        REAL(dp), INTENT(in) :: gasTemperature,gasDensity,size
        REAL(dp), INTENT(in) :: abundances(:)
        INTEGER ::i,j
        CALL READ_COOLANTS()
        DO i=1,ncool
            DO j=1,nspec
                if (coolantNames(i) .eq. specName(j)) coolantIndices(i)=j
            END DO
        END DO

        CALL UPDATE_COOLANT_LINEWIDTHS(gasTemperature)
        DO i=1,NCOOL
            coolants(i)%DENSITY=abundances(coolantIndices(i))*gasDensity
            CALL CALCULATE_LTE_POPULATIONS(coolants(i)%NLEVEL,coolants(i)%ENERGY,coolants(i)%WEIGHT, &
                                          & coolants(i)%DENSITY,gasTemperature, &
                                          & coolants(i)%POPULATION)
        END DO
        CLOUD_SIZE=size
        ! coolantIndices(ncool+1)=nspec
        !Here I could find the reactions that give chemical heating and cooling
        !and store species indices in an array
    END SUBROUTINE initializeHeating


    REAL(dp) FUNCTION getTempDot(gasTemperature,gasDensity,habingField,abundances,h2dis,zeta,cIonRate,dustAbundance&
                                &,exoReactants1,exoReactants2,exoRates,exothermicities)
        !Habing field is radfield diminished by Av
        REAL(dp), INTENT(in) :: gasTemperature,gasDensity,habingField,h2dis,zeta,cIonRate,dustAbundance
        REAL(dp), INTENT(in) :: abundances(:),exoReactants1(:),exoReactants2(:),exoRates(:),exothermicities(:)
        REAL(dp) adiabaticIdx,heating,cooling

        !First calculate adiabatic index - should use number density but that's just an additional common factor
        adiabaticIdx=5.0*(abundances(nh)+abundances(nhe)+abundances(nelec)+abundances(nh2))+2.0*abundances(nh2)
        adiabaticIdx=adiabaticIdx/(3.0*(abundances(nh)+abundances(nhe)+abundances(nelec)+abundances(nh2))+2.0*abundances(nh2))

        !then calculate overall heating/cooling rate
        heating=getHeatingRate(gasTemperature,gasDensity,habingField,abundances,h2dis,zeta,cIonRate,dustAbundance&
            &,exoReactants1,exoReactants2,exoRates,exothermicities)
        !write(*,*) "Total Heating", heating

        IF (gasTemperature .gt. 5.0) THEN
            cooling=getCoolingRate(gasTemperature,gasDensity,abundances,h2dis)
        ELSE
            Cooling=0.0
        END IF

        !write(*,*) "Total Cooling", Cooling

        getTempDot=heating-cooling
        !write(*,*) "Temp Dot",getTempDot
         !and convert to dT/dt
         getTempDot=((adiabaticIdx-1)*getTempDot)/(K_BOLTZ*gasDensity)


    END FUNCTION getTempDot


    REAL(dp) FUNCTION getHeatingRate(gasTemperature,gasDensity,habingField,abundances,h2dis,zeta,cIonRate,dustAbundance&
                                    &,exoReactants1,exoReactants2,exoRates,exothermicities)
        REAL(dp), INTENT(in) :: gasTemperature,gasDensity,habingField,h2dis,zeta,cIonRate,dustAbundance
        REAL(dp), INTENT(IN):: abundances(:),exoReactants1(:),exoReactants2(:),exoRates(:),exothermicities(:)
        REAL(dp) :: heatingMode

            heatingMode=photoelectricHeating(gasTemperature,gasDensity,habingField,abundances(nelec))
            !write(*,*) heatingMode, "photoelectric"
            getHeatingRate=heatingMode

            heatingMode=H2FormationHeating(gasTemperature,gasDensity,abundances(nh))
            !write(*,*) heatingMode, "H2 Form"
            getHeatingRate=getHeatingRate+heatingMode

            heatingMode=h2FUVPumpHeating(abundances(nh),abundances(nh2),gasTemperature,gasDensity,h2dis)
            !write(*,*) heatingMode, "H2 FUV pump"
            getHeatingRate=getHeatingRate+heatingMode

            heatingMode=H2PhotodisHeating(gasDensity,abundances(nh2),h2dis)
            !write(*,*) heatingMode, "H2 Dis"
            getHeatingRate=getHeatingRate+heatingMode
          
            heatingMode=CarbonIonizationHeating(cIonRate,abundances(nc),gasDensity)
            !write(*,*) heatingMode, "C ionization"
            getHeatingRate=getHeatingRate+heatingMode

            heatingMode=cosmicRayHeating(zeta,gasDensity,abundances(nh2))
            !write(*,*) heatingMode, "CR Heating"
            getHeatingRate=getHeatingRate+heatingMode

            heatingMode=chemicalHeating(gasDensity,exoReactants1,exoReactants2,exoRates,exothermicities)
            getHeatingRate=getHeatingRate+heatingMode
          
            heatingMode=gasGrainCollisions(gasTemperature,gasDensity,dustAbundance)
            getHeatingRate=getHeatingRate+heatingMode
    END FUNCTION getHeatingRate

    REAL(dp) FUNCTION getCoolingRate(gasTemperature,gasDensity,abundances,h2dis)
        REAL(dp), INTENT(IN) :: gasTemperature,gasDensity,h2dis
        REAL(dp), INTENT(IN) :: abundances(:)
        real(dp) :: coolingMode, coolings(5)
        INTEGER :: ti


        ! coolingMode=atomicCooling(gasTemperature,gasDensity,abundances(nh),abundances(nhe),&
        !                         &abundances(nelec),abundances(nhx),abundances(nhex))
        ! write(*,*) coolingMode, "atomic"
        ! getCoolingRate=coolingMode

        ! coolingMode=collionallyInducedEmission(gasTemperature,gasDensity,abundances(nh2))
        ! write(*,*) coolingMode, "CIE"
        ! getCoolingRate=getCoolingRate+coolingMode
        
        ! coolingMode=comptonCooling(gasTemperature,gasDensity,abundances(nelec))
        ! write(*,*) coolingMode, "Compton"
        ! getCoolingRate=getCoolingRate+coolingMode

        ! coolingMode=continuumEmission(gasTemperature,gasDensity)
        ! write(*,*) coolingMode, "Continuum"
        ! getCoolingRate=getCoolingRate+coolingMode

        ! coolingMode=H2VibrationalCooling(gasTemperature,gasDensity,abundances(nh2),h2dis)
        ! write(*,*) coolingMode, "H2 Vibrational"
        ! getCoolingRate=getCoolingRate+coolingMode
        getCoolingRate=0.0D0
        DO ti=1,5
            coolings(ti)=lineCooling(gasTemperature,gasDensity,abundances)
        END DO
        call pair_insertion_sort(coolings)
        coolingMode=coolings(3)
        ! write(*,*) "Line cooling"
        ! write(*,*) coolings, coolingMode
        getCoolingRate=getCoolingRate+coolingMode
    END FUNCTION getCoolingRate


    !Put fits from Neufeld here.
    REAL(dp) FUNCTION simpleLineCooling(gasTemperature,gasDensity,abundances)
        REAL(dp), INTENT(IN) :: gasTemperature,gasDensity,abundances(:)
        simpleLineCooling=0.0
    END FUNCTION simpleLineCooling 


    REAL(dp) FUNCTION lineCooling(gasTemperature,gasDensity,abundances)
        REAL(dp), INTENT(IN) :: gasTemperature,gasDensity,abundances(:)
        INTEGER :: N,I!, collisionalIndices(5)=(/nh,nhx,nh2,nhe,nelec/)
        real(dp) :: moleculeCooling=0.0

        CALL UPDATE_COOLANT_LINEWIDTHS(gasTemperature)
        DO N=1,NCOOL
            coolants(N)%DENSITY=abundances(coolantIndices(N))*gasDensity
            CALL CALCULATE_LTE_POPULATIONS(coolants(N)%NLEVEL,coolants(N)%ENERGY,coolants(N)%WEIGHT, &
                                          & coolants(N)%DENSITY,gasTemperature, &
                                          & coolants(N)%POPULATION)
        END DO
        CALL CALCULATE_LINE_OPACITIES()
        CALL CALCULATE_LAMBDA_OPERATOR()

        !!write(*,*)  "lte done"
        !I should then do LVG interactions
         DO I=1,100!while not converged and less than 100 tries:
            DO N=1,NCOOL
                CALL CALCULATE_LEVEL_POPULATIONS(coolants(N),gasTemperature,gasDensity,&
                    &abundances)
            END DO
         !!write(*,*) "after Lvg",abundances(nh)

            !!write(*,*) I

            CALL CALCULATE_LINE_OPACITIES()
            CALL CALCULATE_LAMBDA_OPERATOR()
            IF (CHECK_CONVERGENCE()) EXIT
        END DO 

        !  Calculate the cooling rate due to the Lyman-alpha emission for each particle
        !  using the analytical expression of Spitzer (1978) neglecting photon trapping
        DO N=1,NCOOL
            IF(coolants(N)%NAME.EQ."H") THEN
               coolants(N)%EMISSIVITY(2,1) = 7.3D-19*(abundances(nelec)*gasDensity) &
                                                          & *(abundances(nH)*gasDensity) &
                                                          & *EXP(-118400.0D0/gasTemperature)
            END IF
        END DO

        !Calculate the cooling rates
        lineCooling=0.0
        DO N=1,NCOOL
            moleculeCooling=SUM(coolants(N)%EMISSIVITY,MASK=.NOT.ISNAN(coolants(N)%EMISSIVITY))

            ! !we get these wild changes in cooling rate so let's force it not to change too much in a timestep
            ! IF ( (abs(moleculeCooling-coolants(N)%previousCooling)/coolants(N)%previousCooling) .gt. 2.0d0) THEN
            !     !unless it's step 1, use old cooling for this time step
            !     IF (coolants(N)%previousCooling .ne. 0.0d0) moleculeCooling=coolants(N)%previousCooling
            ! ELSE
            !     coolants(N)%previousCooling=moleculeCooling
            ! END IF
            IF (moleculeCooling .gt. -1.0d-30 .and. abundances(coolantIndices(N)) .gt. 1.0d-20) lineCooling= lineCooling+moleculeCooling
        END DO
        !!write(*,*) lineCooling
    END FUNCTION lineCooling

    ! !-----------------------------------------------------------------------
    ! !  Atomic and ionic cooling rates
    ! !  from Neal et al. 1995 based on Cen (1992) via Grassi et al. (2014)
    ! !-----------------------------------------------------------------------
    REAL(dp) FUNCTION atomicCooling(gasT,gasDensity,hAbund,heAbund,electronAbund,hxAbund,hexAbund)
        REAL(dp), INTENT(IN) :: gasT,gasDensity,hAbund,heAbund,electronAbund,hxAbund,hexAbund
        REAL(dp) :: t5,invT,rootT,collTFactor !temp/10^5, 1/T and a weird factor from the table
        REAL(dp) :: hDens,elecDens,heDens,hxDens,hexDens,gauntFactor
        hDens=gasDensity*hAbund
        elecDens=gasDensity*electronAbund
        heDens=gasDensity*heAbund
        hxDens=gasDensity*hxAbund
        hexDens=gasDensity*hexAbund
        t5=1.0d-5*gasT
        invT=1.0/gasT
        rootT=SQRT(gasT)
        collTFactor=1.0/(1.0+SQRT(t5))

        !gauntFactor from Neal et al. 1995
        gauntFactor=1.1+(0.34*EXP(-((5.5-LOG10(gasT))**2.0)/3.0))
        !Neal et al. 1995 lists several fits to cooling each in ergs/cm3/s so we'll just sum them
        !see table 1 of that paper
        !These are just numerical fits so there's loads of magic numbers
        !I've shorted variable names to make it easier to write/read (tn is temperature/10^n)

        !collisional excitation and ionization
        atomicCooling=(7.5d-19*collTFactor*EXP(-118348.0*invT)*elecDens*hDens) &
            &+(5.54d-17*(gasT**-0.397)*collTFactor*EXP(-473638.0*invT)*elecDens*hexDens)&
            &+(1.27d-21*rootT*EXP(-157809.1*invT)*elecDens*hDens*collTFactor)&
            &+(9.38d-22*rootT*EXP(-285335.4*invT)*elecDens*heDens*collTFactor)&
            &+(4.95d-22*rootT*EXP(-631515.0*invT)*elecDens*hexDens*collTFactor)&
            !dielectric
            &+(1.24d-13*(gasT**-1.5)*EXP(-470000.0*invT)*(1.0+0.3*EXP(-94000.0*invT))*elecDens*hexDens)
        IF (gasT .gt. 1.0d5) THEN
        !recombination
            atomicCooling=atomicCooling&
            &+(8.7d-27*rootT*((1.0d-3*gasT)**-0.2)*elecDens*hxDens/(1.0+((0.1*t5)**0.7)))&
            &+(1.55d-26*(gasT**0.3647)*elecDens*hexDens)&
            !&+(3.48d-26*rootT*((0.001*gasT)**-0.2)*nelec*nhexI/(1+(0.1*t5)**0.7))&
            !Free-free emission
            &+(1.42d-27*rootT*nelec*(nhex+nhx)*gauntFactor)
        END IF

    END FUNCTION atomicCooling

    !!-----------------------------------------------------------------------
    !! Collisionally Induced Emission
    !! Hirano & Yoshida (2013) and Ripamonti & Abel 2004 via Grassi 2012
    !!-----------------------------------------------------------------------
    REAL(dp) FUNCTION collionallyInducedEmission(gasTemperature,gasDensity,h2Abund)
        REAL(dp), INTENT(IN) :: gasTemperature,gasDensity,h2Abund
        REAL(dp), PARAMETER :: aConsts(6)=(/-30.3314216559651,19.0004016698518,-17.1507937874082&
                                                            &,9.49499574218739,-2.54768404538229,0.265382965410969/)
        REAL(dp), PARAMETER :: bConsts(6)=(/-180.992524120965,168.471004362887,-67.499549702687,&
                                                            &13.5075841245848,-1.31983368963974,0.0500087685129987/)
        REAL(dp),PARAMETER :: c=3.0,d=21.2968837223113
        REAL(dp) :: tau,logt
        INTEGER :: i        

        logt=LOG10(gasTemperature)

        tau=(gasDensity*h2Abund/7.0d15)**2.8
        !if (tau.lt.0.2) THEN
            !avoid numerical problems, tau tends to 1 for low tau but fortran can't do it
            !Taylor series is fine for tau<0.2 and that's well above the area we get an issue
         !   tau=1.0-(0.5*tau)+((tau**2.0)/6.0)-((tau**3.0)/24.0)+((tau**4)/120.0)
        !ELSE
            tau=(1.0-dexp(-tau))/tau
        !END IF
        tau=min(1.0,tau)


        collionallyInducedEmission=0.0

        IF (gasTemperature .ge. 1.0d5) THEN
            collionallyInducedEmission=(c*logt)-d            
        ELSE IF (gasTemperature .ge. 891.0) THEN
            DO i=1,SIZE(bConsts)
                collionallyInducedEmission=collionallyInducedEmission+(bConsts(i)*(logt**(i-1)))
            END DO
        !technically fit below is ok down to 100 K but bad fit seems better than no cooling at 70 K
        ELSE IF (gasTemperature .ge. 100.0) THEN
            DO i=1,SIZE(aConsts)
                collionallyInducedEmission=collionallyInducedEmission+(aConsts(i)*(logt**(i-1)))
            END DO
        END IF

        if (gasTemperature .ge. 100.0) collionallyInducedEmission=(10.0**collionallyInducedEmission)*tau
        !collionallyInducedEmission=(10.0**collionallyInducedEmission)*tau
    END FUNCTION collionallyInducedEmission


    !!-----------------------------------------------------------------------
    !! Continuum Emission
    !! Hirano & Yoshida (2013) and Ripamonti & Abel 2004 via Grassi 2012
    !!-----------------------------------------------------------------------
    REAL(dp) FUNCTION continuumEmission(gasTemperature,gasDensity)
        REAL(dp), INTENT(IN) :: gasTemperature,gasDensity
        REAL(dp) :: massDensity,opacity,opticalDepth
        massDensity=min(0.5,gasDensity*MH*1.22) !assume mean molecular weight 1.22 and give max
        opacity=10.0**(1.000042*log(massDensity)+2.14989) !Lenzuni opacity fit

        opticalDepth=SQRT(3.14159*K_BOLTZ*gasTemperature/(massDensity*MH*1.22*GRAV_G))
        opticalDepth=opticalDepth*opacity*massDensity+1.0d-40 ! stop it going to zero

        continuumEmission=4.0*SB_CONST*(gasTemperature**4.0)*opacity*massDensity*min((opticalDepth**(-2.0)),1.0)
    END FUNCTION continuumEmission

    !!-----------------------------------------------------------------------
    !! Compton cooling
    !! Cen 1992 via Grassi 2012
    !! Cooling due to compton scattering of CMB photons.
    !! Shouldn't be important in near universe but include for completion
    !!-----------------------------------------------------------------------
    REAL(dp) FUNCTION comptonCooling(gasTemperature,gasDensity,elecAbund)
        REAL(dp), INTENT(IN) :: gasTemperature,gasDensity,elecAbund
        REAL(dp), PARAMETER :: cmbTemp=2.73
            comptonCooling=1.017d-37*(cmbTemp**4.0)*(gasTemperature-cmbTemp)*elecAbund*gasDensity
    END FUNCTION comptonCooling

    ! !-----------------------------------------------------------------------
    ! !  Grain + PAH photoelectric heating (with graphitic and silicate grains)

    ! !  Weingartner & Draine (2001, ApJS, 134, 263)
    ! !
    ! !  Includes photoelectric heating due to PAHs, VSGs and larger grains
    ! !  Assumes a gas-to-dust mass ratio of 100:1
    ! !-----------------------------------------------------------------------
    REAL(dp) FUNCTION photoelectricHeating(gasTemperature,gasDensity,habingField,electronAbund)
        REAL(dp), INTENT(IN) :: gasTemperature,gasDensity,habingField,electronAbund
        REAL(dp), PARAMETER :: PHI_PAH=0.4d0,ALPHA=0.944D0
        REAL(dp) :: beta,delta,epsilon,nElec,PAH_HEATING_RATE,PAH_COOLING_RATE
        ! REAL(dp), PARAMETER ::C0=5.72D+0,C1=3.45D-2,C2=7.08D-3
        ! REAL(dp), PARAMETER ::C3=1.98D-2, C4=4.95D-1,C5=6.92D-1
        ! REAL(dp), PARAMETER ::C6=5.20D-1

        ! !Weingartner & Draine2001
        ! photoelectricHeating=1.0D-26*(habingField*gasDensity)*(C0+C1*gasTemperature**C4) &
        !         & /(1.0D0+C2*(habingField*SQRT(gasTemperature)/(gasDensity*electronAbund))**C5  &
        !         & *(1.0D0+C3*(habingField*SQRT(gasTemperature)/(gasDensity*electronAbund))**C6))

        !Bakes & Tielens 1994 with updates from Wolfire 2008
        !  Adopt the PAH rate scaling factor of Wolfire et al. (2008, ApJ, 680, 384)
        !  Setting this factor to 1.0 gives the standard Bakes & Tielens expression

        nElec=electronAbund*gasDensity
        BETA=0.735D0/gasTemperature**0.068
        DELTA=habingField*SQRT(gasTemperature)/(nElec*PHI_PAH)
        EPSILON=4.87D-2/(1.0D0+4.0D-3*DELTA**0.73) + 3.65D-2*(gasTemperature/1.0D4)**0.7/(1.0D0+2.0D-4*DELTA)

        PAH_HEATING_RATE=1.30D-24*EPSILON*habingField*gasDensity
        PAH_COOLING_RATE=4.65D-30*gasTemperature**ALPHA*(DELTA**BETA)*nElec*PHI_PAH*gasDensity
        photoelectricHeating=PAH_HEATING_RATE-PAH_COOLING_RATE
    ! !  Assume the PE heating rate scales linearly with metallicity
    !    photoelectricHeating=photoelectricHeating*METALLICITY
    END FUNCTION photoelectricHeating

    ! !-----------------------------------------------------------------------
    ! !  H2 formation heating
    ! !
    ! !  Assume that 1.5 eV is liberated as heat during H2 formation
    ! !
    ! !  See: Hollenbach & Tielens (Review of Modern Physics, 1999, 71, 173)
    ! !  Use the H2 formation rate determined by subroutine H2_FORMATION_RATE
    ! !  and stored as REACTION_RATE(nRGR) (cm^3 s^-1)
    !! JH: I replaced REACTION_RATE(nRGR) with chemistry.f90's h2form=1.0d-17*dsqrt(Y(NEQ-1))

    ! !-----------------------------------------------------------------------
    REAL(dp) FUNCTION H2FormationHeating(gasTemperature,gasDensity,hAbund)
        REAL(dp), INTENT(IN) :: gasTemperature,gasDensity,hAbund
        H2FormationHeating=(1.5*eV)*1.0d-17*dsqrt(gasTemperature)*hAbund*gasDensity*gasDensity
    END FUNCTION H2FormationHeating


    !-----------------------------------------------------------------------
    !  H2 vibrational heating/cooling
    !
    !  Treat the vibrationally excited levels of H2 as a single pseudo level
    !  with effective rates of spontaneous emission, collisional excitation,
    !  FUV pumping and photodissociation that describe the behaviour of all
    !  the vibrational levels combined.
    !
    !  Use the treatment of Rollig et al. (2006, A&A, 451, 917)
    !-----------------------------------------------------------------------
    REAL(dp) FUNCTION H2VibrationalCooling(gasTemperature,gasDensity,h2Abund,h2dis)
        REAL(dp), INTENT(IN) :: gasTemperature,gasDensity,h2Abund,h2dis
        REAL(dp) :: photoDisRate,DELTA_E_10,A_COEFF_10,C_COEFF_10
        REAL(dp) ::DELTA_E_EFF,A_COEFF_EFF,R_PUMP_EFF,R_PHOTO_EFF
        DELTA_E_10=6587.0 ! Energy gap (K) between the v=1 and v=0 levels of H2
        A_COEFF_10=8.6D-7 ! Einstein A-coefficient (s^-1) for emission from the v=1 to v=0 level of H2
        C_COEFF_10=5.4D-13*SQRT(gasTemperature) ! Collisional rate coefficient (cm^3 s^-1) for v=0 to v=1
        photoDisRate=h2dis ! Photodissociation rate (s^-1) for the v=1 level of H2

        DELTA_E_EFF=23500.0 ! Characteristic vibrational level energy (K)
        A_COEFF_EFF=1.9D-6  ! Effective Einstein A-coefficient (s^-1)

        R_PUMP_EFF=11.2*h2dis ! Effective vibrational pumping rate (s^-1)
        R_PHOTO_EFF=18.0*h2dis ! Effective photodissociation rate (s^-1)

        H2VibrationalCooling=K_BOLTZ*DELTA_E_10*C_COEFF_10*gasDensity*EXP(-DELTA_E_10/gasTemperature)*h2Abund*gasDensity&
                           & *(A_COEFF_10+h2dis)/(C_COEFF_10*gasDensity+A_COEFF_10+h2dis)

        !Some heating from H2 vibrational interactions so subtract from cooling rate
        H2VibrationalCooling=H2VibrationalCooling-(h2Abund*gasDensity*(R_PUMP_EFF*K_BOLTZ*DELTA_E_EFF) &
                           & /(1.0D0+(A_COEFF_10+R_PHOTO_EFF)/(C_COEFF_10*gasDensity)))
    END FUNCTION H2VibrationalCooling

    ! !-----------------------------------------------------------------------
    ! !  H2 photodissociation heating
    ! !
    ! !  On average, 0.4 eV of kinetic energy per photodissociated molecule
    ! !
    ! !  Use the H2 photodissociation rate determined by the subroutine
    ! !  CALCULATE_REACTION_RATES and stored as REACTION_RATE(nRH2) (s^-1)
    ! !  JH: again, grabbed h2dis from chemistry.f90 for consistency.
    ! !-----------------------------------------------------------------------
    REAL(dp) FUNCTION H2PhotodisHeating(gasDensity,h2Abund,h2dis)
        REAL(dp), INTENT(IN) :: gasDensity,h2Abund,h2dis
        H2PhotodisHeating=(0.4*eV)*h2dis*gasDensity*h2Abund
    END FUNCTION H2PhotodisHeating

    ! !-----------------------------------------------------------------------
    ! !  Cosmic-ray ionization heating
    ! !
    ! !  20.0 eV of kinetic energy deposited per H2 ionization,
    ! !  based on the estimate of Goldsmith (2001, ApJ, 557,736)
    ! !
    ! !  See also:
    ! !  Clavel et al. (1978, A&A, 65, 435)
    ! !  Tielens & Hollenbach (1985, ApJ, 291, 722)
    ! !  Shull & van Steenberg (1985, ApJ, 298, 268)
    ! !  Kamp & van Zadelhoff (2001)
    ! !-----------------------------------------------------------------------
    REAL(dp) FUNCTION cosmicRayHeating(zeta,gasDensity,h2Abund)
        REAL(dp), INTENT(IN) :: zeta,gasDensity,h2Abund
        cosmicRayHeating=(20.0*eV)*(1.3D-17*zeta)*h2Abund*gasDensity
    END FUNCTION cosmicRayHeating


    ! !-----------------------------------------------------------------------
    ! !  H2 FUV pumping heating
    ! !
    ! !  On average, 2.2 eV released per vibrationally excited H2* molecule
    ! !
    ! !  Use the treatment of Hollenbach & McKee (1979, ApJ)

    ! !  Use the H2 photodissociation rate determined by the subroutine
    ! !  CALCULATE_REACTION_RATES and stored as REACTION_RATE(nRH2) (s^-1)
    ! !
    ! !  Use the H2 critical density expression from Hollenbach & McKee (1979)
    ! !  NOTE: The equation for the collisional de-excitation rate coefficient
    ! !  for the v=2-1 and v=1-0 transitions by collisions with H2 was wrongly
    ! !  stated in the Hollenbach & McKee (1979) paper (equation 6.29), but is
    ! !  corrected in Hollenbach & McKee (1989, equation 2.8) and used below.

    ! ! JH: h2dis instead of Rate(nRH2) again
    ! !-----------------------------------------------------------------------
    REAL(dp) FUNCTION h2FUVPumpHeating(hAbund,h2Abund,gasTemperature,gasDensity,h2dis)
        REAL(dp), INTENT(IN) :: hAbund,h2Abund,gasTemperature,gasDensity,h2dis
        REAL(dp) :: NCRIT_H2
        NCRIT_H2=1.0D6/SQRT(gasTemperature)/(1.6D0*hAbund*EXP(-((400.0D0/gasTemperature)**2)) &
                                      & + 1.4D0*h2Abund*EXP(-(18100.0D0/(gasTemperature+1200.0D0))))

        h2FUVPumpHeating=(2.2*eV)*9.0D0*h2dis*gasDensity*h2Abund/(1.0D0+NCRIT_H2/gasDensity)

        ! !  If vibrationally excited H2 (H2*) is included in the chemical network,
        ! !  then use the treatment of Tielens & Hollenbach (1985, ApJ, 291, 722)
        !    IF(nH2v.NE.0) THEN
        !       H2_FUV_PUMPING_HEATING_RATE=(DENSITY(nH)*1.0D-12*SQRT(gasTemperature)*EXP(-1000.0D0/gasTemperature) &
        !                                & +DENSITY(nH2)*1.4D-12*SQRT(gasTemperature)*EXP(-18100.0D0/(gasTemperature+1200.0D0))) &
        !                                  *(2.6*eV)*DENSITY(nH2v)
        !    END IF

    END FUNCTION h2FUVPumpHeating   



! !-----------------------------------------------------------------------
! !  Carbon photoionization heating
! !
! !  On average, 1 eV of kinetic energy released per carbon ionization
! !  Use the carbon photoionization rate determined by the subroutine
! !  CALCULATE_REACTION_RATES and stored as REACTION_RATE(nRCI) (s^-1)
! !-----------------------------------------------------------------------

FUNCTION CarbonIonizationHeating(cIonizationRate,carbonAbund,gasDensity)
    real(dp), intent(in) :: cIonizationRate,carbonAbund,gasDensity
    real(dp) :: CarbonIonizationHeating
    CarbonIonizationHeating=(1.0*eV)*cIonizationRate*carbonAbund*gasDensity
END FUNCTION

! !-----------------------------------------------------------------------
! !  Exothermic chemical reaction heating
! !
! !  See:
! !  Clavel et al. (1978, A&A,  65, 435)
! !  Meijerink & Spaans (2005, A&A, 436, 397)
! !  Glassgold & Langer (1973, ApJ, 179, L147)
! !
! !  Recombination reactions:
! !     H2+ (10.9 eV); H3+ (9.23+4.76 eV); H3O+ (1.16+5.63+6.27 eV); HCO+ (7.51 eV)
! !
! !  Ion-neutral reactions:
! !     H2+ + H (0.94 eV); He+ + H2 (6.51 eV); He+ + CO (2.22 eV)
! !
! !  For each reaction, the heating rate is given by: n(1) * n(2) * K * E
! !  where n(1) and n(2) are the number densities, K the rate coefficient
! !  (cm^3 s^-1), and E the energy released (erg).
! !-----------------------------------------------------------------------

Function chemicalHeating(gasDensity,exoReactants1,exoReactants2,exoRates,exothermicities)
REAL(dp), INTENT(IN) :: gasDensity,exoReactants1(:),exoReactants2(:),exoRates(:),exothermicities(:)
REAL(dp) :: chemicalHeating

    chemicalHeating=SUM(exoReactants1*exoReactants2*exoRates*exothermicities)
    chemicalHeating=chemicalHeating*gasDensity*gasDensity*EV !each abundance should be a number dnesity to multiply through
  END FUNCTION chemicalHeating
! !-----------------------------------------------------------------------
! !  Gas-grain collisional heating
! !
! !  Use the treatment of Burke & Hollenbach (1983, ApJ, 265, 223)
! !
! !  Other relevant references:
! !  Hollenbach & McKee (1979, ApJS, 41, 555)
! !  Tielens & Hollenbach (1985, ApJ, 291, 722)
! !  Goldsmith (2001, ApJ, 557, 736)
! !  Young et al. (2004, ApJ, 614, 252)
! !
! !  This process is insignificant for the energy balance of the dust,
! !  but can influence the gas temperature. If the dust temperature is
! !  lower than the gas temperature, this becomes a cooling mechanism.
! !-----------------------------------------------------------------------

FUNCTION gasGrainCollisions(gasTemperature,gasDensity,dustAbundance)
    real(dp), intent(in) :: gasTemperature,gasDensity,dustAbundance
    REAL(dp) :: gasGrainCollisions
    REAL(dp) :: nGrain,accommodation,dustTemp=20.0,C_GRAIN
    nGrain=dustAbundance*gasDensity
    C_GRAIN=PI*1.d-5**2

    !!$!  Accommodation fitting formula of Groenewegen (1994, A&A, 290, 531)
    !!$   ACCOMMODATION=0.35D0*EXP(-SQRT((DUST_TEMPERATURE+gasTemperature)/5.0D2))+0.1D0

    !  Accommodation coefficient of Burke & Hollenbach (1983, ApJ, 265, 223)
    accommodation=0.37D0*(1.0D0-0.8D0*EXP(-75.0D0/gasTemperature))

    gasGrainCollisions=nGrain*C_GRAIN*gasDensity*SQRT(8*K_BOLTZ*gasTemperature/(PI*MH)) &
                       & *accommodation*(2*K_BOLTZ*dustTemp-2*K_BOLTZ*gasTemperature)
END FUNCTION gasGrainCollisions



! !-----------------------------------------------------------------------
! !  Coulomb heating
! !
! !  Use the treatment of Meijerink & Spaans (2005, A&A, 436, 397)
! !
! !  Other relevant references:
! !  Dalgarno et al. (1999, ApJS, 125, 237)
! !
! !  This is an X-ray heating mechanism. When X-rays are absorbed, fast
! !  electrons are produced. These fast electrons lose part of their
! !  energy through Coulomb interactions with thermal electrons.
! !-----------------------------------------------------------------------

!    R=DENSITY(nH2)/DENSITY(nH) ! n(H2)/n(H) ratio
!    X_PRIME=1.83D0*ABUNDANCE(nelect)/(1.0D0+0.83D0*ABUNDANCE(nelect)) ! Correction to the electron abundance for a pure H2-He mixture

!    ETA_H2_He=1.0D0+(0.055D0-1.0D0)/(1.0D0+2.17D0*X_PRIME**0.366) ! Heating efficiency for a pure H2-He mixture
!    ETA_H_He =1.0D0+(0.117D0-1.0D0)/(1.0D0+7.95D0*ABUNDANCE(nelect)**0.678) ! Heating efficiency for a pure H-He mixture
!    ETA=(10.0D0*R*ETA_H2_He+ETA_H_He)/(10.0D0*R+1.0D0) ! Total heating efficiency for mixed atomic and molecular gas
!    H_X=PARTICLE%XRAY_ENERGY_DEPOSITION_RATE ! X-ray energy deposition rate per hydrogen nucleus (erg s^-1)

!    COULOMB_HEATING_RATE=ETA*GAS_DENSITY*H_X




! !-----------------------------------------------------------------------
! !  Supersonic turbulent decay heating
! !
! !  Most relevant for the inner parsecs of galaxies (Black)
! !  Black, in Interstellar Processes, 1987, p731
! !  See also: Rodriguez-Fernandez et al., 2001, A&A, 365, 174
! !
! !  V_TURB = turbulent velocity (km/s); Galactic center ~ 15 km/s
! !  L_TURB = turbulent scale length (pc); typically 5 pc
! !-----------------------------------------------------------------------

!    L_TURB=5.0D0
!    TURBULENT_HEATING_RATE=3.5D-28*((V_TURB/1.0D5)**3)*(1.0D0/L_TURB)*GAS_DENSITY




END MODULE heating

!Abandoned heating and cooling mechanisms

! ! !-----------------------------------------------------------------------
! ! !  Molecular Hydrogen cooling 
! ! !  from Galli & Palla (1998) via Grassi et al. (2014)
! ! !  Acceptable up to 10^5 K can easily use  Glover & Abel (2008) fits instead
! ! !  but they're more complex.
! ! ! JH: Dropped in favour of UCLPDR treatment of H2 rotational cooling + h2 Vibrational
! ! !-----------------------------------------------------------------------
! REAL(dp) FUNCTION h2Cooling(gasDensity,gasTemperature,hAbund,h2Abund)
!     REAL(dp), INTENT(IN) :: gasTemperature,gasDensity,hAbund,h2Abund
!     REAL(dp) :: highDensLimit,lowDensLimit,T1000,logT
!     ! Temperature in kilokelvin
!     T1000=gasTemperature*0.001

!     !high density limit is same in all models: sum of vibrational and rotational cooling
!     highDensLimit=((9.5d-22*T1000**3.76)/(1+0.12*(T1000**2.1)))*dexp(-(0.13/T1000)**3.0)
!     highDensLimit=highDensLimit+(3.0d-24*dexp(-0.51/T1000))
    
!     highDensLimit=highDensLimit+(6.7d-19*dexp(-5.86/T1000))+(1.6d-18*dexp(-11.7/T1000))

!     !I'm using Galli & Palli limit here which is ok up to 10^5.
!     !I think Glover and Abel is more accurate but is many fits so hard to imlpement
!     logT=log10(gasTemperature)
!     lowDensLimit=-103.0+(97.59*logT)-(48.05*logT*logT)+(10.8*logT**3.0)-(0.9032*logT**4.0)
!     lowDensLimit=10.0**(lowDensLimit*gasDensity*hAbund)

!     !Combine together for final cooling rate
!     h2Cooling=gasDensity*h2Abund*highDensLimit
!     h2Cooling=h2Cooling/(1.0+(highDensLimit/lowDensLimit))
! END FUNCTION h2Cooling




!=========================================
!     JH: Alternative photoelectric heating rates. The correct one depends on dust distribution
!     The one coded as photoelectricHeating() is the current one in ucl in UCLPDR
!-----------------------------------------------------------------------
!  Grain photoelectric heating (large grains only; r ~ 100 Å)
!
!  Use the treatment of 

!  which in turn follows de Jong (1977, 1980)
!
!  The charge of a dust grain can be found by equating the rate of
!  photo-ejection of electrons from the dust grain to the rate of
!  recombination of electrons with the dust grain (Spitzer)
!
!  The various parameter values are taken from Table 2 of the paper
!-----------------------------------------------------------------------
!     REAL(dp) FUNCTION photoelectricHeating(gasTemperature,gasDensity,habingField,electronAbund)
!         REAL(dp), INTENT(IN) :: gasTemperature,gasDensity,habingField,electronAbund
!         REAL(dp), PARAMETER :: DELTA_D=1.0D0
!         REAL(dp), PARAMETER :: DELTA_UV=1.8D0
!         REAL(dp), PARAMETER :: Y=0.1D0
!         REAL(dp), PARAMETER :: HNU_D=6.0D0
!         REAL(dp), PARAMETER :: HNU_H=13.6D0
!         REAL(dp) :: delta,gamma,XK,XD,X,XX
!         INTEGER :: ITERATION


!         XK=K_BOLTZ*gasTemperature/(HNU_H*eV)
!         XD=HNU_D/HNU_H
!         gamma=2.9D-4*Y*DSQRT(gasTemperature)*habingField/(gasDensity*electronAbund)
!         delta=XK-XD+gamma

!         !  Iterate to determine X by finding the zero of the function F
!         X=0.5D0
!         DO ITERATION=1,100
!           XX=X-(grainChargeFunc(X,DELTA,GAMMA)/deltaGrainChargeFunc(X,DELTA))
!           IF(ABS(XX-X).LT.1.0D-2) EXIT
!           X=XX
!         END DO
!         X=XX

!         IF(ITERATION.GE.100) THEN
!           WRITE(10,*)'WARNING! Grain parameter X not found in PE heating'
!           WRITE(10,*)'Using final value from interation loop: X =',X
!         END IF

!         photoelectricHeating=2.7D-25*DELTA_UV*DELTA_D*gasDensity*Y*habingField &
!                                      & *(((1.0D0-X)**2)/X + XK*((X**2)-1.0D0)/(X**2))

!         !  Assume the PE heating rate scales linearly with metallicity
!         !TH85_PHOTOELECTRIC_HEATING_RATE=TH85_PHOTOELECTRIC_HEATING_RATE*METALLICITY
!     END FUNCTION  photoelectricHeating

! !=======================================================================
! !  X is the grain charge parameter and is the solution to F(X)=0
! !-----------------------------------------------------------------------
!    FUNCTION grainChargeFunc(X,DELTA,GAMMA)
!       IMPLICIT NONE
!       REAL(dp) :: grainChargeFunc
!       REAL(dp), INTENT(IN) :: X,DELTA,GAMMA
!       grainChargeFunc=(X**3)+DELTA*(X**2)-GAMMA
!    END FUNCTION grainChargeFunc

! !=======================================================================
! !  FF(X) is the derivative of F(X) with respect to X
! !-----------------------------------------------------------------------
!    FUNCTION deltaGrainChargeFunc(X,DELTA)
!       IMPLICIT NONE
!       REAL(dp) :: deltaGrainChargeFunc
!       REAL(dp), INTENT(IN) :: X,DELTA
!       deltaGrainChargeFunc=3*(X**2)+DELTA*(2*X)
!    END FUNCTION deltaGrainChargeFunc

! !-----------------------------------------------------------------------
! !  Grain + PAH photoelectric heating (MRN size distribution; r = 3-100 Å)
! !
! !  Use the treatment of Bakes & Tielens (1994, ApJ, 427, 822) with the
! !  modifications suggested by Wolfire et al. (2003, ApJ, 587, 278) to
! !  account for the revised PAH abundance estimate from Spitzer data.
! !
! !  See also:
! !  Wolfire et al. (1995, ApJ, 443, 152)
! !  Le Page, Snow & Bierbaum (2001, ApJS, 132, 233)
! !-----------------------------------------------------------------------

! !  Adopt the PAH rate scaling factor of Wolfire et al. (2008, ApJ, 680, 384)
! !  Setting this factor to 1.0 gives the standard Bakes & Tielens expression
!    PHI_PAH=0.4D0

!    ALPHA=0.944D0
!    BETA=0.735D0/gasTemperature**0.068
!    DELTA=HABING_FIELD*SQRT(gasTemperature)/(DENSITY(nelect)*PHI_PAH)
!    EPSILON=4.87D-2/(1.0D0+4.0D-3*DELTA**0.73) + 3.65D-2*(gasTemperature/1.0D4)**0.7/(1.0D0+2.0D-4*DELTA)

!    PAH_HEATING_RATE=1.30D-24*EPSILON*HABING_FIELD*GAS_DENSITY
!    PAH_COOLING_RATE=4.65D-30*gasTemperature**ALPHA*(DELTA**BETA)*DENSITY(nelect)*PHI_PAH*GAS_DENSITY

!    BT94_PHOTOELECTRIC_HEATING_RATE=PAH_HEATING_RATE - PAH_COOLING_RATE

! !  Assume the PE heating rate scales linearly with metallicity
!    BT94_PHOTOELECTRIC_HEATING_RATE=BT94_PHOTOELECTRIC_HEATING_RATE*METALLICITY

    ! FUNCTION getEquilibriumTemp(gasTemperature,gasDensity,habingField,abundances,h2dis,zeta,cIonRate,dustAbundance&
    !     &,exoReactants1,exoReactants2,exoRates,exothermicities)
    !     REAL(dp), INTENT(in) :: gasTemperature,gasDensity,habingField,h2dis,zeta,cIonRate,dustAbundance
    !     REAL(dp) :: getEquilibriumTemp
    !     REAL(dp), INTENT(in) :: abundances(:),exoReactants1(:),exoReactants2(:),exoRates(:),exothermicities(:)
    !     REAL(dp) :: adiabaticIdx,heating,cooling,fMax=5.0d-6,fRatio,maxTemp=1.0d6,minTemp=5.0d0
    !     INTEGER :: maxTempIter=100,ti

    !     maxTemp=1.0d6
    !     minTemp=5.0d0
    !     !initial guess is current temperature
    !     getEquilibriumTemp=gasTemperature

    !     !find heating and cooling rate for current temp
    !     heating=getHeatingRate(gasTemperature,gasDensity,habingField,abundances,h2dis,zeta,cIonRate,dustAbundance&
    !         &,exoReactants1,exoReactants2,exoRates,exothermicities)
    !     cooling=getCoolingRate(gasTemperature,gasDensity,abundances,h2dis)

    !      !then get measure of their differences
    !     fRatio=2.0D0*abs(heating-cooling)/abs(heating+cooling)
    !     !write(*,*) ti,getEquilibriumTemp,heating,Cooling,fRatio

    !     !whilst that difference is large, change temp to reduce it
    !     IF (fRatio .gt. fMax) THEN
    !         DO ti=1,maxTempIter
    !             !if heating>cooling gas is too cold so minTemp is at least current temp
    !             IF (heating .gt. cooling) minTemp=getEquilibriumTemp
    !             !equally cooling implies temp can't be larger than current.
    !             IF (cooling .gt. heating) maxTemp=getEquilibriumTemp
    !             !take average of min and max to get new trial
    !             getEquilibriumTemp=(minTemp+maxTemp)*0.5d0

    !             !recalculate heatnig,cooling and ratio
    !             heating=getHeatingRate(getEquilibriumTemp,gasDensity,habingField,abundances,h2dis,zeta,cIonRate,dustAbundance&
    !                 &,exoReactants1,exoReactants2,exoRates,exothermicities)
    !             cooling=getCoolingRate(getEquilibriumTemp,gasDensity,abundances,h2dis)
    !             fRatio=2.0D0*abs(heating-cooling)/abs(heating+cooling)
    !             IF (fRatio .lt. fMax) exit
    !             IF (abs(minTemp-maxTemp).lt.0.5) exit
    !             !write(*,*) ti,getEquilibriumTemp,heating,Cooling,minTemp,maxTemp
    !         END DO
    !     ELSE
    !         getEquilibriumTemp=gasTemperature
    !     END IF
    ! END FUNCTION getEquilibriumTemp