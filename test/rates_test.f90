PROGRAM RATES_TEST    
double precision, parameter ::  minTemps (18)=(/2.0000e+03&
    &,1.0000e+01,3.0000e+02,1.0000e+01,7.9500e+03,2.1140e+04,1.0000e+01&
    &,3.0000e+02,1.0000e+01,2.0000e+03,1.0000e+01,1.0000e+03,1.0000e+01&
    &,1.0000e+04,2.0000e+02,4.0000e+03,1.0000e+01,2.0000e+03/)
double precision, parameter ::  maxTemps (18)=(/1.0000e+04&
    &,3.0000e+02,1.3900e+04,7.9500e+03,2.1140e+04,4.1000e+04,3.0000e+02&
    &,1.3900e+04,2.0000e+03,6.0000e+03,1.0000e+03,4.1000e+04,1.0000e+04&
    &,4.1000e+04,4.0000e+03,3.2000e+04,2.0000e+03,1.0000e+04/)

    double precision :: rates(18),gasTemp
    integer :: i=0
    rates=52.0

    gasTemp=10.0
    !this multiplies rate by 0 or 1 depending on whether gastemp>mintemp of a reaction
    rates=rates*min(real(floor(gasTemp/minTemps)),1.0)
    !and this multiplies by 0,1 if gastemp>max temp
    rates=rates*min(real(floor(maxTemps/gasTemp)),1.0)

    DO i=1,18
        write(*,*) minTemps(i),maxTemps(i),rates(i)
    END DO

        rates=52.0

    gasTemp=1000.0
    write(*,*)"*****"
    write(*,*)gasTemp
    write(*,*)"*****"

    !this multiplies rate by 0 or 1 depending on whether gastemp>mintemp of a reaction
    rates=rates*min(real(floor(gasTemp/minTemps)),1.0)
    !and this multiplies by 0,1 if gastemp>max temp
    rates=rates*min(real(floor(maxTemps/gasTemp)),1.0)

    DO i=1,18
        write(*,*) minTemps(i),maxTemps(i),rates(i)
    END DO

END PROGRAM RATES_TEST