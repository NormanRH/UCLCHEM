!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Module that provides heating and cooling rates		  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE heating
IMPLICIT NONE
CONTAINS
    DOUBLE PRECISION FUNCTION getHeatingRate(currentTemp)
        DOUBLE PRECISION, INTENT(in) :: currentTemp
        DOUBLE PRECISION :: heatingRate
        IF (currentTemp .lt. 100) THEN
            getHeatingRate=5.0*3.16d-8
        ELSE
            getHeatingRate=0.0
        END IF
    END FUNCTION getHeatingRate

    DOUBLE PRECISION FUNCTION getCoolingRate()
        DOUBLE PRECISION :: coolingRate
        getCoolingRate=0.0 
    END FUNCTION getCoolingRate
END MODULE heating

