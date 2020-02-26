!WORK IN PROGRESS
!THIS IS A SERIES OF SUBROUTINES THAT CAN BE COMPILED WITH F2PY TO PRODUCE A PYTHON MODULE
!EACH SUBROUTINE BECOMES A PYTHON FUNCTION

SUBROUTINE temperatureTest(initDens,intemp,uvIn,finTime,outFile)
    USE physics
    USE chemistry
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: initDens,finTime,intemp,uvIn
    CHARACTER(LEN=*), INTENT(IN) :: outFile
    !f2py intent(in) initDens,finTime,intemp,outFile,uvIn
    CHARACTER (LEN=100):: abundFile,outputFile,columnFile

    INCLUDE 'defaultparameters.f90'
    close(10)
    close(11)
    close(7)


    open(11,file=outFile,status='unknown') 
    columnFlag=.True. !we've opened a column file not a full output stream


    initialDens=initDens
    initialTemp=intemp
    finalTime=finTime
    radfield=uvIn

    write(*,*) initialDens,initialTemp, radfield,finalTime

    heatingFlag=.True.


    dstep=1
    currentTime=0.0
    timeInYears=0.0

    CALL initializePhysics
    CALL initializeChemistry
    
    !loop until the end condition of the model is reached 
    DO WHILE ((switch .eq. 1 .and. density(1) < finalDens) .or. (switch .eq. 0 .and. timeInYears < finalTime))
        !store current time as starting point for each depth step
        currentTimeold=currentTime

        !Each physics module has a subroutine to set the target time from the current time
        CALL updateTargetTime

        !loop over parcels, counting from centre out to edge of cloud
        DO dstep=1,points
            !update chemistry from currentTime to targetTime
            CALL updateChemistry
            currentTime=targetTime
            !get time in years for output, currentTime is now equal to targetTime
            timeInYears= currentTime/SECONDS_PER_YEAR

            !Update physics so it's correct for new currentTime and start of next time step
            CALL updatePhysics
            !Sublimation checks if Sublimation should happen this time step and does it
            CALL sublimation(abund)
            !write this depth step now time, chemistry and physics are consistent
            CALL output

            !reset time for next depth point
            if (points .gt. 1)currentTime=currentTimeold
        END DO
    END DO 
END SUBROUTINE temperatureTest