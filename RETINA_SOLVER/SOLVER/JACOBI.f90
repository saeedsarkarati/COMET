MODULE JACOBI
    
    USE CRS
    IMPLICIT NONE
    
    TYPE :: J_DATA
        DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: COPY
    END TYPE
    
CONTAINS

    FUNCTION JACOBI_SOLVE_(SIZE, LENGTH, STARTS, INDICES, DIAG_INDICES, A, B, X, MAX_ITER, COPY, W) RESULT(ANS)    
        INTEGER, INTENT(IN) :: SIZE, LENGTH, MAX_ITER
        INTEGER, DIMENSION(0:SIZE), INTENT(IN) :: STARTS
        INTEGER, DIMENSION(0:LENGTH - 1), INTENT(IN) :: INDICES
        INTEGER, DIMENSION(0:SIZE - 1), INTENT(IN) :: DIAG_INDICES
        DOUBLE PRECISION, DIMENSION(0:LENGTH - 1), INTENT(IN) :: A
        DOUBLE PRECISION, DIMENSION(0:SIZE - 1), INTENT(IN) :: B
        DOUBLE PRECISION, DIMENSION(0:SIZE - 1), INTENT(INOUT) :: X
        DOUBLE PRECISION, DIMENSION(0:SIZE - 1), INTENT(INOUT) :: COPY
        DOUBLE PRECISION, INTENT(IN) :: W
        INTEGER :: ANS
        INTEGER :: I, J_INDEX, J, DIAG
        DOUBLE PRECISION :: SUM
        DO ANS = 0, MAX_ITER - 1      
            CALL CRS_COPY_(SIZE, X, COPY)
            !$OMP PARALLEL
                !$OMP DO SCHEDULE(STATIC) PRIVATE(I, DIAG, SUM, J_INDEX, J)
            DO I = 0, SIZE - 1
                DIAG = DIAG_INDICES(I)
                SUM = B(I)
                DO J_INDEX = STARTS(I), DIAG - 1
                    J = INDICES(J_INDEX)
                    SUM = SUM - A(J_INDEX) * COPY(J)
                END DO
                DO J_INDEX = DIAG + 1, STARTS(I + 1) - 1                    
                    J = INDICES(J_INDEX)
                    SUM = SUM - A(J_INDEX) * COPY(J)
                END DO
                X(I) = SUM / A(DIAG) * W + (1.D0 - W) * COPY(I)
            END DO
                !$OMP END DO
            !$OMP END PARALLEL
            IF (ANS .LT. (MAX_ITER - 1)) THEN
                IF (CRS_CONVERGED_(SIZE, X, COPY, 1.D-7, 1.D-12)) THEN
                    EXIT
                END IF
            END IF            
        END DO
        ANS = ANS + 1
    END FUNCTION

    
END MODULE
    