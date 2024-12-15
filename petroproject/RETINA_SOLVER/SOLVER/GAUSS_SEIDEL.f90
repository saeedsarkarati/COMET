MODULE GAUSS_SEIDEL
    
    USE CRS
    IMPLICIT NONE
    
    TYPE :: GS_DATA
        DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: COPY
    END TYPE
    
CONTAINS

    FUNCTION GS_SOLVE_PARALLEL_FORWARD_(SIZE, W_SIZE, LEVEL_STARTS, LENGTH, THREADS, STARTS, INDICES, DIAG_INDICES, INDEPENDENT_ROW_STARTS, INDEPENDENT_ROWS, A, B, Y, X, MAX_ITER, PARAMS) RESULT(ANS)
        INTEGER, INTENT(IN) :: SIZE, W_SIZE, LENGTH, LEVEL_STARTS, THREADS, MAX_ITER
        INTEGER, DIMENSION(0:SIZE), INTENT(IN) :: STARTS
        INTEGER, DIMENSION(0:LENGTH - 1), INTENT(IN) :: INDICES
        INTEGER, DIMENSION(0:SIZE - 1), INTENT(IN) :: DIAG_INDICES
        INTEGER, DIMENSION(0:LEVEL_STARTS - 1), INTENT(IN) :: INDEPENDENT_ROW_STARTS
        INTEGER, DIMENSION(0:SIZE - W_SIZE), INTENT(IN) :: INDEPENDENT_ROWS
        DOUBLE PRECISION, DIMENSION(0:LENGTH - 1), INTENT(IN) :: A
        DOUBLE PRECISION, DIMENSION(0:SIZE - 1), INTENT(IN) :: B
        DOUBLE PRECISION, DIMENSION(0:SIZE - 1), INTENT(INOUT) :: Y
        DOUBLE PRECISION, DIMENSION(0:SIZE - 1), INTENT(OUT) :: X
        TYPE(GS_DATA), INTENT(INOUT) :: PARAMS
        INTEGER :: ANS
        INTEGER :: LEVEL, LEVEL_SIZE, LEVEL_INDEX
        INTEGER :: I, J_INDEX, J
        INTEGER :: S, E, INC
        DOUBLE PRECISION :: SUM
        
        DO ANS = 0, MAX_ITER - 1      
            CALL CRS_COPY_(SIZE, X, PARAMS%COPY)
            !$OMP PARALLEL
                !$OMP DO SCHEDULE(STATIC) PRIVATE(I, SUM, J_INDEX, J)
            DO I = 0, SIZE - 1
                SUM = B(I)
                DO J_INDEX = DIAG_INDICES(I) + 1, STARTS(I + 1) - 1
                    J = INDICES(J_INDEX)
                    SUM = SUM - A(J_INDEX) * X(J)
                END DO
                Y(I) = SUM
            END DO
                !$OMP END DO
            !$OMP END PARALLEL            
            DO LEVEL = 0, LEVEL_STARTS - 2
                LEVEL_SIZE = INDEPENDENT_ROW_STARTS(LEVEL + 1) - INDEPENDENT_ROW_STARTS(LEVEL)
                IF (LEVEL_SIZE .LE. 10) THEN
                    CALL OMP_SET_NUM_THREADS(1)
                ELSE
                    CALL OMP_SET_NUM_THREADS(THREADS)
                END IF
                !$OMP PARALLEL
                    !$OMP DO SCHEDULE(STATIC) PRIVATE(LEVEL_INDEX, I, SUM, J_INDEX, J)
                DO LEVEL_INDEX = INDEPENDENT_ROW_STARTS(LEVEL), INDEPENDENT_ROW_STARTS(LEVEL + 1) - 1
                    I = INDEPENDENT_ROWS(LEVEL_INDEX)
                    SUM = Y(I)
                    DO J_INDEX = STARTS(I), DIAG_INDICES(I) - 1
                        J = INDICES(J_INDEX)
                        SUM = SUM - A(J_INDEX) * X(J)
                    END DO
                    X(I) = SUM / A(J_INDEX)
                END DO
                    !$OMP END DO
                !$OMP END PARALLEL
            END DO
            
            DO I = SIZE - W_SIZE, SIZE - 1
                SUM = Y(I)
                DO J_INDEX = STARTS(I), DIAG_INDICES(I) - 1
                    J = INDICES(J_INDEX)
                    SUM = SUM - A(J_INDEX) * X(J)
                END DO
                X(I) = SUM / A(J_INDEX)
            END DO
            IF (ANS .LT. (MAX_ITER - 1)) THEN
                IF (CRS_CONVERGED_(SIZE, X, PARAMS%COPY, 1.D-7, 1.D-12)) THEN
                    EXIT
                END IF
            END IF            
        END DO
        CALL OMP_SET_NUM_THREADS(THREADS)
        ANS = ANS + 1        
    END FUNCTION
    
    FUNCTION GS_SOLVE_PARALLEL_BACKWARD1_(SIZE, W_SIZE, LEVEL_STARTS, LENGTH, THREADS, STARTS, INDICES, DIAG_INDICES, INDEPENDENT_ROW_STARTS, INDEPENDENT_ROWS, A, B, Y, X, MAX_ITER, PARAMS) RESULT(ANS)
        INTEGER, INTENT(IN) :: SIZE, W_SIZE, LENGTH, LEVEL_STARTS, THREADS, MAX_ITER
        INTEGER, DIMENSION(0:SIZE), INTENT(IN) :: STARTS
        INTEGER, DIMENSION(0:LENGTH - 1), INTENT(IN) :: INDICES
        INTEGER, DIMENSION(0:SIZE - 1), INTENT(IN) :: DIAG_INDICES
        INTEGER, DIMENSION(0:LEVEL_STARTS - 1), INTENT(IN) :: INDEPENDENT_ROW_STARTS
        INTEGER, DIMENSION(0:SIZE - W_SIZE - 1), INTENT(IN) :: INDEPENDENT_ROWS
        DOUBLE PRECISION, DIMENSION(0:LENGTH - 1), INTENT(IN) :: A
        DOUBLE PRECISION, DIMENSION(0:SIZE - 1), INTENT(IN) :: B
        DOUBLE PRECISION, DIMENSION(0:SIZE - 1), INTENT(INOUT) :: Y
        DOUBLE PRECISION, DIMENSION(0:SIZE - 1), INTENT(OUT) :: X
        TYPE(GS_DATA), INTENT(INOUT) :: PARAMS
        INTEGER :: ANS
        INTEGER :: LEVEL, LEVEL_SIZE, LEVEL_INDEX
        INTEGER :: I, J_INDEX, J
        INTEGER :: S, E, INC
        DOUBLE PRECISION :: SUM
        
        DO ANS = 0, MAX_ITER - 1      
            CALL CRS_COPY_(SIZE, X, PARAMS%COPY)
            !$OMP PARALLEL
                !$OMP DO SCHEDULE(STATIC) PRIVATE(I, SUM, J_INDEX, J)
            DO I = 0, SIZE - 1
                SUM = B(I)
                DO J_INDEX = STARTS(I), DIAG_INDICES(I) - 1
                    J = INDICES(J_INDEX)
                    SUM = SUM - A(J_INDEX) * X(J)
                END DO
                Y(I) = SUM
            END DO
                !$OMP END DO
            !$OMP END PARALLEL            
            
            DO I = SIZE - 1, SIZE - W_SIZE, -1
                SUM = Y(I)
                DO J_INDEX = STARTS(I + 1) - 1, DIAG_INDICES(I) + 1, -1
                    J = INDICES(J_INDEX)
                    SUM = SUM - A(J_INDEX) * X(J)
                END DO
                X(I) = SUM / A(J_INDEX)
            END DO

            DO LEVEL = LEVEL_STARTS - 2, 0, -1
                LEVEL_SIZE = INDEPENDENT_ROW_STARTS(LEVEL + 1) - INDEPENDENT_ROW_STARTS(LEVEL)
                IF (LEVEL_SIZE .LE. 10) THEN
                    CALL OMP_SET_NUM_THREADS(1)
                ELSE
                    CALL OMP_SET_NUM_THREADS(THREADS)
                END IF
                !$OMP PARALLEL
                    !$OMP DO SCHEDULE(STATIC) PRIVATE(LEVEL_INDEX, I, SUM, J_INDEX, J)
                DO LEVEL_INDEX = INDEPENDENT_ROW_STARTS(LEVEL + 1) - 1, INDEPENDENT_ROW_STARTS(LEVEL), -1
                    I = INDEPENDENT_ROWS(LEVEL_INDEX)
                !DO I = SIZE - W_SIZE - 1, 0, -1
                    SUM = Y(I)
                    DO J_INDEX = STARTS(I + 1) - 1, DIAG_INDICES(I) + 1, -1
                        J = INDICES(J_INDEX)
                        SUM = SUM - A(J_INDEX) * X(J)
                    END DO
                    X(I) = SUM / A(J_INDEX)
                END DO
                    !$OMP END DO
                !$OMP END PARALLEL
            END DO
            
            IF (ANS .LT. (MAX_ITER - 1)) THEN
                IF (CRS_CONVERGED_(SIZE, X, PARAMS%COPY, 1.D-7, 1.D-12)) THEN
                    EXIT
                END IF
            END IF            
        END DO
        CALL OMP_SET_NUM_THREADS(THREADS)
        ANS = ANS + 1        
    END FUNCTION

    FUNCTION GS_SOLVE_PARALLEL_BACKWARD_(SIZE, W_SIZE, LEVEL_STARTS, LENGTH, THREADS, STARTS, INDICES, DIAG_INDICES, INDEPENDENT_ROW_STARTS, INDEPENDENT_ROWS, A, B, Y, X, MAX_ITER, PARAMS) RESULT(ANS)
        INTEGER, INTENT(IN) :: SIZE, W_SIZE, LENGTH, LEVEL_STARTS, THREADS, MAX_ITER
        INTEGER, DIMENSION(0:SIZE), INTENT(IN) :: STARTS
        INTEGER, DIMENSION(0:LENGTH - 1), INTENT(IN) :: INDICES
        INTEGER, DIMENSION(0:SIZE - 1), INTENT(IN) :: DIAG_INDICES
        INTEGER, DIMENSION(0:LEVEL_STARTS - 1), INTENT(IN) :: INDEPENDENT_ROW_STARTS
        INTEGER, DIMENSION(0:SIZE - W_SIZE - 1), INTENT(IN) :: INDEPENDENT_ROWS
        DOUBLE PRECISION, DIMENSION(0:LENGTH - 1), INTENT(IN) :: A
        DOUBLE PRECISION, DIMENSION(0:SIZE - 1), INTENT(IN) :: B
        DOUBLE PRECISION, DIMENSION(0:SIZE - 1), INTENT(INOUT) :: Y
        DOUBLE PRECISION, DIMENSION(0:SIZE - 1), INTENT(OUT) :: X
        TYPE(GS_DATA), INTENT(INOUT) :: PARAMS
        INTEGER :: ANS
        INTEGER :: LEVEL, LEVEL_SIZE, LEVEL_INDEX
        INTEGER :: I, J_INDEX, J, DIAG
        INTEGER :: S, E, INC
        DOUBLE PRECISION :: SUM
        
        DO ANS = 0, MAX_ITER - 1      
            CALL CRS_COPY_(SIZE, X, PARAMS%COPY)            
            DO I = SIZE - 1, SIZE - W_SIZE, -1
                DIAG = DIAG_INDICES(I)
                SUM = B(I)
                DO J_INDEX = STARTS(I + 1) - 1, DIAG + 1, -1
                    J = INDICES(J_INDEX)
                    SUM = SUM - A(J_INDEX) * X(J)
                END DO
                DO J_INDEX = DIAG - 1, STARTS(I), -1
                    J = INDICES(J_INDEX)
                    SUM = SUM - A(J_INDEX) * X(J)
                END DO                
                X(I) = SUM / A(DIAG)
            END DO

            DO LEVEL = LEVEL_STARTS - 2, 0, -1
                LEVEL_SIZE = INDEPENDENT_ROW_STARTS(LEVEL + 1) - INDEPENDENT_ROW_STARTS(LEVEL)
                IF (LEVEL_SIZE .LE. 10) THEN
                    CALL OMP_SET_NUM_THREADS(1)
                ELSE
                    CALL OMP_SET_NUM_THREADS(THREADS)
                END IF
                !$OMP PARALLEL
                    !$OMP DO SCHEDULE(STATIC) PRIVATE(LEVEL_INDEX, I, DIAG, SUM, J_INDEX, J)
                DO LEVEL_INDEX = INDEPENDENT_ROW_STARTS(LEVEL + 1) - 1, INDEPENDENT_ROW_STARTS(LEVEL), -1
                    I = INDEPENDENT_ROWS(LEVEL_INDEX)
                    DIAG = DIAG_INDICES(I)
                    SUM = B(I)
                    DO J_INDEX = STARTS(I + 1) - 1, DIAG + 1, -1
                        J = INDICES(J_INDEX)
                        SUM = SUM - A(J_INDEX) * X(J)
                    END DO
                    DO J_INDEX = DIAG - 1, STARTS(I), -1
                        J = INDICES(J_INDEX)
                        SUM = SUM - A(J_INDEX) * X(J)
                    END DO
                    X(I) = SUM / A(DIAG)
                END DO
                    !$OMP END DO
                !$OMP END PARALLEL
            END DO
            
            IF (ANS .LT. (MAX_ITER - 1)) THEN
                IF (CRS_CONVERGED_(SIZE, X, PARAMS%COPY, 1.D-7, 1.D-12)) THEN
                    EXIT
                END IF
            END IF            
        END DO
        CALL OMP_SET_NUM_THREADS(THREADS)
        ANS = ANS + 1        
    END FUNCTION
    
    FUNCTION GS_SOLVE_(SIZE, LENGTH, STARTS, INDICES, DIAG_INDICES, A, B, X, MAX_ITER, PARAMS, IS_FORWARD, FORCED_C) RESULT(ANS)    
        INTEGER, INTENT(IN) :: SIZE, LENGTH, MAX_ITER
        INTEGER, DIMENSION(0:SIZE), INTENT(IN) :: STARTS
        INTEGER, DIMENSION(0:LENGTH - 1), INTENT(IN) :: INDICES
        INTEGER, DIMENSION(0:SIZE - 1), INTENT(IN) :: DIAG_INDICES
        DOUBLE PRECISION, DIMENSION(0:LENGTH - 1), INTENT(IN) :: A
        DOUBLE PRECISION, DIMENSION(0:SIZE - 1), INTENT(IN) :: B
        DOUBLE PRECISION, DIMENSION(0:SIZE - 1), INTENT(INOUT) :: X
        TYPE(GS_DATA), INTENT(INOUT) :: PARAMS
        LOGICAL, INTENT(IN) :: IS_FORWARD
        LOGICAL, DIMENSION(0:SIZE - 1), INTENT(IN) :: FORCED_C
        INTEGER :: ANS
        INTEGER :: I, J_INDEX, J, DIAG
        INTEGER :: S, E, INC
        DOUBLE PRECISION :: SUM
        IF (IS_FORWARD) THEN
            S = 0
            E = SIZE - 1
            INC = 1
        ELSE
            S = SIZE - 1
            E = 0
            INC = -1
        END IF        
        DO ANS = 0, MAX_ITER - 1      
            CALL CRS_COPY_(SIZE, X, PARAMS%COPY)
            DO I = S, E, INC
                IF (.NOT. FORCED_C(I)) THEN
                    DIAG = DIAG_INDICES(I)
                    SUM = B(I)
                    DO J_INDEX = STARTS(I), DIAG - 1
                        J = INDICES(J_INDEX)
                        SUM = SUM - A(J_INDEX) * X(J)
                    END DO
                    DO J_INDEX = DIAG + 1, STARTS(I + 1) - 1                    
                        J = INDICES(J_INDEX)
                        SUM = SUM - A(J_INDEX) * X(J)
                    END DO
                    X(I) = SUM / A(DIAG)                
                END IF
            END DO
            IF (ANS .LT. (MAX_ITER - 1)) THEN
                IF (CRS_CONVERGED_(SIZE, X, PARAMS%COPY, 1.D-7, 1.D-12)) THEN
                    EXIT
                END IF
            END IF            
        END DO
        ANS = ANS + 1
    END FUNCTION

    FUNCTION GS_CF_SOLVE_(SIZE, LENGTH, C_SIZE, F_SIZE, STARTS, INDICES, DIAG_INDICES, A, B, X, MAX_ITER, PARAMS, C, F, IS_FORWARD) RESULT(ANS)    
        INTEGER, INTENT(IN) :: SIZE, LENGTH, MAX_ITER, C_SIZE, F_SIZE
        INTEGER, DIMENSION(0:SIZE), INTENT(IN) :: STARTS
        INTEGER, DIMENSION(0:LENGTH - 1), INTENT(IN) :: INDICES
        INTEGER, DIMENSION(0:SIZE - 1), INTENT(IN) :: DIAG_INDICES
        DOUBLE PRECISION, DIMENSION(0:LENGTH - 1), INTENT(IN) :: A
        DOUBLE PRECISION, DIMENSION(0:SIZE - 1), INTENT(IN) :: B
        DOUBLE PRECISION, DIMENSION(0:SIZE - 1), INTENT(INOUT) :: X
        TYPE(GS_DATA), INTENT(INOUT) :: PARAMS
        LOGICAL, INTENT(IN) :: IS_FORWARD
        INTEGER, DIMENSION(0:C_SIZE - 1), INTENT(IN) :: C
        INTEGER, DIMENSION(0:F_SIZE - 1), INTENT(IN) :: F
        INTEGER :: ANS
        INTEGER :: I_INDEX, I, J_INDEX, J, DIAG
        INTEGER :: S, E, INC, CS, CE, FS, FE
        DOUBLE PRECISION :: SUM
        IF (IS_FORWARD) THEN
            S = 0
            CS = 0
            FS = 0
            E = SIZE - 1
            CE = C_SIZE - 1
            FE = F_SIZE - 1
            INC = 1
        ELSE
            S = SIZE - 1
            CS = C_SIZE - 1
            FS = F_SIZE - 1
            E = 0
            CE = 0
            FE = 0
            INC = -1
        END IF        
        DO ANS = 0, MAX_ITER - 1      
            CALL CRS_COPY_(SIZE, X, PARAMS%COPY)
            !$OMP PARALLEL SHARED(X, B, A)
                !$OMP DO SCHEDULE(STATIC) PRIVATE(I_INDEX, I, J_INDEX, J, SUM, DIAG)            
            DO I_INDEX = CS, CE, INC
                I = C(I_INDEX)
                DIAG = DIAG_INDICES(I)
                SUM = B(I)
                DO J_INDEX = STARTS(I), DIAG - 1
                    J = INDICES(J_INDEX)
                    SUM = SUM - A(J_INDEX) * X(J)
                END DO
                DO J_INDEX = DIAG + 1, STARTS(I + 1) - 1                    
                    J = INDICES(J_INDEX)
                    SUM = SUM - A(J_INDEX) * X(J)
                END DO
                X(I) = SUM / A(DIAG)
            END DO
                !$OMP END DO
            !$OMP END PARALLEL
            DO I_INDEX = FS, FE, INC
                I = F(I_INDEX)
                DIAG = DIAG_INDICES(I)
                SUM = B(I)
                DO J_INDEX = STARTS(I), DIAG - 1
                    J = INDICES(J_INDEX)
                    SUM = SUM - A(J_INDEX) * X(J)
                END DO
                DO J_INDEX = DIAG + 1, STARTS(I + 1) - 1                    
                    J = INDICES(J_INDEX)
                    SUM = SUM - A(J_INDEX) * X(J)
                END DO
                X(I) = SUM / A(DIAG)
            END DO            
            IF (ANS .LT. (MAX_ITER - 1)) THEN
                IF (CRS_CONVERGED_(SIZE, X, PARAMS%COPY, 1.D-7, 1.D-12)) THEN
                    EXIT
                END IF
            END IF            
        END DO
        ANS = ANS + 1
    END FUNCTION
    
END MODULE
    