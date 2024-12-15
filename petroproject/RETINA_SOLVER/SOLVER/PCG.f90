MODULE PCG
    
    USE CRS
    USE ILU0
    IMPLICIT NONE
    
    TYPE :: PCG_DATA
        DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: R, Q, P, Z, COPY
    END TYPE
    
CONTAINS

    SUBROUTINE PCG_UPDATE_P_(SIZE, P, Z, BETA_I_1)
        INTEGER, INTENT(IN) :: SIZE
        DOUBLE PRECISION, DIMENSION(0:SIZE - 1), INTENT(IN) :: Z
        DOUBLE PRECISION, DIMENSION(0:SIZE - 1), INTENT(INOUT) :: P
        DOUBLE PRECISION :: BETA_I_1
        INTEGER :: I
        !DEC$ SIMD
        DO I = 0, SIZE - 1
            P(I) = Z(I) + BETA_I_1 * P(I)
        END DO
    END SUBROUTINE

    SUBROUTINE PCG_UPDATE_X_(SIZE, X, P, ALPHA_I)
        INTEGER, INTENT(IN) :: SIZE
        DOUBLE PRECISION, DIMENSION(0:SIZE - 1), INTENT(IN) :: P
        DOUBLE PRECISION, DIMENSION(0:SIZE - 1), INTENT(INOUT) :: X
        DOUBLE PRECISION :: ALPHA_I
        INTEGER :: I
        !DEC$ SIMD
        DO I = 0, SIZE - 1
            X(I) = X(I) + ALPHA_I * P(I)
        END DO
    END SUBROUTINE

    SUBROUTINE PCG_UPDATE_R_(SIZE, R, Q, ALPHA_I)
        INTEGER, INTENT(IN) :: SIZE
        DOUBLE PRECISION, DIMENSION(0:SIZE - 1), INTENT(IN) :: Q
        DOUBLE PRECISION, DIMENSION(0:SIZE - 1), INTENT(INOUT) :: R
        DOUBLE PRECISION :: ALPHA_I
        INTEGER :: I
        !DEC$ SIMD
        DO I = 0, SIZE - 1
            R(I) = R(I) - ALPHA_I * Q(I)
        END DO
    END SUBROUTINE
    
    FUNCTION PCG_SOLVE_(SIZE, LENGTH, M_LENGTH, STARTS, INDICES, M_STARTS, M_INDICES, M_DIAG_INDICES, A, M, B, Y, X, MAX_ITER, PARAMS) RESULT(ANS)    
        INTEGER, INTENT(IN) :: SIZE, LENGTH, MAX_ITER, M_LENGTH
        INTEGER, DIMENSION(0:SIZE), INTENT(IN) :: STARTS, M_STARTS
        INTEGER, DIMENSION(0:LENGTH - 1), INTENT(IN) :: INDICES
        INTEGER, DIMENSION(0:M_LENGTH - 1), INTENT(IN) :: M_INDICES
        INTEGER, DIMENSION(0:SIZE - 1), INTENT(IN) :: M_DIAG_INDICES
        DOUBLE PRECISION, DIMENSION(0:LENGTH - 1), INTENT(IN) :: A, M
        DOUBLE PRECISION, DIMENSION(0:SIZE - 1), INTENT(IN) :: B
        DOUBLE PRECISION, DIMENSION(0:SIZE - 1), INTENT(INOUT) :: Y, X
        TYPE(PCG_DATA), INTENT(INOUT) :: PARAMS
        INTEGER :: ANS
        DOUBLE PRECISION :: RO_I_1, RO_I_2, ALPHA_I, BETA_I_1
        INTEGER :: I, II
        RO_I_1 = 0.D0
        CALL CRS_RESIDUAL_(SIZE, SIZE, LENGTH, STARTS, INDICES, A, B, X, PARAMS%R)
        if (size .EQ. 9026) then
            write (*,*) 'PCG_STAB'
            write (*,*) CRS_VECTOR_NORM_(size, params%r)
        end if
        DO I = 0, MAX_ITER - 1      
            CALL CRS_COPY_(SIZE, X, PARAMS%COPY)
            CALL ILU0_TRI_SOLVE_(SIZE, M_LENGTH, M_STARTS, M_INDICES, M_DIAG_INDICES, M, PARAMS%R, PARAMS%Z, Y)
            RO_I_2 = RO_I_1
            RO_I_1 = CRS_INNER_PRODUCT_(SIZE, PARAMS%R, PARAMS%Z)
            IF (I .EQ. 0) THEN
                CALL CRS_COPY_(SIZE, PARAMS%Z, PARAMS%P)            
            ELSE
                BETA_I_1 = RO_I_1 / RO_I_2
                CALL PCG_UPDATE_P_(SIZE, PARAMS%P, PARAMS%Z, BETA_I_1)
            END IF
            CALL CRS_MULT_MV_(SIZE, SIZE, LENGTH, STARTS, INDICES, A, PARAMS%P, PARAMS%Q)
            ALPHA_I = RO_I_1 / CRS_INNER_PRODUCT_(SIZE, PARAMS%P, PARAMS%Q)
            CALL PCG_UPDATE_X_(SIZE, X, PARAMS%P, ALPHA_I)        
            CALL PCG_UPDATE_R_(SIZE, PARAMS%R, PARAMS%Q, ALPHA_I)    
            if (size .EQ. 9026) then
                write (*,*) 'PCG_STAB'
                write (*,*) CRS_VECTOR_NORM_(size, params%r)
            end if
            IF (I .LT. (MAX_ITER - 1)) THEN
                IF (CRS_CONVERGED_(SIZE, X, PARAMS%COPY, 1.D-7, 1.D-12)) THEN
                    EXIT
                END IF
            END IF
        END DO
        ANS = I + 1
        if (size .EQ. 9026) then
            write (*,*) 'PCG_STAB FINISHED'
        end if
    END FUNCTION
        
END MODULE
    