!****************************************************************************************************************!
!*                                                                                                              *!    
!*                                              GENERAL UTILITIES                                               *!
!*                                                                                                              *!    
!*                                                                                                              *!    
!****************************************************************************************************************!

    SUBROUTINE INV22__(A)
        !DEC$ ATTRIBUTES DLLEXPORT, ALIAS : 'INV22__':: RETINA_SOLVER
        DOUBLE PRECISION, DIMENSION(0:3), INTENT(INOUT) :: A
        DOUBLE PRECISION :: A11, A12, A21, A22, DET
        A11 = A(0)
        A12 = A(1)
        A21 = A(2)
        A22 = A(3)
        DET = A11 * A22 - A12 * A21
        A(0) = A22 / DET
        A(1) = -A12 / DET
        A(2) = -A21 / DET
        A(3) = A11 / DET
    END SUBROUTINE
    
    SUBROUTINE INV__(A)
        !DEC$ ATTRIBUTES DLLEXPORT, ALIAS : 'INV__':: RETINA_SOLVER
        DOUBLE PRECISION, DIMENSION(0:8), INTENT(INOUT) :: A
        DOUBLE PRECISION, DIMENSION(0:8) :: LU
        DOUBLE PRECISION :: A11, A12, A13, A21, A22, A23, A31, A32, A33, DET
        INTEGER :: I, J, PIVOT, I_MAX, P_VAL, COL
        INTEGER, DIMENSION(0:8) :: P
        DOUBLE PRECISION :: PIVOT_VAL, VAL, SUM, MULT
        INTEGER :: C0, C1, C2
        !DEC$ SIMD
        DO I = 0, 2
            C0 = 3 * I
            DO J = 0, 2
                P(C0 + J) = 0.D0
                LU(C0 + J) = A(C0 + J)
                A(C0 + J) = 0.D0
            END DO
            P(C0 + I) = 1.D0
        END DO
        DO PIVOT = 0, 1
            C0 = 3 * PIVOT
            PIVOT_VAL = LU(C0 + PIVOT)
            I_MAX = PIVOT
            DO I = PIVOT + 1, 2
                C1 = 3 * I
                VAL = LU(C1 + PIVOT)
                IF (ABS(VAL) .GT. ABS(PIVOT_VAL)) THEN
                    PIVOT_VAL = VAL
                    I_MAX = I
                END IF
            END DO
            !DEC$ SIMD
            DO J = 0, 2
                C0 = 3 * PIVOT
                C1 = 3 * I_MAX
                VAL = LU(C0 + J)
                P_VAL = P(C0 + J)
                
                LU(C0 + J) = LU(C1 + J)
                P(C0 + J) = P(C1 + J)
                LU(C1 + J) = VAL
                P(C1 + J) = P_VAL
            END DO
            DO I = PIVOT + 1, 2
                C0 = 3 * I
                MULT = LU(C0 + PIVOT) / PIVOT_VAL
                LU(C0 + PIVOT) = MULT
                C1 = 3 * PIVOT
                !DEC$ SIMD
                DO J = PIVOT + 1, 2
                    LU(C0 + J) = LU(C0 + J) - MULT * LU(C1 + J)
                END DO
            END DO
        END DO
        DO COL = 0, 2
            !DEC$ SIMD
            DO I = 0, 2
                C0 = 3 * I
                SUM = P(C0 + COL)                
                DO J = 0, I - 1
                    C1 = 3 * J
                    SUM = SUM - LU(C0 + J) * A(C1 + COL)
                END DO
                A(C0 + COL) = SUM
            END DO
            !DEC$ SIMD
            DO I = 2, 0, -1
                C0 = 3 * I
                SUM = A(C0 + COL)
                DO J = I + 1, 2
                    C1 = 3 * J
                    SUM = SUM - LU(C0 + J) * A(C1 + COL)
                END DO
                A(3 * I + COL) = SUM / LU(3 * I + I)
            END DO
        END DO
        !WRITE (*,*) 'DEBUG'        
        !A11 = A(3 * 0 + 0)
        !A12 = A(3 * 0 + 1)
        !A13 = A(3 * 0 + 2)
        !A21 = A(3 * 1 + 0)
        !A22 = A(3 * 1 + 1)
        !A23 = A(3 * 1 + 2)
        !A31 = A(3 * 2 + 0)
        !A32 = A(3 * 2 + 1)
        !A33 = A(3 * 2 + 2)
        !
        !DET = A11 * (A22 * A33 - A23 * A32) - A12 * (A21 * A33 - A23 * A31) + A13 * (A21 * A32 - A22 * A31)
        !
        !A(3 * 0 + 0) = (A22 * A33 - A23 * A32) / DET
        !A(3 * 0 + 1) = (-A12 * A33 + A13 * A32) / DET
        !A(3 * 0 + 2) = (A12 * A23 - A13 * A22) / DET
        !A(3 * 1 + 0) = (-A21 * A33 + A23 * A31) / DET
        !A(3 * 1 + 1) = (A11 * A33 - A31 * A13) / DET
        !A(3 * 1 + 2) = (-A11 * A23 + A13 * A21) / DET
        !A(3 * 2 + 0) = (A21 * A32 - A22 * A31) / DET
        !A(3 * 2 + 1) = (-A11 * A32 + A31 * A12) / DET
        !A(3 * 2 + 2) = (A11 * A22 - A21 * A12) / DET
        
    END SUBROUTINE

    SUBROUTINE COPY_MULT_MTM__(A, B, C, SCALER)
        !DEC$ ATTRIBUTES DLLEXPORT, ALIAS : 'COPY_MULT_MTM__':: RETINA_SOLVER
        DOUBLE PRECISION, INTENT(IN) :: SCALER
        DOUBLE PRECISION, DIMENSION(0:8), INTENT(IN) :: A, B
        DOUBLE PRECISION, DIMENSION(0:8), INTENT(INOUT) :: C
        DOUBLE PRECISION, DIMENSION(0:8) :: TEMP
        INTEGER :: II, JJ, KK, C1_A, C1_C, C2_B
        !DEC$ SIMD
        DO II = 0, 8
            TEMP(II) = 0.D0
        END DO
        !DEC$ SIMD
        DO II = 0, 2
            C1_C = 3 * II
            DO JJ = 0, 2
                DO KK = 0, 2
                    C1_A = 3 * KK
                    C2_B = 3 * KK
                    TEMP(C1_C + JJ) = TEMP(C1_C + JJ) + SCALER * A(C1_A + II) * B(C2_B + JJ)
                END DO
            END DO
        END DO                    
        !DEC$ SIMD
        DO II = 0, 8
            C(II) = TEMP(II)
        END DO
    END SUBROUTINE
    
    SUBROUTINE COPY_MULT_MM__(A, B, C, SCALER)
        !DEC$ ATTRIBUTES DLLEXPORT, ALIAS : 'COPY_MULT_MM__':: RETINA_SOLVER
        DOUBLE PRECISION, INTENT(IN) :: SCALER
        DOUBLE PRECISION, DIMENSION(0:8), INTENT(IN) :: A, B
        DOUBLE PRECISION, DIMENSION(0:8), INTENT(INOUT) :: C
        DOUBLE PRECISION, DIMENSION(0:8) :: TEMP
        INTEGER :: II, JJ, KK, C1_A, C1_C, C2_B
        !DEC$ SIMD
        DO II = 0, 8
            TEMP(II) = 0.D0
        END DO
        !DEC$ SIMD
        DO II = 0, 2
            C1_C = 3 * II
            C1_A = 3 * II
            DO JJ = 0, 2
                DO KK = 0, 2
                    C2_B = 3 * KK
                    TEMP(C1_C + JJ) = TEMP(C1_C + JJ) + SCALER * A(C1_A + KK) * B(C2_B + JJ)
                END DO
            END DO
        END DO                    
        !DEC$ SIMD
        DO II = 0, 8
            C(II) = TEMP(II)
        END DO
    END SUBROUTINE

    SUBROUTINE COPY_MULT_MV__(A, B, C, SCALER)
        !DEC$ ATTRIBUTES DLLEXPORT, ALIAS : 'COPY_MULT_MV__':: RETINA_SOLVER
        DOUBLE PRECISION, INTENT(IN) :: SCALER
        DOUBLE PRECISION, DIMENSION(0:8), INTENT(IN) :: A
        DOUBLE PRECISION, DIMENSION(0:2), INTENT(IN) :: B
        DOUBLE PRECISION, DIMENSION(0:2), INTENT(INOUT) :: C
        DOUBLE PRECISION, DIMENSIOn(0:2) :: TEMP
        INTEGER :: II, JJ, KK, C1_A
        !DEC$ SIMD
        DO II = 0, 2
            TEMP(II) = 0.D0
        END DO
        !DEC$ SIMD
        DO II = 0, 2
            C1_A = 3 * II
            DO JJ = 0, 2
                TEMP(II) = TEMP(II) + SCALER * A(C1_A + JJ) * B(JJ)
            END DO
        END DO               
        !DEC$ SIMD
        DO II = 0, 2
            C(II) = TEMP(II)
        END DO
    END SUBROUTINE

    SUBROUTINE MULT_MTM_TRANSPOSE__(A, B, C, SCALER)
        !DEC$ ATTRIBUTES DLLEXPORT, ALIAS : 'MULT_MTM_TRANSPOSE__':: RETINA_SOLVER
        DOUBLE PRECISION, INTENT(IN) :: SCALER
        DOUBLE PRECISION, DIMENSION(0:8), INTENT(IN) :: A, B
        DOUBLE PRECISION, DIMENSION(0:8), INTENT(INOUT) :: C
        INTEGER :: II, JJ, KK, C1_A, C1_C, C2_B
        !DEC$ SIMD
        DO II = 0, 2
            DO JJ = 0, 2
                C1_C = 3 * JJ
                DO KK = 0, 2
                    C2_B = 3 * KK
                    C1_A = 3 * KK
                    C(C1_C + II) = C(C1_C + II) + SCALER * A(C1_A + II) * B(C2_B + JJ)
                END DO
            END DO
        END DO                    
    END SUBROUTINE
    
    SUBROUTINE MULT_MTM__(A, B, C, SCALER)
        !DEC$ ATTRIBUTES DLLEXPORT, ALIAS : 'MULT_MTM__':: RETINA_SOLVER
        DOUBLE PRECISION, INTENT(IN) :: SCALER
        DOUBLE PRECISION, DIMENSION(0:8), INTENT(IN) :: A, B
        DOUBLE PRECISION, DIMENSION(0:8), INTENT(INOUT) :: C
        INTEGER :: II, JJ, KK, C1_A, C1_C, C2_B
        !DEC$ SIMD
        DO II = 0, 2
            C1_C = 3 * II
            DO JJ = 0, 2
                DO KK = 0, 2
                    C2_B = 3 * KK
                    C1_A = 3 * KK
                    C(C1_C + JJ) = C(C1_C + JJ) + SCALER * A(C1_A + II) * B(C2_B + JJ)
                END DO
            END DO
        END DO                    
    END SUBROUTINE
    
    SUBROUTINE MULT_MM__(A, B, C, SCALER, MIN_II, MIN_JJ, MIN_KK)
        !DEC$ ATTRIBUTES DLLEXPORT, ALIAS : 'MULT_MM__':: RETINA_SOLVER
        DOUBLE PRECISION, INTENT(IN) :: SCALER
        DOUBLE PRECISION, DIMENSION(0:8), INTENT(IN) :: A, B
        DOUBLE PRECISION, DIMENSION(0:8), INTENT(INOUT) :: C
        INTEGER , INTENT(IN) :: MIN_II, MIN_JJ, MIN_KK
        INTEGER :: II, JJ, KK, C1_A, C1_C, C2_B
        !DEC$ SIMD
        DO II = MIN_II, 2
            C1_C = 3 * II
            C1_A = 3 * II
            DO JJ = MIN_JJ, 2
                DO KK = MIN_KK, 2
                    C2_B = 3 * KK
                    C(C1_C + JJ) = C(C1_C + JJ) + SCALER * A(C1_A + KK) * B(C2_B + JJ)
                END DO
            END DO
        END DO                    
    END SUBROUTINE

    SUBROUTINE MULT_MM_DIAG__(A, B, C, SCALER)
        !DEC$ ATTRIBUTES DLLEXPORT, ALIAS : 'MULT_MM_DIAG__':: RETINA_SOLVER
        DOUBLE PRECISION, INTENT(IN) :: SCALER
        DOUBLE PRECISION, DIMENSION(0:8), INTENT(IN) :: A, B
        DOUBLE PRECISION, DIMENSION(0:8), INTENT(INOUT) :: C
        INTEGER :: II, JJ, KK, C1_A, C1_C, C2_B
        !DEC$ SIMD
        DO II = 0, 2
            C1_C = 3 * II
            C1_A = 3 * II
            DO JJ = 0, 2
                DO KK = 0, 2
                    C2_B = 3 * KK
                    C(C1_C + II) = C(C1_C + II) + SCALER * A(C1_A + KK) * B(C2_B + JJ)
                END DO
            END DO
        END DO                    
    END SUBROUTINE
    
    SUBROUTINE MULT_MV__(A, B, C, SCALER, MIN_II, MIN_JJ)
        !DEC$ ATTRIBUTES DLLEXPORT, ALIAS : 'MULT_MV__':: RETINA_SOLVER    
        DOUBLE PRECISION, INTENT(IN) :: SCALER
        DOUBLE PRECISION, DIMENSION(0:8), INTENT(IN) :: A
        DOUBLE PRECISION, DIMENSION(0:2), INTENT(IN) :: B
        DOUBLE PRECISION, DIMENSION(0:2), INTENT(INOUT) :: C
        INTEGER, INTENT(IN) :: MIN_II, MIN_JJ
        INTEGER :: II, JJ, KK, C1_A
        !DEC$ SIMD
        DO II = MIN_II, 2
            C1_A = 3 * II
            DO JJ = MIN_JJ, 2
                C(II) = C(II) + SCALER * A(C1_A + JJ) * B(JJ)
            END DO
        END DO                    
    END SUBROUTINE
    
    SUBROUTINE COPY_RANGE__(I1, I2, SIZE, SOURCE, DEST)
        !DEC$ ATTRIBUTES DLLEXPORT, ALIAS : 'COPY_RANGE__':: RETINA_SOLVER
        INTEGER, INTENT(IN) :: I1, I2, SIZE
        DOUBLE PRECISION, DIMENSION(0:9 * SIZE - 1), INTENT(IN) :: SOURCE
        DOUBLE PRECISION, DIMENSION(0:9 * SIZE - 1), INTENT(OUT) :: DEST
        INTEGER :: I, C0, II
        !DEC$ SIMD
        DO I = I1, I2
            C0 = 9 * I
            DO II = 0, 8
                DEST(C0 + II) = SOURCE(C0 + II)
            END DO
        END DO
    END SUBROUTINE
    
    SUBROUTINE CLEAR_RANGE_MATRIX__(I1, I2, SIZE, DEST)
        !DEC$ ATTRIBUTES DLLEXPORT, ALIAS : 'CLEAR_RANGE_MATRIX__':: RETINA_SOLVER
        INTEGER, INTENT(IN) :: I1, I2, SIZE
        DOUBLE PRECISION, DIMENSION(0:9 * SIZE - 1), INTENT(OUT) :: DEST
        INTEGER :: II
        !DEC$ SIMD
        DO II = 9 * I1, 9 * I2 + 8
            DEST(II) = 0.D0
        END DO
    END SUBROUTINE

    SUBROUTINE CLEAR_RANGE_VECTOR__(I1, I2, SIZE, DEST)
        !DEC$ ATTRIBUTES DLLEXPORT, ALIAS : 'CLEAR_RANGE_VECTOR__':: RETINA_SOLVER
        INTEGER, INTENT(IN) :: I1, I2, SIZE
        DOUBLE PRECISION, DIMENSION(0:3 * SIZE - 1), INTENT(OUT) :: DEST
        INTEGER :: II
        !DEC$ SIMD
        DO II = 3 * I1, 3 * I2 + 2
            DEST(II) = 0.D0
        END DO
    END SUBROUTINE
    
    SUBROUTINE CLONE__(SIZE, DEST, SOURCE)
        !DEC$ ATTRIBUTES DLLEXPORT, ALIAS : 'CLONE__' :: RETINA_SOLVER
        USE OMP_LIB
        INTEGER, INTENT(IN) :: SIZE
        DOUBLE PRECISION, DIMENSION(0:SIZE - 1), INTENT(IN) :: SOURCE
        DOUBLE PRECISION, DIMENSION(0:SIZE - 1), INTENT(OUT) :: DEST
        INTEGER :: I
        !DEC$ SIMD
        DO I = 0, SIZE - 1
            DEST(I) = SOURCE(I)
        END DO
    END SUBROUTINE
        
    SUBROUTINE ROW_SUM_(I1, I2, SOURCE_SIZE, DEST_SIZE, ROW_STARTS, SOURCE, DEST)
        INTEGER, INTENT(IN) :: I1, I2, SOURCE_SIZE, DEST_SIZE
        INTEGER, DIMENSION(0:DEST_SIZE), INTENT(IN) :: ROW_STARTS
        DOUBLE PRECISION, DIMENSION(0:9 * SOURCE_SIZE - 1), INTENT(IN) :: SOURCE
        DOUBLE PRECISION, DIMENSION(0:9 * DEST_SIZE - 1), INTENT(OUT) :: DEST
        INTEGER :: I, C0, C1, II, J_INDEX
        !DEC$ SIMD
        DO I = I1, I2
            C0 = 9 * I
            DO II = 0, 8
                DEST(C0 + II) = 0.D0
            END DO
        END DO
        DO I = I1, I2
            C0 = 9 * I
            DO J_INDEX = ROW_STARTS(I), ROW_STARTS(I + 1) - 1
                C1 = 9 * J_INDEX
                !DEC$ SIMD
                DO II = 0, 8
                    DEST(C0 + II) = DEST(C0 + II) + SOURCE(C1 + II)
                END DO                
            END DO
        END DO        
    END SUBROUTINE
    
    SUBROUTINE ROW_SUM_TRANSPOSE_(DEST1, DEST2, ROW_START1, ROW_START2, SOURCE_SIZE, DEST_SIZE, ROW_START_SIZE, ROW_STARTS, CONNECTIVITY, SOURCE, DEST)        
        INTEGER, INTENT(IN) :: DEST1, DEST2, ROW_START1, ROW_START2, SOURCE_SIZE, DEST_SIZE, ROW_START_SIZE
        INTEGER, DIMENSION(0:ROW_START_SIZE), INTENT(IN) :: ROW_STARTS
        INTEGER, DIMENSION(0:SOURCE_SIZE - 1), INTENT(IN) :: CONNECTIVITY
        DOUBLE PRECISION, DIMENSION(0:9 * SOURCE_SIZE - 1), INTENT(IN) :: SOURCE
        DOUBLE PRECISION, DIMENSION(0:9 * DEST_SIZE - 1), INTENT(OUT) :: DEST
        INTEGER :: I, C0_DEST, C1_DEST, C0_SOURCE, C1_SOURCE, II, JJ, J_INDEX, J
        !DEC$ SIMD
        DO I = DEST1, DEST2
            C0_DEST = 9 * (I)
            DO II = 0, 8
                DEST(C0_DEST + II) = 0.D0
            END DO
        END DO
        DO I = ROW_START1, ROW_START2
            DO J_INDEX = ROW_STARTS(I), ROW_STARTS(I + 1) - 1
                C0_SOURCE = 9 * J_INDEX
                J = CONNECTIVITY(J_INDEX)
                C0_DEST = 9 * J
                !DEC$ SIMD
                DO II = 0, 2
                    C1_DEST = C0_DEST + 3 * II
                    DO JJ = 0, 2
                        C1_SOURCE = C0_SOURCE + 3 * JJ
                        DEST(C1_DEST + JJ) = DEST(C1_DEST + JJ) + SOURCE(C1_SOURCE + II)
                    END DO
                END DO
            END DO
        END DO
    END SUBROUTINE
    
        