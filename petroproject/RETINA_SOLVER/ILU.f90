!****************************************************************************************************************!
!*                                                                                                              *!    
!*                                              ILU PRECONDITIONER                                              *!
!*                                                                                                              *!    
!*                                                                                                              *!    
!****************************************************************************************************************!

    FUNCTION GET_INDEX_(I_REVERSE_INDEX, J, ACTIVE_SIZE, ILU_PARALLEL_SIZE, ILU_PARALLEL_CONNECTIVITY, ILU_PARALLEL_ROW_STARTS) RESULT(ANS)
        INTEGER, INTENT(IN) :: I_REVERSE_INDEX, J, ACTIVE_SIZE, ILU_PARALLEL_SIZE
        INTEGER, DIMENSION(0:ACTIVE_SIZE), INTENT(IN) :: ILU_PARALLEL_ROW_STARTS
        INTEGER, DIMENSION(0:ILU_PARALLEL_SIZE - 1), INTENT(IN) :: ILU_PARALLEL_CONNECTIVITY
        INTEGER :: ANS
        INTEGER :: START_INDEX, END_INDEX, ID, MED_INDEX
        START_INDEX = ILU_PARALLEL_ROW_STARTS(I_REVERSE_INDEX)
        END_INDEX = ILU_PARALLEL_ROW_STARTS(I_REVERSE_INDEX + 1) - 1
        ANS = -1
        DO WHILE (START_INDEX .LE. END_INDEX)
            MED_INDEX = (START_INDEX + END_INDEX) / 2
            ID = ILU_PARALLEL_CONNECTIVITY(MED_INDEX)
            IF (ID .EQ. J) THEN
                ANS = MED_INDEX
                RETURN
            ELSE IF (ID .GT. J) THEN
                END_INDEX = MED_INDEX - 1
            ELSE
                START_INDEX = MED_INDEX + 1
            END IF                            
        END DO                	
    END FUNCTION
                                      !, CONV_CV_SIZE                              
    SUBROUTINE AIM_ILU_SOLVE__(CV_SIZE, ACTIVE_SIZE, W_SIZE, CON_SIZE, ILU_PARALLEL_SIZE, ILU_LEVEL_STARTS, THREADS, X, B, Y, ILU_FLOW_VALS, ILU_C2W_VALS, ILU_W2C_VALS, ILU_W_VALS, ILU_PARALLEL_CONNECTIVITY, ILU_PARALLEL_ROW_STARTS, ILU_DIAG_INDICES, INDEPENDENT_ROWS, INDEPENDENT_ROW_STARTS, DOF_WELL_CONNECTIVITY, DOF_WELL_ROW_STARTS, WELL_DOF_CONNECTIVITY, WELL_DOF_ROW_STARTS, IMPLICITNESS, NATURAL_ID)
        !DEC$ ATTRIBUTES DLLEXPORT, ALIAS : 'AIM_ILU_SOLVE__':: RETINA_SOLVER
        USE OMP_LIB        
        INTEGER, INTENT(IN) :: CV_SIZE, ACTIVE_SIZE, W_SIZE, CON_SIZE, ILU_PARALLEL_SIZE, ILU_LEVEL_STARTS, THREADS       !, CONV_CV_SIZE
        INTEGER, DIMENSION(0: CV_SIZE- 1), INTENT(IN) :: IMPLICITNESS, NATURAL_ID
        DOUBLE PRECISION, DIMENSION(0:3 * ACTIVE_SIZE + 3 * W_SIZE - 1), INTENT(IN) :: B
        DOUBLE PRECISION, DIMENSION(0:3 * ACTIVE_SIZE + 3 * W_SIZE - 1), INTENT(INOUT) :: Y
        DOUBLE PRECISION, DIMENSION(0:3 * ACTIVE_SIZE + 3 * W_SIZE - 1), INTENT(INOUT) :: X
        DOUBLE PRECISION, DIMENSION(0:9 * ILU_PARALLEL_SIZE - 1), INTENT(IN) :: ILU_FLOW_VALS
        INTEGER, DIMENSION(0:ILU_PARALLEL_SIZE - 1), INTENT(IN) :: ILU_PARALLEL_CONNECTIVITY
        INTEGER, DIMENSION(0:ACTIVE_SIZE), INTENT(IN) :: ILU_PARALLEL_ROW_STARTS
        INTEGER, DIMENSION(0:ACTIVE_SIZE - 1), INTENT(IN) :: ILU_DIAG_INDICES
        DOUBLE PRECISION, DIMENSION(0:9 * CON_SIZE - 1), INTENT(IN) :: ILU_C2W_VALS, ILU_W2C_VALS
        INTEGER, DIMENSION(0:CON_SIZE - 1), INTENT(IN) :: DOF_WELL_CONNECTIVITY, WELL_DOF_CONNECTIVITY
        INTEGER, DIMENSION(0:ACTIVE_SIZE), INTENT(IN) :: DOF_WELL_ROW_STARTS
        INTEGER, DIMENSION(0:W_SIZE), INTENT(IN) :: WELL_DOF_ROW_STARTS
        
        DOUBLE PRECISION, DIMENSION(0:9 * W_SIZE - 1), INTENt(IN) :: ILU_W_VALS
        INTEGER, DIMENSION(0:ILU_LEVEL_STARTS - 1), INTENT(IN) :: INDEPENDENT_ROW_STARTS
        INTEGER, DIMENSION(0:ACTIVE_SIZE - 1), INTENT(IN) :: INDEPENDENT_ROWS
        
        INTEGER :: I, J, K, LEVEL_INDEX, LEVEL, LEVEL_SIZE
        DOUBLE PRECISION, DIMENSION(0:8) :: DIAG_VAL_33
        DOUBLE PRECISION, DIMENSION(0:2) :: SUM_3, X_3        
        DOUBLE PRECISION :: DIAG_VAL_1, SUM_1, VALUE_1, X_1
        INTEGER :: J_INDEX, INDEX
        INTEGER :: C0, C1, C2, C3, II, JJ
        INTEGER :: SIZE, MIN_II, MIN_JJ
        
        !DIR$ ASSUME_ALIGNED ILU_FLOW_VALS: 64
        !DIR$ ASSUME_ALIGNED ILU_PARALLEL_CONNECTIVITY: 64
        !DIR$ ASSUME_ALIGNED ILU_PARALLEL_ROW_STARTS: 64
        !DIR$ ASSUME_ALIGNED X: 64
        !DIR$ ASSUME_ALIGNED B: 64
        !DIR$ ASSUME_ALIGNED Y: 64        
        
        !$OMP PARALLEL
        DO LEVEL = 0, ILU_LEVEL_STARTS - 2
            LEVEL_SIZE = INDEPENDENT_ROW_STARTS(LEVEL + 1) - INDEPENDENT_ROW_STARTS(LEVEL)
            !IF (LEVEL_SIZE .LE. 10) THEN
            !    CALL OMP_SET_NUM_THREADS(1)
            !ELSE
            !    CALL OMP_SET_NUM_THREADS(THREADS)
            !END IF
            !$OMP DO SCHEDULE(STATIC) PRIVATE(LEVEL_INDEX, I, MIN_II, C0, K, SUM_3, J_INDEX, J, MIN_JJ, C1, C2, II, C3, JJ)
            DO LEVEL_INDEX = INDEPENDENT_ROW_STARTS(LEVEL), INDEPENDENT_ROW_STARTS(LEVEL + 1) - 1
                I = INDEPENDENT_ROWS(LEVEL_INDEX)
                IF(IMPLICITNESS(NATURAL_ID(I)) .EQ. 1) THEN
                    MIN_II = 0
                ELSE
                    MIN_II = 2
                END IF
                C0 = 3 * I
                !DEC$ SIMD   
                DO K = MIN_II, 2
                    SUM_3(K) = B(C0 + K)
                END DO                
                DO J_INDEX = ILU_PARALLEL_ROW_STARTS(LEVEL_INDEX), ILU_DIAG_INDICES(LEVEL_INDEX) - 1
                    J = ILU_PARALLEL_CONNECTIVITY(J_INDEX)   
                    IF(IMPLICITNESS(NATURAL_ID(J)) .EQ. 1) THEN
                        MIN_JJ = 0
                    ELSE
                        MIN_JJ = 2
                    END IF
                    C1 = 3 * J
                    C2 = 9 * J_INDEX
                    !DEC$ SIMD
                    DO II = MIN_II, 2
                        C3 = C2 + 3 * II
                        DO JJ = MIN_JJ, 2
                            SUM_3(II) = SUM_3(II) - ILU_FLOW_VALS(C3 + JJ) * X(C1 + JJ)                             
                        END DO
                    END DO
                END DO   
                !DEC$ SIMD
                DO K = MIN_II, 2
                    X(C0 + K) = SUM_3(K)
                END DO
            END DO
            !$OMP END DO        
        END DO      
        !$OMP END PARALLEL
        DO I = 0, W_SIZE - 1
            C0 = 3 * ACTIVE_SIZE + 3 * I
            !DEC$ SIMD
            DO K = 0, 2
                SUM_3(K) = B(C0 + K)
            END DO
            DO J_INDEX = WELL_DOF_ROW_STARTS(I), WELL_DOF_ROW_STARTS(I + 1) - 1
                J = WELL_DOF_CONNECTIVITY(J_INDEX)                
                C1 = 3 * J
                C2 = 9 * J_INDEX
                !DEC$ SIMD
                DO II = 0, 2
                    C3 = C2 + 3 * II
                    DO JJ = 0, 2
                        SUM_3(II) = SUM_3(II) - ILU_W2C_VALS(C3 + JJ) * X(C1 + JJ)
                    END DO
                END DO
            END DO
            !DEC$ SIMD
            DO K = 0, 2
                X(C0 + K) = SUM_3(K)
            END DO
        END DO
        
        DO I = W_SIZE - 1, 0, -1
            C0 = 3 * ACTIVE_SIZE + 3 * I
            C1 = 9 * I
            !DEC$ SIMD
            DO K = 0, 8
                DIAG_VAL_33(K) = ILU_W_VALS(C1 + K)
            END DO
            CALL INV__(DIAG_VAL_33)
            !DEC$ SIMD 
            DO II = 0, 2
                X_3(II) = 0.D0
                C2 = 3 * II
                DO JJ = 0, 2
                    X_3(II) = X_3(II) + DIAG_VAL_33(C2 + JJ) * X(C0 + JJ)
                END DO
                Y(C0 + II) = X_3(II)
            END DO                
        END DO        
        !$OMP PARALLEL
        DO LEVEL = ILU_LEVEL_STARTS - 2, 0, -1
            LEVEL_SIZE = INDEPENDENT_ROW_STARTS(LEVEL + 1) - INDEPENDENT_ROW_STARTS(LEVEL)
            !IF (LEVEL_SIZE .LE. 10) THEN
            !    CALL OMP_SET_NUM_THREADS(1)
            !ELSE
            !    CALL OMP_SET_NUM_THREADS(THREADS)
            !END IF
            !$OMP DO SCHEDULE(STATIC) PRIVATE(LEVEL_INDEX, I, MIN_II, C0, K, SUM_3, J_INDEX, J, MIN_JJ, C1, VALUE_1, C2, II, C3, JJ, DIAG_VAL_33, X_3)
            DO LEVEL_INDEX = INDEPENDENT_ROW_STARTS(LEVEL + 1) - 1, INDEPENDENT_ROW_STARTS(LEVEL), -1
                I = INDEPENDENT_ROWS(LEVEL_INDEX)
                IF(IMPLICITNESS(NATURAL_ID(I)) .EQ. 1) THEN
                    MIN_II = 0
                ELSE
                    MIN_II = 2
                END IF
                C0 = 3 * I
                !DEC$ SIMD
                DO K = MIN_II, 2
                    SUM_3(K) = X(C0 + K)
                END DO                
                DO J_INDEX = DOF_WELL_ROW_STARTS(I + 1) - 1, DOF_WELL_ROW_STARTS(I), -1
                    J = DOF_WELL_CONNECTIVITY(J_INDEX)
                    C1 = 3 * ACTIVE_SIZE + 3 * J
                    C2 = 9 * J_INDEX
                    VALUE_1 = Y(3 * ACTIVE_SIZE + J)
                    !DEC$ SIMD
                    DO II = 0, 2
                        C3 = C2 + 3 * II
                        DO JJ = 0, 2
                            SUM_3(II) = SUM_3(II) - ILU_C2W_VALS(C3 + JJ) * Y(C1 + JJ)
                        END DO
                    END DO
                END DO                
                DO J_INDEX = ILU_PARALLEL_ROW_STARTS(LEVEL_INDEX + 1) - 1, ILU_DIAG_INDICES(LEVEL_INDEX) + 1, -1
                    J = ILU_PARALLEL_CONNECTIVITY(J_INDEX)                    
                    IF(IMPLICITNESS(NATURAL_ID(J)) .EQ. 1) THEN
                        MIN_JJ = 0
                    ELSE
                        MIN_JJ = 2
                    END IF
                    C1 = 3 * J
                    C2 = 9 * J_INDEX
                    !DEC$ SIMD
                    DO II = MIN_II, 2
                        C3 = C2 + 3 * II
                        DO JJ = MIN_JJ, 2
                            SUM_3(II) = SUM_3(II) - ILU_FLOW_VALS(C3 + JJ) * Y(C1 + JJ)                            
                        END DO
                    END DO
                END DO
                C2 = 9 * J_INDEX
                !DEC$ SIMD
                DO K = 0, 8
                    DIAG_VAL_33(K) = ILU_FLOW_VALS(C2 + K)
                END DO
                IF(IMPLICITNESS(NATURAL_ID(I)) .EQ. 1) THEN
                    CALL INV__(DIAG_VAL_33)
                ELSE
                    DIAG_VAL_33(8) = 1.D0/DIAG_VAL_33(8)
                END IF
                !DEC$ SIMD 
                DO II = MIN_II, 2
                    X_3(II) = 0.D0
                    C1 = 3 * II
                    DO JJ = MIN_II, 2
                        X_3(II) = X_3(II) + DIAG_VAL_33(C1 + JJ) * SUM_3(JJ)
                    END DO
                    Y(C0 + II) = X_3(II)
                END DO                
            END DO            
            !$OMP END DO            
        END DO
        !$OMP END PARALLEL
        !DEC$ SIMD
        DO I = 0, W_SIZE - 1
            C0 = 3 * ACTIVE_SIZE + 3 * I
            DO II = 0, 2
                X(C0 + II) = Y(C0 + II)
            END DO
        END DO
        
        !DEC$ SIMD
        DO I = 0, ACTIVE_SIZE - 1
            IF(IMPLICITNESS(NATURAL_ID(I)) .EQ. 1) THEN
                MIN_II = 0
            ELSE
                MIN_II = 2
            END IF    
            C1 = 3 * I
            DO K = MIN_II, 2
                X(C1 + K) = Y(C1 + K)
            END DO
        END DO
        !CALL OMP_SET_NUM_THREADS(THREADS)   
        
    END SUBROUTINE
    
    
    
                            !CONV_CV_SIZE, 
    SUBROUTINE ILU_SOLVE__(ACTIVE_SIZE, W_SIZE, CON_SIZE, ILU_PARALLEL_SIZE, ILU_LEVEL_STARTS, THREADS, X, B, Y, ILU_FLOW_VALS, ILU_C2W_VALS, ILU_W2C_VALS, ILU_W_VALS, ILU_PARALLEL_CONNECTIVITY, ILU_PARALLEL_ROW_STARTS, ILU_DIAG_INDICES, INDEPENDENT_ROWS, INDEPENDENT_ROW_STARTS, DOF_WELL_CONNECTIVITY, DOF_WELL_ROW_STARTS, WELL_DOF_CONNECTIVITY, WELL_DOF_ROW_STARTS)
        !DEC$ ATTRIBUTES DLLEXPORT, ALIAS : 'ILU_SOLVE__':: RETINA_SOLVER
        USE OMP_LIB        
        INTEGER, INTENT(IN) :: ACTIVE_SIZE, W_SIZE, CON_SIZE, ILU_PARALLEL_SIZE, ILU_LEVEL_STARTS, THREADS    !CONV_CV_SIZE, 
        
        DOUBLE PRECISION, DIMENSION(0:3 * ACTIVE_SIZE + 3 * W_SIZE - 1), INTENT(IN) :: B
        DOUBLE PRECISION, DIMENSION(0:3 * ACTIVE_SIZE + 3 * W_SIZE - 1), INTENT(INOUT) :: Y
        DOUBLE PRECISION, DIMENSION(0:3 * ACTIVE_SIZE + 3 * W_SIZE - 1), INTENT(INOUT) :: X
        DOUBLE PRECISION, DIMENSION(0:9 * ILU_PARALLEL_SIZE - 1), INTENT(IN) :: ILU_FLOW_VALS
        INTEGER, DIMENSION(0:ILU_PARALLEL_SIZE - 1), INTENT(IN) :: ILU_PARALLEL_CONNECTIVITY
        INTEGER, DIMENSION(0:ACTIVE_SIZE), INTENT(IN) :: ILU_PARALLEL_ROW_STARTS
        INTEGER, DIMENSION(0:ACTIVE_SIZE - 1), INTENT(IN) :: ILU_DIAG_INDICES
        DOUBLE PRECISION, DIMENSION(0:9 * CON_SIZE - 1), INTENT(IN) :: ILU_C2W_VALS, ILU_W2C_VALS
        INTEGER, DIMENSION(0:CON_SIZE - 1), INTENT(IN) :: DOF_WELL_CONNECTIVITY, WELL_DOF_CONNECTIVITY
        INTEGER, DIMENSION(0:ACTIVE_SIZE), INTENT(IN) :: DOF_WELL_ROW_STARTS
        INTEGER, DIMENSION(0:W_SIZE), INTENT(IN) :: WELL_DOF_ROW_STARTS
        
        DOUBLE PRECISION, DIMENSION(0:9 * W_SIZE - 1), INTENt(IN) :: ILU_W_VALS
        INTEGER, DIMENSION(0:ILU_LEVEL_STARTS - 1), INTENT(IN) :: INDEPENDENT_ROW_STARTS
        INTEGER, DIMENSION(0:ACTIVE_SIZE - 1), INTENT(IN) :: INDEPENDENT_ROWS
        
        INTEGER :: I, J, K, LEVEL_INDEX, LEVEL, LEVEL_SIZE
        DOUBLE PRECISION, DIMENSION(0:8) :: DIAG_VAL_33
        DOUBLE PRECISION, DIMENSION(0:2) :: SUM_3, X_3        
        DOUBLE PRECISION :: DIAG_VAL_1, SUM_1, VALUE_1, X_1
        INTEGER :: J_INDEX, INDEX
        INTEGER :: C0, C1, C2, C3, II, JJ
        INTEGER :: SIZE
        
        !DIR$ ASSUME_ALIGNED ILU_FLOW_VALS: 64
        !DIR$ ASSUME_ALIGNED ILU_PARALLEL_CONNECTIVITY: 64
        !DIR$ ASSUME_ALIGNED ILU_PARALLEL_ROW_STARTS: 64
        !DIR$ ASSUME_ALIGNED X: 64
        !DIR$ ASSUME_ALIGNED B: 64
        !DIR$ ASSUME_ALIGNED Y: 64        
        
        !$OMP PARALLEL
        DO LEVEL = 0, ILU_LEVEL_STARTS - 2
            LEVEL_SIZE = INDEPENDENT_ROW_STARTS(LEVEL + 1) - INDEPENDENT_ROW_STARTS(LEVEL)
            !IF (LEVEL_SIZE .LE. 10) THEN
            !    CALL OMP_SET_NUM_THREADS(1)
            !ELSE
            !    CALL OMP_SET_NUM_THREADS(THREADS)
            !END IF
            !$OMP DO SCHEDULE(STATIC) PRIVATE(LEVEL_INDEX, I, C0, K, SUM_3, J_INDEX, J, C1, C2, II, C3, JJ)
            DO LEVEL_INDEX = INDEPENDENT_ROW_STARTS(LEVEL), INDEPENDENT_ROW_STARTS(LEVEL + 1) - 1
                I = INDEPENDENT_ROWS(LEVEL_INDEX)
                C0 = 3 * I
                !DEC$ SIMD   
                DO K = 0, 2
                    SUM_3(K) = B(C0 + K)
                END DO                
                DO J_INDEX = ILU_PARALLEL_ROW_STARTS(LEVEL_INDEX), ILU_DIAG_INDICES(LEVEL_INDEX) - 1
                    J = ILU_PARALLEL_CONNECTIVITY(J_INDEX)                    
                    C1 = 3 * J
                    C2 = 9 * J_INDEX
                    !DEC$ SIMD
                    DO II = 0, 2
                        C3 = C2 + 3 * II
                        DO JJ = 0, 2
                            SUM_3(II) = SUM_3(II) - ILU_FLOW_VALS(C3 + JJ) * X(C1 + JJ)                             
                        END DO
                    END DO
                END DO   
                !DEC$ SIMD
                DO K = 0, 2
                    X(C0 + K) = SUM_3(K)
                END DO
            END DO
            !$OMP END DO        
        END DO        
        !$OMP END PARALLEL
        DO I = 0, W_SIZE - 1
            C0 = 3 * ACTIVE_SIZE + 3 * I
            !DEC$ SIMD
            DO K = 0, 2
                SUM_3(K) = B(C0 + K)
            END DO
            DO J_INDEX = WELL_DOF_ROW_STARTS(I), WELL_DOF_ROW_STARTS(I + 1) - 1
                J = WELL_DOF_CONNECTIVITY(J_INDEX)                
                C1 = 3 * J
                C2 = 9 * J_INDEX
                !DEC$ SIMD
                DO II = 0, 2
                    C3 = C2 + 3 * II
                    DO JJ = 0, 2
                        SUM_3(II) = SUM_3(II) - ILU_W2C_VALS(C3 + JJ) * X(C1 + JJ)
                    END DO
                END DO
            END DO
            !DEC$ SIMD
            DO K = 0, 2
                X(C0 + K) = SUM_3(K)
            END DO
        END DO
        
        DO I = W_SIZE - 1, 0, -1
            C0 = 3 * ACTIVE_SIZE + 3 * I
            C1 = 9 * I
            !DEC$ SIMD
            DO K = 0, 8
                DIAG_VAL_33(K) = ILU_W_VALS(C1 + K)
            END DO
            CALL INV__(DIAG_VAL_33)
            !DEC$ SIMD 
            DO II = 0, 2
                X_3(II) = 0.D0
                C2 = 3 * II
                DO JJ = 0, 2
                    X_3(II) = X_3(II) + DIAG_VAL_33(C2 + JJ) * X(C0 + JJ)
                END DO
                Y(C0 + II) = X_3(II)
            END DO                
        END DO        
        
        !$OMP PARALLEL
        DO LEVEL = ILU_LEVEL_STARTS - 2, 0, -1
            LEVEL_SIZE = INDEPENDENT_ROW_STARTS(LEVEL + 1) - INDEPENDENT_ROW_STARTS(LEVEL)
            !IF (LEVEL_SIZE .LE. 10) THEN
            !    CALL OMP_SET_NUM_THREADS(1)
            !ELSE
            !    CALL OMP_SET_NUM_THREADS(THREADS)
            !END IF
            !$OMP DO SCHEDULE(STATIC) PRIVATE(LEVEL_INDEX, I, C0, K, SUM_3, J_INDEX, J, C1, VALUE_1, C2, II, C3, JJ, DIAG_VAL_33, X_3)
            DO LEVEL_INDEX = INDEPENDENT_ROW_STARTS(LEVEL + 1) - 1, INDEPENDENT_ROW_STARTS(LEVEL), -1
                I = INDEPENDENT_ROWS(LEVEL_INDEX)
                C0 = 3 * I
                !DEC$ SIMD
                DO K = 0, 2
                    SUM_3(K) = X(C0 + K)
                END DO                
                DO J_INDEX = DOF_WELL_ROW_STARTS(I + 1) - 1, DOF_WELL_ROW_STARTS(I), -1
                    J = DOF_WELL_CONNECTIVITY(J_INDEX)
                    C1 = 3 * ACTIVE_SIZE + 3 * J
                    C2 = 9 * J_INDEX
                    VALUE_1 = Y(3 * ACTIVE_SIZE + J)
                    !DEC$ SIMD
                    DO II = 0, 2
                        C3 = C2 + 3 * II
                        DO JJ = 0, 2
                            SUM_3(II) = SUM_3(II) - ILU_C2W_VALS(C3 + JJ) * Y(C1 + JJ)
                        END DO
                    END DO
                END DO                
                DO J_INDEX = ILU_PARALLEL_ROW_STARTS(LEVEL_INDEX + 1) - 1, ILU_DIAG_INDICES(LEVEL_INDEX) + 1, -1
                    J = ILU_PARALLEL_CONNECTIVITY(J_INDEX)                    
                    C1 = 3 * J
                    C2 = 9 * J_INDEX
                    !DEC$ SIMD
                    DO II = 0, 2
                        C3 = C2 + 3 * II
                        DO JJ = 0, 2
                            SUM_3(II) = SUM_3(II) - ILU_FLOW_VALS(C3 + JJ) * Y(C1 + JJ)                            
                        END DO
                    END DO
                END DO
                C2 = 9 * J_INDEX
                !DEC$ SIMD
                DO K = 0, 8
                    DIAG_VAL_33(K) = ILU_FLOW_VALS(C2 + K)
                END DO
                CALL INV__(DIAG_VAL_33)
                !DEC$ SIMD 
                DO II = 0, 2
                    X_3(II) = 0.D0
                    C1 = 3 * II
                    DO JJ = 0, 2
                        X_3(II) = X_3(II) + DIAG_VAL_33(C1 + JJ) * SUM_3(JJ)
                    END DO
                    Y(C0 + II) = X_3(II)
                END DO                
            END DO            
            !$OMP END DO            
        END DO
        !$OMP END PARALLEL
        
        !DEC$ SIMD
        DO I = 0, W_SIZE - 1
            C0 = 3 * ACTIVE_SIZE + 3 * I
            DO II = 0, 2
                X(C0 + II) = Y(C0 + II)
            END DO
        END DO
        
        !DEC$ SIMD
        DO I = 0, ACTIVE_SIZE - 1
            C1 = 3 * I
            DO K = 0, 2
                X(C1 + K) = Y(C1 + K)
            END DO
        END DO
        !CALL OMP_SET_NUM_THREADS(THREADS)        
    END SUBROUTINE

    SUBROUTINE CONDENSE_ILU_APPROXIMATE_WELL_(CONV_CV_SIZE, ACTIVE_SIZE, W_SIZE, CON_SIZE, ILU_PARALLEL_SIZE, ILU_FLOW_VALS, MATRIX_C2W_VALS, MATRIX_W2C_VALS, MATRIX_W_TEMP, ILU_PARALLEL_CONNECTIVITY, ILU_PARALLEL_ROW_STARTS, DOF_WELL_CONNECTIVITY, DOF_WELL_ROW_STARTS, WELL_DOF_CONNECTIVITY, WELL_DOF_ROW_STARTS, INDEPENDENT_REVERSE_ROWS)
        !DEC$ ATTRIBUTES DLLEXPORT, ALIAS : 'CONDENSE_ILU_APPROXIMATE_WELL_':: RETINA_SOLVER    
        INTEGER, INTENT(IN) :: CONV_CV_SIZE, ACTIVE_SIZE, W_SIZE, CON_SIZE, ILU_PARALLEL_SIZE
        DOUBLE PRECISION, DIMENSION(0:9 * ILU_PARALLEL_SIZE - 1), INTENT(INOUT) :: ILU_FLOW_VALS
        DOUBLE PRECISION, DIMENSION(0:9 * CON_SIZE - 1), INTENT(IN) :: MATRIX_C2W_VALS, MATRIX_W2C_VALS
        DOUBLE PRECISION, DIMENSION(0:9 * W_SIZE - 1), INTENT(INOUT) :: MATRIX_W_TEMP
        INTEGER, DIMENSION(0:ILU_PARALLEL_SIZE - 1), INTENT(IN) :: ILU_PARALLEL_CONNECTIVITY
        INTEGER, DIMENSION(0:ACTIVE_SIZE), INTENT(IN) :: ILU_PARALLEL_ROW_STARTS
        INTEGER, DIMENSION(0:CON_SIZE - 1), INTENT(IN) :: DOF_WELL_CONNECTIVITY, WELL_DOF_CONNECTIVITY
        INTEGER, DIMENSION(0:CONV_CV_SIZE), INTENT(IN) :: DOF_WELL_ROW_STARTS
        INTEGER, DIMENSION(0:W_SIZE), INTENT(IN) :: WELL_DOF_ROW_STARTS
        INTEGER, DIMENSION(0:CONV_CV_SIZE - 1), INTENT(IN) :: INDEPENDENT_REVERSE_ROWS
        INTEGER :: I, J, J_INDEX, J_REVERSE_INDEX, INDEX, ILU_INDEX, GET_INDEX_
        
        ! ROW SUM
        !CALL ROW_SUM_(0, W_SIZE - 1, CON_SIZE, W_SIZE, WELL_DOF_ROW_STARTS, MATRIX_W2C_VALS, MATRIX_W_TEMP)
        !DO I = 0, W_SIZE - 1
        !    DO J_INDEX = WELL_DOF_ROW_STARTS(I), WELL_DOF_ROW_STARTS(I + 1) - 1
        !        J = WELL_DOF_CONNECTIVITY(J_INDEX)
        !        ! ASSUMING EACH CELL CONNECTS TO ONLY ONE WELL AT MOST
        !        INDEX = DOF_WELL_ROW_STARTS(J)    
        !        ILU_INDEX = GET_INDEX_(J, J, CONV_CV_SIZE, ILU_SIZE, ILU_CONNECTIVITY, ILU_ROW_STARTS)
        !        CALL MULT_MM__(MATRIX_C2W_VALS(9 * INDEX:9 * INDEX + 8), MATRIX_W_TEMP(9 * I:9 * I + 8), ILU_FLOW_VALS(9 * ILU_INDEX:9 * ILU_INDEX + 8), -1.D0)
        !    END DO
        !END DO
        
        ! COL SUM
        CALL ROW_SUM_TRANSPOSE_(0, W_SIZE - 1, 0, CONV_CV_SIZE - 1, CON_SIZE, W_SIZE, CONV_CV_SIZE, DOF_WELL_ROW_STARTS, DOF_WELL_CONNECTIVITY, MATRIX_C2W_VALS, MATRIX_W_TEMP)
        DO I = 0, W_SIZE - 1
            DO J_INDEX = WELL_DOF_ROW_STARTS(I), WELL_DOF_ROW_STARTS(I + 1) - 1
                J = WELL_DOF_CONNECTIVITY(J_INDEX)
                J_REVERSE_INDEX = INDEPENDENT_REVERSE_ROWS(J)
                ILU_INDEX = GET_INDEX_(J_REVERSE_INDEX, J, ACTIVE_SIZE, ILU_PARALLEL_SIZE, ILU_PARALLEL_CONNECTIVITY, ILU_PARALLEL_ROW_STARTS)
                CALL MULT_MTM_TRANSPOSE__(MATRIX_W2C_VALS(9 * J_INDEX:9 * J_INDEX + 8), MATRIX_W_TEMP(9 * I:9 * I + 8), ILU_FLOW_VALS(9 * ILU_INDEX:9 * ILU_INDEX + 8), -1.D0)                
            END DO
        END DO
    END SUBROUTINE
    