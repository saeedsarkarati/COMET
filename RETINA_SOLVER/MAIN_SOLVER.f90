
!****************************************************************************************************************!
!*                                                                                                              *!    
!*                                                  MAIN SOLVER                                                 *!
!*                                                                                                              *!    
!*                                                                                                              *!    
!****************************************************************************************************************!
    
    FUNCTION SOLVER_CONVERGED__(SIZE, X_RESULT, X_COPY) RESULT(ANS)
        !DEC$ ATTRIBUTES DLLEXPORT, ALIAS : 'SOLVER_CONVERGED__' :: RETINA_SOLVER
        USE CRS
        INTEGER, INTENT(IN) :: SIZE
        DOUBLE PRECISION, DIMENSION(0:SIZE - 1), INTENT(IN) :: X_RESULT, X_COPY
        LOGICAL :: ANS
        ANS = CRS_CONVERGED_(SIZE, X_RESULT, X_COPY, 1.D-6, 1.D-20)
    END FUNCTION
    
    FUNCTION CALC_VECTOR_NORM__(SIZE, V) RESULT(ANS)
        !DEC$ ATTRIBUTES DLLEXPORT, ALIAS : 'CALC_VECTOR_NORM__' :: RETINA_SOLVER
        USE CRS
        INTEGER, INTENT(IN) :: SIZE
        DOUBLE PRECISION, DIMENSION(0:SIZE - 1), INTENT(IN) :: V
        DOUBLE PRECISION :: ANS
        ANS = CRS_VECTOR_NORM_(SIZE, V)
    END FUNCTION
    
    
    FUNCTION DOT__(SIZE, V1, V2) RESULT(ANS)
        !DEC$ ATTRIBUTES DLLEXPORT, ALIAS : 'DOT__' :: RETINA_SOLVER
        USE CRS
        INTEGER, INTENT(IN) :: SIZE
        DOUBLE PRECISION, DIMENSION(0:SIZE - 1), INTENT(IN) :: V1, V2
        DOUBLE PRECISION :: ANS
        ANS = CRS_INNER_PRODUCT_(SIZE, V1, V2)        
    END FUNCTION
    
    SUBROUTINE SCALER_MULT__(SIZE, V, C)
        !DEC$ ATTRIBUTES DLLEXPORT, ALIAS : 'SCALER_MULT__' :: RETINA_SOLVER
        INTEGER, INTENT(IN) :: SIZE
        DOUBLE PRECISION, DIMENSION(0:SIZE - 1), INTENT(INOUT) :: V
        DOUBLE PRECISION :: C
        INTEGER :: I
        !DEC$ SIMD
        DO I = 0, SIZE - 1
            V(I) = V(I) * C
        END DO
    END SUBROUTINE
                                                          !, CONV_CV_SIZE                         !, ACTIVE_DOF      
    SUBROUTINE AIM_CREATE_R__(CV_SIZE, DOF_SHIFT, W_SHIFT, ACTIVE_SIZE, W_SIZE, CON_SIZE, DOF_SIZE, X_RESULT, VECTOR_VALS, R, MATRIX_FLOW_VALS, MATRIX_C2W_VALS, MATRIX_W2C_VALS, MATRIX_W_VALS, MATRIX_CONNECTIVITY, MATRIX_ROW_STARTS, DOF_WELL_CONNECTIVITY, DOF_WELL_ROW_STARTS, WELL_DOF_CONNECTIVITY, WELL_DOF_ROW_STARTS, IMPLICITNESS, NATURAL_ID)
        ! R = VECTOR_VALS - A * X_RESULT
        !DEC$ ATTRIBUTES DLLEXPORT, ALIAS : 'AIM_CREATE_R__':: RETINA_SOLVER
        USE OMP_LIB    
        INTEGER, INTENT(IN) :: CV_SIZE, DOF_SHIFT, W_SHIFT, ACTIVE_SIZE, W_SIZE, CON_SIZE, DOF_SIZE       !, CONV_CV_SIZE
        !INTEGER, DIMENSION(0:ACTIVE_SIZE - 1), INTENT(IN) :: ACTIVE_DOF
        DOUBLE PRECISION, DIMENSION(0:ACTIVE_SIZE * 3 + 3 * W_SIZE - 1), INTENT(IN) :: X_RESULT, VECTOR_VALS
        DOUBLE PRECISION, DIMENSION(0:ACTIVE_SIZE * 3 + 3 * W_SIZE - 1), INTENT(INOUT) :: R
        DOUBLE PRECISION, DIMENSION(0:9 * DOF_SIZE - 1), INTENT(IN) :: MATRIX_FLOW_VALS
        DOUBLE PRECISION, DIMENSION(0:9 * CON_SIZE - 1), INTENT(IN) :: MATRIX_C2W_VALS, MATRIX_W2C_VALS
        INTEGER, DIMENSION(0:CV_SIZE - 1), INTENT(IN) :: IMPLICITNESS, NATURAL_ID
        DOUBLE PRECISION, DIMENSION(0:9 * W_SIZE - 1), INTENT(IN) :: MATRIX_W_VALS
        INTEGER, DIMENSION(0:DOF_SIZE - 1), INTENT(IN) :: MATRIX_CONNECTIVITY
        INTEGER, DIMENSION(0:ACTIVE_SIZE), INTENT(IN) :: MATRIX_ROW_STARTS
        INTEGER, DIMENSION(0:CON_SIZE - 1), INTENT(IN) :: DOF_WELL_CONNECTIVITY, WELL_DOF_CONNECTIVITY
        INTEGER, DIMENSION(0:ACTIVE_SIZE), INTENT(IN) :: DOF_WELL_ROW_STARTS
        INTEGER, DIMENSION(0:W_SIZE), INTENT(IN) :: WELL_DOF_ROW_STARTS
        !CALL MKL_CSPBLAS_DBSRGEMV('N', CONV_CV_SIZE, 3, MATRIX_FLOW_VALS, MATRIX_ROW_STARTS, MATRIX_CONNECTIVITY, X_RESULT, VECTOR_VALS)
    
        INTEGER :: I, AI, INDEX, COLUMN, II, JJ, C0, C1, C2, C3, MIN_II, MIN_JJ
                                                          !, CONV_CV_SIZE    !, ACTIVE_DOF
        !$OMP PARALLEL SHARED(CV_SIZE, DOF_SHIFT, W_SHIFT, ACTIVE_SIZE, W_SIZE, X_RESULT, VECTOR_VALS, R, MATRIX_FLOW_VALS, MATRIX_C2W_VALS, MATRIX_CONNECTIVITY, MATRIX_ROW_STARTS, DOF_WELL_CONNECTIVITY, DOF_WELL_ROW_STARTS)
            !$OMP DO SCHEDULE(STATIC) PRIVATE(AI, I, MIN_II, C0, II, INDEX, COLUMN, MIN_JJ, C1, C2, C3, JJ)
        DO AI = 0, ACTIVE_SIZE - 1
            I = AI  !ACTIVE_DOF(AI)
            IF(IMPLICITNESS(NATURAL_ID(I)) .EQ. 1) THEN
                MIN_II = 0
            ELSE
                MIN_II = 2
            END IF
            C0 = 3 * I + DOF_SHIFT
            !DEC$ SIMD
            DO II = MIN_II, 2
                R(C0 + II) = VECTOR_VALS(C0 + II)
            END DO
            DO INDEX = MATRIX_ROW_STARTS(I), MATRIX_ROW_STARTS(I + 1) - 1
                COLUMN = MATRIX_CONNECTIVITY(INDEX)
                IF(IMPLICITNESS(NATURAL_ID(COLUMN)) .EQ. 1) THEN
                    MIN_JJ = 0
                ELSE
                    MIN_JJ = 2
                END IF
                C1 = 3 * COLUMN + DOF_SHIFT
                C2 = 9 * INDEX
                !CALL MM_PREFETCH(MATRIX_FLOW_VALS(C2 + 6400), 1)
                !DEC$ SIMD            
                DO II = MIN_II, 2
                    C3 = C2 + 3 * II
                    DO JJ = MIN_JJ, 2
                        R(C0 + II) = R(C0 + II) - MATRIX_FLOW_VALS(C3 + JJ) * X_RESULT(C1 + JJ)
                    END DO
                END DO
            END DO
            DO INDEX = DOF_WELL_ROW_STARTS(I), DOF_WELL_ROW_STARTS(I + 1) - 1
                COLUMN = DOF_WELL_CONNECTIVITY(INDEX)
                C1 = 3 * COLUMN + W_SHIFT
                C2 = 9 * INDEX
                !DEC$ SIMD            
                DO II = 0, 2
                    C3 = C2 + 3 * II
                    DO JJ = 0, 2
                        R(C0 + II) = R(C0 + II) - MATRIX_C2W_VALS(C3 + JJ) * X_RESULT(C1 + JJ)
                    END DO
                END DO
            END DO
        END DO
            !$OMP END DO
        !$OMP END PARALLEL
                                                        !, CONV_CV_SIZE                                
        !$OMP PARALLEL SHARED(DOF_SHIFT, W_SHIFT, W_SIZE, X_RESULT, VECTOR_VALS, R, MATRIX_W_VALS, MATRIX_W2C_VALS, WELL_DOF_CONNECTIVITY, WELL_DOF_ROW_STARTS)
            !$OMP DO SCHEDULE(STATIC) PRIVATE(I, C0, II, INDEX, COLUMN, C1, C2, C3, JJ)
        DO I = 0, W_SIZE - 1
            C0 = 3 * I + W_SHIFT
            !DEC$ SIMD
            DO II = 0, 2
                R(C0 + II) = VECTOR_VALS(C0 + II)
            END DO
            DO INDEX = WELL_DOF_ROW_STARTS(I), WELL_DOF_ROW_STARTS(I + 1) - 1               
                COLUMN = WELL_DOF_CONNECTIVITY(INDEX)
                C1 = 3 * COLUMN + DOF_SHIFT
                C2 = 9 * INDEX
                !DEC$ SIMD            
                DO II = 0, 2
                    C3 = C2 + 3 * II
                    DO JJ = 0, 2
                        R(C0 + II) = R(C0 + II) - MATRIX_W2C_VALS(C3 + JJ) * X_RESULT(C1 + JJ)
                    END DO
                END DO
            END DO
            C2 = 9 * I
            !DEC$ SIMD
            DO II = 0, 2
                C3 = C2 + 3 * II
                DO JJ = 0, 2
                    R(C0 + II) = R(C0 + II) - MATRIX_W_VALS(C3 + JJ) * X_RESULT(C0 + JJ)
                END DO
            END DO
        END DO
            !$OMP END DO
        !$OMP END PARALLEL        
    END SUBROUTINE
    
                                              !, CONV_CV_SIZE                                                      !, ACTIVE_DOF
    SUBROUTINE CREATE_R__(DOF_SHIFT, W_SHIFT, ACTIVE_SIZE, W_SIZE, CON_SIZE, DOF_SIZE, X_RESULT, VECTOR_VALS, R, MATRIX_FLOW_VALS, MATRIX_C2W_VALS, MATRIX_W2C_VALS, MATRIX_W_VALS, MATRIX_CONNECTIVITY, MATRIX_ROW_STARTS, DOF_WELL_CONNECTIVITY, DOF_WELL_ROW_STARTS, WELL_DOF_CONNECTIVITY, WELL_DOF_ROW_STARTS)
        ! R = VECTOR_VALS - A * X_RESULT
        !DEC$ ATTRIBUTES DLLEXPORT, ALIAS : 'CREATE_R__':: RETINA_SOLVER
        USE OMP_LIB    
        INTEGER, INTENT(IN) :: DOF_SHIFT, W_SHIFT, ACTIVE_SIZE, W_SIZE, CON_SIZE, DOF_SIZE    !, CONV_CV_SIZE
        !INTEGER, DIMENSION(0:ACTIVE_SIZE - 1), INTENT(IN) :: ACTIVE_DOF
        DOUBLE PRECISION, DIMENSION(0:ACTIVE_SIZE * 3 + 3 * W_SIZE - 1), INTENT(IN) :: X_RESULT, VECTOR_VALS
        DOUBLE PRECISION, DIMENSION(0:ACTIVE_SIZE * 3 + 3 * W_SIZE - 1), INTENT(INOUT) :: R
        DOUBLE PRECISION, DIMENSION(0:9 * DOF_SIZE - 1), INTENT(IN) :: MATRIX_FLOW_VALS
        DOUBLE PRECISION, DIMENSION(0:9 * CON_sIZE - 1), INTENT(IN) :: MATRIX_C2W_VALS, MATRIX_W2C_VALS
        DOUBLE PRECISION, DIMENSION(0:9 * W_SIZE - 1), INTENT(IN) :: MATRIX_W_VALS
        INTEGER, DIMENSION(0:DOF_SIZE - 1), INTENT(IN) :: MATRIX_CONNECTIVITY
        INTEGER, DIMENSION(0:ACTIVE_SIZE), INTENT(IN) :: MATRIX_ROW_STARTS
        INTEGER, DIMENSION(0:CON_SIZE - 1), INTENT(IN) :: DOF_WELL_CONNECTIVITY, WELL_DOF_CONNECTIVITY
        INTEGER, DIMENSION(0:ACTIVE_SIZE), INTENT(IN) :: DOF_WELL_ROW_STARTS
        INTEGER, DIMENSION(0:W_SIZE), INTENT(IN) :: WELL_DOF_ROW_STARTS
        !CALL MKL_CSPBLAS_DBSRGEMV('N', CONV_CV_SIZE, 3, MATRIX_FLOW_VALS, MATRIX_ROW_STARTS, MATRIX_CONNECTIVITY, X_RESULT, VECTOR_VALS)
    
        INTEGER :: I, AI, INDEX, COLUMN, II, JJ, C0, C1, C2, C3
                                                !, CONV_CV_SIZE                                   !, ACTIVE_DOF                 
        !$OMP PARALLEL SHARED(DOF_SHIFT, W_SHIFT, ACTIVE_SIZE, W_SIZE, X_RESULT, VECTOR_VALS, R, MATRIX_FLOW_VALS, MATRIX_C2W_VALS, MATRIX_CONNECTIVITY, MATRIX_ROW_STARTS, DOF_WELL_CONNECTIVITY, DOF_WELL_ROW_STARTS)
            !$OMP DO SCHEDULE(STATIC) PRIVATE(AI, I, C0, II, INDEX, COLUMN, C1, C2, C3, JJ)
        DO AI = 0, ACTIVE_SIZE - 1
            I = AI  !ACTIVE_DOF(AI)
            C0 = 3 * I + DOF_SHIFT
            !DEC$ SIMD
            DO II = 0, 2
                R(C0 + II) = VECTOR_VALS(C0 + II)
            END DO
            DO INDEX = MATRIX_ROW_STARTS(I), MATRIX_ROW_STARTS(I + 1) - 1
                COLUMN = MATRIX_CONNECTIVITY(INDEX)
                C1 = 3 * COLUMN + DOF_SHIFT
                C2 = 9 * INDEX
                !CALL MM_PREFETCH(MATRIX_FLOW_VALS(C2 + 6400), 1)
                !DEC$ SIMD            
                DO II = 0, 2
                    C3 = C2 + 3 * II
                    DO JJ = 0, 2
                        R(C0 + II) = R(C0 + II) - MATRIX_FLOW_VALS(C3 + JJ) * X_RESULT(C1 + JJ)
                    END DO
                END DO
            END DO
            DO INDEX = DOF_WELL_ROW_STARTS(I), DOF_WELL_ROW_STARTS(I + 1) - 1
                COLUMN = DOF_WELL_CONNECTIVITY(INDEX)
                C1 = 3 * COLUMN + W_SHIFT
                C2 = 9 * INDEX
                !DEC$ SIMD            
                DO II = 0, 2
                    C3 = C2 + 3 * II
                    DO JJ = 0, 2
                        R(C0 + II) = R(C0 + II) - MATRIX_C2W_VALS(C3 + JJ) * X_RESULT(C1 + JJ)
                    END DO
                END DO
            END DO
        END DO
            !$OMP END DO
        !$OMP END PARALLEL
                                                        !, CONV_CV_SIZE                    
        !$OMP PARALLEL SHARED(DOF_SHIFT, W_SHIFT, W_SIZE, X_RESULT, VECTOR_VALS, R, MATRIX_W_VALS, MATRIX_W2C_VALS, WELL_DOF_CONNECTIVITY, WELL_DOF_ROW_STARTS)
            !$OMP DO SCHEDULE(STATIC) PRIVATE(I, C0, II, INDEX, COLUMN, C1, C2, C3, JJ)
        DO I = 0, W_SIZE - 1
            C0 = 3 * I + W_SHIFT
            !DEC$ SIMD
            DO II = 0, 2
                R(C0 + II) = VECTOR_VALS(C0 + II)
            END DO
            DO INDEX = WELL_DOF_ROW_STARTS(I), WELL_DOF_ROW_STARTS(I + 1) - 1               
                COLUMN = WELL_DOF_CONNECTIVITY(INDEX)
                C1 = 3 * COLUMN + DOF_SHIFT
                C2 = 9 * INDEX
                !DEC$ SIMD            
                DO II = 0, 2
                    C3 = C2 + 3 * II
                    DO JJ = 0, 2
                        R(C0 + II) = R(C0 + II) - MATRIX_W2C_VALS(C3 + JJ) * X_RESULT(C1 + JJ)
                    END DO
                END DO
            END DO
            C2 = 9 * I
            !DEC$ SIMD
            DO II = 0, 2
                C3 = C2 + 3 * II
                DO JJ = 0, 2
                    R(C0 + II) = R(C0 + II) - MATRIX_W_VALS(C3 + JJ) * X_RESULT(C0 + JJ)
                END DO
            END DO
        END DO
            !$OMP END DO
        !$OMP END PARALLEL        
    END SUBROUTINE
                                                !, CONV_CV_SIZE                          !, ACTIVE_DOF   
    SUBROUTINE CREATE_CPR_R__(DOF_SHIFT, W_SHIFT, ACTIVE_SIZE, W_SIZE, CON_SIZE, DOF_SIZE, X_RESULT, VECTOR_VALS, R, MATRIX_FLOW_VALS, MATRIX_C2W_VALS, MATRIX_W2C_VALS, MATRIX_W_VALS, MATRIX_CONNECTIVITY, MATRIX_ROW_STARTS, DOF_WELL_CONNECTIVITY, DOF_WELL_ROW_STARTS, WELL_DOF_CONNECTIVITY, WELL_DOF_ROW_STARTS)
        ! R = VECTOR_VALS - A * X_RESULT    
        !DEC$ ATTRIBUTES DLLEXPORT, ALIAS : 'CREATE_CPR_R__':: RETINA_SOLVER
        USE OMP_LIB    
        INTEGER, INTENT(IN) :: DOF_SHIFT, W_SHIFT, ACTIVE_SIZE, W_SIZE, CON_SIZE, DOF_SIZE    !, CONV_CV_SIZE
        !INTEGER, DIMENSION(0:ACTIVE_SIZE - 1), INTENT(IN) :: ACTIVE_DOF
        DOUBLE PRECISION, DIMENSION(0:ACTIVE_SIZE * 3 + 3 * W_SIZE - 1), INTENT(IN) :: X_RESULT, VECTOR_VALS
        DOUBLE PRECISION, DIMENSION(0:ACTIVE_SIZE * 3 + 3 * W_SIZE - 1), INTENT(INOUT) :: R
        DOUBLE PRECISION, DIMENSION(0:9 * DOF_SIZE - 1), INTENT(IN) :: MATRIX_FLOW_VALS
        DOUBLE PRECISION, DIMENSION(0:9 * CON_sIZE - 1), INTENT(IN) :: MATRIX_C2W_VALS, MATRIX_W2C_VALS
        DOUBLE PRECISION, DIMENSION(0:9 * W_SIZE - 1), INTENT(IN) :: MATRIX_W_VALS
        INTEGER, DIMENSION(0:DOF_SIZE - 1), INTENT(IN) :: MATRIX_CONNECTIVITY
        INTEGER, DIMENSION(0:ACTIVE_SIZE), INTENT(IN) :: MATRIX_ROW_STARTS
        INTEGER, DIMENSION(0:CON_SIZE - 1), INTENT(IN) :: DOF_WELL_CONNECTIVITY, WELL_DOF_CONNECTIVITY
        INTEGER, DIMENSION(0:ACTIVE_SIZE), INTENT(IN) :: DOF_WELL_ROW_STARTS
        INTEGER, DIMENSION(0:W_SIZE), INTENT(IN) :: WELL_DOF_ROW_STARTS
        !CALL MKL_CSPBLAS_DBSRGEMV('N', CONV_CV_SIZE, 3, MATRIX_FLOW_VALS, MATRIX_ROW_STARTS, MATRIX_CONNECTIVITY, X_RESULT, VECTOR_VALS)
    
        INTEGER :: I, AI, INDEX, COLUMN, II, C0, C1, C2, C3
                                                !, CONV_CV_SIZE      !, ACTIVE_DOF       
        !$OMP PARALLEL SHARED(DOF_SHIFT, W_SHIFT, ACTIVE_SIZE, W_SIZE, X_RESULT, VECTOR_VALS, R, MATRIX_FLOW_VALS, MATRIX_C2W_VALS, MATRIX_CONNECTIVITY, MATRIX_ROW_STARTS, DOF_WELL_CONNECTIVITY, DOF_WELL_ROW_STARTS)
            !$OMP DO SCHEDULE(STATIC) PRIVATE(AI, I, C0, II, INDEX, COLUMN, C1, C2, C3)
        DO AI = 0, ACTIVE_SIZE - 1
            I = AI         !ACTIVE_DOF(AI)
            C0 = 3 * I + DOF_SHIFT
            !DEC$ SIMD
            DO II = 0, 2
                R(C0 + II) = VECTOR_VALS(C0 + II)
            END DO
            DO INDEX = MATRIX_ROW_STARTS(I), MATRIX_ROW_STARTS(I + 1) - 1
                COLUMN = MATRIX_CONNECTIVITY(INDEX)
                C1 = 3 * COLUMN + DOF_SHIFT
                C2 = 9 * INDEX
                !CALL MM_PREFETCH(MATRIX_FLOW_VALS(C2 + 6400), 1)
                !DEC$ SIMD            
                DO II = 0, 2
                    C3 = C2 + 3 * II
                    R(C0 + II) = R(C0 + II) - MATRIX_FLOW_VALS(C3 + 2) * X_RESULT(C1 + 2)
                END DO
            END DO
            DO INDEX = DOF_WELL_ROW_STARTS(I), DOF_WELL_ROW_STARTS(I + 1) - 1
                COLUMN = DOF_WELL_CONNECTIVITY(INDEX)
                C1 = 3 * COLUMN + W_SHIFT
                C2 = 9 * INDEX
                !DEC$ SIMD            
                DO II = 0, 2
                    C3 = C2 + 3 * II
                    R(C0 + II) = R(C0 + II) - MATRIX_C2W_VALS(C3 + 2) * X_RESULT(C1 + 2)
                END DO
            END DO
        END DO
            !$OMP END DO NOWAIT
        !$OMP END PARALLEL
                                                        !, CONV_CV_SIZE                
        !$OMP PARALLEL SHARED(DOF_SHIFT, W_SHIFT, W_SIZE, X_RESULT, VECTOR_VALS, R, MATRIX_W_VALS, MATRIX_W2C_VALS, WELL_DOF_CONNECTIVITY, WELL_DOF_ROW_STARTS)
            !$OMP DO SCHEDULE(STATIC) PRIVATE(I, C0, II, INDEX, COLUMN, C1, C2, C3)
        DO I = 0, W_SIZE - 1
            C0 = 3 * I + W_SHIFT
            !DEC$ SIMD
            DO II = 0, 2
                R(C0 + II) = VECTOR_VALS(C0 + II)
            END DO
            DO INDEX = WELL_DOF_ROW_STARTS(I), WELL_DOF_ROW_STARTS(I + 1) - 1               
                COLUMN = WELL_DOF_CONNECTIVITY(INDEX)
                C1 = 3 * COLUMN + DOF_SHIFT
                C2 = 9 * INDEX
                !DEC$ SIMD            
                DO II = 0, 2
                    C3 = C2 + 3 * II
                    R(C0 + II) = R(C0 + II) - MATRIX_W2C_VALS(C3 + 2) * X_RESULT(C1 + 2)
                END DO
            END DO
            C2 = 9 * I
            !DEC$ SIMD
            DO II = 0, 2
                C3 = C2 + 3 * II
                R(C0 + II) = R(C0 + II) - MATRIX_W_VALS(C3 + 2) * X_RESULT(C0 + 2)
            END DO
        END DO
            !$OMP END DO NOWAIT
        !$OMP END PARALLEL        
    END SUBROUTINE
                                                             !, CONV_CV_SIZE                           !, ACTIVE_DOF   
    SUBROUTINE AIM_MAIN_MULT_MV__(CV_SIZE, DOF_SHIFT, W_SHIFT, ACTIVE_SIZE, W_SIZE, CON_SIZE, DOF_SIZE, SOURCE, DEST, MATRIX_FLOW_VALS, MATRIX_C2W_VALS, MATRIX_W2C_VALS, MATRIX_W_VALS, MATRIX_CONNECTIVITY, MATRIX_ROW_STARTS, DOF_WELL_CONNECTIVITY, DOF_WELL_ROW_STARTS, WELL_DOF_CONNECTIVITY, WELL_DOF_ROW_STARTS, IMPLICITNESS, NATURAL_ID)
        ! DEST = A * SOURCE
        !DEC$ ATTRIBUTES DLLEXPORT, ALIAS : 'AIM_MAIN_MULT_MV__':: RETINA_SOLVER
        USE OMP_LIB    
        INTEGER, INTENT(IN) :: CV_SIZE, DOF_SHIFT, W_SHIFT, ACTIVE_SIZE, W_SIZE, CON_SIZE, DOF_SIZE   !, CONV_CV_SIZE
        !INTEGER, DIMENSION(0:ACTIVE_SIZE - 1), INTENT(IN) :: ACTIVE_DOF
        DOUBLE PRECISION, DIMENSION(0:ACTIVE_SIZE * 3 + 3 * W_SIZE - 1), INTENT(IN) :: SOURCE
        DOUBLE PRECISION, DIMENSION(0:ACTIVE_SIZE * 3 + 3 * W_SIZE - 1), INTENT(INOUT) :: DEST
        DOUBLE PRECISION, DIMENSION(0:9 * DOF_SIZE - 1), INTENT(IN) :: MATRIX_FLOW_VALS
        DOUBLE PRECISION, DIMENSION(0:9 * CON_sIZE - 1), INTENT(IN) :: MATRIX_C2W_VALS, MATRIX_W2C_VALS
        DOUBLE PRECISION, DIMENSION(0:9 * W_SIZE - 1), INTENT(IN) :: MATRIX_W_VALS
        INTEGER, DIMENSION(0:DOF_SIZE - 1), INTENT(IN) :: MATRIX_CONNECTIVITY
        INTEGER, DIMENSION(0:CV_SIZE - 1), INTENT(IN) :: IMPLICITNESS, NATURAL_ID
        INTEGER, DIMENSION(0:ACTIVE_SIZE), INTENT(IN) :: MATRIX_ROW_STARTS
        INTEGER, DIMENSION(0:CON_SIZE - 1), INTENT(IN) :: DOF_WELL_CONNECTIVITY, WELL_DOF_CONNECTIVITY
        INTEGER, DIMENSION(0:ACTIVE_SIZE), INTENT(IN) :: DOF_WELL_ROW_STARTS
        INTEGER, DIMENSION(0:W_SIZE), INTENT(IN) :: WELL_DOF_ROW_STARTS
        !CALL MKL_CSPBLAS_DBSRGEMV('N', CONV_CV_SIZE, 3, MATRIX_FLOW_VALS, MATRIX_ROW_STARTS, MATRIX_CONNECTIVITY, X_RESULT, VECTOR_VALS)
    
        INTEGER :: AI, I, INDEX, COLUMN, II, JJ, C0, C1, C2, C3, MIN_II, MIN_JJ
                                                        !, CONV_CV_SIZE       !, ACTIVE_DOF 
        !$OMP PARALLEL SHARED(CV_SIZE, DOF_SHIFT, W_SHIFT, ACTIVE_SIZE, W_SIZE, SOURCE, DEST, MATRIX_FLOW_VALS, MATRIX_C2W_VALS, MATRIX_CONNECTIVITY, MATRIX_ROW_STARTS, DOF_WELL_CONNECTIVITY, DOF_WELL_ROW_STARTS)
            !$OMP DO SCHEDULE(STATIC) PRIVATE(AI, I, MIN_II, C0, II, INDEX, COLUMN, MIN_JJ, C1, C2, C3, JJ)
        DO AI = 0, ACTIVE_SIZE - 1
            I = AI      !ACTIVE_DOF(AI)
            IF(IMPLICITNESS(NATURAL_ID(I)) .EQ. 1) THEN
                MIN_II = 0
            ELSE
                MIN_II = 2
            END IF
            C0 = 3 * I + DOF_SHIFT
            !DEC$ SIMD
            DO II = MIN_II, 2
                DEST(C0 + II) = 0.D0
            END DO
            DO INDEX = MATRIX_ROW_STARTS(I), MATRIX_ROW_STARTS(I + 1) - 1
                COLUMN = MATRIX_CONNECTIVITY(INDEX)
                IF(IMPLICITNESS(NATURAL_ID(COLUMN)) .EQ. 1) THEN
                    MIN_JJ = 0
                ELSE
                    MIN_JJ = 2
                END IF
                C1 = 3 * COLUMN + DOF_SHIFT
                C2 = 9 * INDEX
                !CALL MM_PREFETCH(MATRIX_FLOW_VALS(C2 + 6400), 1)
                !DEC$ SIMD            
                DO II = MIN_II, 2
                    C3 = C2 + 3 * II
                    DO JJ = MIN_JJ, 2
                        DEST(C0 + II) = DEST(C0 + II) + MATRIX_FLOW_VALS(C3 + JJ) * SOURCE(C1 + JJ)
                    END DO
                END DO
            END DO
            DO INDEX = DOF_WELL_ROW_STARTS(I), DOF_WELL_ROW_STARTS(I + 1) - 1
                COLUMN = DOF_WELL_CONNECTIVITY(INDEX)
                C1 = 3 * COLUMN + W_SHIFT
                C2 = 9 * INDEX
                !DEC$ SIMD            
                DO II = 0, 2
                    C3 = C2 + 3 * II
                    DO JJ = 0, 2
                        DEST(C0 + II) = DEST(C0 + II) + MATRIX_C2W_VALS(C3 + JJ) * SOURCE(C1 + JJ)
                    END DO
                END DO
            END DO
        END DO
            !$OMP END DO NOWAIT
        !$OMP END PARALLEL
                                                        ! , CONV_CV_SIZE       
        !$OMP PARALLEL SHARED(DOF_SHIFT, W_SHIFT, W_SIZE, SOURCE, DEST, MATRIX_W_VALS, MATRIX_W2C_VALS, WELL_DOF_CONNECTIVITY, WELL_DOF_ROW_STARTS)
            !$OMP DO SCHEDULE(STATIC) PRIVATE(I, C0, II, INDEX, COLUMN, C1, C2, C3, JJ)
        DO I = 0, W_SIZE - 1
            C0 = 3 * I + W_SHIFT
            !DEC$ SIMD
            DO II = 0, 2
                DEST(C0 + II) = 0.D0
            END DO
            DO INDEX = WELL_DOF_ROW_STARTS(I), WELL_DOF_ROW_STARTS(I + 1) - 1
                COLUMN = WELL_DOF_CONNECTIVITY(INDEX)
                C1 = 3 * COLUMN + DOF_SHIFT
                C2 = 9 * INDEX
                !DEC$ SIMD            
                DO II = 0, 2
                    C3 = C2 + 3 * II
                    DO JJ = 0, 2
                        DEST(C0 + II) = DEST(C0 + II) + MATRIX_W2C_VALS(C3 + JJ) * SOURCE(C1 + JJ)
                    END DO
                END DO
            END DO
            C2 = 9 * I
            !DEC$ SIMD
            DO II = 0, 2
                C3 = C2 + 3 * II
                DO JJ = 0, 2
                    DEST(C0 + II) = DEST(C0 + II) + MATRIX_W_VALS(C3 + JJ) * SOURCE(C0 + JJ)
                END DO
            END DO
        END DO
            !$OMP END DO NOWAIT
        !$OMP END PARALLEL
        
        END SUBROUTINE

                                                !, CONV_CV_SIZE                                                        !, ACTIVE_DOF
    SUBROUTINE MAIN_MULT_MV__(DOF_SHIFT, W_SHIFT, ACTIVE_SIZE, W_SIZE, CON_SIZE, DOF_SIZE, SOURCE, DEST, MATRIX_FLOW_VALS, MATRIX_C2W_VALS, MATRIX_W2C_VALS, MATRIX_W_VALS, MATRIX_CONNECTIVITY, MATRIX_ROW_STARTS, DOF_WELL_CONNECTIVITY, DOF_WELL_ROW_STARTS, WELL_DOF_CONNECTIVITY, WELL_DOF_ROW_STARTS)
        ! DEST = A * SOURCE
        !DEC$ ATTRIBUTES DLLEXPORT, ALIAS : 'MAIN_MULT_MV__':: RETINA_SOLVER
        USE OMP_LIB    
        INTEGER, INTENT(IN) :: DOF_SHIFT, W_SHIFT, ACTIVE_SIZE, W_SIZE, CON_SIZE, DOF_SIZE        !, CONV_CV_SIZE
        !INTEGER, DIMENSION(0:ACTIVE_SIZE - 1), INTENT(IN) :: ACTIVE_DOF
        DOUBLE PRECISION, DIMENSION(0:ACTIVE_SIZE * 3 + 3 * W_SIZE - 1), INTENT(IN) :: SOURCE
        DOUBLE PRECISION, DIMENSION(0:ACTIVE_SIZE * 3 + 3 * W_SIZE - 1), INTENT(INOUT) :: DEST
        DOUBLE PRECISION, DIMENSION(0:9 * DOF_SIZE - 1), INTENT(IN) :: MATRIX_FLOW_VALS
        DOUBLE PRECISION, DIMENSION(0:9 * CON_sIZE - 1), INTENT(IN) :: MATRIX_C2W_VALS, MATRIX_W2C_VALS
        DOUBLE PRECISION, DIMENSION(0:9 * W_SIZE - 1), INTENT(IN) :: MATRIX_W_VALS
        INTEGER, DIMENSION(0:DOF_SIZE - 1), INTENT(IN) :: MATRIX_CONNECTIVITY
        INTEGER, DIMENSION(0:ACTIVE_SIZE), INTENT(IN) :: MATRIX_ROW_STARTS
        INTEGER, DIMENSION(0:CON_SIZE - 1), INTENT(IN) :: DOF_WELL_CONNECTIVITY, WELL_DOF_CONNECTIVITY
        INTEGER, DIMENSION(0:ACTIVE_SIZE), INTENT(IN) :: DOF_WELL_ROW_STARTS
        INTEGER, DIMENSION(0:W_SIZE), INTENT(IN) :: WELL_DOF_ROW_STARTS
        !CALL MKL_CSPBLAS_DBSRGEMV('N', CONV_CV_SIZE, 3, MATRIX_FLOW_VALS, MATRIX_ROW_STARTS, MATRIX_CONNECTIVITY, X_RESULT, VECTOR_VALS)
    
        INTEGER :: AI, I, INDEX, COLUMN, II, JJ, C0, C1, C2, C3
                                                !, CONV_CV_SIZE                                    !, ACTIVE_DOF
        !$OMP PARALLEL SHARED(DOF_SHIFT, W_SHIFT, ACTIVE_SIZE, W_SIZE, SOURCE, DEST, MATRIX_FLOW_VALS, MATRIX_C2W_VALS, MATRIX_CONNECTIVITY, MATRIX_ROW_STARTS, DOF_WELL_CONNECTIVITY, DOF_WELL_ROW_STARTS)
            !$OMP DO SCHEDULE(STATIC) PRIVATE(AI, I, C0, II, INDEX, COLUMN, C1, C2, C3, JJ)
        DO AI = 0, ACTIVE_SIZE - 1
            I = AI      !ACTIVE_DOF(AI)
            C0 = 3 * I + DOF_SHIFT
            !DEC$ SIMD
            DO II = 0, 2
                DEST(C0 + II) = 0.D0
            END DO
            DO INDEX = MATRIX_ROW_STARTS(I), MATRIX_ROW_STARTS(I + 1) - 1
                COLUMN = MATRIX_CONNECTIVITY(INDEX)
                C1 = 3 * COLUMN + DOF_SHIFT
                C2 = 9 * INDEX
                !CALL MM_PREFETCH(MATRIX_FLOW_VALS(C2 + 6400), 1)
                !DEC$ SIMD            
                DO II = 0, 2
                    C3 = C2 + 3 * II
                    DO JJ = 0, 2
                        DEST(C0 + II) = DEST(C0 + II) + MATRIX_FLOW_VALS(C3 + JJ) * SOURCE(C1 + JJ)
                    END DO
                END DO
            END DO
            DO INDEX = DOF_WELL_ROW_STARTS(I), DOF_WELL_ROW_STARTS(I + 1) - 1
                COLUMN = DOF_WELL_CONNECTIVITY(INDEX)
                C1 = 3 * COLUMN + W_SHIFT
                C2 = 9 * INDEX
                !DEC$ SIMD            
                DO II = 0, 2
                    C3 = C2 + 3 * II
                    DO JJ = 0, 2
                        DEST(C0 + II) = DEST(C0 + II) + MATRIX_C2W_VALS(C3 + JJ) * SOURCE(C1 + JJ)
                    END DO
                END DO
            END DO
        END DO
            !$OMP END DO NOWAIT
        !$OMP END PARALLEL
                                                        !, CONV_CV_SIZE                        
        !$OMP PARALLEL SHARED(DOF_SHIFT, W_SHIFT, W_SIZE, SOURCE, DEST, MATRIX_W_VALS, MATRIX_W2C_VALS, WELL_DOF_CONNECTIVITY, WELL_DOF_ROW_STARTS)
            !$OMP DO SCHEDULE(STATIC) PRIVATE(I, C0, II, INDEX, COLUMN, C1, C2, C3, JJ)
        DO I = 0, W_SIZE - 1
            C0 = 3 * I + W_SHIFT
            !DEC$ SIMD
            DO II = 0, 2
                DEST(C0 + II) = 0.D0
            END DO
            DO INDEX = WELL_DOF_ROW_STARTS(I), WELL_DOF_ROW_STARTS(I + 1) - 1
                COLUMN = WELL_DOF_CONNECTIVITY(INDEX)
                C1 = 3 * COLUMN + DOF_SHIFT
                C2 = 9 * INDEX
                !DEC$ SIMD            
                DO II = 0, 2
                    C3 = C2 + 3 * II
                    DO JJ = 0, 2
                        DEST(C0 + II) = DEST(C0 + II) + MATRIX_W2C_VALS(C3 + JJ) * SOURCE(C1 + JJ)
                    END DO
                END DO
            END DO
            C2 = 9 * I
            !DEC$ SIMD
            DO II = 0, 2
                C3 = C2 + 3 * II
                DO JJ = 0, 2
                    DEST(C0 + II) = DEST(C0 + II) + MATRIX_W_VALS(C3 + JJ) * SOURCE(C0 + JJ)
                END DO
            END DO
        END DO
            !$OMP END DO NOWAIT
        !$OMP END PARALLEL
        
    END SUBROUTINE

