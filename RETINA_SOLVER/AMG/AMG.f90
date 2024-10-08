    
    
    SUBROUTINE AMG_CREATE__(CONVENTIONAL_CV)
        !DEC$ ATTRIBUTES DLLEXPORT, ALIAS : 'AMG_CREATE__':: RETINA_SOLVER
        USE AMG_MOD
        INTEGER, INTENT(IN) :: CONVENTIONAL_CV
        MAX_LEVEL = 0        
        !IF (CONVENTIONAL_CV .LE. DIRECT_SIZE) THEN
        !    MAX_LEVEL = 0
        !    ALLOCATE(MATS(0:0))
        !ELSE IF (CONVENTIONAL_CV .LE. 8 * DIRECT_SIZE) THEN
        !    MAX_LEVEL = 1
        !    ALLOCATE(MATS(0:1))        
        !    ALLOCATE(LEVELS(0:0))
        !ELSE
        !    MAX_LEVEL = CEILING(LOG(DIRECT_SIZE * 8.D0 / CONVENTIONAL_CV) / LOG(5.D-1))
        !    ALLOCATE(MATS(0:MAX_LEVEL))
        !    ALLOCATE(LEVELS(0:MAX_LEVEL - 1))
        !END IF
    END SUBROUTINE
    
    SUBROUTINE AMG_DESTRUCT__()
        !DEC$ ATTRIBUTES DLLEXPORT, ALIAS : 'AMG_DESTRUCT__':: RETINA_SOLVER
        USE AMG_MOD
        !DEALLOCATE(MATS)
        !IF (MAX_LEVEL .GT. 0) THEN
        !    DEALLOCATE(LEVELS)
        !END IF
    END SUBROUTINE
                                               !, CONVENTIONAL_CV                               !, ACTIVE_MAP                                  !, ACTIVE_DOFS 
    SUBROUTINE AMG_SETUP__(COARSENING_TYPE, CV, ACTIVE_DOF, WELL, CONNECTION, DOF_INDICES_SIZE, DOF_STARTS, DOF_INDICES, DOF_VALUES, DOF_WELL_STARTS, DOF_WELL_INDICES, DOF_WELL_VALUES, WELL_DOF_STARTS, WELL_DOF_INDICES, WELL_DOF_VALUES, WELL_VALUES, DIAG_INDICES, IS_PBH_SOLVE, NATURAL_IDS, ACTIVES)
        !DEC$ ATTRIBUTES DLLEXPORT, ALIAS : 'AMG_SETUP__':: RETINA_SOLVER
        USE AMG_MOD
        INTEGER :: COARSENING_TYPE
        INTEGER, INTENT(IN) :: CV, ACTIVE_DOF, WELL, CONNECTION, DOF_INDICES_SIZE  !CONVENTIONAL_CV, 
        !INTEGER, DIMENSION(0:ACTIVE_DOF - 1), INTENT(IN) :: ACTIVE_DOFS
        !INTEGER, DIMENSION(0:ACTIVE_DOF - 1), INTENT(IN) :: ACTIVE_MAP
        INTEGER, DIMENSION(0:ACTIVE_DOF), INTENT(IN) :: DOF_STARTS
        INTEGER, DIMENSION(0:DOF_INDICES_SIZE - 1), INTENT(IN) :: DOF_INDICES
        DOUBLE PRECISION, DIMENSION(0:9 * DOF_INDICES_SIZE - 1), INTENT(IN) :: DOF_VALUES
        INTEGER, DIMENSION(0:ACTIVE_DOF), INTENT(IN) :: DOF_WELL_STARTS
        INTEGER, DIMENSION(0:CONNECTION - 1), INTENT(IN) :: DOF_WELL_INDICES
        DOUBLE PRECISION, DIMENSION(0:9 * CONNECTION - 1), INTENT(IN) :: DOF_WELL_VALUES
        INTEGER, DIMENSION(0:WELL), INTENT(IN) :: WELL_DOF_STARTS
        INTEGER, DIMENSION(0:CONNECTION - 1), INTENT(IN) :: WELL_DOF_INDICES
        DOUBLE PRECISION, DIMENSION(0:9 * CONNECTION - 1), INTENT(IN) :: WELL_DOF_VALUES
        DOUBLE PRECISION, DIMENSIOn(0:9 * WELL - 1), INTENT(IN) :: WELL_VALUES
        INTEGER, DIMENSION(0:ACTIVE_DOF-1), INTENT(IN) :: DIAG_INDICES
        LOGICAL, DIMENSION(0:WELL-1), INTENT(IN) :: IS_PBH_SOLVE
        INTEGER, DIMENSION(0:CV - 1), INTENT(IN) :: NATURAL_IDS, ACTIVES
        LOGICAL :: OK
        INTEGER :: LEVEL, MAT_NO, C_SIZE, CRITICAL_SIZE
        INTEGER :: I
        OK = .TRUE.                         !, CONVENTIONAL_CV                             !, ACTIVE_MAP                                   !, ACTIVE_DOFS
        CALL AMG_CREATE_PRESSURE_MATRIX_(CV, ACTIVE_DOF, WELL, CONNECTION, DOF_INDICES_SIZE, DOF_STARTS, DOF_INDICES, DOF_VALUES, DOF_WELL_STARTS, DOF_WELL_INDICES, DOF_WELL_VALUES, WELL_DOF_STARTS, WELL_DOF_INDICES, WELL_DOF_VALUES, WELL_VALUES, DIAG_INDICES, IS_PBH_SOLVE, NATURAL_IDS, ACTIVES)        
        !OK = AMG_TEST_NEGATIVE_STRONG_CONNECTIVITY_(MATS(0)%ROWS, MATS(0)%LENGTH, MATS(0)%STARTS, MATS(0)%INDICES, MATS(0)%VALUES)
        !IF (.NOT. OK) THEN
        !    WRITE (*,*) 'UNFORTUNATELY THERE IS A POSITIVE STRONG CONNECTIVITY IN THE MATRIX.'
        !END IF
        !WRITE (*,*) 'IMPES TRANSFORMATION COMPLETED'
        WRITE (*,*) 'MATS(0)%ROWS = ', MATS(0)%ROWS
        MAX_LEVEL = 0
        LEVEL = -1
        C_SIZE = MATS(0)%ROWS        
        CRITICAL_SIZE = DIRECT_SIZE        
        DO WHILE ((C_SIZE .GT. CRITICAL_SIZE) .AND. OK)
            LEVEL = LEVEL + 1
            MAX_LEVEL = MAX_LEVEL + 1
            MAT_NO = LEVEL
            !WRITE (*,*) 'MAT_NORM = ', GET_MATRIX_NORM_(MATS(MAT_NO)%LENGTH, MATS(MAT_NO)%VALUES)
            !WRITE (*,*) 'NNZ = ', MATS(MAT_NO)%LENGTH
            !WRITE (*,*) 'MAX = ', CRS_GET_MAT_MAX_(MATS(MAT_NO)%ROWS, MATS(MAT_NO)%LENGTH, MATS(MAT_NO)%STARTS, MATS(MAT_NO)%VALUES)
            WRITE (*,*) 'CREATING RELAXER'
            CALL AMG_CREATE_RELAXER_(MAT_NO, .FALSE., .TRUE., .FALSE., .FALSE., .FALSE.)
            !CALL AMG_CREATE_ILU0_(MAT_NO, .FALSE., .FALSE., .FALSE., .FALSE., .TRUE.)            
            IF (LEVEL .EQ. 0) THEN
                MATS(MAT_NO)%W_SIZE = WELL
                LEVELS(LEVEL)%COARSENING_TYPE = COARSENING_TYPE
            ELSE
                MATS(MAT_NO)%W_SIZE = 0
                LEVELS(LEVEL)%COARSENING_TYPE = 0
            END IF
            WRITE (*,*) 'COARSENING: METHOD = ', LEVELS(LEVEL)%COARSENING_TYPE
            CALL AMG_GENERATE_COARSE_SET_(LEVEL)            
            C_SIZE = LEVELS(LEVEL)%C_SIZE
            CRITICAL_SIZE = DIRECT_SIZE + LEVELS(LEVEL)%MIN_C
            WRITE (*,*) 'COARSE SET WAS GENERATED SUCCESSFULLY: C_SIZE = ', C_SIZE
            !WRITE (*,*) 'PRINTING FORCED CORASE CELLS...'
            !DO I = 0, MATS(MAT_NO)%ROWS - 1
            !    IF (LEVELS(LEVEL)%FORCED_C(I)) THEN
            !        WRITE (*,*) '@ ROW = ', I
            !    END IF
            !END DO
            !WRITE (*,*) 'GENERATING INTERPOLATOR'
            CALL AMG_GENERATE_INTERPOLATOR_(LEVEL)
            !WRITE (*,*) 'INTERPOLATOR WAS GENERATED SUCCESSFULLY: NNZ = ', LEVELS(LEVEL)%IFC%LENGTH
            !IF (LEVELS(LEVEL)%COARSENING_TYPE .EQ. 1) THEN
            !    !WRITE (*,*) 'SMOOTHING'
            !    CALL AMG_JACOBI_F_SMOOTHING_(LEVEL)
            !    !WRITE (*,*) 'SMOOTHED INTERPOLATOR WAS GENERATED SUCCESSFULLY: NNZ = ', LEVELS(LEVEL)%IFC%LENGTH
            !END IF
            WRITE (*,*) 'PROJECTING MATRIX'
            CALL AMG_GENERATE_MATRIX_(LEVEL)
            WRITE (*,*) 'COARSE MATRIX WAS GENERATED SUCCESSFULLY: NNZ = ', MATS(MAT_NO + 1)%LENGTH
            !OK = AMG_TEST_NEGATIVE_STRONG_CONNECTIVITY_(MATS(MAT_NO + 1)%ROWS, MATS(MAT_NO + 1)%LENGTH, MATS(MAT_NO + 1)%STARTS, MATS(MAT_NO + 1)%INDICES, MATS(MAT_NO + 1)%VALUES)
            !IF (.NOT. OK) THEN
            !    WRITE (*,*) 'UNFORTUNATELY THERE IS A POSITIVE STRONG CONNECTIVITY IN THE MATRIX.'
            !END IF         
            !OK = .FALSE.
        END DO
        MAT_NO = MAX_LEVEL
        WRITE (*,*) 'MAT_NORM = ', GET_MATRIX_NORM_(MATS(MAT_NO)%LENGTH, MATS(MAT_NO)%VALUES)
        !WRITE (*,*) 'MAX = ', CRS_GET_MAT_MAX_(MATS(MAT_NO)%ROWS, MATS(MAT_NO)%LENGTH, MATS(MAT_NO)%STARTS, MATS(MAT_NO)%VALUES)
        CALL AMG_CREATE_ILU0_(MAT_NO, .FALSE., .FALSE., .FALSE., .FALSE., .TRUE.)
        !WRITE (*,*) 'ILU0_NORM = ', GET_MATRIX_NORM_(MATS(MAT_NO)%ILU%LENGTH, MATS(MAT_NO)%ILU%VALUES)
    END SUBROUTINE
                                                                        !, ACTIVE_DOFS
    SUBROUTINE AMG_V_CYCLE__(ROWS, VECTOR, U, N1, N2, ACTIVE_DOF, W_SIZE, LEVEL_STARTS, THREADS, INDEPENDENT_ROW_STARTS, INDEPENDENT_ROWS)
        !DEC$ ATTRIBUTES DLLEXPORT, ALIAS : 'AMG_V_CYCLE__':: RETINA_SOLVER        
        USE AMG_MOD
        USE OMP_LIB
        INTEGER, INTENT(IN) :: ROWS, N1, N2, ACTIVE_DOF, W_SIZE
        !INTEGER, DIMENSION(0:ACTIVE_DOF - 1), INTENT(IN) :: ACTIVE_DOFS
        DOUBLE PRECISION, DIMENSION(0:3 * ROWS - 1), INTENT(IN) :: VECTOR
        DOUBLE PRECISION, DIMENSION(0:3 * ROWS - 1), INTENT(OUT) :: U
        INTEGER, INTENT(IN) ::  LEVEL_STARTS, THREADS
        INTEGER, DIMENSION(0:LEVEL_STARTS - 1), INTENT(IN) :: INDEPENDENT_ROW_STARTS
        INTEGER, DIMENSION(0:ACTIVE_DOF - 1), INTENT(IN) :: INDEPENDENT_ROWS
                                                             !, ACTIVE_DOFS           
        CALL AMG_CREATE_PRESSURE_VECTOR_(ROWS, ACTIVE_DOF, W_SIZE, MATS(0)%VEC, VECTOR)
        CALL AMG_ZERO_INITIAL_GUESS_(MATS(0)%ROWS, MATS(0)%U)
        !call CRS_RESIDUAL_(ROWS, ROWS, MATS(0)%LENGTH, MATS(0)%STARTS, MATS(0)%INDICES, MATS(0)%VALUES, MATS(0)%VEC, MATS(0)%U, MATS(0)%Y)    
        !WRITE (*,*) 'INITIAL RESIDUAL = ', CRS_VECTOR_NORM_(ROWS, MATS(0)%Y)
        !CALL AMG_V_CYCLE_(0, N1, N2)   
        !call CRS_RESIDUAL_(ROWS, ROWS, MATS(0)%LENGTH, MATS(0)%STARTS, MATS(0)%INDICES, MATS(0)%VALUES, MATS(0)%VEC, MATS(0)%U, MATS(0)%Y)    
        !WRITE (*,*) 'MIDDLE RESIDUAL = ', CRS_VECTOR_NORM_(ROWS, MATS(0)%Y)
        CALL AMG_V_CYCLE_(0, N1, N2, LEVEL_STARTS, THREADS, INDEPENDENT_ROW_STARTS, INDEPENDENT_ROWS)
        !call CRS_RESIDUAL_(ROWS, ROWS, MATS(0)%LENGTH, MATS(0)%STARTS, MATS(0)%INDICES, MATS(0)%VALUES, MATS(0)%VEC, MATS(0)%U, MATS(0)%Y)    
        !WRITE (*,*) 'FINAL RESIDUAL = ', CRS_VECTOR_NORM_(ROWS, MATS(0)%Y)
                                                       !, ACTIVE_DOFS 
        CALL AMG_FILL_RESULT_(ROWS, ACTIVE_DOF, W_SIZE, MATS(0)%U, U)
    END SUBROUTINE
    
    SUBROUTINE AMG_DEALLOCATE__()
        !DEC$ ATTRIBUTES DLLEXPORT, ALIAS : 'AMG_DEALLOCATE__':: RETINA_SOLVER
        USE AMG_MOD
        INTEGER :: I
        DEALLOCATE(CPR_LP)
        DO I = 0, MAX_LEVEL
            IF (I .LT. MAX_LEVEL) THEN
                DEALLOCATE(LEVELS(I)%FORCED_C)
                DEALLOCATE(LEVELS(I)%C)
                DEALLOCATE(LEVELS(I)%F)
                DEALLOCATE(LEVELS(I)%IS_C)
                DEALLOCATE(LEVELS(I)%IS_F)
                DEALLOCATE(LEVELS(I)%IS_U)
                DEALLOCATE(LEVELS(I)%CF_INDEX)
                LEVELS(I)%F_SIZE = 0
                LEVELS(I)%C_SIZE = 0
                        
                DEALLOCATE(LEVELS(I)%IFC%STARTS)
                DEALLOCATE(LEVELS(I)%IFC%INDICES)
                DEALLOCATE(LEVELS(I)%IFC%VALUES)
            
                DEALLOCATE(LEVELS(I)%IFCT%STARTS)
                DEALLOCATE(LEVELS(I)%IFCT%INDICES)
                DEALLOCATE(LEVELS(I)%IFCT%VALUES)            
            END IF            
            DEALLOCATE(MATS(I)%STARTS)
            DEALLOCATE(MATS(I)%INDICES)
            DEALLOCATE(MATS(I)%VALUES)
            DEALLOCATE(MATS(I)%VEC)
            DEALLOCATE(MATS(I)%U)
            DEALLOCATE(MATS(I)%DIAG_INDICES)
            
            IF (ALLOCATED(MATS(I)%ILU%STARTS)) THEN
                DEALLOCATE(MATS(I)%ILU%STARTS)
                DEALLOCATE(MATS(I)%ILU%INDICES)
                DEALLOCATE(MATS(I)%ILU%DIAG_INDICES)
                DEALLOCATE(MATS(I)%ILU%VALUES)
            END IF
            DEALLOCATE(MATS(I)%Y)
            IF (ALLOCATED(MATS(I)%J_PARAMS%COPY)) THEN
                DEALLOCATE(MATS(I)%J_PARAMS%COPY)
            END IF
            IF (ALLOCATED(MATS(I)%GS_PARAMS%COPY)) THEN
                DEALLOCATE(MATS(I)%GS_PARAMS%COPY)
            END IF
            IF (ALLOCATED(MATS(I)%ILU_PARAMS%COPY)) THEN
                DEALLOCATE(MATS(I)%ILU_PARAMS%COPY)
                DEALLOCATE(MATS(I)%ILU_PARAMS%RHS)
            END IF
            IF (ALLOCATED(MATS(I)%PCG_PARAMS%COPY)) THEN            
                DEALLOCATE(MATS(I)%PCG_PARAMS%R)
                DEALLOCATE(MATS(I)%PCG_PARAMS%Z)
                DEALLOCATE(MATS(I)%PCG_PARAMS%P)
                DEALLOCATE(MATS(I)%PCG_PARAMS%Q)
                DEALLOCATE(MATS(I)%PCG_PARAMS%COPY)
            END IF
            IF (ALLOCATED(MATS(I)%PBICG_PARAMS%COPY)) THEN            
                DEALLOCATE(MATS(I)%PBICG_PARAMS%R)
                DEALLOCATE(MATS(I)%PBICG_PARAMS%P)
                DEALLOCATE(MATS(I)%PBICG_PARAMS%P_HAT)
                DEALLOCATE(MATS(I)%PBICG_PARAMS%S)
                DEALLOCATE(MATS(I)%PBICG_PARAMS%S_HAT)
                DEALLOCATE(MATS(I)%PBICG_PARAMS%T)
                DEALLOCATE(MATS(I)%PBICG_PARAMS%V)
                DEALLOCATE(MATS(I)%PBICG_PARAMS%RR)
                DEALLOCATE(MATS(I)%PBICG_PARAMS%RMIN)
                DEALLOCATE(MATS(I)%PBICG_PARAMS%COPY)
                DEALLOCATE(MATS(I)%PBICG_PARAMS%XMIN)                
            END IF
            
            IF (ALLOCATED(MATS(I)%S_STARTS)) THEN
                DEALLOCATE(MATS(I)%S_STARTS)
                DEALLOCATE(MATS(I)%S_INDICES)
            END IF
            IF (ALLOCATED(MATS(I)%ST_STARTS)) THEN
                DEALLOCATE(MATS(I)%ST_STARTS)
                DEALLOCATE(MATS(I)%ST_INDICES)
            END IF
        END DO
    END SUBROUTINE
    
    SUBROUTINE AMG_PRINT_MATRIX__(MAT_NO)
        !DEC$ ATTRIBUTES DLLEXPORT, ALIAS : 'AMG_PRINT_MATRIX__':: RETINA_SOLVER
        USE AMG_MOD
        INTEGER, INTENT(IN) :: MAT_NO
        CALL CRS_PRINT_(MATS(MAT_NO)%ROWS, MATS(MAT_NO)%LENGTH, MATS(MAT_NO)%STARTS, MATS(MAT_NO)%INDICES, MATS(MAT_NO)%VALUES)
        CALL CRS_PRINT_VECS_(MATS(MAT_NO)%ROWS, MATS(MAT_NO)%VEC, MATS(MAT_NO)%U, CPR_LP)
    END SUBROUTINE
    
    SUBROUTINE AMG_PRINT_ARRAY__(SIZE, ARR)
        !DEC$ ATTRIBUTES DLLEXPORT, ALIAS : 'AMG_PRINT_ARRAY__':: RETINA_SOLVER
        USE CRS
        INTEGER, INTENT(IN) :: SIZE
        INTEGER, DIMENSION(0:SIZE - 1), INTENT(IN) :: ARR
        CALL CRS_PRINT_ARRAY_(SIZE, ARR)    
    END SUBROUTINE
    
    FUNCTION AMG_GET_MAX_LEVEL__() RESULT(ANS)
        !DEC$ ATTRIBUTES DLLEXPORT, ALIAS : 'AMG_GET_MAX_LEVEL__':: RETINA_SOLVER
        USE AMG_MOD
        INTEGER :: ANS
        ANS = MAX_LEVEL
    END FUNCTION    