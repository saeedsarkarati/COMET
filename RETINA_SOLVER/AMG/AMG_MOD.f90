MODULE AMG_MOD

    USE RED_BLACK_TREE
    USE LINKED_LIST
    USE HEAP_MOD
    USE CRS
    USE OMP_LIB
    USE ILU0
    USE PCG
    USE GAUSS_SEIDEL
    USE JACOBI
    USE PBICG_STAB
    IMPLICIT NONE
    
    TYPE :: LS_TYPE
        INTEGER :: ROWS, LENGTH, S_LENGTH, MAX_CONNECTION
        INTEGER :: W_SIZE
        INTEGER, DIMENSION(:), ALLOCATABLE :: STARTS, INDICES, DIAG_INDICES
        INTEGER, DIMENSION(:), ALLOCATABLE :: S_STARTS, S_INDICES
        INTEGER, DIMENSION(:), ALLOCATABLE :: ST_STARTS, ST_INDICES
        DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: VALUES, VEC, U
        DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Y
        TYPE(CRS_TYPE) :: ILU
        TYPE(PCG_DATA) :: PCG_PARAMS
        TYPE(PBICG_DATA) :: PBICG_PARAMS
        TYPE(GS_DATA) :: GS_PARAMS
        TYPE(ILU_DATA) :: ILU_PARAMS
        TYPE(J_DATA) :: J_PARAMS
    END TYPE

    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: CPR_LP
    
    TYPE :: LEVEL_TYPE
        INTEGER :: F_SIZE, C_SIZE, COARSENING_TYPE
        INTEGER, DIMENSION(:), ALLOCATABLE :: C, F, CF_INDEX
        LOGICAL, DIMENSION(:), ALLOCATABLE :: IS_C, IS_F, IS_U, FORCED_C        
        INTEGER :: MIN_C
        TYPE(CRS_TYPE) :: IFC, IFCT        
    END TYPE
    
    TYPE(LS_TYPE), DIMENSION(0:19) :: MATS
    TYPE(LEVEL_TYPE), DIMENSION(0:19) :: LEVELS
    
    INTEGER :: MAX_LEVEL
    
    INTEGER, PARAMETER :: LU_ = 0, PCG_ = 1, PBICG_ = 2, GS_ = 3
    
    INTEGER, PARAMETER :: DIRECT_SIZE = 200
    DOUBLE PRECISION, PARAMETER :: EPS_STR = 2.5D-1
    
CONTAINS
    
    !******************************************   QUASI-IMPES PART   ******************************************!
                                             !, CONVENTIONAL_CV                              !, ACTIVE_MAP                                   !, ACTIVE_DOFS  
    SUBROUTINE AMG_CREATE_PRESSURE_MATRIX_(CV, ACTIVE_DOF, WELL, CONNECTION, DOF_INDICES_SIZE, DOF_STARTS, DOF_INDICES, DOF_VALUES, DOF_WELL_STARTS, DOF_WELL_INDICES, DOF_WELL_VALUES, WELL_DOF_STARTS, WELL_DOF_INDICES, WELL_DOF_VALUES, WELL_VALUES, DIAG_INDICES, IS_PBH_SOLVE, NATURAL_IDS, ACTIVES)
        INTEGER, INTENT(IN) :: CV, ACTIVE_DOF, WELL, CONNECTION, DOF_INDICES_SIZE      !, CONVENTIONAL_CV
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
        INTEGER, DIMENSION(0:CV-1), INTENT(IN) :: NATURAL_IDS, ACTIVES        
        INTEGER :: ROWS, LENGTH
        ROWS = ACTIVE_DOF + WELL
        LENGTH = DOF_INDICES_SIZE + 2 * CONNECTION + WELL   !- (CONVENTIONAL_CV - ACTIVE_DOF) 
        MATS(0)%ROWS = ROWS
        MATS(0)%LENGTH = LENGTH
        ALLOCATE(MATS(0)%STARTS(0:ROWS))
        ALLOCATE(MATS(0)%INDICES(0:LENGTH - 1))
        ALLOCATE(MATS(0)%VALUES(0:LENGTH - 1))
        ALLOCATE(MATS(0)%VEC(0:ROWS- 1))
        ALLOCATE(MATS(0)%U(0:ROWS- 1))
        ALLOCATE(MATS(0)%DIAG_INDICES(0:ROWS - 1))
        ALLOCATE(CPR_LP(0:3 * ROWS - 1))
        !CALL AMG_QUASI_IMPES_CREATE_PRESSURE_MATRIX_IMPL_(CV, ACTIVE_DOF, CONNECTION, WELL, DOF_INDICES_SIZE, ROWS, LENGTH, MATS(0)%STARTS, MATS(0)%INDICES, MATS(0)%VALUES, DIAG_INDICES, DOF_VALUES, DOF_STARTS, DOF_INDICES, DOF_WELL_STARTS, DOF_WELL_INDICES, DOF_WELL_VALUES, WELL_VALUES, WELL_DOF_STARTS, WELL_DOF_INDICES, WELL_DOF_VALUES, IS_PBH_SOLVE, MATS(0)%DIAG_INDICES)
        !CALL AMG_MIN_SQUARE_QUASI_IMPES_CREATE_PRESSURE_MATRIX_IMPL_(CV, ACTIVE_DOF, CONNECTION, WELL, DOF_INDICES_SIZE, ROWS, LENGTH, MATS(0)%STARTS, MATS(0)%INDICES, MATS(0)%VALUES, DIAG_INDICES, DOF_VALUES, DOF_STARTS, DOF_INDICES, DOF_WELL_STARTS, DOF_WELL_INDICES, DOF_WELL_VALUES, WELL_VALUES, WELL_DOF_STARTS, WELL_DOF_INDICES, WELL_DOF_VALUES, IS_PBH_SOLVE, MATS(0)%DIAG_INDICES)
                                                      !, CONVENTIONAL_CV                                            !, ACTIVE_MAP                                  !, ACTIVE_DOFS        
        CALL AMG_NATURAL_IMPES_CREATE_PRESSURE_MATRIX_IMPL_(CV, ACTIVE_DOF, CONNECTION, WELL, DOF_INDICES_SIZE, ROWS, LENGTH, MATS(0)%STARTS, MATS(0)%INDICES, MATS(0)%VALUES, DIAG_INDICES, DOF_VALUES, DOF_STARTS, DOF_INDICES, DOF_WELL_STARTS, DOF_WELL_INDICES, DOF_WELL_VALUES, WELL_VALUES, WELL_DOF_STARTS, WELL_DOF_INDICES, WELL_DOF_VALUES, IS_PBH_SOLVE, MATS(0)%DIAG_INDICES)
        MATS(0)%MAX_CONNECTION = CRS_GET_MAX_ROW_SIZE_(ROWS, MATS(0)%STARTS)
    END SUBROUTINE
                                                         !, CONV_CV                                                       !, ACTIVE_MAP                   !, ACTIVE_DOFS
    SUBROUTINE AMG_NATURAL_IMPES_CREATE_PRESSURE_MATRIX_IMPL_(CV, ACTIVE_DOF, CONNECTION, WELL_SIZE, DOF_INDICES_SIZE, ROWS, LENGTH, CPR_STARTS, CPR_INDICES, CPR_VALUES, DIAG_IND, MATRIX_FLOW_VALS, MATRIX_ROW_STARTS, MATRIX_CONNECTIVITY, DOF_WELL_ROW_STARTS, DOF_WELL_CONNECTIVITY, MATRIX_C2W_VALS, MATRIX_W_VALS, WELL_DOF_ROW_STARTS, WELL_DOF_CONNECTIVITY, MATRIX_W2C_VALS, IS_PBH_SOLVE, DIAG_INDICES_P)
        INTEGER, INTENT(IN) :: CV, ACTIVE_DOF, CONNECTION, WELL_SIZE, ROWS, LENGTH, DOF_INDICES_SIZE   !, CONV_CV
        !INTEGER, DIMENSION(0:ACTIVE_DOF - 1), INTENT(IN) :: ACTIVE_DOFS
        !INTEGER, DIMENSION(0:ACTIVE_DOF - 1), INTENT(IN) :: ACTIVE_MAP
        INTEGER, DIMENSION(0:ROWS), INTENT(INOUT) :: CPR_STARTS
        INTEGER, DIMENSION(0:LENGTH-1), INTENT(INOUT) :: CPR_INDICES
        DOUBLE PRECISION, DIMENSION(0:LENGTH-1), INTENT(INOUT) :: CPR_VALUES
        INTEGER, DIMENSION(0:ACTIVE_DOF-1), INTENT(IN) :: DIAG_IND
        INTEGER, DIMENSION(0:ROWS - 1), INTENT(OUT) :: DIAG_INDICES_P
        INTEGER, DIMENSION(0:ACTIVE_DOF), INTENT(IN) :: MATRIX_ROW_STARTS, DOF_WELL_ROW_STARTS
        INTEGER, DIMENSION(0:DOF_INDICES_SIZE - 1), INTENT(IN) :: MATRIX_CONNECTIVITY
        INTEGER, DIMENSION(0:CONNECTION - 1), INTENT(IN) :: DOF_WELL_CONNECTIVITY, WELL_DOF_CONNECTIVITY
        INTEGER, DIMENSION(0:WELL_SIZE), INTENT(IN) :: WELL_DOF_ROW_STARTS
        DOUBLE PRECISION, DIMENSION(0:9 * DOF_INDICES_SIZE - 1), INTENT(IN) :: MATRIX_FLOW_VALS
        DOUBLE PRECISION, DIMENSION(0:9 * WELL_SIZE - 1), INTENT(IN) :: MATRIX_W_VALS
        DOUBLE PRECISION, DIMENSION(0:9 * CONNECTION - 1), INTENT(IN) :: MATRIX_W2C_VALS, MATRIX_C2W_VALS
        LOGICAL, DIMENSION(0:WELL_SIZE-1), INTENT(IN) :: IS_PBH_SOLVE        
        INTEGER :: ROW_INDEX, ROW, COLUMN, II, JJ, CPR_PTR, INDEX, WELL, CRS_PTR_OLD, C0, C1, C2
        DOUBLE PRECISION, DIMENSION(0:8) :: A
        DOUBLE PRECISION :: DET
        LOGICAL :: DIAG_FLAG
        ! CREATING CPR_STARTS
        CPR_PTR = 0        
        DO ROW_INDEX = 0, ACTIVE_DOF - 1
            ROW = ROW_INDEX   !ACTIVE_DOFS(ROW_INDEX)
            CPR_STARTS(ROW_INDEX) = CPR_PTR
            CPR_PTR = CPR_PTR + MATRIX_ROW_STARTS(ROW + 1) - MATRIX_ROW_STARTS(ROW) + DOF_WELL_ROW_STARTS(ROW + 1) - DOF_WELL_ROW_STARTS(ROW)
        END DO
        DO ROW = 0, WELL_SIZE - 1
            CPR_STARTS(ROW + ACTIVE_DOF) = CPR_PTR
            CPR_PTR = CPR_PTR + WELL_DOF_ROW_STARTS(ROW + 1) - WELL_DOF_ROW_STARTS(ROW) + 1
        END DO
        CPR_STARTS(ROW + ACTIVE_DOF) = CPR_PTR
        !DEC$ SIMD
        DO II = 0, LENGTH - 1
            CPR_VALUES(II) = 0.D0
        END DO
        DO II = 0, ROWS - 1
            INDEX = 3 * II
            CPR_LP(INDEX) = 0.D0
            CPR_LP(INDEX + 1) = 0.D0
            CPR_LP(INDEX + 2) = 1.D0    
        END DO
        !$OMP PARALLEL 
            !$OMP DO SCHEDULE(STATIC) PRIVATE(ROW_INDEX, ROW, CPR_PTR, INDEX, COLUMN, C1, WELL)
        DO ROW_INDEX = 0, ACTIVE_DOF - 1
            ROW = ROW_INDEX !ACTIVE_DOFS(ROW_INDEX)
            CPR_PTR = CPR_STARTS(ROW_INDEX)
            !DOF2DOF
            DO INDEX = MATRIX_ROW_STARTS(ROW), MATRIX_ROW_STARTS(ROW + 1) - 1                
                COLUMN = MATRIX_CONNECTIVITY(INDEX) 
                IF (ROW .EQ. COLUMN) THEN
                    DIAG_INDICES_P(ROW_INDEX) = CPR_PTR
                END IF
                CPR_INDICES(CPR_PTR) = COLUMN !ACTIVE_MAP(COLUMN)
                C1 = 9 * INDEX
                CPR_VALUES(CPR_PTR) = CPR_VALUES(CPR_PTR) + MATRIX_FLOW_VALS(C1 + 8)
                CPR_PTR = CPR_PTR + 1
            END DO
            !DOF2WELL
            DO INDEX = DOF_WELL_ROW_STARTS(ROW), DOF_WELL_ROW_STARTS(ROW + 1) - 1
                WELL = DOF_WELL_CONNECTIVITY(INDEX)
                COLUMN = WELL + ACTIVE_DOF
                CPR_INDICES(CPR_PTR) = COLUMN 
                IF (IS_PBH_SOLVE(WELL)) THEN
                    C1 = 9 * INDEX
                    CPR_VALUES(CPR_PTR) = CPR_VALUES(CPR_PTR) + MATRIX_C2W_VALS(C1 + 8)
                END IF
                CPR_PTR = CPR_PTR + 1
            END DO
        END DO
            !$OMP END DO
        !$OMP END PARALLEL
        
        !ROWS OF WELL
        !$OMP PARALLEL 
            !$OMP DO SCHEDULE(STATIC) PRIVATE(WELL, C0, ROW_INDEX, ROW, CPR_PTR, INDEX, COLUMN, C1)
        DO WELL = 0, WELL_SIZE - 1
            C0 = 9 * WELL
            ROW_INDEX = WELL + ACTIVE_DOF
            ROW = WELL + ACTIVE_DOF
            CPR_PTR = CPR_STARTS(ROW_INDEX)
            IF((.NOT. IS_PBH_SOLVE(WELL)) .OR. (MATRIX_W_VALS(C0 + 8) .EQ. 1.D0)) THEN
                DO INDEX = WELL_DOF_ROW_STARTS(WELL), WELL_DOF_ROW_STARTS(WELL + 1) - 1
                    COLUMN = WELL_DOF_CONNECTIVITY(INDEX)
                    CPR_INDICES(CPR_PTR) = COLUMN   !ACTIVE_MAP(COLUMN)
                    CPR_PTR = CPR_PTR + 1
                END DO
                CPR_INDICES(CPR_PTR) = ROW_INDEX
                CPR_VALUES(CPR_PTR) = 1.D0
                DIAG_INDICES_P(ROW_INDEX) = CPR_PTR
                CPR_PTR = CPR_PTR + 1
                CPR_LP(3 * ROW_INDEX + 2) = 0.D0
                CYCLE
            END IF
            !WELL2DOF
            DO INDEX = WELL_DOF_ROW_STARTS(WELL), WELL_DOF_ROW_STARTS(WELL + 1) - 1
                COLUMN = WELL_DOF_CONNECTIVITY(INDEX)
                CPR_INDICES(CPR_PTR) = COLUMN   !ACTIVE_MAP(COLUMN)
                C1 = 9 * INDEX
                CPR_VALUES(CPR_PTR) = CPR_VALUES(CPR_PTR) + MATRIX_W2C_VALS(C1 + 8)
                CPR_PTR = CPR_PTR + 1                
            END DO
            !WELL2WELL
            CPR_INDICES(CPR_PTR) = ROW_INDEX
            CPR_VALUES(CPR_PTR) = CPR_VALUES(CPR_PTR) + MATRIX_W_VALS(C0 + 8)
            DIAG_INDICES_P(ROW_INDEX) = CPR_PTR
            CPR_PTR = CPR_PTR + 1
        END DO
            !$OMP END DO
        !$OMP END PARALLEL
    END SUBROUTINE
        
    SUBROUTINE AMG_QUASI_IMPES_CREATE_PRESSURE_MATRIX_IMPL_(CV, ACTIVE_DOF, CONNECTION, WELL_SIZE, DOF_INDICES_SIZE, ROWS, LENGTH, CPR_STARTS, CPR_INDICES, CPR_VALUES, DIAG_IND, MATRIX_FLOW_VALS, MATRIX_ROW_STARTS, MATRIX_CONNECTIVITY, DOF_WELL_ROW_STARTS, DOF_WELL_CONNECTIVITY, MATRIX_C2W_VALS, MATRIX_W_VALS, WELL_DOF_ROW_STARTS, WELL_DOF_CONNECTIVITY, MATRIX_W2C_VALS, IS_PBH_SOLVE, DIAG_INDICES_P)
        INTEGER, INTENT(IN) :: CV, ACTIVE_DOF, CONNECTION, WELL_SIZE, ROWS, LENGTH, DOF_INDICES_SIZE        
        INTEGER, DIMENSION(0:ROWS), INTENT(INOUT) :: CPR_STARTS
        INTEGER, DIMENSION(0:LENGTH-1), INTENT(INOUT) :: CPR_INDICES
        DOUBLE PRECISION, DIMENSION(0:LENGTH-1), INTENT(INOUT) :: CPR_VALUES
        INTEGER, DIMENSION(0:ACTIVE_DOF-1), INTENT(IN) :: DIAG_IND
        INTEGER, DIMENSION(0:ROWS - 1), INTENT(OUT) :: DIAG_INDICES_P
        INTEGER, DIMENSION(0:ACTIVE_DOF), INTENT(IN) :: MATRIX_ROW_STARTS, DOF_WELL_ROW_STARTS
        INTEGER, DIMENSION(0:DOF_INDICES_SIZE - 1), INTENT(IN) :: MATRIX_CONNECTIVITY
        INTEGER, DIMENSION(0:CONNECTION - 1), INTENT(IN) :: DOF_WELL_CONNECTIVITY, WELL_DOF_CONNECTIVITY
        INTEGER, DIMENSION(0:WELL_SIZE), INTENT(IN) :: WELL_DOF_ROW_STARTS
        DOUBLE PRECISION, DIMENSION(0:9 * DOF_INDICES_SIZE - 1), INTENT(IN) :: MATRIX_FLOW_VALS
        DOUBLE PRECISION, DIMENSION(0:9 * WELL_SIZE - 1), INTENT(IN) :: MATRIX_W_VALS
        DOUBLE PRECISION, DIMENSION(0:9 * CONNECTION - 1), INTENT(IN) :: MATRIX_W2C_VALS, MATRIX_C2W_VALS
        LOGICAL, DIMENSION(0:WELL_SIZE-1), INTENT(IN) :: IS_PBH_SOLVE
        INTEGER :: ROW_INDEX, ROW, COLUMN, II, JJ, KK, CPR_PTR, INDEX, LP_INDEX, WELL, CRS_PTR_OLD, C0
        DOUBLE PRECISION, DIMENSION(0:8) :: DIAG, DIAG_INV
        DOUBLE PRECISION :: DET
        LOGICAL :: DIAG_FLAG
        DOUBLE PRECISION :: P_VAL
        ! CREATING CPR_STARTS        
        CPR_PTR = 0        
        DO ROW_INDEX = 0, ACTIVE_DOF - 1
            ROW = ROW_INDEX
            CPR_STARTS(ROW_INDEX) = CPR_PTR
            CPR_PTR = CPR_PTR + MATRIX_ROW_STARTS(ROW + 1) - MATRIX_ROW_STARTS(ROW) + DOF_WELL_ROW_STARTS(ROW + 1) - DOF_WELL_ROW_STARTS(ROW)
        END DO
        DO ROW = 0, WELL_SIZE - 1
            CPR_STARTS(ROW + ACTIVE_DOF) = CPR_PTR
            CPR_PTR = CPR_PTR + WELL_DOF_ROW_STARTS(ROW + 1) - WELL_DOF_ROW_STARTS(ROW) + 1
        END DO
        CPR_STARTS(ROW + ACTIVE_DOF) = CPR_PTR
        !DEC$ SIMD
        DO II = 0, LENGTH - 1
            CPR_VALUES(II) = 0.D0
        END DO
        !$OMP PARALLEL 
            !$OMP DO SCHEDULE(STATIC) PRIVATE(ROW_INDEX, ROW, CPR_PTR, LP_INDEX, C0, DIAG, DIAG_INV, II, JJ, KK, INDEX, COLUMN, P_VAL, WELL)
        DO ROW_INDEX = 0, ACTIVE_DOF - 1
            ROW = ROW_INDEX
            CPR_PTR = CPR_STARTS(ROW_INDEX)
            LP_INDEX = 3 * ROW_INDEX
            ! CPR_LP SHOULD BE CREATED
            C0 = 9 * DIAG_IND(ROW)
            DIAG = MATRIX_FLOW_VALS(C0:C0 + 8)
            DIAG_INV = 0.D0
            DIAG_INV(8) = 1.D0
            DO II = 0, 1
                DO JJ = 0, 1
                    DIAG_INV(3 * II + JJ) = DIAG(3 * II + JJ)
                END DO                
            END DO
            CALL INV__(DIAG_INV)
            DO JJ = 0, 1
                CPR_LP(LP_INDEX + JJ) = 0.D0
                DO KK = 0, 1
                    CPR_LP(LP_INDEX + JJ) = CPR_LP(LP_INDEX + JJ) - DIAG(6 + KK) * DIAG_INV(3 * KK + JJ)
                END DO
            END DO
            IF ((DIAG(8) + CPR_LP(LP_INDEX + 0) * DIAG(2) + CPR_LP(LP_INDEX + 1) * DIAG(5)) .LT. 0.D0) THEN
                CPR_LP(LP_INDEX + 2) = -1.D0
            ELSE
                CPR_LP(LP_INDEX + 2) = 1.D0
            END IF
            !CPR_LP(LP_INDEX + 0) = 0.D0
            !CPR_LP(LP_INDEX + 1) = 0.D0
            !CPR_LP(LP_INDEX + 2) = 1.D0
            !DOF2DOF
            DO INDEX = MATRIX_ROW_STARTS(ROW), MATRIX_ROW_STARTS(ROW + 1) - 1                
                COLUMN = MATRIX_CONNECTIVITY(INDEX) 
                IF (ROW .EQ. COLUMN) THEN
                    DIAG_INDICES_P(ROW_INDEX) = CPR_PTR
                END IF
                CPR_INDICES(CPR_PTR) = COLUMN
                C0 = 9 * INDEX
                P_VAL = MATRIX_FLOW_VALS(C0 + 8)
                DO JJ = 0, 1
                    P_VAL = P_VAL + CPR_LP(LP_INDEX + JJ) * MATRIX_FLOW_VALS(C0 + 3 * JJ + 2)
                END DO                
                CPR_VALUES(CPR_PTR) = P_VAL * CPR_LP(LP_INDEX + 2)
                CPR_PTR = CPR_PTR + 1
            END DO
            !DOF2WELL
            DO INDEX = DOF_WELL_ROW_STARTS(ROW), DOF_WELL_ROW_STARTS(ROW + 1) - 1
                WELL = DOF_WELL_CONNECTIVITY(INDEX)
                COLUMN = WELL + ACTIVE_DOF
                CPR_INDICES(CPR_PTR) = COLUMN 
                IF (IS_PBH_SOLVE(WELL)) THEN
                    C0 = 9 * INDEX
                    P_VAL = MATRIX_C2W_VALS(C0 + 8)
                    DO JJ = 0, 1
                        P_VAL = P_VAL + CPR_LP(LP_INDEX + JJ) * MATRIX_C2W_VALS(C0 + 3 * JJ + 2)
                    END DO
                    CPR_VALUES(CPR_PTR) = P_VAL * CPR_LP(LP_INDEX + 2)
                END IF
                CPR_PTR = CPR_PTR + 1
            END DO
        END DO
            !$OMP END DO
        !$OMP END PARALLEL
        
        !ROWS OF WELL
        !$OMP PARALLEL 
            !$OMP DO SCHEDULE(STATIC) PRIVATE(WELL, C0, ROW_INDEX, ROW, CPR_PTR, LP_INDEX, INDEX, COLUMN, JJ, DIAG, DIAG_INV, II, KK, P_VAL)
        DO WELL = 0, WELL_SIZE - 1
            C0 = 9 * WELL
            ROW_INDEX = WELL + ACTIVE_DOF            
            ROW = WELL + ACTIVE_DOF
            CPR_PTR = CPR_STARTS(ROW_INDEX)
            LP_INDEX = 3 * ROW_INDEX
            IF((.NOT. IS_PBH_SOLVE(WELL)) .OR. (MATRIX_W_VALS(C0 + 8) .EQ. 1.D0)) THEN
                DO INDEX = WELL_DOF_ROW_STARTS(WELL), WELL_DOF_ROW_STARTS(WELL + 1) - 1
                    COLUMN = WELL_DOF_CONNECTIVITY(INDEX)
                    CPR_INDICES(CPR_PTR) = COLUMN
                    CPR_PTR = CPR_PTR + 1
                END DO
                CPR_INDICES(CPR_PTR) = ROW_INDEX
                CPR_VALUES(CPR_PTR) = 1.D0
                DIAG_INDICES_P(ROW_INDEX) = CPR_PTR
                CPR_PTR = CPR_PTR + 1
                DO JJ = 0, 2                    
                    CPR_LP(LP_INDEX + JJ) = 0.D0
                END DO
                CYCLE
            END IF
            ! CPR_LP SHOULD BE CREATED            
            DIAG = MATRIX_W_VALS(C0:C0 + 8)
            DIAG_INV = 0.D0
            DIAG_INV(8) = 1.D0
            DO II = 0, 1
                DO JJ = 0, 1
                    DIAG_INV(3 * II + JJ) = DIAG(3 * II + JJ)
                END DO                
            END DO
            CALL INV__(DIAG_INV)            
            DO JJ = 0, 1
                CPR_LP(LP_INDEX + JJ) = 0.D0
                DO KK = 0, 1
                    CPR_LP(LP_INDEX + JJ) = CPR_LP(LP_INDEX + JJ) - DIAG(6 + KK) * DIAG_INV(3 * KK + JJ)
                END DO
            END DO
            IF ((DIAG(8) + CPR_LP(LP_INDEX + 0) * DIAG(2) + CPR_LP(LP_INDEX + 1) * DIAG(5)) .LT. 0.D0) THEN
                CPR_LP(LP_INDEX + 2) = -1.D0
                CPR_LP(LP_INDEX + 1) = -CPR_LP(LP_INDEX + 1)
                CPR_LP(LP_INDEX + 0) = -CPR_LP(LP_INDEX + 0)
            ELSE
                CPR_LP(LP_INDEX + 2) = 1.D0
            END IF
            !CPR_LP(LP_INDEX + 0) = 0.D0
            !CPR_LP(LP_INDEX + 1) = 0.D0
            !CPR_LP(LP_INDEX + 2) = 1.D0
            
            
            !WELL2DOF
            DO INDEX = WELL_DOF_ROW_STARTS(WELL), WELL_DOF_ROW_STARTS(WELL + 1) - 1
                COLUMN = WELL_DOF_CONNECTIVITY(INDEX)
                CPR_INDICES(CPR_PTR) = COLUMN
                C0 = 9 * INDEX
                P_VAL = MATRIX_W2C_VALS(C0 + 8)
                DO JJ = 0, 1
                    P_VAL = P_VAL + CPR_LP(LP_INDEX + JJ) * MATRIX_W2C_VALS(C0 + 3 * JJ + 2)
                END DO
                CPR_VALUES(CPR_PTR) = P_VAL * CPR_LP(LP_INDEX + 2)
                CPR_PTR = CPR_PTR + 1                
            END DO
            !WELL2WELL
            CPR_INDICES(CPR_PTR) = ROW_INDEX
            C0 = 9 * WELL
            P_VAL = MATRIX_W_VALS(C0 + 8)
            DO JJ = 0, 1
                P_VAL = P_VAL + CPR_LP(LP_INDEX + JJ) * MATRIX_W_VALS(C0 + 3 * JJ + 2)
            END DO
            CPR_VALUES(CPR_PTR) = P_VAL * CPR_LP(LP_INDEX + 2)
            DIAG_INDICES_P(ROW_INDEX) = CPR_PTR
            CPR_PTR = CPR_PTR + 1
        END DO
            !$OMP END DO
        !$OMP END PARALLEL
    END SUBROUTINE

    SUBROUTINE AMG_MIN_SQUARE_QUASI_IMPES_CREATE_PRESSURE_MATRIX_IMPL_(CV, ACTIVE_DOF, CONNECTION, WELL_SIZE, DOF_INDICES_SIZE, ROWS, LENGTH, CPR_STARTS, CPR_INDICES, CPR_VALUES, DIAG_IND, MATRIX_FLOW_VALS, MATRIX_ROW_STARTS, MATRIX_CONNECTIVITY, DOF_WELL_ROW_STARTS, DOF_WELL_CONNECTIVITY, MATRIX_C2W_VALS, MATRIX_W_VALS, WELL_DOF_ROW_STARTS, WELL_DOF_CONNECTIVITY, MATRIX_W2C_VALS, IS_PBH_SOLVE, DIAG_INDICES_P)
        INTEGER, INTENT(IN) :: CV, ACTIVE_DOF, CONNECTION, WELL_SIZE, ROWS, LENGTH, DOF_INDICES_SIZE        
        INTEGER, DIMENSION(0:ROWS), INTENT(INOUT) :: CPR_STARTS
        INTEGER, DIMENSION(0:LENGTH-1), INTENT(INOUT) :: CPR_INDICES
        DOUBLE PRECISION, DIMENSION(0:LENGTH-1), INTENT(INOUT) :: CPR_VALUES
        INTEGER, DIMENSION(0:ACTIVE_DOF-1), INTENT(IN) :: DIAG_IND
        INTEGER, DIMENSION(0:ROWS - 1), INTENT(OUT) :: DIAG_INDICES_P
        INTEGER, DIMENSION(0:ACTIVE_DOF), INTENT(IN) :: MATRIX_ROW_STARTS, DOF_WELL_ROW_STARTS
        INTEGER, DIMENSION(0:DOF_INDICES_SIZE - 1), INTENT(IN) :: MATRIX_CONNECTIVITY
        INTEGER, DIMENSION(0:CONNECTION - 1), INTENT(IN) :: DOF_WELL_CONNECTIVITY, WELL_DOF_CONNECTIVITY
        INTEGER, DIMENSION(0:WELL_SIZE), INTENT(IN) :: WELL_DOF_ROW_STARTS
        DOUBLE PRECISION, DIMENSION(0:9 * DOF_INDICES_SIZE - 1), INTENT(IN) :: MATRIX_FLOW_VALS
        DOUBLE PRECISION, DIMENSION(0:9 * WELL_SIZE - 1), INTENT(IN) :: MATRIX_W_VALS
        DOUBLE PRECISION, DIMENSION(0:9 * CONNECTION - 1), INTENT(IN) :: MATRIX_W2C_VALS, MATRIX_C2W_VALS
        LOGICAL, DIMENSION(0:WELL_SIZE-1), INTENT(IN) :: IS_PBH_SOLVE
        INTEGER :: ROW_INDEX, ROW, COLUMN, II, JJ, KK, CPR_PTR, INDEX, LP_INDEX, WELL, CRS_PTR_OLD, C0
        DOUBLE PRECISION, DIMENSION(0:8) :: DIAG, FV
        DOUBLE PRECISION, DIMENSION(0:2) :: A1, A2, A3
        DOUBLE PRECISION :: A11, A12, A21, A22, B1, B2, NOM_DET, DENOM_DET
        DOUBLE PRECISION :: P_VAL
        
        ! CREATING CPR_STARTS        
        CPR_PTR = 0        
        DO ROW_INDEX = 0, ACTIVE_DOF - 1
            ROW = ROW_INDEX
            CPR_STARTS(ROW_INDEX) = CPR_PTR
            CPR_PTR = CPR_PTR + MATRIX_ROW_STARTS(ROW + 1) - MATRIX_ROW_STARTS(ROW) + DOF_WELL_ROW_STARTS(ROW + 1) - DOF_WELL_ROW_STARTS(ROW)
        END DO
        DO ROW = 0, WELL_SIZE - 1
            CPR_STARTS(ROW + ACTIVE_DOF) = CPR_PTR
            CPR_PTR = CPR_PTR + WELL_DOF_ROW_STARTS(ROW + 1) - WELL_DOF_ROW_STARTS(ROW) + 1
        END DO
        CPR_STARTS(ROW + ACTIVE_DOF) = CPR_PTR
        !DEC$ SIMD
        DO II = 0, LENGTH - 1
            CPR_VALUES(II) = 0.D0
        END DO
        !$OMP PARALLEL 
            !$OMP DO SCHEDULE(STATIC) PRIVATE(ROW_INDEX, ROW, CPR_PTR, LP_INDEX, C0, DIAG, A11, A12, A21, A22, B1, B2, INDEX, FV, A1, A2, A3, DENOM_DET, NOM_DET, COLUMN, P_VAL, JJ, WELL)
        DO ROW_INDEX = 0, ACTIVE_DOF - 1
            ROW = ROW_INDEX
            CPR_PTR = CPR_STARTS(ROW_INDEX)
            LP_INDEX = 3 * ROW_INDEX
            ! CPR_LP SHOULD BE CREATED
            C0 = 9 * DIAG_IND(ROW)
            DIAG = MATRIX_FLOW_VALS(C0:C0 + 8)
            A11 = 0.D0
            A12 = 0.D0
            A21 = 0.D0
            A22 = 0.D0
            B1 = 0.D0
            B2 = 0.D0
            DO INDEX = MATRIX_ROW_STARTS(ROW), MATRIX_ROW_STARTS(ROW + 1) - 1
                C0 = 9 * INDEX
                FV = MATRIX_FLOW_VALS(C0:C0 + 8)
                A1 = FV(0:2)
                A2 = FV(3:5)
                A3 = FV(6:8)
                A11 = A11 + A1(0) * A1(0) + A1(1) * A1(1)
                A12 = A12 + A1(0) * A2(0) + A1(1) * A2(1)                
                A22 = A22 + A2(0) * A2(0) + A2(1) * A2(1)
                B1 = B1 - A1(0) * A3(0) - A1(1) * A3(1)
                B2 = B2 - A2(0) * A3(0) - A2(1) * A3(1)
            END DO
            DO INDEX = DOF_WELL_ROW_STARTS(ROW), DOF_WELL_ROW_STARTS(ROW + 1) - 1
                C0 = 9 * INDEX
                FV = MATRIX_C2W_VALS(C0:C0 + 8)
                A1 = FV(0:2)
                A2 = FV(3:5)
                A3 = FV(6:8)
                A11 = A11 + A1(0) * A1(0) + A1(1) * A1(1)
                A12 = A12 + A1(0) * A2(0) + A1(1) * A2(1)                
                A22 = A22 + A2(0) * A2(0) + A2(1) * A2(1)
                B1 = B1 - A1(0) * A3(0) - A1(1) * A3(1)
                B2 = B2 - A2(0) * A3(0) - A2(1) * A3(1)
            END DO
            A21 = A12
            DENOM_DET = A11 * A22 - A12 * A21
            NOM_DET = B1 * A22 - A12 * B2
            CPR_LP(LP_INDEX + 0) = NOM_DET / DENOM_DET
            NOM_DET = A11 * B2 - B1 * A21
            CPR_LP(LP_INDEX + 1) = NOM_DET / DENOM_DET
            IF ((DIAG(8) + CPR_LP(LP_INDEX + 0) * DIAG(2) + CPR_LP(LP_INDEX + 1) * DIAG(5)) .LT. 0.D0) THEN
                CPR_LP(LP_INDEX + 2) = -1.D0
            ELSE
                CPR_LP(LP_INDEX + 2) = 1.D0
            END IF
            
            !DOF2DOF
            DO INDEX = MATRIX_ROW_STARTS(ROW), MATRIX_ROW_STARTS(ROW + 1) - 1                
                COLUMN = MATRIX_CONNECTIVITY(INDEX) 
                IF (ROW .EQ. COLUMN) THEN
                    DIAG_INDICES_P(ROW_INDEX) = CPR_PTR
                END IF
                CPR_INDICES(CPR_PTR) = COLUMN
                C0 = 9 * INDEX
                P_VAL = MATRIX_FLOW_VALS(C0 + 8)
                DO JJ = 0, 1
                    P_VAL = P_VAL + CPR_LP(LP_INDEX + JJ) * MATRIX_FLOW_VALS(C0 + 3 * JJ + 2)
                END DO                
                CPR_VALUES(CPR_PTR) = P_VAL * CPR_LP(LP_INDEX + 2)
                CPR_PTR = CPR_PTR + 1
            END DO
            !DOF2WELL
            DO INDEX = DOF_WELL_ROW_STARTS(ROW), DOF_WELL_ROW_STARTS(ROW + 1) - 1
                WELL = DOF_WELL_CONNECTIVITY(INDEX)
                COLUMN = WELL + ACTIVE_DOF
                CPR_INDICES(CPR_PTR) = COLUMN 
                IF (IS_PBH_SOLVE(WELL)) THEN
                    C0 = 9 * INDEX
                    P_VAL = MATRIX_C2W_VALS(C0 + 8)
                    DO JJ = 0, 1
                        P_VAL = P_VAL + CPR_LP(LP_INDEX + JJ) * MATRIX_C2W_VALS(C0 + 3 * JJ + 2)
                    END DO
                    CPR_VALUES(CPR_PTR) = P_VAL * CPR_LP(LP_INDEX + 2)
                END IF
                CPR_PTR = CPR_PTR + 1
            END DO
        END DO
            !$OMP END DO
        !$OMP END PARALLEL
        
        !ROWS OF WELL
        !$OMP PARALLEL 
            !$OMP DO SCHEDULE(STATIC) PRIVATE(WELL, C0, ROW_INDEX, ROW, CPR_PTR, LP_INDEX, INDEX, COLUMN, JJ, DIAG, FV, A1, A2, A3, A11, A12, A22, B1, B2, A21, DENOM_DET, NOM_DET, P_VAL)
        DO WELL = 0, WELL_SIZE - 1
            C0 = 9 * WELL
            ROW_INDEX = WELL + ACTIVE_DOF            
            ROW = WELL + ACTIVE_DOF
            CPR_PTR = ROW_INDEX
            LP_INDEX = 3 * ROW_INDEX
            IF((.NOT. IS_PBH_SOLVE(WELL)) .OR. (MATRIX_W_VALS(C0 + 8) .EQ. 1.D0)) THEN
                DO INDEX = WELL_DOF_ROW_STARTS(WELL), WELL_DOF_ROW_STARTS(WELL + 1) - 1
                    COLUMN = WELL_DOF_CONNECTIVITY(INDEX)
                    CPR_INDICES(CPR_PTR) = COLUMN
                    CPR_PTR = CPR_PTR + 1
                END DO
                CPR_INDICES(CPR_PTR) = ROW_INDEX
                CPR_VALUES(CPR_PTR) = 1.D0
                DIAG_INDICES_P(ROW_INDEX) = CPR_PTR
                CPR_PTR = CPR_PTR + 1
                DO JJ = 0, 2                    
                    CPR_LP(LP_INDEX + JJ) = 0.D0
                END DO
                CYCLE
            END IF
            ! CPR_LP SHOULD BE CREATED            
            DIAG = MATRIX_W_VALS(C0:C0 + 8)
            FV = DIAG
            A1 = FV(0:2)
            A2 = FV(3:5)
            A3 = FV(6:8)
            A11 = A1(0) * A1(0) + A1(1) * A1(1)
            A12 = A1(0) * A2(0) + A1(1) * A2(1)                
            A22 = A2(0) * A2(0) + A2(1) * A2(1)
            B1 = - A1(0) * A3(0) - A1(1) * A3(1)
            B2 = - A2(0) * A3(0) - A2(1) * A3(1)
            DO INDEX = WELL_DOF_ROW_STARTS(WELL), WELL_DOF_ROW_STARTS(WELL + 1) - 1
                C0 = 9 * INDEX
                FV = MATRIX_W2C_VALS(C0:C0 + 8)
                A1 = FV(0:2)
                A2 = FV(3:5)
                A3 = FV(6:8)
                A11 = A11 + A1(0) * A1(0) + A1(1) * A1(1)
                A12 = A12 + A1(0) * A2(0) + A1(1) * A2(1)                
                A22 = A22 + A2(0) * A2(0) + A2(1) * A2(1)
                B1 = B1 - A1(0) * A3(0) - A1(1) * A3(1)
                B2 = B2 - A2(0) * A3(0) - A2(1) * A3(1)
            END DO
            A21 = A12
            DENOM_DET = A11 * A22 - A12 * A21
            NOM_DET = B1 * A22 - A12 * B2
            CPR_LP(LP_INDEX + 0) = NOM_DET / DENOM_DET
            NOM_DET = A11 * B2 - B1 * A21
            CPR_LP(LP_INDEX + 1) = NOM_DET / DENOM_DET
            IF ((DIAG(8) + CPR_LP(LP_INDEX + 0) * DIAG(2) + CPR_LP(LP_INDEX + 1) * DIAG(5)) .LT. 0.D0) THEN
                CPR_LP(LP_INDEX + 2) = -1.D0
            ELSE
                CPR_LP(LP_INDEX + 2) = 1.D0
            END IF
            
            
            !WELL2DOF
            DO INDEX = WELL_DOF_ROW_STARTS(WELL), WELL_DOF_ROW_STARTS(WELL + 1) - 1
                COLUMN = WELL_DOF_CONNECTIVITY(INDEX)
                CPR_INDICES(CPR_PTR) = COLUMN
                C0 = 9 * INDEX
                P_VAL = MATRIX_W2C_VALS(C0 + 8)
                DO JJ = 0, 1
                    P_VAL = P_VAL + CPR_LP(LP_INDEX + JJ) * MATRIX_W2C_VALS(C0 + 3 * JJ + 2)
                END DO
                CPR_VALUES(CPR_PTR) = P_VAL * CPR_LP(LP_INDEX + 2)
                CPR_PTR = CPR_PTR + 1                
            END DO
            !WELL2WELL
            C0 = 9 * WELL
            CPR_INDICES(CPR_PTR) = ROW_INDEX
            P_VAL = MATRIX_W_VALS(C0 + 8)
            DO JJ = 0, 1
                P_VAL = P_VAL + CPR_LP(LP_INDEX + JJ) * MATRIX_W_VALS(C0 + 3 * JJ + 2)
            END DO
            CPR_VALUES(CPR_PTR) = P_VAL * CPR_LP(LP_INDEX + 2)
            DIAG_INDICES_P(ROW_INDEX) = CPR_PTR
            CPR_PTR = CPR_PTR + 1
        END DO
            !$OMP END DO
        !$OMP END PARALLEL
    END SUBROUTINE
                                                                   !, ACTIVE_DOFS         
    SUBROUTINE AMG_CREATE_PRESSURE_VECTOR_(SIZE, ACTIVE_DOF, W_SIZE, CPR_VECTOR, VECTOR)        
        INTEGER, INTENT(IN) :: SIZE, ACTIVE_DOF, W_SIZE
        !INTEGER, DIMENSION(0:ACTIVE_DOF - 1), INTENT(IN) :: ACTIVE_DOFS
        DOUBLE PRECISION, DIMENSION(0:ACTIVE_DOF + W_SIZE - 1), INTENT(OUT) :: CPR_VECTOR
        DOUBLE PRECISION, DIMENSION(0:3 * SIZE - 1), INTENT(IN) :: VECTOR
        INTEGER :: I_INDEX, I, II, C0, INDEX
        !$OMP PARALLEL
            !$OMP DO SCHEDULE(STATIC) PRIVATE(I_INDEX, I, C0, II)        
        DO I_INDEX = 0, ACTIVE_DOF - 1
            I = I_INDEX !ACTIVE_DOFS(I_INDEX)
            C0 = 3 * I
            CPR_VECTOR(I_INDEX) = VECTOR(C0 + 2)
            DO II = 0, 1
                CPR_VECTOR(I_INDEX) = CPR_VECTOR(I_INDEX) + CPR_LP(3 * I_INDEX + II) * VECTOR(C0 + II)
            END DO      
            CPR_VECTOR(I_INDEX) = CPR_VECTOR(I_INDEX) * CPR_LP(3 * I_INDEX + 2)
        END DO
            !$OMP END DO
        !$OMP END PARALLEL
        
        !$OMP PARALLEL
            !$OMP DO SCHEDULE(STATIC) PRIVATE(I_INDEX, I, INDEX, C0, II)        
        DO I_INDEX = 0, W_SIZE - 1
            I = I_INDEX + SIZE - W_SIZE
            INDEX = I_INDEX + ACTIVE_DOF
            C0 = 3 * I
            CPR_VECTOR(INDEX) = VECTOR(C0 + 2)
            DO II = 0, 1
                CPR_VECTOR(INDEX) = CPR_VECTOR(INDEX) + CPR_LP(3 * INDEX + II) * VECTOR(C0 + II)
            END DO
            CPR_VECTOR(INDEX) = CPR_VECTOR(INDEX) * CPR_LP(3 * INDEX + 2)
        END DO
            !$OMP END DO
        !$OMP END PARALLEL
    END SUBROUTINE
                                                        !, ACTIVE_DOFS
    SUBROUTINE AMG_FILL_RESULT_(SIZE, ACTIVE_DOF, W_SIZE, SOURCE, DEST)
        INTEGER, INTENT(IN) :: SIZE, ACTIVE_DOF, W_SIZE
        !INTEGER, DIMENSION(0:ACTIVE_DOF - 1), INTENT(IN) :: ACTIVE_DOFS
        DOUBLE PRECISION, DIMENSION(0:ACTIVE_DOF + W_SIZE - 1), INTENT(IN) :: SOURCE
        DOUBLE PRECISION, DIMENSION(0:3 * SIZE - 1), INTENT(OUT) :: DEST
        INTEGER :: I_INDEX, I, C0, II
        !$OMP PARALLEL
            !$OMP DO SCHEDULE(STATIC) PRIVATE(I_INDEX, I, C0, II)
        DO I_INDEX = 0, ACTIVE_DOF - 1
            I = I_INDEX     !ACTIVE_DOFS(I_INDEX)
            C0 = 3 * I
            !DEC$ SIMD
            DO II = 0, 1
                DEST(C0 + II) = 0.D0
            END DO
            DEST(C0 + 2) = SOURCE(I_INDEX)
        END DO
            !$OMP END DO
        !$OMP END PARALLEL
        
        !$OMP PARALLEL
            !$OMP DO SCHEDULE(STATIC) PRIVATE(I_INDEX, I, C0, II)
        DO I_INDEX = 0, W_SIZE - 1
            I = I_INDEX + SIZE - W_SIZE
            C0 = 3 * I
            !DEC$ SIMD
            DO II = 0, 1
                DEST(C0 + II) = 0.D0
            END DO
            DEST(C0 + 2) = SOURCE(I_INDEX + ACTIVE_DOF)
        END DO
            !$OMP END DO
        !$OMP END PARALLEL
    END SUBROUTINE
    

    !******************************************   CONNECTIVITY PART   ******************************************!
    
    SUBROUTINE AMG_ALLOCATE_ST_CON_(MAT_NO, LENGTH, SIZE)
        INTEGER, INTENT(IN) :: MAT_NO, LENGTH, SIZE
        ALLOCATE(MATS(MAT_NO)%ST_STARTS(0:SIZE))
        ALLOCATE(MATS(MAT_NO)%ST_INDICES(0:LENGTH - 1))
    END SUBROUTINE
    
    SUBROUTINE AMG_GENERATE_ST_CON_IMPL_(SIZE, LENGTH, TEMP_STARTS, S_STARTS, S_INDICES, ST_STARTS, ST_INDICES)
        ! TEMP_STARTS CONTAINS ST_ROW_SIZES
        INTEGER, INTENT(IN) :: LENGTH, SIZE
        INTEGER, DIMENSION(0:SIZE), INTENT(IN) :: S_STARTS
        INTEGER, DIMENSION(0:LENGTH - 1), INTENT(IN) :: S_INDICES
        INTEGER, DIMENSION(0:SIZE), INTENT(OUT) :: ST_STARTS
        INTEGER, DIMENSION(0:LENGTH - 1), INTENT(OUT) :: ST_INDICES
        INTEGER, DIMENSION(0:SIZE - 1), INTENT(INOUT) :: TEMP_STARTS
        INTEGER :: I, J_INDEX, J
        ! FILLING ST_CON
        ST_STARTS(SIZE) = LENGTH
        DO I = SIZE - 1, 0, -1
            ST_STARTS(I) = ST_STARTS(I + 1) - TEMP_STARTS(I)
            ! TEMP_STARTS IS CHANGING TO ROW STARTS
            TEMP_STARTS(I) = ST_STARTS(I)
        END DO
        
        DO I = 0, SIZE - 1
            DO J_INDEX = S_STARTS(I), S_STARTS(I + 1) - 1
                J = S_INDICES(J_INDEX)
                ST_INDICES(TEMP_STARTS(J)) = I
                ! TEMP_STARTS IS CHANGING TO ROW POINTERS
                TEMP_STARTS(J) = TEMP_STARTS(J) + 1
            END DO
        END DO
    END SUBROUTINE
    
    FUNCTION AMG_TEST_NEGATIVE_STRONG_CONNECTIVITY_(SIZE, LENGTH, STARTS, INDICES, VALUES) RESULT(OK)
        INTEGER, INTENT(IN) :: SIZE, LENGTH
        INTEGER, DIMENSION(0:SIZE), INTENT(IN) :: STARTS
        INTEGER, DIMENSION(0:LENGTH - 1), INTENT(IN) :: INDICES
        DOUBLE PRECISION, DIMENSION(0:LENGTH - 1), INTENT(IN) :: VALUES
        LOGICAL :: OK
        INTEGER :: I, J_INDEX, J
        DOUBLE PRECISION :: ROW_MAX, VALUE, DIAG
        OK = .TRUE.
        DO I = 0, SIZE - 1
            ! ROW MAX CALCULATION
            ROW_MAX = 0.D0
            DO J_INDEX = STARTS(I), STARTS(I + 1) - 1
                J = INDICES(J_INDEX)
                VALUE = VALUES(J_INDEX)
                IF ((I .NE. J) .AND. (VALUE .LT. 0.D0) .AND. (-VALUE .GT. ROW_MAX)) THEN
                    ROW_MAX = -VALUE
                END IF
                IF (I .EQ. J) THEN
                    DIAG = VALUE
                END IF
            END DO               
            IF (ROW_MAX .GT. DIAG) THEN
                WRITE (*,*) 'ROW_MAX = ', ROW_MAX, ', DIAG = ', DIAG
            END IF
            DO J_INDEX = STARTS(I), STARTS(I + 1) - 1
                J = INDICES(J_INDEX)
                VALUE = VALUES(J_INDEX)
                IF ((I .NE. J) .AND. (VALUE .GT. 0.D0) .AND. (VALUE .GE. EPS_STR * ROW_MAX)) THEN
                    WRITE (*,*) 'POSITIVE = ', VALUE, ', ROW_MAX = ', ROW_MAX, ', DIAG = ', DIAG
                    OK = .FALSE.                    
                END IF          
            END DO
        END DO        
    END FUNCTION
    
    SUBROUTINE AMG_GENERATE_S_CON_IMPL_(SIZE, LENGTH, STARTS, INDICES, VALUES, S_STARTS, S_INDICES, TEMP_STARTS, SPTR, FORCED_C, MIN_C)
        INTEGER, INTENT(IN) :: SIZE, LENGTH
        INTEGER, INTENT(OUT) :: SPTR
        INTEGER, DIMENSION(0:SIZE), INTENT(IN) :: STARTS
        INTEGER, DIMENSION(0:LENGTH - 1), INTENT(IN) :: INDICES
        INTEGER, DIMENSION(0:SIZE), INTENT(OUT) :: S_STARTS
        INTEGER, DIMENSION(0:LENGTH - 1), INTENT(OUT) :: S_INDICES
        ! NEW BUG INTENT(IN)
        DOUBLE PRECISION, DIMENSION(0:LENGTH - 1), INTENT(INOUT) :: VALUES
        INTEGER, DIMENSION(0:SIZE - 1), INTENT(INOUT) :: TEMP_STARTS        
        LOGICAL, DIMENSION(0:SIZE - 1), INTENT(INOUT) :: FORCED_C
        INTEGER, INTENT(INOUT) :: MIN_C
        INTEGER :: I, J_INDEX, J, DIAG_INDEX
        DOUBLE PRECISION :: ROW_MAX, VALUE, DIAG, ROW_SUM
        SPTR = 0
        DO I = 0, SIZE - 1
            ! ROW MAX CALCULATION
            ROW_MAX = 0.D0
            ROW_SUM = 0.D0
            DO J_INDEX = STARTS(I), STARTS(I + 1) - 1
                J = INDICES(J_INDEX)
                VALUE = VALUES(J_INDEX)
                IF (I .EQ. J) THEN
                    DIAG = VALUE
                    DIAG_INDEX = J_INDEX
                ELSE
                    ROW_SUM = ROW_SUM + ABS(VALUE)
                    IF ((VALUE .LT. 0.D0) .AND. (-VALUE .GT. ROW_MAX)) THEN
                        ROW_MAX = -VALUE
                    END IF
                END IF
            END DO                    
            ! NEW BUG
            IF ((DIAG .LE. 0.D0) .OR. (DIAG .LE. (ROW_MAX * 0.995D0)) .OR. (DIAG .LE. (ROW_SUM * 0.9D0))) THEN                            
                FORCED_C(I) = .TRUE.
                MIN_C = MIN_C + 1
            ELSE
                FORCED_C(I) = .FALSE.
            END IF
            !IF ((DIAG .LE. (ROW_MAX * 1.0000001)) .OR. (DIAG .LE. (ROW_SUM * 1.0000001))) THEN                
            !    VALUES(DIAG_INDEX) = 1.0000001 * MAX(ROW_MAX, ROW_SUM)
            !    FORCED_C(I) = .FALSE.
            !    !FORCED_C(I) = .TRUE.
            !    !MIN_C = MIN_C + 1
            !ELSE
            !    FORCED_C(I) = .FALSE.
            !END IF
            
            ! FILLING S_CON
            S_STARTS(I) = SPTR
            DO J_INDEX = STARTS(I), STARTS(I + 1) - 1
                J = INDICES(J_INDEX)
                VALUE = VALUES(J_INDEX)
                IF ((I .NE. J) .AND. (VALUE .LT. 0.D0) .AND. (-VALUE .GE. EPS_STR * ROW_MAX)) THEN
                    S_INDICES(SPTR) = J
                    SPTR = SPTR + 1
                    ! TEMP_STARTS IS CHANGING TO ROW SIZES
                    TEMP_STARTS(J) = TEMP_STARTS(J) + 1
                END IF                
            END DO
        END DO        
        S_STARTS(SIZE) = SPTR        
    END SUBROUTINE
    
    SUBROUTINE AMG_GENERATE_S_CON_(MAT_NO)
        INTEGER, INTENT(IN) :: MAT_NO
        INTEGER :: SIZE, LENGTH, SPTR, I, LEVEL_NO
        INTEGER, DIMENSION(:), ALLOCATABLE :: TEMP_STARTS
        SIZE = MATS(MAT_NO)%ROWS
        LENGTH = MATS(MAT_NO)%LENGTH
        ALLOCATE(MATS(MAT_NO)%S_STARTS(0:SIZE))
        ALLOCATE(MATS(MAT_NO)%S_INDICES(0:LENGTH - 1))
        ALLOCATE(TEMP_STARTS(0:SIZE - 1))
        !DEC$ SIMD
        DO I = 0, SIZE - 1
            TEMP_STARTS(I) = 0
        END DO
        LEVEL_NO = MAT_NO
        ALLOCATE(LEVELS(LEVEL_NO)%FORCED_C(0:SIZE - 1))
        LEVELS(LEVEL_NO)%MIN_C = 0
        CALL AMG_GENERATE_S_CON_IMPL_(SIZE, LENGTH, MATS(MAT_NO)%STARTS, MATS(MAT_NO)%INDICES, MATS(MAT_NO)%VALUES, MATS(MAT_NO)%S_STARTS, MATS(MAT_NO)%S_INDICES, TEMP_STARTS, SPTR, LEVELS(LEVEL_NO)%FORCED_C, LEVELS(LEVEL_NO)%MIN_C)                        
        CALL AMG_ALLOCATE_ST_CON_(MAT_NO, SPTR, SIZE)
        MATS(MAT_NO)%S_LENGTH = SPTR
        CALL AMG_GENERATE_ST_CON_IMPL_(SIZE, SPTR, TEMP_STARTS, MATS(MAT_NO)%S_STARTS, MATS(MAT_NO)%S_INDICES, MATS(MAT_NO)%ST_STARTS, MATS(MAT_NO)%ST_INDICES)
        DEALLOCATE(TEMP_STARTS)
    END SUBROUTINE
    
    !******************************************   COARSE SET GENERATORS PART   ******************************************!
    
    FUNCTION AMG_CALC_POTENTIAL_(ID, SIZE, LENGTH, ST_STARTS, ST_INDICES, IS_U, IS_F) RESULT(ANS)
        INTEGER, INTENT(IN) :: ID, SIZE, LENGTH        
        INTEGER, DIMENSION(0:SIZE), INTENT(IN) :: ST_STARTS
        INTEGER, DIMENSION(0:LENGTH - 1), INTENT(IN) :: ST_INDICES
        LOGICAL, DIMENSION(0:SIZE - 1), INTENT(IN) :: IS_U, IS_F
        INTEGER :: ANS
        INTEGER :: J_INDEX, J
        ANS = 0
        DO J_INDEX = ST_STARTS(ID), ST_STARTS(ID + 1) - 1
            J = ST_INDICES(J_INDEX)
            IF (IS_U(J)) THEN
                ANS = ANS + 1
            ELSEIF (IS_F(J)) THEN
                ANS = ANS + 2                
            END IF
        END DO
    END FUNCTION
    
    SUBROUTINE AMG_INIT_COARSE_SET_GENERATOR_IMPL_(SIZE, LENGTH, STARTS, IS_C, IS_F, IS_U, ST_STARTS, ST_INDICES, FORCED_C)
        INTEGER, INTENT(IN) :: SIZE, LENGTH
        INTEGER, DIMENSION(0:SIZE), INTENT(IN) :: STARTS, ST_STARTS
        INTEGER, DIMENSION(0:LENGTH - 1), INTENT(IN) :: ST_INDICES
        LOGICAL, DIMENSION(0:SIZE - 1), INTENT(OUT) :: IS_C, IS_F, IS_U
        LOGICAL, DIMENSION(0:SIZE - 1), INTENT(IN) :: FORCED_C
        LOGICAL :: FLAG, NOT_FLAG
        INTEGER :: I, VALUE
        !$OMP PARALLEL SHARED(IS_F, IS_U, IS_C, STARTS)
            !$OMP DO SCHEDULE(STATIC) PRIVATE(I, FLAG, NOT_FLAG)
        DO I = 0, SIZE - 1            
            IF (FORCED_C(I)) THEN
                IS_C(I) = .TRUE.
                IS_F(I) = .FALSE.
                IS_U(I) = .FALSE.
            ELSE
                IS_C(I) = .FALSE.
                FLAG = (STARTS(I + 1) - STARTS(I)) .GT. 1
                NOT_FLAG = .NOT. FLAG
                IS_F(I) = NOT_FLAG
                IS_U(I) = FLAG
            END IF
        END DO
            !$OMP END DO
        !$OMP END PARALLEL
        CALL RB_INIT_(SIZE)
        DO I = 0, SIZE - 1
            IF (.NOT. FORCED_C(I)) THEN
                VALUE = AMG_CALC_POTENTIAL_(I, SIZE, LENGTH, ST_STARTS, ST_INDICES, IS_U, IS_F)
                CALL RB_INSERT_(I, VALUE)
            END IF
        END DO
    END SUBROUTINE        

    SUBROUTINE AMG_INIT_COARSE_SET_GENERATOR_(LEVEL_NO)
        INTEGER, INTENT(IN) :: LEVEL_NO                
        INTEGER :: SIZE, MAT_NO, I
        MAT_NO = LEVEL_NO
        SIZE = MATS(MAT_NO)%ROWS
        ALLOCATE(LEVELS(LEVEL_NO)%IS_U(0:SIZE - 1))
        ALLOCATE(LEVELS(LEVEL_NO)%IS_F(0:SIZE - 1))
        ALLOCATE(LEVELS(LEVEL_NO)%IS_C(0:SIZE - 1))
        LEVELS(LEVEL_NO)%F_SIZE = 0
        CALL AMG_INIT_COARSE_SET_GENERATOR_IMPL_(SIZE, MATS(MAT_NO)%S_LENGTH, MATS(MAT_NO)%STARTS, LEVELS(LEVEL_NO)%IS_C, LEVELS(LEVEL_NO)%IS_F, LEVELS(LEVEL_NO)%IS_U, MATS(MAT_NO)%ST_STARTS, MATS(MAT_NO)%ST_INDICES, LEVELS(LEVEL_NO)%FORCED_C)        
    END SUBROUTINE
    
    SUBROUTINE AMG_FIRST_PATH_COARSEN_(FCOUNT, SIZE, LENGTH, ST_STARTS, ST_INDICES, S_STARTS, S_INDICES, IS_C, IS_F, IS_U)
        INTEGER, INTENT(INOUT) :: FCOUNT
        INTEGER, INTENT(IN) :: SIZE, LENGTH        
        INTEGER, DIMENSION(0:SIZE), INTENT(IN) :: ST_STARTS, S_STARTS
        INTEGER, DIMENSION(0:LENGTH - 1), INTENT(IN) :: ST_INDICES, S_INDICES
        LOGICAL, DIMENSION(0:SIZE - 1), INTENT(INOUT) :: IS_C, IS_F, IS_U
        INTEGER :: J, J_INDEX, K, K_INDEX, VALUE, ID
        DO WHILE (RB_SIZE .GT. 0)
            ID = RB_TREE_MAX_()            
            IF (POTENTIAL(ID) .LE. 0) THEN
                EXIT
            END IF
            IS_U(ID) = .FALSE.
            IS_F(ID) = .FALSE.
            IS_C(ID) = .TRUE.
            CALL RB_DELETE_(ID)
            POTENTIAL(ID) = -1
            CALL LIST_CLEAR_()
            DO J_INDEX = ST_STARTS(ID), ST_STARTS(ID + 1) - 1
                J = ST_INDICES(J_INDEX)
                IF (IS_U(J)) THEN
                    IS_F(J) = .TRUE.
                    IS_C(J) = .FALSE.
                    FCOUNT = FCOUNT + 1
                    IS_U(J) = .FALSE.
                    CALL RB_DELETE_(J)
                    POTENTIAL(J) = -1
                    CALL LIST_INSERT_(J)
                END IF
            END DO
            DO J_INDEX = 0, LIST_SIZE - 1
                J = LIST(J_INDEX)
                DO K_INDEX = S_STARTS(J), S_STARTS(J + 1) - 1
                    K = S_INDICES(K_INDEX)
                    IF (IS_U(K)) THEN
                        VALUE = AMG_CALC_POTENTIAL_(K, SIZE, LENGTH, ST_STARTS, ST_INDICES, IS_U, IS_F)
                        CALL RB_DELETE_(K)
                        ! THE FOLLOWING WILL AUTOMATICALLY BE DONE
                        ! POTENTIAL(K) = VALUE
                        CALL RB_INSERT_(K, VALUE)
                    END IF
                END DO
            END DO
        END DO
        ID = ROOT
        DO WHILE (ID .NE. -1)
            CALL RB_DELETE_(ID)
            POTENTIAL(ID) = -1
            IS_U(ID) = .FALSE.
            IF (.NOT. IS_C(ID)) THEN
                IS_F(ID) = .TRUE.
                FCOUNT = FCOUNT + 1
            END IF
            ID = ROOT
        END DO        
    END SUBROUTINE

    SUBROUTINE AMG_FILL_CF_(SIZE, C_SIZE, F_SIZE, C, F, IS_F, CF_INDEX)
        INTEGER, INTENT(IN) :: SIZE, C_SIZE, F_SIZE
        INTEGER, DIMENSION(0:C_SIZE - 1), INTENT(OUT) :: C
        INTEGER, DIMENSION(0:F_SIZE - 1), INTENT(OUT) :: F
        LOGICAL, DIMENSION(0:SIZE - 1), INTENT(IN) :: IS_F
        INTEGER, DIMENSION(0:SIZE - 1), INTENT(OUT) :: CF_INDEX
        INTEGER :: I, CPTR, FPTR
        CPTR = 0
        FPTR = 0
        DO I = 0, SIZE - 1
            IF (IS_F(I)) THEN
                F(FPTR) = I
                CF_INDEX(I) = FPTR
                FPTR = FPTR + 1
            ELSE
                C(CPTR) = I
                CF_INDEX(I) = CPTR
                CPTR = CPTR + 1
            END IF
        END DO        
    END SUBROUTINE
    
    FUNCTION AMG_IS_VITAL_COARSE_CELL_(ID, SIZE, LENGTH, S_STARTS, S_INDICES, ST_STARTS, ST_INDICES, FORCED_C) RESULT(ANS)
        INTEGER, INTENT(IN) :: ID, SIZE, LENGTH
        INTEGER, DIMENSION(0:SIZE), INTENT(IN) :: S_STARTS, ST_STARTS
        INTEGER, DIMENSION(0:LENGTH - 1), INTENT(IN) :: S_INDICES, ST_INDICES
        LOGICAL, DIMENSION(0:SIZE - 1), INTENT(IN) :: FORCED_C
        LOGICAL :: ANS
        INTEGER :: J_INDEX, J, CON_SIZE, J_START
        ANS = .FALSE.        
        IF (FORCED_C(ID)) THEN
            ANS = .TRUE.
            RETURN
        END IF
        ! BIG BUG
        !RETURN
        DO J_INDEX = ST_STARTS(ID), ST_STARTS(ID + 1) - 1
            J = ST_INDICES(J_INDEX)
            J_START = S_STARTS(J)
            CON_SIZE = S_STARTS(J + 1) - J_START
            IF ((CON_SIZE .EQ. 1) .AND. (S_INDICES(J_START) .EQ. ID)) THEN
                ANS = .TRUE.
                RETURN
            END IF
        END DO
    END FUNCTION
    
    SUBROUTINE AMG_FILL_CON_USING_HEAP_(SIZE, LENGTH, ST_SIZES, S_INDICES, PTR)
        INTEGER, INTENT(IN) :: SIZE, LENGTH
        INTEGER, INTENT(INOUT) :: PTR
        INTEGER, DIMENSION(0:SIZE - 1), INTENT(INOUT) :: ST_SIZES
        INTEGER, DIMENSIOn(0:LENGTH - 1), INTENT(INOUT) :: S_INDICES
        INTEGER :: I, ID, VALUE
        ID = -1
        DO WHILE (HEAP_SIZE .GT. 0)
            VALUE = HEAP_EXTRACT_MIN_()
            IF (VALUE .NE. ID) THEN
                ID = VALUE
                S_INDICES(PTR) = ID
                PTR = PTR + 1
                ST_SIZES(ID) = ST_SIZES(ID) + 1
            END IF
        END DO
    END SUBROUTINE
    
    SUBROUTINE AMG_FILL_A1_CON_SET_(I, SIZE, LENGTH, S_STARTS, S_INDICES, IS_C)
        INTEGER, INTENT(IN) :: I, SIZE, LENGTH
        INTEGER, DIMENSION(0:SIZE), INTENT(IN) :: S_STARTS
        INTEGER, DIMENSION(0:LENGTH - 1), INTENT(IN) :: S_INDICES
        LOGICAL, DIMENSION(0:SIZE - 1), INTENT(IN) :: IS_C
        INTEGER :: J_INDEX, J, K_INDEX, K
        DO J_INDEX = S_STARTS(I), S_STARTS(I + 1) - 1
            J = S_INDICES(J_INDEX)
            IF (IS_C(J)) THEN
                CALL HEAP_INSERT_(J)
            END IF
            DO K_INDEX = S_STARTS(J), S_STARTS(J + 1) - 1
                K = S_INDICES(K_INDEX)
                IF ((K .NE. I) .AND. IS_C(K)) THEN
                    CALL HEAP_INSERT_(K)
                END IF
            END DO
        END DO
    END SUBROUTINE
    
    SUBROUTINE AMG_SECOND_PATH_COARSEN_(LEVEL_NO, SIZE, F_SIZE, IS_C, IS_F, IS_U)
        INTEGER, INTENT(IN) :: LEVEL_NO, SIZE
        INTEGER, INTENT(INOUT) :: F_SIZE
        LOGICAL, DIMENSION(0:SIZE - 1), INTENT(INOUT) :: IS_C, IS_F, IS_U
        INTEGER :: COARSENING_TYPE
        INTEGER, DIMENSION(:), ALLOCATABLE :: S_STARTS, S_INDICES, ST_SIZES, ST_STARTS, ST_INDICES
        LOGICAL, DIMENSION(:), ALLOCATABLE :: IS_IN_SECOND_PATH
        INTEGER :: PTR, MAT_NO, S_LENGTH, NEW_S_LENGTH, I, VALUE        
        MAT_NO = LEVEL_NO
        S_LENGTH = MATS(MAT_NO)%S_LENGTH
        COARSENING_TYPE = LEVELS(LEVEL_NO)%COARSENING_TYPE
        ALLOCATE(S_STARTS(0:SIZE))
        ALLOCATE(S_INDICES(0:S_LENGTH - 1))
        ALLOCATE(ST_SIZES(0:SIZE - 1))
        ALLOCATE(IS_IN_SECOND_PATH(0:SIZE - 1))
        !DEC$ SIMD
        DO I = 0, SIZE - 1
            ST_SIZES(I) = 0
        END DO
        PTR = 0        
        DO I = 0, SIZE - 1
            S_STARTS(I) = PTR
            IF (IS_C(I) .AND. (.NOT. AMG_IS_VITAL_COARSE_CELL_(I, SIZE, S_LENGTH, MATS(MAT_NO)%S_STARTS, MATS(MAT_NO)%S_INDICES, MATS(MAT_NO)%ST_STARTS, MATS(MAT_NO)%ST_INDICES, LEVELS(LEVEL_NO)%FORCED_C))) THEN
                IS_IN_SECOND_PATH(I) = .TRUE.
                CALL HEAP_CLEAR_()
                IF (COARSENING_TYPE .EQ. 1) THEN
                    CALL AMG_FILL_A1_CON_SET_(I, SIZE, S_LENGTH, MATS(MAT_NO)%S_STARTS, MATS(MAT_NO)%S_INDICES, IS_C)
                ELSEIF (COARSENING_TYPE .EQ. 2) THEN
                    WRITE (*,*) 'A2 COARSENING NOT SUPPORTED'
                    ! CALL AMG_FILL_A2_CON_SET_()
                END IF
                CALL AMG_FILL_CON_USING_HEAP_(SIZE, S_LENGTH, ST_SIZES, S_INDICES, PTR)                
            ELSE
                IS_IN_SECOND_PATH(I) = .FALSE.
            END IF
        END DO
        S_STARTS(I) = PTR
        NEW_S_LENGTH = PTR
        ALLOCATE(ST_STARTS(0:SIZE))
        ALLOCATE(ST_INDICES(0:NEW_S_LENGTH - 1))
        CALL AMG_GENERATE_ST_CON_IMPL_(SIZE, NEW_S_LENGTH, ST_SIZES, S_STARTS, S_INDICES, ST_STARTS, ST_INDICES)
        DO I = 0, SIZE - 1
            IF (IS_IN_SECOND_PATH(I)) THEN
                IS_U(I) = .TRUE.
                VALUE = AMG_CALC_POTENTIAL_(I, SIZE, NEW_S_LENGTH, ST_STARTS, ST_INDICES, IS_U, IS_F)
                CALL RB_INSERT_(I, VALUE)
            END IF
        END DO
        CALL AMG_FIRST_PATH_COARSEN_(F_SIZE, SIZE, NEW_S_LENGTH, ST_STARTS, ST_INDICES, S_STARTS, S_INDICES, IS_C, IS_F, IS_U)
        DEALLOCATE(S_STARTS)
        DEALLOCATE(S_INDICES)
        DEALLOCATE(ST_SIZES)
        DEALLOCATE(IS_IN_SECOND_PATH)
        DEALLOCATE(ST_STARTS)
        DEALLOCATE(ST_INDICES)
    END SUBROUTINE
    
    SUBROUTINE AMG_GENERATE_COARSE_SET_(LEVEL_NO)
        INTEGER, INTENT(IN) :: LEVEL_NO
        INTEGER :: MAT_NO, F_SIZE, C_SIZE, SIZE, FPTR, CPTR, I
        MAT_NO = LEVEL_NO
        SIZE = MATS(MAT_NO)%ROWS
        CALL AMG_GENERATE_S_CON_(MAT_NO)
        CALL AMG_INIT_COARSE_SET_GENERATOR_(LEVEL_NO)
        F_SIZE = LEVELS(LEVEL_NO)%F_SIZE        
        CALL AMG_FIRST_PATH_COARSEN_(F_SIZE, MATS(MAT_NO)%ROWS, MATS(MAT_NO)%S_LENGTH, MATS(MAT_NO)%ST_STARTS, MATS(MAT_NO)%ST_INDICES, MATS(MAT_NO)%S_STARTS, MATS(MAT_NO)%S_INDICES, LEVELS(LEVEL_NO)%IS_C, LEVELS(LEVEL_NO)%IS_F, LEVELS(LEVEL_NO)%IS_U)
        IF (LEVELS(LEVEL_NO)%COARSENING_TYPE .GT. 0) THEN
            CALL AMG_SECOND_PATH_COARSEN_(LEVEL_NO, SIZE, F_SIZE, LEVELS(LEVEL_NO)%IS_C, LEVELS(LEVEL_NO)%IS_F, LEVELS(LEVEL_NO)%IS_U)
        END IF
        C_SIZE = SIZE - F_SIZE        
        LEVELS(LEVEL_NO)%F_SIZE = F_SIZE
        LEVELS(LEVEL_NO)%C_SIZE  = C_SIZE
        CALL RB_FINALIZE_()
        ALLOCATE(LEVELS(LEVEL_NO)%F(0:F_SIZE - 1))
        ALLOCATE(LEVELS(LEVEL_NO)%C(0:C_SIZE - 1))
        ALLOCATE(LEVELS(LEVEL_NO)%CF_INDEX(0:SIZE - 1))
        CALL AMG_FILL_CF_(SIZE, C_SIZE, F_SIZE, LEVELS(LEVEL_NO)%C, LEVELS(LEVEL_NO)%F, LEVELS(LEVEL_NO)%IS_F, LEVELS(LEVEL_NO)%CF_INDEX)        
    END SUBROUTINE
    
    !******************************************   INTERPOLATION PART   ******************************************!

    SUBROUTINE AMG_GENERATE_INTERPOLATOR_(LEVEL_NO)
        INTEGER, INTENT(IN) :: LEVEL_NO
        INTEGER, DIMENSION (:), ALLOCATABLE:: IS_FSTAR
        INTEGER, DIMENSION(:), ALLOCATABLE :: NO_STRONG
        INTEGER :: CPTR, INTPL_COUNTER, N, CSIZE, I, FSIZE, DS, STRONG_SIZE, T, STRONG_TOTAL, NO_STRONG_SIZE, FPRIME_SIZE, FPRIME_PTR, NEW_FPRIME_SIZE,INTPL_PTR , PTR, NEW_PTR, INDEX, PTR_NEW
        INTEGER, DIMENSION(:), ALLOCATABLE :: FPRIME, NEW_FPRIME, INTERPOLATED, NEW_FPRIME2, NEW_FPRIME3
        INTEGER :: INDICES_LENGTH, SINDEX, SI
        TYPE(CRS_TYPE) :: IFC1, IFC2, IFC3, IFC4

        CSIZE= LEVELS(LEVEL_NO)%C_SIZE
        FSIZE= LEVELS(LEVEL_NO)%F_SIZE
        ALLOCATE(NO_STRONG(0: FSIZE-1))
        ALLOCATE(IS_FSTAR(0: MATS(LEVEL_NO)%ROWS-1))
        INTPL_COUNTER= 0
        IS_FSTAR= 0
        NO_STRONG= 0
        !ISFSTAR OF CORSE POINTS = 1
        DO I=0, CSIZE-1
            IS_FSTAR(LEVELS(LEVEL_NO)%C(I)) = 1
        END DO
        !EVALUATE FPRIME AND 1ST IFC SIZES
        NO_STRONG_SIZE= 0
        INDICES_LENGTH= 0
        DO T=0, FSIZE-1
            I= LEVELS(LEVEL_NO)%F(T)
            DS= MATS(LEVEL_NO)%S_STARTS(I+1)- MATS(LEVEL_NO)%S_STARTS(I)
            IF(DS==0) THEN
                NO_STRONG(T)= 1
                NO_STRONG_SIZE= NO_STRONG_SIZE+1 
            ELSE
                DO SINDEX= MATS(LEVEL_NO)%S_STARTS(I), MATS(LEVEL_NO)%S_STARTS(I+1)-1
                    SI=  MATS(LEVEL_NO)%S_INDICES(SINDEX)
                    IF(IS_FSTAR(SI) .EQ. 1) THEN
                        INDICES_LENGTH = INDICES_LENGTH +1
                    END IF    
                END DO
            END IF
        END DO
        FPRIME_SIZE= FSIZE-NO_STRONG_SIZE
      
        ALLOCATE(IFC1%STARTS(0: FSIZE+CSIZE))
        ALLOCATE(IFC2%STARTS(0: FSIZE+CSIZE))
        ALLOCATE(IFC3%STARTS(0: FSIZE+CSIZE))
        ALLOCATE(IFC4%STARTS(0: FSIZE+CSIZE))
        IFC1%ROWS= FSIZE+CSIZE
        IFC2%ROWS= FSIZE+CSIZE
        IFC3%ROWS= FSIZE+CSIZE
        IFC4%ROWS= FSIZE+CSIZE
        IFC1%STARTS= 0
        IFC2%STARTS= 0
        IFC3%STARTS= 0
        IFC4%STARTS= 0
        
        ALLOCATE(IFC1%INDICES(0: INDICES_LENGTH-1))
        ALLOCATE(IFC1%VALUES (0: INDICES_LENGTH-1))
        IFC1%LENGTH= INDICES_LENGTH
        ALLOCATE(FPRIME(0:FPRIME_SIZE-1))
        !REMOVE ALL WEAK CONNECTIONS FROM FPRIME
        FPRIME_PTR=0
        DO I= 0, FSIZE-1
            IF(NO_STRONG(I) .NE. 1) THEN
                FPRIME(FPRIME_PTR)= LEVELS(LEVEL_NO)%F(I) 
                FPRIME_PTR= FPRIME_PTR+1
            END IF
        END DO
        ALLOCATE(INTERPOLATED(0:FPRIME_SIZE-1))
        INTPL_PTR =0
        
        !1ST PASS
        !WRITE(*,*) "FIRST PASS OF INTERPOLATION"
        CALL INTERPOLATE_(FPRIME_SIZE, MATS(LEVEL_NO)%ROWS, MATS(LEVEL_NO)%LENGTH, MATS(LEVEL_NO)%S_LENGTH, MATS(LEVEL_NO)%ROWS, FPRIME, MATS(LEVEL_NO)%MAX_CONNECTION, IS_FSTAR, MATS(LEVEL_NO)%STARTS, MATS(LEVEL_NO)%INDICES, MATS(LEVEL_NO)%VALUES, MATS(LEVEL_NO)%S_STARTS, MATS(LEVEL_NO)%S_INDICES, INTERPOLATED, INTPL_PTR, IFC1%STARTS, IFC1%INDICES, IFC1%VALUES, IFC1%STARTS, IFC1%INDICES, IFC1%VALUES, IFC1%ROWS, IFC1%LENGTH, IFC1%ROWS, IFC1%LENGTH, LEVELS(LEVEL_NO)%CF_INDEX , 0)       
        NEW_FPRIME_SIZE = FPRIME_SIZE- INTPL_PTR
        IF(NEW_FPRIME_SIZE .NE. 0) THEN
            !WRITE(*,*) "SECOND PASS OF INTERPOLATION"
            ALLOCATE(NEW_FPRIME(0:NEW_FPRIME_SIZE-1))
            CALL DETERMINE_NEW_FPRIME_(FPRIME_SIZE, NEW_FPRIME_SIZE, INTPL_PTR, CSIZE+FSIZE, INTERPOLATED, FPRIME, IS_FSTAR, NEW_FPRIME)
            DEALLOCATE(INTERPOLATED)
            DEALLOCATE(FPRIME)
            INTPL_PTR = 0
            ALLOCATE(INTERPOLATED(0:NEW_FPRIME_SIZE-1))
                                    
            INDICES_LENGTH=  EVALUATE_IFC_SIZES_(NEW_FPRIME_SIZE, MATS(LEVEL_NO)%ROWS, MATS(LEVEL_NO)%LENGTH ,MATS(LEVEL_NO)%ROWS, NEW_FPRIME, MATS(LEVEL_NO)%S_STARTS, MATS(LEVEL_NO)%S_INDICES, IS_FSTAR, IFC1%ROWS, IFC1%STARTS)
            ALLOCATE(IFC2%INDICES(0: INDICES_LENGTH-1))
            ALLOCATE(IFC2%VALUES (0: INDICES_LENGTH-1))
            IFC2%LENGTH= INDICES_LENGTH
            CALL INTERPOLATE_(NEW_FPRIME_SIZE, MATS(LEVEL_NO)%ROWS, MATS(LEVEL_NO)%LENGTH, MATS(LEVEL_NO)%S_LENGTH, MATS(LEVEL_NO)%ROWS, NEW_FPRIME, MATS(LEVEL_NO)%MAX_CONNECTION, IS_FSTAR, MATS(LEVEL_NO)%STARTS, MATS(LEVEL_NO)%INDICES, MATS(LEVEL_NO)%VALUES, MATS(LEVEL_NO)%S_STARTS, MATS(LEVEL_NO)%S_INDICES, INTERPOLATED, INTPL_PTR, IFC1%STARTS, IFC1%INDICES, IFC1%VALUES, IFC2%STARTS, IFC2%INDICES, IFC2%VALUES, IFC1%ROWS, IFC1%LENGTH, IFC2%ROWS, IFC2%LENGTH, LEVELS(LEVEL_NO)%CF_INDEX , 1)              
            FPRIME_SIZE= NEW_FPRIME_SIZE
            NEW_FPRIME_SIZE = NEW_FPRIME_SIZE- INTPL_PTR
            IF(NEW_FPRIME_SIZE .NE. 0) THEN
                !WRITE(*,*) "THIRD PASS OF INTERPOLATION"
                ALLOCATE(NEW_FPRIME2(0:NEW_FPRIME_SIZE-1))
                CALL DETERMINE_NEW_FPRIME_(FPRIME_SIZE, NEW_FPRIME_SIZE, INTPL_PTR, CSIZE+FSIZE, INTERPOLATED, NEW_FPRIME, IS_FSTAR, NEW_FPRIME2)
                DEALLOCATE(INTERPOLATED)
                DEALLOCATE (NEW_FPRIME)
                INTPL_PTR = 0
                ALLOCATE(INTERPOLATED(0:NEW_FPRIME_SIZE-1))
                                        
                INDICES_LENGTH=  EVALUATE_IFC_SIZES_(NEW_FPRIME_SIZE, MATS(LEVEL_NO)%ROWS, MATS(LEVEL_NO)%LENGTH ,MATS(LEVEL_NO)%ROWS, NEW_FPRIME2, MATS(LEVEL_NO)%S_STARTS, MATS(LEVEL_NO)%S_INDICES, IS_FSTAR, IFC2%ROWS, IFC2%STARTS)
                ALLOCATE(IFC3%INDICES(0: INDICES_LENGTH-1))
                ALLOCATE(IFC3%VALUES (0: INDICES_LENGTH-1))
                IFC3%LENGTH= INDICES_LENGTH
                CALL INTERPOLATE_(NEW_FPRIME_SIZE, MATS(LEVEL_NO)%ROWS, MATS(LEVEL_NO)%LENGTH, MATS(LEVEL_NO)%S_LENGTH, MATS(LEVEL_NO)%ROWS, NEW_FPRIME2, MATS(LEVEL_NO)%MAX_CONNECTION, IS_FSTAR, MATS(LEVEL_NO)%STARTS, MATS(LEVEL_NO)%INDICES, MATS(LEVEL_NO)%VALUES, MATS(LEVEL_NO)%S_STARTS, MATS(LEVEL_NO)%S_INDICES, INTERPOLATED, INTPL_PTR, IFC2%STARTS, IFC2%INDICES, IFC2%VALUES, IFC3%STARTS, IFC3%INDICES, IFC3%VALUES, IFC2%ROWS, IFC2%LENGTH, IFC3%ROWS, IFC3%LENGTH, LEVELS(LEVEL_NO)%CF_INDEX , 2)              
                
                FPRIME_SIZE= NEW_FPRIME_SIZE
                NEW_FPRIME_SIZE = NEW_FPRIME_SIZE- INTPL_PTR
                IF(NEW_FPRIME_SIZE .NE. 0) THEN
                    !WRITE(*,*) "FORTH PASS OF INTERPOLATION"
                    ALLOCATE(NEW_FPRIME3(0:NEW_FPRIME_SIZE-1))
                    CALL DETERMINE_NEW_FPRIME_(FPRIME_SIZE, NEW_FPRIME_SIZE, INTPL_PTR, CSIZE+FSIZE, INTERPOLATED, NEW_FPRIME2, IS_FSTAR, NEW_FPRIME3)
                    DEALLOCATE(INTERPOLATED)
                    DEALLOCATE (NEW_FPRIME2)
                    INTPL_PTR = 0
                    ALLOCATE(INTERPOLATED(0:NEW_FPRIME_SIZE-1))
                    INDICES_LENGTH=  EVALUATE_IFC_SIZES_(NEW_FPRIME_SIZE, MATS(LEVEL_NO)%ROWS, MATS(LEVEL_NO)%LENGTH ,MATS(LEVEL_NO)%ROWS, NEW_FPRIME3, MATS(LEVEL_NO)%S_STARTS, MATS(LEVEL_NO)%S_INDICES, IS_FSTAR, IFC3%ROWS, IFC3%STARTS)
                    ALLOCATE(IFC4%INDICES(0: INDICES_LENGTH-1))
                    ALLOCATE(IFC4%VALUES (0: INDICES_LENGTH-1))
                    IFC4%LENGTH= INDICES_LENGTH
                    CALL INTERPOLATE_(NEW_FPRIME_SIZE, MATS(LEVEL_NO)%ROWS, MATS(LEVEL_NO)%LENGTH, MATS(LEVEL_NO)%S_LENGTH, MATS(LEVEL_NO)%ROWS, NEW_FPRIME3, MATS(LEVEL_NO)%MAX_CONNECTION, IS_FSTAR, MATS(LEVEL_NO)%STARTS, MATS(LEVEL_NO)%INDICES, MATS(LEVEL_NO)%VALUES, MATS(LEVEL_NO)%S_STARTS, MATS(LEVEL_NO)%S_INDICES, INTERPOLATED, INTPL_PTR, IFC3%STARTS, IFC3%INDICES, IFC3%VALUES, IFC4%STARTS, IFC4%INDICES, IFC4%VALUES, IFC3%ROWS, IFC3%LENGTH, IFC4%ROWS, IFC4%LENGTH, LEVELS(LEVEL_NO)%CF_INDEX , 3)              
                    FPRIME_SIZE= NEW_FPRIME_SIZE
                    NEW_FPRIME_SIZE = NEW_FPRIME_SIZE- INTPL_PTR
                    IF(NEW_FPRIME_SIZE .NE. 0) THEN
                        WRITE(*,*) "INTERPLATION PASSES ARE MORE THAN 4!"
                        WRITE(*,*) "number of remaind Fs!" , NEW_FPRIME_SIZE
                       !if(LEVELS(LEVEL_NO)%IS_C(3363)) THEN  
                       !     write(*,*) '3363 IS C'
                       ! end if
                       ! if(LEVELS(LEVEL_NO)%IS_C(4603)) THEN  
                       !     write(*,*) '4603 IS C'
                       ! end if
                       ! if(LEVELS(LEVEL_NO)%IS_C(3523)) THEN  
                       !     write(*,*) '3523 IS C'
                       ! end if
                       ! if(LEVELS(LEVEL_NO)%IS_C(3483)) THEN  
                       !     write(*,*) '3483 IS C'
                       ! end if
                       ! if(LEVELS(LEVEL_NO)%IS_C(4643)) THEN  
                       !     write(*,*) '4643 IS C'
                       ! end if
                        WRITE(*,*) "LEVEL_NO=" , LEVEL_NO
                       !WRITE(*,*) "LEVEL_NO=" , MAX_LEVEL
                       !WRITE(*,*) "fprime" , NEW_FPRIME2
                       !WRITE(*,*) "interplated" , INTERPOLATED
                       PAUSE
                    END IF
                    DEALLOCATE (NEW_FPRIME3)
                END IF
            END IF
        END IF 
        LEVELS(LEVEL_NO)%IFC%ROWS= IFC1%ROWS
        LEVELS(LEVEL_NO)%IFC%LENGTH= IFC1%LENGTH + IFC2%LENGTH + IFC3%LENGTH + IFC4%LENGTH + CSIZE
        ALLOCATE (LEVELS(LEVEL_NO)%IFC%STARTS(0: LEVELS(LEVEL_NO)%IFC%ROWS))
        ALLOCATE (LEVELS(LEVEL_NO)%IFC%INDICES(0: LEVELS(LEVEL_NO)%IFC%LENGTH-1))
        ALLOCATE (LEVELS(LEVEL_NO)%IFC%VALUES(0: LEVELS(LEVEL_NO)%IFC%LENGTH-1))
        CALL MERGE_IFCS_(LEVELS(LEVEL_NO)%IFC%ROWS, LEVELS(LEVEL_NO)%IFC%LENGTH, LEVELS(LEVEL_NO)%IFC%STARTS, LEVELS(LEVEL_NO)%IFC%INDICES, LEVELS(LEVEL_NO)%IFC%VALUES, IFC1%LENGTH, IFC1%STARTS, IFC1%INDICES, IFC1%VALUES, IFC2%LENGTH, IFC2%STARTS, IFC2%INDICES, IFC2%VALUES, IFC3%LENGTH, IFC3%STARTS, IFC3%INDICES, IFC3%VALUES, IFC4%LENGTH, IFC4%STARTS, IFC4%INDICES, IFC4%VALUES, CSIZE, LEVELS(LEVEL_NO)%C, LEVELS(LEVEL_NO)%IS_C)
        LEVELS(LEVEL_NO)%IFCT%ROWS= CSIZE
        LEVELS(LEVEL_NO)%IFCT%LENGTH= LEVELS(LEVEL_NO)%IFC%LENGTH 
        ALLOCATE (LEVELS(LEVEL_NO)%IFCT%STARTS(0: LEVELS(LEVEL_NO)%IFCT%ROWS))
        ALLOCATE (LEVELS(LEVEL_NO)%IFCT%INDICES(0: LEVELS(LEVEL_NO)%IFCT%LENGTH-1))
        ALLOCATE (LEVELS(LEVEL_NO)%IFCT%VALUES(0: LEVELS(LEVEL_NO)%IFCT%LENGTH-1))
        CALL CRS_TRANSPOSE_(LEVELS(LEVEL_NO)%IFC%ROWS, CSIZE, LEVELS(LEVEL_NO)%IFC%LENGTH, LEVELS(LEVEL_NO)%IFC%STARTS, LEVELS(LEVEL_NO)%IFC%INDICES, LEVELS(LEVEL_NO)%IFC%VALUES, LEVELS(LEVEL_NO)%IFCT%STARTS, LEVELS(LEVEL_NO)%IFCT%INDICES, LEVELS(LEVEL_NO)%IFCT%VALUES)
        DEALLOCATE (IFC1%STARTS)
        DEALLOCATE (IFC1%INDICES)
        DEALLOCATE (IFC1%VALUES)
        IF(ALLOCATED(IFC2%INDICES)) THEN
            DEALLOCATE (IFC2%STARTS)
            DEALLOCATE (IFC2%INDICES)
            DEALLOCATE (IFC2%VALUES)
        END IF
        IF(ALLOCATED(IFC3%INDICES)) THEN
            DEALLOCATE (IFC3%STARTS)
            DEALLOCATE (IFC3%INDICES)
            DEALLOCATE (IFC3%VALUES)
        END IF
        IF(ALLOCATED(IFC4%INDICES)) THEN
            DEALLOCATE (IFC4%STARTS)
            DEALLOCATE (IFC4%INDICES)
            DEALLOCATE (IFC4%VALUES)
        END IF
        DEALLOCATE(IS_FSTAR)
        DEALLOCATE(NO_STRONG)
        DEALLOCATE(INTERPOLATED)
        !write(*,*) "INTERPOLATION FINISHED"
        !call crs_print_(LEVELS(LEVEL_NO)%IFC%ROWS, LEVELS(LEVEL_NO)%IFC%length, LEVELS(LEVEL_NO)%IFC%starts, LEVELS(LEVEL_NO)%IFC%indices, LEVELS(LEVEL_NO)%IFC%values)
    END SUBROUTINE
     
    FUNCTION EVALUATE_IFC_SIZES_(FPRIME_SIZE, SSTART_SIZE, SINDICES_SIZE, ISFSTAR_SIZE, FPRIME, SSTARTS, SINDICES, IS_FSTAR, IFC_ROWSIZE, IFC_START) RESULT(INDICES_LENGTH)
        INTEGER, INTENT(IN) :: FPRIME_SIZE, SSTART_SIZE, SINDICES_SIZE, ISFSTAR_SIZE, IFC_ROWSIZE
        INTEGER :: INDICES_LENGTH, T, I, DS, SINDEX, SI
        INTEGER, DIMENSION(0:FPRIME_SIZE-1), INTENT(IN) :: FPRIME
        INTEGER, DIMENSION(0:ISFSTAR_SIZE-1), INTENT(IN) :: IS_FSTAR
        INTEGER, DIMENSION(0:SSTART_SIZE), INTENT(IN) :: SSTARTS
        INTEGER, DIMENSION(0:SINDICES_SIZE-1), INTENT(IN) :: SINDICES
        INTEGER, DIMENSION(0:IFC_ROWSIZE), INTENT(IN) :: IFC_START
        INDICES_LENGTH= 0
        DO T=0, FPRIME_SIZE-1
            I= FPRIME(T)
            DS= SSTARTS(I+1)- SSTARTS(I)
            DO SINDEX= SSTARTS(I), SSTARTS(I+1)-1
                SI=  SINDICES(SINDEX)
                IF(IS_FSTAR(SI) .EQ. 1) THEN
                    INDICES_LENGTH = INDICES_LENGTH + IFC_START(SI+1) - IFC_START(SI)
                END IF    
            END DO
        END DO
    END FUNCTION
    
    
    SUBROUTINE DETERMINE_NEW_FPRIME_(FPRIME_SIZE, NEW_FPRIME_SIZE, INTPL_PTR, ISFSTAR_SIZE, INTERPOLATED, FPRIME, IS_FSTAR, NEW_FPRIME)
        INTEGER, INTENT(IN) :: FPRIME_SIZE, NEW_FPRIME_SIZE, INTPL_PTR, ISFSTAR_SIZE
        INTEGER, DIMENSION(0:FPRIME_SIZE-1), INTENT(IN) :: INTERPOLATED
        INTEGER, DIMENSION(0:FPRIME_SIZE-1), INTENT(IN) :: FPRIME
        INTEGER, DIMENSION(0:ISFSTAR_SIZE-1), INTENT(OUT) :: IS_FSTAR
        INTEGER, DIMENSION(0:NEW_FPRIME_SIZE-1), INTENT(OUT) :: NEW_FPRIME
        INTEGER :: PTR_NEW, PTR, INDEX, T, I
        DO INDEX= 0, INTPL_PTR-1
            T= INTERPOLATED(INDEX)
            IS_FSTAR(T) =1
        END DO
        PTR_NEW =0
        PTR=0
        DO I= 0, FPRIME_SIZE-1
            IF(INTERPOLATED(PTR) .NE. FPRIME(I)) THEN 
                NEW_FPRIME(PTR_NEW) = FPRIME(I) 
                PTR_NEW = PTR_NEW +1
            ELSE
                PTR = PTR +1
            END IF
        END DO
    END SUBROUTINE
    
    SUBROUTINE INTERPOLATE_(FPRIME_SIZE, ROW, LENGTH, S_LENGTH, SIZE, FPRIME, MAX_CONNECTION, IS_FSTAR, STARTS, INDICES, VALUES, S_STARTS, S_INDICES, INTERPOLATED, INTPL_PTR, IFC_STARTS, IFC_INDICES, IFC_VALUES, NEW_STARTS, NEW_INDICES, NEW_VALUES , START_SIZE, INDICES_SIZE, NEW_START_SIZE, NEW_INDICES_SIZE, CF_INDEX, COUNTER)     
        INTEGER, INTENT(IN) :: FPRIME_SIZE, MAX_CONNECTION, ROW, LENGTH, S_LENGTH, COUNTER, SIZE, START_SIZE, INDICES_SIZE, NEW_START_SIZE, NEW_INDICES_SIZE
        INTEGER, DIMENSION(0:FPRIME_SIZE-1) , INTENT(IN) :: FPRIME
        INTEGER, DIMENSION(0: ROW-1) , INTENT(IN) :: IS_FSTAR
        INTEGER, DIMENSION(0:ROW) , INTENT(IN) :: STARTS
        INTEGER, DIMENSION(0:LENGTH-1) , INTENT(IN) :: INDICES
        DOUBLE PRECISION, DIMENSION(0:LENGTH-1) , INTENT(IN) :: VALUES
        INTEGER, DIMENSION(0:ROW) , INTENT(IN) :: S_STARTS
        INTEGER, DIMENSION(0:S_LENGTH-1) , INTENT(IN) :: S_INDICES
        INTEGER, DIMENSION(0:FPRIME_SIZE-1) :: INTERPOLATED
        INTEGER, DIMENSION(0:SIZE) , INTENT(IN) :: CF_INDEX
        INTEGER, DIMENSION(0:START_SIZE), INTENT(IN) :: IFC_STARTS
        INTEGER, DIMENSION(0:INDICES_SIZE-1), INTENT(IN) :: IFC_INDICES
        DOUBLE PRECISION, DIMENSION(0:INDICES_SIZE-1), INTENT(IN) :: IFC_VALUES
        INTEGER, DIMENSION(0:NEW_START_SIZE), INTENT(INOUT) :: NEW_STARTS
        INTEGER, DIMENSION(0:NEW_INDICES_SIZE-1), INTENT(INOUT) :: NEW_INDICES
        DOUBLE PRECISION, DIMENSION(0:NEW_INDICES_SIZE-1), INTENT(INOUT) :: NEW_VALUES
        
        INTEGER, DIMENSION(0:MAX_CONNECTION-1) :: PI
        DOUBLE PRECISION, DIMENSION(0:MAX_CONNECTION - 1) :: PI_VALS
        INTEGER, DIMENSION(0:MAX_CONNECTION-1) :: PNTR
        INTEGER, DIMENSION(0:MAX_CONNECTION-1) :: ROW_END
        LOGICAL, DIMENSION(0:MAX_CONNECTION-1) :: MEMBER
        DOUBLE PRECISION, DIMENSION(0:MAX_CONNECTION-1) :: MULT
        INTEGER :: PI_PTR, INTPL_PTR, T, I, SINDEX, SI, PTR, JINDEX, J
        INTEGER :: P, COL, KK, REMAIND_ROW, MINI, MIN_COL, PINDEX_MIN, K, IFCI_PTR, IFCS_PTR, II, IPREV
        DOUBLE PRECISION :: SUM_AIJ_POS , SUM_AIJ_NEG, SUM_AIK_POS, SUM_AIK_NEG, AIJ, AIK, AII, ALPHAI, BETAI, M
        INTERPOLATED= -1
        NEW_VALUES = 0.D0
        NEW_STARTS(0)= 0
        IFCS_PTR=1
        IFCI_PTR=0
        IPREV= 0
        DO T= 0, FPRIME_SIZE -1
            PI= -1
            PNTR=-1
            ROW_END=-1
            PI_PTR=0 
            I= FPRIME(T)
            ! NEWLY ADDED
            SINDEX = S_STARTS(I)
            DO JINDEX = STARTS(I), STARTS(I + 1) - 1
                IF (SINDEX .EQ. S_STARTS(I + 1)) THEN
                    EXIT
                END IF
                J = INDICES(JINDEX)
                SI = S_INDICES(SINDEX)
                IF (SI .EQ. J) THEN                    
                    IF (IS_FSTAR(SI) .EQ. 1) THEN
                        PI(PI_PTR) = SI
                        PI_VALS(PI_PTR) = VALUES(JINDEX)
                        PI_PTR = PI_PTR + 1
                    END IF
                    SINDEX = SINDEX + 1
                END IF
            END DO
            ! TO BE REMOVED
            !DO SINDEX= S_STARTS(I), S_STARTS(I+1)-1
            !    SI=  S_INDICES(SINDEX)
            !    IF(IS_FSTAR(SI) .EQ. 1) THEN
            !        PI(PI_PTR) = SI
            !        PI_PTR = PI_PTR +1
            !    END IF    
            !END DO
            IF(PI_PTR .NE. 0) THEN
                DO II = IPREV+1, I-1
                    NEW_STARTS(II)= IFCI_PTR
                END DO
                NEW_STARTS(II)= IFCI_PTR
                IPREV= I
                INTERPOLATED(INTPL_PTR)= I
                !if(I .eq. 2163 .OR. I .EQ. 3403 .OR. I .EQ. 2323 .OR. I .EQ. 3483 .OR. I .EQ. 4643) then
                !    write(*,*) 'I=', I
                !    write(*,*) 'counter= ', counter
                !    write(*,*) 'PI=' , PI(0:pi_ptr)
                !end if
                !if(I .eq. 963 .OR. I .EQ. 2203 .OR. I .EQ. 1123 .OR. I .EQ. 2283 .OR. I .EQ. 3443) then
                !    write(*,*) 'I=', I
                !    write(*,*) 'counter= ', counter
                !    write(*,*) 'PI=' , PI(0:pi_ptr)
                !end if
                !if(I .eq. 1003 .OR. I .EQ. 1083 .OR. I .EQ. 2243) then
                !    write(*,*) 'I=', I
                !    write(*,*) 'counter= ', counter
                !    write(*,*) 'PI=' , PI(0:pi_ptr)
                !end if
                INTPL_PTR= INTPL_PTR +1
                SUM_AIJ_POS = 0;
                SUM_AIJ_NEG = 0;
                SUM_AIK_POS = 0;
                SUM_AIK_NEG = 0;
                PTR=0
                DO JINDEX= STARTS(I), STARTS(I+1) -1 
                    J= INDICES(JINDEX);
                    AIJ= VALUES(JINDEX)
                    IF(J .NE. I) THEN
                        SUM_AIJ_POS = SUM_AIJ_POS + MAX(AIJ, 0.D0);
                        SUM_AIJ_NEG = SUM_AIJ_NEG + MIN(AIJ, 0.D0);
                        IF(J .EQ. PI(PTR)) THEN
                            SUM_AIK_POS = SUM_AIK_POS + MAX(AIJ, 0.D0);
                            SUM_AIK_NEG = SUM_AIK_NEG + MIN(AIJ, 0.D0);
                            PTR= PTR+ 1
                        END IF
                    ELSE
                        AII= AIJ
                    END IF
                END DO
                ALPHAI= SUM_AIJ_NEG/ SUM_AIK_NEG;
                BETAI= SUM_AIJ_POS/ SUM_AIK_POS;
                PTR=0
                IF(COUNTER .EQ. 0) THEN
                    ! NEWLY ADDED
                    DO P = 0, PI_PTR - 1
                        K = PI(P)
                        AIK = PI_VALS(P)
                        IF(AIK .LE. 0.D0) THEN
                            M = -ALPHAI * AIK / AII
                        ELSE
                            M = -BETAI * AIK / AII
                        END IF
                        NEW_VALUES(IFCI_PTR) = M
                        NEW_INDICES(IFCI_PTR) = CF_INDEX(K)
                           
                        IFCI_PTR = IFCI_PTR + 1
                        PTR = PTR + 1
                    END DO
                    ! TO BE REMOVED
                    !DO JINDEX= STARTS(I), STARTS(I+1) -1 
                    !    J= INDICES(JINDEX);
                    !    K=PI(PTR)
                    !    IF(J .EQ. K) THEN
                    !        AIK= VALUES(JINDEX)
                    !        IF(AIK <= 0.D0) THEN
                    !            M= -ALPHAI*AIK /AII
                    !        ELSE
                    !            M= -BETAI *AIK /AII
                    !        END IF
                    !        NEW_VALUES(IFCI_PTR)= M
                    !        NEW_INDICES(IFCI_PTR)= CF_INDEX(K)
                    !       
                    !        IFCI_PTR = IFCI_PTR+1
                    !        PTR= PTR+ 1
                    !    END IF
                    !    IF(PTR .EQ. PI_PTR) THEN
                    !        EXIT
                    !    END IF
                    !END DO
                ELSE
                    DO P=0, PI_PTR-1
                        K= PI(P)
                        PNTR(P)= IFC_STARTS(K)
                        ROW_END(P)= IFC_STARTS(K+1)-1
                        MEMBER(P)= .TRUE.
                        COL= IFC_INDICES(PNTR(P))
                        ! NEWLY ADDED
                        AIK = PI_VALS(P)
                        ! TO BE REMOVED
                        !DO JINDEX= STARTS(I), STARTS(I+1)-1
                        !    IF(INDICES(JINDEX) .EQ. K) THEN
                        !        AIK= VALUES(JINDEX)
                        !        EXIT
                        !    END IF
                        !END DO
                        IF(AIK .LT. 0.D0) THEN
                            MULT(P)= -ALPHAI*AIK /AII
                        ELSE
                            MULT(P)= -BETAI *AIK /AII
                        END IF
                    END DO
                    REMAIND_ROW= PI_PTR
                    MINI= -1
                    DO WHILE(REMAIND_ROW .NE. 0)
                        MIN_COL=1000000000
                        DO P= 0, PI_PTR-1 
                            IF(MEMBER(P)) THEN
                                MIN_COL= MIN(MIN_COL, IFC_INDICES(PNTR(P)))
                                IF(MIN_COL .EQ. IFC_INDICES(PNTR(P))) THEN
                                    PINDEX_MIN= P
                                END IF
                            END IF
                        END DO
                        IF((MINI .EQ. MIN_COL)) THEN
                            NEW_VALUES(IFCI_PTR-1) = NEW_VALUES(IFCI_PTR-1) + MULT(PINDEX_MIN) * IFC_VALUES(PNTR(PINDEX_MIN)) 
                        ELSE
                            NEW_INDICES(IFCI_PTR)= MIN_COL
                            NEW_VALUES(IFCI_PTR) = MULT(PINDEX_MIN) * IFC_VALUES(PNTR(PINDEX_MIN)) 
                            IFCI_PTR= IFCI_PTR+1
                            MINI= MIN_COL
                        END IF
                        IF(PNTR(PINDEX_MIN)+1 .LE. ROW_END(PINDEX_MIN)) THEN
                            PNTR(PINDEX_MIN) = PNTR(PINDEX_MIN) +1
                        ELSE
                            REMAIND_ROW = REMAIND_ROW -1
                            MEMBER(PINDEX_MIN) = .FALSE.
                        END IF
                    END DO
                END IF
            !else
                !if(counter .eq. 2) then
                !    write(*,*) "not_interpolated", I
                !    write(*,*) "I no of strong connections= ",  S_STARTS(I+1)-1 - S_STARTS(I)
                !    DO SINDEX= S_STARTS(I), S_STARTS(I+1)-1
                !    write(*,*) "strong connected to = " ,S_INDICES(SINDEX)
                !    write(*,*) "is_fstar of its connection = " ,IS_FSTAR(S_INDICES(SINDEX))
                !    END DO
                !end if
            END IF
        END DO 
        DO II = IPREV+1, NEW_START_SIZE-1
            NEW_STARTS(II)= IFCI_PTR
        END DO
        NEW_STARTS(II)= IFCI_PTR
        
    END SUBROUTINE
   
   SUBROUTINE MERGE_IFCS_(ROWSIZE, LENGTH, S, IND, VALS, LENGTH1, S1, IND1, VALS1, LENGTH2, S2, IND2, VALS2, LENGTH3, S3, IND3, VALS3, LENGTH4, S4, IND4, VALS4, CSIZE, C, IS_C)
        INTEGER, INTENT(IN) :: CSIZE, ROWSIZE, LENGTH, LENGTH1, LENGTH2, LENGTH3, LENGTH4
        INTEGER, DIMENSION (0: ROWSIZE), INTENT(IN) :: S1, S2, S3, S4
        INTEGER, DIMENSION (0: ROWSIZE), INTENT(OUT) :: S
        INTEGER, DIMENSION (0: LENGTH-1), INTENT(OUT) :: IND
        DOUBLE PRECISION, DIMENSION (0: LENGTH-1), INTENT(OUT) :: VALS
        INTEGER, DIMENSION (0: LENGTH1-1), INTENT(IN) :: IND1
        DOUBLE PRECISION, DIMENSION (0: LENGTH1-1), INTENT(IN) :: VALS1
        INTEGER, DIMENSION (0: LENGTH2-1), INTENT(IN) :: IND2
        DOUBLE PRECISION, DIMENSION (0: LENGTH2-1), INTENT(IN) :: VALS2
        INTEGER, DIMENSION (0: LENGTH3-1), INTENT(IN) :: IND3
        DOUBLE PRECISION, DIMENSION (0: LENGTH3-1), INTENT(IN) :: VALS3
        INTEGER, DIMENSION (0: LENGTH4-1), INTENT(IN) :: IND4
        DOUBLE PRECISION, DIMENSION (0: LENGTH4-1), INTENT(IN) :: VALS4
        INTEGER, DIMENSION (0: CSIZE-1) :: C
        LOGICAL, DIMENSION (0: ROWSIZE-1) :: IS_C
        INTEGER :: ROW, JINDEX, CPTR, IFC_PTR
        CPTR = 0
        IFC_PTR= 0
        DO ROW = 0, ROWSIZE-1
            S(ROW) = IFC_PTR
            IF(IS_C(ROW)) THEN
                IND(IFC_PTR) = CPTR
                VALS(IFC_PTR)= 1.D0
                IFC_PTR= IFC_PTR +1
                CPTR= CPTR+1
            ELSE
                DO JINDEX= S1(ROW), S1(ROW+1)-1
                    IND(IFC_PTR) = IND1(JINDEX)
                    VALS(IFC_PTR)= VALS1(JINDEX)
                    IFC_PTR= IFC_PTR +1
                END DO
                DO JINDEX= S2(ROW), S2(ROW+1)-1
                    IND(IFC_PTR) = IND2(JINDEX)
                    VALS(IFC_PTR)= VALS2(JINDEX)
                    IFC_PTR= IFC_PTR +1
                END DO
                DO JINDEX= S3(ROW), S3(ROW+1)-1
                    IND(IFC_PTR)  = IND3(JINDEX)
                    VALS(IFC_PTR) = VALS3(JINDEX)
                    IFC_PTR= IFC_PTR +1
                END DO
                DO JINDEX= S4(ROW), S4(ROW+1)-1
                    IND(IFC_PTR)  = IND4(JINDEX)
                    VALS(IFC_PTR) = VALS4(JINDEX)
                    IFC_PTR= IFC_PTR +1
                END DO
            END IF
        END DO
        S(ROW) = IFC_PTR
   END SUBROUTINE
   

    SUBROUTINE AMG_GENERATE_MATRIX_(LEVEL_NO)
        INTEGER, INTENT(IN) :: LEVEL_NO
        INTEGER :: MAT_NO, SIZE, C_SIZE, LENGTH, REAL_LENGTH, IFCT_A_LENGTH, RES_LENGTH
        DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: AT
        TYPE(CRS_TYPE) :: IFCT_A
        MAT_NO = LEVEL_NO
        SIZE = MATS(MAT_NO)%ROWS                
        C_SIZE = LEVELS(LEVEL_NO)%C_SIZE
        LENGTH = MATS(MAT_NO)%LENGTH                
        RES_LENGTH = 3 * LENGTH
        IFCT_A_LENGTH = 3 * LENGTH
        ALLOCATE(IFCT_A%STARTS(0:C_SIZE))
        ALLOCATE(IFCT_A%INDICES(0:IFCT_A_LENGTH - 1))
        ALLOCATE(IFCT_A%VALUES(0:IFCT_A_LENGTH - 1))                                    
        ! GENERATING IFCT_A = IFCT * A
        ALLOCATE(AT(0:LENGTH - 1))
        CALL CRS_TRANSPOSE_SYMMETRIC_STRUCTURE_(SIZE, LENGTH, MATS(MAT_NO)%STARTS, MATS(MAT_NO)%INDICES, MATS(MAT_NO)%VALUES, AT)
        CALL CRS_MULT_MM_(C_SIZE, SIZE, SIZE, IFCT_A_LENGTH, LEVELS(LEVEL_NO)%IFCT%LENGTH, LENGTH, LEVELS(LEVEL_NO)%IFCT%STARTS, LEVELS(LEVEL_NO)%IFCT%INDICES, LEVELS(LEVEL_NO)%IFCT%VALUES, MATS(MAT_NO)%STARTS, MATS(MAT_NO)%INDICES, MATS(MAT_NO)%VALUES, MATS(MAT_NO)%STARTS, MATS(MAT_NO)%INDICES, AT, IFCT_A%STARTS, IFCT_A%INDICES, IFCT_A%VALUES)

        ALLOCATE(MATS(MAT_NO + 1)%STARTS(0:C_SIZE))
        ALLOCATE(MATS(MAT_NO + 1)%INDICES(0:RES_LENGTH - 1))
        ALLOCATE(MATS(MAT_NO + 1)%VALUES(0:RES_LENGTH - 1))                                
        ALLOCATE(MATS(MAT_NO + 1)%DIAG_INDICES(0:C_SIZE - 1))
        ALLOCATE(MATS(MAT_NO + 1)%VEC(0:C_SIZE- 1))
        ALLOCATE(MATS(MAT_NO + 1)%U(0:C_SIZE - 1))
        CALL CRS_MULT_MM_(C_SIZE, SIZE, C_SIZE, RES_LENGTH, IFCT_A_LENGTH, LEVELS(LEVEL_NO)%IFC%LENGTH, IFCT_A%STARTS, IFCT_A%INDICES, IFCT_A%VALUES, LEVELS(LEVEL_NO)%IFC%STARTS, LEVELS(LEVEL_NO)%IFC%INDICES, LEVELS(LEVEL_NO)%IFC%VALUES, LEVELS(LEVEL_NO)%IFCT%STARTS, LEVELS(LEVEL_NO)%IFCT%INDICES, LEVELS(LEVEL_NO)%IFCT%VALUES, MATS(MAT_NO + 1)%STARTS, MATS(MAT_NO + 1)%INDICES, MATS(MAT_NO + 1)%VALUES)
        REAL_LENGTH = MATS(MAT_NO + 1)%STARTS(C_SIZE)
        CALL CRS_FILL_DIAG_INDICES_(C_SIZE, REAL_LENGTH, MATS(MAT_NO + 1)%STARTS, MATS(MAT_NO + 1)%INDICES, MATS(MAT_NO + 1)%VALUES, MATS(MAT_NO + 1)%DIAG_INDICES)
        MATS(MAT_NO + 1)%ROWS = C_SIZE
        MATS(MAT_NO + 1)%LENGTH = REAL_LENGTH
        MATS(MAT_NO + 1)%MAX_CONNECTION = CRS_GET_MAX_ROW_SIZE_(C_SIZE, MATS(MAT_NO + 1)%STARTS)
        DEALLOCATE(IFCT_A%STARTS)
        DEALLOCATE(IFCT_A%INDICES)
        DEALLOCATE(IFCT_A%VALUES)

        DEALLOCATE(AT)
    END SUBROUTINE
    
    
    !******************************************   JACOBI F-SMOOTHING PART   ******************************************!

    SUBROUTINE AMG_PREPARE_JACOBI_F_SMOOTHING_(SIZE, C_SIZE, F_SIZE, A_LENGTH, INC_LENGTH, A_STARTS, A_DIAG, A_INDICES, A, INC_STARTS, INC_INDICES, INC, TFF_STARTS, TFF_INDICES, TFF, IFC_STARTS, IFC_INDICES, IFC, AFC_STARTS, AFC_INDICES, AFC, F, IS_F, CF_INDEX)
        INTEGER, INTENT(IN) :: SIZE, C_SIZE, F_SIZE, A_LENGTH, INC_LENGTH
        INTEGER, DIMENSION(0:SIZE), INTENT(IN) :: A_STARTS
        INTEGER, DIMENSION(0:SIZE - 1), INTENT(IN) :: A_DIAG
        INTEGER, DIMENSION(0:A_LENGTH - 1), INTENT(IN) :: A_INDICES
        DOUBLE PRECISION, DIMENSION(0:A_LENGTH - 1), INTENT(IN) :: A
        INTEGER, DIMENSION(0:SIZE), INTENT(IN) :: INC_STARTS
        INTEGER, DIMENSION(0:INC_LENGTH - 1), INTENT(IN) :: INC_INDICES
        DOUBLE PRECISION, DIMENSION(0:INC_LENGTH - 1), INTENT(IN) :: INC
        INTEGER, DIMENSION(0:F_SIZE), INTENT(OUT) :: TFF_STARTS
        INTEGER, DIMENSION(0:A_LENGTH - 1), INTENT(OUT) :: TFF_INDICES
        DOUBLE PRECISION, DIMENSION(0:A_LENGTH - 1), INTENT(OUT) :: TFF
        INTEGER, DIMENSION(0:F_SIZE), INTENT(OUT) :: AFC_STARTS
        INTEGER, DIMENSION(0:A_LENGTH - 1), INTENT(OUT) :: AFC_INDICES
        DOUBLE PRECISION, DIMENSION(0:A_LENGTH - 1), INTENT(OUT) :: AFC
        INTEGER, DIMENSION(0:F_SIZE), INTENT(OUT) :: IFC_STARTS
        INTEGER, DIMENSION(0:INC_LENGTH - 1), INTENT(OUT) :: IFC_INDICES
        DOUBLE PRECISION, DIMENSION(0:INC_LENGTH - 1), INTENT(OUT) :: IFC
        INTEGER, DIMENSION(0:F_SIZE - 1), INTENT(IN) :: F
        LOGICAL, DIMENSION(0:SIZE - 1), INTENT(IN) :: IS_F
        INTEGER, DIMENSION(0:SIZE - 1), INTENT(IN) :: CF_INDEX
        INTEGER :: I, ID, TFF_PTR, IFC_PTR, AFC_PTR, DIAG_INDEX, J_INDEX, J
        DOUBLE PRECISION :: S
        TFF_PTR = 0
        IFC_PTR = 0
        AFC_PTR = 0        
        DO I = 0, F_SIZE - 1            
            ID = F(I)
            DIAG_INDEX = A_DIAG(ID)
            S = A(DIAG_INDEX)
            TFF_STARTS(I) = TFF_PTR
            IFC_STARTS(I) = IFC_PTR
            AFC_STARTS(I) = AFC_PTR
            DO J_INDEX = A_STARTS(ID), A_STARTS(ID + 1) - 1
                J = A_INDICES(J_INDEX)
                IF (IS_F(J)) THEN
                    IF (ID .NE. J) THEN
                        TFF_INDICES(TFF_PTR) = CF_INDEX(J)
                        TFF(TFF_PTR) = -A(J_INDEX) / S
                        TFF_PTR = TFF_PTR + 1
                    END IF
                ELSE
                    AFC_INDICES(AFC_PTR) = CF_INDEX(J)
                    AFC(AFC_PTR) = -A(J_INDEX) / S
                    AFC_PTR = AFC_PTR + 1
                END IF
            END DO
            DO J_INDEX = INC_STARTS(ID), INC_STARTS(ID + 1) - 1
                J = INC_INDICES(J_INDEX)
                IFC_INDICES(IFC_PTR) = J
                IFC(IFC_PTR) = INC(J_INDEX)
                IFC_PTR = IFC_PTR + 1
            END DO
        END DO
        TFF_STARTS(I) = TFF_PTR
        IFC_STARTS(I) = IFC_PTR
        AFC_STARTS(I) = AFC_PTR
    END SUBROUTINE
    
    SUBROUTINE AMG_FINALIZE_JACOBI_F_SMOOTHING_(SIZE, F_SIZE, IFC_LENGTH, INC_LENGTH, IFC_STARTS, IFC_INDICES, IFC, INC_STARTS, INC_INDICES, INC, F, IS_F, CF_INDEX)
        INTEGER, INTENT(IN) :: SIZE, F_SIZE, IFC_LENGTH, INC_LENGTH
        INTEGER, DIMENSION(0:F_SIZE), INTENT(IN) :: IFC_STARTS
        INTEGER, DIMENSION(0:IFC_LENGTH - 1), INTENT(IN) :: IFC_INDICES
        DOUBLE PRECISION, DIMENSION(0:IFC_LENGTH - 1), INTENT(IN) :: IFC
        INTEGER, DIMENSION(0:SIZE), INTENT(OUT) :: INC_STARTS
        INTEGER, DIMENSION(0:INC_LENGTH - 1), INTENT(OUT) :: INC_INDICES
        DOUBLE PRECISION, DIMENSION(0:INC_LENGTH - 1), INTENT(OUT) :: INC
        INTEGER, DIMENSION(0:F_SIZE - 1), INTENT(IN) :: F
        LOGICAL, DIMENSION(0:SIZE - 1), INTENT(IN) :: IS_F
        INTEGER, DIMENSION(0:SIZE - 1), INTENT(IN) :: CF_INDEX
        INTEGER :: I, J_INDEX, J, PTR, F_INDEX
        PTR = 0
        DO I = 0, SIZE - 1
            INC_STARTS(I) = PTR
            IF (IS_F(I)) THEN
                F_INDEX = CF_INDEX(I)
                DO J_INDEX = IFC_STARTS(F_INDEX), IFC_STARTS(F_INDEX + 1) - 1
                    J = IFC_INDICES(J_INDEX)
                    INC_INDICES(PTR) = J
                    INC(PTR) = IFC(J_INDEX)
                    PTR = PTR + 1
                END DO
            ELSE
                INC_INDICES(PTR) = CF_INDEX(I)
                INC(PTR) = 1.D0
                PTR = PTR + 1
            END IF
        END DO
        INC_STARTS(I) = PTR
    END SUBROUTINE
    
    SUBROUTINE AMG_JACOBI_F_SMOOTHING_(LEVEL_NO)
        INTEGER, INTENT(IN) :: LEVEL_NO
        TYPE(CRS_TYPE) :: IFC, TFF, AFC, TFF_IFC, IFCT, RES
        ! S = DIAG(AFF)
        ! TFF = INV(S) * (S - AFF)
        ! AFC = INV(S) * (-AFC)
        INTEGER :: MAT_NO, SIZE, C_SIZE, F_SIZE, INC_LENGTH, A_LENGTH, TFF_IFC_LENGTH
        MAT_NO = LEVEL_NO
        C_SIZE = LEVELS(LEVEL_NO)%C_SIZE
        F_SIZE = LEVELS(LEVEL_NO)%F_SIZE
        SIZE = C_SIZE + F_SIZE
        A_LENGTH = MATS(MAT_NO)%LENGTH
        INC_LENGTH = LEVELS(LEVEL_NO)%IFC%LENGTH
        TFF_IFC_LENGTH = A_LENGTH
        ALLOCATE(IFC%STARTS(0:F_SIZE))
        ALLOCATE(IFC%INDICES(0:INC_LENGTH - 1))
        ALLOCATE(IFC%VALUES(0:INC_LENGTH - 1))

        ALLOCATE(IFCT%STARTS(0:C_SIZE))
        ALLOCATE(IFCT%INDICES(0:INC_LENGTH - 1))
        ALLOCATE(IFCT%VALUES(0:INC_LENGTH - 1))
        
        ALLOCATE(TFF%STARTS(0:F_SIZE))
        ALLOCATE(TFF%INDICES(0:A_LENGTH - 1))
        ALLOCATE(TFF%VALUES(0:A_LENGTH - 1))
        
        ALLOCATE(AFC%STARTS(0:F_SIZE))
        ALLOCATE(AFC%INDICES(0:A_LENGTH - 1))
        ALLOCATE(AFC%VALUES(0:A_LENGTH - 1))
        
        
        CALL AMG_PREPARE_JACOBI_F_SMOOTHING_(SIZE, C_SIZE, F_SIZE, A_LENGTH, INC_LENGTH, MATS(MAT_NO)%STARTS, MATS(MAT_NO)%DIAG_INDICES, MATS(MAT_NO)%INDICES, MATS(MAT_NO)%VALUES, LEVELS(LEVEL_NO)%IFC%STARTS, LEVELS(LEVEL_NO)%IFC%INDICES, LEVELS(LEVEL_NO)%IFC%VALUES, TFF%STARTS, TFF%INDICES, TFF%VALUES, IFC%STARTS, IFC%INDICES, IFC%VALUES, AFC%STARTS, AFC%INDICES, AFC%VALUES, LEVELS(LEVEL_NO)%F, LEVELS(LEVEL_NO)%IS_F, LEVELS(LEVEL_NO)%CF_INDEX)
        
        ALLOCATE(TFF_IFC%STARTS(0:F_SIZE))
        ALLOCATE(TFF_IFC%INDICES(0:TFF_IFC_LENGTH - 1))
        ALLOCATE(TFF_IFC%VALUES(0:TFF_IFC_LENGTH - 1))
        
        CALL CRS_TRANSPOSE_(F_SIZE, C_SIZE, INC_LENGTH, IFC%STARTS, IFC%INDICES, IFC%VALUES, IFCT%STARTS, IFCT%INDICES, IFCT%VALUES)
        CALL CRS_MULT_MM_(F_SIZE, F_SIZE, C_SIZE, TFF_IFC_LENGTH, A_LENGTH, INC_LENGTH, TFF%STARTS, TFF%INDICES, TFF%VALUES, IFC%STARTS, IFC%INDICES, IFC%VALUES, IFCT%STARTS, IFCT%INDICES, IFCT%VALUES, TFF_IFC%STARTS, TFF_IFC%INDICES, TFF_IFC%VALUES)
        ALLOCATE(RES%STARTS(0:F_SIZE))
        ALLOCATE(RES%INDICES(0:A_LENGTH - 1))
        ALLOCATE(RES%VALUES(0:A_LENGTH - 1))
        CALL CRS_ADD_MM_(F_SIZE, TFF_IFC_LENGTH, A_LENGTH, TFF_IFC_LENGTH, TFF_IFC%STARTS, TFF_IFC%INDICES, TFF_IFC%VALUES, AFC%STARTS, AFC%INDICES, AFC%VALUES, RES%STARTS, RES%INDICES, RES%VALUES)
        DEALLOCATE(LEVELS(LEVEL_NO)%IFC%INDICES)
        DEALLOCATE(LEVELS(LEVEL_NO)%IFC%VALUES)
        DEALLOCATE(LEVELS(LEVEL_NO)%IFCT%INDICES)
        DEALLOCATE(LEVELS(LEVEL_NO)%IFCT%VALUES)
        
        
        ALLOCATE(LEVELS(LEVEL_NO)%IFC%INDICES(0:A_LENGTH - 1))
        ALLOCATE(LEVELS(LEVEL_NO)%IFC%VALUES(0:A_LENGTH - 1))
        ALLOCATE(LEVELS(LEVEL_NO)%IFCT%INDICES(0:A_LENGTH - 1))
        ALLOCATE(LEVELS(LEVEL_NO)%IFCT%VALUES(0:A_LENGTH - 1))
        CALL AMG_FINALIZE_JACOBI_F_SMOOTHING_(SIZE, F_SIZE, A_LENGTH, A_LENGTH, RES%STARTS, RES%INDICES, RES%VALUES, LEVELS(LEVEL_NO)%IFC%STARTS, LEVELS(LEVEL_NO)%IFC%INDICES, LEVELS(LEVEL_NO)%IFC%VALUES, LEVELS(LEVEL_NO)%F, LEVELS(LEVEL_NO)%IS_F, LEVELS(LEVEL_NO)%CF_INDEX)
        INC_LENGTH = LEVELS(LEVEL_NO)%IFC%STARTS(SIZE)
        CALL CRS_TRANSPOSE_(SIZE, C_SIZE, INC_LENGTH, LEVELS(LEVEL_NO)%IFC%STARTS, LEVELS(LEVEL_NO)%IFC%INDICES, LEVELS(LEVEL_NO)%IFC%VALUES, LEVELS(LEVEL_NO)%IFCT%STARTS, LEVELS(LEVEL_NO)%IFCT%INDICES, LEVELS(LEVEL_NO)%IFCT%VALUES)
        LEVELS(LEVEL_NO)%IFC%LENGTH = INC_LENGTH
        LEVELS(LEVEL_NO)%IFCT%LENGTH = INC_LENGTH
        
        DEALLOCATE(IFC%STARTS)
        DEALLOCATE(IFC%INDICES)
        DEALLOCATE(IFC%VALUES)
        
        DEALLOCATE(IFCT%STARTS)
        DEALLOCATE(IFCT%INDICES)
        DEALLOCATE(IFCT%VALUES)
        
        
        DEALLOCATE(TFF%STARTS)
        DEALLOCATE(TFF%INDICES)
        DEALLOCATE(TFF%VALUES)
        
        DEALLOCATE(AFC%STARTS)
        DEALLOCATE(AFC%INDICES)
        DEALLOCATE(AFC%VALUES)
        
        DEALLOCATE(TFF_IFC%STARTS)
        DEALLOCATE(TFF_IFC%INDICES)
        DEALLOCATE(TFF_IFC%VALUES)
        DEALLOCATE(RES%STARTS)
        DEALLOCATE(RES%INDICES)
        DEALLOCATE(RES%VALUES)
        
    END SUBROUTINE


    !******************************************   RELAXATION PART   ******************************************!
    
    SUBROUTINE AMG_CREATE_RELAXER_(MAT_NO, JACOBI, GS, ILU, PCG, PBICG)
        INTEGER, INTENT(IN) :: MAT_NO
        LOGICAL, INTENT(IN) :: GS, ILU, PCG, PBICG, JACOBI
        INTEGER :: SIZE
        SIZE = MATS(MAT_NO)%ROWS
        ALLOCATE(MATS(MAT_NO)%Y(0:SIZE - 1))        
        IF (JACOBI) THEN
            ALLOCATE(MATS(MAT_NO)%J_PARAMS%COPY(0:SIZE - 1))
        END IF
        IF (GS) THEN
            ALLOCATE(MATS(MAT_NO)%GS_PARAMS%COPY(0:SIZE - 1))
        END IF
        IF (ILU) THEN
            ALLOCATE(MATS(MAT_NO)%ILU_PARAMS%COPY(0:SIZE - 1))
            ALLOCATE(MATS(MAT_NO)%ILU_PARAMS%RHS(0:SIZE - 1))
        END IF
        IF (PCG) THEN
            ALLOCATE(MATS(MAT_NO)%PCG_PARAMS%R(0:SIZE - 1))
            ALLOCATE(MATS(MAT_NO)%PCG_PARAMS%Z(0:SIZE - 1))
            ALLOCATE(MATS(MAT_NO)%PCG_PARAMS%P(0:SIZE - 1))
            ALLOCATE(MATS(MAT_NO)%PCG_PARAMS%Q(0:SIZE - 1))
            ALLOCATE(MATS(MAT_NO)%PCG_PARAMS%COPY(0:SIZE - 1))
        END IF
        IF (PBICG) THEN
            ALLOCATE(MATS(MAT_NO)%PBICG_PARAMS%R(0:SIZE - 1))
            ALLOCATE(MATS(MAT_NO)%PBICG_PARAMS%P(0:SIZE - 1))
            ALLOCATE(MATS(MAT_NO)%PBICG_PARAMS%P_HAT(0:SIZE - 1))
            ALLOCATE(MATS(MAT_NO)%PBICG_PARAMS%S(0:SIZE - 1))
            ALLOCATE(MATS(MAT_NO)%PBICG_PARAMS%S_HAT(0:SIZE - 1))
            ALLOCATE(MATS(MAT_NO)%PBICG_PARAMS%T(0:SIZE - 1))
            ALLOCATE(MATS(MAT_NO)%PBICG_PARAMS%V(0:SIZE - 1))
            ALLOCATE(MATS(MAT_NO)%PBICG_PARAMS%RR(0:SIZE - 1))
            ALLOCATE(MATS(MAT_NO)%PBICG_PARAMS%RMIN(0:SIZE - 1))
            ALLOCATE(MATS(MAT_NO)%PBICG_PARAMS%COPY(0:SIZE - 1))
            ALLOCATE(MATS(MAT_NO)%PBICG_PARAMS%XMIN(0:SIZE - 1))
        END IF
    END SUBROUTINE
    
    SUBROUTINE AMG_CREATE_ILU0_(MAT_NO, JACOBI, GS, ILU, PCG, PBICG)        
        INTEGER, INTENT(IN) :: MAT_NO
        LOGICAL, INTENT(IN) :: GS, ILU, PCG, PBICG, JACOBI
        INTEGER :: SIZE, LENGTH
        CALL AMG_CREATE_RELAXER_(MAT_NO, JACOBI, GS, ILU, PCG, PBICG)
        SIZE = MATS(MAT_NO)%ROWS
        LENGTH = MATS(MAT_NO)%LENGTH
        MATS(MAT_NO)%ILU%ROWS = SIZE
        MATS(MAT_NO)%ILU%LENGTH = LENGTH
        ALLOCATE(MATS(MAT_NO)%ILU%STARTS(0:SIZE))
        ALLOCATE(MATS(MAT_NO)%ILU%DIAG_INDICES(0:SIZE - 1))
        ALLOCATE(MATS(MAT_NO)%ILU%INDICES(0:LENGTH - 1))
        ALLOCATE(MATS(MAT_NO)%ILU%VALUES(0:LENGTH - 1))
        CALL ILU0_INITIALIZE_(SIZE, LENGTH, MATS(MAT_NO)%STARTS, MATS(MAT_NO)%INDICES, MATS(MAT_NO)%DIAG_INDICES, MATS(MAT_NO)%VALUES, MATS(MAT_NO)%ILU%STARTS, MATS(MAT_NO)%ILU%INDICES, MATS(MAT_NO)%ILU%DIAG_INDICES, MATS(MAT_NO)%ILU%VALUES)        
        CALL ILU0_CREATE_PRECONDITIONER_(SIZE, LENGTH, MATS(MAT_NO)%ILU%STARTS, MATS(MAT_NO)%ILU%INDICES, MATS(MAT_NO)%ILU%DIAG_INDICES, MATS(MAT_NO)%ILU%VALUES)        
    END SUBROUTINE
    
    
    !******************************************   V-CYCLE PART   ******************************************!

    SUBROUTINE AMG_ZERO_INITIAL_GUESS_(SIZE, U)
        INTEGER, INTENT(IN) :: SIZE
        DOUBLE PRECISION, DIMENSION(0:SIZE - 1), INTENT(OUT) :: U
        INTEGER :: I
        !DEC$ SIMD
        DO I = 0, SIZE - 1
            U(I) = 0.D0
        END DO
    END SUBROUTINE
    
    SUBROUTINE AMG_DIRECT_SOLVE_(METHOD)
        INTEGER, INTENT(IN) :: METHOD
        INTEGER :: LEVEL_NO, MAT_NO, ITER
        MAT_NO = MAX_LEVEL
        IF (METHOD .EQ. LU_) THEN
            ! TODO
        ELSEIF (METHOD .EQ. PCG_) THEN
            ITER = PCG_SOLVE_(MATS(MAT_NO)%ROWS, MATS(MAT_NO)%LENGTH, MATS(MAT_NO)%ILU%LENGTH, MATS(MAT_NO)%STARTS, MATS(MAT_NO)%INDICES, MATS(MAT_NO)%ILU%STARTS, MATS(MAT_NO)%ILU%INDICES, MATS(MAT_NO)%ILU%DIAG_INDICES, MATS(MAT_NO)%VALUES, MATS(MAT_NO)%ILU%VALUES, MATS(MAT_NO)%VEC, MATS(MAT_NO)%Y, MATS(MAT_NO)%U, 200, MATS(MAT_NO)%PCG_PARAMS)        
        ELSEIF (METHOD .EQ. PBICG_) THEN
            CALL ILU0_TRI_SOLVE_(MATS(MAT_NO)%ROWS, MATS(MAT_NO)%ILU%LENGTH, MATS(MAT_NO)%ILU%STARTS, MATS(MAT_NO)%ILU%INDICES, MATS(MAT_NO)%ILU%DIAG_INDICES, MATS(MAT_NO)%ILU%VALUES, MATS(MAT_NO)%VEC, MATS(MAT_NO)%U, MATS(MAT_NO)%Y)
            ITER = PBICG_SOLVE_(MATS(MAT_NO)%ROWS, MATS(MAT_NO)%LENGTH, MATS(MAT_NO)%ILU%LENGTH, MATS(MAT_NO)%STARTS, MATS(MAT_NO)%INDICES, MATS(MAT_NO)%ILU%STARTS, MATS(MAT_NO)%ILU%INDICES, MATS(MAT_NO)%ILU%DIAG_INDICES, MATS(MAT_NO)%VALUES, MATS(MAT_NO)%ILU%VALUES, MATS(MAT_NO)%VEC, MATS(MAT_NO)%Y, MATS(MAT_NO)%U, 10000, MATS(MAT_NO)%PBICG_PARAMS)        
            
            ! COMMENT
            !WRITE (*, 2) ITER, CRS_VECTOR_NORM_(MATS(MAT_NO)%ROWS, MATS(MAT_NO)%U)

            ! TODO
        ELSEIF (METHOD .EQ. GS_) THEN
            ! BUG
            !ITER = GS_SOLVE_(MATS(MAT_NO)%ROWS, MATS(MAT_NO)%LENGTH, MATS(MAT_NO)%STARTS, MATS(MAT_NO)%INDICES, MATS(MAT_NO)%DIAG_INDICES, MATS(MAT_NO)%VALUES, MATS(MAT_NO)%VEC, MATS(MAT_NO)%U, 10, MATS(MAT_NO)%GS_PARAMS, .TRUE.)
        END IF
        !CALL CRS_PRINT_(MATS(MAT_NO)%ROWS, MATS(MAT_NO)%LENGTH, MATS(MAT_NO)%STARTS, MATS(MAT_NO)%INDICES, MATS(MAT_NO)%VALUES)
        !CALL CRS_PRINT_VECS_(MATS(MAT_NO)%ROWS, MATS(MAT_NO)%VEC, MATS(MAT_NO)%U)
2       FORMAT ('PBICG CONVERGED AFTER ', I6, ' ITERATIONS. WITE RESULT NORM = ', E22.15)
    END SUBROUTINE
    
    RECURSIVE SUBROUTINE AMG_V_CYCLE_(LEVEL_NO, N1, N2, LEVEL_STARTS, THREADS, INDEPENDENT_ROW_STARTS, INDEPENDENT_ROWS)
        ! NEW BUG INTENT(IN)
        INTEGER, INTENT(IN) :: LEVEL_NO, N1, N2   
        INTEGER, INTENT(IN) :: LEVEL_STARTS, THREADS
        INTEGER, DIMENSION(0:LEVEL_STARTS - 1), INTENT(IN) :: INDEPENDENT_ROW_STARTS
        INTEGER, DIMENSION(0:MATS(LEVEL_NO)%ROWS - MATS(LEVEL_NO)%W_SIZE - 1), INTENT(IN) :: INDEPENDENT_ROWS
        
        INTEGER :: MAT_NO, SIZE, C_SIZE, F_SIZE, ITER
        LOGICAL :: NAN_FLAG
        REAL :: T11, T22
        INTEGER :: I
        MAT_NO = LEVEL_NO
        SIZE = MATS(MAT_NO)%ROWS
        IF (LEVEL_NO .LT. MAX_LEVEL) THEN
            ! COMMENT
            !WRITE (*,*) '======================  @ LEVEL = ', LEVEL_NO
            C_SIZE = LEVELS(LEVEL_NO)%C_SIZE
            F_SIZE = LEVELS(LEVEL_NO)%F_SIZE
            
            ! COMMENT
            !call CRS_RESIDUAL_(MATS(MAT_NO)%ROWS, MATS(MAT_NO)%ROWS, MATS(MAT_NO)%LENGTH, MATS(MAT_NO)%STARTS, MATS(MAT_NO)%INDICES, MATS(MAT_NO)%VALUES, MATS(MAT_NO)%VEC, MATS(MAT_NO)%U, MATS(MAT_NO)%Y)    
            !WRITE (*,*) 'BEFORE INNER INITIAL RESIDUAL = ', CRS_VECTOR_NORM_(MATS(MAT_NO)%ROWS, MATS(MAT_NO)%Y)
            
            !ITER = JACOBI_SOLVE_(MATS(MAT_NO)%ROWS, MATS(MAT_NO)%LENGTH, MATS(MAT_NO)%STARTS, MATS(MAT_NO)%INDICES, MATS(MAT_NO)%DIAG_INDICES, MATS(MAT_NO)%VALUES, MATS(MAT_NO)%VEC, MATS(MAT_NO)%U, N1, MATS(MAT_NO)%J_PARAMS%COPY, 2.D0 / 3.D0)            
            !IF (LEVEL_NO .EQ. 0) THEN
            !    ITER = GS_SOLVE_PARALLEL_FORWARD_(MATS(MAT_NO)%ROWS, MATS(MAT_NO)%W_SIZE, LEVEL_STARTS, MATS(MAT_NO)%LENGTH, THREADS, MATS(MAT_NO)%STARTS, MATS(MAT_NO)%INDICES, MATS(MAT_NO)%DIAG_INDICES, INDEPENDENT_ROW_STARTS, INDEPENDENT_ROWS, MATS(MAT_NO)%VALUES, MATS(MAT_NO)%VEC, MATS(MAT_NO)%Y, MATS(MAT_NO)%U, N1, MATS(MAT_NO)%GS_PARAMS)                        
            !ELSE
                ITER = GS_SOLVE_(MATS(MAT_NO)%ROWS, MATS(MAT_NO)%LENGTH, MATS(MAT_NO)%STARTS, MATS(MAT_NO)%INDICES, MATS(MAT_NO)%DIAG_INDICES, MATS(MAT_NO)%VALUES, MATS(MAT_NO)%VEC, MATS(MAT_NO)%U, N1, MATS(MAT_NO)%GS_PARAMS, .TRUE., LEVELS(LEVEL_NO)%FORCED_C)            
            !WRITE (*,*) 'GS_SOLVE_RESULT_NORM (', LEVEL_NO, ') = ', CRS_VECTOR_NORM_(MATS(MAT_NO)%ROWS, MATS(MAT_NO)%U)
            !END IF
            
            !ITER = GS_CF_SOLVE_(MATS(MAT_NO)%ROWS, MATS(MAT_NO)%LENGTH, LEVELS(LEVEL_NO)%C_SIZE, LEVELS(LEVEL_NO)%F_SIZE, MATS(MAT_NO)%STARTS, MATS(MAT_NO)%INDICES, MATS(MAT_NO)%DIAG_INDICES, MATS(MAT_NO)%VALUES, MATS(MAT_NO)%VEC, MATS(MAT_NO)%U, N1, MATS(MAT_NO)%GS_PARAMS, LEVELS(LEVEL_NO)%C, LEVELS(LEVEL_NO)%F, .TRUE.)            
            !ITER = ILU0_SOLVE_(MATS(MAT_NO)%ROWS, MATS(MAT_NO)%LENGTH, MATS(MAT_NO)%ILU%LENGTH, MATS(MAT_NO)%STARTS, MATS(MAT_NO)%INDICES, MATS(MAT_NO)%ILU%STARTS, MATS(MAT_NO)%ILU%INDICES, MATS(MAT_NO)%ILU%DIAG_INDICES, MATS(MAT_NO)%VALUES, MATS(MAT_NO)%ILU%VALUES, MATS(MAT_NO)%VEC, MATS(MAT_NO)%Y, MATS(MAT_NO)%U, N1, MATS(MAT_NO)%ILU_PARAMS)
            !CALL PCG_SOLVE_(MATS(MAT_NO)%ROWS, MATS(MAT_NO)%LENGTH, MATS(MAT_NO)%ILU%LENGTH, MATS(MAT_NO)%STARTS, MATS(MAT_NO)%INDICES, MATS(MAT_NO)%ILU%STARTS, MATS(MAT_NO)%ILU%INDICES, MATS(MAT_NO)%ILU%DIAG_INDICES, MATS(MAT_NO)%VALUES, MATS(MAT_NO)%ILU%VALUES, MATS(MAT_NO)%VEC, MATS(MAT_NO)%Y, MATS(MAT_NO)%U, N1, MATS(MAT_NO)%PCG_PARAMS)        
            !ITER = PBICG_SOLVE_(MATS(MAT_NO)%ROWS, MATS(MAT_NO)%LENGTH, MATS(MAT_NO)%ILU%LENGTH, MATS(MAT_NO)%STARTS, MATS(MAT_NO)%INDICES, MATS(MAT_NO)%ILU%STARTS, MATS(MAT_NO)%ILU%INDICES, MATS(MAT_NO)%ILU%DIAG_INDICES, MATS(MAT_NO)%VALUES, MATS(MAT_NO)%ILU%VALUES, MATS(MAT_NO)%VEC, MATS(MAT_NO)%Y, MATS(MAT_NO)%U, N1, MATS(MAT_NO)%PBICG_PARAMS)        
                
            ! COMMENT                
            !call CRS_RESIDUAL_(MATS(MAT_NO)%ROWS, MATS(MAT_NO)%ROWS, MATS(MAT_NO)%LENGTH, MATS(MAT_NO)%STARTS, MATS(MAT_NO)%INDICES, MATS(MAT_NO)%VALUES, MATS(MAT_NO)%VEC, MATS(MAT_NO)%U, MATS(MAT_NO)%Y)    
            !WRITE (*,*) 'BEFORE INNER FINAL RESIDUAL = ', CRS_VECTOR_NORM_(MATS(MAT_NO)%ROWS, MATS(MAT_NO)%Y)
            
            ! GENERATING B_C                                                            
            CALL CRS_RESIDUAL_(SIZE, SIZE, MATS(MAT_NO)%LENGTH, MATS(MAT_NO)%STARTS, MATS(MAT_NO)%INDICES, MATS(MAT_NO)%VALUES, MATS(MAT_NO)%VEC, MATS(MAT_NO)%U, MATS(MAT_NO)%Y)
            CALL CRS_MULT_MV_(C_SIZE, SIZE, LEVELS(LEVEL_NO)%IFCT%LENGTH, LEVELS(LEVEL_NO)%IFCT%STARTS, LEVELS(LEVEL_NO)%IFCT%INDICES, LEVELS(LEVEL_NO)%IFCT%VALUES, MATS(MAT_NO)%Y, MATS(MAT_NO + 1)%VEC)
            
            CALL AMG_ZERO_INITIAL_GUESS_(C_SIZE, MATS(MAT_NO + 1)%U)
            CALL AMG_V_CYCLE_(LEVEL_NO + 1, N1, N2, LEVEL_STARTS, THREADS, INDEPENDENT_ROW_STARTS, INDEPENDENT_ROWS)
            ! GENERATING U
            CALL CRS_MULT_UPDATE_MV_(SIZE, C_SIZE, LEVELS(LEVEL_NO)%IFC%LENGTH, LEVELS(LEVEL_NO)%IFC%STARTS, LEVELS(LEVEL_NO)%IFC%INDICES, LEVELS(LEVEL_NO)%IFC%VALUES, MATS(MAT_NO + 1)%U, MATS(MAT_NO)%U, 1.D0)                       
            
            ! COMMENT
            !call CRS_RESIDUAL_(MATS(MAT_NO)%ROWS, MATS(MAT_NO)%ROWS, MATS(MAT_NO)%LENGTH, MATS(MAT_NO)%STARTS, MATS(MAT_NO)%INDICES, MATS(MAT_NO)%VALUES, MATS(MAT_NO)%VEC, MATS(MAT_NO)%U, MATS(MAT_NO)%Y)    
            !WRITE (*,*) 'AFTER INNER INITIAL RESIDUAL = ', CRS_VECTOR_NORM_(MATS(MAT_NO)%ROWS, MATS(MAT_NO)%Y)            
            
            !ITER = JACOBI_SOLVE_(MATS(MAT_NO)%ROWS, MATS(MAT_NO)%LENGTH, MATS(MAT_NO)%STARTS, MATS(MAT_NO)%INDICES, MATS(MAT_NO)%DIAG_INDICES, MATS(MAT_NO)%VALUES, MATS(MAT_NO)%VEC, MATS(MAT_NO)%U, N1, MATS(MAT_NO)%J_PARAMS%COPY, 2.D0 / 3.D0)            
            !IF (LEVEL_NO .EQ. 0) THEN
            !    ITER = GS_SOLVE_PARALLEL_BACKWARD_(MATS(MAT_NO)%ROWS, MATS(MAT_NO)%W_SIZE, LEVEL_STARTS, MATS(MAT_NO)%LENGTH, THREADS, MATS(MAT_NO)%STARTS, MATS(MAT_NO)%INDICES, MATS(MAT_NO)%DIAG_INDICES, INDEPENDENT_ROW_STARTS, INDEPENDENT_ROWS, MATS(MAT_NO)%VALUES, MATS(MAT_NO)%VEC, MATS(MAT_NO)%Y, MATS(MAT_NO)%U, N1, MATS(MAT_NO)%GS_PARAMS)                        
            !ELSE
                ITER = GS_SOLVE_(MATS(MAT_NO)%ROWS, MATS(MAT_NO)%LENGTH, MATS(MAT_NO)%STARTS, MATS(MAT_NO)%INDICES, MATS(MAT_NO)%DIAG_INDICES, MATS(MAT_NO)%VALUES, MATS(MAT_NO)%VEC, MATS(MAT_NO)%U, N1, MATS(MAT_NO)%GS_PARAMS, .FALSE., LEVELS(LEVEL_NO)%FORCED_C)
            !WRITE (*,*) 'GS_SOLVE_RESULT_NORM (', LEVEL_NO, ') = ', CRS_VECTOR_NORM_(MATS(MAT_NO)%ROWS, MATS(MAT_NO)%U)
            !END IF
            
            !ITER = GS_CF_SOLVE_(MATS(MAT_NO)%ROWS, MATS(MAT_NO)%LENGTH, LEVELS(LEVEL_NO)%C_SIZE, LEVELS(LEVEL_NO)%F_SIZE, MATS(MAT_NO)%STARTS, MATS(MAT_NO)%INDICES, MATS(MAT_NO)%DIAG_INDICES, MATS(MAT_NO)%VALUES, MATS(MAT_NO)%VEC, MATS(MAT_NO)%U, N1, MATS(MAT_NO)%GS_PARAMS, LEVELS(LEVEL_NO)%C, LEVELS(LEVEL_NO)%F, .FALSE.)
            !ITER = ILU0_SOLVE_(MATS(MAT_NO)%ROWS, MATS(MAT_NO)%LENGTH, MATS(MAT_NO)%ILU%LENGTH, MATS(MAT_NO)%STARTS, MATS(MAT_NO)%INDICES, MATS(MAT_NO)%ILU%STARTS, MATS(MAT_NO)%ILU%INDICES, MATS(MAT_NO)%ILU%DIAG_INDICES, MATS(MAT_NO)%VALUES, MATS(MAT_NO)%ILU%VALUES, MATS(MAT_NO)%VEC, MATS(MAT_NO)%Y, MATS(MAT_NO)%U, N1, MATS(MAT_NO)%ILU_PARAMS)
            !ITER = PCG_SOLVE_(MATS(MAT_NO)%ROWS, MATS(MAT_NO)%LENGTH, MATS(MAT_NO)%ILU%LENGTH, MATS(MAT_NO)%STARTS, MATS(MAT_NO)%INDICES, MATS(MAT_NO)%ILU%STARTS, MATS(MAT_NO)%ILU%INDICES, MATS(MAT_NO)%ILU%DIAG_INDICES, MATS(MAT_NO)%VALUES, MATS(MAT_NO)%ILU%VALUES, MATS(MAT_NO)%VEC, MATS(MAT_NO)%Y, MATS(MAT_NO)%U, N2, MATS(MAT_NO)%PCG_PARAMS)
            !ITER = PBICG_SOLVE_(MATS(MAT_NO)%ROWS, MATS(MAT_NO)%LENGTH, MATS(MAT_NO)%ILU%LENGTH, MATS(MAT_NO)%STARTS, MATS(MAT_NO)%INDICES, MATS(MAT_NO)%ILU%STARTS, MATS(MAT_NO)%ILU%INDICES, MATS(MAT_NO)%ILU%DIAG_INDICES, MATS(MAT_NO)%VALUES, MATS(MAT_NO)%ILU%VALUES, MATS(MAT_NO)%VEC, MATS(MAT_NO)%Y, MATS(MAT_NO)%U, N2, MATS(MAT_NO)%PBICG_PARAMS)        
                
            ! COMMENT                
            !call CRS_RESIDUAL_(MATS(MAT_NO)%ROWS, MATS(MAT_NO)%ROWS, MATS(MAT_NO)%LENGTH, MATS(MAT_NO)%STARTS, MATS(MAT_NO)%INDICES, MATS(MAT_NO)%VALUES, MATS(MAT_NO)%VEC, MATS(MAT_NO)%U, MATS(MAT_NO)%Y)    
            !WRITE (*,*) 'AFTER INNER FINAL RESIDUAL = ', CRS_VECTOR_NORM_(MATS(MAT_NO)%ROWS, MATS(MAT_NO)%Y)
        ELSE
            CALL AMG_DIRECT_SOLVE_(PBICG_)
        END IF
    END SUBROUTINE
    
    
END MODULE