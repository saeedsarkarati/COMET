MODULE CRS

    USE OMP_LIB
    USE HEAP_MOD
    IMPLICIT NONE
    
    TYPE :: CRS_TYPE
        INTEGER :: ROWS, LENGTH
        INTEGER, DIMENSION(:), ALLOCATABLE :: STARTS, INDICES, DIAG_INDICES
        DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: VALUES
    END TYPE
        
    CONTAINS
    
    FUNCTION GET_MATRIX_NORM_(LENGTH, VALUES) RESULT(ANS)
        INTEGER, INTENT(IN) :: LENGTH
        DOUBLE PRECISION, DIMENSION(0:LENGTH - 1), INTENT(IN) :: VALUES
        INTEGER :: INDEX
        DOUBLE PRECISION :: ANS
        DOUBLE PRECISION :: VALUE
        ANS = 0.D0
        DO INDEX = 0, LENGTH - 1
            VALUE = VALUES(INDEX)
            ANS = ANS + VALUE * VALUE
        END DO
        ANS = SQRT(ANS / LENGTH)
    END FUNCTION
    
    FUNCTION HAS_SYMMETRIC_STRUCTURE_(SIZE, LENGTH, STARTS, INDICES) RESULT(ANS)
        INTEGER, INTENT(IN) :: SIZE, LENGTH
        INTEGER, DIMENSION(0:SIZE), INTENT(IN) :: STARTS
        INTEGER, DIMENSION(0:LENGTH - 1), INTENT(IN) :: INDICES
        LOGICAL :: ANS
        INTEGER :: I, J_INDEX, J
        ANS = .TRUE.
        DO I = 0, SIZE - 1
            DO J_INDEX = STARTS(I), STARTS(I + 1) - 1
                J = INDICES(J_INDEX)
                IF (CRS_GET_INDEX_(J, I, SIZE, LENGTH, STARTS, INDICES) .EQ. -1) THEN
                    WRITE (*,*) I, J
                    ANS = .FALSE.
                END IF
            END DO
        END DO
    END FUNCTION
    
    FUNCTION CRS_GET_INDEX_(I, J, SIZE, LENGTH, STARTS, INDICES) RESULT(INDEX)
        INTEGER, INTENT(IN) :: I, J, SIZE, LENGTH
        INTEGER, DIMENSION(0:SIZE), INTENT(IN) :: STARTS
        INTEGER, DIMENSION(0:LENGTH - 1), INTENT(IN) :: INDICES        
        INTEGER :: INDEX
        INTEGER :: START_INDEX, END_INDEX, ID
        START_INDEX = STARTS(I)
        END_INDEX = STARTS(I + 1) - 1
        DO WHILE (START_INDEX .LE. END_INDEX)
            INDEX = (START_INDEX + END_INDEX) / 2
            ID = INDICES(INDEX)
            IF (ID .EQ. J) THEN
                RETURN
            ELSE IF (ID .GT. J) THEN
                END_INDEX = INDEX - 1
            ELSE
                START_INDEX = INDEX + 1
            END IF                            
        END DO                	
        INDEX = -1
    END FUNCTION
    
    FUNCTION CRS_GET_(I, J, INDEX, SIZE, LENGTH, STARTS, INDICES, VALUES) RESULT(ANS)
        INTEGER, INTENT(IN) :: I, J, SIZE, LENGTH
        INTEGER, INTENT(OUT) :: INDEX
        INTEGER, DIMENSION(0:SIZE), INTENT(IN) :: STARTS
        INTEGER, DIMENSION(0:LENGTH - 1), INTENT(IN) :: INDICES
        DOUBLE PRECISION, DIMENSION(0:LENGTH - 1), INTENT(IN) :: VALUES
        DOUBLE PRECISION :: ANS
        INTEGER :: II, JJ, C0, C1, C2
        INTEGER :: START_INDEX, END_INDEX, ID
        START_INDEX = STARTS(I)
        END_INDEX = STARTS(I + 1) - 1
        DO WHILE (START_INDEX .LE. END_INDEX)
            INDEX = (START_INDEX + END_INDEX) / 2
            ID = INDICES(INDEX)
            IF (ID .EQ. J) THEN
                ANS = VALUES(INDEX)
                RETURN
            ELSE IF (ID .GT. J) THEN
                END_INDEX = INDEX - 1
            ELSE
                START_INDEX = INDEX + 1
            END IF                            
        END DO                	
        INDEX = -1
        ANS = 0.D0
    END FUNCTION
    
    FUNCTION CRS_GET_NEXT_INDEX_(I, J, SIZE, LENGTH, STARTS, INDICES) RESULT(ANS)
        INTEGER, INTENT(IN) :: I, J, SIZE, LENGTH
        INTEGER, DIMENSION(0:SIZE), INTENT(IN) :: STARTS
        INTEGER, DIMENSION(0:LENGTH - 1), INTENT(IN) :: INDICES
        INTEGER :: ANS
        INTEGER :: START_INDEX, END_INDEX, ID, INDEX
        START_INDEX = STARTS(I)
        END_INDEX = STARTS(I + 1) - 1
        DO WHILE (START_INDEX .LE. END_INDEX)
            INDEX = (START_INDEX + END_INDEX) / 2
            ID = INDICES(INDEX)
            IF (ID .EQ. J) THEN
                END_INDEX = INDEX
                EXIT
            ELSE IF (ID .GT. J) THEN
                END_INDEX = INDEX - 1
            ELSE
                START_INDEX = INDEX + 1
            END IF                            
        END DO                	
        ANS = END_INDEX + 1
    END FUNCTION
    
    SUBROUTINE CRS_PRINT_(SIZE, LENGTH, STARTS, INDICES, VALUES)
        INTEGER, INTENT(IN) :: SIZE, LENGTH 
        INTEGER, DIMENSION(0:SIZE), INTENT(IN) :: STARTS
        INTEGER, DIMENSION(0:LENGTH - 1), INTENT(IN) :: INDICES
        DOUBLE PRECISION, DIMENSION(0:LENGTH - 1), INTENT(IN) :: VALUES
        INTEGER :: FILE_UNIT = 300000, IERROR, I, J_INDEX, J
        OPEN(UNIT = FILE_UNIT, FILE = 'FILES/CRS.dat', STATUS = 'REPLACE', ACTION = 'WRITE', ACCESS = 'SEQUENTIAL', IOSTAT = IERROR)            
        DO I = 0, SIZE - 1
            DO J_INDEX = STARTS(I), STARTS(I + 1) - 1
                J = INDICES(J_INDEX)
                WRITE (FILE_UNIT, 30) I + 1, J + 1, VALUES(J_INDEX)
            END DO
        END DO
30      FORMAT (I10, I10, E27.15)
        CLOSE(FILE_UNIT)
    END SUBROUTINE
    
    SUBROUTINE CRS_PRINT_VECS_(SIZE, B, U, CPR_LP)
        INTEGER, INTENT(IN) :: SIZE
        DOUBLE PRECISION, DIMENSION(0:SIZE - 1), INTENT(IN) :: B, U
        DOUBLE PRECISION, DIMENSION(0:3 * SIZE - 1), INTENT(IN) :: CPR_LP
        INTEGER :: FILE_UNIT = 300001, IERROR, I
        OPEN(UNIT = FILE_UNIT, FILE = 'FILES/CRS_VECS.dat', STATUS = 'REPLACE', ACTION = 'WRITE', ACCESS = 'SEQUENTIAL', IOSTAT = IERROR)            
        DO I = 0, SIZE - 1
            WRITE (FILE_UNIT, 40) B(I), U(I), CPR_LP(3 * I + 2)
        END DO
        CLOSE(FILE_UNIT)
40      FORMAT (E27.15, E27.15, E27.15)
    END SUBROUTINE

    SUBROUTINE CRS_PRINT_ARRAY_(SIZE, ARR)
        INTEGER, INTENT(IN) :: SIZE
        INTEGER, DIMENSION(0:SIZE - 1), INTENT(IN) :: ARR
        INTEGER :: FILE_UNIT = 300002, IERROR, I
        OPEN(UNIT = FILE_UNIT, FILE = 'FILES/ARRAY.dat', STATUS = 'REPLACE', ACTION = 'WRITE', ACCESS = 'SEQUENTIAL', IOSTAT = IERROR)            
        DO I = 0, SIZE - 1
            WRITE (FILE_UNIT, 50) ARR(I)
        END DO
        CLOSE(FILE_UNIT)
50      FORMAT (I10)
    END SUBROUTINE
    
    SUBROUTINE CRS_PRINT_VEC_(SIZE, V)
        INTEGER, INTENT(IN) :: SIZE
        DOUBLE PRECISION, DIMENSION(0:SIZE - 1), INTENT(IN) :: V
        INTEGER :: FILE_UNIT = 300002, IERROR, I
        OPEN(UNIT = FILE_UNIT, FILE = 'FILES/CRS_VEC.dat', STATUS = 'REPLACE', ACTION = 'WRITE', ACCESS = 'SEQUENTIAL', IOSTAT = IERROR)            
        DO I = 0, SIZE - 1
            WRITE (FILE_UNIT, 50) V(I)
        END DO
        CLOSE(FILE_UNIT)
50      FORMAT (E27.15)
    END SUBROUTINE
    
    SUBROUTINE CRS_PRINT_INT_VEC_(SIZE, v)
        INTEGER, INTENT(IN) :: SIZE
        INTEGER, DIMENSION(0:SIZE - 1), INTENT(IN) :: v
        INTEGER :: FILE_UNIT = 300002, IERROR, I
        OPEN(UNIT = FILE_UNIT, FILE = 'FILES/VEC.dat', STATUS = 'REPLACE', ACTION = 'WRITE', ACCESS = 'SEQUENTIAL', IOSTAT = IERROR)            
        DO I = 0, SIZE - 1
            WRITE (FILE_UNIT, 50) V(I)
        END DO
        CLOSE(FILE_UNIT)
50      FORMAT (I10)
    END SUBROUTINE
    
    SUBROUTINE CRS_MULT_MV_(C_SIZE, V_SIZE, LENGTH, STARTS, INDICES, A, V, C)
        ! C = SCALER * C + A*V
        INTEGER, INTENT(IN) :: C_SIZE, V_SIZE, LENGTH
        INTEGER, DIMENSION(0:C_SIZE) :: STARTS
        INTEGER, DIMENSION(0:LENGTH - 1), INTENT(IN) :: INDICES
        DOUBLE PRECISION, DIMENSION(0:LENGTH - 1), INTENT(IN) :: A
        DOUBLE PRECISION, DIMENSION(0:V_SIZE - 1), INTENT(IN) :: V
        DOUBLE PRECISION, DIMENSION(0:C_SIZE - 1), INTENT(OUT) :: C
        INTEGER :: I, J_INDEX, J
        DOUBLE PRECISION :: SUM
        !$OMP PARALLEL SHARED(STARTS, INDICES, A, V, C)
            !$OMP DO SCHEDULE(STATIC) PRIVATE(I, J_INDEX, J, SUM)            
        DO I = 0, C_SIZE - 1
            SUM = 0.D0            
            DO J_INDEX = STARTS(I), STARTS(I + 1) - 1
                J = INDICES(J_INDEX)
                SUM = SUM + A(J_INDEX) * V(J)
            END DO
            C(I) = SUM
        END DO
            !$OMP END DO NOWAIT
        !$OMP END PARALLEL        
    END SUBROUTINE

    
    SUBROUTINE CRS_MULT_UPDATE_MV_(C_SIZE, V_SIZE, LENGTH, STARTS, INDICES, A, V, C, SCALER)
        ! C = SCALER * C + A*V
        INTEGER, INTENT(IN) :: C_SIZE, V_SIZE, LENGTH
        INTEGER, DIMENSION(0:C_SIZE) :: STARTS
        INTEGER, DIMENSION(0:LENGTH - 1), INTENT(IN) :: INDICES
        DOUBLE PRECISION, DIMENSION(0:LENGTH - 1), INTENT(IN) :: A
        DOUBLE PRECISION, DIMENSION(0:V_SIZE - 1), INTENT(IN) :: V
        DOUBLE PRECISION, DIMENSION(0:C_SIZE - 1), INTENT(OUT) :: C
        DOUBLE PRECISION, INTENT(IN) :: SCALER
        INTEGER :: I, J_INDEX, J
        DOUBLE PRECISION :: SUM
        !$OMP PARALLEL SHARED(STARTS, INDICES, A, V, C)
            !$OMP DO SCHEDULE(STATIC) PRIVATE(I, J_INDEX, J, SUM)            
        DO I = 0, C_SIZE - 1
            SUM = 0.D0            
            DO J_INDEX = STARTS(I), STARTS(I + 1) - 1
                J = INDICES(J_INDEX)
                SUM = SUM + A(J_INDEX) * V(J)
            END DO
            C(I) = SCALER * C(I) + SUM
        END DO
            !$OMP END DO NOWAIT
        !$OMP END PARALLEL        
    END SUBROUTINE
    
    SUBROUTINE CRS_RESIDUAL_(C_SIZE, V_SIZE, LENGTH, STARTS, INDICES, A, B, V, C)    
        ! C = B - A*V
        INTEGER, INTENT(IN) :: C_SIZE, V_SIZE, LENGTH
        INTEGER, DIMENSION(0:C_SIZE) :: STARTS
        INTEGER, DIMENSION(0:LENGTH - 1), INTENT(IN) :: INDICES
        DOUBLE PRECISION, DIMENSION(0:LENGTH - 1), INTENT(IN) :: A
        DOUBLE PRECISION, DIMENSION(0:C_SIZE - 1), INTENT(IN) :: B
        DOUBLE PRECISION, DIMENSION(0:V_SIZE - 1), INTENT(IN) :: V
        DOUBLE PRECISION, DIMENSION(0:C_SIZE - 1), INTENT(OUT) :: C
        INTEGER :: I, J_INDEX, J
        DOUBLE PRECISION :: SUM
        !$OMP PARALLEL SHARED(STARTS, INDICES, A, B, V, C)
            !$OMP DO SCHEDULE(STATIC) PRIVATE(I, J_INDEX, J, SUM)            
        DO I = 0, C_SIZE - 1
            SUM = B(I)            
            DO J_INDEX = STARTS(I), STARTS(I + 1) - 1
                J = INDICES(J_INDEX)
                SUM = SUM - A(J_INDEX) * V(J)
            END DO
            C(I) = SUM
        END DO
            !$OMP END DO NOWAIT 
        !$OMP END PARALLEL        
    END SUBROUTINE

    SUBROUTINE CRS_MULT_COUNT_PTR_(SIZE, LENGTH, STARTS, INDICES, BASE, PTR, PTR_VALUE, GOAL, FLAG)
        INTEGER, INTENT(IN) :: LENGTH, SIZE, BASE
        INTEGER, DIMENSION(0:SIZE), INTENT(IN) :: STARTS
        INTEGER, DIMENSION(0:LENGTH - 1), INTENT(IN) :: INDICES
        INTEGER, INTENT(INOUT) :: PTR, PTR_VALUE
        INTEGER, INTENT(IN) :: GOAL
        LOGICAL, INTENT(OUT) :: FLAG
        DO WHILE (PTR_VALUE .LT. GOAL)
            PTR = PTR + 1
            IF (PTR .EQ. STARTS(BASE + 1)) THEN
                FLAG = .FALSE.
                EXIT
            END IF
            PTR_VALUE = INDICES(PTR)
        END DO        
    END SUBROUTINE
        
    
    SUBROUTINE CRS_MULT_MM_(C_SIZE, B_SIZE, BT_SIZE, C_LENGTH, A_LENGTH, B_LENGTH, A_STARTS, A_INDICES, A, B_STARTS, B_INDICES, B, BT_STARTS, BT_INDICES, BT, C_STARTS, C_INDICES, C)
        ! C = A * B
        INTEGER, INTENT(IN) :: C_SIZE, B_SIZE, BT_SIZE, C_LENGTH, A_LENGTH, B_LENGTH
        INTEGER, DIMENSION(0:C_SIZE), INTENT(IN) :: A_STARTS
        INTEGER, DIMENSION(0:B_SIZE), INTENT(IN) :: B_STARTS
        INTEGER, DIMENSION(0:BT_SIZE), INTENT(IN) :: BT_STARTS
        INTEGER, DIMENSION(0:C_SIZE), INTENT(OUT) :: C_STARTS
        INTEGER, DIMENSION(0:A_LENGTH - 1), INTENT(IN) :: A_INDICES
        INTEGER, DIMENSION(0:B_LENGTH - 1), INTENT(IN) ::  B_INDICES, BT_INDICES
        INTEGER, DIMENSION(0:C_LENGTH - 1), INTENT(OUT) :: C_INDICES
        DOUBLE PRECISION, DIMENSION(0:A_LENGTH - 1), INTENT(IN) :: A
        DOUBLE PRECISION, DIMENSION(0:B_LENGTH - 1), INTENT(IN) :: B, BT
        DOUBLE PRECISION, DIMENSION(0:C_LENGTH - 1), INTENT(OUT) :: C        
        INTEGER :: I, J_INDEX, J, K_INDEX, K
        INTEGER :: A_PTR, BT_PTR, A_COL, BT_ROW, CUR, COLUMN, ACCESSOR, COUNT
        LOGICAL :: A_FLAG, BT_FLAG
        DOUBLE PRECISION :: SUM
        LOGICAL, DIMENSION(:), ALLOCATABLE :: FLAG
        ALLOCATE(FLAG(0:BT_SIZE - 1))
        !DIR$ SIMD
        DO I = 0, BT_SIZE - 1
            FLAG(I) = .FALSE.
        END DO
        ACCESSOR = 0
        DO I = 0, C_SIZE - 1
            C_STARTS(I) = ACCESSOR
            CALL HEAP_CLEAR_()
            DO J_INDEX = A_STARTS(I), A_STARTS(I + 1) - 1
                J = A_INDICES(J_INDEX)
                DO K_INDEX = B_STARTS(J), B_STARTS(J + 1) - 1
                    K = B_INDICES(K_INDEX)
                    IF (.NOT. FLAG(K)) THEN
                        CALL HEAP_INSERT_(K)
                        FLAG(K) = .TRUE.
                    END IF
                END DO
            END DO            
            !CUR = -1
            !COUNT = 0
            DO WHILE (HEAP_SIZE .GT. 0) 
                COLUMN = HEAP_EXTRACT_MIN_()
                C_INDICES(ACCESSOR) = COLUMN
                ACCESSOR = ACCESSOR + 1
                FLAG(COLUMN) = .FALSE.
                !IF (COLUMN .NE. CUR) THEN
                !    C_INDICES(ACCESSOR) = COLUMN
                !    ACCESSOR = ACCESSOR + 1
                !    CUR = COLUMN                    
                !ELSE
                !    COUNT = COUNT + 1
                !END IF
            END DO
        END DO
        C_STARTS(I) = ACCESSOR
        DEALLOCATE(FLAG)
        !$OMP PARALLEL
            !$OMP DO SCHEDULE(STATIC) PRIVATE(I, J_INDEX, J, A_PTR, BT_PTR, SUM, A_FLAG, BT_FLAG, A_COL, BT_ROW)            
        DO I = 0, C_SIZE - 1
            DO J_INDEX = C_STARTS(I), C_STARTS(I + 1) - 1
                J = C_INDICES(J_INDEX)
                A_PTR = A_STARTS(I)
                BT_PTR = BT_STARTS(J)
                SUM = 0.D0
                A_FLAG = (A_STARTS(I + 1) - A_STARTS(I)) .GT. 0
                BT_FLAG = (BT_STARTS(J + 1) - BT_STARTS(J)) .GT. 0
                DO WHILE (A_FLAG .AND. BT_FLAG)
                    A_COL = A_INDICES(A_PTR)
                    BT_ROW = BT_INDICES(BT_PTR)                    
                    IF (A_COL .LT. BT_ROW) THEN
                        CALL CRS_MULT_COUNT_PTR_(C_SIZE, A_LENGTH, A_STARTS, A_INDICES, I, A_PTR, A_COL, BT_ROW, A_FLAG)    
                    ELSEIF (A_COL .GT. BT_ROW) THEN
                        CALL CRS_MULT_COUNT_PTR_(BT_SIZE, B_LENGTH, BT_STARTS, BT_INDICES, J, BT_PTR, BT_ROW, A_COL, BT_FLAG)
                    ELSE
                        SUM = SUM + A(A_PTR) * BT(BT_PTR)
                        A_PTR = A_PTR + 1
                        BT_PTR = BT_PTR + 1
                        A_FLAG = A_PTR .LT. A_STARTS(I + 1)
                        BT_FLAG = BT_PTR .LT. BT_STARTS(J + 1)
                    END IF
                END DO
                C(J_INDEX) = SUM
            END DO
        END DO
            !$OMP END DO NOWAIT
        !$OMP END PARALLEL
    END SUBROUTINE

    FUNCTION CRS_GET_MAX_ROW_SIZE_(SIZE, STARTS) RESULT(ANS)
        INTEGER, INTENT(IN) :: SIZE
        INTEGER, DIMENSION(0:SIZE), INTENT(IN) :: STARTS
        INTEGER :: ANS, ROW_SIZE
        INTEGER :: I
        ANS = 0
        DO I = 0, SIZE - 1
            ROW_SIZE = STARTS(I + 1) - STARTS(I)
            IF (ROW_SIZE .GT. ANS) THEN
                ANS = ROW_SIZE
            END IF
        END DO
    END FUNCTION
    
    FUNCTION CRS_GET_MAT_MAX_(SIZE, LENGTH, STARTS, VALUES) RESULT(ANS)
        INTEGER, INTENT(IN) :: SIZE, LENGTH
        INTEGER, DIMENSION(0:SIZE), INTENT(IN) :: STARTS
        DOUBLE PRECISION, DIMENSION(0:LENGTH - 1), INTENT(IN) :: VALUES
        DOUBLE PRECISION :: ANS, VALUE
        INTEGER :: I, J_INDEX
        ANS = 0.D0
        DO I = 0, SIZE - 1
            DO J_INDEX = STARTS(I), STARTS(I + 1) - 1
                VALUE = ABS(VALUES(J_INDEX))
                IF (VALUE .GT. ANS) THEN
                    ANS = VALUE
                END IF
            END DO
        END DO
    END FUNCTION
    
    FUNCTION CRS_GET_VEC_MAX_(SIZE, VALUES) RESULT(ANS)
        INTEGER, INTENT(IN) :: SIZE
        DOUBLE PRECISION, DIMENSION(0:SIZE - 1), INTENT(IN) :: VALUES
        DOUBLE PRECISION :: ANS, VALUE
        INTEGER :: I
        ANS = 0.D0
        DO I = 0, SIZE - 1
            VALUE = ABS(VALUES(I))
            IF (VALUE .GT. ANS) THEN
                ANS = VALUE
            END IF
        END DO
    END FUNCTION
    
    SUBROUTINE CRS_ADD_COUNT_PTR_(SIZE, LENGTH, RESULT_LENGTH, STARTS, INDICES, VALUES, BASE, PTR, PTR_VALUE, GOAL, FLAG, RESULT_PTR, RESULT_INDICES, RESULT_VALUES)
        INTEGER, INTENT(IN) :: LENGTH, SIZE, BASE, RESULT_LENGTH
        INTEGER, DIMENSION(0:SIZE), INTENT(IN) :: STARTS
        INTEGER, DIMENSION(0:LENGTH - 1), INTENT(IN) :: INDICES
        DOUBLE PRECISION, DIMENSION(0:LENGTH - 1), INTENT(IN) :: VALUES
        INTEGER, DIMENSION(0:RESULT_LENGTH - 1), INTENT(INOUT) :: RESULT_INDICES
        DOUBLE PRECISION, DIMENSION(0:RESULT_LENGTH - 1), INTENT(INOUT) :: RESULT_VALUES
        INTEGER, INTENT(INOUT) :: PTR, PTR_VALUE, RESULT_PTR
        INTEGER, INTENT(IN) :: GOAL
        LOGICAL, INTENT(OUT) :: FLAG
        DO WHILE (PTR_VALUE .LT. GOAL)
            RESULT_INDICES(RESULT_PTR) = PTR_VALUE
            RESULT_VALUES(RESULT_PTR) = VALUES(PTR)
            RESULT_PTR = RESULT_PTR + 1
            PTR = PTR + 1
            IF (PTR .EQ. STARTS(BASE + 1)) THEN
                FLAG = .FALSE.
                EXIT
            END IF
            PTR_VALUE = INDICES(PTR)
        END DO        
    END SUBROUTINE

    SUBROUTINE CRS_ADD_FILL_REMAINING_(SIZE, LENGTH, RESULT_LENGTH, STARTS, INDICES, VALUES, BASE, PTR, RESULT_PTR, RESULT_INDICES, RESULT_VALUES)
        INTEGER, INTENT(IN) :: LENGTH, SIZE, BASE, RESULT_LENGTH
        INTEGER, DIMENSION(0:SIZE), INTENT(IN) :: STARTS
        INTEGER, DIMENSION(0:LENGTH - 1), INTENT(IN) :: INDICES
        DOUBLE PRECISION, DIMENSION(0:LENGTH - 1), INTENT(IN) :: VALUES
        INTEGER, DIMENSION(0:RESULT_LENGTH - 1), INTENT(INOUT) :: RESULT_INDICES
        DOUBLE PRECISION, DIMENSION(0:RESULT_LENGTH - 1), INTENT(INOUT) :: RESULT_VALUES
        INTEGER, INTENT(INOUT) :: PTR, RESULT_PTR
        INTEGER :: J_INDEX
        DO WHILE (PTR .LT. STARTS(BASE + 1))
            RESULT_INDICES(RESULT_PTR) = INDICES(PTR)
            RESULT_VALUES(RESULT_PTR) = VALUES(PTR)
            RESULT_PTR = RESULT_PTR + 1
            PTR = PTR + 1
        END DO
    END SUBROUTINE
    
    SUBROUTINE CRS_ADD_MM_(SIZE, A_LENGTH, B_LENGTH, C_LENGTH, A_STARTS, A_INDICES, A, B_STARTS, B_INDICES, B, C_STARTS, C_INDICES, C)
        ! C = A + B
        INTEGER, INTENT(IN) :: SIZE, A_LENGTH, B_LENGTH, C_LENGTH
        INTEGER, DIMENSION(0:SIZE), INTENT(IN) :: A_STARTS, B_STARTS
        INTEGER, DIMENSION(0:SIZE), INTENT(OUT) :: C_STARTS
        INTEGER, DIMENSION(0:A_LENGTH - 1), INTENT(IN) :: A_INDICES
        INTEGER, DIMENSION(0:B_LENGTH - 1), INTENT(IN) :: B_INDICES
        INTEGER, DIMENSION(0:C_LENGTH - 1), INTENT(OUT) :: C_INDICES        
        DOUBLE PRECISION, DIMENSION(0:A_LENGTH - 1), INTENT(IN) :: A
        DOUBLE PRECISION, DIMENSION(0:B_LENGTH - 1), INTENT(IN) :: B
        DOUBLE PRECISION, DIMENSION(0:C_LENGTH - 1), INTENT(OUT) :: C
        INTEGER :: I, C_PTR, A_PTR, B_PTR, A_COL, B_COL
        LOGICAL :: A_FLAG, B_FLAG
        C_PTR = 0
        DO I = 0, SIZE - 1
            C_STARTS(I) = C_PTR
            A_PTR = A_STARTS(I)
            B_PTR = B_STARTS(I)
            A_FLAG = (A_STARTS(I + 1) - A_STARTS(I)) .GT. 0
            B_FLAG = (B_STARTS(I + 1) - B_STARTS(I)) .GT. 0
            DO WHILE (A_FLAG .AND. B_FLAG)
                A_COL = A_INDICES(A_PTR)
                B_COL = B_INDICES(B_PTR)
                IF (A_COL .LT. B_COL) THEN
                    CALL CRS_ADD_COUNT_PTR_(SIZE, A_LENGTH, C_LENGTH, A_STARTS, A_INDICES, A, I, A_PTR, A_COL, B_COL, A_FLAG, C_PTR, C_INDICES, C)
                ELSEIF (A_COL .GT. B_COL) THEN
                    CALL CRS_ADD_COUNT_PTR_(SIZE, B_LENGTH, C_LENGTH, B_STARTS, B_INDICES, B, I, B_PTR, B_COL, A_COL, B_FLAG, C_PTR, C_INDICES, C)
                ELSE
                    C_INDICES(C_PTR) = A_COL
                    C(C_PTR) = A(A_PTR) + B(B_PTR)
                    A_PTR = A_PTR + 1
                    B_PTR = B_PTR + 1
                    A_FLAG = A_PTR .LT. A_STARTS(I + 1)
                    B_FLAG = B_PTR .LT. B_STARTS(I + 1)
                    C_PTR = C_PTR + 1
                END IF
            END DO
            CALL CRS_ADD_FILL_REMAINING_(SIZE, A_LENGTH, C_LENGTH, A_STARTS, A_INDICES, A, I, A_PTR, C_PTR, C_INDICES, C)
            CALL CRS_ADD_FILL_REMAINING_(SIZE, B_LENGTH, C_LENGTH, B_STARTS, B_INDICES, B, I, B_PTR, C_PTR, C_INDICES, C)
        END DO
        C_STARTS(I) = C_PTR
    END SUBROUTINE
    
    SUBROUTINE CRS_TRANSPOSE_SYMMETRIC_STRUCTURE_(SIZE, LENGTH, STARTS, INDICES, A, AT)
        INTEGER, INTENT(IN) :: SIZE, LENGTH
        INTEGER, DIMENSION(0:SIZE), INTENT(IN) :: STARTS
        INTEGER, DIMENSION(0:LENGTH - 1), INTENT(IN) :: INDICES
        DOUBLE PRECISION, DIMENSION(0:LENGTH - 1), INTENT(IN) :: A
        DOUBLE PRECISION, DIMENSION(0:LENGTH - 1), INTENT(OUT) :: AT
        INTEGER, DIMENSION(:), ALLOCATABLE :: TEMP_STARTS
        INTEGER :: I, J_INDEX, J
        ALLOCATE(TEMP_STARTS(0:SIZE - 1))
        !DEC$ SIMD
        DO I = 0, SIZE - 1
            TEMP_STARTS(I) = STARTS(I)
        END DO
        DO I = 0, SIZE - 1
            DO J_INDEX = STARTS(I), STARTS(I + 1) - 1
                J = INDICES(J_INDEX)
                AT(TEMP_STARTS(J)) = A(J_INDEX)
                TEMP_STARTS(J) = TEMP_STARTS(J) + 1
            END DO
        END DO
        DEALLOCATE(TEMP_STARTS)
    END SUBROUTINE

    SUBROUTINE CRS_TRANSPOSE_(A_SIZE, AT_SIZE, LENGTH, A_STARTS, A_INDICES, A, AT_STARTS, AT_INDICES, AT)
        INTEGER, INTENT(IN) :: A_SIZE, AT_SIZE, LENGTH
        INTEGER, DIMENSION(0:A_SIZE), INTENT(IN) :: A_STARTS
        INTEGER, DIMENSION(0:AT_SIZE), INTENT(OUT) :: AT_STARTS
        INTEGER, DIMENSION(0:LENGTH - 1), INTENT(IN) :: A_INDICES
        INTEGER, DIMENSION(0:LENGTH - 1), INTENT(OUT) :: AT_INDICES
        DOUBLE PRECISION, DIMENSION(0:LENGTH - 1), INTENT(IN) :: A
        DOUBLE PRECISION, DIMENSION(0:LENGTH - 1), INTENT(OUT) :: AT
        INTEGER, DIMENSION(:), ALLOCATABLE :: TEMP_STARTS
        INTEGER :: I, J_INDEX, J
        ALLOCATE(TEMP_STARTS(0:AT_SIZE - 1))
        !DEC$ SIMD
        DO I = 0, AT_SIZE - 1
            TEMP_STARTS(I) = 0
        END DO
        DO I = 0, A_SIZE - 1
            DO J_INDEX = A_STARTS(I), A_STARTS(I + 1) - 1
                J = A_INDICES(J_INDEX)
                TEMP_STARTS(J) = TEMP_STARTS(J) + 1
            END DO
        END DO   
        AT_STARTS(AT_SIZE) = LENGTH
        DO I = AT_SIZE - 1, 0, -1
            AT_STARTS(I) = AT_STARTS(I + 1) - TEMP_STARTS(I)
            TEMP_STARTS(I) = AT_STARTS(I)
        END DO
        DO I = 0, A_SIZE - 1
            DO J_INDEX = A_STARTS(I), A_STARTS(I + 1) - 1
                J = A_INDICES(J_INDEX)
                AT_INDICES(TEMP_STARTS(J)) = I
                AT(TEMP_STARTS(J)) = A(J_INDEX)
                TEMP_STARTS(J) = TEMP_STARTS(J) + 1
            END DO
        END DO
        DEALLOCATE(TEMP_STARTS)
    END SUBROUTINE
    
    FUNCTION CRS_INNER_PRODUCT_(SIZE, V1, V2) RESULT(ANS)
        INTEGER, INTENT(IN) :: SIZE
        DOUBLE PRECISION, DIMENSION(0:SIZE - 1), INTENT(IN) :: V1, V2
        DOUBLE PRECISION :: ANS
        INTEGER :: I
        ANS = 0.D0        
        !DEC$ SIMD
        DO I = 0, SIZE - 1
            ANS = ANS + V1(I) * V2(I)
        END DO
    END FUNCTION

    SUBROUTINE CRS_FILL_DIAG_INDICES_(SIZE, LENGTH, STARTS, INDICES, VALUES, DIAG_INDICES)
        INTEGER, INTENT(IN) :: SIZE, LENGTH
        INTEGER, DIMENSION(0:SIZE), INTENT(IN) :: STARTS
        INTEGER, DIMENSION(0:LENGTH - 1), INTENT(IN) :: INDICES
        DOUBLE PRECISION, DIMENSION(0:LENGTH - 1), INTENT(IN) :: VALUES
        INTEGER, DIMENSION(0:SIZE - 1), INTENT(OUT) :: DIAG_INDICES
        INTEGER :: I, J_INDEX, J
        !$OMP PARALLEL
            !$OMP DO SCHEDULE(STATIC) PRIVATE(I, J_INDEX, J)            
        DO I = 0, SIZE - 1
            DO J_INDEX = STARTS(I), STARTS(I + 1) - 1
                J = INDICES(J_INDEX)
                IF (I .EQ. J) THEN
                    DIAG_INDICES(I) = J_INDEX
                END IF
            END DO
        END DO
            !$OMP END DO NOWAIT
        !$OMP END PARALLEL
    END SUBROUTINE
    
    FUNCTION CRS_VECTOR_NORM_(SIZE, V) RESULT(ANS)
        INTEGER, INTENT(IN) :: SIZE
        DOUBLE PRECISION, DIMENSION(0:SIZE - 1), INTENT(IN) :: V
        DOUBLE PRECISION :: ANS, VALUE
        INTEGER :: I
        ANS = 0.D0        
        !DEC$ SIMD
        DO I = 0, SIZE - 1
            VALUE = V(I)
            ANS = ANS + VALUE * VALUE
        END DO
        ANS = SQRT(ANS / SIZE)        
    END FUNCTION
    
    FUNCTION CRS_VEC_ISNAN_(SIZE, V) RESULT(ANS)
        INTEGER, INTENT(IN) :: SIZE
        DOUBLE PRECISION, DIMENSION(0:SIZE - 1), INTENT(IN) :: V
        LOGICAL :: ANS
        INTEGER :: I
        ANS = .FALSE.
        DO I = 0, SIZE - 1
            IF (ISNAN(V(I))) THEN
                WRITE (*,*) I
                ANS = .TRUE.
                PAUSE
            END IF
        END DO
    END FUNCTION
    
    SUBROUTINE CRS_COPY_(SIZE, SOURCE, DEST)
        INTEGER, INTENT(IN) :: SIZE
        DOUBLE PRECISION, DIMENSION(0:SIZE - 1), INTENT(IN) :: SOURCE
        DOUBLE PRECISION, DIMENSION(0:SIZE - 1), INTENT(OUT) :: DEST
        INTEGER :: I
        !DEC$ SIMD
        DO I = 0, SIZE - 1
            DEST(I) = SOURCE(I)
        END DO
    END SUBROUTINE
    
    FUNCTION CRS_CONVERGED_(SIZE, X, COPY, REL_TOL, ABS_TOL) RESULT(ANS)
        INTEGER, INTENT(IN) :: SIZE
        DOUBLE PRECISION, DIMENSION(0:SIZE - 1), INTENT(IN) :: X, COPY
        DOUBLE PRECISION, INTENT(IN) :: REL_TOL, ABS_TOL
        LOGICAL :: ANS
        DOUBLE PRECISION :: DIFF, VALUE, TEMP_DIFF, TEMP_VALUE
        INTEGER :: I
        DIFF = 0.D0
        VALUE = 0.D0
        !DEC$ SIMD
        DO I = 0, SIZE - 1
            TEMP_DIFF = (X(I) - COPY(I))
            DIFF = DIFF + TEMP_DIFF * TEMP_DIFF
            TEMP_VALUE = X(I)
            VALUE = VALUE + TEMP_VALUE * TEMP_VALUE
        END DO
        ANS = SQRT(DIFF / SIZE) .LE. (REL_TOL * SQRT(VALUE / SIZE) + ABS_TOL)
    END FUNCTION
    
END MODULE
    
