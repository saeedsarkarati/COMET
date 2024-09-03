MODULE RED_BLACK_TREE
    
    IMPLICIT NONE
    
    INTEGER, DIMENSION(:), ALLOCATABLE :: POTENTIAL
    INTEGER, DIMENSION(:), ALLOCATABLE :: PARENT, LEFT, RIGHT
    LOGICAL, DIMENSION(:), ALLOCATABLE :: COLOR
    INTEGER :: RB_SIZE
    INTEGER :: ROOT
    
CONTAINS
    
    SUBROUTINE RB_INIT_(MAX_SIZE)
        INTEGER, INTENT(IN) :: MAX_SIZE
        INTEGER :: I
        RB_SIZE = MAX_SIZE        
        ALLOCATE(POTENTIAL(0:RB_SIZE - 1))
        ALLOCATE(PARENT(-1:RB_SIZE - 1))
        ALLOCATE(LEFT(0:RB_SIZE - 1))
        ALLOCATE(RIGHT(0:RB_SIZE - 1))
        ALLOCATE(COLOR(-1:RB_SIZE - 1))
        ROOT = -1
        !DEC$ SIMD
        DO I = 0, RB_SIZE - 1
            COLOR(I) = .FALSE.
            PARENT(I) = -1
            LEFT(I) = -1
            RIGHT(I) = -1
            POTENTIAL(I) = -1
        END DO
        COLOR(-1) = .TRUE.
        RB_SIZE = 0
    END SUBROUTINE
    
    SUBROUTINE RB_FINALIZE_()
        RB_SIZE = 0
        ROOT = -1
        DEALLOCATE(POTENTIAL)
        DEALLOCATE(PARENT)
        DEALLOCATE(LEFT)
        DEALLOCATE(RIGHT)
        DEALLOCATE(COLOR)
    END SUBROUTINE
    
    SUBROUTINE RB_LEFT_ROTATE_(ID)
        INTEGER, INTENT(IN) :: ID
        INTEGER :: Y, X
        X = ID
        Y = RIGHT(X)
        RIGHT(X) = LEFT(Y)
        IF (LEFT(Y) .NE. -1) THEN
            PARENT(LEFT(Y)) = X
        END IF
        PARENT(Y) = PARENT(X)
        IF (PARENT(X) .EQ. -1) THEN
            ROOT = Y
        ELSEIF (X .EQ. LEFT(PARENT(X))) THEN
            LEFT(PARENT(X)) = Y
        ELSE
            RIGHT(PARENT(X)) = Y
        END IF
        LEFT(Y) = X
        PARENT(X) = Y
    END SUBROUTINE
    
    SUBROUTINE RB_RIGHT_ROTATE_(ID)
        INTEGER, INTENT(IN) :: ID
        INTEGER :: Y, X
        X = ID
        Y = LEFT(X)
        LEFT(X) = RIGHT(Y)
        IF (RIGHT(Y) .NE. -1) THEN
            PARENT(RIGHT(Y)) = X
        END IF
        PARENT(Y) = PARENT(X)
        IF (PARENT(X) .EQ. -1) THEN
            ROOT = Y
        ELSEIF (X .EQ. RIGHT(PARENT(X))) THEN
            RIGHT(PARENT(X)) = Y
        ELSE
            LEFT(PARENT(X)) = Y
        END IF
        RIGHT(Y) = X
        PARENT(X) = Y
    END SUBROUTINE
    
    FUNCTION RB_COMPARE_(ID1, VALUE1, ID2, VALUE2) RESULT(ANS)
        INTEGER, INTENT(IN) :: ID1, VALUE1, ID2, VALUE2
        INTEGER :: ANS
        ANS = VALUE1 - VALUE2
        IF (ANS .EQ. 0) THEN
            ANS = ID1 - ID2
        END IF        
    END FUNCTION
    
    SUBROUTINE RB_INSERT_FIXUP_(ID)
        INTEGER, INTENT(IN) :: ID
        INTEGER :: Y, Z
        Z = ID
        DO WHILE (.NOT. COLOR(PARENT(Z)))
            IF (PARENT(Z) .EQ. LEFT(PARENT(PARENT(Z)))) THEN
                Y = RIGHT(PARENT(PARENT(Z)))
                IF (.NOT. COLOR(Y)) THEN
                    COLOR(PARENT(Z)) = .TRUE.
                    COLOR(Y) = .TRUE.
                    COLOR(PARENT(PARENT(Z))) = .FALSE.
                    Z = PARENT(PARENT(Z))
                ELSE 
                    IF (Z .EQ. RIGHT(PARENT(Z))) THEN
                        Z = PARENT(Z)
                        CALL RB_LEFT_ROTATE_(Z)
                    END IF
                    COLOR(PARENT(Z)) = .TRUE.
                    COLOR(PARENT(PARENT(Z))) = .FALSE.
                    CALL RB_RIGHT_ROTATE_(PARENT(PARENT(Z)))
                END IF
            ELSE
                Y = LEFT(PARENT(PARENT(Z)))
                IF (.NOT. COLOR(Y)) THEN
                    COLOR(PARENT(Z)) = .TRUE.
                    COLOR(Y) = .TRUE.
                    COLOR(PARENT(PARENT(Z))) = .FALSE.
                    Z = PARENT(PARENT(Z))
                ELSE 
                    IF (Z .EQ. LEFT(PARENT(Z))) THEN
                        Z = PARENT(Z)
                        CALL RB_RIGHT_ROTATE_(Z)
                    END IF
                    COLOR(PARENT(Z)) = .TRUE.
                    COLOR(PARENT(PARENT(Z))) = .FALSE.
                    CALL RB_LEFT_ROTATE_(PARENT(PARENT(Z)))
                END IF                
            END IF
        END DO
        COLOR(ROOT) = .TRUE.
    END SUBROUTINE
    
    SUBROUTINE RB_INSERT_(ID, VALUE)
        INTEGER, INTENT(IN) :: ID, VALUE
        INTEGER :: X, Y
        Y = -1
        X = ROOT
        DO WHILE (X .NE. -1)
            Y = X
            IF (RB_COMPARE_(ID, VALUE, X, POTENTIAL(X)) .LT. 0) THEN
                X = LEFT(X)
            ELSE
                X = RIGHT(X)
            END IF            
        END DO
        PARENT(ID) = Y
        IF (Y .EQ. -1) THEN
            ROOT = ID
        ELSEIF  (RB_COMPARE_(ID, VALUE, Y, POTENTIAL(Y)) .LT. 0) THEN
            LEFT(Y) = ID
        ELSE
            RIGHT(Y) = ID
        END IF
        LEFT(ID) = -1
        RIGHT(ID) = -1
        COLOR(ID) = .FALSE.
        POTENTIAL(ID) = VALUE
        CALL RB_INSERT_FIXUP_(ID)
        RB_SIZE = RB_SIZE + 1
    END SUBROUTINE
    
    SUBROUTINE RB_TRANSPLANT_(ID1, ID2)
        INTEGER, INTENT(IN) :: ID1, ID2
        INTEGER :: U, V
        U = ID1
        V = ID2
        IF (PARENT(U) .EQ. -1) THEN
            ROOT = V
        ELSEIF (U .EQ. LEFT(PARENT(U))) THEN
            LEFT(PARENT(U)) = V
        ELSE
            RIGHT(PARENT(U)) = V
        END IF
        PARENT(V) = PARENT(U)
    END SUBROUTINE
    
    SUBROUTINE RB_DELETE_FIXUP_(ID)
        INTEGER, INTENT(IN) :: ID
        INTEGER :: W, X, Y
        X = ID
        DO WHILE ((X .NE. ROOT) .AND. COLOR(X)) 
            IF (X .EQ. LEFT(PARENT(X))) THEN
                W = RIGHT(PARENT(X))
                IF (.NOT. COLOR(W)) THEN
                    COLOR(W) = .TRUE.
                    COLOR(PARENT(X)) = .FALSE.
                    CALL RB_LEFT_ROTATE_(PARENT(X))
                    W = RIGHT(PARENT(X))
                END IF
                IF (COLOR(LEFT(W)) .AND. COLOR(RIGHT(W))) THEN
                    COLOR(W) = .FALSE.
                    X = PARENT(X)
                ELSE
                    IF (COLOR(RIGHT(W))) THEN
                        COLOR(LEFT(W)) = .TRUE.
                        COLOR(W) = .FALSE.
                        CALL RB_RIGHT_ROTATE_(W)
                        W = RIGHT(PARENT(X))
                    END IF
                    COLOR(W) = COLOR(PARENT(X))
                    COLOR(PARENT(X)) = .TRUE.
                    COLOR(RIGHT(W)) = .TRUE.
                    CALL RB_LEFT_ROTATE_(PARENT(X))
                    X = ROOT
                END IF
            ELSE
                W = LEFT(PARENT(X))
                IF (.NOT. COLOR(W)) THEN
                    COLOR(W) = .TRUE.
                    COLOR(PARENT(X)) = .FALSE.
                    CALL RB_RIGHT_ROTATE_(PARENT(X))
                    W = LEFT(PARENT(X))
                END IF
                IF (COLOR(RIGHT(W)) .AND. COLOR(LEFT(W))) THEN
                    COLOR(W) = .FALSE.
                    X = PARENT(X)
                ELSE
                    IF (COLOR(LEFT(W))) THEN
                        COLOR(RIGHT(W)) = .TRUE.
                        COLOR(W) = .FALSE.
                        CALL RB_LEFT_ROTATE_(W)
                        W = LEFT(PARENT(X))
                    END IF
                    COLOR(W) = COLOR(PARENT(X))
                    COLOR(PARENT(X)) = .TRUE.
                    COLOR(LEFT(W)) = .TRUE.
                    CALL RB_RIGHT_ROTATE_(PARENT(X))
                    X = ROOT
                END IF
            END IF
        END DO
        COLOR(X) = .TRUE.
    END SUBROUTINE
    
    FUNCTION RB_MIN_(ID) RESULT(ANS)
        INTEGER, INTENT(IN) :: ID
        INTEGER :: ANS
        ANS = ID
        DO WHILE (LEFT(ANS) .NE. -1)
            ANS = LEFT(ANS)
        END DO
    END FUNCTION
    
    SUBROUTINE RB_DELETE_(ID)
        INTEGER, INTENT(IN) :: ID
        INTEGER :: X, Y, Z
        LOGICAL :: Y_ORIGINAL_COLOR
        Z = ID
        Y = Z
        Y_ORIGINAL_COLOR = COLOR(Y)
        IF (LEFT(Z) .EQ. -1) THEN
            X = RIGHT(Z)
            CALL RB_TRANSPLANT_(Z, RIGHT(Z))
        ELSEIF (RIGHT(Z) .EQ. -1) THEN
            X = LEFT(Z)
            CALL RB_TRANSPLANT_(Z, LEFT(Z))
        ELSE
            Y = RB_MIN_(RIGHT(Z))
            Y_ORIGINAL_COLOR = COLOR(Y)
            X = RIGHT(Y)
            IF (PARENT(Y) .EQ. Z) THEN
                PARENT(X) = Y
            ELSE
                CALL RB_TRANSPLANT_(Y, RIGHT(Y))
                RIGHT(Y) = RIGHT(Z)
                PARENT(RIGHT(Y)) = Y
            END IF
            CALL RB_TRANSPLANT_(Z, Y)
            LEFT(Y) = LEFT(Z)
            PARENT(LEFT(Y)) = Y
            COLOR(Y) = COLOR(Z)
        END IF
        IF (Y_ORIGINAL_COLOR) THEN
            CALL RB_DELETE_FIXUP_(X)
        END IF        
        RB_SIZE = RB_SIZE - 1
    END SUBROUTINE
        
    FUNCTION RB_MAX_(ID) RESULT(ANS)
        INTEGER, INTENT(IN) :: ID
        INTEGER :: ANS
        ANS = ID
        DO WHILE (RIGHT(ANS) .NE. -1)
            ANS = RIGHT(ANS)
        END DO
    END FUNCTION
    
    FUNCTION RB_TREE_MAX_() RESULT(ANS)
        INTEGER :: ANS
        ANS = RB_MAX_(ROOT)
    END FUNCTION
    
    FUNCTION RB_SUCCESSOR_(ID) RESULT(ANS)
        INTEGER, INTENT(IN) :: ID
        INTEGER :: ANS
        INTEGER :: X, Y
        X = ID
        IF (RIGHT(X) .NE. -1) THEN
            ANS = RB_MIN_(RIGHT(X))
            RETURN
        END IF
        Y = PARENT(X)
        DO WHILE ((Y .NE. -1) .AND. (X .EQ. RIGHT(Y))) 
            X = Y
            Y = PARENT(Y)
        END DO
        ANS = Y
    END FUNCTION
    
    RECURSIVE SUBROUTINE RB_ORDERED_PRINT_(ID) 
        INTEGER, INTENT(IN) :: ID
        INTEGER :: X
        X = ID
        IF (LEFT(X) .NE. -1) THEN
            CALL RB_ORDERED_PRINT_(LEFT(X))
        END IF
        WRITE (*, 11) X, PARENT(X), LEFT(X), RIGHT(X), COLOR(X), POTENTIAL(X)
        IF (RIGHT(X) .NE. -1) THEN
            CALL RB_ORDERED_PRINT_(RIGHT(X))
        END IF        
11      FORMAT(I5,':', I5, I5, I5, L5, I5)                
    END SUBROUTINE 
    
    
END MODULE