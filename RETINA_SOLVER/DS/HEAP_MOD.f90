MODULE HEAP_MOD
        
    IMPLICIT NONE
    
    INTEGER :: HEAP_SIZE
    INTEGER, DIMENSION(1:99999) :: HEAP
    
    CONTAINS
    
    SUBROUTINE HEAP_CLEAR_()
        HEAP_SIZE = 0
    END SUBROUTINE
    
    RECURSIVE SUBROUTINE MIN_HEAPIFY_(ID)
        INTEGER, INTENT(IN) :: ID
        INTEGER :: I, L, R, SMALLEST, TEMP
        I = ID
        L = 2 * I
        R = 2 * I + 1
        IF ((L .LE. HEAP_SIZE) .AND. (HEAP(L) .LT. HEAP(I))) THEN
            SMALLEST = L
        ELSE
            SMALLEST = I
        END IF
        IF ((R .LE. HEAP_SIZE) .AND. (HEAP(R) .LT. HEAP(SMALLEST))) THEN
            SMALLEST = R
        END IF
        IF (SMALLEST .NE. I) THEN
            TEMP = HEAP(I)
            HEAP(I) = HEAP(SMALLEST)
            HEAP(SMALLEST) = TEMP
            CALL MIN_HEAPIFY_(SMALLEST)
        END IF
    END SUBROUTINE
    
    FUNCTION HEAP_EXTRACT_MIN_() RESULT(ANS)
        INTEGER :: ANS                
        IF (HEAP_SIZE .LT. 1) THEN
            ANS = -1
            RETURN
        END IF
        ANS = HEAP(1)
        HEAP(1) = HEAP(HEAP_SIZE)
        HEAP_SIZE = HEAP_SIZE - 1
        CALL MIN_HEAPIFY_(1)
    END FUNCTION
    
    SUBROUTINE HEAP_DECREASE_(ID, VALUE)
        INTEGER, INTENT(IN) :: ID, VALUE
        INTEGER :: I, PI, TEMP
        I = ID
        IF (VALUE .LT. HEAP(I)) THEN
            HEAP(I) = VALUE
            DO WHILE (I .GT. 1)
                PI = I / 2
                IF (HEAP(PI) .LE. HEAP(I)) THEN
                    EXIT
                END IF
                TEMP = HEAP(PI)
                HEAP(PI) = HEAP(I)
                HEAP(I) = TEMP
                I = PI
            END DO
        END IF
    END SUBROUTINE
    
    SUBROUTINE HEAP_INSERT_(VALUE)
        INTEGER, INTENT(IN) :: VALUE
        HEAP_SIZE = HEAP_SIZE + 1
        HEAP(HEAP_SIZE) = 2000000000
        CALL HEAP_DECREASE_(HEAP_SIZE, VALUE)
    END SUBROUTINE
    
    SUBROUTINE HEAP_PRINT_()
        DO WHILE (HEAP_SIZE .GT. 0) 
            WRITE (*,*) HEAP_EXTRACT_MIN_()
        END DO
    END SUBROUTINE
    
    
END MODULE