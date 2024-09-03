
!****************************************************************************************************************!
!*                                                                                                              *!    
!*                                                ORTHOMIN SOLVER                                               *!
!*                                                                                                              *!    
!*                                                                                                              *!    
!****************************************************************************************************************!

    SUBROUTINE ORTHOMIN_UPDATE_R__(SIZE, ALPHA, R, AD)
        !DEC$ ATTRIBUTES DLLEXPORT, ALIAS : 'ORTHOMIN_UPDATE_R__' :: RETINA_SOLVER
        INTEGER, INTENT(IN) :: SIZE
        DOUBLE PRECISION, INTENT(IN) :: ALPHA
        DOUBLE PRECISION, DIMENSION(0:SIZE - 1), INTENT(IN) :: AD
        DOUBLE PRECISION, DIMENSION(0:SIZE - 1), INTENT(OUT) :: R
        INTEGER :: I
        !DEC$ SIMD
        DO I = 0, SIZE - 1
            R(I) = R(I) - ALPHA * AD(I)
        END DO        
    END SUBROUTINE
    
    SUBROUTINE ORTHOMIN_UPDATE_X__(SIZE, ALPHA, X_RESULT, D)
        !DEC$ ATTRIBUTES DLLEXPORT, ALIAS : 'ORTHOMIN_UPDATE_X__' :: RETINA_SOLVER
        INTEGER, INTENT(IN) :: SIZE
        DOUBLE PRECISION, INTENT(IN) :: ALPHA
        DOUBLE PRECISION, DIMENSION(0:SIZE - 1), INTENT(IN) :: D
        DOUBLE PRECISION, DIMENSION(0:SIZE - 1), INTENT(INOUT) :: X_RESULT
        INTEGER :: I
        !DEC$ SIMD
        DO I = 0, SIZE - 1
            X_RESULT(I) = X_RESULT(I) + ALPHA * D(I)
        END DO
    END SUBROUTINE

    SUBROUTINE ORTHOMIN_UPDATE_BETAS__(SIZE, M, ITER_COUNT, BETAS, ADS, AZ)
        !DEC$ ATTRIBUTES DLLEXPORT, ALIAS : 'ORTHOMIN_UPDATE_BETAS__' :: RETINA_SOLVER
        USE CRS
        INTEGER, INTENT(IN) :: SIZE, M, ITER_COUNT
        DOUBLE PRECISION, DIMENSION(0:M - 1), INTENT(OUT) :: BETAS
        DOUBLE PRECISION, DIMENSION(0:SIZE - 1), INTENT(IN) :: AZ
        DOUBLE PRECISION, DIMENSION(0:M * SIZE - 1), INTENT(IN) :: ADS
        INTEGER :: M_INDEX, C0, M_NUM
        DOUBLE PRECISION :: DENOM
        IF (ITER_COUNT .LT. M) THEN
            M_NUM = ITER_COUNT
        ELSE
            M_NUM = M
        END IF
        !$OMP PARALLEL 
            !$OMP DO SCHEDULE(STATIC) PRIVATE(M_INDEX, C0, DENOM)
        DO M_INDEX = 0, M_NUM - 1
            C0 = M_INDEX * SIZE
            DENOM = CRS_INNER_PRODUCT_(SIZE, ADS(C0:C0 + SIZE - 1), ADS(C0:C0 + SIZE - 1))
            IF (DENOM .LT. 1.D-20) THEN
                BETAS(M_INDEX) = 0.D0
            ELSE
                BETAS(M_INDEX) = -CRS_INNER_PRODUCT_(SIZE, AZ, ADS(C0:C0 + SIZE - 1)) / DENOM
            END IF
        END DO
            !$OMP END DO
        !$OMP END PARALLEL
    END SUBROUTINE

    SUBROUTINE ORTHOMIN_UPDATE_DIRECTION_VECS__(SIZE, M, ITER_COUNT, PTR, Z, AZ, BETAS, DS, ADS)
        !DEC$ ATTRIBUTES DLLEXPORT, ALIAS : 'ORTHOMIN_UPDATE_DIRECTION_VECS__' :: RETINA_SOLVER
        INTEGER, INTENT(IN) :: SIZE, M, PTR, ITER_COUNT
        DOUBLE PRECISION, DIMENSION(0:M - 1), INTENT(IN) :: BETAS
        DOUBLE PRECISION, DIMENSION(0:SIZE - 1), INTENT(IN) :: Z, AZ
        DOUBLE PRECISION, DIMENSION(0:M * SIZE - 1), INTENT(INOUT) :: DS, ADS
        INTEGER :: M_INDEX, I, II, C0, M_NUM
        DOUBLE PRECISION :: SUM_D, SUM_AD
        DOUBLE PRECISION :: BETA
        IF (ITER_COUNT .LT. M) THEN
            M_NUM = ITER_COUNT
        ELSE
            M_NUM = M
        END IF
        !$OMP PARALLEL 
            !$OMP DO SCHEDULE(STATIC) PRIVATE(I, SUM_D, SUM_AD, M_INDEX, C0, BETA, II)
        DO I = 0, SIZE - 1
            SUM_D = Z(I)
            SUM_AD = AZ(I)
            !DEC$ SIMD
            DO M_INDEX = 0, M_NUM - 1
                C0 = M_INDEX * SIZE
                BETA = BETAS(M_INDEX)
                SUM_D = SUM_D + BETA * DS(C0 + I)
                SUM_AD = SUM_AD + BETA * ADS(C0 + I)
            END DO
            C0 = PTR * SIZE
            DS(C0 + I) = SUM_D
            ADS(C0 + I) = SUM_AD
        END DO
            !$OMP END DO
        !$OMP END PARALLEL
    END SUBROUTINE

    SUBROUTINE ORTHOMIN_UPDATE_DIRECTION_VEC__(SIZE, M, ITER_COUNT, PTR, INIT_VEC, BETAS, VECS)
        !DEC$ ATTRIBUTES DLLEXPORT, ALIAS : 'ORTHOMIN_UPDATE_DIRECTION_VEC__' :: RETINA_SOLVER
        INTEGER, INTENT(IN) :: SIZE, M, PTR, ITER_COUNT
        DOUBLE PRECISION, DIMENSION(0:M - 1), INTENT(IN) :: BETAS
        DOUBLE PRECISION, DIMENSION(0:SIZE - 1), INTENT(IN) :: INIT_VEC
        DOUBLE PRECISION, DIMENSION(0:M * SIZE - 1), INTENT(INOUT) :: VECS
        INTEGER :: M_INDEX, I, C0, M_NUM
        DOUBLE PRECISION :: SUM
        IF (ITER_COUNT .LT. M) THEN
            M_NUM = ITER_COUNT
        ELSE
            M_NUM = M
        END IF
        !$OMP PARALLEL 
            !$OMP DO SCHEDULE(STATIC) PRIVATE(I, C0, SUM, M_INDEX)
        DO I = 0, SIZE - 1
            SUM = INIT_VEC(I)
            !DEC$ SIMD
            DO M_INDEX = 0, M_NUM - 1
                C0 = M_INDEX * SIZE
                SUM = SUM + BETAS(M_INDEX) * VECS(C0 + I)
            END DO
            C0 = PTR * SIZE
            VECS(C0 + I) = SUM
        END DO
            !$OMP END DO
        !$OMP END PARALLEL

    END SUBROUTINE

    