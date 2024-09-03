
!****************************************************************************************************************!
!*                                                                                                              *!    
!*                                                  BICG SOLVER                                                 *!
!*                                                                                                              *!    
!*                                                                                                              *!    
!****************************************************************************************************************!

    SUBROUTINE BICG_UPDATE_P__(SIZE, BETHA_I_1, OMEGA_I_1, P, R, V)
        !DEC$ ATTRIBUTES DLLEXPORT, ALIAS : 'BICG_UPDATE_P__' :: RETINA_SOLVER
        INTEGER, INTENT(IN) :: SIZE
        DOUBLE PRECISION, INTENT(IN) :: BETHA_I_1, OMEGA_I_1
        DOUBLE PRECISION, DIMENSION(0:SIZE - 1), INTENT(IN) :: R, V
        DOUBLE PRECISION, DIMENSION(0:SIZE - 1), INTENT(INOUT) :: P
        INTEGER :: I
        !DEC$ SIMD
        DO I = 0, SIZE - 1
            P(I) = R(I) + BETHA_I_1 * (P(I) - OMEGA_I_1 * V(I))
        END DO        
    END SUBROUTINE

    SUBROUTINE BICG_UPDATE_R__(SIZE, OMEGA_I, R, S, T)
        !DEC$ ATTRIBUTES DLLEXPORT, ALIAS : 'BICG_UPDATE_R__' :: RETINA_SOLVER
        INTEGER, INTENT(IN) :: SIZE
        DOUBLE PRECISION, INTENT(IN) :: OMEGA_I
        DOUBLE PRECISION, DIMENSION(0:SIZE - 1), INTENT(IN) :: S, T
        DOUBLE PRECISION, DIMENSION(0:SIZE - 1), INTENT(OUT) :: R
        INTEGER :: I
        !DEC$ SIMD
        DO I=0, SIZE - 1
            R(I) = S(I) - OMEGA_I * T(I)
        END DO        
    END SUBROUTINE
    
    SUBROUTINE BICG_CALC_X__(SIZE, ALPHA_I, OMEGA_I, X_RESULT, P_HAT, S_HAT)
        !DEC$ ATTRIBUTES DLLEXPORT, ALIAS : 'BICG_CALC_X__' :: RETINA_SOLVER
        INTEGER, INTENT(IN) :: SIZE
        DOUBLE PRECISION, INTENT(IN) :: ALPHA_I, OMEGA_I
        DOUBLE PRECISION, DIMENSION(0:SIZE - 1), INTENT(IN) :: P_HAT, S_HAT
        DOUBLE PRECISION, DIMENSION(0:SIZE - 1), INTENT(INOUT) :: X_RESULT
        INTEGER :: I
        !DEC$ SIMD
        DO I=0, SIZE - 1
            X_RESULT(I) = X_RESULT(I) + ALPHA_I * P_HAT(I) + OMEGA_I * S_HAT(I)
        END DO
    END SUBROUTINE

    SUBROUTINE BICG_CALC_X_WITHOUT_S__(SIZE, ALPHA_I, X_RESULT, P_HAT)
        !DEC$ ATTRIBUTES DLLEXPORT, ALIAS : 'BICG_CALC_X_WITHOUT_S__' :: RETINA_SOLVER
        INTEGER, INTENT(IN) :: SIZE
        DOUBLE PRECISION, INTENT(IN) :: ALPHA_I
        DOUBLE PRECISION, DIMENSION(0:SIZE - 1), INTENT(IN) :: P_HAT
        DOUBLE PRECISION, DIMENSION(0:SIZE - 1), INTENT(INOUT) :: X_RESULT
        INTEGER :: I
        !DEC$ SIMD
        DO I=0, SIZE - 1
            X_RESULT(I) = X_RESULT(I) + ALPHA_I * P_HAT(I)
        END DO
    END SUBROUTINE

    SUBROUTINE BICG_UPDATE_X__(SIZE, RESULT_VALS, X_RESULT)
        !DEC$ ATTRIBUTES DLLEXPORT, ALIAS : 'BICG_UPDATE_X__' :: RETINA_SOLVER
        INTEGER, INTENT(IN) :: SIZE
        DOUBLE PRECISION, DIMENSION(0:SIZE - 1), INTENT(IN) :: X_RESULT
        DOUBLE PRECISION, DIMENSION(0:SIZE - 1), INTENT(OUT) :: RESULT_VALS
        INTEGER :: I
        !DEC$ SIMD
        DO I=0, SIZE - 1
            RESULT_VALS(I) = X_RESULT(I)
        END DO
    END SUBROUTINE
    
    SUBROUTINE BICG_UPDATE_S__(SIZE, ALPHA_I, S, R, V)
        !DEC$ ATTRIBUTES DLLEXPORT, ALIAS : 'BICG_UPDATE_S__' :: RETINA_SOLVER
        INTEGER, INTENT(IN) :: SIZE
        DOUBLE PRECISION, INTENT(IN) :: ALPHA_I
        DOUBLE PRECISION, DIMENSION(0:SIZE - 1), INTENT(IN) :: R, V
        DOUBLE PRECISION, DIMENSION(0:SIZE - 1), INTENT(OUT) :: S
        INTEGER :: I
        !DEC$ SIMD
        DO I=0, SIZE - 1
            S(I) = R(I) - ALPHA_I * V(I)
        END DO        
    END SUBROUTINE

    