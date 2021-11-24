!***********************************************************!
!                                                           !
!   NOT USED IN CURRENT PROGRAM, MODIFY FOR STATESPACE EQ.  !
!   OR LCO ANALYSIS ON FREEPLAY STRUCTURES LATER ON         !
!                                                           !
!***********************************************************!
    
    SUBROUTINE STATESPACE
    
    use statespace_var
    use aerodynamic_force
    USE DLM_var
    USE input_var
    USE panel
    
    IMPLICIT NONE
    
    INTEGER W2, W3,W4,W5,PP
    DOUBLE PRECISION :: DESIRED_DAMPING = 0
    
    allocate(Ans(5,10,10), SW(FEM_X, FEM_Y-1), SC(FEM_X,1))
    allocate(Kw(FEM_y-1,FEM_y-1))
    ALLOCATE(RAYLEIGH(FEM_Y-1,2), FREQ(FEM_Y-1), DAMPINGMR(FEM_Y-1,1))
    ALLOCATE(MW(FEM_Y-1,FEM_Y-1), KC(1,1), MC(1,1),K(FEM_Y,FEM_Y))
    
    ANS = AN_
    
    ! GET NORMAL MODES (WINGS)
    
    DO W2 = 1, FEM_Y-1
        SW(:,W2) = S(:,W2)
    END DO
    
    SC(:,1) = S(:,FEM_Y)
    
    ! COMPUTE MODAL STIFFNESS FOR NORMAL MODES
    
    !CALCULATION OF STIFFNESS MATRIX (GSTIFF) AND MASS MATRIX, (GMASS)
    !NEEDS TO BE PERFORMED BEFORE HAND SEE ANALYSIS.F90
    
    KW = MATMUL(MATMUL(TRANSPOSE(SW),GSTIFF1),SW)
    KC = MATMUL(MATMUL(TRANSPOSE(SC),GSTIFF2),SC)
    
    MW = MATMUL(MATMUL(TRANSPOSE(SW),GMASS1),SW)
    MC = MATMUL(MATMUL(TRANSPOSE(SC),GMASS2),SC)
    
    ! MODAL DAMPING COMPUTATIONS
    !BASED ON RALEIGH FORMULATION
    
    DO PP = 1, FEM_Y-1
        FREQ(PP) = SQRT(K(PP,PP))
        RAYLEIGH(PP, 1) = 1/2/FREQ(PP)
        RAYLEIGH(PP,2) = FREQ(PP)/2
        DAMPINGMR(PP,1) = DESIRED_DAMPING;
    END DO
    
    !COEFFICIENTS = 
    
    END SUBROUTINE STATESPACE
    