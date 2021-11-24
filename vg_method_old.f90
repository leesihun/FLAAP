!***********************************************************!
!                                                           !
!   NOT USED IN CURRENT PROGRAM, USED FOR OLDER....         !
!                                                           !
!***********************************************************!
    
    subroutine vgplot
    
        USE aerodynamic_force
    USE ATMOS
    USE CLU
    USE CLU_D
    USE DLM_var
    USE input_var
    USE panel
    USE statespace_Var
    ! USE STATEMENTS.
    
    implicit none
    
    DOUBLE PRECISION :: KRLOW, KRHIGH, INCRE_KR, V_LOW, V_HIGH, INCRE_V
    DOUBLE COMPLEX :: SMALL_S1, small_s2, cmplx, temp_c, small_s
    integer :: errorcode, w3, ITER, w1, w2, marker
    double precision :: kr_
    double precision,allocatable :: KT_inv(:,:)
    double complex, allocatable :: Q_1(:,:), A_1(:,:)
    double precision, allocatable ::A1_re(:,:), A1_im(:,:)
    double precision, allocatable :: VR(:,:), Vi(:,:),  V(:,:)
    double precision, allocatable :: D_(:,:),Dr(:,:), Di(:,:)
    LOGICAL :: AVEC
    double complex, allocatable :: a1_(:,:)
    double complex :: F1, F2 !det(A-pI)
    double complex :: temp1(10,10), EIG_P_TEMP
    double complex :: unity(10, 10)
    double complex :: EIG_PTEMP1, EIG_PTEMP2    ! eigen value, p '
    INTEGER PRINT_, natf!natural frequency
    double precision :: eps = 1.E-10
    integer :: ten = 10
    double precision :: velocity_old
    double precision :: damping(10, 1000)
    double precision :: flutt_freq(10, 1000)
    
    integer :: count
    !logical avec
    real*8:: v_r(10,10), V_i(10,10),d_r(10,10),d_i(10,10),work(20,20)
    ! DECLARIATIONS.
    
    
    DOUBLE PRECISION :: WORKINGSPACE(1000000)
    ! WORKING SPACE FOR CEGIEN, MIGHT NEED TO HAVE MORE SPACE....
    n_mode = 10;
    
    allocate(kt_inv(fem_y, fem_y))
    ALLOCATE (KT(FEM_Y,FEM_Y), MT(FEM_Y,FEM_Y))
    allocate (Q_1(FEM_y,fem_y),A_1(FEM_y,fem_y)       )
    ALLOCATE(A1_(FEM_Y, FEM_Y),A1_re(FEM_Y, FEM_Y), A1_im(FEM_Y, FEM_Y))
    allocate(VR(FEM_Y, FEM_Y), Vi(FEM_Y, FEM_Y), Dr(FEM_Y, FEM_Y), Di(FEM_Y, FEM_Y), V(FEM_Y, FEM_Y), D_(FEM_Y, FEM_Y))
    ALLOCATE(Q_1_1(FEM_Y, FEM_Y),Q_1_2(FEM_Y, FEM_Y))
    allocate(A_eig1(20,20),A_eig2(20,20))  ! A matrix, refer to "Computing Sequences, P-K method fluttr"
    allocate(P_K_results(10, 100, 100,100))
    allocate(eig_pc(1000))
    P_K_results = 0;
    !ALLOCATION DONE
    v_low = 30
    incre_v = 0.5
    v_high = 100
    
    eig_p1 = 0
    eig_p2 = -100;
    MAC = length_y_total!Mean Aerodynamic Chord
    
    
    unity = 0
    do w1 = 1, 10
        unity(w1, w1) = 1;
    end do
    
    ! REDUCED FREQENCY SPECIFICATION
    !****************************************************************************************!
    !do w3 = 1, nmax
    !    beta(w3) = 1.7*(reduced_Freq_max)*w3/((nmax+1)**2);
    !end do
    ! BETA, beta=1.7*kmax*n/(nmax+1)^2, where kmax is the max frequency of interest
    KRLOW = 0.0001
    KRHIGH = 30
    INCRE_KR = 0.0001
    ! FREQUENCY LIMITS, FROM 0.0001 TO 30HZ
    !****************************************************************************************!
    
    
    !MAKE GENERALIZED MODAL MATRICES
    !****************************************************************************************!
    KT = MATMUL(MATMUL(TRANSPOSE(S), GSTIFF1),S) 
    ! GSTIFF IS GLOBAL NODAL STIFFNESS MATRIX FROM NASTRAN
    ! MATRIX SIZE : 10X10 matrix(10MODES)
    
    MT = MATMUL(MATMUL(TRANSPOSE(S), GMASS1),S)
    ! GMASS IS GLOBAL NODAL MASS MATRIX FROM NASTRAN
    ! MATRIX SIZE : 10X10 matrix(10MODES)
    !****************************************************************************************!
    
    
    ITER = 1;
    MACH = 0;
    
    KT(:,:) = GSTIFF(1,:,:)
    MT(:,:) = GMASS(1,:,:)
    
    iter = 1;
    Mach = 0;
    
    do kr_ = krlow, krhigh, incre_kr
        pause;
        
        small_s = cmplx(0, 1)*kr_
        pause;
        Q_1 = (AN_(1,:,:)+AN_(2,:,:)*SMALL_S+AN_(3,:,:)*(SMALL_S**2)+small_s*AN_(4,:,:)/(SMALL_S+BETA(1))+small_s*AN_(5,:,:)/(SMALL_S+BETA(2))) ! ONLY FOR NMAX == 2
        
        call CINVERT_LU_D(KT, ten, KT_INV, errorcode)
        
        A_1 = MATMUL(KT_inv, MATMUL(Mt-rho*(MAC/2)**2/(kr_**2)/2, Q_1))
        
        print*, 'KT : ', KT;pause;
        print*, 'KT_inv : ', KT_inv;pause;
        print*, 'KM : ', mt;pause;
        print*, 'rho : ', rho;pause;
        print*, 'kr_ : ', kr_;pause;
        print*, 'Q1 : ', q_1;pause;
        
        V_r = real(A_1)
        V_i = imag(A_1)
        
        call COMEIG(10,10,100000,v_r,V_i,d_r,d_i,errorcode,work)
        
        D = cmplx(d_r, d_i)
        
        if(iter>1) then
            do w2 = 1, 10
                
                
                
            end do
            !lambda = d1
        else
            !lambda = D1
        end if
        
        !freq = 1/sqrt(real(lambda))/2/3.141592
        !damping = 
        !call COMEIG(real(A_1),imag(A_1),10,10,AVEC,V_r,V_i,d_r,d_i,work,errorcode)
    end do
    
    end subroutine vgplot
    