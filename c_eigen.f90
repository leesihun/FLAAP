
      REAL*8 FUNCTION DABSCD(AR,AI)
!     ABSOLUTE VALUE OF A COMPLEX NUMBER C=AR+I*AI
!     DABSCD=DSQRT(AR**2+AI**2)
      REAL*8 AR,AI,XR,XI,W
      XR=DABS(AR)
      XI=DABS(AI)
      IF(XR.LE.XI) THEN
      W=XR
      XR=XI
      XI=W
      ENDIF
      IF(XI.EQ.0.D0) THEN
      DABSCD=XR
      ELSE
      DABSCD=XR*DSQRT(1.D0+(XI/XR)**2)
      ENDIF
      RETURN
      END

      SUBROUTINE DIVCD(AR,AI,BR,BI,ZR,ZI)
!     COMPLEX DIVISION Z=ZR+I*ZI=(AR+I*AI)/(BR+I*BI)
!     DO NOT USE IF BR=BI=0.
      REAL*8 AR,AI,BR,BI,YR,YI,ZR,ZI,W
      YR=BR
      YI=BI
      IF(DABS(YR).GT.DABS(YI)) THEN
      W=YI/YR
      YR=W*YI+YR
      ZR=(AR+W*AI)/YR
      ZI=(AI-W*AR)/YR
      ELSE
      W=YR/YI
      YI=W*YR+YI
      ZR=(W*AR+AI)/YI
      ZI=(W*AI-AR)/YI
      ENDIF
      RETURN
      END

      SUBROUTINE COMEIG (NDIM,N,NN,A,Z,T,U,IER,EN)
!--------------------------------------------------------------------------------
!     THIS SUBROUTINE CALCULATES THE EIGENVALUES/EIGENVECTORS OF A COMPLEX MATRIX 
!     C = A+I*Z BY THE JACOBI METHOD, THIS METHOD IS ALSO RECOMMANDED TO DEAL WITH 
!     REAL MATRICES THE EIGENVALUES OF WHICH ARE COMPLEX.

!     DATA:
!     NDIM    1ST DIMENSION OF TABLES A, Z, T, U IN MAIN PROGRAM (HERE NDIM=N)
!     N       REAL SIZE OF COMPLEX MATRIX C = A+I*Z
!     NN      MAXIMUM NUMBER OF ITERATIONS
!     A       TABLE STORING THE REAL PART OF GIVEN MATRIX
!     Z       TABLE STORING THE IMAGINARY PART OF GIVEN MATRIX

!     OUTPUTS:
!     A(J,J),Z(J,J),J=1,N   IN MAIN DIAGONALS OF TABLES A AND Z, YOU
!             HAVE NOW RESPECTIVELY THE REAL AND IMAGINARY PARTS OF EIGENVALUES.
!     T,U     THESE TABLES CONTAIN NOW RESPECTIVELY THE REAL AND IMAGINARY PARTS
!             OF THE EIGENVECTORS MATRIX X = T+I*U (STORED IN COLUMNS).
!     IER     ERROR CODE
!             = 0  CONVERGENCE OK
!             = 1  NO CONVERGENCE AFTER NN ITERATIONS

!     WORKING ZONE:
!     EN      TABLE OF SIZE 2*N

!     NOTES:
!     1/      IN CASE OF CONVERGENCE (IER = 0), THE MATRIX EQUATION  C*X = LAMDA*X
!             IS VERIFIED TO THE MACHINE PRECISION.
!
!     REFERENCE:
!     P.J.EBERLEIN, NUMER.MATH 14, PP 232-245 (1970)
!---------------------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(NDIM,*),Z(NDIM,*),T(NDIM,*),U(NDIM,*),EN(*)
      LOGICAL MARK

!     CHECK MACHINE EPSILON (AROUND 1.2E-16 FOR PC)

      EPS = 1.0
   10 EPS = 0.5*EPS
      EPS1 = EPS+1.0
      IF (EPS1.GT.1.0) GO TO 10
      MARK = .FALSE.

!     INITIALIZE EIGENVECTORS

      DO I = 1,N
      T(I,I) = 1.0
      U(I,I) = 0.0
        DO J = I+1,N
        T(I,J) = 0.0
        T(J,I) = 0.0
        U(I,J) = 0.0
        U(J,I) = 0.0
        ENDDO
      ENDDO
      IT = 0
   20 IT = IT+1

!     SAFETY TEST IN CASE OF NO CONVERGENCE

      IF (IT.GT.NN) GO TO 90
      IF (MARK)     GO TO 95

!     DEFINE CONVERGENCE CRITERIUM

      TAU = 0.0
      DO K = 1,N
      W1 = 0.0
        DO I = 1,N
        IF (I.NE.K) W1 = W1+DABS(A(I,K))+DABS(Z(I,K))
        ENDDO
      TAU = TAU+W1
      EN(K) = W1+DABS(A(K,K))+DABS(Z(K,K))
      ENDDO

!     PERMUTE  LINES AND COLUMNS

      DO K = 1,N-1
      EMAX = EN(K)
      I = K
        DO J = K+1,N
        IF (EN(J).GT.EMAX) THEN
        EMAX = EN(J)
        I = J
        ENDIF
        ENDDO
      IF (I.NE.K) THEN
      EN(I) = EN(K)
        DO J = 1,N
        W2 = A(K,J)
        A(K,J) = A(I,J)
        A(I,J) = W2
        W2 = Z(K,J)
        Z(K,J) = Z(I,J)
        Z(I,J) = W2
        ENDDO
        DO J = 1,N
        W2 = A(J,K)
        A(J,K) = A(J,I)
        A(J,I) = W2
        W2 = Z(J,K)
        Z(J,K) = Z(J,I)
        Z(J,I) = W2
        W2 = T(J,K)
        T(J,K) = T(J,I)
        T(J,I) = W2
        W2 = U(J,K)
        U(J,K) = U(J,I)
        U(J,I) = W2
        ENDDO
      END IF
      ENDDO

!     CONVERGENCE IF TAU < 100*EPS

      IF (TAU.LT.100.0*EPS) GO TO 95

!     BEGIN ITERATIONS

      MARK = .TRUE.
      DO K = 1,N-1
      DO M = K+1,N
      G = 0.0
      HR = 0.0
      HJ = 0.0
      HI = 0.0
        DO I = 1,N

        IF (I.NE.K.AND.I.NE.M) THEN

        HR = HR+A(K,I)*A(M,I)+Z(K,I)*Z(M,I)-A(I,K)*A(I,M)-Z(I,K)*Z(I,M)
        HI = HI+Z(K,I)*A(M,I)-A(K,I)*Z(M,I)-A(I,K)*Z(I,M)+Z(I,K)*A(I,M)
        T1 = A(I,K)*A(I,K)+Z(I,K)*Z(I,K)+A(M,I)*A(M,I)+Z(M,I)*Z(M,I)
        T2 = A(I,M)*A(I,M)+Z(I,M)*Z(I,M)+A(K,I)*A(K,I)+Z(K,I)*Z(K,I)
        G = G+T1+T2
        HJ = HJ-T1+T2
        ENDIF
        ENDDO
      BR = A(K,M)+A(M,K)
      BI = Z(K,M)+Z(M,K)
      ER = A(K,M)-A(M,K)
      EI = Z(K,M)-Z(M,K)
      DR = A(K,K)-A(M,M)
      DI = Z(K,K)-Z(M,M)
      T1 = BR*BR+EI*EI+DR*DR
      T2 = BI*BI+ER*ER+DI*DI

      IF (T1.GE.T2) THEN

      SW = 1.0
      C = BR
      S = EI
      D = DR
      DE = DI
      ROOT2 = SQRT(T1)
      ELSE
      SW =-1.0
      C = BI
      S =-ER
      D = DI
      DE = DR
      ROOT2 = SQRT(T2)
      ENDIF
      ROOT1 = SQRT(S*S+C*C)
      SIG = 1.0
      IF (D.LT.0.0) SIG = -1.0
      CA = 1.0
      IF (C.LT.0.0) CA = -1.0
      SA = 0.0

      IF (ROOT1.LT.EPS) THEN

      SX = 0.0
      SA = 0.0
      CX = 1.0
      CA = 1.0

      IF (SW.GT.0.0) THEN
      E = ER
      B = BI
      ELSE
      E = EI
      B =-BR
      ENDIF
      DN = D*D+DE*DE
      GO TO 65
      ENDIF

      IF (DABS(S).GT.EPS) THEN

      CA = C/ROOT1
      SA = S/ROOT1
      ENDIF
      COT2X = D/ROOT1
      COTX = COT2X+SIG*SQRT(1.0+COT2X*COT2X)
      SX = SIG/SQRT(1.0+COTX*COTX)
      CX = SX*COTX
      ETA = (ER*BR+BI*EI)/ROOT1
      TSE = (BR*BI-ER*EI)/ROOT1
      T1 = SIG*(TSE*D-ROOT1*DE)/ROOT2
      T2 = (D*DE+ROOT1*TSE)/ROOT2
      DN = ROOT2*ROOT2+T2*T2
      T2 = HJ*CX*SX
      COS2A = CA*CA-SA*SA
      SIN2A = 2.0*CA*SA
      W1 = HR*COS2A+HI*SIN2A
      W2 = HI*COS2A-HR*SIN2A
      HR = CX*CX*HR-SX*SX*W1-CA*T2
      HI = CX*CX*HI+SX*SX*W2-SA*T2
      B = SW*T1*CA+ETA*SA
      E = CA*ETA-SW*T1*SA

!     ROOT1 < EPS

   65 S = HR-SIG*ROOT2*E
      C = HI-SIG*ROOT2*B
      ROOT = DSQRT(C*C+S*S)

      IF (ROOT.LT.EPS) THEN

      CB = 1.0
      CH = 1.0
      SB = 0.0
      SH = 0.0
      GO TO 70
      END IF
      CB = -C/ROOT
      SB =  S/ROOT
      T2 = CB*B-E*SB
      CN = T2*T2
      TANH = ROOT/(G+2.0*(CN+DN))
      CH = 1.0/DSQRT(1.0-TANH*TANH)
      SH = CH*TANH

!     ROOT < EPS

   70 W1 = SX*SH*(SA*CB-SB*CA)
      C1R = CX*CH-W1
      C2R = CX*CH+W1
      C1I =-SX*SH*(CA*CB+SA*SB)
      C2I = C1I
      W2 = SX*CH*CA
      W1 = CX*SH*SB
      S1R = W2-W1
      S2R =-W2-W1
      W2 = SX*CH*SA
      W1 = CX*SH*CB
      S1I = W2+W1
      S2I = W2-W1
      W1 = SQRT(S1R*S1R+S1I*S1I)
      W2 = SQRT(S2R*S2R+S2I*S2I)

      IF (W1.GT.EPS.OR.W2.GT.EPS) THEN

!     BEGIN TRANSFORMATIONS

      MARK = .FALSE.
      DO I = 1,N
      AKI = A(K,I)
      AMI = A(M,I)
      ZKI = Z(K,I)
      ZMI = Z(M,I)
      A(K,I) = C1R*AKI-C1I*ZKI+S1R*AMI-S1I*ZMI
      Z(K,I) = C1R*ZKI+C1I*AKI+S1R*ZMI+S1I*AMI
      A(M,I) = S2R*AKI-S2I*ZKI+C2R*AMI-C2I*ZMI
      Z(M,I) = S2R*ZKI+S2I*AKI+C2R*ZMI+C2I*AMI
      ENDDO
      DO I = 1,N
      AIK = A(I,K)
      AIM = A(I,M)
      ZIK = Z(I,K)
      ZIM = Z(I,M)
      TIK = T(I,K)
      TIM = T(I,M)
      UIK = U(I,K)
      UIM = U(I,M)
      A(I,K) = C2R*AIK-C2I*ZIK-S2R*AIM+S2I*ZIM
      Z(I,K) = C2R*ZIK+C2I*AIK-S2R*ZIM-S2I*AIM
      A(I,M) = C1R*AIM-C1I*ZIM-S1R*AIK+S1I*ZIK
      Z(I,M) = C1R*ZIM+C1I*AIM-S1R*ZIK-S1I*AIK
      T(I,K) = C2R*TIK-C2I*UIK-S2R*TIM+S2I*UIM
      U(I,K) = C2R*UIK+C2I*TIK-S2R*UIM-S2I*TIM
      T(I,M) = C1R*TIM-C1I*UIM-S1R*TIK+S1I*UIK
      U(I,M) = C1R*UIM+C1I*TIM-S1R*UIK-S1I*TIK
      ENDDO
      ENDIF

!     END TRANSFORMATIONS

      ENDDO
      ENDDO

!     GO TO NEXT ITERATION

      GO TO 20

!     NO CONVERGENCE !

   90 IER = 1
      RETURN

!     CONVERGENCE OK

   95 IER = 0

!     SORT SOLUTIONS IN INCREASING ORDER
      DO J=2,N
      VR=A(J,J)
      VI=Z(J,J)
      DO K=1,N
      EN(K)=T(K,J)
      EN(K+N)=U(K,J)
      ENDDO
         DO I=J-1,1,-1
         IF(DABSCD(A(I,I),Z(I,I)).LE.DABSCD(VR,VI)) GO TO 97
         A(I+1,I+1)=A(I,I)
         Z(I+1,I+1)=Z(I,I)
            DO K=1,N
            T(K,I+1)=T(K,I)
            U(K,I+1)=U(K,I)
            ENDDO
         ENDDO
            I=0
   97       A(I+1,I+1)=VR
            Z(I+1,I+1)=VI
            DO K=1,N
            T(K,I+1)=EN(K)
            U(K,I+1)=EN(K+N)
            ENDDO
      ENDDO

!     NORMALIZE VECTORS (BIGGEST COMPONENT TO UNITY)
      DO J=1,N
      ZM=0.D0
        DO I=1,N
        ZI=DABS(T(I,J))+DABS(U(I,J))
        IF(ZI.GE.ZM) THEN
        IM=I
        ZM=ZI
        ENDIF
        ENDDO
        ZM=T(IM,J)
        ZI=U(IM,J)
        DO I=1,N
        CALL DIVCD(T(I,J),U(I,J),ZM,ZI,TR,TI)
        T(I,J)=TR
        U(I,J)=TI
        ENDDO
      ENDDO
      RETURN
    END

!end of file tcomeig.f90
    
    