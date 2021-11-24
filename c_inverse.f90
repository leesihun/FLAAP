!---------------------------------------------------------------!
!                                                               !
!   Inversion of complex matrix, C_inverse.f90                  !
!   used to inverse AIC matrix                                  !
!                                                               !
!---------------------------------------------------------------!
    
    
!*******************************************************
!*   LU decomposition routines used by test_clu.f90    *
!*                                                     *
!*                 F90 version by J-P Moreau, Paris    *
!*                (with extension to complex domain)   *
!*                        (www.jpmoreau.fr)            *
!* --------------------------------------------------- *
!* Reference:                                          *
!*                                                     *
!* "Numerical Recipes by W.H. Press, B. P. Flannery,   *
!*  S.A. Teukolsky and W.T. Vetterling, Cambridge      *
!*  University Press, 1986".                           *
!*                                                     * 
!*******************************************************
MODULE CLU

CONTAINS

!*****************************************************************
!* Given an N x N complex matrix A, this routine replaces it by  *
!* the LU decomposition of a rowwise permutation of itself. A    *
!* and N are input. INDX is an output vector which records the   *
!* row permutation effected by the partial pivoting; D is output *
!* as -1 or 1,  depending on whether the number of row inter-    *
!* changes was even or odd, respectively. This routine is used in*
!* combination with LUBKSB to solve linear equations or to invert*
!* a matrix. Return code is 1, if matrix is singular.            *
!*****************************************************************
 Subroutine CLUDCMP(A,N,INDX,D,CODE)
 REAL*8, PARAMETER :: TINY=2.2D-16
 double COMPLEX  AMAX,DUM, SUM, A(N,N)
 INTEGER CODE, D, INDX(N)
 double COMPLEX, pointer ::  VV(:)     !complex vector (n)

 allocate(VV(N),stat=ialloc)
 
 D=1; CODE=0

 DO I=1,N
   AMAX=CMPLX(0.,0.)
   DO J=1,N
     IF (ABS(A(I,J))>ABS(AMAX)) then
         AMAX=A(I,J)
     end if
   END DO ! j loop
   IF(ABS(AMAX)<TINY) THEN
     CODE = 1
     print*, 'ERROR !', AMAX
     RETURN
   END IF
   VV(I) = 1. / AMAX
 END DO ! i loop

 DO J=1,N
   DO I=1,J-1
     SUM = A(I,J)
     DO K=1,I-1
       SUM = SUM - A(I,K)*A(K,J) 
     END DO ! k loop
     A(I,J) = SUM
   END DO ! i loop
   AMAX = CMPLX(0.,0.)
   DO I=J,N
     SUM = A(I,J)
     DO K=1,J-1
       SUM = SUM - A(I,K)*A(K,J) 
     END DO ! k loop
     A(I,J) = SUM
     DUM = VV(I)*ABS(SUM)
     IF(CDABS(DUM).GE.CDABS(AMAX)) THEN
       IMAX = I
       AMAX = DUM
     END IF
   END DO ! i loop  
   
   IF(J.NE.IMAX) THEN
     DO K=1,N
       DUM = A(IMAX,K)
       A(IMAX,K) = A(J,K)
       A(J,K) = DUM
     END DO ! k loop
     D = -D
     VV(IMAX) = VV(J)
   END IF

   INDX(J) = IMAX
   IF(ABS(A(J,J)) < TINY)  A(J,J) = CMPLX(TINY,0.)

   IF(J.NE.N) THEN
     DUM = 1. / A(J,J)
     DO I=J+1,N
       A(I,J) = A(I,J)*DUM
     END DO ! i loop
   END IF 
 END DO ! j loop

 RETURN
 END subroutine CLUDCMP


!********************************************************************
!* Solves the set of N complex linear equations A . X = B.  Here A  *
!* is input, not as the matrix A but rather as its LU decomposition,*
!* determined by the routine CLUDCMP. INDX is input as the permuta- *
!* tion vector returned by CLUDCMP. B is input as the right-hand    *
!* side complex vector B, and returns with the solution vector X. A,*
!* N and INDX are not modified by this routine and can be used for  *
!* successive calls with different right-hand sides. This routine is*
!* also efficient for plain complex matrix inversion.               *
!********************************************************************
 Subroutine CLUBKSB(A,N,INDX,B)
 double COMPLEX  SUM, A(N,N),B(N)
 INTEGER INDX(N)

 II = 0
 
 DO I=1,N
   LL = INDX(I)
   SUM = B(LL)
   B(LL) = B(I)
   
   IF(II.NE.0) THEN
     DO J=II,I-1
       SUM = SUM - A(I,J)*B(J)
       
     END DO ! j loop
   ELSE IF(CDABS(SUM).NE.0.) THEN
     II = I
   END IF
   B(I) = SUM
 END DO ! i loop

 DO I=N,1,-1
   SUM = B(I)
   
   IF(I < N) THEN
     DO J=I+1,N
       SUM = SUM - A(I,J)*B(J)
     END DO ! j loop
   END IF
   B(I) = SUM / A(I,I)
 END DO ! i loop
 
 !print*, SUM;pause;

 RETURN
 END subroutine CLUBKSB

!************************************************************
!* Inversion of a complex square matrix by LU decomposition *
!* -------------------------------------------------------- *
!* INPUTS:                                                  *
!*          A: complex square matrix of size n x n          *
!*          n: size of matrix A                             *
!* OUTPUTS:                                                 *
!*          Y; Inverse of A of size n x n                   *
!*         rc: error code (must be zero)                    *
!************************************************************
Subroutine CINVERT_LU(A,n,Y,rc)

  double COMPLEX A(n,n), Y(n,n)

  integer rc, D, i,j

  integer,pointer ::  INDX(:)   !integer vector (n)

  do i = 1, n
              Y(i,i) =(1,0)
              !print*, Y(i,i)
  end do
  
  !dynamic allocations
  allocate(INDX(n),stat=ialloc)

!call LU decomposition routine
  call CLUDCMP(A,n,INDX,D,rc)

!call solver if previous return code is ok
!to obtain inverse of A one column at a time
  if (rc.eq.0) then
    do j=1, n
      call CLUBKSB(A,n,INDX,Y(1,j))
    end do
  end if
!the inverse matrix is now in matrix Y
!the original matrix A is destroyed

  return

END Subroutine CInvert_lu

    END MODULE CLU

! end of file clu.f90
    
    
    
    SUBROUTINE CMATMUL(A,B,C,N,M)                                              
!*******************************************                                     
!*     MULTIPLY TWO COMPLEX MATRICES       *
!* --------------------------------------- *                                     
!* INPUTS:    A  MATRIX N*N                *                                     
!*            B  MATRIX N*M                *                                     
!*            N  INTEGER                   *                                     
!*            M  INTEGER                   *                                     
!* --------------------------------------- *                                     
!* OUTPUT:    C  MATRIX N*M, PRODUCT A*B   *                                     
!*                                         *                                     
!******************************************* 
  double COMPLEX A(N,N),B(N,M),C(N,M),SUM                                           
  DO I=1,N                                                                  
    DO J=1,M                                                                
      SUM=0.                                                                
      DO K=1,N                                                              
        SUM=SUM+A(I,K)*B(K,J)                                               
      ENDDO                                                                 
      C(I,J)=SUM                                                            
    ENDDO                                                                   
  ENDDO                                                                     
  RETURN                                                                    
    END