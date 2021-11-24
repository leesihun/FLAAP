!*************************************************************
!*         Determinant of a complex square matrix            *
!*          By Gauss Method with full pivoting               *
!* --------------------------------------------------------- *
!* SAMPLE RUN:                                               *
!* Calculate the determinant of complex matrix:              *
!*  ( 47,-15) ( 62,5) (  0,-72) (61, 20)                     *
!*  (  6, 14) (-17,3) (-102,91) ( 7,-12)                     *
!*  ( 13, 55) ( 32,8) (  41, 7) (25,  1)                     *
!*  (111,25)  ( 40,0) ( 12,-82) (58,-30)                     *
!*                                                           *
!*  Det =  (1.7416564E+07,-1.0598320E+07)                    *
!*                                                           *
!* --------------------------------------------------------- *
!* Ref.: "Algèbre, Algorithmes et programmes en Pascal       *
!*        By Jean-Louis Jardrin, DUNOD Paris, 1988".         *
!*                                                           *
!*                       F90 Release By J-P Moreau, Paris.   *
!*                              (www.jpmoreau.fr)            *
!************************************************************* 
!Program CDetMat

!integer, parameter :: NMAX = 20

!  Complex A(1:NMAX,1:NMAX), det

! Example #1
! N = 3   (size of complex matrix A)

! A(1,1) = CMPLX(1.0,0.0)
! A(1,2) = CMPLX(0.0,1.0)
! A(1,3) = CMPLX(0.0,1.0)

! A(2,1) = CMPLX(0.0,1.0)
! A(2,2) = CMPLX(1.0,0.0)
! A(2,3) = CMPLX(1.0,0.0)

! A(3,1) = CMPLX( 0.0,1.0)
! A(3,2) = CMPLX( 1.0,0.0)
! A(3,3) = CMPLX(-1.0,0.0)

! ( Det = -4.0 + 0 I )

! Example #2
 ! N = 4

  !A(1,1) = CMPLX(47.0,-15.0)
  !A(1,2) = CMPLX(62.0,  5.0)
  !A(1,3) = CMPLX( 0.0,-72.0)
  !A(1,4) = CMPLX(61.0, 20.0)

  !A(2,1) = CMPLX(   6.0, 14.0)
  !A(2,2) = CMPLX( -17.0,  3.0)
  !A(2,3) = CMPLX(-102.0, 91.0)
  !A(2,4) = CMPLX(   7.0,-12.0)

  !A(3,1) = CMPLX(13.0,-55.0)
  !A(3,2) = CMPLX(32.0,  8.0)
  !A(3,3) = CMPLX(41.0,  7.0)
  !A(3,4) = CMPLX(25.0,  1.0)

  !A(4,1) = CMPLX(111.0, 25.0)
  !A(4,2) = CMPLX( 40.0,  0.0)
  !A(4,3) = CMPLX( 12.0,-82.0)
  !A(4,4) = CMPLX( 58.0,-30.0)

  !eps=1.E-10

  !call DCGT(eps, N, A, det)

  !print *,' '
  !write(*,10,advance='no'); print *, det
  !print *,' '

  !stop

!10 format(' Det = ')

!END

 !*****************************************************************
 !* TSCGT procedure implements the triangularization algorithm of *
 !* Gauss with full pivoting at each step for a complex matrix, A *
 !* and saves the made transformations in KP and LP.              *
 !* ------------------------------------------------------------- *
 !* INPUTS:                                                       *
 !*          N:   size of complex matrix A                        *
 !*          A:   complex matrix of size N x N                    *
 !* OUTPUTS;                                                      *
 !*          it:  =0 if A is singular, else =1.                   *
 !*           C:  contains the upper triangular matrix and the    *
 !*               multipliers used during the process.            *
 !*          KP:  contains the column exchanges.                  *
 !*          LP:  contains the line exchanges.                    *
 !*****************************************************************
  Subroutine TSCGT(eps, N, A, it, C, KP, LP)
    integer, parameter :: NMAX = 10
    double Complex A(NMAX,NMAX), C(NMAX,NMAX)
	Integer KP(NMAX), LP(NMAX)
    double Complex C0,C1,P0,T0
    double precision :: eps

    C=A
    it=1; K=1
    do while (it==1.and.k<N)
      P0=C(k,k); l0=k; k0=k
      do i=k, N
        do j=1, N
          if (ABS(C(i,j)) > ABS(P0)) then
            P0=C(i,j)
            l0=i; k0=j
          end if
        end do
      end do 
      LP(k)=l0; KP(k)=k0
      if (ABS(P0) < eps) then
        it=0
      else
        if (l0.ne.k) then
          do j=k, N
            T0=C(k,j)
            C(k,j)=C(l0,j)
            C(l0,j)=T0
          end do
        end if
        if (k0.ne.k) then
          do i=1, N
            T0=C(i,k)
            C(i,k)=C(i,k0)
            C(i,k0)=T0
          end do
        end if
        do i=k+1, N
          C0=C(i,k)
		  C(i,k) = C0/P0
          do j=k+1, N
            C0=C(i,j)
			C1 = C(i,k)*C(k,j)
			C(i,j) = C0 - C1
          end do
        end do
        k=k+1
      end if
    end do
    if (it==1.and.ABS(C(N,N)) < eps)  it=0
	return
  End  !TSCGT

 !****************************************************************
 !* The DCGT procedure calculates the complex determinant of a   *
 !* complex square matrix by the Gauss method with full pivoting *
 !* ------------------------------------------------------------ *
 !* INPUTS:                                                      *
 !*         eps:  required precision                             *
 !*          N :  size of matrix A                               *
 !*          A :  complex matrix of size N x N                   *
 !* OUTPUT:                                                      *
 !*         det:  complex determinant.                           *                  
 !****************************************************************
  Subroutine DCGT(eps, N, A, det)
    integer, parameter :: NMAX = 10
    double Complex :: A(10,10), det
    double Complex :: C0,Z0,Z1
    double Complex :: C(NMAX,NMAX)
    double precision :: eps
    integer :: KP(NMAX), LP(NMAX)

    Z0 = 0.; Z1 = 1.
    
	call TSCGT(eps,N,A,it,C,KP,LP)
    
	if (it==0) then
      det = Z0
    else
      det=Z1
      do k=1, N
        C0=det 
		det = C0 * C(k,k)
      end do
      l=0
      do K=1, N-1
        if (LP(k).ne.k) l=l+1
        if (KP(k).ne.k) l=l+1
      end do
      if (Mod(l,2).ne.0)  det = -det
    end if
    !print*, 'A : ', A;pause;
    !print*, det;pause;
	return
  end  !DCGT

!end of file cdetmat.f90