!***************************************************
!* Sorting an array with straight insertion method *
!* ----------------------------------------------- *
!* REFERENCE:                                      *
!* "NUMERICAL RECIPES By W.H. Press, B.P. Flannery,*
!*  S.A. Teukolsky and W.T. Vetterling, Cambridge  *
!*  University Press, 1986" [BIBLI 08].            *
!* ----------------------------------------------- *
!* SAMPLE RUN:                                     *
!*                                                 *
!* Table to be sorted:                             *
!*                                                 *
!*  113.4  930.5   73.7  590.9   94.0              *
!*  663.0   82.2   85.7  809.6  615.4              *
!*  448.9  964.8  608.6  297.7  433.1              *
!*  369.6   32.9   24.0  891.4  929.7              *
!*  888.6  295.2  106.8  965.7  791.5              *
!*                                                 *
!* Sorted table (straight insertion):              *
!*                                                 *
!*   24.0   32.9   73.7   82.2   85.7              *
!*   94.0  106.8  113.4  295.2  297.7              *
!*  369.6  433.1  448.9  590.9  608.6              *
!*  615.4  663.0  791.5  809.6  888.6              *
!*  891.4  929.7  930.5  964.8  965.7              *
!*                                                 *
!***************************************************
subroutine C_matrix_sort_imag(A, N)

integer N
double complex    A(N)     !Table to be sorted
double complex    MAX_VALUE  !Maximum value of table

  !initialize random number generator
  
  max_value = cmplx(10**6, 10**6) 

  !N=25  !initialize size of table
  !MAX_VALUE = 1000

  !call TWRIT(N,A)
 
  !call sorting subroutine
  call PIKSRT_imag(N,A)

  !print *,' '
  !print *,'Sorted table (straight insertion):'
  !call TWRIT(N,A)

END subroutine C_matrix_sort_imag


!*****************************************************
!* Sorts an array ARR of length N in ascending order *
!* by straight insertion.                            *
!* ------------------------------------------------- *
!* INPUTS:                                           *
!*	    N	  size of table ARR                  *
!*          ARR	  table to be sorted                 *
!* OUTPUT:                                           *
!*	    ARR   table sorted in ascending order    *
!*                                                   *
!* NOTE: Straight insertion is a N©÷ routine and      *
!*       should only be used for relatively small    *
!*       arrays (N<100).                             *
!*****************************************************         
SUBROUTINE PIKSRT_imag(N,ARR)
implicit none
  double complex ARR(N)
  integer :: N, j, i
  double complex :: a
  do j=2, N
    a=ARR(j)
    do i=j-1,1,-1
      if (imag(ARR(i))<=imag(a)) goto 10
      ARR(i+1)=ARR(i)
    end do
	i=0
10  ARR(i+1)=a
  end do
  return
END

!write table of size N to standard output
SUBROUTINE TWRIT_imag(N,ARR)
double complex :: ARR(N)
  print *,' '
  WRITE(*,10) (ARR(I),I=1,N)
  return
10 FORMAT(5(' ',F6.1))
    END

!end of file sort1.f90
    
    
    
    
    
    
    
    
    !***************************************************
!* Sorting an array with straight insertion method *
!* ----------------------------------------------- *
!* REFERENCE:                                      *
!* "NUMERICAL RECIPES By W.H. Press, B.P. Flannery,*
!*  S.A. Teukolsky and W.T. Vetterling, Cambridge  *
!*  University Press, 1986" [BIBLI 08].            *
!* ----------------------------------------------- *
!* SAMPLE RUN:                                     *
!*                                                 *
!* Table to be sorted:                             *
!*                                                 *
!*  113.4  930.5   73.7  590.9   94.0              *
!*  663.0   82.2   85.7  809.6  615.4              *
!*  448.9  964.8  608.6  297.7  433.1              *
!*  369.6   32.9   24.0  891.4  929.7              *
!*  888.6  295.2  106.8  965.7  791.5              *
!*                                                 *
!* Sorted table (straight insertion):              *
!*                                                 *
!*   24.0   32.9   73.7   82.2   85.7              *
!*   94.0  106.8  113.4  295.2  297.7              *
!*  369.6  433.1  448.9  590.9  608.6              *
!*  615.4  663.0  791.5  809.6  888.6              *
!*  891.4  929.7  930.5  964.8  965.7              *
!*                                                 *
!***************************************************
subroutine C_matrix_sort_real(A, N)

integer N
double complex    A(N)     !Table to be sorted
double complex    MAX_VALUE  !Maximum value of table

  !initialize random number generator
  
  max_value = cmplx(10**6, 10**6) 

  !N=25  !initialize size of table
  !MAX_VALUE = 1000

  !call TWRIT(N,A)
 
  !call sorting subroutine
  call PIKSRT_real(N,A)

  !print *,' '
  !print *,'Sorted table (straight insertion):'
  !call TWRIT(N,A)

END subroutine C_matrix_sort_real


!*****************************************************
!* Sorts an array ARR of length N in ascending order *
!* by straight insertion.                            *
!* ------------------------------------------------- *
!* INPUTS:                                           *
!*	    N	  size of table ARR                  *
!*          ARR	  table to be sorted                 *
!* OUTPUT:                                           *
!*	    ARR   table sorted in ascending order    *
!*                                                   *
!* NOTE: Straight insertion is a N©÷ routine and      *
!*       should only be used for relatively small    *
!*       arrays (N<100).                             *
!*****************************************************         
SUBROUTINE PIKSRT_real(N,ARR)
implicit none
  double complex ARR(N)
  integer :: N, j, i
  double complex :: a
  do j=2, N
    a=ARR(j)
    do i=j-1,1,-1
      if (real(ARR(i))<=real(a)) goto 10
      ARR(i+1)=ARR(i)
    end do
	i=0
10  ARR(i+1)=a
  end do
  return
END

!write table of size N to standard output
SUBROUTINE TWRIT_real(N,ARR)
double complex :: ARR(N)
  print *,' '
  WRITE(*,10) (ARR(I),I=1,N)
  return
10 FORMAT(5(' ',F6.1))
    END

!end of file sort1.f90