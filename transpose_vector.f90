subroutine transpose_Vector(A, N, B)
    
    implicit none
    
    integer N, holo
    double precision :: A(1,N), B(N,1)
    
    do holo = 1,N
        B(holo,1) = A(1, holo)
    end do
    
    end subroutine transpose_vector
    