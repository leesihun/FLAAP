!***********************************************************!
!                                                           !
!   NOT USED IN CURRENT PROGRAM, USE MATLAB                 !
!   .m FILE TO READ .F06 AND PRODUCE 'MODESHAPE_input.txt'  !
!                                                           !
!***********************************************************!
    
    subroutine read_F06
    
    use input_var
    
    character :: tmp_cha
    integer :: number_mode
    double precision :: tmp_double
    integer :: j
    double precision, allocatable :: mode_S(:,:)
    integer :: i, eigen_length, blank
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    eigen_length=45
    number_mode = 10
    blank = 7
    FEM_y = 10
    number_mode = 10
    FEM_x = 45
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    allocate (mode_S(10,45))
    
    open(unit = 180,file ='20.dat', status = 'old')
    
    do i=1, 470
        read(180,*)
    enddo
    do j = 1, 10
        do i = 1, blank
            read(180,*)
        end do
        
        do i=1, eigen_length
            read(180,*), tmp_double, tmp_cha, tmp_double, tmp_double, Mode_S(j, i), tmp_double,tmp_double,tmp_double
        enddo
    end do
    
    close(180)
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                                                       !
    !            MODE_S IS THE MODESHAPE, Z DEFLECT         !
    !                                                       !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
end subroutine read_F06