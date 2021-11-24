
    
    !***************************************************!
    !                                                   !
    !   SAVE DAMPING AND FREQUENCY ACCORDING TO VEL     !
    !   OCITY, 'VG_SAVE.DAT', 'VF_SAVE.DAT'             !
    !                                                   !
    !***************************************************!
    
    subroutine save_vg
    
    use statespace_var
    
    print*, 'saving V-g graph data to Vg_save.dat';
    
    print*, ' ';
    
    open(22, FILE='Vg_save.dat', status = 'new');
    
    write (22, *) 'For each mode, velocity, damping'
    
    do i = 1, 10 ! for each modes
        write (22, *) 'Mode number : ', i
        do j = 1, r_freq_count
            write (22,*) store_v(i, j), store_d(i, j)
        end do
    end do
    
    close (22)
    
    open (33, FILE='Vf_save.dat', status = 'new');
    
    write (33, *) 'For each mode, velocity, frequency'
    
    do i = 1, 10 ! for each modes
        write (33, *) 'Mode number : ', i
        do j = 1, r_freq_count
            write (33,*) store_v(i, j), store_f(i, j)
        end do
    end do
    
    close (33)
    
    end subroutine save_vg