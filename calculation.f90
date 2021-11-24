!---------------------------------------------------------------!
!                                                               !
!   Calculate GAF(Genenralized aerodynamic forces) and others   !
!   corresponds to gnenralized_forces in EZASE                  !
!                                                               !
!---------------------------------------------------------------!
    
subroutine calculation

    use aerodynamic_force
    use DLM_var
    use input_var
    use panel
    use CLU
    use CLU_D
    
    implicit none
    
    !-------------------!
    !                   !
    !   DECLARATIONS    !
    !                   !
    !-------------------!
    
    double precision d_1, xll, xul, xlr, xur, yll, yul,ylr, yur
    integer t1,N_modeshape, tl, p, l,h, tol,dl,m,k
    integer next, w2,w3,w4, newdof, bl, errorcode
    integer experiment
    double precision, allocatable :: arrays(:,:), z_act_temp(:,:), S_z(:,:)
    double precision zll, zul, zlr, zur
    DOUBLE PRECISION BRDD(1,1), BIDD(1,1)
    double precision rll,rul,rlr,rur,ra, BRD, BID
    double precision yp
    double precision, allocatable :: Zxcp(:,:,:), Zf(:,:), Zcp(:,:,:), Z_cp(:,:,:)
    double complex, allocatable :: w5(:), Zf_t(:,:)
    double complex test1(2,2), test2(2,2), test3(2,2)
    double complex, allocatable :: DC2_temp(:,:)
    DOUBLE PRECISION QTG_RRm(1,1), qtg_IIM(1,1)
    DOUBLE PRECISION, ALLOCATABLE :: Z_ACT_T(:,:)
    double precision, allocatable :: Brm(:,:), Bim(:,:), Brm_T(:,:), Bim_T(:,:)
    double complex, allocatable :: An(:,:,:)
    double complex, allocatable :: Qtg_temp1(:,:)
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !---------------------------------------------------!
    !                                                   !   
    !   variables for use in example (DLM by  Blair)    !
    !                                                   !
    !---------------------------------------------------!
    
    double complex, allocatable :: Pressure_coeff(:,:,:), P_temp(:,:)
    double complex :: normalwash_ex(9,1)
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!
    !                   !
    !    ALLOCATIONS    !
    !                   !
    !!!!!!!!!!!!!!!!!!!!!
    
    experiment =1   !EXPERIMENT == 1, COMPARISON WITH DLM_BLAIR
    
    allocate(normalwash1(t, FEM_y)) ! Normalwash = upwash = downwash
    normalwash = (0,0);
    ALLOCATE(QTG(num_reduced_freq, FEM_Y,1))
    QTG = (0,0);
    ALLOCATE(QTG_temp(FEM_Y,num_reduced_freq))
    QTG_temp = (0,0);
    allocate( DC2_temp(t,t))    ! AIC matrix, DC2, where Drs is the element of Dc2, AIC matrix
    DC2_temp = (0,0);
    allocate(An(3+NMAX,FEM_Y,FEM_Y))
    An = (0,0);
    ALLOCATE(AN_(3+NMAX,FEM_Y,FEM_Y))
    An_ = (0,0);
    allocate(Pressure_coeff(num_reduced_freq, t,FEM_y), P_temp(t,FEM_y))    ! 1st : reduced_Freq, 2nd : panel No. 3rd : modeshape No,
    Pressure_coeff = (0,0);
    w3 = FEM_x/3
    
    allocate(AT(3+nmax,FEM_y),  Ata_temp(3+nmax,3+nmax), ATA_inv(3+nmax,3+nmax))
    AT = (0,0);ATA_temp = (0,0);ATA_inv = (0,0);
    allocate(ATA(3+nmax, 3+nmax))
    ATA = (0,0);
    allocate(arrays(w3,1), w5(t*t))
    arrays = (0,0);w5 = (0,0);
    allocate(z_act(num_y_panel+1, num_x_panel+1),z_act_T(num_X_panel+1, num_Y_panel+1), z_act_temp(num_y_panel+1, num_x_panel+1))
    z_act =0.d0;z_act_T =0.d0;z_act_temp =0.d0;
    allocate(z_c(num_X_panel+1, num_Y_panel+1),x_c(num_y_panel+1, num_x_panel+1),y_c(num_y_panel+1, num_x_panel+1),X_w(num_y_panel+1, num_x_panel+1),y_w(num_y_panel+1, num_x_panel+1),z_w(num_X_panel+1, num_Y_panel+1))
    z_c =0.d0;x_c =0.d0;y_c =0.d0;X_w =0.d0;y_w =0.d0;z_w =0.d0;
    allocate(S_z(FEM_x/3, FEM_y))
    S_z =0.d0;
    allocate(ZCp(num_reduced_freq,t,num_reduced_freq), Z_cp(num_reduced_freq,t, FEM_y), Zxcp(num_reduced_freq,t,num_reduced_freq), Zf(t, FEM_y))
    Zcp = 0.d0;z_cp = 0.d0;zxcp = 0.d0;zf = 0.d0;
    allocate(deflect(t), slope(t))  ! Deflection and slope of individual aerodynamic panel
    deflect = 0.d0;slope = 0.d0;
    allocate(normalwash(num_reduced_freq,t, FEM_y))
    normalwash = 0.d0;
    allocate(phaseangle(t), gwash(1,t)) ! used at gust analysis, phase angle and gust wash
    Phaseangle = 0.d0;gwash = 0.d0;
    allocate(Z_f(t, FEM_y), Zf_T(FEM_y, t))
    z_f = 0.d0; zf_t = 0.d0;
    allocate(Qt_temp(FEM_y, FEM_y), Qt(num_reduced_freq,FEM_y, FEM_y),DC2_inv(t,t), QOUT(num_reduced_freq,FEM_y, FEM_y))    ! Dc2 inverse used since Dc2 is destroyed when inversing matrix
    Qt_temp = 0.d0;qt = 0.d0;DC2_inv = 0.d0;qout = 0.d0;
    allocate(C(3+nmax,FEM_y), Cw2(1, FEM_y), ANG(3+nmax,FEM_y),Qsave(num_reduced_freq,FEM_y, FEM_y))
    C= 0.d0;Cw2 = 0.d0;ANG = 0.d0;QSAVE = 0.d0;
    ALLOCATE(QTG_RR(FEM_Y,num_reduced_freq), QTG_II(FEM_Y,num_reduced_freq))
    qtg_rr = 0.d0;qtg_ii = 0.d0;
    allocate(Qtg_r(FEM_Y,num_reduced_freq), Qtg_i(FEM_Y,num_reduced_freq))
    qtg_r = 0.d0;qtg_i = 0.d0;
    allocate(Qt_R(FEM_y, 1),Qt_Rr(1, FEM_y),Qt_i(FEM_y, 1),Qt_ii(1, FEM_y))
    qt_r = 0.d0;qt_rr = 0.d0;qt_i = 0.d0;qt_ii = 0.d0;
    Allocate(Qtg_temp1(t,FEM_Y))
    Qtg_temp = (0,0);
    Zf = (0,0);
    Zf_t = (0,0);
    
    
    newdof=0;
    
    do w2 = 1, FEM_x/3
        S_z(w2,:) = S(w2*3-2,:) ! S is the whole modeshape matrix, S_z is the modeshape matrix in z direction (deflection upward)
    end do
    
    t1 = 0
    
    do w2 = 2, num_X_panel+1
        do w3 = 1, num_y_panel+1
            
            t1 = t1+1
            
            IF (T1<=FEM_X/3) THEN
                Xv(t1) = x_3d(w3,w2)
                Yv(t1) = y_3d(w3,w2)
            END IF
            
        end do
    end do
    
    N_modeshape = FEM_y
    
    do t1 = 1, N_modeshape  ! For each modeshape
        
        print*, 'Calculating modeshapes for ', t1
        
        w3 = FEM_x/3 ! number of modeshape, z direction
        
        do w2 = 1, w3
            arrays(w2,1) = S(w2*3-2-newdof,t1)  ! save z_displacement to vector arrays.
        end do
        
        next = 0;
        
        z_act_temp = reshape(arrays, (/num_y_panel+1, num_x_panel/))    ! reshape arrays vector
        
        do w3 = 1, num_y_panel+1
            do w2 = 1, num_x_panel
                z_act(w3,w2+1) = z_act_temp(w3,w2)
            end do
        end do
        
        z_act(:, 1) = 0
        
        z_act_T = transpose(z_act)
        z_c = z_act_T
        x_c = x_3d
        y_c = y_3d
        x_w = x_3d
        y_w = y_3d
        z_w = z_act_T
 
        do tl = 1, t
            
            xll = y1(tl)
            xul = y4(tl)
            xlr = y2(tl)
            xur = y3(tl)
            yll = x1(tl)
            yul = x4(tl)
            ylr = x2(tl)
            yur = x3(tl)
            
            !!!!!!!!!!!!!!!!line241
            
            zll= 0;zul=0;zlr=0;zur=0
            
            do w2 = 1, FEM_x/3
                if (abs(Xv(w2)-xll)<0.00001 .and. abs(Yv(w2)-yll)<0.00001) then ! close enough!, numerical errors
                    zll = S_z(w2,t1)
                end if
            end do
            
            do w2 = 1, FEM_x/3
                if (abs(Xv(w2)-xul)<0.00001 .and. abs(Yv(w2)-yul)<0.00001) then
                    zul = S_z(w2,t1)
                end if
            end do
            
            do w2 = 1, FEM_x/3
                if (abs(Xv(w2)-xlr)<0.00001 .and. abs(Yv(w2)-ylr)<0.00001) then
                    zlr = S_z(w2,t1)
                end if
            end do
            
            do w2 = 1, FEM_x/3
                if (abs(Xv(w2)-xur)<0.00001 .and. abs(Yv(w2)-yur)<0.00001) then
                    zur = S_z(w2,t1)
                end if
            end do
            
            p=1
            !---------------------------------------------------!
            !                                                   !
            !   Following is for the aeropanel Cps inside       !
            !   calculate Cp according to each reudced freq.    !
            !                                                   !
            !---------------------------------------------------!
            
            do bl = 1, num_reduced_Freq
            
            tol = 1
            
            do w2 = 1, t
                if(Xcp(w2,1) >=xll .and. xcp(w2,1)<=xlr .and. Ycp(w2,1) <= yul .and. Ycp(w2,1)>=yll) then
                    h = w2;
                end if
            end do
            
            !---------------------------------------------------!
            !                                                   !
            !   Find all the aeropanel cps inside of structural !
            !   panel so that we can find cp for each aeropanel !
            !                                                   !
            !---------------------------------------------------!
            
            rll = sqrt(tol*(Xcp(h,1)-xll)**2+(YCP(h,1)-yll)**2)
            rul = sqrt(tol*(Xcp(h,1)-xul)**2+(YCP(h,1)-yul)**2)
            rlr = sqrt(tol*(Xcp(h,1)-xlr)**2+(YCP(h,1)-ylr)**2)
            rur = sqrt(tol*(Xcp(h,1)-xur)**2+(YCP(h,1)-yur)**2)
            ra = (1/rll+1/rul+1/rlr+1/rur)
            
            Zcp(bl,h,1) = (zll/rll+zlr/rlr+zul/rul+zur/rur)/ra
            Z_cp(bl,h,t1) = (zll/rll+zlr/rlr+zul/rul+zur/rur)/ra
            zxcp(bl,h,1) = -((zur+zul)/2-(zlr+zll)/2)/(yur-ylr)
            
            !---------------------------------------------!
            !                                             !
            !   Zcp : average deflection @ control point  !
            !   Zxcp : slope in the direction of the      !
            !   Leading edge (negative in sign convention !
            !                                             !
            !---------------------------------------------!
            
            do w2 = 1,t
                if(fx(w2) ==xll .and. fx(w2)<=xlr .and. fy(w2) <= yul .and. fy(w2)>=yll) then
                    h = w2;
                end if
            end do
            
            ! find distance of the structural node from the
            ! aero force point.
            
            rll = sqrt(tol*(fy(h)-xll)**2+(fx(h)-yll)**2)
            rul = sqrt(tol*(fy(h)-xul)**2+(fx(h)-yul)**2)
            rlr = sqrt(tol*(fY(h)-xlr)**2+(fx(h)-ylr)**2)
            rur = sqrt(tol*(fy(h)-xur)**2+(fx(h)-yur)**2)
            ra = (1/rll+1/rul+1/rlr+1/rur)
            
            Zf(h,t1) = (zll/rll+zlr/rlr+zul/rul+zur/rur)/ra
            if (t1 == 6)then
                !print*, Xv, Yv;pause;
                !print*, Zf(:, t1);pause;
            end if
            !   AVERAGE DEFLECTION AT THE FORCE POINT
            end do
        end do
        
        !---------------------------------------------------!
        !                                                   !
        !   Following calculate non-dimensional pressure    !
        !   coeff., is Cp/(dynamic pressure) for system at  !
        !   1/4 chord points from the leading edge of each  !
        !   aero panels.                                    !
        !                                                   !
        !---------------------------------------------------!
        
        dl = 1; ! For number of airspeed ; 1, not needed to be more than 1
        do m = 1, num_reduced_freq  ! number of reduced_frequency
        
        deflect = Zcp(m,:,1)
        slope = Zxcp(m,:,1)
        
        normalwash(m, :,t1) = slope(:) + cmplx(0,1)*reduced_freq(m)*2/length_y_total*deflect(:)
        ! Column vector of the downwash
        
        end do
    end do
    
    y_gust = 0
    do m = 1, num_reduced_freq
    
        do w2 = 1, t
            yp = CP(w2,1)
            phaseangle (w2) = -2*reduced_freq(m)/length_y_total*((length_y_total-yp)-y_gust)
            gwash(1,w2) = -exp(-cmplx(0,1)*2*reduced_freq(m)/length_y_total*((length_y_total-yp)-y_gust))
        end do
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! GWASH and PHASE ANGLE RIGHT   !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
        Z_f = Zf;
        DO W2 = 1, T
            DO W3 = 1, T
                Dc2(w2,w3)=Drs_rf(M,w2,w3) 
            END DO
        END DO
    
        dl = 1
        DC2_temp = DC2
        
        do w2 = 1, t
            do w3 = 1, FEM_y
                normalwash1(w2,w3) = normalwash(m, w2,w3)
            end do
        end do
    
        DC2_INV = 0.D0
       
        call CINVERT_LU(DC2_temp, t, DC2_inv,errorcode)
    
    ! FOR THE INVERSE TEST....
    !PRINT*, 'INVERSE TEST : ', MATMUL(DC2, DC2_INV);PAUSE;
        DC2_TEMP = MATMUL(DC2, DC2_INV)
    !PRINT*, DC2_INV;PAUSE;
    
        Zf_t = transpose(Zf)
    
        Qt_temp = -matmul(matmul(Zf_t, DC2_inv),normalwash1(:,:))*A_p(1)
        QT(M,:,:) = QT_TEMP(:,:)
        
        !!!!!!!!!!!!!!!!!!!!!!!Added Code
        
        P_temp=matmul(Dc2_inv, normalwash1)
        Pressure_coeff(M,:,:) = P_temp(:,:)
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        Zf_t = zf_t
        dc2_inv = dc2_inv
        gwash = gwash
        a_p = a_p
        
        !print*, FEM_Y,num_reduced_freq;print*, size(Zf_t);print*, size(DC2_inv);print*, size(gwash);print*, sizeof(QTG_temp);pause;
        
        Qtg_temp1= -matmul(Zf_t, DC2_inv)
        
        Qtg_temp = matmul(qtg_temp1,gwash)*A_p(1)
        QTG(M,:,:) = QTG_TEMP(:,:)
    
    end do
    
    !-----------------------------------------------!
    !                                               !
    !   I believe Qt is nondimensional Cp, which    !
    !    P = Cp*(dynamic pressure), therefore,      !
    ! Pressure(per each panel) = Qt*dynamic pressure!
    !                                               !
    !-----------------------------------------------!
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                                      !
    !           Br and Bi verified.        !
    !                                      !   
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    Qsave = Qt;
    
    AT = 0.d0
    
    ATA = 0.d0
    AT = 0.d0
    
    do k = 1, FEM_y
        do dl = 1, 1
            
            ATA = 0.d0
            AT = 0.d0
            
            do m = 1,NUM_REDUCED_FREQ
                do w3 = 1, FEM_y
                Qt_r(w3,1) =  real(Qt(M,w3,k))
                Qt_rr = transpose(Qt_r)
                Qt_i(w3,1) = imag(Qt(M,w3,k))
                Qt_ii = transpose(Qt_i)
                end do
                
                !print*, Qt;pause;
                
                allocate(Brm(1, 3+nmax), Bim(1, 3+nmax), Brm_T(3+nmax,1), Bim_T(3+nmax,1))
                
                Brm(1, :) = Br(m,:)
                Bim(1, :) = Bi(m,:)
                
                !print*, 'nmax : ', Br
                !print*, Br(1,1),Br(1,2),Br(1,3),Br(1,4),Br(1,5),Br(2,1),Br(2,2),Br(2,3),Br(2,4),Br(2,5);pause;
                
                call transpose_vector(Brm, 3+nmax, Brm_T)
                call transpose_vector(Bim, 3+nmax, Bim_T)
                
                AT = AT+MATMUL(BRM_t, QT_RR)+MATMUL(BIM_T, QT_II)
                
                !print*, MATMUL(BIM_T, QT_II);pause;
                
                AtA = (AtA+matmul(Brm_T,BRM)+matmul(Bim_T,Bim))
                
                do w2 = 1, 3+nmax
                    do w3 = 1, 3+nmax
                        BRDD = matmul(transpose(Br),Br)
                        BIDD = matmul(transpose(Bi(:,:)), Bi(:,:))
                        !AtA(w2,w3) = AtA(w2,w3) + BRDD(1,1) + BIDD(1,1)
                    end do
                end do
                
                deallocate(Brm, Bim,Brm_t,Bim_t)
                
                ATA_temp = ATa
            !PRINT*, ATA_TEMP;PAUSE;
                ATA_INV = 0.D0
                ATA = ATA_TEMP
            
                call CINVERT_LU_D(ATA_temp, nmax+3, ATA_inv, errorcode)
            
            !PRINT*, ATA_INV;PAUSE;
            !PRINT*, ATA_TEMP;PAUSE;
            
                if (errorcode == 0) then
                    C = matmul(ATA_inv, AT)
                
                !ORIGINAL :C = matmul(ATA_inv*10, AT)
                !WHY? MAGIC NUMBER??
                else
                    print*, ERRORCODE,'something wrong with CINVERT_LU at calculating C SEE CALCULATION.F90 @ LINENO.412';pause;
                end if
            
                do w2 = 1, 3+nmax
                    An(w2, :, k) = C(w2,:)
                end do
                
                An_temp(m,:,:,:) = An(:,:,:)
                
            end do
            
            
        end do
    end do
    
    !print*, 'C:', At;pause;
    !print*, 'ata:', Ata;pause;
    
    !---------------------------!
    !   VERIFICATION MODULE     !
    !---------------------------!
    
    !print*, 'Xv : ', Xv
    !print*, 'Yv : ', Yv
    !print*, 'xll : ', xll
!    print*, 'yll : ', xul
!    print*, 'yll : ', xlr
!    print*, 'yll : ', xur
!    print*, 'yll : ', yll
!    print*, 'yll : ', yu
!print*, 'yll : ', ylr
!print*, 'yll : ', yur;pause;
!    PRINT*, 'H : ', H;PAUSE;
!    PRINT*, ZLL;PAUSE;
!    PRINT*, ZLR;PAUSE;
    !print*, 'holo', ZCP(10,:,1);pause;
!print*, 'deflect', deflect;pause;
!print*, 'slope', slope;pause;
!    print*, 'normalwash', normalwash;pause;
    !PRINT*, 'Zf : ', s_Z(1,1),s_Z(1,2),s_Z(1,3),s_Z(1,4),s_Z(1,5),s_Z(1,6),s_Z(1,7),s_Z(1,8),s_Z(1,9),s_Z(1,10);PAUSE;
    !print*, normalwash1;pause;
    !PRINT*, 'Qt : ', QT_RR;PAUSE;
    !PRINT*, 'Zf : ', Zf_T;PAUSE;
    !print*, 'AT : ', AT;pause;
    !print*, 'ATa : ', ATa;pause;
    !PRINT*, 'ZCP : ', ZCP;PAUSE;
    !PRINT*, 'ZXCP : ', ZXCP;PAUSE;
    
    !---------------------------------------------------------------------------!
    !                                                                           !
    !                            *STATUS REPORT*                                !
    !                                                                           !
    !   1. DC2, DC2_INV ARE CORRECT AS OF 2017-10-08                            !
    !   2. Zcp, Zxcp doesn't have all t-elements (t= num_x_panel*num_y_panel)   !
    !   3. Something wrong with zll, zul, zlr, zlr                              !
    !   4. xll, xul, xlr, xur and y series match                                !
    !   5. Xv, Yv match as well. then what's the problem?                       !
    !   6. S_z checks out as well                                               !
    !   7. SOLVED!! - 2017-10-09 Trunication error caused, solved by using abs  !
    !      Function and 0.00000001...                                           !
    !   8. ZCP correct, At incorrect, Qt incorrect                              !
    !   9. Qt corrected, double complex problem in normalwash1                  !
    !   10. Fixed with transpose_vector                                         !
    !   11. Qt wrong again........                                              !
    !   12. Solved Qt, ZF correct. what next?                                   !
    !   13. possible Qt error                                                   !
    !   14. Qt is finalized as Cp                                               !
    !   15. Normalwash is the specified upwash                                  !
    !   16. Qt = Dc_inv*upwash*A_p, Pressure Coeff. per box.                    !
    !   17. DONE!!                                                              !
    !                                                                           !
    !---------------------------------------------------------------------------!
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                                       !
    !       GUST MODULE, NOT DONE YET       !
    !                                       !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
!    DO K = 1,1 !number of reduced_freq
!        DO DL = 1,1
!            DO M = 1,NUM_REDUCED_FREQ
!                
!                allocate(Brm(1, 3+nmax), Bim(1, 3+nmax), Brm_T(3+nmax,1), Bim_T(3+nmax,1))
!                
!                Brm(1, :) = Br(m,:)
!                Bim(1, :) = Bi(m,:)
!                
!                call transpose_vector(Brm, 3+nmax, Brm_T)
!                call transpose_vector(Bim, 3+nmax, Bim_T)
!                
!                QtG_r(:,k) =  real(QtG(M,:,k))
!                QtG_rr = transpose(QtG_r)
!                QtG_i(:,k) = imag(QtG(M,:,k))
!                QtG_ii = transpose(QtG_i)
!                
!                QTG_RRM = MATMUL(BRM_T, QTG_RR)
!                
!                QTG_IIM = MATMUL(BIM_T, QTG_II)
!                
!                AT = AT+QTG_RRM+QTG_IIM
!                ATA = ATA+MATMUL(BRM_T,BRM)+MATMUL(BIM_T,BIM)
!                
!                deallocate(Brm, Bim,Brm_t,Bim_t)
!                
!            END DO
!            
!            CALL CINVERT_LU(ATA, 5,ATA_INV, ERRORCODE)
!            
!            IF (ERRORCODE == 1) THEN
!            PRINT*, 'ERROR WITH ATA INVERSE'
!            PAUSE;
!            END IF
!            
!            C = MATMUL(ATA_INV, AT)
!            
!            DO W2 = 1, 3+NMAX
!                DO W4 = 1, FEM_Y
!                    ANG(W2,W4) = C(W2,W4)
!                END DO
!            END DO
!        END DO
!    END DO
    
    AN_ = AN
    
    QOUT = QT;
    
    !QOUT FORM => (M,:,:)
    
    !print*, 'Qout : ', Qout(1,:,1);pause;
    
    deallocate(fx,fy,An)
    deallocate(x_3d, y_3d)
    
end subroutine calculation
