subroutine aero_paneling_vlm
    
    !---------------------------------------------------------------!
    !                                                               !
    !   Subroutine aero_paneling_vlm.f90, perform aeropaneling for  !
    !   Vortex Lattice Method                                       !
    !   Basically same as aeropaneling.f90                          !
    !                                                               !
    !---------------------------------------------------------------!
    
    use panel
    use DLM_var
    use input_var
    use aerodynamic_force
    
    implicit none
    
    double precision x_step, temp_coord, temp_coord_d
    integer rx, ry, i, j, bl, g, k, l,w2, token
    
    allocate(a_p(num_y_panel* num_x_panel))
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      Area of each panels
    
    allocate(Cpy(aero_y_plates*num_x_panel), Cpx(aero_y_plates*num_x_panel))
    cpx = 0.d0; cpy = 0.d0;
    allocate(Y1(aero_y_plates*num_x_panel),Y2(aero_y_plates*num_X_panel),Y3(aero_y_plates*num_X_panel),Y4(aero_y_plates*num_X_panel))
    Y1 = 0.d0;Y2 = 0.d0;Y3 = 0.d0;Y4 = 0.d0
    allocate(x1(aero_y_plates*num_X_panel),x2(aero_y_plates*num_x_panel),X3(aero_y_plates*num_X_panel),x4(aero_y_plates*num_x_panel))
    X1 = 0.d0;x2 = 0.d0;x3 = 0.d0;x4 = 0.d0
    allocate(XCP(aero_y_plates*num_x_panel,2),yCP(aero_y_plates*num_x_panel,2))
    xCP = 0.d0;yCP = 0.d0;
    allocate(CP(aero_y_plates*num_x_panel,2))
    CP = 0.d0;
    allocate(lvpx(aero_y_plates*num_x_panel),lvpy(aero_y_plates*num_x_panel),pln1(aero_y_plates*num_x_panel,2),pln2(aero_y_plates*num_x_panel,2),plf(aero_y_plates*num_x_panel,2))
    lvpx = 0.d0;lvpy = 0.d0;pln1 = 0.d0;pln2 = 0.d0; plf = 0.d0;
    allocate(rvpx(aero_y_plates*num_x_panel),rvpy(aero_y_plates*num_x_panel),fx(aero_y_plates*num_x_panel),fy(aero_y_plates*num_x_panel),nc(aero_y_plates*num_x_panel,2))
    rvpx = 0.d0;Rvpy = 0.d0;fx = 0.d0;fy = 0.d0; nc = 0.d0;
    allocate(N1(aero_y_plates*num_x_panel,2), N2(aero_y_plates*num_x_panel,2))
    allocate(Xv(FEM_x/3),Yv(FEM_x/3))
    
    !xible = length_x_total
    !xibte = 0
    !xoble = length_x_total
    !xobte = 0
    
    dxib = abs(yible-yibte)/aero_y_plates   !!!!!!!!!!!!Changed, 2019-05-24 used to be abs(yible-yoble)/aero_y_plates
    dxob = abs(yoble-yobte)/aero_y_plates   !!!!!!!!!!!!Changed, 2019-05-24 used to be abs(yoble-yibte)/aero_y_plates
    !DYIB = ABS(XIBLE-XOBLE)/NUM_X_PANEL
    !DYOB = ABS(XIBTE-XOBTE)/NUM_X_PANEL
    
    !yoble = length_x_total
    !yible = 0
    
    dy = (length_x_total)/num_x_panel
    
    t =0
    print*, ' '
    print*, 'dxib: ', dxib
    print*, 'dxob: ', dxob
    print*, 'dy: ', dy
    
    if(dxib == dxob) then
        token =1
    else
        token = 0
    end if
    
            PRINT*, 'STILL SOMETHING WRONG WITH X1,X2,X3,X4, MAY NEED ADJUSTING';
            
            
    do rx =1, aero_y_plates
        do ry =1, num_x_panel
            t=t+1
            ix = rx
            iy = ry/aero_y_plates
            y1(t) = dy*(ry-1)
            y2(t) = dy*(ry)
            y3(t) = dy*ry
            y4(t) = dy*(ry-1)
            
            
            x1(t) = yibte-((yibte-yobte)/aero_y_plates)*(ry-1)+&
                ((yibte-yobte)/aero_y_plates)*(ry-1)*(rx-1)&
                +(-ABS(dxib-dxob)*(ry-1)/aero_y_plates+dxib)*(rX-1)
            x2(t) = yibte-((yibte-yobte)/aero_y_plates)*(ry)+&
                ((yibte-yobte)/aero_y_plates)*(ry)*(rx-1)&
                +(-ABS(dxib-dxob)*(ry)/aero_y_plates+dxib)*(rX-1)
            x3(t) = yibte-((yibte-yobte)/aero_y_plates)*(ry)+&
                ((yibte-yobte)/aero_y_plates)*(ry)*(rx)&
                +(-ABS(dxib-dxob)*(ry)/aero_y_plates+dxib)*(rX)
            x4(t) = yibte-((yibte-yobte)/aero_y_plates)*(ry-1)&
                +((yibte-yobte)/aero_y_plates)*(ry-1)*(rx)&
                +(-ABS(dxib-dxob)*(ry-1)/aero_y_plates+dxib)*(rX)
            
            temp_coord = mod(t, num_x_panel);
            temp_coord_d=(t-temp_coord)/num_x_panel
            
            if(rx .eq. 1 .and. ry.eq. 1) then
                x1(t) = yibte;
                x2(t) = yibte+yobte/num_x_panel;
                x3(t) = x2(t)+(dxib-(dxib-dxob)/num_x_panel)*(ry);
                x4(t) = x1(t)+dxib;
            else if (temp_coord_d .eq. 0 .and. temp_coord .ne. 0) then
                x1(t) = x2(t-1)
                x4(t) = x3(t-1)
                x2(t) = x1(t)+yobte/num_x_panel;
                x3(t) = x2(t)+dxib-(dxib-dxob)/num_x_panel*(temp_coord+1)
            endif
            
            if(temp_coord_d .ne. 0) then
                x1(t) = x4(num_x_panel*(temp_coord_d-1)+temp_coord)
                x2(t) = x3(num_x_panel*(temp_coord_d-1)+temp_coord)
                x3(t) = x2(t)+dxib-(dxib-dxob)/num_x_panel*(temp_coord-1)
                x4(t) = x1(t)+dxib-(dxib-dxob)/num_x_panel*(temp_coord)
            endif
            
            
            
            
            Cpy(t) = y1(t)+(y2(t)-y1(t))/2;
            Cpx(t) = x1(t)+(x4(t)-x1(t))/4;
            !print*, t
            
            Xcp(t,1) = Cpy(t)
            Ycp(t,1) = Cpx(t)
            
            Cp(t,1) = Cpx(t)
            Cp(t,2) = Cpy(t)
            
            lvpx(t) = x1(t)+dble(0.75)*(x4(t)-x1(t))
            lvpy(t) = y1(t)
            
            N1(:,1) = lvpx(:)
            N1(:,2) = lvpy(:)
            
            rvpx(t) = lvpx(t)
            rvpy(t) = y2(t)
            
            N2(:,1) = rvpx(:)
            
            N2(:,2) = rvpy(:)
            
            fx(t) = lvpx(t)
            fy(t) = (y2(t)+y1(t))/2
            
            NC(t,1) = fx(t)
            
            NC(t,2) = fy(t)
        end do
        
    end do
    
    aero_x_plates = num_x_panel
    
    x_step = real(real(1)/real(num_x_panel));
    
    allocate(x_2d(num_x_panel+1), increment(num_x_panel+1))
    x_2d = 0.d0; increment = 0.d0
    
    do i =1, (num_x_panel+1)
        
        x_2d(i) = length_x_total*x_step*i !! 1*17matrix
        increment(i) = 0
        !print*, 'x_2d(i) : ',x_2d(i)
    
    end do
    
    !print*, 'i : ',i
    
    allocate(X(num_y_panel+1, num_x_panel+1), Y(num_y_panel+1, num_x_panel+1))
    X = 0.d0; Y = 0.d0
    
    do i = 1,(num_y_panel+1)
        do j = 1, (num_x_panel+1)
            if (i == 1) then
                Y(i,j) = 0
            else
                Y(i, j) = x_2d(j)
            end if
            
        X(i, j) = increment(j)
        increment(j) = increment(j)+length_y_total/num_y_panel
        end do
    end do
    
    do t = 1,(aero_y_plates-1)
        do rx =1, num_y_panel-1
            !print*, 'x1 : ',x1(t)
        end do
    end do
    
    print*, 'aeropaneling for VLM done....'
    
    nc(:,2) = fy(:)
    
    deallocate(x_2d, increment, x,y,cpx,cpy)
    deallocate(rvpx,rvpy,pln1,pln2,plf)
    
end subroutine aero_paneling_vlm