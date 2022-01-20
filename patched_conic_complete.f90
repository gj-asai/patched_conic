! Patched conic para quando o par (v0,lambda1) e conhecido
! Retorna os parametros para pos processamento

subroutine patched_conic_complete(clockwise, v0, lambda1, deltav, deltav1, deltav2, deltat, gamma0, energym, em, trajectory)
    implicit none
    
    logical, intent(in) :: clockwise, trajectory
    double precision, intent(in) :: v0, lambda1
    double precision, intent(out) :: deltav, deltav1, deltav2, deltat, gamma0, energym, em
    
    common /angles/ deg2rad, rad2deg
    double precision :: deg2rad, rad2deg
    
    common /physical/ D, mu_e, Re, mu_m, Rm, vm, wm, Rs, r0, rf, phi0
    double precision :: D, mu_e, Re, mu_m, Rm, vm, wm, Rs, r0, rf, phi0
    
    double precision :: energye, he, pe, ae, ee, r1, r2, gamma1, v1, phi1
    double precision :: phi2, v2
    double precision :: hm, pm, am, vpm, rpm
    double precision :: E0, E1, E2, H2, M2, deltat_e, deltat_m
    double precision :: f1, f2, deltaE, deltaH, w_pm, Ei, Hi, ri, fi, xi, yi, Mi, ti
    
    integer :: i, npontos
    
    ! Orbita geocentrica
    energye = v0**2/2 - mu_e/r0
    he = r0*v0*cos(phi0)
    
    pe = he**2/mu_e
    ee = sqrt(1+2*energye*he**2/mu_e**2)
    ae = -mu_e/2/energye
    
    ! Transicao
    r2 = rs
    
    r1 = sqrt(D**2 + r2**2 - 2*D*r2*cos(lambda1))
    gamma1 = asin(r2*sin(lambda1)/r1)
    
    v1 = sqrt(2*(energye + mu_e/r1))
    phi1 = acos(he/r1/v1)
    
    v2 = sqrt(v1**2 + vm**2 - 2*v1*vm*cos(phi1-gamma1))
    
    if (clockwise) then
        phi2 = -lambda1 + atan(-v1*sin(phi1-gamma1)/(vm-v1*cos(phi1-gamma1)))
    else
        phi2 = lambda1 - atan(-v1*sin(phi1-gamma1)/(vm-v1*cos(phi1-gamma1)))
    end if
    
    ! Orbita selenocentrica
    energym = v2**2/2 - mu_m/r2
    hm = r2*v2*cos(phi2)
    
    pm = hm**2/mu_m
    em = sqrt(1+2*energym*hm**2/mu_m**2)
    am = -mu_m/2/energym
    
    rpm = pm/(1+em)
    vpm = sqrt(2*(energym + mu_m/rpm))
    
    ! Variacao de velocidade nos impulsos
    deltav1 = v0 - sqrt(mu_e/r0)
    deltav2 = vpm - sqrt(mu_m/rf)
    deltav = deltav1 + deltav2
    
    ! Tempo de viagem
    E1 = acos(1/ee * (1-r1/ae))
    deltat_e = sqrt(ae**3/mu_e) * (E1-ee*sin(E1))
    
    if (energym.lt.0) then ! Orbita selenocentrica eliptica
        E2 = -acos(1/em * (1-r2/am))
        deltat_m = sqrt(am**3/mu_m) * (-E2+em*sin(E2))
    else                   ! Orbita selenocentrica hiperbolica
        H2 = -acosh(1/em * (1-r2/am))
        deltat_m = sqrt(-am**3/mu_m) * (-em*sinh(H2)+H2)
    end if
    
    deltat = deltat_e + deltat_m
    
    ! Angulo de saida da terra
    f1 = 2*atan(sqrt((1+ee)/(1-ee))*tan(E1/2)) 
    gamma0 = -f1 + gamma1 + wm*deltat_e
    
    if (.not.trajectory) return
    
    ! Salva parametros importantes
    open(1, file='resultados.dat')
    write(1,*), 'v0      =', v0, 'km/s'
    write(1,*), 'lambda1 =', lambda1*rad2deg, 'graus'
    write(1,*), 'gamma0  =', gamma0*rad2deg, 'graus'
    write(1,*), ''
    write(1,*), 'deltav  =', deltav, 'km/s'
    write(1,*), 'deltav1 =', deltav1, 'km/s'
    write(1,*), 'deltav2 =', deltav2, 'km/s'
    write(1,*), ''
    write(1,*), 'deltat  =', deltat/3600/24, 'dias'
    write(1,*), ''
    write(1,*), 'energym =', energym, 'km^2/s^2'
    write(1,*), 'em      =', em
    close(1)
    
    ! Plot da trajetoria
    ! Referencias:
    ! 0: Inicio da trajetoria geocentrica
    ! 1: Fim da trajetoria geocentrica
    ! 2: Inicio da trajetoria selenocentrica
    open(2, file = 'trajetoria.dat')
    write(2,'(3a25)'), 'X [km]', 'Y [km]', 'T [dias]'
    npontos = 1000
    
    ! Trajetoria geocentrica
    deltaE = E1/(npontos-1)
    do i = 0,npontos-1
        Ei = deltaE*i
        ri = ae*(1-ee*cos(Ei))
        fi = 2*atan(sqrt((1+ee)/(1-ee))*tan(Ei/2))
        
        xi = ri*cos(gamma0 + fi)
        yi = ri*sin(gamma0 + fi)
        
        Mi = Ei - ee*sin(Ei)
        ti = Mi*sqrt(ae**3/mu_e)
        
        write(2,'(3f)'), xi, yi, ti/3600/24
    end do
    
    ! Trajetoria selenocentrica
    if (energym.lt.0) then
        M2 = E2 - em*sin(E2)
        f2 = 2*atan(sqrt((1+em)/(1-em))*tan(E2/2))
        if (clockwise) then
            w_pm = lambda1 - f2 + 180*deg2rad - wm*deltat_e
        else
            w_pm = -lambda1 - f2 + 180*deg2rad + wm*deltat_e
        end if
        
        deltaE = -E2/npontos
        
        do i = 1,npontos
            Ei = E2 + deltaE*i
            ri = am*(1-em*cos(Ei))
            fi = 2*atan(sqrt((1+em)/(1-em))*tan(Ei/2))
            
            Mi = Ei - em*sin(Ei)
            ti = deltat_e + (Mi-M2)*sqrt(am**3/mu_m)
            
            if (clockwise) then
                xi = ri*cos(-w_pm - fi) + D*cos(wm*ti)
                yi = ri*sin(-w_pm - fi) + D*sin(wm*ti)
            else
                xi = ri*cos(w_pm + fi) + D*cos(wm*ti)
                yi = ri*sin(w_pm + fi) + D*sin(wm*ti)
            end if
            
            write(2,'(3f)'), xi, yi, ti/3600/24
        end do
    else
        M2 = em*sinh(H2) - H2
        f2 = 2*atan(sqrt((1+em)/(em-1))*tanh(H2/2))
        if (clockwise) then
            w_pm = lambda1 - f2 + 180*deg2rad - wm*deltat_e
        else
            w_pm = -lambda1 - f2 + 180*deg2rad + wm*deltat_e
        end if
        
        deltaH = -H2/npontos
        
        do i = 1,npontos
            Hi = H2 + deltaH*i
            ri = am*(1-em*cosh(Hi))
            fi = 2*atan(sqrt((1+em)/(em-1))*tanh(Hi/2))
            
            Mi = em*sinh(Hi) - Hi
            ti = deltat_e + (Mi-M2)*sqrt(-am**3/mu_m)
            
            if (clockwise) then
                xi = ri*cos(-w_pm - fi) + D*cos(wm*ti)
                yi = ri*sin(-w_pm - fi) + D*sin(wm*ti)
            else
                xi = ri*cos(w_pm + fi) + D*cos(wm*ti)
                yi = ri*sin(w_pm + fi) + D*sin(wm*ti)
            end if
            
            write(2,'(3f)'), xi, yi, ti/3600/24
        end do
    end if
    
    close(2)
end subroutine