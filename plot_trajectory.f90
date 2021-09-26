! Salva parametros da trajetoria e os pontos para plot
    
subroutine plot_trajectory(v0, lambda1, deltav, deltav1, deltav2, deltat, gamma0)
    implicit none
    
    double precision, intent(in) :: v0, lambda1
    double precision, intent(out) :: deltav, deltav1, deltav2, deltat, gamma0
    
    common /physical/ D, mu_e, Re, mu_m, Rm, vm, wm, Rs, r0, rf, phi0
    double precision :: D, mu_e, Re, mu_m, Rm, vm, wm, Rs, r0, rf, phi0
    
    double precision :: energye, he, pe, ae, ee, r1, r2, gamma1, v1, phi1
    double precision :: epsilon2, v2
    double precision :: energym, hm, pm, am, em, vpm, rpm
    double precision :: E1, F2, deltat_e, deltat_m
    double precision :: f1
    
    integer :: i
    
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
    
    epsilon2 = asin(vm/v2*cos(lambda1) - v1/v2*cos(lambda1+gamma1-phi1))
    
    ! Orbita selenocentrica
    energym = v2**2/2 - mu_m/r2
    hm = r2*v2*sin(epsilon2)
    
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
    E1 = acos(1/ee * (1 - r1/ae))
    F2 = acosh(1/em * (1 - r2/am))
    
    deltat_e = sqrt(ae**3/mu_e) * (E1 - ee*sin(E1))
    deltat_m = sqrt((-am)**3/mu_m) * (em*sinh(F2) - F2)
    deltat = deltat_e + deltat_m
    
    ! Angulo de saida da terra
    f1 = acos((cos(E1) - ee)/(1 - ee*cos(E1)))
    gamma0 = f1 - gamma1 - wm*deltat_e
    
    ! Plot da trajetoria
    open(1, file = 'trajetoria.dat')
    write(1,'(3a25)'), 'X [km]', 'Y [km]', 'T [dias]'
    
    close(1)
end subroutine