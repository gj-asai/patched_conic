! Parametros fisicos do problema

subroutine physical_param(he, hm)
    implicit none
    
    double precision, intent(in) :: he, hm
    
    common /physical/ D, mu_e, Re, mu_m, Rm, vm, wm, Rs, r0, rf, phi0
    double precision :: D, mu_e, Re, mu_m, Rm, vm, wm, Rs, r0, rf, phi0
    
    ! Terra - Lua
    D = 3.84400d5 ! Distancia entre os centros
        
    mu_e = 3.98600d5 ! Parametro gravitacional Terra
    Re = 6.3781366d3 ! Raio Terra
    
    mu_m = 4.902801d3 ! Parametro gravitacional Lua
    Rm = 1.7375d3 ! Raio Lua
    vm = sqrt(mu_e/D) ! Velocidade de translacao Lua
    wm = vm/D ! Velocidade angular Lua
    Rs = D*(mu_m/mu_e)**(0.4) ! Raio da esfera de influencia Lua
    !--------------------------
    
    ! Plutao - Caronte
    !D = 19591. ! distancia entre os centros
    !
    !mu_e = 0.013030d24 * 6.67408d-20 ! parametro gravitacional plutao
    !re = 1.1883d3 ! raio plutao
    !
    !mu_m = 106.17 ! parametro gravitacional caronte
    !rm = 606. ! raio caronte
    !vm = sqrt(mu_e/d) ! velocidade de translacao caronte
    !wm = vm/d ! velocidade angular caronte
    !rs = d*(mu_m/mu_e)**(0.4) ! Raio da esfera de influencia Caronte
    !--------------------------
    
    phi0 = 0.   ! Angulo na saida
    r0 = Re + he ! Raio da orbita inicial
    rf = Rm + hm ! Raio da orbita final
end subroutine