! Parametros fisicos do problema
module physical_param
    implicit none
    
    ! Terra - Lua
    double precision, parameter :: D = 3.84400d5 ! Distancia entre os centros
        
    ! Terra
    double precision, parameter :: mu_e = 3.98600d5 ! Parametro gravitacional
    double precision, parameter :: Re = 6.3781366d3 ! Raio
    
    ! Lua
    double precision, parameter :: mu_m = 4.902801d3 ! Parametro gravitacional
    double precision, parameter :: Rm = 1.7375d3 ! Raio
    double precision, parameter :: vm = sqrt(mu_e/D) ! Velocidade de translacao
    double precision, parameter :: wm = vm/D ! Velocidade angular
    double precision, parameter :: Rs = D*(mu_m/mu_e)**(0.4d0) ! Raio da esfera de influencia
    
    ! Trajetoria a ser otimizada
    double precision, parameter :: r0 = 463d0 + Re, rf = 100d0 + Rm
    double precision, parameter :: phi0 = 0d0
    
    ! Plutao - Caronte
    !double precision, parameter :: D = 19591d0 ! Distancia entre os centros
    !
    !! E: Plutao
    !double precision, parameter :: mu_e = 0.013030d24 * 6.67408d-20 ! Parametro gravitacional
    !double precision, parameter :: Re = 1.1883d3 ! Raio
    !
    !! M: Caronte
    !double precision, parameter :: mu_m = 102.3d0 ! Parametro gravitacional
    !double precision, parameter :: Rm = 603.6d0 ! Raio
    !double precision, parameter :: vm = sqrt(mu_e/D) ! Velocidade de translacao
    !double precision, parameter :: wm = vm/D ! Velocidade angular
    !double precision, parameter :: Rs = D*(mu_m/mu_e)**(0.4d0) ! Raio da esfera de influencia
    !
    !! Trajetoria a ser otimizada
    !double precision, parameter :: r0 = 400d0 + Re, rf = 200d0 + Rm
    !double precision, parameter :: phi0 = 0d0
end module    