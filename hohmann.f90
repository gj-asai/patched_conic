! Calcula o v0 necessario para orbita de transferencia de Hohmann
! Utilizado como chute inicial para Newton-Raphson do ponto inicial da otimizacao

subroutine hohmann(v0)
    implicit none
    
    double precision, intent(out) :: v0
    
    common /physical/ D, mu_e, Re, mu_m, Rm, vm, wm, Rs, r0, rf, phi0
    double precision :: D, mu_e, Re, mu_m, Rm, vm, wm, Rs, r0, rf, phi0
    
    v0 = sqrt(mu_e/r0)*sqrt(2*D/(r0+D))
end subroutine