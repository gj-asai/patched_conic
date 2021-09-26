! Calcula o v0 necessario para orbita de transferencia de Hohmann
! Utilizado como chute inicial para Newton-Raphson do ponto inicial da otimizacao

subroutine hohmann(v0)
    use physical_param
    implicit none
    
    double precision, intent(out) :: v0
    
    v0 = sqrt(mu_e/r0) + sqrt(mu_e/r0)*(sqrt(2*D/(r0+D))-1)
end subroutine