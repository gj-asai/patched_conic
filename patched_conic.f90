! Retorna r_pm, a distancia do periselenio para os inputs dados
! Retorna tambem o deltav total para a transferencia
! Depois deve ser verificado se e igual a orbita final desejada

subroutine patched_conic(v0, lambda1, rpm, deltav)
    use physical_param
    implicit none

    double precision, intent(in) :: v0, lambda1
    double precision, intent(out) :: rpm, deltav
    
    double precision :: energye, he, r1, r2, gamma1, v1, phi1
    double precision :: epsilon2, v2
    double precision :: energym, hm, pm, em, vpm
    double precision :: deltav1, deltav2
    
    ! Orbita geocentrica
    energye = v0**2/2 - mu_e/r0
    he = r0*v0*cos(phi0)
    
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
    
    rpm = pm/(1+em)
    vpm = sqrt(2*(energym + mu_m/rpm))
    
    ! Variacao de velocidade nos impulsos
    deltav1 = v0 - sqrt(mu_e/r0)
    deltav2 = vpm - sqrt(mu_m/rf)
    deltav = deltav1 + deltav2
end subroutine