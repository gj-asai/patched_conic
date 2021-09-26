! Encontra o lambda de minimo delta v
    
subroutine optimize()
    use physical_param
    use optim_config
    implicit none
    
    double precision :: v0, lambda1, deltav, deltav1, deltav2, deltat, gamma0
    double precision :: x(n)
    
    lambda1 = 0.5 ! Valor para inicio da otimizacao
    call hohmann(v0) ! Chute para v0
    call newton_raphson(r0, rf, v0, lambda1) ! Calcula v0 correto para lambda1
    
    ! Otimizacao
    x(1) = v0
    x(2) = lambda1
    call sgra_conjugate(x)
    
    ! Recalcula e apresenta resultado da otimizacao
    call patched_conic_post(x(1), x(2), deltav, deltav1, deltav2, deltat, gamma0)
    print *, 'v0      =', x(1), 'km/s'
    print *, 'lambda1 =', x(2)*180/3.141592653589, 'graus'
    print *, 'gamma0  =', gamma0*180/3.141592653589, 'graus'
    print *, ''
    print *, 'deltav  =', deltav, 'km/s'
    print *, 'deltav1 =', deltav1, 'km/s'
    print *, 'deltav2 =', deltav2, 'km/s'
    print *, ''
    print *, 'deltat  =', deltat/3600/24, 'dias'
end subroutine