! Encontra o lambda de minimo delta v
    
subroutine optimize()
    implicit none
    
    common /angles/ deg2rad, rad2deg
    double precision :: deg2rad, rad2deg
    
    common /physical/ D, mu_e, Re, mu_m, Rm, vm, wm, Rs, r0, rf, phi0
    double precision :: D, mu_e, Re, mu_m, Rm, vm, wm, Rs, r0, rf, phi0
    
    common /sgra/ n, q, maxiter, e1, e2, theta2, theta3
    integer :: n, q, maxiter
    double precision :: e1, e2, theta2, theta3
    
    common /optim/ clockwise, v0, lambda1
    logical :: clockwise
    double precision :: v0, lambda1
    
    double precision :: deltav, deltav1, deltav2, deltat, gamma0, energym, em
    double precision :: x(n)
    
    ! Aumenta precisao do v0
    call newton_raphson(clockwise, r0, rf, v0, lambda1)
    
    ! Otimizacao
    x(1) = v0
    x(2) = lambda1
    call sgra_conjugate(x)
    
    ! Recalcula e apresenta resultado da otimizacao
    ! Tambem gera arquivo com a trajetoria otima
    call patched_conic_complete(clockwise, x(1), x(2), deltav, deltav1, deltav2, deltat, gamma0, energym, em, .true.)
    print *, 'v0      =', v0, 'km/s'
    print *, 'lambda1 =', lambda1*rad2deg, 'graus'
    print *, 'gamma0  =', gamma0*rad2deg, 'graus'
    print *, ''
    print *, 'deltav  =', deltav, 'km/s'
    print *, 'deltav1 =', deltav1, 'km/s'
    print *, 'deltav2 =', deltav2, 'km/s'
    print *, ''
    print *, 'deltat  =', deltat/3600/24, 'dias'
    print *, ''
    print *, 'energym =', energym, 'km^2/s^2'
    print *, 'em      =', em

end subroutine