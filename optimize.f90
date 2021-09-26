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
    
    double precision :: v0, lambda1, deltav, deltav1, deltav2, deltat, gamma0
    double precision :: x(n)
    
    lambda1 = 60*deg2rad ! Valor para inicio da otimizacao
    call hohmann(v0) ! Chute para v0
    call newton_raphson(r0, rf, v0, lambda1) ! Calcula v0 correto para lambda1
    
    ! Otimizacao
    x(1) = v0
    x(2) = lambda1
    call sgra_conjugate(x)
    
    ! Recalcula e apresenta resultado da otimizacao
    call plot_trajectory(x(1), x(2), deltav, deltav1, deltav2, deltat, gamma0)
    print *, 'v0      =', x(1), 'km/s'
    print *, 'lambda1 =', x(2)*rad2deg, 'graus'
    print *, 'gamma0  =', gamma0*rad2deg, 'graus'
    print *, ''
    print *, 'deltav  =', deltav, 'km/s'
    print *, 'deltav1 =', deltav1, 'km/s'
    print *, 'deltav2 =', deltav2, 'km/s'
    print *, ''
    print *, 'deltat  =', deltat/3600/24, 'dias'
end subroutine