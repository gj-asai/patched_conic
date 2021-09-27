program optimization
    implicit none
    
    common /angles/ deg2rad, rad2deg
    double precision :: deg2rad = 3.141592653589/180, rad2deg = 180/3.141592653589
    
    common /bvp/ dv0, tol
    double precision :: dv0 = 1d-6, tol = 1d-12 ! Passo para derivacao, tolerancia para Newton-Raphson
    
    common /sgra/ n, q, maxiter, e1, e2, theta2, theta3
    integer :: n = 2, q = 1, maxiter = 100 ! Dimensoes, restricoes e maximo de iteracoes
    double precision :: e1 = 1d-6          ! Criterio de parada no passo da fase do gradiente: psi_a^2(a)/psi_a^2(0) <= e1
    double precision :: e2 = 1e-6          ! Passo para derivada por diferencas centradas (gradientes e calculo do passo)
    double precision :: theta2 = 1e-15     ! Criterio de parada da fase de restauracao: P(x) <= theta2, P(x) = phi_T(x) * phi(x)
    double precision :: theta3 = 1e-12     ! Criterio de parada do algoritmo: Q(x) <= theta3, Q(x) = F_x_T(x,lambda) * F_x(x,lambda)
    
    call physical_param(167d0, 100d0) ! Altura das orbitas inicial e final
    
    call boundary_value() ! Constroi curva lambda1 vs deltav
    call optimize()       ! Encontra lambda1 de deltav minimo
end program