program optimization
    implicit none
    
    common /angles/ deg2rad, rad2deg
    double precision :: deg2rad = 3.141592653589/180., rad2deg = 180./3.141592653589
    
    common /bvp/ dv0, tol
    double precision :: dv0 = 1d-6, tol = 1d-12 ! Passo para derivacao, tolerancia para Newton-Raphson
    
    common /sgra/ n, q, maxiter, e1, e2, theta2, theta3
    integer :: n = 2, q = 1, maxiter = 100 ! Dimensoes, restricoes e maximo de iteracoes
    double precision :: e1 = 1d-6          ! Criterio de parada no passo da fase do gradiente: psi_a^2(a)/psi_a^2(0) <= e1
    double precision :: e2 = 1e-6          ! Passo para derivada por diferencas centradas (gradientes e calculo do passo)
    double precision :: theta2 = 1d-15     ! Criterio de parada da fase de restauracao: P(x) <= theta2, P(x) = phi_T(x) * phi(x)
    double precision :: theta3 = 1d-12     ! Criterio de parada do algoritmo: Q(x) <= theta3, Q(x) = F_x_T(x,lambda) * F_x(x,lambda)
    
    common /optim/ clockwise, v0, lambda1
    logical :: clockwise
    double precision :: v0, lambda1
    
    clockwise = .true.         ! Sentido da chegada na otimizacao
    lambda1 = 25*deg2rad       ! lambda1 inicial para otimizacao
    v0 = 1.1622238425590812    ! velocidade correspondente ao lambda1 (pode ser retirado do resultado do problema de valor de contorno)
    
    call physical_param(463d0, 500d0) ! Alturas das orbitas inicial e final
    !call physical_param(167d0/6.3781366d3*1.1883d3, 100d0/1.7375d3*606.) ! Alturas das orbitas inicial e final
    
    call boundary_value()                   ! constroi curvas em funcao de lambda1
    !call optimize()                         ! encontra lambda1 de deltav minimo
    !call trajectory(clockwise, lambda1, v0) ! plota uma trajetoria desejada
end program