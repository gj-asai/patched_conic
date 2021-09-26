! Dados do problema para o SGRA:
! n: integer, numero de dimensoes
! q: integer, numero de restricoes
! maxiter: integer, numero maximo de iteracoes
! x0: double precision(n), ponto inicial, nao precisa satisfazer as restricoes
! e1: double precision, criterio de parada no passo da fase do gradiente, psi_a^2(a)/psi_a^2(0) <= e1
! e2: double precision, valor pequeno para estimar derivadas (gradientes e calculo do passo)
! theta2: double precision, criterio de parada da fase de restauracao, P(x) <= theta2, P(x) = phi_T(x) * phi(x)
! theta3: double precision, criterio de parada, Q(x) <= theta3, Q(x) = F_x_T(x,lambda) * F_x(x,lambda)

! f: function (x), funcao a ser otimizada - arquivo function.f90
! phi, dimension(q): function(x), restricoes - arquivo constraints.f90

module optim_config
    implicit none

    integer, parameter ::  n=2, q=1, maxiter=100
    double precision, parameter :: e1=1d-6, e2=1d-6, theta2=1d-15, theta3=1d-12

    contains
        include 'function.f90'
        include 'constraints.f90'

        ! Gradiente da funcao
        function grad_f(x)
            double precision, intent(in) :: x(n)
            double precision :: grad_f(n)
            
            double precision :: x1(n), x2(n)
            integer :: i
            
            do i = 1,n
                x1 = x
                x2 = x
    
                x1(i) = x(i) - e2
                x2(i) = x(i) + e2
    
                grad_f(i) = (f(x2) - f(x1))/(2*e2)
            end do
        end function

        ! Gradiente das restricoes
        function grad_phi(x)
            double precision, intent(in) :: x(n)
            double precision :: grad_phi(n,q)
            
            double precision :: x1(n), x2(n), y1(q), y2(q)
            integer :: i, j
            
            do i = 1,n
                x1 = x
                x2 = x
                    
                x1(i) = x(i) - e2
                x2(i) = x(i) + e2
                    
                y1 = phi(x1)
                y2 = phi(x2)
                
                do j = 1,q
                    grad_phi(i,j) = (y2(j) - y1(j))/(2*e2)
                end do
            end do
        end function

        ! Funcao aumentada para calculo do passo do gradiente
        double precision function f_aug(x,lambda)
            double precision, intent(in) :: x(n), lambda(q)
            f_aug = f(x) + dot_product(lambda,phi(x))
        end function

        ! Gradiente da funcao aumentada
        function grad_f_aug(x,lambda)
            double precision, intent(in) :: x(n), lambda(q)
            double precision :: grad_f_aug(n)
            grad_f_aug = grad_f(x) + matmul(grad_phi(x),lambda)
        end function
end module