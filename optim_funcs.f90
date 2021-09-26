! Funcoes para o SGRA

! f: function (x), funcao a ser otimizada - arquivo function.f90
! phi, dimension(q): function(x), restricoes - arquivo constraints.f90

module optim_funcs
    implicit none
    
    common /sgra/ n, q, maxiter, e1, e2, theta2, theta3
    integer :: n, q, maxiter
    double precision :: e1, e2, theta2, theta3

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