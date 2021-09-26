! Calculo do alpha que minimiza a funcao Psi(alpha) = F(y,lambda), com y = x - alpha*p
! Usa quasilinearizacao, com derivadas numericas

subroutine stepsize(alpha,x,p,lambda)
    use optim_config
    implicit none

    double precision, intent(in) :: x(n), lambda(q), p(n)
    double precision, intent(out) :: alpha

    double precision :: delta, eta, theta, psi_dot, psi_ddot, psi1, psi2
    double precision, dimension(n) :: y, dy, F_x

    F_x = grad_f_aug(x,lambda)
    
    eta = e2/norm2(f_x)

    alpha = 0d0
    theta = e1*dot_product(F_x,p)**2
    do
        y = x - alpha*p
        dy = eta*p

        ! Psi_alpha = -F_x_T(y) * p
        psi_dot = -dot_product(grad_f_aug(y,lambda), p)
        
        ! Psi_alpha_alpha = 1/(2*eta) {F_x_T(y+eta*p) - F_x_T(y-eta*p)} * p
        psi_ddot = 1/(2*eta)*dot_product(grad_f_aug(y+dy,lambda) - grad_f_aug(y-dy,lambda), p)

        if (psi_dot**2.le.theta) exit

        psi1 = f_aug(y,lambda)
        delta = -psi_dot/abs(psi_ddot)
        do
            psi2 = f_aug(x - (alpha+delta)*p, lambda)
            if (psi2.lt.psi1) exit
            delta = delta/2
        end do
        alpha = alpha + delta
    end do
end subroutine