subroutine restoration(y)
    use optim_config
    implicit none

    double precision, intent(inout) :: y(n)
    
    double precision :: P1, P2
    double precision :: grad_phi_y(n,q), grad_phi_y_T(q,n), phi_y(q), sigma(q), deltay(n)
    
    do
        grad_phi_y = grad_phi(y)
        grad_phi_y_T = transpose(grad_phi_y)
        phi_y = phi(y)
        P1 = dot_product(phi_y,phi_y)

        if (P1.le.theta2) exit

        ! sigma = k[phi_x_T(y) * phi_x(y)]^(-1) * phi(y), k=1
        call linear_system(matmul(grad_phi_y_T,grad_phi_y), phi_y, sigma, q)
        do
            deltay = -matmul(grad_phi_y,sigma)
            P2 = dot_product(phi(y+deltay),phi(y+deltay))
            if (P2.lt.P1) exit
            sigma = sigma/2
        end do

        y = y + deltay
    end do
end subroutine