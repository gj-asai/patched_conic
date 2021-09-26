! Atualiza o valor de lambda

subroutine lagrange_mult(lambda, lambda_old, x)
    use optim_funcs
    implicit none
    
    double precision, intent(in) :: x(n)
    double precision, intent(inout) :: lambda(q), lambda_old(q)

    double precision :: grad_phi_x(n,q), grad_phi_x_T(q,n), grad_f_x(n)

    lambda_old = lambda

    grad_phi_x = grad_phi(x)
    grad_phi_x_T = transpose(grad_phi_x)
    grad_f_x = grad_f(x)

    ! lambda = -(phi_x_T * phi_x)^(-1) * phi_x_T * f_x
    call linear_system(matmul(grad_phi_x_T,grad_phi_x), matmul(grad_phi_x_T,grad_f_x), lambda, q)
    lambda = -lambda
end subroutine