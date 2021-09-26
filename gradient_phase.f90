! Apenas atualiza o valor de alpha
! O valor de x e atualizado em sgra_conjugate, somente se representar descida em f apos a restauracao

subroutine gradient(xold, x, p, alpha, lambda, lambda_old)
    use optim_funcs
    implicit none
    
    double precision, intent(in) :: x(n), xold(n), lambda(q), lambda_old(q)
    double precision, intent(inout) :: p(n)
    double precision, intent(out) :: alpha

    double precision :: grad_f_aug_x_old(n), grad_f_aug_x(n)
    double precision :: gamma

    grad_f_aug_x_old = grad_f_aug(xold,lambda_old)
    grad_f_aug_x = grad_f_aug(x,lambda)

    gamma = dot_product(grad_f_aug_x,grad_f_aug_x)/dot_product(grad_f_aug_x_old,grad_f_aug_x_old)
    if (norm2(grad_f_aug_x_old).eq.0) gamma = 0
    p = grad_f_aug_x + gamma*p

    call stepsize(alpha,x,p,lambda)
end subroutine