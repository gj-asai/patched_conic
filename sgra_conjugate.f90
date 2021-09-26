subroutine sgra_conjugate(x)
    use optim_config
    implicit none
    
    double precision, intent(inout) :: x(n)

    integer :: iter, deltaN
    double precision :: alpha, f1, f2, err
    double precision :: xold(n), start(n), p(n), lambda(q), lambda_old(q)
    logical :: restart

    print *, '***** SGRA *****'
    
    xold = x
    lambda_old = 0
    
    p = 0
    lambda = 0
    
    ! O algoritmo sera reiniciado a cada deltaN iteracoes
    deltaN = n - q

    do iter = 1,maxiter
        start = x
        f1 = f(x)

        restart = mod(iter,deltaN).eq.0

        ! Calcula o valor de lambda para a posicao atual de x
        call lagrange_mult(lambda,lambda_old,x)
        err = dot_product(grad_f_aug(x,lambda),grad_f_aug(x,lambda))
        if (err.le.theta3) exit

        ! Calcula a direcao p e o valor de alpha otimo
        call gradient(xold,x,p,alpha,lambda,lambda_old)

        ! Se alpha otimo nao funciona (avanca demais e a funcao aumenta), bissecta alpha ate funcionar, mas depois tem que reiniciar o algoritmo
        do
            x = start - alpha*p
            call restoration(x)
            
            f2 = f(x)
            if (f2.lt.f1) exit
            
            alpha = alpha/2
            restart = .true.
        end do
        xold = start

        if (restart) p = 0
    end do

    print *, 'Iteracoes:', iter-1
    print *, ''
end subroutine