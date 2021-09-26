! Resolve o problema para diversos lambdas
    
subroutine boundary_value_problem()
    use physical_param
    implicit none
    
    double precision :: v0, lambda1, rpm, deltav
    integer :: i
    
    call hohmann(v0) ! Chute inicial para v0, proximas iteracoes usam o valor da iteracao anterior
    
    open(1, file = 'antihorario.dat')
    write(1,'(2a25)'), 'lambda1 [graus]', 'deltav [km/s]'
    
    do i = 30,83
        lambda1 = i * 3.141592653589/180
        call newton_raphson(r0, rf, v0, lambda1)
        call patched_conic(v0, lambda1, rpm, deltav)
        
        write(1,'(2f)'), lambda1*180/3.141592653589, deltav
    end do
    
    close(1)
end subroutine