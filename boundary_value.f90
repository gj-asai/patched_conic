! Resolve o problema para diversos lambdas
    
subroutine boundary_value()
    implicit none
    
    common /angles/ deg2rad, rad2deg
    double precision :: deg2rad, rad2deg
    
    common /physical/ D, mu_e, Re, mu_m, Rm, vm, wm, Rs, r0, rf, phi0
    double precision :: D, mu_e, Re, mu_m, Rm, vm, wm, Rs, r0, rf, phi0
    
    double precision :: v0, lambda1, rpm, deltav
    integer :: i
    
    call hohmann(v0) ! Chute inicial para v0, proximas iteracoes usam o valor da iteracao anterior
    
    open(1, file = 'antihorario.dat')
    write(1,'(2a25)'), 'lambda1 [graus]', 'deltav [km/s]'
    
    do i = 30,83
        lambda1 = i*deg2rad
        call newton_raphson(r0, rf, v0, lambda1)
        call patched_conic(v0, lambda1, rpm, deltav)
        
        write(1,'(2f)'), lambda1*rad2deg, deltav
    end do
    
    close(1)
end subroutine