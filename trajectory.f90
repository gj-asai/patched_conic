! Encontra o lambda de minimo delta v
    
subroutine trajectory(clockwise,lambda1,v0)
    implicit none
    
    common /physical/ D, mu_e, Re, mu_m, Rm, vm, wm, Rs, r0, rf, phi0
    double precision :: D, mu_e, Re, mu_m, Rm, vm, wm, Rs, r0, rf, phi0
    
    logical, intent(in) :: clockwise
    double precision, intent(in) :: lambda1
    double precision :: v0
    
    double precision :: deltav, deltav1, deltav2, deltat, gamma0, energym, em
    
    ! Aumenta precisao do v0
    call newton_raphson(clockwise, r0, rf, v0, lambda1)

    call patched_conic_complete(clockwise, v0, lambda1, deltav, deltav1, deltav2, deltat, gamma0, energym, em, .true.)
end subroutine