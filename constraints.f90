! Restricoes do problema
function phi(x)
    implicit none
    
    common /physical/ D, mu_e, Re, mu_m, Rm, vm, wm, Rs, r0, rf, phi0
    double precision :: D, mu_e, Re, mu_m, Rm, vm, wm, Rs, r0, rf, phi0
    
    common /arrival/ clockwise
    logical :: clockwise
    
    double precision, intent(in) :: x(n)
    double precision :: phi(q)
    
    double precision :: v0, lambda1, rpm, deltav
    
    v0 = x(1)
    lambda1 = x(2)

    call patched_conic(clockwise, v0, lambda1, rpm, deltav)
    phi(1) = rpm - rf
end function