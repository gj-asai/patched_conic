subroutine newton_raphson(r0, rf, v0, lambda1)
    implicit none
    
    double precision, intent(in) :: r0, rf, lambda1
    double precision, intent(inout) :: v0
    
    double precision :: rpm, deltav
    double precision :: v1, v2, rpm1, rpm2, d_rpm, deltav0
    
    common /bvp/ dv0, tol
    double precision :: dv0, tol
    
    do
        call patched_conic(v0, lambda1, rpm, deltav)
        
        ! Calculo da derivada d(rpm)/d(v0)
        ! Se a velocidade testada nao tem energia suficiente, d_rpm = NaN
        v1 = v0 - dv0
        v2 = v0 + dv0
        
        call patched_conic(v1, lambda1, rpm1, deltav)
        call patched_conic(v2, lambda1, rpm2, deltav)
        d_rpm = (rpm2 - rpm1)/(2*dv0)
        
        ! Calculo do novo v0
        ! Se d_rpm = NaN, somente aumenta a velocidade
        if (isnan(d_rpm)) then
            v0 = 1.05*v0
        else
            deltav0 = (rf - rpm)/d_rpm
            v0 = v0 + deltav0
        end if
        
        if (abs(deltav0).lt.tol) exit
    end do
end subroutine