! Funcao a ser otimizada
double precision function f(x) result(deltav)
    implicit none
    
    double precision, intent(in) :: x(n)
    
    common /arrival/ clockwise
    logical :: clockwise
    
    double precision :: v0, lambda1, rpm
    
    v0 = x(1)
    lambda1 = x(2)
    
    call patched_conic(clockwise, v0, lambda1, rpm, deltav)
end function