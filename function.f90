double precision function f(x) result(deltav)
    use physical_param
    implicit none
    
    double precision, intent(in) :: x(n)
    
    double precision :: v0, lambda1, rpm
    
    v0 = x(1)
    lambda1 = x(2)
    
    call patched_conic(v0, lambda1, rpm, deltav)
end function