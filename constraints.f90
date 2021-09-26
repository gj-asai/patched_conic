function phi(x)
    use physical_param
    implicit none
    
    double precision, intent(in) :: x(n)
    double precision :: phi(q)
    
    double precision :: v0, lambda1, rpm, deltav
    
    v0 = x(1)
    lambda1 = x(2)

    call patched_conic(v0, lambda1, rpm, deltav)
    phi(1) = rpm - rf
end function