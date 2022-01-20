! Resolve o problema para diversos lambdas
    
subroutine boundary_value()
    implicit none
    
    common /angles/ deg2rad, rad2deg
    double precision :: deg2rad, rad2deg
    
    common /physical/ D, mu_e, Re, mu_m, Rm, vm, wm, Rs, r0, rf, phi0
    double precision :: D, mu_e, Re, mu_m, Rm, vm, wm, Rs, r0, rf, phi0
    
    double precision :: v0, lambda1, deltav, deltav1, deltav2, deltat, gamma0, energym, em
    double precision :: i
    
    ! Chegada horaria
    call hohmann(v0)
    v0 = v0 + 0.25
    
    open(10, file = 'horario.dat')
    write(10,'(a15,8a25)'), 'lambda1 [graus]', 'v0 [km/s]', 'deltav [km/s]', 'deltav1 [km/s]', 'deltav2 [km/s]', 'deltat [dias]', 'gamma0 [graus]', 'energym [km^2/s^2]', 'em'
    
    do i = 20,90,0.01
        lambda1 = i*deg2rad
        call newton_raphson(.true., r0, rf, v0, lambda1)
        call patched_conic_complete(.true., v0, lambda1, deltav, deltav1, deltav2, deltat, gamma0, energym, em, .false.)
        
        if (deltav .lt. 0) exit
        
        write(10,'(f15.2,8f)'), i, v0, deltav, deltav1, deltav2, deltat/3600/24, gamma0*rad2deg, energym, em
        print '(a15,f15.2)', 'HORARIO', i
    end do
    close(10)
    
    ! Chegada antihoraria
    call hohmann(v0)
    
    open(11, file = 'antihorario.dat')
    write(11,'(a15,8a25)'), 'lambda1 [graus]', 'v0 [km/s]', 'deltav [km/s]', 'deltav1 [km/s]', 'deltav2 [km/s]', 'deltat [dias]', 'gamma0 [graus]', 'energym [km^2/s^2]', 'em'

    do i = 20,90,0.01
        lambda1 = i*deg2rad
        call newton_raphson(.false., r0, rf, v0, lambda1)
        call patched_conic_complete(.false., v0, lambda1, deltav, deltav1, deltav2, deltat, gamma0, energym, em, .false.)
        
        if (deltav .lt. 0) exit
        
        write(11,'(f15.2,8f)'), i, v0, deltav, deltav1, deltav2, deltat/3600/24, gamma0*rad2deg, energym, em
        print '(a15,f15.2)', 'ANTIHORARIO', i
    end do
    
    close(11)
end subroutine