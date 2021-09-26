! Resolve sistema linear Ax=b, com A(n,n), b(n,1)

subroutine linear_system(A,b,x,n)
    implicit none

    integer, intent(in) :: n
    double precision, intent(in) :: A(n,n), b(n)
    double precision, intent(out) :: x(n)

    integer :: i, j, k
    double precision :: m, S, A_tmp(n,n), b_tmp(n)

    A_tmp = A
    b_tmp = b

    ! Triangularizacao
    do k = 1,n-1
        do i = k+1,n
            m = A_tmp(i,k)/A_tmp(k,k)
            A_tmp(i,k) = 0

            do j = k+1,n
                A_tmp(i,j) = A_tmp(i,j) - m*A_tmp(k,j)
            end do
            b_tmp(i) = b_tmp(i) - m*b_tmp(k)
        end do
    end do

    ! Substituicao inversa
    do i = n,1,-1
        S = 0
        do j = i+1,n
            S = S + A_tmp(i,j)*x(j)
        end do
        x(i) = (b_tmp(i) - S)/A_tmp(i,i)
    end do
end subroutine