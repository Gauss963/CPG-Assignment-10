recursive subroutine ifft(x)
    complex, intent(inout) :: x(:)
    complex, allocatable :: even(:), odd(:)
    complex :: w
    integer :: n, k

    n = size(x)
    if (n <= 1) return

    allocate(even(n/2), odd(n/2))

    ! Split into even and odd
    even = x(1:n:2)
    odd = x(2:n:2)

    ! Recursive IFFT
    call ifft(even)
    call ifft(odd)

    ! Combine
    do k = 1, n/2
        w = cmplx(cos(2.0*pi*(k-1)/n), sin(2.0*pi*(k-1)/n))
        x(k) = even(k) + w * odd(k)
        x(k + n/2) = even(k) - w * odd(k)
    end do

    deallocate(even, odd)
end subroutine ifft