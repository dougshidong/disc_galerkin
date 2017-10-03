subroutine legendreGLNodesWeights(N, x, w)
use glob
use prec

implicit none

integer, intent(in) :: N
real(dp), intent(out) :: x(0:N), w(0:N)

integer :: j
integer :: k, kmax = 20000
real(dp) :: dx, q, qd, LN

if(N.eq.1) then
    x(0) = -1.0d0
    w(0) = 1.0d0
    x(1) = 1.0d0
    w(1) = w(0)
else
    x(0) = -1.0d0
    w(0) = 2.0d0 / (N*(N+1.0d0))
    x(N) = 1.0d0
    w(N) = w(0)

    do j = 1, (N+1)/2 - 1
        ! Initial guess
        x(j) = -cos( (j+0.25)*PI/N - 3.0d0 / (8.0d0*N*PI*(j+0.25d0)) )
        ! Newton's method for root finding
        do k = 0, kmax
            call legendreNM(N, x(j), q, qd, LN)
            dx = -q/qd
            x(j) = x(j) + dx
            if(abs(dx) .lt. 1e-18*abs(x(j))) exit
        end do
        call legendreNM(N, x(j), q, qd, LN)
        x(N-j) = -x(j)
        w(j) = 2.0d0 / (N*(N+1.0d0)*LN*LN)
        w(N-j) = w(j)
    end do
end if
if(MOD(N,2).eq.0) then
    j = N/2
    x(j) = 0.0d0
    call legendreNM(N, x(j), q, qd, LN)
    w(j) = 2.0d0 / (N*(N+1.0d0)*LN*LN)
end if

return

end subroutine legendreGLNodesWeights
