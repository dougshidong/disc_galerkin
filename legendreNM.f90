subroutine legendreNM(N,x,q,qd,LN)
! Newton's method to find roots of poly q(x) = L_(N+1) - L_(N-1)
! Returns the poly q, its derivative qd and the Legendre poly L_(N)
! Korpriva, Algorithm 24
use prec
implicit none

integer, intent(in) :: N
real(dp), intent(in) :: x
real(dp), intent(out) :: q, qd, LN

integer :: k
real(dp) :: L(N-2:N+1), LD(N-2:N+1)

L(N-2) = 1.0d0
L(N-1) = x
LD(N-2) = 0.0d0
LD(N-1) = 1.0d0
do k = 2, N
    L(N) = (2.0d0*k-1.0d0) / k * x * L(N-1) - (k-1.0d0)/k * L(N-2)
    LD(N) = LD(N-2) + (2.0d0*k-1.0d0) * L(N-1)
    if(k.eq.N) exit
    L(N-2) = L(N-1)
    L(N-1) = L(N)
    LD(N-2) = LD(N-1)
    LD(N-1) = LD(N)
end do

k = N + 1
L(N+1) = (2.0d0*k-1.0d0) / k * x * L(N) - (k-1.0d0)/k * L(N-1)
LD(N+1) = LD(N-1) + (2.0d0*k-1.0d0) * L(N)

q = L(N+1) - L(N-1)
qd = LD(N+1) - LD(N-1)
LN = L(N)

return
end subroutine legendreNM
