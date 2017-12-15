module jacobi
use prec
use glob

real(dp) :: jacobi_alpha = 0.0, jacobi_beta = 0.0

contains

subroutine JacobiP(x, alpha, beta, order, P)
implicit none
real(dp),   intent(out) :: P(:)
real(dp),   intent(in)  :: x(:), alpha, beta
integer,    intent(in)  :: order

real(dp) :: PL(order+1, size(x))
integer :: i
real(dp) :: gamma0, gamma1, aold, anew, bnew, h1

if(size(x).ne.size(P)) print*,'size wrong,', size(x), size(P)
!P0
gamma0 = 2.0d0**(alpha+beta+1) * gamma(alpha+1) * gamma(beta+1.0d0) / gamma(alpha+beta+2)
PL(1,:) = 1.0d0 / sqrt(gamma0)
if(order == 0) P(:) = PL(order+1,:)
if(order == 0) return

!P1
gamma1 = (alpha+1) * (beta+1) * gamma0 / (alpha+beta+3.0d0)
PL(2,:) = ( (alpha+beta+2) * x(:)/2.0d0 + (alpha-beta) / 2.0d0 ) / sqrt(gamma1)
if(order == 1) P(:) = PL(order+1,:)
if(order == 1) return

aold = 2.0d0 * sqrt((alpha+1.0d0) * (beta+1.0d0) / (alpha+beta+3.0d0)) / (2.0d0+alpha+beta)

do i = 1, order - 1
    h1 = 2*i + alpha + beta
    anew = 2.0d0 / (h1+2) * &
        sqrt( (i+1) * (i+1+alpha+beta) &
            * (i+1+alpha) * (i+1+beta) &
            / (h1+1) / (h1+3) )
    bnew = - (alpha**2.0d0-beta**2.0d0) / h1 / (h1+2.0d0)
    PL(i+2,:) = 1.0d0 / anew * ( -aold*PL(i,:) + (x-bnew) * PL(i+1,:))
    aold = anew;
end do

P(:) = PL(order+1,:)

return

end subroutine JacobiP

subroutine gradJacobiP(r, alpha, beta, order, Pd)
implicit none
integer, intent(in) :: order
real(dp), intent(in) :: r(:), alpha, beta
real(dp), intent(out) :: Pd(:)!size(r))
real(dp) :: P(size(r))

if(order.eq.0) then
    Pd = 0.0d0
else
    call JacobiP(r, alpha+1, beta+1, order-1, P)
    Pd = sqrt(order*(order+alpha+beta+1))*P
end if

return
end subroutine gradJacobiP

end module jacobi
