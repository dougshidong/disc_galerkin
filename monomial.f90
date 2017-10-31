module monomial
use prec
use glob

contains
subroutine monomialP(x, order, P)
implicit none
real(dp),   intent(out) :: P(:)
real(dp),   intent(in)  :: x(:)
integer,    intent(in)  :: order
P = x**order
return
end subroutine monomialP

subroutine gradMonomialP(x, order, P)
implicit none
real(dp),   intent(out) :: P(:)
real(dp),   intent(in)  :: x(:)
integer,    intent(in)  :: order
P = 0.0d0
if(order.eq.0) return
P = order*x**(order-1)
return
end subroutine gradMonomialP


end module monomial
