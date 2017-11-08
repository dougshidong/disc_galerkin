module bspline

use prec
use glob
use quadrature
implicit none

integer :: bs_order, bs_nbasis
logical :: bs_init = .false.
integer :: nknots
integer :: openknot = 1
real(dp), allocatable :: knots(:)

contains
subroutine bspline_init
implicit none
integer :: i
real(dp) :: dk

if(bs_init) return
bs_init = .true.
bs_nbasis = order+1
bs_order = max(bs_nbasis-14,2)
bs_order = order
!bs_nbasis = order+1
!bs_order = bs_nbasis-5

nknots = bs_order+bs_nbasis+1
!if(bs_nbasis.gt.1) dk = (refb-refa) / (bs_nbasis-1.0d0)
if(bs_nbasis.gt.1) dk = (refb-refa) / (nknots-2*(bs_order+1)+1)
!if(bs_nbasis.gt.1) dk = (refb-refa) / (nknots-2*(bs_order+1)+3)

allocate(knots(nknots))
! Initialize knot vector
knots = -99999
knots(1:bs_order+1) = refa
knots(nknots-bs_order:nknots) = refb+1e-15
do i = 1, nknots-2*(bs_order+1)
    knots(i+bs_order+1) = refa + dk*i
end do
!print*,knots

end subroutine bspline_init

subroutine bsplineP(x, ith, bs_order, B)
real(dp),   intent(out) :: B(:)
real(dp),   intent(in)  :: x(:)
integer,    intent(in)  :: ith, bs_order

B = 0
if(ith.lt.1 .or. bs_order.lt.0) return
call coxdeBoor(x, ith, bs_order, B)

return
end subroutine bsplineP

subroutine gradbsplineP(x, ith, bs_order, B)
real(dp),   intent(out) :: B(:)
real(dp),   intent(in)  :: x(:)
integer,    intent(in)  :: ith, bs_order
real(dp), dimension(size(x)) :: P1, P2
real(dp)  :: den1, den2

B = 0
call bsplineP(x, ith, bs_order-1, P1)
call bsplineP(x, ith+1, bs_order-1, P2)
den1 = knots(ith+bs_order) - knots(ith)
den2 = knots(ith+bs_order+1) - knots(ith+1)
if(den1.ne.0) B = B+P1/den1
if(den2.ne.0) B = B-P2/den2
B = bs_order*B

return
end subroutine gradbsplineP

recursive subroutine coxdeBoor(x, ith, bs_order, B)
real(dp),   intent(out) :: B(:)
real(dp),   intent(in)  :: x(:)
integer,    intent(in)  :: ith, bs_order
integer :: i, nx
real(dp), allocatable :: P1(:), P2(:)
real(dp)  :: den1, den2

nx = size(x)
allocate(P1(nx),P2(nx))

if(bs_order.eq.0) then
    do i = 1, nx
!       if(ith.ge.nknots-bs_order) then
!           B(i) = 0
!           cycle
!       endif
       
        if(x(i).ge.knots(ith) .and. x(i).lt.knots(ith+1)) then
            B(i) = 1
        else
            B(i) = 0
        end if
    end do
    return
end if

call coxdeBoor(x,ith  ,bs_order-1,P1)
call coxdeBoor(x,ith+1,bs_order-1,P2)
B = 0
den1 = knots(ith+bs_order) - knots(ith)
den2 = knots(ith+bs_order+1) - knots(ith+1)
if(den1.ne.0) B = B+P1/den1*(x-knots(ith))
if(den2.ne.0) B = B+P2/den2*(knots(ith+bs_order+1)-x)

return
end subroutine coxdeBoor



end module bspline

