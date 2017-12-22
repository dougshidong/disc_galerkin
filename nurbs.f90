module nurbs

use prec
use glob
use bspline
implicit none

integer :: nurbs_order, nurbs_nbasis
logical :: ns_init = .false.
!integer :: nknots
!real(dp), allocatable :: knots(:)
real(dp), allocatable :: nurbs_weights(:)

real(dp) :: nurbs_denominator

contains

subroutine nurbs_init
implicit none
!real(dp) :: dk

if(ns_init) return
ns_init = .true.
call bspline_init
nurbs_nbasis = bs_nbasis
nurbs_order = bs_order
!nknots = nurbs_order+nurbs_nbasis+1
!if(nurbs_nbasis.gt.1) dk = (refb-refa) / (nknots-2*(nurbs_order+1)+1)
!
!allocate(knots(nknots))
!! Initialize knot vector
!knots = -99999
!knots(1:nurbs_order+1) = refa
!knots(nknots-nurbs_order:nknots) = refb+1e-14
!do i = 1, nknots-2*(nurbs_order+1)
!    knots(i+nurbs_order+1) = refa + dk*i
!end do

allocate(nurbs_weights(nurbs_nbasis))
nurbs_weights = 1.0d0
!do i = 1, nurbs_nbasis
!    nurbs_weights(i) = i*0.25
!end do

end subroutine nurbs_init

subroutine nurbsP(x, ith, nurbs_order, B)
implicit none
real(dp),   intent(out) :: B(:)
real(dp),   intent(in)  :: x(:)
integer,    intent(in)  :: ith, nurbs_order
real(dp), dimension(size(x)) :: P1

if(ith.lt.1 .or. nurbs_order.lt.0) return
call bsplineP(x, ith, nurbs_order, P1)
B = P1 * nurbs_weights(ith)

call nurbs_weighting(x, nurbs_nbasis, nurbs_order, P1)

B = B / P1
return
end subroutine nurbsP

subroutine gradnurbsP(x, ith, nurbs_order, B)
real(dp),   intent(out) :: B(:)
real(dp),   intent(in)  :: x(:)
integer,    intent(in)  :: ith, nurbs_order
real(dp), dimension(size(x)) :: W, P2, P3

call nurbs_weighting(x, nurbs_nbasis, nurbs_order, W)
call gradbsplineP(x, ith, nurbs_order, P2)
B = W*P2
call nurbs_weighting_grad(x, nurbs_nbasis, nurbs_order, P2)
call nurbsP(x, ith, nurbs_order, P3)
B = B + P2*P3
B = B / (W**2) * nurbs_weights(ith)

return
end subroutine gradnurbsP

subroutine nurbs_weighting(x, nurbs_nbasis, nurbs_order, B)
! Weighting function W(x) = SUM[ w_{i} * N_{i,p} ]
implicit none
real(dp),   intent(out) :: B(:)
real(dp),   intent(in)  :: x(:)
integer,    intent(in)  :: nurbs_nbasis, nurbs_order
integer :: ibasis
real(dp), dimension(size(x)) :: P1

B = 0.0d0
do ibasis = 1, nurbs_nbasis
    call bsplineP(x, ibasis, nurbs_order, P1)
    B = B + nurbs_weights(ibasis) * P1
end do

return
end subroutine nurbs_weighting
subroutine nurbs_weighting_grad(x, nurbs_nbasis, nurbs_order, B)
! Weighting function W(x) = SUM[ w_{i} * dN_{i,p} ]
implicit none
real(dp),   intent(out) :: B(:)
real(dp),   intent(in)  :: x(:)
integer,    intent(in)  :: nurbs_nbasis, nurbs_order
integer :: ibasis
real(dp), dimension(size(x)) :: P1

B = 0.0d0
do ibasis = 1, nurbs_nbasis
    call gradbsplineP(x, ibasis, nurbs_order, P1)
    B = B + nurbs_weights(ibasis) * P1
end do

return
end subroutine nurbs_weighting_grad

end module nurbs


