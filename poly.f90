module poly
use prec
use glob
use jacobi
use monomial
use bezier
use bspline

real(dp) :: poly_alpha, poly_beta

contains

subroutine poly_eval(polytype, nth, polyx, polyVal) 
implicit none
integer :: polytype, nth
real(dp) :: polyx(:), polyVal(:)

!nth = 0, ..., nbasis
select case(polytype)
    case(1)
        call binomial_eval
        call BezierP(polyx, nth, order, polyVal)
    case(2)
        call bspline_init
        call bsplineP(polyx, nth+1, bs_order, polyVal)
    case(8)
        call monomialP(polyx, nth, polyVal)
    case default
        call JacobiP(polyx, jacobi_alpha, jacobi_beta, nth, polyVal)
end select
end subroutine poly_eval

subroutine polyGrad_eval(polytype, nth, polyx, polyGrad) 
implicit none
integer :: polytype, nth
real(dp) :: polyx(:), polyGrad(:)

select case(polytype)
    case(1)
        call binomial_eval
        call gradBezierP(polyx, nth, order, polyGrad)
    case(2)
        call bspline_init
        call gradbsplineP(polyx, nth+1, bs_order, polyGrad)
    case(8)
        call gradMonomialP(polyx, nth, polyGrad)
    case default
        call gradJacobiP(polyx, jacobi_alpha, jacobi_beta, nth, polyGrad)
end select
end subroutine polyGrad_eval


subroutine gradModalToNodal(r, alpha, beta, order, Vr)
implicit none
integer, intent(in) :: order
real(dp), intent(in) :: r(order+1), alpha, beta
real(dp), intent(out) :: Vr(order+1, order+1)
real(dp) :: Pd(order+1)
integer :: j

do j = 0, order
    call gradJacobiP(r, alpha, beta, j, Pd)
    Vr(:,j+1) = Pd
end do
return
end subroutine gradModalToNodal

subroutine massToStiff(r, alpha, beta, order, Vinv, Dr)
implicit none
integer, intent(in) :: order
real(dp), intent(in) :: r(order+1), alpha, beta
real(dp), intent(in) :: Vinv(order+1, order+1)
real(dp), intent(out) :: Dr(order+1, order+1)
real(dp) :: Vr(order+1, order+1)

! Evaluate Vr, but store in Dr
call gradModalToNodal(r, alpha, beta, order, Vr)
! Evaluate Vtran Drtran = Vrtran
Dr = matmul(Vr,Vinv)
return
end subroutine massToStiff


end module poly
