module poly
use prec
use glob
use jacobi
use monomial
use bezier
use bspline
use nurbs
use lagrange

real(dp) :: poly_alpha, poly_beta

contains

subroutine poly_eval(polytype, ibasis, order, ndim, npts, polyx, polyVal) 
implicit none
integer :: order, ndim, npts
integer :: polytype, ibasis, ib, jb
real(dp) :: polyx(npts, ndim), polyVal(npts)
real(dp) :: polytemp(npts, ndim)

ib = mod(ibasis-1,order+1) + 1
!print*, polyx(2,:)
call poly1D_eval(polytype, ib, polyx(:, 1), polytemp(:, 1))
polyVal = polytemp(:, 1)

! Quads
if(ndim.eq.2) then
    jb = (ibasis-1)/(order+1) + 1
    call poly1D_eval(polytype, jb, polyx(:, 2), polytemp(:, 2))
    polyVal = polytemp(:, 1)*polytemp(:, 2)
end if

return
end subroutine

subroutine polyGrad_eval(polytype, ibasis, order, direction, ndim, npts, polyx, polygrad) 
implicit none
integer :: order, direction, ndim, npts
integer :: polytype, ibasis, ib, jb
real(dp) :: polyx(npts, ndim), polygrad(npts)
real(dp) :: polytemp(npts, ndim)

if(ndim.eq.1) then
    ib = mod(ibasis-1,order+1) + 1
    call poly1DGrad_eval(polytype, ibasis, polyx(:,1), polytemp(:,1)) 
    polygrad = polytemp(:,1)
else if(ndim.eq.2) then
    ib = mod(ibasis-1,order+1) + 1
    jb = (ibasis-1)/(order+1) + 1
    if(direction.eq.1) then
        call poly1DGrad_eval(polytype, ib, polyx(:,1), polytemp(:,1)) 
        call poly1D_eval(polytype, jb, polyx(:,2), polytemp(:,2))
    else if(direction.eq.2) then
        call poly1D_eval(polytype, ib, polyx(:,1), polytemp(:,1))
        call poly1DGrad_eval(polytype, jb, polyx(:,2), polytemp(:,2)) 
    end if
    polygrad = polytemp(:,1)*polytemp(:,2)
end if

return
        
end subroutine

subroutine poly1D_eval(polytype, ibasis, polyx, polyVal) 
implicit none
integer :: polytype, ibasis
real(dp) :: polyx(:), polyVal(:)

!ibasis = 1, ..., nbasis
select case(polytype)
    case(0)
        call LagrangeP(polyx, ibasis, order, polyVal)
    case(1)
        call binomial_eval
        call BezierP(polyx, ibasis-1, order, polyVal)
    case(2)
        call bspline_init
        call bsplineP(polyx, ibasis, bs_order, polyVal)
    case(3)
        call nurbs_init
        call nurbsP(polyx, ibasis, bs_order, polyVal)
    case(8)
        call monomialP(polyx, ibasis-1, polyVal)
    case default
        call JacobiP(polyx, jacobi_alpha, jacobi_beta, ibasis-1, polyVal)
end select
end subroutine poly1D_eval

subroutine poly1DGrad_eval(polytype, ibasis, polyx, polygrad) 
implicit none
integer :: polytype, ibasis
real(dp) :: polyx(:), polyGrad(:)

!ibasis = 1, ..., nbasis
select case(polytype)
    case(0)
        call gradLagrangeP(polyx, ibasis, order, polyGrad)
    case(1)
        call binomial_eval
        call gradBezierP(polyx, ibasis-1, order, polyGrad)
    case(2)
        call bspline_init
        call gradbsplineP(polyx, ibasis, bs_order, polyGrad)
    case(3)
        call nurbs_init
        call gradnurbsP(polyx, ibasis, bs_order, polyGrad)
    case(8)
        call gradMonomialP(polyx, ibasis-1, polyGrad)
    case default
        call gradJacobiP(polyx, jacobi_alpha, jacobi_beta, ibasis-1, polyGrad)
end select
end subroutine poly1DGrad_eval


!subroutine gradModalToNodal(r, alpha, beta, order, Vr)
!implicit none
!integer, intent(in) :: order
!real(dp), intent(in) :: r(order+1), alpha, beta
!real(dp), intent(out) :: Vr(order+1, order+1)
!real(dp) :: Pd(order+1)
!integer :: j
!
!do j = 0, order
!    call gradJacobiP(r, alpha, beta, j, Pd)
!    Vr(:,j+1) = Pd
!end do
!return
!end subroutine gradModalToNodal
!
!subroutine massToStiff(r, alpha, beta, order, Vinv, Dr)
!implicit none
!integer, intent(in) :: order
!real(dp), intent(in) :: r(order+1), alpha, beta
!real(dp), intent(in) :: Vinv(order+1, order+1)
!real(dp), intent(out) :: Dr(order+1, order+1)
!real(dp) :: Vr(order+1, order+1)
!
!! Evaluate Vr, but store in Dr
!call gradModalToNodal(r, alpha, beta, order, Vr)
!! Evaluate Vtran Drtran = Vrtran
!Dr = matmul(Vr,Vinv)
!return
!end subroutine massToStiff


end module poly
