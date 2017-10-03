module poly
use prec
use glob

real(dp) :: poly_alpha, poly_beta
!real(dp), allocatable :: V(:,:), invV(:,:), Dr(:,:)

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

subroutine JacobiGQ(ref, weights, alpha, beta, order)
! Gauss-Quadrature using Golub's tridiagonal eigendecomposition to recover the
! nodes and weights (ref and weights)
implicit none
real(dp),   intent(in)  :: alpha, beta
integer,    intent(in)  :: order
real(dp),   intent(out)  :: ref(order+1), weights(order+1)

real(dp) :: Jac(order+1, order+1)
integer :: i,j
real(dp) :: h, eps = 1.0e-10

integer :: lda, info, lwork
integer,parameter :: lwmax = 100
real(dp) :: eigval(order+1), eigvec(order+1, order+1), work(lwmax)

Jac = 0.0d0
if(order == 0) then
    ref(1) = -(alpha-beta)/(alpha+beta+2)
    weights(1) = 2
    return
end if

! Form tridiagonal symmetric matrix from recurrence
! Diagonal
do i = 1, order+1
    h = 2.0d0*(i-1) + alpha + beta
    Jac(i,i) = -(alpha**2.0d0-beta**2.0d0) / (h*(h+2.0d0))
end do
! Off-diagonal
do i = 1, order
    h = 2.0d0*(i-1) + alpha + beta
    Jac(i,i+1) = 2.0d0 / (h+2.0d0) &
        * sqrt( i * (i+alpha+beta) * (i+alpha) * (i+beta) / (h+1.0d0) / (h+3.0d0) )
    Jac(i+1,i) = Jac(i,i+1)
end do

if (alpha+beta<10*eps) Jac(1,1)=0.0d0

! Eigenvalue decomposition
eigvec = Jac
lda = order+1
lwork = -1
!do, i=1,order+1
!    write(*,"(100E15.5)") ( Jac(i,j), j=1,order+1 )
!enddo
! Query workspace size
call dsyev('V', 'U', order+1, eigvec, lda, eigval, &
            work, lwork, info)
lwork = min( lwmax, int( work( 1 ) ) )
! Evaluate eigenvalues and eigenvectors
call dsyev('V', 'U', order+1, eigvec, lda, eigval, &
            work, lwork, info )
! Check convergence
if( info.gt.0 ) then
    write(*,*) 'DSYEV ERROR IN GETEIG'
    stop
end if

! Compute quadrature by eigenvalue solve
ref = eigval
weights = (eigvec(1,:))**2.0d0 * 2**(alpha+beta+1.0d0) &
    / (alpha+beta+1) * gamma(alpha+1) * gamma(beta+1) / gamma(alpha+beta+1)

return

end subroutine JacobiGQ

subroutine gradJacobiP(r, alpha, beta, order, Pd)
implicit none
integer, intent(in) :: order
real(dp), intent(in) :: r(:), alpha, beta
real(dp), intent(out) :: Pd(size(r))
real(dp) :: P(size(r))

if(order.eq.0) then
    Pd = 0.0d0
else
    call JacobiP(r, alpha+1, beta+1, order-1, P)
    Pd = sqrt(order*(order+alpha+beta+1))*P
end if

return
end subroutine gradJacobiP

subroutine modalToNodal(r, alpha, beta, order, V, invV)
implicit none
integer, intent(in) :: order
real(dp), intent(in) :: r(order+1), alpha, beta
real(dp), intent(out) :: V(order+1, order+1)
real(dp), intent(out) :: invV(order+1, order+1)
real(dp) :: tempV(order+1, order+1)
real(dp) :: P(order+1)
integer :: i, j
integer :: info, ipiv(order+1)

do j = 0, order
    call JacobiP(r, alpha, beta, j, P)
    V(:,j+1) = P
end do
tempV = V
invV = 0.0d0
! Create V inverse
do i = 1, order+1
    invV(i,i) = 1.0d0
end do
call DGESV(order+1, order+1, tempV, order+1, ipiv, invV, order+1, info )
return
end subroutine modalToNodal

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
