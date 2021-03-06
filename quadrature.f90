module quadrature
use prec
use glob, only: refa, refb

integer :: iquad, quad_npts
real(dp) :: quad_alpha = 0.0d0, quad_beta = 0.0d0
real(dp), allocatable :: xquad(:), wquad(:)

contains

subroutine initialize_quad
use glob, only: nref,select_node
quad_npts = nref
quad_npts = nref
allocate(xquad(quad_npts), wquad(quad_npts))
select case(select_node)
    case(11) ! Legendre-Gauss-Lobatto
        call legendreGLNodesWeights(quad_npts-1, xquad, wquad)
    case(21) ! Gauss-Legendre
        call JacobiGQ(xquad, wquad, quad_alpha, quad_beta, quad_npts-1)
    case default
        call JacobiGQ(xquad, wquad, quad_alpha, quad_beta, quad_npts-1)
        !call legendreGLNodesWeights(quad_npts-1, xquad, wquad)
end select
xquad = 0.5d0*((refb-refa)*xquad+(refa+refb))
end subroutine initialize_quad
subroutine finalize_quad
deallocate(xquad, wquad)
end subroutine finalize_quad

real(dp) function integrate(f)
implicit none
real(dp) :: f(quad_npts)
integrate = 0.5d0*(refb-refa)*dot_product(wquad, f)
end function integrate

subroutine JacobiGQ(ref, weights, alpha, beta, order)
! Gauss-Quadrature using Golub's tridiagonal eigendecomposition to recover the
! nodes and weights (ref and weights)
implicit none
real(dp),   intent(in)  :: alpha, beta
integer,    intent(in)  :: order
real(dp),   intent(out)  :: ref(order+1), weights(order+1)

real(dp) :: Jac(order+1, order+1)
real(dp) :: h, eps = 1.0e-10
integer ::i
integer :: lda, info, lwork
integer,parameter :: lwmax = 100
real(dp) :: eigval(order+1), eigvec(order+1, order+1), work(lwmax)

Jac = 0.0d0
if(order == 0) then
    ref(1) = -(alpha-beta)/(alpha+beta+2)
    weights(1) = 2.0d0
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
weights = (eigvec(1,:))**2.0d0 * 2.0d0**(alpha+beta+1.0d0) &
    / (alpha+beta+1) * gamma(alpha+1.0d0) * gamma(beta+1.0d0) / gamma(alpha+beta+1.0d0)

return

end subroutine JacobiGQ

end module quadrature
