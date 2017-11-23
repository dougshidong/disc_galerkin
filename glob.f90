module glob
use prec
implicit none

integer :: order
integer :: nele, nvertex, nref, nface, nfpoint
real(dp), parameter :: PI=4.D0*DATAN(1.D0)
!logical :: inodal = .true.
logical :: inodal = (.not. .true.)
integer :: iequation = 1
! 0 = sine
! 1 = Bump
! 2 = x/abs(x)
! 3 = Differentiation of sine
! 4 = Plot B-Splines
integer :: icase = 0
! 1 = Evaluate eig[semi-discrete]
! 2 = Evaluate eig[semi-discrete] for all wavenumbers
! 3 = Find maximum dt for various RK
integer :: istab = 1
! 1 = Bezier
! 2 = BSpline
! 3 = NURBS
! 8 = Monomial
! Default = Legendre
integer :: polytype = 9
! Node distribution selection
! 1 = Legendre-Gauss-Lobatto
! 2 = Gauss-Legendre
! 3 = Uniform distribution
integer :: select_node = 2
! Reference element Coordinates
real(dp) :: refa = -1.0d0
real(dp) :: refb = 1.0d0

contains 

subroutine printMatrix(A)
real(dp) :: A(:,:)
integer :: i, j
do, i=1,size(A,1)
    write(*,"(100E15.5)") ( A(i,j), j=1,size(A,2) )
enddo
return
end subroutine printMatrix
subroutine printVec(A)
real(dp) :: A(:)
integer :: j
do j = 1, size(A)
write(*,"(E15.5)") A(j)
end do
return
end subroutine printVec

end module glob
