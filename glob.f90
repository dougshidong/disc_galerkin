module glob
use prec
implicit none

integer :: nbasis
integer :: nele, nvertex

integer :: order

real(dp), parameter :: PI=4.D0*DATAN(1.D0)
!logical :: inodal = .true.
!logical :: inodal = (.not. .true.)
logical :: inodal = 0==0
! 0 = sine
! 1 = Bump
! 2 = x/abs(x)
! 3 = Differentiation of sine
! 4 = With source term = -sin(x)
integer :: icase = 0
integer :: iequation = 1
integer :: swform = 0
integer, parameter :: ndim = 2
integer, parameter :: nstate = 1
real(dp) :: finalTime = 2.5d0
real(dp), dimension(2) :: wavespeed = (/ 1.0d0, 1.0d0 /)
!real(dp), dimension(2) :: wavespeed = (/ 0.5d0, 0.5d0 /)
! 1 = Evaluate eig[semi-discrete]
! 2 = Evaluate eig[semi-discrete] for all wavenumbers
! 3 = Find maximum dt for various RK
integer :: istab = -1
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
!real(dp), dimension(ndim) :: refa = -1.0d0
!real(dp), dimension(ndim) :: refb = 1.0d0
real(dp) :: refa = -1.0d0
real(dp) :: refb = 1.0d0
!real(dp) :: refa = 0.0d0
!real(dp) :: refb = 1.0d0

! Flow Solver Settings
integer :: tscheme = 1 ! Time-stepping scheme
real(dp) :: fscheme = 0.0d0 ! Numerical flux scheme

contains 

subroutine printMatrix(A)
real(dp) :: A(:,:)
integer :: i, j
do, i=1,size(A,1)
do, j=1,size(A,2)
    if(abs(A(i,j)).le.1e-13) A(i,j) = 0.0d0
end do
end do
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
