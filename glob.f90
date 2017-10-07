module glob
use prec
implicit none

integer :: order
integer :: nele, nvertex, nref, nface, nfpoint
real(dp), parameter :: PI=4.D0*DATAN(1.D0)
real(dp),allocatable :: L1(:), L2(:)
integer :: iequation = 1, icase
integer :: polytype = 999

contains 

subroutine printMatrix(A)

real(dp) :: A(:,:)
integer :: i, j
do, i=1,size(A,1)
    write(*,"(100E15.5)") ( A(i,j), j=1,size(A,2) )
enddo
return

end subroutine printMatrix

end module glob
