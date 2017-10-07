module bezier
use prec
use glob
implicit none

integer :: bezier_order, bezier_npts
real(dp), allocatable :: bin_array(:,:)
logical :: bezier_init = .false.
contains
subroutine BezierP(x, k, n, P)
implicit none
real(dp),   intent(out) :: P(:)
real(dp),   intent(in)  :: x(:)
integer,    intent(in)  :: k, n

P = 0.0d0
if(k.lt.0 .or. n.lt.0 .or. k.gt.n) return
P = bin_array(n+(1),k+(1)) * (1.0d0-x)**(n-k) * x**(k)
return
end subroutine BezierP

subroutine gradBezierP(x, k, n, P)
implicit none
real(dp),   intent(out) :: P(:)
real(dp),   intent(in)  :: x(:)
integer,    intent(in)  :: k, n
real(dp), dimension(size(x)) :: P1, P2
call BezierP(x, k-1, n-1, P1)
call BezierP(x, k, n-1, P2)
P = n*(P1-P2)
!P = bin_array(n-1+(1),k+(1)) * (1.0d0-x)**(n-1-k) * x**(k) * n
return
end subroutine gradBezierP

subroutine binomial_eval
    integer :: i, j
    if(bezier_init) return
    bezier_order = order
    bezier_npts  = order + 1
    bezier_init = .true.
    allocate(bin_array(bezier_npts, bezier_npts))
    bin_array = 0.0d0
    ! Note, Fortran array index start at 1.
    ! Index 1 is associated with 0th poly/order
    bin_array(:,1) = 1.0d0
    do i = 2, bezier_npts
    do j = 2, i
        bin_array(i,j) = bin_array(i-1,j-1) + bin_array(i-1,j) 
    end do
    end do
    return
end subroutine binomial_eval



end module bezier
