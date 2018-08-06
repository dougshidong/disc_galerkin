module lagrange
! NOTE: Changing the typical interval of [0,1] to [-1,1] to conform with the
!       other bases. Done by making x = (x+1)/2
use prec
use glob
use quadrature
implicit none

contains
subroutine LagrangeP(x, k, n, P)
implicit none
real(dp),   intent(out) :: P(:)
real(dp),   intent(in)  :: x(:)
integer,    intent(in)  :: k, n
integer  j

P = 1.0d0
do j = 1, quad_npts
    if(j.ne.k) then
        P = P * (x - quadR(j)) / (quadR(k) - quadR(j))
    end if
end do

return
end subroutine LagrangeP

subroutine gradLagrangeP(x, k, n, P)
implicit none
real(dp),   intent(out) :: P(:)
real(dp),   intent(in)  :: x(:)
integer,    intent(in)  :: k, n
real(dp), dimension(size(x)) :: P1, P2
integer  j,m,i

P2 = 0.0d0
do i = 1, quad_npts
    if(i.ne.k) then
        P1 = 1.0d0
        do m = 1, quad_npts
            if(m.ne.k .and. m.ne.i) then
                P1 = P1 * (x - quadR(m)) / (quadR(k) - quadR(m))
            end if
        end do
        P2 = P2 + P1/(quadR(k) - quadR(i))
    end if
end do
P = P2
!
!
!call LagrangeP(x, k, n, P1)
!P2 = 0.0d0
!do j = 1, quad_npts
!    if(j.ne.k) then
!        P2 = P2 + 1.0d0/(x(k) - quadR(j))
!    end if
!end do
!
!P = P1*P2
return
end subroutine gradLagrangeP


end module lagrange

