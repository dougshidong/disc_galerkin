subroutine lift1d(V, lift)
! Used to extract surface terms
! (uh - u*) on right face - (uh - u*) on left face
! Minv * [l(x) * term]^xr_xl
use prec
use glob
implicit none

real(dp), intent(in) :: V(nref, nref)
real(dp), intent(out) :: lift(nref, nface*nfpoint)
real(dp) :: Emat(nref, nface*nfpoint)
Emat = 0.0d0
! If LGL points are used, 1s at diagonal corners, and 0s elsewhere
!Emat(:,1) = L1
!Emat(:,2) = L2

! Minv * E
lift = matmul( V, matmul( transpose(V), Emat ) )

return

end subroutine lift1d
