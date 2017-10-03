subroutine advecRHS1D_G(u, time, wavespeed, alpha, rhsu)
use prec
use glob
!use poly, only: V, Dr
use matrices
use grid, only: rx, normals, Fscale
implicit none

real(dp), intent(in) :: u(nref,nele)
real(dp), intent(in) :: time, wavespeed, alpha
real(dp), intent(out) :: rhsu(nref,nele)
real(dp) :: uface(2,nele)
real(dp) :: du(nface*nfpoint,nele)
real(dp) :: f(nele+1)
real(dp) :: uin
integer :: iele, iref

uface = 0.0d0
do iele = 1, nele
    uface(1,iele) = dot_product(u(:,iele), L1(:))
    uface(2,iele) = dot_product(u(:,iele), L2(:))
end do

call flux_G(uface, wavespeed, alpha, f)
! Boundary conditions
! Inlet
uin = 0.0d0
!uin = -sin(wavespeed*time)
f(1) = uin + uface(1,1)

f(1) = 0.5d0*wavespeed*(uin + uface(1,1))
f(1) = f(1) + wavespeed * 0.5d0 *(uin + uface(1,1)*normals(1,1))
! Field differences at faces
do iele = 1, nele
    du(1,iele) = wavespeed * uface(1,iele) - f(iele)
    du(2,iele) = wavespeed * uface(2,iele) - f(iele+1)
end do
du(1,:) = -du(1,:)
du(1,1) = 0.5d0*(uface(1,1) - uin) &
&   * (-wavespeed - (1.0d0-alpha) * wavespeed)
! Outlet
du(2,nele) = 0.0d0

! Evaluate RHS of semi-discrete PDE
! dudt = -wavespeed * Minv * S * u + Minv[ l(x) * du ]^xr_xl
rhsu = -wavespeed * rx * matmul(Differentiation,u)
rhsu = rhsu + matmul(Lift, Fscale*du)

return

end subroutine advecRHS1D_G

subroutine flux_G(uface, wavespeed, alpha, f)
use prec
use glob
use grid, only: normals
implicit none

real(dp), intent(in) :: uface(2,nele)
real(dp), intent(in) :: wavespeed, alpha
real(dp), intent(out) :: f(nele+1)
integer :: iface, iele

!! Upwind flux
f = 0.0
do iface = 2, nele
    iele = iface-1
    f(iface) = uface(2,iele) + uface(1,iele+1) 
    f(iface) = f(iface) + (1.0-alpha) &
        * (uface(2,iele) * normals(2,iele) + uface(1,iele+1)*normals(1,iele+1))
end do
f = f * wavespeed * 0.5d0

! Lax-Friedrichs
do iface = 2, nele
    iele = iface-1
    f(iface) = 0.5d0*wavespeed*(uface(2,iele) + uface(1,iele+1))
    f(iface) = f(iface) + wavespeed * 0.5d0 * &
        (uface(2,iele) * normals(2,iele) + uface(1,iele+1)*normals(1,iele+1))
end do

return
end subroutine flux_G

