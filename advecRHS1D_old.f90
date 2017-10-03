subroutine advecRHS1D(u, time, wavespeed, alpha, rhsu)
use prec
use glob
use poly, only: V, Dr
use grid, only: rx, normals, Fscale
implicit none

real(dp), intent(in) :: u(nref,nele)
real(dp), intent(in) :: time, wavespeed, alpha
real(dp), intent(out) :: rhsu(nref,nele)
real(dp) :: du(nface*nfpoint,nele)
real(dp) :: f(nele+1)
real(dp) :: lift(nref, nface*nfpoint), uin
integer :: iele

call flux(u, wavespeed, alpha, f)
! Boundary conditions
! Inlet
!uin = -sin(wavespeed*time) + -sin(100*wavespeed*time)
!uin = -sin(wavespeed*time)
uin = 0.0d0
!uin = -sin(wavespeed*time)
f(1) = uin + u(1,1)
f(1) = f(1) + (1.0d0-alpha) * (uin + u(1,1)*normals(1,1))
f(1) = f(1) * wavespeed * 0.5d0

f(1) = 0.5d0*wavespeed*(uin + u(1,1))
f(1) = f(1) + wavespeed * 0.5d0 *(uin + u(1,1)*normals(1,1))
!! Field differences at faces
do iele = 1, nele
    du(1,iele) = wavespeed * u(1,iele) - f(iele)
    du(2,iele) = wavespeed * u(nref,iele) - f(iele+1)
end do
du(1,:) = -du(1,:)
du(1,1) = 0.5d0*(u(1,1) - uin) &
&   * (-wavespeed - (1.0d0-alpha) * wavespeed)
! Outlet
du(2,nele) = 0.0d0

call lift1d(V, lift)
! Evaluate RHS of semi-discrete PDE
! dudt = -wavespeed * Minv * S * u + Minv[ l(x) * du ]^xr_xl
rhsu = -wavespeed * rx * matmul(Dr,u)
rhsu = rhsu + matmul(lift, Fscale*du)

!print *, 'rhsu'
!call printMatrix(rhsu)
!print *, 'lift'
!call printMatrix(lift)
!print *, 'Fscale'
!call printMatrix(Fscale)
!print *, 'du'
!call printMatrix(du)

return

end subroutine advecRHS1D

subroutine flux(u, wavespeed, alpha, f)
use prec
use glob
use grid, only: normals
implicit none

real(dp), intent(in) :: u(nref,nele)
real(dp), intent(in) :: wavespeed, alpha
real(dp), intent(out) :: f(nele+1)
integer :: iface, iele

!! Upwind flux
f = 0.0
do iface = 2, nele
    iele = iface-1
    f(iface) = u(nref,iele) + u(1,iele+1) 
    f(iface) = f(iface) + (1.0-alpha) &
        * (u(nref,iele) * normals(2,iele) + u(1,iele+1)*normals(1,iele+1))
end do
f = f * wavespeed * 0.5d0

! Lax-Friedrichs
do iface = 2, nele
    iele = iface-1
    f(iface) = 0.5d0*wavespeed*(u(nref,iele) + u(1,iele+1))
    f(iface) = f(iface) + wavespeed * 0.5d0 * &
        (u(nref,iele) * normals(2,iele) + u(1,iele+1)*normals(1,iele+1))
end do

return
end subroutine flux
