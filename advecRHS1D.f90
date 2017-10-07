subroutine advecRHS1D(u, time, wavespeed, fscheme, rhsu)
use prec
use glob
!use poly, only: V, Dr
use matrices
use grid, only: rx, normals, Fscale
implicit none

real(dp), intent(in) :: u(nref,nele)
real(dp), intent(in) :: time, wavespeed, fscheme
real(dp), intent(out) :: rhsu(nref,nele)
real(dp) :: uface(2,nele)
real(dp) :: du(nface*nfpoint,nele)
real(dp) :: flux_numerical(nvertex)
real(dp) :: uin, u1, u2, n1, n2
integer :: iele, iface

uface = 0.0d0
do iele = 1, nele
    uface(1,iele) = dot_product(u(:,iele), L1(:))
    uface(2,iele) = dot_product(u(:,iele), L2(:))
end do
uface = matmul(transpose(VanderSurf), u)

! Inlet
select case(icase)
    case(0)
        uin = -sin(wavespeed*time)
    case(1)
        uin = 0.0d0
    case default
        print *, 'Invalid case'
end select

!call flux_numerical_eval(uface, wavespeed, fscheme, flux_numerical)
flux_numerical = 0.0
do iface = 2, nele
    iele = iface-1
    u1 = uface(2,iele)
    u2 = uface(1,iele+1)
    n1 = normals(2,iele)
    n2 = normals(1,iele+1)
    flux_numerical(iface) = flux_eval(u1, u2, n1, n2, wavespeed, fscheme)
end do
! Field differences at faces
do iele = 1, nele
    du(1,iele) = flux_analytical(uface(1,iele),wavespeed) - flux_numerical(iele)
    du(2,iele) = flux_analytical(uface(2,iele),wavespeed) - flux_numerical(iele+1)
    !du(2,iele) = wavespeed * uface(2,iele) - flux_numerical(iele+1)
end do
u1 = uin
n1 = 1.0d0
u2 = uface(1,1)
n2 = normals(1,1)
du(1,1) = flux_analytical(uface(1,1),wavespeed) - flux_eval(u1, u2, n1, n2, wavespeed, fscheme)
! du_right - du_left
du(1,:) = -du(1,:)
!du(1,1) = 0.5d0*(uface(1,1) - uin) * (-wavespeed - (1.0d0-fscheme) * wavespeed)
!du(1,1) = 0.5d0*wavespeed*(u1 - u2) + (1.0d0-fscheme) * wavespeed * 0.5d0 * (u1*n1 + u2*n2)

! Outlet
du(2,nele) = 0.0d0

! Evaluate RHS of semi-discrete PDE
! dudt = -wavespeed * Minv * S * u + Minv[ l(x) * du ]^xr_xl
rhsu = -wavespeed * rx * matmul(Differentiation,u)
rhsu = rhsu + matmul(Lift, Fscale*du)

return

contains

subroutine flux_numerical_eval(uface, wavespeed, fscheme, flux_numerical)
use prec
use glob
use grid, only: normals
implicit none

real(dp), intent(in) :: uface(2,nele)
real(dp), intent(in) :: wavespeed, fscheme
real(dp), intent(out) :: flux_numerical(nele+1)
integer :: iface, iele

!! Upwind flux
flux_numerical = 0.0
do iface = 2, nele
    iele = iface-1
    flux_numerical(iface) = uface(2,iele) + uface(1,iele+1) 
    flux_numerical(iface) = flux_numerical(iface) + (1.0-fscheme) &
        * (uface(2,iele) * normals(2,iele) + uface(1,iele+1)*normals(1,iele+1))
end do
flux_numerical = flux_numerical * wavespeed * 0.5d0

! Lax-Friedrichs
do iface = 2, nele
    iele = iface-1
    flux_numerical(iface) = 0.5d0*wavespeed*(uface(2,iele) + uface(1,iele+1))
    flux_numerical(iface) = flux_numerical(iface) + wavespeed * 0.5d0 * &
        (uface(2,iele) * normals(2,iele) + uface(1,iele+1)*normals(1,iele+1))
end do

do iface = 2, nele
    iele = iface-1
    flux_numerical(iface) = 0.5d0*wavespeed*(uface(2,iele) + uface(1,iele+1))
    flux_numerical(iface) = flux_numerical(iface) + wavespeed * 0.5d0 * &
        (uface(2,iele) * normals(2,iele) + uface(1,iele+1)*normals(1,iele+1))
end do

return
end subroutine flux_numerical_eval
real(dp) function flux_eval(u1,u2,n1,n2,A,fscheme)
    implicit none
    real(dp), intent(in) :: u1, u2, n1, n2, A, fscheme
    flux_eval = 0.5d0*A*(u2 + u1) + (1.0d0-fscheme) * A * 0.5d0 * (u1*n1 + u2*n2)
end function flux_eval

real(dp) function flux_analytical(u, wavespeed)
use prec
use glob
implicit none

real(dp), intent(in) :: u, wavespeed

select case(iequation)
    case(1)
        flux_analytical = wavespeed*u
    case default
end select

end function flux_analytical

end subroutine advecRHS1D

