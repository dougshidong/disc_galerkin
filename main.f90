program adv

use prec
use grid
use glob
use poly
use matrices
use quadrature
use stability

implicit none

real(dp) :: xmin, xmax
real(dp), allocatable :: ref(:), weights(:)
real(dp), allocatable :: temp_nref(:)
integer :: tscheme
real(dp) :: fscheme, wavespeed
integer :: i, j, k, l

real(dp), allocatable :: u(:,:)
real(dp), allocatable :: xsol(:,:), usol(:,:)
real(dp), allocatable :: VanderQuad(:,:)
real(dp), allocatable :: polyVal(:)
real(dp) :: finalTime

real(dp) :: wavenumber ! For stability analysis

order   = 6
nele    = 50
open(unit=1,file="input.in",form='formatted')
read(1,*) order, nele
close(1)
tscheme = 1
fscheme = 0.0d0
!fscheme = 1.0d0
wavespeed = 1.0d0
finalTime = 0.50d0

nref    = order + 1
nvertex = nele + 1
nface   = 2
nfpoint = 1

! Setup mass, Vandermonde, differentiation and lift matrices
call allocate_matrices

! Setup Quadrature
call initialize_quad
! Select points in reference element
allocate(ref(nref), weights(nref))
allocate(temp_nref(nref))
weights = 1.0d0
poly_alpha = 0.0d0
poly_beta = 0.0d0
! Node Selection
select case(select_node)
    case(1) ! Legendre-Gauss-Lobatto
        call legendreGLNodesWeights(order, ref, weights)
    case(2) ! Gauss-Legendre
        call JacobiGQ(ref, weights, poly_alpha, poly_beta, order)
    case(3) ! Uniform distribution
        do i = 1, nref
           ref(i) = refa + (i-1)*(refb-refa)/nref
        end do
    case default
        print*, 'Invalid node distribution. Using GL.'
        call JacobiGQ(ref, weights, poly_alpha, poly_beta, order)
end select
! Node distribution normalization to reference element
ref = 0.5d0*((refb-refa)*ref+(refa+refb))

! Evaluate transformation matrices
call buildVandermonde(ref)
call buildMass
call buildStiffness
call buildDifferentiation
call buildLift

j = 10
do i = 0, j
wavenumber = 2.0d0*PI/j * i - PI
call check_stability(wavenumber)
end do
stop
!call printmatrix(Vander)
!print*,
!call printmatrix(Mass)
!print*,
!call printmatrix(Stiff)
!print*,
!call printmatrix(Differentiation)
!print*,
!call printmatrix(Lift)

!call printmatrix(matmul(Mass,MassInv))
!print*,
!call printmatrix(matmul(Vander,VanderInv))
!print*,

! Old subroutines replaced by more general formulations
!call modalToNodal(ref, poly_alpha, poly_beta, order, Vandermonde, VandermondeInv)
!call massToStiff(ref, poly_alpha, poly_beta, order, VanderInv, Differentiation)
!call lift1d(Vander, Lift)

call genGrid(ref, xmin, xmax)

allocate(u(nref,nele))
call initialize_u(u, x)

! umodal = VanderInv * u
allocate(xsol(quad_npts,nele))
allocate(usol(quad_npts,nele))
allocate(VanderQuad(quad_npts, nref))
allocate(polyVal(quad_npts))
! Evaluate Vandermonde
do j = 0, order
    call poly_eval(polytype, j, xquad, polyVal)
    VanderQuad(1:quad_npts,j+1) = polyVal
end do
u = matmul(VanderInv,u)
usol = matmul(VanderQuad,u)
u = matmul(Vander,u)
do i = 1, quad_npts
do j = 1, nele
    xsol(i,j) = vertices(EtoV(1,j)) &
        + 1.0d0/(refb-refa) &!xquad(quad_npts)-xquad(1)) &
        * (xquad(i)-refa) * (vertices(EtoV(2,j)) - vertices(EtoV(1,j)))
end do
end do

open(unit=7, file='output.dat', form='formatted')
! do j = 1, nele
! do i = 1, nref
!  if(abs(u(i,j)).le.1e-30) u(i,j) = 0.0d0
!  write(7,11) x(i,j), u(i,j)
! end do
! end do
do j = 1, nele
do i = 1, quad_npts
 if(abs(usol(i,j)).le.1e-30) usol(i,j) = 0.0d0
 write(7,11) xsol(i,j), usol(i,j)
end do
end do

if(icase.lt.3) then
    call advec1D(u, tscheme, fscheme, wavespeed, finalTime)
    ! umodal = VanderInv * u
    u = matmul(VanderInv,u)
    usol = matmul(VanderQuad,u)
    u = matmul(Vander,u)
else if(icase.eq.3) then
    u = rx*matmul(Differentiation,u)
    u = matmul(VanderInv,u)
    usol = matmul(VanderQuad,u)
else if(icase.eq.4) then
    call bspline_init
    do k = 1, bs_nbasis
        do l = 1, nele
        call bsplineP(x(:,l), k, bs_order, u(:,l))
        end do
        u = matmul(VanderInv,u)
        usol = matmul(VanderQuad,u)
        u = matmul(Vander,u)
        do j = 1, nele
        do i = 1, quad_npts
         if(abs(usol(i,j)).le.1e-30) usol(i,j) = 0.0d0
         write(7,11) xsol(i,j), usol(i,j)
        end do
        end do
    end do
end if

do j = 1, nele
do i = 1, quad_npts
 if(abs(usol(i,j)).le.1e-30) usol(i,j) = 0.0d0
 write(7,11) xsol(i,j), usol(i,j)
end do
end do
do i = 1, quad_npts
 write(7,11) xquad(i), wquad(i)
end do

!do j = 1, nele
!do i = 1, nref
! if(abs(u(i,j)).le.1e-30) u(i,j) = 0.0d0
! write(7,11) x(i,j), u(i,j)
!end do
!end do
!! weights
!do i = 1, nref
! write(7,11) xquad(i), wquad(i)
!end do



close(7)


! Deallocate matrices
call deallocate_matrices
call finalize_quad
call finalize_grid

11 FORMAT(' ',4E24.15)
return
contains

subroutine advec1D(u, tscheme, fscheme, wavespeed, finalTime)
use prec
use glob
use grid
implicit none
integer, intent(in) :: tscheme
real(dp), intent(in) :: finalTime, wavespeed, fscheme
real(dp), intent(inout) :: u(:,:)

integer :: tstep, nstep, intrk
real(dp) :: timelocal, time, dxmin, CFL, dt
real(dp), allocatable :: resiu(:,:), rhsu(:,:)
real(dp) :: rk4a(5), rk4b(5), rk4c(5)

allocate(resiu(nref,nele))
allocate(rhsu(nref,nele))

CFL = 0.001d0
time = 0.0d0
resiu = 0.0d0

dxmin = 1e9
do tstep = 1, nref-1
    timelocal = minval(abs(x(tstep,:)-x(tstep+1,:)))
    dxmin = min(timelocal,dxmin)
end do
dt = 0.5d0 * CFL / wavespeed * dxmin
!dt = 1.0000000000000001e-5

nstep = ceiling(finalTime/dt)
dt = finalTime/nstep

rk4a(1:5) = (/  0.0d0, &
                -567301805773.0d0/1357537059087.0d0, &
                -2404267990393.0d0/2016746695238.0d0, &
                -3550918686646.0d0/2091501179385.0d0, &
                -1275806237668.0d0/842570457699.0d0 /)
rk4b(1:5) = (/  1432997174477.0d0/9575080441755.0d0, &
                5161836677717.0d0/13612068292357.0d0, &
                1720146321549.0d0/2090206949498.0d0, &
                3134564353537.0d0/4481467310338.0d0, &
                2277821191437.0d0/14882151754819.0d0 /)
rk4c(1:5) = (/  0.0d0, &
                1432997174477.0d0/9575080441755.0d0, &
                2526269341429.0d0/6820363962896.0d0, &
                2006345519317.0d0/3224310063776.0d0, &
                2802321613138.0d0/2924317926251.0d0 /)

do tstep = 1, nstep
    if(tscheme.eq.1) then
        do intrk = 1, 5
            timelocal = time + rk4c(intrk)*dt
            call advecRHS1D(u, timelocal, wavespeed, fscheme, rhsu)
            resiu = rk4a(intrk)*resiu + dt*rhsu
            u = u + rk4b(intrk)*resiu
        end do
    else if(tscheme.eq.0) then
        call advecRHS1D(u, time, wavespeed, fscheme, rhsu)
        u = u + dt*rhsu
    end if
    time = time + dt
end do

deallocate(rhsu, resiu)

return

end subroutine advec1D

subroutine initialize_u(u, x)
use glob
use prec
implicit none
real(dp), intent(in) :: x(nref,nele)
real(dp), intent(out) :: u(:,:)
integer :: i, j
u = 0.0d0
select case(icase)
    case(0)
        u = sin(x)
    case(1)
        do j = 1, nele
        do i = 1, nref
         if(abs(x(i,j)).lt.1.0) then
            u(i,j) = exp(-1.0d0/(1.0d0-x(i,j)**2))
         end if
        end do
        end do
    case(2)
        u = x/abs(x)
    case default
        print *, 'Invalid case'
end select

end subroutine initialize_u


end program adv
