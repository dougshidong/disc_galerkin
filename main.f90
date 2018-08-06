program adv

use prec
use glob
use grid
!use poly
use matrices
use quadrature
use cubature2D
use ref_ele_mod
use testcases
use rhs
!use stability

implicit none

integer :: i, j, k, l
integer :: iele, idir, istate

real(dp) :: output_var(1000000,5)
real(dp), allocatable :: utemp(:,:)
real(dp) :: error_abs, error_rel

open(unit=1,file="input.in",form='formatted')
read(1,*) order, nivertex, njvertex
close(1)

if(ndim.eq.1) nbasis = order + 1
if(ndim.eq.2) nbasis = (order + 1)**2

! Setup Quadrature
call initialize_quad
call initialize_cub2D
do i = 1, quad_npts
    l = 0
    output_var(i, 1) = quadR(i)
    output_var(i, 2) = quadW(i)
end do
call output(quad_npts, 2, 'quad.x')
do i = 1, cub2D_npts
    l = 0
    output_var(i, 1) = cubRS(i,1)
    output_var(i, 2) = cubRS(i,2)
    output_var(i, 3) = cubW(i)
end do
call output(cub2D_npts, 3, 'cub2D.x')

! Generate reference element, grid and connectivity
call initialize_reference_elements

! Setup mass, Vandermonde, Mass, Differentiation and Lift matrices
call allocate_matrices
call buildVandermonde
call buildVandermondeGrad
call buildMass
call buildStiffness
call buildDifferentiation
call buildLift
!!!
call genGrid2D

call initialize_u
! umodal = VanderInv * u
do i = 1, nele
    if(.not. inodal) elements(i)%u = matmul(VanderInv,elements(i)%u)
end do

k = 0
do i = 1, nele
    if(.not. inodal) elements(i)%u = matmul(Vander,elements(i)%u)
    do j = 1, ref_ele(elements(i)%typ)%p%nnode_cub
        l = 0
        k = k + 1
        do idir = 1, ndim
            l = l+1
            output_var(k,l) = elements(i)%x_cub(j,idir)
        end do
        do istate = 1, nstate
            l = l+1
            output_var(k,l) = elements(i)%u(j,1)
        end do
    end do
    if(.not. inodal) elements(i)%u = matmul(VanderInv,elements(i)%u)
end do
call output(k, ndim+1, 'init.x')

if(icase.ne.3) then
    call initialize_rhs
    call advec1D
else
    do i = 1, nele
        utemp = matmul(Dr,elements(i)%u)
        utemp(:,1) = elements(i)%drdx_cub(:) * utemp(:,1)
        if(ndim.eq.2) utemp(:,1) = utemp(:,1) + elements(i)%drdy_cub(:) * matmul(Dr,elements(i)%u(:,1))
        if(ndim.eq.2) utemp(:,1) = utemp(:,1) + elements(i)%dsdx_cub(:) * matmul(Ds,elements(i)%u(:,1))
        if(ndim.eq.2) utemp(:,1) = utemp(:,1) + elements(i)%dsdy_cub(:) * matmul(Ds,elements(i)%u(:,1))
        elements(i)%u = utemp
    end do
end if

k = 0
do i = 1, nele
    if(.not. inodal) elements(i)%u = matmul(Vander,elements(i)%u)
    do j = 1, ref_ele(elements(i)%typ)%p%nnode_cub
        l = 0
        k = k + 1
        do idir = 1, ndim
            l = l+1
            output_var(k,l) = elements(i)%x_cub(j,idir)
        end do
        do istate = 1, nstate
            l = l+1
            output_var(k,l) = elements(i)%u(j,1)
        end do
    end do
    if(.not. inodal) elements(i)%u = matmul(VanderInv,elements(i)%u)
end do
call output(k, ndim+1, 'final.x')

!call error_evaluate(error_abs, error_rel)
!print*, 'Absolute error: ', error_abs, 'Relative Error: ', error_rel
!print*, 'Relative Error: ', error_rel

! Deallocate matrices
call deallocate_matrices
call finalize_quad
!!!call finalize_grid


!11 FORMAT(' ',4E24.15)

contains

subroutine advec1D
use prec
use glob
use grid
implicit none
integer :: iele, inode, idir
integer :: tstep, nstep, intrk
real(dp) :: timelocal, time, dxmin, CFL, dt
real(dp) :: rk4a(5), rk4b(5), rk4c(5)
real(dp) :: residual


CFL = 0.1d0
time = 0.0d0

dxmin = 1e9
do iele = 1, nele
do inode = 1, elements(iele)%nnode_cub-1
    dt = abs(elements(iele)%x_cub(inode,1) - elements(iele)%x_cub(inode+1,1))
    dxmin = min(dt,dxmin)
end do
end do
dt = 0.5d0 * CFL / maxval(wavespeed) * dxmin

dt= 1.0e-3

nstep = ceiling(finalTime/dt)
dt = finalTime/dble(nstep)
!dt = 0.01d0
!nstep = 1

!dt = 1.000000000000000e-4
!nstep = 10000
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
!print*, 'a', sum(rk4a(1:5))
!print*, 'b', sum(rk4b(1:5))
!print*, 'c', sum(rk4c(1:5))
!do intrk = 1, 5
!   print*,rk4a(intrk), rk4b(intrk), rk4c(intrk)
!end do

do tstep = 1, nstep
    if(tscheme.eq.1) then
        do intrk = 1, 5
            timelocal = time + rk4c(intrk)*dt

            do iele = 1, nele
                istate = 1
                elements(iele)%f(:,istate,1) = wavespeed(1)*elements(iele)%u(:,istate)
                ! Update F
                if(ndim.ge.2) elements(iele)%f(:,istate,2) = wavespeed(2)*elements(iele)%u(:,istate)

                ! Update source term
                elements(iele)%g(:,istate) = 0.0d0
                if(icase.eq.-1) then
                    elements(iele)%g(:,istate) = 0.5d0*2.15d0*cos(2.15d0*elements(iele)%x_cub(:,1)+0.23d0)
                end if
                if(icase.eq.4) then
                    if(ndim.eq.1) elements(iele)%g(:,istate) = -sin(elements(iele)%x_cub(:,1)) + cos(elements(iele)%x_cub(:,1))
                    if(ndim.eq.2) elements(iele)%g(:,istate) = -sin(elements(iele)%x_cub(:,1)) - cos(elements(iele)%x_cub(:,2))
                end if
            end do

            call advecRHS2D
            residual = 0
            do iele = 1, nele
                elements(iele)%resiu = rk4a(intrk)*elements(iele)%resiu 
                elements(iele)%resiu = elements(iele)%resiu + dt*elements(iele)%rhs
                elements(iele)%u = elements(iele)%u + rk4b(intrk)*elements(iele)%resiu
                residual = residual + norm2(elements(iele)%resiu(:,1))
            end do
        end do
    else if(tscheme.eq.2) then
        do iele = 1, nele
            elements(iele)%u_error = elements(iele)%u
            istate = 1
            elements(iele)%f(:,istate,1) = wavespeed(1)*elements(iele)%u(:,istate)
            ! Update F
            if(ndim.ge.2) elements(iele)%f(:,istate,2) = wavespeed(2)*elements(iele)%u(:,istate)

            ! Update source term
            elements(iele)%g(:,istate) = 0.0d0
            if(icase.eq.-1) then
                elements(iele)%g(:,istate) = 0.5d0*2.15d0*cos(2.15d0*elements(iele)%x_cub(:,1)+0.23d0)
            end if
            if(icase.eq.4) then
                if(ndim.eq.1) elements(iele)%g(:,istate) = -sin(elements(iele)%x_cub(:,1)) + cos(elements(iele)%x_cub(:,1))
                if(ndim.eq.2) elements(iele)%g(:,istate) = -sin(elements(iele)%x_cub(:,1)) - cos(elements(iele)%x_cub(:,2))
            end if
        end do
        call advecRHS2D
        do iele = 1, nele
            elements(iele)%resiu = dt*elements(iele)%rhs
            elements(iele)%u = elements(iele)%u + dt*elements(iele)%rhs
        end do
        do iele = 1, nele
            istate = 1
            elements(iele)%f(:,istate,1) = wavespeed(1)*elements(iele)%u(:,istate)
            ! Update F
            if(ndim.ge.2) elements(iele)%f(:,istate,2) = wavespeed(2)*elements(iele)%u(:,istate)

            ! Update source term
            elements(iele)%g(:,istate) = 0.0d0
            if(icase.eq.-1) then
                elements(iele)%g(:,istate) = 0.5d0*2.15d0*cos(2.15d0*elements(iele)%x_cub(:,1)+0.23d0)
            end if
            if(icase.eq.4) then
                if(ndim.eq.1) elements(iele)%g(:,istate) = -sin(elements(iele)%x_cub(:,1)) + cos(elements(iele)%x_cub(:,1))
                if(ndim.eq.2) elements(iele)%g(:,istate) = -sin(elements(iele)%x_cub(:,1)) - cos(elements(iele)%x_cub(:,2))
            end if
        end do
        call advecRHS2D
        residual = 0
        do iele = 1, nele
            elements(iele)%resiu = dt*elements(iele)%rhs
            elements(iele)%u = 0.5*(elements(iele)%u_error + elements(iele)%u + dt*elements(iele)%rhs)
        end do
        !residual = residual + norm2(elements(iele)%resiu(:,1))
    else if(tscheme.eq.0) then
        do iele = 1, nele
            istate = 1
            elements(iele)%f(:,istate,1) = wavespeed(1)*elements(iele)%u(:,istate)
            ! Update F
            if(ndim.ge.2) elements(iele)%f(:,istate,2) = wavespeed(2)*elements(iele)%u(:,istate)

            ! Update source term
            elements(iele)%g(:,istate) = 0.0d0
            if(icase.eq.-1) then
                elements(iele)%g(:,istate) = 0.5d0*2.15d0*cos(2.15d0*elements(iele)%x_cub(:,1)+0.23d0)
            end if
            if(icase.eq.4) then
                if(ndim.eq.1) elements(iele)%g(:,istate) = -sin(elements(iele)%x_cub(:,1)) + cos(elements(iele)%x_cub(:,1))
                if(ndim.eq.2) elements(iele)%g(:,istate) = -sin(elements(iele)%x_cub(:,1)) - cos(elements(iele)%x_cub(:,2))
            end if
        end do
        !call advecRHS2D_weak
        call advecRHS2D
        residual = 0
        do iele = 1, nele
            elements(iele)%resiu = dt*elements(iele)%rhs
            elements(iele)%u = elements(iele)%u + elements(iele)%resiu
            residual = residual + norm2(elements(iele)%resiu(:,1))
        end do

    else if(tscheme.eq.4) then
        do iele = 1, nele
            elements(iele)%u_error = elements(iele)%u
            istate = 1
            elements(iele)%f(:,istate,1) = wavespeed(1)*elements(iele)%u(:,istate)
            ! Update F
            if(ndim.ge.2) elements(iele)%f(:,istate,2) = wavespeed(2)*elements(iele)%u(:,istate)

            ! Update source term
            elements(iele)%g(:,istate) = 0.0d0
            if(icase.eq.-1) then
                elements(iele)%g(:,istate) = 0.5d0*2.15d0*cos(2.15d0*elements(iele)%x_cub(:,1)+0.23d0)
            end if
            if(icase.eq.4) then
                if(ndim.eq.1) elements(iele)%g(:,istate) = -sin(elements(iele)%x_cub(:,1)) + cos(elements(iele)%x_cub(:,1))
                if(ndim.eq.2) elements(iele)%g(:,istate) = -sin(elements(iele)%x_cub(:,1)) - cos(elements(iele)%x_cub(:,2))
            end if
        end do
        call advecRHS2D
        do iele = 1, nele
            elements(iele)%u1 = dt*elements(iele)%rhs
            elements(iele)%u = elements(iele)%u_error+0.5d0*dt*elements(iele)%rhs
        end do
        do iele = 1, nele
            istate = 1
            elements(iele)%f(:,istate,1) = wavespeed(1)*elements(iele)%u(:,istate)
            ! Update F
            if(ndim.ge.2) elements(iele)%f(:,istate,2) = wavespeed(2)*elements(iele)%u(:,istate)

            ! Update source term
            elements(iele)%g(:,istate) = 0.0d0
            if(icase.eq.-1) then
                elements(iele)%g(:,istate) = 0.5d0*2.15d0*cos(2.15d0*elements(iele)%x_cub(:,1)+0.23d0)
            end if
            if(icase.eq.4) then
                if(ndim.eq.1) elements(iele)%g(:,istate) = -sin(elements(iele)%x_cub(:,1)) + cos(elements(iele)%x_cub(:,1))
                if(ndim.eq.2) elements(iele)%g(:,istate) = -sin(elements(iele)%x_cub(:,1)) - cos(elements(iele)%x_cub(:,2))
            end if
        end do
        call advecRHS2D
        do iele = 1, nele
            elements(iele)%u2 = dt*elements(iele)%rhs
            elements(iele)%u = elements(iele)%u_error+0.5d0*dt*elements(iele)%rhs
        end do
        do iele = 1, nele
            istate = 1
            elements(iele)%f(:,istate,1) = wavespeed(1)*elements(iele)%u(:,istate)
            ! Update F
            if(ndim.ge.2) elements(iele)%f(:,istate,2) = wavespeed(2)*elements(iele)%u(:,istate)

            ! Update source term
            elements(iele)%g(:,istate) = 0.0d0
            if(icase.eq.-1) then
                elements(iele)%g(:,istate) = 0.5d0*2.15d0*cos(2.15d0*elements(iele)%x_cub(:,1)+0.23d0)
            end if
            if(icase.eq.4) then
                if(ndim.eq.1) elements(iele)%g(:,istate) = -sin(elements(iele)%x_cub(:,1)) + cos(elements(iele)%x_cub(:,1))
                if(ndim.eq.2) elements(iele)%g(:,istate) = -sin(elements(iele)%x_cub(:,1)) - cos(elements(iele)%x_cub(:,2))
            end if
        end do
        call advecRHS2D
        do iele = 1, nele
            elements(iele)%u3 = dt*elements(iele)%rhs
            elements(iele)%u = elements(iele)%u_error+dt*elements(iele)%rhs
        end do
        do iele = 1, nele
            istate = 1
            elements(iele)%f(:,istate,1) = wavespeed(1)*elements(iele)%u(:,istate)
            ! Update F
            if(ndim.ge.2) elements(iele)%f(:,istate,2) = wavespeed(2)*elements(iele)%u(:,istate)

            ! Update source term
            elements(iele)%g(:,istate) = 0.0d0
            if(icase.eq.-1) then
                elements(iele)%g(:,istate) = 0.5d0*2.15d0*cos(2.15d0*elements(iele)%x_cub(:,1)+0.23d0)
            end if
            if(icase.eq.4) then
                if(ndim.eq.1) elements(iele)%g(:,istate) = -sin(elements(iele)%x_cub(:,1)) + cos(elements(iele)%x_cub(:,1))
                if(ndim.eq.2) elements(iele)%g(:,istate) = -sin(elements(iele)%x_cub(:,1)) - cos(elements(iele)%x_cub(:,2))
            end if
        end do
        call advecRHS2D
        do iele = 1, nele
            elements(iele)%u4 = dt*elements(iele)%rhs
            elements(iele)%u = elements(iele)%u_error+(elements(iele)%u1+2.0d0*elements(iele)%u2+2.0d0*elements(iele)%u3+elements(iele)%u4)/6.0d0
        end do

        !residual = residual + norm2(elements(iele)%resiu(:,1))
    end if
    if(icase.eq.-1 .or. icase.eq.4) then
        if(modulo(tstep,5000).eq.0) print*, tstep, residual
        if(residual.le.5e-15) exit
        if(tstep .ge. 2e7) exit
    end if

    time = time + dt
end do
!print*,dt,nstep,timelocal
if(icase.eq.4) print*, tstep, residual

return

end subroutine advec1D

subroutine output(nrows, ncols, filename)
integer :: i, j, nrows, ncols
character(len=*) :: filename
open(unit=7, file=filename, form='formatted')
do i = 1, nrows
    do j = 1, ncols
        if(abs(output_var(i, j)).lt.1e-30) output_var(i, j) = 0.0d0
    end do
    write(7,11) output_var(i, 1:ncols)
end do
close(7)

return
11 FORMAT(' ',44(E26.17))
end subroutine output

end program adv
