program adv

use prec
use glob
use grid
!use poly
use matrices
use quadrature
use cubature2D
use ref_ele_mod
!use stability

implicit none

real(dp) :: xmin, xmax, ymin, ymax
!real(dp), allocatable :: xref(:,:), weights(:)
integer :: tscheme
real(dp) :: fscheme, wavespeed(ndim)
real(dp) :: finalTime
integer :: i, j, k, l
logical :: checkDiff

real(dp) :: output_var(10000,5)
!real(dp), allocatable :: u(:,:,:)
real(dp), allocatable :: utemp(:,:), du(:,:)
real(dp) :: error
!real(dp), allocatable :: xsol(:,:), usol(:,:)
!real(dp), allocatable :: VanderQuad(:,:)
!real(dp), allocatable :: polyVal(:)

!Stability analysis
!real(dp) :: wavenumber
!complex(dp), allocatable :: eig(:)

order   = 6
nele    = 50
open(unit=1,file="input.in",form='formatted')
read(1,*) order, nivertex, njvertex
close(1)
tscheme = 1
fscheme = 0.0d0
!fscheme = 1.0d0
wavespeed = 1.0d0
finalTime = 0.50d0

if(ndim.eq.1) nbasis = order + 1
if(ndim.eq.2) nbasis = (order + 1)**2
!nvertex = nele + 1
!nface   = 2
!nfpoint = 1

! Setup Quadrature
call initialize_quad
do i = 1, quad_npts
    output_var(i, 1) = quadR(i)
    output_var(i, 2) = quadW(i)
end do
call output(quad_npts, 2, 'quad.x')
call initialize_cub2D

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

!do i = 1, nbasis
!!print*, elements(1)%x_cub(1,i), elements(1)%x_cub(1,i)
!if(ndim.eq.1) print*, ref_ele(ndim)%p%nodes_cub(1,i)
!if(ndim.eq.1) print*, elements(1)%x_cub(i,1)
!if(ndim.eq.2) print*, ref_ele(ndim)%p%nodes_cub(1,i), ref_ele(elements(1)%typ)%p%nodes_cub(2,i)
!end do

!print*,'dxdr_cub'
!print*, elements(1)%dxdr_cub
!if(ndim.ge.2) then
!    print*, elements(1)%dxds_cub
!    print*, elements(1)%dydr_cub
!    print*, elements(1)%dyds_cub
!end if
!print*,
!print*, elements(1)%drdx_cub
!if(ndim.ge.2) then
!    print*, elements(1)%dsdx_cub
!    print*, elements(1)%drdy_cub
!    print*, elements(1)%dsdy_cub
!end if
!print*,
!print*,

!print*, 'Vandermonde'
!call printmatrix(Vander)
!
!do i=1,ndim; call printmatrix(VanderGrad(:,:,i)); print*,' ' ; enddo;
!print*,
!print*, 'M = (V*Vt)^-1'
!if(inodal) call printmatrix(matmul(transpose(VanderInv),VanderInv))
!if(.not. inodal) call printmatrix(matmul(transpose(VanderInv),VanderInv))
!print*,
!print*, 'Mass'
!call printmatrix(Mass)
!print*,
!print*, 'Sr'
!call printmatrix(Sr)
!print*, 'Dr'
!call printmatrix(Dr)
!print*,
!!call printmatrix(matmul(VanderInv,VanderGrad(:,:,1)))
!call printmatrix(matmul(VanderGrad(:,:,1),VanderInv))
!print*,
!if(ndim.eq.2) then
!    call printmatrix(Ss)
!    print*,
!    call printmatrix(Ds)
!    print*,
!end if
!print*, 'Lift'
!do i = 1, ref_ele(ndim)%p%nface
!    call printmatrix(Lift(:,:,i))
!    print*,
!end do

!Dr = matmul(VanderGrad(:,:,1),VanderInv)
!Ds = matmul(VanderGrad(:,:,2),VanderInv)
do i = 1, nele
    call initialize_u(elements(i)%u, elements(i)%x_cub, nbasis)
    ! umodal = VanderInv * u
    if(.not. inodal) elements(i)%u = matmul(VanderInv,elements(i)%u)
end do

checkDiff = 1==0
if(checkDiff) call printDiff

k = 0
do i = 1, nele
    if(.not. inodal) elements(i)%u = matmul(Vander,elements(i)%u)
    do j = 1, ref_ele(elements(i)%typ)%p%nnode_cub
        k = k + 1
        if(ndim.eq.1) then
            output_var(k,1) = elements(i)%x_cub(j,1)
            output_var(k,2) = elements(i)%u(j,1)
        end if
    end do
    if(.not. inodal) elements(i)%u = matmul(VanderInv,elements(i)%u)
end do
call output(k, 2, 'init.x')

call advec1D(tscheme, fscheme, wavespeed, finalTime)

k = 0
do i = 1, nele
    if(.not. inodal) elements(i)%u = matmul(Vander,elements(i)%u)
    do j = 1, ref_ele(elements(i)%typ)%p%nnode_cub
        k = k + 1
        if(ndim.eq.1) then
            output_var(k,1) = elements(i)%x_cub(j,1)
            output_var(k,2) = elements(i)%u(j,1)
        end if
    end do
    if(.not. inodal) elements(i)%u = matmul(VanderInv,elements(i)%u)
end do
call output(k, 2, 'final.x')

!!!!call printmatrix(matmul(Mass,MassInv))
!!!!print*,
!!!!call printmatrix(matmul(Vander,VanderInv))
!!!!print*,
!!!
!!!! Old subroutines replaced by more general formulations
!!!!call modalToNodal(xref, quad_alpha, quad_beta, order, Vandermonde, VandermondeInv)
!!!!call massToStiff(xref, quad_alpha, quad_beta, order, VanderInv, Differentiation)
!!!!call lift1d(Vander, Lift)
!!!
!!!call genGrid(xref, xmin, xmax)
!!!
!!!! Stability analysis
!!!select case(istab)
!!!    case(1)
!!!        allocate(eig(nbasis))
!!!        wavenumber = 1.0
!!!        call check_stability(wavenumber, eig)
!!!        stop
!!!    case(2)
!!!        allocate(eig(nbasis))
!!!        j = 100
!!!        do i = 0, j
!!!            wavenumber = 2.0d0*PI/j * i - PI
!!!            call check_stability(wavenumber, eig)
!!!        end do
!!!        stop
!!!    case(3)
!!!        call stable_dt
!!!    case default
!!!end select
!!!

!!!allocate(xsol(quad_npts,nele))
!!!allocate(usol(quad_npts,nele))
!!!allocate(VanderQuad(quad_npts, nbasis))
!!!allocate(polyVal(quad_npts))
!!!! Evaluate Vandermonde matrix at quadrature (output) points
!!!do j = 0, order
!!!    call poly_eval(polytype, j, quadR, polyVal)
!!!    VanderQuad(1:quad_npts,j+1) = polyVal
!!!end do
!!!if(inodal) u = matmul(VanderInv,u)
!!!usol = matmul(VanderQuad,u)
!!!if(inodal) u = matmul(Vander,u)
!!!do i = 1, quad_npts
!!!do j = 1, nele
!!!    xsol(i,j) = vertices(EtoV(1,j)) &
!!!        + 1.0d0/(refb-refa) &!quadR(quad_npts)-quadR(1)) &
!!!        * (quadR(i)-refa) * (vertices(EtoV(2,j)) - vertices(EtoV(1,j)))
!!!end do
!!!end do
!!!
!!!open(unit=7, file='output.dat', form='formatted')
!!!do j = 1, nele
!!!do i = 1, quad_npts
!!! if(abs(usol(i,j)).le.1e-30) usol(i,j) = 0.0d0
!!! write(7,11) xsol(i,j), usol(i,j)
!!!end do
!!!end do
!!!
!!!!if(inodal) call printmatrix(matmul(Differentiation,u))
!!!!if(.not.inodal) call printmatrix(matmul(Vander,matmul(Differentiation,u)))
!!!!stop
!!!if(icase.lt.3) then
!!!    call advec1D(u, tscheme, fscheme, wavespeed, finalTime)
!!!    ! umodal = VanderInv * u
!!!    if(inodal) u = matmul(VanderInv,u)
!!!    usol = matmul(VanderQuad,u)
!!!    if(inodal) u = matmul(Vander,u)
!!!else if(icase.eq.3) then
!!!    u = drdx_cub*matmul(Differentiation,u)
!!!    if(inodal) u = matmul(VanderInv,u)
!!!    usol = matmul(VanderQuad,u)
!!!    if(inodal) u = matmul(Vander,u)
!!!else if(icase.eq.4) then
!!!    call bspline_init
!!!    do k = 1, bs_nbasis
!!!        do l = 1, nele
!!!        call bsplineP(x_cub(:,:,l), k, bs_order, u(:,:,l))
!!!        end do
!!!        if(inodal) u = matmul(VanderInv,u)
!!!        usol = matmul(VanderQuad,u)
!!!        if(inodal) u = matmul(Vander,u)
!!!        do j = 1, nele
!!!        do i = 1, quad_npts
!!!         if(abs(usol(i,j)).le.1e-30) usol(i,j) = 0.0d0
!!!         write(7,11) xsol(i,j), usol(i,j)
!!!        end do
!!!        end do
!!!    end do
!!!else if(icase.eq.5) then
!!!    call nurbs_init
!!!    do k = 1, nurbs_nbasis
!!!        do l = 1, nele
!!!        call nurbsP(x_cub(:,:,l), k, nurbs_order, u(:,:,l))
!!!        end do
!!!        if(inodal) u = matmul(VanderInv,u)
!!!        usol = matmul(VanderQuad,u)
!!!        if(inodal) u = matmul(Vander,u)
!!!        do j = 1, nele
!!!        do i = 1, quad_npts
!!!         if(abs(usol(i,j)).le.1e-30) usol(i,j) = 0.0d0
!!!         write(7,11) xsol(i,j), usol(i,j)
!!!        end do
!!!        end do
!!!    end do
!!!end if
!!!
!!!do j = 1, nele
!!!do i = 1, quad_npts
!!! if(abs(usol(i,j)).le.1e-30) usol(i,j) = 0.0d0
!!! write(7,11) xsol(i,j), usol(i,j)
!!!end do
!!!end do
!!!do i = 1, quad_npts
!!! write(7,11) quadR(i), quadW(i)
!!!end do
!!!
!!!!do j = 1, nele
!!!!do i = 1, nbasis
!!!! if(abs(u(i,j)).le.1e-30) u(i,j) = 0.0d0
!!!! write(7,11) x(i,j), u(i,j)
!!!!end do
!!!!end do
!!!!! weights
!!!!do i = 1, nbasis
!!!! write(7,11) quadR(i), quadW(i)
!!!!end do
!!!
!!!
!!!
!!!close(7)
!!!
!!!
!!!! Deallocate matrices
call deallocate_matrices
call finalize_quad
!!!call finalize_grid


!11 FORMAT(' ',4E24.15)

contains

subroutine advec1D(tscheme, fscheme, wavespeed, finalTime)
use prec
use glob
use grid
implicit none
integer, intent(in) :: tscheme
real(dp), intent(in) :: finalTime, wavespeed(ndim), fscheme

integer :: iele, inode, idir
integer :: tstep, nstep, intrk
real(dp) :: timelocal, time, dxmin, CFL, dt
real(dp) :: rk4a(5), rk4b(5), rk4c(5)


CFL = 0.1d0
time = 0.0d0

dxmin = 1e9
do iele = 1, nele
do inode = 1, elements(iele)%nnode_cub-1
    dt = abs(elements(iele)%x_cub(inode,1) - elements(iele)%x_cub(inode+1,1))
    dxmin = min(dt,dxmin)
end do
end do
dt = 0.5d0 * CFL / wavespeed(1) * dxmin
!dt = 1.0000000000000001e-5

nstep = ceiling(finalTime/dt)
dt = finalTime/nstep
!nstep = 0
!do i = 1, nele
!    utemp = matmul(Dr,elements(i)%u)
!    if(.not.inodal) utemp = matmul(Vander, utemp)
!    utemp(:,1) = elements(i)%drdx_cub(:) * utemp(:,1)
!    elements(i)%u = utemp
!end do

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
    do intrk = 1, 5
        timelocal = time + rk4c(intrk)*dt
        ! Update F
        do iele = 1, nele
            do idir = 1, ndim
                elements(iele)%f(:,:,idir) = elements(iele)%u
            end do
        end do
        call advecRHS1D(elements, wavespeed, fscheme)
        do iele = 1, nele
            elements(iele)%resiu = rk4a(intrk)*elements(iele)%resiu + dt*elements(iele)%rhs
            elements(iele)%u = elements(iele)%u + rk4b(intrk)*elements(iele)%resiu
        end do
    end do
    time = time + dt
end do

return

end subroutine advec1D

subroutine initialize_u(u, x, npts)
    use glob
    use prec
    implicit none
    integer, intent(in) :: npts
    real(dp), intent(in) :: x(npts,ndim)
    real(dp), intent(out) :: u(npts,nstate)
    integer :: i, j, istate
    u = 0.0d0
    select case(icase)
        case(0)
            istate = 1
            if(ndim.eq.1) u(:,1) = sin(x(:,1))
            if(ndim.eq.2) u(:,istate) = sin(x(:,1)) + cos(x(:,2)) + 1
        case(1)
            istate = 1
            if(ndim.eq.1) u(:,istate) = sin(x(:,1)) + 1
            if(ndim.eq.2) then
                istate = 1
                u(:,istate) = sin(x(:,1))
            end if
!           if(ndim.eq.2) u(:,1) = sin(x(:,1))
    !   case(1)
    !       do j = 1, nele
    !       do i = 1, nbasis
    !        if(abs(x(i,j)).lt.1.0) then
    !           u(:,i,j) = exp(-1.0d0/(1.0d0-x(:,i,j)**2))
    !        end if
    !       end do
    !       end do
    !   case(2)
    !       u = x/abs(x)
        case default
            print *, 'Invalid case'
    end select

end subroutine initialize_u

subroutine output(nrows, ncols, filename)
integer :: i, j, nrows, ncols
character(len=*) :: filename
open(unit=7, file=filename, form='formatted')
do i = 1, nrows
    write(7,11) output_var(i, 1:ncols)
end do
close(7)

return
11 FORMAT(' ',44E24.15)
end subroutine output
!subroutine output(filename)
!integer :: iele, inode
!type(element) :: ele
!character(len=*) :: filename
!open(unit=7, file=filename, form='formatted')
!do iele = 1, nele
!    ele = elements(iele)
!    do inode = 1, ref_ele(ele%typ)%p%nnode_cub
!        if(ndim.eq.1) write(7,11) ele%x_cub(inode,1), ele%u(inode,1)
!        if(ndim.eq.2) write(7,11) ele%x_cub(inode,1), ele%x_cub(inode,2)
!    end do
!end do
!
!return
!11 FORMAT(' ',4E24.15)
!end subroutine output


subroutine printDiff
print*,'Differentiation solution'
allocate(utemp(nbasis,nstate))
allocate(du(nbasis,nstate))
do i = 1, nele
    if(ndim.eq.1) then
        utemp = matmul(Dr,elements(i)%u)
        if(.not.inodal) utemp = matmul(Vander, utemp)
        utemp(:,1) = elements(i)%drdx_cub(:) * utemp(:,1)
        do j = 1, nbasis
            print*, elements(i)%x_cub(j,1), &
                    'App:', utemp(j,:), &
                    'Exact:', cos(elements(i)%x_cub(j,1)) &
                   ,'Err:', abs(utemp(j,1)-cos(elements(i)%x_cub(j,1)))
        end do
    else if(ndim.eq.2) then
        error = 0.0d0
        ! du/dx = du/dr dr/dx + du/ds ds/dx
        du = matmul(Dr,elements(i)%u)
        if(.not.inodal) du = matmul(Vander, du)
        utemp(:,1) = elements(i)%drdx_cub(:) * du(:,1)

        du = matmul(Ds,elements(i)%u)
        if(.not.inodal) du = matmul(Vander, du)
        utemp(:,1) = utemp(:,1) + elements(i)%dsdx_cub(:) * du(:,1)

        do j = 1, nbasis
!           print*, elements(i)%x_cub(j,1), elements(i)%x_cub(j,2) &
!                   ,'App:', utemp(j,1) &
!                   ,'Ext:', cos(elements(i)%x_cub(j,1)) &
!                   ,'Err:', abs(utemp(j,1)-cos(elements(i)%x_cub(j,1)))
!           print*, elements(i)%x_cub(j,1), elements(i)%x_cub(j,2) &
!                   ,'fun:', sin(elements(i)%x_cub(j,1))+cos(elements(i)%x_cub(j,2)) + 1 &
!                   ,'app:', elements(i)%u(j,:) &
!                   ,'Err:', sin(elements(i)%x_cub(j,1))+cos(elements(i)%x_cub(j,2)) + 1 - elements(i)%u(j,:)
           error = error + abs(utemp(j,1)-cos(elements(i)%x_cub(j,1)))
        end do
        print*, 'Element', i, 'du/dx error: ',error

        error = 0.0d0
        ! du/dy = du/dr dr/dy + du/ds ds/dy
        du = matmul(Dr,elements(i)%u)
        if(.not.inodal) du = matmul(Vander, du)
        utemp(:,1) = elements(i)%drdy_cub(:) * du(:,1)

        du = matmul(Ds,elements(i)%u)
        if(.not.inodal) du = matmul(Vander, du)
        utemp(:,1) = utemp(:,1) + elements(i)%dsdy_cub(:) * du(:,1)

        do j = 1, nbasis
!           print*, elements(i)%x_cub(j,1), elements(i)%x_cub(j,2) &
!                   ,'App:', utemp(j,1) &
!                   ,'Ext:', cos(elements(i)%x_cub(j,1)) &
!                   ,'Err:', abs(utemp(j,1)-cos(elements(i)%x_cub(j,1)))
           error = error + abs(utemp(j,1)-(-sin(elements(i)%x_cub(j,2))))
        end do
        print*, 'Element', i, 'du/dy error: ',error
    end if !ndim
end do
end subroutine printDiff
end program adv
