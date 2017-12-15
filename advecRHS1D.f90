subroutine advecRHS1D(wavespeed, fscheme)
use prec
use glob
!use poly, only: V, Dr
use matrices
use grid
use element_mod
implicit none

integer, parameter :: npts_face_max = 200
integer, parameter :: nface_max = 10
real(dp), intent(in) :: wavespeed(ndim), fscheme
!real(dp) :: uface(2,nele)
!real(dp) :: du(nface_element*nnode_face,nele)
real(dp) :: flux_numerical(npts_face_max, nstate, ndim, nface_global)
real(dp), dimension(ndim,npts_face_max,nface_global) :: flux_approx1, flux_approx2
real(dp) :: uin, u2, n1, n2

integer :: nnode_face

integer :: iele, iface, idir, istate
integer :: iele1, iele2, iface1, iface2
real(dp) :: ufaces(npts_face_max, nstate, nface_global, nele)
real(dp) :: ffaces(npts_face_max, nstate, ndim, nface_global, nele)
real(dp) :: normal1(npts_face_max, ndim), normal2(npts_face_max, ndim)
real(dp), dimension(npts_face_max, nstate) :: uface1, uface2 
real(dp), dimension(npts_face_max, nstate, ndim) :: fface1, fface2 
real(dp), dimension(npts_face_max, nstate, nface_max, nele) :: df_n !strong
real(dp), dimension(npts_face_max, nstate, ndim) :: flux_num, flux_app

! Interpolate state variables at faces
do iele = 1, nele
    nnode_face = ref_ele(elements(iele)%typ)%p%nnode_face
    do iface = 1, nface_element
        ufaces(1:nnode_face, 1:nstate, iface, iele) = &
            matmul( transpose(VanderSurf(1:nbasis, 1:nnode_face, iface)), &
                    elements(iele)%u(1:nbasis, 1:nstate))
        do idir = 1, ndim
!           call printmatrix(transpose(VanderSurf(1:nbasis, 1:nnode_face,iface)))
!           call printmatrix(elements(iele)%f(1:nbasis, 1:nstate, idir))

            ffaces(1:nnode_face, 1:nstate, idir, iface, iele) = &
                matmul( transpose(VanderSurf(1:nbasis, 1:nnode_face, iface)), &
                        elements(iele)%f(1:nbasis, 1:nstate, idir))
        end do
    end do
end do

flux_numerical = 0.0
do iface = 1, nface_global
    iele1 = global_faces(iface)%ele_conn(1) ! Element 1
    iele2 = global_faces(iface)%ele_conn(2) ! Element 2
    nnode_face = ref_ele(elements(iele1)%typ)%p%nnode_face
    if(nnode_face.ne.ref_ele(elements(iele2)%typ)%p%nnode_face) print*,'Error879'
    iface1 = global_faces(iface)%ele_face(1) ! Corresponding face on Element 1
    iface2 = global_faces(iface)%ele_face(2) ! Corresponding face on Element 2
!   print*, 'gface: ', iface, 'elements: ', iele1, iele2
!   print*, 'faces: ', iface1, iface2

    fface1 = 0.0d0; fface2 = 0.0d0
    normal1(1:nnode_face,1:ndim) = elements(iele1)%normals(1:nnode_face,iface1,1:ndim)
    uface1(1:nnode_face,1:nstate) = ufaces(1:nnode_face,1:nstate,iface1,iele1)
    fface1(1:nnode_face,1:nstate,1:ndim) = ffaces(1:nnode_face,1:nstate,1:ndim,iface1,iele1)
    if(global_faces(iface)%reversed) then
        fface2(nnode_face:1:-1,1:nstate,1:ndim) = ffaces(1:nnode_face,1:nstate,1:ndim,iface2,iele2)
        uface2(nnode_face:1:-1,1:nstate) = ufaces(1:nnode_face,1:nstate,iface2,iele2)
        normal2(nnode_face:1:-1,1:ndim) = elements(iele2)%normals(1:nnode_face,iface2,1:ndim)
    else
        fface2(1:nnode_face,1:nstate,1:ndim) = ffaces(1:nnode_face,1:nstate,1:ndim,iface2,iele2)
        uface2(1:nnode_face,1:nstate) = ufaces(1:nnode_face,1:nstate,iface2,iele2)
        normal2(1:nnode_face,1:ndim) = elements(iele2)%normals(1:nnode_face,iface2,1:ndim)
    end if

    ! Linear advection
    ! U = [u, v], F = [ax*u, 0], F2 = [0, ay*v]
!   F1 = 0.0d0; F2 = 0.0d0
    ! Linear advection
!   istate = 1; idir = 1;
!   F1(istate, idir, 1:nnode_face) = wavespeed(idir) * uface1(istate,1:nnode_face)
!   istate = 2; idir = 1;
!   F1(istate, idir, 1:nnode_face) = 0.0d0

!   istate = 1; idir = 2;
!   F2(istate, idir, 1:nnode_face) = 0.0d0
!   istate = 2; idir = 2;
!   F2(istate, idir, 1:nnode_face) = wavespeed(idir) * uface2(istate,1:nnode_face)

    ! Flux numerical in the same orientation as flux from face 1
    do istate = 1, nstate
    do idir = 1, ndim
        !flux_numerical(1:nnode_face, istate, idir, iface) = flux_numerical(1:nnode_face, istate, idir, iface) &
        !   + 0.5d0*( fface1(1:nnode_face, istate, idir) + fface2(1:nnode_face, istate, idir) ) &
        !   + 0.5d0*( uface1(1:nnode_face, istate) * normal1(1:nnode_face, idir) &
        !           + uface2(1:nnode_face, istate) * normal2(1:nnode_face, idir) )
        if(normal1(1,1).eq.1.0d0) flux_numerical(1:nnode_face, istate, idir, iface) = fface1(1:nnode_face, istate, idir)
        if(normal2(1,1).eq.1.0d0) flux_numerical(1:nnode_face, istate, idir, iface) = fface2(1:nnode_face, istate, idir)
    end do
    end do
    !flux_numerical(1:ndim,1:nnode_face,iface) = fluxnum_eval(F1, F2, normal1, normal2, wavespeed, fscheme)
    !print*, iface, nface_global, flux_numerical(1:nnode_face, 1, 1, iface)
end do
!
df_n = 0.0d0
do iele = 1, nele
nnode_face = ref_ele(elements(iele)%typ)%p%nnode_face
do iface = 1, ref_ele(elements(iele)%typ)%p%nface

    normal1(1:nnode_face,1:ndim) = elements(iele)%normals(1:nnode_face, iface, 1:ndim)
    flux_app(1:nnode_face, 1:nstate, 1:ndim) = ffaces(1:nnode_face, 1:nstate, 1:ndim, iface, iele)
    if(global_faces(elements(iele)%face_pointer(iface))%reversed) then
        flux_num(nnode_face:1:-1, 1:nstate, 1:ndim) = flux_numerical(1:nnode_face, 1:nstate, 1:ndim, elements(iele)%face_pointer(iface))
    else
        flux_num(1:nnode_face, 1:nstate, 1:ndim) = flux_numerical(1:nnode_face, 1:nstate, 1:ndim, elements(iele)%face_pointer(iface))
    end if

    do idir = 1, ndim
    do istate = 1, nstate
        df_n(1:nnode_face, istate, iface, iele) = df_n(1:nnode_face, istate, iface, iele) &
            + (flux_app(1:nnode_face, istate, idir) - flux_num(1:nnode_face, istate, idir)) * normal1(1:nnode_face, idir)
    end do
    end do
    !if(iface.eq.1) df_n(1:nnode_face, 1,1, iele) = -(flux_app(1:nnode_face, 1,1) - flux_num(1:nnode_face,1,1))
    !if(iface.eq.2) df_n(1:nnode_face, 1,2, iele) = 0.0d0
    !if(iface.eq.2) print*,df_n(1:nnode_face, 1, iface, iele)
end do
end do

do iele = 1, nele
    elements(iele)%rhs = 0.0d0
    nnode_face = ref_ele(elements(iele)%typ)%p%nnode_face
    do istate = 1, nstate
        if(ndim.eq.1) then
            elements(iele)%rhs(:,istate) = elements(iele)%rhs(:,istate) - elements(iele)%drdx_cub * matmul(Dr,elements(iele)%f(:,istate,1)) 
        else if(ndim.eq.2) then
            elements(iele)%rhs(:,istate) = elements(iele)%rhs(:,istate) &
                - elements(iele)%drdx_cub * matmul(Dr,elements(iele)%f(:,istate,1)) &
                - elements(iele)%dsdx_cub * matmul(Ds,elements(iele)%f(:,istate,1)) &
                - elements(iele)%drdy_cub * matmul(Dr,elements(iele)%f(:,istate,2)) &
                - elements(iele)%dsdy_cub * matmul(Ds,elements(iele)%f(:,istate,2))

        end if
        do iface = 1, elements(iele)%nface
            df_n(1:nnode_face, istate, iface, iele) = df_n(1:nnode_face, istate, iface, iele) * elements(iele)%Fscale(1:nnode_face, iface)
            elements(iele)%rhs(:,istate) = elements(iele)%rhs(:,istate) &
                + matmul(Lift(1:nbasis, 1:nnode_face, iface), df_n(1:nnode_face, istate, iface, iele))
        end do
    end do
end do

! Evaluate RHS of semi-discrete PDE
! dudt = -wavespeed * Minv * S * u + Minv[ l(x) * du ]^xr_xl
!rhsu = -wavespeed * drdx_cub * matmul(Differentiation,elements%u)
!rhsu = rhsu + matmul(Lift, Fscale*du)
!
!
!
 return
!
!contains
!
!subroutine flux_numerical_eval(uface, wavespeed, fscheme, flux_numerical)
!use prec
!use glob
!use grid, only: normals
!implicit none
!
!real(dp), intent(in) :: uface(2,nele)
!real(dp), intent(in) :: wavespeed, fscheme
!real(dp), intent(out) :: flux_numerical(nele+1)
!integer :: iface, iele
!
!!! Upwind flux
!flux_numerical = 0.0
!do iface = 2, nele
!    iele = iface-1
!    flux_numerical(iface) = uface(2,iele) + uface(1,iele+1) 
!    flux_numerical(iface) = flux_numerical(iface) + (1.0-fscheme) &
!        * (uface(2,iele) * normals(2,iele) + uface(1,iele+1)*normals(1,iele+1))
!end do
!flux_numerical = flux_numerical * wavespeed * 0.5d0
!
!! Lax-Friedrichs
!do iface = 2, nele
!    iele = iface-1
!    flux_numerical(iface) = 0.5d0*wavespeed*(uface(2,iele) + uface(1,iele+1))
!    flux_numerical(iface) = flux_numerical(iface) + wavespeed * 0.5d0 * &
!        (uface(2,iele) * normals(2,iele) + uface(1,iele+1)*normals(1,iele+1))
!end do
!
!do iface = 2, nele
!    iele = iface-1
!    flux_numerical(iface) = 0.5d0*wavespeed*(uface(2,iele) + uface(1,iele+1))
!    flux_numerical(iface) = flux_numerical(iface) + wavespeed * 0.5d0 * &
!        (uface(2,iele) * normals(2,iele) + uface(1,iele+1)*normals(1,iele+1))
!end do
!
!return
!end subroutine flux_numerical_eval
!real(dp) function flux_eval(u1,u2,n1,n2,A,fscheme)
!    implicit none
!    real(dp), intent(in) :: u1, u2, n1, n2, A, fscheme
!    flux_eval = 0.5d0*A*(u2 + u1) + (1.0d0-fscheme) * A * 0.5d0 * (u1*n1 + u2*n2)
!end function flux_eval
!
!real(dp) function eval_flux_analytical(u, wavespeed)
!use prec
!use glob
!implicit none
!
!real(dp), intent(in) :: u, wavespeed
!
!select case(iequation)
!    case(1)
!        flux_analytical = wavespeed*u
!    case default
!end select
!
!end function eval_flux_analytical
!
 end subroutine advecRHS1D
!
