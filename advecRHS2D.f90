module rhs
use prec
use glob
use matrices
use grid
use element_mod

integer :: npts_face_max
integer :: nface_max
!integer, parameter :: npts_face_max = 100
!integer, parameter :: nface_max = 10
!real(dp), dimension(ndim,npts_face_max,nface_global) :: flux_approx1, flux_approx2

integer :: nnode_face
!real(dp) :: ufaces(npts_face_max, nstate, nface_global, nele)
!real(dp) :: ffaces(npts_face_max, nstate, ndim, nface_global, nele)
!real(dp) :: flux_numerical(npts_face_max, nstate, ndim, nface_global)
!!real(dp) :: normal1(npts_face_max, ndim), normal2(npts_face_max, ndim)
!!real(dp), dimension(npts_face_max, nstate) :: uface1, uface2 
!!real(dp), dimension(npts_face_max, nstate, ndim) :: fface1, fface2 
!real(dp), dimension(npts_face_max, nstate, nface_max, nele) :: df_n !strong
!!real(dp), dimension(npts_face_max, nstate, ndim) :: flux_num, flux_app

real(dp), allocatable :: ufaces(:, :, :,:)
real(dp), allocatable :: ffaces(:, :, :,:,:)
real(dp), allocatable :: flux_numerical(:,:,:,:)
real(dp), allocatable :: df_n(:,:,:,:)
!real(dp) :: normal1(npts_face_max, ndim), normal2(npts_face_max, ndim)
!real(dp), dimension(npts_face_max, nstate) :: uface1, uface2 
!real(dp), dimension(npts_face_max, nstate, ndim) :: fface1, fface2 
!real(dp), dimension(npts_face_max, nstate, ndim) :: flux_num, flux_app
contains 

subroutine initialize_rhs
integer :: iele
npts_face_max = -1
nface_max = -1
do iele = 1, nele
    npts_face_max = max(ref_ele(elements(iele)%typ)%p%nnode_face, npts_face_max)
    nface_max = max(ref_ele(elements(iele)%typ)%p%nface, nface_max)
end do
allocate(ufaces(npts_face_max, nstate, nface_global, nele))
allocate(ffaces(npts_face_max, nstate, ndim, nface_global, nele))
allocate(flux_numerical(npts_face_max, nstate, ndim, nface_global))
allocate(df_n(npts_face_max, nstate, nface_max, nele))
end subroutine initialize_rhs

subroutine advecRHS2D
implicit none

integer :: iele, iface, idir, istate
integer :: iele1, iele2, iface1, iface2

! Interpolate state variables at faces
do iele = 1, nele
    nnode_face = ref_ele(elements(iele)%typ)%p%nnode_face
    do iface = 1, ref_ele(elements(iele)%typ)%p%nface
        ufaces(1:nnode_face, 1:nstate, iface, iele) = &
            matmul( transpose(VanderSurf(1:nbasis, 1:nnode_face, iface)), &
                    elements(iele)%u(1:nbasis, 1:nstate))
        do idir = 1, ndim
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
    iface1 = global_faces(iface)%ele_face(1) ! Corresponding face on Element 1
    iface2 = global_faces(iface)%ele_face(2) ! Corresponding face on Element 2

    ! Flux numerical in the same orientation as flux from face 1
    do istate = 1, nstate
    do idir = 1, ndim
        if(global_faces(iface)%reversed) then
            flux_numerical(1:nnode_face, istate, idir, iface) = flux_numerical(1:nnode_face, istate, idir, iface) &
               + 0.5d0*(ffaces(1:nnode_face,istate,idir,iface1,iele1) + ffaces(nnode_face:1:-1,istate,idir,iface2,iele2)) &
               + 0.5d0*wavespeed(idir)*(ufaces(1:nnode_face,istate,iface1,iele1) &
                                        * elements(iele1)%normals(1:nnode_face,iface1,idir) &
                                      + ufaces(nnode_face:1:-1,istate,iface2,iele2) &
                                        * elements(iele2)%normals(nnode_face:1:-1,iface2,idir) )
        else
            flux_numerical(1:nnode_face, istate, idir, iface) = flux_numerical(1:nnode_face, istate, idir, iface) &
               + 0.5d0*(ffaces(1:nnode_face,istate,idir,iface1,iele1) + ffaces(1:nnode_face,istate,idir,iface2,iele2)) &
               + 0.5d0*wavespeed(idir)*(ufaces(1:nnode_face,istate,iface1,iele1) &
                                        * elements(iele1)%normals(1:nnode_face,iface1,idir) &
                                      + ufaces(1:nnode_face,istate,iface2,iele2) &
                                        * elements(iele2)%normals(1:nnode_face,iface2,idir) )
        end if
    end do
    end do
    !print*, 'new num', iface, 1, flux_numerical(1:nnode_face, 1, 1, iface)
    !print*, 'new num', iface, 2, flux_numerical(1:nnode_face, 1, 2, iface)
end do
!
df_n = 0.0d0
do iele = 1, nele
nnode_face = ref_ele(elements(iele)%typ)%p%nnode_face
do iface = 1, ref_ele(elements(iele)%typ)%p%nface
    if(global_faces(elements(iele)%face_pointer(iface))%reversed &
       .and. global_faces(elements(iele)%face_pointer(iface))%ele_conn(2) .eq. iele) then
        do idir = 1, ndim
        do istate = 1, nstate
            df_n(1:nnode_face, istate, iface, iele) = df_n(1:nnode_face, istate, iface, iele) &
                + ( ffaces(1:nnode_face, istate, idir, iface, iele) &
                  - flux_numerical(nnode_face:1:-1, istate, idir, elements(iele)%face_pointer(iface)) ) &
                * elements(iele)%normals(1:nnode_face, iface, idir)
        end do
        end do
    else
        do idir = 1, ndim
        do istate = 1, nstate
            df_n(1:nnode_face, istate, iface, iele) = df_n(1:nnode_face, istate, iface, iele) &
                + ( ffaces(1:nnode_face, istate, idir, iface, iele) &
                  - flux_numerical(1:nnode_face, istate, idir, elements(iele)%face_pointer(iface)) ) &
                * elements(iele)%normals(1:nnode_face, iface, idir)
        end do
        end do
    end if
    !print*, 'new df', iele, iface, df_n(1:nnode_face, 1, iface, iele)
end do
end do

do iele = 1, nele
    elements(iele)%rhs = 0.0d0
    nnode_face = ref_ele(elements(iele)%typ)%p%nnode_face
    do istate = 1, nstate
        if(ndim.eq.1) then
            elements(iele)%rhs(:,istate) = elements(iele)%rhs(:,istate) - elements(iele)%drdx_cub * matmul(Dr,elements(iele)%f(:,istate,1)) 
        else if(ndim.eq.2) then
!           print*
!           print*,elements(iele)%Jac_cub
!           print*,-matmul(matmul(MassInv,Sr),elements(iele)%f(:,istate,1))
            elements(iele)%rhs(:,istate) = elements(iele)%rhs(:,istate) &
                - elements(iele)%drdx_cub * matmul(matmul(MassInv,Sr),elements(iele)%f(:,istate,1)) &
                - elements(iele)%dsdx_cub * matmul(Ds,elements(iele)%f(:,istate,1)) &
                - elements(iele)%drdy_cub * matmul(Dr,elements(iele)%f(:,istate,2)) &
                - elements(iele)%dsdy_cub * matmul(Ds,elements(iele)%f(:,istate,2))
            !print*, elements(iele)%rhs(:,istate)
        end if

        !print*,
        do iface = 1, elements(iele)%nface
            df_n(1:nnode_face, istate, iface, iele) = df_n(1:nnode_face, istate, iface, iele) * elements(iele)%Fscale(1:nnode_face, iface)
            elements(iele)%rhs(:,istate) = elements(iele)%rhs(:,istate) &
                + matmul(matmul(MassInv,Lift(1:nbasis, 1:nnode_face, iface)), df_n(1:nnode_face, istate, iface, iele))
        end do

        !elements(iele)%rhs(:,istate) = elements(iele)%rhs(:,istate) + elements(iele)%g(:,istate) 

    end do
end do

return

end subroutine advecRHS2D

subroutine advecRHS2D_weak
implicit none

integer :: iele, iface, idir, istate
integer :: iele1, iele2, iface1, iface2

! Interpolate state variables at faces
do iele = 1, nele
    nnode_face = ref_ele(elements(iele)%typ)%p%nnode_face
    do iface = 1, ref_ele(elements(iele)%typ)%p%nface
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
    iface1 = global_faces(iface)%ele_face(1) ! Corresponding face on Element 1
    iface2 = global_faces(iface)%ele_face(2) ! Corresponding face on Element 2

    ! Flux numerical in the same orientation as flux from face 1
    do istate = 1, nstate
    do idir = 1, ndim
        if(global_faces(iface)%reversed) then
            flux_numerical(1:nnode_face, istate, idir, iface) = flux_numerical(1:nnode_face, istate, idir, iface) &
               + 0.5d0*(ffaces(1:nnode_face,istate,idir,iface1,iele1) + ffaces(nnode_face:1:-1,istate,idir,iface2,iele2)) &
               + 0.5d0*wavespeed(idir)*(ufaces(1:nnode_face,istate,iface1,iele1) &
                                        * elements(iele1)%normals(1:nnode_face,iface1,idir) &
                                      + ufaces(nnode_face:1:-1,istate,iface2,iele2) &
                                        * elements(iele2)%normals(nnode_face:1:-1,iface2,idir) )
        else
            flux_numerical(1:nnode_face, istate, idir, iface) = flux_numerical(1:nnode_face, istate, idir, iface) &
               + 0.5d0*(ffaces(1:nnode_face,istate,idir,iface1,iele1) + ffaces(1:nnode_face,istate,idir,iface2,iele2)) &
               + 0.5d0*wavespeed(idir)*(ufaces(1:nnode_face,istate,iface1,iele1) &
                                        * elements(iele1)%normals(1:nnode_face,iface1,idir) &
                                      + ufaces(1:nnode_face,istate,iface2,iele2) &
                                        * elements(iele2)%normals(1:nnode_face,iface2,idir) )
        end if
    end do
    end do
    !print*, 'new num', iface, 1, flux_numerical(1:nnode_face, 1, 1, iface)
    !print*, 'new num', iface, 2, flux_numerical(1:nnode_face, 1, 2, iface)
end do
!flux_numerical(1:nnode_face, 1, 1, nface_global) = -0.5d0*9.39645473685325e-01
!
df_n = 0.0d0
do iele = 1, nele
nnode_face = ref_ele(elements(iele)%typ)%p%nnode_face
do iface = 1, ref_ele(elements(iele)%typ)%p%nface
    if(global_faces(elements(iele)%face_pointer(iface))%reversed &
       .and. global_faces(elements(iele)%face_pointer(iface))%ele_conn(2) .eq. iele) then
        do idir = 1, ndim
        do istate = 1, nstate
            df_n(1:nnode_face, istate, iface, iele) = df_n(1:nnode_face, istate, iface, iele) &
                - flux_numerical(nnode_face:1:-1, istate, idir, elements(iele)%face_pointer(iface))  &
                * elements(iele)%normals(1:nnode_face, iface, idir)
        end do
        end do
    else
        do idir = 1, ndim
        do istate = 1, nstate
            df_n(1:nnode_face, istate, iface, iele) = df_n(1:nnode_face, istate, iface, iele) &
                - flux_numerical(1:nnode_face, istate, idir, elements(iele)%face_pointer(iface))  &
                * elements(iele)%normals(1:nnode_face, iface, idir)
        end do
        end do
    end if
    !print*, 'new df', iele, iface, df_n(1:nnode_face, 1, iface, iele)
end do
end do
!df_n(1:nnode_face, 1, 2, nele) = -ffaces(1:nnode_face,1,1,2,nele)

do iele = 1, nele
    elements(iele)%rhs = 0.0d0
    nnode_face = ref_ele(elements(iele)%typ)%p%nnode_face
    do istate = 1, nstate
        if(ndim.eq.1) then
            !elements(iele)%rhs(:,istate) = elements(iele)%rhs(:,istate) + elements(iele)%drdx_cub * matmul(matmul(MassInv,transpose(Sr)),elements(iele)%f(:,istate,1)) 
            !elements(iele)%rhs(:,istate) = elements(iele)%rhs(:,istate) + elements(iele)%drdx_cub * matmul(transpose(Sr),elements(iele)%f(:,istate,1)) 
            elements(iele)%rhs(:,istate) = elements(iele)%rhs(:,istate) + matmul(transpose(Sr),elements(iele)%f(:,istate,1)) 
        else if(ndim.eq.2) then
!           print*,
!           print*,elements(iele)%Jac_cub
!           print*,matmul(matmul(MassInv,transpose(Sr)),elements(iele)%f(:,istate,1))
            elements(iele)%rhs(:,istate) = elements(iele)%rhs(:,istate) &
!               + elements(iele)%drdx_cub*elements(iele)%Jac_cub * matmul(transpose(Sr),elements(iele)%f(:,istate,1))! &
!               + elements(iele)%dsdx_cub*elements(iele)%Jac_cub * matmul(transpose(Ss),elements(iele)%f(:,istate,1)) &
!               + elements(iele)%drdy_cub*elements(iele)%Jac_cub * matmul(transpose(Sr),elements(iele)%f(:,istate,2)) &
!               + elements(iele)%dsdy_cub*elements(iele)%Jac_cub * matmul(transpose(Ss),elements(iele)%f(:,istate,2))
!               + elements(iele)%drdx_cub * matmul(matmul(MassInv,transpose(Sr)),elements(iele)%f(:,istate,1)) !&
!               + 1.0d0/elements(iele)%Jac_cub*matmul(matmul(MassInv,transpose(Sr)),elements(iele)%f(:,istate,1)) !&
                + elements(iele)%drdx_cub * matmul(matmul(MassInv,transpose(Sr)),elements(iele)%f(:,istate,1)) &
                + elements(iele)%dsdx_cub * matmul(matmul(MassInv,transpose(Ss)),elements(iele)%f(:,istate,1)) &
                + elements(iele)%drdy_cub * matmul(matmul(MassInv,transpose(Sr)),elements(iele)%f(:,istate,2)) &
                + elements(iele)%dsdy_cub * matmul(matmul(MassInv,transpose(Ss)),elements(iele)%f(:,istate,2))
        end if
        !print*, iele, elements(iele)%rhs(:,istate)

        !print*, iele, elements(iele)%rhs(:,istate)
        !print*,
        do iface = 1, elements(iele)%nface
!           df_n(1:nnode_face, istate, iface, iele) = df_n(1:nnode_face, istate, iface, iele) * elements(iele)%Fscale(1:nnode_face, iface)
!           elements(iele)%rhs(:,istate) = elements(iele)%rhs(:,istate) &
!               + matmul(Lift(1:nbasis, 1:nnode_face, iface), df_n(1:nnode_face, istate, iface, iele))

            df_n(1:nnode_face, istate, iface, iele) = df_n(1:nnode_face, istate, iface, iele) * elements(iele)%Fscale(1:nnode_face, iface)
            elements(iele)%rhs(:,istate) = elements(iele)%rhs(:,istate) &
                + matmul(matmul(MassInv,Lift(1:nbasis, 1:nnode_face, iface)), df_n(1:nnode_face, istate, iface, iele))
        end do

        !elements(iele)%rhs(:,istate) = elements(iele)%rhs(:,istate) + elements(iele)%g(:,istate) 

!        elements(iele)%rhs(:,istate) = elements(iele)%rhs(:,istate) + matmul(phiTw, elements(iele)%g(:,istate)) * elements(iele)%Jac_cub

!        elements(iele)%rhs(:,istate) = 1.0d0/elements(iele)%Jac_cub * matmul(MassInv,elements(iele)%rhs(:,istate))
        !print*,
        !print*, iele
        !write(*,'(1(E25.17))') elements(iele)%rhs(:,istate)
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
end subroutine advecRHS2D_weak

end module rhs
