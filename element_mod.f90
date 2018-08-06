module element_mod
use glob
use ref_ele_mod
!use quadrature, only: quad_npts
implicit none

type element
    integer :: typ
    integer :: nface
    integer :: nnode_cub, nnode_face
    integer, allocatable :: vertices(:)
    integer, allocatable :: face_pointer(:)
    real(dp), allocatable :: u(:,:), f(:,:,:), g(:,:)
    real(dp), allocatable :: u1(:,:)
    real(dp), allocatable :: u2(:,:)
    real(dp), allocatable :: u3(:,:)
    real(dp), allocatable :: u4(:,:)
    real(dp), allocatable :: u_exact(:,:), u_error(:,:)
    real(dp), allocatable :: x_cub(:,:), x_face(:,:,:)
    real(dp), allocatable :: resiu(:,:), rhs(:,:)
    real(dp), allocatable :: normals(:,:,:)
    real(dp), allocatable, dimension(:) :: drdx_cub, drdy_cub, dsdx_cub, dsdy_cub, dxdr_cub, dydr_cub, dxds_cub, dyds_cub
    real(dp), allocatable, dimension(:,:) :: drdx_face, drdy_face, dsdx_face, dsdy_face, dxdr_face, dydr_face, dxds_face, dyds_face
    real(dp), allocatable :: Jac_cub(:), Jac_face(:,:), Fscale(:,:)
end type element

contains

    type(element) function init_ele(typ) result(ele)
    implicit none

    integer, intent(in) :: typ
    integer :: nbasis
    type(ref) :: refe

    ele%typ = typ
    refe = ref_ele(ele%typ)%p
    nbasis = refe%nbasis
    ele%nface = refe%nface
!   ele%nvertex = ref(ele%typ)%nvertex
    ele%nnode_cub = refe%nnode_cub
    ele%nnode_face = refe%nnode_face

    allocate(ele%face_pointer(refe%nface))
    allocate(ele%vertices(refe%nvertex))
    ele%face_pointer(:) = -1
    ele%vertices(:) = -1

    !allocate(ele%normals(refe%nface, ndim))
    allocate(ele%normals(ele%nnode_face, refe%nface, ndim))
    ! Hardcoded normals for now
    if(ndim.eq.1) ele%normals(:, 1, 1) = -1.0d0
    if(ndim.eq.1) ele%normals(:, 2, 1) = 1.0d0

    ! Hardcoded 2D grid. Only accept square grid
    if(ndim.eq.2) ele%normals(:, :, :) = 0.0d0
    if(ndim.eq.2) ele%normals(:, 1, 2) = -1.0d0
    if(ndim.eq.2) ele%normals(:, 2, 1) = 1.0d0
    if(ndim.eq.2) ele%normals(:, 3, 2) = 1.0d0
    if(ndim.eq.2) ele%normals(:, 4, 1) = -1.0d0

    allocate(ele%u(nbasis, nstate))
    allocate(ele%u1(nbasis, nstate))
    allocate(ele%u2(nbasis, nstate))
    allocate(ele%u3(nbasis, nstate))
    allocate(ele%u4(nbasis, nstate))
    allocate(ele%f(nbasis, nstate, ndim))
    allocate(ele%g(nbasis, nstate))

    allocate(ele%u_exact(nbasis, nstate))
    allocate(ele%u_error(nbasis, nstate))

    allocate(ele%x_cub(ele%nnode_cub, ndim))
    allocate(ele%x_face(ele%nnode_face, ndim, refe%nface))

    allocate(ele%resiu(nbasis, nstate))
    allocate(ele%rhs(nbasis, nstate))

    allocate(ele%drdx_cub(ele%nnode_cub))
    allocate(ele%drdy_cub(ele%nnode_cub))
    allocate(ele%dsdx_cub(ele%nnode_cub))
    allocate(ele%dsdy_cub(ele%nnode_cub))

    allocate(ele%dxdr_cub(ele%nnode_cub))
    allocate(ele%dxds_cub(ele%nnode_cub))
    allocate(ele%dydr_cub(ele%nnode_cub))
    allocate(ele%dyds_cub(ele%nnode_cub))

    allocate(ele%Jac_cub(ele%nnode_cub))


    allocate(ele%drdx_face(ele%nnode_face, ele%nface))
    allocate(ele%drdy_face(ele%nnode_face, ele%nface))
    allocate(ele%dsdx_face(ele%nnode_face, ele%nface))
    allocate(ele%dsdy_face(ele%nnode_face, ele%nface))

    allocate(ele%dxdr_face(ele%nnode_face, ele%nface))
    allocate(ele%dxds_face(ele%nnode_face, ele%nface))
    allocate(ele%dydr_face(ele%nnode_face, ele%nface))
    allocate(ele%dyds_face(ele%nnode_face, ele%nface))

    allocate(ele%Jac_face(ele%nnode_face, ele%nface))

    allocate(ele%Fscale(ele%nnode_face, ele%nface))


    return
    end function init_ele

!   subroutine element_jacobian_eval(ele)
!   use matrices, only: Dr, Ds, Vander, VanderInv
!   implicit none
!   type(element) :: ele
!   type(ref) :: refe


!   refe = ref_ele(ele%typ)%p

!   if(inodal) then
!       ele%dxdr_cub = matmul(Dr, ele%x_cub(:,1))
!       if(ndim.ge.2) then
!           ele%dydr_cub = matmul(Dr, ele%x_cub(:,2))
!           ele%dxds_cub = matmul(Ds, ele%x_cub(:,1))
!           ele%dyds_cub = matmul(Ds, ele%x_cub(:,2))
!       end if
!   else
!       ele%dxdr_cub = matmul(Vander,matmul(Dr,matmul(VanderInv, ele%x_cub(:,1))))
!       if(ndim.ge.2) then
!           ele%dxds_cub = matmul(Vander,matmul(Ds,matmul(VanderInv, ele%x_cub(:,1))))
!           ele%dydr_cub = matmul(Vander,matmul(Dr,matmul(VanderInv, ele%x_cub(:,2))))
!           ele%dyds_cub = matmul(Vander,matmul(Ds,matmul(VanderInv, ele%x_cub(:,2))))
!       end if
!   end if
!   if(ndim.eq.1) then
!       ele%drdx_cub = 1.0d0/ele%dxdr_cub
!   else if(ndim.ge.2) then
!       ele%Jac_cub = ele%dxdr_cub*ele%dyds_cub - ele%dxds_cub*ele%dydr_cub
!       ele%drdx_cub = ele%dyds_cub/ele%Jac_cub
!       ele%dsdx_cub = -ele%dydr_cub/ele%Jac_cub
!       ele%drdy_cub = -ele%dxds_cub/ele%Jac_cub
!       ele%dsdy_cub = ele%dxdr_cub/ele%Jac_cub
!   end if

!   if(inodal) then
!       ele%dxdr_face = matmul(Dr, ele%x_face(:,1))
!       if(ndim.ge.2) then
!           ele%dydr_face = matmul(Dr, ele%x_face(:,2))
!           ele%dxds_face = matmul(Ds, ele%x_face(:,1))
!           ele%dyds_face = matmul(Ds, ele%x_face(:,2))
!       end if
!   else
!       ele%dxdr_face = matmul(Vander,matmul(Dr,matmul(VanderInv, ele%x_face(:,1))))
!       if(ndim.ge.2) then
!           ele%dxds_face = matmul(Vander,matmul(Ds,matmul(VanderInv, ele%x_face(:,1))))
!           ele%dydr_face = matmul(Vander,matmul(Dr,matmul(VanderInv, ele%x_face(:,2))))
!           ele%dyds_face = matmul(Vander,matmul(Ds,matmul(VanderInv, ele%x_face(:,2))))
!       end if
!   end if
!   if(ndim.eq.1) then
!       ele%drdx_face = 1.0d0/ele%dxdr_face
!   else if(ndim.ge.2) then
!       ele%Jac_face = ele%dxdr_face*ele%dyds_face - ele%dxds_face*ele%dydr_face
!       ele%drdx_face = ele%dyds_face/ele%Jac_face
!       ele%dsdx_face = -ele%dydr_face/ele%Jac_face
!       ele%drdy_face = -ele%dxds_face/ele%Jac_face
!       ele%dsdy_face = ele%dxdr_face/ele%Jac_face
!   end if

!   end subroutine

end module element_mod
