module ref_ele_mod
use glob
use quadrature, only: quadR, quad_npts
use cubature2D, only: cubRS, cub2D_npts
implicit none

type ref
    integer :: nbasis, order
    integer :: nvertex, nface, nsurfnode
    integer :: nnode_cub, nnode_face
    real(dp), allocatable :: nodes_cub(:,:)
    real(dp), allocatable :: nodes_face(:,:,:)
    real(dp), allocatable :: vertices(:,:)
end type ref


type(ref), target :: line, quad

type refptr
    type(ref), pointer :: p
end type refptr
type(refptr), dimension(3) :: ref_ele

! pointer to array, not array of pointers
! type(ref), pointer :: ref_ele(3)

contains

subroutine initialize_reference_elements
    implicit none
    integer :: i 

    ref_ele(1)%p => line
    ref_ele(2)%p => quad
    ! 1D Line
    line%nface = 2
    line%order = order
    line%nbasis = order+1

    line%nvertex = 2
    allocate(line%vertices(line%nvertex, ndim))
    line%vertices(1, 1) = -1.0d0
    line%vertices(2, 1) = 1.0d0

    line%nnode_cub = quad_npts
    allocate(line%nodes_cub(line%nnode_cub, 1))
    line%nodes_cub(:, 1) = quadR

    line%nnode_face = 1
    line%nsurfnode = line%nnode_face * line%nface
    allocate(line%nodes_face(line%nnode_face, ndim, line%nface))
    line%nodes_face(1, 1, 1) = -1.0d0
    line%nodes_face(1, 1, 2) = 1.0d0

    if(ndim.eq.2) then
    ! Quadrilateral
    quad%nface = 4
    quad%order = order
    quad%nbasis = (order+1)*(order+1)

    quad%nvertex = 4
    allocate(quad%vertices(ndim, quad%nvertex))
    quad%vertices(:, 1) = (/-1.0d0, -1.0d0/)
    quad%vertices(:, 2) = (/ 1.0d0, -1.0d0/)
    quad%vertices(:, 3) = (/ 1.0d0,  1.0d0/)
    quad%vertices(:, 4) = (/-1.0d0,  1.0d0/)

    quad%nnode_cub = cub2D_npts
    allocate(quad%nodes_cub(quad%nnode_cub, 2))
    quad%nodes_cub = cubRS

    quad%nnode_face = quad_npts
    quad%nsurfnode = quad%nnode_face * quad%nface
    allocate(quad%nodes_face(quad%nnode_face, ndim, quad%nface))
    i = 1
    quad%nodes_face(1:quad_npts, 1, i) = quadR
    quad%nodes_face(1:quad_npts, 2, i) = -1.0d0
    i = 2                           
    quad%nodes_face(1:quad_npts, 1, i) = 1.0d0
    quad%nodes_face(1:quad_npts, 2, i) = quadR
    i = 3                           
    quad%nodes_face(1:quad_npts, 1, i) = quadR(size(quadR):1:-1)
    quad%nodes_face(1:quad_npts, 2, i) = 1.0d0
    i = 4                           
    quad%nodes_face(1:quad_npts, 1, i) = -1.0d0
    quad%nodes_face(1:quad_npts, 2, i) = quadR(size(quadR):1:-1)


    ! Triangle
!   ref_ele(3)%nface = 3
!   ref_ele(3)%nvertex = 3
!   ref_ele(3)%order = order
!   ref_ele(3)%nbasis = (order+1)*(order+2) / 2
!   allocate(ref_ele(3)%vertices(ndim, ref_ele(1)%nvertex))
!   ref_ele(3)%vertices(:, 1) = (/ 0.0d0,  0.0d0/)
!   ref_ele(3)%vertices(:, 2) = (/ 1.0d0,  0.0d0/)
!   ref_ele(3)%vertices(:, 3) = (/ 0.0d0,  1.0d0/)
!   allocate(ref_ele(3)%nodes_cub(2, ref_ele(3)%nbasis))

    end if !ndim.eq.2
end subroutine initialize_reference_elements


end module ref_ele_mod
