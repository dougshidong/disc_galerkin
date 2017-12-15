module face_mod
use glob

type face
    integer :: facetype
    integer :: ele_conn(2)
    integer :: ele_face(2)
    logical :: reversed
    integer :: nvertex
    integer, allocatable :: vertices(:)
end type face

! node_face_list(nvertex per face, number of faces, facetype)
! Gives a list of nodes for face number X
integer :: node_face_list(2,4,2)

contains
subroutine initialize_face_node
implicit none
node_face_list = 0
! 1D line element
node_face_list(1,1,1) = 1
node_face_list(1,2,1) = 2
! 2D quad element
node_face_list(1:2,1,2) = (/1, 2/)
node_face_list(1:2,2,2) = (/2, 3/)
node_face_list(1:2,3,2) = (/3, 4/)
node_face_list(1:2,4,2) = (/4, 1/)

return
end subroutine initialize_face_node

end module face_mod
