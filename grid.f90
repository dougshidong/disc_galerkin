module grid
    use glob
    use prec
    use matrices, only: Differentiation

    integer, allocatable :: EtoV(:,:)
    real(dp), allocatable :: vertices(:)
    real(dp), allocatable :: x(:,:)
    real(dp), allocatable :: Fscale(:,:)
    real(dp), allocatable :: Jac(:,:)
    real(dp), allocatable :: rx(:,:)
    real(dp), allocatable :: normals(:,:)

    contains
    subroutine finalize_grid
    deallocate(EtoV)
    deallocate(vertices)
    deallocate(x, rx)
    deallocate(normals)
    deallocate(Fscale)
    end subroutine finalize_grid
    ! Generate an equally spaced grid of elements
    subroutine genGrid(ref, xmin, xmax)
    use glob
    use prec
    use poly
    implicit none

    integer :: i, j
    real(dp) :: xmin, xmax, L
    real(dp) :: ref(:)

    ! Generate 1D grid and map reference elements
    allocate(EtoV(2,nele))
    allocate(vertices(nvertex))
    allocate(x(nref,nele), rx(nref,nele))
    allocate(normals(nfpoint*nface,nele))
    allocate(Fscale(2,nele))
    select case(icase)
        case(0)
            xmin = 0.0d0
            xmax = PI
        case(1)
            xmin = -2.0d0
            xmax = 2.0d0
        case default
            print *, 'Invalid case'
    end select

    L = xmax - xmin

    ! Create nodes at element ends
    do i = 1, nvertex
        vertices(i) = L*(i-1)/(nvertex-1) + xmin
    end do

    ! Create element-to-node connectivity
    do i = 1, nele
        EtoV(1,i) = i
        EtoV(2,i) = i+1
    end do

    ! Create array X of (Np x K) that contains the physical coordinates
    do i = 1, nref
    do j = 1, nele
        x(i,j) = vertices(EtoV(1,j)) &
            + 0.5 * (ref(i)+1.0d0) * (vertices(EtoV(2,j)) - vertices(EtoV(1,j)))
    end do
    end do

    ! Create surface normals
    normals(1,:) = -1.0d0
    normals(2,:) = 1.0d0

    Jac = matmul(Differentiation,x)
    ! rx is the geometric factor (2 / h^k)
    rx = 1.0d0 / Jac
    do i = 1, nele
        Fscale(1,i) = 1/Jac(1,i)
        Fscale(2,i) = 1/Jac(nref,i)
    end do

 
    return
    end subroutine gengrid
end module grid
