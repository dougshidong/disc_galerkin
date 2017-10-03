module grid
    use glob
    use prec
    use matrices, only: Differentiation

    integer, allocatable :: EtoV(:,:)
    real(dp), allocatable :: vertices(:)
    real(dp), allocatable :: x(:,:)
    real(dp), allocatable :: Fx(:,:), Fscale(:,:)
    real(dp), allocatable :: Jac(:,:)
    real(dp), allocatable :: rx(:,:)
    real(dp), allocatable :: normals(:,:)

    contains
!   Read grid from grid.dat file
    subroutine readGrid
    implicit none
    integer :: gfile = 1
    integer :: i

    open(unit=gfile, file="grid.dat", action="read", form="formatted")
    read(*,*)
    read(*,*) order
    read(*,*)
    read(*,*) nele
    nvertex = nele+1
    read(*,*)

    allocate(EtoV(2,nele))
    allocate(vertices(nvertex))

    do i = 1, nvertex
        read(*,*) vertices(i)
    end do
    do i = 1, nele
        read(*,*) EtoV(1,i), EtoV(2,i)
    end do
    close(gfile)
    return

    end subroutine readGrid

    ! Generate an equally spaced grid
    subroutine genGrid(ref, xmin, xmax)
    use glob
    use prec
    use poly
    implicit none

    integer :: i, j
    real(dp) :: xmin, xmax, L
    real(dp) :: ref(:)

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
        Fx(1,i) = x(1,i)
        Fx(2,i) = x(nref,i)
        Fscale(1,i) = 1/Jac(1,i)
        Fscale(2,i) = 1/Jac(nref,i)
    end do

 
    return
    end subroutine gengrid
end module grid
