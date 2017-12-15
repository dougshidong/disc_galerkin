module cubature2D
use prec
use quadrature

integer :: icub, cub2D_npts
real(dp), allocatable :: cubRS(:,:), cubW(:)

contains

subroutine initialize_cub2D
integer :: i, j

call initialize_quad
cub2D_npts = quad_npts*quad_npts
allocate(cubRS(cub2D_npts, 2), cubW(cub2D_npts))
do j = 1, quad_npts
do i = 1, quad_npts
    cubRS(i+(j-1)*quad_npts, 1) = quadR(i)
    cubRS(i+(j-1)*quad_npts, 2) = quadR(j)
end do
end do
do j = 1, quad_npts    
do i = 1, quad_npts
    cubW(i+(j-1)*quad_npts) = quadW(i) * quadW(j)
end do
end do
!cubW = reshape(spread(quadW(1:quad_npts),dim=2,ncopies=quad_npts) &
!              *spread(quadW(1:quad_npts),dim=1,ncopies=quad_npts),&
!               (/quad_npts*quad_npts/), order = (/2 1/))

end subroutine initialize_cub2D

subroutine finalize_cub2D
deallocate(cubRS, cubW)
end subroutine finalize_cub2D

real(dp) function integrate2D(f)
implicit none
real(dp) :: f(cub2D_npts)
integrate2D = dot_product(cubW, f)
end function integrate2D


end module cubature2D
