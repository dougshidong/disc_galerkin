module testcases
use prec
use glob
use grid
use quadrature
use cubature2D

contains
subroutine initialize_u
    use glob
    use prec
    implicit none
    real(dp) :: xshift, yshift
    integer :: i, iele, istate
    xshift = wavespeed(1)*finalTime;
    if(ndim.eq.2) then
        yshift = wavespeed(2)*finalTime;
    end if

    do iele = 1, nele
        elements(iele)%u = 0.0d0
        select case(icase)
            case(-1)
                wavespeed(1) = 0.5d0
                finalTime = 1e5
                xshift = wavespeed(1)*finalTime;
                elements(iele)%u(:,1) = sin(2.15d0*elements(iele)%x_cub(:,1)+0.23d0)
                elements(iele)%g(:,1) = 0.5d0*2.15d0*cos(2.15d0*elements(iele)%x_cub(:,1)+0.23d0)
                elements(iele)%u_exact(:,1) = sin(2.15d0*(elements(iele)%x_cub(:,1))+0.23d0)
            case(0)
                if(nstate.ne.1) print*, 'icase = 0 should have nstate = 1'
                if(ndim.eq.1) then
                    elements(iele)%u(:,1) = sin(elements(iele)%x_cub(:,1))   &
                                            + cos(elements(iele)%x_cub(:,1))
                    elements(iele)%u_exact(:,1) = sin(elements(iele)%x_cub(:,1)-xshift) &
                                                  + cos(elements(iele)%x_cub(:,1)-xshift)
                else if(ndim.eq.2) then
                    elements(iele)%u(:,1) = sin(elements(iele)%x_cub(:,1)) &
                                          + cos(elements(iele)%x_cub(:,1))  &
                                          + sin(elements(iele)%x_cub(:,2)) &
                                          + cos(elements(iele)%x_cub(:,2))
                    elements(iele)%u_exact(:,1) = sin(elements(iele)%x_cub(:,1)-xshift) &
                                                + cos(elements(iele)%x_cub(:,1)-xshift) &
                                                + sin(elements(iele)%x_cub(:,2)-yshift) &
                                                + cos(elements(iele)%x_cub(:,2)-yshift)
                end if
            case(1)
                do i = 1, elements(iele)%nnode_cub
                    if(ndim.eq.1) then
                        elements(iele)%u(i,1) = exp(-10.0d0*elements(iele)%x_cub(i,1)**2)
                        if(elements(iele)%u(i,1).lt.1e-15) elements(iele)%u(i,1) = 0.0d0
                        elements(iele)%u_exact(i,1) = exp(-10.0d0*(elements(iele)%x_cub(i,1)-xshift)**2)
                        if(elements(iele)%u_exact(i,1).lt.1e-15) elements(iele)%u_exact(i,1) = 0.0d0
                    else if(ndim.eq.2) then
                        elements(iele)%u(i,1) = exp(-5.0d0*elements(iele)%x_cub(i,1)**2 &
                                                    -5.0d0*elements(iele)%x_cub(i,2)**2)
                        if(elements(iele)%u(i,1).lt.1e-15) elements(iele)%u(i,1) = 0.0d0
                        elements(iele)%u_exact(i,1) = exp(-5.0d0*(elements(iele)%x_cub(i,1)-xshift)**2 &
                                                          -5.0d0*(elements(iele)%x_cub(i,2)-yshift)**2)
                        if(elements(iele)%u_exact(i,1).lt.1e-15) elements(iele)%u_exact(i,1) = 0.0d0
                    end if
                end do
            case(3)
                if(ndim.eq.1) then
                    elements(iele)%u(:,1) = sin(elements(iele)%x_cub(:,1))
                    elements(iele)%u_exact(:,1) = cos(elements(iele)%x_cub(:,1))
                else if(ndim.eq.2) then
                    ! u_solution = du/dx + du/dy
                    elements(iele)%u(:,1) = sin(elements(iele)%x_cub(:,1)) &
                                          + cos(elements(iele)%x_cub(:,2)) + 1.0d0
                    elements(iele)%u_exact(:,1) = cos(elements(iele)%x_cub(:,1)) &
                                                - sin(elements(iele)%x_cub(:,1))
                 end if
            case(4)
                finalTime = 1000000
                if(nstate.ne.1) print*, 'icase = 0 should have nstate = 1'
                if(ndim.eq.1) then
                    elements(iele)%u(:,1) = 0.0d0!cos(elements(iele)%x_cub(:,1))! * (abs(elements(iele)%x_cub(:,1) / 20))
                    elements(iele)%g(:,1) = -sin(elements(iele)%x_cub(:,1)) + cos(elements(iele)%x_cub(:,1))
                    elements(iele)%u_exact(:,1) = sin(elements(iele)%x_cub(:,1))+cos(elements(iele)%x_cub(:,1))
                else if(ndim.eq.2) then
                    elements(iele)%u(:,1) = 0.0d0
                    elements(iele)%u(:,1) = cos(elements(iele)%x_cub(:,1)) &
                                            - sin(elements(iele)%x_cub(:,2))
                    elements(iele)%g(:,istate) = -sin(elements(iele)%x_cub(:,1)) - cos(elements(iele)%x_cub(:,2))
                    elements(iele)%u_exact(:,1) = cos(elements(iele)%x_cub(:,1)) &
                                                - sin(elements(iele)%x_cub(:,2))
                end if
            case default
                print*, 'Invalid test case'
        end select
    end do
end subroutine initialize_u

subroutine error_evaluate(error_abs, error_rel)
use prec
use glob
real(dp) :: error_abs, error_rel, unorm
integer :: iele, istate

error_abs = 0.0d0
error_rel = 0.0d0
unorm = 0.0d0
do iele = 1, nele
    do istate = 1, nstate
        elements(iele)%u_error(:,istate) = (elements(iele)%u(:,istate) - elements(iele)%u_exact(:,istate))
        if(ndim.eq.1) error_abs = error_abs + abs(integrate(elements(iele)%u_error(:,istate)**2))
        if(ndim.eq.2) error_abs = error_abs + abs(integrate2D(elements(iele)%u_error(:,istate)**2))
        if(ndim.eq.1) unorm = unorm + abs(integrate(elements(iele)%u_exact(:,istate)**2))
        if(ndim.eq.2) unorm = unorm + abs(integrate2D(elements(iele)%u_exact(:,istate)**2))
        !print*, iele, integrate(elements(iele)%u_error(:,istate)), integrate(elements(iele)%u_exact(:,istate))
    end do
end do
error_abs = sqrt(error_abs)
unorm = sqrt(unorm)
error_rel = error_abs / unorm

end subroutine error_evaluate


end module testcases
