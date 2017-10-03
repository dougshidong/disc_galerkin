module matrices
    use prec
    ! Lift matrix contains the Lagrange function evaluated at the surface points
    real(dp), allocatable :: Lift(:,:)
    ! Vandermonde matrix comes from the change of basis
    real(dp), allocatable :: Vandermonde(:,:), VandermondeInv(:,:)
    ! Differentiation matrix defines u'=Dr*u. Dr = Minv*S
    real(dp), allocatable :: Differentiation(:,:)

    contains

    subroutine allocate_matrices
        use glob
        implicit none
        allocate(Lift(nref, nface*nfpoint))
        allocate(Vandermonde(nref,nref), VandermondeInv(nref,nref))
        allocate(Differentiation(nref,nref))
        return
    end subroutine allocate_matrices

    subroutine deallocate_matrices
        deallocate(Lift)
        deallocate(Vandermonde, VandermondeInv)
        deallocate(Differentiation)
        return
    end subroutine deallocate_matrices

end module matrices
