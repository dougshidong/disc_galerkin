module matrices
    use prec
    use poly
    ! Mass = inv(V*V^T) = inv(V^T)*inv(V)
    ! MassInv = V*V^T
    real(dp), allocatable :: Mass(:,:), MassInv(:,:)
    real(dp), allocatable :: Stiff(:,:)
    ! Lift matrix contains the Lagrange function evaluated at the surface points
    real(dp), allocatable :: Lift(:,:)
    ! Vandermonde matrix comes from the change of basis
    ! V^T * lagrange = P
    ! u_nodal = V * umodal
    real(dp), allocatable :: Vander(:,:), VanderInv(:,:)
    real(dp), allocatable :: VanderSurf(:,:)
    ! Differentiation matrix defines u'=Dr*u. Dr = Minv*S
    real(dp), allocatable :: Differentiation(:,:)

    contains
    
    subroutine buildVandermonde(ref)
        use glob
        implicit none
        real(dp), intent(in) :: ref(nref)
        real(dp) :: tempV(nref, nref)
        real(dp) :: polyVal(nref)
        integer :: i, j
        integer :: info, ipiv(nref)

        ! Evaluate Vandermonde
        do j = 0, order
            call poly_eval(polytype, j, ref, polyVal) 
            Vander(1:nref,j+1) = polyVal
        end do
        ! Evaluate inverse Vandermonde
        tempV = Vander
        VanderInv = 0.0d0
        do i = 1, nref
            VanderInv(i,i) = 1.0d0
        end do
        call dgesv(nref, nref, tempV, nref, ipiv, VanderInv, nref, info )
        return
    end subroutine buildVandermonde

    subroutine buildMass
        use glob
        use quadrature
        implicit none
        real(dp) :: ipolyVal(quad_npts), jpolyVal(quad_npts)
        real(dp) :: tempMass(nref,nref)
        integer :: iref, jref
        integer :: info, ipiv(nref)

        ! Evaluate Mass matrix
        do iref = 1, nref
        call poly_eval(polytype, iref-1, xquad, ipolyVal) 
        do jref = 1, nref
            call poly_eval(polytype, jref-1, xquad, jpolyVal) 
            tempMass(iref,jref) = integrate(ipolyVal*jpolyVal)
        end do
        end do
        ! Change of basis to Lagrange polynomial
        Mass = matmul(matmul(transpose(VanderInv), tempMass), VanderInv)
        ! Evaluate inverse Mass matrix
        tempMass = Mass
        MassInv = 0.0d0
        do iref = 1, nref
            MassInv(iref,iref) = 1.0d0
        end do
        call dgesv(nref, nref, tempMass, nref, ipiv, MassInv, nref, info )
        if(info.ne.0) print *, 'Failed to evaluate inverse Mass'
        return
    end subroutine buildMass

    subroutine buildStiffness
        use glob
        use quadrature
        implicit none
        real(dp) :: polyGrad(quad_npts), polyVal(quad_npts)
        real(dp) :: tempStiff(nref,nref)
        integer :: iref, jref

        do iref = 1, nref
        call polyGrad_eval(polytype, iref-1, xquad, polyGrad) 
        do jref = 1, nref
            call poly_eval(polytype, jref-1, xquad, polyVal) 
            ! Weak form dl*l
            !tempStiff(iref,jref) = integrate(polyGrad*polyVal)
            ! Strong form l*dl
            tempStiff(jref,iref) = integrate(polyGrad*polyVal)
        end do
        end do
        !Stiff = matmul(matmul(VanderInv, tempStiff), transpose(VanderInv))
        Stiff = matmul(matmul(transpose(VanderInv), tempStiff), VanderInv)
    end subroutine buildStiffness

    subroutine buildDifferentiation
    Differentiation = matmul(MassInv,Stiff)
    end subroutine buildDifferentiation

    subroutine buildLift
    implicit none
    integer :: iref
    real(dp) :: refSurf(2)
    real(dp) :: polyValSurf(2)
    real(dp) :: basisSurf(nref,2)
    refSurf(1) = -1.0d0
    refSurf(2) = 1.0d0
    do iref = 1, nref
        call poly_eval(polytype, iref-1, refSurf, polyValSurf) 
        basisSurf(iref,:) = polyValSurf
    end do
    VanderSurf = matmul(transpose(VanderInv),basisSurf)
    Lift = matmul(MassInv,VanderSurf)
    end subroutine buildLift

    subroutine allocate_matrices
        use glob
        implicit none
        allocate(Lift(nref, nface*nfpoint))
        allocate(Vander(nref,nref), VanderInv(nref,nref))
        allocate(VanderSurf(nref,nface*nfpoint))
        allocate(Mass(nref,nref), MassInv(nref,nref))
        allocate(Stiff(nref,nref))
        allocate(Differentiation(nref,nref))
        return
    end subroutine allocate_matrices

    subroutine deallocate_matrices
        deallocate(Lift)
        deallocate(Vander, VanderInv)
        deallocate(VanderSurf)
        deallocate(Mass, MassInv)
        deallocate(Stiff)
        deallocate(Differentiation)
        return
    end subroutine deallocate_matrices

end module matrices
