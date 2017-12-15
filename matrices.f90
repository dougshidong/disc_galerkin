module matrices
    use prec
    use glob, only: inodal
    use poly
    ! Mass = inv(V*V^T) = inv(V^T)*inv(V)
    ! MassInv = V*V^T
    real(dp), allocatable :: Mass(:,:), MassInv(:,:)
    real(dp), allocatable :: Sr(:,:), Ss(:,:)
    ! Lift matrix contains the Lagrange function evaluated at the surface points
    real(dp), allocatable :: Lift(:,:,:)
    ! Vandermonde matrix comes from the change of basis
    ! V^T * lagrange = P
    ! u_nodal = V * umodal
    real(dp), allocatable :: Vander(:,:), VanderInv(:,:)
    real(dp), allocatable :: VanderSurf(:,:,:)
    real(dp), allocatable :: VanderGrad(:,:,:)
    ! Dr = Minv*Sr
    real(dp), allocatable :: Dr(:,:), Ds(:,:)

    contains
    
    subroutine buildVandermondeGrad
        use glob
        use ref_ele_mod
        implicit none
        real(dp) :: polyGrad(nbasis), polyVal(nbasis)
        integer :: kdir, ibasis, jbasis

        ! Evaluate Vandermonde
        do kdir = 1, ndim
        do ibasis = 1, nbasis
            call polyGrad_eval(polytype, ibasis, order, kdir, ndim, ref_ele(ndim)%p%nnode_cub, ref_ele(ndim)%p%nodes_cub, polyGrad(:)) 
            VanderGrad(1:nbasis,ibasis,kdir) = polyGrad
        end do
        end do
        return
    end subroutine buildVandermondeGrad

    subroutine VandermondeGrad(VG, refnodes, nnode, porder)
        use glob
        implicit none
        integer :: nnode, porder
        real(dp) :: VG(nnode, nbasis, ndim), refnodes(nnode,ndim)
        real(dp) :: polyGrad(nnode)
        integer :: kdir, ibasis

        ! Evaluate Vandermonde Grad
        do kdir = 1, ndim
        do ibasis = 1, nbasis
            call polyGrad_eval(polytype, ibasis, porder, kdir, ndim, nnode, refnodes, polyGrad) 
            VG(1:nnode, ibasis, kdir) = polyGrad
        end do
        end do
        return
    end subroutine VandermondeGrad
    subroutine Vandermonde(V, refnodes, nnode, porder)
        use glob
        implicit none
        integer :: nnode, porder
        real(dp) :: V(nbasis, nnode), refnodes(nnode,ndim)
        real(dp) :: polyVal(nnode)
        integer :: ibasis

        ! Evaluate Vandermonde
        do ibasis = 1, nbasis
            call poly_eval(polytype, ibasis, porder, ndim, nnode, refnodes, polyVal)
            V(1:nnode,ibasis) = polyVal
        end do
        return
    end subroutine Vandermonde

    subroutine buildVandermonde
        use glob
        use ref_ele_mod
        implicit none
        real(dp) :: tempV(nbasis, nbasis)
        real(dp) :: polyVal(nbasis)
        integer :: info, ipiv(nbasis)
        integer :: ibasis

        ! Evaluate Vandermonde
        do ibasis = 1, nbasis
            call poly_eval(polytype, ibasis, order, ndim, ref_ele(ndim)%p%nnode_cub, ref_ele(ndim)%p%nodes_cub, polyVal) 
            Vander(1:nbasis,ibasis) = polyVal
        end do
        ! Evaluate inverse Vandermonde
        call invertMatrix(Vander, VanderInv, nbasis, info)
        if(info.ne.0) print *, 'Failed to evaluate inverse Vandermonde'
        return
    end subroutine buildVandermonde

    subroutine invertMatrix(A, Ainv, n, info)
        implicit none
        integer :: i, n, info, ipiv(n)
        real(dp), dimension(n,n) :: A, Ainv, tempA

        tempA = A
        Ainv = 0.0d0
        do i = 1, n
            Ainv(i, i) = 1.0d0
        end do
        call dgesv(n, n, tempA, n, ipiv, Ainv, n, info )
        return
    end subroutine

    subroutine buildMass
        use glob
        use quadrature
        use cubature2D
        use ref_ele_mod
        implicit none
        real(dp) :: ipolyVal(nbasis), jpolyVal(nbasis)
        real(dp) :: tempMass(nbasis,nbasis)
        integer :: ibasis, jbasis
        integer :: info, ipiv(nbasis)

        ! Evaluate Mass matrix
        do ibasis = 1, nbasis
        call poly_eval(polytype, ibasis, order, ndim, ref_ele(ndim)%p%nnode_cub, ref_ele(ndim)%p%nodes_cub, ipolyVal) 
        do jbasis = 1, nbasis
            call poly_eval(polytype, jbasis, order, ndim, ref_ele(ndim)%p%nnode_cub, ref_ele(ndim)%p%nodes_cub, jpolyVal) 
            if(ndim.eq.1) tempMass(ibasis,jbasis) = integrate(ipolyVal*jpolyVal)
            if(ndim.eq.2) tempMass(ibasis,jbasis) = integrate2D(ipolyVal*jpolyVal)
        end do
        end do
        ! Change of basis to Lagrange polynomial
        Mass = tempMass
        if(inodal) Mass = matmul(matmul(transpose(VanderInv), tempMass), VanderInv)
        ! Evaluate inverse Mass matrix
        tempMass = Mass
        MassInv = 0.0d0
        do ibasis = 1, nbasis
            MassInv(ibasis,ibasis) = 1.0d0
        end do
        call dgesv(nbasis, nbasis, tempMass, nbasis, ipiv, MassInv, nbasis, info )
        if(info.ne.0) print *, 'Failed to evaluate inverse Mass'
        return
    end subroutine buildMass

    subroutine buildStiffness
        use glob
        use quadrature
        use cubature2D
        use ref_ele_mod
        implicit none
        real(dp) :: polyGrad(nbasis,ndim), polyVal(nbasis)
        real(dp) :: tempStiff(nbasis,nbasis, ndim)
        integer :: idir, ibasis, jbasis

        do idir = 1, ndim
            do ibasis = 1, nbasis
                call polyGrad_eval(polytype, ibasis, order, idir, ndim, ref_ele(ndim)%p%nnode_cub, ref_ele(ndim)%p%nodes_cub, polyGrad(:, idir)) 
                do jbasis = 1, nbasis
                    call poly_eval(polytype, jbasis, order, ndim, ref_ele(ndim)%p%nnode_cub, ref_ele(ndim)%p%nodes_cub, polyVal) 
                    ! Strong form l*dl
                    if(ndim.eq.1) tempStiff(jbasis,ibasis,idir) = integrate(polyGrad(:, idir)*polyVal)
                    if(ndim.eq.2) tempStiff(jbasis,ibasis,idir) = integrate2D(polyGrad(:, idir)*polyVal)
                end do
            end do
        end do
        if(.not.inodal) Sr = tempStiff(:,:,1)
        if(inodal) Sr = matmul(matmul(transpose(VanderInv), tempStiff(:,:,1)), VanderInv)
        if(ndim.ge.2) then
            if(.not.inodal) Ss = tempStiff(:,:,2)
            if(inodal) Ss = matmul(matmul(transpose(VanderInv), tempStiff(:,:,2)), VanderInv)
        end if
    end subroutine buildStiffness

    subroutine buildDifferentiation
        Dr = matmul(MassInv,Sr)
        if(ndim.ge.2) Ds = matmul(MassInv,Ss)
!       Dr = matmul(Mass,Sr)
!       if(ndim.ge.2) Ds = matmul(Mass,Ss)
    end subroutine buildDifferentiation

    subroutine buildLift
        use ref_ele_mod
        implicit none
        integer :: ibasis, iface
        real(dp), allocatable :: polyValSurf(:)
        real(dp), allocatable :: basisSurf(:, :, :)
        allocate(polyValSurf(ref_ele(ndim)%p%nnode_face))
        allocate(basisSurf(nbasis, ref_ele(ndim)%p%nnode_face, ref_ele(ndim)%p%nface))
        do ibasis = 1, nbasis
        do iface = 1, ref_ele(ndim)%p%nface
            call poly_eval(polytype, ibasis, order, ndim, ref_ele(ndim)%p%nnode_face, ref_ele(ndim)%p%nodes_face(:,:,iface), polyValSurf) 
            basisSurf(ibasis,:,iface) = polyValSurf
        end do
        end do
        VanderSurf = basisSurf
        do iface = 1, ref_ele(ndim)%p%nface
            if(inodal) VanderSurf(:, :, iface) = matmul(transpose(VanderInv), basisSurf(:, :, iface))
            VanderSurf(:, :, iface) = matmul(transpose(VanderInv), basisSurf(:, :, iface))
            do ibasis = 1, nbasis
                if(ndim.eq.1) Lift(ibasis, : , iface) = VanderSurf(ibasis, :, iface)
                if(ndim.eq.2) Lift(ibasis, : , iface) = VanderSurf(ibasis, :, iface) * quadW
            end do
            !Lift(:, : , iface) = matmul(MassInv, VanderSurf(:, :, iface))
            Lift(:, : , iface) = matmul(MassInv, Lift(:, : , iface))
        end do
    end subroutine buildLift

    subroutine allocate_matrices
        use glob
        use ref_ele_mod
        implicit none
!       integer :: nface, nfpoint
!       ! Assume only 1 type of element
!       if(ndim.eq.1) then
!           nface = 2
!           nfpoint = 1
!       else if(ndim.eq.2) then !quad
!           nface = 4
!           nfpoint = quad_npts
!       end if
        allocate(Lift(nbasis, ref_ele(ndim)%p%nnode_face, ref_ele(ndim)%p%nface))
        allocate(VanderSurf(nbasis, ref_ele(ndim)%p%nnode_face, ref_ele(ndim)%p%nface))

        allocate(Vander(nbasis,nbasis), VanderInv(nbasis,nbasis))
        allocate(VanderGrad(nbasis, nbasis, ndim))

        allocate(Mass(nbasis,nbasis), MassInv(nbasis,nbasis))
        allocate(Sr(nbasis,nbasis), Dr(nbasis,nbasis))
        if(ndim.ge.2) allocate(Ss(nbasis,nbasis), Ds(nbasis,nbasis))
        return
    end subroutine allocate_matrices

    subroutine deallocate_matrices
        deallocate(Lift)
        deallocate(Vander, VanderInv)
        deallocate(VanderSurf)
        deallocate(Mass, MassInv)
        deallocate(Sr, Dr)
        if(ndim.ge.2) deallocate(Ss, Ds)
        return
    end subroutine deallocate_matrices

end module matrices
