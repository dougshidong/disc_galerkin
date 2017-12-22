module grid
    use glob
    use prec
!   use matrices!, only: Dr, Vander
    use quadrature, only: quad_npts, quadR
    use element_mod
    use face_mod

    integer :: nface_global_max
    integer :: nface_element, npts_face, nvertex_element
    integer :: nface_global, nvertex_global
    integer :: nivertex, njvertex
    integer, allocatable :: EtoV(:,:)
    real(dp), allocatable :: global_vertices(:,:)
    type(face), allocatable :: global_faces(:)
    type(element), allocatable :: elements(:)

    real(dp) :: xmin, xmax, Lx
    real(dp) :: ymin, ymax, Ly

    contains
    subroutine finalize_grid
        deallocate(EtoV)
        deallocate(global_faces)
        deallocate(global_vertices)

        ! Allocatable components are deallocated automatically, but not pointers
        ! Deallocating entire grid will deallocate all its components
        deallocate(elements)
    end subroutine finalize_grid

    ! Generate an equally spaced grid of elements
    subroutine genGrid2D
        use glob
        use prec
        !use poly
        implicit none

        integer :: i, j
        integer :: iele, iele1, iele2
        integer :: iface
        integer :: iface_element1, iface_element2, iface_global
        integer :: ivertex
    !   real(dp) :: ref(ndim, nbasis)

    !   type(ref) :: ref
        integer :: vertices1(ndim), vertices2(ndim)
        integer :: vertices1_sorted(ndim), vertices2_sorted(ndim)
        logical :: found_face

        if(ndim.eq.1) then
            nface_element = 2
            npts_face = 1
            nvertex_element = 2
        else if(ndim.eq.2) then
            nface_element = 4 ! Quadrilaterals
            npts_face = quad_npts
            nvertex_element = nface_element
        end if

        if(ndim.eq.1) njvertex = 1
        nvertex_global = nivertex*njvertex

        if(ndim.eq.1) nele = (nivertex-1)
        if(ndim.eq.2) nele = (njvertex-1)*(nivertex-1)

        allocate(global_vertices(nvertex_global, ndim))
    !   allocate(normals(npts_face,nface_element,nele))
    !   allocate(Fscale(npts_face,nface_element,nele))

        select case(icase)
            case(-1)
                xmin = -1.0d0; xmax = 1.0d0;
                !print *, 'Sine wave case'
            case(0)
                xmin = -1.0d0*PI; xmax = 1.0d0*PI;
                ymin = -1.0d0*PI; ymax = 1.0d0*PI;
                !print *, 'Sine wave case'
            case(1)
                xmin = -3.0d0; xmax = 3.0d0;
                ymin = -3.0d0; ymax = 3.0d0;
                !print *, 'Bump wave case'
            case(3)
                xmin = -1.0d0; xmax = 1.0d0;
                ymin = -1.0d0; ymax = 1.0d0;
                !print *, 'Sine differentiation case'
            case(4)
                xmin = -1.0d0*PI; xmax = 1.0d0*PI;
                ymin = -1.0d0*PI; ymax = 1.0d0*PI;
                !print *, 'Source term sine case'
            case default
                xmin = -0.9d0; xmax = 1.0d0;
                ymin = -0.9d0; ymax = 1.0d0;
                !print *, 'Default Case'
        end select
        Lx = xmax - xmin
        Ly = ymax - ymin

        ! Create nodes at element global_vertices
        ivertex = 0
        do j = 1, njvertex
        do i = 1, nivertex
            ivertex = ivertex + 1
            ! Line grid
            if(ndim.ge.1) global_vertices(ivertex, 1) = Lx*(i-1.0d0)/(nivertex-1) + xmin
            if(ndim.ge.2) then
                ! Square grid
                global_vertices(ivertex,1) = Lx*(i-1.0d0)/(nivertex-1.0d0) + xmin
                global_vertices(ivertex,2) = Ly*(j-1.0d0)/(njvertex-1.0d0) + ymin
                ! Parallelogram grid
                !global_vertices(ivertex,1) = Lx*(i-1.0d0)/(nivertex-1) + xmin
                !global_vertices(ivertex,2) = Ly*(j-1.0d0)/(njvertex-1.0d0) + ymin + global_vertices(ivertex,1)/Lx
            end if
        end do
        end do

        allocate(elements(nele))
        do iele = 1, nele
            if(ndim.eq.1) elements(iele)%typ = 1
            if(ndim.eq.2) elements(iele)%typ = 2
            elements(iele) = init_ele(elements(iele)%typ)
        end do

        allocate(EtoV(nvertex_element,nele))
        call setEtoV

        ! Setup elements vertices, map the reference elements and evaluate the Jacobians
        do iele = 1, nele
            call element_set_vertices(elements(iele), EtoV(:,iele))
            call map_ref_to_phys(elements(iele))
            !call element_jacobian_eval(elements(iele))
            call GeometricFactors(elements(iele)%x_cub, elements(iele)%nnode_cub, elements(iele), &
                elements(iele)%drdx_cub, elements(iele)%dsdx_cub, &
                elements(iele)%drdy_cub, elements(iele)%dsdy_cub, elements(iele)%Jac_cub)
            do iface = 1, elements(iele)%nface
                call GeometricFactors(elements(iele)%x_face(:,:,iface), elements(iele)%nnode_face, elements(iele), &
                    elements(iele)%drdx_face(:,iface), elements(iele)%dsdx_face(:,iface), &
                    elements(iele)%drdy_face(:,iface), elements(iele)%dsdy_face(:,iface), elements(iele)%Jac_face(:,iface))
                ! hardcoded normals of 1.0
                !elements(iele)%Fscale(:,iface) = 1.0d0 / elements(iele)%Jac_face(:,iface)
                elements(iele)%Fscale(:,iface) = elements(iele)%drdx_face(:,iface)! / elements(iele)%Jac_face(:,iface)
            end do
        end do

        ! Create global faces list
        nface_global_max = 4*nele
        allocate(global_faces(nface_global_max))
        if(ndim.eq.1) then ! Point
            global_faces%facetype = 1
            global_faces%nvertex = 1
        else if(ndim.eq.2) then ! Line
            global_faces%facetype = 2
            global_faces%nvertex = 2
        end if
        do iface = 1, nface_global_max
            global_faces(iface)%ele_conn(:) = -1
            allocate(global_faces(iface)%vertices(global_faces(iface)%nvertex))
        end do

        call initialize_face_node
        ! Setup global faces list
        iface_global = 1
        do iele1 = 1, nele
        do iface_element1 = 1, ref_ele(elements(iele1)%typ)%p%nface

            if(elements(iele1)%face_pointer(iface_element1).ne.-1) cycle
            found_face = .false.
            ! New global face found
            elements(iele1)%face_pointer(iface_element1) = iface_global
            global_faces(iface_global)%ele_conn(1) = iele1
            global_faces(iface_global)%ele_face(1) = iface_element1

            ! Set face vertices
            do ivertex = 1, global_faces(iface_global)%nvertex
                global_faces(iface_global)%vertices(ivertex) &
                    = EtoV(node_face_list(ivertex,iface_element1,global_faces(iface_global)%facetype), iele1)
            end do
        
            ! Set vertices to compare
            do ivertex = 1, global_faces(iface_global)%nvertex
                vertices1(ivertex) = global_faces(iface_global)%vertices(ivertex)
            end do
            vertices1_sorted = vertices1
            if(ndim.ge.2) call quicksort(vertices1_sorted, 1, 2)

            ! Find a corresponding face with same vertices
            do iele2 = iele1+1, nele
            do iface_element2 = 1, ref_ele(elements(iele2)%typ)%p%nface
                if(found_face) cycle
                ! Get comparison vertices
                do ivertex = 1,  global_faces(iface_global)%nvertex
                    vertices2(ivertex) = EtoV(node_face_list(ivertex,iface_element2,global_faces(iface_global)%facetype), iele2)
                end do
                vertices2_sorted = vertices2
                if(ndim.ge.2) call quicksort(vertices2_sorted, 1, 2)

                ! Compare faces
                if(allequal(vertices1_sorted, vertices2_sorted)) then
                    found_face = .true.
                    elements(iele2)%face_pointer(iface_element2) = iface_global
                    global_faces(iface_global)%ele_conn(2) = iele2
                    global_faces(iface_global)%ele_face(2) = iface_element2
                    global_faces(iface_global)%reversed = .not.(allequal(vertices1, vertices2))
                end if
            end do
            end do

            if(.not.found_face) cycle
            iface_global = iface_global + 1
        end do
        end do ! Element loop
        nface_global = iface_global - 1

        ! Setup periodic grid (hard-coded)
        if(ndim.eq.1) then
            nface_global = nface_global+1
            iface_global = nface_global
            global_faces(iface_global)%ele_conn(1) = 1
            global_faces(iface_global)%ele_conn(2) = nele
            global_faces(iface_global)%ele_face(1) = 1
            global_faces(iface_global)%ele_face(2) = 2
            global_faces(iface_global)%reversed = .false.
            elements(1)%face_pointer(1) = iface_global
            elements(nele)%face_pointer(2) = iface_global

        else if(ndim.eq.2) then
            do iele1 = 1, nivertex-1
                nface_global = nface_global+1
                iface_global = nface_global
                iele2 = nele - (nivertex-1) + iele1
                global_faces(iface_global)%ele_conn(1) = iele1
                global_faces(iface_global)%ele_conn(2) = iele2
                global_faces(iface_global)%ele_face(1) = 1
                global_faces(iface_global)%ele_face(2) = 3
                global_faces(iface_global)%reversed = .true.
                elements(iele1)%face_pointer(1) = iface_global
                elements(iele2)%face_pointer(3) = iface_global
            end do

            do i = 1, njvertex-1
                nface_global = nface_global+1
                iface_global = nface_global
                iele1 = (i-1) * (nivertex-1) + 1
                iele2 = i*(nivertex-1)
                global_faces(iface_global)%ele_conn(1) = iele1
                global_faces(iface_global)%ele_conn(2) = iele2
                global_faces(iface_global)%ele_face(1) = 4
                global_faces(iface_global)%ele_face(2) = 2
                global_faces(iface_global)%reversed = .true.
                elements(iele1)%face_pointer(4) = iface_global
                elements(iele2)%face_pointer(2) = iface_global
            end do
        end if


    !   dxdr = matmul(Dr,x)
    !   if(.not. inodal) dxdr = matmul(Vander,matmul(Dr,matmul(VanderInv,x)))
    !   ! drdx_cub is the geometric factor (2 / h^k)
    !   drdx_cub = 1.0d0 / dxdr
    !   do i = 1, niele
    !       Fscale(1,i) = 1/dxdr(1,i)
    !       Fscale(2,i) = 1/dxdr(nbasis,i)
    !   end do
    !   print*,dxdr
    !   print*,
    !   print*,drdx_cub

        return
    end subroutine genGrid2D

    subroutine setEtoV
        implicit none
        integer :: ivertex, iele, iv, jv
        ! Create element-to-node connectivity
        ! http://gmsh.info/doc/texinfo/gmsh.html#Node-ordering

        if(ndim.eq.1) then
            ! One-dimensional structured grid
            do iele = 1, nivertex-1
                EtoV(1,iele) = iele
                EtoV(2,iele) = iele+1
            end do
        else if(ndim.eq.2) then
            ! Special case of structured quadrilateral grid
            ivertex = 0; iele = 0
            do jv = 1, njvertex; do iv = 1, nivertex
                ivertex = ivertex + 1
                if(iv.eq.nivertex .or. jv.eq.njvertex) cycle
                iele = iele + 1
                EtoV(1,iele) = ivertex
            end do; end do
            ivertex = 0; iele = 0
            do jv = 1, njvertex; do iv = 1, nivertex
                ivertex = ivertex + 1
                if(iv.eq.1 .or. jv.eq.njvertex) cycle
                iele = iele + 1
                EtoV(2,iele) = ivertex
            end do; end do
            ivertex = 0; iele = 0
            do jv = 1, njvertex; do iv = 1, nivertex
                ivertex = ivertex + 1
                if(iv.eq.1 .or. jv.eq.1) cycle
                iele = iele + 1
                EtoV(3,iele) = ivertex
            end do; end do
            ivertex = 0; iele = 0
            do jv = 1, njvertex; do iv = 1, nivertex
                ivertex = ivertex + 1
                if(iv.eq.nivertex .or. jv.eq.1) cycle
                iele = iele + 1
                EtoV(4,iele) = ivertex
            end do; end do
        end if
    end subroutine setEtoV

    logical function allequal(a, b)
        implicit none
        integer, dimension(:) :: a, b
        integer :: i
        allequal = ( size(a) == size(b) )
        if(allequal) then
        do i = 1, size(a)
            allequal = ( a(i) .eq. b(i) )
            if(.not.allequal) exit
        end do
        end if
        return
    end function

    recursive subroutine quicksort(a, first, last)
        implicit none
        integer :: a(*), x, t
        integer :: first, last
        integer :: i, j

        x = a( (first+last) / 2 )
        i = first
        j = last
        do
         do while (a(i) < x)
            i=i+1
         end do
         do while (x < a(j))
            j=j-1
         end do
         if (i >= j) exit
         t = a(i);  a(i) = a(j);  a(j) = t
         i=i+1
         j=j-1
        end do
        if (first < i-1) call quicksort(a, first, i-1)
        if (j+1 < last)  call quicksort(a, j+1, last)
    end subroutine quicksort


    subroutine element_set_vertices(ele, vertices)
        implicit none
        type(element) :: ele
        integer :: vertices(ref_ele(ele%typ)%p%nvertex)
        ele%vertices(:) = vertices
        return
    end subroutine element_set_vertices

    subroutine map_ref_to_phys(ele)
        implicit none
        type(element) :: ele
        type(ref) :: refe
        real(dp) :: dref(3), x0(3), dx(3), xrel(3)
        real(dp) :: quadref1D(2), quadphys(2,2,2)
        integer :: i, j, edim
        real(dp), allocatable :: ref_v(:,:), phys_v(:,:), ref_n(:,:), phys_n(:,:)
        ! Map reference element onto physical element
        ! Uniform grid
        refe = ref_ele(ele%typ)%p
        if(ele%typ.eq.1) then ! Line
            edim = 1
            dref(1:edim) = refe%vertices(2,1:edim) - refe%vertices(1,1:edim)
            x0(1:edim)   = global_vertices(ele%vertices(1),1:edim)
            dx(1:edim)   = global_vertices(ele%vertices(2),1:edim) - global_vertices(ele%vertices(1),1:edim)
            do i = 1, nbasis
                xrel(1:edim) = refe%nodes_cub(i, 1:edim) - refe%vertices(1, 1:edim)
                ele%x_cub(i, 1:edim) = x0(1:edim) + xrel(1:edim) /dref(1:edim) * dx(1:edim)
            end do
            ele%x_face(1,1:edim,1) = refe%vertices(1,1:edim)
            ele%x_face(1,1:edim,2) = refe%vertices(2,1:edim)
        else if(ele%typ.eq.2) then
            edim = 2
            allocate(phys_n(refe%nnode_cub, edim))
            quadref1D(1) = -1.0d0; quadref1D(2) = 1.0d0;
            !quadref1D(1) = 0.0d0; quadref1D(2) = 1.0d0;
            quadphys(1,1,1:edim) = global_vertices(ele%vertices(1), 1:edim)
            quadphys(2,1,1:edim) = global_vertices(ele%vertices(2), 1:edim)
            quadphys(2,2,1:edim) = global_vertices(ele%vertices(3), 1:edim)
            quadphys(1,2,1:edim) = global_vertices(ele%vertices(4), 1:edim)
            call lagrange2Dquad(quadref1D, quadref1D, quadphys, 2, 2, &
                                refe%nodes_cub(:, 1), refe%nodes_cub(:,2), phys_n, refe%nnode_cub)
            do i = 1, refe%nnode_cub
                ele%x_cub(i, 1:edim) = phys_n(i, 1:edim)
            end do

            deallocate(phys_n)
            allocate(phys_n(refe%nnode_face, edim))
            do i = 1, refe%nface
                call lagrange2Dquad(quadref1D, quadref1D, quadphys, 2, 2, &
                                    refe%nodes_face(:,1,i), refe%nodes_face(:,2,i), phys_n, refe%nnode_face)
                do j = 1, refe%nnode_face
                    ele%x_face(j,1:edim,i) = phys_n(j, 1:edim)
                end do
            end do

        end if
        return
    end subroutine map_ref_to_phys


    subroutine lagrange2Dquad(xr, yr, zr, nxr, nyr, xi, yi, zi, ni)
    implicit none
    integer :: nxr, nyr, ni
    real(dp) :: xr(nxr), yr(nyr), zr(nxr,nyr,2)
    real(dp) :: xi(ni), yi(ni), zi(ni,2)
    integer :: i, j, k
    real(dp) :: Lx, Ly

    do k = 1, ni
        zi(k,:) = 0.0d0
        do i = 1, nxr
        do j = 1, nyr
            call lagrange_basis_function_1D(nxr, xr, i, xi(k), Lx)
            call lagrange_basis_function_1D(nyr, yr, j, yi(k), Ly)
            zi(k,:) = zi(k,:) + zr(i,j,:) * Lx * Ly
        end do
        end do
    end do

    end subroutine lagrange2Dquad

    subroutine lagrange_basis_function_1D(nx, xd, i, xi, yi)
    ! 1D Lagrange basis function evaluation of order nx-1
    implicit none
    integer :: nx, i, j
    real(dp) :: xd(nx), xi, yi
    yi = 0.0d0
    do j = 1, nx
        if(j.ne.i) then
            yi = yi + (xi - xd(j)) / (xd(i) - xd(j))
        end if
    end do
    end subroutine lagrange_basis_function_1D


    subroutine GeometricFactors(xloc, npts, ele, drdx, dsdx, drdy, dsdy, Jac)
    use matrices
    use element_mod
    implicit none
    integer :: npts
    real(dp), dimension(npts, ndim) :: xloc
    type(element) :: ele
    real(dp), dimension(npts, nbasis, ndim) :: Vgrad
    real(dp), dimension(npts) :: drdx, dsdx, drdy, dsdy
    real(dp), dimension(npts) :: dxdr, dxds, dydr, dyds
    real(dp), dimension(npts) :: Jac

    call VandermondeGrad(Vgrad, xloc, npts, ref_ele(ele%typ)%p%order)
    if(.not. inodal) then
        ! Same for nodal and modal????
        dxdr = matmul(Vgrad(:,:,1),matmul(VanderInv, ele%x_cub(:,1)))
        if(ndim.ge.2) then
            dxds = matmul(Vgrad(:,:,2),matmul(VanderInv, ele%x_cub(:,1)))
            dydr = matmul(Vgrad(:,:,1),matmul(VanderInv, ele%x_cub(:,2)))
            dyds = matmul(Vgrad(:,:,2),matmul(VanderInv, ele%x_cub(:,2)))
        end if
    else
        dxdr = matmul(Vgrad(:,:,1),matmul(VanderInv, ele%x_cub(:,1)))
        if(ndim.ge.2) then
            dxds = matmul(Vgrad(:,:,2),matmul(VanderInv, ele%x_cub(:,1)))
            dydr = matmul(Vgrad(:,:,1),matmul(VanderInv, ele%x_cub(:,2)))
            dyds = matmul(Vgrad(:,:,2),matmul(VanderInv, ele%x_cub(:,2)))
        end if
    end if

    if(ndim.eq.1) then
        Jac = dxdr
        drdx = 1.0d0/Jac
    else if(ndim.eq.2) then
        Jac = dxdr*dyds - dxds*dydr
        drdx = dyds/Jac
        dsdx = -dydr/Jac
        drdy = -dxds/Jac
        dsdy = dxdr/Jac
    end if

    end subroutine


end module grid
