module stability

use glob
use prec
use matrices
use poly

implicit none

contains
subroutine stable_dt
implicit none

real(dp) :: dt, dx, wavenumber
complex(dp) :: SD(Np,Np)
complex(dp) :: SD_powers(Np,Np)
complex(dp) :: D(Np,Np) ! Time integration and space discretization
complex(dp) :: eig(Np), maxeig
complex(dp), allocatable :: eigarr(:,:)
integer :: i
real(dp) :: fac, identity(Np,Np), gain(Np)

integer :: it, ik, nt, nk, nrk

nk = 100
nt = 60000

allocate(eigarr(nk,Np))

identity = 0.0d0
do i = 1, Np
    identity(i,i) = 1.0d0
end do

dx = 1.0d0
do nrk = 2, 4
do it = 1, nt
    dt = 10.0d0 / (1.1**(.01*it))
    do ik = 1, nk
        wavenumber = 2.0d0*PI/(nk-1) * (ik-1)
        call build_semidiscrete(wavenumber,SD)

        D = identity
        fac = 1.0d0
        SD_powers = identity
        do i = 1, nrk
            if(i.eq.5) cycle
            fac = fac * i
            SD_powers = matmul(SD_powers,SD)
            D = D + dt**i / fac * SD_powers
        end do
        if(nrk.eq.5) then
            i = 5
            fac = 200.0d0
            SD_powers = matmul(SD,SD_powers)
            D = D + dt**i / fac * SD_powers
        end if

        call complex_eig(D,eig,Np)
        gain = real(eig)**2 + aimag(eig)**2
        if(maxval(gain) .ge. 1.0d0+1e-12) exit
        call complex_eig(SD,eig,Np)
        eigarr(ik,:) = eig
    end do
    if(it.eq.nt) then
        print*, 'Didnt find stable dt', dt
        stop
    end if
    if(ik.gt.nk) then
        !print*, dt
        do ik = 1, nk
            do i = 1, Np
                print*, dt*real(eigarr(ik,i)), dt*aimag(eigarr(ik,i))
            end do
        end do
        exit
    end if
end do ! dt
end do ! nrk

return
end subroutine stable_dt

subroutine check_stability(wavenumber, eig)
implicit none

real(dp) :: wavenumber
complex(dp) :: SD(Np,Np) ! Semi-discrete matrix [d/dt] [u] = [S] [u]
complex(dp) :: eig(Np)

call build_semidiscrete(wavenumber,SD)
if(istab.eq.1) call print_cmatrix(SD)

call complex_eig(SD,eig,Np)
call print_cvector(eig)

end subroutine check_stability

subroutine build_semidiscrete(wavenumber, SD)
implicit none

real(dp) :: wavenumber
real(dp) :: C_m(Np,Np)
real(dp) :: C_0(Np,Np)
complex(dp) :: SD(Np,Np) ! Semi-discrete matrix [d/dt] [u] = [S] [u]
complex(dp) :: du(2,Np)
complex(dp) :: im
real(dp) :: h = 1.0d0 ! Cell length
im = (0.0d0, 1.0d0)

! Alternate form. S = Cm*exp(-I*w) + C0
!du(1,:) = - VanderSurf(:,1)
!du(2,:) = 0.0d0
!C_0 = (refb-refa)/h * (-Differentiation + matmul(Lift, du))
!
!du(1,:) = VanderSurf(:,2)
!du(2,:) = 0.0d0
!C_m = (refb-refa)/h * (matmul(Lift, du))
!
!SD = C_m*exp(-im*wavenumber) + C_0

! VanderSurf(:,1) = poly_eval(X_L)
! VanderSurf(:,2) = poly_eval(X_R)
du(1,:) = VanderSurf(:,2) * exp(-im*wavenumber) - VanderSurf(:,1)
du(2,:) = 0.0d0

SD = (refb-refa)/h * (-Differentiation + matmul(Lift, du))

!call print_cmatrix(SD)

return
end subroutine build_semidiscrete

subroutine print_cmatrix(A)
complex(dp) :: A(:,:)
integer :: i, j
do i=1,size(A,1)
    do j=1,size(A,2)
        if(abs(aimag(A(i,j))).gt.0.0d0) then
            write(*,'(SP,100(E10.2,E9.2,"i,"))',advance="no") real(A(i,j)), aimag(A(i,j))
        else
            !write(*,'(SP,100(E10.2,10X,","))',advance="no") real(A(i,j))
            write(*,'(SP,100(E23.15,10X,","))',advance="no") real(A(i,j))
        end if
    enddo
    write(*,*)
enddo
return
end subroutine print_cmatrix
subroutine print_cvector(v)
complex(dp) :: v(:)
integer :: i
do i=1,size(v,1)
    write(*,*) real(v(i)), aimag(v(i))
!   if(abs(aimag(v(i))).gt.0.0d0) then
!       write(*,'(SP,100(E10.3,E10.3,"i"))') real(v(i)), aimag(v(i))
!   else
!       write(*,'(SP,100(E10.3))') real(v(i))
!       !write(*,'(SP,100(E13.6))') real(v(i))
!   end if
enddo
return
end subroutine print_cvector

subroutine complex_eig(A,eig,n)
implicit none
integer :: n
complex(dp) :: A(n,n)
complex(dp) :: eig(n)
integer, parameter :: lwmax = 9999
integer :: info, lwork
real(dp) :: rwork(2*n)
complex(dp) :: VL(n,n), VR(n,n) ! Eigenvectors
complex(dp) :: Atemp(n,n)
complex(dp) :: work(lwmax)

integer :: matz = 0 ! Want eigenvectors
character(1) :: crv, clv

crv = 'N'; clv = 'N'
if(matz.eq.1) crv ='V'; if(matz.eq.2) clv ='V'

Atemp = A!cmplx(real(A),aimag(A))

lwork = -1
work = 0.0d0
call zgeev(clv, crv, n, Atemp, n, eig, VL, n, VR, n, work, lwork, rwork, info)
lwork = min( lwmax, int( real(work( 1 )) ) )
call zgeev(clv, crv, n, Atemp, n, eig, VL, n, VR, n, work, lwork, rwork, info)
if(info.ne.0) print*,'INFO.NE.0, INFO', INFO

return
end subroutine complex_eig

end module stability
