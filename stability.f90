module stability

use glob
use prec
use matrices
use poly

implicit none

contains
subroutine check_stability(wavenumber)
implicit none

real(dp) :: wavenumber
!real(dp) :: u_ele(nref)
complex(dp) :: SD(nref,nref) ! Semi-discrete matrix [d/dt] [u] = [S] [u]
complex(dp) :: eig(nref)
complex(dp) :: du(2,nref)
real(dp) :: identity(nref,nref)
integer :: i, j
complex(dp) :: im
im = (0.0d0, 1.0d0)

identity = 0.0d0
do i = 1, nref
    identity(i,i) = 1
end do
identity(1,1) = 2

SD = identity + cmplx(identity,2*identity)

! VanderSurf(:,1) = poly_eval(X_L)
! VanderSurf(:,2) = poly_eval(X_R)
du(1,:) = VanderSurf(:,2) * exp(-im*wavenumber) - VanderSurf(:,1)
du(2,:) = 0.0d0
SD = -Differentiation + matmul(Lift, du)

!call print_cmatrix(SD)
print*, 'Wavenumber:',wavenumber
call complex_eig(SD,eig,nref)
call print_cvector(eig)

end subroutine check_stability

subroutine print_cmatrix(A)
complex(dp) :: A(:,:)
integer :: i, j
do i=1,size(A,1)
    do j=1,size(A,2)
        if(abs(aimag(A(i,j))).gt.0.0d0) then
            write(*,'(100(E9.2,"+i",E8.2))',advance="no") real(A(i,j)), aimag(A(i,j))
        else
            write(*,'(100(E20.3))',advance="no") real(A(i,j))
        end if
    enddo
    write(*,*)
enddo
return
end subroutine print_cmatrix
subroutine print_cvector(v)
complex(dp) :: v(:)
integer :: i, j
do i=1,size(v,1)
    if(abs(aimag(v(i))).gt.0.0d0) then
        write(*,'(100(E9.2,"+i",E8.2))') real(v(i)), aimag(v(i))
    else
        write(*,'(100(E20.3))') real(v(i))
    end if
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
complex(dp) :: Atemp(n,n)! Eigenvectors
complex(dp) :: w(n), work(lwmax)

integer :: matz = 0 ! Want eigenvectors
character(1) :: crv, clv

crv = 'N'; clv = 'N'
if(matz.eq.1) crv ='V'; if(matz.eq.2) clv ='V'

Atemp = cmplx(real(A),aimag(A))

lwork = -1
work = 0.0d0
call zgeev(clv, crv, n, Atemp, n, eig, VL, n, VR, n, work, lwork, rwork, info)
lwork = min( lwmax, int( real(work( 1 )) ) )
call zgeev(clv, crv, n, Atemp, n, eig, VL, n, VR, n, work, lwork, rwork, info)
!print*,'info',info

return
end subroutine complex_eig

end module stability
