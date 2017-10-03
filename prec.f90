!     real precision definition
module prec
    integer, parameter :: p1 = selected_real_kind(15,307)
    integer, parameter :: p2 = selected_real_kind(15,307)
!   integer, parameter :: p1 = selected_real_kind(15,307)
!   integer, parameter :: sp = selected_real_kind(6,37)
    integer, parameter :: dp = selected_real_kind(15,307)
!   integer, parameter :: dp = selected_real_kind(33,4931)
end module prec
