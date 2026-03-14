module random_mod
use kind_mod, only: dp
implicit none
private
public :: randn_vec, set_seed
contains

subroutine randn_vec(z)
real(kind=dp), intent(out) :: z(:)
integer :: i, m
real(kind=dp) :: u1, u2, r, theta, pi

pi = acos(-1.0_dp)
m = size(z)
i = 1

do while (i <= m)
 call random_number(u1)
 call random_number(u2)
 u1 = max(u1, 1.0e-12_dp)
 r = sqrt(-2.0_dp * log(u1))
 theta = 2.0_dp * pi * u2
 z(i) = r * cos(theta)
 if (i + 1 <= m) then
  z(i + 1) = r * sin(theta)
 end if
 i = i + 2
end do
end subroutine randn_vec

subroutine set_seed(seed0)
integer, intent(in) :: seed0
integer :: i, nseed
integer, allocatable :: seed(:)
call random_seed(size=nseed)
allocate(seed(nseed))
do i = 1, nseed
 seed(i) = seed0 + 37 * (i - 1)
end do
call random_seed(put=seed)
deallocate(seed)
end subroutine set_seed
end module random_mod
