module vecm_johansen_mod
use kind_mod, only: dp
use linear_mod, only: crossprod1, crossprod2, invert_matrix, symmetric_inv_sqrt_spd, &
   jacobi_eigen_sym, sort_eigen_desc, residualize_on_w, eye, solve_linear, ols_fit, &
   cholesky_lower
use random_mod, only: randn_vec, set_seed
implicit none
private
public :: simulate_vecm, johansen_fit, print_fit_vs_truth, &
   determine_rank_trace, determine_rank_maxeig, johansen_rank_stats, &
   determine_rank_gap, print_rank_gap_report

contains

subroutine simulate_vecm(t_keep, alpha, beta, gamma, sigma_u, burn, y)
! simulate
!
! d y_t = alpha * transpose(beta) * y_{t-1}
!       + sum gamma_i * d y_{t-i}
!       + u_t
!
! with u_t ~ n(0, sigma_u)
integer, intent(in) :: t_keep, burn
real(kind=dp), intent(in) :: alpha(:,:), beta(:,:), gamma(:,:,:), sigma_u(:,:)
real(kind=dp), allocatable, intent(out) :: y(:,:)
integer :: n, r, p, t_total, t, i, info
real(kind=dp), allocatable :: pi(:,:), l(:,:), z(:), u(:), dy(:), yall(:,:)

n = size(alpha,1)
r = size(alpha,2)
if (size(beta,1) /= n .or. size(beta,2) /= r) stop 'bad beta dimensions'
if (size(sigma_u,1) /= n .or. size(sigma_u,2) /= n) stop 'bad sigma_u dimensions'
p = size(gamma,3) + 1

t_total = t_keep + burn
allocate(pi(n,n), l(n,n), z(n), u(n), dy(n), yall(t_total,n))
pi = matmul(alpha, transpose(beta))
yall = 0.0_dp

call cholesky_lower(sigma_u, l, info)
if (info /= 0) stop 'sigma_u is not positive definite'

do t = p + 1, t_total
 dy = matmul(pi, yall(t - 1,:))
 do i = 1, p - 1
  dy = dy + matmul(gamma(:,:,i), yall(t - i,:) - yall(t - i - 1,:))
 end do
 call randn_vec(z)
 u = matmul(l, z)
 dy = dy + u
 yall(t,:) = yall(t - 1,:) + dy
end do

allocate(y(t_keep,n))
y = yall(burn + 1:t_total,:)
deallocate(pi, l, z, u, dy, yall)
end subroutine simulate_vecm

subroutine make_vecm_matrices(y, p, dy, ylag1, w)
! build the usual VECM regression matrices:
!
! dy    = d y_t
! ylag1 = y_{t-1}
! w     = [ d y_{t-1}, d y_{t-2}, ..., d y_{t-p+1} ]
real(kind=dp), intent(in) :: y(:,:)
integer, intent(in) :: p
real(kind=dp), allocatable, intent(out) :: dy(:,:), ylag1(:,:), w(:,:)
integer :: tobs, n, m, row, i, j0, j1

n = size(y,2)
tobs = size(y,1)
m = tobs - p
if (m <= 0) stop 'not enough observations in make_vecm_matrices'

allocate(dy(m,n), ylag1(m,n), w(m, n * max(0, p - 1)))

do row = 1, m
 dy(row,:) = y(row + p,:) - y(row + p - 1,:)
 ylag1(row,:) = y(row + p - 1,:)
end do

do i = 1, p - 1
 j0 = (i - 1) * n + 1
 j1 = i * n
 do row = 1, m
  w(row,j0:j1) = y(row + p - i,:) - y(row + p - i - 1,:)
 end do
end do
end subroutine make_vecm_matrices

subroutine copy_selected_rows(a, rows, b)
real(kind=dp), intent(in) :: a(:,:)
integer, intent(in) :: rows(:)
real(kind=dp), intent(out) :: b(size(rows), size(a,2))
integer :: i

do i = 1, size(rows)
 b(i,:) = a(rows(i),:)
end do
end subroutine copy_selected_rows

subroutine johansen_fit(y, p, r, norm_rows, &
 alpha_hat, beta_hat, gamma_hat, pi_hat, sigma_u_hat, &
 lambda, trace_stat, maxeig_stat)
!
! base-Fortran Johansen estimator for a VECM
!
! steps:
! 1. partial out lagged differences w from d y_t and y_{t-1}
! 2. form s00, s01, s10, s11
! 3. solve the symmetric eigenproblem
!
!    k = s11^(-1/2) * s10 * s00^(-1) * s01 * s11^(-1/2)
!
! 4. recover beta from the eigenvectors of k
! 5. normalize beta so selected rows equal the identity
! 6. regress d y_t on beta' y_{t-1} and lagged differences to get
!    alpha and gamma_i
!
real(kind=dp), intent(in) :: y(:,:)
integer, intent(in) :: p, r
integer, intent(in) :: norm_rows(:)
real(kind=dp), allocatable, intent(out) :: alpha_hat(:,:), beta_hat(:,:)
real(kind=dp), allocatable, intent(out) :: gamma_hat(:,:,:), pi_hat(:,:)
real(kind=dp), allocatable, intent(out) :: sigma_u_hat(:,:), lambda(:)
real(kind=dp), allocatable, intent(out) :: trace_stat(:), maxeig_stat(:)
integer :: n, m, i, j, info, kx, j0, j1
real(kind=dp), allocatable :: dy(:,:), ylag1(:,:), w(:,:), r0(:,:), r1(:,:)
real(kind=dp), allocatable :: s00(:,:), s01(:,:), s10(:,:), s11(:,:)
real(kind=dp), allocatable :: s00_inv(:,:), s11_invsqrt(:,:)
real(kind=dp), allocatable :: eval11(:), evec11(:,:), kmat(:,:)
real(kind=dp), allocatable :: evalk(:), eveck(:,:), beta_raw(:,:)
real(kind=dp), allocatable :: top(:,:), ident_r(:,:), xnorm(:,:)
real(kind=dp), allocatable :: ect(:,:), x(:,:), b(:,:), resid(:,:), sigma(:,:)
real(kind=dp), allocatable :: logterm(:)

n = size(y,2)
call make_vecm_matrices(y, p, dy, ylag1, w)
m = size(dy,1)

call residualize_on_w(dy, w, r0)
call residualize_on_w(ylag1, w, r1)

call crossprod1(r0, s00)
call crossprod2(r0, r1, s01)
allocate(s10(size(s01,2), size(s01,1)))
s10 = transpose(s01)
call crossprod1(r1, s11)

s00 = s00 / real(m, kind=dp)
s01 = s01 / real(m, kind=dp)
s10 = s10 / real(m, kind=dp)
s11 = s11 / real(m, kind=dp)

call invert_matrix(s00, s00_inv, info)
if (info /= 0) stop 'invert_matrix failed for s00'

call symmetric_inv_sqrt_spd(s11, s11_invsqrt, eval11, evec11, info)
if (info /= 0) stop 'symmetric_inv_sqrt_spd failed for s11'

allocate(kmat(n,n))
kmat = matmul(s11_invsqrt, matmul(s10, matmul(s00_inv, matmul(s01, s11_invsqrt))))
kmat = 0.5_dp * (kmat + transpose(kmat))

call jacobi_eigen_sym(kmat, evalk, eveck, info)
if (info /= 0) stop 'jacobi_eigen_sym failed for johansen matrix'
call sort_eigen_desc(evalk, eveck)

allocate(lambda(n), trace_stat(n), maxeig_stat(n))
do i = 1, n
 lambda(i) = min(1.0_dp, max(0.0_dp, evalk(i)))
end do

allocate(beta_raw(n,r))
beta_raw = matmul(s11_invsqrt, eveck(:,1:r))

allocate(top(r,r), ident_r(r,r))
call copy_selected_rows(beta_raw, norm_rows, top)
call eye(r, ident_r)
call solve_linear(top, ident_r, xnorm, info)
if (info /= 0) stop 'normalization of beta failed'

allocate(beta_hat(n,r))
beta_hat = matmul(beta_raw, xnorm)

allocate(ect(m,r))
ect = matmul(ylag1, beta_hat)

kx = r + n * max(0, p - 1)
allocate(x(m,kx))
x(:,1:r) = ect
if (p > 1) then
 x(:,r + 1:kx) = w
end if

call ols_fit(x, dy, b, resid, sigma)

allocate(alpha_hat(n,r))
alpha_hat = transpose(b(1:r,:))

allocate(gamma_hat(n,n,max(0, p - 1)))
gamma_hat = 0.0_dp
do i = 1, p - 1
 j0 = r + (i - 1) * n + 1
 j1 = r + i * n
 gamma_hat(:,:,i) = transpose(b(j0:j1,:))
end do

allocate(pi_hat(n,n), sigma_u_hat(n,n))
pi_hat = matmul(alpha_hat, transpose(beta_hat))
sigma_u_hat = sigma

do i = 1, n
 allocate(logterm(n - i + 1))
 do j = i, n
  logterm(j - i + 1) = log(max(1.0e-12_dp, 1.0_dp - lambda(j)))
 end do
 trace_stat(i) = -real(m, kind=dp) * sum(logterm)
 deallocate(logterm)
 maxeig_stat(i) = -real(m, kind=dp) * log(max(1.0e-12_dp, 1.0_dp - lambda(i)))
end do
end subroutine johansen_fit

function itoa(i) result(s)
integer, intent(in) :: i
character(len=32) :: s
write(s,'(i0)') i
end function itoa

subroutine print_matrix_diff(name, est, truth)
character(len=*), intent(in) :: name
real(kind=dp), intent(in) :: est(:,:), truth(:,:)
integer :: i

print *
print '(a)', trim(name)//' : estimated - true'
do i = 1, size(est,1)
 write(*,'(100(f12.6,1x))') est(i,:) - truth(i,:)
end do
end subroutine print_matrix_diff

subroutine print_fit_vs_truth(alpha_hat, beta_hat, gamma_hat, pi_hat, sigma_u_hat, &
 lambda, trace_stat, maxeig_stat, alpha_true, beta_true, gamma_true, &
 pi_true, sigma_u_true)
real(kind=dp), intent(in) :: alpha_hat(:,:), beta_hat(:,:), gamma_hat(:,:,:)
real(kind=dp), intent(in) :: pi_hat(:,:), sigma_u_hat(:,:), lambda(:)
real(kind=dp), intent(in) :: trace_stat(:), maxeig_stat(:)
real(kind=dp), intent(in) :: alpha_true(:,:), beta_true(:,:), gamma_true(:,:,:)
real(kind=dp), intent(in) :: pi_true(:,:), sigma_u_true(:,:)
integer :: i

print *
print '(a)', 'johansen eigenvalues'
write(*,'(100(f12.6,1x))') lambda

print *
print '(a)', 'johansen trace statistics for h0: rank <= 0, 1, ..., n-1'
write(*,'(100(f12.6,1x))') trace_stat

print *
print '(a)', 'johansen max-eigenvalue statistics'
write(*,'(100(f12.6,1x))') maxeig_stat

call print_matrix_diff('alpha', alpha_hat, alpha_true)
call print_matrix_diff('beta', beta_hat, beta_true)
call print_matrix_diff('pi', pi_hat, pi_true)

do i = 1, size(gamma_hat,3)
 call print_matrix_diff('gamma_'//trim(adjustl(itoa(i))), &
  gamma_hat(:,:,i), gamma_true(:,:,i))
end do

call print_matrix_diff('sigma_u', sigma_u_hat, sigma_u_true)
end subroutine print_fit_vs_truth

subroutine johansen_rank_stats(y, p, lambda, trace_stat, maxeig_stat)
! compute the Johansen eigenvalues and rank-test statistics
! these do not depend on the chosen cointegration rank r

real(kind=dp), intent(in) :: y(:,:)
integer, intent(in) :: p
real(kind=dp), allocatable, intent(out) :: lambda(:), trace_stat(:), maxeig_stat(:)

integer :: n, m, i, info
real(kind=dp), allocatable :: dy(:,:), ylag1(:,:), w(:,:), r0(:,:), r1(:,:)
real(kind=dp), allocatable :: s00(:,:), s01(:,:), s10(:,:), s11(:,:)
real(kind=dp), allocatable :: s00_inv(:,:), s11_invsqrt(:,:)
real(kind=dp), allocatable :: eval11(:), evec11(:,:), kmat(:,:)
real(kind=dp), allocatable :: evalk(:), eveck(:,:)

n = size(y,2)

call make_vecm_matrices(y, p, dy, ylag1, w)
m = size(dy,1)

call residualize_on_w(dy, w, r0)
call residualize_on_w(ylag1, w, r1)

call crossprod1(r0, s00)
call crossprod2(r0, r1, s01)
allocate(s10(size(s01,2), size(s01,1)))
s10 = transpose(s01)
call crossprod1(r1, s11)

s00 = s00 / real(m, kind=dp)
s01 = s01 / real(m, kind=dp)
s10 = s10 / real(m, kind=dp)
s11 = s11 / real(m, kind=dp)

call invert_matrix(s00, s00_inv, info)
if (info /= 0) stop 'invert_matrix failed for s00 in johansen_rank_stats'

call symmetric_inv_sqrt_spd(s11, s11_invsqrt, eval11, evec11, info)
if (info /= 0) stop 'symmetric_inv_sqrt_spd failed for s11 in johansen_rank_stats'

allocate(kmat(n,n))
kmat = matmul(s11_invsqrt, matmul(s10, matmul(s00_inv, matmul(s01, s11_invsqrt))))
kmat = 0.5_dp * (kmat + transpose(kmat))

call jacobi_eigen_sym(kmat, evalk, eveck, info)
if (info /= 0) stop 'jacobi_eigen_sym failed in johansen_rank_stats'

call sort_eigen_desc(evalk, eveck)

allocate(lambda(n), trace_stat(n), maxeig_stat(n))

do i = 1, n
 lambda(i) = min(1.0_dp, max(0.0_dp, evalk(i)))
end do

do i = 1, n
 trace_stat(i) = -real(m, kind=dp) * sum(log(max(1.0e-12_dp, 1.0_dp - lambda(i:n))))
 maxeig_stat(i) = -real(m, kind=dp) * log(max(1.0e-12_dp, 1.0_dp - lambda(i)))
end do

end subroutine johansen_rank_stats

integer function determine_rank_trace(trace_stat, trace_cv) result(rank_est)
! sequential trace test:
! h0: rank <= 0, 1, ..., n-1
! estimated rank is the first r0 for which h0 is not rejected
! if all nulls are rejected, this returns n

real(kind=dp), intent(in) :: trace_stat(:), trace_cv(:)

integer :: n, r0

n = size(trace_stat)
if (size(trace_cv) /= n) stop 'trace_cv has wrong length'

rank_est = n
do r0 = 0, n - 1
 if (trace_stat(r0 + 1) <= trace_cv(r0 + 1)) then
  rank_est = r0
  return
 end if
end do

end function determine_rank_trace

integer function determine_rank_maxeig(maxeig_stat, maxeig_cv) result(rank_est)
! sequential max-eigenvalue test:
! h0: rank = 0, 1, ..., n-1 against h1: rank = r0 + 1
! estimated rank is the first r0 for which h0 is not rejected
! if all nulls are rejected, this returns n

real(kind=dp), intent(in) :: maxeig_stat(:), maxeig_cv(:)

integer :: n, r0

n = size(maxeig_stat)
if (size(maxeig_cv) /= n) stop 'maxeig_cv has wrong length'

rank_est = n
do r0 = 0, n - 1
 if (maxeig_stat(r0 + 1) <= maxeig_cv(r0 + 1)) then
  rank_est = r0
  return
 end if
end do

end function determine_rank_maxeig

integer function determine_rank_gap(lambda, tol_small) result(rank_est)
! determine the cointegration rank from the ordered Johansen eigenvalues
! by finding the largest adjacent drop:
!
!   rank = argmax_i lambda(i) / lambda(i+1)
!
! this is a heuristic, not a formal hypothesis test

real(kind=dp), intent(in) :: lambda(:)
real(kind=dp), intent(in), optional :: tol_small

integer :: i, n, imax
real(kind=dp) :: eps, best_ratio, ratio

n = size(lambda)
eps = 1.0e-8_dp
if (present(tol_small)) eps = tol_small

if (n <= 1) then
 rank_est = 0
 return
end if

if (lambda(1) <= eps) then
 rank_est = 0
 return
end if

best_ratio = -1.0_dp
imax = 0

do i = 1, n - 1
 ratio = lambda(i) / max(lambda(i + 1), eps)
 if (ratio > best_ratio) then
  best_ratio = ratio
  imax = i
 end if
end do

rank_est = imax

end function determine_rank_gap

subroutine print_rank_gap_report(lambda, tol_small)
! print the adjacent eigenvalue ratios used by determine_rank_gap

real(kind=dp), intent(in) :: lambda(:)
real(kind=dp), intent(in), optional :: tol_small

integer :: i, n
real(kind=dp) :: eps, ratio

n = size(lambda)
eps = 1.0e-8_dp
if (present(tol_small)) eps = tol_small

print *
print '(a)', 'eigenvalue-gap report'
print '(a)', 'i        lambda(i)      lambda(i+1)      ratio'

do i = 1, n - 1
 ratio = lambda(i) / max(lambda(i + 1), eps)
 write(*,'(i3,3x,f12.6,3x,f12.6,3x,f12.6)') i, lambda(i), lambda(i + 1), ratio
end do

end subroutine print_rank_gap_report

end module vecm_johansen_mod