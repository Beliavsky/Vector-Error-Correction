program xvecm_rank_select
! Simulate a VECM, compute the Johansen trace and max-eigenvalue rank
! statistics, estimate critical values internally by Monte Carlo for
! ranks 0 through n-1, and report the estimated cointegration rank.
use kind_mod, only: dp
use vecm_johansen_mod, only: simulate_vecm, johansen_rank_stats, &
 determine_rank_trace, determine_rank_maxeig, monte_carlo_rank_cv
use random_mod, only: set_seed
implicit none

integer, parameter :: n = 4           ! number of time series
integer, parameter :: r = 2           ! true number of cointegrating vectors used in simulation
integer, parameter :: p = 2           ! VAR lag order in levels
integer, parameter :: t_keep = 1500   ! number of kept observations
integer, parameter :: burn = 300      ! number of burn-in observations discarded
integer, parameter :: nsim_cv = 2000  ! Monte Carlo replications per null rank

real(kind=dp), parameter :: prob_cv = 0.95_dp ! Monte Carlo critical-value quantile

integer :: r0, rank_trace, rank_maxeig

real(kind=dp) :: alpha_true(n,r)      ! loading matrix for the simulated data
real(kind=dp) :: beta_true(n,r)       ! cointegration matrix for the simulated data
real(kind=dp) :: gamma_true(n,n,p-1)  ! short-run matrices
real(kind=dp) :: sigma_u_true(n,n)    ! innovation covariance
real(kind=dp) :: pi_true(n,n)         ! long-run matrix alpha * beta'

real(kind=dp) :: alpha_cv(n,n-1)      ! template alpha columns used in Monte Carlo null simulations
real(kind=dp) :: beta_cv(n,n-1)       ! template beta columns used in Monte Carlo null simulations

real(kind=dp), allocatable :: y(:,:)         ! simulated level data
real(kind=dp), allocatable :: lambda(:)      ! Johansen eigenvalues
real(kind=dp), allocatable :: trace_stat(:)  ! trace test statistics
real(kind=dp), allocatable :: maxeig_stat(:) ! max-eigenvalue statistics

real(kind=dp) :: trace_cv(n)                 ! Monte Carlo trace critical values
real(kind=dp) :: maxeig_cv(n)                ! Monte Carlo maxeig critical values

! true beta, normalized so rows 1:r are the identity
beta_true(1,:) = [ 1.0_dp,  0.0_dp]
beta_true(2,:) = [ 0.0_dp,  1.0_dp]
beta_true(3,:) = [-1.2_dp,  0.35_dp]
beta_true(4,:) = [ 0.8_dp, -0.9_dp]

! true alpha
alpha_true(1,:) = [-0.25_dp,  0.05_dp]
alpha_true(2,:) = [ 0.12_dp, -0.18_dp]
alpha_true(3,:) = [ 0.06_dp,  0.09_dp]
alpha_true(4,:) = [-0.04_dp,  0.14_dp]

! one short-run matrix because p = 2
gamma_true(:,:,1) = 0.0_dp
gamma_true(1,:,1) = [0.20_dp, 0.02_dp, 0.00_dp, 0.00_dp]
gamma_true(2,:,1) = [0.01_dp, 0.15_dp, 0.01_dp, 0.00_dp]
gamma_true(3,:,1) = [0.00_dp, 0.02_dp, 0.10_dp, 0.01_dp]
gamma_true(4,:,1) = [0.00_dp, 0.00_dp, 0.03_dp, 0.12_dp]

! innovation covariance
sigma_u_true(1,:) = [1.00_dp, 0.30_dp, 0.20_dp, 0.10_dp]
sigma_u_true(2,:) = [0.30_dp, 1.20_dp, 0.25_dp, 0.15_dp]
sigma_u_true(3,:) = [0.20_dp, 0.25_dp, 0.90_dp, 0.35_dp]
sigma_u_true(4,:) = [0.10_dp, 0.15_dp, 0.35_dp, 1.10_dp]

pi_true = matmul(alpha_true, transpose(beta_true))

! template columns used to simulate null ranks 0, 1, 2, 3
! the first r columns match the true simulation model
! the extra column is only used for the Monte Carlo null rank r0 = 3
alpha_cv = 0.0_dp
beta_cv = 0.0_dp

alpha_cv(:,1:r) = alpha_true
beta_cv(:,1:r) = beta_true

alpha_cv(:,3) = [ 0.08_dp, -0.06_dp, 0.10_dp, 0.04_dp]
beta_cv(:,3) = [ 0.00_dp,  0.00_dp, 1.00_dp, -0.40_dp]

call set_seed(12345)

call simulate_vecm(t_keep, alpha_true, beta_true, gamma_true, sigma_u_true, burn, y)
print "('#obs #col:',*(1x,i0))", shape(y)

! compute the Johansen eigenvalues and observed statistics
call johansen_rank_stats(y, p, lambda, trace_stat, maxeig_stat)

! estimate critical values internally by Monte Carlo
call monte_carlo_rank_cv(t_keep, p, burn, alpha_cv, beta_cv, gamma_true, &
 sigma_u_true, nsim_cv, prob_cv, trace_cv, maxeig_cv)

print *
write(*,'(a,i0)') 'number of Monte Carlo replications per null rank = ', nsim_cv
write(*,'(a,f8.4)') 'critical-value quantile = ', prob_cv

print *
print '(a)', 'johansen eigenvalues'
write(*,'(100(f12.6,1x))') lambda

print *
print '(a)', 'Monte Carlo trace critical values for h0: rank <= 0, 1, ..., n-1'
write(*,'(100(f12.6,1x))') trace_cv

print *
print '(a)', 'Monte Carlo max-eigenvalue critical values for h0: rank = 0, 1, ..., n-1'
write(*,'(100(f12.6,1x))') maxeig_cv

print *
print '(a)', 'trace test decisions'
print '(a)', 'r0        statistic      crit_value     decision'
do r0 = 0, n - 1
 if (trace_stat(r0 + 1) > trace_cv(r0 + 1)) then
  write(*,'(i3,3x,f12.6,3x,f12.6,3x,a)') r0, trace_stat(r0 + 1), trace_cv(r0 + 1), 'reject'
 else
  write(*,'(i3,3x,f12.6,3x,f12.6,3x,a)') r0, trace_stat(r0 + 1), trace_cv(r0 + 1), 'fail to reject'
 end if
end do

print *
print '(a)', 'max-eigenvalue test decisions'
print '(a)', 'r0        statistic      crit_value     decision'
do r0 = 0, n - 1
 if (maxeig_stat(r0 + 1) > maxeig_cv(r0 + 1)) then
  write(*,'(i3,3x,f12.6,3x,f12.6,3x,a)') r0, maxeig_stat(r0 + 1), maxeig_cv(r0 + 1), 'reject'
 else
  write(*,'(i3,3x,f12.6,3x,f12.6,3x,a)') r0, maxeig_stat(r0 + 1), maxeig_cv(r0 + 1), 'fail to reject'
 end if
end do

rank_trace = determine_rank_trace(trace_stat, trace_cv)
rank_maxeig = determine_rank_maxeig(maxeig_stat, maxeig_cv)

print *
write(*,'(a,i0)') 'estimated rank by trace test = ', rank_trace
write(*,'(a,i0)') 'estimated rank by max-eigenvalue test = ', rank_maxeig
write(*,'(a,i0)') 'true rank used in simulation = ', r

end program xvecm_rank_select
