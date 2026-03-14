program xvecm_rank_gap
! Simulate a VECM, compute the Johansen eigenvalues, and estimate the
! cointegration rank without critical values by choosing the largest drop
! between successive ordered eigenvalues.
use kind_mod, only: dp
use vecm_johansen_mod, only: simulate_vecm, johansen_rank_stats, &
   determine_rank_gap, print_rank_gap_report
use random_mod, only: randn_vec, set_seed
implicit none

integer, parameter :: n = 4         ! number of time series
integer, parameter :: r = 2         ! true number of cointegrating vectors used in simulation
integer, parameter :: p = 2         ! VAR lag order in levels
integer, parameter :: t_keep = 1500 ! number of simulated observations kept after burn-in
integer, parameter :: burn = 300    ! number of initial observations discarded as burn-in

integer :: rank_gap

real(kind=dp) :: alpha_true(n,r)      ! loading matrix
real(kind=dp) :: beta_true(n,r)       ! cointegration matrix
real(kind=dp) :: gamma_true(n,n,p-1)  ! short-run coefficient matrices
real(kind=dp) :: sigma_u_true(n,n)    ! covariance matrix of innovations
real(kind=dp) :: pi_true(n,n)         ! long-run matrix

real(kind=dp), allocatable :: y(:,:)         ! simulated data
real(kind=dp), allocatable :: lambda(:)      ! Johansen eigenvalues
real(kind=dp), allocatable :: trace_stat(:)  ! trace statistics
real(kind=dp), allocatable :: maxeig_stat(:) ! max-eigenvalue statistics

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

! innovation covariance matrix
sigma_u_true(1,:) = [1.00_dp, 0.30_dp, 0.20_dp, 0.10_dp]
sigma_u_true(2,:) = [0.30_dp, 1.20_dp, 0.25_dp, 0.15_dp]
sigma_u_true(3,:) = [0.20_dp, 0.25_dp, 0.90_dp, 0.35_dp]
sigma_u_true(4,:) = [0.10_dp, 0.15_dp, 0.35_dp, 1.10_dp]

pi_true = matmul(alpha_true, transpose(beta_true))

call set_seed(12345)

call simulate_vecm(t_keep, alpha_true, beta_true, gamma_true, sigma_u_true, burn, y)
print "('#obs #col:',*(1x,i0))", shape(y)
! compute the Johansen eigenvalues and statistics
call johansen_rank_stats(y, p, lambda, trace_stat, maxeig_stat)

! choose rank from the largest eigenvalue gap
rank_gap = determine_rank_gap(lambda)

print *
print '(a)', 'johansen eigenvalues'
write(*,'(100(f12.6,1x))') lambda

call print_rank_gap_report(lambda)

print *
write(*,'(a,i0)') 'estimated rank by largest eigenvalue gap = ', rank_gap
write(*,'(a,i0)') 'true rank used in simulation = ', r

end program xvecm_rank_gap
