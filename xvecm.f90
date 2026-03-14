program xvecm_johansen
! Simulate a VECM, estimate its parameters by the Johansen method for a
! specified cointegration rank, and print estimated minus true values of
! alpha, beta, pi, gamma_i, and sigma_u.
use kind_mod, only: dp
use vecm_johansen_mod, only: simulate_vecm, johansen_fit, print_fit_vs_truth
use random_mod, only: set_seed
implicit none

integer, parameter :: n = 4 ! number of time series (variables) in the VECM
integer, parameter :: r = 2 ! # of cointegrating vectors
integer, parameter :: p = 2 ! lag order of the VAR in levels; VECM has p-1 lagged differences
integer, parameter :: t_keep = 1500 ! number of kept observations
integer, parameter :: burn = 300    ! number of burn-in observations discarded

integer :: norm_rows(r)               ! rows used to normalize beta_hat to an identity matrix

real(kind=dp) :: alpha_true(n,r)      ! loading matrix in the error-correction term
real(kind=dp) :: beta_true(n,r)       ! cointegration matrix; columns are cointegrating vectors
real(kind=dp) :: gamma_true(n,n,p-1)  ! short-run coefficient matrices for lagged differences
real(kind=dp) :: sigma_u_true(n,n)    ! covariance matrix of the innovations
real(kind=dp) :: pi_true(n,n)         ! long-run matrix pi = alpha * transpose(beta)

real(kind=dp), allocatable :: y(:,:)
real(kind=dp), allocatable :: alpha_hat(:,:), beta_hat(:,:)
real(kind=dp), allocatable :: gamma_hat(:,:,:), pi_hat(:,:), sigma_u_hat(:,:)
real(kind=dp), allocatable :: lambda(:), trace_stat(:), maxeig_stat(:)

norm_rows = [1, 2]

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

call set_seed(12345)

call simulate_vecm(t_keep, alpha_true, beta_true, gamma_true, sigma_u_true, burn, y)

call johansen_fit(y, p, r, norm_rows, alpha_hat, beta_hat, gamma_hat, pi_hat, &
 sigma_u_hat, lambda, trace_stat, maxeig_stat)

call print_fit_vs_truth(alpha_hat, beta_hat, gamma_hat, pi_hat, sigma_u_hat, &
 lambda, trace_stat, maxeig_stat, alpha_true, beta_true, gamma_true, &
 pi_true, sigma_u_true)

end program xvecm_johansen
