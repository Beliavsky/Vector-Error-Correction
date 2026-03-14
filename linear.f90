module linear_mod
use kind_mod, only: dp
implicit none
private
public :: cholesky_lower, invert_matrix, symmetric_inv_sqrt_spd, &
residualize_on_w, ols_fit, crossprod1, crossprod2, jacobi_eigen_sym, &
sort_eigen_desc, eye, solve_linear
contains

subroutine eye(n, a)
integer, intent(in) :: n
real(kind=dp), intent(out) :: a(n,n)
integer :: i

a = 0.0_dp
do i = 1, n
 a(i,i) = 1.0_dp
end do
end subroutine eye

subroutine cholesky_lower(a, l, info)
real(kind=dp), intent(in) :: a(:,:)
real(kind=dp), intent(out) :: l(size(a,1),size(a,2))
integer, intent(out) :: info
integer :: i, j, k, n
real(kind=dp) :: s

n = size(a,1)
if (size(a,2) /= n) then
 info = 1
 return
end if

l = 0.0_dp
info = 0

do j = 1, n
 s = a(j,j)
 do k = 1, j - 1
  s = s - l(j,k) * l(j,k)
 end do
 if (s <= 0.0_dp) then
  info = j
  return
 end if
 l(j,j) = sqrt(s)
 do i = j + 1, n
  s = a(i,j)
  do k = 1, j - 1
   s = s - l(i,k) * l(j,k)
  end do
  l(i,j) = s / l(j,j)
 end do
end do
end subroutine cholesky_lower

subroutine solve_linear(a, b, x, info)
! solve a * x = b by Gaussian elimination with partial pivoting
! a is n x n
! b is n x nrhs
real(kind=dp), intent(in) :: a(:,:), b(:,:)
real(kind=dp), allocatable, intent(out) :: x(:,:)
integer, intent(out) :: info
integer :: i, j, k, n, nrhs, pivot
real(kind=dp) :: maxabs, factor, tmp
real(kind=dp), allocatable :: aa(:,:), bb(:,:), rowa(:), rowb(:)

n = size(a,1)
if (size(a,2) /= n) then
 info = 1
 allocate(x(0,0))
 return
end if
if (size(b,1) /= n) then
 info = 2
 allocate(x(0,0))
 return
end if

nrhs = size(b,2)
allocate(aa(n,n), bb(n,nrhs), rowa(n), rowb(nrhs), x(n,nrhs))
aa = a
bb = b
info = 0

do k = 1, n - 1
 pivot = k
 maxabs = abs(aa(k,k))
 do i = k + 1, n
  if (abs(aa(i,k)) > maxabs) then
   maxabs = abs(aa(i,k))
   pivot = i
  end if
 end do
 if (maxabs <= 1.0e-14_dp) then
  info = k
  return
 end if
 if (pivot /= k) then
  rowa = aa(k,:)
  aa(k,:) = aa(pivot,:)
  aa(pivot,:) = rowa
  rowb = bb(k,:)
  bb(k,:) = bb(pivot,:)
  bb(pivot,:) = rowb
 end if
 do i = k + 1, n
  factor = aa(i,k) / aa(k,k)
  aa(i,k) = 0.0_dp
  do j = k + 1, n
   aa(i,j) = aa(i,j) - factor * aa(k,j)
  end do
  do j = 1, nrhs
   bb(i,j) = bb(i,j) - factor * bb(k,j)
  end do
 end do
end do

if (abs(aa(n,n)) <= 1.0e-14_dp) then
 info = n
 return
end if

x = 0.0_dp
do j = 1, nrhs
 x(n,j) = bb(n,j) / aa(n,n)
 do i = n - 1, 1, -1
  tmp = bb(i,j)
  do k = i + 1, n
   tmp = tmp - aa(i,k) * x(k,j)
  end do
  x(i,j) = tmp / aa(i,i)
 end do
end do
end subroutine solve_linear

subroutine invert_matrix(a, ainv, info)
real(kind=dp), intent(in) :: a(:,:)
real(kind=dp), allocatable, intent(out) :: ainv(:,:)
integer, intent(out) :: info
integer :: n
real(kind=dp), allocatable :: ident(:,:), x(:,:)

n = size(a,1)
allocate(ident(n,n))
call eye(n, ident)
call solve_linear(a, ident, x, info)
if (info /= 0) then
 allocate(ainv(0,0))
 return
end if
allocate(ainv(n,n))
ainv = x
deallocate(ident, x)
end subroutine invert_matrix

subroutine crossprod1(a, c)
! c = transpose(a) * a
real(kind=dp), intent(in) :: a(:,:)
real(kind=dp), allocatable, intent(out) :: c(:,:)
integer :: ncol

ncol = size(a,2)
allocate(c(ncol,ncol))
c = matmul(transpose(a), a)
end subroutine crossprod1

subroutine crossprod2(a, b, c)
! c = transpose(a) * b
real(kind=dp), intent(in) :: a(:,:), b(:,:)
real(kind=dp), allocatable, intent(out) :: c(:,:)
integer :: ncol_a, ncol_b

ncol_a = size(a,2)
ncol_b = size(b,2)
allocate(c(ncol_a,ncol_b))
c = matmul(transpose(a), b)
end subroutine crossprod2

subroutine jacobi_eigen_sym(a, eval, evec, info, tol, max_iter)
! Jacobi eigensolver for a real symmetric matrix
real(kind=dp), intent(in) :: a(:,:)
real(kind=dp), allocatable, intent(out) :: eval(:), evec(:,:)
integer, intent(out) :: info
real(kind=dp), intent(in), optional :: tol
integer, intent(in), optional :: max_iter
integer :: i, j, p, q, n, iter, itmax
real(kind=dp) :: threshold, app, aqq, apq, tau, t, c, s, temp
real(kind=dp) :: maxoff
real(kind=dp), allocatable :: b(:,:)

n = size(a,1)
if (size(a,2) /= n) then
 info = 1
 allocate(eval(0), evec(0,0))
 return
end if

threshold = 1.0e-10_dp
if (present(tol)) threshold = tol
itmax = 100 * n * n
if (present(max_iter)) itmax = max_iter

allocate(b(n,n), eval(n), evec(n,n))
b = a
call eye(n, evec)
info = 0

do iter = 1, itmax
 maxoff = 0.0_dp
 p = 1
 q = 1
 do i = 1, n - 1
  do j = i + 1, n
   if (abs(b(i,j)) > maxoff) then
    maxoff = abs(b(i,j))
    p = i
    q = j
   end if
  end do
 end do
 if (maxoff <= threshold) exit

 app = b(p,p)
 aqq = b(q,q)
 apq = b(p,q)
 if (abs(apq) <= threshold) cycle

 tau = (aqq - app) / (2.0_dp * apq)
 if (tau >= 0.0_dp) then
  t = 1.0_dp / (tau + sqrt(1.0_dp + tau * tau))
 else
  t = -1.0_dp / (-tau + sqrt(1.0_dp + tau * tau))
 end if
 c = 1.0_dp / sqrt(1.0_dp + t * t)
 s = t * c

 do j = 1, n
  if (j /= p .and. j /= q) then
   temp = b(j,p)
   b(j,p) = c * temp - s * b(j,q)
   b(p,j) = b(j,p)
   b(j,q) = s * temp + c * b(j,q)
   b(q,j) = b(j,q)
  end if
 end do

 b(p,p) = c * c * app - 2.0_dp * s * c * apq + s * s * aqq
 b(q,q) = s * s * app + 2.0_dp * s * c * apq + c * c * aqq
 b(p,q) = 0.0_dp
 b(q,p) = 0.0_dp

 do j = 1, n
  temp = evec(j,p)
  evec(j,p) = c * temp - s * evec(j,q)
  evec(j,q) = s * temp + c * evec(j,q)
 end do
end do

if (iter > itmax) info = 2

do i = 1, n
 eval(i) = b(i,i)
end do
end subroutine jacobi_eigen_sym

subroutine sort_eigen_desc(eval, evec)
real(kind=dp), intent(inout) :: eval(:), evec(:,:)
integer :: i, j, k, n
real(kind=dp) :: tmp
real(kind=dp), allocatable :: vtmp(:)

n = size(eval)
allocate(vtmp(size(evec,1)))
do i = 1, n - 1
 k = i
 do j = i + 1, n
  if (eval(j) > eval(k)) k = j
 end do
 if (k /= i) then
  tmp = eval(i)
  eval(i) = eval(k)
  eval(k) = tmp
  vtmp = evec(:,i)
  evec(:,i) = evec(:,k)
  evec(:,k) = vtmp
 end if
end do
deallocate(vtmp)
end subroutine sort_eigen_desc

subroutine symmetric_inv_sqrt_spd(a, ainvsqrt, eval, evec, info)
! if a = v * diag(lam) * transpose(v), return
! a^(-1/2) = v * diag(1/sqrt(lam)) * transpose(v)
real(kind=dp), intent(in) :: a(:,:)
real(kind=dp), allocatable, intent(out) :: ainvsqrt(:,:), eval(:), evec(:,:)
integer, intent(out) :: info
integer :: i, n
real(kind=dp), allocatable :: d(:,:), lam(:), v(:,:)

call jacobi_eigen_sym(a, lam, v, info)
if (info /= 0) then
 allocate(ainvsqrt(0,0), eval(0), evec(0,0))
 return
end if
call sort_eigen_desc(lam, v)

n = size(a,1)
allocate(d(n,n), ainvsqrt(n,n), eval(n), evec(n,n))
d = 0.0_dp
do i = 1, n
 if (lam(i) <= 1.0e-14_dp) then
  info = 3
  return
 end if
 d(i,i) = 1.0_dp / sqrt(lam(i))
end do
ainvsqrt = matmul(v, matmul(d, transpose(v)))
eval = lam
evec = v
info = 0
end subroutine symmetric_inv_sqrt_spd

subroutine residualize_on_w(a, w, resid)
! residuals from regressing columns of a on w
real(kind=dp), intent(in) :: a(:,:), w(:,:)
real(kind=dp), allocatable, intent(out) :: resid(:,:)
integer :: k, info
real(kind=dp), allocatable :: wtw(:,:), wta(:,:), coef(:,:)

k = size(w,2)
allocate(resid(size(a,1),size(a,2)))
if (k == 0) then
 resid = a
 return
end if

call crossprod1(w, wtw)
call crossprod2(w, a, wta)
call solve_linear(wtw, wta, coef, info)
if (info /= 0) stop 'solve_linear failed in residualize_on_w'
resid = a - matmul(w, coef)
deallocate(wtw, wta, coef)
end subroutine residualize_on_w

subroutine ols_fit(x, y, b, resid, sigma)
! b = argmin || y - x b ||^2
! sigma = transpose(resid) * resid / nobs
real(kind=dp), intent(in) :: x(:,:), y(:,:)
real(kind=dp), allocatable, intent(out) :: b(:,:), resid(:,:), sigma(:,:)
integer :: info, m
real(kind=dp), allocatable :: xtx(:,:), xty(:,:)

call crossprod1(x, xtx)
call crossprod2(x, y, xty)
call solve_linear(xtx, xty, b, info)
if (info /= 0) stop 'solve_linear failed in ols_fit'
allocate(resid(size(y,1),size(y,2)))
resid = y - matmul(x, b)
m = size(y,1)
call crossprod1(resid, sigma)
sigma = sigma / real(m, kind=dp)
deallocate(xtx, xty)
end subroutine ols_fit

end module linear_mod