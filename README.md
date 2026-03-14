# Vector-Error-Correction
Simulate and fit from Vector Error Correction (VECM) models for cointegrated time series using the [Johansen method](https://academic.oup.com/book/27916).

```fortran
program xvecm_rank_select
! Simulate a VECM, compute the Johansen trace and max-eigenvalue rank
! statistics, estimate critical values internally by Monte Carlo for
! ranks 0 through n-1, and report the estimated cointegration rank.
```
sample output is
```
#obs #col: 1500 4

number of Monte Carlo replications per null rank = 2000
critical-value quantile =   0.9500

johansen eigenvalues
    0.287478     0.129449     0.001799     0.000222

Monte Carlo trace critical values for h0: rank <= 0, 1, ..., n-1
   41.033639    24.300203    12.765328   158.320284

Monte Carlo max-eigenvalue critical values for h0: rank = 0, 1, ..., n-1
   24.603392    17.420267    11.609917   158.320284

trace test decisions
r0        statistic      crit_value     decision
  0     718.434880      41.033639   reject
  1     210.695064      24.300203   reject
  2       3.029339      12.765328   fail to reject
  3       0.331884     158.320284   fail to reject

max-eigenvalue test decisions
r0        statistic      crit_value     decision
  0     507.739816      24.603392   reject
  1     207.665725      17.420267   reject
  2       2.697455      11.609917   fail to reject
  3       0.331884     158.320284   fail to reject

estimated rank by trace test = 2
estimated rank by max-eigenvalue test = 2
true rank used in simulation = 2
```
