
<!-- README.md is generated from README.Rmd. Please edit that file -->

# PoiMA

Jiahui Zou

<!-- badges: start -->
<!-- badges: end -->

Model averaging is a popular methodology to promote prediction and deal
with the uncertainty of covariate choice when we use a statistical
model. This package, PoiMA, is designed to implement model averaging for
Poisson regressions, in which the criterion to calculate model-averaging
weight is unbiased.

## Model

Suppose we observe a dataset $\{(y_i,\mathbf{x}_i),i=1,\ldots,n\}$,
where $y_i$ is the realization of a count variable $Y_i$ measuring the
number of occurrences of an even in a given interval for the $i$th
individual, $\mathbf{x}_i$ is a $p-$dimensional vector, and $n$ denotes
the number of observations. A Poisson distribution is commonly used,
which associates the probability of $Y_i$ events with $\mathbf{x}_i$ as
follows: $$
\begin{align}
P(Y_i=y_i|\mu_i)=\frac{\exp(-\mu_i)\mu_i^{y_i}}{y_i!},
\end{align}
$$ where $\mu_i$ is the Poisson incidence rate that depends on
$\mathbf{x}_i$ as $\mu_i=\exp(\mathbf{x}_i^\top\boldsymbol{\beta})$ and
$\boldsymbol{\beta}$ is the associated $p\times 1$ coefficient vector.

In practice, not all covariates are \`\`useful’’ in predicting $y_i$,
and researchers are typically uncertain regarding which should be
included in the model . Hence, a number of candidate models with
different specifications of $\mathbf{x}_i$ are considered. Let $S$ be
the number of candidate models and let $\Pi_s$ be a $p\times p_s$
selection matrix that consists of 0’s and 1’s and selects $p_s$
covariates ($p_s\le p$) for the $s^{th}$ model. Denote the
$p_s$-dimensional covariate vector in the $s$th model by
$\mathbf{x}^\top_{(s),i}=\mathbf{x}_i^\top\Pi_s$. We model the Poisson
incidence rate of this candidate model as $$
\begin{align}
\mu_{(s),i}=\exp(\mathbf{x}_{(s),i}^\top\boldsymbol{\beta}_{(s)}).
\end{align}
$$ After solving each candidate model, we need to unify the dimension of
the coefficient estimates of each candidate model by
$\hat{\boldsymbol{\beta}}_s=\Pi_s\hat{\boldsymbol{\beta}}_{(s)}$. Then,
we can compute the model averaging estimates of coefficient
$\boldsymbol{\beta}$ by $$
\begin{align}
\hat{\boldsymbol{\beta}}(\mathbf{w})=\sum_{s=1}^S {w}_s\hat{\boldsymbol{\beta}}_s,
\end{align}
$$ where $\mathbf{w}=(w_1,w_2,...,w_S)^\top$ is the weight vector
belonging to the set
$\mathcal{W}=\{\mathbf{w}\in[0,1]^S:\sum_{s=1}^S w_s=1\}$.

The model-averaging weight $\mathbf{w}$ is calculated by minimizing the
following criterion, $$
\begin{align}
C(\mathbf{w})=\log f(\mathbf{y})+\sum_{i=1}^n\left[\hat{\boldsymbol{\mu}}_i(\mathbf{w},\mathbf{y})+\log(\mathbf{y_i!})-y_i\log\{\hat{\boldsymbol{\mu}}_i(\mathbf{w},\mathbf{y}^{(y_i-1)})\}
\right],
\end{align}
$$ where
$\mathbf{y}^{(y_i-1)}=(y_1,\ldots, y_{i-1},y_i-1,y_{i+1},\ldots,y_n)\top$,
which replace the $i$th element of $\mathbf{y}$ by $y_i-1$,
$\hat{\boldsymbol{\mu}}(\mathbf{w},\mathbf{y}^{(y_i-1)})=\left(\hat{\mu}_1(\mathbf{w},\mathbf{y}^{(y_i-1)}),\ldots,\hat{\mu}_n(\mathbf{w},\mathbf{y}^{(y_i-1)}\right)\top$.
Specially, when $y_i=0$,
$y_i\log\{\hat{\mu}_i(\mathbf{w},\mathbf{y}^{(y_i-1)})\}=0$.

## Installation

You can install the development version of PoiMA from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("zoujiahuibin/PoissonMA")
```

## Example

There are basic examples which show you how to solve a common problem:

### Example 1

``` r
options(warn =-1)
library(PoiMA)
library(mvtnorm)
library(magrittr)
beta0=rep(0.1,5)
X=rmvnorm(100,mean=rep(0,5))
lambda=exp(X%*%beta0)
y=rpois(100,lambda)
ModelSetup(y,X,'nested',intercept=FALSE)
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,]    5    0    0    0    0
#> [2,]    5    1    0    0    0
#> [3,]    5    1    4    0    0
#> [4,]    5    1    4    2    0
#> [5,]    5    1    4    2    3
```

``` r
options(warn =-1)
library(PoiMA)
library(mvtnorm)
library(magrittr)
beta0=rep(0.1,5)
X=cbind(1,rmvnorm(100,mean=rep(0,4)))
lambda=exp(X%*%beta0)
y=rpois(100,lambda)
ModelSetup(y,X,'nested',intercept=TRUE)
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,]    1    5    0    0    0
#> [2,]    1    5    3    0    0
#> [3,]    1    5    3    2    0
#> [4,]    1    5    3    2    4
```

### Example 2

``` r
options(warn =-1)
library(PoiMA)
library(mvtnorm)
library(magrittr)
library(Rsolnp)
beta0=rep(0.1,5)
X=rmvnorm(100,mean=rep(0,5))
lambda=exp(X%*%beta0)
y=rpois(100,lambda)
PoissonMA(y,X, modeltype='nested',intercept=FALSE)
#> 
#> Iter: 1 fn: 99.2959   Pars:  0.5367024765 0.0000005306 0.3471932065 0.1161035354 0.0000002510
#> Iter: 2 fn: 99.2959   Pars:  0.53670276091 0.00000034716 0.34719291135 0.11610388819 0.00000009238
#> solnp--> Completed in 2 iterations
#> $Index
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,]    2    0    0    0    0
#> [2,]    2    1    0    0    0
#> [3,]    2    1    5    0    0
#> [4,]    2    1    5    3    0
#> [5,]    2    1    5    3    4
#> 
#> $weight.MA
#> [1] 5.367028e-01 3.471595e-07 3.471929e-01 1.161039e-01 9.238057e-08
#> 
#> $weight.SAIC
#> [1] 0.28072487 0.24673015 0.25219834 0.15823410 0.06211255
#> 
#> $weight.SBIC
#> [1] 0.7588276515 0.1812921316 0.0503724923 0.0085910409 0.0009166837
#> 
#> $OptimalModel.AIC
#> [1] 1
#> 
#> $OptimalModel.BIC
#> [1] 1
#> 
#> $MSE.MA
#> [1] 1.221382
#> 
#> $MSE.AIC
#> [1] 1.24796
#> 
#> $MSE.BIC
#> [1] 1.24796
#> 
#> $MSE.SAIC
#> [1] 1.217641
#> 
#> $MSE.SBIC
#> [1] 1.238622
#> 
#> $MSE.FULL
#> [1] 1.216431
#> 
#> $KLloss.MA
#> [1] 1.316388
#> 
#> $KLloss.AIC
#> [1] 1.331037
#> 
#> $KLloss.BIC
#> [1] 1.326028
#> 
#> $KLloss.SAIC
#> [1] 1.313518
#> 
#> $KLloss.SBIC
#> [1] 1.326028
#> 
#> $KLloss.Full
#> [1] 1.306121
```

``` r
options(warn =-1)
library(PoiMA)
library(mvtnorm)
library(magrittr)
library(Rsolnp)
beta0=rep(0.1,5)
X=cbind(1,rmvnorm(100,mean=rep(0,4)))
lambda=exp(X%*%beta0)
y=rpois(100,lambda)
PoissonMA(y,X, modeltype='nested',intercept=TRUE)
#> 
#> Iter: 1 fn: 94.7682   Pars:  0.153197900 0.364905009 0.481895795 0.000001297
#> Iter: 2 fn: 94.7682   Pars:  0.1531979028 0.3649050081 0.4818965668 0.0000005223
#> solnp--> Completed in 2 iterations
#> $Index
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,]    1    3    0    0    0
#> [2,]    1    3    5    0    0
#> [3,]    1    3    5    2    0
#> [4,]    1    3    5    2    4
#> 
#> $weight.MA
#> [1] 1.531979e-01 3.649050e-01 4.818966e-01 5.222952e-07
#> 
#> $weight.SAIC
#> [1] 0.04194698 0.43389211 0.38028590 0.14387502
#> 
#> $weight.SBIC
#> [1] 0.21975528 0.61789521 0.14721016 0.01513936
#> 
#> $OptimalModel.AIC
#> [1] 2
#> 
#> $OptimalModel.BIC
#> [1] 2
#> 
#> $MSE.MA
#> [1] 0.9674536
#> 
#> $MSE.AIC
#> [1] 0.978013
#> 
#> $MSE.BIC
#> [1] 0.978013
#> 
#> $MSE.SAIC
#> [1] 0.9627208
#> 
#> $MSE.SBIC
#> [1] 0.9792943
#> 
#> $MSE.FULL
#> [1] 0.9579946
#> 
#> $KLloss.MA
#> [1] 1.315458
#> 
#> $KLloss.AIC
#> [1] 1.321035
#> 
#> $KLloss.BIC
#> [1] 1.320033
#> 
#> $KLloss.SAIC
#> [1] 1.314291
#> 
#> $KLloss.SBIC
#> [1] 1.320033
#> 
#> $KLloss.Full
#> [1] 1.312074
```

# Reference

Zou, J., Wang, W., Zhang, X. and Zou, G. (2022). Optimal model averaging
for divergent-dimensional Poisson regressions. Econometric Reviews
41(7): 775–805.
