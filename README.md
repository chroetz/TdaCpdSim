TdaCpdSim – Package Description
================
Christof
02 October, 2021

-   [1 Introduction](#introduction)
    -   [1.1 Installation](#installation)
    -   [1.2 Rough Outline of a
        Simulation](#rough-outline-of-a-simulation)
-   [2 Simulation](#simulation)
    -   [2.1 Samplers](#samplers)
        -   [2.1.1 Distributions](#distributions)
        -   [2.1.2 Kind of Sampling Method](#kind-of-sampling-method)
    -   [2.2 Estimators](#estimators)
        -   [2.2.1 Cusum Post-Processor](#cusum-post-processor)
        -   [2.2.2 Implemented
            Pre-Processor](#implemented-pre-processor)
    -   [2.3 Required Dists](#required-dists)
    -   [2.4 Save and Load](#save-and-load)
-   [3 Evaluate](#evaluate)
    -   [3.1 Relative Mean Absolute
        Error](#relative-mean-absolute-error)
    -   [3.2 Visualization](#visualization)
        -   [3.2.1 Plot Summary](#plot-summary)
        -   [3.2.2 Plot All Runs](#plot-all-runs)
    -   [3.3 Pretty Tables](#pretty-tables)
-   [4 Further Notes](#further-notes)
    -   [4.1 Parallel](#parallel)
    -   [4.2 Reproducibility](#reproducibility)
-   [5 Extending the Package](#extending-the-package)
    -   [5.1 New Pre-Processors](#new-pre-processors)
    -   [5.2 New Post-Processors](#new-post-processors)
    -   [5.3 New Detectors](#new-detectors)
    -   [5.4 New Image Distributions](#new-image-distributions)
    -   [5.5 Other Distributions](#other-distributions)
    -   [5.6 New Kinds of Samplers](#new-kinds-of-samplers)

# 1 Introduction

The R-package `TdaCpdSim` <https://github.com/chroetz/TdaCpdSim>
provides functions to execute and evaluate simulations for the
evaluation of different change point detection (CP/CPD) methods in time
series (TS) using topological data analysis (TDA).

## 1.1 Installation

To install the package from github, execute following commands:

``` r
install.packages("remotes")
remotes::install_github("chroetz/TdaCpdSim")
```

Call user functions of the package by `TdaCpdSim::user_function()` or by
calling `library(TdaCpdSim)` once and then `user_function()`. Call
internal functions of the package by `TdaCpdSim:::internal_function()`.
The package makes use of the `tidyverse`, a collection of packages for
data manipulation. So it is usually convenient to load the `tidyverse`.

``` r
options(tidyverse.quiet = TRUE)
library(tidyverse)
library(TdaCpdSim)
```

## 1.2 Rough Outline of a Simulation

1.  Sample *n* sets
    **X**<sub>1</sub>, …, **X**<sub>*n*</sub> ⊂ ℝ<sup>*p*</sup> of
    points in ℝ<sup>*p*</sup> (usually *p* = 2) to create a time series
    of length *n*. Typically each set has the same size
    \|**X**<sub>*i*</sub>\| = *m*, but in principle different sizes are
    allowed.
2.  Apply a filtration to **X**<sub>1</sub>, …, **X**<sub>*n*</sub> to
    obtain persistence diagrams *P*<sub>1</sub>, …, *P*<sub>*n*</sub>.
3.  Apply estimators for a change point in the time series.
4.  Repeat 1 to 3 several times.
5.  Evaluate the relative mean absolute error of the estimators.

# 2 Simulation

To run a simulation, call the function `TdaCpdSim::simulate_ts()`.

``` r
# samplers <- ...
# estimators <- ...
# required_dists <- ...
results <- simulate_ts(samplers, estimators, required_dists)
```

All three arguments are tables (R-objects of class `tibble`).

The rows of `samplers` describe the distributions of the time series
sampled during the simulation, as well as the number of repetitions and
the filtration applied to the time series of sets of points.

The rows of `required_dists` describe for which distances *d* of
persistence diagrams the values
*d*(*P*<sub>*i*</sub>, *P*<sub>*j*</sub>) are calculated for use in the
estimators. This is not done by individual estimators as it may be
computationally intensive and should not be calculated multiple times.

The rows of `estimators` describe the estimators applied to the time
series. They take the time series of sets of points, the time series of
persistence diagrams, and the distance matrices as arguments and
ultimately produce an element in \[1, *n*\] as estimated value of the
change point.

## 2.1 Samplers

Here is an example of a table of samplers.

``` r
# s_name must be unique
samplers <- tribble(
  ~s_name, ~distri1, ~distri2,
  "o-.", "circ", "dot",
  "o-n", "circ", "horseshoe") %>%
  mutate(m=20, n=300, n1=200, rate=0.05,
         kind = "birth_death",
         filt = "Rips",
         reps=3)
samplers
## # A tibble: 2 x 10
##   s_name distri1 distri2       m     n    n1  rate kind        filt   reps
##   <chr>  <chr>   <chr>     <dbl> <dbl> <dbl> <dbl> <chr>       <chr> <dbl>
## 1 o-.    circ    dot          20   300   200  0.05 birth_death Rips      3
## 2 o-n    circ    horseshoe    20   300   200  0.05 birth_death Rips      3
```

-   `s_name`: a unique identifier of the sampler.
-   `n`: the length of the TS to be produced by the sampler.
-   `kind`: describes which sampling method to use. Currently
    `"birth_death"` or `"iid"`.
-   `filt`: filtration to create persistence diagram as supported by the
    `TDA` R-package. Currently `"Rips"` or `"Alpha"`.
-   `reps`: number of repetitions of one experiment (sample TS, apply
    estimators)
-   `n1`: The TS has two parts of length *n*<sub>1</sub> and
    *n*<sub>2</sub> respectively such that
    *n*<sub>1</sub> + *n*<sub>2</sub> = *n* is the total length of the
    time series. `n1` is *n*<sub>1</sub>. *n*<sub>2</sub> will be
    calculated form `n1` and `n`. The true location of the CP is defined
    as $n\_1 + \\frac12$.
-   `distri1`, `distri2`: The distributions used to sample points in the
    two parts of the TS.
-   only `kind="birth_death"`
    -   `rate`: rate for the exponential distribution of the lifetime of
        points
    -   `m`: fix number of point at every time point
-   only `kind="iid"`
    -   `m1`, `m2`: number of points for the two parts of the TS

### 2.1.1 Distributions

The argument `distri1` and `distri2` name a distribution that is created
from the images in a subfolder named `img` of the working directory
(this path can be changed via `set_image_folder_path("new/path/")`, see
also `get_image_folder_path()`). A distribution (on ℝ<sup>2</sup>) is
created from an image by a mixture of Gaussians on each pixel, where the
individual Gaussian is weighted by the lightness of the pixel. The
variance of each Gaussian is set to (numberOfPixels)<sup> − 1</sup>. As
we do not want to just find changes in mean or variance, these
distributions are then normalized to have zero mean and an identity
covariance matrix.

``` r
distri_names <- c("dot", "8")
distributions <- create_distris(distri_names)
n <- 1e4
for (distri_name in distri_names) {
  points <- distributions[[distri_name]](n)
  rng <- max(abs(points))[c(1,1)] * c(-1,1)
  par(mar=c(2,2,2,2))
  plot(NA, xlim=rng, ylim=rng, asp=1)
  points(points, col=2, pch=".")
  grid()
  title(distri_name)
}
```

![](README_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->![](README_files/figure-gfm/unnamed-chunk-5-2.png)<!-- -->

### 2.1.2 Kind of Sampling Method

-   `"iid"`

Sample `n1` times `m1` points form `distri1` (independently) to obtain
**X**<sub>1</sub>, …, **X**<sub>*n*<sub>1</sub></sub>. Then sample
`n-n1` times `m2` points form `distri2` (independently) to obtain
**X**<sub>*n*<sub>1</sub> + 1</sub>, …, **X**<sub>*n*</sub>. All
**X**<sub>*i*</sub> are independent. Inside each of the two parts of the
time series, the **X**<sub>*i*</sub> are iid. Thus, also the resulting
persistence diagrams are independent / iid.

-   `"birth_death"`

Each **X**<sub>*i*</sub> consists of `m` points which were sampled from
`distri1` or `distri2`. To create **X**<sub>1</sub>, sample `m` points
from `distri1` (independently) and for each point sample a lifetime from
an exponential distribution (independently) with rate `rate`. Given
**X**<sub>*i*</sub> with lifetimes associated to each point, create
**X**<sub>*i* + 1</sub> as follows: Reduce the lifetimes of each point
by 1. Delete all points which now have a negative lifetime. Sample new
points to again have `m` points in total. Reset the lifetimes of the new
points to newly sampled values from the exponential distribution with
rate `rate`. For *i* ≤ *n*<sub>1</sub> + 0.5 − 0.5 \* (1/rate) `distri1`
is used. If *i* is larger, `distri2` is used. This should make the “mean
change point” *n*<sub>1</sub> + 0.5.

## 2.2 Estimators

``` r
# e_name must be unique
estimators <- tribble(
  ~e_name, ~prepro, ~select, ~trim, ~postpro, ~detector,
  "IYG_PC1_D0", "IYG_D0", "1",   0, "cusum", "argmax",
  "IYG_PC3_D0", "IYG_D0", "1:3", 0, "cusum", "argmax",
  "IYG_PC1_D1", "IYG_D1", "1",   0, "cusum", "argmax",
  "IYG_PC3_D1", "IYG_D1", "1:3", 0, "cusum", "argmax",
  "LT_D0",      "LT_D0",  "",    0, "cusum", "argmax",
  "LT_D1",      "LT_D1",  "",    0, "cusum", "argmax",
  "FR_M_D0",    "FR_RM_D0", "1", 20, "id", "argmax",
  "FR_MV_D0",   "FR_RM_D0", "2", 20, "id", "argmax",
  "FR_V_D0",    "FR_RM_D0", "3", 20, "id", "argmax",
  "FR_VV_D0",   "FR_VV_D0", "",  20, "id", "argmax",
  "FR_M_D1",    "FR_RM_D1", "1", 20, "id", "argmax",
  "FR_MV_D1",   "FR_RM_D1", "2", 20, "id", "argmax",
  "FR_V_D1",    "FR_RM_D1", "3", 20, "id", "argmax",
  "FR_VV_D1",   "FR_VV_D1", "",  20, "id", "argmax")
estimators
## # A tibble: 14 x 6
##    e_name     prepro   select  trim postpro detector
##    <chr>      <chr>    <chr>  <dbl> <chr>   <chr>   
##  1 IYG_PC1_D0 IYG_D0   "1"        0 cusum   argmax  
##  2 IYG_PC3_D0 IYG_D0   "1:3"      0 cusum   argmax  
##  3 IYG_PC1_D1 IYG_D1   "1"        0 cusum   argmax  
##  4 IYG_PC3_D1 IYG_D1   "1:3"      0 cusum   argmax  
##  5 LT_D0      LT_D0    ""         0 cusum   argmax  
##  6 LT_D1      LT_D1    ""         0 cusum   argmax  
##  7 FR_M_D0    FR_RM_D0 "1"       20 id      argmax  
##  8 FR_MV_D0   FR_RM_D0 "2"       20 id      argmax  
##  9 FR_V_D0    FR_RM_D0 "3"       20 id      argmax  
## 10 FR_VV_D0   FR_VV_D0 ""        20 id      argmax  
## 11 FR_M_D1    FR_RM_D1 "1"       20 id      argmax  
## 12 FR_MV_D1   FR_RM_D1 "2"       20 id      argmax  
## 13 FR_V_D1    FR_RM_D1 "3"       20 id      argmax  
## 14 FR_VV_D1   FR_VV_D1 ""        20 id      argmax
```

Estimators consist of following parts:

0.  A unique name `e_name`.
1.  A Pre-Processor `prepro`.
2.  Selection `select` and Trimming `trim`.
3.  A Post-Processor `postpro`.
4.  A Detector `detector`.

The Pre-Processor maps from the three lists of a) point sets, b)
persistence diagrams, and c) distance matrices to a matrix
ℝ<sup>ℓ × *q*</sup> which represents a TS of length ℓ (typically ℓ = *n*
or ℓ is close to *n*) of points of dimension *q*, e.g., a
*q*-dimensional statistic calculated form each persistence diagram. The
entry in the table is a character value that refers to an implemented
pre-processor, see below.

Via a selection, one can select *q*<sup>′</sup> columns of the *q*
columns of the output of the pre-processor. The entry in the table is a
character value that can be evaluated as R-code. An empty string selects
all *q* columns.

With `trim`, one can trim the output of the pre-processor from both ends
of the TS. For a trim value *t*, the TS is reduced to length
ℓ<sup>′</sup> = ℓ − 2*t*.

The Post-Processor takes the (ℓ<sup>′</sup> × *q*<sup>′</sup>)-matrix
and maps it to a one-dimensional TS of length *k*. The trivial
post-processor `id` maps
ℝ<sup>ℓ<sup>′</sup> × 1</sup> → ℝ<sup>ℓ<sup>′</sup></sup> via the
identity such that *k* = ℓ<sup>′</sup>. The `cusum` post-processor
applies a *cusum*-like transformation to each column and takes the
squared Euclidean norm in each row. Here, we have
*k* = ℓ<sup>′</sup> − 1.

The Detector takes the one-dimensional TS of length *k* and the original
length *n* and maps them to a value in \[1, *n*\], the estimated
location of the change point. This is usually done by an arg max  in
combination with centering, i.e., adding $\\frac12(n-k)$, to account for
*k* ≠ *n*.

Following functions list the names of all available pre-/post-processors
and detectors.

``` r
get_prepro_names()
## [1] "FR_RM_D0" "FR_RM_D1" "FR_VV_D0" "FR_VV_D1" "IYG_D0"   "IYG_D1"   "LT_D0"   
## [8] "LT_D1"
get_postpro_names()
## [1] "id"    "cusum"
get_detector_names()
## [1] "argmax"
```

### 2.2.1 Cusum Post-Processor

*f*: ℝ<sup>*k* + 1</sup> → ℝ<sup>*k*</sup>, *x* ↦ *z* is applied to each
column of the input *X* ∈ ℝ<sup>ℓ<sup>′</sup> × *q*<sup>′</sup></sup>.
Here *y*<sub>*i*</sub> = *x*<sub>*i*</sub> − *x̄*, for *i* = 1, …, *k*,
where *x̄* is the arithmetic mean of the values of the vector *x*. And
$z\_j = \\sqrt{\\frac{k+1}{j(k+1-j)}}\\sum\_{i=1}^j y\_i$, for
*j* = 1, …, *k*. This results in a matrix
*Z* ∈ ℝ<sup>*k* × *q*<sup>′</sup></sup>. The output of the
`cusum`-post-processor is the vector in ℝ<sup>*k*</sup> of squared
Euclidean norms of the rows of *Z*.

### 2.2.2 Implemented Pre-Processor

-   `IYG_D0`, `IYG_D1`: Method by [Islambekov, Yuvaraj,
    Gel](https://arxiv.org/abs/1910.12939) applied to dimension 0 or
    dimension 1 features
-   `LT_D0`, `LT_D1`: extract the lifetimes (`death - birth`) of
    features from persistence diagram (dim 0 and 1) and order them by
    length; fill with 0s if number of features is not constant in TS.
-   `FR_RM_D0`, `FR_RM_D1`: 3 statistics based on Fréchet change point
    detection that require the calculation of a Fréchet mean, but are
    implemented using the restricted Fréchet mean, i.e., the sample
    element that minimizes the Fréchet objective function (instead of
    the minimizer among all elements in the whole metric space of
    persistence diagrams).
    -   column 1: squared distance of Fréchet means of PDs
        *P*<sub>1</sub>, …, *P*<sub>*i*</sub> and
        *P*<sub>*i* + 1</sub>, …, *P*<sub>*n*</sub>
    -   column 2: mean statistic form [Dubey,
        Müller](https://arxiv.org/abs/1911.11864)
    -   column 3: variance statistic form [Dubey,
        Müller](https://arxiv.org/abs/1911.11864)
-   `FR_VV_D1`, `FR_VV_D1`: A Fréchet change point detection method
    based on the Fréchet variance. Similar to column 3 of `FR_RM_Dx`,
    but does not require calculation of Fréchet means, as it is based on
    $\\mathbf{Var}(X) = \\frac12\\mathbf E\[(X-X^\\prime)^2\]$ for
    independent copy *X*<sup>′</sup> of *X*, instead of
    **V****a****r**(*X*) = **E**\[(*X* − **E**\[*X*\])<sup>2</sup>\].

## 2.3 Required Dists

``` r
# d_name must be unique
required_dists <- tribble(
  ~d_name, ~power, ~dim,
  "0-Inf", Inf, 0,
  "1-Inf", Inf, 1)
required_dists
## # A tibble: 2 x 3
##   d_name power   dim
##   <chr>  <dbl> <dbl>
## 1 0-Inf    Inf     0
## 2 1-Inf    Inf     1
```

This table must contain all distances that are used by one of the
estimators. The estimators may refer to a distance matrix via the name
`d_name`. The available distances are Wasserstein distances of power
`power` (1 to ∞ `Inf`) applied to features of dimension `dim` of the
persistence diagrams.

## 2.4 Save and Load

The three arguments are all tables that can be saved and loaded as
CSV-files. See, e.g., the functions `readr::read_csv()` and
`readr::write_csv()`.

Simulations may take hours to complete. Thus, it recommended to always
save the results (and how the results were created, i.e., the arguments
of `simulate_ts()`) in a file. For example, following code saves results
and arguments as an RDS-file with a file name that contains a time stamp
(requires `library(tidyverse)`).

``` r
# results <- simulate_ts(samplers, estimators, required_dists)
write_rds(
  list(results = results,
       samplers = samplers,
       estimators = estimators,
       required_dists = required_dists),
  paste0("sim_", format(Sys.time(), "%Y%m%d-%H%M%S"), ".RDS"))
```

To load these elements and add them to the global environment, commands
similar to the following can be used.

``` r
lst <- read_rds("sim_example.RDS")
list2env(lst, rlang::global_env())
## <environment: R_GlobalEnv>
```

# 3 Evaluate

`simulate_ts()` returns a table with one row for each repetition of each
sampler and each estimator. The columns are

-   `e_name`: name of the estimator
-   `s_name`: name of the sampler
-   `rep_nr`: number of repetition
-   `estimation`: the estimated location of the change point
-   `postpro_ts`: the time series of length *k* returned by the
    post-processor
-   `len`: *k*

``` r
results
## # A tibble: 84 x 6
##    s_name e_name     rep_nr estimation postpro_ts    len
##    <chr>  <chr>       <int>      <dbl> <list>      <int>
##  1 o-.    IYG_PC1_D0      1      192.  <dbl [299]>   299
##  2 o-.    IYG_PC3_D0      1      190.  <dbl [299]>   299
##  3 o-.    IYG_PC1_D1      1       10.5 <dbl [299]>   299
##  4 o-.    IYG_PC3_D1      1       10.5 <dbl [299]>   299
##  5 o-.    LT_D0           1      196.  <dbl [299]>   299
##  6 o-.    LT_D1           1      164.  <dbl [299]>   299
##  7 o-.    FR_M_D0         1      192.  <dbl [259]>   259
##  8 o-.    FR_MV_D0        1      194.  <dbl [259]>   259
##  9 o-.    FR_V_D0         1       59.5 <dbl [259]>   259
## 10 o-.    FR_VV_D0        1      216.  <dbl [259]>   259
## # ... with 74 more rows
```

## 3.1 Relative Mean Absolute Error

The *Mean Absolute Error (MAE)* of a simulation with *J* repetitions, a
change point at *s* = *n*<sub>1</sub> + 0.5, and estimated CPs
*ŝ*<sub>*j*</sub> is
$\\mathsf{MAE} := \\frac1J \\sum\_{j=1}^J \|\\hat s\_j - s\|$. If one
guesses a location of the CP uniformly at random on a support
{1, …, *k*} with centering for a original time series of length *n*, the
MAE is
$\\mathsf{MAE}\_{\\mathsf{guess}} := \\frac1k \\sum\_{j=1}^k \|(j + \\frac12(n-k)) - s\|$.
The *Relative Mean Absolute Error (RMAE)* is defined as
RMAE := MAE/MAE<sub>guess</sub>.

The function `rmae_eval()` calculates the RMAE values from the
simulation results.

``` r
rmae <- rmae_eval(results, samplers, estimators)
rmae
## # A tibble: 28 x 3
##    s_name e_name       rmae
##    <chr>  <chr>       <dbl>
##  1 o-.    FR_M_D0    0.157 
##  2 o-.    FR_M_D1    0.533 
##  3 o-.    FR_MV_D0   0.0806
##  4 o-.    FR_MV_D1   0.524 
##  5 o-.    FR_V_D0    1.50  
##  6 o-.    FR_V_D1    0.909 
##  7 o-.    FR_VV_D0   0.448 
##  8 o-.    FR_VV_D1   0.228 
##  9 o-.    IYG_PC1_D0 0.0481
## 10 o-.    IYG_PC1_D1 1.16  
## # ... with 18 more rows
```

## 3.2 Visualization

### 3.2.1 Plot Summary

``` r
plot_result_summary(results, samplers, "o-.", "FR_M_D0")
```

![](README_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

-   x-axis: times *t* ∈ {1, …, *n*} of the original time series
-   y-axis: values of the output of the post-processor
-   red: mean post-processor values across all repetitions with
    quartiles (25% and 75% quantiles)
-   green: mean estimated CP with quartiles
-   blue: true CP
-   dashed black lines indicate the support of the post-processed time
    series

### 3.2.2 Plot All Runs

``` r
plot_result_all(results, samplers, "o-.", "FR_M_D0")
```

![](README_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

## 3.3 Pretty Tables

The functions `pretty_table()` and `pretty_color_table()` allow to view
the results in pretty tables.

``` r
rmae %>% 
  group_by(e_name) %>% 
  summarise(sum_rmae = sum(rmae), .groups="drop") %>% 
  mutate(rank = rank(sum_rmae, ties.method="random")) %>% 
  arrange(rank) %>% 
  select(rank, everything()) %>% 
  pretty_table()
```

<table class="table table-striped table-hover table-condensed" style="width: auto !important; margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:right;">
rank
</th>
<th style="text-align:left;">
e\_name
</th>
<th style="text-align:right;">
sum\_rmae
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:right;">
1
</td>
<td style="text-align:left;">
LT\_D0
</td>
<td style="text-align:right;">
0.4732663
</td>
</tr>
<tr>
<td style="text-align:right;">
2
</td>
<td style="text-align:left;">
FR\_M\_D1
</td>
<td style="text-align:right;">
0.7347518
</td>
</tr>
<tr>
<td style="text-align:right;">
3
</td>
<td style="text-align:left;">
FR\_M\_D0
</td>
<td style="text-align:right;">
0.8467566
</td>
</tr>
<tr>
<td style="text-align:right;">
4
</td>
<td style="text-align:left;">
FR\_VV\_D1
</td>
<td style="text-align:right;">
0.8467566
</td>
</tr>
<tr>
<td style="text-align:right;">
5
</td>
<td style="text-align:left;">
FR\_MV\_D0
</td>
<td style="text-align:right;">
0.8915586
</td>
</tr>
<tr>
<td style="text-align:right;">
6
</td>
<td style="text-align:left;">
LT\_D1
</td>
<td style="text-align:right;">
1.1270154
</td>
</tr>
<tr>
<td style="text-align:right;">
7
</td>
<td style="text-align:left;">
IYG\_PC1\_D0
</td>
<td style="text-align:right;">
1.5601744
</td>
</tr>
<tr>
<td style="text-align:right;">
8
</td>
<td style="text-align:left;">
FR\_VV\_D0
</td>
<td style="text-align:right;">
1.5770282
</td>
</tr>
<tr>
<td style="text-align:right;">
9
</td>
<td style="text-align:left;">
IYG\_PC3\_D0
</td>
<td style="text-align:right;">
1.5882495
</td>
</tr>
<tr>
<td style="text-align:right;">
10
</td>
<td style="text-align:left;">
FR\_MV\_D1
</td>
<td style="text-align:right;">
1.5904688
</td>
</tr>
<tr>
<td style="text-align:right;">
11
</td>
<td style="text-align:left;">
FR\_V\_D0
</td>
<td style="text-align:right;">
1.8727210
</td>
</tr>
<tr>
<td style="text-align:right;">
12
</td>
<td style="text-align:left;">
FR\_V\_D1
</td>
<td style="text-align:right;">
2.3521017
</td>
</tr>
<tr>
<td style="text-align:right;">
13
</td>
<td style="text-align:left;">
IYG\_PC1\_D1
</td>
<td style="text-align:right;">
2.3863850
</td>
</tr>
<tr>
<td style="text-align:right;">
14
</td>
<td style="text-align:left;">
IYG\_PC3\_D1
</td>
<td style="text-align:right;">
2.8155332
</td>
</tr>
</tbody>
</table>

``` r
rmae %>% 
  mutate(v = sprintf("%.2f", rmae)) %>% 
  select(s_name, e_name, v) %>% 
  pivot_wider(names_from=s_name, values_from=v) ->
  tb
rmae %>% 
  select(s_name, e_name, rmae) %>%
  pivot_wider(names_from=s_name, values_from=rmae) %>% 
  select(-e_name)->
  color_tb
pretty_color_table(tb, 2:ncol(tb), color_tb)
```

<table class="table table-striped table-hover table-condensed" style="width: auto !important; margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;">
e\_name
</th>
<th style="text-align:left;">
o-.
</th>
<th style="text-align:left;">
o-n
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
FR\_M\_D0
</td>
<td style="text-align:left;background-color: #A9FF80 !important;">
0.16
</td>
<td style="text-align:left;background-color: #FFCC80 !important;">
0.69
</td>
</tr>
<tr>
<td style="text-align:left;">
FR\_M\_D1
</td>
<td style="text-align:left;background-color: #FFF680 !important;">
0.53
</td>
<td style="text-align:left;background-color: #B3FF80 !important;">
0.20
</td>
</tr>
<tr>
<td style="text-align:left;">
FR\_MV\_D0
</td>
<td style="text-align:left;background-color: #94FF80 !important;">
0.08
</td>
<td style="text-align:left;background-color: #FFAE80 !important;">
0.81
</td>
</tr>
<tr>
<td style="text-align:left;">
FR\_MV\_D1
</td>
<td style="text-align:left;background-color: #FFF880 !important;">
0.52
</td>
<td style="text-align:left;background-color: #FF8080 !important;">
1.07
</td>
</tr>
<tr>
<td style="text-align:left;">
FR\_V\_D0
</td>
<td style="text-align:left;background-color: #FF8080 !important;">
1.50
</td>
<td style="text-align:left;background-color: #E1FF80 !important;">
0.38
</td>
</tr>
<tr>
<td style="text-align:left;">
FR\_V\_D1
</td>
<td style="text-align:left;background-color: #FF9480 !important;">
0.91
</td>
<td style="text-align:left;background-color: #FF8080 !important;">
1.44
</td>
</tr>
<tr>
<td style="text-align:left;">
FR\_VV\_D0
</td>
<td style="text-align:left;background-color: #F3FF80 !important;">
0.45
</td>
<td style="text-align:left;background-color: #FF8080 !important;">
1.13
</td>
</tr>
<tr>
<td style="text-align:left;">
FR\_VV\_D1
</td>
<td style="text-align:left;background-color: #BBFF80 !important;">
0.23
</td>
<td style="text-align:left;background-color: #FFDE80 !important;">
0.62
</td>
</tr>
<tr>
<td style="text-align:left;">
IYG\_PC1\_D0
</td>
<td style="text-align:left;background-color: #8CFF80 !important;">
0.05
</td>
<td style="text-align:left;background-color: #FF8080 !important;">
1.51
</td>
</tr>
<tr>
<td style="text-align:left;">
IYG\_PC1\_D1
</td>
<td style="text-align:left;background-color: #FF8080 !important;">
1.16
</td>
<td style="text-align:left;background-color: #FF8080 !important;">
1.23
</td>
</tr>
<tr>
<td style="text-align:left;">
IYG\_PC3\_D0
</td>
<td style="text-align:left;background-color: #97FF80 !important;">
0.09
</td>
<td style="text-align:left;background-color: #FF8080 !important;">
1.50
</td>
</tr>
<tr>
<td style="text-align:left;">
IYG\_PC3\_D1
</td>
<td style="text-align:left;background-color: #FF8080 !important;">
2.15
</td>
<td style="text-align:left;background-color: #FFD280 !important;">
0.67
</td>
</tr>
<tr>
<td style="text-align:left;">
LT\_D0
</td>
<td style="text-align:left;background-color: #8FFF80 !important;">
0.06
</td>
<td style="text-align:left;background-color: #E9FF80 !important;">
0.41
</td>
</tr>
<tr>
<td style="text-align:left;">
LT\_D1
</td>
<td style="text-align:left;background-color: #FFC580 !important;">
0.72
</td>
<td style="text-align:left;background-color: #E9FF80 !important;">
0.41
</td>
</tr>
</tbody>
</table>

``` r
rmae %>% 
  group_by(s_name) %>% 
  mutate(rank = rank(rmae, ties.method="random")) %>% 
  arrange(rank) %>% 
  ungroup() ->
  x
x %>% 
  mutate(v = sprintf("%.2f %s", rmae, e_name)) %>% 
  select(rank, s_name, v) %>% 
  pivot_wider(
    names_from = s_name, 
    values_from = v) ->
  tb
x %>% 
  select(rank, s_name, rmae) %>% 
  pivot_wider(
    names_from = s_name, 
    values_from = rmae) %>% 
  select(-rank) ->
  color_tb
pretty_color_table(tb, 2:ncol(tb), color_tb)
```

<table class="table table-striped table-hover table-condensed" style="width: auto !important; margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:right;">
rank
</th>
<th style="text-align:left;">
o-.
</th>
<th style="text-align:left;">
o-n
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:right;">
1
</td>
<td style="text-align:left;background-color: #8CFF80 !important;">
0.05 IYG\_PC1\_D0
</td>
<td style="text-align:left;background-color: #B3FF80 !important;">
0.20 FR\_M\_D1
</td>
</tr>
<tr>
<td style="text-align:right;">
2
</td>
<td style="text-align:left;background-color: #8FFF80 !important;">
0.06 LT\_D0
</td>
<td style="text-align:left;background-color: #E1FF80 !important;">
0.38 FR\_V\_D0
</td>
</tr>
<tr>
<td style="text-align:right;">
3
</td>
<td style="text-align:left;background-color: #94FF80 !important;">
0.08 FR\_MV\_D0
</td>
<td style="text-align:left;background-color: #E9FF80 !important;">
0.41 LT\_D1
</td>
</tr>
<tr>
<td style="text-align:right;">
4
</td>
<td style="text-align:left;background-color: #97FF80 !important;">
0.09 IYG\_PC3\_D0
</td>
<td style="text-align:left;background-color: #E9FF80 !important;">
0.41 LT\_D0
</td>
</tr>
<tr>
<td style="text-align:right;">
5
</td>
<td style="text-align:left;background-color: #A9FF80 !important;">
0.16 FR\_M\_D0
</td>
<td style="text-align:left;background-color: #FFDE80 !important;">
0.62 FR\_VV\_D1
</td>
</tr>
<tr>
<td style="text-align:right;">
6
</td>
<td style="text-align:left;background-color: #BBFF80 !important;">
0.23 FR\_VV\_D1
</td>
<td style="text-align:left;background-color: #FFD280 !important;">
0.67 IYG\_PC3\_D1
</td>
</tr>
<tr>
<td style="text-align:right;">
7
</td>
<td style="text-align:left;background-color: #F3FF80 !important;">
0.45 FR\_VV\_D0
</td>
<td style="text-align:left;background-color: #FFCC80 !important;">
0.69 FR\_M\_D0
</td>
</tr>
<tr>
<td style="text-align:right;">
8
</td>
<td style="text-align:left;background-color: #FFF880 !important;">
0.52 FR\_MV\_D1
</td>
<td style="text-align:left;background-color: #FFAE80 !important;">
0.81 FR\_MV\_D0
</td>
</tr>
<tr>
<td style="text-align:right;">
9
</td>
<td style="text-align:left;background-color: #FFF680 !important;">
0.53 FR\_M\_D1
</td>
<td style="text-align:left;background-color: #FF8080 !important;">
1.07 FR\_MV\_D1
</td>
</tr>
<tr>
<td style="text-align:right;">
10
</td>
<td style="text-align:left;background-color: #FFC580 !important;">
0.72 LT\_D1
</td>
<td style="text-align:left;background-color: #FF8080 !important;">
1.13 FR\_VV\_D0
</td>
</tr>
<tr>
<td style="text-align:right;">
11
</td>
<td style="text-align:left;background-color: #FF9480 !important;">
0.91 FR\_V\_D1
</td>
<td style="text-align:left;background-color: #FF8080 !important;">
1.23 IYG\_PC1\_D1
</td>
</tr>
<tr>
<td style="text-align:right;">
12
</td>
<td style="text-align:left;background-color: #FF8080 !important;">
1.16 IYG\_PC1\_D1
</td>
<td style="text-align:left;background-color: #FF8080 !important;">
1.44 FR\_V\_D1
</td>
</tr>
<tr>
<td style="text-align:right;">
13
</td>
<td style="text-align:left;background-color: #FF8080 !important;">
1.50 FR\_V\_D0
</td>
<td style="text-align:left;background-color: #FF8080 !important;">
1.50 IYG\_PC3\_D0
</td>
</tr>
<tr>
<td style="text-align:right;">
14
</td>
<td style="text-align:left;background-color: #FF8080 !important;">
2.15 IYG\_PC3\_D1
</td>
<td style="text-align:left;background-color: #FF8080 !important;">
1.51 IYG\_PC1\_D0
</td>
</tr>
</tbody>
</table>

# 4 Further Notes

## 4.1 Parallel

Use the argument `n_cores` to enable parallel runs of repetitions. The
default is `n_cores=1`, which runs the code non-parallel. There is some
overhead when using multiple cores. Thus, `n_cores=2` might not always
be faster than `n_cores=1` and is always less than twice as fast.

``` r
results <- simulate_ts(
  samplers, estimators, required_dists, 
  n_cores=parallel::detectCores())
```

## 4.2 Reproducibility

Set a random seed via `set.seed()` to make simulations reproducible.
This should also work with parallel execution, i.e., with an argument
`n_cores` greater than 1, but was not extensively tested.

``` r
set.seed(1)
results <- simulate_ts(samplers, estimators, required_dists)
```

# 5 Extending the Package

## 5.1 New Pre-Processors

In file `def_prepros.R` add a named function object to the list
`prepros`. The function must take the three arguments
`point_ts, pd_ts, dist_mats` and return a matrix. Alternatively, one can
add Pre-Processors without changing the package code using the function
`register_prepro(name, fun)`, e.g.:

``` r
f <- function(point_ts, pd_ts, dist_mats) {
  t(sapply(point_ts, rowMeans))
}
register_prepro("EUCL_MEAN", f)
```

## 5.2 New Post-Processors

In file `def_postpros.R` add a named function object to the list
`postpros`. The function must one argument – a
(ℓ<sup>′</sup> × *q*<sup>′</sup>)-matrix and return a vector.
Alternatively, one can add Post-Processors without changing the package
code using the function `register_postpro(name, fun)`, e.g.:

``` r
f <- function(X) {
  rowMeans(X^2)
}
register_postpro("norm_2_sq", f)
```

## 5.3 New Detectors

In file `def_dectectors.R` add a named function object to the list
`detectors`. The function must 2 argument – a vector representing the
post-processed TS and the length of the original TS *n* as a single
integer – and return a single `double` (real value) – the estimated
location of the CP (should be in \[1, *n*\]. Alternatively, one can add
Detectors without changing the package code using the function
`register_detector(name, fun)`, e.g.:

``` r
f <- function(y, n) {
  which.min(y) + (n-length(y))/2
}
register_detector("argmin", f)
```

## 5.4 New Image Distributions

Just create new PNGs and put them in a subfolder `img` of your working
directory. They can be referred to by their filename (without ending
`.png`).

## 5.5 Other Distributions

Currently, no other types of distributions are implemented. But the
function `create_distris()` (file `distris.R`) has an argument `type`
which allows for such extensions in the future. Then the call of
`create_distris()` inside of `sample_ts()` might also need change.

## 5.6 New Kinds of Samplers

Change `sample_point_ts()` in file `simulate.R`.
