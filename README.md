# MDEqp: Estimating the Sample Mean and Standard Deviation from Commonly Reported Quantiles

`MDEqp` estimates the sample mean and standard deviation from commonly reported quantiles and sample size.

The package is intended for clinicians, healthcare researchers, and evidence synthesis analysts who need to convert reported summary statistics into estimated sample means and standard deviations. This problem commonly arises in meta-analysis, where many studies report medians, interquartile ranges (IQRs), or other percentiles. `MDEqp` addresses the same general problem setting considered by McGrath et al. (2020), while implementing minimum-distance and weighted quantile/probability estimation approaches.

Specifically, these methods can be applied to studies that report one of the following sets of summary statistics:

-   `S1`: median, minimum and maximum values, and sample size
-   `S2`: median, first and third quartiles, and sample size
-   `S3`: median, minimum and maximum values, first and third quartiles, and sample size

## Installation

You can install the development version of `MDEqp` from GitHub with:

``` r
install.packages("remotes")  # Only needed once
remotes::install_github("luvuu1/MDEqp")
```

After the package is available on CRAN, it can be installed using:

``` r
install.packages("MDEqp")
```

## Usage

`MDEqp` provides four main estimation functions:

-   `MDE()`: quantile-scale minimum distance estimation
-   `MDEp()`: probability-scale minimum distance estimation
-   `wQE()`: weighted quantile estimation
-   `wPE()`: weighted probability estimation

Candidate distributions include the normal, log-normal, gamma, Weibull, and beta distributions. For all four methods, candidate distributions are fitted using the method-specific objective. The final selected distribution is then chosen using the smallest unweighted quantile-scale squared discrepancy between the reported quantiles and fitted quantiles.

Because different reporting scenarios require different summary statistics, users should supply arguments by name. For example, studies reporting only the first quartile, median, third quartile, and sample size can be entered as `MDEp(q1 = 3, med = 5, q3 = 7, n = 100)`; `min` and `max` are not required for this scenario.

## Examples

### Scenario S1: median, minimum, maximum, and sample size

``` r
library(MDEqp)

res_s1 <- MDE(
  min = 1,
  med = 5,
  max = 10,
  n = 50
)

res_s1
summary(res_s1)
```

### Scenario S2: median, first and third quartiles, and sample size

``` r
res_s2 <- MDEp(
  q1 = 3,
  med = 5,
  q3 = 7,
  n = 100
)

res_s2
summary(res_s2)
```

### Scenario S3: median, minimum, maximum, first and third quartiles, and sample size

``` r
res_s3 <- wQE(
  min = 1,
  q1 = 3,
  med = 5,
  q3 = 7,
  max = 10,
  n = 50
)

res_s3
summary(res_s3)
```

## Zero-based bounded outcomes

For zero-based bounded outcomes with a known theoretical upper bound, use the `range` argument. In the current implementation, the lower bound is assumed to be 0, so `range` represents the upper endpoint of the theoretical scale, not the observed sample range.

When `range` is supplied, all reported summary statistics used for fitting must fall within `[0, range]`. These summaries are divided by `range` before model fitting, so all candidate distributions are fitted on the normalized scale. The estimated mean and standard deviation are then transformed back to the original scale. For example, for PHQ-9 scores ranging from 0 to 27, use `range = 27`. For an outcome scale ranging from 0 to 100:

``` r
res_bounded <- wPE(
  min = 5,
  q1 = 20,
  med = 35,
  q3 = 50,
  max = 80,
  n = 100,
  range = 100
)

res_bounded
summary(res_bounded)
```

The beta distribution is fitted only when the quantiles used for fitting are compatible with the `[0, 1]` scale. This can happen either because the original outcome is already on this scale, or because a zero-based theoretical `range` is supplied and the summaries are normalized.

## Output

The returned object includes:

-   `est.mean`: estimated sample mean from the selected distribution
-   `est.sd`: estimated sample standard deviation from the selected distribution
-   `selected.dist`: selected candidate distribution
-   `valuesls`: selection statistic values for candidate distributions
-   `scenario`: detected reporting scenario
-   `fitted.dists`: distribution-specific fitting results

Use `summary()` to compare the estimated mean and standard deviation across all candidate distributions.

``` r
res <- MDE(
  min = 1,
  q1 = 3,
  med = 5,
  q3 = 7,
  max = 10,
  n = 50
)

summary(res)
```

## Handling invalid estimates

If a candidate distribution cannot be evaluated, produces non-finite mean or SD estimates, or gives a negative SD, both estimates are reported as `NA`.

The estimated mean must fall within `[0, range]` when `range` is supplied, or within the reported sample minimum and maximum for scenarios `S1` and `S3` without `range`. If the mean falls outside the applicable interval, both estimates are reported as `NA`. No interval check is applied to the mean for scenario `S2` without `range`.

Excluded candidates have their selection statistics set to `Inf` and are not selected. If no candidate remains eligible, the selected distribution, estimated mean, and estimated SD are reported as `NA`.

## References

Tang X, Tong T, Zhang X, Chu H. (2026). Minimum distance estimation of mean and standard deviation from reported quantiles. *Research Synthesis Methods*, 17(5), 958–975. https://doi.org/10.1017/rsm.2026.10090

Tang X, Wan Y, Zhu C, Zhang X, Tong T, Pu H, Chu H. (under review). Estimating means and standard deviations from reported quantiles for meta-analysis using probability-scale minimum distance methods.

McGrath S, Zhao X, Steele R, et al. (2020). Estimating the sample mean and standard deviation from commonly reported quantiles in meta-analysis. Statistical Methods in Medical Research, 29(9), 2520–2537. https://doi.org/10.1177/0962280219889080


## Implementation reference

The implementation was developed with reference to the fitting structure and helper functions in the [estmeansd](https://github.com/stmcg/estmeansd) package. MDEqp implements the four methods described above and does not require estmeansd to run.

## License

This package is released under the GNU General Public License,
version 3 or later (GPL >= 3).
