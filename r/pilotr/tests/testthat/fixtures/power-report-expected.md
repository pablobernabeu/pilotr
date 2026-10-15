# Sample-size justification

## Primary analysis

- Model and test: lmerTest::lmer, Satterthwaite t
- Significance level: alpha = 0.05, two-sided
- Simplification rule: maximal, then zcp, then intercepts, moving on from singular fits
- Simulated units: subject 24, item 18

## Desired power and sample size

- Desired power: 0.8
- Sample size: 52 (95% confidence interval 41 to 68), the exact solution 51.37 rounded up
- Fit: a probit regression of `power` on the square root of `n_subject`, over 7 points from 8 to 56
- Heterogeneity factor: 1.31, the dispersion that scales the variance behind the interval, which is 1 where the curve fits the probit model

## Focal effects

The rates below are those of the simulated design, each over the replicates that returned an estimate for the effect.

### `effect`

- True value: 0.06
- Justification: Half the priming effect in the pilot (Müller, 2021), the smallest effect of theoretical interest.
- Power: 0.705 (Monte Carlo standard error 0.032; 95% Wilson interval 0.638 to 0.764), from 141 of 200 replicates
- Type S: 0.000 (95% Wilson interval 0.000 to 0.027)
- Type M: 1.273 (Monte Carlo standard error 0.021)

### `age`

- True value: -0.02
- Justification: not given
- Power: 0.187 (Monte Carlo standard error 0.028; 95% Wilson interval 0.139 to 0.247), from 37 of 198 replicates
- Type S: 0.054 (95% Wilson interval 0.015 to 0.177)
- Type M: 2.100 (Monte Carlo standard error 0.190)

## Power when the null hypothesis is true

With every focal effect set to zero, at alpha = 0.05, two-sided:

- `effect`: 0.055 (Monte Carlo standard error 0.016; 95% Wilson interval 0.031 to 0.096), from 11 of 199 replicates
- `age`: 0.045 (Monte Carlo standard error 0.015; 95% Wilson interval 0.024 to 0.084), from 9 of 199 replicates

## Monte Carlo accounting

- Replicates: 200 simulated, 200 returned an estimate
- Fits: 151 converged cleanly, 46 singular, 3 with a warning
- Analysed at each rung: maximal 0.240, zcp 0.115, intercepts 0.645
- Null condition: 200 simulated, 199 returned an estimate

Each rate's Monte Carlo standard error is sqrt(p(1 - p) / n), over the n replicates that returned an estimate, and its interval is the 95% Wilson score interval. The standard error of Type M is the standard deviation of the exaggeration ratios over the square root of their number.

## Reproducibility

- Software: pilotr 0.3.1, R 4.6.1, lme4 1.1.37, lmerTest 3.1.3
- Specification: priming.json, MD5 checksum 5d41402abc4b2a76b9719d911017c592
- Seed: 1
