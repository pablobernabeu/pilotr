# The no-code app

pilotr ships a point-and-click application over the same design
specification that the R and Python packages consume. It is a thin
client. Every control writes into the portable JSON specification, which
you can download and run unchanged in either package to obtain identical
data.

## Two ways to run it

One way is to run it on your own machine. With the package installed,
launch the app locally:

``` r

pilotr::run_app()
```

Each user runs their own R process, so a heavy power run never blocks
anyone else, and simulations use all of your local cores.

The other way needs nothing installed. A serverless build runs entirely
in the browser through WebAssembly, with no server and nothing uploaded.
Try it at the [no-code
app](https://pablobernabeu.github.io/pilotr/app/). It covers the light
path: build a design, simulate, inspect the data, estimate two-group
Gaussian power and export.

## What the app does

Describe the design in the left panel: the sample sizes, the factor and
its two levels, the fixed intercept and effect, a response family, and
for crossed designs the random-effect standard deviations. Then read the
tabs. ‘Design spec’ holds the exact JSON specification, the portable
source of truth, and ‘Data’ the simulated data set. ‘Summary & plot’
gives group summaries and a plot. ‘Power’ reports simulation-based
two-group Gaussian power with the Type S and Type M errors, together
with a power curve over sample size. ‘R script’ emits a self-contained,
reproducible R script, which the installed app can also verify by
running it in a clean R session and confirming bit-for-bit reproduction.

Changing the response family resets the intercept and effect to sensible
values for that family’s scale, so a point-and-click design stays valid.
The advanced ‘paste a JSON spec’ box accepts designs beyond the
point-and-click controls, such as continuous predictors, interactions
and nesting.

## One specification, three interfaces

The app does not own the design. The specification does. Download the
specification (`.json`) or the data (`.csv`), and run the specification
unchanged in R with
[`simulate_design()`](https://pablobernabeu.github.io/pilotr/r/reference/simulate_design.md)
or in Python with `pilotr.simulate()` for identical data. Crossed
mixed-effects power runs in the installed packages rather than in the
browser (in R via `lme4`, the reference backend), and precision/ROPE
design analysis runs in the R package. Whichever interface you use, the
specification you build in the app drives all three.
