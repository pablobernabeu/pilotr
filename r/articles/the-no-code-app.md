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

Each user runs their own R process, so a power run never holds up anyone
else. Within the app, each power analysis or power curve runs in one R
process, and the app waits for its result. With the future and promises
packages installed,
[`run_app()`](https://pablobernabeu.github.io/pilotr/r/reference/run_app.md)
runs it in a background process, leaving the app free for other
sessions, such as a second browser tab. For a heavy run, download the
specification and call
[`power_design()`](https://pablobernabeu.github.io/pilotr/r/reference/power_design.md)
or
[`power_mixed()`](https://pablobernabeu.github.io/pilotr/r/reference/power_mixed.md)
from the console, whose `workers` argument spreads the replicates over
your local cores.

The other way needs nothing installed. A serverless build runs entirely
in the browser through WebAssembly, with no server and nothing uploaded.
Try it at the [no-code
app](https://pablobernabeu.github.io/pilotr/app/). It covers the light
path: build a design, simulate, inspect the data, estimate two-group
Gaussian power with the sample size its curve solves to, and export.

## What the app does

Describe the design in the left panel: the sample sizes, the factor and
its two levels, the fixed intercept and effect, a response family, and
for crossed designs the random-effect standard deviations. A within
design crossed with items can use counterbalanced lists, which show each
subject each item in one condition and rotate the conditions over
subjects and items, and its specification then declares version 0.4.
Then read the tabs. Both apps open on ‘Guide’, which describes the
controls, the response families and the analyses. Three of the other
tabs have shorter names in the browser build, given here in brackets
after the installed app’s. ‘Design spec (JSON)’ (‘Design spec’) holds
the exact JSON specification, the portable source of truth, and ‘Data’
the simulated data set. ‘Summary & plot’ gives group summaries and a
plot. ‘Power & design analysis’ (‘Power’) reports simulation-based
two-group Gaussian power with the Type S and Type M errors, together
with a power curve over sample size and the sample size at which that
curve reaches 0.80, with a confidence interval on it. The curve adds
sample sizes until it brackets a power of 0.80, and for a design without
random effects or continuous predictors the tab also gives the sample
size [`power.t.test()`](https://rdrr.io/r/stats/power.t.test.html)
computes, as a check. The tab runs this analysis for a two-group
Gaussian design with one row per subject, and for any other design it
explains how to run it in the installed package. ‘Reproducible R script’
(‘R script’) emits a self-contained, reproducible R script, which ends
by suggesting the analysis that fits the design. In the installed app,
Verify runs that script in a clean R session and checks that its data
are identical, bit for bit, to the data the current specification
produces. The comparison is with the specification in force when you
select Verify, so it does not depend on when you last selected Simulate.

Changing the response family resets the intercept and effect to sensible
values for that family’s scale, so a point-and-click design stays valid.
The advanced ‘paste a JSON spec’ box accepts designs beyond the
point-and-click controls, such as continuous predictors, interactions
and nesting.

‘Data’, ‘Summary & plot’ and the data download show the last simulation,
and the ‘Data’ tab and the name of the downloaded file give the seed it
used. The specification, power and script tabs and the specification
download follow the controls as they stand. When the controls change the
design after a simulation, a note above the tabs says so until you
select Simulate again.

## One specification, three interfaces

The app does not own the design. The specification does. Download the
specification (`.json`) or the data (`.csv`), and run the specification
unchanged in R with
[`simulate_design()`](https://pablobernabeu.github.io/pilotr/r/reference/simulate_design.md)
or in Python with `pilotr.simulate()` for the same data, identical to
the last bit wherever R and Python share a maths library (see the
[cross-language
article](https://pablobernabeu.github.io/pilotr/r/articles/cross-language.md)).
The app rounds Gaussian and reaction-time responses to four decimals and
leaves a Beta response unrounded, for the reason
[`?build_spec`](https://pablobernabeu.github.io/pilotr/r/reference/build_spec.md)
gives. The heavier analyses stay with the installed packages: crossed
mixed-effects power in R, through `lme4` as the reference backend, and
in Python, through `statsmodels`, and precision/ROPE design analysis in
R. Whichever interface you use, the specification you build in the app
drives all three.
