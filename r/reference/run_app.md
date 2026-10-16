# Launch the pilotr no-code app

Runs the bundled Shiny app locally. When the package is installed, the
app calls the package's functions directly. Running locally gives each
user a private R process, so no one else's work waits on theirs.

## Usage

``` r
run_app(..., async = NULL)
```

## Arguments

- ...:

  passed to
  [`shiny::runApp()`](https://rdrr.io/pkg/shiny/man/runApp.html) (e.g.
  `port`, `launch.browser`).

- async:

  If TRUE (the default when 'future' and 'promises' are installed), the
  app runs each power analysis and power curve in a background R
  process. That leaves the app's own R process free for other sessions,
  such as a second browser tab, while the session that started the run
  waits for it. A parallel future plan you have set is used as it
  stands. Over a sequential plan, the default, `run_app()` sets a
  multisession plan of two workers while the app runs. Either way, the
  plan in force before the call is in force again once the app stops.
  With FALSE, power runs in the app's own R process and no plan is set.

## Value

No return value; called for its side effect of launching the Shiny
application, which blocks the R session until the app is closed.

## Details

Each power analysis and power curve in the app runs in one R process,
and the session that asked for it waits for the result. Its controls and
outputs do not update until the run ends. For a heavy run, download the
specification and call
[`power_design()`](https://pablobernabeu.github.io/pilotr/r/reference/power_design.md)
or
[`power_mixed()`](https://pablobernabeu.github.io/pilotr/r/reference/power_mixed.md)
from the console, where `workers` spreads the replicates over your local
cores.

## Examples

``` r
if (FALSE) { # \dontrun{
run_app()
} # }
```
