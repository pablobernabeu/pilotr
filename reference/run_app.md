# Launch the pilotr no-code app.

Runs the bundled Shiny app locally. When the package is installed, the
app calls the package's functions directly. Running locally gives each
user a private R process, so work is not blocked by a shared process,
and power simulations can be parallelised across the user's own cores.

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

  If TRUE (default when 'future' and 'promises' are installed), set a
  multisession future plan so power runs execute in a background worker
  and do not block the UI. This keeps the app responsive on your own
  machine while a long power run completes in the background.

## Value

No return value; called for its side effect of launching the Shiny
application, which blocks the R session until the app is closed.

## Examples

``` r
if (FALSE) { # \dontrun{
run_app()
} # }
```
