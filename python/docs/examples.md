# Worked examples

Eight ready-to-run specifications ship with the package, covering six of the eight response
families and between-subjects, crossed, nested and partially crossed designs. `pilotr_example()`
lists them and returns the path to each for `load_spec()`, and the same JSON files drive the R
package and the no-code app. Each is simulated below, showing the design it describes and the
first rows of the data it produces. The specification format is documented on the
[Specification](specification.md) page.

```python exec="true" session="ex"
import sys; sys.path.insert(0, "docs")
from _exec import table
```

<!-- The block is given an explicit id because markdown-exec otherwise numbers
     executed blocks with a counter that runs across the whole site, so adding a
     block to any earlier page renumbers the anchors on this one. Each example's
     heading anchor is a link a reader might keep, so it is pinned here instead. -->
```python exec="true" session="ex" id="example"
from pilotr import simulate, load_spec, pilotr_example

DESC = {
    "between_2group_gaussian": "Two-group between-subjects Gaussian.",
    "crossed_mixed_rt":
        "Crossed by-subject and by-item reaction times (shifted lognormal).",
    "beta_proportion": "Bounded proportions through the Beta family.",
    "ordinal_likert_between":
        "Five-point Likert responses via a cumulative-logit model.",
    "poisson_counts_between": "Count outcomes through a log link.",
    "reading_time_continuous":
        "Continuous item and subject predictors and a within factor, with"
        " interactions among them, and a lognormal reading-time outcome.",
    "nested_clusters":
        "Subjects nested in higher-level clusters (an extra grouping factor).",
    "partial_crossing": "Each subject sees a sampled subset of items.",
}

out = []
for name in pilotr_example():
    spec = load_spec(pilotr_example(name))
    fam = spec["response"]["family"]
    d = simulate(spec)
    out.append(f"### {name}")
    out.append(f"*{DESC.get(name, '')}* The `{fam}` family, {len(d)} rows.")
    out.append(table(d.head(4)))
    out.append("")
print("\n".join(out))
```
