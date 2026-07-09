# Worked examples

One ready-to-run specification per design family ships with the package, in the repository's
[`spec/examples/`](https://github.com/pablobernabeu/pilotr/tree/main/spec/examples) directory. The same JSON files drive the R package and the no-code app. Each is
simulated below, showing the design it describes and the first rows of the data it produces. The
specification format is documented on the [Specification](specification.md) page.

```python exec="true" session="ex"
import sys; sys.path.insert(0, "docs")
from _exec import table
```

```python exec="true" session="ex"
import glob
from pilotr import simulate, load_spec

DESC = {
    "between_2group_gaussian": "Two-group between-subjects Gaussian.",
    "crossed_mixed_rt": "Crossed by-subject and by-item reaction times (shifted lognormal).",
    "beta_proportion": "Bounded proportions through the Beta family.",
    "ordinal_likert_between": "Five-point Likert responses via a cumulative-logit model.",
    "poisson_counts_between": "Count outcomes through a log link.",
    "reading_time_continuous": "A continuous predictor with a lognormal reading-time outcome.",
    "nested_clusters": "Subjects nested in higher-level clusters (an extra grouping factor).",
    "partial_crossing": "Each subject sees a sampled subset of items.",
}

out = []
for path in sorted(glob.glob("../spec/examples/*.json")):
    spec = load_spec(path)
    name = spec["name"]
    fam = spec["response"]["family"]
    d = simulate(spec)
    out.append(f"### {name}")
    out.append(f"*{DESC.get(name, '')}* &mdash; `{fam}` family, {len(d)} rows.")
    out.append(table(d.head(4)))
    out.append("")
print("\n".join(out))
```
