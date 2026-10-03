"""Crossed mixed-effects simulation-based power in Python (statsmodels backend).
Demonstrates that the same spec yields mixed-effects power in either ecosystem."""
import os, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
SPEC = os.path.join(HERE, "..", "..", "spec", "examples", "crossed_mixed_rt.json")
from pilotr import power_mixed, load_spec

n_sims = int(os.environ.get("N_SIMS", "40"))
print(f"Fitting {n_sims} crossed MixedLM models in Python (statsmodels)...")
t = time.time()
res = power_mixed(load_spec(SPEC), n_sims=n_sims)
print(f"elapsed {time.time() - t:.0f}s\n")
for k, v in res.items():
    print(f"  {k:14s}: {v}")
print("\nNOTE: the p-values are Wald z tests. With few subjects or items, these reject more")
print("readily than the Satterthwaite tests of the R package's power_mixed(), so in small designs")
print("this power can exceed R's on the same data. R's power_mixed() also fits the random-effect")
print("correlations this specification declares, which the components fitted here leave out.")
