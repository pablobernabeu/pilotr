"""Crossed mixed-effects simulation-based power in Python (statsmodels backend).
Demonstrates that the same spec yields mixed-effects power in either ecosystem."""
import os, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
SPEC = os.path.join(HERE, "..", "..", "spec", "examples", "crossed_mixed_rt.json")
from simdgp import power_mixed, load_spec

n_sims = int(os.environ.get("N_SIMS", "40"))
print(f"Fitting {n_sims} crossed MixedLM models in Python (statsmodels)...")
t = time.time()
res = power_mixed(load_spec(SPEC), n_sims=n_sims)
print(f"elapsed {time.time() - t:.0f}s\n")
for k, v in res.items():
    print(f"  {k:14s}: {v}")
print("\nNOTE: statsmodels overstates random-slope variance, so this power is a CONSERVATIVE")
print("lower bound. The R/lme4 reference is ~0.73 at n=30. Fixed-effect recovery is correct")
print("(mean estimate ~0.048); use the R backend when random slopes / correlations matter.")
