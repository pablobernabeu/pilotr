"""pilotr, for simulating experimental and behavioural data from a portable design specification."""
from .core import RNG, as241, inv_logit, replicate_seeds
from .simulate import simulate, load_spec, Dataset
from .validate import validate_spec, SPEC_VERSION
from .power import power, power_curve, power_mixed
from .examples import pilotr_example

__version__ = "0.3.0"
__all__ = ["RNG", "as241", "inv_logit", "replicate_seeds", "simulate", "load_spec", "Dataset",
           "validate_spec", "SPEC_VERSION", "power", "power_curve", "power_mixed",
           "pilotr_example"]
