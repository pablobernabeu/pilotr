"""pilotr, for simulating experimental and behavioural data from a portable design specification."""
from .core import RNG, as241, inv_logit
from .simulate import simulate, load_spec, Dataset
from .power import power, power_curve, power_mixed
from .examples import pilotr_example

__version__ = "0.1.0"
__all__ = ["RNG", "as241", "inv_logit", "simulate", "load_spec", "Dataset",
           "power", "power_curve", "power_mixed", "pilotr_example"]
