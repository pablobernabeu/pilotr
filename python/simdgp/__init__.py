"""simdgp -- simulate experimental & behavioral data from a portable design spec."""
from .core import RNG, as241, inv_logit
from .simulate import simulate, load_spec, Dataset
from .power import power

__version__ = "0.1.0"
__all__ = ["RNG", "as241", "inv_logit", "simulate", "load_spec", "Dataset", "power"]
