from pystics.simulation import run_pystics_simulation
import importlib.metadata

__version__ = importlib.metadata.version("mypackage")
__all__ = ["run_pystics_simulation"]