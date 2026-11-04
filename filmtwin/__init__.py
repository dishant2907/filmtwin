"""FilmTwin: a physics-based PECVD chamber twin and process engineering toolkit.

Modules: 1 chamber twin, 2 metrology, 3 detection and diagnosis, 4 DOE, SPC and run-to-run control.
"""
from .chamber import ChamberSimulator, run_scenario
from .faults import FaultEvent
from .physics import Condition, Hardware, Recipe, deposit
from .scenarios import SCENARIOS

__version__ = "0.4.0"
__all__ = ["ChamberSimulator", "run_scenario", "FaultEvent", "Condition", "Hardware",
           "Recipe", "deposit", "SCENARIOS"]
