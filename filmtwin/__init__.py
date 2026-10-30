"""FilmTwin: a physics-based PECVD chamber twin and process engineering toolkit.

Module 1 (this release): chamber simulator, faults, scenarios, datasets.
"""
from .chamber import ChamberSimulator, run_scenario
from .faults import FaultEvent
from .physics import Condition, Hardware, Recipe, deposit
from .scenarios import SCENARIOS

__version__ = "0.1.0"
__all__ = ["ChamberSimulator", "run_scenario", "FaultEvent", "Condition", "Hardware",
           "Recipe", "deposit", "SCENARIOS"]
