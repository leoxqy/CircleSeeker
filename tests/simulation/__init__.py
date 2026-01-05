"""
CircleSeeker Simulation and Validation Module

This module provides tools for:
- Generating simulated eccDNA data (reference genome + reads)
- Running validation against ground truth
- Calculating recall, precision, F1, FP, FN metrics

Usage:
    # Generate simulation data
    python -m tests.simulation.eccdna_simulator -o simulation_data

    # Run validation
    python -m tests.simulation.validation_metrics -s simulation_data -r results

    # End-to-end pipeline
    python -m tests.simulation.run_validation
"""

from .eccdna_simulator import (
    SimulationConfig,
    EccDNARecord,
    GenomicRegion,
    run_simulation,
)

from .validation_metrics import (
    ValidationMetrics,
    ValidationCalculator,
    run_validation,
)

__all__ = [
    "SimulationConfig",
    "EccDNARecord",
    "GenomicRegion",
    "run_simulation",
    "ValidationMetrics",
    "ValidationCalculator",
    "run_validation",
]
