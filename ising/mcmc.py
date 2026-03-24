import numpy as np
from .model import IsingModel


def run_step(model: IsingModel, beta: float, rng: np.random.Generator) -> None:
    """Single Metropolis-Hastings step: propose a random spin flip."""
    x = rng.integers(0, model.L)
    y = rng.integers(0, model.L)
    dE = model.delta_energy(x, y)
    if rng.random() < np.exp(-beta * dE):
        model.flip(x, y)


def run_sweep(model: IsingModel, beta: float, rng: np.random.Generator) -> None:
    """One sweep = N Metropolis-Hastings steps."""
    for _ in range(model.N):
        run_step(model, beta, rng)
