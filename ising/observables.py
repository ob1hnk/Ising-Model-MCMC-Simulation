import numpy as np
from .model import IsingModel


def energy_per_spin(model: IsingModel) -> float:
    return model.energy() / model.N


def magnetization_squared(model: IsingModel) -> float:
    m = np.sum(model.spins) / model.N
    return float(m ** 2)
