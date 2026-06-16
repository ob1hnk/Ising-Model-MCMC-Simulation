import numpy as np
import pytest
from ising.model import IsingModel


def test_energy_all_up():
    """All-up config: each spin has 4 neighbors, each bond counted once."""
    L = 3
    model = IsingModel(L, J=1.0)
    model.spins = np.ones((L, L), dtype=int)
    # 2 * L^2 bonds (right + down, periodic), each contributing -J
    expected = -2.0 * L * L
    assert np.isclose(model.energy(), expected)


def test_delta_energy_consistency():
    """delta_energy should match brute-force energy difference."""
    model = IsingModel(5, seed=42)
    x, y = 2, 3
    E_before = model.energy()
    dE = model.delta_energy(x, y)
    model.flip(x, y)
    E_after = model.energy()
    assert np.isclose(dE, E_after - E_before)
