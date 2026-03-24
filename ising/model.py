import numpy as np


class IsingModel:
    """2D Ising model on an L x L lattice with periodic boundary conditions."""

    def __init__(self, L: int, J: float = 1.0, seed: int | None = None):
        self.L = L
        self.N = L * L
        self.J = J
        rng = np.random.default_rng(seed)   # random number generator
        self.spins = rng.choice([-1, 1], size=(L, L))   # LxL array, each index is -1 or 1 randomly

    def energy(self) -> float:
        """Total energy of the current configuration. O(N)."""
        s = self.spins
        # consider only right and down neighbors
        # periodic boundary condition is handled with np.roll()
        interaction = (
            np.roll(s, -1, axis=1)  # right neighbor
            + np.roll(s, -1, axis=0)  # down neighbor
        )
        return -self.J * float(np.sum(s * interaction))

    def delta_energy(self, x: int, y: int) -> float:
        """Energy change from flipping spin at (x, y). O(1)."""
        s = self.spins
        L = self.L
        neighbors_sum = (
            s[(x + 1) % L, y]
            + s[(x - 1) % L, y]
            + s[x, (y + 1) % L]
            + s[x, (y - 1) % L]
        )
        # multiply 2.0 because { unflipped = - flipped } and we are solving for delta value
        return 2.0 * self.J * s[x, y] * neighbors_sum

    def flip(self, x: int, y: int) -> None:
        self.spins[x, y] *= -1
